#include "../gradient_embedder.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::vector;

int main(int argc, char **argv) {
  if (argc < 3 || argc > 6)
    return 2;
  size_t limit = argc > 3 ? std::stoull(argv[3]) : 50000;
  double tolerance = argc > 4 ? std::stod(argv[4]) : 1e-9;
  double factor = argc > 5 ? std::stod(argv[5]) : 1.5;
  std::ifstream in(argv[1]);
  std::ofstream out(argv[2]);
  out << std::setprecision(17);
  size_t cases;
  in >> cases;
  out << "[\n";
  for (size_t test = 0; test < cases; test++) {
    size_t V, F;
    in >> V >> F;
    vector<vector<size_t>> polygons(F), triangles;
    std::map<std::pair<size_t, size_t>, double> lengths;
    for (auto &f : polygons) {
      size_t m;
      in >> m;
      f.resize(m);
      for (auto &v : f)
        in >> v;
      for (size_t i = 0; i < m; i++)
        for (size_t j = i + 1; j < m; j++) {
          size_t d = std::min(j - i, m - (j - i));
          double l = std::sin(M_PI * d / m) / std::sin(M_PI / m);
          lengths[std::minmax(f[i], f[j])] = l;
        }
      for (size_t i = 1; i + 1 < m; i++)
        triangles.push_back({f[0], f[i], f[i + 1]});
    }
    size_t numberOfOrbits;
    in >> numberOfOrbits;
    vector<vector<size_t>> orbits(numberOfOrbits);
    for (auto &orbit : orbits) {
      size_t n;
      in >> n;
      orbit.resize(n);
      for (auto &v : orbit)
        in >> v;
    }
    if (!in)
      return 3;
    auto start = std::chrono::steady_clock::now();
    ManifoldSurfaceMesh mesh(triangles);
    EdgeData<double> L(mesh);
    for (Edge e : mesh.edges())
      L[e] = lengths.at(std::minmax(e.halfedge().vertex().getIndex(),
                                    e.halfedge().twin().vertex().getIndex()));

    crp::OrbitGradient optimizer(mesh, L, orbits);
    if (!optimizer.initialize(1.5))
      return 4;
    auto q = optimizer.q;
    auto g = optimizer.gradient(optimizer.curvature());
    double error = 0;
    for (size_t a = 0; a < orbits.size(); a++) {
      double h = 1e-6 * std::max(1., std::abs(q[a]));
      crp::Vec plus = q, minus = q;
      plus[a] += h;
      minus[a] -= h;
      optimizer.assign(plus);
      double ep = .5 * optimizer.curvature().squaredNorm();
      optimizer.assign(minus);
      double em = .5 * optimizer.curvature().squaredNorm();
      double fd = (ep - em) / (2 * h);
      error = std::max(error, std::abs(fd - g[a]) /
                                  std::max({1., std::abs(fd), std::abs(g[a])}));
    }
    optimizer.assign(q);
    // Rebuilding from the accepted flip history must preserve the geometry.
    auto history = optimizer.geometry->history;
    auto before = optimizer.curvature();
    optimizer.restore(history, q);
    double replay = (before - optimizer.curvature()).cwiseAbs().maxCoeff();
    if (test)
      out << ",\n";
    out << "{\"gradient_relative_error\":" << error
        << ",\"replay_error\":" << replay << "}";
  }
  out << "\n]\n";
}
