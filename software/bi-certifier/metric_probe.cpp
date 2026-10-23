#include "gradient_embedder.h"
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
  if (argc < 3 || argc > 7)
    return 2;
  size_t limit = argc > 3 ? std::stoull(argv[3]) : 50000;
  double tolerance = argc > 4 ? std::stod(argv[4]) : 1e-9;
  double factor = argc > 5 ? std::stod(argv[5]) : 1.5;
  size_t starts = argc > 6 ? std::stoull(argv[6]) : 3;
  if (starts < 1 || starts > 3 || tolerance <= 0 || !std::isfinite(tolerance) ||
      factor <= 1 || !std::isfinite(factor))
    return 2;
  std::ifstream in(argv[1]);
  std::ofstream out(argv[2]);
  out << std::setprecision(17);
  size_t cases = 0;
  if (!(in >> cases) || !out)
    return 3;
  out << "[\n";
  for (size_t test = 0; test < cases; test++) {
    size_t V = 0, F = 0;
    if (!(in >> V >> F))
      return 3;
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
    size_t numberOfOrbits = 0;
    if (!(in >> numberOfOrbits))
      return 3;
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
    std::unique_ptr<crp::OrbitGradient> selected;
    crp::GradientResult result;
    vector<crp::GradientResult> attempts;
    vector<double> factors = {factor, 2., 1.2};
    size_t selectedAttempt = 0;
    for (size_t attempt = 0; attempt < starts; attempt++) {
      auto trial = std::make_unique<crp::OrbitGradient>(mesh, L, orbits);
      crp::GradientResult state;
      try {
        state = trial->solve(limit, tolerance, factors[attempt]);
      } catch (const std::exception &e) {
        state = trial->result;
        state.success = false;
        state.reason = e.what();
      }
      attempts.push_back(state);
      if (!selected || state.success ||
          (state.initialized &&
           (!result.initialized || state.energy < result.energy))) {
        selected = std::move(trial);
        result = state;
        selectedAttempt = attempt;
      }
      if (result.success)
        break;
    }
    auto &optimizer = *selected;
    auto &embedder = *optimizer.geometry;
    bool init = result.initialized, success = result.success;
    if (init)
      embedder.refreshVertexCoordinates();
    vector<vector<size_t>> flips;
    for (auto f : embedder.history)
      flips.emplace_back(f.vertices.begin(), f.vertices.end());
    if (test)
      out << ",\n";
    out << "{\"case\":" << test << ",\"V\":" << V
        << ",\"initialized\":" << (init ? "true" : "false")
        << ",\"success\":" << (success ? "true" : "false") << ",\"seconds\":"
        << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         start)
               .count();
    out << ",\"optimizer\":\"orbit-squared-radius-gradient-descent\","
           "\"iterations\":"
        << result.iterations << ",\"rejected_steps\":" << result.rejected
        << ",\"flip_rollbacks\":" << result.rollbacks
        << ",\"initial_radius\":" << result.initialRadius
        << ",\"initial_curvature_energy\":" << result.initialEnergy
        << ",\"curvature_energy\":" << result.energy
        << ",\"max_radial_defect\":" << result.maxDefect
        << ",\"max_accepted_energy_increase\":" << result.maxIncrease
        << ",\"reason\":\"" << result.reason << "\"";
    out << ",\"selected_attempt\":" << selectedAttempt << ",\"attempts\":[";
    for (size_t a = 0; a < attempts.size(); a++) {
      auto t = attempts[a];
      if (a)
        out << ",";
      out << "{\"factor\":" << factors[a]
          << ",\"initial_radius\":" << t.initialRadius
          << ",\"initial_energy\":" << t.initialEnergy
          << ",\"final_energy\":" << t.energy
          << ",\"max_radial_defect\":" << t.maxDefect
          << ",\"iterations\":" << t.iterations
          << ",\"success\":" << (t.success ? "true" : "false")
          << ",\"reason\":\"" << t.reason << "\"}";
    }
    out << "],\"vertex_orbits\":[";
    for (size_t a = 0; a < orbits.size(); a++) {
      if (a)
        out << ",";
      out << "[";
      for (size_t j = 0; j < orbits[a].size(); j++) {
        if (j)
          out << ",";
        out << orbits[a][j];
      }
      out << "]";
    }
    out << "],\"energy_trace\":[";
    for (size_t a = 0; a < result.trace.size(); a++) {
      auto s = result.trace[a];
      if (a)
        out << ",";
      out << "[" << s.iteration << "," << s.energy << "," << s.defect << ","
          << s.step << "]";
    }
    out << "]";
    if (init) {
      out << ",\"radii\":[";
      for (size_t i = 0; i < V; i++) {
        if (i)
          out << ",";
        out << embedder.r[embedder.currentMesh->vertex(i)];
      }
      out << "],\"squared_radii\":[";
      for (size_t i = 0; i < V; i++) {
        if (i)
          out << ",";
        double r = embedder.r[embedder.currentMesh->vertex(i)];
        out << r * r;
      }
      out << "],\"radial_curvatures\":[";
      for (size_t i = 0; i < V; i++) {
        if (i)
          out << ",";
        out << embedder.kappa(embedder.currentMesh->vertex(i));
      }
      auto apex = embedder.apex;
      out << "],\"apex\":[" << apex.x << "," << apex.y << "," << apex.z << "]";
    }
    out << ",\"flips\":[";
    for (size_t i = 0; i < flips.size(); i++) {
      if (i)
        out << ",";
      out << "[";
      for (size_t j = 0; j < 4; j++) {
        if (j)
          out << ",";
        out << flips[i][j];
      }
      out << "]";
    }
    out << "]";
    if (init) {
      double mismatch = 0, metricError = 0;
      out << ",\"positions\":[";
      for (size_t i = 0; i < V; i++) {
        if (i)
          out << ",";
        auto p = embedder.embedding[embedder.currentMesh->vertex(i)];
        out << "[" << p.x << "," << p.y << "," << p.z << "]";
      }
      out << "],\"triangles\":[";
      bool first = true;
      for (Face f : embedder.currentMesh->faces()) {
        if (!first)
          out << ",";
        first = false;
        out << "[";
        bool f2 = true;
        for (Vertex v : f.adjacentVertices()) {
          if (!f2)
            out << ",";
          f2 = false;
          out << v.getIndex();
        }
        out << "]";
      }
      out << "],\"edges\":[";
      first = true;
      for (Edge e : embedder.currentMesh->edges()) {
        auto a = e.halfedge().vertex(), b = e.halfedge().twin().vertex();
        double l = embedder.intrinsicTriangulation->edgeLengths[e];
        metricError = std::max(
            metricError,
            std::abs((embedder.embedding[a] - embedder.embedding[b]).norm() -
                     l));
        if (!first)
          out << ",";
        first = false;
        out << "[" << a.getIndex() << "," << b.getIndex() << "," << l << "]";
      }
      for (auto item : lengths) {
        auto a = embedder.currentMesh->vertex(item.first.first),
             b = embedder.currentMesh->vertex(item.first.second);
        mismatch = std::max(
            mismatch,
            std::abs((embedder.embedding[a] - embedder.embedding[b]).norm() -
                     item.second));
      }
      out << "],\"max_output_edge_error\":" << metricError
          << ",\"max_original_face_pair_error\":" << mismatch;
    }
    out << "}";
    out.flush();
  }
  out << "\n]\n";
}
