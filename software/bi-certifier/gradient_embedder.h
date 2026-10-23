#pragma once
#include "geometrycentral/surface/embed_convex.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace crp {
using namespace geometrycentral;
using namespace geometrycentral::surface;
using Vec = Eigen::VectorXd;

struct Flip {
  size_t edge;
  std::array<size_t, 4> vertices;
};

// Reuse BI's intrinsic metric, curvature, and analytic curvature derivative.
// Neither embed(), stepEmbedding(), nor solveNewton() is called.
class RadialGeometry : public ConvexEmbedder {
public:
  std::vector<Flip> history;
  RadialGeometry(ManifoldSurfaceMesh &mesh, EdgeData<double> &lengths,
                 EmbedConvexOptions options)
      : ConvexEmbedder(mesh, lengths, options) {
    r = VertexData<double>(*currentMesh);
    intrinsicTriangulation->edgeFlipCallbackList.push_back([this](Edge e) {
      auto h = e.halfedge(), t = h.twin();
      history.push_back({e.getIndex(),
                         {h.next().next().vertex().getIndex(),
                          t.next().next().vertex().getIndex(),
                          h.vertex().getIndex(), t.vertex().getIndex()}});
    });
  }
  using ConvexEmbedder::buildHessian;
  using ConvexEmbedder::kappa;

  bool validTetrahedra() {
    for (Vertex v : currentMesh->vertices())
      if (!std::isfinite(r[v]) || r[v] <= 0)
        return false;
    for (Face f : currentMesh->faces()) {
      auto h = f.halfedge();
      auto i = h.vertex(), j = h.next().vertex(), k = h.next().next().vertex();
      long double x = (long double)r[i] * r[i], y = (long double)r[j] * r[j],
                  z = (long double)r[k] * r[k];
      long double a = l(h), b = l(h.next()), c = l(h.next().next());
      long double u = (x + y - a * a) / 2, v = (y + z - b * b) / 2,
                  w = (z + x - c * c) / 2;
      long double minor = x * y - u * u, det = x * y * z + 2 * u * v * w -
                                               x * v * v - y * w * w -
                                               z * u * u;
      long double scale = std::max({x, y, z, a * a, b * b, c * c});
      if (!(minor > 1e-16L * scale * scale &&
            det > 1e-16L * scale * scale * scale))
        return false;
    }
    return true;
  }

  double initialize(double factor) {
    if (currentMesh->hasBoundary() || !currentMesh->isTriangular() ||
        currentMesh->genus() != 0 || !isMetricConvex())
      return 0;
    intrinsicTriangulation->flipToDelaunay();
    double R = initialRadius() * factor;
    for (int attempt = 0; attempt < 40; attempt++, R *= 1.5) {
      for (Vertex v : currentMesh->vertices())
        r[v] = R;
      if (!validTetrahedra())
        continue;
      bool ok = isEmbeddingConvex();
      for (Vertex v : currentMesh->vertices()) {
        double k = kappa(v), d = delta(v);
        if (!std::isfinite(k) || k < 0 || k > d)
          ok = false;
      }
      if (ok)
        return R;
    }
    return 0;
  }

  bool convexify() {
    std::queue<Edge> queue;
    EdgeData<bool> queued(*currentMesh, true);
    for (Edge e : currentMesh->edges())
      queue.push(e);
    size_t flips = 0, limit = 1000 + 20 * currentMesh->nEdges();
    while (!queue.empty()) {
      Edge e = queue.front();
      queue.pop();
      queued[e] = false;
      double angle = theta(e);
      if (!std::isfinite(angle))
        return false;
      if (angle <= M_PI + options.edgeConvexityTolerance)
        continue;
      if (++flips > limit || !intrinsicTriangulation->flipEdgeIfPossible(e))
        return false;
      auto h = e.halfedge(), t = h.twin();
      for (Edge neighbor : {h.next().edge(), h.next().next().edge(),
                            t.next().edge(), t.next().next().edge()})
        if (!queued[neighbor]) {
          queue.push(neighbor);
          queued[neighbor] = true;
        }
    }
    return validTetrahedra() && isEmbeddingConvex();
  }
};

struct Sample {
  size_t iteration;
  double energy, defect, step;
};
struct GradientResult {
  bool initialized = false, success = false;
  size_t iterations = 0, rejected = 0, rollbacks = 0;
  double initialRadius = 0, initialEnergy = 0, energy = 0, maxDefect = 0,
         maxIncrease = 0;
  std::string reason = "not initialized";
  std::vector<Sample> trace;
};

class OrbitGradient {
public:
  ManifoldSurfaceMesh &original;
  EdgeData<double> &lengths;
  EmbedConvexOptions options;
  std::vector<std::vector<size_t>> orbits;
  std::unique_ptr<RadialGeometry> geometry;
  Vec q, weights;
  double scale = 1;
  GradientResult result;

  OrbitGradient(ManifoldSurfaceMesh &mesh, EdgeData<double> &L,
                const std::vector<std::vector<size_t>> &groups)
      : original(mesh), lengths(L), orbits(groups) {
    options.edgeConvexityTolerance = 1e-10;
    std::vector<bool> seen(mesh.nVertices(), false);
    for (auto &orbit : orbits) {
      if (orbit.empty())
        throw std::invalid_argument("empty vertex orbit");
      for (size_t v : orbit) {
        if (v >= seen.size() || seen[v])
          throw std::invalid_argument("invalid orbit partition");
        seen[v] = true;
      }
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end())
      throw std::invalid_argument("incomplete orbit partition");
    geometry = std::make_unique<RadialGeometry>(mesh, L, options);
    weights = Vec(orbits.size());
    for (size_t a = 0; a < orbits.size(); a++)
      weights[a] = std::sqrt((double)orbits[a].size());
  }

  void assign(const Vec &values) {
    // Compute each radius once, then assign that very same double to its orbit.
    for (size_t a = 0; a < orbits.size(); a++) {
      double radius = std::sqrt(scale * values[a] / weights[a]);
      for (size_t v : orbits[a])
        geometry->r[geometry->currentMesh->vertex(v)] = radius;
    }
  }
  void restore(const std::vector<Flip> &history, const Vec &values) {
    geometry = std::make_unique<RadialGeometry>(original, lengths, options);
    for (const Flip &expected : history) {
      if (!geometry->intrinsicTriangulation->flipEdgeIfPossible(
              geometry->currentMesh->edge(expected.edge)) ||
          geometry->history.back().vertices != expected.vertices)
        throw std::runtime_error(
            "cannot replay accepted intrinsic triangulation");
    }
    assign(values);
  }
  Vec curvature() {
    Vec k(original.nVertices());
    for (Vertex v : geometry->currentMesh->vertices())
      k[v.getIndex()] = geometry->kappa(v);
    return k;
  }
  Vec gradient(const Vec &k) {
    SparseMatrix<double> J;
    bool failed;
    std::tie(J, failed) = geometry->buildHessian();
    if (failed)
      throw std::runtime_error("invalid curvature derivative");
    Vec full = J.transpose() * k, g = Vec::Zero(orbits.size());
    for (size_t a = 0; a < orbits.size(); a++)
      for (size_t v : orbits[a])
        g[a] +=
            full[v] * scale /
            (2 * geometry->r[geometry->currentMesh->vertex(v)] * weights[a]);
    if (!g.allFinite())
      throw std::runtime_error("nonfinite curvature gradient");
    return g;
  }
  bool initialize(double factor) {
    result.initialRadius = geometry->initialize(factor);
    if (!(result.initialRadius > 0))
      return false;
    scale = result.initialRadius * result.initialRadius;
    q = weights;
    assign(q);
    result.initialized = true;
    Vec k = curvature();
    result.initialEnergy = .5 * k.squaredNorm();
    return true;
  }
  GradientResult solve(size_t limit, double tolerance, double factor) {
    if (!initialize(factor))
      return result;
    std::cerr << "Initialized " << original.nVertices() << " vertices in "
              << orbits.size() << " orbits; common radius "
              << result.initialRadius << std::endl;
    Vec oldQ, oldG;
    double previousEnergy = result.initialEnergy;
    for (size_t iteration = 0; iteration <= limit; iteration++) {
      Vec k = curvature();
      double energy = .5 * k.squaredNorm(), defect = k.cwiseAbs().maxCoeff();
      if (iteration % 1000 == 0)
        std::cerr << "Gradient step " << iteration << ", E=" << energy
                  << ", max curvature=" << defect << std::endl;
      result.iterations = iteration;
      result.energy = energy;
      result.maxDefect = defect;
      if (iteration % 100 == 0 || defect <= tolerance || iteration == limit)
        result.trace.push_back({iteration, energy, defect, 0});
      if (!k.allFinite()) {
        result.reason = "nonfinite radial curvature";
        break;
      }
      if (defect <= tolerance) {
        result.success = true;
        result.reason = "radial curvature tolerance reached";
        break;
      }
      if (iteration == limit) {
        result.reason = "gradient iteration budget exhausted";
        break;
      }
      Vec g = gradient(k);
      if (g.norm() < 1e-20) {
        result.reason = "stationary with nonzero radial curvature";
        break;
      }
      double alpha = std::min(1., .1 * std::max(1., q.norm()) / g.norm());
      if (iteration && oldQ.size() == q.size()) {
        Vec s = q - oldQ, y = g - oldG;
        double sy = s.dot(y);
        if (sy > 0 && y.squaredNorm() > 0) {
          // Alternating Barzilai--Borwein scalar step lengths; direction stays
          // exactly the negative gradient. No Newton or quasi-Newton solve.
          alpha = (iteration % 2) ? s.squaredNorm() / sy : sy / y.squaredNorm();
          alpha = std::clamp(alpha, 1e-12, 1e8);
        }
      }
      auto history = geometry->history;
      bool accepted = false;
      for (int backtrack = 0; backtrack < 64; backtrack++, alpha *= .5) {
        Vec candidate = q - alpha * g;
        if ((candidate - q).squaredNorm() == 0)
          break;
        if (!candidate.allFinite() || candidate.minCoeff() <= 0) {
          result.rejected++;
          continue;
        }
        assign(candidate);
        bool valid = geometry->validTetrahedra() && geometry->convexify();
        double candidateEnergy =
            valid ? .5 * curvature().squaredNorm() : INFINITY;
        if (std::isfinite(candidateEnergy) && candidateEnergy < energy &&
            candidateEnergy <= energy - 1e-4 * alpha * g.squaredNorm()) {
          oldQ = q;
          oldG = g;
          q = candidate;
          accepted = true;
          result.maxIncrease =
              std::max(result.maxIncrease, candidateEnergy - previousEnergy);
          previousEnergy = candidateEnergy;
          if (iteration % 100 == 0)
            result.trace.back().step = alpha;
          break;
        }
        result.rejected++;
        if (geometry->history.size() != history.size()) {
          restore(history, q);
          result.rollbacks++;
        } else
          assign(q);
      }
      if (!accepted) {
        result.reason = "line search stalled with nonzero radial curvature";
        break;
      }
    }
    // Values in the final state, including any unsuccessful last iteration.
    Vec k = curvature();
    result.energy = .5 * k.squaredNorm();
    result.maxDefect = k.cwiseAbs().maxCoeff();
    if (result.trace.empty() ||
        result.trace.back().iteration != result.iterations)
      result.trace.push_back(
          {result.iterations, result.energy, result.maxDefect, 0});
    return result;
  }
};
} // namespace crp
