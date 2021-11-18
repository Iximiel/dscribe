#include "weighting.h"
#include <functional>
namespace py = pybind11;
using std::function;
using std::string;

/**
 * Used to calculate the Gaussian weights for each neighbouring atom. Provide
 * either r1s (=r) or r2s (=r^2) and use the boolean "squared" to indicate if
 * r1s should be calculated from r2s.
 */
void getWeights(const int size, double *r1s, const double *r2s,
                const bool squared, const py::dict &weighting,
                double *weights) {
  // No weighting specified
  if (!weighting.contains("function") && !weighting.contains("w0")) {
    for (int i = 0; i < size; i++) {
      weights[i] = 1;
    }
  } else {
    // No weighting function, only w0
    if (squared) {
      for (int i = 0; i < size; i++) {
        r1s[i] = sqrt(r2s[i]);
      }
    }
    getWeights(size, r1s, weighting, weights);
  }
}

/**
 * Used to calculate the Gaussian weights for each neighbouring atom. Provide
 * either r1s (=r) or r2s (=r^2) and use the boolean "squared" to indicate if
 * r1s should be calculated from r2s.
 */
void getWeights(const int size, const double *r1s, const py::dict &weighting,
                double *weights) {
  // No weighting specified
  if (!weighting.contains("function") && !weighting.contains("w0")) {
    for (int i = 0; i < size; i++) {
      weights[i] = 1;
    }
  } else {
    // No weighting function, only w0
    if (!weighting.contains("function") && weighting.contains("w0")) {
      double w0 = weighting["w0"].cast<double>();
      for (int i = 0; i < size; i++) {
        double r = r1s[i];
        if (r == 0) {
          weights[i] = w0;
        } else {
          weights[i] = 1;
        }
      }
    } else {
      // func is initialized by this lambda that is called a it is declared
      const function<double(double)> func =
          [](const py::dict &weightingData) -> function<double(double)> {
        string fname = weightingData["function"].cast<string>();
        if (fname == "poly") {
          double r0 = weightingData["r0"].cast<double>();
          double c = weightingData["c"].cast<double>();
          double m = weightingData["m"].cast<double>();
          return [r0, c, m](double r) -> double {
            return weightPoly(r, r0, c, m);
          };
        } else if (fname == "pow") {
          double r0 = weightingData["r0"].cast<double>();
          double c = weightingData["c"].cast<double>();
          double d = weightingData["d"].cast<double>();
          double m = weightingData["m"].cast<double>();
          return [r0, c, d, m](double r) -> double {
            return weightPow(r, r0, c, d, m);
          };
        } else if (fname == "exp") {
          double r0 = weightingData["r0"].cast<double>();
          double c = weightingData["c"].cast<double>();
          double d = weightingData["d"].cast<double>();
          return
              [r0, c, d](double r) -> double { return weightExp(r, r0, c, d); };
        }
      }(weighting);
      // Weighting function and w0
      if (weighting.contains("w0")) {
        double w0 = weighting["w0"].cast<double>();
        for (int i = 0; i < size; i++) {
          if (r1s[i] == 0) {
            weights[i] = w0;
          } else {
            weights[i] = func(r1s[i]);
          }
        }
        // Weighting function only
      } else {
        for (int i = 0; i < size; i++) {
          weights[i] = func(r1s[i]);
        }
      }
    }
  }
}
