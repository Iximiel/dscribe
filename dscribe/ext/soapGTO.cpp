/*Copyright 2019 DScribe developers

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#include "soapGTO.h"
#include "celllist.h"
#include "soapGTO_getFactors.h"
#include "weighting.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <set>
#include <string>

#define PI2 9.86960440108936
#define PI 3.141592653589793238
#define PI3 31.00627668029982
#define PIHalf 1.57079632679490
/*===========================================================*/
inline int getCrosNumD(int n) { return n * (n + 1) / 2; }
/*================================================================*/
inline void getDeltasD(double *x, double *y, double *z,
                       const py::array_t<double> &positions, const double ix,
                       const double iy, const double iz,
                       const vector<int> &indices) {

  int count = 0;
  auto pos = positions.unchecked<2>();
  for (const int &idx : indices) {
    x[count] = pos(idx, 0) - ix;
    y[count] = pos(idx, 1) - iy;
    z[count] = pos(idx, 2) - iz;
    count++;
  };
}
/*================================================================*/
inline void getRsZsD(double *x, double *x2, double *x4, double *x6, double *x8,
                     double *x10, double *x12, double *x14, double *x16,
                     double *x18, double *y, double *y2, double *y4, double *y6,
                     double *y8, double *y10, double *y12, double *y14,
                     double *y16, double *y18, double *z, double *r2,
                     double *r4, double *r6, double *r8, double *r10,
                     double *r12, double *r14, double *r16, double *r18,
                     double *z2, double *z4, double *z6, double *z8,
                     double *z10, double *z12, double *z14, double *z16,
                     double *z18, double *r20, double *x20, double *y20,
                     double *z20, int size, int lMax) {
  double xx;
  double yy;
  double zz;
  double rr;

  for (int i = 0; i < size; i++) {
    xx = x[i] * x[i];
    yy = y[i] * y[i];
    zz = z[i] * z[i];
    rr = xx + yy + zz;
    x2[i] = xx;
    y2[i] = yy;
    z2[i] = zz;
    r2[i] = rr;

    if (lMax > 3) {
      r4[i] = r2[i] * r2[i];
      z4[i] = z2[i] * z2[i];
      x4[i] = x2[i] * x2[i];
      y4[i] = y2[i] * y2[i];
      if (lMax > 5) {
        r6[i] = r2[i] * r4[i];
        z6[i] = z2[i] * z4[i];
        x6[i] = x2[i] * x4[i];
        y6[i] = y2[i] * y4[i];
        if (lMax > 7) {
          r8[i] = r4[i] * r4[i];
          z8[i] = z4[i] * z4[i];
          x8[i] = x4[i] * x4[i];
          y8[i] = y4[i] * y4[i];
          if (lMax > 9) {
            x10[i] = x6[i] * x4[i];
            y10[i] = y6[i] * y4[i];
            z10[i] = z6[i] * z4[i];
            r10[i] = r6[i] * r4[i];
            if (lMax > 11) {
              x12[i] = x6[i] * x6[i];
              y12[i] = y6[i] * y6[i];
              r12[i] = r6[i] * r6[i];
              z12[i] = z6[i] * z6[i];
              if (lMax > 13) {
                x14[i] = x6[i] * x8[i];
                y14[i] = y6[i] * y8[i];
                r14[i] = r6[i] * r8[i];
                z14[i] = z6[i] * z8[i];
                if (lMax > 15) {
                  x16[i] = x8[i] * x8[i];
                  y16[i] = y8[i] * y8[i];
                  r16[i] = r8[i] * r8[i];
                  z16[i] = z8[i] * z8[i];
                  if (lMax > 17) {
                    x18[i] = x10[i] * x8[i];
                    y18[i] = y10[i] * y8[i];
                    r18[i] = r10[i] * r8[i];
                    z18[i] = z10[i] * z8[i];
                    if (lMax > 19) {
                      x20[i] = x10[i] * x10[i];
                      z20[i] = z10[i] * z10[i];
                      y20[i] = y10[i] * y10[i];
                      r20[i] = r10[i] * r10[i];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
/*================================================================*/
void getAlphaBetaD(double *aOa, double *bOa, double *alphas, double *betas,
                   int Ns, int lMax, double oOeta, double oOeta3O2) {

  int NsNs = Ns * Ns;
  double oneO1alpha;
  double oneO1alpha2;
  double oneO1alphaSqrt;
  double oneO1alphaSqrtX;

  for (int myL = 0; myL < lMax + 1; myL++) {
    for (int k = 0; k < Ns; k++) {
      oneO1alpha = 1.0 / (1.0 + oOeta * alphas[myL * Ns + k]);
      oneO1alphaSqrt = sqrt(oneO1alpha);
      aOa[myL * Ns + k] = -alphas[myL * Ns + k] * oneO1alpha;
      oneO1alpha2 = pow(oneO1alpha, myL + 1);
      oneO1alphaSqrtX = oneO1alphaSqrt * oneO1alpha2;
      for (int n = 0; n < Ns; n++) {
        bOa[myL * NsNs + n * Ns + k] =
            oOeta3O2 * betas[myL * NsNs + n * Ns + k] * oneO1alphaSqrtX;
      }
    }
  }
}
/*================================================================*/
/*==============================================================================================================================*/
void getCD(py::detail::unchecked_mutable_reference<double, 5> &CDevX_mu,
           py::detail::unchecked_mutable_reference<double, 5> &CDevY_mu,
           py::detail::unchecked_mutable_reference<double, 5> &CDevZ_mu,
           double *prCofDX, double *prCofDY, double *prCofDZ,
           py::detail::unchecked_mutable_reference<double, 4> &C_mu,
           double *preCoef, double *x, double *y, double *z, double *r2,
           double *weights, double *bOa, double *aOa, double *exes, int totalAN,
           int Asize, int Ns, int Ntypes, int lMax, int posI, int typeJ,
           const vector<int> &indices, bool return_derivatives) {
  if (Asize == 0) {
    return;
  }
  double sumMe = 0;
  int NsNs = Ns * Ns;
  int NsJ = ((lMax + 1) * (lMax + 1)) * Ns * typeJ;
  int LNsNs;
  int LNs;
  double preExp;
  double preVal;
  double preValX;
  double preValY;
  double preValZ;

  double preVal1;
  double preVal2;
  double preVal3;

  double preValX1;
  double preValY1;
  double preValZ1;

  double preValX2;
  double preValY2;
  double preValZ2;

  double preValX3;
  double preValY3;
  double preValZ3;
  double *preExponentArrya = (double *)malloc(Ns * Asize * sizeof(double));
  // l=0-------------------------------------------------------------------------------------------------
  int shift = 0;
  for (int k = 0; k < Ns; k++) {
    for (int i = 0; i < Asize; i++) {
      preExponentArrya[shift] =
          weights[i] * 1.5707963267948966 * exp(aOa[k] * r2[i]);
      shift++;
    }
  }

  shift = 0;
  for (int k = 0; k < Ns; k++) {
    sumMe = 0;
    for (int i = 0; i < Asize; i++) {
      preExp = preExponentArrya[shift];
      sumMe += preExp;
      shift++;
    }
    for (int n = 0; n < Ns; n++) {
      C_mu(posI, typeJ, n, 0) += bOa[n * Ns + k] * sumMe;
    }
  }
  shift = 0;
  if (return_derivatives) {
    for (int k = 0; k < Ns; k++) {
      for (int i = 0; i < Asize; i++) {
        preExp = preExponentArrya[shift];
        shift++;
        preVal = 2.0 * aOa[k] * preExp;
        preValX = preVal * x[i];
        preValY = preVal * y[i];
        preValZ = preVal * z[i];
        for (int n = 0; n < Ns; n++) {
          CDevX_mu(indices[i], posI, typeJ, n, 0) += bOa[n * Ns + k] * preValX;
          CDevY_mu(indices[i], posI, typeJ, n, 0) += bOa[n * Ns + k] * preValY;
          CDevZ_mu(indices[i], posI, typeJ, n, 0) += bOa[n * Ns + k] * preValZ;
        }
      }
    }
  }
  // l=1-------------------------------------------------------------------------------------------------
  if (lMax > 0) {
    LNsNs = NsNs;
    LNs = Ns;
    shift = 0;

    for (int k = 0; k < Ns; k++) {
      for (int i = 0; i < Asize; i++) {
        preExponentArrya[shift] =
            weights[i] * 2.7206990463849543 * exp(aOa[LNs + k] * r2[i]);
        shift++;
      }
    }

    shift = 0;
    double sumMe1;
    double sumMe2;
    double sumMe3;
    for (int k = 0; k < Ns; k++) {
      sumMe1 = 0;
      sumMe2 = 0;
      sumMe3 = 0;
      for (int i = 0; i < Asize; i++) {
        preExp = preExponentArrya[shift];
        sumMe1 += preExp * z[i];
        sumMe2 += preExp * x[i];
        sumMe3 += preExp * y[i];
        shift++;
      }
      for (int n = 0; n < Ns; n++) {
        C_mu(posI, typeJ, n, 1) += bOa[LNsNs + n * Ns + k] * sumMe1;
        C_mu(posI, typeJ, n, 2) += bOa[LNsNs + n * Ns + k] * sumMe2;
        C_mu(posI, typeJ, n, 3) += bOa[LNsNs + n * Ns + k] * sumMe3;
      }
    }

    if (return_derivatives) {
      shift = 0;
      for (int k = 0; k < Ns; k++) {
        for (int i = 0; i < Asize; i++) {
          preExp = preExponentArrya[shift];
          shift++;
          preVal = 2.0 * aOa[LNs + k] * preExp;
          if (return_derivatives) {
            preVal1 = preVal * z[i];
            preVal2 = preVal * x[i];
            preVal3 = preVal * y[i];

            preValX1 = preVal1 * x[i];
            preValY1 = preVal1 * y[i];
            preValZ1 = preVal1 * z[i] + preExp;

            preValX2 = preVal2 * x[i] + preExp;
            preValY2 = preVal2 * y[i];
            preValZ2 = preVal2 * z[i];

            preValX3 = preVal3 * x[i];
            preValY3 = preVal3 * y[i] + preExp;
            preValZ3 = preVal3 * z[i];
          }
          for (int n = 0; n < Ns; n++) {
            CDevX_mu(indices[i], posI, typeJ, n, 1) +=
                bOa[LNsNs + n * Ns + k] * preValX1;
            CDevY_mu(indices[i], posI, typeJ, n, 1) +=
                bOa[LNsNs + n * Ns + k] * preValY1;
            CDevZ_mu(indices[i], posI, typeJ, n, 1) +=
                bOa[LNsNs + n * Ns + k] * preValZ1;

            CDevX_mu(indices[i], posI, typeJ, n, 2) +=
                bOa[LNsNs + n * Ns + k] * preValX2;
            CDevY_mu(indices[i], posI, typeJ, n, 2) +=
                bOa[LNsNs + n * Ns + k] * preValY2;
            CDevZ_mu(indices[i], posI, typeJ, n, 2) +=
                bOa[LNsNs + n * Ns + k] * preValZ2;

            CDevX_mu(indices[i], posI, typeJ, n, 3) +=
                bOa[LNsNs + n * Ns + k] * preValX3;
            CDevY_mu(indices[i], posI, typeJ, n, 3) +=
                bOa[LNsNs + n * Ns + k] * preValY3;
            CDevZ_mu(indices[i], posI, typeJ, n, 3) +=
                bOa[LNsNs + n * Ns + k] * preValZ3;
          }
        }
      }
    }
  }
  // l>2------------------------------------------------------------------------------------------------------
  //
  //[NsTs100*i_center*totalAN + Ns100*jd*totalAN + buffShift*totalAN*Ns +
  // kd*totalAN + i_atom]
  //
  if (lMax > 1) {
    for (int restOfLs = 2; restOfLs <= lMax; restOfLs++) {
      LNsNs = restOfLs * NsNs;
      LNs = restOfLs * Ns;
      shift = 0;
      for (int k = 0; k < Ns; k++) {
        for (int i = 0; i < Asize; i++) {
          double expSholder = aOa[LNs + k] * r2[i];
          preExponentArrya[shift] = weights[i] * exp(expSholder);
          shift++;
        }
      }

      // double*  sumS = (double*)
      // malloc(sizeof(double)*(restOfLs+1)*(restOfLs+1))
      shift = 0;
      for (int k = 0; k < Ns; k++) {
        for (int m = restOfLs * restOfLs; m < (restOfLs + 1) * (restOfLs + 1);
             m++) {
          sumMe = 0;
          for (int i = 0; i < Asize; i++) {
            preExp = preExponentArrya[Asize * k + i];
            sumMe += preExp * preCoef[totalAN * (m - 4) + i];
            shift++;
          }
          for (int n = 0; n < Ns; n++) {
            C_mu(posI, typeJ, n, m) += bOa[LNsNs + n * Ns + k] * sumMe;
          }
        }
      }

      if (return_derivatives) {
        shift = 0;
        for (int k = 0; k < Ns; k++) {
          for (int i = 0; i < Asize; i++) {
            preExp = preExponentArrya[shift];
            shift++;
            for (int m = restOfLs * restOfLs;
                 m < (restOfLs + 1) * (restOfLs + 1); m++) {
              preVal =
                  2.0 * aOa[LNs + k] * preExp * preCoef[totalAN * (m - 4) + i];
              preValX = x[i] * preVal + preExp * prCofDX[totalAN * (m - 4) + i];
              preValY = y[i] * preVal + preExp * prCofDY[totalAN * (m - 4) + i];
              preValZ = z[i] * preVal + preExp * prCofDZ[totalAN * (m - 4) + i];
              for (int n = 0; n < Ns; n++) {
                CDevX_mu(indices[i], posI, typeJ, n, m) +=
                    bOa[LNsNs + n * Ns + k] * preValX;
                CDevY_mu(indices[i], posI, typeJ, n, m) +=
                    bOa[LNsNs + n * Ns + k] * preValY;
                CDevZ_mu(indices[i], posI, typeJ, n, m) +=
                    bOa[LNsNs + n * Ns + k] * preValZ;
              }
            }
          }
        }
      }
    }
  }
}
/*================================================================================================*/
/**
 * Used to calculate the partial power spectrum.
 */
void getPD(py::detail::unchecked_mutable_reference<double, 2> &descriptor_mu,
           py::detail::unchecked_reference<double, 4> &Cnnd_u, int Ns, int Ts,
           int nCenters, int lMax, bool crossover) {

  int NsTs100 = Ns * Ts * ((lMax + 1) * (lMax + 1));

  // The power spectrum is multiplied by an l-dependent prefactor that comes
  // from the normalization of the Wigner D matrices. This prefactor is
  // mentioned in the arrata of the original SOAP paper: On representing
  // chemical environments, Phys. Rev. B 87, 184115 (2013). Here the square
  // root of the prefactor in the dot-product kernel is used, so that after a
  // possible dot-product the full prefactor is recovered.
  for (int i = 0; i < nCenters; i++) {
    int shiftAll = 0;
    for (int j = 0; j < Ts; j++) {
      int jdLimit = crossover ? Ts : j + 1;
      for (int jd = j; jd < jdLimit; jd++) {
        for (int m = 0; m <= lMax; m++) {
          double prel;
          if (m > 1) {
            prel = PI * sqrt(8.0 / (2.0 * m + 1.0)) * PI3;
          } else {
            prel = PI * sqrt(8.0 / (2.0 * m + 1.0));
          }
          if (j == jd) {
            for (int k = 0; k < Ns; k++) {
              for (int kd = k; kd < Ns; kd++) {
                double buffDouble = 0;
                for (int buffShift = m * m; buffShift < (m + 1) * (m + 1);
                     buffShift++) {
                  buffDouble +=
                      Cnnd_u(i, j, k, buffShift) * Cnnd_u(i, jd, kd, buffShift);
                }
                descriptor_mu(i, shiftAll) = prel * buffDouble;
                shiftAll++;
              }
            }
          } else {
            for (int k = 0; k < Ns; k++) {
              for (int kd = 0; kd < Ns; kd++) {
                double buffDouble = 0;
                for (int buffShift = m * m; buffShift < (m + 1) * (m + 1);
                     buffShift++) {
                  buffDouble +=
                      Cnnd_u(i, j, k, buffShift) * Cnnd_u(i, jd, kd, buffShift);
                }
                descriptor_mu(i, shiftAll) = prel * buffDouble;
                shiftAll++;
              }
            }
          }
        } // end ifelse
      }
    }
  }
}
/*===========================================================================================*/
/**
 * Used to calculate the partial power spectrum derivatives.
 */
void getPDev(py::detail::unchecked_mutable_reference<double, 4> &derivatives_mu,
             py::detail::unchecked_reference<double, 2> &positions_u,
             py::detail::unchecked_reference<int, 1> &indices_u,
             CellList &cell_list,
             py::detail::unchecked_reference<double, 5> &CdevX_u,
             py::detail::unchecked_reference<double, 5> &CdevY_u,
             py::detail::unchecked_reference<double, 5> &CdevZ_u,
             py::detail::unchecked_reference<double, 4> &Cnnd_u, int Ns, int Ts,
             int nCenters, int lMax, bool crossover) {

  // Loop over all given atomic indices for which the derivatives should be
  // calculated for.
  for (int i_idx = 0; i_idx < indices_u.size(); ++i_idx) {
    int i_atom = indices_u(i_idx);

    // Get all neighbouring centers for the current atom
    double ix = positions_u(i_atom, 0);
    double iy = positions_u(i_atom, 1);
    double iz = positions_u(i_atom, 2);
    CellListResult result = cell_list.getNeighboursForPosition(ix, iy, iz);
    vector<int> indices = result.indices;

    // Loop through all neighbouring centers
    for (int j_idx = 0; j_idx < indices.size(); ++j_idx) {
      int i_center = indices[j_idx];
      int shiftAll = 0;
      for (int j = 0; j < Ts; j++) {
        int jdLimit = crossover ? Ts : j + 1;
        for (int jd = j; jd < jdLimit; jd++) {
          for (int m = 0; m <= lMax; m++) {
            double prel = m > 1 ? PI * sqrt(8.0 / (2.0 * m + 1.0)) * PI3
                                : PI * sqrt(8.0 / (2.0 * m + 1.0));
            if (j == jd) {
              for (int k = 0; k < Ns; k++) {
                for (int kd = k; kd < Ns; kd++) {
                  for (int buffShift = m * m; buffShift < (m + 1) * (m + 1);
                       buffShift++) {
                    if (abs(Cnnd_u(i_center, j, k, buffShift)) > 1e-8 ||
                        abs(Cnnd_u(i_center, j, kd, buffShift)) > 1e-8) {
                      derivatives_mu(i_center, i_idx, 0, shiftAll) +=
                          prel *
                          (Cnnd_u(i_center, j, k, buffShift) *
                               CdevX_u(i_atom, i_center, jd, kd, buffShift) +
                           Cnnd_u(i_center, j, kd, buffShift) *
                               CdevX_u(i_atom, i_center, jd, k, buffShift));
                      derivatives_mu(i_center, i_idx, 1, shiftAll) +=
                          prel *
                          (Cnnd_u(i_center, j, k, buffShift) *
                               CdevY_u(i_atom, i_center, jd, kd, buffShift) +
                           Cnnd_u(i_center, j, kd, buffShift) *
                               CdevY_u(i_atom, i_center, jd, k, buffShift));
                      derivatives_mu(i_center, i_idx, 2, shiftAll) +=
                          prel *
                          (Cnnd_u(i_center, j, k, buffShift) *
                               CdevZ_u(i_atom, i_center, jd, kd, buffShift) +
                           Cnnd_u(i_center, j, kd, buffShift) *
                               CdevZ_u(i_atom, i_center, jd, k, buffShift));
                    }
                  }
                  shiftAll++;
                }
              }
            } else {
              for (int k = 0; k < Ns; k++) {
                for (int kd = 0; kd < Ns; kd++) {
                  for (int buffShift = m * m; buffShift < (m + 1) * (m + 1);
                       buffShift++) {
                    if (abs(Cnnd_u(i_center, j, k, buffShift)) > 1e-8 ||
                        abs(Cnnd_u(i_center, jd, kd, buffShift)) > 1e-8) {
                      derivatives_mu(i_center, i_idx, 0, shiftAll) +=
                          prel *
                          (Cnnd_u(i_center, j, k, buffShift) *
                               CdevX_u(i_atom, i_center, jd, kd, buffShift) +
                           Cnnd_u(i_center, jd, kd, buffShift) *
                               CdevX_u(i_atom, i_center, j, k, buffShift));
                      derivatives_mu(i_center, i_idx, 1, shiftAll) +=
                          prel *
                          (Cnnd_u(i_center, j, k, buffShift) *
                               CdevY_u(i_atom, i_center, jd, kd, buffShift) +
                           Cnnd_u(i_center, jd, kd, buffShift) *
                               CdevY_u(i_atom, i_center, j, k, buffShift));
                      derivatives_mu(i_center, i_idx, 2, shiftAll) +=
                          prel *
                          (Cnnd_u(i_center, j, k, buffShift) *
                               CdevZ_u(i_atom, i_center, jd, kd, buffShift) +
                           Cnnd_u(i_center, jd, kd, buffShift) *
                               CdevZ_u(i_atom, i_center, j, k, buffShift));
                    }
                  }
                  shiftAll++;
                }
              }
            }
          }
        }
      }
    }
  }
}
/*=================================================================================================================================================================*/
void soapGTO(py::array_t<double> derivatives, py::array_t<double> descriptor,
             py::array_t<double> cdevX, py::array_t<double> cdevY,
             py::array_t<double> cdevZ, py::array_t<double> positions,
             py::array_t<double> centers, py::array_t<double> alphasArr,
             py::array_t<double> betasArr, py::array_t<int> atomicNumbersArr,
             py::array_t<int> orderedSpeciesArr, const double rCut,
             const double cutoffPadding, const int nMax, const int lMax,
             const double eta, py::dict weighting, const bool crossover,
             string average, py::array_t<int> indices,
             const bool return_descriptor, const bool return_derivatives,
             CellList cell_list_atoms) {
  const int totalAN = atomicNumbersArr.shape(0);
  const int nCenters = centers.shape(0);
  auto derivatives_mu = derivatives.mutable_unchecked<4>();
  auto atomicNumbers = atomicNumbersArr.unchecked<1>();
  auto species = orderedSpeciesArr.unchecked<1>();
  int nSpecies = orderedSpeciesArr.shape(0);
  auto indices_u = indices.unchecked<1>();
  double *alphas = (double *)alphasArr.request().ptr;
  double *betas = (double *)betasArr.request().ptr;
  double oOeta = 1.0 / eta;
  double oOeta3O2 = sqrt(oOeta * oOeta * oOeta);
  double nMax2 = nMax * nMax;
  auto centers_u = centers.unchecked<2>();
  auto positions_u = positions.unchecked<2>();
  const int nFeatures =
      crossover ? (nSpecies * nMax) * (nSpecies * nMax + 1) / 2 * (lMax + 1)
                : nSpecies * (lMax + 1) * ((nMax + 1) * nMax) / 2;
  double *weights = (double *)malloc(sizeof(double) * totalAN);
  constexpr int to2 = 0;
  constexpr int to4 = 1;
  constexpr int to6 = 2;
  constexpr int to8 = 3;
  constexpr int to10 = 4;
  constexpr int to12 = 5;
  constexpr int to14 = 6;
  constexpr int to16 = 7;
  constexpr int to18 = 8;
  constexpr int to20 = 9;
  constexpr int to1 = 10;
  double *dx = (double *)malloc(sizeof(double) * totalAN);
  double *dy = (double *)malloc(sizeof(double) * totalAN);
  double *dz = (double *)malloc(sizeof(double) * totalAN);

  // vectors should help with the memory leaks
  std::array<std::vector<double>, 10> Xpow = {
      std::vector<double>(totalAN),                    // x2
      std::vector<double>((lMax > 3) ? totalAN : 0),   // x4
      std::vector<double>((lMax > 5) ? totalAN : 0),   // x6
      std::vector<double>((lMax > 7) ? totalAN : 0),   // x8
      std::vector<double>((lMax > 9) ? totalAN : 0),   // x10
      std::vector<double>((lMax > 11) ? totalAN : 0),  // x12
      std::vector<double>((lMax > 13) ? totalAN : 0),  // x14
      std::vector<double>((lMax > 15) ? totalAN : 0),  // x16
      std::vector<double>((lMax > 17) ? totalAN : 0),  // x18
      std::vector<double>((lMax > 19) ? totalAN : 0)}; // x20
  std::array<std::vector<double>, 10> Ypow = {
      std::vector<double>(totalAN),                    // y2
      std::vector<double>((lMax > 3) ? totalAN : 0),   // y4
      std::vector<double>((lMax > 5) ? totalAN : 0),   // y6
      std::vector<double>((lMax > 7) ? totalAN : 0),   // y8
      std::vector<double>((lMax > 9) ? totalAN : 0),   // y10
      std::vector<double>((lMax > 11) ? totalAN : 0),  // y12
      std::vector<double>((lMax > 13) ? totalAN : 0),  // y14
      std::vector<double>((lMax > 15) ? totalAN : 0),  // y16
      std::vector<double>((lMax > 17) ? totalAN : 0),  // y18
      std::vector<double>((lMax > 19) ? totalAN : 0)}; // y20
  std::array<std::vector<double>, 10> Zpow = {
      std::vector<double>(totalAN),                    // z2
      std::vector<double>((lMax > 3) ? totalAN : 0),   // z4
      std::vector<double>((lMax > 5) ? totalAN : 0),   // z6
      std::vector<double>((lMax > 7) ? totalAN : 0),   // z8
      std::vector<double>((lMax > 9) ? totalAN : 0),   // z10
      std::vector<double>((lMax > 11) ? totalAN : 0),  // z12
      std::vector<double>((lMax > 13) ? totalAN : 0),  // z14
      std::vector<double>((lMax > 15) ? totalAN : 0),  // z16
      std::vector<double>((lMax > 17) ? totalAN : 0),  // z18
      std::vector<double>((lMax > 19) ? totalAN : 0)}; // z20

  std::array<std::vector<double>, 11> Rpow = {
      std::vector<double>(totalAN),                   // r2
      std::vector<double>((lMax > 3) ? totalAN : 0),  // r4
      std::vector<double>((lMax > 5) ? totalAN : 0),  // r6
      std::vector<double>((lMax > 7) ? totalAN : 0),  // r8
      std::vector<double>((lMax > 9) ? totalAN : 0),  // r10
      std::vector<double>((lMax > 11) ? totalAN : 0), // r12
      std::vector<double>((lMax > 13) ? totalAN : 0), // r14
      std::vector<double>((lMax > 15) ? totalAN : 0), // r16
      std::vector<double>((lMax > 17) ? totalAN : 0), // r18
      std::vector<double>((lMax > 19) ? totalAN : 0), // r20
      std::vector<double>(totalAN)};                  // r1

  double *exes = (double *)malloc(sizeof(double) * totalAN);
  // -4 -> no need for l=0, l=1.
  double *preCoef = (double *)malloc(((lMax + 1) * (lMax + 1) - 4) *
                                     sizeof(double) * totalAN);
  double *prCofDX;
  double *prCofDY;
  double *prCofDZ;
  if (return_derivatives) {
    prCofDX = (double *)malloc(((lMax + 1) * (lMax + 1) - 4) * sizeof(double) *
                               totalAN);
    prCofDY = (double *)malloc(((lMax + 1) * (lMax + 1) - 4) * sizeof(double) *
                               totalAN);
    prCofDZ = (double *)malloc(((lMax + 1) * (lMax + 1) - 4) * sizeof(double) *
                               totalAN);
  }

  double *bOa = (double *)malloc((lMax + 1) * nMax2 * sizeof(double));
  double *aOa = (double *)malloc((lMax + 1) * nMax * sizeof(double));

  // Initialize temporary numpy array for storing the coefficients and the
  // averaged coefficients if inner averaging was requested.
  const int n_coeffs = nSpecies * nMax * (lMax + 1) * (lMax + 1);
  double *cnnd_raw = new double[nCenters * n_coeffs]();
  py::array_t<double> cnnd({nCenters, nSpecies, nMax, (lMax + 1) * (lMax + 1)},
                           cnnd_raw);

  auto cnnd_u = cnnd.unchecked<4>();
  auto cdevX_u = cdevX.unchecked<5>();
  auto cdevY_u = cdevY.unchecked<5>();
  auto cdevZ_u = cdevZ.unchecked<5>();

  auto cnnd_mu = cnnd.mutable_unchecked<4>();
  auto cdevX_mu = cdevX.mutable_unchecked<5>();
  auto cdevY_mu = cdevY.mutable_unchecked<5>();
  auto cdevZ_mu = cdevZ.mutable_unchecked<5>();

  // Initialize binning for atoms and centers
  CellList cell_list_centers(centers, rCut + cutoffPadding);

  // Create a mapping between an atomic index and its internal index in the
  // output. The list of species is already ordered.
  map<int, int> ZIndexMap;
  for (int i = 0; i < species.size(); ++i) {
    ZIndexMap[species(i)] = i;
  }

  getAlphaBetaD(aOa, bOa, alphas, betas, nMax, lMax, oOeta, oOeta3O2);

  // Loop through the centers
  for (int i = 0; i < nCenters; i++) {

    // Get all neighbouring atoms for the center i
    double ix = centers_u(i, 0);
    double iy = centers_u(i, 1);
    double iz = centers_u(i, 2);
    CellListResult result =
        cell_list_atoms.getNeighboursForPosition(ix, iy, iz);

    // Sort the neighbours by type
    map<int, vector<int>> atomicTypeMap;
    for (const int &idx : result.indices) {
      int Z = atomicNumbers(idx);
      atomicTypeMap[Z].push_back(idx);
    };

    // Loop through neighbours sorted by type
    for (const auto &ZIndexPair : atomicTypeMap) {

      // j is the internal index for this atomic number
      int j = ZIndexMap[ZIndexPair.first];
      int n_neighbours = ZIndexPair.second.size();

      // Save the neighbour distances into the arrays dx, dy and dz
      getDeltasD(dx, dy, dz, positions, ix, iy, iz, ZIndexPair.second);
      getRsZsD(dx, Xpow[to2].data(), Xpow[to4].data(), Xpow[to6].data(),
               Xpow[to8].data(), Xpow[to10].data(), Xpow[to12].data(),
               Xpow[to14].data(), Xpow[to16].data(), Xpow[to18].data(), dy,
               Ypow[to2].data(), Ypow[to4].data(), Ypow[to6].data(),
               Ypow[to8].data(), Ypow[to10].data(), Ypow[to12].data(),
               Ypow[to14].data(), Ypow[to16].data(), Ypow[to18].data(), dz,
               Rpow[to2].data(), Rpow[to4].data(), Rpow[to6].data(),
               Rpow[to8].data(), Rpow[to10].data(), Rpow[to12].data(),
               Rpow[to14].data(), Rpow[to16].data(), Rpow[to18].data(),
               Zpow[to2].data(), Zpow[to4].data(), Zpow[to6].data(),
               Zpow[to8].data(), Zpow[to10].data(), Zpow[to12].data(),
               Zpow[to14].data(), Zpow[to16].data(), Zpow[to18].data(),
               Rpow[to20].data(), Xpow[to20].data(), Ypow[to20].data(),
               Zpow[to20].data(), n_neighbours, lMax);
      getWeights(n_neighbours, Rpow[to1].data(), Rpow[to2].data(), true,
                 weighting, weights);
      getCfactorsD(preCoef, prCofDX, prCofDY, prCofDZ, n_neighbours, dx,
                   Xpow[to2].data(), Xpow[to4].data(), Xpow[to6].data(),
                   Xpow[to8].data(), Xpow[to10].data(), Xpow[to12].data(),
                   Xpow[to14].data(), Xpow[to16].data(), Xpow[to18].data(), dy,
                   Ypow[to2].data(), Ypow[to4].data(), Ypow[to6].data(),
                   Ypow[to8].data(), Ypow[to10].data(), Ypow[to12].data(),
                   Ypow[to14].data(), Ypow[to16].data(), Ypow[to18].data(), dz,
                   Zpow[to2].data(), Zpow[to4].data(), Zpow[to6].data(),
                   Zpow[to8].data(), Zpow[to10].data(), Zpow[to12].data(),
                   Zpow[to14].data(), Zpow[to16].data(), Zpow[to18].data(),
                   Rpow[to2].data(), Rpow[to4].data(), Rpow[to6].data(),
                   Rpow[to8].data(), Rpow[to10].data(), Rpow[to12].data(),
                   Rpow[to14].data(), Rpow[to16].data(), Rpow[to18].data(),
                   Rpow[to20].data(), Xpow[to20].data(), Ypow[to20].data(),
                   Zpow[to20].data(), totalAN, lMax, return_derivatives);
      getCD(cdevX_mu, cdevY_mu, cdevZ_mu, prCofDX, prCofDY, prCofDZ, cnnd_mu,
            preCoef, dx, dy, dz, Rpow[to2].data(), weights, bOa, aOa, exes,
            totalAN, n_neighbours, nMax, nSpecies, lMax, i, j,
            ZIndexPair.second, return_derivatives);
    }
  }
  free(dx);
  /*free(x2);
  free(x4);
  free(x6);
  free(x8);
  free(x10);
  free(x12);
  free(x14);
  free(x16);
  free(x18);*/
  free(dy);
  /*
  free(y2);
  free(y4);
  free(y6);
  free(y8);
  free(y10);
  free(y12);
  free(y14);
  free(y16);
  free(y18);
  */
  free(dz);
  /*
  free(z2);
  free(z4);
  free(z6);
  free(z8);
  free(z10);
  free(z12);
  free(z14);
  free(z16);
  free(z18);
  free(r2);
  free(r4);
  free(r6);
  free(r8);
  free(r10);
  free(r12);
  free(r14);
  free(r16);
  free(r18);
  free(r20);
  free(x20);
  free(y20);
  free(z20);*/
  free(exes);
  free(preCoef);
  free(bOa);
  free(aOa);
  free(cnnd_raw);

  // Calculate the descriptor value if requested
  if (return_descriptor) {
    auto descriptor_mu = descriptor.mutable_unchecked<2>();

    // If inner averaging is requested, average the coefficients over the
    // centers (axis 0) before calculating the power spectrum.
    if (average == "inner") {
      double *cnnd_ave_raw = new double[n_coeffs]();
      py::array_t<double> cnnd_ave = py::array_t<double>(
          {1, nSpecies, nMax, (lMax + 1) * (lMax + 1)}, cnnd_ave_raw);

      auto cnnd_ave_mu = cnnd_ave.mutable_unchecked<4>();
      auto cnnd_ave_u = cnnd_ave.unchecked<4>();
      for (int i = 0; i < nCenters; i++) {
        for (int j = 0; j < nSpecies; j++) {
          for (int k = 0; k < nMax; k++) {
            for (int l = 0; l < (lMax + 1) * (lMax + 1); l++) {
              cnnd_ave_mu(0, j, k, l) += cnnd_u(i, j, k, l);
            }
          }
        }
      }
      for (int j = 0; j < nSpecies; j++) {
        for (int k = 0; k < nMax; k++) {
          for (int l = 0; l < (lMax + 1) * (lMax + 1); l++) {
            cnnd_ave_mu(0, j, k, l) =
                cnnd_ave_mu(0, j, k, l) / (double)nCenters;
          }
        }
      }
      getPD(descriptor_mu, cnnd_ave_u, nMax, nSpecies, 1, lMax, crossover);
      delete[] cnnd_ave_raw;
      // If outer averaging is requested, average the power spectrum across the
      // centers.
    } else if (average == "outer") {
      // We allocate the memory and give array_t a pointer to it. This way
      // the memory is owned and freed by C++.
      double *ps_temp_raw = new double[nCenters * nFeatures]();
      py::array_t<double> ps_temp({nCenters, nFeatures}, ps_temp_raw);
      auto ps_temp_mu = ps_temp.mutable_unchecked<2>();
      getPD(ps_temp_mu, cnnd_u, nMax, nSpecies, nCenters, lMax, crossover);
      for (int i = 0; i < nCenters; i++) {
        for (int j = 0; j < nFeatures; j++) {
          descriptor_mu(0, j) += ps_temp_mu(i, j);
        }
      }
      for (int j = 0; j < nFeatures; j++) {
        descriptor_mu(0, j) = descriptor_mu(0, j) / (double)nCenters;
      }
      delete[] ps_temp_raw;
      // Regular power spectrum without averaging
    } else {
      getPD(descriptor_mu, cnnd_u, nMax, nSpecies, nCenters, lMax, crossover);
    }
  }

  // Calculate the derivatives
  if (return_derivatives) {
    getPDev(derivatives_mu, positions_u, indices_u, cell_list_centers, cdevX_u,
            cdevY_u, cdevZ_u, cnnd_u, nMax, nSpecies, nCenters, lMax,
            crossover);
  }

  return;
}
