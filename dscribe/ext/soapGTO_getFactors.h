#ifndef SOAPGTO_GETFACTOR_H
#define SOAPGTO_GETFACTOR_H
void getCfactorsD(double *preCoef, double *prCofDX, double *prCofDY,
                  double *prCofDZ, int Asize, double *x, double *x2, double *x4,
                  double *x6, double *x8, double *x10, double *x12, double *x14,
                  double *x16, double *x18, double *y, double *y2, double *y4,
                  double *y6, double *y8, double *y10, double *y12, double *y14,
                  double *y16, double *y18, double *z, double *z2, double *z4,
                  double *z6, double *z8, double *z10, double *z12, double *z14,
                  double *z16, double *z18, double *r2, double *r4, double *r6,
                  double *r8, double *r10, double *r12, double *r14,
                  double *r16, double *r18, double *r20, double *x20,
                  double *y20, double *z20, int totalAN, int lMax,
                  const bool return_derivatives);

#endif // SOAPGTO_GETFACTOR_H