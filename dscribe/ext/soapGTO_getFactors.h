#ifndef SOAPGTO_GETFACTORS_H
#define SOAPGTO_GETFACTORS_H
void getCfactorsD(
    double *preCoef, double *prCofDX, double *prCofDY, double *prCofDZ,
    const int Asize, const double *x, const double *x2, const double *x4,
    const double *x6, const double *x8, const double *x10, const double *x12,
    const double *x14, const double *x16, const double *x18, const double *y,
    const double *y2, const double *y4, const double *y6, const double *y8,
    const double *y10, const double *y12, const double *y14, const double *y16,
    const double *y18, const double *z, const double *z2, const double *z4,
    const double *z6, const double *z8, const double *z10, const double *z12,
    const double *z14, const double *z16, const double *z18, const double *r2,
    const double *r4, const double *r6, const double *r8, const double *r10,
    const double *r12, const double *r14, const double *r16, const double *r18,
    const double *r20, const double *x20, const double *y20, const double *z20,
    const int totalAN, const int lMax, const bool return_derivatives);

#endif // SOAPGTO_GETFACTORS_H