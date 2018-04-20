#include "float.h"

#define SIGN(x) ((x) >= 0 ? 1 : -1)

#ifdef USE_MIXED_PRECISION
#   define NaN nan("0")
#   define EPSILON DBL_EPSILON
#   define MAXREAL DBL_MAX
#   define MINREAL DBL_MIN
#   define errtol 0.001
#else
#   define NaN nanf("0")
#   define EPSILON FLT_EPSILON
#   define MAXREAL FLT_MAX
#   define MINREAL FLT_MIN
#   define errtol 0.03
#endif

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

inline __device__ void jacobi(mixed u, mixed m, mixed& sn, mixed& cn, mixed& dn)
{
  if (fabs(m) > 1.0)
    sn = cn = dn = NaN;
  else if(fabs(m) < 2.0*EPSILON) {
    sn = sin(u);
    cn = cos(u);
    dn = 1.0;
  }
  else if (fabs(m - 1.0) < 2.0*EPSILON) {
    sn = tanh(u);
    cn = 1.0/cosh(u);
    dn = cn;
  }
  else {
    const int N = 16;
    mixed mu[16];
    mixed nu[16];
    mixed c[16];
    mixed d[16];
    mixed sin_umu, cos_umu, t, r;
    int n = 0;
    mu[0] = 1.0;
    nu[0] = sqrt(1.0 - m);
    while ( fabs(mu[n] - nu[n]) > 4.0 * EPSILON * fabs(mu[n]+nu[n])) {
      mu[n+1] = 0.5 * (mu[n] + nu[n]);
      nu[n+1] = sqrt(mu[n] * nu[n]);
      ++n;
      if (n >= N - 1) {
          sn = cn = dn = NaN;
          return;
      }
    }
    sin_umu = sin(u * mu[n]);
    cos_umu = cos(u * mu[n]);
    /* Since sin(u*mu(n)) can be zero we switch to computing sn(K-u),
       cn(K-u), dn(K-u) when |sin| < |cos| */
    if (fabs(sin_umu) < fabs(cos_umu)) {
        t = sin_umu / cos_umu;
        c[n] = mu[n] * t;
        d[n] = 1.0;
        while(n > 0) {
          n--;
          c[n] = d[n+1] * c[n+1];
          r = (c[n+1] * c[n+1]) / mu[n+1];
          d[n] = (r + nu[n]) / (r + mu[n]);
        }
        dn = sqrt(1.0-m) / d[n];
        cn = (dn) * SIGN(cos_umu) / hypot(1.0, c[n]);
        sn = (cn) * c[n] /sqrt(1.0-m);
    }
    else {
        t = cos_umu / sin_umu;
        c[n] = mu[n] * t;
        d[n] = 1.0;
        while (n > 0) {
          --n;
          c[n] = d[n+1] * c[n+1];
          r = (c[n+1] * c[n+1]) / mu[n+1];
          d[n] = (r + nu[n]) / (r + mu[n]);
        }
        dn = d[n];
        sn = SIGN(sin_umu) / hypot(1.0, c[n]);
        cn = c[n] * (sn);
    }
  }
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed carlsonRC(mixed x, mixed y)
{
  const mixed lolim = 5.0 * MINREAL;
  const mixed uplim = 0.2 * MAXREAL;
  const int nmax = 10000;
  if (x < 0.0 || y < 0.0 || x + y < lolim)
    return NaN;
  else if (max(x, y) < uplim) {
    const mixed c1 = 1.0 / 7.0;
    const mixed c2 = 9.0 / 22.0;
    mixed xn = x;
    mixed yn = y;
    mixed mu, sn, lamda, s;
    int n = 0;
    while (1) {
      mu = (xn + yn + yn) / 3.0;
      sn = (yn + mu) / mu - 2.0;
      if (fabs(sn) < errtol) break;
      lamda = 2.0 * sqrt(xn) * sqrt(yn) + yn;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
      n++;
      if (n==nmax) return NaN;
    }
    s = sn * sn * (0.3 + sn * (c1 + sn * (0.375 + sn * c2)));
    return (1.0 + s) / sqrt(mu);
  }
  else
    return NaN;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed carlsonRF(mixed x, mixed y, mixed z)
{
  const mixed lolim = 5.0 * MINREAL;
  const mixed uplim = 0.2 * MAXREAL;
  const int nmax = 10000;
  if (x < 0.0 || y < 0.0 || z < 0.0)
    return NaN;
  else if (x+y < lolim || x+z < lolim || y+z < lolim)
    return NaN;
  else if (max(x,max(y,z)) < uplim) {
    const mixed c1 = 1.0 / 24.0;
    const mixed c2 = 3.0 / 44.0;
    const mixed c3 = 1.0 / 14.0;
    mixed xn = x;
    mixed yn = y;
    mixed zn = z;
    mixed mu, xndev, yndev, zndev, e2, e3, s;
    int n = 0;
    while (1) {
      mixed epslon, lamda;
      mixed xnroot, ynroot, znroot;
      mu = (xn + yn + zn) / 3.0;
      xndev = 2.0 - (mu + xn) / mu;
      yndev = 2.0 - (mu + yn) / mu;
      zndev = 2.0 - (mu + zn) / mu;
      epslon = max(fabs(xndev), max(fabs(yndev), fabs(zndev)));
      if (epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
      zn = (zn + lamda) * 0.25;
      n++;
      if (n == nmax) return NaN;
    }
    e2 = xndev * yndev - zndev * zndev;
    e3 = xndev * yndev * zndev;
    s = 1.0 + (c1 * e2 - 0.1 - c2 * e3) * e2 + c3 * e3;
    return s / sqrt(mu);
  }
  else
    return NaN;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed carlsonRJ(mixed x, mixed y, mixed z, mixed p)
{
  const mixed lolim =     pow(5.0*MINREAL, 1.0/3.0);
  const mixed uplim = 0.3*pow(0.2*MAXREAL, 1.0/3.0);
  const int nmax = 10000;
  if (x < 0.0 || y < 0.0 || z < 0.0)
    return NaN;
  else if (x + y < lolim || x + z < lolim || y + z < lolim || p < lolim)
    return NaN;
  else if (max(max(x,y),max(z,p)) < uplim) {
    const mixed c1 = 3.0 / 14.0;
    const mixed c2 = 1.0 /  3.0;
    const mixed c3 = 3.0 / 22.0;
    const mixed c4 = 3.0 / 26.0;
    mixed xn = x;
    mixed yn = y;
    mixed zn = z;
    mixed pn = p;
    mixed sigma = 0.0;
    mixed power4 = 1.0;
    mixed mu, xndev, yndev, zndev, pndev;
    mixed ea, eb, ec, e2, e3, s1, s2, s3;
    int n = 0;
    while(1) {
      mixed xnroot, ynroot, znroot;
      mixed lamda, alfa, beta;
      mixed epslon;
      mu = (xn + yn + zn + pn + pn) * 0.2;
      xndev = (mu - xn) / mu;
      yndev = (mu - yn) / mu;
      zndev = (mu - zn) / mu;
      pndev = (mu - pn) / mu;
      epslon = max(max(fabs(xndev), fabs(yndev)), max(fabs(zndev), fabs(pndev)));
      if (epslon < errtol) break;
      xnroot = sqrt(xn);
      ynroot = sqrt(yn);
      znroot = sqrt(zn);
      lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
      alfa = pn * (xnroot + ynroot + znroot) + xnroot * ynroot * znroot;
      alfa = alfa * alfa;
      beta = pn * (pn + lamda) * (pn + lamda);
      mixed rc = carlsonRC(alfa, beta);
      if (isnan(rc)) return NaN;
      sigma  += power4 * rc;
      power4 *= 0.25;
      xn = (xn + lamda) * 0.25;
      yn = (yn + lamda) * 0.25;
      zn = (zn + lamda) * 0.25;
      pn = (pn + lamda) * 0.25;
      n++;
      if (n == nmax) return NaN;
    }
    ea = xndev * (yndev + zndev) + yndev * zndev;
    eb = xndev * yndev * zndev;
    ec = pndev * pndev;
    e2 = ea - 3.0 * ec;
    e3 = eb + 2.0 * pndev * (ea - ec);
    s1 = 1.0 + e2 * (- c1 + 0.75 * c3 * e2 - 1.5 * c4 * e3);
    s2 = eb * (0.5 * c2 + pndev * (- c3 - c3 + pndev * c4));
    s3 = pndev * ea * (c2 - pndev * c3) - c2 * pndev * ec;
    return 3.0 * sigma + power4 * (s1 + s2 + s3) / (mu * sqrt(mu));
  }
  else
    return NaN;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed Omega(mixed x, mixed n, mixed m) {
    mixed x2 = x*x;
    return (-1.0/3.0)*n*x*x2*carlsonRJ(1.0 - x2, 1.0 - m*x2, 1.0, 1.0 + n*x2);
}
