// lambda2.hpp -- lambda2 vortex criterion (Jeong & Hussain 1995) for NekRS v26.
//
// lambda2 = middle eigenvalue of A = S^2 + Omega^2, S/Omega = sym/antisym part of grad(u).
// Same definition as Nek5000 core/postpro.f `lambda2`, but without its post-filter.
// S and Omega come from nrs->strainRotationRate (as for nrs->Qcriterion); a pointwise
// kernel solves the 3x3 symmetric eigenproblem in closed form.
//
// Usage in the .udf (header is included TWICE so that the kernel reaches okl.cpp):
//
//   #ifdef __okl__
//   #include "lambda2.hpp"
//   #endif
//   #include "lambda2.hpp"
//
//   deviceMemory<dfloat> o_lambda2;
//   // UDF_Setup:        o_lambda2.resize(nrs->meshV->Nlocal);
//   //                   nrs->addUserCheckpointField("scalar00", std::vector<deviceMemory<dfloat>>{o_lambda2});
//   // UDF_ExecuteStep:  if (nrs->checkpointStep) lambda2(o_lambda2);
//
// Interface mirrors nrs->Qcriterion (smooth = dssum + inverse-mass average of the gradient;
// smooth = false: element-local gradient, closer to Nek5000's unsmoothed lambda2).

#ifndef NEKRS_LAMBDA2_HPP_
#define NEKRS_LAMBDA2_HPP_

#ifndef __okl__
//#include "nrs.hpp"
//#include "platform.hpp"

inline void lambda2(dlong offset, const occa::memory &o_U, occa::memory &o_L2, bool smooth = true)
{
  // o_SO: S11 S12 S13 S22 S23 S33 | w_z w_y w_x (each field has stride offset)
  auto o_SO = nrs->strainRotationRate(offset, o_U, smooth);
  lambda2Kernel(nrs->meshV->Nlocal, offset, o_SO, o_L2);
}

inline void lambda2(occa::memory &o_L2, bool smooth = true)
{
  lambda2(nrs->fluid->fieldOffset, nrs->fluid->o_U, o_L2, smooth);
}

inline occa::memory lambda2(dlong offset, const occa::memory &o_U, bool smooth = true)
{
  auto o_L2 = platform->deviceMemoryPool.reserve<dfloat>(nrs->meshV->Nlocal);
  lambda2(offset, o_U, o_L2, smooth);
  return o_L2;
}

inline occa::memory lambda2(bool smooth = true)
{
  return lambda2(nrs->fluid->fieldOffset, nrs->fluid->o_U, smooth);
}
#endif

#ifdef __okl__
// Middle eigenvalue of A = S^2 + Omega^2.
// Input: strain entries and w = (SO[8], SO[7], SO[6]) = half the vorticity (wx,wy,wz), so that
//   Omega^2 = w w^T - |w|^2 I   (independent of the sign convention of Omega).
dfloat lambda2Eig(const dfloat s11,
                  const dfloat s12,
                  const dfloat s13,
                  const dfloat s22,
                  const dfloat s23,
                  const dfloat s33,
                  const dfloat wx,
                  const dfloat wy,
                  const dfloat wz)
{
  const dfloat ww = wx * wx + wy * wy + wz * wz;

  const dfloat a11 = s11 * s11 + s12 * s12 + s13 * s13 + wx * wx - ww;
  const dfloat a22 = s12 * s12 + s22 * s22 + s23 * s23 + wy * wy - ww;
  const dfloat a33 = s13 * s13 + s23 * s23 + s33 * s33 + wz * wz - ww;
  const dfloat a12 = s11 * s12 + s12 * s22 + s13 * s23 + wx * wy;
  const dfloat a13 = s11 * s13 + s12 * s23 + s13 * s33 + wx * wz;
  const dfloat a23 = s12 * s13 + s22 * s23 + s23 * s33 + wy * wz;

  // shifted/scaled trigonometric solution of the characteristic cubic (A symmetric => real roots)
  const dfloat q = (a11 + a22 + a33) / 3.0;
  const dfloat d1 = a11 - q;
  const dfloat d2 = a22 - q;
  const dfloat d3 = a33 - q;
  const dfloat p2 = d1 * d1 + d2 * d2 + d3 * d3 + 2.0 * (a12 * a12 + a13 * a13 + a23 * a23);
  const dfloat p = sqrt(p2 / 6.0);

  const dfloat scale =
      fabs(a11) + fabs(a22) + fabs(a33) + fabs(a12) + fabs(a13) + fabs(a23);

  // A ~ q I (zero gradient, A = 0, isotropic): all three eigenvalues equal q
  if (p <= 1e-14 * scale) {
    return q;
  }

  const dfloat b11 = d1 / p;
  const dfloat b22 = d2 / p;
  const dfloat b33 = d3 / p;
  const dfloat b12 = a12 / p;
  const dfloat b13 = a13 / p;
  const dfloat b23 = a23 / p;

  const dfloat detB = b11 * (b22 * b33 - b23 * b23) - b12 * (b12 * b33 - b23 * b13) +
                            b13 * (b12 * b23 - b22 * b13);

  // keep the sign of r (no abs) and clamp: acos stays in its domain
  dfloat r = 0.5 * detB;
  r = (r < -1.0) ? -1.0 : r;
  r = (r > 1.0) ? 1.0 : r;

  const dfloat phi = acos(r) / 3.0;
  const dfloat lmax = q + 2.0 * p * cos(phi);
  const dfloat lmin = q + 2.0 * p * cos(phi + 2.0943951023931953); // + 2 pi/3

  return 3.0 * q - lmax - lmin; // trace(A) = 3q
}

@kernel void lambda2Kernel(const dlong N,
                           const dlong offset,
                           @ restrict const dfloat *SO,
                           @ restrict dfloat *L2)
{
  for (dlong n = 0; n < N; ++n; @tile(256, @outer, @inner)) {
    if (n < N) {
      L2[n] = lambda2Eig(SO[n + 0 * offset],
                         SO[n + 1 * offset],
                         SO[n + 2 * offset],
                         SO[n + 3 * offset],
                         SO[n + 4 * offset],
                         SO[n + 5 * offset],
                         SO[n + 8 * offset],
                         SO[n + 7 * offset],
                         SO[n + 6 * offset]);
    }
  }
}
#endif // __okl__

#endif // NEKRS_LAMBDA2_HPP_
