// fs-vkplate's time step (crates/fs-vkplate/src/lib.rs, Instrument::step) as ONE persistent CUDA kernel.
//
// The host prepends `struct P` (every pointer and size, generated from cuda.rs's one field list), and the
// compile-time sizes GB (simulations per pass over an Airy inverse) and AIRY_F32. All blocks stay
// resident for the whole launch and meet at a grid barrier between dependent phases:
//
//   P1  nodes: L w, linear forces, s1 term, per-body sums;  core nodes: Airy right-hand side q, curvatures;
//       link forces
//   P2  (plates only) every warp: F = A^-1 q, the dense inverse read once for up to GB simulations
//   P3  nodes: links, cavity pressure, radiation damping, drives (each block computes the scalars it needs;
//       one item advances each cavity port), nonlinear gradient g gathered from F;  core: sum F q
//   P4  nodes: source curvature-rate viscosity, SAV auxiliary variable, explicit update wn, Woodbury sums
//   P5  nodes away from contacts: Woodbury correction;  one block per simulation: the same correction on
//       contact nodes and strikes; then contacts in node-disjoint groups and frictions in coloured groups,
//       each group over the whole GPU when it is large (a barrier before it), else on its simulation's block
//   P6  nodes: energy-budget sums;  microphone taps in balanced segments
//   P7  block 0: auxiliary-variable advance, segment sums and port taps into the rings, microphone outputs
//
// A step is latency-bound (a few thousand nodes, 5-6 barriers), so the layout avoids dependent load chains:
// stencils are fixed-width column-major tables (every load of a row issues at once), per-body sums load all
// their partials before any butterfly, and small scalars are recomputed where they are needed instead of
// costing a barrier. Arithmetic follows the CPU code's order wherever no reduction is involved (with
// fused multiply-adds); sums use fixed trees (a butterfly gives every lane, and every block, the same bits),
// so a run is deterministic but rounds differently from the CPU's sequential sums. The CPU step stays the
// reference (render-gpu --check compares the two).

typedef unsigned int u32;

// Same variable-mass edge flux as Csr::mass_flux_correction. The shared
// Laplacian includes fixed-boundary diagonal terms; those remain unchanged.
__device__ __forceinline__ double mass_flux_correction(const P& p, int i, ll o,
                                                       const double* w, const double* wp) {
  const double xi = w[o + i] - wp[o + i], mi = __ldg(&p.mass[i]);
  double correction = 0.0;
  for (int q = 0; q < p.KL; q++) {
    const int j = __ldg(&p.lapE_col[(ll)q * p.NT + i]);
    correction += 0.5 * __ldg(&p.lapE_val[(ll)q * p.NT + i]) * (__ldg(&p.mass[j]) - mi)
                  * ((w[o + j] - wp[o + j]) - xi);
  }
  return correction;
}

#define TPB 256
#ifndef SPIN_NS
#define SPIN_NS 32
#endif

__device__ __forceinline__ void gsync(u32* bar, u32 nb) {
  __syncthreads();
  if (threadIdx.x == 0) {
    volatile u32* gen = bar + 32;
    u32 g = *gen;
    __threadfence();
    if (atomicAdd(bar, 1u) == nb - 1) { atomicExch(bar, 0u); __threadfence(); atomicAdd(bar + 32, 1u); }
    else {
      while (*gen == g) {
#if __CUDA_ARCH__ >= 700
        __nanosleep(SPIN_NS);
#endif
      }
    }
    __threadfence();
  }
  __syncthreads();
}

// Fixed-tree block sum of NV values per thread (a shuffle tree per warp, then the warps' sums in warp
// order); thread 0 gets the sums in v.
template <int NV>
__device__ __forceinline__ void block_sum(double (&v)[NV], double (*sh)[TPB]) {
  const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
#pragma unroll
  for (int k = 0; k < NV; k++) {
    double a = v[k];
    for (int off = 16; off; off >>= 1) a = a + __shfl_down_sync(0xffffffffu, a, off);
    if (lane == 0) sh[k][warp] = a;
  }
  __syncthreads();
  if (warp == 0) {
#pragma unroll
    for (int k = 0; k < NV; k++) {
      double a = lane < TPB / 32 ? sh[k][lane] : 0.0;
      for (int off = 4; off; off >>= 1) a = a + __shfl_down_sync(0xffffffffu, a, off);
      v[k] = a;
    }
  }
  __syncthreads();
}

__device__ __forceinline__ double bfly(double a) {
  for (int off = 16; off; off >>= 1) a = a + __shfl_xor_sync(0xffffffffu, a, off);
  return a;
}

// Sums of partials k0 .. k0+NK-1 over a body's work items [i0, i1) for simulation s, by a whole warp:
// every lane loads its items' partials first (lane, lane + 32, ...), then one butterfly per sum. Every lane,
// and every block, gets the same bits.
template <int NK>
__device__ __forceinline__ void wsums(double (&out)[NK], const double* part, int stride, int k0, int i0, int i1,
                                      int s, ll NI) {
  const int lane = threadIdx.x & 31;
  double a[NK];
#pragma unroll
  for (int k = 0; k < NK; k++) a[k] = 0.0;
  for (int it = i0 + lane; it < i1; it += 32) {
    const double* r = part + ((ll)s * NI + it) * stride + k0;
#pragma unroll
    for (int k = 0; k < NK; k++) a[k] = a[k] + r[k];
  }
#pragma unroll
  for (int k = 0; k < NK; k++) out[k] = bfly(a[k]);
}

// Secant of d_+^(alpha+1)/(alpha+1), matching fs-vkplate power_secant.
// Centered leapfrog uses the previous and next positions, not the current one.
// Bilbao et al., Numerical Modeling of Collisions, section 2.3, equation (11):
// https://arxiv.org/abs/1405.2589. Whole coupled-system passivity is separate.
__device__ __forceinline__ void power_secant(double d, double previous, double alpha,
                                             double& value, double& derivative) {
  if (d <= 0.0 && previous <= 0.0) { value = 0.0; derivative = 0.0; return; }
  if (previous <= 0.0) {
    const double den = d - previous;
    if (alpha == 1.5) {
      const double root = sqrt(d), ratio = d / den;
      value = 0.4 * d * root * ratio;
      derivative = 0.4 * root * ratio * (2.5 - ratio);
      return;
    }
    const double p = pow(d, alpha);
    value = d * p / ((alpha + 1.0) * den);
    derivative = p * (alpha * d - (alpha + 1.0) * previous) / ((alpha + 1.0) * den * den);
    return;
  }
  if (d <= 0.0) {
    const double den = previous - d;
    if (alpha == 1.5) {
      const double root = sqrt(previous), ratio = previous / den;
      value = 0.4 * previous * root * ratio;
      derivative = 0.4 * root * ratio * ratio;
      return;
    }
    const double potential = pow(previous, alpha + 1.0) / (alpha + 1.0);
    value = potential / den;
    derivative = potential / (den * den);
    return;
  }
  if (alpha == 1.5) {
    double a = sqrt(d), b = sqrt(previous);
    const double m = fmax(a, b);
    a = a / m; b = b / m;
    const double a2 = a * a, b2 = b * b, sum = a + b;
    value = 0.4 * m * m * m * (a2*a2 + a2*a*b + a2*b2 + a*b2*b + b2*b2) / sum;
    derivative = 0.2 * m * (3.0*a2*a + 6.0*a2*b + 4.0*a*b2 + 2.0*b2*b) / (sum * sum);
    return;
  }
  const double z = (d - previous) / previous, base = pow(previous, alpha);
  if (fabs(z) * fmax(fabs(alpha), 1.0) < 1e-3) {
    const double a2 = alpha * (alpha - 1.0), a3 = a2 * (alpha - 2.0);
    const double a4 = a3 * (alpha - 3.0), a5 = a4 * (alpha - 4.0);
    value = base * (1.0 + z * (0.5*alpha + z * (a2/6.0 + z * (a3/24.0 + z * (a4/120.0 + z*a5/720.0)))));
    derivative = base / previous * (0.5*alpha + z * (a2/3.0 + z * (a3/8.0 + z * (a4/30.0 + z*a5/144.0))));
  } else {
    const double delta = expm1((alpha + 1.0) * log1p(z));
    value = base * delta / ((alpha + 1.0) * z);
    derivative = base / previous * ((alpha + 1.0) * z * pow(1.0 + z, alpha) - delta) / ((alpha + 1.0) * z * z);
  }
}

__device__ __forceinline__ void contact_phi(double d, double previous, double kk, double alpha,
                                            double mu, double k, double& value, double& derivative) {
  double elastic, elastic_derivative;
  power_secant(d, previous, alpha, elastic, elastic_derivative);
  const double damping_rate = mu / (2.0 * k), damp = 1.0 + damping_rate * (d - previous);
  value = kk * elastic * damp;
  if (value <= 0.0) { value = 0.0; derivative = 0.0; }
  else derivative = kk * (elastic_derivative * damp + elastic * damping_rate);
}

// f=phi(d_pred-comp*f). During unloading either endpoint can be separated;
// phi(d_pred) brackets the force even when d_pred/comp would not.
__device__ double contact_force(double d_pred, double comp, double force, double previous, double kk, double alpha,
                                double mu, int fixed_exit, double k) {
  if (isnan(d_pred) || isnan(previous) || (d_pred <= 0.0 && previous <= 0.0)) return 0.0;
  double lo = 0.0, hi, unused;
  contact_phi(d_pred, previous, kk, alpha, mu, k, hi, unused);
  if (hi <= 0.0) return 0.0;
  double x = (force > 0.0 && force < hi) ? force : 0.5 * hi;
  for (int it = 0; it < 60; it++) {
    double v, dv;
    contact_phi(d_pred - comp * x, previous, kk, alpha, mu, k, v, dv);
    double g = x - v;
    // Keep an exact root before the bracket safeguard can step away from it.
    if (g == 0.0) return x;
    if (g > 0.0) hi = x; else lo = x;
    double next = x - g / (1.0 + comp * dv);
    if (!((next > lo && next < hi) || (fixed_exit && next == x))) next = 0.5 * (lo + hi);
    bool done = fabs(next - x) <= 1e-13 * fmax(fabs(next), 1e-30) || hi - lo <= 1e-13 * fmax(hi, 1e-30);
    x = next;
    if (done) break;
  }
  return x;
}

__device__ __forceinline__ double* wbuf(const P& p, ll r) { return r == 0 ? p.W0 : (r == 1 ? p.W1 : p.W2); }

// Stencil weight of core node c's curvature product on its slot (fs-vkplate nonlinear_gradient's `put`).
__device__ __forceinline__ double slot_val(int slot, double xx, double yy, double xy) {
  double gx = -0.5 * xy;
  switch (slot) {
    case 1: case 2: return yy;
    case 0: return -2.0 * yy - 2.0 * xx;
    case 3: case 4: return xx;
    case 5: case 8: return gx;
    default: return -gx;  // 6, 7
  }
}

// Cavity c's pressure for simulation s from this step's head sums (by a whole warp; every lane gets it).
__device__ __forceinline__ double cavity_pressure(const P& p, int c, int s, ll NI, ll par) {
  double dv = 0.0;
  for (int q = __ldg(&p.cav_hptr[c]); q < __ldg(&p.cav_hptr[c + 1]); q++) {
    const int hb = __ldg(&p.cav_head[q]);
    double sw[1];
    wsums<1>(sw, p.p1, 4, 1, __ldg(&p.bitem0[hb]), __ldg(&p.bitem1[hb]), s, NI);
    dv = dv + sw[0] * __ldg(&p.bh[hb]) * __ldg(&p.bh[hb]);
  }
  const double xport = __ldg(&p.cav_port[c]) ? p.port_x[(par * p.B + s) * p.NCAV + c] : 0.0;
  return -p.gas_bulk_modulus / __ldg(&p.cav_vol[c]) * (dv + xport);
}

// The radiation-damping force density -R q of group r (by a whole warp).
__device__ __forceinline__ double rad_value(const P& p, int r, int s, ll NI) {
  double q = 0.0;
  for (int e = __ldg(&p.rad_ptr[r]); e < __ldg(&p.rad_ptr[r + 1]); e++) {
    const int bi = __ldg(&p.rad_body[e]);
    double sd[1];
    wsums<1>(sd, p.p1, 4, 2, __ldg(&p.bitem0[bi]), __ldg(&p.bitem1[bi]), s, NI);
    q = q + sd[0] / p.k * __ldg(&p.bh[bi]) * __ldg(&p.bh[bi]);
  }
  return -__ldg(&p.rad_res[r]) * q;
}

__device__ __forceinline__ double drive_at(const P& p, int d, double t) {
  const int k0 = __ldg(&p.dv_ptr[d]), k1 = __ldg(&p.dv_ptr[d + 1]);
  if (k1 == k0 || t <= __ldg(&p.dv_t[k0])) return k1 > k0 ? __ldg(&p.dv_v[k0]) : 0.0;
  for (int q = k0; q + 1 < k1; q++) {
    if (t <= __ldg(&p.dv_t[q + 1])) {
      const double uu = (t - __ldg(&p.dv_t[q])) / (__ldg(&p.dv_t[q + 1]) - __ldg(&p.dv_t[q]));
      return __ldg(&p.dv_v[q]) + uu * (__ldg(&p.dv_v[q + 1]) - __ldg(&p.dv_v[q]));
    }
  }
  return __ldg(&p.dv_v[k1 - 1]);
}

// Body scalars of the predict step for (sim s, body b); called by a whole warp, identical in every block.
struct Pred { double ir, sc, extra, mu, h2, psi1, rprev1; int scale; };

__device__ Pred predict_scalars(const P& p, int s, int b) {
  Pred r;
  const double k = p.k, k2 = k * k, h = __ldg(&p.bh[b]), h2 = h * h;
  r.h2 = h2;
  const int i0 = __ldg(&p.bitem0[b]), i1 = __ldg(&p.bitem1[b]);
  const bool vk = __ldg(&p.bvk[b]) != 0, ber = __ldg(&p.bber[b]) > 0.0, nl = vk || ber;
  const ll ob = (ll)s * p.NB + b;
  double psi = p.psiA[ob];
  const double rprev = p.rprevA[ob], budget = p.budgetA[ob];
  double s1[4], fq[1], gw[1];
  wsums<4>(s1, p.p1, 4, 0, i0, i1, s, p.NI);  // w.Lw, sum w, sum(w - wp), sum(2w - wp)
  wsums<1>(gw, p.p3, 1, 0, i0, i1, s, p.NI);
  wsums<1>(fq, p.pc3, 1, 0, __ldg(&p.bcitem0[b]), __ldg(&p.bcitem1[b]), s, p.NCI);
  r.scale = 0; r.ir = 0.0; r.rprev1 = rprev;
  if (nl) {
    double v = 0.0;
    if (ber) v = v + 0.25 * __ldg(&p.bber[b]) * s1[0] * s1[0];
    if (vk) v = v + h2 / (2.0 * __ldg(&p.beh[b])) * fq[0];
    v = fmax(v, 0.0);
    double root = sqrt(2.0 * v);
    if (vk) {
      double target = 0.5 * (root + rprev); // zero is a valid previous potential
      if (psi > target) psi = target;
      else psi = fmin(sqrt(psi * psi + 2.0 * budget), target);
      r.rprev1 = root;
    }
    if (root > 0.0) { r.ir = 1.0 / root; r.scale = 1; }
  }
  r.psi1 = psi;
  double gwp = r.scale ? gw[0] * r.ir : 0.0;
  r.sc = psi - 0.25 * gwp;
  r.mu = __ldg(&p.bairm[b]) > 0.0 ? __ldg(&p.bairm[b]) / (__ldg(&p.barea[b]) * __ldg(&p.barea[b]) * k2) : 0.0;
  double vol_prev = r.mu > 0.0 ? s1[3] * h2 : 0.0;
  r.extra = r.mu > 0.0 ? r.mu * h2 * vol_prev : 0.0;
  return r;
}

// Woodbury coefficients (cg, cu) of (sim s, body b) from this step's predict sums (by a whole warp).
__device__ __forceinline__ void woodbury(const P& p, int s, int b, double& cg, double& cu) {
  const double k = p.k;
  const bool scale = p.savf[(ll)s * p.NB + b] > 0.0;
  const double mu = __ldg(&p.bairm[b]) > 0.0 ? __ldg(&p.bairm[b]) / (__ldg(&p.barea[b]) * __ldg(&p.barea[b]) * (k * k)) : 0.0;
  double a[5];
  wsums<5>(a, p.p4, 5, 0, __ldg(&p.bitem0[b]), __ldg(&p.bitem1[b]), s, p.NI);
  const double gy = a[0], gzg = a[1], gzu = a[2], uy = a[3], uzu = a[4];
  cg = 0.0; cu = 0.0;
  if (scale && !(mu > 0.0)) cg = gy / (4.0 + gzg);
  else if (!scale && mu > 0.0) cu = uy / (1.0 / mu + uzu);
  else if (scale && mu > 0.0) {
    const double a11 = 4.0 + gzg, a12 = gzu, a22 = 1.0 / mu + uzu;
    const double det = a11 * a22 - a12 * a12;
    cg = (a22 * gy - a12 * uy) / det; cu = (a11 * uy - a12 * gy) / det;
  }
}

// F = A^-1 q for every plate and simulation: one warp per row, simulations NG at a time (the batch is
// padded to a multiple of GB in Q, QF and AF), so each pass reads the dense inverse once for NG of them.
template <int NG>
__device__ void airy_rows(const P& p, ll gwarp, ll nwarps) {
  const int lane = threadIdx.x & 31;
  const ll B = p.B, NC = p.NC;
  for (ll r = gwarp; r < p.NROWS; r += nwarps) {
    int v = 0; ll row = r;
    while (row >= __ldg(&p.bnc[__ldg(&p.vk_body[v])])) { row -= __ldg(&p.bnc[__ldg(&p.vk_body[v])]); v++; }
    const int b = __ldg(&p.vk_body[v]);
    const ll nc = __ldg(&p.bnc[b]), coff = __ldg(&p.bcoff[b]);
    for (ll s0 = 0; s0 < B; s0 += NG) {
#if AIRY_F32
      const float* M = p.airy32 + __ldg(&p.airy_off[v]) + row * nc;
      float part[NG];
#pragma unroll
      for (int g = 0; g < NG; g++) part[g] = 0.f;
      for (ll j = lane; j < nc; j += 32) {
        const float a = __ldg(&M[j]);
#pragma unroll
        for (int g = 0; g < NG; g++) part[g] = fmaf(a, p.QF[(s0 + g) * NC + coff + j], part[g]);
      }
      double acc[NG];
#pragma unroll
      for (int g = 0; g < NG; g++) acc[g] = (double)part[g];
#else
      const double* M = p.airy64 + __ldg(&p.airy_off[v]) + row * nc;
      double acc[NG];
#pragma unroll
      for (int g = 0; g < NG; g++) acc[g] = 0.0;
      for (ll j = lane; j < nc; j += 32) {
        const double a = __ldg(&M[j]);
#pragma unroll
        for (int g = 0; g < NG; g++) acc[g] = fma(a, p.Q[(s0 + g) * NC + coff + j], acc[g]);
      }
#endif
#pragma unroll
      for (int g = 0; g < NG; g++) {
        double x = acc[g];
        for (int off = 16; off; off >>= 1) x = x + __shfl_down_sync(0xffffffffu, x, off);
        if (lane == 0) p.AF[(s0 + g) * NC + coff + row] = x;
      }
    }
  }
}

// Contact j of simulation s: its force from the predicted gap (implicit, one scalar Newton solve), applied.
__device__ __forceinline__ void solve_contact(const P& p, ll s, int j, const double* Wp, double* Wn, double k) {
  const ll o = s * p.NT, oj = s * p.NCON + j;
  if (!p.ct_active[oj]) return;
  double d_pred = -__ldg(&p.ct_gap[j]);
  double previous = d_pred;
  for (int e = __ldg(&p.ct_a0[j]); e < __ldg(&p.ct_b0[j]); e++) {
    const int nd = __ldg(&p.ct_node[e]);
    const double w = __ldg(&p.ct_w[e]);
    d_pred = d_pred + w * Wn[o + nd];
    previous = previous + w * Wp[o + nd];
  }
  for (int e = __ldg(&p.ct_b0[j]); e < __ldg(&p.ct_end[j]); e++) {
    const int nd = __ldg(&p.ct_node[e]);
    const double w = __ldg(&p.ct_w[e]);
    d_pred = d_pred - w * Wn[o + nd];
    previous = previous - w * Wp[o + nd];
  }
  const double comp = __ldg(&p.ct_comp[j]);
  // Read the previous trajectory directly; cached current compression belongs
  // to a different time level and can also have stale shared-node updates.
  const double f = contact_force(d_pred, comp, p.ct_force[oj], previous, __ldg(&p.ct_kk[j]), __ldg(&p.ct_alpha[j]),
                                 __ldg(&p.ct_mu[j]), __ldg(&p.ct_fx[j]), k);
  if (f > 0.0) {
    for (int e = __ldg(&p.ct_a0[j]); e < __ldg(&p.ct_b0[j]); e++) {
      const int nd = __ldg(&p.ct_node[e]);
      Wn[o + nd] = Wn[o + nd] - __ldg(&p.ct_w[e]) * __ldg(&p.ia[nd]) * f;
    }
    for (int e = __ldg(&p.ct_b0[j]); e < __ldg(&p.ct_end[j]); e++) {
      const int nd = __ldg(&p.ct_node[e]);
      Wn[o + nd] = Wn[o + nd] + __ldg(&p.ct_w[e]) * __ldg(&p.ia[nd]) * f;
    }
  }
  p.ct_peak[oj] = fmax(p.ct_peak[oj], f);
  p.ct_force[oj] = f;
  p.ct_dnow[oj] = d_pred - comp * f;
  // Coulomb friction has zero support without a compressive normal force.
  // This is only an activity witness: never reorder or approximate an active solve.
  if (f > 0.0 && __ldg(&p.ct_has_friction[j])) atomicExch(p.any_friction, 1);
}

// Friction r of simulation s: constrained quadratic solve, as in fs-vkplate.
__device__ __forceinline__ void solve_friction(const P& p, ll s, int r, const double* W, double* Wn) {
  const ll o = s * p.NT;
  const int c = __ldg(&p.fr_contact[r]);
  const double nf = p.ct_active[s * p.NCON + c] ? p.ct_force[s * p.NCON + c] : 0.0;
  if (!(nf > 0.0) || __ldg(&p.fr_mu[r]) <= 0.0) return;
  double dx = 0.0, dy = 0.0;
  for (int e = __ldg(&p.fr_x0[r]); e < __ldg(&p.fr_y0[r]); e++) dx = dx + __ldg(&p.fr_w[e]) * (Wn[o + __ldg(&p.fr_node[e])] - W[o + __ldg(&p.fr_node[e])]);
  for (int e = __ldg(&p.fr_y0[r]); e < __ldg(&p.fr_end[r]); e++) dy = dy + __ldg(&p.fr_w[e]) * (Wn[o + __ldg(&p.fr_node[e])] - W[o + __ldg(&p.fr_node[e])]);
  const double cxx = __ldg(&p.fr_cxx[r]), cyy = __ldg(&p.fr_cyy[r]), cxy = __ldg(&p.fr_cxy[r]);
  const double det = cxx * cyy - cxy * cxy;
  if (!(det > 0.0)) return;
  double fx = -(cyy * dx - cxy * dy) / det, fy = -(cxx * dy - cxy * dx) / det;
  const double mag = sqrt(fx * fx + fy * fy), cap = __ldg(&p.fr_mu[r]) * nf;
  if (cap == 0.0) return;
  if (!(mag <= cap)) {
    // (C + lambda I) f = -d and |f| = cap. Unlike radial clipping,
    // this opposes final slip even when the compliance is anisotropic.
    double lo = 0.0, hi = sqrt(dx * dx + dy * dy) / cap;
    double lambda = isfinite(mag) ? hi * ((mag - cap) / mag) : 0.5 * hi;
    for (int it = 0; it < 64; it++) {
      const double a = cxx + lambda, b = cyy + lambda;
      const double determinant = a * b - cxy * cxy;
      fx = -(b * dx - cxy * dy) / determinant;
      fy = -(a * dy - cxy * dx) / determinant;
      const double norm = sqrt(fx * fx + fy * fy), residual = norm - cap;
      if (fabs(residual) <= 1e-12 * cap) break;
      if (residual > 0.0) lo = lambda; else hi = lambda;
      const double derivative = -(b * fx * fx - 2.0 * cxy * fx * fy + a * fy * fy) / (determinant * norm);
      const double next = lambda - residual / derivative;
      lambda = next > lo && next < hi ? next : 0.5 * (lo + hi);
    }
  }
  for (int e = __ldg(&p.fr_x0[r]); e < __ldg(&p.fr_y0[r]); e++) {
    const int nd = __ldg(&p.fr_node[e]);
    Wn[o + nd] = Wn[o + nd] + __ldg(&p.ia[nd]) * __ldg(&p.fr_w[e]) * fx;
  }
  for (int e = __ldg(&p.fr_y0[r]); e < __ldg(&p.fr_end[r]); e++) {
    const int nd = __ldg(&p.fr_node[e]);
    Wn[o + nd] = Wn[o + nd] + __ldg(&p.ia[nd]) * __ldg(&p.fr_w[e]) * fy;
  }
}

/// A contact or friction group at least this large (times the batch) is spread over the whole GPU.
#define WIDE 32

extern "C" __global__ void __launch_bounds__(TPB) vk_step(P p) {
  __shared__ double sh[5][TPB];
  __shared__ double sv[8];
  __shared__ double sbody[3 * 64];  // per body: cg, cu, savf (contact block)
  __shared__ Pred pr;
  const int tid = threadIdx.x;
  const int lane = tid & 31, warp = tid >> 5;
  const u32 nbk = gridDim.x;
  const ll B = p.B, NT = p.NT, NI = p.NI, NCI = p.NCI, NC = p.NC;
  const double k = p.k;
  const ll gwarp = ((ll)blockIdx.x * TPB + tid) >> 5;
  const ll nwarps = ((ll)nbk * TPB) >> 5;
  // phase timing on block 0 (cycles from one barrier exit to the next, summed in registers, written at exit)
  long long tprev = clock64(), pt[7] = {0, 0, 0, 0, 0, 0, 0};
  long long active_friction_steps = 0;
#define MARK(ph) if (blockIdx.x == 0 && tid == 0) { long long tn = clock64(); pt[ph] += tn - tprev; tprev = tn; }

  for (ll it = 0; it < p.steps; it++) {
    const ll n = p.n0 + it;
    const ll gstep = p.gstep0 + it;
    const ll par = n & 1;  // port state: read [par], write [1 - par]
    const double t = (double)gstep / p.rate;
    double* Wp = wbuf(p, n % 3);
    double* W = wbuf(p, (n + 1) % 3);
    double* Wn = wbuf(p, (n + 2) % 3);
    double* LW = (n % 2 == 0) ? p.LW0 : p.LW1;
    double* LWP = (n % 2 == 0) ? p.LW1 : p.LW0;

    // ------------------------------------------------------------------ P1
    // The P1 grid barrier publishes this reset before any contact is solved.
    if (blockIdx.x == 0 && tid == 0) *(volatile int*)p.any_friction = 0;
    for (ll wi = blockIdx.x; wi < B * NI; wi += nbk) {
      const int s = (int)(wi / NI), item = (int)(wi % NI);
      const int b = __ldg(&p.item_body[item]);
      const int i = __ldg(&p.item_start[item]) + tid;
      double v[4] = {0.0, 0.0, 0.0, 0.0};
      if (tid < __ldg(&p.item_len[item])) {
        const ll o = (ll)s * NT;
        double acc = 0.0;
#pragma unroll 5
        for (int q = 0; q < p.KL; q++) acc = acc + __ldg(&p.lapE_val[(ll)q * NT + i]) * W[o + __ldg(&p.lapE_col[(ll)q * NT + i])];
        const double lw = acc;
        LW[o + i] = lw;
        double f = __ldg(&p.bias[i]) - __ldg(&p.bten[b]) * lw;
        if (__ldg(&p.bstiff[b])) {
          double a2 = 0.0;
#pragma unroll 8
          for (int q = 0; q < p.KS; q++) a2 = a2 + __ldg(&p.stfE_val[(ll)q * NT + i]) * W[o + __ldg(&p.stfE_col[(ll)q * NT + i])];
          f = f + -1.0 * a2;
        }
        p.F[o + i] = f;
        if (__ldg(&p.bs1[b]) > 0.0) {
          const double correction = mass_flux_correction(p, i, o, W, Wp);
          p.TMP[o + i] = -((lw - LWP[o + i]) + correction / __ldg(&p.mass[i])) * __ldg(&p.bih2[b]);
        }
        const double w = W[o + i], wp = Wp[o + i];
        v[0] = w * lw; v[1] = w; v[2] = w - wp; v[3] = 2.0 * w - wp;
      }
      block_sum<4>(v, sh);
      if (tid == 0) for (int kk = 0; kk < 4; kk++) p.p1[((ll)s * NI + item) * 4 + kk] = v[kk];
    }
    for (ll wi = blockIdx.x; wi < B * NCI; wi += nbk) {
      const int s = (int)(wi / NCI), item = (int)(wi % NCI);
      const int c = __ldg(&p.citem_start[item]) + tid;
      if (tid < __ldg(&p.citem_len[item])) {
        const int b = __ldg(&p.citem_body[item]);
        const ll o = (ll)s * NT;
        const double ih2 = __ldg(&p.bih2[b]);
        if (__ldg(&p.bweak[b])) {
          double sum = 0.0;
          for (int j = __ldg(&p.q_ptr[c]); j < __ldg(&p.q_ptr[c+1]); j++) {
            const int a = __ldg(&p.q_a[j]), d = __ldg(&p.q_b[j]);
            const double wa = W[o+a], wd = W[o+d];
            sum += __ldg(&p.q_weight[j]) * (wa*wd + __ldg(&p.z0[a])*wd + wa*__ldg(&p.z0[d]));
          }
          const double qc = __ldg(&p.beh[b]) * sum;
          const ll oc = (ll)s * NC + c;
          p.Q[oc] = qc; p.QF[oc] = (float)qc;
          p.CXX[oc] = p.CYY[oc] = p.CXY[oc] = 0.0;
        } else {
        int nb[9]; double u[9];
#pragma unroll
        for (int j = 0; j < 9; j++) nb[j] = __ldg(&p.cnb[9ll * c + j]);
#pragma unroll
        for (int j = 0; j < 9; j++) u[j] = nb[j] < 0 ? 0.0 : W[o + nb[j]];
        const double xx = (u[1] - 2.0 * u[0] + u[2]) * ih2;
        const double yy = (u[3] - 2.0 * u[0] + u[4]) * ih2;
        const double xy = (u[5] - u[6] - u[7] + u[8]) * (0.25 * ih2);
        const double zx = __ldg(&p.crest[3ll * c]), zy = __ldg(&p.crest[3ll * c + 1]);
        const double zxy = __ldg(&p.crest[3ll * c + 2]);
        // Same expanded compatibility and rest tangent as the native solver.
        const double qc = -__ldg(&p.beh[b]) * (xx * yy - xy * xy + xx * zy + yy * zx - 2.0 * xy * zxy);
        const ll oc = (ll)s * NC + c;
        p.Q[oc] = qc; p.QF[oc] = (float)qc;
        p.CXX[oc] = xx + zx; p.CYY[oc] = yy + zy; p.CXY[oc] = xy + zxy;
        }
      }
    }
    // link forces, one thread per (simulation, link); the nodes add them in link order in P3
    for (ll e = (ll)tid * nbk + blockIdx.x; e < B * p.NLINK; e += (ll)nbk * TPB) {
      const ll s = e / p.NLINK, o = s * NT;
      const int l = (int)(e % p.NLINK);
      const int a = __ldg(&p.lk_a[l]), bb = __ldg(&p.lk_b[l]);
      const double wa = W[o + a], va = (W[o + a] - Wp[o + a]) / k;
      const double wb = bb >= 0 ? W[o + bb] : 0.0, vb = bb >= 0 ? (W[o + bb] - Wp[o + bb]) / k : 0.0;
      p.LF[e] = -__ldg(&p.lk_k[l]) * (wa - wb) - __ldg(&p.lk_c[l]) * (va - vb);
    }
    gsync(p.bar, nbk);
    MARK(0)

    // ------------------------------------------------------------------ P2
    if (p.NVK > 0) {
      if (B == 1) airy_rows<1>(p, gwarp, nwarps); else airy_rows<GB>(p, gwarp, nwarps);
      gsync(p.bar, nbk);
    }
    MARK(1)

    // ------------------------------------------------------------------ P3
    for (ll wi = blockIdx.x; wi < B * NI; wi += nbk) {
      const int s = (int)(wi / NI), item = (int)(wi % NI);
      const int b = __ldg(&p.item_body[item]);
      const int i = __ldg(&p.item_start[item]) + tid;
      if (warp == 0) {
        double sq[1] = {0.0};
        if (__ldg(&p.bber[b]) > 0.0) wsums<1>(sq, p.p1, 4, 0, __ldg(&p.bitem0[b]), __ldg(&p.bitem1[b]), s, NI);
        const int c = __ldg(&p.bcav[b]), r = __ldg(&p.brad[b]);
        const double cp = c >= 0 ? cavity_pressure(p, c, s, NI, par) : 0.0;
        const double rv = r >= 0 ? rad_value(p, r, s, NI) : 0.0;
        if (lane == 0) {
          sv[0] = __ldg(&p.bber[b]) * sq[0]; sv[1] = cp; sv[2] = rv;
          // the cavity's first head's first item advances its port (from the pressure at the old opening)
          if (c >= 0 && __ldg(&p.cav_port[c]) && b == __ldg(&p.cav_head[__ldg(&p.cav_hptr[c])]) && item == __ldg(&p.bitem0[b])) {
            const ll o0 = (par * B + s) * p.NCAV + c, o1 = ((1 - par) * B + s) * p.NCAV + c;
            const double du = k * (cp - __ldg(&p.port_res[c]) * p.port_u[o0]) * __ldg(&p.port_area[c]) / (p.gas_density * __ldg(&p.port_len[c]));
            const double u = p.port_u[o0] + du;
            p.port_u[o1] = u;
            p.port_du[o1] = du / k;
            p.port_x[o1] = p.port_x[o0] + k * u;
          }
        }
      }
      __syncthreads();
      double v[1] = {0.0};
      if (tid < __ldg(&p.item_len[item])) {
        const ll o = (ll)s * NT;
        const double h = __ldg(&p.bh[b]);
        double f = p.F[o + i];
#pragma unroll 4
        for (int q = __ldg(&p.lnk_ptr[i]); q < __ldg(&p.lnk_ptr[i + 1]); q++) {
          const int l = __ldg(&p.lnk_idx[q]);
          const double fl = p.LF[(ll)s * p.NLINK + l];
          f = (__ldg(&p.lk_a[l]) == i) ? f + fl : f - fl;
        }
        if (__ldg(&p.bcav[b]) >= 0) f = f + sv[1] * h * h;
        if (__ldg(&p.brad[b]) >= 0) f = f + sv[2] * h * h;
        for (int q = __ldg(&p.drv_ptr[i]); q < __ldg(&p.drv_ptr[i + 1]); q++) f = f + drive_at(p, __ldg(&p.drv_idx[q]), t);
        p.F[o + i] = f;
        if (__ldg(&p.bber[b]) > 0.0 || __ldg(&p.bvk[b])) {
          double g = 0.0;
          if (__ldg(&p.bber[b]) > 0.0) g = g + sv[0] * LW[o + i];
          if (__ldg(&p.bvk[b])) {
            if (__ldg(&p.bweak[b])) {
              for (int q = __ldg(&p.wg_ptr[i]); q < __ldg(&p.wg_ptr[i+1]); q++) {
                const int other = __ldg(&p.wg_other[q]);
                const ll oc = (ll)s * NC + __ldg(&p.wg_core[q]);
                g += (h*h*p.AF[oc]*__ldg(&p.wg_weight[q])) * (W[o+other]+__ldg(&p.z0[other]));
              }
            } else {
#pragma unroll 9
            for (int q = 0; q < p.KG; q++) {
              const int gc = __ldg(&p.gtE_core[(ll)q * NT + i]);
              if (gc >= 0) {
                const ll oc = (ll)s * NC + gc;
                g = g - p.AF[oc] * slot_val(__ldg(&p.gtE_slot[(ll)q * NT + i]), p.CXX[oc], p.CYY[oc], p.CXY[oc]);
              }
            }
            }
          }
          p.G[o + i] = g;
          v[0] = g * Wp[o + i];
        }
      }
      block_sum<1>(v, sh);
      if (tid == 0) p.p3[(ll)s * NI + item] = v[0];
    }
    for (ll wi = blockIdx.x; wi < B * NCI; wi += nbk) {
      const int s = (int)(wi / NCI), item = (int)(wi % NCI);
      const int c = __ldg(&p.citem_start[item]) + tid;
      double v[1] = {0.0};
      if (tid < __ldg(&p.citem_len[item])) { const ll oc = (ll)s * NC + c; v[0] = p.AF[oc] * p.Q[oc]; }
      block_sum<1>(v, sh);
      if (tid == 0) p.pc3[(ll)s * NCI + item] = v[0];
    }
    gsync(p.bar, nbk);
    MARK(2)

    // ------------------------------------------------------------------ P4
    for (ll wi = blockIdx.x; wi < B * NI; wi += nbk) {
      const int s = (int)(wi / NI), item = (int)(wi % NI);
      const int b = __ldg(&p.item_body[item]);
      const int i = __ldg(&p.item_start[item]) + tid;
      if (warp == 0) {
        Pred r = predict_scalars(p, s, b);
        if (lane == 0) {
          pr = r;
          if (item == __ldg(&p.bitem0[b])) {
            const ll ob = (ll)s * p.NB + b;
            p.psi1[ob] = r.psi1; p.rprev1[ob] = r.rprev1; p.savf[ob] = r.scale ? 1.0 : 0.0;
          }
        }
      }
      __syncthreads();
      double v[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      if (tid < __ldg(&p.item_len[item])) {
        const ll o = (ll)s * NT + i;
        const bool nl = __ldg(&p.bber[b]) > 0.0 || __ldg(&p.bvk[b]);
        double g = 0.0;
        if (nl) { g = pr.scale ? p.G[o] * pr.ir : 0.0; p.G[o] = g; }
        double force = p.F[o];
        if (__ldg(&p.bviscous[b])) {
          double cv = 0.0;
          const ll base = (ll)s * NT;
          for (int q = __ldg(&p.visc_ptr[i]); q < __ldg(&p.visc_ptr[i + 1]); q++) {
            const int j = __ldg(&p.visc_col[q]);
            cv = cv + __ldg(&p.visc_val[q]) * (W[base + j] - Wp[base + j]);
          }
          force = force + (-1.0 / k) * cv;
          p.F[o] = force;
        }
        double rhs = __ldg(&p.cw[i]) * W[o] - __ldg(&p.cwp[i]) * Wp[o] + force + pr.extra;
        if (__ldg(&p.bs1[b]) > 0.0) rhs = rhs + __ldg(&p.cs1[i]) * p.TMP[o];
        if (pr.scale) rhs = rhs - g * pr.sc;
        const double ia = __ldg(&p.ia[i]);
        const double x = rhs * ia;
        Wn[o] = x;
        if (pr.scale) {
          const double z = g * ia;
          p.TMP[o] = z;
          v[0] = g * x; v[1] = g * z; v[2] = pr.h2 * z;
        }
        if (pr.mu > 0.0) { v[3] = pr.h2 * x; v[4] = pr.h2 * pr.h2 * ia; }
      }
      block_sum<5>(v, sh);
      if (tid == 0) for (int kk = 0; kk < 5; kk++) p.p4[((ll)s * NI + item) * 5 + kk] = v[kk];
    }
    gsync(p.bar, nbk);
    MARK(3)

    // ------------------------------------------------------------------ P5
    // (a) the Woodbury correction on nodes no contact or friction touches, by the nodes' blocks
    for (ll wi = blockIdx.x; wi < B * NI; wi += nbk) {
      const int s = (int)(wi / NI), item = (int)(wi % NI);
      const int b = __ldg(&p.item_body[item]);
      const int i = __ldg(&p.item_start[item]) + tid;
      if (warp == 0) {
        double cg, cu;
        woodbury(p, s, b, cg, cu);
        if (lane == 0) { sv[3] = cg; sv[4] = cu; sv[5] = __ldg(&p.bh[b]) * __ldg(&p.bh[b]) * cu; sv[6] = p.savf[(ll)s * p.NB + b]; }
      }
      __syncthreads();
      if (tid < __ldg(&p.item_len[item]) && !__ldg(&p.cnode[i])) {
        const double cg = sv[3], cu = sv[4], hc = sv[5];
        if (cg != 0.0 || cu != 0.0) {
          const ll o = (ll)s * NT + i;
          const double zg = sv[6] > 0.0 ? p.TMP[o] : 0.0;
          Wn[o] = Wn[o] - (zg * cg + hc * __ldg(&p.ia[i]));
        }
      }
      __syncthreads();
    }
    // (b) one block per simulation: the same correction on contact nodes, strikes, contacts, frictions
    for (ll s = nbk - 1 - blockIdx.x; s < B; s += nbk) {
      const ll o = s * NT;
      for (int b = warp; b < p.NB; b += TPB / 32) {
        double cg, cu;
        woodbury(p, (int)s, b, cg, cu);
        if (lane == 0) { sbody[3 * b] = cg; sbody[3 * b + 1] = cu; sbody[3 * b + 2] = p.savf[s * p.NB + b]; }
      }
      __syncthreads();
      for (int q = tid; q < p.NCN; q += TPB) {
        const int i = __ldg(&p.cn_node[q]), b = __ldg(&p.cn_body[q]);
        const double cg = sbody[3 * b], cu = sbody[3 * b + 1];
        if (cg != 0.0 || cu != 0.0) {
          const double hc = __ldg(&p.bh[b]) * __ldg(&p.bh[b]) * cu;
          const double zg = sbody[3 * b + 2] > 0.0 ? p.TMP[o + i] : 0.0;
          Wn[o + i] = Wn[o + i] - (zg * cg + hc * __ldg(&p.ia[i]));
        }
      }
      __syncthreads();
      if (tid == 0) {
        int cur = p.st_cursor[s];
        while (cur < __ldg(&p.st_ptr[s + 1]) && !(t < __ldg(&p.st_time[cur]))) {
          const int sb = __ldg(&p.st_body[cur]), j = __ldg(&p.st_contact[cur]);
          for (int q = __ldg(&p.bc_ptr[sb]); q < __ldg(&p.bc_ptr[sb + 1]); q++) {
            const int c = __ldg(&p.bc_idx[q]);
            if (c != j) p.ct_active[s * p.NCON + c] = 0;
          }
          double surf = 0.0;
          for (int e = __ldg(&p.ct_a0[j]); e < __ldg(&p.ct_b0[j]); e++) surf = surf + __ldg(&p.ct_w[e]) * W[o + __ldg(&p.ct_node[e])];
          const double z0 = surf - __ldg(&p.ct_gap[j]) + __ldg(&p.st_clear[cur]);
          const int gs = __ldg(&p.item_start[__ldg(&p.bitem0[sb])]);
          W[o + gs] = z0;
          Wp[o + gs] = z0 + __ldg(&p.st_speed[cur]) * k;
          Wn[o + gs] = z0 - __ldg(&p.st_speed[cur]) * k;
          p.ct_active[s * p.NCON + j] = 1;
          p.ct_dnow[s * p.NCON + j] = surf - z0 - __ldg(&p.ct_gap[j]);
          cur++;
        }
        p.st_cursor[s] = cur;
      }
      __syncthreads();
    }
    // contacts, then frictions, group by group. A group large enough to fill the GPU is spread over every
    // block after a grid barrier; a small one runs on its simulation's block (a block barrier between groups).
    {
      bool wide_before = false;
      for (int bt = 0; bt < p.NCB + p.NFB; bt++) {
        if (bt == p.NCB) {
          // Include contacts solved on other blocks before reading the witness.
          // If every friction call would return immediately, its hundreds of
          // ordered batch barriers have no work to protect. Otherwise retain
          // the original batches and every floating-point operation unchanged.
          gsync(p.bar, nbk);
          // Volatile bypasses per-SM L1 caching of this cross-block witness.
          if (!*(volatile int*)p.any_friction) break;
          if (blockIdx.x == 0 && tid == 0) active_friction_steps++;
        }
        const bool fr = bt >= p.NCB;
        const int g = fr ? bt - (int)p.NCB : bt;
        const int q0 = fr ? __ldg(&p.fb_ptr[g]) : __ldg(&p.cb_ptr[g]);
        const int q1 = fr ? __ldg(&p.fb_ptr[g + 1]) : __ldg(&p.cb_ptr[g + 1]);
        const ll size = q1 - q0;
        const bool wide = size * B >= WIDE;
        if (wide || wide_before) gsync(p.bar, nbk);
        if (wide) {
          // block-major: item e on block e % nbk, so a group spreads over every multiprocessor's FP64 units
          for (ll e = (ll)tid * nbk + blockIdx.x; e < B * size; e += (ll)nbk * TPB) {
            const ll s = e / size;
            const int q = q0 + (int)(e % size);
            if (fr) solve_friction(p, s, __ldg(&p.fb_idx[q]), W, Wn); else solve_contact(p, s, __ldg(&p.cb_idx[q]), Wp, Wn, k);
          }
        } else {
          for (ll s = nbk - 1 - blockIdx.x; s < B; s += nbk)
            for (int q = q0 + tid; q < q1; q += TPB) {
              if (fr) solve_friction(p, s, __ldg(&p.fb_idx[q]), W, Wn); else solve_contact(p, s, __ldg(&p.cb_idx[q]), Wp, Wn, k);
            }
          __syncthreads();
        }
        wide_before = wide;
      }
    }
    gsync(p.bar, nbk);
    MARK(4)

    // ------------------------------------------------------------------ P6
    for (ll wi = blockIdx.x; wi < B * NI; wi += nbk) {
      const int s = (int)(wi / NI), item = (int)(wi % NI);
      const int b = __ldg(&p.item_body[item]);
      const int i = __ldg(&p.item_start[item]) + tid;
      double v[4] = {0.0, 0.0, 0.0, 0.0};
      const bool sav = p.savf[(ll)s * p.NB + b] > 0.0;
      const bool vk = __ldg(&p.bvk[b]) != 0;
      const bool viscous = __ldg(&p.bviscous[b]) != 0;
      if ((sav || vk || viscous) && tid < __ldg(&p.item_len[item])) {
        const ll o = (ll)s * NT + i;
        const double dv = Wn[o] - Wp[o];
        if (vk) {
          v[0] = __ldg(&p.ms0[i]) * dv * dv;
          const double correction = mass_flux_correction(p, i, (ll)s * NT, W, Wp);
          v[1] = (__ldg(&p.mass[i]) * (LW[o] - LWP[o]) + correction) * dv;
        }
        if (sav) v[2] = p.G[o] * dv;
        if (viscous) {
          const ll base = (ll)s * NT;
          // Same quarter-cell tensor-root rows as the native solver: squares of
          // G(wn-wp), not a clipped signed matrix quadratic or preset T60.
          for (int r = BENDING_ROOT_ROWS_PER_NODE * i; r < BENDING_ROOT_ROWS_PER_NODE * (i + 1); r++) {
            double strain = 0.0;
            for (int q = __ldg(&p.vroot_ptr[r]); q < __ldg(&p.vroot_ptr[r + 1]); q++) {
              const int j = __ldg(&p.vroot_col[q]);
              strain = strain + __ldg(&p.vroot_val[q]) * (Wn[base + j] - Wp[base + j]);
            }
            v[3] = v[3] + strain * strain;
          }
        }
      }
      if (sav || vk || viscous) {
        block_sum<4>(v, sh);
        if (tid == 0) for (int kk = 0; kk < 4; kk++) p.p7[((ll)s * NI + item) * 4 + kk] = v[kk];
      }
    }
    // microphone taps: one warp per (simulation, segment of at most 64 taps of one delay)
    for (ll e = gwarp; e < B * p.NSEG; e += nwarps) {
      const int s = (int)(e / p.NSEG), sg = (int)(e % p.NSEG);
      const ll o = (ll)s * NT;
      double acc = 0.0;
      for (int q = __ldg(&p.seg_q0[sg]) + lane; q < __ldg(&p.seg_q0[sg + 1]); q += 32) {
        const int nd = __ldg(&p.md_node[q]);
        acc = acc + __ldg(&p.md_w[q]) * (Wn[o + nd] - 2.0 * W[o + nd] + Wp[o + nd]);
      }
      acc = bfly(acc);
      if (lane == 0) p.seg_part[(ll)s * p.NSEG + sg] = acc;
    }
    gsync(p.bar, nbk);
    MARK(5)

    // ------------------------------------------------------------------ P7
    if (blockIdx.x == 0) {
      for (int sb = warp; sb < B * p.NB; sb += TPB / 32) {
        const int s = sb / (int)p.NB, b = sb % (int)p.NB;
        const ll ob = (ll)s * p.NB + b;
        const bool sav = p.savf[ob] > 0.0;
        const bool vk = __ldg(&p.bvk[b]) != 0;
        const bool viscous = __ldg(&p.bviscous[b]) != 0;
        double d[4] = {0.0, 0.0, 0.0, 0.0};
        if (sav || vk || viscous) wsums<4>(d, p.p7, 4, 0, __ldg(&p.bitem0[b]), __ldg(&p.bitem1[b]), s, NI);
        if (lane == 0) {
          const double bending_loss = d[3] / (4.0 * k);
          if (viscous) p.viscous_dissipated[ob] = p.viscous_dissipated[ob] + bending_loss;
          double psi = p.psi1[ob];
          // Refresh shell credit even when the nonlinear gradient vanishes.
          if (vk) {
            const double ih2 = 1.0 / (__ldg(&p.bh[b]) * __ldg(&p.bh[b]));
            const double loss = d[0] / (2.0 * k) + __ldg(&p.bs1[b]) * ih2 / k * d[1] + bending_loss;
            p.budgetA[ob] = loss < 0.0 ? 0.0 : loss; // signed net loss; preserve NaN
          }
          if (sav) psi = psi + 0.5 * d[2];
          p.psiA[ob] = psi;
          p.rprevA[ob] = p.rprev1[ob];
        }
      }
      // each delay's segments, in order, into its ring slot
      for (ll e = tid; e < B * p.NMD; e += TPB) {
        const int s = (int)(e / p.NMD), md = (int)(e % p.NMD);
        double acc = 0.0;
        for (int sg = __ldg(&p.md_seg[md]); sg < __ldg(&p.md_seg[md + 1]); sg++) acc = acc + p.seg_part[(ll)s * p.NSEG + sg];
        const int m = __ldg(&p.md_mic[md]);
        double* ring = p.ring + (ll)s * p.RINGTOT + __ldg(&p.mic_roff[m]);
        const ll slot = (gstep + __ldg(&p.md_delay[md])) & __ldg(&p.mic_rmask[m]);
        ring[slot] = ring[slot] + acc;
      }
      __syncthreads();
      for (int sm = tid; sm < B * p.NM; sm += TPB) {
        const int s = sm / (int)p.NM, m = sm % (int)p.NM;
        double* ring = p.ring + (ll)s * p.RINGTOT + __ldg(&p.mic_roff[m]);
        const ll mask = __ldg(&p.mic_rmask[m]);
        for (int q = __ldg(&p.pt_ptr[m]); q < __ldg(&p.pt_ptr[m + 1]); q++) {
          const double a = p.port_du[((1 - par) * B + s) * p.NCAV + __ldg(&p.pt_cav[q])];
          ring[(gstep + __ldg(&p.pt_delay[q])) & mask] = ring[(gstep + __ldg(&p.pt_delay[q])) & mask] + __ldg(&p.pt_lo[q]) * a;
          ring[(gstep + __ldg(&p.pt_delay[q]) + 1) & mask] = ring[(gstep + __ldg(&p.pt_delay[q]) + 1) & mask] + __ldg(&p.pt_hi[q]) * a;
        }
        const ll slot = gstep & mask;
        double y = ring[slot];
        if (__ldg(&p.mic_tau[m]) >= 0.0) {
          const double g = 2.0 * __ldg(&p.mic_tau[m]) * p.rate;
          double* hp = p.mic_hp + ((ll)s * p.NM + m) * 2;
          const double yn = (g * (y - hp[0]) - (1.0 - g) * hp[1]) / (1.0 + g);
          hp[0] = y; hp[1] = yn;
          y = yn;
        }
        p.out[(it * B + s) * p.NM + m] = y;
        ring[slot] = 0.0;
      }
    }
    MARK(6)
  }
  if (blockIdx.x == 0 && tid == 0) {
    for (int ph = 0; ph < 7; ph++) p.prof[ph] += pt[ph];
    p.prof[7] += active_friction_steps;
  }
}
