//! Acquisition functions: closed-form Expected Improvement and
//! batched q-EI through the Cholesky reparameterization f = μ + L·z
//! with FIXED scrambled-Sobol normal samples — the acquisition
//! surface is DETERMINISTIC, so its optimization replays bitwise.
//! Normal CDF/quantile via deterministic polynomial approximations
//! (Abramowitz–Stegun 7.1.26 for Φ, ~1.5e−7 absolute accuracy;
//! Acklam for Φ⁻¹, ~1e−9 relative) — documented accuracy, no
//! platform libm anywhere in the pipeline.

use crate::gp::Gp;

/// Standard normal PDF.
#[must_use]
pub fn phi_pdf(z: f64) -> f64 {
    fs_math::det::exp(-0.5 * z * z) / fs_math::det::sqrt(2.0 * core::f64::consts::PI)
}

/// Standard normal CDF via the Abramowitz–Stegun 7.1.26 polynomial
/// (max absolute error ≈ 1.5e−7 — plenty for acquisition ranking;
/// DOCUMENTED, deterministic).
#[must_use]
pub fn phi_cdf(z: f64) -> f64 {
    let x = z / fs_math::det::sqrt(2.0);
    let (sign, ax): (f64, f64) = if x < 0.0 { (-1.0, -x) } else { (1.0, x) };
    let t = 1.0 / 0.3275911f64.mul_add(ax, 1.0);
    // Horner: ((((a5·t + a4)·t + a3)·t + a2)·t + a1)·t.
    let mut poly = 1.061405429f64;
    poly = poly.mul_add(t, -1.453152027);
    poly = poly.mul_add(t, 1.421413741);
    poly = poly.mul_add(t, -0.284496736);
    poly = poly.mul_add(t, 0.254829592);
    poly *= t;
    let erf = 1.0 - poly * fs_math::det::exp(-ax * ax);
    0.5 * sign.mul_add(erf, 1.0)
}

/// Standard normal quantile Φ⁻¹ (Acklam's algorithm, relative error
/// ≈ 1e−9): turns Sobol uniforms into deterministic normal samples.
/// Constants are Acklam's verbatim (the lint about excess digits is
/// suppressed: truncating published constants invites transcription
/// bugs; the compiler rounds to the nearest f64 deterministically).
#[must_use]
#[allow(clippy::excessive_precision)]
pub fn phi_inv(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_690e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    assert!(p > 0.0 && p < 1.0, "quantile needs p in (0,1), got {p}");
    let p_low = 0.02425;
    if p < p_low {
        let q = fs_math::det::sqrt(-2.0 * fs_math::det::ln(p));
        return (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
    }
    if p > 1.0 - p_low {
        let q = fs_math::det::sqrt(-2.0 * fs_math::det::ln(1.0 - p));
        return -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
    }
    let q = p - 0.5;
    let r = q * q;
    (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
        / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
}

/// Closed-form Expected Improvement for MINIMIZATION at exploration
/// margin ξ: EI(x) = (f* − μ − ξ)Φ(z) + σφ(z), z = (f* − μ − ξ)/σ.
#[must_use]
pub fn expected_improvement(gp: &Gp, x: &[f64], f_best: f64, xi: f64) -> f64 {
    let (mu, var) = gp.predict(x);
    let sigma = fs_math::det::sqrt(var.max(1e-18));
    let delta = f_best - mu - xi;
    let z = delta / sigma;
    delta.mul_add(phi_cdf(z), sigma * phi_pdf(z)).max(0.0)
}

/// Deterministic q-point normal sample bank: `samples × q` z-values
/// from scrambled Sobol through Φ⁻¹ (fixed per seed — the
/// reparameterization's common random numbers).
#[must_use]
pub fn normal_bank(samples: usize, q: usize, seed: u64) -> Vec<f64> {
    let sobol = fs_rand::qmc::Sobol::scrambled(q, seed);
    let mut bank = vec![0.0f64; samples * q];
    let mut pt = vec![0.0f64; q];
    for s in 0..samples {
        // Skip the first point (all-zeros before scrambling can sit
        // at the box corner); 1-based indexing keeps u in (0,1).
        sobol.point(u32::try_from(s + 1).expect("bank fits u32"), &mut pt);
        for (j, &u) in pt.iter().enumerate() {
            bank[s * q + j] = phi_inv(u.clamp(1e-12, 1.0 - 1e-12));
        }
    }
    bank
}

/// Batched q-EI for MINIMIZATION via the reparameterization
/// f = μ + L·z over the fixed bank: mean of
/// max(0, f* − minᵢ fᵢ(z)). Deterministic given the bank.
#[must_use]
pub fn q_expected_improvement(gp: &Gp, xs: &[Vec<f64>], f_best: f64, bank: &[f64]) -> f64 {
    let q = xs.len();
    let (mu, lflat) = gp.predict_joint(xs);
    let samples = bank.len() / q;
    let mut acc = 0.0f64;
    for s in 0..samples {
        let z = &bank[s * q..(s + 1) * q];
        let mut best = f64::INFINITY;
        for i in 0..q {
            let mut f = mu[i];
            for j in 0..=i {
                f = lflat[i * q + j].mul_add(z[j], f);
            }
            best = best.min(f);
        }
        acc += (f_best - best).max(0.0);
    }
    acc / samples as f64
}
