//! fs-render — unbiased spectral path-tracing core. Layer: L5.
//!
//! Rendering is just another budgeted, cancellable, REPLAYABLE study. This v0
//! is the verifiable Monte-Carlo foundation the whole tracer rests on — the
//! parts with exact or analytic checks, so correctness is a certificate rather
//! than a look:
//!
//! - deterministic LOW-DISCREPANCY sampling ([`radical_inverse`] / [`halton`]) —
//!   an image is as replayable as a solve;
//! - cosine-weighted hemisphere sampling + the LAMBERTIAN FURNACE test
//!   ([`Lambertian::furnace_radiance`]) — a closed uniform-emission scene returns
//!   EXACTLY `albedo · radiance` (energy conservation, zero variance);
//! - MIS balance/power heuristics ([`balance_heuristic`], [`power_heuristic`])
//!   with a WEIGHT-SUM audit ([`mis_weight_sum`]) — no energy lost or gained at
//!   strategy boundaries;
//! - HERO-WAVELENGTH spectral integration ([`spectral_integral`]) — an unbiased
//!   estimate of a spectral integral.
//!
//! Deterministic; no dependencies.

use core::f64::consts::PI;

#[cfg(feature = "chart-backends")]
pub mod charts;

#[cfg(feature = "volumes")]
pub mod volumes;

/// The van der Corput / Halton radical inverse of `i` in the given `base`
/// (a deterministic low-discrepancy coordinate in `[0, 1)`).
#[must_use]
pub fn radical_inverse(base: u32, mut i: u64) -> f64 {
    debug_assert!(base >= 2);
    let b = u64::from(base);
    let bf = f64::from(base);
    let mut inv = 1.0 / bf;
    let mut result = 0.0;
    while i > 0 {
        let digit = (i % b) as f64;
        result += digit * inv;
        inv /= bf;
        i /= b;
    }
    result
}

const PRIMES: [u32; 8] = [2, 3, 5, 7, 11, 13, 17, 19];

/// The `dim`-th Halton coordinate of sample `i` (using the first primes).
///
/// # Panics
/// If `dim >= 8`.
#[must_use]
pub fn halton(dim: usize, i: u64) -> f64 {
    radical_inverse(PRIMES[dim], i)
}

/// Cosine-weighted hemisphere sample from `(u1, u2)` in the local frame
/// (`z` up): returns `(direction, pdf)` with `pdf = cosθ/π`.
#[must_use]
pub fn cosine_sample_hemisphere(u1: f64, u2: f64) -> ([f64; 3], f64) {
    let r = u1.sqrt();
    let phi = 2.0 * PI * u2;
    let x = r * phi.cos();
    let y = r * phi.sin();
    let z = (1.0 - u1).max(0.0).sqrt(); // cosθ
    ([x, y, z], z / PI)
}

/// A grayscale Lambertian BRDF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lambertian {
    /// The albedo (reflectance in `[0, 1]`).
    pub albedo: f64,
}

impl Lambertian {
    /// The BRDF value `ρ/π` (constant for a Lambertian).
    #[must_use]
    pub fn brdf(&self) -> f64 {
        self.albedo / PI
    }

    /// The FURNACE test: Monte-Carlo estimate of the reflected radiance under a
    /// uniform incident radiance. Energy conservation demands exactly
    /// `albedo · incident`, and cosine-weighted importance sampling delivers it
    /// with ZERO variance (every sample equals `albedo · incident`).
    #[must_use]
    pub fn furnace_radiance(&self, incident: f64, samples: usize) -> f64 {
        let mut acc = 0.0;
        for i in 1..=samples as u64 {
            let (dir, pdf) = cosine_sample_hemisphere(radical_inverse(2, i), radical_inverse(3, i));
            let cos_theta = dir[2];
            acc += self.brdf() * incident * cos_theta / pdf;
        }
        acc / samples as f64
    }
}

/// The MIS BALANCE-heuristic weight for strategy `f` given sample counts and the
/// two strategies' pdfs at the sample.
#[must_use]
pub fn balance_heuristic(nf: u32, pf: f64, ng: u32, pg: f64) -> f64 {
    let (a, b) = (f64::from(nf) * pf, f64::from(ng) * pg);
    let denom = a + b;
    if denom <= 0.0 { 0.0 } else { a / denom }
}

/// The MIS POWER-heuristic weight (β = 2) for strategy `f`.
#[must_use]
pub fn power_heuristic(nf: u32, pf: f64, ng: u32, pg: f64) -> f64 {
    let (a, b) = (f64::from(nf) * pf, f64::from(ng) * pg);
    let (a2, b2) = (a * a, b * b);
    let denom = a2 + b2;
    if denom <= 0.0 { 0.0 } else { a2 / denom }
}

/// The MIS WEIGHT-SUM audit: the two balance-heuristic weights at a sample must
/// sum to `1` (no energy lost or gained at the strategy boundary). Returns the
/// sum (nominally `1.0`).
#[must_use]
pub fn mis_weight_sum(pf: f64, pg: f64) -> f64 {
    balance_heuristic(1, pf, 1, pg) + balance_heuristic(1, pg, 1, pf)
}

/// A set of `count` hero wavelengths: the `hero` plus evenly rotated companions
/// wrapped into `[min_nm, max_nm)` (the stratified hero-wavelength set).
#[must_use]
pub fn hero_wavelengths(hero_nm: f64, count: usize, min_nm: f64, max_nm: f64) -> Vec<f64> {
    let range = max_nm - min_nm;
    (0..count)
        .map(|k| {
            let off = range * k as f64 / count as f64;
            min_nm + (hero_nm - min_nm + off).rem_euclid(range)
        })
        .collect()
}

/// An unbiased hero-wavelength estimate of `∫ spectrum(λ) dλ` over
/// `[min_nm, max_nm]`, using `samples` stratified hero draws of 4 wavelengths.
#[must_use]
pub fn spectral_integral(
    spectrum: impl Fn(f64) -> f64,
    min_nm: f64,
    max_nm: f64,
    samples: usize,
) -> f64 {
    let range = max_nm - min_nm;
    let count = 4;
    let mut acc = 0.0;
    for i in 1..=samples as u64 {
        let hero = min_nm + radical_inverse(2, i) * range;
        let set = hero_wavelengths(hero, count, min_nm, max_nm);
        let avg: f64 = set.iter().map(|&l| spectrum(l)).sum::<f64>() / count as f64;
        acc += avg * range;
    }
    acc / samples as f64
}

/// A tiny deterministic MIS demonstrator: estimate `∫₀¹ f(x) dx` by combining
/// UNIFORM sampling (`pdf = 1`) with LINEAR importance sampling (`pdf(x) = 2x`,
/// drawn by inverse CDF `x = √u`), weighted by the balance heuristic — unbiased.
#[must_use]
pub fn mis_integrate_unit(f: impl Fn(f64) -> f64, n: usize) -> f64 {
    let nu = n as u32;
    let mut est = 0.0;
    for i in 1..=n as u64 {
        // uniform strategy.
        let xu = radical_inverse(2, i);
        let (pf_u, pg_u) = (1.0, 2.0 * xu);
        est += balance_heuristic(nu, pf_u, nu, pg_u) * f(xu) / pf_u;
        // linear-importance strategy.
        let xl = radical_inverse(3, i).sqrt();
        let (pf_l, pg_l) = (1.0, 2.0 * xl);
        est += balance_heuristic(nu, pg_l, nu, pf_l) * f(xl) / pg_l;
    }
    est / n as f64
}
