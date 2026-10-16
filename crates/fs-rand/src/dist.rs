//! Extended distributions (plan §6.7, bead 6ys.19): gamma, beta,
//! Dirichlet, categorical (alias method), von Mises–Fisher, truncated
//! variants — all on strict-mode arithmetic with DETERMINISTIC
//! CONSUMPTION CONTRACTS: every rejection advances the stream index like
//! any other draw, so the consumed count is a pure function of stream
//! content (replay-safe, the `next_below` doctrine extended to
//! continuous distributions). Fixed-draw samplers are documented as such.
//!
//! The ziggurat normal and batched bulk Philox fills live in their own
//! explicit fast paths; these extended strict samplers continue to use
//! Box–Muller and scalar draws where that preserves their documented
//! consumption contracts.

use crate::Stream;
use fs_math::det;

impl Stream {
    /// Gamma(α, 1) via Marsaglia–Tsang squeeze (α ≥ 1) with the
    /// deterministic rejection contract; α < 1 uses the Ahrens–Dieter
    /// boost Γ(α) = Γ(α+1)·U^(1/α). Consumption: variable but
    /// content-determined (replay-tested). Panics (structured) for
    /// α ≤ 0.
    #[must_use]
    pub fn next_gamma(&mut self, alpha: f64) -> f64 {
        assert!(
            alpha > 0.0 && alpha.is_finite(),
            "gamma shape must be positive: {alpha}"
        );
        if alpha < 1.0 {
            // Boost: one extra uniform, then the α+1 sampler.
            let u = self.next_f64();
            let g = self.next_gamma(alpha + 1.0);
            // U^(1/α) via strict pow; U ∈ [0,1) → guard the 0 corner.
            return g * det::pow(u.max(f64::MIN_POSITIVE), 1.0 / alpha);
        }
        let (d, v3) = self.next_gamma_base(alpha);
        d * v3
    }

    // Retain the positive α ≥ 1 factors separately for log-domain replay.
    fn next_gamma_base(&mut self, alpha: f64) -> (f64, f64) {
        let d = alpha - 1.0 / 3.0;
        let c = 1.0 / det::sqrt(9.0 * d);
        loop {
            let x = self.next_normal();
            let v = 1.0 + c * x;
            if v <= 0.0 {
                continue; // rejected: index already advanced (2 draws)
            }
            let v3 = v * v * v;
            let u = self.next_f64();
            let x2 = x * x;
            // Squeeze acceptance (cheap path).
            if u < 1.0 - 0.0331 * x2 * x2 {
                return (d, v3);
            }
            // Log acceptance (exact path).
            if det::ln(u) < 0.5f64.mul_add(x2, d * (1.0 - v3 + det::ln(v3))) {
                return (d, v3);
            }
        }
    }

    // `scale <= min(1, alpha)` keeps this log score finite.
    fn next_scaled_log_gamma(&mut self, alpha: f64, scale: f64) -> f64 {
        if alpha < 1.0 {
            let u = self.next_f64();
            let base = self.next_scaled_log_gamma(alpha + 1.0, scale);
            return det::ln(u.max(f64::MIN_POSITIVE)).mul_add(scale / alpha, base);
        }
        let (d, v3) = self.next_gamma_base(alpha);
        scale * (det::ln(d) + det::ln(v3))
    }

    /// Beta(α, β) as Gα/(Gα + Gβ) (inherits the rejection contract).
    #[must_use]
    pub fn next_beta(&mut self, alpha: f64, beta: f64) -> f64 {
        let checkpoint = *self;
        let ga = self.next_gamma(alpha);
        let gb = self.next_gamma(beta);
        let total = ga + gb;
        if total > 0.0 && total.is_finite() {
            ga / total
        } else {
            let mut normalized = [0.0; 2];
            normalize_gamma_replay(checkpoint, &[alpha, beta], &mut normalized);
            normalized[0]
        }
    }

    /// Dirichlet(α₁..α_k): k gammas normalized, written into `out`.
    pub fn next_dirichlet(&mut self, alphas: &[f64], out: &mut [f64]) {
        assert_eq!(alphas.len(), out.len(), "dirichlet output length mismatch");
        assert!(!alphas.is_empty(), "dirichlet needs at least one component");
        let checkpoint = *self;
        let mut total = 0.0f64;
        for (slot, &a) in out.iter_mut().zip(alphas) {
            let g = self.next_gamma(a);
            *slot = g;
            total += g;
        }
        if total > 0.0 && total.is_finite() {
            for slot in out.iter_mut() {
                *slot /= total;
            }
        } else {
            normalize_gamma_replay(checkpoint, alphas, out);
        }
    }

    /// Truncated exponential(1) on [0, cap] via pure inversion — exactly
    /// ONE draw consumed (fixed-consumption contract).
    #[must_use]
    pub fn next_truncated_exponential(&mut self, cap: f64) -> f64 {
        assert!(cap > 0.0, "cap must be positive");
        let u = self.next_f64();
        // Inverse CDF of the truncated exponential:
        // x = −ln(1 − u·(1 − e^(−cap))).
        -det::ln(1.0 - u * (1.0 - det::exp(-cap)))
    }

    /// Standard normal truncated to [lo, ∞) via Robert's exponential
    /// rejection for lo > 0 (deterministic rejection contract), plain
    /// resampling for lo ≤ 0. Consumption content-determined.
    #[must_use]
    pub fn next_truncated_normal(&mut self, lo: f64) -> f64 {
        assert!(lo.is_finite(), "truncation bound must be finite");
        if lo <= 0.0 {
            // Acceptance probability ≥ ½: simple resampling is fine.
            loop {
                let z = self.next_normal();
                if z >= lo {
                    return z;
                }
            }
        }
        // Robert (1995): proposal Exp(λ*) shifted to lo, λ* the optimal rate.
        let lam = f64::midpoint(lo, det::sqrt(lo.mul_add(lo, 4.0)));
        loop {
            let e = self.next_exponential();
            let z = lo + e / lam;
            let diff = z - lam;
            let rho = det::exp(-0.5 * diff * diff);
            let u = self.next_f64();
            if u <= rho {
                return z;
            }
        }
    }

    /// Unit vector on S² from the von Mises–Fisher distribution with mean
    /// direction `mu` (unit) and concentration κ > 0. FIXED consumption:
    /// exactly 2 draws (Ulrich/Wood inversion for the polar component —
    /// no rejection). κ = 0 callers should use a uniform-sphere sampler.
    #[must_use]
    pub fn next_vmf3(&mut self, mu: [f64; 3], kappa: f64) -> [f64; 3] {
        assert!(
            kappa > 0.0 && kappa.is_finite(),
            "vMF needs kappa > 0: {kappa}"
        );
        let u = self.next_f64();
        let phi_u = self.next_f64();
        // Polar component by inversion (Ulrich 1984, S² closed form):
        // w = 1 + ln(u + (1−u)·e^(−2κ))/κ.
        let w =
            1.0 + det::ln(u.max(f64::MIN_POSITIVE) + (1.0 - u) * det::exp(-2.0 * kappa)) / kappa;
        let w = w.clamp(-1.0, 1.0);
        let r = det::sqrt((1.0 - w * w).max(0.0));
        let phi = 2.0 * std::f64::consts::PI * phi_u;
        let local = [r * det::cos(phi), r * det::sin(phi), w];
        rotate_z_to(mu, local)
    }
}

// The checkpoint replay consumes no additional caller-visible draws.
fn normalize_gamma_replay(mut replay: Stream, alphas: &[f64], out: &mut [f64]) {
    let scale = alphas.iter().fold(1.0f64, |acc, &alpha| acc.min(alpha));
    let mut max_score = f64::NEG_INFINITY;
    for (slot, &alpha) in out.iter_mut().zip(alphas) {
        let score = replay.next_scaled_log_gamma(alpha, scale);
        *slot = score;
        max_score = max_score.max(score);
    }

    let total: f64 = out
        .iter_mut()
        .map(|slot| {
            *slot = det::exp((*slot - max_score) / scale);
            *slot
        })
        .sum();
    for slot in out.iter_mut() {
        *slot /= total;
    }
}

/// Rotate `v` (expressed with +z as the pole) so the pole maps to `mu`:
/// deterministic Rodrigues construction with the antipodal branch handled
/// explicitly (no normalization drift — the output norm equals `v`'s).
fn rotate_z_to(mu: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    let nz = mu[2];
    if nz > 0.999_999 {
        return v;
    }
    if nz < -0.999_999 {
        return [v[0], -v[1], -v[2]]; // 180° about x
    }
    // axis = z × mu (unnormalized), angle θ with cosθ = μz.
    let ax = -mu[1];
    let ay = mu[0];
    let s2 = ax.mul_add(ax, ay * ay); // sin²θ
    let c = nz;
    // Rodrigues with unnormalized axis a (|a| = sinθ):
    // R·v = v·cosθ + (â×v)·sinθ + â(â·v)(1−cosθ)
    //     = v·c + (a×v) + a·(a·v)·(1−c)/sin²θ.
    let cross = [ay * v[2], -ax * v[2], ax.mul_add(v[1], -(ay * v[0]))];
    let adotv = ax.mul_add(v[0], ay * v[1]);
    let k = adotv * (1.0 - c) / s2;
    [
        v[0].mul_add(c, cross[0]) + ax * k,
        v[1].mul_add(c, cross[1]) + ay * k,
        v[2].mul_add(c, cross[2]),
    ]
}

/// Vose alias table for O(1) categorical sampling. Construction is
/// DETERMINISTIC: index-order worklists with lowest-index-first
/// processing (P2 applied to setup) — same weights, same table, bitwise.
#[derive(Debug, Clone)]
pub struct AliasTable {
    prob: Vec<f64>,
    alias: Vec<usize>,
}

impl AliasTable {
    /// Build from non-negative weights (at least one positive). Panics
    /// (structured) otherwise.
    #[must_use]
    pub fn new(weights: &[f64]) -> AliasTable {
        let n = weights.len();
        assert!(n > 0, "alias table needs at least one weight");
        let total: f64 = weights.iter().sum();
        assert!(
            total > 0.0 && weights.iter().all(|&w| w >= 0.0 && w.is_finite()),
            "weights must be non-negative and finite with positive sum"
        );
        let scale = n as f64 / total;
        let scaled: Vec<f64> = weights.iter().map(|&w| w * scale).collect();
        // Index-order worklists (VecDeque-free: two vecs used as stacks,
        // filled in ascending index order, popped from the BACK — the
        // fixed processing order that makes construction bitwise).
        let mut small: Vec<usize> = Vec::new();
        let mut large: Vec<usize> = Vec::new();
        for (i, &s) in scaled.iter().enumerate() {
            if s < 1.0 {
                small.push(i);
            } else {
                large.push(i);
            }
        }
        let mut prob = vec![1.0f64; n];
        let mut alias: Vec<usize> = (0..n).collect();
        let mut rem = scaled;
        while let (Some(&s), Some(&l)) = (small.last(), large.last()) {
            small.pop();
            prob[s] = rem[s];
            alias[s] = l;
            rem[l] = (rem[l] + rem[s]) - 1.0;
            if rem[l] < 1.0 {
                large.pop();
                small.push(l);
            }
        }
        // Leftovers (numerical residue) saturate to probability 1.
        for &l in &large {
            prob[l] = 1.0;
        }
        for &s in &small {
            prob[s] = 1.0;
        }
        AliasTable { prob, alias }
    }

    /// Sample a category: exactly ONE uniform draw, split into bucket
    /// and coin (fixed-consumption contract).
    #[must_use]
    pub fn sample(&self, stream: &mut Stream) -> usize {
        let n = self.prob.len();
        let u = stream.next_f64();
        let scaled = u * n as f64;
        let bucket = (scaled as usize).min(n - 1);
        let coin = scaled - bucket as f64;
        if coin < self.prob[bucket] {
            bucket
        } else {
            self.alias[bucket]
        }
    }

    /// Number of categories.
    #[must_use]
    pub fn len(&self) -> usize {
        self.prob.len()
    }

    /// True if the table has no categories (unreachable by construction).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.prob.is_empty()
    }
}
