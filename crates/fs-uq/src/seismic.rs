//! SEISMIC STACK (bead o5kc, the frame flagship's UQ lane): stochastic
//! ground motion from Kanai–Tajimi-class spectra, the CQC fast path
//! for modal combination, and FRAGILITY via incremental dynamic
//! analysis on a nonlinear (bilinear-hysteretic) SDOF — probability-of-
//! exceedance curves whose confidence machinery lives in
//! [`crate::anytime`].
//!
//! Determinism: ground motions are synthesized by spectral
//! representation with counter-based phases (caller-seeded LCG stream,
//! logical identity — no thread dependence).

/// The Kanai–Tajimi power spectral density: filtered white noise
/// through the soil layer (`wg`, `zg`) at intensity `s0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KanaiTajimi {
    /// White-noise intensity (m²/s³).
    pub s0: f64,
    /// Ground filter frequency (rad/s).
    pub wg: f64,
    /// Ground filter damping ratio.
    pub zg: f64,
}

impl KanaiTajimi {
    /// The one-sided PSD at circular frequency `w`.
    #[must_use]
    pub fn psd(&self, w: f64) -> f64 {
        let r = w / self.wg;
        let num = 1.0 + 4.0 * self.zg * self.zg * r * r;
        let den = (1.0 - r * r).powi(2) + 4.0 * self.zg * self.zg * r * r;
        self.s0 * num / den
    }

    /// Synthesize one ground-acceleration record by spectral
    /// representation: `a(t) = Σ √(2·S(ωk)·Δω)·cos(ωk t + φk)` with
    /// deterministic phases from `seed`.
    #[must_use]
    pub fn synthesize(&self, seed: u64, n_freq: usize, dt: f64, n_steps: usize) -> Vec<f64> {
        let w_max = 4.0 * self.wg.max(1.0) * 2.0;
        #[allow(clippy::cast_precision_loss)]
        let dw = w_max / n_freq as f64;
        // Counter-based phases (splitmix-style hash per frequency line).
        let phase = |k: usize| -> f64 {
            let mut z = seed ^ (0x9e37_79b9_7f4a_7c15u64.wrapping_mul(k as u64 + 1));
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^= z >> 31;
            (z >> 11) as f64 / (1u64 << 53) as f64 * std::f64::consts::TAU
        };
        let lines: Vec<(f64, f64, f64)> = (0..n_freq)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let w = (k as f64 + 0.5) * dw;
                ((2.0 * self.psd(w) * dw).sqrt(), w, phase(k))
            })
            .collect();
        (0..n_steps)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let t = i as f64 * dt;
                lines
                    .iter()
                    .map(|&(amp, w, ph)| amp * (w * t + ph).cos())
                    .sum()
            })
            .collect()
    }
}

/// Peak absolute displacement of a LINEAR SDOF (`wn`, `zeta`) under a
/// ground record (Newmark average-acceleration; the response-spectrum
/// ordinate).
#[must_use]
pub fn sdof_peak(record: &[f64], dt: f64, wn: f64, zeta: f64) -> f64 {
    let (mut u, mut v) = (0.0f64, 0.0f64);
    let mut a = -record.first().copied().unwrap_or(0.0);
    let (k, c) = (wn * wn, 2.0 * zeta * wn);
    let mut peak = 0.0f64;
    for &ag in &record[1..] {
        // Newmark beta = 1/4, gamma = 1/2 on u'' + c u' + k u = -ag.
        let rhs = -ag - c * (v + 0.5 * dt * a) - k * (u + dt * v + 0.25 * dt * dt * a);
        let denom = 1.0 + 0.5 * dt * c + 0.25 * dt * dt * k;
        let a_new = rhs / denom;
        u += dt * v + 0.25 * dt * dt * (a + a_new);
        v += 0.5 * dt * (a + a_new);
        a = a_new;
        peak = peak.max(u.abs());
    }
    peak
}

/// CQC (complete quadratic combination) of modal peaks with the
/// Der Kiureghian correlation coefficients — the response-spectrum
/// fast path. Falls back to SRSS exactly when modes are uncorrelated.
#[must_use]
pub fn cqc(peaks: &[f64], freqs: &[f64], damps: &[f64]) -> f64 {
    let n = peaks.len();
    let mut acc = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let (zi, zj) = (damps[i], damps[j]);
            let r = freqs[j] / freqs[i];
            let num = 8.0 * (zi * zj).sqrt() * (zi + r * zj) * r.powf(1.5);
            let den = (1.0 - r * r).powi(2)
                + 4.0 * zi * zj * r * (1.0 + r * r)
                + 4.0 * (zi * zi + zj * zj) * r * r;
            let rho = num / den;
            acc += rho * peaks[i] * peaks[j];
        }
    }
    acc.sqrt()
}

/// SRSS combination (the uncorrelated-modes baseline CQC reduces to).
#[must_use]
pub fn srss(peaks: &[f64]) -> f64 {
    peaks.iter().map(|p| p * p).sum::<f64>().sqrt()
}

/// Peak DUCTILITY of a BILINEAR-hysteretic SDOF (yield displacement
/// `uy`, post-yield stiffness ratio `alpha`) under a scaled record —
/// the nonlinear time-history IDA workhorse.
#[must_use]
pub fn bilinear_peak_ductility(
    record: &[f64],
    dt: f64,
    wn: f64,
    zeta: f64,
    uy: f64,
    alpha: f64,
    scale: f64,
) -> f64 {
    let (k0, c) = (wn * wn, 2.0 * zeta * wn);
    let (mut u, mut v) = (0.0f64, 0.0f64);
    // Hysteretic state: plastic offset of the yield surface.
    let mut up = 0.0f64;
    let mut peak = 0.0f64;
    // Explicit central-difference march with the bilinear force
    // (simple, robust for the fixture's dt).
    let force = |u: f64, up: &mut f64| -> f64 {
        let elastic = u - *up;
        if elastic > uy {
            *up += elastic - uy;
        } else if elastic < -uy {
            *up += elastic + uy;
        }
        k0 * (alpha * u + (1.0 - alpha) * (u - *up))
    };
    for &ag in record {
        let f = force(u, &mut up);
        let acc = -scale * ag - c * v - f;
        v += dt * acc;
        u += dt * v;
        peak = peak.max(u.abs());
    }
    peak / uy
}

/// One fragility point: the empirical exceedance probability at an
/// intensity level over a motion suite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragilityPoint {
    /// Intensity measure (record scale factor).
    pub im: f64,
    /// Empirical exceedance probability.
    pub p: f64,
    /// Suite size behind it.
    pub n: usize,
}

/// INCREMENTAL DYNAMIC ANALYSIS: scale every motion in the suite
/// through the IM ladder, march the bilinear SDOF, and report the
/// exceedance fraction (`ductility > threshold`) per level. The
/// confidence machinery (anytime CS on each point) is
/// [`crate::anytime::estimate_probability_anytime`]'s job.
#[must_use]
#[allow(clippy::too_many_arguments, clippy::cast_precision_loss)]
pub fn ida_fragility(
    kt: &KanaiTajimi,
    suite_seeds: &[u64],
    im_ladder: &[f64],
    wn: f64,
    zeta: f64,
    uy: f64,
    alpha: f64,
    ductility_limit: f64,
) -> Vec<FragilityPoint> {
    let (dt, n_steps, n_freq) = (0.01, 1500, 96);
    let records: Vec<Vec<f64>> = suite_seeds
        .iter()
        .map(|&s| kt.synthesize(s, n_freq, dt, n_steps))
        .collect();
    im_ladder
        .iter()
        .map(|&im| {
            let hits = records
                .iter()
                .filter(|r| {
                    bilinear_peak_ductility(r, dt, wn, zeta, uy, alpha, im) > ductility_limit
                })
                .count();
            FragilityPoint {
                im,
                p: hits as f64 / records.len().max(1) as f64,
                n: records.len(),
            }
        })
        .collect()
}
