//! SEISMIC STACK (bead o5kc, the frame flagship's UQ lane): stochastic
//! ground motion from Kanai–Tajimi-class spectra, the CQC fast path
//! for modal combination, and FRAGILITY via incremental dynamic
//! analysis on a nonlinear (bilinear-hysteretic) SDOF — probability-of-
//! exceedance curves whose confidence machinery lives in
//! [`crate::anytime`].
//!
//! Determinism: ground motions are synthesized by spectral
//! representation with counter-based phases (caller-seeded frequency-line identities,
//! logical identity — no thread dependence).

/// Recorded acceleration histories and physically normalized motion suites.
pub mod recorded;
pub use recorded::{
    AccelerationUnit, BilinearIdaModel, GroundMotionRecord, GroundMotionRecordError, IdaIntensity,
    RecordedIda, ida_fragility_from_records,
};

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

/// Sampling and physical bandwidth of a finite-line ground-motion realization.
///
/// The target PSD is one-sided in angular frequency, with
/// variance = integral S(omega) d(omega). Midpoint lines partition
/// [0, omega_max]; each line contributes S(omega_k) * delta_omega to
/// ensemble variance. Power outside this band is omitted, never aliased or
/// redistributed. Increase n_freq to resolve a narrow soil resonance.
///
/// The frequencies and seeded phases do not depend on dt_s or n_steps.
/// Keeping omega_max and n_freq fixed therefore preserves the same continuous
/// excitation when changing the sample clock or extending a record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundMotionSampling {
    /// Sample interval [s].
    pub dt_s: f64,
    /// Number of samples; N samples span N-1 response-integration intervals.
    pub n_steps: usize,
    /// Number of midpoint frequency lines.
    pub n_freq: usize,
    /// Physical low-pass cutoff [rad/s], at most pi/dt_s.
    pub omega_max: f64,
}

impl GroundMotionSampling {
    #[allow(clippy::cast_precision_loss)]
    fn validate(&self) {
        assert!(
            self.dt_s.is_finite() && self.dt_s > 0.0,
            "ground-motion sample interval must be positive and finite"
        );
        assert!(
            self.n_freq > 0,
            "ground-motion frequency count must be positive"
        );
        assert!(
            self.omega_max.is_finite() && self.omega_max > 0.0,
            "ground-motion bandwidth must be positive and finite"
        );
        assert!(
            self.omega_max <= std::f64::consts::PI / self.dt_s,
            "ground-motion bandwidth exceeds the sampling Nyquist frequency"
        );
        let dw = self.omega_max / self.n_freq as f64;
        assert!(
            dw.is_finite() && dw > 0.0,
            "ground-motion frequency spacing must be representable"
        );
        let final_time = self.n_steps.saturating_sub(1) as f64 * self.dt_s;
        assert!(
            final_time.is_finite() && (self.omega_max * final_time).is_finite(),
            "ground-motion sample times and phases must be representable"
        );
    }
}

impl KanaiTajimi {
    fn validate(&self) {
        assert!(
            self.s0.is_finite() && self.s0 >= 0.0,
            "ground-motion PSD intensity must be nonnegative and finite"
        );
        assert!(
            self.wg.is_finite() && self.wg > 0.0,
            "ground-filter frequency must be positive and finite"
        );
        assert!(
            self.zg.is_finite() && self.zg > 0.0,
            "stationary ground-filter damping must be positive and finite"
        );
    }

    /// One-sided acceleration PSD per rad/s at nonnegative angular frequency w.
    ///
    /// Its integral over an angular-frequency band is acceleration variance;
    /// there is no additional factor of 2*pi. A density per Hz would instead
    /// be 2*pi*S(2*pi*f).
    ///
    /// # Panics
    /// Invalid model parameters or frequency, or an unrepresentable PSD.
    #[must_use]
    pub fn psd(&self, w: f64) -> f64 {
        self.validate();
        assert!(
            w.is_finite() && w >= 0.0,
            "one-sided PSD frequency must be nonnegative and finite"
        );
        let r = w / self.wg;
        let num = 1.0 + 4.0 * self.zg * self.zg * r * r;
        let den = (1.0 - r * r).powi(2) + 4.0 * self.zg * self.zg * r * r;
        let value = self.s0 * num / den;
        assert!(
            value.is_finite() && value >= 0.0,
            "ground-motion PSD must be representable"
        );
        value
    }

    /// Return the legacy convenience API's actual sampling specification.
    ///
    /// The physical cutoff is min(8*wg, pi/dt). This is a sampled, band-limited
    /// approximation to the Kanai-Tajimi process: the omitted high-frequency
    /// power is not restored by rescaling the retained band. Use an explicit
    /// GroundMotionSampling to keep bandwidth fixed across clock refinements.
    ///
    /// # Panics
    /// Invalid model or sampling parameters.
    #[must_use]
    pub fn default_sampling(&self, n_freq: usize, dt: f64, n_steps: usize) -> GroundMotionSampling {
        self.validate();
        let sampling = GroundMotionSampling {
            dt_s: dt,
            n_steps,
            n_freq,
            omega_max: (8.0 * self.wg).min(std::f64::consts::PI / dt),
        };
        sampling.validate();
        sampling
    }

    /// Synthesize a record using the documented default low-pass bandwidth.
    ///
    /// See [Self::default_sampling] for the actual cutoff and
    /// [Self::synthesize_with_sampling] for explicit bandwidth control.
    ///
    /// # Panics
    /// Invalid model/sampling parameters or unrepresentable arithmetic.
    #[must_use]
    pub fn synthesize(&self, seed: u64, n_freq: usize, dt: f64, n_steps: usize) -> Vec<f64> {
        self.synthesize_with_sampling(seed, self.default_sampling(n_freq, dt, n_steps))
    }

    /// Synthesize an explicitly band-limited ground-acceleration record.
    ///
    /// a(t) = sum sqrt(2*S(omega_k)*delta_omega)*cos(omega_k*t + phi_k).
    /// Counter-addressed uniform phases make this finite random-phase
    /// approximation reproducible, with covariance
    /// sum S(omega_k)*delta_omega*cos(omega_k*tau). Every line lies strictly
    /// below Nyquist; this method refuses an incompatible clock rather than
    /// folding high-frequency excitation into lower structural modes.
    ///
    /// # Panics
    /// Invalid model/sampling parameters or unrepresentable arithmetic.
    #[must_use]
    pub fn synthesize_with_sampling(&self, seed: u64, sampling: GroundMotionSampling) -> Vec<f64> {
        let lines = self.spectral_lines(sampling);
        synthesize_lines(seed, sampling, &lines)
    }

    #[allow(clippy::cast_precision_loss)]
    fn spectral_lines(&self, sampling: GroundMotionSampling) -> Vec<(f64, f64)> {
        self.validate();
        sampling.validate();
        let dw = sampling.omega_max / sampling.n_freq as f64;
        let root_dw = fs_math::det::sqrt(dw);
        (0..sampling.n_freq)
            .map(|k| {
                let w = (k as f64 + 0.5) * dw;
                assert!(
                    w.is_finite() && w > 0.0 && w * sampling.dt_s < std::f64::consts::PI,
                    "ground-motion frequency line must lie strictly below Nyquist"
                );
                // Separate square roots avoid overflowing or underflowing
                // the dimensional product S*dw before taking its square root.
                let amp = std::f64::consts::SQRT_2 * fs_math::det::sqrt(self.psd(w)) * root_dw;
                assert!(
                    amp.is_finite(),
                    "ground-motion amplitude must be representable"
                );
                (amp, w)
            })
            .collect()
    }
}

#[allow(clippy::cast_precision_loss)]
fn synthesize_lines(seed: u64, sampling: GroundMotionSampling, lines: &[(f64, f64)]) -> Vec<f64> {
    // Counter-based phases (splitmix-style hash per frequency line). Neither
    // sample count nor clock participates in a line's random identity.
    let phases: Vec<_> = (0..lines.len())
        .map(|k| {
            let mut z = seed ^ (0x9e37_79b9_7f4a_7c15u64.wrapping_mul((k as u64).wrapping_add(1)));
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^= z >> 31;
            (z >> 11) as f64 / (1u64 << 53) as f64 * std::f64::consts::TAU
        })
        .collect();
    (0..sampling.n_steps)
        .map(|i| {
            let t = i as f64 * sampling.dt_s;
            let value: f64 = lines
                .iter()
                .zip(&phases)
                .map(|(&(amp, w), &phase)| amp * fs_math::det::cos(w * t + phase))
                .sum();
            assert!(
                value.is_finite(),
                "ground-acceleration sample must be representable"
            );
            value
        })
        .collect()
}

/// Maximum absolute RELATIVE displacement [m] of a linear SDOF under a
/// ground-acceleration record [m/s²], using Newmark average acceleration.
/// `wn` is the undamped angular frequency [rad/s], `zeta` its damping ratio.
/// The oscillator starts at rest at the first sample; N samples span N−1
/// intervals. Empty and single-sample records therefore have zero response.
///
/// # Panics
/// Nonfinite samples, nonpositive/nonfinite `dt` or `wn`, negative/nonfinite
/// damping, or an unrepresentable response refuse rather than entering a peak.
#[must_use]
pub fn sdof_peak(record: &[f64], dt: f64, wn: f64, zeta: f64) -> f64 {
    let (k, c) = sdof_coefficients(dt, wn, zeta);
    let (mut u, mut v) = (0.0f64, 0.0f64);
    let mut a = -record.first().copied().unwrap_or(0.0);
    assert!(a.is_finite(), "ground acceleration must be finite");
    let denom = 1.0 + 0.5 * dt * c + 0.25 * dt * dt * k;
    assert!(denom.is_finite(), "SDOF effective stiffness must be finite");
    let mut peak = 0.0f64;
    for &ag in record.iter().skip(1) {
        assert!(ag.is_finite(), "ground acceleration must be finite");
        // Newmark beta = 1/4, gamma = 1/2 on u'' + c u' + k u = -ag.
        let rhs = -ag - c * (v + 0.5 * dt * a) - k * (u + dt * v + 0.25 * dt * dt * a);
        let a_new = rhs / denom;
        u += dt * v + 0.25 * dt * dt * (a + a_new);
        v += 0.5 * dt * (a + a_new);
        a = a_new;
        assert!(
            u.is_finite() && v.is_finite() && a.is_finite(),
            "SDOF response must be finite"
        );
        peak = peak.max(u.abs());
    }
    peak
}

/// CQC (complete quadratic combination) of modal peaks with the
/// Der Kiureghian correlation coefficients — the response-spectrum
/// fast path. Falls back to SRSS exactly when modes are uncorrelated.
/// Peaks retain their modal signs. Zero-damped distinct frequencies have
/// zero correlation; equal-frequency, zero-damped modes are identical linear
/// responses and combine coherently. This limiting convention is not a
/// stationary-response claim for undamped white-noise excitation.
///
/// # Panics
/// Array lengths must agree, frequencies must be positive and finite, damping
/// ratios nonnegative and finite, and modal peaks finite.
#[must_use]
pub fn cqc(peaks: &[f64], freqs: &[f64], damps: &[f64]) -> f64 {
    let n = peaks.len();
    assert_eq!(freqs.len(), n, "CQC frequency count must match modal peaks");
    assert_eq!(damps.len(), n, "CQC damping count must match modal peaks");
    let mut scale = 0.0f64;
    for i in 0..n {
        assert!(peaks[i].is_finite(), "CQC modal peaks must be finite");
        assert!(
            freqs[i].is_finite() && freqs[i] > 0.0,
            "CQC frequencies must be positive and finite"
        );
        assert!(
            damps[i].is_finite() && damps[i] >= 0.0,
            "CQC damping must be nonnegative and finite"
        );
        scale = scale.max(peaks[i].abs());
    }
    if scale == 0.0 {
        return 0.0;
    }
    let mut acc = 0.0f64;
    for i in 0..n {
        let pi = peaks[i] / scale;
        // The diagonal is exactly one, including a zero-damped mode; the
        // literal coefficient formula is 0/0 there. Symmetry halves the work.
        acc += pi * pi;
        for j in 0..i {
            let rho = cqc_correlation(freqs[i], freqs[j], damps[i], damps[j]);
            acc += 2.0 * rho * pi * (peaks[j] / scale);
        }
    }
    // Coincident modes of opposite sign can cancel exactly. The correlation
    // form is nonnegative; only summation roundoff may cross zero.
    scale * acc.max(0.0).sqrt()
}

fn cqc_correlation(wi: f64, wj: f64, zi: f64, zj: f64) -> f64 {
    if zi == 0.0 || zj == 0.0 {
        return if zi == zj && wi == wj { 1.0 } else { 0.0 };
    }
    // Choose r <= 1 so that wide modal bands cannot overflow r^4. Scale
    // damping before products, preserving the unequal-damping formula.
    let (r, zi, zj) = if wj <= wi {
        (wj / wi, zi, zj)
    } else {
        (wi / wj, zj, zi)
    };
    if r == 0.0 {
        return 0.0;
    }
    let damping_scale = zi.max(zj);
    let zi = zi / damping_scale;
    let zj = zj / damping_scale;
    if r == 1.0 {
        return 2.0 * zi.sqrt() * zj.sqrt() / (zi + zj);
    }
    let detuning = ((1.0 - r) * (1.0 + r)) / damping_scale;
    let num = 8.0 * zi.sqrt() * zj.sqrt() * (zi + r * zj) * r * r.sqrt();
    let den =
        detuning * detuning + 4.0 * zi * zj * r * (1.0 + r * r) + 4.0 * (zi * zi + zj * zj) * r * r;
    let rho = num / den;
    assert!(rho.is_finite(), "CQC correlation must be representable");
    rho.clamp(0.0, 1.0)
}

fn sdof_coefficients(dt: f64, wn: f64, zeta: f64) -> (f64, f64) {
    assert!(
        dt.is_finite() && dt > 0.0,
        "SDOF time step must be positive and finite"
    );
    assert!(
        wn.is_finite() && wn > 0.0,
        "SDOF frequency must be positive and finite"
    );
    assert!(
        zeta.is_finite() && zeta >= 0.0,
        "SDOF damping must be nonnegative and finite"
    );
    let (k, c) = (wn * wn, 2.0 * zeta * wn);
    assert!(
        k.is_finite() && c.is_finite(),
        "SDOF stiffness and damping must be finite"
    );
    (k, c)
}

/// SRSS combination (the uncorrelated-modes baseline CQC reduces to).
#[must_use]
pub fn srss(peaks: &[f64]) -> f64 {
    peaks.iter().map(|p| p * p).sum::<f64>().sqrt()
}

/// Peak DUCTILITY of a BILINEAR-hysteretic SDOF (yield displacement
/// `uy`, post-yield stiffness ratio `alpha`) under a scaled record —
/// the nonlinear time-history IDA workhorse.
///
/// Uses implicit Newmark average acceleration, with the same sample clock
/// and initial rest state as [`sdof_peak`]. The restoring force is a parallel
/// elastic spring and an elastic-perfectly-plastic spring:
/// `R = alpha*k0*u + (1-alpha)*k0*clamp(u-up, -uy, uy)`.
/// The monotone, piecewise-linear equilibrium is solved exactly in at most
/// two scalar solves; plastic history is committed only after branch selection.
/// This removes the explicit method's stability restriction, but does not
/// remove the need to resolve excitation, oscillation and yield reversals in time.
/// See the Newmark kinematics in
/// <https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/integrator/Newmark.html>.
///
/// # Panics
/// As for [`sdof_peak`], plus nonpositive/nonfinite yield displacement,
/// `alpha` outside `[0, 1]`, nonfinite record scale, or unrepresentable response.
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
    let (k0, c) = sdof_coefficients(dt, wn, zeta);
    assert!(
        uy.is_finite() && uy > 0.0,
        "yield displacement must be positive and finite"
    );
    assert!(
        (0.0..=1.0).contains(&alpha),
        "post-yield stiffness ratio must be in [0, 1]"
    );
    assert!(scale.is_finite(), "ground-motion scale must be finite");
    let (mut u, mut v) = (0.0f64, 0.0f64);
    let mut a = -scale * record.first().copied().unwrap_or(0.0);
    assert!(a.is_finite(), "scaled ground acceleration must be finite");
    // Hysteretic state: plastic offset of the yield surface.
    let mut up = 0.0f64;
    let beta_dt2 = 0.25 * dt * dt;
    let gamma_dt = 0.5 * dt;
    let k_post = alpha * k0;
    let k_plastic = (1.0 - alpha) * k0;
    let elastic_denom = 1.0 + gamma_dt * c + beta_dt2 * k0;
    let yielded_denom = 1.0 + gamma_dt * c + beta_dt2 * k_post;
    assert!(
        elastic_denom.is_finite(),
        "SDOF effective stiffness must be finite"
    );
    let mut peak = 0.0f64;
    for &ag in record.iter().skip(1) {
        let load = -scale * ag;
        assert!(
            load.is_finite(),
            "scaled ground acceleration must be finite"
        );
        let u_predict = u + dt * v + beta_dt2 * a;
        let v_predict = v + gamma_dt * a;
        let rhs = load - c * v_predict;
        let mut a_next = (rhs - k0 * u_predict + k_plastic * up) / elastic_denom;
        let trial = u_predict + beta_dt2 * a_next - up;
        if trial.abs() > uy && alpha < 1.0 {
            let elastic = uy.copysign(trial);
            a_next = (rhs - k_post * u_predict - k_plastic * elastic) / yielded_denom;
            // The accepted displacement returns the plastic spring to its
            // yield surface; failed elastic trials never mutate the history.
            up = u_predict + beta_dt2 * a_next - elastic;
        }
        u = u_predict + beta_dt2 * a_next;
        v = v_predict + gamma_dt * a_next;
        a = a_next;
        assert!(
            u.is_finite() && v.is_finite() && a.is_finite() && up.is_finite(),
            "SDOF response must be finite"
        );
        peak = peak.max(u.abs());
    }
    peak / uy
}

/// One fragility point: the empirical exceedance probability at an
/// intensity level over a motion suite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragilityPoint {
    /// Intensity ordinate; a record scale factor for the synthetic-suite APIs.
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
    ida_fragility_with_sampling(
        kt,
        kt.default_sampling(96, 0.01, 1500),
        suite_seeds,
        im_ladder,
        wn,
        zeta,
        uy,
        alpha,
        ductility_limit,
    )
}

/// Incremental dynamic analysis with caller-controlled excitation sampling.
///
/// The same physical bandwidth and frequency lines can be retained while
/// refining dt_s; the oscillator uses that exact sample interval. Spectral
/// amplitudes are prepared once for the entire seeded suite, and every
/// intensity level reuses the same realized records. The confidence machinery
/// remains [crate::anytime::estimate_probability_anytime]'s responsibility.
///
/// # Panics
/// Invalid sampling/model parameters or unrepresentable synthesis/response.
#[must_use]
#[allow(clippy::too_many_arguments, clippy::cast_precision_loss)]
pub fn ida_fragility_with_sampling(
    kt: &KanaiTajimi,
    sampling: GroundMotionSampling,
    suite_seeds: &[u64],
    im_ladder: &[f64],
    wn: f64,
    zeta: f64,
    uy: f64,
    alpha: f64,
    ductility_limit: f64,
) -> Vec<FragilityPoint> {
    let lines = kt.spectral_lines(sampling);
    let records: Vec<Vec<f64>> = suite_seeds
        .iter()
        .map(|&s| synthesize_lines(s, sampling, &lines))
        .collect();
    im_ladder
        .iter()
        .map(|&im| {
            let hits = records
                .iter()
                .filter(|r| {
                    bilinear_peak_ductility(r, sampling.dt_s, wn, zeta, uy, alpha, im)
                        > ductility_limit
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

#[cfg(test)]
mod tests {
    use super::{
        GroundMotionSampling, KanaiTajimi, bilinear_peak_ductility, cqc,
        ida_fragility_with_sampling, sdof_peak,
    };

    #[test]
    fn g0_empty_and_single_sample_records_do_not_advance_time() {
        for record in [&[][..], &[3.0][..]] {
            assert_eq!(sdof_peak(record, 0.01, 12.0, 0.05), 0.0);
            assert_eq!(
                bilinear_peak_ductility(record, 0.01, 12.0, 0.05, 0.01, 0.1, 2.0),
                0.0
            );
        }
    }

    #[test]
    fn g1_bilinear_elastic_limit_matches_the_linear_response_spectrum() {
        // Nonzero first acceleration exposes an extra initial step or an
        // inconsistent acceleration initialization. Stiff cases exceed the
        // old explicit method's stability limit by an order of magnitude.
        let record: Vec<_> = (0..1000)
            .map(|i| {
                2.0 * fs_math::det::cos(0.017 * f64::from(i))
                    + 0.5 * fs_math::det::sin(0.11 * f64::from(i))
            })
            .collect();
        let scaled: Vec<_> = record.iter().map(|ag| -1.7 * ag).collect();
        for wn in [2.0, 20.0, 1000.0] {
            for zeta in [0.0, 0.05, 1.0, 3.0] {
                let want = sdof_peak(&scaled, 0.01, wn, zeta);
                for (alpha, uy) in [(1.0, 0.001), (0.1, 100.0)] {
                    let got =
                        bilinear_peak_ductility(&record, 0.01, wn, zeta, uy, alpha, -1.7) * uy;
                    assert!(
                        (got - want).abs() <= 2.0e-11 * want.max(1.0e-12),
                        "wn={wn} zeta={zeta} alpha={alpha}: {got} vs {want}"
                    );
                }
            }
        }
    }

    #[test]
    fn g1_short_period_oscillator_remains_bounded_under_constant_acceleration() {
        let record = vec![1.0; 2000];
        let wn = 500.0;
        let uy = 0.001;
        let peak = bilinear_peak_ductility(&record, 0.01, wn, 0.0, uy, 0.1, 1.0) * uy;
        // The initial state is rest; undamped constant loading oscillates
        // about -1/wn² with amplitude 1/wn². Average acceleration preserves
        // that quadratic energy even at wn*dt=5; phase still requires refinement.
        assert!(peak.is_finite() && peak > 1.9 / (wn * wn));
        assert!(peak <= 2.0 / (wn * wn) * (1.0 + 1.0e-10));
    }

    #[test]
    fn g1_yielded_step_solves_the_accepted_restoring_force() {
        // One interval from zero load into yield. beta*h²=.01, k0=100,
        // k_post=10, plastic yield force=.09, so a1=(1-.09)/(1+.1).
        let got = bilinear_peak_ductility(&[0.0, -1.0], 0.2, 10.0, 0.0, 0.001, 0.1, 1.0);
        let want = 0.01 * (1.0 - 0.09) / 1.1 / 0.001;
        assert!((got - want).abs() < 2.0e-14);
        // Reflection reverses the plastic branch but not the absolute peak.
        let reflected = bilinear_peak_ductility(&[0.0, 1.0], 0.2, 10.0, 0.0, 0.001, 0.1, 1.0);
        assert_eq!(got.to_bits(), reflected.to_bits());
    }

    fn cyclic_ground(t: f64) -> f64 {
        3.0 * fs_math::det::sin(7.0 * t) + 1.1 * fs_math::det::sin(19.0 * t)
    }

    fn fine_explicit_hysteresis_peak(n: u32) -> f64 {
        // Independent first-order reference: resolve elastic/plastic reversals
        // with 300k/600k explicit steps, no Newmark prediction or branch solve.
        let dt = 3.0 / f64::from(n);
        let (mut u, mut v, mut up, mut peak) = (0.0_f64, 0.0, 0.0, 0.0_f64);
        for i in 0..n {
            let elastic = (u - up).clamp(-0.01, 0.01);
            up = u - elastic;
            let force = 144.0 * (0.1 * u + 0.9 * elastic);
            let a = -cyclic_ground(f64::from(i) * dt) - 1.2 * v - force;
            v += dt * a;
            u += dt * v;
            peak = peak.max(u.abs());
        }
        peak / 0.01
    }

    #[test]
    fn g1_cyclic_plastic_response_converges_to_independent_fine_history() {
        let reference = fine_explicit_hysteresis_peak(600_000);
        let reference_coarse = fine_explicit_hysteresis_peak(300_000);
        assert!((reference - reference_coarse).abs() < 3.0e-5);
        assert!(
            reference > 10.0,
            "the comparison must exercise substantial yield"
        );
        let mut previous_error = f64::INFINITY;
        for n in [300, 600, 1200] {
            let dt = 3.0 / f64::from(n);
            let record: Vec<_> = (0..=n).map(|i| cyclic_ground(f64::from(i) * dt)).collect();
            let got = bilinear_peak_ductility(&record, dt, 12.0, 0.05, 0.01, 0.1, 1.0);
            let error = (got - reference).abs();
            assert!(
                error < 0.4 * previous_error,
                "n={n}: refinement error {error} vs {previous_error}"
            );
            previous_error = error;
        }
        assert!(previous_error / reference < 1.2e-4);
    }

    #[test]
    fn g0_cqc_handles_undamped_modes_and_coherent_degeneracy() {
        assert_eq!(cqc(&[], &[], &[]), 0.0);
        assert_eq!(cqc(&[-3.0], &[10.0], &[0.0]), 3.0);
        assert_eq!(cqc(&[3.0, 4.0], &[10.0, 20.0], &[0.0, 0.0]), 5.0);
        assert_eq!(cqc(&[3.0, 4.0], &[10.0, 10.0], &[0.0, 0.0]), 7.0);
        assert_eq!(cqc(&[2.0, -1.0], &[10.0, 10.0], &[0.0, 0.0]), 1.0);
        assert_eq!(cqc(&[1.0e200, -1.0e200], &[10.0, 10.0], &[0.0, 0.0]), 0.0);
        assert_eq!(cqc(&[3.0, 4.0], &[1.0e-200, 1.0e200], &[0.05, 0.05]), 5.0);
    }

    #[test]
    fn g3_cqc_unequal_damping_preserves_symmetry_units_and_zero_loss_limit() {
        let peaks = [1.0, -0.6];
        let (zi, zj) = (0.02_f64, 0.08_f64);
        let r = 18.0_f64 / 17.0;
        // Literal published expression, evaluated independently without the
        // implementation's frequency ordering or damping scaling.
        let rho = 8.0 * (zi * zj).sqrt() * (zi + r * zj) * r * r.sqrt()
            / ((1.0 - r * r).powi(2)
                + 4.0 * zi * zj * r * (1.0 + r * r)
                + 4.0 * (zi * zi + zj * zj) * r * r);
        let want = (1.36 - 1.2 * rho).sqrt();
        for scale in [1.0e-100, 1.0, 1.0e100] {
            let got = cqc(&peaks, &[17.0 * scale, 18.0 * scale], &[zi, zj]);
            let reversed = cqc(&[-0.6, 1.0], &[18.0 * scale, 17.0 * scale], &[zj, zi]);
            assert!((got - want).abs() < 5.0e-15);
            assert_eq!(got.to_bits(), reversed.to_bits());
        }
        let lossless = cqc(&peaks, &[17.0, 18.0], &[0.0, 0.0]);
        let small_loss = cqc(&peaks, &[17.0, 18.0], &[1.0e-12, 1.0e-12]);
        assert_eq!(lossless.to_bits(), small_loss.to_bits());
        let degenerate = cqc(&[1.0, 1.0], &[17.0, 17.0], &[zi, zj]);
        assert!((degenerate - (2.0 + 4.0 * (zi * zj).sqrt() / (zi + zj)).sqrt()).abs() < 1.0e-15);
    }

    #[test]
    fn g0_nonphysical_inputs_cannot_turn_into_plausible_peak_values() {
        assert!(std::panic::catch_unwind(|| sdof_peak(&[], 0.0, 10.0, 0.05)).is_err());
        assert!(
            std::panic::catch_unwind(|| sdof_peak(&[0.0, f64::NAN], 0.01, 10.0, 0.05)).is_err()
        );
        assert!(
            std::panic::catch_unwind(|| bilinear_peak_ductility(
                &[0.0],
                0.01,
                10.0,
                -0.1,
                0.01,
                0.1,
                1.0
            ))
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| bilinear_peak_ductility(
                &[0.0],
                0.01,
                10.0,
                0.05,
                0.01,
                -0.1,
                1.0
            ))
            .is_err()
        );
        assert!(std::panic::catch_unwind(|| cqc(&[1.0], &[], &[0.05])).is_err());
        assert!(std::panic::catch_unwind(|| cqc(&[1.0], &[0.0], &[0.05])).is_err());
    }
    fn spectral_model() -> KanaiTajimi {
        KanaiTajimi {
            s0: 0.03,
            wg: 15.0,
            zg: 0.6,
        }
    }

    #[test]
    fn g1_sampled_spectrum_has_the_declared_one_sided_power_and_covariance() {
        let model = spectral_model();
        let sampling = model.default_sampling(16, 0.05, 64);
        let tau = std::f64::consts::TAU;
        assert_eq!(sampling.omega_max, std::f64::consts::PI / sampling.dt_s);
        let dw = sampling.omega_max / sampling.n_freq as f64;
        // These midpoint lines are the odd Fourier bins of a 4*n_freq
        // sample period. Orthogonality gives an exact finite-record oracle,
        // independent of a statistical tolerance on an ensemble of seeds.
        for seed in [0, 7, 101] {
            let record = model.synthesize_with_sampling(seed, sampling);
            let n = record.len() as f64;
            let mut expected_variance = 0.0;
            let mut expected_covariance = 0.0;
            for k in 0..sampling.n_freq {
                let omega = (k as f64 + 0.5) * dw;
                // Independent dimensional form of the target soil transfer.
                let wg2 = model.wg * model.wg;
                let w2 = omega * omega;
                let damping = 4.0 * model.zg * model.zg * wg2 * w2;
                let density = model.s0 * (wg2 * wg2 + damping) / ((wg2 - w2).powi(2) + damping);
                let mut real = 0.0;
                let mut imag = 0.0;
                for (i, value) in record.iter().enumerate() {
                    let phase = tau * (2 * k + 1) as f64 * i as f64 / n;
                    real += value * phase.cos();
                    imag -= value * phase.sin();
                }
                let line_power = 2.0 * ((real / n).powi(2) + (imag / n).powi(2));
                // Convert the one-sided periodogram from power/bin to /Hz.
                let density_per_hz = line_power / (dw / tau);
                assert!(
                    (density_per_hz - tau * density).abs() < 2e-11 * tau * density,
                    "seed={seed}, line={k}: PSD normalization"
                );
                expected_variance += density * dw;
                expected_covariance += density * dw * (3.0 * omega * sampling.dt_s).cos();
            }
            let mean = record.iter().sum::<f64>() / n;
            let variance = record.iter().map(|x| x * x).sum::<f64>() / n;
            let covariance = record
                .iter()
                .enumerate()
                .map(|(i, x)| x * record[(i + 3) % record.len()])
                .sum::<f64>()
                / n;
            assert!(mean.abs() < 1e-13);
            assert!((variance - expected_variance).abs() < 2e-12 * expected_variance);
            assert!((covariance - expected_covariance).abs() < 2e-12 * expected_variance);
            assert_eq!(record, model.synthesize(seed, 16, 0.05, 64));
        }
    }

    #[test]
    fn g1_fixed_bandwidth_preserves_excitation_under_clock_refinement_and_extension() {
        let model = spectral_model();
        let sampling = GroundMotionSampling {
            dt_s: 0.02,
            n_steps: 257,
            n_freq: 96,
            omega_max: 100.0,
        };
        let coarse = model.synthesize_with_sampling(37, sampling);
        let fine = model.synthesize_with_sampling(
            37,
            GroundMotionSampling {
                dt_s: sampling.dt_s / 2.0,
                n_steps: 2 * (sampling.n_steps - 1) + 1,
                ..sampling
            },
        );
        let extended = model.synthesize_with_sampling(
            37,
            GroundMotionSampling {
                n_steps: 313,
                ..sampling
            },
        );
        for i in 0..coarse.len() {
            assert_eq!(coarse[i].to_bits(), fine[2 * i].to_bits());
            assert_eq!(coarse[i].to_bits(), extended[i].to_bits());
        }
    }

    #[test]
    fn g3_ground_motion_preserves_physical_time_frequency_and_intensity_scaling() {
        let model = spectral_model();
        let reference = model.synthesize(29, 48, 0.01, 97);
        for scale in [1.0 / 64.0, 1.0 / 8.0, 8.0, 64.0] {
            // omega -> scale*omega and t -> t/scale preserve phase.
            // S_omega -> S_omega/scale preserves acceleration variance.
            // In particular wg<1 must not trigger an arbitrary 1-rad/s floor.
            let scaled_model = KanaiTajimi {
                wg: scale * model.wg,
                s0: model.s0 / scale,
                ..model
            };
            let actual = scaled_model.synthesize(29, 48, 0.01 / scale, 97);
            for (a, b) in actual.iter().zip(&reference) {
                assert!(
                    (a - b).abs() < 2e-13 * b.abs().max(1.0),
                    "scale={scale}: {a} vs {b}"
                );
            }
        }
        let stronger = KanaiTajimi {
            s0: 4.0 * model.s0,
            ..model
        };
        let actual = stronger.synthesize(29, 48, 0.01, 97);
        for (a, b) in actual.iter().zip(&reference) {
            assert_eq!(a.to_bits(), (2.0 * b).to_bits());
        }
    }

    #[test]
    fn g0_ground_motion_refuses_aliasing_and_invalid_stationary_models() {
        let model = spectral_model();
        let valid = model.default_sampling(8, 0.1, 8);
        let above_nyquist = GroundMotionSampling {
            omega_max: 1.001 * std::f64::consts::PI / valid.dt_s,
            ..valid
        };
        for invalid in [
            above_nyquist,
            GroundMotionSampling { dt_s: 0.0, ..valid },
            GroundMotionSampling { n_freq: 0, ..valid },
            GroundMotionSampling {
                omega_max: 0.0,
                ..valid
            },
        ] {
            assert!(
                std::panic::catch_unwind(|| model.synthesize_with_sampling(1, invalid)).is_err()
            );
        }
        for invalid in [
            KanaiTajimi { s0: -1.0, ..model },
            KanaiTajimi { wg: 0.0, ..model },
            KanaiTajimi { zg: 0.0, ..model },
        ] {
            assert!(
                std::panic::catch_unwind(|| invalid.synthesize_with_sampling(1, valid)).is_err()
            );
        }
        assert_eq!(
            KanaiTajimi { s0: 0.0, ..model }.synthesize_with_sampling(1, valid),
            vec![0.0; valid.n_steps]
        );
        assert!(
            model
                .synthesize_with_sampling(
                    1,
                    GroundMotionSampling {
                        n_steps: 0,
                        ..valid
                    }
                )
                .is_empty()
        );
    }

    #[test]
    fn g1_ida_uses_the_requested_bandwidth_record_length_and_sample_clock() {
        let model = spectral_model();
        let sampling = GroundMotionSampling {
            dt_s: 0.04,
            n_steps: 161,
            n_freq: 24,
            omega_max: 60.0,
        };
        let seeds = [2, 7, 13, 19, 29, 31, 43, 59];
        let records: Vec<_> = seeds
            .iter()
            .map(|&seed| model.synthesize_with_sampling(seed, sampling))
            .collect();
        let response = |record: &[f64], scale| {
            bilinear_peak_ductility(record, sampling.dt_s, 8.0, 0.05, 0.01, 0.1, scale)
        };
        let mut unit_peaks: Vec<_> = records.iter().map(|r| response(r, 1.0)).collect();
        unit_peaks.sort_by(f64::total_cmp);
        assert!(unit_peaks[4] > unit_peaks[3]);
        let limit = 0.5 * (unit_peaks[3] + unit_peaks[4]);
        let intensities = [0.25, 1.0, 4.0];
        let actual = ida_fragility_with_sampling(
            &model,
            sampling,
            &seeds,
            &intensities,
            8.0,
            0.05,
            0.01,
            0.1,
            limit,
        );
        for (point, &im) in actual.iter().zip(&intensities) {
            let hits = records.iter().filter(|r| response(r, im) > limit).count();
            assert_eq!(point.im, im);
            assert_eq!(point.n, seeds.len());
            assert_eq!(point.p, hits as f64 / seeds.len() as f64);
        }
        assert_eq!(
            actual[1].p, 0.5,
            "routing comparison must include mixed outcomes"
        );
        assert!(actual[0].p < actual[2].p);
    }
}
