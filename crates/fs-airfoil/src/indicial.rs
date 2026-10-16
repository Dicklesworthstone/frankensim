//! Indicial kernel definitions + exact-reference state machinery (bead
//! wf-root-guzez.5.1.3, E4.0c). Plan §5.2.1 + the Round-2 reduced-time
//! correction: the motion-memory clock is CHORDWISE reduced time,
//! `ds/dt = 2·U_conv/c` with `U_conv` the POSITIVE chordwise relative-flow
//! component — never the 3-D speed norm. `U_conv = 0` freezes the states;
//! reversed chordwise flow REFUSES rather than hiding behind `|U|`.
//!
//! Kernels are two-pole rational (exponential) approximations
//! `φ(s) = 1 − a₁·e^(−b₁ s) − a₂·e^(−b₂ s)` with the classical constant
//! sets (R.T. Jones' Wagner approximation; Sears/Küssner two-pole set)
//! as REGISTERED DEFAULTS. The per-step transition is the EXACT scalar
//! exponential for each pole (the 2-state matrix exponential is diagonal),
//! so step-size refinement changes nothing — an executable-exactness
//! battery pins that.

use crate::Refusal;
use fs_math::det;

/// Admitted reduced-time increment per update (a 120 Hz tick at Wright
/// speeds gives Δs ≈ 0.1–0.3; 64 is an absurd-input guard, not physics).
pub const MAX_DS: f64 = 64.0;

/// A two-pole indicial kernel: φ(s) = 1 − a₁e^(−b₁s) − a₂e^(−b₂s).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndicialKernel {
    /// Stable kernel id (registered; enters ModelId via the owner record).
    pub kernel_id: &'static str,
    /// Pole amplitudes (a₁, a₂), each ≥ 0, a₁ + a₂ ≤ 1.
    pub a: [f64; 2],
    /// Pole rates (b₁, b₂) in reduced time, each > 0.
    pub b: [f64; 2],
}

/// R.T. Jones' classical two-pole Wagner approximation:
/// φ(0) = 0.5, φ(∞) = 1.
pub const WAGNER_JONES: IndicialKernel = IndicialKernel {
    kernel_id: "wagner-jones-2pole-v1",
    a: [0.165, 0.335],
    b: [0.0455, 0.3],
};

/// Classical two-pole Küssner approximation: ψ(0) = 0, ψ(∞) = 1.
pub const KUSSNER_2POLE: IndicialKernel = IndicialKernel {
    kernel_id: "kussner-2pole-v1",
    a: [0.5, 0.5],
    b: [0.13, 1.0],
};

impl IndicialKernel {
    /// Validate the kernel parameters.
    ///
    /// # Errors
    /// `kernel-params-invalid` (non-finite, negative a, a₁+a₂ > 1, b ≤ 0).
    pub fn admit(&self) -> Result<(), Refusal> {
        let finite = self.a.iter().chain(self.b.iter()).all(|v| v.is_finite());
        let a_ok = self.a.iter().all(|&v| v >= 0.0) && self.a[0] + self.a[1] <= 1.0;
        let b_ok = self.b.iter().all(|&v| v > 0.0);
        if !(finite && a_ok && b_ok) {
            return Err(Refusal {
                code: "kernel-params-invalid",
                message: format!(
                    "kernel {}: a = {:?} (need aᵢ ≥ 0, Σa ≤ 1), b = {:?} (need bᵢ > 0)",
                    self.kernel_id, self.a, self.b
                ),
                ranked_repairs: vec![
                    "use a registered kernel constant (WAGNER_JONES, KUSSNER_2POLE)".into(),
                ],
            });
        }
        Ok(())
    }

    /// Closed-form φ(s) — the exact reference the state update must match.
    #[must_use]
    pub fn phi(&self, s: f64) -> f64 {
        1.0 - self.a[0] * det::exp(-self.b[0] * s) - self.a[1] * det::exp(-self.b[1] * s)
    }
}

/// Deficiency states for one kernel under a unit-step input:
/// x_i(s) = a_i·e^(−b_i·s), so the step response is φ(s) = 1 − x₁ − x₂.
/// States initialize to the TRIM steady state (x = 0, fully developed)
/// unless an impulsive start is explicitly requested (plan §5.1.5).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndicialState {
    /// Deficiency components (x₁, x₂).
    pub x: [f64; 2],
}

impl IndicialState {
    /// Trim steady state: memory fully developed, zero deficiency.
    #[must_use]
    pub fn trim() -> Self {
        IndicialState { x: [0.0, 0.0] }
    }

    /// Impulsive (Wagner) start for a unit step: x_i = a_i, φ = φ(0).
    #[must_use]
    pub fn impulsive_start(kernel: &IndicialKernel) -> Self {
        IndicialState {
            x: [kernel.a[0], kernel.a[1]],
        }
    }

    /// Current step-response value φ = 1 − x₁ − x₂.
    #[must_use]
    pub fn response(&self) -> f64 {
        1.0 - self.x[0] - self.x[1]
    }

    /// Advance the states by the EXACT diagonal matrix exponential over a
    /// reduced-time increment `ds`: x_i ← x_i·e^(−b_i·ds). Exactness means
    /// N sub-steps compose to the identical result as one step (pinned by
    /// the battery to 1e-14).
    ///
    /// # Errors
    /// `reduced-time-increment-invalid` (non-finite, negative, or above
    /// [`MAX_DS`] — tested at cap and cap+1).
    pub fn advance(&mut self, kernel: &IndicialKernel, ds: f64) -> Result<(), Refusal> {
        if !ds.is_finite() || ds < 0.0 || ds > MAX_DS {
            return Err(Refusal {
                code: "reduced-time-increment-invalid",
                message: format!("ds = {ds:?} outside admitted [0, {MAX_DS}]"),
                ranked_repairs: vec![
                    "reduced time never runs backwards; check the U_conv clock".into(),
                ],
            });
        }
        self.x[0] *= det::exp(-kernel.b[0] * ds);
        self.x[1] *= det::exp(-kernel.b[1] * ds);
        Ok(())
    }
}

/// The chordwise reduced-time increment (Round-2 correction, plan §5.2):
/// `Δs = 2·U_conv·Δt / c` where `U_conv` is the POSITIVE chordwise
/// relative-flow component in the section frame.
///
/// - `U_conv = 0` → `Δs = 0` (states FREEZE; a pure vertical gust must not
///   advance the motion-memory clock).
/// - `U_conv < 0` (reversed chordwise flow) → typed REFUSAL: the indicial
///   owner is out of its domain and the caller must switch owners, never
///   hide behind `|U|`.
///
/// # Errors
/// `non-finite-input`, `chord-outside-domain` (chord ≤ 0),
/// `timestep-invalid` (dt < 0), `indicial-flow-reversed`.
pub fn reduced_time_increment(u_conv_mps: f64, chord_m: f64, dt_s: f64) -> Result<f64, Refusal> {
    if !(u_conv_mps.is_finite() && chord_m.is_finite() && dt_s.is_finite()) {
        return Err(Refusal {
            code: "non-finite-input",
            message: format!("u_conv {u_conv_mps:?}, chord {chord_m:?}, dt {dt_s:?}"),
            ranked_repairs: vec!["check the section-frame projection upstream".into()],
        });
    }
    if chord_m <= 0.0 {
        return Err(Refusal {
            code: "chord-outside-domain",
            message: format!("chord {chord_m} m must be positive"),
            ranked_repairs: vec!["pass the section chord in metres".into()],
        });
    }
    if dt_s < 0.0 {
        return Err(Refusal {
            code: "timestep-invalid",
            message: format!("dt {dt_s} s must be non-negative"),
            ranked_repairs: vec!["the simulation clock never runs backwards".into()],
        });
    }
    if u_conv_mps < 0.0 {
        return Err(Refusal {
            code: "indicial-flow-reversed",
            message: format!(
                "chordwise relative flow U_conv = {u_conv_mps} m/s is reversed — the 2-D \
                 indicial owner is outside its domain"
            ),
            ranked_repairs: vec![
                "switch the motion_circulatory owner (AeroEffectOwners) for this strip".into(),
                "never absolute-value the clock: a reversed-flow state has different physics"
                    .into(),
            ],
        });
    }
    Ok(2.0 * u_conv_mps * dt_s / chord_m)
}
