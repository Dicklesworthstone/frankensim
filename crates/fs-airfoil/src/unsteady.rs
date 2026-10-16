//! Two-channel unsteady section response + lagged separation (bead
//! wf-root-guzez.5.6.1, E4.3-i). Plan §5.2.1/V-08a: the LIFT and MOMENT
//! channels are assembled from the registered indicial kernels on the
//! chordwise reduced-time clock, with a first-order lagged separation
//! coordinate modulating both through the Kirchhoff factor. The
//! NONCIRCULATORY (added-mass) contribution is deliberately ABSENT here —
//! its owner is AddedMassOnly (fs-flyer), and the AeroEffectOwners
//! admission (E4.3-ii) is what makes that single-ownership law checkable.
//!
//! Numerics: per step the input is piecewise-constant, so the Duhamel
//! deficiency recurrence x_i ← x_i·e^(−b_i·Δs) + a_i·Δα and the
//! separation lag f ← f_st + (f − f_st)·e^(−Δs/T_f) are EXACT scalar
//! exponentials — substeps with the same input compose bitwise-identically
//! to one step, and the battery pins that.

use crate::Refusal;
use crate::indicial::{IndicialKernel, MAX_DS};
use fs_math::det;

/// Registered separation-lag model (Kirchhoff-class; ceiling `Estimated`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeparationLagModel {
    /// Stable model id (enters ModelId via the owner record).
    pub model_id: &'static str,
    /// Separation lag constant in reduced time (semichords).
    pub t_f: f64,
    /// |α| below which the static separation coordinate is 1 (attached).
    pub alpha_break_rad: f64,
    /// Exponential decay scale of f_static past the break [rad].
    pub decay_rad: f64,
}

/// The registered default (classical T_f ≈ 3 semichords; break/decay from
/// the E4.0 flat-plate blend window 0.30→0.45 rad).
pub const SEPARATION_LAG_V1: SeparationLagModel = SeparationLagModel {
    model_id: "separation-lag-kirchhoff-v1",
    t_f: 3.0,
    alpha_break_rad: 0.30,
    decay_rad: 0.10,
};

impl SeparationLagModel {
    /// Validate the model parameters.
    ///
    /// # Errors
    /// `separation-model-invalid` (non-finite or non-positive constants).
    pub fn admit(&self) -> Result<(), Refusal> {
        let ok = self.t_f.is_finite()
            && self.t_f > 0.0
            && self.alpha_break_rad.is_finite()
            && self.alpha_break_rad > 0.0
            && self.decay_rad.is_finite()
            && self.decay_rad > 0.0;
        if !ok {
            return Err(Refusal {
                code: "separation-model-invalid",
                message: format!(
                    "{}: t_f {:?}, break {:?}, decay {:?} (all must be finite and positive)",
                    self.model_id, self.t_f, self.alpha_break_rad, self.decay_rad
                ),
                ranked_repairs: vec!["use SEPARATION_LAG_V1".into()],
            });
        }
        Ok(())
    }

    /// Static separation coordinate f_st(α) ∈ (0, 1]: 1 while attached,
    /// exponential decay past the break. Even in α (separation does not
    /// care about sign).
    #[must_use]
    pub fn f_static(&self, alpha_rad: f64) -> f64 {
        let a = alpha_rad.abs();
        if a <= self.alpha_break_rad {
            1.0
        } else {
            det::exp(-(a - self.alpha_break_rad) / self.decay_rad)
        }
    }
}

/// Lagged separation coordinate f (the dynamic-stall delay state).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeparationLagState {
    /// Current separation coordinate f ∈ [0, 1].
    pub f: f64,
}

impl SeparationLagState {
    /// Trim initialization: f at its static value for the trim α.
    #[must_use]
    pub fn trim(model: &SeparationLagModel, alpha_rad: f64) -> Self {
        SeparationLagState {
            f: model.f_static(alpha_rad),
        }
    }

    /// Advance by the EXACT scalar exponential over Δs toward f_st(α):
    /// f ← f_st + (f − f_st)·e^(−Δs/T_f). `ds = 0` freezes (`U_conv = 0`
    /// must not relax separation memory).
    ///
    /// # Errors
    /// `reduced-time-increment-invalid` (non-finite, negative, > `MAX_DS`
    /// — same admitted window as the kernel clock, tested at cap and
    /// cap+1).
    pub fn advance(
        &mut self,
        model: &SeparationLagModel,
        ds: f64,
        alpha_rad: f64,
    ) -> Result<(), Refusal> {
        if !ds.is_finite() || ds < 0.0 || ds > MAX_DS {
            return Err(Refusal {
                code: "reduced-time-increment-invalid",
                message: format!("ds = {ds:?} outside admitted [0, {MAX_DS}]"),
                ranked_repairs: vec![
                    "reduced time never runs backwards; check the U_conv clock".into(),
                ],
            });
        }
        let f_st = model.f_static(alpha_rad);
        self.f = f_st + (self.f - f_st) * det::exp(-ds / model.t_f);
        Ok(())
    }
}

/// Duhamel deficiency states for one kernel under ARBITRARY
/// piecewise-constant input: α_eff = α − x₁ − x₂ with
/// x_i ← x_i·e^(−b_i·Δs) + a_i·Δα per step (exact for the held input).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuhamelState {
    /// Deficiency components.
    pub x: [f64; 2],
    /// The input value the deficiencies are measured against.
    pub input_prev: f64,
}

impl DuhamelState {
    /// Trim initialization: memory fully developed at the trim input
    /// (zero deficiency — no startup transient, plan §5.1.5).
    #[must_use]
    pub fn trim(input: f64) -> Self {
        DuhamelState {
            x: [0.0, 0.0],
            input_prev: input,
        }
    }

    /// Effective (lag-filtered) input.
    #[must_use]
    pub fn effective(&self) -> f64 {
        self.input_prev - self.x[0] - self.x[1]
    }

    /// Advance over `ds` with the new held input value.
    ///
    /// # Errors
    /// `reduced-time-increment-invalid` (same admitted window as the
    /// kernel clock); `non-finite-input`.
    pub fn advance(&mut self, kernel: &IndicialKernel, ds: f64, input: f64) -> Result<(), Refusal> {
        if !ds.is_finite() || ds < 0.0 || ds > MAX_DS {
            return Err(Refusal {
                code: "reduced-time-increment-invalid",
                message: format!("ds = {ds:?} outside admitted [0, {MAX_DS}]"),
                ranked_repairs: vec![
                    "reduced time never runs backwards; check the U_conv clock".into(),
                ],
            });
        }
        if !input.is_finite() {
            return Err(Refusal {
                code: "non-finite-input",
                message: format!("held input {input:?}"),
                ranked_repairs: vec!["check the quasi-steady α assembly upstream".into()],
            });
        }
        // MIDPOINT delta-s (plan §5.2.1): the input change is booked at
        // mid-interval, so the new deficiency decays over ds/2 — exact
        // for a step at the midpoint, second-order for smooth input.
        let d = input - self.input_prev;
        self.x[0] = self.x[0] * det::exp(-kernel.b[0] * ds)
            + kernel.a[0] * d * det::exp(-kernel.b[0] * ds / 2.0);
        self.x[1] = self.x[1] * det::exp(-kernel.b[1] * ds)
            + kernel.a[1] * d * det::exp(-kernel.b[1] * ds / 2.0);
        self.input_prev = input;
        Ok(())
    }
}

/// Two-channel unsteady section spec: kernels + separation model +
/// linear-range constants. NO noncirculatory term (AddedMassOnly owner).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnsteadySectionSpec {
    /// Motion-circulatory kernel (Wagner class).
    pub motion_kernel: IndicialKernel,
    /// Incident-gust kernel (Küssner class).
    pub gust_kernel: IndicialKernel,
    /// Separation lag model.
    pub separation: SeparationLagModel,
    /// Attached lift slope [1/rad].
    pub cl_alpha: f64,
    /// Zero-lift moment about quarter chord.
    pub cm0: f64,
    /// Aft aerodynamic-center shift at full separation (fraction of
    /// chord): x_ac(f) = 0.25 + shift·(1 − f).
    pub ac_shift: f64,
}

/// The registered default section (thin-airfoil slope; modest aft shift).
pub const UNSTEADY_SECTION_V1: UnsteadySectionSpec = UnsteadySectionSpec {
    motion_kernel: crate::indicial::WAGNER_JONES,
    gust_kernel: crate::indicial::KUSSNER_2POLE,
    separation: SEPARATION_LAG_V1,
    cl_alpha: core::f64::consts::TAU,
    cm0: 0.0,
    ac_shift: 0.075,
};

/// Full unsteady section state (motion + gust + separation).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnsteadySectionState {
    /// Motion-circulatory Duhamel states on quasi-steady α.
    pub motion: DuhamelState,
    /// Gust Duhamel states on the gust upwash angle.
    pub gust: DuhamelState,
    /// Lagged separation coordinate.
    pub separation: SeparationLagState,
}

/// One tick's two-channel output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionChannels {
    /// Circulatory lift coefficient (Kirchhoff-modulated).
    pub cl: f64,
    /// Quarter-chord moment coefficient.
    pub cm_quarter: f64,
    /// The effective circulatory α the channels were built from [rad].
    pub alpha_eff_rad: f64,
    /// Current separation coordinate.
    pub f: f64,
}

impl UnsteadySectionState {
    /// Trim initialization at (α_trim, gust 0): all memories fully
    /// developed — the first tick's output equals the static answer.
    #[must_use]
    pub fn trim(spec: &UnsteadySectionSpec, alpha_trim_rad: f64) -> Self {
        UnsteadySectionState {
            motion: DuhamelState::trim(alpha_trim_rad),
            gust: DuhamelState::trim(0.0),
            separation: SeparationLagState::trim(&spec.separation, alpha_trim_rad),
        }
    }

    /// Advance all states over `ds` with held inputs, then evaluate both
    /// channels. The separation lag is driven by the EFFECTIVE α (the lag
    /// cascade of dynamic stall: circulation lag first, then f lag).
    ///
    /// # Errors
    /// Propagates the state refusals (`reduced-time-increment-invalid`,
    /// `non-finite-input`).
    pub fn advance(
        &mut self,
        spec: &UnsteadySectionSpec,
        ds: f64,
        alpha_qs_rad: f64,
        alpha_gust_rad: f64,
    ) -> Result<SectionChannels, Refusal> {
        self.motion.advance(&spec.motion_kernel, ds, alpha_qs_rad)?;
        self.gust.advance(&spec.gust_kernel, ds, alpha_gust_rad)?;
        let alpha_eff = self.motion.effective() + self.gust.effective();
        self.separation.advance(&spec.separation, ds, alpha_eff)?;
        let f = self.separation.f;
        // Kirchhoff: CL = CLα·α_eff·((1+√f)/2)².
        let kirchhoff = {
            let k = (1.0 + det::sqrt(f.max(0.0))) / 2.0;
            k * k
        };
        let cl = spec.cl_alpha * alpha_eff * kirchhoff;
        // Moment: cm0 + CL·(0.25 − x_ac(f)) with x_ac aft-shifting as
        // separation grows (both channels carry the kernel dynamics).
        let x_ac = 0.25 + spec.ac_shift * (1.0 - f);
        let cm_quarter = spec.cm0 + cl * (0.25 - x_ac);
        Ok(SectionChannels {
            cl,
            cm_quarter,
            alpha_eff_rad: alpha_eff,
            f,
        })
    }
}
