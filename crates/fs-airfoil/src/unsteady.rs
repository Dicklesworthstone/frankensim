//! Two-channel unsteady section response + lagged separation (bead
//! wf-root-guzez.5.6.1, E4.3-i). Plan §5.2.1/V-08a: the LIFT and MOMENT
//! channels are assembled from the registered indicial kernels on the
//! chordwise reduced-time clock, with a first-order lagged separation
//! coordinate modulating both through the Kirchhoff factor. The
//! NONCIRCULATORY (added-mass) contribution is deliberately ABSENT here —
//! its owner is AddedMassOnly (fs-flyer), and the AeroEffectOwners
//! admission (E4.3-ii) is what makes that single-ownership law checkable.
//!
//! Numerics: `advance` retains the midpoint-step convention.
//! `advance_linear` exactly integrates each circulatory deficiency pole
//! for an input linear in REDUCED time between the accepted endpoints;
//! a zero reduced-time increment freezes its complete state. `apply_jump`
//! explicitly applies an instantaneous input change without relaxation.
//! Separation retains the scalar exponential toward the endpoint effective
//! incidence's static target, not an exact integration of a varying target.
//! All section updates are transactional; composition is subject to roundoff.

use crate::Refusal;
use crate::indicial::{IndicialKernel, MAX_DS};
use fs_math::det;

fn admit_reduced_time(ds: f64) -> Result<(), Refusal> {
    if !ds.is_finite() || ds < 0.0 || ds > MAX_DS {
        return Err(Refusal {
            code: "reduced-time-increment-invalid",
            message: format!("ds = {ds:?} outside admitted [0, {MAX_DS}]"),
            ranked_repairs: vec![
                "reduced time never runs backwards; check the U_conv clock".into(),
            ],
        });
    }
    Ok(())
}

fn admit_input(input: f64) -> Result<(), Refusal> {
    if !input.is_finite() {
        return Err(Refusal {
            code: "non-finite-input",
            message: format!("section input {input:?}"),
            ranked_repairs: vec!["check the quasi-steady alpha assembly upstream".into()],
        });
    }
    Ok(())
}

fn invalid_state(message: &str) -> Refusal {
    Refusal {
        code: "unsteady-state-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "use finite admitted section states and representable inputs and coefficients".into(),
        ],
    }
}

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
    /// cap+1); `separation-model-invalid`, `non-finite-input`, or
    /// `unsteady-state-invalid` (f outside [0, 1]). Refusal leaves f unchanged.
    pub fn advance(
        &mut self,
        model: &SeparationLagModel,
        ds: f64,
        alpha_rad: f64,
    ) -> Result<(), Refusal> {
        admit_reduced_time(ds)?;
        model.admit()?;
        admit_input(alpha_rad)?;
        if !self.f.is_finite() || !(0.0..=1.0).contains(&self.f) {
            return Err(invalid_state("separation coordinate must lie in [0, 1]"));
        }
        // Avoid cancellation in f_st + (f - f_st) at zero elapsed time.
        if ds == 0.0 {
            return Ok(());
        }
        let f_st = model.f_static(alpha_rad);
        self.f = f_st + (self.f - f_st) * det::exp(-ds / model.t_f);
        Ok(())
    }
}

/// Duhamel deficiency states: alpha_eff = alpha - x_1 - x_2.
/// Each pole obeys dx_i/ds = -b_i*x_i + a_i*d(alpha)/ds.
/// Select midpoint steps, endpoint-linear ramps, or explicit jumps.
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

    /// Advance with an input step at the MIDPOINT of `ds`. This retains
    /// the original arithmetic and convention: at ds = 0 a changed input
    /// is an instantaneous jump, not a frozen ramp endpoint.
    ///
    /// # Errors
    /// `reduced-time-increment-invalid`, `kernel-params-invalid`,
    /// `non-finite-input`, or `unsteady-state-invalid`. All refusals
    /// leave the complete state unchanged.
    pub fn advance(&mut self, kernel: &IndicialKernel, ds: f64, input: f64) -> Result<(), Refusal> {
        self.admit_update(kernel, ds, input)?;
        // MIDPOINT delta-s (plan §5.2.1): the input change is booked at
        // mid-interval, so the new deficiency decays over ds/2 — exact
        // for a step at the midpoint, second-order for smooth input.
        let d = input - self.input_prev;
        let x = [
            self.x[0] * det::exp(-kernel.b[0] * ds)
                + kernel.a[0] * d * det::exp(-kernel.b[0] * ds / 2.0),
            self.x[1] * det::exp(-kernel.b[1] * ds)
                + kernel.a[1] * d * det::exp(-kernel.b[1] * ds / 2.0),
        ];
        self.accept(x, input)
    }

    /// Exactly integrate an input linear in REDUCED time from
    /// `input_prev` to `input` over `ds`. For z_i = b_i*ds:
    /// x_i' = exp(-z_i)*x_i + a_i*(input-input_prev)*(1-exp(-z_i))/z_i.
    /// The ramp factor uses expm1 and its z = 0 limit, avoiding cancellation.
    ///
    /// At ds = 0 the complete state, including input_prev, is unchanged;
    /// use [`Self::apply_jump`] to declare an instantaneous change instead.
    /// Linearity in physical time is equivalent only for constant convection
    /// speed. This update allocates no memory on the admitted path.
    ///
    /// # Errors
    /// The same refusals as [`Self::advance`], without partial mutation.
    pub fn advance_linear(
        &mut self,
        kernel: &IndicialKernel,
        ds: f64,
        input: f64,
    ) -> Result<(), Refusal> {
        self.admit_update(kernel, ds, input)?;
        if ds == 0.0 {
            return Ok(());
        }
        let d = input - self.input_prev;
        let mut x = [0.0; 2];
        for i in 0..2 {
            let z = kernel.b[i] * ds;
            let ramp = if z == 0.0 {
                1.0
            } else if z.is_infinite() {
                // Finite admitted b and ds can overflow their product while
                // 1/(b*ds) is still a representable subnormal.
                (1.0 / kernel.b[i]) / ds
            } else {
                -det::expm1(-z) / z
            };
            let carried = self.x[i] * det::exp(-z);
            let driven = kernel.a[i] * d * ramp;
            x[i] = carried + driven;
        }
        self.accept(x, input)
    }

    /// Apply a discontinuous input change with no elapsed reduced time.
    /// x_i' = x_i + a_i*delta, so the instantaneous effective-input change
    /// is (1-sum(a))*delta: half for Wagner, zero for Kussner.
    ///
    /// # Errors
    /// `kernel-params-invalid`, `non-finite-input`, or
    /// `unsteady-state-invalid`; no state changes on refusal.
    pub fn apply_jump(&mut self, kernel: &IndicialKernel, input: f64) -> Result<(), Refusal> {
        self.admit_update(kernel, 0.0, input)?;
        let d = input - self.input_prev;
        let x = [self.x[0] + kernel.a[0] * d, self.x[1] + kernel.a[1] * d];
        self.accept(x, input)
    }

    fn admit_update(&self, kernel: &IndicialKernel, ds: f64, input: f64) -> Result<(), Refusal> {
        admit_reduced_time(ds)?;
        kernel.admit()?;
        admit_input(input)?;
        if !self.input_prev.is_finite()
            || !self.x.iter().all(|x| x.is_finite())
            || !self.effective().is_finite()
        {
            return Err(invalid_state("non-finite Duhamel state"));
        }
        Ok(())
    }

    fn accept(&mut self, x: [f64; 2], input: f64) -> Result<(), Refusal> {
        let next = Self {
            x,
            input_prev: input,
        };
        if !x.iter().all(|value| value.is_finite()) || !next.effective().is_finite() {
            return Err(invalid_state(
                "Duhamel update produced a non-finite response",
            ));
        }
        *self = next;
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

impl UnsteadySectionSpec {
    /// Admit both kernels, separation parameters, and finite load constants.
    ///
    /// # Errors
    /// `kernel-params-invalid`, `separation-model-invalid`, or
    /// `unsteady-section-spec-invalid` for non-finite load constants.
    pub fn admit(&self) -> Result<(), Refusal> {
        self.motion_kernel.admit()?;
        self.gust_kernel.admit()?;
        self.separation.admit()?;
        if ![self.cl_alpha, self.cm0, self.ac_shift]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(Refusal {
                code: "unsteady-section-spec-invalid",
                message: "lift slope, baseline moment and center shift must be finite".into(),
                ranked_repairs: vec!["use finite admitted section load constants".into()],
            });
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum InputEvolution {
    MidpointStep,
    LinearRamp,
    Jump,
}

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

    /// Advance both inputs with the original midpoint-step convention,
    /// then relax separation toward the endpoint EFFECTIVE incidence.
    /// All states commit together only after finite lift and moment result.
    ///
    /// # Errors
    /// Propagates spec/state refusals; every refusal preserves all memories.
    pub fn advance(
        &mut self,
        spec: &UnsteadySectionSpec,
        ds: f64,
        alpha_qs_rad: f64,
        alpha_gust_rad: f64,
    ) -> Result<SectionChannels, Refusal> {
        self.advance_with(
            spec,
            ds,
            alpha_qs_rad,
            alpha_gust_rad,
            InputEvolution::MidpointStep,
        )
    }

    /// Integrate both circulatory inputs exactly as endpoint-linear ramps in
    /// REDUCED time. Separation retains its endpoint-held target exponential;
    /// a varying separation target is not integrated exactly by this method.
    /// At ds = 0 every state bit freezes, including both previous inputs.
    ///
    /// # Errors
    /// Propagates spec/state refusals; every refusal preserves all memories.
    pub fn advance_linear(
        &mut self,
        spec: &UnsteadySectionSpec,
        ds: f64,
        alpha_qs_rad: f64,
        alpha_gust_rad: f64,
    ) -> Result<SectionChannels, Refusal> {
        self.advance_with(
            spec,
            ds,
            alpha_qs_rad,
            alpha_gust_rad,
            InputEvolution::LinearRamp,
        )
    }

    /// Apply explicit instantaneous motion/gust jumps. Circulation has each
    /// kernel's instantaneous gain; the separation coordinate does not relax.
    ///
    /// # Errors
    /// Propagates spec/state refusals; every refusal preserves all memories.
    pub fn apply_jump(
        &mut self,
        spec: &UnsteadySectionSpec,
        alpha_qs_rad: f64,
        alpha_gust_rad: f64,
    ) -> Result<SectionChannels, Refusal> {
        self.advance_with(
            spec,
            0.0,
            alpha_qs_rad,
            alpha_gust_rad,
            InputEvolution::Jump,
        )
    }

    fn advance_with(
        &mut self,
        spec: &UnsteadySectionSpec,
        ds: f64,
        alpha_qs_rad: f64,
        alpha_gust_rad: f64,
        evolution: InputEvolution,
    ) -> Result<SectionChannels, Refusal> {
        spec.admit()?;
        let mut next = *self;
        for (state, kernel, input) in [
            (&mut next.motion, &spec.motion_kernel, alpha_qs_rad),
            (&mut next.gust, &spec.gust_kernel, alpha_gust_rad),
        ] {
            match evolution {
                InputEvolution::MidpointStep => state.advance(kernel, ds, input)?,
                InputEvolution::LinearRamp => state.advance_linear(kernel, ds, input)?,
                InputEvolution::Jump => state.apply_jump(kernel, input)?,
            }
        }
        let alpha_eff = next.motion.effective() + next.gust.effective();
        next.separation.advance(&spec.separation, ds, alpha_eff)?;
        let f = next.separation.f;
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
        if !cl.is_finite() || !cm_quarter.is_finite() {
            return Err(invalid_state(
                "section update produced non-finite lift or moment",
            ));
        }
        *self = next;
        Ok(SectionChannels {
            cl,
            cm_quarter,
            alpha_eff_rad: alpha_eff,
            f,
        })
    }
}
