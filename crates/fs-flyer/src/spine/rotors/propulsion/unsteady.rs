//! Native carrier/shaft dynamics with motion, gust and separation memory.
//!
//! This explicit model choice replaces the wing section lift/moment closure
//! inside the existing two-way wing/BEMT solve. Each nonlinear trial samples a
//! COPY of the same physical history; iterations never advance its clock.
//! The two real mechanical stages advance history in chordwise semichords.
//!
//! The registered Kirchhoff model is an Estimated attached/stall-onset model,
//! admitted for both instantaneous and memory-effective geometric incidence
//! magnitudes <= 0.45 rad at every force sample. It is not a
//! deep-stall, dynamic-stall-vortex or resolved three-dimensional wake model.
//! Profile/induced/parasite drag and axial equilibrium propeller wakes retain
//! their existing owners. No added-mass term or historical calibration is added.

use super::{
    AerodynamicFields, AerodynamicSolution, AerodynamicStage, CoupledStep, GearedLoads,
    GearedPropulsion, GearedState, GearedStep, PropulsionInputs, Refusal, geared, reject,
    wing_error,
};
use fs_airfoil::unsteady::{
    DuhamelState, SectionChannels, SectionMidpointInputs, SeparationLagState, UnsteadySectionSpec,
    UnsteadySectionState,
};
use fs_wing::nonlinear::{SectionFrame, StripRegime};
use std::cell::RefCell;

#[path = "unsteady/gust.rs"]
pub(super) mod gust;
use gust::GustInput;

/// Conservative validity boundary of this attached/stall-onset runtime tier.
pub const MAX_EFFECTIVE_ALPHA_RAD: f64 = 0.45;

/// Complete accepted mechanical and section state, borrowed immutably by steps.
#[derive(Clone, Debug, PartialEq)]
pub struct UnsteadyState {
    pub motion: GearedState,
    /// One motion/gust/separation history per physical lifting strip.
    pub sections: Vec<UnsteadySectionState>,
}

/// Coupled wing/propeller solution and the section channels that generated it.
#[derive(Clone, Debug, PartialEq)]
pub struct UnsteadyAerodynamics {
    pub coupled: CoupledStep,
    /// Lift and pitching moment evaluated from identical frozen stage memory.
    pub sections: Vec<SectionChannels>,
    /// Declared gust component of total geometric incidence [rad], either
    /// explicit or derived from the same-point total/reference incident air.
    pub gust_angle_rad: Vec<f64>,
    /// Actual washed forward chordwise speed used by the physical clock [m/s].
    pub chordwise_speed_mps: Vec<f64>,
}

impl AerodynamicSolution for UnsteadyAerodynamics {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        self.coupled.load_fields()
    }

    fn section_channels(&self) -> Option<&[SectionChannels]> {
        Some(&self.sections)
    }
}

/// One real native force evaluation; no physical history is advanced.
#[derive(Clone, Debug, PartialEq)]
pub struct UnsteadyEvaluation<const N: usize> {
    pub loads: GearedLoads<N>,
    pub aerodynamic_force_body_n: [f64; 3],
    pub propeller_absolute_speed_rad_s: [f64; 2],
    pub aerodynamics: UnsteadyAerodynamics,
}

/// Atomically accepted motion/history and both real force evaluations.
#[derive(Clone, Debug, PartialEq)]
pub struct UnsteadyStep<const N: usize> {
    pub state: UnsteadyState,
    pub motion: GearedStep<N>,
    pub initial: UnsteadyEvaluation<N>,
    pub midpoint: UnsteadyEvaluation<N>,
}

#[derive(Clone, Copy)]
struct SectionModel {
    spec: UnsteadySectionSpec,
    zero_lift_angle_rad: f64,
    chord_m: f64,
    frame: SectionFrame,
}

/// Opt-in unsteady wing closure over the original native geared force owner.
pub struct UnsteadyPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p GearedPropulsion<'a, N>,
    sections: Vec<SectionModel>,
}

fn airfoil_error(error: fs_airfoil::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

fn invalid(message: &str) -> Refusal {
    reject("unsteady-propulsion-invalid", message)
}

fn admit_incidence(alpha: f64) -> Result<(), Refusal> {
    if !alpha.is_finite() || alpha.abs() > MAX_EFFECTIVE_ALPHA_RAD {
        return Err(invalid(
            "effective incidence exceeds the attached/stall-onset domain",
        ));
    }
    Ok(())
}

impl SectionModel {
    fn developed(self, alpha: f64, gust: f64) -> Result<UnsteadySectionState, Refusal> {
        admit_incidence(alpha)?;
        let state = UnsteadySectionState {
            motion: DuhamelState::trim(alpha - gust),
            gust: DuhamelState::trim(gust),
            separation: SeparationLagState::trim(&self.spec.separation, alpha),
        };
        let response = state
            .response(&self.spec, self.zero_lift_angle_rad)
            .map_err(airfoil_error)?;
        admit_incidence(response.alpha_eff_rad)?;
        Ok(state)
    }

    fn sample(
        self,
        history: UnsteadySectionState,
        alpha: f64,
        gust: f64,
    ) -> Result<SectionChannels, Refusal> {
        admit_incidence(alpha)?;
        let mut trial = history;
        trial
            .apply_jump(&self.spec, alpha - gust, gust)
            .map_err(airfoil_error)?;
        let response = trial
            .response(&self.spec, self.zero_lift_angle_rad)
            .map_err(airfoil_error)?;
        admit_incidence(response.alpha_eff_rad)?;
        Ok(response)
    }

    fn regime(self, alpha: f64) -> StripRegime {
        if alpha.abs() <= self.spec.separation.alpha_break_rad {
            StripRegime::Attached
        } else {
            StripRegime::Blended
        }
    }
}

impl<'p, 'a, const N: usize> UnsteadyPropulsion<'p, 'a, N> {
    /// Prepare immutable chord/frame/camber metadata once. The selected spec's
    /// baseline moment is augmented by each strip's native camber moment;
    /// its zero-lift angle is -2*camber. Camber never drives separation memory.
    /// Use [`Self::trimmed_state`] to initialize a developed coupled flight.
    ///
    /// # Errors
    /// Native geometry/spec admission, or a nonpositive selected lift slope.
    pub fn new(
        propulsion: &'p GearedPropulsion<'a, N>,
        spec: UnsteadySectionSpec,
    ) -> Result<Self, Refusal> {
        spec.admit().map_err(airfoil_error)?;
        if spec.cl_alpha <= 0.0 {
            return Err(invalid("the attached lift slope must be positive"));
        }
        let strips = match propulsion.surfaces {
            Some(model) => model.controlled_strips(0.0)?.into_owned(),
            None => propulsion.design.layout(0.0)?.1,
        };
        let sections = strips
            .iter()
            .enumerate()
            .map(|(s, strip)| {
                let (camber, frame) = propulsion.surfaces.map_or(
                    (propulsion.design.camber_ratio, SectionFrame::BODY_XZ),
                    |model| (model.properties[s].camber_ratio, model.frames[s]),
                );
                let section_spec = UnsteadySectionSpec {
                    cm0: spec.cm0 - core::f64::consts::PI * camber,
                    ..spec
                };
                section_spec.admit().map_err(airfoil_error)?;
                Ok(SectionModel {
                    spec: section_spec,
                    zero_lift_angle_rad: -2.0 * camber,
                    chord_m: strip.chord_m,
                    frame,
                })
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        Ok(Self {
            propulsion,
            sections,
        })
    }

    #[must_use]
    pub fn section_count(&self) -> usize {
        self.sections.len()
    }

    fn admit_inputs(&self, gust: GustInput<'_>, warm: Option<[f64; 2]>) -> Result<(), Refusal> {
        if gust.invalid_angles(self.sections.len())
            || warm.is_some_and(|w| !w.iter().all(|value| value.is_finite()))
        {
            return Err(invalid(
                "one finite explicit gust angle per strip and finite warm induction are required",
            ));
        }
        Ok(())
    }

    fn admit_history(&self, history: &[UnsteadySectionState]) -> Result<(), Refusal> {
        if history.len() != self.sections.len() {
            return Err(invalid(
                "one complete physical history per strip is required",
            ));
        }
        for (state, model) in history.iter().zip(&self.sections) {
            // A stored ramp endpoint is only a reference. Physical angle
            // bounds apply after sample rebases it to the actual stage input.
            state
                .response(&model.spec, model.zero_lift_angle_rad)
                .map_err(airfoil_error)?;
        }
        Ok(())
    }

    /// Initialize developed memory from a real steady solve of the SELECTED
    /// section polar, including twist, mutual wing induction and prop wash.
    /// Nonzero gust is already developed in its own channel; motion receives
    /// total effective geometric incidence minus that same gust component.
    /// This initialization performs one coupled force solve.
    ///
    /// `wind_world` is the TOTAL material-air velocity. `gust_angle_rad` is an
    /// explicit per-strip decomposition supplied by that field's owner; it is
    /// not additional wind. Zero means all incidence uses the motion channel.
    ///
    /// # Errors
    /// Native coupled-force admission or invalid/out-of-domain section data.
    pub fn trimmed_state(
        &self,
        motion: GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust_angle_rad: &[f64],
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyState, Refusal> {
        self.trimmed_state_with_gust(
            motion,
            inputs,
            wind_world,
            GustInput::Angles(gust_angle_rad),
            warm,
        )
    }

    /// Initialize developed memory using the actual external gust incidence:
    /// total minus reference normal-section angle at each native strip point.
    /// Both fields use the identical pose and rigid-point velocity. Only total
    /// wind enters loads; the declared reference partitions the two memories.
    /// Invalid/reversed reference chordwise flow refuses without a polar floor.
    pub fn trimmed_state_with_reference_wind(
        &self,
        motion: GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        reference_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyState, Refusal> {
        self.trimmed_state_with_gust(
            motion,
            inputs,
            wind_world,
            GustInput::Reference(reference_world),
            warm,
        )
    }

    fn trimmed_state_with_gust(
        &self,
        motion: GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust: GustInput<'_>,
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyState, Refusal> {
        self.admit_inputs(gust, warm)?;
        let result = self
            .propulsion
            .evaluate_with(&motion, inputs, wind_world, |stage| {
                let gust = gust.resolve(&stage, self.sections.len(), |s| self.sections[s].frame)?;
                self.solve(stage, None, &gust, warm)
            })?;
        let history = result
            .aerodynamics
            .coupled
            .strip_alpha_eff_rad
            .iter()
            .enumerate()
            .map(|(s, &alpha)| {
                self.sections[s].developed(alpha, result.aerodynamics.gust_angle_rad[s])
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(UnsteadyState {
            motion,
            sections: history,
        })
    }

    fn solve(
        &self,
        stage: AerodynamicStage<'_>,
        history: Option<&[UnsteadySectionState]>,
        gust: &[f64],
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyAerodynamics, Refusal> {
        if stage.strips.len() != self.sections.len() {
            return Err(invalid(
                "the bound geometry no longer matches its section history",
            ));
        }
        let failure = RefCell::new(None);
        let response = |s: usize, alpha: f64| {
            let model = self.sections[s];
            match history {
                Some(states) => model.sample(states[s], alpha, gust[s]),
                None => model
                    .developed(alpha, gust[s])?
                    .response(&model.spec, model.zero_lift_angle_rad)
                    .map_err(airfoil_error),
            }
        };
        let closure = |s: usize, alpha: f64| match response(s, alpha) {
            Ok(channels) => (channels.cl, self.sections[s].regime(alpha)),
            Err(error) => {
                if failure.borrow().is_none() {
                    *failure.borrow_mut() = Some(error);
                }
                (f64::NAN, StripRegime::Blended)
            }
        };
        let result = self.propulsion.geometry().equilibrium(
            AerodynamicStage {
                closure: &closure,
                ..stage
            },
            warm,
        );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let coupled = result?;
        let local = coupled
            .strip_velocity_mps
            .as_deref()
            .ok_or_else(|| invalid("coupled force solve omitted its accepted local velocities"))?;
        if local.len() != self.sections.len()
            || coupled.strip_alpha_eff_rad.len() != self.sections.len()
        {
            return Err(invalid("incomplete accepted section fields"));
        }
        let mut channels = Vec::with_capacity(self.sections.len());
        let mut chordwise = Vec::with_capacity(self.sections.len());
        for (s, model) in self.sections.iter().enumerate() {
            channels.push(response(s, coupled.strip_alpha_eff_rad[s])?);
            chordwise.push(
                model
                    .frame
                    .flow(local[s])
                    .map_err(wing_error)?
                    .chordwise_mps,
            );
        }
        Ok(UnsteadyAerodynamics {
            coupled,
            sections: channels,
            gust_angle_rad: gust.to_vec(),
            chordwise_speed_mps: chordwise,
        })
    }

    fn evaluate_history(
        &self,
        motion: &GearedState,
        history: &[UnsteadySectionState],
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust: GustInput<'_>,
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyEvaluation<N>, Refusal> {
        self.admit_inputs(gust, warm)?;
        self.admit_history(history)?;
        let result = self
            .propulsion
            .evaluate_with(motion, inputs, wind_world, |stage| {
                let gust = gust.resolve(&stage, self.sections.len(), |s| self.sections[s].frame)?;
                self.solve(stage, Some(history), &gust, warm)
            })?;
        Ok(UnsteadyEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            aerodynamics: result.aerodynamics,
        })
    }

    /// Sample a fixed physical state. Repeated nonlinear solves do not age
    /// motion/gust/separation memory. Gust angles decompose the TOTAL wind
    /// already sampled by `wind_world`; they never inject a second velocity.
    ///
    /// # Errors
    /// Native force/history admission; the caller's state remains unchanged.
    pub fn evaluate(
        &self,
        state: &UnsteadyState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust_angle_rad: &[f64],
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyEvaluation<N>, Refusal> {
        self.evaluate_history(
            &state.motion,
            &state.sections,
            inputs,
            wind_world,
            GustInput::Angles(gust_angle_rad),
            warm,
        )
    }

    /// Evaluate native lift/moment and shaft loads with a field-derived gust
    /// partition. Each reference point is sampled once before nonlinear wing
    /// iteration; all iterations share the same frozen section history.
    pub fn evaluate_with_reference_wind(
        &self,
        state: &UnsteadyState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        reference_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyEvaluation<N>, Refusal> {
        self.evaluate_history(
            &state.motion,
            &state.sections,
            inputs,
            wind_world,
            GustInput::Reference(reference_world),
            warm,
        )
    }

    fn advance_history(
        &self,
        accepted: &[UnsteadySectionState],
        initial: &UnsteadyAerodynamics,
        midpoint: &UnsteadyAerodynamics,
        dt_s: f64,
    ) -> Result<Vec<UnsteadySectionState>, Refusal> {
        let mut next = accepted.to_vec();
        for (s, model) in self.sections.iter().enumerate() {
            let alpha_start = initial.coupled.strip_alpha_eff_rad[s];
            let alpha_midpoint = midpoint.coupled.strip_alpha_eff_rad[s];
            let ds = 2.0 * midpoint.chordwise_speed_mps[s] * dt_s / model.chord_m;
            next[s]
                .advance_midpoint(
                    &model.spec,
                    ds,
                    SectionMidpointInputs {
                        motion_start_rad: alpha_start - initial.gust_angle_rad[s],
                        gust_start_rad: initial.gust_angle_rad[s],
                        motion_midpoint_rad: alpha_midpoint - midpoint.gust_angle_rad[s],
                        gust_midpoint_rad: midpoint.gust_angle_rad[s],
                        effective_midpoint_rad: midpoint.sections[s].alpha_eff_rad,
                    },
                )
                .map_err(airfoil_error)?;
            next[s]
                .response(&model.spec, model.zero_lift_angle_rad)
                .map_err(airfoil_error)?;
        }
        Ok(next)
    }

    /// Advance native body/shaft motion and per-strip history with exactly two
    /// coupled force solves. The half predictor holds the accepted start input
    /// on the start chordwise clock. Full advancement uses midpoint chordwise
    /// speed, exact Duhamel ramps through start/midpoint inputs, and separation
    /// driven by the accepted midpoint's memory-effective geometric incidence,
    /// after both circulatory channels and before the camber zero-lift offset.
    ///
    /// The reconstructed ramp endpoint is a history reference, not a third
    /// solved force. Next evaluation corrects that reference with a jump on a
    /// copy. Smooth-input temporal accuracy depends on resolving flow variation;
    /// callers should split steps at prescribed wind/control discontinuities.
    /// `gust_angle` is called once per strip at each real stage, after controls.
    ///
    /// # Errors
    /// Every refusal leaves all caller-owned mechanical/history bits intact,
    /// including a late full-clock refusal after both force solves. Callback
    /// side effects remain the caller's responsibility.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W, G>(
        &self,
        state: &UnsteadyState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        inputs: C,
        wind_world: W,
        mut gust_angle: G,
    ) -> Result<UnsteadyStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
        G: FnMut(f64, &GearedState, usize) -> Result<f64, Refusal>,
    {
        self.step_with_evaluation(
            state,
            t_s,
            dt_s,
            angular_origin_world,
            warm,
            inputs,
            |time, mechanical, history, controls, stage_warm| {
                let gust = (0..self.sections.len())
                    .map(|s| gust_angle(time, mechanical, s))
                    .collect::<Result<Vec<_>, _>>()?;
                self.evaluate_history(
                    mechanical,
                    history,
                    controls,
                    &|point| wind_world(time, point),
                    GustInput::Angles(&gust),
                    stage_warm,
                )
            },
        )
    }

    /// Advance real body/shaft and section histories using total/reference
    /// fields at both actual mechanical stages. Control-stage geometry, rigid
    /// angular velocity and surface frames are shared by the decomposition and
    /// native load solve. Only the total field drives aerodynamic velocities.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reference_wind<C, W, R>(
        &self,
        state: &UnsteadyState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        inputs: C,
        wind_world: W,
        reference_world: R,
    ) -> Result<UnsteadyStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
        R: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        self.step_with_evaluation(
            state,
            t_s,
            dt_s,
            angular_origin_world,
            warm,
            inputs,
            |time, mechanical, history, controls, stage_warm| {
                self.evaluate_history(
                    mechanical,
                    history,
                    controls,
                    &|point| wind_world(time, point),
                    GustInput::Reference(&|point| reference_world(time, point)),
                    stage_warm,
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_evaluation<C, E>(
        &self,
        state: &UnsteadyState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        mut inputs: C,
        mut evaluate: E,
    ) -> Result<UnsteadyStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        E: FnMut(
            f64,
            &GearedState,
            &[UnsteadySectionState],
            PropulsionInputs,
            Option<[f64; 2]>,
        ) -> Result<UnsteadyEvaluation<N>, Refusal>,
    {
        self.admit_history(&state.sections)?;
        let mut initial = None;
        let mut midpoint = None;
        let mut stage_history = state.sections.clone();
        let mut stage_warm = warm;
        let motion = geared::step(
            self.propulsion.assembly,
            &state.motion,
            t_s,
            dt_s,
            angular_origin_world,
            |time, mechanical| {
                let controls = inputs(time, mechanical)?;
                let evaluation = evaluate(time, mechanical, &stage_history, controls, stage_warm)?;
                let loads = evaluation.loads;
                stage_warm = Some(evaluation.aerodynamics.coupled.w_slip);
                if initial.is_none() {
                    stage_history = self.advance_history(
                        &state.sections,
                        &evaluation.aerodynamics,
                        &evaluation.aerodynamics,
                        0.5 * dt_s,
                    )?;
                    initial = Some(evaluation);
                } else {
                    midpoint = Some(evaluation);
                }
                Ok(loads)
            },
        )?;
        let initial = initial.ok_or_else(|| invalid("missing initial force evaluation"))?;
        let midpoint = midpoint.ok_or_else(|| invalid("missing midpoint force evaluation"))?;
        let sections = self.advance_history(
            &state.sections,
            &initial.aerodynamics,
            &midpoint.aerodynamics,
            dt_s,
        )?;
        Ok(UnsteadyStep {
            state: UnsteadyState {
                motion: motion.state,
                sections,
            },
            motion,
            initial,
            midpoint,
        })
    }
}
