//! Evolving propeller induction composed with per-strip unsteady lift memory.
//!
//! This is the native composition of the dynamic-inflow owner with the
//! Wagner/Kussner/separation section owner. The accepted state contains all
//! three physical memories: body/shaft motion, disk induction, and section
//! histories. Nonlinear wing iterations sample copies of section history;
//! only the two real midpoint stages advance wake and section state.

use super::{DynamicAerodynamics, DynamicPropulsion, DynamicState};
use super::super::{
    AerodynamicFields, AerodynamicSolution, AerodynamicStage, GearedLoads, GearedStep,
    PropulsionInputs, Refusal, geared, reject, wing_error,
};
use fs_airfoil::unsteady::{
    DuhamelState, SectionChannels, SectionMidpointInputs, SeparationLagState, UnsteadySectionSpec,
    UnsteadySectionState,
};
use fs_airscrew::dynamic_inflow::InflowStep;
use fs_wing::nonlinear::{SectionFrame, StripRegime};
use std::cell::RefCell;

use super::super::unsteady::MAX_EFFECTIVE_ALPHA_RAD;

/// Joint accepted state for body/shaft, dynamic inflow, and wing-section memory.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicUnsteadyState {
    pub dynamic: DynamicState,
    /// One motion/gust/separation history per physical lifting strip.
    pub sections: Vec<UnsteadySectionState>,
}

/// Dynamic-wake loads evaluated with one frozen set of section histories.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicUnsteadyAerodynamics {
    pub dynamic: DynamicAerodynamics,
    pub sections: Vec<SectionChannels>,
    pub gust_angle_rad: Vec<f64>,
    pub chordwise_speed_mps: Vec<f64>,
}

impl AerodynamicSolution for DynamicUnsteadyAerodynamics {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        self.dynamic.load_fields()
    }

    fn section_channels(&self) -> Option<&[SectionChannels]> {
        Some(&self.sections)
    }
}

/// One real force evaluation without publishing any physical-memory update.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicUnsteadyEvaluation<const N: usize> {
    pub loads: GearedLoads<N>,
    pub aerodynamic_force_body_n: [f64; 3],
    pub propeller_absolute_speed_rad_s: [f64; 2],
    pub aerodynamics: DynamicUnsteadyAerodynamics,
}

/// Atomically accepted body/shaft, wake, and section histories.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicUnsteadyStep<const N: usize> {
    pub state: DynamicUnsteadyState,
    pub motion: GearedStep<N>,
    pub inflow: [InflowStep; 2],
    pub initial: DynamicUnsteadyEvaluation<N>,
    pub midpoint: DynamicUnsteadyEvaluation<N>,
}

#[derive(Clone, Copy)]
struct SectionModel {
    spec: UnsteadySectionSpec,
    zero_lift_angle_rad: f64,
    chord_m: f64,
    frame: SectionFrame,
}

fn airfoil_error(error: fs_airfoil::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

fn invalid(message: &str) -> Refusal {
    reject("dynamic-unsteady-propulsion-invalid", message)
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

/// Opt-in composition of native dynamic inflow and native unsteady wing lift.
pub struct DynamicUnsteadyPropulsion<'d, 'p, 'a, const N: usize> {
    dynamic: &'d DynamicPropulsion<'p, 'a, N>,
    sections: Vec<SectionModel>,
}

impl<'d, 'p, 'a, const N: usize> DynamicUnsteadyPropulsion<'d, 'p, 'a, N> {
    /// Prepare immutable section metadata over an already admitted dynamic wake.
    ///
    /// # Errors
    /// Native geometry/spec admission or a nonpositive attached lift slope.
    pub fn new(
        dynamic: &'d DynamicPropulsion<'p, 'a, N>,
        spec: UnsteadySectionSpec,
    ) -> Result<Self, Refusal> {
        spec.admit().map_err(airfoil_error)?;
        if spec.cl_alpha <= 0.0 {
            return Err(invalid("the attached lift slope must be positive"));
        }
        let propulsion = dynamic.propulsion;
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
        Ok(Self { dynamic, sections })
    }

    #[must_use]
    pub fn section_count(&self) -> usize {
        self.sections.len()
    }

    fn admit_inputs(&self, gust: &[f64]) -> Result<(), Refusal> {
        if gust.len() != self.sections.len() || !gust.iter().all(|a| a.is_finite()) {
            return Err(invalid("one finite explicit gust angle per strip is required"));
        }
        Ok(())
    }

    fn admit_history(&self, history: &[UnsteadySectionState]) -> Result<(), Refusal> {
        if history.len() != self.sections.len() {
            return Err(invalid("one complete physical history per strip is required"));
        }
        for (state, model) in history.iter().zip(&self.sections) {
            state
                .response(&model.spec, model.zero_lift_angle_rad)
                .map_err(airfoil_error)?;
        }
        Ok(())
    }

    /// Initialize developed section memory at the actual stored dynamic wake.
    pub fn trimmed_state(
        &self,
        dynamic: DynamicState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust_angle_rad: &[f64],
    ) -> Result<DynamicUnsteadyState, Refusal> {
        self.admit_inputs(gust_angle_rad)?;
        self.dynamic
            .wake
            .admit_state(inputs.density_kg_m3, dynamic.induced_mps)?;
        let result = self.dynamic.propulsion.evaluate_with(
            &dynamic.motion,
            inputs,
            wind_world,
            |stage| self.solve(stage, dynamic.induced_mps, None, gust_angle_rad),
        )?;
        let history = result
            .aerodynamics
            .dynamic
            .strip_alpha_eff_rad
            .iter()
            .enumerate()
            .map(|(s, &alpha)| self.sections[s].developed(alpha, gust_angle_rad[s]))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(DynamicUnsteadyState {
            dynamic,
            sections: history,
        })
    }

    fn solve(
        &self,
        stage: AerodynamicStage<'_>,
        induced_mps: [f64; 2],
        history: Option<&[UnsteadySectionState]>,
        gust: &[f64],
    ) -> Result<DynamicUnsteadyAerodynamics, Refusal> {
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
        let dynamic = self.dynamic.wake.solve(
            AerodynamicStage {
                closure: &closure,
                ..stage
            },
            induced_mps,
        );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let dynamic = dynamic?;
        if dynamic.strip_alpha_eff_rad.len() != self.sections.len()
            || dynamic.strip_velocity_mps.len() != self.sections.len()
        {
            return Err(invalid("incomplete accepted section fields"));
        }
        let mut channels = Vec::with_capacity(self.sections.len());
        let mut chordwise = Vec::with_capacity(self.sections.len());
        for (s, model) in self.sections.iter().enumerate() {
            channels.push(response(s, dynamic.strip_alpha_eff_rad[s])?);
            chordwise.push(
                model
                    .frame
                    .flow(dynamic.strip_velocity_mps[s])
                    .map_err(wing_error)?
                    .chordwise_mps,
            );
        }
        Ok(DynamicUnsteadyAerodynamics {
            dynamic,
            sections: channels,
            gust_angle_rad: gust.to_vec(),
            chordwise_speed_mps: chordwise,
        })
    }

    fn evaluate_history(
        &self,
        dynamic: &DynamicState,
        history: &[UnsteadySectionState],
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust: &[f64],
    ) -> Result<DynamicUnsteadyEvaluation<N>, Refusal> {
        self.admit_inputs(gust)?;
        self.admit_history(history)?;
        self.dynamic
            .wake
            .admit_state(inputs.density_kg_m3, dynamic.induced_mps)?;
        let result = self.dynamic.propulsion.evaluate_with(
            &dynamic.motion,
            inputs,
            wind_world,
            |stage| self.solve(stage, dynamic.induced_mps, Some(history), gust),
        )?;
        Ok(DynamicUnsteadyEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            aerodynamics: result.aerodynamics,
        })
    }

    /// Sample the joint accepted state without advancing either wake or section memory.
    pub fn evaluate(
        &self,
        state: &DynamicUnsteadyState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        gust_angle_rad: &[f64],
    ) -> Result<DynamicUnsteadyEvaluation<N>, Refusal> {
        self.evaluate_history(
            &state.dynamic,
            &state.sections,
            inputs,
            wind_world,
            gust_angle_rad,
        )
    }

    fn advance_history(
        &self,
        accepted: &[UnsteadySectionState],
        initial: &DynamicUnsteadyAerodynamics,
        midpoint: &DynamicUnsteadyAerodynamics,
        dt_s: f64,
    ) -> Result<Vec<UnsteadySectionState>, Refusal> {
        let mut next = accepted.to_vec();
        for (s, model) in self.sections.iter().enumerate() {
            let alpha_start = initial.dynamic.strip_alpha_eff_rad[s];
            let alpha_midpoint = midpoint.dynamic.strip_alpha_eff_rad[s];
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

    /// Advance body/shaft, dynamic inflow, and section histories on one midpoint clock.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W, G>(
        &self,
        state: &DynamicUnsteadyState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
        mut gust_angle: G,
    ) -> Result<DynamicUnsteadyStep<N>, Refusal>
    where
        C: FnMut(f64, &DynamicUnsteadyState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
        G: FnMut(f64, &DynamicUnsteadyState, usize) -> Result<f64, Refusal>,
    {
        self.admit_history(&state.sections)?;
        let mut initial = None;
        let mut midpoint = None;
        let mut stage_induction = state.dynamic.induced_mps;
        let mut stage_history = state.sections.clone();
        let motion = geared::step(
            self.dynamic.propulsion.assembly,
            &state.dynamic.motion,
            t_s,
            dt_s,
            angular_origin_world,
            |time, mechanical| {
                let stage_state = DynamicUnsteadyState {
                    dynamic: DynamicState {
                        motion: *mechanical,
                        induced_mps: stage_induction,
                    },
                    sections: stage_history.clone(),
                };
                let controls = inputs(time, &stage_state)?;
                let gust = (0..self.sections.len())
                    .map(|s| gust_angle(time, &stage_state, s))
                    .collect::<Result<Vec<_>, _>>()?;
                let evaluation = self.evaluate_history(
                    &stage_state.dynamic,
                    &stage_history,
                    controls,
                    &|point| wind_world(time, point),
                    &gust,
                )?;
                let loads = evaluation.loads;
                if initial.is_none() {
                    stage_induction = self
                        .dynamic
                        .wake
                        .advance_wakes(
                            state.dynamic.induced_mps,
                            &evaluation.aerodynamics.dynamic,
                            0.5 * dt_s,
                        )?
                        .map(|wake| wake.induced_mps);
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
        let inflow = self.dynamic.wake.advance_wakes(
            state.dynamic.induced_mps,
            &midpoint.aerodynamics.dynamic,
            dt_s,
        )?;
        let sections = self.advance_history(
            &state.sections,
            &initial.aerodynamics,
            &midpoint.aerodynamics,
            dt_s,
        )?;
        Ok(DynamicUnsteadyStep {
            state: DynamicUnsteadyState {
                dynamic: DynamicState {
                    motion: motion.state,
                    induced_mps: inflow.map(|wake| wake.induced_mps),
                },
                sections,
            },
            motion,
            inflow,
            initial,
            midpoint,
        })
    }
}
