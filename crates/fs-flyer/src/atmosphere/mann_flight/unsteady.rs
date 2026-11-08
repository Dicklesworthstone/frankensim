//! The original Mann stages drive native motion/gust/separation load memory.

use super::{MannFlightBoundary, MannFlightWind, invalid};
use crate::Refusal;
use crate::spine::rotors::geared::GearedState;
use crate::spine::rotors::propulsion::PropulsionInputs;
use crate::spine::rotors::propulsion::dynamic::DynamicState;
use crate::spine::rotors::propulsion::dynamic::unsteady::{
    DynamicUnsteadyPropulsion, DynamicUnsteadyState, DynamicUnsteadyStep,
};
use crate::spine::rotors::propulsion::unsteady::{UnsteadyPropulsion, UnsteadyState, UnsteadyStep};
use fs_atmo::{FlatSiteLogLaw, FlowSample};

fn mean_world(
    mean: FlatSiteLogLaw,
    boundary: MannFlightBoundary,
    ground_down_m: f64,
    point: [f64; 3],
) -> Result<FlowSample, Refusal> {
    let height = ground_down_m - point[2];
    let minimum = match boundary {
        MannFlightBoundary::HomogeneousAbove { minimum_height_m } => minimum_height_m,
        MannFlightBoundary::GroundPlane => 0.0,
    };
    if !point.iter().all(|x| x.is_finite()) || !height.is_finite() || height < minimum {
        return Err(invalid(
            "mann-flight-applicability",
            "reference point lies outside the wind's admitted height domain",
        ));
    }
    let speed = mean.speed(height);
    let gradient = mean.dspeed_dh(height);
    if !speed.is_finite() || !gradient.is_finite() {
        return Err(invalid(
            "mann-flight-invalid",
            "reference wind or gradient is not representable",
        ));
    }
    Ok(FlowSample {
        u: [-speed, 0.0, 0.0],
        grad: [[0.0, 0.0, gradient], [0.0; 3], [0.0; 3]],
    })
}

impl MannFlightWind {
    /// Original stationary log-law reference, in the same world frame and
    /// height domain as `sample`. No OU quadratures are evaluated or advanced.
    /// The unsteady native load owners use its incidence to partition gust
    /// and motion, while `sample` remains the complete physical wind.
    pub fn sample_mean(&self, position_world_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        mean_world(
            self.mean,
            self.boundary,
            self.ground_down_m,
            position_world_m,
        )
    }

    /// Initialize developed native wing memory from this accepted realization.
    /// The actual same-point turbulent incidence occupies the gust channel;
    /// reference mean, rigid motion, controls and induced flow retain theirs.
    pub fn trimmed_unsteady_state<const N: usize>(
        &self,
        model: &UnsteadyPropulsion<'_, '_, N>,
        motion: GearedState,
        inputs: PropulsionInputs,
        warm: Option<[f64; 2]>,
    ) -> Result<UnsteadyState, Refusal> {
        model.trimmed_state_with_reference_wind(
            motion,
            inputs,
            &|point| self.sample(point).map(|flow| flow.u),
            &|point| self.sample_mean(point).map(|flow| flow.u),
            warm,
        )
    }

    /// Advance the original body/shaft and unsteady wing-force owner with this
    /// realization's start and midpoint fields. Every strip uses the identical
    /// total/reference geometry and stage; no scalar gust is inferred at the CG
    /// or added a second time. Histories and wind publish together after success.
    /// The native attached/stall-onset validity and chordwise clocks are retained.
    #[allow(clippy::too_many_arguments)]
    pub fn step_unsteady<const N: usize, C>(
        &mut self,
        model: &UnsteadyPropulsion<'_, '_, N>,
        state: &UnsteadyState,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        maximum_mode_transitions: usize,
        cancelled: impl FnMut() -> bool,
        inputs: C,
    ) -> Result<UnsteadyStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
    {
        self.with_stages(maximum_mode_transitions, cancelled, |stages| {
            let result = model.step_with_reference_wind(
                state,
                stages.times[0],
                stages.dt_s,
                angular_origin_world,
                warm,
                inputs,
                |time, point| stages.sample(time, point).map(|flow| flow.u),
                |_, point| {
                    mean_world(stages.mean, stages.boundary, stages.ground_down_m, point)
                        .map(|flow| flow.u)
                },
            )?;
            stages.sample(stages.times[2], result.state.motion.rigid.pos_m)?;
            Ok(result)
        })
    }

    /// Initialize the same developed wing memory around an already selected
    /// native dynamic propeller-wake state. No induction equilibrium is imposed.
    pub fn trimmed_dynamic_unsteady_state<const N: usize>(
        &self,
        model: &DynamicUnsteadyPropulsion<'_, '_, '_, N>,
        dynamic: DynamicState,
        inputs: PropulsionInputs,
    ) -> Result<DynamicUnsteadyState, Refusal> {
        model.trimmed_state_with_reference_wind(
            dynamic,
            inputs,
            &|point| self.sample(point).map(|flow| flow.u),
            &|point| self.sample_mean(point).map(|flow| flow.u),
        )
    }

    /// Coherent evolving turbulence through real body/shaft, disk-induction and
    /// unsteady wing histories. The original dynamic blade sampling selection
    /// is retained, so every selected blade point shares these same OU stages.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dynamic_unsteady<const N: usize, C>(
        &mut self,
        model: &DynamicUnsteadyPropulsion<'_, '_, '_, N>,
        state: &DynamicUnsteadyState,
        angular_origin_world: [f64; 3],
        maximum_mode_transitions: usize,
        cancelled: impl FnMut() -> bool,
        inputs: C,
    ) -> Result<DynamicUnsteadyStep<N>, Refusal>
    where
        C: FnMut(f64, &DynamicUnsteadyState) -> Result<PropulsionInputs, Refusal>,
    {
        self.with_stages(maximum_mode_transitions, cancelled, |stages| {
            let result = model.step_with_reference_wind(
                state,
                stages.times[0],
                stages.dt_s,
                angular_origin_world,
                inputs,
                |time, point| stages.sample(time, point).map(|flow| flow.u),
                |_, point| {
                    mean_world(stages.mean, stages.boundary, stages.ground_down_m, point)
                        .map(|flow| flow.u)
                },
            )?;
            stages.sample(stages.times[2], result.state.dynamic.motion.rigid.pos_m)?;
            Ok(result)
        })
    }
}
