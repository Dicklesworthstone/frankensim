//! Steady swirling rotors driving the actual geared six-DOF carrier.
//!
//! Reuses the admitted mount geometry, material-air sampling, wing force
//! assembly and mechanical shaft ports of `GearedPropulsion`. The selected
//! native rotor solve changes thrust, aerodynamic shaft torque and their
//! resulting translation, spin and gyroscopic motion. No separate rotor
//! reaction torque is added to the carrier a second time.

use super::{
    AerodynamicFields, AerodynamicSolution, GearedLoads, GearedPropulsion, GearedState, GearedStep,
    PropulsionInputs, Refusal, geared, reject,
};
use crate::propcoupling::swirl::{self as coupling, CoupledSwirlStep, WarmStart};
use fs_airscrew::swirl::Controls;

/// Real accepted wing/swirl loads on the carrier and physical shaft ports.
#[derive(Clone, Debug, PartialEq)]
pub struct SwirlEvaluation<const N: usize> {
    /// Actual forces/moments and signed aerodynamic torques for the gear owner.
    pub loads: GearedLoads<N>,
    /// Body aerodynamic force, before world gravity is added [N].
    pub aerodynamic_force_body_n: [f64; 3],
    /// Signed absolute physical L/R rotor speeds [rad/s].
    pub propeller_absolute_speed_rad_s: [f64; 2],
    /// Same accepted native annular momentum and two-way wing solution.
    pub coupled: CoupledSwirlStep,
}

/// Accepted carrier/shaft motion and its real aerodynamic stage evaluations.
#[derive(Clone, Debug, PartialEq)]
pub struct SwirlStep<const N: usize> {
    /// Existing geared owner, including actual global momentum/work residuals.
    pub motion: GearedStep<N>,
    /// Start-of-step aerodynamic evaluation for the midpoint predictor.
    pub initial: SwirlEvaluation<N>,
    /// Accepted midpoint aerodynamic loads used by the mechanical step.
    pub midpoint: SwirlEvaluation<N>,
}

/// Explicit swirl selection on an existing admitted physical assembly.
/// It does not change an ordinary propulsion model or a historical preset.
pub struct SwirlPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p GearedPropulsion<'a, N>,
    controls: Controls,
}

impl AerodynamicSolution for CoupledSwirlStep {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        self.coupled.load_fields()
    }
}

impl<'p, 'a, const N: usize> SwirlPropulsion<'p, 'a, N> {
    /// Select the native bounded swirl solve. Its controls are admitted by
    /// the native solver before any load evaluation can be returned.
    #[must_use]
    pub fn new(propulsion: &'p GearedPropulsion<'a, N>, controls: Controls) -> Self {
        Self {
            propulsion,
            controls,
        }
    }

    /// Evaluate the actual swirling-rotor/wing solution at one physical state.
    /// Both handedness and carrier spin enter the existing absolute shaft
    /// speed. Required BEMT torque becomes the external aerodynamic shaft
    /// torque exactly once, retaining regeneration and gyroscopic effects.
    ///
    /// Wing/inter-disk wash remains the existing axial actuator-disk field;
    /// this steady selection does not combine swirl with dynamic wake modes.
    ///
    /// # Errors
    /// Existing assembly, material-air, wing, native swirl and two-way
    /// coupling refusals pass through without a partial load state.
    pub fn evaluate(
        &self,
        state: &GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<&WarmStart>,
    ) -> Result<SwirlEvaluation<N>, Refusal> {
        let p = self.propulsion;
        let result = p.evaluate_with(state, inputs, wind_world, |stage| {
            coupling::coupled_step(
                p.spatial_axes,
                stage.op,
                stage.panels,
                stage.strips,
                stage.closure,
                &p.design.rotor,
                stage.disks,
                stage.freestream,
                stage.rho,
                &p.design.coupling,
                self.controls,
                warm,
                stage.incident,
            )
        })?;
        Ok(SwirlEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            coupled: result.aerodynamics,
        })
    }

    /// Advance translation, carrier orientation and geared rotor spin using
    /// real initial and midpoint swirl/wing solves. Numerical branch seeds
    /// come only from the previous accepted aerodynamic stage. Physical
    /// angular momentum remains exclusively in the geared state.
    ///
    /// # Errors
    /// All model/integrator/callback refusals pass through; no partial motion
    /// is published. Caller callback side effects are outside that guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        state: &GearedState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<&WarmStart>,
        mut inputs: C,
        wind_world: W,
    ) -> Result<SwirlStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut initial = None;
        let mut midpoint = None;
        let mut stage_warm = warm.cloned();
        let motion = geared::step(
            self.propulsion.assembly,
            state,
            t_s,
            dt_s,
            angular_origin_world,
            |time, stage| {
                let evaluation = self.evaluate(
                    stage,
                    inputs(time, stage)?,
                    &|point| wind_world(time, point),
                    stage_warm.as_ref(),
                )?;
                stage_warm = Some(evaluation.coupled.warm_start());
                let loads = evaluation.loads;
                if initial.is_none() {
                    initial = Some(evaluation);
                } else {
                    midpoint = Some(evaluation);
                }
                Ok(loads)
            },
        )?;
        Ok(SwirlStep {
            motion,
            initial: initial
                .ok_or_else(|| reject("swirl-coupling-invalid", "missing initial stage"))?,
            midpoint: midpoint
                .ok_or_else(|| reject("swirl-coupling-invalid", "missing midpoint stage"))?,
        })
    }
}
