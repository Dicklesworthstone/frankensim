//! Explicit angular-wake/wing loads driving the existing geared integrator.

use super::{
    AerodynamicFields, AerodynamicSolution, Controls, GearedPropulsion, GearedState,
    PropulsionInputs, Refusal, SwirlEvaluation, SwirlStep, geared, reject,
};
use crate::propcoupling::spatial::DiskAxes;
use crate::propcoupling::swirl::wake::{
    self as coupling, CoupledWakeStep, WakeControls, WakeTransport, WakeWarmStart,
};
use fs_ga::Vec3;

impl AerodynamicSolution for CoupledWakeStep {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        self.solution.load_fields()
    }
}

/// Opt-in downstream angular transport, with admitted energy discrepancy.
/// Select material annuli to include native radial loading and continuity-driven
/// contraction/expansion; the fixed-cylinder projection remains the default.
/// Reuses the same actual stage atmosphere, prepared wing operator, section
/// frames, carrier wrench and signed shaft ports as ordinary geared propulsion.
/// See the coupling module for observation-window and noninteracting-disk limits.
pub struct WakePropulsion<'p, 'a, const N: usize> {
    propulsion: &'p GearedPropulsion<'a, N>,
    native: Controls,
    wake: WakeControls,
    transport: WakeTransport,
}
impl<'p, 'a, const N: usize> WakePropulsion<'p, 'a, N> {
    #[must_use]
    pub fn new(
        propulsion: &'p GearedPropulsion<'a, N>,
        native: Controls,
        wake: WakeControls,
    ) -> Self {
        Self {
            propulsion,
            native,
            wake,
            transport: WakeTransport::FixedCylinder,
        }
    }

    /// Select the source-field physics without replacing the mechanical owner.
    #[must_use]
    pub fn with_transport(mut self, transport: WakeTransport) -> Self {
        self.transport = transport;
        self
    }

    /// Evaluate real downstream angular wash before nonlinear circulation and
    /// complete force/moment assembly. Spin is derived from physical mounted
    /// shaft axes AND handedness; a negative gear ratio is not a thrust axis.
    ///
    /// # Errors
    /// Native stage, swirl, field-coupling, energy-budget and geometry refusals.
    pub fn evaluate(
        &self,
        state: &GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<&WakeWarmStart>,
    ) -> Result<SwirlEvaluation<N, CoupledWakeStep>, Refusal> {
        let p = self.propulsion;
        let axes = p
            .spatial_axes
            .unwrap_or(DiskAxes::new([[1.0, 0.0, 0.0]; 2])?);
        let spin = core::array::from_fn(|k| {
            let mount = p.mounts[k];
            let mechanical = p.assembly.rotors()[mount.rotor_index].axis_body;
            let thrust = mount.thrust_direction_body;
            let alignment = Vec3::new(mechanical[0], mechanical[1], mechanical[2])
                .dot(Vec3::new(thrust[0], thrust[1], thrust[2]));
            mount.handedness.sign() * alignment.signum()
        });
        let result = p.evaluate_with(state, inputs, wind_world, |stage| {
            coupling::coupled_step_with_transport(
                self.transport,
                axes,
                spin,
                stage.op,
                stage.panels,
                stage.strips,
                stage.closure,
                stage.section_frames,
                &p.design.rotor,
                stage.disks,
                stage.freestream,
                stage.rho,
                &p.design.coupling,
                self.native,
                self.wake,
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

    /// Advance actual carrier/shaft momentum using two complete angular-field
    /// solves at the initial and midpoint states. Internal motor torque and
    /// external aerodynamic shaft torque retain their existing exclusive owner.
    /// The field affects the wing wrench; it is never added as a second
    /// artificial shaft reaction on the carrier.
    ///
    /// # Errors
    /// No motion is published if either aerodynamic or mechanical stage refuses.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        state: &GearedState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<&WakeWarmStart>,
        mut inputs: C,
        wind_world: W,
    ) -> Result<SwirlStep<N, CoupledWakeStep>, Refusal>
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
                .ok_or_else(|| reject("swirl-wake-invalid", "missing initial stage"))?,
            midpoint: midpoint
                .ok_or_else(|| reject("swirl-wake-invalid", "missing midpoint stage"))?,
        })
    }
}
