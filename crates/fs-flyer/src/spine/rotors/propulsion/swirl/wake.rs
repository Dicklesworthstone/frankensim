//! Explicit angular-wake/wing loads sharing native aerodynamic evaluation.

use super::super::{AerodynamicStage, PropulsionGeometry};
use super::{
    AerodynamicFields, AerodynamicSolution, Controls, GearedPropulsion, GearedState,
    PropulsionInputs, Refusal, SwirlEvaluation, SwirlStep, geared, reject,
};
use crate::propcoupling::spatial::DiskAxes;
use crate::propcoupling::swirl::wake::{
    self as coupling, ContinuedWakeStep, CoupledWakeStep, WakeControls, WakeTransport,
    WakeWarmStart,
};
use fs_ga::Vec3;
use fs_wing::nonlinear::continuation::{ContinuationControls, WingContinuation};

impl AerodynamicSolution for CoupledWakeStep {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        self.solution.load_fields()
    }
}

impl AerodynamicSolution for ContinuedWakeStep {
    fn load_fields(&self) -> AerodynamicFields<'_> {
        self.solution.load_fields()
    }
}

// Both mechanical consumers share one admitted source-field solve. Geometry
// and physical spin signs are prepared once; no gear coordinate enters here.
pub(in super::super) struct WakeModel<'a> {
    geometry: PropulsionGeometry<'a>,
    axes: DiskAxes,
    spin: [f64; 2],
    native: Controls,
    wake: WakeControls,
    transport: WakeTransport,
}

impl<'a> WakeModel<'a> {
    pub(in super::super) fn new(
        geometry: PropulsionGeometry<'a>,
        mechanical_axes: [[f64; 3]; 2],
        native: Controls,
        wake: WakeControls,
    ) -> Self {
        let axes = geometry
            .spatial_axes
            .unwrap_or_else(|| DiskAxes::new([[1.0, 0.0, 0.0]; 2]).expect("constant unit axes"));
        let spin = core::array::from_fn(|k| {
            let mount = geometry.mounts[k];
            let mechanical = mechanical_axes[k];
            let thrust = mount.thrust_direction_body;
            let alignment = Vec3::new(mechanical[0], mechanical[1], mechanical[2])
                .dot(Vec3::new(thrust[0], thrust[1], thrust[2]));
            mount.handedness.sign() * alignment.signum()
        });
        Self {
            geometry,
            axes,
            spin,
            native,
            wake,
            transport: WakeTransport::FixedCylinder,
        }
    }

    pub(in super::super) fn with_transport(mut self, transport: WakeTransport) -> Self {
        self.transport = transport;
        self
    }

    pub(in super::super) fn solve(
        &self,
        stage: AerodynamicStage<'_>,
        warm: Option<&WakeWarmStart>,
    ) -> Result<CoupledWakeStep, Refusal> {
        coupling::coupled_step_with_transport(
            self.transport,
            self.axes,
            self.spin,
            stage.op,
            stage.panels,
            stage.strips,
            stage.closure,
            stage.section_frames,
            &self.geometry.design.rotor,
            stage.disks,
            stage.freestream,
            stage.rho,
            &self.geometry.design.coupling,
            self.native,
            self.wake,
            warm,
            stage.incident,
        )
    }
    pub(in super::super) fn solve_continued(
        &self,
        stage: AerodynamicStage<'_>,
        warm: Option<&WakeWarmStart>,
        previous: Option<&WingContinuation>,
        domain: ContinuationControls,
    ) -> Result<ContinuedWakeStep, Refusal> {
        coupling::coupled_step_continued(
            self.transport,
            self.axes,
            self.spin,
            stage.op,
            stage.panels,
            stage.strips,
            stage.closure,
            stage.section_frames,
            &self.geometry.design.rotor,
            stage.disks,
            stage.freestream,
            stage.rho,
            &self.geometry.design.coupling,
            self.native,
            self.wake,
            warm,
            stage.incident,
            previous,
            domain,
        )
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
    model: WakeModel<'a>,
}
impl<'p, 'a, const N: usize> WakePropulsion<'p, 'a, N> {
    #[must_use]
    pub fn new(
        propulsion: &'p GearedPropulsion<'a, N>,
        native: Controls,
        wake: WakeControls,
    ) -> Self {
        let mechanical_axes = propulsion
            .mounts
            .map(|mount| propulsion.assembly.rotors()[mount.rotor_index].axis_body);
        Self {
            propulsion,
            model: WakeModel::new(propulsion.geometry(), mechanical_axes, native, wake),
        }
    }

    /// Select the source-field physics without replacing the mechanical owner.
    #[must_use]
    pub fn with_transport(mut self, transport: WakeTransport) -> Self {
        self.model = self.model.with_transport(transport);
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
        let result = self
            .propulsion
            .evaluate_with(state, inputs, wind_world, |stage| {
                self.model.solve(stage, warm)
            })?;
        Ok(SwirlEvaluation {
            loads: result.loads,
            aerodynamic_force_body_n: result.aerodynamic_force_body_n,
            propeller_absolute_speed_rad_s: result.propeller_absolute_speed_rad_s,
            coupled: result.aerodynamics,
        })
    }

    /// Evaluate the native wing branch under the complete annular field.
    /// `previous=None` establishes an initial complete accepted branch.
    /// Numerical field seeds and the caller's fixed wing anchor are distinct.
    ///
    /// # Errors
    /// Native field and wing-continuation refusals leave both inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_continued(
        &self,
        state: &GearedState,
        inputs: PropulsionInputs,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<&WakeWarmStart>,
        previous: Option<&WingContinuation>,
        domain: ContinuationControls,
    ) -> Result<SwirlEvaluation<N, ContinuedWakeStep>, Refusal> {
        let result = self
            .propulsion
            .evaluate_with(state, inputs, wind_world, |stage| {
                self.model.solve_continued(stage, warm, previous, domain)
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
        inputs: C,
        wind_world: W,
    ) -> Result<SwirlStep<N, CoupledWakeStep>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut stage_warm = warm.cloned();
        self.step_with(
            state,
            t_s,
            dt_s,
            angular_origin_world,
            inputs,
            wind_world,
            |stage, inputs, wind| {
                let evaluation = self.evaluate(stage, inputs, wind, stage_warm.as_ref())?;
                stage_warm = Some(evaluation.coupled.warm_start());
                Ok(evaluation)
            },
        )
    }

    /// Advance from accepted initial/midpoint wing branches. The complete
    /// midpoint anchor is returned with its actual loads only after mechanics
    /// succeeds; retain it with the returned geared state for the next step.
    ///
    /// # Errors
    /// Any native field, wing or mechanical refusal publishes neither state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_continued<C, W>(
        &self,
        state: &GearedState,
        wing: &WingContinuation,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<&WakeWarmStart>,
        domain: ContinuationControls,
        inputs: C,
        wind_world: W,
    ) -> Result<SwirlStep<N, ContinuedWakeStep>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut stage_warm = warm.cloned();
        let mut stage_wing = wing.clone();
        self.step_with(
            state,
            t_s,
            dt_s,
            angular_origin_world,
            inputs,
            wind_world,
            |stage, inputs, wind| {
                let evaluation = self.evaluate_continued(
                    stage,
                    inputs,
                    wind,
                    stage_warm.as_ref(),
                    Some(&stage_wing),
                    domain,
                )?;
                stage_warm = Some(evaluation.coupled.solution.warm_start());
                stage_wing = evaluation.coupled.wing_continuation.clone();
                Ok(evaluation)
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with<S, C, W, E>(
        &self,
        state: &GearedState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        mut inputs: C,
        wind_world: W,
        mut evaluate: E,
    ) -> Result<SwirlStep<N, S>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<PropulsionInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
        E: FnMut(
            &GearedState,
            PropulsionInputs,
            &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        ) -> Result<SwirlEvaluation<N, S>, Refusal>,
    {
        let mut initial = None;
        let mut midpoint = None;
        let motion = geared::step(
            self.propulsion.assembly,
            state,
            t_s,
            dt_s,
            angular_origin_world,
            |time, stage| {
                let evaluation = evaluate(stage, inputs(time, stage)?, &|point| {
                    wind_world(time, point)
                })?;
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
