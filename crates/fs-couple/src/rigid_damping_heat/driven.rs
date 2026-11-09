//! Driven spatial damping on the original finite thermal timeline.
//!
//! This pair entry adds simultaneous held actuation to the original native
//! damper. Only its independently integrated positive viscous loss enters the
//! finite caloric receivers. Signed actuator supply/recovery remains work.
//! The same thermal events, transport, observers and whole-window energy gates
//! as the passive network entry remain authoritative.
//!
//! Exactly two original bodies and one canonical [0, 1] connector are admitted.
//! This exposes the native driven-pair owner; it does not approximate a driven
//! network by sequentially applying independent actuators.

/// Simultaneous arbitrary driven connector graphs with finite heat storage.
pub mod network;

use super::{
    FrozenSource, MechanicalAdvance, MechanicalStep, Result, RigidDampingThermalCandidate,
    RigidDampingThermalError, RigidDampingThermalInput, ThermalControls, advance_with_source,
    temperature::{TemperatureDampingLaw, TemperatureDampingThermalCandidate, advance_with_motion},
};
use fs_mbd::{
    RigidBodyState,
    spatial_loads::{
        internal_impulse::{
            InternalWrenchImpulseError,
            viscous::{
                ViscousConnectorControls, ViscousConnectorError,
                driven::DrivenViscousConnectorError, network::ViscousNetworkLink,
            },
        },
        internal_motion::{
            InternalWrenchStepError,
            viscous::{
                ViscousMotionLaw,
                driven::{DrivenViscousMotionError, step_driven_viscous_connector_with_loads},
            },
        },
    },
};
use fs_thermal::enthalpy::{CellId, DepositionObserver, DepositionSimulation, StepConfig};

pub use fs_mbd::spatial_loads::internal_motion::viscous::driven::HeldInternalEffort;

/// Native fixed-size pair controls and original thermal timeline limits.
///
/// Each thermal interval performs exactly two six-coordinate source
/// preparations and two native body drifts. There is no graph subdivision or
/// growing mechanical state; native source scratch is fixed by that dimension.
/// Thermal cells, events and step count retain their original owner's budget.
#[derive(Clone, Copy, Debug)]
pub struct DrivenDampingThermalConfig {
    pub mechanical: ViscousConnectorControls,
    /// Also bounds the moving-geometry and lagged-temperature intervals.
    pub thermal: StepConfig,
    pub maximum_heat_bindings: usize,
    /// One allowance for the complete accepted window, independently applied
    /// to mechanical balance, combined balance and absolute heat-handoff error.
    pub energy_tolerance_j: f64,
}
impl From<DrivenDampingThermalConfig> for ThermalControls {
    fn from(config: DrivenDampingThermalConfig) -> Self {
        Self {
            thermal: config.thermal,
            maximum_heat_bindings: config.maximum_heat_bindings,
            energy_tolerance_j: config.energy_tolerance_j,
        }
    }
}

impl<E, S> From<DrivenViscousMotionError> for RigidDampingThermalError<E, S> {
    fn from(error: DrivenViscousMotionError) -> Self {
        match error {
            DrivenViscousMotionError::Motion(InternalWrenchStepError::Cancelled)
            | DrivenViscousMotionError::Source(DrivenViscousConnectorError::Source(
                ViscousConnectorError::Impulse(InternalWrenchImpulseError::Cancelled),
            )) => Self::Cancelled,
            error => Self::DrivenMotion(error),
        }
    }
}

/// Advance simultaneous held actuation and damping with finite thermal storage.
///
/// Input banks contain the original two bodies and a single link whose body
/// indices are [0, 1]. The effort is applied to attachment B in its declared
/// world or first-material frame; the original native owner transports the
/// opposite wrench and re-expresses world effort after finite motion. It is
/// held over the entire call, including event-clipped thermal intervals.
///
/// Each interval's native actuator work is reported independently from external
/// applied work and physical damping heat. Recovery may make actuator work
/// negative while all deposited damper heat remains nonnegative. Neither work
/// nor a numerical defect is deposited into the thermal owner.
///
/// Inputs are borrowed and a complete mechanical/thermal/observer candidate is
/// returned only after all native and whole-window gates pass. Refine time for
/// changing geometry, thermal transport and source-averaging accuracy.
pub fn advance_candidate<'m, O: DepositionObserver>(
    input: RigidDampingThermalInput<'_>,
    effort: HeldInternalEffort,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    config: DrivenDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<RigidDampingThermalCandidate<'m, O>, O::Error> {
    advance_with_source(
        input,
        thermal,
        observer,
        FrozenSource(input.links),
        DrivenMotion {
            effort,
            controls: config.mechanical,
        },
        target_s,
        config.into(),
        cancelled,
    )
    .map(|(candidate, _)| candidate)
}

/// The same driven pair with original finite-cell temperature feedback.
///
/// The single supplied temperature coordinate is independent of the heat
/// partition. The existing TemperatureDampingLaw receives the original
/// reference matrix and actual caloric temperature at each interval start and
/// at the accepted endpoint. CalibratedKelvinDamping is directly usable here;
/// no calibration, interface temperature or reversible parameter work is
/// invented. The effort itself remains held, independent of temperature.
///
/// This shares the original first-order lagged-temperature discretization.
/// Constant-coefficient source exactness does not certify thermal accuracy.
#[allow(clippy::too_many_arguments)]
pub fn advance_temperature_candidate<'m, O: DepositionObserver, L: TemperatureDampingLaw>(
    input: RigidDampingThermalInput<'_>,
    effort: HeldInternalEffort,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    target_s: f64,
    config: DrivenDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<TemperatureDampingThermalCandidate<'m, O>, O::Error, L::Error> {
    advance_with_motion(
        input,
        thermal,
        observer,
        law,
        temperature_cells,
        DrivenMotion {
            effort,
            controls: config.mechanical,
        },
        target_s,
        config.into(),
        cancelled,
    )
}

#[derive(Clone, Copy)]
struct DrivenMotion {
    effort: HeldInternalEffort,
    controls: ViscousConnectorControls,
}
impl MechanicalAdvance for DrivenMotion {
    fn admit<E, S>(&self, input: RigidDampingThermalInput<'_>) -> Result<(), E, S> {
        if input.states.len() != 2 || input.links.len() != 1 || input.links[0].bodies != [0, 1] {
            return Err(RigidDampingThermalError::Invalid(
                "driven thermal motion needs two bodies and one original [0, 1] connector",
            ));
        }
        Ok(())
    }
    fn step<E, S>(
        &self,
        input: RigidDampingThermalInput<'_>,
        states: &[RigidBodyState],
        links: &[ViscousNetworkLink],
        dt_s: f64,
        requested: &dyn Fn() -> bool,
    ) -> Result<MechanicalStep, E, S> {
        let step = step_driven_viscous_connector_with_loads(
            [states[0], states[1]],
            [input.properties[0], input.properties[1]],
            &ViscousMotionLaw {
                arms_body_m: links[0].arms_body_m,
                damping: links[0].damping.clone(),
            },
            self.effort,
            input.loads,
            dt_s,
            self.controls,
            requested,
        )?;
        Ok(MechanicalStep {
            states_after: step.states_after.to_vec(),
            connector_dissipation_j: vec![step.dissipated_energy_j],
            connector_actuator_work_j: Some(vec![step.actuator_work_j]),
            dissipated_energy_j: step.dissipated_energy_j,
            applied_work_j: step.external_loads.applied_work_j,
            actuator_work_j: step.actuator_work_j,
            work_heat_residual_j: step.work_heat_residual_j,
            native_drift_energy_defect_j: step.native_drift_energy_defect_j,
            pair_relaxations: 2,
        })
    }
}
