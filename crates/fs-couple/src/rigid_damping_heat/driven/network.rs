//! Simultaneous driven damping graphs on the original finite thermal timeline.
//!
//! Every connector acts on the shared body bank during each native source
//! stage. Cycles and parallel links retain independent actuator work and
//! positive heat. No pairwise ordering or source-substep approximation is used.
//! Moving geometry and lagged temperature still need temporal refinement.

use super::super::{
    FrozenSource, MechanicalAdvance, MechanicalStep, Result, RigidDampingThermalCandidate,
    RigidDampingImplicitConfig, RigidDampingThermalError, RigidDampingThermalInput, ThermalControls,
    advance_with_source, advance_with_source_driver, require_implicit_report,
    temperature::{
        TemperatureDampingLaw, TemperatureDampingThermalCandidate, advance_with_motion,
        advance_with_motion_driver,
    },
};
use fs_mbd::{
    RigidBodyState,
    spatial_loads::{
        internal_impulse::viscous::network::ViscousNetworkLink,
        internal_motion::viscous::network::driven::{
            DrivenViscousNetworkMotionControls, DrivenViscousNetworkMotionError,
            step_driven_viscous_network_with_loads,
        },
    },
};
use fs_thermal::enthalpy::{CellId, DepositionObserver, DepositionSimulation, StepConfig};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionReport, ImplicitConductionWorkspace,
};

pub use super::HeldInternalEffort;

/// Complete simultaneous source preparation and native thermal interval limits.
///
/// The native mechanical limits bound the complete dense graph preparation,
/// including both moving-geometry endpoints. They are independent of the sparse
/// passive network's pair-sweep limits. The thermal owner bounds actual event
/// intervals, transport and caloric storage.
#[derive(Clone, Copy, Debug)]
pub struct DrivenNetworkDampingThermalConfig {
    pub mechanical: DrivenViscousNetworkMotionControls,
    /// Also bounds moving geometry and first-order temperature feedback.
    pub thermal: StepConfig,
    pub maximum_heat_bindings: usize,
    pub energy_tolerance_j: f64,
}

/// Original simultaneous graph controls with native BE/SDIRK thermal transport.
pub type DrivenNetworkDampingImplicitConfig =
    RigidDampingImplicitConfig<DrivenViscousNetworkMotionControls>;

impl From<DrivenNetworkDampingThermalConfig> for ThermalControls {
    fn from(config: DrivenNetworkDampingThermalConfig) -> Self {
        Self {
            thermal: config.thermal,
            maximum_heat_bindings: config.maximum_heat_bindings,
            energy_tolerance_j: config.energy_tolerance_j,
        }
    }
}
impl<E, S> From<DrivenViscousNetworkMotionError> for RigidDampingThermalError<E, S> {
    fn from(error: DrivenViscousNetworkMotionError) -> Self {
        match error {
            DrivenViscousNetworkMotionError::Cancelled => Self::Cancelled,
            error => Self::DrivenNetworkMotion(error),
        }
    }
}

/// Advance every original connector's held effort and damping simultaneously.
///
/// There must be one explicit effort for every original link, including zero
/// effort for unactuated links. Each effort is applied at attachment B in its
/// declared world or first-material frame; the native owner transports the
/// reciprocal wrench and refreshes geometry at both source stages.
///
/// Original connector order also identifies heat routes and signed actuator
/// work. Recovery at one actuator cannot erase another connector's positive
/// heat or be deposited as thermal energy. Shared receivers and unallocated
/// outflow use the same native thermal transaction as the passive and pair
/// entries. All physical candidates are discarded on refusal or cancellation.
pub fn advance_candidate<'m, O: DepositionObserver>(
    input: RigidDampingThermalInput<'_>,
    efforts: &[HeldInternalEffort],
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    config: DrivenNetworkDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<RigidDampingThermalCandidate<'m, O>, O::Error> {
    advance_with_source(
        input,
        thermal,
        observer,
        FrozenSource(input.links),
        DrivenNetworkMotion {
            efforts,
            controls: config.mechanical,
        },
        target_s,
        config.into(),
        cancelled,
    )
    .map(|(candidate, _, _)| candidate)
}

/// Advance a simultaneous driven graph with native implicit thermal transport.
///
/// The original whole-graph mechanical owner supplies each connector's signed
/// actuator work and positive heat. One prepared whole-interval heat proposal
/// feeds both SDIRK stages; only accepted thermal intervals advance bodies and
/// histories. Nonlinear retries preserve the accepted mechanical bank.
///
/// Native graph, source and event admission follow
/// [`super::super::advance_implicit_candidate`]. The complete thermal report is
/// returned alongside the joint candidate. Refine maximum_dt_s for moving
/// geometry and held-source accuracy independently of conductive stability.
#[allow(clippy::too_many_arguments)]
pub fn advance_implicit_candidate<'m, O: DepositionObserver>(
    input: RigidDampingThermalInput<'_>,
    efforts: &[HeldInternalEffort],
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: DrivenNetworkDampingImplicitConfig,
    workspace: &mut ImplicitConductionWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<(RigidDampingThermalCandidate<'m, O>, ImplicitConductionReport), O::Error> {
    let (candidate, _, _, report) = advance_with_source_driver(
        input,
        thermal,
        observer,
        FrozenSource(input.links),
        DrivenNetworkMotion {
            efforts,
            controls: config.mechanical,
        },
        target_s,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
}

/// Simultaneous driven graph with original finite-cell temperature feedback.
///
/// Each connector's reference damping form is selected from its explicit
/// temperature coordinate on each actual thermal interval. The original source
/// law and caloric endpoint must admit the accepted result before publication.
/// Reuse the original reference links on continuation, not the shifted endpoint
/// forms. Held effort acquires no invented temperature dependence.
///
/// Feedback is first-order and explicit. Exact held-coefficient graph evolution
/// establishes neither resolved thermal accuracy nor moving-geometry accuracy.
#[allow(clippy::too_many_arguments)]
pub fn advance_temperature_candidate<'m, O: DepositionObserver, L: TemperatureDampingLaw>(
    input: RigidDampingThermalInput<'_>,
    efforts: &[HeldInternalEffort],
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    target_s: f64,
    config: DrivenNetworkDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<TemperatureDampingThermalCandidate<'m, O>, O::Error, L::Error> {
    advance_with_motion(
        input,
        thermal,
        observer,
        law,
        temperature_cells,
        DrivenNetworkMotion {
            efforts,
            controls: config.mechanical,
        },
        target_s,
        config.into(),
        cancelled,
    )
}

/// Simultaneous driven-graph feedback with native implicit thermal conduction.
///
/// Every original effort, reference damping form and explicit temperature
/// binding is retained as in [`advance_temperature_candidate`]. Endpoint law
/// admission follows all native thermal events; refusal publishes no physical
/// candidate. Source preparation and candidate construction share the passive
/// temperature path, including accepted-history and retry semantics.
///
/// Temperature feedback remains first-order and lagged. Both SDIRK stages hold
/// one complete mechanical interval's mean heat; the thermal scheme does not
/// raise the feedback or moving-geometry order. Refine maximum_dt_s separately.
#[allow(clippy::too_many_arguments)]
pub fn advance_temperature_implicit_candidate<
    'm,
    O: DepositionObserver,
    L: TemperatureDampingLaw,
>(
    input: RigidDampingThermalInput<'_>,
    efforts: &[HeldInternalEffort],
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: DrivenNetworkDampingImplicitConfig,
    workspace: &mut ImplicitConductionWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (TemperatureDampingThermalCandidate<'m, O>, ImplicitConductionReport),
    O::Error,
    L::Error,
> {
    let (candidate, report) = advance_with_motion_driver(
        input,
        thermal,
        observer,
        law,
        temperature_cells,
        DrivenNetworkMotion {
            efforts,
            controls: config.mechanical,
        },
        target_s,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
}

#[derive(Clone, Copy)]
struct DrivenNetworkMotion<'a> {
    efforts: &'a [HeldInternalEffort],
    controls: DrivenViscousNetworkMotionControls,
}
impl MechanicalAdvance for DrivenNetworkMotion<'_> {
    fn admit<E, S>(&self, input: RigidDampingThermalInput<'_>) -> Result<(), E, S> {
        if self.efforts.len() != input.links.len()
            || input.states.len() > self.controls.limits.maximum_bodies
            || input.links.len() > self.controls.limits.maximum_connectors
        {
            return Err(RigidDampingThermalError::Invalid(
                "one effort per original link and complete body/link banks within native graph limits required",
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
        let mut step = step_driven_viscous_network_with_loads(
            states,
            input.properties,
            links,
            self.efforts,
            input.loads,
            dt_s,
            self.controls,
            requested,
        )?;
        Ok(MechanicalStep {
            relaxation_memory: None,
            connector_dissipation_j: core::mem::take(&mut step.connector_dissipation_j),
            connector_actuator_work_j: Some(core::mem::take(&mut step.connector_actuator_work_j)),
            dissipated_energy_j: step.dissipated_energy_j,
            applied_work_j: step.external_loads.applied_work_j,
            actuator_work_j: step.actuator_work_j,
            work_heat_residual_j: step.work_heat_residual_j,
            native_drift_energy_defect_j: step.native_drift_energy_defect_j,
            // This native graph uses simultaneous source preparations, not pair sweeps.
            pair_relaxations: 0,
            states_after: step.into_states_after(),
        })
    }
}
