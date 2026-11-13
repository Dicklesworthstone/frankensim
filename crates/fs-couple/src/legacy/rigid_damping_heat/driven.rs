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
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionReport, ImplicitConductionWorkspace,
};

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

/// Original driven-pair controls with native BE/SDIRK thermal transport.
pub type DrivenDampingImplicitConfig = RigidDampingImplicitConfig<ViscousConnectorControls>;

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
    .map(|(candidate, _, _)| candidate)
}

/// Advance the original driven pair with native implicit thermal conduction.
///
/// Signed actuator supply/recovery remains work; only the native positive
/// viscous loss enters finite heat storage. Each attempted thermal interval
/// prepares the same simultaneous mechanical trajectory as [`advance_candidate`].
/// Both SDIRK stages hold that one interval's mean heat, and thermal retries
/// recompute from the accepted bodies without publishing a partial candidate.
///
/// The graph, source restrictions and event admission follow
/// [`super::advance_implicit_candidate`]. Removing the conductive CFL does not
/// resolve moving geometry or held-source accuracy; refine maximum_dt_s
/// independently. The returned report is the original thermal solver's report.
#[allow(clippy::too_many_arguments)]
pub fn advance_implicit_candidate<'m, O: DepositionObserver>(
    input: RigidDampingThermalInput<'_>,
    effort: HeldInternalEffort,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: DrivenDampingImplicitConfig,
    workspace: &mut ImplicitConductionWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<(RigidDampingThermalCandidate<'m, O>, ImplicitConductionReport), O::Error> {
    let (candidate, _, _, report) = advance_with_source_driver(
        input,
        thermal,
        observer,
        FrozenSource(input.links),
        DrivenMotion {
            effort,
            controls: config.mechanical,
        },
        target_s,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
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

/// Driven-pair temperature feedback with native implicit thermal transport.
///
/// Original reference laws, explicit caloric bindings, endpoint source-domain
/// admission and signed-work accounting are those of
/// [`advance_temperature_candidate`]. Thermal retries retain accepted physical
/// history and refresh only the proposed interval. The candidate and complete
/// native thermal report are returned together after every gate accepts.
///
/// BE/SDIRK conduction retains first-order lagged temperature feedback and held
/// whole-interval mechanical heat. SDIRK thermal order does not certify the
/// feedback or moving-geometry split; independently refine maximum_dt_s.
#[allow(clippy::too_many_arguments)]
pub fn advance_temperature_implicit_candidate<
    'm,
    O: DepositionObserver,
    L: TemperatureDampingLaw,
>(
    input: RigidDampingThermalInput<'_>,
    effort: HeldInternalEffort,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: DrivenDampingImplicitConfig,
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
        DrivenMotion {
            effort,
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
            relaxation_memory: None,
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
