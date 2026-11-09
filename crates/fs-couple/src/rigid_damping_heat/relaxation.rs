//! Original Maxwell memory, native rigid motion and finite caloric storage.
//!
//! The fixed spatial source evolves every retained elastic memory coordinate
//! with the shared body bank. Each native moving step supplies independently
//! integrated Kelvin and Maxwell heat, signed actuator work and endpoint memory.
//! The existing thermal transaction commits those histories only after its
//! matching caloric interval accepts. No equilibrium spring or temperature
//! dependence of the retained memory law is inferred.

/// Calibrated WLF/Arrhenius kinetics with retained history and finite-cell feedback.
pub mod temperature;

use super::{
    DampingSource, FrozenSource, MechanicalAdvance, MechanicalStep, Result,
    RigidDampingImplicitConfig, RigidDampingThermalCandidate, RigidDampingThermalError,
    RigidDampingThermalInput, RigidDampingThermalReport, ThermalControls, ThermalDriver,
    advance_with_source_driver, require_implicit_report,
};
use fs_mbd::{
    RigidBodyState,
    spatial_loads::{
        internal_impulse::viscous::network::{
            ViscousNetworkLink,
            relaxation::{
                RelaxationNetworkError, RelaxationNetworkLink, relaxation_memory_energy_change,
            },
        },
        internal_motion::viscous::network::relaxation::{
            RelaxationNetworkMotionError, step_relaxation_network_with_loads,
            step_relaxation_network_with_loads_and_kinetics,
        },
    },
};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionReport, ImplicitConductionWorkspace,
};
use std::borrow::Cow;

pub use fs_mbd::spatial_loads::{
    internal_impulse::viscous::network::relaxation::{
        HeldRelaxationKinetics, RelaxationNetworkLimits, RelaxationNetworkMemory,
        SpatialRelaxationBranch,
    },
    internal_motion::viscous::{
        driven::HeldInternalEffort, network::relaxation::RelaxationNetworkMotionControls,
    },
};

/// Original rigid/viscous graph plus explicit elastic history and held effort.
///
/// Every original connector has one branch list, one memory list and one effort.
/// A memory entry is the native six-component energy amplitude in the original
/// branch's material coordinates [sqrt(J)], with storage |q|²/2. Zero denotes
/// relaxed memory; the current pose does not determine a preloaded history.
/// The original branch order, factors and relaxation times stay fixed.
#[derive(Clone, Copy, Debug)]
pub struct RelaxationDampingThermalInput<'a> {
    pub mechanical: RigidDampingThermalInput<'a>,
    pub branches: &'a [Vec<SpatialRelaxationBranch>],
    pub memory: &'a [Vec<[f64; 6]>],
    pub efforts: &'a [HeldInternalEffort],
}

/// Complete native graph preparation and finite thermal timeline limits.
#[derive(Clone, Copy, Debug)]
pub struct RelaxationDampingThermalConfig {
    pub mechanical: RelaxationNetworkMotionControls,
    /// Bounds native moving-geometry steps as well as caloric transport.
    pub thermal: StepConfig,
    pub maximum_heat_bindings: usize,
    pub energy_tolerance_j: f64,
}

/// Original Maxwell graph controls with native BE/SDIRK thermal transport.
pub type RelaxationDampingImplicitConfig =
    RigidDampingImplicitConfig<RelaxationNetworkMotionControls>;

impl From<RelaxationDampingThermalConfig> for ThermalControls {
    fn from(config: RelaxationDampingThermalConfig) -> Self {
        Self {
            thermal: config.thermal,
            maximum_heat_bindings: config.maximum_heat_bindings,
            energy_tolerance_j: config.energy_tolerance_j,
        }
    }
}
impl<E, S> From<RelaxationNetworkError> for RigidDampingThermalError<E, S> {
    fn from(error: RelaxationNetworkError) -> Self {
        RelaxationNetworkMotionError::from(error).into()
    }
}
impl<E, S> From<RelaxationNetworkMotionError> for RigidDampingThermalError<E, S> {
    fn from(error: RelaxationNetworkMotionError) -> Self {
        match error {
            RelaxationNetworkMotionError::Cancelled => Self::Cancelled,
            error => Self::RelaxationMotion(error),
        }
    }
}

/// Jointly accepted native mechanics, complete Maxwell memory and finite heat.
#[must_use]
pub struct RelaxationDampingThermalCandidate<'m, O> {
    native: RigidDampingThermalCandidate<'m, O>,
    memory: RelaxationNetworkMemory,
}
impl<'m, O> RelaxationDampingThermalCandidate<'m, O> {
    pub fn states(&self) -> &[RigidBodyState] {
        self.native.states()
    }
    pub fn memory(&self) -> &[Vec<[f64; 6]>] {
        &self.memory
    }
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        self.native.thermal()
    }
    pub fn observer(&self) -> &O {
        self.native.observer()
    }
    pub fn report(&self) -> &RigidDampingThermalReport {
        self.native.report()
    }
    pub fn into_parts(self) -> (RigidDampingThermalCandidate<'m, O>, RelaxationNetworkMemory) {
        (self.native, self.memory)
    }
}

/// Advance original body, material-memory and caloric histories together.
///
/// Heat routes refer to original connectors and receive each connector's Kelvin
/// plus Maxwell loss once. Elastic recovery can increase body kinetic energy
/// while still heating its receiver; only independently integrated loss is heat.
/// The final energy gate includes the original-to-final elastic memory change.
///
/// Branch factors/times, instantaneous viscosity and held loads remain fixed.
/// An application's calibrated frequency/strain domain and any thermal change
/// of the memory law remain the caller's responsibility; this function performs
/// no parameter rebind. Native geometry splitting and thermal transport still
/// need timestep refinement independently of the conservation tolerance.
pub fn advance_candidate<'m, O: DepositionObserver>(
    input: RelaxationDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    config: RelaxationDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationDampingThermalCandidate<'m, O>, O::Error> {
    advance_with_driver(
        input,
        thermal,
        observer,
        target_s,
        config.mechanical,
        config.into(),
        ThermalDriver::Explicit,
        cancelled,
    )
    .map(|(candidate, _)| candidate)
}

/// Advance original Maxwell history with native implicit thermal conduction.
///
/// Each attempted physical interval starts from the accepted body and complete
/// memory banks. Its independently integrated Kelvin/Maxwell loss feeds both
/// SDIRK stages as one whole-interval mean; elastic storage and signed actuator
/// work stay in their original accounts. Thermal acceptance commits bodies,
/// memory and caloric storage together. A retry replaces the pending proposal,
/// and any refusal returns no partial candidate.
///
/// Native graph, source and event admission follow
/// [`super::advance_implicit_candidate`]. Fixed branch laws and the independent
/// temporal-refinement obligations of [`advance_candidate`] remain in force:
/// SDIRK thermal order does not establish the moving-mechanics split's order.
#[allow(clippy::too_many_arguments)]
pub fn advance_implicit_candidate<'m, O: DepositionObserver>(
    input: RelaxationDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: RelaxationDampingImplicitConfig,
    workspace: &mut ImplicitConductionWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<(RelaxationDampingThermalCandidate<'m, O>, ImplicitConductionReport), O::Error> {
    let (candidate, report) = advance_with_driver(
        input,
        thermal,
        observer,
        target_s,
        config.mechanical,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
}

#[allow(clippy::too_many_arguments)]
fn advance_with_driver<'m, O: DepositionObserver>(
    input: RelaxationDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    mechanical: RelaxationNetworkMotionControls,
    config: ThermalControls,
    driver: ThermalDriver<'_, '_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<
    (RelaxationDampingThermalCandidate<'m, O>, Option<ImplicitConductionReport>),
    O::Error,
> {
    let links = prepare_links(input, mechanical, &mut cancelled)?;
    let mechanics = RelaxationMotion {
        links: &links,
        initial_memory: input.memory,
        memory: Cow::Borrowed(input.memory),
        efforts: input.efforts,
        controls: mechanical,
    };
    let (native, _, mechanics, report) = advance_with_source_driver(
        input.mechanical,
        thermal,
        observer,
        FrozenSource(input.mechanical.links),
        mechanics,
        target_s,
        config,
        driver,
        cancelled,
    )?;
    Ok((
        RelaxationDampingThermalCandidate {
            native,
            memory: mechanics.memory.into_owned(),
        },
        report,
    ))
}

fn prepare_links<E, S>(
    input: RelaxationDampingThermalInput<'_>,
    mechanical: RelaxationNetworkMotionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<RelaxationNetworkLink>, E, S> {
    if cancelled() {
        return Err(RigidDampingThermalError::Cancelled);
    }
    let n = input.mechanical.links.len();
    if input.branches.len() != n
        || input.memory.len() != n
        || input.efforts.len() != n
        || n > mechanical.limits.source.maximum_connectors
        || input.mechanical.states.len() > mechanical.limits.source.maximum_bodies
    {
        return Err(RigidDampingThermalError::Invalid(
            "complete original graph, branch, history and effort banks within native limits required",
        ));
    }
    let mut branches = 0_usize;
    for (law, history) in input.branches.iter().zip(input.memory) {
        if cancelled() {
            return Err(RigidDampingThermalError::Cancelled);
        }
        branches = branches
            .checked_add(law.len())
            .ok_or(RigidDampingThermalError::Invalid(
                "complete Maxwell branch count overflows",
            ))?;
        if law.len() != history.len() || branches > mechanical.limits.maximum_branches {
            return Err(RigidDampingThermalError::Invalid(
                "every original branch needs one complete history within native limits",
            ));
        }
    }
    let mut links = Vec::with_capacity(n);
    for (viscous, branches) in input.mechanical.links.iter().zip(input.branches) {
        if cancelled() {
            return Err(RigidDampingThermalError::Cancelled);
        }
        links.push(RelaxationNetworkLink {
            viscous: viscous.clone(),
            branches: branches.clone(),
        });
    }
    Ok(links)
}

#[derive(Clone)]
struct RelaxationMotion<'a> {
    links: &'a [RelaxationNetworkLink],
    initial_memory: &'a [Vec<[f64; 6]>],
    memory: Cow<'a, [Vec<[f64; 6]>]>,
    efforts: &'a [HeldInternalEffort],
    controls: RelaxationNetworkMotionControls,
}
impl MechanicalAdvance for RelaxationMotion<'_> {
    fn admit<E, S>(&self, _input: RigidDampingThermalInput<'_>) -> Result<(), E, S> {
        // Validate the complete initial energy coordinates before thermal events.
        relaxation_memory_energy_change(self.initial_memory, self.initial_memory)?;
        Ok(())
    }
    fn step<E, S>(
        &self,
        input: RigidDampingThermalInput<'_>,
        states: &[RigidBodyState],
        _links: &[ViscousNetworkLink],
        dt_s: f64,
        requested: &dyn Fn() -> bool,
    ) -> Result<MechanicalStep, E, S> {
        self.step_native(input, states, None, dt_s, requested)
    }
    fn step_with_source<E, S: DampingSource>(
        &self,
        input: RigidDampingThermalInput<'_>,
        states: &[RigidBodyState],
        source: &S,
        dt_s: f64,
        requested: &dyn Fn() -> bool,
    ) -> Result<MechanicalStep, E, S::Error> {
        self.step_native(input, states, source.relaxation_kinetics(), dt_s, requested)
    }
    fn accept<E, S>(&mut self, step: &mut MechanicalStep) -> Result<(), E, S> {
        self.memory = Cow::Owned(step.relaxation_memory.take().ok_or(
            RigidDampingThermalError::Invalid("native Maxwell step lost its complete memory bank"),
        )?);
        Ok(())
    }
    fn memory_energy_change<E, S>(&self) -> Result<f64, E, S> {
        Ok(relaxation_memory_energy_change(
            self.initial_memory,
            &self.memory,
        )?)
    }
}

impl RelaxationMotion<'_> {
    fn step_native<E, S>(
        &self,
        input: RigidDampingThermalInput<'_>,
        states: &[RigidBodyState],
        kinetics: Option<HeldRelaxationKinetics<'_>>,
        dt_s: f64,
        requested: &dyn Fn() -> bool,
    ) -> Result<MechanicalStep, E, S> {
        let mut native = match kinetics {
            Some(kinetics) => step_relaxation_network_with_loads_and_kinetics(
                states,
                input.properties,
                self.links,
                &self.memory,
                self.efforts,
                kinetics,
                input.loads,
                dt_s,
                self.controls,
                requested,
            )?,
            None => step_relaxation_network_with_loads(
                states,
                input.properties,
                self.links,
                &self.memory,
                self.efforts,
                input.loads,
                dt_s,
                self.controls,
                requested,
            )?,
        };
        let connector_dissipation_j = core::mem::take(&mut native.connector_dissipation_j);
        let connector_actuator_work_j =
            Some(core::mem::take(&mut native.connector_actuator_work_j));
        let dissipated_energy_j = native.dissipated_energy_j;
        let applied_work_j = native.external_loads.applied_work_j;
        let actuator_work_j = native.actuator_work_j;
        let work_heat_residual_j = native.work_heat_residual_j;
        let native_drift_energy_defect_j = native.native_drift_energy_defect_j;
        let (states_after, memory) = native.into_state_and_memory();
        Ok(MechanicalStep {
            states_after,
            relaxation_memory: Some(memory),
            connector_dissipation_j,
            connector_actuator_work_j,
            dissipated_energy_j,
            applied_work_j,
            actuator_work_j,
            work_heat_residual_j,
            native_drift_energy_defect_j,
            pair_relaxations: 0,
        })
    }
}
