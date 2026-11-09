//! Moving spatial damper networks on the original finite thermal timeline.
//!
//! The native mechanical owner supplies each connector's integrated physical
//! loss. Explicit cell fractions route that loss as interval-average watts on
//! the SAME interval selected by the native thermal transport/event owner.
//! Original caloric laws convert energy to temperature; numerical mechanical
//! defects never heat a receiver. Unrouted loss remains an explicit outflow.
//!
//! The baseline entry holds damping forms fixed; [`temperature`] feeds original
//! caloric temperatures back into explicitly authored passive damping laws.
//! Inertias, material arms and external loads remain held. The caller owns their
//! spatial correspondence and material validity. Thermal geometry remains
//! caller-owned; no reversible thermoelastic work or heat partition is inferred.
//! Native transport order and moving-geometry splitting require temporal and
//! spatial refinement independently of conservation tolerances.

/// Temperature-dependent passive laws on the same accepted native timeline.
pub mod temperature;

/// Simultaneous internal actuation, damping heat and finite caloric feedback.
pub mod driven;

/// Original Maxwell memory and native rigid motion on the finite thermal timeline.
pub mod relaxation;

use fs_math::dd::Dd;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::ViscousNetworkLink;
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads,
    viscous::network::{
        ViscousNetworkMotionControls, ViscousNetworkMotionError, ViscousNetworkMotionStep,
        step_viscous_network_with_loads,
    },
};
use fs_mbd::{MassProperties, RigidBodyState, RigidKineticEnergyChange};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionError, ImplicitConductionGraph,
    ImplicitConductionObservedError, ImplicitConductionReport, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    HeatTransfer, StepConfig, collect_source_rates,
};
use std::{borrow::Cow, cell::RefCell, collections::BTreeMap, convert::Infallible};

/// Borrowed original mechanical state, laws and explicitly authored heat routes.
#[derive(Clone, Copy, Debug)]
pub struct RigidDampingThermalInput<'a> {
    /// Accepted mechanical time; it must equal the original thermal clock.
    pub time_s: f64,
    pub states: &'a [RigidBodyState],
    pub properties: &'a [MassProperties],
    pub links: &'a [ViscousNetworkLink],
    /// Held for the complete call, including every thermal event/substep.
    pub loads: HeldExternalLoads<'a>,
    /// One map per ORIGINAL connector, in its original order. Each map assigns
    /// nonnegative fractions totaling at most one to unique physical CellIds.
    /// Different connectors may share a receiver. The remainder is unallocated
    /// loss, not an invented bath. Receiver IDs must be active on each actual
    /// interval; same-time native birth events run before that check.
    pub heat_fractions: &'a [BTreeMap<CellId, f64>],
}

/// Original integration/resource controls and one whole-call energy allowance.
#[derive(Clone, Copy, Debug)]
pub struct RigidDampingThermalConfig {
    /// Per-native-motion-step limits and independent physical admission gates.
    pub mechanical: ViscousNetworkMotionControls,
    /// Native CFL/event/substep controls. max_dt_s also bounds mechanical steps.
    pub thermal: StepConfig,
    /// Total authored connector-to-cell entries, including zero fractions.
    pub maximum_heat_bindings: usize,
    /// Bounds complete residuals, accumulated absolute interval mechanical
    /// defects and absolute cellwise heat-handoff defects independently [J].
    pub energy_tolerance_j: f64,
}

/// Native thermal quadrature for the same accepted mechanical intervals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImplicitThermalScheme {
    BackwardEuler,
    Sdirk2,
}

/// Original mechanical limits and native implicit thermal controls.
#[derive(Clone, Copy, Debug)]
pub struct RigidDampingImplicitConfig<M = ViscousNetworkMotionControls> {
    /// Per-attempt controls for the original mechanical owner. Successfully
    /// prepared rejected attempts retain their reported native work counts.
    pub mechanical: M,
    /// maximum_dt_s also bounds the physical mechanical interval. Source
    /// restrictions and exact native events may shorten it further.
    pub thermal: ImplicitConductionControls,
    pub scheme: ImplicitThermalScheme,
    pub maximum_heat_bindings: usize,
    /// Independent whole-call mechanical, handoff and combined energy gates.
    pub energy_tolerance_j: f64,
}
impl<M> RigidDampingImplicitConfig<M> {
    fn thermal_controls(&self) -> ThermalControls {
        ThermalControls {
            thermal: StepConfig {
                max_dt_s: self.thermal.maximum_dt_s,
                energy_tolerance_j: self.thermal.total_tolerance_j,
            },
            maximum_heat_bindings: self.maximum_heat_bindings,
            energy_tolerance_j: self.energy_tolerance_j,
        }
    }

    fn driver<'a, 'law>(
        &self,
        conduction: &'a ImplicitConductionGraph<'law>,
        workspace: &'a mut ImplicitConductionWorkspace,
    ) -> ThermalDriver<'a, 'law> {
        ThermalDriver::Implicit {
            conduction,
            controls: self.thermal,
            scheme: self.scheme,
            workspace,
        }
    }
}

/// Physical heat, external work and numerical defects retain separate accounts.
#[derive(Clone, Debug)]
pub struct RigidDampingThermalReport {
    pub time_s: f64,
    pub thermal: EnergyAudit,
    /// Original edge order, summed over every accepted interval [J].
    pub connector_heat_j: Vec<f64>,
    /// Signed supplied/recovered work for each original connector [J].
    /// Passive links retain exact zero; heat remains independently integrated.
    pub connector_actuator_work_j: Vec<f64>,
    pub dissipated_energy_j: f64,
    /// Represented generated source energy; included once in the raw thermal audit.
    pub deposited_heat_j: f64,
    pub unallocated_energy_j: f64,
    /// Raw thermal heat less the represented damper source [J].
    pub external_thermal_heat_j: f64,
    /// Applied force/couple work; gravity is in mechanical potential storage.
    pub mechanical_work_j: f64,
    /// Signed internal actuator supply/recovery, independent of passive heat [J].
    pub actuator_work_j: f64,
    /// Elastic Maxwell storage change; exact zero for memoryless entries [J].
    pub memory_energy_change_j: f64,
    /// Kinetic, gravitational and retained elastic-memory storage change [J].
    pub mechanical_energy_change_j: f64,
    pub mechanical_residual_j: f64,
    pub native_drift_energy_defect_j: f64,
    pub absolute_mechanical_defect_j: f64,
    /// Deposited + unallocated - independently integrated mechanical heat [J].
    pub heat_handoff_residual_j: f64,
    pub absolute_heat_handoff_defect_j: f64,
    /// Delta mechanical + delta original thermal storage + unallocated loss
    /// - original incoming inventory - outside mechanical/thermal supply [J].
    pub combined_energy_residual_j: f64,
    /// Actual native work, including successfully prepared but rejected intervals.
    pub pair_relaxations: usize,
}

/// No independently publishable mechanical or thermal half-step on refusal.
#[must_use]
pub struct RigidDampingThermalCandidate<'m, O> {
    states: Vec<RigidBodyState>,
    thermal: DepositionSimulation<'m>,
    observer: O,
    report: RigidDampingThermalReport,
}
impl<'m, O> RigidDampingThermalCandidate<'m, O> {
    pub fn states(&self) -> &[RigidBodyState] {
        &self.states
    }
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub fn observer(&self) -> &O {
        &self.observer
    }
    pub fn report(&self) -> &RigidDampingThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        Vec<RigidBodyState>,
        DepositionSimulation<'m>,
        O,
        RigidDampingThermalReport,
    ) {
        (self.states, self.thermal, self.observer, self.report)
    }
}

#[derive(Debug)]
pub enum RigidDampingThermalError<E, S = Infallible> {
    Invalid(&'static str),
    Motion(ViscousNetworkMotionError),
    DrivenMotion(fs_mbd::spatial_loads::internal_motion::viscous::driven::DrivenViscousMotionError),
    DrivenNetworkMotion(
        fs_mbd::spatial_loads::internal_motion::viscous::network::driven::DrivenViscousNetworkMotionError,
    ),
    RelaxationMotion(
        fs_mbd::spatial_loads::internal_motion::viscous::network::relaxation::RelaxationNetworkMotionError,
    ),
    Thermal(DepositionError),
    ImplicitThermal(ImplicitConductionError),
    Observer(E),
    /// Original connector index and the authored source's own refusal.
    Source {
        connector: usize,
        error: S,
    },
    Energy {
        residual_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}
impl<E, S> From<DepositionError> for RigidDampingThermalError<E, S> {
    fn from(error: DepositionError) -> Self {
        match error {
            DepositionError::Cancelled => Self::Cancelled,
            error => Self::Thermal(error),
        }
    }
}
impl<E, S> From<ViscousNetworkMotionError> for RigidDampingThermalError<E, S> {
    fn from(error: ViscousNetworkMotionError) -> Self {
        match error {
            ViscousNetworkMotionError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}
impl<E, S> From<ImplicitConductionError> for RigidDampingThermalError<E, S> {
    fn from(error: ImplicitConductionError) -> Self {
        match error {
            ImplicitConductionError::Thermal(DepositionError::Cancelled) => Self::Cancelled,
            error => Self::ImplicitThermal(error),
        }
    }
}
impl<E: core::fmt::Debug, S: core::fmt::Debug> core::fmt::Display
    for RigidDampingThermalError<E, S>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid damping/thermal advance refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static, S: std::error::Error + 'static> std::error::Error
    for RigidDampingThermalError<E, S>
{
}
type Result<T, E, S = Infallible> = core::result::Result<T, RigidDampingThermalError<E, S>>;

fn number<E, S>(value: f64) -> Result<f64, E, S> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn value<E, S>(wide: Dd) -> Result<f64, E, S> {
    number(wide.hi)?;
    number(wide.lo)?;
    number(wide.to_f64())
}
fn balance<E, S>(residual: f64, allowance: f64) -> Result<(), E, S> {
    if number(residual)?.abs() > allowance {
        Err(RigidDampingThermalError::Energy {
            residual_j: residual,
            allowance_j: allowance,
        })
    } else {
        Ok(())
    }
}

/// Advance the original body bank, thermal timeline and observer transactionally.
///
/// Every actual thermal interval gets exactly one native moving-network call.
/// Its receipt stays pending until the matching thermal step accepts. Native
/// event clipping, sources, transport, boundary hooks and observer transitions
/// remain active. A rejected heat update, observer callback, energy check or
/// cancellation discards the complete thermal/observer candidate and its bodies.
/// Observer Clone must isolate mutable physical history, as for advance_observed.
///
/// Work is the native mechanical/thermal work plus O(bindings * log(receivers))
/// heat routing per interval. The borrowed initial bank is not cloned; accepted
/// endpoint banks are moved from native receipts. No trajectory is retained.
/// Native per-step allocations remain, plus one bounded receiver map and the
/// current pending receipt. Heat fractions do not infer physical correspondence,
/// damper temperature, moving thermal geometry or material validity.
pub fn advance_candidate<'m, O: DepositionObserver>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    config: RigidDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<RigidDampingThermalCandidate<'m, O>, O::Error> {
    advance_with_source(
        input,
        thermal,
        observer,
        FrozenSource(input.links),
        NetworkMotion(config.mechanical),
        target_s,
        config.into(),
        cancelled,
    )
    .map(|(candidate, _, _)| candidate)
}

/// Advance native rigid damping and stiff thermal transport on one shared clock.
///
/// The caller supplies the original prepared graph and its reusable workspace.
/// One native mechanical proposal supplies whole-interval physical loss; both
/// SDIRK stages use that same mean heat. Only the final accepted thermal step
/// publishes the pending body states and observer history. Nonlinear halving
/// recomputes mechanical loss from the unchanged accepted body bank, retaining
/// attempted work counts. No artificial conductive CFL limits the interval.
///
/// The graph owns all conduction; inner observer conductance overrides refuse
/// through its native admission. Original source timestep restrictions and
/// event/transition clipping remain active. Frozen damping forms, held loads,
/// caloric/geometry correspondence and temporal-refinement obligations remain
/// those of [`advance_candidate`]. SDIRK thermal order does not establish the
/// order of held observed forcing or the moving-mechanics split.
#[allow(clippy::too_many_arguments)]
pub fn advance_implicit_candidate<'m, O: DepositionObserver>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: RigidDampingImplicitConfig,
    workspace: &mut ImplicitConductionWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (
        RigidDampingThermalCandidate<'m, O>,
        ImplicitConductionReport,
    ),
    O::Error,
> {
    let (candidate, _, _, report) = advance_with_source_driver(
        input,
        thermal,
        observer,
        FrozenSource(input.links),
        NetworkMotion(config.mechanical),
        target_s,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
}

fn require_implicit_report<E, S>(
    report: Option<ImplicitConductionReport>,
) -> Result<ImplicitConductionReport, E, S> {
    report.ok_or(RigidDampingThermalError::Invalid(
        "implicit thermal advance omitted its native report",
    ))
}

// Only the source binding varies. The original timeline, native mechanics,
// per-cell heat handoff and independent final energy admission remain shared.
trait DampingSource: Clone {
    type Error;
    fn links(&self) -> &[ViscousNetworkLink];
    /// Explicit rate-only changes in the original Maxwell/Kelvin source basis.
    fn relaxation_kinetics(&self) -> Option<relaxation::HeldRelaxationKinetics<'_>> {
        None
    }
    fn synchronize<E>(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        requested: &dyn Fn() -> bool,
    ) -> Result<(), E, Self::Error>;
}

#[derive(Clone, Copy)]
struct FrozenSource<'a>(&'a [ViscousNetworkLink]);
impl DampingSource for FrozenSource<'_> {
    type Error = Infallible;
    fn links(&self) -> &[ViscousNetworkLink] {
        self.0
    }
    fn synchronize<E>(
        &mut self,
        _thermal: &DepositionSimulation<'_>,
        _requested: &dyn Fn() -> bool,
    ) -> Result<(), E> {
        Ok(())
    }
}

// A small private adapter shares the original thermal transaction among the
// native passive, driven-pair and simultaneous driven-graph solvers.
#[derive(Clone, Copy)]
struct ThermalControls {
    thermal: StepConfig,
    maximum_heat_bindings: usize,
    energy_tolerance_j: f64,
}
impl From<RigidDampingThermalConfig> for ThermalControls {
    fn from(config: RigidDampingThermalConfig) -> Self {
        Self {
            thermal: config.thermal,
            maximum_heat_bindings: config.maximum_heat_bindings,
            energy_tolerance_j: config.energy_tolerance_j,
        }
    }
}

enum ThermalDriver<'a, 'law> {
    Explicit,
    Implicit {
        conduction: &'a ImplicitConductionGraph<'law>,
        controls: ImplicitConductionControls,
        scheme: ImplicitThermalScheme,
        workspace: &'a mut ImplicitConductionWorkspace,
    },
}

#[derive(Clone)]
struct MechanicalStep {
    states_after: Vec<RigidBodyState>,
    // A native candidate stays pending until its matching caloric step accepts.
    relaxation_memory: Option<relaxation::RelaxationNetworkMemory>,
    connector_dissipation_j: Vec<f64>,
    // Passive motion has no actuator bank and needs no per-step zero allocation.
    connector_actuator_work_j: Option<Vec<f64>>,
    dissipated_energy_j: f64,
    applied_work_j: f64,
    actuator_work_j: f64,
    work_heat_residual_j: f64,
    native_drift_energy_defect_j: f64,
    pair_relaxations: usize,
}
impl From<ViscousNetworkMotionStep> for MechanicalStep {
    fn from(step: ViscousNetworkMotionStep) -> Self {
        let [_, endpoint] = step.half_sources;
        Self {
            states_after: endpoint.states_after,
            relaxation_memory: None,
            connector_actuator_work_j: None,
            connector_dissipation_j: step.connector_dissipation_j,
            dissipated_energy_j: step.dissipated_energy_j,
            applied_work_j: step.external_loads.applied_work_j,
            actuator_work_j: 0.0,
            work_heat_residual_j: step.work_heat_residual_j,
            native_drift_energy_defect_j: step.native_drift_energy_defect_j,
            pair_relaxations: step.pair_relaxations,
        }
    }
}
trait MechanicalAdvance: Clone {
    fn admit<E, S>(&self, input: RigidDampingThermalInput<'_>) -> Result<(), E, S>;
    fn step<E, S>(
        &self,
        input: RigidDampingThermalInput<'_>,
        states: &[RigidBodyState],
        links: &[ViscousNetworkLink],
        dt_s: f64,
        requested: &dyn Fn() -> bool,
    ) -> Result<MechanicalStep, E, S>;
    fn step_with_source<E, S: DampingSource>(
        &self,
        input: RigidDampingThermalInput<'_>,
        states: &[RigidBodyState],
        source: &S,
        dt_s: f64,
        requested: &dyn Fn() -> bool,
    ) -> Result<MechanicalStep, E, S::Error> {
        self.step(input, states, source.links(), dt_s, requested)
    }
    fn accept<E, S>(&mut self, step: &mut MechanicalStep) -> Result<(), E, S> {
        if step.relaxation_memory.is_some() {
            return Err(RigidDampingThermalError::Invalid(
                "memoryless motion received an unexpected elastic history",
            ));
        }
        Ok(())
    }
    /// Independent original-to-final storage change after accepted intervals.
    fn memory_energy_change<E, S>(&self) -> Result<f64, E, S> {
        Ok(0.0)
    }
}
#[derive(Clone, Copy)]
struct NetworkMotion(ViscousNetworkMotionControls);
impl MechanicalAdvance for NetworkMotion {
    fn admit<E, S>(&self, input: RigidDampingThermalInput<'_>) -> Result<(), E, S> {
        if input.states.len() > self.0.limits.maximum_bodies
            || input.links.len() > self.0.limits.maximum_connectors
        {
            return Err(RigidDampingThermalError::Invalid(
                "complete original body and connector banks must fit mechanical limits",
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
        Ok(step_viscous_network_with_loads(
            states,
            input.properties,
            links,
            input.loads,
            dt_s,
            self.0,
            requested,
        )?
        .into())
    }
}

fn advance_with_source<'m, O: DepositionObserver, S: DampingSource, M: MechanicalAdvance>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    source: S,
    mechanics: M,
    target_s: f64,
    config: ThermalControls,
    cancelled: impl FnMut() -> bool,
) -> Result<(RigidDampingThermalCandidate<'m, O>, S, M), O::Error, S::Error> {
    advance_with_source_driver(
        input,
        thermal,
        observer,
        source,
        mechanics,
        target_s,
        config,
        ThermalDriver::Explicit,
        cancelled,
    )
    .map(|(candidate, source, mechanics, _)| (candidate, source, mechanics))
}

#[allow(clippy::too_many_arguments)]
fn advance_with_source_driver<'m, O: DepositionObserver, S: DampingSource, M: MechanicalAdvance>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    source: S,
    mechanics: M,
    target_s: f64,
    config: ThermalControls,
    driver: ThermalDriver<'_, '_>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (
        RigidDampingThermalCandidate<'m, O>,
        S,
        M,
        Option<ImplicitConductionReport>,
    ),
    O::Error,
    S::Error,
> {
    let cancellation = RefCell::new(cancelled);
    let requested = || (*cancellation.borrow_mut())();
    let checkpoint = || {
        if requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    if input.time_s.to_bits() != thermal.time_s().to_bits()
        || !target_s.is_finite()
        || target_s <= input.time_s
        || input.states.len() != input.properties.len()
        || input.states.len() != input.loads.wrenches.len()
        || input.links.len() != input.heat_fractions.len()
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j < 0.0
        || !config.thermal.energy_tolerance_j.is_finite()
        || config.thermal.energy_tolerance_j < 0.0
        || config.thermal.energy_tolerance_j > config.energy_tolerance_j
    {
        return Err(RigidDampingThermalError::Invalid(
            "synchronized forward clocks, complete original banks/routes and finite energy limits required",
        ));
    }
    mechanics.admit::<O::Error, S::Error>(input)?;
    let mut binding_count = 0_usize;
    let mut remainders = Vec::with_capacity(input.links.len());
    for route in input.heat_fractions {
        checkpoint()?;
        binding_count = binding_count
            .checked_add(route.len())
            .ok_or(DepositionError::Budget("rigid damper heat bindings"))?;
        if binding_count > config.maximum_heat_bindings {
            return Err(DepositionError::Budget("rigid damper heat bindings").into());
        }
        let mut remaining = Dd::ONE;
        for &fraction in route.values() {
            checkpoint()?;
            if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                return Err(RigidDampingThermalError::Invalid(
                    "finite nonnegative heat fractions required",
                ));
            }
            remaining = remaining - Dd::from_f64(fraction);
            // Retain sub-ulp overspend: rounded sum == 1 is insufficient.
            if remaining.hi < 0.0 || (remaining.hi == 0.0 && remaining.lo < 0.0) {
                return Err(RigidDampingThermalError::Invalid(
                    "connector heat fractions exceed one",
                ));
            }
        }
        remainders.push(remaining);
    }
    let coupled = HeatObserver {
        input,
        source,
        mechanics,
        states: Cow::Borrowed(input.states),
        time_s: input.time_s,
        inner: Cow::Borrowed(observer),
        remainders: &remainders,
        config,
        requested: &requested,
        pending: None,
        accounts: Accounts::default(),
        connector_heat: vec![Dd::ZERO; input.links.len()],
        connector_actuator_work: vec![Dd::ZERO; input.links.len()],
        pair_relaxations: 0,
    };
    let (thermal, mut coupled, audit, implicit_report) = match driver {
        ThermalDriver::Explicit => {
            let (thermal, coupled, audit) =
                thermal.advance_observed(&coupled, target_s, config.thermal, checkpoint)?;
            (thermal, coupled, audit, None)
        }
        ThermalDriver::Implicit {
            conduction,
            controls,
            scheme,
            workspace,
        } => {
            let advance = match scheme {
                ImplicitThermalScheme::BackwardEuler => conduction.advance_observed_to(
                    thermal, &coupled, target_s, controls, workspace, checkpoint,
                ),
                ImplicitThermalScheme::Sdirk2 => conduction.advance_sdirk2_observed_to(
                    thermal, &coupled, target_s, controls, workspace, checkpoint,
                ),
            };
            let (thermal, coupled, report) = advance.map_err(|error| match error {
                ImplicitConductionObservedError::Implicit(error) => error.into(),
                ImplicitConductionObservedError::Observer(error) => error,
            })?;
            (thermal, coupled, report.energy, Some(report))
        }
    };
    if coupled.pending.is_some() || coupled.time_s.to_bits() != target_s.to_bits() {
        return Err(RigidDampingThermalError::Invalid(
            "mechanical and thermal endpoints must agree",
        ));
    }
    // Native events at target_s have all finished. Admit the source kinetics
    // at this exact endpoint before publishing either physical owner.
    coupled
        .source
        .synchronize::<O::Error>(&thermal, &requested)?;
    let mut kinetic = RigidKineticEnergyChange::new();
    let mut potential = Dd::ZERO;
    let g = input.loads.gravity.acceleration_world();
    for ((&before, &after), &properties) in input
        .states
        .iter()
        .zip(coupled.states.iter())
        .zip(input.properties)
    {
        if requested() {
            return Err(RigidDampingThermalError::Cancelled);
        }
        kinetic
            .add_body(properties, before, after)
            .map_err(ViscousNetworkMotionError::from)?;
        let a = before.pose().position_world();
        let b = after.pose().position_world();
        for (initial, final_, acceleration) in [(a.x, b.x, g.x), (a.y, b.y, g.y), (a.z, b.z, g.z)] {
            potential = potential
                - Dd::from_f64(properties.mass())
                    * Dd::from_f64(acceleration)
                    * (Dd::from_f64(final_) - Dd::from_f64(initial));
        }
    }
    let memory_change = number(
        coupled
            .mechanics
            .memory_energy_change::<O::Error, S::Error>()?,
    )?;
    let mechanical_change = kinetic
        .finish_retained()
        .map_err(ViscousNetworkMotionError::from)?
        + potential
        + Dd::from_f64(memory_change);
    let accounts = &coupled.accounts;
    let mechanical = mechanical_change + accounts.heat - accounts.work - accounts.actuator;
    let handoff = accounts.deposited + accounts.unallocated - accounts.heat;
    let combined = mechanical + Dd::from_f64(audit.residual_j) + handoff;
    for residual in [
        mechanical,
        mechanical - accounts.defect,
        handoff,
        combined,
        accounts.absolute_defect,
        accounts.absolute_handoff,
    ] {
        balance(value(residual)?, config.energy_tolerance_j)?;
    }
    let report = RigidDampingThermalReport {
        time_s: target_s,
        thermal: audit,
        connector_heat_j: coupled
            .connector_heat
            .iter()
            .copied()
            .map(value)
            .collect::<Result<Vec<_>, O::Error, S::Error>>()?,
        connector_actuator_work_j: coupled
            .connector_actuator_work
            .iter()
            .copied()
            .map(value)
            .collect::<Result<Vec<_>, O::Error, S::Error>>()?,
        dissipated_energy_j: value(accounts.heat)?,
        deposited_heat_j: value(accounts.deposited)?,
        unallocated_energy_j: value(accounts.unallocated)?,
        external_thermal_heat_j: value(Dd::from_f64(audit.external_heat_j) - accounts.deposited)?,
        mechanical_work_j: value(accounts.work)?,
        actuator_work_j: value(accounts.actuator)?,
        memory_energy_change_j: memory_change,
        mechanical_energy_change_j: value(mechanical_change)?,
        mechanical_residual_j: value(mechanical)?,
        native_drift_energy_defect_j: value(accounts.drift)?,
        absolute_mechanical_defect_j: value(accounts.absolute_defect)?,
        heat_handoff_residual_j: value(handoff)?,
        absolute_heat_handoff_defect_j: value(accounts.absolute_handoff)?,
        combined_energy_residual_j: value(combined)?,
        pair_relaxations: coupled.pair_relaxations,
    };
    if requested() {
        return Err(RigidDampingThermalError::Cancelled);
    }
    Ok((
        RigidDampingThermalCandidate {
            states: coupled.states.into_owned(),
            thermal,
            observer: coupled.inner.into_owned(),
            report,
        },
        coupled.source,
        coupled.mechanics,
        implicit_report,
    ))
}

#[derive(Clone, Default)]
struct Accounts {
    heat: Dd,
    deposited: Dd,
    unallocated: Dd,
    work: Dd,
    actuator: Dd,
    defect: Dd,
    drift: Dd,
    absolute_defect: Dd,
    absolute_handoff: Dd,
}
#[derive(Clone)]
struct Pending {
    start_s: f64,
    dt_s: f64,
    end_s: f64,
    step: MechanicalStep,
    deposited: Dd,
    unallocated: Dd,
    absolute_handoff: Dd,
}
#[derive(Clone)]
struct HeatObserver<'a, 'g, O: Clone, S: DampingSource, M: MechanicalAdvance> {
    input: RigidDampingThermalInput<'a>,
    source: S,
    mechanics: M,
    states: Cow<'a, [RigidBodyState]>,
    time_s: f64,
    inner: Cow<'g, O>,
    remainders: &'g [Dd],
    config: ThermalControls,
    requested: &'g dyn Fn() -> bool,
    pending: Option<Pending>,
    accounts: Accounts,
    connector_heat: Vec<Dd>,
    connector_actuator_work: Vec<Dd>,
    pair_relaxations: usize,
}
impl<O: DepositionObserver, S: DampingSource, M: MechanicalAdvance> DepositionObserver
    for HeatObserver<'_, '_, O, S, M>
{
    type Error = RigidDampingThermalError<O::Error, S::Error>;
    fn begin(&mut self) -> Result<(), O::Error, S::Error> {
        self.inner
            .to_mut()
            .begin()
            .map_err(RigidDampingThermalError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Option<f64>, O::Error, S::Error> {
        self.inner
            .to_mut()
            .next_transition_time_s(thermal, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error, S::Error> {
        collect_source_rates(thermal, |emit| self.visit_sources(thermal, emit, check))
    }
    fn visit_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> core::result::Result<(), DepositionError>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error, S::Error> {
        self.inner
            .to_mut()
            .visit_sources(thermal, emit, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, O::Error, S::Error> {
        self.inner
            .to_mut()
            .conductances(thermal, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error, S::Error> {
        self.inner
            .to_mut()
            .ambient_conductances(thermal, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error, S::Error> {
        self.inner
            .to_mut()
            .reservoir_conductances(thermal, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, O::Error, S::Error> {
        self.inner
            .to_mut()
            .heat_transfers(thermal, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<f64, O::Error, S::Error> {
        self.inner
            .to_mut()
            .limit_timestep_s(thermal, proposed_s, check)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error, S::Error> {
        collect_source_rates(thermal, |emit| {
            self.visit_step_sources(thermal, start_s, dt_s, end_s, emit, check)
        })
    }
    fn visit_step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        emit: &mut dyn FnMut(CellId, f64) -> core::result::Result<(), DepositionError>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error, S::Error> {
        check()?;
        if self.time_s.to_bits() != start_s.to_bits()
            || self
                .pending
                .as_ref()
                .is_some_and(|pending| pending.start_s.to_bits() != start_s.to_bits())
        {
            return Err(RigidDampingThermalError::Invalid(
                "damper preparation must start at the accepted shared clock",
            ));
        }
        // A smaller nonlinear retry replaces only a proposed receipt. The
        // native body bank and constitutive memory still belong to the last
        // accepted step; cumulative preparation work remains charged below.
        self.pending = None;
        for route in self.input.heat_fractions {
            for id in route.keys() {
                check()?;
                if !thermal.cells().contains_key(id) {
                    return Err(RigidDampingThermalError::Invalid(
                        "damper heat receiver is inactive on this interval",
                    ));
                }
            }
        }
        self.source
            .synchronize::<O::Error>(thermal, self.requested)?;
        self.inner
            .to_mut()
            .visit_step_sources(thermal, start_s, dt_s, end_s, emit, check)
            .map_err(RigidDampingThermalError::Observer)?;
        let step = self.mechanics.step_with_source::<O::Error, S>(
            self.input,
            &self.states,
            &self.source,
            dt_s,
            self.requested,
        )?;
        self.pair_relaxations = self
            .pair_relaxations
            .checked_add(step.pair_relaxations)
            .ok_or(DepositionError::Budget("complete damper pair relaxations"))?;
        if step.connector_dissipation_j.len() != self.input.links.len()
            || step
                .connector_actuator_work_j
                .as_ref()
                .is_some_and(|work| work.len() != self.input.links.len())
        {
            return Err(RigidDampingThermalError::Invalid(
                "native work and heat must retain every original connector",
            ));
        }
        let mut edge_work = Dd::ZERO;
        for &work in step.connector_actuator_work_j.iter().flatten() {
            check()?;
            edge_work = edge_work + Dd::from_f64(number(work)?);
        }
        balance(
            value(edge_work - Dd::from_f64(step.actuator_work_j))?,
            self.config.energy_tolerance_j,
        )?;
        let mut allocated = BTreeMap::<CellId, Dd>::new();
        let mut unallocated = Dd::ZERO;
        let mut edge_heat = Dd::ZERO;
        for ((route, &remaining), &heat) in self
            .input
            .heat_fractions
            .iter()
            .zip(self.remainders)
            .zip(&step.connector_dissipation_j)
        {
            check()?;
            if !heat.is_finite() || heat < 0.0 {
                return Err(RigidDampingThermalError::Invalid(
                    "nonnegative native connector heat required",
                ));
            }
            edge_heat = edge_heat + Dd::from_f64(heat);
            let rest = Dd::from_f64(heat) * remaining;
            if value::<O::Error, S::Error>(rest)? == 0.0
                && heat != 0.0
                && (remaining.hi != 0.0 || remaining.lo != 0.0)
            {
                return Err(DepositionError::Arithmetic.into());
            }
            unallocated = unallocated + rest;
            for (&id, &fraction) in route {
                check()?;
                let part = Dd::from_f64(heat) * Dd::from_f64(fraction);
                if value::<O::Error, S::Error>(part)? == 0.0 && heat != 0.0 && fraction != 0.0 {
                    return Err(DepositionError::Arithmetic.into());
                }
                let total = allocated.entry(id).or_insert(Dd::ZERO);
                *total = *total + part;
            }
        }
        let mut deposited = Dd::ZERO;
        let mut ideal = Dd::ZERO;
        let mut absolute = Dd::from_f64(
            value::<O::Error, S::Error>(edge_heat - Dd::from_f64(step.dissipated_energy_j))?.abs(),
        );
        for (id, heat) in allocated {
            check()?;
            ideal = ideal + heat;
            let watts = number(value::<O::Error, S::Error>(heat)? / dt_s)?;
            if watts == 0.0 && (heat.hi != 0.0 || heat.lo != 0.0) {
                return Err(DepositionError::Arithmetic.into());
            }
            emit(id, watts)?;
            let received = Dd::from_f64(watts) * Dd::from_f64(dt_s);
            deposited = deposited + received;
            absolute = absolute + Dd::from_f64(value::<O::Error, S::Error>(received - heat)?.abs());
        }
        absolute = absolute
            + Dd::from_f64(value::<O::Error, S::Error>(ideal + unallocated - edge_heat)?.abs());
        balance(value(absolute)?, self.config.energy_tolerance_j)?;
        self.pending = Some(Pending {
            start_s,
            dt_s,
            end_s,
            step,
            deposited,
            unallocated,
            absolute_handoff: absolute,
        });
        Ok(())
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), O::Error, S::Error> {
        if self.pending.is_some() || self.time_s.to_bits() != thermal.time_s().to_bits() {
            return Err(RigidDampingThermalError::Invalid(
                "event outside accepted shared clock",
            ));
        }
        // All same-time births must finish before requiring every receiver.
        self.inner
            .to_mut()
            .event(event, thermal)
            .map_err(RigidDampingThermalError::Observer)
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error, S::Error> {
        check()?;
        let pending = self
            .pending
            .take()
            .ok_or(RigidDampingThermalError::Invalid(
                "thermal step lacks damper interval",
            ))?;
        if pending.start_s.to_bits() != start_s.to_bits()
            || pending.dt_s.to_bits() != dt_s.to_bits()
            || pending.end_s.to_bits() != thermal.time_s().to_bits()
        {
            return Err(RigidDampingThermalError::Invalid(
                "thermal step changed the mechanical interval",
            ));
        }
        self.inner
            .to_mut()
            .step(start_s, dt_s, old, thermal, check)
            .map_err(RigidDampingThermalError::Observer)?;
        let mut step = pending.step;
        for (total, &heat) in self
            .connector_heat
            .iter_mut()
            .zip(&step.connector_dissipation_j)
        {
            *total = *total + Dd::from_f64(heat);
        }
        if let Some(work) = &step.connector_actuator_work_j {
            for (total, &value) in self.connector_actuator_work.iter_mut().zip(work) {
                *total = *total + Dd::from_f64(value);
            }
        }
        let a = &mut self.accounts;
        a.heat = a.heat + Dd::from_f64(step.dissipated_energy_j);
        a.deposited = a.deposited + pending.deposited;
        a.unallocated = a.unallocated + pending.unallocated;
        a.work = a.work + Dd::from_f64(step.applied_work_j);
        a.actuator = a.actuator + Dd::from_f64(step.actuator_work_j);
        a.defect = a.defect + Dd::from_f64(step.work_heat_residual_j);
        a.drift = a.drift + Dd::from_f64(step.native_drift_energy_defect_j);
        a.absolute_defect = a.absolute_defect + Dd::from_f64(step.work_heat_residual_j.abs());
        a.absolute_handoff = a.absolute_handoff + pending.absolute_handoff;
        balance(value(a.absolute_defect)?, self.config.energy_tolerance_j)?;
        balance(value(a.absolute_handoff)?, self.config.energy_tolerance_j)?;
        self.mechanics.accept::<O::Error, S::Error>(&mut step)?;
        self.states = Cow::Owned(step.states_after);
        self.time_s = pending.end_s;
        check()?;
        Ok(())
    }
}
