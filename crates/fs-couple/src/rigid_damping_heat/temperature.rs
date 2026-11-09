//! Two-way temperature feedback for moving passive rigid damper networks.
//!
//! One explicitly authored caloric CellId supplies each connector's physical
//! temperature. The source resolves its complete passive form from the ORIGINAL
//! reference matrix at the start of every actual thermal interval. Native
//! mechanics integrates physical heat, and the original thermal owner receives
//! it on that same interval. The accepted endpoint law is then synchronized to
//! the actual stored caloric temperature, after all same-time native events.
//!
//! This is first-order, lagged-temperature feedback. Exact frozen-coefficient
//! damping and tight conservation do not resolve rapidly changing temperature;
//! refine the physical maximum interval independently. Native implicit thermal
//! conduction retains this lagged feedback; no simultaneous thermal/mechanical
//! endpoint equation, elastic parameter work, expansion, moving thermal geometry
//! or measured material validity is implied. Pure viscous coefficient refresh
//! changes no stored energy and therefore generates no extra work or heat.

use super::{
    DampingSource, MechanicalAdvance, NetworkMotion, Result, RigidDampingThermalCandidate,
    RigidDampingImplicitConfig, RigidDampingThermalConfig, RigidDampingThermalError,
    RigidDampingThermalInput, RigidDampingThermalReport, ThermalControls, ThermalDriver,
    advance_with_source_driver, require_implicit_report,
};
use fs_mbd::RigidBodyState;
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    PassiveBushingMatrix, network::ViscousNetworkLink,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionObserver, DepositionSimulation};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionGraph, ImplicitConductionReport, ImplicitConductionWorkspace,
};
use std::{borrow::Cow, collections::BTreeMap};

/// Original calibrated temperature shifts for a passive reference graph.
pub mod shift;

/// Stateless, deterministic temperature law for each original connector.
///
/// `reference` is always the original input form, never an already shifted
/// accepted result. The return type preserves the native strict passive-matrix
/// admission; active, nonsymmetric or unresolved forms cannot be supplied.
/// Source validity and any calibration belong to the implementation. A source
/// may refuse its own temperature domain; the whole coupled candidate then
/// refuses, including during the final endpoint refresh.
///
/// This is a memoryless viscous law. Rate/history-dependent materials, elastic
/// storage and thermal expansion require their original physical owners.
/// Evaluation must not mutate external physical history: callback side effects
/// cannot be rolled back, just as for native thermal observer callbacks.
pub trait TemperatureDampingLaw {
    type Error;
    fn at_temperature(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> core::result::Result<PassiveBushingMatrix, Self::Error>;
}

/// Joint accepted mechanics, caloric state and source at one physical endpoint.
#[must_use]
pub struct TemperatureDampingThermalCandidate<'m, O> {
    native: RigidDampingThermalCandidate<'m, O>,
    links: Vec<ViscousNetworkLink>,
    temperatures_k: Vec<f64>,
    law_evaluations: usize,
}
impl<'m, O> TemperatureDampingThermalCandidate<'m, O> {
    pub fn states(&self) -> &[RigidBodyState] {
        self.native.states()
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
    /// Actual endpoint forms in the original edge order. For continuation,
    /// retain the ORIGINAL input links as references; do not adopt these shifted
    /// outputs as fresh reference calibration.
    pub fn endpoint_links(&self) -> &[ViscousNetworkLink] {
        &self.links
    }
    /// Actual original caloric coordinate for each endpoint connector [K].
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    /// Successful source evaluations; unchanged temperature entries are reused.
    pub const fn law_evaluations(&self) -> usize {
        self.law_evaluations
    }
    pub fn into_parts(
        self,
    ) -> (
        RigidDampingThermalCandidate<'m, O>,
        Vec<ViscousNetworkLink>,
        Vec<f64>,
    ) {
        (self.native, self.links, self.temperatures_k)
    }
}

/// Advance the coupled native owners with actual finite-cell temperature feedback.
///
/// `temperature_cells` has one explicit physical coordinate per original link;
/// it is independent of heat-partition maps. Shared source cells are allowed and
/// sampled once per synchronization. Every coordinate must exist after same-time
/// births and at the final accepted endpoint. Topology, arms, inertias and held
/// loads retain the baseline contract. The original input damping forms remain
/// immutable reference calibration across all substeps and continuation calls.
///
/// Uses the same native mechanics, event/transport grid, observer forwarding,
/// compensated heat routing and independent energy admission as the frozen
/// entry. Failed source/endpoint evaluation publishes neither physical owner.
/// Additional live storage is O(connectors + distinct temperature cells), with
/// one cached graph. Source evaluation occurs only when that connector's actual
/// temperature changes, at most connectors * (accepted intervals + 1) times.
pub fn advance_candidate<'m, O: DepositionObserver, L: TemperatureDampingLaw>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    target_s: f64,
    config: RigidDampingThermalConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<TemperatureDampingThermalCandidate<'m, O>, O::Error, L::Error> {
    advance_with_motion(
        input,
        thermal,
        observer,
        law,
        temperature_cells,
        NetworkMotion(config.mechanical),
        target_s,
        config.into(),
        cancelled,
    )
}

/// Advance temperature-dependent native damping with stiff thermal transport.
///
/// The original source law is sampled at each accepted interval start and at
/// the final accepted endpoint, using the same reference forms and explicit
/// caloric coordinates as [`advance_candidate`]. A thermal retry retains the
/// accepted bodies and source cache and recomputes only the pending mechanical
/// interval. The complete native thermal report accompanies the joint candidate.
///
/// Conduction uses the caller's prepared graph and BE or SDIRK2 solver. Both
/// SDIRK stages receive one whole-interval mechanical heat proposal; this does
/// not raise the order of the first-order lagged temperature feedback or the
/// moving-mechanics split. Independently refine maximum_dt_s. Native source
/// restrictions, event admission and refusal of conductance overrides follow
/// [`super::advance_implicit_candidate`].
#[allow(clippy::too_many_arguments)]
pub fn advance_implicit_candidate<'m, O: DepositionObserver, L: TemperatureDampingLaw>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    target_s: f64,
    conduction: &ImplicitConductionGraph<'_>,
    config: RigidDampingImplicitConfig,
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
        NetworkMotion(config.mechanical),
        target_s,
        config.thermal_controls(),
        config.driver(conduction, workspace),
        cancelled,
    )?;
    Ok((candidate, require_implicit_report(report)?))
}

pub(super) fn advance_with_motion<
    'm,
    O: DepositionObserver,
    L: TemperatureDampingLaw,
    M: MechanicalAdvance,
>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    mechanics: M,
    target_s: f64,
    config: ThermalControls,
    cancelled: impl FnMut() -> bool,
) -> Result<TemperatureDampingThermalCandidate<'m, O>, O::Error, L::Error> {
    advance_with_motion_driver(
        input,
        thermal,
        observer,
        law,
        temperature_cells,
        mechanics,
        target_s,
        config,
        ThermalDriver::Explicit,
        cancelled,
    )
    .map(|(candidate, _)| candidate)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance_with_motion_driver<
    'm,
    O: DepositionObserver,
    L: TemperatureDampingLaw,
    M: MechanicalAdvance,
>(
    input: RigidDampingThermalInput<'_>,
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    law: &L,
    temperature_cells: &[CellId],
    mechanics: M,
    target_s: f64,
    config: ThermalControls,
    driver: ThermalDriver<'_, '_>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (TemperatureDampingThermalCandidate<'m, O>, Option<ImplicitConductionReport>),
    O::Error,
    L::Error,
> {
    if temperature_cells.len() != input.links.len() {
        return Err(RigidDampingThermalError::Invalid(
            "one explicit temperature coordinate per original connector required",
        ));
    }
    let source = TemperatureSource {
        reference: input.links,
        links: Cow::Borrowed(input.links),
        cells: temperature_cells,
        law,
        temperatures_k: Vec::new(),
        evaluations: 0,
    };
    let (native, source, _, report) = advance_with_source_driver(
        input, thermal, observer, source, mechanics, target_s, config, driver, cancelled,
    )?;
    Ok((
        TemperatureDampingThermalCandidate {
            native,
            links: source.links.into_owned(),
            temperatures_k: source.temperatures_k,
            law_evaluations: source.evaluations,
        },
        report,
    ))
}

struct TemperatureSource<'a, L> {
    reference: &'a [ViscousNetworkLink],
    links: Cow<'a, [ViscousNetworkLink]>,
    cells: &'a [CellId],
    law: &'a L,
    temperatures_k: Vec<f64>,
    evaluations: usize,
}
impl<L> Clone for TemperatureSource<'_, L> {
    fn clone(&self) -> Self {
        Self {
            reference: self.reference,
            links: self.links.clone(),
            cells: self.cells,
            law: self.law,
            temperatures_k: self.temperatures_k.clone(),
            evaluations: self.evaluations,
        }
    }
}
impl<L: TemperatureDampingLaw> DampingSource for TemperatureSource<'_, L> {
    type Error = L::Error;
    fn links(&self) -> &[ViscousNetworkLink] {
        &self.links
    }
    fn synchronize<E>(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        requested: &dyn Fn() -> bool,
    ) -> Result<(), E, Self::Error> {
        let mut sampled = BTreeMap::new();
        for (connector, &cell) in self.cells.iter().enumerate() {
            if requested() {
                return Err(RigidDampingThermalError::Cancelled);
            }
            let temperature = match sampled.entry(cell) {
                std::collections::btree_map::Entry::Occupied(entry) => *entry.get(),
                std::collections::btree_map::Entry::Vacant(entry) => {
                    let temperature = thermal.temperature_k(cell)?;
                    if !temperature.is_finite() || temperature <= 0.0 {
                        return Err(RigidDampingThermalError::Invalid(
                            "finite positive original caloric temperature required",
                        ));
                    }
                    *entry.insert(temperature)
                }
            };
            if self.temperatures_k.get(connector) == Some(&temperature) {
                continue;
            }
            let damping = self
                .law
                .at_temperature(connector, &self.reference[connector].damping, temperature)
                .map_err(|error| RigidDampingThermalError::Source { connector, error })?;
            self.links.to_mut()[connector].damping = damping;
            if let Some(previous) = self.temperatures_k.get_mut(connector) {
                *previous = temperature;
            } else {
                self.temperatures_k.push(temperature);
            }
            self.evaluations = self
                .evaluations
                .checked_add(1)
                .ok_or(DepositionError::Budget(
                    "temperature-dependent damper evaluations",
                ))?;
        }
        Ok(())
    }
}
