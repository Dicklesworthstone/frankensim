//! Native sparse Kelvin/Maxwell motion and finite-cell self-heating.
//!
//! One original prepared mechanical interval supplies independently integrated
//! heat from every original element. The existing enthalpy timeline distributes
//! that interval's mean power over its own CFL/event/source-limited substeps.
//! No modal basis, dense state exponential or replacement material law enters.
//! The whole call is a candidate: motion, quadrature memory, thermal cells and
//! the caller's observer publish together after the physical balance checks.
//!
//! Mechanical coefficients, geometry and reference mass are frozen. Explicit
//! receiver temperature intervals admit that approximation; they do not supply
//! temperature-dependent relaxation or reciprocal thermoelasticity. Refine the
//! native mechanical duration as well as thermal substeps for time accuracy.

use super::{
    ElementThermalCell, Receiver, Result, SolidThermalConfig, SolidThermalError,
    SolidThermalReport, balance, number,
};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_solid::linear3::visco::{ViscoDynamicState, ViscoDynamicStep, ViscoElementError};
use fs_solid::linear3::{QuadraticThermalStrainField, TetElasticError};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event, HeatTransfer,
    collect_source_rates,
};
use std::collections::{BTreeMap, BTreeSet};

/// Original native mechanical inputs over one prepared interval.
#[derive(Clone, Copy)]
pub struct SpatialViscoLoading<'a> {
    /// Actual interval-mean force in full node-major xyz order [N].
    pub mean_force_n: &'a [f64],
    /// Optional prescribed endpoint velocities, zero on free DOFs [m/s].
    pub support_velocity_mps: Option<&'a [f64]>,
    /// Optional externally prescribed eigenstrain at the original material
    /// points. Its work remains a separate input, never a heat deposit.
    pub eigenstrain: Option<QuadraticThermalStrainField<'a>>,
}

/// Joint result with the original sparse-solve and constitutive observations.
pub struct SpatialViscoThermalReport {
    pub mechanical: ViscoDynamicStep,
    pub energy: SolidThermalReport,
}

/// No independently published mechanical or thermal half-step.
pub struct SpatialViscoThermalCandidate<'p, 'a, 't, O> {
    solid: ViscoDynamicState<'p, 'a>,
    thermal: DepositionSimulation<'t>,
    observer: O,
    report: SpatialViscoThermalReport,
}
impl<'p, 'a, 't, O> SpatialViscoThermalCandidate<'p, 'a, 't, O> {
    pub const fn solid(&self) -> &ViscoDynamicState<'p, 'a> {
        &self.solid
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &SpatialViscoThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        ViscoDynamicState<'p, 'a>,
        DepositionSimulation<'t>,
        O,
        SpatialViscoThermalReport,
    ) {
        (self.solid, self.thermal, self.observer, self.report)
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), E> {
    cx.checkpoint().map_err(|_| SolidThermalError::Cancelled)
}
fn native<E>(error: ViscoElementError) -> SolidThermalError<E> {
    match error {
        ViscoElementError::Elastic(TetElasticError::Cancelled) => SolidThermalError::Cancelled,
        other => SolidThermalError::Spatial(other),
    }
}
pub(super) fn same_inventory(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

/// Advance one fixed native mechanical interval and its original finite cells.
///
/// Tet4, Tet10, HEX8 (including B-bar), HEX20 and HEX27 use the same native
/// element inventories and quadrature histories. This path inherits the sparse
/// operator's work/coordinate limits and introduces no dense state-count cap.
/// Every element needs one unique, complete reference-mass/volume receiver;
/// source material names never infer a heat destination.
///
/// The thermal observer keeps all source, transport, event and step hooks. Its
/// Clone implementation must isolate mutation, as required by the original
/// enthalpy API. Mechanical forcing changes belong at mechanical call ends.
/// Recovered spring energy, support/eigenstrain work and sparse-solve error are
/// kept separate from the original constitutive heat. A late refusal discards
/// the complete candidate, so retry cannot consume physical heat twice.
#[allow(clippy::too_many_arguments)]
pub fn step_thermal_candidate<'p, 'a, 't, O: DepositionObserver>(
    solid: &ViscoDynamicState<'p, 'a>,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    loading: SpatialViscoLoading<'_>,
    cells: &[ElementThermalCell<'_>],
    config: SolidThermalConfig,
    cx: &Cx<'_>,
) -> Result<SpatialViscoThermalCandidate<'p, 'a, 't, O>, O::Error> {
    poll(cx)?;
    let system = solid.operator().system();
    if solid.time_s().to_bits() != thermal.time_s().to_bits()
        || cells.len() != system.element_count()
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j <= 0.0
        || !config.thermal.energy_tolerance_j.is_finite()
        || config.thermal.energy_tolerance_j < 0.0
        || config.thermal.energy_tolerance_j > config.energy_tolerance_j
    {
        return Err(SolidThermalError::Invalid(
            "one shared native/thermal clock, one receiver per element and bounded energy controls required",
        ));
    }
    let mut used = BTreeSet::new();
    let mut receivers = Vec::with_capacity(cells.len());
    for (index, &binding) in cells.iter().enumerate() {
        poll(cx)?;
        let [low, high] = binding.temperature_interval_k;
        if !used.insert(binding.cell)
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high < low
        {
            return Err(SolidThermalError::Invalid(
                "unique caloric receivers with explicit positive temperature intervals required",
            ));
        }
        let original = thermal
            .cells()
            .get(&binding.cell)
            .ok_or(SolidThermalError::Invalid(
                "inactive spatial material-heat destination",
            ))?
            .clone();
        let mass = system
            .element_reference_mass_kg(index)
            .ok_or(SolidThermalError::Invalid("missing original element mass"))?;
        let volume =
            system
                .element_reference_volume_m3(index)
                .ok_or(SolidThermalError::Invalid(
                    "missing original element volume",
                ))?;
        if !same_inventory(original.mass_kg, mass)
            || !same_inventory(original.geometry.volume_m3, volume)
        {
            return Err(SolidThermalError::Invalid(
                "each caloric cell must contain its original whole-element reference mass and volume",
            ));
        }
        let receiver = Receiver { binding, original };
        receiver.admit(thermal)?;
        receivers.push(receiver);
    }

    // The original candidate advance owns sparse inertia, all material memory,
    // support/eigenstrain work and its independent momentum/energy gates.
    let mut next = solid.clone();
    let step = match (loading.support_velocity_mps, loading.eigenstrain) {
        (None, None) => next.advance(loading.mean_force_n, cx),
        (Some(v), None) => next.advance_with_support_velocity(loading.mean_force_n, v, cx),
        (None, Some(e)) => next.advance_with_thermal_strain(loading.mean_force_n, e, cx),
        (Some(v), Some(e)) => {
            next.advance_with_support_velocity_and_thermal_strain(loading.mean_force_n, v, e, cx)
        }
    }
    .map_err(native)?;
    let target = next.time_s();
    let duration = number(target - thermal.time_s())?;
    if duration <= 0.0 || step.constitutive.element_steps.len() != receivers.len() {
        return Err(SolidThermalError::Invalid(
            "native element loss inventory or endpoint clock changed",
        ));
    }
    let mut element_heat_j = Vec::with_capacity(receivers.len());
    let mut total = 0.0;
    for element in &step.constitutive.element_steps {
        poll(cx)?;
        let heat = number(element.step_dissipated_j)?;
        if heat < 0.0 {
            return Err(SolidThermalError::Invalid(
                "negative native constitutive heat",
            ));
        }
        element_heat_j.push(heat);
        total = number(total + heat)?;
    }
    let heat = number(step.constitutive.step_dissipated_j)?;
    if heat < 0.0 {
        return Err(SolidThermalError::Invalid("negative total material heat"));
    }
    balance(total - heat, config.energy_tolerance_j)?;
    let change = step.retained_energy_change_j();
    let work = step.retained_work_j();
    let mechanical = number((change + Dd::from_f64(heat) - work).to_f64())?;
    balance(mechanical, config.energy_tolerance_j)?;

    let (next_thermal, next_observer, audit, deposited_j) = deposit_element_heat(
        thermal, observer, receivers, &element_heat_j, target, config,
        || cx.checkpoint().is_err(),
    )?;
    let handoff = number(
        (Dd::from_f64(deposited_j) - Dd::from_f64(heat)).to_f64(),
    )?;
    let external = Dd::from_f64(audit.external_heat_j) - Dd::from_f64(deposited_j);
    let combined = number(
        (change + Dd::from_f64(audit.stored_energy_change_j)
            - Dd::from_f64(audit.incoming_energy_j) - external - work)
            .to_f64(),
    )?;
    balance(handoff, config.energy_tolerance_j)?;
    balance(combined, config.energy_tolerance_j)?;
    poll(cx)?;
    Ok(SpatialViscoThermalCandidate {
        solid: next,
        thermal: next_thermal,
        observer: next_observer,
        report: SpatialViscoThermalReport {
            mechanical: step,
            energy: SolidThermalReport {
                thermal: audit,
                element_heat_j,
                material_heat_j: heat,
                deposited_heat_j: deposited_j,
                external_thermal_heat_j: number(external.to_f64())?,
                mechanical_work_j: number(work.to_f64())?,
                mechanical_energy_change_j: number(change.to_f64())?,
                mechanical_residual_j: mechanical,
                heat_handoff_residual_j: handoff,
                combined_residual_j: combined,
            },
        },
    })
}

// Both original motion owners deposit already integrated element loss. The
// caloric owner alone chooses thermal cuts; none of those cuts reruns motion.
#[allow(clippy::too_many_arguments)]
pub(super) fn deposit_element_heat<'t, O: DepositionObserver>(
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    receivers: Vec<Receiver<'_>>,
    element_heat_j: &[f64],
    target_s: f64,
    config: SolidThermalConfig,
    cancelled: impl Fn() -> bool,
) -> Result<(DepositionSimulation<'t>, O, EnergyAudit, f64), O::Error> {
    poll_heat(&cancelled)?;
    let duration = number(target_s - thermal.time_s())?;
    if duration <= 0.0 || receivers.len() != element_heat_j.len() {
        return Err(SolidThermalError::Invalid(
            "one complete positive material-heat interval and receiver inventory required",
        ));
    }
    let mut powers = Vec::with_capacity(element_heat_j.len());
    for &heat in element_heat_j {
        poll_heat(&cancelled)?;
        let heat = number(heat)?;
        let power = number(heat / duration)?;
        if heat < 0.0 || heat > 0.0 && power == 0.0 {
            return Err(SolidThermalError::Invalid(
                "unrepresentable nonnegative element heat rate",
            ));
        }
        powers.push(power);
    }
    let heating = NativeHeat {
        inner: observer.clone(),
        receivers,
        powers,
        cancelled: &cancelled,
        accepted_time_s: thermal.time_s(),
        target_s,
        pending: None,
        deposited_j: 0.0,
    };
    let (next, heating, audit) =
        thermal.advance_observed(&heating, target_s, config.thermal, || {
            if cancelled() { Err(DepositionError::Cancelled) } else { Ok(()) }
        })?;
    if heating.pending.is_some()
        || heating.accepted_time_s.to_bits() != target_s.to_bits()
        || next.time_s().to_bits() != target_s.to_bits()
    {
        return Err(SolidThermalError::Invalid(
            "material and caloric candidates ended at different times",
        ));
    }
    poll_heat(&cancelled)?;
    Ok((next, heating.inner, audit, heating.deposited_j))
}

fn poll_heat<E>(cancelled: &dyn Fn() -> bool) -> Result<(), E> {
    if cancelled() { Err(SolidThermalError::Cancelled) } else { Ok(()) }
}

#[derive(Clone, Copy)]
struct Pending {
    start_s: f64,
    dt_s: f64,
    end_s: f64,
    deposited_j: f64,
}

#[derive(Clone)]
struct NativeHeat<'m, 'c, O> {
    inner: O,
    receivers: Vec<Receiver<'m>>,
    powers: Vec<f64>,
    cancelled: &'c dyn Fn() -> bool,
    accepted_time_s: f64,
    target_s: f64,
    pending: Option<Pending>,
    deposited_j: f64,
}
impl<O: DepositionObserver> NativeHeat<'_, '_, O> {
    fn admit(&self, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        for receiver in &self.receivers {
            poll_heat(self.cancelled)?;
            receiver.admit(thermal)?;
        }
        Ok(())
    }
}
impl<O: DepositionObserver> DepositionObserver for NativeHeat<'_, '_, O> {
    type Error = SolidThermalError<O::Error>;
    fn begin(&mut self) -> Result<(), O::Error> {
        self.inner.begin().map_err(SolidThermalError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Option<f64>, O::Error> {
        self.inner
            .next_transition_time_s(thermal, check)
            .map_err(SolidThermalError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        collect_source_rates(thermal, |emit| self.visit_sources(thermal, emit, check))
    }
    fn visit_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> core::result::Result<(), DepositionError>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error> {
        self.admit(thermal)?;
        self.inner
            .visit_sources(thermal, emit, check)
            .map_err(SolidThermalError::Observer)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, O::Error> {
        self.inner
            .conductances(thermal, check)
            .map_err(SolidThermalError::Observer)
    }
    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        self.inner
            .ambient_conductances(thermal, check)
            .map_err(SolidThermalError::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        self.inner
            .reservoir_conductances(thermal, check)
            .map_err(SolidThermalError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, O::Error> {
        self.inner
            .heat_transfers(thermal, check)
            .map_err(SolidThermalError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        dt_s: f64,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<f64, O::Error> {
        self.inner
            .limit_timestep_s(thermal, dt_s, check)
            .map_err(SolidThermalError::Observer)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
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
    ) -> Result<(), O::Error> {
        check()?;
        self.admit(thermal)?;
        if self.pending.is_some()
            || self.accepted_time_s.to_bits() != start_s.to_bits()
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || end_s <= start_s
            || end_s > self.target_s
        {
            return Err(SolidThermalError::Invalid(
                "stale or invalid finite-cell material-heat substep",
            ));
        }
        self.inner
            .visit_step_sources(thermal, start_s, dt_s, end_s, emit, check)
            .map_err(SolidThermalError::Observer)?;
        let mut deposited = 0.0;
        for (receiver, &watts) in self.receivers.iter().zip(&self.powers) {
            check()?;
            if watts == 0.0 {
                continue;
            }
            emit(receiver.binding.cell, watts)?;
            deposited = number(deposited + number(watts * dt_s)?)?;
        }
        self.pending = Some(Pending {
            start_s,
            dt_s,
            end_s,
            deposited_j: deposited,
        });
        Ok(())
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        if self.pending.is_some() || thermal.time_s().to_bits() != self.accepted_time_s.to_bits() {
            return Err(SolidThermalError::Invalid(
                "event outside the accepted thermal substep",
            ));
        }
        self.inner
            .event(event, thermal)
            .map_err(SolidThermalError::Observer)?;
        self.admit(thermal)
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error> {
        let pending = self.pending.take().ok_or(SolidThermalError::Invalid(
            "thermal step lacks original element-heat proposal",
        ))?;
        if pending.start_s.to_bits() != start_s.to_bits()
            || pending.dt_s.to_bits() != dt_s.to_bits()
            || pending.end_s.to_bits() != thermal.time_s().to_bits()
        {
            return Err(SolidThermalError::Invalid(
                "thermal step changed the admitted element-heat interval",
            ));
        }
        self.inner
            .step(start_s, dt_s, old, thermal, check)
            .map_err(SolidThermalError::Observer)?;
        self.admit(thermal)?;
        self.deposited_j = number(self.deposited_j + pending.deposited_j)?;
        self.accepted_time_s = pending.end_s;
        poll_heat(self.cancelled)
    }
}
