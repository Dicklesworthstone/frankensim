//! Material dashpot heat enters the original finite enthalpy timeline.
//!
//! The thermal owner chooses each actual substep, including event/CFL cuts.
//! The original tensor-solid propagator advances once over that interval and
//! supplies its independently integrated element heat, never an energy defect.
//! Full constitutive memory, thermal state and the caller's observer publish
//! together. Mechanical coefficients/geometry remain frozen: this is the
//! explicitly admitted constant-coefficient self-heating approximation, not
//! temperature-dependent viscoelasticity or reciprocal thermoelastic coupling.

/// Closed source-backed relaxation/temperature feedback on the existing
/// implicit thermal residual; stored elastic tensors remain constant.
pub mod relaxation;

use super::{SolidError, SolidStep, TensorSolidHistory};
use fs_exec::CancelGate;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, Event, HeatTransfer, StepConfig, ThermalLaw,
};
use std::collections::{BTreeMap, BTreeSet};

/// One physical heat destination per original element, in connectivity order.
/// The caller owns the spatial correspondence and reference-density agreement;
/// a thermal cell ID is not inferred from a mode number or material name.
#[derive(Clone, Copy)]
pub struct ElementThermalCell<'m> {
    pub cell: CellId,
    /// The EXACT caloric law already borrowed by this thermal cell.
    pub curve: &'m EquilibriumEnthalpyPhaseCurve,
    /// Explicit applicability of the frozen mechanical coefficients [K].
    /// This is a caller admission, not a new measured material-validity claim.
    pub temperature_interval_k: [f64; 2],
}

#[derive(Clone, Copy, Debug)]
pub struct SolidThermalConfig {
    /// Existing explicit enthalpy/CFL owner and whole-call thermal tolerance.
    pub thermal: StepConfig,
    /// Absolute bound for native element heat attribution AND whole-call
    /// mechanical, heat-handoff and combined energy residuals [J].
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum SolidThermalError<E> {
    Invalid(&'static str),
    Solid(SolidError),
    Thermal(DepositionError),
    Phase(PhaseStateError),
    Observer(E),
    Energy { residual_j: f64, allowance_j: f64 },
    Cancelled,
}
impl<E> From<DepositionError> for SolidThermalError<E> {
    fn from(error: DepositionError) -> Self {
        if matches!(error, DepositionError::Cancelled) { Self::Cancelled }
        else { Self::Thermal(error) }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SolidThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tensor solid/thermal coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for SolidThermalError<E> {}
type Result<T, E> = core::result::Result<T, SolidThermalError<E>>;

/// Interval accounts. The raw thermal audit includes generated heat as an
/// incoming source; external_thermal_heat_j removes that INTERNAL transfer.
#[derive(Clone, Debug)]
pub struct SolidThermalReport {
    pub thermal: EnergyAudit,
    pub element_heat_j: Vec<f64>,
    pub material_heat_j: f64,
    pub deposited_heat_j: f64,
    pub external_thermal_heat_j: f64,
    pub mechanical_work_j: f64,
    pub mechanical_energy_change_j: f64,
    pub mechanical_residual_j: f64,
    pub heat_handoff_residual_j: f64,
    pub combined_residual_j: f64,
}

/// No independently publishable thermal or material-memory half-step.
pub struct SolidThermalCandidate<'a, 's, 't, O> {
    solid: TensorSolidHistory<'a, 's>,
    thermal: DepositionSimulation<'t>,
    observer: O,
    report: SolidThermalReport,
}
impl<'a, 's, 't, O> SolidThermalCandidate<'a, 's, 't, O> {
    pub const fn solid(&self) -> &TensorSolidHistory<'a, 's> { &self.solid }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub const fn observer(&self) -> &O { &self.observer }
    pub const fn report(&self) -> &SolidThermalReport { &self.report }
    pub fn into_parts(self) -> (TensorSolidHistory<'a, 's>, DepositionSimulation<'t>, O, SolidThermalReport) {
        (self.solid, self.thermal, self.observer, self.report)
    }
}

fn check<E>(gate: &CancelGate) -> Result<(), E> {
    if gate.is_requested() { Err(SolidThermalError::Cancelled) } else { Ok(()) }
}
fn number<E>(value: f64) -> Result<f64, E> {
    if value.is_finite() { Ok(value) }
    else { Err(SolidThermalError::Invalid("nonfinite coupled heat/work arithmetic")) }
}
fn balance<E>(residual: f64, tolerance: f64) -> Result<(), E> {
    number(residual)?;
    if residual.abs() > tolerance {
        Err(SolidThermalError::Energy { residual_j: residual, allowance_j: tolerance })
    } else { Ok(()) }
}

#[derive(Clone)]
struct Receiver<'m> {
    binding: ElementThermalCell<'m>,
    original: DepositedCell,
}
impl Receiver<'_> {
    fn admit<E>(&self, thermal: &DepositionSimulation<'_>) -> Result<(), E> {
        let cell = thermal.cells().get(&self.binding.cell)
            .ok_or(SolidThermalError::Invalid("material-heat destination is inactive"))?;
        let old = &self.original;
        if cell.id != old.id || cell.material_id != old.material_id
            || cell.segment_id != old.segment_id || cell.body_id != old.body_id
            || cell.mass_kg.to_bits() != old.mass_kg.to_bits() || cell.geometry != old.geometry
            || !thermal.uses_material_law(cell.material_id, self.binding.curve as &dyn ThermalLaw) {
            return Err(SolidThermalError::Invalid("material-heat destination changed its caloric law or reference inventory"));
        }
        let phase = self.binding.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(SolidThermalError::Phase)?;
        let [low, high] = self.binding.temperature_interval_k;
        if !phase.is_fully_solid() || !(low..=high).contains(&phase.temperature_k()) {
            return Err(SolidThermalError::Invalid("material heat left the admitted fully-solid frozen-coefficient interval"));
        }
        Ok(())
    }
}

impl<'a, 's> TensorSolidHistory<'a, 's> {
    /// Advance retained inertia AND all Maxwell history while depositing each
    /// original element's Kelvin/Maxwell loss in its finite caloric receiver.
    /// No scalar modal damping, fitted loss, thermal solver or new material law
    /// is introduced. This works with the existing Tet4/Tet10 projections.
    ///
    /// Forces are held over this call; put force discontinuities at call ends.
    /// The original thermal timeline owns substeps and events. The observer's
    /// source, conductance, pairwise-transfer, timestep, interval and event hooks
    /// remain active, allowing source-backed transport and other finite owners.
    /// Clone must isolate all observer mutation, as required by its native API.
    ///
    /// The exact mechanical interval heat becomes an average thermal power on
    /// THAT interval. Thermal transport retains its original explicit order/CFL;
    /// refine max_dt_s for time accuracy. A whole-call energy check encloses
    /// kinetic, elastic/Maxwell and caloric storage, without counting recovered
    /// branch energy or numerical error as heat. Partial melting, a receiver
    /// replacement, source refusal or cancellation publishes no part of the call.
    ///
    /// Material coefficients, stress-free strain, geometry and supports are
    /// frozen. Temperature bounds must admit that approximation; this does not
    /// update relaxation times, stiffness, expansion or material entropy with T.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermal_candidate<'t, O: DepositionObserver>(
        &self, thermal: &DepositionSimulation<'t>, observer: &O, target_s: f64,
        force_n: &[[f64; 3]], cells: &[ElementThermalCell<'_>],
        config: SolidThermalConfig, gate: &CancelGate,
    ) -> Result<SolidThermalCandidate<'a, 's, 't, O>, O::Error> {
        check(gate)?;
        if self.inner.time_s().to_bits() != thermal.time_s().to_bits()
            || !target_s.is_finite() || target_s <= thermal.time_s()
            || !config.energy_tolerance_j.is_finite() || config.energy_tolerance_j <= 0.0
            || !config.thermal.energy_tolerance_j.is_finite()
            || config.thermal.energy_tolerance_j < 0.0
            || config.thermal.energy_tolerance_j > config.energy_tolerance_j
            || cells.len() != self.solid.projection.element_count() {
            return Err(SolidThermalError::Invalid("one forward shared clock, one receiver per element and bounded energy controls required"));
        }
        // Reuse the original nodal/source dimension and finite-value admission.
        self.solid.projection.project_nodal_force(force_n, gate)
            .map_err(|e| SolidThermalError::Solid(SolidError::Source(e)))?;
        let mut used = BTreeSet::new();
        let mut receivers = Vec::with_capacity(cells.len());
        for &binding in cells {
            check(gate)?;
            let [low, high] = binding.temperature_interval_k;
            if !used.insert(binding.cell) || !low.is_finite() || !high.is_finite()
                || low <= 0.0 || high < low {
                return Err(SolidThermalError::Invalid("unique caloric receivers with explicit finite positive temperature intervals required"));
            }
            let original = thermal.cells().get(&binding.cell)
                .ok_or(SolidThermalError::Invalid("inactive initial material-heat destination"))?.clone();
            let receiver = Receiver { binding, original };
            receiver.admit(thermal)?;
            receivers.push(receiver);
        }
        let coupled = HeatObserver {
            solid: self.clone(), inner: observer.clone(), receivers, force_n, gate,
            tolerance_j: config.energy_tolerance_j, pending: None,
            element_heat_j: vec![0.0; cells.len()], material_heat_j: 0.0,
            deposited_heat_j: 0.0, supplied_work_j: 0.0,
        };
        let (next_thermal, coupled, audit) = thermal.advance_observed(
            &coupled, target_s, config.thermal,
            || if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) },
        )?;
        if coupled.pending.is_some() || coupled.solid.inner.time_s().to_bits() != target_s.to_bits() {
            return Err(SolidThermalError::Invalid("thermal and material-memory endpoints did not publish one clock"));
        }
        let change = number(coupled.solid.inner.stored_energy_j()-self.inner.stored_energy_j())?;
        let external = number(audit.external_heat_j-coupled.deposited_heat_j)?;
        let mechanical = number(change+coupled.material_heat_j-coupled.supplied_work_j)?;
        let handoff = number(coupled.deposited_heat_j-coupled.material_heat_j)?;
        let combined = number(change+audit.stored_energy_change_j-audit.incoming_energy_j
            -external-coupled.supplied_work_j)?;
        for residual in [mechanical, handoff, combined] { balance(residual, config.energy_tolerance_j)?; }
        check(gate)?;
        Ok(SolidThermalCandidate {
            solid: coupled.solid, thermal: next_thermal, observer: coupled.inner,
            report: SolidThermalReport {
                thermal: audit, element_heat_j: coupled.element_heat_j,
                material_heat_j: coupled.material_heat_j, deposited_heat_j: coupled.deposited_heat_j,
                external_thermal_heat_j: external, mechanical_work_j: coupled.supplied_work_j,
                mechanical_energy_change_j: change, mechanical_residual_j: mechanical,
                heat_handoff_residual_j: handoff, combined_residual_j: combined,
            },
        })
    }
}

#[derive(Clone)]
struct Pending<'a, 's> {
    solid: TensorSolidHistory<'a, 's>,
    step: SolidStep,
    start_s: f64,
    dt_s: f64,
    end_s: f64,
    deposited_j: f64,
}
#[derive(Clone)]
struct HeatObserver<'a, 's, 'm, 'g, O> {
    solid: TensorSolidHistory<'a, 's>,
    inner: O,
    receivers: Vec<Receiver<'m>>,
    force_n: &'g [[f64; 3]],
    gate: &'g CancelGate,
    tolerance_j: f64,
    pending: Option<Pending<'a, 's>>,
    element_heat_j: Vec<f64>,
    material_heat_j: f64,
    deposited_heat_j: f64,
    supplied_work_j: f64,
}
impl<O: DepositionObserver> HeatObserver<'_, '_, '_, '_, O> {
    fn admit(&self, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        for receiver in &self.receivers { check(self.gate)?; receiver.admit(thermal)?; }
        Ok(())
    }
}
impl<O: DepositionObserver> DepositionObserver for HeatObserver<'_, '_, '_, '_, O> {
    type Error = SolidThermalError<O::Error>;
    fn begin(&mut self) -> Result<(), O::Error> {
        self.inner.begin().map_err(SolidThermalError::Observer)
    }
    fn requires_step_interval(&self) -> bool { true }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, O::Error> {
        self.admit(thermal)?;
        self.inner.sources(thermal, checkpoint).map_err(SolidThermalError::Observer)
    }
    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, O::Error> {
        self.inner.conductances(thermal, checkpoint).map_err(SolidThermalError::Observer)
    }
    fn heat_transfers(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>)
        -> Result<Vec<HeatTransfer>, O::Error> {
        self.inner.heat_transfers(thermal, checkpoint).map_err(SolidThermalError::Observer)
    }
    fn limit_timestep_s(&mut self, thermal: &DepositionSimulation<'_>, dt_s: f64,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>) -> Result<f64, O::Error> {
        self.inner.limit_timestep_s(thermal, dt_s, checkpoint).map_err(SolidThermalError::Observer)
    }
    fn step_sources(&mut self, thermal: &DepositionSimulation<'_>, start_s: f64, dt_s: f64, end_s: f64,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, O::Error> {
        checkpoint()?;
        self.admit(thermal)?;
        if self.pending.is_some() || self.solid.inner.time_s().to_bits() != start_s.to_bits() {
            return Err(SolidThermalError::Invalid("material-memory interval already pending or stale"));
        }
        let mut rates = self.inner.step_sources(thermal, start_s, dt_s, end_s, checkpoint)
            .map_err(SolidThermalError::Observer)?;
        let mut solid = self.solid.clone();
        let step = solid.advance(self.force_n, dt_s, self.tolerance_j, self.gate)
            .map_err(SolidThermalError::Solid)?;
        if solid.inner.time_s().to_bits() != end_s.to_bits() {
            return Err(SolidThermalError::Invalid("material/thermal interval clock is not exactly representable"));
        }
        let mut deposited = 0.0;
        for (receiver, &heat) in self.receivers.iter().zip(&step.element_heat_j) {
            checkpoint()?;
            if !heat.is_finite() || heat < 0.0 {
                return Err(SolidThermalError::Invalid("negative or nonfinite physical material heat"));
            }
            let watts = number(heat/dt_s)?;
            if heat != 0.0 && watts == 0.0 {
                return Err(SolidThermalError::Invalid("material heat rate underflow"));
            }
            let rate = rates.entry(receiver.binding.cell).or_insert(0.0);
            let old = *rate;
            *rate = number(old+watts)?;
            // Account for the actual representable increment in the merged
            // source map, not an ideal rate lost beside another large source.
            deposited = number(deposited+number((*rate-old)*dt_s)?)?;
        }
        balance(deposited-step.dynamics.dissipated, self.tolerance_j)?;
        self.pending = Some(Pending { solid, step, start_s, dt_s, end_s, deposited_j: deposited });
        Ok(rates)
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        if self.pending.is_some() || self.solid.inner.time_s().to_bits() != thermal.time_s().to_bits() {
            return Err(SolidThermalError::Invalid("event observed outside the shared accepted clock"));
        }
        self.inner.event(event, thermal).map_err(SolidThermalError::Observer)?;
        self.admit(thermal)
    }
    fn step(&mut self, start_s: f64, dt_s: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>) -> Result<(), O::Error> {
        let pending = self.pending.take()
            .ok_or(SolidThermalError::Invalid("thermal step lacks its original material-memory proposal"))?;
        if pending.start_s.to_bits() != start_s.to_bits() || pending.dt_s.to_bits() != dt_s.to_bits()
            || pending.end_s.to_bits() != thermal.time_s().to_bits() {
            return Err(SolidThermalError::Invalid("thermal step changed the admitted material-memory interval"));
        }
        self.inner.step(start_s, dt_s, old, thermal, checkpoint).map_err(SolidThermalError::Observer)?;
        self.admit(thermal)?;
        for (total, heat) in self.element_heat_j.iter_mut().zip(&pending.step.element_heat_j) {
            checkpoint()?; *total = number(*total+heat)?;
        }
        self.material_heat_j = number(self.material_heat_j+pending.step.dynamics.dissipated)?;
        self.deposited_heat_j = number(self.deposited_heat_j+pending.deposited_j)?;
        self.supplied_work_j = number(self.supplied_work_j+pending.step.dynamics.supplied)?;
        check(self.gate)?;
        self.solid = pending.solid;
        Ok(())
    }
}
