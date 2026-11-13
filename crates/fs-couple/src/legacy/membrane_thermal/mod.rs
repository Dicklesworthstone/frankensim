//! Nonlinear membrane damping heat on the original finite enthalpy timeline.
//!
//! Original triangles, reference mass, Kelvin--Voigt source rows or Maxwell
//! histories, inertia and nonlinear mechanical stepping remain with fs-solid.
//! The native thermal owner chooses every conductive/event/source interval;
//! the membrane supplies its native physical loss on that same interval. The original
//! SheetCaloric law converts stored energy to temperature, including variable
//! heat capacity. Numerical mechanical defects never become thermal sources.
//!
//! The baseline entry points hold mechanical coefficients fixed. The
//! `visco::temperature` composition feeds actual receiver temperatures back into
//! declared Maxwell relaxation times while retaining the same energy parameters
//! and natural reference. Every entry requires an explicitly admitted fully solid
//! temperature interval. Changing stress-free geometry, reversible thermoelastic
//! exchange and phase-dependent mechanics need their own thermodynamic closure.
//! Thermal transport retains its native explicit order and CFL limits.

/// Native finite-strain Maxwell relaxation heat on the same enthalpy timeline.
pub mod visco;

use fs_exec::CancelGate;
use fs_material::phase::PhaseStateError;
use fs_material::sheet::thermal::SheetCaloric;
use fs_solid::membrane::{
    MembraneError, MembranePatch,
    dynamics::viscosity::{DampedMembraneMotion, MembraneDampingError},
    inertia::p1_reference_row_masses,
    visco::FiniteViscoMembraneError,
};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, HeatTransfer, StepConfig, ThermalLaw, collect_source_rates,
};
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};

/// One receiver for each ORIGINAL reference triangle, in connectivity order.
/// The caller declares the spatial correspondence and transport surfaces. Mass
/// and volume are checked against the original triangle inventory, never a
/// reduced/free-coordinate mass, deformed area or hot density.
#[derive(Clone, Copy, Debug)]
pub struct ElementThermalCell<'m> {
    pub cell: CellId,
    /// Exact original sheet source and the caloric law borrowed by this cell.
    pub caloric: &'m SheetCaloric,
    /// Explicit applicability of the frozen mechanical coefficients [K].
    /// This is a caller admission, not a measured material-validity claim.
    pub temperature_interval_k: [f64; 2],
}

#[derive(Clone, Copy, Debug)]
pub struct MembraneThermalConfig {
    /// Native transport/event limits; the mechanical step ceiling also applies.
    pub thermal: StepConfig,
    /// Whole-call joule allowance for absolute mechanical numerical defects,
    /// absolute heat handoff error and the combined energy residual.
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum MembraneThermalError<E> {
    Invalid(&'static str),
    Membrane(MembraneError),
    Damping(MembraneDampingError),
    Viscoelastic(FiniteViscoMembraneError),
    Temperature(visco::temperature::ViscoMembraneTemperatureError),
    Thermal(DepositionError),
    Phase(PhaseStateError),
    Observer(E),
    Energy { residual_j: f64, allowance_j: f64 },
    Cancelled,
}
impl<E> From<DepositionError> for MembraneThermalError<E> {
    fn from(error: DepositionError) -> Self {
        if matches!(error, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for MembraneThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "membrane/thermal coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MembraneThermalError<E> {}
type Result<T, E> = core::result::Result<T, MembraneThermalError<E>>;

// Both native mechanical owners supply their own physics and accepted history.
// This common interface only carries the work/heat needed by the original
// thermal timeline, receiver admission and conservative handoff below.
trait ThermalMotion: Clone {
    fn patch(&self) -> &MembranePatch;
    fn time_s(&self) -> f64;
    fn strain_energy_j(&self) -> f64;
    fn kinetic_energy_j(&self) -> f64;
    fn max_step_s(&self) -> f64;
    fn source_temperature_k(&self, element: usize) -> Option<f64>;
    fn element_dissipated_j(&self) -> &[f64];
    fn uses_temperature_feedback(&self) -> bool {
        false
    }
    // Only kinetics may change here: physical storage, reference, momentum and
    // heat are unchanged. Owners with reversible thermal work need a different
    // thermodynamic composition, not this zero-work synchronization hook.
    fn set_temperatures<E>(
        &mut self,
        _temperatures_k: &[f64],
        _cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), E> {
        Ok(())
    }
    fn advance<E>(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<MechanicalInterval, E>;
}

#[derive(Clone)]
struct MechanicalInterval {
    time_s: f64,
    dissipated_j: f64,
    mechanical_work_j: f64,
    // Delta physical mechanical storage + physical heat - external work.
    // It includes signed integration defects, never an extra thermal source.
    energy_defect_j: f64,
    // Sum of absolute constituent defects; opposite numerical terms may not
    // cancel before the whole-call error allowance is checked.
    absolute_energy_defect_j: f64,
    numerical_kinetic_loss_j: f64,
    support_impulse_n_s: [f64; 3],
}

impl ThermalMotion for DampedMembraneMotion<'_, '_> {
    fn patch(&self) -> &MembranePatch {
        self.viscosity().patch()
    }
    fn time_s(&self) -> f64 {
        self.motion().time_s()
    }
    fn strain_energy_j(&self) -> f64 {
        self.motion().strain_energy_j()
    }
    fn kinetic_energy_j(&self) -> f64 {
        self.motion().kinetic_energy_j()
    }
    fn max_step_s(&self) -> f64 {
        self.budget().max_step_s
    }
    fn source_temperature_k(&self, element: usize) -> Option<f64> {
        self.viscosity().laws()[element]
            .source()
            .state_coordinate("T")
    }
    fn element_dissipated_j(&self) -> &[f64] {
        DampedMembraneMotion::element_dissipated_j(self)
    }
    fn advance<E>(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<MechanicalInterval, E> {
        let step = self
            .step(dt_s, dead_load_n, cancelled)
            .map_err(MembraneThermalError::Damping)?;
        Ok(MechanicalInterval {
            time_s: step.time_s,
            dissipated_j: step.dissipated_j,
            mechanical_work_j: step.dead_load_work_j,
            energy_defect_j: step.energy_defect_j,
            absolute_energy_defect_j: step.energy_defect_j.abs(),
            numerical_kinetic_loss_j: 0.0,
            support_impulse_n_s: step.support_impulse_n_s,
        })
    }
}

struct NativeCandidate<'thermal, M, O> {
    motion: M,
    thermal: DepositionSimulation<'thermal>,
    observer: O,
    report: MembraneThermalReport,
}

/// Generated heat is INTERNAL to the combined system. The raw thermal audit
/// counts it as incoming heat; external_thermal_heat_j removes that transfer.
#[derive(Clone, Debug)]
pub struct MembraneThermalReport {
    pub thermal: EnergyAudit,
    pub element_heat_j: Vec<f64>,
    pub material_heat_j: f64,
    /// Represented interval source energy; the native thermal allowance also
    /// bounds rate/storage roundoff after the source is handed off.
    pub deposited_heat_j: f64,
    /// Native aggregate heat less the represented material source. This is
    /// derived from a rounded thermal audit, not a separate external integral.
    pub external_thermal_heat_j: f64,
    pub mechanical_work_j: f64,
    pub mechanical_energy_change_j: f64,
    /// Independently measured mechanical storage change + heat - held-load work.
    pub mechanical_residual_j: f64,
    /// Sum of the native signed numerical defects; never deposited as heat.
    pub mechanical_numerical_defect_j: f64,
    /// Nonnegative time-discretization loss, never deposited into thermal cells.
    /// Zero for the original Verlet/damping composition; explicit for Maxwell BE.
    pub numerical_kinetic_loss_j: f64,
    /// Sum of absolute substep defects, so opposite errors cannot cancel.
    pub absolute_mechanical_defect_j: f64,
    pub heat_handoff_residual_j: f64,
    pub absolute_heat_handoff_defect_j: f64,
    /// Mechanical storage + thermal storage - external work/heat/inventory.
    pub combined_residual_j: f64,
    /// Original elastic/viscous support impulses over the complete call.
    pub support_impulse_n_s: [f64; 3],
}

/// Mechanical history, thermal inventory and observer publish together.
pub struct MembraneThermalCandidate<'model, 'patch, 'thermal, O> {
    motion: DampedMembraneMotion<'model, 'patch>,
    thermal: DepositionSimulation<'thermal>,
    observer: O,
    report: MembraneThermalReport,
}
impl<'model, 'patch, 'thermal, O> MembraneThermalCandidate<'model, 'patch, 'thermal, O> {
    pub const fn motion(&self) -> &DampedMembraneMotion<'model, 'patch> {
        &self.motion
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &MembraneThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DampedMembraneMotion<'model, 'patch>,
        DepositionSimulation<'thermal>,
        O,
        MembraneThermalReport,
    ) {
        (self.motion, self.thermal, self.observer, self.report)
    }
}

fn check<E>(gate: &CancelGate) -> Result<(), E> {
    if gate.is_requested() {
        Err(MembraneThermalError::Cancelled)
    } else {
        Ok(())
    }
}
fn number<E>(value: f64) -> Result<f64, E> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MembraneThermalError::Invalid(
            "nonfinite membrane/thermal energy arithmetic",
        ))
    }
}
fn balance<E>(residual: f64, tolerance: f64) -> Result<(), E> {
    number(residual)?;
    if residual.abs() > tolerance {
        Err(MembraneThermalError::Energy {
            residual_j: residual,
            allowance_j: tolerance,
        })
    } else {
        Ok(())
    }
}

// This is only an arithmetic admission allowance, not uncertain mass/geometry.
fn same_measure(actual: f64, expected: f64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && actual > 0.0
        && expected > 0.0
        && (actual - expected).abs() / actual.max(expected) <= 512.0 * f64::EPSILON
}

// Accumulate physical increments, not large absolute formation-energy datums.
#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), E> {
        number(value)?;
        let next = number(self.sum + value)?;
        let lost = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction = number(self.correction + lost)?;
        self.sum = next;
        Ok(())
    }
    // Keep the low part until all cancelling storage/work terms have entered.
    // Reducing an account to one f64 before combining can erase physical heat.
    fn add_account<E>(&mut self, account: &Self, sign: f64) -> Result<(), E> {
        self.add(sign * account.sum)?;
        self.add(sign * account.correction)
    }
    fn value<E>(&self) -> Result<f64, E> {
        number(self.sum + self.correction)
    }
}

#[derive(Clone)]
struct Receiver<'m> {
    binding: ElementThermalCell<'m>,
    original: DepositedCell,
}
impl Receiver<'_> {
    fn admit<E>(&self, thermal: &DepositionSimulation<'_>) -> Result<(), E> {
        let cell = thermal
            .cells()
            .get(&self.binding.cell)
            .ok_or(MembraneThermalError::Invalid(
                "membrane heat receiver is inactive",
            ))?;
        let old = &self.original;
        if cell.id != old.id
            || cell.material_id != old.material_id
            || cell.segment_id != old.segment_id
            || cell.body_id != old.body_id
            || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
            || cell.geometry != old.geometry
            || !thermal.uses_material_law(
                cell.material_id,
                self.binding.caloric.curve() as &dyn ThermalLaw,
            )
        {
            return Err(MembraneThermalError::Invalid(
                "membrane receiver changed its reference inventory or original caloric law",
            ));
        }
        let phase = self
            .binding
            .caloric
            .curve()
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(MembraneThermalError::Phase)?;
        let [low, high] = self.binding.temperature_interval_k;
        if !phase.is_fully_solid() || !(low..=high).contains(&phase.temperature_k()) {
            return Err(MembraneThermalError::Invalid(
                "membrane heat left the admitted fully-solid frozen-coefficient interval",
            ));
        }
        Ok(())
    }
}

/// Advance finite-amplitude damped motion and every thermal cell on ONE clock.
///
/// Loads are held over this call. The native thermal timeline cuts intervals at
/// events, source transitions, conductive limits and the mechanical step ceiling.
/// Each triangle's analytic damping heat is supplied as its average power on the
/// actual interval. Refinement is still required for nonlinear splitting and
/// thermal transport accuracy; exact local loss flows are not a global exact solve.
///
/// Every receiver uses the SAME complete original sheet source as its triangle
/// and the exact SheetCaloric law already installed in the thermal simulation.
/// Reference mass is the original P1 density integral; reference volume is A*h.
/// The 512-epsilon relative inventory allowance covers floating arithmetic only.
/// Subsequent receiver inventory/law checks are exact.
///
/// The mechanical checkpoint is cloned once on its first actual substep, then
/// its prepared damping/inertia buffers are reused. The native thermal owner
/// clones its state and observer once. Observer Clone must isolate mutable state,
/// as required by DepositionObserver. Cancellation, invalid heat/storage, an
/// observer error or any energy refusal returns no independently publishable
/// mechanical, thermal or observer half-step.
///
/// All mechanical material coefficients and reference geometry remain frozen.
/// A bounded temperature interval admits that approximation, including subsequent
/// calls with the accepted candidate. This does not rebind stiffness/viscosity
/// to the evolving temperature or account for reversible thermoelastic work.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'model, 'patch, 'thermal, O: DepositionObserver>(
    motion: &DampedMembraneMotion<'model, 'patch>,
    thermal: &DepositionSimulation<'thermal>,
    observer: &O,
    target_s: f64,
    dead_load_n: &[[f64; 3]],
    cells: &[ElementThermalCell<'_>],
    config: MembraneThermalConfig,
    gate: &CancelGate,
) -> Result<MembraneThermalCandidate<'model, 'patch, 'thermal, O>, O::Error> {
    let candidate = advance_native(
        motion,
        thermal,
        observer,
        target_s,
        dead_load_n,
        cells,
        config,
        gate,
    )?;
    Ok(MembraneThermalCandidate {
        motion: candidate.motion,
        thermal: candidate.thermal,
        observer: candidate.observer,
        report: candidate.report,
    })
}

#[allow(clippy::too_many_arguments)]
fn advance_native<'thermal, M: ThermalMotion, O: DepositionObserver>(
    motion: &M,
    thermal: &DepositionSimulation<'thermal>,
    observer: &O,
    target_s: f64,
    dead_load_n: &[[f64; 3]],
    cells: &[ElementThermalCell<'_>],
    config: MembraneThermalConfig,
    gate: &CancelGate,
) -> Result<NativeCandidate<'thermal, M, O>, O::Error> {
    check(gate)?;
    let patch = motion.patch();
    if motion.time_s().to_bits() != thermal.time_s().to_bits()
        || !target_s.is_finite()
        || target_s <= thermal.time_s()
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j <= 0.0
        || !config.thermal.energy_tolerance_j.is_finite()
        || config.thermal.energy_tolerance_j < 0.0
        || config.thermal.energy_tolerance_j > config.energy_tolerance_j
        || cells.len() != patch.sheets().len()
        || dead_load_n.len() != patch.mesh().nodes.len()
        || dead_load_n.iter().flatten().any(|x| !x.is_finite())
    {
        return Err(MembraneThermalError::Invalid(
            "one forward shared clock, complete loads/receivers and bounded energy controls required",
        ));
    }
    let mut used = BTreeSet::new();
    let mut receivers = Vec::with_capacity(cells.len());
    for (index, (&binding, area)) in cells
        .iter()
        .zip(patch.reference_triangle_areas_m2())
        .enumerate()
    {
        check(gate)?;
        let [low, high] = binding.temperature_interval_k;
        let sheet = &patch.sheets()[index];
        if !used.insert(binding.cell)
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high < low
            || binding.caloric.sheet().material() != sheet.material()
            || binding.caloric.sheet().thickness_m().to_bits() != sheet.thickness_m().to_bits()
            || binding.caloric.sheet().areal_density_kg_m2().to_bits()
                != sheet.areal_density_kg_m2().to_bits()
        {
            return Err(MembraneThermalError::Invalid(
                "unique receivers, original complete sheet sources and positive temperature intervals required",
            ));
        }
        if let Some(source_t) = motion.source_temperature_k(index) {
            if !(low..=high).contains(&source_t) {
                return Err(MembraneThermalError::Invalid(
                    "frozen mechanical source temperature lies outside the admitted interval",
                ));
            }
        }
        let original = thermal
            .cells()
            .get(&binding.cell)
            .ok_or(MembraneThermalError::Invalid(
                "activate original membrane heat receivers before coupling",
            ))?
            .clone();
        let row_masses =
            p1_reference_row_masses(area, patch.reference_corner_areal_density_kg_m2()[index])
                .map_err(MembraneThermalError::Membrane)?;
        let mass: f64 = row_masses.iter().sum();
        if !same_measure(original.mass_kg, mass)
            || !same_measure(original.geometry.volume_m3, area * sheet.thickness_m())
        {
            return Err(MembraneThermalError::Invalid(
                "membrane receiver must retain the full original triangle mass and reference volume",
            ));
        }
        let receiver = Receiver { binding, original };
        receiver.admit(thermal)?;
        receivers.push(receiver);
    }
    let coupled = HeatObserver {
        motion: Cow::Borrowed(motion),
        inner: observer.clone(),
        receivers,
        dead_load_n,
        gate,
        tolerance_j: config.energy_tolerance_j,
        pending: None,
        element_heat_j: vec![Sum::default(); cells.len()],
        material_heat_j: Sum::default(),
        deposited_heat_j: Sum::default(),
        supplied_work_j: Sum::default(),
        mechanical_defect_j: Sum::default(),
        numerical_kinetic_loss_j: Sum::default(),
        absolute_mechanical_defect_j: Sum::default(),
        absolute_handoff_defect_j: Sum::default(),
        support_impulse_n_s: [Sum::default(); 3],
    };
    let (next_thermal, coupled, audit) =
        thermal.advance_observed(&coupled, target_s, config.thermal, || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })?;
    if coupled.pending.is_some() || coupled.motion.time_s().to_bits() != target_s.to_bits() {
        return Err(MembraneThermalError::Invalid(
            "membrane and thermal endpoints did not publish one clock",
        ));
    }
    let material_heat_j = coupled.material_heat_j.value()?;
    let deposited_heat_j = coupled.deposited_heat_j.value()?;
    let mechanical_work_j = coupled.supplied_work_j.value()?;
    let mut change = Sum::default();
    for value in [
        coupled.motion.strain_energy_j(),
        -motion.strain_energy_j(),
        coupled.motion.kinetic_energy_j(),
        -motion.kinetic_energy_j(),
    ] {
        change.add(value)?;
    }
    let mechanical_energy_change_j = change.value()?;
    let mut mechanical = change;
    mechanical.add_account(&coupled.material_heat_j, 1.0)?;
    mechanical.add_account(&coupled.supplied_work_j, -1.0)?;
    let mechanical_residual_j = mechanical.value()?;
    let mechanical_numerical_defect_j = coupled.mechanical_defect_j.value()?;
    let numerical_kinetic_loss_j = coupled.numerical_kinetic_loss_j.value()?;
    let absolute_mechanical_defect_j = coupled.absolute_mechanical_defect_j.value()?;
    let mut defect_agreement = mechanical;
    defect_agreement.add_account(&coupled.mechanical_defect_j, -1.0)?;
    let mut handoff = coupled.deposited_heat_j;
    handoff.add_account(&coupled.material_heat_j, -1.0)?;
    let heat_handoff_residual_j = handoff.value()?;
    let absolute_heat_handoff_defect_j = coupled.absolute_handoff_defect_j.value()?;
    let mut external = Sum::default();
    external.add(audit.external_heat_j)?;
    external.add_account(&coupled.deposited_heat_j, -1.0)?;
    let external_thermal_heat_j = external.value()?;
    // The native thermal residual subtracts original cell enthalpies before
    // reduction. Preserve our account corrections until after that addition.
    let mut combined = mechanical;
    combined.add(audit.residual_j)?;
    combined.add_account(&handoff, 1.0)?;
    let combined_residual_j = combined.value()?;
    for residual in [
        mechanical_residual_j,
        defect_agreement.value()?,
        absolute_mechanical_defect_j,
        heat_handoff_residual_j,
        absolute_heat_handoff_defect_j,
        combined_residual_j,
    ] {
        balance(residual, config.energy_tolerance_j)?;
    }
    let element_heat_j = coupled
        .element_heat_j
        .iter()
        .map(Sum::value)
        .collect::<Result<Vec<_>, O::Error>>()?;
    let mut support_impulse_n_s = [0.0; 3];
    for (out, sum) in support_impulse_n_s
        .iter_mut()
        .zip(&coupled.support_impulse_n_s)
    {
        *out = sum.value()?;
    }
    check(gate)?;
    Ok(NativeCandidate {
        motion: coupled.motion.into_owned(),
        thermal: next_thermal,
        observer: coupled.inner,
        report: MembraneThermalReport {
            thermal: audit,
            element_heat_j,
            material_heat_j,
            deposited_heat_j,
            external_thermal_heat_j,
            mechanical_work_j,
            mechanical_energy_change_j,
            mechanical_residual_j,
            mechanical_numerical_defect_j,
            numerical_kinetic_loss_j,
            absolute_mechanical_defect_j,
            heat_handoff_residual_j,
            absolute_heat_handoff_defect_j,
            combined_residual_j,
            support_impulse_n_s,
        },
    })
}

#[derive(Clone)]
struct Pending {
    start_s: f64,
    dt_s: f64,
    end_s: f64,
    step: MechanicalInterval,
    deposited_j: f64,
    handoff_defect_j: f64,
}

#[derive(Clone)]
struct HeatObserver<'call, 'm, 'g, M: ThermalMotion, O> {
    motion: Cow<'call, M>,
    inner: O,
    receivers: Vec<Receiver<'m>>,
    dead_load_n: &'g [[f64; 3]],
    gate: &'g CancelGate,
    tolerance_j: f64,
    pending: Option<Pending>,
    element_heat_j: Vec<Sum>,
    material_heat_j: Sum,
    deposited_heat_j: Sum,
    supplied_work_j: Sum,
    mechanical_defect_j: Sum,
    numerical_kinetic_loss_j: Sum,
    absolute_mechanical_defect_j: Sum,
    absolute_handoff_defect_j: Sum,
    support_impulse_n_s: [Sum; 3],
}
impl<M: ThermalMotion, O: DepositionObserver> HeatObserver<'_, '_, '_, M, O> {
    fn admit(&self, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        for receiver in &self.receivers {
            check(self.gate)?;
            receiver.admit(thermal)?;
        }
        Ok(())
    }

    fn update_temperatures(&mut self, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        if !self.motion.uses_temperature_feedback() {
            return Ok(());
        }
        let mut temperatures = Vec::with_capacity(self.receivers.len());
        for receiver in &self.receivers {
            check(self.gate)?;
            temperatures.push(thermal.temperature_k(receiver.binding.cell)?);
        }
        let gate = self.gate;
        self.motion
            .to_mut()
            .set_temperatures(&temperatures, &mut || gate.is_requested())
    }
}
impl<M: ThermalMotion, O: DepositionObserver> DepositionObserver
    for HeatObserver<'_, '_, '_, M, O>
{
    type Error = MembraneThermalError<O::Error>;
    fn begin(&mut self) -> Result<(), O::Error> {
        self.inner.begin().map_err(MembraneThermalError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Option<f64>, O::Error> {
        self.inner
            .next_transition_time_s(thermal, checkpoint)
            .map_err(MembraneThermalError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        collect_source_rates(thermal, |emit| self.visit_sources(thermal, emit, checkpoint))
    }
    fn visit_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> core::result::Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error> {
        self.admit(thermal)?;
        self.inner
            .visit_sources(thermal, emit, checkpoint)
            .map_err(MembraneThermalError::Observer)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, O::Error> {
        self.inner
            .conductances(thermal, checkpoint)
            .map_err(MembraneThermalError::Observer)
    }
    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        self.inner
            .ambient_conductances(thermal, checkpoint)
            .map_err(MembraneThermalError::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        self.inner
            .reservoir_conductances(thermal, checkpoint)
            .map_err(MembraneThermalError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, O::Error> {
        self.inner
            .heat_transfers(thermal, checkpoint)
            .map_err(MembraneThermalError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<f64, O::Error> {
        let capped = proposed_s.min(self.motion.max_step_s());
        let dt = self
            .inner
            .limit_timestep_s(thermal, capped, checkpoint)
            .map_err(MembraneThermalError::Observer)?;
        if !dt.is_finite() || dt <= 0.0 || dt > capped {
            return Err(MembraneThermalError::Invalid(
                "observer must retain the mechanical timestep ceiling",
            ));
        }
        Ok(dt)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, O::Error> {
        collect_source_rates(thermal, |emit| {
            self.visit_step_sources(thermal, start_s, dt_s, end_s, emit, checkpoint)
        })
    }
    fn visit_step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        emit: &mut dyn FnMut(CellId, f64) -> core::result::Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error> {
        checkpoint()?;
        self.admit(thermal)?;
        if self.pending.is_some() || self.motion.time_s().to_bits() != start_s.to_bits() {
            return Err(MembraneThermalError::Invalid(
                "membrane interval already pending or stale",
            ));
        }
        self.inner
            .visit_step_sources(thermal, start_s, dt_s, end_s, emit, checkpoint)
            .map_err(MembraneThermalError::Observer)?;
        self.update_temperatures(thermal)?;
        // This observer is private to the whole native thermal candidate.
        // A later refusal discards it, so no additional per-substep clone is needed.
        let step = self
            .motion
            .to_mut()
            .advance(dt_s, self.dead_load_n, &mut || self.gate.is_requested())?;
        if step.time_s.to_bits() != end_s.to_bits() {
            return Err(MembraneThermalError::Invalid(
                "membrane/thermal interval clock is not exactly representable",
            ));
        }
        let mut heat_total = Sum::default();
        let mut deposited = Sum::default();
        let mut absolute_error = Sum::default();
        for (receiver, &heat) in self
            .receivers
            .iter()
            .zip(self.motion.element_dissipated_j())
        {
            checkpoint()?;
            if !heat.is_finite() || heat < 0.0 {
                return Err(MembraneThermalError::Invalid(
                    "negative or nonfinite membrane material heat",
                ));
            }
            let watts = number(heat / dt_s)?;
            if heat != 0.0 && watts == 0.0 {
                return Err(MembraneThermalError::Invalid(
                    "membrane heat rate underflow",
                ));
            }
            // The native thermal owner combines each independent contribution
            // without rounding it beside the caller's existing sources.
            emit(receiver.binding.cell, watts)?;
            let received = number(watts * dt_s)?;
            heat_total.add(heat)?;
            deposited.add(received)?;
            absolute_error.add(number(received - heat)?.abs())?;
        }
        let attribution_error = number(heat_total.value()? - step.dissipated_j)?.abs();
        absolute_error.add(attribution_error)?;
        let handoff_defect_j = absolute_error.value()?;
        balance(handoff_defect_j, self.tolerance_j)?;
        self.pending = Some(Pending {
            start_s,
            dt_s,
            end_s,
            step,
            deposited_j: deposited.value()?,
            handoff_defect_j,
        });
        Ok(())
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), O::Error> {
        if self.pending.is_some() || self.motion.time_s().to_bits() != thermal.time_s().to_bits() {
            return Err(MembraneThermalError::Invalid(
                "event outside the shared accepted clock",
            ));
        }
        self.inner
            .event(event, thermal)
            .map_err(MembraneThermalError::Observer)?;
        self.admit(thermal)?;
        self.update_temperatures(thermal)
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(), O::Error> {
        let pending = self.pending.take().ok_or(MembraneThermalError::Invalid(
            "thermal step lacks its original membrane interval",
        ))?;
        if pending.start_s.to_bits() != start_s.to_bits()
            || pending.dt_s.to_bits() != dt_s.to_bits()
            || pending.end_s.to_bits() != thermal.time_s().to_bits()
        {
            return Err(MembraneThermalError::Invalid(
                "thermal step changed the admitted membrane interval",
            ));
        }
        self.inner
            .step(start_s, dt_s, old, thermal, checkpoint)
            .map_err(MembraneThermalError::Observer)?;
        self.admit(thermal)?;
        self.update_temperatures(thermal)?;
        for (total, &heat) in self
            .element_heat_j
            .iter_mut()
            .zip(self.motion.element_dissipated_j())
        {
            checkpoint()?;
            total.add(heat)?;
        }
        self.material_heat_j.add(pending.step.dissipated_j)?;
        self.deposited_heat_j.add(pending.deposited_j)?;
        self.supplied_work_j.add(pending.step.mechanical_work_j)?;
        self.mechanical_defect_j.add(pending.step.energy_defect_j)?;
        self.numerical_kinetic_loss_j
            .add(pending.step.numerical_kinetic_loss_j)?;
        self.absolute_mechanical_defect_j
            .add(pending.step.absolute_energy_defect_j)?;
        self.absolute_handoff_defect_j
            .add(pending.handoff_defect_j)?;
        balance(self.absolute_mechanical_defect_j.value()?, self.tolerance_j)?;
        balance(self.absolute_handoff_defect_j.value()?, self.tolerance_j)?;
        for (total, impulse) in self
            .support_impulse_n_s
            .iter_mut()
            .zip(pending.step.support_impulse_n_s)
        {
            total.add(impulse)?;
        }
        check(self.gate)
    }
}

#[cfg(test)]
mod tests {
    use super::{MembraneThermalError, Sum, balance};

    #[test]
    fn heat_account_keeps_small_loss_when_large_stored_energy_or_work_cancels() {
        let mut heat = Sum::default();
        heat.add::<()>(1e16).unwrap();
        heat.add::<()>(1.0).unwrap();
        // A displayed f64 cannot retain this joule beside the large total.
        assert_eq!(heat.value::<()>().unwrap(), 1e16);
        let mut residual = Sum::default();
        residual.add::<()>(-1e16).unwrap();
        residual.add_account::<()>(&heat, 1.0).unwrap();
        assert_eq!(residual.value::<()>().unwrap(), 1.0);
        assert!(matches!(
            balance::<()>(residual.value::<()>().unwrap(), 0.5),
            Err(MembraneThermalError::Energy { .. })
        ));

        // The opposite account sign must retain the same missing work.
        let mut stored = Sum::default();
        stored.add::<()>(1e16).unwrap();
        stored.add_account::<()>(&heat, -1.0).unwrap();
        assert_eq!(stored.value::<()>().unwrap(), -1.0);
    }
}
