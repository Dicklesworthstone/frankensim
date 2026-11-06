//! Nuclear source histories on the existing conservative thermal grid.
//!
//! Only the nuclear source's locally deposited energy is credited to its bound
//! cell. Transportable radiation and nonrecoverable export remain separate in
//! the source account. The original thermal law, contact exchange, boundaries,
//! enthalpy/internal-energy basis and thermal audit retain their existing owners.
//!
//! Prescribed-event sources and multigroup neutron irradiation share one clock,
//! material binding and heat-credit transaction. Their authored transitions clip
//! the native step exactly. Nuclear reservoirs advance over the FINAL accepted
//! interval. A spatial diffusion cohort or compartments exchanging nuclides
//! advance one shared physical history before per-cell heat is credited. Thermal transport
//! retains its original explicit timestep and uniform-cell limitations.

use std::{collections::BTreeMap, sync::Arc};

use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalStorageBasis,
};

use crate::burnup::BurnupError;
use crate::decay::DecayError;
use crate::decay::activation::CaptureActivationError;
use crate::diffusion::burnup::DiffusionBurnupError;
use crate::kinetics::KineticsError;
use crate::{EnergyChannels, FissionError, FissionState};

mod source;
pub use source::{FluxDrivenHeatSource, NuclearHeatSource};
mod irradiation;
pub use irradiation::{
    IrradiationCellSource, IrradiationCellStep, IrradiationHeating, IrradiationHeatingBudget,
    NeutronFluxChange, NeutronFluxSchedule,
};
mod burnup;
pub use burnup::{BurnupCellSource, BurnupHeating, BurnupHeatingBudget};
mod burnup_production;
pub use burnup_production::{BurnupProductionCellSource, BurnupProductionHeating};

mod diffusion;
pub use diffusion::{
    DiffusionCellSource, DiffusionCellStep, DiffusionHeating, DiffusionSourceChange,
    DiffusionSourceSchedule,
};

mod transport;
pub use transport::{
    BurnupTransportCellSource, BurnupTransportCellStep, BurnupTransportHeating,
    NuclideTransportChange, NuclideTransportSchedule,
};

mod production;
pub use production::{DecayProductionCellSource, DecayProductionHeating};
mod activation;
pub use activation::{CaptureActivationCellSource, CaptureActivationHeating};

mod kinetics;
pub use kinetics::{
    KineticCellSource, KineticCellStep, KineticHeating, KineticHeatingBudget, KineticsControl,
    KineticsControlChange, KineticsControlSchedule,
};

mod adaptive;
pub use adaptive::{
    AdaptiveNuclearHeatSource, NuclearAdaptiveControls, NuclearAdaptiveError,
    NuclearAdaptiveReport, NuclearTimeError, NuclearTimeTolerances,
};

/// An absolute-time change to a cell's total prescribed fission rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FissionRateChange {
    /// Finite time after zero [s]; changes must be strictly ordered.
    pub time_s: f64,
    /// Finite, nonnegative total number of fissions per second in this cell.
    pub fissions_per_s: f64,
}

/// Immutable right-continuous rate history with binary-search lookup.
///
/// The initial rate applies from zero until the first change. A source bound
/// at a later thermal time uses its explicitly supplied current fuel/reservoirs;
/// it does not replay or infer an earlier irradiation history.
#[derive(Clone, Debug)]
pub struct FissionRateSchedule {
    initial_fissions_per_s: f64,
    changes: Arc<[FissionRateChange]>,
}

impl FissionRateSchedule {
    /// Validate the complete caller-owned schedule before retaining it.
    pub fn new(
        initial_fissions_per_s: f64,
        changes: Vec<FissionRateChange>,
    ) -> Result<Self, FissionThermalError> {
        if !nonnegative(initial_fissions_per_s) {
            return Err(FissionThermalError::Invalid("invalid initial fission rate"));
        }
        let mut previous = 0.0;
        for change in &changes {
            if !change.time_s.is_finite()
                || change.time_s <= previous
                || !nonnegative(change.fissions_per_s)
            {
                return Err(FissionThermalError::Invalid(
                    "rate changes need strictly increasing positive times and nonnegative rates",
                ));
            }
            previous = change.time_s;
        }
        Ok(Self {
            initial_fissions_per_s,
            changes: changes.into(),
        })
    }

    /// Prescribed rate at the given nonnegative finite absolute time [s].
    pub fn rate_at(&self, time_s: f64) -> Result<f64, FissionThermalError> {
        if !nonnegative(time_s) {
            return Err(FissionThermalError::Invalid(
                "invalid fission schedule time",
            ));
        }
        let index = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        Ok(if index == 0 {
            self.initial_fissions_per_s
        } else {
            self.changes[index - 1].fissions_per_s
        })
    }

    /// Authored changes; no mutable schedule can escape into a running source.
    pub fn changes(&self) -> &[FissionRateChange] {
        &self.changes
    }

    fn next_after(&self, time_s: f64) -> Option<f64> {
        let index = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        self.changes.get(index).map(|change| change.time_s)
    }
}

/// One nuclear inventory and rate history bound to one already active cell.
/// A cell has at most one owner; combine its nuclear inventory before binding.
#[derive(Clone, Debug)]
pub struct FissionCellSource {
    /// Original active thermal cell receiving locally deposited heat.
    pub cell_id: CellId,
    /// Original thermal material identity expected at that cell.
    pub material_id: MaterialId,
    /// Complete nuclear history at the binding's current thermal time.
    pub state: FissionState,
    /// Caller-prescribed total fission rate history for this inventory.
    pub rates: FissionRateSchedule,
}

/// Nuclear release and its thermal credit during the last accepted advance.
/// This is an accounting observation, not a physical model-error bound.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FissionHeatingAudit {
    /// Released channels, including channels deliberately not sent to heat.
    pub released: EnergyChannels,
    /// Actual `average_source_w * accepted_dt_s` delivered to the thermal owner.
    pub thermal_credit_j: f64,
    /// Sum of absolute differences between source deposition and thermal credit.
    pub coupling_roundoff_j: f64,
    /// Sum of nuclear numerical energy defects and declared source roundoff.
    pub nuclear_energy_defect_j: f64,
    /// Total accepted expected fission count over this advance.
    pub fissions: f64,
    /// Total accepted competing captures over this advance.
    pub captures: f64,
    /// Integrated radioactive transformations, including newly produced daughters.
    pub decays: f64,
    /// Supplied nuclei integrated over all accepted intervals and cells.
    pub produced_nuclei: f64,
    /// Noncancelling inventory, production and counter-reduction defects.
    pub population_roundoff_nuclei: f64,
    /// Number of common accepted thermal/source intervals.
    pub substeps: usize,
}

/// Refusal preserving the original thermal and nuclear source histories.
#[derive(Clone, Debug, PartialEq)]
pub enum FissionThermalError {
    /// Original thermal-domain, work, material or cancellation refusal.
    Thermal(DepositionError),
    /// Original nuclear source-model or integration refusal.
    Fission(FissionError),
    /// Original neutron-kinetics model or numerical-accuracy refusal.
    Kinetics(KineticsError),
    /// Original product inventory, reaction-data or propagation refusal.
    Burnup(BurnupError),
    /// Original radioactive inventory, production or propagation refusal.
    Decay(DecayError),
    /// Original capture binding, coupled history or propagation refusal.
    Activation(CaptureActivationError),
    /// Spatial neutron, finite-fuel or shared reaction-account refusal.
    Diffusion(DiffusionBurnupError),
    /// Malformed schedule, duplicate owner or inconsistent hook sequence.
    Invalid(&'static str),
    /// A source names a cell that has not been activated.
    InactiveCell(CellId),
    /// A physical cell has a different material than the source binding.
    MaterialMismatch {
        /// Affected cell.
        cell_id: CellId,
        /// Authored binding identity.
        expected: MaterialId,
        /// Current thermal material identity.
        actual: MaterialId,
    },
    /// Material calorics, mass, geometry or physical correspondence changed.
    BindingChanged(CellId),
    /// Independently advanced owners no longer represent the same physical time.
    ClockMismatch {
        /// Current thermal time [s].
        thermal_time_s: f64,
        /// Expected source time [s].
        fission_time_s: f64,
    },
    /// Noncancelling nuclear and coupling defects exceed the advance allowance.
    EnergyBudget {
        /// Noncancelling accumulated numerical defects [J].
        defect_j: f64,
        /// Caller-supplied advance allowance [J].
        tolerance_j: f64,
    },
    /// Accumulated population defects exceed the source's advance allowance.
    NucleiBudget {
        /// Noncancelling numerical population defects [expected nuclei].
        defect_nuclei: f64,
        /// Caller-supplied advance allowance [expected nuclei].
        tolerance_nuclei: f64,
    },
    /// Finite inputs produced nonfinite arithmetic or an exhausted counter.
    Arithmetic,
}

impl From<DepositionError> for FissionThermalError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<FissionError> for FissionThermalError {
    fn from(error: FissionError) -> Self {
        Self::Fission(error)
    }
}
impl From<KineticsError> for FissionThermalError {
    fn from(error: KineticsError) -> Self {
        Self::Kinetics(error)
    }
}
impl From<DecayError> for FissionThermalError {
    fn from(error: DecayError) -> Self {
        Self::Decay(error)
    }
}
impl From<CaptureActivationError> for FissionThermalError {
    fn from(error: CaptureActivationError) -> Self {
        Self::Activation(error)
    }
}
impl From<BurnupError> for FissionThermalError {
    fn from(error: BurnupError) -> Self {
        match error {
            // Preserve the original callback error for the native observer's
            // cancellation/error recovery rather than wrapping it twice.
            BurnupError::Fission(error) => Self::Fission(error),
            other => Self::Burnup(other),
        }
    }
}
impl From<DiffusionBurnupError> for FissionThermalError {
    fn from(error: DiffusionBurnupError) -> Self {
        match error {
            DiffusionBurnupError::Fission(error) => Self::Fission(error),
            other => Self::Diffusion(other),
        }
    }
}
impl std::fmt::Display for FissionThermalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "fission/thermal coupling refused: {self:?}")
    }
}
impl std::error::Error for FissionThermalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(error) => Some(error),
            Self::Fission(error) => Some(error),
            Self::Kinetics(error) => Some(error),
            Self::Burnup(error) => Some(error),
            Self::Decay(error) => Some(error),
            Self::Activation(error) => Some(error),
            Self::Diffusion(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
struct BoundSource<S: NuclearHeatSource> {
    source: S,
    mass_kg: f64,
    geometry: CellGeometry,
    body_id: u64,
    segment_id: u64,
    caloric_identity: Option<[u8; 32]>,
    storage_basis: ThermalStorageBasis,
    last_step: Option<S::Step>,
}

#[derive(Clone, Debug)]
struct PendingStep<S: NuclearHeatSource> {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    states: Vec<(CellId, S, S::Step)>,
    account: Account,
    audit: FissionHeatingAudit,
}

/// Transactional source owner for [`DepositionSimulation::advance_observed`].
///
/// Clone copies all mutable physical state and shares only immutable source
/// models/rate schedules. Use [`Self::advance_to`] for direct operation; compose
/// this observer with other domains using `fs_thermal::observers::CoupledObservers`.
/// Bindings require active cells and fixed reference mass/geometry. Nuclear
/// burnup does not silently change the original material law or thermal mass.
#[derive(Clone, Debug)]
pub struct NuclearHeating<S: NuclearHeatSource> {
    time_s: f64,
    sources: BTreeMap<CellId, BoundSource<S>>,
    budget: S::Budget,
    pending: Option<PendingStep<S>>,
    account: Account,
    audit: FissionHeatingAudit,
}

/// The shared thermal owner specialized for prescribed fission histories.
pub type FissionHeating = NuclearHeating<FissionCellSource>;

impl<S: NuclearHeatSource> NuclearHeating<S> {
    /// Bind current nuclear inventories to the current thermal clock/materials.
    /// Birth events for all source cells must already have been processed.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        sources: Vec<S>,
        budget: S::Budget,
    ) -> Result<Self, FissionThermalError> {
        if sources.is_empty() || sources.len() > thermal.cells().len() {
            return Err(FissionThermalError::Invalid(
                "source bindings need distinct active cells",
            ));
        }
        S::validate_sources(sources.iter(), &mut || Ok(()))?;
        let mut bound = BTreeMap::new();
        for source in sources {
            let id = source.cell_id();
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(FissionThermalError::InactiveCell(id))?;
            if cell.material_id != source.material_id() {
                return Err(FissionThermalError::MaterialMismatch {
                    cell_id: id,
                    expected: source.material_id(),
                    actual: cell.material_id,
                });
            }
            if source
                .expected_volume_m3()
                .is_some_and(|volume| volume != cell.geometry.volume_m3)
            {
                return Err(FissionThermalError::BindingChanged(id));
            }
            // The same kernel admits its work/tolerance budget, even when the
            // first requested thermal advance will have zero duration.
            source.admit(thermal.time_s(), thermal.temperature_k(id)?, budget)?;
            let value = BoundSource {
                source,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                body_id: cell.body_id,
                segment_id: cell.segment_id,
                caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
                storage_basis: thermal.material_storage_basis(cell.material_id)?,
                last_step: None,
            };
            if bound.insert(id, value).is_some() {
                return Err(FissionThermalError::Invalid(
                    "duplicate fission source cell",
                ));
            }
        }
        Ok(Self {
            time_s: thermal.time_s(),
            sources: bound,
            budget,
            pending: None,
            account: Account::default(),
            audit: FissionHeatingAudit::default(),
        })
    }

    /// Current accepted source/thermal clock [s].
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Accepted source inventory and immutable schedule for one cell.
    pub fn source(&self, cell: CellId) -> Option<&S> {
        self.sources.get(&cell).map(|bound| &bound.source)
    }
    /// Accepted source states in deterministic physical-cell order.
    pub fn sources(&self) -> impl ExactSizeIterator<Item = &S> {
        self.sources.values().map(|bound| &bound.source)
    }
    /// Last accepted nuclear interval at a cell, absent before the first step.
    pub fn last_step(&self, cell: CellId) -> Option<&S::Step> {
        self.sources
            .get(&cell)
            .and_then(|bound| bound.last_step.as_ref())
    }
    /// Source/coupling account for the last accepted advance, excluding other heat.
    pub fn audit(&self) -> FissionHeatingAudit {
        self.audit
    }

    /// Sum the original source owners' lifetime release accounts.
    pub fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        let mut totals = [Sum::default(); 3];
        for source in self.sources() {
            let released = source.cumulative_released()?;
            totals[0].add(released.deposited_j)?;
            totals[1].add(released.radiation_j)?;
            totals[2].add(released.escaped_j)?;
        }
        Ok(EnergyChannels {
            deposited_j: totals[0].value()?,
            radiation_j: totals[1].value()?,
            escaped_j: totals[2].value()?,
        })
    }

    /// Advance both owners atomically, including admission for an idle call.
    pub fn advance_to<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'thermal>, Self, EnergyAudit), FissionThermalError> {
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        let (thermal, source, audit) =
            thermal.advance_observed(self, target_s, config, &mut checkpoint)?;
        source.validate(&thermal, source.time_s, &mut checkpoint)?;
        source.admit_audit(&source.audit)?;
        Ok((thermal, source, audit))
    }

    fn admit_audit(&self, audit: &FissionHeatingAudit) -> Result<(), FissionThermalError> {
        let defect = finite(audit.coupling_roundoff_j + audit.nuclear_energy_defect_j)?;
        if defect > S::energy_tolerance_j(self.budget) {
            return Err(FissionThermalError::EnergyBudget {
                defect_j: defect,
                tolerance_j: S::energy_tolerance_j(self.budget),
            });
        }
        if let Some(tolerance_nuclei) = S::nuclei_tolerance(self.budget) {
            if !nonnegative(tolerance_nuclei) {
                return Err(FissionThermalError::Invalid(
                    "population allowance must be finite and nonnegative",
                ));
            }
            if audit.population_roundoff_nuclei > tolerance_nuclei {
                return Err(FissionThermalError::NucleiBudget {
                    defect_nuclei: audit.population_roundoff_nuclei,
                    tolerance_nuclei,
                });
            }
        }
        Ok(())
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        expected_time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        checkpoint()?;
        if thermal.time_s() != expected_time_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s: thermal.time_s(),
                fission_time_s: expected_time_s,
            });
        }
        S::validate_sources(
            self.sources.values().map(|bound| &bound.source),
            &mut || checkpoint().map_err(Into::into),
        )?;
        for (index, (&id, bound)) in self.sources.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(FissionThermalError::InactiveCell(id))?;
            if cell.material_id != bound.source.material_id() {
                return Err(FissionThermalError::MaterialMismatch {
                    cell_id: id,
                    expected: bound.source.material_id(),
                    actual: cell.material_id,
                });
            }
            if cell.mass_kg != bound.mass_kg
                || cell.geometry != bound.geometry
                || cell.body_id != bound.body_id
                || cell.segment_id != bound.segment_id
                || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
                || thermal.material_storage_basis(cell.material_id)? != bound.storage_basis
            {
                return Err(FissionThermalError::BindingChanged(id));
            }
            bound
                .source
                .validate_temperature(thermal.temperature_k(id)?)?;
        }
        Ok(())
    }
}

impl<S: NuclearHeatSource> DepositionObserver for NuclearHeating<S> {
    type Error = FissionThermalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        if self.pending.is_some() {
            return Err(FissionThermalError::Invalid("uncommitted nuclear interval"));
        }
        self.account = Account::default();
        self.audit = FissionHeatingAudit::default();
        Ok(())
    }

    fn requires_step_interval(&self) -> bool {
        true
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.validate(thermal, self.time_s, checkpoint)?;
        Ok(NuclearHeating::sources(self)
            .filter_map(|source| source.next_transition(self.time_s))
            .min_by(f64::total_cmp))
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.validate(thermal, self.time_s, checkpoint)?;
        if self.pending.is_some()
            || start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !end_s.is_finite()
            || end_s <= start_s
        {
            return Err(FissionThermalError::Invalid(
                "invalid nuclear/thermal interval",
            ));
        }
        let mut source_refs = Vec::with_capacity(self.sources.len());
        let mut temperatures = Vec::with_capacity(self.sources.len());
        for (&id, bound) in &self.sources {
            checkpoint()?;
            if bound
                .source
                .next_transition(start_s)
                .is_some_and(|time| time < end_s)
            {
                return Err(FissionThermalError::Invalid(
                    "nuclear interval crosses a rate change",
                ));
            }
            source_refs.push(&bound.source);
            temperatures.push(thermal.temperature_k(id)?);
        }
        let mut checkpoint_error = None;
        let trials = S::trial_sources(
            &source_refs,
            start_s,
            duration_s,
            &temperatures,
            self.budget,
            &mut || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    FissionError::Cancelled
                })
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let trials = trials?;
        if trials.len() != self.sources.len() {
            return Err(FissionThermalError::Invalid(
                "nuclear collection changed cell coverage",
            ));
        }
        let mut states = Vec::with_capacity(self.sources.len());
        let mut watts = BTreeMap::new();
        let mut account = self.account.clone();
        for ((&id, bound), (state, report)) in self.sources.iter().zip(trials) {
            checkpoint()?;
            if state.cell_id() != id || state.material_id() != bound.source.material_id() {
                return Err(FissionThermalError::Invalid(
                    "nuclear collection changed material correspondence",
                ));
            }
            let power = finite(S::released(&report).deposited_j / duration_s)?;
            let credit = finite(power * duration_s)?;
            account.record::<S>(&report, credit)?;
            watts.insert(id, power);
            states.push((id, state, report));
        }
        account.substeps = account
            .substeps
            .checked_add(1)
            .ok_or(FissionThermalError::Arithmetic)?;
        let audit = account.snapshot()?;
        self.admit_audit(&audit)?;
        checkpoint()?;
        self.pending = Some(PendingStep {
            start_s,
            duration_s,
            end_s,
            states,
            account,
            audit,
        });
        Ok(watts)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.validate(thermal, self.time_s, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let pending = self.pending.as_ref().ok_or(FissionThermalError::Invalid(
            "thermal step has no nuclear candidate",
        ))?;
        if pending.start_s != start_s || pending.duration_s != duration_s {
            return Err(FissionThermalError::Invalid(
                "thermal step changed the nuclear interval",
            ));
        }
        self.validate(thermal, pending.end_s, checkpoint)?;
        self.admit_audit(&pending.audit)?;
        checkpoint()?;
        // All fallible work precedes publication, even for direct hook callers.
        let pending = self.pending.take().expect("validated pending interval");
        for (id, state, report) in pending.states {
            let source = self
                .sources
                .get_mut(&id)
                .expect("candidate retains source identities");
            source.source = state;
            source.last_step = Some(report);
        }
        self.time_s = pending.end_s;
        self.account = pending.account;
        self.audit = pending.audit;
        Ok(())
    }
}

// Compensated scalar totals keep many cells/intervals from dropping small
// released channels. Defect sums contain absolute values, so errors cannot cancel.
#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), FissionThermalError> {
        finite(value)?;
        let next = finite(self.sum + value)?;
        let lost = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction = finite(self.correction + lost)?;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> Result<f64, FissionThermalError> {
        finite(self.sum + self.correction)
    }
    fn export_roundoff(&self) -> Result<f64, FissionThermalError> {
        finite(crate::sum(&[self.sum, self.correction, -self.value()?]).abs())
    }
    // Used only for nonnegative defect sums. Round upward when a retained low
    // part would otherwise disappear at an exact caller allowance boundary.
    fn nonnegative_upper_bound(&self) -> Result<f64, FissionThermalError> {
        let rounded = self.value()?;
        let remainder = finite(crate::sum(&[self.sum, self.correction, -rounded]))?;
        if remainder > 0.0 {
            finite(rounded.next_up())
        } else {
            Ok(rounded)
        }
    }
}

#[derive(Clone, Debug, Default)]
struct Account {
    released: [Sum; 3],
    credit: Sum,
    coupling_roundoff: Sum,
    nuclear_defect: Sum,
    fissions: Sum,
    captures: Sum,
    decays: Sum,
    produced: Sum,
    population_roundoff: Sum,
    substeps: usize,
}
impl Account {
    fn record<S: NuclearHeatSource>(
        &mut self,
        report: &S::Step,
        credit: f64,
    ) -> Result<(), FissionThermalError> {
        let released = S::released(report);
        self.released[0].add(released.deposited_j)?;
        self.released[1].add(released.radiation_j)?;
        self.released[2].add(released.escaped_j)?;
        self.credit.add(credit)?;
        self.coupling_roundoff
            .add((credit - released.deposited_j).abs())?;
        // The core's noncancelling roundoff already includes its residual.
        self.nuclear_defect.add(S::energy_roundoff_j(report))?;
        self.fissions.add(S::fissions(report))?;
        self.captures.add(S::captures(report))?;
        self.decays.add(S::decays(report))?;
        self.produced.add(S::produced_nuclei(report))?;
        self.population_roundoff
            .add(S::population_roundoff_nuclei(report))?;
        Ok(())
    }
    fn snapshot(&self) -> Result<FissionHeatingAudit, FissionThermalError> {
        let mut population_roundoff = self.population_roundoff;
        // The retained low parts must not disappear for free when the new
        // count observations are exported as one scalar per advance.
        population_roundoff.add(self.decays.export_roundoff()?)?;
        population_roundoff.add(self.produced.export_roundoff()?)?;
        population_roundoff.add(self.fissions.export_roundoff()?)?;
        population_roundoff.add(self.captures.export_roundoff()?)?;
        let population_roundoff_nuclei = population_roundoff.nonnegative_upper_bound()?;
        Ok(FissionHeatingAudit {
            released: EnergyChannels {
                deposited_j: self.released[0].value()?,
                radiation_j: self.released[1].value()?,
                escaped_j: self.released[2].value()?,
            },
            thermal_credit_j: self.credit.value()?,
            coupling_roundoff_j: self.coupling_roundoff.value()?,
            nuclear_energy_defect_j: self.nuclear_defect.value()?,
            fissions: self.fissions.value()?,
            captures: self.captures.value()?,
            decays: self.decays.value()?,
            produced_nuclei: self.produced.value()?,
            population_roundoff_nuclei,
            substeps: self.substeps,
        })
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, FissionThermalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FissionThermalError::Arithmetic)
    }
}

#[cfg(test)]
mod tests;
