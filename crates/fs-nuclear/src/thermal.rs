//! Prescribed fission and decay heat on the existing conservative thermal grid.
//!
//! Only the nuclear source's locally deposited energy is credited to its bound
//! cell. Transportable radiation and nonrecoverable export remain separate in
//! the source account. The original thermal law, contact exchange, boundaries,
//! enthalpy/internal-energy basis and thermal audit retain their existing owners.
//!
//! Rates are right-continuous, piecewise constant histories, not a neutron or
//! temperature-feedback calculation. Their transitions clip the native thermal
//! step exactly. Nuclear reservoirs advance analytically over the FINAL accepted
//! interval, including any fuel exhaustion within it. Thermal transport retains
//! its original explicit timestep and uniform-cell limitations.

use std::{collections::BTreeMap, sync::Arc};

use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, Event, MaterialId, StepConfig, ThermalStorageBasis,
};

use crate::{EnergyChannels, FissionBudget, FissionError, FissionState, FissionStep};

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
        Ok(Self { initial_fissions_per_s, changes: changes.into() })
    }

    /// Prescribed rate at the given nonnegative finite absolute time [s].
    pub fn rate_at(&self, time_s: f64) -> Result<f64, FissionThermalError> {
        if !nonnegative(time_s) {
            return Err(FissionThermalError::Invalid("invalid fission schedule time"));
        }
        let index = self.changes.partition_point(|change| change.time_s <= time_s);
        Ok(if index == 0 {
            self.initial_fissions_per_s
        } else {
            self.changes[index - 1].fissions_per_s
        })
    }

    /// Authored changes; no mutable schedule can escape into a running source.
    pub fn changes(&self) -> &[FissionRateChange] { &self.changes }

    fn next_after(&self, time_s: f64) -> Option<f64> {
        let index = self.changes.partition_point(|change| change.time_s <= time_s);
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
    /// Finite inputs produced nonfinite arithmetic or an exhausted counter.
    Arithmetic,
}

impl From<DepositionError> for FissionThermalError {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl From<FissionError> for FissionThermalError {
    fn from(error: FissionError) -> Self { Self::Fission(error) }
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
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
struct BoundSource {
    source: FissionCellSource,
    mass_kg: f64,
    geometry: CellGeometry,
    body_id: u64,
    segment_id: u64,
    caloric_identity: Option<[u8; 32]>,
    storage_basis: ThermalStorageBasis,
    last_step: Option<FissionStep>,
}

#[derive(Clone, Debug)]
struct PendingStep {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    states: Vec<(CellId, FissionState, FissionStep)>,
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
pub struct FissionHeating {
    time_s: f64,
    sources: BTreeMap<CellId, BoundSource>,
    budget: FissionBudget,
    pending: Option<PendingStep>,
    account: Account,
    audit: FissionHeatingAudit,
}

impl FissionHeating {
    /// Bind current nuclear inventories to the current thermal clock/materials.
    /// Birth events for all source cells must already have been processed.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        sources: Vec<FissionCellSource>,
        budget: FissionBudget,
    ) -> Result<Self, FissionThermalError> {
        if sources.is_empty() || sources.len() > thermal.cells().len() {
            return Err(FissionThermalError::Invalid("source bindings need distinct active cells"));
        }
        let mut bound = BTreeMap::new();
        for source in sources {
            let id = source.cell_id;
            let cell = thermal.cells().get(&id).ok_or(FissionThermalError::InactiveCell(id))?;
            if cell.material_id != source.material_id {
                return Err(FissionThermalError::MaterialMismatch {
                    cell_id: id, expected: source.material_id, actual: cell.material_id,
                });
            }
            // The same kernel admits its work/tolerance budget, even when the
            // first requested thermal advance will have zero duration.
            source.state.trial(0.0, 0.0, budget, || Ok(()))?;
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
                return Err(FissionThermalError::Invalid("duplicate fission source cell"));
            }
        }
        Ok(Self {
            time_s: thermal.time_s(), sources: bound, budget, pending: None,
            account: Account::default(), audit: FissionHeatingAudit::default(),
        })
    }

    /// Current accepted source/thermal clock [s].
    pub fn time_s(&self) -> f64 { self.time_s }
    /// Accepted source inventory and immutable schedule for one cell.
    pub fn source(&self, cell: CellId) -> Option<&FissionCellSource> {
        self.sources.get(&cell).map(|bound| &bound.source)
    }
    /// Accepted source states in deterministic physical-cell order.
    pub fn sources(&self) -> impl ExactSizeIterator<Item = &FissionCellSource> {
        self.sources.values().map(|bound| &bound.source)
    }
    /// Last accepted nuclear interval at a cell, absent before the first step.
    pub fn last_step(&self, cell: CellId) -> Option<&FissionStep> {
        self.sources.get(&cell).and_then(|bound| bound.last_step.as_ref())
    }
    /// Source/coupling account for the last accepted advance, excluding other heat.
    pub fn audit(&self) -> FissionHeatingAudit { self.audit }

    /// Sum the original source owners' lifetime release accounts.
    pub fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        let mut totals = [Sum::default(); 3];
        for source in self.sources() {
            let released = source.state.cumulative_released();
            totals[0].add(released.deposited_j)?;
            totals[1].add(released.radiation_j)?;
            totals[2].add(released.escaped_j)?;
        }
        Ok(EnergyChannels {
            deposited_j: totals[0].value()?, radiation_j: totals[1].value()?,
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
        let (thermal, source, audit) = thermal.advance_observed(self, target_s, config, &mut checkpoint)?;
        source.validate(&thermal, source.time_s, &mut checkpoint)?;
        Ok((thermal, source, audit))
    }

    fn validate(
        &self, thermal: &DepositionSimulation<'_>, expected_time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        checkpoint()?;
        if thermal.time_s() != expected_time_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s: thermal.time_s(), fission_time_s: expected_time_s,
            });
        }
        for (index, (&id, bound)) in self.sources.iter().enumerate() {
            if index % 64 == 0 { checkpoint()?; }
            let cell = thermal.cells().get(&id).ok_or(FissionThermalError::InactiveCell(id))?;
            if cell.material_id != bound.source.material_id {
                return Err(FissionThermalError::MaterialMismatch {
                    cell_id: id, expected: bound.source.material_id, actual: cell.material_id,
                });
            }
            if cell.mass_kg != bound.mass_kg || cell.geometry != bound.geometry
                || cell.body_id != bound.body_id || cell.segment_id != bound.segment_id
                || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
                || thermal.material_storage_basis(cell.material_id)? != bound.storage_basis
            {
                return Err(FissionThermalError::BindingChanged(id));
            }
        }
        Ok(())
    }
}

impl DepositionObserver for FissionHeating {
    type Error = FissionThermalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        if self.pending.is_some() {
            return Err(FissionThermalError::Invalid("uncommitted nuclear interval"));
        }
        self.account = Account::default();
        self.audit = FissionHeatingAudit::default();
        Ok(())
    }

    fn requires_step_interval(&self) -> bool { true }

    fn next_transition_time_s(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.validate(thermal, self.time_s, checkpoint)?;
        Ok(self.sources().filter_map(|source| source.rates.next_after(self.time_s))
            .min_by(f64::total_cmp))
    }

    fn step_sources(
        &mut self, thermal: &DepositionSimulation<'_>, start_s: f64, duration_s: f64,
        end_s: f64, checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.validate(thermal, self.time_s, checkpoint)?;
        if self.pending.is_some() || start_s != self.time_s || !duration_s.is_finite()
            || duration_s <= 0.0 || !end_s.is_finite() || end_s <= start_s
        {
            return Err(FissionThermalError::Invalid("invalid nuclear/thermal interval"));
        }
        let mut states = Vec::with_capacity(self.sources.len());
        let mut watts = BTreeMap::new();
        let mut account = self.account.clone();
        for (&id, bound) in &self.sources {
            checkpoint()?;
            if bound.source.rates.next_after(start_s).is_some_and(|time| time < end_s) {
                return Err(FissionThermalError::Invalid("nuclear interval crosses a rate change"));
            }
            let rate = bound.source.rates.rate_at(start_s)?;
            let mut checkpoint_error = None;
            let trial = bound.source.state.trial(duration_s, rate, self.budget, || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    FissionError::Cancelled
                })
            });
            if let Some(error) = checkpoint_error { return Err(error.into()); }
            let (state, report) = trial?;
            let power = finite(report.released.deposited_j / duration_s)?;
            let credit = finite(power * duration_s)?;
            account.record(&report, credit)?;
            watts.insert(id, power);
            states.push((id, state, report));
        }
        account.substeps = account.substeps.checked_add(1).ok_or(FissionThermalError::Arithmetic)?;
        let audit = account.snapshot()?;
        let defect = finite(audit.coupling_roundoff_j + audit.nuclear_energy_defect_j)?;
        if defect > self.budget.energy_tolerance_j {
            return Err(FissionThermalError::EnergyBudget {
                defect_j: defect, tolerance_j: self.budget.energy_tolerance_j,
            });
        }
        checkpoint()?;
        self.pending = Some(PendingStep { start_s, duration_s, end_s, states, account, audit });
        Ok(watts)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.validate(thermal, self.time_s, &mut || Ok(()))
    }

    fn step(
        &mut self, start_s: f64, duration_s: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let pending = self.pending.as_ref()
            .ok_or(FissionThermalError::Invalid("thermal step has no nuclear candidate"))?;
        if pending.start_s != start_s || pending.duration_s != duration_s {
            return Err(FissionThermalError::Invalid("thermal step changed the nuclear interval"));
        }
        self.validate(thermal, pending.end_s, checkpoint)?;
        checkpoint()?;
        // All fallible work precedes publication, even for direct hook callers.
        let pending = self.pending.take().expect("validated pending interval");
        for (id, state, report) in pending.states {
            let source = self.sources.get_mut(&id).expect("candidate retains source identities");
            source.source.state = state;
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
struct Sum { sum: f64, correction: f64 }
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
    fn value(&self) -> Result<f64, FissionThermalError> { finite(self.sum + self.correction) }
}

#[derive(Clone, Debug, Default)]
struct Account {
    released: [Sum; 3],
    credit: Sum,
    coupling_roundoff: Sum,
    nuclear_defect: Sum,
    fissions: Sum,
    substeps: usize,
}
impl Account {
    fn record(&mut self, report: &FissionStep, credit: f64) -> Result<(), FissionThermalError> {
        self.released[0].add(report.released.deposited_j)?;
        self.released[1].add(report.released.radiation_j)?;
        self.released[2].add(report.released.escaped_j)?;
        self.credit.add(credit)?;
        self.coupling_roundoff.add((credit - report.released.deposited_j).abs())?;
        // The core's noncancelling roundoff already includes its residual.
        self.nuclear_defect.add(report.energy_roundoff_j)?;
        self.fissions.add(report.fissions)?;
        Ok(())
    }
    fn snapshot(&self) -> Result<FissionHeatingAudit, FissionThermalError> {
        Ok(FissionHeatingAudit {
            released: EnergyChannels {
                deposited_j: self.released[0].value()?, radiation_j: self.released[1].value()?,
                escaped_j: self.released[2].value()?,
            },
            thermal_credit_j: self.credit.value()?,
            coupling_roundoff_j: self.coupling_roundoff.value()?,
            nuclear_energy_defect_j: self.nuclear_defect.value()?,
            fissions: self.fissions.value()?, substeps: self.substeps,
        })
    }
}

fn nonnegative(value: f64) -> bool { value.is_finite() && value >= 0.0 }
fn finite(value: f64) -> Result<f64, FissionThermalError> {
    if value.is_finite() { Ok(value) } else { Err(FissionThermalError::Arithmetic) }
}

#[cfg(test)]
mod tests;
