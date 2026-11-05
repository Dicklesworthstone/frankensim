//! Radioactive daughter inventories and integrated heat on the native thermal grid.
//!
//! Each cell retains an explicitly supplied decay inventory at the binding time.
//! Its existing prepared decay owner advances all daughters over the FINAL
//! accepted thermal interval. Only integrated locally deposited joules enter
//! the original heat-source path; radiation and nonrecoverable export remain
//! separate accounts. No irradiation history, deposition fraction, nuclear
//! material identity or change in the original thermal mass is inferred.
//!
//! Prepared maps are reused for repeated identical durations. Source evolution
//! is independent of temperature; thermal transport retains the original
//! uniform-cell, fixed-reference-geometry and timestep-refinement assumptions.

use std::{collections::BTreeMap, sync::Arc};

use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalStorageBasis,
};

use super::{DecayBudget, DecayError, DecayInventory, DecayStep, PreparedDecay};
use crate::EnergyChannels;

/// One existing decay inventory bound to one already active thermal cell.
/// A cell has at most one inventory in this observer; combine its inventory
/// before binding. Other independent heat-source owners can be composed through
/// `fs_thermal::observers::CoupledObservers`.
#[derive(Clone, Debug)]
pub struct DecayCellSource {
    /// Original active cell receiving the locally deposited energy.
    pub cell_id: CellId,
    /// Original thermal material identity expected at that cell.
    pub material_id: MaterialId,
    /// Explicit inventory at the current thermal clock; no earlier history is inferred.
    pub inventory: DecayInventory,
}

/// Nuclear and thermal-credit accounts for the last accepted advance.
/// These bounds concern numerical accounting, not nuclear-data uncertainty.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DecayHeatingAudit {
    /// Released energy in all three disjoint destinations.
    pub released: EnergyChannels,
    /// Actual `average_source_w * accepted_duration_s` credited to the thermal owner.
    pub thermal_credit_j: f64,
    /// Sum of absolute differences between integrated deposition and thermal credit.
    pub coupling_roundoff_j: f64,
    /// Noncancelling numerical energy defects reported by the decay owner.
    pub nuclear_energy_defect_j: f64,
    /// Noncancelling per-node inventory defects over all cells and intervals.
    pub population_roundoff_nuclei: f64,
    /// Integrated transformations, including decays of newly produced daughters.
    pub decays: f64,
    /// Number of common accepted thermal/decay intervals.
    pub substeps: usize,
}

/// Refusal preserving both original histories.
#[derive(Clone, Debug, PartialEq)]
pub enum DecayThermalError {
    /// Original thermal-domain, numerical, work or checkpoint refusal.
    Thermal(DepositionError),
    /// Original inventory, propagation or nuclear-account refusal.
    Decay(DecayError),
    /// Duplicate owner or inconsistent observer hook sequence.
    Invalid(&'static str),
    /// A source names a cell that has not been activated.
    InactiveCell(CellId),
    /// The source binding and active thermal cell name different materials.
    MaterialMismatch {
        /// Affected cell.
        cell_id: CellId,
        /// Authored binding identity.
        expected: MaterialId,
        /// Active thermal material identity.
        actual: MaterialId,
    },
    /// The original mass, caloric law, storage basis or physical geometry changed.
    BindingChanged(CellId),
    /// Independently advanced owners no longer represent the same time.
    ClockMismatch {
        /// Current thermal time [s].
        thermal_time_s: f64,
        /// Expected decay-inventory time [s].
        decay_time_s: f64,
    },
    /// Accumulated nuclear and thermal-credit defects exceed the advance allowance.
    EnergyBudget {
        /// Noncancelling numerical defects [J].
        defect_j: f64,
        /// Caller-supplied advance allowance [J].
        tolerance_j: f64,
    },
    /// Accumulated per-node inventory defects exceed the advance allowance.
    NucleiBudget {
        /// Noncancelling numerical population defects [expected nuclei].
        defect_nuclei: f64,
        /// Caller-supplied advance allowance [expected nuclei].
        tolerance_nuclei: f64,
    },
    /// A required finite quantity or counter cannot be represented.
    Arithmetic,
}

impl From<DepositionError> for DecayThermalError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<DecayError> for DecayThermalError {
    fn from(error: DecayError) -> Self {
        Self::Decay(error)
    }
}
impl std::fmt::Display for DecayThermalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "decay/thermal coupling refused: {self:?}")
    }
}
impl std::error::Error for DecayThermalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(error) => Some(error),
            Self::Decay(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
struct PreparedInterval {
    duration_s: f64,
    plan: Arc<PreparedDecay>,
}

#[derive(Clone, Debug)]
struct BoundSource {
    source: DecayCellSource,
    mass_kg: f64,
    geometry: CellGeometry,
    body_id: u64,
    segment_id: u64,
    caloric_identity: Option<[u8; 32]>,
    storage_basis: ThermalStorageBasis,
    prepared: Option<PreparedInterval>,
    last_step: Option<DecayStep>,
}

#[derive(Clone, Debug)]
struct PendingSource {
    cell_id: CellId,
    inventory: DecayInventory,
    report: DecayStep,
    prepared: PreparedInterval,
}

#[derive(Clone, Debug)]
struct PendingStep {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    sources: Vec<PendingSource>,
    account: Account,
    audit: DecayHeatingAudit,
}

/// Transactional decay-heat owner for `DepositionSimulation::advance_observed`.
///
/// Clone copies mutable inventories and shares only immutable models/prepared
/// maps. Use `advance_to` directly, or compose this observer with other physical
/// owners through `fs_thermal::observers::CoupledObservers`. Cell mass, geometry,
/// caloric identity and storage basis remain those of the original thermal law.
/// Nuclear transformations do not silently debit the cell's bulk thermal mass.
#[derive(Clone, Debug)]
pub struct DecayHeating {
    time_s: f64,
    sources: BTreeMap<CellId, BoundSource>,
    budget: DecayBudget,
    pending: Option<PendingStep>,
    account: Account,
    audit: DecayHeatingAudit,
}

impl DecayHeating {
    /// Bind explicit current inventories to the active cells and thermal clock.
    /// Birth events for all source cells must already have been processed.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        sources: Vec<DecayCellSource>,
        budget: DecayBudget,
    ) -> Result<Self, DecayThermalError> {
        if sources.is_empty() || sources.len() > thermal.cells().len() {
            return Err(DecayThermalError::Invalid(
                "source bindings need distinct active cells",
            ));
        }
        let mut bound = BTreeMap::new();
        for source in sources {
            let id = source.cell_id;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(DecayThermalError::InactiveCell(id))?;
            if cell.material_id != source.material_id {
                return Err(DecayThermalError::MaterialMismatch {
                    cell_id: id,
                    expected: source.material_id,
                    actual: cell.material_id,
                });
            }
            // Use the physical owner's admission even for a subsequent idle call.
            source.inventory.trial(0.0, budget, || Ok(()))?;
            let value = BoundSource {
                source,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                body_id: cell.body_id,
                segment_id: cell.segment_id,
                caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
                storage_basis: thermal.material_storage_basis(cell.material_id)?,
                prepared: None,
                last_step: None,
            };
            if bound.insert(id, value).is_some() {
                return Err(DecayThermalError::Invalid("duplicate decay source cell"));
            }
        }
        Ok(Self {
            time_s: thermal.time_s(),
            sources: bound,
            budget,
            pending: None,
            account: Account::default(),
            audit: DecayHeatingAudit::default(),
        })
    }

    /// Current common accepted clock [s].
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Accepted inventory for one physical cell.
    pub fn source(&self, cell: CellId) -> Option<&DecayCellSource> {
        self.sources.get(&cell).map(|bound| &bound.source)
    }

    /// Accepted inventories in deterministic physical-cell order.
    pub fn sources(&self) -> impl ExactSizeIterator<Item = &DecayCellSource> {
        self.sources.values().map(|bound| &bound.source)
    }

    /// Last accepted interval at one cell, absent before its first step.
    pub fn last_step(&self, cell: CellId) -> Option<&DecayStep> {
        self.sources
            .get(&cell)
            .and_then(|bound| bound.last_step.as_ref())
    }

    /// Accounts for the last accepted advance, excluding independently authored heat.
    pub fn audit(&self) -> DecayHeatingAudit {
        self.audit
    }

    /// Sum the original inventory owners' lifetime release accounts.
    pub fn cumulative_released(&self) -> Result<EnergyChannels, DecayThermalError> {
        let mut totals = [Sum::default(); 3];
        for source in self.sources() {
            let released = source.inventory.cumulative_released();
            totals[0].add(released.deposited_j)?;
            totals[1].add(released.radiation_j)?;
            totals[2].add(released.escaped_j)?;
        }
        channels(&totals)
    }

    /// Advance thermal state and all inventories atomically, including idle-call admission.
    pub fn advance_to<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'thermal>, Self, EnergyAudit), DecayThermalError> {
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        let (thermal, source, audit) =
            thermal.advance_observed(self, target_s, config, &mut checkpoint)?;
        source.validate(&thermal, source.time_s, &mut checkpoint)?;
        Ok((thermal, source, audit))
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        expected_time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DecayThermalError> {
        checkpoint()?;
        if thermal.time_s() != expected_time_s {
            return Err(DecayThermalError::ClockMismatch {
                thermal_time_s: thermal.time_s(),
                decay_time_s: expected_time_s,
            });
        }
        for (index, (&id, bound)) in self.sources.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(DecayThermalError::InactiveCell(id))?;
            if cell.material_id != bound.source.material_id {
                return Err(DecayThermalError::MaterialMismatch {
                    cell_id: id,
                    expected: bound.source.material_id,
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
                return Err(DecayThermalError::BindingChanged(id));
            }
        }
        Ok(())
    }
}

impl DepositionObserver for DecayHeating {
    type Error = DecayThermalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        if self.pending.is_some() {
            return Err(DecayThermalError::Invalid("uncommitted decay interval"));
        }
        self.account = Account::default();
        self.audit = DecayHeatingAudit::default();
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
        Ok(None)
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
            return Err(DecayThermalError::Invalid("invalid decay/thermal interval"));
        }
        let mut sources = Vec::with_capacity(self.sources.len());
        let mut watts = BTreeMap::new();
        let mut account = self.account.clone();
        for (&id, bound) in &self.sources {
            checkpoint()?;
            let mut checkpoint_error = None;
            let trial = (|| -> Result<_, DecayError> {
                let mut decay_checkpoint = || {
                    checkpoint().map_err(|error| {
                        checkpoint_error = Some(error);
                        DecayError::Cancelled
                    })
                };
                let prepared = match &bound.prepared {
                    Some(cached) if cached.duration_s == duration_s => cached.clone(),
                    _ => PreparedInterval {
                        duration_s,
                        plan: Arc::new(bound.source.inventory.model().prepare(
                            duration_s,
                            self.budget.propagation,
                            &mut decay_checkpoint,
                        )?),
                    },
                };
                let (inventory, report) = bound.source.inventory.trial_prepared(
                    &prepared.plan,
                    self.budget,
                    &mut decay_checkpoint,
                )?;
                Ok((inventory, report, prepared))
            })();
            if let Some(error) = checkpoint_error {
                return Err(error.into());
            }
            let (inventory, report, prepared) = trial?;
            let power = finite(report.released.deposited_j / duration_s)?;
            let credit = finite(power * duration_s)?;
            account.record(&report, credit)?;
            watts.insert(id, power);
            sources.push(PendingSource {
                cell_id: id,
                inventory,
                report,
                prepared,
            });
        }
        account.substeps = account
            .substeps
            .checked_add(1)
            .ok_or(DecayThermalError::Arithmetic)?;
        let audit = account.snapshot()?;
        let defect = finite(audit.coupling_roundoff_j + audit.nuclear_energy_defect_j)?;
        if defect > self.budget.energy_tolerance_j {
            return Err(DecayThermalError::EnergyBudget {
                defect_j: defect,
                tolerance_j: self.budget.energy_tolerance_j,
            });
        }
        if audit.population_roundoff_nuclei > self.budget.nuclei_tolerance {
            return Err(DecayThermalError::NucleiBudget {
                defect_nuclei: audit.population_roundoff_nuclei,
                tolerance_nuclei: self.budget.nuclei_tolerance,
            });
        }
        checkpoint()?;
        self.pending = Some(PendingStep {
            start_s,
            duration_s,
            end_s,
            sources,
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
        let pending = self.pending.as_ref().ok_or(DecayThermalError::Invalid(
            "thermal step has no decay candidate",
        ))?;
        if pending.start_s != start_s || pending.duration_s != duration_s {
            return Err(DecayThermalError::Invalid(
                "thermal step changed the decay interval",
            ));
        }
        self.validate(thermal, pending.end_s, checkpoint)?;
        checkpoint()?;
        // Publication is infallible after all common-state and budget checks.
        let pending = self.pending.take().expect("validated pending interval");
        for candidate in pending.sources {
            let bound = self
                .sources
                .get_mut(&candidate.cell_id)
                .expect("candidate retains source identities");
            bound.source.inventory = candidate.inventory;
            bound.prepared = Some(candidate.prepared);
            bound.last_step = Some(candidate.report);
        }
        self.time_s = pending.end_s;
        self.account = pending.account;
        self.audit = pending.audit;
        Ok(())
    }
}

// Small interval releases and noncancelling defects retain a compensated account.
#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DecayThermalError> {
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

    fn value(&self) -> Result<f64, DecayThermalError> {
        finite(self.sum + self.correction)
    }
}

#[derive(Clone, Debug, Default)]
struct Account {
    released: [Sum; 3],
    credit: Sum,
    coupling_roundoff: Sum,
    nuclear_defect: Sum,
    population_roundoff: Sum,
    decays: Sum,
    substeps: usize,
}

impl Account {
    fn record(&mut self, report: &DecayStep, credit: f64) -> Result<(), DecayThermalError> {
        self.released[0].add(report.released.deposited_j)?;
        self.released[1].add(report.released.radiation_j)?;
        self.released[2].add(report.released.escaped_j)?;
        self.credit.add(credit)?;
        self.coupling_roundoff
            .add((credit - report.released.deposited_j).abs())?;
        self.nuclear_defect.add(report.energy_roundoff_j)?;
        self.population_roundoff
            .add(report.population_roundoff_nuclei)?;
        for &count in &report.decays {
            self.decays.add(count)?;
        }
        Ok(())
    }

    fn snapshot(&self) -> Result<DecayHeatingAudit, DecayThermalError> {
        Ok(DecayHeatingAudit {
            released: channels(&self.released)?,
            thermal_credit_j: self.credit.value()?,
            coupling_roundoff_j: self.coupling_roundoff.value()?,
            nuclear_energy_defect_j: self.nuclear_defect.value()?,
            population_roundoff_nuclei: self.population_roundoff.value()?,
            decays: self.decays.value()?,
            substeps: self.substeps,
        })
    }
}

fn channels(values: &[Sum; 3]) -> Result<EnergyChannels, DecayThermalError> {
    let result = EnergyChannels {
        deposited_j: values[0].value()?,
        radiation_j: values[1].value()?,
        escaped_j: values[2].value()?,
    };
    finite(result.total_j())?;
    Ok(result)
}

fn finite(value: f64) -> Result<f64, DecayThermalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DecayThermalError::Arithmetic)
    }
}

#[cfg(test)]
mod tests;
