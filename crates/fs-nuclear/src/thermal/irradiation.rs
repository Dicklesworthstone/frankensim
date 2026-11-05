//! Multiple explicit parent-isotope inventories on the original thermal grid.

use std::sync::Arc;

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating, Sum, finite};
use crate::irradiation::{GroupFlux, IrradiationBudget, IrradiationState, IrradiationStep};
use crate::{EnergyChannels, FissionError};
use fs_thermal::enthalpy::{CellId, MaterialId};

/// An absolute-time change in the incident multigroup neutron flux.
#[derive(Clone, Debug)]
pub struct NeutronFluxChange {
    /// Strictly increasing positive absolute time [s].
    pub time_s: f64,
    /// Incident group-integrated scalar flux [neutrons / m² / s].
    pub flux: GroupFlux,
}

/// Immutable piecewise-constant neutron spectrum on one explicitly named grid.
#[derive(Clone, Debug)]
pub struct NeutronFluxSchedule {
    initial: GroupFlux,
    changes: Arc<[NeutronFluxChange]>,
}

impl NeutronFluxSchedule {
    /// Admit times and the exact same energy-group identity and edges throughout.
    pub fn new(
        initial: GroupFlux,
        changes: Vec<NeutronFluxChange>,
    ) -> Result<Self, FissionThermalError> {
        let mut previous = 0.0;
        for change in &changes {
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(FissionThermalError::Invalid(
                    "neutron flux changes need strictly increasing positive times",
                ));
            }
            if change.flux.groups() != initial.groups() {
                return Err(FissionThermalError::Invalid(
                    "neutron flux schedule changes its energy-group identity or edges",
                ));
            }
            previous = change.time_s;
        }
        Ok(Self {
            initial,
            changes: changes.into(),
        })
    }

    /// Right-continuous incident spectrum at a finite nonnegative time [s].
    pub fn flux_at(&self, time_s: f64) -> Result<&GroupFlux, FissionThermalError> {
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(FissionThermalError::Invalid(
                "invalid neutron flux schedule time",
            ));
        }
        let index = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        Ok(if index == 0 {
            &self.initial
        } else {
            &self.changes[index - 1].flux
        })
    }

    /// Authored spectrum changes in increasing absolute-time order.
    pub fn changes(&self) -> &[NeutronFluxChange] {
        &self.changes
    }

    fn next_after(&self, time_s: f64) -> Option<f64> {
        let index = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        self.changes.get(index).map(|change| change.time_s)
    }
}

/// Work limits for one cell's isotope mixture and its common thermal advance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IrradiationHeatingBudget {
    /// Maximum distinct tracked parent isotopes in one thermal cell.
    pub max_isotopes_per_cell: usize,
    /// Original group, fuel and numerical energy allowances.
    /// Its energy allowance also bounds the aggregate across cells and substeps.
    pub source: IrradiationBudget,
}

/// Explicit parent-isotope inventories sharing one incident spectrum and cell.
///
/// The neutron flux is prescribed. Each isotope competes between fission and
/// capture; daughters and fission products are outside this parent-only model.
/// Their material mass is retained by the original caloric owner. A supplied
/// temperature domain is checked at both ends of each thermal interval; the
/// cross sections remain the explicitly supplied frozen coefficients.
#[derive(Clone, Debug)]
pub struct IrradiationCellSource {
    /// Existing active thermal receiver.
    pub cell_id: CellId,
    /// Original thermal material identity.
    pub material_id: MaterialId,
    /// One explicit inventory per distinct parent nuclide, in deterministic order.
    pub inventories: Vec<IrradiationState>,
    /// Shared incident neutron spectrum history.
    pub fluxes: NeutronFluxSchedule,
}

/// One cell interval, retaining each isotope's own accepted reaction account.
#[derive(Clone, Debug, PartialEq)]
pub struct IrradiationCellStep {
    /// Reports in the same order as the bound parent-isotope inventories.
    pub isotopes: Vec<IrradiationStep>,
    /// Sum of physically separate release destinations [J].
    pub released: EnergyChannels,
    /// Total accepted fission events.
    pub fissions: f64,
    /// Total accepted capture events, never credited as fissions.
    pub captures: f64,
    /// Sum of absolute source energy-account defects [J].
    pub energy_roundoff_j: f64,
}

/// The shared conservative thermal owner specialized for neutron irradiation.
pub type IrradiationHeating = NuclearHeating<IrradiationCellSource>;

impl sealed::Sealed for IrradiationCellSource {}
impl NuclearHeatSource for IrradiationCellSource {
    type Budget = IrradiationHeatingBudget;
    type Step = IrradiationCellStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.source.source.energy_tolerance_j
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.fluxes.next_after(time_s)
    }

    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        if self.inventories.is_empty() || self.inventories.len() > budget.max_isotopes_per_cell {
            return Err(FissionThermalError::Invalid(
                "invalid number of parent-isotope inventories",
            ));
        }
        for (index, inventory) in self.inventories.iter().enumerate() {
            if self.inventories[..index]
                .iter()
                .any(|prior| prior.model().nuclide() == inventory.model().nuclide())
            {
                return Err(FissionThermalError::Invalid(
                    "duplicate parent nuclide in one cell",
                ));
            }
            inventory.trial(
                0.0,
                self.fluxes.flux_at(time_s)?,
                temperature_k,
                budget.source,
                || Ok(()),
            )?;
        }
        Ok(())
    }

    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        for inventory in &self.inventories {
            let model = inventory.model();
            if !temperature_k.is_finite()
                || temperature_k < model.temperature_min_k()
                || temperature_k > model.temperature_max_k()
            {
                return Err(FissionError::InvalidInput(
                    "thermal temperature is outside nuclear-data applicability",
                )
                .into());
            }
        }
        Ok(())
    }

    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        self.trial_flux(
            self.fluxes.flux_at(time_s)?,
            duration_s,
            temperature_k,
            budget,
            checkpoint,
        )
    }

    fn released(step: &Self::Step) -> EnergyChannels {
        step.released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.energy_roundoff_j
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.fissions
    }
    fn captures(step: &Self::Step) -> f64 {
        step.captures
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        let mut totals = [Sum::default(); 3];
        for inventory in &self.inventories {
            let channels = inventory.cumulative_released();
            totals[0].add(channels.deposited_j)?;
            totals[1].add(channels.radiation_j)?;
            totals[2].add(channels.escaped_j)?;
        }
        let total = EnergyChannels {
            deposited_j: totals[0].value()?,
            radiation_j: totals[1].value()?,
            escaped_j: totals[2].value()?,
        };
        finite(total.total_j())?;
        Ok(total)
    }
}

impl IrradiationCellSource {
    pub(super) fn trial_flux(
        &self,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: IrradiationHeatingBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, IrradiationCellStep), FissionThermalError> {
        let mut inventories = Vec::with_capacity(self.inventories.len());
        let mut isotopes = Vec::with_capacity(self.inventories.len());
        let mut released = [Sum::default(); 3];
        let mut fissions = Sum::default();
        let mut captures = Sum::default();
        let mut defects = Sum::default();
        for inventory in &self.inventories {
            checkpoint()?;
            let (state, report) = inventory.trial(
                duration_s,
                flux,
                temperature_k,
                budget.source,
                &mut *checkpoint,
            )?;
            released[0].add(report.released.deposited_j)?;
            released[1].add(report.released.radiation_j)?;
            released[2].add(report.released.escaped_j)?;
            fissions.add(report.fissions)?;
            captures.add(report.captures)?;
            defects.add(report.energy_roundoff_j)?;
            inventories.push(state);
            isotopes.push(report);
        }
        // A cell emits one scalar heat rate. Charge the final export from each
        // compensated isotope sum so a small isotope cannot vanish for free.
        for channel in &released {
            defects.add(channel.export_roundoff()?)?;
        }
        let report = IrradiationCellStep {
            isotopes,
            released: EnergyChannels {
                deposited_j: released[0].value()?,
                radiation_j: released[1].value()?,
                escaped_j: released[2].value()?,
            },
            fissions: fissions.value()?,
            captures: captures.value()?,
            energy_roundoff_j: defects.value()?,
        };
        finite(report.released.total_j())?;
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                inventories,
                fluxes: self.fluxes.clone(),
            },
            report,
        ))
    }
}
