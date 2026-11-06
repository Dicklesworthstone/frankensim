//! Explicit global-loss partition; no mechanical or caloric constitutive model.

use fs_math::eft::two_sum;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use std::collections::BTreeMap;

/// Same-time original thermal response and the declared unconverted loss.
#[derive(Clone, Copy, Debug)]
pub struct ImpactHeatReceipt {
    pub time_s: f64,
    /// The original caloric owner's audit, with heat external to that subsystem.
    pub thermal: EnergyAudit,
    pub unallocated_energy_j: f64,
    /// Rounded allocated heat + declared remainder - native mechanical loss.
    /// This numerical defect is gated separately and never credited as heat.
    pub partition_energy_residual_j: f64,
}

pub(super) struct ImpactHeatPartition<'a> {
    fractions: &'a BTreeMap<CellId, f64>,
    remaining_fraction: [f64; 2],
    cell_tolerance: f64,
    tolerance: f64,
}

impl<'a> ImpactHeatPartition<'a> {
    pub(super) fn new(
        thermal: &DepositionSimulation<'_>,
        fractions: &'a BTreeMap<CellId, f64>,
        cell_tolerance: f64,
        tolerance: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !nonnegative(cell_tolerance)
            || !nonnegative(tolerance)
            || fractions.len() > thermal.cells().len()
        {
            return Err(DepositionError::InvalidInput(
                "bounded impact heat partition and finite tolerances required",
            ));
        }
        let mut remaining = [1.0, 0.0];
        for (index, (&id, &fraction)) in fractions.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if !nonnegative(fraction) || fraction > 1.0 || !thermal.cells().contains_key(&id) {
                return Err(DepositionError::InvalidInput(
                    "impact heat needs active cells and nonnegative fractions",
                ));
            }
            add(&mut remaining, -fraction);
            if remaining[0] < 0.0 || (remaining[0] == 0.0 && remaining[1] < 0.0) {
                return Err(DepositionError::InvalidInput(
                    "impact heat fractions exceed the native loss",
                ));
            }
        }
        checkpoint()?;
        Ok(Self {
            fractions,
            remaining_fraction: remaining,
            cell_tolerance,
            tolerance,
        })
    }

    pub(super) fn deposit<'m>(
        self,
        thermal: &DepositionSimulation<'m>,
        loss_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, ImpactHeatReceipt), DepositionError> {
        checkpoint()?;
        if !nonnegative(loss_j) {
            return Err(DepositionError::InvalidInput(
                "finite nonnegative native impact loss required",
            ));
        }
        let mut heat = BTreeMap::new();
        let mut balance = [-loss_j, 0.0];
        for (index, (&id, &fraction)) in self.fractions.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let value = loss_j * fraction;
            if !value.is_finite() || (value == 0.0 && loss_j != 0.0 && fraction != 0.0) {
                return Err(DepositionError::Arithmetic);
            }
            if value != 0.0 {
                heat.insert(id, value);
            }
            add(&mut balance, value);
        }
        let [high, low] = self.remaining_fraction;
        let unallocated = loss_j.mul_add(high, loss_j * low);
        if !nonnegative(unallocated)
            || (unallocated == 0.0 && loss_j != 0.0 && (high != 0.0 || low != 0.0))
        {
            return Err(DepositionError::Arithmetic);
        }
        add(&mut balance, unallocated);
        let residual = balance[0] + balance[1];
        if !residual.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        if residual.abs() > self.tolerance {
            return Err(DepositionError::EnergyBalance(residual));
        }
        let (next, audit) =
            thermal.apply_heat_j(&heat, self.cell_tolerance, self.tolerance, &mut checkpoint)?;
        checkpoint()?;
        Ok((
            next,
            ImpactHeatReceipt {
                time_s: thermal.time_s(),
                thermal: audit,
                unallocated_energy_j: unallocated,
                partition_energy_residual_j: residual,
            },
        ))
    }
}

// Original EFT primitive retains sub-ulp fractions and small energy defects.
fn add(value: &mut [f64; 2], term: f64) {
    let (sum, error) = two_sum(value[0], term);
    let (high, low) = two_sum(sum, value[1] + error);
    *value = [high, low];
}
pub(super) fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = [0.0; 2];
    for value in values {
        add(&mut sum, value);
    }
    sum[0] + sum[1]
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

#[cfg(test)]
#[path = "partition/tests.rs"]
mod tests;
