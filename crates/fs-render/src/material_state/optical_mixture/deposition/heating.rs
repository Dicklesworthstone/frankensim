//! Constituent heat recipients on the existing conservative thermal observer.
//!
//! Each positively absorbing `(optical cell, original component)` is assigned
//! exactly once. Different constituents may heat different native cells; any
//! number of optical contributions may feed one thermal cell. No duplicate
//! caloric law, integrator, event handling or optical clock is introduced.

use super::{ConstituentBeamDeposition, PositiveSum};
use crate::volumes::heating::{OpticalHeating, OpticalHeatingError, OpticalPowerBudget};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

/// Invalid constituent mapping or an original optical/thermal refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstituentHeatingError {
    DuplicateBinding {
        optical_cell: usize,
        component: usize,
    },
    MissingBinding {
        optical_cell: usize,
        component: usize,
    },
    UnknownBinding {
        optical_cell: usize,
        component: usize,
    },
    Heating(OpticalHeatingError),
}

impl From<OpticalHeatingError> for ConstituentHeatingError {
    fn from(error: OpticalHeatingError) -> Self {
        Self::Heating(error)
    }
}
impl From<DepositionError> for ConstituentHeatingError {
    fn from(error: DepositionError) -> Self {
        Self::Heating(OpticalHeatingError::Thermal(error))
    }
}
impl core::fmt::Display for ConstituentHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "constituent optical heating refused: {self:?}")
    }
}
impl core::error::Error for ConstituentHeatingError {}

impl<const N: usize> ConstituentBeamDeposition<N> {
    /// Bind retained constituent watts to active cells of the original solver.
    ///
    /// Bindings are `(optical_cell_index, component_index, thermal_cell_id)`.
    /// Every retained row needs exactly one binding. An omitted nonabsorber
    /// needs none; bindings to omitted rows refuse. The caller owns the
    /// physical mass/control-volume assignment and constituent thermal contacts.
    ///
    /// The returned native [`OpticalHeating`] preserves its existing fixed-
    /// snapshot semantics, clock, energy accounts, transactional observer
    /// composition and authored boundary sources. Recompute optics/deposition
    /// when the accepted constituent state or incident field changes.
    pub fn prepare_heating(
        &self,
        cell_bindings: &[(usize, usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<OpticalHeating<N>, ConstituentHeatingError> {
        checkpoint()?;
        if N == 0
            || thermal.cells().is_empty()
            || !power_tolerance_w.is_finite()
            || power_tolerance_w < 0.0
        {
            return Err(invalid(
                "active thermal cells, nonempty packet and finite nonnegative tolerance required",
            ));
        }
        if cell_bindings.len() > self.cells.len() {
            return Err(invalid(
                "more thermal bindings than absorbing constituent rows",
            ));
        }
        let mut bindings = BTreeMap::new();
        for &(optical_cell, component, recipient) in cell_bindings {
            checkpoint()?;
            if bindings
                .insert((optical_cell, component), recipient)
                .is_some()
            {
                return Err(ConstituentHeatingError::DuplicateBinding {
                    optical_cell,
                    component,
                });
            }
            if !thermal.cells().contains_key(&recipient) {
                return Err(OpticalHeatingError::InactiveCell(recipient).into());
            }
        }
        let mut source_sums = BTreeMap::<CellId, PositiveSum>::new();
        for row in &self.cells {
            checkpoint()?;
            let key = (row.optical_cell_index, row.component_index);
            let recipient =
                bindings
                    .remove(&key)
                    .ok_or(ConstituentHeatingError::MissingBinding {
                        optical_cell: key.0,
                        component: key.1,
                    })?;
            let source = source_sums.entry(recipient).or_default();
            for &power in &row.absorbed_power_w {
                source
                    .add(power)
                    .map_err(|_| invalid("constituent source sum is not representable"))?;
            }
        }
        if let Some((&(optical_cell, component), _)) = bindings.first_key_value() {
            return Err(ConstituentHeatingError::UnknownBinding {
                optical_cell,
                component,
            });
        }
        let mut sources_w = BTreeMap::new();
        let mut total_source = PositiveSum::default();
        for (recipient, sum) in source_sums {
            checkpoint()?;
            let power = sum
                .value()
                .map_err(|_| invalid("constituent source sum is not representable"))?;
            total_source
                .add(power)
                .map_err(|_| invalid("mapped thermal input is not representable"))?;
            sources_w.insert(recipient, power);
        }
        let thermal_input_w = total_source
            .value()
            .map_err(|_| invalid("mapped thermal input is not representable"))?;
        let mut absorbed_total = PositiveSum::default();
        let mut residual_w = [0.0; N];
        for lane in 0..N {
            checkpoint()?;
            let mut accounted = PositiveSum::default();
            for power in [
                self.beam.transmitted_power_w[lane],
                self.beam.absorbed_power_w[lane],
                self.beam.scattered_power_w[lane],
            ] {
                accounted
                    .add(power)
                    .map_err(|_| invalid("native beam balance is not representable"))?;
            }
            residual_w[lane] = self.beam.incident_power_w[lane]
                - accounted
                    .value()
                    .map_err(|_| invalid("native beam balance is not representable"))?;
            if residual_w[lane].abs() > power_tolerance_w
                || self.partition_residual_w[lane].abs() > power_tolerance_w
            {
                return Err(invalid(
                    "beam or constituent partition exceeds power tolerance",
                ));
            }
            absorbed_total
                .add(self.beam.absorbed_power_w[lane])
                .map_err(|_| invalid("absorbed packet power is not representable"))?;
        }
        let absorbed_total = absorbed_total
            .value()
            .map_err(|_| invalid("absorbed packet power is not representable"))?;
        if (thermal_input_w - absorbed_total).abs() > power_tolerance_w {
            return Err(invalid(
                "mapped constituent heat disagrees with true absorption",
            ));
        }
        let power = OpticalPowerBudget {
            medium_identity: self.beam.medium_identity,
            wavelength_bounds_nm: self.beam.wavelength_bounds_nm,
            wavelengths_nm: self.beam.wavelengths_nm,
            incident_w: self.beam.incident_power_w,
            transmitted_w: self.beam.transmitted_power_w,
            absorbed_w: self.beam.absorbed_power_w,
            scattered_w: self.beam.scattered_power_w,
            residual_w,
            thermal_input_w,
        };
        checkpoint()?;
        OpticalHeating::from_admitted_sources(power, sources_w, thermal).map_err(Into::into)
    }
}

fn invalid(reason: &'static str) -> ConstituentHeatingError {
    OpticalHeatingError::Invalid(reason).into()
}
