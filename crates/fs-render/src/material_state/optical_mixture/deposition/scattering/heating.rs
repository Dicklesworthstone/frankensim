//! Constituent-resolved multiple-scattering heat on the native thermal clock.
//!
//! This bridge consumes one already sampled optical report. It does not rerun
//! photon histories, reinterpret escaped radiation, or create a second thermal
//! integrator. Every retained constituent row is assigned exactly once.

use super::super::PositiveSum;
use super::ConstituentScatteringDeposition;
use crate::material_state::optical_mixture::deposition::heating::ConstituentHeatingError;
use crate::volumes::heating::OpticalHeatingError;
use crate::volumes::heating::scattering::{
    ScatteringHeating, ScatteringHeatingTolerance,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

impl<const N: usize> ConstituentScatteringDeposition<N> {
    /// Consume this immutable sampled report and assign each constituent row to
    /// one active native thermal cell.
    ///
    /// Bindings are `(optical_cell_index, component_index, thermal_cell_id)`.
    /// Several optical/component rows may feed the same thermal cell. Means are
    /// summed exactly as heat rates; empirical per-row standard errors are summed
    /// conservatively because constituent shares from one native cell are
    /// correlated. If any sampled absorbing lane was unobserved, its row carries
    /// no numerical error estimate and heating refuses rather than treating that
    /// missing evidence as zero uncertainty.
    ///
    /// This method consumes the report so the large native history summary moves
    /// into `ScatteringHeating` without another clone.
    pub fn into_heating(
        self,
        cell_bindings: &[(usize, usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        tolerance: ScatteringHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ScatteringHeating<N>, ConstituentHeatingError> {
        checkpoint()?;
        if N == 0 || thermal.cells().is_empty() || cell_bindings.len() > self.cells.len() {
            return Err(invalid(
                "active thermal cells, nonempty packet and no excess constituent bindings required",
            ));
        }
        for residual in self.partition_residual_w {
            if !residual.is_finite() || residual.abs() > tolerance.power_roundoff_w {
                return Err(invalid(
                    "constituent partition residual exceeds the scattering heating allowance",
                ));
            }
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

        let mut sources = BTreeMap::<CellId, PositiveSum>::new();
        let mut errors = BTreeMap::<CellId, PositiveSum>::new();
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
            let source = sources.entry(recipient).or_default();
            let error = errors.entry(recipient).or_default();
            for lane in 0..N {
                source
                    .add(row.absorbed_power_w[lane])
                    .map_err(|_| invalid("constituent scattered heat is unrepresentable"))?;
                let standard_error = row.standard_error_w[lane].ok_or_else(|| {
                    invalid(
                        "unobserved constituent absorption cannot be admitted as a thermal source",
                    )
                })?;
                error
                    .add(standard_error)
                    .map_err(|_| invalid("constituent sampling error is unrepresentable"))?;
            }
        }
        if let Some((&(optical_cell, component), _)) = bindings.first_key_value() {
            return Err(ConstituentHeatingError::UnknownBinding {
                optical_cell,
                component,
            });
        }

        let mut sources_w = BTreeMap::new();
        let mut receiving_error_w = BTreeMap::new();
        for (recipient, source) in sources {
            checkpoint()?;
            let power = source
                .value()
                .map_err(|_| invalid("constituent scattered heat is unrepresentable"))?;
            let error = errors
                .remove(&recipient)
                .ok_or_else(|| invalid("missing constituent sampling-error account"))?
                .value()
                .map_err(|_| invalid("constituent sampling error is unrepresentable"))?;
            if power > 0.0 {
                sources_w.insert(recipient, power);
                receiving_error_w.insert(recipient, error);
            } else if error != 0.0 {
                return Err(invalid(
                    "zero constituent heat cannot carry positive receiving uncertainty",
                ));
            }
        }
        if !errors.is_empty() {
            return Err(invalid("orphan constituent sampling-error account"));
        }
        checkpoint()?;
        ScatteringHeating::from_constituent_sources(
            self.scattering,
            sources_w,
            receiving_error_w,
            thermal,
            tolerance,
            checkpoint,
        )
        .map_err(Into::into)
    }
}

fn invalid(reason: &'static str) -> ConstituentHeatingError {
    OpticalHeatingError::Invalid(reason).into()
}
