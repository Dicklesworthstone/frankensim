//! Secondary luminescent-light absorption by explicitly passive receivers.
//!
//! The prepared directed photon operator is shared across native thermal
//! clones and reused for every interval. Integrated emitted energy is routed
//! once: resolved absorption heats its original receiver, escape leaves the
//! grid, and an admitted numerical tail remains unresolved radiation. The
//! receiving absorption channel is declared nonluminescent; it does not
//! silently reuse the pump source outside its excitation domain. Elastic
//! multiple scattering is already included by the original transport owner.

use super::super::exchange::RadiativeThermalCell;
use super::{CheckedSum, OpticalHeatingError, finite, invalid, nonnegative};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::ScatteringGrid;
use crate::volumes::transport::scattering::luminescence::LuminescenceTransferOperator;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Explicit admission of numerical transport diagnostics. A fixed-sample
/// standard error is not a confidence interval or a bound on unseen paths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuminescenceAbsorptionTolerance {
    /// Absolute allowance for each dimensionless probability-row identity.
    pub fraction_roundoff: f64,
    /// Maximum empirical standard error of any retained fraction tally.
    pub max_fraction_standard_error: f64,
    /// Maximum cumulative unresolved radiant energy, in joules. Zero requires
    /// complete histories. A permitted tail is never assigned to heat/escape.
    pub max_unresolved_energy_j: f64,
}

#[derive(Debug)]
struct Model<const M: usize> {
    operator: LuminescenceTransferOperator<M>,
    bindings: BTreeMap<usize, CellId>,
    /// Receiver IDs in exactly the original row/link order, prepared once.
    row_receivers: Vec<Vec<CellId>>,
    temperature_domains: Vec<(CellId, f64, f64)>,
    tolerance: LuminescenceAbsorptionTolerance,
}

/// An open, index-one, fixed optical environment whose emission-band
/// absorption is explicitly passive. This can model host/substrate absorption
/// and nonradiative spectral channels in the emitting material. It is not a
/// declaration that all absorption by a fluorophore is passive.
#[derive(Clone, Debug)]
pub struct AbsorbingLuminescenceEnvironment<const M: usize> {
    model: Arc<Model<M>>,
}

pub(super) struct RoutedEnergy<const M: usize> {
    pub escaped_j: [f64; M],
    pub unresolved_j: [f64; M],
    pub heat_j: f64,
    pub absolute_roundoff_j: f64,
}

impl<const M: usize> AbsorbingLuminescenceEnvironment<M> {
    /// Bind all original positive-density cells, including pure scatterers,
    /// to native thermal cells and explicitly authored optical temperature
    /// domains. The operator must describe this exact absorption, scattering
    /// and phase snapshot; unrelated replacement coefficient fields refuse.
    /// Photon flight is quasistatic relative to the resolved thermal clock.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        medium: &ScatteringGrid<'_>,
        operator: LuminescenceTransferOperator<M>,
        bindings: &[RadiativeThermalCell],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        tolerance: LuminescenceAbsorptionTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        let absorption = medium.absorption_grid();
        if operator.medium_identity != medium.content_identity()
            || operator.absorption_identity != absorption.content_identity()
            || operator.wavelength_bounds_nm != medium.wavelength_bounds_nm()
            || operator.grid_cell_count != absorption.cell_count()
            || absorption.cell_count() > max_grid_cells
            || bindings.len() > max_grid_cells
            || operator.material_cells.len() > max_grid_cells
            || !nonnegative(tolerance.fraction_roundoff)
            || tolerance.fraction_roundoff >= 1.0
            || !nonnegative(tolerance.max_fraction_standard_error)
            || !nonnegative(tolerance.max_unresolved_energy_j)
        {
            return Err(invalid(
                "luminescence transport snapshot, budget or tolerance mismatch",
            ));
        }
        let mut refused = None;
        let validation = operator.validate_with_checkpoint(
            &mut || match checkpoint() {
                Ok(()) => Ok(()),
                Err(error) => {
                    refused = Some(error);
                    Err(VolumeTransportError::Cancelled)
                }
            },
            tolerance.fraction_roundoff,
        );
        if let Some(error) = refused {
            return Err(error.into());
        }
        validation.map_err(|_| invalid("invalid luminescence photon-transfer report"))?;
        // A public report cannot omit an unilluminated material cell whose
        // optical validity may change under another native heat source.
        let mut material = operator.material_cells.iter().copied();
        for cell in 0..absorption.cell_count() {
            checkpoint()?;
            if absorption
                .cell_density_kg_m3(cell)
                .is_some_and(|rho| rho > 0.0)
                && material.next() != Some(cell)
            {
                return Err(invalid(
                    "luminescence transport omitted an original material cell",
                ));
            }
        }
        if material.next().is_some() || bindings.len() != operator.material_cells.len() {
            return Err(invalid(
                "luminescence requires exact original material bindings",
            ));
        }
        let mut receivers = BTreeMap::new();
        let mut domains = BTreeMap::<CellId, (f64, f64)>::new();
        for binding in bindings {
            checkpoint()?;
            let (lo, hi) = (binding.minimum_temperature_k, binding.maximum_temperature_k);
            if !lo.is_finite()
                || lo <= 0.0
                || !hi.is_finite()
                || hi < lo
                || operator
                    .material_cells
                    .binary_search(&binding.grid_cell)
                    .is_err()
            {
                return Err(invalid(
                    "invalid luminescence receiver or optical temperature domain",
                ));
            }
            if receivers
                .insert(binding.grid_cell, binding.thermal_cell)
                .is_some()
            {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell));
            }
            if !thermal.cells().contains_key(&binding.thermal_cell) {
                return Err(OpticalHeatingError::InactiveCell(binding.thermal_cell));
            }
            let domain = domains.entry(binding.thermal_cell).or_insert((lo, hi));
            domain.0 = domain.0.max(lo);
            domain.1 = domain.1.min(hi);
            if domain.0 > domain.1 {
                return Err(invalid(
                    "one thermal receiver has disjoint optical temperature domains",
                ));
            }
        }
        let mut row_receivers = Vec::with_capacity(operator.rows.len());
        for row in &operator.rows {
            checkpoint()?;
            for estimate in [&row.absorbed, &row.escape, &row.unresolved]
                .into_iter()
                .chain(row.absorption.iter().map(|link| &link.fraction))
            {
                if estimate
                    .standard_error
                    .iter()
                    .any(|se| *se > tolerance.max_fraction_standard_error)
                {
                    return Err(invalid(
                        "luminescence transfer sampling error exceeds the supplied allowance",
                    ));
                }
            }
            let mut ids = Vec::with_capacity(row.absorption.len());
            for link in &row.absorption {
                checkpoint()?;
                ids.push(
                    *receivers
                        .get(&link.cell_index)
                        .ok_or(OpticalHeatingError::UnmappedCell(link.cell_index))?,
                );
            }
            row_receivers.push(ids);
        }
        let result = Self {
            model: Arc::new(Model {
                operator,
                bindings: receivers,
                row_receivers,
                temperature_domains: domains
                    .into_iter()
                    .map(|(id, (lo, hi))| (id, lo, hi))
                    .collect(),
                tolerance,
            }),
        };
        result.check_state(thermal, &mut checkpoint)?;
        Ok(result)
    }

    /// Original directed probabilities and sampling diagnostics, including
    /// self-return. They are never renormalized after admission.
    pub fn operator(&self) -> &LuminescenceTransferOperator<M> {
        &self.model.operator
    }

    pub fn tolerance(&self) -> LuminescenceAbsorptionTolerance {
        self.model.tolerance
    }

    pub(super) fn check_emitter(
        &self,
        grid_cell: usize,
        thermal_cell: CellId,
    ) -> Result<(), OpticalHeatingError> {
        if self.model.bindings.get(&grid_cell) != Some(&thermal_cell)
            || self
                .model
                .operator
                .rows
                .binary_search_by_key(&grid_cell, |r| r.emitter_cell)
                .is_err()
        {
            return Err(invalid(
                "luminescence emitter lacks its original transport row or thermal binding",
            ));
        }
        Ok(())
    }

    pub(super) fn check_state(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalHeatingError> {
        for &(id, lo, hi) in &self.model.temperature_domains {
            checkpoint()?;
            let temperature = thermal.temperature_k(id)?;
            if !temperature.is_finite() || !(lo..=hi).contains(&temperature) {
                return Err(invalid(
                    "native temperature exceeds the luminescence transport snapshot domain",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn route(
        &self,
        grid_cell: usize,
        emitted_j: [f64; M],
        duration_s: f64,
        sources: &mut BTreeMap<CellId, CheckedSum>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<RoutedEnergy<M>, OpticalHeatingError> {
        let index = self
            .model
            .operator
            .rows
            .binary_search_by_key(&grid_cell, |r| r.emitter_cell)
            .map_err(|_| invalid("missing luminescence transfer row"))?;
        let row = &self.model.operator.rows[index];
        let mut received = [CheckedSum::default(); M];
        for (link, &receiver) in row.absorption.iter().zip(&self.model.row_receivers[index]) {
            checkpoint()?;
            let mut heat = CheckedSum::default();
            for lane in 0..M {
                let energy = finite(emitted_j[lane] * link.fraction.mean[lane])?;
                received[lane].add(energy)?;
                heat.add(energy)?;
            }
            sources
                .entry(receiver)
                .or_default()
                .add(finite(heat.value()? / duration_s)?)?;
        }
        let mut result = RoutedEnergy {
            escaped_j: [0.0; M],
            unresolved_j: [0.0; M],
            heat_j: 0.0,
            absolute_roundoff_j: 0.0,
        };
        let mut heat = CheckedSum::default();
        let mut defects = CheckedSum::default();
        for lane in 0..M {
            result.escaped_j[lane] = finite(emitted_j[lane] * row.escape.mean[lane])?;
            result.unresolved_j[lane] = finite(emitted_j[lane] * row.unresolved.mean[lane])?;
            let received = received[lane].value()?;
            heat.add(received)?;
            let mut balance = CheckedSum::default();
            for energy in [
                emitted_j[lane],
                -received,
                -result.escaped_j[lane],
                -result.unresolved_j[lane],
            ] {
                balance.add(energy)?;
            }
            defects.add(balance.value()?.abs())?;
        }
        result.heat_j = heat.value()?;
        result.absolute_roundoff_j = defects.value()?;
        Ok(result)
    }
}
