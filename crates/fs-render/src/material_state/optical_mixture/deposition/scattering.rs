//! Constituent absorption of the original native multiple-scattering estimates.

#[cfg(feature = "optical-heating")]
/// Native thermal observer for constituent-resolved scattered absorption.
pub mod heating;

use super::{
    ConstituentDepositionBudget, MaterialOpticalComponent, MaterialOpticalMixtureError,
    MaterialOpticalMixtureGrid, PositiveSum, VolumeTransportError, absorbed_share, admit_closure,
    checkpoint,
};
use crate::volumes::transport::scattering::deposition::{
    ScatteringDepositionConfig, ScatteringPowerDeposition,
};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Fixed-sample absorbed watts for an original constituent in one native cell.
///
/// The constituent share is constant for every history at this cell and
/// wavelength. Both the cell mean and its estimated standard error therefore
/// scale by that same share; no independent constituent histories are invented.
#[derive(Clone, Debug, PartialEq)]
pub struct CellConstituentScatteringAbsorption<const N: usize> {
    optical_cell_index: usize,
    component_index: usize,
    source_identity: ContentHash,
    absorbed_power_w: [f64; N],
    standard_error_w: [Option<f64>; N],
}

impl<const N: usize> CellConstituentScatteringAbsorption<N> {
    pub const fn optical_cell_index(&self) -> usize {
        self.optical_cell_index
    }
    pub const fn component_index(&self) -> usize {
        self.component_index
    }
    /// Original material/particle binding, never an effective blended material.
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    /// Mean absorbed watts per spectral lane, including redirected photons.
    pub const fn absorbed_power_w(&self) -> &[f64; N] {
        &self.absorbed_power_w
    }
    /// Fixed-sample standard error of each constituent mean, in W.
    ///
    /// `None` means a sampled absorbing lane has no observed native deposition
    /// in this cell. `Some(0)` denotes deterministic or structurally zero heat,
    /// or an observed estimate with zero empirical variance. These estimates
    /// are not confidence bounds or guarantees against unobserved rare events.
    pub const fn standard_error_w(&self) -> &[Option<f64>; N] {
        &self.standard_error_w
    }
}

/// One unchanged native multiple-scattering report and its constituent shares.
///
/// Rows preserve the native sparse ascending cell order, followed by original
/// component order. An absent sampled cell is UNOBSERVED, not demonstrated to
/// receive zero heat. Each retained component absorbs at some packet wavelength;
/// molecular Rayleigh-only constituents receive no rows and no heat.
///
/// No aggregate constituent standard error is reported across cells: the native
/// report retains correlated whole-history totals, but does not retain their
/// inter-cell covariance. Its escape, unresolved-power, sample/work counts,
/// source identity and statistical semantics remain unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstituentScatteringDeposition<const N: usize> {
    mixture_identity: ContentHash,
    scattering: ScatteringPowerDeposition<N>,
    cells: Vec<CellConstituentScatteringAbsorption<N>>,
    component_evaluations: usize,
    partition_residual_w: [f64; N],
}

impl<const N: usize> ConstituentScatteringDeposition<N> {
    pub const fn mixture_identity(&self) -> ContentHash {
        self.mixture_identity
    }
    /// Complete original native report, including all statistical diagnostics.
    pub const fn scattering(&self) -> &ScatteringPowerDeposition<N> {
        &self.scattering
    }
    pub fn cells(&self) -> &[CellConstituentScatteringAbsorption<N>] {
        &self.cells
    }
    pub const fn component_evaluations(&self) -> usize {
        self.component_evaluations
    }
    /// Native mean absorbed watts minus the compensated constituent sum.
    /// Retained roundoff is not reassigned as an additional thermal source.
    pub const fn partition_residual_w(&self) -> &[f64; N] {
        &self.partition_residual_w
    }
}

impl MaterialOpticalMixtureGrid<'_> {
    /// Run the native multiple-scattering estimator ONCE and partition its heat.
    ///
    /// At each retained cell and original wavelength, constituent `i` receives
    /// `P_abs * w_i*kappa_ai / sum_j(w_j*kappa_aj)`. Its cell standard error
    /// receives the same fixed multiplier. Existing incident-source/history
    /// budgets and counter streams are passed through unchanged. No extra
    /// tracing, RNG draws, averaged optical coefficients, or thermal debit occurs.
    ///
    /// Additional output/evaluation work obeys `budget`; cancellation, budget
    /// exhaustion, or an unrepresentable positive share refuses the complete
    /// report. The sparse native output is not expanded to the full grid.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_multiple_scattering_power_by_component<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
        budget: ConstituentDepositionBudget,
    ) -> Result<ConstituentScatteringDeposition<N>, MaterialOpticalMixtureError> {
        checkpoint(cx)?;
        let scattering = self.medium().deposit_multiple_scattering_power(
            cx,
            origin,
            unit_direction,
            wavelengths_nm,
            incident_power_w,
            config,
        )?;
        if scattering.medium_identity != self.medium().content_identity() {
            return Err(VolumeTransportError::SourceIdentity.into());
        }
        let mut cells = Vec::new();
        let mut component_evaluations = 0_usize;
        let mut totals = [PositiveSum::default(); N];
        for cell in &scattering.cells {
            checkpoint(cx)?;
            let mixture = self
                .mixtures()
                .at_cell(cell.cell_index)
                .ok_or(VolumeTransportError::Palette)?;
            let work = mixture
                .components()
                .len()
                .checked_mul(N)
                .and_then(|work| component_evaluations.checked_add(work))
                .ok_or(VolumeTransportError::Budget)?;
            if work > budget.max_component_evaluations {
                return Err(VolumeTransportError::Budget.into());
            }
            component_evaluations = work;
            let first_row = cells.len();
            let mut coefficients = [PositiveSum::default(); N];
            for (component_index, (component, &weight)) in mixture
                .components()
                .iter()
                .zip(mixture.composition().fractions())
                .enumerate()
            {
                checkpoint(cx)?;
                let MaterialOpticalComponent::Particles(particles) = component else {
                    continue;
                };
                let mut absorption = [0.0; N];
                for lane in 0..N {
                    absorption[lane] = particles.spectrum().absorption.at(wavelengths_nm[lane])?;
                    coefficients[lane].add(weight * absorption[lane])?;
                }
                if absorption.iter().all(|coefficient| *coefficient == 0.0) {
                    continue;
                }
                if cells.len() == budget.max_rows {
                    return Err(VolumeTransportError::Budget.into());
                }
                if cells.len() == cells.capacity() {
                    cells
                        .try_reserve_exact(
                            cells.capacity().max(16).min(budget.max_rows - cells.len()),
                        )
                        .map_err(|_| VolumeTransportError::Budget)?;
                }
                // Hold original coefficients privately until the denominator
                // and every output share have been evaluated and admitted.
                cells.push(CellConstituentScatteringAbsorption {
                    optical_cell_index: cell.cell_index,
                    component_index,
                    source_identity: particles.binding_identity(),
                    absorbed_power_w: absorption,
                    standard_error_w: [Some(0.0); N],
                });
            }
            let mut partition = [PositiveSum::default(); N];
            for row in &mut cells[first_row..] {
                checkpoint(cx)?;
                let weight = mixture.composition().fractions()[row.component_index];
                for lane in 0..N {
                    let coefficient = row.absorbed_power_w[lane];
                    let total = coefficients[lane].value()?;
                    let mean = cell.absorbed.mean_w[lane];
                    let absorbed = absorbed_share(mean, weight, coefficient, total)?;
                    row.absorbed_power_w[lane] = absorbed;
                    row.standard_error_w[lane] =
                        if coefficient == 0.0 || scattering.sampled_histories[lane] == 0 {
                            Some(0.0)
                        } else if mean == 0.0 {
                            None
                        } else {
                            Some(absorbed_share(
                                cell.absorbed.standard_error_w[lane],
                                weight,
                                coefficient,
                                total,
                            )?)
                        };
                    partition[lane].add(absorbed)?;
                    totals[lane].add(absorbed)?;
                }
            }
            for lane in 0..N {
                admit_closure(
                    cell.absorbed.mean_w[lane],
                    partition[lane].value()?,
                    cells.len() - first_row,
                )?;
            }
        }
        let mut partition_residual_w = [0.0; N];
        for lane in 0..N {
            let total = totals[lane].value()?;
            admit_closure(scattering.absorbed.mean_w[lane], total, cells.len())?;
            partition_residual_w[lane] = scattering.absorbed.mean_w[lane] - total;
        }
        checkpoint(cx)?;
        Ok(ConstituentScatteringDeposition {
            mixture_identity: self.content_identity(),
            scattering,
            cells,
            component_evaluations,
            partition_residual_w,
        })
    }
}
