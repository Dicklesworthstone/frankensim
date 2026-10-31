//! Original-constituent absorption from one native beam transport operation.
//!
//! Transport uses the exact compiled mixture. This owner partitions only its
//! true absorbed power, at every original wavelength, using the constituent
//! mass fractions and absorption laws retained by that same accepted mixture.
//! Scattered radiation remains in the unchanged native beam account.

/// Constituent means and per-cell errors from the original scattering histories.
pub mod scattering;

use super::{
    MaterialOpticalComponent, MaterialOpticalMixtureError, MaterialOpticalMixtureGrid, checkpoint,
};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::deposition::BeamPowerDeposition;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Additional output and constitutive work limits, beyond native ray traversal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstituentDepositionBudget {
    /// Maximum retained positive-absorption `(cell, component)` rows.
    pub max_rows: usize,
    /// Maximum component/wavelength slots examined in absorbing cells.
    /// Zero-power lanes and nonabsorbing components still consume one slot.
    pub max_component_evaluations: usize,
}

/// Absorption by one original constituent in one physical optical cell.
/// Indices remain the original cell and recipe-component ordinals.
#[derive(Clone, Debug, PartialEq)]
pub struct CellConstituentAbsorption<const N: usize> {
    optical_cell_index: usize,
    component_index: usize,
    source_identity: ContentHash,
    absorbed_power_w: [f64; N],
}

impl<const N: usize> CellConstituentAbsorption<N> {
    pub const fn optical_cell_index(&self) -> usize {
        self.optical_cell_index
    }
    pub const fn component_index(&self) -> usize {
        self.component_index
    }
    /// Identity of the original material/particle binding, not an invented blend.
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    pub const fn absorbed_power_w(&self) -> &[f64; N] {
        &self.absorbed_power_w
    }
}

/// Trusted constituent partition retaining the complete original native report.
///
/// Private fields prevent substitution of another medium's report or mutated
/// constituent powers before thermal binding. Rows follow native ray traversal,
/// then original component order. Omitted components/cells have exactly zero
/// absorbed power at every requested wavelength; no positive underflowed share
/// is silently published as zero.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstituentBeamDeposition<const N: usize> {
    mixture_identity: ContentHash,
    beam: BeamPowerDeposition<N>,
    cells: Vec<CellConstituentAbsorption<N>>,
    component_evaluations: usize,
    partition_residual_w: [f64; N],
}

impl<const N: usize> ConstituentBeamDeposition<N> {
    pub const fn mixture_identity(&self) -> ContentHash {
        self.mixture_identity
    }
    /// Unchanged native incident/transmitted/absorbed/scattered energy balance.
    pub const fn beam(&self) -> &BeamPowerDeposition<N> {
        &self.beam
    }
    pub fn cells(&self) -> &[CellConstituentAbsorption<N>] {
        &self.cells
    }
    pub const fn component_evaluations(&self) -> usize {
        self.component_evaluations
    }
    /// Native absorbed watts minus the compensated sum of constituent watts.
    /// This is retained roundoff, not an extra heat source or redistributed heat.
    pub const fn partition_residual_w(&self) -> &[f64; N] {
        &self.partition_residual_w
    }
}

impl MaterialOpticalMixtureGrid<'_> {
    /// Trace one collimated beam ONCE, then resolve absorption by original source.
    ///
    /// For each wavelength and cell, constituent `i` receives
    /// `P_abs * w_i*kappa_ai / sum_j(w_j*kappa_aj)`. This conditions on actual
    /// absorption; it never normalizes the declared mixture mass fractions.
    /// Purely scattering gas/particles receive exactly zero heat. Zero incident
    /// power and cells with no absorption need no partition rows or evaluations.
    ///
    /// The native beam budget bounds ray output independently. Constituent
    /// admission and allocation are additionally bounded by `budget`; failure
    /// or cancellation returns no partial report. Optics remain this accepted
    /// snapshot, and emitted radiation is not part of incident-beam deposition.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_beam_power_by_component<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        max_deposition_cells: usize,
        budget: ConstituentDepositionBudget,
    ) -> Result<ConstituentBeamDeposition<N>, MaterialOpticalMixtureError> {
        checkpoint(cx)?;
        let beam = self.medium().deposit_beam_power(
            cx,
            origin,
            unit_direction,
            max_distance_m,
            wavelengths_nm,
            incident_power_w,
            max_deposition_cells,
        )?;
        if beam.medium_identity != self.medium().content_identity() {
            return Err(VolumeTransportError::SourceIdentity.into());
        }
        let mut cells = Vec::new();
        let mut component_evaluations = 0_usize;
        let mut totals = [PositiveSum::default(); N];
        for cell in &beam.cells {
            checkpoint(cx)?;
            if cell.absorbed_power_w.iter().all(|power| *power == 0.0) {
                continue;
            }
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
                    // This gas optical owner supplies elastic Rayleigh scattering
                    // only; no absorption is inferred from extinction.
                    continue;
                };
                let mut absorption = [0.0; N];
                for lane in 0..N {
                    if cell.absorbed_power_w[lane] == 0.0 {
                        continue;
                    }
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
                // Private construction temporarily retains the original kappa;
                // no row is published until all shares and balances validate.
                cells.push(CellConstituentAbsorption {
                    optical_cell_index: cell.cell_index,
                    component_index,
                    source_identity: particles.binding_identity(),
                    absorbed_power_w: absorption,
                });
            }
            let mut partition = [PositiveSum::default(); N];
            for row in &mut cells[first_row..] {
                checkpoint(cx)?;
                let weight = mixture.composition().fractions()[row.component_index];
                for lane in 0..N {
                    let absorbed = absorbed_share(
                        cell.absorbed_power_w[lane],
                        weight,
                        row.absorbed_power_w[lane],
                        coefficients[lane].value()?,
                    )?;
                    row.absorbed_power_w[lane] = absorbed;
                    partition[lane].add(absorbed)?;
                    totals[lane].add(absorbed)?;
                }
            }
            for lane in 0..N {
                admit_closure(
                    cell.absorbed_power_w[lane],
                    partition[lane].value()?,
                    cells.len() - first_row,
                )?;
            }
        }
        let mut partition_residual_w = [0.0; N];
        for lane in 0..N {
            let total = totals[lane].value()?;
            admit_closure(beam.absorbed_power_w[lane], total, cells.len())?;
            partition_residual_w[lane] = beam.absorbed_power_w[lane] - total;
        }
        checkpoint(cx)?;
        Ok(ConstituentBeamDeposition {
            mixture_identity: self.content_identity(),
            beam,
            cells,
            component_evaluations,
            partition_residual_w,
        })
    }
}

fn absorbed_share(
    power: f64,
    weight: f64,
    coefficient: f64,
    total: f64,
) -> Result<f64, VolumeTransportError> {
    if power == 0.0 || weight == 0.0 || coefficient == 0.0 {
        return Ok(0.0);
    }
    if !total.is_finite() || total <= 0.0 {
        return Err(VolumeTransportError::Unrepresentable);
    }
    let fraction = weight * coefficient / total;
    let direct = power * fraction;
    let value = if fraction >= f64::MIN_POSITIVE && direct > 0.0 && direct.is_finite() {
        direct
    } else {
        // A fraction below f64 range can still carry representable watts.
        // Keep original factors instead of taking ln of an underflowed product.
        det::exp(det::ln(power) + det::ln(weight) + det::ln(coefficient) - det::ln(total))
    };
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(VolumeTransportError::Unrepresentable)
    }
}

fn admit_closure(expected: f64, actual: f64, terms: usize) -> Result<(), VolumeTransportError> {
    let allowance =
        expected * (256.0 * f64::EPSILON * terms.max(1) as f64) + 8.0 * f64::from_bits(1);
    if !actual.is_finite() || (expected - actual).abs() > allowance {
        Err(VolumeTransportError::Unrepresentable)
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct PositiveSum {
    sum: f64,
    correction: f64,
}
impl PositiveSum {
    fn add(&mut self, value: f64) -> Result<(), VolumeTransportError> {
        let next = self.sum + value;
        if !value.is_finite() || value < 0.0 || !next.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        Ok(())
    }
    fn value(self) -> Result<f64, VolumeTransportError> {
        let value = self.sum + self.correction;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(VolumeTransportError::Unrepresentable)
        }
    }
}

/// Explicit constituent-to-cell sources on the existing thermal timeline.
#[cfg(feature = "optical-heating")]
pub mod heating;

#[cfg(test)]
mod arithmetic_tests {
    use super::*;

    #[test]
    fn representable_absorbed_watts_survive_fraction_underflow() {
        let actual = absorbed_share(1e300, 0.5, 2e-300, 1e100).unwrap();
        assert!((actual / 1e-100 - 1.0).abs() < 3e-13);
        assert_eq!(absorbed_share(1e300, 0.5, 0.0, 1e100).unwrap(), 0.0);
        assert_eq!(absorbed_share(0.0, 0.5, 2e-300, 0.0).unwrap(), 0.0);
    }
}
