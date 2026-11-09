//! Frozen-coefficient refinement on the original reservoir caloric owner.
//!
//! The underlying transaction retains the true cell and face transfer low parts.
//! Comparing those accounts prevents a common formation datum or opposing
//! boundaries from hiding temporal error in rounded public readbacks.

use super::super::super::{
    AdaptiveSpeciesDiffusionControls, SpeciesBoundaryDiffusionReport, SpeciesBoundaryTransfer,
    SpeciesEnthalpyDiffusionReport, positive,
};
use super::{
    BoundaryEnthalpyAccount, BoundaryEnthalpySums, BoundarySpeciesEnthalpyDiffusionStep,
    CancelGate, ChtError, EnergySum, MaterialAdvectionStep, SpeciesBoundaryEnthalpyDiffusionReport,
    SpeciesBoundaryEnthalpyTransfer, SpeciesDiffusionBoundary, SpeciesEnthalpyDiffusionStep,
    VoxelDomain, nonnegative, poll, poll_every, refuse,
};

/// Independent absolute scales for the original frozen boundary caloric model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundarySpeciesEnthalpyAccuracy<const N: usize> {
    pub absolute_boundary_transfer_tolerance_kg: [f64; N],
    /// Maximum per-cell error [J/kg] from retained increments / carrier mass.
    /// A formation datum neither scales nor hides this error.
    pub absolute_specific_enthalpy_tolerance_j_kg: f64,
    pub absolute_boundary_enthalpy_tolerance_j: f64,
}

/// Accepted original fine candidate plus complete-interval accuracy/work.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveBoundarySpeciesEnthalpyDiffusionStep<const N: usize> {
    pub step: BoundarySpeciesEnthalpyDiffusionStep<N>,
    pub estimated_mass_fraction_error: [f64; N],
    pub estimated_specific_enthalpy_error_j_kg: f64,
    pub estimated_boundary_transfer_error_kg: [f64; N],
    pub estimated_boundary_enthalpy_error_j: f64,
    pub maximum_error_ratio: f64,
    pub refinements: usize,
    /// Includes every attempted coarse/fine solve, including rejected ones.
    pub total_pair_updates: usize,
    pub total_component_updates: usize,
}

struct Trial<const N: usize> {
    h: Vec<f64>,
    heat: Vec<f64>,
    report: SpeciesEnthalpyDiffusionReport,
    boundary: SpeciesBoundaryDiffusionReport<N>,
    thermal_boundary: SpeciesBoundaryEnthalpyDiffusionReport,
    account: BoundaryEnthalpyAccount<N>,
}

fn difference(mut a: EnergySum, b: EnergySum) -> f64 {
    a.add(-b.sum);
    a.add(-b.correction);
    a.value()
}

impl MaterialAdvectionStep {
    /// Refine the SAME complete frozen boundary caloric interval from original stores.
    ///
    /// Source callbacks have the original fixed owner's arguments and must
    /// return the same face laws at every resolution. They own the constituent
    /// ordering and common caloric convention; neither an EOS nor a physical
    /// surface temperature is inferred. D, reservoir compositions and films
    /// remain fixed throughout this call.
    ///
    /// The original controller checks each constituent concentration. This
    /// extension also checks compensated cell-energy differences / immutable
    /// carrier mass and per-face signed/gross constituent and energy integrals.
    /// For each face take the largest signed/incoming/outgoing difference;
    /// sum those magnitudes across faces, then divide by three. Opposing
    /// reservoirs cannot cancel unresolved exchange. Retained low parts, not
    /// rounded absolute h or scalar heat receipts, supply every caloric estimate.
    ///
    /// The accepted fine inventories/receipts are unmodified. Original work
    /// limits cover all attempted resolutions; storage includes retained coarse
    /// fields and original per-cell/per-face accounts. These asymptotic
    /// second-order estimates describe the frozen model only, without a spatial,
    /// coupled-source, EOS-refresh or phase-transition accuracy guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_with_enthalpy_and_boundaries_adaptive<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        diffusivity_m2_s: &[f64],
        boundaries: &[SpeciesDiffusionBoundary<N>],
        duration_s: f64,
        face_partial_enthalpy: impl Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
        boundary_partial_enthalpy: impl Fn(usize) -> Result<[f64; N], ChtError>,
        controls: AdaptiveSpeciesDiffusionControls<N>,
        accuracy: BoundarySpeciesEnthalpyAccuracy<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<AdaptiveBoundarySpeciesEnthalpyDiffusionStep<N>, ChtError> {
        poll(gate)?;
        positive(accuracy.absolute_specific_enthalpy_tolerance_j_kg)?;
        positive(accuracy.absolute_boundary_enthalpy_tolerance_j)?;
        for tolerance in accuracy.absolute_boundary_transfer_tolerance_kg {
            positive(tolerance)?;
        }
        let extra = domain
            .cell_count()
            .checked_mul(2 * core::mem::size_of::<f64>() + core::mem::size_of::<EnergySum>())
            .and_then(|bytes| {
                boundaries
                    .len()
                    .checked_mul(
                        core::mem::size_of::<SpeciesBoundaryTransfer<N>>()
                            + core::mem::size_of::<SpeciesBoundaryEnthalpyTransfer>()
                            + core::mem::size_of::<BoundaryEnthalpySums<N>>(),
                    )
                    .and_then(|faces| bytes.checked_add(faces))
            })
            .and_then(|bytes| {
                bytes.checked_add(
                    core::mem::size_of::<Trial<N>>()
                        + core::mem::size_of::<AdaptiveBoundarySpeciesEnthalpyDiffusionStep<N>>()
                        + 4 * core::mem::size_of::<[f64; N]>()
                        + 1024,
                )
            })
            .ok_or_else(|| refuse("adaptive boundary enthalpy storage extent overflows"))?;
        let mut estimated_h = 0.0_f64;
        let mut estimated_mass = [0.0; N];
        let mut estimated_energy = 0.0;
        let (species, fine) = self.refine_species(
            domain,
            species_mass_kg,
            controls,
            gate,
            |_| poll(gate),
            extra,
            |solver| {
                let (result, account) = self
                    .diffuse_species_with_enthalpy_and_boundaries_accounted(
                        domain,
                        species_mass_kg,
                        specific_enthalpy_j_kg,
                        diffusivity_m2_s,
                        boundaries,
                        duration_s,
                        &face_partial_enthalpy,
                        &boundary_partial_enthalpy,
                        solver,
                        maximum_enthalpy_defect_j,
                        gate,
                    )?;
                let SpeciesEnthalpyDiffusionStep {
                    species,
                    specific_enthalpy_j_kg,
                    diffusive_enthalpy_change_j,
                    report,
                } = result.step;
                Ok((
                    species,
                    Trial {
                        h: specific_enthalpy_j_kg,
                        heat: diffusive_enthalpy_change_j,
                        report,
                        boundary: result.boundary,
                        thermal_boundary: result.thermal_boundary,
                        account,
                    },
                ))
            },
            |coarse, fine| {
                estimated_h = 0.0;
                for (cell, (&a, &b)) in coarse
                    .account
                    .increments
                    .iter()
                    .zip(&fine.account.increments)
                    .enumerate()
                {
                    poll_every(gate, cell)?;
                    if !domain.is_fluid(cell) {
                        continue;
                    }
                    let error = difference(b, a).abs() / self.mass_kg[cell] / 3.0;
                    nonnegative(error)?;
                    estimated_h = estimated_h.max(error);
                }
                let mut maximum_ratio =
                    estimated_h / accuracy.absolute_specific_enthalpy_tolerance_j_kg;
                let mut mass = [EnergySum::default(); N];
                let mut energy = EnergySum::default();
                for (index, (a, b)) in coarse
                    .account
                    .faces
                    .iter()
                    .zip(&fine.account.faces)
                    .enumerate()
                {
                    poll_every(gate, index)?;
                    for k in 0..N {
                        mass[k].add(
                            difference(b.mass_signed[k], a.mass_signed[k])
                                .abs()
                                .max(difference(b.mass_incoming[k], a.mass_incoming[k]).abs())
                                .max(difference(b.mass_outgoing[k], a.mass_outgoing[k]).abs())
                                / 3.0,
                        );
                    }
                    energy.add(
                        difference(b.signed, a.signed)
                            .abs()
                            .max(difference(b.positive, a.positive).abs())
                            .max(difference(b.negative, a.negative).abs())
                            / 3.0,
                    );
                }
                for k in 0..N {
                    estimated_mass[k] = mass[k].value();
                    nonnegative(estimated_mass[k])?;
                    maximum_ratio = maximum_ratio.max(
                        estimated_mass[k] / accuracy.absolute_boundary_transfer_tolerance_kg[k],
                    );
                }
                estimated_energy = energy.value();
                nonnegative(estimated_energy)?;
                Ok(maximum_ratio
                    .max(estimated_energy / accuracy.absolute_boundary_enthalpy_tolerance_j))
            },
        )?;
        Ok(AdaptiveBoundarySpeciesEnthalpyDiffusionStep {
            step: BoundarySpeciesEnthalpyDiffusionStep {
                step: SpeciesEnthalpyDiffusionStep {
                    species: species.step,
                    specific_enthalpy_j_kg: fine.h,
                    diffusive_enthalpy_change_j: fine.heat,
                    report: fine.report,
                },
                boundary: fine.boundary,
                thermal_boundary: fine.thermal_boundary,
            },
            estimated_mass_fraction_error: species.estimated_mass_fraction_error,
            estimated_specific_enthalpy_error_j_kg: estimated_h,
            estimated_boundary_transfer_error_kg: estimated_mass,
            estimated_boundary_enthalpy_error_j: estimated_energy,
            maximum_error_ratio: species.maximum_error_ratio,
            refinements: species.refinements,
            total_pair_updates: species.total_pair_updates,
            total_component_updates: species.total_component_updates,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_compensated_error_resolves_changes_below_a_large_readback() {
        let mut coarse = EnergySum::default();
        coarse.add(1e16);
        let mut fine = coarse;
        fine.add(1.0);
        assert_eq!(coarse.value(), fine.value());
        assert_eq!(difference(fine, coarse), 1.0);
        assert_eq!(difference(coarse, fine), -1.0);
        // Opposing large transfers leave an original small integral.
        fine.add(-1e16);
        coarse.add(-1e16);
        assert_eq!(difference(fine, coarse), 1.0);
    }
}
