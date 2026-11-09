//! Neutral finite-rate crystals on the original solvent-frame diffusion path.
//!
//! The kinetic material's ideal particle pressure and the exactly zero-virial
//! diffusion view are the SAME thermodynamics. Nonzero interactions would
//! change saturation through activities and are therefore refused, not ignored.
//! Crystal stocks stay local; original phase half-stages surround transport of
//! dissolved species. Their solid chemical reference closes the whole-window
//! free energy without treating dissolution as a source of component mass.

use super::{
    AqueousSolution, OsmoticCompartment, OsmoticNetwork, OsmoticNetworkError,
    SolutionDiffusionControls, SolutionDiffusionLink, compensated_add, finite, poll,
    prepare_solution_diffusion,
};
use fs_material::osmosis::nonideal::VirialSolutionModel;
use fs_qty::{Amount, Temperature, Volume};
use std::sync::Arc;

impl OsmoticNetwork {
    /// Couple original finite-rate neutral solids to solvent-frame molecular
    /// diffusion at fixed temperature and solvent volumes. Only dissolved
    /// molecules cross faces; fixed-area or particle-area crystals remain in
    /// their original material compartments. The native full/two-half timeline
    /// controls the complete first-order phase/transport splitting error.
    ///
    /// The explicitly supplied common thermodynamic model MUST have exactly
    /// zero virial interactions and the same temperature, solvent molar volume
    /// and solute ordering as the original phase owners. Its concentration box
    /// remains binding. All kinetic cells retain the SAME original ordered
    /// solubility laws and pure-solid reference, while rate coefficients,
    /// particle geometry and growth policies can differ. Ordinary ideal cells
    /// can share this network. At least one kinetic compartment is required.
    ///
    /// Existing neutral explicit, implicit binary and implicit multisolute
    /// stages consume the same original Maxwell--Stefan friction data. This is
    /// not an inferred Fick coefficient or a change of reference frame.
    /// Source-backed faces retain their original qualification and transient
    /// query-allocation costs. Constant-data explicit and implicit binary
    /// advancement allocates no heap; the existing implicit multisolute solver
    /// retains its declared transient-vector allowance.
    ///
    /// The closed free-energy report includes the solid-relative Gibbs datum
    /// from each original solubility law; it is not a thermal energy account.
    /// Instantaneous equilibrium solids, ionic assignments, nonideal saturation,
    /// coupled water exchange, nucleation and crystal transport remain outside
    /// this constructor. Refusal/cancellation publishes no state or time.
    pub fn with_kinetic_solution_diffusion(
        species: Vec<String>,
        mut cells: Vec<OsmoticCompartment>,
        mut links: Vec<SolutionDiffusionLink>,
        thermodynamics: Arc<VirialSolutionModel>,
        controls: SolutionDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if thermodynamics.species_count() != species.len()
            || thermodynamics
                .interaction_matrix_m3_mol()
                .iter()
                .any(|&entry| entry != 0.0)
            || cells.iter().any(|cell| cell.solubility.is_some())
        {
            return Err(OsmoticNetworkError::Invalid(
                "kinetic molecular diffusion requires neutral finite-rate stocks and exactly zero virial interactions",
            ));
        }
        let reference = cells
            .iter()
            .find_map(|cell| cell.kinetic_solubility.as_ref())
            .ok_or(OsmoticNetworkError::Invalid(
                "kinetic molecular diffusion requires an original kinetic material owner",
            ))?;
        if links.iter().any(|link| {
            link.electrical.is_some() || link.material().solvent() != reference.solvent()
        }) {
            return Err(OsmoticNetworkError::Invalid(
                "neutral diffusion faces must retain the original phase solvent identity",
            ));
        }
        for cell in &mut cells {
            if cell.solution.temperature_k().to_bits() != thermodynamics.temperature_k().to_bits()
                || cell.solution.solvent_molar_volume_m3_mol().to_bits()
                    != thermodynamics.solvent_molar_volume_m3_mol().to_bits()
                || (!cell.solution.is_ideal()
                    && cell.solution.nonideal_model() != Some(thermodynamics.as_ref()))
            {
                return Err(OsmoticNetworkError::Invalid(
                    "kinetic diffusion must retain the original ideal temperature and solvent convention",
                ));
            }
            // This view changes neither physical amounts nor the phase owner.
            // Preparation allocates once; all subsequent updates reuse it.
            cell.solution = AqueousSolution::nonideal(
                Arc::clone(&thermodynamics),
                cell.solution.solvent_volume_m3(),
                cell.solution.solute_moles().to_vec(),
            )?;
        }
        let required = prepare_solution_diffusion(
            &species,
            &cells,
            &mut links,
            controls,
            maximum_workspace_bytes,
            true,
        )?;
        let mut network = Self::new(
            species,
            cells,
            Vec::new(),
            maximum_workspace_bytes - required,
        )?;
        network.workspace_bytes += required;
        network.diffusion = links;
        network.diffusion_controls = Some(controls);
        Ok(network)
    }
}

pub(super) fn free_energy_change(
    original: &[OsmoticCompartment],
    candidate: &[OsmoticCompartment],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, OsmoticNetworkError> {
    let phases = original
        .iter()
        .find_map(|cell| cell.kinetic_solubility.as_ref())
        .ok_or(OsmoticNetworkError::Invalid(
            "missing original solid chemical reference",
        ))?;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for (before, after) in original.iter().zip(candidate) {
        poll(cancelled)?;
        let (before, after) = (&before.solution, &after.solution);
        if before.solvent_volume_m3().to_bits() != after.solvent_volume_m3().to_bits()
            || before.nonideal_model() != after.nonideal_model()
        {
            return Err(OsmoticNetworkError::Invalid(
                "closed kinetic diffusion changed its solvent volume or original thermodynamics",
            ));
        }
        for ((model, &old), &next) in phases
            .models()
            .iter()
            .zip(before.solute_moles())
            .zip(after.solute_moles())
        {
            // Subtract the SAME pure-solid chemical potential from every
            // dissolved and solid mole. Its linear reference cancels because
            // complete component inventories are conserved across this window.
            // The remaining solid contribution is exactly zero by that datum,
            // not an omitted phase energy or an assumption of fixed crystals.
            let change = model
                .law()
                .solid_relative_free_energy_change_j(
                    Temperature::new(before.temperature_k()),
                    Volume::new(before.solvent_volume_m3()),
                    Amount::new(old),
                    Amount::new(next),
                )
                .map_err(|error| OsmoticNetworkError::Solubility(error.into()))?;
            compensated_add(&mut sum, &mut correction, finite(change)?);
        }
    }
    finite(sum)
}
