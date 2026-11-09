//! Neutral finite-rate crystals on the original solvent-frame diffusion path.
//!
//! Ideal phase owners require an exactly zero-virial diffusion view. Nonzero
//! interactions require each original kinetic material to be explicitly bound
//! to that same nonideal model, so activities change its conditional saturation.
//! Crystal stocks stay local; original phase half-stages surround transport of
//! dissolved species. Their solid chemical reference closes the whole-window
//! free energy without treating dissolution as a source of component mass.

use super::{
    AqueousSolution, OsmoticCompartment, OsmoticNetwork, OsmoticNetworkError,
    SolutionDiffusionControls, SolutionDiffusionLink, compensated_add, finite, poll,
    prepare_solution_diffusion,
};
use fs_material::osmosis::nonideal::VirialSolutionModel;
use std::sync::Arc;

impl OsmoticNetwork {
    /// Couple original finite-rate neutral solids to solvent-frame molecular
    /// diffusion at fixed temperature and solvent volumes. Only dissolved
    /// molecules cross faces; fixed-area or particle-area crystals remain in
    /// their original material compartments. The native full/two-half timeline
    /// controls the complete first-order phase/transport splitting error.
    ///
    /// The explicitly supplied common thermodynamic model must have the same
    /// temperature, solvent molar volume and solute ordering as the original
    /// phase owners. An ideal phase owner requires exactly zero virial entries;
    /// nonzero entries require every kinetic owner to have explicitly admitted
    /// this same model through `with_nonideal_thermodynamics`. Its conditional
    /// activity equilibrium and concentration box remain binding. All kinetic
    /// cells retain the SAME original ordered
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
    /// Instantaneous equilibrium solids, ionic assignments, coupled water
    /// exchange, nucleation and crystal transport remain outside
    /// this constructor. Refusal/cancellation publishes no state or time.
    pub fn with_kinetic_solution_diffusion(
        species: Vec<String>,
        mut cells: Vec<OsmoticCompartment>,
        mut links: Vec<SolutionDiffusionLink>,
        thermodynamics: Arc<VirialSolutionModel>,
        controls: SolutionDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        prepare_kinetic_thermodynamics(&species, &mut cells, &links, &thermodynamics)?;
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

/// Shared phase/diffusion admission for fixed-volume and coupled-water
/// networks. Only the aqueous view changes representation; the original phase
/// owners, actual stocks and source/model identities remain intact.
pub(in crate::osmosis) fn prepare_kinetic_thermodynamics(
    species: &[String],
    cells: &mut [OsmoticCompartment],
    links: &[SolutionDiffusionLink],
    thermodynamics: &Arc<VirialSolutionModel>,
) -> Result<(), OsmoticNetworkError> {
    if thermodynamics.species_count() != species.len()
        || cells.iter().any(|cell| cell.solubility.is_some())
    {
        return Err(OsmoticNetworkError::Invalid(
            "kinetic molecular diffusion requires neutral finite-rate stocks and the declared species count",
        ));
    }
    let interacting = thermodynamics
        .interaction_matrix_m3_mol()
        .iter()
        .any(|&entry| entry != 0.0);
    if interacting
        && cells
            .iter()
            .filter_map(|cell| cell.kinetic_solubility.as_ref())
            .any(|phases| phases.solution().nonideal_model() != Some(thermodynamics.as_ref()))
    {
        return Err(OsmoticNetworkError::Invalid(
            "nonideal molecular diffusion requires every original kinetic owner to admit the same activity model",
        ));
    }
    let reference = cells
        .iter()
        .find_map(|cell| cell.kinetic_solubility.as_ref())
        .ok_or(OsmoticNetworkError::Invalid(
            "kinetic molecular diffusion requires an original kinetic material owner",
        ))?;
    if links
        .iter()
        .any(|link| link.electrical.is_some() || link.material().solvent() != reference.solvent())
    {
        return Err(OsmoticNetworkError::Invalid(
            "neutral diffusion faces must retain the original phase solvent identity",
        ));
    }
    for cell in cells {
        if cell.solution.temperature_k().to_bits() != thermodynamics.temperature_k().to_bits()
            || cell.solution.solvent_molar_volume_m3_mol().to_bits()
                != thermodynamics.solvent_molar_volume_m3_mol().to_bits()
            || (!cell.solution.is_ideal()
                && cell.solution.nonideal_model() != Some(thermodynamics.as_ref()))
        {
            return Err(OsmoticNetworkError::Invalid(
                "kinetic diffusion must retain the original temperature, solvent convention and activity model",
            ));
        }
        // Preparation allocates once; subsequent updates reuse this view.
        cell.solution = AqueousSolution::nonideal(
            Arc::clone(thermodynamics),
            cell.solution.solvent_volume_m3(),
            cell.solution.solute_moles().to_vec(),
        )?;
    }
    Ok(())
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
        // The original phase owner retains the same pure-solid reference and,
        // when explicitly selected, its calibrated nonideal excess potential.
        // Its linear reference cancels only for the conserved COMPLETE stocks.
        let change = phases.solid_relative_inventory_free_energy_change_j(
            before.solvent_volume_m3(),
            before.solute_moles(),
            after.solvent_volume_m3(),
            after.solute_moles(),
        )?;
        compensated_add(&mut sum, &mut correction, finite(change)?);
    }
    finite(sum)
}
