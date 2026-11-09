//! Retained finite-rate crystals in the existing water and molecular pathways.

use super::*;
use fs_material::osmosis::{AqueousSolution, nonideal::VirialSolutionModel};
use std::sync::Arc;

impl OsmoticNetwork {
    /// Evolve finite-rate neutral crystals, retaining water membranes and native
    /// stationary-solvent Maxwell--Stefan diffusion in one closed timeline.
    /// Only dissolved molecules cross the independent molecular links; water
    /// membranes retain every solute. Actual dissolved amounts drive osmosis.
    ///
    /// The supplied thermodynamic model must retain the original temperature,
    /// solvent molar volume and species ordering. Ideal phase owners require
    /// exactly zero virial entries. Nonzero interactions require every kinetic
    /// owner to have explicitly bound this SAME model through
    /// `with_nonideal_thermodynamics`; no ideal source is silently projected.
    /// Its concentration domain remains binding after water motion and phase
    /// exchange. Every crystal-bearing cell retains the same ordered
    /// pair laws and pure-solid references; film coefficients, particle stock
    /// references and growth policies may differ. At least one original kinetic
    /// owner is required. Ordinary aqueous receiving cells are also admitted.
    ///
    /// Native phase half-stages surround the existing water-then-molecular
    /// sweep. The original full/two-half controller therefore measures the
    /// complete first-order splitting error. Particle area evolves from actual
    /// local solid stock, while its reference geometry is never reset by water
    /// or solute transport. No crystal advection or new nucleation law is added.
    ///
    /// The closed potential includes variable-volume solid-relative chemical
    /// work and the exact pressure primitive, including prescribed pressure
    /// load work. Solvent and complete dissolved-plus-solid component balances
    /// use the original independent absolute budgets. This fixed-temperature
    /// account is not a heat balance. Cancellation or late failure publishes no
    /// partial phase stocks, solution volumes or elapsed time.
    ///
    /// Preparation and runtime allocation costs are those of the existing
    /// kinetic diffusion and coupled-water owners. Either path set may be
    /// empty. Ionic phases, caloric coupling and solvent-drag solute transport
    /// require their own constitutive owners.
    #[allow(clippy::too_many_arguments)]
    pub fn with_kinetic_coupled_solution_transport(
        species: Vec<String>,
        mut cells: Vec<OsmoticCompartment>,
        membranes: Vec<OsmoticMembrane>,
        links: Vec<SolutionDiffusionLink>,
        thermodynamics: Arc<VirialSolutionModel>,
        controls: CoupledSolutionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        diffusion::solubility::prepare_kinetic_thermodynamics(
            &species,
            &mut cells,
            &links,
            &thermodynamics,
        )?;
        Self::prepare_coupled_solution_transport(
            species,
            cells,
            membranes,
            links,
            controls,
            maximum_workspace_bytes,
            true,
        )
    }
}

pub(super) fn chemical_change(
    reference: &KineticSolubilitySolution,
    before: &AqueousSolution,
    after_volume: f64,
    after_amounts: &[f64],
) -> Result<f64, OsmoticNetworkError> {
    finite(reference.solid_relative_inventory_free_energy_change_j(
        before.solvent_volume_m3(),
        before.solute_moles(),
        after_volume,
        after_amounts,
    )?)
}
