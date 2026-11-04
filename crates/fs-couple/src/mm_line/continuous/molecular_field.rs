//! Original spatial molecular reservoirs on the existing exact shared clock.
use super::*;
use fs_duct::modal::time_domain::ModalMolecularBudget;
use fs_material::relaxation::RelaxationAbsorption;

impl ContinuousModalLine {
    /// Prepare the full source-backed modal field and every original local
    /// molecular reservoir on one exact flow. The actual gas at each retained
    /// Gauss point determines its relaxation strengths and rates. Local
    /// absorption and dispersion follow the same original material law.
    ///
    /// The prescribed background remains frozen. This adiabatic molecular
    /// model supplies no wall thermal diffusion or evolving gas temperature.
    /// All acoustic and memory states count against both preparation budgets.
    ///
    /// # Errors
    /// Original field, material, complete-capacity, cancellation and exact-flow refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new_field_molecular(
        profile: &RadiusProfile,
        field: &AxialGasField,
        law: RelaxationAbsorption,
        modes: ModeCount,
        elements: usize,
        end: ModalEnd,
        duration: f64,
        molecular_budget: ModalMolecularBudget,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless_field(
            profile,
            field,
            modes,
            elements,
            end,
            budget.max_modes,
            &mut cancelled,
        )?
        .with_field_molecular_relaxation(law, molecular_budget, &mut cancelled)?;
        let (source, system) = image.into_parts();
        Self::prepare(source, system, None, duration, budget, &mut cancelled)
    }

    /// Couple the complete spatial molecular field to every original mouth
    /// load state before preparing the shared exact flow. The mouth uses its
    /// actual source gas and the unchanged physical pressure/flow connection.
    /// Molecular heat and load dissipation both enter the original energy
    /// account; no attenuation is converted into a second heat source.
    ///
    /// # Errors
    /// The refusals of [`Self::new_field_molecular`], plus the original mouth
    /// binding and complete duct/reservoir/load capacity checks.
    #[allow(clippy::too_many_arguments)]
    pub fn new_field_molecular_loaded(
        profile: &RadiusProfile,
        field: &AxialGasField,
        law: RelaxationAbsorption,
        modes: ModeCount,
        elements: usize,
        load: ModalAdmittanceLoad,
        duration: f64,
        molecular_budget: ModalMolecularBudget,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless_field(
            profile,
            field,
            modes,
            elements,
            ModalEnd::Ports,
            budget.max_modes,
            &mut cancelled,
        )?
        .with_field_molecular_relaxation(law, molecular_budget, &mut cancelled)?;
        let (source, system) = image.into_parts();
        Self::validate_load_binding(profile, source.outlet_gas(), modes, &load)?;
        Self::prepare_loaded(source, system, load, duration, budget, &mut cancelled)
    }
}
