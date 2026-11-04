//! Signed modal wall traces and original diffusion histories on one exact clock.
use super::*;
use fs_duct::modal::time_domain::{ModalViscothermalBudget, ModalViscothermalResolution};

impl ContinuousModalLine {
    /// Prepare the complete uniform-cylinder modal field with physical shear
    /// and thermal wall memories. Every wall trace couples the original radial
    /// basis; all acoustic and diffusion coordinates share the exact flow.
    ///
    /// The source retains the declared wall band and thin-layer domain. A
    /// time history does not establish accuracy outside that physical band.
    /// This excludes molecular attachment pending a shared-temperature model.
    ///
    /// # Errors
    /// Original geometry/gas, wall-domain, complete-capacity and exact-flow refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new_viscothermal(
        profile: &RadiusProfile,
        gas: &GasState,
        modes: ModeCount,
        elements: usize,
        end: ModalEnd,
        resolution: ModalViscothermalResolution,
        duration: f64,
        wall_budget: ModalViscothermalBudget,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless(
            profile,
            gas,
            modes,
            elements,
            end,
            budget.max_modes,
            &mut cancelled,
        )?
        .with_viscothermal(resolution, wall_budget, &mut cancelled)?;
        let (source, system) = image.into_parts();
        Self::prepare(source, system, None, duration, budget, &mut cancelled)
    }

    /// Interconnect original modal mouth-load memory with the full wall-loss
    /// operator before preparing the common physical clock. Wall and mouth
    /// heat remain part of the original passive storage/work account.
    ///
    /// # Errors
    /// The refusals of [`Self::new_viscothermal`], plus original load binding
    /// and total duct/diffusion/load capacity.
    #[allow(clippy::too_many_arguments)]
    pub fn new_viscothermal_loaded(
        profile: &RadiusProfile,
        gas: &GasState,
        modes: ModeCount,
        elements: usize,
        load: ModalAdmittanceLoad,
        resolution: ModalViscothermalResolution,
        duration: f64,
        wall_budget: ModalViscothermalBudget,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless(
            profile,
            gas,
            modes,
            elements,
            ModalEnd::Ports,
            budget.max_modes,
            &mut cancelled,
        )?
        .with_viscothermal(resolution, wall_budget, &mut cancelled)?;
        let (source, system) = image.into_parts();
        Self::validate_load_binding(profile, source.outlet_gas(), modes, &load)?;
        Self::prepare_loaded(source, system, load, duration, budget, &mut cancelled)
    }
}
