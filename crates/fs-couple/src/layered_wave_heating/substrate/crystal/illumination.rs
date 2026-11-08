//! Partially polarized and spectral crystal illumination on native thermal cells.
//!
//! Each factor is a coherent Maxwell field, including its ordinary/extraordinary
//! interference. Explicitly uncorrelated factors contribute additive expected
//! dielectric heat. Their physical peak electric fields already include any
//! square-root statistical weights; this adapter neither normalizes them nor
//! interprets discrete factors as a power spectral density.

use super::{
    BoundLayeredThermalBody, CrystalColumn, CrystalHeatingSource, DepositionSimulation,
    LayerConductionConfig, LayeredThermalSolveControls, LayeredThermalWorkspace,
    LayeredWaveHeatingError, OrientedCrystalFieldSide, ResolvedMaterialStatePoint,
    SubstrateCellPower, SubstrateWaveHeatingLimits, SubstrateWavePower, SubstrateWaveStep,
    advance_cells_to, prepare_column,
};
use fs_em::plane_wave::uniaxial::oriented::absorbing::{
    AbsorbingCrystalIllumination, CrystalIlluminationError,
};

impl From<CrystalIlluminationError> for LayeredWaveHeatingError {
    fn from(error: CrystalIlluminationError) -> Self {
        match error {
            CrystalIlluminationError::Cancelled => Self::Cancelled,
            CrystalIlluminationError::Electromagnetic(error) => error.into(),
        }
    }
}

impl CrystalHeatingSource for AbsorbingCrystalIllumination<'_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source().crystal_response().material()
    }

    fn factor_count(&self) -> usize {
        self.factor_count()
    }

    fn absorbed_w(
        &self,
        area: f64,
        start: f64,
        end: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError> {
        Ok(self
            .slab_absorption(area, start, end, cancelled)?
            .absorbed_power_w())
    }

    fn transmitted_flux(
        &self,
        depth: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError> {
        Ok(self.poynting_w_m2(
            OrientedCrystalFieldSide::Transmitted,
            [0.0, 0.0, depth],
            cancelled,
        )?[2])
    }
}

/// Cached expected crystal heat for a bounded, explicitly uncorrelated factor bank.
///
/// Preparation costs O(cells × factors), retaining coherent mode interference
/// within each factor. Each accepted native thermal interval then reuses the
/// summed cell watts without repeating electromagnetic solves or spectral sums.
/// The optical geometry, material context and illumination remain held. For a
/// nonlinear thermal law this is the trajectory driven by expected heat, not a
/// claim that temperature has been averaged over stochastic realizations.
pub struct PreparedCrystalIlluminationHeating<'p, 'f, 's, 'm, 'b, 'a> {
    illumination: &'p AbsorbingCrystalIllumination<'f, 's, 'm>,
    column: CrystalColumn<'b, 'a>,
}

impl<'p, 'f, 's, 'm, 'b, 'a> PreparedCrystalIlluminationHeating<'p, 'f, 's, 'm, 'b, 'a> {
    /// Bind the original illumination to a planar control column inside its
    /// homogeneous crystal half-space. The original thermal body owns all cell
    /// geometry, capacity, conduction, authored heat and accepted thermal state.
    /// Optical and thermal material cards and fixed query context must agree.
    ///
    /// The work budget counts four coherent modal-pair terms per factor for
    /// every cell and both boundary queries. An empty bank gives zero optical
    /// heat while retaining the ordinary native thermal evolution.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        illumination: &'p AbsorbingCrystalIllumination<'f, 's, 'm>,
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        let column = prepare_column(
            thermal,
            illumination,
            depth_start_m,
            limits,
            core::mem::size_of::<Self>(),
            cancelled,
        )?;
        Ok(Self {
            illumination,
            column,
        })
    }

    pub const fn source(&self) -> &'p AbsorbingCrystalIllumination<'f, 's, 'm> {
        self.illumination
    }

    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.column.thermal
    }

    pub fn cells(&self) -> &[SubstrateCellPower] {
        &self.column.cells
    }

    pub const fn power(&self) -> SubstrateWavePower {
        self.column.power
    }

    pub const fn preparation_mode_cell_pairs(&self) -> usize {
        self.column.mode_cell_pairs
    }

    pub const fn workspace_bytes(&self) -> usize {
        self.column.workspace_bytes
    }

    /// Deposit the summed cycle-mean optical watts once over the actual native
    /// thermal interval. Cancellation and failed solves leave `initial` intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SubstrateWaveStep<'t>, LayeredWaveHeatingError> {
        advance_cells_to(
            self.column.thermal,
            &self.column.cells,
            self.column.power,
            initial,
            target_s,
            conduction,
            controls,
            workspace,
            cancelled,
        )
    }
}
