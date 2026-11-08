//! Coherent finite-crystal film heat on the original thermal timeline.
//!
//! All four internal Maxwell waves contribute to each cell's integrated
//! tensor work, including the exit-face reflection. The supplied original
//! optical laws and drive are held; capacity, conduction, authored heat and
//! atomic thermal acceptance remain owned by the native thermal body.
//! The planar source is uniform in transverse intensity. This adapter does
//! not average spatial interference between distinct transverse wavevectors.

use super::{
    BoundLayeredThermalBody, CrystalColumn, CrystalHeatingSource, DepositionSimulation,
    LayerConductionConfig, LayeredThermalSolveControls, LayeredThermalWorkspace,
    LayeredWaveHeatingError, ResolvedMaterialStatePoint, SubstrateCellPower,
    SubstrateWaveHeatingLimits, SubstrateWavePower, SubstrateWaveStep, advance_cells_to,
    prepare_column,
};
use fs_em::plane_wave::uniaxial::oriented::{
    absorbing::slab::PreparedAbsorbingSlabField, slab::SlabRegion,
};

/// Independent polarization and spectral factors through the finite film.
pub mod illumination;
pub use illumination::PreparedCrystalSlabIlluminationHeating;

impl CrystalHeatingSource for PreparedAbsorbingSlabField<'_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source()
            .entrance_interface()
            .crystal_response()
            .material()
    }

    fn factor_count(&self) -> usize {
        1
    }

    fn mode_pair_count(&self) -> Option<usize> {
        Some(16)
    }

    fn absorbed_w(
        &self,
        area: f64,
        start: f64,
        end: f64,
        _: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError> {
        Ok(self.slab_absorption(area, start, end)?.absorbed_power_w())
    }

    fn transmitted_flux(
        &self,
        depth: f64,
        _: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError> {
        Ok(self
            .field_at(SlabRegion::Crystal, [0.0, 0.0, depth])?
            .field
            .poynting_w_m2[2])
    }
}

/// Cached coherent film watts in original thermal-cell order.
///
/// Preparation integrates the complete internal standing wave in each
/// physical cell. Subsequent accepted thermal intervals reuse those watts
/// without repeating the Maxwell boundary solve or the depth integration.
/// Optical coefficients, thickness, optic axis and illumination remain held;
/// this owner does not infer thermo-optic feedback or carrier transients.
pub struct PreparedCrystalSlabHeating<'p, 'f, 's, 'm, 'b, 'a> {
    field: &'p PreparedAbsorbingSlabField<'f, 's, 'm>,
    column: CrystalColumn<'b, 'a>,
}

impl<'p, 'f, 's, 'm, 'b, 'a> PreparedCrystalSlabHeating<'p, 'f, 's, 'm, 'b, 'a> {
    /// Bind a finite crystal field to its original planar thermal column.
    ///
    /// Every conductivity and capacity law must retain the optical material
    /// card and fixed non-temperature query context. The thermal geometry
    /// supplies area, widths, subdivisions and cell identities; all depth
    /// cells must lie inside the original optical film. The work budget
    /// counts sixteen coherent mode pairs per cell and boundary query.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        field: &'p PreparedAbsorbingSlabField<'f, 's, 'm>,
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        let column = prepare_column(
            thermal,
            field,
            depth_start_m,
            limits,
            core::mem::size_of::<Self>(),
            cancelled,
        )?;
        Ok(Self { field, column })
    }

    pub const fn source(&self) -> &'p PreparedAbsorbingSlabField<'f, 's, 'm> {
        self.field
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

    /// Deposit the integrated optical watts exactly once over the accepted
    /// native interval. Refusal or cancellation leaves the initial state intact.
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
