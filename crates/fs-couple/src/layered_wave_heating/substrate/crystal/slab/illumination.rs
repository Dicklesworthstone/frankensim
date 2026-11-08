//! Independent polarization and spectral film factors on native thermal cells.
//!
//! Each admitted group contains at most one prepared incident plane wave.
//! Its four internal eigenwaves remain coherent, including exit reflections,
//! while independent factors contribute additive expected tensor heat. The
//! physical peak fields already contain their square-root statistical weights.
//! Coherent multi-beam groups refuse: a point sample times an area cannot
//! integrate their possible transverse interference fringes.

use super::super::{
    BoundLayeredThermalBody, CrystalColumn, CrystalHeatingSource, DepositionSimulation,
    LayerConductionConfig, LayeredThermalSolveControls, LayeredThermalWorkspace,
    LayeredWaveHeatingError, ResolvedMaterialStatePoint, SubstrateCellPower,
    SubstrateWaveHeatingLimits, SubstrateWavePower, SubstrateWaveStep, advance_cells_to, invalid,
    poll, prepare_column,
};
use fs_em::plane_wave::uniaxial::oriented::{
    absorbing::slab::AbsorbingSlabIllumination, slab::SlabRegion,
};
use fs_math::dd::Dd;

impl CrystalHeatingSource for AbsorbingSlabIllumination<'_, '_, '_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source()
            .entrance_interface()
            .crystal_response()
            .material()
    }

    fn factor_count(&self) -> usize {
        self.component_count()
    }

    fn mode_pair_count(&self) -> Option<usize> {
        self.component_count().checked_mul(16)
    }

    fn absorbed_w(
        &self,
        area: f64,
        start: f64,
        end: f64,
        mut cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        // prepare_column owns positive area and ordered finite cell depths.
        // The original optical extent still applies when every group is dark.
        if end > self.source().thickness_m() {
            return Err(invalid(
                "film illumination thermal cells extend beyond the original crystal",
            ));
        }
        let mut sum = Dd::ZERO;
        for group in self.groups() {
            poll(&mut cancelled)?;
            match group.components() {
                [] => {}
                [field] => {
                    sum = sum
                        + Dd::from_f64(field.slab_absorption(area, start, end)?.absorbed_power_w());
                }
                _ => {
                    return Err(invalid(
                        "planar film heating requires at most one plane wave per independent group",
                    ));
                }
            }
        }
        let watts = sum.to_f64();
        if !watts.is_finite() || watts < 0.0 {
            return Err(invalid(
                "film illumination absorbed watts are unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(watts)
    }

    fn transmitted_flux(
        &self,
        depth: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError> {
        Ok(self.poynting_w_m2(SlabRegion::Crystal, [0.0, 0.0, depth], cancelled)?[2])
    }
}

/// Cached expected heat from independent finite-film polarization factors.
///
/// Original material/card identity, fixed non-temperature context, film
/// extent and thermal cell identities remain authoritative. Optical laws,
/// geometry and illumination are held. Capacity, conduction, authored heat
/// and atomic acceptance remain owned by the native thermal body. For a
/// nonlinear thermal law this is the trajectory driven by expected heat,
/// not an average temperature over random optical realizations.
pub struct PreparedCrystalSlabIlluminationHeating<'i, 'g, 'p, 'f, 's, 'm, 'b, 'a> {
    illumination: &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm>,
    column: CrystalColumn<'b, 'a>,
}

impl<'i, 'g, 'p, 'f, 's, 'm, 'b, 'a>
    PreparedCrystalSlabIlluminationHeating<'i, 'g, 'p, 'f, 's, 'm, 'b, 'a>
{
    /// Bind explicitly independent plane-wave groups to the original planar
    /// thermal volume. Incidence direction, polarization and frequency may
    /// differ across groups. Each field already includes its statistical weight.
    ///
    /// Multiple coherent incident beams in one group require an authored
    /// transverse volume integral and refuse here, even if a particular point
    /// is dark. Work counts sixteen modal pairs per component for every cell
    /// and both boundaries. Empty illumination retains native thermal evolution.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        illumination: &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm>,
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        for group in illumination.groups() {
            poll(&mut cancelled)?;
            if group.component_count() > 1 {
                return Err(invalid(
                    "planar film heating requires at most one plane wave per independent group",
                ));
            }
        }
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

    pub const fn source(&self) -> &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm> {
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

    /// Deposit summed cycle-mean watts once over the accepted native interval.
    /// Cancellation and failed solves leave the supplied simulation unchanged.
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
