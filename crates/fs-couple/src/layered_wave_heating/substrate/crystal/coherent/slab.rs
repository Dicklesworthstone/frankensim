//! Finite-footprint heat from original coherent and partially coherent films.
//!
//! Every beam retains all four internal Maxwell waves, both authored optical
//! faces and any coherent drives from those faces. Only volume-integrated
//! dielectric work heats the original thermal cells. Six-face optical flow
//! remains an independent account, including lateral transport.
//!
//! The shared cell cache owns geometry/material admission, work and storage
//! bounds, and the single thermal debit. These consumers do not reconstruct
//! an interface, substitute a half-space, or repeat optical integration during
//! thermal advancement. Optics and source statistics remain held.

use super::{
    AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, BoundLayeredThermalBody,
    CoherentCrystalWavePower, CoherentCrystalWaveStep, CrystalBoxQuery, CrystalBoxSource,
    DepositionSimulation, LayerConductionConfig, LayeredThermalSolveControls,
    LayeredThermalWorkspace, LayeredWaveHeatingError, PreparedCrystalBoxCells,
    ResolvedMaterialStatePoint, SubstrateCellPower, SubstrateWaveHeatingLimits,
};
use fs_em::plane_wave::uniaxial::oriented::absorbing::slab::{
    AbsorbingSlabCoherentIllumination, AbsorbingSlabIllumination,
};

/// Cached native thermal-cell heat from one coherent angular film illumination.
///
/// Every original depth cell lies inside the same optical film. The footprint
/// must match the native planar area; material card, fixed query axes, capacity
/// and conduction remain those of the supplied body. Four modes per beam retain
/// angular fringes and internal standing waves before the dielectric loss square.
/// No lateral thermal resolution or thermo-optic feedback is inferred.
pub struct PreparedCoherentCrystalSlabHeating<'i, 'p, 'f, 's, 'm, 'b, 'a> {
    illumination: &'i AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm>,
    thermal: &'b BoundLayeredThermalBody<'a>,
    prepared: PreparedCrystalBoxCells,
}

impl<'i, 'p, 'f, 's, 'm, 'b, 'a>
    PreparedCoherentCrystalSlabHeating<'i, 'p, 'f, 's, 'm, 'b, 'a>
{
    /// Integrate the original film's physical boxes before publishing any cache.
    /// Work is cells * (4 * components)^2. Storage covers this owner, cell powers,
    /// and peak optical-query or interval scratch. Native thermal workspace and
    /// borrowed material/source owners remain separate. Invalid film extent,
    /// geometry, material, budget or cancellation publishes no candidate.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        illumination: &'i AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm>,
        footprint_origin_m: [f64; 2],
        footprint_extent_m: [f64; 2],
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        let prepared = PreparedCrystalBoxCells::try_new(
            thermal, illumination, footprint_origin_m, footprint_extent_m,
            depth_start_m, limits, core::mem::size_of::<Self>(), cancelled,
        )?;
        Ok(Self { illumination, thermal, prepared })
    }

    pub const fn source(&self) -> &'i AbsorbingSlabCoherentIllumination<'p, 'f, 's, 'm> {
        self.illumination
    }

    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }

    pub const fn region(&self) -> &AbsorbingCrystalBox {
        &self.prepared.region
    }

    pub fn cells(&self) -> &[SubstrateCellPower] {
        &self.prepared.cells
    }

    pub const fn power(&self) -> CoherentCrystalWavePower {
        self.prepared.power
    }

    pub const fn preparation_mode_cell_pairs(&self) -> usize {
        self.prepared.mode_cell_pairs
    }

    pub const fn workspace_bytes(&self) -> usize {
        self.prepared.workspace_bytes
    }

    /// Add cached dielectric heat once over the accepted native thermal interval.
    /// Six signed face energies use that same duration and are not added as heat.
    /// Failure or cancellation leaves the supplied initial state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCrystalWaveStep<'t>, LayeredWaveHeatingError> {
        self.prepared.advance_to(
            self.thermal, initial, target_s, conduction, controls, workspace, cancelled,
        )
    }
}

/// Expected native cell heat from independent coherent angular/spectral groups.
///
/// All four waves interfere within each original beam and coherent group;
/// independent groups add their positive dielectric powers. Physical source
/// amplitudes already contain statistical weights. This drives a deterministic
/// heat solve with expected optical watts; for a nonlinear thermal body it is
/// not an ensemble-mean temperature prediction. Film boundaries and material
/// ownership remain exact even for an empty illumination.
pub struct PreparedPartiallyCoherentCrystalSlabHeating<'i, 'g, 'p, 'f, 's, 'm, 'b, 'a> {
    illumination: &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm>,
    thermal: &'b BoundLayeredThermalBody<'a>,
    prepared: PreparedCrystalBoxCells,
}

impl<'i, 'g, 'p, 'f, 's, 'm, 'b, 'a>
    PreparedPartiallyCoherentCrystalSlabHeating<'i, 'g, 'p, 'f, 's, 'm, 'b, 'a>
{
    /// Bind every cell to an original finite-film box and sum native group heat.
    /// Work is cells * sum_g (4 * components_g)^2; the native optical owner also
    /// bounds group traversal. Geometry, material, memory and atomic acceptance
    /// use exactly the same shared cache as the coherent consumer.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        illumination: &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm>,
        footprint_origin_m: [f64; 2],
        footprint_extent_m: [f64; 2],
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        let prepared = PreparedCrystalBoxCells::try_new(
            thermal, illumination, footprint_origin_m, footprint_extent_m,
            depth_start_m, limits, core::mem::size_of::<Self>(), cancelled,
        )?;
        Ok(Self { illumination, thermal, prepared })
    }

    pub const fn source(&self) -> &'i AbsorbingSlabIllumination<'g, 'p, 'f, 's, 'm> {
        self.illumination
    }

    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }

    pub const fn region(&self) -> &AbsorbingCrystalBox {
        &self.prepared.region
    }

    pub fn cells(&self) -> &[SubstrateCellPower] {
        &self.prepared.cells
    }

    pub const fn power(&self) -> CoherentCrystalWavePower {
        self.prepared.power
    }

    pub const fn preparation_mode_cell_pairs(&self) -> usize {
        self.prepared.mode_cell_pairs
    }

    pub const fn workspace_bytes(&self) -> usize {
        self.prepared.workspace_bytes
    }

    /// Advance the native body with each expected cell watt added exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCrystalWaveStep<'t>, LayeredWaveHeatingError> {
        self.prepared.advance_to(
            self.thermal, initial, target_s, conduction, controls, workspace, cancelled,
        )
    }
}

impl CrystalBoxSource for AbsorbingSlabCoherentIllumination<'_, '_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source().entrance_interface().crystal_response().material()
    }

    fn mode_pairs(&self) -> Result<usize, LayeredWaveHeatingError> {
        Ok(self.mode_pair_count())
    }

    fn integrated_power(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CrystalBoxQuery, LayeredWaveHeatingError> {
        let result = self.box_absorption(region, limits, cancelled)?;
        Ok(CrystalBoxQuery {
            absorbed_w: result.absorbed_power_w(),
            outward_face_w: result.outward_face_power_w(),
            workspace_bytes: result.workspace_bytes(),
        })
    }
}

impl CrystalBoxSource for AbsorbingSlabIllumination<'_, '_, '_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source().entrance_interface().crystal_response().material()
    }

    fn mode_pairs(&self) -> Result<usize, LayeredWaveHeatingError> {
        Ok(self.mode_pair_count())
    }

    fn integrated_power(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CrystalBoxQuery, LayeredWaveHeatingError> {
        let result = self.box_absorption(region, limits, cancelled)?;
        Ok(CrystalBoxQuery {
            absorbed_w: result.absorbed_power_w(),
            outward_face_w: result.outward_face_power_w(),
            workspace_bytes: result.workspace_bytes(),
        })
    }
}
