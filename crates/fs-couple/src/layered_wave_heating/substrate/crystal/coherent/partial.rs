//! Expected native cell heat from a partially coherent angular illumination.
//!
//! Coherent beam and eigenmode interference stays inside each original factor.
//! Only explicitly independent factors add powers, including factors at different
//! frequencies. Their physical amplitudes already contain any statistical weights.
//! The resulting mean cell watts drive the native deterministic heat solve. For a
//! nonlinear thermal body, this is not an ensemble-mean temperature prediction.

use super::{
    AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, BoundLayeredThermalBody,
    CoherentCrystalWavePower, CoherentCrystalWaveStep, CrystalBoxQuery, CrystalBoxSource,
    DepositionSimulation, LayerConductionConfig, LayeredThermalSolveControls,
    LayeredThermalWorkspace, LayeredWaveHeatingError, PreparedCrystalBoxCells,
    ResolvedMaterialStatePoint, SubstrateCellPower, SubstrateWaveHeatingLimits,
};
use fs_em::plane_wave::uniaxial::oriented::absorbing::
    AbsorbingCrystalPartiallyCoherentIllumination;

/// Cached finite-footprint absorption of independent coherent source factors.
///
/// Every factor retains its original material, frequency, angular fields and
/// phases. Preparation sums the native EM volume heat in every original depth
/// cell; six signed exterior face powers form the optical account. There is no
/// coherence between independent factors, no extra statistical weight, and no
/// conversion of face flux or balance residual into heat.
///
/// Each accepted interval reuses the cached watts. Material response, source
/// statistics and geometry remain held; there is no thermo-optic feedback or
/// lateral temperature resolution. Source and thermal ownership stay borrowed.
pub struct PreparedPartiallyCoherentCrystalHeating<'i, 'f, 'p, 's, 'm, 'b, 'a> {
    illumination: &'i AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm>,
    thermal: &'b BoundLayeredThermalBody<'a>,
    prepared: PreparedCrystalBoxCells,
}

impl<'i, 'f, 'p, 's, 'm, 'b, 'a>
    PreparedPartiallyCoherentCrystalHeating<'i, 'f, 'p, 's, 'm, 'b, 'a>
{
    /// Bind expected dielectric heat to the original planar thermal cells.
    ///
    /// The rectangle and all cell volumes must match the native geometry.
    /// Work counts cells * sum_f (2 * components_f)^2. The original EM owner
    /// separately bounds factor traversal. Storage includes this owner, cell
    /// powers, and peak EM-query or interval scratch, whichever is larger.
    /// Native thermal state/workspace and borrowed factors remain separate.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        illumination: &'i AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm>,
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

    pub const fn source(
        &self,
    ) -> &'i AbsorbingCrystalPartiallyCoherentIllumination<'f, 'p, 's, 'm> {
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

    /// Advance the original native body with expected optical heat added once.
    /// Optical face energies use the same accepted thermal duration. Existing
    /// authored heat remains intact; failures preserve the supplied state.
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

impl CrystalBoxSource for AbsorbingCrystalPartiallyCoherentIllumination<'_, '_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source().crystal_response().material()
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
