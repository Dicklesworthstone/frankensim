//! Polarization-dependent crystal absorption on the native thermal timeline.
//!
//! The thermal body is a finite planar control column inside the original
//! homogeneous uniaxial crystal. The base adapter uses a half-space; the slab
//! adapter retains both optical faces of a finite absorbing film.
//! Actual area and depth cells receive analytic integrals of
//! the original coherent ordinary/extraordinary field's dielectric loss.
//! Interference is retained before heating; summing independent mode powers
//! or sampling heat only at cell centres would describe a different source.
//!
//! Geometry, optic axis, optical coefficients and monochromatic illumination
//! are held. The existing thermal owner supplies capacity, conduction, authored
//! heat, boundaries and atomic acceptance. This is cycle-mean heating, without
//! optical feedback, carrier transients, phase change or lateral temperature
//! gradients. Equal material cards retain provenance, not measured validation.

use super::{
    BoundLayeredThermalBody, CellId, DepositionSimulation, LayerConductionConfig,
    LayeredThermalSolveControls, LayeredThermalWorkspace, LayeredWaveHeatingError,
    SubstrateCellPower, SubstrateWaveHeatingLimits, SubstrateWavePower, SubstrateWaveStep,
    admit_context, advance_cells_to, invalid, poll, reserve,
};
use fs_em::plane_wave::uniaxial::oriented::{
    OrientedCrystalFieldSide, absorbing::PreparedAbsorbingCrystalField,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_math::dd::Dd;

/// Independent polarization and spectral factors on the same thermal owner.
pub mod illumination;
pub use illumination::PreparedCrystalIlluminationHeating;

/// Coherent finite-footprint crystal heat, including lateral optical flux.
pub mod coherent;
pub use coherent::{
    CoherentCrystalWaveEnergy, CoherentCrystalWavePower, CoherentCrystalWaveStep,
    PreparedCoherentCrystalHeating, PreparedPartiallyCoherentCrystalHeating,
    PreparedCoherentCrystalSlabHeating, PreparedPartiallyCoherentCrystalSlabHeating,
};

/// Finite absorbing crystal films on the same thermal-cell owner.
pub mod slab;
pub use slab::{PreparedCrystalSlabHeating, PreparedCrystalSlabIlluminationHeating};

/// Independently oriented absorbing layers on their original thermal cells.
pub mod stack;
pub use stack::PreparedCrystalStackHeating;

// Only native electromagnetic owners implement this internal source seam.
// Geometry, material admission, cell integration and thermal accounting stay
// identical for a coherent field and an explicitly uncorrelated illumination.
trait CrystalHeatingSource {
    fn material(&self) -> &ResolvedMaterialStatePoint;
    fn factor_count(&self) -> usize;
    fn mode_pair_count(&self) -> Option<usize> {
        self.factor_count().checked_mul(4)
    }
    fn absorbed_w(
        &self,
        area: f64,
        start: f64,
        end: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError>;
    fn transmitted_flux(
        &self,
        depth: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, LayeredWaveHeatingError>;
}

impl CrystalHeatingSource for PreparedAbsorbingCrystalField<'_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source().crystal_response().material()
    }
    fn factor_count(&self) -> usize {
        1
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
            .field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, depth])?
            .field
            .poynting_w_m2[2])
    }
}

struct CrystalColumn<'b, 'a> {
    thermal: &'b BoundLayeredThermalBody<'a>,
    cells: Vec<SubstrateCellPower>,
    power: SubstrateWavePower,
    mode_cell_pairs: usize,
    workspace_bytes: usize,
}

/// Cached coherent crystal watts in the original thermal-cell order.
/// Preparation is O(cells) with four modal-pair terms per interval. Accepted
/// thermal intervals reuse the same watts without solving Maxwell again.
pub struct PreparedCrystalWaveHeating<'p, 's, 'm, 'b, 'a> {
    field: &'p PreparedAbsorbingCrystalField<'s, 'm>,
    thermal: &'b BoundLayeredThermalBody<'a>,
    cells: Vec<SubstrateCellPower>,
    power: SubstrateWavePower,
    mode_cell_pairs: usize,
    workspace_bytes: usize,
}

impl<'p, 's, 'm, 'b, 'a> PreparedCrystalWaveHeating<'p, 's, 'm, 'b, 'a> {
    /// Bind the original crystal field to one planar thermal column. All
    /// conductivity and capacity sources must retain the crystal's material
    /// card and shared fixed query context. Optical temperature stays explicit
    /// and held, as in the existing substrate adapter.
    ///
    /// Depth is measured from the original optical interface, z=0. Area,
    /// widths, subdivision and cell IDs come only from the thermal body. The
    /// work budget counts four coherent mode pairs per cell and boundary query.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        field: &'p PreparedAbsorbingCrystalField<'s, 'm>,
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
        Ok(Self {
            field,
            thermal,
            cells: column.cells,
            power: column.power,
            mode_cell_pairs: column.mode_cell_pairs,
            workspace_bytes: column.workspace_bytes,
        })
    }

    pub const fn source(&self) -> &'p PreparedAbsorbingCrystalField<'s, 'm> {
        self.field
    }
    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }
    pub fn cells(&self) -> &[SubstrateCellPower] {
        &self.cells
    }
    pub const fn power(&self) -> SubstrateWavePower {
        self.power
    }
    pub const fn preparation_mode_cell_pairs(&self) -> usize {
        self.mode_cell_pairs
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Deposit these held cycle-mean watts once over the actual thermal
    /// interval, sharing the original nonlinear thermal solve and heat account.
    /// All failure/cancellation paths leave the supplied simulation unchanged.
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
            self.thermal,
            &self.cells,
            self.power,
            initial,
            target_s,
            conduction,
            controls,
            workspace,
            cancelled,
        )
    }
}

fn prepare_column<'b, 'a>(
    thermal: &'b BoundLayeredThermalBody<'a>,
    source: &impl CrystalHeatingSource,
    depth_start_m: f64,
    limits: SubstrateWaveHeatingLimits,
    owner_bytes: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CrystalColumn<'b, 'a>, LayeredWaveHeatingError> {
    poll(&mut cancelled)?;
    let body = thermal.prepared_body();
    let count = body.cell_count();
    if !depth_start_m.is_finite()
        || depth_start_m < 0.0
        || !body.conduction().geometry().is_planar()
        || count == 0
        || count > limits.maximum_cells
    {
        return Err(invalid(
            "crystal heat requires finite nonnegative depth and admitted planar thermal cells",
        ));
    }
    let mode_pairs = source
        .mode_pair_count()
        .ok_or(invalid("crystal mode-pair work overflows"))?;
    let mode_cell_pairs = count
        .checked_add(2)
        .and_then(|n| n.checked_mul(mode_pairs))
        .ok_or(invalid("crystal cell-mode work overflows"))?;
    let per_cell = core::mem::size_of::<SubstrateCellPower>()
        + core::mem::size_of::<CellId>()
        + core::mem::size_of::<f64>()
        + core::mem::size_of::<(CellId, f64)>();
    let workspace_bytes = count
        .checked_mul(per_cell)
        .and_then(|n| n.checked_add(owner_bytes))
        .and_then(|n| n.checked_add(core::mem::size_of::<SubstrateWaveStep<'_>>()))
        .ok_or(invalid("crystal heat storage overflows"))?;
    if mode_cell_pairs > limits.maximum_mode_cell_pairs
        || workspace_bytes > limits.maximum_workspace_bytes
    {
        return Err(invalid("crystal heat exceeds its work or storage budget"));
    }
    admit_crystal_material(thermal, source.material(), &mut cancelled)?;
    let area = body.conduction().area_m2();
    let ids = thermal.cell_ids();
    let mut cells = reserve(count)?;
    let mut layer_start = Dd::from_f64(depth_start_m);
    let mut left = depth_start_m;
    let mut absorbed = Dd::ZERO;
    for (layer, &n) in body
        .conduction()
        .layers()
        .iter()
        .zip(body.layer_cell_counts())
    {
        let width = Dd::from_f64(layer.thickness_m());
        for local in 0..n {
            poll(&mut cancelled)?;
            let right = (layer_start
                + width * Dd::from_f64((local + 1) as f64) / Dd::from_f64(n as f64))
            .to_f64();
            if !right.is_finite() || right <= left {
                return Err(invalid("crystal cell depth interval is unrepresentable"));
            }
            let watts = source.absorbed_w(area, left, right, &mut cancelled)?;
            if !watts.is_finite() || watts < 0.0 {
                return Err(invalid("crystal cell watts are unrepresentable"));
            }
            cells.push(SubstrateCellPower {
                cell_id: ids[cells.len()],
                depth_m: [left, right],
                absorbed_w: watts,
            });
            absorbed = absorbed + Dd::from_f64(watts);
            left = right;
        }
        layer_start = layer_start + width;
    }
    let mut boundary_power = |depth| -> Result<f64, LayeredWaveHeatingError> {
        let flux = source.transmitted_flux(depth, &mut cancelled)?;
        let watts = area * flux;
        if !watts.is_finite() || (flux != 0.0 && watts == 0.0) {
            return Err(invalid("crystal boundary watts are unrepresentable"));
        }
        Ok(watts)
    };
    // Positive z points into the top face and out of the bottom face. A
    // coherently driven face contributes its total net Poynting flux.
    let top_w = boundary_power(depth_start_m)?;
    let bottom_w = boundary_power(left)?;
    let entering = Dd::from_f64(top_w.max(0.0)) + Dd::from_f64((-bottom_w).max(0.0));
    let leaving = Dd::from_f64((-top_w).max(0.0)) + Dd::from_f64(bottom_w.max(0.0));
    let entering_w = entering.to_f64();
    let leaving_w = leaving.to_f64();
    let absorbed_w = absorbed.to_f64();
    let residual_w = (entering - leaving - absorbed).to_f64();
    if !entering_w.is_finite()
        || !leaving_w.is_finite()
        || !absorbed_w.is_finite()
        || !residual_w.is_finite()
    {
        return Err(invalid("crystal power account is unrepresentable"));
    }
    poll(&mut cancelled)?;
    Ok(CrystalColumn {
        thermal,
        cells,
        power: SubstrateWavePower {
            entering_w,
            leaving_w,
            absorbed_w,
            residual_w,
        },
        mode_cell_pairs,
        workspace_bytes,
    })
}

fn admit_crystal_material(
    thermal: &BoundLayeredThermalBody<'_>,
    original: &ResolvedMaterialStatePoint,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LayeredWaveHeatingError> {
    admit_crystal_materials(thermal, |_| Ok(original), cancelled)
}

fn admit_crystal_materials<'m>(
    thermal: &BoundLayeredThermalBody<'_>,
    mut original: impl FnMut(usize) -> Result<&'m ResolvedMaterialStatePoint, LayeredWaveHeatingError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LayeredWaveHeatingError> {
    let body = thermal.prepared_body();
    for (index, (layer, storage)) in body
        .conduction()
        .layers()
        .iter()
        .zip(body.storages())
        .enumerate()
    {
        poll(cancelled)?;
        let original = original(index)?;
        if original.card_identity() != layer.material_card_identity() {
            return Err(invalid(
                "crystal conductivity must retain the original optical material card",
            ));
        }
        admit_context(
            original,
            &layer.query_endpoints()[0],
            layer.temperature_axis(),
        )?;
        for component in storage.components() {
            poll(cancelled)?;
            if original.card_identity() != component.material_card_identity() {
                return Err(invalid(
                    "crystal capacity must retain the original optical material card",
                ));
            }
            admit_context(
                original,
                &component.query_endpoints()[0],
                component.temperature_axis(),
            )?;
        }
    }
    Ok(())
}
