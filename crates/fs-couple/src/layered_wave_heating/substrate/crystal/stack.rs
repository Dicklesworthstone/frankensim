//! Layer-resolved polarized optical heat on the original thermal timeline.
//!
//! Each optical layer maps to the corresponding original thermal layer.
//! Its cells receive analytic integrals of all four coherent Maxwell waves;
//! independently oriented material tensors never become an averaged crystal.
//! The existing thermal body still owns capacity, conduction, authored heat,
//! boundaries and atomic acceptance. Optical coefficients and drive are held.

use super::{
    BoundLayeredThermalBody, CellId, CrystalColumn, DepositionSimulation, LayerConductionConfig,
    LayeredThermalSolveControls, LayeredThermalWorkspace, LayeredWaveHeatingError,
    SubstrateCellPower, SubstrateWaveHeatingLimits, SubstrateWavePower, SubstrateWaveStep,
    admit_crystal_materials, advance_cells_to, invalid, poll, reserve,
};
use fs_em::plane_wave::uniaxial::oriented::absorbing::stack::{
    PreparedAbsorbingStackField, StackRegion,
};
use fs_math::dd::Dd;

/// Cached coherent multilayer heat in native thermal-cell order.
/// Preparation is O(cells + layers). Accepted thermal intervals reuse cell
/// watts without a Maxwell solve, material query or depth quadrature.
pub struct PreparedCrystalStackHeating<'p, 'f, 's, 'm, 'b, 'a> {
    field: &'p PreparedAbsorbingStackField<'f, 's, 'm>,
    column: CrystalColumn<'b, 'a>,
}

impl<'p, 'f, 's, 'm, 'b, 'a> PreparedCrystalStackHeating<'p, 'f, 's, 'm, 'b, 'a> {
    /// Bind the complete original optical stack to a planar thermal body.
    /// Ordered layer counts, exact authored thicknesses, original material
    /// cards and fixed non-temperature query context must match. Each thermal
    /// layer may contain any admitted native cell subdivision.
    ///
    /// Temperature evolves on the native timeline; the supplied optical
    /// material snapshots remain held. Temperature-dependent optical feedback
    /// requires refreshing every affected layer and preparing a new field.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        field: &'p PreparedAbsorbingStackField<'f, 's, 'm>,
        limits: SubstrateWaveHeatingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let body = thermal.prepared_body();
        let optical = field.source().layers();
        let thermal_layers = body.conduction().layers();
        let count = body.cell_count();
        if !body.conduction().geometry().is_planar()
            || count == 0
            || count > limits.maximum_cells
            || optical.len() != thermal_layers.len()
        {
            return Err(invalid(
                "crystal stack heat requires matching ordered planar thermal layers",
            ));
        }
        for (source, layer) in optical.iter().zip(thermal_layers) {
            poll(&mut cancelled)?;
            if source.thickness_m() != layer.thickness_m() {
                return Err(invalid(
                    "crystal optical and thermal layer thicknesses must match exactly",
                ));
            }
        }
        let mode_cell_pairs = count
            .checked_add(2)
            .and_then(|n| n.checked_mul(16))
            .ok_or(invalid("crystal stack cell-mode work overflows"))?;
        let per_cell = core::mem::size_of::<SubstrateCellPower>()
            + core::mem::size_of::<CellId>()
            + core::mem::size_of::<f64>()
            + core::mem::size_of::<(CellId, f64)>();
        let workspace_bytes = count
            .checked_mul(per_cell)
            .and_then(|n| n.checked_add(core::mem::size_of::<Self>()))
            .and_then(|n| n.checked_add(core::mem::size_of::<SubstrateWaveStep<'_>>()))
            .ok_or(invalid("crystal stack heat storage overflows"))?;
        if mode_cell_pairs > limits.maximum_mode_cell_pairs
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(invalid(
                "crystal stack heat exceeds its work or storage budget",
            ));
        }
        admit_crystal_materials(
            thermal,
            |index| {
                optical
                    .get(index)
                    .map(|layer| layer.crystal_response().material())
                    .ok_or(invalid("unknown original crystal stack material"))
            },
            &mut cancelled,
        )?;

        let area = body.conduction().area_m2();
        let ids = thermal.cell_ids();
        let mut cells = reserve(count)?;
        let mut absorbed = Dd::ZERO;
        for (index, (layer, &n)) in field
            .response()
            .layers
            .iter()
            .zip(body.layer_cell_counts())
            .enumerate()
        {
            poll(&mut cancelled)?;
            let width = Dd::from_f64(layer.thickness_m);
            let entrance = Dd::from_f64(layer.entrance_z_m);
            let mut start = 0.0;
            let mut global_start = layer.entrance_z_m;
            for cell in 0..n {
                poll(&mut cancelled)?;
                let local_end = if cell + 1 == n {
                    width
                } else {
                    width * Dd::from_f64((cell + 1) as f64) / Dd::from_f64(n as f64)
                };
                let end = local_end.to_f64();
                let global_end = if cell + 1 == n {
                    layer.exit_z_m
                } else {
                    (entrance + local_end).to_f64()
                };
                if !global_end.is_finite() || global_end <= global_start || end <= start {
                    return Err(invalid(
                        "crystal stack cell depth interval is unrepresentable",
                    ));
                }
                // Integrate LOCAL coordinates directly. Subtracting two large
                // global positions would erase phase and thin-cell precision.
                let watts = field
                    .layer_absorption(index, area, start, end)?
                    .absorbed_power_w();
                if !watts.is_finite() || watts < 0.0 {
                    return Err(invalid("crystal stack cell watts are unrepresentable"));
                }
                cells.push(SubstrateCellPower {
                    cell_id: ids[cells.len()],
                    depth_m: [global_start, global_end],
                    absorbed_w: watts,
                });
                absorbed = absorbed + Dd::from_f64(watts);
                start = end;
                global_start = global_end;
            }
        }
        let mut boundary = |region, depth| -> Result<f64, LayeredWaveHeatingError> {
            poll(&mut cancelled)?;
            let flux = field
                .field_at(region, [0.0, 0.0, depth])?
                .field
                .poynting_w_m2[2];
            let watts = area * flux;
            if !watts.is_finite() || (flux != 0.0 && watts == 0.0) {
                return Err(invalid("crystal stack boundary watts are unrepresentable"));
            }
            Ok(watts)
        };
        let top = boundary(StackRegion::TotalIncidentSide, 0.0)?;
        let bottom = boundary(StackRegion::Transmitted, field.source().thickness_m())?;
        let entering = Dd::from_f64(top.max(0.0)) + Dd::from_f64((-bottom).max(0.0));
        let leaving = Dd::from_f64((-top).max(0.0)) + Dd::from_f64(bottom.max(0.0));
        let power = SubstrateWavePower {
            entering_w: entering.to_f64(),
            leaving_w: leaving.to_f64(),
            absorbed_w: absorbed.to_f64(),
            residual_w: (entering - leaving - absorbed).to_f64(),
        };
        if [
            power.entering_w,
            power.leaving_w,
            power.absorbed_w,
            power.residual_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("crystal stack power account is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            field,
            column: CrystalColumn {
                thermal,
                cells,
                power,
                mode_cell_pairs,
                workspace_bytes,
            },
        })
    }

    pub const fn source(&self) -> &'p PreparedAbsorbingStackField<'f, 's, 'm> {
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

    /// Deposit these held optical watts once over the accepted native step.
    /// A failed gate or cancellation leaves the supplied state unchanged.
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
