//! Coherent crystal absorption over a real rectangular thermal footprint.
//!
//! Every original depth cell receives the volume integral of the complete
//! transmitted field. Crossed-beam fringes and ordinary/extraordinary cross
//! terms enter the original dielectric loss before heat is deposited. The
//! optical account includes all six faces; lateral Poynting flux is not heat.
//!
//! The native thermal model has one temperature per cross section and resolves
//! conduction through depth. It does not resolve lateral temperature gradients,
//! lateral thermal conduction, patch-edge scattering, thermo-optic feedback or
//! carrier transients. The rectangular support selects a control volume inside
//! the original homogeneous crystal; it adds no optical boundary. Finite-film
//! consumers retain both original optical faces and admit only boxes inside them.

use super::{
    BoundLayeredThermalBody, CellId, DepositionSimulation, LayerConductionConfig,
    LayeredThermalSolveControls, LayeredThermalWorkspace, LayeredWaveHeatingError,
    SubstrateCellPower, SubstrateWaveHeatingLimits, admit_crystal_material, invalid, poll,
    reserve, ResolvedMaterialStatePoint,
};
use super::super::advance_cell_heat_to;
use fs_em::plane_wave::uniaxial::oriented::absorbing::{
    AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, AbsorbingCrystalCoherentIllumination,
};
use fs_material::thermal::transport::stack::transient::implicit::LayeredThermalSolveReport;
use fs_math::dd::Dd;

/// Signed mean optical flow through the actual finite control volume.
/// Independent coherent factors contribute their expected powers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentCrystalWavePower {
    /// Outward watts at [x_min, x_max, y_min, y_max, z_min, z_max].
    /// A negative value is net entry through that face. These are face
    /// integrals, not a split of incoming and outgoing rays within each face.
    pub outward_face_w: [f64; 6],
    /// Sum of the original thermal cells' volume-integrated dielectric heat.
    pub absorbed_w: f64,
    /// Minus the outward face sum minus heat; never reassigned to cell heat.
    pub residual_w: f64,
}

/// Original six-face optical account over the accepted native thermal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentCrystalWaveEnergy {
    pub outward_face_j: [f64; 6],
    /// Sum of the actual cell-energy additions supplied once to the heat owner.
    pub absorbed_j: f64,
    pub residual_j: f64,
}

/// One atomic native thermal candidate and the source's finite-volume account.
pub struct CoherentCrystalWaveStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    pub heat_report: LayeredThermalSolveReport,
    pub optical: CoherentCrystalWaveEnergy,
    pub cell_energy_j: Vec<(CellId, f64)>,
}

/// Expected heat from independent coherent angular/spectral factors.
pub mod partial;
pub use partial::PreparedPartiallyCoherentCrystalHeating;

/// Angular interference and independent coherent groups in finite films.
pub mod slab;
pub use slab::{PreparedCoherentCrystalSlabHeating, PreparedPartiallyCoherentCrystalSlabHeating};

/// Cached coherent cell watts from the original crystal, illumination and body.
///
/// Preparation integrates all (2 * components)^2 transmitted-mode pairs in
/// every native cell. Repeated thermal intervals reuse those watts without
/// another Maxwell solve or interference integration. Optical material state,
/// geometry and source amplitudes remain held over the thermal evolution.
pub struct PreparedCoherentCrystalHeating<'p, 'f, 's, 'm, 'b, 'a> {
    illumination: &'p AbsorbingCrystalCoherentIllumination<'f, 's, 'm>,
    thermal: &'b BoundLayeredThermalBody<'a>,
    prepared: PreparedCrystalBoxCells,
}

impl<'p, 'f, 's, 'm, 'b, 'a> PreparedCoherentCrystalHeating<'p, 'f, 's, 'm, 'b, 'a> {
    /// Bind the physical x/y rectangle to the original planar thermal area.
    ///
    /// Depth subdivision, conductivity, capacity, material cards and fixed
    /// query axes belong to the supplied native thermal body. The rectangle's
    /// area and each integrated volume must match that body's geometry.
    ///
    /// The work budget counts every modal pair in every cell; each pair also
    /// supplies the six face integrals. Storage includes this owner, its cell
    /// powers, and the larger of EM query scratch or interval-adapter scratch.
    /// Borrowed source/body, native thermal state/workspace and allocator
    /// overhead remain separate, as in the existing substrate adapters.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        illumination: &'p AbsorbingCrystalCoherentIllumination<'f, 's, 'm>,
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

    pub const fn source(&self) -> &'p AbsorbingCrystalCoherentIllumination<'f, 's, 'm> {
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

    /// Add cached cell heat once on the original accepted thermal clock.
    /// Face energies remain an optical account, separate from dielectric heat.
    /// Failure/cancellation preserves the supplied initial state.
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

// This private seam only accepts original EM owners. Every public adapter uses
// the same native cell geometry, material admission and single thermal debit.
trait CrystalBoxSource {
    fn material(&self) -> &ResolvedMaterialStatePoint;
    fn mode_pairs(&self) -> Result<usize, LayeredWaveHeatingError>;
    fn integrated_power(
        &self,
        region: AbsorbingCrystalBox,
        limits: AbsorbingCrystalBoxLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CrystalBoxQuery, LayeredWaveHeatingError>;
}

struct CrystalBoxQuery {
    absorbed_w: f64,
    outward_face_w: [f64; 6],
    workspace_bytes: usize,
}

impl CrystalBoxSource for AbsorbingCrystalCoherentIllumination<'_, '_, '_> {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.source().crystal_response().material()
    }

    fn mode_pairs(&self) -> Result<usize, LayeredWaveHeatingError> {
        self.component_count().checked_mul(2)
            .and_then(|n| n.checked_mul(n))
            .ok_or(invalid("coherent crystal modal-pair work overflows"))
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

struct PreparedCrystalBoxCells {
    region: AbsorbingCrystalBox,
    cells: Vec<SubstrateCellPower>,
    power: CoherentCrystalWavePower,
    mode_cell_pairs: usize,
    workspace_bytes: usize,
}

impl PreparedCrystalBoxCells {
    #[allow(clippy::too_many_arguments)]
    fn try_new(
        thermal: &BoundLayeredThermalBody<'_>,
        source: &impl CrystalBoxSource,
        footprint_origin_m: [f64; 2],
        footprint_extent_m: [f64; 2],
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        owner_header_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let body = thermal.prepared_body();
        let count = body.cell_count();
        if count == 0
            || count > limits.maximum_cells
            || !body.conduction().geometry().is_planar()
            || !depth_start_m.is_finite()
            || depth_start_m < 0.0
            || footprint_origin_m.iter().any(|v| !v.is_finite())
            || footprint_extent_m.iter().any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid("coherent crystal heat requires a finite planar control volume"));
        }
        let area = body.conduction().area_m2();
        match_measure(
            Dd::from_f64(footprint_extent_m[0]) * Dd::from_f64(footprint_extent_m[1]),
            Dd::from_f64(area),
        )?;
        let mode_pairs = source.mode_pairs()?;
        let mode_cell_pairs = mode_pairs
            .checked_mul(count)
            .ok_or(invalid("coherent crystal cell work overflows"))?;
        let retained_bytes = count
            .checked_mul(core::mem::size_of::<SubstrateCellPower>())
            .and_then(|n| n.checked_add(owner_header_bytes))
            .ok_or(invalid("coherent crystal storage overflows"))?;
        let interval_bytes = count
            .checked_mul(core::mem::size_of::<f64>() + core::mem::size_of::<(CellId, f64)>())
            .and_then(|n| n.checked_add(core::mem::size_of::<CoherentCrystalWaveStep<'_>>()))
            .ok_or(invalid("coherent crystal interval storage overflows"))?;
        let available = limits
            .maximum_workspace_bytes
            .checked_sub(retained_bytes)
            .ok_or(invalid("coherent crystal storage exceeds its budget"))?;
        if mode_cell_pairs > limits.maximum_mode_cell_pairs || interval_bytes > available {
            return Err(invalid("coherent crystal heat exceeds its work or storage budget"));
        }
        admit_crystal_material(
            thermal,
            source.material(),
            &mut cancelled,
        )?;

        let query_limits = AbsorbingCrystalBoxLimits {
            maximum_mode_pairs: mode_pairs,
            maximum_workspace_bytes: available,
        };
        let ids = thermal.cell_ids();
        let mut cells = reserve(count)?;
        let mut layer_start = Dd::from_f64(depth_start_m);
        let mut left = depth_start_m;
        let mut absorbed = Dd::ZERO;
        let mut outward = [Dd::ZERO; 6];
        let mut query_bytes = 0;
        for (layer, &n) in body.conduction().layers().iter().zip(body.layer_cell_counts()) {
            let width = Dd::from_f64(layer.thickness_m());
            for local in 0..n {
                poll(&mut cancelled)?;
                let right = (layer_start
                    + width * Dd::from_f64((local + 1) as f64) / Dd::from_f64(n as f64))
                .to_f64();
                if !right.is_finite() || right <= left {
                    return Err(invalid("coherent crystal cell depth is unrepresentable"));
                }
                let region = AbsorbingCrystalBox::try_new(
                    [footprint_origin_m[0], footprint_origin_m[1], left],
                    [footprint_extent_m[0], footprint_extent_m[1], right - left],
                )?;
                match_measure(
                    Dd::from_f64(region.volume_m3()),
                    Dd::from_f64(area) * width / Dd::from_f64(n as f64),
                )?;
                let result = source.integrated_power(region, query_limits, &mut cancelled)?;
                query_bytes = query_bytes.max(result.workspace_bytes);
                let watts = result.absorbed_w;
                let faces = result.outward_face_w;
                if !watts.is_finite()
                    || watts < 0.0
                    || faces.iter().any(|value| !value.is_finite())
                {
                    return Err(invalid("coherent crystal cell power is unrepresentable"));
                }
                // Lateral faces partition along depth. Only the original
                // first/last z faces bound the assembled control volume;
                // internal cell faces are not counted as external power.
                for face in 0..4 {
                    outward[face] = outward[face] + Dd::from_f64(faces[face]);
                }
                if cells.is_empty() {
                    outward[4] = Dd::from_f64(faces[4]);
                }
                outward[5] = Dd::from_f64(faces[5]);
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
        if !finite(absorbed) || outward.iter().any(|value| !finite(*value)) {
            return Err(invalid("coherent crystal power reduction is unrepresentable"));
        }
        let outward_face_w = outward.map(Dd::to_f64);
        let power = CoherentCrystalWavePower {
            outward_face_w,
            absorbed_w: absorbed.to_f64(),
            residual_w: residual(outward_face_w, absorbed)?,
        };
        let region = AbsorbingCrystalBox::try_new(
            [footprint_origin_m[0], footprint_origin_m[1], depth_start_m],
            [footprint_extent_m[0], footprint_extent_m[1], left - depth_start_m],
        )?;
        let workspace_bytes = retained_bytes
            .checked_add(query_bytes.max(interval_bytes))
            .ok_or(invalid("coherent crystal workspace overflows"))?;
        if workspace_bytes > limits.maximum_workspace_bytes {
            return Err(invalid("coherent crystal workspace exceeds its budget"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            region,
            cells,
            power,
            mode_cell_pairs,
            workspace_bytes,
        })
    }

    /// Add the cached cell heat once on the original accepted thermal clock.
    /// Six signed face energies describe the optical account; only dielectric
    /// heat enters the native solve. Failure/cancellation preserves initial.
    #[allow(clippy::too_many_arguments)]
    fn advance_to<'t>(
        &self,
        thermal: &BoundLayeredThermalBody<'_>,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCrystalWaveStep<'t>, LayeredWaveHeatingError> {
        let result = advance_cell_heat_to(
            thermal,
            &self.cells,
            initial,
            target_s,
            conduction,
            controls,
            workspace,
            |duration, absorbed| {
                let mut faces = [0.0; 6];
                for (joules, watts) in faces.iter_mut().zip(self.power.outward_face_w) {
                    *joules = watts * duration;
                    if !joules.is_finite() || (watts != 0.0 && duration > 0.0 && *joules == 0.0) {
                        return Err(invalid("coherent crystal signed face energy is unrepresentable"));
                    }
                }
                if !finite(absorbed) {
                    return Err(invalid("coherent crystal absorbed energy is unrepresentable"));
                }
                Ok(CoherentCrystalWaveEnergy {
                    outward_face_j: faces,
                    absorbed_j: absorbed.to_f64(),
                    residual_j: residual(faces, absorbed)?,
                })
            },
            cancelled,
        )?;
        Ok(CoherentCrystalWaveStep {
            thermal: result.thermal,
            heat_report: result.heat_report,
            optical: result.account,
            cell_energy_j: result.cell_energy_j,
        })
    }
}

fn match_measure(actual: Dd, expected: Dd) -> Result<(), LayeredWaveHeatingError> {
    if !finite(actual)
        || !finite(expected)
        || actual.to_f64() <= 0.0
        || expected.to_f64() <= 0.0
        || ((actual - expected) / expected).to_f64().abs() > 64.0 * f64::EPSILON
    {
        return Err(invalid(
            "coherent optical footprint and cell volume must match the original thermal geometry",
        ));
    }
    Ok(())
}

fn residual(faces: [f64; 6], absorbed: Dd) -> Result<f64, LayeredWaveHeatingError> {
    let value = faces
        .into_iter()
        .fold(-absorbed, |sum, face| sum - Dd::from_f64(face));
    if finite(value) {
        Ok(value.to_f64())
    } else {
        Err(invalid("coherent crystal finite-volume account is unrepresentable"))
    }
}

fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
