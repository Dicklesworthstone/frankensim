//! Physical terminal faces of the original layered thermal control volumes.
//!
//! Selected channels conduct from the cell's volume-centre trace through its
//! original material and an authored film. Partial patches are independent
//! normal streamtubes; this one-dimensional model does not infer lateral heat
//! spreading. The original thermal owner retains coverage, baths and accounts.

use super::{BoundLayeredThermalBody, LayerStackGeometry};
use crate::thermal::conduction::{
    LayerConductionConfig, ReservoirLayerGeometry, ReservoirLayerResponse,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, backward_euler::BoundaryConductances,
};
use std::collections::BTreeMap;

/// Physical lower/inner or upper/outer terminal face, in original layer order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum LayeredFace {
    Inner,
    Outer,
}

/// An existing prescribed bath; its original temperature and coverage survive.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum LayeredBoundaryChannel {
    Ambient,
    Reservoir,
}

/// Positive physical area assigned to one existing bath on one terminal face.
/// Each selected cell/channel must be covered completely by its declared
/// patches. At most one patch per (face, channel) is admitted, hence at most four.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredSurfacePatch {
    pub face: LayeredFace,
    pub channel: LayeredBoundaryChannel,
    pub area_m2: f64,
}

/// The actual material-face state and outward heat at one immutable candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredSurfaceResponse {
    pub patch: LayeredSurfacePatch,
    pub cell_id: CellId,
    pub cell_temperature_k: f64,
    pub bath_temperature_k: f64,
    /// Positive power leaves the cell; derivatives are with respect to
    /// [cell temperature, bath temperature]. Insulation has exactly zero heat.
    pub response: ReservoirLayerResponse,
}

struct BoundPatch {
    patch: LayeredSurfacePatch,
    cell_id: CellId,
    layer: usize,
    geometry: ReservoirLayerGeometry,
    bath_k: f64,
    film_resistance_k_w: Option<f64>,
}

/// Geometry and original boundary data are prepared once, not at every trial.
#[derive(Default)]
pub(super) struct SurfaceBoundaries {
    patches: Vec<BoundPatch>,
}

#[derive(Default)]
pub(super) struct ResolvedSurfaceBoundaries {
    pub boundaries: BoundaryConductances,
    pub derivatives: BoundaryConductances,
    pub responses: Vec<LayeredSurfaceResponse>,
}

impl SurfaceBoundaries {
    pub fn new(
        binding: &BoundLayeredThermalBody<'_>,
        initial: &DepositionSimulation<'_>,
        patches: &[LayeredSurfacePatch],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if patches.len() > 4 {
            return Err(bad("at most four layered surface patches are supported"));
        }
        let mut ordered = patches.to_vec();
        ordered.sort_by_key(|patch| (patch.face, patch.channel));
        if ordered
            .windows(2)
            .any(|pair| pair[0].face == pair[1].face && pair[0].channel == pair[1].channel)
        {
            return Err(bad("duplicate layered surface face/channel"));
        }
        let body = binding.body;
        let mut prepared = Vec::with_capacity(ordered.len());
        let mut coverage = BTreeMap::<(CellId, LayeredBoundaryChannel), f64>::new();
        let mut face_coverage = [0.0; 2];
        for patch in ordered {
            checkpoint()?;
            if !patch.area_m2.is_finite() || patch.area_m2 <= 0.0 {
                return Err(bad(
                    "layered surface patch area must be positive and finite",
                ));
            }
            let (face_index, cell_index) = match patch.face {
                LayeredFace::Inner => (0, 0),
                LayeredFace::Outer => (1, body.cells.len() - 1),
            };
            let cell = &body.cells[cell_index];
            let id = binding.cell_id(cell_index);
            let boundary = initial.cells()[&id].boundary;
            let (native_area, film_g, bath_k) = match patch.channel {
                LayeredBoundaryChannel::Ambient => (
                    initial.exposed_area_m2(id)?,
                    patch.area_m2 * boundary.htc_w_m2_k,
                    boundary.ambient_k,
                ),
                LayeredBoundaryChannel::Reservoir => {
                    let reservoir = boundary.reservoir.ok_or_else(|| {
                        bad("surface reservoir patch needs an original reservoir")
                    })?;
                    (
                        reservoir.area_m2,
                        reservoir.conductance_w_k * (patch.area_m2 / reservoir.area_m2),
                        reservoir.temperature_k,
                    )
                }
            };
            if native_area <= 0.0 || !native_area.is_finite() {
                return Err(bad("selected surface channel has no original area"));
            }
            let film_resistance_k_w = if film_g == 0.0 {
                None
            } else {
                let resistance = 1.0 / film_g;
                if !film_g.is_finite()
                    || film_g < 0.0
                    || !resistance.is_finite()
                    || resistance <= 0.0
                {
                    return Err(DepositionError::Arithmetic);
                }
                Some(resistance)
            };
            let surface = match patch.face {
                LayeredFace::Inner => cell.geometry.start_m[0],
                LayeredFace::Outer => cell.geometry.end_m[0],
            };
            let geometry = body.conduction.geometry();
            let area = geometry.face_area_m2(surface).map_err(bad)?;
            face_coverage[face_index] += patch.area_m2;
            if !face_coverage[face_index].is_finite()
                || face_coverage[face_index] - area > area_roundoff(area)
            {
                return Err(bad(
                    "layered surface patches exceed their physical face area",
                ));
            }
            let full_face = match geometry {
                LayerStackGeometry::Planar { area_m2 } => {
                    ReservoirLayerGeometry::planar(area_m2, (cell.centroid - surface).abs())?
                }
                LayerStackGeometry::Cylindrical { length_m, .. } => {
                    ReservoirLayerGeometry::cylindrical_face(surface, cell.centroid, length_m)?
                }
                LayerStackGeometry::Spherical { .. } => {
                    ReservoirLayerGeometry::spherical_face(surface, cell.centroid)?
                }
            };
            // Original radial surface sums may carry one outward-rounded ulp.
            // Never grow the material streamtube beyond its actual full face.
            let geometry = full_face.with_area_fraction((patch.area_m2 / area).min(1.0))?;
            *coverage.entry((id, patch.channel)).or_default() += patch.area_m2;
            prepared.push(BoundPatch {
                patch,
                cell_id: id,
                layer: cell.layer,
                geometry,
                bath_k,
                film_resistance_k_w,
            });
        }
        for ((id, channel), covered) in coverage {
            checkpoint()?;
            let cell = &initial.cells()[&id];
            let expected = match channel {
                LayeredBoundaryChannel::Ambient => initial.exposed_area_m2(id)?,
                LayeredBoundaryChannel::Reservoir => {
                    cell.boundary
                        .reservoir
                        .ok_or_else(|| bad("surface reservoir is absent"))?
                        .area_m2
                }
            };
            // Coverage subtracts original contacts and reservoir patches. Its
            // rounding scale is the cell's physical surface, not a tiny remainder.
            if !covered.is_finite()
                || (covered - expected).abs() > area_roundoff(cell.geometry.surface_area_m2)
            {
                return Err(bad(
                    "surface patches must cover their entire original cell/channel",
                ));
            }
        }
        checkpoint()?;
        Ok(Self { patches: prepared })
    }

    pub fn resolve(
        &self,
        binding: &BoundLayeredThermalBody<'_>,
        candidate: &DepositionSimulation<'_>,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedSurfaceBoundaries, DepositionError> {
        let mut resolved = ResolvedSurfaceBoundaries::default();
        for bound in &self.patches {
            checkpoint()?;
            let cell_k = candidate.temperature_k(bound.cell_id)?;
            let response = if let Some(film) = bound.film_resistance_k_w {
                bound.geometry.evaluate_material_to_bath(
                    &binding.body.conduction.curves[bound.layer],
                    [cell_k, bound.bath_k],
                    film,
                    config,
                    &mut checkpoint,
                )?
            } else {
                ReservoirLayerResponse {
                    power_to_reservoir_w: 0.0,
                    conductance_w_k: 0.0,
                    endpoint_derivative_w_k: [0.0, 0.0],
                    material_face_temperature_k: cell_k,
                    heat_residual_w: 0.0,
                    iterations: 0,
                }
            };
            let (conductances, derivatives) = match bound.patch.channel {
                LayeredBoundaryChannel::Ambient => (
                    &mut resolved.boundaries.ambient_w_k,
                    &mut resolved.derivatives.ambient_w_k,
                ),
                LayeredBoundaryChannel::Reservoir => (
                    &mut resolved.boundaries.reservoir_w_k,
                    &mut resolved.derivatives.reservoir_w_k,
                ),
            };
            let conductance = conductances.entry(bound.cell_id).or_default();
            *conductance += response.conductance_w_k;
            let derivative = derivatives.entry(bound.cell_id).or_default();
            *derivative += response.endpoint_derivative_w_k[0];
            if !conductance.is_finite()
                || *conductance < 0.0
                || !derivative.is_finite()
                || *derivative < 0.0
            {
                return Err(DepositionError::Arithmetic);
            }
            resolved.responses.push(LayeredSurfaceResponse {
                patch: bound.patch,
                cell_id: bound.cell_id,
                cell_temperature_k: cell_k,
                bath_temperature_k: bound.bath_k,
                response,
            });
        }
        checkpoint()?;
        Ok(resolved)
    }
}

fn area_roundoff(scale: f64) -> f64 {
    8.0 * f64::EPSILON * scale
}
fn bad(reason: &'static str) -> DepositionError {
    DepositionError::InvalidInput(reason)
}
