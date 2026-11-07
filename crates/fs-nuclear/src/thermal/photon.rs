//! Planar photon transport derived from the original material cells.
//!
//! The declared reduction treats ordered, contiguous, collinear cell segments
//! as homogeneous layers of a common cross-sectional area. Original volume and
//! mass determine thickness and reference density; no thermal emissivity,
//! equilibrium density or replacement material law supplies photon physics.
//! Photon mass kernels are constant over their declared material applicability.

use std::collections::BTreeSet;

use super::*;
use crate::endf::decay::spectra::DecayRadiation;
use crate::radiation::transport::material::PhotonMassMaterial;
use crate::radiation::transport::refinement::{
    PhotonSlabRefinement, RefinedPhotonSlabResponse, refine_slab_response,
};
use crate::radiation::transport::{
    PhotonOrdinate, PhotonSlabCell, PhotonSlabLimits, PhotonSlabResponse, admit_slab_dimensions,
    prepare_slab_response,
};

/// One original cell in geometric order along the planar stack's normal.
#[derive(Clone, Debug)]
pub struct PhotonSlabLayer {
    pub cell_id: CellId,
    /// Fixed mass coefficients for this cell's original material identity.
    pub material: Arc<PhotonMassMaterial>,
}

/// Explicit geometric reduction to laterally infinite planar layers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonSlabGeometry {
    /// Common physical cross-sectional area [m²]. Each original cell must have
    /// volume / area equal to its centreline depth within the position bound.
    pub area_m2: f64,
    /// Absolute coordinate tolerance [m] for endpoints, common direction and
    /// layer contiguity. This is a geometry admission, not a transport error.
    pub position_tolerance_m: f64,
}

struct PhotonSlabInput {
    radiation: DecayRadiation,
    edges_j: Vec<f64>,
    cells: Vec<PhotonSlabCell>,
    source_ids: Vec<u64>,
    receivers: Vec<RadiationReceiver>,
}

impl<S: DecayRadiationHeatSource> NuclearHeating<S> {
    /// Solve steady photon transport through the original material stack and
    /// attach its response to the existing radioactive heating transaction.
    ///
    /// Every source must occupy one supplied layer. Source emission is uniform
    /// and isotropic in that layer. Mass/volume set reference density, and
    /// volume/area sets depth; located centreline endpoints must describe that
    /// same ordered planar stack. All layers, including transparent or purely
    /// scattering layers, retain the existing material/geometry binding.
    ///
    /// This is a laterally infinite, vacuum-ended, fixed-property slab model.
    /// It does not infer three-dimensional shielding, thermal expansion,
    /// changing opacity, secondary-electron transport or photon flight time.
    /// The supplied kernels must apply over the subsequent material states.
    /// Source, transport and receiver admission precede any accepted-history
    /// mutation. The returned diagnostics belong to the prepared transport solve.
    #[allow(clippy::too_many_arguments)]
    pub fn with_photon_slab(
        self,
        thermal: &DepositionSimulation<'_>,
        data_id: String,
        geometry: PhotonSlabGeometry,
        layers: Vec<PhotonSlabLayer>,
        emitters: Vec<DecayRadiationEmitter>,
        ordinates: Vec<PhotonOrdinate>,
        limits: PhotonSlabLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, PhotonSlabResponse), FissionThermalError> {
        let input = self.photon_slab_input(
            thermal,
            geometry,
            layers,
            &emitters,
            ordinates.len(),
            limits,
            &mut checkpoint,
        )?;
        let prepared = photon_call(&mut checkpoint, |radiation_checkpoint| {
            prepare_slab_response(
                data_id,
                input.radiation,
                input.edges_j,
                input.cells,
                input.source_ids,
                ordinates,
                limits,
                radiation_checkpoint,
            )
        })?;
        let coupled = self.with_grouped_decay_radiation(
            thermal,
            Arc::clone(&prepared.response),
            emitters,
            input.receivers,
            checkpoint,
        )?;
        Ok((coupled, prepared))
    }

    /// Resolve the flat-source spatial approximation before binding radioactive
    /// heating to the same original material cells as [`Self::with_photon_slab`].
    ///
    /// Refinement divides only the photon solve: nuclear inventories, original
    /// material masses, thermal cells and source identities remain unchanged.
    /// Accepted sublayer heat is reduced back to each original receiver before
    /// the existing grouped radiation transaction can credit it. A refused or
    /// interrupted refinement leaves the original heating binding unchanged.
    ///
    /// The returned grid change is an observed spatial convergence diagnostic.
    /// Angular/energy resolution, fixed material applicability and the planar
    /// transport assumptions still require independent justification.
    #[allow(clippy::too_many_arguments)]
    pub fn with_refined_photon_slab(
        self,
        thermal: &DepositionSimulation<'_>,
        data_id: String,
        geometry: PhotonSlabGeometry,
        layers: Vec<PhotonSlabLayer>,
        emitters: Vec<DecayRadiationEmitter>,
        ordinates: Vec<PhotonOrdinate>,
        limits: PhotonSlabLimits,
        refinement: PhotonSlabRefinement,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, RefinedPhotonSlabResponse), FissionThermalError> {
        let input = self.photon_slab_input(
            thermal,
            geometry,
            layers,
            &emitters,
            ordinates.len(),
            limits,
            &mut checkpoint,
        )?;
        let prepared = photon_call(&mut checkpoint, |radiation_checkpoint| {
            refine_slab_response(
                data_id,
                input.radiation,
                input.edges_j,
                input.cells,
                input.source_ids,
                ordinates,
                limits,
                refinement,
                radiation_checkpoint,
            )
        })?;
        let coupled = self.with_grouped_decay_radiation(
            thermal,
            Arc::clone(&prepared.transport.response),
            emitters,
            input.receivers,
            checkpoint,
        )?;
        Ok((coupled, prepared))
    }

    fn photon_slab_input(
        &self,
        thermal: &DepositionSimulation<'_>,
        geometry: PhotonSlabGeometry,
        layers: Vec<PhotonSlabLayer>,
        emitters: &[DecayRadiationEmitter],
        ordinates: usize,
        limits: PhotonSlabLimits,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PhotonSlabInput, FissionThermalError> {
        checkpoint()?;
        self.validate(thermal, self.time_s, &mut || Ok(()))?;
        if self.pending.is_some() {
            return Err(FissionThermalError::Invalid(
                "cannot change photon transport during a pending interval",
            ));
        }
        let groups = layers
            .first()
            .map_or(0, |layer| layer.material.group_edges_j().len() - 1);
        // Apply the transport owner's actual caps before copying any grid or
        // expanding immutable mass coefficients into macroscopic cell rows.
        admit_slab_dimensions(layers.len(), groups, self.sources.len(), ordinates, limits)?;
        if !geometry.area_m2.is_finite()
            || geometry.area_m2 <= 0.0
            || !nonnegative(geometry.position_tolerance_m)
        {
            return Err(FissionThermalError::Invalid(
                "positive common photon area and finite nonnegative geometry tolerance required",
            ));
        }
        let source_ids: Vec<_> = self.sources.keys().map(|cell| cell.0).collect();
        let declared: BTreeSet<_> = emitters.iter().map(|emitter| emitter.cell_id.0).collect();
        if source_ids.is_empty()
            || emitters.len() != source_ids.len()
            || !declared.iter().copied().eq(source_ids.iter().copied())
        {
            return Err(FissionThermalError::Invalid(
                "photon emitters must cover the nonempty original source cohort exactly",
            ));
        }
        let radiation = layers[0].material.radiation();
        let edges = layers[0].material.group_edges_j().to_vec();
        let first = thermal
            .cells()
            .get(&layers[0].cell_id)
            .ok_or(FissionThermalError::InactiveCell(layers[0].cell_id))?;
        let direction = normal(first.geometry)?;
        let mut expected_start = first.geometry.start_m;
        let mut previous_end = expected_start;
        let mut seen = BTreeSet::new();
        let mut receivers = Vec::with_capacity(layers.len());
        let mut physical = Vec::with_capacity(layers.len());
        for layer in layers {
            checkpoint()?;
            if !seen.insert(layer.cell_id) {
                return Err(FissionThermalError::Invalid(
                    "duplicate photon material cell",
                ));
            }
            if layer.material.radiation() != radiation || layer.material.group_edges_j() != edges {
                return Err(FissionThermalError::Invalid(
                    "photon materials must share the same radiation kind and energy grid",
                ));
            }
            let cell = thermal
                .cells()
                .get(&layer.cell_id)
                .ok_or(FissionThermalError::InactiveCell(layer.cell_id))?;
            let material_id = MaterialId(layer.material.material_id());
            if cell.material_id != material_id {
                return Err(FissionThermalError::MaterialMismatch {
                    cell_id: layer.cell_id,
                    expected: material_id,
                    actual: cell.material_id,
                });
            }
            let thickness = cell.geometry.volume_m3 / geometry.area_m2;
            let density = cell.mass_kg / cell.geometry.volume_m3;
            if !thickness.is_finite() || thickness <= 0.0 || !density.is_finite() || density <= 0.0
            {
                return Err(FissionThermalError::Arithmetic);
            }
            // A tolerance at least half the layer depth cannot establish even
            // its ordered extent. Compare by division to avoid half underflow.
            if geometry.position_tolerance_m / thickness >= 0.5 {
                return Err(FissionThermalError::Invalid(
                    "photon geometry tolerance must resolve every material layer",
                ));
            }
            let alignment: f64 = normal(cell.geometry)?
                .into_iter()
                .zip(direction)
                .map(|(actual, reference)| actual * reference)
                .sum();
            if !alignment.is_finite() || alignment <= 0.0 {
                return Err(FissionThermalError::Invalid(
                    "original photon layers must point forward along the stack normal",
                ));
            }
            let expected_end = std::array::from_fn(|axis| {
                direction[axis].mul_add(thickness, expected_start[axis])
            });
            if !matches_position(
                cell.geometry.start_m,
                expected_start,
                geometry.position_tolerance_m,
            ) || !matches_position(
                cell.geometry.start_m,
                previous_end,
                geometry.position_tolerance_m,
            ) || !matches_position(
                cell.geometry.end_m,
                expected_end,
                geometry.position_tolerance_m,
            ) || cell.geometry.start_m == cell.geometry.end_m
                || expected_start == expected_end
            {
                return Err(FissionThermalError::Invalid(
                    "original cells must form a contiguous collinear photon stack with volume/area depths",
                ));
            }
            thermal.temperature_k(layer.cell_id)?;
            receivers.push(RadiationReceiver {
                cell_id: layer.cell_id,
                material_id,
            });
            physical.push((layer, thickness, density));
            expected_start = expected_end;
            previous_end = cell.geometry.end_m;
        }
        if source_ids.iter().any(|&id| !seen.contains(&CellId(id))) {
            return Err(FissionThermalError::Invalid(
                "every photon source must occupy an original material layer",
            ));
        }
        let cells = photon_call(checkpoint, |radiation_checkpoint| {
            let mut cells = Vec::with_capacity(physical.len());
            for (layer, thickness, density) in physical {
                cells.push(layer.material.macroscopic_cell(
                    layer.cell_id.0,
                    thickness,
                    density,
                    &mut *radiation_checkpoint,
                )?);
            }
            Ok(cells)
        })?;
        Ok(PhotonSlabInput {
            radiation,
            edges_j: edges,
            cells,
            source_ids,
            receivers,
        })
    }
}

// Material conversion and both transport routes share the same exact callback
// translation; no wrapped cancellation can replace the original source error.
fn photon_call<R>(
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    run: impl FnOnce(&mut dyn FnMut() -> Result<(), RadiationError>) -> Result<R, RadiationError>,
) -> Result<R, FissionThermalError> {
    let mut checkpoint_error = None;
    let result = run(&mut || {
        checkpoint().map_err(|error| {
            checkpoint_error = Some(error);
            RadiationError::Cancelled
        })
    });
    if let Some(error) = checkpoint_error {
        Err(error.into())
    } else {
        result.map_err(Into::into)
    }
}

fn matches_position(actual: [f64; 3], expected: [f64; 3], tolerance: f64) -> bool {
    actual.into_iter().zip(expected).all(|(actual, expected)| {
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance
    })
}

fn normal(geometry: CellGeometry) -> Result<[f64; 3], FissionThermalError> {
    let delta: [f64; 3] = std::array::from_fn(|axis| geometry.end_m[axis] - geometry.start_m[axis]);
    let scale = delta.iter().copied().map(f64::abs).fold(0.0, f64::max);
    if !delta.iter().all(|value| value.is_finite()) || scale == 0.0 {
        return Err(FissionThermalError::Invalid(
            "positive finite original photon-layer centreline extent required",
        ));
    }
    let scaled = delta.map(|value| value / scale);
    let length = fs_math::det::sqrt(scaled.iter().map(|value| value * value).sum());
    Ok(scaled.map(|value| value / length))
}
