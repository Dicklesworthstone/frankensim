//! Storage of thermally thin, perfectly bonded boundary coatings.
//!
//! A coating shares the substrate's original P1 boundary temperature. Its
//! areal heat capacity adds storage without introducing another temperature,
//! contact resistance, tangential sheet-conduction operator, or thin volume
//! mesh. Do not add a layer whose volume storage is already represented.

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_sparse::Coo;

use super::{BackwardEuler, poll};
use crate::ConductionError;
use crate::material::TemperatureSpan;
use crate::transient::MaterialHeatCapacity;

/// Frozen, nonnegative areal heat capacity, J/(m²·K).
///
/// Zero means no added surface storage; the underlying volume still retains
/// its positive heat capacity. Temperature-dependent calorics and latent heat
/// require a corresponding enthalpy model, not a frozen value in this lane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArealHeatCapacity {
    value_j_per_m2_k: f64,
}

impl ArealHeatCapacity {
    /// Declare a finite, nonnegative areal heat capacity in J/(m²·K).
    ///
    /// # Errors
    /// Non-finite or negative values are refused.
    pub fn declared(value_j_per_m2_k: f64) -> Result<Self, ConductionError> {
        crate::require_finite("areal heat capacity", value_j_per_m2_k)?;
        if value_j_per_m2_k < 0.0 {
            return Err(invalid("areal heat capacity must be nonnegative"));
        }
        Ok(Self { value_j_per_m2_k })
    }

    /// Areal heat capacity, J/(m²·K).
    #[must_use]
    pub const fn value_j_per_m2_k(self) -> f64 {
        self.value_j_per_m2_k
    }
}

/// A frozen P1 areal-capacity field on one original boundary triangle.
///
/// `boundary_face` is a slot in `ConductionMesh::boundary()`, and the three
/// coefficients follow that face's original `vertices` order. Adjacent faces
/// may represent different coatings. Faces must appear at most once in one
/// addition; combine coincident layers' areal capacities explicitly.
///
/// The layer is thermally thin and perfectly bonded to the solid: its
/// temperature is the solid boundary trace, including at prescribed nodes.
/// This adds no contact resistance or in-plane film conduction. A separately
/// meshed coating must not be counted again here.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceHeatCapacity<'s> {
    boundary_face: usize,
    nodal: [ArealHeatCapacity; 3],
    material: Option<&'s MaterialHeatCapacity>,
}

impl<'s> SurfaceHeatCapacity<'s> {
    /// A uniform declared coating on one original boundary face.
    #[must_use]
    pub const fn uniform(boundary_face: usize, capacity: ArealHeatCapacity) -> Self {
        Self {
            boundary_face,
            nodal: [capacity; 3],
            material: None,
        }
    }

    /// A declared P1 field in the original boundary face's vertex order.
    #[must_use]
    pub const fn nodal(boundary_face: usize, capacity: [ArealHeatCapacity; 3]) -> Self {
        Self {
            boundary_face,
            nodal: capacity,
            material: None,
        }
    }

    /// Bind one admitted density/Cp source to three coating thicknesses, m.
    ///
    /// `C_s = thickness * rho * c_p`. The source owner already requires
    /// validity-wide constant scalar calorics; a temperature-dependent material
    /// is not sampled and silently frozen. Preparation retains the original
    /// source bundle and intersects its temperature domain at every vertex
    /// receiving storage from the face.
    ///
    /// Thickness follows the original face's vertex order, is nonnegative,
    /// and is frozen in the reference geometry. This does not add expansion,
    /// ablation or a new material state. Use `[thickness_m; 3]` for a uniform film.
    ///
    /// # Errors
    /// Non-finite or negative thickness, or unrepresentable areal capacity.
    pub fn from_material(
        boundary_face: usize,
        material: &'s MaterialHeatCapacity,
        thickness_m: [f64; 3],
    ) -> Result<Self, ConductionError> {
        let mut nodal = [ArealHeatCapacity {
            value_j_per_m2_k: 0.0,
        }; 3];
        for (capacity, thickness) in nodal.iter_mut().zip(thickness_m) {
            crate::require_finite("surface layer thickness", thickness)?;
            if thickness < 0.0 {
                return Err(invalid("surface layer thickness must be nonnegative"));
            }
            let value = thickness * material.value_j_per_m3_k();
            *capacity = ArealHeatCapacity::declared(value)?;
            if thickness > 0.0 && value == 0.0 {
                return Err(invalid(
                    "positive coating thickness underflows its areal capacity",
                ));
            }
        }
        Ok(Self {
            boundary_face,
            nodal,
            material: Some(material),
        })
    }

    /// Slot in the prepared mesh's original boundary-face table.
    #[must_use]
    pub const fn boundary_face(self) -> usize {
        self.boundary_face
    }

    /// P1 coefficients in the original boundary face's vertex order.
    #[must_use]
    pub const fn nodal_capacity(self) -> [ArealHeatCapacity; 3] {
        self.nodal
    }

    /// The optional original material source; declared fields have no source.
    #[must_use]
    pub const fn material(self) -> Option<&'s MaterialHeatCapacity> {
        self.material
    }

    // One exact row-moment owner for temperature and primary-enthalpy storage.
    // Keep the original P1 thickness/capacity correlation and arithmetic order.
    pub(crate) fn nodal_moments(
        self,
        mesh: &crate::ConductionMesh,
    ) -> Result<([usize; 3], [f64; 3]), ConductionError> {
        let face = mesh.boundary().get(self.boundary_face).ok_or_else(|| {
            invalid("surface capacity face is outside the original boundary table")
        })?;
        let vertices = face.vertices.map(|vertex| vertex as usize);
        let values = self.nodal.map(ArealHeatCapacity::value_j_per_m2_k);
        if values.iter().all(|&value| value == 0.0) {
            return Ok((vertices, [0.0; 3]));
        }
        let mut moments = [0.0; 3];
        for (local, entry) in moments.iter_mut().enumerate() {
            *entry = face.area
                * (values[local] / 6.0
                    + values[(local + 1) % 3] / 12.0
                    + values[(local + 2) % 3] / 12.0);
            crate::require_finite("surface heat-capacity moment", *entry)?;
            if *entry <= 0.0 {
                return Err(invalid(
                    "positive surface storage has an unrepresentable face moment",
                ));
            }
        }
        Ok((vertices, moments))
    }
}

impl<'m> BackwardEuler<'m> {
    /// Add perfectly bonded surface storage to the prepared volume capacity.
    ///
    /// For a triangle of area `A`, the lumped row moment is exactly
    /// `integral(C_s N_i dA) = A * (2 C_i + C_j + C_k) / 12`.
    /// Thus a P1 thickness/capacity field preserves its spatial correlation
    /// with a P1 temperature; replacing it with a face average does not.
    ///
    /// All existing backward-Euler endpoint, nonlinear-conductivity,
    /// contact/radiation, reaction and adjoint paths use this combined storage.
    /// The coating adds neither a second temperature nor a heat-transfer
    /// resistance. Its capacity remains frozen. Omit storage already present
    /// in a separately meshed layer.
    ///
    /// The operation is additive. A face may occur once per call; repeated
    /// calls intentionally add further layers. Empty or all-zero additions
    /// preserve the original operator exactly and retain no inactive sources.
    ///
    /// # Errors
    /// Invalid or duplicate original face slots, non-finite/unrepresentable
    /// integrated storage, disjoint active source domains, or cancellation.
    /// No partly prepared model is returned.
    pub fn with_surface_heat_capacity(
        mut self,
        cx: &Cx<'_>,
        surfaces: &[SurfaceHeatCapacity<'_>],
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let n = self.mesh.vertex_count();
        let mut added = vec![0.0; n];
        let mut seen = BTreeSet::new();
        let mut active = false;
        for (slot, surface) in surfaces.iter().enumerate() {
            if slot % 512 == 0 {
                poll(cx, slot)?;
            }
            let (vertices, moments) = surface.nodal_moments(self.mesh)?;
            if !seen.insert(surface.boundary_face) {
                return Err(invalid("surface capacity face is declared more than once"));
            }
            if moments.iter().all(|&value| value == 0.0) {
                continue;
            }
            active = true;
            for (&vertex, entry) in vertices.iter().zip(moments) {
                added[vertex] += entry;
                crate::require_finite("assembled surface heat capacity", added[vertex])?;
            }
            if let Some(material) = surface.material {
                let spans = self
                    .capacity_temperature_spans
                    .get_or_insert_with(|| vec![TemperatureSpan::Unbounded; n]);
                // A zero coefficient at one corner does not remove its
                // positive basis moment when the rest of the face is active.
                for vertex in vertices {
                    let span = spans[vertex].intersect(material.temperature_span());
                    if matches!(span, TemperatureSpan::Sampled { low, high } if low > high) {
                        return Err(invalid(
                            "volume and surface capacity sources have disjoint temperature support",
                        ));
                    }
                    spans[vertex] = span;
                }
                if !self.material_capacities.contains(material) {
                    self.material_capacities.push(material.clone());
                }
            }
        }
        if !active {
            poll(cx, surfaces.len())?;
            return Ok(self);
        }
        let mut combined = Coo::new(n, n);
        for (row, &surface) in added.iter().enumerate() {
            if row % 512 == 0 {
                poll(cx, row)?;
            }
            let (columns, values) = self.capacity.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if column == row {
                    crate::require_finite(
                        "combined volume and surface heat capacity",
                        value + surface,
                    )?;
                }
                combined.push(row, column, value);
            }
            if surface > 0.0 {
                combined.push(row, row, surface);
            }
        }
        let combined = combined.assemble();
        if self.volume_capacity.is_none() {
            // Move the original bulk operator into its retained owner before
            // replacing physical storage. Subtracting the coating afterwards
            // could lose a small bulk contribution to floating-point rounding.
            self.volume_capacity = Some(std::mem::replace(&mut self.capacity, combined));
        } else {
            // Further bonded layers change physical storage, never the original
            // volume-capacity coordinate differentiated by the product.
            self.capacity = combined;
        }
        poll(cx, n)?;
        Ok(self)
    }
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "surface heat capacity",
        what: what.to_string(),
    }
}

