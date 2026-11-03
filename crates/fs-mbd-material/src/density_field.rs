//! Source-bound P1 density, exact mass moments and conserved material transport.
//!
//! Density is resolved once at each reference vertex and interpolated linearly
//! inside each tetrahedron. The native inertia owner integrates that field;
//! neither a centroid sample nor a uniform-density inertia replaces it.
//! An affine cell deformation transports the original density by its Jacobian.
//! It does not query new material density or add/remove matter.

use std::collections::BTreeSet;

use fs_ga::Vec3;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_mbd::articulated::SpatialInertia;
use fs_qty::Density;

use crate::principal::{PrincipalFrame, PrincipalFrameError};
use crate::{
    MAX_BODY_REGIONS, MaterialInertiaError, VolumeMoments, merge_inertias, poll, positive,
};

/// Density-weighted nodal motion projected into the existing rigid-body owner.
#[path = "density_field/motion.rs"]
pub mod motion;

mod partition;
pub use partition::{TetrahedralDensityPart, TetrahedralDensitySplit};

/// Maximum reference vertices; the original body-region budget caps cells.
pub const MAX_DENSITY_VERTICES: usize = 4 * MAX_BODY_REGIONS;

/// One vertex's explicit position and complete source-backed density query.
#[derive(Clone, Debug)]
pub struct DensityVertexBinding<'a> {
    /// Position in the reference body's shared coordinate frame [m].
    pub position_m: Vec3,
    /// Immutable material card, without material-name density inference.
    pub material: &'a MaterialCard,
    /// Full physical state at this reference vertex.
    pub point: &'a QueryPoint,
    /// Original claim selection for this vertex's density.
    pub selection: MaterialPropertySelection,
}

/// Original density source and reference position of a shared mesh vertex.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDensityVertex {
    position_m: Vec3,
    material: ResolvedMaterialStatePoint,
}
impl ResolvedDensityVertex {
    /// Reference position [m].
    pub const fn position_m(&self) -> Vec3 {
        self.position_m
    }
    /// Complete source state and selected density receipt, retained unchanged.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Positive reference density at this vertex [kg/m3].
    pub fn density_kg_m3(&self) -> f64 {
        self.material
            .property(DENSITY_PROPERTY)
            .expect("density is the required source property")
            .value_si()
    }
}

/// Geometry and exact P1 mass distribution of one reference or deformed cell.
#[derive(Clone, Debug, PartialEq)]
pub struct DensityCell {
    geometry: VolumeMoments,
    vertex_density_kg_m3: [f64; 4],
    inertia: SpatialInertia,
}
impl DensityCell {
    /// Geometric volume and unweighted moments of this current tetrahedron.
    pub const fn geometry(&self) -> VolumeMoments {
        self.geometry
    }
    /// Density at each vertex, in the cell's published topology order [kg/m3].
    /// After deformation these are transported values, not new measurements.
    pub const fn vertex_density_kg_m3(&self) -> [f64; 4] {
        self.vertex_density_kg_m3
    }
    /// Native mass, density-weighted centroid and full centroidal inertia.
    pub const fn inertia(&self) -> SpatialInertia {
        self.inertia
    }
    /// Cell-average density, consistent with the original retained mass [kg/m3].
    pub fn mean_density_kg_m3(&self) -> f64 {
        self.inertia.mass() / self.geometry.volume_m3()
    }
}

/// Immutable reference body with continuous P1 density on shared mesh vertices.
///
/// Positive nodal densities make the field positive throughout each cell.
/// The interpolation and tetrahedral partition are explicit model assumptions,
/// not new material measurements or a mixing rule. All cells must be disjoint
/// except at shared faces/edges; coverage and overlap remain caller obligations.
#[derive(Clone, Debug, PartialEq)]
pub struct TetrahedralDensityBody {
    vertices: Vec<ResolvedDensityVertex>,
    tetrahedra: Vec<[usize; 4]>,
    determinants: Vec<f64>,
    cells: Vec<DensityCell>,
    inertia: SpatialInertia,
}

/// Current coordinates and the same reference mass distribution after a P1 map.
///
/// Each cell carries rho_current = rho_reference / J. Different cell Jacobians
/// can make current densities discontinuous at a shared vertex; the per-cell
/// values retain that physical distinction instead of averaging matter away.
#[derive(Clone, Debug, PartialEq)]
pub struct DeformedDensityBody<'a> {
    reference: &'a TetrahedralDensityBody,
    vertices_m: Vec<Vec3>,
    cells: Vec<DensityCell>,
    inertia: SpatialInertia,
}

impl TetrahedralDensityBody {
    /// Resolve the same density requirement at each vertex, then integrate the
    /// complete P1 field through the original native spatial-inertia owner.
    ///
    /// Cell vertex indices and cells are canonically sorted so input cell order
    /// does not change the mass reduction. Source queries occur once per vertex.
    /// Construction uses O(vertices + cells log cells) work and bounded storage;
    /// cancellation is checked before work, per source/cell and at publication.
    pub fn resolve(
        vertices: &[DensityVertexBinding<'_>],
        tetrahedra: &[[usize; 4]],
        max_vertices: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialInertiaError> {
        poll(&mut cancelled)?;
        if max_vertices == 0
            || max_vertices > MAX_DENSITY_VERTICES
            || vertices.len() < 4
            || vertices.len() > max_vertices
            || tetrahedra.is_empty()
            || tetrahedra.len() > MAX_BODY_REGIONS
        {
            return Err(invalid("bounded vertices and 1..=4096 tetrahedra required"));
        }
        let mut ordered = BTreeSet::new();
        for &tetrahedron in tetrahedra {
            poll(&mut cancelled)?;
            let mut cell = tetrahedron;
            cell.sort_unstable();
            if cell[3] >= vertices.len()
                || cell.windows(2).any(|pair| pair[0] == pair[1])
                || !ordered.insert(cell)
            {
                return Err(invalid(
                    "tetrahedra need four distinct valid vertices and no duplicates",
                ));
            }
        }
        let requirement = ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let mut resolved = Vec::with_capacity(vertices.len());
        for vertex in vertices {
            poll(&mut cancelled)?;
            finite_vertex(vertex.position_m)?;
            let material = resolve_material_state_point(
                vertex.material,
                vertex.point,
                core::slice::from_ref(&requirement),
                vertex.selection.clone(),
            )?;
            resolved.push(ResolvedDensityVertex {
                position_m: vertex.position_m,
                material,
            });
        }
        let tetrahedra: Vec<_> = ordered.into_iter().collect();
        let mut cells = Vec::with_capacity(tetrahedra.len());
        let mut determinants = Vec::with_capacity(tetrahedra.len());
        let mut combined = None;
        for &tet in &tetrahedra {
            poll(&mut cancelled)?;
            let positions = tet.map(|v| resolved[v].position_m);
            let geometry = VolumeMoments::tetrahedron(positions)?;
            let densities = tet.map(|v| resolved[v].density_kg_m3());
            let mean_density = mean_density(densities);
            let mass = geometry.volume_m3() * mean_density;
            let inertia = SpatialInertia::from_weighted_tetrahedron(mass, positions, densities)?;
            combined = Some(combine(combined, inertia)?);
            determinants.push(determinant(positions));
            cells.push(DensityCell {
                geometry,
                vertex_density_kg_m3: densities,
                inertia,
            });
        }
        let inertia = combined.expect("the tetrahedral partition is nonempty");
        inertia.matrix()?;
        poll(&mut cancelled)?;
        Ok(Self {
            vertices: resolved,
            tetrahedra,
            determinants,
            cells,
            inertia,
        })
    }

    /// Every original position, material state and density source receipt.
    pub fn vertices(&self) -> &[ResolvedDensityVertex] {
        &self.vertices
    }
    /// Canonical shared-vertex topology, also defining each cell's density order.
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        &self.tetrahedra
    }
    /// Exact per-cell mass, centroid and inertia, in canonical topology order.
    pub fn cells(&self) -> &[DensityCell] {
        &self.cells
    }
    /// Existing spatial inertia consumed by the actual articulated dynamics.
    pub const fn inertia(&self) -> SpatialInertia {
        self.inertia
    }
    /// Existing free-body, contact and hydrostatic frame; no second motion model.
    pub fn principal_frame(
        &self,
        relative_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PrincipalFrame, PrincipalFrameError> {
        PrincipalFrame::from_inertia(self.inertia, relative_tolerance, cancelled)
    }

    /// Transport reference matter through the actual shared current vertices.
    ///
    /// Every cell's total mass is copied exactly. Its relative barycentric
    /// density weights remain unchanged, while geometry, centroid and complete
    /// inertia follow the affine cell map. Repeated calls start from reference
    /// matter; no density refresh, cumulative mass drift or isotropic substitute
    /// is introduced. This publishes a configuration, not live-state momentum
    /// transfer, elastic stresses or deformation velocities.
    pub fn deform(
        &self,
        vertices_m: &[Vec3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DeformedDensityBody<'_>, MaterialInertiaError> {
        poll(&mut cancelled)?;
        if vertices_m.len() != self.vertices.len() {
            return Err(invalid("deformation must retain every reference vertex"));
        }
        for &vertex in vertices_m {
            poll(&mut cancelled)?;
            finite_vertex(vertex)?;
        }
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut combined = None;
        for (index, &tet) in self.tetrahedra.iter().enumerate() {
            poll(&mut cancelled)?;
            let original = &self.cells[index];
            let cell = if tet
                .iter()
                .all(|&v| vertices_m[v] == self.vertices[v].position_m)
            {
                original.clone()
            } else {
                let positions = tet.map(|v| vertices_m[v]);
                let current_determinant = determinant(positions);
                if !current_determinant.is_finite()
                    || current_determinant == 0.0
                    || current_determinant.is_sign_positive()
                        != self.determinants[index].is_sign_positive()
                {
                    return Err(invalid(
                        "deformation inverts, collapses or overflows a cell",
                    ));
                }
                let geometry = VolumeMoments::tetrahedron(positions)?;
                let current_mean = original.inertia.mass() / geometry.volume_m3();
                let reference_mean = mean_density(original.vertex_density_kg_m3);
                let densities = original
                    .vertex_density_kg_m3
                    .map(|rho| current_mean * (rho / reference_mean));
                if !positive(current_mean) || !densities.into_iter().all(positive) {
                    return Err(invalid("transported nodal density is unrepresentable"));
                }
                let inertia = SpatialInertia::from_weighted_tetrahedron(
                    original.inertia.mass(),
                    positions,
                    original.vertex_density_kg_m3,
                )?;
                DensityCell {
                    geometry,
                    vertex_density_kg_m3: densities,
                    inertia,
                }
            };
            combined = Some(combine(combined, cell.inertia)?);
            cells.push(cell);
        }
        let inertia = combined.expect("the reference partition is nonempty");
        inertia.matrix()?;
        poll(&mut cancelled)?;
        Ok(DeformedDensityBody {
            reference: self,
            vertices_m: vertices_m.to_vec(),
            cells,
            inertia,
        })
    }
}

impl DeformedDensityBody<'_> {
    /// Original material cards, selected density receipts and interpolation model.
    pub const fn reference(&self) -> &TetrahedralDensityBody {
        self.reference
    }
    /// Actual current shared vertices [m].
    pub fn vertices_m(&self) -> &[Vec3] {
        &self.vertices_m
    }
    /// Unchanged canonical topology and cell order.
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        self.reference.tetrahedra()
    }
    /// Current cell densities and full inertias with exact reference masses.
    pub fn cells(&self) -> &[DensityCell] {
        &self.cells
    }
    /// Existing native spatial inertia of the deformed material distribution.
    pub const fn inertia(&self) -> SpatialInertia {
        self.inertia
    }
    /// Same principal-frame adapter as a resolved uniform-region body.
    pub fn principal_frame(
        &self,
        relative_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PrincipalFrame, PrincipalFrameError> {
        PrincipalFrame::from_inertia(self.inertia, relative_tolerance, cancelled)
    }
}

fn mean_density(values: [f64; 4]) -> f64 {
    let scale = values.into_iter().fold(0.0_f64, f64::max);
    scale
        * values
            .into_iter()
            .map(|rho| 0.25 * (rho / scale))
            .sum::<f64>()
}
fn determinant(vertices: [Vec3; 4]) -> f64 {
    (vertices[1] - vertices[0]).dot((vertices[2] - vertices[0]).cross(vertices[3] - vertices[0]))
}
fn combine(
    previous: Option<SpatialInertia>,
    cell: SpatialInertia,
) -> Result<SpatialInertia, MaterialInertiaError> {
    match previous {
        Some(previous) => merge_inertias(previous, cell),
        None => Ok(cell),
    }
}
fn finite_vertex(vertex: Vec3) -> Result<(), MaterialInertiaError> {
    if [vertex.x, vertex.y, vertex.z]
        .into_iter()
        .all(f64::is_finite)
    {
        Ok(())
    } else {
        Err(invalid("density mesh vertices must be finite"))
    }
}
fn invalid(message: &'static str) -> MaterialInertiaError {
    MaterialInertiaError::Invalid(message)
}
