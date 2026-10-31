//! Automatic rigid-sphere contact with the CURRENT exterior of a P1 tet body.
//!
//! The original pressure/load geometry owner admits the closed reference
//! topology. Preparation binds each exterior triangle to its original tet and
//! caches the query owner's BVH. Queries refit current coordinates, retain
//! outward winding and map triangle barycentrics into the original consistent
//! mass/velocity trace. Shared edges/vertices produce one canonical witness.
//!
//! A sphere's actual material arm and the solid's actual material point remain
//! distinct at positive gap. Geometry never moves either point to hide the
//! angular couple from a tangential impulse. Interior sphere centres, collapsed
//! or inverted elements, unresolvable geometry and exhausted bounds refuse.
//! Raw separation queries may report a negative sphere gap for event bracketing;
//! contact admission refuses it. This is discrete-configuration geometry, not
//! a claim of continuous collision detection or global mesh injectivity. The
//! original volume mesh must be conforming and globally non-overlapping.

use std::collections::BTreeMap;

use super::super::impulse::FiniteSolidPoint;
use super::super::{FiniteSolidBody, FiniteSolidError};
use super::FiniteSolidRigidContact;
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::Cx;
use fs_mbd::{DynamicsError, RigidBodyState, Vec3};
use fs_query::{PreparedTriangleSurface, TriangleFeature, TrianglePoint, TriangleQueryError};
use fs_solid::hyper3::{HyperTetError, pressure::ClosedTetPressureSurface};

/// One geometric sphere carried by an existing rigid body's material frame.
/// Multiple offset spheres may form a declared compound rigid geometry; their
/// mass/inertia remain those supplied to the original rigid mechanics owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidSphere {
    pub id: u64,
    pub rigid_body: usize,
    pub center_body_m: Vec3,
    pub radius_m: f64,
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Explicit geometric activation, quality, work and output bounds.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidSurfaceControls {
    pub maximum_spheres: usize,
    /// Aggregate BVH visits, triangle tests and tetrahedron containment tests.
    /// Linear refit work is additionally bounded by the prepared mesh size.
    pub maximum_query_visits: usize,
    /// Aggregate nearby triangles before shared-feature deduplication.
    pub maximum_candidates: usize,
    pub maximum_contacts: usize,
    /// Retains separated potential contacts without relocating their points.
    pub activation_tolerance_m: f64,
    /// Minimum |det(edges / longest_edge)| for EACH current tetrahedron.
    pub minimum_scaled_jacobian: f64,
}
impl Default for FiniteSolidSurfaceControls {
    fn default() -> Self {
        Self {
            maximum_spheres: 64,
            maximum_query_visits: 1 << 20,
            maximum_candidates: 1 << 16,
            maximum_contacts: 21,
            activation_tolerance_m: 1e-10,
            minimum_scaled_jacobian: 1e-12,
        }
    }
}

/// An actual sphere witness paired with an original tetrahedral material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidSphereWitness {
    pub sphere_id: u64,
    /// Stable canonical exterior face used to address the material point.
    pub face: usize,
    pub contact: FiniteSolidRigidContact,
    pub rigid_point_world_m: Vec3,
    pub solid_point_world_m: Vec3,
    /// Signed normal separation of the two ACTUAL reconstructed witnesses.
    pub gap_m: f64,
}

/// Complete discrete-geometry result, ordered by stable sphere/contact IDs.
#[derive(Clone, Debug)]
pub struct FiniteSolidSurfaceQuery {
    /// One closest witness per sphere. Negative gaps remain negative.
    pub closest: Vec<FiniteSolidSphereWitness>,
    /// All near exterior features, with shared edge/vertex duplicates removed.
    /// These are raw geometry; call `admitted_contacts` before impact response.
    pub contacts: Vec<FiniteSolidSphereWitness>,
    pub visited_nodes: usize,
    pub tested_triangles: usize,
    pub tested_tetrahedra: usize,
    pub candidate_triangles: usize,
}
impl FiniteSolidSurfaceQuery {
    /// Admit nonpenetrating actual witnesses without tolerance-clamping gaps.
    pub fn admitted_contacts(
        &self,
    ) -> Result<Vec<FiniteSolidRigidContact>, FiniteSolidSurfaceError> {
        for witness in self.closest.iter().chain(&self.contacts) {
            if witness.gap_m < 0.0 {
                return Err(FiniteSolidSurfaceError::Penetration {
                    sphere_id: witness.sphere_id,
                    gap_m: witness.gap_m,
                });
            }
        }
        Ok(self.contacts.iter().map(|w| w.contact).collect())
    }
}

#[derive(Debug)]
pub enum FiniteSolidSurfaceError {
    Invalid(&'static str),
    ReferenceMismatch,
    Surface(HyperTetError),
    Query(TriangleQueryError),
    Solid(FiniteSolidError),
    Rigid(DynamicsError),
    DegenerateTetrahedron { tetrahedron: usize },
    InvertedTetrahedron { tetrahedron: usize },
    SphereCenterInside { sphere_id: u64, tetrahedron: usize },
    UnresolvedExterior { sphere_id: u64 },
    Penetration { sphere_id: u64, gap_m: f64 },
    WorkBudgetExhausted { maximum_query_visits: usize },
    ContactBudgetExhausted { maximum_contacts: usize },
    ContactIdOverflow { sphere_id: u64 },
    Cancelled,
}
impl core::fmt::Display for FiniteSolidSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solid sphere geometry refused: {self:?}")
    }
}
impl std::error::Error for FiniteSolidSurfaceError {}
impl From<HyperTetError> for FiniteSolidSurfaceError {
    fn from(error: HyperTetError) -> Self {
        Self::Surface(error)
    }
}
impl From<TriangleQueryError> for FiniteSolidSurfaceError {
    fn from(error: TriangleQueryError) -> Self {
        Self::Query(error)
    }
}
impl From<FiniteSolidError> for FiniteSolidSurfaceError {
    fn from(error: FiniteSolidError) -> Self {
        Self::Solid(error)
    }
}
impl From<DynamicsError> for FiniteSolidSurfaceError {
    fn from(error: DynamicsError) -> Self {
        Self::Rigid(error)
    }
}

#[derive(Debug)]
struct MaterialFace {
    tetrahedron: usize,
    local_nodes: [usize; 3],
    /// Canonical IDs for local edges (0,1), (0,2), (1,2).
    edge_ids: [usize; 3],
}

/// Reusable topology, material incidence and BVH for one exact reference mesh.
#[derive(Debug)]
pub struct FiniteSolidContactSurface {
    exterior: ClosedTetPressureSurface,
    triangle_surface: PreparedTriangleSurface,
    tetrahedra: Vec<[usize; 4]>,
    reference_orientation: Vec<f64>,
    material_faces: Vec<MaterialFace>,
    node_count: usize,
    edge_count: usize,
    feature_count: u64,
}

impl FiniteSolidContactSurface {
    /// Reuse native manifold/quality admission and cache each exterior face's
    /// original tet owner. `maximum_face_entries` bounds all incidence work by
    /// four entries per tet before allocation, through the original owner.
    pub fn new(
        solid: &FiniteSolidBody<'_>,
        maximum_face_entries: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteSolidSurfaceError> {
        poll(cx)?;
        let problem = solid.problem();
        let exterior = ClosedTetPressureSurface::new(
            problem.nodes_m,
            problem.tetrahedra,
            problem.budget,
            maximum_face_entries,
            cx,
        )?;
        let mut owners = BTreeMap::new();
        for (tetrahedron, tet) in problem.tetrahedra.iter().enumerate() {
            poll(cx)?;
            for opposite in 0..4 {
                let mut face = [0usize; 3];
                let mut slot = 0;
                for (local, &node) in tet.iter().enumerate() {
                    if local != opposite {
                        face[slot] = node;
                        slot += 1;
                    }
                }
                face.sort_unstable();
                owners.insert(face, tetrahedron);
            }
        }
        let mut edges = BTreeMap::new();
        for &nodes in exterior.faces() {
            poll(cx)?;
            for (a, b) in [(0, 1), (0, 2), (1, 2)] {
                edges.insert(edge(nodes[a], nodes[b]), 0usize);
            }
        }
        for (id, value) in edges.values_mut().enumerate() {
            *value = id;
        }
        let mut material_faces = Vec::with_capacity(exterior.faces().len());
        let mut triangles = Vec::with_capacity(exterior.faces().len());
        for &nodes in exterior.faces() {
            poll(cx)?;
            let mut key = nodes;
            key.sort_unstable();
            let tetrahedron = *owners.get(&key).ok_or(FiniteSolidSurfaceError::Invalid(
                "exterior lost its original tet owner",
            ))?;
            let tet = problem.tetrahedra[tetrahedron];
            let mut local_nodes = [0usize; 3];
            for (slot, node) in nodes.into_iter().enumerate() {
                local_nodes[slot] = tet.iter().position(|&id| id == node).ok_or(
                    FiniteSolidSurfaceError::Invalid("exterior node is absent from original tet"),
                )?;
            }
            material_faces.push(MaterialFace {
                tetrahedron,
                local_nodes,
                edge_ids: [(0, 1), (0, 2), (1, 2)].map(|(a, b)| edges[&edge(nodes[a], nodes[b])]),
            });
            let mut tri = [0u32; 3];
            for (slot, node) in nodes.into_iter().enumerate() {
                tri[slot] = u32::try_from(node).map_err(|_| {
                    FiniteSolidSurfaceError::Invalid("triangle vertex exceeds query index range")
                })?;
            }
            triangles.push(tri);
        }
        let triangle_surface = PreparedTriangleSurface::new(problem.nodes_m, &triangles, cx)?;
        let mut reference_orientation = Vec::with_capacity(problem.tetrahedra.len());
        for (tetrahedron, tet) in problem.tetrahedra.iter().enumerate() {
            poll(cx)?;
            let cell = CurrentTet::new(
                tet.map(|i| problem.nodes_m[i]),
                tetrahedron,
                None,
                problem.budget.minimum_scaled_jacobian,
            )?;
            reference_orientation.push(cell.determinant.signum());
        }
        let feature_count = u64::try_from(
            problem
                .nodes_m
                .len()
                .saturating_add(edges.len())
                .saturating_add(material_faces.len()),
        )
        .map_err(|_| FiniteSolidSurfaceError::Invalid("feature ID range exceeded"))?;
        Ok(Self {
            exterior,
            triangle_surface,
            tetrahedra: problem.tetrahedra.to_vec(),
            reference_orientation,
            material_faces,
            node_count: problem.nodes_m.len(),
            edge_count: edges.len(),
            feature_count,
        })
    }

    /// Original outward-wound exterior topology; indices address original nodes.
    #[must_use]
    pub fn faces(&self) -> &[[usize; 3]] {
        self.exterior.faces()
    }

    pub(crate) fn triangle_surface(&self) -> &PreparedTriangleSurface {
        &self.triangle_surface
    }

    pub(crate) fn reference_orientations(&self) -> &[f64] {
        &self.reference_orientation
    }

    pub(crate) fn current_positions(
        &self,
        solid: &FiniteSolidBody<'_>,
        cx: &Cx<'_>,
    ) -> Result<Vec<[f64; 3]>, FiniteSolidSurfaceError> {
        poll(cx)?;
        if !self
            .exterior
            .matches_reference(solid.problem().nodes_m, solid.problem().tetrahedra)
        {
            return Err(FiniteSolidSurfaceError::ReferenceMismatch);
        }
        Ok(solid.current_positions(cx)?)
    }

    /// Query current native coordinates. Raw negative gaps are retained so an
    /// event driver can bracket a crossing without accepting the overlap.
    pub fn sphere_separations(
        &self,
        solid: &FiniteSolidBody<'_>,
        rigid_bodies: &[RigidBodyState],
        spheres: &[FiniteSolidSphere],
        controls: FiniteSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteSolidSurfaceQuery, FiniteSolidSurfaceError> {
        let positions = self.current_positions(solid, cx)?;
        self.sphere_separations_at_positions(&positions, rigid_bodies, spheres, controls, cx)
    }

    /// Complete nonpenetrating event contacts for the existing full-mass/full-
    /// inertia response adapter. Positive activation never repairs penetration.
    pub fn sphere_contacts(
        &self,
        solid: &FiniteSolidBody<'_>,
        rigid_bodies: &[RigidBodyState],
        spheres: &[FiniteSolidSphere],
        controls: FiniteSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<Vec<FiniteSolidRigidContact>, FiniteSolidSurfaceError> {
        self.sphere_separations(solid, rigid_bodies, spheres, controls, cx)?
            .admitted_contacts()
    }

    /// A coordinate query for an explicitly declared trajectory image. It does
    /// not construct or claim a native accepted mechanical state.
    pub(crate) fn sphere_separations_at_positions(
        &self,
        positions: &[[f64; 3]],
        rigid_bodies: &[RigidBodyState],
        spheres: &[FiniteSolidSphere],
        controls: FiniteSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteSolidSurfaceQuery, FiniteSolidSurfaceError> {
        poll(cx)?;
        if positions.len() != self.node_count
            || spheres.len() > controls.maximum_spheres
            || !controls.activation_tolerance_m.is_finite()
            || controls.activation_tolerance_m < 0.0
            || !controls.minimum_scaled_jacobian.is_finite()
            || controls.minimum_scaled_jacobian <= 0.0
        {
            return Err(FiniteSolidSurfaceError::Invalid(
                "matching bounded geometry and positive finite quality required",
            ));
        }
        let cells = self.current_tetrahedra(positions, controls.minimum_scaled_jacobian, cx)?;
        let surface = self.triangle_surface.refit(positions, cx)?;
        let [volume_min, volume_max] = surface.bounds();
        let mut order: Vec<_> = (0..spheres.len()).collect();
        order.sort_unstable_by_key(|&i| spheres[i].id);
        let mut report = FiniteSolidSurfaceQuery {
            closest: Vec::with_capacity(spheres.len()),
            contacts: Vec::new(),
            visited_nodes: 0,
            tested_triangles: 0,
            tested_tetrahedra: 0,
            candidate_triangles: 0,
        };
        let mut contacts = BTreeMap::new();
        for (slot, &index) in order.iter().enumerate() {
            poll(cx)?;
            let sphere = spheres[index];
            if sphere.rigid_body >= rigid_bodies.len()
                || !sphere.center_body_m.is_finite()
                || !sphere.radius_m.is_finite()
                || sphere.radius_m <= 0.0
                || !sphere.friction_coefficient.is_finite()
                || sphere.friction_coefficient < 0.0
                || (slot > 0 && spheres[order[slot - 1]].id == sphere.id)
            {
                return Err(FiniteSolidSurfaceError::Invalid(
                    "distinct sphere IDs, valid rigid body, positive radius and nonnegative friction required",
                ));
            }
            sphere
                .id
                .checked_mul(self.feature_count)
                .and_then(|base| base.checked_add(self.feature_count - 1))
                .ok_or(FiniteSolidSurfaceError::ContactIdOverflow {
                    sphere_id: sphere.id,
                })?;
            let rigid = rigid_bodies[sphere.rigid_body];
            let center = rigid.pose().point_world_from_body(sphere.center_body_m)?;
            let center_array = array(center);
            // Most spheres approaching an exterior face are outside the whole
            // volume's box. They need no per-element containment scan. Centres
            // in that box still receive the complete bounded volume test.
            let outside_volume_box = (0..3).any(|axis| {
                center_array[axis] < volume_min[axis] || center_array[axis] > volume_max[axis]
            });
            for (tetrahedron, cell) in
                cells
                    .iter()
                    .enumerate()
                    .take(if outside_volume_box { 0 } else { cells.len() })
            {
                poll(cx)?;
                if report.work() >= controls.maximum_query_visits {
                    return Err(FiniteSolidSurfaceError::WorkBudgetExhausted {
                        maximum_query_visits: controls.maximum_query_visits,
                    });
                }
                report.tested_tetrahedra += 1;
                if cell.contains(center_array)? {
                    return Err(FiniteSolidSurfaceError::SphereCenterInside {
                        sphere_id: sphere.id,
                        tetrahedron,
                    });
                }
            }
            let within = finite(sphere.radius_m + controls.activation_tolerance_m)?.next_up();
            let query = surface.point_query(
                center_array,
                within,
                controls.maximum_query_visits.saturating_sub(report.work()),
                controls
                    .maximum_candidates
                    .saturating_sub(report.candidate_triangles),
                cx,
            )?;
            report.visited_nodes += query.visited_nodes;
            report.tested_triangles += query.tested_triangles;
            report.candidate_triangles += query.nearby.len();
            let closest = self.witness(
                positions,
                sphere,
                rigid,
                center,
                query.nearest_triangle,
                query.nearest,
                controls.activation_tolerance_m,
            )?;
            let side = closest
                .contact
                .normal_world
                .dot(vec(query.nearest.unit_normal));
            if side < -128.0 * f64::EPSILON {
                return Err(FiniteSolidSurfaceError::UnresolvedExterior {
                    sphere_id: sphere.id,
                });
            }
            report.closest.push(closest);
            for (face, point) in query.nearby {
                let witness = self.witness(
                    positions,
                    sphere,
                    rigid,
                    center,
                    face,
                    point,
                    controls.activation_tolerance_m,
                )?;
                if witness.gap_m > controls.activation_tolerance_m
                    || witness.contact.normal_world.dot(vec(point.unit_normal))
                        < -128.0 * f64::EPSILON
                {
                    continue;
                }
                if !contacts.contains_key(&witness.contact.id) {
                    if contacts.len() == controls.maximum_contacts {
                        return Err(FiniteSolidSurfaceError::ContactBudgetExhausted {
                            maximum_contacts: controls.maximum_contacts,
                        });
                    }
                    contacts.insert(witness.contact.id, witness);
                }
            }
        }
        report.contacts = contacts.into_values().collect();
        poll(cx)?;
        Ok(report)
    }

    fn current_tetrahedra(
        &self,
        positions: &[[f64; 3]],
        minimum_scaled_jacobian: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<CurrentTet>, FiniteSolidSurfaceError> {
        let mut cells = Vec::with_capacity(self.tetrahedra.len());
        for (tetrahedron, (&tet, &orientation)) in self
            .tetrahedra
            .iter()
            .zip(&self.reference_orientation)
            .enumerate()
        {
            poll(cx)?;
            cells.push(CurrentTet::new(
                tet.map(|i| positions[i]),
                tetrahedron,
                Some(orientation),
                minimum_scaled_jacobian,
            )?);
        }
        Ok(cells)
    }

    #[allow(clippy::too_many_arguments)]
    fn witness(
        &self,
        positions: &[[f64; 3]],
        sphere: FiniteSolidSphere,
        rigid: RigidBodyState,
        center: Vec3,
        face: usize,
        point: TrianglePoint,
        activation_tolerance_m: f64,
    ) -> Result<FiniteSolidSphereWitness, FiniteSolidSurfaceError> {
        let material = &self.material_faces[face];
        let mut weights = [0.0; 4];
        for (local, weight) in material.local_nodes.into_iter().zip(point.barycentric) {
            weights[local] = weight;
        }
        let solid_point = FiniteSolidPoint::new(material.tetrahedron, weights)?;
        let mut position = [0.0; 3];
        // Match the original material-point trace's node order and multiply-add.
        for (node, weight) in self.tetrahedra[material.tetrahedron]
            .into_iter()
            .zip(solid_point.weights())
        {
            for axis in 0..3 {
                position[axis] = weight.mul_add(positions[node][axis], position[axis]);
            }
        }
        let solid_point_world_m = vec(vector(position)?);
        let radial = center.sub(solid_point_world_m);
        let distance = magnitude(radial)?;
        if distance == 0.0 {
            return Err(FiniteSolidSurfaceError::UnresolvedExterior {
                sphere_id: sphere.id,
            });
        }
        let normal_world = Vec3::new(
            radial.x / distance,
            radial.y / distance,
            radial.z / distance,
        );
        let arm_body_m = sphere.center_body_m.sub(
            rigid
                .pose()
                .orientation()
                .rotate_world_to_body(normal_world.scale(sphere.radius_m)),
        );
        let rigid_point_world_m = rigid.pose().point_world_from_body(arm_body_m)?;
        let gap_m = finite(normal_world.dot(rigid_point_world_m.sub(solid_point_world_m)))?;
        let feature = match point.feature {
            TriangleFeature::Vertex(local) => self.faces()[face][local],
            TriangleFeature::Edge(a, b) => {
                let local_edge = match (a, b) {
                    (0, 1) => 0,
                    (0, 2) => 1,
                    (1, 2) => 2,
                    _ => {
                        return Err(FiniteSolidSurfaceError::Invalid(
                            "unknown local triangle edge",
                        ));
                    }
                };
                self.node_count + material.edge_ids[local_edge]
            }
            TriangleFeature::Face => self.node_count + self.edge_count + face,
        };
        let feature =
            u64::try_from(feature).map_err(|_| FiniteSolidSurfaceError::ContactIdOverflow {
                sphere_id: sphere.id,
            })?;
        let id = sphere
            .id
            .checked_mul(self.feature_count)
            .and_then(|base| base.checked_add(feature))
            .ok_or(FiniteSolidSurfaceError::ContactIdOverflow {
                sphere_id: sphere.id,
            })?;
        let components = [
            normal_world.x.abs(),
            normal_world.y.abs(),
            normal_world.z.abs(),
        ];
        let axis = if components[0] <= components[1] && components[0] <= components[2] {
            Vec3::new(1.0, 0.0, 0.0)
        } else if components[1] <= components[2] {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let tangent = normal_world.cross(axis);
        let tangent_world = tangent.scale(magnitude(tangent)?.recip());
        Ok(FiniteSolidSphereWitness {
            sphere_id: sphere.id,
            face,
            rigid_point_world_m,
            solid_point_world_m,
            gap_m,
            contact: FiniteSolidRigidContact {
                id,
                rigid_body: sphere.rigid_body,
                arm_body_m,
                solid_point,
                normal_world,
                tangent_world,
                activation_tolerance_m,
                friction_coefficient: sphere.friction_coefficient,
                normal_law: sphere.normal_law,
            },
        })
    }
}

impl FiniteSolidSurfaceQuery {
    fn work(&self) -> usize {
        self.visited_nodes
            .saturating_add(self.tested_triangles)
            .saturating_add(self.tested_tetrahedra)
    }
}

struct CurrentTet {
    origin: [f64; 3],
    min: [f64; 3],
    max: [f64; 3],
    scale: f64,
    cofactors: [[f64; 3]; 3],
    determinant: f64,
}
impl CurrentTet {
    fn new(
        points: [[f64; 3]; 4],
        tetrahedron: usize,
        orientation: Option<f64>,
        floor: f64,
    ) -> Result<Self, FiniteSolidSurfaceError> {
        let mut scale = 0.0_f64;
        for p in points {
            vector(p)?;
        }
        for a in 0..4 {
            for b in a + 1..4 {
                scale = scale.max(magnitude(vec(difference(points[a], points[b])?))?);
            }
        }
        if scale == 0.0 {
            return Err(FiniteSolidSurfaceError::DegenerateTetrahedron { tetrahedron });
        }
        let a = difference(points[1], points[0])?.map(|v| v / scale);
        let b = difference(points[2], points[0])?.map(|v| v / scale);
        let c = difference(points[3], points[0])?.map(|v| v / scale);
        let cofactors = [cross(b, c), cross(c, a), cross(a, b)];
        let determinant = finite(dot(a, cofactors[0]))?;
        if orientation.is_some_and(|sign| sign * determinant < 0.0) {
            return Err(FiniteSolidSurfaceError::InvertedTetrahedron { tetrahedron });
        }
        if determinant.abs() < floor {
            return Err(FiniteSolidSurfaceError::DegenerateTetrahedron { tetrahedron });
        }
        let min = core::array::from_fn(|axis| {
            points
                .iter()
                .fold(f64::INFINITY, |v, p| v.min(p[axis]))
                .next_down()
        });
        let max = core::array::from_fn(|axis| {
            points
                .iter()
                .fold(f64::NEG_INFINITY, |v, p| v.max(p[axis]))
                .next_up()
        });
        vector(min)?;
        vector(max)?;
        Ok(Self {
            origin: points[0],
            min,
            max,
            scale,
            cofactors,
            determinant,
        })
    }
    fn contains(&self, point: [f64; 3]) -> Result<bool, FiniteSolidSurfaceError> {
        if (0..3).any(|axis| point[axis] < self.min[axis] || point[axis] > self.max[axis]) {
            return Ok(false);
        }
        let offset = difference(point, self.origin)?.map(|v| v / self.scale);
        let mut weights = [0.0; 4];
        for i in 0..3 {
            weights[i + 1] = finite(dot(self.cofactors[i], offset) / self.determinant)?;
        }
        weights[0] = finite(1.0 - weights[1] - weights[2] - weights[3])?;
        // A centre on an unresolved tet boundary is unsafe for a positive-radius
        // sphere too. Do not turn small negative barycentric roundoff into an
        // asserted exterior clearance.
        let coordinate_scale = offset.iter().fold(1.0_f64, |s, x| s.max(x.abs()));
        let guard = finite(128.0 * f64::EPSILON * coordinate_scale / self.determinant.abs())?;
        Ok(weights.into_iter().all(|w| w >= -guard && w <= 1.0 + guard))
    }
}

fn edge(a: usize, b: usize) -> [usize; 2] {
    [a.min(b), a.max(b)]
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vec(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn difference(a: [f64; 3], b: [f64; 3]) -> Result<[f64; 3], FiniteSolidSurfaceError> {
    vector(core::array::from_fn(|i| a[i] - b[i]))
}
fn vector(v: [f64; 3]) -> Result<[f64; 3], FiniteSolidSurfaceError> {
    for x in v {
        finite(x)?;
    }
    Ok(v)
}
fn finite(v: f64) -> Result<f64, FiniteSolidSurfaceError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(FiniteSolidSurfaceError::Invalid(
            "geometry arithmetic is not representable",
        ))
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}
fn magnitude(v: Vec3) -> Result<f64, FiniteSolidSurfaceError> {
    finite(v.x.hypot(v.y).hypot(v.z))
}
fn poll(cx: &Cx<'_>) -> Result<(), FiniteSolidSurfaceError> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidSurfaceError::Cancelled)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn current_tet_geometry_refuses_inversion_collapse_and_ambiguous_interior() {
        let points = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let cell = CurrentTet::new(points, 7, Some(1.0), 1e-12).unwrap();
        assert!(cell.contains([0.125, 0.25, 0.5]).unwrap());
        assert!(cell.contains([0.25, 0.25, 0.0]).unwrap());
        assert!(!cell.contains([0.5, 0.5, 0.5]).unwrap());
        assert!(!cell.contains([0.25, 0.25, -0.125]).unwrap());
        let reversed = [points[0], points[2], points[1], points[3]];
        assert!(matches!(
            CurrentTet::new(reversed, 7, Some(1.0), 1e-12),
            Err(FiniteSolidSurfaceError::InvertedTetrahedron { tetrahedron: 7 })
        ));
        assert!(CurrentTet::new(reversed, 7, Some(-1.0), 1e-12).is_ok());
        let flat = [points[0], points[1], points[2], [0.0, 0.0, 1e-14]];
        assert!(matches!(
            CurrentTet::new(flat, 7, Some(1.0), 1e-12),
            Err(FiniteSolidSurfaceError::DegenerateTetrahedron { tetrahedron: 7 })
        ));
    }
}
