//! Actual quadratic tetrahedral boundaries feeding native simultaneous impacts.
//!
//! The original solid topology owner admits conforming TET10 incidence. A compact
//! corner mesh supplies only closed-boundary incidence and reference winding;
//! every current distance, normal and material witness uses the six-node curved
//! trace. Signed TET10 shape weights then feed the original native mass/contact
//! owner, including its pressure constraints and retained material histories.
//!
//! Detection is a bounded discrete query. Distance intervals and localized
//! minimizers must resolve nonpenetration and the authored activation band before
//! mechanical admission. Backside or unresolved witnesses refuse. A positive
//! Jacobian at the selected material point is a local admission check, not a
//! certificate of global injectivity: the source mesh must be a nonoverlapping,
//! consistently oriented solid. There is no curved-volume corner containment
//! approximation, penetration repair or time integrator. The path module
//! separately bounds declared affine coordinate motion; it does not certify a
//! material trajectory between integration endpoints.

use crate::finite_solid::contact::surface::{
    FiniteSolidSphere, FiniteSolidSurfaceControls, FiniteSolidSurfaceError,
};
use crate::inelastic_impact::multiple::InelasticMixedContact;
use fs_mbd::{DynamicsError, Vec3};
use fs_query::triangle::quadratic::{
    PreparedQuadraticTriangleSurface, QuadraticQueryControls, QuadraticQueryError,
};
use fs_solid::{
    hyper3::{HyperTetError, mixed::visco::dynamics::Tet10ImpulsePoint},
    linear3::TetElasticError,
};

mod geometry;
mod topology;

/// Native circular and pressure-constrained mechanical adapters.
pub mod mechanical;

/// Whole-path bounds for explicitly declared affine nodal and sphere motion.
pub mod path;

/// Bounded geometry controls; all distances are in metres.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticSolidSurfaceControls {
    /// Sphere, traversal, candidate, contact and local Jacobian admission limits.
    /// Traversal work counts BVH visits plus tested patches, across all spheres.
    pub surface: FiniteSolidSurfaceControls,
    /// Curved subdivision limits. `maximum_cells` is aggregate across spheres;
    /// the active-cell limit bounds one live branch-and-bound frontier.
    pub query: QuadraticQueryControls,
    /// Maximum physical distance from the returned witness to any possible
    /// minimizing point. The query's witness tolerance is tightened to this
    /// value; distinct unresolved minima cannot become a claimed single contact.
    pub maximum_witness_span_m: f64,
}

impl Default for QuadraticSolidSurfaceControls {
    fn default() -> Self {
        Self {
            surface: FiniteSolidSurfaceControls::default(),
            query: QuadraticQueryControls {
                distance_tolerance: 1.0e-10,
                witness_tolerance: 1.0e-8,
                maximum_cells: 1 << 20,
                maximum_active_cells: 1 << 16,
            },
            maximum_witness_span_m: 1.0e-8,
        }
    }
}

/// Real radial sphere point and original quadratic material trace.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticSolidSphereWitness {
    pub sphere_id: u64,
    pub face: usize,
    pub contact: InelasticMixedContact<Tet10ImpulsePoint>,
    pub rigid_point_world_m: Vec3,
    pub solid_point_world_m: Vec3,
    /// Separation of the actual reconstructed witnesses; may be negative.
    pub gap_m: f64,
    /// Conservative curved-distance gap interval, including reconstruction error.
    pub gap_lower_m: f64,
    pub gap_upper_m: f64,
}

/// Raw bounded detection receipt. Admission is deliberately a separate operation.
#[derive(Debug)]
pub struct QuadraticSolidSurfaceQuery {
    /// One globally nearest exterior witness per sphere, sorted by sphere ID.
    pub closest: Vec<QuadraticSolidSphereWitness>,
    /// Canonical resolved active features, sorted by contact ID and deduplicated.
    pub contacts: Vec<QuadraticSolidSphereWitness>,
    /// First conservative candidate whose distance interval crosses the authored
    /// activation boundary. No potentially active feature is silently discarded.
    pub unresolved_activation: Option<QuadraticSolidSphereWitness>,
    pub visited_nodes: usize,
    pub tested_triangles: usize,
    pub examined_cells: usize,
    pub candidate_triangles: usize,
}

impl QuadraticSolidSurfaceQuery {
    /// Admit only resolved nonpenetrating geometry for the unchanged mechanical
    /// owner. A raw negative gap remains available to an external event locator.
    pub fn admitted_contacts(
        &self,
    ) -> Result<Vec<InelasticMixedContact<Tet10ImpulsePoint>>, QuadraticSolidSurfaceError> {
        for witness in self.closest.iter().chain(&self.contacts) {
            if witness.gap_m < 0.0 {
                return Err(FiniteSolidSurfaceError::Penetration {
                    sphere_id: witness.sphere_id,
                    gap_m: witness.gap_m,
                }
                .into());
            }
            if witness.gap_lower_m < 0.0 {
                return Err(QuadraticSolidSurfaceError::UnresolvedGap {
                    sphere_id: witness.sphere_id,
                    lower_gap_m: witness.gap_lower_m,
                    upper_gap_m: witness.gap_upper_m,
                });
            }
        }
        if let Some(witness) = self.unresolved_activation {
            return Err(QuadraticSolidSurfaceError::UnresolvedActivation {
                sphere_id: witness.sphere_id,
                face: witness.face,
                lower_gap_m: witness.gap_lower_m,
                upper_gap_m: witness.gap_upper_m,
                activation_tolerance_m: witness.contact.activation_tolerance_m,
            });
        }
        Ok(self
            .contacts
            .iter()
            .map(|witness| witness.contact)
            .collect())
    }
}

/// Preserve original topology/solid/rigid refusals and bounded curved-query evidence.
#[derive(Debug)]
pub enum QuadraticSolidSurfaceError {
    Surface(FiniteSolidSurfaceError),
    SolidGeometry(TetElasticError),
    Query(QuadraticQueryError),
    UnresolvedGap {
        sphere_id: u64,
        lower_gap_m: f64,
        upper_gap_m: f64,
    },
    UnresolvedActivation {
        sphere_id: u64,
        face: usize,
        lower_gap_m: f64,
        upper_gap_m: f64,
        activation_tolerance_m: f64,
    },
    AmbiguousWitness {
        sphere_id: u64,
        face: usize,
        error_bound_m: f64,
        maximum_error_m: f64,
    },
}

impl core::fmt::Display for QuadraticSolidSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quadratic solid surface refused: {self:?}")
    }
}
impl std::error::Error for QuadraticSolidSurfaceError {}
impl From<FiniteSolidSurfaceError> for QuadraticSolidSurfaceError {
    fn from(error: FiniteSolidSurfaceError) -> Self {
        Self::Surface(error)
    }
}
impl From<TetElasticError> for QuadraticSolidSurfaceError {
    fn from(error: TetElasticError) -> Self {
        Self::SolidGeometry(error)
    }
}
impl From<QuadraticQueryError> for QuadraticSolidSurfaceError {
    fn from(error: QuadraticQueryError) -> Self {
        Self::Query(error)
    }
}
impl From<HyperTetError> for QuadraticSolidSurfaceError {
    fn from(error: HyperTetError) -> Self {
        Self::Surface(FiniteSolidSurfaceError::Surface(error))
    }
}
impl From<DynamicsError> for QuadraticSolidSurfaceError {
    fn from(error: DynamicsError) -> Self {
        Self::Surface(FiniteSolidSurfaceError::Rigid(error))
    }
}

#[derive(Debug)]
struct MaterialFace {
    tetrahedron: usize,
    local_corners: [usize; 3],
    /// Canonical exterior edges in local `(0,1),(0,2),(1,2)` order.
    edge_ids: [usize; 3],
}

/// Immutable original TET10 identity, exterior incidence and reusable curved BVH.
#[derive(Debug)]
pub struct QuadraticSolidContactSurface {
    reference_nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 10]>,
    triangle_surface: PreparedQuadraticTriangleSurface,
    material_faces: Vec<MaterialFace>,
    edge_count: usize,
    feature_count: u64,
}

fn checkpoint(cx: &fs_exec::Cx<'_>) -> Result<(), QuadraticSolidSurfaceError> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidSurfaceError::Cancelled.into())
}
fn invalid(message: &'static str) -> QuadraticSolidSurfaceError {
    FiniteSolidSurfaceError::Invalid(message).into()
}
fn finite(value: f64) -> Result<f64, QuadraticSolidSurfaceError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("finite quadratic contact geometry required"))
    }
}
fn vector(value: [f64; 3]) -> Result<[f64; 3], QuadraticSolidSurfaceError> {
    for coordinate in value {
        finite(coordinate)?;
    }
    Ok(value)
}
fn edge(a: usize, b: usize) -> [usize; 2] {
    [a.min(b), a.max(b)]
}

#[cfg(test)]
mod tests;
