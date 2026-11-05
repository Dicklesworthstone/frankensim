//! Explicit tetrahedral through-seams and conservative physical-state transfer.
//!
//! Selected internal faces become paired exterior traces. Node duplication
//! follows connected tetrahedral stars, so every unselected face keeps its
//! original continuity. A selected face that remains tied at a crack-front
//! vertex refuses: the existing matching P1 cohesive owner requires distinct
//! degrees of freedom on both traces. This is a through-seam preparation path,
//! not automatic crack discovery, remeshing, or a bulk-damage-to-surface-law
//! conversion. In particular, inserting a fresh cohesive law into an already
//! loaded solid needs a separate traction/history/energy transfer admission.
//!
//! Original elements, reference coordinates and integration-point order do not
//! change. Motion is duplicated exactly, and its consistent P1 mass, momentum,
//! angular momentum and kinetic energy are reconstructed on both meshes before
//! publication. Source bindings and element histories can therefore retain
//! their original ordered slices; this does not rebind opaque solver-owned
//! checkpoints to a new system. Existing cohesive dynamics owns all subsequent
//! fracture loss, motion and compressive reclosure.

use crate::cohesive::p1::FacePair;
use crate::loads3::TetLoadGeometry3;
use fs_exec::Cx;
use fs_material::cohesive::BilinearCohesive;

mod mesh;
pub use mesh::{
    FractureMotionTransfer, FractureTopologyBudget, FractureTopologyError, InteriorSeam,
    P1PhysicalInventory, SplitFace, TetFractureTopology,
};

impl<'g> TetFractureTopology<'g> {
    /// Split explicitly selected internal faces of admitted reference geometry.
    /// All allocations and face work are bounded before construction; the
    /// original geometry and all caller state remain unchanged on refusal.
    pub fn new(
        geometry: &'g TetLoadGeometry3<'_>,
        seams: &[InteriorSeam],
        budget: FractureTopologyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FractureTopologyError> {
        Self::new_controlled(geometry, seams, budget, &mut || cx.checkpoint().is_err())
    }

    /// Predicate form with identical topology, work bounds and cancellation.
    pub fn new_controlled<C: FnMut() -> bool>(
        geometry: &'g TetLoadGeometry3<'_>,
        seams: &[InteriorSeam],
        budget: FractureTopologyBudget,
        cancelled: &mut C,
    ) -> Result<Self, FractureTopologyError> {
        Self::prepare(
            geometry.nodes_m(),
            geometry.tetrahedra(),
            geometry.element_volumes_m3(),
            seams,
            budget,
            cancelled,
        )
    }

    /// Duplicate intensive nodal motion and reconstruct consistent P1 inertia.
    /// Densities follow the unchanged element order and have units kg/m³.
    pub fn transfer_motion(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        density_kg_m3: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FractureMotionTransfer, FractureTopologyError> {
        self.transfer_motion_controlled(displacement_m, velocity_m_s, density_kg_m3, &mut || {
            cx.checkpoint().is_err()
        })
    }
}

impl SplitFace {
    /// Use the existing matching-trace quadrature and cohesive mechanics owner.
    /// The caller's explicit A-side and interface identity survive splitting.
    /// This only binds geometry; it does not create or reset material history.
    #[must_use]
    pub fn face_pair<'a>(&self, law: &'a BilinearCohesive) -> FacePair<'a> {
        FacePair {
            id: self.id,
            side_a: self.side_a,
            side_b: self.side_b,
            law,
        }
    }
}
