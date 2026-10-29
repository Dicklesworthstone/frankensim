//! Original P1 exterior geometry retained beside a complete material image.
use super::{ReduceError, ValveCardRequest, ValveElementOrder, invalid};
use crate::linear3::TetAssemblyBudget;
use crate::loads3::{LoadError3, LoadFace3, TetLoadGeometry3};
use fs_exec::Cx;

/// Immutable source nodes and original exterior triangles. These are the same
/// reference geometry as the volume material, supports and pressure ports.
/// Retaining them once lets downstream contact use actual physical areas and
/// nodal traces without inventing an effective face or rebuilding incidence.
#[derive(Debug)]
pub struct MaterialValveSurface {
    nodes_m: Vec<[f64; 3]>,
    faces: Vec<LoadFace3>,
}
impl MaterialValveSurface {
    pub fn nodes_m(&self) -> &[[f64; 3]] {
        &self.nodes_m
    }
    /// Canonical boundary-face order from `TetLoadGeometry3`; never interior
    /// faces. Original outward normals are independent of tet winding.
    pub fn faces(&self) -> &[LoadFace3] {
        &self.faces
    }

    pub(super) fn from_request<M>(
        request: &ValveCardRequest<'_, M>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Option<Self>, ReduceError> {
        if !matches!(request.element_order, ValveElementOrder::Tet4) {
            // The P1 boundary owner cannot represent curved quadratic traces.
            return Ok(None);
        }
        let cap = budget
            .maximum_tetrahedra
            .checked_mul(4)
            .ok_or_else(|| invalid("original exterior face budget overflow"))?;
        let geometry = TetLoadGeometry3::new(request.nodes_m, request.tetrahedra, budget, cap, cx)
            .map_err(|error| match error {
                LoadError3::Cancelled => crate::linear3::TetElasticError::Cancelled.into(),
                LoadError3::Solver(error) => error.into(),
                _ => invalid("original material valve exterior geometry is inadmissible"),
            })?;
        super::poll(cx)?;
        Ok(Some(Self {
            nodes_m: request.nodes_m.to_vec(),
            faces: geometry.faces().to_vec(),
        }))
    }
}
