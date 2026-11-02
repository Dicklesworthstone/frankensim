//! Immutable geometry for continuous P2 displacement / P1 pressure tetrahedra.
//!
//! The original Tet10 owner supplies the basis, map, Jacobian and shared-edge
//! admission. Duffy quadrature is explicit: at least four points per collapsed
//! direction integrate pressure-test times the cubic affine-element det(F).
//! Nonlinear energies and curved geometry still require a quadrature ladder.
//! Positive sampled Jacobians do not certify injectivity between those samples.

use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;

use crate::linear3::tet10::{MAX_QUADRATURE_POINTS_PER_DIR, tet10_point, validate_tet10_topology};
use crate::linear3::{TetAssemblyBudget, TetElasticError};

#[derive(Clone, Debug)]
pub(super) struct MixedPoint {
    pub(super) weight_m3: f64,
    /// Row zero is balanced against the other rows so the exact derivative of
    /// relative-displacement kinematics also annihilates rigid translation.
    pub(super) gradients: [[f64; 3]; 10],
    pub(super) shapes: [f64; 10],
    pub(super) pressure_shapes: [f64; 4],
}

#[derive(Clone, Debug)]
pub(super) struct PreparedMixedElement {
    pub(super) connectivity: [usize; 10],
    pub(super) pressure_indices: [usize; 4],
    pub(super) points: Vec<MixedPoint>,
    pub(super) reference_volume_m3: f64,
}

#[derive(Clone, Debug)]
pub(super) struct PreparedMixedMesh {
    pub(super) nodes_count: usize,
    pub(super) elements: Vec<PreparedMixedElement>,
    /// Sorted original corner-node IDs; each owns ONE shared pressure unknown.
    pub(super) pressure_nodes: Vec<usize>,
    pub(super) points_per_element: usize,
    pub(super) total_volume_m3: f64,
}

impl PreparedMixedMesh {
    /// Prepare only reference geometry. The mixed owner separately admits its
    /// constrained displacement/pressure DOFs and complete KKT storage budget.
    pub(super) fn new(
        nodes_m: &[[f64; 3]],
        elements: &[[usize; 10]],
        budget: TetAssemblyBudget,
        quadrature_points_per_dir: usize,
        maximum_quadrature_points: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        if nodes_m.is_empty() || elements.is_empty() {
            return Err(TetElasticError::EmptyMesh);
        }
        if !(4..=MAX_QUADRATURE_POINTS_PER_DIR).contains(&quadrature_points_per_dir) {
            return Err(TetElasticError::QuadratureInadmissible {
                what: "mixed P2/P1 Duffy points per direction must be in 4..=12",
            });
        }
        let points_per_element = quadrature_points_per_dir.pow(3);
        let total_points = elements
            .len()
            .checked_mul(points_per_element)
            .ok_or(TetElasticError::SizeOverflow)?;
        for (what, requested, maximum) in [
            ("nodes", nodes_m.len(), budget.maximum_nodes),
            ("tetrahedra", elements.len(), budget.maximum_tetrahedra),
            (
                "mixed_quadrature_points",
                total_points,
                maximum_quadrature_points,
            ),
        ] {
            if requested > maximum {
                return Err(TetElasticError::BudgetExceeded {
                    what,
                    requested,
                    maximum,
                });
            }
        }
        for (node, coordinates) in nodes_m.iter().enumerate() {
            if node % 256 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            if coordinates.iter().any(|value| !value.is_finite()) {
                return Err(TetElasticError::NonFiniteNode { node });
            }
        }
        let pressure_nodes: Vec<_> = validate_tet10_topology(nodes_m.len(), elements, cx)?
            .into_iter()
            .collect();
        nodes_m
            .len()
            .checked_mul(3)
            .and_then(|dofs| dofs.checked_add(pressure_nodes.len()))
            .ok_or(TetElasticError::SizeOverflow)?;
        let quadrature = duffy_quadrature(quadrature_points_per_dir);
        let mut prepared = Vec::with_capacity(elements.len());
        let mut total_volume_m3 = 0.0;
        for (element, &connectivity) in elements.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let positions = connectivity.map(|node| nodes_m[node]);
            let pressure_indices = core::array::from_fn(|a| {
                pressure_nodes
                    .binary_search(&connectivity[a])
                    .expect("admitted corner pressure node")
            });
            let mut points = Vec::with_capacity(points_per_element);
            let mut reference_volume_m3 = 0.0;
            for (q, &(lambda, reference_weight)) in quadrature.iter().enumerate() {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                let point = tet10_point(
                    positions,
                    lambda,
                    budget.minimum_scaled_jacobian,
                    element,
                    q,
                )?;
                let weight_m3 = reference_weight * point.determinant_m3;
                if !weight_m3.is_finite() || weight_m3 <= 0.0 {
                    return Err(TetElasticError::QuadratureInadmissible {
                        what: "mixed reference quadrature weight must be finite and positive",
                    });
                }
                let mut gradients = point.gradients;
                gradients[0] =
                    core::array::from_fn(|c| -gradients[1..].iter().map(|g| g[c]).sum::<f64>());
                if gradients.iter().flatten().any(|value| !value.is_finite()) {
                    return Err(TetElasticError::QuadratureInadmissible {
                        what: "mixed reference gradients are unrepresentable",
                    });
                }
                reference_volume_m3 += weight_m3;
                points.push(MixedPoint {
                    weight_m3,
                    gradients,
                    shapes: point.shapes,
                    pressure_shapes: lambda,
                });
            }
            total_volume_m3 += reference_volume_m3;
            if !total_volume_m3.is_finite() || total_volume_m3 <= 0.0 {
                return Err(TetElasticError::QuadratureInadmissible {
                    what: "mixed total reference volume must be finite and positive",
                });
            }
            prepared.push(PreparedMixedElement {
                connectivity,
                pressure_indices,
                points,
                reference_volume_m3,
            });
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(Self {
            nodes_count: nodes_m.len(),
            elements: prepared,
            pressure_nodes,
            points_per_element,
            total_volume_m3,
        })
    }
}
