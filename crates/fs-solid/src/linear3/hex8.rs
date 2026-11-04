//! Trilinear three-dimensional elasticity with consistent reference inertia.
//!
//! The original element-constant Mandel tensor enters the shared matrix
//! assembler unchanged. B-bar replaces only the volumetric part of each
//! strain column by its physical-volume-weighted element mean:
//! epsilon_bar = dev(epsilon) + mean(trace(epsilon)) I / 3.
//! Forces and storage therefore derive from the same full B_bar^T C B_bar,
//! including anisotropic normal/shear coupling. This is a small-strain
//! assumed-strain brick, not an exact incompressibility constraint.
//!
//! Full deviatoric quadrature is retained, without hourglass stiffness.
//! Stiffness uses a caller-selected tensor Gauss rule of order 2..=12;
//! distorted maps still require quadrature/mesh convergence. Consistent mass
//! uses at least three points per direction: N_a N_b det(J) can be quartic
//! in each coordinate on a trilinear map, so two points are insufficient.
//! Mass is independent of the strain formulation and is never lumped.
//!
//! Source-integrated elementwise thermal expansion enters the same projected
//! weak form. Thermal equilibrium retains full nodal support reactions and
//! integrates stress and stored energy from that same assumed-strain field.
//! Temperature gradients within an element require a separately admitted
//! pressure/eigenstrain projection and are not silently averaged here.
//!
//! The returned ordinary stiffness/mass pencil feeds the existing spectral
//! and linear-dynamics owners. Thermal or dead-load equilibrium also feeds the
//! shared signed initial-stress tangent, retaining the same strain formulation
//! and reference mass. Supports are homogeneous displacement DOFs.
//! Geometry admission samples corners, centre and integration points; this
//! does not certify global injectivity of a warped brick. No finite-strain,
//! general shear-locking or continuum-convergence claim is made.

use super::hex::{HEX27_REFERENCE_NODES, hex8_point};
use super::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetMaterialField, isoparametric,
};
use crate::linear::Formulation;
use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
use std::collections::BTreeSet;

pub use super::{
    QuadraticElementThermalResponse as Hex8ElementThermalResponse,
    QuadraticPrestressedAssembly as Hex8PrestressedAssembly,
    QuadraticThermalResponse as Hex8ThermalResponse,
};

mod thermal;

/// A conforming eight-node small-strain solid on its reference geometry.
#[derive(Debug, Clone, Copy)]
pub struct Hex8ElasticProblem<'a> {
    /// Reference coordinates in metres.
    pub nodes_m: &'a [[f64; 3]],
    /// Canonical corners: four around z-, then corresponding corners at z+.
    /// This is the first eight entries of HEX27_REFERENCE_NODES.
    pub hexahedra: &'a [[usize; 8]],
    /// Original full stiffness and positive density, uniform or per element.
    pub materials: TetMaterialField<'a>,
    /// Standard compatible strain or three-dimensional volumetric B-bar.
    pub formulation: Formulation,
    /// Homogeneous supported displacement DOFs, in full node-major xyz order.
    pub fixed_dofs: &'a [usize],
    /// Existing work/geometry envelope; maximum_tetrahedra bounds brick count.
    pub budget: TetAssemblyBudget,
    /// Stiffness Gauss points per direction, in 2..=12. Use 2 for the standard
    /// fully integrated brick and increase it to check distorted-map accuracy.
    /// Mass uses max(3, stiffness_points_per_dir).
    pub stiffness_points_per_dir: usize,
}

impl Hex8ElasticProblem<'_> {
    /// Assemble the original material stiffness and consistent mass. All six
    /// free-body rigid modes are allowed; downstream owners decide which
    /// constraints a particular static or dynamic operation requires.
    pub fn assemble(&self, cx: &Cx<'_>) -> Result<TetElasticAssembly, TetElasticError> {
        let minimum_probe = self.validate(true, cx)?;
        let stiffness_rule = tensor_rule(self.stiffness_points_per_dir);
        let extra_mass_rule = (self.stiffness_points_per_dir < 3).then(|| tensor_rule(3));
        let mass_rule = extra_mass_rule.as_deref().unwrap_or(&stiffness_rule);
        let mut assembly = isoparametric::assemble_with_formulation::<8, 24, _>(
            self.nodes_m,
            self.hexahedra,
            self.materials,
            None,
            self.fixed_dofs,
            self.budget,
            self.formulation,
            &stiffness_rule,
            mass_rule,
            |nodes, reference, element, point| self.point(nodes, reference, element, point),
            |_, _| None,
            cx,
        )?;
        assembly.minimum_scaled_jacobian = assembly.minimum_scaled_jacobian.min(minimum_probe);
        Ok(assembly)
    }

    pub(super) fn point(
        &self,
        nodes: [[f64; 3]; 8],
        reference: [f64; 3],
        element: usize,
        point: usize,
    ) -> Result<isoparametric::Point<8>, TetElasticError> {
        let p = hex8_point(
            nodes,
            reference,
            self.budget.minimum_scaled_jacobian,
            element,
            point,
        )?;
        Ok(isoparametric::Point {
            shapes: p.shapes,
            gradients: p.gradients,
            determinant: p.determinant_m3,
            scaled_jacobian: p.scaled_jacobian,
        })
    }

    fn validate(&self, require_free_dofs: bool, cx: &Cx<'_>) -> Result<f64, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        if self.nodes_m.is_empty() || self.hexahedra.is_empty() {
            return Err(TetElasticError::EmptyMesh);
        }
        validate_stiffness_order(self.stiffness_points_per_dir)?;
        for (what, requested, maximum) in [
            ("nodes", self.nodes_m.len(), self.budget.maximum_nodes),
            (
                "hexahedra",
                self.hexahedra.len(),
                self.budget.maximum_tetrahedra,
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
        let minimum = self.budget.minimum_scaled_jacobian;
        if !(minimum.is_finite() && minimum > 0.0) {
            return Err(TetElasticError::DegenerateTet {
                element: 0,
                scaled_jacobian: f64::NAN,
                minimum,
            });
        }
        let full = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        let mut fixed = BTreeSet::new();
        for (index, &dof) in self.fixed_dofs.iter().enumerate() {
            if index % 256 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            if dof >= full {
                return Err(TetElasticError::FixedDofOutOfRange {
                    dof,
                    full_dofs: full,
                });
            }
            if !fixed.insert(dof) {
                return Err(TetElasticError::DuplicateFixedDof { dof });
            }
        }
        let free = full - fixed.len();
        if require_free_dofs && free == 0 {
            return Err(TetElasticError::NoFreeDofs);
        }
        if free > self.budget.maximum_free_dofs {
            return Err(TetElasticError::BudgetExceeded {
                what: "free_dofs",
                requested: free,
                maximum: self.budget.maximum_free_dofs,
            });
        }
        for (node, p) in self.nodes_m.iter().enumerate() {
            if node % 256 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            if p.iter().any(|v| !v.is_finite()) {
                return Err(TetElasticError::NonFiniteNode { node });
            }
        }
        self.materials.validate(self.hexahedra.len())?;
        let mut incident = vec![false; self.nodes_m.len()];
        let mut cells = BTreeSet::new();
        let mut minimum_probe = f64::INFINITY;
        for (element, connectivity) in self.hexahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for (local, &node) in connectivity.iter().enumerate() {
                if node >= self.nodes_m.len() {
                    return Err(TetElasticError::NodeIndexOutOfRange { element, node });
                }
                if connectivity[..local].contains(&node) {
                    return Err(TetElasticError::RepeatedTetNode { element });
                }
                incident[node] = true;
            }
            let mut sorted = *connectivity;
            sorted.sort_unstable();
            if !cells.insert(sorted) {
                return Err(TetElasticError::InvalidGeometryUpdate {
                    what: "duplicate hexahedron",
                });
            }
            let nodes = connectivity.map(|i| self.nodes_m[i]);
            // Integration points alone can miss a folded corner. Probe the
            // complete corner set and centre before allocating matrix storage.
            for (point, reference) in HEX27_REFERENCE_NODES[..8]
                .iter()
                .map(|p| p.map(f64::from))
                .chain(std::iter::once([0.0; 3]))
                .enumerate()
            {
                let p = hex8_point(nodes, reference, minimum, element, point)?;
                minimum_probe = minimum_probe.min(p.scaled_jacobian);
            }
        }
        if let Some(node) = incident.iter().position(|v| !*v) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        Ok(minimum_probe)
    }
}

pub(super) fn validate_stiffness_order(order: usize) -> Result<(), TetElasticError> {
    if order < 2 {
        return Err(TetElasticError::UnderIntegrated {
            points_per_dir: order,
        });
    }
    if order > 12 {
        return Err(TetElasticError::QuadratureInadmissible {
            what: "HEX8 stiffness Gauss points per direction must be in 2..=12",
        });
    }
    Ok(())
}

pub(super) fn tensor_rule(order: usize) -> Vec<([f64; 3], f64)> {
    let (points, weights) = gauss_legendre(order);
    let mut rule = Vec::with_capacity(order * order * order);
    for i in 0..order {
        for j in 0..order {
            for k in 0..order {
                rule.push((
                    [points[i], points[j], points[k]],
                    weights[i] * weights[j] * weights[k],
                ));
            }
        }
    }
    rule
}
