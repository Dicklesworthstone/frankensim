//! Prescribed initial strain in the original plane CutFEM equilibrium.
//!
//! The physical stress is `C (epsilon(u) - epsilon_free)`. Its weak load is
//! `integral B^T C epsilon_free dA` and, on the Nitsche displacement boundary,
//! `-integral N^T (C epsilon_free) n ds`. The latter sign follows by substituting
//! the physical stress in the consistency traction; omitting it gives spurious
//! forces even under completely restrained uniform expansion.
//!
//! Fields are explicitly cellwise Q1, with four supplied corner tensors. They
//! may be discontinuous across material cells. A callback supplies coefficients,
//! not an arbitrary continuous function whose interpolation is claimed exact.
//! Thermal callers supply integrated expansion from their actual stress-free
//! reference, including the out-of-plane component. No temperature law, source
//! extrapolation, heat balance or updated geometry is inferred here.

use super::{CellKey, CutElasticityOperator, CutFemError, ElasticPlane};
use crate::fem::q1;
use crate::grid::{cell_corner_nodes, cell_rect};
use crate::quad::tensor_gauss;
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt, tensor};
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use std::sync::Arc;

/// Dimensionless three-dimensional free strain on one Q1 background cell.
/// Tensors use `[xx, yy, zz, xy, yz, zx]` with physical tensor shear.
/// Corners are lower-left, lower-right, upper-right, upper-left.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Q1Eigenstrain {
    corners: [Voigt; 4],
}

impl Q1Eigenstrain {
    /// Admit finite tensor coefficients for an explicitly bilinear field.
    /// Constitutive and total-strain limits are checked on the solved field;
    /// an eigenstrain itself is neither an elastic strain nor a displacement.
    pub fn try_new(corners: [Voigt; 4]) -> Result<Self, CutFemError> {
        if corners
            .iter()
            .any(|e| e.iter().any(|v| !v.is_finite()) || !tensor::frobenius_norm(e).is_finite())
        {
            return Err(invalid(
                "Q1 eigenstrain requires finite tensor coefficients",
            ));
        }
        Ok(Self { corners })
    }

    /// One constant, fully three-dimensional initial strain throughout a cell.
    pub fn try_uniform(strain: Voigt) -> Result<Self, CutFemError> {
        Self::try_new([strain; 4])
    }

    /// Isotropic integrated linear expansion at each original corner.
    /// Values are dimensionless `integral alpha(T) dT`, not kelvin or an
    /// endpoint coefficient. This retains expansion in all three directions.
    pub fn try_isotropic(expansion: [f64; 4]) -> Result<Self, CutFemError> {
        Self::try_new(expansion.map(|e| [e, e, e, 0.0, 0.0, 0.0]))
    }

    /// The original coefficients, including their out-of-plane components.
    #[must_use]
    pub const fn corner_strains(&self) -> &[Voigt; 4] {
        &self.corners
    }

    pub(super) fn at(&self, shape: [f64; 4]) -> Result<Voigt, CutFemError> {
        // A constant source remains exactly constant even when the four shape
        // weights do not sum to one bit-for-bit. Translation of this field
        // similarly does not amplify cancellation between individual products.
        let strain = core::array::from_fn(|component| {
            (1..4).fold(self.corners[0][component], |sum, corner| {
                shape[corner].mul_add(
                    self.corners[corner][component] - self.corners[0][component],
                    sum,
                )
            })
        });
        if strain.iter().any(|v| !v.is_finite()) {
            return Err(invalid("Q1 eigenstrain interpolation is unrepresentable"));
        }
        Ok(strain)
    }
}

#[derive(Debug)]
pub(super) struct EigenstrainField {
    maximum_kinematic_strain: f64,
    cells: BTreeMap<CellKey, Q1Eigenstrain>,
}

pub(super) type RetainedEigenstrain = Arc<EigenstrainField>;

fn invalid(what: &str) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn initial_stress(
    material: &IsotropicElastic,
    plane: ElasticPlane,
    free: Voigt,
) -> Result<[f64; 3], CutFemError> {
    let stress = match plane {
        ElasticPlane::Strain => {
            // epsilon_zz = 0 constrains TOTAL strain. The initial zz strain
            // still contributes through the three-dimensional bulk modulus.
            let full = material.stress(&free, &());
            [full[0], full[1], full[3]]
        }
        ElasticPlane::Stress => {
            // Eliminate all transverse elastic stresses analytically before
            // evaluating the volumetric law, as in the original plane owner.
            let (lambda, mu) = super::reduced_lame(material, plane);
            super::constitutive_mul(lambda, mu, [free[0], free[1], 2.0 * free[3]])
        }
    };
    if stress.iter().any(|v| !v.is_finite()) {
        return Err(invalid("initial elastic stress is unrepresentable"));
    }
    Ok(stress)
}

impl CutElasticityOperator {
    /// Add prescribed initial strain to the original assembled equilibrium.
    ///
    /// The callback is called once per active cell with that cell's physical
    /// corners. It returns the four coefficients of an explicitly Q1 free-strain
    /// field. Supply original source-integrated thermal, hygroscopic or residual
    /// strain; no material-name inference or temperature approximation occurs.
    /// All loads use unit reference thickness, like the original plane operator.
    /// `maximum_kinematic_strain` is an explicit Frobenius-norm ceiling in (0,1)
    /// for TOTAL strain. The original material's calibration limit continues to
    /// bound ELASTIC strain. Free expansion does not consume the elastic domain,
    /// and that domain is never enlarged to make a thermal state pass. The
    /// kinematic ceiling is a small-strain model budget, not a finite-rotation
    /// or deformation-injectivity proof.
    ///
    /// This consumes the admitted operator and reuses its stiffness, cut rules,
    /// supports, hanging-node transform and original plane reduction. The field
    /// and complete RHS are retained by the returned candidate; later solves and
    /// observations never rerun the callback. A second attachment refuses so
    /// source refresh cannot silently apply the old expansion twice. Keep an
    /// unmodified operator clone when reusing one geometry for thermal updates.
    ///
    /// Use [`Self::solve_controlled`] on the returned operator. That same solver
    /// now admits both total and elastic tensor strain at every cell corner;
    /// [`Self::strain_stress`] reports the actual physical stress and storage.
    /// There is no change to ghost stabilization or stiffness VJP semantics.
    ///
    /// Cancellation is polled per cell, every 64 integration sites, and before
    /// publication. `Break` returns the caller's reason without a partial field.
    /// Work is linear in retained quadrature and sparse nodal expansions. No cut
    /// geometry, global stiffness or external source is reconstructed.
    pub fn with_eigenstrain<B>(
        mut self,
        maximum_kinematic_strain: f64,
        mut field: impl FnMut(CellKey, [[f64; 2]; 4]) -> Result<Q1Eigenstrain, CutFemError>,
        mut control: impl FnMut() -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, Self>, CutFemError> {
        if let ControlFlow::Break(reason) = control() {
            return Ok(ControlFlow::Break(reason));
        }
        if self.eigenstrain.is_some() {
            return Err(invalid(
                "an initial-strain field is already attached to this operator",
            ));
        }
        if !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(invalid(
                "the total-strain kinematic ceiling must lie in (0,1)",
            ));
        }
        let mut retained = BTreeMap::new();
        for &cell in &self.active {
            if let ControlFlow::Break(reason) = control() {
                return Ok(ControlFlow::Break(reason));
            }
            let (lo, hi) = cell_rect(cell);
            let corners = [lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]];
            let eigenstrain = field(cell, corners)?;
            let mut local = [0.0; 8];
            let full_rule;
            let bulk = if let Some(rule) = self.rules.get(&cell) {
                rule.bulk.as_slice()
            } else {
                full_rule = {
                    let mut points = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut points);
                    points
                };
                full_rule.as_slice()
            };
            for (index, &(point, weight)) in bulk.iter().enumerate() {
                if index % 64 == 0 {
                    if let ControlFlow::Break(reason) = control() {
                        return Ok(ControlFlow::Break(reason));
                    }
                }
                let (shape, gradients) = q1(lo, hi, point);
                let stress = initial_stress(&self.material, self.plane, eigenstrain.at(shape)?)?;
                for node in 0..4 {
                    for component in 0..2 {
                        local[2 * node + component] += weight
                            * super::dot3(super::strain_row(gradients[node], component), stress);
                    }
                }
            }
            if self.embedded_dirichlet {
                if let Some(rule) = self.rules.get(&cell) {
                    for (index, &(point, weight, normal)) in rule.iface.iter().enumerate() {
                        if index % 64 == 0 {
                            if let ControlFlow::Break(reason) = control() {
                                return Ok(ControlFlow::Break(reason));
                            }
                        }
                        let (shape, _) = q1(lo, hi, point);
                        let [xx, yy, xy] =
                            initial_stress(&self.material, self.plane, eigenstrain.at(shape)?)?;
                        let traction = [
                            xx * normal[0] + xy * normal[1],
                            xy * normal[0] + yy * normal[1],
                        ];
                        for node in 0..4 {
                            for component in 0..2 {
                                local[2 * node + component] -=
                                    weight * shape[node] * traction[component];
                            }
                        }
                    }
                }
            }
            if local.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "initial-strain load integration is unrepresentable",
                ));
            }
            for (local_node, node) in cell_corner_nodes(cell, self.lattice_level)
                .iter()
                .enumerate()
            {
                if let Some(expansions) = &self.node_expansions {
                    for &(id, weight) in &expansions[node] {
                        for component in 0..2 {
                            let dof = 2 * id + component;
                            if !self.clamped[dof] {
                                let value = weight * local[2 * local_node + component];
                                if value != 0.0 {
                                    self.rhs[dof] += value;
                                    if !self.rhs[dof].is_finite() {
                                        return Err(invalid(
                                            "initial-strain load scatter is unrepresentable",
                                        ));
                                    }
                                }
                            }
                        }
                    }
                } else {
                    let id = self.node_ids[node];
                    for component in 0..2 {
                        let dof = 2 * id + component;
                        let value = local[2 * local_node + component];
                        if !self.clamped[dof] && value != 0.0 {
                            self.rhs[dof] += value;
                            if !self.rhs[dof].is_finite() {
                                return Err(invalid(
                                    "initial-strain load scatter is unrepresentable",
                                ));
                            }
                        }
                    }
                }
            }
            retained.insert(cell, eigenstrain);
        }
        if let ControlFlow::Break(reason) = control() {
            return Ok(ControlFlow::Break(reason));
        }
        self.eigenstrain = Some(Arc::new(EigenstrainField {
            maximum_kinematic_strain,
            cells: retained,
        }));
        Ok(ControlFlow::Continue(self))
    }

    /// Original Q1 free-strain coefficients for an active cell, when attached.
    #[must_use]
    pub fn eigenstrain_cell(&self, cell: CellKey) -> Option<&Q1Eigenstrain> {
        self.eigenstrain.as_ref()?.cells.get(&cell)
    }

    /// Total-strain model ceiling, independently declared when initial strain
    /// is attached. Without it, the original material-strain ceiling applies.
    #[must_use]
    pub fn maximum_kinematic_strain(&self) -> f64 {
        self.eigenstrain
            .as_ref()
            .map_or(self.material.strain_limit, |field| {
                field.maximum_kinematic_strain
            })
    }

    pub(super) fn eigenstrain_at(
        &self,
        cell: CellKey,
        point: [f64; 2],
    ) -> Result<Option<Voigt>, CutFemError> {
        let Some(fields) = &self.eigenstrain else {
            return Ok(None);
        };
        let field = fields
            .cells
            .get(&cell)
            .ok_or_else(|| invalid("missing retained cell eigenstrain"))?;
        let (lo, hi) = cell_rect(cell);
        let (shape, _) = q1(lo, hi, point);
        field.at(shape).map(Some)
    }
}
