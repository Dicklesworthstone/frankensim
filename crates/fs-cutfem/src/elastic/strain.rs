//! Three-dimensional physical observations of the original plane Q1 field.
//!
//! The plane-stress reduction retains its relaxed thickness strain. All tensor
//! components use fs-material's `[xx, yy, zz, xy, yz, zx]` tensor-shear convention;
//! `gamma_xy = 2 epsilon_xy`. Energy is physical bulk energy, excluding Nitsche,
//! ghost stabilization and the work of prescribed boundary displacements.

use super::{CellKey, CutElasticityOperator, CutElasticitySolution, CutFemError, ElasticPlane};
use crate::Quadtree;
use crate::fem::q1;
use crate::grid::{cell_corner_nodes, cell_rect};
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt, tensor};
use std::ops::ControlFlow;

/// Local response of the immutable material and plane reduction used to solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElasticStrainStress {
    /// Dimensionless total symmetric strain, including relaxed epsilon_zz.
    pub strain: Voigt,
    /// Actual constitutive strain, `strain - free_strain`, with tensor shear.
    pub elastic_strain: Voigt,
    /// Original prescribed free strain at the observation point, or zero.
    pub free_strain: Voigt,
    /// Cauchy stress [Pa] in the same tensor order, including sigma_zz.
    pub stress: Voigt,
    /// Physical storage per reference volume [J/m^3], `sigma:elastic_strain/2`.
    pub energy_density: f64,
    /// Rotation-invariant total-strain norm, including both shear entries.
    pub strain_magnitude: f64,
    /// Rotation-invariant norm of the actual elastic strain.
    pub elastic_strain_magnitude: f64,
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn strain(material: &IsotropicElastic, plane: ElasticPlane, gradient: [[f64; 2]; 2]) -> Voigt {
    let xx = gradient[0][0];
    let yy = gradient[1][1];
    let zz = match plane {
        ElasticPlane::Strain => 0.0,
        ElasticPlane::Stress => -material.poisson / (1.0 - material.poisson) * (xx + yy),
    };
    [
        xx,
        yy,
        zz,
        f64::midpoint(gradient[0][1], gradient[1][0]),
        0.0,
        0.0,
    ]
}

fn strain_pair(
    material: &IsotropicElastic,
    plane: ElasticPlane,
    gradient: [[f64; 2]; 2],
    free: Option<Voigt>,
) -> (Voigt, Voigt) {
    let mut total = strain(material, plane, gradient);
    let Some(free) = free else {
        return (total, total);
    };
    let mut elastic = core::array::from_fn(|i| total[i] - free[i]);
    if plane == ElasticPlane::Stress {
        // The transverse equilibrium constrains ELASTIC stress/strain. In
        // particular free isotropic expansion relaxes to epsilon_zz=epsilon_free,
        // rather than acquiring the cold law's artificial Poisson contraction.
        elastic[2] = -material.poisson / (1.0 - material.poisson) * (elastic[0] + elastic[1]);
        elastic[4] = 0.0;
        elastic[5] = 0.0;
        total[2] = free[2] + elastic[2];
        total[4] = free[4];
        total[5] = free[5];
    }
    (total, elastic)
}

fn admit(material: &IsotropicElastic, strain: &Voigt) -> Result<f64, CutFemError> {
    let magnitude = tensor::frobenius_norm(strain);
    if !magnitude.is_finite() || magnitude > material.strain_limit {
        return Err(invalid(format!(
            "three-dimensional tensor strain magnitude {magnitude} exceeds the original material limit {}",
            material.strain_limit
        )));
    }
    Ok(magnitude)
}

fn admit_total(strain: &Voigt, maximum: f64) -> Result<f64, CutFemError> {
    let magnitude = tensor::frobenius_norm(strain);
    if !magnitude.is_finite() || magnitude > maximum {
        return Err(invalid(format!(
            "total tensor strain magnitude {magnitude} exceeds the declared kinematic limit {maximum}"
        )));
    }
    Ok(magnitude)
}

fn response(
    material: &IsotropicElastic,
    plane: ElasticPlane,
    gradient: [[f64; 2]; 2],
    free: Option<Voigt>,
    maximum_kinematic_strain: f64,
) -> Result<ElasticStrainStress, CutFemError> {
    let (strain, elastic_strain) = strain_pair(material, plane, gradient, free);
    let strain_magnitude = if free.is_some() {
        admit_total(&strain, maximum_kinematic_strain)?
    } else {
        admit(material, &strain)?
    };
    let elastic_strain_magnitude = if free.is_some() {
        admit(material, &elastic_strain)?
    } else {
        strain_magnitude
    };
    // Use the SAME reduced coefficients as assembly. Evaluating 3-D lambda on
    // a nearly incompressible plane-stress state would lose stress by cancelling
    // its large hydrostatic terms, despite the bounded physical plane response.
    let (lambda, mu) = super::reduced_lame(material, plane);
    let in_plane = super::constitutive_mul(
        lambda,
        mu,
        [
            elastic_strain[0],
            elastic_strain[1],
            2.0 * elastic_strain[3],
        ],
    );
    let zz = match plane {
        ElasticPlane::Strain => lambda * (strain[0] + strain[1]),
        ElasticPlane::Stress => 0.0,
    };
    let stress = if free.is_some() && plane == ElasticPlane::Strain {
        // The total transverse strain vanishes, but its elastic part need not.
        // Retain all native transverse stresses and their stored energy.
        material.stress(&elastic_strain, &())
    } else {
        [in_plane[0], in_plane[1], zz, in_plane[2], 0.0, 0.0]
    };
    let energy_density = 0.5 * tensor::contract(&elastic_strain, &stress);
    if stress.iter().any(|v| !v.is_finite()) || !energy_density.is_finite() {
        return Err(invalid("unrepresentable physical elastic stress or energy"));
    }
    Ok(ElasticStrainStress {
        strain,
        elastic_strain,
        free_strain: free.unwrap_or([0.0; 6]),
        stress,
        energy_density,
        strain_magnitude,
        elastic_strain_magnitude,
    })
}

impl CutElasticitySolution {
    /// Recover the original three-dimensional material response at one point.
    /// Uses the same displacement gradient as [`Self::value_gradient`].
    /// No new material or plane assumption can be supplied after the solve.
    ///
    /// An active cut cell includes a fictitious extension outside the body.
    /// This method evaluates that polynomial; it does not classify the point
    /// against the SDF. The caller must select a point in the physical domain
    /// when interpreting this as a specimen observation.
    ///
    /// # Errors
    /// All field-evaluation refusals and the original material strain ceiling.
    /// Algebraic solver convergence alone is not material-domain admission.
    pub fn strain_stress(
        &self,
        grid: &Quadtree,
        cell: CellKey,
        point: [f64; 2],
    ) -> Result<ElasticStrainStress, CutFemError> {
        let (_, gradient) = self.value_gradient(grid, cell, point)?;
        response(
            &self.material,
            self.plane,
            gradient,
            None,
            self.material.strain_limit,
        )
    }
}

impl CutElasticityOperator {
    fn observation_values(
        &self,
        coefficients: &[f64],
        cell: CellKey,
    ) -> Result<[[f64; 2]; 4], CutFemError> {
        if coefficients.len() != self.dof_count() {
            return Err(invalid("strain observation requires complete coefficients"));
        }
        if self.active.binary_search(&cell).is_err() {
            return Err(invalid(format!(
                "cannot observe inactive elasticity cell {cell:?}"
            )));
        }
        let mut values = [[0.0; 2]; 4];
        for (value, node) in values
            .iter_mut()
            .zip(cell_corner_nodes(cell, self.lattice_level))
        {
            if let Some(expansions) = &self.node_expansions {
                for &(id, weight) in &expansions[&node] {
                    for c in 0..2 {
                        value[c] += weight * coefficients[2 * id + c];
                    }
                }
            } else {
                let id = self.node_ids[&node];
                *value = [coefficients[2 * id], coefficients[2 * id + 1]];
            }
            if value.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "nonfinite physical nodal displacement in strain observation",
                ));
            }
        }
        Ok(values)
    }

    /// Observe any complete displacement vector, including a controlled solve's
    /// coefficients, on the original material and hanging-node trial space.
    /// This bounded point query allocates nothing on success and reconstructs
    /// only the four required physical nodes, rather than a whole nodal field.
    /// It observes the active-cell extension with the same physical-domain
    /// caveat as [`CutElasticitySolution::strain_stress`].
    ///
    /// # Errors
    /// Inactive cell, wrong lattice/vector, nonfinite data, point outside the
    /// cell, arithmetic overflow or the original full tensor strain limit.
    pub fn strain_stress(
        &self,
        grid: &Quadtree,
        coefficients: &[f64],
        cell: CellKey,
        point: [f64; 2],
    ) -> Result<ElasticStrainStress, CutFemError> {
        if grid.max_level() != self.lattice_level {
            return Err(invalid(
                "strain observation requires the original node lattice",
            ));
        }
        let values = self.observation_values(coefficients, cell)?;
        let gradient = gradient(cell, point, &values)?;
        response(
            &self.material,
            self.plane,
            gradient,
            self.eigenstrain_at(cell, point)?,
            self.maximum_kinematic_strain(),
        )
    }

    /// Admit the full Q1 strain field, rather than only sampled quadrature sites.
    /// Returns the maximum total or elastic tensor norm over every active cell.
    /// The original gradient is affine; an attached Q1 free-strain field makes
    /// the relaxed total and elastic strains bilinear. Both are affine in each
    /// coordinate separately, so their convex norms attain a maximum at a vertex.
    /// Four observations cover each complete polynomial in real arithmetic.
    /// With attached eigenstrain, total strain uses its explicit kinematic
    /// ceiling while elastic strain retains the original material calibration.
    ///
    /// For cut cells this bounds the fictitious extension too; refusal can be
    /// conservative for the physical subdomain. This is a floating-point model
    /// bound, not an interval enclosure, finite-strain or continuum-error proof.
    /// No FEM basis vector or unit adjoint is automatically treated as a physical
    /// displacement. The supplied vector must be the actual state being assessed.
    ///
    /// Cancellation is polled per cell and before publication; no partial bound
    /// is returned. Work is O(active cells * log(nodes)), with no heap allocation.
    ///
    /// # Errors
    /// The same field and original material-strain refusals as point recovery.
    pub fn admit_strain_domain<B>(
        &self,
        grid: &Quadtree,
        coefficients: &[f64],
        control: impl FnMut() -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, f64>, CutFemError> {
        if grid.max_level() != self.lattice_level {
            return Err(invalid(
                "strain observation requires the original node lattice",
            ));
        }
        self.admit_retained_strain_domain(coefficients, control)
    }

    // The solve owns all active cells, node maps and material data already.
    // Reuse the physical-field gate without storing or reconstructing a tree.
    pub(crate) fn admit_retained_strain_domain<B>(
        &self,
        coefficients: &[f64],
        mut control: impl FnMut() -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, f64>, CutFemError> {
        let mut maximum = 0.0_f64;
        for &cell in &self.active {
            if let ControlFlow::Break(reason) = control() {
                return Ok(ControlFlow::Break(reason));
            }
            let values = self.observation_values(coefficients, cell)?;
            let (lo, hi) = cell_rect(cell);
            for point in [lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]] {
                let gradient = gradient(cell, point, &values)?;
                let free = self.eigenstrain_at(cell, point)?;
                let (total, elastic) = strain_pair(&self.material, self.plane, gradient, free);
                if free.is_some() {
                    maximum = maximum.max(admit_total(&total, self.maximum_kinematic_strain())?);
                    maximum = maximum.max(admit(&self.material, &elastic)?);
                } else {
                    maximum = maximum.max(admit(&self.material, &total)?);
                }
            }
        }
        if let ControlFlow::Break(reason) = control() {
            return Ok(ControlFlow::Break(reason));
        }
        Ok(ControlFlow::Continue(maximum))
    }
}

fn gradient(
    cell: CellKey,
    point: [f64; 2],
    values: &[[f64; 2]; 4],
) -> Result<[[f64; 2]; 2], CutFemError> {
    let (lo, hi) = cell_rect(cell);
    if (0..2).any(|a| !point[a].is_finite() || point[a] < lo[a] || point[a] > hi[a]) {
        return Err(invalid(
            "strain observation point must be finite and inside its cell",
        ));
    }
    let (_, gradients) = q1(lo, hi, point);
    let gradient = gradient_from_values(&gradients, values);
    if gradient.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid(
            "unrepresentable displacement gradient in strain observation",
        ));
    }
    Ok(gradient)
}

// One kinematic owner for solved-field and operator observations. A rigid
// translation has zero strain: remove it before differentiating so the final
// sum does not have to cancel large, individually rounded products.
pub(super) fn gradient_from_values(
    gradients: &[[f64; 2]; 4],
    values: &[[f64; 2]; 4],
) -> [[f64; 2]; 2] {
    let mut gradient = [[0.0; 2]; 2];
    for corner in 0..4 {
        for component in 0..2 {
            for axis in 0..2 {
                gradient[component][axis] +=
                    gradients[corner][axis] * (values[corner][component] - values[0][component]);
            }
        }
    }
    gradient
}
