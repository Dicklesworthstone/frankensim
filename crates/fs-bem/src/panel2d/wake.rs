//! Retained Hess–Smith solves for nonuniform incident flow and Kelvin shedding.
//!
//! The geometry-only matrix and its LU are shared by every wake stage. The
//! newborn vortex enters both the no-penetration and Kutta rows; eliminating
//! its strength with Kelvin's theorem adds one scalar Schur complement, not
//! an empirical circulation relaxation. Velocities and geometry must use the
//! same frame and units. This solves an inviscid boundary condition, not the
//! unsteady pressure or viscous separation problem.

use super::{
    Airfoil2d, Geometry, assemble, dense_dimension, geometry, source_velocity, zeroed_vec,
};
use crate::BemError;
use fs_la::factor::{FactorError, Lu, lu};

/// The body field satisfying one set of incident velocities and the Kutta row.
#[derive(Debug, Clone)]
pub struct PanelState2d {
    /// Constant source density on each panel, in velocity units.
    pub sources: Vec<f64>,
    /// Shared clockwise-positive vortex density, in velocity units.
    pub sheet_density: f64,
    /// Clockwise-positive circulation: sheet density times section perimeter.
    /// Units are velocity times length; this is not pressure-integrated Cl/2.
    pub bound_circulation: f64,
    /// Total tangential velocity at every panel midpoint, including incidence.
    pub tangential_velocity: Vec<f64>,
    /// Largest absolute total normal velocity at the collocation points.
    pub max_normal_residual: f64,
    /// Signed sum of the two trailing-edge tangential velocities.
    pub kutta_residual: f64,
}

/// One immutable geometry and retained factorization for repeated wake stages.
///
/// Preparation retains the original boundary matrix, its LU, and tangential
/// influence rows in O(n²) storage. Each RHS solve and surface diagnostic is
/// O(n²) arithmetic without reevaluating the logarithmic panel integrals.
#[derive(Debug, Clone)]
pub struct PreparedPanel2d {
    foil: Airfoil2d,
    geometry: Geometry,
    factor: Lu,
    boundary: Vec<f64>,
    tangential: Vec<f64>,
    perimeter: f64,
}

impl PreparedPanel2d {
    /// Prepare the same matrix and principal values used by [`super::solve`].
    /// Geometry and dense-work admission are unchanged from the steady solver.
    pub fn new(foil: &Airfoil2d) -> Result<Self, BemError> {
        let geometry = geometry(foil)?;
        let n = foil.nodes.len();
        let (dim, _) = dense_dimension(n)?;
        let (boundary, _) = assemble(foil, &geometry, 0.0)?;
        let perimeter: f64 = geometry.len.iter().sum();
        if !perimeter.is_finite() || boundary.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "prepared 2D panel assembly",
            });
        }
        let factor = lu(&boundary, dim).map_err(|source| BemError::LinearSolve {
            stage: "prepared 2D panel",
            source,
        })?;
        let mut tangential = zeroed_vec(n * dim, "2D tangential influence rows")?;
        for (i, row) in tangential.chunks_exact_mut(dim).enumerate() {
            let tangent = geometry.tangent[i];
            for j in 0..n {
                if i == j {
                    row[n] += 0.5 * tangent[0] * tangent[0] + 0.5 * tangent[1] * tangent[1];
                } else {
                    let source = source_velocity(foil, &geometry, j, geometry.mid[i]);
                    row[j] = source[0] * tangent[0] + source[1] * tangent[1];
                    row[n] += source[1] * tangent[0] - source[0] * tangent[1];
                }
            }
        }
        if tangential.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "2D tangential influence rows",
            });
        }
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(n)
            .map_err(|_| BemError::AllocationFailed {
                operation: "prepared 2D panel geometry",
            })?;
        nodes.extend_from_slice(&foil.nodes);
        Ok(Self {
            foil: Airfoil2d { nodes },
            geometry,
            factor,
            boundary,
            tangential,
            perimeter,
        })
    }

    /// Collocation points, in the supplied section's frame and length units.
    #[must_use]
    pub fn midpoints(&self) -> &[[f64; 2]] {
        &self.geometry.mid
    }

    /// Clockwise section polygon, available for exterior-point admission.
    #[must_use]
    pub fn nodes(&self) -> &[[f64; 2]] {
        &self.foil.nodes
    }

    /// Solve no penetration and Kutta for the incident velocity at each midpoint.
    /// The incident field includes the freestream and any already-known wake.
    pub fn solve_incident(&self, incident: &[[f64; 2]]) -> Result<PanelState2d, BemError> {
        let coefficients = self.solve_coefficients(incident, "panel incident velocities")?;
        self.state(coefficients, |i| incident[i])
    }

    /// Solve the body and a newborn vortex with Kelvin circulation conservation.
    ///
    /// `incident` excludes the newborn; `unit_newborn_incident` is its velocity
    /// per unit clockwise circulation at the same midpoints. Both fields enter
    /// the normal and Kutta equations. `previous_bound` is the circulation to
    /// conserve between the body and newborn while the older wake is fixed.
    /// The returned scalar is the newborn circulation in velocity-length units.
    ///
    /// Two solves with the retained LU give x₀ and x₁. For perimeter P,
    /// Γ₀ = P x₀[n], Γ₁ = P x₁[n], and γ = (Γ_previous − Γ₀)/(1 + Γ₁).
    /// A numerically singular final Schur pivot is refused without clamping.
    pub fn solve_kelvin(
        &self,
        incident: &[[f64; 2]],
        unit_newborn_incident: &[[f64; 2]],
        previous_bound: f64,
    ) -> Result<(PanelState2d, f64), BemError> {
        if !previous_bound.is_finite() {
            return Err(BemError::InvalidScalar {
                name: "previous bound circulation",
                value: previous_bound,
                requirement: "finite circulation",
            });
        }
        let mut coefficients = self.solve_coefficients(incident, "panel incident velocities")?;
        let response =
            self.solve_coefficients(unit_newborn_incident, "newborn incident velocities")?;
        let n = self.foil.nodes.len();
        let base_circulation = self.perimeter * coefficients[n];
        let circulation_response = self.perimeter * response[n];
        if !base_circulation.is_finite() || !circulation_response.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "2D panel Kelvin response",
            });
        }
        let denominator = 1.0 + circulation_response;
        if denominator.abs() <= 64.0 * f64::EPSILON * (1.0 + circulation_response.abs()) {
            return Err(BemError::LinearSolve {
                stage: "2D panel Kelvin coupling",
                source: FactorError::Singular { index: n + 1 },
            });
        }
        let newborn = (previous_bound - base_circulation) / denominator;
        if !newborn.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "2D panel Kelvin circulation",
            });
        }
        for (value, response) in coefficients.iter_mut().zip(response) {
            *value = newborn.mul_add(response, *value);
        }
        let state = self.state(coefficients, |i| {
            [
                newborn.mul_add(unit_newborn_incident[i][0], incident[i][0]),
                newborn.mul_add(unit_newborn_incident[i][1], incident[i][1]),
            ]
        })?;
        Ok((state, newborn))
    }

    /// Evaluate the induced body field; the caller adds freestream and wake.
    ///
    /// Intended for exterior points away from panel endpoints. On-surface
    /// queries do not replace the principal-value collocation diagnostics.
    pub fn body_velocity(
        &self,
        state: &PanelState2d,
        point: [f64; 2],
    ) -> Result<[f64; 2], BemError> {
        let n = self.foil.nodes.len();
        if state.sources.len() != n {
            return Err(BemError::VectorLength {
                field: "body source densities",
                expected: n,
                actual: state.sources.len(),
            });
        }
        if !state.sheet_density.is_finite() || state.sources.iter().any(|value| !value.is_finite())
        {
            return Err(BemError::NonFiniteResult {
                operation: "2D body field coefficients",
            });
        }
        for value in point {
            if !value.is_finite() {
                return Err(BemError::InvalidScalar {
                    name: "panel field coordinate",
                    value,
                    requirement: "finite coordinate",
                });
            }
        }
        let mut velocity = [0.0, 0.0];
        for (j, &density) in state.sources.iter().enumerate() {
            let source = source_velocity(&self.foil, &self.geometry, j, point);
            velocity[0] += source[0] * density + source[1] * state.sheet_density;
            velocity[1] += source[1] * density - source[0] * state.sheet_density;
        }
        if velocity.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "2D body field evaluation",
            });
        }
        Ok(velocity)
    }

    fn solve_coefficients(
        &self,
        incident: &[[f64; 2]],
        field: &'static str,
    ) -> Result<Vec<f64>, BemError> {
        let n = self.foil.nodes.len();
        if incident.len() != n {
            return Err(BemError::VectorLength {
                field,
                expected: n,
                actual: incident.len(),
            });
        }
        if let Some(&value) = incident.iter().flatten().find(|value| !value.is_finite()) {
            return Err(BemError::InvalidScalar {
                name: field,
                value,
                requirement: "finite velocity components",
            });
        }
        let mut rhs = zeroed_vec(n + 1, "2D incident velocity RHS")?;
        for (i, velocity) in incident.iter().enumerate() {
            let normal = self.geometry.normal[i];
            rhs[i] = -(velocity[0] * normal[0] + velocity[1] * normal[1]);
        }
        for &i in &[0usize, n - 1] {
            let tangent = self.geometry.tangent[i];
            rhs[n] -= incident[i][0] * tangent[0] + incident[i][1] * tangent[1];
        }
        self.factor.solve(&mut rhs);
        if rhs.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "2D incident velocity solution",
            });
        }
        Ok(rhs)
    }

    fn state(
        &self,
        mut coefficients: Vec<f64>,
        incident: impl Fn(usize) -> [f64; 2],
    ) -> Result<PanelState2d, BemError> {
        let n = self.foil.nodes.len();
        let dim = n + 1;
        let sheet_density = coefficients[n];
        let bound_circulation = self.perimeter * sheet_density;
        if !bound_circulation.is_finite() || coefficients.iter().any(|value| !value.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "2D coupled panel solution",
            });
        }
        let mut tangential_velocity = zeroed_vec(n, "2D coupled tangential velocities")?;
        let mut max_normal_residual = 0.0f64;
        for (i, tangential_velocity) in tangential_velocity.iter_mut().enumerate() {
            let velocity = incident(i);
            let normal = self.geometry.normal[i];
            let tangent = self.geometry.tangent[i];
            let mut normal_velocity = velocity[0] * normal[0] + velocity[1] * normal[1];
            *tangential_velocity = velocity[0] * tangent[0] + velocity[1] * tangent[1];
            for (j, &coefficient) in coefficients.iter().enumerate() {
                normal_velocity = self.boundary[i * dim + j].mul_add(coefficient, normal_velocity);
                *tangential_velocity =
                    self.tangential[i * dim + j].mul_add(coefficient, *tangential_velocity);
            }
            if !normal_velocity.is_finite() || !tangential_velocity.is_finite() {
                return Err(BemError::NonFiniteResult {
                    operation: "2D coupled surface velocity",
                });
            }
            max_normal_residual = max_normal_residual.max(normal_velocity.abs());
        }
        let kutta_residual = tangential_velocity[0] + tangential_velocity[n - 1];
        if !kutta_residual.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "2D coupled Kutta residual",
            });
        }
        coefficients.truncate(n);
        Ok(PanelState2d {
            sources: coefficients,
            sheet_density,
            bound_circulation,
            tangential_velocity,
            max_normal_residual,
            kutta_residual,
        })
    }
}
