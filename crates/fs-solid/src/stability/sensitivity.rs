//! Total buckling sensitivities for element-wise elastic sizing.
//!
//! For fixed geometry, Poisson ratio, dead loads and prescribed displacements,
//! `K(s) u(s) = f` and `(K(s) + lambda G(s, u(s))) phi = 0`. Differentiating
//! only K freezes a stress state that generally redistributes with the design.
//! This owner includes both the direct material derivative of G and its
//! displacement derivative. One static adjoint per accepted mode supplies all
//! element derivatives; no perturbed forward solves or dense third-order
//! derivative tensors are formed. P1/Q1 and Standard/B-bar use the original
//! equilibrium quadrature, projection and homogeneous incremental supports.

use super::{
    BucklingResult, buckling_loads_cancellable, check_cancelled, expand_mode,
    reduced_pencil_scaled, reference_strain,
};
use crate::SolidError;
use crate::linear::{Formulation, Jacobi, LinearProblem, lame};
use crate::mesh2::{quad_points, shapes_at};
use fs_solver::checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg};
use fs_solver::op::CsrOp;

/// Accepted modes and their total derivatives with respect to the original
/// element Young's-modulus multipliers. Fields cannot be replaced by unrelated
/// modes, gradients or scale vectors after analysis.
#[derive(Debug, Clone)]
pub struct ElementScaleSensitivity {
    buckling: BucklingResult,
    scales: Vec<f64>,
    gradients: Vec<Vec<f64>>,
    adjoint_relative_residuals: Vec<f64>,
}

impl ElementScaleSensitivity {
    /// Physical critical loads, homogeneous-support modes and pencil residuals.
    #[must_use]
    pub fn buckling(&self) -> &BucklingResult {
        &self.buckling
    }

    /// The dimensionless multipliers used in this complete forward/adjoint run.
    #[must_use]
    pub fn element_scales(&self) -> &[f64] {
        &self.scales
    }

    /// `gradients[mode][element] = d lambda_mode / d scale_element`.
    /// Only spectrally separated modes are admitted by this owner.
    #[must_use]
    pub fn gradients(&self) -> &[Vec<f64>] {
        &self.gradients
    }

    /// Recomputed true relative residual for each static adjoint, at most 1e-8.
    #[must_use]
    pub fn adjoint_relative_residuals(&self) -> &[f64] {
        &self.adjoint_relative_residuals
    }
}

/// Solve the scaled elastic equilibrium and compute total element sensitivities
/// of its lowest positive proportional buckling loads.
///
/// Loads and imposed displacement values are independent of the design. For a
/// density or other nonlinear parameterization, apply its scalar chain rule to
/// the returned derivatives. This is infinitesimal initial-stress buckling; it
/// does not differentiate load-following geometry or a nonlinear equilibrium.
///
/// One additional eigenpair is requested to check separation at the retained
/// spectrum boundary. Relative gaps <= 1e-6 refuse: individual branches do not
/// have a unique gradient at a repeated eigenvalue. Use a complete-cluster
/// spectral formulation for that case. The existing 4096-free-DOF factor cap
/// and bounded small-pencil fallback also apply here.
///
/// # Errors
/// Invalid scales, unresolved clusters, unrepresentable derivatives, failed
/// equilibrium, eigenpair or adjoint residual gates return `SolidError`.
pub fn element_scale_sensitivities(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    modes: usize,
    iterations: usize,
) -> Result<ElementScaleSensitivity, SolidError> {
    element_scale_sensitivities_cancellable(problem, scales, modes, iterations, &|| false)
}

/// The same analysis with cancellation at element boundaries and within the
/// native equilibrium, eigenvalue and static-adjoint iterations. Cancellation
/// publishes no partial sensitivity result.
///
/// # Errors
/// As [`element_scale_sensitivities`], plus `SolidError::Cancelled`.
pub fn element_scale_sensitivities_cancellable(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    modes: usize,
    iterations: usize,
    cancelled: &impl Fn() -> bool,
) -> Result<ElementScaleSensitivity, SolidError> {
    check_cancelled(cancelled)?;
    if modes == 0 || iterations == 0 {
        return Err(invalid(
            "buckling sensitivities require positive mode and iteration allowances",
        ));
    }
    let (k, kg, map, reference) = reduced_pencil_scaled(problem, Some(scales), cancelled)?;
    let mut buckling = buckling_loads_cancellable(
        &k,
        &kg,
        &map,
        modes.saturating_add(1),
        iterations,
        cancelled,
    )?;
    let retained = modes.min(buckling.loads.len());
    for i in 0..retained {
        if let Some(&next) = buckling.loads.get(i + 1) {
            if (next - buckling.loads[i]) / next <= 1e-6 {
                return Err(invalid(
                    "buckling sensitivity requires separated eigenvalues, including the retained spectrum boundary",
                ));
            }
        }
    }
    buckling.loads.truncate(retained);
    buckling.modes.truncate(retained);
    buckling.relative_residuals.truncate(retained);
    let preconditioner = Jacobi::new(&k);
    let op = CsrOp::symmetric(k);
    let mut gradients = Vec::with_capacity(retained);
    let mut residuals = Vec::with_capacity(retained);
    for (mode, &load) in buckling.modes.iter().zip(&buckling.loads) {
        check_cancelled(cancelled)?;
        let full_mode = expand_mode(mode, &map);
        let mut geometric_mode = vec![0.0; mode.len()];
        kg.spmv(mode, &mut geometric_mode);
        let denominator = -mode
            .iter()
            .zip(&geometric_mode)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(invalid(
                "buckling sensitivity requires positive finite geometric modal work",
            ));
        }
        let (rhs, direct, reference_forces) = contractions(
            problem, scales, &map, &reference, &full_mode, load, cancelled,
        )?;
        let adjoint = checked_cg(
            &op,
            &preconditioner,
            &rhs,
            CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 40_000,
                max_corrections: 2,
            },
            |_| check_cancelled(cancelled),
        )
        .map_err(|error| match error {
            CheckedCgError::InvalidInput(what) => invalid(what),
            CheckedCgError::Interrupted(error) => error,
        })?;
        let residual = adjoint.report.rel_residual;
        if !residual.is_finite() || residual > 1e-8 {
            return Err(SolidError::SolveFailed {
                iters: adjoint.report.iters,
                rel_residual: residual,
            });
        }
        let mut gradient = Vec::with_capacity(scales.len());
        for (element, conn) in problem.mesh.elems.iter().enumerate() {
            check_cancelled(cancelled)?;
            // The full K'_e u includes prescribed columns. The adjoint and
            // mode are homogeneous only on the incremental fixed DOFs.
            let mut redistribution = 0.0;
            for (a, &node) in conn.iter().enumerate() {
                for c in 0..2 {
                    if let Some(free) = map[2 * node + c] {
                        redistribution += adjoint.x[free] * reference_forces[element][2 * a + c];
                    }
                }
            }
            let value = (-load).mul_add(redistribution, direct[element]) / denominator;
            if !value.is_finite() {
                return Err(invalid("total buckling derivative is not representable"));
            }
            gradient.push(value);
        }
        gradients.push(gradient);
        residuals.push(residual);
    }
    check_cancelled(cancelled)?;
    Ok(ElementScaleSensitivity {
        buckling,
        scales: scales.to_vec(),
        gradients,
        adjoint_relative_residuals: residuals,
    })
}

/// q = d(phi^T G phi)/du, direct = phi^T(K'_e + lambda G'_e)phi,
/// and unscaled element internal forces K'_e u. Storage is O(elements + DOFs);
/// quadrature work is linear in element arity rather than a dense Hessian.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn contractions(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    map: &[Option<usize>],
    reference: &[[f64; 2]],
    mode: &[[f64; 2]],
    load: f64,
    cancelled: &impl Fn() -> bool,
) -> Result<(Vec<f64>, Vec<f64>, Vec<[f64; 8]>), SolidError> {
    let (lambda, mu) = lame(problem.youngs, problem.poisson, problem.plane);
    let mut rhs = vec![0.0; map.iter().flatten().count()];
    let mut direct = vec![0.0; scales.len()];
    let mut forces = vec![[0.0; 8]; scales.len()];
    for (element, conn) in problem.mesh.elems.iter().enumerate() {
        check_cancelled(cancelled)?;
        let points = quad_points(conn.len());
        let mut mean_grad = [[0.0; 2]; 4];
        let mut mean_reference = 0.0;
        let mut mean_mode = 0.0;
        if problem.formulation == Formulation::BBar {
            let mut volume = 0.0;
            for &(xi, eta, weight) in &points {
                let (_, grads, det) = shapes_at(&problem.mesh.nodes, conn, xi, eta);
                let w = weight * det;
                volume += w;
                let u = reference_strain(conn, &grads, reference);
                let p = reference_strain(conn, &grads, mode);
                mean_reference += w * (u[0] + u[1]);
                mean_mode += w * (p[0] + p[1]);
                for (a, g) in grads.iter().enumerate() {
                    for c in 0..2 {
                        mean_grad[a][c] += w * g[c];
                    }
                }
            }
            mean_reference /= volume;
            mean_mode /= volume;
            for g in &mut mean_grad[..conn.len()] {
                for component in g {
                    *component /= volume;
                }
            }
        }
        for &(xi, eta, weight) in &points {
            let (_, grads, det) = shapes_at(&problem.mesh.nodes, conn, xi, eta);
            let w = weight * det;
            let mut eu = reference_strain(conn, &grads, reference);
            let mut ep = reference_strain(conn, &grads, mode);
            if problem.formulation == Formulation::BBar {
                for (strain, mean) in [(&mut eu, mean_reference), (&mut ep, mean_mode)] {
                    let correction = 0.5 * (mean - strain[0] - strain[1]);
                    strain[0] += correction;
                    strain[1] += correction;
                }
            }
            let stress_u = stress(lambda, mu, eu);
            let stress_p = stress(lambda, mu, ep);
            // A = sum_c grad(phi_c) tensor grad(phi_c); its engineering
            // contraction [Axx, Ayy, 2 Axy] pairs with physical shear stress.
            let mut gradient = [[0.0; 2]; 2];
            for (a, &node) in conn.iter().enumerate().skip(1) {
                for c in 0..2 {
                    let delta = mode[node][c] - mode[conn[0]][c];
                    for j in 0..2 {
                        gradient[c][j] += grads[a][j] * delta;
                    }
                }
            }
            let a = [
                gradient[0][0].powi(2) + gradient[1][0].powi(2),
                gradient[0][1].powi(2) + gradient[1][1].powi(2),
                2.0 * (gradient[0][0] * gradient[0][1] + gradient[1][0] * gradient[1][1]),
            ];
            let pullback = stress(lambda * scales[element], mu * scales[element], a);
            direct[element] += w * load.mul_add(dot3(stress_u, a), dot3(stress_p, ep));
            // The transpose of the same B-bar operator is essential here;
            // a pointwise B would reintroduce the hydrostatic locking term.
            for (local, &node) in conn.iter().enumerate() {
                for c in 0..2 {
                    let g = grads[local];
                    let mut b = if c == 0 {
                        [g[0], 0.0, g[1]]
                    } else {
                        [0.0, g[1], g[0]]
                    };
                    if problem.formulation == Formulation::BBar {
                        let correction = 0.5 * (mean_grad[local][c] - g[c]);
                        b[0] += correction;
                        b[1] += correction;
                    }
                    forces[element][2 * local + c] += w * dot3(b, stress_u);
                    if let Some(free) = map[2 * node + c] {
                        rhs[free] += w * dot3(b, pullback);
                    }
                }
            }
        }
    }
    if rhs
        .iter()
        .chain(&direct)
        .chain(forces.iter().flatten())
        .any(|x| !x.is_finite())
    {
        return Err(invalid("buckling sensitivity contractions must be finite"));
    }
    Ok((rhs, direct, forces))
}

fn stress(lambda: f64, mu: f64, strain: [f64; 3]) -> [f64; 3] {
    let mean = (lambda + mu) * (strain[0] + strain[1]);
    let deviator = mu * (strain[0] - strain[1]);
    [mean + deviator, mean - deviator, mu * strain[2]]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
