//! Spatial load control on the original Kelvin/Maxwell constitutive histories.
//!
//! A prescribed affine displacement ramp gives an affine endpoint force map.
//! Its exact consistent tangent is SPD after the existing rigid-mode admission,
//! so one bounded sparse CG solve enforces the requested endpoint nodal loads.
//! This is a quasistatic time discretization: the material ramp is integrated
//! exactly, while a continuously prescribed load history requires time-step
//! convergence. No inertia, numerical damping or replacement loss law enters.

use super::*;
use crate::linear::Jacobi;
use crate::linear3::{TetStaticSolveConfig, euclidean_norm};
use fs_solver::krylov::CgState;
use fs_solver::op::{CsrOp, LinearOp};

/// Reusable exact-ramp stiffness for one system and one positive duration.
/// Preparation retains all source tensors, quadrature and off-diagonal terms.
/// Reuse this operator for repeated creep steps with the same duration; changes
/// in load, support displacement or accepted memory require no reassembly.
pub struct ViscoQuasistaticOperator<'a> {
    system: &'a ViscoElements,
    dt_s: f64,
    tangents: Vec<Mandel6>,
    reduced: Vec<Option<usize>>,
    operator: CsrOp,
    preconditioner: Jacobi,
}

impl ViscoQuasistaticOperator<'_> {
    /// Positive interval duration used by the constitutive tangent [s].
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }

    /// Full sparse consistent tangent in canonical free-DOF order [N/m].
    pub fn stiffness(&self) -> &Csr {
        self.operator.matrix()
    }
}

/// Accepted load-controlled endpoint, with original constitutive work and heat.
pub struct ViscoQuasistaticStep {
    /// Full node-major displacement, including moving prescribed supports [m].
    pub displacement_m: Vec<f64>,
    /// Internal force minus applied force on every DOF [N]. Fixed rows are
    /// support reactions; free rows retain the actual equilibrium residual.
    pub reaction_force_n: Vec<f64>,
    /// Exact material-ramp forces, physical storage, work and dashpot heat.
    pub constitutive: ViscoForceStep,
    /// Sparse CG iterations used in this step.
    pub iterations: usize,
    /// Independently recomputed endpoint force residual, divided by the largest
    /// free-force norm among applied load, history forcing and solve RHS.
    pub true_relative_residual: f64,
}

impl ViscoElements {
    /// Prepare the consistent tangent for load-controlled quasistatic creep.
    /// Every disconnected component must have its six rigid motions fixed,
    /// using the existing geometric solid admission. Arbitrary nonzero support
    /// displacements can subsequently be prescribed on those same fixed DOFs.
    /// Every original material band must include DC for sustained-load creep;
    /// this does not certify the transient spectrum of a supplied load history.
    pub fn quasistatic_operator(
        &self,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ViscoQuasistaticOperator<'_>, ViscoElementError> {
        checkpoint(cx)?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(ViscoElementError::InvalidHistory {
                what: "quasistatic creep requires finite positive duration",
            });
        }
        if let Some(error) = &self.static_constraint_error {
            return Err(error.clone().into());
        }
        let n = self.equilibrium.free_dofs.len();
        let mut reduced = vec![None; self.full_dofs];
        for (i, &dof) in self.equilibrium.free_dofs.iter().enumerate() {
            reduced[dof] = Some(i);
        }
        let mut tangents = Vec::with_capacity(self.elements.len());
        let mut stiffness = Coo::new(n, n);
        for element in &self.elements {
            checkpoint(cx)?;
            // The original spectral owner admits the equilibrium limit. A
            // narrow-band fit cannot silently become a long-term creep law.
            let _ = element.law.complex_modulus(0.0)?;
            let tangent = element
                .law
                .tangent_mandel(dt_s)
                .map_err(ViscoElementError::Projection)?;
            let d = element.dofs.len();
            let mut local = vec![0.0; d * d];
            for point in &element.points {
                checkpoint(cx)?;
                for i in 0..d {
                    for j in i..d {
                        local[i * d + j] +=
                            point.volume_m3 * bilinear(&point.b[i], &tangent, &point.b[j]);
                    }
                }
            }
            for i in 0..d {
                for j in i..d {
                    let value = local[i * d + j];
                    if !value.is_finite() {
                        return Err(ViscoElementError::NonFiniteResult);
                    }
                    if let (Some(row), Some(col)) =
                        (reduced[element.dofs[i]], reduced[element.dofs[j]])
                    {
                        stiffness.push(row, col, value);
                        if row != col {
                            stiffness.push(col, row, value);
                        }
                    }
                }
            }
            tangents.push(tangent);
        }
        let stiffness = stiffness.assemble();
        if (0..n).any(|i| stiffness.row(i).1.iter().any(|v| !v.is_finite())) {
            return Err(ViscoElementError::NonFiniteResult);
        }
        checkpoint(cx)?;
        let preconditioner = Jacobi::new(&stiffness);
        Ok(ViscoQuasistaticOperator {
            system: self,
            dt_s,
            tangents,
            reduced,
            operator: CsrOp::symmetric(stiffness),
            preconditioner,
        })
    }
}

impl<'a> ViscoElementState<'a> {
    /// Solve complete endpoint nodal loads with prescribed support motion.
    /// `support_displacement_m` is a full node-major vector with zero entries
    /// on free DOFs. Loads on fixed DOFs are retained in reaction recovery.
    /// The operator must borrow this exact original system.
    ///
    /// Each accepted displacement interpolates affinely from the prior state.
    /// Consequently reported work is the material's integrated stress work,
    /// including support motion, not a trapezoid of endpoint nodal loads.
    /// Endpoint force balance does not certify the continuous load trajectory
    /// between endpoints. Refine the interval to resolve creep/load evolution.
    /// No state is published until both sparse and original-force residuals
    /// pass, every source strain is admitted, and cancellation is checked.
    pub fn step_quasistatic(
        &mut self,
        prepared: &ViscoQuasistaticOperator<'a>,
        applied_force_n: &[f64],
        support_displacement_m: &[f64],
        config: TetStaticSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<ViscoQuasistaticStep, ViscoElementError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.system, prepared.system)
            || applied_force_n.len() != self.system.full_dofs
            || support_displacement_m.len() != self.system.full_dofs
            || applied_force_n
                .iter()
                .chain(support_displacement_m)
                .any(|v| !v.is_finite())
            || self
                .system
                .equilibrium
                .free_dofs
                .iter()
                .any(|&d| support_displacement_m[d] != 0.0)
        {
            return Err(ViscoElementError::InvalidHistory {
                what: "same-system operator, complete finite loads and fixed-DOF-only support motion required",
            });
        }
        if !config.relative_residual_tolerance.is_finite()
            || config.relative_residual_tolerance <= 0.0
            || config.relative_residual_tolerance >= 1.0
            || config.maximum_iterations == 0
        {
            return Err(TetElasticError::InvalidStaticSolveConfig {
                what: "positive iteration budget and relative residual tolerance in (0,1) required",
            }
            .into());
        }
        let free = &self.system.equilibrium.free_dofs;
        let mut history_force = vec![0.0; free.len()];
        for ((element, histories), tangent) in self
            .system
            .elements
            .iter()
            .zip(&self.histories)
            .zip(&prepared.tangents)
        {
            for (point, accepted) in element.points.iter().zip(histories) {
                checkpoint(cx)?;
                // Hold the admitted previous strain through the original
                // constitutive owner. The tangent then supplies the exact
                // affine offset for moving supports. Solve free displacement
                // increments around the accepted configuration so the large
                // relaxed elastic load need not be solved again at each step.
                let mut hold = accepted.clone();
                let mut stress = element
                    .law
                    .step(&mut hold, accepted.strain, prepared.dt_s)?;
                let delta: [f64; 6] = core::array::from_fn(|i| {
                    point
                        .b
                        .iter()
                        .zip(&element.dofs)
                        .filter(|(_, d)| prepared.reduced[**d].is_none())
                        .map(|(b, &d)| b[i] * (support_displacement_m[d] - self.displacement_m[d]))
                        .sum::<f64>()
                });
                for i in 0..6 {
                    stress[i] += (0..6).map(|j| tangent[i][j] * delta[j]).sum::<f64>();
                }
                for (b, &dof) in point.b.iter().zip(&element.dofs) {
                    if let Some(row) = prepared.reduced[dof] {
                        history_force[row] += point.volume_m3
                            * b.iter()
                                .zip(stress)
                                .map(|(b, stress)| b * stress)
                                .sum::<f64>();
                    }
                }
            }
        }
        let rhs: Vec<_> = free
            .iter()
            .zip(&history_force)
            .map(|(&d, &history)| applied_force_n[d] - history)
            .collect();
        let force_scale = euclidean_norm(rhs.iter().copied())
            .max(euclidean_norm(free.iter().map(|&d| applied_force_n[d])))
            .max(euclidean_norm(history_force.iter().copied()))
            .max(f64::MIN_POSITIVE);
        if !force_scale.is_finite() || rhs.iter().chain(&history_force).any(|v| !v.is_finite()) {
            return Err(ViscoElementError::NonFiniteResult);
        }
        let mut solve = CgState::new(&prepared.operator, &prepared.preconditioner, &rhs);
        while solve.iters < config.maximum_iterations
            && solve.rel_residual() >= config.relative_residual_tolerance
        {
            checkpoint(cx)?;
            let remaining = config.maximum_iterations - solve.iters;
            let _ = solve.run(
                &prepared.operator,
                &prepared.preconditioner,
                config.relative_residual_tolerance,
                remaining.min(64),
            );
        }
        let mut applied = vec![0.0; free.len()];
        prepared.operator.apply(&solve.x, &mut applied);
        let residual = euclidean_norm(rhs.iter().zip(&applied).map(|(b, ax)| b - ax)) / force_scale;
        check_residual(residual, solve.iters, config)?;

        let mut displacement_m = support_displacement_m.to_vec();
        for (&dof, &value) in free.iter().zip(&solve.x) {
            displacement_m[dof] = self.displacement_m[dof] + value;
        }
        let mut candidate = self.clone();
        let constitutive = candidate.step_internal_force(&displacement_m, prepared.dt_s, cx)?;
        let reaction_force_n: Vec<_> = constitutive
            .full_force_n
            .iter()
            .zip(applied_force_n)
            .map(|(internal, external)| internal - external)
            .collect();
        if reaction_force_n.iter().any(|v| !v.is_finite()) {
            return Err(ViscoElementError::NonFiniteResult);
        }
        let true_relative_residual =
            euclidean_norm(free.iter().map(|&d| reaction_force_n[d])) / force_scale;
        check_residual(true_relative_residual, solve.iters, config)?;
        checkpoint(cx)?;
        *self = candidate;
        Ok(ViscoQuasistaticStep {
            displacement_m,
            reaction_force_n,
            constitutive,
            iterations: solve.iters,
            true_relative_residual,
        })
    }
}

fn check_residual(
    residual: f64,
    iterations: usize,
    config: TetStaticSolveConfig,
) -> Result<(), ViscoElementError> {
    if !residual.is_finite() || residual > config.relative_residual_tolerance {
        Err(TetElasticError::StaticSolveFailed {
            iterations,
            true_relative_residual: residual,
        }
        .into())
    } else {
        Ok(())
    }
}
