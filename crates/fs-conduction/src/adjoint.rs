//! Conductivity-design gradients and discrete thermal goal comparisons.
//!
//! [`ConductivityDesign`] computes `dJ/dρ` for a LINEAR conduction solve
//! whose design parameters are per-element conductivity multipliers. Fixed
//! matching or nonmatching thermal contacts can be retained through
//! [`ConductivityDesign::with_interfaces`]; their resistance is not scaled by ρ.
//! [`compare_discrete_goal`] also supports nonlinear k(T) and matching
//! contacts, retaining the finite linearization remainder explicitly.
//!
//! # What this establishes
//!
//! The design map is `K_e ← ρ_e · K` on the conduction block only.
//! With `R(T; ρ) = A(ρ)T − b(ρ) = 0` on the free dofs (`b` carries the
//! Dirichlet lift, which is itself `ρ`-dependent — hence the pullback
//! uses the FULL nodal temperature, not the free part), the implicit
//! function theorem gives
//!
//! ```text
//!   dJ/dρ_e = ∂J/∂ρ_e − λᵀ (∂R/∂ρ_e),   (∂R/∂T)ᵀ λ = ∂J/∂T
//!   ∂R/∂ρ_e = K_e · T_full |_free       (exact: A is AFFINE in ρ)
//! ```
//!
//! First- and second-order goals share one residual-checked, preconditioned
//! solve on `fs-solver`'s stack. The linear conductivity operator is symmetric,
//! so its adjoint uses the same PCG operator, never differentiated iterations.
//!
//! # What this does NOT establish
//!
//! - The conductivity DESIGN hook refuses nonlinear `k(T)`. The separate
//!   Robin/load responses and two-field goal comparison use
//!   [`crate::assemble_jacobian`]'s nonsymmetric operator, but supply no
//!   nonlinear conductivity-parameter gradient.
//! - Nothing about SHAPE derivatives. `ρ` is a coefficient, not a
//!   geometry; mesh-motion sensitivity is `fs-adjoint`'s Hadamard path
//!   and is not wired here.
//! - No continuum goal-error bound. [`ConductivityDesign::discrete_goal_error`]
//!   measures a linear goal discrepancy against this mesh's discrete solution.
//!   Mesh error requires a separately constructed enriched space; a residual
//!   on the original solved space is not a discretization estimate.

use fs_adjoint::AdjointReport;
use fs_exec::Cx;

/// Existing sparse-proof controls and outcomes used by the public Robin
/// analyzer. Consumers need no direct solver dependency to request this proof.
pub use fs_solver::goal::inverse::spectral::{
    SpectralInverseLimits, SpectralPreparation, SpectralStop,
};

use crate::ConductionError;
use crate::assemble::{
    DofMap, assemble_operator_scaled_with_interfaces, reduce,
};
use crate::interface::ThermalInterfaces;
use crate::material::{ElementMaterials, TemperatureSpan};
use crate::solve::{ConductionProblem, LinearConfig};

/// A conduction problem parameterized by per-element conductivity
/// multipliers, with an adjoint gradient.
pub struct ConductivityDesign<'m> {
    problem: ConductionProblem<'m>,
    dofs: DofMap,
    linear: LinearConfig,
    element_materials: Option<&'m ElementMaterials>,
    /// Fixed contact conductances, independent of the bulk design multipliers.
    interfaces: Option<&'m ThermalInterfaces>,
    /// Per-element LINEAR conductivity tensor at construction. A
    /// temperature-dependent model is refused, so this is not a
    /// linearization of `k(T)`.
    element_tensors: Vec<[[f64; 3]; 3]>,
    /// A material-admissible assembly field. A constant sampled table still
    /// has a bounded validity interval and must not be sampled at 0 K.
    reference_temperature: Vec<f64>,
}

/// A primal solve at one design point.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignSolution {
    /// Free-dof temperature, K.
    pub free_temperature: Vec<f64>,
    /// Full nodal temperature, K.
    pub temperature: Vec<f64>,
    /// Recomputed `‖b − Ax‖₂/‖b‖₂` of the primal solve.
    pub primal_relative_residual: f64,
    /// Krylov iterations spent on the primal solve.
    pub primal_iterations: usize,
}

/// Dual-weighted discrepancy to this mesh's linear discrete solution.
/// This is not a continuum error certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscreteGoalError {
    /// `J(T_discrete) - J(T_approximate)`, up to dual solver error.
    pub signed_error: f64,
    /// Signed `z_i (b - A T_approximate)_i`, in free-dof order and goal units.
    /// These are algebraic contributions, not elementwise refinement indicators.
    pub free_dof_contributions: Vec<f64>,
    /// Recomputed relative residual of `Aᵀ z = weights`.
    pub dual_relative_residual: f64,
    /// Krylov iterations used by the dual solve.
    pub dual_iterations: usize,
}

impl<'m> ConductivityDesign<'m> {
    /// Bind a design problem. A heterogeneous base tensor comes from
    /// the problem's own [`ConductionProblem::element_materials`] root;
    /// `ρ_e` still multiplies that element's own linear `K_e`. Every
    /// assigned model must be temperature-independent.
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] when the effective material is
    /// temperature dependent (this hook is the LINEAR case only);
    /// assignment validation refusals; [`ConductionError::NoFreeDofs`] /
    /// [`ConductionError::SingularPureNeumann`] from the degree-of-freedom
    /// map.
    pub fn new(
        problem: ConductionProblem<'m>,
        linear: LinearConfig,
    ) -> Result<ConductivityDesign<'m>, ConductionError> {
        if problem.boundary.has_temperature_dependent_surface_conduction() {
            return Err(ConductionError::Conductivity {
                what: "the linear IFT hook requires temperature-independent sheet conductivity; use the nonlinear Jacobian response".into(),
            });
        }
        let materials = problem.element_materials;
        if let Some(assigned) = materials {
            assigned.validate_for(problem.mesh)?;
        } else if problem.material.is_temperature_dependent() {
            return Err(ConductionError::Conductivity {
                what: "the IFT hook in this crate covers the LINEAR case only; a k(T) model \
                       needs the nonsymmetric Newton Jacobian and its own verification"
                    .to_string(),
            });
        }
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        let ne = problem.mesh.element_count();
        let mut element_tensors = Vec::with_capacity(ne);
        let mut vertex_spans = vec![TemperatureSpan::Unbounded; problem.mesh.vertex_count()];
        let reference = |span: TemperatureSpan| -> Result<f64, ConductionError> {
            let temperature = match span {
                TemperatureSpan::Unbounded => 0.0,
                TemperatureSpan::Sampled { low, high } => f64::midpoint(low, high),
            };
            span.check(temperature)?;
            Ok(temperature)
        };
        for e in 0..ne {
            let model = match materials {
                Some(assigned) => assigned.model_for(e)?,
                None => problem.material,
            };
            if model.is_temperature_dependent() {
                return Err(ConductionError::Conductivity {
                    what: "the IFT hook in this crate covers the LINEAR case only; a k(T) model \
                           needs the nonsymmetric Newton Jacobian and its own verification"
                        .to_string(),
                });
            }
            let span = model.temperature_span();
            element_tensors.push(model.tensor_at(reference(span)?)?);
            for &vertex in &problem.mesh.complex().tets[e] {
                vertex_spans[vertex as usize] = vertex_spans[vertex as usize].intersect(span);
            }
        }
        for sheet in problem.boundary.surface_heat_conduction() {
            for vertex in problem.mesh.boundary()[sheet.boundary_face()].vertices {
                vertex_spans[vertex as usize] = vertex_spans[vertex as usize].intersect(sheet.temperature_span());
            }
        }
        let reference_temperature = vertex_spans
            .into_iter()
            .map(reference)
            .collect::<Result<Vec<_>, _>>()?;
        Ok(ConductivityDesign {
            problem,
            dofs,
            linear,
            element_materials: materials,
            interfaces: None,
            element_tensors,
            reference_temperature,
        })
    }

    /// Bind fixed matching or nonmatching thermal contacts to the design.
    ///
    /// Their currently bound resistance coefficients enter every primal,
    /// adjoint, discrete-goal and Hessian solve, including the Dirichlet lift.
    /// Only the bulk tensors are multiplied by `ρ_e`; the contact operator
    /// and its geometry remain fixed. Temperature-, pressure- or
    /// design-dependent contact-law updates are not differentiated here.
    ///
    /// # Errors
    /// Refuses an interface bound to another mesh or incompatible boundary,
    /// using the same validation as the native contact assembler.
    pub fn with_interfaces(
        mut self,
        interfaces: &'m ThermalInterfaces,
    ) -> Result<Self, ConductionError> {
        interfaces.validate_for(self.problem.mesh, self.problem.boundary)?;
        self.interfaces = Some(interfaces);
        Ok(self)
    }

    /// One design parameter per element.
    #[must_use]
    pub fn parameter_count(&self) -> usize {
        self.problem.mesh.element_count()
    }

    /// The degree-of-freedom map.
    #[must_use]
    pub const fn dofs(&self) -> &DofMap {
        &self.dofs
    }

    fn check_rho(&self, rho: &[f64]) -> Result<(), ConductionError> {
        if rho.len() != self.parameter_count() {
            return Err(ConductionError::FieldLength {
                field: "conductivity design vector",
                expected: self.parameter_count(),
                found: rho.len(),
            });
        }
        for &r in rho {
            if !(r.is_finite() && r > 0.0) {
                return Err(ConductionError::Conductivity {
                    what: format!("design multiplier {r} must be finite and positive"),
                });
            }
        }
        Ok(())
    }

    /// Solve the primal problem at a design point.
    ///
    /// # Errors
    /// Assembly and linear-solve refusals.
    pub fn solve(&self, cx: &Cx<'_>, rho: &[f64]) -> Result<DesignSolution, ConductionError> {
        self.check_rho(rho)?;
        let (a_ff, b_f) = self.system(cx, rho)?;
        let (x, report) = second_order::solve(cx, &a_ff, &b_f, self.linear)?;
        let temperature = self.dofs.scatter(&x);
        self.check_temperature_support(cx, &temperature)?;
        Ok(DesignSolution {
            temperature,
            free_temperature: x,
            primal_relative_residual: report.relative_residual,
            primal_iterations: report.iterations,
        })
    }

    fn system(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
    ) -> Result<(fs_sparse::Csr, Vec<f64>), ConductionError> {
        let system = assemble_operator_scaled_with_interfaces(
            cx,
            self.problem.mesh,
            self.problem.boundary,
            self.problem.material,
            self.problem.source,
            &self.reference_temperature,
            Some(rho),
            self.interfaces,
            self.element_materials,
        )?;
        crate::solve::require_anchored_components(cx, self.problem, &self.dofs, &system.operator)?;
        Ok(reduce(&system, &self.dofs))
    }

    fn check_temperature_support(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<(), ConductionError> {
        // A P1 field reaches its extrema at vertices. A supported element
        // average cannot license extrapolation at a hot/cold corner. Check the
        // actual material of every incident element, including fixed vertices.
        for e in 0..self.problem.mesh.element_count() {
            if e % crate::ASSEMBLY_TILE == 0 {
                cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                    stage: "conductivity-design-material-support",
                    at: e,
                })?;
            }
            let model = match self.element_materials {
                Some(assigned) => assigned.model_for(e)?,
                None => self.problem.material,
            };
            for vertex in self.problem.mesh.complex().tets[e] {
                let t = temperature[vertex as usize];
                crate::require_finite("conductivity design temperature", t)?;
                model.temperature_span().check(t)?;
            }
        }
        // A constant sheet can have bounded source support independently of
        // the bulk. Its real trace temperatures must remain supported too.
        for (index, sheet) in self.problem.boundary.surface_heat_conduction().iter().enumerate() {
            if index % crate::ASSEMBLY_TILE == 0 {
                cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                    stage: "conductivity-design-sheet-support",
                    at: index,
                })?;
            }
            let span = sheet.temperature_span();
            for vertex in self.problem.mesh.boundary()[sheet.boundary_face()].vertices {
                span.check(temperature[vertex as usize])?;
            }
        }
        Ok(())
    }

    /// Evaluate a linear goal's discrepancy without another primal solve.
    ///
    /// Both vectors use this mesh's free-dof order. Prescribed temperatures
    /// enter through the primal's Dirichlet lift. For `Aᵀ z = weights`,
    /// `J(T_discrete) - J(T_approximate) = zᵀ(b - A T_approximate)`.
    /// The SPD linear operator makes the transposed solve identical to PCG.
    ///
    /// An enriched-mesh caller must prolong the coarse field and assemble the
    /// SAME goal there. This method neither constructs that enrichment nor
    /// bounds the remaining continuum error or dual error.
    ///
    /// # Errors
    /// Refuses invalid lengths, nonfinite inputs/results, a stalled dual,
    /// cancellation, and every assembly refusal of [`Self::solve`].
    pub fn discrete_goal_error(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
        approximate_free_temperature: &[f64],
        weights: &[f64],
    ) -> Result<DiscreteGoalError, ConductionError> {
        if !self.linear.tolerance.is_finite() || self.linear.tolerance <= 0.0 {
            return Err(ConductionError::Config {
                parameter: "linear.tolerance",
                what: "dual relative tolerance must be finite and positive".to_string(),
            });
        }
        self.check_rho(rho)?;
        for (field, values) in [
            ("approximate free temperature", approximate_free_temperature),
            ("goal weights", weights),
        ] {
            if values.len() != self.dofs.n() {
                return Err(ConductionError::FieldLength {
                    field,
                    expected: self.dofs.n(),
                    found: values.len(),
                });
            }
            if let Some(value) = values.iter().find(|value| !value.is_finite()) {
                return Err(ConductionError::NonFinite {
                    field,
                    bits: value.to_bits(),
                });
            }
        }
        self.check_temperature_support(cx, &self.dofs.scatter(approximate_free_temperature))?;
        let (matrix, rhs) = self.system(cx, rho)?;
        let checkpoint = |at| {
            cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                stage: "discrete-goal-dual",
                at,
            })
        };
        checkpoint(0)?;
        // Scale the goal before Krylov norms to avoid weight under/overflow.
        let scale = weights.iter().map(|w| w.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Ok(DiscreteGoalError {
                signed_error: 0.0,
                free_dof_contributions: vec![0.0; self.dofs.n()],
                dual_relative_residual: 0.0,
                dual_iterations: 0,
            });
        }
        let dual_rhs: Vec<f64> = weights.iter().map(|w| w / scale).collect();
        let (dual, report) = second_order::solve(cx, &matrix, &dual_rhs, self.linear)?;
        let mut applied = vec![0.0; self.dofs.n()];
        matrix.spmv(approximate_free_temperature, &mut applied);
        let free_dof_contributions: Vec<f64> = dual
            .iter()
            .zip(rhs.iter().zip(&applied))
            .map(|(z, (b, a))| (z * (b - a)) * scale)
            .collect();
        let signed_error: f64 = free_dof_contributions.iter().sum();
        if let Some(value) = free_dof_contributions
            .iter()
            .chain(std::iter::once(&signed_error))
            .find(|value| !value.is_finite())
        {
            return Err(ConductionError::NonFinite {
                field: "discrete goal error",
                bits: value.to_bits(),
            });
        }
        checkpoint(report.iterations)?;
        Ok(DiscreteGoalError {
            signed_error,
            free_dof_contributions,
            dual_relative_residual: report.relative_residual,
            dual_iterations: report.iterations,
        })
    }

    /// The linear quantity of interest `J(ρ) = Σ_i w_i T_i` over FREE
    /// dofs.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`] for a mis-sized weight vector;
    /// every refusal [`ConductivityDesign::solve`] can produce.
    pub fn objective(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
        weights: &[f64],
    ) -> Result<f64, ConductionError> {
        if weights.len() != self.dofs.n() {
            return Err(ConductionError::FieldLength {
                field: "objective weights",
                expected: self.dofs.n(),
                found: weights.len(),
            });
        }
        for &weight in weights {
            crate::require_finite("objective weights", weight)?;
        }
        let solution = self.solve(cx, rho)?;
        let objective = solution
            .free_temperature
            .iter()
            .zip(weights)
            .map(|(t, w)| t * w)
            .sum();
        crate::require_finite("conductivity objective", objective)?;
        Ok(objective)
    }

    /// `dJ/dρ` by the IFT adjoint, plus the adjoint-solve report.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`] for a mis-sized weight vector;
    /// every refusal [`ConductivityDesign::solve`] can produce.
    pub fn gradient(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
        weights: &[f64],
    ) -> Result<(Vec<f64>, AdjointReport), ConductionError> {
        if weights.len() != self.dofs.n() {
            return Err(ConductionError::FieldLength {
                field: "objective weights",
                expected: self.dofs.n(),
                found: weights.len(),
            });
        }
        // Reuse the retained goal's primal, transpose and material pullback.
        // This assembles once and enforces the actual iteration allowance,
        // cancellation, finite arithmetic and remeasured residual for BOTH
        // solves; a failed adjoint can never publish a usable gradient.
        let goal = self.linearize_goal(cx, rho, weights)?;
        let dual = goal.dual_solve();
        Ok((
            goal.gradient().to_vec(),
            AdjointReport {
                iters: dual.iterations,
                adjoint_residual: dual.relative_residual,
                converged: true,
            },
        ))
    }

    /// `(∂R/∂ρ)ᵀ λ`, element by element: `λᵀ (K_e · T_full)|_free`.
    /// Exact, because the bulk operator is LINEAR in `ρ`; fixed contact and
    /// boundary operators have zero parameter derivative.
    #[must_use]
    pub fn parameter_pullback(&self, lambda: &[f64], temperature: &[f64]) -> Vec<f64> {
        let mesh = self.problem.mesh;
        let mut out = vec![0.0f64; mesh.element_count()];
        for (e, slot) in out.iter_mut().enumerate() {
            let tet = mesh.complex().tets[e];
            let load = crate::assemble::element_conduction_action(
                mesh,
                e,
                &self.element_tensors[e],
                temperature,
            );
            let mut acc = 0.0f64;
            for a in 0..4 {
                let ia = self.dofs.slot_of(tet[a] as usize);
                let Some(ia) = ia else { continue };
                acc = lambda[ia].mul_add(load[a], acc);
            }
            *slot = acc;
        }
        out
    }
}

/// Robin-boundary and assembled-load tangents and adjoints, including k(T).
pub mod robin;

// Keep the public response types beside the goal analyzers used by native
// radiation, natural-convection and nonlinear-solid consumers. These are the
// same types as adjoint::robin, not alternative numerical implementations.
pub use robin::{MeanFluxFeedback, PrescribedFluxControls, PrescribedFluxGradient,
    RobinGradient, RobinResponse};

pub use robin::goal::{
    DiscreteGoalComparison, LinearGoalAnalysis, LinearGoalAnalysisConfig, LinearGoalAnalyzer,
    LinearGoalSolve, LinearGoalSolveConfig, LinearGoalStop, LinearMaximumAnalysis, LinearMaximumSolve,
    LinearRobinFeedbackAnalyzer, LinearRobinMaximumAnalysis, RobinFeedbackAnalysisConfig,
    RobinGoalFeedback, RobinGoalLinearization, analyze_linear_goal, analyze_linear_maximum,
    compare_discrete_goal,
};

/// Retained linear thermal goals and implicit second-order design derivatives.
pub mod second_order;

