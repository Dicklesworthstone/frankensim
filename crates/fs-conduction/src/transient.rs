//! Transient conduction by the method of lines: a declared volumetric heat
//! capacity, the exact P1 capacitance matrix, and theta-method stepping over
//! the crate's existing steady operator.
//!
//! Bead `frankensim-extreal-program-f85xj.5.11`, the first item staged by
//! `f85xj.5.9`. The steady solve answers "where does it settle"; this answers
//! "how long does it take to get there, and how hot does it get on the way",
//! which is the question a duty cycle actually asks.
//!
//! # The discrete problem
//!
//! Semi-discretizing in space leaves an ODE system in the nodal temperatures,
//!
//! ```text
//!   C dT/dt + K(T) T = b
//! ```
//!
//! with `K` and `b` exactly the operator and load the steady path assembles,
//! and `C` the capacitance matrix built from the same P1 mass matrix the
//! source load already uses, `∫ λ_a λ_b dV = V(1+δ_ab)/20`, scaled by the
//! declared volumetric heat capacity `rho·c_p`. Reusing the steady assembly
//! is the point: a transient run cannot drift from the steady physics,
//! because it IS the steady physics with one extra term.
//!
//! The theta method advances it as
//!
//! ```text
//!   (C/dt + θK) T^{n+1} = (C/dt − (1−θ)K) T^n + b
//! ```
//!
//! # Only unconditionally stable schemes are admitted
//!
//! `theta` is restricted to `[0.5, 1]`. Below one half the scheme is only
//! conditionally stable, and the step size at which it stops being stable
//! depends on the mesh — so an explicit-leaning theta would hand callers a
//! configuration that silently produces oscillating garbage on a refined
//! mesh while looking fine on a coarse one. Callers who want that tradeoff
//! should get it from a scheme that states its own CFL condition, not from a
//! parameter that looks continuous.
//!
//! `theta = 1` is backward Euler (first order, strongly damping) and
//! `theta = 0.5` is Crank-Nicolson (second order, non-damping). Both are
//! offered because the choice is a real modelling decision: Crank-Nicolson
//! is more accurate and will ring on a step load that backward Euler smooths.
//!
//! # Linear conduction only
//!
//! A temperature-dependent `k(T)` is REFUSED rather than linearized at the
//! old step. Freezing conductivity across a step is a defensible scheme, but
//! it is a different one with its own error behaviour, and adopting it
//! silently would make the observed time order depend on how strongly `k`
//! varies. This mirrors `fs-adjoint`, which refuses a `k(T)` material rather
//! than linearizing it.

use crate::ConductionError;
use crate::assemble::{DofMap, assemble_operator, reduce_matrix_and_lift};
use crate::bc::ThermalBoundary;
use crate::field::ScalarField;
use crate::material::ConductivityModel;
use crate::mesh::ConductionMesh;
use crate::solve::{LinearConfig, spd_preconditioner};
use fs_exec::Cx;
use fs_solver::krylov::CgState;
use fs_solver::{CsrOp, norm2};
use fs_sparse::{Coo, Csr};

/// SI exponents of volumetric heat capacity, J/(m³·K).
pub const VOLUMETRIC_HEAT_CAPACITY_DIMS: fs_qty::Dims = fs_qty::Dims([-1, 1, -2, -1, 0, 0]);

/// Smallest admitted theta. Below one half the scheme is only conditionally
/// stable and its limit depends on the mesh.
pub const MIN_THETA: f64 = 0.5;

/// Declared volumetric heat capacity `rho · c_p`, J/(m³·K).
///
/// Declared rather than looked up, for now: a card-backed capacity with the
/// receipt discipline `ConductivityModel` already has is the natural next
/// step, and this type is shaped so that becoming card-backed does not change
/// its callers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumetricHeatCapacity {
    value_j_per_m3_k: f64,
}

impl VolumetricHeatCapacity {
    /// Declare a positive volumetric heat capacity.
    ///
    /// # Errors
    /// [`ConductionError::NonFinite`] for a non-finite value, and
    /// [`ConductionError::Config`] for a non-positive one — a zero capacity
    /// is not a fast material, it is a steady problem, and silently solving a
    /// different equation is the failure this refuses.
    pub fn declared(value_j_per_m3_k: f64) -> Result<Self, ConductionError> {
        crate::require_finite("volumetric heat capacity", value_j_per_m3_k)?;
        if value_j_per_m3_k <= 0.0 {
            return Err(ConductionError::Config {
                parameter: "volumetric heat capacity",
                what: format!(
                    "{value_j_per_m3_k} J/(m^3 K) is not positive; declare a positive rho*c_p, and for a zero-capacity model solve the steady problem instead"
                ),
            });
        }
        Ok(Self { value_j_per_m3_k })
    }

    /// The declared value, J/(m³·K).
    #[must_use]
    pub const fn value_j_per_m3_k(self) -> f64 {
        self.value_j_per_m3_k
    }
}

/// Assemble the P1 capacitance matrix `C_ab = ∫ rho c_p λ_a λ_b dV`.
///
/// Uses the same exact tet mass matrix `V(1+δ_ab)/20` the source load uses,
/// so the discrete heat content of a uniform temperature field is exactly
/// `rho c_p V T` — the property the uniform-heating test pins.
///
/// # Errors
/// [`ConductionError::Cancelled`] at a tile boundary.
pub fn assemble_capacitance(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    capacity: VolumetricHeatCapacity,
) -> Result<Csr, ConductionError> {
    let n = mesh.vertex_count();
    let mut coo = Coo::new(n, n);
    for element in 0..mesh.element_count() {
        if element % crate::assemble::ASSEMBLY_TILE == 0 {
            cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                stage: "assemble-capacitance",
                at: element,
            })?;
        }
        let volume = mesh.element_volume(element);
        let vertices = mesh.complex().tets[element];
        for (a, &va) in vertices.iter().enumerate() {
            for (b, &vb) in vertices.iter().enumerate() {
                let kronecker = if a == b { 1.0 } else { 0.0 };
                let entry = capacity.value_j_per_m3_k * volume * (1.0 + kronecker) / 20.0;
                coo.push(va as usize, vb as usize, entry);
            }
        }
    }
    Ok(coo.assemble())
}

/// Time-integration configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientConfig {
    theta: f64,
    dt_s: f64,
    linear: LinearConfig,
}

impl TransientConfig {
    /// Declare the scheme and step.
    ///
    /// # Errors
    /// [`ConductionError::Config`] for a theta outside `[MIN_THETA, 1]` or a
    /// non-positive/non-finite step.
    pub fn new(theta: f64, dt_s: f64, linear: LinearConfig) -> Result<Self, ConductionError> {
        if !theta.is_finite() || !(MIN_THETA..=1.0).contains(&theta) {
            return Err(ConductionError::Config {
                parameter: "theta",
                what: format!(
                    "{theta} is outside the admitted [{MIN_THETA}, 1]; use 1 (backward Euler) or 0.5 (Crank-Nicolson), because below 0.5 the scheme is only conditionally stable and its limit depends on the mesh"
                ),
            });
        }
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(ConductionError::Config {
                parameter: "time step",
                what: format!("{dt_s} s is not finite and positive"),
            });
        }
        Ok(Self {
            theta,
            dt_s,
            linear,
        })
    }

    /// Backward Euler at the given step.
    ///
    /// # Errors
    /// As [`TransientConfig::new`].
    pub fn backward_euler(dt_s: f64, linear: LinearConfig) -> Result<Self, ConductionError> {
        Self::new(1.0, dt_s, linear)
    }

    /// Crank-Nicolson at the given step.
    ///
    /// # Errors
    /// As [`TransientConfig::new`].
    pub fn crank_nicolson(dt_s: f64, linear: LinearConfig) -> Result<Self, ConductionError> {
        Self::new(MIN_THETA, dt_s, linear)
    }

    /// The declared theta.
    #[must_use]
    pub const fn theta(&self) -> f64 {
        self.theta
    }

    /// The declared step, s.
    #[must_use]
    pub const fn dt_s(&self) -> f64 {
        self.dt_s
    }

    /// The scheme's nominal temporal order: 2 at Crank-Nicolson, else 1.
    ///
    /// Nominal, not measured. `tests/transient.rs` measures the observed
    /// order and gates it; this is only what the scheme claims.
    #[must_use]
    pub fn nominal_order(&self) -> u32 {
        if (self.theta - 0.5).abs() < 1.0e-12 {
            2
        } else {
            1
        }
    }
}

/// One completed step's evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientStep {
    time_s: f64,
    krylov_iterations: usize,
    true_relative_residual: f64,
}

impl TransientStep {
    /// Simulated time at the END of this step, s.
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Krylov iterations the step consumed.
    #[must_use]
    pub const fn krylov_iterations(&self) -> usize {
        self.krylov_iterations
    }

    /// The crate's own recomputed Euclidean relative residual, not the
    /// method's recursive estimate.
    #[must_use]
    pub const fn true_relative_residual(&self) -> f64 {
        self.true_relative_residual
    }
}

/// A completed march.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientSolution {
    /// Nodal temperature over ALL vertices at the final time, K.
    pub temperature: Vec<f64>,
    /// Final simulated time, s.
    pub time_s: f64,
    /// Per-step evidence, in order.
    pub steps: Vec<TransientStep>,
}

/// The spatial problem a transient march advances.
///
/// Mirrors [`crate::solve::ConductionProblem`] and adds the capacity, so the
/// transient entry point takes the same shape as the steady one rather than a
/// long positional argument list.
#[derive(Debug, Clone, Copy)]
pub struct TransientProblem<'m> {
    /// The prepared mesh.
    pub mesh: &'m ConductionMesh,
    /// The boundary partition.
    pub boundary: &'m ThermalBoundary,
    /// The conductivity model. Must be temperature independent.
    pub material: &'m ConductivityModel,
    /// The volumetric source, W/m³.
    pub source: &'m ScalarField,
    /// Declared volumetric heat capacity.
    pub capacity: VolumetricHeatCapacity,
}

/// March a linear transient conduction problem for `steps` steps.
///
/// `initial` is the nodal temperature at `t = 0` over all vertices. Dirichlet
/// vertices are lifted to their prescribed values at every step, so an
/// initial condition inconsistent with the boundary is corrected at the first
/// step rather than carried.
///
/// # Errors
/// [`ConductionError::Config`] for a temperature-dependent conductivity or a
/// mismatched initial vector, [`ConductionError::LinearSolveFailed`] when a
/// step does not converge, and [`ConductionError::Cancelled`] at a tile or
/// step boundary.
pub fn march(
    cx: &Cx<'_>,
    problem: TransientProblem<'_>,
    config: &TransientConfig,
    initial: &[f64],
    steps: usize,
) -> Result<TransientSolution, ConductionError> {
    let TransientProblem {
        mesh,
        boundary,
        material,
        source,
        capacity,
    } = problem;
    if material.is_temperature_dependent() {
        return Err(ConductionError::Config {
            parameter: "conductivity",
            what: "transient marching received a temperature-dependent conductivity; supply a constant one, because freezing k(T) across a step is a different scheme with its own error behaviour and is not adopted silently".to_string(),
        });
    }
    let n = mesh.vertex_count();
    if initial.len() != n {
        return Err(ConductionError::FieldLength {
            field: "initial temperature",
            expected: n,
            found: initial.len(),
        });
    }
    for &value in initial {
        crate::require_finite("initial temperature", value)?;
    }

    let dofs = DofMap::new(boundary, n)?;
    let capacitance = assemble_capacitance(cx, mesh, capacity)?;

    // K and b are temperature independent here, so one assembly serves every
    // step. That is exactly why k(T) is refused rather than re-assembled.
    let system = assemble_operator(cx, mesh, boundary, material, source, initial)?;
    let stiffness = &system.operator;
    let load = &system.load;

    let inverse_dt = 1.0 / config.dt_s;
    let theta = config.theta;

    // A = C/dt + theta*K, assembled once.
    let lhs = axpy_csr(&capacitance, inverse_dt, stiffness, theta);
    let (a_ff, lift) = reduce_matrix_and_lift(&lhs, &dofs);
    let precond = spd_preconditioner(&a_ff);

    let mut temperature = initial.to_vec();
    for (index, &vertex) in dofs.fixed().iter().enumerate() {
        let _ = index;
        temperature[vertex] = dofs.prescribed()[vertex];
    }

    let mut recorded = Vec::with_capacity(steps);
    let mut scratch_c = vec![0.0f64; n];
    let mut scratch_k = vec![0.0f64; n];

    for step in 0..steps {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "transient-step",
            at: step,
        })?;

        // r = C/dt T^n − (1−θ) K T^n + b, over all vertices.
        capacitance.spmv(&temperature, &mut scratch_c);
        stiffness.spmv(&temperature, &mut scratch_k);
        let mut rhs = Vec::with_capacity(dofs.n());
        for (slot, &vertex) in dofs.free().iter().enumerate() {
            let value = scratch_c[vertex].mul_add(
                inverse_dt,
                (1.0 - theta).mul_add(-scratch_k[vertex], load[vertex]),
            ) + lift[slot];
            if !value.is_finite() {
                return Err(ConductionError::NonFinite {
                    field: "transient right-hand side",
                    bits: value.to_bits(),
                });
            }
            rhs.push(value);
        }

        let op = CsrOp::symmetric(a_ff.clone());
        let mut cg = CgState::new(&op, &precond, &rhs);
        let report = cg.run(
            &op,
            &precond,
            config.linear.tolerance,
            config.linear.max_iterations,
        );
        let truth = true_relative_residual(&a_ff, &rhs, &cg.x);
        if truth.is_nan() || truth >= config.linear.tolerance {
            return Err(ConductionError::LinearSolveFailed {
                iteration: step,
                krylov_iterations: report.iters,
                true_relative_residual: truth,
                tolerance: config.linear.tolerance,
            });
        }

        for (slot, &vertex) in dofs.free().iter().enumerate() {
            temperature[vertex] = cg.x[slot];
        }
        recorded.push(TransientStep {
            time_s: config.dt_s * ((step + 1) as f64),
            krylov_iterations: report.iters,
            true_relative_residual: truth,
        });
    }

    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "transient-publish",
        at: steps,
    })?;

    Ok(TransientSolution {
        temperature,
        time_s: config.dt_s * (steps as f64),
        steps: recorded,
    })
}

/// `alpha·A + beta·B` for two same-shaped CSR matrices, via the crate's
/// deterministic COO assembly.
fn axpy_csr(a: &Csr, alpha: f64, b: &Csr, beta: f64) -> Csr {
    let n = a.nrows();
    let mut coo = Coo::new(n, n);
    for row in 0..n {
        let (columns, values) = a.row(row);
        for (&column, &value) in columns.iter().zip(values.iter()) {
            coo.push(row, column, alpha * value);
        }
        let (columns, values) = b.row(row);
        for (&column, &value) in columns.iter().zip(values.iter()) {
            coo.push(row, column, beta * value);
        }
    }
    coo.assemble()
}

/// The recomputed Euclidean relative residual, not the method's recursive
/// estimate.
fn true_relative_residual(a: &Csr, b: &[f64], x: &[f64]) -> f64 {
    let mut ax = vec![0.0f64; b.len()];
    a.spmv(x, &mut ax);
    let mut residual = Vec::with_capacity(b.len());
    for (value, target) in ax.iter().zip(b.iter()) {
        residual.push(target - value);
    }
    let denominator = norm2(b);
    if denominator > 0.0 {
        norm2(&residual) / denominator
    } else {
        norm2(&residual)
    }
}
