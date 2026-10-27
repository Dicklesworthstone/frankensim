//! Sparse average-acceleration Newmark dynamics on the existing P1 operators.
//!
//! M a + (eta M + zeta K) v + K u = f_external(t) + f_thermal.
//! Topology and restraint components are fixed for the run; `thermal` admits
//! source-refreshed eigenstrain and stiffness at unchanged reference density,
//! while `motion` additionally admits prescribed support trajectories.
//! Free bodies are admitted: mass, not artificial clamps,
//! makes the acceleration solve nonsingular. Forces retain their original mesh.
//!
//! This is the sparse acceleration-form backend of the average-acceleration
//! scheme also used by `shell_time`: beta=1/4, gamma=1/2. No dense matrix or
//! effective-matrix assembly is needed. Element operators and PCG stay in their
//! existing owners. Every candidate passes a fresh linear residual and a discrete
//! work/energy check before replacing the accepted state. The latter is a
//! numerical identity, NOT a temporal, spatial or continuum error certificate.

use super::*;
use crate::linear3::TetElasticAssembly;
use crate::shell_time::DynamicState;
use fs_solver::krylov::CgState;
use fs_solver::op::LinearOp;
use fs_sparse::{Csr, precond::Precond};
use std::sync::Arc;

/// Physical stress, inertia and support reactions at accepted endpoints.
pub mod response;
/// Finite load jumps and bounded material reconfiguration at accepted times.
pub mod events;
/// Discrete state and load pullbacks through the original inertial equations.
pub mod adjoint;
/// Absolute prescribed support motion with full inertial coupling and work.
pub mod motion;
/// Source-refreshed thermal strain with physical elastic energy and input work.
pub mod thermal;
mod thermal_damping;

/// Explicit resource, solve, motion and energy allowances, all in coherent SI.
#[derive(Debug, Clone, Copy)]
pub struct TransientConfig3 {
    pub maximum_steps: usize,
    pub maximum_free_dofs: usize,
    /// Preflight bound on 144 scalar element entries per tetrahedron, per matrix.
    pub maximum_element_entries: usize,
    pub solve: TetStaticSolveConfig,
    pub absolute_force_tolerance_n: f64,
    pub absolute_energy_tolerance_j: f64,
    pub relative_energy_tolerance: f64,
    pub maximum_displacement_m: f64,
    /// eta >= 0 in C = eta M + zeta K, units s^-1.
    pub mass_damping_per_s: f64,
    /// zeta >= 0 in C = eta M + zeta K, units s.
    pub stiffness_damping_s: f64,
}
impl Default for TransientConfig3 {
    fn default() -> Self {
        Self {
            maximum_steps: 10_000, maximum_free_dofs: 50_000,
            maximum_element_entries: 14_400_000,
            solve: TetStaticSolveConfig::default(), absolute_force_tolerance_n: 1e-12,
            absolute_energy_tolerance_j: 1e-12, relative_energy_tolerance: 1e-8,
            maximum_displacement_m: 0.01, mass_damping_per_s: 0.0,
            stiffness_damping_s: 0.0,
        }
    }
}
impl TransientConfig3 {
    fn validate(self) -> Result<(), TransientError3> {
        let positive = [self.absolute_force_tolerance_n, self.absolute_energy_tolerance_j,
            self.relative_energy_tolerance, self.maximum_displacement_m,
            self.solve.relative_residual_tolerance];
        if positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || self.solve.relative_residual_tolerance >= 1.0
            || self.relative_energy_tolerance >= 1.0
            || self.maximum_steps == 0 || self.maximum_free_dofs == 0
            || self.maximum_element_entries == 0 || self.solve.maximum_iterations == 0
            || [self.mass_damping_per_s, self.stiffness_damping_s].iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        { return Err(TransientError3::Invalid("invalid transient resource/solve/energy/damping allowance")); }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum TransientError3 {
    Invalid(&'static str),
    Loads(LoadError3),
    Budget,
    Solve { iterations: usize, residual_n: f64, allowance_n: f64 },
    Momentum { residual_kg_m_s: f64, allowance_kg_m_s: f64 },
    Energy { defect_j: f64, allowance_j: f64 },
    Cancelled,
}
impl std::fmt::Display for TransientError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "tetrahedral transient refused: {self:?}")
    }
}
impl std::error::Error for TransientError3 {}
impl From<LoadError3> for TransientError3 {
    fn from(e: LoadError3) -> Self {
        if e == LoadError3::Cancelled { Self::Cancelled } else { Self::Loads(e) }
    }
}
impl From<TetElasticError> for TransientError3 {
    fn from(e: TetElasticError) -> Self { LoadError3::from(e).into() }
}
fn check(cx: &Cx<'_>) -> Result<(), TransientError3> {
    cx.checkpoint().map_err(|_| TransientError3::Cancelled)
}
fn number(v: f64) -> Result<f64, TransientError3> {
    if v.is_finite() { Ok(v) } else { Err(TransientError3::Invalid("nonfinite transient arithmetic")) }
}
fn inner(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (a, b)| a.mul_add(*b, sum))
}
fn length(a: &[f64]) -> f64 { a.iter().fold(0.0, |n, x| n.hypot(*x)) }

/// A linear combination applied without constructing another sparse matrix.
struct Operator<'a, 'cx> { m: &'a Csr, k: &'a Csr, cm: f64, ck: f64, cx: &'a Cx<'cx> }
impl LinearOp for Operator<'_, '_> {
    fn n(&self) -> usize { self.m.nrows() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut work = 0usize;
        for i in 0..self.n() {
            let mut value = 0.0;
            for (matrix, coefficient) in [(self.m, self.cm), (self.k, self.ck)] {
                if coefficient == 0.0 { continue; }
                let (columns, values) = matrix.row(i);
                let mut applied = 0.0;
                for (&column, &entry) in columns.iter().zip(values) {
                    if work % 256 == 0 && self.cx.checkpoint().is_err() {
                        // LinearOp has no error return. Poison this private
                        // candidate; the owner polls again before publication.
                        y.fill(f64::NAN); return;
                    }
                    applied = entry.mul_add(x[column], applied);
                    work += 1;
                }
                value = coefficient.mul_add(applied, value);
            }
            y[i] = value;
        }
    }
}
struct Diagonal(Vec<f64>);
impl Precond for Diagonal {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((z, r), inv) in z.iter_mut().zip(r).zip(&self.0) { *z = r * inv; }
    }
}
fn solve(op: &Operator<'_, '_>, rhs: &[f64], config: TransientConfig3)
    -> Result<(Vec<f64>, usize, f64), TransientError3> {
    solve_with_diagonal(op, rhs, config, op.cx,
        |i| op.cm * op.m.get(i, i) + op.ck * op.k.get(i, i))
}

// The same PCG, preconditioner and fresh residual gate also admit a combination
// of the two endpoint stiffness operators. No parallel numerical solver exists.
fn solve_with_diagonal(op: &impl LinearOp, rhs: &[f64], config: TransientConfig3,
    cx: &Cx<'_>, diagonal_entry: impl Fn(usize) -> f64)
    -> Result<(Vec<f64>, usize, f64), TransientError3> {
    check(cx)?;
    let rhs_norm = number(length(rhs))?;
    let allowance = number(config.absolute_force_tolerance_n
        + config.solve.relative_residual_tolerance * rhs_norm)?;
    if rhs_norm == 0.0 { return Ok((vec![0.0; rhs.len()], 0, 0.0)); }
    let mut inverse = Vec::with_capacity(rhs.len());
    for i in 0..rhs.len() {
        if i % 256 == 0 { check(cx)?; }
        let d = number(diagonal_entry(i))?;
        if d <= 0.0 { return Err(TransientError3::Invalid("effective mass diagonal is not positive")); }
        inverse.push(number(d.recip())?);
    }
    let diagonal = Diagonal(inverse);
    let mut state = CgState::new(op, &diagonal, rhs);
    // Solve more tightly than the acceptance gate, avoiding recursive-residual
    // drift at the boundary. Absolute allowance still handles physical zero.
    let target = (0.1 * allowance / rhs_norm).min(0.1);
    while state.iters < config.solve.maximum_iterations && state.rel_residual() > target {
        check(cx)?;
        let before = state.iters;
        let remaining = config.solve.maximum_iterations - before;
        let _ = state.run(op, &diagonal, target, remaining.min(32));
        if state.iters == before { break; }
    }
    check(cx)?;
    let mut applied = vec![0.0; rhs.len()]; op.apply(&state.x, &mut applied);
    check(cx)?;
    let residual = rhs.iter().zip(&applied).fold(0.0_f64, |n, (f, a)| n.hypot(f-a));
    if !residual.is_finite() || residual > allowance || state.x.iter().any(|v| !v.is_finite()) {
        return Err(TransientError3::Solve { iterations: state.iters, residual_n: residual, allowance_n: allowance });
    }
    Ok((state.x, state.iters, residual))
}

/// Mechanical source work includes the equivalent thermal load.
/// The reported energy is 1/2 v'Mv + 1/2 u'Ku, not the shifted thermodynamic
/// elastic energy of (strain - thermal_strain). That distinction is explicit.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientStep3 {
    pub time_s: f64,
    pub step_size_s: f64,
    pub iterations: usize,
    pub force_residual_n: f64,
    pub mechanical_work_j: f64,
    /// Work of support reactions ON the solid; zero for stationary zero supports.
    pub support_work_j: f64,
    pub thermal_equivalent_work_j: f64,
    /// Signed parameter work 1/2 q_old' (K_next-K_old) q_next [J].
    /// Zero for unchanged stiffness; this is not dissipated heat.
    pub stiffness_parameter_work_j: f64,
    pub damping_loss_j: f64,
    pub quadratic_energy_j: f64,
    pub balance_error_j: f64,
}

/// Immutable operators plus the last accepted state. Clone is a same-process
/// checkpoint/fork; matrices are shared through Arc, state vectors are copied.
/// No history grows with step count. Failed/cancelled steps change nothing.
#[derive(Clone)]
pub struct TetTransientRun3<'g, 'm, 'p, 's> {
    pub(super) geometry: &'g TetLoadGeometry3<'m>,
    pub(super) problem: &'p TetLinearElasticProblem<'s>,
    operators: Arc<TetElasticAssembly>,
    state: DynamicState,
    thermal: Vec<f64>,
    thermal_strains: Option<response::FrozenThermal3>,
    force: Vec<f64>,
    full_force: Vec<(usize, f64)>,
    config: TransientConfig3,
    prescribed: Option<motion::PrescribedMotion3>,
}
impl<'g, 'm> TetLoadCase3<'g, 'm> {
    /// Assemble K and consistent M ONCE and solve M a0=f0-C v0-K u0.
    /// Initial displacement/velocity are full node-major xyz SI vectors; fixed
    /// components must be exactly zero. Thermal strain is optional and frozen.
    /// Materials remain borrowed for density-aware history evaluation.
    pub fn begin_transient<'p, 's>(&self, problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>, time_s: f64,
        displacement_m: &[[f64; 3]], velocity_m_s: &[[f64; 3]], config: TransientConfig3,
        cx: &Cx<'_>) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientError3> {
        self.begin_transient_impl(problem, thermal, time_s, displacement_m, velocity_m_s,
            None, None, config, cx)
    }

    #[allow(clippy::too_many_arguments)]
    fn begin_transient_impl<'p, 's>(&self, problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>, time_s: f64,
        displacement_m: &[[f64; 3]], velocity_m_s: &[[f64; 3]],
        prescribed_acceleration_m_s2: Option<&[[f64; 3]]>,
        initial_reduced_momentum: Option<&[f64]>, config: TransientConfig3,
        cx: &Cx<'_>) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientError3> {
        check(cx)?; config.validate()?; number(time_s)?;
        if initial_reduced_momentum.is_some() && prescribed_acceleration_m_s2.is_some() {
            return Err(TransientError3::Invalid("material reconfiguration requires stationary zero supports"));
        }
        if problem.nodes_m != self.geometry.nodes || problem.tetrahedra != self.geometry.tets {
            return Err(TransientError3::Invalid("initial load belongs to a different reference mesh"));
        }
        let full = problem.nodes_m.len().checked_mul(3).ok_or(TransientError3::Budget)?;
        if full.saturating_sub(problem.fixed_dofs.len()) > config.maximum_free_dofs
            || problem.tetrahedra.len().checked_mul(144)
                .is_none_or(|n| n > config.maximum_element_entries) { return Err(TransientError3::Budget); }
        if displacement_m.len() != problem.nodes_m.len() || velocity_m_s.len() != problem.nodes_m.len() {
            return Err(TransientError3::Invalid("initial motion needs one xyz vector per node"));
        }
        for (u, v) in displacement_m.iter().zip(velocity_m_s) {
            check(cx)?;
            if !u.iter().chain(v).all(|x| x.is_finite()) || norm(*u) > config.maximum_displacement_m {
                return Err(TransientError3::Invalid("nonfinite or excessive initial motion"));
            }
        }
        let mut request = TetLinearElasticProblem { nodes_m: problem.nodes_m,
            tetrahedra: problem.tetrahedra, materials: problem.materials,
            fixed_dofs: problem.fixed_dofs, budget: problem.budget };
        request.budget.maximum_free_dofs = request.budget.maximum_free_dofs.min(config.maximum_free_dofs);
        // A fully clamped solid still stores thermal stress energy. Reuse the
        // native zero-free-DOF assembly, which keeps all geometry/mass checks.
        let operators = Arc::new(request.assemble_prescribed_dynamics(cx)?);
        for &d in problem.fixed_dofs {
            check(cx)?;
            if prescribed_acceleration_m_s2.is_none()
                && (displacement_m[d/3][d%3] != 0.0 || velocity_m_s[d/3][d%3] != 0.0) {
                return Err(TransientError3::Invalid("transient component supports must be stationary and zero"));
            }
        }
        for matrix in [&operators.mass, &operators.stiffness] {
            for row in 0..matrix.nrows() {
                check(cx)?;
                if matrix.row(row).1.iter().any(|v| !v.is_finite()) {
                    return Err(TransientError3::Invalid("nonfinite assembled mass or stiffness"));
                }
            }
        }
        let n = operators.free_dofs.len();
        let thermal_field = thermal;
        let (thermal, full_thermal) = match thermal {
            Some(field) => {
                let load = request.assemble_thermal_load(field, cx)?;
                (load.reduced_force_n, load.full_equivalent_force_n)
            }
            None => (vec![0.0; n], Vec::new()),
        };
        let prescribed = prescribed_acceleration_m_s2.map(|a|
            motion::PrescribedMotion3::initial(problem, &operators.free_dofs,
                displacement_m, velocity_m_s, a, full_thermal, cx)).transpose()?;
        let thermal_strains = thermal_field.map(|field| response::FrozenThermal3::capture(field, cx)).transpose()?;
        let q: Vec<_> = operators.free_dofs.iter().map(|&d| displacement_m[d/3][d%3]).collect();
        let v: Vec<_> = operators.free_dofs.iter().map(|&d| velocity_m_s[d/3][d%3]).collect();
        let mut run = TetTransientRun3 { geometry: self.geometry, problem, operators,
            state: DynamicState { time_s, step: 0, displacement: q, velocity: v,
                acceleration: vec![0.0; n], kinetic_energy_j: 0.0, strain_energy_j: 0.0,
                total_energy_j: 0.0 }, thermal, thermal_strains, force: vec![0.0; n],
            full_force: self.nodal_forces_n().to_vec(), config, prescribed };
        if let Some(momentum) = initial_reduced_momentum {
            if momentum.len() != n || momentum.iter().any(|v| !v.is_finite()) {
                return Err(TransientError3::Invalid("one finite canonical momentum per free DOF required"));
            }
            let current = run.apply(1.0, 0.0, &run.state.velocity, cx)?;
            // Identical mass preserves the accepted velocity bit-for-bit.
            if current.iter().zip(momentum).any(|(a,b)| a != b) {
                let scale = number(length(momentum))?;
                run.state.velocity = if scale == 0.0 { vec![0.0; n] } else {
                    let normalized: Vec<_> = momentum.iter().map(|v| v/scale).collect();
                    // Solve a unit-norm equation: the primal force tolerance
                    // has no units appropriate to a momentum conservation gate.
                    let controls = TransientConfig3 {
                        absolute_force_tolerance_n: f64::MIN_POSITIVE, ..config };
                    let mut velocity = match solve(&run.operator(1.0, 0.0, cx), &normalized, controls) {
                        Ok((velocity, _, _)) => velocity,
                        Err(TransientError3::Solve { residual_n, allowance_n, .. }) => {
                            return Err(TransientError3::Momentum {
                                residual_kg_m_s: number(scale*residual_n)?,
                                allowance_kg_m_s: number(scale*allowance_n)? });
                        }
                        Err(error) => return Err(error),
                    };
                    for value in &mut velocity { check(cx)?; *value = number(*value*scale)?; }
                    velocity
                };
            }
        }
        run.force = run.reduce_force(self, cx)?;
        let ku = run.apply(0.0, 1.0, &run.state.displacement, cx)?;
        let cv = run.apply(config.mass_damping_per_s, config.stiffness_damping_s, &run.state.velocity, cx)?;
        let mut rhs: Vec<_> = (0..n).map(|i| run.force[i] + run.thermal[i] - ku[i] - cv[i]).collect();
        if let Some(supports) = &run.prescribed { run.subtract_support_force(supports, &mut rhs, cx)?; }
        let op = run.operator(1.0, 0.0, cx);
        run.state.acceleration = solve(&op, &rhs, config)?.0;
        let (kinetic, strain) = run.energy(&run.state.displacement, &run.state.velocity, cx)?;
        run.state.kinetic_energy_j = kinetic; run.state.strain_energy_j = strain;
        run.state.total_energy_j = number(kinetic + strain)?;
        run.initialize_support_energy(cx)?;
        check(cx)?; Ok(run)
    }
}
impl TetTransientRun3<'_, '_, '_, '_> {
    /// Reduced state, ordered by `free_dofs`; no mutable checkpoint escape.
    pub fn state(&self) -> &DynamicState { &self.state }
    pub fn free_dofs(&self) -> &[usize] { &self.operators.free_dofs }
    /// The admitted bound on full nodal displacement, including support motion.
    pub fn maximum_displacement_m(&self) -> f64 { self.config.maximum_displacement_m }
    /// Exact material/source problem retained by the accepted trajectory.
    pub fn material_problem(&self) -> &TetLinearElasticProblem<'_> { self.problem }
    /// Expand one accepted state to full node-major xyz displacement/velocity/acceleration.
    pub fn nodal_motion(&self, cx: &Cx<'_>) -> Result<[Vec<[f64; 3]>; 3], TransientError3> {
        check(cx)?;
        let mut out = std::array::from_fn(|_| vec![[0.0; 3]; self.geometry.nodes.len()]);
        for (i, &d) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 { check(cx)?; }
            out[0][d/3][d%3] = self.state.displacement[i];
            out[1][d/3][d%3] = self.state.velocity[i];
            out[2][d/3][d%3] = self.state.acceleration[i];
        }
        if let Some(supports) = &self.prescribed {
            for &d in self.problem.fixed_dofs {
                check(cx)?;
                for field in 0..3 { out[field][d/3][d%3] = supports.motion[field][d]; }
            }
        }
        check(cx)?; Ok(out)
    }
    fn operator<'a, 'cx>(&'a self, cm: f64, ck: f64, cx: &'a Cx<'cx>) -> Operator<'a, 'cx> {
        Operator { m: &self.operators.mass, k: &self.operators.stiffness, cm, ck, cx }
    }
    fn apply(&self, cm: f64, ck: f64, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, TransientError3> {
        let mut out = vec![0.0; x.len()]; self.operator(cm, ck, cx).apply(x, &mut out); check(cx)?;
        if out.iter().any(|v| !v.is_finite()) { return Err(TransientError3::Invalid("nonfinite operator application")); }
        Ok(out)
    }
    fn energy(&self, q: &[f64], v: &[f64], cx: &Cx<'_>) -> Result<(f64, f64), TransientError3> {
        Ok((number(0.5 * inner(v, &self.apply(1.0, 0.0, v, cx)?))?,
            number(0.5 * inner(q, &self.apply(0.0, 1.0, q, cx)?))?))
    }
    fn reduce_force(&self, load: &TetLoadCase3<'_, '_>, cx: &Cx<'_>) -> Result<Vec<f64>, TransientError3> {
        if load.geometry.nodes != self.geometry.nodes || load.geometry.tets != self.geometry.tets {
            return Err(TransientError3::Invalid("transient load belongs to another reference mesh"));
        }
        let mut out = vec![0.0; self.free_dofs().len()];
        for &(d, f) in load.nodal_forces_n() {
            check(cx)?;
            if let Ok(i) = self.free_dofs().binary_search(&d) { out[i] = f; }
        }
        Ok(out)
    }
    /// Advance one absolute force endpoint. M,C,K stay frozen; loads are sampled
    /// at t+h, NOT at an alpha stage or a wrapped/clamped frame clock. Newmark
    /// integrates the piecewise-linear endpoint force; callers must resolve
    /// discontinuities and short pulses explicitly. No hidden event detection.
    pub fn advance(&mut self, load_next: &TetLoadCase3<'_, '_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientStep3, TransientError3> {
        check(cx)?;
        if self.prescribed.is_some() {
            return Err(TransientError3::Invalid("prescribed-support runs require an explicit next support endpoint"));
        }
        self.advance_impl(load_next, None, h_s, cx)
    }

    fn advance_impl(&mut self, load_next: &TetLoadCase3<'_, '_>,
        next_supports: Option<motion::PrescribedMotion3>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientStep3, TransientError3> {
        self.advance_with_thermal_impl(load_next, next_supports, None, None, h_s, cx)
    }

    fn advance_with_thermal_impl(&mut self, load_next: &TetLoadCase3<'_, '_>,
        mut next_supports: Option<motion::PrescribedMotion3>, thermal_next: Option<&[f64]>,
        stiffness_before: Option<&Csr>, h_s: f64, cx: &Cx<'_>) -> Result<TransientStep3, TransientError3> {
        check(cx)?;
        if thermal_next.is_some_and(|field| field.len() != self.thermal.len())
            || (thermal_next.is_some() && next_supports.is_some()) {
            return Err(TransientError3::Invalid("thermal endpoint requires matching fields and stationary zero supports"));
        }
        if stiffness_before.is_some_and(|matrix| matrix.nrows() != self.free_dofs().len())
            || (stiffness_before.is_some() && next_supports.is_some()) {
            return Err(TransientError3::Invalid("changing stiffness requires fixed mass and stationary zero supports"));
        }
        if self.state.step >= self.config.maximum_steps { return Err(TransientError3::Budget); }
        let next_time = self.state.time_s + h_s;
        if !h_s.is_finite() || h_s <= 0.0 || !next_time.is_finite() || next_time <= self.state.time_s {
            return Err(TransientError3::Invalid("step must reach a distinct finite later time"));
        }
        let b = number(0.25 * h_s * h_s)?;
        if b == 0.0 { return Err(TransientError3::Invalid("squared time step underflows")); }
        let force = self.reduce_force(load_next, cx)?;
        let n = force.len();
        let eta = self.config.mass_damping_per_s; let zeta = self.config.stiffness_damping_s;
        let thermal_force = thermal_next.unwrap_or(&self.thermal);
        let changing_damping = stiffness_before.filter(|_| zeta != 0.0);
        let (q, v, a, iterations, force_residual_n) = if let Some(previous) = changing_damping {
            let endpoint = self.damped_material_endpoint(previous, &force, thermal_force, h_s, b, cx)?;
            (endpoint.displacement, endpoint.velocity, endpoint.acceleration,
                endpoint.iterations, endpoint.force_residual_n)
        } else {
            // Preserve the original arithmetic for all previously admitted runs.
            let mut qp = Vec::with_capacity(n); let mut vp = Vec::with_capacity(n);
            for i in 0..n {
                if i % 256 == 0 { check(cx)?; }
                qp.push(number(self.state.displacement[i] + h_s*self.state.velocity[i] + b*self.state.acceleration[i])?);
                vp.push(number(self.state.velocity[i] + 0.5*h_s*self.state.acceleration[i])?);
            }
            let ku = self.apply(0.0, 1.0, &qp, cx)?; let cv = self.apply(eta, zeta, &vp, cx)?;
            let mut rhs: Vec<_> = (0..n).map(|i| force[i]+thermal_force[i]-ku[i]-cv[i]).collect();
            if let Some(supports) = &next_supports { self.subtract_support_force(supports, &mut rhs, cx)?; }
            let op = self.operator(number(1.0+0.5*h_s*eta)?, number(b+0.5*h_s*zeta)?, cx);
            let (a, iterations, force_residual_n) = solve(&op, &rhs, self.config)?;
            let mut q = Vec::with_capacity(n); let mut v = Vec::with_capacity(n);
            for i in 0..n {
                if i % 256 == 0 { check(cx)?; }
                q.push(number(qp[i]+b*a[i])?); v.push(number(vp[i]+0.5*h_s*a[i])?);
            }
            (q, v, a, iterations, force_residual_n)
        };
        let mut dq = Vec::with_capacity(n); let mut average_v = Vec::with_capacity(n);
        for i in 0..n {
            if i % 256 == 0 { check(cx)?; }
            dq.push(number(q[i]-self.state.displacement[i])?);
            average_v.push(0.5*v[i]+0.5*self.state.velocity[i]);
        }
        // A nodal bound applies to the vector norm, not individual components.
        let mut node = usize::MAX; let mut magnitude = 0.0_f64;
        for (&d, &u) in self.free_dofs().iter().zip(&q) {
            if d/3 != node { node = d/3; magnitude = 0.0; }
            magnitude = magnitude.hypot(u);
            if magnitude > self.config.maximum_displacement_m {
                return Err(TransientError3::Invalid("candidate exceeds displacement allowance"));
            }
        }
        let (kinetic, strain, mechanical_work_j, support_work_j,
            thermal_equivalent_work_j, damping_loss_j) = if let Some(supports) = &mut next_supports {
            self.supported_accounting(supports, [&q, &v, &a], load_next, h_s, cx)?
        } else {
            let (kinetic, strain) = self.energy(&q, &v, cx)?;
            let average_force: Vec<_> = force.iter().zip(&self.force).map(|(a,b)| 0.5*a+0.5*b).collect();
            let thermal_work = if let Some(next) = thermal_next {
                let average: Vec<_> = self.thermal.iter().zip(next).map(|(a,b)| 0.5*a+0.5*b).collect();
                number(inner(&average, &dq))?
            } else { number(inner(&self.thermal, &dq))? };
            let damping_loss = if let Some(previous) = changing_damping {
                self.material_damping_loss(previous, &average_v, h_s, cx)?
            } else {
                number(h_s*inner(&average_v, &self.apply(eta,zeta,&average_v,cx)?))?
            };
            (kinetic, strain, number(inner(&average_force, &dq))?, 0.0, thermal_work, damping_loss)
        };
        let stiffness_parameter_work_j = if let Some(previous) = stiffness_before {
            let mut old_action = vec![0.0; n];
            Operator { m: &self.operators.mass, k: previous, cm: 0.0, ck: 1.0, cx }
                .apply(&q, &mut old_action);
            check(cx)?;
            let new_action = self.apply(0.0, 1.0, &q, cx)?;
            let mut work = 0.0;
            for i in 0..n {
                if i % 256 == 0 { check(cx)?; }
                work = number(work + 0.5*self.state.displacement[i]
                    * number(new_action[i]-old_action[i])?)?;
            }
            work
        } else { 0.0 };
        let quadratic_energy_j = number(kinetic+strain)?;
        let balance_error_j = number(quadratic_energy_j - self.state.total_energy_j
            + damping_loss_j - mechanical_work_j - thermal_equivalent_work_j - support_work_j
            - stiffness_parameter_work_j)?;
        let scale = self.state.total_energy_j.abs().max(quadratic_energy_j.abs())
            .max(mechanical_work_j.abs()).max(thermal_equivalent_work_j.abs()).max(damping_loss_j.abs())
            .max(support_work_j.abs()).max(stiffness_parameter_work_j.abs());
        let allowance = number(self.config.absolute_energy_tolerance_j + self.config.relative_energy_tolerance*scale)?;
        if balance_error_j.abs() > allowance || kinetic < -allowance || strain < -allowance || damping_loss_j < -allowance {
            return Err(TransientError3::Energy { defect_j: balance_error_j, allowance_j: allowance });
        }
        let report = TransientStep3 { time_s: next_time, step_size_s: h_s, iterations, force_residual_n,
            mechanical_work_j, support_work_j, thermal_equivalent_work_j, stiffness_parameter_work_j,
            damping_loss_j, quadratic_energy_j, balance_error_j };
        let full_force = load_next.nodal_forces_n().to_vec();
        let thermal = thermal_next.map(<[f64]>::to_vec);
        check(cx)?;
        self.state = DynamicState { time_s: next_time, step: self.state.step.checked_add(1).ok_or(TransientError3::Budget)?,
            displacement: q, velocity: v, acceleration: a, kinetic_energy_j: kinetic,
            strain_energy_j: strain, total_energy_j: quadratic_energy_j };
        self.force = force;
        self.full_force = full_force;
        self.prescribed = next_supports;
        if let Some(thermal) = thermal { self.thermal = thermal; }
        Ok(report)
    }
}
