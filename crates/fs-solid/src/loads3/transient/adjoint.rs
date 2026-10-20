//! Discrete load/state adjoints of the admitted acceleration-form Newmark map.
//!
//! For b=h^2/4 and d=h/2, A=M+d*C+b*K and A*a1=f1+fT-K*qp-C*vp.
//! Solve A^T*lambda = b*q_bar+d*v_bar+a_bar, then reverse the original
//! predictors. The same mass solve reverses initialization and finite jumps.
//! These are derivatives of the discrete equations at FIXED h, M, C, K and
//! restraints, not derivatives of PCG iterations, acceptance tests or roundoff.
//! The load/state map is linear, so its pullback needs no stored primal states.
//! Material-scale contractions are supplied by `material`. Geometry, event-time
//! and direct stress/reaction/energy-observation derivatives are not supplied.

use super::*;

/// Independent transposed-solve and memory limits. Cotangents have objective-
/// dependent units; no primal newton-valued absolute tolerance is reused.
#[derive(Debug, Clone, Copy)]
pub struct TransientAdjointConfig3 {
    pub solve: TetStaticSolveConfig,
    pub maximum_free_dofs: usize,
    /// Combined retained nonzero entries of M^T and K^T.
    pub maximum_transpose_entries: usize,
}
impl Default for TransientAdjointConfig3 {
    fn default() -> Self {
        Self { solve: TetStaticSolveConfig { relative_residual_tolerance: 1e-11,
            maximum_iterations: 40_000 }, maximum_free_dofs: 50_000,
            maximum_transpose_entries: 28_800_000 }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum TransientAdjointError3 {
    Invalid(&'static str),
    Budget,
    /// Residual of the normalized cotangent equation, not a force in newtons.
    Solve { iterations: usize, relative_residual: f64, tolerance: f64 },
    Dynamics(TransientError3),
    History(crate::loads3::history::LoadHistoryError3),
    Cancelled,
}
impl std::fmt::Display for TransientAdjointError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "tetrahedral transient adjoint refused: {self:?}")
    }
}
impl std::error::Error for TransientAdjointError3 {}
impl From<TransientError3> for TransientAdjointError3 {
    fn from(e: TransientError3) -> Self {
        if e == TransientError3::Cancelled { Self::Cancelled } else { Self::Dynamics(e) }
    }
}

/// Derivatives of a scalar objective with respect to reduced q, v and a.
/// Order is `TransientAdjoint3::free_dofs`; units are objective/m,
/// objective/(m/s), and objective/(m/s^2), NOT displacement/velocity/acceleration.
#[derive(Debug, Clone, PartialEq)]
pub struct MotionCotangent3 {
    pub q: Vec<f64>,
    pub v: Vec<f64>,
    pub a: Vec<f64>,
}
impl MotionCotangent3 {
    fn validate(&self, n: usize, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
        for values in [&self.q, &self.v, &self.a] { validate_values(values, n, cx)?; }
        Ok(())
    }
}

/// Pullback of one interval or acceleration reset. `force` is also the
/// derivative with respect to an independently varied equivalent thermal load.
/// It excludes direct objective dependence on forces/reactions/energy reports.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientPullback3 {
    pub previous: MotionCotangent3,
    /// One derivative per free force component, objective/N.
    pub force: Vec<f64>,
    pub iterations: usize,
    pub true_relative_residual: f64,
}

#[derive(Debug)]
struct Transposes3 { mass: Csr, stiffness: Csr }

/// Frozen, shared transposed operators for the exact run that created them.
/// Constructing this optional object leaves the original trajectory unchanged.
#[derive(Debug, Clone)]
pub struct TransientAdjoint3<'g, 'm> {
    geometry: &'g TetLoadGeometry3<'m>,
    operators: Arc<TetElasticAssembly>,
    transposes: Arc<Transposes3>,
    primal_config: TransientConfig3,
    config: TransientAdjointConfig3,
}

fn validate_values(values: &[f64], n: usize, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
    check(cx)?;
    if values.len() != n { return Err(TransientAdjointError3::Invalid("cotangent dimension differs from free DOFs")); }
    for (i, &value) in values.iter().enumerate() {
        if i % 256 == 0 { check(cx)?; }
        if !value.is_finite() { return Err(TransientAdjointError3::Invalid("nonfinite cotangent")); }
    }
    Ok(())
}

impl<'g, 'm> TetTransientRun3<'g, 'm, '_, '_> {
    /// Prepare actual sparse transposes once. Do not assume bitwise symmetry of
    /// independently rounded element entries. Sparse transposition remains its
    /// existing bounded-by-entry-count indivisible owner operation; applications
    /// and PCG retain the original inner cancellation checkpoints.
    pub fn adjoint(&self, config: TransientAdjointConfig3, cx: &Cx<'_>)
        -> Result<TransientAdjoint3<'g, 'm>, TransientAdjointError3> {
        check(cx)?;
        if self.prescribed.is_some() {
            return Err(TransientAdjointError3::Invalid("reduced-only adjoint cannot reverse supports; use prescribed_adjoint"));
        }
        self.prepare_adjoint(config, cx)
    }

    // Shared sparse owner. The public reduced-only entry point retains its
    // homogeneous-support contract; `support` adds the omitted cross terms.
    fn prepare_adjoint(&self, config: TransientAdjointConfig3, cx: &Cx<'_>)
        -> Result<TransientAdjoint3<'g, 'm>, TransientAdjointError3> {
        check(cx)?;
        let tol = config.solve.relative_residual_tolerance;
        if !tol.is_finite() || tol <= 0.0 || tol >= 1.0 || config.solve.maximum_iterations == 0 {
            return Err(TransientAdjointError3::Invalid("positive iteration limit and relative tolerance in (0,1) required"));
        }
        if self.free_dofs().len() > config.maximum_free_dofs { return Err(TransientAdjointError3::Budget); }
        let mut entries = 0usize;
        for matrix in [&self.operators.mass, &self.operators.stiffness] {
            for row in 0..matrix.nrows() {
                check(cx)?;
                entries = entries.checked_add(matrix.row(row).1.len()).ok_or(TransientAdjointError3::Budget)?;
                if entries > config.maximum_transpose_entries { return Err(TransientAdjointError3::Budget); }
            }
        }
        let mass = fs_sparse::ops::transpose(&self.operators.mass);
        check(cx)?;
        let stiffness = fs_sparse::ops::transpose(&self.operators.stiffness);
        check(cx)?;
        Ok(TransientAdjoint3 { geometry: self.geometry, operators: self.operators.clone(),
            transposes: Arc::new(Transposes3 { mass, stiffness }), primal_config: self.config, config })
    }
}

impl TransientAdjoint3<'_, '_> {
    pub fn free_dofs(&self) -> &[usize] { &self.operators.free_dofs }
    fn operator<'a, 'cx>(&'a self, cm: f64, ck: f64, cx: &'a Cx<'cx>) -> Operator<'a, 'cx> {
        Operator { m: &self.transposes.mass, k: &self.transposes.stiffness, cm, ck, cx }
    }
    fn apply(&self, cm: f64, ck: f64, x: &[f64], cx: &Cx<'_>)
        -> Result<Vec<f64>, TransientAdjointError3> {
        let mut out = vec![0.0; x.len()];
        self.operator(cm, ck, cx).apply(x, &mut out);
        validate_values(&out, x.len(), cx)?;
        Ok(out)
    }
    fn solve(&self, cm: f64, ck: f64, rhs: &[f64], cx: &Cx<'_>)
        -> Result<(Vec<f64>, usize, f64), TransientAdjointError3> {
        check(cx)?;
        let scale = number(length(rhs))?;
        if scale == 0.0 { return Ok((vec![0.0; rhs.len()], 0, 0.0)); }
        let normalized: Vec<_> = rhs.iter().map(|v| v/scale).collect();
        // Reuse native PCG on a unit-norm RHS. This makes its stopping rule
        // independent of objective scaling, including very small seeds.
        let config = TransientConfig3 { solve: self.config.solve,
            absolute_force_tolerance_n: f64::MIN_POSITIVE, ..self.primal_config };
        let op = self.operator(cm, ck, cx);
        let (mut lambda, iterations, residual) = match super::solve(&op, &normalized, config) {
            Ok(result) => result,
            Err(TransientError3::Solve { iterations, residual_n, .. }) => {
                return Err(TransientAdjointError3::Solve { iterations,
                    relative_residual: residual_n/length(&normalized),
                    tolerance: config.solve.relative_residual_tolerance });
            }
            Err(error) => return Err(error.into()),
        };
        let relative = number(residual/length(&normalized))?;
        if relative > config.solve.relative_residual_tolerance {
            return Err(TransientAdjointError3::Solve { iterations, relative_residual: relative,
                tolerance: config.solve.relative_residual_tolerance });
        }
        for (i, value) in lambda.iter_mut().enumerate() {
            if i % 256 == 0 { check(cx)?; }
            *value = number(*value * scale)?;
        }
        check(cx)?;
        Ok((lambda, iterations, relative))
    }

    /// Reverse one accepted smooth interval at fixed positive h. Initial q/v/a
    /// are independent here; chain `pullback_acceleration` at initialization to
    /// account for the actual initial mass solve. No primal trajectory is run.
    pub fn pullback_step(&self, h_s: f64, seed: &MotionCotangent3, cx: &Cx<'_>)
        -> Result<TransientPullback3, TransientAdjointError3> {
        seed.validate(self.free_dofs().len(), cx)?;
        if !h_s.is_finite() || h_s <= 0.0 { return Err(TransientAdjointError3::Invalid("positive finite adjoint step required")); }
        let b = number(0.25*h_s*h_s)?; let d = 0.5*h_s;
        if b == 0.0 { return Err(TransientAdjointError3::Invalid("squared adjoint step underflows")); }
        let mut rhs = Vec::with_capacity(seed.q.len());
        for i in 0..seed.q.len() {
            if i % 256 == 0 { check(cx)?; }
            rhs.push(number(b*seed.q[i] + d*seed.v[i] + seed.a[i])?);
        }
        let eta = self.primal_config.mass_damping_per_s;
        let zeta = self.primal_config.stiffness_damping_s;
        let (force, iterations, true_relative_residual) =
            self.solve(number(1.0+d*eta)?, number(b+d*zeta)?, &rhs, cx)?;
        let kt = self.apply(0.0, 1.0, &force, cx)?;
        let ct = self.apply(eta, zeta, &force, cx)?;
        let mut previous = MotionCotangent3 { q: Vec::with_capacity(rhs.len()),
            v: Vec::with_capacity(rhs.len()), a: Vec::with_capacity(rhs.len()) };
        for i in 0..rhs.len() {
            if i % 256 == 0 { check(cx)?; }
            let q = number(seed.q[i]-kt[i])?;
            let v = number(seed.v[i]-ct[i])?;
            previous.q.push(q);
            previous.v.push(number(h_s*q+v)?);
            previous.a.push(number(b*q+d*v)?);
        }
        check(cx)?;
        Ok(TransientPullback3 { previous, force, iterations, true_relative_residual })
    }

    /// Reverse a=M^-1(f+fT-C*v-K*q). Used for initialization AND a finite load
    /// jump, which preserves q/v and overwrites the previous acceleration.
    /// Consequently the returned previous acceleration cotangent is exactly zero.
    pub fn pullback_acceleration(&self, seed: &MotionCotangent3, cx: &Cx<'_>)
        -> Result<TransientPullback3, TransientAdjointError3> {
        seed.validate(self.free_dofs().len(), cx)?;
        let (force, iterations, true_relative_residual) = self.solve(1.0, 0.0, &seed.a, cx)?;
        let kt = self.apply(0.0, 1.0, &force, cx)?;
        let ct = self.apply(self.primal_config.mass_damping_per_s,
            self.primal_config.stiffness_damping_s, &force, cx)?;
        let mut previous = MotionCotangent3 { q: Vec::with_capacity(force.len()),
            v: Vec::with_capacity(force.len()), a: vec![0.0; force.len()] };
        for i in 0..force.len() {
            if i % 256 == 0 { check(cx)?; }
            previous.q.push(number(seed.q[i]-kt[i])?);
            previous.v.push(number(seed.v[i]-ct[i])?);
        }
        check(cx)?;
        Ok(TransientPullback3 { previous, force, iterations, true_relative_residual })
    }

    /// Chain a reduced force cotangent through an already integrated load
    /// direction, e.g. unit pressure/traction or unit acceleration at the same
    /// fixed material. Constrained-node forces have zero effect on motion, not
    /// zero effect on an objective with direct reaction dependence.
    pub fn contract_load_direction(&self, force: &[f64], direction: &TetLoadCase3<'_, '_>, cx: &Cx<'_>)
        -> Result<f64, TransientAdjointError3> {
        validate_values(force, self.free_dofs().len(), cx)?;
        if direction.geometry.nodes != self.geometry.nodes || direction.geometry.tets != self.geometry.tets {
            return Err(TransientAdjointError3::Invalid("load direction belongs to another reference mesh"));
        }
        let mut result = 0.0;
        for &(d, value) in direction.nodal_forces_n() {
            check(cx)?;
            if let Ok(i) = self.free_dofs().binary_search(&d) { result = number(value.mul_add(force[i], result))?; }
        }
        check(cx)?; Ok(result)
    }
}

/// End-to-end load-control gradients on a bounded recorded trajectory.
pub mod trajectory;

impl From<crate::loads3::history::LoadHistoryError3> for TransientAdjointError3 {
    fn from(error: crate::loads3::history::LoadHistoryError3) -> Self {
        if error == crate::loads3::history::LoadHistoryError3::Cancelled { Self::Cancelled }
        else { Self::History(error) }
    }
}

/// Endpoint stiffness, density and Rayleigh-damping derivatives.
pub mod material;

/// Full nodal state and prescribed-support control pullbacks.
pub mod support;
