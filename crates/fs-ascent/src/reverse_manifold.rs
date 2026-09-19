//! Reverse-mode L-BFGS on heterogeneous products of the live manifold types.
//! Geometry and line search delegate to fs-opt and the existing Wolfe engine.

use crate::runner::Packing;
use crate::stop::{StopObservation, StopReason, StopRule};
use crate::wolfe::try_strong_wolfe_with_budget;
use fs_exec::Cx;
use fs_opt::reverse::ReverseError;
use fs_opt::{
    Manifold, OptError, ProductDifferentialError, ProductFactor, ProductFactorId,
    ProductManifold, ProductManifoldError, ReverseProblem, ReverseProblemError,
};
use std::collections::VecDeque;

type Pair = (Vec<f64>, Vec<f64>, f64);

/// Structural, geometry, numerical or cancellation refusal. No panic adapter.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseManifoldError {
    /// This driver does not solve additional equality/inequality constraints.
    ConstraintsUnsupported,
    /// Invalid stopping rule or initial packed length.
    InvalidInput(&'static str),
    /// Original factor/shape/allocation refusal.
    Geometry(ProductDifferentialError),
    /// Original reverse-tape refusal.
    Evaluation(ReverseProblemError),
    /// Unrepresentable numerical work, not convergence.
    Numerical(&'static str),
    /// Cancellation; accepted state remains usable and spent work is retained.
    Cancelled,
}

impl From<ProductDifferentialError> for ReverseManifoldError {
    fn from(error: ProductDifferentialError) -> Self { Self::Geometry(error) }
}
impl From<ProductManifoldError> for ReverseManifoldError {
    fn from(error: ProductManifoldError) -> Self { Self::Geometry(error.into()) }
}
impl From<ReverseProblemError> for ReverseManifoldError {
    fn from(error: ReverseProblemError) -> Self {
        match error {
            ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::Cancelled)) => Self::Cancelled,
            other => Self::Evaluation(other),
        }
    }
}
impl core::fmt::Display for ReverseManifoldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ConstraintsUnsupported => write!(f, "product-manifold L-BFGS cannot drop declared constraints"),
            Self::InvalidInput(what) | Self::Numerical(what) => write!(f, "{what}"),
            Self::Geometry(error) => write!(f, "{error}"),
            Self::Evaluation(error) => write!(f, "{error}"),
            Self::Cancelled => write!(f, "product-manifold study cancelled at its accepted checkpoint"),
        }
    }
}
impl std::error::Error for ReverseManifoldError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Geometry(e) => Some(e), Self::Evaluation(e) => Some(e), _ => None }
    }
}

/// Cached endpoint diagnostics; obtaining these never evaluates the problem.
#[derive(Debug, Clone)]
pub struct ReverseManifoldReport {
    /// Stop attribution; budget has priority over simultaneous convergence.
    pub reason: StopReason,
    /// Signed, weighted objective at the last accepted point.
    pub f: f64,
    /// Infinity norm in the concatenated retraction-parameter coordinates.
    pub grad_norm: f64,
    /// Euclidean norm in those same parameter coordinates, not physical units.
    pub grad_l2_norm: f64,
    /// Cumulative accepted steps.
    pub iters: usize,
    /// Cumulative sample attempts, including rejected geometry and arithmetic.
    pub evals: usize,
    /// Cumulative unavailable-domain trials.
    pub rejected_trials: usize,
    /// Accepted steps that restarted memory after a numerical transport refusal.
    pub memory_restarts: usize,
}

fn poll(cx: Option<&Cx<'_>>) -> Result<(), ReverseManifoldError> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| ReverseManifoldError::Cancelled)?; }
    Ok(())
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    x.iter().zip(y).map(|(x, y)| x*y).sum()
}
fn inf(x: &[f64]) -> f64 { x.iter().map(|x| x.abs()).fold(0.0, f64::max) }
fn norm(x: &[f64]) -> f64 {
    let scale = inf(x);
    if scale == 0.0 { return 0.0; }
    scale * fs_math::det::sqrt(x.iter().map(|x| (x/scale)*(x/scale)).sum())
}
fn rho(s: &[f64], y: &[f64]) -> Option<f64> {
    let (sn, yn, sy) = (norm(s), norm(y), dot(s, y));
    if !sn.is_finite() || !yn.is_finite() || sn == 0.0 || yn == 0.0 || !sy.is_finite() || sy <= 0.0 { return None; }
    let cosine: f64 = s.iter().zip(y).map(|(s, y)| (s/sn)*(y/yn)).sum();
    let reciprocal = 1.0/sy;
    (cosine.is_finite() && cosine > 1e-14 && reciprocal.is_finite()).then_some(reciprocal)
}
fn rule_budget(rule: &StopRule) -> Result<Option<usize>, ReverseManifoldError> {
    let valid = match rule {
        StopRule::GradNorm(t) => t.is_finite() && *t >= 0.0,
        StopRule::ObjectiveBelow(t) => t.is_finite(),
        StopRule::Budget(b) => return Ok(Some(*b)),
        StopRule::Stall { rel, window } => rel.is_finite() && *rel >= 0.0 && *window > 0 && *window < usize::MAX,
        StopRule::Any(children) | StopRule::All(children) => {
            if children.is_empty() { return Err(ReverseManifoldError::InvalidInput("empty stopping rule")); }
            let mut bound: Option<usize> = None;
            for child in children {
                if let Some(b) = rule_budget(child)? { bound = Some(bound.map_or(b, |old| old.min(b))); }
            }
            return Ok(bound);
        }
    };
    if valid { Ok(None) } else { Err(ReverseManifoldError::InvalidInput("invalid stopping rule")) }
}
fn recoverable(error: &ReverseManifoldError) -> bool {
    matches!(error,
        ReverseManifoldError::Numerical(_)
        | ReverseManifoldError::Evaluation(ReverseProblemError::NonFiniteObjective { .. })
        | ReverseManifoldError::Evaluation(ReverseProblemError::Reverse(ReverseError::NonFiniteAdjoint { .. }))
        | ReverseManifoldError::Evaluation(ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::EvalNonFinite { .. })))
        | ReverseManifoldError::Geometry(ProductDifferentialError::Geometry(ProductManifoldError::FactorOperation {
            source: OptError::RetractionNonFinite { .. } | OptError::RetractionDomain { .. }, ..
        }))
    )
}

/// A problem-bound checkpoint over any ordered mixture of Rn, Sphere, SO(3)
/// and Stiefel variables. The oracle and declaration order cannot change on
/// resume. Additional constraints are refused, never ignored. Zero memory
/// selects retracted steepest descent. No global-convergence claim is made.
#[derive(Debug, Clone)]
pub struct ReverseManifoldStudy<'oracle, 'problem> {
    oracle: &'oracle ReverseProblem<'problem>,
    product: ProductManifold,
    packing: Packing,
    x: Vec<f64>,
    f: f64,
    g: Vec<f64>,
    pairs: VecDeque<Pair>,
    memory: usize,
    history: Vec<f64>,
    iters: usize,
    evals: usize,
    limit: usize,
    rejected: usize,
    memory_restarts: usize,
    last_rejection: Option<ReverseManifoldError>,
}

impl<'oracle, 'problem> ReverseManifoldStudy<'oracle, 'problem> {
    /// Validate and bind the product, then spend one initial primal/gradient
    /// attempt. Only SO(3) starts adopt their canonical antipodal representative;
    /// already valid Sphere and Stiefel initial bits are not renormalized.
    pub fn new(
        oracle: &'oracle ReverseProblem<'problem>, point: &[f64], memory: usize,
        cx: Option<&Cx<'_>>,
    ) -> Result<Self, ReverseManifoldError> {
        poll(cx)?;
        if !oracle.problem().constraints().is_empty() { return Err(ReverseManifoldError::ConstraintsUnsupported); }
        let factors = oracle.problem().vars().iter().enumerate().map(|(i, v)| {
            ProductFactor::new(ProductFactorId::new(i as u32), v.manifold)
        }).collect();
        let product = ProductManifold::new(factors)?;
        let packing = Packing::new(oracle.problem());
        if point.len() != packing.dim { return Err(ReverseManifoldError::InvalidInput("initial product point length mismatch")); }
        product.validate_point(point)?;
        let mut x = point.to_vec();
        for block in product.layout().factors() {
            if matches!(block.factor().manifold(), Manifold::So3) {
                let start = block.point_offset().get() as usize;
                let end = start + 4;
                let canonical = Manifold::So3.retract(&x[start..end], &[0.0; 3])
                    .map_err(|source| ProductManifoldError::FactorOperation {
                        id: block.factor().id(), index: block.index(), source,
                    })?;
                x[start..end].copy_from_slice(&canonical);
            }
        }
        let (f, g) = sample(oracle, &packing, &product, &x, cx)?;
        poll(cx)?;
        Ok(Self {
            oracle, product, packing, x, f, g, pairs: VecDeque::new(), memory,
            history: vec![f], iters: 0, evals: 1, rejected: 0, memory_restarts: 0,
            limit: oracle.problem().budget().limit.maximum().map_or(usize::MAX,
                |b| usize::try_from(b.get()).unwrap_or(usize::MAX)),
            last_rejection: None,
        })
    }

    /// Last accepted packed point (read-only, including quaternion storage).
    #[must_use]
    pub fn point(&self) -> &[f64] { &self.x }
    /// Current gradient in retraction-parameter storage.
    #[must_use]
    pub fn gradient(&self) -> &[f64] { &self.g }
    /// Initial value followed by accepted values only.
    #[must_use]
    pub fn history(&self) -> &[f64] { &self.history }
    /// Immutable geometry associated with the problem's declaration order.
    #[must_use]
    pub fn manifold(&self) -> &ProductManifold { &self.product }
    /// Attempted samples, including failed and cancelled work.
    #[must_use]
    pub fn evaluations(&self) -> usize { self.evals }
    /// Original latest numerical trial/transport refusal, not a placeholder.
    #[must_use]
    pub fn last_rejection(&self) -> Option<&ReverseManifoldError> { self.last_rejection.as_ref() }

    fn report(&self, reason: StopReason) -> ReverseManifoldReport {
        ReverseManifoldReport { reason, f: self.f, grad_norm: inf(&self.g), grad_l2_norm: norm(&self.g),
            iters: self.iters, evals: self.evals, rejected_trials: self.rejected,
            memory_restarts: self.memory_restarts }
    }

    fn direction(&self) -> Result<Vec<f64>, ReverseManifoldError> {
        let mut q = self.g.clone();
        let mut alphas = Vec::with_capacity(self.pairs.len());
        for (s, y, rho) in self.pairs.iter().rev() {
            let a = rho*dot(s, &q);
            for (q, y) in q.iter_mut().zip(y) { *q = a.mul_add(-y, *q); }
            alphas.push(a);
        }
        if let Some((s, y, _)) = self.pairs.back() {
            let gamma = dot(s, y)/dot(y, y);
            if gamma.is_finite() && gamma > 0.0 { for q in &mut q { *q *= gamma; } }
        }
        for ((s, y, rho), a) in self.pairs.iter().zip(alphas.iter().rev()) {
            let b = rho*dot(y, &q);
            for (q, s) in q.iter_mut().zip(s) { *q = (a-b).mul_add(*s, *q); }
        }
        for q in &mut q { *q = -*q; }
        if q.iter().any(|v| !v.is_finite()) { return Err(ReverseManifoldError::Numerical("nonfinite quasi-Newton direction")); }
        // Parameter vectors are NOT ambient quaternion gradients. Only the
        // embedded Sphere/Stiefel blocks need reprojection after roundoff.
        for block in self.product.layout().factors() {
            let man = block.factor().manifold();
            if matches!(man, Manifold::Sphere { .. } | Manifold::Stiefel { .. }) {
                let p = block.point_offset().get() as usize;
                let t = block.param_offset().get() as usize;
                let n = block.manifold_layout().param_dim().get() as usize;
                let projected = man.parameter_gradient(&self.x[p..p+n], &q[t..t+n])
                    .map_err(|source| ProductManifoldError::FactorOperation {
                        id: block.factor().id(), index: block.index(), source,
                    })?;
                q[t..t+n].copy_from_slice(&projected);
            }
        }
        self.product.validate_parameter_tangent(&self.x, &q)?;
        Ok(q)
    }

    fn moved_pairs(
        &self, direction: &[f64], alpha: f64, to: &[f64], gradient: &[f64],
        restart: bool, cx: Option<&Cx<'_>>,
    ) -> Result<VecDeque<Pair>, ReverseManifoldError> {
        let mut pairs = VecDeque::new();
        if self.memory == 0 { return Ok(pairs); }
        let step: Vec<f64> = direction.iter().map(|d| alpha*d).collect();
        let transport = |v: &[f64]| -> Result<Vec<f64>, ReverseManifoldError> {
            poll(cx)?;
            Ok(self.product.transport_parameter(&self.x, &step, to, v)?)
        };
        if !restart {
            for (s, y, _) in &self.pairs {
                let (s, y) = (transport(s)?, transport(y)?);
                if let Some(r) = rho(&s, &y) { pairs.push_back((s, y, r)); }
            }
        }
        let td = transport(direction)?;
        let tg = transport(&self.g)?;
        let stretch = (norm(&td)/norm(direction)).max(1.0);
        if !stretch.is_finite() { return Err(ReverseManifoldError::Numerical("nonfinite transport stretch")); }
        let s: Vec<f64> = td.iter().map(|d| alpha*d).collect();
        // Same cautious scaling used by the single-manifold Riemannian engine.
        let y: Vec<f64> = gradient.iter().zip(tg).map(|(g, old)| g/stretch-old).collect();
        if let Some(r) = rho(&s, &y) {
            if pairs.len() == self.memory { pairs.pop_front(); }
            pairs.push_back((s, y, r));
        }
        Ok(pairs)
    }

    /// Run at most `additional_steps` accepted retractions. Every Budget leaf
    /// is a hard cumulative ceiling even inside All, alongside the problem cap.
    /// Each attempted trial costs at most one primal and one reverse sweep.
    /// Geometry failures are conservatively charged as sample attempts too.
    ///
    /// State publication follows line-search success and complete transport.
    /// Errors/cancellation keep the old accepted state but never refund work.
    /// Numerical transport failures restart optional curvature memory rather
    /// than discarding a valid Wolfe step. Mid-search resumes may repeat work;
    /// splits at completed steps preserve trajectory and accounting.
    pub fn run(
        &mut self, rule: &StopRule, additional_steps: usize, cx: Option<&Cx<'_>>,
    ) -> Result<ReverseManifoldReport, ReverseManifoldError> {
        let cap = rule_budget(rule)?.unwrap_or(usize::MAX).min(self.limit);
        let mut completed = 0usize;
        loop {
            poll(cx)?;
            if self.evals >= cap { return Ok(self.report(StopReason::Budget)); }
            let observation = StopObservation { grad_norm: inf(&self.g), objective: self.f,
                evals: self.evals, history: &self.history };
            if let Some(reason) = rule.check(&observation) { return Ok(self.report(reason)); }
            if completed == additional_steps { return Ok(self.report(StopReason::IterationCap)); }
            let (mut d, mut restart) = match self.direction() {
                Ok(d) => (d, false),
                Err(error) if recoverable(&error) => (self.g.iter().map(|g| -g).collect(), true),
                Err(error) => return Err(error),
            };
            let mut slope = dot(&self.g, &d);
            if !slope.is_finite() || slope >= 0.0 {
                d = self.g.iter().map(|g| -g).collect();
                slope = dot(&self.g, &d);
                restart = true;
            }
            if !slope.is_finite() || slope >= 0.0 {
                return Ok(self.report(if self.g.iter().all(|g| *g == 0.0) { StopReason::GradNorm } else { StopReason::Stall }));
            }
            let mut spent = 0;
            let mut rejected = 0;
            let mut last_error = None;
            let mut last = None;
            let result = {
                let mut phi = |alpha: f64| -> Result<(f64, f64), ReverseManifoldError> {
                    poll(cx)?;
                    spent += 1;
                    let trial = (|| {
                        let curve = self.product.retract_curve(&self.x, &d, alpha)?;
                        if curve.point == self.x { return Err(ReverseManifoldError::Numerical("retraction did not move point")); }
                        let (f, g) = sample(self.oracle, &self.packing, &self.product, &curve.point, cx)?;
                        let slope = dot(&g, &curve.velocity);
                        if !slope.is_finite() { return Err(ReverseManifoldError::Numerical("nonfinite curve derivative")); }
                        Ok((curve.point, f, g, slope))
                    })();
                    match trial {
                        Ok((point, f, g, slope)) => { last = Some((point, f, g)); Ok((f, slope)) }
                        Err(error) if recoverable(&error) => {
                            last = None; rejected += 1; last_error = Some(error);
                            Ok((f64::INFINITY, 0.0))
                        }
                        Err(error) => Err(error),
                    }
                };
                try_strong_wolfe_with_budget(&mut phi, self.f, slope, 1.0, 1e-4, 0.9, cap-self.evals)
            };
            self.evals += spent;
            self.rejected += rejected;
            if last_error.is_some() { self.last_rejection = last_error; }
            let outcome = result?;
            poll(cx)?;
            if !outcome.success {
                return Ok(self.report(if self.evals >= cap { StopReason::Budget } else { StopReason::Stall }));
            }
            let (point, f, g) = last.ok_or(ReverseManifoldError::Numerical("Wolfe accepted without a valid trial"))?;
            let mut transport_error = None;
            let pairs = match self.moved_pairs(&d, outcome.alpha, &point, &g, restart, cx) {
                Ok(pairs) => pairs,
                Err(error) if recoverable(&error) => { transport_error = Some(error); VecDeque::new() }
                Err(error) => return Err(error),
            };
            poll(cx)?;
            if let Some(error) = transport_error {
                self.last_rejection = Some(error); self.memory_restarts += 1;
            }
            self.x = point; self.f = f; self.g = g; self.pairs = pairs;
            self.iters += 1; self.history.push(f); completed += 1;
        }
    }
}

fn sample(
    oracle: &ReverseProblem<'_>, packing: &Packing, product: &ProductManifold,
    point: &[f64], cx: Option<&Cx<'_>>,
) -> Result<(f64, Vec<f64>), ReverseManifoldError> {
    poll(cx)?;
    let bindings = packing.unpack(point);
    let tape = match cx {
        Some(cx) => oracle.evaluate_cancellable(&bindings, cx)?,
        None => oracle.evaluate(&bindings)?,
    };
    let blocks = match cx {
        Some(cx) => tape.objective_gradient_cancellable(cx)?,
        None => tape.objective_gradient()?,
    };
    let ambient: Vec<f64> = blocks.into_iter().flatten().collect();
    let g = product.parameter_gradient(point, &ambient)?;
    poll(cx)?;
    Ok((tape.objective_value(), g))
}
