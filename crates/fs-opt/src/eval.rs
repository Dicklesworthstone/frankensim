//! Evaluation + the TOY Riemannian descent that proves the manifold
//! metadata is consumable: retractions keep iterates ON their
//! manifolds ("optimize an orientation" never becomes "optimize 9
//! numbers and renormalize when it explodes"). Gradients here are
//! finite differences through retractions — a deliberately simple
//! consumer; exact adjoints are the gradient-stack bead's.

use crate::admission::AdmissionCaps;
use crate::ir::{
    EvalLimit, Expr, Manifold, NodeId, ObjectiveEvalSite, OptError, ProbeDirection, Problem, Shape,
    VarId,
};
use fs_exec::Cx;

/// Squared-norm/dot tolerance for deciding whether a finite stored
/// point is already on its declared manifold.
const MANIFOLD_MEMBERSHIP_TOL: f64 = 1e-10;
/// Candidates below this squared norm are numerically rank-deficient;
/// normalizing them would fabricate a direction.
const RETRACTION_MIN_NORM_SQ: f64 = 1e-24;
/// Maximum traversed scalar elements between retraction/domain cancellation polls.
const RETRACTION_CHECKPOINT_STRIDE: usize = 256;

fn checkpoint_retraction_work<F>(index: usize, checkpoint: &mut F) -> Result<(), OptError>
where
    F: FnMut() -> Result<(), OptError>,
{
    if index % RETRACTION_CHECKPOINT_STRIDE == 0 {
        checkpoint()?;
    }
    Ok(())
}

/// An evaluated node value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Scalar.
    S(f64),
    /// Vector.
    V(Vec<f64>),
}

impl Value {
    /// The scalar, if it is one.
    #[must_use]
    pub fn scalar(&self) -> Option<f64> {
        match self {
            Value::S(v) => Some(*v),
            Value::V(_) => None,
        }
    }
}

fn runtime_vector_len(values: &[f64]) -> u64 {
    u64::try_from(values.len()).unwrap_or(u64::MAX)
}

fn validate_runtime_shape(problem: &Problem, node: NodeId, value: &Value) -> Result<(), OptError> {
    let expected = problem.shape(node)?;
    let matches = match (expected, value) {
        (Shape::Scalar, Value::S(_)) => true,
        (Shape::Vector(expected_len), Value::V(values)) => {
            u64::from(expected_len) == runtime_vector_len(values)
        }
        _ => false,
    };
    if matches {
        Ok(())
    } else {
        Err(OptError::EvalShape {
            node: node.0,
            expected,
            actual_vector_len: match value {
                Value::S(_) => None,
                Value::V(values) => Some(runtime_vector_len(values)),
            },
        })
    }
}

fn checked_component(node: NodeId, index: u32, values: &[f64]) -> Result<f64, OptError> {
    values
        .get(usize::try_from(index).unwrap_or(usize::MAX))
        .copied()
        .ok_or(OptError::EvalIndexOut {
            node: node.0,
            index,
            len: runtime_vector_len(values),
        })
}

/// A complete, validated runtime point for every variable in one
/// sealed [`Problem`]. Construction is keyed by [`VarId`], so caller
/// entry order carries no semantics; the retained problem reference
/// prevents applying the frame to a different graph.
#[derive(Debug)]
pub struct BindingFrame<'problem, 'value> {
    problem: &'problem Problem,
    ordered: Vec<&'value [f64]>,
}

impl<'problem, 'value> BindingFrame<'problem, 'value> {
    /// Validate and canonicalize one keyed runtime frame.
    ///
    /// # Errors
    /// [`OptError::UnknownVar`], [`OptError::BindingDuplicate`],
    /// [`OptError::BindingMissing`], [`OptError::BindingLen`],
    /// [`OptError::BindingNonFinite`], [`OptError::BindingDomain`], or
    /// [`OptError::CapExceeded`].
    pub fn new<I>(problem: &'problem Problem, bindings: I) -> Result<Self, OptError>
    where
        I: IntoIterator<Item = (VarId, &'value [f64])>,
    {
        Self::new_with_caps(problem, bindings, &AdmissionCaps::default())
    }

    fn new_with_caps<I>(
        problem: &'problem Problem,
        bindings: I,
        caps: &AdmissionCaps,
    ) -> Result<Self, OptError>
    where
        I: IntoIterator<Item = (VarId, &'value [f64])>,
    {
        validate_binding_frame_envelope(problem, caps)?;
        let mut slots = vec![None; problem.vars.len()];
        for (var, binding) in bindings {
            let slot = slots
                .get_mut(var.0 as usize)
                .ok_or(OptError::UnknownVar { id: var.0 })?;
            if slot.replace(binding).is_some() {
                return Err(OptError::BindingDuplicate { var: var.0 });
            }
        }
        let mut ordered = Vec::with_capacity(slots.len());
        for (var, binding) in slots.into_iter().enumerate() {
            ordered.push(binding.ok_or(OptError::BindingMissing { var: var as u32 })?);
        }
        validate_ordered_bindings(problem, &ordered, caps)?;
        Ok(Self { problem, ordered })
    }

    /// Evaluate an arbitrary node using this already validated frame.
    ///
    /// # Errors
    /// [`OptError::UnknownNode`], [`OptError::Unevaluable`],
    /// [`OptError::EvalShape`], [`OptError::EvalIndexOut`],
    /// [`OptError::EvalNonFinite`], or [`OptError::CapExceeded`].
    pub fn eval(&self, node: NodeId) -> Result<Value, OptError> {
        validate_eval_envelope(self.problem, node)?;
        eval_validated(self.problem, node, &self.ordered)
    }
}

/// Evaluate `node` with variable `bindings` (exactly one point per
/// declared variable, stored per its manifold and indexed by `VarId`).
/// Arbitrary subgraph roots remain evaluable, but their runtime frame
/// is still complete so an unused declaration cannot be silently
/// omitted. PDE and stochastic nodes are NOT
/// evaluable here — they refuse with a teaching error (their execution
/// belongs to FLUX/UQ runners; the IR only carries them).
///
/// # Errors
/// [`OptError::Unevaluable`] / [`OptError::UnknownVar`] /
/// [`OptError::UnknownNode`] / [`OptError::BindingCount`] /
/// [`OptError::BindingLen`] / [`OptError::BindingNonFinite`] /
/// [`OptError::BindingDomain`] / [`OptError::EvalShape`] /
/// [`OptError::EvalIndexOut`] / [`OptError::EvalNonFinite`] /
/// [`OptError::CapExceeded`].
pub fn eval(problem: &Problem, node: NodeId, bindings: &[Vec<f64>]) -> Result<Value, OptError> {
    let caps = validate_eval_envelope(problem, node)?;
    if bindings.len() != problem.vars.len() {
        return Err(OptError::BindingCount {
            vars: problem.vars.len() as u32,
            got: bindings.len() as u64,
        });
    }
    let ordered: Vec<&[f64]> = bindings.iter().map(Vec::as_slice).collect();
    validate_ordered_bindings(problem, &ordered, &caps)?;
    eval_validated(problem, node, &ordered)
}

/// Evaluate `node` with a complete binding frame keyed by [`VarId`].
/// Entry order is irrelevant. Unknown, duplicate, and missing ids are
/// refused before any graph arithmetic, and each accepted value
/// inherits the units and manifold declared by its variable.
///
/// # Errors
/// The graph, payload, and arithmetic errors documented by [`eval`];
/// keyed structural failures use [`OptError::BindingDuplicate`] and
/// [`OptError::BindingMissing`] instead of positional
/// [`OptError::BindingCount`].
pub fn eval_keyed<'a, I>(problem: &Problem, node: NodeId, bindings: I) -> Result<Value, OptError>
where
    I: IntoIterator<Item = (VarId, &'a [f64])>,
{
    let caps = validate_eval_envelope(problem, node)?;
    let frame = BindingFrame::new_with_caps(problem, bindings, &caps)?;
    eval_validated(problem, node, &frame.ordered)
}

fn validate_eval_envelope(problem: &Problem, node: NodeId) -> Result<AdmissionCaps, OptError> {
    if node.0 as usize >= problem.exprs.len() {
        return Err(OptError::UnknownNode { id: node.0 });
    }
    let caps = AdmissionCaps::default();
    let depth = problem.node_depth(node)?;
    if depth > caps.max_graph_depth {
        return Err(OptError::CapExceeded {
            what: "graph depth",
            count: u64::from(depth),
            cap: u64::from(caps.max_graph_depth),
        });
    }
    validate_binding_frame_envelope(problem, &caps)?;
    Ok(caps)
}

fn validate_binding_frame_envelope(
    problem: &Problem,
    caps: &AdmissionCaps,
) -> Result<(), OptError> {
    let variable_count = problem.vars.len() as u64;
    if variable_count > u64::from(caps.max_vars) {
        return Err(OptError::CapExceeded {
            what: "runtime binding variables",
            count: variable_count,
            cap: u64::from(caps.max_vars),
        });
    }
    let work = problem.total_admission_work();
    if work > caps.max_total_work {
        return Err(OptError::CapExceeded {
            what: "total admission work",
            count: work,
            cap: caps.max_total_work,
        });
    }
    let mut point_storage = 0u64;
    let mut validation_work = 0u64;
    for variable in &problem.vars {
        let point_dim = variable
            .manifold
            .point_dim()
            .ok_or_else(|| OptError::ManifoldInvalid {
                what: format!(
                    "{:?} has no representable runtime point dimension",
                    variable.manifold
                ),
            })?;
        if point_dim > caps.max_point_dim {
            return Err(OptError::CapExceeded {
                what: "runtime binding point dimension",
                count: u64::from(point_dim),
                cap: u64::from(caps.max_point_dim),
            });
        }
        point_storage = point_storage.saturating_add(u64::from(point_dim));
        if point_storage > caps.max_total_point_storage {
            return Err(OptError::CapExceeded {
                what: "runtime binding point storage",
                count: point_storage,
                cap: caps.max_total_point_storage,
            });
        }
        let domain_work = match variable.manifold {
            Manifold::Rn { .. } => 0,
            Manifold::Sphere { .. } | Manifold::So3 => u64::from(point_dim),
            Manifold::Stiefel { n, p } => {
                let gram_entries = u64::from(p).saturating_mul(u64::from(p).saturating_add(1)) / 2;
                u64::from(n).saturating_mul(gram_entries)
            }
        };
        validation_work = validation_work
            .saturating_add(u64::from(point_dim))
            .saturating_add(domain_work);
        if validation_work > caps.max_total_work {
            return Err(OptError::CapExceeded {
                what: "runtime binding validation work",
                count: validation_work,
                cap: caps.max_total_work,
            });
        }
    }
    Ok(())
}

fn validate_ordered_bindings(
    problem: &Problem,
    bindings: &[&[f64]],
    caps: &AdmissionCaps,
) -> Result<(), OptError> {
    let mut total_binding_components = 0u64;
    for (i, (binding, var)) in bindings.iter().zip(&problem.vars).enumerate() {
        let expected = var
            .manifold
            .point_dim()
            .expect("sealed problems carry validated manifolds");
        if binding.len() as u64 != u64::from(expected) {
            return Err(OptError::BindingLen {
                var: i as u32,
                expected,
                got: binding.len() as u64,
            });
        }
        total_binding_components = total_binding_components.saturating_add(binding.len() as u64);
        if total_binding_components > caps.max_total_point_storage {
            return Err(OptError::CapExceeded {
                what: "runtime binding components",
                count: total_binding_components,
                cap: caps.max_total_point_storage,
            });
        }
        for (component_index, component) in binding.iter().enumerate() {
            if !component.is_finite() {
                return Err(OptError::BindingNonFinite {
                    var: i as u32,
                    component: component_index as u32,
                    bits: component.to_bits(),
                });
            }
        }
        var.manifold.validate_binding_domain(binding, i as u32)?;
    }
    Ok(())
}

fn eval_validated(problem: &Problem, node: NodeId, bindings: &[&[f64]]) -> Result<Value, OptError> {
    // Arena order guarantees every dependency has a lower id, so a
    // root-bounded memo prefix is sufficient; unrelated later nodes
    // cannot force memo allocation for this evaluation.
    //
    // EXPLICIT-STACK EVALUATION (bead frankensim-xf8v7): the walk is
    // physically iterative so no admitted graph — at the depth cap or
    // otherwise — can overflow the call stack. Reachability first (a
    // worklist; children always have lower ids), then one bottom-up
    // sweep in arena order, which is a topological order by
    // construction. Only reachable nodes are evaluated, so unreachable
    // Unevaluable nodes still cannot poison an evaluation, exactly as
    // under the recursive walk.
    let root = node.0 as usize;
    let mut memo: Vec<Option<Value>> = vec![None; root + 1];
    let mut reachable = vec![false; root + 1];
    let mut worklist = vec![node];
    while let Some(n) = worklist.pop() {
        let i = n.0 as usize;
        if !reachable[i] {
            reachable[i] = true;
            worklist.extend(crate::ir::children(&problem.exprs[i]));
        }
    }
    for i in 0..=root {
        if reachable[i] {
            let value = eval_node(problem, NodeId(i as u32), bindings, &memo)?;
            memo[i] = Some(value);
        }
    }
    Ok(memo[root]
        .take()
        .expect("the root is reachable from itself"))
}

/// Evaluate ONE node whose children are already in `memo` (guaranteed
/// by the bottom-up arena-order sweep in [`eval`]). Never recurses.
#[allow(clippy::too_many_lines)] // one arm per node kind: the evaluator IS the semantics
fn eval_node(
    problem: &Problem,
    node: NodeId,
    bindings: &[&[f64]],
    memo: &[Option<Value>],
) -> Result<Value, OptError> {
    let ev = |n: NodeId| -> Value {
        memo[n.0 as usize]
            .clone()
            .expect("arena order: children are evaluated before their parents")
    };
    let scalar = |v: Value| -> f64 {
        match v {
            Value::S(x) => x,
            Value::V(_) => unreachable!("builder enforced scalar shape"),
        }
    };
    let out = match &problem.exprs[node.0 as usize] {
        Expr::Var(v) => {
            let x = bindings
                .get(v.0 as usize)
                .ok_or(OptError::UnknownVar { id: v.0 })?;
            Value::V((*x).to_vec())
        }
        Expr::Component { of, index } => {
            let v = ev(*of);
            match v {
                Value::V(xs) => Value::S(checked_component(node, *index, &xs)?),
                Value::S(_) => unreachable!("builder enforced vector shape"),
            }
        }
        Expr::Const { value, .. } => Value::S(*value),
        Expr::Add(a, b) => match (ev(*a), ev(*b)) {
            (Value::S(x), Value::S(y)) => Value::S(x + y),
            (Value::V(x), Value::V(y)) => Value::V(x.iter().zip(&y).map(|(p, q)| p + q).collect()),
            _ => unreachable!("builder enforced matching shapes"),
        },
        Expr::Sub(a, b) => match (ev(*a), ev(*b)) {
            (Value::S(x), Value::S(y)) => Value::S(x - y),
            (Value::V(x), Value::V(y)) => Value::V(x.iter().zip(&y).map(|(p, q)| p - q).collect()),
            _ => unreachable!("builder enforced matching shapes"),
        },
        Expr::Mul(a, b) => match (ev(*a), ev(*b)) {
            (Value::S(x), Value::S(y)) => Value::S(x * y),
            (Value::S(s), Value::V(v)) | (Value::V(v), Value::S(s)) => {
                Value::V(v.iter().map(|p| p * s).collect())
            }
            _ => unreachable!("builder rejected vector*vector"),
        },
        Expr::Div(a, b) => {
            let (x, y) = (scalar(ev(*a)), scalar(ev(*b)));
            Value::S(x / y)
        }
        Expr::Neg(a) => match ev(*a) {
            Value::S(x) => Value::S(-x),
            Value::V(v) => Value::V(v.iter().map(|p| -p).collect()),
        },
        Expr::Powi { base, exp } => Value::S(fs_math::det::powi(scalar(ev(*base)), *exp)),
        Expr::Sqrt(a) => Value::S(fs_math::det::sqrt(scalar(ev(*a)))),
        Expr::Exp(a) => Value::S(fs_math::det::exp(scalar(ev(*a)))),
        Expr::Ln(a) => Value::S(fs_math::det::ln(scalar(ev(*a)))),
        Expr::Tanh(a) => Value::S(fs_math::det::tanh(scalar(ev(*a)))),
        Expr::Dot(a, b) => match (ev(*a), ev(*b)) {
            (Value::V(x), Value::V(y)) => Value::S(x.iter().zip(&y).map(|(p, q)| p * q).sum()),
            _ => unreachable!("builder enforced vectors"),
        },
        Expr::NormSq(a) => match ev(*a) {
            Value::V(x) => Value::S(x.iter().map(|p| p * p).sum()),
            Value::S(_) => unreachable!("builder enforced vector"),
        },
        Expr::Min(a, b) => Value::S(scalar(ev(*a)).min(scalar(ev(*b)))),
        Expr::Max(a, b) => Value::S(scalar(ev(*a)).max(scalar(ev(*b)))),
        Expr::Abs(a) => Value::S(scalar(ev(*a)).abs()),
        Expr::PdeResidual { study, .. } => {
            return Err(OptError::Unevaluable {
                node: node.0,
                kind: format!("pde_residual `{study}` (FLUX executes physics, not the IR)"),
            });
        }
        Expr::Expectation { uq_config, .. }
        | Expr::Cvar { uq_config, .. }
        | Expr::Quantile { uq_config, .. } => {
            return Err(OptError::Unevaluable {
                node: node.0,
                kind: format!("stochastic node over `{uq_config}` (UQ runners execute these)"),
            });
        }
    };
    validate_runtime_shape(problem, node, &out)?;
    match &out {
        Value::S(value) if !value.is_finite() => {
            return Err(OptError::EvalNonFinite {
                node: node.0,
                component: None,
                bits: value.to_bits(),
            });
        }
        Value::V(values) => {
            if let Some((component, value)) = values
                .iter()
                .enumerate()
                .find(|(_, value)| !value.is_finite())
            {
                return Err(OptError::EvalNonFinite {
                    node: node.0,
                    component: Some(component as u32),
                    bits: value.to_bits(),
                });
            }
        }
        Value::S(_) => {}
    }
    Ok(out)
}

impl Manifold {
    /// Descent parameter dimension (what the FD gradient has): ambient
    /// storage for Rn/Sphere/Stiefel (projection happens inside the
    /// retraction), axis-angle 3 for SO(3). CHECKED like
    /// [`Manifold::point_dim`]; `None` only for descriptors a sealed
    /// problem can never contain.
    #[must_use]
    pub fn param_dim(&self) -> Option<u32> {
        match *self {
            Manifold::So3 => Some(3),
            m => m.point_dim(),
        }
    }

    fn label(&self) -> &'static str {
        match self {
            Manifold::Rn { .. } => "Rn",
            Manifold::Sphere { .. } => "Sphere",
            Manifold::So3 => "SO(3)",
            Manifold::Stiefel { .. } => "Stiefel",
        }
    }

    fn domain_error(
        &self,
        what: &'static str,
        location: Option<(u32, u32)>,
        measurement: f64,
    ) -> OptError {
        OptError::RetractionDomain {
            manifold: self.label(),
            what,
            location,
            measurement_bits: measurement.to_bits(),
        }
    }

    fn validate_retraction_finite_with_checkpoint<F>(
        input: &'static str,
        values: &[f64],
        checkpoint: &mut F,
    ) -> Result<(), OptError>
    where
        F: FnMut() -> Result<(), OptError>,
    {
        for (component, value) in values.iter().enumerate() {
            checkpoint_retraction_work(component, checkpoint)?;
            if !value.is_finite() {
                return Err(OptError::RetractionNonFinite {
                    input,
                    component: component as u32,
                    bits: value.to_bits(),
                });
            }
        }
        checkpoint()
    }

    fn validate_point_domain_with_checkpoint<F>(
        &self,
        x: &[f64],
        checkpoint: &mut F,
    ) -> Result<(), OptError>
    where
        F: FnMut() -> Result<(), OptError>,
    {
        match *self {
            Manifold::Rn { .. } => checkpoint(),
            Manifold::Sphere { .. } | Manifold::So3 => {
                let mut norm_sq = 0.0;
                for (component, value) in x.iter().enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    norm_sq += value * value;
                }
                checkpoint()?;
                if !norm_sq.is_finite() {
                    return Err(self.domain_error(
                        "point norm squared must be finite",
                        None,
                        norm_sq,
                    ));
                }
                if (norm_sq - 1.0).abs() > MANIFOLD_MEMBERSHIP_TOL {
                    return Err(self.domain_error(
                        "point must have unit norm before retraction",
                        None,
                        norm_sq,
                    ));
                }
                Ok(())
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                for column in 0..p {
                    for against in 0..=column {
                        checkpoint()?;
                        let mut dot = 0.0;
                        for row in 0..n {
                            checkpoint_retraction_work(row, checkpoint)?;
                            dot += x[column * n + row] * x[against * n + row];
                        }
                        checkpoint()?;
                        let location = Some((column as u32, against as u32));
                        if !dot.is_finite() {
                            return Err(self.domain_error(
                                "point Gram entry must be finite",
                                location,
                                dot,
                            ));
                        }
                        let expected = if column == against { 1.0 } else { 0.0 };
                        if (dot - expected).abs() > MANIFOLD_MEMBERSHIP_TOL {
                            return Err(self.domain_error(
                                "point columns must be orthonormal",
                                location,
                                dot,
                            ));
                        }
                    }
                }
                checkpoint()
            }
        }
    }

    fn validate_point_domain(&self, x: &[f64]) -> Result<(), OptError> {
        self.validate_point_domain_with_checkpoint(x, &mut || Ok(()))
    }

    fn validate_binding_domain(&self, x: &[f64], var: u32) -> Result<(), OptError> {
        match self.validate_point_domain(x) {
            Ok(()) => Ok(()),
            Err(OptError::RetractionDomain {
                manifold,
                what,
                location,
                measurement_bits,
            }) => Err(OptError::BindingDomain {
                var,
                manifold,
                what,
                location,
                measurement_bits,
            }),
            Err(other) => Err(other),
        }
    }

    fn validate_retraction_output_with_checkpoint<F>(
        &self,
        output: Vec<f64>,
        checkpoint: &mut F,
    ) -> Result<Vec<f64>, OptError>
    where
        F: FnMut() -> Result<(), OptError>,
    {
        Self::validate_retraction_finite_with_checkpoint("retraction output", &output, checkpoint)?;
        self.validate_point_domain_with_checkpoint(&output, checkpoint)?;
        Ok(output)
    }

    /// Retract: move from `x` along parameter vector `t`, landing ON
    /// the manifold. Rn: translation. Sphere: normalize(x+t). SO(3):
    /// right-multiply by `exp(ω/2)` (unit quaternion). Stiefel:
    /// Gram-Schmidt of `X+T` (QR retraction). Raw point and parameter
    /// storage must exactly match this manifold and contain finite
    /// components. The base point must already belong to the manifold;
    /// zero-norm and rank-deficient candidates are refused rather than
    /// normalized into fabricated points.
    ///
    /// # Errors
    /// [`OptError::ManifoldInvalid`], [`OptError::RetractionLen`],
    /// [`OptError::RetractionNonFinite`], or [`OptError::RetractionDomain`].
    #[must_use]
    pub fn retract(&self, x: &[f64], t: &[f64]) -> Result<Vec<f64>, OptError> {
        self.retract_with_checkpoint(x, t, &mut || Ok(()))
    }

    #[allow(clippy::too_many_lines)] // one auditable branch per manifold and refusal phase
    fn retract_with_checkpoint<F>(
        &self,
        x: &[f64],
        t: &[f64],
        checkpoint: &mut F,
    ) -> Result<Vec<f64>, OptError>
    where
        F: FnMut() -> Result<(), OptError>,
    {
        checkpoint()?;
        self.validate(&AdmissionCaps::default())?;
        let point_dim = self.point_dim().ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{self:?} has no representable point dimension"),
        })?;
        if x.len() as u64 != u64::from(point_dim) {
            return Err(OptError::RetractionLen {
                input: "retraction point",
                expected: point_dim,
                got: x.len() as u64,
            });
        }
        Self::validate_retraction_finite_with_checkpoint("retraction point", x, checkpoint)?;
        self.validate_point_domain_with_checkpoint(x, checkpoint)?;
        let param_dim = self.param_dim().ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{self:?} has no representable retraction parameter dimension"),
        })?;
        if t.len() as u64 != u64::from(param_dim) {
            return Err(OptError::RetractionLen {
                input: "retraction step",
                expected: param_dim,
                got: t.len() as u64,
            });
        }
        Self::validate_retraction_finite_with_checkpoint("retraction step", t, checkpoint)?;
        match *self {
            Manifold::Rn { .. } => {
                let mut output = Vec::with_capacity(x.len());
                for (component, (a, b)) in x.iter().zip(t).enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    output.push(a + b);
                }
                self.validate_retraction_output_with_checkpoint(output, checkpoint)
            }
            Manifold::Sphere { .. } => {
                let mut y = Vec::with_capacity(x.len());
                for (component, (a, b)) in x.iter().zip(t).enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    y.push(a + b);
                }
                Self::validate_retraction_finite_with_checkpoint(
                    "retraction candidate",
                    &y,
                    checkpoint,
                )?;
                let mut norm_sq = 0.0;
                for (component, value) in y.iter().enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    norm_sq += value * value;
                }
                checkpoint()?;
                if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                    return Err(self.domain_error(
                        "candidate norm squared must be finite and nonsingular",
                        None,
                        norm_sq,
                    ));
                }
                let norm = norm_sq.sqrt();
                let mut output = Vec::with_capacity(y.len());
                for (component, value) in y.iter().enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    output.push(value / norm);
                }
                self.validate_retraction_output_with_checkpoint(output, checkpoint)
            }
            Manifold::So3 => {
                checkpoint()?;
                let half = [t[0] * 0.5, t[1] * 0.5, t[2] * 0.5];
                let angle_sq = half[0] * half[0] + half[1] * half[1] + half[2] * half[2];
                if !angle_sq.is_finite() {
                    return Err(self.domain_error(
                        "step half-angle norm squared must be finite",
                        None,
                        angle_sq,
                    ));
                }
                let ang = angle_sq.sqrt();
                let (s, c) = if ang > 1e-12 {
                    (ang.sin() / ang, ang.cos())
                } else {
                    (1.0, 1.0)
                };
                let dq = [c, half[0] * s, half[1] * s, half[2] * s];
                let q = [x[0], x[1], x[2], x[3]];
                let mut out = [
                    q[0] * dq[0] - q[1] * dq[1] - q[2] * dq[2] - q[3] * dq[3],
                    q[0] * dq[1] + q[1] * dq[0] + q[2] * dq[3] - q[3] * dq[2],
                    q[0] * dq[2] - q[1] * dq[3] + q[2] * dq[0] + q[3] * dq[1],
                    q[0] * dq[3] + q[1] * dq[2] - q[2] * dq[1] + q[3] * dq[0],
                ];
                Self::validate_retraction_finite_with_checkpoint(
                    "retraction candidate",
                    &out,
                    checkpoint,
                )?;
                let mut norm_sq = 0.0;
                for (component, value) in out.iter().enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    norm_sq += value * value;
                }
                checkpoint()?;
                if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                    return Err(self.domain_error(
                        "candidate norm squared must be finite and nonsingular",
                        None,
                        norm_sq,
                    ));
                }
                let norm = norm_sq.sqrt();
                for (component, v) in out.iter_mut().enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    *v /= norm;
                }
                self.validate_retraction_output_with_checkpoint(out.to_vec(), checkpoint)
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                let mut cols = Vec::with_capacity(p);
                for column in 0..p {
                    checkpoint()?;
                    let mut candidate = Vec::with_capacity(n);
                    for row in 0..n {
                        checkpoint_retraction_work(row, checkpoint)?;
                        candidate.push(x[column * n + row] + t[column * n + row]);
                    }
                    Self::validate_retraction_finite_with_checkpoint(
                        "retraction candidate",
                        &candidate,
                        checkpoint,
                    )?;
                    cols.push(candidate);
                }
                // Deterministic Gram-Schmidt (QR retraction).
                for j in 0..p {
                    checkpoint()?;
                    for k in 0..j {
                        checkpoint()?;
                        let mut d = 0.0;
                        for i in 0..n {
                            checkpoint_retraction_work(i, checkpoint)?;
                            d += cols[j][i] * cols[k][i];
                        }
                        checkpoint()?;
                        if !d.is_finite() {
                            return Err(self.domain_error(
                                "candidate Gram projection must be finite",
                                Some((j as u32, k as u32)),
                                d,
                            ));
                        }
                        let (prior_columns, current_and_later) = cols.split_at_mut(j);
                        let prior = &prior_columns[k];
                        let current = &mut current_and_later[0];
                        for i in 0..n {
                            checkpoint_retraction_work(i, checkpoint)?;
                            current[i] -= d * prior[i];
                        }
                        Self::validate_retraction_finite_with_checkpoint(
                            "retraction candidate",
                            current,
                            checkpoint,
                        )?;
                    }
                    let current = &mut cols[j];
                    let mut norm_sq = 0.0;
                    for (component, value) in current.iter().enumerate() {
                        checkpoint_retraction_work(component, checkpoint)?;
                        norm_sq += value * value;
                    }
                    checkpoint()?;
                    if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                        return Err(self.domain_error(
                            "candidate column is rank-deficient",
                            Some((j as u32, j as u32)),
                            norm_sq,
                        ));
                    }
                    let norm = norm_sq.sqrt();
                    for (component, v) in current.iter_mut().enumerate() {
                        checkpoint_retraction_work(component, checkpoint)?;
                        *v /= norm;
                    }
                }
                let mut output = Vec::with_capacity(x.len());
                for column in &cols {
                    for (row, value) in column.iter().enumerate() {
                        checkpoint_retraction_work(row, checkpoint)?;
                        output.push(*value);
                    }
                }
                self.validate_retraction_output_with_checkpoint(output, checkpoint)
            }
        }
    }

    fn retract_with_cx(&self, x: &[f64], t: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, OptError> {
        self.retract_with_checkpoint(x, t, &mut || descent_checkpoint(cx))
    }
}

/// Descent policy.
#[derive(Debug, Clone, Copy)]
pub struct DescentOptions {
    /// Maximum steps.
    pub steps: u32,
    /// Learning rate.
    pub lr: f64,
    /// Finite-difference step.
    pub fd_h: f64,
}

impl Default for DescentOptions {
    fn default() -> Self {
        DescentOptions {
            steps: 200,
            lr: 0.2,
            fd_h: 1e-6,
        }
    }
}

/// What the descent did.
#[derive(Debug, Clone)]
pub struct DescentReport {
    /// Final point.
    pub x: Vec<f64>,
    /// Initial objective value.
    pub f0: f64,
    /// Final objective value.
    pub f_final: f64,
    /// Objective evaluations spent.
    pub evals: u64,
    /// Steps taken.
    pub steps_taken: u32,
    /// True when the P4 budget stopped the run (a receipt, not an
    /// error — the point is still valid).
    pub budget_stopped: bool,
}

/// Toy Riemannian gradient descent of a closure over ONE manifold
/// variable: FD gradient through the retraction, fixed step. Proves
/// retraction metadata is consumable; polls cancellation before objective
/// work and at bounded intervals inside domain/retraction loops, and honors
/// an explicit [`EvalLimit`]. The manifold, start point
/// (length AND finite components), and step policy (`fd_h`/`lr`
/// finite, positive) are leaf-gated BEFORE any descent arithmetic.
/// With an unwind-capable panic strategy, an ordinary raw-objective
/// panic whose hook and payload cleanup both complete normally returns
/// [`OptError::ObjectivePanicked`] with bounded path attribution. The
/// active panic hook still runs first and may emit the original payload
/// and location. `panic=abort`, a panicking/aborting hook, and a panic
/// payload whose destructor panics are process-level no-claim boundaries;
/// caller and hook side effects are not rolled back.
///
/// # Errors
/// [`OptError::Cancelled`] / [`OptError::ManifoldInvalid`] /
/// [`OptError::BindingLen`] / [`OptError::NonFinite`] /
/// [`OptError::ObjectivePanicked`] / [`OptError::BadParam`].
pub fn descend_fn(
    manifold: Manifold,
    f: &dyn Fn(&[f64]) -> f64,
    x0: &[f64],
    opts: DescentOptions,
    eval_limit: EvalLimit,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    let evaluation = std::cell::Cell::new(0u64);
    let checked_f = |x: &[f64], site: ObjectiveEvalSite| {
        let ordinal = evaluation.get().saturating_add(1);
        evaluation.set(ordinal);
        let value =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| f(x))).map_err(|_| {
                OptError::ObjectivePanicked {
                    evaluation: ordinal,
                    site,
                }
            })?;
        finite_descent_value(value, "descent objective result")
    };
    descend_fn_checked(manifold, &checked_f, x0, opts, eval_limit, cx)
}

fn finite_descent_value(value: f64, what: &'static str) -> Result<f64, OptError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OptError::NonFinite {
            what,
            bits: value.to_bits(),
        })
    }
}

fn descent_checkpoint(cx: &Cx<'_>) -> Result<(), OptError> {
    cx.checkpoint().map_err(|_| OptError::Cancelled)
}

fn descend_fn_checked(
    manifold: Manifold,
    f: &dyn Fn(&[f64], ObjectiveEvalSite) -> Result<f64, OptError>,
    x0: &[f64],
    opts: DescentOptions,
    eval_limit: EvalLimit,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    manifold.validate(&AdmissionCaps::default())?;
    let point_dim = manifold
        .point_dim()
        .ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{manifold:?} has no representable point dimension"),
        })?;
    if x0.len() as u64 != u64::from(point_dim) {
        return Err(OptError::BindingLen {
            var: 0,
            expected: point_dim,
            got: x0.len() as u64,
        });
    }
    if !(opts.fd_h.is_finite() && opts.fd_h > 0.0) {
        return Err(OptError::BadParam {
            what: "descent finite-difference step fd_h (finite, > 0)",
            value: opts.fd_h,
        });
    }
    if !(opts.lr.is_finite() && opts.lr > 0.0) {
        return Err(OptError::BadParam {
            what: "descent learning rate lr (finite, > 0; descent, not ascent)",
            value: opts.lr,
        });
    }
    let param_dim = manifold
        .param_dim()
        .ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{manifold:?} has no representable descent parameter dimension"),
        })?;
    let pd = usize::try_from(param_dim).map_err(|_| OptError::ManifoldInvalid {
        what: format!("{manifold:?} descent parameter dimension does not fit usize"),
    })?;
    let atomic_step_evals = u64::from(param_dim).saturating_mul(2).saturating_add(1);
    // Preserve deterministic cheap-refusal precedence for malformed manifold,
    // length, and policy metadata. Long point scans begin only after those O(1)
    // gates, and are cancellation-aware from their first component onward.
    descent_checkpoint(cx)?;
    // Leaf gating (review High #6, bead j3vb5): a non-finite start
    // point or degenerate step policy must refuse BEFORE any descent
    // arithmetic — NaN would otherwise propagate through retractions
    // and finite differences as plausible-looking garbage.
    for (index, component) in x0.iter().enumerate() {
        if index % RETRACTION_CHECKPOINT_STRIDE == 0 {
            descent_checkpoint(cx)?;
        }
        if !component.is_finite() {
            return Err(OptError::NonFinite {
                what: "descent initial point component",
                bits: component.to_bits(),
            });
        }
    }
    descent_checkpoint(cx)?;
    manifold.validate_point_domain_with_checkpoint(x0, &mut || descent_checkpoint(cx))?;
    descent_checkpoint(cx)?;
    let mut x = x0.to_vec();
    let mut evals = 0u64;
    let mut budget_stopped = false;
    let f0 = f(&x, ObjectiveEvalSite::Initial)?;
    evals += 1;
    descent_checkpoint(cx)?;
    let mut steps_taken = 0;
    'outer: for step_index in 0..opts.steps {
        descent_checkpoint(cx)?;
        // One landed step is atomic: reserve its complete central-FD
        // gradient (two probes per parameter) plus the terminal value
        // before spending the first probe. A one-short cap therefore
        // leaves both x and the evaluation ledger at the prior receipt.
        if matches!(
            eval_limit,
            EvalLimit::Limited(maximum)
                if evals.saturating_add(atomic_step_evals) > maximum.get()
        ) {
            budget_stopped = true;
            break 'outer;
        }
        let mut g = vec![0.0; pd];
        for (i, gi) in g.iter_mut().enumerate() {
            descent_checkpoint(cx)?;
            let mut t = vec![0.0; pd];
            t[i] = opts.fd_h;
            descent_checkpoint(cx)?;
            let xp = manifold.retract_with_cx(&x, &t, cx)?;
            descent_checkpoint(cx)?;
            let fp = f(
                &xp,
                ObjectiveEvalSite::Probe {
                    step: step_index,
                    parameter: i as u32,
                    direction: ProbeDirection::Positive,
                },
            )?;
            evals += 1;
            descent_checkpoint(cx)?;
            t[i] = -opts.fd_h;
            let xm = manifold.retract_with_cx(&x, &t, cx)?;
            descent_checkpoint(cx)?;
            let fm = f(
                &xm,
                ObjectiveEvalSite::Probe {
                    step: step_index,
                    parameter: i as u32,
                    direction: ProbeDirection::Negative,
                },
            )?;
            evals += 1;
            descent_checkpoint(cx)?;
            *gi = finite_descent_value(
                (fp - fm) / (2.0 * opts.fd_h),
                "finite-difference gradient component",
            )?;
        }
        descent_checkpoint(cx)?;
        let step: Vec<f64> = g.iter().map(|v| -opts.lr * v).collect();
        descent_checkpoint(cx)?;
        x = manifold.retract_with_cx(&x, &step, cx)?;
        descent_checkpoint(cx)?;
        steps_taken += 1;
    }
    let f_final = if steps_taken == 0 {
        f0
    } else {
        descent_checkpoint(cx)?;
        let value = f(&x, ObjectiveEvalSite::Final { steps_taken })?;
        evals += 1;
        descent_checkpoint(cx)?;
        value
    };
    descent_checkpoint(cx)?;
    Ok(DescentReport {
        x,
        f0,
        f_final,
        evals,
        steps_taken,
        budget_stopped,
    })
}

/// Toy descent of a problem's FIRST objective over its FIRST variable
/// (the IR-driven variant; enforces `problem.budget` per P4).
///
/// # Errors
/// [`OptError::Cancelled`] / evaluation teaching errors.
pub fn descend_ir(
    problem: &Problem,
    x0: &[f64],
    opts: DescentOptions,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    // A problem with no objective or no variable is unsolvable — return a
    // structured error, never an index panic (`ProblemBuilder` does not
    // require either, and `descend_ir` is public).
    let obj = *problem
        .objectives
        .first()
        .ok_or(OptError::IndexOut { index: 0, len: 0 })?;
    let sign = match obj.sense {
        crate::ir::Sense::Minimize => 1.0,
        crate::ir::Sense::Maximize => -1.0,
    };
    let manifold = problem
        .vars
        .first()
        .ok_or(OptError::IndexOut { index: 0, len: 0 })?
        .manifold;
    // The fallible objective runs only inside the shared descent seam,
    // after manifold/start/options leaf validation. Its first successful
    // invocation is f0 and is counted exactly once, including under a
    // one-evaluation budget.
    let f = |x: &[f64], _site: ObjectiveEvalSite| -> Result<f64, OptError> {
        let scalar = eval(problem, obj.node, &[x.to_vec()])?
            .scalar()
            .ok_or(OptError::NotScalar { node: obj.node.0 })?;
        finite_descent_value(sign * obj.weight * scalar, "weighted IR objective")
    };
    descend_fn_checked(manifold, &f, x0, opts, problem.budget.limit, cx)
}

#[cfg(test)]
mod runtime_shape_tests {
    use super::*;
    use crate::ir::ProblemBuilder;
    use fs_qty::Dims;

    #[test]
    fn runtime_values_must_match_sealed_shape_receipts() {
        let mut builder = ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
            .expect("variable");
        let vector = builder.var_ref(variable).expect("vector node");
        let scalar = builder.component(vector, 0).expect("scalar node");
        let problem = builder.finish();

        assert_eq!(
            validate_runtime_shape(&problem, vector, &Value::V(vec![1.0, 2.0])),
            Ok(())
        );
        assert_eq!(
            validate_runtime_shape(&problem, scalar, &Value::S(1.0)),
            Ok(())
        );
        for (actual_len, actual_len_u64) in [(1usize, 1u64), (3, 3)] {
            assert_eq!(
                validate_runtime_shape(&problem, vector, &Value::V(vec![0.0; actual_len])),
                Err(OptError::EvalShape {
                    node: vector.0,
                    expected: Shape::Vector(2),
                    actual_vector_len: Some(actual_len_u64),
                })
            );
        }
        assert_eq!(
            validate_runtime_shape(&problem, vector, &Value::S(0.0)),
            Err(OptError::EvalShape {
                node: vector.0,
                expected: Shape::Vector(2),
                actual_vector_len: None,
            })
        );
        assert_eq!(
            validate_runtime_shape(&problem, scalar, &Value::V(vec![0.0])),
            Err(OptError::EvalShape {
                node: scalar.0,
                expected: Shape::Scalar,
                actual_vector_len: Some(1),
            })
        );
    }

    #[test]
    fn runtime_component_access_is_checked_and_attributed() {
        let node = NodeId(17);
        assert_eq!(checked_component(node, 0, &[4.0]), Ok(4.0));
        assert_eq!(
            checked_component(node, 1, &[4.0]),
            Err(OptError::EvalIndexOut {
                node: 17,
                index: 1,
                len: 1,
            })
        );
    }

    #[test]
    fn retraction_checkpoint_adapter_preserves_bits_errors_and_poll_replay() {
        let cases = [
            (Manifold::Rn { dim: 2 }, vec![1.0, -2.0], vec![0.5, 0.25]),
            (
                Manifold::Sphere { ambient: 3 },
                vec![1.0, 0.0, 0.0],
                vec![0.0, 0.25, -0.125],
            ),
            (
                Manifold::So3,
                vec![1.0, 0.0, 0.0, 0.0],
                vec![0.1, -0.2, 0.3],
            ),
            (
                Manifold::Stiefel { n: 3, p: 2 },
                vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                vec![0.0, 0.125, 0.0, -0.125, 0.0, 0.0],
            ),
        ];

        for (manifold, point, step) in cases {
            let public = manifold.retract(&point, &step).expect("public retraction");
            let mut first_polls = 0usize;
            let first = manifold
                .retract_with_checkpoint(&point, &step, &mut || {
                    first_polls += 1;
                    Ok(())
                })
                .expect("checkpointed retraction");
            let mut replay_polls = 0usize;
            let replay = manifold
                .retract_with_checkpoint(&point, &step, &mut || {
                    replay_polls += 1;
                    Ok(())
                })
                .expect("checkpointed replay");

            let bits = |values: &[f64]| {
                values
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            };
            assert_eq!(
                bits(&first),
                bits(&public),
                "{manifold:?} changed output bits"
            );
            assert_eq!(bits(&replay), bits(&first), "{manifold:?} replay drifted");
            assert_eq!(replay_polls, first_polls, "{manifold:?} poll count drifted");

            let mut invalid_step = step;
            invalid_step[0] = f64::NAN;
            assert_eq!(
                manifold.retract_with_checkpoint(&point, &invalid_step, &mut || Ok(())),
                manifold.retract(&point, &invalid_step),
                "{manifold:?} changed typed error attribution"
            );
        }
    }

    #[test]
    fn stiefel_retraction_work_is_checkpointed_before_and_inside_long_loops() {
        let manifold = Manifold::Stiefel { n: 1024, p: 2 };
        let mut point = vec![0.0; 2048];
        point[0] = 1.0;
        point[1024 + 1] = 1.0;
        let step = vec![0.0; 2048];

        // These deterministic ordinals land at first contact, point-domain
        // validation, QR projection/update, and output revalidation. Each
        // injected cancellation must abort without publishing a point.
        for cancel_at in [1usize, 20, 85, 129] {
            let calls = std::cell::Cell::new(0usize);
            let result = manifold.retract_with_checkpoint(&point, &step, &mut || {
                let next = calls.get() + 1;
                calls.set(next);
                if next == cancel_at {
                    Err(OptError::Cancelled)
                } else {
                    Ok(())
                }
            });
            assert_eq!(result, Err(OptError::Cancelled));
            assert_eq!(calls.get(), cancel_at);
        }
    }
}
