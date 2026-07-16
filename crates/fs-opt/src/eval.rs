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
use core::num::NonZeroU64;
use fs_exec::Cx;

/// Squared-norm/dot tolerance for deciding whether a finite stored
/// point is already on its declared manifold.
const MANIFOLD_MEMBERSHIP_TOL: f64 = 1e-10;
/// Candidates below this squared norm are numerically rank-deficient;
/// normalizing them would fabricate a direction.
const RETRACTION_MIN_NORM_SQ: f64 = 1e-24;
/// Maximum traversed scalar elements between retraction/domain cancellation polls.
const RETRACTION_CHECKPOINT_STRIDE: usize = 256;
/// Maximum evaluator work items between cancellation polls inside one phase.
const EVAL_CHECKPOINT_STRIDE: usize = 256;
const DEFAULT_DESCENT_MAX_WORK_UNITS: u64 = 1 << 24;
const DEFAULT_DESCENT_MAX_WORKSPACE_BYTES: u64 = 1 << 30;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EvalPhase {
    BindingEnvelope,
    BindingValues,
    BindingDomain,
    StorageInitialization,
    Reachability,
    NodeSweep,
    VectorConstruction,
    VectorReduction,
    OutputValidation,
    Finalize,
}

fn no_eval_checkpoint(_phase: EvalPhase) -> Result<(), OptError> {
    Ok(())
}

fn checkpoint_eval_work<F>(
    index: usize,
    phase: EvalPhase,
    checkpoint: &mut F,
) -> Result<(), OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    if index % EVAL_CHECKPOINT_STRIDE == 0 {
        checkpoint(phase)?;
    }
    Ok(())
}

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

fn allocation_len(len: usize) -> u64 {
    u64::try_from(len).unwrap_or(u64::MAX)
}

fn try_vec_capacity<T>(
    path: &'static str,
    node: Option<NodeId>,
    variable: Option<VarId>,
    len: usize,
) -> Result<Vec<T>, OptError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| OptError::RuntimeAllocationRefused {
            path,
            node: node.map(|node| node.0),
            variable: variable.map(|variable| variable.0),
            elements: allocation_len(len),
            element_bytes: u64::try_from(core::mem::size_of::<T>()).unwrap_or(u64::MAX),
        })?;
    Ok(values)
}

fn try_clone_vector<F>(
    path: &'static str,
    node: Option<NodeId>,
    variable: Option<VarId>,
    values: &[f64],
    checkpoint: &mut F,
) -> Result<Vec<f64>, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    checkpoint(EvalPhase::VectorConstruction)?;
    let mut output = try_vec_capacity(path, node, variable, values.len())?;
    for (component, value) in values.iter().enumerate() {
        checkpoint_eval_work(component, EvalPhase::VectorConstruction, checkpoint)?;
        output.push(*value);
    }
    checkpoint(EvalPhase::VectorConstruction)?;
    Ok(output)
}

fn try_owned_diagnostic(
    path: &'static str,
    node: Option<NodeId>,
    variable: Option<VarId>,
    text: &'static str,
) -> Result<String, OptError> {
    let mut output = String::new();
    output
        .try_reserve_exact(text.len())
        .map_err(|_| OptError::RuntimeAllocationRefused {
            path,
            node: node.map(|node| node.0),
            variable: variable.map(|variable| variable.0),
            elements: allocation_len(text.len()),
            element_bytes: 1,
        })?;
    output.push_str(text);
    Ok(output)
}

fn try_map_vector<F>(
    path: &'static str,
    node: Option<NodeId>,
    values: &[f64],
    mut map: impl FnMut(f64) -> f64,
    checkpoint: &mut F,
) -> Result<Vec<f64>, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    checkpoint(EvalPhase::VectorConstruction)?;
    let mut output = try_vec_capacity(path, node, None, values.len())?;
    for (component, value) in values.iter().enumerate() {
        checkpoint_eval_work(component, EvalPhase::VectorConstruction, checkpoint)?;
        output.push(map(*value));
    }
    checkpoint(EvalPhase::VectorConstruction)?;
    Ok(output)
}

fn try_zip_vectors<F>(
    path: &'static str,
    node: Option<NodeId>,
    left: &[f64],
    right: &[f64],
    mut map: impl FnMut(f64, f64) -> f64,
    checkpoint: &mut F,
) -> Result<Vec<f64>, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    checkpoint(EvalPhase::VectorConstruction)?;
    let mut output = try_vec_capacity(path, node, None, left.len())?;
    for (component, (left, right)) in left.iter().zip(right).enumerate() {
        checkpoint_eval_work(component, EvalPhase::VectorConstruction, checkpoint)?;
        output.push(map(*left, *right));
    }
    checkpoint(EvalPhase::VectorConstruction)?;
    Ok(output)
}

fn reduce_dot<F>(left: &[f64], right: &[f64], checkpoint: &mut F) -> Result<f64, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    let mut sum = 0.0;
    for (component, (left, right)) in left.iter().zip(right).enumerate() {
        checkpoint_eval_work(component, EvalPhase::VectorReduction, checkpoint)?;
        sum += *left * *right;
    }
    checkpoint(EvalPhase::VectorReduction)?;
    Ok(sum)
}

fn reduce_norm_sq<F>(values: &[f64], checkpoint: &mut F) -> Result<f64, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    let mut sum = 0.0;
    for (component, value) in values.iter().enumerate() {
        checkpoint_eval_work(component, EvalPhase::VectorReduction, checkpoint)?;
        sum += *value * *value;
    }
    checkpoint(EvalPhase::VectorReduction)?;
    Ok(sum)
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
    /// [`OptError::BindingNonFinite`], [`OptError::BindingDomain`],
    /// [`OptError::CapExceeded`], or [`OptError::RuntimeAllocationRefused`].
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
        let mut slots = try_vec_capacity("binding-frame/slots", None, None, problem.vars.len())?;
        slots.resize(problem.vars.len(), None);
        for (var, binding) in bindings {
            let slot = slots
                .get_mut(var.0 as usize)
                .ok_or(OptError::UnknownVar { id: var.0 })?;
            if slot.replace(binding).is_some() {
                return Err(OptError::BindingDuplicate { var: var.0 });
            }
        }
        if let Some((var, _)) = slots
            .iter()
            .enumerate()
            .find(|(_, binding)| binding.is_none())
        {
            return Err(OptError::BindingMissing { var: var as u32 });
        }
        let mut ordered = try_vec_capacity("binding-frame/order", None, None, slots.len())?;
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
    /// [`OptError::EvalNonFinite`], [`OptError::CapExceeded`], or
    /// [`OptError::RuntimeAllocationRefused`].
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
/// [`OptError::CapExceeded`] / [`OptError::RuntimeAllocationRefused`].
pub fn eval(problem: &Problem, node: NodeId, bindings: &[Vec<f64>]) -> Result<Value, OptError> {
    let caps = validate_eval_envelope(problem, node)?;
    if bindings.len() != problem.vars.len() {
        return Err(OptError::BindingCount {
            vars: problem.vars.len() as u32,
            got: bindings.len() as u64,
        });
    }
    let mut ordered = try_vec_capacity("eval/positional-frame", Some(node), None, bindings.len())?;
    for binding in bindings {
        ordered.push(binding.as_slice());
    }
    eval_ordered_with_caps(problem, node, &ordered, &caps)
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

fn eval_borrowed_with_checkpoint<F>(
    problem: &Problem,
    node: NodeId,
    bindings: &[&[f64]],
    checkpoint: &mut F,
) -> Result<Value, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    let caps = validate_eval_header(problem, node)?;
    validate_binding_frame_header(problem, &caps)?;
    if bindings.len() != problem.vars.len() {
        return Err(OptError::BindingCount {
            vars: problem.vars.len() as u32,
            got: allocation_len(bindings.len()),
        });
    }
    validate_ordered_binding_headers(problem, bindings, &caps)?;
    validate_binding_frame_envelope_with_checkpoint(problem, &caps, checkpoint)?;
    validate_ordered_binding_values_with_checkpoint(problem, bindings, checkpoint)?;
    eval_validated_with_checkpoint(problem, node, bindings, checkpoint)
}

fn eval_ordered_with_caps(
    problem: &Problem,
    node: NodeId,
    bindings: &[&[f64]],
    caps: &AdmissionCaps,
) -> Result<Value, OptError> {
    eval_ordered_with_caps_and_checkpoint(problem, node, bindings, caps, &mut no_eval_checkpoint)
}

fn eval_ordered_with_caps_and_checkpoint<F>(
    problem: &Problem,
    node: NodeId,
    bindings: &[&[f64]],
    caps: &AdmissionCaps,
    checkpoint: &mut F,
) -> Result<Value, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    validate_ordered_bindings_with_checkpoint(problem, bindings, caps, checkpoint)?;
    eval_validated_with_checkpoint(problem, node, bindings, checkpoint)
}

fn validate_eval_envelope(problem: &Problem, node: NodeId) -> Result<AdmissionCaps, OptError> {
    validate_eval_envelope_with_checkpoint(problem, node, &mut no_eval_checkpoint)
}

fn validate_eval_envelope_with_checkpoint<F>(
    problem: &Problem,
    node: NodeId,
    checkpoint: &mut F,
) -> Result<AdmissionCaps, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    let caps = validate_eval_header(problem, node)?;
    validate_binding_frame_envelope_with_checkpoint(problem, &caps, checkpoint)?;
    Ok(caps)
}

fn validate_eval_header(problem: &Problem, node: NodeId) -> Result<AdmissionCaps, OptError> {
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
    Ok(caps)
}

fn validate_binding_frame_envelope(
    problem: &Problem,
    caps: &AdmissionCaps,
) -> Result<(), OptError> {
    validate_binding_frame_envelope_with_checkpoint(problem, caps, &mut no_eval_checkpoint)
}

fn validate_binding_frame_envelope_with_checkpoint<F>(
    problem: &Problem,
    caps: &AdmissionCaps,
    checkpoint: &mut F,
) -> Result<(), OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    validate_binding_frame_header(problem, caps)?;
    let mut point_storage = 0u64;
    let mut validation_work = 0u64;
    for (variable_index, variable) in problem.vars.iter().enumerate() {
        let Some(point_dim) = variable.manifold.point_dim() else {
            return Err(OptError::ManifoldInvalid {
                what: try_owned_diagnostic(
                    "binding-frame/manifold-diagnostic",
                    None,
                    Some(VarId(variable_index as u32)),
                    "sealed manifold has no representable runtime point dimension",
                )?,
            });
        };
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
        checkpoint_eval_work(variable_index, EvalPhase::BindingEnvelope, checkpoint)?;
    }
    checkpoint(EvalPhase::BindingEnvelope)
}

fn validate_binding_frame_header(problem: &Problem, caps: &AdmissionCaps) -> Result<(), OptError> {
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
    Ok(())
}

fn validate_ordered_bindings(
    problem: &Problem,
    bindings: &[&[f64]],
    caps: &AdmissionCaps,
) -> Result<(), OptError> {
    validate_ordered_bindings_with_checkpoint(problem, bindings, caps, &mut no_eval_checkpoint)
}

fn validate_ordered_bindings_with_checkpoint<F>(
    problem: &Problem,
    bindings: &[&[f64]],
    caps: &AdmissionCaps,
    checkpoint: &mut F,
) -> Result<(), OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
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
        validate_binding_value_with_checkpoint(binding, var.manifold, i as u32, checkpoint)?;
    }
    checkpoint(EvalPhase::BindingValues)
}

fn validate_ordered_binding_values_with_checkpoint<F>(
    problem: &Problem,
    bindings: &[&[f64]],
    checkpoint: &mut F,
) -> Result<(), OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    for (i, (binding, var)) in bindings.iter().zip(&problem.vars).enumerate() {
        validate_binding_value_with_checkpoint(binding, var.manifold, i as u32, checkpoint)?;
    }
    checkpoint(EvalPhase::BindingValues)
}

fn validate_binding_value_with_checkpoint<F>(
    binding: &[f64],
    manifold: Manifold,
    variable: u32,
    checkpoint: &mut F,
) -> Result<(), OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    for (component_index, component) in binding.iter().enumerate() {
        checkpoint_eval_work(component_index, EvalPhase::BindingValues, checkpoint)?;
        if !component.is_finite() {
            return Err(OptError::BindingNonFinite {
                var: variable,
                component: component_index as u32,
                bits: component.to_bits(),
            });
        }
    }
    checkpoint(EvalPhase::BindingValues)?;
    manifold.validate_binding_domain_with_checkpoint(binding, variable, &mut || {
        checkpoint(EvalPhase::BindingDomain)
    })
}

fn validate_ordered_binding_headers(
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
    }
    Ok(())
}

fn eval_validated(problem: &Problem, node: NodeId, bindings: &[&[f64]]) -> Result<Value, OptError> {
    eval_validated_with_checkpoint(problem, node, bindings, &mut no_eval_checkpoint)
}

fn eval_validated_with_checkpoint<F>(
    problem: &Problem,
    node: NodeId,
    bindings: &[&[f64]],
    checkpoint: &mut F,
) -> Result<Value, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
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
    let prefix_len = root + 1;
    checkpoint(EvalPhase::StorageInitialization)?;
    let mut memo = try_vec_capacity("eval/memo", Some(node), None, prefix_len)?;
    for index in 0..prefix_len {
        checkpoint_eval_work(index, EvalPhase::StorageInitialization, checkpoint)?;
        memo.push(None);
    }
    checkpoint(EvalPhase::StorageInitialization)?;
    let mut reachable = try_vec_capacity("eval/reachability", Some(node), None, prefix_len)?;
    for index in 0..prefix_len {
        checkpoint_eval_work(index, EvalPhase::StorageInitialization, checkpoint)?;
        reachable.push(false);
    }
    checkpoint(EvalPhase::StorageInitialization)?;
    let mut worklist = try_vec_capacity("eval/worklist", Some(node), None, prefix_len)?;
    reachable[root] = true;
    worklist.push(node);
    let mut reachability_work = 0usize;
    while let Some(n) = worklist.pop() {
        checkpoint_eval_work(reachability_work, EvalPhase::Reachability, checkpoint)?;
        reachability_work = reachability_work.saturating_add(1);
        let i = n.0 as usize;
        for child in crate::ir::children(&problem.exprs[i]) {
            checkpoint_eval_work(reachability_work, EvalPhase::Reachability, checkpoint)?;
            reachability_work = reachability_work.saturating_add(1);
            let child_index = child.0 as usize;
            if !reachable[child_index] {
                reachable[child_index] = true;
                worklist.push(child);
            }
        }
    }
    checkpoint(EvalPhase::Reachability)?;
    for i in 0..=root {
        checkpoint_eval_work(i, EvalPhase::NodeSweep, checkpoint)?;
        if reachable[i] {
            let value =
                eval_node_with_checkpoint(problem, NodeId(i as u32), bindings, &memo, checkpoint)?;
            memo[i] = Some(value);
        }
    }
    checkpoint(EvalPhase::NodeSweep)?;
    let value = memo[root]
        .take()
        .expect("the root is reachable from itself");
    checkpoint(EvalPhase::Finalize)?;
    Ok(value)
}

/// Evaluate ONE node whose children are already in `memo` (guaranteed
/// by the bottom-up arena-order sweep in [`eval`]). Never recurses.
#[allow(clippy::too_many_lines)] // one arm per node kind: the evaluator IS the semantics
fn eval_node_with_checkpoint<F>(
    problem: &Problem,
    node: NodeId,
    bindings: &[&[f64]],
    memo: &[Option<Value>],
    checkpoint: &mut F,
) -> Result<Value, OptError>
where
    F: FnMut(EvalPhase) -> Result<(), OptError>,
{
    let ev = |n: NodeId| -> &Value {
        memo[n.0 as usize]
            .as_ref()
            .expect("arena order: children are evaluated before their parents")
    };
    let scalar = |v: &Value| -> f64 {
        match v {
            Value::S(x) => *x,
            Value::V(_) => unreachable!("builder enforced scalar shape"),
        }
    };
    let out =
        match &problem.exprs[node.0 as usize] {
            Expr::Var(v) => {
                let x = bindings
                    .get(v.0 as usize)
                    .ok_or(OptError::UnknownVar { id: v.0 })?;
                Value::V(try_clone_vector(
                    "eval/variable-value",
                    Some(node),
                    Some(*v),
                    x,
                    checkpoint,
                )?)
            }
            Expr::Component { of, index } => {
                let v = ev(*of);
                match v {
                    Value::V(xs) => Value::S(checked_component(node, *index, xs)?),
                    Value::S(_) => unreachable!("builder enforced vector shape"),
                }
            }
            Expr::Const { value, .. } => Value::S(*value),
            Expr::Add(a, b) => match (ev(*a), ev(*b)) {
                (Value::S(x), Value::S(y)) => Value::S(*x + *y),
                (Value::V(x), Value::V(y)) => Value::V(try_zip_vectors(
                    "eval/vector-add",
                    Some(node),
                    x,
                    y,
                    |p, q| p + q,
                    checkpoint,
                )?),
                _ => unreachable!("builder enforced matching shapes"),
            },
            Expr::Sub(a, b) => match (ev(*a), ev(*b)) {
                (Value::S(x), Value::S(y)) => Value::S(*x - *y),
                (Value::V(x), Value::V(y)) => Value::V(try_zip_vectors(
                    "eval/vector-sub",
                    Some(node),
                    x,
                    y,
                    |p, q| p - q,
                    checkpoint,
                )?),
                _ => unreachable!("builder enforced matching shapes"),
            },
            Expr::Mul(a, b) => match (ev(*a), ev(*b)) {
                (Value::S(x), Value::S(y)) => Value::S(*x * *y),
                (Value::S(s), Value::V(v)) | (Value::V(v), Value::S(s)) => Value::V(
                    try_map_vector("eval/vector-scale", Some(node), v, |p| p * *s, checkpoint)?,
                ),
                _ => unreachable!("builder rejected vector*vector"),
            },
            Expr::Div(a, b) => {
                let (x, y) = (scalar(ev(*a)), scalar(ev(*b)));
                Value::S(x / y)
            }
            Expr::Neg(a) => match ev(*a) {
                Value::S(x) => Value::S(-*x),
                Value::V(v) => Value::V(try_map_vector(
                    "eval/vector-negate",
                    Some(node),
                    v,
                    |p| -p,
                    checkpoint,
                )?),
            },
            Expr::Powi { base, exp } => Value::S(fs_math::det::powi(scalar(ev(*base)), *exp)),
            Expr::Sqrt(a) => Value::S(fs_math::det::sqrt(scalar(ev(*a)))),
            Expr::Exp(a) => Value::S(fs_math::det::exp(scalar(ev(*a)))),
            Expr::Ln(a) => Value::S(fs_math::det::ln(scalar(ev(*a)))),
            Expr::Tanh(a) => Value::S(fs_math::det::tanh(scalar(ev(*a)))),
            Expr::Dot(a, b) => match (ev(*a), ev(*b)) {
                (Value::V(x), Value::V(y)) => Value::S(reduce_dot(x, y, checkpoint)?),
                _ => unreachable!("builder enforced vectors"),
            },
            Expr::NormSq(a) => match ev(*a) {
                Value::V(x) => Value::S(reduce_norm_sq(x, checkpoint)?),
                Value::S(_) => unreachable!("builder enforced vector"),
            },
            Expr::Min(a, b) => Value::S(scalar(ev(*a)).min(scalar(ev(*b)))),
            Expr::Max(a, b) => Value::S(scalar(ev(*a)).max(scalar(ev(*b)))),
            Expr::Abs(a) => Value::S(scalar(ev(*a)).abs()),
            Expr::PdeResidual { .. } => {
                return Err(OptError::Unevaluable {
                    node: node.0,
                    kind: "pde_residual (FLUX executes physics, not the IR)",
                });
            }
            Expr::Expectation { .. } | Expr::Cvar { .. } | Expr::Quantile { .. } => {
                return Err(OptError::Unevaluable {
                    node: node.0,
                    kind: "stochastic node (UQ runners execute these, not the IR)",
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
            for (component, value) in values.iter().enumerate() {
                checkpoint_eval_work(component, EvalPhase::OutputValidation, checkpoint)?;
                if !value.is_finite() {
                    return Err(OptError::EvalNonFinite {
                        node: node.0,
                        component: Some(component as u32),
                        bits: value.to_bits(),
                    });
                }
            }
            checkpoint(EvalPhase::OutputValidation)?;
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

    fn validate_binding_domain_with_checkpoint<F>(
        &self,
        x: &[f64],
        var: u32,
        checkpoint: &mut F,
    ) -> Result<(), OptError>
    where
        F: FnMut() -> Result<(), OptError>,
    {
        match self.validate_point_domain_with_checkpoint(x, checkpoint) {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DescentEnvelope {
    owned_work: u64,
    scratch_bytes: u64,
}

fn maximum_landed_steps(steps: u32, param_dim: u32, eval_limit: EvalLimit) -> u64 {
    let requested = u64::from(steps);
    match eval_limit {
        EvalLimit::Unlimited => requested,
        EvalLimit::Limited(maximum) => {
            // N landed steps consume 1 initial + 2*param_dim*N probes +
            // 1 terminal evaluation. The one terminal slot is reserved before
            // any probe, so limits below that complete envelope admit no step.
            let probe_evals = u64::from(param_dim) * 2;
            maximum
                .get()
                .saturating_sub(2)
                .checked_div(probe_evals)
                .unwrap_or(0)
                .min(requested)
        }
    }
}

fn plan_add(resource: &'static str, left: u64, right: u64) -> Result<u64, OptError> {
    left.checked_add(right)
        .ok_or(OptError::DescentPlanOverflow { resource })
}

fn plan_mul(resource: &'static str, left: u64, right: u64) -> Result<u64, OptError> {
    left.checked_mul(right)
        .ok_or(OptError::DescentPlanOverflow { resource })
}

fn descent_envelope(
    manifold: Manifold,
    point_dim: u32,
    param_dim: u32,
    steps: u32,
    eval_limit: EvalLimit,
) -> Result<DescentEnvelope, OptError> {
    let point = u64::from(point_dim);
    let param = u64::from(param_dim);
    let active_steps = maximum_landed_steps(steps, param_dim, eval_limit);

    // One unit conservatively represents one scalar slot initialized or one
    // scalar value visited by fs-opt's own descent/retraction plumbing. The
    // Stiefel multiplier covers input/output Gram scans and deterministic QR.
    // Caller-owned objective work is intentionally outside this envelope.
    let retraction_factor = match manifold {
        Manifold::Stiefel { p, .. } => {
            plan_add("work units", plan_mul("work units", u64::from(p), 8)?, 32)?
        }
        Manifold::Rn { .. } | Manifold::Sphere { .. } | Manifold::So3 => 24,
    };
    let retraction_work = plan_add(
        "work units",
        plan_mul("work units", point, retraction_factor)?,
        param,
    )?;
    let retractions_per_step = plan_add("work units", plan_mul("work units", param, 2)?, 1)?;
    let probe_initialization = plan_mul("work units", param, param)?;
    let retraction_visits = plan_mul("work units", retractions_per_step, retraction_work)?;
    let gradient_plumbing = plan_add(
        "work units",
        plan_mul("work units", param, 4)?,
        plan_mul("work units", point, 2)?, // point scale + landed displacement
    )?;
    let per_step = plan_add(
        "work units",
        plan_add("work units", probe_initialization, retraction_visits)?,
        gradient_plumbing,
    )?;
    let probe_preflight = if active_steps == 0 {
        0
    } else {
        plan_add(
            "work units",
            param,
            plan_mul(
                "work units",
                plan_mul("work units", param, 2)?,
                retraction_work,
            )?,
        )?
    };
    let owned_work = plan_add(
        "work units",
        plan_add(
            "work units",
            plan_add("work units", retraction_work, point)?,
            probe_preflight,
        )?,
        plan_mul("work units", active_steps, per_step)?,
    )?;

    // Peak retained scalar storage is deliberately conservative: current
    // iterate, positive/negative probes, retraction scratch/output, gradient,
    // coordinate probe, and landed step. Include Vec headers, including one
    // per Stiefel column. Caller inputs and objective-owned storage are not
    // charged because fs-opt neither allocates nor controls them.
    let active_scalar_slots = if active_steps == 0 {
        point
    } else {
        plan_add(
            "workspace bytes",
            plan_mul("workspace bytes", point, 6)?,
            plan_mul("workspace bytes", param, 3)?,
        )?
    };
    let column_headers = match manifold {
        Manifold::Stiefel { p, .. } if active_steps > 0 => u64::from(p),
        _ => 0,
    };
    let vec_headers = plan_add("workspace bytes", 16, column_headers)?;
    let scratch_bytes = plan_add(
        "workspace bytes",
        plan_mul(
            "workspace bytes",
            active_scalar_slots,
            core::mem::size_of::<f64>() as u64,
        )?,
        plan_mul(
            "workspace bytes",
            vec_headers,
            core::mem::size_of::<Vec<f64>>() as u64,
        )?,
    )?;

    Ok(DescentEnvelope {
        owned_work,
        scratch_bytes,
    })
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
    /// Unitless relative retracted-point displacement at or below which descent stops.
    pub closure_threshold: f64,
    /// Hard upper bound for descent-engine/retraction scalar work units admitted
    /// up front (default `2^24`).
    /// Objective implementation work, including IR evaluation, has its own
    /// admission boundary and is not included.
    pub max_work_units: NonZeroU64,
    /// Hard upper bound for peak descent-engine and retraction workspace, in
    /// bytes (default 1 GiB). Objective implementation workspace, including
    /// IR evaluation, is not included.
    pub max_workspace_bytes: NonZeroU64,
}

impl Default for DescentOptions {
    fn default() -> Self {
        DescentOptions {
            steps: 200,
            lr: 0.2,
            fd_h: 1e-6,
            closure_threshold: 1e-12,
            max_work_units: NonZeroU64::new(DEFAULT_DESCENT_MAX_WORK_UNITS)
                .expect("positive default descent work cap"),
            max_workspace_bytes: NonZeroU64::new(DEFAULT_DESCENT_MAX_WORKSPACE_BYTES)
                .expect("positive default descent workspace cap"),
        }
    }
}

/// Why a successful descent report stopped.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DescentStop {
    /// The configured iteration limit was reached.
    StepLimit,
    /// The relative retracted-candidate displacement met the configured threshold.
    ClosureThreshold,
    /// The explicit evaluation limit could not fund the next atomic step.
    EvaluationLimit,
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
    /// Typed successful terminal state.
    pub stop: DescentStop,
    /// True when the P4 budget stopped the run (a receipt, not an
    /// error — the point is still valid). Exactly equivalent to
    /// `stop == DescentStop::EvaluationLimit`.
    pub budget_stopped: bool,
    /// Conservative descent-engine/retraction work bound admitted before the
    /// first objective. Objective implementation work is excluded.
    pub work_upper_bound: u64,
    /// Conservative peak descent-engine/retraction workspace bound admitted
    /// before the first objective. Objective implementation workspace is excluded.
    pub workspace_upper_bound_bytes: u64,
}

/// Toy Riemannian gradient descent of a closure over ONE manifold
/// variable: FD gradient through the retraction, fixed step. Proves
/// retraction metadata is consumable; polls cancellation before objective
/// work and at bounded intervals inside domain/retraction loops, and honors
/// an explicit [`EvalLimit`]. The manifold, start point
/// (length AND finite components), step/closure policy, initial coordinate
/// retractions, and descent-engine/retraction work/workspace envelope are gated
/// BEFORE f0.
/// Closure uses the unitless ratio
/// `max_abs(retract(x, step) - x) / max(max_abs(x), fd_h)`.
/// The resource envelope excludes objective work/allocation. For
/// [`descend_ir`], this currently includes the separately admitted, fallibly
/// allocated IR evaluator.
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
/// [`OptError::ObjectivePanicked`] / [`OptError::BadParam`] /
/// [`OptError::DescentCapExceeded`] / [`OptError::DescentPlanOverflow`].
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

fn validate_initial_probe_retractions(
    manifold: Manifold,
    x0: &[f64],
    param_dim: usize,
    fd_h: f64,
    cx: &Cx<'_>,
) -> Result<(), OptError> {
    let mut probe = vec![0.0; param_dim];
    for parameter in 0..param_dim {
        descent_checkpoint(cx)?;
        probe[parameter] = fd_h;
        drop(manifold.retract_with_cx(x0, &probe, cx)?);
        descent_checkpoint(cx)?;
        probe[parameter] = -fd_h;
        drop(manifold.retract_with_cx(x0, &probe, cx)?);
        probe[parameter] = 0.0;
    }
    descent_checkpoint(cx)
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
    if !(opts.fd_h.is_finite() && opts.fd_h > 0.0 && (2.0 * opts.fd_h).is_finite()) {
        return Err(OptError::BadParam {
            what: "descent finite-difference step fd_h (finite, > 0, finite doubled denominator)",
            value: opts.fd_h,
        });
    }
    if !(opts.lr.is_finite() && opts.lr > 0.0 && opts.lr <= 1.0) {
        return Err(OptError::BadParam {
            what: "descent learning rate lr (finite, 0 < lr <= 1; descent, not ascent)",
            value: opts.lr,
        });
    }
    if !(opts.closure_threshold.is_finite()
        && opts.closure_threshold > 0.0
        && opts.closure_threshold <= 1.0)
    {
        return Err(OptError::BadParam {
            what: "descent closure threshold (finite, 0 < threshold <= 1)",
            value: opts.closure_threshold,
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
    let reachable_steps = maximum_landed_steps(opts.steps, param_dim, eval_limit);
    let envelope = descent_envelope(manifold, point_dim, param_dim, opts.steps, eval_limit)?;
    if envelope.owned_work > opts.max_work_units.get() {
        return Err(OptError::DescentCapExceeded {
            resource: "work units",
            required: envelope.owned_work,
            cap: opts.max_work_units.get(),
        });
    }
    if envelope.scratch_bytes > opts.max_workspace_bytes.get() {
        return Err(OptError::DescentCapExceeded {
            resource: "workspace bytes",
            required: envelope.scratch_bytes,
            cap: opts.max_workspace_bytes.get(),
        });
    }
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
    if reachable_steps > 0 {
        validate_initial_probe_retractions(manifold, x0, pd, opts.fd_h, cx)?;
    }
    descent_checkpoint(cx)?;
    let mut x = x0.to_vec();
    let mut evals = 0u64;
    let mut stop = DescentStop::StepLimit;
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
            stop = DescentStop::EvaluationLimit;
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
        let mut step = Vec::with_capacity(g.len());
        for (component, gradient) in g.iter().enumerate() {
            checkpoint_retraction_work(component, &mut || descent_checkpoint(cx))?;
            let value =
                finite_descent_value(-opts.lr * gradient, "descent parameter-step component")?;
            step.push(value);
        }
        descent_checkpoint(cx)?;
        let candidate = manifold.retract_with_cx(&x, &step, cx)?;
        descent_checkpoint(cx)?;
        if candidate.len() != x.len() {
            return Err(OptError::RetractionLen {
                input: "retraction output",
                expected: point_dim,
                got: candidate.len() as u64,
            });
        }
        let mut point_max_abs = 0.0f64;
        let mut landed_step_max_abs = 0.0f64;
        for (component, (candidate_value, value)) in candidate.iter().zip(&x).enumerate() {
            checkpoint_retraction_work(component, &mut || descent_checkpoint(cx))?;
            point_max_abs = point_max_abs.max(value.abs());
            let displacement = finite_descent_value(
                *candidate_value - *value,
                "descent landed-point displacement component",
            )?;
            landed_step_max_abs = landed_step_max_abs.max(displacement.abs());
        }
        let closure_scale = point_max_abs.max(opts.fd_h);
        let relative_step = finite_descent_value(
            landed_step_max_abs / closure_scale,
            "descent relative-step closure ratio",
        )?;
        if relative_step <= opts.closure_threshold {
            stop = DescentStop::ClosureThreshold;
            break 'outer;
        }
        x = candidate;
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
        stop,
        budget_stopped: stop == DescentStop::EvaluationLimit,
        work_upper_bound: envelope.owned_work,
        workspace_upper_bound_bytes: envelope.scratch_bytes,
    })
}

/// Toy descent of a problem's FIRST objective over its FIRST variable
/// (the IR-driven variant; enforces `problem.budget` per P4). The caller's
/// [`Cx`] is threaded through binding validation and every proportional graph
/// evaluation phase after the capped cheap node/count/length preflight, with at
/// most 256 evaluator work items between polls.
///
/// # Errors
/// [`OptError::Cancelled`] / evaluation teaching errors.
pub fn descend_ir(
    problem: &Problem,
    x0: &[f64],
    opts: DescentOptions,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    descend_ir_with_eval_checkpoint(problem, x0, opts, cx, &|_phase| descent_checkpoint(cx))
}

fn descend_ir_with_eval_checkpoint<F>(
    problem: &Problem,
    x0: &[f64],
    opts: DescentOptions,
    cx: &Cx<'_>,
    eval_checkpoint: &F,
) -> Result<DescentReport, OptError>
where
    F: Fn(EvalPhase) -> Result<(), OptError>,
{
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
        let mut checkpoint = eval_checkpoint;
        let scalar = eval_borrowed_with_checkpoint(problem, obj.node, &[x], &mut checkpoint)?
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

    fn checkpoint_problem() -> (Problem, NodeId, Vec<f64>) {
        let dimension = u32::try_from(EVAL_CHECKPOINT_STRIDE * 2 + 1)
            .expect("checkpoint fixture dimension fits u32");
        let mut builder = ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: dimension }, Dims::NONE)
            .expect("variable");
        let point = builder.var_ref(variable).expect("point");
        let negated = builder.neg(point).expect("negation");
        let cancelled = builder.add(point, negated).expect("vector addition");
        let scale = builder.konst(0.5, Dims::NONE).expect("scale");
        let scaled = builder.mul(scale, cancelled).expect("scaled vector");
        let norm = builder.norm_sq(scaled).expect("norm squared");
        let mut leaves = vec![norm];
        for index in 1..=512u32 {
            let value = builder
                .konst(f64::from(index), Dims::NONE)
                .expect("unique padding constant");
            let negative = builder.neg(value).expect("padding negation");
            leaves.push(builder.add(value, negative).expect("exact zero leaf"));
        }
        while leaves.len() > 1 {
            let mut next = Vec::with_capacity(leaves.len().div_ceil(2));
            let mut pairs = leaves.chunks_exact(2);
            for pair in &mut pairs {
                next.push(builder.add(pair[0], pair[1]).expect("balanced scalar sum"));
            }
            if let Some(remainder) = pairs.remainder().first() {
                next.push(*remainder);
            }
            leaves = next;
        }
        let objective = leaves[0];
        builder
            .objective(objective, crate::ir::Sense::Minimize, 1.0)
            .expect("objective");
        let problem = builder.finish();
        let mut binding: Vec<f64> = (0..dimension)
            .map(|index| f64::from(index % 17) - 8.0)
            .collect();
        binding[0] = -0.0;
        (problem, objective, binding)
    }

    fn with_test_cx<R>(gate: &fs_exec::CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                fs_exec::StreamKey {
                    seed: 0xE7A1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                asupersync::types::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn assert_value_bits_eq(left: &Value, right: &Value) {
        match (left, right) {
            (Value::S(left), Value::S(right)) => assert_eq!(left.to_bits(), right.to_bits()),
            (Value::V(left), Value::V(right)) => assert_eq!(
                left.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
                right
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            ),
            _ => panic!("value shape changed under checkpointing"),
        }
    }

    #[test]
    fn checkpointed_eval_replays_and_every_poll_aborts_atomically() {
        let (problem, objective, binding) = checkpoint_problem();
        let public =
            eval(&problem, objective, std::slice::from_ref(&binding)).expect("public evaluation");

        let mut phases = Vec::new();
        let checkpointed = eval_borrowed_with_checkpoint(
            &problem,
            objective,
            &[binding.as_slice()],
            &mut |phase| {
                phases.push(phase);
                Ok(())
            },
        )
        .expect("checkpointed evaluation");
        assert_value_bits_eq(&checkpointed, &public);
        assert_eq!(
            checkpointed.scalar().expect("scalar objective").to_bits(),
            0.0f64.to_bits(),
            "balanced exact-zero padding must retain the expected objective bits"
        );

        let mut replay_phases = Vec::new();
        let replay = eval_borrowed_with_checkpoint(
            &problem,
            objective,
            &[binding.as_slice()],
            &mut |phase| {
                replay_phases.push(phase);
                Ok(())
            },
        )
        .expect("checkpoint replay");
        assert_value_bits_eq(&replay, &checkpointed);
        assert_eq!(replay_phases, phases, "phase trace must replay exactly");

        let expected_histogram = [
            (EvalPhase::BindingEnvelope, 2usize),
            (EvalPhase::BindingValues, 5),
            (EvalPhase::BindingDomain, 1),
            (EvalPhase::StorageInitialization, 21),
            (EvalPhase::Reachability, 20),
            (EvalPhase::NodeSweep, 10),
            (EvalPhase::VectorConstruction, 20),
            (EvalPhase::VectorReduction, 4),
            (EvalPhase::OutputValidation, 16),
            (EvalPhase::Finalize, 1),
        ];
        for (phase, expected) in expected_histogram {
            let actual = phases.iter().filter(|seen| **seen == phase).count();
            assert_eq!(actual, expected, "checkpoint cadence drifted for {phase:?}");
        }
        assert_eq!(phases.len(), 100, "complete phase trace is a G4 receipt");

        for target in 1..=phases.len() {
            let mut calls = 0usize;
            let mut observed = Vec::new();
            let result = eval_borrowed_with_checkpoint(
                &problem,
                objective,
                &[binding.as_slice()],
                &mut |phase| {
                    calls += 1;
                    observed.push(phase);
                    if calls == target {
                        Err(OptError::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            );
            assert_eq!(result, Err(OptError::Cancelled), "target poll {target}");
            assert_eq!(calls, target, "work continued after target poll {target}");
            assert_eq!(
                observed.as_slice(),
                &phases[..target],
                "trace drift at poll {target}"
            );
        }
    }

    #[test]
    fn checkpointed_eval_preserves_cheap_metadata_error_precedence() {
        let (problem, objective, binding) = checkpoint_problem();
        let mut polls = 0usize;
        let missing = eval_borrowed_with_checkpoint(&problem, objective, &[], &mut |_phase| {
            polls += 1;
            Err(OptError::Cancelled)
        });
        assert!(matches!(
            missing,
            Err(OptError::BindingCount { got: 0, .. })
        ));
        assert_eq!(
            polls, 0,
            "binding count is checked before proportional work"
        );

        let short = &binding[..binding.len() - 1];
        let malformed =
            eval_borrowed_with_checkpoint(&problem, objective, &[short], &mut |_phase| {
                polls += 1;
                Err(OptError::Cancelled)
            });
        assert!(matches!(malformed, Err(OptError::BindingLen { .. })));
        assert_eq!(
            polls, 0,
            "binding length is checked before proportional work"
        );

        let unknown = eval_borrowed_with_checkpoint(
            &problem,
            NodeId(u32::MAX),
            &[binding.as_slice()],
            &mut |_phase| {
                polls += 1;
                Err(OptError::Cancelled)
            },
        );
        assert_eq!(unknown, Err(OptError::UnknownNode { id: u32::MAX }));
        assert_eq!(polls, 0, "unknown node is checked before proportional work");
    }

    #[test]
    fn binding_envelope_checkpoint_stride_is_pinned() {
        let mut builder = ProblemBuilder::new();
        let mut first = None;
        for index in 0..=EVAL_CHECKPOINT_STRIDE {
            let variable = builder
                .var(&format!("x{index}"), Manifold::Rn { dim: 1 }, Dims::NONE)
                .expect("bounded variable");
            first.get_or_insert(variable);
        }
        let point = builder
            .var_ref(first.expect("at least one variable"))
            .expect("point");
        let root = builder.component(point, 0).expect("scalar root");
        let problem = builder.finish();
        let owned = vec![vec![0.0]; EVAL_CHECKPOINT_STRIDE + 1];
        let bindings: Vec<&[f64]> = owned.iter().map(Vec::as_slice).collect();
        let mut phases = Vec::new();
        let value = eval_borrowed_with_checkpoint(&problem, root, &bindings, &mut |phase| {
            phases.push(phase);
            Ok(())
        })
        .expect("wide binding frame");
        assert_eq!(value, Value::S(0.0));
        assert_eq!(
            phases
                .iter()
                .filter(|phase| **phase == EvalPhase::BindingEnvelope)
                .count(),
            3,
            "indices 0 and 256 plus the terminal boundary must poll"
        );

        let mut envelope_polls = 0usize;
        let cancelled = eval_borrowed_with_checkpoint(&problem, root, &bindings, &mut |phase| {
            if phase == EvalPhase::BindingEnvelope {
                envelope_polls += 1;
                if envelope_polls == 2 {
                    return Err(OptError::Cancelled);
                }
            }
            Ok(())
        });
        assert_eq!(cancelled, Err(OptError::Cancelled));
        assert_eq!(envelope_polls, 2);
    }

    #[test]
    fn descend_ir_routes_every_evaluator_poll_through_cx() {
        let (problem, _objective, binding) = checkpoint_problem();
        let options = DescentOptions {
            steps: 0,
            ..DescentOptions::default()
        };
        let baseline_gate = fs_exec::CancelGate::new_clock_free();
        let baseline_phases = std::cell::RefCell::new(Vec::new());
        let baseline = with_test_cx(&baseline_gate, |cx| {
            descend_ir_with_eval_checkpoint(&problem, &binding, options, cx, &|phase| {
                baseline_phases.borrow_mut().push(phase);
                Ok(())
            })
            .expect("baseline IR descent")
        });
        let phases = baseline_phases.into_inner();
        assert!(!phases.is_empty());
        assert_eq!(
            phases.len(),
            100,
            "f0 must traverse one complete evaluator trace"
        );
        assert_eq!(baseline.evals, 1);
        assert_eq!(baseline.steps_taken, 0);
        assert_eq!(baseline.stop, DescentStop::StepLimit);
        assert!(!baseline.budget_stopped);
        assert_eq!(baseline.f0.to_bits(), 0.0f64.to_bits());
        assert_eq!(baseline.f_final.to_bits(), baseline.f0.to_bits());
        assert_eq!(
            baseline
                .x
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            binding
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            "steps=0 must return the exact input point bits"
        );
        assert_eq!(baseline.x[0].to_bits(), (-0.0f64).to_bits());

        let replay_gate = fs_exec::CancelGate::new_clock_free();
        let replay_phases = std::cell::RefCell::new(Vec::new());
        let replay = with_test_cx(&replay_gate, |cx| {
            descend_ir_with_eval_checkpoint(&problem, &binding, options, cx, &|phase| {
                replay_phases.borrow_mut().push(phase);
                Ok(())
            })
            .expect("replayed IR descent")
        });
        assert_eq!(replay_phases.into_inner(), phases);
        assert_eq!(
            replay
                .x
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            baseline
                .x
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(replay.f0.to_bits(), baseline.f0.to_bits());
        assert_eq!(replay.f_final.to_bits(), baseline.f_final.to_bits());
        assert_eq!(replay.evals, baseline.evals);
        assert_eq!(replay.steps_taken, baseline.steps_taken);
        assert_eq!(replay.stop, baseline.stop);
        assert_eq!(replay.budget_stopped, baseline.budget_stopped);
        assert_eq!(replay.work_upper_bound, baseline.work_upper_bound);
        assert_eq!(
            replay.workspace_upper_bound_bytes,
            baseline.workspace_upper_bound_bytes
        );

        for target in 1..=phases.len() {
            let gate = fs_exec::CancelGate::new_clock_free();
            let calls = std::cell::Cell::new(0usize);
            let observed = std::cell::RefCell::new(Vec::new());
            with_test_cx(&gate, |cx| {
                let result =
                    descend_ir_with_eval_checkpoint(&problem, &binding, options, cx, &|phase| {
                        let next = calls.get() + 1;
                        calls.set(next);
                        observed.borrow_mut().push(phase);
                        if next == target {
                            gate.request();
                        }
                        descent_checkpoint(cx)
                    });
                assert!(
                    matches!(result, Err(OptError::Cancelled)),
                    "target poll {target}"
                );
            });
            assert!(gate.is_requested());
            assert_eq!(calls.get(), target, "work continued after poll {target}");
            let observed = observed.into_inner();
            assert_eq!(
                observed.as_slice(),
                &phases[..target],
                "IR descent trace drift at poll {target}"
            );
        }
    }

    #[test]
    fn runtime_vector_reservations_refuse_capacity_overflow_without_partial_output() {
        let len = usize::MAX;
        let error = try_vec_capacity::<u64>("test/capacity-overflow", Some(NodeId(17)), None, len)
            .expect_err("an address-space-sized u64 vector cannot be reserved");
        assert_eq!(
            error,
            OptError::RuntimeAllocationRefused {
                path: "test/capacity-overflow",
                node: Some(17),
                variable: None,
                elements: allocation_len(len),
                element_bytes: u64::try_from(core::mem::size_of::<u64>())
                    .expect("u64 layout fits the diagnostic domain"),
            }
        );
        assert!(
            error
                .to_string()
                .contains("no partial result was published"),
            "the refusal must teach the atomic-publication boundary"
        );
    }

    #[test]
    fn fallible_runtime_vector_builders_preserve_arithmetic_bits() {
        let left = [1.0, -0.0, f64::MIN_POSITIVE];
        let right = [2.0, 0.0, -f64::MIN_POSITIVE];
        let node = Some(NodeId(3));
        let mut checkpoint = no_eval_checkpoint;
        let cloned = try_clone_vector("test/clone", node, Some(VarId(2)), &left, &mut checkpoint)
            .expect("bounded clone");
        let negated = try_map_vector("test/map", node, &left, |value| -value, &mut checkpoint)
            .expect("bounded map");
        let added = try_zip_vectors(
            "test/zip",
            node,
            &left,
            &right,
            |a, b| a + b,
            &mut checkpoint,
        )
        .expect("bounded zip");

        assert_eq!(
            cloned
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            left.iter().map(|value| value.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(
            negated
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            left.iter()
                .map(|value| (-*value).to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            added
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            left.iter()
                .zip(&right)
                .map(|(a, b)| (*a + *b).to_bits())
                .collect::<Vec<_>>()
        );
    }

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
    fn high_p_stiefel_work_envelope_dominates_scalar_visit_lower_bound() {
        let n = 64u32;
        let p = 64u32;
        let point_dim = n * p;
        let envelope = descent_envelope(
            Manifold::Stiefel { n, p },
            point_dim,
            point_dim,
            0,
            EvalLimit::Unlimited,
        )
        .expect("small exact envelope");

        // One current Stiefel retraction necessarily visits at least
        // N*(4.5*p + 10.5) scalar slots across input/output Gram scans,
        // deterministic QR dot/update/revalidation, and copies. Compare
        // doubled integers so the tripwire itself has no rounding.
        let lower_bound_twice = u64::from(point_dim) * (9 * u64::from(p) + 21);
        assert!(
            envelope.owned_work * 2 >= lower_bound_twice,
            "high-p Stiefel plans must never understate scalar visits"
        );
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
