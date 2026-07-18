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
/// Scale-free tolerance for validating a parameter vector as tangent.
const MANIFOLD_TANGENCY_TOL: f64 = 64.0 * MANIFOLD_MEMBERSHIP_TOL;
/// Bounded Richardson refinements for the near-identity Stiefel Sylvester
/// equation used by the embedded-metric tangent projection.
///
/// Admission gives `p <= n` and `n * p <= u32::MAX`, hence `p <= 65_535`.
/// Every entry of `E = X^T X - I` is at most `1e-10`, so the exact-real
/// max-entry residual recurrence
///
/// `R_next = -(E R + R E) / 2`
///
/// contracts by at most `p * 1e-10 < 2^-17`. One hundred twenty-eight
/// contractions take even an initially finite `f64::MAX` residual below the
/// minimum subnormal magnitude. Floating-point roundoff can still prevent that
/// ideal contraction, so the measured compensated residual below remains the
/// authority and the operation fails closed at this explicit bound.
const STIEFEL_PROJECTION_REFINEMENT_LIMIT: usize = 128;
/// Candidates below this squared norm are numerically rank-deficient;
/// normalizing them would fabricate a direction.
const RETRACTION_MIN_NORM_SQ: f64 = 1e-24;

/// Add one term with deterministic Neumaier compensation.
///
/// Projection and tangent checks scale unconstrained vectors to unit max norm;
/// point-membership checks apply it to the declared unit/orthonormal domain.
/// A non-finite accumulator is therefore an explicit domain failure, not an
/// intermediate projection overflow that later cancellation should repair.
fn compensated_add(sum: &mut f64, correction: &mut f64, term: f64) {
    let next = *sum + term;
    if !next.is_finite() {
        *sum = next;
        *correction = 0.0;
        return;
    }
    let roundoff = if sum.abs() >= term.abs() {
        (*sum - next) + term
    } else {
        (term - next) + *sum
    };
    let corrected = *correction + roundoff;
    if !corrected.is_finite() {
        *sum = corrected;
        *correction = 0.0;
        return;
    }
    *sum = next;
    *correction = corrected;
}

fn compensated_sum(terms: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for term in terms {
        compensated_add(&mut sum, &mut correction, term);
    }
    sum + correction
}

fn max_abs(values: &[f64]) -> f64 {
    values
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max)
}

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

/// Select the lexicographically positive representative of an SO(3) antipodal
/// pair and give every exact-zero lane one deterministic bit spelling.
fn canonicalize_so3_representative(quaternion: &mut [f64; 4]) {
    let negate = quaternion
        .iter()
        .copied()
        .find(|component| *component != 0.0)
        .is_some_and(f64::is_sign_negative);
    if negate {
        for component in quaternion.iter_mut() {
            *component = -*component;
        }
    }
    for component in quaternion.iter_mut() {
        if *component == 0.0 {
            *component = 0.0;
        }
    }
}

/// An evaluated node value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Scalar.
    S(f64),
    /// Vector.
    V(Vec<f64>),
}

/// One point on a manifold retraction curve and its optimizer-parameter
/// velocity at that point.
///
/// `point` always has [`Manifold::point_dim`] entries and `velocity` always has
/// [`Manifold::param_dim`] entries. For Rn, Sphere, and Stiefel the current
/// parameter representation is ambient, so `velocity` is also the derivative
/// of the stored point coordinates. SO(3) is deliberately different: its point
/// is a four-coordinate quaternion while `velocity` is the constant
/// three-coordinate right/body Lie velocity of the declared curve.
#[derive(Debug, Clone, PartialEq)]
pub struct RetractionCurve {
    /// Retraction output at the requested scalar curve parameter.
    pub point: Vec<f64>,
    /// Curve velocity in the manifold's optimizer-parameter representation.
    pub velocity: Vec<f64>,
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

fn try_clone_runtime_slice<F>(
    path: &'static str,
    values: &[f64],
    checkpoint: &mut F,
) -> Result<Vec<f64>, OptError>
where
    F: FnMut() -> Result<(), OptError>,
{
    let mut output = try_vec_capacity(path, None, None, values.len())?;
    for (component, value) in values.iter().enumerate() {
        checkpoint_retraction_work(component, checkpoint)?;
        output.push(*value);
    }
    checkpoint()?;
    Ok(output)
}

fn try_filled_runtime_vector<F>(
    path: &'static str,
    len: usize,
    value: f64,
    checkpoint: &mut F,
) -> Result<Vec<f64>, OptError>
where
    F: FnMut() -> Result<(), OptError>,
{
    let mut output = try_vec_capacity(path, None, None, len)?;
    for component in 0..len {
        checkpoint_retraction_work(component, checkpoint)?;
        output.push(value);
    }
    checkpoint()?;
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
                let mut norm_correction = 0.0;
                for (component, value) in x.iter().enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    compensated_add(&mut norm_sq, &mut norm_correction, value * value);
                }
                norm_sq += norm_correction;
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
                        let mut dot_correction = 0.0;
                        for row in 0..n {
                            checkpoint_retraction_work(row, checkpoint)?;
                            compensated_add(
                                &mut dot,
                                &mut dot_correction,
                                x[column * n + row] * x[against * n + row],
                            );
                        }
                        dot += dot_correction;
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

    fn validate_operation_point(&self, point: &[f64]) -> Result<(), OptError> {
        self.validate(&AdmissionCaps::default())?;
        let point_dim = match self.point_dim() {
            Some(point_dim) => point_dim,
            None => {
                return Err(OptError::ManifoldInvalid {
                    what: try_owned_diagnostic(
                        "manifold-ops/point-diagnostic",
                        None,
                        None,
                        "manifold operation has no representable point dimension",
                    )?,
                });
            }
        };
        if point.len() as u64 != u64::from(point_dim) {
            return Err(OptError::RetractionLen {
                input: "manifold operation point",
                expected: point_dim,
                got: point.len() as u64,
            });
        }
        let mut checkpoint = || Ok(());
        Self::validate_retraction_finite_with_checkpoint(
            "manifold operation point",
            point,
            &mut checkpoint,
        )?;
        self.validate_point_domain_with_checkpoint(point, &mut checkpoint)
    }

    fn validate_parameter_payload(
        &self,
        input: &'static str,
        parameter: &[f64],
    ) -> Result<(), OptError> {
        let param_dim = match self.param_dim() {
            Some(param_dim) => param_dim,
            None => {
                return Err(OptError::ManifoldInvalid {
                    what: try_owned_diagnostic(
                        "manifold-ops/parameter-diagnostic",
                        None,
                        None,
                        "manifold operation has no representable parameter dimension",
                    )?,
                });
            }
        };
        if parameter.len() as u64 != u64::from(param_dim) {
            return Err(OptError::RetractionLen {
                input,
                expected: param_dim,
                got: parameter.len() as u64,
            });
        }
        Self::validate_retraction_finite_with_checkpoint(input, parameter, &mut || Ok(()))
    }

    fn validate_ambient_payload(
        &self,
        input: &'static str,
        ambient: &[f64],
    ) -> Result<(), OptError> {
        let point_dim = match self.point_dim() {
            Some(point_dim) => point_dim,
            None => {
                return Err(OptError::ManifoldInvalid {
                    what: try_owned_diagnostic(
                        "manifold-ops/ambient-diagnostic",
                        None,
                        None,
                        "manifold operation has no representable point dimension",
                    )?,
                });
            }
        };
        if ambient.len() as u64 != u64::from(point_dim) {
            return Err(OptError::RetractionLen {
                input,
                expected: point_dim,
                got: ambient.len() as u64,
            });
        }
        Self::validate_retraction_finite_with_checkpoint(input, ambient, &mut || Ok(()))
    }

    fn parameter_tangency_residual(&self, point: &[f64], parameter: &[f64]) -> f64 {
        let scale = max_abs(parameter);
        if scale == 0.0 {
            return 0.0;
        }
        match *self {
            Manifold::Rn { .. } | Manifold::So3 => 0.0,
            Manifold::Sphere { .. } => compensated_sum(
                point
                    .iter()
                    .zip(parameter)
                    .map(|(coordinate, tangent)| coordinate * (tangent / scale)),
            )
            .abs(),
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                let mut residual = 0.0_f64;
                for row in 0..p {
                    for column in 0..p {
                        let left = compensated_sum((0..n).map(|index| {
                            point[row * n + index] * (parameter[column * n + index] / scale)
                        }));
                        let right = compensated_sum((0..n).map(|index| {
                            (parameter[row * n + index] / scale) * point[column * n + index]
                        }));
                        let symmetric = compensated_sum([left, right]);
                        residual = residual.max(symmetric.abs());
                    }
                }
                residual
            }
        }
    }

    fn checked_manifold_work_len(
        &self,
        left: usize,
        right: usize,
        rule: &'static str,
    ) -> Result<usize, OptError> {
        left.checked_mul(right)
            .ok_or_else(|| self.domain_error(rule, None, (left as f64) * (right as f64)))
    }

    fn normalized_ambient_component(
        &self,
        value: f64,
        ambient_scale: f64,
        component: usize,
    ) -> Result<f64, OptError> {
        let normalized = value / ambient_scale;
        if value != 0.0 && normalized == 0.0 {
            return Err(self.domain_error(
                "ambient gradient component is lost during max-abs normalization",
                Some((component as u32, 0)),
                value,
            ));
        }
        Ok(normalized)
    }

    fn rescale_parameter_gradient(
        &self,
        output: &mut [f64],
        ambient_scale: f64,
    ) -> Result<(), OptError> {
        for (component, value) in output.iter_mut().enumerate() {
            let rescaled = *value * ambient_scale;
            if !rescaled.is_finite() {
                return Err(OptError::RetractionNonFinite {
                    input: "parameter gradient output",
                    component: component as u32,
                    bits: rescaled.to_bits(),
                });
            }
            if *value != 0.0 && rescaled == 0.0 {
                return Err(self.domain_error(
                    "parameter gradient component is lost while restoring caller scale",
                    Some((component as u32, 0)),
                    *value,
                ));
            }
            *value = rescaled;
        }
        Ok(())
    }

    fn reconstruct_sphere_tangent_pivot(
        &self,
        point: &[f64],
        tangent: &mut [f64],
    ) -> Result<(), OptError> {
        let mut pivot = 0_usize;
        for component in 1..point.len() {
            if point[component].abs() > point[pivot].abs() {
                pivot = component;
            }
        }
        let pivot_coordinate = point[pivot];
        if pivot_coordinate == 0.0 {
            return Err(self.domain_error(
                "Sphere tangent reconstruction requires a nonzero point pivot",
                Some((pivot as u32, 0)),
                pivot_coordinate,
            ));
        }
        let other_contraction =
            compensated_sum(point.iter().zip(tangent.iter()).enumerate().filter_map(
                |(component, (coordinate, value))| {
                    (component != pivot).then_some(coordinate * value)
                },
            ));
        let reconstructed = -other_contraction / pivot_coordinate;
        if !reconstructed.is_finite() {
            return Err(self.domain_error(
                "Sphere tangent pivot reconstruction must remain finite",
                Some((pivot as u32, 0)),
                reconstructed,
            ));
        }
        tangent[pivot] = reconstructed;
        let residual = self.parameter_tangency_residual(point, tangent);
        if !residual.is_finite() || residual > MANIFOLD_TANGENCY_TOL {
            return Err(self.domain_error(
                "Sphere tangent pivot reconstruction did not reach the tangent space",
                Some((pivot as u32, 0)),
                residual,
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)] // explicit Sylvester residual/refusal phases stay auditable
    fn refine_stiefel_projection(
        &self,
        point: &[f64],
        projected: &mut [f64],
    ) -> Result<(), OptError> {
        let Manifold::Stiefel { n, p } = *self else {
            unreachable!("Stiefel projection refinement is Stiefel-only");
        };
        let (n, p) = (n as usize, p as usize);
        let residual_storage = self.checked_manifold_work_len(
            p,
            p,
            "Stiefel projection Sylvester workspace must be representable",
        )?;
        let mut residual = try_vec_capacity(
            "manifold-ops/stiefel-gradient-sylvester",
            None,
            None,
            residual_storage,
        )?;

        // For an admitted frame M = X^T X differs from I by at most the
        // membership tolerance entrywise. Starting from the usual
        // X-orthonormal formula, each update V <- V - X R/2 applies a
        // deterministic Richardson correction to the Sylvester equation
        // M S + S M = X^T G + G^T X. In exact arithmetic the residual obeys
        // R_next = -(E R + R E)/2 for E=M-I. The descriptor constraints imply
        // p^2 <= n*p <= u32::MAX, and therefore
        // ||R_next||_max <= p*1e-10*||R||_max < 2^-17*||R||_max. The explicit
        // 128-step bound is consequently enough to drive any initially finite
        // f64 residual below the minimum subnormal in the ideal recurrence.
        // Floating arithmetic is accepted only by the measured scale-free
        // postcondition; stagnation or other unresolved roundoff fails closed.
        for _ in 0..STIEFEL_PROJECTION_REFINEMENT_LIMIT {
            if projected.iter().all(|value| *value == 0.0) {
                return Ok(());
            }
            let tangency_residual = self.parameter_tangency_residual(point, projected);
            if tangency_residual <= MANIFOLD_TANGENCY_TOL {
                return Ok(());
            }

            residual.clear();
            for left in 0..p {
                for right in 0..p {
                    let first = compensated_sum(
                        (0..n).map(|row| point[left * n + row] * projected[right * n + row]),
                    );
                    let second = compensated_sum(
                        (0..n).map(|row| projected[left * n + row] * point[right * n + row]),
                    );
                    let symmetric = compensated_sum([first, second]);
                    if !symmetric.is_finite() {
                        return Err(self.domain_error(
                            "Stiefel projection Sylvester residual must be finite",
                            Some((left as u32, right as u32)),
                            symmetric,
                        ));
                    }
                    residual.push(symmetric);
                }
            }

            for column in 0..p {
                for row in 0..n {
                    let correction = 0.5
                        * compensated_sum(
                            (0..p)
                                .map(|basis| point[basis * n + row] * residual[basis * p + column]),
                        );
                    projected[column * n + row] -= correction;
                }
            }
            Self::validate_retraction_finite_with_checkpoint(
                "Stiefel projection refinement",
                projected,
                &mut || Ok(()),
            )?;
            let relative_residual = self.parameter_tangency_residual(point, projected);
            if relative_residual <= MANIFOLD_TANGENCY_TOL {
                return Ok(());
            }
        }

        let tangency_residual = self.parameter_tangency_residual(point, projected);
        if !tangency_residual.is_finite() || tangency_residual > MANIFOLD_TANGENCY_TOL {
            return Err(self.domain_error(
                "Stiefel projection Sylvester refinement did not reach the tangent space",
                None,
                tangency_residual,
            ));
        }
        Ok(())
    }

    /// Validate a vector in the current retraction-parameter representation as
    /// tangent at `point`.
    ///
    /// Rn parameters and SO(3) right/body Lie coordinates are unconstrained.
    /// Sphere parameters must be orthogonal to the stored point. Stiefel
    /// ambient parameters must satisfy `X^T V + V^T X = 0` under the embedded
    /// Frobenius metric.
    ///
    /// # Errors
    /// Returns the typed manifold/retraction error family for malformed
    /// descriptors, lengths, non-finite values, off-manifold points, or a
    /// finite parameter vector outside the declared tangent space.
    pub fn validate_parameter_tangent(
        &self,
        point: &[f64],
        parameter: &[f64],
    ) -> Result<(), OptError> {
        self.validate_operation_point(point)?;
        self.validate_parameter_payload("manifold tangent parameter", parameter)?;
        let residual = self.parameter_tangency_residual(point, parameter);
        if !residual.is_finite() || residual > MANIFOLD_TANGENCY_TOL {
            return Err(self.domain_error(
                "parameter vector must belong to the point tangent space",
                None,
                residual,
            ));
        }
        Ok(())
    }

    /// Pull an ambient point-storage gradient back to the current retraction
    /// parameterization.
    ///
    /// Rn is the identity. Sphere and Stiefel use embedded-metric tangent
    /// projection. SO(3) applies the differential of the declared
    /// right-composed quaternion exponential, producing three body-frame
    /// components from a four-coordinate ambient quaternion gradient. Every
    /// reduction-bearing branch projects a max-abs-scaled ambient gradient,
    /// then checks the component-attributed rescale back to caller units.
    ///
    /// # Errors
    /// Returns the typed manifold/retraction error family for invalid points,
    /// lengths, domains, non-finite inputs, or non-finite projected output.
    #[must_use]
    #[allow(clippy::too_many_lines)] // one explicit scale/project/rescale branch per manifold
    pub fn parameter_gradient(
        &self,
        point: &[f64],
        ambient_gradient: &[f64],
    ) -> Result<Vec<f64>, OptError> {
        self.validate_operation_point(point)?;
        self.validate_ambient_payload("ambient manifold gradient", ambient_gradient)?;
        let mut output = try_vec_capacity(
            "manifold-ops/parameter-gradient",
            None,
            None,
            usize::try_from(self.param_dim().expect("validated parameter dimension"))
                .expect("u32 parameter dimension fits usize"),
        )?;
        match *self {
            Manifold::Rn { .. } => output.extend_from_slice(ambient_gradient),
            Manifold::Sphere { .. } => {
                let ambient_scale = max_abs(ambient_gradient);
                if ambient_scale == 0.0 {
                    output.resize(ambient_gradient.len(), 0.0);
                } else {
                    for (component, value) in ambient_gradient.iter().copied().enumerate() {
                        output.push(self.normalized_ambient_component(
                            value,
                            ambient_scale,
                            component,
                        )?);
                    }
                    let norm_sq =
                        compensated_sum(point.iter().map(|coordinate| coordinate * coordinate));
                    let normal = compensated_sum(
                        point
                            .iter()
                            .zip(&output)
                            .map(|(coordinate, gradient)| coordinate * gradient),
                    ) / norm_sq;
                    for (coordinate, gradient) in point.iter().zip(&mut output) {
                        *gradient = (-normal).mul_add(*coordinate, *gradient);
                    }

                    // Membership deliberately admits tiny norm error. Remove
                    // the resulting roundoff-normal component using the same
                    // x.x geometry, rather than rejecting an otherwise valid
                    // radial gradient after scale-free tangent validation.
                    for _ in 0..2 {
                        let residual = compensated_sum(
                            point
                                .iter()
                                .zip(&output)
                                .map(|(coordinate, gradient)| coordinate * gradient),
                        ) / norm_sq;
                        if residual == 0.0 {
                            break;
                        }
                        for (coordinate, gradient) in point.iter().zip(&mut output) {
                            *gradient = (-residual).mul_add(*coordinate, *gradient);
                        }
                    }
                    // Membership permits a small norm defect, so cancellation
                    // in the normal correction need not make x.v bit-exact.
                    // Reconstruct the lane with the largest |x_i| from all
                    // others. That pivot minimizes division amplification,
                    // preserves every non-pivot tangent component, and makes a
                    // purely axial normal projection exactly zero without any
                    // magnitude-based whole-vector canonicalization.
                    self.reconstruct_sphere_tangent_pivot(point, &mut output)?;
                    self.rescale_parameter_gradient(&mut output, ambient_scale)?;
                }
            }
            Manifold::So3 => {
                let [w, x, y, z] = [point[0], point[1], point[2], point[3]];
                let ambient_scale = max_abs(ambient_gradient);
                if ambient_scale == 0.0 {
                    output.extend([0.0; 3]);
                } else {
                    let [gw, gx, gy, gz] = [
                        self.normalized_ambient_component(ambient_gradient[0], ambient_scale, 0)?,
                        self.normalized_ambient_component(ambient_gradient[1], ambient_scale, 1)?,
                        self.normalized_ambient_component(ambient_gradient[2], ambient_scale, 2)?,
                        self.normalized_ambient_component(ambient_gradient[3], ambient_scale, 3)?,
                    ];
                    output.extend([
                        0.5 * compensated_sum([-x * gw, w * gx, z * gy, -y * gz]),
                        0.5 * compensated_sum([-y * gw, -z * gx, w * gy, x * gz]),
                        0.5 * compensated_sum([-z * gw, y * gx, -x * gy, w * gz]),
                    ]);
                    self.rescale_parameter_gradient(&mut output, ambient_scale)?;
                }
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                let ambient_scale = max_abs(ambient_gradient);
                if ambient_scale == 0.0 {
                    output.resize(ambient_gradient.len(), 0.0);
                } else {
                    for (component, value) in ambient_gradient.iter().copied().enumerate() {
                        output.push(self.normalized_ambient_component(
                            value,
                            ambient_scale,
                            component,
                        )?);
                    }
                    let gram_storage = self.checked_manifold_work_len(
                        p,
                        p,
                        "Stiefel projection Gram workspace must be representable",
                    )?;
                    let mut gram = try_vec_capacity(
                        "manifold-ops/stiefel-gradient-gram",
                        None,
                        None,
                        gram_storage,
                    )?;
                    for left in 0..p {
                        for right in 0..p {
                            let entry = compensated_sum(
                                (0..n).map(|row| point[left * n + row] * output[right * n + row]),
                            );
                            if !entry.is_finite() {
                                return Err(self.domain_error(
                                    "Stiefel projection Gram entry must be finite",
                                    Some((left as u32, right as u32)),
                                    entry,
                                ));
                            }
                            gram.push(entry);
                        }
                    }
                    for column in 0..p {
                        for row in 0..n {
                            let correction = compensated_sum((0..p).map(|basis| {
                                let symmetric = 0.5
                                    * compensated_sum([
                                        gram[basis * p + column],
                                        gram[column * p + basis],
                                    ]);
                                point[basis * n + row] * symmetric
                            }));
                            output[column * n + row] -= correction;
                        }
                    }
                    Self::validate_retraction_finite_with_checkpoint(
                        "scaled Stiefel parameter gradient",
                        &output,
                        &mut || Ok(()),
                    )?;
                    self.refine_stiefel_projection(point, &mut output)?;
                    self.rescale_parameter_gradient(&mut output, ambient_scale)?;
                }
            }
        }
        Self::validate_retraction_finite_with_checkpoint(
            "parameter gradient output",
            &output,
            &mut || Ok(()),
        )?;
        self.validate_parameter_tangent(point, &output)?;
        Ok(output)
    }

    #[allow(clippy::too_many_lines)] // differentiated audit trail mirrors each QR phase
    fn stiefel_qr_differential(
        &self,
        point: &[f64],
        step: &[f64],
        variation: &[f64],
    ) -> Result<RetractionCurve, OptError> {
        let Manifold::Stiefel { n, p } = *self else {
            unreachable!("QR differential is Stiefel-only");
        };
        let (n, p) = (n as usize, p as usize);
        let mut columns = try_vec_capacity("manifold-ops/stiefel-qr-columns", None, None, p)?;
        let mut derivatives =
            try_vec_capacity("manifold-ops/stiefel-qr-derivatives", None, None, p)?;
        for column in 0..p {
            let mut candidate = try_vec_capacity("manifold-ops/stiefel-qr-column", None, None, n)?;
            let mut derivative =
                try_vec_capacity("manifold-ops/stiefel-qr-derivative", None, None, n)?;
            for row in 0..n {
                candidate.push(point[column * n + row] + step[column * n + row]);
                derivative.push(variation[column * n + row]);
            }
            columns.push(candidate);
            derivatives.push(derivative);
        }

        for current_index in 0..p {
            for prior_index in 0..current_index {
                let projection: f64 = columns[current_index]
                    .iter()
                    .zip(&columns[prior_index])
                    .map(|(current, prior)| current * prior)
                    .sum();
                let projection_velocity: f64 = derivatives[current_index]
                    .iter()
                    .zip(&columns[prior_index])
                    .map(|(current, prior)| current * prior)
                    .sum::<f64>()
                    + columns[current_index]
                        .iter()
                        .zip(&derivatives[prior_index])
                        .map(|(current, prior)| current * prior)
                        .sum::<f64>();
                let (prior_columns, current_columns) = columns.split_at_mut(current_index);
                let prior = &prior_columns[prior_index];
                let current = &mut current_columns[0];
                let (prior_derivatives, current_derivatives) =
                    derivatives.split_at_mut(current_index);
                let prior_derivative = &prior_derivatives[prior_index];
                let current_derivative = &mut current_derivatives[0];
                for row in 0..n {
                    current_derivative[row] -=
                        projection_velocity * prior[row] + projection * prior_derivative[row];
                    current[row] -= projection * prior[row];
                }
            }
            let current = &mut columns[current_index];
            let current_derivative = &mut derivatives[current_index];
            let norm_sq: f64 = current.iter().map(|value| value * value).sum();
            if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                return Err(self.domain_error(
                    "candidate column is rank-deficient",
                    Some((current_index as u32, current_index as u32)),
                    norm_sq,
                ));
            }
            let norm = norm_sq.sqrt();
            let norm_velocity = current
                .iter()
                .zip(current_derivative.iter())
                .map(|(candidate, derivative)| candidate * derivative)
                .sum::<f64>()
                / norm;
            for row in 0..n {
                let normalized = current[row] / norm;
                current_derivative[row] =
                    (current_derivative[row] - normalized * norm_velocity) / norm;
                current[row] = normalized;
            }
        }

        let storage = self.checked_manifold_work_len(
            n,
            p,
            "Stiefel differentiated QR storage must be representable",
        )?;
        let mut output = try_vec_capacity("manifold-ops/stiefel-qr-output", None, None, storage)?;
        let mut velocity =
            try_vec_capacity("manifold-ops/stiefel-qr-velocity", None, None, storage)?;
        for column in 0..p {
            output.extend_from_slice(&columns[column]);
            velocity.extend_from_slice(&derivatives[column]);
        }
        Self::validate_retraction_finite_with_checkpoint(
            "Stiefel curve output",
            &output,
            &mut || Ok(()),
        )?;
        Self::validate_retraction_finite_with_checkpoint(
            "Stiefel curve velocity",
            &velocity,
            &mut || Ok(()),
        )?;
        Ok(RetractionCurve {
            point: output,
            velocity,
        })
    }

    fn require_same_curve_point(
        &self,
        declared: &[f64],
        reconstructed: &[f64],
        rule: &'static str,
    ) -> Result<(), OptError> {
        if declared.len() != reconstructed.len() {
            return Err(OptError::RetractionLen {
                input: rule,
                expected: u32::try_from(reconstructed.len()).unwrap_or(u32::MAX),
                got: declared.len() as u64,
            });
        }
        if declared
            .iter()
            .zip(reconstructed)
            .all(|(left, right)| left.to_bits() == right.to_bits())
        {
            return Ok(());
        }
        let mismatch = declared
            .iter()
            .zip(reconstructed)
            .map(|(left, right)| (left - right).abs())
            .fold(0.0_f64, f64::max);
        Err(self.domain_error(rule, None, mismatch))
    }

    /// Evaluate the authoritative retraction curve `R_point(alpha * direction)`
    /// and its optimizer-parameter velocity at the landed point.
    ///
    /// The SO(3) branch returns the constant three-coordinate right/body Lie
    /// velocity, distinct from the four-coordinate derivative of its stored
    /// quaternion point. The Stiefel branch differentiates the same
    /// deterministic modified-Gram-Schmidt QR program used by
    /// [`Manifold::retract`].
    ///
    /// # Errors
    /// Returns the typed manifold/retraction error family for invalid points,
    /// parameter directions, scalar parameters, singular curves, or non-finite
    /// outputs.
    #[must_use]
    #[allow(clippy::too_many_lines)] // one explicit velocity branch per manifold
    pub fn retract_curve(
        &self,
        point: &[f64],
        direction: &[f64],
        alpha: f64,
    ) -> Result<RetractionCurve, OptError> {
        self.validate_operation_point(point)?;
        self.validate_parameter_payload("retraction curve direction", direction)?;
        if !alpha.is_finite() {
            return Err(OptError::RetractionNonFinite {
                input: "retraction curve alpha",
                component: 0,
                bits: alpha.to_bits(),
            });
        }
        let mut step = try_vec_capacity(
            "manifold-ops/retraction-curve-step",
            None,
            None,
            direction.len(),
        )?;
        step.extend(direction.iter().map(|value| alpha * value));
        Self::validate_retraction_finite_with_checkpoint(
            "retraction curve step",
            &step,
            &mut || Ok(()),
        )?;
        let landed = self.retract(point, &step)?;
        let velocity = match *self {
            Manifold::Rn { .. } => try_clone_runtime_slice(
                "manifold-ops/rn-curve-velocity",
                direction,
                &mut || Ok(()),
            )?,
            Manifold::Sphere { .. } => {
                let mut moved = try_vec_capacity(
                    "manifold-ops/sphere-curve-candidate",
                    None,
                    None,
                    point.len(),
                )?;
                for (coordinate, increment) in point.iter().zip(&step) {
                    // Preserve the public retraction's multiply-then-add
                    // operation order: `step = alpha * direction`, followed
                    // by `candidate = point + step`.
                    moved.push(coordinate + increment);
                }
                let norm = moved.iter().map(|value| value * value).sum::<f64>().sqrt();
                let radial: f64 = landed
                    .iter()
                    .zip(direction)
                    .map(|(coordinate, tangent)| coordinate * tangent)
                    .sum();
                let mut velocity = try_vec_capacity(
                    "manifold-ops/sphere-curve-velocity",
                    None,
                    None,
                    direction.len(),
                )?;
                for (tangent, coordinate) in direction.iter().zip(&landed) {
                    velocity.push(((-radial).mul_add(*coordinate, *tangent)) / norm);
                }
                velocity
            }
            Manifold::So3 => {
                try_clone_runtime_slice("manifold-ops/so3-curve-velocity", direction, &mut || {
                    Ok(())
                })?
            }
            Manifold::Stiefel { .. } => {
                let differentiated = self.stiefel_qr_differential(point, &step, direction)?;
                self.require_same_curve_point(
                    &landed,
                    &differentiated.point,
                    "differentiated QR point must match the authoritative retraction",
                )?;
                differentiated.velocity
            }
        };
        Self::validate_retraction_finite_with_checkpoint(
            "retraction curve velocity",
            &velocity,
            &mut || Ok(()),
        )?;
        let tangency_residual = self.parameter_tangency_residual(&landed, &velocity);
        if !tangency_residual.is_finite() || tangency_residual > MANIFOLD_TANGENCY_TOL {
            return Err(self.domain_error(
                "retraction curve velocity must be tangent at the landed point",
                None,
                tangency_residual,
            ));
        }
        Ok(RetractionCurve {
            point: landed,
            velocity,
        })
    }

    /// Transport a tangent parameter vector along one authoritative retraction
    /// step.
    ///
    /// Rn is the identity. Sphere uses shortest-geodesic isometric transport.
    /// SO(3) keeps right/body components unchanged, an isometric trivialized
    /// transport (not a Levi-Civita-parallel claim). Stiefel applies the
    /// differential of the production QR retraction; it is a deterministic
    /// tangent transport but no metric-isometry or superlinear-BFGS claim is
    /// made.
    ///
    /// `to` must be the exact result of `self.retract(from, step)`, preventing
    /// persisted or caller-spliced endpoints from silently changing transport
    /// semantics. `step` follows the existing retraction-parameter contract:
    /// Sphere and Stiefel accept ambient-coordinate increments, so callers that
    /// need a textbook tangent-bundle curve must validate the step separately
    /// with [`Manifold::validate_parameter_tangent`].
    ///
    /// # Errors
    /// Returns the typed manifold/retraction error family for malformed points,
    /// step/vector dimensions, non-finite inputs, non-tangent vectors,
    /// antipodal Sphere endpoints, singular QR paths, or an incompatible
    /// destination.
    #[must_use]
    pub fn transport_parameter(
        &self,
        from: &[f64],
        step: &[f64],
        to: &[f64],
        vector: &[f64],
    ) -> Result<Vec<f64>, OptError> {
        self.validate_operation_point(from)?;
        self.validate_operation_point(to)?;
        self.validate_parameter_payload("transport retraction step", step)?;
        self.validate_parameter_tangent(from, vector)?;
        let expected_to = self.retract(from, step)?;
        self.require_same_curve_point(
            to,
            &expected_to,
            "transport destination must equal the authoritative retraction",
        )?;
        let output = match *self {
            Manifold::Rn { .. } | Manifold::So3 => {
                try_clone_runtime_slice("manifold-ops/trivialized-transport", vector, &mut || {
                    Ok(())
                })?
            }
            Manifold::Sphere { .. } => {
                let from_norm = from.iter().map(|value| value * value).sum::<f64>().sqrt();
                let to_norm = to.iter().map(|value| value * value).sum::<f64>().sqrt();
                let vector_scale = vector
                    .iter()
                    .map(|value| value.abs())
                    .fold(0.0_f64, f64::max);
                let denominator = 1.0
                    + from
                        .iter()
                        .zip(to)
                        .map(|(left, right)| (left / from_norm) * (right / to_norm))
                        .sum::<f64>();
                let dimension_roundoff = 64.0 * f64::EPSILON * from.len() as f64;
                let denominator_floor = (MANIFOLD_MEMBERSHIP_TOL + dimension_roundoff).sqrt();
                if !denominator.is_finite() || denominator <= denominator_floor {
                    return Err(self.domain_error(
                        "Sphere transport is undefined or ill-conditioned near antipodal points",
                        None,
                        denominator,
                    ));
                }
                let scaled_coefficient = if vector_scale == 0.0 {
                    0.0
                } else {
                    vector
                        .iter()
                        .zip(to)
                        .map(|(tangent, coordinate)| {
                            (tangent / vector_scale) * (coordinate / to_norm)
                        })
                        .sum::<f64>()
                        / denominator
                };
                let mut transported =
                    try_vec_capacity("manifold-ops/sphere-transport", None, None, vector.len())?;
                for (value, (old, new)) in vector.iter().zip(from.iter().zip(to)) {
                    let unit_sum = old / from_norm + new / to_norm;
                    let scaled_value = if vector_scale == 0.0 {
                        0.0
                    } else {
                        value / vector_scale
                    };
                    transported
                        .push(vector_scale * (-scaled_coefficient).mul_add(unit_sum, scaled_value));
                }

                let comparison_scale = vector
                    .iter()
                    .chain(&transported)
                    .map(|value| value.abs())
                    .fold(0.0_f64, f64::max);
                let isometry_residual = if comparison_scale == 0.0 {
                    0.0
                } else {
                    let before = vector
                        .iter()
                        .map(|value| {
                            let scaled = value / comparison_scale;
                            scaled * scaled
                        })
                        .sum::<f64>();
                    let after = transported
                        .iter()
                        .map(|value| {
                            let scaled = value / comparison_scale;
                            scaled * scaled
                        })
                        .sum::<f64>();
                    (before - after).abs() / before.max(after)
                };
                let isometry_tolerance = 8.0 * MANIFOLD_TANGENCY_TOL + dimension_roundoff;
                if !isometry_residual.is_finite() || isometry_residual > isometry_tolerance {
                    return Err(self.domain_error(
                        "Sphere transport must preserve tangent norm within conditioning tolerance",
                        None,
                        isometry_residual,
                    ));
                }
                transported
            }
            Manifold::Stiefel { .. } => {
                let differentiated = self.stiefel_qr_differential(from, step, vector)?;
                self.require_same_curve_point(
                    to,
                    &differentiated.point,
                    "transport QR point must match the authoritative retraction",
                )?;
                differentiated.velocity
            }
        };
        Self::validate_retraction_finite_with_checkpoint(
            "transported parameter vector",
            &output,
            &mut || Ok(()),
        )?;
        self.validate_parameter_tangent(to, &output)?;
        Ok(output)
    }

    /// Retract: move from `x` along parameter vector `t`, landing ON
    /// the manifold. Rn: translation. Sphere: normalize(x+t). SO(3):
    /// right-multiply by `exp(ω/2)` (unit quaternion), then choose the sign
    /// whose first exact nonzero `(w, x, y, z)` component is positive and
    /// spell every zero lane as `+0.0`. Stiefel: Gram-Schmidt of `X+T` (QR
    /// retraction). Raw point and parameter storage must exactly match this
    /// manifold and contain finite components. The base point must already
    /// belong to the manifold; zero-norm and rank-deficient candidates are
    /// refused rather than normalized into fabricated points.
    ///
    /// # Errors
    /// [`OptError::ManifoldInvalid`], [`OptError::RetractionLen`],
    /// [`OptError::RetractionNonFinite`], [`OptError::RetractionDomain`], or
    /// [`OptError::RuntimeAllocationRefused`].
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
        let point_dim = match self.point_dim() {
            Some(point_dim) => point_dim,
            None => {
                return Err(OptError::ManifoldInvalid {
                    what: try_owned_diagnostic(
                        "retract/manifold-diagnostic",
                        None,
                        None,
                        "retraction manifold has no representable point dimension",
                    )?,
                });
            }
        };
        if x.len() as u64 != u64::from(point_dim) {
            return Err(OptError::RetractionLen {
                input: "retraction point",
                expected: point_dim,
                got: x.len() as u64,
            });
        }
        Self::validate_retraction_finite_with_checkpoint("retraction point", x, checkpoint)?;
        self.validate_point_domain_with_checkpoint(x, checkpoint)?;
        let param_dim = match self.param_dim() {
            Some(param_dim) => param_dim,
            None => {
                return Err(OptError::ManifoldInvalid {
                    what: try_owned_diagnostic(
                        "retract/manifold-diagnostic",
                        None,
                        None,
                        "retraction manifold has no representable parameter dimension",
                    )?,
                });
            }
        };
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
                let mut output = try_vec_capacity("retract/rn-output", None, None, x.len())?;
                for (component, (a, b)) in x.iter().zip(t).enumerate() {
                    checkpoint_retraction_work(component, checkpoint)?;
                    output.push(a + b);
                }
                self.validate_retraction_output_with_checkpoint(output, checkpoint)
            }
            Manifold::Sphere { .. } => {
                let mut y = try_vec_capacity("retract/sphere-candidate", None, None, x.len())?;
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
                let mut output = try_vec_capacity("retract/sphere-output", None, None, y.len())?;
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
                canonicalize_so3_representative(&mut out);
                let mut output = try_vec_capacity("retract/so3-output", None, None, out.len())?;
                output.extend_from_slice(&out);
                self.validate_retraction_output_with_checkpoint(output, checkpoint)
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                let mut cols = try_vec_capacity("retract/stiefel-columns", None, None, p)?;
                for column in 0..p {
                    checkpoint()?;
                    let mut candidate = try_vec_capacity("retract/stiefel-column", None, None, n)?;
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
                let mut output = try_vec_capacity("retract/stiefel-output", None, None, x.len())?;
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

fn plan_layout_bytes<T>(resource: &'static str) -> Result<u64, OptError> {
    u64::try_from(core::mem::size_of::<T>()).map_err(|_| OptError::DescentPlanOverflow { resource })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct IrEvalEnvelope {
    planning_work: u64,
    work_per_evaluation: u64,
    workspace_bytes: u64,
    #[cfg(test)]
    variables: u64,
    #[cfg(test)]
    binding_components: u64,
    #[cfg(test)]
    domain_work: u64,
    #[cfg(test)]
    prefix_nodes: u64,
    #[cfg(test)]
    prefix_edges: u64,
    #[cfg(test)]
    vector_components: u64,
    #[cfg(test)]
    reduction_components: u64,
}

#[derive(Clone, Copy)]
struct IrObjectivePlan<'problem> {
    problem: &'problem Problem,
    root: NodeId,
}

#[derive(Debug, Clone, Copy)]
struct IrPlanningLimit {
    base_work: u64,
    cap: u64,
}

fn charge_ir_planning_work<F>(
    work: &mut u64,
    limit: Option<IrPlanningLimit>,
    checkpoint: &mut F,
) -> Result<(), OptError>
where
    F: FnMut() -> Result<(), OptError>,
{
    let next = plan_add("IR evaluator planning work units", *work, 1)?;
    if let Some(limit) = limit {
        let required = plan_add("IR evaluator work units", limit.base_work, next)?;
        if required > limit.cap {
            return Err(OptError::DescentCapExceeded {
                resource: "work units",
                required,
                cap: limit.cap,
            });
        }
    }
    let stride =
        u64::try_from(EVAL_CHECKPOINT_STRIDE).map_err(|_| OptError::DescentPlanOverflow {
            resource: "IR evaluator planning work units",
        })?;
    if *work % stride == 0 {
        checkpoint()?;
    }
    *work = next;
    Ok(())
}

fn ir_domain_work(manifold: Manifold, point_dim: u32) -> Result<u64, OptError> {
    match manifold {
        Manifold::Rn { .. } => Ok(0),
        Manifold::Sphere { .. } | Manifold::So3 => Ok(u64::from(point_dim)),
        Manifold::Stiefel { n, p } => {
            let p = u64::from(p);
            let gram_entries = plan_mul(
                "IR evaluator work units",
                p,
                plan_add("IR evaluator work units", p, 1)?,
            )? / 2;
            plan_mul("IR evaluator work units", u64::from(n), gram_entries)
        }
    }
}

#[cfg(test)]
fn ir_eval_envelope_with_checkpoint<F>(
    problem: &Problem,
    root: NodeId,
    planning_limit: Option<IrPlanningLimit>,
    checkpoint: &mut F,
) -> Result<IrEvalEnvelope, OptError>
where
    F: FnMut() -> Result<(), OptError>,
{
    let caps = validate_eval_header(problem, root)?;
    validate_binding_frame_header(problem, &caps)?;
    ir_eval_envelope_with_caps(problem, root, &caps, planning_limit, checkpoint)
}

#[allow(clippy::too_many_lines)] // each term is an auditable part of one exact receipt
fn ir_eval_envelope_with_caps<F>(
    problem: &Problem,
    root: NodeId,
    caps: &AdmissionCaps,
    planning_limit: Option<IrPlanningLimit>,
    checkpoint: &mut F,
) -> Result<IrEvalEnvelope, OptError>
where
    F: FnMut() -> Result<(), OptError>,
{
    let root_index = usize::try_from(root.0).map_err(|_| OptError::DescentPlanOverflow {
        resource: "IR evaluator prefix nodes",
    })?;
    if root_index >= problem.exprs.len() {
        return Err(OptError::UnknownNode { id: root.0 });
    }
    let prefix_len = root_index
        .checked_add(1)
        .ok_or(OptError::DescentPlanOverflow {
            resource: "IR evaluator prefix nodes",
        })?;
    let prefix_nodes = u64::try_from(prefix_len).map_err(|_| OptError::DescentPlanOverflow {
        resource: "IR evaluator prefix nodes",
    })?;
    let variables =
        u64::try_from(problem.vars.len()).map_err(|_| OptError::DescentPlanOverflow {
            resource: "IR evaluator binding variables",
        })?;

    let mut planning_work = 0u64;
    let mut binding_components = 0u64;
    let mut domain_work = 0u64;
    for (variable_index, variable) in problem.vars.iter().enumerate() {
        let point_dim = match variable.manifold.point_dim() {
            Some(point_dim) => point_dim,
            None => {
                return Err(OptError::ManifoldInvalid {
                    what: try_owned_diagnostic(
                        "descent/ir-envelope-manifold-diagnostic",
                        None,
                        Some(VarId(u32::try_from(variable_index).map_err(|_| {
                            OptError::DescentPlanOverflow {
                                resource: "IR evaluator binding variables",
                            }
                        })?)),
                        "sealed IR variable has no representable runtime point dimension",
                    )?,
                });
            }
        };
        if point_dim > caps.max_point_dim {
            return Err(OptError::CapExceeded {
                what: "runtime binding point dimension",
                count: u64::from(point_dim),
                cap: u64::from(caps.max_point_dim),
            });
        }
        binding_components = plan_add(
            "IR evaluator work units",
            binding_components,
            u64::from(point_dim),
        )?;
        if binding_components > caps.max_total_point_storage {
            return Err(OptError::CapExceeded {
                what: "runtime binding point storage",
                count: binding_components,
                cap: caps.max_total_point_storage,
            });
        }
        domain_work = plan_add(
            "IR evaluator work units",
            domain_work,
            ir_domain_work(variable.manifold, point_dim)?,
        )?;
        let validation_work = plan_add("IR evaluator work units", binding_components, domain_work)?;
        if validation_work > caps.max_total_work {
            return Err(OptError::CapExceeded {
                what: "runtime binding validation work",
                count: validation_work,
                cap: caps.max_total_work,
            });
        }
        charge_ir_planning_work(&mut planning_work, planning_limit, checkpoint)?;
    }

    let mut prefix_edges = 0u64;
    let mut vector_components = 0u64;
    let mut reduction_components = 0u64;
    for (index, expr) in problem.exprs[..prefix_len].iter().enumerate() {
        charge_ir_planning_work(&mut planning_work, planning_limit, checkpoint)?;
        let node = NodeId(
            u32::try_from(index).map_err(|_| OptError::DescentPlanOverflow {
                resource: "IR evaluator prefix nodes",
            })?,
        );
        match problem.shape(node)? {
            Shape::Scalar => {}
            Shape::Vector(len) => {
                vector_components =
                    plan_add("IR evaluator work units", vector_components, u64::from(len))?;
            }
        }
        for _child in crate::ir::children(expr) {
            charge_ir_planning_work(&mut planning_work, planning_limit, checkpoint)?;
            prefix_edges = plan_add("IR evaluator work units", prefix_edges, 1)?;
        }
        let reduction_input = match expr {
            Expr::Dot(left, _) | Expr::NormSq(left) => Some(*left),
            _ => None,
        };
        if let Some(input) = reduction_input {
            match problem.shape(input)? {
                Shape::Vector(len) => {
                    reduction_components = plan_add(
                        "IR evaluator work units",
                        reduction_components,
                        u64::from(len),
                    )?;
                }
                actual => {
                    return Err(OptError::ShapeMismatch {
                        op: "IR evaluator reduction input",
                        left: actual,
                        right: Shape::Vector(1),
                    });
                }
            }
        }
    }
    checkpoint()?;

    // Conservative logical visits for one complete evaluator invocation:
    // three variable-level scans; binding/domain scalar visits; two table
    // initializations; reachability and child-edge traversal; the arena-prefix
    // sweep; node/result validation and memo publication; a second child read
    // pass during execution; vector construction/output scans; reductions; and
    // root/final publication. Prefix accounting intentionally overcharges
    // unreachable nodes before the objective root.
    let mut work_per_evaluation = 2u64;
    for term in [
        plan_mul("IR evaluator work units", variables, 3)?,
        binding_components,
        domain_work,
        plan_mul("IR evaluator work units", prefix_nodes, 6)?,
        plan_mul("IR evaluator work units", prefix_edges, 2)?,
        plan_mul("IR evaluator work units", vector_components, 2)?,
        reduction_components,
    ] {
        work_per_evaluation = plan_add("IR evaluator work units", work_per_evaluation, term)?;
    }

    // Logical requested Rust storage. Allocator metadata and any capacity in
    // excess of try_reserve_exact's request remain explicit no-claims.
    let mut workspace_bytes = plan_layout_bytes::<Value>("IR evaluator workspace bytes")?;
    for term in [
        plan_mul(
            "IR evaluator workspace bytes",
            prefix_nodes,
            plan_layout_bytes::<Option<Value>>("IR evaluator workspace bytes")?,
        )?,
        plan_mul(
            "IR evaluator workspace bytes",
            prefix_nodes,
            plan_layout_bytes::<bool>("IR evaluator workspace bytes")?,
        )?,
        plan_mul(
            "IR evaluator workspace bytes",
            prefix_nodes,
            plan_layout_bytes::<NodeId>("IR evaluator workspace bytes")?,
        )?,
        plan_mul(
            "IR evaluator workspace bytes",
            vector_components,
            plan_layout_bytes::<f64>("IR evaluator workspace bytes")?,
        )?,
        plan_layout_bytes::<Vec<Option<Value>>>("IR evaluator workspace bytes")?,
        plan_layout_bytes::<Vec<bool>>("IR evaluator workspace bytes")?,
        plan_layout_bytes::<Vec<NodeId>>("IR evaluator workspace bytes")?,
    ] {
        workspace_bytes = plan_add("IR evaluator workspace bytes", workspace_bytes, term)?;
    }

    Ok(IrEvalEnvelope {
        planning_work,
        work_per_evaluation,
        workspace_bytes,
        #[cfg(test)]
        variables,
        #[cfg(test)]
        binding_components,
        #[cfg(test)]
        domain_work,
        #[cfg(test)]
        prefix_nodes,
        #[cfg(test)]
        prefix_edges,
        #[cfg(test)]
        vector_components,
        #[cfg(test)]
        reduction_components,
    })
}

fn ir_objective_evaluation_count(param_dim: u32, active_steps: u64) -> Result<u64, OptError> {
    if active_steps == 0 {
        return Ok(1);
    }
    let probes_per_step = plan_mul("IR objective evaluation count", u64::from(param_dim), 2)?;
    plan_add(
        "IR objective evaluation count",
        plan_mul(
            "IR objective evaluation count",
            probes_per_step,
            active_steps,
        )?,
        2,
    )
}

fn compose_ir_descent_envelope(
    descent: DescentEnvelope,
    ir: IrEvalEnvelope,
    param_dim: u32,
    active_steps: u64,
) -> Result<DescentEnvelope, OptError> {
    let evaluations = ir_objective_evaluation_count(param_dim, active_steps)?;
    let evaluator_work = plan_mul(
        "IR evaluator work units",
        evaluations,
        ir.work_per_evaluation,
    )?;
    let owned_work = plan_add(
        "IR evaluator work units",
        plan_add(
            "IR evaluator work units",
            descent.owned_work,
            ir.planning_work,
        )?,
        evaluator_work,
    )?;
    let scratch_bytes = plan_add(
        "IR evaluator workspace bytes",
        descent.scratch_bytes,
        ir.workspace_bytes,
    )?;
    Ok(DescentEnvelope {
        owned_work,
        scratch_bytes,
    })
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
    /// Hard upper bound for fs-opt-owned scalar/logical work admitted up front
    /// (default `2^24`). [`descend_ir`] includes its one-time prefix planner and
    /// the maximum budget-reachable evaluator invocations; an opaque
    /// [`descend_fn`] closure remains caller-owned and excluded.
    pub max_work_units: NonZeroU64,
    /// Hard upper bound for peak fs-opt-owned logical requested workspace, in
    /// bytes (default 1 GiB). [`descend_ir`] adds its evaluator tables/vector
    /// payload; an opaque [`descend_fn`] closure remains caller-owned and
    /// excluded.
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
    /// The explicit evaluation limit could not fund the next atomic step. The
    /// report retains the last completely landed iterate as a restart point;
    /// a fresh invocation has a new segment-local evaluation/site ledger.
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
    /// Objective evaluations spent by this invocation/segment.
    pub evals: u64,
    /// Steps taken.
    pub steps_taken: u32,
    /// Typed successful terminal state.
    pub stop: DescentStop,
    /// True when the P4 budget stopped the run (a receipt, not an
    /// error — the point is still valid). Exactly equivalent to
    /// `stop == DescentStop::EvaluationLimit`.
    pub budget_stopped: bool,
    /// Conservative fs-opt-owned work bound admitted before the first
    /// objective. IR-driven reports include evaluator planning/execution;
    /// opaque raw-closure work is excluded.
    pub work_upper_bound: u64,
    /// Conservative peak fs-opt-owned logical requested workspace admitted
    /// before the first objective. IR-driven reports include evaluator
    /// workspace; opaque raw-closure workspace is excluded.
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
/// The raw-closure resource envelope excludes caller objective work/allocation;
/// [`descend_ir`] composes its owned evaluator plan, maximum invocation work,
/// and peak logical workspace into the same caps and report fields.
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
/// [`OptError::DescentCapExceeded`] / [`OptError::DescentPlanOverflow`] /
/// [`OptError::RuntimeAllocationRefused`].
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
    descend_fn_checked(manifold, &checked_f, x0, opts, eval_limit, None, cx)
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
    probe: &mut [f64],
    fd_h: f64,
    cx: &Cx<'_>,
) -> Result<(), OptError> {
    for parameter in 0..probe.len() {
        descent_checkpoint(cx)?;
        probe[parameter] = fd_h;
        drop(manifold.retract_with_cx(x0, probe, cx)?);
        descent_checkpoint(cx)?;
        probe[parameter] = -fd_h;
        drop(manifold.retract_with_cx(x0, probe, cx)?);
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
    ir_plan: Option<IrObjectivePlan<'_>>,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    manifold.validate(&AdmissionCaps::default())?;
    let point_dim = match manifold.point_dim() {
        Some(point_dim) => point_dim,
        None => {
            return Err(OptError::ManifoldInvalid {
                what: try_owned_diagnostic(
                    "descent/manifold-diagnostic",
                    None,
                    None,
                    "descent manifold has no representable point dimension",
                )?,
            });
        }
    };
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
    let param_dim = match manifold.param_dim() {
        Some(param_dim) => param_dim,
        None => {
            return Err(OptError::ManifoldInvalid {
                what: try_owned_diagnostic(
                    "descent/manifold-diagnostic",
                    None,
                    None,
                    "descent manifold has no representable parameter dimension",
                )?,
            });
        }
    };
    let pd = match usize::try_from(param_dim) {
        Ok(pd) => pd,
        Err(_) => {
            return Err(OptError::ManifoldInvalid {
                what: try_owned_diagnostic(
                    "descent/manifold-diagnostic",
                    None,
                    None,
                    "descent parameter dimension does not fit target usize",
                )?,
            });
        }
    };
    let atomic_step_evals = u64::from(param_dim).saturating_mul(2).saturating_add(1);
    let reachable_steps = maximum_landed_steps(opts.steps, param_dim, eval_limit);
    let base_envelope = descent_envelope(manifold, point_dim, param_dim, opts.steps, eval_limit)?;
    if base_envelope.owned_work > opts.max_work_units.get() {
        return Err(OptError::DescentCapExceeded {
            resource: "work units",
            required: base_envelope.owned_work,
            cap: opts.max_work_units.get(),
        });
    }
    if base_envelope.scratch_bytes > opts.max_workspace_bytes.get() {
        return Err(OptError::DescentCapExceeded {
            resource: "workspace bytes",
            required: base_envelope.scratch_bytes,
            cap: opts.max_workspace_bytes.get(),
        });
    }
    let envelope = if let Some(plan) = ir_plan {
        let caps = validate_eval_header(plan.problem, plan.root)?;
        validate_binding_frame_header(plan.problem, &caps)?;
        if plan.problem.vars.len() != 1 {
            return Err(OptError::BindingCount {
                vars: u32::try_from(plan.problem.vars.len()).map_err(|_| {
                    OptError::DescentPlanOverflow {
                        resource: "IR evaluator binding variables",
                    }
                })?,
                got: 1,
            });
        }
        let planning_limit = IrPlanningLimit {
            base_work: base_envelope.owned_work,
            cap: opts.max_work_units.get(),
        };
        let ir = ir_eval_envelope_with_caps(
            plan.problem,
            plan.root,
            &caps,
            Some(planning_limit),
            &mut || descent_checkpoint(cx),
        )?;
        compose_ir_descent_envelope(base_envelope, ir, param_dim, reachable_steps)?
    } else {
        base_envelope
    };
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
    let mut tangent = if reachable_steps > 0 {
        try_filled_runtime_vector("descent/tangent", pd, 0.0, &mut || descent_checkpoint(cx))?
    } else {
        Vec::new()
    };
    let mut gradient = if reachable_steps > 0 {
        try_filled_runtime_vector("descent/gradient", pd, 0.0, &mut || descent_checkpoint(cx))?
    } else {
        Vec::new()
    };
    let mut step = if reachable_steps > 0 {
        try_vec_capacity("descent/step", None, None, pd)?
    } else {
        Vec::new()
    };
    if reachable_steps > 0 {
        validate_initial_probe_retractions(manifold, x0, &mut tangent, opts.fd_h, cx)?;
    }
    descent_checkpoint(cx)?;
    let mut x = try_clone_runtime_slice("descent/point", x0, &mut || descent_checkpoint(cx))?;
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
        for (i, gi) in gradient.iter_mut().enumerate() {
            descent_checkpoint(cx)?;
            tangent[i] = opts.fd_h;
            descent_checkpoint(cx)?;
            let xp = manifold.retract_with_cx(&x, &tangent, cx)?;
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
            tangent[i] = -opts.fd_h;
            let xm = manifold.retract_with_cx(&x, &tangent, cx)?;
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
            tangent[i] = 0.0;
        }
        descent_checkpoint(cx)?;
        step.clear();
        for (component, gradient) in gradient.iter().enumerate() {
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

/// Toy descent of a problem's FIRST objective over its sole variable
/// (the IR-driven variant; enforces `problem.budget` per P4). Multi-variable
/// problems refuse with [`OptError::BindingCount`] because this API accepts one
/// point and cannot fabricate a complete runtime frame. The caller's
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
    descend_fn_checked(
        manifold,
        &f,
        x0,
        opts,
        problem.budget.limit,
        Some(IrObjectivePlan {
            problem,
            root: obj.node,
        }),
        cx,
    )
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

    fn layout_bytes<T>() -> u64 {
        u64::try_from(core::mem::size_of::<T>()).expect("Rust layout fits the receipt domain")
    }

    fn ir_prefix_problem(with_suffix: bool) -> (Problem, NodeId) {
        let mut builder = ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: 3 }, Dims::NONE)
            .expect("variable");
        let point = builder.var_ref(variable).expect("point");
        let negated = builder.neg(point).expect("negated point");
        let cancelled = builder.add(point, negated).expect("vector cancellation");
        let root = builder.norm_sq(cancelled).expect("scalar objective");
        assert_eq!(builder.var_ref(variable), Ok(point), "var_ref hash-consing");
        assert_eq!(builder.neg(point), Ok(negated), "negation hash-consing");
        assert_eq!(
            builder.add(point, negated),
            Ok(cancelled),
            "addition hash-consing"
        );
        assert_eq!(
            builder.norm_sq(cancelled),
            Ok(root),
            "reduction hash-consing"
        );
        builder
            .objective(root, crate::ir::Sense::Minimize, 1.0)
            .expect("objective");
        if with_suffix {
            let scale = builder.konst(7.0, Dims::NONE).expect("suffix scalar");
            builder
                .mul(scale, point)
                .expect("unrelated vector-valued suffix");
        }
        (builder.finish(), root)
    }

    fn ir_norm_problem() -> (Problem, NodeId) {
        let mut builder = ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: 3 }, Dims::NONE)
            .expect("variable");
        let point = builder.var_ref(variable).expect("point");
        let root = builder.norm_sq(point).expect("squared norm");
        builder
            .objective(root, crate::ir::Sense::Minimize, 1.0)
            .expect("objective");
        (builder.finish(), root)
    }

    #[test]
    fn ir_evaluator_envelope_is_exact_prefix_scoped_and_checkpointed() {
        let (problem, root) = ir_prefix_problem(false);
        assert_eq!(root, NodeId(3));
        let mut polls = 0usize;
        let envelope = ir_eval_envelope_with_checkpoint(&problem, root, None, &mut || {
            polls += 1;
            Ok(())
        })
        .expect("bounded evaluator envelope");
        let workspace_bytes = layout_bytes::<Value>()
            + 4 * (layout_bytes::<Option<Value>>()
                + layout_bytes::<bool>()
                + layout_bytes::<NodeId>())
            + 9 * layout_bytes::<f64>()
            + layout_bytes::<Vec<Option<Value>>>()
            + layout_bytes::<Vec<bool>>()
            + layout_bytes::<Vec<NodeId>>();
        assert_eq!(
            envelope,
            IrEvalEnvelope {
                planning_work: 9,
                work_per_evaluation: 61,
                workspace_bytes,
                variables: 1,
                binding_components: 3,
                domain_work: 0,
                prefix_nodes: 4,
                prefix_edges: 4,
                vector_components: 9,
                reduction_components: 3,
            }
        );
        assert_eq!(polls, 2, "work item zero plus terminal publication");

        let (with_suffix, same_root) = ir_prefix_problem(true);
        let suffix_envelope =
            ir_eval_envelope_with_checkpoint(&with_suffix, same_root, None, &mut || Ok(()))
                .expect("suffix-insensitive envelope");
        assert_eq!(
            suffix_envelope, envelope,
            "nodes after the objective root cannot inflate its evaluator receipt"
        );

        let mut capped_polls = 0usize;
        let capped = ir_eval_envelope_with_checkpoint(
            &problem,
            root,
            Some(IrPlanningLimit {
                base_work: 0,
                cap: 8,
            }),
            &mut || {
                capped_polls += 1;
                Ok(())
            },
        );
        assert_eq!(
            capped,
            Err(OptError::DescentCapExceeded {
                resource: "work units",
                required: 9,
                cap: 8,
            })
        );
        assert_eq!(
            capped_polls, 1,
            "the ninth item refuses before work or a terminal checkpoint"
        );
    }

    #[test]
    fn ir_evaluator_planner_replays_and_every_poll_aborts_atomically() {
        let (problem, root, _binding) = checkpoint_problem();
        let mut polls = 0usize;
        let envelope = ir_eval_envelope_with_checkpoint(&problem, root, None, &mut || {
            polls += 1;
            Ok(())
        })
        .expect("large evaluator plan");
        let stride = u64::try_from(EVAL_CHECKPOINT_STRIDE).expect("checkpoint stride fits u64");
        let expected_polls = usize::try_from(envelope.planning_work.div_ceil(stride))
            .expect("admitted planner poll count fits usize")
            + 1;
        assert_eq!(polls, expected_polls, "work-stride polls plus publication");

        let mut replay_polls = 0usize;
        let replay = ir_eval_envelope_with_checkpoint(&problem, root, None, &mut || {
            replay_polls += 1;
            Ok(())
        })
        .expect("replayed evaluator plan");
        assert_eq!(replay, envelope);
        assert_eq!(replay_polls, polls, "planner cadence must replay exactly");

        for target in 1..=polls {
            let mut observed = 0usize;
            let result = ir_eval_envelope_with_checkpoint(&problem, root, None, &mut || {
                observed += 1;
                if observed == target {
                    Err(OptError::Cancelled)
                } else {
                    Ok(())
                }
            });
            assert_eq!(result, Err(OptError::Cancelled), "target poll {target}");
            assert_eq!(
                observed, target,
                "planner continued after cancellation at poll {target}"
            );
        }
    }

    #[test]
    fn ir_evaluator_domain_and_composition_receipts_are_exact() {
        assert_eq!(ir_domain_work(Manifold::Rn { dim: 3 }, 3), Ok(0));
        assert_eq!(ir_domain_work(Manifold::Sphere { ambient: 3 }, 3), Ok(3));
        assert_eq!(ir_domain_work(Manifold::So3, 4), Ok(4));
        assert_eq!(ir_domain_work(Manifold::Stiefel { n: 3, p: 2 }, 6), Ok(9));
        assert_eq!(
            ir_domain_work(
                Manifold::Stiefel {
                    n: u32::MAX,
                    p: u32::MAX,
                },
                u32::MAX,
            ),
            Err(OptError::DescentPlanOverflow {
                resource: "IR evaluator work units",
            })
        );

        let (problem, root) = ir_norm_problem();
        let ir = ir_eval_envelope_with_checkpoint(&problem, root, None, &mut || Ok(()))
            .expect("small evaluator envelope");
        assert_eq!(ir.planning_work, 4);
        assert_eq!(ir.work_per_evaluation, 31);
        assert_eq!(ir_objective_evaluation_count(3, 0), Ok(1));
        assert_eq!(ir_objective_evaluation_count(3, 1), Ok(8));

        let zero_base = descent_envelope(Manifold::Rn { dim: 3 }, 3, 3, 0, EvalLimit::Unlimited)
            .expect("zero-step descent envelope");
        let zero =
            compose_ir_descent_envelope(zero_base, ir, 3, 0).expect("zero-step composed envelope");
        assert_eq!(zero_base.owned_work, 78);
        assert_eq!(zero.owned_work, 113);
        assert_eq!(
            zero.scratch_bytes,
            zero_base.scratch_bytes + ir.workspace_bytes
        );

        let one_base = descent_envelope(Manifold::Rn { dim: 3 }, 3, 3, 1, EvalLimit::Unlimited)
            .expect("one-step descent envelope");
        let one =
            compose_ir_descent_envelope(one_base, ir, 3, 1).expect("one-step composed envelope");
        assert_eq!(one_base.owned_work, 1_083);
        assert_eq!(one.owned_work, 1_335);
        assert_eq!(
            one.scratch_bytes,
            one_base.scratch_bytes + ir.workspace_bytes
        );

        let limit_7 = EvalLimit::Limited(NonZeroU64::new(7).expect("positive limit"));
        let limit_8 = EvalLimit::Limited(NonZeroU64::new(8).expect("positive limit"));
        assert_eq!(maximum_landed_steps(1, 3, limit_7), 0);
        assert_eq!(maximum_landed_steps(1, 3, limit_8), 1);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // exact cap/call receipts stay in one end-to-end tripwire
    fn ir_descent_combined_caps_and_evaluator_call_receipts_are_exact() {
        let (problem, _root) = ir_norm_problem();
        let start = [1.0, -2.0, 0.5];
        let zero_options = DescentOptions {
            steps: 0,
            ..DescentOptions::default()
        };
        let zero_gate = fs_exec::CancelGate::new_clock_free();
        let zero_finalizations = std::cell::Cell::new(0usize);
        let zero_report = with_test_cx(&zero_gate, |cx| {
            descend_ir_with_eval_checkpoint(&problem, &start, zero_options, cx, &|phase| {
                if phase == EvalPhase::Finalize {
                    zero_finalizations.set(zero_finalizations.get() + 1);
                }
                Ok(())
            })
            .expect("zero-step IR descent")
        });
        assert_eq!(zero_report.work_upper_bound, 113);
        assert_eq!(zero_report.evals, 1);
        assert_eq!(zero_report.steps_taken, 0);
        assert_eq!(zero_finalizations.get(), 1);

        let exact = DescentOptions {
            max_work_units: NonZeroU64::new(zero_report.work_upper_bound)
                .expect("positive work receipt"),
            max_workspace_bytes: NonZeroU64::new(zero_report.workspace_upper_bound_bytes)
                .expect("positive workspace receipt"),
            ..zero_options
        };
        let exact_gate = fs_exec::CancelGate::new_clock_free();
        with_test_cx(&exact_gate, |cx| {
            descend_ir(&problem, &start, exact, cx).expect("exact combined caps admit");
        });

        for options in [
            DescentOptions {
                max_work_units: NonZeroU64::new(zero_report.work_upper_bound - 1)
                    .expect("one-short work cap"),
                ..zero_options
            },
            DescentOptions {
                max_workspace_bytes: NonZeroU64::new(zero_report.workspace_upper_bound_bytes - 1)
                    .expect("one-short workspace cap"),
                ..zero_options
            },
        ] {
            let gate = fs_exec::CancelGate::new_clock_free();
            let evaluator_polls = std::cell::Cell::new(0usize);
            let error = with_test_cx(&gate, |cx| {
                descend_ir_with_eval_checkpoint(&problem, &start, options, cx, &|_phase| {
                    evaluator_polls.set(evaluator_polls.get() + 1);
                    Ok(())
                })
                .expect_err("one-short combined cap must refuse")
            });
            match error {
                OptError::DescentCapExceeded {
                    resource: "work units",
                    required,
                    cap,
                } => {
                    assert_eq!(required, zero_report.work_upper_bound);
                    assert_eq!(cap + 1, required);
                }
                OptError::DescentCapExceeded {
                    resource: "workspace bytes",
                    required,
                    cap,
                } => {
                    assert_eq!(required, zero_report.workspace_upper_bound_bytes);
                    assert_eq!(cap + 1, required);
                }
                other => panic!("unexpected combined-cap refusal: {other:?}"),
            }
            assert_eq!(evaluator_polls.get(), 0, "refusal must dominate f0");
        }

        let one_options = DescentOptions {
            steps: 1,
            ..DescentOptions::default()
        };
        let one_gate = fs_exec::CancelGate::new_clock_free();
        let one_finalizations = std::cell::Cell::new(0usize);
        let one_report = with_test_cx(&one_gate, |cx| {
            descend_ir_with_eval_checkpoint(&problem, &start, one_options, cx, &|phase| {
                if phase == EvalPhase::Finalize {
                    one_finalizations.set(one_finalizations.get() + 1);
                }
                Ok(())
            })
            .expect("one-step IR descent")
        });
        assert_eq!(one_report.work_upper_bound, 1_335);
        assert_eq!(one_report.evals, 8);
        assert_eq!(one_report.steps_taken, 1);
        assert_eq!(one_finalizations.get(), 8);
    }

    #[test]
    fn ir_planning_overflow_and_multi_binding_refusal_are_typed() {
        let mut work = u64::MAX;
        let mut checkpoints = 0usize;
        assert_eq!(
            charge_ir_planning_work(&mut work, None, &mut || {
                checkpoints += 1;
                Ok(())
            }),
            Err(OptError::DescentPlanOverflow {
                resource: "IR evaluator planning work units",
            })
        );
        assert_eq!(work, u64::MAX);
        assert_eq!(checkpoints, 0, "overflow refuses before a checkpoint");
        assert_eq!(
            ir_objective_evaluation_count(1, u64::MAX),
            Err(OptError::DescentPlanOverflow {
                resource: "IR objective evaluation count",
            })
        );

        let mut builder = ProblemBuilder::new();
        let first = builder
            .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("first variable");
        builder
            .var("y", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("second variable");
        let first_point = builder.var_ref(first).expect("first point");
        let root = builder.norm_sq(first_point).expect("objective root");
        builder
            .objective(root, crate::ir::Sense::Minimize, 1.0)
            .expect("objective");
        let problem = builder.finish();
        let gate = fs_exec::CancelGate::new_clock_free();
        with_test_cx(&gate, |cx| {
            let invalid_options = DescentOptions {
                lr: 0.0,
                ..DescentOptions::default()
            };
            assert!(matches!(
                descend_ir(&problem, &[1.0], invalid_options, cx),
                Err(OptError::BadParam { .. })
            ));
            assert!(matches!(
                descend_ir(&problem, &[1.0], DescentOptions::default(), cx),
                Err(OptError::BindingCount { vars: 2, got: 1 })
            ));
        });
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
    fn runtime_scratch_builders_are_fallible_checkpointed_and_bit_exact() {
        let mut source = vec![1.0; RETRACTION_CHECKPOINT_STRIDE * 2 + 1];
        source[0] = -0.0;
        source[RETRACTION_CHECKPOINT_STRIDE] = f64::MIN_POSITIVE;
        source[RETRACTION_CHECKPOINT_STRIDE * 2] = -f64::MIN_POSITIVE;

        let mut clone_polls = 0usize;
        let cloned = try_clone_runtime_slice("test/runtime-clone", &source, &mut || {
            clone_polls += 1;
            Ok(())
        })
        .expect("bounded runtime clone");
        assert_eq!(clone_polls, 4, "indices 0/256/512 plus terminal");
        assert_eq!(
            cloned
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            source
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );

        let mut fill_polls = 0usize;
        let filled =
            try_filled_runtime_vector("test/runtime-fill", source.len(), -0.0, &mut || {
                fill_polls += 1;
                Ok(())
            })
            .expect("bounded runtime fill");
        assert_eq!(fill_polls, 4, "indices 0/256/512 plus terminal");
        assert!(
            filled
                .iter()
                .all(|value| value.to_bits() == (-0.0f64).to_bits())
        );

        for target in 1..=clone_polls {
            let mut calls = 0usize;
            let result = try_clone_runtime_slice("test/runtime-clone", &source, &mut || {
                calls += 1;
                if calls == target {
                    Err(OptError::Cancelled)
                } else {
                    Ok(())
                }
            });
            assert_eq!(result, Err(OptError::Cancelled), "target poll {target}");
            assert_eq!(calls, target, "copy continued after target poll {target}");
        }

        for target in 1..=fill_polls {
            let mut calls = 0usize;
            let result =
                try_filled_runtime_vector("test/runtime-fill", source.len(), -0.0, &mut || {
                    calls += 1;
                    if calls == target {
                        Err(OptError::Cancelled)
                    } else {
                        Ok(())
                    }
                });
            assert_eq!(result, Err(OptError::Cancelled), "target poll {target}");
            assert_eq!(calls, target, "fill continued after target poll {target}");
        }

        let mut checkpoint = || Ok(());
        assert_eq!(
            try_filled_runtime_vector(
                "test/runtime-fill-overflow",
                usize::MAX,
                0.0,
                &mut checkpoint,
            ),
            Err(OptError::RuntimeAllocationRefused {
                path: "test/runtime-fill-overflow",
                node: None,
                variable: None,
                elements: allocation_len(usize::MAX),
                element_bytes: core::mem::size_of::<f64>() as u64,
            })
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
