//! Native reverse-mode execution of the live, sealed optimization IR.
//!
//! Compile once, evaluate every reachable node once, then obtain a weighted
//! vector-Jacobian product in one reverse sweep. Vector lanes stay contiguous;
//! shared subexpressions accumulate all parent contributions. Elementary
//! primals use the same `fs_math::det` functions and reduction order as `eval`.
//! Gradients are chain-rule derivatives of the mathematical operators, not
//! derivatives of their floating-point implementations or error certificates.
//!
//! Smooth algebra is executable by default. Explicit [`physics`] bindings may
//! supply solved scalar PDE residuals and their total first-order derivatives.
//! Unbound physics, UQ and kinks still refuse; no PDE adjoint is fabricated.
//! Existing `eval` and finite-difference descent are unchanged.

use crate::{BindingFrame, Expr, NodeId, OptError, Problem, Shape, VarId, children};
use fs_exec::Cx;

pub mod physics;
use physics::{CompiledPhysics, PhysicsBinding, PhysicsError};

/// Explicit storage/work envelope for a compiled reverse program.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReverseLimits {
    /// Maximum graph-node, root, or variable count (each checked separately).
    pub max_nodes: usize,
    /// Maximum scalar lanes, including every declared variable's point.
    /// Each evaluation stores this many primals; a pullback also needs this
    /// many adjoints. Output gradients add at most one variable-point copy.
    pub max_scalar_slots: usize,
}

/// A checked refusal from compilation, primal execution, or a pullback.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseError {
    /// Existing graph, binding, domain, resource, or cancellation refusal.
    Evaluation(OptError),
    /// A supplied physics binding disagrees with its node or declaration.
    PhysicsBinding { node: NodeId, what: &'static str },
    /// Preserve a model's original typed refusal and the responsible node.
    Physics { node: NodeId, source: PhysicsError },
    /// A provider omitted or invented derivative components.
    PhysicsGradientLength { node: NodeId, expected: usize, actual: usize },
    /// No smooth derivative is claimed for this reachable operation.
    Nonsmooth {
        /// Offending live-IR node.
        node: NodeId,
    },
    /// A pullback needs exactly one seed per compiled scalar root.
    SeedCount {
        /// Compiled scalar roots.
        expected: usize,
        /// Supplied seeds.
        actual: usize,
    },
    /// A derivative contribution or its accumulated sum became non-finite.
    NonFiniteAdjoint {
        /// Operation producing the contribution (or seeded root).
        node: NodeId,
        /// Component within the destination slot.
        component: usize,
        /// Exact offending IEEE-754 bits.
        bits: u64,
    },
}

impl From<OptError> for ReverseError {
    fn from(error: OptError) -> Self { Self::Evaluation(error) }
}

impl core::fmt::Display for ReverseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Evaluation(error) => write!(f, "{error}"),
            Self::PhysicsBinding { node, what } => write!(f, "physics node {}: {what}", node.0),
            Self::Physics { node, source } => write!(f, "physics node {}: {source}", node.0),
            Self::PhysicsGradientLength { node, expected, actual } => write!(f, "physics node {} needs {expected} gradient components, received {actual}", node.0),
            Self::Nonsmooth { node } => write!(f, "node {} has no smooth reverse rule", node.0),
            Self::SeedCount { expected, actual } => write!(f, "reverse program needs {expected} seeds, received {actual}"),
            Self::NonFiniteAdjoint { node, component, bits } => write!(f, "node {} produced a non-finite adjoint at component {component}: {bits:#018x}", node.0),
        }
    }
}

impl std::error::Error for ReverseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Evaluation(error) => Some(error),
            Self::Physics { source, .. } => Some(source),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct Slot { start: usize, len: usize }

/// A compiled smooth subgraph bound by reference to its immutable problem.
/// Scalar roots retain caller order, including duplicates.
#[derive(Debug)]
pub struct ReverseProgram<'problem> {
    problem: &'problem Problem,
    roots: Vec<NodeId>,
    slots: Vec<Slot>,
    variables: Vec<Slot>,
    order: Vec<NodeId>,
    scalar_slots: usize,
    physics: Vec<CompiledPhysics<'problem>>,
}

fn cap(what: &'static str, count: usize, limit: usize) -> Result<(), OptError> {
    if count > limit {
        Err(OptError::CapExceeded { what, count: count as u64, cap: limit as u64 })
    } else { Ok(()) }
}

fn capacity<T>(len: usize) -> Result<Vec<T>, OptError> {
    let mut out = Vec::new();
    out.try_reserve_exact(len).map_err(|_| OptError::RuntimeAllocationRefused {
        path: "reverse/workspace", node: None, variable: None,
        elements: len as u64, element_bytes: core::mem::size_of::<T>() as u64,
    })?;
    Ok(out)
}

fn poll(cx: Option<&Cx>) -> Result<(), OptError> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| OptError::Cancelled)?; }
    Ok(())
}

fn tick(index: usize, cx: Option<&Cx>) -> Result<(), OptError> {
    if index % 256 == 0 { poll(cx)?; }
    Ok(())
}

fn zeros(len: usize, cx: Option<&Cx>) -> Result<Vec<f64>, OptError> {
    let mut out = capacity(len)?;
    for i in 0..len { tick(i, cx)?; out.push(0.0); }
    Ok(out)
}

impl<'problem> ReverseProgram<'problem> {
    /// Compile the union of the requested scalar subgraphs without recursion.
    /// Bounds are checked before allocating their corresponding storage.
    /// Compilation performs no objective evaluations and consumes no problem
    /// evaluation budget; the calling optimizer owns that accounting.
    pub fn new(problem: &'problem Problem, roots: &[NodeId], limits: ReverseLimits) -> Result<Self, ReverseError> {
        Self::new_with_physics(problem, roots, limits, &[])
    }

    /// Bind explicitly declared scalar physics studies to their live executors.
    /// Every supplied signature is checked, including unreachable bindings.
    /// Retained physics gradients are charged to `max_scalar_slots`. Provider
    /// scratch space, solves and their work budgets remain provider-owned.
    pub fn new_with_physics(
        problem: &'problem Problem,
        roots: &[NodeId],
        limits: ReverseLimits,
        bindings: &[PhysicsBinding<'problem>],
    ) -> Result<Self, ReverseError> {
        cap("reverse graph nodes", problem.exprs().len(), limits.max_nodes)?;
        cap("reverse roots", roots.len(), limits.max_nodes)?;
        cap("reverse variables", problem.vars().len(), limits.max_nodes)?;
        cap("reverse physics bindings", bindings.len(), limits.max_nodes)?;
        let bindings = physics::admitted_bindings(problem, bindings)?;
        let mut needed = capacity(problem.exprs().len())?;
        needed.resize(problem.exprs().len(), false);
        let mut retained_roots = capacity(roots.len())?;
        for &root in roots {
            if problem.shape(root)? != Shape::Scalar { return Err(OptError::NotScalar { node: root.0 }.into()); }
            needed[root.0 as usize] = true;
            retained_roots.push(root);
        }
        // Every child precedes its parent in a sealed Problem.
        for i in (0..needed.len()).rev() {
            if needed[i] {
                match &problem.exprs()[i] {
                    Expr::Min(..) | Expr::Max(..) | Expr::Abs(..) => return Err(ReverseError::Nonsmooth { node: NodeId(i as u32) }),
                    Expr::PdeResidual { .. } if bindings.binary_search_by_key(&NodeId(i as u32), |binding| binding.node).is_ok() => {},
                    Expr::PdeResidual { .. } | Expr::Expectation { .. } | Expr::Cvar { .. } | Expr::Quantile { .. } => return Err(OptError::Unevaluable {
                        node: i as u32, kind: "reverse program requires a live algebraic derivative; physics/UQ needs its own executor",
                    }.into()),
                    Expr::Var(_) | Expr::Component { .. } | Expr::Const { .. }
                    | Expr::Add(..) | Expr::Sub(..) | Expr::Mul(..) | Expr::Div(..)
                    | Expr::Neg(_) | Expr::Powi { .. } | Expr::Sqrt(_) | Expr::Exp(_)
                    | Expr::Ln(_) | Expr::Tanh(_) | Expr::Dot(..) | Expr::NormSq(_) => {}
                }
                for child in children(&problem.exprs()[i]) { needed[child.0 as usize] = true; }
            }
        }
        let mut scalar_slots = 0usize;
        let mut allocate = |len: usize| -> Result<Slot, OptError> {
            let end = scalar_slots.checked_add(len).ok_or(OptError::CapExceeded {
                what: "reverse scalar slots", count: u64::MAX, cap: limits.max_scalar_slots as u64,
            })?;
            cap("reverse scalar slots", end, limits.max_scalar_slots)?;
            let slot = Slot { start: scalar_slots, len };
            scalar_slots = end;
            Ok(slot)
        };
        let mut variables = capacity(problem.vars().len())?;
        for variable in problem.vars() {
            variables.push(allocate(variable.manifold.point_dim().expect("sealed manifold") as usize)?);
        }
        let mut slots = capacity(needed.len())?;
        slots.resize(needed.len(), Slot::default());
        let mut order = capacity(needed.iter().filter(|&&value| value).count())?;
        for (i, &used) in needed.iter().enumerate() {
            if !used { continue; }
            let node = NodeId(i as u32);
            slots[i] = match problem.expr(node)? {
                Expr::Var(var) => variables[var.0 as usize],
                _ => allocate(match problem.shape(node)? { Shape::Scalar => 1, Shape::Vector(n) => n as usize })?,
            };
            order.push(node);
        }
        let mut physics = capacity(bindings.len())?;
        for binding in bindings {
            if !needed[binding.node.0 as usize] { continue; }
            let Expr::PdeResidual { over, .. } = problem.expr(binding.node)? else {
                unreachable!("binding admission requires a PDE residual");
            };
            physics.push(CompiledPhysics {
                node: binding.node, model: binding.model, variable: *over,
                gradient: allocate(variables[over.0 as usize].len)?,
            });
        }
        Ok(Self { problem, roots: retained_roots, slots, variables, order, scalar_slots, physics })
    }

    /// Number of reachable IR nodes; independent of the number of variables.
    #[must_use]
    pub fn node_count(&self) -> usize { self.order.len() }

    /// Scalar lanes retained by each primal or adjoint workspace.
    #[must_use]
    pub fn scalar_slots(&self) -> usize { self.scalar_slots }

    /// Execute one primal sweep. `None` explicitly opts out of cancellation.
    /// Binding validation delegates to the existing BindingFrame authority;
    /// it is bracketed by polls. Numerical sweeps and workspace initialization
    /// additionally poll every 256 scalar lanes and at every node boundary.
    /// A failed/cancelled call publishes no partial evaluation.
    pub fn evaluate(&self, bindings: &[Vec<f64>], cx: Option<&Cx>) -> Result<ReverseEvaluation<'_, 'problem>, ReverseError> {
        poll(cx)?;
        if bindings.len() != self.variables.len() {
            return Err(OptError::BindingCount { vars: self.variables.len() as u32, got: bindings.len() as u64 }.into());
        }
        let _frame = BindingFrame::new(self.problem, bindings.iter().enumerate().map(|(i, values)| (VarId(i as u32), values.as_slice())))?;
        poll(cx)?;
        let mut values = zeros(self.scalar_slots, cx)?;
        for (slot, binding) in self.variables.iter().zip(bindings) {
            for (j, &value) in binding.iter().enumerate() { tick(j, cx)?; values[slot.start + j] = value; }
        }
        for &node in &self.order {
            poll(cx)?;
            let out = self.slots[node.0 as usize];
            if let Ok(index) = self.physics.binary_search_by_key(&node, |binding| binding.node) {
                let binding = &self.physics[index];
                let input = self.variables[binding.variable.0 as usize];
                let sample = binding.model.value_gradient(
                    &values[input.start..input.start + input.len], cx,
                ).map_err(|source| match source {
                    PhysicsError::Cancelled => ReverseError::Evaluation(OptError::Cancelled),
                    source => ReverseError::Physics { node, source },
                })?;
                poll(cx)?;
                if sample.gradient.len() != input.len {
                    return Err(ReverseError::PhysicsGradientLength {
                        node, expected: input.len, actual: sample.gradient.len(),
                    });
                }
                if !sample.value.is_finite() {
                    return Err(OptError::EvalNonFinite {
                        node: node.0, component: None, bits: sample.value.to_bits(),
                    }.into());
                }
                values[out.start] = sample.value;
                for (i, &derivative) in sample.gradient.iter().enumerate() {
                    tick(i, cx)?;
                    if !derivative.is_finite() {
                        return Err(ReverseError::NonFiniteAdjoint {
                            node, component: i, bits: derivative.to_bits(),
                        });
                    }
                    values[binding.gradient.start + i] = derivative;
                }
                continue;
            }
            let at = |id: NodeId| self.slots[id.0 as usize];
            for j in 0..out.len {
                tick(j, cx)?;
                let value = match self.problem.expr(node)? {
                    Expr::Var(_) => continue,
                    Expr::Const { value, .. } => *value,
                    Expr::Component { of, index } => values[at(*of).start + *index as usize],
                    Expr::Add(a, b) => values[at(*a).start + j] + values[at(*b).start + j],
                    Expr::Sub(a, b) => values[at(*a).start + j] - values[at(*b).start + j],
                    Expr::Mul(a, b) => {
                        let (a, b) = (at(*a), at(*b));
                        values[a.start + if a.len == 1 { 0 } else { j }] * values[b.start + if b.len == 1 { 0 } else { j }]
                    }
                    Expr::Div(a, b) => values[at(*a).start] / values[at(*b).start],
                    Expr::Neg(a) => -values[at(*a).start + j],
                    Expr::Powi { base, exp } => fs_math::det::powi(values[at(*base).start], *exp),
                    Expr::Sqrt(a) => fs_math::det::sqrt(values[at(*a).start]),
                    Expr::Exp(a) => fs_math::det::exp(values[at(*a).start]),
                    Expr::Ln(a) => fs_math::det::ln(values[at(*a).start]),
                    Expr::Tanh(a) => fs_math::det::tanh(values[at(*a).start]),
                    Expr::Dot(a, b) => {
                        let (a, b) = (at(*a), at(*b));
                        let mut sum = 0.0;
                        for k in 0..a.len { tick(k, cx)?; sum += values[a.start + k] * values[b.start + k]; }
                        sum
                    }
                    Expr::NormSq(a) => {
                        let a = at(*a);
                        let mut sum = 0.0;
                        for k in 0..a.len { tick(k, cx)?; let x = values[a.start + k]; sum += x * x; }
                        sum
                    }
                    _ => unreachable!("compiler refused non-algebraic or nonsmooth nodes"),
                };
                if !value.is_finite() {
                    return Err(OptError::EvalNonFinite { node: node.0, component: match self.problem.shape(node)? { Shape::Scalar => None, Shape::Vector(_) => Some(j as u32) }, bits: value.to_bits() }.into());
                }
                values[out.start + j] = value;
            }
        }
        let mut outputs = capacity(self.roots.len())?;
        for (i, root) in self.roots.iter().enumerate() { tick(i, cx)?; outputs.push(values[self.slots[root.0 as usize].start]); }
        poll(cx)?;
        Ok(ReverseEvaluation { program: self, values, outputs })
    }
}

/// Immutable primal evaluation. Multiple pullbacks reuse the same primal
/// values without re-running graph arithmetic or perturbing any variable.
#[derive(Debug)]
pub struct ReverseEvaluation<'program, 'problem> {
    program: &'program ReverseProgram<'problem>,
    values: Vec<f64>,
    outputs: Vec<f64>,
}

fn accumulate(adjoint: &mut [f64], slot: Slot, component: usize, contribution: f64, node: NodeId) -> Result<(), ReverseError> {
    let index = slot.start + component;
    let next = adjoint[index] + contribution;
    if !contribution.is_finite() || !next.is_finite() {
        return Err(ReverseError::NonFiniteAdjoint { node, component, bits: if contribution.is_finite() { next.to_bits() } else { contribution.to_bits() } });
    }
    adjoint[index] = next;
    Ok(())
}

impl ReverseEvaluation<'_, '_> {
    /// Scalar root values in the exact order passed to compilation.
    #[must_use]
    pub fn values(&self) -> &[f64] { &self.outputs }

    /// Return J^T * seeds in each variable's ambient point coordinates.
    /// Contributions from duplicate roots and shared children accumulate.
    /// Zero cotangents skip local derivative arithmetic, but do not suppress
    /// compile-time refusals or invalid primal values. Seeds must be finite.
    pub fn pullback(&self, seeds: &[f64], cx: Option<&Cx>) -> Result<Vec<Vec<f64>>, ReverseError> {
        poll(cx)?;
        if seeds.len() != self.program.roots.len() { return Err(ReverseError::SeedCount { expected: self.program.roots.len(), actual: seeds.len() }); }
        let mut adjoint = zeros(self.program.scalar_slots, cx)?;
        for (i, (&node, &seed)) in self.program.roots.iter().zip(seeds).enumerate() {
            tick(i, cx)?;
            accumulate(&mut adjoint, self.program.slots[node.0 as usize], 0, seed, node)?;
        }
        let values = &self.values;
        let at = |id: NodeId| self.program.slots[id.0 as usize];
        for &node in self.program.order.iter().rev() {
            poll(cx)?;
            let out = at(node);
            for j in 0..out.len {
                tick(j, cx)?;
                let g = adjoint[out.start + j];
                if g == 0.0 { continue; }
                match self.program.problem.expr(node)? {
                    Expr::Var(_) | Expr::Const { .. } => {},
                    Expr::PdeResidual { over, .. } => {
                        let index = self.program.physics.binary_search_by_key(&node, |binding| binding.node)
                            .expect("compiled PDE node has a bound provider");
                        let binding = &self.program.physics[index];
                        let input = self.program.variables[over.0 as usize];
                        for k in 0..input.len {
                            tick(k, cx)?;
                            accumulate(&mut adjoint, input, k, g * values[binding.gradient.start + k], node)?;
                        }
                    }
                    Expr::Component { of, index } => accumulate(&mut adjoint, at(*of), *index as usize, g, node)?,
                    Expr::Add(a, b) | Expr::Sub(a, b) => {
                        accumulate(&mut adjoint, at(*a), j, g, node)?;
                        let sign = if matches!(self.program.problem.expr(node)?, Expr::Sub(..)) { -1.0 } else { 1.0 };
                        accumulate(&mut adjoint, at(*b), j, sign * g, node)?;
                    }
                    Expr::Mul(a, b) => {
                        let (a, b) = (at(*a), at(*b));
                        let (ia, ib) = (if a.len == 1 { 0 } else { j }, if b.len == 1 { 0 } else { j });
                        accumulate(&mut adjoint, a, ia, g * values[b.start + ib], node)?;
                        accumulate(&mut adjoint, b, ib, g * values[a.start + ia], node)?;
                    }
                    Expr::Div(a, b) => {
                        let (a, b) = (at(*a), at(*b));
                        accumulate(&mut adjoint, a, 0, g / values[b.start], node)?;
                        accumulate(&mut adjoint, b, 0, -g * values[out.start] / values[b.start], node)?;
                    }
                    Expr::Neg(a) => accumulate(&mut adjoint, at(*a), j, -g, node)?,
                    Expr::Powi { base, exp } => {
                        let a = at(*base);
                        let derivative = match *exp {
                            0 => 0.0,
                            1 => 1.0,
                            i32::MIN => (*exp as f64) * (values[out.start] / values[a.start]),
                            _ => (*exp as f64) * fs_math::det::powi(values[a.start], *exp - 1),
                        };
                        accumulate(&mut adjoint, a, 0, g * derivative, node)?;
                    }
                    Expr::Sqrt(a) => accumulate(&mut adjoint, at(*a), 0, (0.5 * g) / values[out.start], node)?,
                    Expr::Exp(a) => accumulate(&mut adjoint, at(*a), 0, g * values[out.start], node)?,
                    Expr::Ln(a) => accumulate(&mut adjoint, at(*a), 0, g / values[at(*a).start], node)?,
                    Expr::Tanh(a) => {
                        // Avoid 1-tanh(x)^2 losing every derivative bit in
                        // saturated tails while the derivative is representable.
                        let t = fs_math::det::exp(-2.0 * values[at(*a).start].abs());
                        accumulate(&mut adjoint, at(*a), 0, g * (4.0 * t / ((1.0 + t) * (1.0 + t))), node)?;
                    }
                    Expr::Dot(a, b) => {
                        let (a, b) = (at(*a), at(*b));
                        for k in 0..a.len {
                            tick(k, cx)?;
                            accumulate(&mut adjoint, a, k, g * values[b.start + k], node)?;
                            accumulate(&mut adjoint, b, k, g * values[a.start + k], node)?;
                        }
                    }
                    Expr::NormSq(a) => {
                        let a = at(*a);
                        for k in 0..a.len { tick(k, cx)?; accumulate(&mut adjoint, a, k, g * (2.0 * values[a.start + k]), node)?; }
                    }
                    _ => unreachable!("compiler refused unsupported reverse rules"),
                }
            }
        }
        let mut result = capacity(self.program.variables.len())?;
        for slot in &self.program.variables {
            let mut gradient = capacity(slot.len)?;
            for j in 0..slot.len { tick(j, cx)?; gradient.push(adjoint[slot.start + j]); }
            result.push(gradient);
        }
        poll(cx)?;
        Ok(result)
    }

    /// Pull back further into each manifold's authoritative retraction
    /// parameters. In particular SO(3) maps four quaternion components to
    /// three body-frame parameters; Sphere/Stiefel use their live projections.
    pub fn parameter_pullback(&self, seeds: &[f64], cx: Option<&Cx>) -> Result<Vec<Vec<f64>>, ReverseError> {
        let ambient = self.pullback(seeds, cx)?;
        let mut result = capacity(ambient.len())?;
        for ((variable, slot), gradient) in self.program.problem.vars().iter().zip(&self.program.variables).zip(ambient) {
            poll(cx)?;
            result.push(variable.manifold.parameter_gradient(&self.values[slot.start..slot.start + slot.len], &gradient)?);
        }
        poll(cx)?;
        Ok(result)
    }
}

mod second_order;
pub use second_order::HessianError;
