//! Directional differentiation of the reverse sweep (forward-over-reverse).
//! The primal tape is reused: no finite differences and no dense Hessian.

use super::{ReverseError, ReverseEvaluation, Slot, accumulate, capacity, poll, tick, zeros};
use crate::{Expr, NodeId, OptError, children};
use fs_exec::Cx;

/// Refusal from a mathematical ambient Hessian-vector product.
#[derive(Debug, Clone, PartialEq)]
pub enum HessianError {
    /// Original graph, seed, allocation, cancellation or adjoint error.
    Reverse(ReverseError),
    /// Directions require one block for every declared variable.
    DirectionCount { expected: usize, actual: usize },
    /// A direction block must match its variable's point storage.
    DirectionLength { variable: usize, expected: usize, actual: usize },
    /// A direction component was not finite.
    NonFiniteDirection { variable: usize, component: usize, bits: u64 },
    /// A local derivative or directional accumulation was not finite.
    NonFiniteDerivative { node: NodeId, component: usize, bits: u64 },
}

impl From<ReverseError> for HessianError {
    fn from(error: ReverseError) -> Self { Self::Reverse(error) }
}
impl From<OptError> for HessianError {
    fn from(error: OptError) -> Self { Self::Reverse(error.into()) }
}
impl core::fmt::Display for HessianError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Reverse(error) => write!(f, "{error}"),
            Self::DirectionCount { expected, actual } => write!(f, "Hessian product needs {expected} direction blocks, received {actual}"),
            Self::DirectionLength { variable, expected, actual } => write!(f, "Hessian direction {variable} needs {expected} coordinates, received {actual}"),
            Self::NonFiniteDirection { variable, component, bits } => write!(f, "Hessian direction {variable}:{component} is non-finite ({bits:#018x})"),
            Self::NonFiniteDerivative { node, component, bits } => write!(f, "node {} has a non-finite second-order derivative at {component} ({bits:#018x})", node.0),
        }
    }
}
impl std::error::Error for HessianError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Reverse(error) => Some(error), _ => None }
    }
}

fn finite(value: f64, node: NodeId, component: usize) -> Result<f64, HessianError> {
    if value.is_finite() { Ok(value) } else {
        Err(HessianError::NonFiniteDerivative { node, component, bits: value.to_bits() })
    }
}

/// Emit each scalar dependency with its first partial and directional partial.
/// One rule is shared by the forward and differentiated reverse sweeps.
fn edges(
    expression: &Expr, lane: usize, slots: &[Slot], values: &[f64], tangent: &[f64],
    output: Slot, cx: Option<&Cx<'_>>,
    mut emit: impl FnMut(Slot, usize, f64, f64) -> Result<(), HessianError>,
) -> Result<(), HessianError> {
    let at = |id: NodeId| slots[id.0 as usize];
    match expression {
        Expr::Var(_) | Expr::Const { .. } => {}
        Expr::Component { of, index } => emit(at(*of), *index as usize, 1.0, 0.0)?,
        Expr::Add(a, b) | Expr::Sub(a, b) => {
            emit(at(*a), lane, 1.0, 0.0)?;
            emit(at(*b), lane, if matches!(expression, Expr::Sub(..)) { -1.0 } else { 1.0 }, 0.0)?;
        }
        Expr::Neg(a) => emit(at(*a), lane, -1.0, 0.0)?,
        Expr::Mul(a, b) => {
            let (a, b) = (at(*a), at(*b));
            let (i, j) = (if a.len == 1 { 0 } else { lane }, if b.len == 1 { 0 } else { lane });
            emit(a, i, values[b.start + j], tangent[b.start + j])?;
            emit(b, j, values[a.start + i], tangent[a.start + i])?;
        }
        Expr::Div(a, b) => {
            let (a, b) = (at(*a), at(*b));
            let inverse = 1.0 / values[b.start];
            let relative = tangent[b.start] / values[b.start];
            let quotient = values[output.start];
            emit(a, 0, inverse, -relative * inverse)?;
            emit(b, 0, -quotient * inverse,
                (2.0 * quotient * relative - tangent[a.start] * inverse) * inverse)?;
        }
        Expr::Dot(a, b) => {
            let (a, b) = (at(*a), at(*b));
            for i in 0..a.len {
                tick(i, cx)?;
                emit(a, i, values[b.start + i], tangent[b.start + i])?;
                emit(b, i, values[a.start + i], tangent[a.start + i])?;
            }
        }
        Expr::NormSq(a) => {
            let a = at(*a);
            for i in 0..a.len {
                tick(i, cx)?;
                emit(a, i, 2.0 * values[a.start + i], 2.0 * tangent[a.start + i])?;
            }
        }
        Expr::Powi { base, exp } => {
            let a = at(*base);
            let x = values[a.start];
            let p = f64::from(*exp);
            // checked_sub keeps i32::MIN and i32::MIN+1 legal exponents.
            let first = match *exp {
                0 => 0.0, 1 => 1.0,
                _ => p * exp.checked_sub(1).map_or_else(
                    || values[output.start] / x, |power| fs_math::det::powi(x, power)),
            };
            let second = match *exp {
                0 | 1 => 0.0, 2 => 2.0,
                _ => p * (p - 1.0) * exp.checked_sub(2).map_or_else(
                    || (values[output.start] / x) / x, |power| fs_math::det::powi(x, power)),
            };
            emit(a, 0, first, second * tangent[a.start])?;
        }
        Expr::Sqrt(a) | Expr::Exp(a) | Expr::Ln(a) | Expr::Tanh(a) => {
            let a = at(*a);
            let x = values[a.start];
            let y = values[output.start];
            let (first, second) = match expression {
                Expr::Sqrt(_) => { let d = 0.5 / y; (d, (-0.5 * d) / x) }
                Expr::Exp(_) => (y, y),
                Expr::Ln(_) => { let d = 1.0 / x; (d, -d / x) }
                Expr::Tanh(_) => {
                    let t = fs_math::det::exp(-2.0 * x.abs());
                    let d = 4.0 * t / ((1.0 + t) * (1.0 + t));
                    (d, -2.0 * y * d)
                }
                _ => unreachable!(),
            };
            emit(a, 0, first, second * tangent[a.start])?;
        }
        _ => unreachable!("reverse compiler refused unsupported operations"),
    }
    Ok(())
}

impl ReverseEvaluation<'_, '_> {
    /// Apply the Hessian of `sum(seeds[i] * root[i])` to an ambient direction.
    /// Seeds are held constant. Directions and results use variable declaration
    /// order and point storage, NOT manifold retraction parameters. This is not
    /// a Riemannian Hessian: differentiating the manifold pullback is separate.
    ///
    /// Reuses this immutable primal tape. One directional forward sweep and one
    /// differentiated reverse sweep require three scalar workspaces of the
    /// program's scalar_slots size, plus a node mask and the output. Neither an
    /// objective re-evaluation nor an n-by-n Hessian is constructed. The caller
    /// owns product-count budgets; numerical loops poll cancellation.
    ///
    /// Active primitives must have finite first and directional second partials.
    /// Unsupported nodes retain compile-time refusals; zero-seeded subgraphs
    /// need no derivative arithmetic. These are chain-rule derivatives of the
    /// mathematical operators, not interval certificates or derivatives of the
    /// floating-point implementation. A refusal publishes no partial product.
    pub fn hessian_vector_product(
        &self, seeds: &[f64], direction: &[Vec<f64>], cx: Option<&Cx<'_>>,
    ) -> Result<Vec<Vec<f64>>, HessianError> {
        poll(cx)?;
        let program = self.program;
        if seeds.len() != program.roots.len() {
            return Err(ReverseError::SeedCount { expected: program.roots.len(), actual: seeds.len() }.into());
        }
        if direction.len() != program.variables.len() {
            return Err(HessianError::DirectionCount { expected: program.variables.len(), actual: direction.len() });
        }
        for (variable, (block, slot)) in direction.iter().zip(&program.variables).enumerate() {
            poll(cx)?;
            if block.len() != slot.len {
                return Err(HessianError::DirectionLength { variable, expected: slot.len, actual: block.len() });
            }
            for (component, value) in block.iter().enumerate() {
                tick(component, cx)?;
                if !value.is_finite() {
                    return Err(HessianError::NonFiniteDirection { variable, component, bits: value.to_bits() });
                }
            }
        }
        let mut adjoint = zeros(program.scalar_slots, cx)?;
        for (i, (&node, &seed)) in program.roots.iter().zip(seeds).enumerate() {
            tick(i, cx)?;
            accumulate(&mut adjoint, program.slots[node.0 as usize], 0, seed, node)?;
        }
        let mut active = capacity(program.slots.len())?;
        for i in 0..program.slots.len() { tick(i, cx)?; active.push(false); }
        for (i, &node) in program.roots.iter().enumerate() {
            tick(i, cx)?;
            active[node.0 as usize] = adjoint[program.slots[node.0 as usize].start] != 0.0;
        }
        for &node in program.order.iter().rev() {
            poll(cx)?;
            if active[node.0 as usize] {
                if matches!(program.problem.expr(node)?, Expr::PdeResidual { .. }) {
                    return Err(OptError::Unevaluable {
                        node: node.0,
                        kind: "bound physics supplies first-order derivatives, not a Hessian action",
                    }.into());
                }
                for child in children(program.problem.expr(node)?) { active[child.0 as usize] = true; }
            }
        }
        let mut tangent = zeros(program.scalar_slots, cx)?;
        let mut product = zeros(program.scalar_slots, cx)?;
        for (slot, block) in program.variables.iter().zip(direction) {
            for (i, &value) in block.iter().enumerate() { tick(i, cx)?; tangent[slot.start + i] = value; }
        }
        for &node in &program.order {
            poll(cx)?;
            if !active[node.0 as usize] || matches!(program.problem.expr(node)?, Expr::Var(_)) { continue; }
            let output = program.slots[node.0 as usize];
            for lane in 0..output.len {
                tick(lane, cx)?;
                let mut sum = 0.0;
                edges(program.problem.expr(node)?, lane, &program.slots, &self.values, &tangent, output, cx,
                    |child, component, first, _| {
                        finite(first, node, lane)?;
                        sum = finite(sum + first * tangent[child.start + component], node, lane)?;
                        Ok(())
                    })?;
                tangent[output.start + lane] = sum;
            }
        }
        for &node in program.order.iter().rev() {
            poll(cx)?;
            if !active[node.0 as usize] { continue; }
            let output = program.slots[node.0 as usize];
            for lane in 0..output.len {
                tick(lane, cx)?;
                let bar = adjoint[output.start + lane];
                let dot_bar = product[output.start + lane];
                if bar == 0.0 && dot_bar == 0.0 { continue; }
                edges(program.problem.expr(node)?, lane, &program.slots, &self.values, &tangent, output, cx,
                    |child, component, first, dot_first| {
                        let first = finite(first, node, lane)?;
                        // Zero adjoints must not manufacture 0 * infinity.
                        let second_term = if bar == 0.0 { 0.0 } else { bar * finite(dot_first, node, lane)? };
                        accumulate(&mut product, child, component,
                            finite(dot_bar * first + second_term, node, lane)?, node)?;
                        accumulate(&mut adjoint, child, component, bar * first, node)?;
                        Ok(())
                    })?;
            }
        }
        let mut result = capacity(program.variables.len())?;
        for slot in &program.variables {
            let mut block = capacity(slot.len)?;
            for lane in 0..slot.len { tick(lane, cx)?; block.push(product[slot.start + lane]); }
            result.push(block);
        }
        poll(cx)?;
        Ok(result)
    }
}