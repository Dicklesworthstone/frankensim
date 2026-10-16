//! First-order directional sweep sharing the HVP primitive rules.
use super::{HessianError, ReverseError, ReverseEvaluation, edges, finite};
use super::super::{capacity, poll, tick, zeros};
use crate::{Expr, children};
use fs_exec::Cx;

impl ReverseEvaluation<'_, '_> {
    /// Apply the Jacobian of all compiled scalar roots to an ambient direction.
    /// One forward sweep reuses the immutable primal values and bound physics
    /// gradients; no primal solve, provider Hessian or dense Jacobian is needed.
    /// Roots (including duplicates) retain their compiled order. Direction and
    /// arithmetic refusals share the existing [`HessianError`] protocol.
    pub fn jacobian_vector_product(
        &self, direction: &[Vec<f64>], cx: Option<&Cx<'_>>,
    ) -> Result<Vec<f64>, HessianError> {
        self.selected_jacobian_vector_product(None, direction, cx)
    }

    pub(crate) fn selected_jacobian_vector_product(
        &self, selected: Option<&[bool]>, direction: &[Vec<f64>], cx: Option<&Cx<'_>>,
    ) -> Result<Vec<f64>, HessianError> {
        poll(cx)?;
        let p = self.program;
        if let Some(mask) = selected {
            if mask.len() != p.roots.len() {
                return Err(ReverseError::SeedCount { expected: p.roots.len(), actual: mask.len() }.into());
            }
        }
        if direction.len() != p.variables.len() {
            return Err(HessianError::DirectionCount { expected: p.variables.len(), actual: direction.len() });
        }
        for (variable, (block, slot)) in direction.iter().zip(&p.variables).enumerate() {
            poll(cx)?;
            if block.len() != slot.len {
                return Err(HessianError::DirectionLength { variable, expected: slot.len, actual: block.len() });
            }
            for (component, &value) in block.iter().enumerate() {
                tick(component, cx)?;
                if !value.is_finite() {
                    return Err(HessianError::NonFiniteDirection { variable, component, bits: value.to_bits() });
                }
            }
        }
        let mut active = capacity(p.slots.len())?;
        for i in 0..p.slots.len() { tick(i, cx)?; active.push(false); }
        for (i, &node) in p.roots.iter().enumerate() {
            tick(i, cx)?;
            if selected.is_none_or(|mask| mask[i]) { active[node.0 as usize] = true; }
        }
        for &node in p.order.iter().rev() {
            poll(cx)?;
            if active[node.0 as usize] {
                for child in children(p.problem.expr(node)?) { active[child.0 as usize] = true; }
            }
        }
        let mut tangent = zeros(p.scalar_slots, cx)?;
        for (slot, block) in p.variables.iter().zip(direction) {
            for (i, &v) in block.iter().enumerate() { tick(i, cx)?; tangent[slot.start + i] = v; }
        }
        for &node in &p.order {
            poll(cx)?;
            let expression = p.problem.expr(node)?;
            if !active[node.0 as usize] || matches!(expression, Expr::Var(_)) { continue; }
            let output = p.slots[node.0 as usize];
            if let Ok(index) = p.physics.binary_search_by_key(&node, |binding| binding.node) {
                let binding = &p.physics[index];
                let input = p.variables[binding.variable.0 as usize];
                let mut sum = 0.0;
                for i in 0..input.len {
                    tick(i, cx)?;
                    sum = finite(sum + self.values[binding.gradient.start + i] * tangent[input.start + i], node, i)?;
                }
                tangent[output.start] = sum;
                continue;
            }
            for lane in 0..output.len {
                tick(lane, cx)?;
                let mut sum = 0.0;
                edges(expression, lane, &p.slots, &self.values, &tangent, output, cx,
                    |child, component, first, _| {
                        let first = finite(first, node, lane)?;
                        sum = finite(sum + first * tangent[child.start + component], node, lane)?;
                        Ok(())
                    })?;
                tangent[output.start + lane] = sum;
            }
        }
        let mut result = capacity(p.roots.len())?;
        for (i, &node) in p.roots.iter().enumerate() {
            tick(i, cx)?;
            result.push(if selected.is_none_or(|mask| mask[i]) { tangent[p.slots[node.0 as usize].start] } else { 0.0 });
        }
        poll(cx)?;
        Ok(result)
    }
}
