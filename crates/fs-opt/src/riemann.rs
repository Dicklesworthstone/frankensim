//! Toy Riemannian gradient descent — NOT an optimizer offering (those are
//! later ASCENT beads); it exists to prove the manifold metadata is
//! consumable end-to-end: ambient gradient → tangent projection →
//! retraction, per variable, deterministically.

use crate::graph::{Problem, Unevaluable, VarId};
use std::collections::HashMap;

/// Fixed-step Riemannian descent on the problem's objective.
/// Deterministic: fixed iteration count, per-variable processing in
/// `VarId` order.
///
/// # Errors
/// Propagates [`Unevaluable`] from structure-only nodes.
pub fn descend(
    problem: &Problem,
    x0: HashMap<VarId, Vec<f64>>,
    step: f64,
    iters: usize,
) -> Result<(HashMap<VarId, Vec<f64>>, f64), Unevaluable> {
    let mut x = x0;
    for _ in 0..iters {
        let grads = problem.gradient(&x)?;
        let mut order: Vec<VarId> = x.keys().copied().collect();
        order.sort();
        for v in order {
            let man = problem.manifold(v).clone();
            let xt = x.get_mut(&v).expect("var present");
            let tangent = man.tangent_project(xt, &grads[&v]);
            let neg_step: Vec<f64> = tangent.iter().map(|g| -step * g).collect();
            *xt = man.retract(xt, &neg_step);
        }
    }
    let f = problem.eval(&x)?;
    Ok((x, f))
}
