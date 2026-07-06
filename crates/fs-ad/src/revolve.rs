//! Checkpointed reverse sweeps (plan §6.6): adjoint a length-L step chain
//! holding only O(log L) states instead of L. v1 implements BINARY
//! TREEVERSE (recursive bisection): peak snapshots ≤ ⌈log₂ L⌉ + 1, total
//! forward re-evaluations ≤ L·⌈log₂ L⌉. The binomially-optimal
//! Griewank–Walther schedule (fewer recomputations at a FIXED small
//! budget) and ledger spill for snapshots are recorded refinements.
//!
//! THE DETERMINISM PAYOFF (headline test): because every forward step is
//! deterministic, recomputed states are BIT-IDENTICAL to stored ones —
//! so the checkpointed adjoint equals the full-storage adjoint bitwise,
//! not just approximately. Checkpointing is a memory policy, invisible
//! to the numbers.

/// Instrumentation returned by a sweep: the resource-bound evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RevolveStats {
    /// Total forward step evaluations (including recomputation).
    pub forward_steps: u64,
    /// Peak simultaneously-held snapshots.
    pub peak_snapshots: usize,
}

/// Minimum snapshot budget binary treeverse needs for `steps` steps.
#[must_use]
pub fn min_budget(steps: usize) -> usize {
    (usize::BITS - steps.next_power_of_two().leading_zeros()) as usize
}

/// Reverse-sweep a deterministic step chain with checkpointing.
///
/// - `x0`: state before step 0.
/// - `steps`: chain length L (steps are indexed 0..L).
/// - `budget`: max snapshots held at once (≥ [`min_budget`]; structured
///   panic otherwise — silently exceeding memory budgets is the failure
///   mode this exists to prevent).
/// - `forward(i, x_i) -> x_{i+1}` — must be PURE and deterministic.
/// - `reverse(i, x_i, bar_{i+1}) -> bar_i` — pulls the adjoint back
///   across step i (the caller's vjp for that step).
/// - `seed`: the adjoint at the END of the chain (∂J/∂x_L).
///
/// Returns the adjoint at the start (∂J/∂x_0 shape, in the caller's `B`)
/// plus the instrumentation.
pub fn checkpointed_adjoint<S, B, F, R>(
    x0: &S,
    steps: usize,
    budget: usize,
    forward: &F,
    reverse: &R,
    seed: B,
) -> (B, RevolveStats)
where
    S: Clone,
    F: Fn(usize, &S) -> S,
    R: Fn(usize, &S, B) -> B,
{
    assert!(
        budget >= min_budget(steps),
        "snapshot budget {budget} < required {} for {steps} steps (binary treeverse)",
        min_budget(steps)
    );
    let mut stats = RevolveStats { forward_steps: 0, peak_snapshots: 0 };
    if steps == 0 {
        return (seed, stats);
    }
    let bar = sweep(x0, 0, steps, seed, forward, reverse, &mut stats, 1);
    (bar, stats)
}

/// Full-storage reference adjoint (stores every state) — the oracle the
/// checkpointed sweep must match BITWISE, and a fine choice when memory
/// is plentiful.
pub fn full_adjoint<S, B, F, R>(
    x0: &S,
    steps: usize,
    forward: &F,
    reverse: &R,
    seed: B,
) -> B
where
    S: Clone,
    F: Fn(usize, &S) -> S,
    R: Fn(usize, &S, B) -> B,
{
    let mut states = Vec::with_capacity(steps + 1);
    states.push(x0.clone());
    for i in 0..steps {
        let next = forward(i, &states[i]);
        states.push(next);
    }
    let mut bar = seed;
    for i in (0..steps).rev() {
        bar = reverse(i, &states[i], bar);
    }
    bar
}

#[allow(clippy::too_many_arguments)]
fn sweep<S, B, F, R>(
    state: &S,
    begin: usize,
    end: usize,
    bar: B,
    forward: &F,
    reverse: &R,
    stats: &mut RevolveStats,
    depth: usize,
) -> B
where
    S: Clone,
    F: Fn(usize, &S) -> S,
    R: Fn(usize, &S, B) -> B,
{
    stats.peak_snapshots = stats.peak_snapshots.max(depth);
    if end - begin == 1 {
        return reverse(begin, state, bar);
    }
    let mid = usize::midpoint(begin, end + 1); // ceil((begin+end)/2)
    // Advance a copy to mid while HOLDING `state` (the live snapshot at
    // this depth); the right half reverses first.
    let mut s_mid = state.clone();
    for i in begin..mid {
        s_mid = forward(i, &s_mid);
        stats.forward_steps += 1;
    }
    let bar = sweep(&s_mid, mid, end, bar, forward, reverse, stats, depth + 1);
    drop(s_mid);
    sweep(state, begin, mid, bar, forward, reverse, stats, depth)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dual::Dual64;

    /// Step: x ← x − 0.1·(x³ − θ). Adjoint state (x̄, θ̄).
    fn fwd(theta: f64) -> impl Fn(usize, &f64) -> f64 {
        move |_i, &x| (-0.1f64).mul_add(x * x * x - theta, x)
    }

    fn rev(theta: f64) -> impl Fn(usize, &f64, (f64, f64)) -> (f64, f64) {
        let _ = theta;
        move |_i, &x, (xbar, tbar)| {
            // x' = x − 0.1(x³ − θ): ∂x'/∂x = 1 − 0.3x², ∂x'/∂θ = 0.1.
            ((-0.3f64 * x * x + 1.0) * xbar, 0.1f64.mul_add(xbar, tbar))
        }
    }

    const STEPS: usize = 100;
    const X0: f64 = 0.3;
    const THETA: f64 = 0.8;

    #[test]
    fn checkpointed_equals_full_storage_bitwise() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let full = full_adjoint(&X0, STEPS, &f, &r, (1.0, 0.0));
        let budget = min_budget(STEPS);
        let (ck, stats) = checkpointed_adjoint(&X0, STEPS, budget, &f, &r, (1.0, 0.0));
        assert_eq!(full.0.to_bits(), ck.0.to_bits(), "xbar must be BIT-identical");
        assert_eq!(full.1.to_bits(), ck.1.to_bits(), "thetabar must be BIT-identical");
        // Resource bounds hold.
        assert!(
            stats.peak_snapshots <= budget,
            "peak {} exceeded budget {budget}",
            stats.peak_snapshots
        );
        let log2 = u64::from(usize::BITS - STEPS.leading_zeros());
        assert!(
            stats.forward_steps <= (STEPS as u64) * log2,
            "forward steps {} above the L*log2(L) bound",
            stats.forward_steps
        );
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"revolve\",\"verdict\":\"pass\",\"detail\":\"bitwise == full storage; {} fwd re-evals, peak {} snaps (budget {budget})\"}}",
            stats.forward_steps, stats.peak_snapshots
        );
    }

    #[test]
    fn gradient_matches_forward_duals() {
        // dx_L/dθ via the adjoint sweep must equal the forward-dual run.
        let f = fwd(THETA);
        let r = rev(THETA);
        let (bar, _) =
            checkpointed_adjoint(&X0, STEPS, min_budget(STEPS), &f, &r, (1.0, 0.0));
        // Forward dual: seed θ.
        let mut x = Dual64::<1>::constant(X0);
        let th = Dual64::<1>::variable(THETA, 0);
        for _ in 0..STEPS {
            x = x - Dual64::constant(0.1) * (x * x * x - th);
        }
        assert!(
            (bar.1 - x.eps[0]).abs() < 1e-13 * x.eps[0].abs().max(1.0),
            "adjoint dθ {} vs forward dual {}",
            bar.1,
            x.eps[0]
        );
        // And the x̄ channel: dx_L/dx_0.
        let mut y = Dual64::<1>::variable(X0, 0);
        let thc = Dual64::<1>::constant(THETA);
        for _ in 0..STEPS {
            y = y - Dual64::constant(0.1) * (y * y * y - thc);
        }
        assert!(
            (bar.0 - y.eps[0]).abs() < 1e-13 * y.eps[0].abs().max(1.0),
            "adjoint dx0 {} vs forward dual {}",
            bar.0,
            y.eps[0]
        );
    }

    #[test]
    fn insufficient_budget_is_refused() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let res = std::panic::catch_unwind(|| {
            checkpointed_adjoint(&X0, 64, 2, &f, &r, (1.0, 0.0))
        });
        assert!(res.is_err(), "budget 2 for 64 steps must be refused loudly");
    }

    #[test]
    fn zero_and_one_step_chains() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let (bar, stats) = checkpointed_adjoint(&X0, 0, 1, &f, &r, (2.0, 3.0));
        assert_eq!(bar, (2.0, 3.0), "empty chain returns the seed");
        assert_eq!(stats.forward_steps, 0);
        let (bar1, _) = checkpointed_adjoint(&X0, 1, 1, &f, &r, (1.0, 0.0));
        let full1 = full_adjoint(&X0, 1, &f, &r, (1.0, 0.0));
        assert_eq!(bar1.0.to_bits(), full1.0.to_bits());
    }
}
