//! Checkpointed reverse sweeps (plan §6.6): adjoint a length-L step chain
//! holding only O(log L) states instead of L. v1 implements BINARY
//! TREEVERSE (recursive bisection): peak snapshots are ⌊log₂ L⌋ + 1 for
//! L > 0 (zero for an empty chain), and total forward re-evaluations are
//! ≤ L·⌈log₂ L⌉. Bead o3ui adds the
//! binomially-OPTIMAL Griewank–Walther schedule
//! ([`checkpointed_adjoint_binomial`]: recomputations exactly
//! r·L − β(s+1, r−1), the proven minimum for a fixed snapshot budget)
//! and the snapshot SPILL hook ([`checkpointed_adjoint_spilling`]:
//! RAM budget exceeded → overflow snapshots go to a caller-provided
//! [`SnapshotStore`] — e.g. an fs-ledger adapter at a higher layer —
//! instead of refusing; the no-store default keeps the structured
//! panic).
//!
//! THE DETERMINISM PAYOFF (headline test): because every forward step is
//! deterministic, recomputed states are BIT-IDENTICAL to stored ones —
//! so the checkpointed adjoint equals the full-storage adjoint bitwise,
//! not just approximately. Checkpointing is a memory policy, invisible
//! to the numbers. The same holds for spilled snapshots, PROVIDED the
//! store round-trips byte-exactly (the trait's contract).

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
    (usize::BITS - steps.leading_zeros()) as usize
}

fn ceil_midpoint(begin: usize, end: usize) -> usize {
    debug_assert!(begin <= end);
    let span = end - begin;
    begin + span / 2 + span % 2
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
    let mut stats = RevolveStats {
        forward_steps: 0,
        peak_snapshots: 0,
    };
    if steps == 0 {
        return (seed, stats);
    }
    let bar = sweep(x0, 0, steps, seed, forward, reverse, &mut stats, 1);
    (bar, stats)
}

/// Full-storage reference adjoint (stores every state) — the oracle the
/// checkpointed sweep must match BITWISE, and a fine choice when memory
/// is plentiful.
pub fn full_adjoint<S, B, F, R>(x0: &S, steps: usize, forward: &F, reverse: &R, seed: B) -> B
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
    let mid = ceil_midpoint(begin, end);
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

// ---------------------------------------------------------------------------
// Binomial (Griewank–Walther) schedule.
// ---------------------------------------------------------------------------

/// Binomial capacity β(s, r) = C(s+r, s): the longest chain reversible
/// with s snapshots and at most r forward sweeps over any step
/// (saturating — capacities beyond u64 are "unbounded" for any real L).
#[must_use]
pub fn beta(s: usize, r: usize) -> u64 {
    // C(s+r, min(s, r)) with saturating arithmetic.
    let k = s.min(r) as u64;
    let n = (s + r) as u64;
    let mut acc = 1u128;
    for i in 1..=k {
        acc = acc.saturating_mul(u128::from(n - k + i)) / u128::from(i);
        if acc > u128::from(u64::MAX) {
            return u64::MAX;
        }
    }
    u64::try_from(acc).unwrap_or(u64::MAX)
}

/// Minimal sweep count r with β(s, r) ≥ l (the "repetition number").
fn reps(l: u64, s: usize) -> usize {
    assert!(s >= 1, "binomial schedule needs at least one snapshot");
    let mut r = 0usize;
    while beta(s, r) < l {
        r += 1;
    }
    r
}

/// The PROVEN minimum forward re-evaluations for reversing `steps`
/// steps with `budget` snapshots: r·L − β(s+1, r−1), where r is minimal
/// with β(s, r) ≥ L (Griewank & Walther, Algorithm 799). The binomial
/// battery asserts the measured count EQUALS this — optimality is a
/// gate, not a hope.
#[must_use]
pub fn binomial_reevals(steps: usize, budget: usize) -> u64 {
    let l = steps as u64;
    if l <= 1 {
        return 0;
    }
    let r = reps(l, budget);
    if r == 0 {
        return 0;
    }
    (r as u64) * l - beta(budget + 1, r - 1)
}

/// Reverse-sweep with the binomially-optimal checkpoint schedule.
///
/// Same contract as [`checkpointed_adjoint`], but the snapshot budget
/// is a FIXED cap independent of L (any `budget ≥ 1` works; smaller
/// budgets trade re-evaluations, down to the s = 1 quadratic sweep),
/// and the forward re-evaluation count achieves the binomial optimum
/// r·L − β(s+1, r−1) exactly.
///
/// `budget` = s counts PARKED checkpoints (β's semantics: states
/// retained while other segments are swept, including x0). The live
/// sweep state is one more, so RAM holds at most `budget + 1` states
/// simultaneously — that is what [`RevolveStats::peak_snapshots`]
/// reports (all alive states, same truth as the treeverse counter).
///
/// # Panics
/// If `budget == 0` with `steps ≥ 1`.
pub fn checkpointed_adjoint_binomial<S, B, F, R>(
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
    let mut stats = RevolveStats {
        forward_steps: 0,
        peak_snapshots: 0,
    };
    if steps == 0 {
        return (seed, stats);
    }
    assert!(budget >= 1, "binomial schedule needs snapshot budget >= 1");
    let bar = binomial_sweep(x0, 0, steps, budget, seed, forward, reverse, &mut stats, 1);
    (bar, stats)
}

/// Theorem cost t(l, s) = r·l − β(s+1, r−1) — the proven minimum
/// (infeasible-sentinel for l ≥ 2 with no snapshots).
fn t_opt(l: u64, s: usize) -> u64 {
    if l <= 1 {
        return 0;
    }
    if s == 0 {
        return u64::MAX / 4;
    }
    let r = reps(l, s);
    (r as u64) * l - beta(s + 1, r - 1)
}

/// The optimal split via the DP recurrence
/// t(l, s) = min over l̂ of [l̂ + t(l̂, s) + t(l−l̂, s−1)]: the cost is
/// piecewise-linear in l̂ with breakpoints exactly at the β values
/// (where a subsegment's repetition class changes), so scanning the
/// O(r) breakpoint candidates against the CLOSED-FORM t recovers the
/// exact optimum. (A naive "l̂ = β(s, r−1)" split was measured 10
/// re-evals above optimal at s = 2, L = 100 — subsegments can land in
/// lower repetition classes than the parent's r−1/r.)
fn binomial_split(l: u64, s: usize) -> u64 {
    let r = reps(l, s);
    let mut best_cost = u64::MAX;
    let mut best = 1u64;
    let consider = |lhat: u64, best_cost: &mut u64, best: &mut u64| {
        if lhat >= 1 && lhat < l {
            let cost = lhat
                .saturating_add(t_opt(lhat, s))
                .saturating_add(t_opt(l - lhat, s - 1));
            if cost < *best_cost || (cost == *best_cost && lhat < *best) {
                *best_cost = cost;
                *best = lhat;
            }
        }
    };
    for rho in 0..=r {
        consider(beta(s, rho).min(l - 1), &mut best_cost, &mut best);
        consider(
            l.saturating_sub(beta(s - 1, rho)).max(1),
            &mut best_cost,
            &mut best,
        );
    }
    best
}

#[allow(clippy::too_many_arguments)]
fn binomial_sweep<S, B, F, R>(
    state: &S,
    begin: usize,
    end: usize,
    slots: usize,
    bar: B,
    forward: &F,
    reverse: &R,
    stats: &mut RevolveStats,
    held: usize,
) -> B
where
    S: Clone,
    F: Fn(usize, &S) -> S,
    R: Fn(usize, &S, B) -> B,
{
    stats.peak_snapshots = stats.peak_snapshots.max(held);
    let l = (end - begin) as u64;
    if l == 1 {
        return reverse(begin, state, bar);
    }
    let lhat = usize::try_from(binomial_split(l, slots)).expect("split fits usize");
    let mid = begin + lhat;
    // Advance a copy to mid while holding `state`; the mid snapshot
    // occupies one more slot, so the right segment sees slots − 1.
    let mut s_mid = state.clone();
    for i in begin..mid {
        s_mid = forward(i, &s_mid);
        stats.forward_steps += 1;
    }
    let bar = binomial_sweep(
        &s_mid,
        mid,
        end,
        slots - 1,
        bar,
        forward,
        reverse,
        stats,
        held + 1,
    );
    drop(s_mid);
    // The left segment reuses this level's slot set (mid was dropped).
    binomial_sweep(state, begin, mid, slots, bar, forward, reverse, stats, held)
}

// ---------------------------------------------------------------------------
// Snapshot spill (ledger hook).
// ---------------------------------------------------------------------------

/// External snapshot storage for sweeps whose RAM budget is smaller
/// than the schedule needs. CONTRACT: `read(key)` must return exactly
/// the bytes passed to `write(key, ..)` — byte-exact round-trip is what
/// keeps the checkpointed adjoint bitwise-equal to full storage. Keys
/// are written once and may be read MULTIPLE times; `evict` is a hint
/// that a key is dead (default no-op). The fsqlite-backed adapter
/// belongs to fs-ledger (L6) — this trait is the L1 seam.
pub trait SnapshotStore<S> {
    /// Persist a snapshot under `key` (keys are unique per sweep).
    fn write(&mut self, key: u64, state: &S);
    /// Recall the snapshot stored under `key`.
    fn read(&mut self, key: u64) -> S;
    /// The snapshot under `key` is no longer needed (hint; no-op default).
    fn evict(&mut self, key: u64) {
        let _ = key;
    }
}

/// Instrumentation for a spilling sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpillStats {
    /// The in-RAM sweep counters.
    pub revolve: RevolveStats,
    /// Peak snapshots held IN RAM (≤ the RAM budget, gated).
    pub peak_ram_snapshots: usize,
    /// Snapshots written to the store.
    pub spills: u64,
    /// Snapshots read back from the store.
    pub restores: u64,
}

/// Where a level's base snapshot lives.
enum Slot<S> {
    Ram(S),
    Spilled(u64),
}

/// Reverse-sweep with binary treeverse, spilling snapshots BEYOND
/// `ram_budget` to `store` instead of refusing. `ram_budget ≥ 1` (the
/// root state itself); every deeper level that exceeds the budget
/// round-trips through the store. Bitwise equality with full storage
/// holds whenever the store honors its byte-exact contract (gated in
/// the battery with an in-memory store).
///
/// Use [`checkpointed_adjoint`] when the budget suffices — the no-spill
/// structured panic remains the default failure mode; this entry point
/// is the OPT-IN escape valve.
///
/// # Panics
/// If `ram_budget == 0` with `steps ≥ 1`.
pub fn checkpointed_adjoint_spilling<S, B, F, R>(
    x0: &S,
    steps: usize,
    ram_budget: usize,
    store: &mut dyn SnapshotStore<S>,
    forward: &F,
    reverse: &R,
    seed: B,
) -> (B, SpillStats)
where
    S: Clone,
    F: Fn(usize, &S) -> S,
    R: Fn(usize, &S, B) -> B,
{
    let mut stats = SpillStats {
        revolve: RevolveStats {
            forward_steps: 0,
            peak_snapshots: 0,
        },
        peak_ram_snapshots: 0,
        spills: 0,
        restores: 0,
    };
    if steps == 0 {
        return (seed, stats);
    }
    assert!(
        ram_budget >= 1,
        "spilling sweep needs RAM for the live state"
    );
    let mut next_key = 0u64;
    let slot = Slot::Ram(x0.clone());
    let bar = spill_sweep(
        &slot,
        0,
        steps,
        seed,
        forward,
        reverse,
        &mut stats,
        1,
        ram_budget,
        store,
        &mut next_key,
    );
    (bar, stats)
}

#[allow(clippy::too_many_arguments)]
fn spill_sweep<S, B, F, R>(
    slot: &Slot<S>,
    begin: usize,
    end: usize,
    bar: B,
    forward: &F,
    reverse: &R,
    stats: &mut SpillStats,
    ram_held: usize,
    ram_budget: usize,
    store: &mut dyn SnapshotStore<S>,
    next_key: &mut u64,
) -> B
where
    S: Clone,
    F: Fn(usize, &S) -> S,
    R: Fn(usize, &S, B) -> B,
{
    stats.revolve.peak_snapshots = stats.revolve.peak_snapshots.max(ram_held);
    stats.peak_ram_snapshots = stats.peak_ram_snapshots.max(ram_held);
    // Materialize this level's base state.
    let state: S = match slot {
        Slot::Ram(s) => s.clone(),
        Slot::Spilled(key) => {
            stats.restores += 1;
            store.read(*key)
        }
    };
    if end - begin == 1 {
        return reverse(begin, &state, bar);
    }
    let mid = ceil_midpoint(begin, end);
    let mut s_mid = state;
    for i in begin..mid {
        s_mid = forward(i, &s_mid);
        stats.revolve.forward_steps += 1;
    }
    // The mid snapshot takes a slot: RAM while within budget, else spill.
    let (mid_slot, child_ram) = if ram_held < ram_budget {
        (Slot::Ram(s_mid), ram_held + 1)
    } else {
        let key = *next_key;
        *next_key += 1;
        stats.spills += 1;
        store.write(key, &s_mid);
        drop(s_mid);
        (Slot::Spilled(key), ram_held)
    };
    stats.peak_ram_snapshots = stats.peak_ram_snapshots.max(child_ram);
    let bar = spill_sweep(
        &mid_slot, mid, end, bar, forward, reverse, stats, child_ram, ram_budget, store, next_key,
    );
    if let Slot::Spilled(key) = mid_slot {
        store.evict(key);
    }
    drop(mid_slot);
    spill_sweep(
        slot, begin, mid, bar, forward, reverse, stats, ram_held, ram_budget, store, next_key,
    )
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
    fn minimum_budget_matches_binary_treeverse_depth() {
        for (steps, expected) in [
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 2),
            (4, 3),
            (5, 3),
            (7, 3),
            (8, 4),
            (9, 4),
            (usize::MAX, usize::BITS as usize),
        ] {
            assert_eq!(min_budget(steps), expected, "steps={steps}");
        }
        assert_eq!(
            ceil_midpoint(0, usize::MAX),
            usize::MAX / 2 + 1,
            "the first split is overflow-free at the full usize boundary"
        );
        assert_eq!(ceil_midpoint(usize::MAX - 2, usize::MAX), usize::MAX - 1);
    }

    #[test]
    fn zero_step_chain_needs_no_budget_or_callbacks() {
        let forward_calls = std::cell::Cell::new(0usize);
        let reverse_calls = std::cell::Cell::new(0usize);
        let forward = |_: usize, state: &f64| {
            forward_calls.set(forward_calls.get() + 1);
            *state
        };
        let reverse = |_: usize, _: &f64, bar: (f64, f64)| {
            reverse_calls.set(reverse_calls.get() + 1);
            bar
        };

        let seed = (2.0, 3.0);
        let (bar, stats) = checkpointed_adjoint(&X0, 0, 0, &forward, &reverse, seed);
        assert_eq!(bar, seed);
        assert_eq!(
            stats,
            RevolveStats {
                forward_steps: 0,
                peak_snapshots: 0,
            }
        );
        assert_eq!(forward_calls.get(), 0);
        assert_eq!(reverse_calls.get(), 0);
    }

    #[test]
    fn non_power_of_two_chains_run_at_the_exact_minimum() {
        let forward = fwd(THETA);
        let reverse = rev(THETA);
        for steps in [3usize, 5] {
            let full = full_adjoint(&X0, steps, &forward, &reverse, (1.0, 0.0));
            let budget = min_budget(steps);
            let (checkpointed, stats) =
                checkpointed_adjoint(&X0, steps, budget, &forward, &reverse, (1.0, 0.0));
            assert_eq!(full.0.to_bits(), checkpointed.0.to_bits());
            assert_eq!(full.1.to_bits(), checkpointed.1.to_bits());
            assert_eq!(stats.peak_snapshots, budget, "steps={steps}");
        }
    }

    #[test]
    fn checkpointed_equals_full_storage_bitwise() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let full = full_adjoint(&X0, STEPS, &f, &r, (1.0, 0.0));
        let budget = min_budget(STEPS);
        let (ck, stats) = checkpointed_adjoint(&X0, STEPS, budget, &f, &r, (1.0, 0.0));
        assert_eq!(
            full.0.to_bits(),
            ck.0.to_bits(),
            "xbar must be BIT-identical"
        );
        assert_eq!(
            full.1.to_bits(),
            ck.1.to_bits(),
            "thetabar must be BIT-identical"
        );
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
        let (bar, _) = checkpointed_adjoint(&X0, STEPS, min_budget(STEPS), &f, &r, (1.0, 0.0));
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
        let res = std::panic::catch_unwind(|| checkpointed_adjoint(&X0, 64, 2, &f, &r, (1.0, 0.0)));
        assert!(res.is_err(), "budget 2 for 64 steps must be refused loudly");
    }

    #[test]
    fn binomial_is_bitwise_and_exactly_optimal() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let full = full_adjoint(&X0, STEPS, &f, &r, (1.0, 0.0));
        // Across a range of FIXED budgets — including ones far below
        // binary treeverse's requirement — the sweep stays bitwise and
        // hits the binomial optimum r·L − β(s+1, r−1) EXACTLY.
        for budget in [1usize, 2, 3, 5, 8] {
            let (ck, stats) = checkpointed_adjoint_binomial(&X0, STEPS, budget, &f, &r, (1.0, 0.0));
            assert_eq!(full.0.to_bits(), ck.0.to_bits(), "xbar at budget {budget}");
            assert_eq!(full.1.to_bits(), ck.1.to_bits(), "tbar at budget {budget}");
            // s parked checkpoints + the live sweep state.
            assert!(
                stats.peak_snapshots <= budget + 1,
                "alive peak {} exceeded s+1 = {}",
                stats.peak_snapshots,
                budget + 1
            );
            let optimal = binomial_reevals(STEPS, budget);
            assert_eq!(
                stats.forward_steps, optimal,
                "budget {budget}: measured {} vs binomial optimum {optimal}",
                stats.forward_steps
            );
        }
        // RAM-fair comparison: treeverse holding `tb` alive states vs
        // binomial with s = tb − 1 parked (+1 live = same RAM). The
        // optimal schedule must not lose.
        let tb = min_budget(STEPS);
        let (_, bin) = checkpointed_adjoint_binomial(&X0, STEPS, tb - 1, &f, &r, (1.0, 0.0));
        let (_, tree) = checkpointed_adjoint(&X0, STEPS, tb, &f, &r, (1.0, 0.0));
        assert_eq!(tb, 7, "100-step treeverse has an exact seven-state peak");
        assert_eq!(tree.peak_snapshots, tb);
        assert_eq!(bin.forward_steps, binomial_reevals(STEPS, tb - 1));
        assert_eq!(bin.forward_steps, 280);
        assert_eq!(tree.forward_steps, 356);
        assert!(
            bin.forward_steps <= tree.forward_steps,
            "binomial {} must be <= treeverse {} at equal RAM {tb}",
            bin.forward_steps,
            tree.forward_steps
        );
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"revolve-binomial\",\"verdict\":\"pass\",\"detail\":\"bitwise across budgets 1..8; counts == r*L - beta(s+1,r-1) exactly; equal-RAM {tb}: binomial {} vs treeverse {}\"}}",
            bin.forward_steps, tree.forward_steps
        );
    }

    #[test]
    fn binomial_beta_table() {
        // Known values: beta(s, r) = C(s+r, s).
        assert_eq!(beta(1, 3), 4);
        assert_eq!(beta(2, 2), 6);
        assert_eq!(beta(3, 4), 35);
        assert_eq!(beta(0, 5), 1);
        assert_eq!(beta(5, 0), 1);
        // s = 1 degenerates to the quadratic single-snapshot sweep.
        assert_eq!(binomial_reevals(10, 1), 45); // L(L-1)/2
    }

    /// Byte-exact in-memory store that COUNTS its traffic (the mock for
    /// the fs-ledger adapter, which must honor the same contract).
    struct MapStore {
        map: std::collections::HashMap<u64, Vec<u8>>,
        evictions: u64,
    }

    impl SnapshotStore<f64> for MapStore {
        fn write(&mut self, key: u64, state: &f64) {
            self.map.insert(key, state.to_le_bytes().to_vec());
        }
        fn read(&mut self, key: u64) -> f64 {
            let b = self.map.get(&key).expect("spilled key present");
            f64::from_le_bytes(b.as_slice().try_into().expect("8 bytes"))
        }
        fn evict(&mut self, key: u64) {
            self.map.remove(&key);
            self.evictions += 1;
        }
    }

    #[test]
    fn spilling_sweep_is_bitwise_under_tiny_ram() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let full = full_adjoint(&X0, STEPS, &f, &r, (1.0, 0.0));
        // RAM budget 2 for a 100-step chain: binary treeverse would
        // refuse (needs 7+); the spill path must complete BITWISE.
        let mut store = MapStore {
            map: std::collections::HashMap::new(),
            evictions: 0,
        };
        let (ck, stats) =
            checkpointed_adjoint_spilling(&X0, STEPS, 2, &mut store, &f, &r, (1.0, 0.0));
        assert_eq!(full.0.to_bits(), ck.0.to_bits(), "xbar bitwise via spill");
        assert_eq!(full.1.to_bits(), ck.1.to_bits(), "tbar bitwise via spill");
        assert!(
            stats.peak_ram_snapshots <= 2,
            "RAM peak {} exceeded budget 2",
            stats.peak_ram_snapshots
        );
        assert!(stats.spills > 0, "tiny RAM must actually spill");
        assert!(
            stats.restores >= stats.spills,
            "every spill is read back at least once"
        );
        assert_eq!(
            store.evictions, stats.spills,
            "every spilled key is evicted when dead"
        );
        assert!(store.map.is_empty(), "store drains by the end of the sweep");
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"revolve-spill\",\"verdict\":\"pass\",\"detail\":\"bitwise == full storage at RAM budget 2 ({} spills, {} restores, {} fwd re-evals); store drained\"}}",
            stats.spills, stats.restores, stats.revolve.forward_steps
        );
    }

    #[test]
    fn zero_and_one_step_chains() {
        let f = fwd(THETA);
        let r = rev(THETA);
        let (bar, stats) = checkpointed_adjoint(&X0, 0, 0, &f, &r, (2.0, 3.0));
        assert_eq!(bar, (2.0, 3.0), "empty chain returns the seed");
        assert_eq!(stats.forward_steps, 0);
        let (bar1, _) = checkpointed_adjoint(&X0, 1, 1, &f, &r, (1.0, 0.0));
        let full1 = full_adjoint(&X0, 1, &f, &r, (1.0, 0.0));
        assert_eq!(bar1.0.to_bits(), full1.0.to_bits());
    }
}
