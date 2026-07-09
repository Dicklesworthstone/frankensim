//! The budget ALLOCATOR (bead gp3.9, plan §11.4 Bet 12): budget
//! allocation as an optimization problem HELM solves about itself —
//! "drag to 2% in 2 hours" becomes: choose a setting per knob (mesh,
//! order, tolerance, surrogate, samples, …) minimizing modeled error
//! subject to a WALL-CLOCK budget, where wall-clock composes
//! TROPICALLY (max across parallel tracks of the sum within each
//! track — §14.3: an upgrade off the critical path spends only slack,
//! and the planner sees that).
//!
//! V1 (this module, always on): GREEDY-PLUS-LOOKUP — start at the
//! cheapest feasible plan, repeatedly take the Pareto-ladder upgrade
//! with the best marginal utility Δerror/Δwall (slack upgrades are
//! free and taken first), re-planned ONLINE as a-posteriori estimates
//! replace a-priori model values. Always yields a plan or a STRUCTURED
//! [`BudgetInfeasible`] with ranked relaxations — never a shrug.
//!
//! V2 (feature `moonshot-planner`, [M]): the co-optimizer in
//! `moonshot.rs` — exact per-track multiple-choice-knapsack DP plus
//! water-filling/CMA-ES on rate-based continuous relaxations. Ships
//! OFF: the promotion gate (beats hand allocation on all three
//! flagships) belongs to the huq.15 Gauntlet.

use std::fmt::Write as _;

/// One discrete setting of a knob: the modeled error contribution and
/// the modeled cost (seconds) at this setting.
#[derive(Debug, Clone)]
pub struct KnobSetting {
    /// Human-readable label ("m=32", "tol=1e-8", …).
    pub label: String,
    /// Modeled error contribution (additive across knobs — the Error
    /// Ledger's attribution shape).
    pub error: f64,
    /// Modeled cost in seconds (additive within a track).
    pub cost: f64,
}

/// A knob: a named discrete dial on one execution track.
#[derive(Debug, Clone)]
pub struct Knob {
    /// Name ("mesh", "order", …).
    pub name: String,
    /// Execution track (wall-clock = max over tracks of within-track
    /// sums — the tropical composition).
    pub track: usize,
    /// Settings; the constructor prunes to the Pareto ladder (cost
    /// strictly increasing, error strictly decreasing).
    pub settings: Vec<KnobSetting>,
}

impl Knob {
    /// Build a knob, pruning dominated settings (anything with both
    /// higher cost and no error improvement) and sorting by cost.
    ///
    /// # Panics
    /// If no settings remain (an empty knob has no plan).
    #[must_use]
    pub fn new(name: &str, track: usize, mut settings: Vec<KnobSetting>) -> Knob {
        settings.sort_by(|a, b| a.cost.total_cmp(&b.cost).then(a.error.total_cmp(&b.error)));
        let mut ladder: Vec<KnobSetting> = Vec::new();
        for s in settings {
            if let Some(last) = ladder.last()
                && s.error >= last.error
            {
                continue; // dominated: costs more, no better
            }
            ladder.push(s);
        }
        assert!(!ladder.is_empty(), "knob {name} has no settings");
        Knob {
            name: name.to_owned(),
            track,
            settings: ladder,
        }
    }
}

/// The allocation problem: knobs, a wall-clock budget, an error target.
#[derive(Debug, Clone)]
pub struct AllocProblem {
    /// The dials.
    pub knobs: Vec<Knob>,
    /// Wall-clock budget in seconds (max over tracks).
    pub budget_s: f64,
    /// Total-error target (sum over knob contributions).
    pub error_target: f64,
}

/// A concrete plan: one setting index per knob, with its modeled
/// totals and the greedy rationale (per-upgrade marginal utilities).
#[derive(Debug, Clone)]
pub struct Plan {
    /// Chosen setting index per knob.
    pub choice: Vec<usize>,
    /// Modeled total error.
    pub total_error: f64,
    /// Modeled wall-clock (tropical: max over track sums).
    pub wall_clock: f64,
    /// Why: one line per accepted upgrade (knob, from→to, Δe, Δwall).
    pub rationale: Vec<String>,
}

/// Structured infeasibility: what was asked, what is possible, and the
/// RANKED relaxations that would make the problem solvable.
#[derive(Debug, Clone)]
pub struct BudgetInfeasible {
    /// Best achievable error within the budget.
    pub best_error_in_budget: f64,
    /// Wall-clock needed to reach the error target (ignoring budget).
    pub budget_needed_for_target: f64,
    /// Relaxations, ranked by relative distance from the request.
    pub relaxations: Vec<String>,
}

/// Wall-clock of a choice vector: tropical max over track sums.
fn wall_clock(knobs: &[Knob], choice: &[usize]) -> f64 {
    let ntracks = knobs.iter().map(|k| k.track).max().unwrap_or(0) + 1;
    let mut track_sum = vec![0.0f64; ntracks];
    for (k, &c) in knobs.iter().zip(choice) {
        track_sum[k.track] += k.settings[c].cost;
    }
    track_sum.iter().fold(0.0f64, |m, &v| m.max(v))
}

fn total_error(knobs: &[Knob], choice: &[usize]) -> f64 {
    knobs
        .iter()
        .zip(choice)
        .map(|(k, &c)| k.settings[c].error)
        .sum()
}

/// The greedy-plus-lookup allocator (V1). Deterministic: ties break by
/// lowest knob index.
///
/// # Errors
/// [`BudgetInfeasible`] when no plan meets both the budget and the
/// error target — with the achievable frontier and ranked relaxations.
pub fn allocate(problem: &AllocProblem) -> Result<Plan, BudgetInfeasible> {
    let knobs = &problem.knobs;
    let mut choice = vec![0usize; knobs.len()];
    let mut rationale = Vec::new();
    // Greedy ascent on the Pareto ladders.
    loop {
        let err = total_error(knobs, &choice);
        if err <= problem.error_target {
            return Ok(Plan {
                total_error: err,
                wall_clock: wall_clock(knobs, &choice),
                choice,
                rationale,
            });
        }
        let wall = wall_clock(knobs, &choice);
        // Candidate upgrades: next ladder step per knob that fits.
        let mut best: Option<(f64, usize, f64, f64)> = None; // (utility, knob, de, dw)
        for (ki, k) in knobs.iter().enumerate() {
            let c = choice[ki];
            if c + 1 >= k.settings.len() {
                continue;
            }
            let mut trial = choice.clone();
            trial[ki] = c + 1;
            let new_wall = wall_clock(knobs, &trial);
            if new_wall > problem.budget_s {
                continue;
            }
            let de = k.settings[c].error - k.settings[c + 1].error;
            let dw = new_wall - wall;
            // Slack upgrades (off the critical path) are free:
            // infinite utility, taken immediately in index order.
            let utility = if dw <= 0.0 { f64::INFINITY } else { de / dw };
            let better = match best {
                None => true,
                Some((u, _, _, _)) => utility > u,
            };
            if better {
                best = Some((utility, ki, de, dw));
            }
        }
        let Some((utility, ki, de, dw)) = best else {
            // Nothing fits: build the structured infeasibility.
            return Err(infeasible(problem, &choice));
        };
        let c = choice[ki];
        let mut line = String::new();
        let _ = write!(
            line,
            "{}: {} -> {} (dE {de:.3e}, dWall {dw:.3e}s, utility {utility:.3e})",
            knobs[ki].name,
            knobs[ki].settings[c].label,
            knobs[ki].settings[c + 1].label
        );
        rationale.push(line);
        choice[ki] = c + 1;
    }
}

/// Build the structured infeasibility record from the stalled state.
fn infeasible(problem: &AllocProblem, stalled: &[usize]) -> BudgetInfeasible {
    let knobs = &problem.knobs;
    // Best error within budget: greedy already maximized error
    // reduction under the budget; report where it stalled.
    let best_error_in_budget = total_error(knobs, stalled);
    // Budget needed for the target: continue the greedy ladder
    // ignoring the budget (max settings are the floor of achievable
    // error; if even that misses the target, say so).
    let mut choice = stalled.to_vec();
    loop {
        let err = total_error(knobs, &choice);
        if err <= problem.error_target {
            break;
        }
        let mut best: Option<(f64, usize)> = None;
        for (ki, k) in knobs.iter().enumerate() {
            let c = choice[ki];
            if c + 1 >= k.settings.len() {
                continue;
            }
            let de = k.settings[c].error - k.settings[c + 1].error;
            let dc = k.settings[c + 1].cost - k.settings[c].cost;
            let u = if dc <= 0.0 { f64::INFINITY } else { de / dc };
            if best.is_none_or(|(bu, _)| u > bu) {
                best = Some((u, ki));
            }
        }
        let Some((_, ki)) = best else { break };
        choice[ki] += 1;
    }
    let budget_needed = wall_clock(knobs, &choice);
    let floor_error = total_error(knobs, &choice);
    let mut relaxations = Vec::new();
    let mut ranked: Vec<(f64, String)> = Vec::new();
    if floor_error <= problem.error_target {
        let rel = (budget_needed - problem.budget_s) / problem.budget_s.max(1e-30);
        ranked.push((
            rel,
            format!(
                "raise budget to {budget_needed:.3e}s (+{:.1}%) to reach the error target",
                rel * 100.0
            ),
        ));
    } else {
        ranked.push((
            f64::INFINITY,
            format!(
                "error target {:.3e} is below the model floor {floor_error:.3e} even at max settings — target unreachable at ANY budget",
                problem.error_target
            ),
        ));
    }
    let rel_t = (best_error_in_budget - problem.error_target) / problem.error_target.max(1e-30);
    ranked.push((
        rel_t,
        format!(
            "relax error target to {best_error_in_budget:.3e} (+{:.1}%) to fit the current budget",
            rel_t * 100.0
        ),
    ));
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
    for (_, r) in ranked {
        relaxations.push(r);
    }
    BudgetInfeasible {
        best_error_in_budget,
        budget_needed_for_target: budget_needed,
        relaxations,
    }
}

/// The ONLINE re-planner: holds the problem plus a-posteriori error
/// estimates that override the a-priori model values (the DWR estimate
/// arrives → the plan updates).
#[derive(Debug, Clone)]
pub struct Allocator {
    problem: AllocProblem,
}

impl Allocator {
    /// Wrap a problem.
    #[must_use]
    pub fn new(problem: AllocProblem) -> Allocator {
        Allocator { problem }
    }

    /// The current problem view (with overrides applied).
    #[must_use]
    pub fn problem(&self) -> &AllocProblem {
        &self.problem
    }

    /// Replace the modeled error of one knob setting with a measured
    /// a-posteriori estimate (e.g. the DWR estimate for the current
    /// mesh level).
    ///
    /// # Panics
    /// On out-of-range indices (programmer error).
    pub fn observe_error(&mut self, knob: usize, setting: usize, measured: f64) {
        self.problem.knobs[knob].settings[setting].error = measured;
    }

    /// Re-plan from the current estimates.
    ///
    /// # Errors
    /// [`BudgetInfeasible`] as in [`allocate`].
    pub fn replan(&self) -> Result<Plan, BudgetInfeasible> {
        allocate(&self.problem)
    }
}

/// Brute-force oracle (exponential; FIXTURE-SCALE ONLY): the exact
/// minimum error subject to the budget — what the batteries compare
/// the planners against.
#[must_use]
pub fn oracle_min_error(problem: &AllocProblem) -> Option<(Vec<usize>, f64)> {
    let knobs = &problem.knobs;
    let mut best: Option<(Vec<usize>, f64)> = None;
    let mut choice = vec![0usize; knobs.len()];
    loop {
        if wall_clock(knobs, &choice) <= problem.budget_s {
            let e = total_error(knobs, &choice);
            if best.as_ref().is_none_or(|(_, be)| e < *be) {
                best = Some((choice.clone(), e));
            }
        }
        // Odometer increment.
        let mut i = 0;
        loop {
            if i == knobs.len() {
                return best;
            }
            choice[i] += 1;
            if choice[i] < knobs[i].settings.len() {
                break;
            }
            choice[i] = 0;
            i += 1;
        }
    }
}

/// Public wall-clock/error evaluators (the batteries and the moonshot
/// co-optimizer share them).
#[must_use]
pub fn plan_wall_clock(knobs: &[Knob], choice: &[usize]) -> f64 {
    wall_clock(knobs, choice)
}

/// Total modeled error of a choice.
#[must_use]
pub fn plan_total_error(knobs: &[Knob], choice: &[usize]) -> f64 {
    total_error(knobs, choice)
}
