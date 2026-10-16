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
//! [`AllocationError`] with typed input refusals or ranked budget
//! relaxations — never a shrug.
//!
//! V2 (feature `moonshot-planner`, [M]): the co-optimizer in
//! `moonshot.rs` — exact per-track multiple-choice-knapsack DP plus
//! water-filling/CMA-ES on rate-based continuous relaxations. Ships
//! OFF: the promotion gate (beats hand allocation on all three
//! flagships) belongs to the huq.15 Gauntlet.

use std::fmt::Write as _;

/// Maximum number of logical execution tracks accepted by the V1
/// allocator. Track IDs index planner metadata, so they are bounded
/// independently of the amount of host memory available.
pub const MAX_EXECUTION_TRACKS: usize = 1024;

/// Maximum number of knobs accepted by one allocation problem.
pub const MAX_ALLOCATION_KNOBS: usize = 256;

/// Maximum settings retained on one knob before Pareto pruning.
pub const MAX_SETTINGS_PER_KNOB: usize = 256;

/// Maximum aggregate setting rows accepted across all knobs.
pub const MAX_TOTAL_SETTINGS: usize = 4096;

/// Maximum Cartesian choices visited by the fixture-only exact oracle.
pub const MAX_ORACLE_COMBINATIONS: usize = 1_000_000;

/// A malformed allocation problem or online observation.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanInputError {
    /// A problem contains too many independent knobs.
    TooManyKnobs {
        /// Supplied knob count.
        count: usize,
        /// Inclusive limit.
        max: usize,
    },
    /// A track ID falls outside the planner's bounded domain.
    TrackOutOfRange {
        /// Knob carrying the invalid track.
        knob: String,
        /// Supplied track ID.
        track: usize,
        /// Exclusive upper bound.
        max_exclusive: usize,
    },
    /// A knob has no selectable settings.
    EmptySettings {
        /// Knob with the empty setting list.
        knob: String,
    },
    /// One knob contains too many setting rows.
    TooManySettings {
        /// Knob carrying the oversized ladder.
        knob: String,
        /// Supplied setting count.
        count: usize,
        /// Inclusive limit.
        max: usize,
    },
    /// The problem contains too many setting rows in aggregate.
    TooManyTotalSettings {
        /// Supplied aggregate setting count.
        count: usize,
        /// Inclusive limit.
        max: usize,
    },
    /// A setting contains a nonfinite or negative scalar.
    InvalidSettingValue {
        /// Knob carrying the setting.
        knob: String,
        /// Setting index.
        setting: usize,
        /// Invalid field (`error`, `cost`, or `observed error`).
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// A problem-level scalar is nonfinite or negative.
    InvalidProblemValue {
        /// Invalid field (`budget_s` or `error_target`).
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// Public struct construction bypassed the Pareto-ladder constructor.
    InvalidParetoLadder {
        /// Malformed knob.
        knob: String,
        /// First index in the invalid adjacent pair.
        lower: usize,
        /// Second index in the invalid adjacent pair.
        upper: usize,
    },
    /// Valid individual scalars overflow when combined into a plan total.
    AggregateOverflow {
        /// Aggregate that overflowed.
        field: &'static str,
    },
    /// An online observation named a missing knob.
    KnobIndexOutOfRange {
        /// Supplied knob index.
        knob: usize,
        /// Number of knobs in the problem.
        knob_count: usize,
    },
    /// An online observation named a missing setting.
    SettingIndexOutOfRange {
        /// Knob index.
        knob: usize,
        /// Supplied setting index.
        setting: usize,
        /// Number of settings on the knob.
        setting_count: usize,
    },
    /// A choice vector does not name exactly one setting per knob.
    ChoiceLengthMismatch {
        /// Number of knobs.
        knob_count: usize,
        /// Number of supplied choices.
        choice_count: usize,
    },
    /// A choice vector names a missing setting.
    ChoiceIndexOutOfRange {
        /// Knob index.
        knob: usize,
        /// Supplied setting index.
        choice: usize,
        /// Number of settings on the knob.
        setting_count: usize,
    },
    /// The exact fixture oracle would exceed its deterministic work cap.
    OracleWorkLimitExceeded {
        /// Maximum combinations the fixture oracle may enumerate.
        max_combinations: usize,
    },
}

impl core::fmt::Display for PlanInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooManyKnobs { count, max } => {
                write!(f, "allocation has {count} knobs; the limit is {max}")
            }
            Self::TrackOutOfRange {
                knob,
                track,
                max_exclusive,
            } => write!(
                f,
                "knob {knob:?} uses track {track}; tracks must be below {max_exclusive}"
            ),
            Self::EmptySettings { knob } => {
                write!(f, "knob {knob:?} has no settings")
            }
            Self::TooManySettings { knob, count, max } => write!(
                f,
                "knob {knob:?} has {count} settings; the per-knob limit is {max}"
            ),
            Self::TooManyTotalSettings { count, max } => write!(
                f,
                "allocation has {count} total settings; the aggregate limit is {max}"
            ),
            Self::InvalidSettingValue {
                knob,
                setting,
                field,
                value,
            } => write!(
                f,
                "knob {knob:?} setting {setting} has invalid {field} {value:?}; values must be finite and non-negative"
            ),
            Self::InvalidProblemValue { field, value } => write!(
                f,
                "allocation {field} is {value:?}; values must be finite and non-negative"
            ),
            Self::InvalidParetoLadder { knob, lower, upper } => write!(
                f,
                "knob {knob:?} settings {lower} and {upper} are not a strict Pareto ladder (cost must increase while error decreases)"
            ),
            Self::AggregateOverflow { field } => {
                write!(f, "allocation {field} overflows finite f64 range")
            }
            Self::KnobIndexOutOfRange { knob, knob_count } => write!(
                f,
                "online observation names knob {knob}, but the problem has {knob_count} knob(s)"
            ),
            Self::SettingIndexOutOfRange {
                knob,
                setting,
                setting_count,
            } => write!(
                f,
                "online observation names setting {setting} on knob {knob}, but it has {setting_count} setting(s)"
            ),
            Self::ChoiceLengthMismatch {
                knob_count,
                choice_count,
            } => write!(
                f,
                "choice vector has {choice_count} entries for {knob_count} knob(s)"
            ),
            Self::ChoiceIndexOutOfRange {
                knob,
                choice,
                setting_count,
            } => write!(
                f,
                "choice {choice} on knob {knob} is out of range for {setting_count} setting(s)"
            ),
            Self::OracleWorkLimitExceeded { max_combinations } => write!(
                f,
                "fixture oracle would exceed its {max_combinations} combination work limit"
            ),
        }
    }
}

impl std::error::Error for PlanInputError {}

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
    /// # Errors
    /// [`PlanInputError`] when the track is out of range, the setting
    /// list is empty/oversized, or a modeled error/cost is nonfinite or
    /// negative.
    pub fn new(
        name: &str,
        track: usize,
        settings: Vec<KnobSetting>,
    ) -> Result<Knob, PlanInputError> {
        validate_knob_parts(name, track, &settings)?;
        Ok(Knob {
            name: name.to_owned(),
            track,
            settings: pareto_ladder(settings),
        })
    }
}

fn pareto_ladder(mut settings: Vec<KnobSetting>) -> Vec<KnobSetting> {
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
    ladder
}

fn validate_knob_parts(
    name: &str,
    track: usize,
    settings: &[KnobSetting],
) -> Result<(), PlanInputError> {
    if track >= MAX_EXECUTION_TRACKS {
        return Err(PlanInputError::TrackOutOfRange {
            knob: name.to_owned(),
            track,
            max_exclusive: MAX_EXECUTION_TRACKS,
        });
    }
    if settings.is_empty() {
        return Err(PlanInputError::EmptySettings {
            knob: name.to_owned(),
        });
    }
    if settings.len() > MAX_SETTINGS_PER_KNOB {
        return Err(PlanInputError::TooManySettings {
            knob: name.to_owned(),
            count: settings.len(),
            max: MAX_SETTINGS_PER_KNOB,
        });
    }
    for (setting, value) in settings.iter().enumerate() {
        for (field, scalar) in [("error", value.error), ("cost", value.cost)] {
            if !scalar.is_finite() || scalar < 0.0 {
                return Err(PlanInputError::InvalidSettingValue {
                    knob: name.to_owned(),
                    setting,
                    field,
                    value: scalar,
                });
            }
        }
    }
    Ok(())
}

/// The allocation problem: knobs, a wall-clock budget, an error target.
/// An empty knob list is the intentional additive/tropical identity:
/// its sole plan has zero error and zero wall cost.
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
    /// Modeled error at the deterministic greedy allocator's stalled
    /// in-budget plan. This is a feasible heuristic result, not the exact
    /// optimum or a lower bound on achievable error.
    pub best_error_in_budget: f64,
    /// Wall-clock sufficient for the same greedy policy to reach the target
    /// when the target is reachable. This is not a minimum-budget claim; when
    /// the target is below the model floor, it is the maximum-setting wall
    /// clock used to diagnose that floor.
    pub budget_needed_for_target: f64,
    /// Relaxations, ranked by relative distance from the request.
    pub relaxations: Vec<String>,
}

/// A typed allocator refusal.
#[derive(Debug, Clone)]
pub enum AllocationError {
    /// The request is malformed and cannot be reasoned about safely.
    InvalidInput(PlanInputError),
    /// Even the cheapest setting on every knob exceeds the wall budget.
    MinimumPlanExceedsBudget {
        /// Requested wall budget.
        budget_s: f64,
        /// Wall time of the cheapest complete plan.
        minimum_wall_s: f64,
    },
    /// At least one complete plan fits, but none meets the error target.
    BudgetInfeasible(BudgetInfeasible),
}

impl core::fmt::Display for AllocationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(error) => write!(f, "invalid allocation input: {error}"),
            Self::MinimumPlanExceedsBudget {
                budget_s,
                minimum_wall_s,
            } => write!(
                f,
                "minimum complete plan requires {minimum_wall_s:.3e}s, exceeding the {budget_s:.3e}s wall budget"
            ),
            Self::BudgetInfeasible(error) => write!(
                f,
                "greedy allocator stalled at in-budget error {}; target-budget diagnostic is {}s (see ranked relaxations)",
                error.best_error_in_budget, error.budget_needed_for_target
            ),
        }
    }
}

impl std::error::Error for AllocationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidInput(error) => Some(error),
            Self::MinimumPlanExceedsBudget { .. } | Self::BudgetInfeasible(_) => None,
        }
    }
}

impl From<PlanInputError> for AllocationError {
    fn from(error: PlanInputError) -> Self {
        Self::InvalidInput(error)
    }
}

fn checked_add(current: f64, increment: f64, field: &'static str) -> Result<f64, PlanInputError> {
    let total = current + increment;
    if total.is_finite() {
        Ok(total)
    } else {
        Err(PlanInputError::AggregateOverflow { field })
    }
}

fn validate_knobs(knobs: &[Knob]) -> Result<(), PlanInputError> {
    if knobs.len() > MAX_ALLOCATION_KNOBS {
        return Err(PlanInputError::TooManyKnobs {
            count: knobs.len(),
            max: MAX_ALLOCATION_KNOBS,
        });
    }

    let mut total_settings = 0usize;
    for knob in knobs {
        validate_knob_parts(&knob.name, knob.track, &knob.settings)?;
        total_settings += knob.settings.len();
        if total_settings > MAX_TOTAL_SETTINGS {
            return Err(PlanInputError::TooManyTotalSettings {
                count: total_settings,
                max: MAX_TOTAL_SETTINGS,
            });
        }
        for (lower, pair) in knob.settings.windows(2).enumerate() {
            if pair[0].cost >= pair[1].cost || pair[0].error <= pair[1].error {
                return Err(PlanInputError::InvalidParetoLadder {
                    knob: knob.name.clone(),
                    lower,
                    upper: lower + 1,
                });
            }
        }
    }

    // Bound all possible choice totals, not just the cheapest plan, so
    // subsequent utility arithmetic cannot silently become infinite.
    let mut maximum_track_cost = [0.0f64; MAX_EXECUTION_TRACKS];
    let mut maximum_total_error = 0.0f64;
    for knob in knobs {
        let maximum_cost = knob.settings.last().expect("validated non-empty").cost;
        maximum_track_cost[knob.track] =
            checked_add(maximum_track_cost[knob.track], maximum_cost, "track cost")?;
        let maximum_error = knob.settings[0].error;
        maximum_total_error = checked_add(maximum_total_error, maximum_error, "total error")?;
    }
    Ok(())
}

/// Validate every invariant assumed by allocation, including public
/// struct values that bypassed [`Knob::new`]. Validation completes
/// before any work is sized by a caller-controlled track ID.
pub(crate) fn validate_problem(problem: &AllocProblem) -> Result<(), PlanInputError> {
    for (field, value) in [
        ("budget_s", problem.budget_s),
        ("error_target", problem.error_target),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(PlanInputError::InvalidProblemValue { field, value });
        }
    }
    validate_knobs(&problem.knobs)
}

fn validate_choice(knobs: &[Knob], choice: &[usize]) -> Result<(), PlanInputError> {
    if choice.len() != knobs.len() {
        return Err(PlanInputError::ChoiceLengthMismatch {
            knob_count: knobs.len(),
            choice_count: choice.len(),
        });
    }
    for (knob, (&selected, settings)) in choice
        .iter()
        .zip(knobs.iter().map(|knob| &knob.settings))
        .enumerate()
    {
        if selected >= settings.len() {
            return Err(PlanInputError::ChoiceIndexOutOfRange {
                knob,
                choice: selected,
                setting_count: settings.len(),
            });
        }
    }
    Ok(())
}

/// Wall-clock of a choice vector: tropical max over track sums.
fn wall_clock(knobs: &[Knob], choice: &[usize]) -> f64 {
    // Store only tracks that actually occur. A sparse or attacker-sized
    // numeric ID can therefore never determine an allocation size.
    let mut track_sum: Vec<(usize, f64)> = Vec::with_capacity(knobs.len());
    for (k, &c) in knobs.iter().zip(choice) {
        if let Some((_, sum)) = track_sum.iter_mut().find(|(track, _)| *track == k.track) {
            *sum += k.settings[c].cost;
        } else {
            track_sum.push((k.track, k.settings[c].cost));
        }
    }
    track_sum
        .iter()
        .fold(0.0f64, |maximum, (_, sum)| maximum.max(*sum))
}

fn total_error(knobs: &[Knob], choice: &[usize]) -> f64 {
    knobs
        .iter()
        .zip(choice)
        .map(|(k, &c)| k.settings[c].error)
        .fold(0.0f64, |total, error| total + error)
}

/// Select the next upgrade under the V1 marginal-utility rule. `None` as the
/// budget means that every finite upgrade is admitted; a finite budget filters
/// out trial choices whose tropical wall clock would exceed it. Iteration order
/// preserves the public lowest-knob-index tie break in both modes.
fn best_greedy_upgrade(
    knobs: &[Knob],
    choice: &[usize],
    budget_s: Option<f64>,
) -> Option<(f64, usize, f64, f64)> {
    let wall = wall_clock(knobs, choice);
    let mut best: Option<(f64, usize, f64, f64)> = None; // (utility, knob, de, dw)
    for (ki, knob) in knobs.iter().enumerate() {
        let current = choice[ki];
        if current + 1 >= knob.settings.len() {
            continue;
        }
        let mut trial = choice.to_vec();
        trial[ki] = current + 1;
        let new_wall = wall_clock(knobs, &trial);
        if budget_s.is_some_and(|budget_s| new_wall > budget_s) {
            continue;
        }
        let de = knob.settings[current].error - knob.settings[current + 1].error;
        let dw = new_wall - wall;
        // Slack upgrades (off the critical path) are free. Equal utilities
        // retain the earlier knob because replacement is strictly greater.
        let utility = if dw <= 0.0 { f64::INFINITY } else { de / dw };
        if best.is_none_or(|(best_utility, _, _, _)| utility > best_utility) {
            best = Some((utility, ki, de, dw));
        }
    }
    best
}

/// The greedy-plus-lookup allocator (V1). Deterministic: ties break by
/// lowest knob index.
///
/// # Errors
/// [`AllocationError::InvalidInput`] for malformed scalar/track/ladder
/// inputs; [`AllocationError::MinimumPlanExceedsBudget`] when no
/// complete plan fits at all; [`AllocationError::BudgetInfeasible`]
/// when a plan fits but cannot meet the error target.
pub fn allocate(problem: &AllocProblem) -> Result<Plan, AllocationError> {
    validate_problem(problem)?;
    let knobs = &problem.knobs;
    let mut choice = vec![0usize; knobs.len()];
    let mut rationale = Vec::new();
    let minimum_wall_s = wall_clock(knobs, &choice);
    if minimum_wall_s > problem.budget_s {
        return Err(AllocationError::MinimumPlanExceedsBudget {
            budget_s: problem.budget_s,
            minimum_wall_s,
        });
    }
    // Greedy ascent on the Pareto ladders.
    loop {
        let err = total_error(knobs, &choice);
        if err <= problem.error_target {
            let wall_clock = wall_clock(knobs, &choice);
            debug_assert!(wall_clock <= problem.budget_s);
            return Ok(Plan {
                total_error: err,
                wall_clock,
                choice,
                rationale,
            });
        }
        let Some((utility, ki, de, dw)) =
            best_greedy_upgrade(knobs, &choice, Some(problem.budget_s))
        else {
            // Nothing fits: build the structured infeasibility.
            return Err(AllocationError::BudgetInfeasible(infeasible(
                problem, &choice,
            )));
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
    // Report the actual greedy stalled point. V1 is deliberately heuristic,
    // so this is neither an exact optimum nor a lower bound.
    let best_error_in_budget = total_error(knobs, stalled);
    // Restart at the minimum plan and replay the exact allocation rule without
    // a cap. The first target-reaching prefix has nondecreasing wall clock, so
    // every one of its selected upgrades fits when allocation is rerun at the
    // final prefix wall. Because the uncapped winner already has maximum
    // marginal utility (with the same index tie break), the capped rerun follows
    // the identical prefix. The result is therefore sufficient for this greedy
    // policy, although it need not be the globally minimum feasible budget.
    let mut choice = vec![0usize; knobs.len()];
    loop {
        let err = total_error(knobs, &choice);
        if err <= problem.error_target {
            break;
        }
        let Some((_, ki, _, _)) = best_greedy_upgrade(knobs, &choice, None) else {
            break;
        };
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
                "raise budget to {budget_needed}s (+{:.1}%) to reach the error target",
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
            "relax error target to {best_error_in_budget} (+{:.1}%) to accept the current greedy plan",
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
    /// Wrap a validated problem.
    ///
    /// # Errors
    /// [`PlanInputError`] when the supplied problem is malformed.
    pub fn new(problem: AllocProblem) -> Result<Allocator, PlanInputError> {
        validate_problem(&problem)?;
        Ok(Allocator { problem })
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
    /// The update is transactional: invalid indices/values leave the
    /// current problem untouched. A valid update re-prunes the knob's
    /// Pareto ladder before it becomes visible to [`replan`](Self::replan).
    ///
    /// # Errors
    /// [`PlanInputError`] for out-of-range indices, invalid measurements,
    /// or a resulting aggregate outside finite `f64` range.
    pub fn observe_error(
        &mut self,
        knob: usize,
        setting: usize,
        measured: f64,
    ) -> Result<(), PlanInputError> {
        let Some(current_knob) = self.problem.knobs.get(knob) else {
            return Err(PlanInputError::KnobIndexOutOfRange {
                knob,
                knob_count: self.problem.knobs.len(),
            });
        };
        if setting >= current_knob.settings.len() {
            return Err(PlanInputError::SettingIndexOutOfRange {
                knob,
                setting,
                setting_count: current_knob.settings.len(),
            });
        }
        if !measured.is_finite() || measured < 0.0 {
            return Err(PlanInputError::InvalidSettingValue {
                knob: current_knob.name.clone(),
                setting,
                field: "observed error",
                value: measured,
            });
        }

        let mut candidate = self.problem.clone();
        candidate.knobs[knob].settings[setting].error = measured;
        let settings = std::mem::take(&mut candidate.knobs[knob].settings);
        candidate.knobs[knob].settings = pareto_ladder(settings);
        validate_problem(&candidate)?;
        self.problem = candidate;
        Ok(())
    }

    /// Re-plan from the current estimates.
    ///
    /// # Errors
    /// [`AllocationError`] as in [`allocate`].
    pub fn replan(&self) -> Result<Plan, AllocationError> {
        allocate(&self.problem)
    }
}

/// Brute-force oracle (exponential; FIXTURE-SCALE ONLY): the exact
/// minimum error subject to the budget — what the batteries compare
/// the planners against.
///
/// # Errors
/// [`PlanInputError`] for a malformed problem or when the Cartesian
/// choice space exceeds [`MAX_ORACLE_COMBINATIONS`].
pub fn oracle_min_error(
    problem: &AllocProblem,
) -> Result<Option<(Vec<usize>, f64)>, PlanInputError> {
    validate_problem(problem)?;
    let knobs = &problem.knobs;
    if knobs.is_empty() {
        return Ok(Some((Vec::new(), 0.0)));
    }
    let mut combinations = 1usize;
    for knob in knobs {
        if combinations > MAX_ORACLE_COMBINATIONS / knob.settings.len() {
            return Err(PlanInputError::OracleWorkLimitExceeded {
                max_combinations: MAX_ORACLE_COMBINATIONS,
            });
        }
        combinations *= knob.settings.len();
    }
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
                return Ok(best);
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
///
/// # Errors
/// [`PlanInputError`] when the knob domain or choice vector is invalid.
pub fn plan_wall_clock(knobs: &[Knob], choice: &[usize]) -> Result<f64, PlanInputError> {
    validate_knobs(knobs)?;
    validate_choice(knobs, choice)?;
    Ok(wall_clock(knobs, choice))
}

/// Total modeled error of a choice.
///
/// # Errors
/// [`PlanInputError`] when the knob domain or choice vector is invalid.
pub fn plan_total_error(knobs: &[Knob], choice: &[usize]) -> Result<f64, PlanInputError> {
    validate_knobs(knobs)?;
    validate_choice(knobs, choice)?;
    Ok(total_error(knobs, choice))
}
