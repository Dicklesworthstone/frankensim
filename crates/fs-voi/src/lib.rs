//! fs-voi — value-of-information and active validation. Layer: L4.
//!
//! The strategic layer deciding WHAT INFORMATION TO ACQUIRE NEXT. The whole
//! point is the ESTIMATOR-vs-DECISION distinction: spend compute where it can
//! CHANGE A DECISION (flip which design ranks first), not merely where an
//! estimator is uncertain. A design that is wildly uncertain but clearly last
//! is not worth another cent.
//!
//! Each design's objective is an Evidence-carrying estimate: a mean plus an
//! [`Uncertainty`] with NUMERICAL, STATISTICAL, and MODEL-form components. The
//! decision is "which design minimizes the objective". Over that:
//! - [`ranking_flip_probability`] — `P(a runner-up is actually better)`;
//! - [`expected_opportunity_loss`] — the FULL multi-alternative expected
//!   opportunity loss of the current decision (the ceiling on what any
//!   information is worth): EVERY alternative participates, so a
//!   worse-mean/large-variance third design cannot be silently dropped;
//! - [`top_two_evpi_surrogate`] — the cheap two-design closed form, kept ONLY
//!   as an action-ranking surrogate; it must never gate a global robustness
//!   claim;
//! - [`action_value`] — how much a menu action (surrogate / simulate / sample /
//!   test / refine) REDUCES the surrogate EVPI, per unit cost. An action on a
//!   decision-IRRELEVANT design is worth ~0;
//! - [`recommend`] — the best action, or STOP when the decision is already
//!   robust (the STOP gate is the FULL evaluator, never the surrogate);
//! - [`heuristic_choice`] — the uncertainty-proportional baseline VOI must beat
//!   ([M] discipline).
//!
//! Deterministic; no dependencies (Gaussian decision algebra with an in-house
//! normal CDF).

/// Semantic version of the decision algebra. Bump this whenever a change can
/// alter opportunity-loss result bits, best/runner-up selection, uncertainty
/// composition, or the Gaussian helper semantics.
// v2 (bead sj31i.5): `evpi`/`evpi_by` became the renamed
// `top_two_evpi_surrogate{,_by}` (mean ordering is not stochastic
// dominance — the two lowest means do NOT bound the decision's true
// opportunity loss); the robustness-bearing evaluator is the new full
// multi-alternative `expected_opportunity_loss{,_by}`. Uncertainty
// composition and pairwise deviation sums moved to overflow-safe
// scaled norms, so finite near-`sqrt(MAX)` components no longer
// overflow into false infinities (this moves result bits at extreme
// scales; ordinary-scale bits move only through the scaled-norm
// rounding order).
pub const EVPI_SEMANTICS_VERSION: u64 = 2;

/// Composite-Simpson panel count for one [`expected_opportunity_loss`]
/// evaluation (even, fixed): the deterministic quadrature reads every
/// design's survival function at `EOL_QUADRATURE_PANELS + 1` nodes.
pub const EOL_QUADRATURE_PANELS: usize = 256;

/// Half-width, in per-design standard deviations, of the integration
/// window. Beyond `12σ` a Gaussian tail carries `< 2e-33` probability,
/// far below the quadrature's own resolution.
const EOL_TAIL_HALF_WIDTHS: f64 = 12.0;

/// The three uncertainty components of an estimate (they compose in quadrature).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Uncertainty {
    /// Numerical / discretization uncertainty.
    pub numerical: f64,
    /// Statistical (sampling) uncertainty.
    pub statistical: f64,
    /// Model-form uncertainty.
    pub model: f64,
}

/// Which uncertainty component an action reduces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Component {
    /// Numerical / discretization.
    Numerical,
    /// Statistical / sampling.
    Statistical,
    /// Model-form.
    Model,
}

impl Uncertainty {
    /// The total standard deviation `√(num² + stat² + model²)`, computed
    /// as a SCALED norm: finite components near `√MAX` compose to their
    /// representable total instead of overflowing the naive variance sum
    /// into a false infinity (bead sj31i.5).
    #[must_use]
    pub fn total_std(&self) -> f64 {
        scaled_norm(&[self.numerical, self.statistical, self.model])
    }

    /// The dominant uncertainty component.
    #[must_use]
    pub fn dominant(&self) -> Component {
        let (n, s, m) = (
            self.numerical.abs(),
            self.statistical.abs(),
            self.model.abs(),
        );
        if m >= n && m >= s {
            Component::Model
        } else if s >= n {
            Component::Statistical
        } else {
            Component::Numerical
        }
    }

    fn reduced(&self, component: Component, factor: f64) -> Uncertainty {
        let f = (1.0 - factor).clamp(0.0, 1.0);
        let mut u = *self;
        match component {
            Component::Numerical => u.numerical *= f,
            Component::Statistical => u.statistical *= f,
            Component::Model => u.model *= f,
        }
        u
    }
}

/// A design with an Evidence-carrying objective estimate (minimizing).
#[derive(Debug, Clone, PartialEq)]
pub struct DesignEstimate {
    /// The design id.
    pub name: String,
    /// The estimated objective (lower is better).
    pub mean: f64,
    /// Its decomposed uncertainty.
    pub uncertainty: Uncertainty,
}

impl DesignEstimate {
    /// A design estimate.
    #[must_use]
    pub fn new(name: impl Into<String>, mean: f64, uncertainty: Uncertainty) -> DesignEstimate {
        DesignEstimate {
            name: name.into(),
            mean,
            uncertainty,
        }
    }
    fn std(&self) -> f64 {
        self.uncertainty.total_std()
    }
}

/// `P(other is actually better than chosen)` — the ranking-flip probability
/// under Gaussian objectives (minimizing).
#[must_use]
pub fn ranking_flip_probability(chosen: &DesignEstimate, other: &DesignEstimate) -> f64 {
    let sigma = scaled_norm(&[chosen.std(), other.std()]);
    if sigma <= 0.0 {
        return f64::from(u8::from(other.mean < chosen.mean));
    }
    // P(obj_other < obj_chosen) = Φ((mean_chosen − mean_other) / σ).
    normal_cdf((chosen.mean - other.mean) / sigma)
}

/// The current decision posture: the best design, the runner-up, and the
/// probability the ranking flips between them.
#[derive(Debug, Clone, PartialEq)]
pub struct DecisionPosture {
    /// The current best design (lowest mean).
    pub best: String,
    /// The closest competitor.
    pub runner_up: String,
    /// `P(runner-up actually better)`.
    pub flip_probability: f64,
}

fn top_two(designs: &[DesignEstimate]) -> Option<(&DesignEstimate, &DesignEstimate)> {
    let (best, runner) = top_two_indices_by(designs.len(), &|idx| designs[idx].mean)?;
    Some((&designs[best], &designs[runner]))
}

/// The allocation-free top-two scan over means-by-index — the shared
/// core `evpi`, `evpi_by`, and `decision_posture` all route through
/// (one code path, so accessor-driven callers are bitwise-identical to
/// slice-driven ones). Non-finite means are skipped; ties break toward
/// the LOWER index, so a caller holding a canonically ordered menu
/// gets canonical tie-breaking with no clone and no sort.
fn top_two_indices_by(len: usize, mean_at: &dyn Fn(usize) -> f64) -> Option<(usize, usize)> {
    if len < 2 {
        return None;
    }
    let mut best: Option<(usize, f64)> = None;
    let mut runner: Option<(usize, f64)> = None;
    for idx in 0..len {
        let mean = mean_at(idx);
        if !mean.is_finite() {
            continue;
        }
        match best {
            None => best = Some((idx, mean)),
            Some(current_best) if estimate_precedes((idx, mean), current_best) => {
                runner = best;
                best = Some((idx, mean));
            }
            Some(_) => match runner {
                None => runner = Some((idx, mean)),
                Some(current_runner) if estimate_precedes((idx, mean), current_runner) => {
                    runner = Some((idx, mean));
                }
                Some(_) => {}
            },
        }
    }
    let (best, _) = best?;
    let (runner, _) = runner?;
    Some((best, runner))
}

fn estimate_precedes((a_idx, a_mean): (usize, f64), (b_idx, b_mean): (usize, f64)) -> bool {
    match a_mean.total_cmp(&b_mean) {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Equal => a_idx < b_idx,
        std::cmp::Ordering::Greater => false,
    }
}

/// The decision posture over a set of designs.
#[must_use]
pub fn decision_posture(designs: &[DesignEstimate]) -> Option<DecisionPosture> {
    let (best, runner) = top_two(designs)?;
    Some(DecisionPosture {
        best: best.name.clone(),
        runner_up: runner.name.clone(),
        flip_probability: ranking_flip_probability(best, runner),
    })
}

/// Pairwise expected opportunity loss `E[(obj_chosen − obj_other)⁺]` for the
/// two closest designs — the expected regret of the current top-two decision.
fn pairwise_evpi_scalar(chosen_mean: f64, chosen_std: f64, other_mean: f64, other_std: f64) -> f64 {
    let sigma = scaled_norm(&[chosen_std, other_std]);
    let delta = chosen_mean - other_mean; // ≤ 0 when `chosen` is best
    if sigma <= 0.0 {
        return (-delta).max(0.0) * 0.0; // no uncertainty → no opportunity loss
    }
    let z = delta / sigma;
    sigma * normal_pdf(z) + delta * normal_cdf(z)
}

/// The TOP-TWO closed-form opportunity-loss SURROGATE (bead sj31i.5:
/// the renamed former `evpi`). It sees only the two lowest posterior
/// means; a third alternative with a worse mean but a large uncertainty
/// can hold material probability of being optimal while this value sits
/// near zero — mean ordering is not stochastic dominance. It exists as
/// a cheap ACTION-RANKING surrogate and MUST NOT gate a global
/// robustness claim; robustness belongs to
/// [`expected_opportunity_loss`], which includes every alternative.
#[must_use]
pub fn top_two_evpi_surrogate(designs: &[DesignEstimate]) -> f64 {
    top_two_evpi_surrogate_by(designs.len(), &|idx| designs[idx].mean, &|idx| {
        designs[idx].std()
    })
}

/// Allocation-free [`top_two_evpi_surrogate`] over indexed accessors —
/// one shared code path, so results are bitwise-identical to the slice
/// form. `std_at` is consulted only for the final top-two pair. Callers
/// own their index order: ties break toward the LOWER index. The same
/// surrogate caveat applies: never a global robustness gate.
#[must_use]
pub fn top_two_evpi_surrogate_by(
    len: usize,
    mean_at: &dyn Fn(usize) -> f64,
    std_at: &dyn Fn(usize) -> f64,
) -> f64 {
    match top_two_indices_by(len, mean_at) {
        Some((best, runner)) => {
            pairwise_evpi_scalar(mean_at(best), std_at(best), mean_at(runner), std_at(runner))
        }
        None => 0.0,
    }
}

/// The FULL expected value of perfect information for the current
/// decision: `E[obj_best − min_j obj_j]` with EVERY alternative
/// participating (independent Gaussian objectives, minimizing). Near
/// zero means the decision is robust against the WHOLE menu — this is
/// the only evaluator in this crate allowed to support a global
/// robustness claim.
#[must_use]
pub fn expected_opportunity_loss(designs: &[DesignEstimate]) -> f64 {
    expected_opportunity_loss_by(designs.len(), &|idx| designs[idx].mean, &|idx| {
        designs[idx].std()
    })
}

/// Allocation-free [`expected_opportunity_loss`] over indexed accessors
/// (one shared code path — bitwise-identical to the slice form).
///
/// Method: with `S_j` the Gaussian survival functions,
/// `E[min_j X_j] = L + ∫_L^U Π_j S_j(x) dx` where `L` truncates below
/// every design's `12σ` lower tail and `U` is the smallest `12σ` upper
/// edge (beyond `U` some survival factor is `< 2e-33`, so the product
/// vanishes; below `L` it is `1` to the same tolerance). The integral
/// is a fixed [`EOL_QUADRATURE_PANELS`]-panel composite Simpson rule —
/// deterministic, allocation-free, and monotone in menu content. The
/// result is an ESTIMATED value with quadrature/tail resolution around
/// `(U−L)/panels` curvature error, not a certified enclosure.
///
/// Domain conventions match the surrogate scan: designs with
/// non-finite means are excluded from the decision; fewer than two
/// finite designs mean no decision and zero loss; a zero-`σ` design
/// truncates the minimum exactly. A non-finite deviation on an
/// included design poisons the result to NaN, exactly as it poisons
/// the surrogate.
#[must_use]
pub fn expected_opportunity_loss_by(
    len: usize,
    mean_at: &dyn Fn(usize) -> f64,
    std_at: &dyn Fn(usize) -> f64,
) -> f64 {
    let Some((best, _)) = top_two_indices_by(len, mean_at) else {
        return 0.0;
    };
    let best_mean = mean_at(best);
    // Integration window over the FINITE designs (the same inclusion
    // rule as the scan), saturated so astronomically scaled menus keep
    // a representable window instead of overflowing to ±inf.
    let mut lower = f64::INFINITY;
    let mut upper = f64::INFINITY;
    for idx in 0..len {
        let mean = mean_at(idx);
        if !mean.is_finite() {
            continue;
        }
        let std = std_at(idx);
        if !(std.is_finite() && std >= 0.0) {
            return f64::NAN;
        }
        let half = (EOL_TAIL_HALF_WIDTHS * std).clamp(0.0, f64::MAX);
        lower = lower.min((mean - half).clamp(-f64::MAX, f64::MAX));
        upper = upper.min((mean + half).clamp(-f64::MAX, f64::MAX));
    }
    if !(lower.is_finite() && upper.is_finite()) || upper <= lower {
        // Degenerate window: every included design is exact (σ = 0) or
        // one exact design caps the minimum below every tail — the
        // minimum is deterministic at `lower`-side resolution.
        return (best_mean - best_mean.min(upper)).max(0.0);
    }
    // Survival product Π_j S_j(x); zero-σ designs are exact step
    // functions (their mean is ≥ `upper` by window construction, so
    // they contribute 1 across the open window).
    let survival_product = |x: f64| -> f64 {
        let mut product = 1.0f64;
        for idx in 0..len {
            let mean = mean_at(idx);
            if !mean.is_finite() {
                continue;
            }
            let std = std_at(idx);
            if std <= 0.0 {
                if x >= mean {
                    return 0.0;
                }
                continue;
            }
            product *= 1.0 - normal_cdf((x - mean) / std);
            if product <= 0.0 {
                return 0.0;
            }
        }
        product
    };
    // Composite Simpson over the fixed panel count.
    let panels = EOL_QUADRATURE_PANELS as f64;
    let step = (upper - lower) / panels;
    let mut integral = survival_product(lower) + survival_product(upper);
    for node in 1..EOL_QUADRATURE_PANELS {
        let weight = if node % 2 == 1 { 4.0 } else { 2.0 };
        integral += weight * survival_product(step.mul_add(node as f64, lower));
    }
    integral *= step / 3.0;
    let expected_min = lower + integral;
    (best_mean - expected_min).max(0.0)
}

/// A menu action's kind (which uncertainty it reduces).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActionKind {
    /// Run a cheaper surrogate — reduces numerical uncertainty.
    Surrogate,
    /// Run a higher-resolution simulation — reduces numerical uncertainty.
    Simulate,
    /// Refine the mesh near the decision boundary — reduces numerical.
    Refine,
    /// Sample more scenarios — reduces statistical uncertainty.
    Sample,
    /// Build / test a prototype — reduces model-form uncertainty.
    Test,
}

impl ActionKind {
    /// The uncertainty component this action reduces.
    #[must_use]
    pub fn component(self) -> Component {
        match self {
            ActionKind::Surrogate | ActionKind::Simulate | ActionKind::Refine => {
                Component::Numerical
            }
            ActionKind::Sample => Component::Statistical,
            ActionKind::Test => Component::Model,
        }
    }
}

/// A candidate information-acquisition action.
#[derive(Debug, Clone, PartialEq)]
pub struct Action {
    /// The action id.
    pub name: String,
    /// What it does.
    pub kind: ActionKind,
    /// Which design it informs.
    pub target_design: String,
    /// The fraction by which it reduces the target component (`0..=1`).
    pub reduction: f64,
    /// Its cost.
    pub cost: f64,
}

/// The decision value of an action.
#[derive(Debug, Clone, PartialEq)]
pub struct ActionValue {
    /// The action id.
    pub action: String,
    /// The EVPI reduction it buys (its decision value).
    pub value: f64,
    /// Its cost.
    pub cost: f64,
    /// Decision value per unit cost.
    pub value_per_cost: f64,
}

/// The decision value of `action`: how much it reduces the TOP-TWO
/// SURROGATE EVPI (per cost) — a deliberate ranking heuristic, cheap
/// enough to evaluate per action. An action on a design outside the
/// decision boundary reduces it by ~0 — worthless however uncertain
/// that design is. Surrogate caveat: this ranks actions; it never
/// claims the decision is globally robust (that gate is
/// [`expected_opportunity_loss`] inside [`recommend`]).
#[must_use]
pub fn action_value(designs: &[DesignEstimate], action: &Action) -> ActionValue {
    let before = top_two_evpi_surrogate(designs);
    let after_designs: Vec<DesignEstimate> = designs
        .iter()
        .map(|d| {
            if d.name == action.target_design {
                DesignEstimate {
                    uncertainty: d
                        .uncertainty
                        .reduced(action.kind.component(), action.reduction),
                    ..d.clone()
                }
            } else {
                d.clone()
            }
        })
        .collect();
    let after = top_two_evpi_surrogate(&after_designs);
    let value = (before - after).max(0.0);
    let value_per_cost = if action.cost.is_finite() && action.cost > 0.0 {
        value / action.cost
    } else if action.cost.is_finite()
        && action.cost >= 0.0
        && action.cost <= f64::EPSILON
        && value > 0.0
    {
        f64::INFINITY
    } else {
        0.0
    };
    ActionValue {
        action: action.name.clone(),
        value,
        cost: action.cost,
        value_per_cost,
    }
}

/// The VOI recommendation.
#[derive(Debug, Clone, PartialEq)]
pub enum Recommendation {
    /// Take this action (best decision value per cost).
    Act {
        /// The action id.
        action: String,
        /// Its decision value per cost.
        value_per_cost: f64,
    },
    /// Stop — the decision is already robust (EVPI below the threshold).
    Stop {
        /// Why.
        reason: String,
    },
}

/// Recommend the highest decision-value-per-cost action, or STOP when the
/// decision is already robust (current EVPI `<= stop_threshold`).
///
/// The STOP gate is the FULL multi-alternative
/// [`expected_opportunity_loss`] (bead sj31i.5): a high-variance
/// third alternative keeps the campaign alive even when the top-two
/// surrogate reads zero. Action ranking below the gate remains
/// surrogate-driven.
#[must_use]
pub fn recommend(
    designs: &[DesignEstimate],
    actions: &[Action],
    stop_threshold: f64,
) -> Recommendation {
    let current = expected_opportunity_loss(designs);
    if current <= stop_threshold {
        return Recommendation::Stop {
            reason: format!("decision robust: EVPI {current:.3e} <= {stop_threshold:.3e}"),
        };
    }
    let best = actions
        .iter()
        .map(|a| action_value(designs, a))
        .filter(|v| v.value > 0.0 && v.value_per_cost > 0.0)
        .max_by(|a, b| {
            a.value_per_cost
                .partial_cmp(&b.value_per_cost)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    match best {
        Some(v) => Recommendation::Act {
            action: v.action,
            value_per_cost: v.value_per_cost,
        },
        None => Recommendation::Stop {
            reason: "no action changes the decision".to_string(),
        },
    }
}

/// The [M] baseline VOI must beat: uncertainty-proportional allocation — pick
/// the action informing the design with the LARGEST total uncertainty,
/// regardless of decision relevance.
#[must_use]
pub fn heuristic_choice<'a>(
    designs: &[DesignEstimate],
    actions: &'a [Action],
) -> Option<&'a Action> {
    let most_uncertain = designs.iter().max_by(|a, b| {
        a.std()
            .partial_cmp(&b.std())
            .unwrap_or(std::cmp::Ordering::Equal)
    })?;
    actions
        .iter()
        .find(|a| a.target_design == most_uncertain.name)
}

// -- Gaussian helpers (in-house; no external special functions) -------------

/// Overflow-safe Euclidean norm (bead sj31i.5): the largest finite
/// magnitude is factored out before squaring, so finite components near
/// `√MAX` compose to their representable norm instead of overflowing
/// the naive square sum. Non-finite components propagate (`inf` wins
/// over `NaN`, matching IEEE hypot semantics closely enough for
/// deviation composition; a lone `NaN` poisons the result).
fn scaled_norm(components: &[f64]) -> f64 {
    let mut scale = 0.0f64;
    for &component in components {
        if component.is_infinite() {
            return f64::INFINITY;
        }
        scale = scale.max(component.abs());
    }
    if scale == 0.0 {
        // All zeros, or zeros mixed with NaN (`max` skips NaN): the
        // fold below yields exact zero or the poisoning NaN.
        return components.iter().fold(0.0f64, |acc, &c| acc + c * c);
    }
    let mut sum = 0.0f64;
    for &component in components {
        let scaled = component / scale;
        sum += scaled * scaled;
    }
    scale * sum.sqrt()
}

fn normal_pdf(z: f64) -> f64 {
    (-0.5 * z * z).exp() / std::f64::consts::TAU.sqrt()
}

fn normal_cdf(z: f64) -> f64 {
    f64::midpoint(1.0, erf(z / std::f64::consts::SQRT_2))
}

/// `erf` via Abramowitz & Stegun 7.1.26 (published constants, ≤ 1.5e-7 error).
#[allow(clippy::unreadable_literal)]
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    sign * y
}
