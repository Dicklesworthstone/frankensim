//! Horizon trigger 4 (bead `frankensim-epic-addendum-xpck.5.4`): the
//! paying-workload splitting-error demand gate for proposal 4 ("Extend the
//! complex into time", [`crate::proposals`]).
//!
//! PURE decision logic over measured coupled-transient error budgets. The
//! activation law is deliberately narrow — adaptive spacetime control may
//! switch on ONLY when, for a PAYING workload, the splitting-error term is
//! simultaneously (a) the LARGEST actionable term and (b) at least 20% of
//! the complete budget, with stable coupling. Everything else measures and
//! stays instrument-only; ties and equality at the threshold do NOT fire.
//!
//! Honesty laws enforced here:
//! - Budget completeness: shares must be present for every required term
//!   (including `splitting` itself), in `[0,1]`, summing to 1 within
//!   tolerance — a dropped larger term is a REFUSAL, never a verdict.
//! - Non-paying workloads are measured but can never activate.
//! - Unstable coupling is nonactivating by definition (the controller's
//!   premise dies with the stability assumption).
//! - An empty retained population yields [`NoData`] — absence of paying
//!   workloads is a measurement gap, not a green.

/// The activation threshold: splitting error must be AT LEAST this share of
/// the complete error budget (equality alone does not activate — see
/// [`splitting_verdict`]: it must ALSO be strictly the largest term).
pub const SPLITTING_SHARE_MIN: f64 = 0.20;

/// Tolerance on budget shares summing to one.
const SUM_TOLERANCE: f64 = 1e-9;

/// The required budget term names (a complete budget binds every
/// independently actionable error source; dropping one would flatter the
/// splitting share).
pub const REQUIRED_TERMS: [&str; 4] = ["splitting", "discretization", "iteration", "model"];

/// One actionable error term: an attributed fraction of the workload's
/// complete error budget.
#[derive(Debug, Clone, PartialEq)]
pub struct ErrorTerm {
    pub name: String,
    /// Fraction of the total budget in `[0, 1]`.
    pub share: f64,
}

impl ErrorTerm {
    #[must_use]
    pub fn new(name: &str, share: f64) -> Self {
        Self { name: name.to_string(), share }
    }
}

/// Retention class of a bound workload. Only [`WorkloadClass::Paying`]
/// workloads can ever satisfy the demand gate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkloadClass {
    Paying,
    NonPaying,
}

/// A bound workload's complete error budget plus its coupling-stability
/// attestation.
#[derive(Debug, Clone)]
pub struct WorkloadBudget {
    pub workload_id: String,
    pub class: WorkloadClass,
    pub terms: Vec<ErrorTerm>,
    /// Attestation that the coupled system was numerically stable over the
    /// measured window (an unstable run invalidates the budget split).
    pub coupling_stable: bool,
}

/// Typed refusals: every malformed budget refuses BY NAME. None of these
/// are verdicts; all are nonactivating.
#[derive(Debug, Clone, PartialEq)]
pub enum BudgetRefusal {
    /// No terms at all — a zero budget is a measurement gap.
    EmptyTerms,
    /// A share outside `[0, 1]`.
    ShareOutOfRange { term: String, share: f64 },
    /// Shares do not sum to one within [`SUM_TOLERANCE`] (double-counted or
    /// dropped error).
    SharesDontSum { sum: f64 },
    /// A [`REQUIRED_TERMS`] name absent from the budget.
    MissingRequiredTerm { term: String },
}

/// Per-workload verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplittingVerdict {
    /// Gate satisfied: paying workload, stable coupling, splitting error
    /// strictly dominant and >= 20% of the complete budget.
    Activate,
    /// Measured, attributed, retained — but not activated.
    InstrumentOnly,
}

impl SplittingVerdict {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Activate => "Activate",
            Self::InstrumentOnly => "InstrumentOnly",
        }
    }
}

/// Adjudicate one workload budget against the demand gate.
///
/// Refusals first (fail closed), then nonactivation classes, then the gate:
/// `splitting >= SPLITTING_SHARE_MIN` AND `splitting > max(other terms)`.
#[must_use]
pub fn splitting_verdict(w: &WorkloadBudget) -> Result<SplittingVerdict, BudgetRefusal> {
    if w.terms.is_empty() {
        return Err(BudgetRefusal::EmptyTerms);
    }
    let mut sum = 0.0;
    for t in &w.terms {
        if !(0.0..=1.0).contains(&t.share) {
            return Err(BudgetRefusal::ShareOutOfRange { term: t.name.clone(), share: t.share });
        }
        sum += t.share;
    }
    if (sum - 1.0).abs() > SUM_TOLERANCE {
        return Err(BudgetRefusal::SharesDontSum { sum });
    }
    for required in REQUIRED_TERMS {
        if !w.terms.iter().any(|t| t.name == required) {
            return Err(BudgetRefusal::MissingRequiredTerm { term: required.to_string() });
        }
    }
    if !w.coupling_stable {
        return Ok(SplittingVerdict::InstrumentOnly);
    }
    if w.class != WorkloadClass::Paying {
        return Ok(SplittingVerdict::InstrumentOnly);
    }
    let splitting = w
        .terms
        .iter()
        .find(|t| t.name == "splitting")
        .map_or(0.0, |t| t.share);
    let largest_other = w
        .terms
        .iter()
        .filter(|t| t.name != "splitting")
        .map(|t| t.share)
        .fold(0.0_f64, f64::max);
    if splitting >= SPLITTING_SHARE_MIN && splitting > largest_other {
        Ok(SplittingVerdict::Activate)
    } else {
        Ok(SplittingVerdict::InstrumentOnly)
    }
}

/// Aggregate disposition over the retained workload population.
#[derive(Debug, Clone, PartialEq)]
pub enum PopulationDisposition {
    /// No retained population (or no paying member): the trigger cannot even
    /// be evaluated. This is today's honest state.
    NoData { reason: String },
    /// Every paying workload passed the gate.
    Activate { verdicts: Vec<(String, SplittingVerdict)> },
    /// Bound and measured, but at least one paying workload stayed below
    /// the gate (or refused — refusals are surfaced, never swallowed).
    InstrumentOnly { verdicts: Vec<(String, SplittingVerdict)> },
}

/// Aggregate a retained population. Weakest link: ANY refusing budget or
/// any non-activating paying workload holds the population at
/// InstrumentOnly. No paying members at all => NoData.
#[must_use]
pub fn population_disposition(workloads: &[WorkloadBudget]) -> PopulationDisposition {
    let mut verdicts = Vec::new();
    let mut saw_paying = false;
    let mut all_activated = true;
    for w in workloads {
        if w.class == WorkloadClass::Paying {
            saw_paying = true;
        }
        match splitting_verdict(w) {
            Ok(v) => {
                if w.class == WorkloadClass::Paying && v != SplittingVerdict::Activate {
                    all_activated = false;
                }
                verdicts.push((w.workload_id.clone(), v));
            }
            Err(r) => {
                // A refusing budget blocks activation for the whole
                // population: incomplete evidence never activates.
                all_activated = false;
                verdicts.push((
                    format!("{} [refused]", w.workload_id),
                    SplittingVerdict::InstrumentOnly,
                ));
                let _ = r;
            }
        }
    }
    if !saw_paying {
        return PopulationDisposition::NoData {
            reason: "no retained paying coupled-transient workload is bound".to_string(),
        };
    }
    if all_activated {
        PopulationDisposition::Activate { verdicts }
    } else {
        PopulationDisposition::InstrumentOnly { verdicts }
    }
}

/// Today's receipt: the retained paying-workload population is EMPTY, so
/// proposal 4's demand gate is NoData by construction. The engine above is
/// live and tested; the controller stays off until real budgets bind.
#[must_use]
pub fn current_receipt() -> PopulationDisposition {
    PopulationDisposition::NoData {
        reason: "no paying coupled-transient workload exists in the program yet \
                 (fs-wedge has no bound production coupled-transient customer)"
            .to_string(),
    }
}
