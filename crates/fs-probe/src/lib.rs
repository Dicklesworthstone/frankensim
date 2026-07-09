//! fs-probe — discrepancy probes + the budget pie (plan addendum, Proposal 3).
//! Layer: L3. Makes the THIRD epistemic color (`estimated`) actually
//! COMPUTABLE, and makes model-form error visible.
//!
//! Two mechanisms:
//!
//! 1. **Adjacent-rung discrepancy probes.** Given the shared fidelity-ladder
//!    registry ([`fs_ladder`]), evaluate the same design on two ADJACENT rungs,
//!    prolongate the coarse solution onto the fine rung, and take the
//!    difference: `|fine − prolongate(coarse)|` per subdomain is a MEASURED,
//!    LOCALIZED model-form error field. It is stored **estimated-color** (an
//!    [`fs_evidence::Color::Estimated`] with the probe as estimator and its
//!    dispersion). Humans never do this systematically because it is tedious;
//!    a swarm does it overnight.
//!
//! 2. **The budget pie.** Aggregate a result's error contributions BY COLOR:
//!    verified (interval-certified numerics ≈ discretization/mesh error),
//!    validated (regime-anchored), estimated (model-form: closure/surrogate).
//!    The pie tells the operator where the error budget is ACTUALLY spent — so
//!    "refine the mesh" is never prescribed when the budget says "your closure
//!    is the problem". This is the single most operator-legible artifact of the
//!    whole epistemic type system.
//!
//! Probe compute is bounded by a [`ProbeBudget`] (a HARD ceiling as a fixed
//! fraction of the fleet budget). Everything here is a pure, deterministic
//! function of its inputs (no RNG, no I/O), so probe fields and budget pies are
//! bit-reproducible on replay.

use fs_evidence::{Color, ColorRank};
use fs_ladder::{Ladder, LadderError};
use std::error::Error;
use std::fmt;

/// A structured probe error (a refusal that teaches).
#[derive(Debug, Clone, PartialEq)]
pub enum ProbeError {
    /// The ladder rejected the prolongation (bad rung, at the top, …).
    Ladder(LadderError),
    /// The prolongated coarse state and the fine state have different lengths
    /// — they are not on the same rung's grid.
    DimMismatch {
        /// Length after prolongation.
        prolongated_len: usize,
        /// Length of the supplied fine state.
        fine_len: usize,
    },
    /// A probe would exceed the fleet-budget cap.
    BudgetExceeded {
        /// The requested cost.
        requested: f64,
        /// The remaining budget under the cap.
        remaining: f64,
    },
    /// A non-finite or negative cost was requested.
    BadCost {
        /// The offending cost.
        cost: f64,
    },
}

impl fmt::Display for ProbeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProbeError::Ladder(e) => write!(f, "ladder rejected the probe: {e}"),
            ProbeError::DimMismatch {
                prolongated_len,
                fine_len,
            } => write!(
                f,
                "prolongated coarse state has {prolongated_len} points but the fine state has \
                 {fine_len}; fix: evaluate both on the same adjacent-rung grid"
            ),
            ProbeError::BudgetExceeded {
                requested,
                remaining,
            } => write!(
                f,
                "probe cost {requested} exceeds the {remaining} remaining under the fleet-budget \
                 cap; fix: raise the cap or skip the probe"
            ),
            ProbeError::BadCost { cost } => {
                write!(f, "probe cost {cost} must be finite and non-negative")
            }
        }
    }
}

impl Error for ProbeError {}

impl From<LadderError> for ProbeError {
    fn from(e: LadderError) -> ProbeError {
        ProbeError::Ladder(e)
    }
}

/// A measured, localized model-form discrepancy between two adjacent rungs.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscrepancyField {
    /// The kernel whose ladder was probed.
    pub kernel: String,
    /// The coarse rung index.
    pub from_rung: u32,
    /// The fine rung index (`from_rung + 1`).
    pub to_rung: u32,
    /// The per-subdomain absolute discrepancy `|fine − prolongate(coarse)|`.
    pub per_subdomain: Vec<f64>,
    /// The largest discrepancy (the headline magnitude).
    pub l_inf: f64,
    /// The mean discrepancy.
    pub mean: f64,
    /// The epistemic color — always `Estimated` (a model-form probe is not a
    /// certified bound): the estimator identity + its dispersion.
    pub color: Color,
}

/// Run an adjacent-rung discrepancy probe: prolongate the `from_rung` coarse
/// solution onto `from_rung + 1` and measure the model-form gap against the
/// fine solution. The result is estimated-color; a near-zero gap gives a
/// near-zero dispersion (no manufactured error).
///
/// # Errors
/// [`ProbeError::Ladder`] if the ladder cannot prolongate at `from_rung`;
/// [`ProbeError::DimMismatch`] if the prolongated and fine states differ in
/// length.
pub fn probe_adjacent(
    ladder: &Ladder,
    from_rung: u32,
    coarse: &[f64],
    fine: &[f64],
) -> Result<DiscrepancyField, ProbeError> {
    let prolongated = ladder.prolongate(from_rung, coarse)?;
    if prolongated.len() != fine.len() {
        return Err(ProbeError::DimMismatch {
            prolongated_len: prolongated.len(),
            fine_len: fine.len(),
        });
    }
    let per_subdomain: Vec<f64> = fine
        .iter()
        .zip(&prolongated)
        .map(|(f, p)| (f - p).abs())
        .collect();
    let l_inf = per_subdomain.iter().copied().fold(0.0_f64, f64::max);
    let mean = if per_subdomain.is_empty() {
        0.0
    } else {
        per_subdomain.iter().sum::<f64>() / per_subdomain.len() as f64
    };
    let color = Color::Estimated {
        estimator: format!(
            "adjacent-rung-probe:{}:{from_rung}->{}",
            ladder.kernel(),
            from_rung + 1
        ),
        dispersion: l_inf,
    };
    Ok(DiscrepancyField {
        kernel: ladder.kernel().to_string(),
        from_rung,
        to_rung: from_rung + 1,
        per_subdomain,
        l_inf,
        mean,
        color,
    })
}

/// One error contribution to a result: its source, its epistemic color, and
/// its magnitude (in the QoI's error units).
#[derive(Debug, Clone, PartialEq)]
pub struct ErrorContribution {
    /// What produced this error (e.g. `"mesh"`, `"turbulence-closure"`).
    pub source: String,
    /// Its epistemic color.
    pub color: Color,
    /// Its magnitude (non-negative).
    pub magnitude: f64,
}

impl ErrorContribution {
    /// A contribution.
    #[must_use]
    pub fn new(source: impl Into<String>, color: Color, magnitude: f64) -> ErrorContribution {
        ErrorContribution {
            source: source.into(),
            color,
            magnitude,
        }
    }
}

/// The budget pie: where the error budget is spent, by color.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BudgetPie {
    /// Total error magnitude.
    pub total: f64,
    /// Verified (interval-certified numerics ≈ discretization/mesh) error.
    pub verified: f64,
    /// Validated (regime-anchored) error.
    pub validated: f64,
    /// Estimated (model-form: closure / surrogate) error.
    pub estimated: f64,
}

impl BudgetPie {
    /// Aggregate contributions into a pie (magnitudes summed by color rank).
    #[must_use]
    pub fn of(contributions: &[ErrorContribution]) -> BudgetPie {
        let mut pie = BudgetPie {
            total: 0.0,
            verified: 0.0,
            validated: 0.0,
            estimated: 0.0,
        };
        for c in contributions {
            let m = c.magnitude.max(0.0);
            pie.total += m;
            match c.color.rank() {
                ColorRank::Verified => pie.verified += m,
                ColorRank::Validated => pie.validated += m,
                ColorRank::Estimated => pie.estimated += m,
            }
        }
        pie
    }

    /// The fraction of the total budget in a given color (0 if the total is 0).
    #[must_use]
    pub fn fraction(&self, rank: ColorRank) -> f64 {
        if self.total <= 0.0 {
            return 0.0;
        }
        let part = match rank {
            ColorRank::Verified => self.verified,
            ColorRank::Validated => self.validated,
            ColorRank::Estimated => self.estimated,
        };
        part / self.total
    }

    /// The color that dominates the budget (`None` if there is no budget).
    /// Ties break toward the WEAKER color (estimated > validated > verified),
    /// so an ambiguous case is reported conservatively as the weaker claim.
    #[must_use]
    pub fn dominant(&self) -> Option<ColorRank> {
        if self.total <= 0.0 {
            return None;
        }
        // start from the weakest so equal magnitudes resolve conservatively.
        let mut best = ColorRank::Estimated;
        let mut best_mag = self.estimated;
        if self.validated > best_mag {
            best = ColorRank::Validated;
            best_mag = self.validated;
        }
        if self.verified > best_mag {
            best = ColorRank::Verified;
        }
        Some(best)
    }

    /// An operator-legible verdict: where the budget says to spend effort.
    #[must_use]
    pub fn verdict(&self) -> &'static str {
        match self.dominant() {
            None => "no error budget recorded",
            Some(ColorRank::Verified) => {
                "numerical/discretization error dominates — refine the mesh or raise the order"
            }
            Some(ColorRank::Validated) => {
                "validated (regime-anchored) error dominates — check the regime and its anchor"
            }
            Some(ColorRank::Estimated) => {
                "MODEL-FORM (closure/surrogate) error dominates — refining the mesh will NOT help; \
                 fix the closure or validate it"
            }
        }
    }
}

/// A HARD ceiling on probe compute, as a fixed fraction of the fleet budget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProbeBudget {
    fleet_budget: f64,
    cap_fraction: f64,
    spent: f64,
}

impl ProbeBudget {
    /// A probe budget capping spend at `cap_fraction` of `fleet_budget`.
    /// `cap_fraction` is clamped to `[0, 1]`; non-finite inputs fail closed
    /// to a zero cap.
    #[must_use]
    pub fn new(fleet_budget: f64, cap_fraction: f64) -> ProbeBudget {
        let fleet_budget = if fleet_budget.is_finite() {
            fleet_budget.max(0.0)
        } else {
            0.0
        };
        let cap_fraction = if cap_fraction.is_finite() {
            cap_fraction.clamp(0.0, 1.0)
        } else {
            0.0
        };
        ProbeBudget {
            fleet_budget,
            cap_fraction,
            spent: 0.0,
        }
    }

    /// The absolute cap (`cap_fraction × fleet_budget`).
    #[must_use]
    pub fn cap(&self) -> f64 {
        self.fleet_budget * self.cap_fraction
    }

    /// How much has been spent on probes so far.
    #[must_use]
    pub fn spent(&self) -> f64 {
        self.spent
    }

    /// The remaining probe budget under the cap.
    #[must_use]
    pub fn remaining(&self) -> f64 {
        (self.cap() - self.spent).max(0.0)
    }

    /// Try to spend `cost` on a probe. The cap is a HARD ceiling: spending up
    /// to EXACTLY the cap is allowed; anything beyond is refused.
    ///
    /// # Errors
    /// [`ProbeError::BadCost`] if `cost` is negative or non-finite;
    /// [`ProbeError::BudgetExceeded`] if it would exceed the cap.
    pub fn try_spend(&mut self, cost: f64) -> Result<(), ProbeError> {
        if !cost.is_finite() || cost < 0.0 {
            return Err(ProbeError::BadCost { cost });
        }
        if self.spent + cost > self.cap() {
            return Err(ProbeError::BudgetExceeded {
                requested: cost,
                remaining: self.remaining(),
            });
        }
        self.spent += cost;
        Ok(())
    }
}
