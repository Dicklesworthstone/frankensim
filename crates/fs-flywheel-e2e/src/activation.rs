//! Activation statistics for the flywheel checkpoint (bead
//! frankensim-sj31i.18).
//!
//! Root cause this module closes: the phase-1 "six-month checkpoint"
//! treated 20 correlated grid points from a deliberately linear θ-family
//! (with the proposer cache drawn from the same family) as
//! customer-realistic activation evidence, and an empty warm-start
//! ratio list collapsed to `+INFINITY` median savings — which passed
//! the ≥ 1.5× gate. That fixture remains valuable as CONFORMANCE (the
//! machinery runs, the hostile control can fail); it is not activation
//! evidence and can no longer masquerade as such.
//!
//! Activation now requires:
//! - a preregistered [`SamplingFrame`] — named strata with minimum
//!   sample counts, a seed, and the thresholds — whose content identity
//!   is retained in the verdict (no post-result threshold editing);
//! - independent holdout evidence only: any development-tagged row is a
//!   typed [`ActivationRefusal::HoldoutLeakage`];
//! - exact denominators: a warm-started sample with a zero warm cost is
//!   refused, an all-outright stratum has UNMEASURED savings and
//!   refuses instead of minting infinite ratios, and duplicate sample
//!   identities cannot inflate `n`;
//! - anytime-valid statistics: each gate is a betting e-process
//!   (Ville's inequality bounds the false-activation rate at the
//!   preregistered δ under optional stopping), so stopping early or
//!   late cannot manufacture a pass;
//! - stratified adjudication: EVERY stratum must reject both null
//!   hypotheses — an adverse stratum cannot hide inside a favorable
//!   pool.
//!
//! `AcceptedOutright` and measured savings are distinct: outright
//! accepts feed the accept-rate gate only, while the savings gate runs
//! exclusively over warm-started samples with exact cold/warm counts.

use std::collections::BTreeMap;

/// One preregistered stratum: a name and the minimum holdout sample
/// count below which adjudication refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StratumSpec {
    /// Stable stratum name (kernel class × regime, for example).
    pub name: &'static str,
    /// Minimum admitted holdout samples; below this the frame refuses.
    pub min_samples: u32,
}

/// Preregistered activation thresholds. All fields are part of the
/// frame identity; editing any of them after results exist produces a
/// different frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActivationThresholds {
    /// H0 for the accept gate: true accept probability ≤ this floor.
    pub accept_rate_floor: f64,
    /// Savings gate: a warm-started sample counts as a win when
    /// `cold / warm ≥ savings_floor`; H0 is P(win) ≤ 1/2 (a median
    /// claim).
    pub savings_floor: f64,
    /// Anytime false-activation bound δ ∈ (0, 1); each e-process must
    /// reach 1/δ.
    pub confidence_delta: f64,
    /// Betting fraction λ ∈ (0, 1) applied within each e-process step
    /// (scaled so wealth stays positive for either outcome).
    pub betting_fraction: f64,
}

/// Whether a sample was drawn for development or held out for
/// adjudication. Only holdout rows are admissible activation evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Partition {
    /// Used to build/tune proposers; never activation evidence.
    Development,
    /// Independent holdout; the only admissible evidence.
    Holdout,
}

/// The outcome of one speculation query, with exact denominators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// The proposal was accepted outright — no warm-start measurement
    /// exists for this sample, so it carries NO savings evidence.
    AcceptedOutright,
    /// A warm start ran: exact cold and warm work counts.
    WarmStarted {
        /// Work units the cold solve would have used.
        cold: u32,
        /// Work units the warm-started solve used.
        warm: u32,
    },
    /// The proposal was rejected and the solve ran cold.
    ColdSolve,
}

/// One admitted evidence row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sample {
    /// Stable sample identity (hash of problem/θ/kernel/regime).
    /// Duplicates are refused rather than double-counted.
    pub key: u64,
    /// The preregistered stratum this sample belongs to.
    pub stratum: &'static str,
    /// Development or holdout.
    pub partition: Partition,
    /// The measured outcome.
    pub outcome: Outcome,
}

/// The preregistered sampling frame: strata, seed, and thresholds.
#[derive(Debug, Clone, PartialEq)]
pub struct SamplingFrame {
    /// Study name binding the frame to its corpus.
    pub study: &'static str,
    /// The corpus draw seed (part of identity; replay uses it).
    pub seed: u64,
    /// Preregistered strata with minimum counts.
    pub strata: Vec<StratumSpec>,
    /// Preregistered thresholds.
    pub thresholds: ActivationThresholds,
}

impl SamplingFrame {
    /// Content identity of the frame: any edit to strata, seed, or
    /// thresholds is a different frame. Retained in every verdict.
    #[must_use]
    pub fn identity(&self) -> fs_ledger::ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"fs-flywheel-e2e.sampling-frame.v1\x00");
        bytes.extend_from_slice(self.study.as_bytes());
        bytes.push(0);
        bytes.extend_from_slice(&self.seed.to_le_bytes());
        bytes.extend_from_slice(
            &u64::try_from(self.strata.len())
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        for stratum in &self.strata {
            bytes.extend_from_slice(stratum.name.as_bytes());
            bytes.push(0);
            bytes.extend_from_slice(&stratum.min_samples.to_le_bytes());
        }
        bytes.extend_from_slice(&self.thresholds.accept_rate_floor.to_le_bytes());
        bytes.extend_from_slice(&self.thresholds.savings_floor.to_le_bytes());
        bytes.extend_from_slice(&self.thresholds.confidence_delta.to_le_bytes());
        bytes.extend_from_slice(&self.thresholds.betting_fraction.to_le_bytes());
        fs_ledger::hash_bytes(&bytes)
    }
}

/// Structural refusals: the evidence cannot even be adjudicated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActivationRefusal {
    /// A threshold field is outside its admitted domain.
    MalformedThresholds,
    /// A frame declares no strata or a zero minimum.
    MalformedFrame,
    /// A sample names a stratum the frame never preregistered.
    UnknownStratum {
        /// The offending sample identity.
        key: u64,
    },
    /// A development-tagged row reached the adjudicator.
    HoldoutLeakage {
        /// The offending sample identity.
        key: u64,
    },
    /// Two rows share one sample identity; correlated duplicates
    /// cannot inflate the denominator.
    DuplicateSample {
        /// The duplicated identity.
        key: u64,
    },
    /// A warm-started row with a zero warm denominator (the +∞ hole).
    ZeroWarmDenominator {
        /// The offending sample identity.
        key: u64,
    },
    /// A stratum has fewer holdout samples than it preregistered.
    InsufficientSamples {
        /// The short stratum.
        stratum: &'static str,
        /// Admitted rows.
        have: u32,
        /// Preregistered minimum.
        need: u32,
    },
    /// A stratum has no warm-started rows at all: its savings are
    /// unmeasured, and unmeasured savings never activate.
    NoWarmStartEvidence {
        /// The unmeasured stratum.
        stratum: &'static str,
    },
}

impl core::fmt::Display for ActivationRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MalformedThresholds => f.write_str("activation thresholds outside domain"),
            Self::MalformedFrame => f.write_str("sampling frame has no usable strata"),
            Self::UnknownStratum { key } => {
                write!(f, "sample {key:#x} names an unregistered stratum")
            }
            Self::HoldoutLeakage { key } => write!(
                f,
                "sample {key:#x} is development-tagged; holdout evidence only"
            ),
            Self::DuplicateSample { key } => {
                write!(f, "sample identity {key:#x} appears more than once")
            }
            Self::ZeroWarmDenominator { key } => {
                write!(
                    f,
                    "sample {key:#x} reports a warm start with zero warm work"
                )
            }
            Self::InsufficientSamples {
                stratum,
                have,
                need,
            } => write!(
                f,
                "stratum '{stratum}' has {have} holdout samples; {need} preregistered"
            ),
            Self::NoWarmStartEvidence { stratum } => write!(
                f,
                "stratum '{stratum}' has no warm-started rows; savings are unmeasured"
            ),
        }
    }
}

impl core::error::Error for ActivationRefusal {}

/// Per-stratum adjudication record retained in the report.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StratumReport {
    /// The stratum name.
    pub stratum: &'static str,
    /// Admitted holdout rows.
    pub samples: u32,
    /// Rows where the proposal was adopted (outright or warm).
    pub accepts: u32,
    /// Warm-started rows (the savings denominator).
    pub warm_samples: u32,
    /// Warm-started rows meeting the savings floor.
    pub savings_wins: u32,
    /// Running-maximum e-value against H0: accept rate ≤ floor.
    pub accept_e_max: f64,
    /// Running-maximum e-value against H0: P(win) ≤ 1/2.
    pub savings_e_max: f64,
    /// Whether the accept gate rejected its null at 1/δ.
    pub accept_rejected: bool,
    /// Whether the savings gate rejected its null at 1/δ.
    pub savings_rejected: bool,
}

/// The replayable activation verdict.
#[derive(Debug, Clone, PartialEq)]
pub struct ActivationReport {
    /// Identity of the preregistered frame this verdict binds to.
    pub frame_identity: fs_ledger::ContentHash,
    /// Per-stratum records, in frame order.
    pub strata: Vec<StratumReport>,
    /// True only when EVERY stratum rejected BOTH nulls.
    pub activated: bool,
}

/// One-sided Bernoulli e-process against H0: p ≤ p0.
///
/// Wealth update per observation x ∈ {0, 1}:
/// `W ← W · (1 + λ · (x − p0) / m)` with `m = max(p0, 1 − p0)` so the
/// per-step factor stays in (0, 2) for λ ∈ (0, 1). Ville's inequality
/// gives `P(sup W ≥ 1/δ) ≤ δ` under H0 at the boundary, uniformly over
/// stopping times. Returns the running-maximum wealth; a caller
/// comparing it against `1/δ` gets an anytime-valid level-δ test.
/// Public so paired comparisons elsewhere in the harness (the phase-2
/// adjoint-vs-baseline gate, bead sj31i.28) share one preregistered
/// statistical primitive.
#[must_use]
pub fn e_process(observations: impl Iterator<Item = bool>, p0: f64, lambda: f64) -> f64 {
    let scale = p0.max(1.0 - p0);
    let mut wealth = 1.0_f64;
    let mut max_wealth = 1.0_f64;
    for x in observations {
        let x = if x { 1.0 } else { 0.0 };
        wealth *= 1.0 + lambda * (x - p0) / scale;
        if wealth > max_wealth {
            max_wealth = wealth;
        }
    }
    max_wealth
}

fn validate(frame: &SamplingFrame) -> Result<(), ActivationRefusal> {
    let t = frame.thresholds;
    let thresholds_ok = t.accept_rate_floor > 0.0
        && t.accept_rate_floor < 1.0
        && t.savings_floor.is_finite()
        && t.savings_floor >= 1.0
        && t.confidence_delta > 0.0
        && t.confidence_delta < 1.0
        && t.betting_fraction > 0.0
        && t.betting_fraction < 1.0;
    if !thresholds_ok {
        return Err(ActivationRefusal::MalformedThresholds);
    }
    if frame.strata.is_empty() || frame.strata.iter().any(|s| s.min_samples == 0) {
        return Err(ActivationRefusal::MalformedFrame);
    }
    Ok(())
}

/// Adjudicate holdout evidence against a preregistered frame.
///
/// Pure and deterministic: samples are canonicalized by identity, so
/// the verdict is a function of the (frame, evidence-set) pair alone —
/// the replayable statistical claim behind every threshold.
///
/// # Errors
/// Typed [`ActivationRefusal`] for structural defects (unknown or
/// short strata, leakage, duplicates, zero denominators, unmeasured
/// savings). Statistical shortfall is NOT an error: it returns
/// `activated: false` with the per-stratum e-values retained.
pub fn adjudicate(
    frame: &SamplingFrame,
    samples: &[Sample],
) -> Result<ActivationReport, ActivationRefusal> {
    validate(frame)?;
    let registered: BTreeMap<&'static str, u32> = frame
        .strata
        .iter()
        .map(|s| (s.name, s.min_samples))
        .collect();

    // Canonicalize by identity; refuse duplicates and structural rot.
    let mut canonical: BTreeMap<u64, Sample> = BTreeMap::new();
    for sample in samples {
        if !registered.contains_key(sample.stratum) {
            return Err(ActivationRefusal::UnknownStratum { key: sample.key });
        }
        if sample.partition == Partition::Development {
            return Err(ActivationRefusal::HoldoutLeakage { key: sample.key });
        }
        if let Outcome::WarmStarted { warm: 0, .. } = sample.outcome {
            return Err(ActivationRefusal::ZeroWarmDenominator { key: sample.key });
        }
        if canonical.insert(sample.key, *sample).is_some() {
            return Err(ActivationRefusal::DuplicateSample { key: sample.key });
        }
    }

    let t = frame.thresholds;
    let target = 1.0 / t.confidence_delta;
    let mut reports = Vec::with_capacity(frame.strata.len());
    let mut activated = true;
    for spec in &frame.strata {
        let rows: Vec<&Sample> = canonical
            .values()
            .filter(|s| s.stratum == spec.name)
            .collect();
        let have = u32::try_from(rows.len()).unwrap_or(u32::MAX);
        if have < spec.min_samples {
            return Err(ActivationRefusal::InsufficientSamples {
                stratum: spec.name,
                have,
                need: spec.min_samples,
            });
        }
        let accepts = rows
            .iter()
            .filter(|s| {
                matches!(
                    s.outcome,
                    Outcome::AcceptedOutright | Outcome::WarmStarted { .. }
                )
            })
            .count();
        let warm_rows: Vec<(u32, u32)> = rows
            .iter()
            .filter_map(|s| match s.outcome {
                Outcome::WarmStarted { cold, warm } => Some((cold, warm)),
                _ => None,
            })
            .collect();
        if warm_rows.is_empty() {
            return Err(ActivationRefusal::NoWarmStartEvidence { stratum: spec.name });
        }
        let savings_wins = warm_rows
            .iter()
            .filter(|(cold, warm)| f64::from(*cold) / f64::from(*warm) >= t.savings_floor)
            .count();

        let accept_e_max = e_process(
            rows.iter().map(|s| {
                matches!(
                    s.outcome,
                    Outcome::AcceptedOutright | Outcome::WarmStarted { .. }
                )
            }),
            t.accept_rate_floor,
            t.betting_fraction,
        );
        let savings_e_max = e_process(
            warm_rows
                .iter()
                .map(|(cold, warm)| f64::from(*cold) / f64::from(*warm) >= t.savings_floor),
            0.5,
            t.betting_fraction,
        );
        let accept_rejected = accept_e_max >= target;
        let savings_rejected = savings_e_max >= target;
        activated &= accept_rejected && savings_rejected;
        reports.push(StratumReport {
            stratum: spec.name,
            samples: have,
            accepts: u32::try_from(accepts).unwrap_or(u32::MAX),
            warm_samples: u32::try_from(warm_rows.len()).unwrap_or(u32::MAX),
            savings_wins: u32::try_from(savings_wins).unwrap_or(u32::MAX),
            accept_e_max,
            savings_e_max,
            accept_rejected,
            savings_rejected,
        });
    }
    Ok(ActivationReport {
        frame_identity: frame.identity(),
        strata: reports,
        activated,
    })
}
