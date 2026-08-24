//! Horizon trigger D (bead `frankensim-epic-addendum-xpck.5.7`): the
//! agent-operator mode and endpoint-targeting evidence gate for Proposal D
//! ("Goodhart guard: optimizer-endpoint escalation", [`crate::proposals`]).
//!
//! Rule-4 condition:
//! - Human-driven mode emits [`GoodhartDisposition::Rule4Defer`].
//! - Agent-operator mode activates ONLY when:
//!   1. All four escalation steps (rung $k+1$, cross-rep re-solve,
//!      $\delta$-perturbation, independent estimator) are verified available
//!      in fixed order or carry explicit provisional status.
//!   2. A preregistered endpoint-versus-random study proves that optimizer
//!      endpoints exhibit a statistically distinguishable bug/infeasibility
//!      catch rate ($p_{\text{endpoint}} > p_{\text{random}}$ with $p \le 0.05$).
//!   3. If indistinguishable, budget returns to general falsification and
//!      activation is deferred.

use fs_blake3::hash_bytes;
pub use crate::horizon_explanation::OperatorMode;

/// Required escalation step kind for Proposal D.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EscalationStep {
    /// Rung $k+1$ higher-fidelity re-solve.
    RungKPlus1,
    /// Cross-representation re-solve (e.g. NURBS vs CSG vs SDF).
    CrossRepresentation,
    /// $\delta$-perturbation neighborhood test.
    DeltaPerturbation,
    /// Independent numerical estimator.
    IndependentEstimator,
}

/// Execution / readiness status of an escalation step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepStatus {
    /// Fully implemented and active.
    Available,
    /// Explicit provisional implementation under monitored envelope.
    Provisional,
    /// Step not yet implemented (blocks full activation).
    Unavailable,
}

/// Statistical catch-rate study comparing optimizer endpoints against random sampling.
#[derive(Debug, Clone, PartialEq)]
pub struct EndpointStudy {
    /// Study identifier / preregistration DOI or artifact hash.
    pub preregistration_ref: String,
    /// Total optimizer endpoint evaluations tested.
    pub endpoint_sample_count: usize,
    /// Catches (defects / discrepancies) found at optimizer endpoints.
    pub endpoint_catches: usize,
    /// Total random baseline points tested.
    pub random_sample_count: usize,
    /// Catches found at random baseline points.
    pub random_catches: usize,
    /// One-sided statistical significance level ($p$-value).
    pub p_value: f64,
}

impl EndpointStudy {
    /// Catch rate at optimizer endpoints.
    #[must_use]
    pub fn endpoint_rate(&self) -> f64 {
        if self.endpoint_sample_count == 0 {
            0.0
        } else {
            self.endpoint_catches as f64 / self.endpoint_sample_count as f64
        }
    }

    /// Catch rate at random exploration baseline.
    #[must_use]
    pub fn random_rate(&self) -> f64 {
        if self.random_sample_count == 0 {
            0.0
        } else {
            self.random_catches as f64 / self.random_sample_count as f64
        }
    }

    /// True if endpoint targeting is statistically distinguishable from random exploration.
    #[must_use]
    pub fn is_distinguishable(&self) -> bool {
        self.endpoint_sample_count >= 30
            && self.random_sample_count >= 30
            && self.endpoint_rate() > self.random_rate()
            && self.p_value.is_finite()
            && self.p_value <= 0.05
    }
}

/// Input premises for Proposal D evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct ProposalDPremises {
    pub step_statuses: [(EscalationStep, StepStatus); 4],
    pub study: Option<EndpointStudy>,
}

/// Typed refusals for malformed Proposal D premises.
#[derive(Debug, Clone, PartialEq)]
pub enum TriggerDRefusal {
    /// Missing endpoint study in agent-operator mode.
    MissingStudy,
    /// Inadmissible sample counts (e.g. zero or catches > samples).
    InadmissibleStudyData { endpoint_catches: usize, endpoint_samples: usize },
    /// Non-finite or negative p-value.
    InvalidPValue { p: f64 },
    /// Escalation steps missing or incomplete.
    IncompleteEscalationSteps,
}

/// Activation verdict for Proposal D.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriggerDVerdict {
    /// Agent-operator mode, 4-step escalation available, and endpoint targeting proven: activate.
    Activate,
    /// Human-driven operator mode: deferred under Rule 4.
    Rule4Defer,
    /// Indistinguishable catch rate: return budget to general falsification, defer activation.
    IndistinguishableDefer,
    /// One or more escalation steps unavailable: provisional defer.
    ProvisionalDefer,
}

/// Overall population disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GoodhartDisposition {
    /// Promoted and activated.
    Activate,
    /// Deferred under Rule-4 human-driven posture.
    Rule4Defer,
    /// Evaluated and deferred (budget returned to general falsification).
    Defer,
    /// No endpoint targeting study conducted yet.
    NoData,
}

/// Immutable receipt of a Trigger D evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct TriggerDReceipt {
    pub proposal: &'static str,
    pub disposition: GoodhartDisposition,
    pub verdict: TriggerDVerdict,
    pub operator_mode: OperatorMode,
    pub receipt_hash: String,
    pub reason: String,
}

/// Evaluate Proposal D activation conditions.
///
/// # Errors
/// Returns [`TriggerDRefusal`] if study data or escalation steps are malformed.
pub fn evaluate_trigger_d(
    mode: OperatorMode,
    premises: &ProposalDPremises,
) -> Result<TriggerDVerdict, TriggerDRefusal> {
    if mode == OperatorMode::HumanDriven {
        return Ok(TriggerDVerdict::Rule4Defer);
    }

    // Verify all 4 escalation steps exist
    let mut has_unavailable = false;
    let mut has_provisional = false;
    for (_step, status) in &premises.step_statuses {
        match status {
            StepStatus::Unavailable => has_unavailable = true,
            StepStatus::Provisional => has_provisional = true,
            StepStatus::Available => {}
        }
    }

    if has_unavailable {
        return Ok(TriggerDVerdict::ProvisionalDefer);
    }

    let Some(study) = &premises.study else {
        return Err(TriggerDRefusal::MissingStudy);
    };

    if study.endpoint_catches > study.endpoint_sample_count
        || study.random_catches > study.random_sample_count
    {
        return Err(TriggerDRefusal::InadmissibleStudyData {
            endpoint_catches: study.endpoint_catches,
            endpoint_samples: study.endpoint_sample_count,
        });
    }

    if !study.p_value.is_finite() || study.p_value < 0.0 || study.p_value > 1.0 {
        return Err(TriggerDRefusal::InvalidPValue { p: study.p_value });
    }

    if !study.is_distinguishable() {
        return Ok(TriggerDVerdict::IndistinguishableDefer);
    }

    if has_provisional {
        Ok(TriggerDVerdict::ProvisionalDefer)
    } else {
        Ok(TriggerDVerdict::Activate)
    }
}

/// Mint an immutable decision receipt for Proposal D.
#[must_use]
pub fn mint_trigger_d_receipt(
    mode: OperatorMode,
    premises_opt: Option<&ProposalDPremises>,
) -> TriggerDReceipt {
    if mode == OperatorMode::HumanDriven {
        let hash = hash_bytes(b"org.frankensim.horizon-trigger-d.rule4defer.v1").to_hex();
        return TriggerDReceipt {
            proposal: "D",
            disposition: GoodhartDisposition::Rule4Defer,
            verdict: TriggerDVerdict::Rule4Defer,
            operator_mode: mode,
            receipt_hash: hash,
            reason: "Rule 4: human-driven operator mode defers optimizer-endpoint Goodhart guard activation".into(),
        };
    }

    let Some(premises) = premises_opt else {
        let hash = hash_bytes(b"org.frankensim.horizon-trigger-d.nodata.v1").to_hex();
        return TriggerDReceipt {
            proposal: "D",
            disposition: GoodhartDisposition::NoData,
            verdict: TriggerDVerdict::IndistinguishableDefer,
            operator_mode: mode,
            receipt_hash: hash,
            reason: "no preregistered endpoint-targeting study evaluated in agent-operator mode yet".into(),
        };
    };

    match evaluate_trigger_d(mode, premises) {
        Ok(TriggerDVerdict::Activate) => {
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-d.activate.v1").to_hex();
            TriggerDReceipt {
                proposal: "D",
                disposition: GoodhartDisposition::Activate,
                verdict: TriggerDVerdict::Activate,
                operator_mode: mode,
                receipt_hash: hash,
                reason: "agent-operator mode active with 4-step escalation available and endpoint targeting statistically proven (p <= 0.05)".into(),
            }
        }
        Ok(TriggerDVerdict::IndistinguishableDefer) => {
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-d.defer-indistinguishable.v1").to_hex();
            TriggerDReceipt {
                proposal: "D",
                disposition: GoodhartDisposition::Defer,
                verdict: TriggerDVerdict::IndistinguishableDefer,
                operator_mode: mode,
                receipt_hash: hash,
                reason: "endpoint targeting statistically indistinguishable from random baseline; budget returned to general falsification pool".into(),
            }
        }
        Ok(TriggerDVerdict::ProvisionalDefer) => {
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-d.defer-provisional.v1").to_hex();
            TriggerDReceipt {
                proposal: "D",
                disposition: GoodhartDisposition::Defer,
                verdict: TriggerDVerdict::ProvisionalDefer,
                operator_mode: mode,
                receipt_hash: hash,
                reason: "one or more Goodhart escalation steps are provisional or unavailable; deferring activation".into(),
            }
        }
        Ok(TriggerDVerdict::Rule4Defer) => {
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-d.rule4defer.v1").to_hex();
            TriggerDReceipt {
                proposal: "D",
                disposition: GoodhartDisposition::Rule4Defer,
                verdict: TriggerDVerdict::Rule4Defer,
                operator_mode: mode,
                receipt_hash: hash,
                reason: "Rule 4: human-driven operator mode defers Goodhart guard activation".into(),
            }
        }
        Err(refusal) => {
            let hash = hash_bytes(format!("{:?}", refusal).as_bytes()).to_hex();
            TriggerDReceipt {
                proposal: "D",
                disposition: GoodhartDisposition::Defer,
                verdict: TriggerDVerdict::IndistinguishableDefer,
                operator_mode: mode,
                receipt_hash: hash,
                reason: format!("inadmissible Goodhart premises ({refusal:?}); deferring activation"),
            }
        }
    }
}
