//! Horizon trigger B (bead `frankensim-epic-addendum-xpck.5.6`): the
//! agent-operator mode and explanation honesty gate for Proposal B
//! ("Explanation objects", [`crate::proposals`]).
//!
//! Rule-4 condition:
//! - Human-driven mode emits [`ExplanationDisposition::Rule4Defer`] (human engineers
//!   inspect raw diagnostics and plots directly; explanation synthesis is not activated).
//! - Agent-operator mode activates ONLY when:
//!   1. Attributed channels + residual reconcile with observed QoI within bounds,
//!      with a failure rate $\le 10\%$.
//!   2. Every high-residual case refuses the honesty gate, and narrative generation
//!      is also strictly refused (no storytelling over unexplained residuals).

use fs_blake3::hash_bytes;

/// Maximum admitted reconciliation failure rate (10%).
pub const MAX_RECONCILIATION_FAILURE_RATE: f64 = 0.10;

/// Default tolerance on attribution reconciliation $|Q - \sum c_i| \le \text{tol}$.
pub const DEFAULT_RECONCILIATION_TOL: f64 = 1e-4;

/// Governed execution operator mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperatorMode {
    /// Human-driven engineer interaction (Rule-4 deferral applies).
    HumanDriven,
    /// Autonomous agent-orchestrated closed-loop design/solve mode.
    AgentOperator,
}

/// One explanation case record in the retained battery.
#[derive(Debug, Clone, PartialEq)]
pub struct ExplanationCase {
    /// Case identifier.
    pub case_id: String,
    /// Observed Quantity of Interest (QoI).
    pub observed_qoi: f64,
    /// Attributed decomposition channels (e.g. skin-friction, induced, wave drag).
    pub attributed_channels: Vec<f64>,
    /// Allowed attribution tolerance.
    pub tolerance: f64,
    /// True if the honesty gate passed; false if refused.
    pub honesty_gate_passed: bool,
    /// True if narrative text was emitted; false if narrative was refused.
    pub narrative_emitted: bool,
}

impl ExplanationCase {
    /// Compute the absolute attribution residual $|Q - \sum c_i|$.
    #[must_use]
    pub fn residual(&self) -> f64 {
        let sum: f64 = self.attributed_channels.iter().sum();
        (self.observed_qoi - sum).abs()
    }

    /// Check if the decomposition reconciles within tolerance.
    #[must_use]
    pub fn is_reconciled(&self) -> bool {
        self.residual() <= self.tolerance
    }
}

/// Typed refusals for malformed explanation batteries.
#[derive(Debug, Clone, PartialEq)]
pub enum TriggerBRefusal {
    /// Empty explanation battery.
    EmptyBattery,
    /// Non-finite QoI or channel value.
    NonFiniteQuantity { case_id: String },
    /// Non-positive or non-finite tolerance.
    InvalidTolerance { case_id: String, tol: f64 },
    /// Narrative emitted despite honesty gate refusal (storytelling violation).
    NarrativeOverRefusal { case_id: String },
    /// Unreconciled residual passed the honesty gate (smear violation).
    UnreconciledPassedGate { case_id: String, residual: f64, tol: f64 },
}

/// Activation verdict for Proposal B.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriggerBVerdict {
    /// Agent-operator mode and reconciliation / honesty pass: activate explanation objects.
    Activate,
    /// Human-driven operator mode: deferred under Rule 4.
    Rule4Defer,
    /// Agent-operator mode but reconciliation or honesty failed: deferred.
    Defer,
}

/// Overall population disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExplanationDisposition {
    /// Promoted and activated.
    Activate,
    /// Deferred under Rule-4 human-driven posture.
    Rule4Defer,
    /// Evaluated and deferred due to failure rate or honesty violations.
    Defer,
    /// No retained explanation battery present.
    NoData,
}

/// Immutable receipt of a Trigger B evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct TriggerBReceipt {
    pub proposal: &'static str,
    pub disposition: ExplanationDisposition,
    pub verdict: TriggerBVerdict,
    pub operator_mode: OperatorMode,
    pub failure_rate: f64,
    pub receipt_hash: String,
    pub reason: String,
}

/// Evaluate Proposal B activation conditions.
///
/// # Errors
/// Returns [`TriggerBRefusal`] if any case in the battery violates honesty invariants.
pub fn evaluate_trigger_b(
    mode: OperatorMode,
    battery: &[ExplanationCase],
) -> Result<TriggerBVerdict, TriggerBRefusal> {
    if mode == OperatorMode::HumanDriven {
        return Ok(TriggerBVerdict::Rule4Defer);
    }

    if battery.is_empty() {
        return Err(TriggerBRefusal::EmptyBattery);
    }

    let mut failures = 0usize;
    for case in battery {
        if !case.observed_qoi.is_finite() || !case.tolerance.is_finite() || case.tolerance <= 0.0 {
            return Err(TriggerBRefusal::InvalidTolerance {
                case_id: case.case_id.clone(),
                tol: case.tolerance,
            });
        }
        for &ch in &case.attributed_channels {
            if !ch.is_finite() {
                return Err(TriggerBRefusal::NonFiniteQuantity { case_id: case.case_id.clone() });
            }
        }

        let reconciled = case.is_reconciled();
        if !reconciled {
            failures += 1;
            // High residual MUST NOT pass the honesty gate
            if case.honesty_gate_passed {
                return Err(TriggerBRefusal::UnreconciledPassedGate {
                    case_id: case.case_id.clone(),
                    residual: case.residual(),
                    tol: case.tolerance,
                });
            }
            // Narrative MUST be refused when honesty gate refuses
            if case.narrative_emitted {
                return Err(TriggerBRefusal::NarrativeOverRefusal {
                    case_id: case.case_id.clone(),
                });
            }
        }
    }

    let failure_rate = failures as f64 / battery.len() as f64;
    if failure_rate <= MAX_RECONCILIATION_FAILURE_RATE {
        Ok(TriggerBVerdict::Activate)
    } else {
        Ok(TriggerBVerdict::Defer)
    }
}

/// Mint an immutable decision receipt for Proposal B.
#[must_use]
pub fn mint_trigger_b_receipt(
    mode: OperatorMode,
    battery_opt: Option<&[ExplanationCase]>,
) -> TriggerBReceipt {
    if mode == OperatorMode::HumanDriven {
        let hash = hash_bytes(b"org.frankensim.horizon-trigger-b.rule4defer.v1").to_hex();
        return TriggerBReceipt {
            proposal: "B",
            disposition: ExplanationDisposition::Rule4Defer,
            verdict: TriggerBVerdict::Rule4Defer,
            operator_mode: mode,
            failure_rate: 0.0,
            receipt_hash: hash,
            reason: "Rule 4: human-driven operator mode defers explanation-object activation".into(),
        };
    }

    let Some(battery) = battery_opt else {
        let hash = hash_bytes(b"org.frankensim.horizon-trigger-b.nodata.v1").to_hex();
        return TriggerBReceipt {
            proposal: "B",
            disposition: ExplanationDisposition::NoData,
            verdict: TriggerBVerdict::Defer,
            operator_mode: mode,
            failure_rate: f64::NAN,
            receipt_hash: hash,
            reason: "no retained explanation battery evaluated in agent-operator mode yet".into(),
        };
    };

    match evaluate_trigger_b(mode, battery) {
        Ok(TriggerBVerdict::Activate) => {
            let failures = battery.iter().filter(|c| !c.is_reconciled()).count();
            let rate = failures as f64 / battery.len() as f64;
            let mut payload = Vec::new();
            payload.extend_from_slice(b"org.frankensim.horizon-trigger-b.activate.v1");
            payload.extend_from_slice(rate.to_le_bytes().as_slice());
            let hash = hash_bytes(&payload).to_hex();
            TriggerBReceipt {
                proposal: "B",
                disposition: ExplanationDisposition::Activate,
                verdict: TriggerBVerdict::Activate,
                operator_mode: mode,
                failure_rate: rate,
                receipt_hash: hash,
                reason: format!("agent-operator mode active with reconciliation failure rate ({:.1}%, {}/{}) <= 10% and honesty gate respected", rate * 100.0, failures, battery.len()),
            }
        }
        Ok(TriggerBVerdict::Defer) => {
            let failures = battery.iter().filter(|c| !c.is_reconciled()).count();
            let rate = failures as f64 / battery.len() as f64;
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-b.defer.v1").to_hex();
            TriggerBReceipt {
                proposal: "B",
                disposition: ExplanationDisposition::Defer,
                verdict: TriggerBVerdict::Defer,
                operator_mode: mode,
                failure_rate: rate,
                receipt_hash: hash,
                reason: format!("reconciliation failure rate ({:.1}%, {}/{}) exceeds 10% ceiling; deferring activation", rate * 100.0, failures, battery.len()),
            }
        }
        Ok(TriggerBVerdict::Rule4Defer) => {
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-b.rule4defer.v1").to_hex();
            TriggerBReceipt {
                proposal: "B",
                disposition: ExplanationDisposition::Rule4Defer,
                verdict: TriggerBVerdict::Rule4Defer,
                operator_mode: mode,
                failure_rate: 0.0,
                receipt_hash: hash,
                reason: "Rule 4: human-driven operator mode defers explanation-object activation".into(),
            }
        }
        Err(refusal) => {
            let hash = hash_bytes(format!("{:?}", refusal).as_bytes()).to_hex();
            TriggerBReceipt {
                proposal: "B",
                disposition: ExplanationDisposition::Defer,
                verdict: TriggerBVerdict::Defer,
                operator_mode: mode,
                failure_rate: f64::NAN,
                receipt_hash: hash,
                reason: format!("inadmissible explanation battery ({refusal:?}); deferring activation"),
            }
        }
    }
}
