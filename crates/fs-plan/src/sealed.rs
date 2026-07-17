//! Sealed cost-model authority (bead 2pmb).
//!
//! [`crate::cost::CostModel`] is deliberately freely constructible —
//! it is provisional MATH. Authority is a different thing: the exact
//! roofline loader spends real validation (receipt scope, provenance
//! edges, dependency digests, build identity) before trusting a tune
//! row, and that work must survive the API boundary instead of being
//! discarded into an indistinguishable struct.
//!
//! [`SealedCostModel`] is the carrier. Its `ExactRooflineReceipt`
//! class is mintable ONLY by [`crate::oracle::cost_model_from_tune`]
//! after full validation; the fields are private, so a caller-fitted
//! model can never impersonate one. Tests and non-authoritative hints
//! use [`SealedCostModel::provisional_unaudited`], which works
//! identically but stamps every prediction — and therefore every
//! downstream admission finding and session receipt — with
//! [`CostEvidenceClass::ProvisionalUnaudited`]. Composition cannot
//! upgrade the class: authority follows the mint, not the arithmetic.

use crate::cost::{CostModel, CostPrediction, CostRefusal};

/// How a sealed model's evidence was established.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CostEvidenceClass {
    /// Minted by the exact roofline loader after full receipt,
    /// provenance, dependency, and build-identity validation, AND
    /// still inside the caller's freshness contract at assessment.
    ExactRooflineReceipt,
    /// Explicitly provisional: caller-fitted observations with NO
    /// receipt validation. Useful for tests and non-authoritative
    /// hints; every consumer sees this label and must surface it.
    ProvisionalUnaudited,
    /// Receipt-backed evidence whose FRESHNESS contract is violated at
    /// assessment time (bead jle3m): aged past the caller's horizon,
    /// recorded on a different machine fingerprint, or carrying a
    /// future timestamp. Ranks BELOW a fresh provisional fit: a
    /// violated contract is a red flag, not neutral seniority — a
    /// year-old sealed roofline must never outrank a fresh fit
    /// forever.
    StaleRooflineReceipt,
}

impl CostEvidenceClass {
    /// Stable name for findings and receipts.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            CostEvidenceClass::ExactRooflineReceipt => "exact-roofline-receipt",
            CostEvidenceClass::ProvisionalUnaudited => "provisional-unaudited",
            CostEvidenceClass::StaleRooflineReceipt => "stale-roofline-receipt",
        }
    }

    /// Total evidence lattice (bead jle3m):
    /// fresh exact (2) > fresh provisional (1) > stale exact (0).
    /// Weakest-wins folds compare ranks; mixing can only degrade.
    #[must_use]
    pub const fn rank(self) -> u8 {
        match self {
            CostEvidenceClass::ExactRooflineReceipt => 2,
            CostEvidenceClass::ProvisionalUnaudited => 1,
            CostEvidenceClass::StaleRooflineReceipt => 0,
        }
    }

    /// The weaker of two classes under [`Self::rank`].
    #[must_use]
    pub const fn weakest(self, other: CostEvidenceClass) -> CostEvidenceClass {
        if self.rank() <= other.rank() {
            self
        } else {
            other
        }
    }
}

/// Preregistered freshness contract for assessing sealed evidence
/// (bead jle3m). Pure data: the caller supplies the observation time
/// and current machine fingerprint at assessment — this crate never
/// reads a clock or the host identity itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FreshnessPolicy {
    horizon_ns: i64,
}

impl FreshnessPolicy {
    /// A policy that ages receipts out after `horizon_ns`.
    ///
    /// # Errors
    /// Refuses a non-positive horizon (a zero horizon would stale
    /// every receipt at its own recording instant).
    pub const fn new(horizon_ns: i64) -> Result<FreshnessPolicy, CostRefusal> {
        if horizon_ns <= 0 {
            return Err(CostRefusal::InvalidFreshnessHorizon { horizon_ns });
        }
        Ok(FreshnessPolicy { horizon_ns })
    }

    /// The admitted horizon in ledger nanoseconds.
    #[must_use]
    pub const fn horizon_ns(self) -> i64 {
        self.horizon_ns
    }
}

/// Why an assessment kept or degraded a sealed model's class
/// (bead jle3m). Retained alongside the assessed class so consumers
/// and routers can distinguish age from drift without re-deriving.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StalenessVerdict {
    /// Inside the freshness contract on the same machine.
    Fresh,
    /// Not receipt-backed; freshness does not apply.
    NotApplicable,
    /// Older than the policy horizon at the supplied observation time.
    AgedOut {
        /// Observed age (observation time minus recording time).
        age_ns: i64,
        /// The policy horizon it exceeded.
        horizon_ns: i64,
    },
    /// Recorded under a different machine fingerprint than the one
    /// supplied at assessment.
    MachineDrift,
    /// Recorded in the future of the supplied observation time —
    /// clock skew or forged provenance; fail closed as stale.
    FutureRecording {
        /// How far in the future the recording claims to be.
        ahead_ns: i64,
    },
}

/// The validated scope a sealed model speaks for. For provisional
/// models every receipt-derived field is the literal string
/// `"provisional"` / zero — visibly not a receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CostModelScope {
    kernel: String,
    shape_class: String,
    machine: Vec<u8>,
    run_receipt: String,
    op: u64,
    build_identity: String,
    recorded_at_ns: i64,
}

impl CostModelScope {
    pub(crate) fn from_validated(
        kernel: String,
        shape_class: String,
        machine: Vec<u8>,
        run_receipt: String,
        op: u64,
        build_identity: String,
        recorded_at_ns: i64,
    ) -> CostModelScope {
        CostModelScope {
            kernel,
            shape_class,
            machine,
            run_receipt,
            op,
            build_identity,
            recorded_at_ns,
        }
    }

    fn provisional(label: &str) -> CostModelScope {
        CostModelScope {
            kernel: format!("provisional:{label}"),
            shape_class: "provisional".to_string(),
            machine: Vec::new(),
            run_receipt: "provisional".to_string(),
            op: 0,
            build_identity: "provisional".to_string(),
            recorded_at_ns: 0,
        }
    }

    /// The exact tune kernel (or `provisional:<label>`).
    #[must_use]
    pub fn kernel(&self) -> &str {
        &self.kernel
    }

    /// The exact tune shape class.
    #[must_use]
    pub fn shape_class(&self) -> &str {
        &self.shape_class
    }

    /// The exact roofline machine key bytes (empty for provisional).
    #[must_use]
    pub fn machine(&self) -> &[u8] {
        &self.machine
    }

    /// The finalized-run receipt digest the row was bound to.
    #[must_use]
    pub fn run_receipt(&self) -> &str {
        &self.run_receipt
    }

    /// The ledger operation id the evidence was recorded under.
    #[must_use]
    pub fn op(&self) -> u64 {
        self.op
    }

    /// The validated build identity.
    #[must_use]
    pub fn build_identity(&self) -> &str {
        &self.build_identity
    }

    /// When the evidence operation finished (ledger nanoseconds).
    #[must_use]
    pub fn recorded_at_ns(&self) -> i64 {
        self.recorded_at_ns
    }
}

/// A cost model whose provenance survived the loader boundary.
/// Private fields: the `ExactRooflineReceipt` class cannot be forged
/// outside this crate.
#[derive(Debug, Clone)]
pub struct SealedCostModel {
    model: CostModel,
    scope: CostModelScope,
    class: CostEvidenceClass,
}

impl SealedCostModel {
    pub(crate) fn mint_exact(model: CostModel, scope: CostModelScope) -> SealedCostModel {
        SealedCostModel {
            model,
            scope,
            class: CostEvidenceClass::ExactRooflineReceipt,
        }
    }

    /// Wrap a caller-fitted model as EXPLICITLY provisional evidence.
    /// The label lands in the scope's kernel field and the class lands
    /// in every prediction; nothing downstream can mistake this for
    /// receipt-backed authority.
    #[must_use]
    pub fn provisional_unaudited(model: CostModel, label: &str) -> SealedCostModel {
        SealedCostModel {
            model,
            scope: CostModelScope::provisional(label),
            class: CostEvidenceClass::ProvisionalUnaudited,
        }
    }

    /// Observation count behind the fit.
    #[must_use]
    pub fn n_obs(&self) -> usize {
        self.model.n_obs()
    }

    /// The validated (or visibly provisional) scope.
    #[must_use]
    pub fn scope(&self) -> &CostModelScope {
        &self.scope
    }

    /// The evidence class stamped into every prediction.
    #[must_use]
    pub fn evidence_class(&self) -> CostEvidenceClass {
        self.class
    }

    /// Assess the minted class against the caller's freshness
    /// contract (bead jle3m): receipt-backed evidence degrades to
    /// [`CostEvidenceClass::StaleRooflineReceipt`] when aged past the
    /// horizon, recorded under a different machine fingerprint, or
    /// carrying a future timestamp. Provisional evidence is untouched
    /// (freshness cannot upgrade anything). Pure: `now_ns` and
    /// `current_machine` are supplied by the caller.
    #[must_use]
    pub fn assess(
        &self,
        now_ns: i64,
        current_machine: &[u8],
        policy: FreshnessPolicy,
    ) -> (CostEvidenceClass, StalenessVerdict) {
        if self.class != CostEvidenceClass::ExactRooflineReceipt {
            return (self.class, StalenessVerdict::NotApplicable);
        }
        if self.scope.machine != current_machine {
            return (
                CostEvidenceClass::StaleRooflineReceipt,
                StalenessVerdict::MachineDrift,
            );
        }
        let recorded = self.scope.recorded_at_ns;
        if recorded > now_ns {
            return (
                CostEvidenceClass::StaleRooflineReceipt,
                StalenessVerdict::FutureRecording {
                    ahead_ns: recorded - now_ns,
                },
            );
        }
        let age_ns = now_ns - recorded;
        if age_ns > policy.horizon_ns() {
            return (
                CostEvidenceClass::StaleRooflineReceipt,
                StalenessVerdict::AgedOut {
                    age_ns,
                    horizon_ns: policy.horizon_ns(),
                },
            );
        }
        (
            CostEvidenceClass::ExactRooflineReceipt,
            StalenessVerdict::Fresh,
        )
    }

    /// Predict wall cost at `size` with the class ASSESSED against the
    /// freshness contract instead of the mint-time stamp.
    ///
    /// # Errors
    /// Exactly [`CostModel::predict`]'s refusals.
    pub fn predict_assessed(
        &self,
        size: f64,
        now_ns: i64,
        current_machine: &[u8],
        policy: FreshnessPolicy,
    ) -> Result<(SealedCostPrediction, StalenessVerdict), CostRefusal> {
        let (class, verdict) = self.assess(now_ns, current_machine, policy);
        let mut sealed = self.predict(size)?;
        sealed.evidence = class;
        Ok((sealed, verdict))
    }

    /// Predict wall cost at `size`, carrying scope and class.
    ///
    /// # Errors
    /// Exactly [`CostModel::predict`]'s refusals.
    pub fn predict(&self, size: f64) -> Result<SealedCostPrediction, CostRefusal> {
        let prediction = self.model.predict(size)?;
        Ok(SealedCostPrediction {
            prediction,
            scope: self.scope.clone(),
            evidence: self.class,
        })
    }
}

/// A prediction that remembers exactly what it speaks for.
///
/// Private fields keep the authority carrier opaque: callers may read the
/// numeric bands, complete validated scope, and evidence class, but cannot
/// forge receipt-backed provenance around an arbitrary [`CostPrediction`].
#[derive(Debug, Clone, PartialEq)]
pub struct SealedCostPrediction {
    prediction: CostPrediction,
    scope: CostModelScope,
    evidence: CostEvidenceClass,
}

impl SealedCostPrediction {
    /// The quantile bands (P10/P50/P90, observation count, and extrapolation
    /// flag). Returning the `Copy` value cannot alter the sealed carrier.
    #[must_use]
    pub const fn prediction(&self) -> CostPrediction {
        self.prediction
    }

    /// The complete validated (or visibly provisional) scope.
    #[must_use]
    pub const fn scope(&self) -> &CostModelScope {
        &self.scope
    }

    /// The mint-time or freshness-assessed evidence class.
    #[must_use]
    pub const fn evidence_class(&self) -> CostEvidenceClass {
        self.evidence
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::CostObservation;

    fn fitted() -> CostModel {
        let obs: Vec<CostObservation> = (1..=6)
            .map(|k| CostObservation {
                size: f64::from(k) * 256.0,
                cost_s: 0.05 * f64::from(k) * 256.0,
            })
            .collect();
        CostModel::fit(&obs).expect("fits")
    }

    fn exact_sealed(recorded_at_ns: i64, machine: &[u8]) -> SealedCostModel {
        SealedCostModel::mint_exact(
            fitted(),
            CostModelScope::from_validated(
                "gemm-f64".to_string(),
                "square".to_string(),
                machine.to_vec(),
                "run-receipt-digest".to_string(),
                7,
                "build-identity".to_string(),
                recorded_at_ns,
            ),
        )
    }

    #[test]
    fn freshness_policy_refuses_nonpositive_horizons() {
        assert!(matches!(
            FreshnessPolicy::new(0),
            Err(CostRefusal::InvalidFreshnessHorizon { horizon_ns: 0 })
        ));
        assert!(matches!(
            FreshnessPolicy::new(-1),
            Err(CostRefusal::InvalidFreshnessHorizon { horizon_ns: -1 })
        ));
        assert_eq!(
            FreshnessPolicy::new(1)
                .expect("one nanosecond")
                .horizon_ns(),
            1
        );
    }

    #[test]
    fn assessment_keeps_fresh_receipts_exact_at_the_boundary() {
        let sealed = exact_sealed(1_000, b"machine-a");
        let policy = FreshnessPolicy::new(500).expect("policy");
        // Exactly at the horizon is still fresh (age == horizon).
        assert_eq!(
            sealed.assess(1_500, b"machine-a", policy),
            (
                CostEvidenceClass::ExactRooflineReceipt,
                StalenessVerdict::Fresh
            )
        );
        // One past the horizon ages out with exact accounting.
        assert_eq!(
            sealed.assess(1_501, b"machine-a", policy),
            (
                CostEvidenceClass::StaleRooflineReceipt,
                StalenessVerdict::AgedOut {
                    age_ns: 501,
                    horizon_ns: 500
                }
            )
        );
    }

    #[test]
    fn machine_drift_and_future_recordings_fail_closed() {
        let sealed = exact_sealed(1_000, b"machine-a");
        let policy = FreshnessPolicy::new(1_000_000).expect("policy");
        assert_eq!(
            sealed.assess(1_100, b"machine-b", policy),
            (
                CostEvidenceClass::StaleRooflineReceipt,
                StalenessVerdict::MachineDrift
            )
        );
        assert_eq!(
            sealed.assess(900, b"machine-a", policy),
            (
                CostEvidenceClass::StaleRooflineReceipt,
                StalenessVerdict::FutureRecording { ahead_ns: 100 }
            )
        );
    }

    #[test]
    fn provisional_evidence_is_untouched_by_assessment() {
        let sealed = SealedCostModel::provisional_unaudited(fitted(), "unit-test");
        let policy = FreshnessPolicy::new(1).expect("policy");
        assert_eq!(
            sealed.assess(i64::MAX, b"any-machine", policy),
            (
                CostEvidenceClass::ProvisionalUnaudited,
                StalenessVerdict::NotApplicable
            )
        );
    }

    #[test]
    fn assessed_predictions_carry_the_degraded_class() {
        let sealed = exact_sealed(1_000, b"machine-a");
        let policy = FreshnessPolicy::new(500).expect("policy");
        let (fresh, verdict) = sealed
            .predict_assessed(512.0, 1_200, b"machine-a", policy)
            .expect("predicts");
        assert_eq!(
            fresh.evidence_class(),
            CostEvidenceClass::ExactRooflineReceipt
        );
        assert_eq!(verdict, StalenessVerdict::Fresh);
        let (stale, verdict) = sealed
            .predict_assessed(512.0, 9_000, b"machine-a", policy)
            .expect("predicts");
        assert_eq!(
            stale.evidence_class(),
            CostEvidenceClass::StaleRooflineReceipt
        );
        assert!(matches!(verdict, StalenessVerdict::AgedOut { .. }));
        // The numeric bands are identical; only the evidence label
        // degrades — staleness never rewrites the science.
        assert_eq!(fresh.prediction(), stale.prediction());
        assert_eq!(fresh.scope(), stale.scope());
        // The mint-time stamp itself is immutable.
        assert_eq!(
            sealed.evidence_class(),
            CostEvidenceClass::ExactRooflineReceipt
        );
    }

    #[test]
    fn the_evidence_lattice_never_upgrades_by_mixing() {
        use CostEvidenceClass::{
            ExactRooflineReceipt as Exact, ProvisionalUnaudited as Provisional,
            StaleRooflineReceipt as Stale,
        };
        assert!(Exact.rank() > Provisional.rank());
        assert!(Provisional.rank() > Stale.rank());
        assert_eq!(Exact.weakest(Provisional), Provisional);
        assert_eq!(Provisional.weakest(Stale), Stale);
        assert_eq!(Exact.weakest(Stale), Stale);
        assert_eq!(Exact.weakest(Exact), Exact);
        // A year-old sealed roofline no longer outranks a fresh
        // provisional fit: the stale class loses the fold.
        assert_eq!(Stale.weakest(Provisional), Stale);
    }

    #[test]
    fn provisional_mint_labels_everything_and_upgrades_nothing() {
        let sealed = SealedCostModel::provisional_unaudited(fitted(), "unit-test");
        assert_eq!(
            sealed.evidence_class(),
            CostEvidenceClass::ProvisionalUnaudited
        );
        assert_eq!(sealed.scope().kernel(), "provisional:unit-test");
        assert_eq!(sealed.scope().run_receipt(), "provisional");
        assert_eq!(sealed.scope().op(), 0);
        assert!(sealed.scope().machine().is_empty());
        let prediction = sealed.predict(512.0).expect("predicts");
        assert_eq!(
            prediction.evidence_class(),
            CostEvidenceClass::ProvisionalUnaudited,
            "the class travels into every prediction"
        );
        assert_eq!(prediction.scope(), sealed.scope());
        assert_eq!(prediction.scope().kernel(), "provisional:unit-test");
        assert_eq!(prediction.scope().recorded_at_ns(), 0);
        // The math is untouched by the seal: bands match the raw model.
        let raw = fitted().predict(512.0).expect("raw predicts");
        assert_eq!(prediction.prediction(), raw);
    }

    #[test]
    fn exact_mint_carries_the_full_validated_scope() {
        let scope = CostModelScope::from_validated(
            "simd-axpy-f64".to_string(),
            "roofline-v1:run=abc:op=41".to_string(),
            vec![7u8; 40],
            "abc".to_string(),
            41,
            "build-xyz".to_string(),
            1_784_000_000_000_000_000,
        );
        let sealed = SealedCostModel::mint_exact(fitted(), scope);
        assert_eq!(
            sealed.evidence_class(),
            CostEvidenceClass::ExactRooflineReceipt
        );
        let prediction = sealed.predict(512.0).expect("predicts");
        assert_eq!(
            prediction.evidence_class(),
            CostEvidenceClass::ExactRooflineReceipt
        );
        assert_eq!(prediction.scope(), sealed.scope());
        assert_eq!(prediction.scope().kernel(), "simd-axpy-f64");
        assert_eq!(
            prediction.scope().shape_class(),
            "roofline-v1:run=abc:op=41"
        );
        assert_eq!(prediction.scope().run_receipt(), "abc");
        assert_eq!(prediction.scope().op(), 41);
        assert_eq!(prediction.scope().machine(), &[7u8; 40][..]);
        assert_eq!(prediction.scope().build_identity(), "build-xyz");
        assert_eq!(
            prediction.scope().recorded_at_ns(),
            1_784_000_000_000_000_000
        );
    }

    #[test]
    fn refusals_pass_through_the_seal() {
        let sealed = SealedCostModel::provisional_unaudited(CostModel::new(), "empty");
        assert!(
            sealed.predict(512.0).is_err(),
            "an empty model refuses through the seal exactly as raw"
        );
        assert_eq!(sealed.n_obs(), 0);
    }

    #[test]
    fn class_names_are_stable_receipt_vocabulary() {
        assert_eq!(
            CostEvidenceClass::ExactRooflineReceipt.name(),
            "exact-roofline-receipt"
        );
        assert_eq!(
            CostEvidenceClass::ProvisionalUnaudited.name(),
            "provisional-unaudited"
        );
    }
}
