//! fs-evidence — `Evidence<T>` / `Certified<T>`: certificates that travel
//! INSIDE values (patch Rev B; plan §11.4's Error Ledger substrate).
//! Layer: UTIL (usable by every layer without violations).
//!
//! The core argument: without model-form evidence, FrankenSim can produce
//! beautifully certified WRONG answers — mesh error 0.7%, residual
//! irrelevant, but turbulence-closure discrepancy 8–15%, so the design
//! ranking is not decision-grade. Being able to SAY that sentence is the
//! difference between a solver and an engineering-intelligence system.
//!
//! An [`Evidence<T>`] carries a value plus FOUR uncertainty slices:
//! - [`NumericalCertificate`] — enclosure/estimate of the scalar quantity
//!   of interest (the plan Appendix B `Certified<T>` fields: value +
//!   interval bound + provenance + adjoint hook — kept intact as the
//!   numerical slice);
//! - [`StatisticalCertificate`] — e-values / confidence half-widths;
//! - [`ModelEvidence`] — which model cards, whose assumptions, what
//!   validity domain, how big the model-form discrepancy band is;
//! - [`SensitivitySummary`] — d(qoi)/d(param) headline entries.
//!
//! Composition ([`Evidence::combine`]) is CONSERVATIVE by construction:
//! numerical enclosures compose with outward rounding, validity domains
//! intersect, assumptions union, discrepancy bands add (first-order),
//! provenance hashes chain — the G0 conservativeness laws in the
//! conformance suite are the contract.
//!
//! Decision honesty: [`Evidence::assess`] splits the total band by source
//! and answers "is this decision-grade at threshold θ, and if not, WHICH
//! uncertainty dominates" — the input to decision-aware escalation
//! ([`EscalationAdvice`]): refine numerics, gather samples, or escalate
//! MODEL fidelity when model-form dominates (the case cheap refinement
//! cannot fix).

use core::fmt;
use std::collections::BTreeMap;
use std::fmt::Write as _;

pub mod action;
pub mod admitted;
pub mod balance;
mod cards;
pub mod color;
mod discrepancy;
pub mod falsify;
pub mod vv;

pub use admitted::{
    AdmissionDecision, AdmissionReceipt, AdmissionRejection, AdmissionVerifier, AdmittedColor,
    NoAdmissionVerifier, no_admission_policy,
};
pub use cards::{Ambition, ModelCard, ModelRegistry, RegistryError};
pub use color::{
    COLOR_ALGEBRA_VERSION, Color, ColorError, ColorPayloadError, ColorRank, Demotion, IntervalOp,
    MAX_COLOR_IDENTITY_BYTES, check_regime, color_identity_reason, color_leaf_identity_reason,
    color_of, compose, demotion_estimator_identity, intersect_domains, regime_demotion,
    validate_color_payload, verified_from,
};
pub use discrepancy::{
    DiscrepancyBand, DiscrepancyModel, FidelityPair, FitError, ModelBracket, OutOfDomain,
};
pub use falsify::{
    AttemptRecord, ClaimContext, EstimatorBugCandidate, FalsifierAttempt, FalsifierHistory,
    FalsifierOutcome, FalsifierRegistry, FalsifierSpec, FalsifyError, TombstoneCandidate,
    allocate_budget,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Content-address of an artifact/operation (FNV-1a until the BLAKE3-class
/// ledger hash supersedes it — same upgrade path as fs-obs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ProvenanceHash(pub u64);

impl ProvenanceHash {
    /// Hash raw bytes into a provenance id.
    #[must_use]
    pub fn of_bytes(bytes: &[u8]) -> Self {
        ProvenanceHash(fs_obs::fnv1a64(bytes))
    }

    /// Chain an operation over operand provenances (deterministic,
    /// order-sensitive — `combine("sub", [a, b])` differs from `[b, a]`).
    ///
    /// Encoded with the canonical replay-identity schema (gp3.14):
    /// typed length-prefixed fields, NOT `op|hash|hash` concatenation —
    /// an op name containing `|` or a hex-looking suffix can no longer
    /// imitate a different chain.
    #[must_use]
    pub fn chain(op: &str, operands: &[ProvenanceHash]) -> Self {
        let mut b = fs_obs::ident::IdentityBuilder::new("provenance-chain").str("op", op);
        for p in operands {
            b = b.child_root64("operand", p.0);
        }
        ProvenanceHash(b.finish().root())
    }
}

/// How strong the numerical bound is (severity-ordered: composing takes
/// the weakest, and float composition never claims `Exact`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum NumericalKind {
    /// The bound is exact (integer-representable results, by-construction
    /// identities).
    Exact,
    /// A rigorous outward-rounded enclosure (fs-ivl class).
    Enclosure,
    /// A reported band without rigor (a-posteriori estimators).
    Estimate,
    /// No numerical claim: the band is infinite and everything downstream
    /// inherits the refusal.
    NoClaim,
}

impl NumericalKind {
    fn name(self) -> &'static str {
        match self {
            NumericalKind::Exact => "exact",
            NumericalKind::Enclosure => "enclosure",
            NumericalKind::Estimate => "estimate",
            NumericalKind::NoClaim => "no-claim",
        }
    }
}

/// Canonical `(lo ≤ hi)` ordering that PROPAGATES NaN. `f64::min`/`max`
/// silently DISCARD a NaN operand — normalizing `(NaN, 1.0)` to `(1.0, 1.0)`
/// would mint razor-thin false precision from a garbage bound (bead wa8i E4);
/// a NaN input yields a NaN interval that fails closed at the color gate.
fn ordered_bounds(lo: f64, hi: f64) -> (f64, f64) {
    if lo.is_nan() || hi.is_nan() {
        (f64::NAN, f64::NAN)
    } else if lo <= hi {
        (lo, hi)
    } else {
        (hi, lo)
    }
}

/// The numerical slice: `[lo, hi]` encloses (or estimates) the scalar QoI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NumericalCertificate {
    /// Bound strength.
    pub kind: NumericalKind,
    /// Lower bound on the QoI.
    pub lo: f64,
    /// Upper bound on the QoI.
    pub hi: f64,
}

impl NumericalCertificate {
    /// An exact value.
    #[must_use]
    pub fn exact(v: f64) -> Self {
        NumericalCertificate {
            kind: NumericalKind::Exact,
            lo: v,
            hi: v,
        }
    }

    /// A rigorous enclosure (callers guarantee `lo <= hi`; out-of-order bounds
    /// are normalized by swapping — a teaching-free total function). A NaN bound
    /// is PROPAGATED, not silently dropped, so it fails closed at the color gate.
    #[must_use]
    pub fn enclosure(lo: f64, hi: f64) -> Self {
        let (lo, hi) = ordered_bounds(lo, hi);
        NumericalCertificate {
            kind: NumericalKind::Enclosure,
            lo,
            hi,
        }
    }

    /// A non-rigorous band.
    #[must_use]
    pub fn estimate(lo: f64, hi: f64) -> Self {
        let (lo, hi) = ordered_bounds(lo, hi);
        NumericalCertificate {
            kind: NumericalKind::Estimate,
            lo,
            hi,
        }
    }

    /// The explicit refusal.
    #[must_use]
    pub fn no_claim() -> Self {
        NumericalCertificate {
            kind: NumericalKind::NoClaim,
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }

    /// Relative half-width against a reference magnitude (∞ for NoClaim).
    #[must_use]
    pub fn rel_half_width(&self, reference: f64) -> f64 {
        if self.kind == NumericalKind::NoClaim
            || !reference.is_finite()
            || self.lo.is_nan()
            || self.hi.is_nan()
            || self.lo > self.hi
        {
            return f64::INFINITY;
        }
        let half = 0.5 * (self.hi - self.lo);
        let relative = half / reference.abs().max(f64::MIN_POSITIVE);
        if half.is_finite() && half >= 0.0 && relative.is_finite() {
            relative
        } else {
            f64::INFINITY
        }
    }

    fn compose(a: &Self, b: &Self, op: Op) -> Self {
        if a.kind == NumericalKind::NoClaim || b.kind == NumericalKind::NoClaim {
            return Self::no_claim();
        }
        // Public fields permit hand-built certificates. Preserve invalidity
        // instead of allowing f64::min/max to discard a NaN operand.
        let kind = a.kind.max(b.kind).max(NumericalKind::Enclosure);
        if a.lo.is_nan() || a.hi.is_nan() || b.lo.is_nan() || b.hi.is_nan() {
            return Self::whole_line(kind);
        }
        if a.lo > a.hi || b.lo > b.hi {
            return Self::whole_line(kind);
        }
        let (lo, hi) = match op {
            Op::Add => (a.lo + b.lo, a.hi + b.hi),
            Op::Sub => (a.lo - b.hi, a.hi - b.lo),
            Op::Mul => {
                let c = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi];
                if c.iter().any(|v| v.is_nan()) {
                    return Self::whole_line(kind);
                }
                (
                    c.iter().copied().fold(f64::INFINITY, f64::min),
                    c.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                )
            }
            Op::Min => (a.lo.min(b.lo), a.hi.min(b.hi)),
            Op::Max => (a.lo.max(b.lo), a.hi.max(b.hi)),
        };
        // Outward rounding: one ulp each way covers the composition
        // arithmetic itself. Float composition never claims Exact.
        if lo.is_nan() || hi.is_nan() || lo > hi {
            Self::whole_line(kind)
        } else {
            NumericalCertificate {
                kind,
                lo: lo.next_down(),
                hi: hi.next_up(),
            }
        }
    }

    fn whole_line(kind: NumericalKind) -> Self {
        NumericalCertificate {
            kind,
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }
}

/// The statistical slice. Composition is CONSERVATIVE-WEAKEST v1 (proper
/// e-value arithmetic arrives with fs-eproc — CONTRACT no-claims).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StatisticalCertificate {
    /// No stochastic component.
    None,
    /// An anytime-valid e-value against the stated null level.
    EValue {
        /// The finite, non-negative e-value.
        e: f64,
        /// The finite design level in `(0, 1)`.
        alpha: f64,
    },
    /// A confidence half-width around the QoI.
    HalfWidth {
        /// Finite, non-negative absolute half-width.
        half_width: f64,
        /// Finite confidence level in `(0, 1)`.
        confidence: f64,
    },
}

impl StatisticalCertificate {
    /// Relative statistical band against a reference magnitude. E-values
    /// contribute zero width (the decision they certify is already
    /// anytime-valid); half-widths scale by the reference.
    #[must_use]
    pub fn rel_width(&self, reference: f64) -> f64 {
        if self.validation_issue().is_some() {
            return f64::INFINITY;
        }
        match self {
            StatisticalCertificate::None | StatisticalCertificate::EValue { .. } => 0.0,
            StatisticalCertificate::HalfWidth { half_width, .. } => {
                if reference.is_finite() {
                    half_width / reference.abs().max(f64::MIN_POSITIVE)
                } else {
                    f64::INFINITY
                }
            }
        }
    }

    fn compose(a: &Self, b: &Self) -> Self {
        use StatisticalCertificate as S;
        if a.validation_issue().is_some() || b.validation_issue().is_some() {
            return S::HalfWidth {
                half_width: f64::INFINITY,
                confidence: 0.0,
            };
        }
        match (*a, *b) {
            (S::None, x) | (x, S::None) => x,
            (
                S::HalfWidth {
                    half_width: w1,
                    confidence: c1,
                },
                S::HalfWidth {
                    half_width: w2,
                    confidence: c2,
                },
            ) => S::HalfWidth {
                half_width: w1 + w2,
                confidence: c1.min(c2),
            },
            (S::EValue { e: e1, alpha: a1 }, S::EValue { e: e2, alpha: a2 }) => S::EValue {
                e: e1.min(e2),
                alpha: a1.max(a2),
            },
            // Mixed kinds: keep the width-bearing certificate (the weaker
            // decision story) — conservative-weakest v1.
            (w @ S::HalfWidth { .. }, S::EValue { .. })
            | (S::EValue { .. }, w @ S::HalfWidth { .. }) => w,
        }
    }

    fn validation_issue(&self) -> Option<(&'static str, f64, &'static str)> {
        match *self {
            StatisticalCertificate::None => None,
            StatisticalCertificate::EValue { e, .. } if !e.is_finite() || e < 0.0 => {
                Some(("e", e, "must be finite and non-negative"))
            }
            StatisticalCertificate::EValue { alpha, .. }
                if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 =>
            {
                Some((
                    "alpha",
                    alpha,
                    "must be finite and strictly between 0 and 1",
                ))
            }
            StatisticalCertificate::HalfWidth { half_width, .. }
                if !half_width.is_finite() || half_width < 0.0 =>
            {
                Some(("half_width", half_width, "must be finite and non-negative"))
            }
            StatisticalCertificate::HalfWidth { confidence, .. }
                if !confidence.is_finite() || confidence <= 0.0 || confidence >= 1.0 =>
            {
                Some((
                    "confidence",
                    confidence,
                    "must be finite and strictly between 0 and 1",
                ))
            }
            StatisticalCertificate::EValue { .. } | StatisticalCertificate::HalfWidth { .. } => {
                None
            }
        }
    }
}

/// A named-parameter validity box: which region of parameter space the
/// evidence is good for. Missing parameters are unconstrained.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ValidityDomain {
    bounds: BTreeMap<String, (f64, f64)>,
}

impl ValidityDomain {
    /// The unconstrained domain.
    #[must_use]
    pub fn unconstrained() -> Self {
        ValidityDomain::default()
    }

    /// The declared axis bounds (read-only view).
    #[must_use]
    pub fn bounds(&self) -> &BTreeMap<String, (f64, f64)> {
        &self.bounds
    }

    /// Constrain one parameter to `[lo, hi]`. Reversed finite endpoints
    /// normalize; a NaN endpoint is preserved as an unusable domain.
    #[must_use]
    pub fn with(mut self, param: impl Into<String>, lo: f64, hi: f64) -> Self {
        self.bounds.insert(param.into(), ordered_bounds(lo, hi));
        self
    }

    /// Constraint bounds for a parameter, if any.
    #[must_use]
    pub fn bound(&self, param: &str) -> Option<(f64, f64)> {
        self.bounds.get(param).copied()
    }

    /// True when the point satisfies every constraint.
    #[must_use]
    pub fn contains(&self, point: &BTreeMap<String, f64>) -> bool {
        self.bounds.iter().all(|(k, &(lo, hi))| {
            lo.is_finite()
                && hi.is_finite()
                && lo <= hi
                && point
                    .get(k)
                    .is_some_and(|&v| v.is_finite() && v >= lo && v <= hi)
        })
    }

    /// Per-parameter intersection (the composition law: composed validity
    /// is EXACTLY the intersection — never wider).
    #[must_use]
    pub fn intersect(&self, other: &Self) -> Self {
        let mut out = self.bounds.clone();
        for (k, &(lo2, hi2)) in &other.bounds {
            out.entry(k.clone())
                .and_modify(|(lo, hi)| {
                    if lo.is_nan() || hi.is_nan() || lo2.is_nan() || hi2.is_nan() {
                        *lo = f64::NAN;
                        *hi = f64::NAN;
                    } else {
                        *lo = lo.max(lo2);
                        *hi = hi.min(hi2);
                    }
                })
                .or_insert((lo2, hi2));
        }
        ValidityDomain { bounds: out }
    }

    /// True when some parameter's interval is empty or unusable.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bounds
            .values()
            .any(|&(lo, hi)| !lo.is_finite() || !hi.is_finite() || lo > hi)
    }

    /// Constrained parameter names, sorted (BTreeMap order — deterministic
    /// renderings and audits key off this).
    #[must_use]
    pub fn param_names(&self) -> Vec<String> {
        self.bounds.keys().cloned().collect()
    }

    pub(crate) fn to_json(&self) -> String {
        let mut s = String::from("{");
        for (i, (k, (lo, hi))) in self.bounds.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{}:[{},{}]", json_string(k), fmt_f64(*lo), fmt_f64(*hi));
        }
        s.push('}');
        s
    }
}

/// Headline sensitivities d(qoi)/d(param). Merging keeps the larger
/// magnitude per parameter (conservative).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SensitivitySummary {
    /// Parameter → d(qoi)/d(param).
    pub d_qoi: BTreeMap<String, f64>,
}

impl SensitivitySummary {
    fn merge(a: &Self, b: &Self) -> Self {
        let mut out = a.d_qoi.clone();
        for (k, &v) in &b.d_qoi {
            out.entry(k.clone())
                .and_modify(|cur| {
                    if v.abs() > cur.abs() {
                        *cur = v;
                    }
                })
                .or_insert(v);
        }
        SensitivitySummary { d_qoi: out }
    }
}

/// The model-form slice: which cards, whose assumptions, what validity,
/// how big the discrepancy band.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelEvidence {
    /// Names of the model cards in play (sorted, deduplicated).
    pub cards: Vec<String>,
    /// Union of stated assumptions (sorted, deduplicated).
    pub assumptions: Vec<String>,
    /// Validity domain (composition intersects).
    pub validity: ValidityDomain,
    /// Non-negative relative model-form discrepancy band (composition adds —
    /// first-order conservative). Positive infinity is an explicit unbounded
    /// claim; NaN and negative values are invalid.
    pub discrepancy_rel: f64,
    /// False the moment any constituent was queried outside its domain.
    pub in_domain: bool,
}

impl ModelEvidence {
    /// No model-form claims (pure numerics).
    #[must_use]
    pub fn none() -> Self {
        ModelEvidence {
            cards: Vec::new(),
            assumptions: Vec::new(),
            validity: ValidityDomain::unconstrained(),
            discrepancy_rel: 0.0,
            in_domain: true,
        }
    }

    /// Evidence contributed by one model card evaluated at `point`.
    #[must_use]
    pub fn from_card(card: &ModelCard, point: &BTreeMap<String, f64>) -> Self {
        ModelEvidence {
            cards: vec![card.name.clone()],
            assumptions: card.assumptions.clone(),
            validity: card.validity.clone(),
            discrepancy_rel: card.discrepancy_rel,
            in_domain: card.validity.contains(point),
        }
    }

    fn compose(a: &Self, b: &Self) -> Self {
        let mut cards = [a.cards.clone(), b.cards.clone()].concat();
        cards.sort_unstable();
        cards.dedup();
        let mut assumptions = [a.assumptions.clone(), b.assumptions.clone()].concat();
        assumptions.sort_unstable();
        assumptions.dedup();
        ModelEvidence {
            cards,
            assumptions,
            validity: a.validity.intersect(&b.validity),
            discrepancy_rel: if Self::valid_discrepancy(a.discrepancy_rel)
                && Self::valid_discrepancy(b.discrepancy_rel)
            {
                a.discrepancy_rel + b.discrepancy_rel
            } else {
                f64::INFINITY
            },
            in_domain: a.in_domain && b.in_domain,
        }
    }

    fn valid_discrepancy(discrepancy_rel: f64) -> bool {
        !discrepancy_rel.is_nan() && discrepancy_rel >= 0.0
    }
}

/// Which uncertainty source dominates. Declaration order is the
/// deterministic tie-break (model-form first: it is the band cheap
/// refinement cannot shrink, so ties escalate the model).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UncertaintySource {
    /// Model-form discrepancy / out-of-domain use.
    ModelForm,
    /// Statistical half-width.
    Statistical,
    /// Numerical band.
    Numerical,
}

impl UncertaintySource {
    fn name(self) -> &'static str {
        match self {
            UncertaintySource::ModelForm => "model-form",
            UncertaintySource::Statistical => "statistical",
            UncertaintySource::Numerical => "numerical",
        }
    }
}

/// Is the evidence good enough for the pending decision?
#[derive(Debug, Clone, PartialEq)]
pub enum DecisionStatus {
    /// Total relative band within the decision threshold.
    DecisionGrade,
    /// Not decision-grade; says WHICH source dominates and why.
    NotDecisionGrade {
        /// The dominant uncertainty source.
        dominant: UncertaintySource,
        /// A teaching sentence for the report.
        detail: String,
    },
}

/// Per-source relative bands (the assess/escalation currency).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UncertaintyBreakdown {
    /// Numerical relative half-width.
    pub numerical_rel: f64,
    /// Statistical relative width.
    pub statistical_rel: f64,
    /// Model-form relative band (∞ when out of domain).
    pub model_rel: f64,
}

impl UncertaintyBreakdown {
    /// First-order conservative total.
    #[must_use]
    pub fn total_rel(&self) -> f64 {
        Self::usable_band(self.numerical_rel)
            + Self::usable_band(self.statistical_rel)
            + Self::usable_band(self.model_rel)
    }

    /// The dominant source (declaration-order tie-break).
    #[must_use]
    pub fn dominant(&self) -> UncertaintySource {
        let entries = [
            (
                UncertaintySource::ModelForm,
                Self::usable_band(self.model_rel),
            ),
            (
                UncertaintySource::Statistical,
                Self::usable_band(self.statistical_rel),
            ),
            (
                UncertaintySource::Numerical,
                Self::usable_band(self.numerical_rel),
            ),
        ];
        let mut best = entries[0];
        for &(src, band) in &entries[1..] {
            if band > best.1 {
                best = (src, band);
            }
        }
        best.0
    }

    fn usable_band(band: f64) -> f64 {
        if band.is_nan() || band < 0.0 {
            f64::INFINITY
        } else {
            band
        }
    }
}

/// Decision-aware escalation: what to spend money on next (the HELM hook).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EscalationAdvice {
    /// Already decision-grade.
    NoneNeeded,
    /// Numerical band dominates: refine mesh/tolerance.
    RefineNumerics,
    /// Statistical width dominates: gather more samples.
    GatherMoreSamples,
    /// Model-form dominates: cheap refinement CANNOT fix this — escalate
    /// model fidelity (or bracket).
    EscalateModelFidelity,
}

/// QoI composition operators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Sum of QoIs.
    Add,
    /// Difference of QoIs.
    Sub,
    /// Product of QoIs.
    Mul,
    /// Pointwise minimum.
    Min,
    /// Pointwise maximum.
    Max,
}

impl Op {
    fn name(self) -> &'static str {
        match self {
            Op::Add => "add",
            Op::Sub => "sub",
            Op::Mul => "mul",
            Op::Min => "min",
            Op::Max => "max",
        }
    }

    fn apply(self, a: f64, b: f64) -> f64 {
        match self {
            Op::Add => a + b,
            Op::Sub => a - b,
            Op::Mul => a * b,
            Op::Min => a.min(b),
            Op::Max => a.max(b),
        }
    }
}

/// A value with its full evidence: the noun that travels through every
/// layer (plan §1's founding move).
#[derive(Debug, Clone, PartialEq)]
pub struct Evidence<T> {
    /// The carried value.
    pub value: T,
    /// The scalar quantity of interest the certificates describe (equals
    /// `value` for scalar evidence).
    pub qoi: f64,
    /// Numerical slice.
    pub numerical: NumericalCertificate,
    /// Statistical slice.
    pub statistical: StatisticalCertificate,
    /// Model-form slice.
    pub model: ModelEvidence,
    /// Sensitivity headline.
    pub sensitivity: SensitivitySummary,
    /// Content-address of the producing operation.
    pub provenance: ProvenanceHash,
    /// Adjoint hook: content-address of the tape/solver artifact that can
    /// produce d(qoi)/d(inputs) (fs-ad wires this; None = no gradient
    /// support claimed).
    pub adjoint_ref: Option<ProvenanceHash>,
}

/// `Certified<T>`: evidence whose certificate slate PASSED the full
/// validation boundary in [`Evidence::certified`] — rigorous numerics
/// (Exact/Enclosure), finite ordered bounds, QoI containment, valid
/// statistical parameters, and in-domain non-negative model evidence.
///
/// UNFORGEABLE BY CONSTRUCTION (bead gp3.2.1): the inner evidence is
/// private, the ONLY constructor is [`Evidence::certified`], and no
/// mutable access exists — so a certificate cannot be forged by literal
/// construction or weakened by post-certification mutation. Reads go
/// through `Deref<Target = Evidence<T>>` (immutable). To modify, call
/// [`Certified::into_evidence`]: the downgrade is explicit, the
/// certification mark is LOST, and re-certifying re-validates. The
/// escape hatches are plain `Evidence<T>` or
/// [`NumericalCertificate::no_claim`] — never a `Certified<T>`.
#[derive(Debug, Clone, PartialEq)]
pub struct Certified<T>(Evidence<T>);

impl<T> Certified<T> {
    /// Read-only view of the certified evidence.
    #[must_use]
    pub fn evidence(&self) -> &Evidence<T> {
        &self.0
    }

    /// EXPLICIT downgrade: surrender the certification mark to get a
    /// mutable `Evidence<T>` back. Any change must pass
    /// [`Evidence::certified`] again to regain `Certified<T>` — this is
    /// the reconstruction/round-trip path, and it re-validates.
    #[must_use]
    pub fn into_evidence(self) -> Evidence<T> {
        self.0
    }
}

impl<T> core::ops::Deref for Certified<T> {
    type Target = Evidence<T>;
    fn deref(&self) -> &Evidence<T> {
        &self.0
    }
}

/// Why a value could not be certified (Decalogue P10: teaching refusals).
/// Carries the offending numbers; `Eq` is deliberately absent (f64
/// payloads), `Copy + PartialEq` retained.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CertifyError {
    /// The numerical certificate is not rigorous.
    NotRigorous {
        /// The offending kind.
        kind: NumericalKind,
    },
    /// A constituent model was used outside its validity domain.
    OutOfDomain,
    /// An Exact certificate over a non-finite value.
    ExactNotFinite {
        /// The rejected QoI.
        qoi: f64,
    },
    /// An Exact certificate whose bounds do not both equal the QoI.
    ExactInconsistent {
        /// The QoI the certificates describe.
        qoi: f64,
        /// The certificate's lower bound.
        lo: f64,
        /// The certificate's upper bound.
        hi: f64,
    },
    /// An enclosure with a NaN or infinite bound.
    NonFiniteBounds {
        /// The certificate's lower bound.
        lo: f64,
        /// The certificate's upper bound.
        hi: f64,
    },
    /// An enclosure with `lo > hi` (encloses nothing).
    InvertedBounds {
        /// The certificate's lower bound.
        lo: f64,
        /// The certificate's upper bound.
        hi: f64,
    },
    /// The QoI lies outside its own claimed enclosure.
    QoiOutsideEnclosure {
        /// The QoI the certificates describe.
        qoi: f64,
        /// The certificate's lower bound.
        lo: f64,
        /// The certificate's upper bound.
        hi: f64,
    },
    /// Scalar evidence whose carried value differs from the certified QoI.
    ScalarValueMismatch {
        /// The carried scalar value.
        value: f64,
        /// The scalar QoI described by the certificate.
        qoi: f64,
    },
    /// A statistical certificate contains an invalid uncertainty value.
    InvalidStatistical {
        /// The invalid field.
        field: &'static str,
        /// The rejected value.
        value: f64,
        /// The field's validity requirement.
        requirement: &'static str,
    },
    /// A model-form discrepancy is NaN or negative.
    InvalidModelDiscrepancy {
        /// The rejected relative discrepancy.
        discrepancy_rel: f64,
    },
    /// The declared model-validity box is empty or has a non-finite bound.
    InvalidModelValidity,
}

impl fmt::Display for CertifyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CertifyError::NotRigorous { kind } => write!(
                f,
                "Certified<T> requires a rigorous numerical certificate; got {} — enclose the \
                 QoI (fs-ivl) or keep it as plain Evidence",
                kind.name()
            ),
            CertifyError::OutOfDomain => write!(
                f,
                "Certified<T> refused: a constituent model was used OUTSIDE its validity \
                 domain; bracket the model or escalate fidelity"
            ),
            CertifyError::ExactNotFinite { qoi } => write!(
                f,
                "Certified<T> refused: an Exact certificate over the non-finite value {qoi} \
                 certifies nothing — use NumericalCertificate::no_claim for missing values"
            ),
            CertifyError::ExactInconsistent { qoi, lo, hi } => write!(
                f,
                "Certified<T> refused: Exact demands lo == qoi == hi, got qoi {qoi} with \
                 bounds [{lo}, {hi}] — the certificate does not describe its own QoI"
            ),
            CertifyError::NonFiniteBounds { lo, hi } => write!(
                f,
                "Certified<T> refused: enclosure [{lo}, {hi}] has a NaN or infinite bound; \
                 a certified enclosure must be finite — keep it as Evidence or no-claim"
            ),
            CertifyError::InvertedBounds { lo, hi } => write!(
                f,
                "Certified<T> refused: enclosure [{lo}, {hi}] is inverted (lo > hi) and \
                 encloses NOTHING — a forged or corrupted certificate fails closed"
            ),
            CertifyError::QoiOutsideEnclosure { qoi, lo, hi } => write!(
                f,
                "Certified<T> refused: qoi {qoi} lies outside its own claimed enclosure \
                 [{lo}, {hi}] — the certificate contradicts the value it travels with"
            ),
            CertifyError::ScalarValueMismatch { value, qoi } => write!(
                f,
                "Certified<f64> refused: carried value {value} differs from certificate QoI \
                 {qoi} — scalar evidence must describe the number it carries"
            ),
            CertifyError::InvalidStatistical {
                field,
                value,
                requirement,
            } => write!(
                f,
                "Certified<T> refused: statistical field `{field}` = {value} is invalid; it \
                 {requirement} — keep the value as plain Evidence until uncertainty is valid"
            ),
            CertifyError::InvalidModelDiscrepancy { discrepancy_rel } => write!(
                f,
                "Certified<T> refused: model discrepancy {discrepancy_rel} is NaN or negative; \
                 relative discrepancy must be non-negative (infinity is an honest unbounded \
                 claim)"
            ),
            CertifyError::InvalidModelValidity => write!(
                f,
                "Certified<T> refused: the model validity domain is empty or contains a \
                 non-finite bound; an impossible regime cannot be asserted in-domain"
            ),
        }
    }
}

impl core::error::Error for CertifyError {}

impl Evidence<f64> {
    /// Scalar evidence with an exact numerical certificate.
    #[must_use]
    pub fn exact(value: f64, provenance: ProvenanceHash) -> Self {
        Evidence {
            value,
            qoi: value,
            numerical: NumericalCertificate::exact(value),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        }
    }

    /// Scalar evidence with a rigorous enclosure.
    #[must_use]
    pub fn enclosed(value: f64, lo: f64, hi: f64, provenance: ProvenanceHash) -> Self {
        Evidence {
            numerical: NumericalCertificate::enclosure(lo, hi),
            ..Evidence::exact(value, provenance)
        }
    }
}

impl<T> Evidence<T> {
    /// The `Certified<T>` constructor discipline — the ONLY door into
    /// [`Certified<T>`], and the FULL validation boundary (bead
    /// gp3.2.1): rigorous numerics (Exact/Enclosure), finite ordered
    /// bounds, the QoI contained in its own enclosure, valid statistical
    /// parameters, and EXPLICIT in-domain non-negative model evidence.
    /// Pure-math values certify with
    /// `ModelEvidence::none()` — that IS the explicit "no model
    /// involved" statement.
    ///
    /// Because `NumericalCertificate` fields are public (composition
    /// ergonomics), a hand-built certificate can hold NaN, infinite,
    /// or inverted bounds — everything here validates the ACTUAL
    /// numbers, not the constructor that claimed them.
    ///
    /// Certified payloads must be `'static`: scalar evidence is identified
    /// at this boundary so its carried `f64` can be bound bit-exactly to the
    /// certified QoI. Borrow-scoped payloads remain plain `Evidence<T>` until
    /// the caller promotes an owned value.
    ///
    /// # Errors
    /// [`CertifyError`] naming what is missing and how to fix it.
    pub fn certified(self) -> Result<Certified<T>, CertifyError>
    where
        T: 'static,
    {
        let (qoi, lo, hi) = (self.qoi, self.numerical.lo, self.numerical.hi);
        if let Some(value) = (&self.value as &dyn core::any::Any).downcast_ref::<f64>()
            && value.to_bits() != qoi.to_bits()
        {
            return Err(CertifyError::ScalarValueMismatch { value: *value, qoi });
        }
        match self.numerical.kind {
            NumericalKind::Estimate | NumericalKind::NoClaim => {
                return Err(CertifyError::NotRigorous {
                    kind: self.numerical.kind,
                });
            }
            NumericalKind::Exact => {
                if !qoi.is_finite() {
                    return Err(CertifyError::ExactNotFinite { qoi });
                }
                // Bit identity, not approximate equality: Exact CLAIMS
                // the bounds are the value (and fails closed on the
                // -0.0 vs 0.0 edge rather than minting a certificate).
                if lo.to_bits() != qoi.to_bits() || hi.to_bits() != qoi.to_bits() {
                    return Err(CertifyError::ExactInconsistent { qoi, lo, hi });
                }
            }
            NumericalKind::Enclosure => {
                if !(lo.is_finite() && hi.is_finite()) {
                    return Err(CertifyError::NonFiniteBounds { lo, hi });
                }
                if lo > hi {
                    return Err(CertifyError::InvertedBounds { lo, hi });
                }
                if !(qoi >= lo && qoi <= hi) {
                    return Err(CertifyError::QoiOutsideEnclosure { qoi, lo, hi });
                }
            }
        }
        if self.model.validity.is_empty() {
            return Err(CertifyError::InvalidModelValidity);
        }
        if !self.model.in_domain {
            return Err(CertifyError::OutOfDomain);
        }
        if let Some((field, value, requirement)) = self.statistical.validation_issue() {
            return Err(CertifyError::InvalidStatistical {
                field,
                value,
                requirement,
            });
        }
        if !ModelEvidence::valid_discrepancy(self.model.discrepancy_rel) {
            return Err(CertifyError::InvalidModelDiscrepancy {
                discrepancy_rel: self.model.discrepancy_rel,
            });
        }
        Ok(Certified(self))
    }

    /// Attach a statistical certificate.
    #[must_use]
    pub fn with_statistical(mut self, s: StatisticalCertificate) -> Self {
        self.statistical = s;
        self
    }

    /// Attach model evidence.
    #[must_use]
    pub fn with_model(mut self, m: ModelEvidence) -> Self {
        self.model = m;
        self
    }

    /// Attach a sensitivity headline.
    #[must_use]
    pub fn with_sensitivity(mut self, s: SensitivitySummary) -> Self {
        self.sensitivity = s;
        self
    }

    /// Attach an adjoint hook.
    #[must_use]
    pub fn with_adjoint(mut self, adjoint: ProvenanceHash) -> Self {
        self.adjoint_ref = Some(adjoint);
        self
    }

    /// Per-source relative bands (out-of-domain model use is an infinite
    /// model band — un-decidable by construction).
    #[must_use]
    pub fn breakdown(&self) -> UncertaintyBreakdown {
        UncertaintyBreakdown {
            numerical_rel: self.numerical.rel_half_width(self.qoi),
            statistical_rel: self.statistical.rel_width(self.qoi),
            model_rel: if self.model.in_domain
                && !self.model.validity.is_empty()
                && ModelEvidence::valid_discrepancy(self.model.discrepancy_rel)
            {
                self.model.discrepancy_rel
            } else {
                f64::INFINITY
            },
        }
    }

    /// Is this decision-grade at relative threshold `threshold_rel`?
    /// The teaching detail names the dominant source and the numbers.
    #[must_use]
    pub fn assess(&self, threshold_rel: f64) -> DecisionStatus {
        let b = self.breakdown();
        let total_rel = b.total_rel();
        let valid_threshold = threshold_rel.is_finite() && threshold_rel >= 0.0;
        if valid_threshold && total_rel.is_finite() && total_rel <= threshold_rel {
            return DecisionStatus::DecisionGrade;
        }
        let dominant = b.dominant();
        let suffix = if dominant == UncertaintySource::ModelForm {
            "; cheap refinement cannot fix this, escalate model fidelity or bracket"
        } else {
            ""
        };
        let detail = if valid_threshold {
            format!(
                "total relative band {:.3} exceeds the decision threshold {:.3}: numerical \
                 {:.3}, statistical {:.3}, model-form {:.3} — {} uncertainty dominates{}",
                total_rel,
                threshold_rel,
                b.numerical_rel,
                b.statistical_rel,
                b.model_rel,
                dominant.name(),
                suffix,
            )
        } else {
            format!(
                "decision threshold {threshold_rel} is invalid (must be finite and \
                 non-negative): numerical {:.3}, statistical {:.3}, model-form {:.3} — {} \
                 uncertainty dominates{}",
                b.numerical_rel,
                b.statistical_rel,
                b.model_rel,
                dominant.name(),
                suffix,
            )
        };
        DecisionStatus::NotDecisionGrade { dominant, detail }
    }

    /// Decision-aware escalation advice at `threshold_rel` (the HELM
    /// governor hook).
    #[must_use]
    pub fn escalation_advice(&self, threshold_rel: f64) -> EscalationAdvice {
        match self.assess(threshold_rel) {
            DecisionStatus::DecisionGrade => EscalationAdvice::NoneNeeded,
            DecisionStatus::NotDecisionGrade { dominant, .. } => match dominant {
                UncertaintySource::ModelForm => EscalationAdvice::EscalateModelFidelity,
                UncertaintySource::Statistical => EscalationAdvice::GatherMoreSamples,
                UncertaintySource::Numerical => EscalationAdvice::RefineNumerics,
            },
        }
    }

    /// Conservative composition: `op` on the QoIs, certificates composed
    /// per the module-level laws, provenance chained. The adjoint hook
    /// does NOT compose (fs-ad owns composed tapes — no-claim).
    #[must_use]
    pub fn combine<U, V>(op: Op, a: &Evidence<T>, b: &Evidence<U>, value: V) -> Evidence<V> {
        Evidence {
            value,
            qoi: op.apply(a.qoi, b.qoi),
            numerical: NumericalCertificate::compose(&a.numerical, &b.numerical, op),
            statistical: StatisticalCertificate::compose(&a.statistical, &b.statistical),
            model: ModelEvidence::compose(&a.model, &b.model),
            sensitivity: SensitivitySummary::merge(&a.sensitivity, &b.sensitivity),
            provenance: ProvenanceHash::chain(op.name(), &[a.provenance, b.provenance]),
            adjoint_ref: None,
        }
    }

    /// The ledger `evidence` table row / `explain()` payload (canonical
    /// field order, deterministic — no clocks, no addresses).
    #[must_use]
    pub fn to_ledger_row_json(&self) -> String {
        let mut s = String::with_capacity(384);
        let _ = write!(
            s,
            "{{\"provenance\":\"{:016x}\",\"qoi\":{},\"numerical\":{{\"kind\":\"{}\",\
             \"lo\":{},\"hi\":{}}},\"statistical\":",
            self.provenance.0,
            fmt_f64(self.qoi),
            self.numerical.kind.name(),
            fmt_f64(self.numerical.lo),
            fmt_f64(self.numerical.hi),
        );
        match self.statistical {
            StatisticalCertificate::None => s.push_str("{\"kind\":\"none\"}"),
            StatisticalCertificate::EValue { e, alpha } => {
                let _ = write!(
                    s,
                    "{{\"kind\":\"e-value\",\"e\":{},\"alpha\":{}}}",
                    fmt_f64(e),
                    fmt_f64(alpha)
                );
            }
            StatisticalCertificate::HalfWidth {
                half_width,
                confidence,
            } => {
                let _ = write!(
                    s,
                    "{{\"kind\":\"half-width\",\"half_width\":{},\"confidence\":{}}}",
                    fmt_f64(half_width),
                    fmt_f64(confidence)
                );
            }
        }
        let _ = write!(
            s,
            ",\"model\":{{\"cards\":[{}],\"assumptions\":[{}],\"discrepancy_rel\":{},\
             \"in_domain\":{},\"validity\":{}}},\"sensitivity\":{},\"adjoint\":{}}}",
            canonical_json_string_list(&self.model.cards),
            canonical_json_string_list(&self.model.assumptions),
            fmt_f64(self.model.discrepancy_rel),
            self.model.in_domain,
            self.model.validity.to_json(),
            sensitivity_json(&self.sensitivity),
            self.adjoint_ref
                .map_or_else(|| "null".to_string(), |a| format!("\"{:016x}\"", a.0)),
        );
        s
    }
}

/// Finite floats print shortest-round-trip; non-finite prints as a tagged
/// string (fs-obs doctrine: JSON has no Inf/NaN).
pub(crate) fn fmt_f64(v: f64) -> String {
    if v.is_finite() {
        format!("{v}")
    } else {
        format!("\"non-finite:{v}\"")
    }
}

pub(crate) fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Render a set-like public string vector in canonical lexical order.
pub(crate) fn canonical_json_string_list(values: &[String]) -> String {
    let mut canonical: Vec<&str> = values.iter().map(String::as_str).collect();
    canonical.sort_unstable();
    canonical.dedup();
    canonical
        .into_iter()
        .map(json_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn sensitivity_json(summary: &SensitivitySummary) -> String {
    let mut out = String::from("{");
    for (index, (parameter, derivative)) in summary.d_qoi.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "{}:{}", json_string(parameter), fmt_f64(*derivative));
    }
    out.push('}');
    out
}

impl fmt::Display for DecisionStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DecisionStatus::DecisionGrade => write!(f, "decision-grade"),
            DecisionStatus::NotDecisionGrade { detail, .. } => {
                write!(f, "NOT decision-grade: {detail}")
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(pairs: &[(&str, f64)]) -> BTreeMap<String, f64> {
        pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn numerical_composition_is_outward_and_never_exact_after_floats() {
        let a = NumericalCertificate::exact(2.0);
        let b = NumericalCertificate::enclosure(3.0, 3.5);
        let c = NumericalCertificate::compose(&a, &b, Op::Mul);
        assert_eq!(c.kind, NumericalKind::Enclosure);
        assert!(c.lo <= 6.0 && c.hi >= 7.0, "{c:?}");
        let nc = NumericalCertificate::compose(&a, &NumericalCertificate::no_claim(), Op::Add);
        assert_eq!(nc.kind, NumericalKind::NoClaim, "NoClaim absorbs");
        assert!(nc.lo.is_infinite() && nc.hi.is_infinite());
    }

    #[test]
    fn validity_intersection_and_containment_laws() {
        let a = ValidityDomain::unconstrained().with("Re", 1e4, 1e5);
        let b = ValidityDomain::unconstrained()
            .with("Re", 5e4, 2e5)
            .with("Ma", 0.0, 0.3);
        let i = a.intersect(&b);
        assert_eq!(i.bound("Re"), Some((5e4, 1e5)));
        assert_eq!(i.bound("Ma"), Some((0.0, 0.3)));
        let p_in = point(&[("Re", 7e4), ("Ma", 0.1)]);
        let p_out = point(&[("Re", 2e4), ("Ma", 0.1)]);
        assert!(i.contains(&p_in));
        assert!(!i.contains(&p_out));
        // Composed ⊆ both constituents (the conservativeness law).
        assert!(a.contains(&p_in) && b.contains(&p_in));
        // Empty intersection is detectable.
        let e = a.intersect(&ValidityDomain::unconstrained().with("Re", 2e5, 3e5));
        assert!(e.is_empty());
        // Missing keys are unconstrained.
        assert!(a.contains(&point(&[("Re", 5e4)])));
        assert!(
            !a.contains(&point(&[("Ma", 0.1)])),
            "missing constrained key fails"
        );
    }

    #[test]
    fn certified_discipline_refuses_estimates_and_out_of_domain_models() {
        let p = ProvenanceHash::of_bytes(b"t");
        let est = Evidence {
            numerical: NumericalCertificate::estimate(1.0, 2.0),
            ..Evidence::exact(1.5, p)
        };
        let err = est.certified().expect_err("estimates are not certified");
        assert!(err.to_string().contains("rigorous"), "{err}");
        let mut out_of_domain = Evidence::exact(1.5, p);
        out_of_domain.model = ModelEvidence {
            cards: vec!["les-smagorinsky".into()],
            assumptions: vec![],
            validity: ValidityDomain::unconstrained().with("Re", 0.0, 1e5),
            discrepancy_rel: 0.1,
            in_domain: false,
        };
        let err = out_of_domain
            .certified()
            .expect_err("out of domain refused");
        assert!(err.to_string().contains("validity"), "{err}");
        assert!(
            Evidence::exact(1.5, p).certified().is_ok(),
            "pure math certifies"
        );
    }

    #[test]
    fn dominant_source_ties_break_to_model_form() {
        let b = UncertaintyBreakdown {
            numerical_rel: 0.1,
            statistical_rel: 0.1,
            model_rel: 0.1,
        };
        assert_eq!(b.dominant(), UncertaintySource::ModelForm);
        let b2 = UncertaintyBreakdown {
            numerical_rel: 0.2,
            statistical_rel: 0.1,
            model_rel: 0.1,
        };
        assert_eq!(b2.dominant(), UncertaintySource::Numerical);
    }

    #[test]
    fn provenance_chains_are_order_sensitive_and_deterministic() {
        let a = ProvenanceHash::of_bytes(b"a");
        let b = ProvenanceHash::of_bytes(b"b");
        assert_eq!(
            ProvenanceHash::chain("add", &[a, b]),
            ProvenanceHash::chain("add", &[a, b])
        );
        assert_ne!(
            ProvenanceHash::chain("sub", &[a, b]),
            ProvenanceHash::chain("sub", &[b, a])
        );
    }

    #[test]
    fn ledger_row_is_canonical_and_single_line() {
        let p = ProvenanceHash::of_bytes(b"row");
        let ev = Evidence::enclosed(12.4, 12.3, 12.5, p)
            .with_statistical(StatisticalCertificate::HalfWidth {
                half_width: 0.05,
                confidence: 0.95,
            })
            .with_adjoint(ProvenanceHash::of_bytes(b"tape"));
        let row = ev.to_ledger_row_json();
        assert_eq!(row, ev.to_ledger_row_json(), "deterministic");
        assert!(!row.contains('\n'));
        for key in [
            "\"provenance\":",
            "\"numerical\":",
            "\"statistical\":",
            "\"model\":",
            "\"sensitivity\":",
            "\"adjoint\":",
        ] {
            assert!(row.contains(key), "missing {key} in {row}");
        }
    }

    #[test]
    fn ledger_row_escapes_all_evidence_metadata_and_keeps_missing_slices() {
        let mut sensitivity = SensitivitySummary::default();
        sensitivity
            .d_qoi
            .insert("gain\"\n\u{0001}".to_string(), f64::NAN);
        let evidence = Evidence::exact(1.0, ProvenanceHash::of_bytes(b"hostile-row"))
            .with_model(ModelEvidence {
                cards: vec!["card\"\n\u{0002}".to_string()],
                assumptions: vec!["assume\\\r\u{0003}".to_string()],
                validity: ValidityDomain::unconstrained().with(
                    "axis\"\n\u{0004}",
                    f64::NEG_INFINITY,
                    f64::INFINITY,
                ),
                discrepancy_rel: f64::INFINITY,
                in_domain: false,
            })
            .with_sensitivity(sensitivity);

        let row = evidence.to_ledger_row_json();
        assert!(row.contains("\"assumptions\":["), "{row}");
        assert!(row.contains("\"sensitivity\":{"), "{row}");
        assert!(row.contains("card\\\"\\n\\u0002"), "{row}");
        assert!(row.contains("axis\\\"\\n\\u0004"), "{row}");
        assert!(row.contains("gain\\\"\\n\\u0001"), "{row}");
        assert!(row.contains("\"non-finite:inf\""), "{row}");
        assert!(row.contains("\"non-finite:NaN\""), "{row}");
        assert!(!row.chars().any(|ch| u32::from(ch) < 0x20), "{row:?}");
    }

    #[test]
    fn ledger_row_canonicalizes_public_model_sets() {
        let provenance = ProvenanceHash::of_bytes(b"canonical-public-model-sets");
        let base = Evidence::exact(1.0, provenance);
        let first = base.clone().with_model(ModelEvidence {
            cards: vec!["zeta".to_string(), "alpha".to_string(), "zeta".to_string()],
            assumptions: vec![
                "steady".to_string(),
                "isothermal".to_string(),
                "steady".to_string(),
            ],
            validity: ValidityDomain::unconstrained(),
            discrepancy_rel: 0.1,
            in_domain: true,
        });
        let second = base.with_model(ModelEvidence {
            cards: vec!["alpha".to_string(), "zeta".to_string()],
            assumptions: vec!["isothermal".to_string(), "steady".to_string()],
            validity: ValidityDomain::unconstrained(),
            discrepancy_rel: 0.1,
            in_domain: true,
        });

        let first_row = first.to_ledger_row_json();
        assert_eq!(
            first_row,
            second.to_ledger_row_json(),
            "caller ordering and duplicates cannot change a set-like durable row"
        );
        assert_eq!(first_row.matches("zeta").count(), 1, "duplicates survive");
        let alpha = first_row.find("alpha").expect("alpha card retained");
        let zeta = first_row.find("zeta").expect("zeta card retained");
        assert!(
            alpha < zeta,
            "sets are not lexically canonical: {first_row}"
        );
    }
}
