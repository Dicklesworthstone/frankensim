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

mod cards;
mod discrepancy;

pub use cards::{Ambition, ModelCard, ModelRegistry, RegistryError};
pub use discrepancy::{
    DiscrepancyBand, DiscrepancyModel, FidelityPair, FitError, ModelBracket, OutOfDomain,
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
    #[must_use]
    pub fn chain(op: &str, operands: &[ProvenanceHash]) -> Self {
        let mut s = String::with_capacity(16 + operands.len() * 17);
        s.push_str(op);
        for p in operands {
            let _ = write!(s, "|{:016x}", p.0);
        }
        Self::of_bytes(s.as_bytes())
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

    /// A rigorous enclosure (callers guarantee `lo <= hi`; violations are
    /// normalized by swapping — a teaching-free total function).
    #[must_use]
    pub fn enclosure(lo: f64, hi: f64) -> Self {
        NumericalCertificate {
            kind: NumericalKind::Enclosure,
            lo: lo.min(hi),
            hi: lo.max(hi),
        }
    }

    /// A non-rigorous band.
    #[must_use]
    pub fn estimate(lo: f64, hi: f64) -> Self {
        NumericalCertificate {
            kind: NumericalKind::Estimate,
            lo: lo.min(hi),
            hi: lo.max(hi),
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
        if self.kind == NumericalKind::NoClaim {
            return f64::INFINITY;
        }
        let half = 0.5 * (self.hi - self.lo);
        half / reference.abs().max(f64::MIN_POSITIVE)
    }

    fn compose(a: &Self, b: &Self, op: Op) -> Self {
        if a.kind == NumericalKind::NoClaim || b.kind == NumericalKind::NoClaim {
            return Self::no_claim();
        }
        let (lo, hi) = match op {
            Op::Add => (a.lo + b.lo, a.hi + b.hi),
            Op::Sub => (a.lo - b.hi, a.hi - b.lo),
            Op::Mul => {
                let c = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi];
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
        let kind = a.kind.max(b.kind).max(NumericalKind::Enclosure);
        NumericalCertificate {
            kind,
            lo: lo.next_down(),
            hi: hi.next_up(),
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
        /// The e-value.
        e: f64,
        /// The design level the e-threshold was set for.
        alpha: f64,
    },
    /// A confidence half-width around the QoI.
    HalfWidth {
        /// Absolute half-width.
        half_width: f64,
        /// Confidence level in (0, 1).
        confidence: f64,
    },
}

impl StatisticalCertificate {
    /// Relative statistical band against a reference magnitude. E-values
    /// contribute zero width (the decision they certify is already
    /// anytime-valid); half-widths scale by the reference.
    #[must_use]
    pub fn rel_width(&self, reference: f64) -> f64 {
        match self {
            StatisticalCertificate::None | StatisticalCertificate::EValue { .. } => 0.0,
            StatisticalCertificate::HalfWidth { half_width, .. } => {
                half_width / reference.abs().max(f64::MIN_POSITIVE)
            }
        }
    }

    fn compose(a: &Self, b: &Self) -> Self {
        use StatisticalCertificate as S;
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

    /// Constrain one parameter to `[lo, hi]`.
    #[must_use]
    pub fn with(mut self, param: impl Into<String>, lo: f64, hi: f64) -> Self {
        self.bounds.insert(param.into(), (lo.min(hi), lo.max(hi)));
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
        self.bounds
            .iter()
            .all(|(k, &(lo, hi))| point.get(k).is_some_and(|&v| v >= lo && v <= hi))
    }

    /// Per-parameter intersection (the composition law: composed validity
    /// is EXACTLY the intersection — never wider).
    #[must_use]
    pub fn intersect(&self, other: &Self) -> Self {
        let mut out = self.bounds.clone();
        for (k, &(lo2, hi2)) in &other.bounds {
            out.entry(k.clone())
                .and_modify(|(lo, hi)| {
                    *lo = lo.max(lo2);
                    *hi = hi.min(hi2);
                })
                .or_insert((lo2, hi2));
        }
        ValidityDomain { bounds: out }
    }

    /// True when some parameter's interval is empty (nothing satisfies it).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bounds.values().any(|&(lo, hi)| lo > hi)
    }

    /// Constrained parameter names, sorted (BTreeMap order — deterministic
    /// renderings and audits key off this).
    #[must_use]
    pub fn param_names(&self) -> Vec<String> {
        self.bounds.keys().cloned().collect()
    }

    fn to_json(&self) -> String {
        let mut s = String::from("{");
        for (i, (k, (lo, hi))) in self.bounds.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "\"{k}\":[{lo},{hi}]");
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
    /// Relative model-form discrepancy band (composition adds —
    /// first-order conservative).
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
            discrepancy_rel: a.discrepancy_rel + b.discrepancy_rel,
            in_domain: a.in_domain && b.in_domain,
        }
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
        self.numerical_rel + self.statistical_rel + self.model_rel
    }

    /// The dominant source (declaration-order tie-break).
    #[must_use]
    pub fn dominant(&self) -> UncertaintySource {
        let entries = [
            (UncertaintySource::ModelForm, self.model_rel),
            (UncertaintySource::Statistical, self.statistical_rel),
            (UncertaintySource::Numerical, self.numerical_rel),
        ];
        let mut best = entries[0];
        for &(src, band) in &entries[1..] {
            if band > best.1 {
                best = (src, band);
            }
        }
        best.0
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

/// `Certified<T>`: evidence whose numerical slice is rigorous (Exact or
/// Enclosure) AND whose model evidence is explicit. The constructor
/// discipline lives in [`Evidence::certified`].
pub type Certified<T> = Evidence<T>;

/// Why a value could not be certified (Decalogue P10: teaching refusals).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CertifyError {
    /// The numerical certificate is not rigorous.
    NotRigorous {
        /// The offending kind.
        kind: NumericalKind,
    },
    /// A constituent model was used outside its validity domain.
    OutOfDomain,
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
    /// The `Certified<T>` constructor discipline: rigorous numerics
    /// (Exact/Enclosure) plus EXPLICIT, in-domain model evidence.
    /// Pure-math values certify with `ModelEvidence::none()` — that IS the
    /// explicit "no model involved" statement.
    ///
    /// # Errors
    /// [`CertifyError`] naming what is missing and how to fix it.
    pub fn certified(self) -> Result<Certified<T>, CertifyError> {
        if matches!(
            self.numerical.kind,
            NumericalKind::Estimate | NumericalKind::NoClaim
        ) {
            return Err(CertifyError::NotRigorous {
                kind: self.numerical.kind,
            });
        }
        if !self.model.in_domain {
            return Err(CertifyError::OutOfDomain);
        }
        Ok(self)
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
            model_rel: if self.model.in_domain {
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
        if b.total_rel() <= threshold_rel {
            return DecisionStatus::DecisionGrade;
        }
        let dominant = b.dominant();
        let detail = format!(
            "total relative band {:.3} exceeds the decision threshold {:.3}: numerical {:.3}, \
             statistical {:.3}, model-form {:.3} — {} uncertainty dominates{}",
            b.total_rel(),
            threshold_rel,
            b.numerical_rel,
            b.statistical_rel,
            b.model_rel,
            dominant.name(),
            if dominant == UncertaintySource::ModelForm {
                "; cheap refinement cannot fix this, escalate model fidelity or bracket"
            } else {
                ""
            }
        );
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
            ",\"model\":{{\"cards\":[{}],\"discrepancy_rel\":{},\"in_domain\":{},\
             \"validity\":{}}},\"adjoint\":{}}}",
            self.model
                .cards
                .iter()
                .map(|c| format!("\"{c}\""))
                .collect::<Vec<_>>()
                .join(","),
            fmt_f64(self.model.discrepancy_rel),
            self.model.in_domain,
            self.model.validity.to_json(),
            self.adjoint_ref
                .map_or_else(|| "null".to_string(), |a| format!("\"{:016x}\"", a.0)),
        );
        s
    }
}

/// Finite floats print shortest-round-trip; non-finite prints as a tagged
/// string (fs-obs doctrine: JSON has no Inf/NaN).
fn fmt_f64(v: f64) -> String {
    if v.is_finite() {
        format!("{v}")
    } else {
        format!("\"non-finite:{v}\"")
    }
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
            "\"adjoint\":",
        ] {
            assert!(row.contains(key), "missing {key} in {row}");
        }
    }
}
