//! I08.1 evidence-action vocabulary: proposals, responses, dependencies,
//! and budgets for buying uncertainty reduction.
//!
//! Compute refinements and physical experiments share one decision-facing
//! vocabulary ([`ActionProposal`]) without pretending their costs, lead
//! times, or evidence colors are interchangeable. Every proposal targets one
//! claim slice, carries explicit cost states per axis (an unknown cost stays
//! [`CostState::Unknown`] with a named authority gap — never a silent zero),
//! declares dependencies, exclusivity, correlation, expiry, and an evidence
//! color CEILING, and has a stable content-addressed identity so identical
//! proposals are idempotent.
//!
//! The planning/execution split is structural: an [`ActionProposal`] can
//! never carry an outcome color; only an [`ExecutionReceipt`] admitted
//! against the proposal can, and its outcome can never exceed the ceiling.
//! A planned physical test therefore cannot raise evidence until executed.
//!
//! No-claims: cost and response "distributions" are interval envelopes, not
//! parametric distributions; price normalization across currencies is a
//! policy input this module refuses to guess; correlation is carried as
//! group identity only, with no joint-distribution claim; portfolio
//! admission checks feasibility structure (cycles, exclusivity, expiry,
//! duplicate identity), not decision optimality.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};

use crate::color::ColorRank;

pub use crate::balance::BoundedId;

/// Current canonical schema version for evidence-action proposals.
pub const ACTION_SCHEMA_VERSION: u32 = 1;
/// Maximum dependencies one proposal may declare.
pub const MAX_ACTION_DEPENDENCIES: usize = 256;
/// Maximum proposals one portfolio may admit.
pub const MAX_PORTFOLIO_ACTIONS: usize = 4096;
/// Maximum capability identities one proposal may demand.
pub const MAX_ACTION_CAPABILITIES: usize = 64;
/// Maximum accepted canonical transport size in bytes.
pub const MAX_ACTION_CANONICAL_BYTES: usize = 256 * 1024;

const IDENTITY_DOMAIN: &str = "org.frankensim.fs-evidence.evidence-action.v1";
const MAGIC: &[u8; 4] = b"FSEA";

/// Stable named rule violated by a refusal. Every refusal names exactly one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ActionRule {
    /// Identifier empty, over budget, or outside the bounded grammar.
    IdBounds,
    /// A numeric field was NaN or infinite.
    NonFiniteValue,
    /// A range lower endpoint exceeds its upper endpoint.
    InvertedRange,
    /// A cost or duration was negative (or negative zero).
    NegativeCost,
    /// A currency code was not exactly three ASCII uppercase letters.
    CurrencyCode,
    /// Cost axes with different units cannot add; normalization is policy.
    UnitMismatch,
    /// An expected response factor left (0, 1].
    ResponseRange,
    /// Two proposals in one portfolio share a content identity.
    DuplicateProposal,
    /// A declared dependency is absent from the portfolio.
    MissingDependency,
    /// The dependency graph contains a cycle.
    DependencyCycle,
    /// Two admitted proposals share an exclusivity group.
    ExclusivityViolation,
    /// The proposal expired at or before the admission instant.
    Expired,
    /// Expiry and admission instant reference different clocks.
    ClockMismatch,
    /// An execution outcome exceeded the proposal's evidence ceiling.
    CeilingExceeded,
    /// A receipt names a proposal identity it was not admitted against.
    ReceiptProposalMismatch,
    /// A collection exceeded its declared budget.
    CollectionBudget,
    /// Schema version unknown to this decoder.
    SchemaVersion,
}

impl ActionRule {
    /// Stable machine-readable slug logged with every refusal.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::IdBounds => "action-id-bounds",
            Self::NonFiniteValue => "action-non-finite-value",
            Self::InvertedRange => "action-inverted-range",
            Self::NegativeCost => "action-negative-cost",
            Self::CurrencyCode => "action-currency-code",
            Self::UnitMismatch => "action-unit-mismatch",
            Self::ResponseRange => "action-response-range",
            Self::DuplicateProposal => "action-duplicate-proposal",
            Self::MissingDependency => "action-missing-dependency",
            Self::DependencyCycle => "action-dependency-cycle",
            Self::ExclusivityViolation => "action-exclusivity-violation",
            Self::Expired => "action-expired",
            Self::ClockMismatch => "action-clock-mismatch",
            Self::CeilingExceeded => "action-ceiling-exceeded",
            Self::ReceiptProposalMismatch => "action-receipt-proposal-mismatch",
            Self::CollectionBudget => "action-collection-budget",
            Self::SchemaVersion => "action-schema-version",
        }
    }
}

/// A typed refusal naming the violated rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ActionError {
    rule: ActionRule,
    detail: String,
}

impl ActionError {
    fn new(rule: ActionRule, detail: impl Into<String>) -> Self {
        Self {
            rule,
            detail: detail.into(),
        }
    }

    /// The violated rule.
    #[must_use]
    pub const fn rule(&self) -> ActionRule {
        self.rule
    }

    /// Human-readable detail; never participates in rule matching.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for ActionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule.slug(), self.detail)
    }
}

impl std::error::Error for ActionError {}

fn refuse<T>(rule: ActionRule, detail: impl Into<String>) -> Result<T, ActionError> {
    Err(ActionError::new(rule, detail))
}

/// What an action buys: the taxonomy is closed and decision-facing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ActionKind {
    /// Tighten a solver tolerance.
    SolverTolerance,
    /// Refine the spatial mesh.
    MeshRefinement,
    /// Refine the time discretization.
    TimeRefinement,
    /// Escalate representation or model fidelity.
    RepresentationEscalation,
    /// Draw additional UQ samples.
    UqSamples,
    /// Run a physical material or coupon test.
    MaterialCouponTest,
    /// Deploy sensors or run a physical measurement campaign.
    SensorCampaign,
    /// Attempt a falsification probe.
    Falsification,
    /// Discharge a standards or process obligation.
    StandardsObligation,
    /// Explicitly refuse: record that no admissible action exists.
    Refusal,
}

impl ActionKind {
    /// Whether executing this action requires the physical world (and so
    /// carries lead time and can never be pre-counted as evidence).
    #[must_use]
    pub const fn is_physical(self) -> bool {
        matches!(self, Self::MaterialCouponTest | Self::SensorCampaign)
    }

    const fn code(self) -> u8 {
        match self {
            Self::SolverTolerance => 1,
            Self::MeshRefinement => 2,
            Self::TimeRefinement => 3,
            Self::RepresentationEscalation => 4,
            Self::UqSamples => 5,
            Self::MaterialCouponTest => 6,
            Self::SensorCampaign => 7,
            Self::Falsification => 8,
            Self::StandardsObligation => 9,
            Self::Refusal => 10,
        }
    }

    fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            1 => Self::SolverTolerance,
            2 => Self::MeshRefinement,
            3 => Self::TimeRefinement,
            4 => Self::RepresentationEscalation,
            5 => Self::UqSamples,
            6 => Self::MaterialCouponTest,
            7 => Self::SensorCampaign,
            8 => Self::Falsification,
            9 => Self::StandardsObligation,
            10 => Self::Refusal,
            _ => return None,
        })
    }
}

/// The uncertainty component one action targets, aligned with the V&V
/// prediction-assessment decomposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum UncertaintyComponent {
    /// Numerical (discretization/solver) uncertainty.
    Numerical,
    /// Material/property data uncertainty.
    Material,
    /// Model-form uncertainty.
    ModelForm,
    /// Measured-data uncertainty.
    Data,
    /// Aleatory variability.
    Aleatory,
    /// Epistemic ignorance.
    Epistemic,
}

const fn component_code(c: UncertaintyComponent) -> u8 {
    match c {
        UncertaintyComponent::Numerical => 1,
        UncertaintyComponent::Material => 2,
        UncertaintyComponent::ModelForm => 3,
        UncertaintyComponent::Data => 4,
        UncertaintyComponent::Aleatory => 5,
        UncertaintyComponent::Epistemic => 6,
    }
}

fn component_from_code(code: u8) -> Option<UncertaintyComponent> {
    Some(match code {
        1 => UncertaintyComponent::Numerical,
        2 => UncertaintyComponent::Material,
        3 => UncertaintyComponent::ModelForm,
        4 => UncertaintyComponent::Data,
        5 => UncertaintyComponent::Aleatory,
        6 => UncertaintyComponent::Epistemic,
        _ => return None,
    })
}

/// The claim slice one action targets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetSlice {
    /// The claim identity whose uncertainty is targeted.
    pub claim: BoundedId,
    /// The targeted uncertainty component.
    pub component: UncertaintyComponent,
}

/// A non-negative bounded range with exact `f64` endpoints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonNegRange {
    lo: f64,
    hi: f64,
}

impl NonNegRange {
    /// Admit a non-negative range or refuse with a named rule.
    pub fn new(lo: f64, hi: f64) -> Result<Self, ActionError> {
        if !lo.is_finite() || !hi.is_finite() {
            return refuse(ActionRule::NonFiniteValue, "range endpoint not finite");
        }
        if lo > hi {
            return refuse(
                ActionRule::InvertedRange,
                format!("range [{lo:e}, {hi:e}] inverted"),
            );
        }
        if lo < 0.0 || lo.is_sign_negative() || hi.is_sign_negative() {
            return refuse(ActionRule::NegativeCost, "range endpoint negative");
        }
        Ok(Self { lo, hi })
    }

    /// Lower endpoint.
    #[must_use]
    pub const fn lo(&self) -> f64 {
        self.lo
    }

    /// Upper endpoint.
    #[must_use]
    pub const fn hi(&self) -> f64 {
        self.hi
    }

    /// Conservative sum: endpoints add, upper widens one ulp outward.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            lo: self.lo + other.lo,
            hi: (self.hi + other.hi).next_up(),
        }
    }
}

/// The unit a cost range is stated in. Addition requires identical units;
/// converting between them (including currencies) is a policy decision this
/// module refuses to make.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum CostUnit {
    /// Physical money in one ISO-4217-style three-letter currency.
    Currency {
        /// Exactly three ASCII uppercase letters.
        code: [u8; 3],
    },
    /// Compute, in CPU core seconds.
    CpuCoreSeconds,
    /// Wall time, in seconds.
    WallSeconds,
    /// Calendar lead time, in days.
    CalendarDays,
    /// Memory, in bytes.
    MemoryBytes,
}

impl CostUnit {
    /// Admit a currency unit from a three-letter code.
    pub fn currency(code: &str) -> Result<Self, ActionError> {
        let bytes = code.as_bytes();
        if bytes.len() != 3 || !bytes.iter().all(u8::is_ascii_uppercase) {
            return refuse(
                ActionRule::CurrencyCode,
                format!("currency code {code:?} is not three ASCII uppercase letters"),
            );
        }
        Ok(Self::Currency {
            code: [bytes[0], bytes[1], bytes[2]],
        })
    }
}

/// One cost axis: a known interval envelope in one unit, or an explicit
/// unknown with a named missing authority. Unknown never silently reads as
/// zero.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum CostState {
    /// Known interval envelope.
    Known {
        /// The envelope.
        range: NonNegRange,
        /// The unit the envelope is stated in.
        unit: CostUnit,
    },
    /// No cost authority exists; the gap is named.
    Unknown {
        /// The named missing authority.
        authority_gap: BoundedId,
    },
}

impl CostState {
    /// Sum two cost states. Known+Known requires identical units;
    /// any Unknown makes the sum Unknown with the (first) named gap —
    /// missing authorities stay explicit through aggregation.
    pub fn add(&self, other: &Self) -> Result<Self, ActionError> {
        match (self, other) {
            (Self::Known { range: a, unit: ua }, Self::Known { range: b, unit: ub }) => {
                if ua != ub {
                    return refuse(
                        ActionRule::UnitMismatch,
                        "cost units differ; normalization is a policy input",
                    );
                }
                Ok(Self::Known {
                    range: a.add(b),
                    unit: ua.clone(),
                })
            }
            (Self::Unknown { authority_gap }, _) | (_, Self::Unknown { authority_gap }) => {
                Ok(Self::Unknown {
                    authority_gap: authority_gap.clone(),
                })
            }
        }
    }
}

/// The demand vector of one proposal across the non-interchangeable axes.
#[derive(Debug, Clone, PartialEq)]
pub struct Budgets {
    /// Physical money.
    pub money: CostState,
    /// Compute time.
    pub compute: CostState,
    /// Peak memory.
    pub memory: CostState,
    /// Calendar lead time.
    pub lead_time: CostState,
    /// Named capabilities (instruments, licenses, facilities) demanded.
    pub capabilities: BTreeSet<BoundedId>,
}

impl Budgets {
    fn validate(&self) -> Result<(), ActionError> {
        if self.capabilities.len() > MAX_ACTION_CAPABILITIES {
            return refuse(
                ActionRule::CollectionBudget,
                format!("{} capabilities above budget", self.capabilities.len()),
            );
        }
        Ok(())
    }
}

/// Expected uncertainty response: the anticipated remaining-uncertainty
/// factor for the targeted slice, or an explicit unknown. This is a
/// PLANNING expectation and never evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum ExpectedResponse {
    /// Anticipated remaining-uncertainty factor interval in (0, 1].
    Factor {
        /// Optimistic (smallest remaining) factor.
        lo: f64,
        /// Pessimistic (largest remaining) factor.
        hi: f64,
    },
    /// No response model exists.
    Unknown,
}

impl ExpectedResponse {
    fn validate(&self) -> Result<(), ActionError> {
        if let Self::Factor { lo, hi } = self {
            if !lo.is_finite() || !hi.is_finite() {
                return refuse(ActionRule::NonFiniteValue, "response factor not finite");
            }
            if !(*lo > 0.0 && lo <= hi && *hi <= 1.0) {
                return refuse(
                    ActionRule::ResponseRange,
                    format!("response factor [{lo}, {hi}] outside 0 < lo <= hi <= 1"),
                );
            }
        }
        Ok(())
    }
}

/// Expiry instant on a named logical clock.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Expiry {
    /// The logical clock identity.
    pub clock: BoundedId,
    /// The tick at which the proposal expires (exclusive validity).
    pub tick: i64,
}

/// Stable content identity of one admitted proposal.
pub type ActionId = ContentHash;

/// Everything needed to admit one proposal.
#[derive(Debug, Clone)]
pub struct ActionDraft {
    /// What the action buys.
    pub kind: ActionKind,
    /// The claim slice it targets.
    pub target: TargetSlice,
    /// The anticipated response (planning-only).
    pub response: ExpectedResponse,
    /// The demand vector.
    pub budgets: Budgets,
    /// Content identities of proposals this one depends on.
    pub dependencies: BTreeSet<ActionId>,
    /// Mutually exclusive alternatives share one group identity.
    pub exclusivity_group: Option<BoundedId>,
    /// Correlated responses share one group identity (identity only; no
    /// joint-distribution claim).
    pub correlation_group: Option<BoundedId>,
    /// Optional expiry.
    pub expiry: Option<Expiry>,
    /// The strongest evidence color executing this action can EVER yield.
    pub ceiling: ColorRank,
    /// The proposing party.
    pub proposer: BoundedId,
}

/// One admitted, content-addressed evidence-action proposal.
#[derive(Debug, Clone, PartialEq)]
pub struct ActionProposal {
    schema_version: u32,
    kind: ActionKind,
    target: TargetSlice,
    response: ExpectedResponse,
    budgets: Budgets,
    dependencies: BTreeSet<ActionId>,
    exclusivity_group: Option<BoundedId>,
    correlation_group: Option<BoundedId>,
    expiry: Option<Expiry>,
    ceiling: ColorRank,
    proposer: BoundedId,
}

impl ActionProposal {
    /// Admit a proposal or refuse with a named rule.
    pub fn admit(draft: ActionDraft) -> Result<Self, ActionError> {
        draft.response.validate()?;
        draft.budgets.validate()?;
        if draft.dependencies.len() > MAX_ACTION_DEPENDENCIES {
            return refuse(
                ActionRule::CollectionBudget,
                format!("{} dependencies above budget", draft.dependencies.len()),
            );
        }
        Ok(Self {
            schema_version: ACTION_SCHEMA_VERSION,
            kind: draft.kind,
            target: draft.target,
            response: draft.response,
            budgets: draft.budgets,
            dependencies: draft.dependencies,
            exclusivity_group: draft.exclusivity_group,
            correlation_group: draft.correlation_group,
            expiry: draft.expiry,
            ceiling: draft.ceiling,
            proposer: draft.proposer,
        })
    }

    /// Schema version the proposal was admitted under.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// What the action buys.
    #[must_use]
    pub const fn kind(&self) -> ActionKind {
        self.kind
    }

    /// The targeted claim slice.
    #[must_use]
    pub fn target(&self) -> &TargetSlice {
        &self.target
    }

    /// The planning-only anticipated response.
    #[must_use]
    pub const fn response(&self) -> ExpectedResponse {
        self.response
    }

    /// The demand vector.
    #[must_use]
    pub fn budgets(&self) -> &Budgets {
        &self.budgets
    }

    /// Declared dependencies.
    #[must_use]
    pub fn dependencies(&self) -> &BTreeSet<ActionId> {
        &self.dependencies
    }

    /// The exclusivity group, if any.
    #[must_use]
    pub fn exclusivity_group(&self) -> Option<&BoundedId> {
        self.exclusivity_group.as_ref()
    }

    /// The correlation group, if any.
    #[must_use]
    pub fn correlation_group(&self) -> Option<&BoundedId> {
        self.correlation_group.as_ref()
    }

    /// The expiry, if any.
    #[must_use]
    pub fn expiry(&self) -> Option<&Expiry> {
        self.expiry.as_ref()
    }

    /// The evidence color ceiling.
    #[must_use]
    pub const fn ceiling(&self) -> ColorRank {
        self.ceiling
    }

    /// The proposing party.
    #[must_use]
    pub fn proposer(&self) -> &BoundedId {
        &self.proposer
    }

    /// Canonical transport bytes (versioned, bounded, deterministic).
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&self.schema_version.to_be_bytes());
        out.push(self.kind.code());
        push_str(&mut out, self.target.claim.as_str());
        out.push(component_code(self.target.component));
        match self.response {
            ExpectedResponse::Factor { lo, hi } => {
                out.push(1);
                out.extend_from_slice(&lo.to_bits().to_be_bytes());
                out.extend_from_slice(&hi.to_bits().to_be_bytes());
            }
            ExpectedResponse::Unknown => out.push(2),
        }
        for axis in [
            &self.budgets.money,
            &self.budgets.compute,
            &self.budgets.memory,
            &self.budgets.lead_time,
        ] {
            encode_cost(&mut out, axis);
        }
        out.extend_from_slice(&u32_len(self.budgets.capabilities.len()).to_be_bytes());
        for capability in &self.budgets.capabilities {
            push_str(&mut out, capability.as_str());
        }
        out.extend_from_slice(&u32_len(self.dependencies.len()).to_be_bytes());
        for dependency in &self.dependencies {
            out.extend_from_slice(dependency.as_bytes());
        }
        encode_opt_id(&mut out, self.exclusivity_group.as_ref());
        encode_opt_id(&mut out, self.correlation_group.as_ref());
        match &self.expiry {
            None => out.push(0),
            Some(expiry) => {
                out.push(1);
                push_str(&mut out, expiry.clock.as_str());
                out.extend_from_slice(&expiry.tick.to_be_bytes());
            }
        }
        out.push(match self.ceiling {
            ColorRank::Estimated => 1,
            ColorRank::Validated => 2,
            ColorRank::Verified => 3,
        });
        push_str(&mut out, self.proposer.as_str());
        out
    }

    /// Stable content identity: identical proposals are idempotent.
    #[must_use]
    pub fn content_id(&self) -> ActionId {
        hash_domain(IDENTITY_DOMAIN, &self.canonical_bytes())
    }

    /// Decode canonical transport bytes. Fail-closed: bounded,
    /// version-gated, validated, and canonical (re-encoding must reproduce
    /// the input bit-for-bit).
    pub fn decode(bytes: &[u8]) -> Result<Self, ActionCodecError> {
        if bytes.len() > MAX_ACTION_CANONICAL_BYTES {
            return Err(ActionCodecError::at(0, "transport above size budget"));
        }
        let mut r = Reader { bytes, pos: 0 };
        if r.take(4)? != MAGIC {
            return Err(ActionCodecError::at(0, "bad magic"));
        }
        let schema_version = r.u32()?;
        if schema_version != ACTION_SCHEMA_VERSION {
            return Err(ActionCodecError::at(
                r.pos,
                format!("unknown schema version {schema_version}"),
            ));
        }
        let kind = ActionKind::from_code(r.u8()?)
            .ok_or_else(|| ActionCodecError::at(r.pos, "bad kind tag"))?;
        let claim = r.bounded_id()?;
        let component = component_from_code(r.u8()?)
            .ok_or_else(|| ActionCodecError::at(r.pos, "bad component tag"))?;
        let response = match r.u8()? {
            1 => ExpectedResponse::Factor {
                lo: f64::from_bits(r.u64()?),
                hi: f64::from_bits(r.u64()?),
            },
            2 => ExpectedResponse::Unknown,
            other => {
                return Err(ActionCodecError::at(
                    r.pos,
                    format!("bad response tag {other}"),
                ));
            }
        };
        let money = decode_cost(&mut r)?;
        let compute = decode_cost(&mut r)?;
        let memory = decode_cost(&mut r)?;
        let lead_time = decode_cost(&mut r)?;
        let capability_count = r.len_u32(MAX_ACTION_CAPABILITIES)?;
        let mut capabilities = BTreeSet::new();
        let mut prev_capability: Option<BoundedId> = None;
        for _ in 0..capability_count {
            let capability = r.bounded_id()?;
            if prev_capability.as_ref().is_some_and(|p| p >= &capability) {
                return Err(ActionCodecError::at(r.pos, "capabilities not sorted"));
            }
            prev_capability = Some(capability.clone());
            capabilities.insert(capability);
        }
        let dependency_count = r.len_u32(MAX_ACTION_DEPENDENCIES)?;
        let mut dependencies = BTreeSet::new();
        let mut prev_dep: Option<ActionId> = None;
        for _ in 0..dependency_count {
            let dependency = ContentHash::from_slice(r.take(32)?)
                .ok_or_else(|| ActionCodecError::at(r.pos, "bad dependency hash"))?;
            if prev_dep.is_some_and(|p| p.as_bytes() >= dependency.as_bytes()) {
                return Err(ActionCodecError::at(r.pos, "dependencies not sorted"));
            }
            prev_dep = Some(dependency);
            dependencies.insert(dependency);
        }
        let exclusivity_group = decode_opt_id(&mut r)?;
        let correlation_group = decode_opt_id(&mut r)?;
        let expiry = match r.u8()? {
            0 => None,
            1 => Some(Expiry {
                clock: r.bounded_id()?,
                tick: r.i64()?,
            }),
            other => {
                return Err(ActionCodecError::at(
                    r.pos,
                    format!("bad expiry tag {other}"),
                ));
            }
        };
        let ceiling = match r.u8()? {
            1 => ColorRank::Estimated,
            2 => ColorRank::Validated,
            3 => ColorRank::Verified,
            other => {
                return Err(ActionCodecError::at(
                    r.pos,
                    format!("bad ceiling tag {other}"),
                ));
            }
        };
        let proposer = r.bounded_id()?;
        if r.pos != bytes.len() {
            return Err(ActionCodecError::at(r.pos, "trailing bytes"));
        }
        let proposal = Self::admit(ActionDraft {
            kind,
            target: TargetSlice { claim, component },
            response,
            budgets: Budgets {
                money,
                compute,
                memory,
                lead_time,
                capabilities,
            },
            dependencies,
            exclusivity_group,
            correlation_group,
            expiry,
            ceiling,
            proposer,
        })
        .map_err(|e| ActionCodecError::at(0, e.to_string()))?;
        if proposal.canonical_bytes() != bytes {
            return Err(ActionCodecError::at(0, "non-canonical encoding"));
        }
        Ok(proposal)
    }
}

/// Portfolio-level totals for one cost axis: a known same-unit envelope or
/// an explicit incomparable state naming every gap.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PortfolioCost {
    /// All contributions were known and in one unit.
    Known {
        /// Summed envelope.
        range: NonNegRange,
        /// Shared unit.
        unit: CostUnit,
    },
    /// At least one contribution was unknown or unit-incomparable; every
    /// contributing gap or unit is listed so nothing hides.
    Incomparable {
        /// Named authority gaps from unknown contributions.
        gaps: Vec<BoundedId>,
        /// Distinct units seen across known contributions.
        units_seen: usize,
    },
    /// No proposal demanded this axis.
    Zero,
}

/// One admitted feasible portfolio of proposals.
#[derive(Debug, Clone, PartialEq)]
pub struct Portfolio {
    actions: BTreeMap<ActionId, ActionProposal>,
    order: Vec<ActionId>,
    unexecuted_physical: bool,
}

impl Portfolio {
    /// Admit a portfolio at one instant, or refuse with a named rule.
    /// Checks duplicate identity, dependency closure and acyclicity,
    /// exclusivity, expiry against `now`, and collection budgets.
    pub fn admit(proposals: Vec<ActionProposal>, now: &Expiry) -> Result<Self, ActionError> {
        if proposals.len() > MAX_PORTFOLIO_ACTIONS {
            return refuse(
                ActionRule::CollectionBudget,
                format!("{} proposals above budget", proposals.len()),
            );
        }
        let mut actions: BTreeMap<ActionId, ActionProposal> = BTreeMap::new();
        let mut exclusivity: BTreeMap<String, usize> = BTreeMap::new();
        let mut unexecuted_physical = false;
        for proposal in &proposals {
            if let Some(expiry) = proposal.expiry() {
                if expiry.clock != now.clock {
                    return refuse(
                        ActionRule::ClockMismatch,
                        "expiry clock differs from the admission clock",
                    );
                }
                if expiry.tick <= now.tick {
                    return refuse(
                        ActionRule::Expired,
                        format!("proposal expired at tick {}", expiry.tick),
                    );
                }
            }
            if proposal.kind().is_physical() {
                unexecuted_physical = true;
            }
        }
        for proposal in proposals {
            let id = proposal.content_id();
            if let Some(group) = proposal.exclusivity_group() {
                let seen = exclusivity.entry(group.as_str().to_owned()).or_insert(0);
                *seen += 1;
                if *seen > 1 {
                    return refuse(
                        ActionRule::ExclusivityViolation,
                        format!("exclusivity group {:?} admitted twice", group.as_str()),
                    );
                }
            }
            if actions.insert(id, proposal).is_some() {
                return refuse(
                    ActionRule::DuplicateProposal,
                    "identical proposal content admitted twice",
                );
            }
        }
        for (id, proposal) in &actions {
            for dependency in proposal.dependencies() {
                if dependency == id {
                    return refuse(ActionRule::DependencyCycle, "self-dependency");
                }
                if !actions.contains_key(dependency) {
                    return refuse(
                        ActionRule::MissingDependency,
                        format!("dependency {} absent", dependency.to_hex()),
                    );
                }
            }
        }
        // Kahn's algorithm over the dependency DAG; deterministic because
        // both the ready set and adjacency iterate in BTreeMap id order.
        let mut in_degree: BTreeMap<ActionId, usize> = actions
            .iter()
            .map(|(id, p)| (*id, p.dependencies().len()))
            .collect();
        let mut order = Vec::with_capacity(actions.len());
        loop {
            let Some(next) = in_degree.iter().find(|(_, d)| **d == 0).map(|(id, _)| *id) else {
                break;
            };
            in_degree.remove(&next);
            order.push(next);
            for (id, proposal) in &actions {
                if proposal.dependencies().contains(&next)
                    && let Some(d) = in_degree.get_mut(id)
                {
                    *d -= 1;
                }
            }
        }
        if !in_degree.is_empty() {
            return refuse(
                ActionRule::DependencyCycle,
                format!("{} proposals stuck in a dependency cycle", in_degree.len()),
            );
        }
        Ok(Self {
            actions,
            order,
            unexecuted_physical,
        })
    }

    /// Admitted proposals, keyed by content identity.
    #[must_use]
    pub fn actions(&self) -> &BTreeMap<ActionId, ActionProposal> {
        &self.actions
    }

    /// One valid execution order (dependencies first, deterministic).
    #[must_use]
    pub fn order(&self) -> &[ActionId] {
        &self.order
    }

    /// Whether the portfolio contains physical actions that have not been
    /// executed. While true, NO evidence effect exists for those actions —
    /// planned physical tests cannot raise evidence.
    #[must_use]
    pub const fn contains_unexecuted_physical(&self) -> bool {
        self.unexecuted_physical
    }

    /// The evidence effect of the PLAN itself: none, structurally. Only
    /// [`ExecutionReceipt`]s admitted per proposal carry outcome colors.
    #[must_use]
    pub const fn planned_evidence_effect(&self) -> Option<ColorRank> {
        None
    }

    fn total(&self, pick: impl Fn(&Budgets) -> &CostState) -> PortfolioCost {
        let mut known: Option<(NonNegRange, CostUnit)> = None;
        let mut gaps: Vec<BoundedId> = Vec::new();
        let mut units: BTreeSet<CostUnit> = BTreeSet::new();
        for proposal in self.actions.values() {
            match pick(proposal.budgets()) {
                CostState::Known { range, unit } => {
                    units.insert(unit.clone());
                    known = Some(match known {
                        None => (*range, unit.clone()),
                        Some((acc, acc_unit)) => {
                            if &acc_unit != unit {
                                (acc, acc_unit)
                            } else {
                                (acc.add(range), acc_unit)
                            }
                        }
                    });
                }
                CostState::Unknown { authority_gap } => gaps.push(authority_gap.clone()),
            }
        }
        if !gaps.is_empty() || units.len() > 1 {
            return PortfolioCost::Incomparable {
                gaps,
                units_seen: units.len(),
            };
        }
        match known {
            None => PortfolioCost::Zero,
            Some((range, unit)) => PortfolioCost::Known { range, unit },
        }
    }

    /// Total money demand (explicitly incomparable across currencies).
    #[must_use]
    pub fn total_money(&self) -> PortfolioCost {
        self.total(|b| &b.money)
    }

    /// Total compute demand.
    #[must_use]
    pub fn total_compute(&self) -> PortfolioCost {
        self.total(|b| &b.compute)
    }

    /// Total memory demand.
    #[must_use]
    pub fn total_memory(&self) -> PortfolioCost {
        self.total(|b| &b.memory)
    }

    /// Total calendar lead-time demand.
    #[must_use]
    pub fn total_lead_time(&self) -> PortfolioCost {
        self.total(|b| &b.lead_time)
    }
}

/// Proof that one proposal was actually executed, with its outcome color.
/// This is the ONLY object that can carry evidence color for an action.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExecutionReceipt {
    proposal: ActionId,
    outcome: ColorRank,
    executed_at_clock: BoundedId,
    executed_at_tick: i64,
    executor: BoundedId,
}

impl ExecutionReceipt {
    /// Admit an execution receipt against the exact proposal it executes.
    /// The outcome can never exceed the proposal's ceiling.
    pub fn admit(
        proposal: &ActionProposal,
        claimed_proposal_id: ActionId,
        outcome: ColorRank,
        executed_at: Expiry,
        executor: BoundedId,
    ) -> Result<Self, ActionError> {
        if proposal.content_id() != claimed_proposal_id {
            return refuse(
                ActionRule::ReceiptProposalMismatch,
                "receipt names a different proposal identity",
            );
        }
        if outcome > proposal.ceiling() {
            return refuse(
                ActionRule::CeilingExceeded,
                "execution outcome outranks the proposal's evidence ceiling",
            );
        }
        Ok(Self {
            proposal: claimed_proposal_id,
            outcome,
            executed_at_clock: executed_at.clock,
            executed_at_tick: executed_at.tick,
            executor,
        })
    }

    /// The executed proposal's identity.
    #[must_use]
    pub const fn proposal(&self) -> ActionId {
        self.proposal
    }

    /// The achieved outcome color (never above the ceiling).
    #[must_use]
    pub const fn outcome(&self) -> ColorRank {
        self.outcome
    }

    /// The executing party.
    #[must_use]
    pub fn executor(&self) -> &BoundedId {
        &self.executor
    }

    /// The clock the execution instant is stated on.
    #[must_use]
    pub fn executed_at_clock(&self) -> &BoundedId {
        &self.executed_at_clock
    }

    /// The execution instant tick.
    #[must_use]
    pub const fn executed_at_tick(&self) -> i64 {
        self.executed_at_tick
    }
}

/// A bounded canonical-transport refusal for evidence actions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ActionCodecError {
    offset: usize,
    detail: String,
}

impl ActionCodecError {
    fn at(offset: usize, detail: impl Into<String>) -> Self {
        Self {
            offset,
            detail: detail.into(),
        }
    }

    /// Stable rule identifier for every wire-level refusal.
    #[must_use]
    pub const fn rule_name(&self) -> &'static str {
        "action-canonical-identity"
    }

    /// Byte offset at which decoding refused.
    #[must_use]
    pub const fn offset(&self) -> usize {
        self.offset
    }

    /// Human-readable detail; never participates in rule matching.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for ActionCodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "action-canonical-identity at byte {}: {}",
            self.offset, self.detail
        )
    }
}

impl std::error::Error for ActionCodecError {}

fn u32_len(len: usize) -> u32 {
    u32::try_from(len).unwrap_or(u32::MAX)
}

fn push_str(out: &mut Vec<u8>, s: &str) {
    out.extend_from_slice(&u32_len(s.len()).to_be_bytes()[2..]);
    out.extend_from_slice(s.as_bytes());
}

fn encode_cost(out: &mut Vec<u8>, cost: &CostState) {
    match cost {
        CostState::Known { range, unit } => {
            out.push(1);
            out.extend_from_slice(&range.lo().to_bits().to_be_bytes());
            out.extend_from_slice(&range.hi().to_bits().to_be_bytes());
            match unit {
                CostUnit::Currency { code } => {
                    out.push(1);
                    out.extend_from_slice(code);
                }
                CostUnit::CpuCoreSeconds => out.push(2),
                CostUnit::WallSeconds => out.push(3),
                CostUnit::CalendarDays => out.push(4),
                CostUnit::MemoryBytes => out.push(5),
            }
        }
        CostState::Unknown { authority_gap } => {
            out.push(2);
            push_str(out, authority_gap.as_str());
        }
    }
}

fn encode_opt_id(out: &mut Vec<u8>, id: Option<&BoundedId>) {
    match id {
        None => out.push(0),
        Some(id) => {
            out.push(1);
            push_str(out, id.as_str());
        }
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], ActionCodecError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.bytes.len())
            .ok_or_else(|| ActionCodecError::at(self.pos, "truncated"))?;
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn u8(&mut self) -> Result<u8, ActionCodecError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, ActionCodecError> {
        let b = self.take(4)?;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn u64(&mut self) -> Result<u64, ActionCodecError> {
        let b = self.take(8)?;
        let mut arr = [0u8; 8];
        arr.copy_from_slice(b);
        Ok(u64::from_be_bytes(arr))
    }

    fn i64(&mut self) -> Result<i64, ActionCodecError> {
        Ok(self.u64()?.cast_signed())
    }

    fn len_u32(&mut self, budget: usize) -> Result<usize, ActionCodecError> {
        let b = self.take(4)?;
        let len = u32::from_be_bytes([b[0], b[1], b[2], b[3]]) as usize;
        if len > budget {
            return Err(ActionCodecError::at(self.pos, "collection above budget"));
        }
        Ok(len)
    }

    fn bounded_id(&mut self) -> Result<BoundedId, ActionCodecError> {
        let b = self.take(2)?;
        let len = usize::from(u16::from_be_bytes([b[0], b[1]]));
        let raw = self.take(len)?;
        let text = core::str::from_utf8(raw)
            .map_err(|_| ActionCodecError::at(self.pos, "identifier not UTF-8"))?;
        BoundedId::new(text).map_err(|e| ActionCodecError::at(self.pos, e.to_string()))
    }
}

fn decode_cost(r: &mut Reader<'_>) -> Result<CostState, ActionCodecError> {
    Ok(match r.u8()? {
        1 => {
            let lo = f64::from_bits(r.u64()?);
            let hi = f64::from_bits(r.u64()?);
            let range =
                NonNegRange::new(lo, hi).map_err(|e| ActionCodecError::at(r.pos, e.to_string()))?;
            let unit = match r.u8()? {
                1 => {
                    let code = r.take(3)?;
                    let text = core::str::from_utf8(code)
                        .map_err(|_| ActionCodecError::at(r.pos, "currency not UTF-8"))?
                        .to_owned();
                    CostUnit::currency(&text)
                        .map_err(|e| ActionCodecError::at(r.pos, e.to_string()))?
                }
                2 => CostUnit::CpuCoreSeconds,
                3 => CostUnit::WallSeconds,
                4 => CostUnit::CalendarDays,
                5 => CostUnit::MemoryBytes,
                other => {
                    return Err(ActionCodecError::at(r.pos, format!("bad unit tag {other}")));
                }
            };
            CostState::Known { range, unit }
        }
        2 => CostState::Unknown {
            authority_gap: r.bounded_id()?,
        },
        other => {
            return Err(ActionCodecError::at(r.pos, format!("bad cost tag {other}")));
        }
    })
}

fn decode_opt_id(r: &mut Reader<'_>) -> Result<Option<BoundedId>, ActionCodecError> {
    match r.u8()? {
        0 => Ok(None),
        1 => Ok(Some(r.bounded_id()?)),
        other => Err(ActionCodecError::at(
            r.pos,
            format!("bad option tag {other}"),
        )),
    }
}
