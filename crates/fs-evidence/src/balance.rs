//! I04.1 conservation-defect ontology: one typed accounting language for
//! balance defects.
//!
//! Cell residuals, port mismatches, remap errors, solver defects, jump work,
//! and unowned losses become comparable and composable only through one typed
//! receipt: [`BalanceDefectReceipt`]. The ontology pins quantity/species
//! semantics, extensive-versus-rate meaning, sign/orientation, spatial
//! chain/cochain support, time/window support, producer identity, the
//! accounting chart, expected closure, the measured defect, split
//! numerical/material/model uncertainty, the evidence color rank, and stable
//! content-addressed lineage.
//!
//! Exact zero, bounded defect, unknown, unowned remainder, and inapplicable
//! are DISTINCT states ([`DefectState`]) and never collapse into each other.
//! Receipts compose over disjoint spatial partitions and adjacent time
//! windows without double counting; every incompatibility refuses with a
//! named rule ([`BalanceRule`]).
//!
//! No-claims: this module accounts for defects it is HANDED — it does not
//! detect physical events, does not convert between unit scales, assumes
//! non-relativistic mass accounting, and adds split uncertainties linearly
//! (the conservative perfect-correlation bound), claiming no independence.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};

use crate::color::ColorRank;

/// Current canonical schema version for balance-defect receipts.
pub const BALANCE_SCHEMA_VERSION: u32 = 1;
/// Maximum UTF-8 bytes for one bounded identifier.
pub const MAX_BALANCE_ID_BYTES: usize = 96;
/// Maximum terms carried by one receipt.
pub const MAX_BALANCE_TERMS: usize = 4096;
/// Maximum spatial support ids carried by one receipt.
pub const MAX_BALANCE_SUPPORT_IDS: usize = 65_536;
/// Maximum accounts in one accounting chart.
pub const MAX_BALANCE_CHART_ACCOUNTS: usize = 256;
/// Maximum lineage parents recorded on one receipt.
pub const MAX_BALANCE_LINEAGE_PARENTS: usize = 64;
/// Maximum accepted canonical transport size in bytes.
pub const MAX_BALANCE_CANONICAL_BYTES: usize = 1024 * 1024;

const IDENTITY_DOMAIN: &str = "org.frankensim.fs-evidence.balance-defect.v1";
const CHART_DOMAIN: &str = "org.frankensim.fs-evidence.balance-chart.v1";
const MAGIC: &[u8; 4] = b"FSBD";

/// Stable named rule violated by a refusal. Every refusal names exactly one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum BalanceRule {
    /// Identifier empty, over budget, or outside `[a-z0-9._:-]`.
    IdBounds,
    /// A numeric field was NaN or infinite.
    NonFiniteValue,
    /// Interval lower endpoint exceeds upper endpoint.
    InvertedInterval,
    /// Negative zero attempted to alias canonical zero.
    SignedZeroAlias,
    /// A split uncertainty magnitude was negative.
    NegativeUncertainty,
    /// Element atomic number outside 1..=118.
    ElementRange,
    /// Momentum axis outside 0..=2 or chain/cochain degree above 15.
    SupportDegreeRange,
    /// Non-global support carried no ids, or global support carried ids.
    SupportShape,
    /// Window start does not precede window end.
    TimeOrder,
    /// Extensive amounts need a window; rates composed across windows.
    SemanticsMismatch,
    /// Receipts reference different charts (id, version, or digest).
    ChartMismatch,
    /// A term names an account absent from the chart.
    MissingChartAccount,
    /// The chart declares one account name twice.
    DuplicateChartAccount,
    /// One receipt carries two terms for the same (account, owner) pair.
    DuplicateTermOwner,
    /// Spatial supports overlap: composing them would double count.
    SupportOverlap,
    /// Supports live on different complexes or support kinds.
    SupportMismatch,
    /// Time supports reference different clocks.
    ClockMismatch,
    /// Windows are not adjacent (gap or overlap) for window composition.
    WindowsNotAdjacent,
    /// Receipts carry opposite sign conventions.
    SignConventionMismatch,
    /// Receipts account different quantities.
    QuantityMismatch,
    /// Receipts carry different decimal unit scales.
    ScaleMismatch,
    /// A term violates the quantity's admissible accounting rule.
    AccountingRuleViolation,
    /// Defect states cannot compose (inapplicable against anything else).
    DefectStateIncompatible,
    /// A collection exceeded its declared budget.
    CollectionBudget,
    /// Schema version unknown to this decoder.
    SchemaVersion,
}

impl BalanceRule {
    /// Stable machine-readable slug logged with every refusal.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::IdBounds => "balance-id-bounds",
            Self::NonFiniteValue => "balance-non-finite-value",
            Self::InvertedInterval => "balance-inverted-interval",
            Self::SignedZeroAlias => "balance-signed-zero-alias",
            Self::NegativeUncertainty => "balance-negative-uncertainty",
            Self::ElementRange => "balance-element-range",
            Self::SupportDegreeRange => "balance-support-degree-range",
            Self::SupportShape => "balance-support-shape",
            Self::TimeOrder => "balance-time-order",
            Self::SemanticsMismatch => "balance-semantics-mismatch",
            Self::ChartMismatch => "balance-chart-mismatch",
            Self::MissingChartAccount => "balance-missing-chart-account",
            Self::DuplicateChartAccount => "balance-duplicate-chart-account",
            Self::DuplicateTermOwner => "balance-duplicate-term-owner",
            Self::SupportOverlap => "balance-support-overlap",
            Self::SupportMismatch => "balance-support-mismatch",
            Self::ClockMismatch => "balance-clock-mismatch",
            Self::WindowsNotAdjacent => "balance-windows-not-adjacent",
            Self::SignConventionMismatch => "balance-sign-convention-mismatch",
            Self::QuantityMismatch => "balance-quantity-mismatch",
            Self::ScaleMismatch => "balance-scale-mismatch",
            Self::AccountingRuleViolation => "balance-accounting-rule-violation",
            Self::DefectStateIncompatible => "balance-defect-state-incompatible",
            Self::CollectionBudget => "balance-collection-budget",
            Self::SchemaVersion => "balance-schema-version",
        }
    }
}

/// A typed refusal naming the violated rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BalanceError {
    rule: BalanceRule,
    detail: String,
}

impl BalanceError {
    fn new(rule: BalanceRule, detail: impl Into<String>) -> Self {
        Self {
            rule,
            detail: detail.into(),
        }
    }

    /// The violated rule.
    #[must_use]
    pub const fn rule(&self) -> BalanceRule {
        self.rule
    }

    /// Human-readable detail; never participates in rule matching.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for BalanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule.slug(), self.detail)
    }
}

impl std::error::Error for BalanceError {}

fn refuse<T>(rule: BalanceRule, detail: impl Into<String>) -> Result<T, BalanceError> {
    Err(BalanceError::new(rule, detail))
}

/// A bounded lower-case identifier: 1..=96 bytes of `[a-z0-9._:-]`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BoundedId(String);

impl BoundedId {
    /// Admit a bounded identifier or refuse with [`BalanceRule::IdBounds`].
    pub fn new(value: &str) -> Result<Self, BalanceError> {
        if value.is_empty() || value.len() > MAX_BALANCE_ID_BYTES {
            return refuse(
                BalanceRule::IdBounds,
                format!("identifier length {} outside 1..=96", value.len()),
            );
        }
        if !value
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b"-_.:".contains(&b))
        {
            return refuse(
                BalanceRule::IdBounds,
                format!("identifier {value:?} outside [a-z0-9._:-]"),
            );
        }
        Ok(Self(value.to_owned()))
    }

    /// The admitted identifier text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Which conserved-or-produced quantity one receipt accounts.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum QuantityKind {
    /// Total energy, in the canonical unit joule.
    Energy,
    /// Exergy (available work), in joule.
    Exergy,
    /// Entropy, in joule per kelvin.
    Entropy,
    /// Electric charge, in coulomb.
    Charge,
    /// Total mass, in kilogram (non-relativistic no-claim applies).
    Mass,
    /// Amount of one chemical element, in mole.
    Element {
        /// Atomic number, 1..=118.
        atomic_number: u8,
    },
    /// Mass of one named species, in kilogram; reactions may produce or
    /// consume it, so its ledger is open.
    SpeciesMass {
        /// Species identity.
        species: BoundedId,
    },
    /// One linear-momentum component, in kilogram metre per second.
    Momentum {
        /// Frame identity the axis is resolved in.
        frame: BoundedId,
        /// Axis index, 0..=2.
        axis: u8,
    },
    /// One angular-momentum component, in kilogram metre squared per second.
    AngularMomentum {
        /// Frame identity the axis is resolved in.
        frame: BoundedId,
        /// Axis index, 0..=2.
        axis: u8,
    },
    /// Information, in bits; its ledger is open.
    Information,
}

/// The admissible accounting rule for one quantity kind. Different
/// quantities obey different balance laws; the ontology never forces one
/// generic equation beyond its semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AccountingRule {
    /// Strictly conserved: production/destruction terms are model-owned
    /// defect terms, not physics.
    Conserved,
    /// Production terms must be non-negative (second law: entropy).
    ProductionNonNegative,
    /// Destruction terms must be non-negative (exergy destruction).
    DestructionNonNegative,
    /// Open ledger: signed production admitted, but every term stays owned.
    OpenLedger,
}

impl QuantityKind {
    /// Validate the kind's own bounds.
    pub fn validate(&self) -> Result<(), BalanceError> {
        match self {
            Self::Element { atomic_number } => {
                if !(1..=118).contains(atomic_number) {
                    return refuse(
                        BalanceRule::ElementRange,
                        format!("atomic number {atomic_number} outside 1..=118"),
                    );
                }
            }
            Self::Momentum { axis, .. } | Self::AngularMomentum { axis, .. } if *axis > 2 => {
                return refuse(
                    BalanceRule::SupportDegreeRange,
                    format!("momentum axis {axis} outside 0..=2"),
                );
            }
            _ => {}
        }
        Ok(())
    }

    /// The admissible accounting rule for this quantity.
    #[must_use]
    pub fn accounting_rule(&self) -> AccountingRule {
        match self {
            Self::Energy
            | Self::Charge
            | Self::Mass
            | Self::Element { .. }
            | Self::Momentum { .. }
            | Self::AngularMomentum { .. } => AccountingRule::Conserved,
            Self::Entropy => AccountingRule::ProductionNonNegative,
            Self::Exergy => AccountingRule::DestructionNonNegative,
            Self::SpeciesMass { .. } | Self::Information => AccountingRule::OpenLedger,
        }
    }

    /// The canonical SI-coherent unit symbol every value is stated in
    /// (before the receipt's decimal scale).
    #[must_use]
    pub fn canonical_unit(&self) -> &'static str {
        match self {
            Self::Energy | Self::Exergy => "J",
            Self::Entropy => "J/K",
            Self::Charge => "C",
            Self::Mass | Self::SpeciesMass { .. } => "kg",
            Self::Element { .. } => "mol",
            Self::Momentum { .. } => "kg.m/s",
            Self::AngularMomentum { .. } => "kg.m^2/s",
            Self::Information => "bit",
        }
    }
}

/// Extensive-versus-rate meaning of every value in one receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Semantics {
    /// Amounts accumulated over the receipt's window.
    ExtensiveAmount,
    /// Flow rates (per unit time) on the receipt's time support.
    Rate,
}

/// Sign/orientation convention for boundary terms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignConvention {
    /// Flow into the support counts positive.
    InflowPositive,
    /// Flow out of the support counts positive.
    OutflowPositive,
}

impl SignConvention {
    /// The opposite orientation.
    #[must_use]
    pub const fn reversed(self) -> Self {
        match self {
            Self::InflowPositive => Self::OutflowPositive,
            Self::OutflowPositive => Self::InflowPositive,
        }
    }
}

/// The spatial support family one receipt lives on.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SupportKind {
    /// The whole named complex; carries no ids.
    Global,
    /// Named volumetric regions.
    Region,
    /// Named interfaces between regions.
    Interface,
    /// Mesh cells (chains) of one dimension.
    Cells {
        /// Cell dimension, 0..=15.
        dim: u8,
    },
    /// Cochain support of one degree.
    Cochain {
        /// Cochain degree, 0..=15.
        degree: u8,
    },
}

/// Spatial support: a bounded id-set on one named complex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpatialSupport {
    complex: BoundedId,
    kind: SupportKind,
    ids: BTreeSet<u64>,
}

impl SpatialSupport {
    /// Admit a spatial support or refuse with a named rule.
    pub fn new(
        complex: BoundedId,
        kind: SupportKind,
        ids: BTreeSet<u64>,
    ) -> Result<Self, BalanceError> {
        match kind {
            SupportKind::Global => {
                if !ids.is_empty() {
                    return refuse(
                        BalanceRule::SupportShape,
                        "global support must not carry ids",
                    );
                }
            }
            SupportKind::Cells { dim } => {
                if dim > 15 {
                    return refuse(
                        BalanceRule::SupportDegreeRange,
                        format!("cell dimension {dim} above 15"),
                    );
                }
                if ids.is_empty() {
                    return refuse(BalanceRule::SupportShape, "cell support needs ids");
                }
            }
            SupportKind::Cochain { degree } => {
                if degree > 15 {
                    return refuse(
                        BalanceRule::SupportDegreeRange,
                        format!("cochain degree {degree} above 15"),
                    );
                }
                if ids.is_empty() {
                    return refuse(BalanceRule::SupportShape, "cochain support needs ids");
                }
            }
            SupportKind::Region | SupportKind::Interface => {
                if ids.is_empty() {
                    return refuse(BalanceRule::SupportShape, "named support needs ids");
                }
            }
        }
        if ids.len() > MAX_BALANCE_SUPPORT_IDS {
            return refuse(
                BalanceRule::CollectionBudget,
                format!("{} support ids above budget", ids.len()),
            );
        }
        Ok(Self { complex, kind, ids })
    }

    /// The named complex.
    #[must_use]
    pub fn complex(&self) -> &BoundedId {
        &self.complex
    }

    /// The support family.
    #[must_use]
    pub fn kind(&self) -> &SupportKind {
        &self.kind
    }

    /// The admitted id-set.
    #[must_use]
    pub fn ids(&self) -> &BTreeSet<u64> {
        &self.ids
    }

    /// Union of two provably disjoint supports on the same complex and kind.
    /// Refuses double counting ([`BalanceRule::SupportOverlap`]) and
    /// cross-kind or cross-complex composition ([`BalanceRule::SupportMismatch`]).
    pub fn disjoint_union(&self, other: &Self) -> Result<Self, BalanceError> {
        if self.complex != other.complex || self.kind != other.kind {
            return refuse(
                BalanceRule::SupportMismatch,
                "supports live on different complexes or kinds",
            );
        }
        if matches!(self.kind, SupportKind::Global) {
            return refuse(
                BalanceRule::SupportOverlap,
                "global supports always overlap",
            );
        }
        if self.ids.intersection(&other.ids).next().is_some() {
            return refuse(
                BalanceRule::SupportOverlap,
                "support id-sets intersect: composition would double count",
            );
        }
        let ids: BTreeSet<u64> = self.ids.union(&other.ids).copied().collect();
        if ids.len() > MAX_BALANCE_SUPPORT_IDS {
            return refuse(BalanceRule::CollectionBudget, "union support above budget");
        }
        Ok(Self {
            complex: self.complex.clone(),
            kind: self.kind.clone(),
            ids,
        })
    }
}

/// Time support: one instant or one half-open window `[start, end)` on a
/// named logical clock.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TimeSupport {
    /// One instant tick.
    Instant {
        /// The logical clock identity.
        clock: BoundedId,
        /// The instant, in clock ticks.
        tick: i64,
    },
    /// One half-open window of clock ticks.
    Window {
        /// The logical clock identity.
        clock: BoundedId,
        /// Inclusive window start tick.
        start: i64,
        /// Exclusive window end tick.
        end: i64,
    },
}

impl TimeSupport {
    /// Validate ordering.
    pub fn validate(&self) -> Result<(), BalanceError> {
        if let Self::Window { start, end, .. } = self
            && start >= end
        {
            return refuse(
                BalanceRule::TimeOrder,
                format!("window start {start} not before end {end}"),
            );
        }
        Ok(())
    }

    /// The clock identity.
    #[must_use]
    pub fn clock(&self) -> &BoundedId {
        match self {
            Self::Instant { clock, .. } | Self::Window { clock, .. } => clock,
        }
    }
}

/// The declared role of one chart account.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AccountRole {
    /// Stored amount change inside the support.
    Storage,
    /// Flux across the support boundary.
    BoundaryFlux,
    /// Interior production.
    Production,
    /// Interior destruction.
    Destruction,
    /// Representation-transfer (remap) contribution.
    Remap,
    /// Solver-introduced defect.
    SolverDefect,
    /// Work done by jumps/impulses at events.
    JumpWork,
    /// Explicit transfer between named parties.
    Transfer,
}

const fn role_code(role: AccountRole) -> u8 {
    match role {
        AccountRole::Storage => 1,
        AccountRole::BoundaryFlux => 2,
        AccountRole::Production => 3,
        AccountRole::Destruction => 4,
        AccountRole::Remap => 5,
        AccountRole::SolverDefect => 6,
        AccountRole::JumpWork => 7,
        AccountRole::Transfer => 8,
    }
}

/// One versioned accounting chart: the closed set of admissible accounts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AccountingChart {
    id: BoundedId,
    version: u32,
    accounts: BTreeMap<BoundedId, AccountRole>,
}

impl AccountingChart {
    /// Admit a chart; duplicate account names refuse.
    pub fn new(
        id: BoundedId,
        version: u32,
        accounts: Vec<(BoundedId, AccountRole)>,
    ) -> Result<Self, BalanceError> {
        if accounts.is_empty() || accounts.len() > MAX_BALANCE_CHART_ACCOUNTS {
            return refuse(
                BalanceRule::CollectionBudget,
                format!("{} chart accounts outside 1..=256", accounts.len()),
            );
        }
        let mut map = BTreeMap::new();
        for (name, role) in accounts {
            if map.insert(name.clone(), role).is_some() {
                return refuse(
                    BalanceRule::DuplicateChartAccount,
                    format!("account {:?} declared twice", name.as_str()),
                );
            }
        }
        Ok(Self {
            id,
            version,
            accounts: map,
        })
    }

    /// The chart identity.
    #[must_use]
    pub fn id(&self) -> &BoundedId {
        &self.id
    }

    /// The chart version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// The declared role of one account, if present.
    #[must_use]
    pub fn role(&self, account: &BoundedId) -> Option<AccountRole> {
        self.accounts.get(account).copied()
    }

    /// Content digest binding the exact account set.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut bytes = Vec::new();
        push_str(&mut bytes, self.id.as_str());
        bytes.extend_from_slice(&self.version.to_be_bytes());
        bytes.extend_from_slice(&u32_len(self.accounts.len()).to_be_bytes());
        for (name, role) in &self.accounts {
            push_str(&mut bytes, name.as_str());
            bytes.push(role_code(*role));
        }
        hash_domain(CHART_DOMAIN, &bytes)
    }

    /// The exact reference a receipt pins.
    #[must_use]
    pub fn reference(&self) -> ChartRef {
        ChartRef {
            id: self.id.clone(),
            version: self.version,
            digest: self.digest(),
        }
    }
}

/// Exact reference to one accounting chart: id, version, and content digest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChartRef {
    /// Chart identity.
    pub id: BoundedId,
    /// Chart version.
    pub version: u32,
    /// Chart content digest.
    pub digest: ContentHash,
}

fn norm_zero(x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x }
}

/// A closed interval with exact `f64` endpoints. Negative zero endpoints
/// refuse: canonical transport admits exactly one zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    /// Admit an interval or refuse with a named rule.
    pub fn new(lo: f64, hi: f64) -> Result<Self, BalanceError> {
        if !lo.is_finite() || !hi.is_finite() {
            return refuse(BalanceRule::NonFiniteValue, "interval endpoint not finite");
        }
        if lo > hi {
            return refuse(
                BalanceRule::InvertedInterval,
                format!("interval [{lo:e}, {hi:e}] inverted"),
            );
        }
        if (lo == 0.0 && lo.is_sign_negative()) || (hi == 0.0 && hi.is_sign_negative()) {
            return refuse(
                BalanceRule::SignedZeroAlias,
                "negative zero endpoint aliases canonical zero",
            );
        }
        Ok(Self { lo, hi })
    }

    /// The exact degenerate point interval.
    pub fn point(value: f64) -> Result<Self, BalanceError> {
        Self::new(value, value)
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

    /// Conservative enclosure of the sum: endpoints add, then widen one ulp
    /// outward each side, so the result contains the exact real sum.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        let lo = norm_zero((self.lo + other.lo).next_down());
        let hi = norm_zero((self.hi + other.hi).next_up());
        Self { lo, hi }
    }

    /// Exact negation (endpoint swap); an involution bit-for-bit.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            lo: norm_zero(-self.hi),
            hi: norm_zero(-self.lo),
        }
    }

    /// Whether `other` lies entirely inside this interval.
    #[must_use]
    pub fn contains(&self, other: &Self) -> bool {
        self.lo <= other.lo && other.hi <= self.hi
    }
}

/// Split uncertainty magnitudes: numerical, material, and model components,
/// each a non-negative half-width in the receipt's unit and scale. They add
/// linearly (perfect-correlation conservative bound); no independence claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UncertaintySplit {
    numerical: f64,
    material: f64,
    model: f64,
}

impl UncertaintySplit {
    /// Admit an uncertainty split or refuse with a named rule.
    pub fn new(numerical: f64, material: f64, model: f64) -> Result<Self, BalanceError> {
        for (name, v) in [
            ("numerical", numerical),
            ("material", material),
            ("model", model),
        ] {
            if !v.is_finite() {
                return refuse(
                    BalanceRule::NonFiniteValue,
                    format!("{name} uncertainty not finite"),
                );
            }
            if v < 0.0 || (v == 0.0 && v.is_sign_negative()) {
                return refuse(
                    BalanceRule::NegativeUncertainty,
                    format!("{name} uncertainty negative"),
                );
            }
        }
        Ok(Self {
            numerical,
            material,
            model,
        })
    }

    /// Numerical component.
    #[must_use]
    pub const fn numerical(&self) -> f64 {
        self.numerical
    }

    /// Material component.
    #[must_use]
    pub const fn material(&self) -> f64 {
        self.material
    }

    /// Model component.
    #[must_use]
    pub const fn model(&self) -> f64 {
        self.model
    }

    /// Linear (conservative) componentwise sum, widened one ulp upward.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            numerical: (self.numerical + other.numerical).next_up(),
            material: (self.material + other.material).next_up(),
            model: (self.model + other.model).next_up(),
        }
    }
}

/// The measured defect. The five states are distinct and never collapse.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum DefectState {
    /// The balance closed exactly (bit-exact accounting identity).
    ExactZero,
    /// The defect is bounded inside this interval.
    Bounded(Interval),
    /// The defect could not be measured; visibility is preserved.
    Unknown,
    /// A remainder no party owns; it stays visible through composition.
    UnownedRemainder(Interval),
    /// This quantity's accounting does not apply here, for a named reason.
    Inapplicable {
        /// The named reason accounting does not apply.
        reason: BoundedId,
    },
}

impl DefectState {
    /// Compose two defect states over a partition or window union.
    /// Unknown absorbs measured states; an unowned remainder never launders
    /// into a plain bound; inapplicable composes only with the identically
    /// reasoned inapplicable.
    pub fn compose(&self, other: &Self) -> Result<Self, BalanceError> {
        match (self, other) {
            (Self::Inapplicable { reason: a }, Self::Inapplicable { reason: b }) => {
                if a == b {
                    Ok(Self::Inapplicable { reason: a.clone() })
                } else {
                    refuse(
                        BalanceRule::DefectStateIncompatible,
                        "inapplicable reasons differ",
                    )
                }
            }
            (Self::Inapplicable { .. }, _) | (_, Self::Inapplicable { .. }) => refuse(
                BalanceRule::DefectStateIncompatible,
                "inapplicable cannot compose with a measured state",
            ),
            (Self::Unknown, _) | (_, Self::Unknown) => Ok(Self::Unknown),
            (Self::UnownedRemainder(a), Self::UnownedRemainder(b)) => {
                Ok(Self::UnownedRemainder(a.add(b)))
            }
            (Self::UnownedRemainder(a), Self::Bounded(b))
            | (Self::Bounded(b), Self::UnownedRemainder(a)) => Ok(Self::UnownedRemainder(a.add(b))),
            (Self::UnownedRemainder(a), Self::ExactZero)
            | (Self::ExactZero, Self::UnownedRemainder(a)) => Ok(Self::UnownedRemainder(*a)),
            (Self::Bounded(a), Self::Bounded(b)) => Ok(Self::Bounded(a.add(b))),
            (Self::Bounded(a), Self::ExactZero) | (Self::ExactZero, Self::Bounded(a)) => {
                Ok(Self::Bounded(*a))
            }
            (Self::ExactZero, Self::ExactZero) => Ok(Self::ExactZero),
        }
    }

    /// Orientation reversal: intervals negate; the state family is preserved.
    #[must_use]
    pub fn reversed(&self) -> Self {
        match self {
            Self::Bounded(i) => Self::Bounded(i.neg()),
            Self::UnownedRemainder(i) => Self::UnownedRemainder(i.neg()),
            other => other.clone(),
        }
    }
}

/// One owned contribution to one chart account.
#[derive(Debug, Clone, PartialEq)]
pub struct BalanceTerm {
    /// The chart account this term posts to.
    pub account: BoundedId,
    /// The party owning this contribution.
    pub owner: BoundedId,
    /// The measured contribution interval.
    pub value: Interval,
    /// Evidence color rank of this contribution.
    pub color: ColorRank,
}

/// One composable conservation-defect receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct BalanceDefectReceipt {
    schema_version: u32,
    quantity: QuantityKind,
    semantics: Semantics,
    sign: SignConvention,
    scale_pow10: i8,
    chart: ChartRef,
    support: SpatialSupport,
    time: TimeSupport,
    producer: BoundedId,
    terms: Vec<BalanceTerm>,
    expected_closure: Interval,
    measured: DefectState,
    uncertainty: UncertaintySplit,
    lineage: Vec<ContentHash>,
}

/// Everything needed to admit one leaf receipt.
#[derive(Debug, Clone)]
pub struct BalanceDraft {
    /// Quantity accounted.
    pub quantity: QuantityKind,
    /// Extensive-versus-rate meaning.
    pub semantics: Semantics,
    /// Sign/orientation convention.
    pub sign: SignConvention,
    /// Decimal power-of-ten scale applied to the canonical unit.
    pub scale_pow10: i8,
    /// Exact chart reference.
    pub chart: ChartRef,
    /// Spatial support.
    pub support: SpatialSupport,
    /// Time support.
    pub time: TimeSupport,
    /// Producing party.
    pub producer: BoundedId,
    /// Owned term contributions.
    pub terms: Vec<BalanceTerm>,
    /// Expected closure interval (what a perfect balance would read).
    pub expected_closure: Interval,
    /// Measured defect state.
    pub measured: DefectState,
    /// Split uncertainty.
    pub uncertainty: UncertaintySplit,
}

impl BalanceDefectReceipt {
    /// Admit a leaf receipt against its chart, or refuse with a named rule.
    pub fn admit(draft: BalanceDraft, chart: &AccountingChart) -> Result<Self, BalanceError> {
        draft.quantity.validate()?;
        draft.time.validate()?;
        if chart.reference() != draft.chart {
            return refuse(
                BalanceRule::ChartMismatch,
                "draft chart reference does not match the supplied chart",
            );
        }
        if matches!(draft.semantics, Semantics::ExtensiveAmount)
            && matches!(draft.time, TimeSupport::Instant { .. })
        {
            return refuse(
                BalanceRule::SemanticsMismatch,
                "extensive amounts need a window, not an instant",
            );
        }
        if draft.terms.len() > MAX_BALANCE_TERMS {
            return refuse(
                BalanceRule::CollectionBudget,
                format!("{} terms above budget", draft.terms.len()),
            );
        }
        let rule = draft.quantity.accounting_rule();
        let mut seen: BTreeSet<(&str, &str)> = BTreeSet::new();
        for term in &draft.terms {
            let Some(role) = chart.role(&term.account) else {
                return refuse(
                    BalanceRule::MissingChartAccount,
                    format!("account {:?} absent from chart", term.account.as_str()),
                );
            };
            if !seen.insert((term.account.as_str(), term.owner.as_str())) {
                return refuse(
                    BalanceRule::DuplicateTermOwner,
                    format!(
                        "duplicate (account, owner) = ({:?}, {:?})",
                        term.account.as_str(),
                        term.owner.as_str()
                    ),
                );
            }
            match (rule, role) {
                (AccountingRule::ProductionNonNegative, AccountRole::Production)
                | (AccountingRule::DestructionNonNegative, AccountRole::Destruction)
                    if term.value.lo() < 0.0 =>
                {
                    return refuse(
                        BalanceRule::AccountingRuleViolation,
                        format!("{:?} admits no negative {:?} term", draft.quantity, role),
                    );
                }
                _ => {}
            }
        }
        let mut terms = draft.terms;
        terms.sort_by(|a, b| (&a.account, &a.owner).cmp(&(&b.account, &b.owner)));
        Ok(Self {
            schema_version: BALANCE_SCHEMA_VERSION,
            quantity: draft.quantity,
            semantics: draft.semantics,
            sign: draft.sign,
            scale_pow10: draft.scale_pow10,
            chart: draft.chart,
            support: draft.support,
            time: draft.time,
            producer: draft.producer,
            terms,
            expected_closure: draft.expected_closure,
            measured: draft.measured,
            uncertainty: draft.uncertainty,
            lineage: Vec::new(),
        })
    }

    /// Schema version the receipt was admitted under.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Quantity accounted.
    #[must_use]
    pub fn quantity(&self) -> &QuantityKind {
        &self.quantity
    }

    /// Extensive-versus-rate meaning.
    #[must_use]
    pub const fn semantics(&self) -> Semantics {
        self.semantics
    }

    /// Sign convention.
    #[must_use]
    pub const fn sign(&self) -> SignConvention {
        self.sign
    }

    /// Decimal scale exponent.
    #[must_use]
    pub const fn scale_pow10(&self) -> i8 {
        self.scale_pow10
    }

    /// Exact chart reference.
    #[must_use]
    pub fn chart(&self) -> &ChartRef {
        &self.chart
    }

    /// Spatial support.
    #[must_use]
    pub fn support(&self) -> &SpatialSupport {
        &self.support
    }

    /// Time support.
    #[must_use]
    pub fn time(&self) -> &TimeSupport {
        &self.time
    }

    /// Producing party.
    #[must_use]
    pub fn producer(&self) -> &BoundedId {
        &self.producer
    }

    /// Owned terms, canonically sorted by (account, owner).
    #[must_use]
    pub fn terms(&self) -> &[BalanceTerm] {
        &self.terms
    }

    /// Expected closure interval.
    #[must_use]
    pub fn expected_closure(&self) -> &Interval {
        &self.expected_closure
    }

    /// Measured defect state.
    #[must_use]
    pub fn measured(&self) -> &DefectState {
        &self.measured
    }

    /// Split uncertainty.
    #[must_use]
    pub fn uncertainty(&self) -> &UncertaintySplit {
        &self.uncertainty
    }

    /// Content ids of the parent receipts this one composed from.
    #[must_use]
    pub fn lineage(&self) -> &[ContentHash] {
        &self.lineage
    }

    /// Weakest evidence rank over all terms (no-laundering law); a receipt
    /// with no terms carries only Estimated authority.
    #[must_use]
    pub fn evidence_rank(&self) -> ColorRank {
        self.terms
            .iter()
            .map(|t| t.color)
            .min()
            .unwrap_or(ColorRank::Estimated)
    }

    /// Reverse the sign/orientation convention. Terms, closure, and the
    /// measured state negate; double reversal restores the receipt
    /// bit-for-bit.
    #[must_use]
    pub fn reverse_orientation(&self) -> Self {
        let mut out = self.clone();
        out.sign = self.sign.reversed();
        for term in &mut out.terms {
            term.value = term.value.neg();
        }
        out.expected_closure = self.expected_closure.neg();
        out.measured = self.measured.reversed();
        out
    }

    fn check_compatible(&self, other: &Self) -> Result<(), BalanceError> {
        if self.schema_version != other.schema_version {
            return refuse(BalanceRule::SchemaVersion, "schema versions differ");
        }
        if self.quantity != other.quantity {
            return refuse(BalanceRule::QuantityMismatch, "quantities differ");
        }
        if self.semantics != other.semantics {
            return refuse(BalanceRule::SemanticsMismatch, "semantics differ");
        }
        if self.sign != other.sign {
            return refuse(
                BalanceRule::SignConventionMismatch,
                "sign conventions differ",
            );
        }
        if self.scale_pow10 != other.scale_pow10 {
            return refuse(BalanceRule::ScaleMismatch, "decimal scales differ");
        }
        if self.chart != other.chart {
            return refuse(BalanceRule::ChartMismatch, "chart references differ");
        }
        if self.time.clock() != other.time.clock() {
            return refuse(BalanceRule::ClockMismatch, "clocks differ");
        }
        Ok(())
    }

    fn merge_terms(a: &[BalanceTerm], b: &[BalanceTerm]) -> Result<Vec<BalanceTerm>, BalanceError> {
        let mut merged: BTreeMap<(BoundedId, BoundedId), (Interval, ColorRank)> = BTreeMap::new();
        for term in a.iter().chain(b) {
            let key = (term.account.clone(), term.owner.clone());
            match merged.get_mut(&key) {
                None => {
                    merged.insert(key, (term.value, term.color));
                }
                Some((value, color)) => {
                    *value = value.add(&term.value);
                    *color = (*color).min(term.color);
                }
            }
        }
        if merged.len() > MAX_BALANCE_TERMS {
            return refuse(BalanceRule::CollectionBudget, "merged terms above budget");
        }
        Ok(merged
            .into_iter()
            .map(|((account, owner), (value, color))| BalanceTerm {
                account,
                owner,
                value,
                color,
            })
            .collect())
    }

    fn composed(
        &self,
        other: &Self,
        support: SpatialSupport,
        time: TimeSupport,
        composer: BoundedId,
    ) -> Result<Self, BalanceError> {
        let terms = Self::merge_terms(&self.terms, &other.terms)?;
        let measured = self.measured.compose(&other.measured)?;
        let lineage = vec![self.content_id(), other.content_id()];
        if lineage.len() > MAX_BALANCE_LINEAGE_PARENTS {
            return refuse(BalanceRule::CollectionBudget, "lineage above budget");
        }
        Ok(Self {
            schema_version: self.schema_version,
            quantity: self.quantity.clone(),
            semantics: self.semantics,
            sign: self.sign,
            scale_pow10: self.scale_pow10,
            chart: self.chart.clone(),
            support,
            time,
            producer: composer,
            terms,
            expected_closure: self.expected_closure.add(&other.expected_closure),
            measured,
            uncertainty: self.uncertainty.add(&other.uncertainty),
            lineage,
        })
    }

    /// Compose two receipts over a disjoint spatial partition. Both must
    /// share quantity, semantics, sign, scale, chart, clock, and the exact
    /// same time support; supports must be provably disjoint.
    pub fn compose_partition(
        &self,
        other: &Self,
        composer: BoundedId,
    ) -> Result<Self, BalanceError> {
        self.check_compatible(other)?;
        if self.time != other.time {
            return refuse(
                BalanceRule::WindowsNotAdjacent,
                "partition composition needs the exact same time support",
            );
        }
        let support = self.support.disjoint_union(&other.support)?;
        self.composed(other, support, self.time.clone(), composer)
    }

    /// Compose two extensive receipts over adjacent windows on the same
    /// support. Rates refuse: averaging rates across windows is not an
    /// accounting identity.
    pub fn compose_windows(&self, other: &Self, composer: BoundedId) -> Result<Self, BalanceError> {
        self.check_compatible(other)?;
        if self.semantics != Semantics::ExtensiveAmount {
            return refuse(
                BalanceRule::SemanticsMismatch,
                "only extensive amounts compose across windows",
            );
        }
        if self.support != other.support {
            return refuse(
                BalanceRule::SupportMismatch,
                "window composition needs the exact same support",
            );
        }
        let (
            TimeSupport::Window {
                clock,
                start: a_start,
                end: a_end,
            },
            TimeSupport::Window {
                start: b_start,
                end: b_end,
                ..
            },
        ) = (&self.time, &other.time)
        else {
            return refuse(
                BalanceRule::SemanticsMismatch,
                "window composition needs window supports",
            );
        };
        if a_end != b_start {
            return refuse(
                BalanceRule::WindowsNotAdjacent,
                format!("window end {a_end} does not meet next start {b_start}"),
            );
        }
        let time = TimeSupport::Window {
            clock: clock.clone(),
            start: *a_start,
            end: *b_end,
        };
        self.composed(other, self.support.clone(), time, composer)
    }

    /// Canonical transport bytes (versioned, bounded, deterministic).
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&self.schema_version.to_be_bytes());
        encode_quantity(&mut out, &self.quantity);
        out.push(match self.semantics {
            Semantics::ExtensiveAmount => 1,
            Semantics::Rate => 2,
        });
        out.push(match self.sign {
            SignConvention::InflowPositive => 1,
            SignConvention::OutflowPositive => 2,
        });
        out.push(self.scale_pow10.cast_unsigned());
        push_str(&mut out, self.chart.id.as_str());
        out.extend_from_slice(&self.chart.version.to_be_bytes());
        out.extend_from_slice(self.chart.digest.as_bytes());
        push_str(&mut out, self.support.complex.as_str());
        match &self.support.kind {
            SupportKind::Global => out.push(1),
            SupportKind::Region => out.push(2),
            SupportKind::Interface => out.push(3),
            SupportKind::Cells { dim } => {
                out.push(4);
                out.push(*dim);
            }
            SupportKind::Cochain { degree } => {
                out.push(5);
                out.push(*degree);
            }
        }
        out.extend_from_slice(&u32_len(self.support.ids.len()).to_be_bytes());
        for id in &self.support.ids {
            out.extend_from_slice(&id.to_be_bytes());
        }
        match &self.time {
            TimeSupport::Instant { clock, tick } => {
                out.push(1);
                push_str(&mut out, clock.as_str());
                out.extend_from_slice(&tick.to_be_bytes());
            }
            TimeSupport::Window { clock, start, end } => {
                out.push(2);
                push_str(&mut out, clock.as_str());
                out.extend_from_slice(&start.to_be_bytes());
                out.extend_from_slice(&end.to_be_bytes());
            }
        }
        push_str(&mut out, self.producer.as_str());
        out.extend_from_slice(&u32_len(self.terms.len()).to_be_bytes());
        for term in &self.terms {
            push_str(&mut out, term.account.as_str());
            push_str(&mut out, term.owner.as_str());
            push_interval(&mut out, &term.value);
            out.push(match term.color {
                ColorRank::Estimated => 1,
                ColorRank::Validated => 2,
                ColorRank::Verified => 3,
            });
        }
        push_interval(&mut out, &self.expected_closure);
        match &self.measured {
            DefectState::ExactZero => out.push(1),
            DefectState::Bounded(i) => {
                out.push(2);
                push_interval(&mut out, i);
            }
            DefectState::Unknown => out.push(3),
            DefectState::UnownedRemainder(i) => {
                out.push(4);
                push_interval(&mut out, i);
            }
            DefectState::Inapplicable { reason } => {
                out.push(5);
                push_str(&mut out, reason.as_str());
            }
        }
        out.extend_from_slice(&self.uncertainty.numerical.to_bits().to_be_bytes());
        out.extend_from_slice(&self.uncertainty.material.to_bits().to_be_bytes());
        out.extend_from_slice(&self.uncertainty.model.to_bits().to_be_bytes());
        out.extend_from_slice(&u32_len(self.lineage.len()).to_be_bytes());
        for parent in &self.lineage {
            out.extend_from_slice(parent.as_bytes());
        }
        out
    }

    /// Domain-separated content identity over the canonical bytes.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        hash_domain(IDENTITY_DOMAIN, &self.canonical_bytes())
    }

    /// Decode canonical transport bytes. Decoding is fail-closed: bounded,
    /// version-checked, exhaustively validated, and canonical (re-encoding
    /// must reproduce the input bit-for-bit, so trailing bytes, unsorted
    /// collections, and aliased encodings all refuse).
    // The decoder deliberately mirrors the complete canonical field order in
    // one place so reviewers can compare it directly with `canonical_bytes`.
    #[allow(clippy::too_many_lines)]
    pub fn decode(bytes: &[u8]) -> Result<Self, BalanceCodecError> {
        if bytes.len() > MAX_BALANCE_CANONICAL_BYTES {
            return Err(BalanceCodecError::at(0, "transport above size budget"));
        }
        let mut r = Reader { bytes, pos: 0 };
        let magic = r.take(4)?;
        if magic != MAGIC {
            return Err(BalanceCodecError::at(0, "bad magic"));
        }
        let schema_version = r.u32()?;
        if schema_version != BALANCE_SCHEMA_VERSION {
            return Err(BalanceCodecError::at(
                r.pos,
                format!("unknown schema version {schema_version}"),
            ));
        }
        let quantity = decode_quantity(&mut r)?;
        let semantics = match r.u8()? {
            1 => Semantics::ExtensiveAmount,
            2 => Semantics::Rate,
            other => {
                return Err(BalanceCodecError::at(
                    r.pos,
                    format!("bad semantics tag {other}"),
                ));
            }
        };
        let sign = match r.u8()? {
            1 => SignConvention::InflowPositive,
            2 => SignConvention::OutflowPositive,
            other => Err(BalanceCodecError::at(
                r.pos,
                format!("bad sign tag {other}"),
            ))?,
        };
        let scale_pow10 = r.u8()?.cast_signed();
        let chart_id = r.bounded_id()?;
        let chart_version = r.u32()?;
        let chart_digest = ContentHash::from_slice(r.take(32)?)
            .ok_or_else(|| BalanceCodecError::at(r.pos, "bad chart digest"))?;
        let complex = r.bounded_id()?;
        let kind = match r.u8()? {
            1 => SupportKind::Global,
            2 => SupportKind::Region,
            3 => SupportKind::Interface,
            4 => SupportKind::Cells { dim: r.u8()? },
            5 => SupportKind::Cochain { degree: r.u8()? },
            other => {
                return Err(BalanceCodecError::at(
                    r.pos,
                    format!("bad support tag {other}"),
                ));
            }
        };
        let id_count = r.len_u32(MAX_BALANCE_SUPPORT_IDS)?;
        let mut ids = BTreeSet::new();
        let mut prev: Option<u64> = None;
        for _ in 0..id_count {
            let id = r.u64()?;
            if prev.is_some_and(|p| p >= id) {
                return Err(BalanceCodecError::at(
                    r.pos,
                    "support ids not strictly sorted",
                ));
            }
            prev = Some(id);
            ids.insert(id);
        }
        let support = SpatialSupport::new(complex, kind, ids)
            .map_err(|e| BalanceCodecError::at(r.pos, e.to_string()))?;
        let time = match r.u8()? {
            1 => TimeSupport::Instant {
                clock: r.bounded_id()?,
                tick: r.i64()?,
            },
            2 => TimeSupport::Window {
                clock: r.bounded_id()?,
                start: r.i64()?,
                end: r.i64()?,
            },
            other => {
                return Err(BalanceCodecError::at(
                    r.pos,
                    format!("bad time tag {other}"),
                ));
            }
        };
        time.validate()
            .map_err(|e| BalanceCodecError::at(r.pos, e.to_string()))?;
        let producer = r.bounded_id()?;
        let term_count = r.len_u32(MAX_BALANCE_TERMS)?;
        let mut terms = Vec::with_capacity(term_count);
        for _ in 0..term_count {
            let account = r.bounded_id()?;
            let owner = r.bounded_id()?;
            let value = r.interval()?;
            let color = match r.u8()? {
                1 => ColorRank::Estimated,
                2 => ColorRank::Validated,
                3 => ColorRank::Verified,
                other => {
                    return Err(BalanceCodecError::at(
                        r.pos,
                        format!("bad color tag {other}"),
                    ));
                }
            };
            terms.push(BalanceTerm {
                account,
                owner,
                value,
                color,
            });
        }
        let expected_closure = r.interval()?;
        let measured = match r.u8()? {
            1 => DefectState::ExactZero,
            2 => DefectState::Bounded(r.interval()?),
            3 => DefectState::Unknown,
            4 => DefectState::UnownedRemainder(r.interval()?),
            5 => DefectState::Inapplicable {
                reason: r.bounded_id()?,
            },
            other => {
                return Err(BalanceCodecError::at(
                    r.pos,
                    format!("bad defect tag {other}"),
                ));
            }
        };
        let uncertainty = UncertaintySplit::new(
            f64::from_bits(r.u64()?),
            f64::from_bits(r.u64()?),
            f64::from_bits(r.u64()?),
        )
        .map_err(|e| BalanceCodecError::at(r.pos, e.to_string()))?;
        let lineage_count = r.len_u32(MAX_BALANCE_LINEAGE_PARENTS)?;
        let mut lineage = Vec::with_capacity(lineage_count);
        for _ in 0..lineage_count {
            lineage.push(
                ContentHash::from_slice(r.take(32)?)
                    .ok_or_else(|| BalanceCodecError::at(r.pos, "bad lineage hash"))?,
            );
        }
        if r.pos != bytes.len() {
            return Err(BalanceCodecError::at(r.pos, "trailing bytes"));
        }
        let receipt = Self {
            schema_version,
            quantity,
            semantics,
            sign,
            scale_pow10,
            chart: ChartRef {
                id: chart_id,
                version: chart_version,
                digest: chart_digest,
            },
            support,
            time,
            producer,
            terms,
            expected_closure,
            measured,
            uncertainty,
            lineage,
        };
        receipt
            .quantity
            .validate()
            .map_err(|e| BalanceCodecError::at(0, e.to_string()))?;
        if receipt.canonical_bytes() != bytes {
            return Err(BalanceCodecError::at(0, "non-canonical encoding"));
        }
        Ok(receipt)
    }
}

/// A bounded canonical-transport refusal for balance receipts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BalanceCodecError {
    offset: usize,
    detail: String,
}

impl BalanceCodecError {
    fn at(offset: usize, detail: impl Into<String>) -> Self {
        Self {
            offset,
            detail: detail.into(),
        }
    }

    /// Stable rule identifier for every wire-level refusal.
    #[must_use]
    pub const fn rule_name(&self) -> &'static str {
        "balance-canonical-identity"
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

impl fmt::Display for BalanceCodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "balance-canonical-identity at byte {}: {}",
            self.offset, self.detail
        )
    }
}

impl std::error::Error for BalanceCodecError {}

fn u32_len(len: usize) -> u32 {
    u32::try_from(len).unwrap_or(u32::MAX)
}

fn push_str(out: &mut Vec<u8>, s: &str) {
    out.extend_from_slice(&u32_len(s.len()).to_be_bytes()[2..]);
    out.extend_from_slice(s.as_bytes());
}

fn push_interval(out: &mut Vec<u8>, i: &Interval) {
    out.extend_from_slice(&i.lo().to_bits().to_be_bytes());
    out.extend_from_slice(&i.hi().to_bits().to_be_bytes());
}

fn encode_quantity(out: &mut Vec<u8>, q: &QuantityKind) {
    match q {
        QuantityKind::Energy => out.push(1),
        QuantityKind::Exergy => out.push(2),
        QuantityKind::Entropy => out.push(3),
        QuantityKind::Charge => out.push(4),
        QuantityKind::Mass => out.push(5),
        QuantityKind::Element { atomic_number } => {
            out.push(6);
            out.push(*atomic_number);
        }
        QuantityKind::SpeciesMass { species } => {
            out.push(7);
            push_str(out, species.as_str());
        }
        QuantityKind::Momentum { frame, axis } => {
            out.push(8);
            push_str(out, frame.as_str());
            out.push(*axis);
        }
        QuantityKind::AngularMomentum { frame, axis } => {
            out.push(9);
            push_str(out, frame.as_str());
            out.push(*axis);
        }
        QuantityKind::Information => out.push(10),
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], BalanceCodecError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.bytes.len())
            .ok_or_else(|| BalanceCodecError::at(self.pos, "truncated"))?;
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn u8(&mut self) -> Result<u8, BalanceCodecError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, BalanceCodecError> {
        let b = self.take(4)?;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn u64(&mut self) -> Result<u64, BalanceCodecError> {
        let b = self.take(8)?;
        let mut arr = [0u8; 8];
        arr.copy_from_slice(b);
        Ok(u64::from_be_bytes(arr))
    }

    fn i64(&mut self) -> Result<i64, BalanceCodecError> {
        Ok(self.u64()?.cast_signed())
    }

    fn len_u32(&mut self, budget: usize) -> Result<usize, BalanceCodecError> {
        let b = self.take(4)?;
        let len = u32::from_be_bytes([b[0], b[1], b[2], b[3]]) as usize;
        if len > budget {
            return Err(BalanceCodecError::at(self.pos, "collection above budget"));
        }
        Ok(len)
    }

    fn bounded_id(&mut self) -> Result<BoundedId, BalanceCodecError> {
        let b = self.take(2)?;
        let len = usize::from(u16::from_be_bytes([b[0], b[1]]));
        if len > MAX_BALANCE_ID_BYTES {
            return Err(BalanceCodecError::at(self.pos, "identifier above budget"));
        }
        let raw = self.take(len)?;
        let text = core::str::from_utf8(raw)
            .map_err(|_| BalanceCodecError::at(self.pos, "identifier not UTF-8"))?;
        BoundedId::new(text).map_err(|e| BalanceCodecError::at(self.pos, e.to_string()))
    }

    fn interval(&mut self) -> Result<Interval, BalanceCodecError> {
        let lo = f64::from_bits(self.u64()?);
        let hi = f64::from_bits(self.u64()?);
        Interval::new(lo, hi).map_err(|e| BalanceCodecError::at(self.pos, e.to_string()))
    }
}

fn decode_quantity(r: &mut Reader<'_>) -> Result<QuantityKind, BalanceCodecError> {
    Ok(match r.u8()? {
        1 => QuantityKind::Energy,
        2 => QuantityKind::Exergy,
        3 => QuantityKind::Entropy,
        4 => QuantityKind::Charge,
        5 => QuantityKind::Mass,
        6 => QuantityKind::Element {
            atomic_number: r.u8()?,
        },
        7 => QuantityKind::SpeciesMass {
            species: r.bounded_id()?,
        },
        8 => QuantityKind::Momentum {
            frame: r.bounded_id()?,
            axis: r.u8()?,
        },
        9 => QuantityKind::AngularMomentum {
            frame: r.bounded_id()?,
            axis: r.u8()?,
        },
        10 => QuantityKind::Information,
        other => {
            return Err(BalanceCodecError::at(
                r.pos,
                format!("bad quantity tag {other}"),
            ));
        }
    })
}
