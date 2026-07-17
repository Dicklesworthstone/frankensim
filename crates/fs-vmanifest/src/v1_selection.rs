//! V.1.5 deterministic selection semantics (bead
//! frankensim-leapfrog-2026-program-i94v.7.1.5): claim strata, campaign
//! profiles, filters, budgets, capabilities, and shard assignment over
//! the stable v1 identities — expanded DETERMINISTICALLY, receipted
//! BEFORE execution, and diffed exactly.
//!
//! Two orthogonal, versioned axes: [`Stratum`] {core, max} selects the
//! declared claim/capability surface; a [`ProfileId`] selects campaign
//! intensity and evidence obligations from the eight atomic built-ins or
//! a versioned manifest-defined composite whose ordered inputs,
//! precedence, and digest are frozen before execution. Repeated profile
//! flags and implicit string composition ("smoke+soak") are INVALID;
//! legacy SMOKE/MID/FULL names refuse here with a ranked migration to the
//! V.4.6 adapters — they can never silently alias a stratum/profile pair.
//!
//! [`expand_selection`] is a pure function of (manifest cases, selection
//! input): the same inputs select the same logical cases independent of
//! enumeration order, worker count, or presentation; capability-blocked
//! and predicate-skipped cases stay VISIBLE with named reasons; shards
//! partition the selection exactly; an empty selection is explicitly
//! non-green. [`semantic_diff`] turns any input change into an exact
//! added/removed/changed report.
//!
//! No-claims: selection is metadata — expanding it runs no production
//! computation and adjudicates nothing; timeout/cancellation policy and
//! budgets are carried and receipted, enforced by the executing campaign
//! runner; the legacy adapters themselves are V.4.6 scope.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};

use crate::v1::{CaseId, V1Error};

const SELECTION_DOMAIN: &str = "org.frankensim.fs-vmanifest.selection-receipt.v1";

/// The claim/capability stratum: scientific scope, orthogonal to
/// campaign intensity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Stratum {
    /// The core declared claim surface.
    Core,
    /// The maximal declared claim surface (a superset of core).
    Max,
}

impl Stratum {
    /// Stable lowercase name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Core => "core",
            Self::Max => "max",
        }
    }
}

/// The eight atomic built-in campaign profiles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum BuiltinProfile {
    /// Fast confidence pass.
    Smoke,
    /// The standard battery.
    Standard,
    /// Adversarial/falsification emphasis.
    Adversarial,
    /// Long-duration soak.
    Soak,
    /// Security emphasis.
    Security,
    /// Chaos/fault injection.
    Chaos,
    /// Cross-ISA reproduction.
    CrossIsa,
    /// Release qualification.
    Release,
}

impl BuiltinProfile {
    /// Stable lowercase id.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Smoke => "smoke",
            Self::Standard => "standard",
            Self::Adversarial => "adversarial",
            Self::Soak => "soak",
            Self::Security => "security",
            Self::Chaos => "chaos",
            Self::CrossIsa => "cross-isa",
            Self::Release => "release",
        }
    }

    fn parse(value: &str) -> Option<BuiltinProfile> {
        Some(match value {
            "smoke" => Self::Smoke,
            "standard" => Self::Standard,
            "adversarial" => Self::Adversarial,
            "soak" => Self::Soak,
            "security" => Self::Security,
            "chaos" => Self::Chaos,
            "cross-isa" => Self::CrossIsa,
            "release" => Self::Release,
            _ => return None,
        })
    }
}

/// A versioned manifest-defined composite profile: ordered inputs,
/// explicit precedence rule, frozen before execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositeProfile {
    /// The manifest-scoped composite id.
    pub id: String,
    /// The composite's frozen version.
    pub version: u32,
    /// Ordered atomic inputs (order IS precedence for conflicts).
    pub inputs: Vec<BuiltinProfile>,
    /// The stated precedence/conflict rule.
    pub precedence_rule: String,
}

impl CompositeProfile {
    /// Validate: at least two distinct ordered inputs, a stated rule.
    pub fn validate(&self) -> Result<(), V1Error> {
        if self.inputs.len() < 2 {
            return Err(v1err(
                "v1-profile-composition",
                "a composite profile needs at least two ordered inputs",
                &["use the atomic built-in profile directly"],
            ));
        }
        let mut seen = BTreeSet::new();
        for input in &self.inputs {
            if !seen.insert(*input) {
                return Err(v1err(
                    "v1-profile-composition",
                    format!("composite input {:?} repeated", input.name()),
                    &["each atomic input appears once; order encodes precedence"],
                ));
            }
        }
        if self.precedence_rule.is_empty() {
            return Err(v1err(
                "v1-profile-composition",
                "composite precedence/conflict rule must be stated",
                &["freeze the precedence rule text before execution"],
            ));
        }
        Ok(())
    }

    /// The frozen semantic digest of the composite definition.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(self.id.as_bytes());
        bytes.push(0);
        bytes.extend_from_slice(&self.version.to_be_bytes());
        for input in &self.inputs {
            bytes.extend_from_slice(input.name().as_bytes());
            bytes.push(0);
        }
        bytes.extend_from_slice(self.precedence_rule.as_bytes());
        hash_domain("org.frankensim.fs-vmanifest.composite-profile.v1", &bytes)
    }
}

/// One atomic, manifest-resolved campaign profile selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProfileId {
    /// A built-in atomic profile.
    Builtin(BuiltinProfile),
    /// A frozen composite.
    Composite(CompositeProfile),
}

impl ProfileId {
    /// The stable rendered id.
    #[must_use]
    pub fn render(&self) -> String {
        match self {
            Self::Builtin(b) => b.name().to_owned(),
            Self::Composite(c) => format!("{}@v{}", c.id, c.version),
        }
    }

    /// Parse ONE atomic profile flag. Implicit composition, repetition,
    /// and legacy tier names refuse with ranked fixes.
    pub fn parse_atomic(value: &str) -> Result<ProfileId, V1Error> {
        if let Some(builtin) = BuiltinProfile::parse(value) {
            return Ok(ProfileId::Builtin(builtin));
        }
        let upper = value.to_ascii_uppercase();
        if ["SMOKE", "MID", "FULL"].contains(&upper.as_str()) && value != "smoke" {
            return Err(v1err(
                "v1-legacy-alias",
                format!("legacy tier name {value:?} cannot alias a stratum/profile pair"),
                &[
                    "route through the V.4.6 adapter, which must reconcile cases, claims, \
                     budgets, capabilities, acceptance rules, and artifacts or refuse",
                    "or select an explicit stratum plus one atomic profile",
                ],
            ));
        }
        if value.contains('+') || value.contains(',') {
            return Err(v1err(
                "v1-profile-composition",
                format!("implicit string composition {value:?} is invalid"),
                &[
                    "define a versioned manifest composite profile with ordered inputs and a \
                     frozen precedence rule",
                    "or expand into separate CampaignRuns with separately retained receipts",
                ],
            ));
        }
        Err(v1err(
            "v1-profile-unknown",
            format!("unknown profile {value:?}"),
            &["use one of smoke|standard|adversarial|soak|security|chaos|cross-isa|release"],
        ))
    }

    /// Parse a full `--profile` flag LIST. More than one flag refuses:
    /// a suite needing several profiles expands into separate runs.
    pub fn parse_flags(flags: &[&str]) -> Result<ProfileId, V1Error> {
        match flags {
            [] => Err(v1err(
                "v1-profile-missing",
                "a CampaignRun selects exactly one profile",
                &["pass exactly one --profile"],
            )),
            [one] => Self::parse_atomic(one),
            many => Err(v1err(
                "v1-profile-composition",
                format!("{} repeated --profile flags are invalid", many.len()),
                &[
                    "expand into separate CampaignRuns with separately retained receipts",
                    "or freeze a versioned composite profile in the manifest",
                ],
            )),
        }
    }
}

fn v1err(rule: &'static str, detail: impl Into<String>, fixes: &[&str]) -> V1Error {
    V1Error::with_fixes(rule, detail, fixes)
}

/// One selectable case row from the manifest: stable identity plus its
/// declared stratum, profile memberships, and capability demands.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectableCase {
    /// The stable case identity.
    pub case: CaseId,
    /// The stratum the case's claim belongs to.
    pub stratum: Stratum,
    /// The atomic profiles that include this case.
    pub profiles: BTreeSet<BuiltinProfile>,
    /// Capabilities the case requires (instruments, ISAs, hosts).
    pub required_capabilities: BTreeSet<String>,
}

/// A prefix filter over case ids, applied in declared order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Filter {
    /// Keep only cases whose id starts with the prefix.
    IncludePrefix(String),
    /// Drop cases whose id starts with the prefix.
    ExcludePrefix(String),
}

/// A named skip predicate: the NAME is receipted so nothing disappears
/// anonymously.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NamedSkip {
    /// The predicate's stable name.
    pub name: String,
    /// Case-id prefix the predicate suppresses.
    pub prefix: String,
    /// The stated reason.
    pub reason: String,
}

/// The complete selection input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionInput {
    /// The requested stratum (max includes core's surface).
    pub stratum: Stratum,
    /// Exactly one atomic, manifest-resolved profile.
    pub profile: ProfileId,
    /// Ordered filters.
    pub filters: Vec<Filter>,
    /// Available capabilities.
    pub capabilities: BTreeSet<String>,
    /// Named budgets (accuracy/time/memory/samples), echoed verbatim.
    pub budgets: BTreeMap<String, u64>,
    /// Named skip predicates.
    pub skips: Vec<NamedSkip>,
    /// Shard count (>= 1).
    pub shards: u32,
}

/// One visibly-skipped case: id plus the NAMED reason.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SkippedCase {
    /// The case.
    pub case: CaseId,
    /// The named reason ("capability:<missing>", "skip:<name>").
    pub reason: String,
}

/// The canonical pre-execution selection receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionReceipt {
    /// The requested stratum.
    pub stratum: Stratum,
    /// The rendered profile id.
    pub profile: String,
    /// Budgets echoed verbatim.
    pub budgets: BTreeMap<String, u64>,
    /// Selected case ids, ascending.
    pub selected: Vec<CaseId>,
    /// Visibly skipped cases with named reasons, ascending by case.
    pub skipped: Vec<SkippedCase>,
    /// Shard assignment: `shards[i]` lists shard i's cases, ascending.
    pub shards: Vec<Vec<CaseId>>,
}

impl SelectionReceipt {
    /// The deterministic selection digest (the semantic-diff anchor).
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(self.stratum.name().as_bytes());
        bytes.push(0);
        bytes.extend_from_slice(self.profile.as_bytes());
        bytes.push(0);
        for (name, value) in &self.budgets {
            bytes.extend_from_slice(name.as_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        for case in &self.selected {
            bytes.extend_from_slice(case.as_str().as_bytes());
            bytes.push(0);
        }
        for skip in &self.skipped {
            bytes.extend_from_slice(skip.case.as_str().as_bytes());
            bytes.push(1);
            bytes.extend_from_slice(skip.reason.as_bytes());
            bytes.push(0);
        }
        bytes.extend_from_slice(&(self.shards.len() as u32).to_be_bytes());
        for shard in &self.shards {
            bytes.extend_from_slice(&(shard.len() as u32).to_be_bytes());
            for case in shard {
                bytes.extend_from_slice(case.as_str().as_bytes());
                bytes.push(0);
            }
        }
        hash_domain(SELECTION_DOMAIN, &bytes)
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Deterministically expand a selection: pure in (cases, input), stable
/// under input reordering, capability- and predicate-skips visible by
/// name, shards an exact partition, empty selection explicitly refused.
pub fn expand_selection(
    cases: &[SelectableCase],
    input: &SelectionInput,
) -> Result<SelectionReceipt, V1Error> {
    if input.shards == 0 {
        return Err(v1err(
            "v1-selection-shards",
            "shard count must be at least one",
            &["pass shards >= 1"],
        ));
    }
    let profile_members: BTreeSet<BuiltinProfile> = match &input.profile {
        ProfileId::Builtin(b) => BTreeSet::from([*b]),
        ProfileId::Composite(c) => {
            c.validate()?;
            c.inputs.iter().copied().collect()
        }
    };

    // Canonical order: ascending case id, independent of input order.
    let mut rows: Vec<&SelectableCase> = cases.iter().collect();
    rows.sort_by(|a, b| a.case.cmp(&b.case));
    let mut ids = BTreeSet::new();
    for row in &rows {
        if !ids.insert(&row.case) {
            return Err(v1err(
                "v1-duplicate-case",
                format!("case id {:?} declared twice", row.case.as_str()),
                &["case identities are stable and unique per manifest"],
            ));
        }
    }

    let mut selected = Vec::new();
    let mut skipped = Vec::new();
    for row in rows {
        // Stratum: core selects core; max selects core + max.
        let stratum_admits = match input.stratum {
            Stratum::Core => row.stratum == Stratum::Core,
            Stratum::Max => true,
        };
        if !stratum_admits {
            continue; // Out of scientific scope: not selected, not "skipped".
        }
        if row.profiles.intersection(&profile_members).next().is_none() {
            continue; // Different campaign intensity: out of this run.
        }
        // Ordered filters: last matching filter wins; default include.
        let mut keep = true;
        for filter in &input.filters {
            match filter {
                Filter::IncludePrefix(p) => {
                    if row.case.as_str().starts_with(p.as_str()) {
                        keep = true;
                    }
                }
                Filter::ExcludePrefix(p) => {
                    if row.case.as_str().starts_with(p.as_str()) {
                        keep = false;
                    }
                }
            }
        }
        if !keep {
            continue;
        }
        // Capability routing: missing capability = VISIBLE named skip.
        if let Some(missing) = row
            .required_capabilities
            .difference(&input.capabilities)
            .next()
        {
            skipped.push(SkippedCase {
                case: row.case.clone(),
                reason: format!("capability:{missing}"),
            });
            continue;
        }
        // Named skip predicates.
        if let Some(skip) = input
            .skips
            .iter()
            .find(|s| row.case.as_str().starts_with(s.prefix.as_str()))
        {
            skipped.push(SkippedCase {
                case: row.case.clone(),
                reason: format!("skip:{}:{}", skip.name, skip.reason),
            });
            continue;
        }
        selected.push(row.case.clone());
    }

    if selected.is_empty() {
        return Err(v1err(
            "v1-empty-selection",
            "the selection is empty: an empty or unsupported selection is non-green",
            &[
                "widen the filters or supply the missing capabilities",
                "an intentionally empty run must be declared as such upstream, never inferred",
            ],
        ));
    }

    // Deterministic shard partition by case-id hash (stable across
    // worker counts and enumeration order).
    let mut shards: Vec<Vec<CaseId>> = vec![Vec::new(); input.shards as usize];
    for case in &selected {
        let shard = (fnv1a64(case.as_str().as_bytes()) % u64::from(input.shards)) as usize;
        shards[shard].push(case.clone());
    }

    Ok(SelectionReceipt {
        stratum: input.stratum,
        profile: input.profile.render(),
        budgets: input.budgets.clone(),
        selected,
        skipped,
        shards,
    })
}

/// The exact semantic diff between two selection receipts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionDiff {
    /// Cases in `after` but not `before`.
    pub added: Vec<CaseId>,
    /// Cases in `before` but not `after` — nothing disappears silently.
    pub removed: Vec<CaseId>,
    /// Budget rows that changed, as (name, before, after).
    pub budget_changes: Vec<(String, Option<u64>, Option<u64>)>,
    /// Whether stratum or profile changed (scientific scope / intensity).
    pub scope_changed: bool,
}

impl SelectionDiff {
    /// Whether the two receipts are semantically identical.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.added.is_empty()
            && self.removed.is_empty()
            && self.budget_changes.is_empty()
            && !self.scope_changed
    }
}

impl fmt::Display for SelectionDiff {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "selection diff: +{} -{} budgets~{} scope_changed={}",
            self.added.len(),
            self.removed.len(),
            self.budget_changes.len(),
            self.scope_changed
        )
    }
}

/// Compute the exact diff between two receipts.
#[must_use]
pub fn semantic_diff(before: &SelectionReceipt, after: &SelectionReceipt) -> SelectionDiff {
    let b: BTreeSet<&CaseId> = before.selected.iter().collect();
    let a: BTreeSet<&CaseId> = after.selected.iter().collect();
    let added = a.difference(&b).map(|c| (*c).clone()).collect();
    let removed = b.difference(&a).map(|c| (*c).clone()).collect();
    let mut budget_changes = Vec::new();
    let names: BTreeSet<&String> = before.budgets.keys().chain(after.budgets.keys()).collect();
    for name in names {
        let old = before.budgets.get(name).copied();
        let new = after.budgets.get(name).copied();
        if old != new {
            budget_changes.push((name.clone(), old, new));
        }
    }
    SelectionDiff {
        added,
        removed,
        budget_changes,
        scope_changed: before.stratum != after.stratum || before.profile != after.profile,
    }
}
