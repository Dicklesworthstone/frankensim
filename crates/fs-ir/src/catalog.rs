//! The machine-readable operator catalog (bead gp3.6, plan §11.1):
//! `catalog(query)` is how an agent DISCOVERS the system. Entries are
//! generated FROM the code surfaces that admission, lowering, and the
//! query menu actually execute — hand-written prose is forbidden for
//! facts the code knows, and the battery proves the remaining declared
//! facts against execution (a deliberately mismatched declaration
//! fails, so documentation rots loudly instead of silently).
//!
//! This module is also the REGISTRY the other surfaces consume:
//! [`ARITH_SAME_DIMS`] and [`COMPARE_FORMS`] are imported by admission's
//! dimensional dimension, and [`SUGAR_VERBS`] is what the drift battery
//! holds against `lower`'s dispatch — one source of truth, no parallel
//! lists.

use crate::query::Qoi;
use fs_qty::Dims;
use std::fmt::Write as _;

/// Catalog schema version: bound into the canonical JSON export so
/// agents can detect representation changes across releases.
/// v2: typed `consumes`/`produces` signature columns (the
/// signature-compatibility search vocabulary).
pub const CATALOG_SCHEMA_VERSION: u32 = 2;

/// Arithmetic forms whose operands must share exact dimensions
/// (admission's dimensional dimension consumes this list).
pub const ARITH_SAME_DIMS: &[&str] = &["+", "-", "min", "max"];

/// Comparison forms: same-dims operands, dimensionless verdict
/// (admission's dimensional dimension consumes this list).
pub const COMPARE_FORMS: &[&str] = &["=", "<", ">", "<=", ">="];

/// Progressive-disclosure sugar verbs `lower` expands. The drift
/// battery executes every one of these through `lower` and refuses if
/// the dispatch and this registry disagree in either direction.
pub const SUGAR_VERBS: &[&str] = &["optimize-shape", "simulate-pour"];

/// Determinism class of an operator's execution (typed here so it is a
/// fact the code knows; CONTRACT.md prose may elaborate, never
/// contradict).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeterminismClass {
    /// Bit-identical across runs, thread counts, and ISAs.
    BitStableCrossIsa,
    /// Deterministic given the same host/ISA; last-ULP movement across
    /// ISAs is possible and documented.
    DeterministicWithinIsa,
    /// A measurement of the physical machine (never replayable).
    Measured,
    /// Pure syntax/structure transformation — no floating point at all.
    Structural,
}

impl DeterminismClass {
    /// Stable token for the canonical export.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::BitStableCrossIsa => "bit-stable-cross-isa",
            Self::DeterministicWithinIsa => "deterministic-within-isa",
            Self::Measured => "measured",
            Self::Structural => "structural",
        }
    }
}

/// Cancellation behavior at this operator's boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CancellationBehavior {
    /// Bounded synchronous work; callers own chunking and polls.
    BoundedSynchronous,
    /// Long work polls an injected `Cx` at bounded checkpoints
    /// (request → drain → finalize).
    CxPolled,
}

impl CancellationBehavior {
    /// Stable token for the canonical export.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::BoundedSynchronous => "bounded-synchronous",
            Self::CxPolled => "cx-polled",
        }
    }
}

/// Ambition tag (plan Ambition-Tag rule).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ambition {
    /// Shipping surface.
    S,
    /// Flagship loop surface.
    F,
    /// Moonshot surface.
    M,
}

impl Ambition {
    /// Stable token for the canonical export.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::S => "S",
            Self::F => "F",
            Self::M => "M",
        }
    }
}

/// A typed value kind in an operator's signature: what the operator
/// CONSUMES and PRODUCES semantically, beyond the surface syntax of
/// its arguments. This is the signature-compatibility search
/// vocabulary ("what consumes a field and produces a probability?"),
/// so tokens are stable and closed — extending the vocabulary is a
/// schema-version event, not an edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SigType {
    /// A named simulation field (temperature, stress, …).
    Field,
    /// A named spatial region.
    Region,
    /// A dimensioned quantity `(qty v unit)`.
    Quantity,
    /// A named uncertainty environment (probabilistic QoIs).
    Environment,
    /// Design levers (`:over` targets).
    Levers,
    /// A scalar objective form.
    ObjectiveScalar,
    /// A vessel/geometry reference.
    Vessel,
    /// A fluid/material reference.
    Fluid,
    /// A pour/motion schedule.
    Schedule,
    /// An explicit lowered IR form (what sugar expansion yields).
    IrForm,
    /// A dimensioned QoI value under tolerance semantics.
    QoiValue,
    /// A dimensionless probability.
    Probability,
    /// A dimensionless comparison verdict.
    Verdict,
    /// A registered long-running campaign handle (ledger-bound).
    Campaign,
}

impl SigType {
    /// Stable token for the canonical export and query surface.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Field => "field",
            Self::Region => "region",
            Self::Quantity => "quantity",
            Self::Environment => "environment",
            Self::Levers => "levers",
            Self::ObjectiveScalar => "objective-scalar",
            Self::Vessel => "vessel",
            Self::Fluid => "fluid",
            Self::Schedule => "schedule",
            Self::IrForm => "ir-form",
            Self::QoiValue => "qoi-value",
            Self::Probability => "probability",
            Self::Verdict => "verdict",
            Self::Campaign => "campaign",
        }
    }

    /// Parse a query token. `None` for tokens outside the vocabulary —
    /// query filters treat those as matching no entry (fail-safe), so a
    /// typo'd search returns the empty set rather than everything.
    #[must_use]
    pub fn from_name(token: &str) -> Option<Self> {
        match token {
            "field" => Some(Self::Field),
            "region" => Some(Self::Region),
            "quantity" => Some(Self::Quantity),
            "environment" => Some(Self::Environment),
            "levers" => Some(Self::Levers),
            "objective-scalar" => Some(Self::ObjectiveScalar),
            "vessel" => Some(Self::Vessel),
            "fluid" => Some(Self::Fluid),
            "schedule" => Some(Self::Schedule),
            "ir-form" => Some(Self::IrForm),
            "qoi-value" => Some(Self::QoiValue),
            "probability" => Some(Self::Probability),
            "verdict" => Some(Self::Verdict),
            "campaign" => Some(Self::Campaign),
            _ => None,
        }
    }
}

/// What kind of operator an entry describes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OperatorKind {
    /// A progressive-disclosure sugar verb `lower` expands; carries the
    /// exact target head and the DEFAULTS the expansion injects when
    /// the caller omits them (drift-checked against the real trace).
    SugarVerb {
        /// The explicit operator this verb lowers to.
        lowers_to: &'static str,
        /// Injected-default trace lines exactly as `lower` records them
        /// for a minimal invocation.
        injected_defaults: &'static [&'static str],
    },
    /// A namespaced core operator (`namespace.name`) executed by a
    /// downstream engine.
    CoreOperator,
    /// A declarative QoI operator from the fixed query menu.
    QoiOperator {
        /// Linear in the field (planner metadata).
        linear: bool,
        /// Adjoint availability (planner metadata).
        adjoint_available: bool,
        /// Accuracy-ladder applicability (planner metadata).
        ladder_applicable: bool,
        /// Probabilistic QoIs require a named environment.
        probabilistic: bool,
    },
    /// An arithmetic/comparison form of the study surface.
    ArithmeticForm {
        /// Dimensional rule token ("same-dims" operands; comparisons
        /// yield a dimensionless verdict).
        dims_rule: &'static str,
    },
}

/// One argument in an operator's surface signature.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArgSpec {
    /// Argument name (keyword name for keyworded args).
    pub name: &'static str,
    /// Whether the argument is required.
    pub required: bool,
    /// Surface syntax description (an IR form, not prose).
    pub syntax: &'static str,
}

/// One catalog entry: every column the plan's catalog contract names.
#[derive(Debug, Clone, PartialEq)]
pub struct OperatorEntry {
    /// Operator head as it appears in IR programs.
    pub name: &'static str,
    /// What kind of operator this is (with kind-specific metadata).
    pub kind: OperatorKind,
    /// Surface signature: arguments in order.
    pub args: &'static [ArgSpec],
    /// Typed signature: semantic kinds this operator consumes, in
    /// argument order (the signature-compatibility search columns).
    pub consumes: &'static [SigType],
    /// Typed signature: the semantic kind this operator produces.
    pub produces: SigType,
    /// What the operator produces (surface description token).
    pub output: &'static str,
    /// Dimensional rule token for the dims dimension.
    pub dims_rule: &'static str,
    /// The key under which `AdmissionContext.cost_models` must register
    /// a sealed cost model to price this operator (None: never priced).
    pub cost_model_key: Option<&'static str>,
    /// CONTRACT.md error-model reference (crate whose `## Error model`
    /// section governs refusals at this boundary).
    pub error_model: &'static str,
    /// Determinism class (typed; see [`DeterminismClass`]).
    pub determinism: DeterminismClass,
    /// Capability grants required to admit this operator (exact names
    /// or `ns.*` globs, admission's capability dimension).
    pub capabilities: &'static [&'static str],
    /// Cancellation behavior at this boundary.
    pub cancellation: CancellationBehavior,
    /// Executable worked examples (IR s-expressions). The battery
    /// parses and executes every example; broken documentation fails
    /// the build.
    pub examples: &'static [&'static str],
    /// Semantic version of this operator's surface.
    pub semver: &'static str,
    /// Evidence model-card reference (evidence-model-ledger id), if
    /// one exists yet.
    pub model_card: Option<&'static str>,
    /// Ambition tag.
    pub ambition: Ambition,
    /// Feature flag gating the operator (None: always available).
    pub feature_flag: Option<&'static str>,
}

/// A structured catalog query. All present filters must match
/// (conjunction); absent filters match everything.
#[derive(Debug, Clone, Default)]
pub struct CatalogQuery {
    /// Exact name or `ns.*` glob.
    pub name: Option<String>,
    /// Namespace prefix ("ascent", "flux", "query", "form", "sugar").
    pub namespace: Option<String>,
    /// Only operators admissible under this capability grant set.
    pub granted_capabilities: Option<Vec<String>>,
    /// Signature compatibility: every listed [`SigType`] token must be
    /// among the entry's consumed kinds ("what can I feed an X?").
    /// Tokens outside the vocabulary match no entry.
    pub consumes: Option<Vec<String>>,
    /// Signature compatibility: the entry's produced kind must equal
    /// this [`SigType`] token ("what gives me a Y?"). Tokens outside
    /// the vocabulary match no entry.
    pub produces: Option<String>,
    /// Case-insensitive token that must appear in name/output/syntax.
    pub text: Option<String>,
}

/// The generated catalog.
#[derive(Debug)]
pub struct Catalog {
    entries: Vec<OperatorEntry>,
}

/// A validation refusal: the catalog and the code disagree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CatalogDrift {
    /// The entry (or registry surface) at fault.
    pub entry: String,
    /// What disagreed.
    pub detail: String,
}

impl Catalog {
    /// Generate the built-in catalog from the live code surfaces.
    ///
    /// QoI entries are DERIVED from [`Qoi`]'s own metadata (not
    /// re-authored); arithmetic forms are derived from the shared
    /// registries admission consumes.
    #[must_use]
    pub fn builtin() -> Self {
        let mut entries = vec![
            OperatorEntry {
                name: "optimize-shape",
                kind: OperatorKind::SugarVerb {
                    lowers_to: "ascent.optimize",
                    injected_defaults: &[
                        ":method (lbfgs :m 17)",
                        ":until (grad-norm 1e-5)",
                        ":emit (ledger report)",
                    ],
                },
                args: &[
                    ArgSpec {
                        name: "min",
                        required: true,
                        syntax: ":min <objective>",
                    },
                    ArgSpec {
                        name: "over",
                        required: true,
                        syntax: ":over <levers>",
                    },
                    ArgSpec {
                        name: "method",
                        required: false,
                        syntax: ":method <method-form>",
                    },
                    ArgSpec {
                        name: "until",
                        required: false,
                        syntax: ":until <stop-form>",
                    },
                ],
                consumes: &[SigType::ObjectiveScalar, SigType::Levers],
                produces: SigType::IrForm,
                output: "explicit (ascent.optimize …) form with every default named in the trace",
                dims_rule: "objective-scalar",
                cost_model_key: Some("ascent.optimize"),
                error_model: "fs-ir",
                determinism: DeterminismClass::Structural,
                capabilities: &["ascent.optimize"],
                cancellation: CancellationBehavior::BoundedSynchronous,
                examples: &["(optimize-shape :min (mass bracket) :over (thickness))"],
                semver: "1.0.0",
                model_card: None,
                ambition: Ambition::S,
                feature_flag: None,
            },
            OperatorEntry {
                name: "simulate-pour",
                kind: OperatorKind::SugarVerb {
                    lowers_to: "flux.free-surface-lbm",
                    injected_defaults: &[],
                },
                args: &[
                    ArgSpec {
                        name: "vessel",
                        required: true,
                        syntax: "<vessel>",
                    },
                    ArgSpec {
                        name: "fluid",
                        required: true,
                        syntax: "<fluid>",
                    },
                    ArgSpec {
                        name: "schedule",
                        required: true,
                        syntax: "<schedule>",
                    },
                ],
                consumes: &[SigType::Vessel, SigType::Fluid, SigType::Schedule],
                produces: SigType::IrForm,
                output: "explicit (flux.free-surface-lbm …) form",
                dims_rule: "engine-defined",
                cost_model_key: Some("flux.free-surface-lbm"),
                error_model: "fs-ir",
                determinism: DeterminismClass::Structural,
                capabilities: &["flux.free-surface-lbm"],
                cancellation: CancellationBehavior::BoundedSynchronous,
                examples: &["(simulate-pour crucible bronze schedule-a)"],
                semver: "1.0.0",
                model_card: None,
                ambition: Ambition::S,
                feature_flag: None,
            },
            OperatorEntry {
                name: "ascent.optimize",
                kind: OperatorKind::CoreOperator,
                args: &[
                    ArgSpec {
                        name: "objective",
                        required: true,
                        syntax: "(min <objective>)",
                    },
                    ArgSpec {
                        name: "over",
                        required: true,
                        syntax: ":over <levers>",
                    },
                    ArgSpec {
                        name: "method",
                        required: true,
                        syntax: ":method <method-form>",
                    },
                    ArgSpec {
                        name: "until",
                        required: true,
                        syntax: ":until <stop-form>",
                    },
                    ArgSpec {
                        name: "emit",
                        required: true,
                        syntax: ":emit (ledger report)",
                    },
                ],
                consumes: &[SigType::ObjectiveScalar, SigType::Levers],
                produces: SigType::Campaign,
                output: "optimization campaign registered against the design ledger",
                dims_rule: "objective-scalar",
                cost_model_key: Some("ascent.optimize"),
                error_model: "fs-ir",
                determinism: DeterminismClass::DeterministicWithinIsa,
                capabilities: &["ascent.optimize"],
                cancellation: CancellationBehavior::CxPolled,
                examples: &[
                    "(ascent.optimize (min (mass bracket)) :over (thickness) :method (lbfgs :m 17) :until (grad-norm 1e-5) :emit (ledger report))",
                ],
                semver: "1.0.0",
                model_card: None,
                ambition: Ambition::S,
                feature_flag: None,
            },
            OperatorEntry {
                name: "flux.free-surface-lbm",
                kind: OperatorKind::CoreOperator,
                args: &[
                    ArgSpec {
                        name: "vessel",
                        required: true,
                        syntax: "<vessel>",
                    },
                    ArgSpec {
                        name: "fluid",
                        required: true,
                        syntax: "<fluid>",
                    },
                    ArgSpec {
                        name: "schedule",
                        required: true,
                        syntax: "<schedule>",
                    },
                ],
                consumes: &[SigType::Vessel, SigType::Fluid, SigType::Schedule],
                produces: SigType::Campaign,
                output: "free-surface LBM campaign",
                dims_rule: "engine-defined",
                cost_model_key: Some("flux.free-surface-lbm"),
                error_model: "fs-ir",
                determinism: DeterminismClass::BitStableCrossIsa,
                capabilities: &["flux.free-surface-lbm"],
                cancellation: CancellationBehavior::CxPolled,
                examples: &["(flux.free-surface-lbm crucible bronze schedule-a)"],
                semver: "1.0.0",
                model_card: None,
                ambition: Ambition::F,
                feature_flag: None,
            },
        ];
        entries.extend(qoi_entries());
        entries.extend(arithmetic_entries());
        entries.sort_by(|a, b| a.name.cmp(b.name));
        Self { entries }
    }

    /// Every entry, in canonical (name) order.
    #[must_use]
    pub fn entries(&self) -> &[OperatorEntry] {
        &self.entries
    }

    /// Battery-only constructor for drift-refusal fixtures.
    #[doc(hidden)]
    #[must_use]
    pub fn from_entries_for_test(mut entries: Vec<OperatorEntry>) -> Self {
        entries.sort_by(|a, b| a.name.cmp(b.name));
        Self { entries }
    }

    /// Structured search: all present filters are conjunctive.
    #[must_use]
    pub fn query(&self, query: &CatalogQuery) -> Vec<&OperatorEntry> {
        self.entries
            .iter()
            .filter(|entry| {
                if let Some(name) = &query.name
                    && !name_matches(name, entry.name)
                {
                    return false;
                }
                if let Some(namespace) = &query.namespace
                    && entry_namespace(entry) != namespace.as_str()
                {
                    return false;
                }
                if let Some(grants) = &query.granted_capabilities
                    && !entry
                        .capabilities
                        .iter()
                        .all(|needed| grants.iter().any(|grant| grant_covers(grant, needed)))
                {
                    return false;
                }
                if let Some(consumes) = &query.consumes
                    && !consumes.iter().all(|token| {
                        SigType::from_name(token).is_some_and(|t| entry.consumes.contains(&t))
                    })
                {
                    return false;
                }
                if let Some(produces) = &query.produces
                    && SigType::from_name(produces) != Some(entry.produces)
                {
                    return false;
                }
                if let Some(text) = &query.text {
                    let needle = text.to_ascii_lowercase();
                    let haystack = format!(
                        "{} {} {}",
                        entry.name,
                        entry.output,
                        entry
                            .args
                            .iter()
                            .map(|a| a.syntax)
                            .collect::<Vec<_>>()
                            .join(" ")
                    )
                    .to_ascii_lowercase();
                    if !haystack.contains(&needle) {
                        return false;
                    }
                }
                true
            })
            .collect()
    }

    /// Canonical, diffable JSON-lines export: one deterministic line
    /// per entry plus a version header. Byte-identical across
    /// generations of the same code (proven in the battery).
    #[must_use]
    pub fn to_canonical_jsonl(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "{{\"schema\":\"fs-ir-operator-catalog\",\"version\":{CATALOG_SCHEMA_VERSION},\"entries\":{}}}",
            self.entries.len()
        );
        for entry in &self.entries {
            let kind = match &entry.kind {
                OperatorKind::SugarVerb {
                    lowers_to,
                    injected_defaults,
                } => format!(
                    "{{\"sugar\":{{\"lowers_to\":\"{}\",\"injected_defaults\":[{}]}}}}",
                    lowers_to,
                    injected_defaults
                        .iter()
                        .map(|d| format!("\"{}\"", escape(d)))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                OperatorKind::CoreOperator => "{\"core\":{}}".to_string(),
                OperatorKind::QoiOperator {
                    linear,
                    adjoint_available,
                    ladder_applicable,
                    probabilistic,
                } => format!(
                    "{{\"qoi\":{{\"linear\":{linear},\"adjoint_available\":{adjoint_available},\"ladder_applicable\":{ladder_applicable},\"probabilistic\":{probabilistic}}}}}"
                ),
                OperatorKind::ArithmeticForm { dims_rule } => {
                    format!("{{\"form\":{{\"dims_rule\":\"{dims_rule}\"}}}}")
                }
            };
            let args = entry
                .args
                .iter()
                .map(|a| {
                    format!(
                        "{{\"name\":\"{}\",\"required\":{},\"syntax\":\"{}\"}}",
                        a.name,
                        a.required,
                        escape(a.syntax)
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            let consumes = entry
                .consumes
                .iter()
                .map(|t| format!("\"{}\"", t.name()))
                .collect::<Vec<_>>()
                .join(",");
            let _ = writeln!(
                out,
                "{{\"name\":\"{}\",\"kind\":{kind},\"args\":[{args}],\"consumes\":[{consumes}],\"produces\":\"{}\",\"output\":\"{}\",\"dims_rule\":\"{}\",\"cost_model_key\":{},\"error_model\":\"{}\",\"determinism\":\"{}\",\"capabilities\":[{}],\"cancellation\":\"{}\",\"examples\":[{}],\"semver\":\"{}\",\"model_card\":{},\"ambition\":\"{}\",\"feature_flag\":{}}}",
                entry.name,
                entry.produces.name(),
                escape(entry.output),
                entry.dims_rule,
                entry
                    .cost_model_key
                    .map_or("null".to_string(), |k| format!("\"{k}\"")),
                entry.error_model,
                entry.determinism.name(),
                entry
                    .capabilities
                    .iter()
                    .map(|c| format!("\"{c}\""))
                    .collect::<Vec<_>>()
                    .join(","),
                entry.cancellation.name(),
                entry
                    .examples
                    .iter()
                    .map(|e| format!("\"{}\"", escape(e)))
                    .collect::<Vec<_>>()
                    .join(","),
                entry.semver,
                entry
                    .model_card
                    .map_or("null".to_string(), |m| format!("\"{m}\"")),
                entry.ambition.name(),
                entry
                    .feature_flag
                    .map_or("null".to_string(), |f| format!("\"{f}\"")),
            );
        }
        out
    }

    /// Diff this catalog against a previously exported canonical JSONL:
    /// which operator names were added, removed, or changed. Agents ask
    /// "what changed since the version I know".
    #[must_use]
    pub fn diff(&self, old_jsonl: &str) -> CatalogDiffReport {
        let old: std::collections::BTreeMap<String, String> = old_jsonl
            .lines()
            .skip(1)
            .filter_map(|line| {
                line.split_once("\"name\":\"")
                    .and_then(|(_, rest)| rest.split_once('"'))
                    .map(|(name, _)| (name.to_string(), line.to_string()))
            })
            .collect();
        let current = self.to_canonical_jsonl();
        let new: std::collections::BTreeMap<String, String> = current
            .lines()
            .skip(1)
            .filter_map(|line| {
                line.split_once("\"name\":\"")
                    .and_then(|(_, rest)| rest.split_once('"'))
                    .map(|(name, _)| (name.to_string(), line.to_string()))
            })
            .collect();
        CatalogDiffReport {
            added: new
                .keys()
                .filter(|k| !old.contains_key(*k))
                .cloned()
                .collect(),
            removed: old
                .keys()
                .filter(|k| !new.contains_key(*k))
                .cloned()
                .collect(),
            changed: new
                .iter()
                .filter(|(k, line)| old.get(*k).is_some_and(|o| o != *line))
                .map(|(k, _)| k.clone())
                .collect(),
        }
    }

    /// Structural validation: unique names, sugar targets present,
    /// examples non-empty, registry agreement.
    ///
    /// # Errors
    /// One [`CatalogDrift`] per disagreement.
    pub fn validate(&self) -> Result<(), Vec<CatalogDrift>> {
        let mut drifts = Vec::new();
        let mut names = std::collections::BTreeSet::new();
        for entry in &self.entries {
            if !names.insert(entry.name) {
                drifts.push(CatalogDrift {
                    entry: entry.name.to_string(),
                    detail: "duplicate catalog name".to_string(),
                });
            }
            if entry.examples.is_empty() {
                drifts.push(CatalogDrift {
                    entry: entry.name.to_string(),
                    detail: "no executable example".to_string(),
                });
            }
            if entry.consumes.is_empty() {
                drifts.push(CatalogDrift {
                    entry: entry.name.to_string(),
                    detail: "no consumed signature types".to_string(),
                });
            }
            if let OperatorKind::QoiOperator { probabilistic, .. } = &entry.kind {
                let produces_probability = entry.produces == SigType::Probability;
                if *probabilistic != produces_probability {
                    drifts.push(CatalogDrift {
                        entry: entry.name.to_string(),
                        detail: format!(
                            "probabilistic={probabilistic} but produces {:?} — a QoI produces \
                             probability exactly when it is probabilistic",
                            entry.produces.name()
                        ),
                    });
                }
            }
            if let OperatorKind::SugarVerb { lowers_to, .. } = &entry.kind {
                if entry.produces != SigType::IrForm {
                    drifts.push(CatalogDrift {
                        entry: entry.name.to_string(),
                        detail: format!(
                            "sugar verbs produce ir-form (structural expansion), not {:?}",
                            entry.produces.name()
                        ),
                    });
                }
                if !self.entries.iter().any(|e| e.name == *lowers_to) {
                    drifts.push(CatalogDrift {
                        entry: entry.name.to_string(),
                        detail: format!("lowers_to target {lowers_to:?} has no catalog entry"),
                    });
                }
                if !SUGAR_VERBS.contains(&entry.name) {
                    drifts.push(CatalogDrift {
                        entry: entry.name.to_string(),
                        detail: "sugar entry missing from the SUGAR_VERBS registry".to_string(),
                    });
                }
            }
        }
        for verb in SUGAR_VERBS {
            if !self.entries.iter().any(|e| e.name == *verb) {
                drifts.push(CatalogDrift {
                    entry: (*verb).to_string(),
                    detail: "registered sugar verb has no catalog entry".to_string(),
                });
            }
        }
        // Growth discipline: lower's dispatch table and the registry
        // must be SET-EQUAL — a new dispatch arm cannot land without a
        // registry entry (and hence a cataloged, example-executed
        // operator), nor a registry entry without a dispatch arm.
        let dispatch: std::collections::BTreeSet<&str> =
            crate::lower::dispatch_heads().into_iter().collect();
        let registry: std::collections::BTreeSet<&str> = SUGAR_VERBS.iter().copied().collect();
        if dispatch != registry {
            drifts.push(CatalogDrift {
                entry: "sugar-dispatch".to_string(),
                detail: format!(
                    "lower dispatch heads {dispatch:?} and the SUGAR_VERBS registry \
                     {registry:?} are not set-equal"
                ),
            });
        }
        for form in ARITH_SAME_DIMS.iter().chain(COMPARE_FORMS) {
            if !self.entries.iter().any(|e| e.name == *form) {
                drifts.push(CatalogDrift {
                    entry: (*form).to_string(),
                    detail: "registered arithmetic form has no catalog entry".to_string(),
                });
            }
        }
        if drifts.is_empty() {
            Ok(())
        } else {
            Err(drifts)
        }
    }
}

/// What changed between two catalog exports.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CatalogDiffReport {
    /// Operator names present now but not before.
    pub added: Vec<String>,
    /// Operator names present before but not now.
    pub removed: Vec<String>,
    /// Operator names whose entry line changed.
    pub changed: Vec<String>,
}

fn entry_namespace(entry: &OperatorEntry) -> &'static str {
    match &entry.kind {
        OperatorKind::SugarVerb { .. } => "sugar",
        OperatorKind::QoiOperator { .. } => "query",
        OperatorKind::ArithmeticForm { .. } => "form",
        OperatorKind::CoreOperator => {
            entry.name.split_once('.').map_or("core", |(namespace, _)| {
                // Namespaces are 'static by construction: core names are
                // 'static string literals.
                match namespace {
                    "ascent" => "ascent",
                    "flux" => "flux",
                    _ => "core",
                }
            })
        }
    }
}

fn name_matches(pattern: &str, name: &str) -> bool {
    if let Some(prefix) = pattern.strip_suffix(".*") {
        name.strip_prefix(prefix)
            .is_some_and(|rest| rest.starts_with('.'))
    } else {
        pattern == name
    }
}

fn grant_covers(grant: &str, needed: &str) -> bool {
    if let Some(prefix) = grant.strip_suffix(".*") {
        needed
            .strip_prefix(prefix)
            .is_some_and(|rest| rest.starts_with('.'))
    } else {
        grant == needed
    }
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\").replace('"', "\\\"")
}

/// QoI entries derived from the live [`Qoi`] menu: kind names, planner
/// metadata, probabilistic flag, and dims behavior all come from the
/// enum's own methods executed here — not re-authored.
fn qoi_entries() -> Vec<OperatorEntry> {
    let samples: [(
        Qoi,
        &'static [ArgSpec],
        &'static [SigType],
        SigType,
        &'static str,
        &'static [&'static str],
    ); 3] = [
        (
            Qoi::MaxOverRegion {
                field: "catalog-probe-field".to_string(),
                region: "everywhere".to_string(),
            },
            &[
                ArgSpec {
                    name: "field",
                    required: true,
                    syntax: "<field>",
                },
                ArgSpec {
                    name: "region",
                    required: true,
                    syntax: "<region>",
                },
            ],
            &[SigType::Field, SigType::Region],
            SigType::QoiValue,
            "field dimensions",
            &[
                "(query (max-over-region temperature hot-corner) (tolerance 0.5 kelvin) 100.0 3600.0)",
            ],
        ),
        (
            Qoi::Integral {
                field: "catalog-probe-field".to_string(),
                region: "everywhere".to_string(),
            },
            &[
                ArgSpec {
                    name: "field",
                    required: true,
                    syntax: "<field>",
                },
                ArgSpec {
                    name: "region",
                    required: true,
                    syntax: "<region>",
                },
            ],
            &[SigType::Field, SigType::Region],
            SigType::QoiValue,
            "field dimensions x volume",
            &["(query (integral heat-flux outlet) (tolerance 1.0 watt) 250.0 7200.0)"],
        ),
        (
            Qoi::Exceedance {
                field: "catalog-probe-field".to_string(),
                region: "everywhere".to_string(),
                threshold: 1.0,
                threshold_dims: Dims([0, 0, 0, 1, 0, 0]),
                environment: "catalog-probe-env".to_string(),
            },
            &[
                ArgSpec {
                    name: "field",
                    required: true,
                    syntax: "<field>",
                },
                ArgSpec {
                    name: "region",
                    required: true,
                    syntax: "<region>",
                },
                ArgSpec {
                    name: "threshold",
                    required: true,
                    syntax: "<threshold qty>",
                },
                ArgSpec {
                    name: "environment",
                    required: true,
                    syntax: "<named environment>",
                },
            ],
            &[
                SigType::Field,
                SigType::Region,
                SigType::Quantity,
                SigType::Environment,
            ],
            SigType::Probability,
            "dimensionless probability",
            &[
                "(query (exceedance stress weld-line (qty 250.0 megapascal) storm-env) (confidence 0.95) 500.0 86400.0)",
            ],
        ),
    ];
    samples
        .into_iter()
        .map(|(qoi, args, consumes, produces, output, examples)| {
            let meta = qoi.meta();
            OperatorEntry {
                name: qoi_static_name(&qoi),
                kind: OperatorKind::QoiOperator {
                    linear: meta.linear,
                    adjoint_available: meta.adjoint_available,
                    ladder_applicable: meta.ladder_applicable,
                    probabilistic: qoi.is_probabilistic(),
                },
                args,
                consumes,
                produces,
                output,
                dims_rule: "qoi-value-dims",
                cost_model_key: None,
                error_model: "fs-ir",
                determinism: DeterminismClass::Structural,
                capabilities: &[],
                cancellation: CancellationBehavior::BoundedSynchronous,
                examples,
                semver: "1.0.0",
                model_card: None,
                ambition: Ambition::S,
                feature_flag: None,
            }
        })
        .collect()
}

fn qoi_static_name(qoi: &Qoi) -> &'static str {
    match qoi {
        Qoi::MaxOverRegion { .. } => "max-over-region",
        Qoi::Integral { .. } => "integral",
        Qoi::Exceedance { .. } => "exceedance",
    }
}

/// Arithmetic/comparison entries derived from the shared registries.
fn arithmetic_entries() -> Vec<OperatorEntry> {
    let arith = ARITH_SAME_DIMS.iter().map(|name| OperatorEntry {
        name,
        kind: OperatorKind::ArithmeticForm {
            dims_rule: "same-dims",
        },
        args: &[
            ArgSpec {
                name: "lhs",
                required: true,
                syntax: "<operand>",
            },
            ArgSpec {
                name: "rhs",
                required: true,
                syntax: "<operand>",
            },
        ],
        consumes: &[SigType::Quantity, SigType::Quantity],
        produces: SigType::Quantity,
        output: "operand dimensions",
        dims_rule: "same-dims",
        cost_model_key: None,
        error_model: "fs-ir",
        determinism: DeterminismClass::Structural,
        capabilities: &[],
        cancellation: CancellationBehavior::BoundedSynchronous,
        examples: &["(+ (qty 1.0 meter) (qty 2.0 meter))"],
        semver: "1.0.0",
        model_card: None,
        ambition: Ambition::S,
        feature_flag: None,
    });
    let compare = COMPARE_FORMS.iter().map(|name| OperatorEntry {
        name,
        kind: OperatorKind::ArithmeticForm {
            dims_rule: "same-dims-comparison",
        },
        args: &[
            ArgSpec {
                name: "lhs",
                required: true,
                syntax: "<operand>",
            },
            ArgSpec {
                name: "rhs",
                required: true,
                syntax: "<operand>",
            },
        ],
        consumes: &[SigType::Quantity, SigType::Quantity],
        produces: SigType::Verdict,
        output: "dimensionless verdict",
        dims_rule: "same-dims-comparison",
        cost_model_key: None,
        error_model: "fs-ir",
        determinism: DeterminismClass::Structural,
        capabilities: &[],
        cancellation: CancellationBehavior::BoundedSynchronous,
        examples: &["(< (qty 1.0 meter) (qty 2.0 meter))"],
        semver: "1.0.0",
        model_card: None,
        ambition: Ambition::S,
        feature_flag: None,
    });
    arith.chain(compare).collect()
}
