//! The `.fsim` semantic model for the ratified thermal-design-assurance
//! vertical (record `frankensim-vertical-ratification-v1`), with the Five
//! Explicits mandatory and every remaining mandatory section named by a typed
//! violation when absent.
//!
//! Sections are `Option`s so recognition can be lenient and validation can
//! name every omission instead of refusing at the first; an ADMISSIBLE
//! project is exactly one whose [`ProjectSpec::validate`] returns no
//! violations. Declared-empty collections are deliberately distinct from
//! omitted sections: `(cooling (fans) (vents))` states "no fans" as a fact,
//! while a missing `cooling` section is a violation.

use std::collections::{BTreeMap, BTreeSet};

use fs_io::MeshSelector;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{EntityDeclaration, EntityId, EntityKind, InterfacePair, Violation};

/// SI dimension vectors this schema checks against.
pub mod dims {
    use fs_qty::Dims;

    /// Kelvin.
    pub const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
    /// Watt.
    pub const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
    /// Pascal.
    pub const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
    /// Second.
    pub const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
    /// Cubic metre per second.
    pub const VOLUMETRIC_FLOW: Dims = Dims([3, 0, -1, 0, 0, 0]);
    /// Square metre.
    pub const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
}

/// Project metadata: who the study is for and what decision it feeds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Metadata {
    /// Project name.
    pub name: String,
    /// Creation date (ISO-8601 day).
    pub created: String,
    /// Context of use: where the answer will be relied on.
    pub context_of_use: String,
    /// The engineering decision the result is intended to inform.
    pub intended_decision: String,
    /// Decision gate controlling whether an indeterminate result may be used.
    pub decision_gate: DecisionGate,
    /// Consequence framing for a wrong or unsupported decision.
    pub consequence: ConsequenceClass,
}

impl Metadata {
    /// Whether this context may proceed with an explicitly indeterminate
    /// assessment. Safety-significant contexts always require a determinate
    /// result even when the intended use is nominally exploratory.
    #[must_use]
    pub const fn permits_indeterminate(&self) -> bool {
        matches!(self.decision_gate, DecisionGate::ScopingEstimate)
            && !matches!(self.consequence, ConsequenceClass::SafetyCritical)
    }
}

/// Intended decision gate for one project context.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DecisionGate {
    /// Early scoping may retain an explicit Indeterminate result.
    ScopingEstimate,
    /// A design choice requires a determinate comparison.
    DesignSelection,
    /// Compliance sign-off requires a determinate, fully sourced verdict.
    ComplianceSignoff,
}

/// Consequence framing for the intended decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConsequenceClass {
    /// Advisory exploration with no direct release authority.
    Advisory,
    /// Reliability or service-life decision.
    Reliability,
    /// Damage, personnel, or regulatory safety consequence.
    SafetyCritical,
}

/// The versions pillar of the Five Explicits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Versions {
    /// `.fsim` schema version (must equal [`crate::FSIM_VERSION`]).
    pub schema: u32,
    /// Constellation lock digest the project was authored against.
    pub constellation: String,
    /// Workspace revision the project was authored against.
    pub workspace: String,
}

/// The seeds pillar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seeds {
    /// Master seed; every derived stream is keyed from it by logical
    /// identity, never by scheduling.
    pub master: u64,
}

/// The budgets pillar.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Budgets {
    /// Wall-clock solve budget (seconds).
    pub solve_time: QtyAny,
    /// Peak memory budget in bytes.
    pub memory_bytes: u64,
    /// Relative accuracy target for the primary quantity of interest.
    pub accuracy_rel: f64,
}

/// The units pillar: an explicit acknowledgment that stored quantities are
/// coherent SI base values, plus the preferred display family.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnitsDoctrine {
    /// Storage convention; the only admitted value is `"si-base"`.
    pub storage: String,
    /// Display preference recorded for report rendering.
    pub display: String,
}

/// One imported geometry artifact, referenced through its quarantine
/// receipt — geometry lives in artifacts, never inline in the project.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometryArtifact {
    /// Role in the assembly ("enclosure", "board", "heatsink", ...).
    pub role: String,
    /// Declared source format from the import receipt.
    pub format: String,
    /// Content hash of the raw source bytes from `fs_io::ImportReceipt`.
    pub source_hash: u64,
    /// Parser version that produced the receipt.
    pub parser_version: String,
}

/// One mesh-index-free assignment declared in the project file.
///
/// `target` is a project-local declared entity name. The L6 assignment
/// adapter resolves it to the actual [`EntityId`] and gives fs-io only that
/// persistent token; a raw mesh ordinal is never an identity.
#[derive(Debug, Clone, PartialEq)]
pub struct GeometryAssignment {
    /// Geometry artifact role this selector operates on.
    pub artifact: String,
    /// Declared region or interface name to bind.
    pub target: String,
    /// Explicit coordinate unit carried into the fs-io receipt.
    pub length_unit: String,
    /// Deterministic named/geometric selector.
    pub selector: MeshSelector,
    /// Whether overlap with other overlap-enabled assignments is intentional.
    pub allow_overlap: bool,
}

/// One entity declaration with a persistent identity. Names are the
/// project-local namespace; identities are recomputed from the declarations
/// through `fs-scenario`'s entity machinery, and an optional expected token
/// pins the identity so silent drift is caught on load.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntityDecl {
    /// Root assembly.
    Assembly {
        /// Declared (identity-bearing) name.
        name: String,
        /// Display name (never part of the identity).
        display: String,
        /// Expected identity token (`kind:hex`), verified when present.
        expect_id: Option<String>,
    },
    /// Part under an assembly.
    Part {
        /// Parent assembly's declared name.
        parent: String,
        /// Declared name.
        name: String,
        /// Display name.
        display: String,
        /// Expected identity token.
        expect_id: Option<String>,
    },
    /// Region under a part.
    Region {
        /// Parent part's declared name.
        parent: String,
        /// Declared name.
        name: String,
        /// Display name.
        display: String,
        /// Expected identity token.
        expect_id: Option<String>,
    },
    /// Interface between two regions.
    Interface {
        /// Parent assembly's declared name.
        parent: String,
        /// Declared name.
        name: String,
        /// Display name.
        display: String,
        /// Region name on the `from` side.
        from: String,
        /// Region name on the `to` side.
        to: String,
        /// Expected identity token.
        expect_id: Option<String>,
    },
}

impl EntityDecl {
    /// The declared (identity-bearing) name.
    #[must_use]
    pub fn name(&self) -> &str {
        match self {
            EntityDecl::Assembly { name, .. }
            | EntityDecl::Part { name, .. }
            | EntityDecl::Region { name, .. }
            | EntityDecl::Interface { name, .. } => name,
        }
    }

    fn expect_id(&self) -> Option<&str> {
        match self {
            EntityDecl::Assembly { expect_id, .. }
            | EntityDecl::Part { expect_id, .. }
            | EntityDecl::Region { expect_id, .. }
            | EntityDecl::Interface { expect_id, .. } => expect_id.as_deref(),
        }
    }
}

/// One material binding: a matdb card reference with state, admitted
/// temperature range, and source channel.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialBinding {
    /// Region (or part) declared name the card is bound to.
    pub region: String,
    /// Hex digest of the matdb card's content hash.
    pub card: String,
    /// Explicit claim pin (hex digest of a claim content hash) selecting
    /// between coexisting conflicting claims on the card. Optional: a card
    /// without conflicts needs no pin, and a conflict without a pin is a
    /// validation refusal — the system never auto-picks between sources.
    pub claim: Option<String>,
    /// Material state identifier on the card.
    pub state: String,
    /// Lower bound of the admitted temperature range.
    pub temp_lo: QtyAny,
    /// Upper bound of the admitted temperature range.
    pub temp_hi: QtyAny,
    /// Source channel the card came from.
    pub source: String,
}

/// One interface-card binding: a TIM/contact system card from matdb bound to
/// a declared interface entity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InterfaceCardBinding {
    /// Interface declared name the card is bound to.
    pub interface: String,
    /// Hex digest of the matdb interface-system card's content hash.
    pub card: String,
    /// Explicit claim pin (hex digest of a claim content hash) selecting
    /// between coexisting conflicting claims on the card; same semantics as
    /// [`MaterialBinding::claim`].
    pub claim: Option<String>,
    /// Source channel the card came from.
    pub source: String,
}

/// One power dissipation row of the power map.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerDissipation {
    /// Region declared name that dissipates.
    pub region: String,
    /// Dissipated power (W).
    pub watts: QtyAny,
    /// Duty factor in `0.0..=1.0`. Omission in the lenient wire spelling is
    /// the schema's ONE receipted power default (`1.0`).
    pub duty: f64,
}

/// One fan declaration (correlation-rung operating point).
#[derive(Debug, Clone, PartialEq)]
pub struct Fan {
    /// Fan name.
    pub name: String,
    /// Volumetric flow at the operating point (m^3/s).
    pub flow: QtyAny,
    /// Static pressure at the operating point (Pa).
    pub static_pressure: QtyAny,
}

/// One vent declaration.
#[derive(Debug, Clone, PartialEq)]
pub struct Vent {
    /// Surface/region declared name the vent occupies.
    pub region: String,
    /// Open area (m^2).
    pub area: QtyAny,
}

/// The cooling section: declared-empty lists are facts, not omissions.
#[derive(Debug, Clone, PartialEq)]
pub struct Cooling {
    /// Fans (may be empty: "no fans" is a declaration).
    pub fans: Vec<Fan>,
    /// Vents (may be empty).
    pub vents: Vec<Vent>,
    /// Non-modeled leakage/background dissipation (W).
    pub leakage: QtyAny,
}

/// The operating envelope.
#[derive(Debug, Clone, PartialEq)]
pub struct Envelope {
    /// Lowest ambient temperature the answer must cover (K).
    pub ambient_lo: QtyAny,
    /// Highest ambient temperature the answer must cover (K).
    pub ambient_hi: QtyAny,
    /// Ambient pressure (Pa).
    pub pressure: QtyAny,
}

/// Kind of retained requirement authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RequirementSourceKind {
    /// Published standard or code clause.
    Standard,
    /// Manufacturer datasheet or product specification.
    Datasheet,
    /// Versioned internal engineering policy.
    InternalPolicy,
    /// Explicit user declaration with no external-document authority.
    UserDeclaration,
}

/// Versioned document and exact locator supporting a requirement or policy.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RequirementSource {
    /// Authority family.
    pub kind: RequirementSourceKind,
    /// Stable document, dataset, or policy identity.
    pub document: String,
    /// Exact edition, revision, or semantic version.
    pub version: String,
    /// Clause, table, section, or declaration locator.
    pub locator: String,
}

/// Direction of a scalar requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RequirementDirection {
    /// The effective quantity must not exceed the limit.
    AtMost,
    /// The effective quantity must not fall below the limit.
    AtLeast,
}

/// Engineering meaning of violating a requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RequirementSeverity {
    /// Preferred derating or reliability margin.
    ReliabilityDerating,
    /// Component or assembly damage limit.
    DamageLimit,
    /// Safety or regulatory limit.
    SafetyCritical,
}

/// Safety-factor authority already reflected in the effective limit.
///
/// No universal multiply/divide rule is inferred. The source owns how the
/// finite factor was applied, and [`ThermalLimit::limit`] stores the effective
/// value consumed by compliance evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct SafetyFactorPolicy {
    /// Applied factor, finite and at least one.
    pub factor: f64,
    /// Versioned policy defining the application rule.
    pub source: RequirementSource,
}

/// One sourced thermal requirement with explicit decision semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalLimit {
    /// Quantity-of-interest identity consumed by the decision layer.
    pub qoi: String,
    /// Component class the limit applies to ("junction", "case", ...).
    pub class: String,
    /// Region declared name the limit is evaluated on.
    pub region: String,
    /// Comparison direction.
    pub direction: RequirementDirection,
    /// Effective limit temperature (K), after the retained factor policy.
    pub limit: QtyAny,
    /// Required margin below the limit (K).
    pub margin: QtyAny,
    /// Versioned authority for the base requirement.
    pub source: RequirementSource,
    /// Already-applied safety-factor authority.
    pub safety_factor: SafetyFactorPolicy,
    /// Consequence of violating this requirement.
    pub severity: RequirementSeverity,
}

/// Requirement-source version change requiring human review.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RequirementSourceReview {
    /// Quantity whose authority changed.
    pub qoi: String,
    /// Region whose requirement changed.
    pub region: String,
    /// `requirement` or `safety-factor`.
    pub role: &'static str,
    /// Stable source document identity.
    pub document: String,
    /// Exact source locator.
    pub locator: String,
    /// Previously reviewed version.
    pub previous_version: String,
    /// Newly declared version.
    pub current_version: String,
}

/// Compare two admitted requirement sets and name every source-version change.
///
/// Matching is by QoI, region, authority role, source kind, document, and
/// locator. New/deleted requirements and entirely different authorities are
/// not mislabeled as version bumps; those belong to ordinary project diff
/// review.
#[must_use]
pub fn requirement_source_reviews(
    previous: &[ThermalLimit],
    current: &[ThermalLimit],
) -> Vec<RequirementSourceReview> {
    type SourceKey = (
        String,
        String,
        &'static str,
        RequirementSourceKind,
        String,
        String,
    );

    fn sources(limit: &ThermalLimit) -> [(&'static str, &RequirementSource); 2] {
        [
            ("requirement", &limit.source),
            ("safety-factor", &limit.safety_factor.source),
        ]
    }

    let mut previous_versions = BTreeMap::<SourceKey, String>::new();
    for limit in previous {
        for (role, source) in sources(limit) {
            previous_versions.insert(
                (
                    limit.qoi.clone(),
                    limit.region.clone(),
                    role,
                    source.kind,
                    source.document.clone(),
                    source.locator.clone(),
                ),
                source.version.clone(),
            );
        }
    }

    let mut reviews = Vec::new();
    for limit in current {
        for (role, source) in sources(limit) {
            let key = (
                limit.qoi.clone(),
                limit.region.clone(),
                role,
                source.kind,
                source.document.clone(),
                source.locator.clone(),
            );
            if let Some(previous_version) = previous_versions.get(&key)
                && previous_version != &source.version
            {
                reviews.push(RequirementSourceReview {
                    qoi: limit.qoi.clone(),
                    region: limit.region.clone(),
                    role,
                    document: source.document.clone(),
                    locator: source.locator.clone(),
                    previous_version: previous_version.clone(),
                    current_version: source.version.clone(),
                });
            }
        }
    }
    reviews.sort();
    reviews
}

/// Solver settings: fidelity selection and the stop rule.
#[derive(Debug, Clone, PartialEq)]
pub struct SolverSettings {
    /// `"auto"` or an explicit registered capability id.
    pub fidelity: String,
    /// Relative stop-rule tolerance.
    pub tolerance_rel: f64,
}

/// One output request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputRequest {
    /// Output name.
    pub name: String,
    /// Kind: `"scalar"`, `"field"`, or `"report"`.
    pub kind: String,
}

/// One receipted default: the lenient wire spelling applied a documented
/// default and says so in the validation output (no silent defaults).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DefaultReceipt {
    /// Dotted path of the defaulted field.
    pub field: String,
    /// The applied value, spelled canonically.
    pub value: String,
    /// Why this default is admissible at all.
    pub rationale: &'static str,
}

/// The `.fsim` project: every section of the cooling vertical's user-facing
/// contract, with omissions representable so validation can name them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ProjectSpec {
    /// Project metadata + context of use + intended decision.
    pub metadata: Option<Metadata>,
    /// Five Explicits: versions.
    pub versions: Option<Versions>,
    /// Five Explicits: seeds.
    pub seeds: Option<Seeds>,
    /// Five Explicits: budgets.
    pub budgets: Option<Budgets>,
    /// Five Explicits: capabilities (required, at least one).
    pub capabilities: Option<Vec<String>>,
    /// Five Explicits: units doctrine.
    pub units: Option<UnitsDoctrine>,
    /// Imported geometry artifact references.
    pub geometry: Option<Vec<GeometryArtifact>>,
    /// Mesh-index-free region/interface assignments.
    pub assignments: Option<Vec<GeometryAssignment>>,
    /// Assembly/part/region/interface declarations.
    pub assembly: Option<Vec<EntityDecl>>,
    /// Material bindings.
    pub materials: Option<Vec<MaterialBinding>>,
    /// Interface (TIM/contact) card bindings.
    pub interface_cards: Option<Vec<InterfaceCardBinding>>,
    /// Power map.
    pub power: Option<Vec<PowerDissipation>>,
    /// Fans, vents, leakage.
    pub cooling: Option<Cooling>,
    /// Operating envelope.
    pub envelope: Option<Envelope>,
    /// Requirements with margins.
    pub requirements: Option<Vec<ThermalLimit>>,
    /// Solver settings.
    pub solver: Option<SolverSettings>,
    /// Output requests.
    pub outputs: Option<Vec<OutputRequest>>,
}

#[allow(clippy::too_many_arguments)]
fn interface_declaration(
    ids: &BTreeMap<String, EntityId>,
    kinds: &BTreeMap<String, &'static str>,
    parent: &str,
    name: &str,
    display: &str,
    from: &str,
    to: &str,
    out: &mut Vec<Violation>,
) -> Option<(EntityDeclaration, &'static str)> {
    let parent_ok = matches!(kinds.get(parent), Some(&"assembly"));
    let from_ok = matches!(kinds.get(from), Some(&"region"));
    let to_ok = matches!(kinds.get(to), Some(&"region"));
    if parent_ok && from_ok && to_ok {
        let pair = InterfacePair::ordered(ids[from], ids[to]);
        Some((
            EntityDeclaration::interface(ids[parent], name, pair).with_display_name(display),
            "interface",
        ))
    } else {
        out.push(violation(
            "project-entity-parent-unknown",
            format!("interface `{name}` references `{parent}`/`{from}`/`{to}`, not all declared above it with the right kinds"),
            "declare the parent assembly and both regions before the interface",
        ));
        None
    }
}

/// A 64-char lowercase-or-uppercase hex content digest, the only admitted
/// spelling for card and claim references.
fn is_hex_digest(text: &str) -> bool {
    text.len() == 64 && text.bytes().all(|b| b.is_ascii_hexdigit())
}

fn violation(code: &'static str, what: impl Into<String>, fix: impl Into<String>) -> Violation {
    Violation {
        code,
        what: what.into(),
        fix: fix.into(),
    }
}

fn check_dims(
    out: &mut Vec<Violation>,
    code: &'static str,
    context: &str,
    got: QtyAny,
    expected: Dims,
) {
    if got.dims != expected {
        out.push(violation(
            code,
            format!(
                "{context} has dimensions {:?}, expected {:?}",
                got.dims, expected
            ),
            format!("spell {context} as a quantity with dimensions {expected:?}"),
        ));
    }
}

impl ProjectSpec {
    /// Validate the project, naming every omission and inconsistency as a
    /// typed violation. Empty output means the project is admissible.
    #[must_use]
    pub fn validate(&self) -> Vec<Violation> {
        let mut out = Vec::new();
        self.check_sections(&mut out);
        self.check_quantities(&mut out);
        let ids = self.resolve_entities(&mut out);
        self.check_references(&ids, &mut out);
        self.check_assignment_coverage(&ids, &mut out);
        out
    }

    fn check_sections(&self, out: &mut Vec<Violation>) {
        let mandatory: [(&'static str, bool, &str); 17] = [
            (
                "project-metadata-missing",
                self.metadata.is_some(),
                "metadata",
            ),
            (
                "project-versions-missing",
                self.versions.is_some(),
                "versions",
            ),
            ("project-seeds-missing", self.seeds.is_some(), "seeds"),
            ("project-budgets-missing", self.budgets.is_some(), "budgets"),
            (
                "project-capabilities-missing",
                self.capabilities.is_some(),
                "capabilities",
            ),
            ("project-units-missing", self.units.is_some(), "units"),
            (
                "project-geometry-missing",
                self.geometry.is_some(),
                "geometry",
            ),
            (
                "project-assignments-missing",
                self.assignments.is_some(),
                "assignments",
            ),
            (
                "project-assembly-missing",
                self.assembly.is_some(),
                "assembly",
            ),
            (
                "project-materials-missing",
                self.materials.is_some(),
                "materials",
            ),
            (
                "project-interface-cards-missing",
                self.interface_cards.is_some(),
                "interface-cards",
            ),
            ("project-power-missing", self.power.is_some(), "power"),
            ("project-cooling-missing", self.cooling.is_some(), "cooling"),
            (
                "project-envelope-missing",
                self.envelope.is_some(),
                "envelope",
            ),
            (
                "project-requirements-missing",
                self.requirements.is_some(),
                "requirements",
            ),
            ("project-solver-missing", self.solver.is_some(), "solver"),
            ("project-outputs-missing", self.outputs.is_some(), "outputs"),
        ];
        for (code, present, section) in mandatory {
            if !present {
                out.push(violation(
                    code,
                    format!("mandatory section `{section}` is absent"),
                    format!("declare `({section} ...)`; the Five Explicits and every cooling section are required, never defaulted"),
                ));
            }
        }

        self.check_nonempty_sections(out);
        self.check_field_content(out);
    }

    fn check_nonempty_sections(&self, out: &mut Vec<Violation>) {
        let nonempty: [(&'static str, Option<usize>, &str); 7] = [
            (
                "project-capabilities-empty",
                self.capabilities.as_ref().map(Vec::len),
                "capabilities",
            ),
            (
                "project-geometry-empty",
                self.geometry.as_ref().map(Vec::len),
                "geometry",
            ),
            (
                "project-assignments-empty",
                self.assignments.as_ref().map(Vec::len),
                "assignments",
            ),
            (
                "project-assembly-empty",
                self.assembly.as_ref().map(Vec::len),
                "assembly",
            ),
            (
                "project-power-empty",
                self.power.as_ref().map(Vec::len),
                "power",
            ),
            (
                "project-requirements-empty",
                self.requirements.as_ref().map(Vec::len),
                "requirements",
            ),
            (
                "project-outputs-empty",
                self.outputs.as_ref().map(Vec::len),
                "outputs",
            ),
        ];
        for (code, len, section) in nonempty {
            if len == Some(0) {
                out.push(violation(
                    code,
                    format!("section `{section}` is declared but empty"),
                    format!("a cooling project needs at least one `{section}` row; an intentionally empty list is only meaningful for cooling fans/vents"),
                ));
            }
        }
    }

    fn check_field_content(&self, out: &mut Vec<Violation>) {
        if let Some(metadata) = &self.metadata {
            for (field, value) in [
                ("metadata.name", &metadata.name),
                ("metadata.created", &metadata.created),
                ("metadata.context-of-use", &metadata.context_of_use),
                ("metadata.intended-decision", &metadata.intended_decision),
            ] {
                if value.trim().is_empty() {
                    out.push(violation(
                        "project-metadata-field-empty",
                        format!("`{field}` is empty"),
                        format!(
                            "state `{field}`: the record of why this study exists is mandatory"
                        ),
                    ));
                }
            }
        }
        if let Some(requirements) = &self.requirements {
            for limit in requirements {
                for (field, value) in [
                    ("qoi", limit.qoi.as_str()),
                    ("class", limit.class.as_str()),
                    ("region", limit.region.as_str()),
                    ("source.document", limit.source.document.as_str()),
                    ("source.version", limit.source.version.as_str()),
                    ("source.locator", limit.source.locator.as_str()),
                    (
                        "safety-factor.source.document",
                        limit.safety_factor.source.document.as_str(),
                    ),
                    (
                        "safety-factor.source.version",
                        limit.safety_factor.source.version.as_str(),
                    ),
                    (
                        "safety-factor.source.locator",
                        limit.safety_factor.source.locator.as_str(),
                    ),
                ] {
                    if value.is_empty()
                        || value.trim() != value
                        || value.chars().any(char::is_control)
                    {
                        out.push(violation(
                            "project-requirement-source-invalid",
                            format!(
                                "requirement `{}` field `{field}` is not a canonical nonempty identity",
                                limit.qoi
                            ),
                            format!(
                                "state a nonempty, trim-canonical, control-free `{field}` value"
                            ),
                        ));
                    }
                }
                if !(limit.safety_factor.factor.is_finite() && limit.safety_factor.factor >= 1.0) {
                    out.push(violation(
                        "project-requirement-safety-factor",
                        format!(
                            "requirement `{}` has safety factor {}",
                            limit.qoi, limit.safety_factor.factor
                        ),
                        "state the finite factor of at least one already reflected in the effective limit",
                    ));
                }
            }
        }
        if let Some(versions) = &self.versions {
            if versions.schema != crate::FSIM_VERSION {
                out.push(violation(
                    "project-schema-version-mismatch",
                    format!(
                        "versions.schema is {} but this reader admits only {}",
                        versions.schema,
                        crate::FSIM_VERSION
                    ),
                    "run the explicit migration path; version semantics are never rewritten implicitly",
                ));
            }
            if versions.constellation.trim().is_empty() || versions.workspace.trim().is_empty() {
                out.push(violation(
                    "project-versions-field-empty",
                    "versions.constellation or versions.workspace is empty",
                    "pin both digests; replay without pins is not replay",
                ));
            }
        }
        if let Some(units) = &self.units
            && units.storage != "si-base"
        {
            out.push(violation(
                "project-units-storage",
                format!("units.storage is `{}`", units.storage),
                "the only admitted storage convention is `si-base`; quantities are stored in coherent SI base units",
            ));
        }
        if let Some(solver) = &self.solver {
            if !(solver.tolerance_rel.is_finite() && solver.tolerance_rel > 0.0) {
                out.push(violation(
                    "project-solver-tolerance",
                    format!("solver.tolerance-rel is {}", solver.tolerance_rel),
                    "state a finite positive relative tolerance",
                ));
            }
            if solver.fidelity.trim().is_empty() {
                out.push(violation(
                    "project-solver-fidelity",
                    "solver.fidelity is empty",
                    "state `auto` or an explicit registered capability id",
                ));
            }
        }
        if let Some(outputs) = &self.outputs {
            for output in outputs {
                if !matches!(output.kind.as_str(), "scalar" | "field" | "report") {
                    out.push(violation(
                        "project-output-kind",
                        format!("output `{}` has kind `{}`", output.name, output.kind),
                        "use one of `scalar`, `field`, `report`",
                    ));
                }
            }
        }
        if let Some(budgets) = &self.budgets
            && !(budgets.accuracy_rel.is_finite() && budgets.accuracy_rel > 0.0)
        {
            out.push(violation(
                "project-budget-accuracy",
                format!("budgets.accuracy-rel is {}", budgets.accuracy_rel),
                "state a finite positive relative accuracy target",
            ));
        }
        if let Some(assignments) = &self.assignments {
            for assignment in assignments {
                for (field, value) in [
                    ("artifact", assignment.artifact.as_str()),
                    ("target", assignment.target.as_str()),
                    ("length-unit", assignment.length_unit.as_str()),
                ] {
                    if value.is_empty()
                        || value.trim() != value
                        || value.chars().any(char::is_control)
                    {
                        out.push(violation(
                            "project-assignment-field-invalid",
                            format!(
                                "assignment field `{field}` value {value:?} is not a canonical nonempty label"
                            ),
                            format!(
                                "state a nonempty, trim-canonical, control-free `{field}` value"
                            ),
                        ));
                    }
                }
            }
        }
    }

    fn check_quantities(&self, out: &mut Vec<Violation>) {
        if let Some(budgets) = &self.budgets {
            check_dims(
                out,
                "project-budget-dims",
                "budgets.solve-time",
                budgets.solve_time,
                dims::TIME,
            );
        }
        if let Some(power) = &self.power {
            for row in power {
                check_dims(
                    out,
                    "project-power-dims",
                    &format!("power for `{}`", row.region),
                    row.watts,
                    dims::POWER,
                );
                if !(row.duty.is_finite() && (0.0..=1.0).contains(&row.duty)) {
                    out.push(violation(
                        "project-duty-range",
                        format!("duty for `{}` is {}", row.region, row.duty),
                        "duty must lie in 0.0..=1.0",
                    ));
                }
            }
        }
        if let Some(cooling) = &self.cooling {
            for fan in &cooling.fans {
                check_dims(
                    out,
                    "project-fan-dims",
                    &format!("fan `{}` flow", fan.name),
                    fan.flow,
                    dims::VOLUMETRIC_FLOW,
                );
                check_dims(
                    out,
                    "project-fan-dims",
                    &format!("fan `{}` static pressure", fan.name),
                    fan.static_pressure,
                    dims::PRESSURE,
                );
            }
            for vent in &cooling.vents {
                check_dims(
                    out,
                    "project-vent-dims",
                    &format!("vent on `{}`", vent.region),
                    vent.area,
                    dims::AREA,
                );
            }
            check_dims(
                out,
                "project-leakage-dims",
                "cooling.leakage",
                cooling.leakage,
                dims::POWER,
            );
        }
        self.check_range_quantities(out);
    }

    fn check_range_quantities(&self, out: &mut Vec<Violation>) {
        if let Some(envelope) = &self.envelope {
            check_dims(
                out,
                "project-envelope-dims",
                "envelope.ambient-lo",
                envelope.ambient_lo,
                dims::TEMPERATURE,
            );
            check_dims(
                out,
                "project-envelope-dims",
                "envelope.ambient-hi",
                envelope.ambient_hi,
                dims::TEMPERATURE,
            );
            check_dims(
                out,
                "project-envelope-dims",
                "envelope.pressure",
                envelope.pressure,
                dims::PRESSURE,
            );
            if envelope.ambient_lo.dims == dims::TEMPERATURE
                && envelope.ambient_hi.dims == dims::TEMPERATURE
                && envelope.ambient_lo.value > envelope.ambient_hi.value
            {
                out.push(violation(
                    "project-envelope-range",
                    "envelope.ambient-lo exceeds envelope.ambient-hi",
                    "state the ambient range low..high",
                ));
            }
        }
        if let Some(requirements) = &self.requirements {
            for limit in requirements {
                check_dims(
                    out,
                    "project-limit-dims",
                    &format!("limit for `{}`", limit.region),
                    limit.limit,
                    dims::TEMPERATURE,
                );
                check_dims(
                    out,
                    "project-limit-dims",
                    &format!("margin for `{}`", limit.region),
                    limit.margin,
                    dims::TEMPERATURE,
                );
                if !limit.limit.value.is_finite()
                    || !limit.margin.value.is_finite()
                    || limit.margin.value < 0.0
                {
                    out.push(violation(
                        "project-limit-range",
                        format!(
                            "requirement `{}` has effective limit {} and margin {}",
                            limit.qoi, limit.limit.value, limit.margin.value
                        ),
                        "state a finite effective limit and a finite nonnegative guard margin",
                    ));
                }
            }
        }
        self.check_card_bindings(out);
    }

    fn check_card_bindings(&self, out: &mut Vec<Violation>) {
        if let Some(materials) = &self.materials {
            for binding in materials {
                check_dims(
                    out,
                    "project-material-dims",
                    &format!("temp-lo for `{}`", binding.region),
                    binding.temp_lo,
                    dims::TEMPERATURE,
                );
                check_dims(
                    out,
                    "project-material-dims",
                    &format!("temp-hi for `{}`", binding.region),
                    binding.temp_hi,
                    dims::TEMPERATURE,
                );
                if binding.temp_lo.dims == dims::TEMPERATURE
                    && binding.temp_hi.dims == dims::TEMPERATURE
                    && binding.temp_lo.value > binding.temp_hi.value
                {
                    out.push(violation(
                        "project-material-range",
                        format!("material range for `{}` is inverted", binding.region),
                        "state the admitted temperature range low..high",
                    ));
                }
                if !is_hex_digest(&binding.card) {
                    out.push(violation(
                        "project-material-card",
                        format!(
                            "card reference for `{}` is not a 64-hex content hash",
                            binding.region
                        ),
                        "reference the matdb card by its full content hash",
                    ));
                }
                if let Some(claim) = &binding.claim
                    && !is_hex_digest(claim)
                {
                    out.push(violation(
                        "project-material-claim",
                        format!(
                            "claim pin for `{}` is not a 64-hex claim content hash",
                            binding.region
                        ),
                        "pin the exact claim by its full content hash, or drop the pin if the card has no conflicting claims",
                    ));
                }
            }
        }
        if let Some(interface_cards) = &self.interface_cards {
            for binding in interface_cards {
                if !is_hex_digest(&binding.card) {
                    out.push(violation(
                        "project-interface-card",
                        format!(
                            "interface card for `{}` is not a 64-hex content hash",
                            binding.interface
                        ),
                        "reference the matdb interface-system card by its full content hash",
                    ));
                }
                if let Some(claim) = &binding.claim
                    && !is_hex_digest(claim)
                {
                    out.push(violation(
                        "project-interface-claim",
                        format!(
                            "claim pin for `{}` is not a 64-hex claim content hash",
                            binding.interface
                        ),
                        "pin the exact claim by its full content hash, or drop the pin if the card has no conflicting claims",
                    ));
                }
            }
        }
    }

    /// Recompute persistent identities from the declarations, verifying any
    /// expected tokens. Returned map is declared-name -> identity.
    pub fn resolve_entities(&self, out: &mut Vec<Violation>) -> BTreeMap<String, EntityId> {
        let mut ids: BTreeMap<String, EntityId> = BTreeMap::new();
        let mut kinds: BTreeMap<String, &'static str> = BTreeMap::new();
        let Some(assembly) = &self.assembly else {
            return ids;
        };
        for decl in assembly {
            if ids.contains_key(decl.name()) {
                out.push(violation(
                    "project-entity-duplicate",
                    format!("entity name `{}` is declared twice", decl.name()),
                    "declared names are the project-local namespace and must be unique",
                ));
                continue;
            }
            let declaration = match decl {
                EntityDecl::Assembly { name, display, .. } => Some((
                    EntityDeclaration::assembly(name).with_display_name(display),
                    "assembly",
                )),
                EntityDecl::Part {
                    parent,
                    name,
                    display,
                    ..
                } => {
                    if let (Some(parent_id), Some(&"assembly")) =
                        (ids.get(parent), kinds.get(parent))
                    {
                        Some((
                            EntityDeclaration::part(*parent_id, name).with_display_name(display),
                            "part",
                        ))
                    } else {
                        out.push(violation(
                            "project-entity-parent-unknown",
                            format!("part `{name}` names parent assembly `{parent}`, which is not declared above it"),
                            "declare parents before children; a part's parent must be an assembly",
                        ));
                        None
                    }
                }
                EntityDecl::Region {
                    parent,
                    name,
                    display,
                    ..
                } => {
                    if let (Some(parent_id), Some(&"part")) = (ids.get(parent), kinds.get(parent)) {
                        Some((
                            EntityDeclaration::region(*parent_id, name).with_display_name(display),
                            "region",
                        ))
                    } else {
                        out.push(violation(
                            "project-entity-parent-unknown",
                            format!("region `{name}` names parent part `{parent}`, which is not declared above it"),
                            "declare parents before children; a region's parent must be a part",
                        ));
                        None
                    }
                }
                EntityDecl::Interface {
                    parent,
                    name,
                    display,
                    from,
                    to,
                    ..
                } => interface_declaration(&ids, &kinds, parent, name, display, from, to, out),
            };
            if let Some((declaration, kind)) = declaration {
                let id = declaration.identity();
                if let Some(expected) = decl.expect_id() {
                    let token = id.token();
                    if token != expected {
                        out.push(violation(
                            "project-entity-id-mismatch",
                            format!(
                                "entity `{}` recomputes identity `{token}` but the project pins `{expected}`",
                                decl.name()
                            ),
                            "an identity pin proves byte-equal derivation inputs; if the rename/re-import was intentional, record the rebind and update the pin",
                        ));
                    }
                }
                ids.insert(decl.name().to_string(), id);
                kinds.insert(decl.name().to_string(), kind);
            }
        }
        ids
    }

    fn check_references(&self, ids: &BTreeMap<String, EntityId>, out: &mut Vec<Violation>) {
        if self.assembly.is_none() {
            return;
        }
        fn check_ref(
            ids: &BTreeMap<String, EntityId>,
            out: &mut Vec<Violation>,
            context: String,
            name: &str,
        ) {
            if !ids.contains_key(name) {
                out.push(violation(
                    "project-ref-unknown",
                    format!("{context} references `{name}`, which resolves to no declared entity"),
                    "reference entities by their declared names from the assembly section",
                ));
            }
        }
        if let Some(materials) = &self.materials {
            for binding in materials {
                check_ref(
                    ids,
                    out,
                    format!("material binding `{}`", binding.card),
                    &binding.region,
                );
            }
        }
        if let Some(power) = &self.power {
            for row in power {
                check_ref(ids, out, "power map row".to_string(), &row.region);
            }
        }
        if let Some(cooling) = &self.cooling {
            for vent in &cooling.vents {
                check_ref(
                    ids,
                    out,
                    format!("vent (area {})", vent.area.value),
                    &vent.region,
                );
            }
        }
        if let Some(requirements) = &self.requirements {
            for limit in requirements {
                check_ref(
                    ids,
                    out,
                    format!("requirement `{}`", limit.class),
                    &limit.region,
                );
            }
        }
        if let Some(interface_cards) = &self.interface_cards {
            for binding in interface_cards {
                check_ref(
                    ids,
                    out,
                    format!("interface card `{}`", binding.card),
                    &binding.interface,
                );
            }
        }
        if let Some(assignments) = &self.assignments {
            let geometry_roles: BTreeSet<&str> = self
                .geometry
                .as_deref()
                .unwrap_or_default()
                .iter()
                .map(|artifact| artifact.role.as_str())
                .collect();
            for assignment in assignments {
                check_ref(
                    ids,
                    out,
                    format!("geometry assignment on `{}`", assignment.artifact),
                    &assignment.target,
                );
                if !geometry_roles.contains(assignment.artifact.as_str()) {
                    out.push(violation(
                        "project-assignment-artifact-unknown",
                        format!(
                            "assignment for `{}` references geometry role `{}`, which is not declared",
                            assignment.target, assignment.artifact
                        ),
                        "reference the exact role of one declared geometry artifact",
                    ));
                }
                if let Some(id) = ids.get(&assignment.target)
                    && !matches!(id.kind(), EntityKind::Region | EntityKind::Interface)
                {
                    out.push(violation(
                        "project-assignment-target-kind",
                        format!(
                            "assignment target `{}` is a {}, not a region or interface",
                            assignment.target,
                            id.kind().label()
                        ),
                        "assign finite mesh faces only to declared regions or interfaces",
                    ));
                }
            }
        }
    }

    fn check_assignment_coverage(
        &self,
        ids: &BTreeMap<String, EntityId>,
        out: &mut Vec<Violation>,
    ) {
        let Some(assignments) = &self.assignments else {
            return;
        };
        let mut targets = BTreeSet::new();
        let mut artifacts = BTreeSet::new();
        for assignment in assignments {
            if !targets.insert(assignment.target.as_str()) {
                out.push(violation(
                    "project-assignment-target-duplicate",
                    format!(
                        "entity `{}` has more than one geometry assignment",
                        assignment.target
                    ),
                    "declare exactly one artifact and selector for each region or interface",
                ));
            }
            artifacts.insert(assignment.artifact.as_str());
        }
        for (name, id) in ids {
            if matches!(id.kind(), EntityKind::Region | EntityKind::Interface)
                && !targets.contains(name.as_str())
            {
                out.push(violation(
                    "project-assignment-target-unbound",
                    format!("{} `{name}` has no geometry assignment", id.kind().label()),
                    "declare one mesh-index-free named or geometric selector for this entity",
                ));
            }
        }
        if let Some(geometry) = &self.geometry {
            let mut roles = BTreeSet::new();
            for artifact in geometry {
                if !roles.insert(artifact.role.as_str()) {
                    out.push(violation(
                        "project-geometry-role-duplicate",
                        format!(
                            "geometry role `{}` is declared more than once",
                            artifact.role
                        ),
                        "give each imported geometry artifact a unique project role",
                    ));
                }
                if !artifacts.contains(artifact.role.as_str()) {
                    out.push(violation(
                        "project-assignment-artifact-unassigned",
                        format!(
                            "geometry role `{}` has no region/interface assignments",
                            artifact.role
                        ),
                        "declare at least one selector for every imported geometry artifact, or remove the unused artifact",
                    ));
                }
            }
        }
    }
}
