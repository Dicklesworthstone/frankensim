//! Wire spellings for the `.fsim` project: one semantic model, two concrete
//! syntaxes (canonical s-expression and isomorphic JSON), both over `fs-ir`'s
//! typed AST, exactly as precedented there. Canonical bytes — the input to
//! [`canonical_hash`] — are the checked s-expression render of the envelope;
//! the JSON spelling parses to the same AST and therefore the same canonical
//! bytes, so the hash is stable across spellings by construction.
//!
//! The strict parsers refuse noncanonical input with a first-difference span.
//! The lenient s-expression parser admits two documented relaxations — an
//! omitted power `:duty` (defaulted to `1.0`) and noncanonical whitespace or
//! clause spelling — and reports both loudly: every applied default becomes a
//! [`DefaultReceipt`] and any byte difference becomes a
//! [`CanonicalizationReceipt`] carrying both hashes. No silent defaults, no
//! silent re-emission.

use fs_blake3::{ContentHash, hash_domain};
use fs_io::{HalfSpaceSide, MeshSelector};
use fs_ir::{Node, NodeKind};
use fs_qty::QtyAny;
use fs_scenario::Violation;

use crate::FSIM_VERSION;
use crate::spec::{
    Budgets, ConsequenceClass, Cooling, DecisionGate, DefaultReceipt, EntityDecl, Envelope, Fan,
    GeometryArtifact, GeometryAssignment, InterfaceCardBinding, MaterialBinding, Metadata,
    OutputRequest, PowerDissipation, ProjectSpec, RequirementDirection, RequirementSeverity,
    RequirementSource, RequirementSourceKind, SafetyFactorPolicy, Seeds, SolverSettings,
    ThermalLimit, UnitsDoctrine, Vent, Versions,
};

/// Domain for canonical `.fsim` byte hashing.
pub const FSIM_CANONICAL_DOMAIN: &str = "org.frankensim.fs-project.canonical.v1";

/// Hash canonical `.fsim` bytes under the schema's domain.
#[must_use]
pub fn canonical_hash(bytes: &[u8]) -> ContentHash {
    hash_domain(FSIM_CANONICAL_DOMAIN, bytes)
}

const fn decision_gate_tag(value: DecisionGate) -> &'static str {
    match value {
        DecisionGate::ScopingEstimate => "scoping-estimate",
        DecisionGate::DesignSelection => "design-selection",
        DecisionGate::ComplianceSignoff => "compliance-signoff",
    }
}

const fn consequence_tag(value: ConsequenceClass) -> &'static str {
    match value {
        ConsequenceClass::Advisory => "advisory",
        ConsequenceClass::Reliability => "reliability",
        ConsequenceClass::SafetyCritical => "safety-critical",
    }
}

const fn requirement_direction_tag(value: RequirementDirection) -> &'static str {
    match value {
        RequirementDirection::AtMost => "at-most",
        RequirementDirection::AtLeast => "at-least",
    }
}

const fn requirement_source_kind_tag(value: RequirementSourceKind) -> &'static str {
    match value {
        RequirementSourceKind::Standard => "standard",
        RequirementSourceKind::Datasheet => "datasheet",
        RequirementSourceKind::InternalPolicy => "internal-policy",
        RequirementSourceKind::UserDeclaration => "user-declaration",
    }
}

const fn requirement_severity_tag(value: RequirementSeverity) -> &'static str {
    match value {
        RequirementSeverity::ReliabilityDerating => "reliability-derating",
        RequirementSeverity::DamageLimit => "damage-limit",
        RequirementSeverity::SafetyCritical => "safety-critical",
    }
}

/// Typed refusal from the wire layer (syntax, envelope, canonicality).
/// Validation findings are [`Violation`]s, not errors: recognition is
/// deliberately lenient so every omission can be named at once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectError {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// What went wrong, with context.
    pub detail: String,
    /// How to fix it.
    pub hint: String,
}

impl ProjectError {
    fn syntax(error: &fs_ir::IrError) -> Self {
        ProjectError {
            code: "fsim-syntax",
            detail: error.to_string(),
            hint: "fix the underlying syntax error; the .fsim wire grammar is fs-ir's".to_string(),
        }
    }
}

impl core::fmt::Display for ProjectError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}: {} ({})", self.code, self.detail, self.hint)
    }
}

/// Receipt issued when accepted input bytes were not the canonical spelling
/// and the reader re-emitted them canonically. Mirrors fs-scenario's
/// `SourceCanonicalizationReceipt`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CanonicalizationReceipt {
    /// Hash of the exact accepted source bytes.
    pub source_hash: ContentHash,
    /// Hash of the canonical re-emission.
    pub canonical_hash: ContentHash,
}

impl CanonicalizationReceipt {
    /// Re-verify this receipt against the exact byte strings it binds.
    #[must_use]
    pub fn verifies(&self, source: &[u8], canonical: &[u8]) -> bool {
        canonical_hash(source) == self.source_hash
            && canonical_hash(canonical) == self.canonical_hash
    }
}

/// A decoded project: the spec, the canonical bytes it hashes to, and every
/// receipt the decode path owed.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedProject {
    /// The recognized project.
    pub spec: ProjectSpec,
    /// Recognition-level findings (unknown fields, malformed clauses).
    pub recognition: Vec<Violation>,
    /// Receipted defaults the lenient path applied.
    pub defaults: Vec<DefaultReceipt>,
    /// Present exactly when the accepted bytes were not canonical.
    pub canonicalization: Option<CanonicalizationReceipt>,
    /// The canonical s-expression bytes of the decoded project.
    pub canonical: String,
}

impl DecodedProject {
    /// Hash of the canonical bytes.
    #[must_use]
    pub fn hash(&self) -> ContentHash {
        canonical_hash(self.canonical.as_bytes())
    }

    /// All findings: recognition-level plus semantic validation.
    #[must_use]
    pub fn findings(&self) -> Vec<Violation> {
        let mut out = self.recognition.clone();
        out.extend(self.spec.validate());
        out
    }
}

fn sym(name: &str) -> Node {
    Node::synthetic(NodeKind::Symbol(name.to_string()))
}

fn kw(name: &str) -> Node {
    Node::synthetic(NodeKind::Keyword(name.to_string()))
}

fn text(value: &str) -> Node {
    Node::synthetic(NodeKind::Str(value.to_string()))
}

fn int(value: i64) -> Node {
    Node::synthetic(NodeKind::Int(value))
}

fn float(value: f64) -> Node {
    Node::synthetic(NodeKind::Float(value))
}

fn seed(value: u64) -> Node {
    Node::synthetic(NodeKind::Seed(value))
}

fn boolean(value: bool) -> Node {
    sym(if value { "true" } else { "false" })
}

fn list(items: Vec<Node>) -> Node {
    Node::synthetic(NodeKind::List(items))
}

fn qty(value: QtyAny) -> Result<Node, ProjectError> {
    Node::quantity(value.value, value.dims).map_err(|error| ProjectError {
        code: "fsim-quantity",
        detail: error.to_string(),
        hint: "quantities must be finite and canonically spellable".to_string(),
    })
}

/// Lower a project to the envelope AST. Absent sections are simply not
/// emitted; lowering is total over whatever is present so a draft can be
/// rendered for inspection.
pub fn lower(spec: &ProjectSpec) -> Result<Node, ProjectError> {
    let mut sections: Vec<Node> = vec![
        sym("fsim-project"),
        kw("version"),
        int(i64::from(FSIM_VERSION)),
    ];
    lower_pillars(spec, &mut sections)?;
    lower_structure(spec, &mut sections)?;
    lower_operations(spec, &mut sections)?;
    Ok(list(sections))
}

fn lower_pillars(spec: &ProjectSpec, sections: &mut Vec<Node>) -> Result<(), ProjectError> {
    if let Some(m) = &spec.metadata {
        sections.push(list(vec![
            sym("metadata"),
            kw("name"),
            text(&m.name),
            kw("created"),
            text(&m.created),
            kw("context-of-use"),
            text(&m.context_of_use),
            kw("intended-decision"),
            text(&m.intended_decision),
            kw("decision-gate"),
            text(decision_gate_tag(m.decision_gate)),
            kw("consequence"),
            text(consequence_tag(m.consequence)),
        ]));
    }
    if let Some(v) = &spec.versions {
        sections.push(list(vec![
            sym("versions"),
            kw("schema"),
            int(i64::from(v.schema)),
            kw("constellation"),
            text(&v.constellation),
            kw("workspace"),
            text(&v.workspace),
        ]));
    }
    if let Some(s) = &spec.seeds {
        sections.push(list(vec![sym("seeds"), kw("master"), seed(s.master)]));
    }
    if let Some(b) = &spec.budgets {
        let memory = i64::try_from(b.memory_bytes).map_err(|_| ProjectError {
            code: "fsim-budget-overflow",
            detail: format!(
                "memory budget {} exceeds the wire integer range",
                b.memory_bytes
            ),
            hint: "state the memory budget in bytes below 2^63".to_string(),
        })?;
        sections.push(list(vec![
            sym("budgets"),
            kw("solve-time"),
            qty(b.solve_time)?,
            kw("memory-bytes"),
            int(memory),
            kw("accuracy-rel"),
            float(b.accuracy_rel),
        ]));
    }
    if let Some(caps) = &spec.capabilities {
        let mut items = vec![sym("capabilities")];
        items.extend(caps.iter().map(|c| text(c)));
        sections.push(list(items));
    }
    if let Some(u) = &spec.units {
        sections.push(list(vec![
            sym("units"),
            kw("storage"),
            text(&u.storage),
            kw("display"),
            text(&u.display),
        ]));
    }
    Ok(())
}

fn lower_structure(spec: &ProjectSpec, sections: &mut Vec<Node>) -> Result<(), ProjectError> {
    if let Some(geometry) = &spec.geometry {
        let mut items = vec![sym("geometry")];
        for artifact in geometry {
            items.push(list(vec![
                sym("artifact"),
                kw("role"),
                text(&artifact.role),
                kw("format"),
                text(&artifact.format),
                kw("source-hash"),
                text(&format!("{:016x}", artifact.source_hash)),
                kw("parser"),
                text(&artifact.parser_version),
            ]));
        }
        sections.push(list(items));
    }
    if let Some(assignments) = &spec.assignments {
        let mut items = vec![sym("assignments")];
        for assignment in assignments {
            items.push(list(vec![
                sym("assign"),
                kw("artifact"),
                text(&assignment.artifact),
                kw("target"),
                text(&assignment.target),
                kw("length-unit"),
                text(&assignment.length_unit),
                kw("selector"),
                lower_selector(&assignment.selector),
                kw("allow-overlap"),
                boolean(assignment.allow_overlap),
            ]));
        }
        sections.push(list(items));
    }
    if let Some(assembly) = &spec.assembly {
        let mut items = vec![sym("assembly")];
        for decl in assembly {
            items.push(lower_entity(decl));
        }
        sections.push(list(items));
    }
    if let Some(materials) = &spec.materials {
        let mut items = vec![sym("materials")];
        for binding in materials {
            let mut row = vec![
                sym("binding"),
                kw("region"),
                text(&binding.region),
                kw("card"),
                text(&binding.card),
            ];
            if let Some(claim) = &binding.claim {
                row.push(kw("claim"));
                row.push(text(claim));
            }
            row.extend([
                kw("state"),
                text(&binding.state),
                kw("temp-lo"),
                qty(binding.temp_lo)?,
                kw("temp-hi"),
                qty(binding.temp_hi)?,
                kw("source"),
                text(&binding.source),
            ]);
            items.push(list(row));
        }
        sections.push(list(items));
    }
    if let Some(interface_cards) = &spec.interface_cards {
        let mut items = vec![sym("interface-cards")];
        for binding in interface_cards {
            let mut row = vec![
                sym("card"),
                kw("interface"),
                text(&binding.interface),
                kw("card"),
                text(&binding.card),
            ];
            if let Some(claim) = &binding.claim {
                row.push(kw("claim"));
                row.push(text(claim));
            }
            row.extend([kw("source"), text(&binding.source)]);
            items.push(list(row));
        }
        sections.push(list(items));
    }
    if let Some(power) = &spec.power {
        let mut items = vec![sym("power")];
        for row in power {
            items.push(list(vec![
                sym("dissipation"),
                kw("region"),
                text(&row.region),
                kw("watts"),
                qty(row.watts)?,
                kw("duty"),
                float(row.duty),
            ]));
        }
        sections.push(list(items));
    }
    Ok(())
}

fn lower_operations(spec: &ProjectSpec, sections: &mut Vec<Node>) -> Result<(), ProjectError> {
    if let Some(cooling) = &spec.cooling {
        let mut fans = vec![sym("fans")];
        for fan in &cooling.fans {
            fans.push(list(vec![
                sym("fan"),
                kw("name"),
                text(&fan.name),
                kw("flow"),
                qty(fan.flow)?,
                kw("static-pressure"),
                qty(fan.static_pressure)?,
            ]));
        }
        let mut vents = vec![sym("vents")];
        for vent in &cooling.vents {
            vents.push(list(vec![
                sym("vent"),
                kw("region"),
                text(&vent.region),
                kw("area"),
                qty(vent.area)?,
            ]));
        }
        sections.push(list(vec![
            sym("cooling"),
            list(fans),
            list(vents),
            list(vec![sym("leakage"), kw("watts"), qty(cooling.leakage)?]),
        ]));
    }
    if let Some(envelope) = &spec.envelope {
        sections.push(list(vec![
            sym("envelope"),
            kw("ambient-lo"),
            qty(envelope.ambient_lo)?,
            kw("ambient-hi"),
            qty(envelope.ambient_hi)?,
            kw("pressure"),
            qty(envelope.pressure)?,
        ]));
    }
    if let Some(requirements) = &spec.requirements {
        let mut items = vec![sym("requirements")];
        for limit in requirements {
            items.push(list(vec![
                sym("t-limit"),
                kw("qoi"),
                text(&limit.qoi),
                kw("class"),
                text(&limit.class),
                kw("region"),
                text(&limit.region),
                kw("direction"),
                text(requirement_direction_tag(limit.direction)),
                kw("limit"),
                qty(limit.limit)?,
                kw("margin"),
                qty(limit.margin)?,
                kw("source-kind"),
                text(requirement_source_kind_tag(limit.source.kind)),
                kw("source-document"),
                text(&limit.source.document),
                kw("source-version"),
                text(&limit.source.version),
                kw("source-locator"),
                text(&limit.source.locator),
                kw("safety-factor"),
                float(limit.safety_factor.factor),
                kw("safety-source-kind"),
                text(requirement_source_kind_tag(limit.safety_factor.source.kind)),
                kw("safety-source-document"),
                text(&limit.safety_factor.source.document),
                kw("safety-source-version"),
                text(&limit.safety_factor.source.version),
                kw("safety-source-locator"),
                text(&limit.safety_factor.source.locator),
                kw("severity"),
                text(requirement_severity_tag(limit.severity)),
            ]));
        }
        sections.push(list(items));
    }
    if let Some(solver) = &spec.solver {
        sections.push(list(vec![
            sym("solver"),
            kw("fidelity"),
            text(&solver.fidelity),
            kw("tolerance-rel"),
            float(solver.tolerance_rel),
        ]));
    }
    if let Some(outputs) = &spec.outputs {
        let mut items = vec![sym("outputs")];
        for output in outputs {
            items.push(list(vec![
                sym("qoi"),
                kw("name"),
                text(&output.name),
                kw("kind"),
                text(&output.kind),
            ]));
        }
        sections.push(list(items));
    }
    Ok(())
}

fn lower_entity(decl: &EntityDecl) -> Node {
    let mut items = match decl {
        EntityDecl::Assembly { name, display, .. } => vec![
            sym("assembly-decl"),
            kw("name"),
            text(name),
            kw("display"),
            text(display),
        ],
        EntityDecl::Part {
            parent,
            name,
            display,
            ..
        } => vec![
            sym("part"),
            kw("parent"),
            text(parent),
            kw("name"),
            text(name),
            kw("display"),
            text(display),
        ],
        EntityDecl::Region {
            parent,
            name,
            display,
            ..
        } => vec![
            sym("region"),
            kw("parent"),
            text(parent),
            kw("name"),
            text(name),
            kw("display"),
            text(display),
        ],
        EntityDecl::Interface {
            parent,
            name,
            display,
            from,
            to,
            ..
        } => vec![
            sym("interface"),
            kw("parent"),
            text(parent),
            kw("name"),
            text(name),
            kw("display"),
            text(display),
            kw("from"),
            text(from),
            kw("to"),
            text(to),
        ],
    };
    let expect = match decl {
        EntityDecl::Assembly { expect_id, .. }
        | EntityDecl::Part { expect_id, .. }
        | EntityDecl::Region { expect_id, .. }
        | EntityDecl::Interface { expect_id, .. } => expect_id,
    };
    if let Some(expected) = expect {
        items.push(kw("id"));
        items.push(text(expected));
    }
    list(items)
}

fn vec3(value: [f64; 3]) -> Node {
    list(vec![
        sym("vec3"),
        float(value[0]),
        float(value[1]),
        float(value[2]),
    ])
}

fn lower_selector(selector: &MeshSelector) -> Node {
    match selector {
        MeshSelector::NamedGroup { name } => list(vec![sym("named-group"), kw("name"), text(name)]),
        MeshSelector::HalfSpace {
            normal,
            offset,
            side,
            tolerance,
        } => list(vec![
            sym("half-space"),
            kw("normal"),
            vec3(*normal),
            kw("offset"),
            float(*offset),
            kw("side"),
            sym(match side {
                HalfSpaceSide::AtMost => "at-most",
                HalfSpaceSide::AtLeast => "at-least",
            }),
            kw("tolerance"),
            float(*tolerance),
        ]),
        MeshSelector::Box {
            min,
            max,
            tolerance,
        } => list(vec![
            sym("box"),
            kw("min"),
            vec3(*min),
            kw("max"),
            vec3(*max),
            kw("tolerance"),
            float(*tolerance),
        ]),
        MeshSelector::Cylinder {
            origin,
            axis,
            radius,
            axial_min,
            axial_max,
            tolerance,
        } => list(vec![
            sym("cylinder"),
            kw("origin"),
            vec3(*origin),
            kw("axis"),
            vec3(*axis),
            kw("radius"),
            float(*radius),
            kw("axial-min"),
            float(*axial_min),
            kw("axial-max"),
            float(*axial_max),
            kw("tolerance"),
            float(*tolerance),
        ]),
        MeshSelector::NearestDatum {
            point,
            max_distance,
            tolerance,
        } => list(vec![
            sym("nearest-datum"),
            kw("point"),
            vec3(*point),
            kw("max-distance"),
            float(*max_distance),
            kw("tolerance"),
            float(*tolerance),
        ]),
        MeshSelector::ExplicitFaceSet {
            faces,
            fragility_acknowledged,
        } => {
            let mut face_items = vec![sym("faces")];
            face_items.extend(faces.iter().map(|face| int(i64::from(*face))));
            list(vec![
                sym("explicit-face-set"),
                kw("faces"),
                list(face_items),
                kw("fragility-acknowledged"),
                boolean(*fragility_acknowledged),
            ])
        }
    }
}

/// Render the canonical s-expression bytes for a project.
pub fn print_sexpr(spec: &ProjectSpec) -> Result<String, ProjectError> {
    let node = lower(spec)?;
    fs_ir::sexpr::print_checked(&node).map_err(|e| ProjectError::syntax(&e))
}

/// Render the isomorphic JSON spelling.
pub fn print_json(spec: &ProjectSpec) -> Result<String, ProjectError> {
    let node = lower(spec)?;
    fs_ir::json::print_checked(&node).map_err(|e| ProjectError::syntax(&e))
}

/// Strictly parse canonical s-expression bytes: any byte difference from the
/// canonical re-render refuses with the first difference position.
pub fn parse_sexpr(src: &str) -> Result<DecodedProject, ProjectError> {
    let decoded = parse_sexpr_lenient(src)?;
    require_canonical(src, &decoded)?;
    Ok(decoded)
}

/// Strictly parse the JSON spelling: the input must be the checked JSON
/// render of the decoded project's AST.
pub fn parse_json(src: &str) -> Result<DecodedProject, ProjectError> {
    let node = fs_ir::json::parse(src).map_err(|e| ProjectError::syntax(&e))?;
    let decoded = recognize(&node)?;
    let canonical_json = print_json(&decoded.spec)?;
    if src != canonical_json {
        return Err(non_canonical(src, &canonical_json, "json"));
    }
    require_no_relaxations(&decoded)?;
    Ok(decoded)
}

/// Leniently parse s-expression bytes: noncanonical spellings are accepted
/// with a [`CanonicalizationReceipt`], and the documented power-duty default
/// is applied with a [`DefaultReceipt`]. Nothing is applied silently.
pub fn parse_sexpr_lenient(src: &str) -> Result<DecodedProject, ProjectError> {
    let node = fs_ir::sexpr::parse(src).map_err(|e| ProjectError::syntax(&e))?;
    let mut decoded = recognize(&node)?;
    if src != decoded.canonical {
        decoded.canonicalization = Some(CanonicalizationReceipt {
            source_hash: canonical_hash(src.as_bytes()),
            canonical_hash: canonical_hash(decoded.canonical.as_bytes()),
        });
    }
    Ok(decoded)
}

fn require_canonical(src: &str, decoded: &DecodedProject) -> Result<(), ProjectError> {
    if src != decoded.canonical {
        return Err(non_canonical(src, &decoded.canonical, "s-expression"));
    }
    require_no_relaxations(decoded)
}

fn require_no_relaxations(decoded: &DecodedProject) -> Result<(), ProjectError> {
    if decoded.defaults.is_empty() {
        Ok(())
    } else {
        Err(ProjectError {
            code: "fsim-default-in-strict-mode",
            detail: format!(
                "the strict parser refuses applied defaults; {} field(s) were omitted",
                decoded.defaults.len()
            ),
            hint: "spell every field explicitly, or use the lenient parser and retain its receipts"
                .to_string(),
        })
    }
}

fn non_canonical(src: &str, canonical: &str, syntax: &str) -> ProjectError {
    let position = src
        .bytes()
        .zip(canonical.bytes())
        .position(|(a, b)| a != b)
        .unwrap_or_else(|| src.len().min(canonical.len()));
    ProjectError {
        code: "fsim-non-canonical",
        detail: format!("{syntax} input differs from the canonical spelling at byte {position}"),
        hint: "persist the canonical bytes; use the lenient parser to canonicalize with a receipt"
            .to_string(),
    }
}

struct Reader<'a> {
    items: &'a [Node],
    at: usize,
}

impl<'a> Reader<'a> {
    fn keyword(&mut self) -> Option<&'a str> {
        if let Some(Node {
            kind: NodeKind::Keyword(k),
            ..
        }) = self.items.get(self.at)
        {
            self.at += 1;
            Some(k)
        } else {
            None
        }
    }

    fn next(&mut self) -> Option<&'a Node> {
        let node = self.items.get(self.at);
        self.at += 1;
        node
    }
}

fn unknown_field(recognition: &mut Vec<Violation>, context: &str, field: &str) {
    recognition.push(Violation {
        code: "project-unknown-field",
        what: format!("`{context}` carries unknown field `{field}`"),
        fix: format!(
            "remove `{field}`; unknown fields are refused so typos cannot silently drop intent"
        ),
    });
}

fn expect_str(node: Option<&Node>, context: &str, out: &mut Vec<Violation>) -> String {
    if let Some(Node {
        kind: NodeKind::Str(s),
        ..
    }) = node
    {
        s.clone()
    } else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: format!("`{context}` expected a string value"),
            fix: format!("spell `{context}` as a double-quoted string"),
        });
        String::new()
    }
}

fn expect_qty(node: Option<&Node>, context: &str, out: &mut Vec<Violation>) -> QtyAny {
    if let Some(Node {
        kind: NodeKind::Qty { value, dims, .. },
        ..
    }) = node
    {
        QtyAny::new(*value, *dims)
    } else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: format!("`{context}` expected a quantity literal"),
            fix: format!("spell `{context}` as a dimensioned quantity (SI base units)"),
        });
        QtyAny::dimensionless(f64::NAN)
    }
}

fn expect_float(node: Option<&Node>, context: &str, out: &mut Vec<Violation>) -> f64 {
    match node {
        Some(Node {
            kind: NodeKind::Float(v),
            ..
        }) => *v,
        Some(Node {
            kind: NodeKind::Int(v),
            ..
        }) => {
            out.push(Violation {
                code: "project-malformed-clause",
                what: format!("`{context}` expected a float literal, got integer {v}"),
                fix: format!("spell `{context}` with a decimal point"),
            });
            #[allow(clippy::cast_precision_loss)]
            {
                *v as f64
            }
        }
        _ => {
            out.push(Violation {
                code: "project-malformed-clause",
                what: format!("`{context}` expected a float literal"),
                fix: format!("spell `{context}` as a decimal number"),
            });
            f64::NAN
        }
    }
}

fn expect_boolean(node: Option<&Node>, context: &str, out: &mut Vec<Violation>) -> bool {
    match node {
        Some(Node {
            kind: NodeKind::Symbol(value),
            ..
        }) if value == "true" => true,
        Some(Node {
            kind: NodeKind::Symbol(value),
            ..
        }) if value == "false" => false,
        _ => {
            out.push(Violation {
                code: "project-malformed-clause",
                what: format!("`{context}` expected symbol `true` or `false`"),
                fix: format!("spell `{context}` as exactly `true` or `false`"),
            });
            false
        }
    }
}

fn expect_vec3(node: Option<&Node>, context: &str, out: &mut Vec<Violation>) -> [f64; 3] {
    let Some(Node {
        kind: NodeKind::List(items),
        ..
    }) = node
    else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: format!("`{context}` expected `(vec3 <x> <y> <z>)`"),
            fix: format!("spell `{context}` as three explicit decimal coordinates"),
        });
        return [f64::NAN; 3];
    };
    if !matches!(
        items.first(),
        Some(Node {
            kind: NodeKind::Symbol(head),
            ..
        }) if head == "vec3"
    ) || items.len() != 4
    {
        out.push(Violation {
            code: "project-malformed-clause",
            what: format!("`{context}` expected exactly `(vec3 <x> <y> <z>)`"),
            fix: format!("spell `{context}` as three explicit decimal coordinates"),
        });
        return [f64::NAN; 3];
    }
    [
        expect_float(items.get(1), &format!("{context}[0]"), out),
        expect_float(items.get(2), &format!("{context}[1]"), out),
        expect_float(items.get(3), &format!("{context}[2]"), out),
    ]
}

fn section_name(node: &Node) -> Option<(&str, &[Node])> {
    if let NodeKind::List(items) = &node.kind
        && let Some(Node {
            kind: NodeKind::Symbol(s),
            ..
        }) = items.first()
    {
        return Some((s.as_str(), &items[1..]));
    }
    None
}

/// Recognize the envelope AST into a project, collecting recognition-level
/// violations rather than refusing at the first problem.
pub fn recognize(node: &Node) -> Result<DecodedProject, ProjectError> {
    let NodeKind::List(items) = &node.kind else {
        return Err(not_a_project("the document root is not a list"));
    };
    match items.first() {
        Some(Node {
            kind: NodeKind::Symbol(s),
            ..
        }) if s == "fsim-project" => {}
        _ => {
            return Err(not_a_project(
                "the document does not open with `fsim-project`",
            ));
        }
    }
    match (items.get(1), items.get(2)) {
        (
            Some(Node {
                kind: NodeKind::Keyword(k),
                ..
            }),
            Some(Node {
                kind: NodeKind::Int(v),
                ..
            }),
        ) if k == "version" => {
            if u32::try_from(*v) != Ok(FSIM_VERSION) {
                return Err(ProjectError {
                    code: "fsim-unsupported-version",
                    detail: format!("envelope declares version {v}; this reader admits only {FSIM_VERSION}"),
                    hint: "use the explicit, receipted migration path; version semantics are never rewritten implicitly".to_string(),
                });
            }
        }
        _ => {
            return Err(not_a_project(
                "the envelope lacks `:version <int>` in position",
            ));
        }
    }

    let mut spec = ProjectSpec::default();
    let mut recognition = Vec::new();
    let mut defaults = Vec::new();
    for section in &items[3..] {
        let Some((name, body)) = section_name(section) else {
            recognition.push(Violation {
                code: "project-malformed-clause",
                what: "a top-level section is not a `(name ...)` list".to_string(),
                fix: "every section is a list opening with its section symbol".to_string(),
            });
            continue;
        };
        match name {
            "metadata" => spec.metadata = Some(read_metadata(body, &mut recognition)),
            "versions" => spec.versions = Some(read_versions(body, &mut recognition)),
            "seeds" => spec.seeds = read_seeds(body, &mut recognition),
            "budgets" => spec.budgets = Some(read_budgets(body, &mut recognition)),
            "capabilities" => spec.capabilities = Some(read_capabilities(body, &mut recognition)),
            "units" => spec.units = Some(read_units(body, &mut recognition)),
            "geometry" => spec.geometry = Some(read_geometry(body, &mut recognition)),
            "assignments" => {
                spec.assignments = Some(read_assignments(body, &mut recognition));
            }
            "assembly" => spec.assembly = Some(read_assembly(body, &mut recognition)),
            "materials" => spec.materials = Some(read_materials(body, &mut recognition)),
            "interface-cards" => {
                spec.interface_cards = Some(read_interface_cards(body, &mut recognition));
            }
            "power" => spec.power = Some(read_power(body, &mut recognition, &mut defaults)),
            "cooling" => spec.cooling = read_cooling(body, &mut recognition),
            "envelope" => spec.envelope = Some(read_envelope(body, &mut recognition)),
            "requirements" => spec.requirements = Some(read_requirements(body, &mut recognition)),
            "solver" => spec.solver = Some(read_solver(body, &mut recognition)),
            "outputs" => spec.outputs = Some(read_outputs(body, &mut recognition)),
            other => unknown_field(&mut recognition, "fsim-project", other),
        }
    }

    let canonical = print_sexpr(&spec)?;
    Ok(DecodedProject {
        spec,
        recognition,
        defaults,
        canonicalization: None,
        canonical,
    })
}

fn not_a_project(detail: &str) -> ProjectError {
    ProjectError {
        code: "fsim-not-a-project",
        detail: detail.to_string(),
        hint: "a .fsim document is `(fsim-project :version N <sections>...)`".to_string(),
    }
}

fn read_pairs<'a>(
    body: &'a [Node],
    context: &str,
    known: &[&str],
    recognition: &mut Vec<Violation>,
) -> Vec<(&'a str, &'a Node)> {
    let mut reader = Reader { items: body, at: 0 };
    let mut pairs = Vec::new();
    while reader.at < body.len() {
        let Some(key) = reader.keyword() else {
            recognition.push(Violation {
                code: "project-malformed-clause",
                what: format!("`{context}` expected a `:keyword value` pair"),
                fix: format!("`{context}` fields are keyword/value pairs in canonical order"),
            });
            break;
        };
        let Some(value) = reader.next() else {
            recognition.push(Violation {
                code: "project-malformed-clause",
                what: format!("`{context}` field `:{key}` has no value"),
                fix: "every keyword takes exactly one value".to_string(),
            });
            break;
        };
        if known.contains(&key) {
            pairs.push((key, value));
        } else {
            unknown_field(recognition, context, key);
        }
    }
    pairs
}

fn field<'a>(pairs: &[(&'a str, &'a Node)], key: &str) -> Option<&'a Node> {
    pairs.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

fn invalid_enum(out: &mut Vec<Violation>, field: &'static str, value: &str, allowed: &str) {
    if value.is_empty() {
        return;
    }
    out.push(Violation {
        code: "project-field-enum",
        what: format!("`{field}` has unknown value {value:?}"),
        fix: format!("use one of {allowed}"),
    });
}

fn parse_decision_gate(value: String, out: &mut Vec<Violation>) -> DecisionGate {
    match value.as_str() {
        "scoping-estimate" => DecisionGate::ScopingEstimate,
        "design-selection" => DecisionGate::DesignSelection,
        "compliance-signoff" => DecisionGate::ComplianceSignoff,
        _ => {
            invalid_enum(
                out,
                "metadata.decision-gate",
                &value,
                "`scoping-estimate`, `design-selection`, `compliance-signoff`",
            );
            DecisionGate::ComplianceSignoff
        }
    }
}

fn parse_consequence(value: String, out: &mut Vec<Violation>) -> ConsequenceClass {
    match value.as_str() {
        "advisory" => ConsequenceClass::Advisory,
        "reliability" => ConsequenceClass::Reliability,
        "safety-critical" => ConsequenceClass::SafetyCritical,
        _ => {
            invalid_enum(
                out,
                "metadata.consequence",
                &value,
                "`advisory`, `reliability`, `safety-critical`",
            );
            ConsequenceClass::SafetyCritical
        }
    }
}

fn parse_requirement_direction(value: String, out: &mut Vec<Violation>) -> RequirementDirection {
    match value.as_str() {
        "at-most" => RequirementDirection::AtMost,
        "at-least" => RequirementDirection::AtLeast,
        _ => {
            invalid_enum(out, "t-limit.direction", &value, "`at-most`, `at-least`");
            RequirementDirection::AtMost
        }
    }
}

fn parse_requirement_source_kind(
    value: String,
    field: &'static str,
    out: &mut Vec<Violation>,
) -> RequirementSourceKind {
    match value.as_str() {
        "standard" => RequirementSourceKind::Standard,
        "datasheet" => RequirementSourceKind::Datasheet,
        "internal-policy" => RequirementSourceKind::InternalPolicy,
        "user-declaration" => RequirementSourceKind::UserDeclaration,
        _ => {
            invalid_enum(
                out,
                field,
                &value,
                "`standard`, `datasheet`, `internal-policy`, `user-declaration`",
            );
            RequirementSourceKind::UserDeclaration
        }
    }
}

fn parse_requirement_severity(value: String, out: &mut Vec<Violation>) -> RequirementSeverity {
    match value.as_str() {
        "reliability-derating" => RequirementSeverity::ReliabilityDerating,
        "damage-limit" => RequirementSeverity::DamageLimit,
        "safety-critical" => RequirementSeverity::SafetyCritical,
        _ => {
            invalid_enum(
                out,
                "t-limit.severity",
                &value,
                "`reliability-derating`, `damage-limit`, `safety-critical`",
            );
            RequirementSeverity::SafetyCritical
        }
    }
}

fn read_metadata(body: &[Node], out: &mut Vec<Violation>) -> Metadata {
    let pairs = read_pairs(
        body,
        "metadata",
        &[
            "name",
            "created",
            "context-of-use",
            "intended-decision",
            "decision-gate",
            "consequence",
        ],
        out,
    );
    Metadata {
        name: expect_str(field(&pairs, "name"), "metadata.name", out),
        created: expect_str(field(&pairs, "created"), "metadata.created", out),
        context_of_use: expect_str(
            field(&pairs, "context-of-use"),
            "metadata.context-of-use",
            out,
        ),
        intended_decision: expect_str(
            field(&pairs, "intended-decision"),
            "metadata.intended-decision",
            out,
        ),
        decision_gate: parse_decision_gate(
            expect_str(
                field(&pairs, "decision-gate"),
                "metadata.decision-gate",
                out,
            ),
            out,
        ),
        consequence: parse_consequence(
            expect_str(field(&pairs, "consequence"), "metadata.consequence", out),
            out,
        ),
    }
}

fn read_versions(body: &[Node], out: &mut Vec<Violation>) -> Versions {
    let pairs = read_pairs(
        body,
        "versions",
        &["schema", "constellation", "workspace"],
        out,
    );
    let schema = if let Some(Node {
        kind: NodeKind::Int(v),
        ..
    }) = field(&pairs, "schema")
    {
        u32::try_from(*v).unwrap_or(0)
    } else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`versions.schema` expected an integer".to_string(),
            fix: "spell `versions.schema` as the integer schema version".to_string(),
        });
        0
    };
    Versions {
        schema,
        constellation: expect_str(
            field(&pairs, "constellation"),
            "versions.constellation",
            out,
        ),
        workspace: expect_str(field(&pairs, "workspace"), "versions.workspace", out),
    }
}

fn read_seeds(body: &[Node], out: &mut Vec<Violation>) -> Option<Seeds> {
    let pairs = read_pairs(body, "seeds", &["master"], out);
    if let Some(Node {
        kind: NodeKind::Seed(v),
        ..
    }) = field(&pairs, "master")
    {
        Some(Seeds { master: *v })
    } else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`seeds.master` expected a 0x-prefixed seed literal".to_string(),
            fix: "spell the master seed as `0x...` hex fitting u64".to_string(),
        });
        None
    }
}

fn read_budgets(body: &[Node], out: &mut Vec<Violation>) -> Budgets {
    let pairs = read_pairs(
        body,
        "budgets",
        &["solve-time", "memory-bytes", "accuracy-rel"],
        out,
    );
    let memory_bytes = match field(&pairs, "memory-bytes") {
        Some(Node {
            kind: NodeKind::Int(v),
            ..
        }) if *v >= 0 => u64::try_from(*v).unwrap_or(0),
        _ => {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`budgets.memory-bytes` expected a non-negative integer".to_string(),
                fix: "state the memory budget in bytes".to_string(),
            });
            0
        }
    };
    Budgets {
        solve_time: expect_qty(field(&pairs, "solve-time"), "budgets.solve-time", out),
        memory_bytes,
        accuracy_rel: expect_float(field(&pairs, "accuracy-rel"), "budgets.accuracy-rel", out),
    }
}

fn read_capabilities(body: &[Node], out: &mut Vec<Violation>) -> Vec<String> {
    let mut capabilities = Vec::new();
    for node in body {
        match &node.kind {
            NodeKind::Str(s) => capabilities.push(s.clone()),
            _ => out.push(Violation {
                code: "project-malformed-clause",
                what: "`capabilities` entries must be strings".to_string(),
                fix: "list required capability ids as strings".to_string(),
            }),
        }
    }
    capabilities
}

fn read_units(body: &[Node], out: &mut Vec<Violation>) -> UnitsDoctrine {
    let pairs = read_pairs(body, "units", &["storage", "display"], out);
    UnitsDoctrine {
        storage: expect_str(field(&pairs, "storage"), "units.storage", out),
        display: expect_str(field(&pairs, "display"), "units.display", out),
    }
}

fn read_geometry(body: &[Node], out: &mut Vec<Violation>) -> Vec<GeometryArtifact> {
    let mut artifacts = Vec::new();
    for node in body {
        let Some(("artifact", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`geometry` rows must be `(artifact ...)`".to_string(),
                fix: "reference imported artifacts through their quarantine receipts".to_string(),
            });
            continue;
        };
        let pairs = read_pairs(
            inner,
            "artifact",
            &["role", "format", "source-hash", "parser"],
            out,
        );
        let source_hash = {
            let raw = expect_str(field(&pairs, "source-hash"), "artifact.source-hash", out);
            match u64::from_str_radix(&raw, 16) {
                Ok(v) if raw.len() == 16 => v,
                _ => {
                    out.push(Violation {
                        code: "project-malformed-clause",
                        what: format!("artifact source-hash `{raw}` is not 16 hex digits"),
                        fix: "copy the import receipt's source hash as 16 lowercase hex digits"
                            .to_string(),
                    });
                    0
                }
            }
        };
        artifacts.push(GeometryArtifact {
            role: expect_str(field(&pairs, "role"), "artifact.role", out),
            format: expect_str(field(&pairs, "format"), "artifact.format", out),
            source_hash,
            parser_version: expect_str(field(&pairs, "parser"), "artifact.parser", out),
        });
    }
    artifacts
}

fn read_assignments(body: &[Node], out: &mut Vec<Violation>) -> Vec<GeometryAssignment> {
    let mut assignments = Vec::new();
    for node in body {
        let Some(("assign", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`assignments` rows must be `(assign ...)`".to_string(),
                fix: "declare one artifact, target, length unit, selector, and overlap policy per row"
                    .to_string(),
            });
            continue;
        };
        let pairs = read_pairs(
            inner,
            "assign",
            &[
                "artifact",
                "target",
                "length-unit",
                "selector",
                "allow-overlap",
            ],
            out,
        );
        assignments.push(GeometryAssignment {
            artifact: expect_str(field(&pairs, "artifact"), "assign.artifact", out),
            target: expect_str(field(&pairs, "target"), "assign.target", out),
            length_unit: expect_str(field(&pairs, "length-unit"), "assign.length-unit", out),
            selector: read_selector(field(&pairs, "selector"), out),
            allow_overlap: expect_boolean(
                field(&pairs, "allow-overlap"),
                "assign.allow-overlap",
                out,
            ),
        });
    }
    assignments
}

fn read_selector(node: Option<&Node>, out: &mut Vec<Violation>) -> MeshSelector {
    let Some(node) = node else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`assign.selector` is missing".to_string(),
            fix: "declare one named, geometric, datum, or acknowledged explicit-face selector"
                .to_string(),
        });
        return MeshSelector::NamedGroup {
            name: String::new(),
        };
    };
    let Some((kind, body)) = section_name(node) else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`assign.selector` must be a selector list".to_string(),
            fix: "use `(named-group ...)`, `(half-space ...)`, `(box ...)`, `(cylinder ...)`, `(nearest-datum ...)`, or `(explicit-face-set ...)`".to_string(),
        });
        return MeshSelector::NamedGroup {
            name: String::new(),
        };
    };
    match kind {
        "named-group" => {
            let pairs = read_pairs(body, kind, &["name"], out);
            MeshSelector::NamedGroup {
                name: expect_str(field(&pairs, "name"), "named-group.name", out),
            }
        }
        "half-space" => {
            let pairs = read_pairs(body, kind, &["normal", "offset", "side", "tolerance"], out);
            let side = match field(&pairs, "side") {
                Some(Node {
                    kind: NodeKind::Symbol(side),
                    ..
                }) if side == "at-most" => HalfSpaceSide::AtMost,
                Some(Node {
                    kind: NodeKind::Symbol(side),
                    ..
                }) if side == "at-least" => HalfSpaceSide::AtLeast,
                _ => {
                    out.push(Violation {
                        code: "project-malformed-clause",
                        what: "`half-space.side` must be `at-most` or `at-least`".to_string(),
                        fix: "state which closed half-space the selector admits".to_string(),
                    });
                    HalfSpaceSide::AtMost
                }
            };
            MeshSelector::HalfSpace {
                normal: expect_vec3(field(&pairs, "normal"), "half-space.normal", out),
                offset: expect_float(field(&pairs, "offset"), "half-space.offset", out),
                side,
                tolerance: expect_float(field(&pairs, "tolerance"), "half-space.tolerance", out),
            }
        }
        "box" => {
            let pairs = read_pairs(body, kind, &["min", "max", "tolerance"], out);
            MeshSelector::Box {
                min: expect_vec3(field(&pairs, "min"), "box.min", out),
                max: expect_vec3(field(&pairs, "max"), "box.max", out),
                tolerance: expect_float(field(&pairs, "tolerance"), "box.tolerance", out),
            }
        }
        "cylinder" => {
            let pairs = read_pairs(
                body,
                kind,
                &[
                    "origin",
                    "axis",
                    "radius",
                    "axial-min",
                    "axial-max",
                    "tolerance",
                ],
                out,
            );
            MeshSelector::Cylinder {
                origin: expect_vec3(field(&pairs, "origin"), "cylinder.origin", out),
                axis: expect_vec3(field(&pairs, "axis"), "cylinder.axis", out),
                radius: expect_float(field(&pairs, "radius"), "cylinder.radius", out),
                axial_min: expect_float(field(&pairs, "axial-min"), "cylinder.axial-min", out),
                axial_max: expect_float(field(&pairs, "axial-max"), "cylinder.axial-max", out),
                tolerance: expect_float(field(&pairs, "tolerance"), "cylinder.tolerance", out),
            }
        }
        "nearest-datum" => {
            let pairs = read_pairs(body, kind, &["point", "max-distance", "tolerance"], out);
            MeshSelector::NearestDatum {
                point: expect_vec3(field(&pairs, "point"), "nearest-datum.point", out),
                max_distance: expect_float(
                    field(&pairs, "max-distance"),
                    "nearest-datum.max-distance",
                    out,
                ),
                tolerance: expect_float(field(&pairs, "tolerance"), "nearest-datum.tolerance", out),
            }
        }
        "explicit-face-set" => {
            let pairs = read_pairs(body, kind, &["faces", "fragility-acknowledged"], out);
            MeshSelector::ExplicitFaceSet {
                faces: read_faces(field(&pairs, "faces"), out),
                fragility_acknowledged: expect_boolean(
                    field(&pairs, "fragility-acknowledged"),
                    "explicit-face-set.fragility-acknowledged",
                    out,
                ),
            }
        }
        other => {
            unknown_field(out, "assign.selector", other);
            MeshSelector::NamedGroup {
                name: String::new(),
            }
        }
    }
}

fn read_faces(node: Option<&Node>, out: &mut Vec<Violation>) -> Vec<u32> {
    let Some(Node {
        kind: NodeKind::List(items),
        ..
    }) = node
    else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`explicit-face-set.faces` expected `(faces <ordinal>...)`".to_string(),
            fix: "list non-negative u32 face ordinals and acknowledge their remeshing fragility"
                .to_string(),
        });
        return Vec::new();
    };
    if !matches!(
        items.first(),
        Some(Node {
            kind: NodeKind::Symbol(head),
            ..
        }) if head == "faces"
    ) {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`explicit-face-set.faces` list must open with `faces`".to_string(),
            fix: "spell the value as `(faces <ordinal>...)`".to_string(),
        });
        return Vec::new();
    }
    let mut faces = Vec::new();
    for item in &items[1..] {
        let Node {
            kind: NodeKind::Int(value),
            ..
        } = item
        else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`explicit-face-set.faces` contains a non-u32 ordinal".to_string(),
                fix: "list only non-negative face ordinals that fit u32".to_string(),
            });
            continue;
        };
        match u32::try_from(*value) {
            Ok(face) => faces.push(face),
            Err(_) => out.push(Violation {
                code: "project-malformed-clause",
                what: "`explicit-face-set.faces` contains a non-u32 ordinal".to_string(),
                fix: "list only non-negative face ordinals that fit u32".to_string(),
            }),
        }
    }
    faces
}

fn read_assembly(body: &[Node], out: &mut Vec<Violation>) -> Vec<EntityDecl> {
    let mut declarations = Vec::new();
    for node in body {
        let Some((kind, inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`assembly` rows must be entity declaration lists".to_string(),
                fix: "declare `(assembly-decl ...)`, `(part ...)`, `(region ...)`, or `(interface ...)`".to_string(),
            });
            continue;
        };
        let known: &[&str] = match kind {
            "assembly-decl" => &["name", "display", "id"],
            "part" | "region" => &["parent", "name", "display", "id"],
            "interface" => &["parent", "name", "display", "from", "to", "id"],
            other => {
                unknown_field(out, "assembly", other);
                continue;
            }
        };
        let pairs = read_pairs(inner, kind, known, out);
        let name = expect_str(field(&pairs, "name"), &format!("{kind}.name"), out);
        let display = expect_str(field(&pairs, "display"), &format!("{kind}.display"), out);
        let expect_id = field(&pairs, "id").map(|node| expect_str(Some(node), "id", out));
        let declaration = match kind {
            "assembly-decl" => EntityDecl::Assembly {
                name,
                display,
                expect_id,
            },
            "part" => EntityDecl::Part {
                parent: expect_str(field(&pairs, "parent"), "part.parent", out),
                name,
                display,
                expect_id,
            },
            "region" => EntityDecl::Region {
                parent: expect_str(field(&pairs, "parent"), "region.parent", out),
                name,
                display,
                expect_id,
            },
            "interface" => EntityDecl::Interface {
                parent: expect_str(field(&pairs, "parent"), "interface.parent", out),
                name,
                display,
                from: expect_str(field(&pairs, "from"), "interface.from", out),
                to: expect_str(field(&pairs, "to"), "interface.to", out),
                expect_id,
            },
            _ => unreachable!("kind was matched above"),
        };
        declarations.push(declaration);
    }
    declarations
}

fn read_materials(body: &[Node], out: &mut Vec<Violation>) -> Vec<MaterialBinding> {
    let mut bindings = Vec::new();
    for node in body {
        let Some(("binding", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`materials` rows must be `(binding ...)`".to_string(),
                fix: "bind matdb cards with `(binding :region ... :card ... :state ... :temp-lo ... :temp-hi ... :source ...)`".to_string(),
            });
            continue;
        };
        let pairs = read_pairs(
            inner,
            "binding",
            &[
                "region", "card", "claim", "state", "temp-lo", "temp-hi", "source",
            ],
            out,
        );
        bindings.push(MaterialBinding {
            region: expect_str(field(&pairs, "region"), "binding.region", out),
            card: expect_str(field(&pairs, "card"), "binding.card", out),
            claim: field(&pairs, "claim").map(|node| expect_str(Some(node), "binding.claim", out)),
            state: expect_str(field(&pairs, "state"), "binding.state", out),
            temp_lo: expect_qty(field(&pairs, "temp-lo"), "binding.temp-lo", out),
            temp_hi: expect_qty(field(&pairs, "temp-hi"), "binding.temp-hi", out),
            source: expect_str(field(&pairs, "source"), "binding.source", out),
        });
    }
    bindings
}

fn read_interface_cards(body: &[Node], out: &mut Vec<Violation>) -> Vec<InterfaceCardBinding> {
    let mut bindings = Vec::new();
    for node in body {
        let Some(("card", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`interface-cards` rows must be `(card ...)`".to_string(),
                fix: "bind TIM/contact cards with `(card :interface ... :card ... :source ...)`"
                    .to_string(),
            });
            continue;
        };
        let pairs = read_pairs(
            inner,
            "card",
            &["interface", "card", "claim", "source"],
            out,
        );
        bindings.push(InterfaceCardBinding {
            interface: expect_str(field(&pairs, "interface"), "card.interface", out),
            card: expect_str(field(&pairs, "card"), "card.card", out),
            claim: field(&pairs, "claim").map(|node| expect_str(Some(node), "card.claim", out)),
            source: expect_str(field(&pairs, "source"), "card.source", out),
        });
    }
    bindings
}

fn read_power(
    body: &[Node],
    out: &mut Vec<Violation>,
    defaults: &mut Vec<DefaultReceipt>,
) -> Vec<PowerDissipation> {
    let mut rows = Vec::new();
    for node in body {
        let Some(("dissipation", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`power` rows must be `(dissipation ...)`".to_string(),
                fix: "declare `(dissipation :region ... :watts ... :duty ...)`".to_string(),
            });
            continue;
        };
        let pairs = read_pairs(inner, "dissipation", &["region", "watts", "duty"], out);
        let region = expect_str(field(&pairs, "region"), "dissipation.region", out);
        let duty = if let Some(node) = field(&pairs, "duty") {
            expect_float(Some(node), "dissipation.duty", out)
        } else {
            defaults.push(DefaultReceipt {
                field: format!("power.dissipation[{region}].duty"),
                value: "1.0".to_string(),
                rationale: "continuous dissipation is the conservative thermal assumption; \
                            any other duty must be stated",
            });
            1.0
        };
        rows.push(PowerDissipation {
            region,
            watts: expect_qty(field(&pairs, "watts"), "dissipation.watts", out),
            duty,
        });
    }
    rows
}

fn read_cooling(body: &[Node], out: &mut Vec<Violation>) -> Option<Cooling> {
    let mut fans = Vec::new();
    let mut vents = Vec::new();
    let mut leakage = None;
    for node in body {
        match section_name(node) {
            Some(("fans", inner)) => {
                for fan_node in inner {
                    let Some(("fan", fan_body)) = section_name(fan_node) else {
                        out.push(Violation {
                            code: "project-malformed-clause",
                            what: "`cooling.fans` rows must be `(fan ...)`".to_string(),
                            fix: "declare `(fan :name ... :flow ... :static-pressure ...)`"
                                .to_string(),
                        });
                        continue;
                    };
                    let pairs =
                        read_pairs(fan_body, "fan", &["name", "flow", "static-pressure"], out);
                    fans.push(Fan {
                        name: expect_str(field(&pairs, "name"), "fan.name", out),
                        flow: expect_qty(field(&pairs, "flow"), "fan.flow", out),
                        static_pressure: expect_qty(
                            field(&pairs, "static-pressure"),
                            "fan.static-pressure",
                            out,
                        ),
                    });
                }
            }
            Some(("vents", inner)) => {
                for vent_node in inner {
                    let Some(("vent", vent_body)) = section_name(vent_node) else {
                        out.push(Violation {
                            code: "project-malformed-clause",
                            what: "`cooling.vents` rows must be `(vent ...)`".to_string(),
                            fix: "declare `(vent :region ... :area ...)`".to_string(),
                        });
                        continue;
                    };
                    let pairs = read_pairs(vent_body, "vent", &["region", "area"], out);
                    vents.push(Vent {
                        region: expect_str(field(&pairs, "region"), "vent.region", out),
                        area: expect_qty(field(&pairs, "area"), "vent.area", out),
                    });
                }
            }
            Some(("leakage", inner)) => {
                let pairs = read_pairs(inner, "leakage", &["watts"], out);
                leakage = Some(expect_qty(field(&pairs, "watts"), "leakage.watts", out));
            }
            _ => {
                out.push(Violation {
                    code: "project-unknown-field",
                    what: "`cooling` carries an unknown subsection".to_string(),
                    fix: "cooling contains exactly `(fans ...)`, `(vents ...)`, `(leakage ...)`"
                        .to_string(),
                });
            }
        }
    }
    let Some(leakage) = leakage else {
        out.push(Violation {
            code: "project-malformed-clause",
            what: "`cooling` lacks its `(leakage :watts ...)` declaration".to_string(),
            fix: "state leakage explicitly, even `0 W`; absence of leakage is a claim".to_string(),
        });
        return None;
    };
    Some(Cooling {
        fans,
        vents,
        leakage,
    })
}

fn read_envelope(body: &[Node], out: &mut Vec<Violation>) -> Envelope {
    let pairs = read_pairs(
        body,
        "envelope",
        &["ambient-lo", "ambient-hi", "pressure"],
        out,
    );
    Envelope {
        ambient_lo: expect_qty(field(&pairs, "ambient-lo"), "envelope.ambient-lo", out),
        ambient_hi: expect_qty(field(&pairs, "ambient-hi"), "envelope.ambient-hi", out),
        pressure: expect_qty(field(&pairs, "pressure"), "envelope.pressure", out),
    }
}

fn read_requirements(body: &[Node], out: &mut Vec<Violation>) -> Vec<ThermalLimit> {
    let mut limits = Vec::new();
    for node in body {
        let Some(("t-limit", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`requirements` rows must be `(t-limit ...)`".to_string(),
                fix: "declare a fully sourced `(t-limit :qoi ... :class ... :region ... \
                      :direction ... :limit ... :margin ... :source-* ... :safety-* ... \
                      :severity ...)`"
                    .to_string(),
            });
            continue;
        };
        let pairs = read_pairs(
            inner,
            "t-limit",
            &[
                "qoi",
                "class",
                "region",
                "direction",
                "limit",
                "margin",
                "source-kind",
                "source-document",
                "source-version",
                "source-locator",
                "safety-factor",
                "safety-source-kind",
                "safety-source-document",
                "safety-source-version",
                "safety-source-locator",
                "severity",
            ],
            out,
        );
        limits.push(ThermalLimit {
            qoi: expect_str(field(&pairs, "qoi"), "t-limit.qoi", out),
            class: expect_str(field(&pairs, "class"), "t-limit.class", out),
            region: expect_str(field(&pairs, "region"), "t-limit.region", out),
            direction: parse_requirement_direction(
                expect_str(field(&pairs, "direction"), "t-limit.direction", out),
                out,
            ),
            limit: expect_qty(field(&pairs, "limit"), "t-limit.limit", out),
            margin: expect_qty(field(&pairs, "margin"), "t-limit.margin", out),
            source: RequirementSource {
                kind: parse_requirement_source_kind(
                    expect_str(field(&pairs, "source-kind"), "t-limit.source-kind", out),
                    "t-limit.source-kind",
                    out,
                ),
                document: expect_str(
                    field(&pairs, "source-document"),
                    "t-limit.source-document",
                    out,
                ),
                version: expect_str(
                    field(&pairs, "source-version"),
                    "t-limit.source-version",
                    out,
                ),
                locator: expect_str(
                    field(&pairs, "source-locator"),
                    "t-limit.source-locator",
                    out,
                ),
            },
            safety_factor: SafetyFactorPolicy {
                factor: expect_float(field(&pairs, "safety-factor"), "t-limit.safety-factor", out),
                source: RequirementSource {
                    kind: parse_requirement_source_kind(
                        expect_str(
                            field(&pairs, "safety-source-kind"),
                            "t-limit.safety-source-kind",
                            out,
                        ),
                        "t-limit.safety-source-kind",
                        out,
                    ),
                    document: expect_str(
                        field(&pairs, "safety-source-document"),
                        "t-limit.safety-source-document",
                        out,
                    ),
                    version: expect_str(
                        field(&pairs, "safety-source-version"),
                        "t-limit.safety-source-version",
                        out,
                    ),
                    locator: expect_str(
                        field(&pairs, "safety-source-locator"),
                        "t-limit.safety-source-locator",
                        out,
                    ),
                },
            },
            severity: parse_requirement_severity(
                expect_str(field(&pairs, "severity"), "t-limit.severity", out),
                out,
            ),
        });
    }
    limits
}

fn read_solver(body: &[Node], out: &mut Vec<Violation>) -> SolverSettings {
    let pairs = read_pairs(body, "solver", &["fidelity", "tolerance-rel"], out);
    SolverSettings {
        fidelity: expect_str(field(&pairs, "fidelity"), "solver.fidelity", out),
        tolerance_rel: expect_float(field(&pairs, "tolerance-rel"), "solver.tolerance-rel", out),
    }
}

fn read_outputs(body: &[Node], out: &mut Vec<Violation>) -> Vec<OutputRequest> {
    let mut outputs = Vec::new();
    for node in body {
        let Some(("qoi", inner)) = section_name(node) else {
            out.push(Violation {
                code: "project-malformed-clause",
                what: "`outputs` rows must be `(qoi ...)`".to_string(),
                fix: "declare `(qoi :name ... :kind scalar|field|report)`".to_string(),
            });
            continue;
        };
        let pairs = read_pairs(inner, "qoi", &["name", "kind"], out);
        outputs.push(OutputRequest {
            name: expect_str(field(&pairs, "name"), "qoi.name", out),
            kind: expect_str(field(&pairs, "kind"), "qoi.kind", out),
        });
    }
    outputs
}
