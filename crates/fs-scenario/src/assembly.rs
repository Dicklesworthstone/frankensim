//! Generic, content-addressed assembly DATA (IRP M6.1-W1a).
//!
//! Route B: the artifact binds canonical scenario IR without changing its
//! bytes. Structural admission checks dimensions, identities and explicit
//! unknowns; module validity, material resolution and rendering are consumers'
//! responsibilities. A source citation does not establish scientific authority.
#![deny(unsafe_code)]

use crate::entity::{EntityCatalog, EntityId, EntityKind, GeometryFingerprint, ToleranceSource};
use crate::{FrameId, Scenario};
use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::{Dims, QtyAny};
use std::collections::BTreeSet;
use std::fmt;

pub mod generators;
pub(crate) mod wire;

/// Frozen `scenario.assembly` wire version; scenario IR remains v2.
pub const ASSEMBLY_DESCRIPTION_SCHEMA_VERSION: u32 = 1;
/// Artifact kind for the content-addressed sidecar.
pub const ASSEMBLY_DESCRIPTION_ARTIFACT_KIND: &str = "scenario-assembly-v1";
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// Generic parameter spelling and its declared six-base SI dimensions.
/// Defined once here for descriptions, gestures and performances.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ParameterKey {
    /// Stable, whitespace-free key, e.g. `open_fraction`.
    pub name: String,
    /// Dimensions expected by every assignment.
    pub dims: Dims,
}
impl ParameterKey {
    /// Construct a dimensioned key. Assembly admission checks its spelling.
    #[must_use]
    pub fn new(name: impl Into<String>, dims: Dims) -> Self {
        Self {
            name: name.into(),
            dims,
        }
    }
}
impl PartialOrd for ParameterKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for ParameterKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (&self.name, self.dims.0).cmp(&(&other.name, other.dims.0))
    }
}
/// A parameter attached to a persistent identity, never a display name.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ParameterPath {
    /// Catalog identity of the parameter owner.
    pub entity: EntityId,
    /// Stable parameter and dimensions.
    pub parameter: ParameterKey,
}
/// Provenance reserved here so anatomy consumers never redefine states.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StateProvenance {
    /// Speech imaging does not imply a playing configuration.
    SpeechImaging {
        ledger_row: String,
        subject: String,
        token: String,
        playing_envelope: Option<ContentHash>,
    },
    /// Measurement of a playing configuration.
    PlayingMeasured { ledger_row: String },
    /// An explicitly declared exploratory state.
    Declared { rationale: String },
}
/// Named assignments shared by all configuration consumers.
#[derive(Debug, Clone, PartialEq)]
pub struct ConfigurationState {
    /// Exact state name; not an instrument-specific Rust type.
    pub name: String,
    /// Ordered assignments; repeated paths refuse.
    pub assignments: Vec<(ParameterPath, QtyAny)>,
    /// Explicit origin and limits of the state.
    pub provenance: StateProvenance,
}
/// Physical reference/load condition of the source geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceState {
    Unloaded,
    Strung,
    Assembled,
    Played,
}
/// Physical part classes, independent of instrument families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartKind {
    String,
    Shell,
    Solid,
    Duct,
    Cavity,
    ValvePattern,
    DrivenBody,
}
/// Imported chart representation; validity is checked by the importing consumer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Representation {
    Mesh,
    Nurbs,
    Sdf,
    Voxel,
}
/// Geometry is generated from data or cites an imported source.
#[derive(Debug, Clone, PartialEq)]
pub enum GeometrySource {
    /// Generator id and dimensioned, ordered parameters.
    Parametric {
        generator: String,
        params: Vec<(ParameterKey, QtyAny)>,
    },
    /// Exact source identity and ledger row; no geometry bytes embedded here.
    Imported {
        fingerprint: GeometryFingerprint,
        representation: Representation,
        ledger_row: String,
    },
}
/// Coordinate-wise reconstruction/material bound, not a joint confidence claim.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionBound {
    pub property: ParameterKey,
    pub lower: QtyAny,
    pub upper: QtyAny,
    pub source: ToleranceSource,
}
/// An unavailable material property cannot carry an invented nominal value.
#[derive(Debug, Clone, PartialEq)]
pub enum UnknownProperty {
    Unknown {
        property: ParameterKey,
        bound: Option<DimensionBound>,
    },
}
/// Material citation plus its physical state, frame and optional layer.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialRef {
    pub card: ContentHash,
    pub query_point: Vec<(ParameterKey, QtyAny)>,
    pub frame: FrameId,
    pub layer: Option<u32>,
    pub unknown: Vec<UnknownProperty>,
}
impl MaterialRef {
    /// Refuse a consumer's attempt to substitute a value for an explicit unknown.
    /// A successful check does not resolve a card or prove that it has a value.
    pub fn require_known(&self, property: &str) -> Result<(), AssemblyError> {
        if self.unknown.iter().any(
            |u| matches!(u, UnknownProperty::Unknown { property: p, .. } if p.name == property),
        ) {
            return Err(AssemblyError::UnknownPropertyDefaulted {
                property: property.to_owned(),
            });
        }
        Ok(())
    }
}
/// One identity-bearing body with geometry and material-state citations.
#[derive(Debug, Clone, PartialEq)]
pub struct PartDecl {
    pub entity: EntityId,
    pub kind: PartKind,
    pub geometry: GeometrySource,
    pub materials: Vec<MaterialRef>,
    pub uncertainty: Vec<DimensionBound>,
}
/// Schema and law content identities, resolved by the coupling consumer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InterfaceKind {
    PowerPort { schema: ContentHash },
    Contact { law: ContentHash },
    Glue { law: ContentHash },
}
/// An interface uses the catalog's ordered or unordered endpoint pairing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InterfaceDecl {
    pub entity: EntityId,
    pub kind: InterfaceKind,
}
/// Dimensioned observation location in an existing scenario frame.
#[derive(Debug, Clone, PartialEq)]
pub struct ObserverDecl {
    pub entity: EntityId,
    pub frame: FrameId,
    pub quantity: ParameterKey,
    pub position: [QtyAny; 3],
}
/// Untrusted authoring input. `None` reference state refuses rather than defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyInput {
    pub root: EntityId,
    pub catalog: EntityCatalog,
    pub parts: Vec<PartDecl>,
    pub interfaces: Vec<InterfaceDecl>,
    pub configurations: Vec<ConfigurationState>,
    pub observers: Vec<ObserverDecl>,
    pub reference_state: Option<ReferenceState>,
}
/// Caller-selected limits applied before traversal/allocation of wire records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssemblyBudget {
    pub max_parts: usize,
    pub max_interfaces: usize,
    pub max_configurations: usize,
    pub max_assignments_per_state: usize,
    pub max_wire_bytes: usize,
    pub max_catalog_entities: usize,
    pub max_catalog_receipts: usize,
    pub max_nested_items: usize,
    pub max_text_bytes: usize,
}
/// Explicit default profile, not an implicit physical value.
pub const DEFAULT_ASSEMBLY_BUDGET: AssemblyBudget = AssemblyBudget {
    max_parts: 4096,
    max_interfaces: 16384,
    max_configurations: 4096,
    max_assignments_per_state: 16384,
    max_wire_bytes: 16 * 1024 * 1024,
    max_catalog_entities: 65536,
    max_catalog_receipts: 262144,
    max_nested_items: 65536,
    max_text_bytes: 4096,
};
/// Structural admission and wire errors. Names retain the offending parameter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssemblyError {
    UnknownPropertyDefaulted {
        property: String,
    },
    MissingReferenceState,
    DimensionMismatch {
        parameter: String,
        expected: Dims,
        found: Dims,
    },
    BudgetExceeded {
        resource: &'static str,
        cap: usize,
    },
    ScenarioHashMismatch {
        expected: ContentHash,
        found: ContentHash,
    },
    Invalid {
        field: &'static str,
        detail: String,
    },
    Entity {
        detail: String,
    },
    Wire {
        offset: usize,
        detail: String,
    },
    UnknownVersion {
        found: u32,
    },
    TrailingBytes,
}
impl fmt::Display for AssemblyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "assembly: {self:?}")
    }
}
impl std::error::Error for AssemblyError {}
impl From<crate::EntityError> for AssemblyError {
    fn from(e: crate::EntityError) -> Self {
        Self::Entity {
            detail: e.to_string(),
        }
    }
}
pub(crate) fn invalid(field: &'static str, detail: impl Into<String>) -> AssemblyError {
    AssemblyError::Invalid {
        field,
        detail: detail.into(),
    }
}
pub(crate) fn cap(n: usize, limit: usize, resource: &'static str) -> Result<(), AssemblyError> {
    if n > limit {
        Err(AssemblyError::BudgetExceeded {
            resource,
            cap: limit,
        })
    } else {
        Ok(())
    }
}
pub(crate) fn text(
    value: &str,
    field: &'static str,
    budget: AssemblyBudget,
) -> Result<(), AssemblyError> {
    cap(value.len(), budget.max_text_bytes, "text bytes")?;
    if value.trim().is_empty() || value.chars().any(char::is_control) {
        return Err(invalid(field, "empty or control-containing text"));
    }
    Ok(())
}
fn key(key: &ParameterKey, budget: AssemblyBudget) -> Result<(), AssemblyError> {
    text(&key.name, "parameter", budget)?;
    if !key
        .name
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-'))
    {
        return Err(invalid("parameter", &key.name));
    }
    Ok(())
}
fn quantity(q: QtyAny, name: &str, dims: Dims) -> Result<(), AssemblyError> {
    if q.dims != dims {
        return Err(AssemblyError::DimensionMismatch {
            parameter: name.to_owned(),
            expected: dims,
            found: q.dims,
        });
    }
    if !q.value.is_finite() {
        return Err(invalid("quantity", name));
    }
    Ok(())
}
fn hash(hash: ContentHash, field: &'static str) -> Result<(), AssemblyError> {
    if hash.0 == [0; 32] {
        return Err(invalid(field, "zero content identity"));
    }
    Ok(())
}
fn bound(value: &DimensionBound, budget: AssemblyBudget) -> Result<(), AssemblyError> {
    key(&value.property, budget)?;
    quantity(value.lower, &value.property.name, value.property.dims)?;
    quantity(value.upper, &value.property.name, value.property.dims)?;
    if value.lower.value > value.upper.value {
        return Err(invalid("bound", &value.property.name));
    }
    for (field, value) in value.source.nonempty_fields() {
        text(value, field, budget)?;
    }
    Ok(())
}
fn parameters(
    values: &[(ParameterKey, QtyAny)],
    budget: AssemblyBudget,
) -> Result<(), AssemblyError> {
    cap(values.len(), budget.max_nested_items, "parameters")?;
    let mut names = BTreeSet::new();
    for (k, q) in values {
        key(k, budget)?;
        quantity(*q, &k.name, k.dims)?;
        if !names.insert(&k.name) {
            return Err(invalid("parameter", format!("duplicate {}", k.name)));
        }
    }
    Ok(())
}
fn entity(
    catalog: &EntityCatalog,
    id: EntityId,
    expected: Option<EntityKind>,
    root: EntityId,
) -> Result<(), AssemblyError> {
    if !catalog.is_active(id) || expected.is_some_and(|kind| id.kind() != kind) {
        return Err(invalid("entity", id.token()));
    }
    let mut current = Some(id);
    for _ in 0..=catalog.entities().len() {
        if current == Some(root) {
            return Ok(());
        }
        current = current
            .and_then(|id| catalog.get(id))
            .and_then(|entry| entry.parent());
        if current.is_none() {
            break;
        }
    }
    Err(invalid("containment", id.token()))
}
fn frame(scenario: &Scenario, frame: FrameId) -> Result<(), AssemblyError> {
    if frame.0 != 0 && !scenario.frames.frames.iter().any(|f| f.id == frame) {
        return Err(invalid("frame", frame.0.to_string()));
    }
    Ok(())
}

pub(crate) fn validate(
    input: &AssemblyInput,
    scenario: &Scenario,
    budget: AssemblyBudget,
) -> Result<(), AssemblyError> {
    cap(input.parts.len(), budget.max_parts, "parts")?;
    cap(input.interfaces.len(), budget.max_interfaces, "interfaces")?;
    cap(
        input.configurations.len(),
        budget.max_configurations,
        "configurations",
    )?;
    cap(input.observers.len(), budget.max_nested_items, "observers")?;
    cap(
        input.catalog.entities().len(),
        budget.max_catalog_entities,
        "catalog entities",
    )?;
    cap(
        input.catalog.receipts().len(),
        budget.max_catalog_receipts,
        "catalog receipts",
    )?;
    for n in [
        input.catalog.datums().len(),
        input.catalog.tolerances().len(),
        input.catalog.placements().len(),
    ] {
        cap(n, budget.max_nested_items, "catalog records")?;
    }
    if input.reference_state.is_none() {
        return Err(AssemblyError::MissingReferenceState);
    }
    entity(
        &input.catalog,
        input.root,
        Some(EntityKind::Assembly),
        input.root,
    )?;
    if !input.catalog.verify_receipts() {
        return Err(invalid("catalog", "invalid receipt chain"));
    }
    let findings = input.catalog.validate();
    if !findings.is_empty() {
        return Err(invalid("catalog", format!("{findings:?}")));
    }
    for p in input.catalog.placements() {
        frame(scenario, p.frame())?;
    }
    let mut seen = BTreeSet::new();
    for part in &input.parts {
        entity(
            &input.catalog,
            part.entity,
            Some(EntityKind::Part),
            input.root,
        )?;
        if !seen.insert(part.entity) {
            return Err(invalid("part", "duplicate identity"));
        }
        match &part.geometry {
            GeometrySource::Parametric { generator, params } => {
                text(generator, "generator", budget)?;
                parameters(params, budget)?;
                generators::validate_source(&part.geometry, budget.max_nested_items)
                    .map_err(|e| invalid("geometry generator", e.to_string()))?;
            }
            GeometrySource::Imported {
                fingerprint,
                ledger_row,
                ..
            } => {
                hash(fingerprint.hash(), "geometry")?;
                text(ledger_row, "ledger row", budget)?;
                if input.catalog.get(part.entity).and_then(|e| e.fingerprint())
                    != Some(*fingerprint)
                {
                    return Err(invalid(
                        "geometry fingerprint",
                        "imported chart does not match its entity declaration",
                    ));
                }
            }
        }
        cap(part.materials.len(), budget.max_nested_items, "materials")?;
        cap(
            part.uncertainty.len(),
            budget.max_nested_items,
            "uncertainty bounds",
        )?;
        for b in &part.uncertainty {
            bound(b, budget)?;
        }
        for material in &part.materials {
            hash(material.card, "material card")?;
            frame(scenario, material.frame)?;
            parameters(&material.query_point, budget)?;
            cap(
                material.unknown.len(),
                budget.max_nested_items,
                "unknown properties",
            )?;
            let mut unknowns = BTreeSet::new();
            for UnknownProperty::Unknown { property, bound: b } in &material.unknown {
                key(property, budget)?;
                if !unknowns.insert(&property.name) {
                    return Err(invalid("unknown property", &property.name));
                }
                if material
                    .query_point
                    .iter()
                    .any(|(k, _)| k.name == property.name)
                {
                    return Err(AssemblyError::UnknownPropertyDefaulted {
                        property: property.name.clone(),
                    });
                }
                if let Some(b) = b {
                    bound(b, budget)?;
                    if b.property != *property {
                        return Err(invalid("unknown property bound", &property.name));
                    }
                }
            }
        }
    }
    seen.clear();
    for interface in &input.interfaces {
        entity(
            &input.catalog,
            interface.entity,
            Some(EntityKind::Interface),
            input.root,
        )?;
        if !seen.insert(interface.entity) {
            return Err(invalid("interface", "duplicate identity"));
        }
        let pair = input
            .catalog
            .get(interface.entity)
            .and_then(|e| e.pair())
            .ok_or_else(|| invalid("interface", "missing endpoints"))?;
        for endpoint in [pair.from(), pair.to()] {
            entity(&input.catalog, endpoint, None, input.root)?;
        }
        let identity = match interface.kind {
            InterfaceKind::PowerPort { schema } => schema,
            InterfaceKind::Contact { law } | InterfaceKind::Glue { law } => law,
        };
        hash(identity, "interface schema/law")?;
    }
    let mut states = BTreeSet::new();
    for state in &input.configurations {
        text(&state.name, "configuration", budget)?;
        if !states.insert(&state.name) {
            return Err(invalid("configuration", "duplicate name"));
        }
        cap(
            state.assignments.len(),
            budget.max_assignments_per_state,
            "assignments per state",
        )?;
        let mut paths = BTreeSet::new();
        for (path, value) in &state.assignments {
            entity(&input.catalog, path.entity, None, input.root)?;
            key(&path.parameter, budget)?;
            quantity(*value, &path.parameter.name, path.parameter.dims)?;
            if !paths.insert((path.entity, &path.parameter.name)) {
                return Err(invalid("assignment", "duplicate parameter path"));
            }
        }
        match &state.provenance {
            StateProvenance::Declared { rationale } => text(rationale, "state rationale", budget)?,
            StateProvenance::PlayingMeasured { ledger_row } => {
                text(ledger_row, "state ledger row", budget)?
            }
            StateProvenance::SpeechImaging {
                ledger_row,
                subject,
                token,
                playing_envelope,
            } => {
                text(ledger_row, "state ledger row", budget)?;
                text(subject, "subject", budget)?;
                text(token, "speech token", budget)?;
                if let Some(h) = playing_envelope {
                    hash(*h, "playing envelope")?;
                }
            }
        }
    }
    for observer in &input.observers {
        entity(&input.catalog, observer.entity, None, input.root)?;
        frame(scenario, observer.frame)?;
        key(&observer.quantity, budget)?;
        for p in observer.position {
            quantity(p, "observer position", LENGTH)?;
        }
    }
    Ok(())
}

/// Structurally admitted artifact; no material/kernel validity is inferred.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyDescription {
    scenario_ir_hash: ContentHash,
    input: AssemblyInput,
}

/// Applied rules for the additive sidecar route; no tolerance records were
/// present in the old positional scenario grammar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyAttachmentRule {
    /// The existing v1-to-v2 dimension crosswalk ran.
    AppendMoleZero,
    /// A noncanonical v2 spelling was re-emitted by the existing IR reader.
    CanonicalizeScenarioV2,
    /// Canonical scenario values/bytes are retained unchanged by attachment.
    PreserveScenarioV2,
    /// The new frozen assembly v1 artifact was attached by its content hash.
    AttachAssemblyV1,
}
/// Parser-issued binding of exact prior source, unchanged canonical scenario
/// and the new assembly artifact. This is migration/attachment evidence, not
/// measurement validation or an old tolerance-IR migration claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyAttachmentReceipt {
    source_hash: ContentHash,
    canonical_hash: ContentHash,
    assembly_hash: ContentHash,
    rules: Vec<AssemblyAttachmentRule>,
    dimension_crosswalk: Option<crate::ir::DimensionCrosswalkReceipt>,
    source_canonicalization: Option<crate::ir::SourceCanonicalizationReceipt>,
}
impl AssemblyAttachmentReceipt {
    /// Ordered applied rules, including preservation and attachment.
    #[must_use]
    pub fn rules(&self) -> &[AssemblyAttachmentRule] {
        &self.rules
    }
    /// Existing five-to-six receipt when the source was legacy v1.
    #[must_use]
    pub const fn dimension_crosswalk(&self) -> Option<&crate::ir::DimensionCrosswalkReceipt> {
        self.dimension_crosswalk.as_ref()
    }
    /// Exact prior bytes and artifact identity must match the retained receipt.
    #[must_use]
    pub fn verifies(&self, source: &[u8], canonical: &[u8], assembly: &[u8]) -> bool {
        self.source_hash == hash_bytes(source)
            && self.canonical_hash == hash_bytes(canonical)
            && self.assembly_hash == hash_bytes(assembly)
            && self
                .dimension_crosswalk
                .as_ref()
                .is_none_or(|r| r.verifies(source, canonical))
            && self
                .source_canonicalization
                .as_ref()
                .is_none_or(|r| r.verifies(source, canonical))
    }
}
impl AssemblyDescription {
    /// Admit any supported prior scenario through its existing reader, attach
    /// the sidecar and retain the rules and exact source/canonical identities.
    pub fn attach_to_ir(
        source: &str,
        input: AssemblyInput,
        budget: AssemblyBudget,
    ) -> Result<(Self, AssemblyAttachmentReceipt), AssemblyError> {
        let decoded =
            crate::ir::parse_ir(source).map_err(|e| invalid("scenario IR", e.to_string()))?;
        let description = Self::new(decoded.scenario(), input, budget)?;
        let mut rules = Vec::new();
        if decoded.migration().is_some() {
            rules.push(AssemblyAttachmentRule::AppendMoleZero);
        }
        if decoded.canonicalization().is_some() {
            rules.push(AssemblyAttachmentRule::CanonicalizeScenarioV2);
        }
        rules.extend([
            AssemblyAttachmentRule::PreserveScenarioV2,
            AssemblyAttachmentRule::AttachAssemblyV1,
        ]);
        let receipt = AssemblyAttachmentReceipt {
            source_hash: hash_bytes(source.as_bytes()),
            canonical_hash: description.scenario_ir_hash,
            assembly_hash: description.content_hash(budget)?,
            rules,
            dimension_crosswalk: decoded.migration().cloned(),
            source_canonicalization: decoded.canonicalization().cloned(),
        };
        Ok((description, receipt))
    }
    /// Admit explicit assembly data and bind its scenario's canonical bytes.
    pub fn new(
        scenario: &Scenario,
        input: AssemblyInput,
        budget: AssemblyBudget,
    ) -> Result<Self, AssemblyError> {
        validate(&input, scenario, budget)?;
        let result = Self {
            scenario_ir_hash: hash_bytes(crate::ir::write_ir(scenario).as_bytes()),
            input,
        };
        result.to_bytes(budget)?;
        Ok(result)
    }
    /// Immutable typed data for validity, coupling and performance consumers.
    #[must_use]
    pub const fn input(&self) -> &AssemblyInput {
        &self.input
    }
    /// Exact canonical scenario binding.
    #[must_use]
    pub const fn scenario_ir_hash(&self) -> ContentHash {
        self.scenario_ir_hash
    }
    /// Refuse a sidecar attached to a different scenario.
    pub fn verify_scenario(&self, scenario: &Scenario) -> Result<(), AssemblyError> {
        let found = hash_bytes(crate::ir::write_ir(scenario).as_bytes());
        if found != self.scenario_ir_hash {
            return Err(AssemblyError::ScenarioHashMismatch {
                expected: self.scenario_ir_hash,
                found,
            });
        }
        Ok(())
    }
    /// Canonical bounded bytes: finite shortest-round-trip quantities, six SI exponents.
    pub fn to_bytes(&self, budget: AssemblyBudget) -> Result<Vec<u8>, AssemblyError> {
        wire::encode(self, budget)
    }
    /// Decode, re-admit and verify an exact scenario-bound artifact.
    pub fn from_bytes(
        bytes: &[u8],
        scenario: &Scenario,
        budget: AssemblyBudget,
    ) -> Result<Self, AssemblyError> {
        wire::decode(bytes, scenario, budget)
    }
    /// Content identity of the exact frozen artifact bytes.
    pub fn content_hash(&self, budget: AssemblyBudget) -> Result<ContentHash, AssemblyError> {
        Ok(hash_bytes(&self.to_bytes(budget)?))
    }
}
