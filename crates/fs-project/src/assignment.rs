//! L6 project adapter for fs-io's deterministic mesh assignments.
//!
//! The `.fsim` document names geometry artifacts and declares selectors by
//! project-local entity name. This module resolves those names through
//! [`ProjectSpec::resolve_entities`], supplies the resulting persistent
//! [`EntityId`] tokens to fs-io, and retains fs-io's complete lower-layer
//! report without changing its meaning. Imported meshes remain
//! caller-supplied: this adapter proves exact project/report binding, not who
//! supplied the geometry or whether it represents the physical part.

use core::fmt::Write as _;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{DomainHasher, hash_domain};
use fs_conduction::{ConductionMesh, InterfaceFacePair, ThermalInterfaces};
use fs_exec::Cx;
use fs_io::{
    AssignmentLimits, AssignmentReport, AssignmentRequest, NamedFaceGroup, resolve_mesh_assignments,
};
use fs_rep_mesh::{Soup, point_triangle_distance};
use fs_scenario::{EntityId, EntityKind, Violation};

use crate::spec::{EntityDecl, GeometryArtifact, GeometryAssignment, ProjectSpec};

/// Domain for the identity hook derived from one exact project geometry row.
pub const GEOMETRY_SOURCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-project.geometry-source.v1";

/// Domain for retained exact fs-io assignment-report JSON bytes.
pub const GEOMETRY_ASSIGNMENT_REPORT_DOMAIN: &str =
    "org.frankensim.fs-project.geometry-assignment-report.v1";

const PROJECT_ASSIGNMENT_NO_CLAIM: &str = "the adapter binds exact project entity identities, one declared selector plan, and one supplied finite tessellation; it does not authenticate the mesh supplier, prove continuum/CAD/physical-region sameness, or make fs-io face ordinals stable across re-import";
const INTERFACE_AUDIT_NO_CLAIM: &str = "the audit reports finite-mesh region pairs whose supplied triangle soups approach within the declared tolerance in one shared coordinate unit; it does not certify continuum contact, infer a physical interface law, authenticate assembly transforms, or prove that a declared interface is complete";
const CONDUCTION_INTERFACE_NO_CLAIM: &str = "the lowering binds exact coordinate-bit equality among one resolved imported triangle soup, its retained face ordinals, declared region/interface selectors, and one exact ConductionMesh boundary; it does not authenticate the importer, prove continuum or CAD identity, infer region ownership when oriented source faces do not establish it, select an interface card, or construct a thermal operator";

/// One caller-supplied promoted mesh and its importer-provided named groups.
#[derive(Debug)]
struct ImportedMesh {
    soup: Soup,
    length_unit: String,
    named_groups: Vec<NamedFaceGroup>,
}

/// Caller-supplied mesh store keyed by the exact project geometry-row
/// identity. The store computes its own keys, so a key cannot lie about the
/// [`GeometryArtifact`] it was inserted for. Authenticity of the supplied
/// mesh remains the caller's trust channel.
#[derive(Debug, Default)]
pub struct ImportedMeshLibrary {
    meshes: BTreeMap<String, ImportedMesh>,
}

impl ImportedMeshLibrary {
    /// An empty imported-mesh library.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a promoted mesh for one exact project geometry row.
    ///
    /// Returns the source identity fs-io will retain in its receipt.
    pub fn insert(
        &mut self,
        artifact: &GeometryArtifact,
        soup: Soup,
        length_unit: impl Into<String>,
        named_groups: Vec<NamedFaceGroup>,
    ) -> String {
        let source_identity = geometry_source_identity(artifact);
        self.meshes.insert(
            source_identity.clone(),
            ImportedMesh {
                soup,
                length_unit: length_unit.into(),
                named_groups,
            },
        );
        source_identity
    }

    fn get(&self, source_identity: &str) -> Option<&ImportedMesh> {
        self.meshes.get(source_identity)
    }
}

/// The persistent entity corresponding to one row of a retained fs-io
/// assignment report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedProjectAssignment {
    /// Project-local declared name.
    pub declared_target: String,
    /// Recomputed persistent identity.
    pub entity_id: EntityId,
}

/// One geometry artifact's exact lower-layer report plus its L6 entity
/// binding and retention material.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedGeometryArtifact {
    /// Project geometry role.
    pub artifact_role: String,
    /// Exact identity hook supplied to fs-io.
    pub source_identity: String,
    /// fs-io's complete success report, unchanged.
    pub report: AssignmentReport,
    /// Exact canonical one-line JSON returned by `AssignmentReport::to_json`.
    pub report_bytes: Vec<u8>,
    /// Domain-separated hash of `report_bytes`.
    pub report_hash: String,
    /// Entity bindings in exactly the same order as `report.assignments`.
    pub entities: Vec<ResolvedProjectAssignment>,
}

/// Atomic result of resolving every declared geometry assignment.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GeometryResolution {
    /// Typed refusals. Empty means every declared assignment resolved.
    pub violations: Vec<Violation>,
    /// Success rows in project geometry declaration order. Empty on any
    /// refusal: partial assignment publication is forbidden.
    pub artifacts: Vec<ResolvedGeometryArtifact>,
}

impl GeometryResolution {
    /// True when every declared assignment resolved.
    #[must_use]
    pub fn admissible(&self) -> bool {
        self.violations.is_empty()
    }

    /// Exact retained assignment reports, ready for ledger storage.
    pub fn receipts(&self) -> impl Iterator<Item = &ResolvedGeometryArtifact> {
        self.artifacts.iter()
    }

    /// Deterministic assignment table for CLI logs and end-to-end evidence.
    #[must_use]
    pub fn render_table(&self) -> String {
        let mut output = String::new();
        for artifact in &self.artifacts {
            for (entity, assignment) in artifact.entities.iter().zip(&artifact.report.assignments) {
                let volume = assignment
                    .stats
                    .enclosed_volume
                    .map_or_else(|| "-".to_string(), |value| value.to_string());
                let _ = writeln!(
                    output,
                    "{} | entity {} | artifact {} | source {} | unit {} | selector {:016x} | faces {} | area {} | volume {} | bounds [{},{},{}]..[{},{},{}] | report {}",
                    entity.declared_target,
                    entity.entity_id,
                    artifact.artifact_role,
                    artifact.source_identity,
                    artifact.report.receipt.length_unit(),
                    assignment.selector_fingerprint,
                    assignment.stats.face_count,
                    assignment.stats.surface_area,
                    volume,
                    assignment.stats.bounds_min[0],
                    assignment.stats.bounds_min[1],
                    assignment.stats.bounds_min[2],
                    assignment.stats.bounds_max[0],
                    assignment.stats.bounds_max[1],
                    assignment.stats.bounds_max[2],
                    artifact.report_hash,
                );
            }
        }
        output
    }

    /// Explicit no-claim boundary shared by every successful row.
    #[must_use]
    pub const fn no_claim() -> &'static str {
        PROJECT_ASSIGNMENT_NO_CLAIM
    }
}

/// Explicit resource and geometric envelope for undeclared-contact detection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceAuditLimits {
    /// Inclusive separation threshold in the common coordinate unit retained
    /// by the compared assignment reports.
    pub proximity_tolerance: f64,
    /// Maximum exact triangle-pair tests across the complete audit.
    pub max_triangle_pair_tests: u64,
}

impl InterfaceAuditLimits {
    /// Conservative default for interactive project validation.
    pub const DEFAULT: Self = Self {
        proximity_tolerance: 1.0e-6,
        max_triangle_pair_tests: 1_000_000,
    };
}

impl Default for InterfaceAuditLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// One finite-mesh proximity that lacks a declared interface entity.
#[derive(Debug, Clone, PartialEq)]
pub struct UndeclaredContactCandidate {
    /// Lexicographically first region name.
    pub first_region: String,
    /// Lexicographically second region name.
    pub second_region: String,
    /// Geometry role supplying `first_region`.
    pub first_artifact: String,
    /// Geometry role supplying `second_region`.
    pub second_artifact: String,
    /// Closest retained face ordinal on `first_region`.
    pub first_face: u32,
    /// Closest retained face ordinal on `second_region`.
    pub second_face: u32,
    /// Computed finite-triangle separation in `length_unit`.
    pub separation: f64,
    /// Common coordinate unit of both supplied meshes.
    pub length_unit: String,
}

/// Fail-closed result of checking resolved region geometry against declared
/// interface pairs.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct InterfaceDeclarationAudit {
    /// Typed precondition, resource, or cancellation refusals.
    pub violations: Vec<Violation>,
    /// Undeclared region pairs at or below the admitted proximity tolerance.
    pub undeclared_contacts: Vec<UndeclaredContactCandidate>,
    /// Exact triangle-pair distance tests consumed.
    pub triangle_pair_tests: u64,
}

impl InterfaceDeclarationAudit {
    /// True only when the audit completed and found no undeclared contact.
    #[must_use]
    pub fn admissible(&self) -> bool {
        self.violations.is_empty() && self.undeclared_contacts.is_empty()
    }

    /// Explicit finite-mesh and coordinate-frame no-claim boundary.
    #[must_use]
    pub const fn no_claim() -> &'static str {
        INTERFACE_AUDIT_NO_CLAIM
    }
}

/// Explicit resource envelope for exact imported-face to conduction-boundary
/// lowering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConductionInterfaceLimits {
    /// Maximum selected source faces indexed across the complete lowering.
    pub max_source_faces: u64,
}

impl ConductionInterfaceLimits {
    /// Conservative default for one project-to-conduction adapter call.
    pub const DEFAULT: Self = Self {
        max_source_faces: 1_000_000,
    };
}

impl Default for ConductionInterfaceLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Exact imported face retained as provenance for one conduction boundary
/// slot.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ConductionSourceFace {
    /// Project geometry role whose immutable assignment report selected it.
    pub artifact: String,
    /// Exact source identity used by the assignment receipt.
    pub source_identity: String,
    /// Face ordinal in that exact imported triangle soup.
    pub face: u32,
}

/// One declared scenario interface lowered to an oriented pair of exact
/// `ConductionMesh::boundary()` slots.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedConductionInterfacePair {
    /// Declared interface entity name.
    pub interface: String,
    /// Declared `from` region name.
    pub from_region: String,
    /// Declared `to` region name.
    pub to_region: String,
    /// Boundary slot whose outward-oriented source face belongs to `from`.
    pub from_boundary_slot: usize,
    /// Boundary slot whose outward-oriented source face belongs to `to`.
    pub to_boundary_slot: usize,
    /// Imported source face proving the `from` slot association.
    pub from_source: ConductionSourceFace,
    /// Imported source face proving the `to` slot association.
    pub to_source: ConductionSourceFace,
    /// Imported faces selected by the interface entity for this exact
    /// geometric triangle. One shared CAD face or both duplicated traces are
    /// admitted; an empty selector match refuses.
    pub interface_sources: Vec<ConductionSourceFace>,
}

impl ResolvedConductionInterfacePair {
    /// The fs-conduction face-pair value, preserving declared scenario
    /// orientation.
    #[must_use]
    pub const fn face_pair(&self) -> InterfaceFacePair {
        InterfaceFacePair {
            side_a: self.from_boundary_slot,
            side_b: self.to_boundary_slot,
        }
    }
}

/// Atomic result of lowering every coincident conduction trace to declared
/// scenario interface semantics.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ConductionInterfaceResolution {
    /// Typed precondition, geometry, resource, or cancellation refusals.
    pub violations: Vec<Violation>,
    /// Lowered pairs in deterministic conduction candidate order. Empty on
    /// any refusal; partial interface publication is forbidden.
    pub pairs: Vec<ResolvedConductionInterfacePair>,
    /// Exact selected source faces indexed before publication.
    pub source_faces_indexed: u64,
}

impl ConductionInterfaceResolution {
    /// True only when every coincident conduction trace has exactly one
    /// declared interface and every declared interface lowered at least once.
    #[must_use]
    pub fn admissible(&self) -> bool {
        self.violations.is_empty()
    }

    /// Deterministic table for end-to-end evidence and operator handoff.
    #[must_use]
    pub fn render_table(&self) -> String {
        let mut output = String::new();
        for pair in &self.pairs {
            let interface_faces = pair
                .interface_sources
                .iter()
                .map(|source| {
                    format!(
                        "{}:{}:{}",
                        source.artifact, source.source_identity, source.face
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            let _ = writeln!(
                output,
                "{} | {} slot {} <- {}:{}:{} | {} slot {} <- {}:{}:{} | interface-faces [{}]",
                pair.interface,
                pair.from_region,
                pair.from_boundary_slot,
                pair.from_source.artifact,
                pair.from_source.source_identity,
                pair.from_source.face,
                pair.to_region,
                pair.to_boundary_slot,
                pair.to_source.artifact,
                pair.to_source.source_identity,
                pair.to_source.face,
                interface_faces,
            );
        }
        output
    }

    /// Explicit authority boundary shared by every successful row.
    #[must_use]
    pub const fn no_claim() -> &'static str {
        CONDUCTION_INTERFACE_NO_CLAIM
    }
}

/// Derive the opaque fs-io source-identity hook from one exact project
/// geometry row. Every variable-width field is length-framed.
#[must_use]
pub fn geometry_source_identity(artifact: &GeometryArtifact) -> String {
    let mut hasher = DomainHasher::new(GEOMETRY_SOURCE_IDENTITY_DOMAIN);
    absorb_bytes(&mut hasher, artifact.role.as_bytes());
    absorb_bytes(&mut hasher, artifact.format.as_bytes());
    hasher.update(&artifact.source_hash.to_le_bytes());
    absorb_bytes(&mut hasher, artifact.parser_version.as_bytes());
    format!("geometry:{}", hasher.finalize().to_hex())
}

fn absorb_bytes(hasher: &mut DomainHasher, bytes: &[u8]) {
    let length = bytes.len() as u128;
    hasher.update(&length.to_le_bytes());
    hasher.update(bytes);
}

fn violation(code: &'static str, what: impl Into<String>, fix: impl Into<String>) -> Violation {
    Violation {
        code,
        what: what.into(),
        fix: fix.into(),
    }
}

/// Resolve every persisted selector plan against the exact imported meshes.
///
/// Entity resolution and plan preflight complete before any geometric work.
/// Reports publish atomically only after every artifact succeeds.
#[must_use]
pub fn resolve_geometry_assignments(
    spec: &ProjectSpec,
    library: &ImportedMeshLibrary,
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> GeometryResolution {
    let mut result = GeometryResolution::default();
    if cx.checkpoint().is_err() {
        result
            .violations
            .push(cancelled("project-assignment-entry"));
        return result;
    }

    let (Some(geometry), Some(assignments), Some(_assembly)) =
        (&spec.geometry, &spec.assignments, &spec.assembly)
    else {
        result.violations.push(violation(
            "project-assignment-preconditions",
            "geometry assignment resolution needs the `geometry`, `assignments`, and `assembly` sections",
            "declare the missing sections and run structural project validation first",
        ));
        return result;
    };

    let ids = spec.resolve_entities(&mut result.violations);
    if !result.violations.is_empty() {
        return result;
    }

    let mut geometry_by_role = BTreeMap::new();
    for artifact in geometry {
        if geometry_by_role
            .insert(artifact.role.as_str(), artifact)
            .is_some()
        {
            result.violations.push(violation(
                "project-geometry-role-duplicate",
                format!(
                    "geometry role `{}` is declared more than once",
                    artifact.role
                ),
                "give every imported geometry artifact one unique project role",
            ));
        }
    }

    let mut expected_targets = BTreeSet::new();
    for (name, id) in &ids {
        if matches!(id.kind(), EntityKind::Region | EntityKind::Interface) {
            expected_targets.insert(name.as_str());
        }
    }

    let mut rows_by_artifact: BTreeMap<&str, Vec<&GeometryAssignment>> = BTreeMap::new();
    let mut assigned_targets = BTreeSet::new();
    for assignment in assignments {
        let Some(id) = ids.get(&assignment.target) else {
            result.violations.push(violation(
                "project-assignment-target-unknown",
                format!(
                    "assignment on `{}` targets undeclared entity `{}`",
                    assignment.artifact, assignment.target
                ),
                "reference a declared region or interface by its exact project-local name",
            ));
            continue;
        };
        if !matches!(id.kind(), EntityKind::Region | EntityKind::Interface) {
            result.violations.push(violation(
                "project-assignment-target-kind",
                format!(
                    "assignment target `{}` is a {}, not a region or interface",
                    assignment.target,
                    id.kind().label()
                ),
                "assign finite mesh faces only to declared regions or interfaces",
            ));
            continue;
        }
        if !geometry_by_role.contains_key(assignment.artifact.as_str()) {
            result.violations.push(violation(
                "project-assignment-artifact-unknown",
                format!(
                    "assignment for `{}` references undeclared geometry role `{}`",
                    assignment.target, assignment.artifact
                ),
                "reference the exact role of one declared geometry artifact",
            ));
            continue;
        }
        if !assigned_targets.insert(assignment.target.as_str()) {
            result.violations.push(violation(
                "project-assignment-target-duplicate",
                format!(
                    "entity `{}` has more than one geometry assignment",
                    assignment.target
                ),
                "declare exactly one artifact and selector for each region or interface",
            ));
            continue;
        }
        rows_by_artifact
            .entry(assignment.artifact.as_str())
            .or_default()
            .push(assignment);
    }

    for target in expected_targets {
        if !assigned_targets.contains(target) {
            result.violations.push(violation(
                "project-assignment-target-unbound",
                format!("region/interface `{target}` has no geometry assignment"),
                "declare one mesh-index-free named or geometric selector for this entity",
            ));
        }
    }
    for artifact in geometry {
        if !rows_by_artifact.contains_key(artifact.role.as_str()) {
            result.violations.push(violation(
                "project-assignment-artifact-unassigned",
                format!(
                    "geometry role `{}` has no region/interface assignments",
                    artifact.role
                ),
                "declare at least one selector for every imported geometry artifact, or remove the unused artifact",
            ));
        }
    }
    if assignments.len() > limits.max_requests {
        result.violations.push(violation(
            "mesh-assignment-resource-bound",
            format!(
                "project declares {} assignment requests, exceeding the admitted total {}",
                assignments.len(),
                limits.max_requests
            ),
            "raise the explicit assignment limit within the run budget or reduce the selector plan",
        ));
    }
    if !result.violations.is_empty() {
        return result;
    }

    let mut pending = Vec::new();
    let mut total_selected = 0usize;
    for (artifact_index, artifact) in geometry.iter().enumerate() {
        if cx.checkpoint().is_err() {
            result
                .violations
                .push(cancelled("project-assignment-artifact"));
            return result;
        }
        let source_identity = geometry_source_identity(artifact);
        let Some(imported) = library.get(&source_identity) else {
            result.violations.push(violation(
                "project-assignment-mesh-unavailable",
                format!(
                    "geometry role `{}` has no supplied promoted mesh for source identity `{source_identity}`",
                    artifact.role
                ),
                "run quarantine/import/promotion for this exact geometry receipt and insert the resulting mesh under the project geometry row",
            ));
            return result;
        };
        let Some(rows) = rows_by_artifact.get(artifact.role.as_str()) else {
            result.violations.push(violation(
                "project-assignment-report-mismatch",
                format!(
                    "geometry role `{}` lost its assignment plan after successful preflight",
                    artifact.role
                ),
                "treat this as an internal adapter defect; do not retain or publish a report",
            ));
            return result;
        };
        if rows
            .iter()
            .any(|row| row.length_unit != imported.length_unit)
        {
            result.violations.push(violation(
                "project-assignment-unit-mismatch",
                format!(
                    "geometry role `{}` was supplied in unit `{}` but one or more selector rows declare another length unit",
                    artifact.role, imported.length_unit
                ),
                "make every selector row state the exact coordinate unit carried by the promoted mesh",
            ));
            return result;
        }

        let mut requests = Vec::with_capacity(rows.len());
        let mut entities = Vec::with_capacity(rows.len());
        for row in rows {
            let Some(entity_id) = ids.get(&row.target).copied() else {
                result.violations.push(violation(
                    "project-assignment-report-mismatch",
                    format!(
                        "geometry role `{}` lost target `{}` after successful preflight",
                        artifact.role, row.target
                    ),
                    "treat this as an internal adapter defect; do not retain or publish a report",
                ));
                return result;
            };
            requests.push(AssignmentRequest {
                subject: entity_id.token(),
                selector: row.selector.clone(),
                allow_overlap: row.allow_overlap,
            });
            entities.push(ResolvedProjectAssignment {
                declared_target: row.target.clone(),
                entity_id,
            });
        }
        let report = match resolve_mesh_assignments(
            &imported.soup,
            &source_identity,
            &imported.length_unit,
            &imported.named_groups,
            &requests,
            limits,
            cx,
        ) {
            Ok(report) => report,
            Err(refusal) => {
                result.violations.push(violation(
                    refusal.code,
                    format!("geometry role `{}`: {}", artifact.role, refusal.what),
                    refusal.fix,
                ));
                return result;
            }
        };

        let selected = report
            .assignments
            .iter()
            .try_fold(0usize, |sum, assignment| {
                sum.checked_add(assignment.faces.len())
            });
        let Some(selected) = selected.and_then(|value| total_selected.checked_add(value)) else {
            result.violations.push(violation(
                "mesh-assignment-resource-bound",
                "aggregate selected-face count overflowed the platform range",
                "reduce the imported mesh or assignment plan",
            ));
            return result;
        };
        total_selected = selected;
        if total_selected > limits.max_selected_faces {
            result.violations.push(violation(
                "mesh-assignment-resource-bound",
                format!(
                    "project-wide selected-face count {total_selected} exceeds admitted total {}",
                    limits.max_selected_faces
                ),
                "raise the explicit assignment limit within the run budget or reduce the imported mesh",
            ));
            return result;
        }

        if entities.len() != report.assignments.len()
            || entities
                .iter()
                .zip(&report.assignments)
                .any(|(entity, assignment)| entity.entity_id.token() != assignment.subject)
        {
            result.violations.push(violation(
                "project-assignment-report-mismatch",
                format!(
                    "geometry role `{}` returned a report whose subject order differs from the compiled persistent-identity plan",
                    artifact.role
                ),
                "treat this as an internal adapter defect; do not retain or publish the report",
            ));
            return result;
        }
        let report_bytes = report.to_json().into_bytes();
        let report_hash = hash_domain(GEOMETRY_ASSIGNMENT_REPORT_DOMAIN, &report_bytes).to_hex();
        pending.push(ResolvedGeometryArtifact {
            artifact_role: artifact.role.clone(),
            source_identity,
            report,
            report_bytes,
            report_hash,
            entities,
        });

        if artifact_index + 1 < geometry.len() && cx.checkpoint().is_err() {
            result
                .violations
                .push(cancelled("project-assignment-between-artifacts"));
            return result;
        }
    }

    result.artifacts = pending;
    result
}

type ConductionPointKey = [u64; 3];
type ConductionFaceKey = [ConductionPointKey; 3];

#[derive(Debug, Clone)]
struct IndexedConductionSourceFace {
    source: ConductionSourceFace,
    outward: [f64; 3],
}

/// Lower exact, oriented imported surface assignments to the boundary-slot
/// identity consumed by fs-conduction.
///
/// This function deliberately resolves geometry assignments itself while
/// borrowing the imported library and conduction mesh for the complete call.
/// A source face associates with a boundary slot only when all three vertex
/// coordinate bit patterns agree and the source triangle orientation agrees
/// with that slot's outward normal. Interface selectors must cover the same
/// exact triangle, every coincident conduction pair must be claimed once, and
/// every declared interface must lower at least once.
#[must_use]
pub fn resolve_conduction_interface_pairs(
    spec: &ProjectSpec,
    library: &ImportedMeshLibrary,
    assignment_limits: AssignmentLimits,
    limits: ConductionInterfaceLimits,
    mesh: &ConductionMesh,
    cx: &Cx<'_>,
) -> ConductionInterfaceResolution {
    let mut result = ConductionInterfaceResolution::default();
    if cx.checkpoint().is_err() {
        result
            .violations
            .push(conduction_interface_cancelled("conduction-interface-entry"));
        return result;
    }

    let geometry = resolve_geometry_assignments(spec, library, assignment_limits, cx);
    if !geometry.admissible() {
        result.violations = geometry.violations;
        return result;
    }
    let Some(assembly) = &spec.assembly else {
        result.violations.push(violation(
            "project-conduction-interface-preconditions",
            "conduction interface lowering requires the `assembly` section",
            "declare the region and interface entities before lowering boundary slots",
        ));
        return result;
    };

    let mut declarations = Vec::new();
    let mut declaration_pairs = BTreeMap::<(&str, &str), &str>::new();
    for entity in assembly {
        if let EntityDecl::Interface { name, from, to, .. } = entity {
            let region_pair = ordered_pair(from, to);
            if let Some(previous) = declaration_pairs.insert(region_pair, name) {
                result.violations.push(violation(
                    "project-conduction-interface-pair-ambiguous",
                    format!(
                        "interfaces `{previous}` and `{name}` both join regions `{}` and `{}`",
                        region_pair.0, region_pair.1
                    ),
                    "declare one thermal interface entity per region pair; staged mechanical joints need a separate typed relationship",
                ));
            }
            declarations.push((name.as_str(), from.as_str(), to.as_str()));
        }
    }
    if !result.violations.is_empty() {
        return result;
    }

    let mut index =
        BTreeMap::<(String, ConductionFaceKey), Vec<IndexedConductionSourceFace>>::new();
    for artifact in &geometry.artifacts {
        let Some(imported) = library.get(&artifact.source_identity) else {
            result.violations.push(violation(
                "project-conduction-interface-mesh-unavailable",
                format!(
                    "resolved geometry role `{}` no longer has its exact imported mesh",
                    artifact.artifact_role
                ),
                "keep the imported mesh library alive through conduction boundary lowering",
            ));
            return result;
        };
        if artifact.entities.len() != artifact.report.assignments.len() {
            result.violations.push(violation(
                "project-conduction-interface-report-mismatch",
                format!(
                    "resolved geometry role `{}` has different entity and assignment row counts",
                    artifact.artifact_role
                ),
                "discard the inconsistent resolution and rerun the project assignment adapter",
            ));
            return result;
        }
        if artifact.report.receipt.length_unit() != "m" {
            result.violations.push(violation(
                "project-conduction-interface-unit",
                format!(
                    "geometry role `{}` retains length unit `{}`, but ConductionMesh coordinates are SI metres",
                    artifact.artifact_role,
                    artifact.report.receipt.length_unit()
                ),
                "promote and assign the geometry in explicit `m` coordinates before exact conduction lowering",
            ));
            return result;
        }
        for (entity, assignment) in artifact.entities.iter().zip(&artifact.report.assignments) {
            for &face in &assignment.faces {
                if result.source_faces_indexed == limits.max_source_faces {
                    result.violations.push(violation(
                        "project-conduction-interface-resource-bound",
                        format!(
                            "conduction interface lowering reached the admitted {} selected source faces",
                            limits.max_source_faces
                        ),
                        "raise the explicit source-face limit within the run budget or reduce the selector plan",
                    ));
                    return result;
                }
                if result.source_faces_indexed % 1_024 == 0 && cx.checkpoint().is_err() {
                    result.violations.push(conduction_interface_cancelled(
                        "conduction-interface-source-index",
                    ));
                    return result;
                }
                result.source_faces_indexed += 1;
                let face_index = face as usize;
                let triangle = imported.soup.tri(face_index);
                let coordinates = triangle.map(|point| [point.x, point.y, point.z]);
                let outward = cross3(
                    subtract(coordinates[1], coordinates[0]),
                    subtract(coordinates[2], coordinates[0]),
                );
                if !(dot3(outward, outward).is_finite() && dot3(outward, outward) > 0.0) {
                    result.violations.push(violation(
                        "project-conduction-interface-source-face",
                        format!(
                            "entity `{}` selected degenerate or non-finite source face {} from `{}`",
                            entity.declared_target, face, artifact.artifact_role
                        ),
                        "repair the promoted triangle soup before assigning it to conduction semantics",
                    ));
                    return result;
                }
                index
                    .entry((
                        entity.declared_target.clone(),
                        coordinate_face_key(coordinates),
                    ))
                    .or_default()
                    .push(IndexedConductionSourceFace {
                        source: ConductionSourceFace {
                            artifact: artifact.artifact_role.clone(),
                            source_identity: artifact.source_identity.clone(),
                            face,
                        },
                        outward,
                    });
            }
        }
    }

    let candidates = match ThermalInterfaces::coincident_face_pairs(mesh) {
        Ok(candidates) => candidates,
        Err(error) => {
            result.violations.push(violation(
                "project-conduction-interface-geometry",
                format!("fs-conduction refused the exact boundary candidate set: {error}"),
                "repair the duplicated matching-P1 traces before lowering scenario interfaces",
            ));
            return result;
        }
    };
    let mut claimed_declarations = BTreeSet::new();
    let mut pending = Vec::new();
    for (candidate_index, candidate) in candidates.into_iter().enumerate() {
        if candidate_index % 1_024 == 0 && cx.checkpoint().is_err() {
            result.violations.push(conduction_interface_cancelled(
                "conduction-interface-candidate",
            ));
            return result;
        }
        let key = conduction_face_key(mesh, candidate.side_a);
        let mut claims = Vec::new();
        for &(interface, from, to) in &declarations {
            let direct_from = match unique_oriented_source(
                &index,
                from,
                key,
                mesh.boundary()[candidate.side_a].outward_normal,
            ) {
                Ok(source) => source,
                Err(violation) => {
                    result.violations.push(violation);
                    return result;
                }
            };
            let direct_to = match unique_oriented_source(
                &index,
                to,
                key,
                mesh.boundary()[candidate.side_b].outward_normal,
            ) {
                Ok(source) => source,
                Err(violation) => {
                    result.violations.push(violation);
                    return result;
                }
            };
            let reverse_from = match unique_oriented_source(
                &index,
                from,
                key,
                mesh.boundary()[candidate.side_b].outward_normal,
            ) {
                Ok(source) => source,
                Err(violation) => {
                    result.violations.push(violation);
                    return result;
                }
            };
            let reverse_to = match unique_oriented_source(
                &index,
                to,
                key,
                mesh.boundary()[candidate.side_a].outward_normal,
            ) {
                Ok(source) => source,
                Err(violation) => {
                    result.violations.push(violation);
                    return result;
                }
            };
            let orientation = match (direct_from, direct_to, reverse_from, reverse_to) {
                (Some(from_source), Some(to_source), None, None) => {
                    Some((candidate.side_a, candidate.side_b, from_source, to_source))
                }
                (None, None, Some(from_source), Some(to_source)) => {
                    Some((candidate.side_b, candidate.side_a, from_source, to_source))
                }
                (None, None, None, None) => None,
                _ => {
                    result.violations.push(violation(
                        "project-conduction-interface-orientation",
                        format!(
                            "declared interface `{interface}` does not uniquely orient coincident boundary slots {} and {} from `{from}` to `{to}`",
                            candidate.side_a, candidate.side_b
                        ),
                        "supply exactly one outward-oriented imported face for each declared region side",
                    ));
                    return result;
                }
            };
            let Some((from_slot, to_slot, from_source, to_source)) = orientation else {
                continue;
            };
            let Some(interface_rows) = index.get(&(interface.to_string(), key)) else {
                result.violations.push(violation(
                    "project-conduction-interface-selector-miss",
                    format!(
                        "interface `{interface}` joins the exact boundary triangle at slots {} and {}, but its geometry selector does not select that triangle",
                        candidate.side_a, candidate.side_b
                    ),
                    "make the interface selector cover the exact shared face while retaining explicit overlap with both region selectors",
                ));
                return result;
            };
            let mut interface_sources = interface_rows
                .iter()
                .map(|row| row.source.clone())
                .collect::<Vec<_>>();
            interface_sources.sort();
            interface_sources.dedup();
            if interface_sources.len() > 2 {
                result.violations.push(violation(
                    "project-conduction-interface-selector-ambiguous",
                    format!(
                        "interface `{interface}` selects {} source faces for one exact geometric triangle",
                        interface_sources.len()
                    ),
                    "select one shared CAD face or the two explicit duplicated traces; model higher-order junctions separately",
                ));
                return result;
            }
            claims.push(ResolvedConductionInterfacePair {
                interface: interface.to_string(),
                from_region: from.to_string(),
                to_region: to.to_string(),
                from_boundary_slot: from_slot,
                to_boundary_slot: to_slot,
                from_source,
                to_source,
                interface_sources,
            });
        }
        match claims.as_slice() {
            [] => {
                result.violations.push(violation(
                    "project-conduction-interface-undeclared",
                    format!(
                        "coincident conduction boundary slots {} and {} have no uniquely oriented declared scenario interface",
                        candidate.side_a, candidate.side_b
                    ),
                    "declare and geometrically assign the interface between the two owning regions; no contact law is inferred",
                ));
                return result;
            }
            [claim] => {
                claimed_declarations.insert(claim.interface.clone());
                pending.push(claim.clone());
            }
            many => {
                result.violations.push(violation(
                    "project-conduction-interface-multiple",
                    format!(
                        "coincident conduction boundary slots {} and {} are claimed by interfaces {}",
                        candidate.side_a,
                        candidate.side_b,
                        many.iter()
                            .map(|claim| claim.interface.as_str())
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    "bind each exact duplicated trace to one declared scenario interface",
                ));
                return result;
            }
        }
    }
    for &(interface, _, _) in &declarations {
        if !claimed_declarations.contains(interface) {
            result.violations.push(violation(
                "project-conduction-interface-no-face-pair",
                format!(
                    "declared interface `{interface}` did not lower to any exact coincident ConductionMesh boundary pair"
                ),
                "supply duplicated matching-P1 traces and outward-oriented region assignments for this interface",
            ));
        }
    }
    if result.violations.is_empty() {
        result.pairs = pending;
    }
    result
}

/// Detect resolved region surfaces that approach within the admitted tolerance
/// but have no interface entity joining their project-local names.
///
/// The audit first resolves assignments itself while holding one immutable
/// borrow of the mesh library. A caller therefore cannot substitute a stale
/// resolution or replace a mesh between selection and proximity testing. It
/// refuses incomparable coordinate units, invalid limits, cancellation, and
/// exhausted work budgets without publishing a partial contact list.
#[must_use]
pub fn audit_interface_declarations(
    spec: &ProjectSpec,
    library: &ImportedMeshLibrary,
    assignment_limits: AssignmentLimits,
    limits: InterfaceAuditLimits,
    cx: &Cx<'_>,
) -> InterfaceDeclarationAudit {
    let mut result = InterfaceDeclarationAudit::default();
    if cx.checkpoint().is_err() {
        result
            .violations
            .push(interface_audit_cancelled("interface-audit-entry"));
        return result;
    }
    let resolution = resolve_geometry_assignments(spec, library, assignment_limits, cx);
    if !resolution.admissible() {
        result.violations = resolution.violations;
        return result;
    }
    if !limits.proximity_tolerance.is_finite() || limits.proximity_tolerance < 0.0 {
        result.violations.push(violation(
            "project-interface-audit-tolerance",
            format!(
                "interface proximity tolerance {} is not finite and non-negative",
                limits.proximity_tolerance
            ),
            "supply one finite non-negative tolerance in the retained mesh coordinate unit",
        ));
        return result;
    }
    let Some(assembly) = &spec.assembly else {
        result.violations.push(violation(
            "project-interface-audit-preconditions",
            "interface proximity auditing requires the `assembly` section",
            "declare the assembly regions and their interface entities before auditing geometry",
        ));
        return result;
    };

    let mut declared_pairs = BTreeSet::new();
    for entity in assembly {
        if let EntityDecl::Interface { from, to, .. } = entity {
            declared_pairs.insert(ordered_pair(from, to));
        }
    }

    struct RegionSurface<'a> {
        name: &'a str,
        artifact: &'a str,
        length_unit: &'a str,
        soup: &'a Soup,
        faces: &'a [u32],
        bounds_min: [f64; 3],
        bounds_max: [f64; 3],
    }

    let mut surfaces = Vec::new();
    for artifact in &resolution.artifacts {
        let Some(imported) = library.get(&artifact.source_identity) else {
            result.violations.push(violation(
                "project-interface-audit-mesh-unavailable",
                format!(
                    "resolved geometry role `{}` no longer has its retained supplied mesh",
                    artifact.artifact_role
                ),
                "keep the exact imported-mesh library alive through assignment resolution and interface auditing",
            ));
            return result;
        };
        if artifact.entities.len() != artifact.report.assignments.len() {
            result.violations.push(violation(
                "project-interface-audit-report-mismatch",
                format!(
                    "resolved geometry role `{}` has different entity and assignment row counts",
                    artifact.artifact_role
                ),
                "discard the inconsistent resolution and rerun the project assignment adapter",
            ));
            return result;
        }
        for (entity, assignment) in artifact.entities.iter().zip(&artifact.report.assignments) {
            if entity.entity_id.kind() == EntityKind::Region {
                surfaces.push(RegionSurface {
                    name: &entity.declared_target,
                    artifact: &artifact.artifact_role,
                    length_unit: artifact.report.receipt.length_unit(),
                    soup: &imported.soup,
                    faces: &assignment.faces,
                    bounds_min: assignment.stats.bounds_min,
                    bounds_max: assignment.stats.bounds_max,
                });
            }
        }
    }
    surfaces.sort_by(|left, right| left.name.cmp(right.name));

    let tolerance_squared = limits.proximity_tolerance * limits.proximity_tolerance;
    let mut pending = Vec::new();
    for first_index in 0..surfaces.len() {
        for second_index in (first_index + 1)..surfaces.len() {
            let first = &surfaces[first_index];
            let second = &surfaces[second_index];
            if declared_pairs.contains(&ordered_pair(first.name, second.name)) {
                continue;
            }
            if first.length_unit != second.length_unit {
                result.violations.push(violation(
                    "project-interface-audit-unit-mismatch",
                    format!(
                        "undeclared region pair `{}` / `{}` uses incomparable coordinate units `{}` and `{}`",
                        first.name, second.name, first.length_unit, second.length_unit
                    ),
                    "promote both meshes into one explicit assembly coordinate unit before proximity auditing",
                ));
                return result;
            }
            if aabb_distance_squared(
                first.bounds_min,
                first.bounds_max,
                second.bounds_min,
                second.bounds_max,
            ) > tolerance_squared
            {
                continue;
            }

            let mut closest = f64::INFINITY;
            let mut closest_faces = (0, 0);
            for &first_face in first.faces {
                let first_triangle = first.soup.tri(first_face as usize);
                let first_coordinates = [
                    [
                        first_triangle[0].x,
                        first_triangle[0].y,
                        first_triangle[0].z,
                    ],
                    [
                        first_triangle[1].x,
                        first_triangle[1].y,
                        first_triangle[1].z,
                    ],
                    [
                        first_triangle[2].x,
                        first_triangle[2].y,
                        first_triangle[2].z,
                    ],
                ];
                for &second_face in second.faces {
                    if result.triangle_pair_tests == limits.max_triangle_pair_tests {
                        result.violations.push(violation(
                            "project-interface-audit-resource-bound",
                            format!(
                                "interface audit exhausted its {} triangle-pair tests before completing",
                                limits.max_triangle_pair_tests
                            ),
                            "raise the explicit interface-audit work budget or reduce the admitted contact surfaces",
                        ));
                        return result;
                    }
                    result.triangle_pair_tests += 1;
                    if result.triangle_pair_tests.is_multiple_of(1024) && cx.checkpoint().is_err() {
                        result
                            .violations
                            .push(interface_audit_cancelled("interface-audit-triangle-pair"));
                        return result;
                    }

                    let second_triangle = second.soup.tri(second_face as usize);
                    let second_coordinates = [
                        [
                            second_triangle[0].x,
                            second_triangle[0].y,
                            second_triangle[0].z,
                        ],
                        [
                            second_triangle[1].x,
                            second_triangle[1].y,
                            second_triangle[1].z,
                        ],
                        [
                            second_triangle[2].x,
                            second_triangle[2].y,
                            second_triangle[2].z,
                        ],
                    ];
                    let mut vertex_distance = f64::INFINITY;
                    for point in first_triangle {
                        vertex_distance = vertex_distance.min(point_triangle_distance(
                            point,
                            second_triangle[0],
                            second_triangle[1],
                            second_triangle[2],
                        ));
                    }
                    for point in second_triangle {
                        vertex_distance = vertex_distance.min(point_triangle_distance(
                            point,
                            first_triangle[0],
                            first_triangle[1],
                            first_triangle[2],
                        ));
                    }
                    let distance =
                        triangle_distance(first_coordinates, second_coordinates, vertex_distance);
                    if distance < closest {
                        closest = distance;
                        closest_faces = (first_face, second_face);
                    }
                }
            }
            if closest <= limits.proximity_tolerance {
                pending.push(UndeclaredContactCandidate {
                    first_region: first.name.to_string(),
                    second_region: second.name.to_string(),
                    first_artifact: first.artifact.to_string(),
                    second_artifact: second.artifact.to_string(),
                    first_face: closest_faces.0,
                    second_face: closest_faces.1,
                    separation: closest,
                    length_unit: first.length_unit.to_string(),
                });
            }
        }
        if first_index + 1 < surfaces.len() && cx.checkpoint().is_err() {
            result
                .violations
                .push(interface_audit_cancelled("interface-audit-region-pair"));
            return result;
        }
    }

    for contact in &pending {
        result.violations.push(violation(
            "project-interface-undeclared-contact",
            format!(
                "regions `{}` and `{}` approach to {} {} on faces {} / {} but have no declared interface",
                contact.first_region,
                contact.second_region,
                contact.separation,
                contact.length_unit,
                contact.first_face,
                contact.second_face
            ),
            "declare an interface entity and its deliberate thermal law; perfect contact must be explicit rather than inferred",
        ));
    }
    result.undeclared_contacts = pending;
    result
}

fn ordered_pair<'a>(first: &'a str, second: &'a str) -> (&'a str, &'a str) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn aabb_distance_squared(
    first_min: [f64; 3],
    first_max: [f64; 3],
    second_min: [f64; 3],
    second_max: [f64; 3],
) -> f64 {
    let mut squared = 0.0;
    for axis in 0..3 {
        let gap = if first_max[axis] < second_min[axis] {
            second_min[axis] - first_max[axis]
        } else if second_max[axis] < first_min[axis] {
            first_min[axis] - second_max[axis]
        } else {
            0.0
        };
        squared = gap.mul_add(gap, squared);
    }
    squared
}

fn triangle_distance(
    first_coordinates: [[f64; 3]; 3],
    second_coordinates: [[f64; 3]; 3],
    mut closest: f64,
) -> f64 {
    if triangles_intersect(first_coordinates, second_coordinates) {
        return 0.0;
    }
    for first_edge in 0..3 {
        for second_edge in 0..3 {
            closest = closest.min(segment_segment_distance(
                first_coordinates[first_edge],
                first_coordinates[(first_edge + 1) % 3],
                second_coordinates[second_edge],
                second_coordinates[(second_edge + 1) % 3],
            ));
        }
    }
    closest
}

fn triangles_intersect(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> bool {
    (0..3).any(|edge| {
        segment_intersects_triangle(first[edge], first[(edge + 1) % 3], second)
            || segment_intersects_triangle(second[edge], second[(edge + 1) % 3], first)
    })
}

fn segment_intersects_triangle(start: [f64; 3], end: [f64; 3], tri: [[f64; 3]; 3]) -> bool {
    let direction = subtract(end, start);
    let edge_1 = subtract(tri[1], tri[0]);
    let edge_2 = subtract(tri[2], tri[0]);
    let cross_direction = cross3(direction, edge_2);
    let determinant = dot3(edge_1, cross_direction);
    let scale = norm3(direction)
        .max(norm3(edge_1))
        .max(norm3(edge_2))
        .max(1.0);
    let determinant_epsilon = 32.0 * f64::EPSILON * scale * scale * scale;
    if determinant.abs() <= determinant_epsilon {
        return false;
    }
    let barycentric_epsilon = 32.0 * f64::EPSILON;
    let inverse = determinant.recip();
    let origin_delta = subtract(start, tri[0]);
    let u = inverse * dot3(origin_delta, cross_direction);
    if u < -barycentric_epsilon || u > 1.0 + barycentric_epsilon {
        return false;
    }
    let cross_origin = cross3(origin_delta, edge_1);
    let v = inverse * dot3(direction, cross_origin);
    if v < -barycentric_epsilon || u + v > 1.0 + barycentric_epsilon {
        return false;
    }
    let along_segment = inverse * dot3(edge_2, cross_origin);
    along_segment >= -barycentric_epsilon && along_segment <= 1.0 + barycentric_epsilon
}

fn segment_segment_distance(
    first_start: [f64; 3],
    first_end: [f64; 3],
    second_start: [f64; 3],
    second_end: [f64; 3],
) -> f64 {
    let first_direction = subtract(first_end, first_start);
    let second_direction = subtract(second_end, second_start);
    let origins = subtract(first_start, second_start);
    let first_length = dot3(first_direction, first_direction);
    let second_length = dot3(second_direction, second_direction);
    let second_projection = dot3(second_direction, origins);
    let epsilon = f64::EPSILON;

    let (mut first_parameter, mut second_parameter);
    if first_length <= epsilon && second_length <= epsilon {
        return norm3(origins);
    }
    if first_length <= epsilon {
        first_parameter = 0.0;
        second_parameter = (second_projection / second_length).clamp(0.0, 1.0);
    } else {
        let first_projection = dot3(first_direction, origins);
        if second_length <= epsilon {
            second_parameter = 0.0;
            first_parameter = (-first_projection / first_length).clamp(0.0, 1.0);
        } else {
            let coupling = dot3(first_direction, second_direction);
            let denominator = first_length * second_length - coupling * coupling;
            first_parameter = if denominator.abs() > epsilon {
                ((coupling * second_projection - first_projection * second_length) / denominator)
                    .clamp(0.0, 1.0)
            } else {
                0.0
            };
            second_parameter = (coupling * first_parameter + second_projection) / second_length;
            if second_parameter < 0.0 {
                second_parameter = 0.0;
                first_parameter = (-first_projection / first_length).clamp(0.0, 1.0);
            } else if second_parameter > 1.0 {
                second_parameter = 1.0;
                first_parameter = ((coupling - first_projection) / first_length).clamp(0.0, 1.0);
            }
        }
    }
    let first_closest = add(first_start, scale3(first_direction, first_parameter));
    let second_closest = add(second_start, scale3(second_direction, second_parameter));
    norm3(subtract(first_closest, second_closest))
}

fn add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn scale3(vector: [f64; 3], scale: f64) -> [f64; 3] {
    [vector[0] * scale, vector[1] * scale, vector[2] * scale]
}

fn dot3(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross3(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn norm3(vector: [f64; 3]) -> f64 {
    dot3(vector, vector).sqrt()
}

fn coordinate_point_key(point: [f64; 3]) -> ConductionPointKey {
    [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()]
}

fn coordinate_face_key(coordinates: [[f64; 3]; 3]) -> ConductionFaceKey {
    let mut key = coordinates.map(coordinate_point_key);
    key.sort_unstable();
    key
}

fn conduction_face_key(mesh: &ConductionMesh, slot: usize) -> ConductionFaceKey {
    let coordinates = mesh.boundary()[slot]
        .vertices
        .map(|vertex| mesh.positions()[vertex as usize]);
    coordinate_face_key(coordinates)
}

fn unique_oriented_source(
    index: &BTreeMap<(String, ConductionFaceKey), Vec<IndexedConductionSourceFace>>,
    entity: &str,
    key: ConductionFaceKey,
    outward_normal: [f64; 3],
) -> Result<Option<ConductionSourceFace>, Violation> {
    let Some(rows) = index.get(&(entity.to_string(), key)) else {
        return Ok(None);
    };
    let matching = rows
        .iter()
        .filter(|row| dot3(row.outward, outward_normal) > 0.0)
        .collect::<Vec<_>>();
    match matching.as_slice() {
        [] => Ok(None),
        [row] => Ok(Some(row.source.clone())),
        many => Err(violation(
            "project-conduction-interface-source-ambiguous",
            format!(
                "entity `{entity}` has {} outward-oriented imported faces for one exact conduction boundary slot",
                many.len()
            ),
            "remove duplicate source faces or make the selector identify one retained oriented trace",
        )),
    }
}

fn conduction_interface_cancelled(stage: &'static str) -> Violation {
    violation(
        "project-conduction-interface-cancelled",
        format!("conduction interface lowering was cancelled at `{stage}`"),
        "retry under a live deterministic cancellation context; no partial boundary-slot mapping was published",
    )
}

fn interface_audit_cancelled(stage: &'static str) -> Violation {
    violation(
        "project-interface-audit-cancelled",
        format!("interface declaration audit was cancelled at `{stage}`"),
        "retry under a live deterministic cancellation context; no partial contact list was published",
    )
}

fn cancelled(stage: &'static str) -> Violation {
    violation(
        "mesh-assignment-cancelled",
        format!("geometry assignment was cancelled at `{stage}`"),
        "retry under a live deterministic cancellation context; no partial report was published",
    )
}
