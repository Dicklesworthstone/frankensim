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
use fs_exec::Cx;
use fs_io::{
    AssignmentLimits, AssignmentReport, AssignmentRequest, NamedFaceGroup, resolve_mesh_assignments,
};
use fs_rep_mesh::Soup;
use fs_scenario::{EntityId, EntityKind, Violation};

use crate::spec::{GeometryArtifact, GeometryAssignment, ProjectSpec};

/// Domain for the identity hook derived from one exact project geometry row.
pub const GEOMETRY_SOURCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-project.geometry-source.v1";

/// Domain for retained exact fs-io assignment-report JSON bytes.
pub const GEOMETRY_ASSIGNMENT_REPORT_DOMAIN: &str =
    "org.frankensim.fs-project.geometry-assignment-report.v1";

const PROJECT_ASSIGNMENT_NO_CLAIM: &str = "the adapter binds exact project entity identities, one declared selector plan, and one supplied finite tessellation; it does not authenticate the mesh supplier, prove continuum/CAD/physical-region sameness, or make fs-io face ordinals stable across re-import";

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

fn cancelled(stage: &'static str) -> Violation {
    violation(
        "mesh-assignment-cancelled",
        format!("geometry assignment was cancelled at `{stage}`"),
        "retry under a live deterministic cancellation context; no partial report was published",
    )
}
