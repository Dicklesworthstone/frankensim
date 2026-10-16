//! Selector-anchored association and semantic revision diffs for promoted
//! triangle soups.
//!
//! This module compares two already-resolved [`AssignmentReport`] values. It
//! deliberately does not invent persistent CAD topology: subjects and source
//! identities remain caller-presented hooks, while geometric fingerprints are
//! deterministic diagnostics over finite tessellations.

use core::fmt;
use fs_exec::{Cx, ExecMode};
use fs_geom::Point3;
use fs_rep_mesh::Soup;
use std::fmt::Write as _;

use crate::{AssignmentReport, ResolvedAssignment, selection};

/// Versioned semantics for mesh association, revision diffs, and migration.
pub const MESH_ASSOCIATION_SEMANTICS_VERSION: &str = "fs-io/mesh-association/v1";

/// Maximum default owned-loop work between cancellation polls.
pub const MESH_ASSOCIATION_POLL_STRIDE: usize = 4096;

const ASSOCIATION_AUTHORITY: &str = "finite-tessellation-association-diagnostic";
const ASSOCIATION_NO_CLAIM: &str = "selector agreement, topology signatures, and geometric fingerprints are deterministic diagnostics over the two supplied finite tessellations; caller identities and the declared frame transform are retained but not authenticated; no persistent topological naming, CAD-semantic equivalence, continuum correspondence, collision-resistant identity, or safe automatic migration under adversarial topology change is certified";

/// A caller-declared proper rigid transform from the source frame into the
/// target frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidTransform3 {
    /// Row-major rotation matrix.
    pub rotation: [[f64; 3]; 3],
    /// Target-frame translation.
    pub translation: [f64; 3],
}

impl RigidTransform3 {
    /// Identity frame alignment.
    pub const IDENTITY: Self = Self {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0, 0.0, 0.0],
    };

    /// Construct an explicit source-to-target rigid transform.
    #[must_use]
    pub const fn new(rotation: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    fn apply(self, point: Point3) -> [f64; 3] {
        let source = [point.x, point.y, point.z];
        [
            dot(self.rotation[0], source) + self.translation[0],
            dot(self.rotation[1], source) + self.translation[1],
            dot(self.rotation[2], source) + self.translation[2],
        ]
    }
}

/// Explicit tolerances and resource limits for one association operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssociationPolicy {
    /// Stable relative surface-area drift.
    pub stable_relative_area: f64,
    /// Stable centroid drift in `length_unit`.
    pub stable_distance: f64,
    /// Stable maximum absolute orientation-moment drift.
    pub stable_orientation: f64,
    /// Stable relative axis-aligned extent drift in the target frame.
    pub stable_relative_extent: f64,
    /// Largest relative area drift admitted as a migration proposal.
    pub moved_relative_area: f64,
    /// Largest centroid drift admitted as a migration proposal.
    pub moved_distance: f64,
    /// Largest orientation-moment drift admitted as a migration proposal.
    pub moved_orientation: f64,
    /// Largest relative extent drift admitted as a migration proposal.
    pub moved_relative_extent: f64,
    /// Maximum normalized-score gap that makes two fallback candidates
    /// ambiguous.
    pub ambiguity_score_gap: f64,
    /// Maximum orthonormality/determinant residual for the declared rotation.
    pub frame_tolerance: f64,
    /// Maximum assignments on either side.
    pub max_assignments: usize,
    /// Maximum vertices in either supplied soup.
    pub max_mesh_vertices: usize,
    /// Maximum faces in either supplied soup.
    pub max_mesh_faces: usize,
    /// Maximum aggregate selected-face references across both sides.
    pub max_face_references: usize,
    /// Maximum aggregate edge records created for topology signatures.
    pub max_edge_records: usize,
    /// Maximum source-target candidate comparisons.
    pub max_candidate_tests: u64,
    /// Maximum owned-loop work between cancellation polls.
    pub poll_stride: usize,
}

impl AssociationPolicy {
    /// Construct a conservative policy around explicit length tolerances.
    ///
    /// The caller still owns the model-scale meaning of these distances. All
    /// relative thresholds and resource limits remain visible in the receipt.
    #[must_use]
    pub const fn engineering(stable_distance: f64, moved_distance: f64) -> Self {
        Self {
            stable_relative_area: 1.0e-9,
            stable_distance,
            stable_orientation: 1.0e-9,
            stable_relative_extent: 1.0e-9,
            moved_relative_area: 0.25,
            moved_distance,
            moved_orientation: 0.25,
            moved_relative_extent: 0.25,
            ambiguity_score_gap: 1.0e-9,
            frame_tolerance: 1.0e-10,
            max_assignments: 4096,
            max_mesh_vertices: 1_000_000,
            max_mesh_faces: 1_000_000,
            max_face_references: 4_000_000,
            max_edge_records: 12_000_000,
            max_candidate_tests: 16_000_000,
            poll_stride: MESH_ASSOCIATION_POLL_STRIDE,
        }
    }
}

/// Coarse topology features that are invariant to ordinary edge subdivision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopologySignature {
    /// Face-connected components, where faces meet across an edge.
    pub component_count: usize,
    /// Connected components of the selected boundary-edge graph.
    pub boundary_component_count: usize,
    /// Boundary components in which every vertex has degree two.
    pub boundary_loop_count: usize,
    /// Undirected edges incident to more than two selected faces.
    pub nonmanifold_edge_count: usize,
    /// Two-face edges with equal rather than opposite orientation.
    pub orientation_conflict_count: usize,
    /// True only when every selected edge is paired exactly once with opposite
    /// orientation.
    pub closed_oriented_boundary: bool,
}

/// A deterministic geometric diagnostic for one resolved assignment.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceFingerprint {
    /// Selected face count. This is receipted but is not a topology invariant.
    pub face_count: usize,
    /// Surface area in squared `length_unit`.
    pub surface_area: f64,
    /// Area-weighted centroid in the target comparison frame.
    pub centroid: [f64; 3],
    /// Axis-aligned extents in the target comparison frame.
    pub extents: [f64; 3],
    /// Normal-sign-invariant surface orientation tensor entries:
    /// `xx, yy, zz, xy, xz, yz`.
    pub orientation_moments: [f64; 6],
    /// Coarse selected-mesh topology.
    pub topology: TopologySignature,
    /// Local FNV-1a replay root over the fields above.
    pub local_fingerprint: u64,
}

/// Measured drift from one source fingerprint to one target candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssociationDrift {
    /// Whether the selector-semantic fingerprints are exactly equal.
    pub selector_agrees: bool,
    /// Whether the coarse topology signatures are exactly equal.
    pub topology_agrees: bool,
    /// Relative surface-area drift.
    pub relative_area: f64,
    /// Centroid distance in `length_unit`.
    pub centroid_distance: f64,
    /// Maximum absolute orientation-moment drift.
    pub orientation: f64,
    /// Maximum relative axis-aligned extent drift.
    pub relative_extent: f64,
    /// Maximum moved-threshold-normalized geometric drift.
    pub normalized_score: f64,
}

/// Association verdict for one source assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssociationVerdict {
    /// Same subject, selector, topology, and stable-threshold fingerprint.
    Stable,
    /// One candidate resolves within moved thresholds but requires review.
    Moved,
    /// Multiple fallback candidates are indistinguishable under the declared
    /// ambiguity gap.
    Ambiguous,
    /// No safe association is available.
    Lost,
}

impl AssociationVerdict {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Stable => "stable",
            Self::Moved => "moved",
            Self::Ambiguous => "ambiguous",
            Self::Lost => "lost",
        }
    }
}

/// Region-level semantic change classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionChange {
    /// No change beyond stable tolerances.
    Unchanged,
    /// Shape diagnostics remain stable but placement moved.
    Moved,
    /// Topology is stable while shape diagnostics changed.
    Deformed,
    /// The coarse topology signature changed.
    TopologyChanged,
    /// Source region disappeared.
    Removed,
    /// Multiple target candidates remain plausible.
    Ambiguous,
}

impl RegionChange {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Unchanged => "unchanged",
            Self::Moved => "moved",
            Self::Deformed => "deformed",
            Self::TopologyChanged => "topology-changed",
            Self::Removed => "removed",
            Self::Ambiguous => "ambiguous",
        }
    }
}

/// Assignment carry-over disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MigrationAction {
    /// Selector and fingerprint agreement permit deterministic carry-over.
    AutoApply,
    /// A unique moved/deformed candidate is proposed for explicit review.
    Propose,
    /// Migration must stop for user action.
    Refuse,
}

impl MigrationAction {
    const fn as_str(self) -> &'static str {
        match self {
            Self::AutoApply => "auto-apply",
            Self::Propose => "propose",
            Self::Refuse => "refuse",
        }
    }
}

/// One per-source association, semantic diff, and migration decision.
#[derive(Debug, Clone, PartialEq)]
pub struct AssociationDecision {
    /// Persistent source subject.
    pub source_subject: String,
    /// Chosen target subject for stable/moved decisions.
    pub target_subject: Option<String>,
    /// Association verdict.
    pub verdict: AssociationVerdict,
    /// Region-level semantic change.
    pub change: RegionChange,
    /// Carry-over action.
    pub migration: MigrationAction,
    /// Source geometric fingerprint.
    pub source_fingerprint: SurfaceFingerprint,
    /// Chosen target fingerprint, if any.
    pub target_fingerprint: Option<SurfaceFingerprint>,
    /// Drift to the chosen target, if any.
    pub drift: Option<AssociationDrift>,
    /// Canonically ordered fallback candidates when association is ambiguous.
    pub candidates: Vec<String>,
    /// Stable explanatory reason.
    pub reason: &'static str,
}

/// One target assignment that was not consumed by a source association.
#[derive(Debug, Clone, PartialEq)]
pub struct AddedRegion {
    /// Target subject.
    pub subject: String,
    /// Target fingerprint.
    pub fingerprint: SurfaceFingerprint,
}

/// Receipt binding all association inputs, policy, and report output.
#[derive(Debug, Clone, PartialEq)]
pub struct AssociationReceipt {
    source_identity: String,
    target_identity: String,
    length_unit: String,
    source_assignments_fingerprint: u64,
    target_assignments_fingerprint: u64,
    source_to_target: RigidTransform3,
    policy: AssociationPolicy,
    report_fingerprint: u64,
}

impl AssociationReceipt {
    /// Caller-presented source identity retained by the source assignment
    /// receipt.
    #[must_use]
    pub fn source_identity(&self) -> &str {
        &self.source_identity
    }

    /// Caller-presented target identity retained by the target assignment
    /// receipt.
    #[must_use]
    pub fn target_identity(&self) -> &str {
        &self.target_identity
    }

    /// Common declared length unit.
    #[must_use]
    pub fn length_unit(&self) -> &str {
        &self.length_unit
    }

    /// Local replay root over all published decisions and additions.
    #[must_use]
    pub const fn report_fingerprint(&self) -> u64 {
        self.report_fingerprint
    }

    /// Declared source-to-target frame transform.
    #[must_use]
    pub const fn source_to_target(&self) -> RigidTransform3 {
        self.source_to_target
    }

    /// Explicit policy used for this report.
    #[must_use]
    pub const fn policy(&self) -> AssociationPolicy {
        self.policy
    }
}

/// Atomic association, semantic-diff, and migration report.
#[derive(Debug, Clone, PartialEq)]
pub struct AssociationReport {
    /// Decisions in source assignment order.
    pub decisions: Vec<AssociationDecision>,
    /// Unconsumed target assignments in canonical subject order.
    pub added: Vec<AddedRegion>,
    /// Receipt binding the complete report.
    pub receipt: AssociationReceipt,
}

impl AssociationReport {
    /// Canonical one-line JSON artifact suitable for HELM-side storage and
    /// `frankensim compare` consumption.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut output = String::from("{\"kind\":\"mesh-association-report\",\"version\":");
        push_json_string(&mut output, MESH_ASSOCIATION_SEMANTICS_VERSION);
        output.push_str(",\"source_identity\":");
        push_json_string(&mut output, self.receipt.source_identity());
        output.push_str(",\"target_identity\":");
        push_json_string(&mut output, self.receipt.target_identity());
        output.push_str(",\"length_unit\":");
        push_json_string(&mut output, self.receipt.length_unit());
        let _ = write!(
            output,
            ",\"source_assignments_fingerprint\":\"{:016x}\",\"target_assignments_fingerprint\":\"{:016x}\"",
            self.receipt.source_assignments_fingerprint,
            self.receipt.target_assignments_fingerprint
        );
        push_transform_json(&mut output, self.receipt.source_to_target);
        push_policy_json(&mut output, self.receipt.policy);
        output.push_str(",\"decisions\":[");
        for (index, decision) in self.decisions.iter().enumerate() {
            if index != 0 {
                output.push(',');
            }
            push_decision_json(&mut output, decision);
        }
        output.push_str("],\"added\":[");
        for (index, added) in self.added.iter().enumerate() {
            if index != 0 {
                output.push(',');
            }
            output.push_str("{\"subject\":");
            push_json_string(&mut output, &added.subject);
            output.push_str(",\"fingerprint\":");
            push_surface_json(&mut output, &added.fingerprint);
            output.push('}');
        }
        let _ = write!(
            output,
            "],\"report_fingerprint\":\"{:016x}\",\"authority\":",
            self.receipt.report_fingerprint
        );
        push_json_string(&mut output, ASSOCIATION_AUTHORITY);
        output.push_str(",\"no_claim\":");
        push_json_string(&mut output, ASSOCIATION_NO_CLAIM);
        output.push('}');
        output
    }

    /// Deterministic human-readable region diff and migration table.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut output = String::from("# Imported geometry revision diff\n\n");
        let _ = writeln!(
            output,
            "- Source: `{}`",
            markdown_code(self.receipt.source_identity())
        );
        let _ = writeln!(
            output,
            "- Target: `{}`",
            markdown_code(self.receipt.target_identity())
        );
        let _ = writeln!(
            output,
            "- Length unit: `{}`",
            markdown_code(self.receipt.length_unit())
        );
        let _ = writeln!(
            output,
            "- Receipt: `{:016x}`\n",
            self.receipt.report_fingerprint
        );
        output.push_str(
            "| Source subject | Target subject | Association | Change | Migration | Centroid drift | Area drift |\n",
        );
        output.push_str("| --- | --- | --- | --- | --- | ---: | ---: |\n");
        for decision in &self.decisions {
            let target = decision.target_subject.as_deref().unwrap_or("—");
            let (distance, area) =
                decision
                    .drift
                    .map_or((String::from("—"), String::from("—")), |drift| {
                        (
                            drift.centroid_distance.to_string(),
                            drift.relative_area.to_string(),
                        )
                    });
            let _ = writeln!(
                output,
                "| {} | {} | {} | {} | {} | {} | {} |",
                markdown_cell(&decision.source_subject),
                markdown_cell(target),
                decision.verdict.as_str(),
                decision.change.as_str(),
                decision.migration.as_str(),
                distance,
                area
            );
        }
        for added in &self.added {
            let _ = writeln!(
                output,
                "| — | {} | added | added | refuse | — | — |",
                markdown_cell(&added.subject)
            );
        }
        output.push_str("\n_No claim: ");
        output.push_str(ASSOCIATION_NO_CLAIM);
        output.push_str("._\n");
        output
    }
}

/// Structured refusal before any report is published.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssociationRefusal {
    /// Stable machine-facing code.
    pub code: &'static str,
    /// Specific diagnosis.
    pub what: String,
    /// Actionable correction.
    pub fix: String,
}

impl fmt::Display for AssociationRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}: {} (fix: {})",
            self.code, self.what, self.fix
        )
    }
}

impl std::error::Error for AssociationRefusal {}

/// Associate resolved assignments across two promoted tessellations.
///
/// Exact subject identity always takes precedence and can never silently
/// rebound to a different target. When an exact subject is absent, only
/// selector-agreeing, topology-agreeing candidates inside moved thresholds are
/// considered. Stable decisions auto-apply; moved/deformed decisions are
/// proposals; ambiguity, loss, and topology change refuse migration.
///
/// # Errors
///
/// Returns [`AssociationRefusal`] for invalid policy/frame/unit inputs,
/// malformed report-to-soup references, resource-bound violations, numerical
/// overflow, non-deterministic execution, or cancellation. Publication is
/// atomic.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn associate_mesh_assignments(
    source_soup: &Soup,
    source: &AssignmentReport,
    target_soup: &Soup,
    target: &AssignmentReport,
    source_to_target: RigidTransform3,
    policy: AssociationPolicy,
    cx: &Cx<'_>,
) -> Result<AssociationReport, AssociationRefusal> {
    checkpoint(cx, "association-entry", 0)?;
    if cx.mode() != ExecMode::Deterministic {
        return Err(refusal(
            "mesh-association-fast-mode",
            "association receipts require deterministic execution mode",
            "retry with ExecMode::Deterministic",
        ));
    }
    validate_policy(policy)?;
    validate_transform(source_to_target, policy.frame_tolerance)?;
    if source.receipt.length_unit() != target.receipt.length_unit() {
        return Err(refusal(
            "mesh-association-unit-mismatch",
            format!(
                "source length unit {:?} differs from target length unit {:?}",
                source.receipt.length_unit(),
                target.receipt.length_unit()
            ),
            "convert both promoted tessellations into one explicit length unit before association",
        ));
    }
    validate_resource_envelope(source_soup, source, target_soup, target, policy)?;
    validate_subjects(&source.assignments, "source", policy, cx)?;
    validate_subjects(&target.assignments, "target", policy, cx)?;
    validate_receipt_binding(source_soup, source, "source", cx)?;
    validate_receipt_binding(target_soup, target, "target", cx)?;

    let mut source_fingerprints = Vec::new();
    source_fingerprints
        .try_reserve_exact(source.assignments.len())
        .map_err(|_| allocation_refusal("source fingerprints", source.assignments.len()))?;
    for (index, assignment) in source.assignments.iter().enumerate() {
        checkpoint(cx, "source-fingerprint", index)?;
        source_fingerprints.push(surface_fingerprint(
            source_soup,
            assignment,
            source_to_target,
            policy,
            cx,
            "source",
        )?);
    }

    let mut target_fingerprints = Vec::new();
    target_fingerprints
        .try_reserve_exact(target.assignments.len())
        .map_err(|_| allocation_refusal("target fingerprints", target.assignments.len()))?;
    for (index, assignment) in target.assignments.iter().enumerate() {
        checkpoint(cx, "target-fingerprint", index)?;
        target_fingerprints.push(surface_fingerprint(
            target_soup,
            assignment,
            RigidTransform3::IDENTITY,
            policy,
            cx,
            "target",
        )?);
    }

    let exact_targets: Vec<Option<usize>> = source
        .assignments
        .iter()
        .map(|source_assignment| {
            target
                .assignments
                .iter()
                .position(|candidate| candidate.subject == source_assignment.subject)
        })
        .collect();
    let mut exact_reserved = vec![false; target.assignments.len()];
    for target_index in exact_targets.iter().flatten().copied() {
        exact_reserved[target_index] = true;
    }
    let mut fallback_candidates = Vec::new();
    fallback_candidates
        .try_reserve_exact(source.assignments.len())
        .map_err(|_| allocation_refusal("fallback candidate rows", source.assignments.len()))?;
    fallback_candidates.resize_with(source.assignments.len(), Vec::new);
    let mut fallback_owner_count = vec![0usize; target.assignments.len()];
    let mut candidate_visits = 0usize;
    for (source_index, source_assignment) in source.assignments.iter().enumerate() {
        if exact_targets[source_index].is_some() {
            continue;
        }
        let source_fingerprint = &source_fingerprints[source_index];
        let candidates = &mut fallback_candidates[source_index];
        candidates
            .try_reserve_exact(target.assignments.len())
            .map_err(|_| allocation_refusal("fallback candidates", target.assignments.len()))?;
        for (target_index, target_assignment) in target.assignments.iter().enumerate() {
            if exact_reserved[target_index] {
                continue;
            }
            poll(
                cx,
                "association-candidate",
                candidate_visits,
                policy.poll_stride,
            )?;
            candidate_visits = candidate_visits.saturating_add(1);
            let drift = measure_drift(
                source_assignment,
                source_fingerprint,
                target_assignment,
                &target_fingerprints[target_index],
                policy,
            );
            if drift.selector_agrees && drift.topology_agrees && within_moved(drift, policy) {
                fallback_owner_count[target_index] =
                    fallback_owner_count[target_index].saturating_add(1);
                candidates.push((target_index, drift));
            }
        }
        candidates.sort_unstable_by(|left, right| {
            left.1
                .normalized_score
                .total_cmp(&right.1.normalized_score)
                .then_with(|| {
                    target.assignments[left.0]
                        .subject
                        .cmp(&target.assignments[right.0].subject)
                })
        });
    }

    let mut target_used = exact_reserved;
    let mut decisions = Vec::new();
    decisions
        .try_reserve_exact(source.assignments.len())
        .map_err(|_| allocation_refusal("association decisions", source.assignments.len()))?;
    for (source_index, source_assignment) in source.assignments.iter().enumerate() {
        checkpoint(cx, "association-decision", source_index)?;
        let source_fingerprint = &source_fingerprints[source_index];
        if let Some(target_index) = exact_targets[source_index] {
            let target_assignment = &target.assignments[target_index];
            let target_fingerprint = &target_fingerprints[target_index];
            let drift = measure_drift(
                source_assignment,
                source_fingerprint,
                target_assignment,
                target_fingerprint,
                policy,
            );
            let (verdict, change, migration, reason) = classify_exact(drift, policy);
            decisions.push(AssociationDecision {
                source_subject: source_assignment.subject.clone(),
                target_subject: Some(target_assignment.subject.clone()),
                verdict,
                change,
                migration,
                source_fingerprint: source_fingerprint.clone(),
                target_fingerprint: Some(target_fingerprint.clone()),
                drift: Some(drift),
                candidates: Vec::new(),
                reason,
            });
            continue;
        }

        let candidates = &fallback_candidates[source_index];
        let Some(&(best_index, best_drift)) = candidates.first() else {
            decisions.push(AssociationDecision {
                source_subject: source_assignment.subject.clone(),
                target_subject: None,
                verdict: AssociationVerdict::Lost,
                change: RegionChange::Removed,
                migration: MigrationAction::Refuse,
                source_fingerprint: source_fingerprint.clone(),
                target_fingerprint: None,
                drift: None,
                candidates: Vec::new(),
                reason: "no exact subject or selector-and-topology-compatible candidate falls within moved thresholds",
            });
            continue;
        };

        let score_ambiguous = candidates.get(1).is_some_and(|(_, second)| {
            second.normalized_score - best_drift.normalized_score <= policy.ambiguity_score_gap
        });
        let ownership_ambiguous = fallback_owner_count[best_index] > 1;
        if score_ambiguous || ownership_ambiguous {
            let threshold = best_drift.normalized_score + policy.ambiguity_score_gap;
            let mut names: Vec<String> = candidates
                .iter()
                .take_while(|(_, drift)| drift.normalized_score <= threshold)
                .map(|(target_index, _)| target.assignments[*target_index].subject.clone())
                .collect();
            names.sort_unstable();
            decisions.push(AssociationDecision {
                source_subject: source_assignment.subject.clone(),
                target_subject: None,
                verdict: AssociationVerdict::Ambiguous,
                change: RegionChange::Ambiguous,
                migration: MigrationAction::Refuse,
                source_fingerprint: source_fingerprint.clone(),
                target_fingerprint: None,
                drift: None,
                candidates: names,
                reason: if ownership_ambiguous {
                    "the best target candidate is compatible with multiple source subjects"
                } else {
                    "multiple selector-and-topology-compatible candidates fall within the declared ambiguity score gap"
                },
            });
            continue;
        }

        target_used[best_index] = true;
        let target_assignment = &target.assignments[best_index];
        let target_fingerprint = &target_fingerprints[best_index];
        decisions.push(AssociationDecision {
            source_subject: source_assignment.subject.clone(),
            target_subject: Some(target_assignment.subject.clone()),
            verdict: AssociationVerdict::Moved,
            change: classify_shape_change(best_drift, policy),
            migration: MigrationAction::Propose,
            source_fingerprint: source_fingerprint.clone(),
            target_fingerprint: Some(target_fingerprint.clone()),
            drift: Some(best_drift),
            candidates: Vec::new(),
            reason: "one fallback candidate agrees in selector and topology and falls within moved thresholds",
        });
    }

    let mut added = Vec::new();
    added
        .try_reserve_exact(target.assignments.len())
        .map_err(|_| allocation_refusal("added regions", target.assignments.len()))?;
    for (index, used) in target_used.iter().copied().enumerate() {
        poll(cx, "association-additions", index, policy.poll_stride)?;
        if !used {
            added.push(AddedRegion {
                subject: target.assignments[index].subject.clone(),
                fingerprint: target_fingerprints[index].clone(),
            });
        }
    }
    added.sort_unstable_by(|left, right| left.subject.cmp(&right.subject));
    checkpoint(cx, "association-publication", decisions.len())?;
    let report_fingerprint = fingerprint_report(
        source,
        target,
        source_to_target,
        policy,
        &decisions,
        &added,
        cx,
    )?;
    Ok(AssociationReport {
        decisions,
        added,
        receipt: AssociationReceipt {
            source_identity: source.receipt.source_identity().to_string(),
            target_identity: target.receipt.source_identity().to_string(),
            length_unit: source.receipt.length_unit().to_string(),
            source_assignments_fingerprint: source.receipt.assignments_fingerprint(),
            target_assignments_fingerprint: target.receipt.assignments_fingerprint(),
            source_to_target,
            policy,
            report_fingerprint,
        },
    })
}

fn validate_policy(policy: AssociationPolicy) -> Result<(), AssociationRefusal> {
    for (name, stable, moved) in [
        (
            "relative area",
            policy.stable_relative_area,
            policy.moved_relative_area,
        ),
        (
            "centroid distance",
            policy.stable_distance,
            policy.moved_distance,
        ),
        (
            "orientation",
            policy.stable_orientation,
            policy.moved_orientation,
        ),
        (
            "relative extent",
            policy.stable_relative_extent,
            policy.moved_relative_extent,
        ),
    ] {
        if !stable.is_finite() || stable < 0.0 || !moved.is_finite() || moved < stable {
            return Err(refusal(
                "mesh-association-invalid-policy",
                format!(
                    "{name} thresholds require finite 0 <= stable <= moved, got {stable} and {moved}"
                ),
                "supply ordered finite association thresholds",
            ));
        }
    }
    for (name, value) in [
        ("ambiguity_score_gap", policy.ambiguity_score_gap),
        ("frame_tolerance", policy.frame_tolerance),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(refusal(
                "mesh-association-invalid-policy",
                format!("{name} must be finite and non-negative, got {value}"),
                "supply a finite non-negative association policy value",
            ));
        }
    }
    if policy.max_assignments == 0
        || policy.max_mesh_vertices == 0
        || policy.max_mesh_faces == 0
        || policy.max_face_references == 0
        || policy.max_edge_records == 0
        || policy.max_candidate_tests == 0
        || policy.poll_stride == 0
        || policy.poll_stride > MESH_ASSOCIATION_POLL_STRIDE
    {
        return Err(refusal(
            "mesh-association-invalid-policy",
            format!(
                "resource limits must be nonzero and poll_stride must be in 1..={MESH_ASSOCIATION_POLL_STRIDE}"
            ),
            "supply a nonzero bounded association resource envelope",
        ));
    }
    Ok(())
}

fn validate_transform(
    transform: RigidTransform3,
    tolerance: f64,
) -> Result<(), AssociationRefusal> {
    if transform
        .rotation
        .iter()
        .flatten()
        .chain(transform.translation.iter())
        .any(|value| !value.is_finite())
    {
        return Err(refusal(
            "mesh-association-invalid-frame",
            "the source-to-target transform contains a non-finite value",
            "supply one finite proper rigid transform",
        ));
    }
    for row in 0..3 {
        for other in row..3 {
            let expected = if row == other { 1.0 } else { 0.0 };
            let residual =
                (dot(transform.rotation[row], transform.rotation[other]) - expected).abs();
            if residual > tolerance {
                return Err(refusal(
                    "mesh-association-invalid-frame",
                    format!(
                        "rotation row dot product ({row},{other}) has residual {residual}, above {tolerance}"
                    ),
                    "orthonormalize the declared rotation or raise the receipted frame tolerance",
                ));
            }
        }
    }
    let determinant = determinant(transform.rotation);
    if (determinant - 1.0).abs() > tolerance {
        return Err(refusal(
            "mesh-association-invalid-frame",
            format!("rotation determinant {determinant} is not +1 within tolerance {tolerance}"),
            "supply a proper rotation without reflection or scale",
        ));
    }
    Ok(())
}

fn validate_resource_envelope(
    source_soup: &Soup,
    source: &AssignmentReport,
    target_soup: &Soup,
    target: &AssignmentReport,
    policy: AssociationPolicy,
) -> Result<(), AssociationRefusal> {
    for (side, count) in [
        ("source assignments", source.assignments.len()),
        ("target assignments", target.assignments.len()),
    ] {
        if count > policy.max_assignments {
            return Err(resource_refusal(side, count, policy.max_assignments));
        }
    }
    for (side, soup) in [("source", source_soup), ("target", target_soup)] {
        if soup.positions.len() > policy.max_mesh_vertices {
            return Err(resource_refusal(
                &format!("{side} mesh vertices"),
                soup.positions.len(),
                policy.max_mesh_vertices,
            ));
        }
        if soup.triangles.len() > policy.max_mesh_faces {
            return Err(resource_refusal(
                &format!("{side} mesh faces"),
                soup.triangles.len(),
                policy.max_mesh_faces,
            ));
        }
    }
    let source_faces = sum_faces(&source.assignments)?;
    let target_faces = sum_faces(&target.assignments)?;
    let face_references = source_faces.checked_add(target_faces).ok_or_else(|| {
        resource_refusal(
            "aggregate selected-face references",
            usize::MAX,
            policy.max_face_references,
        )
    })?;
    if face_references > policy.max_face_references {
        return Err(resource_refusal(
            "aggregate selected-face references",
            face_references,
            policy.max_face_references,
        ));
    }
    let edge_records = face_references.checked_mul(3).ok_or_else(|| {
        resource_refusal(
            "aggregate topology edge records",
            usize::MAX,
            policy.max_edge_records,
        )
    })?;
    if edge_records > policy.max_edge_records {
        return Err(resource_refusal(
            "aggregate topology edge records",
            edge_records,
            policy.max_edge_records,
        ));
    }
    let candidate_tests = u64::try_from(source.assignments.len())
        .ok()
        .and_then(|left| {
            u64::try_from(target.assignments.len())
                .ok()
                .and_then(|right| left.checked_mul(right))
        })
        .ok_or_else(|| {
            resource_refusal(
                "source-target candidate tests",
                usize::MAX,
                usize::try_from(policy.max_candidate_tests).unwrap_or(usize::MAX),
            )
        })?;
    if candidate_tests > policy.max_candidate_tests {
        return Err(refusal(
            "mesh-association-work-limit",
            format!(
                "association requires {candidate_tests} candidate tests; cap is {}",
                policy.max_candidate_tests
            ),
            "reduce the assignment sets or explicitly raise max_candidate_tests",
        ));
    }
    Ok(())
}

fn sum_faces(assignments: &[ResolvedAssignment]) -> Result<usize, AssociationRefusal> {
    assignments.iter().try_fold(0usize, |sum, assignment| {
        sum.checked_add(assignment.faces.len()).ok_or_else(|| {
            refusal(
                "mesh-association-work-overflow",
                "selected-face reference count overflowed",
                "reduce the association input",
            )
        })
    })
}

fn validate_subjects(
    assignments: &[ResolvedAssignment],
    side: &str,
    policy: AssociationPolicy,
    cx: &Cx<'_>,
) -> Result<(), AssociationRefusal> {
    let mut subjects = Vec::new();
    subjects
        .try_reserve_exact(assignments.len())
        .map_err(|_| allocation_refusal("association subject index", assignments.len()))?;
    for (index, assignment) in assignments.iter().enumerate() {
        poll(cx, "association-subjects", index, policy.poll_stride)?;
        if assignment.subject.is_empty()
            || assignment.subject.trim() != assignment.subject
            || assignment.subject.chars().any(char::is_control)
        {
            return Err(refusal(
                "mesh-association-invalid-subject",
                format!("{side} subject {:?} is not canonical", assignment.subject),
                "re-resolve assignments with nonempty trim-canonical control-free subjects",
            ));
        }
        subjects.push(assignment.subject.as_str());
    }
    subjects.sort_unstable();
    if let Some(duplicate) = subjects.windows(2).find(|pair| pair[0] == pair[1]) {
        return Err(refusal(
            "mesh-association-duplicate-subject",
            format!("{side} subject {:?} occurs more than once", duplicate[0]),
            "publish exactly one resolved assignment per persistent subject",
        ));
    }
    Ok(())
}

fn validate_receipt_binding(
    soup: &Soup,
    report: &AssignmentReport,
    side: &str,
    cx: &Cx<'_>,
) -> Result<(), AssociationRefusal> {
    let soup_fingerprint = selection::fingerprint_soup(soup, cx)
        .map_err(|error| map_assignment_refusal(error, side))?;
    if soup_fingerprint != report.receipt.source_mesh_fingerprint() {
        return Err(refusal(
            "mesh-association-soup-receipt-mismatch",
            format!(
                "{side} soup fingerprint {soup_fingerprint:016x} differs from assignment receipt {:016x}",
                report.receipt.source_mesh_fingerprint()
            ),
            "pair each assignment report with the exact promoted soup it resolved",
        ));
    }
    let assignments_fingerprint = selection::fingerprint_assignments(&report.assignments, cx)
        .map_err(|error| map_assignment_refusal(error, side))?;
    if assignments_fingerprint != report.receipt.assignments_fingerprint() {
        return Err(refusal(
            "mesh-association-assignment-receipt-mismatch",
            format!(
                "{side} assignment rows fingerprint {assignments_fingerprint:016x} differs from receipt {:016x}",
                report.receipt.assignments_fingerprint()
            ),
            "use the immutable assignment rows originally published with the receipt",
        ));
    }
    Ok(())
}

fn map_assignment_refusal(error: crate::AssignmentRefusal, side: &str) -> AssociationRefusal {
    if error.code == "mesh-assignment-cancelled" {
        refusal(
            "mesh-association-cancelled",
            format!("cancellation observed while binding the {side} assignment receipt"),
            "retry under a fresh scope from the same immutable inputs",
        )
    } else {
        refusal(
            "mesh-association-assignment-receipt",
            format!("{side} assignment receipt could not be rebound: {error}"),
            "re-run assignment resolution on the exact promoted soup before association",
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeRecord {
    low: u32,
    high: u32,
    face: usize,
    forward: bool,
}

fn surface_fingerprint(
    soup: &Soup,
    assignment: &ResolvedAssignment,
    transform: RigidTransform3,
    policy: AssociationPolicy,
    cx: &Cx<'_>,
    side: &str,
) -> Result<SurfaceFingerprint, AssociationRefusal> {
    if assignment.faces.is_empty() {
        return Err(refusal(
            "mesh-association-empty-assignment",
            format!("{side} subject {:?} selects no faces", assignment.subject),
            "re-run assignment resolution and refuse empty selections",
        ));
    }
    if assignment.faces.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(refusal(
            "mesh-association-face-order",
            format!(
                "{side} subject {:?} faces are not sorted and unique",
                assignment.subject
            ),
            "re-run deterministic assignment resolution before association",
        ));
    }
    let edge_capacity = assignment.faces.len().checked_mul(3).ok_or_else(|| {
        resource_refusal(
            "assignment topology edge records",
            usize::MAX,
            policy.max_edge_records,
        )
    })?;
    let mut edges = Vec::new();
    edges
        .try_reserve_exact(edge_capacity)
        .map_err(|_| allocation_refusal("assignment topology edge records", edge_capacity))?;
    let mut area = 0.0;
    let mut centroid_numerator = [0.0; 3];
    let mut orientation = [0.0; 6];
    let mut bounds_min = [f64::INFINITY; 3];
    let mut bounds_max = [f64::NEG_INFINITY; 3];

    for (local_face, face) in assignment.faces.iter().copied().enumerate() {
        poll(cx, "association-face", local_face, policy.poll_stride)?;
        let face_index = usize::try_from(face).map_err(|_| {
            refusal(
                "mesh-association-face-range",
                format!("{side} face {face} is not representable"),
                "repair the assignment report",
            )
        })?;
        let Some(indices) = soup.triangles.get(face_index).copied() else {
            return Err(refusal(
                "mesh-association-face-range",
                format!(
                    "{side} subject {:?} references missing face {face}",
                    assignment.subject
                ),
                "pair the assignment report with the exact promoted soup it resolved",
            ));
        };
        let mut points = [[0.0; 3]; 3];
        for (corner, vertex) in indices.iter().copied().enumerate() {
            let vertex_index = usize::try_from(vertex).map_err(|_| {
                refusal(
                    "mesh-association-vertex-range",
                    format!("{side} face {face} references unrepresentable vertex {vertex}"),
                    "repair the promoted soup before association",
                )
            })?;
            let Some(point) = soup.positions.get(vertex_index).copied() else {
                return Err(refusal(
                    "mesh-association-vertex-range",
                    format!("{side} face {face} references missing vertex {vertex}"),
                    "repair the promoted soup before association",
                ));
            };
            if !finite3([point.x, point.y, point.z]) {
                return Err(refusal(
                    "mesh-association-nonfinite-vertex",
                    format!("{side} face {face} references a non-finite vertex"),
                    "refuse or repair non-finite geometry before association",
                ));
            }
            let transformed = transform.apply(point);
            if !finite3(transformed) {
                return Err(refusal(
                    "mesh-association-frame-overflow",
                    format!(
                        "{side} face {face} overflows under the declared source-to-target transform"
                    ),
                    "rescale the geometry or supply a numerically representable frame transform",
                ));
            }
            points[corner] = transformed;
            for axis in 0..3 {
                bounds_min[axis] = bounds_min[axis].min(transformed[axis]);
                bounds_max[axis] = bounds_max[axis].max(transformed[axis]);
            }
        }
        let normal = cross(sub(points[1], points[0]), sub(points[2], points[0]));
        let twice_area = norm(normal);
        let triangle_area = 0.5 * twice_area;
        if !triangle_area.is_finite() || triangle_area <= 0.0 {
            return Err(refusal(
                "mesh-association-degenerate-face",
                format!("{side} face {face} has non-positive or overflowing area"),
                "repair degenerate geometry before association",
            ));
        }
        area += triangle_area;
        let triangle_centroid = scale(add(add(points[0], points[1]), points[2]), 1.0 / 3.0);
        for axis in 0..3 {
            centroid_numerator[axis] += triangle_area * triangle_centroid[axis];
        }
        let unit = scale(normal, 1.0 / twice_area);
        for (slot, value) in [
            unit[0] * unit[0],
            unit[1] * unit[1],
            unit[2] * unit[2],
            unit[0] * unit[1],
            unit[0] * unit[2],
            unit[1] * unit[2],
        ]
        .into_iter()
        .enumerate()
        {
            orientation[slot] += triangle_area * value;
        }
        for (from, to) in [
            (indices[0], indices[1]),
            (indices[1], indices[2]),
            (indices[2], indices[0]),
        ] {
            edges.push(EdgeRecord {
                low: from.min(to),
                high: from.max(to),
                face: local_face,
                forward: from <= to,
            });
        }
    }
    if !area.is_finite()
        || !finite3(centroid_numerator)
        || orientation.iter().any(|value| !value.is_finite())
    {
        return Err(refusal(
            "mesh-association-statistics-overflow",
            format!(
                "{side} subject {:?} fingerprint statistics overflowed",
                assignment.subject
            ),
            "rescale or split the selected geometry",
        ));
    }
    let mut centroid = scale(centroid_numerator, 1.0 / area);
    let mut extents = sub(bounds_max, bounds_min);
    for value in &mut orientation {
        *value = canonical_zero(*value / area);
    }
    for value in centroid.iter_mut().chain(extents.iter_mut()) {
        *value = canonical_zero(*value);
    }
    let topology = topology_signature(
        &mut edges,
        soup.positions.len(),
        assignment.faces.len(),
        policy,
        cx,
    )?;
    let mut fingerprint = Fingerprint::new();
    fingerprint.absorb_usize(assignment.faces.len());
    fingerprint.absorb_f64(area);
    fingerprint.absorb_f64s(&centroid);
    fingerprint.absorb_f64s(&extents);
    fingerprint.absorb_f64s(&orientation);
    fingerprint.absorb_usize(topology.component_count);
    fingerprint.absorb_usize(topology.boundary_component_count);
    fingerprint.absorb_usize(topology.boundary_loop_count);
    fingerprint.absorb_usize(topology.nonmanifold_edge_count);
    fingerprint.absorb_usize(topology.orientation_conflict_count);
    fingerprint.absorb_bool(topology.closed_oriented_boundary);
    Ok(SurfaceFingerprint {
        face_count: assignment.faces.len(),
        surface_area: canonical_zero(area),
        centroid,
        extents,
        orientation_moments: orientation,
        topology,
        local_fingerprint: fingerprint.finish(),
    })
}

fn topology_signature(
    edges: &mut [EdgeRecord],
    vertex_count: usize,
    face_count: usize,
    policy: AssociationPolicy,
    cx: &Cx<'_>,
) -> Result<TopologySignature, AssociationRefusal> {
    edges.sort_unstable();
    let mut face_sets = DisjointSet::new(face_count)?;
    let mut boundary_edges = Vec::new();
    boundary_edges
        .try_reserve_exact(edges.len())
        .map_err(|_| allocation_refusal("boundary edge candidates", edges.len()))?;
    let mut nonmanifold_edge_count = 0usize;
    let mut orientation_conflict_count = 0usize;
    let mut start = 0usize;
    while start < edges.len() {
        poll(cx, "association-edge-runs", start, policy.poll_stride)?;
        let key = (edges[start].low, edges[start].high);
        let mut end = start + 1;
        while end < edges.len() && (edges[end].low, edges[end].high) == key {
            end += 1;
        }
        for position in start + 1..end {
            face_sets.union(edges[start].face, edges[position].face);
        }
        match end - start {
            1 => boundary_edges.push(key),
            2 => {
                if edges[start].forward == edges[start + 1].forward {
                    orientation_conflict_count += 1;
                }
            }
            _ => nonmanifold_edge_count += 1,
        }
        start = end;
    }
    let component_count = face_sets.root_count();
    let (boundary_component_count, boundary_loop_count) =
        boundary_signature(&boundary_edges, vertex_count, policy, cx)?;
    Ok(TopologySignature {
        component_count,
        boundary_component_count,
        boundary_loop_count,
        nonmanifold_edge_count,
        orientation_conflict_count,
        closed_oriented_boundary: boundary_edges.is_empty()
            && nonmanifold_edge_count == 0
            && orientation_conflict_count == 0,
    })
}

fn boundary_signature(
    edges: &[(u32, u32)],
    vertex_count: usize,
    policy: AssociationPolicy,
    cx: &Cx<'_>,
) -> Result<(usize, usize), AssociationRefusal> {
    if edges.is_empty() {
        return Ok((0, 0));
    }
    let endpoint_capacity = edges.len().checked_mul(2).ok_or_else(|| {
        resource_refusal(
            "boundary endpoints",
            usize::MAX,
            policy.max_edge_records.saturating_mul(2),
        )
    })?;
    let mut vertices = Vec::new();
    vertices
        .try_reserve_exact(endpoint_capacity)
        .map_err(|_| allocation_refusal("boundary endpoints", endpoint_capacity))?;
    for (index, &(left, right)) in edges.iter().enumerate() {
        poll(
            cx,
            "association-boundary-endpoints",
            index,
            policy.poll_stride,
        )?;
        for vertex in [left, right] {
            if usize::try_from(vertex).map_or(true, |vertex| vertex >= vertex_count) {
                return Err(refusal(
                    "mesh-association-vertex-range",
                    format!("boundary edge references missing vertex {vertex}"),
                    "repair the promoted soup before association",
                ));
            }
            vertices.push(vertex);
        }
    }
    vertices.sort_unstable();
    vertices.dedup();
    let mut sets = DisjointSet::new(vertices.len())?;
    let mut degrees = vec![0usize; vertices.len()];
    for (index, &(left, right)) in edges.iter().enumerate() {
        poll(cx, "association-boundary-graph", index, policy.poll_stride)?;
        let left = vertices.binary_search(&left).map_err(|_| {
            refusal(
                "mesh-association-internal-boundary",
                "a boundary endpoint was absent from its canonical vertex index",
                "retry from the exact promoted soup and report",
            )
        })?;
        let right = vertices.binary_search(&right).map_err(|_| {
            refusal(
                "mesh-association-internal-boundary",
                "a boundary endpoint was absent from its canonical vertex index",
                "retry from the exact promoted soup and report",
            )
        })?;
        degrees[left] = degrees[left].saturating_add(1);
        degrees[right] = degrees[right].saturating_add(1);
        sets.union(left, right);
    }
    let mut roots = Vec::new();
    roots
        .try_reserve_exact(vertices.len())
        .map_err(|_| allocation_refusal("boundary component roots", vertices.len()))?;
    for vertex in 0..vertices.len() {
        roots.push(sets.find(vertex));
    }
    let mut component_roots = roots.clone();
    component_roots.sort_unstable();
    component_roots.dedup();
    let mut all_degree_two = vec![true; vertices.len()];
    for (vertex, root) in roots.into_iter().enumerate() {
        if degrees[vertex] != 2 {
            all_degree_two[root] = false;
        }
    }
    let boundary_loop_count = component_roots
        .iter()
        .filter(|root| all_degree_two[**root])
        .count();
    Ok((component_roots.len(), boundary_loop_count))
}

fn measure_drift(
    source_assignment: &ResolvedAssignment,
    source: &SurfaceFingerprint,
    target_assignment: &ResolvedAssignment,
    target: &SurfaceFingerprint,
    policy: AssociationPolicy,
) -> AssociationDrift {
    let relative_area = relative_difference(source.surface_area, target.surface_area);
    let centroid_distance = norm(sub(target.centroid, source.centroid));
    let orientation = source
        .orientation_moments
        .iter()
        .zip(target.orientation_moments)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max);
    let relative_extent = source
        .extents
        .iter()
        .zip(target.extents)
        .map(|(left, right)| relative_difference(*left, right))
        .fold(0.0, f64::max);
    let normalized_score = [
        threshold_ratio(relative_area, policy.moved_relative_area),
        threshold_ratio(centroid_distance, policy.moved_distance),
        threshold_ratio(orientation, policy.moved_orientation),
        threshold_ratio(relative_extent, policy.moved_relative_extent),
    ]
    .into_iter()
    .fold(0.0, f64::max);
    AssociationDrift {
        selector_agrees: source_assignment.selector_fingerprint
            == target_assignment.selector_fingerprint,
        topology_agrees: source.topology == target.topology,
        relative_area,
        centroid_distance,
        orientation,
        relative_extent,
        normalized_score,
    }
}

fn classify_exact(
    drift: AssociationDrift,
    policy: AssociationPolicy,
) -> (
    AssociationVerdict,
    RegionChange,
    MigrationAction,
    &'static str,
) {
    if !drift.topology_agrees {
        return (
            AssociationVerdict::Lost,
            RegionChange::TopologyChanged,
            MigrationAction::Refuse,
            "the exact subject resolves but its coarse topology signature changed",
        );
    }
    if drift.selector_agrees && within_stable(drift, policy) {
        return (
            AssociationVerdict::Stable,
            RegionChange::Unchanged,
            MigrationAction::AutoApply,
            "exact subject, selector, topology, and geometric fingerprint agree within stable thresholds",
        );
    }
    if within_moved(drift, policy) {
        return (
            AssociationVerdict::Moved,
            classify_shape_change(drift, policy),
            MigrationAction::Propose,
            "the exact subject remains unique but selector or geometric drift requires review",
        );
    }
    (
        AssociationVerdict::Lost,
        RegionChange::Deformed,
        MigrationAction::Refuse,
        "the exact subject resolves outside the declared moved thresholds",
    )
}

fn classify_shape_change(drift: AssociationDrift, policy: AssociationPolicy) -> RegionChange {
    if drift.relative_area <= policy.stable_relative_area
        && drift.orientation <= policy.stable_orientation
        && drift.relative_extent <= policy.stable_relative_extent
    {
        RegionChange::Moved
    } else {
        RegionChange::Deformed
    }
}

fn within_stable(drift: AssociationDrift, policy: AssociationPolicy) -> bool {
    drift.topology_agrees
        && drift.relative_area <= policy.stable_relative_area
        && drift.centroid_distance <= policy.stable_distance
        && drift.orientation <= policy.stable_orientation
        && drift.relative_extent <= policy.stable_relative_extent
}

fn within_moved(drift: AssociationDrift, policy: AssociationPolicy) -> bool {
    drift.topology_agrees
        && drift.relative_area <= policy.moved_relative_area
        && drift.centroid_distance <= policy.moved_distance
        && drift.orientation <= policy.moved_orientation
        && drift.relative_extent <= policy.moved_relative_extent
}

fn relative_difference(left: f64, right: f64) -> f64 {
    let scale = left.abs().max(right.abs());
    if scale == 0.0 {
        0.0
    } else {
        (left - right).abs() / scale
    }
}

fn threshold_ratio(value: f64, threshold: f64) -> f64 {
    if threshold == 0.0 {
        if value == 0.0 { 0.0 } else { f64::MAX }
    } else {
        value / threshold
    }
}

struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(length: usize) -> Result<Self, AssociationRefusal> {
        let mut parent = Vec::new();
        parent
            .try_reserve_exact(length)
            .map_err(|_| allocation_refusal("disjoint-set parents", length))?;
        parent.extend(0..length);
        let mut rank = Vec::new();
        rank.try_reserve_exact(length)
            .map_err(|_| allocation_refusal("disjoint-set ranks", length))?;
        rank.resize(length, 0);
        Ok(Self { parent, rank })
    }

    fn find(&mut self, mut value: usize) -> usize {
        let mut root = value;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        while self.parent[value] != value {
            let next = self.parent[value];
            self.parent[value] = root;
            value = next;
        }
        root
    }

    fn union(&mut self, left: usize, right: usize) {
        let mut left = self.find(left);
        let mut right = self.find(right);
        if left == right {
            return;
        }
        if self.rank[left] < self.rank[right] {
            core::mem::swap(&mut left, &mut right);
        }
        self.parent[right] = left;
        if self.rank[left] == self.rank[right] {
            self.rank[left] = self.rank[left].saturating_add(1);
        }
    }

    fn root_count(&mut self) -> usize {
        (0..self.parent.len())
            .filter(|index| self.find(*index) == *index)
            .count()
    }
}

fn fingerprint_report(
    source: &AssignmentReport,
    target: &AssignmentReport,
    source_to_target: RigidTransform3,
    policy: AssociationPolicy,
    decisions: &[AssociationDecision],
    added: &[AddedRegion],
    cx: &Cx<'_>,
) -> Result<u64, AssociationRefusal> {
    let mut fingerprint = Fingerprint::new();
    fingerprint.absorb_bytes(MESH_ASSOCIATION_SEMANTICS_VERSION.as_bytes());
    fingerprint.absorb_bytes(source.receipt.source_identity().as_bytes());
    fingerprint.absorb_bytes(target.receipt.source_identity().as_bytes());
    fingerprint.absorb_bytes(source.receipt.length_unit().as_bytes());
    fingerprint.absorb_u64(source.receipt.assignments_fingerprint());
    fingerprint.absorb_u64(target.receipt.assignments_fingerprint());
    for row in source_to_target.rotation {
        fingerprint.absorb_f64s(&row);
    }
    fingerprint.absorb_f64s(&source_to_target.translation);
    fingerprint.absorb_f64(policy.stable_relative_area);
    fingerprint.absorb_f64(policy.stable_distance);
    fingerprint.absorb_f64(policy.stable_orientation);
    fingerprint.absorb_f64(policy.stable_relative_extent);
    fingerprint.absorb_f64(policy.moved_relative_area);
    fingerprint.absorb_f64(policy.moved_distance);
    fingerprint.absorb_f64(policy.moved_orientation);
    fingerprint.absorb_f64(policy.moved_relative_extent);
    fingerprint.absorb_f64(policy.ambiguity_score_gap);
    fingerprint.absorb_f64(policy.frame_tolerance);
    fingerprint.absorb_usize(policy.max_assignments);
    fingerprint.absorb_usize(policy.max_mesh_vertices);
    fingerprint.absorb_usize(policy.max_mesh_faces);
    fingerprint.absorb_usize(policy.max_face_references);
    fingerprint.absorb_usize(policy.max_edge_records);
    fingerprint.absorb_u64(policy.max_candidate_tests);
    fingerprint.absorb_usize(policy.poll_stride);
    fingerprint.absorb_usize(decisions.len());
    for (index, decision) in decisions.iter().enumerate() {
        poll(
            cx,
            "association-output-fingerprint",
            index,
            policy.poll_stride,
        )?;
        fingerprint.absorb_bytes(decision.source_subject.as_bytes());
        fingerprint.absorb_bool(decision.target_subject.is_some());
        fingerprint.absorb_bytes(
            decision
                .target_subject
                .as_deref()
                .unwrap_or_default()
                .as_bytes(),
        );
        fingerprint.absorb_bytes(decision.verdict.as_str().as_bytes());
        fingerprint.absorb_bytes(decision.change.as_str().as_bytes());
        fingerprint.absorb_bytes(decision.migration.as_str().as_bytes());
        fingerprint.absorb_u64(decision.source_fingerprint.local_fingerprint);
        fingerprint.absorb_u64(
            decision
                .target_fingerprint
                .as_ref()
                .map_or(0, |surface| surface.local_fingerprint),
        );
        fingerprint.absorb_bool(decision.drift.is_some());
        if let Some(drift) = decision.drift {
            fingerprint.absorb_bool(drift.selector_agrees);
            fingerprint.absorb_bool(drift.topology_agrees);
            fingerprint.absorb_f64(drift.relative_area);
            fingerprint.absorb_f64(drift.centroid_distance);
            fingerprint.absorb_f64(drift.orientation);
            fingerprint.absorb_f64(drift.relative_extent);
            fingerprint.absorb_f64(drift.normalized_score);
        }
        fingerprint.absorb_usize(decision.candidates.len());
        for candidate in &decision.candidates {
            fingerprint.absorb_bytes(candidate.as_bytes());
        }
        fingerprint.absorb_bytes(decision.reason.as_bytes());
    }
    fingerprint.absorb_usize(added.len());
    for (index, added) in added.iter().enumerate() {
        poll(
            cx,
            "association-added-fingerprint",
            index,
            policy.poll_stride,
        )?;
        fingerprint.absorb_bytes(added.subject.as_bytes());
        fingerprint.absorb_u64(added.fingerprint.local_fingerprint);
    }
    Ok(fingerprint.finish())
}

struct Fingerprint(u64);

impl Fingerprint {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    const fn new() -> Self {
        Self(Self::OFFSET)
    }

    fn absorb_bytes(&mut self, bytes: &[u8]) {
        self.absorb_usize(bytes.len());
        for byte in bytes {
            self.0 ^= u64::from(*byte);
            self.0 = self.0.wrapping_mul(Self::PRIME);
        }
    }

    fn absorb_u64(&mut self, value: u64) {
        for byte in value.to_le_bytes() {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(Self::PRIME);
        }
    }

    fn absorb_usize(&mut self, value: usize) {
        self.absorb_u64(u64::try_from(value).unwrap_or(u64::MAX));
    }

    fn absorb_f64(&mut self, value: f64) {
        self.absorb_u64(canonical_zero(value).to_bits());
    }

    fn absorb_f64s(&mut self, values: &[f64]) {
        for value in values {
            self.absorb_f64(*value);
        }
    }

    fn absorb_bool(&mut self, value: bool) {
        self.absorb_u64(u64::from(value));
    }

    const fn finish(self) -> u64 {
        self.0
    }
}

fn push_decision_json(output: &mut String, decision: &AssociationDecision) {
    output.push_str("{\"source_subject\":");
    push_json_string(output, &decision.source_subject);
    output.push_str(",\"target_subject\":");
    if let Some(target) = &decision.target_subject {
        push_json_string(output, target);
    } else {
        output.push_str("null");
    }
    output.push_str(",\"association\":");
    push_json_string(output, decision.verdict.as_str());
    output.push_str(",\"change\":");
    push_json_string(output, decision.change.as_str());
    output.push_str(",\"migration\":");
    push_json_string(output, decision.migration.as_str());
    output.push_str(",\"source_fingerprint\":");
    push_surface_json(output, &decision.source_fingerprint);
    output.push_str(",\"target_fingerprint\":");
    if let Some(target) = &decision.target_fingerprint {
        push_surface_json(output, target);
    } else {
        output.push_str("null");
    }
    output.push_str(",\"drift\":");
    if let Some(drift) = decision.drift {
        let _ = write!(
            output,
            "{{\"selector_agrees\":{},\"topology_agrees\":{},\"relative_area\":{},\"centroid_distance\":{},\"orientation\":{},\"relative_extent\":{},\"normalized_score\":{}}}",
            drift.selector_agrees,
            drift.topology_agrees,
            drift.relative_area,
            drift.centroid_distance,
            drift.orientation,
            drift.relative_extent,
            drift.normalized_score
        );
    } else {
        output.push_str("null");
    }
    output.push_str(",\"candidates\":[");
    for (index, candidate) in decision.candidates.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        push_json_string(output, candidate);
    }
    output.push_str("],\"reason\":");
    push_json_string(output, decision.reason);
    output.push('}');
}

fn push_surface_json(output: &mut String, surface: &SurfaceFingerprint) {
    let topology = surface.topology;
    let _ = write!(
        output,
        "{{\"local_fingerprint\":\"{:016x}\",\"face_count\":{},\"surface_area\":{},\"centroid\":[{},{},{}],\"extents\":[{},{},{}],\"orientation_moments\":[{},{},{},{},{},{}],\"topology\":{{\"components\":{},\"boundary_components\":{},\"boundary_loops\":{},\"nonmanifold_edges\":{},\"orientation_conflicts\":{},\"closed_oriented_boundary\":{}}}}}",
        surface.local_fingerprint,
        surface.face_count,
        surface.surface_area,
        surface.centroid[0],
        surface.centroid[1],
        surface.centroid[2],
        surface.extents[0],
        surface.extents[1],
        surface.extents[2],
        surface.orientation_moments[0],
        surface.orientation_moments[1],
        surface.orientation_moments[2],
        surface.orientation_moments[3],
        surface.orientation_moments[4],
        surface.orientation_moments[5],
        topology.component_count,
        topology.boundary_component_count,
        topology.boundary_loop_count,
        topology.nonmanifold_edge_count,
        topology.orientation_conflict_count,
        topology.closed_oriented_boundary
    );
}

fn push_transform_json(output: &mut String, transform: RigidTransform3) {
    let _ = write!(
        output,
        ",\"source_to_target\":{{\"rotation\":[[{},{},{}],[{},{},{}],[{},{},{}]],\"translation\":[{},{},{}]}}",
        transform.rotation[0][0],
        transform.rotation[0][1],
        transform.rotation[0][2],
        transform.rotation[1][0],
        transform.rotation[1][1],
        transform.rotation[1][2],
        transform.rotation[2][0],
        transform.rotation[2][1],
        transform.rotation[2][2],
        transform.translation[0],
        transform.translation[1],
        transform.translation[2]
    );
}

fn push_policy_json(output: &mut String, policy: AssociationPolicy) {
    let _ = write!(
        output,
        ",\"policy\":{{\"stable_relative_area\":{},\"stable_distance\":{},\"stable_orientation\":{},\"stable_relative_extent\":{},\"moved_relative_area\":{},\"moved_distance\":{},\"moved_orientation\":{},\"moved_relative_extent\":{},\"ambiguity_score_gap\":{},\"frame_tolerance\":{},\"max_assignments\":{},\"max_mesh_vertices\":{},\"max_mesh_faces\":{},\"max_face_references\":{},\"max_edge_records\":{},\"max_candidate_tests\":{},\"poll_stride\":{}}}",
        policy.stable_relative_area,
        policy.stable_distance,
        policy.stable_orientation,
        policy.stable_relative_extent,
        policy.moved_relative_area,
        policy.moved_distance,
        policy.moved_orientation,
        policy.moved_relative_extent,
        policy.ambiguity_score_gap,
        policy.frame_tolerance,
        policy.max_assignments,
        policy.max_mesh_vertices,
        policy.max_mesh_faces,
        policy.max_face_references,
        policy.max_edge_records,
        policy.max_candidate_tests,
        policy.poll_stride
    );
}

fn push_json_string(output: &mut String, value: &str) {
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\u{08}' => output.push_str("\\b"),
            '\u{0c}' => output.push_str("\\f"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character <= '\u{1f}' => {
                let _ = write!(output, "\\u{:04x}", u32::from(character));
            }
            character => output.push(character),
        }
    }
    output.push('"');
}

fn markdown_code(value: &str) -> String {
    value.replace('`', "\\`")
}

fn markdown_cell(value: &str) -> String {
    value.replace('|', "\\|").replace('\n', " ")
}

fn checkpoint(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), AssociationRefusal> {
    cx.checkpoint().map_err(|_| {
        refusal(
            "mesh-association-cancelled",
            format!("cancellation observed during {stage} at work item {at}"),
            "retry under a fresh scope from the same immutable inputs",
        )
    })
}

fn poll(
    cx: &Cx<'_>,
    stage: &'static str,
    at: usize,
    stride: usize,
) -> Result<(), AssociationRefusal> {
    if at.is_multiple_of(stride) {
        checkpoint(cx, stage, at)?;
    }
    Ok(())
}

fn resource_refusal(resource: &str, requested: usize, limit: usize) -> AssociationRefusal {
    refusal(
        "mesh-association-resource-bound",
        format!("{resource} requires {requested}; cap is {limit}"),
        "reduce the revision pair or explicitly raise the receipted association limit",
    )
}

fn allocation_refusal(resource: &str, requested: usize) -> AssociationRefusal {
    refusal(
        "mesh-association-allocation",
        format!("could not reserve {requested} {resource}"),
        "reduce the revision pair or provide a larger explicit memory budget",
    )
}

fn refusal(
    code: &'static str,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> AssociationRefusal {
    AssociationRefusal {
        code,
        what: what.into(),
        fix: fix.into(),
    }
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}

fn finite3(value: [f64; 3]) -> bool {
    value.into_iter().all(f64::is_finite)
}

fn add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn scale(value: [f64; 3], factor: f64) -> [f64; 3] {
    [value[0] * factor, value[1] * factor, value[2] * factor]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn norm(value: [f64; 3]) -> f64 {
    dot(value, value).sqrt()
}
