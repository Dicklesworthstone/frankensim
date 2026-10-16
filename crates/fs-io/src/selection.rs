//! Deterministic, mesh-index-free surface assignment for promoted triangle
//! soups.
//!
//! Higher layers supply persistent entity-identity tokens and source-artifact
//! identities. This L2 module deliberately treats both as opaque bytes: it
//! resolves geometry and emits replay evidence, but it does not derive or
//! authenticate L3/L6 identities.

use core::fmt;
use fs_exec::{Cx, ExecMode};
use fs_geom::Point3;
use fs_rep_mesh::Soup;
use std::fmt::Write as _;

/// Versioned meaning of geometric selection and assignment.
pub const MESH_ASSIGNMENT_SEMANTICS_VERSION: &str = "fs-io/mesh-assignment/v1";

/// Maximum owned loop work between explicit cancellation polls.
pub const MESH_ASSIGNMENT_POLL_STRIDE: usize = 4096;

const RECEIPT_NO_CLAIM: &str = "selectors classify the supplied finite tessellation only; caller-supplied source and subject identities are retained but not authenticated; no between-facet, continuum-topology, self-intersection, CAD-semantic, or physical-region-sameness claim is made";

/// Which closed side of a half-space inequality is selected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HalfSpaceSide {
    /// Select triangles whose vertices satisfy `normal · point <= offset`.
    AtMost,
    /// Select triangles whose vertices satisfy `normal · point >= offset`.
    AtLeast,
}

/// A deterministic selector over triangle faces.
///
/// Geometric selectors require every vertex of a triangle to satisfy the
/// predicate. This avoids centroid-only selections that can change merely
/// because a face is subdivided.
#[derive(Debug, Clone, PartialEq)]
pub enum MeshSelector {
    /// Use a format- or adapter-provided named group.
    NamedGroup {
        /// Exact group name.
        name: String,
    },
    /// Select one closed half-space, expanded by `tolerance`.
    HalfSpace {
        /// Plane normal; need not be normalized.
        normal: [f64; 3],
        /// Plane offset in `normal · point = offset`.
        offset: f64,
        /// Selected inequality side.
        side: HalfSpaceSide,
        /// Non-negative tolerance in dot-product units.
        tolerance: f64,
    },
    /// Select triangles fully contained by an axis-aligned box.
    Box {
        /// Inclusive lower corner.
        min: [f64; 3],
        /// Inclusive upper corner.
        max: [f64; 3],
        /// Non-negative coordinate tolerance.
        tolerance: f64,
    },
    /// Select triangles fully contained by a finite cylinder.
    Cylinder {
        /// A point on the cylinder axis at axial coordinate zero.
        origin: [f64; 3],
        /// Nonzero axis direction; normalized during admission.
        axis: [f64; 3],
        /// Non-negative cylinder radius.
        radius: f64,
        /// Inclusive lower axial coordinate.
        axial_min: f64,
        /// Inclusive upper axial coordinate.
        axial_max: f64,
        /// Non-negative spatial tolerance.
        tolerance: f64,
    },
    /// Select every face tied for nearest to a datum within `tolerance`.
    NearestDatum {
        /// Datum point.
        point: [f64; 3],
        /// Maximum admitted point-to-triangle distance.
        max_distance: f64,
        /// Non-negative tie tolerance.
        tolerance: f64,
    },
    /// Escape hatch for explicitly enumerated face ordinals.
    ExplicitFaceSet {
        /// Exact face ordinals in the supplied tessellation.
        faces: Vec<u32>,
        /// Must be true: callers explicitly acknowledge remeshing fragility.
        fragility_acknowledged: bool,
    },
}

/// A name-to-face mapping retained from an importer or external adapter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NamedFaceGroup {
    /// Exact group name.
    pub name: String,
    /// Face ordinals in the supplied tessellation.
    pub faces: Vec<u32>,
}

/// One persistent-identity assignment request.
#[derive(Debug, Clone, PartialEq)]
pub struct AssignmentRequest {
    /// Opaque persistent entity token supplied by the higher layer.
    pub subject: String,
    /// How faces are selected.
    pub selector: MeshSelector,
    /// Whether overlap with other overlap-enabled assignments is intentional.
    pub allow_overlap: bool,
}

/// Explicit resource envelope for one assignment operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssignmentLimits {
    /// Maximum input vertices.
    pub max_mesh_vertices: usize,
    /// Maximum input faces.
    pub max_mesh_faces: usize,
    /// Maximum assignment requests.
    pub max_requests: usize,
    /// Maximum named groups.
    pub max_named_groups: usize,
    /// Maximum aggregate face references across named groups.
    pub max_group_faces: usize,
    /// Maximum aggregate published face references.
    pub max_selected_faces: usize,
    /// Maximum face-predicate evaluations.
    pub max_predicate_tests: u64,
    /// Maximum bytes in any caller-supplied identity or name.
    pub max_label_bytes: usize,
}

impl AssignmentLimits {
    /// Production default aligned with the bounded STEP tessellation handoff.
    pub const DEFAULT: Self = Self {
        max_mesh_vertices: 1_000_000,
        max_mesh_faces: 1_000_000,
        max_requests: 4096,
        max_named_groups: 4096,
        max_group_faces: 4_000_000,
        max_selected_faces: 4_000_000,
        max_predicate_tests: 100_000_000,
        max_label_bytes: 4096,
    };
}

impl Default for AssignmentLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Geometric sanity statistics for one resolved assignment.
#[derive(Debug, Clone, PartialEq)]
pub struct AssignmentStats {
    /// Number of selected faces.
    pub face_count: usize,
    /// Sum of selected triangle areas, in squared `length_unit`.
    pub surface_area: f64,
    /// Enclosed oriented volume when the selected subset is a closed,
    /// consistently oriented triangle boundary; otherwise `None`.
    pub enclosed_volume: Option<f64>,
    /// Inclusive lower bound of selected vertices.
    pub bounds_min: [f64; 3],
    /// Inclusive upper bound of selected vertices.
    pub bounds_max: [f64; 3],
}

/// One successfully resolved persistent subject.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedAssignment {
    /// Opaque persistent entity token supplied by the caller.
    pub subject: String,
    /// Selected face ordinals, sorted and unique.
    pub faces: Vec<u32>,
    /// Geometric sanity statistics.
    pub stats: AssignmentStats,
    /// Versioned local replay fingerprint of the selector semantics.
    pub selector_fingerprint: u64,
    /// Whether overlap was explicitly admitted for this assignment.
    pub allow_overlap: bool,
}

/// Success-only receipt for one complete assignment operation.
#[derive(Debug, Clone, PartialEq)]
pub struct AssignmentReceipt {
    source_identity: String,
    length_unit: String,
    source_mesh_fingerprint: u64,
    named_groups_fingerprint: u64,
    requests_fingerprint: u64,
    assignments_fingerprint: u64,
}

impl AssignmentReceipt {
    /// Caller-supplied source-artifact identity hook.
    #[must_use]
    pub fn source_identity(&self) -> &str {
        &self.source_identity
    }

    /// Caller-supplied length-unit identifier.
    #[must_use]
    pub fn length_unit(&self) -> &str {
        &self.length_unit
    }

    /// Local exact-bit replay fingerprint of the supplied soup.
    #[must_use]
    pub const fn source_mesh_fingerprint(&self) -> u64 {
        self.source_mesh_fingerprint
    }

    /// Local replay fingerprint of named-group mappings.
    #[must_use]
    pub const fn named_groups_fingerprint(&self) -> u64 {
        self.named_groups_fingerprint
    }

    /// Local replay fingerprint of ordered requests.
    #[must_use]
    pub const fn requests_fingerprint(&self) -> u64 {
        self.requests_fingerprint
    }

    /// Local replay fingerprint of the published assignments.
    #[must_use]
    pub const fn assignments_fingerprint(&self) -> u64 {
        self.assignments_fingerprint
    }

    fn to_json(&self, assignments: &[ResolvedAssignment]) -> String {
        let mut output = String::from("{\"kind\":\"mesh-assignment-receipt\",\"version\":");
        push_json_string(&mut output, MESH_ASSIGNMENT_SEMANTICS_VERSION);
        output.push_str(",\"source_identity\":");
        push_json_string(&mut output, &self.source_identity);
        output.push_str(",\"length_unit\":");
        push_json_string(&mut output, &self.length_unit);
        let _ = write!(
            output,
            ",\"source_mesh_fingerprint\":\"{:016x}\",\"named_groups_fingerprint\":\"{:016x}\",\"requests_fingerprint\":\"{:016x}\",\"assignments_fingerprint\":\"{:016x}\",\"assignments\":[",
            self.source_mesh_fingerprint,
            self.named_groups_fingerprint,
            self.requests_fingerprint,
            self.assignments_fingerprint
        );
        for (index, assignment) in assignments.iter().enumerate() {
            if index != 0 {
                output.push(',');
            }
            output.push_str("{\"subject\":");
            push_json_string(&mut output, &assignment.subject);
            let _ = write!(
                output,
                ",\"selector_fingerprint\":\"{:016x}\",\"face_count\":{},\"surface_area\":{},\"enclosed_volume\":",
                assignment.selector_fingerprint,
                assignment.stats.face_count,
                assignment.stats.surface_area
            );
            match assignment.stats.enclosed_volume {
                Some(volume) => {
                    let _ = write!(output, "{volume}");
                }
                None => output.push_str("null"),
            }
            let _ = write!(
                output,
                ",\"bounds_min\":[{},{},{}],\"bounds_max\":[{},{},{}],\"allow_overlap\":{}}}",
                assignment.stats.bounds_min[0],
                assignment.stats.bounds_min[1],
                assignment.stats.bounds_min[2],
                assignment.stats.bounds_max[0],
                assignment.stats.bounds_max[1],
                assignment.stats.bounds_max[2],
                assignment.allow_overlap
            );
        }
        output.push_str("],\"authority\":\"finite-tessellation-selection\",\"no_claim\":");
        push_json_string(&mut output, RECEIPT_NO_CLAIM);
        output.push('}');
        output
    }
}

/// Atomic success value: no partial assignments or success receipt escape a
/// refusal.
#[derive(Debug, Clone, PartialEq)]
pub struct AssignmentReport {
    /// Assignments in request order.
    pub assignments: Vec<ResolvedAssignment>,
    /// Receipt binding the complete operation.
    pub receipt: AssignmentReceipt,
}

impl AssignmentReport {
    /// Canonical, one-line JSON binding this report's receipt to its
    /// assignments.
    #[must_use]
    pub fn to_json(&self) -> String {
        self.receipt.to_json(&self.assignments)
    }
}

/// Structured, actionable assignment refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssignmentRefusal {
    /// Stable machine-facing refusal code.
    pub code: &'static str,
    /// Specific diagnosis.
    pub what: String,
    /// Actionable correction.
    pub fix: String,
}

impl fmt::Display for AssignmentRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}: {} (fix: {})",
            self.code, self.what, self.fix
        )
    }
}

impl std::error::Error for AssignmentRefusal {}

/// Resolve named, geometric, or explicitly acknowledged face selections.
///
/// # Errors
///
/// Returns [`AssignmentRefusal`] for malformed geometry or selector inputs,
/// resource-bound violations, cancellation, empty selections, undeclared
/// named groups, fragile face sets without acknowledgement, and unintended
/// overlap. Publication is atomic.
#[allow(clippy::too_many_arguments)]
pub fn resolve_mesh_assignments(
    soup: &Soup,
    source_identity: &str,
    length_unit: &str,
    named_groups: &[NamedFaceGroup],
    requests: &[AssignmentRequest],
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> Result<AssignmentReport, AssignmentRefusal> {
    checkpoint(cx, "assignment-entry", 0)?;
    if cx.mode() != ExecMode::Deterministic {
        return Err(refusal(
            "mesh-assignment-fast-mode",
            "mesh assignment receipts require deterministic execution mode",
            "retry with ExecMode::Deterministic",
        ));
    }
    validate_limits(limits)?;
    validate_label(source_identity, "source identity", limits.max_label_bytes)?;
    validate_label(length_unit, "length unit", limits.max_label_bytes)?;
    validate_soup(soup, limits, cx)?;
    if requests.is_empty() {
        return Err(refusal(
            "mesh-assignment-empty-plan",
            "the assignment plan contains no requests",
            "declare at least one persistent subject and selector",
        ));
    }
    if requests.len() > limits.max_requests {
        return Err(resource_refusal(
            "assignment requests",
            requests.len(),
            limits.max_requests,
        ));
    }
    if named_groups.len() > limits.max_named_groups {
        return Err(resource_refusal(
            "named groups",
            named_groups.len(),
            limits.max_named_groups,
        ));
    }

    let group_order = validate_groups(soup, named_groups, limits, cx)?;
    validate_requests(requests, limits, cx)?;
    preflight_predicate_work(soup, requests, limits)?;

    let source_mesh_fingerprint = fingerprint_soup(soup, cx)?;
    let named_groups_fingerprint = fingerprint_groups(named_groups, cx)?;
    let requests_fingerprint = fingerprint_requests(requests, cx)?;

    let mut assignments = Vec::new();
    assignments
        .try_reserve_exact(requests.len())
        .map_err(|_| allocation_refusal("assignment rows", requests.len()))?;
    let mut total_selected = 0usize;
    for (request_index, request) in requests.iter().enumerate() {
        checkpoint(cx, "assignment-request", request_index)?;
        let mut faces = select_faces(
            soup,
            named_groups,
            &group_order,
            &request.selector,
            limits,
            cx,
        )?;
        faces.sort_unstable();
        if faces.windows(2).any(|window| window[0] == window[1]) {
            return Err(refusal(
                "mesh-assignment-duplicate-face",
                format!(
                    "subject {:?} selected the same face more than once",
                    request.subject
                ),
                "deduplicate the named or explicit face set at its source",
            ));
        }
        if faces.is_empty() {
            return Err(refusal(
                "mesh-assignment-empty-selection",
                format!(
                    "selector for subject {:?} resolved to zero faces",
                    request.subject
                ),
                "correct the selector, units, tolerance, source group, or geometry artifact",
            ));
        }
        total_selected = total_selected.checked_add(faces.len()).ok_or_else(|| {
            resource_refusal(
                "aggregate selected faces",
                usize::MAX,
                limits.max_selected_faces,
            )
        })?;
        if total_selected > limits.max_selected_faces {
            return Err(resource_refusal(
                "aggregate selected faces",
                total_selected,
                limits.max_selected_faces,
            ));
        }
        let stats = assignment_stats(soup, &faces, cx)?;
        assignments.push(ResolvedAssignment {
            subject: request.subject.clone(),
            selector_fingerprint: fingerprint_selector(&request.selector, cx)?,
            allow_overlap: request.allow_overlap,
            faces,
            stats,
        });
    }
    validate_overlap(soup.triangles.len(), &assignments, cx)?;
    checkpoint(cx, "assignment-publication", assignments.len())?;
    let assignments_fingerprint = fingerprint_assignments(&assignments, cx)?;
    Ok(AssignmentReport {
        receipt: AssignmentReceipt {
            source_identity: source_identity.to_string(),
            length_unit: length_unit.to_string(),
            source_mesh_fingerprint,
            named_groups_fingerprint,
            requests_fingerprint,
            assignments_fingerprint,
        },
        assignments,
    })
}

fn validate_limits(limits: AssignmentLimits) -> Result<(), AssignmentRefusal> {
    for (name, value) in [
        ("max_mesh_vertices", limits.max_mesh_vertices),
        ("max_mesh_faces", limits.max_mesh_faces),
        ("max_requests", limits.max_requests),
        ("max_named_groups", limits.max_named_groups),
        ("max_group_faces", limits.max_group_faces),
        ("max_selected_faces", limits.max_selected_faces),
        ("max_label_bytes", limits.max_label_bytes),
    ] {
        if value == 0 {
            return Err(refusal(
                "mesh-assignment-invalid-limits",
                format!("{name} must be nonzero"),
                "provide a nonzero explicit assignment resource envelope",
            ));
        }
    }
    if limits.max_predicate_tests == 0 {
        return Err(refusal(
            "mesh-assignment-invalid-limits",
            "max_predicate_tests must be nonzero",
            "provide a nonzero explicit assignment resource envelope",
        ));
    }
    Ok(())
}

fn validate_label(value: &str, label: &str, max_bytes: usize) -> Result<(), AssignmentRefusal> {
    if value.is_empty() || value.trim() != value || value.chars().any(char::is_control) {
        return Err(refusal(
            "mesh-assignment-invalid-label",
            format!("{label} must be nonempty, trim-canonical, and control-free"),
            format!("supply a canonical {label}"),
        ));
    }
    if value.len() > max_bytes {
        return Err(resource_refusal(label, value.len(), max_bytes));
    }
    Ok(())
}

fn validate_soup(
    soup: &Soup,
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> Result<(), AssignmentRefusal> {
    if soup.positions.is_empty() || soup.triangles.is_empty() {
        return Err(refusal(
            "mesh-assignment-empty-geometry",
            "the supplied mesh must contain vertices and faces",
            "run the import quarantine and promotion path before assignment",
        ));
    }
    if soup.positions.len() > limits.max_mesh_vertices {
        return Err(resource_refusal(
            "mesh vertices",
            soup.positions.len(),
            limits.max_mesh_vertices,
        ));
    }
    if soup.triangles.len() > limits.max_mesh_faces {
        return Err(resource_refusal(
            "mesh faces",
            soup.triangles.len(),
            limits.max_mesh_faces,
        ));
    }
    for (index, point) in soup.positions.iter().enumerate() {
        poll(cx, "assignment-vertices", index)?;
        if !finite3([point.x, point.y, point.z]) {
            return Err(refusal(
                "mesh-assignment-nonfinite-vertex",
                format!("vertex {index} has non-finite coordinates"),
                "refuse or repair the geometry in quarantine before assignment",
            ));
        }
    }
    for (face, triangle) in soup.triangles.iter().enumerate() {
        poll(cx, "assignment-faces", face)?;
        let mut points = [Point3::new(0.0, 0.0, 0.0); 3];
        for (corner, vertex) in triangle.iter().copied().enumerate() {
            let vertex = usize::try_from(vertex).map_err(|_| {
                refusal(
                    "mesh-assignment-missing-vertex",
                    format!("face {face} references unrepresentable vertex {vertex}"),
                    "repair the importer so every face references an admitted vertex",
                )
            })?;
            let Some(point) = soup.positions.get(vertex).copied() else {
                return Err(refusal(
                    "mesh-assignment-missing-vertex",
                    format!("face {face} references missing vertex {vertex}"),
                    "repair the importer so every face references an admitted vertex",
                ));
            };
            points[corner] = point;
        }
        let twice_area = norm(cross(
            delta(points[1], points[0]),
            delta(points[2], points[0]),
        ));
        if !(twice_area.is_finite() && twice_area > 0.0) {
            return Err(refusal(
                "mesh-assignment-degenerate-face",
                format!("face {face} is degenerate or has non-finite area"),
                "run repair and promotion before resolving assignments",
            ));
        }
    }
    Ok(())
}

fn validate_groups(
    soup: &Soup,
    groups: &[NamedFaceGroup],
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> Result<Vec<usize>, AssignmentRefusal> {
    let mut total_faces = 0usize;
    let mut order = Vec::new();
    order
        .try_reserve_exact(groups.len())
        .map_err(|_| allocation_refusal("named-group index", groups.len()))?;
    for (index, group) in groups.iter().enumerate() {
        poll(cx, "assignment-groups", index)?;
        validate_label(&group.name, "named group", limits.max_label_bytes)?;
        if group.faces.is_empty() {
            return Err(refusal(
                "mesh-assignment-empty-group",
                format!("named group {:?} contains no faces", group.name),
                "remove the group or populate it from the importer",
            ));
        }
        total_faces = total_faces.checked_add(group.faces.len()).ok_or_else(|| {
            resource_refusal("named-group faces", usize::MAX, limits.max_group_faces)
        })?;
        if total_faces > limits.max_group_faces {
            return Err(resource_refusal(
                "named-group faces",
                total_faces,
                limits.max_group_faces,
            ));
        }
        let mut sorted = clone_faces(&group.faces)?;
        sorted.sort_unstable();
        if sorted.windows(2).any(|window| window[0] == window[1]) {
            return Err(refusal(
                "mesh-assignment-duplicate-group-face",
                format!("named group {:?} repeats a face", group.name),
                "deduplicate the importer-provided named group",
            ));
        }
        if let Some(face) = sorted
            .iter()
            .copied()
            .find(|face| usize::try_from(*face).map_or(true, |face| face >= soup.triangles.len()))
        {
            return Err(refusal(
                "mesh-assignment-group-face-range",
                format!(
                    "named group {:?} references missing face {face}",
                    group.name
                ),
                "repair the importer-provided group mapping",
            ));
        }
        order.push(index);
    }
    order.sort_unstable_by(|left, right| groups[*left].name.cmp(&groups[*right].name));
    for duplicate in order.windows(2) {
        if groups[duplicate[0]].name == groups[duplicate[1]].name {
            return Err(refusal(
                "mesh-assignment-duplicate-group",
                format!(
                    "named group {:?} is declared more than once",
                    groups[duplicate[0]].name
                ),
                "retain exactly one mapping for each exact group name",
            ));
        }
    }
    Ok(order)
}

fn validate_requests(
    requests: &[AssignmentRequest],
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> Result<(), AssignmentRefusal> {
    let mut subjects = Vec::new();
    subjects
        .try_reserve_exact(requests.len())
        .map_err(|_| allocation_refusal("subject index", requests.len()))?;
    for (index, request) in requests.iter().enumerate() {
        poll(cx, "assignment-request-admission", index)?;
        validate_label(
            &request.subject,
            "persistent subject token",
            limits.max_label_bytes,
        )?;
        validate_selector(&request.selector, limits)?;
        subjects.push((request.subject.as_str(), index));
    }
    subjects.sort_unstable_by(|left, right| left.0.cmp(right.0));
    if let Some(duplicate) = subjects
        .windows(2)
        .find(|window| window[0].0 == window[1].0)
    {
        return Err(refusal(
            "mesh-assignment-duplicate-subject",
            format!(
                "persistent subject {:?} has multiple assignment requests",
                duplicate[0].0
            ),
            "compose one selector per persistent subject",
        ));
    }
    Ok(())
}

fn validate_selector(
    selector: &MeshSelector,
    limits: AssignmentLimits,
) -> Result<(), AssignmentRefusal> {
    match selector {
        MeshSelector::NamedGroup { name } => {
            validate_label(name, "named-group selector", limits.max_label_bytes)
        }
        MeshSelector::HalfSpace {
            normal,
            offset,
            tolerance,
            ..
        } => {
            validate_finite_vector(*normal, "half-space normal")?;
            validate_direction(*normal, "half-space normal")?;
            validate_finite(*offset, "half-space offset")?;
            validate_nonnegative(*tolerance, "half-space tolerance")?;
            validate_finite(*offset - *tolerance, "half-space lower threshold")?;
            validate_finite(*offset + *tolerance, "half-space upper threshold")
        }
        MeshSelector::Box {
            min,
            max,
            tolerance,
        } => {
            validate_finite_vector(*min, "box minimum")?;
            validate_finite_vector(*max, "box maximum")?;
            if (0..3).any(|axis| min[axis] > max[axis]) {
                return Err(selector_refusal(
                    "box minimum must not exceed maximum on any axis",
                ));
            }
            validate_nonnegative(*tolerance, "box tolerance")?;
            for axis in 0..3 {
                validate_finite(min[axis] - *tolerance, "box admitted minimum")?;
                validate_finite(max[axis] + *tolerance, "box admitted maximum")?;
            }
            Ok(())
        }
        MeshSelector::Cylinder {
            origin,
            axis,
            radius,
            axial_min,
            axial_max,
            tolerance,
        } => {
            validate_finite_vector(*origin, "cylinder origin")?;
            validate_finite_vector(*axis, "cylinder axis")?;
            validate_direction(*axis, "cylinder axis")?;
            validate_nonnegative(*radius, "cylinder radius")?;
            validate_finite(*axial_min, "cylinder axial minimum")?;
            validate_finite(*axial_max, "cylinder axial maximum")?;
            if axial_min > axial_max {
                return Err(selector_refusal(
                    "cylinder axial minimum must not exceed maximum",
                ));
            }
            validate_nonnegative(*tolerance, "cylinder tolerance")?;
            validate_finite(*radius + *tolerance, "cylinder admitted radius")?;
            validate_finite(*axial_min - *tolerance, "cylinder admitted axial minimum")?;
            validate_finite(*axial_max + *tolerance, "cylinder admitted axial maximum")
        }
        MeshSelector::NearestDatum {
            point,
            max_distance,
            tolerance,
        } => {
            validate_finite_vector(*point, "datum point")?;
            validate_nonnegative(*max_distance, "datum maximum distance")?;
            validate_nonnegative(*tolerance, "datum tie tolerance")
        }
        MeshSelector::ExplicitFaceSet {
            faces,
            fragility_acknowledged,
        } => {
            if !fragility_acknowledged {
                return Err(refusal(
                    "mesh-assignment-fragility-unacknowledged",
                    "an explicit face-set selector did not acknowledge remeshing fragility",
                    "set fragility_acknowledged=true or use a named/geometric selector",
                ));
            }
            if faces.is_empty() {
                return Err(selector_refusal(
                    "an explicit face-set selector must contain at least one face",
                ));
            }
            if faces.len() > limits.max_selected_faces {
                return Err(resource_refusal(
                    "explicit face-set faces",
                    faces.len(),
                    limits.max_selected_faces,
                ));
            }
            Ok(())
        }
    }
}

fn preflight_predicate_work(
    soup: &Soup,
    requests: &[AssignmentRequest],
    limits: AssignmentLimits,
) -> Result<(), AssignmentRefusal> {
    let geometric = requests
        .iter()
        .filter(|request| {
            matches!(
                &request.selector,
                MeshSelector::HalfSpace { .. }
                    | MeshSelector::Box { .. }
                    | MeshSelector::Cylinder { .. }
                    | MeshSelector::NearestDatum { .. }
            )
        })
        .count();
    let tests = u64::try_from(geometric)
        .ok()
        .and_then(|geometric| {
            u64::try_from(soup.triangles.len())
                .ok()
                .and_then(|faces| geometric.checked_mul(faces))
        })
        .ok_or_else(|| {
            refusal(
                "mesh-assignment-work-overflow",
                "face-predicate work overflowed",
                "reduce the mesh or assignment plan",
            )
        })?;
    if tests > limits.max_predicate_tests {
        return Err(refusal(
            "mesh-assignment-work-limit",
            format!(
                "assignment requires {tests} face-predicate tests; cap is {}",
                limits.max_predicate_tests
            ),
            "split the plan, reduce the mesh, or explicitly raise max_predicate_tests",
        ));
    }
    Ok(())
}

fn select_faces(
    soup: &Soup,
    groups: &[NamedFaceGroup],
    group_order: &[usize],
    selector: &MeshSelector,
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> Result<Vec<u32>, AssignmentRefusal> {
    match selector {
        MeshSelector::NamedGroup { name } => {
            let found =
                group_order.binary_search_by(|index| groups[*index].name.as_str().cmp(name));
            let Ok(found) = found else {
                return Err(refusal(
                    "mesh-assignment-group-unknown",
                    format!("named group {name:?} was not supplied by the importer"),
                    "choose a retained group name or use a geometric selector",
                ));
            };
            clone_faces(&groups[group_order[found]].faces)
        }
        MeshSelector::ExplicitFaceSet { faces, .. } => {
            let selected = clone_faces(faces)?;
            if let Some(face) = selected.iter().copied().find(|face| {
                usize::try_from(*face).map_or(true, |face| face >= soup.triangles.len())
            }) {
                return Err(refusal(
                    "mesh-assignment-face-range",
                    format!("explicit face set references missing face {face}"),
                    "update the fragile face set for this exact tessellation",
                ));
            }
            Ok(selected)
        }
        MeshSelector::NearestDatum {
            point,
            max_distance,
            tolerance,
        } => select_nearest(
            soup,
            Point3::new(point[0], point[1], point[2]),
            *max_distance,
            *tolerance,
            limits,
            cx,
        ),
        _ => {
            let mut selected = Vec::new();
            selected
                .try_reserve_exact(soup.triangles.len())
                .map_err(|_| allocation_refusal("selected faces", soup.triangles.len()))?;
            for face in 0..soup.triangles.len() {
                poll(cx, "assignment-predicate", face)?;
                let triangle = soup_triangle(soup, face);
                if triangle_matches(&triangle, selector) {
                    selected.push(u32::try_from(face).map_err(|_| {
                        refusal(
                            "mesh-assignment-face-ordinal",
                            format!("face ordinal {face} does not fit u32"),
                            "reduce the mesh below the admitted face-index representation",
                        )
                    })?);
                }
            }
            Ok(selected)
        }
    }
}

fn select_nearest(
    soup: &Soup,
    point: Point3,
    max_distance: f64,
    tolerance: f64,
    limits: AssignmentLimits,
    cx: &Cx<'_>,
) -> Result<Vec<u32>, AssignmentRefusal> {
    let mut distances = Vec::new();
    distances
        .try_reserve_exact(soup.triangles.len())
        .map_err(|_| allocation_refusal("nearest-datum distances", soup.triangles.len()))?;
    let mut nearest = f64::INFINITY;
    for face in 0..soup.triangles.len() {
        poll(cx, "assignment-nearest", face)?;
        let [a, b, c] = soup_triangle(soup, face);
        let distance = point_triangle_distance(point, a, b, c);
        if !distance.is_finite() {
            return Err(refusal(
                "mesh-assignment-distance-overflow",
                format!("datum distance overflowed at face {face}"),
                "rescale coordinates or repair the geometry",
            ));
        }
        nearest = nearest.min(distance);
        distances.push(distance);
    }
    if nearest > max_distance {
        return Ok(Vec::new());
    }
    let threshold = nearest + tolerance;
    if !threshold.is_finite() {
        return Err(selector_refusal(
            "datum nearest-distance plus tolerance must stay finite",
        ));
    }
    let count = distances
        .iter()
        .filter(|distance| **distance <= threshold && **distance <= max_distance)
        .count();
    if count > limits.max_selected_faces {
        return Err(resource_refusal(
            "nearest-datum selected faces",
            count,
            limits.max_selected_faces,
        ));
    }
    let mut selected = Vec::new();
    selected
        .try_reserve_exact(count)
        .map_err(|_| allocation_refusal("nearest-datum selected faces", count))?;
    for (face, distance) in distances.into_iter().enumerate() {
        if distance <= threshold && distance <= max_distance {
            selected.push(u32::try_from(face).map_err(|_| {
                refusal(
                    "mesh-assignment-face-ordinal",
                    format!("face ordinal {face} does not fit u32"),
                    "reduce the mesh below the admitted face-index representation",
                )
            })?);
        }
    }
    Ok(selected)
}

fn triangle_matches(triangle: &[Point3; 3], selector: &MeshSelector) -> bool {
    match selector {
        MeshSelector::HalfSpace {
            normal,
            offset,
            side,
            tolerance,
        } => triangle.iter().all(|point| {
            let projection = dot(*normal, [point.x, point.y, point.z]);
            match side {
                HalfSpaceSide::AtMost => projection <= offset + tolerance,
                HalfSpaceSide::AtLeast => projection >= offset - tolerance,
            }
        }),
        MeshSelector::Box {
            min,
            max,
            tolerance,
        } => triangle.iter().all(|point| {
            let point = [point.x, point.y, point.z];
            (0..3).all(|axis| {
                point[axis] >= min[axis] - tolerance && point[axis] <= max[axis] + tolerance
            })
        }),
        MeshSelector::Cylinder {
            origin,
            axis,
            radius,
            axial_min,
            axial_max,
            tolerance,
        } => {
            let axis_norm = norm(*axis);
            let unit = [
                axis[0] / axis_norm,
                axis[1] / axis_norm,
                axis[2] / axis_norm,
            ];
            let admitted_radius = radius + tolerance;
            triangle.iter().all(|point| {
                let relative = [
                    point.x - origin[0],
                    point.y - origin[1],
                    point.z - origin[2],
                ];
                let axial = dot(relative, unit);
                let radial = [
                    relative[0] - axial * unit[0],
                    relative[1] - axial * unit[1],
                    relative[2] - axial * unit[2],
                ];
                axial >= axial_min - tolerance
                    && axial <= axial_max + tolerance
                    && norm(radial) <= admitted_radius
            })
        }
        MeshSelector::NamedGroup { .. }
        | MeshSelector::NearestDatum { .. }
        | MeshSelector::ExplicitFaceSet { .. } => false,
    }
}

fn assignment_stats(
    soup: &Soup,
    faces: &[u32],
    cx: &Cx<'_>,
) -> Result<AssignmentStats, AssignmentRefusal> {
    let mut surface_area = 0.0;
    let mut signed_volume = 0.0;
    let mut bounds_min = [f64::INFINITY; 3];
    let mut bounds_max = [f64::NEG_INFINITY; 3];
    let edge_count = faces.len().checked_mul(3).ok_or_else(|| {
        refusal(
            "mesh-assignment-edge-overflow",
            "selected edge count overflowed",
            "reduce the selected face set",
        )
    })?;
    let mut edges = Vec::new();
    edges
        .try_reserve_exact(edge_count)
        .map_err(|_| allocation_refusal("selected oriented edges", edge_count))?;
    for (position, face) in faces.iter().copied().enumerate() {
        poll(cx, "assignment-statistics", position)?;
        let face_index = usize::try_from(face).map_err(|_| {
            refusal(
                "mesh-assignment-face-ordinal",
                format!("face ordinal {face} is not representable"),
                "repair the selector output",
            )
        })?;
        let triangle_indices = soup.triangles[face_index];
        let [a, b, c] = soup_triangle(soup, face_index);
        let twice_area = norm(cross(delta(b, a), delta(c, a)));
        surface_area += 0.5 * twice_area;
        signed_volume += dot([a.x, a.y, a.z], cross([b.x, b.y, b.z], [c.x, c.y, c.z])) / 6.0;
        for point in [a, b, c] {
            let point = [point.x, point.y, point.z];
            for axis in 0..3 {
                bounds_min[axis] = bounds_min[axis].min(point[axis]);
                bounds_max[axis] = bounds_max[axis].max(point[axis]);
            }
        }
        for (from, to) in [
            (triangle_indices[0], triangle_indices[1]),
            (triangle_indices[1], triangle_indices[2]),
            (triangle_indices[2], triangle_indices[0]),
        ] {
            let key = (from.min(to), from.max(to));
            let orientation = if from <= to { 1i8 } else { -1i8 };
            edges.push((key.0, key.1, orientation));
        }
    }
    if !surface_area.is_finite()
        || !signed_volume.is_finite()
        || !finite3(bounds_min)
        || !finite3(bounds_max)
    {
        return Err(refusal(
            "mesh-assignment-statistics-overflow",
            "assignment area, volume, or bounds overflowed",
            "rescale the geometry or split the assignment",
        ));
    }
    edges.sort_unstable();
    let mut closed = true;
    let mut start = 0usize;
    while start < edges.len() {
        let key = (edges[start].0, edges[start].1);
        let mut end = start;
        let mut balance = 0i32;
        while end < edges.len() && (edges[end].0, edges[end].1) == key {
            balance += i32::from(edges[end].2);
            end += 1;
        }
        if end - start != 2 || balance != 0 {
            closed = false;
            break;
        }
        start = end;
    }
    Ok(AssignmentStats {
        face_count: faces.len(),
        surface_area,
        enclosed_volume: closed.then_some(signed_volume.abs()),
        bounds_min,
        bounds_max,
    })
}

fn validate_overlap(
    face_count: usize,
    assignments: &[ResolvedAssignment],
    cx: &Cx<'_>,
) -> Result<(), AssignmentRefusal> {
    let mut owner = Vec::new();
    owner
        .try_reserve_exact(face_count)
        .map_err(|_| allocation_refusal("face ownership table", face_count))?;
    owner.resize(face_count, None::<usize>);
    let mut visited = 0usize;
    for (assignment_index, assignment) in assignments.iter().enumerate() {
        for face in &assignment.faces {
            poll(cx, "assignment-overlap", visited)?;
            visited += 1;
            let face = usize::try_from(*face).map_err(|_| {
                refusal(
                    "mesh-assignment-face-ordinal",
                    "selected face ordinal is not representable",
                    "repair the selector output",
                )
            })?;
            if let Some(previous) = owner[face] {
                if !(assignments[previous].allow_overlap && assignment.allow_overlap) {
                    return Err(refusal(
                        "mesh-assignment-overlap",
                        format!(
                            "face {face} is assigned to both {:?} and {:?}",
                            assignments[previous].subject, assignment.subject
                        ),
                        "make selectors disjoint or explicitly enable overlap on every participating assignment",
                    ));
                }
            } else {
                owner[face] = Some(assignment_index);
            }
        }
    }
    Ok(())
}

fn point_triangle_distance(point: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    // Real-Time Collision Detection, Christer Ericson, closest-point regions.
    let ab = delta(b, a);
    let ac = delta(c, a);
    let ap = delta(point, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return norm(ap);
    }
    let bp = delta(point, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return norm(bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return norm(sub(ap, scale(ab, v)));
    }
    let cp = delta(point, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return norm(cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return norm(sub(ap, scale(ac, w)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let edge = delta(c, b);
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return norm(sub(bp, scale(edge, w)));
    }
    let denominator = 1.0 / (va + vb + vc);
    let v = vb * denominator;
    let w = vc * denominator;
    norm(sub(ap, add(scale(ab, v), scale(ac, w))))
}

pub(crate) fn fingerprint_soup(soup: &Soup, cx: &Cx<'_>) -> Result<u64, AssignmentRefusal> {
    let mut hash = Fingerprint::new();
    hash.bytes(MESH_ASSIGNMENT_SEMANTICS_VERSION.as_bytes());
    hash.usize(soup.positions.len());
    for (index, point) in soup.positions.iter().enumerate() {
        poll(cx, "assignment-source-fingerprint", index)?;
        hash.u64(point.x.to_bits());
        hash.u64(point.y.to_bits());
        hash.u64(point.z.to_bits());
    }
    hash.usize(soup.triangles.len());
    for (index, triangle) in soup.triangles.iter().enumerate() {
        poll(
            cx,
            "assignment-source-fingerprint",
            soup.positions.len() + index,
        )?;
        for vertex in triangle {
            hash.u32(*vertex);
        }
    }
    Ok(hash.finish())
}

fn fingerprint_groups(groups: &[NamedFaceGroup], cx: &Cx<'_>) -> Result<u64, AssignmentRefusal> {
    let mut hash = Fingerprint::new();
    hash.bytes(MESH_ASSIGNMENT_SEMANTICS_VERSION.as_bytes());
    hash.usize(groups.len());
    let mut visited = 0usize;
    for group in groups {
        poll(cx, "assignment-group-fingerprint", visited)?;
        visited += 1;
        hash.framed(group.name.as_bytes());
        hash.usize(group.faces.len());
        for face in &group.faces {
            poll(cx, "assignment-group-fingerprint", visited)?;
            visited += 1;
            hash.u32(*face);
        }
    }
    Ok(hash.finish())
}

fn fingerprint_requests(
    requests: &[AssignmentRequest],
    cx: &Cx<'_>,
) -> Result<u64, AssignmentRefusal> {
    let mut hash = Fingerprint::new();
    hash.bytes(MESH_ASSIGNMENT_SEMANTICS_VERSION.as_bytes());
    hash.usize(requests.len());
    for (index, request) in requests.iter().enumerate() {
        poll(cx, "assignment-request-fingerprint", index)?;
        hash.framed(request.subject.as_bytes());
        hash.u8(u8::from(request.allow_overlap));
        absorb_selector(
            &mut hash,
            &request.selector,
            cx,
            "assignment-request-fingerprint",
        )?;
    }
    Ok(hash.finish())
}

fn fingerprint_selector(selector: &MeshSelector, cx: &Cx<'_>) -> Result<u64, AssignmentRefusal> {
    let mut hash = Fingerprint::new();
    hash.bytes(MESH_ASSIGNMENT_SEMANTICS_VERSION.as_bytes());
    absorb_selector(&mut hash, selector, cx, "assignment-selector-fingerprint")?;
    Ok(hash.finish())
}

fn absorb_selector(
    hash: &mut Fingerprint,
    selector: &MeshSelector,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<(), AssignmentRefusal> {
    match selector {
        MeshSelector::NamedGroup { name } => {
            hash.u8(1);
            hash.framed(name.as_bytes());
        }
        MeshSelector::HalfSpace {
            normal,
            offset,
            side,
            tolerance,
        } => {
            hash.u8(2);
            for value in normal {
                hash.u64(value.to_bits());
            }
            hash.u64(offset.to_bits());
            hash.u8(match side {
                HalfSpaceSide::AtMost => 0,
                HalfSpaceSide::AtLeast => 1,
            });
            hash.u64(tolerance.to_bits());
        }
        MeshSelector::Box {
            min,
            max,
            tolerance,
        } => {
            hash.u8(3);
            for value in min.iter().chain(max) {
                hash.u64(value.to_bits());
            }
            hash.u64(tolerance.to_bits());
        }
        MeshSelector::Cylinder {
            origin,
            axis,
            radius,
            axial_min,
            axial_max,
            tolerance,
        } => {
            hash.u8(4);
            for value in origin.iter().chain(axis) {
                hash.u64(value.to_bits());
            }
            for value in [radius, axial_min, axial_max, tolerance] {
                hash.u64(value.to_bits());
            }
        }
        MeshSelector::NearestDatum {
            point,
            max_distance,
            tolerance,
        } => {
            hash.u8(5);
            for value in point {
                hash.u64(value.to_bits());
            }
            hash.u64(max_distance.to_bits());
            hash.u64(tolerance.to_bits());
        }
        MeshSelector::ExplicitFaceSet {
            faces,
            fragility_acknowledged,
        } => {
            hash.u8(6);
            hash.u8(u8::from(*fragility_acknowledged));
            hash.usize(faces.len());
            for (index, face) in faces.iter().enumerate() {
                poll(cx, stage, index)?;
                hash.u32(*face);
            }
        }
    }
    Ok(())
}

pub(crate) fn fingerprint_assignments(
    assignments: &[ResolvedAssignment],
    cx: &Cx<'_>,
) -> Result<u64, AssignmentRefusal> {
    let mut hash = Fingerprint::new();
    hash.bytes(MESH_ASSIGNMENT_SEMANTICS_VERSION.as_bytes());
    hash.usize(assignments.len());
    let mut visited = 0usize;
    for assignment in assignments {
        poll(cx, "assignment-output-fingerprint", visited)?;
        visited += 1;
        hash.framed(assignment.subject.as_bytes());
        hash.u64(assignment.selector_fingerprint);
        hash.u8(u8::from(assignment.allow_overlap));
        hash.usize(assignment.faces.len());
        for face in &assignment.faces {
            poll(cx, "assignment-output-fingerprint", visited)?;
            visited += 1;
            hash.u32(*face);
        }
    }
    Ok(hash.finish())
}

struct Fingerprint(u64);

impl Fingerprint {
    const fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    fn bytes(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 ^= u64::from(*byte);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    fn framed(&mut self, bytes: &[u8]) {
        self.usize(bytes.len());
        self.bytes(bytes);
    }

    fn u8(&mut self, value: u8) {
        self.bytes(&[value]);
    }

    fn u32(&mut self, value: u32) {
        self.bytes(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }

    fn usize(&mut self, value: usize) {
        self.bytes(&u64::try_from(value).unwrap_or(u64::MAX).to_le_bytes());
    }

    const fn finish(self) -> u64 {
        self.0
    }
}

fn soup_triangle(soup: &Soup, face: usize) -> [Point3; 3] {
    let [a, b, c] = soup.triangles[face];
    [
        soup.positions[a as usize],
        soup.positions[b as usize],
        soup.positions[c as usize],
    ]
}

fn clone_faces(faces: &[u32]) -> Result<Vec<u32>, AssignmentRefusal> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(faces.len())
        .map_err(|_| allocation_refusal("selected faces", faces.len()))?;
    output.extend_from_slice(faces);
    Ok(output)
}

fn checkpoint(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), AssignmentRefusal> {
    cx.checkpoint().map_err(|_| AssignmentRefusal {
        code: "mesh-assignment-cancelled",
        what: format!("cancellation observed at {stage} position {at}"),
        fix: "retry the atomic assignment operation under a live deterministic context".to_string(),
    })
}

fn poll(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), AssignmentRefusal> {
    if at % MESH_ASSIGNMENT_POLL_STRIDE == 0 {
        checkpoint(cx, stage, at)?;
    }
    Ok(())
}

fn selector_refusal(what: impl Into<String>) -> AssignmentRefusal {
    refusal(
        "mesh-assignment-invalid-selector",
        what,
        "supply finite, ordered selector parameters in the declared length unit",
    )
}

fn resource_refusal(resource: &str, requested: usize, limit: usize) -> AssignmentRefusal {
    refusal(
        "mesh-assignment-resource-limit",
        format!("{resource} request {requested} exceeds limit {limit}"),
        "reduce the input/plan or explicitly raise the bounded assignment envelope",
    )
}

fn allocation_refusal(resource: &str, requested: usize) -> AssignmentRefusal {
    refusal(
        "mesh-assignment-allocation",
        format!("could not reserve {requested} {resource}"),
        "reduce the input or available-memory demand and retry",
    )
}

fn refusal(
    code: &'static str,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> AssignmentRefusal {
    AssignmentRefusal {
        code,
        what: what.into(),
        fix: fix.into(),
    }
}

fn validate_finite(value: f64, name: &str) -> Result<(), AssignmentRefusal> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(selector_refusal(format!("{name} must be finite")))
    }
}

fn validate_nonnegative(value: f64, name: &str) -> Result<(), AssignmentRefusal> {
    validate_finite(value, name)?;
    if value < 0.0 {
        Err(selector_refusal(format!("{name} must be non-negative")))
    } else {
        Ok(())
    }
}

fn validate_finite_vector(value: [f64; 3], name: &str) -> Result<(), AssignmentRefusal> {
    if finite3(value) {
        Ok(())
    } else {
        Err(selector_refusal(format!("{name} must be finite")))
    }
}

fn validate_direction(value: [f64; 3], name: &str) -> Result<(), AssignmentRefusal> {
    let magnitude = norm(value);
    if magnitude.is_finite() && magnitude > 0.0 {
        Ok(())
    } else {
        Err(selector_refusal(format!(
            "{name} must have a finite, nonzero magnitude"
        )))
    }
}

fn finite3(value: [f64; 3]) -> bool {
    value.into_iter().all(f64::is_finite)
}

fn add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn delta(left: Point3, right: Point3) -> [f64; 3] {
    [left.x - right.x, left.y - right.y, left.z - right.z]
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

fn push_json_string(output: &mut String, value: &str) {
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
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
