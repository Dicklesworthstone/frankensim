//! Strict native decoding for a triangle-only ISO 10303 faceted-B-rep subset.
//!
//! This module interprets only one caller-selected, root-reachable chain:
//! `FACETED_BREP -> CLOSED_SHELL -> (FACE | FACE_SURFACE) -> FACE_OUTER_BOUND
//! -> POLY_LOOP -> CARTESIAN_POINT`. A `FACE_SURFACE` must additionally resolve
//! through `PLANE -> AXIS2_PLACEMENT_3D` to a 3-D location and optional 3-D
//! directions. Every face must contain one triangular outer loop. Plane-backed
//! faces are checked for coplanarity and orientation before the resulting soup
//! is routed through [`crate::import_step_tessellation`], which owns repair,
//! bounded mesh-integrity checks, evidence composition, and SDF sampling.
//! A declared AP203/AP214-family schema is an admission label only; this module
//! does not validate an EXPRESS schema, representation context, units, product
//! linkage, application-protocol global rules, or general surface geometry.

use crate::step::{ParsedStep, StepEntity, StepInstance, StepValue};
use crate::step_import::{
    MAX_STEP_TESSELLATION_TRIANGLES, MAX_STEP_TESSELLATION_VERTICES,
    MAX_STEP_TOPOLOGY_AUXILIARY_BYTES, StepImportOutcome, StepImportRefusal,
    StepTessellatorIdentity, import_step_tessellation,
};
use fs_evidence::{NumericalCertificate, vv::UnitId};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_rep_mesh::Soup;
use std::{cmp::Reverse, collections::BinaryHeap};

/// Versioned meaning of the strict triangular faceted-B-rep decoder.
pub const STEP_FACETED_DECODER_VERSION: &str = "step-triangular-faceted-brep-v2";
/// Stable materializer name retained by the downstream STEP import receipt.
pub const STEP_FACETED_MATERIALIZER_NAME: &str = "fs-io-native-faceted-brep";

const CANCELLATION_POLL_STRIDE: usize = 4_096;
const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Exact schema declaration admitted by the resource-entity decoder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepFacetedProfile {
    /// `CONFIG_CONTROL_DESIGN`, commonly associated with AP203.
    ConfigControlDesign,
    /// `AUTOMOTIVE_DESIGN`, commonly associated with AP214.
    AutomotiveDesign,
}

impl StepFacetedProfile {
    /// Exact admitted `FILE_SCHEMA` string.
    #[must_use]
    pub const fn schema_identifier(self) -> &'static str {
        match self {
            Self::ConfigControlDesign => "CONFIG_CONTROL_DESIGN",
            Self::AutomotiveDesign => "AUTOMOTIVE_DESIGN",
        }
    }
}

/// Defensive bounds for native faceted-B-rep semantic materialization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StepFacetedLimits {
    /// Maximum distinct reachable polygon `CARTESIAN_POINT` instances.
    pub max_vertices: usize,
    /// Maximum reachable triangular faces.
    pub max_triangles: usize,
    /// Maximum admitted logical element-payload bytes.
    ///
    /// Container headers and allocator rounding are outside this portable
    /// estimate; allocation failure still maps to [`StepFacetedRefusal::Resource`].
    pub max_auxiliary_bytes: usize,
}

impl Default for StepFacetedLimits {
    fn default() -> Self {
        Self {
            max_vertices: MAX_STEP_TESSELLATION_VERTICES,
            max_triangles: MAX_STEP_TESSELLATION_TRIANGLES,
            max_auxiliary_bytes: MAX_STEP_TOPOLOGY_AUXILIARY_BYTES,
        }
    }
}

/// Fail-closed native faceted-B-rep decoding error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepFacetedRefusal {
    /// The declared schema set is missing, ambiguous, or outside the subset.
    Schema {
        /// Exact-source syntax fingerprint.
        source_fingerprint: u64,
        /// Actionable diagnosis.
        what: String,
    },
    /// One root-reachable instance violates the admitted entity shape.
    Entity {
        /// Exact-source syntax fingerprint.
        source_fingerprint: u64,
        /// Offending or expected instance ID.
        instance_id: u64,
        /// Stable relationship from the selected root.
        relationship: &'static str,
        /// Actionable diagnosis.
        what: String,
    },
    /// A semantic-work or allocation bound refused the request.
    Resource {
        /// Exact-source syntax fingerprint.
        source_fingerprint: u64,
        /// Decoder stage.
        stage: &'static str,
        /// Actionable diagnosis.
        what: String,
    },
    /// The caller context observed cancellation before publication.
    Cancelled {
        /// Exact-source syntax fingerprint.
        source_fingerprint: u64,
        /// Decoder stage that observed cancellation.
        stage: &'static str,
    },
}

impl core::fmt::Display for StepFacetedRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Schema { what, .. } => write!(f, "STEP faceted schema refused: {what}"),
            Self::Entity {
                instance_id,
                relationship,
                what,
                ..
            } => write!(
                f,
                "STEP faceted entity #{instance_id} at {relationship} refused: {what}"
            ),
            Self::Resource { stage, what, .. } => {
                write!(
                    f,
                    "STEP faceted resource admission refused at {stage}: {what}"
                )
            }
            Self::Cancelled { stage, .. } => {
                write!(f, "STEP faceted decoding cancelled at {stage}")
            }
        }
    }
}

impl std::error::Error for StepFacetedRefusal {}

/// Native decoding or downstream mesh/SDF refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum StepFacetedImportRefusal {
    /// The root-reachable faceted closure could not be materialized.
    Decode(StepFacetedRefusal),
    /// The materialized soup failed the existing topology/evidence/SDF handoff.
    Import {
        /// Successful native materialization retained for refusal provenance.
        decoder_receipt: StepFacetedReceipt,
        /// Downstream structured refusal.
        error: StepImportRefusal,
    },
}

impl core::fmt::Display for StepFacetedImportRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Decode(error) => error.fmt(f),
            Self::Import { error, .. } => error.fmt(f),
        }
    }
}

impl std::error::Error for StepFacetedImportRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Decode(error) => Some(error),
            Self::Import { error, .. } => Some(error),
        }
    }
}

/// Source-bound receipt for one strict triangular faceted-B-rep closure.
#[derive(Debug, Clone, PartialEq)]
pub struct StepFacetedReceipt {
    source_fingerprint: u64,
    canonical_layout_fingerprint: u64,
    profile: StepFacetedProfile,
    root_id: u64,
    shell_id: u64,
    vertex_count: usize,
    triangle_count: usize,
    bare_face_count: usize,
    plane_face_count: usize,
    reversed_bounds: usize,
    coordinate_conversion: NumericalCertificate,
    plane_consistency: NumericalCertificate,
    materialization_deviation: NumericalCertificate,
    semantic_fingerprint: u64,
    limits: StepFacetedLimits,
}

impl StepFacetedReceipt {
    /// Exact-source syntax fingerprint.
    #[must_use]
    pub const fn source_fingerprint(&self) -> u64 {
        self.source_fingerprint
    }

    /// Canonical-layout syntax fingerprint.
    #[must_use]
    pub const fn canonical_layout_fingerprint(&self) -> u64 {
        self.canonical_layout_fingerprint
    }

    /// Exact admitted schema declaration.
    #[must_use]
    pub const fn profile(&self) -> StepFacetedProfile {
        self.profile
    }

    /// Caller-selected `FACETED_BREP` instance ID.
    #[must_use]
    pub const fn root_id(&self) -> u64 {
        self.root_id
    }

    /// Reachable `CLOSED_SHELL` instance ID.
    #[must_use]
    pub const fn shell_id(&self) -> u64 {
        self.shell_id
    }

    /// Distinct reachable Cartesian-point count.
    #[must_use]
    pub const fn vertex_count(&self) -> usize {
        self.vertex_count
    }

    /// Reachable triangular face count.
    #[must_use]
    pub const fn triangle_count(&self) -> usize {
        self.triangle_count
    }

    /// Reachable resource-level `FACE` count.
    #[must_use]
    pub const fn bare_face_count(&self) -> usize {
        self.bare_face_count
    }

    /// Reachable plane-backed `FACE_SURFACE` count.
    #[must_use]
    pub const fn plane_face_count(&self) -> usize {
        self.plane_face_count
    }

    /// Face-bound orientations reversed while materializing loops.
    #[must_use]
    pub const fn reversed_bounds(&self) -> usize {
        self.reversed_bounds
    }

    /// Conservative estimated spatial displacement from decimal-to-f64 conversion.
    #[must_use]
    pub const fn coordinate_conversion(&self) -> NumericalCertificate {
        self.coordinate_conversion
    }

    /// Estimated maximum accepted distance between a face vertex and its plane.
    #[must_use]
    pub const fn plane_consistency(&self) -> NumericalCertificate {
        self.plane_consistency
    }

    /// Estimated spatial deviation handed to the downstream STEP importer.
    #[must_use]
    pub const fn materialization_deviation(&self) -> NumericalCertificate {
        self.materialization_deviation
    }

    /// Deterministic non-cryptographic identity of the admitted closure and soup.
    #[must_use]
    pub const fn semantic_fingerprint(&self) -> u64 {
        self.semantic_fingerprint
    }

    /// Exact semantic limits used for admission.
    #[must_use]
    pub const fn limits(&self) -> StepFacetedLimits {
        self.limits
    }

    /// Canonical JSON suitable for a schema-declaration-gated import event.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"kind\":\"step-triangular-faceted-brep-receipt\",\
             \"authority\":\"schema-declaration-gated-resource-decoding\",\
             \"decoder_version\":\"{}\",\"fs_io_version\":\"{}\",\
             \"source_fingerprint_fnv1a64\":\"{:016x}\",\
             \"canonical_layout_fingerprint_fnv1a64\":\"{:016x}\",\
             \"semantic_fingerprint_fnv1a64\":\"{:016x}\",\
             \"schema\":\"{}\",\"root_id\":{},\"shell_id\":{},\
             \"vertices\":{},\"triangles\":{},\"bare_faces\":{},\"plane_faces\":{},\
             \"reversed_bounds\":{},\
             \"coordinate_conversion\":{{\"kind\":\"estimate\",\"lo\":{},\"hi\":{}}},\
             \"plane_consistency\":{{\"kind\":\"estimate\",\"lo\":{},\"hi\":{}}},\
             \"materialization_deviation\":{{\"kind\":\"estimate\",\"lo\":{},\"hi\":{}}},\
             \"limits\":{{\"vertices\":{},\"triangles\":{},\"auxiliary_bytes\":{}}},\
             \"no_claim\":\"pinned FACE and plane-backed FACE_SURFACE resource subset only; no full EXPRESS or AP conformance, representation/unit context, product linkage, general surfaces, self-intersection certificate, or topology authority\"}}",
            STEP_FACETED_DECODER_VERSION,
            crate::VERSION,
            self.source_fingerprint,
            self.canonical_layout_fingerprint,
            self.semantic_fingerprint,
            self.profile.schema_identifier(),
            self.root_id,
            self.shell_id,
            self.vertex_count,
            self.triangle_count,
            self.bare_face_count,
            self.plane_face_count,
            self.reversed_bounds,
            self.coordinate_conversion.lo,
            self.coordinate_conversion.hi,
            self.plane_consistency.lo,
            self.plane_consistency.hi,
            self.materialization_deviation.lo,
            self.materialization_deviation.hi,
            self.limits.max_vertices,
            self.limits.max_triangles,
            self.limits.max_auxiliary_bytes,
        )
    }
}

/// Sealed triangle soup and receipt produced by native faceted decoding.
#[derive(Debug)]
pub struct DecodedFacetedBrep {
    soup: Soup,
    receipt: StepFacetedReceipt,
}

impl DecodedFacetedBrep {
    /// Read-only materialized triangle soup.
    #[must_use]
    pub const fn soup(&self) -> &Soup {
        &self.soup
    }

    /// Read-only source-bound decoder receipt.
    #[must_use]
    pub const fn receipt(&self) -> &StepFacetedReceipt {
        &self.receipt
    }

    /// Consume the sealed result for the downstream topology handoff.
    #[must_use]
    pub fn into_parts(self) -> (Soup, StepFacetedReceipt) {
        (self.soup, self.receipt)
    }
}

/// Estimated SDF import paired with the native semantic-decoder receipt.
#[derive(Debug)]
pub struct StepFacetedImportOutcome {
    decoder_receipt: StepFacetedReceipt,
    import: StepImportOutcome,
}

impl StepFacetedImportOutcome {
    /// Root-reachable semantic materialization receipt.
    #[must_use]
    pub const fn decoder_receipt(&self) -> &StepFacetedReceipt {
        &self.decoder_receipt
    }

    /// Existing topology/evidence/SDF import result.
    #[must_use]
    pub const fn import(&self) -> &StepImportOutcome {
        &self.import
    }
}

#[derive(Debug, Clone, Copy)]
struct RawFace {
    face_id: u64,
    bound_id: u64,
    loop_id: u64,
    point_ids: [u64; 3],
    plane: Option<RawPlane>,
}

#[derive(Debug, Clone, Copy)]
struct RawPlane {
    plane_id: u64,
    placement_id: u64,
    location_id: u64,
    axis_id: Option<u64>,
    ref_direction_id: Option<u64>,
    location: Point3,
    axis: [f64; 3],
    ref_direction: Option<[f64; 3]>,
    same_sense: bool,
}

/// Decode one caller-selected triangular `FACETED_BREP` closure.
///
/// Reachable instances must use the strict entity closure documented at module
/// scope. Unknown instances outside the selected closure are ignored.
/// The result is a deterministic soup, not topology or application-protocol
/// authority.
///
/// # Errors
/// [`StepFacetedRefusal`] for schema drift, malformed reachable entities,
/// bounds/allocation refusal, non-finite coordinate conversion, or cancellation.
#[allow(clippy::too_many_lines)]
pub fn decode_faceted_brep_with_limits(
    parsed: &ParsedStep,
    root_id: u64,
    limits: StepFacetedLimits,
    cx: &Cx<'_>,
) -> Result<DecodedFacetedBrep, StepFacetedRefusal> {
    let source_fingerprint = parsed.receipt().source_fingerprint();
    checkpoint(cx, source_fingerprint, "entry")?;
    validate_limits(limits, source_fingerprint)?;
    let profile = admitted_profile(parsed, source_fingerprint)?;
    if root_id == 0 {
        return Err(entity_refusal(
            source_fingerprint,
            root_id,
            "root",
            "FACETED_BREP instance IDs must be positive",
        ));
    }

    let instances = &parsed.document().instances;
    let index_bytes = instances
        .len()
        .checked_mul(core::mem::size_of::<(u64, usize)>())
        .ok_or_else(|| {
            resource_refusal(
                source_fingerprint,
                "instance-index",
                "byte estimate overflowed",
            )
        })?;
    let index_sort_heap_bytes = instances
        .len()
        .div_ceil(CANCELLATION_POLL_STRIDE)
        .checked_mul(core::mem::size_of::<Reverse<(u64, usize, usize)>>())
        .ok_or_else(|| {
            resource_refusal(
                source_fingerprint,
                "instance-index",
                "sort-heap byte estimate overflowed",
            )
        })?;
    let index_sort_bytes = index_bytes
        .checked_mul(2)
        .and_then(|bytes| bytes.checked_add(index_sort_heap_bytes))
        .ok_or_else(|| {
            resource_refusal(
                source_fingerprint,
                "instance-index",
                "sort byte estimate overflowed",
            )
        })?;
    ensure_auxiliary_bound(
        index_sort_bytes,
        limits,
        source_fingerprint,
        "instance-index",
    )?;
    let mut index = Vec::new();
    reserve_exact(
        &mut index,
        instances.len(),
        source_fingerprint,
        "instance-index",
    )?;
    for (offset, instance) in instances.iter().enumerate() {
        poll(cx, source_fingerprint, "instance-index", offset)?;
        index.push((instance.id, offset));
    }
    cancellable_sort_by_key(
        &mut index,
        |entry| entry.0,
        cx,
        source_fingerprint,
        "instance-index-sort",
    )?;
    for (offset, pair) in index.windows(2).enumerate() {
        poll(
            cx,
            source_fingerprint,
            "instance-index-duplicate-scan",
            offset,
        )?;
        if pair[0].0 == pair[1].0 {
            return Err(entity_refusal(
                source_fingerprint,
                pair[0].0,
                "instance-index",
                "duplicate instance ID survived syntax admission",
            ));
        }
    }

    let root = simple_entity(
        instances,
        &index,
        root_id,
        "FACETED_BREP",
        "root",
        source_fingerprint,
    )?;
    require_shape(root, root_id, 2, "root", source_fingerprint)?;
    let shell_id = require_reference(
        &root.parameters[1],
        root_id,
        "root.outer",
        source_fingerprint,
    )?;
    let shell = simple_entity(
        instances,
        &index,
        shell_id,
        "CLOSED_SHELL",
        "root.outer",
        source_fingerprint,
    )?;
    require_shape(shell, shell_id, 2, "root.outer", source_fingerprint)?;
    let shell_faces = require_aggregate(
        &shell.parameters[1],
        shell_id,
        "root.outer.faces",
        source_fingerprint,
    )?;
    if shell_faces.is_empty() {
        return Err(entity_refusal(
            source_fingerprint,
            shell_id,
            "root.outer.faces",
            "CLOSED_SHELL must contain at least one face",
        ));
    }
    if shell_faces.len() > limits.max_triangles {
        return Err(resource_refusal(
            source_fingerprint,
            "shell-faces",
            format!(
                "CLOSED_SHELL has {} faces; cap is {}",
                shell_faces.len(),
                limits.max_triangles
            ),
        ));
    }

    let face_count = shell_faces.len();
    let per_face_bytes = core::mem::size_of::<RawFace>()
        + 3 * core::mem::size_of::<u64>()
        + 3 * core::mem::size_of::<u64>()
        + core::mem::size_of::<[u32; 3]>()
        + 3 * core::mem::size_of::<Point3>();
    let semantic_bytes = face_count
        .checked_mul(per_face_bytes)
        .and_then(|bytes| bytes.checked_add(index_bytes))
        .ok_or_else(|| {
            resource_refusal(
                source_fingerprint,
                "semantic-plan",
                "byte estimate overflowed",
            )
        })?;
    ensure_auxiliary_bound(semantic_bytes, limits, source_fingerprint, "semantic-plan")?;

    let mut face_ids = Vec::new();
    reserve_exact(&mut face_ids, face_count, source_fingerprint, "shell-faces")?;
    for (face_offset, value) in shell_faces.iter().enumerate() {
        poll(
            cx,
            source_fingerprint,
            "shell-face-reference-plan",
            face_offset,
        )?;
        face_ids.push(require_reference(
            value,
            shell_id,
            "root.outer.faces",
            source_fingerprint,
        )?);
    }
    cancellable_sort_by_key(
        &mut face_ids,
        |face_id| *face_id,
        cx,
        source_fingerprint,
        "shell-face-sort",
    )?;
    if let Some(duplicate) = adjacent_duplicate(
        &face_ids,
        cx,
        source_fingerprint,
        "shell-face-duplicate-scan",
    )? {
        return Err(entity_refusal(
            source_fingerprint,
            duplicate,
            "root.outer.faces",
            "CLOSED_SHELL face SET contains a duplicate reference",
        ));
    }

    let mut raw_faces = Vec::new();
    let mut used_bounds = Vec::new();
    let mut used_loops = Vec::new();
    reserve_exact(&mut raw_faces, face_count, source_fingerprint, "face-plan")?;
    reserve_exact(
        &mut used_bounds,
        face_count,
        source_fingerprint,
        "bound-uniqueness",
    )?;
    reserve_exact(
        &mut used_loops,
        face_count,
        source_fingerprint,
        "loop-uniqueness",
    )?;
    let mut reversed_bounds = 0usize;
    let mut bare_face_count = 0usize;
    let mut plane_face_count = 0usize;
    let mut coordinate_upper = 0.0f64;
    for (face_offset, &face_id) in face_ids.iter().enumerate() {
        poll(cx, source_fingerprint, "face-traversal", face_offset)?;
        let face = simple_entity_any(
            instances,
            &index,
            face_id,
            "root.outer.faces[]",
            source_fingerprint,
        )?;
        let plane = match face.name.as_str() {
            "FACE" => {
                require_shape(face, face_id, 2, "root.outer.faces[]", source_fingerprint)?;
                bare_face_count += 1;
                None
            }
            "FACE_SURFACE" => {
                require_shape(face, face_id, 4, "root.outer.faces[]", source_fingerprint)?;
                let plane_id = require_reference(
                    &face.parameters[2],
                    face_id,
                    "face_surface.face_geometry",
                    source_fingerprint,
                )?;
                let same_sense = require_boolean(
                    &face.parameters[3],
                    face_id,
                    "face_surface.same_sense",
                    source_fingerprint,
                )?;
                let (plane, location_upper) =
                    parse_raw_plane(instances, &index, plane_id, same_sense, source_fingerprint)?;
                coordinate_upper = coordinate_upper.max(location_upper);
                plane_face_count += 1;
                Some(plane)
            }
            other => {
                return Err(entity_refusal(
                    source_fingerprint,
                    face_id,
                    "root.outer.faces[]",
                    format!("expected FACE or FACE_SURFACE, found {other}"),
                ));
            }
        };
        let bounds = require_aggregate(
            &face.parameters[1],
            face_id,
            "face.bounds",
            source_fingerprint,
        )?;
        if bounds.len() != 1 {
            return Err(entity_refusal(
                source_fingerprint,
                face_id,
                "face.bounds",
                format!(
                    "triangle-only FACE requires exactly one outer bound; got {}",
                    bounds.len()
                ),
            ));
        }
        let bound_id =
            require_reference(&bounds[0], face_id, "face.bounds[0]", source_fingerprint)?;
        let bound = simple_entity(
            instances,
            &index,
            bound_id,
            "FACE_OUTER_BOUND",
            "face.bounds[0]",
            source_fingerprint,
        )?;
        require_shape(bound, bound_id, 3, "face.bounds[0]", source_fingerprint)?;
        let loop_id = require_reference(
            &bound.parameters[1],
            bound_id,
            "face.outer_bound.loop",
            source_fingerprint,
        )?;
        let orientation = require_boolean(
            &bound.parameters[2],
            bound_id,
            "face.outer_bound.orientation",
            source_fingerprint,
        )?;
        let reverse = !orientation;
        let poly_loop = simple_entity(
            instances,
            &index,
            loop_id,
            "POLY_LOOP",
            "face.outer_bound.loop",
            source_fingerprint,
        )?;
        require_shape(
            poly_loop,
            loop_id,
            2,
            "face.outer_bound.loop",
            source_fingerprint,
        )?;
        let polygon = require_aggregate(
            &poly_loop.parameters[1],
            loop_id,
            "poly_loop.polygon",
            source_fingerprint,
        )?;
        if polygon.len() != 3 {
            return Err(entity_refusal(
                source_fingerprint,
                loop_id,
                "poly_loop.polygon",
                format!(
                    "triangle-only POLY_LOOP requires exactly three points; got {}",
                    polygon.len()
                ),
            ));
        }
        let mut point_ids = [0u64; 3];
        for (corner, value) in polygon.iter().enumerate() {
            point_ids[corner] =
                require_reference(value, loop_id, "poly_loop.polygon[]", source_fingerprint)?;
        }
        if point_ids[0] == point_ids[1]
            || point_ids[0] == point_ids[2]
            || point_ids[1] == point_ids[2]
        {
            return Err(entity_refusal(
                source_fingerprint,
                loop_id,
                "poly_loop.polygon",
                "POLY_LOOP point references must be unique",
            ));
        }
        if reverse {
            point_ids.swap(1, 2);
            reversed_bounds += 1;
        }
        used_bounds.push(bound_id);
        used_loops.push(loop_id);
        raw_faces.push(RawFace {
            face_id,
            bound_id,
            loop_id,
            point_ids,
            plane,
        });
    }

    cancellable_sort_by_key(
        &mut used_bounds,
        |bound_id| *bound_id,
        cx,
        source_fingerprint,
        "bound-uniqueness-sort",
    )?;
    if let Some(duplicate) =
        adjacent_duplicate(&used_bounds, cx, source_fingerprint, "bound-duplicate-scan")?
    {
        return Err(entity_refusal(
            source_fingerprint,
            duplicate,
            "face.bounds",
            "FACE_OUTER_BOUND is reused by more than one face",
        ));
    }
    cancellable_sort_by_key(
        &mut used_loops,
        |loop_id| *loop_id,
        cx,
        source_fingerprint,
        "loop-uniqueness-sort",
    )?;
    if let Some(duplicate) =
        adjacent_duplicate(&used_loops, cx, source_fingerprint, "loop-duplicate-scan")?
    {
        return Err(entity_refusal(
            source_fingerprint,
            duplicate,
            "face.outer_bound.loop",
            "POLY_LOOP is reused by more than one face",
        ));
    }

    let point_reference_count = face_count.checked_mul(3).ok_or_else(|| {
        resource_refusal(
            source_fingerprint,
            "point-plan",
            "point-reference count overflowed",
        )
    })?;
    let mut point_ids = Vec::new();
    reserve_exact(
        &mut point_ids,
        point_reference_count,
        source_fingerprint,
        "point-plan",
    )?;
    for (face_offset, face) in raw_faces.iter().enumerate() {
        poll(cx, source_fingerprint, "point-reference-plan", face_offset)?;
        point_ids.extend_from_slice(&face.point_ids);
    }
    cancellable_sort_by_key(
        &mut point_ids,
        |point_id| *point_id,
        cx,
        source_fingerprint,
        "point-reference-sort",
    )?;
    dedup_sorted(
        &mut point_ids,
        cx,
        source_fingerprint,
        "point-reference-dedup",
    )?;
    checkpoint(cx, source_fingerprint, "point-reference-sort")?;
    if point_ids.len() > limits.max_vertices {
        return Err(resource_refusal(
            source_fingerprint,
            "point-plan",
            format!(
                "reachable closure has {} points; cap is {}",
                point_ids.len(),
                limits.max_vertices
            ),
        ));
    }

    let mut positions = Vec::new();
    reserve_exact(
        &mut positions,
        point_ids.len(),
        source_fingerprint,
        "point-materialization",
    )?;
    for (point_offset, &point_id) in point_ids.iter().enumerate() {
        poll(
            cx,
            source_fingerprint,
            "point-materialization",
            point_offset,
        )?;
        let (position, spatial_estimate) = parse_cartesian_point(
            instances,
            &index,
            point_id,
            "poly_loop.polygon[]",
            source_fingerprint,
        )?;
        coordinate_upper = coordinate_upper.max(spatial_estimate);
        positions.push(position);
    }

    let mut triangles = Vec::new();
    reserve_exact(
        &mut triangles,
        raw_faces.len(),
        source_fingerprint,
        "triangle-materialization",
    )?;
    let mut plane_consistency_upper = 0.0f64;
    for (face_offset, face) in raw_faces.iter().enumerate() {
        poll(
            cx,
            source_fingerprint,
            "triangle-materialization",
            face_offset,
        )?;
        let mut triangle = [0u32; 3];
        for (corner, point_id) in face.point_ids.iter().enumerate() {
            let ordinal = point_ids.binary_search(point_id).map_err(|_| {
                entity_refusal(
                    source_fingerprint,
                    face.face_id,
                    "poly_loop.polygon[]",
                    format!("point reference #{point_id} disappeared from the semantic plan"),
                )
            })?;
            triangle[corner] = u32::try_from(ordinal).map_err(|_| {
                resource_refusal(
                    source_fingerprint,
                    "triangle-materialization",
                    format!("compact point ordinal {ordinal} exceeds u32"),
                )
            })?;
        }
        if let Some(plane) = face.plane {
            plane_consistency_upper = plane_consistency_upper.max(validate_plane_face(
                face,
                plane,
                triangle,
                &positions,
                coordinate_upper,
                source_fingerprint,
            )?);
        }
        triangles.push(triangle);
    }

    let materialization_upper = if plane_face_count == 0 {
        coordinate_upper
    } else {
        outward_nonnegative_sum(coordinate_upper, plane_consistency_upper).ok_or_else(|| {
            entity_refusal(
                source_fingerprint,
                root_id,
                "materialization-deviation",
                "coordinate conversion and plane-consistency estimates do not have a finite sum",
            )
        })?
    };

    let semantic_fingerprint = semantic_fingerprint(
        profile,
        root_id,
        shell_id,
        &point_ids,
        &positions,
        &raw_faces,
        &triangles,
        source_fingerprint,
        cx,
    )?;
    checkpoint(cx, source_fingerprint, "publication")?;
    let receipt = StepFacetedReceipt {
        source_fingerprint,
        canonical_layout_fingerprint: parsed.receipt().canonical_layout_fingerprint(),
        profile,
        root_id,
        shell_id,
        vertex_count: positions.len(),
        triangle_count: triangles.len(),
        bare_face_count,
        plane_face_count,
        reversed_bounds,
        coordinate_conversion: NumericalCertificate::estimate(0.0, coordinate_upper),
        plane_consistency: NumericalCertificate::estimate(0.0, plane_consistency_upper),
        materialization_deviation: NumericalCertificate::estimate(0.0, materialization_upper),
        semantic_fingerprint,
        limits,
    };
    Ok(DecodedFacetedBrep {
        soup: Soup {
            positions,
            triangles,
        },
        receipt,
    })
}

/// Decode a strict triangular faceted B-rep and run the existing topology and
/// estimated-SDF handoff.
///
/// The caller supplies the unit because this closure deliberately ignores
/// representation context. Decimal-coordinate conversion and accepted
/// face-plane residual are retained as an estimated materialization deviation.
/// The returned decoder receipt remains separate from the existing topology/SDF
/// receipt.
///
/// # Errors
/// [`StepFacetedImportRefusal`] for native decoding or downstream import refusal.
pub fn import_faceted_brep(
    parsed: &ParsedStep,
    root_id: u64,
    length_unit: UnitId,
    target_h: f64,
    cx: &Cx<'_>,
) -> Result<StepFacetedImportOutcome, StepFacetedImportRefusal> {
    let decoded =
        decode_faceted_brep_with_limits(parsed, root_id, StepFacetedLimits::default(), cx)
            .map_err(StepFacetedImportRefusal::Decode)?;
    let (soup, decoder_receipt) = decoded.into_parts();
    let configuration_fingerprint = fnv1a64_update(
        fnv1a64_update(FNV_OFFSET, STEP_FACETED_DECODER_VERSION.as_bytes()),
        &decoder_receipt.semantic_fingerprint().to_le_bytes(),
    );
    let import = import_step_tessellation(
        parsed,
        soup,
        StepTessellatorIdentity {
            name: STEP_FACETED_MATERIALIZER_NAME.to_string(),
            version: STEP_FACETED_DECODER_VERSION.to_string(),
            configuration_fingerprint,
        },
        decoder_receipt.materialization_deviation(),
        length_unit,
        target_h,
        cx,
    )
    .map_err(|error| StepFacetedImportRefusal::Import {
        decoder_receipt: decoder_receipt.clone(),
        error,
    })?;
    Ok(StepFacetedImportOutcome {
        decoder_receipt,
        import,
    })
}

fn admitted_profile(
    parsed: &ParsedStep,
    source_fingerprint: u64,
) -> Result<StepFacetedProfile, StepFacetedRefusal> {
    let schemas = parsed.receipt().schema_identifiers();
    if schemas.len() != 1 {
        return Err(StepFacetedRefusal::Schema {
            source_fingerprint,
            what: format!(
                "exactly one unambiguous schema declaration is required; got {:?}",
                schemas
            ),
        });
    }
    match schemas[0].as_str() {
        "CONFIG_CONTROL_DESIGN" => Ok(StepFacetedProfile::ConfigControlDesign),
        "AUTOMOTIVE_DESIGN" => Ok(StepFacetedProfile::AutomotiveDesign),
        schema => Err(StepFacetedRefusal::Schema {
            source_fingerprint,
            what: format!(
                "schema {schema:?} is outside the admitted CONFIG_CONTROL_DESIGN/AUTOMOTIVE_DESIGN subset"
            ),
        }),
    }
}

fn validate_limits(
    limits: StepFacetedLimits,
    source_fingerprint: u64,
) -> Result<(), StepFacetedRefusal> {
    for (name, offered, hard) in [
        (
            "max_vertices",
            limits.max_vertices,
            MAX_STEP_TESSELLATION_VERTICES,
        ),
        (
            "max_triangles",
            limits.max_triangles,
            MAX_STEP_TESSELLATION_TRIANGLES,
        ),
        (
            "max_auxiliary_bytes",
            limits.max_auxiliary_bytes,
            MAX_STEP_TOPOLOGY_AUXILIARY_BYTES,
        ),
    ] {
        if offered == 0 || offered > hard {
            return Err(resource_refusal(
                source_fingerprint,
                "limits",
                format!("{name} {offered} is outside 1..={hard}"),
            ));
        }
    }
    Ok(())
}

fn simple_entity<'a>(
    instances: &'a [StepInstance],
    index: &[(u64, usize)],
    id: u64,
    expected: &'static str,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<&'a StepEntity, StepFacetedRefusal> {
    let entity = simple_entity_any(instances, index, id, relationship, source_fingerprint)?;
    if entity.name != expected {
        return Err(entity_refusal(
            source_fingerprint,
            id,
            relationship,
            format!("expected {expected}, found {}", entity.name),
        ));
    }
    Ok(entity)
}

fn simple_entity_any<'a>(
    instances: &'a [StepInstance],
    index: &[(u64, usize)],
    id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<&'a StepEntity, StepFacetedRefusal> {
    let offset = index
        .binary_search_by_key(&id, |entry| entry.0)
        .map(|position| index[position].1)
        .map_err(|_| {
            entity_refusal(
                source_fingerprint,
                id,
                relationship,
                "referenced instance does not exist",
            )
        })?;
    let instance = &instances[offset];
    if instance.components.len() != 1 {
        return Err(entity_refusal(
            source_fingerprint,
            id,
            relationship,
            format!(
                "reachable instance must be simple; found {} complex components",
                instance.components.len()
            ),
        ));
    }
    let entity = &instance.components[0];
    Ok(entity)
}

fn require_shape(
    entity: &StepEntity,
    id: u64,
    arity: usize,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<(), StepFacetedRefusal> {
    if entity.parameters.len() != arity {
        return Err(entity_refusal(
            source_fingerprint,
            id,
            relationship,
            format!(
                "{} requires {arity} physical parameters including inherited name; got {}",
                entity.name,
                entity.parameters.len()
            ),
        ));
    }
    if !matches!(&entity.parameters[0], StepValue::String(_)) {
        return Err(entity_refusal(
            source_fingerprint,
            id,
            relationship,
            format!("{} inherited name must be a STEP string", entity.name),
        ));
    }
    Ok(())
}

fn require_reference(
    value: &StepValue,
    owner_id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<u64, StepFacetedRefusal> {
    if let StepValue::Reference(id) = value {
        Ok(*id)
    } else {
        Err(entity_refusal(
            source_fingerprint,
            owner_id,
            relationship,
            "expected an entity reference",
        ))
    }
}

fn require_optional_reference(
    value: &StepValue,
    owner_id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<Option<u64>, StepFacetedRefusal> {
    match value {
        StepValue::Omitted => Ok(None),
        StepValue::Reference(id) => Ok(Some(*id)),
        _ => Err(entity_refusal(
            source_fingerprint,
            owner_id,
            relationship,
            "expected an entity reference or omitted value ($)",
        )),
    }
}

fn require_aggregate<'a>(
    value: &'a StepValue,
    owner_id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<&'a [StepValue], StepFacetedRefusal> {
    if let StepValue::Aggregate(values) = value {
        Ok(values)
    } else {
        Err(entity_refusal(
            source_fingerprint,
            owner_id,
            relationship,
            "expected an aggregate",
        ))
    }
}

fn require_boolean(
    value: &StepValue,
    owner_id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<bool, StepFacetedRefusal> {
    match value {
        StepValue::Enumeration(value) if value == "T" => Ok(true),
        StepValue::Enumeration(value) if value == "F" => Ok(false),
        _ => Err(entity_refusal(
            source_fingerprint,
            owner_id,
            relationship,
            "boolean value must be .T. or .F.",
        )),
    }
}

fn parse_raw_plane(
    instances: &[StepInstance],
    index: &[(u64, usize)],
    plane_id: u64,
    same_sense: bool,
    source_fingerprint: u64,
) -> Result<(RawPlane, f64), StepFacetedRefusal> {
    let plane = simple_entity(
        instances,
        index,
        plane_id,
        "PLANE",
        "face_surface.face_geometry",
        source_fingerprint,
    )?;
    require_shape(
        plane,
        plane_id,
        2,
        "face_surface.face_geometry",
        source_fingerprint,
    )?;
    let placement_id = require_reference(
        &plane.parameters[1],
        plane_id,
        "plane.position",
        source_fingerprint,
    )?;
    let placement = simple_entity(
        instances,
        index,
        placement_id,
        "AXIS2_PLACEMENT_3D",
        "plane.position",
        source_fingerprint,
    )?;
    require_shape(
        placement,
        placement_id,
        4,
        "plane.position",
        source_fingerprint,
    )?;
    let location_id = require_reference(
        &placement.parameters[1],
        placement_id,
        "plane.position.location",
        source_fingerprint,
    )?;
    let axis_id = require_optional_reference(
        &placement.parameters[2],
        placement_id,
        "plane.position.axis",
        source_fingerprint,
    )?;
    let ref_direction_id = require_optional_reference(
        &placement.parameters[3],
        placement_id,
        "plane.position.ref_direction",
        source_fingerprint,
    )?;
    let (location, location_upper) = parse_cartesian_point(
        instances,
        index,
        location_id,
        "plane.position.location",
        source_fingerprint,
    )?;
    let axis = match axis_id {
        Some(id) => parse_direction(
            instances,
            index,
            id,
            "plane.position.axis",
            source_fingerprint,
        )?,
        None => [0.0, 0.0, 1.0],
    };
    let ref_direction = match ref_direction_id {
        Some(id) => Some(parse_direction(
            instances,
            index,
            id,
            "plane.position.ref_direction",
            source_fingerprint,
        )?),
        None => None,
    };
    if axis_id.is_some() {
        if let Some(reference) = ref_direction {
            let separation = cross_norm(axis, reference);
            if !separation.is_finite() || separation <= 256.0 * f64::EPSILON {
                return Err(entity_refusal(
                    source_fingerprint,
                    placement_id,
                    "plane.position",
                    "explicit axis and ref_direction must be numerically non-parallel",
                ));
            }
        }
    }
    Ok((
        RawPlane {
            plane_id,
            placement_id,
            location_id,
            axis_id,
            ref_direction_id,
            location,
            axis,
            ref_direction,
            same_sense,
        },
        location_upper,
    ))
}

fn parse_cartesian_point(
    instances: &[StepInstance],
    index: &[(u64, usize)],
    point_id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<(Point3, f64), StepFacetedRefusal> {
    let point = simple_entity(
        instances,
        index,
        point_id,
        "CARTESIAN_POINT",
        relationship,
        source_fingerprint,
    )?;
    require_shape(point, point_id, 2, relationship, source_fingerprint)?;
    let coordinates = require_aggregate(
        &point.parameters[1],
        point_id,
        "cartesian_point.coordinates",
        source_fingerprint,
    )?;
    let values = parse_finite_triplet(
        coordinates,
        point_id,
        "cartesian_point.coordinates",
        "coordinate",
        source_fingerprint,
    )?;
    let point_ulp = values
        .iter()
        .fold(0.0f64, |upper, value| upper.max(ulp_magnitude(*value)));
    let spatial_estimate = point_ulp * 2.0;
    if !spatial_estimate.is_finite() {
        return Err(entity_refusal(
            source_fingerprint,
            point_id,
            "cartesian_point.coordinates",
            "finite coordinate has no finite conservative spatial conversion estimate",
        ));
    }
    Ok((
        Point3::new(values[0], values[1], values[2]),
        spatial_estimate,
    ))
}

fn parse_direction(
    instances: &[StepInstance],
    index: &[(u64, usize)],
    direction_id: u64,
    relationship: &'static str,
    source_fingerprint: u64,
) -> Result<[f64; 3], StepFacetedRefusal> {
    let direction = simple_entity(
        instances,
        index,
        direction_id,
        "DIRECTION",
        relationship,
        source_fingerprint,
    )?;
    require_shape(direction, direction_id, 2, relationship, source_fingerprint)?;
    let ratios = require_aggregate(
        &direction.parameters[1],
        direction_id,
        "direction.direction_ratios",
        source_fingerprint,
    )?;
    let values = parse_finite_triplet(
        ratios,
        direction_id,
        "direction.direction_ratios",
        "direction ratio",
        source_fingerprint,
    )?;
    normalize_direction(values).ok_or_else(|| {
        entity_refusal(
            source_fingerprint,
            direction_id,
            "direction.direction_ratios",
            "3-D direction ratios must define a finite nonzero direction",
        )
    })
}

fn parse_finite_triplet(
    values: &[StepValue],
    owner_id: u64,
    relationship: &'static str,
    noun: &'static str,
    source_fingerprint: u64,
) -> Result<[f64; 3], StepFacetedRefusal> {
    if values.len() != 3 {
        return Err(entity_refusal(
            source_fingerprint,
            owner_id,
            relationship,
            format!(
                "3-D faceted subset requires exactly three {noun}s; got {}",
                values.len()
            ),
        ));
    }
    let mut parsed = [0.0f64; 3];
    for (axis, value) in values.iter().enumerate() {
        let StepValue::Number(token) = value else {
            return Err(entity_refusal(
                source_fingerprint,
                owner_id,
                relationship,
                format!("{noun}s must be decimal number tokens"),
            ));
        };
        let number = token.parse::<f64>().map_err(|error| {
            entity_refusal(
                source_fingerprint,
                owner_id,
                relationship,
                format!("{noun} {token:?} does not convert to f64: {error}"),
            )
        })?;
        if !number.is_finite() {
            return Err(entity_refusal(
                source_fingerprint,
                owner_id,
                relationship,
                format!("{noun} {token:?} converts to non-finite f64"),
            ));
        }
        parsed[axis] = number;
    }
    Ok(parsed)
}

fn normalize_direction(values: [f64; 3]) -> Option<[f64; 3]> {
    let scale = values
        .iter()
        .fold(0.0f64, |largest, value| largest.max(value.abs()));
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let scaled = [values[0] / scale, values[1] / scale, values[2] / scale];
    let magnitude = (scaled[0] * scaled[0] + scaled[1] * scaled[1] + scaled[2] * scaled[2]).sqrt();
    if !magnitude.is_finite() || magnitude == 0.0 {
        return None;
    }
    Some([
        scaled[0] / magnitude,
        scaled[1] / magnitude,
        scaled[2] / magnitude,
    ])
}

fn validate_plane_face(
    face: &RawFace,
    plane: RawPlane,
    triangle: [u32; 3],
    positions: &[Point3],
    coordinate_upper: f64,
    source_fingerprint: u64,
) -> Result<f64, StepFacetedRefusal> {
    let mut points = [Point3::new(0.0, 0.0, 0.0); 3];
    for (corner, ordinal) in triangle.into_iter().enumerate() {
        points[corner] = *positions.get(ordinal as usize).ok_or_else(|| {
            entity_refusal(
                source_fingerprint,
                face.face_id,
                "face_surface.bounds",
                format!("materialized point ordinal {ordinal} is unavailable"),
            )
        })?;
    }

    let edge_a = [
        points[1].x - points[0].x,
        points[1].y - points[0].y,
        points[1].z - points[0].z,
    ];
    let edge_b = [
        points[2].x - points[0].x,
        points[2].y - points[0].y,
        points[2].z - points[0].z,
    ];
    let edge_scale = edge_a
        .iter()
        .chain(&edge_b)
        .fold(0.0f64, |largest, value| largest.max(value.abs()));
    if !edge_scale.is_finite() || edge_scale == 0.0 {
        return Err(entity_refusal(
            source_fingerprint,
            face.face_id,
            "face_surface.bounds",
            "plane-backed triangle has no finite nonzero edge scale",
        ));
    }
    let a = [
        edge_a[0] / edge_scale,
        edge_a[1] / edge_scale,
        edge_a[2] / edge_scale,
    ];
    let b = [
        edge_b[0] / edge_scale,
        edge_b[1] / edge_scale,
        edge_b[2] / edge_scale,
    ];
    let triangle_normal = cross(a, b);
    let triangle_normal_magnitude = vector_norm(triangle_normal);
    if !triangle_normal_magnitude.is_finite() || triangle_normal_magnitude <= 256.0 * f64::EPSILON {
        return Err(entity_refusal(
            source_fingerprint,
            face.face_id,
            "face_surface.bounds",
            "plane-backed triangle is numerically degenerate",
        ));
    }
    let expected_normal = if plane.same_sense {
        plane.axis
    } else {
        [-plane.axis[0], -plane.axis[1], -plane.axis[2]]
    };
    let alignment = dot(triangle_normal, expected_normal);
    let alignment_floor = 256.0 * f64::EPSILON * triangle_normal_magnitude;
    if !alignment.is_finite() || alignment <= alignment_floor {
        return Err(entity_refusal(
            source_fingerprint,
            face.face_id,
            "face_surface.same_sense",
            format!(
                "triangle winding disagrees with PLANE #{} and same_sense",
                plane.plane_id
            ),
        ));
    }

    let coordinate_scale = points.iter().fold(
        plane
            .location
            .x
            .abs()
            .max(plane.location.y.abs())
            .max(plane.location.z.abs())
            .max(1.0),
        |largest, point| {
            largest
                .max(point.x.abs())
                .max(point.y.abs())
                .max(point.z.abs())
        },
    );
    let arithmetic_tolerance = 256.0 * f64::EPSILON * coordinate_scale;
    let conversion_tolerance = coordinate_upper * 8.0;
    let tolerance = outward_nonnegative_sum(arithmetic_tolerance, conversion_tolerance)
        .ok_or_else(|| {
            entity_refusal(
                source_fingerprint,
                plane.plane_id,
                "plane.position",
                "plane-consistency tolerance is not finite",
            )
        })?;
    let mut accepted_upper = 0.0f64;
    for point in points {
        let residual =
            point_plane_residual(point, plane.location, plane.axis).ok_or_else(|| {
                entity_refusal(
                    source_fingerprint,
                    face.face_id,
                    "face_surface.face_geometry",
                    "point-to-plane residual is not finite",
                )
            })?;
        if residual > tolerance {
            return Err(entity_refusal(
                source_fingerprint,
                face.face_id,
                "face_surface.face_geometry",
                format!(
                    "triangle vertex is {residual:e} from PLANE #{}; admitted tolerance is {tolerance:e}",
                    plane.plane_id
                ),
            ));
        }
        let residual_upper =
            outward_nonnegative_sum(residual, arithmetic_tolerance).ok_or_else(|| {
                entity_refusal(
                    source_fingerprint,
                    face.face_id,
                    "face_surface.face_geometry",
                    "accepted point-to-plane estimate is not finite",
                )
            })?;
        accepted_upper = accepted_upper.max(residual_upper);
    }
    Ok(accepted_upper)
}

fn point_plane_residual(point: Point3, location: Point3, normal: [f64; 3]) -> Option<f64> {
    let delta = [
        point.x - location.x,
        point.y - location.y,
        point.z - location.z,
    ];
    let scale = delta
        .iter()
        .fold(0.0f64, |largest, value| largest.max(value.abs()));
    if !scale.is_finite() {
        return None;
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let scaled = [delta[0] / scale, delta[1] / scale, delta[2] / scale];
    let residual = dot(scaled, normal).abs() * scale;
    residual.is_finite().then_some(residual)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vector_norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn cross_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    vector_norm(cross(a, b))
}

fn outward_nonnegative_sum(a: f64, b: f64) -> Option<f64> {
    if !a.is_finite() || !b.is_finite() || a < 0.0 || b < 0.0 {
        return None;
    }
    let sum = a + b;
    if !sum.is_finite() {
        return None;
    }
    let rounded = if sum == 0.0 {
        0.0
    } else {
        f64::from_bits(sum.to_bits() + 1)
    };
    rounded.is_finite().then_some(rounded)
}

fn reserve_exact<T>(
    values: &mut Vec<T>,
    additional: usize,
    source_fingerprint: u64,
    stage: &'static str,
) -> Result<(), StepFacetedRefusal> {
    values
        .try_reserve_exact(additional)
        .map_err(|error| resource_refusal(source_fingerprint, stage, error.to_string()))
}

fn cancellable_sort_by_key<T, F>(
    values: &mut Vec<T>,
    key: F,
    cx: &Cx<'_>,
    source_fingerprint: u64,
    stage: &'static str,
) -> Result<(), StepFacetedRefusal>
where
    T: Copy,
    F: Fn(&T) -> u64,
{
    checkpoint(cx, source_fingerprint, stage)?;
    for chunk in values.chunks_mut(CANCELLATION_POLL_STRIDE) {
        chunk.sort_unstable_by_key(&key);
        checkpoint(cx, source_fingerprint, stage)?;
    }
    if values.len() <= CANCELLATION_POLL_STRIDE {
        return Ok(());
    }

    let chunk_count = values.len().div_ceil(CANCELLATION_POLL_STRIDE);
    let mut heap = BinaryHeap::new();
    heap.try_reserve_exact(chunk_count)
        .map_err(|error| resource_refusal(source_fingerprint, stage, error.to_string()))?;
    for chunk_index in 0..chunk_count {
        let index = chunk_index * CANCELLATION_POLL_STRIDE;
        heap.push(Reverse((key(&values[index]), chunk_index, index)));
    }

    let mut sorted = Vec::new();
    reserve_exact(&mut sorted, values.len(), source_fingerprint, stage)?;
    while let Some(Reverse((_entry_key, chunk_index, index))) = heap.pop() {
        poll(cx, source_fingerprint, stage, sorted.len())?;
        sorted.push(values[index]);
        let next = index + 1;
        let chunk_end = ((chunk_index + 1) * CANCELLATION_POLL_STRIDE).min(values.len());
        if next < chunk_end {
            heap.push(Reverse((key(&values[next]), chunk_index, next)));
        }
    }
    *values = sorted;
    checkpoint(cx, source_fingerprint, stage)
}

fn ensure_auxiliary_bound(
    estimated: usize,
    limits: StepFacetedLimits,
    source_fingerprint: u64,
    stage: &'static str,
) -> Result<(), StepFacetedRefusal> {
    if estimated > limits.max_auxiliary_bytes {
        return Err(resource_refusal(
            source_fingerprint,
            stage,
            format!(
                "estimated {estimated} logical element-payload bytes exceed cap {}",
                limits.max_auxiliary_bytes
            ),
        ));
    }
    Ok(())
}

fn adjacent_duplicate(
    values: &[u64],
    cx: &Cx<'_>,
    source_fingerprint: u64,
    stage: &'static str,
) -> Result<Option<u64>, StepFacetedRefusal> {
    for (offset, pair) in values.windows(2).enumerate() {
        poll(cx, source_fingerprint, stage, offset)?;
        if pair[0] == pair[1] {
            return Ok(Some(pair[0]));
        }
    }
    Ok(None)
}

fn dedup_sorted(
    values: &mut Vec<u64>,
    cx: &Cx<'_>,
    source_fingerprint: u64,
    stage: &'static str,
) -> Result<(), StepFacetedRefusal> {
    if values.len() < 2 {
        return Ok(());
    }
    let mut write = 1usize;
    for read in 1..values.len() {
        poll(cx, source_fingerprint, stage, read)?;
        let value = values[read];
        if value != values[write - 1] {
            values[write] = value;
            write += 1;
        }
    }
    values.truncate(write);
    Ok(())
}

fn checkpoint(
    cx: &Cx<'_>,
    source_fingerprint: u64,
    stage: &'static str,
) -> Result<(), StepFacetedRefusal> {
    cx.checkpoint().map_err(|_| StepFacetedRefusal::Cancelled {
        source_fingerprint,
        stage,
    })
}

fn poll(
    cx: &Cx<'_>,
    source_fingerprint: u64,
    stage: &'static str,
    offset: usize,
) -> Result<(), StepFacetedRefusal> {
    if offset % CANCELLATION_POLL_STRIDE == 0 {
        checkpoint(cx, source_fingerprint, stage)?;
    }
    Ok(())
}

fn entity_refusal(
    source_fingerprint: u64,
    instance_id: u64,
    relationship: &'static str,
    what: impl Into<String>,
) -> StepFacetedRefusal {
    StepFacetedRefusal::Entity {
        source_fingerprint,
        instance_id,
        relationship,
        what: what.into(),
    }
}

fn resource_refusal(
    source_fingerprint: u64,
    stage: &'static str,
    what: impl Into<String>,
) -> StepFacetedRefusal {
    StepFacetedRefusal::Resource {
        source_fingerprint,
        stage,
        what: what.into(),
    }
}

fn ulp_magnitude(value: f64) -> f64 {
    let magnitude = value.abs();
    if magnitude == 0.0 {
        return f64::from_bits(1);
    }
    if magnitude == f64::MAX {
        return f64::INFINITY;
    }
    f64::from_bits(magnitude.to_bits() + 1) - magnitude
}

fn semantic_fingerprint(
    profile: StepFacetedProfile,
    root_id: u64,
    shell_id: u64,
    point_ids: &[u64],
    positions: &[Point3],
    raw_faces: &[RawFace],
    triangles: &[[u32; 3]],
    source_fingerprint: u64,
    cx: &Cx<'_>,
) -> Result<u64, StepFacetedRefusal> {
    let mut hash = fnv1a64_update(FNV_OFFSET, STEP_FACETED_DECODER_VERSION.as_bytes());
    hash = fnv1a64_update(hash, profile.schema_identifier().as_bytes());
    hash = fnv1a64_update(hash, &root_id.to_le_bytes());
    hash = fnv1a64_update(hash, &shell_id.to_le_bytes());
    for (point_offset, (&point_id, point)) in point_ids.iter().zip(positions).enumerate() {
        poll(
            cx,
            source_fingerprint,
            "semantic-point-fingerprint",
            point_offset,
        )?;
        hash = fnv1a64_update(hash, &point_id.to_le_bytes());
        hash = fnv1a64_update(hash, &point.x.to_bits().to_le_bytes());
        hash = fnv1a64_update(hash, &point.y.to_bits().to_le_bytes());
        hash = fnv1a64_update(hash, &point.z.to_bits().to_le_bytes());
    }
    for (face_offset, (face, triangle)) in raw_faces.iter().zip(triangles).enumerate() {
        poll(
            cx,
            source_fingerprint,
            "semantic-face-fingerprint",
            face_offset,
        )?;
        hash = fnv1a64_update(hash, &face.face_id.to_le_bytes());
        hash = fnv1a64_update(hash, &face.bound_id.to_le_bytes());
        hash = fnv1a64_update(hash, &face.loop_id.to_le_bytes());
        for vertex in triangle {
            hash = fnv1a64_update(hash, &vertex.to_le_bytes());
        }
        match face.plane {
            None => hash = fnv1a64_update(hash, &[0]),
            Some(plane) => {
                hash = fnv1a64_update(hash, &[1]);
                hash = fnv1a64_update(hash, &plane.plane_id.to_le_bytes());
                hash = fnv1a64_update(hash, &plane.placement_id.to_le_bytes());
                hash = fnv1a64_update(hash, &plane.location_id.to_le_bytes());
                hash = fingerprint_optional_id(hash, plane.axis_id);
                hash = fingerprint_optional_id(hash, plane.ref_direction_id);
                hash = fnv1a64_update(hash, &[u8::from(plane.same_sense)]);
                for value in [plane.location.x, plane.location.y, plane.location.z] {
                    hash = fnv1a64_update(hash, &value.to_bits().to_le_bytes());
                }
                for value in plane.axis {
                    hash = fnv1a64_update(hash, &value.to_bits().to_le_bytes());
                }
                match plane.ref_direction {
                    None => hash = fnv1a64_update(hash, &[0]),
                    Some(reference) => {
                        hash = fnv1a64_update(hash, &[1]);
                        for value in reference {
                            hash = fnv1a64_update(hash, &value.to_bits().to_le_bytes());
                        }
                    }
                }
            }
        }
    }
    Ok(hash)
}

fn fingerprint_optional_id(mut hash: u64, id: Option<u64>) -> u64 {
    match id {
        None => fnv1a64_update(hash, &[0]),
        Some(id) => {
            hash = fnv1a64_update(hash, &[1]);
            fnv1a64_update(hash, &id.to_le_bytes())
        }
    }
}

fn fnv1a64_update(mut hash: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}
