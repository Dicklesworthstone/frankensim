//! Strict native decoding for a triangle-only ISO 10303 faceted-B-rep subset.
//!
//! This module interprets only one caller-selected, root-reachable chain:
//! `FACETED_BREP -> CLOSED_SHELL -> FACE -> FACE_OUTER_BOUND -> POLY_LOOP ->
//! CARTESIAN_POINT`. Every face must contain one triangular outer loop. The
//! resulting soup is still routed through [`crate::import_step_tessellation`],
//! which owns repair, bounded mesh-integrity checks, evidence composition, and
//! SDF sampling.
//! A declared AP203/AP214-family schema is an admission label only; this module
//! does not validate an EXPRESS schema, representation context, units, product
//! linkage, or application-protocol global rules. In particular, v1 refuses
//! plane-backed `FACE_SURFACE` records used by ordinary AP203/AP214 exchange;
//! those require a pinned `PLANE`/placement/direction subset rather than being
//! silently treated as bare resource `FACE` records.

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
pub const STEP_FACETED_DECODER_VERSION: &str = "step-triangular-faceted-brep-v1";
/// Stable materializer name retained by the downstream STEP import receipt.
pub const STEP_FACETED_MATERIALIZER_NAME: &str = "fs-io-native-faceted-brep";

const CANCELLATION_POLL_STRIDE: usize = 4_096;
const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Exact schema declaration admitted by the v1 resource-entity decoder.
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
    /// Maximum distinct reachable Cartesian points.
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
    reversed_bounds: usize,
    coordinate_conversion: NumericalCertificate,
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
             \"vertices\":{},\"triangles\":{},\"reversed_bounds\":{},\
             \"coordinate_conversion\":{{\"kind\":\"estimate\",\"lo\":{},\"hi\":{}}},\
             \"limits\":{{\"vertices\":{},\"triangles\":{},\"auxiliary_bytes\":{}}},\
             \"no_claim\":\"bare FACE resource subset only; no FACE_SURFACE or PLANE support; no full EXPRESS or AP conformance, unit-context correspondence, product linkage, self-intersection certificate, or topology authority\"}}",
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
            self.reversed_bounds,
            self.coordinate_conversion.lo,
            self.coordinate_conversion.hi,
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
}

/// Decode one caller-selected triangular `FACETED_BREP` closure.
///
/// Reachable instances must use the strict six-entity subset documented at
/// module scope. Unknown instances outside the selected closure are ignored.
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
    for (face_offset, &face_id) in face_ids.iter().enumerate() {
        poll(cx, source_fingerprint, "face-traversal", face_offset)?;
        let face = simple_entity(
            instances,
            &index,
            face_id,
            "FACE",
            "root.outer.faces[]",
            source_fingerprint,
        )?;
        require_shape(face, face_id, 2, "root.outer.faces[]", source_fingerprint)?;
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
        let reverse = require_boolean(
            &bound.parameters[2],
            bound_id,
            "face.outer_bound.orientation",
            source_fingerprint,
        )?;
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
    let mut coordinate_upper = 0.0f64;
    for (point_offset, &point_id) in point_ids.iter().enumerate() {
        poll(
            cx,
            source_fingerprint,
            "point-materialization",
            point_offset,
        )?;
        let point = simple_entity(
            instances,
            &index,
            point_id,
            "CARTESIAN_POINT",
            "poly_loop.polygon[]",
            source_fingerprint,
        )?;
        require_shape(point, point_id, 2, "cartesian_point", source_fingerprint)?;
        let coordinates = require_aggregate(
            &point.parameters[1],
            point_id,
            "cartesian_point.coordinates",
            source_fingerprint,
        )?;
        if coordinates.len() != 3 {
            return Err(entity_refusal(
                source_fingerprint,
                point_id,
                "cartesian_point.coordinates",
                format!(
                    "3-D faceted subset requires exactly three coordinates; got {}",
                    coordinates.len()
                ),
            ));
        }
        let mut values = [0.0f64; 3];
        let mut point_ulp = 0.0f64;
        for (axis, coordinate) in coordinates.iter().enumerate() {
            let StepValue::Number(token) = coordinate else {
                return Err(entity_refusal(
                    source_fingerprint,
                    point_id,
                    "cartesian_point.coordinates[]",
                    "coordinates must be decimal number tokens",
                ));
            };
            let value = token.parse::<f64>().map_err(|error| {
                entity_refusal(
                    source_fingerprint,
                    point_id,
                    "cartesian_point.coordinates[]",
                    format!("coordinate {token:?} does not convert to f64: {error}"),
                )
            })?;
            if !value.is_finite() {
                return Err(entity_refusal(
                    source_fingerprint,
                    point_id,
                    "cartesian_point.coordinates[]",
                    format!("coordinate {token:?} converts to non-finite f64"),
                ));
            }
            values[axis] = value;
            point_ulp = point_ulp.max(ulp_magnitude(value));
        }
        let spatial_estimate = point_ulp * 2.0;
        if !spatial_estimate.is_finite() {
            return Err(entity_refusal(
                source_fingerprint,
                point_id,
                "cartesian_point.coordinates",
                "finite coordinate has no finite conservative spatial conversion estimate",
            ));
        }
        coordinate_upper = coordinate_upper.max(spatial_estimate);
        positions.push(Point3::new(values[0], values[1], values[2]));
    }

    let mut triangles = Vec::new();
    reserve_exact(
        &mut triangles,
        raw_faces.len(),
        source_fingerprint,
        "triangle-materialization",
    )?;
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
        triangles.push(triangle);
    }

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
        reversed_bounds,
        coordinate_conversion: NumericalCertificate::estimate(0.0, coordinate_upper),
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
/// The caller supplies the unit because this v1 closure deliberately ignores
/// representation context. Decimal-coordinate conversion is retained as an
/// estimated materialization deviation. The returned decoder receipt remains
/// separate from the existing topology/SDF receipt.
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
        decoder_receipt.coordinate_conversion(),
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
    let offset = index
        .binary_search_by_key(&id, |entry| entry.0)
        .map(|position| index[position].1)
        .map_err(|_| {
            entity_refusal(
                source_fingerprint,
                id,
                relationship,
                format!("referenced {expected} instance does not exist"),
            )
        })?;
    let instance = &instances[offset];
    if instance.components.len() != 1 {
        return Err(entity_refusal(
            source_fingerprint,
            id,
            relationship,
            format!(
                "reachable {expected} must be a simple instance; found {} complex components",
                instance.components.len()
            ),
        ));
    }
    let entity = &instance.components[0];
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
        StepValue::Enumeration(value) if value == "T" => Ok(false),
        StepValue::Enumeration(value) if value == "F" => Ok(true),
        _ => Err(entity_refusal(
            source_fingerprint,
            owner_id,
            relationship,
            "orientation must be .T. or .F.",
        )),
    }
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
    }
    Ok(hash)
}

fn fnv1a64_update(mut hash: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}
