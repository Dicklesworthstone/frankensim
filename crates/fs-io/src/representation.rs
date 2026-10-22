//! Representation-aware geometry admission (IRP M6.3-W1).
//!
//! [`crate::promote`] repairs whatever it finds: it fan-fills every boundary
//! loop up to `max_hole_edges`, so a tone hole or a sound hole would be closed
//! silently. Instrument geometry instead DECLARES what it is:
//! - a closed solid (watertight 2-manifold, outward);
//! - a midsurface (open edges declared loop by loop, thickness per vertex or
//!   per triangle);
//! - an opened solid (each real opening declared).
//!
//! [`admit_representation`] checks the soup against that declaration and
//! never repairs. The declaration also supplies the length unit, the frame, the
//! reference state, and the part/material map. File labels are never units or
//! materials. The ledger row is an opaque token: fs-io (L2) does not read the
//! data ledger.
//!
//! The admitted value has positions `origin + R (s x)` in metres. Its receipt
//! records every declared fact and every count it was checked against.

use crate::quarantine::{
    CensusRefusal, ImportCensusPolicy, IntersectionInspection, Quarantined, census_with_policy,
    simple_boundary_loops,
};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_rep_mesh::Soup;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Semantic identity of [`RepresentationReceipt`].
pub const REPRESENTATION_SEMANTICS_VERSION: &str = "fs-io-representation-v1";
/// Largest admitted number of declared loops (edges or openings).
pub const MAX_DECLARED_LOOPS: usize = 256;
/// Largest admitted number of part map entries.
pub const MAX_PARTS: usize = 256;
/// Frame admission tolerance on `|det R - 1|` and `||R^T R - I||_inf`.
pub const FRAME_TOLERANCE: f64 = 1.0e-12;
/// Cancellation poll stride of the admission loops.
const POLL_STRIDE: usize = 1024;
/// Intersection pairs the embedded census may inspect before it samples.
const INTERSECTION_PAIRS: usize = 1 << 20;

/// What the source geometry represents.
#[derive(Debug, Clone, PartialEq)]
pub enum Representation {
    /// A watertight, outward-oriented 2-manifold bounding a volume.
    ClosedSolid,
    /// A zero-thickness midsurface with declared free edges and thickness.
    Midsurface {
        /// Thickness, in the declared length unit.
        thickness: ThicknessField,
        /// One declaration per free-edge loop.
        edges: Vec<DeclaredLoop>,
    },
    /// A solid with real openings (tone holes, sound holes, mouths).
    OpenedSolid {
        /// One declaration per opening.
        openings: Vec<DeclaredLoop>,
    },
}

/// A declared boundary loop, located by a point on it.
#[derive(Debug, Clone, PartialEq)]
pub struct DeclaredLoop {
    /// Caller label, unique within the declaration.
    pub label: String,
    /// A point on the loop, in source coordinates and units.
    pub probe: [f64; 3],
}

/// Midsurface thickness, in the declared length unit.
#[derive(Debug, Clone, PartialEq)]
pub enum ThicknessField {
    /// One value per vertex.
    PerVertex(Vec<f64>),
    /// One value per triangle.
    PerTriangle(Vec<f64>),
}

/// A declared length unit. It is never guessed from the file.
#[derive(Debug, Clone, PartialEq)]
pub struct LengthUnit {
    /// Unit label recorded in the receipt ("m", "mm", "in").
    pub label: String,
    /// Metres per source unit.
    pub metres_per_unit: f64,
}

/// Placement of the source frame in the assembly frame: `x = origin + R x_src`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame {
    /// Origin of the source frame [m].
    pub origin_m: [f64; 3],
    /// Rows of `R`, which must be right-handed and orthonormal.
    pub axes: [[f64; 3]; 3],
}

impl Frame {
    /// The identity placement.
    pub const IDENTITY: Self = Self {
        origin_m: [0.0; 3],
        axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };
}

/// The mechanical state the geometry describes; owned by fs-geom so imports
/// and generated parts declare the same type.
pub use fs_geom::ReferenceState;

/// A part's material binding.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialRef {
    /// A material card from a declared pack.
    Card {
        /// Pack identity.
        pack_id: String,
        /// Card identity within the pack.
        card_id: String,
    },
    /// Declared as unknown, with the reason.
    DeclaredUnknown {
        /// Why no card applies.
        reason: String,
    },
    /// A label carried by the source file (OBJ `usemtl`, group or object name).
    /// A converter forwarding it lands here. It is always refused: a file label
    /// names a region, not a material.
    SourceLabel(String),
}

/// Which triangles one part owns.
#[derive(Debug, Clone, PartialEq)]
pub enum PartSelector {
    /// Explicit triangle indices.
    Triangles(Vec<u32>),
    /// The half-open triangle range `start..end`.
    Range {
        /// First triangle.
        start: u32,
        /// One past the last triangle.
        end: u32,
    },
}

/// One part of the source geometry and its material.
#[derive(Debug, Clone, PartialEq)]
pub struct PartMap {
    /// Owned triangles.
    pub selector: PartSelector,
    /// Opaque part token (the scenario binds it).
    pub part_token: String,
    /// Material binding.
    pub material: MaterialRef,
}

/// Everything the caller declares about one geometry source.
#[derive(Debug, Clone, PartialEq)]
pub struct RepresentationDeclaration {
    /// What the source represents.
    pub representation: Representation,
    /// Source length unit; `None` refuses.
    pub unit: Option<LengthUnit>,
    /// Placement in the assembly frame.
    pub frame: Frame,
    /// Reference state; `None` refuses.
    pub reference_state: Option<ReferenceState>,
    /// Part/material map covering every triangle exactly once.
    pub parts: Vec<PartMap>,
    /// Data-ledger row id (opaque here; the ledger checker resolves it).
    pub ledger_row: String,
    /// Probe-to-loop matching and census tolerance [m].
    pub tolerance_m: f64,
}

/// Typed refusals of [`admit_representation`].
#[derive(Debug, Clone, PartialEq)]
pub enum RepresentationRefusal {
    /// A closed solid has open edges.
    NotClosed {
        /// Edges used by exactly one triangle.
        open_edges: usize,
    },
    /// An edge is used by three or more triangles.
    NonManifold {
        /// Such edges.
        edges: usize,
    },
    /// Open edges that do not form simple loops.
    NonSimpleBoundary {
        /// Edges used by exactly one triangle.
        open_edges: usize,
        /// Of those, edges on simple loops.
        looped_edges: usize,
    },
    /// Adjacent triangles disagree on orientation.
    Orientation {
        /// Manifold edges traversed twice in the same direction.
        inconsistent_edges: usize,
    },
    /// A closed solid encloses non-positive signed volume (inside out).
    InsideOut {
        /// Signed volume [m^3].
        signed_volume_m3: f64,
    },
    /// The embedded census found a blocking defect; nothing is repaired here.
    Defect {
        /// Census class.
        class: &'static str,
        /// Instances.
        count: usize,
    },
    /// A boundary loop has no declaration.
    UndeclaredLoop {
        /// The loop's smallest vertex index.
        first_vertex: u32,
        /// Label of the nearest declared probe, if any.
        nearest_probe: Option<String>,
    },
    /// A declaration matches no boundary loop within tolerance.
    MissingDeclaredLoop {
        /// Its label.
        label: String,
    },
    /// Two declarations match the same loop, or a label repeats.
    AmbiguousLoop {
        /// The second declaration's label.
        label: String,
    },
    /// The thickness array does not have one value per vertex/triangle.
    ThicknessLength {
        /// Required length.
        expected: usize,
        /// Supplied length.
        got: usize,
    },
    /// A thickness value is not finite and positive.
    Thickness {
        /// Its index.
        index: usize,
        /// Its value.
        value: f64,
    },
    /// No (or an unusable) length unit was declared.
    UnitMissing,
    /// The frame is not right-handed orthonormal.
    Frame {
        /// `det R`.
        det: f64,
        /// `||R^T R - I||_inf`.
        orthogonality: f64,
    },
    /// Triangles no part owns.
    Unmapped {
        /// Count.
        triangles: usize,
    },
    /// Triangles owned by more than one part.
    DoublyMapped {
        /// Count.
        triangles: usize,
    },
    /// A part selector names a triangle outside the mesh.
    SelectorOutOfRange {
        /// Part token.
        part: String,
    },
    /// A material was taken from a file label.
    MaterialFromLabel {
        /// The label.
        label: String,
    },
    /// A required text field is empty.
    MissingText {
        /// Which field.
        field: &'static str,
    },
    /// No reference state was declared.
    ReferenceStateMissing,
    /// A declared count exceeds its cap.
    Cap {
        /// What was counted.
        what: &'static str,
        /// The cap.
        cap: usize,
        /// The count.
        got: usize,
    },
    /// The embedded census refused or was cancelled.
    Census(CensusRefusal),
}

impl core::fmt::Display for RepresentationRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "representation admission refused: {self:?}")
    }
}

impl core::error::Error for RepresentationRefusal {}

/// One declared loop as matched.
#[derive(Debug, Clone, PartialEq)]
pub struct MatchedLoop {
    /// Declaration label.
    pub label: String,
    /// The loop's vertex indices, ascending.
    pub vertices: Vec<u32>,
    /// Probe-to-loop distance [m].
    pub probe_distance_m: f64,
}

/// What [`admit_representation`] checked and admitted.
#[derive(Debug, Clone, PartialEq)]
pub struct RepresentationReceipt {
    /// [`REPRESENTATION_SEMANTICS_VERSION`].
    pub semantics_version: &'static str,
    /// "closed-solid", "midsurface" or "opened-solid".
    pub representation: &'static str,
    /// Data-ledger row id.
    pub ledger_row: String,
    /// Source format and hash, from the quarantine receipt.
    pub source_format: &'static str,
    /// Source content hash.
    pub source_hash: u64,
    /// Declared unit label.
    pub unit: String,
    /// Metres per source unit, applied exactly once.
    pub metres_per_unit: f64,
    /// Frame origin [m].
    pub frame_origin_m: [f64; 3],
    /// `det R`.
    pub frame_det: f64,
    /// `||R^T R - I||_inf`.
    pub frame_orthogonality: f64,
    /// Declared reference state.
    pub reference_state: &'static str,
    /// Vertex and triangle counts.
    pub vertices: usize,
    /// Triangle count.
    pub triangles: usize,
    /// Simple boundary loops found.
    pub boundary_loops: usize,
    /// Declared loops, as matched.
    pub declared_loops: Vec<MatchedLoop>,
    /// Edges used by exactly one triangle.
    pub open_edges: usize,
    /// Edges used by three or more triangles.
    pub nonmanifold_edges: usize,
    /// Repair actions applied (always zero: admission never repairs).
    pub repairs: usize,
    /// Signed enclosed volume [m^3] (closed solids only).
    pub signed_volume_m3: Option<f64>,
    /// Thickness range [m] (midsurfaces only).
    pub thickness_range_m: Option<(f64, f64)>,
    /// Thickness layout: "per-vertex" or "per-triangle" (midsurfaces only).
    pub thickness_layout: Option<&'static str>,
    /// Part tokens in declaration order.
    pub parts: Vec<String>,
}

impl RepresentationReceipt {
    /// Canonical JSON (one line per admitted fixture).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = format!(
            "{{\"semantics\":\"{}\",\"representation\":\"{}\",\"ledger_row\":\"{}\",\
             \"source_format\":\"{}\",\"source_hash\":\"{:016x}\",\"unit\":\"{}\",\
             \"metres_per_unit\":{:e},\"frame_origin_m\":[{:e},{:e},{:e}],\"frame_det\":{:e},\
             \"frame_orthogonality\":{:e},\"reference_state\":\"{}\",\"vertices\":{},\
             \"triangles\":{},\"boundary_loops\":{},\"declared_loops\":[",
            self.semantics_version,
            self.representation,
            escape(&self.ledger_row),
            self.source_format,
            self.source_hash,
            escape(&self.unit),
            self.metres_per_unit,
            self.frame_origin_m[0],
            self.frame_origin_m[1],
            self.frame_origin_m[2],
            self.frame_det,
            self.frame_orthogonality,
            self.reference_state,
            self.vertices,
            self.triangles,
            self.boundary_loops,
        );
        for (i, matched) in self.declared_loops.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(
                s,
                "{{\"label\":\"{}\",\"vertices\":{},\"probe_distance_m\":{:e}}}",
                escape(&matched.label),
                matched.vertices.len(),
                matched.probe_distance_m
            );
        }
        let _ = write!(
            s,
            "],\"open_edges\":{},\"nonmanifold_edges\":{},\"repairs\":{},\"signed_volume_m3\":{},\
             \"thickness_layout\":{},\"thickness_range_m\":{},\"parts\":[",
            self.open_edges,
            self.nonmanifold_edges,
            self.repairs,
            self.signed_volume_m3
                .map_or_else(|| "null".to_string(), |v| format!("{v:e}")),
            self.thickness_layout
                .map_or_else(|| "null".to_string(), |v| format!("\"{v}\"")),
            self.thickness_range_m
                .map_or_else(|| "null".to_string(), |(lo, hi)| format!("[{lo:e},{hi:e}]")),
        );
        for (i, part) in self.parts.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "\"{}\"", escape(part));
        }
        s.push_str("]}");
        s
    }
}

fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out
}

fn cancelled(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), RepresentationRefusal> {
    if at.is_multiple_of(POLL_STRIDE) {
        cx.checkpoint().map_err(|_| {
            RepresentationRefusal::Census(CensusRefusal {
                stage,
                at,
                reason: "cancellation requested before admission".to_string(),
                cancelled: true,
            })
        })?;
    }
    Ok(())
}

fn require_text(field: &'static str, text: &str) -> Result<(), RepresentationRefusal> {
    if text.trim().is_empty() {
        Err(RepresentationRefusal::MissingText { field })
    } else {
        Ok(())
    }
}

/// `(det R, ||R^T R - I||_inf)`.
fn frame_quality(axes: &[[f64; 3]; 3]) -> (f64, f64) {
    let r = axes;
    let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    let mut worst = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let gram: f64 = (0..3).map(|k| r[k][i] * r[k][j]).sum();
            let identity = if i == j { 1.0 } else { 0.0 };
            worst = worst.max((gram - identity).abs());
        }
    }
    (det, worst)
}

fn point_segment_distance(p: [f64; 3], a: Point3, b: Point3) -> f64 {
    let ab = [b.x - a.x, b.y - a.y, b.z - a.z];
    let ap = [p[0] - a.x, p[1] - a.y, p[2] - a.z];
    let len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let t = if len2 > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let d = [ap[0] - t * ab[0], ap[1] - t * ab[1], ap[2] - t * ab[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// Census classes that block admission (it never repairs). A near-loop gap is
/// an exporter crack (near-duplicate vertices): refused, never welded.
const BLOCKING_CLASSES: [&str; 6] = [
    "invalid-face-index",
    "non-finite-vertex",
    "degenerate-face",
    "duplicate-face",
    "unreferenced-vertex",
    "near-boundary-loop-gap",
];

/// Admit a quarantined soup against its representation declaration.
///
/// Checks run in this order, each refusing by name:
/// 1. declared facts (unit, frame, reference state, ledger row, caps, parts,
///    materials, thickness);
/// 2. the embedded census (blocking classes; no repair);
/// 3. edge use, orientation, boundary loops vs declarations;
/// 4. closed-solid signed volume.
///
/// # Errors
/// [`RepresentationRefusal`]; nothing is published on any refusal, including
/// cancellation.
#[allow(clippy::too_many_lines)] // one admission checklist, checked in its stated order
pub fn admit_representation(
    q: Quarantined<Soup>,
    d: &RepresentationDeclaration,
    cx: &Cx<'_>,
) -> Result<(Evidence<Soup>, RepresentationReceipt), RepresentationRefusal> {
    // 1. Declared facts.
    let unit = match &d.unit {
        Some(unit)
            if !unit.label.trim().is_empty()
                && unit.metres_per_unit.is_finite()
                && unit.metres_per_unit > 0.0 =>
        {
            unit
        }
        _ => return Err(RepresentationRefusal::UnitMissing),
    };
    let (det, orthogonality) = frame_quality(&d.frame.axes);
    if !((det - 1.0).abs() <= FRAME_TOLERANCE && orthogonality <= FRAME_TOLERANCE)
        || d.frame.origin_m.iter().any(|v| !v.is_finite())
    {
        return Err(RepresentationRefusal::Frame { det, orthogonality });
    }
    let reference_state = d
        .reference_state
        .ok_or(RepresentationRefusal::ReferenceStateMissing)?;
    require_text("ledger_row", &d.ledger_row)?;
    if !(d.tolerance_m.is_finite() && d.tolerance_m > 0.0) {
        return Err(RepresentationRefusal::MissingText {
            field: "tolerance_m (finite and positive)",
        });
    }
    let (kind, declared): (&'static str, &[DeclaredLoop]) = match &d.representation {
        Representation::ClosedSolid => ("closed-solid", &[]),
        Representation::Midsurface { edges, .. } => ("midsurface", edges),
        Representation::OpenedSolid { openings } => ("opened-solid", openings),
    };
    if declared.len() > MAX_DECLARED_LOOPS {
        return Err(RepresentationRefusal::Cap {
            what: "declared_loops",
            cap: MAX_DECLARED_LOOPS,
            got: declared.len(),
        });
    }
    if d.parts.len() > MAX_PARTS {
        return Err(RepresentationRefusal::Cap {
            what: "parts",
            cap: MAX_PARTS,
            got: d.parts.len(),
        });
    }
    for (i, loop_) in declared.iter().enumerate() {
        require_text("declared_loop.label", &loop_.label)?;
        if declared[..i].iter().any(|other| other.label == loop_.label) {
            return Err(RepresentationRefusal::AmbiguousLoop {
                label: loop_.label.clone(),
            });
        }
        if loop_.probe.iter().any(|v| !v.is_finite()) {
            return Err(RepresentationRefusal::MissingDeclaredLoop {
                label: loop_.label.clone(),
            });
        }
    }
    let soup = &q.raw;
    let n_vertices = soup.positions.len();
    let n_triangles = soup.triangles.len();

    // Parts and materials: every triangle exactly once; no material from a label.
    let mut owners = vec![0u32; n_triangles];
    for part in &d.parts {
        require_text("part_token", &part.part_token)?;
        match &part.material {
            MaterialRef::Card { pack_id, card_id } => {
                require_text("material.pack_id", pack_id)?;
                require_text("material.card_id", card_id)?;
            }
            MaterialRef::DeclaredUnknown { reason } => require_text("material.reason", reason)?,
            MaterialRef::SourceLabel(label) => {
                return Err(RepresentationRefusal::MaterialFromLabel {
                    label: label.clone(),
                });
            }
        }
        let out_of_range = || RepresentationRefusal::SelectorOutOfRange {
            part: part.part_token.clone(),
        };
        match &part.selector {
            PartSelector::Triangles(list) => {
                for (i, &t) in list.iter().enumerate() {
                    cancelled(cx, "part-map", i)?;
                    let slot = owners.get_mut(t as usize).ok_or_else(out_of_range)?;
                    *slot = slot.saturating_add(1);
                }
            }
            PartSelector::Range { start, end } => {
                if start > end || *end as usize > n_triangles {
                    return Err(out_of_range());
                }
                for (i, slot) in owners[*start as usize..*end as usize]
                    .iter_mut()
                    .enumerate()
                {
                    cancelled(cx, "part-map", i)?;
                    *slot = slot.saturating_add(1);
                }
            }
        }
    }
    let unmapped = owners.iter().filter(|&&o| o == 0).count();
    if unmapped > 0 {
        return Err(RepresentationRefusal::Unmapped {
            triangles: unmapped,
        });
    }
    let doubly = owners.iter().filter(|&&o| o > 1).count();
    if doubly > 0 {
        return Err(RepresentationRefusal::DoublyMapped { triangles: doubly });
    }

    // Thickness: one finite positive value per vertex or per triangle.
    let mut thickness_range_m = None;
    let mut thickness_layout = None;
    if let Representation::Midsurface { thickness, .. } = &d.representation {
        let (values, expected, layout) = match thickness {
            ThicknessField::PerVertex(v) => (v, n_vertices, "per-vertex"),
            ThicknessField::PerTriangle(v) => (v, n_triangles, "per-triangle"),
        };
        if values.len() != expected {
            return Err(RepresentationRefusal::ThicknessLength {
                expected,
                got: values.len(),
            });
        }
        let mut range = (f64::INFINITY, 0.0_f64);
        for (index, &value) in values.iter().enumerate() {
            let metres = value * unit.metres_per_unit;
            if !(value.is_finite() && value > 0.0 && metres.is_finite() && metres > 0.0) {
                return Err(RepresentationRefusal::Thickness { index, value });
            }
            range = (range.0.min(metres), range.1.max(metres));
        }
        thickness_range_m = Some(range);
        thickness_layout = Some(layout);
    }

    // 2. Embedded census: blocking classes refuse; nothing is repaired.
    let policy = ImportCensusPolicy::try_new(
        d.tolerance_m / unit.metres_per_unit,
        IntersectionInspection::DeterministicSampleF64 {
            sample_count: INTERSECTION_PAIRS,
        },
        POLL_STRIDE,
    )
    .map_err(RepresentationRefusal::Census)?;
    let census = census_with_policy(soup, policy, cx).map_err(RepresentationRefusal::Census)?;
    for class in BLOCKING_CLASSES {
        let count = census.count(class);
        if count > 0 {
            return Err(RepresentationRefusal::Defect { class, count });
        }
    }

    // 3. Edge use (split the census's lumped open/non-manifold class) and orientation.
    let mut directed: BTreeMap<(u32, u32), (u32, u32)> = BTreeMap::new();
    for (face, t) in soup.triangles.iter().enumerate() {
        cancelled(cx, "edge-use", face)?;
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let entry = directed.entry((a.min(b), a.max(b))).or_insert((0, 0));
            if a < b {
                entry.0 += 1;
            } else {
                entry.1 += 1;
            }
        }
    }
    let mut open_edges = Vec::new();
    let mut nonmanifold_edges = 0usize;
    let mut inconsistent_edges = 0usize;
    for (&edge, &(forward, backward)) in &directed {
        match forward + backward {
            1 => open_edges.push(edge),
            2 if forward != 1 => inconsistent_edges += 1,
            2 => {}
            _ => nonmanifold_edges += 1,
        }
    }
    if nonmanifold_edges > 0 {
        return Err(RepresentationRefusal::NonManifold {
            edges: nonmanifold_edges,
        });
    }
    if inconsistent_edges > 0 {
        return Err(RepresentationRefusal::Orientation { inconsistent_edges });
    }
    if matches!(d.representation, Representation::ClosedSolid) && !open_edges.is_empty() {
        return Err(RepresentationRefusal::NotClosed {
            open_edges: open_edges.len(),
        });
    }
    let loops =
        simple_boundary_loops(soup, POLL_STRIDE, cx).map_err(RepresentationRefusal::Census)?;
    let looped_edges: usize = loops.iter().map(Vec::len).sum();
    if looped_edges != open_edges.len() {
        return Err(RepresentationRefusal::NonSimpleBoundary {
            open_edges: open_edges.len(),
            looped_edges,
        });
    }

    // Match each declaration to the loop nearest its probe; every loop exactly once.
    let mut loop_of_vertex: BTreeMap<u32, usize> = BTreeMap::new();
    for (index, vertices) in loops.iter().enumerate() {
        for &v in vertices {
            loop_of_vertex.insert(v, index);
        }
    }
    let probe_distances = |probe: [f64; 3]| -> Vec<f64> {
        let mut best = vec![f64::INFINITY; loops.len()];
        for &(a, b) in &open_edges {
            if let Some(&index) = loop_of_vertex.get(&a) {
                let distance = point_segment_distance(
                    probe,
                    soup.positions[a as usize],
                    soup.positions[b as usize],
                ) * unit.metres_per_unit;
                best[index] = best[index].min(distance);
            }
        }
        best
    };
    let mut claimed: Vec<Option<usize>> = vec![None; loops.len()];
    let mut matched = Vec::with_capacity(declared.len());
    let mut all_distances = Vec::with_capacity(declared.len());
    for (i, loop_) in declared.iter().enumerate() {
        cancelled(cx, "loop-match", i)?;
        let distances = probe_distances(loop_.probe);
        let nearest = distances
            .iter()
            .enumerate()
            .filter(|&(_, &dist)| dist <= d.tolerance_m)
            .min_by(|a, b| a.1.total_cmp(b.1));
        let Some((index, &distance)) = nearest else {
            return Err(RepresentationRefusal::MissingDeclaredLoop {
                label: loop_.label.clone(),
            });
        };
        if claimed[index].is_some() {
            return Err(RepresentationRefusal::AmbiguousLoop {
                label: loop_.label.clone(),
            });
        }
        claimed[index] = Some(i);
        matched.push(MatchedLoop {
            label: loop_.label.clone(),
            vertices: loops[index].clone(),
            probe_distance_m: distance,
        });
        all_distances.push(distances);
    }
    if let Some(index) = claimed.iter().position(Option::is_none) {
        let nearest_probe = all_distances
            .iter()
            .enumerate()
            .min_by(|a, b| a.1[index].total_cmp(&b.1[index]))
            .map(|(i, _)| declared[i].label.clone());
        return Err(RepresentationRefusal::UndeclaredLoop {
            first_vertex: loops[index][0],
            nearest_probe,
        });
    }

    // Place in metres: x = origin + R (s x_src), applied exactly once.
    let s = unit.metres_per_unit;
    let r = d.frame.axes;
    let o = d.frame.origin_m;
    let positions: Vec<Point3> = soup
        .positions
        .iter()
        .map(|p| {
            let x = [p.x * s, p.y * s, p.z * s];
            let row = |k: usize| o[k] + (r[k][0] * x[0] + r[k][1] * x[1] + r[k][2] * x[2]);
            Point3::new(row(0), row(1), row(2))
        })
        .collect();

    // 4. Closed solids enclose positive volume.
    let mut signed_volume_m3 = None;
    if matches!(d.representation, Representation::ClosedSolid) {
        let volume: f64 = soup
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| positions[i as usize]);
                let cross = [
                    b.y * c.z - b.z * c.y,
                    b.z * c.x - b.x * c.z,
                    b.x * c.y - b.y * c.x,
                ];
                (a.x * cross[0] + a.y * cross[1] + a.z * cross[2]) / 6.0
            })
            .sum();
        if !(volume.is_finite() && volume > 0.0) {
            return Err(RepresentationRefusal::InsideOut {
                signed_volume_m3: volume,
            });
        }
        signed_volume_m3 = Some(volume);
    }

    let receipt = RepresentationReceipt {
        semantics_version: REPRESENTATION_SEMANTICS_VERSION,
        representation: kind,
        ledger_row: d.ledger_row.clone(),
        source_format: q.source_receipt.format,
        source_hash: q.source_receipt.source_hash,
        unit: unit.label.clone(),
        metres_per_unit: s,
        frame_origin_m: o,
        frame_det: det,
        frame_orthogonality: orthogonality,
        reference_state: reference_state.name(),
        vertices: n_vertices,
        triangles: n_triangles,
        boundary_loops: loops.len(),
        declared_loops: matched,
        open_edges: open_edges.len(),
        nonmanifold_edges,
        repairs: 0,
        signed_volume_m3,
        thickness_range_m,
        thickness_layout,
        parts: d.parts.iter().map(|p| p.part_token.clone()).collect(),
    };
    let json = receipt.to_json();
    #[allow(clippy::cast_precision_loss)]
    let qoi = n_triangles as f64;
    let triangles = q.raw.triangles;
    Ok((
        Evidence {
            value: Soup {
                positions,
                triangles,
            },
            qoi,
            numerical: NumericalCertificate::exact(qoi),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::of_bytes(json.as_bytes()),
            adjoint_ref: None,
        },
        receipt,
    ))
}
