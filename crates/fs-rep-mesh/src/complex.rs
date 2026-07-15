//! Oriented cell complexes (plan §7.2): genuine 2-D triangle complexes and
//! tet element storage with complete signed incidence chains — and the δδ = 0 identity
//! verified EXACTLY (integer arithmetic), because fs-feec's discrete
//! exterior calculus is only as good as the orientation bookkeeping here
//! (plan §8.1: "the discrete de Rham sequence is exact... purely
//! combinatorial"). Hex elements get SoA storage now; their incidence
//! joins when fs-feec's tensor-product families land (CONTRACT no-claim).

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EntityId, Field, FieldSpec,
    NeverCancel, StrongIdentity, WireType,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// An oriented tetrahedral complex with derived, canonically-ordered edge
/// and face tables and signed incidence operators d0, d1, d2.
#[derive(Debug, Clone)]
pub struct TetComplex {
    /// Vertex count.
    pub vertex_count: usize,
    /// Tets as vertex quadruples (orientation = even permutations of the
    /// stored order).
    pub tets: Vec<[u32; 4]>,
    /// Canonical edges (sorted vertex pairs), sorted.
    pub edges: Vec<[u32; 2]>,
    /// Canonical faces (sorted vertex triples), sorted.
    pub faces: Vec<[u32; 3]>,
    edge_index: BTreeMap<[u32; 2], usize>,
    face_index: BTreeMap<[u32; 3], usize>,
}

/// A signed sparse incidence operator: rows of `(column, ±1)` pairs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Incidence {
    /// `rows[r]` lists `(column, sign)` entries of row r.
    pub rows: Vec<Vec<(usize, i8)>>,
    /// Column count.
    pub cols: usize,
}

impl Incidence {
    /// Apply to an integer cochain (exact arithmetic — this is the δδ = 0
    /// verification path).
    #[must_use]
    pub fn apply(&self, x: &[i64]) -> Vec<i64> {
        self.rows
            .iter()
            .map(|row| row.iter().map(|&(c, s)| i64::from(s) * x[c]).sum())
            .collect()
    }
}

impl TetComplex {
    /// Build from tets (vertex indices; orientation as stored).
    #[must_use]
    pub fn from_tets(vertex_count: usize, tets: Vec<[u32; 4]>) -> Self {
        let mut edge_set: BTreeMap<[u32; 2], usize> = BTreeMap::new();
        let mut face_set: BTreeMap<[u32; 3], usize> = BTreeMap::new();
        for t in &tets {
            for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                let mut e = [t[i], t[j]];
                e.sort_unstable();
                let next = edge_set.len();
                edge_set.entry(e).or_insert(next);
            }
            for skip in 0..4 {
                let mut f = [0u32; 3];
                let mut k = 0;
                for (v, &vert) in t.iter().enumerate() {
                    if v != skip {
                        f[k] = vert;
                        k += 1;
                    }
                }
                f.sort_unstable();
                let next = face_set.len();
                face_set.entry(f).or_insert(next);
            }
        }
        // Re-key to sorted order (BTreeMap iteration is sorted; indices
        // follow it for determinism).
        let edges: Vec<[u32; 2]> = edge_set.keys().copied().collect();
        let faces: Vec<[u32; 3]> = face_set.keys().copied().collect();
        let edge_index = edges.iter().enumerate().map(|(i, &e)| (e, i)).collect();
        let face_index = faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();
        TetComplex {
            vertex_count,
            tets,
            edges,
            faces,
            edge_index,
            face_index,
        }
    }

    /// d0: vertices → edges (gradient shape). Row per edge `[a, b]` with
    /// a < b: `x(b) − x(a)`.
    #[must_use]
    pub fn d0(&self) -> Incidence {
        Incidence {
            rows: self
                .edges
                .iter()
                .map(|&[a, b]| vec![(a as usize, -1i8), (b as usize, 1i8)])
                .collect(),
            cols: self.vertex_count,
        }
    }

    /// d1: edges → faces (curl shape). Face `[a, b, c]` with a < b < c
    /// has boundary edges (a,b) +, (b,c) +, (a,c) −.
    #[must_use]
    pub fn d1(&self) -> Incidence {
        Incidence {
            rows: self
                .faces
                .iter()
                .map(|&[a, b, c]| {
                    vec![
                        (self.edge_index[&[a, b]], 1i8),
                        (self.edge_index[&[b, c]], 1i8),
                        (self.edge_index[&[a, c]], -1i8),
                    ]
                })
                .collect(),
            cols: self.edges.len(),
        }
    }

    /// d2: faces → tets (divergence shape). Tet `[v0, v1, v2, v3]`
    /// (canonically sorted per row) has face `omit i` with sign (−1)^i.
    #[must_use]
    pub fn d2(&self) -> Incidence {
        Incidence {
            rows: self
                .tets
                .iter()
                .map(|t| {
                    let mut sorted = *t;
                    sorted.sort_unstable();
                    (0..4)
                        .map(|skip| {
                            let mut f = [0u32; 3];
                            let mut k = 0;
                            for (v, &vert) in sorted.iter().enumerate() {
                                if v != skip {
                                    f[k] = vert;
                                    k += 1;
                                }
                            }
                            let sign = if skip % 2 == 0 { 1i8 } else { -1i8 };
                            (self.face_index[&f], sign)
                        })
                        .collect()
                })
                .collect(),
            cols: self.faces.len(),
        }
    }
}

/// Typed identity schema for one mesh lineage supplied by a machine graph or
/// another durable owner.
pub enum TriComplex2LineageSchema {}

impl CanonicalSchema for TriComplex2LineageSchema {
    const DOMAIN: &'static str = "org.frankensim.fs-rep-mesh.tri-complex2-lineage.v1";
    const NAME: &'static str = "tri-complex2-lineage";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str =
        "durable lineage namespace shared by TriComplex2 revisions and machine-IR crosswalks";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("namespace", WireType::Utf8)];
}

/// Typed identity schema for a vertex, edge, or face within one durable mesh
/// lineage.
pub enum TriFeatureSchema {}

/// bead sj31i.52.10: the feature schema admits EXACTLY a
/// `TriComplex2LineageSchema` entity as its lineage child — any other
/// role/schema/version refuses at encode time, and this binding is part
/// of `TriFeatureSchema`'s own schema identity.
static TRI_FEATURE_LINEAGE_CHILD: fs_blake3::identity::ChildSpec =
    fs_blake3::identity::ChildSpec::for_identity::<
        fs_blake3::identity::EntityId<TriComplex2LineageSchema>,
    >();

impl CanonicalSchema for TriFeatureSchema {
    const DOMAIN: &'static str = "org.frankensim.fs-rep-mesh.tri-feature.v1";
    const NAME: &'static str = "tri-complex2-feature";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str =
        "topological feature identity from lineage, dimension, and canonical vertex keys";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("lineage", &TRI_FEATURE_LINEAGE_CHILD),
        FieldSpec::required("topological-dimension", WireType::U64),
        FieldSpec::required("vertex-keys", WireType::CanonicalSet),
    ];
}

const TRI_IDENTITY_LIMITS: CanonicalLimits = CanonicalLimits::new(4096, 1024, 8, 8, 256);

/// Durable identity of a related family of [`TriComplex2`] revisions.
///
/// Equality is a typed digest comparison, not an authority or provenance
/// claim. Callers retain the namespace-to-machine-entity crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TriComplex2LineageId(EntityId<TriComplex2LineageSchema>);

impl TriComplex2LineageId {
    /// Strictly parse retained digest bytes. Parsing does not authenticate the
    /// lineage.
    #[must_use]
    pub fn parse_hex(value: &str) -> Option<Self> {
        EntityId::parse_hex(value).map(Self)
    }

    /// Exact typed digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

/// Stable typed identity of one topological feature within a lineage.
///
/// Feature identity is independent of storage order, embedding coordinates,
/// and face orientation. It moves when the lineage, topological dimension, or
/// canonical set of caller-owned vertex keys moves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TriFeatureId(EntityId<TriFeatureSchema>);

impl TriFeatureId {
    /// Strictly parse retained digest bytes. Parsing does not authenticate the
    /// feature or prove that it belongs to a live complex.
    #[must_use]
    pub fn parse_hex(value: &str) -> Option<Self> {
        EntityId::parse_hex(value).map(Self)
    }

    /// Exact typed digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

/// Mint a typed mesh-lineage identity from an exact caller-owned namespace.
///
/// No Unicode or path normalization is performed. The namespace is semantic
/// data and must remain stable across refinements that preserve feature
/// lineage.
pub fn tri_complex2_lineage_id(namespace: &str) -> Result<TriComplex2LineageId, TriComplex2Error> {
    if namespace.is_empty() {
        return Err(TriComplex2Error::EmptyLineageNamespace);
    }
    CanonicalEncoder::<EntityId<TriComplex2LineageSchema>, _>::new(TRI_IDENTITY_LIMITS, NeverCancel)
        .and_then(|encoder| encoder.utf8(Field::new(0, "namespace"), namespace))
        .and_then(CanonicalEncoder::finish)
        .map(|receipt| TriComplex2LineageId(receipt.id()))
        .map_err(TriComplex2Error::Identity)
}

fn tri_feature_id(
    lineage: TriComplex2LineageId,
    topological_dimension: u64,
    vertex_keys: &[u64],
) -> Result<TriFeatureId, TriComplex2Error> {
    let mut keys: Vec<[u8; 8]> = vertex_keys.iter().map(|key| key.to_le_bytes()).collect();
    keys.sort_unstable();
    CanonicalEncoder::<EntityId<TriFeatureSchema>, _>::new(TRI_IDENTITY_LIMITS, NeverCancel)
        .and_then(|encoder| encoder.child(Field::new(0, "lineage"), lineage.0))
        .and_then(|encoder| {
            encoder.u64(
                Field::new(1, "topological-dimension"),
                topological_dimension,
            )
        })
        .and_then(|encoder| {
            encoder.canonical_set(
                Field::new(2, "vertex-keys"),
                keys.len() as u64,
                keys.iter().map(|key| &key[..]),
            )
        })
        .and_then(CanonicalEncoder::finish)
        .map(|receipt| TriFeatureId(receipt.id()))
        .map_err(TriComplex2Error::Identity)
}

/// Coordinate measure attached to a two-dimensional topological complex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Metric2 {
    /// Cartesian planar measure. Face measure is area times `thickness`.
    Planar {
        /// Out-of-plane thickness.
        thickness: f64,
    },
    /// Axisymmetric meridian measure. Coordinate zero is radius and face
    /// measure is `angular_span * integral(radius dA)`.
    Axisymmetric {
        /// Swept angle in radians, in `(0, 2π]`.
        angular_span: f64,
    },
}

impl Metric2 {
    /// Construct a finite positive planar-thickness measure.
    pub fn planar(thickness: f64) -> Result<Self, Metric2Error> {
        validate_positive_metric("thickness", thickness)?;
        Ok(Self::Planar { thickness })
    }

    /// Construct a finite positive axisymmetric sweep of at most one turn.
    pub fn axisymmetric(angular_span: f64) -> Result<Self, Metric2Error> {
        validate_positive_metric("angular span", angular_span)?;
        if angular_span > core::f64::consts::TAU {
            return Err(Metric2Error::AngularSpanExceedsTurn {
                bits: angular_span.to_bits(),
            });
        }
        Ok(Self::Axisymmetric { angular_span })
    }

    /// True for the axisymmetric weighted measure.
    #[must_use]
    pub const fn is_axisymmetric(self) -> bool {
        matches!(self, Self::Axisymmetric { .. })
    }
}

fn validate_positive_metric(name: &'static str, value: f64) -> Result<(), Metric2Error> {
    if !value.is_finite() {
        return Err(Metric2Error::NonFinite {
            name,
            bits: value.to_bits(),
        });
    }
    if value <= 0.0 {
        return Err(Metric2Error::NonPositive {
            name,
            bits: value.to_bits(),
        });
    }
    Ok(())
}

/// Refusal from metric-metadata construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Metric2Error {
    /// A scale was NaN or infinite.
    NonFinite {
        /// Scale name.
        name: &'static str,
        /// Exact refused bits.
        bits: u64,
    },
    /// A scale was zero or negative.
    NonPositive {
        /// Scale name.
        name: &'static str,
        /// Exact refused bits.
        bits: u64,
    },
    /// An axisymmetric sweep exceeded one turn.
    AngularSpanExceedsTurn {
        /// Exact refused bits.
        bits: u64,
    },
}

impl fmt::Display for Metric2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite { name, bits } => {
                write!(f, "{name} is non-finite (bits 0x{bits:016x})")
            }
            Self::NonPositive { name, bits } => {
                write!(f, "{name} is not positive (bits 0x{bits:016x})")
            }
            Self::AngularSpanExceedsTurn { bits } => write!(
                f,
                "axisymmetric angular span exceeds one turn (bits 0x{bits:016x})"
            ),
        }
    }
}

impl core::error::Error for Metric2Error {}

/// One oriented edge in a trace subcomplex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TraceEdge2 {
    /// Index in the parent complex's canonical edge table.
    pub global_edge: usize,
    /// Parent face whose selected-side orientation induced this trace edge.
    pub source_face: usize,
    /// Global vertices in the selected face's boundary orientation.
    pub oriented_vertices: [u32; 2],
    /// Typed identity of the parent edge.
    pub feature_id: TriFeatureId,
}

/// Deterministic boundary/trace extraction for a selected face subcomplex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraceMap2 {
    /// Canonically sorted parent vertex indices used by the trace.
    pub vertices: Vec<u32>,
    /// Parent edges on the selected subcomplex's boundary, in canonical parent
    /// edge order and oriented by the selected side.
    pub edges: Vec<TraceEdge2>,
    /// Exact trace incidence from trace-local vertices to trace edges.
    pub d0: Incidence,
}

/// A genuine oriented two-dimensional triangle cell complex.
#[derive(Debug, Clone)]
pub struct TriComplex2 {
    lineage: TriComplex2LineageId,
    metric: Metric2,
    /// Embedding coordinates `(x, y)` for planar measure or `(radius, z)` for
    /// axisymmetric measure.
    pub vertices: Vec<[f64; 2]>,
    /// Stable caller-owned keys, one per vertex.
    pub vertex_keys: Vec<u64>,
    /// Oriented triangles. Even permutations preserve orientation.
    pub faces: Vec<[u32; 3]>,
    /// Canonical unoriented edges, sorted by storage vertex index.
    pub edges: Vec<[u32; 2]>,
    vertex_ids: Vec<TriFeatureId>,
    edge_ids: Vec<TriFeatureId>,
    face_ids: Vec<TriFeatureId>,
    edge_index: BTreeMap<[u32; 2], usize>,
    d1_rows: Vec<Vec<(usize, i8)>>,
    edge_measures: Vec<f64>,
    face_measures: Vec<f64>,
}

impl TriComplex2 {
    /// Build an admissible complex under an explicit durable lineage and
    /// caller-owned vertex keys.
    ///
    /// Construction refuses payloads without a two-cell, malformed indices,
    /// repeated keys/cells, non-manifold or incoherently oriented edges,
    /// non-finite embedding data, negative axisymmetric radii, and
    /// non-representable measures.
    pub fn from_triangles(
        lineage: TriComplex2LineageId,
        vertices: Vec<[f64; 2]>,
        vertex_keys: Vec<u64>,
        faces: Vec<[u32; 3]>,
        metric: Metric2,
    ) -> Result<Self, TriComplex2Error> {
        match metric {
            Metric2::Planar { thickness } => {
                validate_positive_metric("thickness", thickness)
                    .map_err(TriComplex2Error::Metric)?;
            }
            Metric2::Axisymmetric { angular_span } => {
                validate_positive_metric("angular span", angular_span)
                    .map_err(TriComplex2Error::Metric)?;
                if angular_span > core::f64::consts::TAU {
                    return Err(TriComplex2Error::Metric(
                        Metric2Error::AngularSpanExceedsTurn {
                            bits: angular_span.to_bits(),
                        },
                    ));
                }
            }
        }
        if vertex_keys.len() != vertices.len() {
            return Err(TriComplex2Error::VertexKeyCount {
                vertices: vertices.len(),
                keys: vertex_keys.len(),
            });
        }
        if faces.is_empty() {
            return Err(TriComplex2Error::NoFaces);
        }

        let mut seen_keys = BTreeMap::new();
        for (vertex, &key) in vertex_keys.iter().enumerate() {
            if let Some(first) = seen_keys.insert(key, vertex) {
                return Err(TriComplex2Error::DuplicateVertexKey {
                    key,
                    first,
                    second: vertex,
                });
            }
        }
        for (vertex, coordinates) in vertices.iter().enumerate() {
            for (axis, value) in coordinates.iter().copied().enumerate() {
                if !value.is_finite() {
                    return Err(TriComplex2Error::NonFiniteCoordinate {
                        vertex,
                        axis,
                        bits: value.to_bits(),
                    });
                }
            }
            if metric.is_axisymmetric() && coordinates[0] < 0.0 {
                return Err(TriComplex2Error::NegativeRadius {
                    vertex,
                    bits: coordinates[0].to_bits(),
                });
            }
        }

        let mut edge_set = BTreeSet::new();
        let mut face_set = BTreeMap::new();
        for (face_index, face) in faces.iter().copied().enumerate() {
            for (local, vertex) in face.iter().copied().enumerate() {
                if vertex as usize >= vertices.len() {
                    return Err(TriComplex2Error::VertexIndexOutOfRange {
                        face: face_index,
                        local,
                        vertex,
                        vertex_count: vertices.len(),
                    });
                }
            }
            if face[0] == face[1] || face[1] == face[2] || face[2] == face[0] {
                return Err(TriComplex2Error::RepeatedFaceVertex {
                    face: face_index,
                    vertices: face,
                });
            }
            let mut canonical_face = face;
            canonical_face.sort_unstable();
            if let Some(first) = face_set.insert(canonical_face, face_index) {
                return Err(TriComplex2Error::DuplicateFace {
                    first,
                    second: face_index,
                    vertices: canonical_face,
                });
            }
            for (from, to) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                edge_set.insert(canonical_edge(from, to));
            }
        }

        let edges: Vec<[u32; 2]> = edge_set.into_iter().collect();
        let edge_index: BTreeMap<[u32; 2], usize> = edges
            .iter()
            .enumerate()
            .map(|(i, &edge)| (edge, i))
            .collect();
        let mut d1_rows = Vec::with_capacity(faces.len());
        let mut edge_uses = vec![Vec::new(); edges.len()];
        for (face_index, face) in faces.iter().copied().enumerate() {
            let mut row = Vec::with_capacity(3);
            for (from, to) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                let edge = canonical_edge(from, to);
                let sign = if [from, to] == edge { 1 } else { -1 };
                let global_edge = edge_index[&edge];
                row.push((global_edge, sign));
                edge_uses[global_edge].push((face_index, sign));
            }
            d1_rows.push(row);
        }
        for (edge, uses) in edges.iter().copied().zip(&edge_uses) {
            if uses.len() > 2 {
                return Err(TriComplex2Error::NonManifoldEdge {
                    edge,
                    incident_faces: uses.len(),
                });
            }
            if let [(first_face, first_sign), (second_face, second_sign)] = uses.as_slice()
                && first_sign == second_sign
            {
                return Err(TriComplex2Error::IncoherentOrientation {
                    edge,
                    first_face: *first_face,
                    second_face: *second_face,
                });
            }
        }

        let vertex_ids = vertex_keys
            .iter()
            .map(|key| tri_feature_id(lineage, 0, core::slice::from_ref(key)))
            .collect::<Result<Vec<_>, _>>()?;
        let edge_ids = edges
            .iter()
            .map(|&[a, b]| {
                tri_feature_id(
                    lineage,
                    1,
                    &[vertex_keys[a as usize], vertex_keys[b as usize]],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let face_ids = faces
            .iter()
            .map(|&[a, b, c]| {
                tri_feature_id(
                    lineage,
                    2,
                    &[
                        vertex_keys[a as usize],
                        vertex_keys[b as usize],
                        vertex_keys[c as usize],
                    ],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;

        let face_measures = faces
            .iter()
            .copied()
            .enumerate()
            .map(|(face, indices)| face_measure(metric, &vertices, face, indices))
            .collect::<Result<Vec<_>, _>>()?;
        let edge_measures = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(edge, indices)| edge_measure(metric, &vertices, edge, indices))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            lineage,
            metric,
            vertices,
            vertex_keys,
            faces,
            edges,
            vertex_ids,
            edge_ids,
            face_ids,
            edge_index,
            d1_rows,
            edge_measures,
            face_measures,
        })
    }

    /// Build with stable index-derived vertex keys under an exact namespace.
    /// Appending refinement vertices preserves all prior keys; reindexing does
    /// not, so durable machine graphs should call [`Self::from_triangles`]
    /// with explicit keys.
    pub fn from_indexed_triangles(
        namespace: &str,
        vertices: Vec<[f64; 2]>,
        faces: Vec<[u32; 3]>,
        metric: Metric2,
    ) -> Result<Self, TriComplex2Error> {
        let vertex_keys = (0..vertices.len())
            .map(|index| u64::try_from(index).map_err(|_| TriComplex2Error::VertexCountOverflow))
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_triangles(
            tri_complex2_lineage_id(namespace)?,
            vertices,
            vertex_keys,
            faces,
            metric,
        )
    }

    /// Topological dimension of this cell complex.
    #[must_use]
    pub const fn topological_dimension(&self) -> u8 {
        2
    }

    /// Coordinate embedding dimension. Axisymmetric weighting does not turn
    /// the stored meridian mesh into a 3-D cell complex.
    #[must_use]
    pub const fn embedding_dimension(&self) -> u8 {
        2
    }

    /// Durable lineage shared across related revisions.
    #[must_use]
    pub const fn lineage(&self) -> TriComplex2LineageId {
        self.lineage
    }

    /// Declared integration measure.
    #[must_use]
    pub const fn metric(&self) -> Metric2 {
        self.metric
    }

    /// Stable typed vertex identities in storage order.
    #[must_use]
    pub fn vertex_ids(&self) -> &[TriFeatureId] {
        &self.vertex_ids
    }

    /// Stable typed edge identities in canonical edge order.
    #[must_use]
    pub fn edge_ids(&self) -> &[TriFeatureId] {
        &self.edge_ids
    }

    /// Stable typed face identities in face storage order.
    #[must_use]
    pub fn face_ids(&self) -> &[TriFeatureId] {
        &self.face_ids
    }

    /// Exact signed incidence d0: vertices → canonical edges.
    #[must_use]
    pub fn d0(&self) -> Incidence {
        Incidence {
            rows: self
                .edges
                .iter()
                .map(|&[a, b]| vec![(a as usize, -1), (b as usize, 1)])
                .collect(),
            cols: self.vertices.len(),
        }
    }

    /// Exact signed incidence d1: canonical edges → oriented faces.
    #[must_use]
    pub fn d1(&self) -> Incidence {
        Incidence {
            rows: self.d1_rows.clone(),
            cols: self.edges.len(),
        }
    }

    /// Prevalidated edge measure under the declared metric.
    #[must_use]
    pub fn edge_measure(&self, edge: usize) -> Option<f64> {
        self.edge_measures.get(edge).copied()
    }

    /// Prevalidated face measure under the declared metric.
    #[must_use]
    pub fn face_measure(&self, face: usize) -> Option<f64> {
        self.face_measures.get(face).copied()
    }

    /// Extract the oriented outer boundary trace of the complete complex.
    pub fn boundary_trace(&self) -> Result<TraceMap2, TriComplex2Error> {
        self.trace_for_faces(0..self.faces.len())
    }

    /// Extract the oriented boundary of a selected face subcomplex.
    ///
    /// Interior selected-selected edges cancel exactly; selected-unselected
    /// edges become deterministic interface traces. Duplicate or out-of-range
    /// face selections refuse.
    pub fn trace_for_faces(
        &self,
        face_indices: impl IntoIterator<Item = usize>,
    ) -> Result<TraceMap2, TriComplex2Error> {
        let mut selected = BTreeSet::new();
        for face in face_indices {
            if face >= self.faces.len() {
                return Err(TriComplex2Error::TraceFaceOutOfRange {
                    face,
                    face_count: self.faces.len(),
                });
            }
            if !selected.insert(face) {
                return Err(TriComplex2Error::DuplicateTraceFace { face });
            }
        }

        let mut sums = vec![0i8; self.edges.len()];
        let mut source_faces = vec![None; self.edges.len()];
        for face in selected {
            for &(edge, sign) in &self.d1_rows[face] {
                sums[edge] += sign;
                source_faces[edge] = Some(face);
            }
        }
        let mut global_vertices = BTreeSet::new();
        let mut trace_edges = Vec::new();
        for (global_edge, sign) in sums.into_iter().enumerate() {
            if sign == 0 {
                continue;
            }
            debug_assert!(sign == -1 || sign == 1);
            let [a, b] = self.edges[global_edge];
            let oriented_vertices = if sign == 1 { [a, b] } else { [b, a] };
            global_vertices.extend(oriented_vertices);
            trace_edges.push(TraceEdge2 {
                global_edge,
                source_face: source_faces[global_edge]
                    .expect("nonzero trace coefficient has source"),
                oriented_vertices,
                feature_id: self.edge_ids[global_edge],
            });
        }
        let vertices: Vec<u32> = global_vertices.into_iter().collect();
        let trace_vertex: BTreeMap<u32, usize> = vertices
            .iter()
            .copied()
            .enumerate()
            .map(|(trace, global)| (global, trace))
            .collect();
        let d0 = Incidence {
            rows: trace_edges
                .iter()
                .map(|edge| {
                    vec![
                        (trace_vertex[&edge.oriented_vertices[0]], -1),
                        (trace_vertex[&edge.oriented_vertices[1]], 1),
                    ]
                })
                .collect(),
            cols: vertices.len(),
        };
        Ok(TraceMap2 {
            vertices,
            edges: trace_edges,
            d0,
        })
    }

    /// Locate a canonical edge without exposing the internal map.
    #[must_use]
    pub fn edge_index(&self, a: u32, b: u32) -> Option<usize> {
        self.edge_index.get(&canonical_edge(a, b)).copied()
    }
}

fn canonical_edge(a: u32, b: u32) -> [u32; 2] {
    if a < b { [a, b] } else { [b, a] }
}

fn signed_area(vertices: &[[f64; 2]], [a, b, c]: [u32; 3]) -> f64 {
    let [ax, ay] = vertices[a as usize];
    let [bx, by] = vertices[b as usize];
    let [cx, cy] = vertices[c as usize];
    0.5 * ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
}

fn face_measure(
    metric: Metric2,
    vertices: &[[f64; 2]],
    face: usize,
    indices: [u32; 3],
) -> Result<f64, TriComplex2Error> {
    let area = signed_area(vertices, indices).abs();
    if area == 0.0 {
        return Err(TriComplex2Error::DegenerateFace { face });
    }
    let measure = match metric {
        Metric2::Planar { thickness } => area * thickness,
        Metric2::Axisymmetric { angular_span } => {
            let [a, b, c] = indices;
            let mean_radius = vertices[a as usize][0] / 3.0
                + vertices[b as usize][0] / 3.0
                + vertices[c as usize][0] / 3.0;
            area * angular_span * mean_radius
        }
    };
    if !area.is_finite() || !measure.is_finite() || measure <= 0.0 {
        return Err(TriComplex2Error::NonRepresentableMeasure {
            topological_dimension: 2,
            cell: face,
        });
    }
    Ok(measure)
}

fn edge_measure(
    metric: Metric2,
    vertices: &[[f64; 2]],
    edge: usize,
    [a, b]: [u32; 2],
) -> Result<f64, TriComplex2Error> {
    let [ax, ay] = vertices[a as usize];
    let [bx, by] = vertices[b as usize];
    let length = (bx - ax).hypot(by - ay);
    let measure = match metric {
        Metric2::Planar { thickness } => length * thickness,
        Metric2::Axisymmetric { angular_span } => length * angular_span * (ax / 2.0 + bx / 2.0),
    };
    if !length.is_finite() || length <= 0.0 || !measure.is_finite() || measure < 0.0 {
        return Err(TriComplex2Error::NonRepresentableMeasure {
            topological_dimension: 1,
            cell: edge,
        });
    }
    Ok(measure)
}

/// Fail-closed construction or trace-extraction error for [`TriComplex2`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriComplex2Error {
    /// A durable lineage namespace was empty.
    EmptyLineageNamespace,
    /// Canonical typed identity construction refused.
    Identity(CanonicalError),
    /// Metric metadata was invalid.
    Metric(Metric2Error),
    /// No two-dimensional cell was supplied.
    NoFaces,
    /// Vertex and key tables differed in length.
    VertexKeyCount {
        /// Coordinate count.
        vertices: usize,
        /// Key count.
        keys: usize,
    },
    /// Two vertices reused one stable key.
    DuplicateVertexKey {
        /// Reused key.
        key: u64,
        /// First storage index.
        first: usize,
        /// Second storage index.
        second: usize,
    },
    /// A face referenced a missing vertex.
    VertexIndexOutOfRange {
        /// Face index.
        face: usize,
        /// Local corner index.
        local: usize,
        /// Refused vertex index.
        vertex: u32,
        /// Available vertex count.
        vertex_count: usize,
    },
    /// A face repeated one vertex.
    RepeatedFaceVertex {
        /// Face index.
        face: usize,
        /// Refused face.
        vertices: [u32; 3],
    },
    /// Two face rows described the same unoriented cell.
    DuplicateFace {
        /// First face index.
        first: usize,
        /// Second face index.
        second: usize,
        /// Canonical vertex triple.
        vertices: [u32; 3],
    },
    /// More than two faces used an edge.
    NonManifoldEdge {
        /// Canonical edge.
        edge: [u32; 2],
        /// Incident face count.
        incident_faces: usize,
    },
    /// Two adjacent faces traversed their shared edge in the same direction.
    IncoherentOrientation {
        /// Canonical edge.
        edge: [u32; 2],
        /// First incident face.
        first_face: usize,
        /// Second incident face.
        second_face: usize,
    },
    /// An embedding coordinate was NaN or infinite.
    NonFiniteCoordinate {
        /// Vertex index.
        vertex: usize,
        /// Coordinate axis.
        axis: usize,
        /// Exact refused bits.
        bits: u64,
    },
    /// An axisymmetric radius was negative.
    NegativeRadius {
        /// Vertex index.
        vertex: usize,
        /// Exact refused bits.
        bits: u64,
    },
    /// A face had zero embedded area.
    DegenerateFace {
        /// Face index.
        face: usize,
    },
    /// A cell's weighted measure underflowed, overflowed, or was invalid.
    NonRepresentableMeasure {
        /// Cell dimension.
        topological_dimension: u8,
        /// Cell index in its canonical table.
        cell: usize,
    },
    /// A trace selected a missing face.
    TraceFaceOutOfRange {
        /// Refused face index.
        face: usize,
        /// Available face count.
        face_count: usize,
    },
    /// A trace selected one face twice.
    DuplicateTraceFace {
        /// Duplicate face index.
        face: usize,
    },
    /// The platform vertex count could not be represented by stable u64 keys.
    VertexCountOverflow,
}

impl fmt::Display for TriComplex2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyLineageNamespace => f.write_str("TriComplex2 lineage namespace is empty"),
            Self::Identity(error) => write!(f, "TriComplex2 identity refused: {error}"),
            Self::Metric(error) => write!(f, "TriComplex2 metric refused: {error}"),
            Self::NoFaces => f.write_str("TriComplex2 requires at least one face"),
            Self::VertexKeyCount { vertices, keys } => {
                write!(
                    f,
                    "TriComplex2 has {vertices} vertices but {keys} stable keys"
                )
            }
            Self::DuplicateVertexKey { key, first, second } => write!(
                f,
                "TriComplex2 stable key {key} is shared by vertices {first} and {second}"
            ),
            Self::VertexIndexOutOfRange {
                face,
                local,
                vertex,
                vertex_count,
            } => write!(
                f,
                "TriComplex2 face {face} corner {local} references vertex {vertex}, count {vertex_count}"
            ),
            Self::RepeatedFaceVertex { face, vertices } => {
                write!(f, "TriComplex2 face {face} repeats a vertex: {vertices:?}")
            }
            Self::DuplicateFace {
                first,
                second,
                vertices,
            } => write!(
                f,
                "TriComplex2 faces {first} and {second} duplicate cell {vertices:?}"
            ),
            Self::NonManifoldEdge {
                edge,
                incident_faces,
            } => write!(
                f,
                "TriComplex2 edge {edge:?} has {incident_faces} incident faces"
            ),
            Self::IncoherentOrientation {
                edge,
                first_face,
                second_face,
            } => write!(
                f,
                "TriComplex2 faces {first_face} and {second_face} traverse edge {edge:?} alike"
            ),
            Self::NonFiniteCoordinate { vertex, axis, bits } => write!(
                f,
                "TriComplex2 vertex {vertex} axis {axis} is non-finite (bits 0x{bits:016x})"
            ),
            Self::NegativeRadius { vertex, bits } => write!(
                f,
                "TriComplex2 vertex {vertex} has negative radius (bits 0x{bits:016x})"
            ),
            Self::DegenerateFace { face } => {
                write!(f, "TriComplex2 face {face} has zero embedded area")
            }
            Self::NonRepresentableMeasure {
                topological_dimension,
                cell,
            } => write!(
                f,
                "TriComplex2 dimension-{topological_dimension} cell {cell} has non-representable measure"
            ),
            Self::TraceFaceOutOfRange { face, face_count } => write!(
                f,
                "TriComplex2 trace selects face {face}, count {face_count}"
            ),
            Self::DuplicateTraceFace { face } => {
                write!(f, "TriComplex2 trace selects face {face} twice")
            }
            Self::VertexCountOverflow => {
                f.write_str("TriComplex2 vertex count does not fit stable u64 keys")
            }
        }
    }
}

impl core::error::Error for TriComplex2Error {}

/// Oriented hex element storage (SoA-friendly quadruple-pair layout).
/// Incidence operators join with fs-feec's tensor-product families —
/// storage lands now so meshing beads have a target (CONTRACT no-claim).
#[derive(Debug, Clone, Default)]
pub struct HexComplex {
    /// Vertex count.
    pub vertex_count: usize,
    /// Hexes as vertex octuples (VTK ordering).
    pub hexes: Vec<[u32; 8]>,
}
