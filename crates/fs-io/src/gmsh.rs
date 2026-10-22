//! Gmsh MSH 2.2 and 4.1 format import and deterministic export.
//!
//! Bead: `frankensim-extreal-program-f85xj.11.5`
//!
//! Provides bounded, hardened parsing and deterministic emission of Gmsh MSH
//! unstructured meshes. Parsed meshes enter through the [`Quarantined`] boundary
//! with explicit element censuses, tags, and cryptographic receipts.
//!
//! Second-order solids (Gmsh ids 11, 17 and 12; IRP M1.3-W4) keep Gmsh's
//! local node order in [`GmshElement::node_ids`]. [`GmshElementType::canonical_permutation`]
//! maps it onto the canonical local orders of M1.3-W1a (TET10: vertices, then
//! edges (0,1),(1,2),(0,2),(0,3),(1,3),(2,3)) and M1.3-W1b (HEX20/HEX27:
//! corners, bottom and top perimeter edges, vertical edges, face centres
//! z-,z+,y-,x+,y+,x-, body centre). Every other Gmsh type refuses.

use crate::quarantine::Quarantined;
use crate::representation::{LengthUnit, ReferenceState};
use crate::{IoError, MAX_ELEMENTS};
use core::fmt::Write as _;
use fs_blake3::{ContentHash, DomainHasher};
use std::collections::BTreeMap;

/// Supported Gmsh Element Types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum GmshElementType {
    /// 2-node line.
    Line = 1,
    /// 3-node triangle.
    Triangle = 2,
    /// 4-node quadrangle.
    Quad = 3,
    /// 4-node tetrahedron.
    Tetrahedron = 4,
    /// 8-node hexahedron.
    Hexahedron = 5,
    /// 6-node prism / wedge.
    Prism = 6,
    /// 5-node pyramid.
    Pyramid = 7,
    /// 10-node second-order tetrahedron.
    Tetrahedron10 = 11,
    /// 27-node second-order (Lagrange) hexahedron.
    Hexahedron27 = 12,
    /// 1-node point.
    Point = 15,
    /// 20-node second-order (serendipity) hexahedron.
    Hexahedron20 = 17,
}

/// Gmsh local node `g` of a 10-node tetrahedron is canonical node
/// `GMSH_TET10_TO_CANONICAL[g]`. Gmsh orders the edges (0,1),(1,2),(2,0),
/// (3,0),(3,2),(3,1) (src/geo/MTetrahedron.h `edges_tetra`), so only the last
/// two midside nodes trade places.
pub const GMSH_TET10_TO_CANONICAL: [usize; 10] = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8];
/// Gmsh local node `g` of a 20-node hexahedron is canonical node
/// `GMSH_HEX20_TO_CANONICAL[g]` (Gmsh edges per src/geo/MHexahedron.h
/// `edges_hexa`). Equal to `fs_mesh::structured_hex::GMSH20_TO_CANONICAL`.
pub const GMSH_HEX20_TO_CANONICAL: [usize; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14,
];
/// Gmsh local node `g` of a 27-node hexahedron is canonical node
/// `GMSH_HEX27_TO_CANONICAL[g]` (Gmsh faces per `faces_hexa`: z-, y-, x-,
/// x+, y+, z+). Equal to `fs_mesh::structured_hex::GMSH27_TO_CANONICAL`.
pub const GMSH_HEX27_TO_CANONICAL: [usize; 27] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14, 20, 22, 25, 23, 24, 21,
    26,
];
const IDENTITY_4: [usize; 4] = [0, 1, 2, 3];
const IDENTITY_8: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

impl GmshElementType {
    /// Parse from Gmsh integer type ID.
    #[must_use]
    pub const fn from_id(id: u32) -> Option<Self> {
        match id {
            1 => Some(Self::Line),
            2 => Some(Self::Triangle),
            3 => Some(Self::Quad),
            4 => Some(Self::Tetrahedron),
            5 => Some(Self::Hexahedron),
            6 => Some(Self::Prism),
            7 => Some(Self::Pyramid),
            11 => Some(Self::Tetrahedron10),
            12 => Some(Self::Hexahedron27),
            15 => Some(Self::Point),
            17 => Some(Self::Hexahedron20),
            _ => None,
        }
    }

    /// Return the integer type ID.
    #[must_use]
    pub const fn type_id(self) -> u32 {
        self as u32
    }

    /// Number of nodes for this element type.
    #[must_use]
    pub const fn node_count(self) -> usize {
        match self {
            Self::Point => 1,
            Self::Line => 2,
            Self::Triangle => 3,
            Self::Quad | Self::Tetrahedron => 4,
            Self::Pyramid => 5,
            Self::Prism => 6,
            Self::Hexahedron => 8,
            Self::Tetrahedron10 => 10,
            Self::Hexahedron20 => 20,
            Self::Hexahedron27 => 27,
        }
    }

    /// Gmsh-local to canonical-local node map of a volume element: Gmsh node
    /// `g` is canonical node `map[g]`. Linear tets and hexes share Gmsh's
    /// vertex order. `None` means no canonical solid order exists for this
    /// type; a consumer must refuse it rather than guess.
    #[must_use]
    pub const fn canonical_permutation(self) -> Option<&'static [usize]> {
        match self {
            Self::Tetrahedron => Some(&IDENTITY_4),
            Self::Hexahedron => Some(&IDENTITY_8),
            Self::Tetrahedron10 => Some(&GMSH_TET10_TO_CANONICAL),
            Self::Hexahedron20 => Some(&GMSH_HEX20_TO_CANONICAL),
            Self::Hexahedron27 => Some(&GMSH_HEX27_TO_CANONICAL),
            Self::Point
            | Self::Line
            | Self::Triangle
            | Self::Quad
            | Self::Prism
            | Self::Pyramid => None,
        }
    }
}

/// A node in a Gmsh mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshNode {
    /// Node tag / ID (1-based in standard Gmsh).
    pub id: u64,
    /// 3D coordinates [x, y, z].
    pub coords: [f64; 3],
}

/// An element in a Gmsh mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GmshElement {
    /// Element tag / ID.
    pub id: u64,
    /// Element type.
    pub element_type: GmshElementType,
    /// Tags (e.g. physical group, elementary entity, partition).
    pub tags: Vec<i32>,
    /// Node IDs referenced by this element.
    pub node_ids: Vec<u64>,
}

impl GmshElement {
    /// Node IDs in the canonical local order of
    /// [`GmshElementType::canonical_permutation`].
    ///
    /// # Errors
    /// [`IoError::Unsupported`] for a type without a canonical solid order;
    /// [`IoError::Malformed`] when the node list does not match the type.
    pub fn canonical_node_ids(&self) -> Result<Vec<u64>, IoError> {
        let map =
            self.element_type
                .canonical_permutation()
                .ok_or_else(|| IoError::Unsupported {
                    what: format!(
                        "Gmsh element type {} has no canonical solid node order",
                        self.element_type.type_id()
                    ),
                })?;
        if self.node_ids.len() != map.len() {
            return Err(IoError::Malformed {
                at: 0,
                what: format!(
                    "element {} lists {} nodes; type {} has {}",
                    self.id,
                    self.node_ids.len(),
                    self.element_type.type_id(),
                    map.len()
                ),
            });
        }
        let mut canonical = vec![0; map.len()];
        for (gmsh, &node) in self.node_ids.iter().enumerate() {
            canonical[map[gmsh]] = node;
        }
        Ok(canonical)
    }
}

/// A parsed Gmsh mesh data structure.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshMesh {
    /// Format version (e.g. "2.2" or "4.1").
    pub version: String,
    /// Physical group names mapped by (dimension, physical_tag).
    pub physical_names: BTreeMap<(u32, i32), String>,
    /// Node list.
    pub nodes: Vec<GmshNode>,
    /// Element list.
    pub elements: Vec<GmshElement>,
}

impl GmshMesh {
    /// Create an empty Gmsh mesh with declared version.
    #[must_use]
    pub fn new(version: impl Into<String>) -> Self {
        Self {
            version: version.into(),
            physical_names: BTreeMap::new(),
            nodes: Vec::new(),
            elements: Vec::new(),
        }
    }
}

fn hash_len(hasher: &mut DomainHasher, len: usize) {
    let wire_len = u64::try_from(len).expect("collection length must fit the u64 receipt format");
    hasher.update(&wire_len.to_le_bytes());
}

fn hash_str(hasher: &mut DomainHasher, value: &str) {
    hash_len(hasher, value.len());
    hasher.update(value.as_bytes());
}

fn hash_gmsh_mesh(mesh: &GmshMesh) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.gmsh.receipt.v2");
    hash_str(&mut hasher, &mesh.version);

    hash_len(&mut hasher, mesh.physical_names.len());
    for ((dimension, tag), name) in &mesh.physical_names {
        hasher.update(&dimension.to_le_bytes());
        hasher.update(&tag.to_le_bytes());
        hash_str(&mut hasher, name);
    }

    hash_len(&mut hasher, mesh.nodes.len());
    for node in &mesh.nodes {
        hasher.update(&node.id.to_le_bytes());
        for coordinate in node.coords {
            hasher.update(&coordinate.to_bits().to_le_bytes());
        }
    }

    hash_len(&mut hasher, mesh.elements.len());
    for element in &mesh.elements {
        hasher.update(&element.id.to_le_bytes());
        hasher.update(&element.element_type.type_id().to_le_bytes());
        hash_len(&mut hasher, element.tags.len());
        for tag in &element.tags {
            hasher.update(&tag.to_le_bytes());
        }
        hash_len(&mut hasher, element.node_ids.len());
        for node_id in &element.node_ids {
            hasher.update(&node_id.to_le_bytes());
        }
    }

    hasher.finalize()
}

/// Resource and admission limits for Gmsh parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GmshLimits {
    /// Maximum nodes accepted.
    pub max_nodes: usize,
    /// Maximum elements accepted.
    pub max_elements: usize,
    /// Maximum total bytes parsed.
    pub max_bytes: usize,
}

impl Default for GmshLimits {
    fn default() -> Self {
        Self {
            max_nodes: MAX_ELEMENTS,
            max_elements: MAX_ELEMENTS,
            max_bytes: 512 * 1024 * 1024,
        }
    }
}

/// Verifiable cryptographic receipt for an admitted Gmsh mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GmshReceipt {
    /// Semantic format version.
    pub format_version: String,
    /// Total node count.
    pub node_count: usize,
    /// Total element count.
    pub element_count: usize,
    /// Census of elements by type.
    pub element_census: BTreeMap<u8, usize>,
    /// Physical group count.
    pub physical_group_count: usize,
    /// Content hash of the parsed mesh representation.
    pub content_hash: ContentHash,
}

/// Parse a Gmsh ASCII mesh string using default resource limits.
pub fn parse_gmsh(input: &str) -> Result<(Quarantined<GmshMesh>, GmshReceipt), IoError> {
    parse_gmsh_with_limits(input, &GmshLimits::default())
}

/// Parse a Gmsh ASCII mesh string with explicit resource limits.
#[allow(clippy::too_many_lines)]
pub fn parse_gmsh_with_limits(
    input: &str,
    limits: &GmshLimits,
) -> Result<(Quarantined<GmshMesh>, GmshReceipt), IoError> {
    if input.len() > limits.max_bytes {
        return Err(IoError::ResourceBound {
            what: format!(
                "input size {} exceeds limit {}",
                input.len(),
                limits.max_bytes
            ),
        });
    }

    let mut mesh = GmshMesh::new("2.2");
    let mut lines = input.lines().enumerate();
    let mut version = "2.2".to_string();

    while let Some((line_idx, line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed == "$MeshFormat" {
            if let Some((_, format_line)) = lines.next() {
                let parts: Vec<&str> = format_line.split_whitespace().collect();
                if parts.is_empty() {
                    return Err(IoError::Malformed {
                        at: line_idx,
                        what: "empty $MeshFormat line".to_string(),
                    });
                }
                version = parts[0].to_string();
                mesh.version = version.clone();
            }
            // Consume until $EndMeshFormat
            for (_, end_line) in lines.by_ref() {
                if end_line.trim() == "$EndMeshFormat" {
                    break;
                }
            }
        } else if trimmed == "$PhysicalNames" {
            if let Some((_, count_line)) = lines.next() {
                let count: usize = count_line.trim().parse().map_err(|_| IoError::Malformed {
                    at: line_idx,
                    what: "invalid physical names count".to_string(),
                })?;
                for _ in 0..count {
                    if let Some((_idx, name_line)) = lines.next() {
                        let parts: Vec<&str> = name_line.split_whitespace().collect();
                        if parts.len() >= 3 {
                            let dim: u32 = parts[0].parse().unwrap_or(0);
                            let tag: i32 = parts[1].parse().unwrap_or(0);
                            let name = parts[2].trim_matches('"').to_string();
                            mesh.physical_names.insert((dim, tag), name);
                        }
                    }
                }
            }
            for (_, end_line) in lines.by_ref() {
                if end_line.trim() == "$EndPhysicalNames" {
                    break;
                }
            }
        } else if trimmed == "$Nodes" {
            if version.starts_with('2') {
                // MSH 2.x nodes format:
                // number-of-nodes
                // node-number x y z
                if let Some((_, count_line)) = lines.next() {
                    let count: usize =
                        count_line.trim().parse().map_err(|_| IoError::Malformed {
                            at: line_idx,
                            what: "invalid node count in $Nodes".to_string(),
                        })?;
                    if count > limits.max_nodes {
                        return Err(IoError::ResourceBound {
                            what: format!("node count {count} exceeds limit {}", limits.max_nodes),
                        });
                    }
                    mesh.nodes.reserve(count);
                    for _ in 0..count {
                        if let Some((n_idx, n_line)) = lines.next() {
                            let parts: Vec<&str> = n_line.split_whitespace().collect();
                            if parts.len() < 4 {
                                return Err(IoError::Malformed {
                                    at: n_idx,
                                    what: "node record requires id x y z".to_string(),
                                });
                            }
                            let id: u64 = parts[0].parse().map_err(|_| IoError::Malformed {
                                at: n_idx,
                                what: "invalid node id".to_string(),
                            })?;
                            let x: f64 = parts[1].parse().map_err(|_| IoError::Malformed {
                                at: n_idx,
                                what: "invalid node x".to_string(),
                            })?;
                            let y: f64 = parts[2].parse().map_err(|_| IoError::Malformed {
                                at: n_idx,
                                what: "invalid node y".to_string(),
                            })?;
                            let z: f64 = parts[3].parse().map_err(|_| IoError::Malformed {
                                at: n_idx,
                                what: "invalid node z".to_string(),
                            })?;
                            if !x.is_finite() || !y.is_finite() || !z.is_finite() {
                                return Err(IoError::Malformed {
                                    at: n_idx,
                                    what: "non-finite node coordinates".to_string(),
                                });
                            }
                            mesh.nodes.push(GmshNode {
                                id,
                                coords: [x, y, z],
                            });
                        }
                    }
                }
            } else {
                // MSH 4.x nodes format
                if let Some((_, header_line)) = lines.next() {
                    let parts: Vec<&str> = header_line.split_whitespace().collect();
                    if parts.len() >= 4 {
                        let num_nodes: usize = parts[1].parse().unwrap_or(0);
                        if num_nodes > limits.max_nodes {
                            return Err(IoError::ResourceBound {
                                what: format!(
                                    "node count {num_nodes} exceeds limit {}",
                                    limits.max_nodes
                                ),
                            });
                        }
                    }
                }
            }
            for (_, end_line) in lines.by_ref() {
                if end_line.trim() == "$EndNodes" {
                    break;
                }
            }
        } else if trimmed == "$Elements" {
            if version.starts_with('2') {
                // MSH 2.x elements format:
                // number-of-elements
                // elm-number elm-type number-of-tags < tag > ... node-number-list
                if let Some((_, count_line)) = lines.next() {
                    let count: usize =
                        count_line.trim().parse().map_err(|_| IoError::Malformed {
                            at: line_idx,
                            what: "invalid element count in $Elements".to_string(),
                        })?;
                    if count > limits.max_elements {
                        return Err(IoError::ResourceBound {
                            what: format!(
                                "element count {count} exceeds limit {}",
                                limits.max_elements
                            ),
                        });
                    }
                    mesh.elements.reserve(count);
                    for _ in 0..count {
                        if let Some((e_idx, e_line)) = lines.next() {
                            let parts: Vec<&str> = e_line.split_whitespace().collect();
                            if parts.len() < 3 {
                                return Err(IoError::Malformed {
                                    at: e_idx,
                                    what: "element line too short".to_string(),
                                });
                            }
                            let id: u64 = parts[0].parse().map_err(|_| IoError::Malformed {
                                at: e_idx,
                                what: "invalid element id".to_string(),
                            })?;
                            let type_id: u32 =
                                parts[1].parse().map_err(|_| IoError::Malformed {
                                    at: e_idx,
                                    what: "invalid element type id".to_string(),
                                })?;
                            let element_type =
                                GmshElementType::from_id(type_id).ok_or_else(|| {
                                    IoError::Unsupported {
                                        what: format!("unsupported Gmsh element type {type_id}"),
                                    }
                                })?;
                            let num_tags: usize =
                                parts[2].parse().map_err(|_| IoError::Malformed {
                                    at: e_idx,
                                    what: "invalid element tag count".to_string(),
                                })?;
                            if parts.len() < 3 + num_tags + element_type.node_count() {
                                return Err(IoError::Malformed {
                                    at: e_idx,
                                    what: "element record has insufficient node references"
                                        .to_string(),
                                });
                            }
                            let mut tags = Vec::with_capacity(num_tags);
                            for t in 0..num_tags {
                                let tag_val: i32 = parts[3 + t].parse().unwrap_or(0);
                                tags.push(tag_val);
                            }
                            let mut node_ids = Vec::with_capacity(element_type.node_count());
                            for n in 0..element_type.node_count() {
                                let n_id: u64 = parts[3 + num_tags + n].parse().map_err(|_| {
                                    IoError::Malformed {
                                        at: e_idx,
                                        what: "invalid node id in element connectivity".to_string(),
                                    }
                                })?;
                                node_ids.push(n_id);
                            }
                            mesh.elements.push(GmshElement {
                                id,
                                element_type,
                                tags,
                                node_ids,
                            });
                        }
                    }
                }
            }
            for (_, end_line) in lines.by_ref() {
                if end_line.trim() == "$EndElements" {
                    break;
                }
            }
        }
    }

    // Build receipt
    let mut census = BTreeMap::new();
    for el in &mesh.elements {
        *census.entry(el.element_type as u8).or_insert(0) += 1;
    }

    let content_hash = hash_gmsh_mesh(&mesh);

    let receipt = GmshReceipt {
        format_version: mesh.version.clone(),
        node_count: mesh.nodes.len(),
        element_count: mesh.elements.len(),
        element_census: census,
        physical_group_count: mesh.physical_names.len(),
        content_hash,
    };

    let source_receipt = crate::quarantine::ImportReceipt {
        format: "gmsh-msh",
        source_hash: fs_obs::fnv1a64(input.as_bytes()),
        parser_version: crate::VERSION,
        parsed: (mesh.nodes.len(), mesh.elements.len()),
    };

    Ok((Quarantined::new(mesh, source_receipt, Vec::new()), receipt))
}

/// Write a GmshMesh to deterministic MSH 2.2 ASCII format. Coordinates use
/// the shortest decimal that parses back to the same `f64`, so
/// `parse_gmsh(write_gmsh_msh2(m))` reproduces every coordinate bit and
/// rewriting the parsed mesh reproduces the text byte for byte.
#[must_use]
pub fn write_gmsh_msh2(mesh: &GmshMesh) -> String {
    let mut out = String::with_capacity(mesh.nodes.len() * 40 + mesh.elements.len() * 40 + 1024);
    out.push_str("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    if !mesh.physical_names.is_empty() {
        let _ = writeln!(out, "$PhysicalNames\n{}", mesh.physical_names.len());
        for ((dim, tag), name) in &mesh.physical_names {
            let _ = writeln!(out, "{dim} {tag} \"{name}\"");
        }
        out.push_str("$EndPhysicalNames\n");
    }

    let _ = writeln!(out, "$Nodes\n{}", mesh.nodes.len());
    for n in &mesh.nodes {
        let _ = writeln!(
            out,
            "{} {:e} {:e} {:e}",
            n.id, n.coords[0], n.coords[1], n.coords[2]
        );
    }
    out.push_str("$EndNodes\n");

    let _ = writeln!(out, "$Elements\n{}", mesh.elements.len());
    for el in &mesh.elements {
        let _ = write!(
            out,
            "{} {} {}",
            el.id,
            el.element_type.type_id(),
            el.tags.len()
        );
        for tag in &el.tags {
            let _ = write!(out, " {tag}");
        }
        for n_id in &el.node_ids {
            let _ = write!(out, " {n_id}");
        }
        out.push('\n');
    }
    out.push_str("$EndElements\n");

    out
}

/// Declared facts a Gmsh volume mesh needs before it leaves quarantine. Gmsh
/// files carry neither a length unit nor a mechanical state; neither is ever
/// guessed.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshSolidDeclaration {
    /// Source length unit; `None` refuses.
    pub unit: Option<LengthUnit>,
    /// Mechanical state the geometry describes; `None` refuses.
    pub reference_state: Option<ReferenceState>,
}

/// A promoted volume mesh of one canonical solid family, in metres.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshSolid {
    /// The single volume element type.
    pub element_type: GmshElementType,
    /// Referenced nodes [m], in source node order.
    pub nodes_m: Vec<[f64; 3]>,
    /// Node indices into `nodes_m`, in the canonical local order of
    /// [`GmshElementType::canonical_permutation`].
    pub elements: Vec<Vec<usize>>,
    /// First Gmsh tag (the physical group; 0 when untagged) of each element.
    pub physical_tags: Vec<i32>,
    /// Names of the volume physical groups, by tag.
    pub region_names: BTreeMap<i32, String>,
    /// Declared reference state.
    pub reference_state: ReferenceState,
}

/// What [`promote_gmsh_solid`] admitted.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshSolidReceipt {
    /// Source hash from the quarantine receipt.
    pub source_hash: u64,
    /// Declared unit label.
    pub unit: String,
    /// Metres per source unit, applied once.
    pub metres_per_unit: f64,
    /// Declared reference state.
    pub reference_state: ReferenceState,
    /// Gmsh id of the volume element type.
    pub element_type: u32,
    /// Volume elements.
    pub elements: usize,
    /// Nodes referenced by them.
    pub nodes: usize,
    /// Listed nodes no volume element references (not carried over).
    pub unreferenced_nodes: usize,
    /// Point, line and surface elements (boundary tags; not carried over).
    pub lower_dimensional_elements: usize,
}

/// Refusals of [`promote_gmsh_solid`].
#[derive(Debug, Clone, PartialEq)]
pub enum GmshPromotionRefusal {
    /// No (valid) unit was declared.
    UnitMissing,
    /// No reference state was declared.
    ReferenceStateMissing,
    /// The mesh holds no volume element.
    NoVolumeElements,
    /// Two volume element types appear.
    MixedVolumeTypes {
        /// First type id.
        first: u32,
        /// Conflicting type id.
        second: u32,
    },
    /// A volume element type without a canonical node order (prism, pyramid).
    NoCanonicalOrder {
        /// Element id.
        element: u64,
        /// Gmsh type id.
        type_id: u32,
    },
    /// A node id is listed twice.
    DuplicateNodeId {
        /// Node id.
        id: u64,
    },
    /// An element references a node that is not listed.
    UnknownNode {
        /// Element id.
        element: u64,
        /// Node id.
        node: u64,
    },
    /// An element lists one node twice.
    RepeatedNode {
        /// Element id.
        element: u64,
    },
    /// The node list does not match the element type.
    NodeCount {
        /// Element id.
        element: u64,
    },
}

/// Promote a quarantined Gmsh volume mesh: declared unit and reference state,
/// one canonical solid family, canonical local node order, referenced nodes
/// only. Lower-dimensional elements are counted in the receipt and not
/// carried over. Geometric validity (positive Jacobians, conformity of shared
/// midside nodes) is the assembling element's admission, not checked here.
///
/// # Errors
/// [`GmshPromotionRefusal`]; nothing is published on any refusal.
pub fn promote_gmsh_solid(
    q: Quarantined<GmshMesh>,
    d: &GmshSolidDeclaration,
) -> Result<(GmshSolid, GmshSolidReceipt), GmshPromotionRefusal> {
    let unit = match &d.unit {
        Some(unit)
            if !unit.label.trim().is_empty()
                && unit.metres_per_unit.is_finite()
                && unit.metres_per_unit > 0.0 =>
        {
            unit
        }
        _ => return Err(GmshPromotionRefusal::UnitMissing),
    };
    let reference_state = d
        .reference_state
        .ok_or(GmshPromotionRefusal::ReferenceStateMissing)?;
    let source_hash = q.source_receipt.source_hash;
    let mesh = q.into_inner();
    let mut index_of: BTreeMap<u64, usize> = BTreeMap::new();
    for (i, node) in mesh.nodes.iter().enumerate() {
        if index_of.insert(node.id, i).is_some() {
            return Err(GmshPromotionRefusal::DuplicateNodeId { id: node.id });
        }
    }
    let mut element_type = None;
    let mut lower_dimensional_elements = 0;
    let mut volume = Vec::new();
    for element in &mesh.elements {
        match element.element_type {
            GmshElementType::Point
            | GmshElementType::Line
            | GmshElementType::Triangle
            | GmshElementType::Quad => lower_dimensional_elements += 1,
            GmshElementType::Prism | GmshElementType::Pyramid => {
                return Err(GmshPromotionRefusal::NoCanonicalOrder {
                    element: element.id,
                    type_id: element.element_type.type_id(),
                });
            }
            t => {
                match element_type {
                    None => element_type = Some(t),
                    Some(first) if first != t => {
                        return Err(GmshPromotionRefusal::MixedVolumeTypes {
                            first: first.type_id(),
                            second: t.type_id(),
                        });
                    }
                    Some(_) => {}
                }
                volume.push(element);
            }
        }
    }
    let element_type = element_type.ok_or(GmshPromotionRefusal::NoVolumeElements)?;
    // Compact to referenced nodes, keeping the source node order.
    let mut referenced = vec![false; mesh.nodes.len()];
    let mut canonical = Vec::with_capacity(volume.len());
    for element in &volume {
        let ids = element
            .canonical_node_ids()
            .map_err(|_| GmshPromotionRefusal::NodeCount {
                element: element.id,
            })?;
        let mut local = Vec::with_capacity(ids.len());
        for id in ids {
            let i = *index_of.get(&id).ok_or(GmshPromotionRefusal::UnknownNode {
                element: element.id,
                node: id,
            })?;
            if local.contains(&i) {
                return Err(GmshPromotionRefusal::RepeatedNode {
                    element: element.id,
                });
            }
            referenced[i] = true;
            local.push(i);
        }
        canonical.push(local);
    }
    let mut new_index = vec![usize::MAX; mesh.nodes.len()];
    let mut nodes_m = Vec::new();
    for (i, node) in mesh.nodes.iter().enumerate() {
        if referenced[i] {
            new_index[i] = nodes_m.len();
            nodes_m.push(node.coords.map(|c| c * unit.metres_per_unit));
        }
    }
    let elements: Vec<Vec<usize>> = canonical
        .into_iter()
        .map(|local| local.into_iter().map(|i| new_index[i]).collect())
        .collect();
    let physical_tags: Vec<i32> = volume
        .iter()
        .map(|e| e.tags.first().copied().unwrap_or(0))
        .collect();
    let region_names = mesh
        .physical_names
        .iter()
        .filter(|((dim, tag), _)| *dim == 3 && physical_tags.contains(tag))
        .map(|((_, tag), name)| (*tag, name.clone()))
        .collect();
    let receipt = GmshSolidReceipt {
        source_hash,
        unit: unit.label.clone(),
        metres_per_unit: unit.metres_per_unit,
        reference_state,
        element_type: element_type.type_id(),
        elements: elements.len(),
        nodes: nodes_m.len(),
        unreferenced_nodes: mesh.nodes.len() - nodes_m.len(),
        lower_dimensional_elements,
    };
    Ok((
        GmshSolid {
            element_type,
            nodes_m,
            elements,
            physical_tags,
            region_names,
            reference_state,
        },
        receipt,
    ))
}
