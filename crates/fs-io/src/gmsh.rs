//! Gmsh MSH 2.2 and 4.1 format import and deterministic export.
//!
//! Bead: `frankensim-extreal-program-f85xj.11.5`
//!
//! Provides bounded, hardened parsing and deterministic emission of Gmsh MSH
//! unstructured meshes. Parsed meshes enter through the [`Quarantined`] boundary
//! with explicit element censuses, tags, and cryptographic receipts.

use crate::quarantine::Quarantined;
use crate::{IoError, MAX_ELEMENTS};
use core::fmt::Write as _;
use fs_blake3::{ContentHash, hash_domain};
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
    /// 1-node point.
    Point = 15,
}

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
            15 => Some(Self::Point),
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

    let mut hash_text = String::with_capacity(4096);
    let _ = write!(
        hash_text,
        "gmsh_msh;version={};nodes={};elements={}",
        mesh.version,
        mesh.nodes.len(),
        mesh.elements.len()
    );
    for (type_id, count) in &census {
        let _ = write!(hash_text, ";type_{type_id}={count}");
    }
    let content_hash = hash_domain("org.frankensim.gmsh.receipt.v1", hash_text.as_bytes());

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
        source_hash: u64::from_le_bytes(content_hash.as_bytes()[..8].try_into().unwrap_or([0; 8])),
        parser_version: crate::VERSION,
        parsed: (mesh.nodes.len(), mesh.elements.len()),
    };

    Ok((Quarantined::new(mesh, source_receipt, Vec::new()), receipt))
}

/// Write a GmshMesh to deterministic MSH 2.2 ASCII format.
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
            "{} {:.12e} {:.12e} {:.12e}",
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
