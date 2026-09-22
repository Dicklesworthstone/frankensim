//! OBJ geometry with an optional, loss-aware physical-region import.
//!
//! `read_obj` retains its geometry-only contract. `read_obj_document` also
//! retains object names, group membership, material assignments and MTL
//! references, without interpreting visual MTL parameters as physical laws.
//! Neither reader opens referenced files. Positions and faces are supported;
//! texture/normal indices are accepted but not retained. Simple planar polygons
//! use the shared concavity-aware triangulator; crossed/nonplanar faces refuse.

use crate::{IoError, MAX_ELEMENTS};
use fs_geom::Point3;
use fs_rep_mesh::{MAX_POLYGON_VERTICES, PolygonError, Soup, triangulate_polygon};
use std::fmt::Write as _;
use std::ops::Range;

const MAX_LABEL_BYTES: usize = 4 * 1024 * 1024;
const MAX_REGIONS: usize = 100_000;
const MAX_GROUPS: usize = 64;

/// Authored OBJ labels on a contiguous range of triangulated faces.
/// A face can belong to several groups. Object and group names are NOT
/// material identities, and MTL names are NOT constitutive material cards.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ObjRegion {
    pub triangles: Range<usize>,
    pub object: Option<String>,
    pub groups: Vec<String>,
    pub material: Option<String>,
}
impl ObjRegion {
    /// Exact, case-sensitive selection; no fuzzy assignment of physical parts.
    #[must_use]
    pub fn has_label(&self, label: &str) -> bool {
        self.object.as_deref() == Some(label) || self.groups.iter().any(|g| g == label)
    }
}

/// Geometry and its authored region partition. Every output triangle belongs
/// to exactly one region, including unlabelled triangles. Source vertices are
/// not welded, scaled, reoriented or silently repaired.
#[derive(Debug)]
pub struct ObjDocument {
    pub soup: Soup,
    pub regions: Vec<ObjRegion>,
    /// Unresolved `mtllib` references, in source order. Import never follows
    /// these paths (in particular, no path traversal or network I/O occurs).
    pub material_libraries: Vec<String>,
}

/// Import the geometry-only OBJ subset (backwards-compatible, documented loss).
///
/// Triangle faces retain their original indices, including degeneracies for
/// downstream repair/quarantine. Polygon faces must be simple and planar and
/// are bounded by [`MAX_POLYGON_VERTICES`]. Negative indices are supported.
/// Texture coordinates, normals, labels and material references are discarded.
///
/// # Errors
/// [`IoError`] for malformed indices/coordinates, inadmissible polygons, or
/// resource bounds. No partially imported soup is returned.
pub fn read_obj(text: &str) -> Result<Soup, IoError> {
    Ok(read_document(text, false)?.soup)
}

/// Import geometry while retaining `o`, `g`, `usemtl` and `mtllib` semantics.
///
/// Triangulating an n-gon does not lose its labels: every generated triangle
/// inherits the active assignment. Adjacent identical assignments coalesce.
/// Empty `o`, `g` and `usemtl` statements clear the corresponding assignment.
/// Coordinates remain in SOURCE units; OBJ has no dependable SI-unit field.
/// This is an asset import, not permission to treat a render mesh as a plate.
///
/// # Errors
/// In addition to geometry errors, refuses oversized metadata (4 MiB of
/// declared/retained label text, 100,000 regions, 64 groups per assignment).
/// The legacy reader does not allocate or validate ignored metadata.
pub fn read_obj_document(text: &str) -> Result<ObjDocument, IoError> {
    read_document(text, true)
}

fn budget_label(bytes: &mut usize, extra: usize) -> Result<(), IoError> {
    *bytes = bytes.checked_add(extra).filter(|&n| n <= MAX_LABEL_BYTES)
        .ok_or_else(|| IoError::ResourceBound { what: "OBJ label text exceeds the metadata cap".into() })?;
    Ok(())
}

fn read_document(text: &str, labels: bool) -> Result<ObjDocument, IoError> {
    let mut positions: Vec<Point3> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut regions: Vec<ObjRegion> = Vec::new();
    let mut material_libraries = Vec::new();
    let mut object: Option<String> = None;
    let mut groups: Vec<String> = Vec::new();
    let mut material: Option<String> = None;
    let mut label_bytes = 0usize;
    for (ln, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() { continue; }
        let mut it = line.split_whitespace();
        let tag = it.next().unwrap_or("");
        match tag {
            "v" => {
                let mut v = [0.0f64; 3];
                for slot in &mut v {
                    let tok = it.next().ok_or(IoError::Malformed {
                        at: ln + 1, what: "v needs three coordinates".to_string(),
                    })?;
                    *slot = tok.parse::<f64>().map_err(|_| IoError::Malformed {
                        at: ln + 1, what: format!("bad coordinate {tok:?}"),
                    })?;
                    if !slot.is_finite() {
                        return Err(IoError::Malformed { at: ln + 1, what: "non-finite coordinate".into() });
                    }
                }
                if positions.len() == MAX_ELEMENTS {
                    return Err(IoError::ResourceBound { what: "vertex count exceeds the element cap".into() });
                }
                positions.push(Point3::new(v[0], v[1], v[2]));
            }
            "f" => {
                let mut idx: Vec<u32> = Vec::new();
                for tok in it {
                    if idx.len() == MAX_POLYGON_VERTICES {
                        return Err(IoError::ResourceBound {
                            what: format!("OBJ face at line {} exceeds the polygon vertex cap", ln + 1),
                        });
                    }
                    let first = tok.split('/').next().unwrap_or("");
                    let signed: i64 = first.parse().map_err(|_| IoError::Malformed {
                        at: ln + 1, what: format!("bad face index {tok:?}"),
                    })?;
                    let resolved = if signed < 0 {
                        i64::try_from(positions.len()).expect("cap").checked_add(signed)
                    } else { signed.checked_sub(1) };
                    let resolved = resolved.filter(|&n| n >= 0 && n < positions.len() as i64)
                        .ok_or_else(|| IoError::Malformed {
                            at: ln + 1, what: format!("face index {signed} out of range"),
                        })?;
                    idx.push(u32::try_from(resolved).expect("range checked"));
                }
                let first = triangles.len();
                append_face(&positions, &idx, &mut triangles, ln + 1)?;
                if labels {
                    let same = regions.last().is_some_and(|r| {
                        r.object == object && r.groups == groups && r.material == material
                    });
                    if same {
                        if let Some(last) = regions.last_mut() { last.triangles.end = triangles.len(); }
                    } else {
                        if regions.len() == MAX_REGIONS {
                            return Err(IoError::ResourceBound { what: "OBJ region count exceeds the metadata cap".into() });
                        }
                        let extra = object.as_ref().map_or(0, String::len)
                            + material.as_ref().map_or(0, String::len)
                            + groups.iter().map(String::len).sum::<usize>();
                        budget_label(&mut label_bytes, extra)?;
                        regions.push(ObjRegion {
                            triangles: first..triangles.len(), object: object.clone(),
                            groups: groups.clone(), material: material.clone(),
                        });
                    }
                }
            }
            "o" | "usemtl" if labels => {
                let name = line[tag.len()..].trim();
                // Unused declarations also consume the metadata budget.
                budget_label(&mut label_bytes, name.len())?;
                let value = (!name.is_empty()).then(|| name.to_owned());
                if tag == "o" { object = value; } else { material = value; }
            }
            "g" if labels => {
                groups.clear();
                for group in it {
                    if groups.len() == MAX_GROUPS {
                        return Err(IoError::ResourceBound { what: "OBJ group membership exceeds the metadata cap".into() });
                    }
                    budget_label(&mut label_bytes, group.len())?;
                    groups.push(group.to_owned());
                }
            }
            "mtllib" if labels => {
                for library in it {
                    if material_libraries.len() == MAX_REGIONS {
                        return Err(IoError::ResourceBound { what: "OBJ material library count exceeds the metadata cap".into() });
                    }
                    budget_label(&mut label_bytes, library.len())?;
                    material_libraries.push(library.to_owned());
                }
            }
            _ => {} // vt/vn/s/... ignored: see the loss contract above.
        }
    }
    if triangles.is_empty() {
        return Err(IoError::Malformed { at: 0, what: "OBJ contains no faces".to_string() });
    }
    Ok(ObjDocument { soup: Soup { positions, triangles }, regions, material_libraries })
}

// Shared by the text mesh readers: never fan-fill a concave polygon, and
// reserve the entire output face before publishing any of its triangles.
pub(crate) fn append_face(
    positions: &[Point3], indices: &[u32], triangles: &mut Vec<[u32; 3]>, at: usize,
) -> Result<(), IoError> {
    let count = indices.len().checked_sub(2).filter(|&n| n > 0)
        .ok_or_else(|| IoError::Malformed { at, what: "face needs at least three vertices".into() })?;
    if triangles.len().checked_add(count).is_none_or(|n| n > MAX_ELEMENTS) {
        return Err(IoError::ResourceBound { what: "triangle count exceeds the element cap".into() });
    }
    if indices.len() > MAX_POLYGON_VERTICES {
        return Err(IoError::ResourceBound { what: "polygon vertex count exceeds the face cap".into() });
    }
    let face = if indices.len() == 3 { vec![[indices[0], indices[1], indices[2]]] } else {
        triangulate_polygon(positions, indices).map_err(|error| match error {
            PolygonError::Resource(what) => IoError::ResourceBound { what: what.into() },
            error => IoError::Malformed { at, what: error.to_string() },
        })?
    };
    triangles.try_reserve(face.len()).map_err(|_| IoError::ResourceBound {
        what: "face triangle allocation failed".into(),
    })?;
    triangles.extend(face);
    Ok(())
}

/// Export geometry only (deterministic; f64 round-trip precision).
#[must_use]
pub fn write_obj(soup: &Soup) -> String {
    let mut out = String::with_capacity(soup.positions.len() * 32);
    for p in &soup.positions { let _ = writeln!(out, "v {} {} {}", p.x, p.y, p.z); }
    for t in &soup.triangles { let _ = writeln!(out, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1); }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    const QUAD: &str = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n";
    #[test]
    fn relative_concave_faces_preserve_area_and_round_trip() {
        let source = "v 0 0 0\nv 3 0 0\nv 3 3 0\nv 2 3 0\nv 2 1 0\nv 1 1 0\nv 1 3 0\nv 0 3 0\nf -8/1/1 -7/2/1 -6/3/1 -5/4/1 -4/5/1 -3/6/1 -2/7/1 -1/8/1\n";
        let soup = read_obj(source).unwrap();
        assert_eq!(soup.triangles.len(), 6);
        let area: f64 = soup.triangles.iter().map(|t| {
            let [a, b, c] = t.map(|i| soup.positions[i as usize]);
            let area = ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)) * 0.5;
            assert!(area > 0.0); area
        }).sum();
        assert_eq!(area, 7.0);
        let round_trip = read_obj(&write_obj(&soup)).unwrap();
        assert_eq!(round_trip.positions, soup.positions);
        assert_eq!(round_trip.triangles, soup.triangles);
    }
    #[test]
    fn bad_polygon_refuses_at_source_line_without_changing_triangle_repair() {
        let error = read_obj("v 0 0 0\nv 1 1 0\nv 0 1 0\nv 1 0 0\nf 1 2 3 4\n").unwrap_err();
        assert!(matches!(error, IoError::Malformed { at: 5, .. }));
        let soup = read_obj("v 0 0 0\nv 1 0 0\nf 1 1 2\n").unwrap();
        assert_eq!(soup.triangles, vec![[0, 0, 1]]);
    }
    #[test]
    fn polygon_triangulation_preserves_groups_and_materials() {
        let text = format!("{QUAD}mtllib woods.mtl hardware.mtl\no Steinway soundboard\ng panel acoustic\nusemtl Sitka spruce\nf -4 -3 -2 -1\nf 1 2 3\ng lid\nusemtl lacquer\nf 1 3 4\n");
        let doc = read_obj_document(&text).unwrap();
        assert_eq!(doc.soup.triangles.len(), 4);
        assert_eq!(doc.regions.len(), 2);
        assert_eq!(doc.regions[0].triangles, 0..3);
        assert!(doc.regions[0].has_label("panel"));
        assert!(doc.regions[0].has_label("Steinway soundboard"));
        assert!(!doc.regions[0].has_label("Panel"));
        assert_eq!(doc.regions[0].material.as_deref(), Some("Sitka spruce"));
        assert_eq!(doc.regions[1].triangles, 3..4);
        assert_eq!(doc.material_libraries, vec!["woods.mtl", "hardware.mtl"]);
        let old = read_obj(&text).unwrap();
        assert_eq!(old.positions, doc.soup.positions);
        assert_eq!(old.triangles, doc.soup.triangles);
    }
    #[test]
    fn empty_assignments_clear_and_every_triangle_is_covered() {
        let text = format!("{QUAD}f 1 2 3\no board\ng top structural\nusemtl spruce\nf 1 3 4\no\ng\nusemtl\nf 1 2 3\nf 1 3 4\n");
        let doc = read_obj_document(&text).unwrap();
        assert_eq!(doc.regions.len(), 3);
        assert_eq!(doc.regions[2].triangles, 2..4);
        assert_eq!(doc.regions[2].object, None);
        assert_eq!(doc.regions[2].material, None);
        assert!(doc.regions[2].groups.is_empty());
        let covered: Vec<_> = doc.regions.iter().flat_map(|r| r.triangles.clone()).collect();
        assert_eq!(covered, (0..4).collect::<Vec<_>>());
    }
    #[test]
    fn hostile_metadata_refuses_without_changing_legacy_geometry() {
        let groups = std::iter::repeat_n("part", MAX_GROUPS + 1).collect::<Vec<_>>().join(" ");
        let text = format!("{QUAD}g {groups}\nf 1 2 3\n");
        assert!(matches!(read_obj_document(&text), Err(IoError::ResourceBound { .. })));
        assert_eq!(read_obj(&text).unwrap().triangles.len(), 1);
        let mut bytes = MAX_LABEL_BYTES;
        assert!(budget_label(&mut bytes, 1).is_err());
        for index in ["0", "-9223372036854775808", "9223372036854775807"] {
            assert!(read_obj_document(&format!("{QUAD}f {index} 2 3\n")).is_err());
        }
    }
    #[test]
    fn library_paths_are_data_not_io_requests() {
        let doc = read_obj_document(&format!("{QUAD}mtllib ../../private.mtl https://example.invalid/m.mtl\nf 1 2 3\n")).unwrap();
        assert_eq!(doc.material_libraries.len(), 2);
        assert_eq!(doc.material_libraries[0], "../../private.mtl");
    }
}
