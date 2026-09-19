//! OBJ import/export (subset): `v` positions and `f` faces. Simple planar
//! polygons are triangulated without filling concave notches; crossed and
//! nonplanar polygon faces refuse. `v/vt/vn` index forms are accepted, texture
//! and normal indices ignored (documented lossy). Negative (relative) indices
//! are supported. Export retains full f64 round-trip precision.

use crate::{IoError, MAX_ELEMENTS};
use fs_geom::Point3;
use fs_rep_mesh::{MAX_POLYGON_VERTICES, PolygonError, Soup, triangulate_polygon};
use std::fmt::Write as _;

/// Import an OBJ subset.
///
/// Triangle faces retain their original indices, including degeneracies for
/// the existing downstream repair/quarantine policy. Faces with more than
/// three vertices must be simple and planar and are bounded by
/// [`MAX_POLYGON_VERTICES`]; their winding and boundary vertices are preserved.
///
/// # Errors
/// [`IoError`] on malformed indices/coordinates, inadmissible polygons, or
/// resource bounds. No partially imported soup is returned.
pub fn read_obj(text: &str) -> Result<Soup, IoError> {
    let mut positions: Vec<Point3> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    for (ln, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        match it.next() {
            Some("v") => {
                let mut v = [0.0f64; 3];
                for slot in &mut v {
                    let tok = it.next().ok_or(IoError::Malformed {
                        at: ln + 1,
                        what: "v needs three coordinates".to_string(),
                    })?;
                    *slot = tok.parse::<f64>().map_err(|_| IoError::Malformed {
                        at: ln + 1,
                        what: format!("bad coordinate {tok:?}"),
                    })?;
                    if !slot.is_finite() {
                        return Err(IoError::Malformed {
                            at: ln + 1,
                            what: "non-finite coordinate".to_string(),
                        });
                    }
                }
                if positions.len() == MAX_ELEMENTS {
                    return Err(IoError::ResourceBound {
                        what: "vertex count exceeds the element cap".to_string(),
                    });
                }
                positions.push(Point3::new(v[0], v[1], v[2]));
            }
            Some("f") => {
                let mut idx: Vec<u32> = Vec::new();
                for tok in it {
                    if idx.len() == MAX_POLYGON_VERTICES {
                        return Err(IoError::ResourceBound {
                            what: format!("OBJ face at line {} exceeds the polygon vertex cap", ln + 1),
                        });
                    }
                    let first = tok.split('/').next().unwrap_or("");
                    let signed: i64 = first.parse().map_err(|_| IoError::Malformed {
                        at: ln + 1,
                        what: format!("bad face index {tok:?}"),
                    })?;
                    let resolved: i64 = if signed < 0 {
                        i64::try_from(positions.len()).expect("cap") + signed
                    } else {
                        signed - 1
                    };
                    if resolved < 0 || resolved >= i64::try_from(positions.len()).expect("cap") {
                        return Err(IoError::Malformed {
                            at: ln + 1,
                            what: format!("face index {signed} out of range"),
                        });
                    }
                    idx.push(u32::try_from(resolved).expect("range checked"));
                }
                append_face(&positions, &idx, &mut triangles, ln + 1)?;
            }
            _ => {} // vt/vn/usemtl/o/g/s… ignored (documented subset)
        }
    }
    if triangles.is_empty() {
        return Err(IoError::Malformed {
            at: 0,
            what: "OBJ contains no faces".to_string(),
        });
    }
    Ok(Soup { positions, triangles })
}

// Shared by the text mesh readers: never fan-fill a concave polygon, and
// reserve the entire output face before publishing any of its triangles.
pub(crate) fn append_face(
    positions: &[Point3],
    indices: &[u32],
    triangles: &mut Vec<[u32; 3]>,
    at: usize,
) -> Result<(), IoError> {
    let count = indices.len().checked_sub(2).filter(|&n| n > 0)
        .ok_or_else(|| IoError::Malformed { at, what: "face needs at least three vertices".into() })?;
    if triangles.len().checked_add(count).is_none_or(|n| n > MAX_ELEMENTS) {
        return Err(IoError::ResourceBound { what: "triangle count exceeds the element cap".into() });
    }
    if indices.len() > MAX_POLYGON_VERTICES {
        return Err(IoError::ResourceBound { what: "polygon vertex count exceeds the face cap".into() });
    }
    // Preserve the pre-existing triangular-soup admission boundary. Repair,
    // not a polygon triangulator, owns degenerate source triangles.
    let face = if indices.len() == 3 {
        vec![[indices[0], indices[1], indices[2]]]
    } else {
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

/// Export as OBJ (deterministic; f64 round-trip precision).
#[must_use]
pub fn write_obj(soup: &Soup) -> String {
    let mut out = String::with_capacity(soup.positions.len() * 32);
    for p in &soup.positions {
        let _ = writeln!(out, "v {} {} {}", p.x, p.y, p.z);
    }
    for t in &soup.triangles {
        let _ = writeln!(out, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_concave_faces_preserve_area_and_round_trip() {
        let source = "v 0 0 0\nv 3 0 0\nv 3 3 0\nv 2 3 0\nv 2 1 0\nv 1 1 0\nv 1 3 0\nv 0 3 0\nf -8/1/1 -7/2/1 -6/3/1 -5/4/1 -4/5/1 -3/6/1 -2/7/1 -1/8/1\n";
        let soup = read_obj(source).unwrap();
        assert_eq!(soup.triangles.len(), 6);
        let area: f64 = soup.triangles.iter().map(|t| {
            let [a, b, c] = t.map(|i| soup.positions[i as usize]);
            let area = ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)) * 0.5;
            assert!(area > 0.0);
            area
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
}
