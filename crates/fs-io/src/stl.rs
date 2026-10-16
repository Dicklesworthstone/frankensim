//! STL import/export: binary (auto-detected) and ASCII. Import is
//! hostile-input hardened (length-checked, element-capped); export is
//! byte-deterministic. Fidelity contract: STL carries positions only —
//! connectivity is reconstructed by exact coordinate matching, normals
//! are recomputed (documented lossy).

use crate::{IoError, MAX_ELEMENTS};
use fs_geom::Point3;
use fs_rep_mesh::Soup;
use std::collections::BTreeMap;

fn f32_le(bytes: &[u8], at: usize) -> Result<f32, IoError> {
    bytes
        .get(at..at + 4)
        .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        .ok_or(IoError::Malformed {
            at,
            what: "truncated float".to_string(),
        })
}

/// Weld exactly-equal vertex positions into an indexed soup (STL stores
/// per-facet corners; bitwise-equal coordinates weld — no tolerance).
fn weld(corners: &[[f64; 3]]) -> Soup {
    let mut index: BTreeMap<[u64; 3], u32> = BTreeMap::new();
    let mut positions = Vec::new();
    let mut triangles = Vec::new();
    let (tris, _rest) = corners.as_chunks::<3>();
    for tri in tris {
        let mut ids = [0u32; 3];
        for (slot, c) in ids.iter_mut().zip(tri) {
            let key = [c[0].to_bits(), c[1].to_bits(), c[2].to_bits()];
            *slot = *index.entry(key).or_insert_with(|| {
                positions.push(Point3::new(c[0], c[1], c[2]));
                u32::try_from(positions.len() - 1).expect("element cap")
            });
        }
        triangles.push(ids);
    }
    Soup {
        positions,
        triangles,
    }
}

fn read_binary(bytes: &[u8]) -> Result<Soup, IoError> {
    if bytes.len() < 84 {
        return Err(IoError::Malformed {
            at: bytes.len(),
            what: "binary STL shorter than header + count".to_string(),
        });
    }
    let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
    if count > MAX_ELEMENTS {
        return Err(IoError::ResourceBound {
            what: format!("{count} facets exceeds the element cap"),
        });
    }
    let need = 84 + count * 50;
    if bytes.len() < need {
        return Err(IoError::Malformed {
            at: bytes.len(),
            what: format!("binary STL declares {count} facets but is truncated"),
        });
    }
    let mut corners = Vec::with_capacity(count * 3);
    for f in 0..count {
        let base = 84 + f * 50 + 12; // skip the stored normal
        for corner in 0..3 {
            let at = base + corner * 12;
            let x = f32_le(bytes, at)?;
            let y = f32_le(bytes, at + 4)?;
            let z = f32_le(bytes, at + 8)?;
            if !(x.is_finite() && y.is_finite() && z.is_finite()) {
                return Err(IoError::Malformed {
                    at,
                    what: format!("non-finite vertex in facet {f}"),
                });
            }
            // f32 -> f64 is exact and injective: binary welding is unchanged.
            corners.push([f64::from(x), f64::from(y), f64::from(z)]);
        }
    }
    Ok(weld(&corners))
}

fn read_ascii(text: &str) -> Result<Soup, IoError> {
    let mut corners: Vec<[f64; 3]> = Vec::new();
    for (ln, line) in text.lines().enumerate() {
        let mut it = line.split_whitespace();
        if it.next() == Some("vertex") {
            let mut v = [0.0f64; 3];
            for slot in &mut v {
                let tok = it.next().ok_or(IoError::Malformed {
                    at: ln + 1,
                    what: "vertex needs three coordinates".to_string(),
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
            corners.push(v);
            if corners.len() > 3 * MAX_ELEMENTS {
                return Err(IoError::ResourceBound {
                    what: "vertex count exceeds the element cap".to_string(),
                });
            }
        }
    }
    if corners.is_empty() || !corners.len().is_multiple_of(3) {
        return Err(IoError::Malformed {
            at: 0,
            what: format!(
                "ASCII STL has {} vertices (must be a positive multiple of 3)",
                corners.len()
            ),
        });
    }
    Ok(weld(&corners))
}

/// Import an STL (binary auto-detected by the `solid` prefix heuristic
/// PLUS a length check — binary files that begin with "solid" are
/// handled correctly).
///
/// # Errors
/// [`IoError`] on malformed/oversized input.
pub fn read_stl(bytes: &[u8]) -> Result<Soup, IoError> {
    if bytes.len() >= 84 {
        let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        if count <= MAX_ELEMENTS && bytes.len() == 84 + count * 50 {
            return read_binary(bytes);
        }
    }
    let text = core::str::from_utf8(bytes).map_err(|e| IoError::Malformed {
        at: e.valid_up_to(),
        what: "neither a well-sized binary STL nor UTF-8 text".to_string(),
    })?;
    if text.trim_start().starts_with("solid") {
        read_ascii(text)
    } else {
        Err(IoError::Malformed {
            at: 0,
            what: "not an STL (no binary sizing, no 'solid' header)".to_string(),
        })
    }
}

/// Export a soup as binary STL (deterministic bytes; normals computed
/// from winding).
#[must_use]
pub fn write_stl(soup: &Soup) -> Vec<u8> {
    let mut out = vec![0u8; 80];
    out.extend_from_slice(
        &u32::try_from(soup.triangles.len())
            .expect("element cap")
            .to_le_bytes(),
    );
    for t in 0..soup.triangles.len() {
        let [a, b, c] = soup.tri(t);
        let u = [b.x - a.x, b.y - a.y, b.z - a.z];
        let v = [c.x - a.x, c.y - a.y, c.z - a.z];
        let mut n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let len = fs_math::det::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if len > 0.0 {
            for x in &mut n {
                *x /= len;
            }
        }
        for x in n {
            #[allow(clippy::cast_possible_truncation)]
            out.extend_from_slice(&(x as f32).to_le_bytes());
        }
        for p in [a, b, c] {
            #[allow(clippy::cast_possible_truncation)]
            for x in [p.x, p.y, p.z] {
                out.extend_from_slice(&(x as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&[0, 0]); // attribute byte count
    }
    out
}
