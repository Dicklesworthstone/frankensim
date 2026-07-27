//! PLY import/export: ASCII and binary_little_endian, vertex x/y/z
//! (float/double) + face vertex-index lists. Other elements/properties
//! are skipped with correct stride accounting (binary) or token counting
//! (ASCII) — documented subset, structured rejection beyond it.

use crate::{IoError, MAX_ELEMENTS};
use fs_geom::Point3;
use fs_rep_mesh::Soup;
use std::fmt::Write as _;

const MAX_PLY_HEADER_LINES: usize = 10_000;

/// Maximum sum of all element row counts declared by one PLY header.
pub const MAX_PLY_TOTAL_ROWS: u64 = MAX_ELEMENTS as u64;

/// Maximum header-derived row/property dispatch work for one PLY body.
///
/// One unit is charged for entering an element row and one for visiting each
/// property declared on that row. List items are separately limited by
/// [`MAX_PLY_LIST_ITEMS`] and must each consume an input token or byte stride.
pub const MAX_PLY_DISPATCH_WORK: u64 = 4 * MAX_PLY_TOTAL_ROWS;

/// Maximum number of items admitted in one PLY list property.
pub const MAX_PLY_LIST_ITEMS: usize = 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    Ascii,
    BinaryLe,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    F32,
    F64,
}

impl Ty {
    fn parse(tok: &str) -> Option<Ty> {
        match tok {
            "char" | "int8" => Some(Ty::I8),
            "uchar" | "uint8" => Some(Ty::U8),
            "short" | "int16" => Some(Ty::I16),
            "ushort" | "uint16" => Some(Ty::U16),
            "int" | "int32" => Some(Ty::I32),
            "uint" | "uint32" => Some(Ty::U32),
            "float" | "float32" => Some(Ty::F32),
            "double" | "float64" => Some(Ty::F64),
            _ => None,
        }
    }

    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::F64 => 8,
        }
    }

    fn read_f64(self, b: &[u8]) -> f64 {
        match self {
            Ty::I8 => f64::from(b[0].cast_signed()),
            Ty::U8 => f64::from(b[0]),
            Ty::I16 => f64::from(i16::from_le_bytes([b[0], b[1]])),
            Ty::U16 => f64::from(u16::from_le_bytes([b[0], b[1]])),
            Ty::I32 => f64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            Ty::U32 => f64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            Ty::F32 => f64::from(f32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            Ty::F64 => f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]),
        }
    }
}

#[derive(Debug)]
enum Prop {
    Scalar(Ty, String),
    List(Ty, Ty, String),
}

#[derive(Debug)]
struct Element {
    name: String,
    count: usize,
    props: Vec<Prop>,
}

struct Header {
    format: Format,
    elements: Vec<Element>,
    body_start: usize,
}

fn parse_header(bytes: &[u8]) -> Result<Header, IoError> {
    let mut pos = 0usize;
    let mut lines = Vec::new();
    loop {
        let end = bytes[pos..]
            .iter()
            .position(|&b| b == b'\n')
            .ok_or(IoError::Malformed {
                at: pos,
                what: "header never ends".to_string(),
            })?;
        let line = core::str::from_utf8(&bytes[pos..pos + end])
            .map_err(|_| IoError::Malformed {
                at: pos,
                what: "non-UTF-8 header line".to_string(),
            })?
            .trim_end_matches('\r')
            .to_string();
        pos += end + 1;
        let stop = line == "end_header";
        lines.push(line);
        if stop {
            break;
        }
        if lines.len() > MAX_PLY_HEADER_LINES {
            return Err(IoError::ResourceBound {
                what: "header line cap".to_string(),
            });
        }
    }
    if lines.first().map(String::as_str) != Some("ply") {
        return Err(IoError::Malformed {
            at: 0,
            what: "missing 'ply' magic".to_string(),
        });
    }
    let mut format = None;
    let mut elements: Vec<Element> = Vec::new();
    for (ln, line) in lines.iter().enumerate() {
        header_line(line, ln, &mut format, &mut elements)?;
    }
    let format = format.ok_or(IoError::Malformed {
        at: 0,
        what: "no format line".to_string(),
    })?;
    preflight_body_dispatch(&elements)?;
    Ok(Header {
        format,
        elements,
        body_start: pos,
    })
}

fn preflight_body_dispatch(elements: &[Element]) -> Result<(), IoError> {
    let mut total_rows = 0u64;
    let mut dispatch_work = 0u64;
    for element in elements {
        if element.count != 0 && element.props.is_empty() {
            return Err(IoError::Unsupported {
                what: format!(
                    "PLY element {:?} declares {} rows but has no properties",
                    element.name, element.count
                ),
            });
        }

        let count = u64::try_from(element.count).map_err(|_| IoError::ResourceBound {
            what: "PLY aggregate element-row count overflow".to_string(),
        })?;
        let properties =
            u64::try_from(element.props.len()).map_err(|_| IoError::ResourceBound {
                what: "PLY aggregate property count overflow".to_string(),
            })?;
        total_rows = total_rows
            .checked_add(count)
            .ok_or(IoError::ResourceBound {
                what: "PLY aggregate element-row count overflow".to_string(),
            })?;
        let per_row = properties.checked_add(1).ok_or(IoError::ResourceBound {
            what: "PLY aggregate dispatch work overflow".to_string(),
        })?;
        let element_work = count.checked_mul(per_row).ok_or(IoError::ResourceBound {
            what: "PLY aggregate dispatch work overflow".to_string(),
        })?;
        dispatch_work = dispatch_work
            .checked_add(element_work)
            .ok_or(IoError::ResourceBound {
                what: "PLY aggregate dispatch work overflow".to_string(),
            })?;
    }

    if total_rows > MAX_PLY_TOTAL_ROWS {
        return Err(IoError::ResourceBound {
            what: format!(
                "PLY aggregate element rows {total_rows} exceed cap {MAX_PLY_TOTAL_ROWS}"
            ),
        });
    }
    if dispatch_work > MAX_PLY_DISPATCH_WORK {
        return Err(IoError::ResourceBound {
            what: format!(
                "PLY aggregate dispatch work {dispatch_work} exceeds cap \
                 {MAX_PLY_DISPATCH_WORK} (row iterations plus property visits)"
            ),
        });
    }
    Ok(())
}

fn header_line(
    line: &str,
    ln: usize,
    format: &mut Option<Format>,
    elements: &mut Vec<Element>,
) -> Result<(), IoError> {
    let mut it = line.split_whitespace();
    match it.next() {
        Some("format") => {
            *format = match it.next() {
                Some("ascii") => Some(Format::Ascii),
                Some("binary_little_endian") => Some(Format::BinaryLe),
                Some(other) => {
                    return Err(IoError::Unsupported {
                        what: format!("PLY format {other} (ascii/binary_little_endian only)"),
                    });
                }
                None => None,
            };
        }
        Some("element") => {
            let name = it.next().unwrap_or("").to_string();
            let count: usize =
                it.next()
                    .and_then(|t| t.parse().ok())
                    .ok_or(IoError::Malformed {
                        at: ln,
                        what: "element needs a count".to_string(),
                    })?;
            if count > MAX_ELEMENTS {
                return Err(IoError::ResourceBound {
                    what: format!("element {name} count {count} exceeds the cap"),
                });
            }
            elements.push(Element {
                name,
                count,
                props: Vec::new(),
            });
        }
        Some("property") => {
            let el = elements.last_mut().ok_or(IoError::Malformed {
                at: ln,
                what: "property before any element".to_string(),
            })?;
            let first = it.next().unwrap_or("");
            if first == "list" {
                let count_ty = Ty::parse(it.next().unwrap_or("")).ok_or(IoError::Unsupported {
                    what: "unknown list count type".to_string(),
                })?;
                let item_ty = Ty::parse(it.next().unwrap_or("")).ok_or(IoError::Unsupported {
                    what: "unknown list item type".to_string(),
                })?;
                let name = it.next().unwrap_or("").to_string();
                el.props.push(Prop::List(count_ty, item_ty, name));
            } else {
                let ty = Ty::parse(first).ok_or(IoError::Unsupported {
                    what: format!("unknown property type {first}"),
                })?;
                let name = it.next().unwrap_or("").to_string();
                el.props.push(Prop::Scalar(ty, name));
            }
        }
        _ => {} // comment / obj_info / ply / end_header
    }
    Ok(())
}

fn take<'a>(bytes: &'a [u8], pos: &mut usize, n: usize) -> Result<&'a [u8], IoError> {
    let out = bytes.get(*pos..*pos + n).ok_or(IoError::Malformed {
        at: *pos,
        what: "truncated body".to_string(),
    })?;
    *pos += n;
    Ok(out)
}

/// Import a PLY subset.
///
/// # Errors
/// [`IoError`] on malformed/unsupported/oversized input.
pub fn read_ply(bytes: &[u8]) -> Result<Soup, IoError> {
    let header = parse_header(bytes)?;
    let mut positions: Vec<Point3> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    match header.format {
        Format::Ascii => read_ascii_body(bytes, &header, &mut positions, &mut triangles)?,
        Format::BinaryLe => read_binary_body(bytes, &header, &mut positions, &mut triangles)?,
    }
    if positions.is_empty() || triangles.is_empty() {
        return Err(IoError::Malformed {
            at: 0,
            what: "PLY has no vertex/face payload in the supported subset".to_string(),
        });
    }
    // DEFERRED index validation (bead wqd.25.1): PLY headers define
    // element order, and a legal file may place faces BEFORE vertices —
    // validating during the parse rejected such files against a vertex
    // count of zero. Indices are checked here, against the FINAL count,
    // once every element has been consumed; `at` is the exact ordinal
    // of the offending triangle.
    let n_vertices = positions.len();
    for (ordinal, tri) in triangles.iter().enumerate() {
        for &i in tri {
            let in_range = usize::try_from(i).is_ok_and(|v| v < n_vertices);
            if !in_range {
                return Err(IoError::Malformed {
                    at: ordinal,
                    what: format!("face index {i} out of range (mesh has {n_vertices} vertices)"),
                });
            }
        }
    }
    Ok(Soup {
        positions,
        triangles,
    })
}

fn read_ascii_body(
    bytes: &[u8],
    header: &Header,
    positions: &mut Vec<Point3>,
    triangles: &mut Vec<[u32; 3]>,
) -> Result<(), IoError> {
    {
        {
            let text = core::str::from_utf8(&bytes[header.body_start..]).map_err(|e| {
                IoError::Malformed {
                    at: header.body_start + e.valid_up_to(),
                    what: "non-UTF-8 ASCII body".to_string(),
                }
            })?;
            let mut tokens = text.split_whitespace();
            for el in &header.elements {
                for _ in 0..el.count {
                    let mut xyz = [f64::NAN; 3];
                    for prop in &el.props {
                        match prop {
                            Prop::Scalar(_, name) => {
                                let v = next_f64_token(next_token(&mut tokens, name)?, name)?;
                                match name.as_str() {
                                    "x" => xyz[0] = v,
                                    "y" => xyz[1] = v,
                                    "z" => xyz[2] = v,
                                    _ => {}
                                }
                            }
                            Prop::List(_, _, name) => {
                                let count_tok = next_token(&mut tokens, "list count")?;
                                let n = parse_usize_token(count_tok, "list count")?;
                                if n > MAX_PLY_LIST_ITEMS {
                                    return Err(IoError::ResourceBound {
                                        what: format!("list longer than {MAX_PLY_LIST_ITEMS}"),
                                    });
                                }
                                if el.name == "face"
                                    && (name == "vertex_indices" || name == "vertex_index")
                                {
                                    let mut idx = Vec::with_capacity(n);
                                    for _ in 0..n {
                                        let item_tok = next_token(&mut tokens, "list item")?;
                                        idx.push(parse_u32_token(item_tok, "list item")?);
                                    }
                                    push_face(triangles, &idx)?;
                                } else {
                                    for _ in 0..n {
                                        let item_tok = next_token(&mut tokens, "list item")?;
                                        let _ = next_f64_token(item_tok, "list item")?;
                                    }
                                }
                            }
                        }
                    }
                    if el.name == "vertex" {
                        push_vertex(positions, xyz)?;
                    }
                }
            }
        }
    }
    Ok(())
}

fn read_binary_body(
    bytes: &[u8],
    header: &Header,
    positions: &mut Vec<Point3>,
    triangles: &mut Vec<[u32; 3]>,
) -> Result<(), IoError> {
    {
        {
            let mut pos = header.body_start;
            for el in &header.elements {
                for _ in 0..el.count {
                    let mut xyz = [f64::NAN; 3];
                    for prop in &el.props {
                        match prop {
                            Prop::Scalar(ty, name) => {
                                let b = take(bytes, &mut pos, ty.size())?;
                                let v = ty.read_f64(b);
                                match name.as_str() {
                                    "x" => xyz[0] = v,
                                    "y" => xyz[1] = v,
                                    "z" => xyz[2] = v,
                                    _ => {}
                                }
                            }
                            Prop::List(count_ty, item_ty, name) => {
                                let cb = take(bytes, &mut pos, count_ty.size())?;
                                let n = parse_usize_value(count_ty.read_f64(cb), "list count")?;
                                if n > MAX_PLY_LIST_ITEMS {
                                    return Err(IoError::ResourceBound {
                                        what: format!("list longer than {MAX_PLY_LIST_ITEMS}"),
                                    });
                                }
                                if el.name == "face"
                                    && (name == "vertex_indices" || name == "vertex_index")
                                {
                                    let mut idx = Vec::with_capacity(n);
                                    for _ in 0..n {
                                        let ib = take(bytes, &mut pos, item_ty.size())?;
                                        idx.push(parse_u32_value(
                                            item_ty.read_f64(ib),
                                            "list item",
                                        )?);
                                    }
                                    push_face(triangles, &idx)?;
                                } else {
                                    for _ in 0..n {
                                        let _ = take(bytes, &mut pos, item_ty.size())?;
                                    }
                                }
                            }
                        }
                    }
                    if el.name == "vertex" {
                        push_vertex(positions, xyz)?;
                    }
                }
            }
        }
    }
    Ok(())
}

fn next_token<'a>(
    tokens: &mut core::str::SplitWhitespace<'a>,
    what: &str,
) -> Result<&'a str, IoError> {
    tokens.next().ok_or(IoError::Malformed {
        at: 0,
        what: format!("body ended early ({what})"),
    })
}

fn next_f64_token(tok: &str, what: &str) -> Result<f64, IoError> {
    tok.parse().map_err(|_| IoError::Malformed {
        at: 0,
        what: format!("bad number {tok:?} ({what})"),
    })
}

fn parse_usize_token(tok: &str, what: &str) -> Result<usize, IoError> {
    let v = next_f64_token(tok, what)?;
    parse_usize_value(v, what)
}

fn parse_u32_token(tok: &str, what: &str) -> Result<u32, IoError> {
    let v = next_f64_token(tok, what)?;
    parse_u32_value(v, what)
}

#[allow(clippy::cast_precision_loss)]
fn parse_usize_value(v: f64, what: &str) -> Result<usize, IoError> {
    if !(v.is_finite() && v.fract() == 0.0 && v >= 0.0) {
        return Err(IoError::Malformed {
            at: 0,
            what: format!("{what} must be a non-negative integer, got {v}"),
        });
    }
    if v > MAX_PLY_LIST_ITEMS as f64 {
        return MAX_PLY_LIST_ITEMS
            .checked_add(1)
            .ok_or_else(|| IoError::ResourceBound {
                what: "PLY list-item cap has no representable refusal sentinel".to_string(),
            });
    }
    Ok(v as usize)
}

fn parse_u32_value(v: f64, what: &str) -> Result<u32, IoError> {
    if !(v.is_finite() && v.fract() == 0.0 && v >= 0.0 && v <= f64::from(u32::MAX)) {
        return Err(IoError::Malformed {
            at: 0,
            what: format!("{what} must be a non-negative u32, got {v}"),
        });
    }
    Ok(v as u32)
}

fn push_vertex(positions: &mut Vec<Point3>, xyz: [f64; 3]) -> Result<(), IoError> {
    if xyz.iter().any(|v| !v.is_finite()) {
        return Err(IoError::Malformed {
            at: positions.len(),
            what: "vertex missing x/y/z or non-finite".to_string(),
        });
    }
    positions.push(Point3::new(xyz[0], xyz[1], xyz[2]));
    Ok(())
}

/// Triangulate one face fan into pending triangles. Index RANGE checks
/// are deferred to [`read_ply`]'s post-parse pass (legal PLY may put
/// faces before vertices); structural checks and the resource cap stay
/// here.
fn push_face(triangles: &mut Vec<[u32; 3]>, idx: &[u32]) -> Result<(), IoError> {
    if idx.len() < 3 {
        return Err(IoError::Malformed {
            at: triangles.len(),
            what: "face with fewer than three indices".to_string(),
        });
    }
    for k in 1..idx.len() - 1 {
        triangles.push([idx[0], idx[k], idx[k + 1]]);
        if triangles.len() > MAX_ELEMENTS {
            return Err(IoError::ResourceBound {
                what: "triangle cap".to_string(),
            });
        }
    }
    Ok(())
}

/// Export as deterministic ASCII PLY with shortest-round-trip f64 text.
///
/// When a soup has at least one vertex and one triangle, every position
/// coordinate is finite, every triangle index is in range, and its counts fit
/// [`read_ply`]'s documented envelopes, re-import preserves position bits and
/// triangle indices exactly. This infallible writer does not itself validate
/// those preconditions.
#[must_use]
pub fn write_ply(soup: &Soup) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "ply\nformat ascii 1.0");
    let _ = writeln!(out, "element vertex {}", soup.positions.len());
    let _ = writeln!(
        out,
        "property double x\nproperty double y\nproperty double z"
    );
    let _ = writeln!(out, "element face {}", soup.triangles.len());
    let _ = writeln!(out, "property list uchar uint vertex_indices\nend_header");
    for p in &soup.positions {
        let _ = writeln!(out, "{} {} {}", p.x, p.y, p.z);
    }
    for t in &soup.triangles {
        let _ = writeln!(out, "3 {} {} {}", t[0], t[1], t[2]);
    }
    out
}

#[cfg(all(test, target_pointer_width = "64"))]
mod tests {
    use super::{Element, Prop, Ty, preflight_body_dispatch};
    use crate::IoError;

    #[test]
    fn aggregate_dispatch_arithmetic_refuses_overflow() {
        let elements = [Element {
            name: "synthetic-overflow-probe".to_string(),
            count: usize::MAX,
            props: vec![Prop::Scalar(Ty::U8, "value".to_string())],
        }];
        assert_eq!(
            preflight_body_dispatch(&elements),
            Err(IoError::ResourceBound {
                what: "PLY aggregate dispatch work overflow".to_string(),
            })
        );
    }
}
