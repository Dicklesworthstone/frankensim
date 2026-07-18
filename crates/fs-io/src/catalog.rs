//! CSV/JSON catalog ingestion with SCHEMA VALIDATION — the AISC-catalog
//! path for the frame flagship. Every cell is checked against a declared
//! column spec; violations are helpful errors naming the row, column,
//! offending text, and the expectation. Quoted CSV fields (RFC-4180
//! subset with escaped quotes) are supported.

use crate::IoError;
use std::collections::BTreeMap;

/// What a column must contain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColumnKind {
    /// Any nonempty string.
    Text,
    /// A finite float, optionally bounded.
    Number {
        /// Inclusive lower bound.
        min: f64,
        /// Inclusive upper bound.
        max: f64,
    },
}

/// One column's contract.
#[derive(Debug, Clone)]
pub struct ColumnSpec {
    /// Column name (must appear in the header).
    pub name: &'static str,
    /// The value contract.
    pub kind: ColumnKind,
    /// Whether empty cells are allowed.
    pub required: bool,
}

/// A catalog schema.
#[derive(Debug, Clone)]
pub struct Schema {
    /// Column contracts (order-independent; matched by header name).
    pub columns: Vec<ColumnSpec>,
}

/// A validated catalog: rows of (column name → text) with numbers
/// pre-parsed where the schema demands them.
#[derive(Debug, Clone, PartialEq)]
pub struct Catalog {
    /// Row-major cells keyed by column name.
    pub rows: Vec<BTreeMap<String, String>>,
    /// Pre-parsed numeric views for Number columns.
    pub numbers: Vec<BTreeMap<String, f64>>,
}

/// Split one CSV record (RFC-4180 subset: quoted fields, `""` escapes).
fn split_csv(line: &str, row: usize) -> Result<Vec<String>, IoError> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut chars = line.chars().peekable();
    let mut quoted = false;
    while let Some(c) = chars.next() {
        if quoted {
            match c {
                '"' if chars.peek() == Some(&'"') => {
                    cur.push('"');
                    chars.next();
                }
                '"' => quoted = false,
                other => cur.push(other),
            }
        } else {
            match c {
                '"' if cur.is_empty() => quoted = true,
                ',' => fields.push(core::mem::take(&mut cur)),
                other => cur.push(other),
            }
        }
    }
    if quoted {
        return Err(IoError::Malformed {
            at: row,
            what: "unterminated quoted CSV field".to_string(),
        });
    }
    fields.push(cur);
    Ok(fields)
}

impl Schema {
    /// Validate one cell against its spec.
    fn check_cell(spec: &ColumnSpec, text: &str, row: usize) -> Result<Option<f64>, IoError> {
        let trimmed = text.trim();
        if trimmed.is_empty() {
            if spec.required {
                return Err(IoError::Schema {
                    row,
                    column: spec.name.to_string(),
                    what: "required cell is empty".to_string(),
                });
            }
            return Ok(None);
        }
        match spec.kind {
            ColumnKind::Text => Ok(None),
            ColumnKind::Number { min, max } => {
                let v: f64 = trimmed.parse().map_err(|_| IoError::Schema {
                    row,
                    column: spec.name.to_string(),
                    what: format!("{trimmed:?} is not a number"),
                })?;
                if !v.is_finite() || v < min || v > max {
                    return Err(IoError::Schema {
                        row,
                        column: spec.name.to_string(),
                        what: format!("{v} outside the declared range [{min}, {max}]"),
                    });
                }
                Ok(Some(v))
            }
        }
    }

    /// Parse + validate a CSV catalog (first record is the header).
    ///
    /// # Errors
    /// [`IoError::Malformed`] for CSV structure; [`IoError::Schema`] with
    /// row/column/expectation for value violations; missing schema
    /// columns are named.
    pub fn parse_csv(&self, text: &str) -> Result<Catalog, IoError> {
        let mut lines = text
            .lines()
            .enumerate()
            .filter(|(_, l)| !l.trim().is_empty());
        let (_, header_line) = lines.next().ok_or(IoError::Malformed {
            at: 0,
            what: "empty catalog".to_string(),
        })?;
        let header = split_csv(header_line, 0)?;
        for spec in &self.columns {
            if !header.iter().any(|h| h.trim() == spec.name) {
                return Err(IoError::Schema {
                    row: 0,
                    column: spec.name.to_string(),
                    what: format!(
                        "column missing from the header (found: {})",
                        header.join(", ")
                    ),
                });
            }
        }
        let mut rows = Vec::new();
        let mut numbers = Vec::new();
        for (data_row, (ln, line)) in lines.enumerate() {
            let row_no = data_row + 1; // 1-based, header excluded
            let fields = split_csv(line, row_no)?;
            if fields.len() != header.len() {
                return Err(IoError::Malformed {
                    at: ln + 1,
                    what: format!(
                        "record has {} fields, header has {}",
                        fields.len(),
                        header.len()
                    ),
                });
            }
            let mut row = BTreeMap::new();
            let mut nums = BTreeMap::new();
            for (h, cell) in header.iter().zip(&fields) {
                row.insert(h.trim().to_string(), cell.clone());
            }
            for spec in &self.columns {
                let cell = row.get(spec.name).cloned().unwrap_or_default();
                if let Some(v) = Self::check_cell(spec, &cell, row_no)? {
                    nums.insert(spec.name.to_string(), v);
                }
            }
            rows.push(row);
            numbers.push(nums);
        }
        Ok(Catalog { rows, numbers })
    }

    /// Parse + validate a JSON catalog: an array of flat objects
    /// (string/number members). Uses a minimal in-house JSON reader —
    /// structured rejection on anything outside that subset.
    ///
    /// # Errors
    /// [`IoError`] for JSON structure or schema violations.
    pub fn parse_json(&self, text: &str) -> Result<Catalog, IoError> {
        let rows_raw = mini_json_array_of_objects(text)?;
        let mut rows = Vec::new();
        let mut numbers = Vec::new();
        for (i, obj) in rows_raw.iter().enumerate() {
            let mut nums = BTreeMap::new();
            for spec in &self.columns {
                let cell = obj.get(spec.name).cloned().unwrap_or_default();
                if let Some(v) = Self::check_cell(spec, &cell, i + 1)? {
                    nums.insert(spec.name.to_string(), v);
                }
            }
            rows.push(obj.clone());
            numbers.push(nums);
        }
        Ok(Catalog { rows, numbers })
    }
}

/// A deliberately tiny JSON reader: `[ {"k": "v" | number, ...}, ... ]`.
fn mini_json_array_of_objects(text: &str) -> Result<Vec<BTreeMap<String, String>>, IoError> {
    let bytes = text.as_bytes();
    let mut pos = 0usize;
    let err = |pos: usize, what: &str| IoError::Malformed {
        at: pos,
        what: what.to_string(),
    };
    let skip_ws = |pos: &mut usize| {
        while *pos < bytes.len() && matches!(bytes[*pos], b' ' | b'\n' | b'\r' | b'\t') {
            *pos += 1;
        }
    };
    skip_ws(&mut pos);
    if bytes.get(pos) != Some(&b'[') {
        return Err(err(pos, "expected a JSON array"));
    }
    pos += 1;
    let mut out = Vec::new();
    loop {
        skip_ws(&mut pos);
        match bytes.get(pos) {
            Some(&b']') => {
                pos += 1;
                skip_ws(&mut pos);
                if pos != bytes.len() {
                    return Err(err(pos, "trailing bytes after the JSON array"));
                }
                break;
            }
            Some(&b'{') => {
                pos += 1;
                let mut obj = BTreeMap::new();
                loop {
                    skip_ws(&mut pos);
                    match bytes.get(pos) {
                        Some(&b'}') => {
                            pos += 1;
                            break;
                        }
                        Some(&b'"') => {
                            let key = read_string(bytes, &mut pos)?;
                            skip_ws(&mut pos);
                            if bytes.get(pos) != Some(&b':') {
                                return Err(err(pos, "expected ':'"));
                            }
                            pos += 1;
                            skip_ws(&mut pos);
                            let value = match bytes.get(pos) {
                                Some(&b'"') => read_string(bytes, &mut pos)?,
                                Some(c) if c.is_ascii_digit() || *c == b'-' || *c == b'+' => {
                                    let start = pos;
                                    while pos < bytes.len()
                                        && (bytes[pos].is_ascii_digit()
                                            || matches!(
                                                bytes[pos],
                                                b'.' | b'e' | b'E' | b'-' | b'+'
                                            ))
                                    {
                                        pos += 1;
                                    }
                                    core::str::from_utf8(&bytes[start..pos])
                                        .map_err(|_| err(start, "bad number"))?
                                        .to_string()
                                }
                                _ => return Err(err(pos, "expected a string or number value")),
                            };
                            obj.insert(key, value);
                            skip_ws(&mut pos);
                            if bytes.get(pos) == Some(&b',') {
                                pos += 1;
                            }
                        }
                        _ => return Err(err(pos, "expected a key or '}'")),
                    }
                }
                out.push(obj);
                skip_ws(&mut pos);
                if bytes.get(pos) == Some(&b',') {
                    pos += 1;
                }
            }
            _ => return Err(err(pos, "expected an object or ']'")),
        }
        if out.len() > crate::MAX_ELEMENTS {
            return Err(IoError::ResourceBound {
                what: "catalog row cap".to_string(),
            });
        }
    }
    Ok(out)
}

fn read_string(bytes: &[u8], pos: &mut usize) -> Result<String, IoError> {
    *pos += 1; // opening quote
    // Accumulate raw bytes and decode as UTF-8 at the close, so multi-byte
    // scalars (e.g. "café") survive instead of being split into Latin-1 chars.
    let mut buf: Vec<u8> = Vec::new();
    loop {
        match bytes.get(*pos) {
            None => {
                return Err(IoError::Malformed {
                    at: *pos,
                    what: "unterminated string".to_string(),
                });
            }
            Some(&b'"') => {
                *pos += 1;
                return String::from_utf8(buf).map_err(|_| IoError::Malformed {
                    at: *pos,
                    what: "invalid UTF-8 in string".to_string(),
                });
            }
            Some(&b'\\') => {
                let escaped = bytes.get(*pos + 1).ok_or(IoError::Malformed {
                    at: *pos,
                    what: "dangling escape".to_string(),
                })?;
                buf.push(match escaped {
                    b'n' => b'\n',
                    b't' => b'\t',
                    other => *other,
                });
                *pos += 2;
            }
            Some(&c) => {
                buf.push(c);
                *pos += 1;
            }
        }
    }
}
