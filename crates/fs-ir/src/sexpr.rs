//! Canonical s-expression syntax (plan §11.1): total hand-rolled reader
//! with byte-exact spans, and the canonical printer whose output reparses
//! to the same shape (round-trip law, conformance-tested).
//!
//! Atom classification (deterministic, no silent fallbacks): `:name` →
//! keyword; `0x…` → seed; number-leading tokens are numeric and MUST parse
//! as int, float, fs-qty quantity, or count — a number with a garbage
//! suffix is a structured error pointing at its span, never a symbol;
//! `"…"` → string (with `\" \\ \n \t` escapes); anything else → symbol.
//! Comments run `;` to end of line. Recursion depth is capped: adversarial
//! nesting is a structured rejection, not a stack overflow (fuzz law).

use crate::ast::{CountUnit, Node, NodeKind, Span};
use crate::{IrError, IrErrorKind};

/// Maximum nesting depth (structured rejection beyond — G0 fuzz law).
pub const MAX_DEPTH: usize = 256;

/// Parse one s-expression program (exactly one top-level form).
///
/// # Errors
/// Structured [`IrError`] with the offending span and a fix hint.
pub fn parse(src: &str) -> Result<Node, IrError> {
    let mut p = Parser {
        src,
        bytes: src.as_bytes(),
        pos: 0,
    };
    p.skip_trivia();
    let node = p.parse_node(0)?;
    p.skip_trivia();
    if p.pos != p.bytes.len() {
        return Err(IrError {
            span: Span::new(p.pos, p.bytes.len()),
            kind: IrErrorKind::TrailingInput,
            detail: "input continues after the top-level form".to_string(),
            hint: "a program is exactly one form; wrap multiple forms in (study ...)".to_string(),
        });
    }
    Ok(node)
}

struct Parser<'a> {
    src: &'a str,
    bytes: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn skip_trivia(&mut self) {
        while self.pos < self.bytes.len() {
            match self.bytes[self.pos] {
                b' ' | b'\t' | b'\r' | b'\n' => self.pos += 1,
                b';' => {
                    while self.pos < self.bytes.len() && self.bytes[self.pos] != b'\n' {
                        self.pos += 1;
                    }
                }
                _ => break,
            }
        }
    }

    fn parse_node(&mut self, depth: usize) -> Result<Node, IrError> {
        if depth > MAX_DEPTH {
            return Err(IrError {
                span: Span::new(self.pos, self.pos + 1),
                kind: IrErrorKind::TooDeep,
                detail: format!("nesting exceeds the {MAX_DEPTH}-level cap"),
                hint: "flatten the program; adversarial nesting is refused by design".to_string(),
            });
        }
        self.skip_trivia();
        let Some(&c) = self.bytes.get(self.pos) else {
            return Err(IrError {
                span: Span::new(self.pos, self.pos),
                kind: IrErrorKind::UnexpectedEnd,
                detail: "expected a form or atom, found end of input".to_string(),
                hint: "supply a complete s-expression, e.g. (study \"name\" ...)".to_string(),
            });
        };
        match c {
            b'(' => self.parse_list(depth),
            b')' => Err(IrError {
                span: Span::new(self.pos, self.pos + 1),
                kind: IrErrorKind::UnexpectedChar,
                detail: "unmatched closing paren".to_string(),
                hint: "remove it or open a matching '(' earlier".to_string(),
            }),
            b'"' => self.parse_string(),
            _ => self.parse_atom(),
        }
    }

    fn parse_list(&mut self, depth: usize) -> Result<Node, IrError> {
        let start = self.pos;
        self.pos += 1; // consume '('
        let mut items = Vec::new();
        loop {
            self.skip_trivia();
            match self.bytes.get(self.pos) {
                None => {
                    return Err(IrError {
                        span: Span::new(start, self.pos),
                        kind: IrErrorKind::UnclosedParen,
                        detail: "this form is never closed".to_string(),
                        hint: "add the matching ')'".to_string(),
                    });
                }
                Some(b')') => {
                    self.pos += 1;
                    return Ok(Node {
                        kind: NodeKind::List(items),
                        span: Span::new(start, self.pos),
                    });
                }
                Some(_) => items.push(self.parse_node(depth + 1)?),
            }
        }
    }

    fn parse_string(&mut self) -> Result<Node, IrError> {
        let start = self.pos;
        self.pos += 1; // consume '"'
        let mut out = String::new();
        while let Some(&c) = self.bytes.get(self.pos) {
            match c {
                b'"' => {
                    self.pos += 1;
                    return Ok(Node {
                        kind: NodeKind::Str(out),
                        span: Span::new(start, self.pos),
                    });
                }
                b'\\' => {
                    let esc = self.bytes.get(self.pos + 1).copied();
                    let ch = match esc {
                        Some(b'"') => '"',
                        Some(b'\\') => '\\',
                        Some(b'n') => '\n',
                        Some(b't') => '\t',
                        _ => {
                            return Err(IrError {
                                span: Span::new(self.pos, (self.pos + 2).min(self.bytes.len())),
                                kind: IrErrorKind::BadEscape,
                                detail: "unknown escape sequence".to_string(),
                                hint: r#"supported escapes: \" \\ \n \t"#.to_string(),
                            });
                        }
                    };
                    out.push(ch);
                    self.pos += 2;
                }
                _ => {
                    // Advance one UTF-8 scalar (source is &str: boundaries ok).
                    let ch_len = self.src[self.pos..]
                        .chars()
                        .next()
                        .map_or(1, char::len_utf8);
                    out.push_str(&self.src[self.pos..self.pos + ch_len]);
                    self.pos += ch_len;
                }
            }
        }
        Err(IrError {
            span: Span::new(start, self.pos),
            kind: IrErrorKind::UnclosedString,
            detail: "string literal is never closed".to_string(),
            hint: "add the closing '\"'".to_string(),
        })
    }

    fn parse_atom(&mut self) -> Result<Node, IrError> {
        let start = self.pos;
        while let Some(&c) = self.bytes.get(self.pos) {
            if c.is_ascii_whitespace() || c == b'(' || c == b')' || c == b'"' || c == b';' {
                break;
            }
            self.pos += 1;
        }
        let span = Span::new(start, self.pos);
        let text = &self.src[start..self.pos];
        classify_atom(text, span)
    }
}

/// Classify one atom token (shared with the JSON side's qty re-parsing).
pub(crate) fn classify_atom(text: &str, span: Span) -> Result<Node, IrError> {
    if text.is_empty() {
        // Reachable when the scanner meets a byte only the atom reader
        // treats as a delimiter (fuzz-found): refuse, never panic (P10).
        return Err(IrError {
            span,
            kind: IrErrorKind::UnexpectedChar,
            detail: "empty token (stray delimiter byte)".to_string(),
            hint: "remove the stray character".to_string(),
        });
    }
    if let Some(name) = text.strip_prefix(':') {
        if name.is_empty() {
            return Err(IrError {
                span,
                kind: IrErrorKind::BadKeyword,
                detail: "bare ':' is not a keyword".to_string(),
                hint: "keywords are :name".to_string(),
            });
        }
        return Ok(Node {
            kind: NodeKind::Keyword(name.to_string()),
            span,
        });
    }
    if let Some(hex) = text.strip_prefix("0x").or_else(|| text.strip_prefix("0X")) {
        return match u64::from_str_radix(hex, 16) {
            Ok(v) => Ok(Node {
                kind: NodeKind::Seed(v),
                span,
            }),
            Err(_) => Err(IrError {
                span,
                kind: IrErrorKind::BadSeed,
                detail: format!("{text:?} is not a valid 64-bit hex seed"),
                hint: "seeds are 0x-prefixed hex fitting u64, e.g. 0x5EED0001".to_string(),
            }),
        };
    }
    let b = text.as_bytes();
    let digit_at = |i: usize| b.get(i).is_some_and(u8::is_ascii_digit);
    // Numeric lead: 5, .5, -5, +5, -.5, +.5 — but NOT "..", "-", "+", ".".
    let numeric_lead = digit_at(0)
        || ((b[0] == b'-' || b[0] == b'+')
            && (digit_at(1) || (b.get(1) == Some(&b'.') && digit_at(2))))
        || (b[0] == b'.' && digit_at(1));
    if numeric_lead {
        return classify_numeric(text, span);
    }
    Ok(Node {
        kind: NodeKind::Symbol(text.to_string()),
        span,
    })
}

fn classify_numeric(text: &str, span: Span) -> Result<Node, IrError> {
    if let Ok(i) = text.parse::<i64>() {
        return Ok(Node {
            kind: NodeKind::Int(i),
            span,
        });
    }
    if let Ok(f) = text.parse::<f64>() {
        if !f.is_finite() {
            return Err(IrError {
                span,
                kind: IrErrorKind::BadNumber,
                detail: format!("{text:?} is not finite"),
                hint: "IR literals must be finite".to_string(),
            });
        }
        return Ok(Node {
            kind: NodeKind::Float(f),
            span,
        });
    }
    // Count units (information/core grants) before SI quantities: fs-qty
    // refuses information units by design.
    let split = text
        .char_indices()
        .find(|(_, c)| c.is_ascii_alphabetic())
        .map(|(i, _)| i);
    if let Some(i) = split
        && let Some(unit) = CountUnit::from_suffix(&text[i..])
        && let Ok(v) = text[..i].parse::<f64>()
        && v.is_finite()
    {
        return Ok(Node {
            kind: NodeKind::Count { value: v, unit },
            span,
        });
    }
    match fs_qty::parse::parse_qty(text) {
        Ok(q) => Ok(Node {
            kind: NodeKind::Qty {
                value: q.value,
                dims: q.dims,
                text: text.to_string(),
            },
            span,
        }),
        Err(e) => Err(IrError {
            span,
            kind: IrErrorKind::BadQuantity,
            detail: format!("{text:?} is not an int, float, quantity, or count: {e}"),
            hint: "numeric tokens must fully parse; e.g. 0.12Pa*s, 65deg, 384GiB, 2e-2".to_string(),
        }),
    }
}

/// Print a node in canonical s-expression form; `parse(print(x))` has the
/// same shape as `x` (round-trip law).
#[must_use]
pub fn print(node: &Node) -> String {
    let mut out = String::new();
    print_into(node, &mut out);
    out
}

fn print_into(node: &Node, out: &mut String) {
    use std::fmt::Write as _;
    match &node.kind {
        NodeKind::Int(i) => {
            let _ = write!(out, "{i}");
        }
        NodeKind::Float(f) => {
            let _ = write!(out, "{f:?}");
        }
        NodeKind::Qty { text, .. } => out.push_str(text),
        NodeKind::Count { value, unit } => {
            let _ = write!(out, "{value:?}{}", unit.suffix());
        }
        NodeKind::Seed(v) => {
            let _ = write!(out, "0x{v:X}");
        }
        NodeKind::Str(s) => {
            out.push('"');
            for c in s.chars() {
                match c {
                    '"' => out.push_str("\\\""),
                    '\\' => out.push_str("\\\\"),
                    '\n' => out.push_str("\\n"),
                    '\t' => out.push_str("\\t"),
                    _ => out.push(c),
                }
            }
            out.push('"');
        }
        NodeKind::Symbol(s) => out.push_str(s),
        NodeKind::Keyword(k) => {
            out.push(':');
            out.push_str(k);
        }
        NodeKind::List(items) => {
            out.push('(');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(' ');
                }
                print_into(item, out);
            }
            out.push(')');
        }
    }
}
