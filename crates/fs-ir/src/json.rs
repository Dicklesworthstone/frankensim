//! The lossless JSON mapping (plan §11.1): the second concrete syntax,
//! isomorphic to the s-expressions — both parse to the same typed AST
//! (property-tested, not aspirational).
//!
//! Mapping (single-key tagged objects; unknown tags are structured
//! rejections):
//!
//! | AST          | JSON                          |
//! |--------------|-------------------------------|
//! | `Int`        | `{"i": 42}`                   |
//! | `Float`      | `{"f": 2e-2}`                 |
//! | `Qty`        | `{"q": "0.12Pa*s"}` (literal) |
//! | `Count`      | `{"c": "384GiB"}` (literal)   |
//! | `Seed`       | `{"seed": "0xF00D0002"}`      |
//! | `Str`        | `{"s": "..."}`                |
//! | `Symbol`     | `{"sym": "..."}`              |
//! | `Keyword`    | `{"kw": "name"}`              |
//! | `List`       | `[...]`                       |
//!
//! Qty/Count/Seed reuse the s-expr literal grammar inside the string, so
//! one classifier owns numeric semantics for both syntaxes. Hand-rolled
//! parser (P1: no serde) with byte spans, depth cap, and finite-only
//! numbers.

use crate::ast::{Node, NodeKind, Span};
use crate::sexpr::{MAX_DEPTH, classify_atom};
use crate::{IrError, IrErrorKind};

fn err(span: Span, kind: IrErrorKind, detail: &str, hint: &str) -> IrError {
    IrError {
        span,
        kind,
        detail: detail.to_string(),
        hint: hint.to_string(),
    }
}

/// Parse one JSON-mapped program.
///
/// # Errors
/// Structured [`IrError`] with the offending span and a fix hint.
pub fn parse(src: &str) -> Result<Node, IrError> {
    let mut p = Parser {
        src,
        bytes: src.as_bytes(),
        pos: 0,
    };
    p.skip_ws();
    let node = p.parse_value(0)?;
    p.skip_ws();
    if p.pos != p.bytes.len() {
        return Err(err(
            Span::new(p.pos, p.bytes.len()),
            IrErrorKind::TrailingInput,
            "input continues after the top-level value",
            "a program is exactly one JSON value",
        ));
    }
    Ok(node)
}

struct Parser<'a> {
    src: &'a str,
    bytes: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn skip_ws(&mut self) {
        while let Some(&c) = self.bytes.get(self.pos) {
            if c == b' ' || c == b'\t' || c == b'\r' || c == b'\n' {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    fn expect(&mut self, c: u8, what: &str) -> Result<(), IrError> {
        if self.bytes.get(self.pos) == Some(&c) {
            self.pos += 1;
            Ok(())
        } else {
            Err(err(
                Span::new(self.pos, (self.pos + 1).min(self.bytes.len())),
                IrErrorKind::JsonSyntax,
                &format!("expected {what}"),
                "check the JSON structure against the fs-ir mapping table",
            ))
        }
    }

    fn parse_value(&mut self, depth: usize) -> Result<Node, IrError> {
        if depth > MAX_DEPTH {
            return Err(err(
                Span::new(self.pos, self.pos + 1),
                IrErrorKind::TooDeep,
                &format!("nesting exceeds the {MAX_DEPTH}-level cap"),
                "flatten the program; adversarial nesting is refused by design",
            ));
        }
        self.skip_ws();
        match self.bytes.get(self.pos) {
            Some(b'[') => self.parse_array(depth),
            Some(b'{') => self.parse_tagged(depth),
            Some(_) => Err(err(
                Span::new(self.pos, (self.pos + 1).min(self.bytes.len())),
                IrErrorKind::JsonSyntax,
                "expected an array (list) or a single-key tagged object (atom)",
                "atoms are {\"i\":..}/{\"f\":..}/{\"q\":..}/{\"c\":..}/{\"seed\":..}/\
                 {\"s\":..}/{\"sym\":..}/{\"kw\":..}",
            )),
            None => Err(err(
                Span::new(self.pos, self.pos),
                IrErrorKind::UnexpectedEnd,
                "expected a value, found end of input",
                "supply a complete JSON value",
            )),
        }
    }

    fn parse_array(&mut self, depth: usize) -> Result<Node, IrError> {
        let start = self.pos;
        self.pos += 1; // '['
        let mut items = Vec::new();
        self.skip_ws();
        if self.bytes.get(self.pos) == Some(&b']') {
            self.pos += 1;
            return Ok(Node {
                kind: NodeKind::List(items),
                span: Span::new(start, self.pos),
            });
        }
        loop {
            items.push(self.parse_value(depth + 1)?);
            self.skip_ws();
            match self.bytes.get(self.pos) {
                Some(b',') => {
                    self.pos += 1;
                }
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Node {
                        kind: NodeKind::List(items),
                        span: Span::new(start, self.pos),
                    });
                }
                _ => {
                    return Err(err(
                        Span::new(self.pos, (self.pos + 1).min(self.bytes.len())),
                        IrErrorKind::JsonSyntax,
                        "expected ',' or ']' in array",
                        "separate list items with commas",
                    ));
                }
            }
        }
    }

    fn parse_tagged(&mut self, _depth: usize) -> Result<Node, IrError> {
        let start = self.pos;
        self.pos += 1; // '{'
        self.skip_ws();
        let (tag, _tag_span) = self.parse_string_token()?;
        self.skip_ws();
        self.expect(b':', "':' after the tag")?;
        self.skip_ws();
        let kind = match tag.as_str() {
            "i" => {
                let (text, span) = self.parse_number_token()?;
                match text.parse::<i64>() {
                    Ok(v) => NodeKind::Int(v),
                    Err(_) => {
                        return Err(err(
                            span,
                            IrErrorKind::BadNumber,
                            &format!("{text:?} is not an i64"),
                            "\"i\" carries integers; use \"f\" for floats",
                        ));
                    }
                }
            }
            "f" => {
                let (text, span) = self.parse_number_token()?;
                match text.parse::<f64>() {
                    Ok(v) if v.is_finite() => NodeKind::Float(v),
                    _ => {
                        return Err(err(
                            span,
                            IrErrorKind::BadNumber,
                            &format!("{text:?} is not a finite f64"),
                            "IR literals must be finite",
                        ));
                    }
                }
            }
            "q" | "c" | "seed" => {
                let (text, span) = self.parse_string_token()?;
                let node = classify_atom(&text, span)?;
                match (tag.as_str(), &node.kind) {
                    ("q", NodeKind::Qty { .. })
                    | ("c", NodeKind::Count { .. })
                    | ("seed", NodeKind::Seed(_)) => node.kind,
                    _ => {
                        return Err(err(
                            span,
                            IrErrorKind::JsonTagMismatch,
                            &format!("literal {text:?} does not match tag {tag:?}"),
                            "\"q\" carries quantities like \"0.12Pa*s\"; \"c\" counts like \
                             \"384GiB\"; \"seed\" hex like \"0xF00D0002\"",
                        ));
                    }
                }
            }
            "s" => NodeKind::Str(self.parse_string_token()?.0),
            "sym" => NodeKind::Symbol(self.parse_string_token()?.0),
            "kw" => NodeKind::Keyword(self.parse_string_token()?.0),
            other => {
                return Err(err(
                    Span::new(start, self.pos),
                    IrErrorKind::JsonUnknownTag,
                    &format!("unknown atom tag {other:?}"),
                    "known tags: i f q c seed s sym kw",
                ));
            }
        };
        self.skip_ws();
        self.expect(b'}', "'}' closing the atom object")?;
        Ok(Node {
            kind,
            span: Span::new(start, self.pos),
        })
    }

    fn parse_string_token(&mut self) -> Result<(String, Span), IrError> {
        let start = self.pos;
        self.expect(b'"', "'\"' opening a string")?;
        let mut out = String::new();
        while let Some(&c) = self.bytes.get(self.pos) {
            match c {
                b'"' => {
                    self.pos += 1;
                    return Ok((out, Span::new(start, self.pos)));
                }
                b'\\' => {
                    let esc = self.bytes.get(self.pos + 1).copied();
                    match esc {
                        Some(b'"') => out.push('"'),
                        Some(b'\\') => out.push('\\'),
                        Some(b'/') => out.push('/'),
                        Some(b'n') => out.push('\n'),
                        Some(b't') => out.push('\t'),
                        Some(b'r') => out.push('\r'),
                        Some(b'b') => out.push('\u{0008}'),
                        Some(b'f') => out.push('\u{000C}'),
                        Some(b'u') => {
                            let hex_span =
                                Span::new(self.pos, (self.pos + 6).min(self.bytes.len()));
                            let hex = self
                                .src
                                .get(self.pos + 2..self.pos + 6)
                                .ok_or_else(|| self.u_escape_err(hex_span))?;
                            let code = u32::from_str_radix(hex, 16)
                                .map_err(|_| self.u_escape_err(hex_span))?;
                            let ch =
                                char::from_u32(code).ok_or_else(|| self.u_escape_err(hex_span))?;
                            out.push(ch);
                            self.pos += 4; // beyond the standard 2 below
                        }
                        _ => {
                            return Err(err(
                                Span::new(self.pos, (self.pos + 2).min(self.bytes.len())),
                                IrErrorKind::BadEscape,
                                "unknown escape sequence",
                                r#"supported: \ " / \n \t \r \b \f \uXXXX (no surrogate pairs)"#,
                            ));
                        }
                    }
                    self.pos += 2;
                }
                _ => {
                    let ch_len = self.src[self.pos..]
                        .chars()
                        .next()
                        .map_or(1, char::len_utf8);
                    out.push_str(&self.src[self.pos..self.pos + ch_len]);
                    self.pos += ch_len;
                }
            }
        }
        Err(err(
            Span::new(start, self.pos),
            IrErrorKind::UnclosedString,
            "string is never closed",
            "add the closing '\"'",
        ))
    }

    fn u_escape_err(&self, span: Span) -> IrError {
        let _ = self;
        err(
            span,
            IrErrorKind::BadEscape,
            "malformed \\uXXXX escape (surrogate halves are rejected)",
            "use 4 hex digits naming a Unicode scalar value",
        )
    }

    fn parse_number_token(&mut self) -> Result<(String, Span), IrError> {
        let start = self.pos;
        if self.bytes.get(self.pos) == Some(&b'-') {
            self.pos += 1;
        }
        while let Some(&c) = self.bytes.get(self.pos) {
            if c.is_ascii_digit() || matches!(c, b'.' | b'e' | b'E' | b'+' | b'-') {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return Err(err(
                Span::new(start, (start + 1).min(self.bytes.len())),
                IrErrorKind::JsonSyntax,
                "expected a number",
                "\"i\" and \"f\" tags carry JSON numbers",
            ));
        }
        Ok((
            self.src[start..self.pos].to_string(),
            Span::new(start, self.pos),
        ))
    }
}

/// Print a node in the canonical JSON mapping; `parse(print(x))` has the
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
            let _ = write!(out, "{{\"i\":{i}}}");
        }
        NodeKind::Float(f) => {
            let _ = write!(out, "{{\"f\":{f:?}}}");
        }
        NodeKind::Qty { text, .. } => {
            let _ = write!(out, "{{\"q\":");
            print_json_string(text, out);
            out.push('}');
        }
        NodeKind::Count { value, unit } => {
            let _ = write!(out, "{{\"c\":");
            let literal = match value {
                crate::ast::CountValue::Exact(v) => format!("{v}{}", unit.suffix()),
                crate::ast::CountValue::Fractional(decimal) => {
                    format!("{}{}", decimal.canonical(), unit.suffix())
                }
            };
            print_json_string(&literal, out);
            out.push('}');
        }
        NodeKind::Seed(v) => {
            let _ = write!(out, "{{\"seed\":\"0x{v:X}\"}}");
        }
        NodeKind::Str(s) => {
            out.push_str("{\"s\":");
            print_json_string(s, out);
            out.push('}');
        }
        NodeKind::Symbol(s) => {
            out.push_str("{\"sym\":");
            print_json_string(s, out);
            out.push('}');
        }
        NodeKind::Keyword(k) => {
            out.push_str("{\"kw\":");
            print_json_string(k, out);
            out.push('}');
        }
        NodeKind::List(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                print_into(item, out);
            }
            out.push(']');
        }
    }
}

fn print_json_string(s: &str, out: &mut String) {
    use std::fmt::Write as _;
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04X}", c as u32);
            }
            _ => out.push(c),
        }
    }
    out.push('"');
}
