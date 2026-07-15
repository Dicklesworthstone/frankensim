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

use crate::ast::{Node, NodeKind, Span, canonical_quantity_text};
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
    node.validate()?;
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
                            let code = u16::from_str_radix(hex, 16)
                                .map_err(|_| self.u_escape_err(hex_span))?;
                            if (0xD800..=0xDBFF).contains(&code) {
                                let pair_span =
                                    Span::new(self.pos, (self.pos + 12).min(self.bytes.len()));
                                if self.bytes.get(self.pos + 6) != Some(&b'\\')
                                    || self.bytes.get(self.pos + 7) != Some(&b'u')
                                {
                                    return Err(self.u_escape_err(pair_span));
                                }
                                let low_hex = self
                                    .src
                                    .get(self.pos + 8..self.pos + 12)
                                    .ok_or_else(|| self.u_escape_err(pair_span))?;
                                let low = u16::from_str_radix(low_hex, 16)
                                    .map_err(|_| self.u_escape_err(pair_span))?;
                                if !(0xDC00..=0xDFFF).contains(&low) {
                                    return Err(self.u_escape_err(pair_span));
                                }
                                let scalar = 0x1_0000
                                    + ((u32::from(code) - 0xD800) << 10)
                                    + (u32::from(low) - 0xDC00);
                                out.push(
                                    char::from_u32(scalar)
                                        .ok_or_else(|| self.u_escape_err(pair_span))?,
                                );
                                self.pos += 10; // plus the common two-byte escape advance below
                            } else if (0xDC00..=0xDFFF).contains(&code) {
                                return Err(self.u_escape_err(hex_span));
                            } else {
                                out.push(
                                    char::from_u32(u32::from(code))
                                        .ok_or_else(|| self.u_escape_err(hex_span))?,
                                );
                                self.pos += 4; // plus the common two-byte escape advance below
                            }
                        }
                        _ => {
                            return Err(err(
                                Span::new(self.pos, (self.pos + 2).min(self.bytes.len())),
                                IrErrorKind::BadEscape,
                                "unknown escape sequence",
                                r#"supported: \ " / \n \t \r \b \f \uXXXX (paired when surrogate-encoded)"#,
                            ));
                        }
                    }
                    self.pos += 2;
                }
                0x00..=0x1F => {
                    return Err(err(
                        Span::new(self.pos, self.pos + 1),
                        IrErrorKind::JsonSyntax,
                        "unescaped control character in JSON string",
                        "encode control characters with a JSON escape such as \\n or \\u000A",
                    ));
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
            "malformed \\uXXXX escape or surrogate pair",
            "use one Unicode scalar or a high-surrogate/low-surrogate pair",
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
        let text = &self.src[start..self.pos];
        if !is_rfc8259_number(text) {
            return Err(err(
                Span::new(start, self.pos),
                IrErrorKind::BadNumber,
                &format!("{text:?} is not an RFC 8259 JSON number"),
                "use -?(0|[1-9][0-9]*)(.[0-9]+)?([eE][+-]?[0-9]+)?",
            ));
        }
        Ok((text.to_string(), Span::new(start, self.pos)))
    }
}

fn is_rfc8259_number(text: &str) -> bool {
    let bytes = text.as_bytes();
    let mut index = 0;
    if bytes.get(index) == Some(&b'-') {
        index += 1;
    }
    match bytes.get(index) {
        Some(b'0') => index += 1,
        Some(b'1'..=b'9') => {
            index += 1;
            while bytes.get(index).is_some_and(u8::is_ascii_digit) {
                index += 1;
            }
        }
        _ => return false,
    }
    if bytes.get(index) == Some(&b'.') {
        index += 1;
        let start = index;
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        if index == start {
            return false;
        }
    }
    if matches!(bytes.get(index), Some(b'e' | b'E')) {
        index += 1;
        if matches!(bytes.get(index), Some(b'+' | b'-')) {
            index += 1;
        }
        let start = index;
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        if index == start {
            return false;
        }
    }
    index == bytes.len()
}

/// Validate and print a node in the canonical JSON mapping.
///
/// # Errors
/// Rejects the first invalid atom with its source span and exact tree path.
pub fn print(node: &Node) -> Result<String, IrError> {
    print_checked(node)
}

/// Validate and print one canonical JSON mapping.
///
/// # Errors
/// Rejects the first invalid atom with its source span and exact tree path.
pub fn print_checked(node: &Node) -> Result<String, IrError> {
    node.validate()?;
    let mut out = String::new();
    print_into(node, &mut out)?;
    Ok(out)
}

fn print_into(node: &Node, out: &mut String) -> Result<(), IrError> {
    use std::fmt::Write as _;
    match &node.kind {
        NodeKind::Int(i) => {
            let _ = write!(out, "{{\"i\":{i}}}");
        }
        NodeKind::Float(f) => {
            let _ = write!(out, "{{\"f\":{f:?}}}");
        }
        NodeKind::Qty { value, dims, .. } => {
            let _ = write!(out, "{{\"q\":");
            print_json_string(&canonical_quantity_text(*value, *dims, node.span)?, out);
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
                print_into(item, out)?;
            }
            out.push(']');
        }
    }
    Ok(())
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
