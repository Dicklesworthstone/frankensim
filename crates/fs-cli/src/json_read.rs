//! Strict, allocation-bounded JSON reader for the CLI's own retained receipts.
//!
//! The solve stages write their receipts with hand-rolled writers (no serde
//! in the Franken-only runtime graph). The report stage and the export verbs
//! must read those receipts back from the ledger without trusting them: a
//! receipt is an artifact, and an artifact may be truncated, forged, or from
//! an older writer. This reader therefore:
//!
//! - accepts exactly RFC 8259 JSON (objects, arrays, strings, numbers,
//!   `true`/`false`/`null`), with no comments, no trailing commas, no NaN;
//! - refuses inputs deeper than [`MAX_DEPTH`] or longer than [`MAX_BYTES`]
//!   so a hostile artifact cannot exhaust the stack or memory;
//! - keeps the original number spelling alongside the parsed `f64`, because
//!   receipt identities are computed over the exact bytes the writer emitted;
//! - preserves object key order (receipt writers emit canonical order, and
//!   the report renders in that order).
//!
//! It makes no semantic claim: a value that parses is merely well-formed.

use std::fmt;

/// Maximum accepted nesting depth (objects and arrays combined).
pub const MAX_DEPTH: usize = 64;
/// Maximum accepted input length in bytes.
pub const MAX_BYTES: usize = 64 * 1024 * 1024;

/// A parsed JSON value.
#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    /// `null`.
    Null,
    /// `true` / `false`.
    Bool(bool),
    /// A number with its exact source spelling.
    Number { value: f64, raw: String },
    /// A string with escapes resolved.
    Str(String),
    /// An array in source order.
    Array(Vec<JsonValue>),
    /// An object in source order (duplicate keys refused).
    Object(Vec<(String, JsonValue)>),
}

/// A parse refusal with a byte offset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JsonReadError {
    /// Byte offset where the problem was detected.
    pub offset: usize,
    /// What was wrong.
    pub what: String,
}

impl fmt::Display for JsonReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid JSON at byte {}: {}", self.offset, self.what)
    }
}

impl JsonValue {
    /// Parse one complete JSON document.
    ///
    /// # Errors
    /// [`JsonReadError`] on any syntax violation, on trailing bytes, on
    /// duplicate object keys, or when the depth or size bounds are exceeded.
    pub fn parse(text: &str) -> Result<JsonValue, JsonReadError> {
        if text.len() > MAX_BYTES {
            return Err(JsonReadError {
                offset: 0,
                what: format!("document is {} bytes, above the {MAX_BYTES}-byte bound", text.len()),
            });
        }
        let mut parser = Parser {
            bytes: text.as_bytes(),
            pos: 0,
        };
        parser.skip_ws();
        let value = parser.value(0)?;
        parser.skip_ws();
        if parser.pos != parser.bytes.len() {
            return Err(parser.error("trailing bytes after the document"));
        }
        Ok(value)
    }

    /// Object member lookup by exact key.
    #[must_use]
    pub fn get(&self, key: &str) -> Option<&JsonValue> {
        match self {
            JsonValue::Object(members) => members
                .iter()
                .find_map(|(k, v)| (k == key).then_some(v)),
            _ => None,
        }
    }

    /// Nested lookup along a key path.
    #[must_use]
    pub fn path(&self, keys: &[&str]) -> Option<&JsonValue> {
        keys.iter().try_fold(self, |node, key| node.get(key))
    }

    /// The string payload, if this is a string.
    #[must_use]
    pub fn as_str(&self) -> Option<&str> {
        match self {
            JsonValue::Str(s) => Some(s.as_str()),
            _ => None,
        }
    }

    /// The numeric value, if this is a number.
    #[must_use]
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            JsonValue::Number { value, .. } => Some(*value),
            _ => None,
        }
    }

    /// The exact number spelling, if this is a number.
    #[must_use]
    pub fn number_raw(&self) -> Option<&str> {
        match self {
            JsonValue::Number { raw, .. } => Some(raw.as_str()),
            _ => None,
        }
    }

    /// The array items, if this is an array.
    #[must_use]
    pub fn as_array(&self) -> Option<&[JsonValue]> {
        match self {
            JsonValue::Array(items) => Some(items.as_slice()),
            _ => None,
        }
    }

    /// The object members in source order, if this is an object.
    #[must_use]
    pub fn as_object(&self) -> Option<&[(String, JsonValue)]> {
        match self {
            JsonValue::Object(members) => Some(members.as_slice()),
            _ => None,
        }
    }

    /// String member by key, or `None` when absent or not a string.
    #[must_use]
    pub fn str_field(&self, key: &str) -> Option<&str> {
        self.get(key).and_then(JsonValue::as_str)
    }

    /// Numeric member by key, or `None` when absent or not a number.
    #[must_use]
    pub fn f64_field(&self, key: &str) -> Option<f64> {
        self.get(key).and_then(JsonValue::as_f64)
    }
}

struct Parser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn error(&self, what: impl Into<String>) -> JsonReadError {
        JsonReadError {
            offset: self.pos,
            what: what.into(),
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.pos += 1;
        }
    }

    fn expect_literal(&mut self, literal: &str) -> Result<(), JsonReadError> {
        if self.bytes[self.pos..].starts_with(literal.as_bytes()) {
            self.pos += literal.len();
            Ok(())
        } else {
            Err(self.error(format!("expected `{literal}`")))
        }
    }

    fn value(&mut self, depth: usize) -> Result<JsonValue, JsonReadError> {
        if depth > MAX_DEPTH {
            return Err(self.error(format!("nesting deeper than {MAX_DEPTH}")));
        }
        match self.peek() {
            None => Err(self.error("unexpected end of input")),
            Some(b'{') => self.object(depth),
            Some(b'[') => self.array(depth),
            Some(b'"') => self.string().map(JsonValue::Str),
            Some(b't') => self.expect_literal("true").map(|()| JsonValue::Bool(true)),
            Some(b'f') => self.expect_literal("false").map(|()| JsonValue::Bool(false)),
            Some(b'n') => self.expect_literal("null").map(|()| JsonValue::Null),
            Some(b'-' | b'0'..=b'9') => self.number(),
            Some(other) => Err(self.error(format!("unexpected byte 0x{other:02x}"))),
        }
    }

    fn object(&mut self, depth: usize) -> Result<JsonValue, JsonReadError> {
        self.pos += 1; // '{'
        let mut members: Vec<(String, JsonValue)> = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(JsonValue::Object(members));
        }
        loop {
            self.skip_ws();
            if self.peek() != Some(b'"') {
                return Err(self.error("expected a string key"));
            }
            let key = self.string()?;
            if members.iter().any(|(existing, _)| *existing == key) {
                return Err(self.error(format!("duplicate key `{key}`")));
            }
            self.skip_ws();
            if self.peek() != Some(b':') {
                return Err(self.error("expected `:` after key"));
            }
            self.pos += 1;
            self.skip_ws();
            let value = self.value(depth + 1)?;
            members.push((key, value));
            self.skip_ws();
            match self.peek() {
                Some(b',') => {
                    self.pos += 1;
                }
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(JsonValue::Object(members));
                }
                _ => return Err(self.error("expected `,` or `}` in object")),
            }
        }
    }

    fn array(&mut self, depth: usize) -> Result<JsonValue, JsonReadError> {
        self.pos += 1; // '['
        let mut items = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(JsonValue::Array(items));
        }
        loop {
            self.skip_ws();
            items.push(self.value(depth + 1)?);
            self.skip_ws();
            match self.peek() {
                Some(b',') => {
                    self.pos += 1;
                }
                Some(b']') => {
                    self.pos += 1;
                    return Ok(JsonValue::Array(items));
                }
                _ => return Err(self.error("expected `,` or `]` in array")),
            }
        }
    }

    fn number(&mut self) -> Result<JsonValue, JsonReadError> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        match self.peek() {
            Some(b'0') => {
                self.pos += 1;
            }
            Some(b'1'..=b'9') => {
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.pos += 1;
                }
            }
            _ => return Err(self.error("expected a digit")),
        }
        if self.peek() == Some(b'.') {
            self.pos += 1;
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err(self.error("expected a digit after `.`"));
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err(self.error("expected a digit in exponent"));
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        let raw = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| self.error("number bytes are not UTF-8"))?;
        let value: f64 = raw
            .parse()
            .map_err(|_| self.error(format!("number `{raw}` does not parse")))?;
        if !value.is_finite() {
            return Err(self.error(format!("number `{raw}` is not finite")));
        }
        Ok(JsonValue::Number {
            value,
            raw: raw.to_string(),
        })
    }

    fn string(&mut self) -> Result<String, JsonReadError> {
        self.pos += 1; // opening quote
        let mut out = String::new();
        loop {
            let Some(byte) = self.peek() else {
                return Err(self.error("unterminated string"));
            };
            match byte {
                b'"' => {
                    self.pos += 1;
                    return Ok(out);
                }
                b'\\' => {
                    self.pos += 1;
                    let Some(escape) = self.peek() else {
                        return Err(self.error("unterminated escape"));
                    };
                    self.pos += 1;
                    match escape {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{0008}'),
                        b'f' => out.push('\u{000C}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let unit = self.hex4()?;
                            let ch = if (0xD800..=0xDBFF).contains(&unit) {
                                // Surrogate pair.
                                if !self.bytes[self.pos..].starts_with(b"\\u") {
                                    return Err(self.error("unpaired high surrogate"));
                                }
                                self.pos += 2;
                                let low = self.hex4()?;
                                if !(0xDC00..=0xDFFF).contains(&low) {
                                    return Err(self.error("invalid low surrogate"));
                                }
                                let code =
                                    0x10000 + ((u32::from(unit) - 0xD800) << 10) + (u32::from(low) - 0xDC00);
                                char::from_u32(code).ok_or_else(|| self.error("invalid code point"))?
                            } else if (0xDC00..=0xDFFF).contains(&unit) {
                                return Err(self.error("unpaired low surrogate"));
                            } else {
                                char::from_u32(u32::from(unit))
                                    .ok_or_else(|| self.error("invalid code point"))?
                            };
                            out.push(ch);
                        }
                        other => {
                            return Err(self.error(format!("invalid escape `\\{}`", other as char)));
                        }
                    }
                }
                0x00..=0x1F => return Err(self.error("raw control character in string")),
                _ => {
                    // Copy one UTF-8 scalar; the input is a `&str`, so the
                    // sequence is valid by construction.
                    let rest = std::str::from_utf8(&self.bytes[self.pos..])
                        .map_err(|_| self.error("string bytes are not UTF-8"))?;
                    let ch = rest.chars().next().expect("non-empty remainder");
                    out.push(ch);
                    self.pos += ch.len_utf8();
                }
            }
        }
    }

    fn hex4(&mut self) -> Result<u16, JsonReadError> {
        let end = self.pos + 4;
        if end > self.bytes.len() {
            return Err(self.error("truncated \\u escape"));
        }
        let hex = std::str::from_utf8(&self.bytes[self.pos..end])
            .map_err(|_| self.error("\\u escape is not UTF-8"))?;
        let unit = u16::from_str_radix(hex, 16).map_err(|_| self.error("invalid \\u escape"))?;
        self.pos = end;
        Ok(unit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_receipt_shaped_documents_and_keeps_number_spelling() {
        let doc = r#"{"schema":"x.v1","qoi":[{"name":"temperature-max","value":312.5,"unit":"kelvin"}],"terms":[],"total":"unknown","ok":true,"none":null,"neg":-1e-3}"#;
        let value = JsonValue::parse(doc).expect("parses");
        assert_eq!(value.str_field("schema"), Some("x.v1"));
        let first = &value.get("qoi").unwrap().as_array().unwrap()[0];
        assert_eq!(first.f64_field("value"), Some(312.5));
        assert_eq!(first.get("value").unwrap().number_raw(), Some("312.5"));
        assert_eq!(value.get("neg").unwrap().number_raw(), Some("-1e-3"));
        assert_eq!(value.get("ok"), Some(&JsonValue::Bool(true)));
        assert_eq!(value.get("none"), Some(&JsonValue::Null));
        assert_eq!(value.path(&["qoi"]).unwrap().as_array().unwrap().len(), 1);
    }

    #[test]
    fn refuses_malformed_inputs() {
        for bad in [
            "",
            "{",
            "[1,]",
            "{\"a\":1,}",
            "{\"a\":1,\"a\":2}",
            "01",
            "1.",
            "NaN",
            "\"\\x\"",
            "{\"a\":1} x",
            "\"unterminated",
            "\"\\ud800\"",
        ] {
            assert!(JsonValue::parse(bad).is_err(), "accepted {bad:?}");
        }
    }

    #[test]
    fn bounds_depth() {
        let deep = "[".repeat(MAX_DEPTH + 2) + &"]".repeat(MAX_DEPTH + 2);
        assert!(JsonValue::parse(&deep).is_err());
        let ok = "[".repeat(MAX_DEPTH) + &"]".repeat(MAX_DEPTH);
        assert!(JsonValue::parse(&ok).is_ok());
    }

    #[test]
    fn resolves_escapes_and_surrogates() {
        let value = JsonValue::parse(r#""a\"b\\c\n\u00e9\ud83d\ude00""#).expect("parses");
        assert_eq!(value.as_str(), Some("a\"b\\c\né😀"));
    }
}
