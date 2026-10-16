//! Minimal strict JSON for the browser boundary: one FLAT object whose values
//! are finite numbers, arrays of finite numbers, booleans or strings. Anything
//! else is refused rather than guessed. No serde: keeps the wasm small and the
//! admission rules explicit.

/// One admitted value of the flat parameter object.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// A finite JSON number.
    Number(f64),
    /// An array of finite JSON numbers.
    Numbers(Vec<f64>),
    /// `true` / `false`.
    Bool(bool),
    /// A string without escape processing beyond `\"` and `\\`.
    Text(String),
}

struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Cursor<'_> {
    fn skip(&mut self) {
        while self.at < self.bytes.len() && self.bytes[self.at].is_ascii_whitespace() {
            self.at += 1;
        }
    }
    fn peek(&mut self) -> Option<u8> {
        self.skip();
        self.bytes.get(self.at).copied()
    }
    fn expect(&mut self, byte: u8) -> Result<(), String> {
        if self.peek() == Some(byte) {
            self.at += 1;
            Ok(())
        } else {
            Err(format!("expected '{}' at byte {}", byte as char, self.at))
        }
    }
    fn string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut out = Vec::new();
        loop {
            let byte = *self.bytes.get(self.at).ok_or("unterminated string")?;
            self.at += 1;
            match byte {
                b'"' => break,
                b'\\' => {
                    let escaped = *self.bytes.get(self.at).ok_or("unterminated escape")?;
                    self.at += 1;
                    match escaped {
                        b'"' | b'\\' | b'/' => out.push(escaped),
                        _ => return Err("unsupported string escape".into()),
                    }
                }
                other => out.push(other),
            }
        }
        String::from_utf8(out).map_err(|_| "string is not UTF-8".to_string())
    }
    fn number(&mut self) -> Result<f64, String> {
        self.skip();
        let start = self.at;
        while self.at < self.bytes.len()
            && matches!(self.bytes[self.at], b'0'..=b'9' | b'+' | b'-' | b'.' | b'e' | b'E')
        {
            self.at += 1;
        }
        let text = core::str::from_utf8(&self.bytes[start..self.at]).map_err(|_| "bad number")?;
        let value: f64 = text
            .parse()
            .map_err(|_| format!("invalid number at byte {start}"))?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(format!("nonfinite number at byte {start}"))
        }
    }
    fn literal(&mut self, word: &str) -> bool {
        if self.bytes[self.at..].starts_with(word.as_bytes()) {
            self.at += word.len();
            true
        } else {
            false
        }
    }
    fn value(&mut self) -> Result<Value, String> {
        match self.peek().ok_or("missing value")? {
            b'"' => Ok(Value::Text(self.string()?)),
            b'[' => {
                self.at += 1;
                let mut items = Vec::new();
                if self.peek() == Some(b']') {
                    self.at += 1;
                    return Ok(Value::Numbers(items));
                }
                loop {
                    items.push(self.number()?);
                    if items.len() > 64 {
                        return Err("numeric array exceeds 64 entries".into());
                    }
                    match self.peek() {
                        Some(b',') => self.at += 1,
                        Some(b']') => {
                            self.at += 1;
                            return Ok(Value::Numbers(items));
                        }
                        _ => return Err("malformed numeric array".into()),
                    }
                }
            }
            b't' if self.literal("true") => Ok(Value::Bool(true)),
            b'f' if self.literal("false") => Ok(Value::Bool(false)),
            b'{' => Err("nested objects are not part of the flat parameter object".into()),
            b'n' => Err("null is not an admitted parameter value; omit the key instead".into()),
            _ => Ok(Value::Number(self.number()?)),
        }
    }
}

/// Parse one flat object. Duplicate keys are refused.
///
/// # Errors
/// Any syntax outside the admitted flat subset, with a byte position.
pub fn parse_flat_object(text: &str) -> Result<Vec<(String, Value)>, String> {
    if text.len() > 65_536 {
        return Err("parameter object exceeds 64 KiB".into());
    }
    let mut cursor = Cursor { bytes: text.as_bytes(), at: 0 };
    let mut out: Vec<(String, Value)> = Vec::new();
    cursor.expect(b'{')?;
    if cursor.peek() == Some(b'}') {
        cursor.at += 1;
    } else {
        loop {
            cursor.skip();
            let key = cursor.string()?;
            cursor.expect(b':')?;
            let value = cursor.value()?;
            if out.iter().any(|(k, _)| *k == key) {
                return Err(format!("duplicate key \"{key}\""));
            }
            out.push((key, value));
            match cursor.peek() {
                Some(b',') => cursor.at += 1,
                Some(b'}') => {
                    cursor.at += 1;
                    break;
                }
                _ => return Err(format!("expected ',' or '}}' at byte {}", cursor.at)),
            }
        }
    }
    if cursor.peek().is_some() {
        return Err("trailing characters after the parameter object".into());
    }
    Ok(out)
}

/// Shortest round-trip decimal for a finite f64 (`null` if nonfinite).
#[must_use]
pub fn number(value: f64) -> String {
    if value.is_finite() {
        let text = format!("{value:?}");
        // Rust prints 1e21 as "1e21" and 0.1 as "0.1": both valid JSON.
        text
    } else {
        "null".to_string()
    }
}

/// `[a,b,c]` for finite f64 values.
#[must_use]
pub fn numbers(values: &[f64]) -> String {
    let mut out = String::from("[");
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&number(*v));
    }
    out.push(']');
    out
}

/// A JSON string literal with the mandatory escapes.
#[must_use]
pub fn string(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for ch in text.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_objects_parse_and_everything_else_is_refused() {
        let parsed = parse_flat_object(
            r#" {"a": 1.5e3, "b":[1, -2.0,3], "c": true, "d":"x\"y"} "#,
        )
        .unwrap();
        assert_eq!(parsed[0], ("a".into(), Value::Number(1500.0)));
        assert_eq!(parsed[1], ("b".into(), Value::Numbers(vec![1.0, -2.0, 3.0])));
        assert_eq!(parsed[2], ("c".into(), Value::Bool(true)));
        assert_eq!(parsed[3], ("d".into(), Value::Text("x\"y".into())));
        assert!(parse_flat_object("{}").unwrap().is_empty());
        for bad in [
            "", "[]", "{\"a\":}", "{\"a\":1,}", "{\"a\":1,\"a\":2}", "{\"a\":{\"b\":1}}",
            "{\"a\":null}", "{\"a\":1} x", "{\"a\":1e999}", "{\"a\":[1,\"x\"]}", "{a:1}",
        ] {
            assert!(parse_flat_object(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn writers_emit_valid_json_fragments() {
        assert_eq!(number(0.1), "0.1");
        assert_eq!(number(f64::NAN), "null");
        assert_eq!(numbers(&[1.0, 2.5]), "[1.0,2.5]");
        assert_eq!(string("a\"b\\c\n"), "\"a\\\"b\\\\c\\n\"");
    }
}
