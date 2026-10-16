//! Minimal s-expression reader for the skeleton's study format — the embryo
//! of fs-ir's canonical syntax (plan §11.1). Hand-rolled, total (never
//! panics), position-carrying errors.

/// An s-expression: an atom or a list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Sexpr {
    /// Bare token (symbol, number literal, quoted string with quotes removed).
    Atom(String),
    /// Parenthesized list.
    List(Vec<Sexpr>),
}

impl Sexpr {
    /// If this is a list whose head is the atom `name`, return the tail.
    #[must_use]
    pub fn as_form(&self, name: &str) -> Option<&[Sexpr]> {
        match self {
            Sexpr::List(items) => match items.split_first() {
                Some((Sexpr::Atom(head), tail)) if head == name => Some(tail),
                _ => None,
            },
            Sexpr::Atom(_) => None,
        }
    }

    /// Atom text, if an atom.
    #[must_use]
    pub fn atom(&self) -> Option<&str> {
        match self {
            Sexpr::Atom(a) => Some(a),
            Sexpr::List(_) => None,
        }
    }

    /// Find the first sub-form `(name ...)` in a list body.
    #[must_use]
    pub fn find_form<'a>(body: &'a [Sexpr], name: &str) -> Option<&'a [Sexpr]> {
        body.iter().find_map(|s| s.as_form(name))
    }
}

/// Parse one top-level s-expression.
///
/// # Errors
/// Returns `(byte position, message)` on malformed input.
pub fn parse(input: &str) -> Result<Sexpr, (usize, String)> {
    let mut pos = 0usize;
    let bytes = input.as_bytes();
    let expr = parse_at(input, bytes, &mut pos)?;
    skip_ws(bytes, &mut pos);
    if pos != bytes.len() {
        return Err((pos, "trailing input after the study form".to_string()));
    }
    Ok(expr)
}

fn skip_ws(bytes: &[u8], pos: &mut usize) {
    while *pos < bytes.len() {
        match bytes[*pos] {
            b' ' | b'\t' | b'\n' | b'\r' => *pos += 1,
            b';' => {
                while *pos < bytes.len() && bytes[*pos] != b'\n' {
                    *pos += 1;
                }
            }
            _ => break,
        }
    }
}

fn parse_at(input: &str, bytes: &[u8], pos: &mut usize) -> Result<Sexpr, (usize, String)> {
    skip_ws(bytes, pos);
    if *pos >= bytes.len() {
        return Err((*pos, "unexpected end of input".to_string()));
    }
    match bytes[*pos] {
        b'(' => {
            *pos += 1;
            let mut items = Vec::new();
            loop {
                skip_ws(bytes, pos);
                if *pos >= bytes.len() {
                    return Err((*pos, "unclosed ( — missing )".to_string()));
                }
                if bytes[*pos] == b')' {
                    *pos += 1;
                    return Ok(Sexpr::List(items));
                }
                items.push(parse_at(input, bytes, pos)?);
            }
        }
        b')' => Err((*pos, "unexpected )".to_string())),
        b'"' => {
            let start = *pos + 1;
            let mut end = start;
            while end < bytes.len() && bytes[end] != b'"' {
                end += 1;
            }
            if end >= bytes.len() {
                return Err((start, "unterminated string".to_string()));
            }
            *pos = end + 1;
            Ok(Sexpr::Atom(input[start..end].to_string()))
        }
        _ => {
            let start = *pos;
            while *pos < bytes.len()
                && !matches!(
                    bytes[*pos],
                    b'(' | b')' | b'"' | b' ' | b'\t' | b'\n' | b'\r' | b';'
                )
            {
                *pos += 1;
            }
            Ok(Sexpr::Atom(input[start..*pos].to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nested_forms_with_comments_and_strings() {
        let s = parse("(study \"pv-1\" ; comment\n (seed 0x5EED) (grid 33))").expect("parses");
        let body = s.as_form("study").expect("study form");
        assert_eq!(body[0].atom(), Some("pv-1"));
        assert_eq!(
            Sexpr::find_form(body, "grid").and_then(|g| g[0].atom()),
            Some("33")
        );
    }

    #[test]
    fn errors_carry_positions_and_never_panic() {
        for bad in ["", "(", ")", "(study \"x", "(a) b"] {
            let e = parse(bad).expect_err(bad);
            assert!(!e.1.is_empty());
        }
        // Garbage battery: totality.
        let mut st: u64 = 7;
        for _ in 0..5_000 {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            let s: String = (0..st % 20)
                .map(|i| b"()\"; ax0\n"[((st >> (i % 50)) % 9) as usize] as char)
                .collect();
            let _ = parse(&s);
        }
    }
}
