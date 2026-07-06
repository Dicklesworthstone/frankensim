//! Minimal in-house JSON round-trip for [`QtyAny`]:
//! `{"value":0.12,"dims":[-1,1,-1,0,0]}` (dims order `[m,kg,s,K,A]`).
//!
//! In-house because the runtime dependency set is std + the Franken
//! constellation only (Decalogue P1) — serde is not on that list. The writer
//! uses Rust's shortest-round-trip float formatting, so `from_json(to_json(q))`
//! is bit-exact for finite values. Non-finite values are rejected at
//! serialization (JSON has no NaN/Infinity; ledger artifacts must not smuggle
//! them through text) — a documented policy, not an accident.

use crate::{Dims, QtyAny};
use core::fmt;

/// JSON encode/decode failures with position and guidance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JsonError {
    /// Byte offset (0 for serialization errors).
    pub at: usize,
    /// Description with fix guidance.
    pub message: String,
}

impl fmt::Display for JsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "QtyAny JSON error at byte {}: {}", self.at, self.message)
    }
}

impl core::error::Error for JsonError {}

/// Serialize to the canonical JSON object.
///
/// # Errors
/// Returns [`JsonError`] for non-finite values (JSON cannot represent them).
pub fn to_json(q: QtyAny) -> Result<String, JsonError> {
    if !q.value.is_finite() {
        return Err(JsonError {
            at: 0,
            message: format!(
                "non-finite value {:?} cannot be encoded as JSON; if this arrived from a \
                 computation, the computation should have reported a structured error instead",
                q.value
            ),
        });
    }
    let [m, kg, s, k, a] = q.dims.0;
    Ok(format!(
        "{{\"value\":{},\"dims\":[{m},{kg},{s},{k},{a}]}}",
        q.value
    ))
}

struct Cursor<'a> {
    s: &'a [u8],
    i: usize,
}

impl Cursor<'_> {
    fn skip_ws(&mut self) {
        while self.i < self.s.len() && matches!(self.s[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn expect(&mut self, tok: &str) -> Result<(), JsonError> {
        self.skip_ws();
        if self.s[self.i..].starts_with(tok.as_bytes()) {
            self.i += tok.len();
            Ok(())
        } else {
            Err(JsonError {
                at: self.i,
                message: format!("expected {tok:?}"),
            })
        }
    }

    fn number(&mut self) -> Result<f64, JsonError> {
        self.skip_ws();
        let start = self.i;
        while self.i < self.s.len()
            && matches!(
                self.s[self.i],
                b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E'
            )
        {
            self.i += 1;
        }
        core::str::from_utf8(&self.s[start..self.i])
            .ok()
            .and_then(|t| t.parse().ok())
            .ok_or(JsonError {
                at: start,
                message: "expected a JSON number".to_string(),
            })
    }

    fn int_i8(&mut self) -> Result<i8, JsonError> {
        let v = self.number()?;
        let i = v as i8;
        if (f64::from(i) - v).abs() > 0.0 {
            return Err(JsonError {
                at: self.i,
                message: format!("dimension exponent {v} is not a small integer"),
            });
        }
        Ok(i)
    }
}

/// Parse the canonical JSON object (field order fixed: value, dims).
///
/// # Errors
/// Returns [`JsonError`] with byte position on any deviation from the
/// canonical shape — this parser is intentionally strict: the writer is ours,
/// so any deviation indicates corruption, not dialect.
pub fn from_json(text: &str) -> Result<QtyAny, JsonError> {
    let mut c = Cursor {
        s: text.as_bytes(),
        i: 0,
    };
    c.expect("{")?;
    c.expect("\"value\"")?;
    c.expect(":")?;
    let value = c.number()?;
    c.expect(",")?;
    c.expect("\"dims\"")?;
    c.expect(":")?;
    c.expect("[")?;
    let m = c.int_i8()?;
    c.expect(",")?;
    let kg = c.int_i8()?;
    c.expect(",")?;
    let s = c.int_i8()?;
    c.expect(",")?;
    let k = c.int_i8()?;
    c.expect(",")?;
    let a = c.int_i8()?;
    c.expect("]")?;
    c.expect("}")?;
    c.skip_ws();
    if c.i != text.len() {
        return Err(JsonError {
            at: c.i,
            message: "trailing input after object".to_string(),
        });
    }
    Ok(QtyAny::new(value, Dims([m, kg, s, k, a])))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DynViscosity, Pressure};

    #[test]
    fn round_trip_is_bit_exact_for_finite_values() {
        // Deterministic grid including awkward values (subnormal-adjacent,
        // negative, high-precision decimals).
        let values = [
            0.0,
            -0.0,
            0.12,
            -3.5e-9,
            1.0 / 3.0,
            6.02214076e23,
            f64::MIN_POSITIVE,
            -f64::MAX / 2.0,
        ];
        for &v in &values {
            let q = QtyAny::new(v, Pressure::DIMS);
            let text = to_json(q).expect("finite");
            let back = from_json(&text).unwrap_or_else(|e| panic!("{text}: {e}"));
            assert_eq!(back.value.to_bits(), v.to_bits(), "value bits for {v}");
            assert_eq!(back.dims, Pressure::DIMS);
        }
    }

    #[test]
    fn canonical_shape_matches_spec() {
        let q = DynViscosity::new(0.12).erase();
        assert_eq!(
            to_json(q).unwrap(),
            r#"{"value":0.12,"dims":[-1,1,-1,0,0]}"#
        );
    }

    #[test]
    fn non_finite_values_are_refused_with_guidance() {
        let e = to_json(QtyAny::dimensionless(f64::NAN)).unwrap_err();
        assert!(e.message.contains("structured error"), "{e}");
        assert!(to_json(QtyAny::dimensionless(f64::INFINITY)).is_err());
    }

    #[test]
    fn corruption_is_rejected_with_position() {
        for bad in [
            "",
            "{}",
            r#"{"value":1}"#,
            r#"{"value":1,"dims":[1,2,3,4]}"#,
            r#"{"value":1,"dims":[1,2,3,4,5]} extra"#,
            r#"{"value":1,"dims":[1.5,0,0,0,0]}"#,
        ] {
            assert!(from_json(bad).is_err(), "should reject {bad:?}");
        }
    }

    #[test]
    fn whitespace_tolerant_parse() {
        let q =
            from_json(" { \"value\" : 2.5 , \"dims\" : [ 1 , 0 , -1 , 0 , 0 ] } ").expect("parses");
        assert!((q.value - 2.5).abs() < 1e-15);
        assert_eq!(q.dims, Dims([1, 0, -1, 0, 0]));
    }
}
