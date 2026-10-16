//! Minimal in-house JSON round-trip for [`QtyAny`]. The current canonical
//! wire is version 2:
//! `{"schema_version":2,"value":0.12,"dims":[-1,1,-1,0,0,0]}` with
//! dimensions ordered `[m,kg,s,K,A,mol]`.
//!
//! In-house because the runtime dependency set is std + the Franken
//! constellation only (Decalogue P1) — serde is not on that list. The writer
//! uses Rust's shortest-round-trip float formatting, so `from_json(to_json(q))`
//! is bit-exact for finite values. Non-finite values are rejected at
//! serialization (JSON has no NaN/Infinity; ledger artifacts must not smuggle
//! them through text) — a documented policy, not an accident. Exact legacy
//! version-1 five-vector bytes remain decodable only through [`decode_json`],
//! which returns an immutable old-hash → new-hash migration receipt.

use crate::{Dims, QtyAny};
use core::fmt;
use fs_blake3::hash_bytes;

pub use fs_blake3::ContentHash;

/// Historical implicit five-base wire version.
pub const LEGACY_WIRE_VERSION: u32 = 1;
/// Current explicit six-base wire version.
pub const WIRE_VERSION: u32 = 2;

/// Quantity JSON schema understood by the decoder.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum QtyWireVersion {
    /// Exact historical `{"value":...,"dims":[m,kg,s,K,A]}` bytes.
    LegacyFive = LEGACY_WIRE_VERSION,
    /// Current `{"schema_version":2,...,"dims":[m,kg,s,K,A,mol]}` bytes.
    SixBase = WIRE_VERSION,
}

/// The only admitted semantic rule for a legacy five-base quantity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FiveToSixRule {
    /// Preserve all five exponents and append an exact zero mole exponent.
    AppendMoleZero,
}

/// Immutable evidence that exact legacy bytes were mapped to canonical v2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DimensionCrosswalkReceipt {
    source_version: QtyWireVersion,
    target_version: QtyWireVersion,
    old_hash: ContentHash,
    new_hash: ContentHash,
    rule: FiveToSixRule,
}

impl DimensionCrosswalkReceipt {
    /// Source schema named by the receipt.
    #[must_use]
    pub const fn source_version(&self) -> QtyWireVersion {
        self.source_version
    }

    /// Target schema named by the receipt.
    #[must_use]
    pub const fn target_version(&self) -> QtyWireVersion {
        self.target_version
    }

    /// BLAKE3 content hash of the exact source bytes.
    #[must_use]
    pub const fn old_hash(&self) -> ContentHash {
        self.old_hash
    }

    /// BLAKE3 content hash of the exact canonical target bytes.
    #[must_use]
    pub const fn new_hash(&self) -> ContentHash {
        self.new_hash
    }

    /// Semantic migration rule applied.
    #[must_use]
    pub const fn rule(&self) -> FiveToSixRule {
        self.rule
    }

    /// Verify this receipt against the exact preserved source and target bytes.
    #[must_use]
    pub fn verifies(&self, old_bytes: &[u8], new_bytes: &[u8]) -> bool {
        self.source_version == QtyWireVersion::LegacyFive
            && self.target_version == QtyWireVersion::SixBase
            && self.rule == FiveToSixRule::AppendMoleZero
            && hash_bytes(old_bytes) == self.old_hash
            && hash_bytes(new_bytes) == self.new_hash
    }
}

/// Version-aware decode outcome. Legacy outcomes always carry a receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedQty {
    qty: QtyAny,
    source_version: QtyWireVersion,
    migration: Option<DimensionCrosswalkReceipt>,
}

impl DecodedQty {
    /// Decoded six-base quantity.
    #[must_use]
    pub const fn qty(&self) -> QtyAny {
        self.qty
    }

    /// Schema of the supplied bytes.
    #[must_use]
    pub const fn source_version(&self) -> QtyWireVersion {
        self.source_version
    }

    /// Migration evidence; present exactly for legacy five-vector input.
    #[must_use]
    pub const fn migration(&self) -> Option<&DimensionCrosswalkReceipt> {
        self.migration.as_ref()
    }
}

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
    let [m, kg, s, k, a, mol] = q.dims.0;
    Ok(format!(
        "{{\"schema_version\":{WIRE_VERSION},\"value\":{},\"dims\":[{m},{kg},{s},{k},{a},{mol}]}}",
        q.value
    ))
}

/// Reproduce the exact historical v1 bytes for a mol-free value.
///
/// This exists for immutable artifact verification and golden fixtures, not
/// as the default writer. New artifacts must use [`to_json`].
///
/// # Errors
/// Returns [`JsonError`] for non-finite values or a nonzero mole exponent.
pub fn to_legacy_json(q: QtyAny) -> Result<String, JsonError> {
    legacy_json(q, false)
}

fn legacy_json(q: QtyAny, explicit_version: bool) -> Result<String, JsonError> {
    if q.dims.0[5] != 0 {
        return Err(JsonError {
            at: 0,
            message: "a nonzero mole exponent cannot be represented by legacy five-base JSON"
                .to_string(),
        });
    }
    if !q.value.is_finite() {
        return Err(JsonError {
            at: 0,
            message: "non-finite values cannot be encoded as legacy JSON".to_string(),
        });
    }
    let [m, kg, s, k, a, _] = q.dims.0;
    if explicit_version {
        Ok(format!(
            "{{\"schema_version\":{LEGACY_WIRE_VERSION},\"value\":{},\"dims\":[{m},{kg},{s},{k},{a}]}}",
            q.value
        ))
    } else {
        Ok(format!(
            "{{\"value\":{},\"dims\":[{m},{kg},{s},{k},{a}]}}",
            q.value
        ))
    }
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

    /// Scan exactly one RFC 8259 JSON number and report whether its spelling
    /// is an integer (no fraction or exponent). The semantic integer readers
    /// use the spelling, not the resulting `f64`, so `1.0` and `1e0` cannot
    /// masquerade as dimension exponents or schema versions.
    fn number_span(&mut self) -> Result<(usize, usize, bool), JsonError> {
        self.skip_ws();
        let start = self.i;

        if self.s.get(self.i) == Some(&b'-') {
            self.i += 1;
        }

        match self.s.get(self.i).copied() {
            Some(b'0') => {
                self.i += 1;
                if self.s.get(self.i).is_some_and(u8::is_ascii_digit) {
                    return Err(JsonError {
                        at: self.i,
                        message: "JSON numbers cannot contain leading zeros".to_string(),
                    });
                }
            }
            Some(b'1'..=b'9') => {
                self.i += 1;
                while self.s.get(self.i).is_some_and(u8::is_ascii_digit) {
                    self.i += 1;
                }
            }
            _ => {
                return Err(JsonError {
                    at: start,
                    message: "expected a JSON number".to_string(),
                });
            }
        }

        let mut integer_syntax = true;
        if self.s.get(self.i) == Some(&b'.') {
            integer_syntax = false;
            self.i += 1;
            let fraction_start = self.i;
            while self.s.get(self.i).is_some_and(u8::is_ascii_digit) {
                self.i += 1;
            }
            if self.i == fraction_start {
                return Err(JsonError {
                    at: self.i,
                    message: "a JSON fraction requires at least one digit after '.'".to_string(),
                });
            }
        }

        if matches!(self.s.get(self.i).copied(), Some(b'e' | b'E')) {
            integer_syntax = false;
            self.i += 1;
            if matches!(self.s.get(self.i).copied(), Some(b'+' | b'-')) {
                self.i += 1;
            }
            let exponent_start = self.i;
            while self.s.get(self.i).is_some_and(u8::is_ascii_digit) {
                self.i += 1;
            }
            if self.i == exponent_start {
                return Err(JsonError {
                    at: self.i,
                    message: "a JSON exponent requires at least one digit".to_string(),
                });
            }
        }

        Ok((start, self.i, integer_syntax))
    }

    fn number(&mut self) -> Result<f64, JsonError> {
        let (start, end, _) = self.number_span()?;
        let value: f64 = core::str::from_utf8(&self.s[start..end])
            .ok()
            .and_then(|t| t.parse().ok())
            .ok_or(JsonError {
                at: start,
                message: "expected a JSON number".to_string(),
            })?;
        if !value.is_finite() {
            return Err(JsonError {
                at: start,
                message: "JSON number is outside the finite f64 domain".to_string(),
            });
        }
        Ok(value)
    }

    fn int_i8(&mut self) -> Result<i8, JsonError> {
        let (start, end, integer_syntax) = self.number_span()?;
        let raw = core::str::from_utf8(&self.s[start..end]).map_err(|_| JsonError {
            at: start,
            message: "dimension exponent is not valid UTF-8".to_string(),
        })?;
        if !integer_syntax {
            return Err(JsonError {
                at: start,
                message: format!(
                    "dimension exponent {raw:?} must use integer JSON syntax without a fraction or exponent"
                ),
            });
        }
        raw.parse::<i8>().map_err(|_| JsonError {
            at: start,
            message: format!("dimension exponent {raw:?} is not a small integer"),
        })
    }

    fn int_u32(&mut self) -> Result<u32, JsonError> {
        let (start, end, integer_syntax) = self.number_span()?;
        let raw = core::str::from_utf8(&self.s[start..end]).map_err(|_| JsonError {
            at: start,
            message: "schema version is not valid UTF-8".to_string(),
        })?;
        if !integer_syntax || raw.starts_with('-') {
            return Err(JsonError {
                at: start,
                message: format!(
                    "schema version {raw:?} must use unsigned integer JSON syntax without a sign, fraction, or exponent"
                ),
            });
        }
        raw.parse::<u32>().map_err(|_| JsonError {
            at: start,
            message: format!("schema version {raw:?} is not an unsigned integer"),
        })
    }
}

fn validate_canonical_source(
    text: &str,
    qty: QtyAny,
    source_version: QtyWireVersion,
    explicit_version: bool,
) -> Result<(), JsonError> {
    let canonical_source = match (source_version, explicit_version) {
        (QtyWireVersion::LegacyFive, false) => to_legacy_json(qty),
        (QtyWireVersion::LegacyFive, true) => legacy_json(qty, true),
        (QtyWireVersion::SixBase, true) => to_json(qty),
        (QtyWireVersion::SixBase, false) => Err(JsonError {
            at: 0,
            message: "six-base quantity JSON requires an explicit schema_version".to_string(),
        }),
    }?;
    if text == canonical_source {
        return Ok(());
    }
    let at = text
        .bytes()
        .zip(canonical_source.bytes())
        .position(|(actual, canonical)| actual != canonical)
        .unwrap_or(text.len().min(canonical_source.len()));
    Err(JsonError {
        at,
        message: "quantity JSON must exactly match its canonical versioned encoding; preserve historical v1 bytes and use to_json for v2"
            .to_string(),
    })
}

/// Decode either exact legacy v1 five-vector JSON or current v2 six-vector
/// JSON. The field order is fixed within each schema.
///
/// Legacy bytes are never silently reinterpreted: their outcome always
/// carries a [`DimensionCrosswalkReceipt`] that binds the exact old bytes to
/// the exact canonical v2 bytes.
///
/// # Errors
/// Returns [`JsonError`] with byte position on any deviation from the
/// canonical shape — this parser is intentionally strict: the writer is ours,
/// so any deviation indicates corruption, not dialect.
pub fn decode_json(text: &str) -> Result<DecodedQty, JsonError> {
    let mut c = Cursor {
        s: text.as_bytes(),
        i: 0,
    };
    c.expect("{")?;

    c.skip_ws();
    let explicit_version = c.s[c.i..].starts_with(b"\"schema_version\"");
    let source_version = if explicit_version {
        c.expect("\"schema_version\"")?;
        c.expect(":")?;
        let raw_version = c.int_u32()?;
        let version = match raw_version {
            LEGACY_WIRE_VERSION => QtyWireVersion::LegacyFive,
            WIRE_VERSION => QtyWireVersion::SixBase,
            _ => {
                return Err(JsonError {
                    at: c.i,
                    message: format!(
                        "unsupported quantity schema version {raw_version}; supported versions are {LEGACY_WIRE_VERSION} and {WIRE_VERSION}"
                    ),
                });
            }
        };
        c.expect(",")?;
        version
    } else {
        // Exact historical bytes carried no explicit tag; their fixed shape
        // is the immutable implicit-v1 schema.
        QtyWireVersion::LegacyFive
    };

    c.expect("\"value\"")?;
    c.expect(":")?;
    let value = c.number()?;
    if !value.is_finite() {
        return Err(JsonError {
            at: c.i,
            message: "non-finite quantity values are not valid ledger JSON".to_string(),
        });
    }
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
    let mol = if source_version == QtyWireVersion::SixBase {
        c.expect(",")?;
        c.int_i8()?
    } else {
        0
    };
    c.expect("]")?;
    c.expect("}")?;
    c.skip_ws();
    if c.i != text.len() {
        return Err(JsonError {
            at: c.i,
            message: "trailing input after object".to_string(),
        });
    }
    let qty = QtyAny::new(value, Dims([m, kg, s, k, a, mol]));
    validate_canonical_source(text, qty, source_version, explicit_version)?;
    let migration = if source_version == QtyWireVersion::LegacyFive {
        let new_bytes = to_json(qty)?;
        Some(DimensionCrosswalkReceipt {
            source_version,
            target_version: QtyWireVersion::SixBase,
            old_hash: hash_bytes(text.as_bytes()),
            new_hash: hash_bytes(new_bytes.as_bytes()),
            rule: FiveToSixRule::AppendMoleZero,
        })
    } else {
        None
    };
    Ok(DecodedQty {
        qty,
        source_version,
        migration,
    })
}

/// Parse only the current canonical v2 shape.
///
/// Legacy input is rejected here because returning only [`QtyAny`] would
/// discard mandatory migration evidence. Use [`decode_json`] for v1 bytes.
///
/// # Errors
/// Returns [`JsonError`] for malformed/unsupported input or for legacy input
/// whose receipt would otherwise be lost.
pub fn from_json(text: &str) -> Result<QtyAny, JsonError> {
    let decoded = decode_json(text)?;
    if decoded.migration.is_some() {
        return Err(JsonError {
            at: 0,
            message: "legacy five-base JSON requires decode_json so its semantic-crosswalk receipt is retained"
                .to_string(),
        });
    }
    Ok(decoded.qty)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AmountConcentration, DynViscosity, Pressure};

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
            r#"{"schema_version":2,"value":0.12,"dims":[-1,1,-1,0,0,0]}"#
        );
    }

    #[test]
    fn nonzero_mole_dimension_round_trips_in_v2() {
        let q = AmountConcentration::new(3.25).erase();
        let text = to_json(q).expect("finite");
        let decoded = decode_json(&text).expect("v2 decodes");
        assert_eq!(decoded.source_version(), QtyWireVersion::SixBase);
        assert!(decoded.migration().is_none());
        assert_eq!(decoded.qty(), q);
    }

    #[test]
    fn legacy_bytes_require_and_verify_immutable_crosswalk_receipt() {
        const OLD: &str = r#"{"value":0.12,"dims":[-1,1,-1,0,0]}"#;
        const NEW: &str = r#"{"schema_version":2,"value":0.12,"dims":[-1,1,-1,0,0,0]}"#;
        let decoded = decode_json(OLD).expect("legacy decodes with evidence");
        assert_eq!(decoded.source_version(), QtyWireVersion::LegacyFive);
        assert_eq!(decoded.qty().dims, Dims([-1, 1, -1, 0, 0, 0]));
        let receipt = decoded.migration().expect("receipt is mandatory");
        assert_eq!(receipt.rule(), FiveToSixRule::AppendMoleZero);
        assert_eq!(
            receipt.old_hash(),
            ContentHash::from_hex(
                "b97ca96f12cf487bc90760adad7257311fed950f95ab834c9107e51bf5f31ef1"
            )
            .expect("pinned old hash")
        );
        assert_eq!(
            receipt.new_hash(),
            ContentHash::from_hex(
                "8353a2a85f0de4a46f8cb31cb1673198c9bae9526b848369be545031d495bbb5"
            )
            .expect("pinned new hash")
        );
        assert!(receipt.verifies(OLD.as_bytes(), NEW.as_bytes()));
        assert!(!receipt.verifies(b"tampered", NEW.as_bytes()));
        assert!(!receipt.verifies(OLD.as_bytes(), b"tampered"));
        assert_eq!(to_legacy_json(decoded.qty()).unwrap(), OLD);
        assert!(from_json(OLD).unwrap_err().message.contains("receipt"));
    }

    #[test]
    fn explicit_v1_is_supported_but_version_arity_mismatches_fail_closed() {
        let v1 = r#"{"schema_version":1,"value":2.5,"dims":[1,0,-1,0,0]}"#;
        let decoded = decode_json(v1).expect("explicit v1");
        assert_eq!(decoded.qty().dims, Dims([1, 0, -1, 0, 0, 0]));
        assert!(decoded.migration().is_some());

        for bad in [
            r#"{"schema_version":1,"value":1,"dims":[1,2,3,4,5,6]}"#,
            r#"{"schema_version":2,"value":1,"dims":[1,2,3,4,5]}"#,
            r#"{"schema_version":3,"value":1,"dims":[1,2,3,4,5,6]}"#,
            r#"{"schema_version":+1,"value":1,"dims":[1,2,3,4,5]}"#,
            r#"{"schema_version":01,"value":1,"dims":[1,2,3,4,5]}"#,
            r#"{"schema_version":1.0,"value":1,"dims":[1,2,3,4,5]}"#,
            r#"{"schema_version":1e0,"value":1,"dims":[1,2,3,4,5]}"#,
            r#"{"schema_version":-1,"value":1,"dims":[1,2,3,4,5]}"#,
            r#"{"schema_version":4294967296,"value":1,"dims":[1,2,3,4,5]}"#,
        ] {
            assert!(decode_json(bad).is_err(), "must reject {bad}");
        }
    }

    #[test]
    fn json_number_grammar_is_strict_and_role_aware() {
        for expected in [0.0_f64, -0.0, 0.125, -3.5e-9, 6.022_140_76e23] {
            let text = to_json(QtyAny::dimensionless(expected)).expect("canonical number");
            let decoded = decode_json(&text).unwrap_or_else(|e| panic!("must accept {text}: {e}"));
            assert_eq!(decoded.qty().value.to_bits(), expected.to_bits(), "{text}");
        }

        for raw in ["+1", "01", "-01", ".5", "1.", "1e", "1e+", "1e999"] {
            let text = format!(r#"{{"schema_version":2,"value":{raw},"dims":[0,0,0,0,0,0]}}"#);
            assert!(decode_json(&text).is_err(), "must reject value {raw}");
        }

        for raw in ["+1", "01", "-01", "1.0", "1e0", "-0.0", "128", "-129"] {
            let text = format!(r#"{{"schema_version":2,"value":1,"dims":[{raw},0,0,0,0,0]}}"#);
            assert!(
                decode_json(&text).is_err(),
                "must reject dimension exponent {raw}"
            );
        }
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
            assert!(decode_json(bad).is_err(), "should reject {bad:?}");
        }
    }

    #[test]
    fn noncanonical_v1_and_v2_mutations_fail_closed() {
        for mutated in [
            " { \"schema_version\" : 2 , \"value\" : 2.5 , \"dims\" : [ 1 , 0 , -1 , 0 , 0 , 0 ] } ",
            r#"{"schema_version":2,"value":1e0,"dims":[0,0,0,0,0,0]}"#,
            r#"{"value":1,"schema_version":2,"dims":[0,0,0,0,0,0]}"#,
            r#"{ "value":1,"dims":[0,0,0,0,0]}"#,
            r#"{"schema_version":1,"value":1.0,"dims":[0,0,0,0,0]}"#,
            r#"{"schema_version":1,"dims":[0,0,0,0,0],"value":1}"#,
        ] {
            assert!(
                decode_json(mutated).is_err(),
                "noncanonical mutation must refuse before issuing a receipt: {mutated}"
            );
        }
    }
}
