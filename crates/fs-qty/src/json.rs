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

    fn int_u32(&mut self) -> Result<u32, JsonError> {
        let at = self.i;
        let v = self.number()?;
        let i = v as u32;
        if !v.is_finite() || f64::from(i) != v {
            return Err(JsonError {
                at,
                message: format!("schema version {v} is not an unsigned integer"),
            });
        }
        Ok(i)
    }
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
    let source_version = if c.s[c.i..].starts_with(b"\"schema_version\"") {
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
        const NEW: &str =
            r#"{"schema_version":2,"value":0.12,"dims":[-1,1,-1,0,0,0]}"#;
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
        ] {
            assert!(decode_json(bad).is_err(), "must reject {bad}");
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
    fn whitespace_tolerant_parse() {
        let q = from_json(
            " { \"schema_version\" : 2 , \"value\" : 2.5 , \"dims\" : [ 1 , 0 , -1 , 0 , 0 , 0 ] } ",
        )
        .expect("parses");
        assert!((q.value - 2.5).abs() < 1e-15);
        assert_eq!(q.dims, Dims([1, 0, -1, 0, 0, 0]));
    }
}
