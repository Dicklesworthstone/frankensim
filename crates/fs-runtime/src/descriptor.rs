//! Assembly-derived descriptors, structured refusals, and canonical JSON writers.

use core::fmt::{self, Write as _};
use fs_qty::Dims;

/// Schema version for exported [`SnapshotDescriptor`] JSON payloads.
pub const SNAPSHOT_DESCRIPTOR_SCHEMA: &str = "frankensim.runtime.snapshot_descriptor.v1";

/// Structured refusal surfaced by `AssemblyRuntime`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeRefusal {
    /// Stable machine-readable refusal code.
    pub code: &'static str,
    /// Lifecycle stage where refusal occurred (`load`, `input`, `step`, `probe`).
    pub stage: &'static str,
    /// Human-readable explanation.
    pub message: String,
}

impl RuntimeRefusal {
    /// Construct a typed refusal.
    #[must_use]
    pub fn new(code: &'static str, stage: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            stage,
            message: message.into(),
        }
    }

    /// Canonical JSON encoding `{"code":...,"stage":...,"message":...}`.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::with_capacity(64 + self.message.len());
        out.push_str("{\"code\":");
        push_json_str(&mut out, self.code);
        out.push_str(",\"stage\":");
        push_json_str(&mut out, self.stage);
        out.push_str(",\"message\":");
        push_json_str(&mut out, &self.message);
        out.push('}');
        out
    }
}

impl fmt::Display for RuntimeRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}/{}] {}", self.stage, self.code, self.message)
    }
}

impl std::error::Error for RuntimeRefusal {}

/// Descriptor of one interactive input channel exposed by the runtime.
#[derive(Debug, Clone, PartialEq)]
pub struct InputDescriptor {
    /// Zero-based channel index passed to `set_input`.
    pub index: usize,
    /// Unique declared input name.
    pub name: String,
    /// Target entity name.
    pub target: String,
    /// Target port name on `target`.
    pub port: String,
    /// SI unit label derived from `dims`.
    pub unit: String,
    /// Six base SI dimension exponents `[L, M, T, Theta, I, N]`.
    pub dims: [i8; 6],
    /// Scalar channel width (`1` for scalar inputs).
    pub width: usize,
    /// Initial default value held before any `set_input` call.
    pub default_value: f64,
}

/// Descriptor of one probe slot in the flat snapshot buffer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProbeDescriptor {
    /// Zero-based probe index.
    pub index: usize,
    /// Unique declared probe name.
    pub name: String,
    /// Receiver/probe class (`field-quantity`, `thermocouple`, `custom`, ...).
    pub kind: String,
    /// Target entity name.
    pub target: String,
    /// Observed physical quantity name.
    pub quantity: String,
    /// SI unit label derived from `dims`.
    pub unit: String,
    /// Six base SI dimension exponents `[L, M, T, Theta, I, N]`.
    pub dims: [i8; 6],
    /// Starting offset in the flat `f64` snapshot buffer.
    pub offset: usize,
    /// Number of contiguous `f64` elements occupied by this probe.
    pub len: usize,
    /// Logical tensor shape (empty for scalar probes, `[n]` for vectors).
    pub shape: Vec<usize>,
}

/// Generated snapshot layout exported so front ends never hard-code offsets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnapshotDescriptor {
    /// Schema identifier ([`SNAPSHOT_DESCRIPTOR_SCHEMA`]).
    pub schema: String,
    /// Assembly name.
    pub assembly: String,
    /// Total `f64` length of the flat snapshot buffer.
    pub total_len: usize,
    /// Ordered probe descriptors covering `0..total_len`.
    pub probes: Vec<ProbeDescriptor>,
}

impl SnapshotDescriptor {
    /// Look up a probe descriptor by declared name.
    #[must_use]
    pub fn probe(&self, name: &str) -> Option<&ProbeDescriptor> {
        self.probes.iter().find(|p| p.name == name)
    }

    /// Encode this descriptor as deterministic canonical JSON.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::with_capacity(128 + self.probes.len() * 140);
        out.push_str("{\"schema\":");
        push_json_str(&mut out, &self.schema);
        out.push_str(",\"assembly\":");
        push_json_str(&mut out, &self.assembly);
        let _ = write!(out, ",\"total_len\":{},\"probes\":[", self.total_len);
        for (i, p) in self.probes.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "{{\"index\":{},\"name\":", p.index);
            push_json_str(&mut out, &p.name);
            out.push_str(",\"kind\":");
            push_json_str(&mut out, &p.kind);
            out.push_str(",\"target\":");
            push_json_str(&mut out, &p.target);
            out.push_str(",\"quantity\":");
            push_json_str(&mut out, &p.quantity);
            out.push_str(",\"unit\":");
            push_json_str(&mut out, &p.unit);
            let d = p.dims;
            let _ = write!(
                out,
                ",\"dims\":[{},{},{},{},{},{}],\"offset\":{},\"len\":{},\"shape\":[",
                d[0], d[1], d[2], d[3], d[4], d[5], p.offset, p.len
            );
            for (j, dim) in p.shape.iter().enumerate() {
                if j > 0 {
                    out.push(',');
                }
                let _ = write!(out, "{dim}");
            }
            out.push_str("]}");
        }
        out.push_str("]}");
        out
    }
}

/// Encode interactive input descriptors as canonical JSON.
#[must_use]
pub fn inputs_to_json(inputs: &[InputDescriptor]) -> String {
    let mut out = String::with_capacity(16 + inputs.len() * 140);
    out.push('[');
    for (i, inp) in inputs.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "{{\"index\":{},\"name\":", inp.index);
        push_json_str(&mut out, &inp.name);
        out.push_str(",\"target\":");
        push_json_str(&mut out, &inp.target);
        out.push_str(",\"port\":");
        push_json_str(&mut out, &inp.port);
        out.push_str(",\"unit\":");
        push_json_str(&mut out, &inp.unit);
        let d = inp.dims;
        let _ = write!(
            out,
            ",\"dims\":[{},{},{},{},{},{}],\"width\":{},\"default_value\":",
            d[0], d[1], d[2], d[3], d[4], d[5], inp.width
        );
        if inp.default_value.is_finite() {
            let v = if inp.default_value == 0.0 {
                0.0
            } else {
                inp.default_value
            };
            let _ = write!(out, "{v:.17e}}}");
        } else {
            out.push_str("null}");
        }
    }
    out.push(']');
    out
}

/// Canonical SI unit label for common physical dimensions (unified into `fs-fsim` in R2b).
#[must_use]
pub const fn dims_to_unit(dims: Dims) -> Option<&'static str> {
    match dims.0 {
        [1, 0, 0, 0, 0, 0] => Some("m"),
        [3, 0, 0, 0, 0, 0] => Some("m^3"),
        [1, 0, -1, 0, 0, 0] => Some("m/s"),
        [0, 0, 1, 0, 0, 0] => Some("s"),
        [0, 0, 0, 1, 0, 0] => Some("K"),
        [2, 1, -3, 0, 0, 0] => Some("W"),
        [0, 0, -1, 0, 0, 0] => Some("rad/s"),
        [0, 0, 0, 0, 0, 0] => Some("1"),
        [0, 1, 0, 0, 0, 0] => Some("kg"),
        [1, 1, -2, 0, 0, 0] => Some("N"),
        [2, 1, -2, 0, 0, 0] => Some("N*m"),
        [-1, 1, -2, 0, 0, 0] => Some("Pa"),
        [2, 1, -3, 0, -1, 0] => Some("V"),
        [0, 0, 0, 0, 1, 0] => Some("A"),
        _ => None,
    }
}

pub(crate) fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_json_writers_pin_exact_bytes() {
        let desc = SnapshotDescriptor {
            schema: SNAPSHOT_DESCRIPTOR_SCHEMA.to_owned(),
            assembly: "osc \"a\"\\1\n\x07".to_owned(),
            total_len: 4,
            probes: vec![
                ProbeDescriptor {
                    index: 0,
                    name: "pos".to_owned(),
                    kind: "field-quantity".to_owned(),
                    target: "mass".to_owned(),
                    quantity: "position".to_owned(),
                    unit: dims_to_unit(Dims([1, 0, 0, 0, 0, 0]))
                        .unwrap_or("SI")
                        .to_owned(),
                    dims: [1, 0, 0, 0, 0, 0],
                    offset: 0,
                    len: 1,
                    shape: Vec::new(),
                },
                ProbeDescriptor {
                    index: 1,
                    name: "vel".to_owned(),
                    kind: "custom".to_owned(),
                    target: "body".to_owned(),
                    quantity: "velocity".to_owned(),
                    unit: "m/s".to_owned(),
                    dims: [1, 0, -1, 0, 0, 0],
                    offset: 1,
                    len: 3,
                    shape: vec![3],
                },
            ],
        };
        assert_eq!(desc.probe("vel").unwrap().offset, 1);
        assert!(desc.probe("missing").is_none());
        let expected = concat!(
            r#"{"schema":"frankensim.runtime.snapshot_descriptor.v1","#,
            r#""assembly":"osc \"a\"\\1\n\u0007","total_len":4,"probes":["#,
            r#"{"index":0,"name":"pos","kind":"field-quantity","target":"mass","#,
            r#""quantity":"position","unit":"m","dims":[1,0,0,0,0,0],"offset":0,"len":1,"shape":[]},"#,
            r#"{"index":1,"name":"vel","kind":"custom","target":"body","#,
            r#""quantity":"velocity","unit":"m/s","dims":[1,0,-1,0,0,0],"offset":1,"len":3,"shape":[3]}]}"#,
        );
        assert_eq!(desc.to_json(), expected);

        let mut inputs = [InputDescriptor {
            index: 0,
            name: "drive".to_owned(),
            target: "mass".to_owned(),
            port: "mechanics".to_owned(),
            unit: dims_to_unit(Dims([1, 1, -2, 0, 0, 0]))
                .unwrap_or("SI")
                .to_owned(),
            dims: [1, 1, -2, 0, 0, 0],
            width: 1,
            default_value: -0.0,
        }];
        assert_eq!(
            inputs_to_json(&inputs),
            concat!(
                r#"[{"index":0,"name":"drive","target":"mass","port":"mechanics","#,
                r#""unit":"N","dims":[1,1,-2,0,0,0],"width":1,"default_value":0.00000000000000000e0}]"#,
            )
        );
        inputs[0].default_value = f64::NAN;
        assert!(inputs_to_json(&inputs).ends_with(r#""default_value":null}]"#));
        let r = RuntimeRefusal::new("invalid-dt", "step", "bad \"dt\"");
        assert_eq!(
            r.to_json(),
            r#"{"code":"invalid-dt","stage":"step","message":"bad \"dt\""}"#
        );
    }
}
