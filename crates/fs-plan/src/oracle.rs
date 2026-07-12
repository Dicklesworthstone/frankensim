//! Consumer wiring (plan §11.4 "consumers wired"): fitted cost models
//! behind fs-geom's `CostOracle` so the Rep Router plans with THIS
//! machine's measured history, and a fs-ledger `tune`-table loader so
//! models are rebuilt deterministically from ledger snapshots.

use std::collections::BTreeMap;

use crate::cost::{CostModel, CostObservation, CostRefusal, MAX_COST_OBSERVATIONS};
use fs_ledger::Ledger;

/// Maximum distinct converter edges tracked by one planner oracle.
pub const MAX_PLAN_ORACLE_EDGES: usize = 4_096;

/// Maximum absolute-error observations retained per converter edge.
pub const MAX_PLAN_ORACLE_ERROR_OBSERVATIONS: usize = MAX_COST_OBSERVATIONS;

/// Maximum UTF-8 bytes in one oracle edge identity.
pub const MAX_PLAN_ORACLE_EDGE_BYTES: usize = 4_096;

/// Why the planner oracle refused configuration or evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PlanOracleError {
    /// Edge identity is empty, oversized, or contains controls.
    InvalidEdge,
    /// A reference size is outside the positive finite model domain.
    InvalidReferenceSize,
    /// Re-registering an observed edge attempted to change its size domain.
    ReferenceSizeConflict,
    /// A reused edge name attempted to cross converter-specification domains.
    SpecificationConflict,
    /// An executed edge was never registered with an explicit reference size.
    UnregisteredEdge,
    /// The measured absolute error is not finite and nonnegative.
    InvalidError,
    /// The bounded distinct-edge budget is exhausted.
    EdgeLimit {
        /// Maximum distinct edges.
        limit: usize,
    },
    /// The bounded error-history budget is exhausted.
    ErrorObservationLimit {
        /// Maximum error observations per edge.
        limit: usize,
    },
    /// The per-edge cost model refused the observation.
    Cost(CostRefusal),
}

impl core::fmt::Display for PlanOracleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidEdge => write!(
                f,
                "edge identity must be nonempty, control-free, and at most {MAX_PLAN_ORACLE_EDGE_BYTES} bytes"
            ),
            Self::InvalidReferenceSize => {
                write!(f, "edge reference size must be positive and finite")
            }
            Self::ReferenceSizeConflict => write!(
                f,
                "cannot change an edge reference size after observations exist"
            ),
            Self::SpecificationConflict => write!(
                f,
                "edge identity is already registered for a different converter specification"
            ),
            Self::UnregisteredEdge => write!(
                f,
                "edge is not registered; register its explicit reference size before recording"
            ),
            Self::InvalidError => {
                write!(f, "measured absolute error must be nonnegative and finite")
            }
            Self::EdgeLimit { limit } => {
                write!(f, "planner oracle edge limit {limit} is exhausted")
            }
            Self::ErrorObservationLimit { limit } => write!(
                f,
                "planner oracle error-observation limit {limit} is exhausted"
            ),
            Self::Cost(error) => write!(f, "cost observation refused: {error}"),
        }
    }
}

impl core::error::Error for PlanOracleError {}

impl From<PlanOracleError> for fs_geom::CostOracleError {
    fn from(error: PlanOracleError) -> Self {
        let problem = error.to_string();
        match error {
            PlanOracleError::InvalidEdge
            | PlanOracleError::UnregisteredEdge
            | PlanOracleError::ReferenceSizeConflict
            | PlanOracleError::SpecificationConflict => Self::InvalidEdge { problem },
            PlanOracleError::InvalidReferenceSize => Self::InvalidMeasurement {
                field: "reference_size",
                problem,
            },
            PlanOracleError::InvalidError => Self::InvalidMeasurement {
                field: "error_abs",
                problem,
            },
            PlanOracleError::EdgeLimit { limit } => Self::CapacityExceeded {
                resource: "edges",
                limit,
            },
            PlanOracleError::ErrorObservationLimit { limit } => Self::CapacityExceeded {
                resource: "error_observations",
                limit,
            },
            PlanOracleError::Cost(CostRefusal::BadInput) => Self::InvalidMeasurement {
                field: "cost_s",
                problem,
            },
            PlanOracleError::Cost(CostRefusal::ObservationLimit { limit }) => {
                Self::CapacityExceeded {
                    resource: "cost_observations",
                    limit,
                }
            }
            PlanOracleError::Cost(_) => Self::Backend { problem },
        }
    }
}

#[derive(Debug, Clone)]
struct RegisteredPlanEdge {
    spec: fs_geom::ConverterSpec,
    reference_size: f64,
    model: CostModel,
}

/// A [`fs_geom::CostOracle`] backed by spec-scoped per-edge quantile cost
/// models. Each converter is registered with the reference problem size its
/// routing requests are quoted at; sealed execution observations feed the
/// online refits.
#[derive(Debug, Clone, Default)]
pub struct PlanCostOracle {
    models: BTreeMap<String, RegisteredPlanEdge>,
    errors: BTreeMap<String, Vec<f64>>,
}

impl PlanCostOracle {
    /// An empty oracle.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register an edge with its reference size (idempotent; keeps
    /// existing observations).
    pub fn register_edge(
        &mut self,
        spec: &fs_geom::ConverterSpec,
        reference_size: f64,
    ) -> Result<(), PlanOracleError> {
        Self::validate_edge(&spec.name)?;
        if !reference_size.is_finite() || reference_size <= 0.0 {
            return Err(PlanOracleError::InvalidReferenceSize);
        }
        if let Some(registered) = self.models.get_mut(&spec.name) {
            if registered.spec != *spec {
                return Err(PlanOracleError::SpecificationConflict);
            }
            if registered.model.n_obs() > 0
                && registered.reference_size.to_bits() != reference_size.to_bits()
            {
                return Err(PlanOracleError::ReferenceSizeConflict);
            }
            registered.reference_size = reference_size;
            return Ok(());
        }
        if self.models.len() >= MAX_PLAN_ORACLE_EDGES {
            return Err(PlanOracleError::EdgeLimit {
                limit: MAX_PLAN_ORACLE_EDGES,
            });
        }
        self.models.insert(
            spec.name.clone(),
            RegisteredPlanEdge {
                spec: spec.clone(),
                reference_size,
                model: CostModel::new(),
            },
        );
        Ok(())
    }

    /// The fitted model for an edge, if any.
    #[must_use]
    pub fn model(&self, edge: &str) -> Option<&CostModel> {
        self.models.get(edge).map(|registered| &registered.model)
    }

    fn validate_edge(edge: &str) -> Result<(), PlanOracleError> {
        if edge.is_empty()
            || edge.len() > MAX_PLAN_ORACLE_EDGE_BYTES
            || edge.chars().any(char::is_control)
        {
            Err(PlanOracleError::InvalidEdge)
        } else {
            Ok(())
        }
    }

    /// Record one edge atomically: cost and error histories either both
    /// advance, or neither does.
    fn record_one(
        &mut self,
        observation: &fs_geom::ValidatedEdgeObservation,
    ) -> Result<(), PlanOracleError> {
        let edge = observation.edge();
        let cost_s = observation.cost_s();
        let error_abs = observation.conservative_error_abs();
        Self::validate_edge(edge)?;
        if !error_abs.is_finite() || error_abs < 0.0 {
            return Err(PlanOracleError::InvalidError);
        }
        let Some(registered) = self.models.get(edge) else {
            return Err(PlanOracleError::UnregisteredEdge);
        };
        if registered.spec != *observation.spec() {
            return Err(PlanOracleError::SpecificationConflict);
        }
        let prior_errors = self.errors.get(edge).map_or(0, Vec::len);
        if prior_errors >= MAX_PLAN_ORACLE_ERROR_OBSERVATIONS {
            return Err(PlanOracleError::ErrorObservationLimit {
                limit: MAX_PLAN_ORACLE_ERROR_OBSERVATIONS,
            });
        }
        let mut candidate_model = registered.model.clone();
        candidate_model
            .observe(CostObservation {
                size: registered.reference_size,
                cost_s,
            })
            .map_err(PlanOracleError::Cost)?;
        let mut candidate_errors = self.errors.get(edge).cloned().unwrap_or_default();
        let insertion = candidate_errors
            .binary_search_by(|value| value.total_cmp(&error_abs))
            .unwrap_or_else(|index| index);
        candidate_errors.insert(insertion, error_abs);

        let entry = self
            .models
            .get_mut(edge)
            .ok_or(PlanOracleError::UnregisteredEdge)?;
        entry.model = candidate_model;
        self.errors.insert(edge.to_string(), candidate_errors);
        Ok(())
    }
}

impl fs_geom::CostOracle for PlanCostOracle {
    fn measured_cost_s(
        &self,
        spec: &fs_geom::ConverterSpec,
    ) -> Result<Option<f64>, fs_geom::CostOracleError> {
        let Some(registered) = self.models.get(&spec.name) else {
            return Ok(None);
        };
        if registered.spec != *spec {
            return Err(PlanOracleError::SpecificationConflict.into());
        }
        match registered.model.predict(registered.reference_size) {
            Ok(prediction) => Ok(Some(prediction.p50)),
            Err(CostRefusal::InsufficientData { .. }) => Ok(None),
            Err(error) => Err(fs_geom::CostOracleError::Backend {
                problem: format!("registered cost model prediction failed: {error}"),
            }),
        }
    }

    fn measured_error_abs(
        &self,
        spec: &fs_geom::ConverterSpec,
    ) -> Result<Option<f64>, fs_geom::CostOracleError> {
        let Some(registered) = self.models.get(&spec.name) else {
            return Ok(None);
        };
        if registered.spec != *spec {
            return Err(PlanOracleError::SpecificationConflict.into());
        }
        Ok(self
            .errors
            .get(&spec.name)
            .and_then(|errs| errs.last().copied()))
    }

    fn record_batch(
        &mut self,
        observations: &[fs_geom::ValidatedEdgeObservation],
    ) -> Result<(), fs_geom::CostOracleError> {
        let mut candidate = self.clone();
        for observation in observations {
            candidate
                .record_one(observation)
                .map_err(fs_geom::CostOracleError::from)?;
        }
        *self = candidate;
        Ok(())
    }
}

/// Receipt schema written by `fs-roofline::Attainment::to_jsonl`.
pub const ROOFLINE_RECEIPT_VERSION: u64 = 3;

/// Production tune-row parameter schema written by fs-roofline.
pub const ROOFLINE_ROW_SCHEMA: &str = "fs-roofline-ledger-row-v4";

/// Production tune shape prefix written by fs-roofline.
pub const ROOFLINE_TUNE_SHAPE_PREFIX: &str = "roofline-v7";

/// Exact byte width of fs-roofline's fingerprint + baseline machine key.
pub const ROOFLINE_MACHINE_KEY_BYTES: usize = 40;

/// Maximum receipt bytes decoded by the planner. This matches the ledger's
/// tune measurement admission bound, so parsing is bounded before allocation.
pub const MAX_ROOFLINE_RECEIPT_BYTES: usize = fs_ledger::MAX_TUNE_MEASURED_BYTES;

const MAX_RECEIPT_JSON_DEPTH: usize = 32;
const MAX_RECEIPT_JSON_NODES: usize = 65_536;
const MAX_RECEIPT_JSON_STRING_BYTES: usize = 256 * 1024;
const MAX_RECEIPT_JSON_CONTAINER_ITEMS: usize = 32_768;
const MAX_RECEIPT_JSON_NUMBER_BYTES: usize = 64;
const MAX_ROOFLINE_REPS: usize = 1_000;

/// Why an exact ledger tune row could not become cost-model evidence.
#[derive(Debug)]
pub enum TuneModelError {
    /// Ledger lookup or storage validation failed.
    Ledger(fs_ledger::LedgerError),
    /// No row exists at the exact `(kernel, shape_class, machine)` key.
    MissingRow,
    /// JSON syntax, bounds, types, or receipt schema were invalid.
    InvalidReceipt {
        /// Receipt/parameter location.
        field: String,
        /// Refusal diagnosis.
        problem: String,
    },
    /// A row's embedded identity disagreed with its exact ledger key.
    ScopeMismatch {
        /// Mismatched identity dimension.
        field: &'static str,
    },
    /// The decoded observation was refused by the numerical model.
    Cost(CostRefusal),
}

impl core::fmt::Display for TuneModelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Ledger(error) => write!(f, "tune ledger read failed: {error}"),
            Self::MissingRow => write!(
                f,
                "no tune row exists at the exact kernel/shape/machine key"
            ),
            Self::InvalidReceipt { field, problem } => {
                write!(f, "invalid roofline receipt at {field}: {problem}")
            }
            Self::ScopeMismatch { field } => write!(
                f,
                "roofline row embedded {field} does not match its exact ledger key"
            ),
            Self::Cost(error) => write!(f, "decoded tune observation refused: {error}"),
        }
    }
}

impl core::error::Error for TuneModelError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Ledger(error) => Some(error),
            Self::Cost(error) => Some(error),
            Self::MissingRow | Self::InvalidReceipt { .. } | Self::ScopeMismatch { .. } => None,
        }
    }
}

impl From<fs_ledger::LedgerError> for TuneModelError {
    fn from(error: fs_ledger::LedgerError) -> Self {
        Self::Ledger(error)
    }
}

impl From<CostRefusal> for TuneModelError {
    fn from(error: CostRefusal) -> Self {
        Self::Cost(error)
    }
}

#[derive(Debug)]
enum JsonValue {
    Null,
    Bool(bool),
    String(String),
    Number(String),
    Array(Vec<JsonValue>),
    Object(Vec<(String, JsonValue)>),
}

impl JsonValue {
    fn kind(&self) -> &'static str {
        match self {
            Self::Null => "null",
            Self::Bool(_) => "boolean",
            Self::String(_) => "string",
            Self::Number(_) => "number",
            Self::Array(_) => "array",
            Self::Object(_) => "object",
        }
    }
}

struct StrictJson<'a> {
    bytes: &'a [u8],
    offset: usize,
    nodes: usize,
}

impl<'a> StrictJson<'a> {
    fn parse(text: &'a str, max_bytes: usize) -> Result<JsonValue, TuneModelError> {
        if text.len() > max_bytes {
            return Err(invalid_receipt(
                "json",
                format!("{} bytes exceeds limit {max_bytes}", text.len()),
            ));
        }
        let mut parser = Self {
            bytes: text.as_bytes(),
            offset: 0,
            nodes: 0,
        };
        let value = parser.value(0)?;
        parser.whitespace();
        if parser.offset != parser.bytes.len() {
            return Err(parser.error("json", "trailing bytes after root value"));
        }
        Ok(value)
    }

    fn error(&self, field: &str, problem: impl Into<String>) -> TuneModelError {
        invalid_receipt(format!("{field} (byte {})", self.offset), problem.into())
    }

    fn whitespace(&mut self) {
        while self
            .bytes
            .get(self.offset)
            .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r'))
        {
            self.offset += 1;
        }
    }

    fn value(&mut self, depth: usize) -> Result<JsonValue, TuneModelError> {
        if depth > MAX_RECEIPT_JSON_DEPTH {
            return Err(self.error(
                "json",
                format!("nesting depth exceeds limit {MAX_RECEIPT_JSON_DEPTH}"),
            ));
        }
        self.nodes = self
            .nodes
            .checked_add(1)
            .ok_or_else(|| self.error("json", "node counter overflow"))?;
        if self.nodes > MAX_RECEIPT_JSON_NODES {
            return Err(self.error(
                "json",
                format!("node count exceeds limit {MAX_RECEIPT_JSON_NODES}"),
            ));
        }
        self.whitespace();
        match self.bytes.get(self.offset).copied() {
            Some(b'{') => self.object(depth),
            Some(b'[') => self.array(depth),
            Some(b'"') => self.string().map(JsonValue::String),
            Some(b'n') if self.consume_literal(b"null") => Ok(JsonValue::Null),
            Some(b't') if self.consume_literal(b"true") => Ok(JsonValue::Bool(true)),
            Some(b'f') if self.consume_literal(b"false") => Ok(JsonValue::Bool(false)),
            Some(byte) if byte.is_ascii_digit() || byte == b'-' => {
                self.number().map(JsonValue::Number)
            }
            _ => Err(self.error("json", "unexpected byte or end of input")),
        }
    }

    fn consume_literal(&mut self, literal: &[u8]) -> bool {
        if self
            .bytes
            .get(self.offset..)
            .is_some_and(|remaining| remaining.starts_with(literal))
        {
            self.offset += literal.len();
            true
        } else {
            false
        }
    }

    fn object(&mut self, depth: usize) -> Result<JsonValue, TuneModelError> {
        self.offset += 1;
        let mut fields = Vec::new();
        let mut keys = std::collections::BTreeSet::new();
        self.whitespace();
        if self.bytes.get(self.offset) == Some(&b'}') {
            self.offset += 1;
            return Ok(JsonValue::Object(fields));
        }
        loop {
            let key = self.string()?;
            self.whitespace();
            if self.bytes.get(self.offset) != Some(&b':') {
                return Err(self.error("object", "expected ':' after key"));
            }
            self.offset += 1;
            let value = self.value(depth + 1)?;
            if !keys.insert(key.clone()) {
                return Err(self.error("object", format!("duplicate key {key:?}")));
            }
            fields.push((key, value));
            if fields.len() > MAX_RECEIPT_JSON_CONTAINER_ITEMS {
                return Err(self.error(
                    "object",
                    format!("member count exceeds limit {MAX_RECEIPT_JSON_CONTAINER_ITEMS}"),
                ));
            }
            self.whitespace();
            match self.bytes.get(self.offset) {
                Some(b',') => {
                    self.offset += 1;
                    self.whitespace();
                }
                Some(b'}') => {
                    self.offset += 1;
                    return Ok(JsonValue::Object(fields));
                }
                _ => return Err(self.error("object", "expected ',' or '}'")),
            }
        }
    }

    fn array(&mut self, depth: usize) -> Result<JsonValue, TuneModelError> {
        self.offset += 1;
        let mut values = Vec::new();
        self.whitespace();
        if self.bytes.get(self.offset) == Some(&b']') {
            self.offset += 1;
            return Ok(JsonValue::Array(values));
        }
        loop {
            values.push(self.value(depth + 1)?);
            if values.len() > MAX_RECEIPT_JSON_CONTAINER_ITEMS {
                return Err(self.error(
                    "array",
                    format!("element count exceeds limit {MAX_RECEIPT_JSON_CONTAINER_ITEMS}"),
                ));
            }
            self.whitespace();
            match self.bytes.get(self.offset) {
                Some(b',') => self.offset += 1,
                Some(b']') => {
                    self.offset += 1;
                    return Ok(JsonValue::Array(values));
                }
                _ => return Err(self.error("array", "expected ',' or ']'")),
            }
        }
    }

    fn number(&mut self) -> Result<String, TuneModelError> {
        let start = self.offset;
        if self.bytes.get(self.offset) == Some(&b'-') {
            self.offset += 1;
        }
        match self.bytes.get(self.offset) {
            Some(b'0') => {
                self.offset += 1;
                if self.bytes.get(self.offset).is_some_and(u8::is_ascii_digit) {
                    return Err(self.error("number", "leading zeros are not canonical JSON"));
                }
            }
            Some(b'1'..=b'9') => {
                while self.bytes.get(self.offset).is_some_and(u8::is_ascii_digit) {
                    self.offset += 1;
                }
            }
            _ => return Err(self.error("number", "missing integer digits")),
        }
        if self.bytes.get(self.offset) == Some(&b'.') {
            self.offset += 1;
            let fraction_start = self.offset;
            while self.bytes.get(self.offset).is_some_and(u8::is_ascii_digit) {
                self.offset += 1;
            }
            if self.offset == fraction_start {
                return Err(self.error("number", "fraction has no digits"));
            }
        }
        if self
            .bytes
            .get(self.offset)
            .is_some_and(|byte| matches!(byte, b'e' | b'E'))
        {
            self.offset += 1;
            if self
                .bytes
                .get(self.offset)
                .is_some_and(|byte| matches!(byte, b'+' | b'-'))
            {
                self.offset += 1;
            }
            let exponent_start = self.offset;
            while self.bytes.get(self.offset).is_some_and(u8::is_ascii_digit) {
                self.offset += 1;
            }
            if self.offset == exponent_start {
                return Err(self.error("number", "exponent has no digits"));
            }
        }
        if self.offset - start > MAX_RECEIPT_JSON_NUMBER_BYTES {
            return Err(self.error(
                "number",
                format!("token exceeds {MAX_RECEIPT_JSON_NUMBER_BYTES} bytes"),
            ));
        }
        let text = self
            .bytes
            .get(start..self.offset)
            .and_then(|bytes| core::str::from_utf8(bytes).ok())
            .ok_or_else(|| self.error("number", "number is not UTF-8"))?;
        let finite = text.parse::<f64>().is_ok_and(f64::is_finite);
        if !finite {
            return Err(self.error("number", "number is not finite"));
        }
        Ok(text.to_string())
    }

    fn string(&mut self) -> Result<String, TuneModelError> {
        self.whitespace();
        if self.bytes.get(self.offset) != Some(&b'"') {
            return Err(self.error("string", "expected opening quote"));
        }
        self.offset += 1;
        let mut output = String::new();
        loop {
            match self.bytes.get(self.offset).copied() {
                None => return Err(self.error("string", "unterminated string")),
                Some(b'"') => {
                    self.offset += 1;
                    return Ok(output);
                }
                Some(b'\\') => {
                    self.offset += 1;
                    let escape = self
                        .bytes
                        .get(self.offset)
                        .copied()
                        .ok_or_else(|| self.error("string", "unterminated escape"))?;
                    self.offset += 1;
                    match escape {
                        b'"' => output.push('"'),
                        b'\\' => output.push('\\'),
                        b'/' => output.push('/'),
                        b'b' => output.push('\u{0008}'),
                        b'f' => output.push('\u{000c}'),
                        b'n' => output.push('\n'),
                        b'r' => output.push('\r'),
                        b't' => output.push('\t'),
                        b'u' => output.push(self.unicode_escape()?),
                        _ => return Err(self.error("string", "invalid escape")),
                    }
                }
                Some(byte) if byte < 0x20 => {
                    return Err(self.error("string", "unescaped control character"));
                }
                Some(byte) => {
                    let width = if byte < 0x80 {
                        1
                    } else if byte >> 5 == 0b110 {
                        2
                    } else if byte >> 4 == 0b1110 {
                        3
                    } else if byte >> 3 == 0b11110 {
                        4
                    } else {
                        return Err(self.error("string", "invalid UTF-8 lead byte"));
                    };
                    let chunk = self
                        .bytes
                        .get(self.offset..self.offset + width)
                        .and_then(|bytes| core::str::from_utf8(bytes).ok())
                        .ok_or_else(|| self.error("string", "invalid UTF-8"))?;
                    output.push_str(chunk);
                    self.offset += width;
                }
            }
            if output.len() > MAX_RECEIPT_JSON_STRING_BYTES {
                return Err(self.error(
                    "string",
                    format!("decoded string exceeds {MAX_RECEIPT_JSON_STRING_BYTES} bytes"),
                ));
            }
        }
    }

    fn hex_escape_unit(&mut self) -> Result<u16, TuneModelError> {
        let end = self
            .offset
            .checked_add(4)
            .ok_or_else(|| self.error("string", "unicode escape offset overflow"))?;
        let digits = self
            .bytes
            .get(self.offset..end)
            .and_then(|bytes| core::str::from_utf8(bytes).ok())
            .ok_or_else(|| self.error("string", "short unicode escape"))?;
        let unit = u16::from_str_radix(digits, 16)
            .map_err(|_| self.error("string", "invalid unicode escape"))?;
        self.offset = end;
        Ok(unit)
    }

    fn unicode_escape(&mut self) -> Result<char, TuneModelError> {
        let first = self.hex_escape_unit()?;
        let scalar = if (0xd800..=0xdbff).contains(&first) {
            if self.bytes.get(self.offset..self.offset + 2) != Some(b"\\u") {
                return Err(self.error("string", "high surrogate lacks low surrogate"));
            }
            self.offset += 2;
            let second = self.hex_escape_unit()?;
            if !(0xdc00..=0xdfff).contains(&second) {
                return Err(self.error("string", "invalid low surrogate"));
            }
            0x1_0000 + ((u32::from(first) - 0xd800) << 10) + (u32::from(second) - 0xdc00)
        } else if (0xdc00..=0xdfff).contains(&first) {
            return Err(self.error("string", "unpaired low surrogate"));
        } else {
            u32::from(first)
        };
        char::from_u32(scalar)
            .ok_or_else(|| self.error("string", "unicode escape is not a scalar value"))
    }
}

fn invalid_receipt(field: impl Into<String>, problem: impl Into<String>) -> TuneModelError {
    TuneModelError::InvalidReceipt {
        field: field.into(),
        problem: problem.into(),
    }
}

struct ObjectFields {
    what: String,
    fields: BTreeMap<String, JsonValue>,
}

impl ObjectFields {
    fn new(value: JsonValue, what: impl Into<String>) -> Result<Self, TuneModelError> {
        let what = what.into();
        let JsonValue::Object(fields) = value else {
            return Err(invalid_receipt(
                &what,
                format!("expected object, got {}", value.kind()),
            ));
        };
        Ok(Self {
            what,
            fields: fields.into_iter().collect(),
        })
    }

    fn take(&mut self, field: &'static str) -> Result<JsonValue, TuneModelError> {
        self.fields
            .remove(field)
            .ok_or_else(|| invalid_receipt(format!("{}.{}", self.what, field), "missing field"))
    }

    fn finish(self) -> Result<(), TuneModelError> {
        if let Some(field) = self.fields.keys().next() {
            Err(invalid_receipt(
                format!("{}.{}", self.what, field),
                "unknown field",
            ))
        } else {
            Ok(())
        }
    }
}

fn expect_string(value: JsonValue, field: &str) -> Result<String, TuneModelError> {
    match value {
        JsonValue::String(value) => Ok(value),
        other => Err(invalid_receipt(
            field,
            format!("expected string, got {}", other.kind()),
        )),
    }
}

fn expect_number(value: JsonValue, field: &str) -> Result<f64, TuneModelError> {
    match value {
        JsonValue::Number(value) => value
            .parse::<f64>()
            .ok()
            .filter(|number| number.is_finite())
            .ok_or_else(|| invalid_receipt(field, "expected finite number")),
        other => Err(invalid_receipt(
            field,
            format!("expected number, got {}", other.kind()),
        )),
    }
}

fn expect_bool(value: JsonValue, field: &str) -> Result<bool, TuneModelError> {
    match value {
        JsonValue::Bool(value) => Ok(value),
        other => Err(invalid_receipt(
            field,
            format!("expected boolean, got {}", other.kind()),
        )),
    }
}

fn expect_u64(value: JsonValue, field: &str) -> Result<u64, TuneModelError> {
    match value {
        JsonValue::Number(value)
            if !value.is_empty() && value.bytes().all(|byte| byte.is_ascii_digit()) =>
        {
            value
                .parse::<u64>()
                .map_err(|_| invalid_receipt(field, "unsigned integer is out of range"))
        }
        JsonValue::Number(_) => Err(invalid_receipt(field, "expected unsigned integer")),
        other => Err(invalid_receipt(
            field,
            format!("expected number, got {}", other.kind()),
        )),
    }
}

fn expect_array(value: JsonValue, field: &str) -> Result<Vec<JsonValue>, TuneModelError> {
    match value {
        JsonValue::Array(values) => Ok(values),
        other => Err(invalid_receipt(
            field,
            format!("expected array, got {}", other.kind()),
        )),
    }
}

fn expect_hex_u64(value: JsonValue, field: &str) -> Result<u64, TuneModelError> {
    let text = expect_string(value, field)?;
    if text.len() != 16
        || !text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(invalid_receipt(
            field,
            "expected exactly 16 lowercase hexadecimal digits",
        ));
    }
    u64::from_str_radix(&text, 16).map_err(|_| invalid_receipt(field, "hex value is out of range"))
}

fn expect_hash(value: JsonValue, field: &str) -> Result<String, TuneModelError> {
    let text = expect_string(value, field)?;
    if text.len() == 64
        && text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        Ok(text)
    } else {
        Err(invalid_receipt(
            field,
            "expected exactly 64 lowercase hexadecimal digits",
        ))
    }
}

fn finite_from_bits(
    value: JsonValue,
    field: &str,
    positive: bool,
) -> Result<(u64, f64), TuneModelError> {
    let bits = expect_hex_u64(value, field)?;
    let decoded = f64::from_bits(bits);
    let admitted = decoded.is_finite()
        && if positive {
            decoded > 0.0
        } else {
            decoded >= 0.0
        };
    if admitted {
        Ok((bits, decoded))
    } else {
        Err(invalid_receipt(
            field,
            if positive {
                "bits must encode a positive finite value"
            } else {
                "bits must encode a nonnegative finite value"
            },
        ))
    }
}

#[derive(Debug)]
struct ReceiptObservation {
    kernel: String,
    version: String,
    machine: u64,
    observations: Vec<CostObservation>,
    reps: u64,
}

#[derive(Debug)]
struct RowBinding {
    op: u64,
    run_receipt: String,
    payload_artifact: String,
    dependency_receipt_artifact: String,
    dependency_receipt_digest: String,
    baseline_hash: String,
    build_identity: String,
    reps: u64,
}

struct DecodedMeasurement {
    elements: u64,
    sample_times: Vec<f64>,
    decision_count: usize,
    median_bits: u64,
    median_seconds: f64,
    p25_bits: u64,
    p75_bits: u64,
    dispersion_bits: u64,
}

struct ReceiptMetrics {
    rate_bits: u64,
    reps: u64,
}

fn decode_row_binding(text: &str) -> Result<RowBinding, TuneModelError> {
    let mut object = ObjectFields::new(
        StrictJson::parse(text, fs_ledger::MAX_TUNE_PARAMS_BYTES)?,
        "params",
    )?;
    let schema = expect_string(object.take("schema")?, "params.schema")?;
    if schema != ROOFLINE_ROW_SCHEMA {
        return Err(invalid_receipt(
            "params.schema",
            format!("unsupported schema {schema:?}"),
        ));
    }
    let op = expect_u64(object.take("op")?, "params.op")?;
    if op == 0 || i64::try_from(op).is_err() {
        return Err(invalid_receipt(
            "params.op",
            "operation id must be in 1..=i64::MAX",
        ));
    }
    let run_receipt = expect_hash(object.take("run_receipt")?, "params.run_receipt")?;
    let payload_artifact =
        expect_hash(object.take("payload_artifact")?, "params.payload_artifact")?;
    let dependency_receipt_artifact = expect_hash(
        object.take("dependency_receipt_artifact")?,
        "params.dependency_receipt_artifact",
    )?;
    let dependency_receipt_digest = expect_hash(
        object.take("dependency_receipt_digest")?,
        "params.dependency_receipt_digest",
    )?;
    let baseline_hash = expect_hash(object.take("baseline_hash")?, "params.baseline_hash")?;
    let build_identity = expect_hash(object.take("build_identity")?, "params.build_identity")?;
    let reps = expect_u64(object.take("reps")?, "params.reps")?;
    if reps == 0 || reps > MAX_ROOFLINE_REPS as u64 {
        return Err(invalid_receipt(
            "params.reps",
            format!("must be in 1..={MAX_ROOFLINE_REPS}"),
        ));
    }
    for field in [
        "post_bandwidth_single_bits",
        "post_bandwidth_all_core_bits",
        "post_peak_single_bits",
        "post_peak_all_core_bits",
    ] {
        let _ = finite_from_bits(object.take(field)?, &format!("params.{field}"), true)?;
    }
    object.finish()?;
    Ok(RowBinding {
        op,
        run_receipt,
        payload_artifact,
        dependency_receipt_artifact,
        dependency_receipt_digest,
        baseline_hash,
        build_identity,
        reps,
    })
}

fn decode_receipt_axes(value: JsonValue) -> Result<(), TuneModelError> {
    let mut axes = ObjectFields::new(value, "receipt.axes")?;
    let logical_cpus = expect_u64(axes.take("logical_cpus")?, "receipt.axes.logical_cpus")?;
    if logical_cpus == 0 || u32::try_from(logical_cpus).is_err() {
        return Err(invalid_receipt(
            "receipt.axes.logical_cpus",
            "must be in 1..=u32::MAX",
        ));
    }
    for field in [
        "bandwidth_single_bits",
        "bandwidth_all_core_bits",
        "peak_single_bits",
        "peak_all_core_bits",
    ] {
        let _ = finite_from_bits(axes.take(field)?, &format!("receipt.axes.{field}"), true)?;
    }
    axes.finish()
}

fn decode_receipt_spec(value: JsonValue) -> Result<(), TuneModelError> {
    let mut spec = ObjectFields::new(value, "receipt.spec")?;
    let _ = finite_from_bits(
        spec.take("bytes_per_elem_bits")?,
        "receipt.spec.bytes_per_elem_bits",
        true,
    )?;
    let _ = finite_from_bits(
        spec.take("flops_per_elem_bits")?,
        "receipt.spec.flops_per_elem_bits",
        false,
    )?;
    let threading = expect_string(spec.take("threading")?, "receipt.spec.threading")?;
    if !matches!(threading.as_str(), "single_thread" | "all_core") {
        return Err(invalid_receipt(
            "receipt.spec.threading",
            "unknown threading class",
        ));
    }
    let target_axis = expect_string(spec.take("target_axis")?, "receipt.spec.target_axis")?;
    if !matches!(
        target_axis.as_str(),
        "binding_roof" | "compute_peak" | "memory_bandwidth"
    ) {
        return Err(invalid_receipt(
            "receipt.spec.target_axis",
            "unknown target axis",
        ));
    }
    match spec.take("target_fraction_bits")? {
        JsonValue::Null => {}
        value => {
            let (_, target) = finite_from_bits(value, "receipt.spec.target_fraction_bits", true)?;
            if target > 1.0 {
                return Err(invalid_receipt(
                    "receipt.spec.target_fraction_bits",
                    "target fraction exceeds 1",
                ));
            }
        }
    }
    spec.finish()
}

fn validate_decision_hashes(values: &[JsonValue]) -> Result<(), TuneModelError> {
    for (index, value) in values.iter().enumerate() {
        match value {
            JsonValue::Null => {}
            JsonValue::String(hash)
                if hash.len() == 64
                    && hash
                        .bytes()
                        .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)) => {}
            other => {
                return Err(invalid_receipt(
                    format!("receipt.measurement.decision_binding_hashes[{index}]"),
                    format!(
                        "expected null or lowercase hash string, got {}",
                        other.kind()
                    ),
                ));
            }
        }
    }
    Ok(())
}

fn decode_receipt_measurement(value: JsonValue) -> Result<DecodedMeasurement, TuneModelError> {
    let mut measurement = ObjectFields::new(value, "receipt.measurement")?;
    let origin = expect_string(measurement.take("origin")?, "receipt.measurement.origin")?;
    if origin != "timed" {
        return Err(invalid_receipt(
            "receipt.measurement.origin",
            "only timed production receipts carry cost evidence",
        ));
    }
    let elements = expect_u64(
        measurement.take("elements")?,
        "receipt.measurement.elements",
    )?;
    if elements == 0 || usize::try_from(elements).is_err() {
        return Err(invalid_receipt(
            "receipt.measurement.elements",
            "must be in 1..=usize::MAX",
        ));
    }
    let warmup_runs = expect_u64(
        measurement.take("warmup_runs")?,
        "receipt.measurement.warmup_runs",
    )?;
    if warmup_runs > MAX_ROOFLINE_REPS as u64 {
        return Err(invalid_receipt(
            "receipt.measurement.warmup_runs",
            format!("exceeds limit {MAX_ROOFLINE_REPS}"),
        ));
    }
    let sample_bits = expect_array(
        measurement.take("sample_seconds_bits")?,
        "receipt.measurement.sample_seconds_bits",
    )?;
    if sample_bits.is_empty() || sample_bits.len() > MAX_ROOFLINE_REPS {
        return Err(invalid_receipt(
            "receipt.measurement.sample_seconds_bits",
            format!("length must be in 1..={MAX_ROOFLINE_REPS}"),
        ));
    }
    let mut sample_times = Vec::with_capacity(sample_bits.len());
    for (index, value) in sample_bits.into_iter().enumerate() {
        let (_, seconds) = finite_from_bits(
            value,
            &format!("receipt.measurement.sample_seconds_bits[{index}]"),
            true,
        )?;
        sample_times.push(seconds);
    }
    let decision_hashes = expect_array(
        measurement.take("decision_binding_hashes")?,
        "receipt.measurement.decision_binding_hashes",
    )?;
    if decision_hashes.len() > MAX_ROOFLINE_REPS {
        return Err(invalid_receipt(
            "receipt.measurement.decision_binding_hashes",
            format!("length exceeds limit {MAX_ROOFLINE_REPS}"),
        ));
    }
    validate_decision_hashes(&decision_hashes)?;
    let (median_bits, median_seconds) = finite_from_bits(
        measurement.take("median_seconds_bits")?,
        "receipt.measurement.median_seconds_bits",
        true,
    )?;
    let (p25_bits, p25_seconds) = finite_from_bits(
        measurement.take("p25_seconds_bits")?,
        "receipt.measurement.p25_seconds_bits",
        true,
    )?;
    let (p75_bits, p75_seconds) = finite_from_bits(
        measurement.take("p75_seconds_bits")?,
        "receipt.measurement.p75_seconds_bits",
        true,
    )?;
    let (dispersion_bits, _) = finite_from_bits(
        measurement.take("dispersion_bits")?,
        "receipt.measurement.dispersion_bits",
        false,
    )?;
    if p25_seconds > median_seconds || median_seconds > p75_seconds {
        return Err(invalid_receipt(
            "receipt.measurement.quantiles",
            "seconds quantiles are not ordered",
        ));
    }
    measurement.finish()?;
    Ok(DecodedMeasurement {
        elements,
        sample_times,
        decision_count: decision_hashes.len(),
        median_bits,
        median_seconds,
        p25_bits,
        p75_bits,
        dispersion_bits,
    })
}

fn decode_receipt_metrics(
    receipt: &mut ObjectFields,
    measurement: &DecodedMeasurement,
) -> Result<ReceiptMetrics, TuneModelError> {
    match receipt.take("execution")? {
        JsonValue::Null | JsonValue::Object(_) => {}
        other => {
            return Err(invalid_receipt(
                "receipt.execution",
                format!("expected null or object, got {}", other.kind()),
            ));
        }
    }
    let (rate_bits, _) = finite_from_bits(
        receipt.take("elems_per_sec_bits")?,
        "receipt.elems_per_sec_bits",
        true,
    )?;
    for (field, positive) in [
        ("gbs_bits", false),
        ("gflops_bits", false),
        ("limit_elems_per_sec_bits", true),
        ("attainment_bits", false),
        ("target_attainment_bits", false),
    ] {
        let _ = finite_from_bits(receipt.take(field)?, &format!("receipt.{field}"), positive)?;
    }
    let (dispersion_bits, _) = finite_from_bits(
        receipt.take("dispersion_bits")?,
        "receipt.dispersion_bits",
        false,
    )?;
    if dispersion_bits != measurement.dispersion_bits {
        return Err(invalid_receipt(
            "receipt.dispersion_bits",
            "does not match measurement dispersion",
        ));
    }
    for (field, positive) in [
        ("elems_per_sec_display", true),
        ("gbs", false),
        ("gflops", false),
        ("limit_elems_per_sec", true),
        ("attainment", false),
        ("target_attainment", false),
        ("dispersion", false),
    ] {
        let value = expect_number(receipt.take(field)?, &format!("receipt.{field}"))?;
        if if positive { value <= 0.0 } else { value < 0.0 } {
            return Err(invalid_receipt(
                format!("receipt.{field}"),
                if positive {
                    "must be positive"
                } else {
                    "must be nonnegative"
                },
            ));
        }
    }
    let roof = expect_string(receipt.take("roof")?, "receipt.roof")?;
    if !matches!(roof.as_str(), "bandwidth" | "compute") {
        return Err(invalid_receipt("receipt.roof", "unknown roof side"));
    }
    let reps = expect_u64(receipt.take("reps")?, "receipt.reps")?;
    if reps == 0 || reps > MAX_ROOFLINE_REPS as u64 {
        return Err(invalid_receipt(
            "receipt.reps",
            format!("must be in 1..={MAX_ROOFLINE_REPS}"),
        ));
    }
    let expected_reps = sample_bits_len(reps)?;
    if measurement.sample_times.len() != expected_reps
        || measurement.decision_count != expected_reps
    {
        return Err(invalid_receipt(
            "receipt.measurement",
            "sample and decision-binding lengths must both match reps",
        ));
    }
    let verdict = expect_string(receipt.take("verdict")?, "receipt.verdict")?;
    if !matches!(verdict.as_str(), "within_band" | "below_band" | "no_target") {
        return Err(invalid_receipt(
            "receipt.verdict",
            "non-citable or unknown verdict",
        ));
    }
    if !matches!(receipt.take("invalid_reason")?, JsonValue::Null) {
        return Err(invalid_receipt(
            "receipt.invalid_reason",
            "citable tune evidence must have null invalid_reason",
        ));
    }
    Ok(ReceiptMetrics { rate_bits, reps })
}

fn validate_rederived_metrics(
    measurement: &DecodedMeasurement,
    rate_bits: u64,
) -> Result<f64, TuneModelError> {
    let mut sorted_times = measurement.sample_times.clone();
    sorted_times.sort_by(f64::total_cmp);
    let recomputed_median = nearest_rank(&sorted_times, 0.5)?;
    let recomputed_p25 = nearest_rank(&sorted_times, 0.25)?;
    let recomputed_p75 = nearest_rank(&sorted_times, 0.75)?;
    let recomputed_dispersion = (recomputed_p75 - recomputed_p25) / recomputed_median;
    if recomputed_median.to_bits() != measurement.median_bits
        || recomputed_p25.to_bits() != measurement.p25_bits
        || recomputed_p75.to_bits() != measurement.p75_bits
        || recomputed_dispersion.to_bits() != measurement.dispersion_bits
    {
        return Err(invalid_receipt(
            "receipt.measurement.statistics",
            "stored median/p25/p75/dispersion do not rederive from sample_seconds_bits",
        ));
    }
    let size = measurement.elements as f64;
    if (size / measurement.median_seconds).to_bits() != rate_bits {
        return Err(invalid_receipt(
            "receipt.elems_per_sec_bits",
            "does not rederive from measurement elements / median seconds",
        ));
    }
    Ok(size)
}

fn decode_receipt(text: &str) -> Result<ReceiptObservation, TuneModelError> {
    let mut receipt = ObjectFields::new(
        StrictJson::parse(text, MAX_ROOFLINE_RECEIPT_BYTES)?,
        "receipt",
    )?;
    let version = expect_u64(receipt.take("receipt_version")?, "receipt.receipt_version")?;
    if version != ROOFLINE_RECEIPT_VERSION {
        return Err(invalid_receipt(
            "receipt.receipt_version",
            format!("unsupported version {version}"),
        ));
    }
    let kernel = expect_string(receipt.take("kernel")?, "receipt.kernel")?;
    let version = expect_string(receipt.take("version")?, "receipt.version")?;
    if kernel.is_empty() || version.is_empty() {
        return Err(invalid_receipt(
            "receipt.identity",
            "kernel and version must be nonempty",
        ));
    }
    let machine = expect_hex_u64(receipt.take("machine")?, "receipt.machine")?;

    decode_receipt_axes(receipt.take("axes")?)?;

    decode_receipt_spec(receipt.take("spec")?)?;

    let measurement = decode_receipt_measurement(receipt.take("measurement")?)?;

    let metrics = decode_receipt_metrics(&mut receipt, &measurement)?;
    receipt.finish()?;
    let size = validate_rederived_metrics(&measurement, metrics.rate_bits)?;
    Ok(ReceiptObservation {
        kernel,
        version,
        machine,
        observations: measurement
            .sample_times
            .into_iter()
            .map(|cost_s| CostObservation { size, cost_s })
            .collect(),
        reps: metrics.reps,
    })
}

fn sample_bits_len(reps: u64) -> Result<usize, TuneModelError> {
    usize::try_from(reps)
        .map_err(|_| invalid_receipt("receipt.reps", "repetition count exceeds usize"))
}

fn nearest_rank(sorted: &[f64], probability: f64) -> Result<f64, TuneModelError> {
    let last = sorted
        .len()
        .checked_sub(1)
        .ok_or_else(|| invalid_receipt("receipt.measurement", "empty timed sample"))?;
    let index = ((sorted.len() as f64 - 1.0) * probability).round() as usize;
    sorted
        .get(index.min(last))
        .copied()
        .ok_or_else(|| invalid_receipt("receipt.measurement", "quantile index out of bounds"))
}

fn validate_result_manifest(
    value: JsonValue,
    expected_kernel: &str,
    expected_version: &str,
    expected_payload: &str,
    expected_count: u64,
) -> Result<(), TuneModelError> {
    let mut manifest = ObjectFields::new(value, "op.ir.result_manifest")?;
    let schema = expect_string(manifest.take("schema")?, "op.ir.result_manifest.schema")?;
    if schema != "fs-roofline-run-manifest-v1" {
        return Err(invalid_receipt(
            "op.ir.result_manifest.schema",
            "unsupported result manifest schema",
        ));
    }
    let entries = expect_array(manifest.take("entries")?, "op.ir.result_manifest.entries")?;
    manifest.finish()?;
    if u64::try_from(entries.len()).ok() != Some(expected_count) || entries.is_empty() {
        return Err(invalid_receipt(
            "op.ir.result_manifest.entries",
            "entry count does not match op.ir.kernels",
        ));
    }
    let mut matching_entries = 0_usize;
    for (expected_ordinal, entry) in entries.into_iter().enumerate() {
        let mut entry = ObjectFields::new(
            entry,
            format!("op.ir.result_manifest.entries[{expected_ordinal}]"),
        )?;
        let ordinal = expect_u64(
            entry.take("ordinal")?,
            &format!("op.ir.result_manifest.entries[{expected_ordinal}].ordinal"),
        )?;
        if usize::try_from(ordinal).ok() != Some(expected_ordinal) {
            return Err(invalid_receipt(
                format!("op.ir.result_manifest.entries[{expected_ordinal}].ordinal"),
                "ordinals must be exactly 0..entry_count in order",
            ));
        }
        let kernel = expect_string(
            entry.take("kernel")?,
            &format!("op.ir.result_manifest.entries[{expected_ordinal}].kernel"),
        )?;
        let version = expect_string(
            entry.take("version")?,
            &format!("op.ir.result_manifest.entries[{expected_ordinal}].version"),
        )?;
        let payload = expect_hash(
            entry.take("payload")?,
            &format!("op.ir.result_manifest.entries[{expected_ordinal}].payload"),
        )?;
        entry.finish()?;
        if kernel == expected_kernel && version == expected_version && payload == expected_payload {
            matching_entries += 1;
        }
    }
    if matching_entries != 1 {
        return Err(invalid_receipt(
            "op.ir.result_manifest.entries",
            "exactly one manifest member must bind this kernel/version/payload",
        ));
    }
    Ok(())
}

fn validate_versions(versions: &str, build_identity: &str) -> Result<(), TuneModelError> {
    let mut versions = ObjectFields::new(
        StrictJson::parse(versions, fs_ledger::MAX_TUNE_PARAMS_BYTES)?,
        "op.versions",
    )?;
    let executable = expect_hash(
        versions.take("frankensim_executable")?,
        "op.versions.frankensim_executable",
    )?;
    if executable != build_identity {
        return Err(TuneModelError::ScopeMismatch {
            field: "build_identity",
        });
    }
    let roofline_version = expect_string(versions.take("fs-roofline")?, "op.versions.fs-roofline")?;
    if roofline_version != env!("CARGO_PKG_VERSION") {
        return Err(invalid_receipt(
            "op.versions.fs-roofline",
            format!("unsupported producer version {roofline_version:?}"),
        ));
    }
    versions.finish()
}

fn validate_op_ir(
    ir: &str,
    binding: &RowBinding,
    receipt: &ReceiptObservation,
) -> Result<(), TuneModelError> {
    let mut ir = ObjectFields::new(
        StrictJson::parse(ir, fs_ledger::MAX_TUNE_PARAMS_BYTES)?,
        "op.ir",
    )?;
    if expect_string(ir.take("op")?, "op.ir.op")? != "perf.roofline" {
        return Err(invalid_receipt("op.ir.op", "expected perf.roofline"));
    }
    let kernel_count = expect_u64(ir.take("kernels")?, "op.ir.kernels")?;
    if kernel_count == 0 {
        return Err(invalid_receipt("op.ir.kernels", "must be positive"));
    }
    let fingerprint = expect_hex_u64(ir.take("fingerprint")?, "op.ir.fingerprint")?;
    let _ = expect_hex_u64(ir.take("post_fingerprint")?, "op.ir.post_fingerprint")?;
    if fingerprint != receipt.machine {
        return Err(TuneModelError::ScopeMismatch {
            field: "op fingerprint",
        });
    }
    if !expect_bool(
        ir.take("measurement_admitted")?,
        "op.ir.measurement_admitted",
    )? || !expect_bool(ir.take("admitted")?, "op.ir.admitted")?
    {
        return Err(invalid_receipt(
            "op.ir.admission",
            "production tune evidence requires admitted=true and measurement_admitted=true",
        ));
    }
    if expect_string(ir.take("protocol")?, "op.ir.protocol")? != "production-v2" {
        return Err(invalid_receipt("op.ir.protocol", "expected production-v2"));
    }
    let _ = expect_hash(ir.take("run_nonce")?, "op.ir.run_nonce")?;
    let _ = expect_hash(ir.take("pre_axes_receipt")?, "op.ir.pre_axes_receipt")?;
    let _ = expect_hash(ir.take("post_axes_receipt")?, "op.ir.post_axes_receipt")?;
    if expect_string(
        ir.take("dependency_graph_evidence")?,
        "op.ir.dependency_graph_evidence",
    )? != "operator-observed-receipt"
    {
        return Err(invalid_receipt(
            "op.ir.dependency_graph_evidence",
            "expected operator-observed-receipt",
        ));
    }
    if expect_hash(
        ir.take("dependency_receipt_digest")?,
        "op.ir.dependency_receipt_digest",
    )? != binding.dependency_receipt_digest
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "dependency_receipt_digest",
        });
    }
    if expect_hash(
        ir.take("dependency_receipt_artifact")?,
        "op.ir.dependency_receipt_artifact",
    )? != binding.dependency_receipt_artifact
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "dependency_receipt_artifact",
        });
    }
    if expect_hash(
        ir.take("finalized_run_receipt")?,
        "op.ir.finalized_run_receipt",
    )? != binding.run_receipt
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "finalized_run_receipt",
        });
    }
    validate_result_manifest(
        ir.take("result_manifest")?,
        &receipt.kernel,
        &receipt.version,
        &binding.payload_artifact,
        kernel_count,
    )?;
    if !matches!(ir.take("baseline_admission")?, JsonValue::Object(_)) {
        return Err(invalid_receipt(
            "op.ir.baseline_admission",
            "expected object",
        ));
    }
    ir.finish()
}

fn content_hash(text: &str, field: &str) -> Result<fs_ledger::ContentHash, TuneModelError> {
    fs_ledger::ContentHash::from_hex(text)
        .ok_or_else(|| invalid_receipt(field, "invalid content hash"))
}

fn validate_provenance(
    ledger: &Ledger,
    row: &fs_ledger::TuneRow,
    binding: &RowBinding,
    receipt: &ReceiptObservation,
    machine: &[u8],
) -> Result<(), TuneModelError> {
    let baseline = content_hash(&binding.baseline_hash, "params.baseline_hash")?;
    if machine.get(8..) != Some(baseline.as_bytes()) {
        return Err(TuneModelError::ScopeMismatch {
            field: "baseline_hash",
        });
    }
    let op_id = i64::try_from(binding.op)
        .map_err(|_| invalid_receipt("params.op", "operation id exceeds i64"))?;
    let op = ledger
        .op(op_id)?
        .ok_or(TuneModelError::ScopeMismatch { field: "operation" })?;
    if op.id != op_id
        || op.session.as_deref() != Some(b"roofline".as_slice())
        || op.seed != b"roofline"
        || op.budget != "{\"wall_s\":600}"
        || op.capability != "{\"ops\":[\"perf.roofline\"]}"
        || op.outcome.as_deref() != Some("ok")
        || op.diag.is_some()
        || op.t_end.is_none_or(|end| end < op.t_start)
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "operation envelope",
        });
    }
    validate_versions(&op.versions, &binding.build_identity)?;
    validate_op_ir(&op.ir, binding, receipt)?;

    let payload_hash = content_hash(&binding.payload_artifact, "params.payload_artifact")?;
    let payload_info =
        ledger
            .artifact_info(&payload_hash)?
            .ok_or(TuneModelError::ScopeMismatch {
                field: "payload_artifact",
            })?;
    if payload_info.kind != "roofline-benchmark-result"
        || payload_info.meta.as_deref() != Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}")
        || u64::try_from(row.measured.len()).ok() != Some(payload_info.len)
        || ledger.get_artifact(&payload_hash)?.as_deref() != Some(row.measured.as_bytes())
        || !ledger.edge_exists(op_id, &payload_hash, fs_ledger::EdgeRole::Out)?
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "payload_artifact",
        });
    }

    let dependency_hash = content_hash(
        &binding.dependency_receipt_artifact,
        "params.dependency_receipt_artifact",
    )?;
    let dependency_info =
        ledger
            .artifact_info(&dependency_hash)?
            .ok_or(TuneModelError::ScopeMismatch {
                field: "dependency_receipt_artifact",
            })?;
    let dependency_bytes =
        ledger
            .get_artifact(&dependency_hash)?
            .ok_or(TuneModelError::ScopeMismatch {
                field: "dependency_receipt_artifact",
            })?;
    if dependency_info.kind != "fs-la-depgraph-receipt"
        || dependency_info.meta.as_deref()
            != Some("{\"schema\":\"fs-la-depgraph-receipt-v1\",\"trust\":\"operator-observed\"}")
        || !ledger.edge_exists(op_id, &dependency_hash, fs_ledger::EdgeRole::In)?
        || fs_blake3::hash_domain(
            "org.frankensim.fs-la.depgraph-receipt.v1",
            &dependency_bytes,
        )
        .to_string()
            != binding.dependency_receipt_digest
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "dependency_receipt",
        });
    }
    Ok(())
}

/// Rebuild one model from one exact production roofline tune key.
///
/// This API deliberately uses `tune_get`, not a per-kernel scan: foreign
/// machines and neighboring shape classes can never contribute. A current
/// production row contributes each bounded timed repetition as a same-size
/// observation. Rows with fewer than [`crate::cost::MIN_OBS`] repetitions
/// therefore continue to refuse prediction honestly.
///
/// # Errors
/// Ledger corruption/absence, malformed or duplicate-key JSON, unsupported
/// schema versions, scope mismatches, nonfinite measurements, and model
/// refusals are all returned explicitly.
pub fn cost_model_from_tune(
    ledger: &Ledger,
    kernel: &str,
    shape_class: &str,
    machine: &[u8],
) -> Result<CostModel, TuneModelError> {
    if machine.len() != ROOFLINE_MACHINE_KEY_BYTES {
        return Err(TuneModelError::ScopeMismatch { field: "machine" });
    }
    let row = ledger
        .tune_get(kernel, shape_class, machine)?
        .ok_or(TuneModelError::MissingRow)?;
    let binding = decode_row_binding(&row.params)?;
    let receipt = decode_receipt(&row.measured)?;
    if receipt.kernel != kernel || row.kernel != kernel {
        return Err(TuneModelError::ScopeMismatch { field: "kernel" });
    }
    let expected_shape = format!(
        "{ROOFLINE_TUNE_SHAPE_PREFIX}:{}:run={}:op={}",
        receipt.version, binding.run_receipt, binding.op
    );
    if row.shape_class != shape_class || expected_shape != shape_class {
        return Err(TuneModelError::ScopeMismatch {
            field: "shape_class",
        });
    }
    let fingerprint_bytes = machine
        .get(..8)
        .and_then(|bytes| bytes.try_into().ok())
        .ok_or(TuneModelError::ScopeMismatch { field: "machine" })?;
    let fingerprint = u64::from_le_bytes(fingerprint_bytes);
    if row.machine != machine || receipt.machine != fingerprint {
        return Err(TuneModelError::ScopeMismatch { field: "machine" });
    }
    if binding.reps != receipt.reps {
        return Err(TuneModelError::ScopeMismatch { field: "reps" });
    }
    if binding.payload_artifact != fs_ledger::hash_bytes(row.measured.as_bytes()).to_string() {
        return Err(TuneModelError::ScopeMismatch {
            field: "payload_artifact",
        });
    }
    validate_provenance(ledger, &row, &binding, &receipt, machine)?;
    CostModel::fit(&receipt.observations).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn oracle_spec(name: &str) -> fs_geom::ConverterSpec {
        fs_geom::ConverterSpec {
            from: "source".to_string(),
            to: "target".to_string(),
            name: name.to_string(),
            base_cost_s: 1.0,
            error: fs_geom::ErrorModel::AdditiveAbs(0.1),
            certified: false,
        }
    }

    fn bits(value: f64) -> String {
        format!("{:016x}", value.to_bits())
    }

    fn production_receipt() -> String {
        let elements = 1_000_u64;
        let sample = [0.003_f64, 0.001, 0.002];
        let median = 0.002_f64;
        let p25 = 0.002_f64;
        let p75 = 0.003_f64;
        let dispersion = (p75 - p25) / median;
        let rate = elements as f64 / median;
        let sample_bits = sample
            .iter()
            .map(|seconds| format!("\"{}\"", bits(*seconds)))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{{\"receipt_version\":3,\"kernel\":\"simd-axpy-f64\",\"version\":\"1\",\"machine\":\"0102030405060708\",\
             \"axes\":{{\"logical_cpus\":8,\"bandwidth_single_bits\":\"{}\",\"bandwidth_all_core_bits\":\"{}\",\"peak_single_bits\":\"{}\",\"peak_all_core_bits\":\"{}\"}},\
             \"spec\":{{\"bytes_per_elem_bits\":\"{}\",\"flops_per_elem_bits\":\"{}\",\"threading\":\"single_thread\",\"target_axis\":\"binding_roof\",\"target_fraction_bits\":\"{}\"}},\
             \"measurement\":{{\"origin\":\"timed\",\"elements\":{elements},\"warmup_runs\":1,\"sample_seconds_bits\":[{sample_bits}],\"decision_binding_hashes\":[null,null,null],\"median_seconds_bits\":\"{}\",\"p25_seconds_bits\":\"{}\",\"p75_seconds_bits\":\"{}\",\"dispersion_bits\":\"{}\"}},\
             \"execution\":null,\"elems_per_sec_bits\":\"{}\",\"gbs_bits\":\"{}\",\"gflops_bits\":\"{}\",\"limit_elems_per_sec_bits\":\"{}\",\"attainment_bits\":\"{}\",\"target_attainment_bits\":\"{}\",\"dispersion_bits\":\"{}\",\
             \"elems_per_sec_display\":{rate},\"gbs\":0.004,\"gflops\":0.001,\"limit_elems_per_sec\":1000000,\"roof\":\"bandwidth\",\"attainment\":0.5,\"target_attainment\":1.0,\"dispersion\":{dispersion},\"reps\":3,\"verdict\":\"within_band\",\"invalid_reason\":null}}",
            bits(100.0),
            bits(200.0),
            bits(1_000.0),
            bits(2_000.0),
            bits(8.0),
            bits(2.0),
            bits(0.5),
            bits(median),
            bits(p25),
            bits(p75),
            bits(dispersion),
            bits(rate),
            bits(0.004),
            bits(0.001),
            bits(1_000_000.0),
            bits(0.5),
            bits(1.0),
            bits(dispersion),
        )
    }

    #[test]
    fn receipt_v3_decodes_every_timed_sample() {
        let decoded = decode_receipt(&production_receipt()).unwrap();
        assert_eq!(decoded.observations.len(), 3);
        assert_eq!(
            decoded.observations[0].size.to_bits(),
            1_000.0_f64.to_bits()
        );
        assert_eq!(
            decoded.observations[0].cost_s.to_bits(),
            0.003_f64.to_bits()
        );
        let model = CostModel::fit(&decoded.observations).unwrap();
        assert_eq!(model.n_obs(), 3);
        assert!(model.predict(1_000.0).is_ok());
    }

    #[test]
    fn receipt_parser_rejects_duplicate_keys_and_forged_statistics() {
        let receipt = production_receipt();
        let duplicate = receipt.replacen(
            "{\"receipt_version\":3,",
            "{\"receipt_version\":3,\"kernel\":\"shadow\",",
            1,
        );
        assert!(matches!(
            decode_receipt(&duplicate),
            Err(TuneModelError::InvalidReceipt { problem, .. })
                if problem.contains("duplicate key")
        ));

        let median = bits(0.002);
        let forged = receipt.replacen(
            &format!("\"median_seconds_bits\":\"{median}\""),
            &format!("\"median_seconds_bits\":\"{}\"", bits(0.0025)),
            1,
        );
        assert!(matches!(
            decode_receipt(&forged),
            Err(TuneModelError::InvalidReceipt { field, .. })
                if field == "receipt.measurement.statistics"
        ));
    }

    #[test]
    fn strict_json_accepts_byte_cap_and_refuses_limit_plus_one() {
        let at_cap = format!("null{}", " ".repeat(MAX_ROOFLINE_RECEIPT_BYTES - 4));
        assert!(StrictJson::parse(&at_cap, MAX_ROOFLINE_RECEIPT_BYTES).is_ok());
        let over_cap = format!("{at_cap} ");
        assert!(matches!(
            StrictJson::parse(&over_cap, MAX_ROOFLINE_RECEIPT_BYTES),
            Err(TuneModelError::InvalidReceipt { problem, .. })
                if problem.contains("exceeds limit")
        ));
    }

    #[test]
    fn planner_oracle_registration_is_spec_scoped() {
        let mut oracle = PlanCostOracle::new();
        let spec = oracle_spec("frep->sdf");
        oracle.register_edge(&spec, 100.0).unwrap();
        assert_eq!(oracle.model("frep->sdf").unwrap().n_obs(), 0);

        let mut changed = spec.clone();
        changed.error = fs_geom::ErrorModel::AdditiveAbs(0.2);
        assert_eq!(
            oracle.register_edge(&changed, 100.0),
            Err(PlanOracleError::SpecificationConflict)
        );
        oracle.register_edge(&spec, 200.0).unwrap();
        assert_eq!(
            oracle.register_edge(&spec, f64::NAN),
            Err(PlanOracleError::InvalidReferenceSize)
        );
    }

    #[test]
    fn planner_oracle_accepts_edge_cap_and_refuses_limit_plus_one() {
        let mut oracle = PlanCostOracle::new();
        for index in 0..MAX_PLAN_ORACLE_EDGES {
            oracle
                .register_edge(&oracle_spec(&format!("edge-{index}")), 1.0)
                .unwrap();
        }
        assert_eq!(
            oracle.register_edge(&oracle_spec("one-too-many"), 1.0),
            Err(PlanOracleError::EdgeLimit {
                limit: MAX_PLAN_ORACLE_EDGES
            })
        );
    }
}
