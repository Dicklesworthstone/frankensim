//! Consumer wiring (plan §11.4 "consumers wired"): fitted cost models
//! behind fs-geom's `CostOracle` so the Rep Router plans with THIS
//! machine's measured history, and a fs-ledger `tune`-table loader so
//! models are rebuilt deterministically from ledger snapshots.

use std::collections::BTreeMap;

use crate::cost::{CostModel, CostObservation, CostRefusal, MAX_COST_OBSERVATIONS};
use crate::sealed::{CostModelScope, SealedCostModel};
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
pub const ROOFLINE_TUNE_SHAPE_PREFIX: &str = "roofline-v8";

/// Exact byte width of fs-roofline's fingerprint + baseline machine key.
pub const ROOFLINE_MACHINE_KEY_BYTES: usize = 40;

/// Maximum receipt bytes decoded by the planner. This matches the ledger's
/// tune measurement admission bound, so parsing is bounded before allocation.
pub const MAX_ROOFLINE_RECEIPT_BYTES: usize = fs_ledger::MAX_TUNE_MEASURED_BYTES;

/// Exact payload ceiling owned by fs-la's dependency-receipt grammar.
/// A source-pin test fails if the producer declaration drifts without this
/// consumer being updated; no L6-to-L1 runtime edge is added solely for a cap.
pub const MAX_DEPGRAPH_RECEIPT_BYTES: u64 = 1_048_576;

const MAX_RECEIPT_JSON_DEPTH: usize = 32;
const MAX_RECEIPT_JSON_NODES: usize = 65_536;
const MAX_RECEIPT_JSON_STRING_BYTES: usize = 256 * 1024;
const MAX_RECEIPT_JSON_CONTAINER_ITEMS: usize = 32_768;
const MAX_RECEIPT_JSON_NUMBER_BYTES: usize = 64;
const MAX_PRODUCTION_ELEMENTS: u64 = 1 << 24;
const MAX_PRODUCTION_KERNEL_RUNS: u64 = 64;
const MAX_PRODUCTION_WARMUP: u64 = MAX_PRODUCTION_KERNEL_RUNS - 1;
const MAX_PRODUCTION_REPS: u64 = MAX_PRODUCTION_KERNEL_RUNS;
const MAX_PRODUCTION_REGISTRY_FLOPS: u128 = 1 << 39;
const MAX_PRODUCTION_REGISTRY_BYTES: u128 = 1 << 33;
const MAX_PRODUCTION_LOGICAL_CPUS: u64 = 4_096;
const SEALED_PRODUCTION_REGISTRY: [(&str, &str); 4] = [
    ("simd-axpy-f64", "1"),
    ("simd-dot-f64", "1"),
    ("simd-sum-f64", "1"),
    ("gemm-f64", "2"),
];
const MAX_BASELINE_STRING_BYTES: usize = 4_096;
const MAX_BASELINE_LINE_BYTES: usize = 16 * 1_024;
const MAX_BASELINE_AGE_DAYS: u64 = 365;
const MIN_PROMOTION_RUNS: usize = 3;
const WALL_NS_PER_DAY: u64 = 86_400 * 1_000_000_000;
const BASELINE_HASH_DOMAIN: &str = "frankensim.fs-roofline.baseline.v1";
const PRODUCTION_AXES_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-roofline.production-axes-receipt.v1";
const RESULT_MANIFEST_DOMAIN: &str = "org.frankensim.fs-roofline.run-result-manifest.v1";
const FINALIZED_RUN_DOMAIN: &str = "org.frankensim.fs-roofline.finalized-run.v3";

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

#[derive(Debug, Clone)]
struct ReceiptObservation {
    kernel: String,
    version: String,
    machine: u64,
    logical_cpus: u64,
    axis_bits: [u64; 4],
    elements: u64,
    warmup_runs: u64,
    observations: Vec<CostObservation>,
    reps: u64,
}

#[derive(Debug, Clone)]
struct RowBinding {
    op: u64,
    run_receipt: String,
    payload_artifact: String,
    dependency_receipt_artifact: String,
    dependency_receipt_digest: String,
    baseline_hash: String,
    build_identity: String,
    reps: u64,
    post_axis_bits: [u64; 4],
}

#[derive(Debug, Clone)]
struct ManifestEntry {
    ordinal: u64,
    kernel: String,
    version: String,
    payload: String,
}

#[derive(Debug)]
struct ValidatedProductionOp {
    baseline_admission: String,
    result_manifest: String,
    manifest: Vec<ManifestEntry>,
    finalized_run_receipt: String,
    admission: ValidatedAxisAdmission,
}

#[derive(Debug)]
struct ValidatedAxisAdmission {
    decision_day: u64,
    baseline_hash: String,
    pre: DecodedAxes,
    post: DecodedAxes,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct DecodedIdentity {
    fingerprint: u64,
    cpu_brand: String,
    logical_cpus: u64,
    os: String,
    arch: String,
    firmware: String,
    canonical: String,
}

#[derive(Debug)]
struct DecodedAxes {
    fingerprint: u64,
    cpu_brand: String,
    logical_cpus: u64,
    bits: [u64; 4],
    canonical: String,
}

#[derive(Debug)]
struct DecodedBaseline {
    identity: DecodedIdentity,
    bits: [u64; 4],
    source_receipts: Vec<String>,
    promoted_day: u64,
    age_policy_days: u64,
    canonical: String,
}

struct DecodedMeasurement {
    elements: u64,
    warmup_runs: u64,
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

fn validate_production_run_counts(warmup_runs: u64, reps: u64) -> Result<u64, TuneModelError> {
    if warmup_runs > MAX_PRODUCTION_WARMUP {
        return Err(invalid_receipt(
            "receipt.measurement.warmup_runs",
            format!("must be in 0..={MAX_PRODUCTION_WARMUP}"),
        ));
    }
    if reps == 0 || reps > MAX_PRODUCTION_REPS {
        return Err(invalid_receipt(
            "receipt.reps",
            format!("must be in 1..={MAX_PRODUCTION_REPS}"),
        ));
    }
    let runs_per_kernel = warmup_runs.checked_add(reps).ok_or_else(|| {
        invalid_receipt(
            "receipt.measurement.warmup_runs",
            "warmup + repetition count overflows u64",
        )
    })?;
    if runs_per_kernel > MAX_PRODUCTION_KERNEL_RUNS {
        return Err(invalid_receipt(
            "receipt.measurement.warmup_runs",
            format!(
                "warmup + reps must be at most {MAX_PRODUCTION_KERNEL_RUNS}, got {runs_per_kernel}"
            ),
        ));
    }
    Ok(runs_per_kernel)
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
    if reps == 0 || reps > MAX_PRODUCTION_REPS {
        return Err(invalid_receipt(
            "params.reps",
            format!("must be in 1..={MAX_PRODUCTION_REPS}"),
        ));
    }
    let mut post_axis_bits = [0_u64; 4];
    for (slot, field) in post_axis_bits.iter_mut().zip([
        "post_bandwidth_single_bits",
        "post_bandwidth_all_core_bits",
        "post_peak_single_bits",
        "post_peak_all_core_bits",
    ]) {
        *slot = finite_from_bits(object.take(field)?, &format!("params.{field}"), true)?.0;
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
        post_axis_bits,
    })
}

fn decode_receipt_axes(value: JsonValue) -> Result<(u64, [u64; 4]), TuneModelError> {
    let mut axes = ObjectFields::new(value, "receipt.axes")?;
    let logical_cpus = expect_u64(axes.take("logical_cpus")?, "receipt.axes.logical_cpus")?;
    if logical_cpus == 0 || logical_cpus > MAX_PRODUCTION_LOGICAL_CPUS {
        return Err(invalid_receipt(
            "receipt.axes.logical_cpus",
            format!("must be in 1..={MAX_PRODUCTION_LOGICAL_CPUS}"),
        ));
    }
    let mut axis_bits = [0_u64; 4];
    for (slot, field) in axis_bits.iter_mut().zip([
        "bandwidth_single_bits",
        "bandwidth_all_core_bits",
        "peak_single_bits",
        "peak_all_core_bits",
    ]) {
        *slot = finite_from_bits(axes.take(field)?, &format!("receipt.axes.{field}"), true)?.0;
    }
    axes.finish()?;
    Ok((logical_cpus, axis_bits))
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
    if elements == 0 || elements > MAX_PRODUCTION_ELEMENTS {
        return Err(invalid_receipt(
            "receipt.measurement.elements",
            format!("must be in 1..={MAX_PRODUCTION_ELEMENTS}"),
        ));
    }
    let warmup_runs = expect_u64(
        measurement.take("warmup_runs")?,
        "receipt.measurement.warmup_runs",
    )?;
    if warmup_runs > MAX_PRODUCTION_WARMUP {
        return Err(invalid_receipt(
            "receipt.measurement.warmup_runs",
            format!("must be in 0..={MAX_PRODUCTION_WARMUP}"),
        ));
    }
    let sample_bits = expect_array(
        measurement.take("sample_seconds_bits")?,
        "receipt.measurement.sample_seconds_bits",
    )?;
    if sample_bits.is_empty()
        || u64::try_from(sample_bits.len())
            .ok()
            .is_none_or(|len| len > MAX_PRODUCTION_REPS)
    {
        return Err(invalid_receipt(
            "receipt.measurement.sample_seconds_bits",
            format!("length must be in 1..={MAX_PRODUCTION_REPS}"),
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
    if u64::try_from(decision_hashes.len())
        .ok()
        .is_none_or(|len| len > MAX_PRODUCTION_REPS)
    {
        return Err(invalid_receipt(
            "receipt.measurement.decision_binding_hashes",
            format!("length exceeds limit {MAX_PRODUCTION_REPS}"),
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
        warmup_runs,
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
    if reps == 0 || reps > MAX_PRODUCTION_REPS {
        return Err(invalid_receipt(
            "receipt.reps",
            format!("must be in 1..={MAX_PRODUCTION_REPS}"),
        ));
    }
    validate_production_run_counts(measurement.warmup_runs, reps)?;
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

    let (logical_cpus, axis_bits) = decode_receipt_axes(receipt.take("axes")?)?;

    decode_receipt_spec(receipt.take("spec")?)?;

    let measurement = decode_receipt_measurement(receipt.take("measurement")?)?;

    let metrics = decode_receipt_metrics(&mut receipt, &measurement)?;
    receipt.finish()?;
    let size = validate_rederived_metrics(&measurement, metrics.rate_bits)?;
    Ok(ReceiptObservation {
        kernel,
        version,
        machine,
        logical_cpus,
        axis_bits,
        elements: measurement.elements,
        warmup_runs: measurement.warmup_runs,
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

fn canonical_json_string(value: &str) -> String {
    use std::fmt::Write as _;

    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if control.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

fn validate_baseline_text(value: &str, field: &str) -> Result<(), TuneModelError> {
    if value.trim().is_empty()
        || value.len() > MAX_BASELINE_STRING_BYTES
        || value.chars().any(char::is_control)
    {
        return Err(invalid_receipt(
            field,
            "expected nonblank, control-free text within the producer bound",
        ));
    }
    Ok(())
}

fn validate_authority_text(value: &str, field: &str) -> Result<(), TuneModelError> {
    if value.is_empty()
        || value.trim() != value
        || value.len() > MAX_BASELINE_STRING_BYTES
        || value.chars().any(char::is_control)
    {
        return Err(invalid_receipt(
            field,
            "expected trimmed, nonblank, control-free text within the authority bound",
        ));
    }
    Ok(())
}

fn decode_hash_array(value: JsonValue, field: &str) -> Result<Vec<String>, TuneModelError> {
    let values = expect_array(value, field)?;
    if values.len() < MIN_PROMOTION_RUNS || values.len() > MAX_RECEIPT_JSON_CONTAINER_ITEMS {
        return Err(invalid_receipt(
            field,
            format!("expected {MIN_PROMOTION_RUNS}..={MAX_RECEIPT_JSON_CONTAINER_ITEMS} receipts"),
        ));
    }
    let mut hashes = Vec::with_capacity(values.len());
    for (index, value) in values.into_iter().enumerate() {
        hashes.push(expect_hash(value, &format!("{field}[{index}]"))?);
    }
    if hashes.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(invalid_receipt(field, "hashes must be sorted and unique"));
    }
    Ok(hashes)
}

fn hash_array_json(hashes: &[String]) -> String {
    format!(
        "[{}]",
        hashes
            .iter()
            .map(|hash| format!("\"{hash}\""))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn identity_json(identity: &DecodedIdentity) -> String {
    format!(
        "{{\"fingerprint\":\"{:016x}\",\"cpu_brand\":{},\"logical_cpus\":{},\"os\":{},\"arch\":{},\"firmware\":{}}}",
        identity.fingerprint,
        canonical_json_string(&identity.cpu_brand),
        identity.logical_cpus,
        canonical_json_string(&identity.os),
        canonical_json_string(&identity.arch),
        canonical_json_string(&identity.firmware),
    )
}

fn decode_identity(value: JsonValue, what: &str) -> Result<DecodedIdentity, TuneModelError> {
    let mut identity = ObjectFields::new(value, what)?;
    let fingerprint = expect_hex_u64(
        identity.take("fingerprint")?,
        &format!("{what}.fingerprint"),
    )?;
    let cpu_brand = expect_string(identity.take("cpu_brand")?, &format!("{what}.cpu_brand"))?;
    let logical_cpus = expect_u64(
        identity.take("logical_cpus")?,
        &format!("{what}.logical_cpus"),
    )?;
    let os = expect_string(identity.take("os")?, &format!("{what}.os"))?;
    let arch = expect_string(identity.take("arch")?, &format!("{what}.arch"))?;
    let firmware = expect_string(identity.take("firmware")?, &format!("{what}.firmware"))?;
    identity.finish()?;
    if logical_cpus == 0 || u32::try_from(logical_cpus).is_err() {
        return Err(invalid_receipt(
            format!("{what}.logical_cpus"),
            "must be in 1..=u32::MAX",
        ));
    }
    for (field, value) in [
        ("cpu_brand", cpu_brand.as_str()),
        ("os", os.as_str()),
        ("arch", arch.as_str()),
        ("firmware", firmware.as_str()),
    ] {
        validate_baseline_text(value, &format!("{what}.{field}"))?;
    }
    let mut decoded = DecodedIdentity {
        fingerprint,
        cpu_brand,
        logical_cpus,
        os,
        arch,
        firmware,
        canonical: String::new(),
    };
    decoded.canonical = identity_json(&decoded);
    Ok(decoded)
}

fn axes_values(bits: [u64; 4]) -> [f64; 4] {
    bits.map(f64::from_bits)
}

fn axes_are_plausible(bits: [u64; 4]) -> bool {
    let values = axes_values(bits);
    values.iter().all(|value| value.is_finite() && *value > 0.0)
        && values[0] >= 5.0
        && values[2] >= 5.0
        && values[1] >= values[0] * 0.5
        && values[3] >= values[2] * 0.5
}

fn decode_axes(value: JsonValue, what: &str) -> Result<DecodedAxes, TuneModelError> {
    let mut axes = ObjectFields::new(value, what)?;
    let fingerprint = expect_hex_u64(axes.take("fingerprint")?, &format!("{what}.fingerprint"))?;
    let cpu_brand = expect_string(axes.take("cpu_brand")?, &format!("{what}.cpu_brand"))?;
    validate_baseline_text(&cpu_brand, &format!("{what}.cpu_brand"))?;
    let logical_cpus = expect_u64(axes.take("logical_cpus")?, &format!("{what}.logical_cpus"))?;
    if logical_cpus == 0 || u32::try_from(logical_cpus).is_err() {
        return Err(invalid_receipt(
            format!("{what}.logical_cpus"),
            "must be in 1..=u32::MAX",
        ));
    }
    let mut bits = [0_u64; 4];
    for (slot, field) in bits.iter_mut().zip([
        "bandwidth_single_bits",
        "bandwidth_all_core_bits",
        "peak_single_bits",
        "peak_all_core_bits",
    ]) {
        *slot = finite_from_bits(axes.take(field)?, &format!("{what}.{field}"), true)?.0;
    }
    axes.finish()?;
    if !axes_are_plausible(bits) {
        return Err(invalid_receipt(
            what,
            "axes fail producer plausibility floors",
        ));
    }
    let canonical = format!(
        "{{\"fingerprint\":\"{fingerprint:016x}\",\"cpu_brand\":{},\"logical_cpus\":{logical_cpus},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}}",
        canonical_json_string(&cpu_brand),
        bits[0],
        bits[1],
        bits[2],
        bits[3],
    );
    Ok(DecodedAxes {
        fingerprint,
        cpu_brand,
        logical_cpus,
        bits,
        canonical,
    })
}

fn decode_baseline(value: JsonValue) -> Result<DecodedBaseline, TuneModelError> {
    let what = "op.ir.baseline_admission.baseline";
    let mut baseline = ObjectFields::new(value, what)?;
    if expect_u64(
        baseline.take("schema_version")?,
        &format!("{what}.schema_version"),
    )? != 1
    {
        return Err(invalid_receipt(
            format!("{what}.schema_version"),
            "expected baseline schema 1",
        ));
    }
    for (field, expected) in [
        ("low_band_bits", 0.70_f64.to_bits()),
        ("high_band_bits", 1.15_f64.to_bits()),
        ("promotion_drift_bits", 0.25_f64.to_bits()),
    ] {
        if expect_hex_u64(baseline.take(field)?, &format!("{what}.{field}"))? != expected {
            return Err(invalid_receipt(
                format!("{what}.{field}"),
                "producer policy constant mismatch",
            ));
        }
    }
    let fingerprint = expect_hex_u64(
        baseline.take("fingerprint")?,
        &format!("{what}.fingerprint"),
    )?;
    let cpu_brand = expect_string(baseline.take("cpu_brand")?, &format!("{what}.cpu_brand"))?;
    let logical_cpus = expect_u64(
        baseline.take("logical_cpus")?,
        &format!("{what}.logical_cpus"),
    )?;
    let os = expect_string(baseline.take("os")?, &format!("{what}.os"))?;
    let arch = expect_string(baseline.take("arch")?, &format!("{what}.arch"))?;
    let firmware = expect_string(baseline.take("firmware")?, &format!("{what}.firmware"))?;
    if logical_cpus == 0 || u32::try_from(logical_cpus).is_err() {
        return Err(invalid_receipt(
            format!("{what}.logical_cpus"),
            "must be in 1..=u32::MAX",
        ));
    }
    for (field, value) in [
        ("cpu_brand", cpu_brand.as_str()),
        ("os", os.as_str()),
        ("arch", arch.as_str()),
        ("firmware", firmware.as_str()),
    ] {
        validate_baseline_text(value, &format!("{what}.{field}"))?;
    }
    let mut bits = [0_u64; 4];
    for (slot, field) in bits.iter_mut().zip([
        "bandwidth_single_bits",
        "bandwidth_all_core_bits",
        "peak_single_bits",
        "peak_all_core_bits",
    ]) {
        *slot = finite_from_bits(baseline.take(field)?, &format!("{what}.{field}"), true)?.0;
    }
    if !axes_are_plausible(bits) {
        return Err(invalid_receipt(
            what,
            "baseline axes fail producer plausibility floors",
        ));
    }
    let source_receipts = decode_hash_array(
        baseline.take("source_receipts")?,
        &format!("{what}.source_receipts"),
    )?;
    let promoted_by = expect_string(
        baseline.take("promoted_by")?,
        &format!("{what}.promoted_by"),
    )?;
    let justification = expect_string(
        baseline.take("justification")?,
        &format!("{what}.justification"),
    )?;
    validate_baseline_text(&promoted_by, &format!("{what}.promoted_by"))?;
    validate_baseline_text(&justification, &format!("{what}.justification"))?;
    let promoted_day = expect_u64(
        baseline.take("promoted_day")?,
        &format!("{what}.promoted_day"),
    )?;
    let source_runs = expect_u64(
        baseline.take("source_runs")?,
        &format!("{what}.source_runs"),
    )?;
    if usize::try_from(source_runs).ok() != Some(source_receipts.len()) {
        return Err(invalid_receipt(
            format!("{what}.source_runs"),
            "does not equal source_receipts length",
        ));
    }
    let age_policy_days = expect_u64(
        baseline.take("age_policy_days")?,
        &format!("{what}.age_policy_days"),
    )?;
    if !(1..=MAX_BASELINE_AGE_DAYS).contains(&age_policy_days) {
        return Err(invalid_receipt(
            format!("{what}.age_policy_days"),
            "outside 1..=365",
        ));
    }
    baseline.finish()?;
    let mut identity = DecodedIdentity {
        fingerprint,
        cpu_brand: cpu_brand.clone(),
        logical_cpus,
        os: os.clone(),
        arch: arch.clone(),
        firmware: firmware.clone(),
        canonical: String::new(),
    };
    identity.canonical = identity_json(&identity);
    let canonical = format!(
        "{{\"schema_version\":1,\"low_band_bits\":\"{:016x}\",\"high_band_bits\":\"{:016x}\",\"promotion_drift_bits\":\"{:016x}\",\"fingerprint\":\"{fingerprint:016x}\",\"cpu_brand\":{},\"logical_cpus\":{logical_cpus},\"os\":{},\"arch\":{},\"firmware\":{},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\",\"source_receipts\":{},\"promoted_by\":{},\"justification\":{},\"promoted_day\":{promoted_day},\"source_runs\":{source_runs},\"age_policy_days\":{age_policy_days}}}",
        0.70_f64.to_bits(),
        1.15_f64.to_bits(),
        0.25_f64.to_bits(),
        canonical_json_string(&cpu_brand),
        canonical_json_string(&os),
        canonical_json_string(&arch),
        canonical_json_string(&firmware),
        bits[0],
        bits[1],
        bits[2],
        bits[3],
        hash_array_json(&source_receipts),
        canonical_json_string(&promoted_by),
        canonical_json_string(&justification),
    );
    if canonical.len() > MAX_BASELINE_LINE_BYTES {
        return Err(invalid_receipt(
            what,
            format!("canonical baseline exceeds {MAX_BASELINE_LINE_BYTES} bytes"),
        ));
    }
    Ok(DecodedBaseline {
        identity,
        bits,
        source_receipts,
        promoted_day,
        age_policy_days,
        canonical,
    })
}

fn validate_trusted_axis_math(
    admission_day: u64,
    baseline: &DecodedBaseline,
    pre: &DecodedAxes,
    post: &DecodedAxes,
) -> Result<(), TuneModelError> {
    if admission_day < baseline.promoted_day
        || admission_day - baseline.promoted_day > baseline.age_policy_days
    {
        return Err(invalid_receipt(
            "op.ir.baseline_admission.verdict",
            "trusted verdict contradicts baseline age policy",
        ));
    }
    let pre_values = axes_values(pre.bits);
    let post_values = axes_values(post.bits);
    let baseline_values = axes_values(baseline.bits);
    for index in 0..4 {
        let drift = (pre_values[index] - post_values[index]).abs()
            / pre_values[index].abs().max(post_values[index].abs());
        let pre_ratio = pre_values[index] / baseline_values[index];
        let post_ratio = post_values[index] / baseline_values[index];
        if drift > 0.25
            || !(0.70..=1.15).contains(&pre_ratio)
            || !(0.70..=1.15).contains(&post_ratio)
        {
            return Err(invalid_receipt(
                "op.ir.baseline_admission.verdict",
                "trusted verdict contradicts axis drift or baseline bands",
            ));
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn validate_axis_admission(
    value: JsonValue,
    binding: &RowBinding,
    receipt: &ReceiptObservation,
    fingerprint: u64,
    post_fingerprint: u64,
    pre_axes_receipt: &str,
    post_axes_receipt: &str,
) -> Result<(ValidatedAxisAdmission, String), TuneModelError> {
    let what = "op.ir.baseline_admission";
    let mut admission = ObjectFields::new(value, what)?;
    if expect_string(admission.take("schema")?, &format!("{what}.schema"))?
        != "fs-roofline-axis-admission-v2"
    {
        return Err(invalid_receipt(
            format!("{what}.schema"),
            "expected fs-roofline-axis-admission-v2",
        ));
    }
    if expect_string(admission.take("tier")?, &format!("{what}.tier"))? != "attested" {
        return Err(invalid_receipt(
            format!("{what}.tier"),
            "production evidence requires attested tier",
        ));
    }
    let now_day = expect_u64(admission.take("now_day")?, &format!("{what}.now_day"))?;
    let decision_day = expect_u64(
        admission.take("decision_day")?,
        &format!("{what}.decision_day"),
    )?;
    if now_day != decision_day {
        return Err(invalid_receipt(
            format!("{what}.decision_day"),
            "must equal now_day",
        ));
    }
    let identity = decode_identity(admission.take("identity")?, &format!("{what}.identity"))?;
    let pre = decode_axes(admission.take("pre")?, &format!("{what}.pre"))?;
    let post = decode_axes(admission.take("post")?, &format!("{what}.post"))?;
    let baseline_hash = expect_hash(
        admission.take("baseline_hash")?,
        &format!("{what}.baseline_hash"),
    )?;
    let baseline = decode_baseline(admission.take("baseline")?)?;
    let mut attestation = ObjectFields::new(
        admission.take("attestation")?,
        format!("{what}.attestation"),
    )?;
    let key_id = expect_string(
        attestation.take("key_id")?,
        &format!("{what}.attestation.key_id"),
    )?;
    let signature = expect_string(
        attestation.take("signature")?,
        &format!("{what}.attestation.signature"),
    )?;
    attestation.finish()?;
    validate_authority_text(&key_id, &format!("{what}.attestation.key_id"))?;
    validate_authority_text(&signature, &format!("{what}.attestation.signature"))?;
    let required_sources = decode_hash_array(
        admission.take("required_source_receipts")?,
        &format!("{what}.required_source_receipts"),
    )?;
    let mut authority =
        ObjectFields::new(admission.take("authority")?, format!("{what}.authority"))?;
    if expect_string(
        authority.take("verdict")?,
        &format!("{what}.authority.verdict"),
    )? != "authorized"
    {
        return Err(invalid_receipt(
            format!("{what}.authority.verdict"),
            "production evidence requires authorized authority verdict",
        ));
    }
    let policy_receipt = expect_hash(
        authority.take("policy_receipt")?,
        &format!("{what}.authority.policy_receipt"),
    )?;
    authority.finish()?;
    let mut verdict = ObjectFields::new(admission.take("verdict")?, format!("{what}.verdict"))?;
    if expect_string(
        verdict.take("baseline")?,
        &format!("{what}.verdict.baseline"),
    )? != "trusted"
    {
        return Err(invalid_receipt(
            format!("{what}.verdict.baseline"),
            "production evidence requires trusted baseline verdict",
        ));
    }
    verdict.finish()?;
    admission.finish()?;

    if identity != baseline.identity
        || identity.fingerprint != pre.fingerprint
        || identity.fingerprint != post.fingerprint
        || identity.fingerprint != fingerprint
        || identity.cpu_brand != pre.cpu_brand
        || identity.cpu_brand != post.cpu_brand
        || identity.logical_cpus != pre.logical_cpus
        || identity.logical_cpus != post.logical_cpus
        || post.fingerprint != post_fingerprint
    {
        return Err(invalid_receipt(
            format!("{what}.identity"),
            "baseline, probe, or operation identity mismatch",
        ));
    }
    if pre.logical_cpus != receipt.logical_cpus
        || pre.bits != receipt.axis_bits
        || post.bits != binding.post_axis_bits
    {
        return Err(invalid_receipt(
            format!("{what}.axes"),
            "probe axes do not bind the measured receipt and row parameters",
        ));
    }
    if required_sources != baseline.source_receipts {
        return Err(invalid_receipt(
            format!("{what}.required_source_receipts"),
            "must equal the canonical baseline provenance",
        ));
    }
    let rederived_baseline_hash =
        fs_blake3::hash_domain(BASELINE_HASH_DOMAIN, baseline.canonical.as_bytes()).to_string();
    if baseline_hash != rederived_baseline_hash || baseline_hash != binding.baseline_hash {
        return Err(TuneModelError::ScopeMismatch {
            field: "baseline_admission.baseline_hash",
        });
    }
    if fs_blake3::hash_domain(PRODUCTION_AXES_RECEIPT_DOMAIN, pre.canonical.as_bytes()).to_string()
        != pre_axes_receipt
        || fs_blake3::hash_domain(PRODUCTION_AXES_RECEIPT_DOMAIN, post.canonical.as_bytes())
            .to_string()
            != post_axes_receipt
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "axis receipt",
        });
    }
    validate_trusted_axis_math(now_day, &baseline, &pre, &post)?;
    let canonical = format!(
        "{{\"schema\":\"fs-roofline-axis-admission-v2\",\"tier\":\"attested\",\"now_day\":{now_day},\"decision_day\":{decision_day},\"identity\":{},\"pre\":{},\"post\":{},\"baseline_hash\":\"{baseline_hash}\",\"baseline\":{},\"attestation\":{{\"key_id\":{},\"signature\":{}}},\"required_source_receipts\":{},\"authority\":{{\"verdict\":\"authorized\",\"policy_receipt\":\"{policy_receipt}\"}},\"verdict\":{{\"baseline\":\"trusted\"}}}}",
        identity.canonical,
        pre.canonical,
        post.canonical,
        baseline.canonical,
        canonical_json_string(&key_id),
        canonical_json_string(&signature),
        hash_array_json(&required_sources),
    );
    Ok((
        ValidatedAxisAdmission {
            decision_day,
            baseline_hash,
            pre,
            post,
        },
        canonical,
    ))
}

fn validate_result_manifest(
    value: JsonValue,
    expected_kernel: &str,
    expected_version: &str,
    expected_payload: &str,
    expected_count: u64,
) -> Result<(String, Vec<ManifestEntry>), TuneModelError> {
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
    let sealed_count = u64::try_from(SEALED_PRODUCTION_REGISTRY.len())
        .expect("sealed production registry length fits u64");
    if expected_count != sealed_count || u64::try_from(entries.len()).ok() != Some(sealed_count) {
        return Err(invalid_receipt(
            "op.ir.result_manifest.entries",
            format!(
                "production-v3 requires exactly {sealed_count} sealed registry entries matching op.ir.kernels"
            ),
        ));
    }
    let mut matching_entries = 0_usize;
    let mut decoded_entries = Vec::with_capacity(entries.len());
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
        let (sealed_kernel, sealed_version) = SEALED_PRODUCTION_REGISTRY[expected_ordinal];
        if kernel != sealed_kernel || version != sealed_version {
            return Err(invalid_receipt(
                format!("op.ir.result_manifest.entries[{expected_ordinal}].identity"),
                format!(
                    "expected sealed production-v3 member {sealed_kernel}/{sealed_version}, got {kernel}/{version}"
                ),
            ));
        }
        if kernel == expected_kernel && version == expected_version && payload == expected_payload {
            matching_entries += 1;
        }
        decoded_entries.push(ManifestEntry {
            ordinal,
            kernel,
            version,
            payload,
        });
    }
    if matching_entries != 1 {
        return Err(invalid_receipt(
            "op.ir.result_manifest.entries",
            "exactly one manifest member must bind this kernel/version/payload",
        ));
    }
    let canonical = format!(
        "{{\"schema\":\"fs-roofline-run-manifest-v1\",\"entries\":[{}]}}",
        decoded_entries
            .iter()
            .map(|entry| format!(
                "{{\"ordinal\":{},\"kernel\":{},\"version\":{},\"payload\":\"{}\"}}",
                entry.ordinal,
                canonical_json_string(&entry.kernel),
                canonical_json_string(&entry.version),
                entry.payload
            ))
            .collect::<Vec<_>>()
            .join(",")
    );
    Ok((canonical, decoded_entries))
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
) -> Result<ValidatedProductionOp, TuneModelError> {
    let original_ir = ir;
    let mut ir = ObjectFields::new(
        StrictJson::parse(ir, fs_ledger::MAX_TUNE_PARAMS_BYTES)?,
        "op.ir",
    )?;
    if expect_string(ir.take("op")?, "op.ir.op")? != "perf.roofline" {
        return Err(invalid_receipt("op.ir.op", "expected perf.roofline"));
    }
    let kernel_count = expect_u64(ir.take("kernels")?, "op.ir.kernels")?;
    let sealed_kernel_count = u64::try_from(SEALED_PRODUCTION_REGISTRY.len())
        .expect("sealed production registry length fits u64");
    if kernel_count != sealed_kernel_count {
        return Err(invalid_receipt(
            "op.ir.kernels",
            format!("production-v3 requires exactly {sealed_kernel_count} kernels"),
        ));
    }
    let fingerprint = expect_hex_u64(ir.take("fingerprint")?, "op.ir.fingerprint")?;
    let post_fingerprint = expect_hex_u64(ir.take("post_fingerprint")?, "op.ir.post_fingerprint")?;
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
    if expect_string(ir.take("protocol")?, "op.ir.protocol")? != "production-v3" {
        return Err(invalid_receipt("op.ir.protocol", "expected production-v3"));
    }
    let run_nonce = expect_hash(ir.take("run_nonce")?, "op.ir.run_nonce")?;
    let pre_axes_receipt = expect_hash(ir.take("pre_axes_receipt")?, "op.ir.pre_axes_receipt")?;
    let post_axes_receipt = expect_hash(ir.take("post_axes_receipt")?, "op.ir.post_axes_receipt")?;
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
    let dependency_receipt_digest = expect_hash(
        ir.take("dependency_receipt_digest")?,
        "op.ir.dependency_receipt_digest",
    )?;
    if dependency_receipt_digest != binding.dependency_receipt_digest {
        return Err(TuneModelError::ScopeMismatch {
            field: "dependency_receipt_digest",
        });
    }
    let dependency_receipt_artifact = expect_hash(
        ir.take("dependency_receipt_artifact")?,
        "op.ir.dependency_receipt_artifact",
    )?;
    if dependency_receipt_artifact != binding.dependency_receipt_artifact {
        return Err(TuneModelError::ScopeMismatch {
            field: "dependency_receipt_artifact",
        });
    }
    let finalized_run_receipt = expect_hash(
        ir.take("finalized_run_receipt")?,
        "op.ir.finalized_run_receipt",
    )?;
    if finalized_run_receipt != binding.run_receipt {
        return Err(TuneModelError::ScopeMismatch {
            field: "finalized_run_receipt",
        });
    }
    let (result_manifest, manifest) = validate_result_manifest(
        ir.take("result_manifest")?,
        &receipt.kernel,
        &receipt.version,
        &binding.payload_artifact,
        kernel_count,
    )?;
    let (admission, baseline_admission) = validate_axis_admission(
        ir.take("baseline_admission")?,
        binding,
        receipt,
        fingerprint,
        post_fingerprint,
        &pre_axes_receipt,
        &post_axes_receipt,
    )?;
    ir.finish()?;
    let canonical = format!(
        "{{\"op\":\"perf.roofline\",\"kernels\":{kernel_count},\"fingerprint\":\"{fingerprint:016x}\",\"post_fingerprint\":\"{post_fingerprint:016x}\",\"measurement_admitted\":true,\"admitted\":true,\"protocol\":\"production-v3\",\"run_nonce\":\"{run_nonce}\",\"pre_axes_receipt\":\"{pre_axes_receipt}\",\"post_axes_receipt\":\"{post_axes_receipt}\",\"dependency_graph_evidence\":\"operator-observed-receipt\",\"dependency_receipt_digest\":\"{dependency_receipt_digest}\",\"dependency_receipt_artifact\":\"{dependency_receipt_artifact}\",\"finalized_run_receipt\":\"{finalized_run_receipt}\",\"result_manifest\":{result_manifest},\"baseline_admission\":{baseline_admission}}}"
    );
    if canonical != original_ir {
        return Err(invalid_receipt(
            "op.ir",
            "production envelope is not byte-identical to canonical serialization",
        ));
    }
    Ok(ValidatedProductionOp {
        baseline_admission,
        result_manifest,
        manifest,
        finalized_run_receipt,
        admission,
    })
}

fn content_hash(text: &str, field: &str) -> Result<fs_ledger::ContentHash, TuneModelError> {
    fs_ledger::ContentHash::from_hex(text)
        .ok_or_else(|| invalid_receipt(field, "invalid content hash"))
}

fn wall_ns_day(wall_ns: i64) -> Result<u64, TuneModelError> {
    let wall_ns = u64::try_from(wall_ns)
        .map_err(|_| invalid_receipt("op.t_end", "must be a nonnegative Unix timestamp"))?;
    Ok(wall_ns / WALL_NS_PER_DAY)
}

fn push_receipt_field(payload: &mut Vec<u8>, bytes: &[u8]) -> Result<(), TuneModelError> {
    let len = u64::try_from(bytes.len())
        .map_err(|_| invalid_receipt("finalized_run_receipt", "field length exceeds u64"))?;
    payload.extend_from_slice(&len.to_le_bytes());
    payload.extend_from_slice(bytes);
    Ok(())
}

fn validate_payload_artifact(
    ledger: &Ledger,
    row: &fs_ledger::TuneRow,
    binding: &RowBinding,
) -> Result<(), TuneModelError> {
    let payload_hash = content_hash(&binding.payload_artifact, "params.payload_artifact")?;
    let payload_info =
        ledger
            .artifact_info(&payload_hash)?
            .ok_or(TuneModelError::ScopeMismatch {
                field: "payload_artifact",
            })?;
    let payload_limit = u64::try_from(row.measured.len())
        .map_err(|_| invalid_receipt("measured", "payload length exceeds u64"))?;
    let op_id = i64::try_from(binding.op)
        .map_err(|_| invalid_receipt("params.op", "operation id exceeds i64"))?;
    if payload_info.kind != "roofline-benchmark-result"
        || payload_info.meta.as_deref() != Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}")
        || u64::try_from(row.measured.len()).ok() != Some(payload_info.len)
        || ledger
            .get_artifact_bounded(&payload_hash, payload_limit)?
            .as_deref()
            != Some(row.measured.as_bytes())
        || !ledger.edge_exists(op_id, &payload_hash, fs_ledger::EdgeRole::Out)?
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "payload_artifact",
        });
    }
    Ok(())
}

fn manifest_shape(version: &str, run_receipt: &str, op: u64) -> String {
    format!("{ROOFLINE_TUNE_SHAPE_PREFIX}:{version}:run={run_receipt}:op={op}")
}

fn checked_production_product(
    left: u128,
    right: u128,
    component: &str,
) -> Result<u128, TuneModelError> {
    left.checked_mul(right).ok_or_else(|| {
        invalid_receipt(
            "result_manifest production profile",
            format!("production-v3 {component} multiplication overflowed u128"),
        )
    })
}

fn checked_production_sum(
    left: u128,
    right: u128,
    component: &str,
) -> Result<u128, TuneModelError> {
    left.checked_add(right).ok_or_else(|| {
        invalid_receipt(
            "result_manifest production profile",
            format!("production-v3 {component} addition overflowed u128"),
        )
    })
}

fn production_registry_work(n: u64, runs_per_kernel: u64) -> Result<(u128, u128), TuneModelError> {
    let side = u128::from(n.isqrt().max(256));
    let n = u128::from(n);
    let gemm_outputs = checked_production_product(side, side, "GEMM output extent")?;
    let gemm_flops_per_output = checked_production_product(side, 2, "GEMM per-output FLOPs")?;
    let vector_flops = checked_production_product(n, 5, "vector FLOPs")?;
    let gemm_flops = checked_production_product(gemm_outputs, gemm_flops_per_output, "GEMM FLOPs")?;
    let flops_per_run = checked_production_sum(vector_flops, gemm_flops, "registry FLOPs")?;
    let vector_bytes = checked_production_product(n, 48, "vector bytes")?;
    let gemm_bytes = checked_production_product(gemm_outputs, 24, "GEMM bytes")?;
    let bytes_per_run = checked_production_sum(vector_bytes, gemm_bytes, "registry bytes")?;
    let runs = u128::from(runs_per_kernel);
    Ok((
        checked_production_product(flops_per_run, runs, "total FLOPs")?,
        checked_production_product(bytes_per_run, runs, "total bytes")?,
    ))
}

fn validate_production_registry_profile(
    receipts: &[ReceiptObservation],
) -> Result<(), TuneModelError> {
    if receipts.len() != SEALED_PRODUCTION_REGISTRY.len() {
        return Err(invalid_receipt(
            "result_manifest production profile",
            format!(
                "production-v3 requires exactly {} result receipts",
                SEALED_PRODUCTION_REGISTRY.len()
            ),
        ));
    }
    let target = &receipts[0];
    let runs_per_kernel = validate_production_run_counts(target.warmup_runs, target.reps)?;
    for (index, (receipt, (expected_kernel, expected_version))) in
        receipts.iter().zip(SEALED_PRODUCTION_REGISTRY).enumerate()
    {
        if receipt.kernel != expected_kernel || receipt.version != expected_version {
            return Err(invalid_receipt(
                format!("result_manifest production profile[{index}].identity"),
                format!(
                    "expected {expected_kernel}/{expected_version}, got {}/{}",
                    receipt.kernel, receipt.version
                ),
            ));
        }
        if receipt.elements == 0 || receipt.elements > MAX_PRODUCTION_ELEMENTS {
            return Err(invalid_receipt(
                format!("result_manifest production profile[{index}].elements"),
                format!("must be in 1..={MAX_PRODUCTION_ELEMENTS}"),
            ));
        }
        if receipt.logical_cpus == 0 || receipt.logical_cpus > MAX_PRODUCTION_LOGICAL_CPUS {
            return Err(invalid_receipt(
                format!("result_manifest production profile[{index}].logical_cpus"),
                format!("must be in 1..={MAX_PRODUCTION_LOGICAL_CPUS}"),
            ));
        }
        if receipt.logical_cpus != target.logical_cpus
            || receipt.warmup_runs != target.warmup_runs
            || receipt.reps != target.reps
        {
            return Err(TuneModelError::ScopeMismatch {
                field: "result_manifest sibling configuration",
            });
        }
    }

    let n = receipts[0].elements;
    if receipts[1].elements != n || receipts[2].elements != n {
        return Err(TuneModelError::ScopeMismatch {
            field: "result_manifest vector elements",
        });
    }
    let gemm_side = n.isqrt().max(256);
    let expected_gemm_elements = gemm_side.checked_mul(gemm_side).ok_or_else(|| {
        invalid_receipt(
            "result_manifest GEMM elements",
            "derived GEMM output extent overflowed u64",
        )
    })?;
    if receipts[3].elements != expected_gemm_elements {
        return Err(TuneModelError::ScopeMismatch {
            field: "result_manifest GEMM elements",
        });
    }

    let (total_flops, total_bytes) = production_registry_work(n, runs_per_kernel)?;
    if total_flops > MAX_PRODUCTION_REGISTRY_FLOPS {
        return Err(invalid_receipt(
            "result_manifest production profile",
            format!(
                "production-v3 requires {total_flops} modeled FLOPs, exceeding {MAX_PRODUCTION_REGISTRY_FLOPS}"
            ),
        ));
    }
    if total_bytes > MAX_PRODUCTION_REGISTRY_BYTES {
        return Err(invalid_receipt(
            "result_manifest production profile",
            format!(
                "production-v3 requires {total_bytes} modeled logical bytes, exceeding {MAX_PRODUCTION_REGISTRY_BYTES}"
            ),
        ));
    }
    Ok(())
}

fn validate_manifest_sibling_configuration(
    binding: &RowBinding,
    receipt: &ReceiptObservation,
    shared: &RowBinding,
    target: &ReceiptObservation,
) -> Result<(), TuneModelError> {
    if binding.reps != shared.reps
        || receipt.reps != binding.reps
        || receipt.reps != target.reps
        || receipt.warmup_runs != target.warmup_runs
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "result_manifest sibling configuration",
        });
    }
    Ok(())
}

fn validate_manifest_sibling(
    ledger: &Ledger,
    entry: &ManifestEntry,
    protocol: &ValidatedProductionOp,
    shared: &RowBinding,
    target: &ReceiptObservation,
    machine: &[u8],
) -> Result<(String, ReceiptObservation), TuneModelError> {
    let shape = manifest_shape(&entry.version, &protocol.finalized_run_receipt, shared.op);
    let row =
        ledger
            .tune_get(&entry.kernel, &shape, machine)?
            .ok_or(TuneModelError::ScopeMismatch {
                field: "result_manifest sibling",
            })?;
    let binding = decode_row_binding(&row.params)?;
    let receipt = decode_receipt(&row.measured)?;
    validate_manifest_sibling_configuration(&binding, &receipt, shared, target)?;
    if row.kernel != entry.kernel
        || row.shape_class != shape
        || row.machine != machine
        || binding.op != shared.op
        || binding.run_receipt != protocol.finalized_run_receipt
        || binding.payload_artifact != entry.payload
        || binding.dependency_receipt_artifact != shared.dependency_receipt_artifact
        || binding.dependency_receipt_digest != shared.dependency_receipt_digest
        || binding.baseline_hash != shared.baseline_hash
        || binding.build_identity != shared.build_identity
        || binding.post_axis_bits != shared.post_axis_bits
        || receipt.kernel != entry.kernel
        || receipt.version != entry.version
        || receipt.machine != protocol.admission.pre.fingerprint
        || receipt.logical_cpus != protocol.admission.pre.logical_cpus
        || receipt.axis_bits != protocol.admission.pre.bits
        || fs_ledger::hash_bytes(row.measured.as_bytes()).to_string() != entry.payload
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "result_manifest sibling",
        });
    }
    validate_payload_artifact(ledger, &row, &binding)?;
    Ok((row.measured, receipt))
}

fn validate_finalized_run_receipt(
    ledger: &Ledger,
    protocol: &ValidatedProductionOp,
    binding: &RowBinding,
    target: &ReceiptObservation,
    machine: &[u8],
) -> Result<(), TuneModelError> {
    if protocol.manifest.len() != SEALED_PRODUCTION_REGISTRY.len() {
        return Err(invalid_receipt(
            "op.ir.result_manifest.entries",
            format!(
                "finalized production-v3 run requires exactly {} manifest entries",
                SEALED_PRODUCTION_REGISTRY.len()
            ),
        ));
    }
    let mut result_payloads = Vec::with_capacity(protocol.manifest.len());
    let mut receipts = Vec::with_capacity(protocol.manifest.len());
    for entry in &protocol.manifest {
        let (payload, receipt) =
            validate_manifest_sibling(ledger, entry, protocol, binding, target, machine)?;
        result_payloads.push(payload);
        receipts.push(receipt);
    }
    validate_production_registry_profile(&receipts)?;
    if finalized_run_receipt_for_payloads(
        &protocol.baseline_admission,
        &result_payloads,
        &protocol.result_manifest,
    )? != protocol.finalized_run_receipt
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "finalized_run_receipt",
        });
    }
    Ok(())
}

fn finalized_run_receipt_for_payloads(
    baseline_admission: &str,
    result_payloads: &[String],
    result_manifest: &str,
) -> Result<String, TuneModelError> {
    let mut payload = Vec::new();
    push_receipt_field(&mut payload, baseline_admission.as_bytes())?;
    let result_count = u64::try_from(result_payloads.len())
        .map_err(|_| invalid_receipt("op.ir.result_manifest", "entry count exceeds u64"))?;
    payload.extend_from_slice(&result_count.to_le_bytes());
    for measured in result_payloads {
        push_receipt_field(&mut payload, measured.as_bytes())?;
    }
    let manifest_hash = fs_blake3::hash_domain(RESULT_MANIFEST_DOMAIN, result_manifest.as_bytes());
    push_receipt_field(&mut payload, manifest_hash.as_bytes())?;
    Ok(fs_blake3::hash_domain(FINALIZED_RUN_DOMAIN, &payload).to_string())
}

/// Full provenance validation; returns the evidence operation's
/// completion time so the sealed scope can retain it.
fn validate_provenance(
    ledger: &Ledger,
    row: &fs_ledger::TuneRow,
    binding: &RowBinding,
    receipt: &ReceiptObservation,
    machine: &[u8],
) -> Result<i64, TuneModelError> {
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
    let recorded_at_ns = op.t_end.ok_or(TuneModelError::ScopeMismatch {
        field: "operation envelope",
    })?;
    if op.id != op_id
        || op.session.as_deref() != Some(b"roofline".as_slice())
        || op.seed != b"roofline"
        || op.budget != "{\"wall_s\":600}"
        || op.capability != "{\"ops\":[\"perf.roofline\"]}"
        || op.outcome.as_deref() != Some("ok")
        || op.diag.is_some()
        || recorded_at_ns < op.t_start
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "operation envelope",
        });
    }
    validate_versions(&op.versions, &binding.build_identity)?;
    let protocol = validate_op_ir(&op.ir, binding, receipt)?;
    if protocol.admission.baseline_hash != binding.baseline_hash
        || protocol.admission.post.bits != binding.post_axis_bits
        || wall_ns_day(recorded_at_ns)? != protocol.admission.decision_day
    {
        return Err(TuneModelError::ScopeMismatch {
            field: "baseline admission",
        });
    }
    validate_payload_artifact(ledger, row, binding)?;
    validate_finalized_run_receipt(ledger, &protocol, binding, receipt, machine)?;

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
    let dependency_bytes = ledger
        .get_artifact_bounded(&dependency_hash, MAX_DEPGRAPH_RECEIPT_BYTES)?
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
    Ok(recorded_at_ns)
}

/// Rebuild one model from one exact production roofline tune key.
///
/// This API deliberately uses `tune_get`, not a per-kernel scan: foreign
/// machines and neighboring shape classes can never contribute. A current
/// production row contributes each bounded timed repetition as a same-size
/// observation. Rows with fewer than [`crate::cost::MIN_OBS`] repetitions
/// therefore continue to refuse prediction honestly.
///
/// The returned model is SEALED (bead 2pmb): the validation this
/// function performs — receipt scope, operation envelope, build
/// identity, payload and dependency digests, finalized-run receipt —
/// travels with the model as a [`CostModelScope`] under
/// [`crate::sealed::CostEvidenceClass::ExactRooflineReceipt`], which
/// only this loader can mint. A caller-fitted [`CostModel`] can enter
/// consumers only as explicitly provisional evidence.
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
) -> Result<SealedCostModel, TuneModelError> {
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
    let recorded_at_ns = validate_provenance(ledger, &row, &binding, &receipt, machine)?;
    let model = CostModel::fit(&receipt.observations)?;
    Ok(SealedCostModel::mint_exact(
        model,
        CostModelScope::from_validated(
            kernel.to_string(),
            shape_class.to_string(),
            machine.to_vec(),
            binding.run_receipt.clone(),
            binding.op,
            binding.build_identity.clone(),
            recorded_at_ns,
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dependency_receipt_cap_matches_the_fs_la_producer_source() {
        let source = include_str!("../../fs-la/depgraph_receipt_format.rs");
        let declaration = source
            .lines()
            .find_map(|line| {
                line.trim()
                    .strip_prefix("pub const MAX_RECEIPT_BYTES: usize = ")
            })
            .and_then(|value| value.strip_suffix(';'))
            .expect("fs-la producer must declare MAX_RECEIPT_BYTES");
        let producer_cap = declaration
            .replace('_', "")
            .parse::<u64>()
            .expect("fs-la producer cap must remain a decimal byte count");
        assert_eq!(producer_cap, MAX_DEPGRAPH_RECEIPT_BYTES);
    }

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

    struct AdmissionFixture {
        admission: String,
        pre_receipt: String,
        post_receipt: String,
        manifest: String,
        payloads: Vec<String>,
        binding: RowBinding,
        receipt: ReceiptObservation,
    }

    fn numbered_hash(value: u64) -> String {
        format!("{value:064x}")
    }

    fn admission_fixture() -> AdmissionFixture {
        admission_fixture_for_target(std::env::consts::OS, std::env::consts::ARCH)
    }

    fn admission_fixture_for_target(os: &str, arch: &str) -> AdmissionFixture {
        let fingerprint = 0x0102_0304_0506_0708_u64;
        let logical_cpus = 8_u64;
        let axis_bits = [
            100.0_f64.to_bits(),
            200.0_f64.to_bits(),
            1_000.0_f64.to_bits(),
            2_000.0_f64.to_bits(),
        ];
        let cpu_brand = "fixture-cpu";
        let firmware = "fixture-firmware";
        let identity = format!(
            "{{\"fingerprint\":\"{fingerprint:016x}\",\"cpu_brand\":\"{cpu_brand}\",\"logical_cpus\":{logical_cpus},\"os\":\"{os}\",\"arch\":\"{arch}\",\"firmware\":\"{firmware}\"}}"
        );
        let axes = format!(
            "{{\"fingerprint\":\"{fingerprint:016x}\",\"cpu_brand\":\"{cpu_brand}\",\"logical_cpus\":{logical_cpus},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}}",
            axis_bits[0], axis_bits[1], axis_bits[2], axis_bits[3]
        );
        let sources = [numbered_hash(1), numbered_hash(2), numbered_hash(3)];
        let sources_json = format!("[\"{}\",\"{}\",\"{}\"]", sources[0], sources[1], sources[2]);
        let baseline = format!(
            "{{\"schema_version\":1,\"low_band_bits\":\"{:016x}\",\"high_band_bits\":\"{:016x}\",\"promotion_drift_bits\":\"{:016x}\",\"fingerprint\":\"{fingerprint:016x}\",\"cpu_brand\":\"{cpu_brand}\",\"logical_cpus\":{logical_cpus},\"os\":\"{os}\",\"arch\":\"{arch}\",\"firmware\":\"{firmware}\",\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\",\"source_receipts\":{sources_json},\"promoted_by\":\"fixture-operator\",\"justification\":\"fixture baseline\",\"promoted_day\":100,\"source_runs\":3,\"age_policy_days\":90}}",
            0.70_f64.to_bits(),
            1.15_f64.to_bits(),
            0.25_f64.to_bits(),
            axis_bits[0],
            axis_bits[1],
            axis_bits[2],
            axis_bits[3],
        );
        let baseline_hash =
            fs_blake3::hash_domain(BASELINE_HASH_DOMAIN, baseline.as_bytes()).to_string();
        let policy_receipt = numbered_hash(9);
        let admission = format!(
            "{{\"schema\":\"fs-roofline-axis-admission-v2\",\"tier\":\"attested\",\"now_day\":100,\"decision_day\":100,\"identity\":{identity},\"pre\":{axes},\"post\":{axes},\"baseline_hash\":\"{baseline_hash}\",\"baseline\":{baseline},\"attestation\":{{\"key_id\":\"fixture-key\",\"signature\":\"fixture-signature\"}},\"required_source_receipts\":{sources_json},\"authority\":{{\"verdict\":\"authorized\",\"policy_receipt\":\"{policy_receipt}\"}},\"verdict\":{{\"baseline\":\"trusted\"}}}}"
        );
        let pre_receipt =
            fs_blake3::hash_domain(PRODUCTION_AXES_RECEIPT_DOMAIN, axes.as_bytes()).to_string();
        let post_receipt = pre_receipt.clone();
        let payloads = SEALED_PRODUCTION_REGISTRY
            .iter()
            .enumerate()
            .map(|(ordinal, _)| format!("fixture-result-payload-{ordinal}"))
            .collect::<Vec<_>>();
        let payload_hashes = payloads
            .iter()
            .map(|payload| fs_ledger::hash_bytes(payload.as_bytes()).to_string())
            .collect::<Vec<_>>();
        let manifest = format!(
            "{{\"schema\":\"fs-roofline-run-manifest-v1\",\"entries\":[{}]}}",
            SEALED_PRODUCTION_REGISTRY
                .iter()
                .zip(&payload_hashes)
                .enumerate()
                .map(|(ordinal, ((kernel, version), payload))| format!(
                    "{{\"ordinal\":{ordinal},\"kernel\":\"{kernel}\",\"version\":\"{version}\",\"payload\":\"{payload}\"}}"
                ))
                .collect::<Vec<_>>()
                .join(",")
        );
        let run_receipt =
            finalized_run_receipt_for_payloads(&admission, &payloads, &manifest).unwrap();
        AdmissionFixture {
            admission,
            pre_receipt,
            post_receipt,
            manifest,
            payloads,
            binding: RowBinding {
                op: 7,
                run_receipt,
                payload_artifact: payload_hashes[0].clone(),
                dependency_receipt_artifact: numbered_hash(10),
                dependency_receipt_digest: numbered_hash(11),
                baseline_hash,
                build_identity: numbered_hash(12),
                reps: 3,
                post_axis_bits: axis_bits,
            },
            receipt: ReceiptObservation {
                kernel: "simd-axpy-f64".to_string(),
                version: "1".to_string(),
                machine: fingerprint,
                logical_cpus,
                axis_bits,
                elements: 1_000,
                warmup_runs: 1,
                observations: Vec::new(),
                reps: 3,
            },
        }
    }

    fn fixture_op_ir(fixture: &AdmissionFixture, admission: &str) -> String {
        format!(
            "{{\"op\":\"perf.roofline\",\"kernels\":{},\"fingerprint\":\"{:016x}\",\"post_fingerprint\":\"{:016x}\",\"measurement_admitted\":true,\"admitted\":true,\"protocol\":\"production-v3\",\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"dependency_graph_evidence\":\"operator-observed-receipt\",\"dependency_receipt_digest\":\"{}\",\"dependency_receipt_artifact\":\"{}\",\"finalized_run_receipt\":\"{}\",\"result_manifest\":{},\"baseline_admission\":{admission}}}",
            SEALED_PRODUCTION_REGISTRY.len(),
            fixture.receipt.machine,
            fixture.receipt.machine,
            numbered_hash(13),
            fixture.pre_receipt,
            fixture.post_receipt,
            fixture.binding.dependency_receipt_digest,
            fixture.binding.dependency_receipt_artifact,
            fixture.binding.run_receipt,
            fixture.manifest,
        )
    }

    #[test]
    fn production_v3_consumer_refuses_empty_candidate_and_untrusted_admission() {
        let fixture = admission_fixture();
        let valid = fixture_op_ir(&fixture, &fixture.admission);
        assert!(validate_op_ir(&valid, &fixture.binding, &fixture.receipt).is_ok());

        for invalid in [
            "{}".to_string(),
            fixture
                .admission
                .replacen("\"tier\":\"attested\"", "\"tier\":\"candidate\"", 1),
            fixture.admission.replacen(
                "\"verdict\":\"authorized\"",
                "\"verdict\":\"revoked-key\"",
                1,
            ),
            fixture
                .admission
                .replacen("\"baseline\":\"trusted\"", "\"baseline\":\"degraded\"", 1),
        ] {
            assert!(
                validate_op_ir(
                    &fixture_op_ir(&fixture, &invalid),
                    &fixture.binding,
                    &fixture.receipt,
                )
                .is_err(),
                "admission unexpectedly accepted: {invalid}"
            );
        }
    }

    #[test]
    fn production_v3_replay_is_independent_of_the_consumer_host_target() {
        let fixture = admission_fixture_for_target("historical-os", "historical-arch");
        let valid = fixture_op_ir(&fixture, &fixture.admission);
        assert!(
            validate_op_ir(&valid, &fixture.binding, &fixture.receipt).is_ok(),
            "a self-consistent historical receipt must replay on a different audit host"
        );

        let mismatched_identity =
            fixture
                .admission
                .replacen("\"os\":\"historical-os\"", "\"os\":\"substituted-os\"", 1);
        assert!(
            validate_op_ir(
                &fixture_op_ir(&fixture, &mismatched_identity),
                &fixture.binding,
                &fixture.receipt,
            )
            .is_err(),
            "identity disagreement inside the historical receipt must still fail closed"
        );
    }

    #[test]
    fn production_v3_manifest_requires_the_exact_ordered_sealed_registry() {
        let fixture = admission_fixture();
        let wrong_count = fixture_op_ir(&fixture, &fixture.admission).replacen(
            &format!("\"kernels\":{}", SEALED_PRODUCTION_REGISTRY.len()),
            "\"kernels\":3",
            1,
        );
        assert!(validate_op_ir(&wrong_count, &fixture.binding, &fixture.receipt).is_err());

        let wrong_first_member = fixture.manifest.replacen(
            "\"kernel\":\"simd-axpy-f64\"",
            "\"kernel\":\"simd-dot-f64\"",
            1,
        );
        let wrong_manifest = fixture_op_ir(&fixture, &fixture.admission).replacen(
            &fixture.manifest,
            &wrong_first_member,
            1,
        );
        assert!(validate_op_ir(&wrong_manifest, &fixture.binding, &fixture.receipt).is_err());
    }

    fn production_profile(n: u64, warmup_runs: u64, reps: u64) -> Vec<ReceiptObservation> {
        let gemm_side = n.isqrt().max(256);
        SEALED_PRODUCTION_REGISTRY
            .iter()
            .enumerate()
            .map(|(index, (kernel, version))| ReceiptObservation {
                kernel: (*kernel).to_string(),
                version: (*version).to_string(),
                machine: 0x0102_0304_0506_0708,
                logical_cpus: 8,
                axis_bits: [1; 4],
                elements: if index == 3 { gemm_side * gemm_side } else { n },
                warmup_runs,
                observations: Vec::new(),
                reps,
            })
            .collect()
    }

    #[test]
    fn production_v3_profile_enforces_exact_run_and_aggregate_work_caps() {
        assert_eq!(
            production_registry_work(1, 1).unwrap(),
            (33_554_437, 1_572_912),
            "consumer work algebra must include the floor GEMM and all vector kernels"
        );
        assert_eq!(validate_production_run_counts(63, 1).unwrap(), 64);
        assert!(validate_production_run_counts(64, 1).is_err());
        assert!(validate_production_run_counts(0, 65).is_err());
        assert!(validate_production_run_counts(63, 2).is_err());

        validate_production_registry_profile(&production_profile(1, 0, 64))
            .expect("minimum shape admits the exact per-kernel run cap");
        validate_production_registry_profile(&production_profile(1 << 24, 2, 1))
            .expect("three maximum-shape runs remain inside the FLOP cap");
        let flops = validate_production_registry_profile(&production_profile(1 << 24, 3, 1))
            .expect_err("four maximum-shape runs exceed the FLOP cap");
        assert!(format!("{flops}").contains("modeled FLOPs"), "{flops}");

        validate_production_registry_profile(&production_profile(1 << 22, 27, 1))
            .expect("twenty-eight default-shape runs remain inside the byte cap");
        let bytes = validate_production_registry_profile(&production_profile(1 << 22, 28, 1))
            .expect_err("twenty-nine default-shape runs exceed the byte cap");
        assert!(
            format!("{bytes}").contains("modeled logical bytes"),
            "{bytes}"
        );
    }

    #[test]
    fn production_v3_profile_and_sibling_configuration_fail_closed() {
        let valid = production_profile(1_000, 1, 3);
        validate_production_registry_profile(&valid).expect("valid sealed profile");

        let mut reordered = valid.clone();
        reordered.swap(0, 1);
        assert!(validate_production_registry_profile(&reordered).is_err());
        let mut changed_vector = valid.clone();
        changed_vector[1].elements += 1;
        assert!(validate_production_registry_profile(&changed_vector).is_err());
        let mut changed_gemm = valid.clone();
        changed_gemm[3].elements += 1;
        assert!(validate_production_registry_profile(&changed_gemm).is_err());
        let mut changed_warmup = valid.clone();
        changed_warmup[2].warmup_runs += 1;
        assert!(validate_production_registry_profile(&changed_warmup).is_err());

        let fixture = admission_fixture();
        let target = fixture.receipt.clone();
        let sibling = target.clone();
        let shared = fixture.binding.clone();
        let mut sibling_binding = shared.clone();
        validate_manifest_sibling_configuration(&sibling_binding, &sibling, &shared, &target)
            .expect("identical sibling configuration");
        sibling_binding.reps += 1;
        assert!(
            validate_manifest_sibling_configuration(&sibling_binding, &sibling, &shared, &target,)
                .is_err(),
            "sibling row-binding reps must equal the target binding"
        );
        let mut changed_receipt = sibling;
        changed_receipt.warmup_runs += 1;
        assert!(
            validate_manifest_sibling_configuration(&shared, &changed_receipt, &shared, &target,)
                .is_err(),
            "sibling receipt warmup must equal the target receipt"
        );
    }

    #[test]
    fn admission_source_baseline_policy_and_day_mutations_fail_closed() {
        let fixture = admission_fixture();
        let source_mutation = fixture
            .admission
            .replacen(&numbered_hash(3), &numbered_hash(4), 1);
        assert!(
            validate_op_ir(
                &fixture_op_ir(&fixture, &source_mutation),
                &fixture.binding,
                &fixture.receipt,
            )
            .is_err()
        );

        let baseline_mutation =
            fixture
                .admission
                .replacen("\"age_policy_days\":90", "\"age_policy_days\":91", 1);
        assert!(
            validate_op_ir(
                &fixture_op_ir(&fixture, &baseline_mutation),
                &fixture.binding,
                &fixture.receipt,
            )
            .is_err()
        );

        let day_mutation =
            fixture
                .admission
                .replacen("\"decision_day\":100", "\"decision_day\":101", 1);
        assert!(
            validate_op_ir(
                &fixture_op_ir(&fixture, &day_mutation),
                &fixture.binding,
                &fixture.receipt,
            )
            .is_err()
        );

        let policy_mutation = fixture
            .admission
            .replacen(&numbered_hash(9), &numbered_hash(8), 1);
        let mutated_receipt = finalized_run_receipt_for_payloads(
            &policy_mutation,
            &fixture.payloads,
            &fixture.manifest,
        )
        .unwrap();
        assert_ne!(fixture.binding.run_receipt, mutated_receipt);
        assert_eq!(wall_ns_day(100 * WALL_NS_PER_DAY as i64).unwrap(), 100);
        assert_ne!(wall_ns_day(101 * WALL_NS_PER_DAY as i64).unwrap(), 100);
    }

    #[test]
    fn finalized_run_v3_binds_every_sibling_payload_and_manifest_byte() {
        let fixture = admission_fixture();
        let payloads = vec!["first-result".to_string(), "second-result".to_string()];
        let first = fs_ledger::hash_bytes(payloads[0].as_bytes());
        let second = fs_ledger::hash_bytes(payloads[1].as_bytes());
        let manifest = format!(
            "{{\"schema\":\"fs-roofline-run-manifest-v1\",\"entries\":[{{\"ordinal\":0,\"kernel\":\"first\",\"version\":\"1\",\"payload\":\"{first}\"}},{{\"ordinal\":1,\"kernel\":\"second\",\"version\":\"1\",\"payload\":\"{second}\"}}]}}"
        );
        let retained =
            finalized_run_receipt_for_payloads(&fixture.admission, &payloads, &manifest).unwrap();

        let mut changed_sibling = payloads.clone();
        changed_sibling[1].push('!');
        assert_ne!(
            retained,
            finalized_run_receipt_for_payloads(&fixture.admission, &changed_sibling, &manifest,)
                .unwrap()
        );
        let changed_manifest = manifest.replacen("\"ordinal\":1", "\"ordinal\":2", 1);
        assert_ne!(
            retained,
            finalized_run_receipt_for_payloads(&fixture.admission, &payloads, &changed_manifest,)
                .unwrap()
        );
    }

    #[test]
    fn roofline_protocol_constants_are_source_pinned() {
        let lib = include_str!("../../fs-roofline/src/lib.rs");
        let baseline = include_str!("../../fs-roofline/src/baseline.rs");
        let axes = include_str!("../../fs-roofline/src/axes.rs");
        let production = include_str!("../../fs-roofline/src/production.rs");
        let kernels = include_str!("../../fs-roofline/src/kernels.rs");
        let authority = include_str!("../../fs-roofline/src/authority.rs");
        for needle in [
            FINALIZED_RUN_DOMAIN,
            RESULT_MANIFEST_DOMAIN,
            "fs-roofline-axis-admission-v2",
            "fs-roofline-run-manifest-v1",
        ] {
            assert!(lib.contains(needle), "fs-roofline lib drifted at {needle}");
        }
        for needle in [
            BASELINE_HASH_DOMAIN,
            "pub const MIN_PROMOTION_RUNS: usize = 3;",
            "pub const MAX_BASELINE_AGE_DAYS: u32 = 365;",
            "pub const BASELINE_LOW_BAND: f64 = 0.70;",
            "pub const BASELINE_HIGH_BAND: f64 = 1.15;",
            "const MAX_BASELINE_LINE_BYTES: usize = 16 * 1024;",
        ] {
            assert!(
                baseline.contains(needle),
                "fs-roofline baseline drifted at {needle}"
            );
        }
        assert!(axes.contains("pub const MAX_AXIS_REPROBE_DRIFT: f64 = 0.25;"));
        for needle in [
            PRODUCTION_AXES_RECEIPT_DOMAIN,
            "pub const MAX_PRODUCTION_ELEMENTS: usize = crate::kernels::MAX_VECTOR_KERNEL_ELEMENTS;",
            "pub const MAX_PRODUCTION_KERNEL_RUNS: usize = 64;",
            "pub const MAX_PRODUCTION_WARMUP: usize = MAX_PRODUCTION_KERNEL_RUNS - 1;",
            "pub const MAX_PRODUCTION_REPS: usize = MAX_PRODUCTION_KERNEL_RUNS;",
            "pub const MAX_PRODUCTION_REGISTRY_FLOPS: u128 = 1 << 39;",
            "pub const MAX_PRODUCTION_REGISTRY_BYTES: u128 = 1 << 33;",
        ] {
            assert!(
                production.contains(needle),
                "fs-roofline production envelope drifted at {needle}"
            );
        }
        for needle in [
            "pub const MAX_VECTOR_KERNEL_ELEMENTS: usize = 1 << 24;",
            "pub const MAX_GEMM_THREADS: usize = 4_096;",
            "n.isqrt().max(256)",
            "ProductionKernelWork::scaled(\"simd-axpy-f64\", n_u128, 2, 24)",
            "ProductionKernelWork::scaled(\"simd-dot-f64\", n_u128, 2, 16)",
            "ProductionKernelWork::scaled(\"simd-sum-f64\", n_u128, 1, 8)",
            "ProductionKernelWork::scaled(\"gemm-f64\", gemm_outputs, gemm_flops_per_output, 24)",
            "pub const GEMM_ROOFLINE_VERSION: &str = \"2\";",
            "let mut registry = default_registry(n)?;",
            "registry.push(Box::new(GemmKernel::new(",
        ] {
            assert!(
                kernels.contains(needle),
                "fs-roofline sealed registry drifted at {needle}"
            );
        }
        let mut prior = None;
        for constructor in [
            "Box::new(AxpyKernel::new(n)?)",
            "Box::new(DotKernel::new(n)?)",
            "Box::new(SumKernel::new(n)?)",
            "registry.push(Box::new(GemmKernel::new(",
        ] {
            let position = kernels
                .find(constructor)
                .unwrap_or_else(|| panic!("missing sealed registry constructor {constructor}"));
            if let Some(prior) = prior {
                assert!(
                    position > prior,
                    "sealed registry constructor order drifted at {constructor}"
                );
            }
            prior = Some(position);
        }
        assert_eq!(MAX_PRODUCTION_ELEMENTS, 1 << 24);
        assert_eq!(MAX_PRODUCTION_KERNEL_RUNS, 64);
        assert_eq!(MAX_PRODUCTION_WARMUP, 63);
        assert_eq!(MAX_PRODUCTION_REPS, 64);
        assert_eq!(MAX_PRODUCTION_REGISTRY_FLOPS, 1 << 39);
        assert_eq!(MAX_PRODUCTION_REGISTRY_BYTES, 1 << 33);
        assert_eq!(MAX_PRODUCTION_LOGICAL_CPUS, 4_096);
        assert!(authority.contains("const MAX_PROMOTION_AUTHORITY_FIELD_BYTES: usize = 4096;"));
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
        assert_eq!(decoded.elements, 1_000);
        assert_eq!(decoded.warmup_runs, 1);
        assert_eq!(decoded.reps, 3);
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
    fn receipt_v3_refuses_values_outside_the_sealed_production_envelope() {
        let receipt = production_receipt();
        for hostile in [
            receipt.replacen(
                "\"elements\":1000",
                &format!("\"elements\":{}", MAX_PRODUCTION_ELEMENTS + 1),
                1,
            ),
            receipt.replacen("\"warmup_runs\":1", "\"warmup_runs\":64", 1),
            receipt.replacen("\"warmup_runs\":1", "\"warmup_runs\":63", 1),
            receipt.replacen("\"reps\":3", "\"reps\":65", 1),
            receipt.replacen(
                "\"logical_cpus\":8",
                &format!("\"logical_cpus\":{}", MAX_PRODUCTION_LOGICAL_CPUS + 1),
                1,
            ),
        ] {
            assert!(
                decode_receipt(&hostile).is_err(),
                "producer-impossible receipt unexpectedly decoded: {hostile}"
            );
        }
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
