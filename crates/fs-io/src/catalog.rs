//! CSV/JSON catalog ingestion with SCHEMA VALIDATION — the AISC-catalog
//! path for the frame flagship. Every cell is checked against a declared
//! column spec; violations are helpful errors naming the row, column,
//! offending text, and the expectation. Quoted CSV fields (RFC-4180
//! subset with escaped quotes) are supported.

use crate::IoError;
use std::collections::BTreeMap;

/// Version of the sealed catalog-schema admission contract.
pub const CATALOG_SCHEMA_VERSION: &str = "fs-io/catalog-schema/v1";

/// Resource envelope for admitting a catalog schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CatalogSchemaLimits {
    /// Maximum number of declared columns.
    pub max_columns: usize,
    /// Maximum UTF-8 bytes in one canonical column name.
    pub max_name_bytes: usize,
    /// Maximum UTF-8 bytes summed over all canonical column names.
    pub max_total_name_bytes: usize,
}

impl CatalogSchemaLimits {
    /// Default schema envelope for [`Schema::admit`].
    pub const DEFAULT: Self = Self {
        max_columns: 4_096,
        max_name_bytes: 256,
        max_total_name_bytes: 64 * 1024,
    };
}

impl Default for CatalogSchemaLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// What a column must contain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColumnKind {
    /// Any nonempty string.
    Text,
    /// A finite float, optionally bounded.
    Number {
        /// Inclusive lower bound.
        min: f64,
        /// Inclusive upper bound.
        max: f64,
    },
}

/// One column's contract.
#[derive(Debug, Clone, PartialEq)]
pub struct ColumnSpec {
    /// Canonical column name. Required columns must appear in every document;
    /// optional columns may be absent.
    pub name: &'static str,
    /// The value contract.
    pub kind: ColumnKind,
    /// Whether empty cells are allowed.
    pub required: bool,
}

/// Deterministic evidence retained by an admitted catalog schema.
///
/// `local_identity_fnv1a64` is a stable replay fingerprint over the version,
/// limits, ordered column contracts, and lookup policies. It is not a
/// collision-resistant content address; HELM must upgrade it before using it
/// as ledger authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CatalogSchemaReceipt {
    /// Versioned identity domain and admission semantics.
    pub schema_version: &'static str,
    /// Caller-selected schema limits used during admission.
    pub limits: CatalogSchemaLimits,
    /// Number of admitted columns.
    pub column_count: usize,
    /// UTF-8 bytes summed over all admitted column names.
    pub total_name_bytes: usize,
    /// Deterministic, non-cryptographic local replay identity.
    pub local_identity_fnv1a64: u64,
}

/// Why an unchecked column declaration could not become a [`Schema`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SchemaDefinitionRefusal {
    /// At least one column is required.
    EmptySchema,
    /// The declaration exceeds the admitted column-count cap.
    ColumnCount {
        /// Supplied columns.
        count: usize,
        /// Admitted maximum.
        limit: usize,
    },
    /// A name is empty after applying the CSV lookup normalization.
    EmptyName {
        /// One-based declaration position.
        column: usize,
    },
    /// A name contains leading or trailing whitespace and would therefore
    /// alias another spelling under CSV header lookup.
    NonCanonicalName {
        /// One-based declaration position.
        column: usize,
    },
    /// One name exceeds the per-name byte cap.
    NameBytes {
        /// One-based declaration position.
        column: usize,
        /// Supplied UTF-8 bytes.
        bytes: usize,
        /// Admitted maximum.
        limit: usize,
    },
    /// Aggregate name bytes exceed the schema envelope.
    TotalNameBytes {
        /// Bytes through the first refusing declaration.
        bytes: usize,
        /// Admitted maximum.
        limit: usize,
    },
    /// Two declarations have the same canonical lookup name.
    DuplicateName {
        /// One-based position of the first declaration.
        first_column: usize,
        /// One-based position of the duplicate declaration.
        duplicate_column: usize,
    },
    /// A numeric lower or upper bound is NaN or infinite.
    NonFiniteNumberBound {
        /// One-based declaration position.
        column: usize,
        /// `true` for the lower bound, `false` for the upper bound.
        lower: bool,
    },
    /// A numeric lower bound is greater than its upper bound.
    InvertedNumberBounds {
        /// One-based declaration position.
        column: usize,
    },
}

impl core::fmt::Display for SchemaDefinitionRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptySchema => write!(f, "catalog schema must declare at least one column"),
            Self::ColumnCount { count, limit } => {
                write!(f, "catalog schema has {count} columns; limit is {limit}")
            }
            Self::EmptyName { column } => {
                write!(f, "catalog schema column {column} has an empty name")
            }
            Self::NonCanonicalName { column } => write!(
                f,
                "catalog schema column {column} has leading or trailing whitespace"
            ),
            Self::NameBytes {
                column,
                bytes,
                limit,
            } => write!(
                f,
                "catalog schema column {column} name has {bytes} bytes; limit is {limit}"
            ),
            Self::TotalNameBytes { bytes, limit } => write!(
                f,
                "catalog schema names total {bytes} bytes; limit is {limit}"
            ),
            Self::DuplicateName {
                first_column,
                duplicate_column,
            } => write!(
                f,
                "catalog schema columns {first_column} and {duplicate_column} have the same name"
            ),
            Self::NonFiniteNumberBound { column, lower } => write!(
                f,
                "catalog schema column {column} has a non-finite {} bound",
                if *lower { "lower" } else { "upper" }
            ),
            Self::InvertedNumberBounds { column } => write!(
                f,
                "catalog schema column {column} has a lower bound greater than its upper bound"
            ),
        }
    }
}

impl std::error::Error for SchemaDefinitionRefusal {}

/// An admitted, immutable catalog schema.
#[derive(Debug, Clone, PartialEq)]
pub struct Schema {
    columns: Vec<ColumnSpec>,
    receipt: CatalogSchemaReceipt,
}

/// A validated catalog: rows of (column name → text) with numbers
/// pre-parsed where the schema demands them.
#[derive(Debug, Clone, PartialEq)]
pub struct Catalog {
    /// Row-major cells keyed by column name.
    pub rows: Vec<BTreeMap<String, String>>,
    /// Pre-parsed numeric views for Number columns.
    pub numbers: Vec<BTreeMap<String, f64>>,
}

/// Shared validation/projection/output envelope for CSV and JSON catalogs.
///
/// Counts are logical retained payload and deterministic work, not allocator
/// metadata or `BTreeMap` node size. Format readers compose this envelope with
/// their syntax-specific input limits before parsing or projecting any row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CatalogProjectionLimits {
    /// Maximum `(row, schema-column)` validation visits.
    pub max_validation_visits: usize,
    /// Maximum numeric entries retained across all projected rows.
    pub max_numeric_entries: usize,
    /// Maximum UTF-8 bytes retained by cloned numeric-projection keys.
    pub max_numeric_key_bytes: usize,
    /// Maximum logical UTF-8 bytes retained by row and numeric maps.
    pub max_output_bytes: usize,
}

impl CatalogProjectionLimits {
    /// Default projection envelope shared by both formats.
    pub const DEFAULT: Self = Self {
        max_validation_visits: 16_000_000,
        max_numeric_entries: 1_000_000,
        max_numeric_key_bytes: 32 * 1024 * 1024,
        max_output_bytes: 256 * 1024 * 1024,
    };
}

impl Default for CatalogProjectionLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Resource envelope for CSV syntax, decoded fields, and shared projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CatalogCsvLimits {
    /// Maximum UTF-8 input bytes, including delimiters and line endings.
    pub max_input_bytes: usize,
    /// Maximum nonempty data records, excluding the header.
    pub max_rows: usize,
    /// Maximum fields in the header or one data record.
    pub max_fields_per_record: usize,
    /// Maximum fields summed over the header and all data records.
    pub max_total_fields: usize,
    /// Maximum decoded UTF-8 bytes in one header or data field.
    pub max_field_bytes: usize,
    /// Maximum decoded UTF-8 bytes summed over raw header fields before
    /// whitespace normalization.
    pub max_header_bytes: usize,
    /// Maximum decoded UTF-8 bytes summed over every field.
    pub max_decoded_bytes: usize,
    /// Shared schema-validation and retained-output envelope.
    pub projection: CatalogProjectionLimits,
}

impl CatalogCsvLimits {
    /// Default bounded CSV world-boundary envelope.
    pub const DEFAULT: Self = Self {
        max_input_bytes: 64 * 1024 * 1024,
        max_rows: 250_000,
        max_fields_per_record: 4_096,
        max_total_fields: 1_000_000,
        max_field_bytes: 1024 * 1024,
        max_header_bytes: 64 * 1024,
        max_decoded_bytes: 32 * 1024 * 1024,
        projection: CatalogProjectionLimits::DEFAULT,
    };
}

impl Default for CatalogCsvLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Resource envelope for the strict catalog-JSON reader.
///
/// The limits count logical payload, not allocator metadata: input bytes include
/// JSON whitespace and delimiters, decoded bytes include every decoded key and
/// value plus every retained number lexeme, and string/number limits apply to
/// one token. All caps are checked before growing an owned payload. `BTreeMap`
/// does not expose a fallible node-reservation API, but its insertions happen
/// only after the row/member/payload caps have admitted the member.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CatalogJsonLimits {
    /// Maximum UTF-8 input bytes, including JSON syntax and whitespace.
    pub max_input_bytes: usize,
    /// Maximum objects in the top-level array.
    pub max_rows: usize,
    /// Maximum members in one object.
    pub max_members_per_object: usize,
    /// Maximum members summed over every object.
    pub max_total_members: usize,
    /// Maximum decoded UTF-8 bytes in one key or string value.
    pub max_string_bytes: usize,
    /// Maximum bytes in one retained JSON number lexeme.
    pub max_number_bytes: usize,
    /// Maximum decoded key/value/number bytes summed over the catalog.
    pub max_decoded_bytes: usize,
    /// Shared schema-validation and retained-output envelope.
    pub projection: CatalogProjectionLimits,
}

impl CatalogJsonLimits {
    /// Default world-boundary envelope for [`Schema::parse_json`].
    pub const DEFAULT: Self = Self {
        max_input_bytes: 64 * 1024 * 1024,
        max_rows: 250_000,
        max_members_per_object: 4_096,
        max_total_members: 1_000_000,
        max_string_bytes: 1024 * 1024,
        max_number_bytes: 256,
        max_decoded_bytes: 32 * 1024 * 1024,
        projection: CatalogProjectionLimits::DEFAULT,
    };
}

impl Default for CatalogJsonLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

#[derive(Debug, Default)]
struct CsvCounters {
    total_fields: usize,
    decoded_bytes: usize,
    header_decoded_bytes: usize,
}

fn csv_cap_refusal(cap: &str, limit: usize, row: usize) -> IoError {
    IoError::ResourceBound {
        what: format!("catalog CSV {cap} cap {limit} exceeded at record {row}"),
    }
}

fn push_csv_char(
    field: &mut String,
    character: char,
    row: usize,
    limits: CatalogCsvLimits,
) -> Result<(), IoError> {
    let width = character.len_utf8();
    let next = field
        .len()
        .checked_add(width)
        .ok_or_else(|| csv_cap_refusal("field decoded-byte", limits.max_field_bytes, row))?;
    if next > limits.max_field_bytes {
        return Err(csv_cap_refusal(
            "field decoded-byte",
            limits.max_field_bytes,
            row,
        ));
    }
    field
        .try_reserve(width)
        .map_err(|_| allocation_refusal("decoded CSV field", row))?;
    field.push(character);
    Ok(())
}

fn finish_csv_field(
    fields: &mut Vec<String>,
    field: &mut String,
    row: usize,
    limits: CatalogCsvLimits,
    counters: &mut CsvCounters,
) -> Result<(), IoError> {
    if fields.len() >= limits.max_fields_per_record {
        return Err(csv_cap_refusal(
            "per-record field",
            limits.max_fields_per_record,
            row,
        ));
    }
    counters.total_fields = counters
        .total_fields
        .checked_add(1)
        .filter(|total| *total <= limits.max_total_fields)
        .ok_or_else(|| csv_cap_refusal("aggregate field", limits.max_total_fields, row))?;
    counters.decoded_bytes = counters
        .decoded_bytes
        .checked_add(field.len())
        .filter(|total| *total <= limits.max_decoded_bytes)
        .ok_or_else(|| csv_cap_refusal("aggregate decoded-byte", limits.max_decoded_bytes, row))?;
    if row == 0 {
        counters.header_decoded_bytes = counters
            .header_decoded_bytes
            .checked_add(field.len())
            .filter(|total| *total <= limits.max_header_bytes)
            .ok_or_else(|| csv_cap_refusal("header decoded-byte", limits.max_header_bytes, row))?;
    }
    fields
        .try_reserve(1)
        .map_err(|_| allocation_refusal("CSV field index", row))?;
    fields.push(core::mem::take(field));
    Ok(())
}

/// Split one bounded CSV record (RFC-4180 subset: quoted fields, `""` escapes).
fn split_csv(
    line: &str,
    row: usize,
    limits: CatalogCsvLimits,
    counters: &mut CsvCounters,
) -> Result<Vec<String>, IoError> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut chars = line.chars().peekable();
    let mut quoted = false;
    while let Some(c) = chars.next() {
        if quoted {
            match c {
                '"' if chars.peek() == Some(&'"') => {
                    push_csv_char(&mut cur, '"', row, limits)?;
                    chars.next();
                }
                '"' => quoted = false,
                other => push_csv_char(&mut cur, other, row, limits)?,
            }
        } else {
            match c {
                '"' if cur.is_empty() => quoted = true,
                ',' => finish_csv_field(&mut fields, &mut cur, row, limits, counters)?,
                other => push_csv_char(&mut cur, other, row, limits)?,
            }
        }
    }
    if quoted {
        return Err(IoError::Malformed {
            at: row,
            what: "unterminated quoted CSV field".to_string(),
        });
    }
    finish_csv_field(&mut fields, &mut cur, row, limits, counters)?;
    Ok(fields)
}

const FNV1A64_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV1A64_PRIME: u64 = 0x0000_0100_0000_01b3;

fn schema_hash_bytes(state: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *state ^= u64::from(*byte);
        *state = state.wrapping_mul(FNV1A64_PRIME);
    }
}

fn schema_hash_usize(state: &mut u64, mut value: usize) {
    loop {
        let mut byte = (value & 0x7f) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        schema_hash_bytes(state, &[byte]);
        if value == 0 {
            return;
        }
    }
}

fn schema_identity(columns: &[ColumnSpec], limits: CatalogSchemaLimits) -> u64 {
    let mut state = FNV1A64_OFFSET;
    schema_hash_bytes(
        &mut state,
        b"fs-io/catalog-schema/v1\0csv-name=trim\0json-name=exact\0unknown=preserve\0optional=may-omit\0validation=declaration-order\0",
    );
    schema_hash_usize(&mut state, limits.max_columns);
    schema_hash_usize(&mut state, limits.max_name_bytes);
    schema_hash_usize(&mut state, limits.max_total_name_bytes);
    schema_hash_usize(&mut state, columns.len());
    for column in columns {
        schema_hash_usize(&mut state, column.name.len());
        schema_hash_bytes(&mut state, column.name.as_bytes());
        match column.kind {
            ColumnKind::Text => schema_hash_bytes(&mut state, &[0]),
            ColumnKind::Number { min, max } => {
                schema_hash_bytes(&mut state, &[1]);
                schema_hash_bytes(&mut state, &min.to_bits().to_le_bytes());
                schema_hash_bytes(&mut state, &max.to_bits().to_le_bytes());
            }
        }
        schema_hash_bytes(&mut state, &[u8::from(column.required)]);
    }
    state
}

const MAX_DIAGNOSTIC_TEXT_BYTES: usize = 96;
const MAX_DIAGNOSTIC_HEADER_NAMES: usize = 8;

fn bounded_diagnostic_text(text: &str) -> String {
    if text.len() <= MAX_DIAGNOSTIC_TEXT_BYTES {
        return text.to_owned();
    }
    let mut end = MAX_DIAGNOSTIC_TEXT_BYTES;
    while !text.is_char_boundary(end) {
        end -= 1;
    }
    format!("{}… ({} UTF-8 bytes)", &text[..end], text.len())
}

fn header_witness(header: &[String]) -> String {
    let mut witness = String::new();
    for (index, name) in header.iter().take(MAX_DIAGNOSTIC_HEADER_NAMES).enumerate() {
        if index != 0 {
            witness.push_str(", ");
        }
        witness.push_str(&bounded_diagnostic_text(name));
    }
    if header.len() > MAX_DIAGNOSTIC_HEADER_NAMES {
        witness.push_str(&format!(
            ", … ({} more columns)",
            header.len() - MAX_DIAGNOSTIC_HEADER_NAMES
        ));
    }
    witness
}

fn fallible_copy(text: &str, payload: &str, at: usize) -> Result<String, IoError> {
    let mut copy = String::new();
    copy.try_reserve_exact(text.len())
        .map_err(|_| allocation_refusal(payload, at))?;
    copy.push_str(text);
    Ok(copy)
}

fn normalize_csv_header(raw_header: Vec<String>) -> Result<Vec<String>, IoError> {
    let mut header = Vec::new();
    header
        .try_reserve_exact(raw_header.len())
        .map_err(|_| allocation_refusal("normalized CSV header", 0))?;
    for (index, raw_name) in raw_header.into_iter().enumerate() {
        let name = raw_name.trim();
        if name.is_empty() {
            return Err(IoError::Schema {
                row: 0,
                column: format!("header column {}", index + 1),
                what: "CSV header name is empty after whitespace normalization".to_string(),
            });
        }
        header.push(fallible_copy(name, "normalized CSV header name", 0)?);
    }

    let mut first_positions = BTreeMap::<&str, usize>::new();
    for (index, name) in header.iter().enumerate() {
        if let Some(first_index) = first_positions.insert(name, index) {
            return Err(IoError::Schema {
                row: 0,
                column: bounded_diagnostic_text(name),
                what: format!(
                    "duplicate CSV header after whitespace normalization at columns {} and {}",
                    first_index + 1,
                    index + 1
                ),
            });
        }
    }
    Ok(header)
}

fn operation_refusal(format: &str, detail: impl Into<String>) -> IoError {
    IoError::ResourceBound {
        what: format!(
            "catalog {format} operation admission refused: {}",
            detail.into()
        ),
    }
}

fn checked_operation_add(
    left: usize,
    right: usize,
    format: &str,
    field: &str,
) -> Result<usize, IoError> {
    left.checked_add(right).ok_or_else(|| {
        operation_refusal(
            format,
            format!("{field} arithmetic overflow in checked addition"),
        )
    })
}

fn checked_operation_product(
    left: usize,
    right: usize,
    format: &str,
    field: &str,
) -> Result<usize, IoError> {
    left.checked_mul(right).ok_or_else(|| {
        operation_refusal(
            format,
            format!("{field} arithmetic overflow in checked multiplication"),
        )
    })
}

fn ensure_operation_cap(
    actual: usize,
    limit: usize,
    format: &str,
    cap: &str,
) -> Result<(), IoError> {
    if actual <= limit {
        Ok(())
    } else {
        Err(operation_refusal(
            format,
            format!("{cap} requires {actual}; declared limit is {limit}"),
        ))
    }
}

#[derive(Debug, Clone, Copy)]
struct ProjectionPlan {
    max_rows: usize,
    numeric_key_bytes: usize,
}

#[derive(Debug)]
struct ProjectionUse {
    format: &'static str,
    limits: CatalogProjectionLimits,
    validation_visits: usize,
    numeric_entries: usize,
    numeric_key_bytes: usize,
    output_bytes: usize,
}

impl ProjectionUse {
    fn new(
        format: &'static str,
        limits: CatalogProjectionLimits,
        initial_output_bytes: usize,
    ) -> Result<Self, IoError> {
        ensure_operation_cap(
            initial_output_bytes,
            limits.max_output_bytes,
            format,
            "logical output bytes",
        )?;
        Ok(Self {
            format,
            limits,
            validation_visits: 0,
            numeric_entries: 0,
            numeric_key_bytes: 0,
            output_bytes: initial_output_bytes,
        })
    }

    fn charge(
        current: &mut usize,
        amount: usize,
        limit: usize,
        format: &str,
        cap: &str,
    ) -> Result<(), IoError> {
        let next = checked_operation_add(*current, amount, format, cap)?;
        ensure_operation_cap(next, limit, format, cap)?;
        *current = next;
        Ok(())
    }

    fn validation_visit(&mut self) -> Result<(), IoError> {
        Self::charge(
            &mut self.validation_visits,
            1,
            self.limits.max_validation_visits,
            self.format,
            "schema validation visits",
        )
    }

    fn output(&mut self, bytes: usize) -> Result<(), IoError> {
        Self::charge(
            &mut self.output_bytes,
            bytes,
            self.limits.max_output_bytes,
            self.format,
            "logical output bytes",
        )
    }

    fn numeric_entry(&mut self, key_bytes: usize) -> Result<(), IoError> {
        Self::charge(
            &mut self.numeric_entries,
            1,
            self.limits.max_numeric_entries,
            self.format,
            "numeric projection entries",
        )?;
        Self::charge(
            &mut self.numeric_key_bytes,
            key_bytes,
            self.limits.max_numeric_key_bytes,
            self.format,
            "numeric projection key bytes",
        )?;
        self.output(key_bytes)
    }
}

impl Schema {
    /// Admit a schema under [`CatalogSchemaLimits::DEFAULT`].
    ///
    /// Declaration order fixes deterministic validation-error priority and is
    /// therefore identity-bearing. Document column/member order is not.
    ///
    /// # Errors
    /// Returns [`SchemaDefinitionRefusal`] before any schema can be used when
    /// the declaration is empty, ambiguous, out of bounds, or has invalid
    /// numeric bounds.
    pub fn admit(columns: Vec<ColumnSpec>) -> Result<Self, SchemaDefinitionRefusal> {
        Self::admit_with_limits(columns, CatalogSchemaLimits::DEFAULT)
    }

    /// Admit a schema under caller-explicit definition limits.
    ///
    /// # Errors
    /// Returns the first refusal in declaration order. Names must already be
    /// in their `str::trim` canonical form so CSV and JSON lookup cannot
    /// disagree about aliases.
    pub fn admit_with_limits(
        columns: Vec<ColumnSpec>,
        limits: CatalogSchemaLimits,
    ) -> Result<Self, SchemaDefinitionRefusal> {
        if columns.is_empty() {
            return Err(SchemaDefinitionRefusal::EmptySchema);
        }
        if columns.len() > limits.max_columns {
            return Err(SchemaDefinitionRefusal::ColumnCount {
                count: columns.len(),
                limit: limits.max_columns,
            });
        }

        let mut total_name_bytes = 0usize;
        let mut names = BTreeMap::<&str, usize>::new();
        for (index, column) in columns.iter().enumerate() {
            let ordinal = index + 1;
            let canonical = column.name.trim();
            if canonical.is_empty() {
                return Err(SchemaDefinitionRefusal::EmptyName { column: ordinal });
            }
            if canonical != column.name {
                return Err(SchemaDefinitionRefusal::NonCanonicalName { column: ordinal });
            }
            if column.name.len() > limits.max_name_bytes {
                return Err(SchemaDefinitionRefusal::NameBytes {
                    column: ordinal,
                    bytes: column.name.len(),
                    limit: limits.max_name_bytes,
                });
            }
            let next_total = total_name_bytes.checked_add(column.name.len()).ok_or(
                SchemaDefinitionRefusal::TotalNameBytes {
                    bytes: usize::MAX,
                    limit: limits.max_total_name_bytes,
                },
            )?;
            if next_total > limits.max_total_name_bytes {
                return Err(SchemaDefinitionRefusal::TotalNameBytes {
                    bytes: next_total,
                    limit: limits.max_total_name_bytes,
                });
            }
            total_name_bytes = next_total;
            if let Some(first_index) = names.insert(column.name, index) {
                return Err(SchemaDefinitionRefusal::DuplicateName {
                    first_column: first_index + 1,
                    duplicate_column: ordinal,
                });
            }
            if let ColumnKind::Number { min, max } = column.kind {
                if !min.is_finite() {
                    return Err(SchemaDefinitionRefusal::NonFiniteNumberBound {
                        column: ordinal,
                        lower: true,
                    });
                }
                if !max.is_finite() {
                    return Err(SchemaDefinitionRefusal::NonFiniteNumberBound {
                        column: ordinal,
                        lower: false,
                    });
                }
                if min > max {
                    return Err(SchemaDefinitionRefusal::InvertedNumberBounds { column: ordinal });
                }
            }
        }

        let receipt = CatalogSchemaReceipt {
            schema_version: CATALOG_SCHEMA_VERSION,
            limits,
            column_count: columns.len(),
            total_name_bytes,
            local_identity_fnv1a64: schema_identity(&columns, limits),
        };
        Ok(Self { columns, receipt })
    }

    /// Admitted column contracts in deterministic validation order.
    #[must_use]
    pub fn columns(&self) -> &[ColumnSpec] {
        &self.columns
    }

    /// Versioned deterministic schema-admission evidence.
    #[must_use]
    pub const fn receipt(&self) -> &CatalogSchemaReceipt {
        &self.receipt
    }

    fn preflight_projection(
        &self,
        format: &'static str,
        max_rows: usize,
        max_fields_per_row: usize,
        max_total_fields: usize,
        csv_header: Option<&[String]>,
        limits: CatalogProjectionLimits,
    ) -> Result<ProjectionPlan, IoError> {
        let required_columns = self.columns.iter().filter(|column| column.required).count();
        ensure_operation_cap(
            required_columns,
            max_fields_per_row,
            format,
            "required schema columns per row",
        )?;

        let header_fields = csv_header.map_or(0, <[String]>::len);
        ensure_operation_cap(
            header_fields,
            max_total_fields,
            format,
            "minimum header fields",
        )?;
        let data_field_budget = max_total_fields - header_fields;
        // A CSV record always has the admitted header width, so its aggregate
        // field cap also bounds row count. JSON admits empty objects; hostile
        // or schema-invalid rows can therefore still consume the full row and
        // validation budget without consuming a member.
        let minimum_fields_per_row = csv_header.map_or(0, <[String]>::len);
        let jointly_admitted_rows = if minimum_fields_per_row == 0 {
            max_rows
        } else {
            max_rows.min(data_field_budget / minimum_fields_per_row)
        };

        let validation_visits = checked_operation_product(
            self.columns.len(),
            jointly_admitted_rows,
            format,
            "schema validation visits",
        )?;
        ensure_operation_cap(
            validation_visits,
            limits.max_validation_visits,
            format,
            "schema validation visits",
        )?;

        let mut numeric_columns = 0usize;
        let mut numeric_name_bytes = 0usize;
        for column in &self.columns {
            let present =
                csv_header.is_none_or(|header| header.iter().any(|name| name == column.name));
            if present && matches!(column.kind, ColumnKind::Number { .. }) {
                numeric_columns =
                    checked_operation_add(numeric_columns, 1, format, "numeric schema columns")?;
                numeric_name_bytes = checked_operation_add(
                    numeric_name_bytes,
                    column.name.len(),
                    format,
                    "numeric schema-name bytes",
                )?;
            }
        }
        let numeric_entries = checked_operation_product(
            numeric_columns.min(max_fields_per_row),
            jointly_admitted_rows,
            format,
            "numeric projection entries",
        )?
        .min(data_field_budget);
        ensure_operation_cap(
            numeric_entries,
            limits.max_numeric_entries,
            format,
            "numeric projection entries",
        )?;
        let numeric_key_bytes = checked_operation_product(
            numeric_name_bytes,
            jointly_admitted_rows,
            format,
            "numeric projection key bytes",
        )?;
        ensure_operation_cap(
            numeric_key_bytes,
            limits.max_numeric_key_bytes,
            format,
            "numeric projection key bytes",
        )?;
        Ok(ProjectionPlan {
            max_rows: jointly_admitted_rows,
            numeric_key_bytes,
        })
    }

    fn preflight_json_output(
        max_decoded_bytes: usize,
        plan: ProjectionPlan,
        limits: CatalogProjectionLimits,
    ) -> Result<(), IoError> {
        let output_bytes = checked_operation_add(
            max_decoded_bytes,
            plan.numeric_key_bytes,
            "JSON",
            "logical output bytes",
        )?;
        ensure_operation_cap(
            output_bytes,
            limits.max_output_bytes,
            "JSON",
            "logical output bytes",
        )
    }

    fn preflight_csv_output(
        header: &[String],
        counters: &CsvCounters,
        limits: CatalogCsvLimits,
        plan: ProjectionPlan,
    ) -> Result<(), IoError> {
        let mut normalized_header_bytes = 0usize;
        for name in header {
            normalized_header_bytes = checked_operation_add(
                normalized_header_bytes,
                name.len(),
                "CSV",
                "normalized header key bytes",
            )?;
        }
        let repeated_header_bytes = checked_operation_product(
            normalized_header_bytes,
            plan.max_rows,
            "CSV",
            "repeated output header-key bytes",
        )?;
        let decoded_value_bytes = limits
            .max_decoded_bytes
            .checked_sub(counters.header_decoded_bytes)
            .ok_or_else(|| {
                operation_refusal(
                    "CSV",
                    "header-byte accounting exceeds the aggregate decoded-byte envelope",
                )
            })?;
        let output_bytes = checked_operation_add(
            decoded_value_bytes,
            repeated_header_bytes,
            "CSV",
            "logical row output bytes",
        )?;
        let output_bytes = checked_operation_add(
            output_bytes,
            plan.numeric_key_bytes,
            "CSV",
            "logical output bytes",
        )?;
        ensure_operation_cap(
            output_bytes,
            limits.projection.max_output_bytes,
            "CSV",
            "logical output bytes",
        )
    }

    /// Validate one cell against its spec.
    fn check_cell(spec: &ColumnSpec, text: &str, row: usize) -> Result<Option<f64>, IoError> {
        let trimmed = text.trim();
        if trimmed.is_empty() {
            if spec.required {
                return Err(IoError::Schema {
                    row,
                    column: spec.name.to_string(),
                    what: "required cell is empty".to_string(),
                });
            }
            return Ok(None);
        }
        match spec.kind {
            ColumnKind::Text => Ok(None),
            ColumnKind::Number { min, max } => {
                let v: f64 = trimmed.parse().map_err(|_| IoError::Schema {
                    row,
                    column: spec.name.to_string(),
                    what: format!("{:?} is not a number", bounded_diagnostic_text(trimmed)),
                })?;
                if !v.is_finite() || v < min || v > max {
                    return Err(IoError::Schema {
                        row,
                        column: spec.name.to_string(),
                        what: format!("{v} outside the declared range [{min}, {max}]"),
                    });
                }
                Ok(Some(v))
            }
        }
    }

    /// Parse + validate a CSV catalog (first record is the header).
    ///
    /// # Errors
    /// [`IoError::Malformed`] for CSV structure; [`IoError::Schema`] with
    /// row/column/expectation for value violations; missing schema
    /// columns are named.
    pub fn parse_csv(&self, text: &str) -> Result<Catalog, IoError> {
        self.parse_csv_with_limits(text, CatalogCsvLimits::DEFAULT)
    }

    /// Parse + validate a CSV catalog under one syntax/decoded/projection
    /// envelope shared with the JSON projection path.
    ///
    /// # Errors
    /// [`IoError::ResourceBound`] is returned before unadmitted field growth,
    /// schema work, numeric projection, or output retention. Syntax and schema
    /// mismatches retain their existing structured variants.
    pub fn parse_csv_with_limits(
        &self,
        text: &str,
        limits: CatalogCsvLimits,
    ) -> Result<Catalog, IoError> {
        if text.len() > limits.max_input_bytes {
            return Err(csv_cap_refusal("input-byte", limits.max_input_bytes, 0));
        }
        let required_columns = self.columns.iter().filter(|column| column.required).count();
        ensure_operation_cap(
            required_columns,
            limits.max_fields_per_record,
            "CSV",
            "required schema columns per row",
        )?;
        ensure_operation_cap(
            required_columns.max(1),
            limits.max_total_fields,
            "CSV",
            "minimum header fields",
        )?;
        let mut counters = CsvCounters::default();
        let mut lines = text
            .lines()
            .enumerate()
            .filter(|(_, l)| !l.trim().is_empty());
        let (_, header_line) = lines.next().ok_or(IoError::Malformed {
            at: 0,
            what: "empty catalog".to_string(),
        })?;
        let header = normalize_csv_header(split_csv(header_line, 0, limits, &mut counters)?)?;
        for spec in &self.columns {
            if spec.required && !header.iter().any(|name| name == spec.name) {
                return Err(IoError::Schema {
                    row: 0,
                    column: spec.name.to_string(),
                    what: format!(
                        "column missing from the header (found: {})",
                        header_witness(&header)
                    ),
                });
            }
        }
        let plan = self.preflight_projection(
            "CSV",
            limits.max_rows,
            limits.max_fields_per_record,
            limits.max_total_fields,
            Some(&header),
            limits.projection,
        )?;
        Self::preflight_csv_output(&header, &counters, limits, plan)?;
        let mut rows = Vec::new();
        let mut numbers = Vec::new();
        let mut projection = ProjectionUse::new("CSV", limits.projection, 0)?;
        for (data_row, (ln, line)) in lines.enumerate() {
            let row_no = data_row + 1; // 1-based, header excluded
            if data_row >= limits.max_rows {
                return Err(csv_cap_refusal("row", limits.max_rows, row_no));
            }
            let fields = split_csv(line, row_no, limits, &mut counters)?;
            if fields.len() != header.len() {
                return Err(IoError::Malformed {
                    at: ln + 1,
                    what: format!(
                        "record has {} fields, header has {}",
                        fields.len(),
                        header.len()
                    ),
                });
            }
            rows.try_reserve(1)
                .map_err(|_| allocation_refusal("CSV output row index", row_no))?;
            numbers
                .try_reserve(1)
                .map_err(|_| allocation_refusal("CSV numeric-row index", row_no))?;
            let mut row = BTreeMap::new();
            let mut nums = BTreeMap::new();
            for (name, cell) in header.iter().zip(fields) {
                let retained_bytes =
                    checked_operation_add(name.len(), cell.len(), "CSV", "row key/value bytes")?;
                projection.output(retained_bytes)?;
                row.insert(fallible_copy(name, "CSV output column key", row_no)?, cell);
            }
            for spec in &self.columns {
                projection.validation_visit()?;
                let cell = row.get(spec.name).map(String::as_str).unwrap_or_default();
                if let Some(v) = Self::check_cell(spec, cell, row_no)? {
                    projection.numeric_entry(spec.name.len())?;
                    nums.insert(
                        fallible_copy(spec.name, "CSV numeric projection key", row_no)?,
                        v,
                    );
                }
            }
            rows.push(row);
            numbers.push(nums);
        }
        Ok(Catalog { rows, numbers })
    }

    /// Parse + validate a JSON catalog: an array of flat objects
    /// (string/number members). The bounded in-house reader implements strict
    /// RFC 8259 grammar and rejects anything outside that declared subset.
    ///
    /// # Errors
    /// [`IoError`] for JSON structure or schema violations.
    pub fn parse_json(&self, text: &str) -> Result<Catalog, IoError> {
        self.parse_json_with_limits(text, CatalogJsonLimits::DEFAULT)
    }

    /// Parse + validate a JSON catalog under a caller-explicit resource
    /// envelope. The accepted language is RFC 8259 JSON restricted to one
    /// top-level array of flat objects whose values are strings or numbers.
    ///
    /// # Errors
    /// [`IoError::Malformed`] identifies the first invalid byte offset;
    /// [`IoError::ResourceBound`] names the cap, limit, and refusal offset;
    /// [`IoError::Schema`] reports row/column value violations.
    pub fn parse_json_with_limits(
        &self,
        text: &str,
        limits: CatalogJsonLimits,
    ) -> Result<Catalog, IoError> {
        if text.len() > limits.max_input_bytes {
            return Err(cap_refusal(
                "input-byte",
                limits.max_input_bytes,
                limits.max_input_bytes,
            ));
        }
        let plan = self.preflight_projection(
            "JSON",
            limits.max_rows,
            limits.max_members_per_object,
            limits.max_total_members,
            None,
            limits.projection,
        )?;
        Self::preflight_json_output(limits.max_decoded_bytes, plan, limits.projection)?;
        let parsed = parse_json_catalog_document(text, limits)?;
        let mut projection = ProjectionUse::new("JSON", limits.projection, parsed.decoded_bytes)?;
        let rows = parsed.rows;
        let mut numbers = Vec::new();
        numbers
            .try_reserve_exact(rows.len())
            .map_err(|_| allocation_refusal("catalog numeric-row index", 0))?;
        for (i, obj) in rows.iter().enumerate() {
            let mut nums = BTreeMap::new();
            for spec in &self.columns {
                projection.validation_visit()?;
                let cell = obj.get(spec.name).map(String::as_str).unwrap_or_default();
                if let Some(v) = Self::check_cell(spec, cell, i + 1)? {
                    projection.numeric_entry(spec.name.len())?;
                    nums.insert(
                        fallible_copy(spec.name, "JSON numeric projection key", i + 1)?,
                        v,
                    );
                }
            }
            numbers.push(nums);
        }
        Ok(Catalog { rows, numbers })
    }
}

fn malformed(at: usize, what: impl Into<String>) -> IoError {
    IoError::Malformed {
        at,
        what: what.into(),
    }
}

fn cap_refusal(cap: &str, limit: usize, at: usize) -> IoError {
    IoError::ResourceBound {
        what: format!("catalog JSON {cap} cap {limit} exceeded at byte offset {at}"),
    }
}

fn allocation_refusal(payload: &str, at: usize) -> IoError {
    IoError::ResourceBound {
        what: format!("allocation failed for {payload} at deterministic position {at}"),
    }
}

/// Strict RFC 8259 reader for `[ {"k": "v" | number, ...}, ... ]`.
#[cfg(test)]
fn mini_json_array_of_objects(
    text: &str,
    limits: CatalogJsonLimits,
) -> Result<Vec<BTreeMap<String, String>>, IoError> {
    parse_json_catalog_document(text, limits).map(|parsed| parsed.rows)
}

struct ParsedJsonCatalog {
    rows: Vec<BTreeMap<String, String>>,
    decoded_bytes: usize,
}

fn parse_json_catalog_document(
    text: &str,
    limits: CatalogJsonLimits,
) -> Result<ParsedJsonCatalog, IoError> {
    if text.len() > limits.max_input_bytes {
        return Err(cap_refusal(
            "input-byte",
            limits.max_input_bytes,
            limits.max_input_bytes,
        ));
    }
    JsonCatalogParser {
        bytes: text.as_bytes(),
        pos: 0,
        limits,
        total_members: 0,
        decoded_bytes: 0,
    }
    .parse()
}

struct JsonCatalogParser<'a> {
    bytes: &'a [u8],
    pos: usize,
    limits: CatalogJsonLimits,
    total_members: usize,
    decoded_bytes: usize,
}

impl JsonCatalogParser<'_> {
    fn parse(mut self) -> Result<ParsedJsonCatalog, IoError> {
        self.skip_ws();
        self.expect(b'[', "expected a JSON array")?;
        self.skip_ws();

        let mut rows = Vec::new();
        if self.peek() == Some(b']') {
            self.pos += 1;
            self.finish_document()?;
            return Ok(ParsedJsonCatalog {
                rows,
                decoded_bytes: self.decoded_bytes,
            });
        }

        loop {
            if self.peek() != Some(b'{') {
                return Err(malformed(self.pos, "expected a JSON object"));
            }
            if rows.len() >= self.limits.max_rows {
                return Err(cap_refusal("row", self.limits.max_rows, self.pos));
            }
            rows.try_reserve(1)
                .map_err(|_| allocation_refusal("catalog row index", self.pos))?;
            rows.push(self.parse_object()?);
            self.skip_ws();
            match self.peek() {
                Some(b',') => {
                    self.pos += 1;
                    self.skip_ws();
                    if self.peek() == Some(b']') {
                        return Err(malformed(self.pos, "trailing comma in JSON array"));
                    }
                }
                Some(b']') => {
                    self.pos += 1;
                    self.finish_document()?;
                    return Ok(ParsedJsonCatalog {
                        rows,
                        decoded_bytes: self.decoded_bytes,
                    });
                }
                _ => {
                    return Err(malformed(self.pos, "expected ',' or ']' after JSON object"));
                }
            }
        }
    }

    fn parse_object(&mut self) -> Result<BTreeMap<String, String>, IoError> {
        self.expect(b'{', "expected a JSON object")?;
        self.skip_ws();
        let mut object = BTreeMap::new();
        let mut object_members = 0usize;
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(object);
        }

        loop {
            let key_at = self.pos;
            if self.peek() != Some(b'"') {
                return Err(malformed(self.pos, "expected a quoted JSON object key"));
            }
            if object_members >= self.limits.max_members_per_object {
                return Err(cap_refusal(
                    "per-object member",
                    self.limits.max_members_per_object,
                    self.pos,
                ));
            }
            if self.total_members >= self.limits.max_total_members {
                return Err(cap_refusal(
                    "aggregate member",
                    self.limits.max_total_members,
                    self.pos,
                ));
            }
            let key = self.read_string()?;
            if object.contains_key(&key) {
                return Err(malformed(
                    key_at,
                    format!(
                        "duplicate JSON object key ({} decoded UTF-8 bytes)",
                        key.len()
                    ),
                ));
            }
            self.charge_decoded(key.len(), key_at)?;

            self.skip_ws();
            self.expect(b':', "expected ':' after JSON object key")?;
            self.skip_ws();
            let value_at = self.pos;
            let value = match self.peek() {
                Some(b'"') => self.read_string()?,
                Some(b'-' | b'0'..=b'9') => self.read_number()?,
                _ => {
                    return Err(malformed(
                        self.pos,
                        "expected a JSON string or number value",
                    ));
                }
            };
            self.charge_decoded(value.len(), value_at)?;

            object.insert(key, value);
            object_members += 1;
            self.total_members += 1;

            self.skip_ws();
            match self.peek() {
                Some(b',') => {
                    self.pos += 1;
                    self.skip_ws();
                    if self.peek() == Some(b'}') {
                        return Err(malformed(self.pos, "trailing comma in JSON object"));
                    }
                }
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(object);
                }
                _ => {
                    return Err(malformed(
                        self.pos,
                        "expected ',' or '}' after JSON object member",
                    ));
                }
            }
        }
    }

    fn read_string(&mut self) -> Result<String, IoError> {
        self.expect(b'"', "expected JSON string")?;
        let mut output = String::new();
        loop {
            let chunk_start = self.pos;
            let string_remaining = self.limits.max_string_bytes.saturating_sub(output.len());
            let aggregate_prefix =
                self.decoded_bytes
                    .checked_add(output.len())
                    .ok_or_else(|| {
                        cap_refusal(
                            "aggregate decoded-byte",
                            self.limits.max_decoded_bytes,
                            self.pos,
                        )
                    })?;
            let aggregate_remaining = self
                .limits
                .max_decoded_bytes
                .saturating_sub(aggregate_prefix);
            let raw_remaining = string_remaining.min(aggregate_remaining);
            while let Some(byte) = self.peek() {
                if byte == b'"' || byte == b'\\' || byte < 0x20 {
                    break;
                }
                if self.pos - chunk_start >= raw_remaining {
                    let (cap, limit) = if string_remaining <= aggregate_remaining {
                        ("string decoded-byte", self.limits.max_string_bytes)
                    } else {
                        ("aggregate decoded-byte", self.limits.max_decoded_bytes)
                    };
                    return Err(cap_refusal(cap, limit, self.pos));
                }
                self.pos += 1;
            }
            if self.pos > chunk_start {
                let chunk = core::str::from_utf8(&self.bytes[chunk_start..self.pos])
                    .map_err(|_| malformed(chunk_start, "invalid UTF-8 in JSON string"))?;
                self.append_string(&mut output, chunk, chunk_start)?;
            }

            match self.peek() {
                None => return Err(malformed(self.pos, "unterminated JSON string")),
                Some(b'"') => {
                    self.pos += 1;
                    return Ok(output);
                }
                Some(byte) if byte < 0x20 => {
                    return Err(malformed(self.pos, "raw C0 control byte in JSON string"));
                }
                Some(b'\\') => {
                    let escape_at = self.pos;
                    let escaped = *self
                        .bytes
                        .get(self.pos + 1)
                        .ok_or_else(|| malformed(escape_at, "dangling JSON string escape"))?;
                    match escaped {
                        b'"' => {
                            self.pos += 2;
                            self.append_char(&mut output, '"', escape_at)?;
                        }
                        b'\\' => {
                            self.pos += 2;
                            self.append_char(&mut output, '\\', escape_at)?;
                        }
                        b'/' => {
                            self.pos += 2;
                            self.append_char(&mut output, '/', escape_at)?;
                        }
                        b'b' => {
                            self.pos += 2;
                            self.append_char(&mut output, '\u{0008}', escape_at)?;
                        }
                        b'f' => {
                            self.pos += 2;
                            self.append_char(&mut output, '\u{000c}', escape_at)?;
                        }
                        b'n' => {
                            self.pos += 2;
                            self.append_char(&mut output, '\n', escape_at)?;
                        }
                        b'r' => {
                            self.pos += 2;
                            self.append_char(&mut output, '\r', escape_at)?;
                        }
                        b't' => {
                            self.pos += 2;
                            self.append_char(&mut output, '\t', escape_at)?;
                        }
                        b'u' => {
                            let first = self.read_hex_quad()?;
                            let scalar = if (0xd800..=0xdbff).contains(&first) {
                                let second_at = self.pos;
                                if self.bytes.get(self.pos..self.pos.saturating_add(2))
                                    != Some(&b"\\u"[..])
                                {
                                    return Err(malformed(
                                        second_at,
                                        "high surrogate must be followed by a low-surrogate escape",
                                    ));
                                }
                                let second = self.read_hex_quad()?;
                                if !(0xdc00..=0xdfff).contains(&second) {
                                    return Err(malformed(
                                        second_at,
                                        "high surrogate followed by a non-low surrogate",
                                    ));
                                }
                                0x1_0000
                                    + (((u32::from(first) - 0xd800) << 10)
                                        | (u32::from(second) - 0xdc00))
                            } else if (0xdc00..=0xdfff).contains(&first) {
                                return Err(malformed(
                                    escape_at,
                                    "unpaired low surrogate in JSON string",
                                ));
                            } else {
                                u32::from(first)
                            };
                            let character = char::from_u32(scalar).ok_or_else(|| {
                                malformed(escape_at, "invalid Unicode scalar in JSON string")
                            })?;
                            self.append_char(&mut output, character, escape_at)?;
                        }
                        _ => {
                            return Err(malformed(
                                self.pos + 1,
                                format!("unknown JSON string escape byte 0x{escaped:02x}"),
                            ));
                        }
                    }
                }
                Some(_) => {
                    return Err(malformed(
                        self.pos,
                        "invalid byte reached the JSON string decoder",
                    ));
                }
            }
        }
    }

    fn read_hex_quad(&mut self) -> Result<u16, IoError> {
        if self.peek() != Some(b'\\') || self.bytes.get(self.pos + 1) != Some(&b'u') {
            return Err(malformed(self.pos, "expected a Unicode escape"));
        }
        self.pos += 2;
        let mut value = 0u16;
        for _ in 0..4 {
            let at = self.pos;
            let byte = *self
                .bytes
                .get(self.pos)
                .ok_or_else(|| malformed(at, "truncated four-digit Unicode escape"))?;
            let digit = match byte {
                b'0'..=b'9' => u16::from(byte - b'0'),
                b'a'..=b'f' => u16::from(byte - b'a' + 10),
                b'A'..=b'F' => u16::from(byte - b'A' + 10),
                _ => return Err(malformed(at, "non-hex digit in Unicode escape")),
            };
            value = (value << 4) | digit;
            self.pos += 1;
        }
        Ok(value)
    }

    fn read_number(&mut self) -> Result<String, IoError> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.advance_number_byte(start)?;
        }
        match self.peek() {
            Some(b'0') => {
                self.advance_number_byte(start)?;
                if matches!(self.peek(), Some(b'0'..=b'9')) {
                    return Err(malformed(
                        self.pos,
                        "leading zero in JSON number integer part",
                    ));
                }
            }
            Some(b'1'..=b'9') => {
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.advance_number_byte(start)?;
                }
            }
            _ => {
                return Err(malformed(
                    self.pos,
                    "expected a digit in JSON number integer part",
                ));
            }
        }
        if self.peek() == Some(b'.') {
            self.advance_number_byte(start)?;
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err(malformed(
                    self.pos,
                    "expected a digit after JSON number decimal point",
                ));
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.advance_number_byte(start)?;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.advance_number_byte(start)?;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.advance_number_byte(start)?;
            }
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err(malformed(
                    self.pos,
                    "expected a digit in JSON number exponent",
                ));
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.advance_number_byte(start)?;
            }
        }
        if let Some(byte) = self.peek()
            && !matches!(byte, b' ' | b'\n' | b'\r' | b'\t' | b',' | b'}')
        {
            return Err(malformed(self.pos, "invalid byte after JSON number"));
        }

        let token = &self.bytes[start..self.pos];
        if self
            .decoded_bytes
            .checked_add(token.len())
            .is_none_or(|total| total > self.limits.max_decoded_bytes)
        {
            return Err(cap_refusal(
                "aggregate decoded-byte",
                self.limits.max_decoded_bytes,
                start,
            ));
        }
        let mut output = String::new();
        output
            .try_reserve_exact(token.len())
            .map_err(|_| allocation_refusal("JSON number token", start))?;
        let token = core::str::from_utf8(token)
            .map_err(|_| malformed(start, "non-ASCII byte in JSON number"))?;
        output.push_str(token);
        Ok(output)
    }

    fn advance_number_byte(&mut self, start: usize) -> Result<(), IoError> {
        if self.pos - start >= self.limits.max_number_bytes {
            return Err(cap_refusal(
                "number-token byte",
                self.limits.max_number_bytes,
                self.pos,
            ));
        }
        self.pos += 1;
        Ok(())
    }

    fn append_char(&self, output: &mut String, character: char, at: usize) -> Result<(), IoError> {
        let mut encoded = [0u8; 4];
        self.append_string(output, character.encode_utf8(&mut encoded), at)
    }

    fn append_string(&self, output: &mut String, text: &str, at: usize) -> Result<(), IoError> {
        let new_len = output
            .len()
            .checked_add(text.len())
            .ok_or_else(|| cap_refusal("string decoded-byte", self.limits.max_string_bytes, at))?;
        if new_len > self.limits.max_string_bytes {
            return Err(cap_refusal(
                "string decoded-byte",
                self.limits.max_string_bytes,
                at,
            ));
        }
        if self
            .decoded_bytes
            .checked_add(new_len)
            .is_none_or(|total| total > self.limits.max_decoded_bytes)
        {
            return Err(cap_refusal(
                "aggregate decoded-byte",
                self.limits.max_decoded_bytes,
                at,
            ));
        }
        output
            .try_reserve(text.len())
            .map_err(|_| allocation_refusal("decoded JSON string", at))?;
        output.push_str(text);
        Ok(())
    }

    fn charge_decoded(&mut self, amount: usize, at: usize) -> Result<(), IoError> {
        self.decoded_bytes = self
            .decoded_bytes
            .checked_add(amount)
            .filter(|total| *total <= self.limits.max_decoded_bytes)
            .ok_or_else(|| {
                cap_refusal("aggregate decoded-byte", self.limits.max_decoded_bytes, at)
            })?;
        Ok(())
    }

    fn finish_document(&mut self) -> Result<(), IoError> {
        self.skip_ws();
        if self.pos == self.bytes.len() {
            Ok(())
        } else {
            Err(malformed(self.pos, "trailing bytes after the JSON array"))
        }
    }

    fn expect(&mut self, expected: u8, what: &str) -> Result<(), IoError> {
        if self.peek() == Some(expected) {
            self.pos += 1;
            Ok(())
        } else {
            Err(malformed(self.pos, what))
        }
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\n' | b'\r' | b'\t')) {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text_column(name: &'static str, required: bool) -> ColumnSpec {
        ColumnSpec {
            name,
            kind: ColumnKind::Text,
            required,
        }
    }

    fn number_column(name: &'static str, min: f64, max: f64) -> ColumnSpec {
        ColumnSpec {
            name,
            kind: ColumnKind::Number { min, max },
            required: true,
        }
    }

    fn parse_rows(input: &str) -> Result<Vec<BTreeMap<String, String>>, IoError> {
        mini_json_array_of_objects(input, CatalogJsonLimits::DEFAULT)
    }

    fn assert_malformed(input: &str, expected_at: usize, expected_detail: &str, case: &str) {
        match parse_rows(input) {
            Err(IoError::Malformed { at, what }) => {
                assert_eq!(at, expected_at, "{case}: wrong refusal offset: {what}");
                assert!(
                    what.contains(expected_detail),
                    "{case}: expected detail {expected_detail:?}, got {what:?} at {at}"
                );
            }
            other => panic!("{case}: expected Malformed, got {other:?}"),
        }
    }

    fn assert_resource(input: &str, limits: CatalogJsonLimits, expected_cap: &str, case: &str) {
        match mini_json_array_of_objects(input, limits) {
            Err(IoError::ResourceBound { what }) => assert!(
                what.contains(expected_cap) && what.contains("byte offset"),
                "{case}: expected cap {expected_cap:?} with offset, got {what:?}"
            ),
            other => panic!("{case}: expected ResourceBound, got {other:?}"),
        }
    }

    /// G0: unchecked declarations cannot become authority across every schema
    /// definition boundary; equal and extreme finite numeric bounds remain
    /// legal inclusive ranges.
    #[test]
    fn g0_schema_admission_refuses_ambiguous_or_invalid_definitions() {
        assert!(matches!(
            Schema::admit(Vec::new()),
            Err(SchemaDefinitionRefusal::EmptySchema)
        ));

        let one_column = CatalogSchemaLimits {
            max_columns: 1,
            max_name_bytes: 8,
            max_total_name_bytes: 8,
        };
        Schema::admit_with_limits(vec![text_column("a", true)], one_column)
            .expect("exact column-count boundary must admit");
        assert!(matches!(
            Schema::admit_with_limits(
                vec![text_column("a", true), text_column("b", true)],
                one_column
            ),
            Err(SchemaDefinitionRefusal::ColumnCount { count: 2, limit: 1 })
        ));

        for (column, expected) in [
            (
                text_column("", true),
                SchemaDefinitionRefusal::EmptyName { column: 1 },
            ),
            (
                text_column(" a", true),
                SchemaDefinitionRefusal::NonCanonicalName { column: 1 },
            ),
            (
                text_column("a ", true),
                SchemaDefinitionRefusal::NonCanonicalName { column: 1 },
            ),
        ] {
            assert_eq!(Schema::admit(vec![column]), Err(expected));
        }

        assert!(matches!(
            Schema::admit_with_limits(
                vec![text_column("abc", true)],
                CatalogSchemaLimits {
                    max_columns: 1,
                    max_name_bytes: 2,
                    max_total_name_bytes: 8,
                }
            ),
            Err(SchemaDefinitionRefusal::NameBytes {
                column: 1,
                bytes: 3,
                limit: 2
            })
        ));
        assert!(matches!(
            Schema::admit_with_limits(
                vec![text_column("ab", true), text_column("cd", true)],
                CatalogSchemaLimits {
                    max_columns: 2,
                    max_name_bytes: 2,
                    max_total_name_bytes: 3,
                }
            ),
            Err(SchemaDefinitionRefusal::TotalNameBytes { bytes: 4, limit: 3 })
        ));
        assert!(matches!(
            Schema::admit(vec![text_column("a", true), text_column("a", false)]),
            Err(SchemaDefinitionRefusal::DuplicateName {
                first_column: 1,
                duplicate_column: 2
            })
        ));

        assert!(matches!(
            Schema::admit(vec![number_column("n", f64::NAN, 1.0)]),
            Err(SchemaDefinitionRefusal::NonFiniteNumberBound {
                column: 1,
                lower: true
            })
        ));
        assert!(matches!(
            Schema::admit(vec![number_column("n", 0.0, f64::INFINITY)]),
            Err(SchemaDefinitionRefusal::NonFiniteNumberBound {
                column: 1,
                lower: false
            })
        ));
        assert!(matches!(
            Schema::admit(vec![number_column("n", 2.0, 1.0)]),
            Err(SchemaDefinitionRefusal::InvertedNumberBounds { column: 1 })
        ));

        Schema::admit(vec![number_column("equal", 1.0, 1.0)])
            .expect("equal inclusive bounds are valid");
        Schema::admit(vec![number_column("finite", f64::MIN, f64::MAX)])
            .expect("extreme finite bounds are valid");
    }

    /// G3/G5: admission is byte-stable and every policy-bearing declaration
    /// input moves the local replay identity.
    #[test]
    fn g3_schema_identity_is_stable_and_policy_sensitive() {
        let columns = vec![text_column("id", true), number_column("value", -1.0, 1.0)];
        let first = Schema::admit(columns.clone()).expect("baseline schema");
        let retry = Schema::admit(columns).expect("identical retry");
        assert_eq!(first.receipt(), retry.receipt());
        assert_eq!(first.receipt().schema_version, CATALOG_SCHEMA_VERSION);

        let variants = [
            Schema::admit(vec![
                text_column("id", false),
                number_column("value", -1.0, 1.0),
            ])
            .expect("required-policy variant"),
            Schema::admit(vec![
                text_column("id", true),
                number_column("value", -2.0, 1.0),
            ])
            .expect("bounds variant"),
            Schema::admit(vec![
                number_column("value", -1.0, 1.0),
                text_column("id", true),
            ])
            .expect("validation-order variant"),
            Schema::admit_with_limits(
                vec![text_column("id", true), number_column("value", -1.0, 1.0)],
                CatalogSchemaLimits {
                    max_columns: 8,
                    ..CatalogSchemaLimits::DEFAULT
                },
            )
            .expect("limit variant"),
        ];
        for variant in variants {
            assert_ne!(
                variant.receipt().local_identity_fnv1a64,
                first.receipt().local_identity_fnv1a64,
                "each authority-bearing schema or limit change must move identity"
            );
        }
    }

    /// G0/G3: CSV header aliases refuse before row-map insertion. Optional
    /// schema columns may be absent and unknown canonical-name columns are
    /// preserved identically by CSV and JSON; document order is immaterial.
    #[test]
    fn g0_g3_csv_header_admission_and_cross_format_projection_policy() {
        let schema = Schema::admit(vec![text_column("id", true), text_column("note", false)])
            .expect("valid projection schema");

        for csv in ["id,id\nA,B\n", "id, id \nA,B\n"] {
            match schema.parse_csv(csv) {
                Err(IoError::Schema {
                    row: 0,
                    column,
                    what,
                }) => {
                    assert_eq!(column, "id");
                    assert!(what.contains("duplicate CSV header"));
                }
                other => panic!("normalized duplicate header must refuse, got {other:?}"),
            }
        }

        let csv = schema
            .parse_csv("extra,id\nkeep,A\n")
            .expect("optional column may be absent and extra column is preserved");
        let permuted_csv = schema
            .parse_csv("id,extra\nA,keep\n")
            .expect("document column permutation must be immaterial");
        let json = schema
            .parse_json(r#"[{"id":"A","extra":"keep"}]"#)
            .expect("JSON follows the same optional/unknown-column policy");
        assert_eq!(csv, permuted_csv);
        assert_eq!(csv, json);
        assert!(!csv.rows[0].contains_key("note"));
        assert_eq!(csv.rows[0]["extra"], "keep");
    }

    /// G0: attacker-sized cell text cannot become an attacker-sized teaching
    /// error even before the shared CSV operation envelope lands.
    #[test]
    fn g0_schema_error_witness_is_bounded() {
        let schema =
            Schema::admit(vec![number_column("n", 0.0, 1.0)]).expect("valid numeric schema");
        let offender = "x".repeat(16 * 1024);
        let csv = format!("n\n{offender}\n");
        match schema.parse_csv(&csv) {
            Err(IoError::Schema { what, .. }) => {
                assert!(
                    what.len() < 192,
                    "diagnostic must stay bounded: {}",
                    what.len()
                );
                assert!(what.contains("UTF-8 bytes"));
            }
            other => panic!("non-number must produce a schema error, got {other:?}"),
        }
    }

    fn assert_catalog_resource(error: IoError, expected: &str, case: &str) {
        match error {
            IoError::ResourceBound { what } => assert!(
                what.contains(expected),
                "{case}: expected {expected:?} in resource refusal, got {what:?}"
            ),
            other => panic!("{case}: expected ResourceBound, got {other:?}"),
        }
    }

    /// G0/G3: CSV and JSON share one exact projection envelope. Every CSV
    /// syntax/payload dimension and every common projection dimension admits
    /// its exact boundary and refuses boundary-minus-one before publication.
    #[test]
    fn g0_g3_common_catalog_operation_budget_is_exact_and_cross_format() {
        let schema = Schema::admit(vec![text_column("id", true), number_column("n", 0.0, 9.0)])
            .expect("valid common catalog schema");
        let projection = CatalogProjectionLimits {
            max_validation_visits: 2,
            max_numeric_entries: 1,
            max_numeric_key_bytes: 1,
            max_output_bytes: 6,
        };
        let csv_text = "id,n\n\"A\",1\n";
        let csv_limits = CatalogCsvLimits {
            max_input_bytes: csv_text.len(),
            max_rows: 1,
            max_fields_per_record: 2,
            max_total_fields: 4,
            max_field_bytes: 2,
            max_header_bytes: 3,
            max_decoded_bytes: 5,
            projection,
        };
        let json_text = r#"[{"id":"A","n":1}]"#;
        let json_limits = CatalogJsonLimits {
            max_input_bytes: json_text.len(),
            max_rows: 1,
            max_members_per_object: 2,
            max_total_members: 2,
            max_string_bytes: 2,
            max_number_bytes: 1,
            max_decoded_bytes: 5,
            projection,
        };

        let csv = schema
            .parse_csv_with_limits(csv_text, csv_limits)
            .expect("CSV exact operation boundary");
        let unquoted_csv = schema
            .parse_csv_with_limits("id,n\nA,1\n", csv_limits)
            .expect("quote-preserving CSV rewrite");
        let json = schema
            .parse_json_with_limits(json_text, json_limits)
            .expect("JSON exact operation boundary");
        assert_eq!(csv, unquoted_csv);
        assert_eq!(csv, json);

        for (limits, expected, case) in [
            (
                CatalogCsvLimits {
                    max_input_bytes: csv_text.len() - 1,
                    ..csv_limits
                },
                "input-byte",
                "CSV input byte boundary",
            ),
            (
                CatalogCsvLimits {
                    max_rows: 0,
                    ..csv_limits
                },
                "row",
                "CSV row boundary",
            ),
            (
                CatalogCsvLimits {
                    max_fields_per_record: 1,
                    ..csv_limits
                },
                "required schema columns per row",
                "CSV per-record field boundary",
            ),
            (
                CatalogCsvLimits {
                    max_total_fields: 3,
                    ..csv_limits
                },
                "aggregate field",
                "CSV aggregate field boundary",
            ),
            (
                CatalogCsvLimits {
                    max_field_bytes: 1,
                    ..csv_limits
                },
                "field decoded-byte",
                "CSV field byte boundary",
            ),
            (
                CatalogCsvLimits {
                    max_header_bytes: 2,
                    ..csv_limits
                },
                "header decoded-byte",
                "CSV header byte boundary",
            ),
            (
                CatalogCsvLimits {
                    max_decoded_bytes: 4,
                    ..csv_limits
                },
                "aggregate decoded-byte",
                "CSV decoded byte boundary",
            ),
            (
                CatalogCsvLimits {
                    projection: CatalogProjectionLimits {
                        max_validation_visits: 1,
                        ..projection
                    },
                    ..csv_limits
                },
                "schema validation visits",
                "CSV validation visit boundary",
            ),
            (
                CatalogCsvLimits {
                    projection: CatalogProjectionLimits {
                        max_numeric_entries: 0,
                        ..projection
                    },
                    ..csv_limits
                },
                "numeric projection entries",
                "CSV numeric entry boundary",
            ),
            (
                CatalogCsvLimits {
                    projection: CatalogProjectionLimits {
                        max_numeric_key_bytes: 0,
                        ..projection
                    },
                    ..csv_limits
                },
                "numeric projection key bytes",
                "CSV numeric key boundary",
            ),
            (
                CatalogCsvLimits {
                    projection: CatalogProjectionLimits {
                        max_output_bytes: 5,
                        ..projection
                    },
                    ..csv_limits
                },
                "logical output bytes",
                "CSV output boundary",
            ),
        ] {
            let error = schema.parse_csv_with_limits(csv_text, limits).unwrap_err();
            assert_catalog_resource(error, expected, case);
        }

        for (limits, expected, case) in [
            (
                CatalogJsonLimits {
                    projection: CatalogProjectionLimits {
                        max_validation_visits: 1,
                        ..projection
                    },
                    ..json_limits
                },
                "schema validation visits",
                "JSON validation visit boundary",
            ),
            (
                CatalogJsonLimits {
                    projection: CatalogProjectionLimits {
                        max_numeric_entries: 0,
                        ..projection
                    },
                    ..json_limits
                },
                "numeric projection entries",
                "JSON numeric entry boundary",
            ),
            (
                CatalogJsonLimits {
                    projection: CatalogProjectionLimits {
                        max_numeric_key_bytes: 0,
                        ..projection
                    },
                    ..json_limits
                },
                "numeric projection key bytes",
                "JSON numeric key boundary",
            ),
            (
                CatalogJsonLimits {
                    projection: CatalogProjectionLimits {
                        max_output_bytes: 5,
                        ..projection
                    },
                    ..json_limits
                },
                "logical output bytes",
                "JSON output boundary",
            ),
        ] {
            let error = schema
                .parse_json_with_limits(json_text, limits)
                .unwrap_err();
            assert_catalog_resource(error, expected, case);
        }

        let optional_schema = Schema::admit(vec![text_column("a", false), text_column("b", false)])
            .expect("valid optional-only schema");
        let empty_object_text = "[{},{}]";
        let empty_object_limits = CatalogJsonLimits {
            max_input_bytes: empty_object_text.len(),
            max_rows: 2,
            max_members_per_object: 0,
            max_total_members: 0,
            max_string_bytes: 0,
            max_number_bytes: 0,
            max_decoded_bytes: 0,
            projection: CatalogProjectionLimits {
                max_validation_visits: 4,
                max_numeric_entries: 0,
                max_numeric_key_bytes: 0,
                max_output_bytes: 0,
            },
        };
        let empty_objects = optional_schema
            .parse_json_with_limits(empty_object_text, empty_object_limits)
            .expect("empty objects still fit the exact validation envelope");
        assert_eq!(empty_objects.rows.len(), 2);
        assert!(empty_objects.rows.iter().all(|row| row.is_empty()));
        let error = optional_schema
            .parse_json_with_limits(
                empty_object_text,
                CatalogJsonLimits {
                    projection: CatalogProjectionLimits {
                        max_validation_visits: 3,
                        ..empty_object_limits.projection
                    },
                    ..empty_object_limits
                },
            )
            .expect_err("member caps cannot hide empty-object validation work");
        assert_catalog_resource(
            error,
            "schema validation visits",
            "JSON empty-object validation boundary",
        );

        let overflow = CatalogCsvLimits {
            max_rows: usize::MAX,
            max_total_fields: usize::MAX,
            projection: CatalogProjectionLimits {
                max_validation_visits: usize::MAX,
                max_numeric_entries: usize::MAX,
                max_numeric_key_bytes: usize::MAX,
                max_output_bytes: usize::MAX,
            },
            ..csv_limits
        };
        let error = schema
            .parse_csv_with_limits(csv_text, overflow)
            .expect_err("checked operation overflow must refuse before row projection");
        assert_catalog_resource(error, "arithmetic overflow", "checked preflight overflow");
    }

    /// G0: every RFC 8259 string escape, BMP escape, surrogate pair, and raw
    /// UTF-8 scalar decodes exactly once into the retained catalog payload.
    #[test]
    fn g0_json_string_escapes_and_surrogates_decode_exactly() {
        let input = r#"[{"simple":"\"\\\/\b\f\n\r\t","nul":"\u0000","bmp":"\u20aC","pair":"\uD834\uDd1E","first":"\uD800\uDC00","last":"\uDBFF\uDFFF","raw":"café–90"}]"#;
        let rows = parse_rows(input).expect("complete RFC string fixture must parse");
        assert_eq!(rows.len(), 1, "one input object must remain one row");
        assert_eq!(
            rows[0]["simple"], "\"\\/\u{0008}\u{000c}\n\r\t",
            "all eight simple escapes must decode with their RFC meaning"
        );
        assert_eq!(rows[0]["nul"], "\0", "escaped NUL is legal JSON text");
        assert_eq!(rows[0]["bmp"], "€", "mixed-case hex digits are legal");
        assert_eq!(
            rows[0]["pair"], "𝄞",
            "UTF-16 surrogate pair must become one scalar"
        );
        assert_eq!(
            rows[0]["first"], "\u{10000}",
            "lowest surrogate pair must map to the first supplementary scalar"
        );
        assert_eq!(
            rows[0]["last"], "\u{10ffff}",
            "highest surrogate pair must map to the last Unicode scalar"
        );
        assert_eq!(rows[0]["raw"], "café–90", "raw UTF-8 must remain exact");
    }

    /// G0: the first malformed escape byte is stable and actionable.
    #[test]
    fn g0_json_malformed_unicode_and_escape_offsets_are_exact() {
        let bad_hex = r#"[{"k":"\u12G4"}]"#;
        assert_malformed(
            bad_hex,
            bad_hex.find('G').expect("fixture has G"),
            "non-hex",
            "bad Unicode hex digit",
        );

        let lone_low = r#"[{"k":"\uDC00"}]"#;
        assert_malformed(
            lone_low,
            lone_low.find("\\uDC00").expect("fixture has low surrogate"),
            "unpaired low surrogate",
            "lone low surrogate",
        );

        let lone_high = r#"[{"k":"\uD800"}]"#;
        assert_malformed(
            lone_high,
            lone_high
                .find("\\uD800")
                .expect("fixture has high surrogate")
                + 6,
            "high surrogate",
            "lone high surrogate",
        );

        let wrong_pair = r#"[{"k":"\uD800\u0041"}]"#;
        assert_malformed(
            wrong_pair,
            wrong_pair
                .find("\\u0041")
                .expect("fixture has second escape"),
            "non-low surrogate",
            "high surrogate followed by BMP scalar",
        );

        let truncated = r#"[{"k":"\u12"#;
        assert_malformed(
            truncated,
            truncated.len(),
            "truncated",
            "truncated Unicode escape",
        );

        let dangling = "[{\"k\":\"\\";
        assert_malformed(
            dangling,
            dangling.len() - 1,
            "dangling",
            "dangling terminal escape",
        );

        let prefix = b"[{\"k\":\"";
        let suffix = b"\"}]";
        for escaped in 0u8..=0x7f {
            if matches!(
                escaped,
                b'"' | b'\\' | b'/' | b'b' | b'f' | b'n' | b'r' | b't' | b'u'
            ) {
                continue;
            }
            let mut bytes = Vec::from(prefix);
            bytes.push(b'\\');
            bytes.push(escaped);
            bytes.extend_from_slice(suffix);
            let input = String::from_utf8(bytes).expect("ASCII escape fixture is UTF-8");
            assert_malformed(
                &input,
                prefix.len() + 1,
                "unknown JSON string escape",
                &format!("unknown escape byte 0x{escaped:02x}"),
            );
        }

        let incomplete_u = r#"[{"k":"\u"}]"#;
        assert_malformed(
            incomplete_u,
            incomplete_u.find("\\u").expect("fixture has escape") + 2,
            "non-hex",
            "Unicode escape without digits",
        );
    }

    /// G0: all 32 raw C0 bytes are forbidden inside a JSON string, including
    /// the four bytes that are legal only as whitespace outside strings.
    #[test]
    fn g0_json_raw_c0_controls_are_exhaustively_rejected() {
        let prefix = b"[{\"k\":\"";
        let suffix = b"\"}]";
        for control in 0u8..=0x1f {
            let mut bytes = Vec::from(prefix);
            bytes.push(control);
            bytes.extend_from_slice(suffix);
            let input = String::from_utf8(bytes).expect("C0 fixture remains valid UTF-8");
            assert_malformed(
                &input,
                prefix.len(),
                "raw C0 control",
                &format!("raw control byte 0x{control:02x}"),
            );
        }
    }

    /// G0: JSON's number production is implemented literally and preserves
    /// each admitted lexeme for downstream schema validation/provenance.
    #[test]
    fn g0_json_number_grammar_accepts_only_rfc_8259_lexemes() {
        for number in [
            "0",
            "-0",
            "10",
            "-10",
            "0.0",
            "-0.125",
            "1e2",
            "1E+2",
            "1e-2",
            "1234567890.0123456789e+123",
        ] {
            let input = format!("[{{\"n\":{number}}}]");
            let rows = parse_rows(&input)
                .unwrap_or_else(|error| panic!("valid number {number:?} refused: {error:?}"));
            assert_eq!(
                rows[0]["n"], number,
                "valid number {number:?} must retain its exact lexeme"
            );
        }

        let prefix = "[{\"n\":";
        for (number, relative_at, detail) in [
            ("+1", 0, "expected a JSON string or number"),
            ("01", 1, "leading zero"),
            ("-01", 2, "leading zero"),
            (".1", 0, "expected a JSON string or number"),
            ("1.", 2, "after JSON number decimal point"),
            ("1e", 2, "JSON number exponent"),
            ("1e+", 3, "JSON number exponent"),
            ("--1", 1, "integer part"),
            ("1_0", 1, "invalid byte after JSON number"),
            ("0x1", 1, "invalid byte after JSON number"),
            ("NaN", 0, "expected a JSON string or number"),
            ("Infinity", 0, "expected a JSON string or number"),
            ("1e--2", 3, "JSON number exponent"),
            ("1..0", 2, "after JSON number decimal point"),
        ] {
            let input = format!("{prefix}{number}}}]");
            assert_malformed(
                &input,
                prefix.len() + relative_at,
                detail,
                &format!("invalid number {number:?}"),
            );
        }
    }

    /// G0: comma/colon grammar is explicit, and duplicate detection operates
    /// on decoded keys so alternate escape spellings cannot overwrite data.
    #[test]
    fn g0_json_delimiters_and_duplicate_keys_are_strict() {
        let missing_object_comma = r#"[{"a":"1" "b":"2"}]"#;
        assert_malformed(
            missing_object_comma,
            missing_object_comma.rfind("\"b\"").expect("second key"),
            "expected ',' or '}'",
            "missing object comma",
        );

        let trailing_object_comma = r#"[{"a":"1",}]"#;
        assert_malformed(
            trailing_object_comma,
            trailing_object_comma.find('}').expect("object close"),
            "trailing comma in JSON object",
            "trailing object comma",
        );

        let missing_array_comma = r#"[{"a":"1"} {"b":"2"}]"#;
        assert_malformed(
            missing_array_comma,
            missing_array_comma.rfind('{').expect("second object"),
            "expected ',' or ']'",
            "missing array comma",
        );

        let trailing_array_comma = r#"[{"a":"1"},]"#;
        assert_malformed(
            trailing_array_comma,
            trailing_array_comma.find(']').expect("array close"),
            "trailing comma in JSON array",
            "trailing array comma",
        );

        let doubled_comma = r#"[{"a":"1",,"b":"2"}]"#;
        assert_malformed(
            doubled_comma,
            doubled_comma.find(",,").expect("double comma") + 1,
            "quoted JSON object key",
            "doubled object comma",
        );

        let missing_colon = r#"[{"a" "1"}]"#;
        assert_malformed(
            missing_colon,
            missing_colon.rfind("\"1\"").expect("value"),
            "expected ':'",
            "missing colon",
        );

        let duplicate = r#"[{"a":"1","a":"2"}]"#;
        assert_malformed(
            duplicate,
            duplicate.rfind("\"a\"").expect("second key"),
            "duplicate JSON object key",
            "literal duplicate key",
        );

        let decoded_duplicate = r#"[{"a":"1","\u0061":"2"}]"#;
        assert_malformed(
            decoded_duplicate,
            decoded_duplicate.find("\"\\u0061\"").expect("escaped key"),
            "duplicate JSON object key",
            "escape-equivalent duplicate key",
        );

        for (input, case) in [
            (r#"[{"a":true}]"#, "boolean value"),
            (r#"[{"a":null}]"#, "null value"),
            (r#"[{"a":[]}]"#, "nested array value"),
            (r#"[{"a":{}}]"#, "nested object value"),
            (r#"[1]"#, "non-object array member"),
        ] {
            assert!(
                matches!(parse_rows(input), Err(IoError::Malformed { .. })),
                "{case}: unsupported flat-catalog value must be Malformed"
            );
        }
    }

    /// G3: insignificant whitespace, key ordering, raw Unicode, and escaped
    /// Unicode are semantic-preserving rewrites of this restricted language.
    #[test]
    fn g3_json_equivalent_rewrites_produce_identical_rows() {
        let compact = r#"[{"a":"café","b":"𝄞","n":-1.25e+2}]"#;
        let whitespace = " \n[ \t{ \"a\" : \"café\" , \"b\" : \"𝄞\" , \"n\" : -1.25e+2 } \r] \t";
        let escaped = r#"[{"n":-1.25e+2,"b":"\uD834\uDD1E","a":"caf\u00e9"}]"#;
        let expected = parse_rows(compact).expect("compact fixture");
        assert_eq!(
            parse_rows(whitespace).expect("whitespace rewrite"),
            expected,
            "RFC whitespace insertion must not move catalog semantics"
        );
        assert_eq!(
            parse_rows(escaped).expect("escape/member-order rewrite"),
            expected,
            "member permutation and equivalent Unicode escaping must agree"
        );
    }

    /// G3: no proper prefix of a valid document may publish a partial row;
    /// each truncation reports a byte offset inside or immediately after the
    /// available prefix.
    #[test]
    fn g3_json_all_truncation_prefixes_refuse_without_partial_results() {
        let complete = r#"[{"a":"\uD834\uDD1E","n":-1.25e+2},{"b":"escaped\ntext"}]"#;
        parse_rows(complete).expect("complete truncation fixture must parse");
        for cut in 0..complete.len() {
            match parse_rows(&complete[..cut]) {
                Err(IoError::Malformed { at, what }) => assert!(
                    at <= cut,
                    "truncation at {cut}: refusal offset {at} is outside prefix; detail={what:?}"
                ),
                other => panic!(
                    "truncation at byte {cut} must not publish a partial catalog; got {other:?}"
                ),
            }
        }
    }

    /// G0/G3: every logical resource dimension accepts its exact boundary and
    /// refuses the first excess before growing the corresponding payload.
    #[test]
    fn g0_json_resource_caps_are_exact_and_compositional() {
        let base = CatalogJsonLimits {
            max_input_bytes: 1024,
            max_rows: 8,
            max_members_per_object: 8,
            max_total_members: 16,
            max_string_bytes: 64,
            max_number_bytes: 64,
            max_decoded_bytes: 256,
            projection: CatalogProjectionLimits::DEFAULT,
        };

        let input = "[]";
        let exact = CatalogJsonLimits {
            max_input_bytes: input.len(),
            ..base
        };
        mini_json_array_of_objects(input, exact).expect("exact input-byte cap");
        assert_resource(
            input,
            CatalogJsonLimits {
                max_input_bytes: input.len() - 1,
                ..exact
            },
            "input-byte",
            "first input byte beyond cap",
        );

        let rows = "[{},{}]";
        mini_json_array_of_objects(
            rows,
            CatalogJsonLimits {
                max_rows: 2,
                ..base
            },
        )
        .expect("exact row cap");
        assert_resource(
            rows,
            CatalogJsonLimits {
                max_rows: 1,
                ..base
            },
            "row cap",
            "first row beyond cap",
        );
        match mini_json_array_of_objects(
            "[{},,]",
            CatalogJsonLimits {
                max_rows: 1,
                ..base
            },
        ) {
            Err(IoError::Malformed { at: 4, .. }) => {}
            other => panic!(
                "malformed token at exhausted row cap must remain a syntax refusal: {other:?}"
            ),
        }

        let members = r#"[{"a":"","b":""}]"#;
        mini_json_array_of_objects(
            members,
            CatalogJsonLimits {
                max_members_per_object: 2,
                ..base
            },
        )
        .expect("exact per-object member cap");
        assert_resource(
            members,
            CatalogJsonLimits {
                max_members_per_object: 1,
                ..base
            },
            "per-object member",
            "first object member beyond cap",
        );
        match mini_json_array_of_objects(
            r#"[{"a":"",,}]"#,
            CatalogJsonLimits {
                max_members_per_object: 1,
                ..base
            },
        ) {
            Err(IoError::Malformed { at, .. })
                if at == r#"[{"a":"",,}]"#.find(",,").expect("double comma") + 1 => {}
            other => panic!(
                "malformed token at exhausted member cap must remain a syntax refusal: {other:?}"
            ),
        }

        let aggregate_members = r#"[{"a":""},{"b":""}]"#;
        mini_json_array_of_objects(
            aggregate_members,
            CatalogJsonLimits {
                max_total_members: 2,
                ..base
            },
        )
        .expect("exact aggregate-member cap");
        assert_resource(
            aggregate_members,
            CatalogJsonLimits {
                max_total_members: 1,
                ..base
            },
            "aggregate member",
            "first aggregate member beyond cap",
        );

        let scalar = r#"[{"k":"\uD834\uDD1E"}]"#;
        for (spelling, case) in [(scalar, "escaped scalar"), (r#"[{"k":"𝄞"}]"#, "raw scalar")] {
            mini_json_array_of_objects(
                spelling,
                CatalogJsonLimits {
                    max_string_bytes: 4,
                    ..base
                },
            )
            .unwrap_or_else(|error| panic!("{case} at exact four-byte cap: {error:?}"));
            assert_resource(
                spelling,
                CatalogJsonLimits {
                    max_string_bytes: 3,
                    ..base
                },
                "string decoded-byte",
                &format!("{case} beyond string cap"),
            );
        }

        let number = r#"[{"n":-1.25e+2}]"#;
        mini_json_array_of_objects(
            number,
            CatalogJsonLimits {
                max_number_bytes: 8,
                ..base
            },
        )
        .expect("eight-byte number at exact token cap");
        assert_resource(
            number,
            CatalogJsonLimits {
                max_number_bytes: 7,
                ..base
            },
            "number-token byte",
            "first number byte beyond token cap",
        );

        let decoded = r#"[{"a":"bc"}]"#;
        mini_json_array_of_objects(
            decoded,
            CatalogJsonLimits {
                max_decoded_bytes: 3,
                ..base
            },
        )
        .expect("key plus value at exact aggregate decoded cap");
        assert_resource(
            decoded,
            CatalogJsonLimits {
                max_decoded_bytes: 2,
                ..base
            },
            "aggregate decoded-byte",
            "first decoded byte beyond aggregate cap",
        );
    }

    /// G3: at document end only the four RFC 8259 whitespace bytes are
    /// semantic-preserving; every other ASCII byte is trailing garbage.
    #[test]
    fn g3_json_trailing_ascii_matrix_only_accepts_json_whitespace() {
        let baseline = parse_rows("[]").expect("empty array baseline");
        for suffix in 0u8..=0x7f {
            let input = format!("[]{}", char::from(suffix));
            if matches!(suffix, b' ' | b'\n' | b'\r' | b'\t') {
                assert_eq!(
                    parse_rows(&input).unwrap_or_else(|error| panic!(
                        "JSON whitespace 0x{suffix:02x}: {error:?}"
                    )),
                    baseline,
                    "JSON whitespace suffix 0x{suffix:02x} must be inert"
                );
            } else {
                assert_malformed(
                    &input,
                    2,
                    "trailing bytes",
                    &format!("trailing ASCII byte 0x{suffix:02x}"),
                );
            }
        }
        let non_ascii_space = "[]\u{00a0}";
        assert_malformed(
            non_ascii_space,
            2,
            "trailing bytes",
            "non-JSON Unicode whitespace suffix",
        );
    }
}
