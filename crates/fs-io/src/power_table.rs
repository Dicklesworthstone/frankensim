//! Tabular component-power ingestion through the quarantine path.
//!
//! Real component power arrives as a spreadsheet export from a power-budget
//! tool. `fs-conduction` deliberately admits only already-parsed typed rows —
//! file parsing is a quarantine concern, not a physics one — so this module
//! turns a declared table into rows that `ComponentPower::new` can consume.
//!
//! Doctrine:
//! - **The unit is declared, checked, and recorded.** A die table in
//!   milliwatts is ordinary, and assuming watts because a header says "power"
//!   is a silent factor-of-1000 error. The file must declare its unit from a
//!   closed set; an undeclared or unrecognized unit REFUSES. The declared unit
//!   and the exact factor applied are both recorded in the receipt, so the
//!   conversion is auditable rather than invisible.
//! - **Refusals name the row AND the physical line.** A spreadsheet export is
//!   edited by humans in a text editor, so a record ordinal alone is not
//!   actionable. Directive lines, comments and blanks all advance the line
//!   counter.
//! - **Nothing is coerced or skipped.** A malformed row, a missing column, a
//!   non-numeric wattage, or a duplicate component name refuses loudly.
//! - **The parser does not resolve names to geometry.** It produces names and
//!   watts; binding a name to a vertex set belongs to the scenario layer that
//!   owns persistent entity identity. Keeping that split is what stops the
//!   parser from needing a mesh.
//! - **The parser does not audit the total.** `PowerMap::new` and its
//!   volumetric lowering already refuse a table that does not balance, and
//!   re-implementing that here would let the two disagree. A declared total is
//!   carried through verbatim so the existing refusal surfaces downstream.
//!
//! No-claim boundary: ingesting a file does not validate its contents. A
//! parsed power table is a caller declaration with quarantine provenance,
//! exactly like the mesh import path — it says the bytes were well-formed and
//! the declared unit was applied, never that the watts are real.

use std::collections::BTreeMap;

use crate::IoError;
use crate::catalog::{CatalogCsvLimits, ColumnKind, ColumnSpec, Schema};

/// Semantics version of the power-table ingestion surface.
pub const POWER_TABLE_SEMANTICS_VERSION: &str = "fs-io-power-table-v1";

/// Directive that declares the table's power unit.
pub const POWER_UNIT_DIRECTIVE: &str = "#power-unit:";

/// Directive that declares the table's total, in the declared unit.
pub const DECLARED_TOTAL_DIRECTIVE: &str = "#declared-total:";

/// Required component-name column.
pub const COMPONENT_COLUMN: &str = "component";

/// Required power column, in the declared unit.
pub const POWER_COLUMN: &str = "power";

/// Optional symmetric half-width column, in the declared unit.
pub const HALF_WIDTH_COLUMN: &str = "half_width";

/// Optional confidence column for the half width.
pub const CONFIDENCE_COLUMN: &str = "confidence";

/// The power unit a table declares. A closed set: an unrecognized unit is a
/// refusal, never a guess.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PowerUnit {
    /// Watts.
    Watts,
    /// Milliwatts.
    Milliwatts,
    /// Kilowatts.
    Kilowatts,
}

impl PowerUnit {
    /// Parse a declared unit token.
    #[must_use]
    pub fn parse(token: &str) -> Option<Self> {
        match token.trim() {
            "W" | "w" | "watt" | "watts" => Some(Self::Watts),
            "mW" | "mw" | "milliwatt" | "milliwatts" => Some(Self::Milliwatts),
            "kW" | "kw" | "kilowatt" | "kilowatts" => Some(Self::Kilowatts),
            _ => None,
        }
    }

    /// Stable label, as recorded in the receipt.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Watts => "W",
            Self::Milliwatts => "mW",
            Self::Kilowatts => "kW",
        }
    }

    /// Multiplier converting a declared value into watts.
    #[must_use]
    pub const fn watts_factor(self) -> f64 {
        match self {
            Self::Watts => 1.0,
            Self::Milliwatts => 1e-3,
            Self::Kilowatts => 1e3,
        }
    }

    /// Every accepted spelling, for diagnostics.
    #[must_use]
    pub const fn accepted() -> &'static str {
        "W, mW, kW"
    }
}

/// A refusal while ingesting a power table.
///
/// Carries BOTH the physical source line and, where one exists, the 1-based
/// data-row ordinal. The catalog reader reports the ordinal; a human editing a
/// spreadsheet export needs the line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PowerTableError {
    /// One-based physical source line, or zero for a whole-document refusal.
    pub line: usize,
    /// One-based data-row ordinal, when the refusal is attributable to a row.
    pub row: Option<usize>,
    /// Stable field or stage token.
    pub field: &'static str,
    /// Actionable diagnosis.
    pub reason: String,
}

impl core::fmt::Display for PowerTableError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match (self.line, self.row) {
            (0, _) => write!(
                formatter,
                "power table refused at {}: {}",
                self.field, self.reason
            ),
            (line, Some(row)) => write!(
                formatter,
                "power table line {line} (data row {row}) field {}: {}",
                self.field, self.reason
            ),
            (line, None) => write!(
                formatter,
                "power table line {line} field {}: {}",
                self.field, self.reason
            ),
        }
    }
}

impl std::error::Error for PowerTableError {}

fn refuse(line: usize, row: Option<usize>, field: &'static str, reason: String) -> PowerTableError {
    PowerTableError {
        line,
        row,
        field,
        reason,
    }
}

/// One admitted component row, in WATTS.
///
/// Deliberately carries no geometry. Binding a name to a vertex set is the
/// scenario layer's responsibility; a parser that resolved names would need a
/// mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerTableRow {
    name: String,
    watts: f64,
    half_width_w: Option<f64>,
    confidence: Option<f64>,
    line: usize,
}

impl PowerTableRow {
    /// Declared component name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Dissipation in watts, after the declared-unit conversion.
    #[must_use]
    pub const fn watts(&self) -> f64 {
        self.watts
    }

    /// Symmetric half width in watts, when declared.
    #[must_use]
    pub const fn half_width_w(&self) -> Option<f64> {
        self.half_width_w
    }

    /// Confidence for the half width, when declared.
    #[must_use]
    pub const fn confidence(&self) -> Option<f64> {
        self.confidence
    }

    /// Physical source line this row came from.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

/// Provenance for one ingested power table.
///
/// Mirrors the import-receipt shape: what was read, from which bytes, by which
/// parser, under which declared unit, and with what trust.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerTableReceipt {
    semantics: &'static str,
    source_hash: u64,
    parser_version: &'static str,
    declared_unit: PowerUnit,
    watts_factor: f64,
    rows: usize,
    declared_total_w: Option<f64>,
}

impl PowerTableReceipt {
    /// Semantics version of the ingestion surface.
    #[must_use]
    pub const fn semantics(&self) -> &'static str {
        self.semantics
    }

    /// Content hash of the exact source bytes.
    #[must_use]
    pub const fn source_hash(&self) -> u64 {
        self.source_hash
    }

    /// Parser version that produced this receipt.
    #[must_use]
    pub const fn parser_version(&self) -> &'static str {
        self.parser_version
    }

    /// The unit the file declared.
    #[must_use]
    pub const fn declared_unit(&self) -> PowerUnit {
        self.declared_unit
    }

    /// The exact factor applied to reach watts.
    #[must_use]
    pub const fn watts_factor(&self) -> f64 {
        self.watts_factor
    }

    /// Admitted row count.
    #[must_use]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    /// Declared total in watts, when the file declared one.
    #[must_use]
    pub const fn declared_total_w(&self) -> Option<f64> {
        self.declared_total_w
    }

    /// Canonical JSON, in the import-receipt house style.
    #[must_use]
    pub fn to_json(&self) -> String {
        let total = self.declared_total_w.map_or_else(
            || "null".to_string(),
            |value| {
                format!(
                    "{{\"watts\":{value},\"bits\":\"{:016x}\"}}",
                    value.to_bits()
                )
            },
        );
        format!(
            "{{\"kind\":\"power-table-receipt\",\"semantics\":\"{}\",\"source_hash\":\"{:016x}\",\
             \"parser\":\"{}\",\"declared_unit\":\"{}\",\"watts_factor\":{},\"rows\":{},\
             \"declared_total_w\":{total},\"trust\":\"quarantined\"}}",
            self.semantics,
            self.source_hash,
            self.parser_version,
            self.declared_unit.label(),
            self.watts_factor,
            self.rows,
        )
    }
}

/// An ingested power table: rows in watts, plus its quarantine provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerTable {
    rows: Vec<PowerTableRow>,
    receipt: PowerTableReceipt,
}

impl PowerTable {
    /// Admitted rows, in declaration order.
    #[must_use]
    pub fn rows(&self) -> &[PowerTableRow] {
        &self.rows
    }

    /// Quarantine provenance.
    #[must_use]
    pub const fn receipt(&self) -> &PowerTableReceipt {
        &self.receipt
    }

    /// Declared total in watts, when the file declared one.
    ///
    /// Passed through verbatim so the existing `PowerMap` balance refusal is
    /// the one that fires. This module never corrects a total.
    #[must_use]
    pub const fn declared_total_w(&self) -> Option<f64> {
        self.receipt.declared_total_w
    }

    /// Sum of admitted row wattages. A diagnostic, NOT an audit: the
    /// authoritative balance check lives in `PowerMap`.
    #[must_use]
    pub fn row_total_w(&self) -> f64 {
        self.rows.iter().map(PowerTableRow::watts).sum()
    }
}

fn parse_directive_value<'a>(line: &'a str, directive: &str) -> Option<&'a str> {
    let trimmed = line.trim_start();
    let compact: String = trimmed.chars().filter(|c| !c.is_whitespace()).collect();
    let key: String = directive.chars().filter(|c| !c.is_whitespace()).collect();
    if compact.starts_with(&key) {
        let idx = trimmed.find(':')?;
        Some(trimmed[idx + 1..].trim())
    } else {
        None
    }
}

fn schema() -> Result<Schema, PowerTableError> {
    Schema::admit(vec![
        ColumnSpec {
            name: COMPONENT_COLUMN,
            kind: ColumnKind::Text,
            required: true,
        },
        ColumnSpec {
            name: POWER_COLUMN,
            kind: ColumnKind::Number {
                min: 0.0,
                max: f64::MAX,
            },
            required: true,
        },
        ColumnSpec {
            name: HALF_WIDTH_COLUMN,
            kind: ColumnKind::Number {
                min: 0.0,
                max: f64::MAX,
            },
            required: false,
        },
        ColumnSpec {
            name: CONFIDENCE_COLUMN,
            kind: ColumnKind::Number { min: 0.0, max: 1.0 },
            required: false,
        },
    ])
    .map_err(|refusal| refuse(0, None, "schema", format!("{refusal:?}")))
}

/// Ingest a declared component-power table.
///
/// The document is a directive prologue followed by a CSV body:
///
/// ```text
/// #power-unit: mW
/// #declared-total: 4200
/// component,power,half_width,confidence
/// die,1200,50,0.95
/// ```
///
/// The unit directive is mandatory. Comment (`#`) and blank lines are skipped
/// but still advance the physical line counter, so a refusal names the line a
/// human would see in an editor.
///
/// # Errors
///
/// [`PowerTableError`] for a missing or unrecognized unit directive, a
/// malformed declared total, an empty table, a duplicate component name, or
/// any refusal the catalog schema raises — remapped onto the physical line.
#[allow(clippy::too_many_lines)] // One linear ingest; splitting it would separate a refusal from its line.
pub fn parse_power_table(input: &str) -> Result<PowerTable, PowerTableError> {
    let source_hash = fs_obs::fnv1a64(input.as_bytes());

    let mut declared_unit: Option<(PowerUnit, usize)> = None;
    let mut declared_total_raw: Option<(f64, usize)> = None;
    let mut body = String::new();
    let mut body_first_line: Option<usize> = None;

    for (index, line) in input.lines().enumerate() {
        let line_number = index + 1;
        let trimmed = line.trim();

        if body_first_line.is_none() {
            if let Some(value) = parse_directive_value(line, POWER_UNIT_DIRECTIVE) {
                if declared_unit.is_some() {
                    return Err(refuse(
                        line_number,
                        None,
                        "power-unit",
                        "the power unit is declared more than once".to_string(),
                    ));
                }
                let unit = PowerUnit::parse(value).ok_or_else(|| {
                    refuse(
                        line_number,
                        None,
                        "power-unit",
                        format!(
                            "declared unit {value:?} is not recognized; accepted units are {}",
                            PowerUnit::accepted()
                        ),
                    )
                })?;
                declared_unit = Some((unit, line_number));
                continue;
            }
            if let Some(value) = parse_directive_value(line, DECLARED_TOTAL_DIRECTIVE) {
                let parsed: f64 = value.trim().parse().map_err(|_| {
                    refuse(
                        line_number,
                        None,
                        "declared-total",
                        format!("declared total {value:?} is not a number"),
                    )
                })?;
                if !parsed.is_finite() || parsed < 0.0 {
                    return Err(refuse(
                        line_number,
                        None,
                        "declared-total",
                        format!("declared total {parsed} is not finite and non-negative"),
                    ));
                }
                declared_total_raw = Some((parsed, line_number));
                continue;
            }
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            body_first_line = Some(line_number);
        }

        body.push_str(line);
        body.push('\n');
    }

    let (unit, _unit_line) = declared_unit.ok_or_else(|| {
        refuse(
            0,
            None,
            "power-unit",
            format!(
                "no {POWER_UNIT_DIRECTIVE} directive; a power table must declare its unit \
                 ({}) rather than let a reader assume watts",
                PowerUnit::accepted()
            ),
        )
    })?;

    let header_line = body_first_line.ok_or_else(|| {
        refuse(
            0,
            None,
            "body",
            "the table declares a unit but has no header row or data".to_string(),
        )
    })?;

    let schema = schema()?;
    let catalog = schema
        .parse_csv_with_limits(&body, CatalogCsvLimits::DEFAULT)
        .map_err(|error| map_catalog_error(error, header_line))?;

    if catalog.rows.is_empty() {
        return Err(refuse(
            header_line,
            None,
            "body",
            "the table has a header but no component rows".to_string(),
        ));
    }

    let factor = unit.watts_factor();
    let mut rows: Vec<PowerTableRow> = Vec::new();
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();

    for (index, text_row) in catalog.rows.iter().enumerate() {
        let row_ordinal = index + 1;
        let line = header_line + row_ordinal;
        let numbers = &catalog.numbers[index];

        let name = text_row
            .get(COMPONENT_COLUMN)
            .map(|value| value.trim().to_string())
            .unwrap_or_default();
        if name.is_empty() {
            return Err(refuse(
                line,
                Some(row_ordinal),
                COMPONENT_COLUMN,
                "component name is empty".to_string(),
            ));
        }
        if let Some(previous) = seen.get(&name) {
            return Err(refuse(
                line,
                Some(row_ordinal),
                COMPONENT_COLUMN,
                format!("component {name:?} is already declared on line {previous}"),
            ));
        }
        seen.insert(name.clone(), line);

        let declared = *numbers.get(POWER_COLUMN).ok_or_else(|| {
            refuse(
                line,
                Some(row_ordinal),
                POWER_COLUMN,
                "power cell is missing".to_string(),
            )
        })?;
        let watts = declared * factor;
        if !watts.is_finite() {
            return Err(refuse(
                line,
                Some(row_ordinal),
                POWER_COLUMN,
                format!(
                    "power {declared} {} does not convert to a finite wattage",
                    unit.label()
                ),
            ));
        }

        let half_width_w = numbers.get(HALF_WIDTH_COLUMN).map(|value| value * factor);
        let confidence = numbers.get(CONFIDENCE_COLUMN).copied();
        if half_width_w.is_some() && confidence.is_none() {
            return Err(refuse(
                line,
                Some(row_ordinal),
                CONFIDENCE_COLUMN,
                format!(
                    "a declared {HALF_WIDTH_COLUMN} needs its {CONFIDENCE_COLUMN}; an interval \
                     without a confidence states no uncertainty"
                ),
            ));
        }
        if let Some(confidence) = confidence
            && !(confidence > 0.0 && confidence < 1.0)
        {
            return Err(refuse(
                line,
                Some(row_ordinal),
                CONFIDENCE_COLUMN,
                format!("confidence {confidence} is not strictly inside (0, 1)"),
            ));
        }

        rows.push(PowerTableRow {
            name,
            watts,
            half_width_w,
            confidence,
            line,
        });
    }

    let declared_total_w = declared_total_raw.map(|(value, _)| value * factor);
    let receipt = PowerTableReceipt {
        semantics: POWER_TABLE_SEMANTICS_VERSION,
        source_hash,
        parser_version: crate::VERSION,
        declared_unit: unit,
        watts_factor: factor,
        rows: rows.len(),
        declared_total_w,
    };

    Ok(PowerTable { rows, receipt })
}

/// Remap a catalog refusal onto the physical source line.
///
/// The catalog reader reports a data-row ordinal with the header excluded.
/// Adding the header's own physical line recovers the line a human sees.
fn map_catalog_error(error: IoError, header_line: usize) -> PowerTableError {
    match error {
        IoError::Schema { row, column, what } => {
            if row == 0 {
                // Row zero is a header-level refusal in the catalog reader.
                PowerTableError {
                    line: header_line,
                    row: None,
                    field: "header",
                    reason: format!("column {column:?}: {what}"),
                }
            } else {
                PowerTableError {
                    line: header_line + row,
                    row: Some(row),
                    field: "cell",
                    reason: format!("column {column:?}: {what}"),
                }
            }
        }
        IoError::Malformed { at, what } => PowerTableError {
            line: header_line + at,
            row: None,
            field: "record",
            reason: what,
        },
        IoError::ResourceBound { what } => PowerTableError {
            line: 0,
            row: None,
            field: "resource-bound",
            reason: what,
        },
        IoError::Unsupported { what } => PowerTableError {
            line: 0,
            row: None,
            field: "unsupported",
            reason: what,
        },
        IoError::Cancelled { stage, at } => PowerTableError {
            line: 0,
            row: None,
            field: "cancelled",
            reason: format!("catalog read cancelled at {stage} position {at}"),
        },
    }
}
