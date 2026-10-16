//! Cross-domain source index for the Instrument Realism Program (M8-W1).
//!
//! Parsing is pure: it checks metadata and decisions, never fetches sources,
//! reads artifacts, or grants scientific authority. Domain records retain the
//! values and their detailed uncertainty. Repository cross-references are
//! exercised by `tests/instrument_data.rs`, outside browser dependency cones.

#![deny(unsafe_code)]

use std::collections::BTreeSet;
use std::fmt;

/// Repository-relative source ledger, read by consumers rather than embedded.
pub const INSTRUMENT_DATA_LEDGER_LOCATOR: &str = "data/instruments/ledger.tsv";
/// Exact first line of the v1 TSV grammar.
pub const INSTRUMENT_DATA_LEDGER_HEADER: &str = "frankensim.instrument-data-ledger.v1";
/// Exact number of fields in every source row.
pub const LEDGER_COLUMNS: usize = 14;
/// Maximum UTF-8 bytes per field, including notes and the artifact list.
pub const LEDGER_FIELD_BYTE_CAP: usize = 4096;
/// Mandatory second line, in column order.
pub const LEDGER_FIELDS: [&str; LEDGER_COLUMNS] = [
    "source_id",
    "item",
    "kind",
    "what",
    "locator",
    "specimen",
    "licence",
    "redistribution",
    "uncertainty",
    "digitization",
    "decision",
    "artifacts",
    "retrieved_utc",
    "notes",
];

macro_rules! vocabulary {
    ($(#[$meta:meta])* $name:ident { $($(#[$doc:meta])* $variant:ident => $tag:literal),+ $(,)? }) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
        pub enum $name { $($(#[$doc])* $variant),+ }
        impl $name {
            /// Stable spelling in the v1 ledger.
            #[must_use]
            pub const fn tag(self) -> &'static str {
                match self { $(Self::$variant => $tag),+ }
            }
            fn parse(value: &str, field: &'static str) -> Result<Self, LedgerErrorKind> {
                match value {
                    $($tag => Ok(Self::$variant)),+,
                    _ => Err(unknown(field, value)),
                }
            }
        }
    };
}

vocabulary! {
    /// Domain of the indexed source, independent of retention permission.
    DataKind {
        /// Constitutive or physical-property data.
        Material => "material",
        /// Instrument or part geometry.
        Geometry => "geometry",
        /// Recording reserved for grading.
        ReferenceAudio => "reference-audio",
        /// Human performance measurements.
        PlayerData => "player-data",
        /// Controlled artificial excitation measurements.
        MachinePlayerData => "machine-player-data",
        /// Impedance, admittance, modes or decay of an identified object.
        PassiveResponse => "passive-response",
        /// Player anatomy or configuration.
        Anatomy => "anatomy",
        /// Acquisition or reconstruction method only.
        Method => "method",
    }
}
vocabulary! {
    /// Recorded retention terms; not an independent legal determination.
    Redistribution {
        /// Retention permitted under the named source terms.
        Permitted => "permitted",
        /// Facts only; no figure or table reproduction permission implied.
        FactsOnly => "facts-only",
        /// Permission is still required.
        PermissionRequired => "permission-required",
        /// Retention refused.
        Refused => "refused",
    }
}
vocabulary! {
    /// Where uncertainty stands; numerical values stay in domain records.
    UncertaintyStatement {
        /// The source states uncertainty.
        Stated => "stated",
        /// Uncertainty is explicitly unknown.
        Unstated => "unstated",
        /// An identified procedure bounds uncertainty.
        Bounded => "bounded",
    }
}
vocabulary! {
    /// How retained facts were acquired, without implying validation.
    Digitization {
        /// No acquisition yet.
        None => "none",
        /// Values explicitly stated in text.
        TextTranscribed => "text-transcribed",
        /// Values transcribed from tables.
        TableTranscribed => "table-transcribed",
        /// Values extracted from a figure, with an error statement.
        FigureDigitized => "figure-digitized",
        /// Owner measurement.
        OwnerMeasured => "owner-measured",
        /// Measurement by a service.
        ServiceMeasured => "service-measured",
        /// Scan acquisition.
        Scan => "scan",
        /// Photo reconstruction.
        Photogrammetry => "photogrammetry",
        /// Authored fixture rather than measured material data.
        Authored => "authored",
    }
}
vocabulary! {
    /// Source disposition. Admission is metadata, never a fidelity claim.
    Decision {
        /// Retain the declared source artifacts.
        Admit => "admit",
        /// Retain only the identified facts.
        AdmitFactsOnly => "admit-facts-only",
        /// Restricted to fixtures, not acceptance renders.
        FixtureOnly => "fixture-only",
        /// Deliberately refused source.
        Refuse => "refuse",
        /// Permission not established.
        PendingPermission => "pending-permission",
        /// Acquisition or source hunt remains open.
        PendingHunt => "pending-hunt",
    }
}

impl Decision {
    /// Whether the row may index retained artifacts (fixtures included).
    #[must_use]
    pub const fn retains_artifacts(self) -> bool {
        matches!(self, Self::Admit | Self::AdmitFactsOnly | Self::FixtureOnly)
    }
}

/// Typed link to an existing domain record. A reference is not a receipt.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ArtifactRef {
    /// Slug under `data/matdb/seed-v1/`.
    Matdb(String),
    /// Registered value-bearing acoustic case.
    Vvreg(String),
    /// Registered absence or deliberate refusal.
    VvregAbsence(String),
    /// Specimen directory under `data/instruments/geometry/`.
    Geometry(String),
    /// Recording id; resolution is delegated to M8-W6a.
    Reference(String),
    /// Existing repository-relative transcription path.
    Table(String),
}

impl ArtifactRef {
    fn parse(text: &str) -> Result<Self, LedgerErrorKind> {
        let (prefix, value) = text
            .split_once(':')
            .ok_or_else(|| unknown("artifacts", text))?;
        if prefix == "table" {
            // Refs stay within the repository, independently of the host OS.
            if value.is_empty()
                || value.contains('\\')
                || value
                    .split('/')
                    .any(|p| p.is_empty() || p == "." || p == "..")
            {
                return Err(unknown("artifacts", text));
            }
            return Ok(Self::Table(value.to_owned()));
        }
        if !slug(value) {
            return Err(unknown("artifacts", text));
        }
        match prefix {
            "matdb" => Ok(Self::Matdb(value.to_owned())),
            "vvreg" => Ok(Self::Vvreg(value.to_owned())),
            "vvreg-absence" => Ok(Self::VvregAbsence(value.to_owned())),
            "geometry" => Ok(Self::Geometry(value.to_owned())),
            "reference" => Ok(Self::Reference(value.to_owned())),
            _ => Err(unknown("artifacts", text)),
        }
    }
}

/// A validated source row. Text fields preserve their exact input bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LedgerRow {
    /// Unique source key.
    pub source_id: String,
    /// Owning plan item, or `pre-IRP` for backfill.
    pub item: String,
    /// Source domain.
    pub kind: DataKind,
    /// Description of the retained facts or open hunt.
    pub what: String,
    /// Source DOI, URL, object number or owner measurement id.
    pub locator: String,
    /// Identified object; `None` encodes the explicit `-` sentinel.
    pub specimen: Option<String>,
    /// Source licence token, preserved for exact manifest comparison.
    pub licence: String,
    /// Declared retention terms.
    pub redistribution: Redistribution,
    /// Uncertainty status, with details in the linked domain record.
    pub uncertainty: UncertaintyStatement,
    /// Acquisition method.
    pub digitization: Digitization,
    /// Retention or pending decision.
    pub decision: Decision,
    /// Sorted typed references; empty encodes `-`.
    pub artifacts: Vec<ArtifactRef>,
    /// ISO calendar date of the ledger's source review, not scan acquisition.
    pub retrieved_utc: String,
    /// Limitations, defect explanation, or explicit `-`.
    pub notes: String,
}

/// Parsed rows sorted by source id, making all lookups row-order independent.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InstrumentDataLedger {
    rows: Vec<LedgerRow>,
}

impl InstrumentDataLedger {
    /// Sources in deterministic source-id order.
    #[must_use]
    pub fn rows(&self) -> &[LedgerRow] {
        &self.rows
    }
    /// Exact source lookup.
    #[must_use]
    pub fn by_source(&self, id: &str) -> Option<&LedgerRow> {
        self.rows
            .binary_search_by(|r| r.source_id.as_str().cmp(id))
            .ok()
            .map(|i| &self.rows[i])
    }
    /// All rows naming this artifact, in source-id order.
    #[must_use]
    pub fn for_artifact(&self, artifact: &ArtifactRef) -> Vec<&LedgerRow> {
        self.rows
            .iter()
            .filter(|r| r.artifacts.contains(artifact))
            .collect()
    }
    /// All sources assigned to one plan item, in source-id order.
    #[must_use]
    pub fn for_item(&self, item: &str) -> Vec<&LedgerRow> {
        self.rows.iter().filter(|r| r.item == item).collect()
    }
}

/// Structured schema or decision refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LedgerErrorKind {
    /// Schema identifier or column heading does not match.
    Header,
    /// Wrong number of row fields.
    ColumnCount { expected: usize, found: usize },
    /// A mandatory column was empty or whitespace-only.
    MissingField(&'static str),
    /// Unknown vocabulary, malformed id, path, licence token or date.
    UnknownValue { field: &'static str, value: String },
    /// Repeated source key; no first-row-wins fallback.
    DuplicateSourceId(String),
    /// Field exceeded the UTF-8 byte cap.
    FieldTooLong { field: &'static str, cap: usize },
    /// Row-local retention or identity rule failed.
    Inconsistent(&'static str),
}

/// Refusal with a one-based physical input line number.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LedgerError {
    /// One-based input line (including the two header lines).
    pub line: usize,
    /// Structured cause.
    pub kind: LedgerErrorKind,
}

impl fmt::Display for LedgerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "instrument ledger line {}: {:?}", self.line, self.kind)
    }
}
impl std::error::Error for LedgerError {}

fn unknown(field: &'static str, value: &str) -> LedgerErrorKind {
    LedgerErrorKind::UnknownValue {
        field,
        value: value.to_owned(),
    }
}
fn slug(value: &str) -> bool {
    !value.is_empty()
        && value != "-"
        && value
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
}
fn iso_date(value: &str) -> bool {
    let b = value.as_bytes();
    if b.len() != 10
        || b[4] != b'-'
        || b[7] != b'-'
        || !b
            .iter()
            .enumerate()
            .all(|(i, b)| i == 4 || i == 7 || b.is_ascii_digit())
    {
        return false;
    }
    let year: u32 = value[..4].parse().unwrap_or(0);
    let month: u32 = value[5..7].parse().unwrap_or(0);
    let day: u32 = value[8..].parse().unwrap_or(0);
    let days = match month {
        2 if year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400)) => {
            29
        }
        2 => 28,
        4 | 6 | 9 | 11 => 30,
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        _ => 0,
    };
    year > 0 && day > 0 && day <= days
}

fn parse_row(fields: &[&str]) -> Result<LedgerRow, LedgerErrorKind> {
    if fields.len() != LEDGER_COLUMNS {
        return Err(LedgerErrorKind::ColumnCount {
            expected: LEDGER_COLUMNS,
            found: fields.len(),
        });
    }
    for (&field, value) in LEDGER_FIELDS.iter().zip(fields) {
        if value.len() > LEDGER_FIELD_BYTE_CAP {
            return Err(LedgerErrorKind::FieldTooLong {
                field,
                cap: LEDGER_FIELD_BYTE_CAP,
            });
        }
        if value.trim().is_empty() {
            return Err(LedgerErrorKind::MissingField(field));
        }
        if value.chars().any(char::is_control) {
            return Err(unknown(field, value));
        }
    }
    if !slug(fields[0]) {
        return Err(unknown("source_id", fields[0]));
    }
    if fields[5] != "-" && !slug(fields[5]) {
        return Err(unknown("specimen", fields[5]));
    }
    if !iso_date(fields[12]) {
        return Err(unknown("retrieved_utc", fields[12]));
    }
    // SPDX-shaped identifiers and legacy matdb tokens are preserved exactly.
    // Catalogue membership and legal authority are not inferred by this parser.
    let licence = fields[6];
    let token = licence.strip_prefix("terms:").unwrap_or(licence);
    if token.is_empty()
        || !token
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'.' | b'+'))
    {
        return Err(unknown("licence", licence));
    }
    let mut artifacts = if fields[11] == "-" {
        Vec::new()
    } else {
        fields[11]
            .split(';')
            .map(ArtifactRef::parse)
            .collect::<Result<Vec<_>, _>>()?
    };
    artifacts.sort();
    artifacts.dedup();
    let row = LedgerRow {
        source_id: fields[0].to_owned(),
        item: fields[1].to_owned(),
        kind: DataKind::parse(fields[2], "kind")?,
        what: fields[3].to_owned(),
        locator: fields[4].to_owned(),
        specimen: (fields[5] != "-").then(|| fields[5].to_owned()),
        licence: licence.to_owned(),
        redistribution: Redistribution::parse(fields[7], "redistribution")?,
        uncertainty: UncertaintyStatement::parse(fields[8], "uncertainty")?,
        digitization: Digitization::parse(fields[9], "digitization")?,
        decision: Decision::parse(fields[10], "decision")?,
        artifacts,
        retrieved_utc: fields[12].to_owned(),
        notes: fields[13].to_owned(),
    };
    let pending = matches!(
        row.decision,
        Decision::PendingPermission | Decision::PendingHunt
    );
    if row.licence == "unknown" && !pending {
        return Err(LedgerErrorKind::Inconsistent(
            "unknown licence requires a pending decision",
        ));
    }
    if row.decision.retains_artifacts()
        && !matches!(
            row.redistribution,
            Redistribution::Permitted | Redistribution::FactsOnly
        )
    {
        return Err(LedgerErrorKind::Inconsistent(
            "retained artifacts require permitted or facts-only redistribution",
        ));
    }
    if row.digitization == Digitization::FigureDigitized
        && (row.redistribution != Redistribution::Permitted
            || row.uncertainty == UncertaintyStatement::Unstated)
    {
        return Err(LedgerErrorKind::Inconsistent(
            "figure digitization requires permission and stated or bounded error",
        ));
    }
    if row.decision == Decision::FixtureOnly
        && row.digitization != Digitization::Authored
        && row.notes == "-"
    {
        return Err(LedgerErrorKind::Inconsistent(
            "non-authored fixture requires a note naming its defect",
        ));
    }
    if !row.decision.retains_artifacts()
        && row
            .artifacts
            .iter()
            .any(|a| !matches!(a, ArtifactRef::VvregAbsence(_)))
    {
        return Err(LedgerErrorKind::Inconsistent(
            "refused or pending sources may reference only absences",
        ));
    }
    if matches!(
        row.kind,
        DataKind::Geometry | DataKind::Anatomy | DataKind::PassiveResponse
    ) && row.specimen.is_none()
    {
        return Err(LedgerErrorKind::Inconsistent(
            "geometry, anatomy and passive response require a specimen",
        ));
    }
    Ok(row)
}

/// Parse the exact two-header-line TSV format, refusing blanks and unknowns.
///
/// Rows and references are sorted for deterministic lookups. There is no file
/// I/O, implicit uncertainty, source retrieval, scientific validation, or claim
/// promotion. A header-only ledger is valid; seed non-vacuity is tested separately.
pub fn parse_instrument_data_ledger(text: &str) -> Result<InstrumentDataLedger, LedgerError> {
    let mut lines = text.lines();
    if lines.next() != Some(INSTRUMENT_DATA_LEDGER_HEADER) {
        return Err(LedgerError {
            line: 1,
            kind: LedgerErrorKind::Header,
        });
    }
    if lines.next() != Some(LEDGER_FIELDS.join("\t").as_str()) {
        return Err(LedgerError {
            line: 2,
            kind: LedgerErrorKind::Header,
        });
    }
    let mut rows = Vec::new();
    let mut ids = BTreeSet::new();
    for (index, line) in lines.enumerate() {
        let line_number = index + 3;
        let row = parse_row(&line.split('\t').collect::<Vec<_>>()).map_err(|kind| LedgerError {
            line: line_number,
            kind,
        })?;
        if !ids.insert(row.source_id.clone()) {
            return Err(LedgerError {
                line: line_number,
                kind: LedgerErrorKind::DuplicateSourceId(row.source_id),
            });
        }
        rows.push(row);
    }
    rows.sort_by(|a, b| a.source_id.cmp(&b.source_id));
    Ok(InstrumentDataLedger { rows })
}
