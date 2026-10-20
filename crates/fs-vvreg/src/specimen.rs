//! Specimen geometry tables for the Instrument Realism Program (M8-W3a..n
//! under the plan §8 geometry completion rule and the M6.3 source
//! hierarchy).
//!
//! A specimen directory `data/instruments/geometry/<id>/` holds small
//! tab-separated tables. Each starts with an exact schema line and an exact
//! column line ([`SpecimenTable`]). Parsing is pure: callers read the bytes,
//! nothing is fetched, and no value is guessed or defaulted. An absent
//! required value refuses with a typed error that names it
//! ([`SpecimenErrorKind::Missing`]); an unknown uncertainty is the explicit
//! word `unstated`, never a blank. Lookups are keyed and sorted, so row order
//! never changes a result (G3). Every parsed table keeps the BLAKE3 digest of
//! the exact bytes it read, so a consumer can record which version it used.
//!
//! No-claim: admission checks form, units and cross-table identity. It does
//! not check that a value is physically right, that a source is authentic, or
//! that a specimen is complete; `unknowns.tsv` carries that boundary.

#![deny(unsafe_code)]

use fs_blake3::{ContentHash, hash_bytes};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Repository-relative root of the specimen directories.
pub const SPECIMEN_GEOMETRY_ROOT: &str = "data/instruments/geometry";
/// Maximum UTF-8 bytes per field.
pub const SPECIMEN_FIELD_BYTE_CAP: usize = 4096;
/// Keys every `specimen.tsv` must carry (further keys are kept verbatim).
pub const SPECIMEN_KEYS: [&str; 11] = [
    "specimen_id",
    "maker",
    "model",
    "serial",
    "year",
    "system",
    "holder",
    "reference_state",
    "environment",
    "source_ids",
    "target_relation",
];

/// The table kinds of one specimen directory.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SpecimenTable {
    /// `specimen.tsv`: identity, rights and reference state (key/value).
    Specimen,
    /// `frames.tsv`: one SI frame per described part.
    Frames,
    /// `dimensions.tsv`: long-format SI values with uncertainty.
    Dimensions,
    /// `parts.tsv`: material, state and joins per part.
    Parts,
    /// `unknowns.tsv`: named gaps and the claims they block.
    Unknowns,
    /// `fingerings.tsv`: hole states per named fingering (winds only).
    Fingerings,
}

impl SpecimenTable {
    /// File name inside the specimen directory.
    #[must_use]
    pub const fn file_name(self) -> &'static str {
        match self {
            Self::Specimen => "specimen.tsv",
            Self::Frames => "frames.tsv",
            Self::Dimensions => "dimensions.tsv",
            Self::Parts => "parts.tsv",
            Self::Unknowns => "unknowns.tsv",
            Self::Fingerings => "fingerings.tsv",
        }
    }

    /// Exact first line.
    #[must_use]
    pub const fn schema(self) -> &'static str {
        match self {
            Self::Specimen => "frankensim.specimen.v1",
            Self::Frames => "frankensim.specimen-frames.v1",
            Self::Dimensions => "frankensim.specimen-dimensions.v1",
            Self::Parts => "frankensim.specimen-parts.v1",
            Self::Unknowns => "frankensim.specimen-unknowns.v1",
            Self::Fingerings => "frankensim.specimen-fingerings.v1",
        }
    }

    /// Exact second line, in column order.
    #[must_use]
    pub const fn fields(self) -> &'static [&'static str] {
        match self {
            Self::Specimen => &["key", "value"],
            Self::Frames => &["part", "origin", "axis_x", "axis_y", "axis_z", "note"],
            Self::Dimensions => &[
                "part",
                "feature",
                "quantity",
                "value_si",
                "unit",
                "sigma_si",
                "route",
                "source_id",
                "note",
            ],
            Self::Parts => &["part", "material", "state", "joins", "note"],
            Self::Unknowns => &[
                "id",
                "what",
                "bound",
                "correlation_group",
                "blocks",
                "route",
            ],
            Self::Fingerings => &[
                "fingering",
                "written",
                "concert",
                "register",
                "open_holes",
                "source_id",
                "note",
            ],
        }
    }
}

/// Structured refusal cause.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpecimenErrorKind {
    /// Schema line or column line does not match.
    Header,
    /// Wrong number of row fields.
    ColumnCount {
        /// Columns the schema declares.
        expected: usize,
        /// Columns the row carries.
        found: usize,
    },
    /// A field was empty or whitespace-only.
    MissingField(&'static str),
    /// A field exceeded [`SPECIMEN_FIELD_BYTE_CAP`].
    FieldTooLong {
        /// Column name.
        field: &'static str,
        /// The byte cap.
        cap: usize,
    },
    /// Unknown vocabulary, malformed identifier or number.
    UnknownValue {
        /// Column name.
        field: &'static str,
        /// The rejected text.
        value: String,
    },
    /// Repeated key; no first-row-wins fallback.
    Duplicate(String),
    /// A required key, row or value is absent; names it.
    Missing(String),
    /// Two tables disagree (a part without a frame, an unlisted source).
    Inconsistent(String),
}

/// Refusal naming the table and the one-based input line (0: not line-bound).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpecimenError {
    /// Table that refused.
    pub table: SpecimenTable,
    /// One-based line including the two header lines; 0 for lookups.
    pub line: usize,
    /// Structured cause.
    pub kind: SpecimenErrorKind,
}

impl fmt::Display for SpecimenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} line {}: {:?}",
            self.table.file_name(),
            self.line,
            self.kind
        )
    }
}
impl std::error::Error for SpecimenError {}

fn err(table: SpecimenTable, line: usize, kind: SpecimenErrorKind) -> SpecimenError {
    SpecimenError { table, line, kind }
}
fn unknown(field: &'static str, value: &str) -> SpecimenErrorKind {
    SpecimenErrorKind::UnknownValue {
        field,
        value: value.to_owned(),
    }
}
fn slug(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
}

/// Header-checked rows as `(line, fields)`.
fn table_rows(table: SpecimenTable, text: &str) -> Result<Vec<(usize, Vec<&str>)>, SpecimenError> {
    let mut lines = text.lines();
    if lines.next() != Some(table.schema()) {
        return Err(err(table, 1, SpecimenErrorKind::Header));
    }
    if lines.next() != Some(table.fields().join("\t").as_str()) {
        return Err(err(table, 2, SpecimenErrorKind::Header));
    }
    let fields = table.fields();
    let mut out = Vec::new();
    for (index, line) in lines.enumerate() {
        let number = index + 3;
        let row: Vec<&str> = line.split('\t').collect();
        if row.len() != fields.len() {
            return Err(err(
                table,
                number,
                SpecimenErrorKind::ColumnCount {
                    expected: fields.len(),
                    found: row.len(),
                },
            ));
        }
        for (&field, value) in fields.iter().zip(&row) {
            if value.len() > SPECIMEN_FIELD_BYTE_CAP {
                return Err(err(
                    table,
                    number,
                    SpecimenErrorKind::FieldTooLong {
                        field,
                        cap: SPECIMEN_FIELD_BYTE_CAP,
                    },
                ));
            }
            if value.trim().is_empty() {
                return Err(err(table, number, SpecimenErrorKind::MissingField(field)));
            }
            if value.chars().any(char::is_control) {
                return Err(err(table, number, unknown(field, value)));
            }
        }
        out.push((number, row));
    }
    Ok(out)
}

/// The §1.4 reference/load state a description was measured in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceState {
    /// No preload (the reference for loaded equilibria).
    Unloaded,
    /// Strung to playing tension.
    Strung,
    /// Assembled as played, unstrung or without player loads.
    Assembled,
    /// Under player loads.
    Played,
}

impl ReferenceState {
    fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "unloaded" => Self::Unloaded,
            "strung" => Self::Strung,
            "assembled" => Self::Assembled,
            "played" => Self::Played,
            _ => return None,
        })
    }
}

/// `specimen.tsv`: identity, rights and reference state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Specimen {
    entries: BTreeMap<String, String>,
    reference_state: ReferenceState,
    source_ids: Vec<String>,
    digest: ContentHash,
}

impl Specimen {
    /// Parse and require every [`SPECIMEN_KEYS`] entry.
    ///
    /// # Errors
    /// Header, field, duplicate-key, missing-key and vocabulary refusals.
    pub fn parse(text: &str) -> Result<Self, SpecimenError> {
        let table = SpecimenTable::Specimen;
        let mut entries = BTreeMap::new();
        for (line, row) in table_rows(table, text)? {
            if !slug(&row[0].replace('_', "-")) {
                return Err(err(table, line, unknown("key", row[0])));
            }
            if entries
                .insert(row[0].to_owned(), row[1].to_owned())
                .is_some()
            {
                return Err(err(
                    table,
                    line,
                    SpecimenErrorKind::Duplicate(row[0].to_owned()),
                ));
            }
        }
        for key in SPECIMEN_KEYS {
            if !entries.contains_key(key) {
                return Err(err(
                    table,
                    0,
                    SpecimenErrorKind::Missing(format!("specimen key {key}")),
                ));
            }
        }
        if !slug(&entries["specimen_id"]) {
            return Err(err(
                table,
                0,
                unknown("specimen_id", &entries["specimen_id"]),
            ));
        }
        let state = &entries["reference_state"];
        let reference_state = ReferenceState::parse(state)
            .ok_or_else(|| err(table, 0, unknown("reference_state", state)))?;
        let source_ids: Vec<String> = entries["source_ids"]
            .split(';')
            .map(str::to_owned)
            .collect();
        let mut seen = BTreeSet::new();
        for id in &source_ids {
            if !slug(id) || !seen.insert(id) {
                return Err(err(table, 0, unknown("source_ids", id)));
            }
        }
        Ok(Self {
            entries,
            reference_state,
            source_ids,
            digest: hash_bytes(text.as_bytes()),
        })
    }
    /// The specimen slug (its directory name).
    #[must_use]
    pub fn id(&self) -> &str {
        &self.entries["specimen_id"]
    }
    /// Any entry, verbatim.
    #[must_use]
    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries.get(key).map(String::as_str)
    }
    /// The declared reference/load state.
    #[must_use]
    pub fn reference_state(&self) -> ReferenceState {
        self.reference_state
    }
    /// Ledger `source_id`s this specimen draws on.
    #[must_use]
    pub fn source_ids(&self) -> &[String] {
        &self.source_ids
    }
    /// BLAKE3 of the parsed bytes.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        self.digest
    }
}

/// SI unit of a dimension value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiUnit {
    /// Metre.
    Metre,
    /// Square metre.
    SquareMetre,
    /// Cubic metre.
    CubicMetre,
    /// Kilogram.
    Kilogram,
    /// Radian.
    Radian,
    /// Dimensionless.
    One,
}

impl SiUnit {
    fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "m" => Self::Metre,
            "m2" => Self::SquareMetre,
            "m3" => Self::CubicMetre,
            "kg" => Self::Kilogram,
            "rad" => Self::Radian,
            "1" => Self::One,
            _ => return None,
        })
    }
}

/// Plan M6.3 source hierarchy, per value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AcquisitionRoute {
    /// 1: direct measurement of the specimen.
    SpecimenMeasured,
    /// 2: published measured geometry (CT studies, theses, maker drawings).
    PublishedMeasured,
    /// 3: online CAD or scans checked against independent dimensions.
    OnlineCadChecked,
    /// 4: reconstruction from calibrated photographs.
    PhotoReconstruction,
}

impl AcquisitionRoute {
    fn parse(value: &str) -> Option<Self> {
        Some(match value {
            "specimen-measured" => Self::SpecimenMeasured,
            "published-measured" => Self::PublishedMeasured,
            "online-cad-checked" => Self::OnlineCadChecked,
            "photo-reconstruction" => Self::PhotoReconstruction,
            _ => return None,
        })
    }
}

/// Uncertainty of one value: a σ or bound in the value's unit, or `unstated`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Sigma {
    /// Stated standard uncertainty or precision bound (finite, ≥ 0).
    Value(f64),
    /// The source states none.
    Unstated,
}

/// One long-format dimension value.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionRow {
    /// Described part (slug).
    pub part: String,
    /// Feature of the part (slug), e.g. `hole-07`.
    pub feature: String,
    /// Measured quantity (slug), e.g. `chimney-height`.
    pub quantity: String,
    /// Value in SI.
    pub value_si: f64,
    /// Unit of `value_si`.
    pub unit: SiUnit,
    /// Uncertainty in the same unit.
    pub sigma: Sigma,
    /// Acquisition route.
    pub route: AcquisitionRoute,
    /// Ledger source.
    pub source_id: String,
    /// Where in the source, and any caveat.
    pub note: String,
}

type DimensionKey = (String, String, String);

/// `dimensions.tsv`, keyed by `(part, feature, quantity)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Dimensions {
    rows: BTreeMap<DimensionKey, DimensionRow>,
    digest: ContentHash,
}

impl Dimensions {
    /// Parse every row; duplicates refuse.
    ///
    /// # Errors
    /// Header, field, slug, unit, route, number and duplicate refusals.
    pub fn parse(text: &str) -> Result<Self, SpecimenError> {
        let table = SpecimenTable::Dimensions;
        let mut rows = BTreeMap::new();
        for (line, r) in table_rows(table, text)? {
            for (field, value) in [
                ("part", r[0]),
                ("feature", r[1]),
                ("quantity", r[2]),
                ("source_id", r[7]),
            ] {
                if !slug(value) {
                    return Err(err(table, line, unknown(field, value)));
                }
            }
            let value_si: f64 = r[3]
                .parse()
                .ok()
                .filter(|v: &f64| v.is_finite())
                .ok_or_else(|| err(table, line, unknown("value_si", r[3])))?;
            let unit =
                SiUnit::parse(r[4]).ok_or_else(|| err(table, line, unknown("unit", r[4])))?;
            let sigma = if r[5] == "unstated" {
                Sigma::Unstated
            } else {
                r[5].parse()
                    .ok()
                    .filter(|s: &f64| s.is_finite() && *s >= 0.0)
                    .map(Sigma::Value)
                    .ok_or_else(|| err(table, line, unknown("sigma_si", r[5])))?
            };
            let route = AcquisitionRoute::parse(r[6])
                .ok_or_else(|| err(table, line, unknown("route", r[6])))?;
            let key = (r[0].to_owned(), r[1].to_owned(), r[2].to_owned());
            let row = DimensionRow {
                part: key.0.clone(),
                feature: key.1.clone(),
                quantity: key.2.clone(),
                value_si,
                unit,
                sigma,
                route,
                source_id: r[7].to_owned(),
                note: r[8].to_owned(),
            };
            if rows.insert(key, row).is_some() {
                return Err(err(
                    table,
                    line,
                    SpecimenErrorKind::Duplicate(format!("{}/{}/{}", r[0], r[1], r[2])),
                ));
            }
        }
        Ok(Self {
            rows,
            digest: hash_bytes(text.as_bytes()),
        })
    }
    /// Rows in key order.
    pub fn rows(&self) -> impl Iterator<Item = &DimensionRow> {
        self.rows.values()
    }
    /// One value, if present.
    #[must_use]
    pub fn get(&self, part: &str, feature: &str, quantity: &str) -> Option<&DimensionRow> {
        self.rows
            .get(&(part.to_owned(), feature.to_owned(), quantity.to_owned()))
    }
    /// One value; absence refuses naming `part/feature/quantity`.
    ///
    /// # Errors
    /// [`SpecimenErrorKind::Missing`].
    pub fn require(
        &self,
        part: &str,
        feature: &str,
        quantity: &str,
    ) -> Result<&DimensionRow, SpecimenError> {
        self.get(part, feature, quantity).ok_or_else(|| {
            err(
                SpecimenTable::Dimensions,
                0,
                SpecimenErrorKind::Missing(format!("dimension {part}/{feature}/{quantity}")),
            )
        })
    }
    /// Distinct features of a part, sorted.
    #[must_use]
    pub fn features(&self, part: &str) -> Vec<&str> {
        let set: BTreeSet<&str> = self
            .rows
            .values()
            .filter(|r| r.part == part)
            .map(|r| r.feature.as_str())
            .collect();
        set.into_iter().collect()
    }
    /// Distinct parts, sorted.
    #[must_use]
    pub fn parts(&self) -> BTreeSet<&str> {
        self.rows.values().map(|r| r.part.as_str()).collect()
    }
    /// BLAKE3 of the parsed bytes.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        self.digest
    }
}

/// A keyed table whose rows are kept as verbatim fields (frames, parts,
/// unknowns): the first column is a unique slug key.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KeyedTable {
    table: SpecimenTable,
    rows: BTreeMap<String, Vec<String>>,
    digest: ContentHash,
}

impl KeyedTable {
    /// Parse `frames.tsv`, `parts.tsv` or `unknowns.tsv`.
    ///
    /// # Errors
    /// Header, field, slug-key and duplicate refusals; any other table kind
    /// refuses as [`SpecimenErrorKind::Header`].
    pub fn parse(table: SpecimenTable, text: &str) -> Result<Self, SpecimenError> {
        if !matches!(
            table,
            SpecimenTable::Frames | SpecimenTable::Parts | SpecimenTable::Unknowns
        ) {
            return Err(err(table, 0, SpecimenErrorKind::Header));
        }
        let mut rows = BTreeMap::new();
        for (line, r) in table_rows(table, text)? {
            if !slug(r[0]) {
                return Err(err(table, line, unknown(table.fields()[0], r[0])));
            }
            let fields = r.iter().map(|s| (*s).to_owned()).collect();
            if rows.insert(r[0].to_owned(), fields).is_some() {
                return Err(err(
                    table,
                    line,
                    SpecimenErrorKind::Duplicate(r[0].to_owned()),
                ));
            }
        }
        Ok(Self {
            table,
            rows,
            digest: hash_bytes(text.as_bytes()),
        })
    }
    /// The row for a key, as fields in column order.
    #[must_use]
    pub fn row(&self, key: &str) -> Option<&[String]> {
        self.rows.get(key).map(Vec::as_slice)
    }
    /// One field of a keyed row, by column name.
    #[must_use]
    pub fn field(&self, key: &str, column: &str) -> Option<&str> {
        let index = self.table.fields().iter().position(|c| *c == column)?;
        self.rows.get(key).map(|r| r[index].as_str())
    }
    /// Keys, sorted.
    pub fn keys(&self) -> impl Iterator<Item = &str> {
        self.rows.keys().map(String::as_str)
    }
    /// BLAKE3 of the parsed bytes.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        self.digest
    }
}

/// One hole's state in a fingering: `fraction` 1 is open, `0 < σ < 1` a vent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoleOpening {
    /// Hole index in the specimen's declared numbering (≥ 1).
    pub hole: u16,
    /// Opening fraction σ ∈ (0, 1].
    pub fraction: f64,
}

/// One named fingering; every hole not listed is closed.
#[derive(Debug, Clone, PartialEq)]
pub struct FingeringRow {
    /// Unique label (slug).
    pub fingering: String,
    /// Written note.
    pub written: String,
    /// Concert (sounding) note.
    pub concert: String,
    /// Register label.
    pub register: String,
    /// Open or vented holes, strictly increasing by index.
    pub open_holes: Vec<HoleOpening>,
    /// Ledger source of the chart.
    pub source_id: String,
    /// Source location and mapping notes.
    pub note: String,
}

/// `fingerings.tsv`, in file order (the chart order); labels are unique.
#[derive(Debug, Clone, PartialEq)]
pub struct Fingerings {
    rows: Vec<FingeringRow>,
    digest: ContentHash,
}

fn parse_openings(value: &str) -> Option<Vec<HoleOpening>> {
    if value == "-" {
        return Some(Vec::new());
    }
    let mut out: Vec<HoleOpening> = Vec::new();
    for item in value.split(',') {
        let (hole, fraction) = match item.split_once('@') {
            Some((h, s)) => (h, s.parse::<f64>().ok().filter(|s| *s > 0.0 && *s < 1.0)?),
            None => (item, 1.0),
        };
        let hole: u16 = hole.parse().ok().filter(|h| *h >= 1)?;
        if out.last().is_some_and(|last| last.hole >= hole) {
            return None;
        }
        out.push(HoleOpening { hole, fraction });
    }
    Some(out)
}

impl Fingerings {
    /// Parse the chart. `open_holes` is `-` (all closed) or comma-separated
    /// strictly increasing indices, each optionally `@σ` for a vent.
    ///
    /// # Errors
    /// Header, field, slug, hole-list and duplicate-label refusals.
    pub fn parse(text: &str) -> Result<Self, SpecimenError> {
        let table = SpecimenTable::Fingerings;
        let mut rows = Vec::new();
        let mut labels = BTreeSet::new();
        for (line, r) in table_rows(table, text)? {
            for (field, value) in [("fingering", r[0]), ("source_id", r[5])] {
                if !slug(value) {
                    return Err(err(table, line, unknown(field, value)));
                }
            }
            let open_holes = parse_openings(r[4])
                .ok_or_else(|| err(table, line, unknown("open_holes", r[4])))?;
            if !labels.insert(r[0]) {
                return Err(err(
                    table,
                    line,
                    SpecimenErrorKind::Duplicate(r[0].to_owned()),
                ));
            }
            rows.push(FingeringRow {
                fingering: r[0].to_owned(),
                written: r[1].to_owned(),
                concert: r[2].to_owned(),
                register: r[3].to_owned(),
                open_holes,
                source_id: r[5].to_owned(),
                note: r[6].to_owned(),
            });
        }
        Ok(Self {
            rows,
            digest: hash_bytes(text.as_bytes()),
        })
    }
    /// Rows in chart order.
    #[must_use]
    pub fn rows(&self) -> &[FingeringRow] {
        &self.rows
    }
    /// BLAKE3 of the parsed bytes.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        self.digest
    }
}

/// A cross-checked specimen directory.
#[derive(Debug, Clone, PartialEq)]
pub struct SpecimenDescription {
    /// `specimen.tsv`.
    pub specimen: Specimen,
    /// `frames.tsv`.
    pub frames: KeyedTable,
    /// `dimensions.tsv`.
    pub dimensions: Dimensions,
    /// `parts.tsv`.
    pub parts: KeyedTable,
    /// `unknowns.tsv`.
    pub unknowns: KeyedTable,
    /// `fingerings.tsv`, when the specimen has one.
    pub fingerings: Option<Fingerings>,
}

impl SpecimenDescription {
    /// Admit parsed tables that agree with each other: every dimensioned part
    /// has a frame and a parts row, and every cited source is listed in
    /// `specimen.tsv`.
    ///
    /// # Errors
    /// [`SpecimenErrorKind::Inconsistent`] naming the first disagreement, or
    /// a table passed under the wrong kind.
    pub fn try_new(
        specimen: Specimen,
        frames: KeyedTable,
        dimensions: Dimensions,
        parts: KeyedTable,
        unknowns: KeyedTable,
        fingerings: Option<Fingerings>,
    ) -> Result<Self, SpecimenError> {
        for (got, want) in [
            (frames.table, SpecimenTable::Frames),
            (parts.table, SpecimenTable::Parts),
            (unknowns.table, SpecimenTable::Unknowns),
        ] {
            if got != want {
                return Err(err(want, 0, SpecimenErrorKind::Header));
            }
        }
        let inconsistent =
            |table, what: String| err(table, 0, SpecimenErrorKind::Inconsistent(what));
        for part in dimensions.parts() {
            if frames.row(part).is_none() {
                return Err(inconsistent(
                    SpecimenTable::Frames,
                    format!("part {part} has no frame"),
                ));
            }
            if parts.row(part).is_none() {
                return Err(inconsistent(
                    SpecimenTable::Parts,
                    format!("part {part} has no parts row"),
                ));
            }
        }
        let listed: BTreeSet<&str> = specimen.source_ids().iter().map(String::as_str).collect();
        let cited = dimensions
            .rows()
            .map(|r| (SpecimenTable::Dimensions, r.source_id.as_str()))
            .chain(fingerings.iter().flat_map(|f| {
                f.rows()
                    .iter()
                    .map(|r| (SpecimenTable::Fingerings, r.source_id.as_str()))
            }));
        for (table, source) in cited {
            if !listed.contains(source) {
                return Err(inconsistent(
                    table,
                    format!("source {source} is not listed in specimen.tsv"),
                ));
            }
        }
        Ok(Self {
            specimen,
            frames,
            dimensions,
            parts,
            unknowns,
            fingerings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc(table: SpecimenTable, rows: &[&str]) -> String {
        let mut out = format!("{}\n{}\n", table.schema(), table.fields().join("\t"));
        for row in rows {
            out.push_str(row);
            out.push('\n');
        }
        out
    }
    const DIM: &str =
        "bore\tstation-01\tdiameter\t0.01465\tm\t1e-05\tpublished-measured\tsrc-a\tTable 3";
    fn specimen_rows(source_ids: &str) -> Vec<String> {
        SPECIMEN_KEYS
            .iter()
            .map(|k| match *k {
                "specimen_id" => "specimen_id\tclarinet-test".to_owned(),
                "reference_state" => "reference_state\tassembled".to_owned(),
                "source_ids" => format!("source_ids\t{source_ids}"),
                key => format!("{key}\tunstated"),
            })
            .collect()
    }
    fn specimen(source_ids: &str) -> Specimen {
        let rows = specimen_rows(source_ids);
        let refs: Vec<&str> = rows.iter().map(String::as_str).collect();
        Specimen::parse(&doc(SpecimenTable::Specimen, &refs)).expect("specimen")
    }
    fn keyed(table: SpecimenTable, rows: &[&str]) -> KeyedTable {
        KeyedTable::parse(table, &doc(table, rows)).expect("keyed")
    }

    #[test]
    fn headers_columns_and_blank_fields_refuse() {
        let t = SpecimenTable::Dimensions;
        let bad_schema = doc(t, &[DIM]).replacen("v1", "v2", 1);
        assert_eq!(
            Dimensions::parse(&bad_schema).unwrap_err().kind,
            SpecimenErrorKind::Header
        );
        let short = doc(t, &["bore\tstation-01"]);
        assert_eq!(
            Dimensions::parse(&short).unwrap_err().kind,
            SpecimenErrorKind::ColumnCount {
                expected: 9,
                found: 2
            }
        );
        let blank = doc(t, &[&DIM.replace("Table 3", " ")]);
        let e = Dimensions::parse(&blank).unwrap_err();
        assert_eq!(
            (e.line, e.kind),
            (3, SpecimenErrorKind::MissingField("note"))
        );
    }

    #[test]
    fn field_byte_cap_admits_cap_and_refuses_cap_plus_one() {
        let t = SpecimenTable::Dimensions;
        let at_cap = DIM.replace("Table 3", &"x".repeat(SPECIMEN_FIELD_BYTE_CAP));
        assert!(Dimensions::parse(&doc(t, &[&at_cap])).is_ok());
        let over = DIM.replace("Table 3", &"x".repeat(SPECIMEN_FIELD_BYTE_CAP + 1));
        assert!(matches!(
            Dimensions::parse(&doc(t, &[&over])).unwrap_err().kind,
            SpecimenErrorKind::FieldTooLong { field: "note", .. }
        ));
    }

    #[test]
    fn units_routes_numbers_and_sigma_refuse_unknowns() {
        let t = SpecimenTable::Dimensions;
        for (from, to, field) in [
            ("\tm\t", "\tmm\t", "unit"),
            ("published-measured", "guessed", "route"),
            ("0.01465", "NaN", "value_si"),
            ("0.01465", "inf", "value_si"),
            ("1e-05", "-1e-05", "sigma_si"),
            ("1e-05", "", "sigma_si"),
            ("station-01", "Station 1", "feature"),
        ] {
            let row = DIM.replacen(from, to, 1);
            let kind = Dimensions::parse(&doc(t, &[&row])).unwrap_err().kind;
            match kind {
                SpecimenErrorKind::UnknownValue { field: f, .. } => assert_eq!(f, field, "{row}"),
                SpecimenErrorKind::MissingField(f) => assert_eq!(f, field, "{row}"),
                other => panic!("{row}: {other:?}"),
            }
        }
        let unstated = DIM.replace("1e-05", "unstated");
        let parsed = Dimensions::parse(&doc(t, &[&unstated])).expect("unstated is explicit");
        assert_eq!(parsed.rows().next().unwrap().sigma, Sigma::Unstated);
    }

    #[test]
    fn duplicates_refuse_and_require_names_the_absent_value() {
        let t = SpecimenTable::Dimensions;
        assert!(matches!(
            Dimensions::parse(&doc(t, &[DIM, DIM])).unwrap_err().kind,
            SpecimenErrorKind::Duplicate(_)
        ));
        let dims = Dimensions::parse(&doc(t, &[DIM])).unwrap();
        assert_eq!(
            dims.require("bore", "station-01", "diameter")
                .unwrap()
                .value_si,
            0.01465
        );
        assert_eq!(
            dims.require("bore", "station-02", "diameter")
                .unwrap_err()
                .kind,
            SpecimenErrorKind::Missing("dimension bore/station-02/diameter".to_owned())
        );
    }

    #[test]
    fn row_permutation_is_lookup_invariant_and_digest_tracks_bytes() {
        let t = SpecimenTable::Dimensions;
        let second = DIM
            .replace("station-01", "station-02")
            .replace("0.01465", "0.01470");
        let a = Dimensions::parse(&doc(t, &[DIM, &second])).unwrap();
        let b = Dimensions::parse(&doc(t, &[&second, DIM])).unwrap();
        assert_eq!(a.rows().collect::<Vec<_>>(), b.rows().collect::<Vec<_>>());
        assert_eq!(a.features("bore"), vec!["station-01", "station-02"]);
        assert_ne!(
            a.digest(),
            b.digest(),
            "the digest identifies bytes, not content"
        );
    }

    #[test]
    fn specimen_requires_every_key_and_known_state() {
        let mut rows = specimen_rows("src-a");
        rows.retain(|r| !r.starts_with("holder\t"));
        let refs: Vec<&str> = rows.iter().map(String::as_str).collect();
        assert_eq!(
            Specimen::parse(&doc(SpecimenTable::Specimen, &refs))
                .unwrap_err()
                .kind,
            SpecimenErrorKind::Missing("specimen key holder".to_owned())
        );
        let rows = specimen_rows("src-a");
        let swapped: Vec<String> = rows
            .iter()
            .map(|r| r.replace("reference_state\tassembled", "reference_state\tloaded"))
            .collect();
        let refs: Vec<&str> = swapped.iter().map(String::as_str).collect();
        assert!(matches!(
            Specimen::parse(&doc(SpecimenTable::Specimen, &refs))
                .unwrap_err()
                .kind,
            SpecimenErrorKind::UnknownValue {
                field: "reference_state",
                ..
            }
        ));
        let s = specimen("src-a;src-b");
        assert_eq!(
            (s.id(), s.reference_state()),
            ("clarinet-test", ReferenceState::Assembled)
        );
        assert_eq!(s.source_ids(), ["src-a", "src-b"]);
    }

    #[test]
    fn fingering_lists_parse_vents_and_refuse_malformed_states() {
        let row = |holes: &str| format!("a3\tA3\tG3\tfirst\t{holes}\tsrc-a\tchart");
        let chart =
            Fingerings::parse(&doc(SpecimenTable::Fingerings, &[&row("1,3,5@0.25,24")])).unwrap();
        let holes = &chart.rows()[0].open_holes;
        assert_eq!(
            holes.iter().map(|h| h.hole).collect::<Vec<_>>(),
            [1, 3, 5, 24]
        );
        assert_eq!((holes[2].fraction, holes[3].fraction), (0.25, 1.0));
        let closed = Fingerings::parse(&doc(SpecimenTable::Fingerings, &[&row("-")])).unwrap();
        assert!(closed.rows()[0].open_holes.is_empty());
        for bad in ["0", "3,1", "1,1", "2@1", "2@0", "2@x", "", "1,,2"] {
            let parsed = Fingerings::parse(&doc(SpecimenTable::Fingerings, &[&row(bad)]));
            assert!(parsed.is_err(), "{bad:?} must refuse");
        }
        let twice = row("1");
        assert!(matches!(
            Fingerings::parse(&doc(SpecimenTable::Fingerings, &[&twice, &twice]))
                .unwrap_err()
                .kind,
            SpecimenErrorKind::Duplicate(_)
        ));
    }

    #[test]
    fn description_refuses_a_part_without_frame_or_an_unlisted_source() {
        let frames = |part: &str| {
            keyed(
                SpecimenTable::Frames,
                &[&format!("{part}\torigin\t+x\tunstated\tunstated\tnote")],
            )
        };
        let parts = keyed(
            SpecimenTable::Parts,
            &["bore\tunstated\tassembled\t-\tnote"],
        );
        let unknowns = keyed(
            SpecimenTable::Unknowns,
            &["u-a\twhat\tunbounded\tg\tclaims\troute"],
        );
        let dims = || Dimensions::parse(&doc(SpecimenTable::Dimensions, &[DIM])).unwrap();
        let ok = SpecimenDescription::try_new(
            specimen("src-a"),
            frames("bore"),
            dims(),
            parts.clone(),
            unknowns.clone(),
            None,
        );
        assert!(ok.is_ok());
        let no_frame = SpecimenDescription::try_new(
            specimen("src-a"),
            frames("barrel"),
            dims(),
            parts.clone(),
            unknowns.clone(),
            None,
        );
        assert_eq!(
            no_frame.unwrap_err().kind,
            SpecimenErrorKind::Inconsistent("part bore has no frame".to_owned())
        );
        let unlisted = SpecimenDescription::try_new(
            specimen("src-b"),
            frames("bore"),
            dims(),
            parts.clone(),
            unknowns.clone(),
            None,
        );
        assert!(matches!(
            unlisted.unwrap_err().kind,
            SpecimenErrorKind::Inconsistent(_)
        ));
        assert!(KeyedTable::parse(SpecimenTable::Dimensions, "").is_err());
    }
}
