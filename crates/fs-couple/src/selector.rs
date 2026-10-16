//! The selector (music bead `frankensim-music-v8-root-3ez8g.14`).
//!
//! Gesture × claim × budget → image set. A filling with a menu needs a
//! chooser: which image(s) run THIS block given the gesture phase (phrase vs
//! held), the claim being served, and the CPU budget. Without it menus are
//! documentation; with it they are a runtime capability no sampler can
//! imitate.
//!
//! DOCTRINE ENFORCED HERE (as admission rules, not prose):
//!
//! - D19: claims are `(filling, image, qoi)` rows in
//!   `instrument-claims.json`. The selector consults a parsed snapshot of
//!   that registry; a menu entry with no admitting row does not exist.
//! - D25: a live-default image requires `gate == green` AND a budget row.
//!   An UNGATED image is never selectable as a live default. A REFUSED
//!   image is never selectable at all.
//! - CPU panic drops to a PRE-GATED cheaper image of the SAME filling,
//!   never to a forbidden or off-menu image (structural: candidates come
//!   only from admitted menu entries).
//! - Mid-note hops require a supported D17 state lift on the target entry
//!   AND the articulation bead's measured settle policy
//!   (`data/claims/wind-hop-policy.tsv`, minted by
//!   `wind_hop_tests::mint_hop_policy_artifact`); otherwise the hop is
//!   deferred to the next phrase boundary.
//! - Peak-location questions route to the FD oracle WITHOUT rendering audio.
//!
//! NO PHYSICS LIVES HERE: the selector chooses, it never computes sound.
//! Determinism contract: identical (registry snapshot, menu, policy,
//! request) tuples produce byte-identical receipts — decisions iterate the
//! caller-declared menu order only, never hash maps.
//!
//! Serving scope note: `Play` and `PeakLocations` are the servings the
//! DONE-WHEN demos name. "Design" and "lock-question" servings join when
//! their consuming beads land; the enum is the extension point, not a stub
//! (every variant here is exercised by tests).

use fs_blake3::{ContentHash, hash_domain};

/// Schema string stamped on every emitted selection receipt.
pub const SELECTION_RECEIPT_SCHEMA: &str = "frankensim-selection-receipt-v1";

/// Registry schema this module parses. Mirrors
/// `xtask::instrument_claims::REGISTRY_SCHEMA`; kept literal because
/// fs-couple (L3) cannot depend on xtask (TOOL layer).
pub const REGISTRY_SCHEMA_V1: &str = "frankensim-instrument-claims-v1";

/// Hash domain for selection-receipt content identities.
const RECEIPT_HASH_DOMAIN: &str = "frankensim-fs-couple-selection-receipt-v1";

/// Typed gate status of a registry row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Gate {
    /// Implemented but no gates-bead evidence yet.
    Ungated,
    /// Gates bead landed evidence; selectable.
    Green,
    /// Deliberately refused (D21: rows transition, never vanish).
    Refused,
}

impl Gate {
    fn parse(value: &str) -> Result<Self, SelectorError> {
        match value {
            "ungated" => Ok(Self::Ungated),
            "green" => Ok(Self::Green),
            "refused" => Ok(Self::Refused),
            other => Err(SelectorError::Registry(format!(
                "invalid gate value {other:?}"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Ungated => "ungated",
            Self::Green => "green",
            Self::Refused => "refused",
        }
    }
}

/// One `(filling, image, qoi)` claim row, reduced to what selection reasons
/// about.
#[derive(Clone, Debug, PartialEq)]
pub struct RegistryRow {
    /// Instrument filling, e.g. `wind-reed`.
    pub filling: String,
    /// Performance image, e.g. `char-line`.
    pub image: String,
    /// Quality of interest served, e.g. `quarter-wave-lock`.
    pub qoi: String,
    /// Gate status.
    pub gate: Gate,
    /// Registry `live_default` flag (`yes`).
    pub live_default: bool,
    /// Budget-row reference when present (D25's second live-default limb).
    pub budget_row: Option<String>,
}

/// A parsed snapshot of `instrument-claims.json`.
///
/// Construction validates the schema tag and every consumed field, refusing
/// malformed input instead of silently skipping rows: one bad row must not
/// hide from an admission decision.
#[derive(Clone, Debug, PartialEq)]
pub struct RegistrySnapshot {
    rows: Vec<RegistryRow>,
}

impl RegistrySnapshot {
    /// Parse the tracked registry source.
    ///
    /// # Errors
    /// [`SelectorError::Registry`] for schema mismatch, structural defects,
    /// or invalid field values.
    pub fn parse(source: &str) -> Result<Self, SelectorError> {
        let rows = RegistryJsonParser::new(source)?.parse_rows()?;
        Ok(Self { rows })
    }

    /// All rows for a `(filling, image)` pair.
    #[must_use]
    pub fn rows_for(&self, filling: &str, image: &str) -> Vec<&RegistryRow> {
        self.rows
            .iter()
            .filter(|r| r.filling == filling && r.image == image)
            .collect()
    }

    /// The row governing one exact claim.
    #[must_use]
    pub fn row(&self, filling: &str, image: &str, qoi: &str) -> Option<&RegistryRow> {
        self.rows
            .iter()
            .find(|r| r.filling == filling && r.image == image && r.qoi == qoi)
    }

    /// Row count (diagnostics).
    #[must_use]
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// Whether the snapshot holds no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}

/// Minimal JSON-subset parser for exactly the registry fields selection
/// consumes. fs-couple has no JSON dependency (Franken-only graph); each
/// crate hand-rolls the narrow parser it needs (fs-ledger `RowJsonParser`
/// precedent). Unknown members are skipped, never interpreted.
struct RegistryJsonParser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> RegistryJsonParser<'a> {
    /// Enter the top object and validate the schema tag; cursor lands just
    /// past `"schema":"..."`.
    fn new(source: &'a str) -> Result<Self, SelectorError> {
        let mut p = Self {
            bytes: source.as_bytes(),
            pos: 0,
        };
        p.skip_ws();
        p.expect_byte(b'{')?;
        let schema = p.string_field("schema")?;
        if schema != REGISTRY_SCHEMA_V1 {
            return Err(SelectorError::Registry(format!(
                "registry schema {schema:?} is not {REGISTRY_SCHEMA_V1:?}"
            )));
        }
        Ok(p)
    }

    fn err_at(&self, what: &str) -> SelectorError {
        SelectorError::Registry(format!("{what} at byte {}", self.pos))
    }

    fn skip_ws(&mut self) {
        while self.pos < self.bytes.len()
            && matches!(self.bytes[self.pos], b' ' | b'\t' | b'\n' | b'\r')
        {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn eat_byte(&mut self, expected: u8) -> bool {
        if self.peek() == Some(expected) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expect_byte(&mut self, expected: u8) -> Result<(), SelectorError> {
        self.skip_ws();
        if self.eat_byte(expected) {
            Ok(())
        } else {
            Err(self.err_at(&format!("expected {:?}", expected as char)))
        }
    }

    fn parse_string(&mut self) -> Result<String, SelectorError> {
        self.skip_ws();
        if !self.eat_byte(b'"') {
            return Err(self.err_at("expected string"));
        }
        let mut out = String::new();
        loop {
            match self.peek() {
                None => return Err(self.err_at("unterminated string")),
                Some(b'"') => {
                    self.pos += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    self.pos += 1;
                    let esc = self.peek().ok_or_else(|| self.err_at("truncated escape"))?;
                    self.pos += 1;
                    match esc {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'n' => out.push('\n'),
                        b't' => out.push('\t'),
                        b'r' => out.push('\r'),
                        // \uXXXX (with surrogate-pair support): the tracked
                        // registry legitimately contains \u2014 em-dashes in
                        // note prose, so this is the real path, not a
                        // theoretical one.
                        b'u' => {
                            let cp = self.parse_hex4()?;
                            let decoded = if (0xD800..=0xDBFF).contains(&cp) {
                                if self.peek() != Some(b'\\') {
                                    return Err(self.err_at("lone high surrogate"));
                                }
                                self.pos += 1;
                                if self.peek() != Some(b'u') {
                                    return Err(self.err_at("lone high surrogate"));
                                }
                                self.pos += 1;
                                let low = self.parse_hex4()?;
                                if !(0xDC00..=0xDFFF).contains(&low) {
                                    return Err(self.err_at("invalid low surrogate"));
                                }
                                0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00)
                            } else {
                                cp
                            };
                            match char::from_u32(decoded) {
                                Some(c) => out.push(c),
                                None => {
                                    return Err(self.err_at("invalid \\u escape"));
                                }
                            }
                        }
                        other => {
                            return Err(self.err_at(&format!("unknown escape {:?}", other as char)));
                        }
                    }
                }
                Some(byte) => {
                    out.push(byte as char);
                    self.pos += 1;
                }
            }
        }
    }

    /// Read exactly four hex digits (the payload of a `\u` escape).
    fn parse_hex4(&mut self) -> Result<u32, SelectorError> {
        let mut value: u32 = 0;
        for _ in 0..4 {
            let byte = self.peek().ok_or_else(|| self.err_at("truncated \\u"))?;
            let digit = (byte as char)
                .to_digit(16)
                .ok_or_else(|| self.err_at("bad hex in \\u"))?;
            value = value * 16 + digit;
            self.pos += 1;
        }
        Ok(value)
    }

    fn skip_value(&mut self) -> Result<(), SelectorError> {
        self.skip_ws();
        match self.peek() {
            Some(b'"') => {
                self.parse_string()?;
                Ok(())
            }
            Some(open @ (b'{' | b'[')) => {
                let close = if open == b'{' { b'}' } else { b']' };
                let mut depth = 0usize;
                loop {
                    match self.peek() {
                        None => return Err(self.err_at("unterminated container")),
                        Some(b'"') => {
                            self.parse_string()?;
                        }
                        Some(b) if b == open => {
                            depth += 1;
                            self.pos += 1;
                        }
                        Some(b) if b == close => {
                            depth -= 1;
                            self.pos += 1;
                            if depth == 0 {
                                return Ok(());
                            }
                        }
                        Some(_) => self.pos += 1,
                    }
                }
            }
            Some(_) => {
                while let Some(b) = self.peek() {
                    if matches!(b, b',' | b'}' | b']') {
                        break;
                    }
                    self.pos += 1;
                }
                Ok(())
            }
            None => Err(self.err_at("expected value")),
        }
    }

    /// Read `"name":"<string>"` at the cursor.
    fn string_field(&mut self, name: &str) -> Result<String, SelectorError> {
        self.skip_ws();
        let key = self.parse_string()?;
        if key != name {
            return Err(self.err_at(&format!("expected key {name:?}, got {key:?}")));
        }
        self.expect_byte(b':')?;
        self.parse_string()
    }

    /// If the next member's key is `name`, consume `key:` and return true
    /// with the cursor on the value. Restores the cursor otherwise.
    fn key_field(&mut self, name: &str) -> Result<bool, SelectorError> {
        let save = self.pos;
        self.skip_ws();
        if self.peek() != Some(b'"') {
            return Ok(false);
        }
        match self.parse_string() {
            Ok(key) if key == name => {
                self.expect_byte(b':')?;
                Ok(true)
            }
            _ => {
                self.pos = save;
                Ok(false)
            }
        }
    }

    fn parse_rows(&mut self) -> Result<Vec<RegistryRow>, SelectorError> {
        loop {
            self.skip_ws();
            if self.eat_byte(b',') {
                continue;
            }
            if self.eat_byte(b'}') {
                return Err(SelectorError::Registry("registry has no rows".into()));
            }
            if self.key_field("rows")? {
                return self.parse_row_array();
            }
            self.skip_ws();
            let _key = self.parse_string()?;
            self.expect_byte(b':')?;
            self.skip_value()?;
        }
    }

    fn parse_row_array(&mut self) -> Result<Vec<RegistryRow>, SelectorError> {
        self.skip_ws();
        if !self.eat_byte(b'[') {
            return Err(self.err_at("rows is not an array"));
        }
        let mut rows = Vec::new();
        loop {
            self.skip_ws();
            if self.eat_byte(b']') {
                return Ok(rows);
            }
            if !rows.is_empty() {
                self.expect_byte(b',')?;
                self.skip_ws();
                if self.eat_byte(b']') {
                    return Ok(rows);
                }
            }
            self.expect_byte(b'{')?;
            rows.push(self.parse_row_object()?);
        }
    }

    fn parse_row_object(&mut self) -> Result<RegistryRow, SelectorError> {
        let mut filling = None;
        let mut image = None;
        let mut qoi = None;
        let mut gate = None;
        let mut live_default = false;
        let mut budget_row = None;
        let mut first = true;
        loop {
            self.skip_ws();
            if self.eat_byte(b'}') {
                break;
            }
            if !first {
                self.expect_byte(b',')?;
                self.skip_ws();
                if self.eat_byte(b'}') {
                    break;
                }
            }
            first = false;
            self.skip_ws();
            let key = self.parse_string()?;
            self.expect_byte(b':')?;
            match key.as_str() {
                "filling" => filling = Some(self.parse_string()?),
                "image" => image = Some(self.parse_string()?),
                "qoi" => qoi = Some(self.parse_string()?),
                "gate" => {
                    let g = self.parse_string()?;
                    gate = Some(Gate::parse(&g)?);
                }
                "live_default" => {
                    let v = self.parse_string()?;
                    match v.as_str() {
                        "yes" => live_default = true,
                        "no" => live_default = false,
                        other => {
                            return Err(SelectorError::Registry(format!(
                                "invalid live_default value {other:?}"
                            )));
                        }
                    }
                }
                "budget_row" => {
                    self.skip_ws();
                    if self.peek() == Some(b'n') {
                        self.pos += 1;
                        if self.bytes.get(self.pos..self.pos + 3) != Some(b"ull") {
                            return Err(self.err_at("malformed null"));
                        }
                        self.pos += 3;
                    } else {
                        budget_row = Some(self.parse_string()?);
                    }
                }
                _ => self.skip_value()?,
            }
        }
        let missing = |f: &str| SelectorError::Registry(format!("row missing field {f:?}"));
        Ok(RegistryRow {
            filling: filling.ok_or_else(|| missing("filling"))?,
            image: image.ok_or_else(|| missing("image"))?,
            qoi: qoi.ok_or_else(|| missing("qoi"))?,
            gate: gate.ok_or_else(|| missing("gate"))?,
            live_default,
            budget_row,
        })
    }
}

/// Which claim is the block serving?
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Serving {
    /// Render audio for performance.
    Play,
    /// Answer "where are the peaks?" — FD oracle, NO audio.
    PeakLocations,
}

impl Serving {
    fn as_str(self) -> &'static str {
        match self {
            Self::Play => "play",
            Self::PeakLocations => "peak-locations",
        }
    }

    fn parse(value: &str) -> Result<Self, SelectorError> {
        match value {
            "play" => Ok(Self::Play),
            "peak-locations" => Ok(Self::PeakLocations),
            other => Err(SelectorError::Receipt(format!("invalid serving {other:?}"))),
        }
    }
}

/// Gesture phase, decided by the MEASURED settle policy (never ad-hoc).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GesturePhase {
    /// Attack / moving gesture — spatial-or-island images.
    Phrase,
    /// Settled — the hold image (cheaper) when gated.
    Held,
}

impl GesturePhase {
    fn as_str(self) -> &'static str {
        match self {
            Self::Phrase => "phrase",
            Self::Held => "held",
        }
    }

    fn parse(value: &str) -> Result<Self, SelectorError> {
        match value {
            "phrase" => Ok(Self::Phrase),
            "held" => Ok(Self::Held),
            other => Err(SelectorError::Receipt(format!(
                "invalid gesture phase {other:?}"
            ))),
        }
    }
}

/// Role a menu entry plays in the routing policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EntryClass {
    /// Spatial/island image (char lines, pHS chains): phrase riding.
    SpatialOrIsland,
    /// Settled-note hold image (VFIT hold, modal 1-port).
    SettledHold,
    /// FD spectral oracle: answers peak questions, renders nothing.
    SpectralOracle,
}

/// One caller-declared candidate image on a filling's menu.
///
/// The menu declares WHICH images exist for this voice and their relative
/// cost; the REGISTRY decides admission. Neither alone suffices.
#[derive(Clone, Debug)]
pub struct MenuEntry {
    /// Registry image key.
    pub image: &'static str,
    /// Registry qoi whose gate admits this use.
    pub qoi: &'static str,
    /// Relative CPU cost tier (lower is cheaper). Budget squeezing picks the
    /// cheapest ADMITTED candidate.
    pub cost_tier: u8,
    /// Routing role.
    pub class: EntryClass,
    /// Does the target expose a D17 state lift (mid-note hop support)?
    pub lift_supported: bool,
}

/// A filling's declared image menu (caller data, e.g. from a voice card).
#[derive(Clone, Debug)]
pub struct VoiceMenu {
    /// Registry filling key.
    pub filling: &'static str,
    /// Candidates, in deterministic preference-scan order.
    pub entries: Vec<MenuEntry>,
}

/// The measured hop policy, parsed from the articulation bead's committed
/// artifact `data/claims/wind-hop-policy.tsv` (header line:
/// `# settle detector: relative windowed-RMS drift < 0.05 for 4 consecutive
/// 25ms blocks; first settled block (control run): 6`).
///
/// PARAMETERS ARE DATA: this type carries them verbatim; changing the
/// measurement means re-minting the artifact, not editing code.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HopPolicy {
    /// Relative windowed-RMS drift below which a block counts as settled.
    pub drift_threshold: f64,
    /// Consecutive settled blocks required to declare Held.
    pub consecutive_blocks: u32,
    /// Block length in milliseconds used by the measurement.
    pub window_ms: f64,
    /// First control-block index eligible for a hop (measured control run).
    pub first_settled_block: u32,
}

impl HopPolicy {
    /// Parse the committed artifact's policy header.
    ///
    /// # Errors
    /// [`SelectorError::Policy`] when the header is absent or malformed.
    pub fn parse(tsv: &str) -> Result<Self, SelectorError> {
        let header = tsv
            .lines()
            .find(|l| l.starts_with("# settle detector:"))
            .ok_or_else(|| SelectorError::Policy("no '# settle detector:' header".into()))?;
        let grab = |pat: &str| -> Result<f64, SelectorError> {
            let idx = header
                .find(pat)
                .ok_or_else(|| SelectorError::Policy(format!("policy header missing {pat:?}")))?;
            let tail = &header[idx + pat.len()..];
            let num: String = tail
                .chars()
                .take_while(|c| c.is_ascii_digit() || *c == '.')
                .collect();
            num.parse::<f64>().map_err(|_| {
                SelectorError::Policy(format!("unparsable number after {pat:?}: {num:?}"))
            })
        };
        let drift_threshold = grab("drift < ")?;
        let consecutive_blocks = grab("for ")?.round() as u32;
        let window_ms = grab("consecutive ")?;
        let first_settled_block = header
            .split("(control run): ")
            .nth(1)
            .and_then(|tail| {
                let digits: String = tail.chars().take_while(char::is_ascii_digit).collect();
                digits.parse::<u32>().ok()
            })
            .ok_or_else(|| SelectorError::Policy("missing '(control run): N'".into()))?;
        Ok(Self {
            drift_threshold,
            consecutive_blocks,
            window_ms,
            first_settled_block,
        })
    }

    /// Classify detector outputs into a gesture phase. Pure and
    /// deterministic: same numbers in, same phase out.
    #[must_use]
    pub fn classify(
        &self,
        rel_drift: f64,
        consecutive_settled: u32,
        block_index: u32,
    ) -> GesturePhase {
        let quiet = rel_drift < self.drift_threshold;
        let streak = consecutive_settled >= self.consecutive_blocks;
        let eligible = block_index >= self.first_settled_block;
        if quiet && streak && eligible {
            GesturePhase::Held
        } else {
            GesturePhase::Phrase
        }
    }
}

/// Per-entry admission verdict recorded in the receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct ConsideredEntry {
    /// Menu entry image.
    pub image: String,
    /// Consulted qoi.
    pub qoi: String,
    /// Registry gate string consulted (`absent` when no row exists).
    pub gate: String,
    /// Cost tier from the menu.
    pub cost_tier: u8,
    /// Why this entry was or was not chosen (stable machine words).
    pub verdict: String,
}

/// Whether a mid-note hop to the chosen image is permitted now.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HopVerdict {
    /// Chosen image equals the current one.
    NoChange,
    /// At a phrase boundary: hop freely.
    AtBoundary,
    /// Mid-note but the D17 lift supports it and the phase policy admits it.
    MidNoteAllowed,
    /// Mid-note and policy defers it to the next phrase boundary.
    DeferredToBoundary,
}

impl HopVerdict {
    fn as_str(self) -> &'static str {
        match self {
            Self::NoChange => "no-change",
            Self::AtBoundary => "at-boundary",
            Self::MidNoteAllowed => "mid-note-allowed",
            Self::DeferredToBoundary => "deferred-to-boundary",
        }
    }

    fn parse(value: &str) -> Result<Self, SelectorError> {
        match value {
            "no-change" => Ok(Self::NoChange),
            "at-boundary" => Ok(Self::AtBoundary),
            "mid-note-allowed" => Ok(Self::MidNoteAllowed),
            "deferred-to-boundary" => Ok(Self::DeferredToBoundary),
            other => Err(SelectorError::Receipt(format!(
                "invalid hop verdict {other:?}"
            ))),
        }
    }
}

/// One decision's full, inspectable record — the artifact that answers
/// "why did it pick that image at t=3.2s" without re-running anything.
#[derive(Clone, Debug, PartialEq)]
pub struct SelectionReceipt {
    /// Schema tag.
    pub schema: String,
    /// Filling consulted.
    pub filling: String,
    /// Claim served.
    pub serving: Serving,
    /// Detector-derived gesture phase.
    pub phase: GesturePhase,
    /// Affordable cost tier from the budget config.
    pub budget_headroom_tier: u8,
    /// Chosen `(image, qoi)` when a candidate admitted.
    pub chosen: Option<(String, String)>,
    /// Machine-readable reason for the choice.
    pub chosen_because: String,
    /// True when the pick came from budget fallback, not class preference.
    pub fallback_used: bool,
    /// Every candidate scanned, in menu order, with verdicts.
    pub considered: Vec<ConsideredEntry>,
    /// Whether the chosen route renders audio (oracle routes do not).
    pub render_audio: bool,
    /// Mid-note hop disposition.
    pub hop: HopVerdict,
}

/// Read one `"key":<value>` field out of a canonical receipt line. String
/// values are unescaped; bare literals (`null`, numbers) come back raw.
fn receipt_field(line: &str, key: &str) -> Result<String, SelectorError> {
    let pat = format!("\"{key}\":");
    let idx = line
        .find(&pat)
        .ok_or_else(|| SelectorError::Receipt(format!("missing key {key:?}")))?;
    let tail = line[idx + pat.len()..].trim_start();
    if let Some(body) = tail.strip_prefix('"') {
        let end = body
            .find('"')
            .ok_or_else(|| SelectorError::Receipt(format!("unterminated {key:?}")))?;
        Ok(body[..end].replace("\\\"", "\"").replace("\\\\", "\\"))
    } else {
        let end = tail.find([',', '}']).unwrap_or(tail.len());
        Ok(tail[..end].trim().to_string())
    }
}

/// Decode the `considered:[...]` payload of a canonical receipt line.
fn decode_considered(blob: &str) -> Result<Vec<ConsideredEntry>, SelectorError> {
    let mut considered = Vec::new();
    for item in blob.split("},{") {
        if item.is_empty() {
            continue;
        }
        let item = item.trim_matches(|c| c == '{' || c == '}');
        let field = |k: &str| -> Result<String, SelectorError> {
            let pat = format!("\"{k}\":");
            let i = item
                .find(&pat)
                .ok_or_else(|| SelectorError::Receipt(format!("considered missing {k:?}")))?;
            let raw = &item[i + pat.len()..];
            if let Some(r) = raw.strip_prefix('"') {
                Ok(r.split('"').next().unwrap_or_default().to_string())
            } else {
                let end = raw.find([',', '}']).unwrap_or(raw.len());
                Ok(raw[..end].trim().to_string())
            }
        };
        considered.push(ConsideredEntry {
            image: field("image")?,
            qoi: field("qoi")?,
            gate: field("gate")?,
            cost_tier: field("cost_tier")?
                .parse::<u8>()
                .map_err(|_| SelectorError::Receipt("bad cost_tier".into()))?,
            verdict: field("verdict")?,
        });
    }
    Ok(considered)
}

impl SelectionReceipt {
    /// Canonical JSONL encoding (single line, fixed field order).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let esc = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"");
        let considered = self
            .considered
            .iter()
            .map(|c| {
                format!(
                    "{{\"image\":\"{}\",\"qoi\":\"{}\",\"gate\":\"{}\",\"cost_tier\":{},\"verdict\":\"{}\"}}",
                    esc(&c.image),
                    esc(&c.qoi),
                    esc(&c.gate),
                    c.cost_tier,
                    esc(&c.verdict)
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let (img, qoi) = match &self.chosen {
            Some((i, q)) => (format!("\"{}\"", esc(i)), format!("\"{}\"", esc(q))),
            None => ("null".to_string(), "null".to_string()),
        };
        format!(
            "{{\"schema\":\"{}\",\"filling\":\"{}\",\"serving\":\"{}\",\"phase\":\"{}\",\"budget_headroom_tier\":{},\"chosen_image\":{},\"chosen_qoi\":{},\"chosen_because\":\"{}\",\"fallback_used\":{},\"render_audio\":{},\"hop\":\"{}\",\"considered\":[{}]}}",
            esc(&self.schema),
            esc(&self.filling),
            self.serving.as_str(),
            self.phase.as_str(),
            self.budget_headroom_tier,
            img,
            qoi,
            esc(&self.chosen_because),
            self.fallback_used,
            self.render_audio,
            self.hop.as_str(),
            considered
        )
    }

    /// Decode [`Self::to_jsonl`] output. Strict on the consumed fields.
    ///
    /// # Errors
    /// [`SelectorError::Receipt`] on any structural mismatch.
    pub fn from_jsonl(line: &str) -> Result<Self, SelectorError> {
        let get = |key: &str| receipt_field(line, key);
        let parse_bool = |key: &str| -> Result<bool, SelectorError> {
            match get(key)?.as_str() {
                "true" => Ok(true),
                "false" => Ok(false),
                other => Err(SelectorError::Receipt(format!(
                    "bad bool {key:?}={other:?}"
                ))),
            }
        };
        let schema = get("schema")?;
        if schema != SELECTION_RECEIPT_SCHEMA {
            return Err(SelectorError::Receipt(format!(
                "schema {schema:?} is not {SELECTION_RECEIPT_SCHEMA:?}"
            )));
        }
        let chosen_image_raw = get("chosen_image")?;
        let chosen_qoi_raw = get("chosen_qoi")?;
        let chosen = if chosen_image_raw == "null" {
            if chosen_qoi_raw != "null" {
                return Err(SelectorError::Receipt("chosen qoi without image".into()));
            }
            None
        } else {
            Some((chosen_image_raw, chosen_qoi_raw))
        };
        let considered_blob = {
            let pat = "\"considered\":[";
            let idx = line
                .rfind(pat)
                .ok_or_else(|| SelectorError::Receipt("missing considered".into()))?;
            let tail = &line[idx + pat.len()..];
            let end = tail
                .rfind(']')
                .ok_or_else(|| SelectorError::Receipt("unterminated considered".into()))?;
            tail[..end].to_string()
        };
        Ok(Self {
            schema,
            filling: get("filling")?,
            serving: Serving::parse(&get("serving")?)?,
            phase: GesturePhase::parse(&get("phase")?)?,
            budget_headroom_tier: get("budget_headroom_tier")?
                .parse::<u8>()
                .map_err(|_| SelectorError::Receipt("bad number".into()))?,
            chosen,
            chosen_because: get("chosen_because")?,
            fallback_used: parse_bool("fallback_used")?,
            render_audio: parse_bool("render_audio")?,
            hop: HopVerdict::parse(&get("hop")?)?,
            considered: decode_considered(&considered_blob)?,
        })
    }

    /// Content identity of the canonical JSONL bytes (provenance citation).
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        hash_domain(RECEIPT_HASH_DOMAIN, self.to_jsonl().as_bytes())
    }
}

/// Typed selection failures. A typed refusal beats a fabricated pick.
#[derive(Clone, Debug, PartialEq)]
pub enum SelectorError {
    /// Registry parse/validation failure.
    Registry(String),
    /// Hop-policy artifact parse failure.
    Policy(String),
    /// Receipt encode/decode failure.
    Receipt(String),
    /// The menu named a filling the registry snapshot does not know at all.
    UnknownFilling(String),
    /// No candidate passed admission; the message lists each verdict.
    NoAdmittedImage(String),
}

impl core::fmt::Display for SelectorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Registry(d) => write!(f, "registry: {d}"),
            Self::Policy(d) => write!(f, "hop policy: {d}"),
            Self::Receipt(d) => write!(f, "receipt: {d}"),
            Self::UnknownFilling(filling) => {
                write!(f, "registry knows no rows for filling {filling:?}")
            }
            Self::NoAdmittedImage(d) => write!(f, "no admitted image: {d}"),
        }
    }
}

impl std::error::Error for SelectorError {}

/// Session posture: does THIS render claim the registry's live-default slot?
///
/// D25 bites here: live-default sessions enforce the full gate+budget-row
/// conjunction; explicitly non-default sessions (experiments, bake-offs,
/// probes) may select ungated-but-admitted rows because they are not being
/// offered to users as THE default.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SessionPosture {
    /// User-facing default: D25 fully enforced.
    LiveDefault,
    /// Explicitly non-default: ungated rows become selectable, still never
    /// `refused` ones.
    NonDefaultDeclared,
}

/// One selection request at one control boundary.
#[derive(Clone, Copy, Debug)]
pub struct SelectionRequest {
    /// Claim served this block.
    pub serving: Serving,
    /// Detector-derived gesture phase.
    pub phase: GesturePhase,
    /// Cheapest affordable tier: candidates with
    /// `cost_tier > headroom` are unaffordable this block.
    pub budget_headroom_tier: u8,
    /// Image currently loaded, when one is (hop reasoning needs it).
    pub current_image: Option<&'static str>,
    /// Is the current control point a phrase boundary (hop permitted)?
    pub at_phrase_boundary: bool,
    /// Session posture (D25 enforcement level).
    pub posture: SessionPosture,
}

/// Scan a filling's menu against the registry, session posture, and budget,
/// recording every candidate's verdict in `considered` (menu order).
/// Returns the admitted `(menu index, entry)` pairs.
fn admit_candidates<'a>(
    registry: &RegistrySnapshot,
    menu: &'a VoiceMenu,
    request: &SelectionRequest,
    considered: &mut Vec<ConsideredEntry>,
) -> Vec<(usize, &'a MenuEntry)> {
    let mut admitted = Vec::new();
    for (idx, entry) in menu.entries.iter().enumerate() {
        // Serving discipline: the oracle ANSWERS questions and never
        // renders audio; audio images never answer peak questions.
        let serving_ok = match entry.class {
            EntryClass::SpectralOracle => request.serving == Serving::PeakLocations,
            _ => request.serving == Serving::Play,
        };
        if !serving_ok {
            let verdict = if request.serving == Serving::Play {
                "oracle-not-for-audio"
            } else {
                "audio-image-not-oracle-route"
            };
            considered.push(ConsideredEntry {
                image: entry.image.to_string(),
                qoi: entry.qoi.to_string(),
                gate: "n/a".to_string(),
                cost_tier: entry.cost_tier,
                verdict: verdict.to_string(),
            });
            continue;
        }
        let Some(row) = registry.row(menu.filling, entry.image, entry.qoi) else {
            considered.push(ConsideredEntry {
                image: entry.image.to_string(),
                qoi: entry.qoi.to_string(),
                gate: "absent".to_string(),
                cost_tier: entry.cost_tier,
                verdict: "no-registry-row".to_string(),
            });
            continue;
        };
        let gate_s = row.gate.as_str();
        let base = |verdict: &str| ConsideredEntry {
            image: entry.image.to_string(),
            qoi: entry.qoi.to_string(),
            gate: gate_s.to_string(),
            cost_tier: entry.cost_tier,
            verdict: verdict.to_string(),
        };
        match row.gate {
            Gate::Refused => considered.push(base("refused-image")),
            Gate::Ungated => {
                if request.posture == SessionPosture::LiveDefault {
                    considered.push(base("ungated-not-live-default"));
                } else if entry.cost_tier > request.budget_headroom_tier {
                    considered.push(base("over-budget"));
                } else {
                    considered.push(base("admitted-non-default"));
                    admitted.push((idx, entry));
                }
            }
            Gate::Green => {
                // D25 full form binds only the live-default slot.
                if request.posture == SessionPosture::LiveDefault
                    && row.live_default
                    && row.budget_row.is_none()
                {
                    considered.push(base("live-default-without-budget-row"));
                } else if entry.cost_tier > request.budget_headroom_tier {
                    considered.push(base("over-budget"));
                } else {
                    considered.push(base("admitted"));
                    admitted.push((idx, entry));
                }
            }
        }
    }
    admitted
}

/// The selector. Stateless: all inputs are explicit, all outputs receipts.
#[derive(Debug, Clone, Copy, Default)]
pub struct Selector;

impl Selector {
    /// Decide the image for one control boundary and emit its receipt.
    ///
    /// # Errors
    /// [`SelectorError::UnknownFilling`] when the registry has no rows for
    /// the filling at all; [`SelectorError::NoAdmittedImage`] when every
    /// candidate refused (there is no honest pick to cite).
    pub fn select(
        registry: &RegistrySnapshot,
        menu: &VoiceMenu,
        request: &SelectionRequest,
    ) -> Result<SelectionReceipt, SelectorError> {
        if registry.rows.iter().all(|r| r.filling != menu.filling) {
            return Err(SelectorError::UnknownFilling(menu.filling.to_string()));
        }
        let mut considered = Vec::new();
        let admitted = admit_candidates(registry, menu, request, &mut considered);
        if admitted.is_empty() {
            let detail = considered
                .iter()
                .map(|c| format!("{}={}", c.image, c.verdict))
                .collect::<Vec<_>>()
                .join(", ");
            return Err(SelectorError::NoAdmittedImage(detail));
        }
        // Preference scan: class preference first (gesture/serving policy),
        // then cost. Menu order breaks ties — deterministic.
        let preferred_class = match request.serving {
            Serving::PeakLocations => EntryClass::SpectralOracle,
            Serving::Play => match request.phase {
                GesturePhase::Phrase => EntryClass::SpatialOrIsland,
                GesturePhase::Held => EntryClass::SettledHold,
            },
        };
        let class_pick = admitted
            .iter()
            .filter(|(_, e)| e.class == preferred_class)
            .min_by_key(|(idx, e)| (e.cost_tier, *idx))
            .copied();
        let cheap_pick = admitted
            .iter()
            .min_by_key(|(idx, e)| (e.cost_tier, *idx))
            .copied();
        let fallback_used = class_pick.is_none();
        // `admitted` was checked non-empty above, so `cheap_pick` is always
        // Some and this last resort is dead in practice — kept total so the
        // decision core contains no panic paths.
        let (_menu_idx, pick_entry) = class_pick.or(cheap_pick).unwrap_or(admitted[0]);
        let render_audio = request.serving != Serving::PeakLocations;
        let hop = if request.current_image == Some(pick_entry.image) {
            HopVerdict::NoChange
        } else if request.at_phrase_boundary {
            HopVerdict::AtBoundary
        } else if pick_entry.lift_supported {
            // The measured policy already gated the PHASE (a Held-phase
            // pick IS the settled state its parameters define), so a
            // mid-note switch is admitted exactly when the target exposes
            // the D17 lift.
            HopVerdict::MidNoteAllowed
        } else {
            HopVerdict::DeferredToBoundary
        };
        let chosen_because = if request.serving == Serving::PeakLocations {
            "spectral-question-routed-to-oracle-no-audio"
        } else if fallback_used {
            "budget-squeeze-to-cheapest-admitted"
        } else {
            match request.phase {
                GesturePhase::Phrase => "phrase-prefers-spatial-or-island",
                GesturePhase::Held => "held-prefers-gated-hold",
            }
        };
        Ok(SelectionReceipt {
            schema: SELECTION_RECEIPT_SCHEMA.to_string(),
            filling: menu.filling.to_string(),
            serving: request.serving,
            phase: request.phase,
            budget_headroom_tier: request.budget_headroom_tier,
            chosen: Some((pick_entry.image.to_string(), pick_entry.qoi.to_string())),
            chosen_because: chosen_because.to_string(),
            fallback_used,
            considered,
            render_audio,
            hop,
        })
    }
}
