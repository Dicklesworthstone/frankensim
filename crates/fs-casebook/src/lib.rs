//! The shared conformance-case harness (bead huq.5, plan §13.3): named
//! cases with structured JSON-lines records — the executable half of the
//! CONTRACT.md discipline.
//!
//! Every crate's conformance suite registers named cases; running the suite
//! executes them in registration order. Legacy cases retain their stable v1
//! JSON line carrying the case id, input digest, tolerance, verdict, and
//! evidence pointers. Replayable cases add a companion with the exact rerun
//! command and canonical input frame; implementation/reference differences add
//! a first-disagreement companion suitable for a bug report.
//!
//! The case format is DATA-FIRST so an fs-ir front end can wrap it additively:
//! [`Case`] execution returns a [`CaseOutcome`] value and [`SuiteReport`] holds
//! typed records. Printing is a thin, separable layer, and an IR-speaking
//! runner can consume the same records without rewriting suites.
//!
//! No-claims: the harness runs and records; it does not decide what a
//! tolerance MEANS physically, does not compare across ISAs by itself
//! (cross-ISA evidence is a Gauntlet lane running the same suite on both
//! hosts), and awards no certification tiers (that is fs-conform's
//! converter scope).

use core::fmt::{self, Write as _};

/// Version stamped into every emitted record.
pub const CASEBOOK_RECORD_VERSION: u32 = 1;

/// Version stamped into every replay companion record.
pub const CASEBOOK_REPLAY_RECORD_VERSION: u32 = 1;

/// Version stamped into every structured disagreement record.
pub const CASEBOOK_DISAGREEMENT_RECORD_VERSION: u32 = 1;

/// FNV-1a 64 over raw bytes — the canonical inputs-digest helper so
/// suites never hand-roll their own.
#[must_use]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn hex_encode(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut out, "{byte:02x}").expect("writing to String cannot fail");
    }
    out
}

fn hex_nibble(byte: u8, offset: usize) -> Result<u8, ReplayError> {
    match byte {
        b'0'..=b'9' => Ok(byte - b'0'),
        b'a'..=b'f' => Ok(byte - b'a' + 10),
        b'A'..=b'F' => Ok(byte - b'A' + 10),
        _ => Err(ReplayError::InvalidHex { offset, byte }),
    }
}

fn hex_decode(encoded: &str) -> Result<Vec<u8>, ReplayError> {
    if encoded.len() % 2 != 0 {
        return Err(ReplayError::OddHexLength {
            length: encoded.len(),
        });
    }
    let encoded = encoded.as_bytes();
    let mut out = Vec::with_capacity(encoded.len() / 2);
    for offset in (0..encoded.len()).step_by(2) {
        let high = hex_nibble(encoded[offset], offset)?;
        let low = hex_nibble(encoded[offset + 1], offset + 1)?;
        out.push((high << 4) | low);
    }
    Ok(out)
}

/// A replay companion record failed canonical-frame verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReplayError {
    /// Hex text must contain exactly two characters per byte.
    OddHexLength {
        /// Observed character count.
        length: usize,
    },
    /// A non-hex character occurred in the canonical frame.
    InvalidHex {
        /// Character offset in the encoded frame.
        offset: usize,
        /// Invalid byte.
        byte: u8,
    },
    /// The frame decoded, but its text was not canonical lowercase hex.
    NonCanonicalHex {
        /// Hex text carried by the replay record.
        declared: String,
        /// Canonical lowercase encoding of the same bytes.
        canonical: String,
    },
    /// The decoded frame length did not match the declared length.
    LengthMismatch {
        /// Length carried by the replay record.
        declared: usize,
        /// Length reconstructed from its hex frame.
        reconstructed: usize,
    },
    /// The reconstructed frame did not hash to the declared digest.
    DigestMismatch {
        /// Digest carried by the replay record.
        declared: String,
        /// Digest recomputed from the reconstructed bytes.
        reconstructed: String,
    },
}

impl fmt::Display for ReplayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OddHexLength { length } => {
                write!(f, "canonical input hex has odd length {length}")
            }
            Self::InvalidHex { offset, byte } => write!(
                f,
                "canonical input hex has invalid byte 0x{byte:02x} at offset {offset}"
            ),
            Self::NonCanonicalHex {
                declared,
                canonical,
            } => write!(
                f,
                "canonical input hex is not lowercase: declared {declared}, canonical {canonical}"
            ),
            Self::LengthMismatch {
                declared,
                reconstructed,
            } => write!(
                f,
                "canonical input length mismatch: declared {declared}, reconstructed {reconstructed}"
            ),
            Self::DigestMismatch {
                declared,
                reconstructed,
            } => write!(
                f,
                "canonical input digest mismatch: declared {declared}, reconstructed {reconstructed}"
            ),
        }
    }
}

impl std::error::Error for ReplayError {}

/// Exact instructions and canonical input bytes for replaying one case.
///
/// [`Suite::case_replayable`] derives the case digest from these bytes;
/// callers cannot supply a digest that silently disagrees with the retained
/// frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplaySpec {
    command: String,
    canonical_inputs: Vec<u8>,
}

impl ReplaySpec {
    /// Construct a replay specification from exact canonical input bytes.
    ///
    /// # Panics
    /// Panics when `command` is empty because such a record would not teach an
    /// operator how to rerun the case.
    #[must_use]
    pub fn new(command: impl Into<String>, canonical_inputs: impl Into<Vec<u8>>) -> Self {
        let command = command.into();
        assert!(!command.is_empty(), "replay command must be non-empty");
        Self {
            command,
            canonical_inputs: canonical_inputs.into(),
        }
    }

    /// Reconstruct a replay specification from a retained canonical hex frame.
    ///
    /// # Errors
    /// [`ReplayError`] if `canonical_inputs_hex` is malformed.
    ///
    /// # Panics
    /// Panics when `command` is empty; see [`ReplaySpec::new`].
    pub fn from_hex(
        command: impl Into<String>,
        canonical_inputs_hex: &str,
    ) -> Result<Self, ReplayError> {
        Ok(Self::new(command, hex_decode(canonical_inputs_hex)?))
    }

    /// Exact command or selector that reruns the owning test case.
    #[must_use]
    pub fn command(&self) -> &str {
        &self.command
    }

    /// Canonical case-input bytes.
    #[must_use]
    pub fn canonical_inputs(&self) -> &[u8] {
        &self.canonical_inputs
    }

    /// Canonical input frame encoded as lowercase hex.
    #[must_use]
    pub fn canonical_inputs_hex(&self) -> String {
        hex_encode(&self.canonical_inputs)
    }

    /// FNV-1a digest derived from the canonical input bytes.
    #[must_use]
    pub fn inputs_digest(&self) -> u64 {
        fnv1a64(&self.canonical_inputs)
    }

    fn record(&self, suite: &str, case: &str) -> ReplayRecord {
        ReplayRecord {
            version: CASEBOOK_REPLAY_RECORD_VERSION,
            suite: suite.to_owned(),
            case: case.to_owned(),
            command: self.command.clone(),
            inputs_digest: format!("{:016x}", self.inputs_digest()),
            inputs_len: self.canonical_inputs.len(),
            canonical_inputs_hex: self.canonical_inputs_hex(),
        }
    }
}

/// Bug-report-ready replay metadata emitted beside a replayable case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplayRecord {
    /// Replay-record schema version.
    pub version: u32,
    /// Owning suite.
    pub suite: String,
    /// Stable case id.
    pub case: String,
    /// Exact command or selector that reruns the case.
    pub command: String,
    /// FNV-1a digest of the canonical input frame.
    pub inputs_digest: String,
    /// Canonical input length in bytes.
    pub inputs_len: usize,
    /// Complete canonical input frame as lowercase hex.
    pub canonical_inputs_hex: String,
}

impl ReplayRecord {
    /// Decode the retained frame and verify its declared length and digest.
    ///
    /// # Errors
    /// [`ReplayError`] if the hex is malformed or either declaration differs
    /// from the reconstructed bytes.
    pub fn verify_and_decode(&self) -> Result<Vec<u8>, ReplayError> {
        let decoded = hex_decode(&self.canonical_inputs_hex)?;
        let canonical = hex_encode(&decoded);
        if canonical != self.canonical_inputs_hex {
            return Err(ReplayError::NonCanonicalHex {
                declared: self.canonical_inputs_hex.clone(),
                canonical,
            });
        }
        if decoded.len() != self.inputs_len {
            return Err(ReplayError::LengthMismatch {
                declared: self.inputs_len,
                reconstructed: decoded.len(),
            });
        }
        let reconstructed = format!("{:016x}", fnv1a64(&decoded));
        if reconstructed != self.inputs_digest {
            return Err(ReplayError::DigestMismatch {
                declared: self.inputs_digest.clone(),
                reconstructed,
            });
        }
        Ok(decoded)
    }
}

/// The first exact-byte disagreement between an implementation and reference.
///
/// A byte mismatch carries both bytes. A length-boundary mismatch carries
/// `None` for the side that ended first. Digests bind the complete frames even
/// though the record localizes only the first disagreement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DisagreementRecord {
    /// Disagreement-record schema version.
    pub version: u32,
    /// Owning suite.
    pub suite: String,
    /// Stable case id.
    pub case: String,
    /// Implementation under test.
    pub implementation: String,
    /// Comparison reference or oracle.
    pub reference: String,
    /// Complete implementation-frame length.
    pub implementation_len: usize,
    /// Complete reference-frame length.
    pub reference_len: usize,
    /// FNV-1a digest of the implementation frame.
    pub implementation_digest: String,
    /// FNV-1a digest of the reference frame.
    pub reference_digest: String,
    /// First byte offset where the frames differ.
    pub first_offset: usize,
    /// Implementation byte at `first_offset`, or `None` at its length boundary.
    pub implementation_byte: Option<u8>,
    /// Reference byte at `first_offset`, or `None` at its length boundary.
    pub reference_byte: Option<u8>,
}

impl DisagreementRecord {
    /// Compare two exact frames and localize their first disagreement.
    /// Returns `None` for identical frames.
    ///
    /// `suite` and `case` provide context when the record is rendered on its
    /// own. [`Suite::run`] always binds attached records to their actual owning
    /// registration before emission and report lookup.
    #[must_use]
    pub fn first(
        suite: impl Into<String>,
        case: impl Into<String>,
        implementation: impl Into<String>,
        reference: impl Into<String>,
        implementation_frame: &[u8],
        reference_frame: &[u8],
    ) -> Option<Self> {
        let shared = implementation_frame.len().min(reference_frame.len());
        let first_offset = (0..shared)
            .find(|&offset| implementation_frame[offset] != reference_frame[offset])
            .or_else(|| (implementation_frame.len() != reference_frame.len()).then_some(shared))?;
        Some(Self {
            version: CASEBOOK_DISAGREEMENT_RECORD_VERSION,
            suite: suite.into(),
            case: case.into(),
            implementation: implementation.into(),
            reference: reference.into(),
            implementation_len: implementation_frame.len(),
            reference_len: reference_frame.len(),
            implementation_digest: format!("{:016x}", fnv1a64(implementation_frame)),
            reference_digest: format!("{:016x}", fnv1a64(reference_frame)),
            first_offset,
            implementation_byte: implementation_frame.get(first_offset).copied(),
            reference_byte: reference_frame.get(first_offset).copied(),
        })
    }

    /// Stable mismatch classification used in the JSON record.
    #[must_use]
    pub fn mismatch_kind(&self) -> &'static str {
        if self.implementation_byte.is_some() && self.reference_byte.is_some() {
            "byte"
        } else {
            "length-boundary"
        }
    }
}

/// The tolerance a case claims to hold. `Structural` marks cases whose
/// verdict is a typed structural fact (a refusal fired, a schema
/// round-tripped) rather than a numeric bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ToleranceSpec {
    /// Bit-exact equality.
    Exact,
    /// At most this many ULPs of deviation.
    Ulps(u32),
    /// Relative deviation at most this bound.
    RelativeLe(f64),
    /// Absolute deviation at most this bound.
    AbsoluteLe(f64),
    /// A typed structural verdict with no numeric bound.
    Structural,
}

impl fmt::Display for ToleranceSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Exact => write!(f, "exact"),
            Self::Ulps(n) => write!(f, "ulps<={n}"),
            Self::RelativeLe(b) => write!(f, "rel<={b:e}"),
            Self::AbsoluteLe(b) => write!(f, "abs<={b:e}"),
            Self::Structural => write!(f, "structural"),
        }
    }
}

/// What one executed case reports back.
#[derive(Debug, Clone)]
pub struct CaseOutcome {
    /// Whether the case held its claim.
    pub pass: bool,
    /// Human-readable measurement detail (goes into the record verbatim).
    pub details: String,
    /// Evidence pointers: fixture hashes, artifact paths, log ids.
    pub evidence: Vec<String>,
    /// Structured implementation/reference disagreements discovered by the case.
    pub disagreements: Vec<DisagreementRecord>,
}

impl CaseOutcome {
    /// A passing outcome.
    #[must_use]
    pub fn pass(details: impl Into<String>) -> CaseOutcome {
        CaseOutcome {
            pass: true,
            details: details.into(),
            evidence: Vec::new(),
            disagreements: Vec::new(),
        }
    }

    /// A failing outcome.
    #[must_use]
    pub fn fail(details: impl Into<String>) -> CaseOutcome {
        CaseOutcome {
            pass: false,
            details: details.into(),
            evidence: Vec::new(),
            disagreements: Vec::new(),
        }
    }

    /// Attach an evidence pointer.
    #[must_use]
    pub fn with_evidence(mut self, pointer: impl Into<String>) -> CaseOutcome {
        self.evidence.push(pointer.into());
        self
    }

    /// Attach a localized disagreement and fail closed.
    ///
    /// Setting `pass = false` here prevents a caller from recording a known
    /// implementation/reference difference while accidentally leaving the
    /// owning case green. [`Suite::run`] also binds the record's suite and case
    /// fields to the registration that returned this outcome.
    #[must_use]
    pub fn with_disagreement(mut self, disagreement: DisagreementRecord) -> CaseOutcome {
        self.pass = false;
        self.disagreements.push(disagreement);
        self
    }
}

/// One executed case's structured record — the unit of the JSON-lines
/// output and of cross-run comparison.
#[derive(Debug, Clone)]
pub struct CaseRecord {
    /// Record schema version.
    pub version: u32,
    /// The owning suite name.
    pub suite: String,
    /// The stable case id.
    pub case: String,
    /// Hex FNV-1a 64 digest of the case's declared inputs.
    pub inputs_digest: String,
    /// The claimed tolerance, rendered stably.
    pub tolerance: String,
    /// The verdict.
    pub pass: bool,
    /// Measurement detail.
    pub details: String,
    /// Evidence pointers.
    pub evidence: Vec<String>,
}

fn escape_json(s: &str, out: &mut String) {
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                write!(out, "\\u{:04x}", c as u32)
                    .expect("writing JSON escape bytes to a String cannot fail");
            }
            c => out.push(c),
        }
    }
}

impl ReplayRecord {
    /// Render one deterministic JSONL replay companion.
    #[must_use]
    pub fn json_line(&self) -> String {
        let mut out = String::with_capacity(self.canonical_inputs_hex.len() + 192);
        out.push_str("{\"casebook_replay\":");
        out.push_str(&self.version.to_string());
        out.push_str(",\"suite\":\"");
        escape_json(&self.suite, &mut out);
        out.push_str("\",\"case\":\"");
        escape_json(&self.case, &mut out);
        out.push_str("\",\"command\":\"");
        escape_json(&self.command, &mut out);
        out.push_str("\",\"inputs_digest\":\"");
        escape_json(&self.inputs_digest, &mut out);
        out.push_str("\",\"inputs_len\":");
        out.push_str(&self.inputs_len.to_string());
        out.push_str(",\"canonical_inputs_hex\":\"");
        escape_json(&self.canonical_inputs_hex, &mut out);
        out.push_str("\"}");
        out
    }
}

fn push_optional_byte_json(out: &mut String, byte: Option<u8>) {
    match byte {
        Some(byte) => {
            out.push('"');
            write!(out, "{byte:02x}").expect("writing to String cannot fail");
            out.push('"');
        }
        None => out.push_str("null"),
    }
}

impl DisagreementRecord {
    /// Render one deterministic, bug-report-ready JSONL disagreement.
    #[must_use]
    pub fn json_line(&self) -> String {
        let mut out = String::with_capacity(320);
        out.push_str("{\"casebook_disagreement\":");
        out.push_str(&self.version.to_string());
        out.push_str(",\"suite\":\"");
        escape_json(&self.suite, &mut out);
        out.push_str("\",\"case\":\"");
        escape_json(&self.case, &mut out);
        out.push_str("\",\"implementation\":\"");
        escape_json(&self.implementation, &mut out);
        out.push_str("\",\"reference\":\"");
        escape_json(&self.reference, &mut out);
        out.push_str("\",\"mismatch_kind\":\"");
        out.push_str(self.mismatch_kind());
        out.push_str("\",\"implementation_len\":");
        out.push_str(&self.implementation_len.to_string());
        out.push_str(",\"reference_len\":");
        out.push_str(&self.reference_len.to_string());
        out.push_str(",\"implementation_digest\":\"");
        escape_json(&self.implementation_digest, &mut out);
        out.push_str("\",\"reference_digest\":\"");
        escape_json(&self.reference_digest, &mut out);
        out.push_str("\",\"first_offset\":");
        out.push_str(&self.first_offset.to_string());
        out.push_str(",\"implementation_byte\":");
        push_optional_byte_json(&mut out, self.implementation_byte);
        out.push_str(",\"reference_byte\":");
        push_optional_byte_json(&mut out, self.reference_byte);
        out.push('}');
        out
    }
}

impl CaseRecord {
    /// The record as one JSON line (hand-rolled, dependency-free,
    /// deterministic field order).
    #[must_use]
    pub fn json_line(&self) -> String {
        let mut out = String::with_capacity(128);
        out.push_str("{\"casebook\":");
        out.push_str(&self.version.to_string());
        out.push_str(",\"suite\":\"");
        escape_json(&self.suite, &mut out);
        out.push_str("\",\"case\":\"");
        escape_json(&self.case, &mut out);
        out.push_str("\",\"inputs_digest\":\"");
        escape_json(&self.inputs_digest, &mut out);
        out.push_str("\",\"tolerance\":\"");
        escape_json(&self.tolerance, &mut out);
        out.push_str("\",\"pass\":");
        out.push_str(if self.pass { "true" } else { "false" });
        out.push_str(",\"details\":\"");
        escape_json(&self.details, &mut out);
        out.push_str("\",\"evidence\":[");
        for (i, e) in self.evidence.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push('"');
            escape_json(e, &mut out);
            out.push('"');
        }
        out.push_str("]}");
        out
    }
}

type CaseFn = Box<dyn FnOnce() -> CaseOutcome>;

struct Case {
    id: &'static str,
    inputs_digest: u64,
    tolerance: ToleranceSpec,
    replay: Option<ReplaySpec>,
    run: CaseFn,
}

/// A named conformance suite: cases registered in a fixed order, executed
/// deterministically in that order.
pub struct Suite {
    name: &'static str,
    cases: Vec<Case>,
}

impl Suite {
    /// A new suite. The name is the stable cross-run identity prefix.
    ///
    /// # Panics
    /// If the name is empty.
    #[must_use]
    pub fn new(name: &'static str) -> Suite {
        assert!(!name.is_empty(), "suite name must be non-empty");
        Suite {
            name,
            cases: Vec::new(),
        }
    }

    /// Register one named case. `inputs_digest` binds the exact inputs
    /// the case runs on (use [`fnv1a64`] over the canonical input bytes);
    /// `tolerance` is the claim the case holds.
    #[must_use]
    pub fn case(
        mut self,
        id: &'static str,
        inputs_digest: u64,
        tolerance: ToleranceSpec,
        run: impl FnOnce() -> CaseOutcome + 'static,
    ) -> Suite {
        self.cases.push(Case {
            id,
            inputs_digest,
            tolerance,
            replay: None,
            run: Box::new(run),
        });
        self
    }

    /// Register a case with a complete replay companion.
    ///
    /// Unlike [`Suite::case`], this derives the input digest from the exact
    /// canonical bytes in `replay`, eliminating caller-supplied digest/frame
    /// drift. The legacy case record remains byte-for-byte unchanged; a
    /// separate replay JSONL record carries the command and full input frame.
    #[must_use]
    pub fn case_replayable(
        mut self,
        id: &'static str,
        replay: ReplaySpec,
        tolerance: ToleranceSpec,
        run: impl FnOnce() -> CaseOutcome + 'static,
    ) -> Suite {
        let inputs_digest = replay.inputs_digest();
        self.cases.push(Case {
            id,
            inputs_digest,
            tolerance,
            replay: Some(replay),
            run: Box::new(run),
        });
        self
    }

    /// Execute every case in registration order, emitting one JSON line
    /// per case to stdout, and return the typed report. Duplicate case
    /// ids and empty ids are recorded as structural FAILURES (fail
    /// closed), never silently accepted.
    #[must_use]
    pub fn run(self) -> SuiteReport {
        let mut seen: Vec<&'static str> = Vec::new();
        let mut records = Vec::with_capacity(self.cases.len());
        let mut replay_records = Vec::new();
        let mut disagreements = Vec::new();
        for case in self.cases {
            let Case {
                id,
                inputs_digest,
                tolerance,
                replay,
                run,
            } = case;
            let structural_defect = if id.is_empty() {
                Some("empty case id".to_owned())
            } else if seen.contains(&id) {
                Some(format!("duplicate case id {id:?}"))
            } else {
                seen.push(id);
                None
            };
            let outcome = match structural_defect {
                Some(defect) => CaseOutcome::fail(defect),
                None => run(),
            };
            let CaseOutcome {
                pass,
                details,
                evidence,
                disagreements: mut case_disagreements,
            } = outcome;
            let pass = pass && case_disagreements.is_empty();
            for disagreement in &mut case_disagreements {
                disagreement.suite = self.name.to_owned();
                disagreement.case = id.to_owned();
            }
            let record = CaseRecord {
                version: CASEBOOK_RECORD_VERSION,
                suite: self.name.to_owned(),
                case: id.to_owned(),
                inputs_digest: format!("{inputs_digest:016x}"),
                tolerance: tolerance.to_string(),
                pass,
                details,
                evidence,
            };
            println!("{}", record.json_line());
            if let Some(replay) = replay {
                let replay_record = replay.record(self.name, id);
                println!("{}", replay_record.json_line());
                replay_records.push(replay_record);
            }
            for disagreement in case_disagreements {
                println!("{}", disagreement.json_line());
                disagreements.push(disagreement);
            }
            records.push(record);
        }
        SuiteReport {
            records,
            replay_records,
            disagreements,
        }
    }
}

/// The typed result of one suite run.
#[derive(Debug, Clone)]
pub struct SuiteReport {
    /// One record per executed case, in execution order.
    pub records: Vec<CaseRecord>,
    /// Replay companions, in owning-case registration order.
    pub replay_records: Vec<ReplayRecord>,
    /// Structured disagreements, in owning-case and discovery order.
    pub disagreements: Vec<DisagreementRecord>,
}

impl SuiteReport {
    /// Whether every case passed (an empty suite is NOT green — a suite
    /// that ran nothing proved nothing).
    #[must_use]
    pub fn all_passed(&self) -> bool {
        !self.records.is_empty() && self.records.iter().all(|r| r.pass)
    }

    /// The failing records.
    #[must_use]
    pub fn failures(&self) -> Vec<&CaseRecord> {
        self.records.iter().filter(|r| !r.pass).collect()
    }

    /// Replay companion for one case, if it registered one.
    #[must_use]
    pub fn replay_for(&self, case: &str) -> Option<&ReplayRecord> {
        self.replay_records
            .iter()
            .find(|record| record.case == case)
    }

    /// Structured disagreements attached to one case.
    #[must_use]
    pub fn disagreements_for(&self, case: &str) -> Vec<&DisagreementRecord> {
        self.disagreements
            .iter()
            .filter(|record| record.case == case)
            .collect()
    }

    /// Assert the suite is green, panicking with every failing record's
    /// JSON line so the test log carries the full structured evidence.
    ///
    /// # Panics
    /// If any case failed or the suite ran zero cases.
    pub fn assert_green(&self) {
        if self.all_passed() {
            return;
        }
        let mut msg = String::from("conformance suite not green:\n");
        if self.records.is_empty() {
            msg.push_str("  (zero cases executed)\n");
        }
        for failure in self.failures() {
            msg.push_str("  ");
            msg.push_str(&failure.json_line());
            msg.push('\n');
            if let Some(replay) = self.replay_for(&failure.case) {
                msg.push_str("  ");
                msg.push_str(&replay.json_line());
                msg.push('\n');
            }
            for disagreement in self.disagreements_for(&failure.case) {
                msg.push_str("  ");
                msg.push_str(&disagreement.json_line());
                msg.push('\n');
            }
        }
        panic!("{msg}");
    }
}
