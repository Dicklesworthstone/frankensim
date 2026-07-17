//! The shared conformance-case harness (bead huq.5, plan §13.3): named
//! cases with structured JSON-lines records — the executable half of the
//! CONTRACT.md discipline.
//!
//! Every crate's conformance suite registers named cases; running the
//! suite executes them in registration order and emits ONE JSON line per
//! case carrying the stable case id, the inputs digest, the tolerance
//! specification, the verdict, and evidence pointers — so any failure is
//! reproducible from its log line alone (CONVENTIONS "structured events")
//! and a reimplementation can be held to its predecessor's suite by
//! replaying the same case ids and digests.
//!
//! The v0 case format is DATA-FIRST so the planned fs-ir front end wraps
//! additively: [`Case`] execution returns a [`CaseOutcome`] value and
//! [`SuiteReport`] holds typed [`CaseRecord`]s — printing is a thin,
//! separable layer, and an IR-speaking runner can consume the same
//! records without rewriting suites.
//!
//! No-claims: the harness runs and records; it does not decide what a
//! tolerance MEANS physically, does not compare across ISAs by itself
//! (cross-ISA evidence is a Gauntlet lane running the same suite on both
//! hosts), and awards no certification tiers (that is fs-conform's
//! converter scope).

use core::fmt::{self, Write as _};

/// Version stamped into every emitted record.
pub const CASEBOOK_RECORD_VERSION: u32 = 1;

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
}

impl CaseOutcome {
    /// A passing outcome.
    #[must_use]
    pub fn pass(details: impl Into<String>) -> CaseOutcome {
        CaseOutcome {
            pass: true,
            details: details.into(),
            evidence: Vec::new(),
        }
    }

    /// A failing outcome.
    #[must_use]
    pub fn fail(details: impl Into<String>) -> CaseOutcome {
        CaseOutcome {
            pass: false,
            details: details.into(),
            evidence: Vec::new(),
        }
    }

    /// Attach an evidence pointer.
    #[must_use]
    pub fn with_evidence(mut self, pointer: impl Into<String>) -> CaseOutcome {
        self.evidence.push(pointer.into());
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
        for case in self.cases {
            let structural_defect = if case.id.is_empty() {
                Some("empty case id".to_owned())
            } else if seen.contains(&case.id) {
                Some(format!("duplicate case id {:?}", case.id))
            } else {
                seen.push(case.id);
                None
            };
            let outcome = match structural_defect {
                Some(defect) => CaseOutcome::fail(defect),
                None => (case.run)(),
            };
            let record = CaseRecord {
                version: CASEBOOK_RECORD_VERSION,
                suite: self.name.to_owned(),
                case: case.id.to_owned(),
                inputs_digest: format!("{:016x}", case.inputs_digest),
                tolerance: case.tolerance.to_string(),
                pass: outcome.pass,
                details: outcome.details,
                evidence: outcome.evidence,
            };
            println!("{}", record.json_line());
            records.push(record);
        }
        SuiteReport { records }
    }
}

/// The typed result of one suite run.
#[derive(Debug, Clone)]
pub struct SuiteReport {
    /// One record per executed case, in execution order.
    pub records: Vec<CaseRecord>,
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
        }
        panic!("{msg}");
    }
}
