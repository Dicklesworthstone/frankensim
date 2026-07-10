//! fs-obs — structured observability: the ONE event schema for kernels,
//! solvers, test suites, and (once it lands) the ledger `events` table.
//!
//! Roughly forty beads specify per-suite logging ("structured JSON records
//! sufficient to reproduce any failure from logs alone"). Without a single
//! owned schema, forty suites invent forty dialects and that promise dies.
//! This crate makes "diagnosable from logs alone" a CHECKABLE property:
//! every emitter produces [`Event`]s, every event serializes to one JSON-line
//! dialect, and the [`validate_line`] / [`lint_failure_record`] gates run in
//! CI.
//!
//! Determinism split (Decalogue P2): an event's CONTENT (kind, payload,
//! scope, seq) is deterministic in deterministic mode and is what
//! [`Event::content_hash`] covers; wall-clock time lives in the envelope
//! only, excluded from the hash, so logs from two runs of the same seed diff
//! cleanly.
//!
//! Serialization is in-house (Decalogue P1: std + constellation only — serde
//! is not on that list). The wire format is JSON-lines with CANONICAL field
//! order; the strict validator treats deviation as corruption, not dialect.
//!
//! The [`ident`] module owns the CANONICAL REPLAY IDENTITY encoding
//! (bead gp3.14): versioned, typed, length-prefixed — the shared
//! replacement for ad hoc delimiter-concatenation identities.

pub mod ident;

use core::fmt;
use std::fmt::Write as _;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Schema version stamped into every line; bump on any non-additive change.
pub const SCHEMA_VERSION: u32 = 1;

/// Severity ladder. `Error` events MUST satisfy [`lint_failure_record`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// High-volume diagnostics (tile completions).
    Trace,
    /// Normal progress (solver iterations, case verdicts).
    Info,
    /// Something degraded but the run continues.
    Warn,
    /// A failure record — must be reproducible from its own payload.
    Error,
}

impl Severity {
    fn name(self) -> &'static str {
        match self {
            Severity::Trace => "trace",
            Severity::Info => "info",
            Severity::Warn => "warn",
            Severity::Error => "error",
        }
    }
}

/// Typed payload registry, v1. Additive evolution only: new kinds may be
/// added; existing fields may never change meaning (validator-enforced by
/// the golden lines in the conformance suite).
#[derive(Debug, Clone, PartialEq)]
pub enum EventKind {
    /// One iteration of an iterative solver.
    SolverResidual {
        /// Solver identity (e.g. "cg", "p-mg").
        solver: String,
        /// Iteration index.
        iter: u64,
        /// Residual norm.
        residual: f64,
    },
    /// A tile finished executing.
    TileComplete {
        /// Logical tile identity (the SAME identity that keys RNG streams
        /// and deterministic reductions — plan §5.2).
        tile: u64,
        /// Kernel name.
        kernel: String,
    },
    /// A scope was cancelled.
    Cancellation {
        /// Why (kill-handle, budget, panic containment, ...).
        reason: String,
    },
    /// Budget accounting delta (P4: budgets first).
    BudgetDelta {
        /// Resource name ("wall_s", "mem_bytes", "core_s", "energy_j").
        resource: String,
        /// Amount spent in this step.
        spent: f64,
        /// Remaining grant.
        remaining: f64,
    },
    /// A gradient verification outcome (the merge-gate evidence, §8.7).
    GradientCheck {
        /// Operator under test.
        op: String,
        /// Max relative error across probed directions.
        max_rel_err: f64,
        /// Verdict.
        pass: bool,
    },
    /// One conformance-suite case verdict (plan §13.3).
    ConformanceCase {
        /// Suite id (e.g. "fs-qty/conformance").
        suite: String,
        /// Case id (e.g. "qty-001/0.12Pa*s").
        case: String,
        /// Verdict.
        pass: bool,
        /// Human/agent-readable detail; REQUIRED non-empty when `pass=false`
        /// (lint-enforced) so failures reproduce from the log alone.
        detail: String,
        /// Replay seed when the case is randomized (lint: required on fail
        /// for randomized cases; 0 = not randomized).
        seed: u64,
    },
    /// A performance measurement (roofline harness rows).
    BenchmarkResult {
        /// Kernel name.
        kernel: String,
        /// Metric name ("glups", "gflops", "bandwidth_gbs", "mrays").
        metric: String,
        /// Measured value.
        value: f64,
        /// Machine fingerprint hash (fs-substrate probe).
        machine: u64,
    },
    /// A chaos/storm assertion outcome (G4).
    StormAssertion {
        /// Assertion name ("no-arena-leak", "cancel-latency-p99").
        name: String,
        /// Verdict.
        pass: bool,
        /// Storm seed for replay.
        seed: u64,
    },
    /// Escape hatch for kinds not yet in the registry; the payload must be a
    /// single pre-serialized JSON object (validated for balance, not schema).
    Custom {
        /// Kind name (kebab-case).
        name: String,
        /// Pre-serialized JSON object.
        json: String,
    },
}

impl EventKind {
    fn kind_name(&self) -> &'static str {
        match self {
            EventKind::SolverResidual { .. } => "solver_residual",
            EventKind::TileComplete { .. } => "tile_complete",
            EventKind::Cancellation { .. } => "cancellation",
            EventKind::BudgetDelta { .. } => "budget_delta",
            EventKind::GradientCheck { .. } => "gradient_check",
            EventKind::ConformanceCase { .. } => "conformance_case",
            EventKind::BenchmarkResult { .. } => "benchmark_result",
            EventKind::StormAssertion { .. } => "storm_assertion",
            EventKind::Custom { .. } => "custom",
        }
    }
}

/// One observability event: envelope + typed payload.
#[derive(Debug, Clone, PartialEq)]
pub struct Event {
    /// Session identity (study/run scope; content-hashed).
    pub session: String,
    /// Slash-separated scope path mirroring the asupersync scope tree
    /// (e.g. "study-x/op-3/kernel-lbm/tile-42"); content-hashed.
    pub scope: String,
    /// Per-emitter monotone sequence number; content-hashed (gives logs a
    /// deterministic total order per scope without wall-clock).
    pub seq: u64,
    /// Severity.
    pub severity: Severity,
    /// Typed payload.
    pub kind: EventKind,
    /// Wall-clock nanoseconds since the unix epoch. ENVELOPE ONLY: excluded
    /// from `content_hash`, always serialized LAST, `None` in deterministic
    /// replay comparisons.
    pub wall_ns: Option<u64>,
}

fn esc(out: &mut String, s: &str) {
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
}

fn push_str_field(out: &mut String, key: &str, val: &str) {
    let _ = write!(out, "\"{key}\":\"");
    esc(out, val);
    out.push('"');
}

/// Serialize a float for the wire: finite → shortest-round-trip; non-finite
/// → tagged string (JSON has no NaN/Inf; readers must handle both shapes).
fn push_f64(out: &mut String, key: &str, v: f64) {
    if v.is_finite() {
        let _ = write!(out, "\"{key}\":{v}");
    } else {
        let _ = write!(out, "\"{key}\":\"non-finite:{v}\"");
    }
}

impl Event {
    /// Serialize the CONTENT portion (everything except `wall_ns`) in
    /// canonical field order.
    fn content_json(&self) -> String {
        let mut s = String::with_capacity(160);
        let _ = write!(s, "{{\"v\":{SCHEMA_VERSION},");
        push_str_field(&mut s, "session", &self.session);
        s.push(',');
        push_str_field(&mut s, "scope", &self.scope);
        let _ = write!(s, ",\"seq\":{},", self.seq);
        push_str_field(&mut s, "severity", self.severity.name());
        s.push(',');
        push_str_field(&mut s, "kind", self.kind.kind_name());
        s.push_str(",\"payload\":{");
        match &self.kind {
            EventKind::SolverResidual {
                solver,
                iter,
                residual,
            } => {
                push_str_field(&mut s, "solver", solver);
                let _ = write!(s, ",\"iter\":{iter},");
                push_f64(&mut s, "residual", *residual);
            }
            EventKind::TileComplete { tile, kernel } => {
                let _ = write!(s, "\"tile\":{tile},");
                push_str_field(&mut s, "kernel", kernel);
            }
            EventKind::Cancellation { reason } => {
                push_str_field(&mut s, "reason", reason);
            }
            EventKind::BudgetDelta {
                resource,
                spent,
                remaining,
            } => {
                push_str_field(&mut s, "resource", resource);
                s.push(',');
                push_f64(&mut s, "spent", *spent);
                s.push(',');
                push_f64(&mut s, "remaining", *remaining);
            }
            EventKind::GradientCheck {
                op,
                max_rel_err,
                pass,
            } => {
                push_str_field(&mut s, "op", op);
                s.push(',');
                push_f64(&mut s, "max_rel_err", *max_rel_err);
                let _ = write!(s, ",\"pass\":{pass}");
            }
            EventKind::ConformanceCase {
                suite,
                case,
                pass,
                detail,
                seed,
            } => {
                push_str_field(&mut s, "suite", suite);
                s.push(',');
                push_str_field(&mut s, "case", case);
                let _ = write!(s, ",\"pass\":{pass},");
                push_str_field(&mut s, "detail", detail);
                let _ = write!(s, ",\"seed\":{seed}");
            }
            EventKind::BenchmarkResult {
                kernel,
                metric,
                value,
                machine,
            } => {
                push_str_field(&mut s, "kernel", kernel);
                s.push(',');
                push_str_field(&mut s, "metric", metric);
                s.push(',');
                push_f64(&mut s, "value", *value);
                let _ = write!(s, ",\"machine\":{machine}");
            }
            EventKind::StormAssertion { name, pass, seed } => {
                push_str_field(&mut s, "name", name);
                let _ = write!(s, ",\"pass\":{pass},\"seed\":{seed}");
            }
            EventKind::Custom { name, json } => {
                push_str_field(&mut s, "name", name);
                let _ = write!(s, ",\"data\":{json}");
            }
        }
        s.push('}');
        s.push('}');
        s
    }

    /// Serialize one JSON-line: content with `wall_ns` spliced in as the
    /// LAST field before the closing brace (envelope position).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let mut content = self.content_json();
        if let Some(ns) = self.wall_ns {
            content.pop(); // strip trailing '}'
            let _ = write!(content, ",\"wall_ns\":{ns}}}");
        }
        content
    }

    /// Deterministic content hash (FNV-1a 64 over the content JSON —
    /// EXCLUDING wall-clock). Not cryptographic; ledger-grade content
    /// addressing (BLAKE3-class) arrives with fs-ledger and this method's
    /// contract permits strengthening the hash at a schema-version bump.
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        fnv1a64(self.content_json().as_bytes())
    }
}

/// FNV-1a 64-bit (in-house, deterministic across platforms).
#[must_use]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// A per-scope emitter handing out monotone sequence numbers.
#[derive(Debug)]
pub struct Emitter {
    session: String,
    scope: String,
    seq: u64,
}

impl Emitter {
    /// Create an emitter for one (session, scope) pair.
    #[must_use]
    pub fn new(session: impl Into<String>, scope: impl Into<String>) -> Self {
        Emitter {
            session: session.into(),
            scope: scope.into(),
            seq: 0,
        }
    }

    /// Build the next event (seq auto-increments). `wall_ns` is supplied by
    /// the caller because THIS crate must stay deterministic and I/O-free;
    /// runtime layers pass real clocks, tests pass `None`.
    pub fn emit(&mut self, severity: Severity, kind: EventKind, wall_ns: Option<u64>) -> Event {
        let e = Event {
            session: self.session.clone(),
            scope: self.scope.clone(),
            seq: self.seq,
            severity,
            kind,
            wall_ns,
        };
        self.seq += 1;
        e
    }
}

/// Validation failure for a JSON-line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaError {
    /// Byte position (approximate for structural errors).
    pub at: usize,
    /// What is wrong and how to fix it.
    pub message: String,
}

impl fmt::Display for SchemaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "event schema violation at byte {}: {}",
            self.at, self.message
        )
    }
}

impl core::error::Error for SchemaError {}

/// Registry of known kinds (validator uses it; keep in sync with EventKind).
pub const KNOWN_KINDS: &[&str] = &[
    "solver_residual",
    "tile_complete",
    "cancellation",
    "budget_delta",
    "gradient_check",
    "conformance_case",
    "benchmark_result",
    "storm_assertion",
    "custom",
];

/// Strict structural validation of one JSON-line: required envelope keys in
/// canonical order, known kind, balanced payload object. (The writer is
/// ours; deviation means corruption, so the checks are cheap and strict
/// rather than a full JSON parser.)
///
/// # Errors
/// Returns [`SchemaError`] naming the first violated requirement.
pub fn validate_line(line: &str) -> Result<(), SchemaError> {
    let need = |cond: bool, at: usize, msg: &str| -> Result<(), SchemaError> {
        if cond {
            Ok(())
        } else {
            Err(SchemaError {
                at,
                message: msg.to_string(),
            })
        }
    };
    need(
        line.starts_with('{') && line.ends_with('}'),
        0,
        "line must be one JSON object",
    )?;
    let ver = format!("{{\"v\":{SCHEMA_VERSION},");
    need(
        line.starts_with(&ver),
        0,
        "first field must be the schema version \"v\"",
    )?;
    for key in [
        "\"session\":",
        "\"scope\":",
        "\"seq\":",
        "\"severity\":",
        "\"kind\":",
        "\"payload\":{",
    ] {
        need(
            line.contains(key),
            0,
            &format!("missing required field {key}"),
        )?;
    }
    let kind_pos = line.find("\"kind\":\"").ok_or(SchemaError {
        at: 0,
        message: "kind must be a string".to_string(),
    })?;
    let kind_rest = &line[kind_pos + 8..];
    let kind_end = kind_rest.find('"').unwrap_or(0);
    let kind = &kind_rest[..kind_end];
    need(
        KNOWN_KINDS.contains(&kind),
        kind_pos,
        &format!("unknown kind {kind:?}; register new kinds in the payload registry"),
    )?;
    // Braces balance (escaped quotes handled by the writer's escaping).
    let mut depth = 0i64;
    let mut in_str = false;
    let mut prev_backslash = false;
    for (i, c) in line.char_indices() {
        match c {
            '"' if !prev_backslash => in_str = !in_str,
            '{' if !in_str => depth += 1,
            '}' if !in_str => depth -= 1,
            _ => {}
        }
        prev_backslash = c == '\\' && !prev_backslash;
        if depth < 0 {
            return Err(SchemaError {
                at: i,
                message: "unbalanced braces".to_string(),
            });
        }
    }
    need(
        depth == 0 && !in_str,
        line.len(),
        "unbalanced braces or unterminated string",
    )?;
    Ok(())
}

/// The failure-record completeness lint, v1 (the "log-replay lint" of the
/// observability bead): failure verdicts must carry enough to reproduce.
///
/// # Errors
/// Returns [`SchemaError`] describing the missing reproduction ingredient.
pub fn lint_failure_record(event: &Event) -> Result<(), SchemaError> {
    match &event.kind {
        EventKind::ConformanceCase {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing conformance case must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::StormAssertion {
            pass: false, seed, ..
        } if *seed == 0 => Err(SchemaError {
            at: 0,
            message: "failing storm assertion must carry its replay seed".to_string(),
        }),
        EventKind::GradientCheck {
            pass: false, op, ..
        } if op.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing gradient check must name its operator".to_string(),
        }),
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_events() -> Vec<Event> {
        let mut em = Emitter::new("study-x", "op-1/kernel-cg");
        vec![
            em.emit(
                Severity::Info,
                EventKind::SolverResidual {
                    solver: "cg".into(),
                    iter: 3,
                    residual: 1.5e-7,
                },
                Some(1_000),
            ),
            em.emit(
                Severity::Trace,
                EventKind::TileComplete {
                    tile: 42,
                    kernel: "lbm_d3q19".into(),
                },
                None,
            ),
            em.emit(
                Severity::Warn,
                EventKind::BudgetDelta {
                    resource: "wall_s".into(),
                    spent: 12.5,
                    remaining: 7187.5,
                },
                Some(2_000),
            ),
            em.emit(
                Severity::Error,
                EventKind::ConformanceCase {
                    suite: "fs-qty/conformance".into(),
                    case: "qty-001".into(),
                    pass: false,
                    detail: "value mismatch: got 0.13, want 0.12".into(),
                    seed: 7,
                },
                Some(3_000),
            ),
            em.emit(
                Severity::Info,
                EventKind::Custom {
                    name: "regime-report".into(),
                    json: r#"{"re":100.5,"we":0.3}"#.into(),
                },
                None,
            ),
        ]
    }

    #[test]
    fn every_kind_serializes_and_validates() {
        for e in sample_events() {
            let line = e.to_jsonl();
            validate_line(&line).unwrap_or_else(|err| panic!("{line}: {err}"));
        }
    }

    #[test]
    fn wall_clock_is_envelope_only() {
        let mut a = sample_events().remove(0);
        let h1 = a.content_hash();
        a.wall_ns = Some(999_999_999);
        let h2 = a.content_hash();
        assert_eq!(h1, h2, "content hash must exclude wall clock");
        // ...but the serialized line DOES carry it, as the last field.
        assert!(a.to_jsonl().ends_with(",\"wall_ns\":999999999}"));
    }

    #[test]
    fn content_hash_is_sensitive_to_content() {
        let events = sample_events();
        let mut hashes: Vec<u64> = events.iter().map(Event::content_hash).collect();
        hashes.sort_unstable();
        hashes.dedup();
        assert_eq!(
            hashes.len(),
            events.len(),
            "distinct events must hash distinctly"
        );
    }

    #[test]
    fn sequence_numbers_are_monotone_per_emitter() {
        let events = sample_events();
        for (i, e) in events.iter().enumerate() {
            assert_eq!(e.seq, i as u64);
        }
    }

    #[test]
    fn golden_line_shape_is_stable() {
        // Schema evolution is additive-only; this golden line is the contract.
        // Changing it requires a SCHEMA_VERSION bump and a semantic justification
        // (golden-evidence policy, AGENTS.md).
        let e = Event {
            session: "s".into(),
            scope: "a/b".into(),
            seq: 5,
            severity: Severity::Info,
            kind: EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 1e-9,
                pass: true,
            },
            wall_ns: None,
        };
        // Note: Rust's shortest-round-trip float Display never uses scientific
        // notation, so 1e-9 serializes as 0.000000001 — that IS the contract.
        assert_eq!(
            e.to_jsonl(),
            "{\"v\":1,\"session\":\"s\",\"scope\":\"a/b\",\"seq\":5,\"severity\":\"info\",\
             \"kind\":\"gradient_check\",\"payload\":{\"op\":\"poisson\",\
             \"max_rel_err\":0.000000001,\"pass\":true}}"
        );
    }

    #[test]
    fn validator_rejects_corruption() {
        let good = sample_events()[0].to_jsonl();
        for bad in [
            String::new(),
            "not json".to_string(),
            good.replace("\"v\":1", "\"v\":99"),
            good.replace("solver_residual", "mystery_kind"),
            good.replace("\"session\"", "\"sesion\""),
            good[..good.len() - 1].to_string(),
        ] {
            assert!(validate_line(&bad).is_err(), "should reject: {bad}");
        }
    }

    #[test]
    fn failure_lint_demands_reproduction_ingredients() {
        let mut em = Emitter::new("s", "x");
        let bad = em.emit(
            Severity::Error,
            EventKind::ConformanceCase {
                suite: "s".into(),
                case: "c".into(),
                pass: false,
                detail: String::new(),
                seed: 0,
            },
            None,
        );
        assert!(lint_failure_record(&bad).is_err());
        let good = em.emit(
            Severity::Error,
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: false,
                seed: 42,
            },
            None,
        );
        assert!(lint_failure_record(&good).is_ok());
        let bad_storm = em.emit(
            Severity::Error,
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: false,
                seed: 0,
            },
            None,
        );
        assert!(lint_failure_record(&bad_storm).is_err());
    }

    #[test]
    fn escaping_handles_hostile_strings() {
        let mut em = Emitter::new("s\"es\\sion\n", "sc\tope");
        let e = em.emit(
            Severity::Error,
            EventKind::Cancellation {
                reason: "quote\" backslash\\ newline\n tab\t".into(),
            },
            None,
        );
        let line = e.to_jsonl();
        validate_line(&line).unwrap_or_else(|err| panic!("{line}: {err}"));
        assert!(!line.contains('\n'), "JSONL lines must be single-line");
    }

    #[test]
    fn non_finite_floats_are_tagged_not_invalid() {
        let mut em = Emitter::new("s", "x");
        let e = em.emit(
            Severity::Info,
            EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 1,
                residual: f64::NAN,
            },
            None,
        );
        let line = e.to_jsonl();
        validate_line(&line).expect("tagged non-finite must stay valid");
        assert!(line.contains("non-finite:NaN"));
    }

    #[test]
    fn fnv_matches_known_answers() {
        // Published FNV-1a 64 test vectors.
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a64(b"foobar"), 0x85944171f73967e8);
    }
}
