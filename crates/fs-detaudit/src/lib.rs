//! G5 determinism-audit harness (bead 6nb.2 sibling, epic 6nb: bead
//! frankensim-epic-gauntlet-6nb.6): bit-identity audits across runs and
//! worker counts on one ISA, first-divergence localization over staged
//! hash traces, the cross-ISA divergence classification report, and
//! ExecMode fast/deterministic delta measurement.
//!
//! The engine is subject-agnostic: a [`Subject`] runs a workload at a
//! requested worker count and returns a [`StagedTrace`] — an ordered list
//! of (stage label, content hash) pairs hashed over exact result bits.
//! [`audit`] executes the worker-count matrix with repeats, requires every
//! trace to be bit-identical to the baseline, and on divergence names the
//! FIRST differing stage — the locator that turns "something drifted"
//! into "this reduction, this stage".
//!
//! Cross-ISA: [`classify_cross_isa`] compares two per-ISA artifact
//! ledgers under an explicit [`DivergencePolicy`]. Every difference must
//! land in a DECLARED category (FMA contraction, bounded libm ULP —
//! verified against value bits when supplied) or the report is not clean:
//! unclassified rows and envelope violations are build failures, and
//! reduction-shape divergence is forbidden outright in deterministic
//! mode. [`CrossIsaReport::render_markdown`] emits the documentation-of-
//! record artifact.
//!
//! No-claims: the harness audits what subjects hand it — it does not
//! schedule threads itself, does not prove the absence of nondeterminism
//! beyond the exercised matrix and repeats, makes no throughput claim
//! (mode deltas are MEASURED per call, never assumed), and the cross-ISA
//! report is evidence of classified divergence, not a cross-ISA equality
//! certificate.

use core::fmt::{self, Write as _};
use std::collections::BTreeMap;
use std::time::Instant;

fn escape_json_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{0008}' => escaped.push_str("\\b"),
            '\u{000c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '\u{0000}'..='\u{001f}' => {
                write!(&mut escaped, "\\u{:04x}", u32::from(ch))
                    .expect("writing to a String cannot fail");
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

/// FNV-1a 64 over raw bytes — the artifact/content hash helper.
#[must_use]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// One stage of a subject run: a label and the content hash of every
/// result bit the stage produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageHash {
    /// Stable stage label (names the kernel/reduction/publication step).
    pub label: String,
    /// Content hash over exact result bits at this stage.
    pub hash: u64,
}

/// The ordered staged hash trace of one subject run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StagedTrace {
    /// Stages in execution order.
    pub stages: Vec<StageHash>,
}

impl StagedTrace {
    /// The final artifact hash (last stage), if any stage ran.
    #[must_use]
    pub fn final_hash(&self) -> Option<u64> {
        self.stages.last().map(|s| s.hash)
    }
}

/// A determinism-audit subject: runs the workload at a worker count and
/// reports its staged trace. A subject claiming ExecMode::Deterministic
/// must produce bit-identical traces for every worker count and repeat.
pub struct Subject {
    /// Stable subject name.
    pub name: &'static str,
    /// Execute at the given worker count.
    pub run: Box<dyn Fn(usize) -> StagedTrace>,
}

/// The worker-count matrix an audit exercises.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WorkerMatrix {
    /// Distinct worker counts, ascending.
    pub counts: Vec<usize>,
}

impl WorkerMatrix {
    /// The bead's canonical matrix for this host: 1, 2, P, P+2, and
    /// oversubscribed 2P (deduplicated, ascending), with P the available
    /// parallelism.
    #[must_use]
    pub fn host_default() -> WorkerMatrix {
        let p = std::thread::available_parallelism().map_or(4, usize::from);
        let mut counts = vec![1, 2, p, p + 2, 2 * p];
        counts.sort_unstable();
        counts.dedup();
        WorkerMatrix { counts }
    }

    /// An explicit matrix.
    ///
    /// # Panics
    /// If `counts` is empty or contains zero.
    #[must_use]
    pub fn explicit(mut counts: Vec<usize>) -> WorkerMatrix {
        assert!(
            !counts.is_empty() && counts.iter().all(|&c| c > 0),
            "worker matrix must be non-empty and positive"
        );
        counts.sort_unstable();
        counts.dedup();
        WorkerMatrix { counts }
    }
}

/// Audit configuration: the matrix and the repeats per worker count.
#[derive(Debug, Clone)]
pub struct AuditConfig {
    /// Worker counts to exercise.
    pub matrix: WorkerMatrix,
    /// Runs per worker count (>= 1).
    pub repeats: usize,
}

impl AuditConfig {
    /// A config over the host-default matrix.
    #[must_use]
    pub fn host_default(repeats: usize) -> AuditConfig {
        AuditConfig {
            matrix: WorkerMatrix::host_default(),
            repeats: repeats.max(1),
        }
    }
}

/// Where a diverging run first left the baseline trace.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DivergenceLocator {
    /// Worker count of the diverging run.
    pub workers: usize,
    /// Repeat index of the diverging run.
    pub repeat: usize,
    /// Index of the FIRST differing stage (or first missing stage on a
    /// length mismatch).
    pub first_stage: usize,
    /// The baseline's label at that stage (empty when the baseline is
    /// shorter than the diverging trace).
    pub stage_label: String,
    /// Baseline hash at that stage (0 when absent).
    pub baseline_hash: u64,
    /// Observed hash at that stage (0 when absent).
    pub observed_hash: u64,
}

/// One executed audit row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditRow {
    /// Worker count.
    pub workers: usize,
    /// Repeat index.
    pub repeat: usize,
    /// Final artifact hash of the run.
    pub final_hash: u64,
    /// Whether the full staged trace matched the baseline bit-for-bit.
    pub identical: bool,
}

/// The audit verdict and evidence.
#[derive(Debug, Clone)]
pub struct AuditReport {
    /// The subject name.
    pub subject: &'static str,
    /// The baseline (first matrix entry, repeat 0) staged trace.
    pub baseline: StagedTrace,
    /// Every executed row (matrix × repeats), in execution order.
    pub rows: Vec<AuditRow>,
    /// Every divergence, localized to its first differing stage.
    pub divergences: Vec<DivergenceLocator>,
}

impl AuditReport {
    /// Whether every run reproduced the baseline bit-for-bit.
    #[must_use]
    pub fn identical(&self) -> bool {
        self.divergences.is_empty() && self.rows.iter().all(|r| r.identical)
    }

    /// One JSON line per row plus one per divergence (structured log).
    #[must_use]
    pub fn json_lines(&self) -> Vec<String> {
        let mut out = Vec::with_capacity(self.rows.len() + self.divergences.len());
        let subject = escape_json_string(self.subject);
        for row in &self.rows {
            out.push(format!(
                "{{\"detaudit\":\"row\",\"subject\":\"{}\",\"workers\":{},\"repeat\":{},\"final_hash\":\"{:016x}\",\"identical\":{}}}",
                subject,
                row.workers,
                row.repeat,
                row.final_hash,
                row.identical
            ));
        }
        for d in &self.divergences {
            out.push(format!(
                "{{\"detaudit\":\"divergence\",\"subject\":\"{}\",\"workers\":{},\"repeat\":{},\"first_stage\":{},\"stage_label\":\"{}\",\"baseline\":\"{:016x}\",\"observed\":\"{:016x}\"}}",
                subject,
                d.workers,
                d.repeat,
                d.first_stage,
                escape_json_string(&d.stage_label),
                d.baseline_hash,
                d.observed_hash
            ));
        }
        out
    }
}

fn locate(
    baseline: &StagedTrace,
    observed: &StagedTrace,
    workers: usize,
    repeat: usize,
) -> Option<DivergenceLocator> {
    let n = baseline.stages.len().max(observed.stages.len());
    for i in 0..n {
        let b = baseline.stages.get(i);
        let o = observed.stages.get(i);
        let same = match (b, o) {
            (Some(b), Some(o)) => b == o,
            _ => false,
        };
        if !same {
            return Some(DivergenceLocator {
                workers,
                repeat,
                first_stage: i,
                stage_label: b.map_or_else(String::new, |s| s.label.clone()),
                baseline_hash: b.map_or(0, |s| s.hash),
                observed_hash: o.map_or(0, |s| s.hash),
            });
        }
    }
    None
}

/// Run the same-ISA determinism audit: every (worker count, repeat) run
/// must reproduce the baseline staged trace bit-for-bit; divergences are
/// localized to their first differing stage and reported, never hidden.
#[must_use]
pub fn audit(subject: &Subject, config: &AuditConfig) -> AuditReport {
    let baseline_workers = config.matrix.counts[0];
    let baseline = (subject.run)(baseline_workers);
    let mut rows = Vec::new();
    let mut divergences = Vec::new();
    for &workers in &config.matrix.counts {
        for repeat in 0..config.repeats {
            let trace = if workers == baseline_workers && repeat == 0 {
                baseline.clone()
            } else {
                (subject.run)(workers)
            };
            let identical = trace == baseline;
            rows.push(AuditRow {
                workers,
                repeat,
                final_hash: trace.final_hash().unwrap_or(0),
                identical,
            });
            if !identical && let Some(loc) = locate(&baseline, &trace, workers, repeat) {
                divergences.push(loc);
            }
        }
    }
    AuditReport {
        subject: subject.name,
        baseline,
        rows,
        divergences,
    }
}

/// The documented cross-ISA divergence categories (plan §5.4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DivergenceClass {
    /// FMA contraction differences (declared and accepted).
    FmaContraction,
    /// libm ULP differences bounded by the fs-math contract envelope;
    /// verified against value bits when both ledgers supply them.
    LibmUlp {
        /// Maximum admitted ULP distance.
        max_ulps: u32,
    },
}

impl fmt::Display for DivergenceClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FmaContraction => write!(f, "fma-contraction"),
            Self::LibmUlp { max_ulps } => write!(f, "libm-ulp<={max_ulps}"),
        }
    }
}

/// One artifact row in a per-ISA ledger.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LedgerRow {
    /// Content hash over the artifact's exact bits.
    pub hash: u64,
    /// The artifact's scalar value bits, when the artifact IS one f64
    /// (enables ULP verification); `None` for composite artifacts.
    pub value_bits: Option<u64>,
}

/// A per-ISA artifact ledger: artifact name → row.
#[derive(Debug, Clone)]
pub struct IsaLedger {
    /// The ISA label (e.g. "aarch64", "x86-64").
    pub isa: String,
    /// Artifact rows.
    pub rows: BTreeMap<String, LedgerRow>,
}

/// The declared admissible divergence per artifact. Artifacts absent from
/// the policy must match bit-for-bit; reduction-shape divergence has NO
/// admissible declaration in deterministic mode (plan: "should be NONE"),
/// so any unexplained mismatch fails the report.
#[derive(Debug, Clone, Default)]
pub struct DivergencePolicy {
    /// Artifact name → declared class.
    pub declared: BTreeMap<String, DivergenceClass>,
}

/// How one artifact classified.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Classification {
    /// Bit-identical across ISAs.
    Identical,
    /// Divergent inside its declared category (detail records the check).
    Classified {
        /// The declared class.
        class: DivergenceClass,
        /// The verification detail.
        detail: String,
    },
    /// Divergent with no admissible declaration, a violated envelope, or
    /// a ledger-shape mismatch — a build failure.
    Unclassified {
        /// The reason.
        reason: String,
    },
}

/// One cross-ISA report row.
#[derive(Debug, Clone, PartialEq)]
pub struct CrossIsaRow {
    /// The artifact name.
    pub artifact: String,
    /// Hash on the first ledger's ISA (0 when absent).
    pub hash_a: u64,
    /// Hash on the second ledger's ISA (0 when absent).
    pub hash_b: u64,
    /// The classification.
    pub classification: Classification,
}

/// The cross-ISA divergence report: the documentation-of-record artifact.
#[derive(Debug, Clone)]
pub struct CrossIsaReport {
    /// First ISA label.
    pub isa_a: String,
    /// Second ISA label.
    pub isa_b: String,
    /// One row per artifact across both ledgers.
    pub rows: Vec<CrossIsaRow>,
}

fn ulp_distance(a_bits: u64, b_bits: u64) -> Option<u64> {
    let a = f64::from_bits(a_bits);
    let b = f64::from_bits(b_bits);
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    if a.is_sign_negative() != b.is_sign_negative() && a != 0.0 && b != 0.0 {
        return None;
    }
    // Monotone integer mapping of finite f64 within one sign.
    let key = |bits: u64| -> i128 {
        let v = f64::from_bits(bits);
        if v.is_sign_negative() {
            -i128::from(bits & 0x7fff_ffff_ffff_ffff)
        } else {
            i128::from(bits)
        }
    };
    Some(
        key(a_bits)
            .abs_diff(key(b_bits))
            .try_into()
            .unwrap_or(u64::MAX),
    )
}

/// Classify every artifact difference between two per-ISA ledgers under
/// the declared policy. Every divergence must land in a declared category
/// with a verified envelope; anything else is `Unclassified` and the
/// report is not clean.
#[must_use]
pub fn classify_cross_isa(
    a: &IsaLedger,
    b: &IsaLedger,
    policy: &DivergencePolicy,
) -> CrossIsaReport {
    let mut names: Vec<&String> = a.rows.keys().chain(b.rows.keys()).collect();
    names.sort();
    names.dedup();
    let mut rows = Vec::with_capacity(names.len());
    for name in names {
        let ra = a.rows.get(name);
        let rb = b.rows.get(name);
        let classification = match (ra, rb) {
            (None, _) | (_, None) => Classification::Unclassified {
                reason: "artifact missing from one ISA ledger".to_owned(),
            },
            (Some(ra), Some(rb)) if ra.hash == rb.hash => Classification::Identical,
            (Some(ra), Some(rb)) => match policy.declared.get(name) {
                None => Classification::Unclassified {
                    reason: "divergence with no declared category (reduction-shape \
                             divergence is forbidden in deterministic mode)"
                        .to_owned(),
                },
                Some(class @ DivergenceClass::FmaContraction) => Classification::Classified {
                    class: *class,
                    detail: "declared FMA-contraction divergence accepted".to_owned(),
                },
                Some(class @ DivergenceClass::LibmUlp { max_ulps }) => {
                    match (ra.value_bits, rb.value_bits) {
                        (Some(va), Some(vb)) => match ulp_distance(va, vb) {
                            Some(d) if d <= u64::from(*max_ulps) => Classification::Classified {
                                class: *class,
                                detail: format!("measured {d} ulps within envelope"),
                            },
                            Some(d) => Classification::Unclassified {
                                reason: format!("libm envelope violated: {d} ulps > {max_ulps}"),
                            },
                            None => Classification::Unclassified {
                                reason: "libm class needs comparable finite same-sign values"
                                    .to_owned(),
                            },
                        },
                        _ => Classification::Unclassified {
                            reason: "libm class declared but value bits missing from a ledger"
                                .to_owned(),
                        },
                    }
                }
            },
        };
        rows.push(CrossIsaRow {
            artifact: name.clone(),
            hash_a: ra.map_or(0, |r| r.hash),
            hash_b: rb.map_or(0, |r| r.hash),
            classification,
        });
    }
    CrossIsaReport {
        isa_a: a.isa.clone(),
        isa_b: b.isa.clone(),
        rows,
    }
}

impl CrossIsaReport {
    /// Whether every row is identical or classified (no unclassified rows
    /// — the bead's acceptance gate).
    #[must_use]
    pub fn clean(&self) -> bool {
        self.rows
            .iter()
            .all(|r| !matches!(r.classification, Classification::Unclassified { .. }))
    }

    /// Render the documentation-of-record markdown artifact.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "# Cross-ISA divergence report: {} vs {}\n\n\
             | artifact | {} | {} | classification |\n|---|---|---|---|\n",
            self.isa_a, self.isa_b, self.isa_a, self.isa_b
        ));
        for row in &self.rows {
            let class = match &row.classification {
                Classification::Identical => "identical".to_owned(),
                Classification::Classified { class, detail } => format!("{class}: {detail}"),
                Classification::Unclassified { reason } => format!("UNCLASSIFIED: {reason}"),
            };
            out.push_str(&format!(
                "| {} | {:016x} | {:016x} | {} |\n",
                row.artifact, row.hash_a, row.hash_b, class
            ));
        }
        out.push_str(&format!(
            "\nverdict: {}\n",
            if self.clean() {
                "clean (every divergence classified)"
            } else {
                "NOT CLEAN (unclassified rows are build failures)"
            }
        ));
        out
    }
}

/// The measured ExecMode fast/deterministic tradeoff for one workload:
/// the throughput gain and the reproducibility loss, OBSERVED per call.
#[derive(Debug, Clone)]
pub struct ModeDeltaReport {
    /// Workload name.
    pub name: &'static str,
    /// Total deterministic-mode wall time over the repeats.
    pub deterministic_ns: u128,
    /// Total fast-mode wall time over the repeats.
    pub fast_ns: u128,
    /// deterministic_ns / fast_ns (> 1 means fast mode is faster).
    pub gain_ratio: f64,
    /// Whether deterministic-mode hashes were identical across repeats.
    pub deterministic_reproducible: bool,
    /// Whether fast-mode hashes were identical across repeats.
    pub fast_reproducible: bool,
}

/// Measure the fast-versus-deterministic tradeoff: times both closures
/// over `repeats` runs and records whether each mode's artifact hash was
/// reproducible — so the tradeoff stays observed rather than assumed.
#[must_use]
pub fn measure_mode_delta(
    name: &'static str,
    repeats: usize,
    deterministic: impl Fn() -> u64,
    fast: impl Fn() -> u64,
) -> ModeDeltaReport {
    let repeats = repeats.max(2);
    let mut det_hashes = Vec::with_capacity(repeats);
    let det_start = Instant::now();
    for _ in 0..repeats {
        det_hashes.push(deterministic());
    }
    let deterministic_ns = det_start.elapsed().as_nanos();
    let mut fast_hashes = Vec::with_capacity(repeats);
    let fast_start = Instant::now();
    for _ in 0..repeats {
        fast_hashes.push(fast());
    }
    let fast_ns = fast_start.elapsed().as_nanos();
    let gain_ratio = if fast_ns == 0 {
        f64::INFINITY
    } else {
        let d = deterministic_ns as f64;
        let f = fast_ns as f64;
        d / f
    };
    ModeDeltaReport {
        name,
        deterministic_ns,
        fast_ns,
        gain_ratio,
        deterministic_reproducible: det_hashes.windows(2).all(|w| w[0] == w[1]),
        fast_reproducible: fast_hashes.windows(2).all(|w| w[0] == w[1]),
    }
}
