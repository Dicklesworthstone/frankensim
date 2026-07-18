//! IN-HOUSE PROPERTY-BASED TESTING with integrated shrinking (bead
//! frankensim-4nh8; Gauntlet G0, plan §13.1). The plan requires
//! "proptest-class shrinking" without external dependencies: this crate
//! is the shared engine the law suites adopt, extracted from the
//! shrinking convention `fs-bisect::compound` proved on the real
//! powi/rand_nla incident.
//!
//! DESIGN CONTRACT:
//! - **Deterministic**: every case derives from `(suite seed, case index)`
//!   via a splitmix64 stream — same seed, same cases, any machine, any
//!   thread count. A failure prints the CASE SEED; setting
//!   `FSIM_PROPCHECK_REPLAY=<case_seed>` reruns exactly that case.
//!   The same JSONL failure row is appended to a process-scoped replay
//!   artifact (or `FSIM_PROPCHECK_REPLAY_FILE` when configured).
//! - **Shrinking**: on failure the input is minimized by greedy
//!   first-failing-candidate descent (the shared `fs-propcheck` convention,
//!   publicly re-exported by `fs-bisect`: candidates
//!   most-aggressive-first, fixed order, bounded steps), so CI failures
//!   arrive as local-minimum counterexamples, not haystacks.
//! - **Zero runtime dependencies** (UTIL layer): the generator stream is
//!   a test-input source, not a simulation RNG — kernels keep using
//!   fs-rand's Philox discipline; this stream never touches ledgers.
//! - **JSONL failure rows**: every failure emits one structured line
//!   (property, case seed, shrink steps, minimized debug form) before
//!   panicking, per the house verdict-logging style.

/// Gauntlet G3 relation declarations and the metamorphic runner.
pub mod metamorphic;

/// Deterministic generator stream for test inputs (splitmix64).
///
/// NOT a simulation RNG: no Philox identity, no ledger coupling — a
/// reproducible source of test cases only.
#[derive(Debug, Clone)]
pub struct Stream {
    state: u64,
}

impl Stream {
    /// A stream seeded for `(suite_seed, case_index)` — the replay key.
    #[must_use]
    pub fn for_case(suite_seed: u64, case_index: u64) -> Stream {
        // Mix the pair so nearby cases decorrelate.
        let mut s = Stream {
            state: suite_seed ^ case_index.wrapping_mul(0x9e37_79b9_7f4a_7c15),
        };
        // Warm the mixer so low-entropy seeds still spread.
        let _ = s.next_u64();
        s
    }

    /// Next raw 64-bit value (splitmix64 step).
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform in `[lo, hi]` (inclusive), using rejection sampling rather
    /// than biased modulo reduction. `lo > hi` is a caller bug and panics
    /// with the offending bounds (test-harness semantics).
    pub fn int_in(&mut self, lo: i64, hi: i64) -> i64 {
        assert!(lo <= hi, "int_in bounds inverted: [{lo}, {hi}]");
        let span = hi.wrapping_sub(lo).cast_unsigned();
        if span == u64::MAX {
            return self.next_u64().cast_signed();
        }
        let width = span + 1;
        // `threshold` is 2^64 mod width. Rejecting raw values below it
        // leaves an admitted population whose size is exactly divisible by
        // width, so every residue has the same number of preimages.
        let threshold = width.wrapping_neg() % width;
        loop {
            let raw = self.next_u64();
            if raw >= threshold {
                return lo.wrapping_add((raw % width).cast_signed());
            }
        }
    }

    /// Uniform-ish finite f64 in `[lo, hi)`, with a 1-in-8 chance of drawing
    /// a SPECIAL value (0.0, -0.0, either endpoint-neighbor, 1.0, -1.0) —
    /// the corners where numeric laws actually break.
    ///
    /// Bounds must be finite. The convex interpolation avoids overflowing
    /// `hi - lo`, and the final clamp prevents rounding up to the excluded
    /// upper endpoint.
    pub fn f64_in(&mut self, lo: f64, hi: f64) -> f64 {
        const TWO_POW_53: f64 = 9_007_199_254_740_992.0;

        assert!(
            lo.is_finite() && hi.is_finite(),
            "f64_in bounds must be finite: [{lo}, {hi})"
        );
        assert!(lo < hi, "f64_in bounds inverted: [{lo}, {hi})");
        if self.next_u64().is_multiple_of(8) {
            let specials = [0.0, -0.0, lo, next_down(hi), 1.0, -1.0];
            let pick = specials[(self.next_u64() % specials.len() as u64) as usize];
            if pick >= lo && pick < hi {
                return pick;
            }
        }
        let unit = (self.next_u64() >> 11) as f64 / TWO_POW_53;
        let drawn = lo * (1.0 - unit) + hi * unit;
        if drawn < lo {
            lo
        } else if drawn >= hi || !drawn.is_finite() {
            next_down(hi)
        } else {
            drawn
        }
    }

    /// A vector of `n` in `[0, max_len]` elements drawn by `elem`.
    ///
    /// `usize::MAX` is refused explicitly because the inclusive-width
    /// calculation cannot represent `max_len + 1` and such an allocation is
    /// not a meaningful test-harness request.
    pub fn vec_of<T>(&mut self, max_len: usize, mut elem: impl FnMut(&mut Stream) -> T) -> Vec<T> {
        assert!(
            max_len != usize::MAX,
            "vec_of max_len must be less than usize::MAX"
        );
        let n = (self.next_u64() % (max_len as u64 + 1)) as usize;
        (0..n).map(|_| elem(self)).collect()
    }
}

fn next_down(value: f64) -> f64 {
    debug_assert!(value.is_finite());
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    if value > 0.0 {
        f64::from_bits(bits - 1)
    } else {
        f64::from_bits(bits + 1)
    }
}

/// Deterministic shrinking: candidates strictly "smaller" than `self`,
/// most aggressive first, in a FIXED order. Empty = fully shrunk.
///
/// This is the canonical shared trait; `fs-bisect::compound` publicly
/// re-exports it while retaining its own failure-predicate minimizer.
pub trait Shrink: Clone {
    /// Smaller candidate inputs, most aggressive first.
    fn shrink_candidates(&self) -> Vec<Self>;
}

impl Shrink for i64 {
    fn shrink_candidates(&self) -> Vec<i64> {
        let x = *self;
        if x == 0 {
            return vec![];
        }
        let mut out = vec![0];
        // Halving ladder toward zero, then the predecessor.
        let half = x / 2;
        if half != x {
            out.push(half);
        }
        out.push(x - x.signum());
        let magnitude = x.unsigned_abs();
        out.retain(|c| c.unsigned_abs() < magnitude || (c.unsigned_abs() == magnitude && *c > x));
        out.dedup();
        out
    }
}

impl Shrink for u64 {
    fn shrink_candidates(&self) -> Vec<u64> {
        let x = *self;
        if x == 0 {
            return vec![];
        }
        let mut out = vec![0, x / 2, x - 1];
        out.retain(|c| *c < x);
        out.dedup();
        out
    }
}

impl Shrink for f64 {
    fn shrink_candidates(&self) -> Vec<f64> {
        let x = *self;
        if x.is_nan() || x == 0.0 {
            return vec![];
        }
        if x.is_infinite() {
            return vec![0.0, x.signum() * f64::MAX];
        }
        let mut out = vec![0.0, x / 2.0, x.trunc()];
        out.retain(|c| c.abs() < x.abs());
        out.dedup_by(|a, b| a.to_bits() == b.to_bits());
        out
    }
}

impl<T: Shrink> Shrink for Vec<T> {
    fn shrink_candidates(&self) -> Vec<Vec<T>> {
        let mut out = Vec::new();
        if self.is_empty() {
            return out;
        }
        // Most aggressive: drop halves, then single elements, then
        // shrink one element in place.
        let midpoint = self.len() / 2;
        out.push(self[..midpoint].to_vec());
        if midpoint > 0 {
            out.push(self[midpoint..].to_vec());
        }
        if self.len() > 1 {
            for i in 0..self.len() {
                let mut fewer = self.clone();
                fewer.remove(i);
                out.push(fewer);
            }
        }
        for (i, item) in self.iter().enumerate() {
            for cand in item.shrink_candidates() {
                let mut smaller = self.clone();
                smaller[i] = cand;
                out.push(smaller);
            }
        }
        out
    }
}

impl<A: Shrink, B: Shrink> Shrink for (A, B) {
    fn shrink_candidates(&self) -> Vec<(A, B)> {
        let a_candidates = self.0.shrink_candidates();
        let b_candidates = self.1.shrink_candidates();
        let mut out: Vec<(A, B)> = a_candidates
            .iter()
            .cloned()
            .zip(b_candidates.iter().cloned())
            .collect();
        out.extend(a_candidates.into_iter().map(|a| (a, self.1.clone())));
        out.extend(b_candidates.into_iter().map(|b| (self.0.clone(), b)));
        out
    }
}

impl<A: Shrink, B: Shrink, C: Shrink> Shrink for (A, B, C) {
    fn shrink_candidates(&self) -> Vec<(A, B, C)> {
        let a_candidates = self.0.shrink_candidates();
        let b_candidates = self.1.shrink_candidates();
        let c_candidates = self.2.shrink_candidates();
        let mut out: Vec<(A, B, C)> = a_candidates
            .iter()
            .cloned()
            .zip(b_candidates.iter().cloned())
            .zip(c_candidates.iter().cloned())
            .map(|((a, b), c)| (a, b, c))
            .collect();
        out.extend(
            a_candidates
                .iter()
                .cloned()
                .zip(b_candidates.iter().cloned())
                .map(|(a, b)| (a, b, self.2.clone())),
        );
        out.extend(
            a_candidates
                .iter()
                .cloned()
                .zip(c_candidates.iter().cloned())
                .map(|(a, c)| (a, self.1.clone(), c)),
        );
        out.extend(
            b_candidates
                .iter()
                .cloned()
                .zip(c_candidates.iter().cloned())
                .map(|(b, c)| (self.0.clone(), b, c)),
        );
        out.extend(
            a_candidates
                .into_iter()
                .map(|a| (a, self.1.clone(), self.2.clone())),
        );
        out.extend(
            b_candidates
                .into_iter()
                .map(|b| (self.0.clone(), b, self.2.clone())),
        );
        out.extend(
            c_candidates
                .into_iter()
                .map(|c| (self.0.clone(), self.1.clone(), c)),
        );
        out
    }
}

/// Outcome of a shrink descent (mirrors `fs-bisect::MinimizeReport`).
#[derive(Debug, Clone)]
pub struct ShrinkReport<T> {
    /// The smallest failing input found.
    pub minimized: T,
    /// Accepted shrink steps.
    pub steps: usize,
    /// Total property evaluations, including the seed input and any
    /// candidate used to prove that the accepted-step budget was exhausted.
    pub tried: usize,
    /// False when a configured step, evaluation, or per-step candidate
    /// ceiling prevents proving that the retained input is a fixpoint.
    pub converged: bool,
}

/// Default ceiling on property evaluations during one shrink descent,
/// including the failing seed evaluation.
pub const DEFAULT_MAX_MINIMIZE_EVALUATIONS: usize = 100_000;

/// Default ceiling on the candidates inspected at one shrink step.
pub const DEFAULT_MAX_SHRINK_CANDIDATES_PER_STEP: usize = 4_096;

/// Explicit work envelope for [`minimize_with_budget`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinimizeBudget {
    /// Maximum accepted failing descents.
    pub max_steps: usize,
    /// Maximum property evaluations, including the failing seed.
    pub max_evaluations: usize,
    /// Maximum candidates inspected from any one candidate surface.
    pub max_candidates_per_step: usize,
}

impl MinimizeBudget {
    /// The bounded default envelope used by [`minimize`].
    #[must_use]
    pub const fn for_steps(max_steps: usize) -> Self {
        Self {
            max_steps,
            max_evaluations: DEFAULT_MAX_MINIMIZE_EVALUATIONS,
            max_candidates_per_step: DEFAULT_MAX_SHRINK_CANDIDATES_PER_STEP,
        }
    }
}

/// Typed refusal from [`minimize`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MinimizeError {
    /// The supplied seed satisfies the property, so there is no failure to
    /// minimize.
    SeedPasses,
    /// The work envelope cannot even evaluate the supplied seed.
    EmptyEvaluationBudget,
}

impl core::fmt::Display for MinimizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SeedPasses => f.write_str("minimize seed satisfies the property"),
            Self::EmptyEvaluationBudget => {
                f.write_str("minimize requires at least one property evaluation")
            }
        }
    }
}

impl std::error::Error for MinimizeError {}

/// Greedy deterministic minimization: `property(input) == true` means PASS;
/// repeatedly take the FIRST candidate for which it returns false until a
/// fixpoint or `max_steps` accepted steps. Same input + property gives the
/// identical trajectory.
///
/// The candidate surface is inspected once more at the exact accepted-step
/// budget: reaching a fixpoint on the final admitted step reports
/// `converged = true`; finding another failure reports `converged = false`
/// without accepting the over-budget candidate.
///
/// # Errors
/// [`MinimizeError::SeedPasses`] when `seed_failure` is not actually failing.
pub fn minimize<T: Shrink>(
    seed_failure: T,
    property: impl Fn(&T) -> bool,
    max_steps: usize,
) -> Result<ShrinkReport<T>, MinimizeError> {
    minimize_with_budget(seed_failure, property, MinimizeBudget::for_steps(max_steps))
}

/// Greedy deterministic minimization under an explicit property-evaluation
/// and candidate-count work envelope. Reaching either ceiling returns the
/// best retained failing input with `converged = false`.
///
/// Candidate construction is delegated to [`Shrink::shrink_candidates`]; the
/// returned vector is count-checked before any candidate property evaluation.
///
/// # Errors
/// [`MinimizeError::SeedPasses`] when `seed_failure` is not actually failing,
/// or [`MinimizeError::EmptyEvaluationBudget`] when the envelope cannot
/// evaluate the seed.
pub fn minimize_with_budget<T: Shrink>(
    seed_failure: T,
    property: impl Fn(&T) -> bool,
    budget: MinimizeBudget,
) -> Result<ShrinkReport<T>, MinimizeError> {
    if budget.max_evaluations == 0 {
        return Err(MinimizeError::EmptyEvaluationBudget);
    }
    if property(&seed_failure) {
        return Err(MinimizeError::SeedPasses);
    }
    let mut current = seed_failure;
    let mut steps = 0usize;
    let mut tried = 1usize;
    loop {
        let mut next_failure = None;
        let candidates = current.shrink_candidates();
        if candidates.len() > budget.max_candidates_per_step {
            return Ok(ShrinkReport {
                minimized: current,
                steps,
                tried,
                converged: false,
            });
        }
        for cand in candidates {
            if tried == budget.max_evaluations {
                return Ok(ShrinkReport {
                    minimized: current,
                    steps,
                    tried,
                    converged: false,
                });
            }
            tried += 1;
            if !property(&cand) {
                next_failure = Some(cand);
                break;
            }
        }
        let Some(candidate) = next_failure else {
            return Ok(ShrinkReport {
                minimized: current,
                steps,
                tried,
                converged: true,
            });
        };
        if steps == budget.max_steps {
            return Ok(ShrinkReport {
                minimized: current,
                steps,
                tried,
                converged: false,
            });
        }
        current = candidate;
        steps += 1;
    }
}

#[derive(Debug)]
pub(crate) enum StructuredVerdict {
    Pass,
    Fail(StructuredFailure),
}

#[derive(Debug)]
pub(crate) struct StructuredFailure {
    kind: &'static str,
    detail: String,
    context: Vec<(&'static str, String)>,
}

impl StructuredFailure {
    pub(crate) fn new(
        kind: &'static str,
        detail: impl Into<String>,
        context: Vec<(&'static str, String)>,
    ) -> Self {
        Self {
            kind,
            detail: detail.into(),
            context,
        }
    }
}

#[derive(Debug)]
enum PropertyOutcome {
    Pass,
    Failed(StructuredFailure),
    Panicked(String),
}

impl PropertyOutcome {
    fn passed(&self) -> bool {
        matches!(self, Self::Pass)
    }

    fn failure_kind(&self) -> &'static str {
        match self {
            Self::Pass => "non-deterministic",
            Self::Failed(failure) => failure.kind,
            Self::Panicked(_) => "panic",
        }
    }

    fn detail(&self) -> &str {
        match self {
            Self::Pass => "the minimized input passed when re-evaluated",
            Self::Failed(failure) => &failure.detail,
            Self::Panicked(message) => message,
        }
    }

    fn context(&self) -> &[(&'static str, String)] {
        match self {
            Self::Failed(failure) => &failure.context,
            Self::Pass | Self::Panicked(_) => &[],
        }
    }
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "non-string panic payload".to_string()
    }
}

fn evaluate_property<T>(input: &T, property: &impl Fn(&T) -> StructuredVerdict) -> PropertyOutcome {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| property(input))) {
        Ok(StructuredVerdict::Pass) => PropertyOutcome::Pass,
        Ok(StructuredVerdict::Fail(failure)) => PropertyOutcome::Failed(failure),
        Err(payload) => PropertyOutcome::Panicked(panic_message(payload.as_ref())),
    }
}

fn json_string(value: &str) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut escaped = String::with_capacity(value.len() + 2);
    escaped.push('"');
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let code = control as usize;
                escaped.push_str("\\u00");
                escaped.push(HEX[code >> 4] as char);
                escaped.push(HEX[code & 0x0f] as char);
            }
            other => escaped.push(other),
        }
    }
    escaped.push('"');
    escaped
}

fn json_string_object(fields: &[(&'static str, String)]) -> String {
    let mut object = String::from("{");
    for (index, (name, value)) in fields.iter().enumerate() {
        if index > 0 {
            object.push(',');
        }
        object.push_str(&json_string(name));
        object.push(':');
        object.push_str(&json_string(value));
    }
    object.push('}');
    object
}

fn parse_replay_value(value: &str) -> Result<u64, &'static str> {
    if value.is_empty() || !value.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err("expected a nonempty unsigned decimal case index");
    }
    value
        .parse::<u64>()
        .map_err(|_| "case index does not fit in u64")
}

fn replay_case_index() -> Option<u64> {
    let raw = std::env::var_os("FSIM_PROPCHECK_REPLAY")?;
    let value = raw
        .into_string()
        .unwrap_or_else(|_| panic!("FSIM_PROPCHECK_REPLAY must be valid UTF-8 unsigned decimal"));
    Some(
        parse_replay_value(&value)
            .unwrap_or_else(|problem| panic!("invalid FSIM_PROPCHECK_REPLAY={value:?}: {problem}")),
    )
}

fn replay_artifact_path_from(
    configured: Option<std::ffi::OsString>,
    cargo_target_dir: Option<std::ffi::OsString>,
    temp_dir: &std::path::Path,
    process_id: u32,
) -> std::path::PathBuf {
    if let Some(path) = configured {
        return path.into();
    }
    let directory =
        cargo_target_dir.map_or_else(|| temp_dir.to_path_buf(), std::path::PathBuf::from);
    directory.join(format!("fs-propcheck-replay-{process_id}.jsonl"))
}

fn replay_artifact_path() -> std::path::PathBuf {
    replay_artifact_path_from(
        std::env::var_os("FSIM_PROPCHECK_REPLAY_FILE"),
        std::env::var_os("CARGO_TARGET_DIR"),
        &std::env::temp_dir(),
        std::process::id(),
    )
}

fn write_replay_row(writer: &mut impl std::io::Write, row: &str) -> std::io::Result<()> {
    writer.write_all(row.as_bytes())?;
    writer.write_all(b"\n")?;
    writer.flush()
}

fn append_replay_artifact(path: &std::path::Path, row: &str) -> Result<(), String> {
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .map_err(|error| format!("open {}: {error}", path.display()))?;
    write_replay_row(&mut file, row).map_err(|error| format!("write {}: {error}", path.display()))
}

/// Run `cases` generated checks of `property`; on the first failure,
/// including a caught property panic, SHRINK to a local-minimum
/// counterexample, emit one valid dependency-free JSONL row to stdout and a
/// replay artifact, and panic with the replay seed. By default the artifact is
/// process-scoped under `CARGO_TARGET_DIR` (or the OS temp directory); set
/// `FSIM_PROPCHECK_REPLAY_FILE=<path>` to select an exact CI artifact path.
/// Set `FSIM_PROPCHECK_REPLAY=<case_seed>` to rerun exactly the failing case
/// (generation and shrink included). A present but malformed replay value is a
/// caller error and fails closed.
///
/// # Panics
/// On malformed replay configuration or the first (shrunk) property failure —
/// that is the test-runner API.
pub fn check<T: Shrink + core::fmt::Debug>(
    property_name: &str,
    suite_seed: u64,
    cases: u64,
    generate: impl Fn(&mut Stream) -> T,
    property: impl Fn(&T) -> bool,
) {
    check_structured(property_name, suite_seed, cases, generate, |input| {
        if property(input) {
            StructuredVerdict::Pass
        } else {
            StructuredVerdict::Fail(StructuredFailure::new(
                "returned-false",
                "property returned false",
                Vec::new(),
            ))
        }
    });
}

pub(crate) fn check_structured<T: Shrink + core::fmt::Debug>(
    property_name: &str,
    suite_seed: u64,
    cases: u64,
    generate: impl Fn(&mut Stream) -> T,
    property: impl Fn(&T) -> StructuredVerdict,
) {
    let run_case = |case_index: u64| {
        let mut stream = Stream::for_case(suite_seed, case_index);
        let input = generate(&mut stream);
        let initial_outcome = evaluate_property(&input, &property);
        if initial_outcome.passed() {
            return;
        }
        let initial_failure_kind = initial_outcome.failure_kind();
        let report = minimize(
            input.clone(),
            |candidate| {
                let outcome = evaluate_property(candidate, &property);
                outcome.passed() || outcome.failure_kind() != initial_failure_kind
            },
            10_000,
        )
        .unwrap_or_else(|error| match error {
            MinimizeError::SeedPasses => ShrinkReport {
                minimized: input,
                steps: 0,
                tried: 1,
                converged: false,
            },
            MinimizeError::EmptyEvaluationBudget => {
                unreachable!("the check runner's default budget evaluates the seed")
            }
        });
        let final_outcome = evaluate_property(&report.minimized, &property);
        let minimized = format!("{:?}", report.minimized);
        let counterexample_status = if report.converged {
            "local-fixpoint"
        } else {
            "budget-limited"
        };
        let replay_path = replay_artifact_path();
        let property_json = json_string(property_name);
        let minimized_json = json_string(&minimized);
        let detail_json = json_string(final_outcome.detail());
        let context_json = json_string_object(final_outcome.context());
        let replay_path_json = json_string(&replay_path.to_string_lossy());
        let failure_row = format!(
            "{{\"suite\":\"fs-propcheck\",\"property\":{property_json},\
             \"verdict\":\"fail\",\"suite_seed\":{suite_seed},\
             \"case_seed\":{case_index},\"shrink_steps\":{},\"shrink_tried\":{},\
             \"converged\":{},\"counterexample_status\":\"{counterexample_status}\",\
             \"failure_kind\":\"{}\",\"failure_detail\":{detail_json},\
             \"failure_context\":{context_json},\"minimized\":{minimized_json},\
             \"replay_file\":{replay_path_json}}}",
            report.steps,
            report.tried,
            report.converged,
            final_outcome.failure_kind(),
        );
        let artifact_result = append_replay_artifact(&replay_path, &failure_row);
        ::std::println!("{failure_row}");
        let counterexample_label = if report.converged {
            "local-minimum counterexample"
        } else {
            "best-known counterexample"
        };
        match artifact_result {
            Ok(()) => panic!(
                "property `{property_name}` failed ({}: {}); {counterexample_label}: {:?} \
                 (replay: FSIM_PROPCHECK_REPLAY={case_index}, suite seed {suite_seed}, \
                 replay artifact: {}, {} shrink steps)",
                final_outcome.failure_kind(),
                final_outcome.detail(),
                report.minimized,
                replay_path.display(),
                report.steps
            ),
            Err(artifact_error) => panic!(
                "property `{property_name}` failed ({}: {}); {counterexample_label}: {:?} \
                 (replay: FSIM_PROPCHECK_REPLAY={case_index}, suite seed {suite_seed}, \
                 replay artifact emission failed: {artifact_error}, {} shrink steps)",
                final_outcome.failure_kind(),
                final_outcome.detail(),
                report.minimized,
                report.steps
            ),
        }
    };

    if let Some(case_index) = replay_case_index() {
        run_case(case_index);
    } else {
        for case_index in 0..cases {
            run_case(case_index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        StructuredFailure, StructuredVerdict, check_structured, json_string, json_string_object,
        parse_replay_value, replay_artifact_path_from, write_replay_row,
    };

    #[test]
    fn json_string_escapes_controls_quotes_and_backslashes() {
        assert_eq!(
            json_string("quote=\" slash=\\ newline=\n nul=\0 snowman=☃"),
            "\"quote=\\\" slash=\\\\ newline=\\n nul=\\u0000 snowman=☃\""
        );
        assert_eq!(
            json_string_object(&[("relation", "rigid\nmove".to_string())]),
            "{\"relation\":\"rigid\\nmove\"}"
        );
    }

    #[test]
    fn replay_parser_refuses_malformed_or_overflowing_values() {
        assert_eq!(parse_replay_value("42"), Ok(42));
        assert!(parse_replay_value("").is_err());
        assert!(parse_replay_value(" 42").is_err());
        assert!(parse_replay_value("+42").is_err());
        assert!(parse_replay_value("18446744073709551616").is_err());
    }

    #[test]
    fn structured_shrinking_preserves_the_initial_failure_kind() {
        let result = std::panic::catch_unwind(|| {
            check_structured(
                "failure-kind-stability",
                7,
                1,
                |_| 9_i64,
                |value| {
                    if *value >= 3 {
                        StructuredVerdict::Fail(StructuredFailure::new(
                            "primary-kind",
                            "primary failure",
                            Vec::new(),
                        ))
                    } else if *value == 0 {
                        StructuredVerdict::Fail(StructuredFailure::new(
                            "incidental-kind",
                            "incidental failure",
                            Vec::new(),
                        ))
                    } else {
                        StructuredVerdict::Pass
                    }
                },
            );
        });
        let payload = result.expect_err("the primary failure remains");
        let message = payload
            .downcast_ref::<String>()
            .expect("string panic from runner");
        assert!(message.contains("primary-kind"), "{message}");
        assert!(
            message.contains("local-minimum counterexample: 3"),
            "{message}"
        );
    }

    #[test]
    fn replay_artifact_path_and_jsonl_writer_are_deterministic() {
        let configured = replay_artifact_path_from(
            Some(std::ffi::OsString::from("/ci/replay.jsonl")),
            Some(std::ffi::OsString::from("/target")),
            std::path::Path::new("/tmp"),
            42,
        );
        assert_eq!(configured, std::path::Path::new("/ci/replay.jsonl"));

        let defaulted = replay_artifact_path_from(
            None,
            Some(std::ffi::OsString::from("/target")),
            std::path::Path::new("/tmp"),
            42,
        );
        assert_eq!(
            defaulted,
            std::path::Path::new("/target/fs-propcheck-replay-42.jsonl")
        );

        let mut bytes = Vec::new();
        write_replay_row(&mut bytes, "{\"case_seed\":7}").expect("in-memory write");
        assert_eq!(bytes, b"{\"case_seed\":7}\n");
    }
}
