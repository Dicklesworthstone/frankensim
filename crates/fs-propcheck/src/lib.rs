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
//! - **Shrinking**: on failure the input is minimized by greedy
//!   first-failing-candidate descent (the `fs-bisect` convention:
//!   candidates most-aggressive-first, fixed order, bounded steps), so
//!   CI failures arrive as minimal counterexamples, not haystacks.
//! - **Zero runtime dependencies** (UTIL layer): the generator stream is
//!   a test-input source, not a simulation RNG — kernels keep using
//!   fs-rand's Philox discipline; this stream never touches ledgers.
//! - **JSONL failure rows**: every failure emits one structured line
//!   (property, case seed, shrink steps, minimized debug form) before
//!   panicking, per the house verdict-logging style.

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

    /// Uniform in `[lo, hi]` (inclusive); `lo > hi` is a caller bug and
    /// panics with the offending bounds (test-harness semantics).
    pub fn int_in(&mut self, lo: i64, hi: i64) -> i64 {
        assert!(lo <= hi, "int_in bounds inverted: [{lo}, {hi}]");
        let span = hi.wrapping_sub(lo).cast_unsigned();
        if span == u64::MAX {
            return self.next_u64().cast_signed();
        }
        lo.wrapping_add((self.next_u64() % (span + 1)).cast_signed())
    }

    /// Uniform-ish f64 in `[lo, hi)`, with a 1-in-8 chance of drawing a
    /// SPECIAL value (0.0, -0.0, lo, hi, 1.0, -1.0) — the corners where
    /// numeric laws actually break.
    pub fn f64_in(&mut self, lo: f64, hi: f64) -> f64 {
        assert!(lo < hi, "f64_in bounds inverted: [{lo}, {hi})");
        if self.next_u64().is_multiple_of(8) {
            let specials = [0.0, -0.0, lo, hi, 1.0, -1.0];
            let pick = specials[(self.next_u64() % specials.len() as u64) as usize];
            if pick >= lo && pick < hi {
                return pick;
            }
        }
        let unit = (self.next_u64() >> 11) as f64 / f64::from(1u32 << 26) / f64::from(1u32 << 27);
        lo + (hi - lo) * unit
    }

    /// A vector of `n` in `[0, max_len]` elements drawn by `elem`.
    pub fn vec_of<T>(&mut self, max_len: usize, mut elem: impl FnMut(&mut Stream) -> T) -> Vec<T> {
        let n = (self.next_u64() % (max_len as u64 + 1)) as usize;
        (0..n).map(|_| elem(self)).collect()
    }
}

/// Deterministic shrinking: candidates strictly "smaller" than `self`,
/// most aggressive first, in a FIXED order. Empty = fully shrunk.
///
/// This is the `fs-bisect::compound` convention verbatim, so the two
/// engines stay interchangeable; `fs-bisect` adoption of this shared
/// trait is a follow-up slice on the bead.
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
        out.retain(|c| c.abs() < x.abs() || (c.abs() == x.abs() && *c > x));
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
        if x == 0.0 || !x.is_finite() {
            return vec![];
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
        out.push(self[..self.len() / 2].to_vec());
        out.push(self[self.len() / 2..].to_vec());
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
        let mut out: Vec<(A, B)> = self
            .0
            .shrink_candidates()
            .into_iter()
            .map(|a| (a, self.1.clone()))
            .collect();
        out.extend(
            self.1
                .shrink_candidates()
                .into_iter()
                .map(|b| (self.0.clone(), b)),
        );
        out
    }
}

impl<A: Shrink, B: Shrink, C: Shrink> Shrink for (A, B, C) {
    fn shrink_candidates(&self) -> Vec<(A, B, C)> {
        let mut out: Vec<(A, B, C)> = self
            .0
            .shrink_candidates()
            .into_iter()
            .map(|a| (a, self.1.clone(), self.2.clone()))
            .collect();
        out.extend(
            self.1
                .shrink_candidates()
                .into_iter()
                .map(|b| (self.0.clone(), b, self.2.clone())),
        );
        out.extend(
            self.2
                .shrink_candidates()
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
    /// Total property evaluations during shrinking.
    pub tried: usize,
    /// False when the step budget ran out before a fixpoint.
    pub converged: bool,
}

/// Greedy deterministic minimization: repeatedly take the FIRST failing
/// candidate until fixpoint or `max_steps`. Same input + property gives
/// the identical trajectory.
pub fn minimize<T: Shrink>(
    seed_failure: T,
    property: impl Fn(&T) -> bool,
    max_steps: usize,
) -> ShrinkReport<T> {
    let mut current = seed_failure;
    let mut steps = 0;
    let mut tried = 0;
    'outer: while steps < max_steps {
        for cand in current.shrink_candidates() {
            tried += 1;
            if !property(&cand) {
                current = cand;
                steps += 1;
                continue 'outer;
            }
        }
        return ShrinkReport {
            minimized: current,
            steps,
            tried,
            converged: true,
        };
    }
    ShrinkReport {
        minimized: current,
        steps,
        tried,
        converged: false,
    }
}

/// Run `cases` generated checks of `property`; on the first failure,
/// SHRINK to a minimal counterexample, emit one JSONL row, and panic
/// with the replay seed. Set `FSIM_PROPCHECK_REPLAY=<case_seed>` to
/// rerun exactly the failing case (generation and shrink included).
///
/// # Panics
/// On the first (shrunk) property failure — that is the point.
pub fn check<T: Shrink + core::fmt::Debug>(
    property_name: &str,
    suite_seed: u64,
    cases: u64,
    generate: impl Fn(&mut Stream) -> T,
    property: impl Fn(&T) -> bool,
) {
    let replay = std::env::var("FSIM_PROPCHECK_REPLAY")
        .ok()
        .and_then(|v| v.parse::<u64>().ok());
    let indices: Vec<u64> = match replay {
        Some(case_seed) => vec![case_seed],
        None => (0..cases).collect(),
    };
    for case_index in indices {
        let mut stream = Stream::for_case(suite_seed, case_index);
        let input = generate(&mut stream);
        if property(&input) {
            continue;
        }
        let report = minimize(input, &property, 10_000);
        println!(
            "{{\"suite\":\"fs-propcheck\",\"property\":\"{property_name}\",\
             \"verdict\":\"fail\",\"suite_seed\":{suite_seed},\
             \"case_seed\":{case_index},\"shrink_steps\":{},\"shrink_tried\":{},\
             \"converged\":{},\"minimized\":\"{:?}\"}}",
            report.steps,
            report.tried,
            report.converged,
            report.minimized
        );
        panic!(
            "property `{property_name}` failed; minimal counterexample: {:?} \
             (replay: FSIM_PROPCHECK_REPLAY={case_index}, suite seed {suite_seed}, \
             {} shrink steps)",
            report.minimized, report.steps
        );
    }
}
