//! compound — the Gauntlet failure-compounding workflow (bead 6nb.9).
//!
//! Every golden break, falsifier hit, guard failure, or property
//! counterexample should STRENGTHEN the permanent test surface instead of
//! being fixed and forgotten. This module is the v1 mechanism:
//!
//! 1. **Capture** the failure as a [`FailureCase`]: seed, typed input, the
//!    violated [`InvariantClass`], and the contract surface it broke.
//! 2. **Minimize** it ([`minimize`]): deterministic greedy descent through
//!    [`Shrink`] candidates, keeping the invariant violated at every step.
//!    A non-failing input is a typed refusal, never a silent no-op.
//! 3. **Probe the neighborhood** ([`probe_neighborhood`]): bounded, labeled
//!    perturbations around the minimum expose whether the failure is a
//!    point or a region.
//! 4. **Land a family** ([`RegressionFamily`]): the minimum plus its failing
//!    neighbors, with tracking-issue references and a recommended admission
//!    rule when the class is general.
//! 5. **Replay** ([`RegressionFamily::replay`]): the family is
//!    content-addressed ([`Canon`] bytes → FNV-64) and re-executable — a
//!    member that stops failing is REPORTED, because a regression family
//!    whose members silently pass is stale evidence.
//!
//! Everything is plain data and deterministic: same case + same predicate
//! ⇒ bitwise-identical minimum, probes, manifest, and content hash, on
//! every ISA and in every build mode (the canon encoding is integer bytes
//! and `f64::to_bits`, never formatted floats).
//!
//! What this module does NOT do (no-claims): it does not write to the
//! ledger or emit fs-obs events (recorded follow-up once the huq.16 schema
//! lands), and it does not itself change admission rules — the family
//! CARRIES the recommendation (as check-powi was born from the powi
//! incident); enacting it is the responding agent's task.

/// Semantic version of the canon encoding + content-hash assembly
/// (golden-couplings surface `fs-bisect:compound-canon`). Changing the
/// [`Canon`] byte layout, the tag values, the FNV constants, or the
/// field order in [`RegressionFamily::content_hash`] changes every
/// family hash — bump this and deliberately re-freeze the dependents
/// listed in golden-couplings.json (docs/GOLDEN_POLICY.md).
pub const COMPOUND_CANON_VERSION: u32 = 1;

/// The invariant a failure violated — the classification axis that decides
/// which sibling surfaces the lesson propagates to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InvariantClass {
    /// Bits differ between debug and release builds (e.g. the `f64::powi`
    /// incident, bead 4xnt).
    BuildModeDeterminism,
    /// Bits differ across ISAs where the contract claims they must not.
    CrossIsaDeterminism,
    /// A frozen golden hash no longer matches the observed bits.
    GoldenDrift,
    /// A certified enclosure excludes the true value.
    EnclosureViolation,
    /// A certificate accepted something its falsifier refutes.
    CertificateForgery,
    /// A conserved quantity drifted beyond its stated band.
    ConservationViolation,
    /// A gradient/adjoint disagrees with its independent check.
    AdjointInconsistency,
    /// Anything else — named, never silent.
    Other(String),
}

impl InvariantClass {
    /// Stable name for manifests and hashes.
    #[must_use]
    pub fn name(&self) -> &str {
        match self {
            InvariantClass::BuildModeDeterminism => "build-mode-determinism",
            InvariantClass::CrossIsaDeterminism => "cross-isa-determinism",
            InvariantClass::GoldenDrift => "golden-drift",
            InvariantClass::EnclosureViolation => "enclosure-violation",
            InvariantClass::CertificateForgery => "certificate-forgery",
            InvariantClass::ConservationViolation => "conservation-violation",
            InvariantClass::AdjointInconsistency => "adjoint-inconsistency",
            InvariantClass::Other(s) => s,
        }
    }
}

/// A captured failure: everything needed to reproduce it deterministically.
#[derive(Debug, Clone)]
pub struct FailureCase<I> {
    /// Stable identifier (used in manifests and issue references).
    pub id: String,
    /// The seed that produced the input (0 when the input is explicit).
    pub seed: u64,
    /// The failing input itself.
    pub input: I,
    /// Which invariant broke.
    pub invariant: InvariantClass,
    /// The contract surface that broke, as `crate::surface` prose.
    pub contract: String,
    /// One-line human detail (observed vs expected).
    pub detail: String,
}

/// Typed refusals — the workflow never silently does nothing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompoundError {
    /// The captured input does not fail the predicate: there is nothing to
    /// minimize, and pretending otherwise would produce a fake regression.
    NotFailing {
        /// The case id that was expected to fail.
        id: String,
    },
}

/// Deterministic shrinking: candidates strictly "smaller" than `self`, in a
/// FIXED order. An empty vector means fully shrunk.
pub trait Shrink: Clone {
    /// Smaller candidate inputs, most aggressive first (convention).
    fn shrink_candidates(&self) -> Vec<Self>;
}

/// The result of [`minimize`].
#[derive(Debug, Clone)]
pub struct MinimizeReport<I> {
    /// The smallest input found that still fails.
    pub minimized: I,
    /// Accepted shrink steps.
    pub steps: usize,
    /// Total candidates evaluated.
    pub tried: usize,
    /// False when the step budget ran out before a fixpoint — the minimum
    /// is honest but possibly not minimal.
    pub converged: bool,
}

/// Greedy deterministic minimization: repeatedly take the FIRST failing
/// shrink candidate until none fails (fixpoint) or `max_steps` accepted
/// steps. Same input + same predicate ⇒ identical trajectory.
///
/// # Errors
/// [`CompoundError::NotFailing`] when `input` does not fail `fails`.
pub fn minimize<I: Shrink>(
    id: &str,
    input: &I,
    fails: &dyn Fn(&I) -> bool,
    max_steps: usize,
) -> Result<MinimizeReport<I>, CompoundError> {
    if !fails(input) {
        return Err(CompoundError::NotFailing { id: id.to_string() });
    }
    let mut current = input.clone();
    let mut steps = 0usize;
    let mut tried = 0usize;
    let mut converged = false;
    'outer: loop {
        if steps == max_steps {
            break;
        }
        for cand in current.shrink_candidates() {
            tried += 1;
            if fails(&cand) {
                current = cand;
                steps += 1;
                continue 'outer;
            }
        }
        converged = true;
        break;
    }
    Ok(MinimizeReport {
        minimized: current,
        steps,
        tried,
        converged,
    })
}

/// One labeled neighborhood probe outcome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NeighborProbe {
    /// The caller's label for this neighbor (e.g. `"k=5"`).
    pub label: String,
    /// Whether the invariant is violated there too.
    pub fails: bool,
}

/// The bounded neighborhood around a minimized failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NeighborhoodReport {
    /// Every probe, in the caller's (deterministic) order.
    pub probes: Vec<NeighborProbe>,
    /// How many neighbors also fail (region vs point evidence).
    pub failing: usize,
}

/// Evaluate a bounded, labeled set of neighbors of a minimized failure.
/// The caller supplies the neighbors, so the probe is bounded by
/// construction and its order is the caller's order.
pub fn probe_neighborhood<I>(
    neighbors: &[(String, I)],
    fails: &dyn Fn(&I) -> bool,
) -> NeighborhoodReport {
    let probes: Vec<NeighborProbe> = neighbors
        .iter()
        .map(|(label, input)| NeighborProbe {
            label: label.clone(),
            fails: fails(input),
        })
        .collect();
    let failing = probes.iter().filter(|p| p.fails).count();
    NeighborhoodReport { probes, failing }
}

/// Canonical bytes for content addressing. Tagged and length-prefixed so
/// distinct structures cannot collide by concatenation; floats canonicalize
/// through `to_bits`, never through formatting.
pub trait Canon {
    /// Append this value's canonical bytes.
    fn canon(&self, out: &mut Vec<u8>);
}

impl Canon for u64 {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(1);
        out.extend_from_slice(&self.to_le_bytes());
    }
}
impl Canon for i64 {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(2);
        out.extend_from_slice(&self.to_le_bytes());
    }
}
impl Canon for i32 {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(3);
        out.extend_from_slice(&self.to_le_bytes());
    }
}
impl Canon for f64 {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(4);
        out.extend_from_slice(&self.to_bits().to_le_bytes());
    }
}
impl Canon for bool {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(5);
        out.push(u8::from(*self));
    }
}
impl Canon for str {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(6);
        out.extend_from_slice(&(self.len() as u64).to_le_bytes());
        out.extend_from_slice(self.as_bytes());
    }
}
impl Canon for String {
    fn canon(&self, out: &mut Vec<u8>) {
        self.as_str().canon(out);
    }
}
impl<T: Canon> Canon for Vec<T> {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(7);
        out.extend_from_slice(&(self.len() as u64).to_le_bytes());
        for item in self {
            item.canon(out);
        }
    }
}
impl<A: Canon, B: Canon> Canon for (A, B) {
    fn canon(&self, out: &mut Vec<u8>) {
        out.push(8);
        self.0.canon(out);
        self.1.canon(out);
    }
}

/// FNV-1a 64 over canonical bytes — the house content-hash idiom.
#[must_use]
pub fn fnv64(bytes: &[u8]) -> u64 {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        acc ^= u64::from(b);
        acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
    acc
}

/// A permanent regression family: the minimized case plus its failing
/// neighbors, with provenance. This is the artifact a failure leaves
/// behind — MORE than one example, linked to its tracking issue, carrying
/// the admission-rule lesson when one generalizes.
#[derive(Debug, Clone)]
pub struct RegressionFamily<I> {
    /// Family name (stable, kebab-case).
    pub name: String,
    /// The invariant every member violates.
    pub invariant: InvariantClass,
    /// Labeled members; `members[0]` is the minimized case by convention.
    pub members: Vec<(String, I)>,
    /// Tracking references (bead ids / issue ids) — never empty for a
    /// landed family; a failure without a paper trail cannot compound.
    pub tracking: Vec<String>,
    /// The generalized lesson, when there is one (e.g. "lint variable-
    /// exponent powi out of deterministic paths").
    pub recommended_admission: Option<String>,
}

/// A replayed family: which members still fail (live) and which now pass
/// (stale — the regression they pinned was fixed or the predicate moved).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplayReport {
    /// Labels of members that still violate the invariant.
    pub still_failing: Vec<String>,
    /// Labels of members that no longer violate it.
    pub now_passing: Vec<String>,
}

impl<I: Canon> RegressionFamily<I> {
    /// Content hash over the full canonical encoding: name, invariant,
    /// members (labels + inputs), tracking, admission. Deterministic across
    /// runs, build modes, and ISAs.
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        let mut bytes = Vec::new();
        self.name.canon(&mut bytes);
        self.invariant.name().canon(&mut bytes);
        bytes.push(7);
        bytes.extend_from_slice(&(self.members.len() as u64).to_le_bytes());
        for (label, input) in &self.members {
            label.canon(&mut bytes);
            input.canon(&mut bytes);
        }
        self.tracking.canon(&mut bytes);
        match &self.recommended_admission {
            Some(a) => {
                bytes.push(9);
                a.canon(&mut bytes);
            }
            None => bytes.push(10),
        }
        fnv64(&bytes)
    }

    /// The replayable manifest: JSON-lines, one header, one line per member
    /// (canonical bytes hex-encoded), one trailer with the content hash.
    #[must_use]
    pub fn manifest(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        let _ = writeln!(
            out,
            "{{\"family\":\"{}\",\"invariant\":\"{}\",\"members\":{},\"tracking\":\"{}\"}}",
            self.name,
            self.invariant.name(),
            self.members.len(),
            self.tracking.join(",")
        );
        for (label, input) in &self.members {
            let mut bytes = Vec::new();
            input.canon(&mut bytes);
            let hex: String = bytes.iter().fold(String::new(), |mut s, b| {
                let _ = write!(s, "{b:02x}");
                s
            });
            let _ = writeln!(out, "{{\"member\":\"{label}\",\"canon\":\"{hex}\"}}");
        }
        let _ = writeln!(
            out,
            "{{\"content_hash\":\"{:#018x}\"}}",
            self.content_hash()
        );
        out
    }

    /// Re-execute every member against the predicate. A live family has
    /// `now_passing` empty; anything else is stale evidence to act on.
    #[must_use]
    pub fn replay(&self, fails: &dyn Fn(&I) -> bool) -> ReplayReport {
        let mut still_failing = Vec::new();
        let mut now_passing = Vec::new();
        for (label, input) in &self.members {
            if fails(input) {
                still_failing.push(label.clone());
            } else {
                now_passing.push(label.clone());
            }
        }
        ReplayReport {
            still_failing,
            now_passing,
        }
    }
}

/// The full workflow output: minimized case, neighborhood, family, hash.
#[derive(Debug, Clone)]
pub struct CompoundReport<I> {
    /// The captured case with its input replaced by the minimum.
    pub case: FailureCase<I>,
    /// Minimization statistics.
    pub steps: usize,
    /// Whether minimization reached a fixpoint.
    pub converged: bool,
    /// The bounded neighborhood around the minimum.
    pub neighborhood: NeighborhoodReport,
    /// The landed family (minimum first, then failing neighbors).
    pub family: RegressionFamily<I>,
    /// The family's content hash (also in the manifest trailer).
    pub content_hash: u64,
}

/// The v1 workflow driver: minimize → probe → land the family.
///
/// `neighbors_of` receives the MINIMIZED input and must return a bounded,
/// deterministically ordered, labeled neighbor set. Failing neighbors join
/// the family behind the minimum.
///
/// # Errors
/// [`CompoundError::NotFailing`] when the captured input does not fail.
pub fn compound<I: Shrink + Canon>(
    case: FailureCase<I>,
    fails: &dyn Fn(&I) -> bool,
    neighbors_of: &dyn Fn(&I) -> Vec<(String, I)>,
    tracking: Vec<String>,
    recommended_admission: Option<String>,
    max_steps: usize,
) -> Result<CompoundReport<I>, CompoundError> {
    let report = minimize(&case.id, &case.input, fails, max_steps)?;
    let neighbors = neighbors_of(&report.minimized);
    let neighborhood = probe_neighborhood(&neighbors, fails);
    let mut members: Vec<(String, I)> = vec![("minimized".to_string(), report.minimized.clone())];
    for ((label, input), probe) in neighbors.into_iter().zip(&neighborhood.probes) {
        if probe.fails {
            members.push((label, input));
        }
    }
    let family = RegressionFamily {
        name: case.id.clone(),
        invariant: case.invariant.clone(),
        members,
        tracking,
        recommended_admission,
    };
    let content_hash = family.content_hash();
    Ok(CompoundReport {
        case: FailureCase {
            input: report.minimized,
            ..case
        },
        steps: report.steps,
        converged: report.converged,
        neighborhood,
        family,
        content_hash,
    })
}
