//! compound — the Gauntlet failure-compounding workflow (bead 6nb.9).
//!
//! Every golden break, falsifier hit, guard failure, or property
//! counterexample should STRENGTHEN the permanent test surface instead of
//! being fixed and forgotten. This module is the v3 mechanism:
//!
//! 1. **Capture** the failure as a [`FailureCase`]: seed, typed input, the
//!    violated [`InvariantClass`], and the contract surface it broke.
//! 2. **Minimize** it ([`minimize`]): deterministic greedy descent through
//!    [`Shrink`] candidates, keeping the invariant violated at every step.
//!    A non-failing input is a typed refusal, never a silent no-op.
//! 3. **Probe the neighborhood** ([`probe_neighborhood`]): bounded, uniquely
//!    labeled perturbations around the minimum expose whether the failure is
//!    a point or a region.
//! 4. **Land a family** ([`RegressionFamily`]): the minimum plus its failing
//!    neighbors, with tracking-issue references and a recommended admission
//!    rule when the class is general.
//! 5. **Replay** ([`RegressionFamily::replay`]): the family is
//!    content-addressed (stable [`Canon`] codec/schema + bounded bytes →
//!    domain-separated BLAKE3) and re-executable. Live typed values must match
//!    the sealed snapshots before any predicate work; an authentic member that
//!    stops failing is REPORTED because silently passing members are stale
//!    evidence.
//!
//! Everything is plain data and deterministic: same case + same predicate +
//! stable member codec ⇒ bitwise-identical minimum, probes, manifest, and
//! content hash. Cross-ISA equality is a member-codec claim: the built-in
//! integer and `f64::to_bits` codecs make it, while caller-defined codecs must
//! state and test their own portability boundary.
//!
//! What this module does NOT do (no-claims): it does not write to the
//! ledger or emit fs-obs events (recorded follow-up once the huq.16 schema
//! lands), and it does not itself change admission rules — the family
//! CARRIES the recommendation (as check-powi was born from the powi
//! incident); enacting it is the responding agent's task.

/// Semantic version of the canon encoding + content-hash assembly
/// (golden-couplings surface `fs-bisect:compound-canon`). Changing the
/// [`Canon`] byte layout, the tag values, the hash domain, or the
/// field order in [`RegressionFamily::content_hash`] changes every
/// family hash — bump this and deliberately re-freeze the dependents
/// listed in golden-couplings.json (docs/GOLDEN_POLICY.md).
pub const COMPOUND_CANON_VERSION: u32 = 3;

/// Domain separating regression-family identities from every other BLAKE3 use.
pub const COMPOUND_FAMILY_HASH_DOMAIN: &str = "org.frankensim.fs-bisect.compound-family.v3";

/// Maximum accepted minimizer descent steps.
pub const MAX_MINIMIZE_STEPS: usize = 65_536;
/// Maximum shrink candidates returned for one descent step.
pub const MAX_SHRINK_CANDIDATES_PER_STEP: usize = 4_096;
/// Maximum predicate evaluations across one minimization.
pub const MAX_MINIMIZE_EVALUATIONS: usize = 1_000_000;
/// Maximum neighboring inputs evaluated and retained for one family.
pub const MAX_NEIGHBOR_PROBES: usize = 4_096;
/// Maximum tracking references attached to one regression family.
pub const MAX_TRACKING_REFS: usize = 64;
/// Maximum bytes in a case/family/member/tracking identifier.
pub const MAX_IDENTIFIER_BYTES: usize = 256;
/// Maximum bytes in a contract or admission-rule description.
pub const MAX_DESCRIPTION_BYTES: usize = 16 * 1024;
/// Maximum canonical payload bytes retained for one regression member.
pub const MAX_CANONICAL_MEMBER_BYTES: usize = 1024 * 1024;
/// Maximum canonical payload bytes retained across one family.
pub const MAX_CANONICAL_FAMILY_BYTES: usize = 16 * 1024 * 1024;
/// Maximum bytes in the stable member codec/schema descriptor.
pub const MAX_CANONICAL_SCHEMA_BYTES: usize = 16 * 1024;

const RESERVED_INVARIANT_NAMES: [&str; 7] = [
    "build-mode-determinism",
    "cross-isa-determinism",
    "golden-drift",
    "enclosure-violation",
    "certificate-forgery",
    "conservation-violation",
    "adjoint-inconsistency",
];

fn visible_identifier(value: &str) -> bool {
    value.bytes().all(|byte| (b'!'..=b'~').contains(&byte))
}

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

    fn validate(&self) -> Result<(), CompoundError> {
        let Self::Other(name) = self else {
            return Ok(());
        };
        validate_identifier("invariant", name)?;
        if RESERVED_INVARIANT_NAMES.contains(&name.as_str()) {
            return Err(CompoundError::InvalidField {
                field: "invariant",
                problem: format!(
                    "custom invariant name {name:?} is reserved by a built-in variant"
                ),
            });
        }
        Ok(())
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

/// Immutable provenance bound into a regression family's identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FamilyProvenance {
    seed: u64,
    contract: String,
    detail: String,
}

impl FamilyProvenance {
    /// Construct bounded provenance for one captured failure.
    ///
    /// # Errors
    /// Empty or oversized contract/detail descriptions are refused.
    pub fn new(seed: u64, contract: String, detail: String) -> Result<Self, CompoundError> {
        validate_description("contract", &contract)?;
        validate_description("detail", &detail)?;
        Ok(Self {
            seed,
            contract,
            detail,
        })
    }

    /// Reproduction seed (`0` when the member was explicit).
    #[must_use]
    pub const fn seed(&self) -> u64 {
        self.seed
    }

    /// Contract surface that was violated.
    #[must_use]
    pub fn contract(&self) -> &str {
        &self.contract
    }

    /// Captured expected/observed diagnosis.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
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
    /// A caller-controlled field violates the canonical family envelope.
    InvalidField {
        /// Stable field name.
        field: &'static str,
        /// Actionable diagnosis.
        problem: String,
    },
    /// A deterministic work or collection bound was exceeded.
    LimitExceeded {
        /// Bounded resource.
        resource: &'static str,
        /// Requested or observed value.
        requested: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// Two labels or tracking references would make the manifest ambiguous.
    DuplicateIdentity {
        /// Collection containing the duplicate.
        field: &'static str,
        /// Repeated value.
        value: String,
    },
    /// Permanent families require a completed minimization, not merely the
    /// best witness found before a caller budget expired.
    MinimizationIncomplete {
        /// The case whose minimization did not reach a fixpoint.
        id: String,
        /// Accepted shrink steps before exhaustion.
        steps: usize,
        /// Predicate evaluations, including the seed input.
        evaluations: usize,
    },
    /// A live typed member or its codec schema no longer matches the identity
    /// sealed at family construction, so replay semantics are unauthenticated.
    ReplayIdentityDrift {
        /// Member label, or `None` when the codec/schema itself drifted.
        member: Option<String>,
    },
    /// An evidence-producing callback mutated the canonical semantics of the
    /// value it inspected.
    CallbackIdentityDrift {
        /// Callback phase (`shrink_candidates`, `minimize`, `neighborhood`, or
        /// `neighbors_of`).
        phase: &'static str,
        /// Stable case/member identity.
        identity: String,
    },
}

impl core::fmt::Display for CompoundError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NotFailing { id } => {
                write!(f, "captured case {id:?} does not fail its predicate")
            }
            Self::InvalidField { field, problem } => {
                write!(f, "invalid failure-family field {field}: {problem}")
            }
            Self::LimitExceeded {
                resource,
                requested,
                max,
            } => write!(
                f,
                "failure compounding {resource} request {requested} exceeds limit {max}"
            ),
            Self::DuplicateIdentity { field, value } => {
                write!(f, "duplicate {field} identity {value:?}")
            }
            Self::MinimizationIncomplete {
                id,
                steps,
                evaluations,
            } => write!(
                f,
                "failure family {id:?} did not reach a minimization fixpoint \
                 ({steps} accepted steps, {evaluations} predicate evaluations)"
            ),
            Self::ReplayIdentityDrift {
                member: Some(member),
            } => write!(
                f,
                "failure-family member {member:?} changed after its identity was sealed"
            ),
            Self::ReplayIdentityDrift { member: None } => write!(
                f,
                "failure-family member codec/schema changed after its identity was sealed"
            ),
            Self::CallbackIdentityDrift { phase, identity } => write!(
                f,
                "failure-compounding callback {phase} mutated canonical identity {identity:?}"
            ),
        }
    }
}

impl std::error::Error for CompoundError {}

fn validate_identifier(field: &'static str, value: &str) -> Result<(), CompoundError> {
    if value.is_empty() {
        return Err(CompoundError::InvalidField {
            field,
            problem: "must not be empty".to_string(),
        });
    }
    if value.len() > MAX_IDENTIFIER_BYTES {
        return Err(CompoundError::LimitExceeded {
            resource: field,
            requested: value.len(),
            max: MAX_IDENTIFIER_BYTES,
        });
    }
    if !visible_identifier(value) {
        return Err(CompoundError::InvalidField {
            field,
            problem: "must contain visible ASCII bytes only".to_string(),
        });
    }
    Ok(())
}

fn validate_family_name(value: &str) -> Result<(), CompoundError> {
    validate_identifier("case_id", value)?;
    let mut bytes = value.bytes();
    let first = bytes.next().expect("validated non-empty");
    let last = value.as_bytes()[value.len() - 1];
    let alphanumeric = |byte: u8| byte.is_ascii_lowercase() || byte.is_ascii_digit();
    if !alphanumeric(first)
        || !alphanumeric(last)
        || !value.bytes().all(|byte| alphanumeric(byte) || byte == b'-')
        || value.as_bytes().windows(2).any(|pair| pair == b"--")
    {
        return Err(CompoundError::InvalidField {
            field: "case_id",
            problem: "must be lowercase kebab-case (ASCII letters/digits separated by '-')"
                .to_string(),
        });
    }
    Ok(())
}

fn validate_description(field: &'static str, value: &str) -> Result<(), CompoundError> {
    if value.trim().is_empty() {
        return Err(CompoundError::InvalidField {
            field,
            problem: "must not be empty or whitespace-only".to_string(),
        });
    }
    if value.len() > MAX_DESCRIPTION_BYTES {
        return Err(CompoundError::LimitExceeded {
            resource: field,
            requested: value.len(),
            max: MAX_DESCRIPTION_BYTES,
        });
    }
    Ok(())
}

fn validate_tracking_refs(tracking: &[String]) -> Result<(), CompoundError> {
    if tracking.is_empty() {
        return Err(CompoundError::InvalidField {
            field: "tracking",
            problem: "at least one Beads or issue reference is required".to_string(),
        });
    }
    if tracking.len() > MAX_TRACKING_REFS {
        return Err(CompoundError::LimitExceeded {
            resource: "tracking_refs",
            requested: tracking.len(),
            max: MAX_TRACKING_REFS,
        });
    }
    let mut tracking_refs = std::collections::BTreeSet::new();
    for reference in tracking {
        validate_identifier("tracking_ref", reference)?;
        if !tracking_refs.insert(reference.as_str()) {
            return Err(CompoundError::DuplicateIdentity {
                field: "tracking_ref",
                value: reference.clone(),
            });
        }
    }
    Ok(())
}

fn validate_admission_rule(rule: Option<&str>) -> Result<(), CompoundError> {
    if let Some(rule) = rule {
        validate_description("recommended_admission", rule)?;
    }
    Ok(())
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
    /// Total predicate evaluations, including the captured seed input.
    pub tried: usize,
    /// False when the step budget ran out before a fixpoint — the minimum
    /// is honest but possibly not minimal.
    pub converged: bool,
}

fn evaluate_stable<I: Canon>(
    phase: &'static str,
    identity: &str,
    input: &I,
    predicate: &dyn Fn(&I) -> bool,
) -> Result<bool, CompoundError> {
    let before = canonical_bytes(input)?;
    let result = predicate(input);
    if canonical_bytes(input)? != before {
        return Err(CompoundError::CallbackIdentityDrift {
            phase,
            identity: identity.to_string(),
        });
    }
    Ok(result)
}

/// Greedy deterministic minimization: repeatedly take the FIRST failing
/// shrink candidate until none fails (fixpoint) or `max_steps` accepted
/// steps. Same input + same predicate ⇒ identical trajectory.
///
/// # Errors
/// [`CompoundError::NotFailing`] when `input` does not fail `fails`.
pub fn minimize<I: Shrink + Canon>(
    id: &str,
    input: &I,
    fails: &dyn Fn(&I) -> bool,
    max_steps: usize,
) -> Result<MinimizeReport<I>, CompoundError> {
    validate_family_name(id)?;
    if max_steps > MAX_MINIMIZE_STEPS {
        return Err(CompoundError::LimitExceeded {
            resource: "minimize_steps",
            requested: max_steps,
            max: MAX_MINIMIZE_STEPS,
        });
    }
    let mut tried = 1usize;
    if !evaluate_stable("minimize", id, input, fails)? {
        return Err(CompoundError::NotFailing { id: id.to_string() });
    }
    let mut current = input.clone();
    let mut steps = 0usize;
    let mut converged = false;
    'outer: loop {
        let current_before_candidates = canonical_bytes(&current)?;
        let candidates = current.shrink_candidates();
        if canonical_bytes(&current)? != current_before_candidates {
            return Err(CompoundError::CallbackIdentityDrift {
                phase: "shrink_candidates",
                identity: id.to_string(),
            });
        }
        if candidates.len() > MAX_SHRINK_CANDIDATES_PER_STEP {
            return Err(CompoundError::LimitExceeded {
                resource: "shrink_candidates_per_step",
                requested: candidates.len(),
                max: MAX_SHRINK_CANDIDATES_PER_STEP,
            });
        }
        for cand in candidates {
            if tried == MAX_MINIMIZE_EVALUATIONS {
                return Err(CompoundError::LimitExceeded {
                    resource: "minimize_evaluations",
                    requested: tried.saturating_add(1),
                    max: MAX_MINIMIZE_EVALUATIONS,
                });
            }
            tried = tried.checked_add(1).ok_or(CompoundError::LimitExceeded {
                resource: "minimize_evaluations",
                requested: usize::MAX,
                max: MAX_MINIMIZE_EVALUATIONS,
            })?;
            let candidate_fails = evaluate_stable("minimize", id, &cand, fails)?;
            if canonical_bytes(&current)? != current_before_candidates {
                return Err(CompoundError::CallbackIdentityDrift {
                    phase: "minimize",
                    identity: id.to_string(),
                });
            }
            if candidate_fails {
                if steps == max_steps {
                    break 'outer;
                }
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
/// The caller supplies the neighbors; this function validates the hard count
/// cap, aggregate canonical-byte cap, and unique canonical labels before making
/// any predicate call. Its order is the caller's order. A final identity pass
/// catches later predicates that mutate an earlier neighbor through shared
/// interior state.
pub fn probe_neighborhood<I: Canon>(
    neighbors: &[(String, I)],
    fails: &dyn Fn(&I) -> bool,
) -> Result<NeighborhoodReport, CompoundError> {
    if neighbors.len() > MAX_NEIGHBOR_PROBES {
        return Err(CompoundError::LimitExceeded {
            resource: "neighbor_probes",
            requested: neighbors.len(),
            max: MAX_NEIGHBOR_PROBES,
        });
    }
    let mut seen = std::collections::BTreeSet::new();
    for (label, _) in neighbors {
        validate_identifier("neighbor_label", label)?;
        if label == "minimized" {
            return Err(CompoundError::DuplicateIdentity {
                field: "neighbor_label",
                value: label.clone(),
            });
        }
        if !seen.insert(label.as_str()) {
            return Err(CompoundError::DuplicateIdentity {
                field: "neighbor_label",
                value: label.clone(),
            });
        }
    }
    let mut snapshots = Vec::with_capacity(neighbors.len());
    let mut canonical_bytes_total = 0usize;
    for (_, input) in neighbors {
        let snapshot = canonical_bytes(input)?;
        canonical_bytes_total = canonical_bytes_total.checked_add(snapshot.len()).ok_or(
            CompoundError::LimitExceeded {
                resource: "canonical_neighborhood_bytes",
                requested: usize::MAX,
                max: MAX_CANONICAL_FAMILY_BYTES,
            },
        )?;
        if canonical_bytes_total > MAX_CANONICAL_FAMILY_BYTES {
            return Err(CompoundError::LimitExceeded {
                resource: "canonical_neighborhood_bytes",
                requested: canonical_bytes_total,
                max: MAX_CANONICAL_FAMILY_BYTES,
            });
        }
        snapshots.push(snapshot);
    }
    let mut probes = Vec::with_capacity(neighbors.len());
    for ((label, input), snapshot) in neighbors.iter().zip(&snapshots) {
        if canonical_bytes(input)? != *snapshot {
            return Err(CompoundError::CallbackIdentityDrift {
                phase: "neighborhood",
                identity: label.clone(),
            });
        }
        probes.push(NeighborProbe {
            label: label.clone(),
            fails: evaluate_stable("neighborhood", label, input, fails)?,
        });
    }
    for ((label, input), snapshot) in neighbors.iter().zip(&snapshots) {
        if canonical_bytes(input)? != *snapshot {
            return Err(CompoundError::CallbackIdentityDrift {
                phase: "neighborhood",
                identity: label.clone(),
            });
        }
    }
    let failing = probes.iter().filter(|p| p.fails).count();
    Ok(NeighborhoodReport { probes, failing })
}

/// Bounded append-only sink supplied to [`Canon`] implementations.
///
/// The sink refuses an oversized append before allocating it. A codec can
/// still perform caller-owned work internally, but it cannot force this crate
/// to retain an unbounded canonical payload.
#[derive(Debug)]
pub struct CanonWriter {
    bytes: Vec<u8>,
    max: usize,
    resource: &'static str,
}

impl CanonWriter {
    fn new(max: usize, resource: &'static str) -> Self {
        Self {
            bytes: Vec::new(),
            max,
            resource,
        }
    }

    /// Append one byte within the configured bound.
    pub fn push(&mut self, byte: u8) -> Result<(), CompoundError> {
        self.extend_from_slice(&[byte])
    }

    /// Append a byte slice, refusing before allocation when it would exceed
    /// the configured bound.
    pub fn extend_from_slice(&mut self, bytes: &[u8]) -> Result<(), CompoundError> {
        let requested =
            self.bytes
                .len()
                .checked_add(bytes.len())
                .ok_or(CompoundError::LimitExceeded {
                    resource: self.resource,
                    requested: usize::MAX,
                    max: self.max,
                })?;
        if requested > self.max {
            return Err(CompoundError::LimitExceeded {
                resource: self.resource,
                requested,
                max: self.max,
            });
        }
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }

    /// Append `count` copies without allocating a temporary buffer.
    pub fn repeat(&mut self, byte: u8, count: usize) -> Result<(), CompoundError> {
        let requested =
            self.bytes
                .len()
                .checked_add(count)
                .ok_or(CompoundError::LimitExceeded {
                    resource: self.resource,
                    requested: usize::MAX,
                    max: self.max,
                })?;
        if requested > self.max {
            return Err(CompoundError::LimitExceeded {
                resource: self.resource,
                requested,
                max: self.max,
            });
        }
        self.bytes.resize(requested, byte);
        Ok(())
    }

    fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

/// Canonical bytes for content addressing. Every implementation declares a
/// globally unique, stable codec id and schema version in addition to its
/// tagged, length-prefixed payload. Floats use `to_bits`, never formatting.
///
/// # Codec trust boundary
/// Implementations must be pure and must encode every field that can affect
/// the failure predicate. The workflow detects persistent mutation around
/// callbacks, but Rust cannot prove that a caller-written codec is complete or
/// that a callback did not transiently mutate and restore hidden state. Use
/// derived/field-by-field codecs and immutable replay inputs for authoritative
/// families; incomplete or impure codecs are explicitly outside the replay
/// authentication claim.
pub trait Canon {
    /// Stable codec id. Changing its meaning requires a new id or schema
    /// version; Rust type names are intentionally not used because they are
    /// not a durable wire contract.
    const TYPE_ID: &'static str;
    /// Version of this type's canonical payload semantics. Zero is reserved.
    const SCHEMA_VERSION: u32 = 1;

    /// Append this value's complete, side-effect-free canonical bytes to the
    /// bounded sink.
    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError>;

    /// Append stable child codec schemas for generic containers. The outer
    /// type's own id/version is always emitted by this crate and cannot be
    /// skipped by an implementation.
    #[doc(hidden)]
    fn canon_child_schemas(_out: &mut CanonWriter) -> Result<(), CompoundError> {
        Ok(())
    }
}

fn append_canon_schema<T: Canon + ?Sized>(out: &mut CanonWriter) -> Result<(), CompoundError> {
    validate_identifier("canon_type_id", T::TYPE_ID)?;
    if T::SCHEMA_VERSION == 0 {
        return Err(CompoundError::InvalidField {
            field: "canon_schema_version",
            problem: format!("codec {:?} uses reserved schema version 0", T::TYPE_ID),
        });
    }
    out.push(13)?;
    out.extend_from_slice(&(T::TYPE_ID.len() as u64).to_le_bytes())?;
    out.extend_from_slice(T::TYPE_ID.as_bytes())?;
    out.extend_from_slice(&T::SCHEMA_VERSION.to_le_bytes())?;
    T::canon_child_schemas(out)
}

fn canon_schema<T: Canon + ?Sized>() -> Result<Vec<u8>, CompoundError> {
    let mut out = CanonWriter::new(MAX_CANONICAL_SCHEMA_BYTES, "canonical_schema_bytes");
    append_canon_schema::<T>(&mut out)?;
    Ok(out.into_bytes())
}

/// Produce one bounded canonical payload using the same path family
/// construction and replay use.
pub fn canonical_bytes<T: Canon + ?Sized>(value: &T) -> Result<Vec<u8>, CompoundError> {
    let mut out = CanonWriter::new(MAX_CANONICAL_MEMBER_BYTES, "canonical_member_bytes");
    value.canon(&mut out)?;
    let bytes = out.into_bytes();
    if bytes.is_empty() {
        return Err(CompoundError::InvalidField {
            field: "member_canon",
            problem: "Canon implementations must emit a non-empty tagged value".to_string(),
        });
    }
    Ok(bytes)
}

impl Canon for u64 {
    const TYPE_ID: &'static str = "org.frankensim.canon.u64";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(1)?;
        out.extend_from_slice(&self.to_le_bytes())
    }
}
impl Canon for u32 {
    const TYPE_ID: &'static str = "org.frankensim.canon.u32";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(11)?;
        out.extend_from_slice(&self.to_le_bytes())
    }
}
impl Canon for i64 {
    const TYPE_ID: &'static str = "org.frankensim.canon.i64";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(2)?;
        out.extend_from_slice(&self.to_le_bytes())
    }
}
impl Canon for i32 {
    const TYPE_ID: &'static str = "org.frankensim.canon.i32";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(3)?;
        out.extend_from_slice(&self.to_le_bytes())
    }
}
impl Canon for f64 {
    const TYPE_ID: &'static str = "org.frankensim.canon.f64-bits";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(4)?;
        out.extend_from_slice(&self.to_bits().to_le_bytes())
    }
}
impl Canon for bool {
    const TYPE_ID: &'static str = "org.frankensim.canon.bool";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(5)?;
        out.push(u8::from(*self))
    }
}
impl Canon for str {
    const TYPE_ID: &'static str = "org.frankensim.canon.str";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(6)?;
        out.extend_from_slice(&(self.len() as u64).to_le_bytes())?;
        out.extend_from_slice(self.as_bytes())
    }
}
impl Canon for String {
    const TYPE_ID: &'static str = "org.frankensim.canon.string";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        self.as_str().canon(out)
    }
}
impl<T: Canon> Canon for Vec<T> {
    const TYPE_ID: &'static str = "org.frankensim.canon.vec";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(7)?;
        out.extend_from_slice(&(self.len() as u64).to_le_bytes())?;
        for item in self {
            item.canon(out)?;
        }
        Ok(())
    }

    fn canon_child_schemas(out: &mut CanonWriter) -> Result<(), CompoundError> {
        append_canon_schema::<T>(out)
    }
}
impl<A: Canon, B: Canon> Canon for (A, B) {
    const TYPE_ID: &'static str = "org.frankensim.canon.tuple2";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        out.push(8)?;
        self.0.canon(out)?;
        self.1.canon(out)
    }

    fn canon_child_schemas(out: &mut CanonWriter) -> Result<(), CompoundError> {
        append_canon_schema::<A>(out)?;
        append_canon_schema::<B>(out)
    }
}

impl Canon for InvariantClass {
    const TYPE_ID: &'static str = "org.frankensim.fs-bisect.invariant-class";

    fn canon(&self, out: &mut CanonWriter) -> Result<(), CompoundError> {
        match self {
            Self::BuildModeDeterminism => out.push(0),
            Self::CrossIsaDeterminism => out.push(1),
            Self::GoldenDrift => out.push(2),
            Self::EnclosureViolation => out.push(3),
            Self::CertificateForgery => out.push(4),
            Self::ConservationViolation => out.push(5),
            Self::AdjointInconsistency => out.push(6),
            Self::Other(name) => {
                out.push(7)?;
                name.canon(out)
            }
        }
    }
}

/// A permanent regression family: the minimized case plus its failing
/// neighbors, with provenance. This is the artifact a failure leaves
/// behind — MORE than one example, linked to its tracking issue, carrying
/// the admission-rule lesson when one generalizes.
#[derive(Debug, Clone)]
pub struct RegressionFamily<I> {
    /// Family name (stable, kebab-case).
    name: String,
    /// The invariant every member violates.
    invariant: InvariantClass,
    /// Reproduction seed and violated contract context.
    provenance: FamilyProvenance,
    /// Labeled members; `members[0]` is the minimized case by convention.
    members: Vec<(String, I)>,
    /// Construction-time canonical snapshots paired one-for-one with members.
    /// Hashes and manifests never re-run a stateful caller implementation.
    member_canon: Vec<Vec<u8>>,
    /// Stable codec/schema domain for `I`, including generic child schemas.
    member_schema: Vec<u8>,
    /// Tracking references (bead ids / issue ids) — never empty for a
    /// landed family; a failure without a paper trail cannot compound.
    tracking: Vec<String>,
    /// The generalized lesson, when there is one (e.g. "lint variable-
    /// exponent powi out of deterministic paths").
    recommended_admission: Option<String>,
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

fn write_json_string(out: &mut String, value: &str) {
    use std::fmt::Write as _;

    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            control if control <= '\u{1f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            ordinary => out.push(ordinary),
        }
    }
    out.push('"');
}

impl<I: Canon> RegressionFamily<I> {
    /// Build a canonical, bounded regression family.
    ///
    /// # Errors
    /// Refuses empty/duplicate/oversized identifiers, an empty tracking set,
    /// an empty member set, or a malformed custom invariant.
    fn new(
        name: String,
        invariant: InvariantClass,
        members: Vec<(String, I)>,
        tracking: Vec<String>,
        recommended_admission: Option<String>,
        provenance: FamilyProvenance,
    ) -> Result<Self, CompoundError> {
        validate_family_name(&name)?;
        invariant.validate()?;
        validate_tracking_refs(&tracking)?;
        validate_admission_rule(recommended_admission.as_deref())?;
        let member_schema = canon_schema::<I>()?;
        if members.is_empty() {
            return Err(CompoundError::InvalidField {
                field: "members",
                problem: "a regression family must retain at least its minimized case".to_string(),
            });
        }
        let max_members = MAX_NEIGHBOR_PROBES + 1;
        if members.len() > max_members {
            return Err(CompoundError::LimitExceeded {
                resource: "family_members",
                requested: members.len(),
                max: max_members,
            });
        }
        if members[0].0 != "minimized" {
            return Err(CompoundError::InvalidField {
                field: "members",
                problem: "the first member must be labeled \"minimized\"".to_string(),
            });
        }
        let mut member_labels = std::collections::BTreeSet::new();
        let mut member_canon = Vec::with_capacity(members.len());
        let mut canonical_family_bytes = 0usize;
        for (label, input) in &members {
            validate_identifier("member_label", label)?;
            if !member_labels.insert(label.as_str()) {
                return Err(CompoundError::DuplicateIdentity {
                    field: "member_label",
                    value: label.clone(),
                });
            }
            let canonical = canonical_bytes(input)?;
            canonical_family_bytes = canonical_family_bytes.checked_add(canonical.len()).ok_or(
                CompoundError::LimitExceeded {
                    resource: "canonical_family_bytes",
                    requested: usize::MAX,
                    max: MAX_CANONICAL_FAMILY_BYTES,
                },
            )?;
            if canonical_family_bytes > MAX_CANONICAL_FAMILY_BYTES {
                return Err(CompoundError::LimitExceeded {
                    resource: "canonical_family_bytes",
                    requested: canonical_family_bytes,
                    max: MAX_CANONICAL_FAMILY_BYTES,
                });
            }
            member_canon.push(canonical);
        }
        Ok(Self {
            name,
            invariant,
            provenance,
            members,
            member_canon,
            member_schema,
            tracking,
            recommended_admission,
        })
    }

    /// Stable family name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Classified invariant.
    #[must_use]
    pub fn invariant(&self) -> &InvariantClass {
        &self.invariant
    }

    /// Captured seed and violated contract context.
    #[must_use]
    pub fn provenance(&self) -> &FamilyProvenance {
        &self.provenance
    }

    /// Number of sealed members (the minimized case followed by failing
    /// neighbors).
    #[must_use]
    pub fn member_count(&self) -> usize {
        self.members.len()
    }

    /// Stable label for one sealed member. Typed member values are not
    /// exposed because mutating interior state must not bypass replay's
    /// identity check.
    #[must_use]
    pub fn member_label(&self, index: usize) -> Option<&str> {
        self.members.get(index).map(|(label, _)| label.as_str())
    }

    /// Beads or issue references that own the family.
    #[must_use]
    pub fn tracking(&self) -> &[String] {
        &self.tracking
    }

    /// Generalized admission recommendation, if one was justified.
    #[must_use]
    pub fn recommended_admission(&self) -> Option<&str> {
        self.recommended_admission.as_deref()
    }

    /// Content hash over the full canonical encoding: name, invariant,
    /// member codec/schema, members (labels + inputs), tracking, and admission.
    /// Deterministic across runs and build modes; cross-ISA equality additionally
    /// requires the member codec and predicate domain to make that claim.
    #[must_use]
    pub fn content_hash(&self) -> fs_blake3::ContentHash {
        let mut bytes = CanonWriter::new(usize::MAX, "sealed_family_hash_bytes");
        COMPOUND_CANON_VERSION
            .canon(&mut bytes)
            .expect("sealed canon version is bounded");
        self.name
            .canon(&mut bytes)
            .expect("validated family name is bounded");
        self.invariant
            .canon(&mut bytes)
            .expect("validated invariant is bounded");
        self.provenance
            .seed
            .canon(&mut bytes)
            .expect("seed is fixed width");
        self.provenance
            .contract
            .canon(&mut bytes)
            .expect("validated contract is bounded");
        self.provenance
            .detail
            .canon(&mut bytes)
            .expect("validated detail is bounded");
        bytes.push(14).expect("unbounded sealed writer");
        bytes
            .extend_from_slice(&(self.member_schema.len() as u64).to_le_bytes())
            .expect("unbounded sealed writer");
        bytes
            .extend_from_slice(&self.member_schema)
            .expect("sealed schema is bounded");
        bytes.push(7).expect("unbounded sealed writer");
        bytes
            .extend_from_slice(&(self.members.len() as u64).to_le_bytes())
            .expect("unbounded sealed writer");
        for ((label, _), canonical) in self.members.iter().zip(&self.member_canon) {
            label
                .canon(&mut bytes)
                .expect("validated member label is bounded");
            bytes.push(12).expect("unbounded sealed writer");
            bytes
                .extend_from_slice(
                    &u64::try_from(canonical.len())
                        .expect("bounded canonical member length fits u64")
                        .to_le_bytes(),
                )
                .expect("unbounded sealed writer");
            bytes
                .extend_from_slice(canonical)
                .expect("sealed canonical member is bounded");
        }
        self.tracking
            .canon(&mut bytes)
            .expect("validated tracking references are bounded");
        match &self.recommended_admission {
            Some(a) => {
                bytes.push(9).expect("unbounded sealed writer");
                a.canon(&mut bytes)
                    .expect("validated admission rule is bounded");
            }
            None => bytes.push(10).expect("unbounded sealed writer"),
        }
        fs_blake3::hash_domain(COMPOUND_FAMILY_HASH_DOMAIN, &bytes.into_bytes())
    }

    /// The canonical capture manifest: JSON-lines, one header, one line per
    /// member (canonical bytes hex-encoded), one trailer with the content hash.
    /// Decoding arbitrary caller-defined member types remains the family
    /// owner's responsibility.
    #[must_use]
    pub fn manifest(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"canon_version\":{COMPOUND_CANON_VERSION},\"family\":"
        );
        write_json_string(&mut out, &self.name);
        let _ = write!(out, ",\"invariant\":");
        write_json_string(&mut out, self.invariant.name());
        let _ = write!(out, ",\"seed\":{},\"contract\":", self.provenance.seed);
        write_json_string(&mut out, &self.provenance.contract);
        let _ = write!(out, ",\"detail\":");
        write_json_string(&mut out, &self.provenance.detail);
        let _ = write!(out, ",\"member_type\":");
        write_json_string(&mut out, I::TYPE_ID);
        let schema_hex: String = self
            .member_schema
            .iter()
            .fold(String::new(), |mut s, byte| {
                let _ = write!(s, "{byte:02x}");
                s
            });
        let _ = write!(
            out,
            ",\"member_schema_version\":{},\"member_schema\":\"{schema_hex}\"",
            I::SCHEMA_VERSION
        );
        let _ = write!(out, ",\"members\":{},\"tracking\":[", self.members.len());
        for (index, reference) in self.tracking.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            write_json_string(&mut out, reference);
        }
        let _ = write!(out, "],\"recommended_admission\":");
        if let Some(rule) = &self.recommended_admission {
            write_json_string(&mut out, rule);
        } else {
            out.push_str("null");
        }
        out.push_str("}\n");
        for ((label, _), canonical) in self.members.iter().zip(&self.member_canon) {
            let hex: String = canonical.iter().fold(String::new(), |mut s, b| {
                let _ = write!(s, "{b:02x}");
                s
            });
            out.push_str("{\"member\":");
            write_json_string(&mut out, label);
            let _ = writeln!(out, ",\"canon\":\"{hex}\"}}");
        }
        let _ = writeln!(out, "{{\"content_hash\":\"{}\"}}", self.content_hash());
        out
    }

    /// Re-execute every member against the predicate. All live typed values
    /// are re-canonicalized and matched to the sealed snapshots before any
    /// predicate work. A live family has `now_passing` empty; anything else is
    /// stale evidence to act on.
    ///
    /// # Errors
    /// Refuses replay if the codec schema or any live member no longer matches
    /// the content-addressed identity.
    pub fn replay(&self, fails: &dyn Fn(&I) -> bool) -> Result<ReplayReport, CompoundError> {
        if canon_schema::<I>()? != self.member_schema {
            return Err(CompoundError::ReplayIdentityDrift { member: None });
        }
        self.verify_live_members()?;
        let mut still_failing = Vec::new();
        let mut now_passing = Vec::new();
        for ((label, input), canonical) in self.members.iter().zip(&self.member_canon) {
            if canonical_bytes(input)? != *canonical {
                return Err(CompoundError::ReplayIdentityDrift {
                    member: Some(label.clone()),
                });
            }
            let failed = fails(input);
            if canonical_bytes(input)? != *canonical {
                return Err(CompoundError::ReplayIdentityDrift {
                    member: Some(label.clone()),
                });
            }
            if failed {
                still_failing.push(label.clone());
            } else {
                now_passing.push(label.clone());
            }
        }
        // A later callback can share interior state with and mutate an earlier
        // member; the final pass closes that persistent cross-member TOCTOU.
        self.verify_live_members()?;
        Ok(ReplayReport {
            still_failing,
            now_passing,
        })
    }

    fn verify_live_members(&self) -> Result<(), CompoundError> {
        for ((label, input), canonical) in self.members.iter().zip(&self.member_canon) {
            if canonical_bytes(input)? != *canonical {
                return Err(CompoundError::ReplayIdentityDrift {
                    member: Some(label.clone()),
                });
            }
        }
        Ok(())
    }
}

/// The full workflow output: minimized case, neighborhood, family, hash.
#[derive(Debug, Clone)]
pub struct CompoundReport<I> {
    /// The captured case with its input replaced by the minimum.
    pub case: FailureCase<I>,
    /// Accepted minimization steps.
    pub steps: usize,
    /// Predicate evaluations during minimization, including the seed input.
    pub tried: usize,
    /// Whether minimization reached a fixpoint.
    pub converged: bool,
    /// The bounded neighborhood around the minimum.
    pub neighborhood: NeighborhoodReport,
    /// The landed family (minimum first, then failing neighbors).
    pub family: RegressionFamily<I>,
    /// The family's content hash (also in the manifest trailer).
    pub content_hash: fs_blake3::ContentHash,
}

/// The v3 workflow driver: validate → minimize → probe → seal the family.
///
/// `neighbors_of` receives the MINIMIZED input and returns a deterministically
/// ordered, labeled neighbor set. Its callback work is caller-owned; the
/// returned set is count/identity-validated before any neighbor predicate is
/// evaluated. Failing neighbors join the family behind the minimum.
///
/// # Errors
/// [`CompoundError::NotFailing`] when the captured input does not fail,
/// [`CompoundError::MinimizationIncomplete`] when the caller budget expires
/// before a fixpoint, plus structured field, identity, and deterministic
/// work-limit refusals.
pub fn compound<I: Shrink + Canon>(
    case: FailureCase<I>,
    fails: &dyn Fn(&I) -> bool,
    neighbors_of: &dyn Fn(&I) -> Vec<(String, I)>,
    tracking: Vec<String>,
    recommended_admission: Option<String>,
    max_steps: usize,
) -> Result<CompoundReport<I>, CompoundError> {
    validate_family_name(&case.id)?;
    case.invariant.validate()?;
    let provenance = FamilyProvenance::new(case.seed, case.contract.clone(), case.detail.clone())?;
    validate_tracking_refs(&tracking)?;
    validate_admission_rule(recommended_admission.as_deref())?;
    let report = minimize(&case.id, &case.input, fails, max_steps)?;
    if !report.converged {
        return Err(CompoundError::MinimizationIncomplete {
            id: case.id,
            steps: report.steps,
            evaluations: report.tried,
        });
    }
    let minimized_before_neighbors = canonical_bytes(&report.minimized)?;
    let neighbors = neighbors_of(&report.minimized);
    if canonical_bytes(&report.minimized)? != minimized_before_neighbors {
        return Err(CompoundError::CallbackIdentityDrift {
            phase: "neighbors_of",
            identity: case.id,
        });
    }
    let neighborhood = probe_neighborhood(&neighbors, fails)?;
    let mut members: Vec<(String, I)> = vec![("minimized".to_string(), report.minimized.clone())];
    for ((label, input), probe) in neighbors.into_iter().zip(&neighborhood.probes) {
        if probe.fails {
            members.push((label, input));
        }
    }
    let family = RegressionFamily::new(
        case.id.clone(),
        case.invariant.clone(),
        members,
        tracking,
        recommended_admission,
        provenance,
    )?;
    let content_hash = family.content_hash();
    Ok(CompoundReport {
        case: FailureCase {
            input: report.minimized,
            ..case
        },
        steps: report.steps,
        tried: report.tried,
        converged: true,
        neighborhood,
        family,
        content_hash,
    })
}
