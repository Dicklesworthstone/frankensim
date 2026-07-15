//! FALSIFIER PAIRING (addendum Proposal 6): a catalog of intended independent
//! checks plus diagnostic consequence × doubt budget telemetry. This module is
//! deliberately NOT release authority: names and method prose cannot prove an
//! executable binding, independence, or that the exact claim instance survived
//! a retained run. Authenticated instance admission lives in the package/checker
//! release boundary.
//! Certificates prove what the model claims; falsifiers probe whether
//! the claims connect to reality; the gap between them is where
//! simulation systems silently rot. Popper as infrastructure.
//!
//! The catalog refuses empty or malformed declarations, the allocator retains
//! honest cold-start boundaries, discrepancy reports produce pending
//! tombstone/bug candidates,
//! and preliminary class-level yield tracking supports review. Typed
//! per-falsifier/window receipts and policy-checked independence remain required
//! before this telemetry can govern release.

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Doubt never reaches zero: even a perfect record keeps a floor of
/// falsification pressure (the record could be luck or blind spots).
pub const DOUBT_FLOOR: f64 = 0.05;

/// Cold-start doubt: a class with NO history in a regime is maximally
/// doubted, never trusted by default.
pub const DOUBT_COLD_START: f64 = 1.0;

/// A claim with no downstream dependents still gets a minimal-but-
/// nonzero consequence weight (someone may read it directly).
pub const CONSEQUENCE_FLOOR: f64 = 0.01;

/// Yield threshold: at or above this many newly observed runs in a review
/// window, a class with zero reported discrepancies has meaningful volume.
pub const RENT_VOLUME: u64 = 100;

/// Budget-share multiplier applied per class-level rent review (never to zero).
pub const RENT_DECAY: f64 = 0.5;

/// Floor on the decayed share multiplier.
pub const RENT_SHARE_FLOOR: f64 = 0.1;

/// Per-row error level used by the time-uniform union-bound heuristic that
/// keeps finite pass histories from collapsing doubt prematurely. It does not
/// control multiplicity across classes/regimes.
pub const DOUBT_ALPHA: f64 = 0.05;

/// Defensive bound on one diagnostic allocation request.
pub const MAX_CLAIMS_PER_ALLOCATION: usize = 65_536;

/// Defensive bound on one declaration-catalog lint request.
pub const MAX_CATALOG_CLASSES_PER_LINT: usize = 262_144;

/// Maximum number of certificate classes in the in-memory declaration catalog.
pub const MAX_CATALOG_CLASSES: usize = 16_384;

/// Maximum declarations admitted for one certificate class.
pub const MAX_FALSIFIERS_PER_CLASS: usize = 256;

/// Maximum aggregate declaration name/method bytes for one class.
pub const MAX_CATALOG_TEXT_BYTES_PER_CLASS: usize = 1_048_576;

/// Maximum aggregate retained identifier/declaration bytes in one catalog.
pub const MAX_CATALOG_TEXT_BYTES_TOTAL: usize = 16_777_216;

/// Maximum aggregate bytes returned by one missing-class lint.
pub const MAX_CATALOG_LINT_OUTPUT_BYTES: usize = 16_777_216;

/// Maximum missing identifiers materialized by one catalog lint.
pub const MAX_CATALOG_LINT_OUTPUT_CLASSES: usize = 65_536;

/// Defensive bound on distinct `(class, regime)` telemetry rows. This legacy
/// in-memory accumulator has no cancellation context, so admission is bounded
/// before a new caller-controlled key can grow it without limit.
pub const MAX_HISTORY_ROWS: usize = 65_536;

/// Maximum distinct idempotent diagnostic attempts retained in memory.
pub const MAX_ATTEMPTS: usize = 65_536;

/// Maximum aggregate retained class/regime key bytes in diagnostic history.
pub const MAX_HISTORY_KEY_BYTES_TOTAL: usize = 16_777_216;

/// Maximum aggregate retained attempt-identity and detail bytes.
pub const MAX_ATTEMPT_BYTES_TOTAL: usize = 16_777_216;

const MAX_IDENTIFIER_BYTES: usize = 1_024;
const MAX_DETAIL_BYTES: usize = 65_536;

/// One catalog declaration for an intended independent check. Public strings
/// carry no executable or release authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FalsifierSpec {
    /// Stable nonempty catalog name.
    pub name: String,
    /// The independent method, stated (audit text).
    pub method: String,
}

/// Validation failure for falsifier-catalog or diagnostic-history data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FalsifyError {
    /// A certificate class tried to register without any falsifier.
    NoFalsifier {
        /// The offending class.
        class: String,
    },
    /// Duplicate class registration.
    Duplicate {
        /// The class.
        class: String,
    },
    /// A query named an unregistered class.
    Unknown {
        /// The class.
        class: String,
    },
    /// An identifier is empty, oversized, has a non-alphanumeric first byte,
    /// or contains a byte outside the stable ASCII-slug alphabet.
    InvalidIdentifier {
        /// Field name.
        field: &'static str,
    },
    /// Caller-provided audit text is empty or exceeds its field-specific cap.
    InvalidText {
        /// Field name.
        field: &'static str,
        /// Public byte cap for this field.
        max_bytes: usize,
    },
    /// A catalog entry has an empty/oversized name or method.
    InvalidSpec {
        /// Certificate class being registered.
        class: String,
        /// Entry index.
        index: usize,
    },
    /// A catalog repeats the same falsifier name within one class.
    DuplicateSpec {
        /// Certificate class.
        class: String,
        /// Duplicate falsifier name.
        name: String,
    },
    /// A compute charge, total budget, consequence, or accumulated total is
    /// negative or non-finite.
    InvalidNumber {
        /// Numeric field.
        field: &'static str,
        /// Exact IEEE-754 bits.
        bits: u64,
    },
    /// A strictly positive finite compute charge is zero, negative, or
    /// non-finite.
    InvalidPositiveNumber {
        /// Numeric field.
        field: &'static str,
        /// Exact IEEE-754 bits.
        bits: u64,
    },
    /// Checked history arithmetic could not represent the next state.
    HistoryOverflow {
        /// Certificate class.
        class: String,
        /// Regime key.
        regime: String,
    },
    /// An allocation vector exceeds its defensive legacy-path ceiling.
    TooManyClaims {
        /// Requested claim count.
        requested: usize,
        /// Public cap.
        cap: usize,
    },
    /// A catalog-lint request exceeds its defensive ceiling.
    TooManyCatalogClasses {
        /// Requested class count.
        requested: usize,
        /// Public cap.
        cap: usize,
    },
    /// A new distinct telemetry row would exceed the in-memory history cap.
    TooManyHistoryRows {
        /// Current row count (also the rejected postcondition boundary).
        current: usize,
        /// Public cap.
        cap: usize,
    },
    /// A caller-controlled collection or byte budget exceeds a named ceiling.
    ResourceLimit {
        /// Limited field or collection.
        field: &'static str,
        /// Requested count or bytes.
        requested: usize,
        /// Public cap.
        cap: usize,
    },
    /// A bounded output allocation could not be reserved.
    ResourceExhausted {
        /// Operation whose output allocation failed.
        operation: &'static str,
    },
    /// An attempt names a falsifier that is not declared for its class.
    UnregisteredFalsifier {
        /// Certificate class.
        class: String,
        /// Falsifier name.
        falsifier: String,
    },
    /// One stable attempt ID was reused for different content.
    AttemptCollision {
        /// Conflicting stable attempt ID.
        attempt_id: String,
    },
}

impl core::fmt::Display for FalsifyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FalsifyError::NoFalsifier { class } => write!(
                f,
                "certificate class {class:?} cannot enter the falsifier catalog \
                 without at least one declaration"
            ),
            FalsifyError::Duplicate { class } => write!(f, "class {class:?} already registered"),
            FalsifyError::Unknown { class } => write!(f, "unknown certificate class {class:?}"),
            FalsifyError::InvalidIdentifier { field } => {
                write!(
                    f,
                    "{field} must be an ASCII slug of at most {MAX_IDENTIFIER_BYTES} bytes, beginning with an alphanumeric byte"
                )
            }
            FalsifyError::InvalidText { field, max_bytes } => {
                write!(f, "{field} must be nonempty and at most {max_bytes} bytes")
            }
            FalsifyError::InvalidSpec { class, index } => write!(
                f,
                "certificate class {class:?} has a malformed falsifier declaration at index {index}"
            ),
            FalsifyError::DuplicateSpec { class, name } => write!(
                f,
                "certificate class {class:?} repeats falsifier declaration {name:?}"
            ),
            FalsifyError::InvalidNumber { field, bits } => write!(
                f,
                "{field} must be finite and non-negative (bits 0x{bits:016x})"
            ),
            FalsifyError::InvalidPositiveNumber { field, bits } => write!(
                f,
                "{field} must be finite and strictly positive (bits 0x{bits:016x})"
            ),
            FalsifyError::HistoryOverflow { class, regime } => write!(
                f,
                "falsifier history overflow for class {class:?} in regime {regime:?}"
            ),
            FalsifyError::TooManyClaims { requested, cap } => write!(
                f,
                "falsifier allocation has {requested} claims, above defensive cap {cap}"
            ),
            FalsifyError::TooManyCatalogClasses { requested, cap } => write!(
                f,
                "falsifier catalog lint has {requested} classes, above defensive cap {cap}"
            ),
            FalsifyError::TooManyHistoryRows { current, cap } => write!(
                f,
                "falsifier telemetry already has {current} distinct rows; defensive cap is {cap}"
            ),
            FalsifyError::ResourceLimit {
                field,
                requested,
                cap,
            } => write!(f, "{field} has size {requested}, above defensive cap {cap}"),
            FalsifyError::ResourceExhausted { operation } => {
                write!(f, "could not reserve bounded output for {operation}")
            }
            FalsifyError::UnregisteredFalsifier { class, falsifier } => write!(
                f,
                "falsifier {falsifier:?} is not declared for certificate class {class:?}"
            ),
            FalsifyError::AttemptCollision { attempt_id } => write!(
                f,
                "falsifier attempt ID {attempt_id:?} was reused for different content"
            ),
        }
    }
}

impl std::error::Error for FalsifyError {}

/// The falsifier registry: certificate class → its independent checks.
#[derive(Debug, Default)]
pub struct FalsifierRegistry {
    classes: BTreeMap<String, Vec<FalsifierSpec>>,
    text_bytes: usize,
}

impl FalsifierRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> Self {
        FalsifierRegistry::default()
    }

    /// The starting declaration catalog from the proposal. Entries state
    /// intended checker families; they do not assert implementations exist or
    /// confer release authority.
    #[must_use]
    pub fn standard() -> Self {
        let mut r = FalsifierRegistry::new();
        let pairs: [(&str, &str, &str); 7] = [
            (
                "sampled-interface-agreement",
                "independent-retained-sample-replay",
                "separately structured checker re-evaluates retained chart samples and interval predicates against source-bound replay artifacts",
            ),
            (
                "continuum-watertightness",
                "certified-oriented-intersection-winding",
                "certified oriented intersections plus winding/degree and coverage-complete interval subdivision on an independently structured checker path",
            ),
            (
                "conservation",
                "global-flux-audit",
                "independent global flux balance on a DIFFERENT quadrature rule",
            ),
            (
                "adjoint-gradient",
                "finite-difference-spot-check",
                "central differences along random directions, independent of the tape",
            ),
            (
                "surrogate-accept",
                "held-out-point-evaluation",
                "full-fidelity evaluation at points the surrogate never saw",
            ),
            (
                "symmetry-block-solve",
                "occasional-full-solve",
                "solve the UNREDUCED system on random instances and compare",
            ),
            (
                "validated-color",
                "held-out-experimental-anchor",
                "compare against experimental anchors withheld from calibration",
            ),
        ];
        for (class, name, method) in pairs {
            r.register(
                class,
                vec![FalsifierSpec {
                    name: name.to_string(),
                    method: method.to_string(),
                }],
            )
            .expect("standard registry is well-formed");
        }
        r
    }

    /// Register a certificate class with intended falsifier declarations.
    /// Refuses empty/malformed identifiers, empty/oversized method text, and
    /// duplicate names. Registration remains catalog metadata, not execution.
    ///
    /// # Errors
    /// Returns a structured error for malformed/duplicate declarations or any
    /// class, per-class declaration, or aggregate text resource ceiling.
    pub fn register(
        &mut self,
        class: &str,
        falsifiers: Vec<FalsifierSpec>,
    ) -> Result<(), FalsifyError> {
        validate_identifier("certificate class", class)?;
        if falsifiers.is_empty() {
            return Err(FalsifyError::NoFalsifier {
                class: class.to_string(),
            });
        }
        if self.classes.contains_key(class) {
            return Err(FalsifyError::Duplicate {
                class: class.to_string(),
            });
        }
        if self.classes.len() >= MAX_CATALOG_CLASSES {
            return Err(FalsifyError::ResourceLimit {
                field: "falsifier catalog classes",
                requested: self.classes.len().saturating_add(1),
                cap: MAX_CATALOG_CLASSES,
            });
        }
        if falsifiers.len() > MAX_FALSIFIERS_PER_CLASS {
            return Err(FalsifyError::ResourceLimit {
                field: "falsifiers per certificate class",
                requested: falsifiers.len(),
                cap: MAX_FALSIFIERS_PER_CLASS,
            });
        }
        let mut names = std::collections::BTreeSet::new();
        let mut text_bytes = class.len();
        for (index, spec) in falsifiers.iter().enumerate() {
            if validate_identifier("falsifier name", &spec.name).is_err()
                || spec.method.trim().is_empty()
                || spec.method.len() > MAX_DETAIL_BYTES
            {
                return Err(FalsifyError::InvalidSpec {
                    class: class.to_string(),
                    index,
                });
            }
            text_bytes = text_bytes
                .checked_add(spec.name.len())
                .and_then(|bytes| bytes.checked_add(spec.method.len()))
                .ok_or(FalsifyError::ResourceLimit {
                    field: "falsifier declaration bytes per class",
                    requested: usize::MAX,
                    cap: MAX_CATALOG_TEXT_BYTES_PER_CLASS,
                })?;
            if text_bytes > MAX_CATALOG_TEXT_BYTES_PER_CLASS {
                return Err(FalsifyError::ResourceLimit {
                    field: "falsifier declaration bytes per class",
                    requested: text_bytes,
                    cap: MAX_CATALOG_TEXT_BYTES_PER_CLASS,
                });
            }
            if !names.insert(spec.name.clone()) {
                return Err(FalsifyError::DuplicateSpec {
                    class: class.to_string(),
                    name: spec.name.clone(),
                });
            }
        }
        let next_total =
            self.text_bytes
                .checked_add(text_bytes)
                .ok_or(FalsifyError::ResourceLimit {
                    field: "aggregate falsifier catalog text bytes",
                    requested: usize::MAX,
                    cap: MAX_CATALOG_TEXT_BYTES_TOTAL,
                })?;
        if next_total > MAX_CATALOG_TEXT_BYTES_TOTAL {
            return Err(FalsifyError::ResourceLimit {
                field: "aggregate falsifier catalog text bytes",
                requested: next_total,
                cap: MAX_CATALOG_TEXT_BYTES_TOTAL,
            });
        }
        self.classes.insert(class.to_string(), falsifiers);
        self.text_bytes = next_total;
        Ok(())
    }

    /// The falsifiers for a class.
    ///
    /// # Errors
    /// [`FalsifyError::Unknown`].
    pub fn falsifiers(&self, class: &str) -> Result<&[FalsifierSpec], FalsifyError> {
        validate_identifier("certificate class", class)?;
        self.classes
            .get(class)
            .map(Vec::as_slice)
            .ok_or_else(|| FalsifyError::Unknown {
                class: class.to_string(),
            })
    }

    fn declares(&self, class: &str, falsifier: &str) -> bool {
        self.classes
            .get(class)
            .is_some_and(|specs| specs.iter().any(|spec| spec.name.as_str() == falsifier))
    }

    /// Catalog-completeness lint: every distinct named class must have a
    /// declaration. Output is unique and sorted by the canonical class slug,
    /// so duplicate caller entries cannot amplify diagnostics.
    /// Empty output means only that catalog metadata is present. It does not
    /// prove executable binding, independence, an exact-instance run, retained
    /// artifacts, or release admissibility.
    ///
    /// # Errors
    /// Returns [`FalsifyError::InvalidIdentifier`] for a malformed class and
    /// [`FalsifyError::TooManyCatalogClasses`] when the bounded request cannot
    /// be admitted.
    pub fn catalog_gate(&self, certificate_classes: &[&str]) -> Result<Vec<String>, FalsifyError> {
        if certificate_classes.len() > MAX_CATALOG_CLASSES_PER_LINT {
            return Err(FalsifyError::TooManyCatalogClasses {
                requested: certificate_classes.len(),
                cap: MAX_CATALOG_CLASSES_PER_LINT,
            });
        }
        let mut unique_classes = Vec::new();
        unique_classes
            .try_reserve_exact(certificate_classes.len())
            .map_err(|_| FalsifyError::ResourceExhausted {
                operation: "falsifier catalog lint input",
            })?;
        for class in certificate_classes {
            validate_identifier("certificate class", class)?;
            unique_classes.push(*class);
        }
        unique_classes.sort_unstable();
        unique_classes.dedup();

        let mut missing = Vec::new();
        missing
            .try_reserve_exact(unique_classes.len().min(MAX_CATALOG_LINT_OUTPUT_CLASSES))
            .map_err(|_| FalsifyError::ResourceExhausted {
                operation: "falsifier catalog lint",
            })?;
        let mut output_bytes = 0usize;
        for class in unique_classes {
            if !self.classes.contains_key(class) {
                if missing.len() >= MAX_CATALOG_LINT_OUTPUT_CLASSES {
                    return Err(FalsifyError::ResourceLimit {
                        field: "falsifier catalog lint missing classes",
                        requested: missing.len().saturating_add(1),
                        cap: MAX_CATALOG_LINT_OUTPUT_CLASSES,
                    });
                }
                output_bytes =
                    output_bytes
                        .checked_add(class.len())
                        .ok_or(FalsifyError::ResourceLimit {
                            field: "falsifier catalog lint output bytes",
                            requested: usize::MAX,
                            cap: MAX_CATALOG_LINT_OUTPUT_BYTES,
                        })?;
                if output_bytes > MAX_CATALOG_LINT_OUTPUT_BYTES {
                    return Err(FalsifyError::ResourceLimit {
                        field: "falsifier catalog lint output bytes",
                        requested: output_bytes,
                        cap: MAX_CATALOG_LINT_OUTPUT_BYTES,
                    });
                }
                missing.push(class.to_string());
            }
        }
        Ok(missing)
    }
}

fn validate_identifier(field: &'static str, value: &str) -> Result<(), FalsifyError> {
    let allowed = |byte: u8| {
        byte.is_ascii_alphanumeric()
            || matches!(byte, b'-' | b'.' | b'_' | b':' | b'/' | b'@' | b'+' | b'~')
    };
    if value.is_empty()
        || value.len() > MAX_IDENTIFIER_BYTES
        || !value.as_bytes()[0].is_ascii_alphanumeric()
        || value.bytes().any(|byte| !allowed(byte))
    {
        return Err(FalsifyError::InvalidIdentifier { field });
    }
    Ok(())
}

fn validate_nonnegative_finite(field: &'static str, value: f64) -> Result<(), FalsifyError> {
    if !value.is_finite() || value < 0.0 {
        return Err(FalsifyError::InvalidNumber {
            field,
            bits: value.to_bits(),
        });
    }
    Ok(())
}

fn validate_positive_finite(field: &'static str, value: f64) -> Result<(), FalsifyError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(FalsifyError::InvalidPositiveNumber {
            field,
            bits: value.to_bits(),
        });
    }
    Ok(())
}

/// Per-(class, regime) falsification history: the DOUBT source.
#[derive(Debug, Default)]
pub struct FalsifierHistory {
    /// (class, regime-key) → (passes, hits, compute spent).
    rows: BTreeMap<(String, String), (u64, u64, f64)>,
    /// Stable attempt ID → immutable content fingerprint. Replays of identical
    /// content are idempotent; conflicting reuse is rejected.
    attempts: BTreeMap<String, AttemptFingerprint>,
    /// Aggregate retained class/regime key bytes.
    history_key_bytes: usize,
    /// Aggregate retained attempt fingerprint bytes.
    attempt_bytes: usize,
    /// Rent-decay multipliers per class (1.0 until decayed).
    share: BTreeMap<String, f64>,
    /// Invocation-schedule-independent fixed-volume rent windows.
    rent_windows: BTreeMap<String, RentWindowState>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AttemptFingerprint {
    class: String,
    regime: String,
    falsifier: String,
    claim_revision: String,
    artifact_id: String,
    seed: u64,
    compute_bits: u64,
    outcome: FalsifierOutcome,
}

#[derive(Debug, Clone, Default)]
struct RentWindowState {
    current_runs: u64,
    current_hits: u64,
    closed_window_hits: Vec<u64>,
}

/// Typed outcome of one exact diagnostic falsifier attempt. A discrepancy is a
/// candidate observation, not an adjudicated proof that an estimator is wrong.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FalsifierOutcome {
    /// The declared checker found no discrepancy in this attempt.
    NoDiscrepancy,
    /// The declared checker reported audit text requiring adjudication.
    Discrepancy {
        /// What disagreed.
        detail: String,
    },
}

/// Immutable identity and accounting envelope for one diagnostic attempt.
/// All identifiers are stable ASCII slugs; `attempt_id` is the idempotency key.
#[derive(Debug, Clone, PartialEq)]
pub struct FalsifierAttempt {
    /// Stable unique attempt identity.
    pub attempt_id: String,
    /// Certificate class.
    pub class: String,
    /// Regime key the claim lived in.
    pub regime: String,
    /// Which falsifier caught it.
    pub falsifier: String,
    /// Caller-asserted claim or claim-set revision reference exercised.
    pub claim_revision: String,
    /// Caller-asserted retained input/output evidence artifact reference.
    pub artifact_id: String,
    /// Explicit deterministic seed.
    pub seed: u64,
    /// Positive finite compute charge.
    pub compute_s: f64,
    /// Observed diagnostic outcome.
    pub outcome: FalsifierOutcome,
}

/// Deterministic JSON for a *candidate* tombstone. A policy authority must bind
/// it to an exact claim and independently adjudicate it before activation.
#[derive(Debug, Clone, PartialEq)]
pub struct TombstoneCandidate {
    json: String,
}

impl TombstoneCandidate {
    /// Schema-versioned, fixed-order JSON payload.
    #[must_use]
    pub fn json(&self) -> &str {
        &self.json
    }
}

/// Deterministic JSON for a candidate estimator bug report.
#[derive(Debug, Clone, PartialEq)]
pub struct EstimatorBugCandidate {
    json: String,
}

impl EstimatorBugCandidate {
    /// Schema-versioned, fixed-order JSON payload.
    #[must_use]
    pub fn json(&self) -> &str {
        &self.json
    }
}

/// Result of ingesting one attempt. An identical replay returns the same
/// deterministic projections with `duplicate=true` and does not mutate state.
#[derive(Debug, Clone, PartialEq)]
pub struct AttemptRecord {
    /// Whether this exact attempt was already present.
    pub duplicate: bool,
    /// Present only for discrepancy outcomes.
    pub tombstone: Option<TombstoneCandidate>,
    /// Present only for discrepancy outcomes.
    pub estimator_bug: Option<EstimatorBugCandidate>,
}

fn attempt_retained_bytes(attempt: &FalsifierAttempt) -> Result<usize, FalsifyError> {
    let detail_bytes = match &attempt.outcome {
        FalsifierOutcome::NoDiscrepancy => 0,
        FalsifierOutcome::Discrepancy { detail } => detail.len(),
    };
    [
        attempt.attempt_id.len(),
        attempt.class.len(),
        attempt.regime.len(),
        attempt.falsifier.len(),
        attempt.claim_revision.len(),
        attempt.artifact_id.len(),
        detail_bytes,
    ]
    .into_iter()
    .try_fold(0usize, |acc, bytes| {
        acc.checked_add(bytes).ok_or(FalsifyError::ResourceLimit {
            field: "aggregate falsifier attempt bytes",
            requested: usize::MAX,
            cap: MAX_ATTEMPT_BYTES_TOTAL,
        })
    })
}

fn candidate_record(attempt: &FalsifierAttempt, duplicate: bool) -> AttemptRecord {
    let FalsifierOutcome::Discrepancy { detail } = &attempt.outcome else {
        return AttemptRecord {
            duplicate,
            tombstone: None,
            estimator_bug: None,
        };
    };
    let attempt_id = crate::json_string(&attempt.attempt_id);
    let class = crate::json_string(&attempt.class);
    let regime = crate::json_string(&attempt.regime);
    let falsifier = crate::json_string(&attempt.falsifier);
    let claim_revision = crate::json_string(&attempt.claim_revision);
    let artifact_id = crate::json_string(&attempt.artifact_id);
    let detail = crate::json_string(detail);
    let seed_bits = crate::json_string(&format!("{:016x}", attempt.seed));
    let compute_bits = crate::json_string(&format!("{:016x}", attempt.compute_s.to_bits()));
    let mut tombstone = String::from(
        "{\"schema_id\":\"fs-evidence/falsifier-candidate\",\"schema_version\":1,\"kind\":\"tombstone-candidate\",\"source\":\"falsifier-attempt-candidate\",\"adjudication\":\"pending\"",
    );
    let _ = write!(
        tombstone,
        ",\"attempt_id\":{attempt_id},\"class\":{class},\"regime\":{regime},\"falsifier\":{falsifier},\"claim_revision\":{claim_revision},\"artifact_id\":{artifact_id},\"seed_bits\":{seed_bits},\"compute_s_bits\":{compute_bits},\"detail\":{detail}}}"
    );
    let mut estimator_bug = String::from(
        "{\"schema_id\":\"fs-evidence/falsifier-candidate\",\"schema_version\":1,\"kind\":\"estimator-bug-candidate\",\"adjudication\":\"pending\"",
    );
    let _ = write!(
        estimator_bug,
        ",\"attempt_id\":{attempt_id},\"class\":{class},\"regime\":{regime},\"caught_by\":{falsifier},\"claim_revision\":{claim_revision},\"artifact_id\":{artifact_id},\"seed_bits\":{seed_bits},\"compute_s_bits\":{compute_bits},\"evidence\":{detail}}}"
    );
    AttemptRecord {
        duplicate,
        tombstone: Some(TombstoneCandidate { json: tombstone }),
        estimator_bug: Some(EstimatorBugCandidate {
            json: estimator_bug,
        }),
    }
}

fn validated_attempt_fingerprint(
    attempt: &FalsifierAttempt,
) -> Result<AttemptFingerprint, FalsifyError> {
    validate_identifier("attempt ID", &attempt.attempt_id)?;
    validate_identifier("certificate class", &attempt.class)?;
    validate_identifier("regime", &attempt.regime)?;
    validate_identifier("falsifier name", &attempt.falsifier)?;
    validate_identifier("claim revision", &attempt.claim_revision)?;
    validate_identifier("artifact ID", &attempt.artifact_id)?;
    validate_positive_finite("compute seconds", attempt.compute_s)?;
    if let FalsifierOutcome::Discrepancy { detail } = &attempt.outcome
        && (detail.trim().is_empty() || detail.len() > MAX_DETAIL_BYTES)
    {
        return Err(FalsifyError::InvalidText {
            field: "falsifier detail",
            max_bytes: MAX_DETAIL_BYTES,
        });
    }

    Ok(AttemptFingerprint {
        class: attempt.class.clone(),
        regime: attempt.regime.clone(),
        falsifier: attempt.falsifier.clone(),
        claim_revision: attempt.claim_revision.clone(),
        artifact_id: attempt.artifact_id.clone(),
        seed: attempt.seed,
        compute_bits: attempt.compute_s.to_bits(),
        outcome: attempt.outcome.clone(),
    })
}

fn validate_attempt_catalog(
    registry: &FalsifierRegistry,
    attempt: &FalsifierAttempt,
) -> Result<(), FalsifyError> {
    if registry.falsifiers(&attempt.class).is_err() {
        return Err(FalsifyError::Unknown {
            class: attempt.class.clone(),
        });
    }
    if !registry.declares(&attempt.class, &attempt.falsifier) {
        return Err(FalsifyError::UnregisteredFalsifier {
            class: attempt.class.clone(),
            falsifier: attempt.falsifier.clone(),
        });
    }
    Ok(())
}

impl FalsifierHistory {
    /// An empty history.
    #[must_use]
    pub fn new() -> Self {
        FalsifierHistory::default()
    }

    /// Ingest one source-referencing diagnostic attempt. The caller-provided
    /// references are grammar-checked here; authenticated content/claim binding
    /// belongs to the package/checker successor. Identical retries are
    /// idempotent; reuse of an attempt ID for different content is rejected.
    /// On first ingestion, the registry check prevents typos/invented checker
    /// names from steering telemetry, but the declaration catalog still does
    /// not prove executable independence or confer release authority. An exact
    /// retry is compared with the immutable accepted fingerprint before the
    /// caller's current registry is consulted, so catalog replacement cannot
    /// break idempotent delivery or rewrite accepted history.
    pub fn record_attempt(
        &mut self,
        registry: &FalsifierRegistry,
        attempt: &FalsifierAttempt,
    ) -> Result<AttemptRecord, FalsifyError> {
        let fingerprint = validated_attempt_fingerprint(attempt)?;
        if let Some(existing) = self.attempts.get(&attempt.attempt_id) {
            if existing != &fingerprint {
                return Err(FalsifyError::AttemptCollision {
                    attempt_id: attempt.attempt_id.clone(),
                });
            }
            return Ok(candidate_record(attempt, true));
        }
        validate_attempt_catalog(registry, attempt)?;
        if self.attempts.len() >= MAX_ATTEMPTS {
            return Err(FalsifyError::ResourceLimit {
                field: "retained falsifier attempts",
                requested: self.attempts.len().saturating_add(1),
                cap: MAX_ATTEMPTS,
            });
        }

        let retained_attempt_bytes = attempt_retained_bytes(attempt)?;
        let next_attempt_bytes = self
            .attempt_bytes
            .checked_add(retained_attempt_bytes)
            .ok_or(FalsifyError::ResourceLimit {
                field: "aggregate falsifier attempt bytes",
                requested: usize::MAX,
                cap: MAX_ATTEMPT_BYTES_TOTAL,
            })?;
        if next_attempt_bytes > MAX_ATTEMPT_BYTES_TOTAL {
            return Err(FalsifyError::ResourceLimit {
                field: "aggregate falsifier attempt bytes",
                requested: next_attempt_bytes,
                cap: MAX_ATTEMPT_BYTES_TOTAL,
            });
        }

        let key = (attempt.class.clone(), attempt.regime.clone());
        let new_key_bytes = self.ensure_row_capacity(&key)?;
        let (passes, hits, compute) = self.rows.get(&key).copied().unwrap_or((0, 0, 0.0));
        let discrepancy = matches!(&attempt.outcome, FalsifierOutcome::Discrepancy { .. });
        let (next_passes, next_hits) = if discrepancy {
            (
                passes,
                hits.checked_add(1)
                    .ok_or_else(|| history_overflow(&attempt.class, &attempt.regime))?,
            )
        } else {
            (
                passes
                    .checked_add(1)
                    .ok_or_else(|| history_overflow(&attempt.class, &attempt.regime))?,
                hits,
            )
        };
        let _ = next_passes
            .checked_add(next_hits)
            .ok_or_else(|| history_overflow(&attempt.class, &attempt.regime))?;
        let next_compute = compute + attempt.compute_s;
        validate_nonnegative_finite("accumulated compute seconds", next_compute)?;

        let mut next_window = self
            .rent_windows
            .get(&attempt.class)
            .cloned()
            .unwrap_or_default();
        next_window.current_runs = next_window
            .current_runs
            .checked_add(1)
            .ok_or_else(|| history_overflow(&attempt.class, "<rent-window>"))?;
        if discrepancy {
            next_window.current_hits = next_window
                .current_hits
                .checked_add(1)
                .ok_or_else(|| history_overflow(&attempt.class, "<rent-window>"))?;
        }
        if next_window.current_runs == RENT_VOLUME {
            next_window.closed_window_hits.try_reserve(1).map_err(|_| {
                FalsifyError::ResourceExhausted {
                    operation: "falsifier rent window",
                }
            })?;
            next_window
                .closed_window_hits
                .push(next_window.current_hits);
            next_window.current_runs = 0;
            next_window.current_hits = 0;
        }

        let record = candidate_record(attempt, false);
        self.rows
            .insert(key, (next_passes, next_hits, next_compute));
        self.attempts
            .insert(attempt.attempt_id.clone(), fingerprint);
        self.rent_windows.insert(attempt.class.clone(), next_window);
        self.history_key_bytes = self
            .history_key_bytes
            .checked_add(new_key_bytes)
            .expect("row-key byte cap checked before insertion");
        self.attempt_bytes = next_attempt_bytes;
        Ok(record)
    }

    /// Conservative diagnostic doubt for one class/regime. The raw discrepancy
    /// rate is widened by a per-row, time-uniform union-bound Hoeffding radius
    /// at [`DOUBT_ALPHA`], then floored. It does not control multiplicity across
    /// rows and assumes an admitted stationary Bernoulli process; this mixed,
    /// unauthenticated legacy history cannot govern release.
    #[must_use]
    pub fn doubt(&self, class: &str, regime: &str) -> f64 {
        if validate_identifier("certificate class", class).is_err()
            || validate_identifier("regime", regime).is_err()
        {
            return DOUBT_COLD_START;
        }
        match self.rows.get(&(class.to_string(), regime.to_string())) {
            None => DOUBT_COLD_START,
            Some((passes, hits, _)) => {
                let total = passes + hits;
                if total == 0 {
                    return DOUBT_COLD_START;
                }
                #[allow(clippy::cast_precision_loss)]
                let n = total as f64;
                #[allow(clippy::cast_precision_loss)]
                let empirical = *hits as f64 / n;
                let alpha_n = DOUBT_ALPHA / (n * (n + 1.0));
                let radius = (-alpha_n.ln() / (2.0 * n)).sqrt();
                (empirical + radius).clamp(DOUBT_FLOOR, 1.0)
            }
        }
    }

    /// Preliminary YIELD telemetry for a class: (reported discrepancies,
    /// compute spent, runs) across all regimes.
    pub fn yield_of(&self, class: &str) -> Result<(u64, f64, u64), FalsifyError> {
        validate_identifier("certificate class", class)?;
        let mut hits = 0u64;
        let mut compute = 0.0f64;
        let mut runs = 0u64;
        for ((c, _), (p, h, s)) in &self.rows {
            if c == class {
                hits = hits
                    .checked_add(*h)
                    .ok_or_else(|| FalsifyError::HistoryOverflow {
                        class: class.to_string(),
                        regime: "<all>".to_string(),
                    })?;
                compute += *s;
                validate_nonnegative_finite("aggregated compute seconds", compute)?;
                let row_runs = p
                    .checked_add(*h)
                    .ok_or_else(|| FalsifyError::HistoryOverflow {
                        class: class.to_string(),
                        regime: "<all>".to_string(),
                    })?;
                runs = runs
                    .checked_add(row_runs)
                    .ok_or_else(|| FalsifyError::HistoryOverflow {
                        class: class.to_string(),
                        regime: "<all>".to_string(),
                    })?;
            }
        }
        Ok((hits, compute, runs))
    }

    /// The current budget-share multiplier for a class.
    #[must_use]
    pub fn share(&self, class: &str) -> f64 {
        self.share.get(class).copied().unwrap_or(1.0)
    }

    /// Preliminary class-level governance-window review. Windows close after
    /// exactly [`RENT_VOLUME`] newly ingested, ordered attempts, independent of
    /// when this method is called. Thus one call after 200 clean attempts has
    /// the same effect as calls after each block of 100. Per-falsifier policy
    /// receipts remain successor work.
    pub fn rent_review(&mut self) -> Result<Vec<(String, f64)>, FalsifyError> {
        let mut decisions = Vec::new();
        decisions
            .try_reserve_exact(self.rent_windows.len())
            .map_err(|_| FalsifyError::ResourceExhausted {
                operation: "falsifier rent review",
            })?;
        for (class, window) in &self.rent_windows {
            if window.closed_window_hits.is_empty() {
                continue;
            }
            let mut next = self.share(class);
            for hits in &window.closed_window_hits {
                if *hits == 0 {
                    next = (next * RENT_DECAY).max(RENT_SHARE_FLOOR);
                }
            }
            decisions.push((class.clone(), next));
        }
        let mut decayed = Vec::new();
        decayed.try_reserve_exact(decisions.len()).map_err(|_| {
            FalsifyError::ResourceExhausted {
                operation: "falsifier rent review output",
            }
        })?;
        for (class, next) in decisions {
            let current = self.share(&class);
            if next < current {
                self.share.insert(class.clone(), next);
                decayed.push((class.clone(), next));
            }
            self.rent_windows
                .get_mut(&class)
                .expect("decision came from rent window")
                .closed_window_hits
                .clear();
        }
        Ok(decayed)
    }

    fn ensure_row_capacity(&self, key: &(String, String)) -> Result<usize, FalsifyError> {
        if self.rows.contains_key(key) {
            return Ok(0);
        }
        if self.rows.len() >= MAX_HISTORY_ROWS {
            return Err(FalsifyError::TooManyHistoryRows {
                current: self.rows.len(),
                cap: MAX_HISTORY_ROWS,
            });
        }
        let key_bytes =
            key.0
                .len()
                .checked_add(key.1.len())
                .ok_or(FalsifyError::ResourceLimit {
                    field: "aggregate falsifier history key bytes",
                    requested: usize::MAX,
                    cap: MAX_HISTORY_KEY_BYTES_TOTAL,
                })?;
        let next_bytes =
            self.history_key_bytes
                .checked_add(key_bytes)
                .ok_or(FalsifyError::ResourceLimit {
                    field: "aggregate falsifier history key bytes",
                    requested: usize::MAX,
                    cap: MAX_HISTORY_KEY_BYTES_TOTAL,
                })?;
        if next_bytes > MAX_HISTORY_KEY_BYTES_TOTAL {
            return Err(FalsifyError::ResourceLimit {
                field: "aggregate falsifier history key bytes",
                requested: next_bytes,
                cap: MAX_HISTORY_KEY_BYTES_TOTAL,
            });
        }
        Ok(key_bytes)
    }
}

fn history_overflow(class: &str, regime: &str) -> FalsifyError {
    FalsifyError::HistoryOverflow {
        class: class.to_string(),
        regime: regime.to_string(),
    }
}

/// One claim awaiting falsification budget.
#[derive(Debug, Clone, PartialEq)]
pub struct ClaimContext {
    /// Certificate class.
    pub class: String,
    /// Regime key.
    pub regime: String,
    /// Downstream decision weight (DAG dependents; the ledger scores it).
    pub consequence: f64,
}

/// Allocate a diagnostic falsification budget across claims by consequence ×
/// doubt × class-review share. Inputs are finite, nonnegative, and bounded;
/// normalization is max-rescaled to avoid overflow. This does not authenticate
/// a falsifier or authorize release.
pub fn allocate_budget(
    total_budget_s: f64,
    claims: &[ClaimContext],
    history: &FalsifierHistory,
) -> Result<Vec<f64>, FalsifyError> {
    validate_nonnegative_finite("total falsification budget", total_budget_s)?;
    if claims.len() > MAX_CLAIMS_PER_ALLOCATION {
        return Err(FalsifyError::TooManyClaims {
            requested: claims.len(),
            cap: MAX_CLAIMS_PER_ALLOCATION,
        });
    }
    for claim in claims {
        validate_identifier("certificate class", &claim.class)?;
        validate_identifier("regime", &claim.regime)?;
        validate_nonnegative_finite("claim consequence", claim.consequence)?;
    }
    if claims.is_empty() {
        return Ok(Vec::new());
    }
    let mut weights = Vec::new();
    weights
        .try_reserve_exact(claims.len())
        .map_err(|_| FalsifyError::ResourceExhausted {
            operation: "falsification budget allocation",
        })?;
    if total_budget_s == 0.0 {
        weights.resize(claims.len(), 0.0);
        return Ok(weights);
    }
    let mut max_weight = 0.0f64;
    for claim in claims {
        let consequence = claim.consequence.max(CONSEQUENCE_FLOOR);
        let weight =
            consequence * history.doubt(&claim.class, &claim.regime) * history.share(&claim.class);
        validate_nonnegative_finite("falsification weight", weight)?;
        max_weight = max_weight.max(weight);
        weights.push(weight);
    }
    if max_weight == 0.0 {
        return Ok(vec![0.0; claims.len()]);
    }
    let scaled_total: f64 = weights.iter().map(|weight| weight / max_weight).sum();
    validate_nonnegative_finite("scaled falsification weight total", scaled_total)?;
    Ok(weights
        .iter()
        .map(|weight| total_budget_s * (weight / max_weight) / scaled_total)
        .collect())
}
