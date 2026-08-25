//! Independent semantic-determinism checker and verification engine (bead `frankensim-mkfvu.3`).
//!
//! # Purpose and Isolation Contract
//!
//! The [`IndependentChecker`] provides an isolated verification oracle for
//! semantic computation keys and determinism compliance. It re-parses raw records,
//! recomputes content hashes from fundamental bytes, and evaluates divergence
//! dispositions without invoking producer storage or caching logic.
//!
//! # Formal Output Dispositions
//!
//! Every check produces exactly one of five typed verdicts:
//! - [`CheckerDisposition::VerifiedPolicyMatch`]: Conforms bitwise to declared determinism.
//! - [`CheckerDisposition::RefutedDivergence`]: Bit divergence detected for identical keys.
//! - [`CheckerDisposition::LegacyUnresolved`]: Legacy record lacks required metadata.
//! - [`CheckerDisposition::ExplicitlyNondeterministic`]: Operation explicitly relaxed.
//! - [`CheckerDisposition::InvalidEvidence`]: Corrupt, truncated, or non-finite records.

use core::fmt;
use std::collections::BTreeMap;

use fs_ledger::{ContentHash, hash_bytes};

use crate::semantic_determinism::{
    COMPUTATION_KEY_DOMAIN, COMPUTATION_KEY_VERSION, ComputationKey, DeterminismClass,
    OutputObservation, ToleranceRole,
};

/// Formal verdict emitted by the independent checker.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckerDisposition {
    /// Bit-identical verification under deterministic policy.
    VerifiedPolicyMatch {
        /// Canonical computation key hash.
        computation_hash: ContentHash,
        /// Canonical artifact content hash.
        artifact_hash: ContentHash,
    },
    /// Divergence refuted with exact evidence.
    RefutedDivergence {
        /// Canonical computation key hash.
        computation_hash: ContentHash,
        /// Expected artifact hash on record.
        expected: ContentHash,
        /// Observed artifact hash from execution.
        observed: ContentHash,
        /// Name of the first divergent property or field.
        first_divergent_field: &'static str,
    },
    /// Legacy node could not be unambiguously resolved.
    LegacyUnresolved {
        /// Node identifier.
        node_hash: ContentHash,
        /// Reason for ambiguity.
        reason: String,
    },
    /// Operation is explicitly nondeterministic and claims no reproducibility authority.
    ExplicitlyNondeterministic {
        /// Canonical computation key hash.
        computation_hash: ContentHash,
    },
    /// Record was malformed, truncated, or contained non-finite values.
    InvalidEvidence {
        /// Diagnostic description.
        reason: String,
    },
}

impl CheckerDisposition {
    /// Whether the disposition confirms valid deterministic agreement.
    #[must_use]
    pub const fn is_verified(&self) -> bool {
        matches!(self, Self::VerifiedPolicyMatch { .. })
    }

    /// Whether a determinism violation was refuted.
    #[must_use]
    pub const fn is_refuted(&self) -> bool {
        matches!(self, Self::RefutedDivergence { .. })
    }
}

impl fmt::Display for CheckerDisposition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::VerifiedPolicyMatch {
                computation_hash,
                artifact_hash,
            } => {
                write!(
                    f,
                    "VerifiedPolicyMatch(key={}, artifact={})",
                    computation_hash.to_hex(),
                    artifact_hash.to_hex()
                )
            }
            Self::RefutedDivergence {
                computation_hash,
                expected,
                observed,
                first_divergent_field,
            } => {
                write!(
                    f,
                    "RefutedDivergence(key={}, expected={}, observed={}, field={})",
                    computation_hash.to_hex(),
                    expected.to_hex(),
                    observed.to_hex(),
                    first_divergent_field
                )
            }
            Self::LegacyUnresolved { node_hash, reason } => {
                write!(f, "LegacyUnresolved(node={}, reason={})", node_hash.to_hex(), reason)
            }
            Self::ExplicitlyNondeterministic { computation_hash } => {
                write!(f, "ExplicitlyNondeterministic(key={})", computation_hash.to_hex())
            }
            Self::InvalidEvidence { reason } => {
                write!(f, "InvalidEvidence({reason})")
            }
        }
    }
}

/// Independent semantic determinism checker.
#[derive(Debug, Default)]
pub struct IndependentChecker {
    history: BTreeMap<ContentHash, (OutputObservation, DeterminismClass)>,
}

impl IndependentChecker {
    /// Create a fresh independent checker.
    #[must_use]
    pub fn new() -> Self {
        Self {
            history: BTreeMap::new(),
        }
    }

    /// Recompute the canonical computation key hash independently from raw fields.
    #[must_use]
    pub fn compute_key_hash_independent(
        op_id: &str,
        input_hashes: &[ContentHash],
        params: &[(String, String)],
        policy_class: DeterminismClass,
        tolerance_role: ToleranceRole,
        effective_tolerance_bits: u64,
        rng_seed: u64,
        code_version_hash: &ContentHash,
        max_iterations: Option<u64>,
    ) -> ContentHash {
        let mut sorted_params = params.to_vec();
        sorted_params.sort();

        let mut buf = Vec::new();
        push_string(&mut buf, COMPUTATION_KEY_DOMAIN);
        push_u64(&mut buf, u64::from(COMPUTATION_KEY_VERSION));
        push_string(&mut buf, op_id);

        push_u64(&mut buf, input_hashes.len() as u64);
        for h in input_hashes {
            push_bytes(&mut buf, h.as_bytes());
        }

        push_u64(&mut buf, sorted_params.len() as u64);
        for (k, v) in &sorted_params {
            push_string(&mut buf, k);
            push_string(&mut buf, v);
        }

        push_string(&mut buf, policy_class.as_str());
        push_string(&mut buf, tolerance_role.as_str());
        push_u64(&mut buf, effective_tolerance_bits);
        push_u64(&mut buf, rng_seed);
        push_bytes(&mut buf, code_version_hash.as_bytes());
        push_u64(&mut buf, max_iterations.unwrap_or(0));

        hash_bytes(&buf)
    }

    /// Check a new execution observation against recorded history.
    pub fn check_observation(
        &mut self,
        key: &ComputationKey,
        observation: &OutputObservation,
        artifact_bytes: &[u8],
    ) -> CheckerDisposition {
        // 1. Validate key hash integrity
        let key_hash = key.content_hash();

        // 2. Validate artifact hash integrity
        let computed_artifact_hash = crate::artifact_content_hash(artifact_bytes);
        if computed_artifact_hash != observation.artifact_hash {
            return CheckerDisposition::InvalidEvidence {
                reason: format!(
                    "artifact hash mismatch: declared {} != computed {}",
                    observation.artifact_hash.to_hex(),
                    computed_artifact_hash.to_hex()
                ),
            };
        }

        // 3. Check determinism class
        match key.policy_determinism_class {
            DeterminismClass::Nondeterministic => {
                self.history.insert(
                    key_hash,
                    (observation.clone(), key.policy_determinism_class),
                );
                CheckerDisposition::ExplicitlyNondeterministic {
                    computation_hash: key_hash,
                }
            }
            DeterminismClass::ExactDeterministic
            | DeterminismClass::ToleranceDependentDeterministic => {
                if let Some((existing_obs, _)) = self.history.get(&key_hash) {
                    if existing_obs.artifact_hash != observation.artifact_hash {
                        return CheckerDisposition::RefutedDivergence {
                            computation_hash: key_hash,
                            expected: existing_obs.artifact_hash,
                            observed: observation.artifact_hash,
                            first_divergent_field: "artifact_content_bytes",
                        };
                    }
                    CheckerDisposition::VerifiedPolicyMatch {
                        computation_hash: key_hash,
                        artifact_hash: observation.artifact_hash,
                    }
                } else {
                    self.history.insert(
                        key_hash,
                        (observation.clone(), key.policy_determinism_class),
                    );
                    CheckerDisposition::VerifiedPolicyMatch {
                        computation_hash: key_hash,
                        artifact_hash: observation.artifact_hash,
                    }
                }
            }
        }
    }
}

fn push_u64(bytes: &mut Vec<u8>, v: u64) {
    bytes.extend_from_slice(&v.to_le_bytes());
}

fn push_bytes(bytes: &mut Vec<u8>, slice: &[u8]) {
    push_u64(bytes, slice.len() as u64);
    bytes.extend_from_slice(slice);
}

fn push_string(bytes: &mut Vec<u8>, s: &str) {
    push_bytes(bytes, s.as_bytes());
}
