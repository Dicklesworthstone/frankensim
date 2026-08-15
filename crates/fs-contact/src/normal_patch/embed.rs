//! Transactional embedding of public normal-patch receipts into generic ports.
//!
//! This module deliberately consumes only [`NormalPatchRequest::evaluate`]. It
//! neither inspects nor reconstructs private constitutive-law state.

use core::fmt::Write as _;

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_tribo::{ExactlyOnceKeyError, ExactlyOnceKeyLedger};

use super::{
    ApplicabilityRatios, InputUncertainty, LineNormalPatchReceipt, NormalPatchError,
    NormalPatchReceipt, NormalPatchRequest, PointNormalPatchReceipt,
};

const EMBED_DOMAIN: &str = "org.frankensim.fs-contact.normal-patch.embed.v1";
const FRAME_TOLERANCE: f64 = 1.0e-12;

/// Explicit solver lane. Only a smooth, fixed branch can publish a residual.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationLane {
    /// Fixed active set: a constitutive tangent is valid for this solve step.
    SmoothFixed,
    /// Event handling is in progress; no smooth residual is published.
    Eventful,
}

/// Stable identifiers supplied by the consuming solver.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalPatchEmbedIdentity {
    pub solver_id: String,
    pub contact_id: String,
    pub feature_id: String,
    pub sample_id: String,
}

/// Solver-owned kinematics and iteration state in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalPatchKinematics {
    /// Signed normal separation; penetration is negative.
    pub declared_gap_m: f64,
    /// Required to equal `max(-declared_gap_m, 0)` on a fixed branch.
    pub approach_m: f64,
    pub approach_rate_m_per_s: f64,
    pub time_s: f64,
    pub step_s: f64,
    pub iteration: u64,
    /// Unit normal from the reacting body toward the acting body.
    pub normal: [f64; 3],
    /// Moment arm from the port origin to the contact feature.
    pub moment_arm_m: [f64; 3],
}

/// One solver-to-law request. The supplied law request is cloned and mapped;
/// the public model/source identity is preserved while its state identity is
/// made specific to this contact/time/iteration sample.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalPatchEmbedRequest {
    pub identity: NormalPatchEmbedIdentity,
    pub lane: IntegrationLane,
    pub converged: bool,
    pub kinematics: NormalPatchKinematics,
    pub law_request: NormalPatchRequest,
}

/// Checkpointable, deterministic ledger state. A state is immutable from the
/// caller's perspective: a successful transition returns its successor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalPatchEmbedState {
    anchor_time_bits: u64,
    max_forward_step_bits: u64,
    last_time_bits: u64,
    last_iteration: u64,
    work_ledger: NormalWorkKeyLedger,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum NormalWorkKeyLedger {
    Retained(ExactlyOnceKeyLedger<String>),
    StrictSequence(ExactlyOnceKeyLedger<u64>),
}

/// A rollback token made from a complete deterministic state snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalPatchEmbedCheckpoint {
    pub checkpoint_id: ContentHash,
    state: NormalPatchEmbedState,
}

/// Point-contact action/reaction wrench in SI `N` and `N m`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointNormalPort {
    pub action_force_n: [f64; 3],
    pub action_moment_n_m: [f64; 3],
    pub reaction_force_n: [f64; 3],
    pub reaction_moment_n_m: [f64; 3],
    pub residual_force_n: [f64; 3],
    pub tangent_n_per_m: f64,
    pub dissipated_power_w: f64,
    pub irreversible_work_j: f64,
}

/// Line-contact action/reaction wrench. Force is `N/m`; its moment about a
/// point is `N` because it is normalized by unit axial length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineNormalPort {
    pub action_force_n_per_m: [f64; 3],
    pub action_moment_n: [f64; 3],
    pub reaction_force_n_per_m: [f64; 3],
    pub reaction_moment_n: [f64; 3],
    pub residual_force_n_per_m: [f64; 3],
    pub tangent_pa: f64,
    pub dissipated_power_w_per_m: f64,
    pub irreversible_work_j_per_m: f64,
}

/// A typed generic port; point and line resultants cannot be mixed.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalPatchPort {
    Point(PointNormalPort),
    Line(LineNormalPort),
}

/// Successful, exactly-once transition with the public constitutive receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalPatchEmbedTransition {
    pub embedding_id: ContentHash,
    pub law_request_id: ContentHash,
    pub receipt_id: ContentHash,
    pub receipt: NormalPatchReceipt,
    pub port: NormalPatchPort,
    pub applicability: ApplicabilityRatios,
    pub uncertainty: InputUncertainty,
    pub next_state: NormalPatchEmbedState,
}

/// Refusal surface before a residual, tangent, or work record is published.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalPatchEmbedError {
    InvalidIdentity {
        field: &'static str,
    },
    InvalidKinematics {
        field: &'static str,
    },
    InvalidFrame,
    ApproachGapMismatch {
        declared_gap_m: f64,
        approach_m: f64,
    },
    EventfulLane,
    Nonconverged,
    StaleState {
        time_s: f64,
        iteration: u64,
    },
    FutureState {
        time_s: f64,
        maximum_time_s: f64,
    },
    DuplicateWorkKey {
        key: String,
    },
    OutOfSequenceWork {
        iteration: u64,
    },
    WorkKeyCapacityExceeded {
        maximum: usize,
    },
    Law(NormalPatchError),
}

impl core::fmt::Display for NormalPatchEmbedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidIdentity { field } => {
                write!(f, "nonblank embed identity required: {field}")
            }
            Self::InvalidKinematics { field } => {
                write!(f, "invalid normal-patch kinematics: {field}")
            }
            Self::InvalidFrame => write!(f, "normal-patch port requires a finite unit normal"),
            Self::ApproachGapMismatch {
                declared_gap_m,
                approach_m,
            } => write!(
                f,
                "approach {approach_m} does not equal active penetration from gap {declared_gap_m}"
            ),
            Self::EventfulLane => {
                write!(f, "eventful contact lane cannot publish a smooth tangent")
            }
            Self::Nonconverged => write!(f, "nonconverged contact iteration cannot publish a port"),
            Self::StaleState { time_s, iteration } => {
                write!(
                    f,
                    "stale contact sample at time {time_s}, iteration {iteration}"
                )
            }
            Self::FutureState {
                time_s,
                maximum_time_s,
            } => {
                write!(f, "future contact sample {time_s} exceeds {maximum_time_s}")
            }
            Self::DuplicateWorkKey { key } => {
                write!(f, "exactly-once work key already committed: {key}")
            }
            Self::OutOfSequenceWork { iteration } => {
                write!(
                    f,
                    "normal-patch work iteration is out of sequence: {iteration}"
                )
            }
            Self::WorkKeyCapacityExceeded { maximum } => {
                write!(f, "normal-patch work-key capacity {maximum} is exhausted")
            }
            Self::Law(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for NormalPatchEmbedError {}

impl From<NormalPatchError> for NormalPatchEmbedError {
    fn from(value: NormalPatchError) -> Self {
        Self::Law(value)
    }
}

impl NormalPatchEmbedState {
    /// Starts an ordered contact ledger at `anchor_time_s`; later samples may
    /// advance by at most `max_forward_step_s` before the solver checkpoints.
    pub fn new(anchor_time_s: f64, max_forward_step_s: f64) -> Result<Self, NormalPatchEmbedError> {
        if !anchor_time_s.is_finite() || anchor_time_s < 0.0 {
            return Err(NormalPatchEmbedError::InvalidKinematics {
                field: "anchor_time_s",
            });
        }
        if !max_forward_step_s.is_finite() || max_forward_step_s <= 0.0 {
            return Err(NormalPatchEmbedError::InvalidKinematics {
                field: "max_forward_step_s",
            });
        }
        Ok(Self {
            anchor_time_bits: anchor_time_s.to_bits(),
            max_forward_step_bits: max_forward_step_s.to_bits(),
            last_time_bits: anchor_time_s.to_bits(),
            last_iteration: 0,
            work_ledger: NormalWorkKeyLedger::Retained(
                ExactlyOnceKeyLedger::retained_set(usize::MAX)
                    .expect("usize::MAX is a nonzero work-key capacity"),
            ),
        })
    }

    /// Starts a fixed-memory ledger accepting iterations `1, 2, ...`.
    pub fn new_strict_sequence(
        anchor_time_s: f64,
        max_forward_step_s: f64,
        maximum_committed_work_keys: usize,
    ) -> Result<Self, NormalPatchEmbedError> {
        let mut state = Self::new(anchor_time_s, max_forward_step_s)?;
        state.work_ledger = NormalWorkKeyLedger::StrictSequence(
            ExactlyOnceKeyLedger::strict_sequence(1, maximum_committed_work_keys)
                .map_err(|error| normal_ledger_error(error, 0, String::new()))?,
        );
        Ok(state)
    }

    /// Captures a rollback token without changing this state.
    pub fn checkpoint(&self) -> NormalPatchEmbedCheckpoint {
        let checkpoint_id = hash_domain(EMBED_DOMAIN, self.canonical().as_bytes());
        NormalPatchEmbedCheckpoint {
            checkpoint_id,
            state: self.clone(),
        }
    }

    /// Restores the exact snapshot. This is the only rollback operation.
    pub fn rollback(checkpoint: &NormalPatchEmbedCheckpoint) -> Self {
        checkpoint.state.clone()
    }

    fn anchor_time_s(&self) -> f64 {
        f64::from_bits(self.anchor_time_bits)
    }

    fn max_forward_step_s(&self) -> f64 {
        f64::from_bits(self.max_forward_step_bits)
    }

    fn last_time_s(&self) -> f64 {
        f64::from_bits(self.last_time_bits)
    }

    fn canonical(&self) -> String {
        let keys = match &self.work_ledger {
            NormalWorkKeyLedger::Retained(ledger) => ledger
                .retained_keys()
                .expect("retained ledger exposes its exact keys")
                .iter()
                .cloned()
                .collect::<Vec<_>>()
                .join(","),
            NormalWorkKeyLedger::StrictSequence(ledger) => format!(
                "strict:{}:{}",
                ledger.committed_count(),
                ledger
                    .strict_next_key()
                    .expect("strict ledger exposes its next key")
            ),
        };
        format!(
            "v1|{:.17e}|{:.17e}|{:.17e}|{}|{keys}",
            self.anchor_time_s(),
            self.max_forward_step_s(),
            self.last_time_s(),
            self.last_iteration,
        )
    }
}

impl NormalPatchEmbedRequest {
    /// Evaluates and stages one valid fixed-branch contact sample. A refusal
    /// returns no port and leaves the caller's state untouched.
    pub fn evaluate(
        &self,
        state: &NormalPatchEmbedState,
    ) -> Result<NormalPatchEmbedTransition, NormalPatchEmbedError> {
        let work_key = self.work_key();
        let next_work_ledger = match &state.work_ledger {
            NormalWorkKeyLedger::Retained(ledger) => {
                NormalWorkKeyLedger::Retained(ledger.advance(&work_key, None).map_err(|error| {
                    normal_ledger_error(error, self.kinematics.iteration, work_key.clone())
                })?)
            }
            NormalWorkKeyLedger::StrictSequence(ledger) => {
                let successor = self
                    .kinematics
                    .iteration
                    .checked_add(1)
                    .ok_or(NormalPatchEmbedError::InvalidKinematics { field: "iteration" })?;
                NormalWorkKeyLedger::StrictSequence(
                    ledger
                        .advance(&self.kinematics.iteration, Some(successor))
                        .map_err(|error| {
                            normal_ledger_error(error, self.kinematics.iteration, work_key.clone())
                        })?,
                )
            }
        };
        self.validate(state)?;
        let mut query = self.law_request.clone();
        query.indentation_m = self.kinematics.approach_m;
        query.indentation_rate_m_per_s = self.kinematics.approach_rate_m_per_s;
        query.step_s = self.kinematics.step_s;
        query.identity.state_id = self.mapped_law_state_id();
        let receipt = query.evaluate()?;
        let (port, applicability, uncertainty) = self.port_from_receipt(&receipt);
        let embedding_id = hash_embedding_id(
            &work_key,
            receipt.request_id(),
            receipt.receipt_id(),
            &query.identity.model_id,
        );
        let mut next_state = state.clone();
        next_state.last_time_bits = self.kinematics.time_s.to_bits();
        next_state.last_iteration = self.kinematics.iteration;
        next_state.work_ledger = next_work_ledger;
        Ok(NormalPatchEmbedTransition {
            embedding_id,
            law_request_id: receipt.request_id(),
            receipt_id: receipt.receipt_id(),
            receipt,
            port,
            applicability,
            uncertainty,
            next_state,
        })
    }

    fn validate(&self, state: &NormalPatchEmbedState) -> Result<(), NormalPatchEmbedError> {
        for (value, field) in [
            (&self.identity.solver_id, "solver_id"),
            (&self.identity.contact_id, "contact_id"),
            (&self.identity.feature_id, "feature_id"),
            (&self.identity.sample_id, "sample_id"),
        ] {
            if value.trim().is_empty() {
                return Err(NormalPatchEmbedError::InvalidIdentity { field });
            }
        }
        if self.lane == IntegrationLane::Eventful {
            return Err(NormalPatchEmbedError::EventfulLane);
        }
        if !self.converged {
            return Err(NormalPatchEmbedError::Nonconverged);
        }
        if !self.kinematics.declared_gap_m.is_finite() {
            return Err(NormalPatchEmbedError::InvalidKinematics {
                field: "declared_gap_m",
            });
        }
        for (value, field) in [
            (self.kinematics.approach_m, "approach_m"),
            (self.kinematics.time_s, "time_s"),
            (self.kinematics.step_s, "step_s"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(NormalPatchEmbedError::InvalidKinematics { field });
            }
        }
        if !self.kinematics.approach_rate_m_per_s.is_finite() {
            return Err(NormalPatchEmbedError::InvalidKinematics {
                field: "approach_rate_m_per_s",
            });
        }
        if self.kinematics.step_s <= 0.0 {
            return Err(NormalPatchEmbedError::InvalidKinematics { field: "step_s" });
        }
        let expected_approach = (-self.kinematics.declared_gap_m).max(0.0);
        if (self.kinematics.approach_m - expected_approach).abs()
            > FRAME_TOLERANCE * expected_approach.max(1.0)
        {
            return Err(NormalPatchEmbedError::ApproachGapMismatch {
                declared_gap_m: self.kinematics.declared_gap_m,
                approach_m: self.kinematics.approach_m,
            });
        }
        let normal_norm_squared = dot(self.kinematics.normal, self.kinematics.normal);
        if !normal_norm_squared.is_finite()
            || (normal_norm_squared.sqrt() - 1.0).abs() > FRAME_TOLERANCE
            || !self
                .kinematics
                .moment_arm_m
                .iter()
                .all(|value| value.is_finite())
        {
            return Err(NormalPatchEmbedError::InvalidFrame);
        }
        let last_time_s = state.last_time_s();
        if self.kinematics.time_s < last_time_s
            || (self.kinematics.time_s == last_time_s
                && self.kinematics.iteration <= state.last_iteration)
        {
            return Err(NormalPatchEmbedError::StaleState {
                time_s: self.kinematics.time_s,
                iteration: self.kinematics.iteration,
            });
        }
        let maximum_time_s = last_time_s + state.max_forward_step_s();
        if self.kinematics.time_s > maximum_time_s {
            return Err(NormalPatchEmbedError::FutureState {
                time_s: self.kinematics.time_s,
                maximum_time_s,
            });
        }
        Ok(())
    }

    fn port_from_receipt(
        &self,
        receipt: &NormalPatchReceipt,
    ) -> (NormalPatchPort, ApplicabilityRatios, InputUncertainty) {
        match receipt {
            NormalPatchReceipt::Point(receipt) => (
                NormalPatchPort::Point(point_port(
                    receipt,
                    self.kinematics.normal,
                    self.kinematics.moment_arm_m,
                )),
                receipt.ratios,
                receipt.uncertainty,
            ),
            NormalPatchReceipt::Line(receipt) => (
                NormalPatchPort::Line(line_port(
                    receipt,
                    self.kinematics.normal,
                    self.kinematics.moment_arm_m,
                )),
                receipt.ratios,
                receipt.uncertainty,
            ),
        }
    }

    fn mapped_law_state_id(&self) -> String {
        format!(
            "embed/v1/{}/{}/{}/{:016x}/{}",
            self.identity.solver_id,
            self.identity.contact_id,
            self.identity.feature_id,
            self.kinematics.time_s.to_bits(),
            self.kinematics.iteration,
        )
    }

    fn work_key(&self) -> String {
        format!(
            "v1|{}|{}|{}|{}|{:016x}|{}",
            self.identity.solver_id,
            self.identity.contact_id,
            self.identity.feature_id,
            self.identity.sample_id,
            self.kinematics.time_s.to_bits(),
            self.kinematics.iteration,
        )
    }
}

fn hash_embedding_id(
    work_key: &str,
    request_id: ContentHash,
    receipt_id: ContentHash,
    model_id: &str,
) -> ContentHash {
    let mut hasher = DomainHasher::new(EMBED_DOMAIN);
    write!(
        &mut hasher,
        "{work_key}|{request_id}|{receipt_id}|{model_id}"
    )
    .expect("writing to DomainHasher cannot fail");
    hasher.finalize()
}

fn normal_ledger_error(
    error: ExactlyOnceKeyError,
    iteration: u64,
    work_key: String,
) -> NormalPatchEmbedError {
    match error {
        ExactlyOnceKeyError::Duplicate => NormalPatchEmbedError::DuplicateWorkKey { key: work_key },
        ExactlyOnceKeyError::OutOfSequence => {
            NormalPatchEmbedError::OutOfSequenceWork { iteration }
        }
        ExactlyOnceKeyError::CapacityExceeded { maximum } => {
            NormalPatchEmbedError::WorkKeyCapacityExceeded { maximum }
        }
        ExactlyOnceKeyError::ZeroCapacity | ExactlyOnceKeyError::MissingSuccessor => {
            NormalPatchEmbedError::InvalidKinematics {
                field: "exactly-once work ledger",
            }
        }
    }
}

fn point_port(
    receipt: &PointNormalPatchReceipt,
    normal: [f64; 3],
    arm: [f64; 3],
) -> PointNormalPort {
    let action_force_n = scale(normal, receipt.normal_force_n);
    let action_moment_n_m = cross(arm, action_force_n);
    PointNormalPort {
        action_force_n,
        action_moment_n_m,
        reaction_force_n: scale(action_force_n, -1.0),
        reaction_moment_n_m: scale(action_moment_n_m, -1.0),
        residual_force_n: [0.0; 3],
        tangent_n_per_m: receipt.tangent_n_per_m,
        dissipated_power_w: receipt.dissipated_power_w,
        irreversible_work_j: receipt.irreversible_work_j,
    }
}

fn line_port(receipt: &LineNormalPatchReceipt, normal: [f64; 3], arm: [f64; 3]) -> LineNormalPort {
    let action_force_n_per_m = scale(normal, receipt.normal_line_load_n_per_m);
    let action_moment_n = cross(arm, action_force_n_per_m);
    LineNormalPort {
        action_force_n_per_m,
        action_moment_n,
        reaction_force_n_per_m: scale(action_force_n_per_m, -1.0),
        reaction_moment_n: scale(action_moment_n, -1.0),
        residual_force_n_per_m: [0.0; 3],
        tangent_pa: receipt.tangent_pa,
        dissipated_power_w_per_m: receipt.dissipated_power_w_per_m,
        irreversible_work_j_per_m: receipt.irreversible_work_j_per_m,
    }
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::hash_embedding_id;
    use fs_blake3::ContentHash;

    #[test]
    fn g0_embedding_id_matches_pre_streaming_golden() {
        assert_eq!(
            hash_embedding_id(
                "work-key",
                ContentHash([0x11; 32]),
                ContentHash([0x22; 32]),
                "model",
            )
            .to_hex(),
            "c76bf66353e457faade444ddbca8e0393d7f42854b6b0e5bf0ae92146a46b9b4"
        );
    }
}
