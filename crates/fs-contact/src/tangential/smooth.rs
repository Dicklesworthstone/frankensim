//! Smooth-solver adapter for the generic finite-patch partial-slip law.
//!
//! This module owns solver transaction semantics only.  It delegates every
//! constitutive update, capacity, storage, and heat calculation to
//! [`fs_tribo::partial_slip`]; it neither reimplements nor promotes that law.
//! The regularization is a declared, bounded kinematic preprocessing map.  It
//! is C-infinity but does **not** remove the constituent law's branch events,
//! so derivatives are available only when all finite-difference probes retain
//! one explicitly reported branch.

use core::fmt;

use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchAuthority, NormalPatchView, PartialSlipCheckpoint,
    PartialSlipError, PartialSlipInterface, PartialSlipKinematics, PartialSlipLaw,
    PartialSlipStateKind, TangentFrame, TangentialWrench,
};

/// Stable identity of this solver embedding.
pub const SMOOTH_TANGENTIAL_ADAPTER_MODEL_ID: &str = "fs-contact/smooth-partial-slip-adapter-v1";

/// Explicit acceptance policy for caller-supplied authority labels.
///
/// This is a refusal boundary, not a ranking of authority kinds.  An accepted
/// input keeps its original label in every receipt; the adapter never upgrades
/// it merely because a smooth solver consumed it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmoothAuthorityPolicy {
    /// Admit caller-declared normal/interface receipts.
    pub allow_caller_declared: bool,
    /// Admit synthetic fixtures.  Production callers normally set this false.
    pub allow_synthetic_fixture: bool,
    /// Admit estimate-only normal/interface receipts.
    pub allow_estimated: bool,
}

impl SmoothAuthorityPolicy {
    /// Policy suitable for numerical fixtures only.
    #[must_use]
    pub const fn test_only() -> Self {
        Self {
            allow_caller_declared: true,
            allow_synthetic_fixture: true,
            allow_estimated: true,
        }
    }

    fn admits(self, authority: NormalPatchAuthority) -> bool {
        match authority {
            NormalPatchAuthority::CallerDeclared => self.allow_caller_declared,
            NormalPatchAuthority::SyntheticFixture => self.allow_synthetic_fixture,
            NormalPatchAuthority::Estimated => self.allow_estimated,
        }
    }
}

/// Units-explicit smoothing and derivative-probe controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothRegularization {
    /// Positive dimensionless deadzone scale for the tangent-plane
    /// creepage MAGNITUDE. The smoothing is isotropic: it rescales the
    /// creepage vector by the regularized magnitude, never the frame
    /// components separately (a componentwise deadzone breaks the SO(2)
    /// frame covariance the receipt claims — executed: the
    /// frame-rotation conformance test measured a 7e-5 relative
    /// world-force violation before this was fixed).
    pub creepage_scale: f64,
    /// Positive scale for relative torsional spin, in rad/s.
    pub torsional_spin_scale_rad_per_s: f64,
    /// Positive central-difference creepage probe, dimensionless.
    pub tangent_probe_creepage: f64,
    /// Positive central-difference spin probe, in rad/s.
    pub tangent_probe_spin_rad_per_s: f64,
}

impl SmoothRegularization {
    fn validate(self) -> Result<(), SmoothTangentialError> {
        positive(self.creepage_scale, "creepage_scale")?;
        positive(
            self.torsional_spin_scale_rad_per_s,
            "torsional_spin_scale_rad_per_s",
        )?;
        positive(self.tangent_probe_creepage, "tangent_probe_creepage")?;
        positive(
            self.tangent_probe_spin_rad_per_s,
            "tangent_probe_spin_rad_per_s",
        )
    }
}

/// Named smooth-lane identity, authority policy, and regularization receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothTangentialAdapter {
    adapter_id: String,
    source_id: String,
    regularization: SmoothRegularization,
    authority_policy: SmoothAuthorityPolicy,
}

impl SmoothTangentialAdapter {
    /// Creates a deterministic smooth-lane adapter.
    pub fn new(
        adapter_id: impl Into<String>,
        source_id: impl Into<String>,
        regularization: SmoothRegularization,
        authority_policy: SmoothAuthorityPolicy,
    ) -> Result<Self, SmoothTangentialError> {
        let value = Self {
            adapter_id: adapter_id.into(),
            source_id: source_id.into(),
            regularization,
            authority_policy,
        };
        value.validate()?;
        Ok(value)
    }

    /// Stable adapter identity.
    #[must_use]
    pub fn adapter_id(&self) -> &str {
        &self.adapter_id
    }

    /// Caller source identity for this adapter configuration.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Declared regularization controls.
    #[must_use]
    pub const fn regularization(&self) -> SmoothRegularization {
        self.regularization
    }

    /// Explicit input-authority admission policy.
    #[must_use]
    pub const fn authority_policy(&self) -> SmoothAuthorityPolicy {
        self.authority_policy
    }

    /// Creates a zero-history checkpoint bound to all upstream identities.
    pub fn initial_state(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        max_committed_work_keys: usize,
    ) -> Result<SmoothTangentialState, SmoothTangentialError> {
        self.validate_inputs(patch, interface)?;
        if max_committed_work_keys == 0 {
            return Err(SmoothTangentialError::InvalidInput {
                field: "max_committed_work_keys",
            });
        }
        let law_checkpoint = PartialSlipCheckpoint::new(
            patch.clone(),
            interface.clone(),
            law.clone(),
            fs_tribo::partial_slip::PartialSlipState::zero(),
        )?;
        Ok(SmoothTangentialState {
            checkpoint: SmoothTangentialCheckpoint {
                adapter_id: self.adapter_id.clone(),
                adapter_source_id: self.source_id.clone(),
                regularization: self.regularization,
                authority_policy: self.authority_policy,
                law_checkpoint,
                committed_version: 0,
                committed_work_keys: Vec::new(),
                max_committed_work_keys,
            },
        })
    }

    /// Prepares one candidate without mutating the supplied state.
    ///
    /// The caller may discard it through [`Self::rollback`] and retry with the
    /// same work key.  Work is consumed only by [`Self::commit`].
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        state: &SmoothTangentialState,
        query: &SmoothTangentialQuery,
    ) -> Result<SmoothTangentialReceipt, SmoothTangentialError> {
        self.validate_state(law, patch, interface, state)?;
        self.validate_query(state, query, true)?;
        let regularized_kinematics = self.regularize_kinematics(query.kinematics)?;
        let law_state =
            law.restore_checkpoint(patch, interface, &state.checkpoint.law_checkpoint)?;
        let step = law.advance(
            patch,
            interface,
            query.frame,
            regularized_kinematics,
            &query.work_ownership,
            &law_state,
        )?;
        let action_reaction = ActionReactionWrench::from_action(step.wrench);
        let next_checkpoint = SmoothTangentialCheckpoint {
            adapter_id: self.adapter_id.clone(),
            adapter_source_id: self.source_id.clone(),
            regularization: self.regularization,
            authority_policy: self.authority_policy,
            law_checkpoint: step.checkpoint.clone(),
            committed_version: state.checkpoint.committed_version.checked_add(1).ok_or(
                SmoothTangentialError::InvalidDerived {
                    field: "committed_version overflow",
                },
            )?,
            committed_work_keys: state.checkpoint.committed_work_keys.clone(),
            max_committed_work_keys: state.checkpoint.max_committed_work_keys,
        };
        Ok(SmoothTangentialReceipt {
            query_id: query.query_id.clone(),
            parent_version: state.checkpoint.committed_version,
            branch: step.state,
            regularized_kinematics,
            action_reaction,
            residual: SmoothResidual {
                tangent_force_n: step.tangent_force_n,
                torsional_moment_nm: step.torsional_moment_nm,
                endpoint_relative_power_w: step.generalized_work.endpoint_body_power_w,
                reconstructed_body_power_w: step.generalized_work.reconstructed_body_power_w,
            },
            partial_slip_step: step,
            prior_state: state.clone(),
            next_checkpoint,
        })
    }

    /// Returns a central finite-difference tangent only within one law branch.
    ///
    /// This deliberately refuses if the base or either probe has a different
    /// [`PartialSlipStateKind`].  It is a local numerical derivative, not a
    /// global differentiability claim across stick/slip events.
    #[allow(clippy::too_many_arguments)]
    pub fn fixed_branch_tangent(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        state: &SmoothTangentialState,
        query: &SmoothTangentialQuery,
    ) -> Result<FixedBranchTangent, SmoothTangentialError> {
        let base = self.prepare_preview(law, patch, interface, state, query)?;
        let mut derivative = [[0.0; 3]; 3];
        for column in 0..3 {
            let step = if column < 2 {
                self.regularization.tangent_probe_creepage
            } else {
                self.regularization.tangent_probe_spin_rad_per_s
            };
            let positive =
                self.preview_perturbed(law, patch, interface, state, query, column, step)?;
            let negative =
                self.preview_perturbed(law, patch, interface, state, query, column, -step)?;
            if positive.branch != base.branch || negative.branch != base.branch {
                return Err(SmoothTangentialError::NoDerivativeOnBranchChange {
                    base: base.branch,
                    positive: positive.branch,
                    negative: negative.branch,
                });
            }
            let plus = positive.residual.components();
            let minus = negative.residual.components();
            for row in 0..3 {
                derivative[row][column] = (plus[row] - minus[row]) / (2.0 * step);
                finite(derivative[row][column], "fixed_branch_tangent")?;
            }
        }
        Ok(FixedBranchTangent {
            branch: base.branch,
            derivative,
            input_units: [
                "dimensionless creepage",
                "dimensionless creepage",
                "rad/s torsional spin",
            ],
            output_units: ["N", "N", "N m"],
        })
    }

    /// Commits a prepared receipt exactly once, producing an immutable new state.
    pub fn commit(
        &self,
        state: &SmoothTangentialState,
        receipt: &SmoothTangentialReceipt,
    ) -> Result<SmoothTangentialState, SmoothTangentialError> {
        self.validate()?;
        if receipt.prior_state != *state {
            return match receipt
                .parent_version
                .cmp(&state.checkpoint.committed_version)
            {
                core::cmp::Ordering::Less => Err(SmoothTangentialError::StaleState {
                    expected: state.checkpoint.committed_version,
                    observed: receipt.parent_version,
                }),
                core::cmp::Ordering::Greater => Err(SmoothTangentialError::FutureState {
                    expected: state.checkpoint.committed_version,
                    observed: receipt.parent_version,
                }),
                core::cmp::Ordering::Equal => Err(SmoothTangentialError::ReceiptDoesNotMatchState),
            };
        }
        if state
            .checkpoint
            .committed_work_keys
            .iter()
            .any(|key| key == &receipt.partial_slip_step.generalized_work.ownership)
        {
            return Err(SmoothTangentialError::DuplicateWorkKey);
        }
        if state.checkpoint.committed_work_keys.len() >= state.checkpoint.max_committed_work_keys {
            return Err(SmoothTangentialError::WorkKeyCapacityExceeded {
                max: state.checkpoint.max_committed_work_keys,
            });
        }
        let mut next = receipt.next_checkpoint.clone();
        next.committed_work_keys
            .push(receipt.partial_slip_step.generalized_work.ownership.clone());
        Ok(SmoothTangentialState { checkpoint: next })
    }

    /// Discards an uncommitted candidate and returns its exact prior state.
    pub fn rollback(
        &self,
        receipt: &SmoothTangentialReceipt,
    ) -> Result<SmoothTangentialState, SmoothTangentialError> {
        self.validate()?;
        if receipt.next_checkpoint.adapter_id != self.adapter_id
            || receipt.next_checkpoint.adapter_source_id != self.source_id
            || receipt.next_checkpoint.regularization != self.regularization
            || receipt.next_checkpoint.authority_policy != self.authority_policy
        {
            return Err(SmoothTangentialError::CheckpointIdentityMismatch {
                field: "adapter receipt",
            });
        }
        Ok(receipt.prior_state.clone())
    }

    /// Restores a checkpoint only if this adapter and all upstream law inputs match.
    #[allow(clippy::too_many_arguments)]
    pub fn restore_checkpoint(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        checkpoint: SmoothTangentialCheckpoint,
    ) -> Result<SmoothTangentialState, SmoothTangentialError> {
        self.validate_inputs(patch, interface)?;
        if checkpoint.adapter_id != self.adapter_id {
            return Err(SmoothTangentialError::CheckpointIdentityMismatch {
                field: "adapter_id",
            });
        }
        if checkpoint.adapter_source_id != self.source_id {
            return Err(SmoothTangentialError::CheckpointIdentityMismatch {
                field: "adapter_source_id",
            });
        }
        if checkpoint.regularization != self.regularization {
            return Err(SmoothTangentialError::CheckpointIdentityMismatch {
                field: "regularization",
            });
        }
        if checkpoint.authority_policy != self.authority_policy {
            return Err(SmoothTangentialError::CheckpointIdentityMismatch {
                field: "authority_policy",
            });
        }
        if checkpoint.max_committed_work_keys == 0
            || checkpoint.committed_work_keys.len() > checkpoint.max_committed_work_keys
        {
            return Err(SmoothTangentialError::InvalidInput {
                field: "checkpoint work-key budget",
            });
        }
        let _ = law.restore_checkpoint(patch, interface, &checkpoint.law_checkpoint)?;
        Ok(SmoothTangentialState { checkpoint })
    }

    fn prepare_preview(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        state: &SmoothTangentialState,
        query: &SmoothTangentialQuery,
    ) -> Result<SmoothTangentialReceipt, SmoothTangentialError> {
        self.validate_state(law, patch, interface, state)?;
        self.validate_query(state, query, false)?;
        let regularized_kinematics = self.regularize_kinematics(query.kinematics)?;
        let law_state =
            law.restore_checkpoint(patch, interface, &state.checkpoint.law_checkpoint)?;
        let step = law.advance(
            patch,
            interface,
            query.frame,
            regularized_kinematics,
            &query.work_ownership,
            &law_state,
        )?;
        Ok(SmoothTangentialReceipt {
            query_id: query.query_id.clone(),
            parent_version: state.checkpoint.committed_version,
            branch: step.state,
            regularized_kinematics,
            action_reaction: ActionReactionWrench::from_action(step.wrench),
            residual: SmoothResidual {
                tangent_force_n: step.tangent_force_n,
                torsional_moment_nm: step.torsional_moment_nm,
                endpoint_relative_power_w: step.generalized_work.endpoint_body_power_w,
                reconstructed_body_power_w: step.generalized_work.reconstructed_body_power_w,
            },
            partial_slip_step: step,
            prior_state: state.clone(),
            next_checkpoint: state.checkpoint.clone(),
        })
    }

    fn preview_perturbed(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        state: &SmoothTangentialState,
        query: &SmoothTangentialQuery,
        coordinate: usize,
        delta: f64,
    ) -> Result<SmoothTangentialReceipt, SmoothTangentialError> {
        let mut perturbed = query.clone();
        match coordinate {
            0 => perturbed.kinematics.creepage[0] += delta,
            1 => perturbed.kinematics.creepage[1] += delta,
            2 => perturbed.kinematics.torsional_spin_rad_per_s += delta,
            _ => {
                return Err(SmoothTangentialError::InvalidInput {
                    field: "tangent coordinate",
                });
            }
        }
        self.prepare_preview(law, patch, interface, state, &perturbed)
    }

    fn regularize_kinematics(
        &self,
        mut kinematics: PartialSlipKinematics,
    ) -> Result<PartialSlipKinematics, SmoothTangentialError> {
        // Isotropic tangent-plane deadzone: regularize the creepage
        // MAGNITUDE and rescale the vector. Componentwise smoothing is
        // not rotation covariant, so it would make the world-frame
        // receipt depend on the caller's tangent-frame choice.
        let magnitude = kinematics.creepage[0].hypot(kinematics.creepage[1]);
        let regular_magnitude = regularize(magnitude, self.regularization.creepage_scale)?;
        if magnitude > 0.0 {
            let factor = regular_magnitude / magnitude;
            kinematics.creepage[0] *= factor;
            kinematics.creepage[1] *= factor;
        }
        // The torsional spin is a scalar about the shared normal and is
        // invariant under tangent-frame rotation, so its scalar deadzone
        // is covariant as is.
        kinematics.torsional_spin_rad_per_s = regularize(
            kinematics.torsional_spin_rad_per_s,
            self.regularization.torsional_spin_scale_rad_per_s,
        )?;
        Ok(kinematics)
    }

    fn validate_query(
        &self,
        state: &SmoothTangentialState,
        query: &SmoothTangentialQuery,
        reject_duplicate: bool,
    ) -> Result<(), SmoothTangentialError> {
        nonblank(&query.query_id, "query_id")?;
        if query.expected_state_version < state.checkpoint.committed_version {
            return Err(SmoothTangentialError::StaleState {
                expected: state.checkpoint.committed_version,
                observed: query.expected_state_version,
            });
        }
        if query.expected_state_version > state.checkpoint.committed_version {
            return Err(SmoothTangentialError::FutureState {
                expected: state.checkpoint.committed_version,
                observed: query.expected_state_version,
            });
        }
        if reject_duplicate
            && state
                .checkpoint
                .committed_work_keys
                .iter()
                .any(|key| key == &query.work_ownership)
        {
            return Err(SmoothTangentialError::DuplicateWorkKey);
        }
        Ok(())
    }

    fn validate_state(
        &self,
        law: &PartialSlipLaw,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        state: &SmoothTangentialState,
    ) -> Result<(), SmoothTangentialError> {
        self.validate_inputs(patch, interface)?;
        self.restore_checkpoint(law, patch, interface, state.checkpoint.clone())
            .map(|_| ())
    }

    fn validate_inputs(
        &self,
        patch: &NormalPatchView,
        interface: &PartialSlipInterface,
    ) -> Result<(), SmoothTangentialError> {
        self.validate()?;
        if !self.authority_policy.admits(patch.authority()) {
            return Err(SmoothTangentialError::AuthorityRefused {
                input: "normal_patch",
                authority: patch.authority(),
            });
        }
        if !self.authority_policy.admits(interface.authority()) {
            return Err(SmoothTangentialError::AuthorityRefused {
                input: "ordered_interface",
                authority: interface.authority(),
            });
        }
        Ok(())
    }

    fn validate(&self) -> Result<(), SmoothTangentialError> {
        nonblank(&self.adapter_id, "adapter_id")?;
        nonblank(&self.source_id, "adapter_source_id")?;
        self.regularization.validate()
    }
}

/// Caller-owned query for one solver attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothTangentialQuery {
    /// Stable query/attempt identity; diagnostics must keep it external to this adapter.
    pub query_id: String,
    /// Version this solver evaluated.  Older/newer values refuse deterministically.
    pub expected_state_version: u64,
    /// Contact frame for the requested kinematics.
    pub frame: TangentFrame,
    /// Relative motion supplied to the constitutive law after declared smoothing.
    pub kinematics: PartialSlipKinematics,
    /// Exactly-once generalized-work ownership for a prospective commit.
    pub work_ownership: GeneralizedWorkOwnership,
}

/// Equal-and-opposite wrench receipt.  No body-velocity split is inferred here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActionReactionWrench {
    /// Wrench on the body whose relative motion was supplied to the law.
    pub action_on_declared_body: TangentialWrench,
    /// Equal-and-opposite wrench on the ordered counterbody.
    pub reaction_on_counterbody: TangentialWrench,
}

impl ActionReactionWrench {
    fn from_action(action: TangentialWrench) -> Self {
        Self {
            action_on_declared_body: action,
            reaction_on_counterbody: TangentialWrench {
                force_n: neg3(action.force_n),
                torque_nm: neg3(action.torque_nm),
            },
        }
    }
}

/// Mixed-unit residual delivered to a smooth contact solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothResidual {
    /// Tangent-frame longitudinal/lateral force residual, N.
    pub tangent_force_n: [f64; 2],
    /// Normal-axis torsional-moment residual, N m.
    pub torsional_moment_nm: f64,
    /// Endpoint relative-wrench power, W, for the regularized kinematics only.
    pub endpoint_relative_power_w: f64,
    /// Exact discrete-work power reconstructed by the constitutive rung, W.
    pub reconstructed_body_power_w: f64,
}

impl SmoothResidual {
    fn components(self) -> [f64; 3] {
        [
            self.tangent_force_n[0],
            self.tangent_force_n[1],
            self.torsional_moment_nm,
        ]
    }
}

/// Fixed-branch numerical tangent of [`SmoothResidual`] components.
#[derive(Debug, Clone, PartialEq)]
pub struct FixedBranchTangent {
    /// Branch shared by the base and all six probes.
    pub branch: PartialSlipStateKind,
    /// Rows `(F_long, F_lat, M_normal)` versus columns `(creep_long, creep_lat, spin)`.
    pub derivative: [[f64; 3]; 3],
    /// Units of the three input columns.
    pub input_units: [&'static str; 3],
    /// Units of the three output rows.
    pub output_units: [&'static str; 3],
}

/// Immutable state that advances only on a successful commit.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothTangentialState {
    checkpoint: SmoothTangentialCheckpoint,
}

impl SmoothTangentialState {
    /// Current committed version.
    #[must_use]
    pub const fn committed_version(&self) -> u64 {
        self.checkpoint.committed_version
    }

    /// Deterministic checkpoint containing adapter and upstream identities.
    #[must_use]
    pub fn checkpoint(&self) -> SmoothTangentialCheckpoint {
        self.checkpoint.clone()
    }
}

/// Full deterministic replay identity for this adapter's state.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothTangentialCheckpoint {
    adapter_id: String,
    adapter_source_id: String,
    regularization: SmoothRegularization,
    authority_policy: SmoothAuthorityPolicy,
    law_checkpoint: PartialSlipCheckpoint,
    committed_version: u64,
    committed_work_keys: Vec<GeneralizedWorkOwnership>,
    max_committed_work_keys: usize,
}

impl SmoothTangentialCheckpoint {
    /// Adapter model identity bound into the checkpoint.
    #[must_use]
    pub fn adapter_id(&self) -> &str {
        &self.adapter_id
    }

    /// Adapter configuration source bound into the checkpoint.
    #[must_use]
    pub fn adapter_source_id(&self) -> &str {
        &self.adapter_source_id
    }

    /// Full upstream partial-slip checkpoint.
    #[must_use]
    pub fn law_checkpoint(&self) -> &PartialSlipCheckpoint {
        &self.law_checkpoint
    }
}

/// Candidate receipt; call [`SmoothTangentialAdapter::commit`] or rollback.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothTangentialReceipt {
    /// Caller query identity.
    pub query_id: String,
    /// State version this candidate was based on.
    pub parent_version: u64,
    /// Constitutive branch selected by the delegated law.
    pub branch: PartialSlipStateKind,
    /// Kinematics after bounded declared smoothing.
    pub regularized_kinematics: PartialSlipKinematics,
    /// Equal-and-opposite physical wrench pair.
    pub action_reaction: ActionReactionWrench,
    /// Solver residual and power receipt.
    pub residual: SmoothResidual,
    /// Delegated constitutive receipt, including storage and heat ownership.
    pub partial_slip_step: fs_tribo::partial_slip::PartialSlipStep,
    prior_state: SmoothTangentialState,
    next_checkpoint: SmoothTangentialCheckpoint,
}

/// Refusal surface for smooth-lane transaction and derivative semantics.
#[derive(Debug, Clone, PartialEq)]
pub enum SmoothTangentialError {
    /// Required adapter/query identity was blank.
    MissingIdentity { field: &'static str },
    /// Invalid finite/budget/regularization input.
    InvalidInput { field: &'static str },
    /// A derived value was non-finite or overflowed.
    InvalidDerived { field: &'static str },
    /// The caller's policy does not admit this upstream authority label.
    AuthorityRefused {
        /// Input receipt that was refused.
        input: &'static str,
        /// Unpromoted authority label that triggered the policy.
        authority: NormalPatchAuthority,
    },
    /// Query/receipt based on an earlier state version.
    StaleState { expected: u64, observed: u64 },
    /// Query/receipt claims a state version not yet committed.
    FutureState { expected: u64, observed: u64 },
    /// A committed work key was offered again.
    DuplicateWorkKey,
    /// The explicit bounded exactly-once ledger cannot accept another key.
    WorkKeyCapacityExceeded { max: usize },
    /// Central probes crossed a constitutive branch, so no derivative is supplied.
    NoDerivativeOnBranchChange {
        /// Branch at the unperturbed point.
        base: PartialSlipStateKind,
        /// Branch at the positive probe.
        positive: PartialSlipStateKind,
        /// Branch at the negative probe.
        negative: PartialSlipStateKind,
    },
    /// A candidate was forged or based on a different same-version state.
    ReceiptDoesNotMatchState,
    /// Adapter receipt identity does not match this adapter.
    CheckpointIdentityMismatch { field: &'static str },
    /// Delegated partial-slip refusal, preserved without reinterpretation.
    PartialSlip(PartialSlipError),
}

impl fmt::Display for SmoothTangentialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingIdentity { field } => write!(f, "nonblank identity required: {field}"),
            Self::InvalidInput { field } => write!(f, "invalid smooth tangential input: {field}"),
            Self::InvalidDerived { field } => {
                write!(f, "invalid smooth tangential result: {field}")
            }
            Self::AuthorityRefused { input, authority } => {
                write!(f, "smooth lane refuses {input} authority {authority:?}")
            }
            Self::StaleState { expected, observed } => {
                write!(
                    f,
                    "stale smooth state: expected version {expected}, observed {observed}"
                )
            }
            Self::FutureState { expected, observed } => {
                write!(
                    f,
                    "future smooth state: expected version {expected}, observed {observed}"
                )
            }
            Self::DuplicateWorkKey => write!(f, "smooth lane work key was already committed"),
            Self::WorkKeyCapacityExceeded { max } => {
                write!(
                    f,
                    "smooth lane exactly-once work ledger reached its {max}-key budget"
                )
            }
            Self::NoDerivativeOnBranchChange { .. } => {
                write!(
                    f,
                    "no fixed-branch derivative across a partial-slip branch change"
                )
            }
            Self::ReceiptDoesNotMatchState => {
                write!(f, "smooth receipt does not match supplied state")
            }
            Self::CheckpointIdentityMismatch { field } => {
                write!(f, "smooth tangential checkpoint identity mismatch: {field}")
            }
            Self::PartialSlip(error) => {
                write!(f, "partial-slip law refused smooth-lane query: {error}")
            }
        }
    }
}

impl std::error::Error for SmoothTangentialError {}

impl From<PartialSlipError> for SmoothTangentialError {
    fn from(value: PartialSlipError) -> Self {
        Self::PartialSlip(value)
    }
}

/// `x - e tanh(x/e)`: odd, C-infinity, bounded departure `<= e`, and tends to `x` as `e -> 0`.
fn regularize(value: f64, scale: f64) -> Result<f64, SmoothTangentialError> {
    finite(value, "kinematics")?;
    positive(scale, "regularization scale")?;
    let output = value - scale * (value / scale).tanh();
    finite(output, "regularized kinematics")?;
    Ok(output)
}

fn neg3(value: [f64; 3]) -> [f64; 3] {
    [-value[0], -value[1], -value[2]]
}

fn nonblank(value: &str, field: &'static str) -> Result<(), SmoothTangentialError> {
    if value.trim().is_empty() {
        Err(SmoothTangentialError::MissingIdentity { field })
    } else {
        Ok(())
    }
}

fn positive(value: f64, field: &'static str) -> Result<(), SmoothTangentialError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(SmoothTangentialError::InvalidInput { field })
    }
}

fn finite(value: f64, field: &'static str) -> Result<(), SmoothTangentialError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(SmoothTangentialError::InvalidDerived { field })
    }
}
