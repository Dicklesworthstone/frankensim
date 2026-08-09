//! Transactional moving-contact response of a resolved rectangular modal base.
//!
//! This port is geometry and material agnostic. An upstream structural solve
//! supplies a mass-normalized [`RectangularPlateModalBasis`], a constitutive
//! damping state supplies per-mode loss, and each accepted mechanics interval
//! supplies an actual contact point and normal force. The same modal state can
//! therefore drive contact kinematics, structural acoustics, and rendering.

use std::sync::Arc;

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel,
};
use fs_material::visco::RayleighDamping;
use fs_math::c64::C64;

use crate::structural_acoustics::{
    PlatePointForceProjection, RectangularPlateModalBasis, StructuralModalBasisError,
};

const MODAL_BASE_LINEAGE_DOMAIN: &str =
    "org.frankensim.euler.rectangular-modal-base-step-lineage.v1";

/// Stable owner identity for one immutable structural base model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RectangularModalBaseIdentity {
    /// Physical base identity.
    pub model_id: String,
    /// Complete caller-selected configuration identity.
    pub configuration_id: String,
}

/// Immutable moving-contact base port.
#[derive(Clone, Debug)]
pub struct RectangularModalBasePort {
    identity: RectangularModalBaseIdentity,
    basis: Arc<RectangularPlateModalBasis>,
    runtime_template: ModalAcousticTimeModel,
    maximum_accepted_steps: u64,
    maximum_contact_distance_m: f64,
}

/// Accepted fixed-size lineage state for a moving-contact modal base.
#[derive(Clone, Debug, PartialEq)]
pub struct RectangularModalBaseCheckpoint {
    identity: RectangularModalBaseIdentity,
    accepted_version: u64,
    elapsed_time_s: f64,
    states: Vec<ModalAcousticState>,
    total_modal_energy_j: f64,
    cumulative_external_work_j: f64,
    cumulative_dissipation_j: f64,
    last_contact_point_base_m: [f64; 3],
    last_surface_state: ModalBaseSurfaceState,
    accepted_step_lineage_root: ContentHash,
}

impl RectangularModalBaseCheckpoint {
    /// Number of committed intervals.
    #[must_use]
    pub const fn accepted_version(&self) -> u64 {
        self.accepted_version
    }

    /// Accepted physical time [s].
    #[must_use]
    pub const fn elapsed_time_s(&self) -> f64 {
        self.elapsed_time_s
    }

    /// Mass-normalized modal coordinates at the accepted boundary.
    #[must_use]
    pub fn states(&self) -> &[ModalAcousticState] {
        &self.states
    }

    /// Total retained structural energy [J].
    #[must_use]
    pub const fn total_modal_energy_j(&self) -> f64 {
        self.total_modal_energy_j
    }

    /// Cumulative work delivered by contact [J].
    #[must_use]
    pub const fn cumulative_external_work_j(&self) -> f64 {
        self.cumulative_external_work_j
    }

    /// Cumulative viscous dissipation [J], including only admitted roundoff.
    #[must_use]
    pub const fn cumulative_dissipation_j(&self) -> f64 {
        self.cumulative_dissipation_j
    }

    /// Last accepted contact location, or the declared reference point at rest.
    #[must_use]
    pub const fn last_contact_point_base_m(&self) -> [f64; 3] {
        self.last_contact_point_base_m
    }

    /// Local motion at the last accepted contact/reference point.
    #[must_use]
    pub const fn last_surface_state(&self) -> ModalBaseSurfaceState {
        self.last_surface_state
    }

    /// Domain-separated fixed-size identity of the accepted step sequence.
    #[must_use]
    pub const fn accepted_step_lineage_root(&self) -> ContentHash {
        self.accepted_step_lineage_root
    }
}

/// Base-frame transverse surface motion at one material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalBaseSurfaceState {
    /// Upward transverse displacement [m].
    pub displacement_m: f64,
    /// Upward transverse velocity [m/s].
    pub velocity_m_per_s: f64,
}

/// One caller-owned moving-contact interval.
#[derive(Clone, Debug, PartialEq)]
pub struct RectangularModalBaseStepInput {
    /// Bounded descriptive label in the version-scoped idempotency tuple.
    pub step_id: String,
    /// Exact accepted version extended by this interval.
    pub expected_version: u64,
    /// Positive interval length [s].
    pub duration_s: f64,
    /// Base-frame contact point at interval start [m].
    pub contact_point_start_base_m: [f64; 3],
    /// Base-frame point at which the held force is projected [m].
    pub contact_point_force_base_m: [f64; 3],
    /// Base-frame contact point at interval end [m].
    pub contact_point_end_base_m: [f64; 3],
    /// Nonnegative compressive force applied into the base [N].
    pub compressive_normal_force_on_base_n: f64,
}

/// Exact accepted-interval accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct RectangularModalBaseStepReceipt {
    /// Consumed checkpoint version.
    pub parent_version: u64,
    /// Produced checkpoint version.
    pub next_version: u64,
    /// Accepted interval [s].
    pub duration_s: f64,
    /// Start contact point [m].
    pub contact_point_start_base_m: [f64; 3],
    /// Held-force projection point [m].
    pub contact_point_force_base_m: [f64; 3],
    /// End contact point [m].
    pub contact_point_end_base_m: [f64; 3],
    /// Surface state at the actual start point.
    pub surface_start: ModalBaseSurfaceState,
    /// Surface state at the actual end point.
    pub surface_end: ModalBaseSurfaceState,
    /// Compressive force acting into the base [N].
    pub compressive_normal_force_on_base_n: f64,
    /// Equal-and-opposite reaction on the disc in the base frame [N].
    pub normal_reaction_on_disc_base_n: [f64; 3],
    /// Work delivered to all retained modes [J].
    pub external_work_j: f64,
    /// Retained modal-energy change [J].
    pub stored_energy_change_j: f64,
    /// Viscous energy loss [J].
    pub viscous_dissipation_j: f64,
    /// `delta_energy - work + dissipation` [J].
    pub energy_closure_residual_j: f64,
}

/// Prepared transition whose state remains private until accepted.
#[derive(Clone, Debug, PartialEq)]
pub struct RectangularModalBaseProposal {
    parent: RectangularModalBaseCheckpoint,
    next: RectangularModalBaseCheckpoint,
    receipt: RectangularModalBaseStepReceipt,
}

impl RectangularModalBaseProposal {
    /// Immutable proposed interval accounting.
    #[must_use]
    pub const fn receipt(&self) -> &RectangularModalBaseStepReceipt {
        &self.receipt
    }
}

/// Typed admission or transactional-step refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum RectangularModalBaseError {
    /// A scalar, identity, point, or work bound is malformed.
    InvalidInput { what: &'static str },
    /// A checkpoint or proposal belongs to another immutable port.
    IdentityMismatch,
    /// The caller attempted to extend a stale or skipped version.
    VersionMismatch { expected: u64, observed: u64 },
    /// The declared trajectory-step budget is exhausted.
    StepBudgetExceeded,
    /// Structural point projection refused the contact location.
    Projection { detail: String },
    /// Exact modal time integration refused the candidate state.
    ModalTime { detail: String },
}

impl core::fmt::Display for RectangularModalBaseError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for RectangularModalBaseError {}

impl RectangularModalBasePort {
    /// Bind a resolved structural basis, damping law, and bounded time runtime.
    pub fn try_new(
        identity: RectangularModalBaseIdentity,
        basis: Arc<RectangularPlateModalBasis>,
        damping: RayleighDamping,
        nominal_sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
        maximum_accepted_steps: u64,
        maximum_contact_distance_m: f64,
    ) -> Result<Self, RectangularModalBaseError> {
        validate_identity(&identity)?;
        if maximum_accepted_steps == 0
            || !(maximum_contact_distance_m.is_finite() && maximum_contact_distance_m >= 0.0)
        {
            return Err(RectangularModalBaseError::InvalidInput {
                what: "modal base work and surface-distance bounds must be finite and nonzero where required",
            });
        }
        let modes = basis
            .modes
            .iter()
            .map(|mode| ModalAcousticMode {
                angular_frequency_rad_s: mode.angular_frequency_rad_s,
                damping_ratio: damping.zeta_at(mode.angular_frequency_rad_s),
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect();
        let runtime_template =
            ModalAcousticTimeModel::try_new(nominal_sample_rate_hz, modes, budget)
                .map_err(modal_time)?;
        Ok(Self {
            identity,
            basis,
            runtime_template,
            maximum_accepted_steps,
            maximum_contact_distance_m,
        })
    }

    /// Exact resolved structural basis shared with acoustics and visualization.
    #[must_use]
    pub fn basis(&self) -> &RectangularPlateModalBasis {
        &self.basis
    }

    /// Start an unloaded base at rest.
    #[must_use]
    pub fn initial_checkpoint(&self) -> RectangularModalBaseCheckpoint {
        self.checkpoint_from_runtime(&self.runtime_template, [0.0; 3], 0.0, 0.0, 0.0)
    }

    /// Start from static equilibrium under one already-established contact load.
    pub fn initial_static_contact_checkpoint(
        &self,
        contact_point_base_m: [f64; 3],
        compressive_normal_force_on_base_n: f64,
    ) -> Result<RectangularModalBaseCheckpoint, RectangularModalBaseError> {
        if !(compressive_normal_force_on_base_n.is_finite()
            && compressive_normal_force_on_base_n >= 0.0)
        {
            return Err(RectangularModalBaseError::InvalidInput {
                what: "initial compressive normal force must be finite and nonnegative",
            });
        }
        let projection = self.project(contact_point_base_m, -compressive_normal_force_on_base_n)?;
        let mut runtime = self.runtime_template.clone();
        runtime
            .initialize_static_equilibrium(&projection.modal_force_n_per_sqrt_kg)
            .map_err(modal_time)?;
        Ok(self.checkpoint_from_runtime(&runtime, contact_point_base_m, 0.0, 0.0, 0.0))
    }

    /// Evaluate actual local plate motion at an arbitrary admitted surface point.
    pub fn surface_state(
        &self,
        checkpoint: &RectangularModalBaseCheckpoint,
        point_base_m: [f64; 3],
    ) -> Result<ModalBaseSurfaceState, RectangularModalBaseError> {
        self.validate_checkpoint(checkpoint)?;
        let unit_projection = self.project(point_base_m, 1.0)?;
        surface_state(&checkpoint.states, &unit_projection)
    }

    /// Prepare one exact-ZOH moving-contact interval without committing it.
    pub fn propose(
        &self,
        checkpoint: &RectangularModalBaseCheckpoint,
        input: &RectangularModalBaseStepInput,
    ) -> Result<RectangularModalBaseProposal, RectangularModalBaseError> {
        self.validate_checkpoint(checkpoint)?;
        validate_step(input, checkpoint, self.maximum_accepted_steps)?;
        let surface_start = self.surface_state(checkpoint, input.contact_point_start_base_m)?;
        let force_projection = self.project(
            input.contact_point_force_base_m,
            -input.compressive_normal_force_on_base_n,
        )?;
        let mut runtime = self.runtime_template.clone();
        runtime
            .restore_states(&checkpoint.states)
            .map_err(modal_time)?;
        let frame = runtime
            .step_duration(
                &force_projection.modal_force_n_per_sqrt_kg,
                input.duration_s,
            )
            .map_err(modal_time)?;
        let states = runtime.states().to_vec();
        let unit_end = self.project(input.contact_point_end_base_m, 1.0)?;
        let surface_end = surface_state(&states, &unit_end)?;
        let stored_energy_change_j = frame.total_modal_energy_j - checkpoint.total_modal_energy_j;
        let energy_closure_residual_j =
            stored_energy_change_j - frame.input_work_j + frame.viscous_dissipation_j;
        let next = RectangularModalBaseCheckpoint {
            identity: checkpoint.identity.clone(),
            accepted_version: checkpoint.accepted_version + 1,
            elapsed_time_s: checkpoint.elapsed_time_s + input.duration_s,
            states,
            total_modal_energy_j: frame.total_modal_energy_j,
            cumulative_external_work_j: checkpoint.cumulative_external_work_j + frame.input_work_j,
            cumulative_dissipation_j: checkpoint.cumulative_dissipation_j
                + frame.viscous_dissipation_j,
            last_contact_point_base_m: input.contact_point_end_base_m,
            last_surface_state: surface_end,
            accepted_step_lineage_root: extend_lineage(
                checkpoint.accepted_step_lineage_root,
                checkpoint.accepted_version,
                &input.step_id,
            ),
        };
        Ok(RectangularModalBaseProposal {
            parent: checkpoint.clone(),
            next,
            receipt: RectangularModalBaseStepReceipt {
                parent_version: checkpoint.accepted_version,
                next_version: checkpoint.accepted_version + 1,
                duration_s: input.duration_s,
                contact_point_start_base_m: input.contact_point_start_base_m,
                contact_point_force_base_m: input.contact_point_force_base_m,
                contact_point_end_base_m: input.contact_point_end_base_m,
                surface_start,
                surface_end,
                compressive_normal_force_on_base_n: input.compressive_normal_force_on_base_n,
                normal_reaction_on_disc_base_n: [
                    0.0,
                    0.0,
                    input.compressive_normal_force_on_base_n,
                ],
                external_work_j: frame.input_work_j,
                stored_energy_change_j,
                viscous_dissipation_j: frame.viscous_dissipation_j,
                energy_closure_residual_j,
            },
        })
    }

    /// Commit a proposal only when it exactly extends the supplied checkpoint.
    pub fn accept(
        &self,
        checkpoint: &RectangularModalBaseCheckpoint,
        proposal: RectangularModalBaseProposal,
    ) -> Result<RectangularModalBaseCheckpoint, RectangularModalBaseError> {
        self.validate_checkpoint(checkpoint)?;
        if proposal.parent != *checkpoint || proposal.next.identity != self.identity {
            return Err(RectangularModalBaseError::IdentityMismatch);
        }
        Ok(proposal.next)
    }

    fn project(
        &self,
        point_base_m: [f64; 3],
        transverse_force_n: f64,
    ) -> Result<PlatePointForceProjection, RectangularModalBaseError> {
        self.basis
            .project_transverse_point_force(
                point_base_m,
                transverse_force_n,
                self.maximum_contact_distance_m,
            )
            .map_err(projection)
    }

    fn checkpoint_from_runtime(
        &self,
        runtime: &ModalAcousticTimeModel,
        last_contact_point_base_m: [f64; 3],
        elapsed_time_s: f64,
        cumulative_external_work_j: f64,
        cumulative_dissipation_j: f64,
    ) -> RectangularModalBaseCheckpoint {
        let last_surface_state = self
            .project(last_contact_point_base_m, 1.0)
            .and_then(|projection| surface_state(runtime.states(), &projection))
            .expect("an admitted modal basis contains its declared reference contact point");
        RectangularModalBaseCheckpoint {
            identity: self.identity.clone(),
            accepted_version: 0,
            elapsed_time_s,
            states: runtime.states().to_vec(),
            total_modal_energy_j: total_energy(runtime.modes(), runtime.states()),
            cumulative_external_work_j,
            cumulative_dissipation_j,
            last_contact_point_base_m,
            last_surface_state,
            accepted_step_lineage_root: hash_domain(MODAL_BASE_LINEAGE_DOMAIN, &[]),
        }
    }

    fn validate_checkpoint(
        &self,
        checkpoint: &RectangularModalBaseCheckpoint,
    ) -> Result<(), RectangularModalBaseError> {
        if checkpoint.identity != self.identity
            || checkpoint.states.len() != self.runtime_template.modes().len()
            || checkpoint.accepted_version > self.maximum_accepted_steps
        {
            return Err(RectangularModalBaseError::IdentityMismatch);
        }
        Ok(())
    }
}

fn validate_identity(
    identity: &RectangularModalBaseIdentity,
) -> Result<(), RectangularModalBaseError> {
    if [&identity.model_id, &identity.configuration_id]
        .into_iter()
        .any(|value| value.is_empty() || value.len() > 256 || !value.is_ascii())
    {
        return Err(RectangularModalBaseError::InvalidInput {
            what: "modal base identity must be bounded nonempty ASCII",
        });
    }
    Ok(())
}

fn validate_step(
    input: &RectangularModalBaseStepInput,
    checkpoint: &RectangularModalBaseCheckpoint,
    maximum_accepted_steps: u64,
) -> Result<(), RectangularModalBaseError> {
    if input.step_id.is_empty() || input.step_id.len() > 256 || !input.step_id.is_ascii() {
        return Err(RectangularModalBaseError::InvalidInput {
            what: "modal base step identity must be bounded nonempty ASCII",
        });
    }
    if input.expected_version != checkpoint.accepted_version {
        return Err(RectangularModalBaseError::VersionMismatch {
            expected: checkpoint.accepted_version,
            observed: input.expected_version,
        });
    }
    if checkpoint.accepted_version >= maximum_accepted_steps {
        return Err(RectangularModalBaseError::StepBudgetExceeded);
    }
    if !(input.duration_s.is_finite() && input.duration_s > 0.0)
        || !(input.compressive_normal_force_on_base_n.is_finite()
            && input.compressive_normal_force_on_base_n >= 0.0)
        || input
            .contact_point_start_base_m
            .iter()
            .chain(&input.contact_point_force_base_m)
            .chain(&input.contact_point_end_base_m)
            .any(|value| !value.is_finite())
    {
        return Err(RectangularModalBaseError::InvalidInput {
            what: "modal base duration, force, and moving contact points must be finite and physical",
        });
    }
    Ok(())
}

fn surface_state(
    states: &[ModalAcousticState],
    unit_force_projection: &PlatePointForceProjection,
) -> Result<ModalBaseSurfaceState, RectangularModalBaseError> {
    if states.len() != unit_force_projection.modal_force_n_per_sqrt_kg.len() {
        return Err(RectangularModalBaseError::InvalidInput {
            what: "modal state and structural projection cardinalities differ",
        });
    }
    let mut displacement_m = 0.0;
    let mut velocity_m_per_s = 0.0;
    for (state, shape_per_sqrt_kg) in states
        .iter()
        .zip(&unit_force_projection.modal_force_n_per_sqrt_kg)
    {
        displacement_m += state.displacement_m_sqrt_kg * shape_per_sqrt_kg;
        velocity_m_per_s += state.velocity_m_sqrt_kg_per_s * shape_per_sqrt_kg;
    }
    if !(displacement_m.is_finite() && velocity_m_per_s.is_finite()) {
        return Err(RectangularModalBaseError::InvalidInput {
            what: "modal surface evaluation became nonfinite",
        });
    }
    Ok(ModalBaseSurfaceState {
        displacement_m,
        velocity_m_per_s,
    })
}

fn total_energy(modes: &[ModalAcousticMode], states: &[ModalAcousticState]) -> f64 {
    modes
        .iter()
        .zip(states)
        .map(|(mode, state)| {
            0.5 * state.velocity_m_sqrt_kg_per_s.powi(2)
                + 0.5 * mode.angular_frequency_rad_s.powi(2) * state.displacement_m_sqrt_kg.powi(2)
        })
        .sum()
}

fn extend_lineage(parent: ContentHash, parent_version: u64, step_id: &str) -> ContentHash {
    let mut hasher = DomainHasher::new(MODAL_BASE_LINEAGE_DOMAIN);
    hasher.update(parent.as_bytes());
    hasher.update(&parent_version.to_le_bytes());
    let length = u64::try_from(step_id.len()).expect("step identities are at most 256 bytes");
    hasher.update(&length.to_le_bytes());
    hasher.update(step_id.as_bytes());
    hasher.finalize()
}

fn projection(error: StructuralModalBasisError) -> RectangularModalBaseError {
    RectangularModalBaseError::Projection {
        detail: error.to_string(),
    }
}

fn modal_time(error: ModalAcousticTimeError) -> RectangularModalBaseError {
    RectangularModalBaseError::ModalTime {
        detail: error.to_string(),
    }
}
