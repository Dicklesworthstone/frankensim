//! Euler-disc parameterization of the generic reduced exterior-air wrench.
//!
//! This is an adapter, not an Euler-disc air law.  It admits only free exterior
//! gas around the disc and delegates force, torque, correlation-domain, and
//! exact-once work rules to `fs_flux::reduced_aero`.  Thin-gap pressure and
//! target-fitted terms have no representation here.

use core::fmt;

use fs_flux::{
    AlternativeWrenchSet, BodyKinematics, CandidateWrench, DiscGeometry, DiscPose, GasProperties,
    GasPropertyCard, ReducedAeroError, ReducedAeroInput, ReducedAeroModel, SurfaceRoughness, Vec3,
    WorkReceipt, WorkWindow,
};
use fs_tribo::{ExactlyOnceKeyError, ExactlyOnceKeyLedger};

const UNIT_TOLERANCE: f64 = 1.0e-12;

/// Explicit spatial ownership for this adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExternalAirDomain {
    /// Free gas surrounding the disc's exterior faces and rim.
    ExteriorFreeGas,
    /// Thin-gap pressure belongs to `fs_flux::gas_film`, never this adapter.
    ThinGap,
}

/// A named, right-handed body frame expressed in the request's world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EulerDiscBodyFrame {
    /// Unit body x direction in the world frame.
    pub x_world: Vec3,
    /// Unit body z direction/away-from-base disc normal in the world frame.
    pub z_world: Vec3,
}

impl EulerDiscBodyFrame {
    fn validate(self) -> Result<(), ExternalAirError> {
        for vector in [self.x_world, self.z_world] {
            if !(vector.x.is_finite() && vector.y.is_finite() && vector.z.is_finite()) {
                return Err(ExternalAirError::InvalidInput {
                    field: "body_frame",
                });
            }
        }
        let x_norm = norm(self.x_world);
        let z_norm = norm(self.z_world);
        if (x_norm - 1.0).abs() > UNIT_TOLERANCE || (z_norm - 1.0).abs() > UNIT_TOLERANCE {
            return Err(ExternalAirError::NonUnitBodyFrame { x_norm, z_norm });
        }
        if dot(self.x_world, self.z_world).abs() > UNIT_TOLERANCE {
            return Err(ExternalAirError::NonOrthogonalBodyFrame);
        }
        Ok(())
    }

    fn y_world(self) -> Vec3 {
        cross(self.z_world, self.x_world)
    }

    fn into_body(self, vector_world: Vec3) -> Vec3 {
        Vec3::new(
            dot(vector_world, self.x_world),
            dot(vector_world, self.y_world()),
            dot(vector_world, self.z_world),
        )
    }
}

/// Euler-disc exterior geometry admitted by the generic circular-disc API.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EulerDiscExteriorGeometry {
    /// Disc radius [m].
    pub radius_m: f64,
    /// Exterior rim thickness [m], not an air-film gap.
    pub exterior_thickness_m: f64,
}

/// Pose and rigid-body rates for one exterior-air request.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EulerDiscExteriorState {
    /// Disc centre/reference point [m] in the named world frame.
    pub center_world_m: Vec3,
    /// Linear velocity [m s^-1] in the named world frame.
    pub center_velocity_world_m_per_s: Vec3,
    /// Angular velocity [rad s^-1] in the named world frame.
    pub angular_velocity_world_rad_per_s: Vec3,
    /// Orientation used to report body-frame wrench components.
    pub body_frame: EulerDiscBodyFrame,
}

/// Identity and declared source routing retained by the adapter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalAirIdentity {
    /// Caller-owned case identity.
    pub case_id: String,
    /// Explicit inertial/world-frame identity.
    pub world_frame_id: String,
    /// Named disc body-frame identity used for body-frame wrench components.
    pub body_frame_id: String,
    /// Source identity for the admitted disc geometry.
    pub geometry_source_id: String,
    /// Source identity for pose and rate inputs.
    pub state_source_id: String,
    /// Source identity for the external-domain declaration.
    pub domain_source_id: String,
}

/// Full Euler-disc request to the generic exterior-air API.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerExternalAirInput {
    /// Domain must be [`ExternalAirDomain::ExteriorFreeGas`].
    pub domain: ExternalAirDomain,
    /// Case, frame, geometry/state, and domain identities.
    pub identity: ExternalAirIdentity,
    /// Exterior geometry only.
    pub geometry: EulerDiscExteriorGeometry,
    /// Pose and velocities.
    pub state: EulerDiscExteriorState,
    /// Complete gas property card, including far-field world velocity.
    pub gas: GasPropertyCard,
    /// Absolute pressure and source identity retained without inventing a scaling law.
    pub pressure: ExteriorAirPressure,
    /// Exterior surface roughness card. It is not a gas-film roughness rule.
    pub exterior_roughness: SurfaceRoughness,
    /// Independently retained generic correlations. They are never averaged.
    pub alternatives: Vec<ReducedAeroModel>,
}

/// Pressure state retained beside the gas-property card.
#[derive(Debug, Clone, PartialEq)]
pub struct ExteriorAirPressure {
    /// Absolute far-field pressure [Pa].
    pub absolute_pressure_pa: f64,
    /// Source identity for this pressure state.
    pub source_id: String,
}

/// A generic exterior wrench with its equivalent components in the disc body frame.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerExternalAirCandidate {
    /// Generic result, force and torque in the named world frame.
    pub world_wrench: CandidateWrench,
    /// Same total force in the named disc body frame [N].
    pub force_body_n: Vec3,
    /// Same total torque in the named disc body frame [N m].
    pub torque_body_n_m: Vec3,
    /// Dissipation has no allocated thermal destination in this mechanical model.
    pub heat: ExteriorAirHeatDisposition,
}

/// Honest thermal boundary for an isothermal mechanical wrench correlation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExteriorAirHeatDisposition {
    /// Relative mechanical dissipation is reported, but no gas/body heat split is modelled.
    UnallocatedNoThermalModel,
}

/// Applied exterior-air output with retained alternative disagreement.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerExternalAirSet {
    /// Case/frame/source identity carried from the admitted request.
    pub identity: ExternalAirIdentity,
    /// Retained pressure state. It does not alter a correlation that lacks pressure dependence.
    pub pressure: ExteriorAirPressure,
    /// Identity of the free exterior-gas domain.
    pub domain: ExternalAirDomain,
    /// The generic correlation applicability envelopes admitted this exact state.
    pub applicability: ExteriorAirApplicability,
    /// The pressure-scaling boundary of the mapped generic API.
    pub pressure_scaling: ExteriorAirPressureScaling,
    /// All candidates in the generic API's deterministic correlation-identity order.
    pub candidates: Vec<EulerExternalAirCandidate>,
    /// True when distinct alternatives produce a different force or torque.
    pub has_force_or_torque_disagreement: bool,
}

/// Pressure treatment that prevents an unimplemented universal scaling claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExteriorAirPressureScaling {
    /// The generic correlation receives density/viscosity/sound speed, not pressure directly.
    /// A caller must supply a gas card evaluated at the pressure state it declares.
    NoDirectScaling,
}

/// Applicability that is established by successful generic model admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExteriorAirApplicability {
    /// This adapter admits only free exterior gas, never a thin gap.
    pub domain: ExternalAirDomain,
    /// Every returned candidate passed its own generic correlation envelope.
    pub generic_correlation_domain_admitted: bool,
    /// Outputs retain only the generic model's Estimate-only authority.
    pub estimate_only: bool,
}

/// Typed refusal from the Euler-disc exterior-air adapter.
#[derive(Debug, Clone, PartialEq)]
pub enum ExternalAirError {
    /// A required source/frame/case identity is empty or non-canonical.
    InvalidIdentity { field: &'static str },
    /// A scalar/vector does not meet the adapter's input contract.
    InvalidInput { field: &'static str },
    /// The declared body-frame axes are not unit directions.
    NonUnitBodyFrame { x_norm: f64, z_norm: f64 },
    /// The declared body-frame axes do not form an orthogonal disc frame.
    NonOrthogonalBodyFrame,
    /// Thin-gap pressure is deliberately routed elsewhere.
    ThinGapDomainRejected,
    /// A generic correlation model refused the mapped exterior request.
    GenericRefusal { detail: ReducedAeroError },
    /// A caller tried to stage exterior work already accepted in this state.
    DuplicateAcceptedWork { exchange_key: u64 },
    /// A fixed-memory transaction supplied a skipped or reordered exchange key.
    OutOfSequenceAcceptedWork { exchange_key: u64 },
    /// The caller-declared bounded exactly-once ledger has no capacity for a
    /// new accepted exterior-work interval.
    AcceptedWorkCapacityExceeded { maximum: usize },
    /// A staged exterior-work proposal did not originate from the state being
    /// accepted, so it cannot be applied safely.
    ProposalDoesNotMatchState,
}

impl fmt::Display for ExternalAirError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ExternalAirError {}

/// Map one admitted Euler-disc exterior-air request to the generic model.
///
/// Each successful candidate remains Estimate-only, including its correlation,
/// gas, roughness, uncertainty, applicability, mechanical-power, and
/// moving-ambient accounting receipt. A body-only power sign is not a total
/// energy claim when the ambient moves.
pub fn evaluate_euler_disc_external_air(
    input: &EulerExternalAirInput,
) -> Result<EulerExternalAirSet, ExternalAirError> {
    if input.domain != ExternalAirDomain::ExteriorFreeGas {
        return Err(ExternalAirError::ThinGapDomainRejected);
    }
    validate_identity(&input.identity.case_id, "identity.case_id")?;
    validate_identity(&input.identity.world_frame_id, "identity.world_frame_id")?;
    validate_identity(&input.identity.body_frame_id, "identity.body_frame_id")?;
    validate_identity(
        &input.identity.geometry_source_id,
        "identity.geometry_source_id",
    )?;
    validate_identity(&input.identity.state_source_id, "identity.state_source_id")?;
    validate_identity(
        &input.identity.domain_source_id,
        "identity.domain_source_id",
    )?;
    validate_identity(&input.pressure.source_id, "pressure.source_id")?;
    if !(input.pressure.absolute_pressure_pa.is_finite()
        && input.pressure.absolute_pressure_pa > 0.0)
    {
        return Err(ExternalAirError::InvalidInput {
            field: "pressure.absolute_pressure_pa",
        });
    }
    input.state.body_frame.validate()?;

    let geometry = DiscGeometry {
        radius_m: input.geometry.radius_m,
        exterior_thickness_m: input.geometry.exterior_thickness_m,
    };
    let pose = DiscPose::try_new(input.state.body_frame.z_world).map_err(generic_refusal)?;
    let gas = GasProperties::try_from(input.gas.clone()).map_err(generic_refusal)?;
    let generic_input = ReducedAeroInput {
        world_frame_id: input.identity.world_frame_id.clone(),
        geometry,
        pose,
        kinematics: BodyKinematics {
            reference_point_world_m: input.state.center_world_m,
            linear_velocity_world_m_per_s: input.state.center_velocity_world_m_per_s,
            angular_velocity_world_rad_per_s: input.state.angular_velocity_world_rad_per_s,
        },
        gas,
        roughness: input.exterior_roughness.clone(),
    };
    let generic_set = AlternativeWrenchSet::evaluate(&input.alternatives, &generic_input)
        .map_err(generic_refusal)?;
    Ok(map_set(
        generic_set,
        input.identity.clone(),
        input.pressure.clone(),
        input.state.body_frame,
    ))
}

/// Exact-once exterior-air work accounting. Thermal allocation is still unavailable.
#[derive(Debug, Default)]
pub struct EulerExternalAirWorkWindow {
    inner: WorkWindow,
}

impl EulerExternalAirWorkWindow {
    /// Record one candidate's generic body work and passive relative dissipation once.
    pub fn record_once(
        &mut self,
        exchange_key: u64,
        duration_s: f64,
        candidate: &EulerExternalAirCandidate,
    ) -> Result<WorkReceipt, ExternalAirError> {
        self.inner
            .record_once(exchange_key, duration_s, &candidate.world_wrench)
            .map_err(generic_refusal)
    }

    /// Accumulated work into the body [J]. It may have either sign with moving ambient gas.
    #[must_use]
    pub fn body_work_j(&self) -> f64 {
        self.inner.body_work_j()
    }

    /// Accumulated passive relative dissipation [J], not a heat-allocation claim.
    #[must_use]
    pub fn relative_dissipation_j(&self) -> f64 {
        self.inner.relative_dissipation_j()
    }
}

/// Cloneable accepted-step accounting state for the exterior free-gas domain.
///
/// This is deliberately separate from [`EulerExternalAirWorkWindow`], whose
/// immediate mutation remains useful for one-shot callers. A coupled solver
/// should use this state to stage a selected correlation candidate, accept it
/// only with the outer mechanics step, and retain the result in its checkpoint.
/// It never selects or averages alternatives and it cannot account for a gas
/// film or any non-exterior channel.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerExternalAirWorkState {
    state_id: String,
    exchange_ledger: ExactlyOnceKeyLedger<u64>,
    body_work_j: f64,
    relative_dissipation_j: f64,
    committed_version: u64,
}

impl EulerExternalAirWorkState {
    /// Creates an identity-bound, bounded exterior-work checkpoint.
    ///
    /// The capacity is an admission budget, not a ring-buffer: accepted keys
    /// are never discarded, because doing so would re-admit duplicate work.
    pub fn new(
        state_id: impl Into<String>,
        maximum_committed_exchange_keys: usize,
    ) -> Result<Self, ExternalAirError> {
        let state_id = state_id.into();
        validate_identity(&state_id, "exterior_work_state_id")?;
        if maximum_committed_exchange_keys == 0 {
            return Err(ExternalAirError::InvalidInput {
                field: "maximum_committed_exchange_keys",
            });
        }
        Ok(Self {
            state_id,
            exchange_ledger: ExactlyOnceKeyLedger::retained_set(maximum_committed_exchange_keys)
                .map_err(exactly_once_error)?,
            body_work_j: 0.0,
            relative_dissipation_j: 0.0,
            committed_version: 0,
        })
    }

    /// Creates a fixed-memory ledger that accepts `1, 2, ...` in order.
    pub fn new_strict_sequence(
        state_id: impl Into<String>,
        maximum_committed_exchange_keys: usize,
    ) -> Result<Self, ExternalAirError> {
        let state_id = state_id.into();
        validate_identity(&state_id, "exterior_work_state_id")?;
        Ok(Self {
            state_id,
            exchange_ledger: ExactlyOnceKeyLedger::strict_sequence(
                1,
                maximum_committed_exchange_keys,
            )
            .map_err(exactly_once_error)?,
            body_work_j: 0.0,
            relative_dissipation_j: 0.0,
            committed_version: 0,
        })
    }

    /// Stable caller-owned identity retained across accept/refuse/checkpoint
    /// transitions for this one exterior-gas work ledger.
    #[must_use]
    pub fn state_id(&self) -> &str {
        &self.state_id
    }

    /// Maximum accepted exchange keys retained by this checkpoint.
    #[must_use]
    pub const fn maximum_committed_exchange_keys(&self) -> usize {
        self.exchange_ledger.maximum_committed()
    }

    /// Total accepted work into the body [J]. Moving ambient gas may make this
    /// positive; it is not a total system-energy closure claim.
    #[must_use]
    pub const fn body_work_j(&self) -> f64 {
        self.body_work_j
    }

    /// Total accepted passive relative dissipation [J].
    #[must_use]
    pub const fn relative_dissipation_j(&self) -> f64 {
        self.relative_dissipation_j
    }

    /// Count of accepted exterior-work proposals.
    #[must_use]
    pub const fn committed_version(&self) -> u64 {
        self.committed_version
    }

    /// Stages the exact work receipt for one already-selected exterior
    /// candidate. The generic `WorkWindow` is used as a non-mutating validator
    /// for the candidate arithmetic; the returned proposal owns no shared
    /// mutation until [`Self::commit`] is called.
    pub fn prepare(
        &self,
        exchange_key: u64,
        duration_s: f64,
        candidate: &EulerExternalAirCandidate,
    ) -> Result<EulerExternalAirWorkProposal, ExternalAirError> {
        let strict_successor = self
            .exchange_ledger
            .strict_next_key()
            .map(|_| {
                exchange_key
                    .checked_add(1)
                    .ok_or(ExternalAirError::InvalidInput {
                        field: "exterior work exchange-key successor",
                    })
            })
            .transpose()?;
        let exchange_ledger = self
            .exchange_ledger
            .advance(&exchange_key, strict_successor)
            .map_err(|error| exactly_once_prepare_error(error, exchange_key))?;
        let mut validator = WorkWindow::default();
        let receipt = validator
            .record_once(exchange_key, duration_s, &candidate.world_wrench)
            .map_err(generic_refusal)?;
        let body_work_j = self.body_work_j + receipt.body_work_j;
        let relative_dissipation_j = self.relative_dissipation_j + receipt.relative_dissipation_j;
        if !(body_work_j.is_finite() && relative_dissipation_j.is_finite()) {
            return Err(ExternalAirError::InvalidInput {
                field: "accepted exterior work accumulation",
            });
        }
        let committed_version =
            self.committed_version
                .checked_add(1)
                .ok_or(ExternalAirError::InvalidInput {
                    field: "accepted exterior work version",
                })?;
        Ok(EulerExternalAirWorkProposal {
            receipt,
            correlation_id: candidate.world_wrench.correlation.id.clone(),
            parent_version: self.committed_version,
            prior_state: self.clone(),
            next_state: Self {
                state_id: self.state_id.clone(),
                exchange_ledger,
                body_work_j,
                relative_dissipation_j,
                committed_version,
            },
        })
    }

    /// Commits a staged exterior-work proposal only against its exact parent.
    pub fn commit(
        &self,
        proposal: &EulerExternalAirWorkProposal,
    ) -> Result<Self, ExternalAirError> {
        if *self != proposal.prior_state || self.committed_version != proposal.parent_version {
            return Err(ExternalAirError::ProposalDoesNotMatchState);
        }
        Ok(proposal.next_state.clone())
    }

    /// Refuses a staged exterior-work proposal and recovers its exact parent.
    #[must_use]
    pub fn rollback(proposal: &EulerExternalAirWorkProposal) -> Self {
        proposal.prior_state.clone()
    }
}

/// Uncommitted exact-once exterior-work candidate. The composition owner
/// chooses the correlation before this API is called; no ranking is encoded.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerExternalAirWorkProposal {
    /// Generic work receipt for the selected candidate and duration.
    pub receipt: WorkReceipt,
    /// Retained selected-correlation identity; alternatives remain available
    /// on the original [`EulerExternalAirSet`].
    pub correlation_id: String,
    /// Accepted-state version observed while staging this proposal.
    pub parent_version: u64,
    prior_state: EulerExternalAirWorkState,
    next_state: EulerExternalAirWorkState,
}

fn map_set(
    generic_set: AlternativeWrenchSet,
    identity: ExternalAirIdentity,
    pressure: ExteriorAirPressure,
    frame: EulerDiscBodyFrame,
) -> EulerExternalAirSet {
    let has_force_or_torque_disagreement = generic_set.has_force_or_torque_disagreement();
    let candidates = generic_set
        .candidates
        .into_iter()
        .map(|world_wrench| EulerExternalAirCandidate {
            force_body_n: frame.into_body(world_wrench.force_world_n),
            torque_body_n_m: frame.into_body(world_wrench.torque_world_n_m),
            world_wrench,
            heat: ExteriorAirHeatDisposition::UnallocatedNoThermalModel,
        })
        .collect();
    EulerExternalAirSet {
        identity,
        pressure,
        domain: ExternalAirDomain::ExteriorFreeGas,
        applicability: ExteriorAirApplicability {
            domain: ExternalAirDomain::ExteriorFreeGas,
            generic_correlation_domain_admitted: true,
            estimate_only: true,
        },
        pressure_scaling: ExteriorAirPressureScaling::NoDirectScaling,
        candidates,
        has_force_or_torque_disagreement,
    }
}

fn generic_refusal(detail: ReducedAeroError) -> ExternalAirError {
    ExternalAirError::GenericRefusal { detail }
}

fn validate_identity(value: &str, field: &'static str) -> Result<(), ExternalAirError> {
    if value.is_empty()
        || !value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(ExternalAirError::InvalidIdentity { field });
    }
    Ok(())
}

fn exactly_once_error(error: ExactlyOnceKeyError) -> ExternalAirError {
    exactly_once_prepare_error(error, 0)
}

fn exactly_once_prepare_error(error: ExactlyOnceKeyError, exchange_key: u64) -> ExternalAirError {
    match error {
        ExactlyOnceKeyError::Duplicate => ExternalAirError::DuplicateAcceptedWork { exchange_key },
        ExactlyOnceKeyError::OutOfSequence => {
            ExternalAirError::OutOfSequenceAcceptedWork { exchange_key }
        }
        ExactlyOnceKeyError::CapacityExceeded { maximum } => {
            ExternalAirError::AcceptedWorkCapacityExceeded { maximum }
        }
        ExactlyOnceKeyError::ZeroCapacity | ExactlyOnceKeyError::MissingSuccessor => {
            ExternalAirError::InvalidInput {
                field: "exterior work exactly-once ledger",
            }
        }
    }
}

fn dot(left: Vec3, right: Vec3) -> f64 {
    left.x
        .mul_add(right.x, left.y.mul_add(right.y, left.z * right.z))
}

fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}

fn cross(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(
        left.y.mul_add(right.z, -(left.z * right.y)),
        left.z.mul_add(right.x, -(left.x * right.z)),
        left.x.mul_add(right.y, -(left.y * right.x)),
    )
}
