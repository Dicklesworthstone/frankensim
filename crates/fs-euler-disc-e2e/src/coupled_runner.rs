//! Closed, deterministic reduced Euler-disc trajectory runner.
//!
//! This runner is deliberately a bounded numerical model, not a validation
//! claim.  It advances the production `fs_mbd` rigid-body state while retaining
//! separately evaluated contact, rolling, base, and gas wrenches.  The new
//! Euler adapter modules are exported by the crate, but this runner does not
//! yet consume their receipts: its channel laws remain deliberately reduced.

use core::fmt;

use crate::specimen::ResolvedDiscProfile;
use crate::{
    ContactDiscGeometry, ContactGeometry, contact_geometry, profile_contact_geometry,
    profile_state_at_ground_contact, state_at_ground_contact,
};
use fs_exec::Cx;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

/// A factor set supplied by a campaign case, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledFactors {
    pub mass_kg: f64,
    pub radius_m: f64,
    pub thickness_m: f64,
    pub transverse_inertia_kg_m2: f64,
    pub axial_inertia_kg_m2: f64,
    pub gravity_m_per_s2: f64,
    /// Saturated gross-sliding Coulomb coefficient; no static-stick solve is performed.
    pub sliding_friction_coefficient: f64,
    pub rolling_resistance_m: f64,
    /// Kelvin-Voigt normal-contact stiffness [N/m].
    pub contact_stiffness_n_per_m: f64,
    /// Kelvin-Voigt normal-contact damping [N s/m].
    pub contact_damping_n_s_per_m: f64,
    /// Effective vertical inertia of the base mode [kg].
    pub base_effective_mass_kg: f64,
    pub base_stiffness_n_per_m: f64,
    pub base_damping_n_s_per_m: f64,
    pub gas_rotational_damping_n_m_s: f64,
    pub gas_translation_damping_n_s_per_m: f64,
}

/// Geometry-independent material and environment factors, in SI units.
///
/// Profile trajectories derive mass, center of mass, principal inertia, outer
/// radius, and thickness from one [`ResolvedDiscProfile`]. Keeping those
/// quantities out of this input prevents a caller from pairing one contact
/// shape with hand-entered inertia from another shape.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledChannelFactors {
    pub gravity_m_per_s2: f64,
    pub sliding_friction_coefficient: f64,
    pub rolling_resistance_m: f64,
    pub contact_stiffness_n_per_m: f64,
    pub contact_damping_n_s_per_m: f64,
    pub base_effective_mass_kg: f64,
    pub base_stiffness_n_per_m: f64,
    pub base_damping_n_s_per_m: f64,
    pub gas_rotational_damping_n_m_s: f64,
    pub gas_translation_damping_n_s_per_m: f64,
}

impl CoupledFactors {
    pub fn channel_factors(self) -> CoupledChannelFactors {
        CoupledChannelFactors {
            gravity_m_per_s2: self.gravity_m_per_s2,
            sliding_friction_coefficient: self.sliding_friction_coefficient,
            rolling_resistance_m: self.rolling_resistance_m,
            contact_stiffness_n_per_m: self.contact_stiffness_n_per_m,
            contact_damping_n_s_per_m: self.contact_damping_n_s_per_m,
            base_effective_mass_kg: self.base_effective_mass_kg,
            base_stiffness_n_per_m: self.base_stiffness_n_per_m,
            base_damping_n_s_per_m: self.base_damping_n_s_per_m,
            gas_rotational_damping_n_m_s: self.gas_rotational_damping_n_m_s,
            gas_translation_damping_n_s_per_m: self.gas_translation_damping_n_s_per_m,
        }
    }
}

/// Fixed deterministic run controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledControls {
    pub timestep_s: f64,
    pub maximum_steps: u32,
    pub terminal_inclination_rad: f64,
    pub reimpact_limit: u32,
}

/// Initial pose/twist.  Inclination is the angle of the disc plane above the base.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledInitialState {
    pub inclination_rad: f64,
    pub precession_rad_per_s: f64,
    pub spin_rad_per_s: f64,
}

/// An independently evaluated channel wrench and its body work over one step.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChannelWrench {
    pub force_world_n: Vec3,
    pub torque_world_nm: Vec3,
    pub work_j: f64,
}

/// Channel ownership is explicit: no channel may borrow another's work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChannelOwnership {
    pub gravity: ChannelWrench,
    pub contact: ChannelWrench,
    pub rolling: ChannelWrench,
    pub base: ChannelWrench,
    pub gas: ChannelWrench,
}

/// One retained, time-evolving sample.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledSample {
    pub time_s: f64,
    /// Complete accepted rigid-body state at `time_s`.
    pub state: RigidBodyState,
    /// World-frame center-of-mass velocity derived from `state` and the exact
    /// mass properties used by the runner [m/s].
    pub center_of_mass_velocity_world_m_per_s: Vec3,
    /// Complete accepted one-mode base displacement [m].
    pub base_deflection_m: f64,
    /// Complete accepted one-mode base velocity [m/s].
    pub base_velocity_m_per_s: f64,
    /// Unilateral branch immediately after the accepted interval.
    pub contact_branch: CoupledContactBranch,
    /// Contact/support geometry evaluated at the accepted endpoint pose. The
    /// geometry is retained even on the open branch for diagnostics; public
    /// render trajectories expose it only for closed contact.
    pub endpoint_contact_geometry: ContactGeometry,
    /// Signed endpoint gap relative to the moving base mode [m].
    pub endpoint_signed_gap_m: f64,
    pub inclination_rad: f64,
    pub precession_rad_per_s: f64,
    pub spin_rad_per_s: f64,
    /// Finite-difference `d(precession)/dt` [rad/s²].
    pub precession_acceleration_rad_per_s2: f64,
    /// Gap evaluated at the beginning of the interval ending at `time_s` [m].
    pub interval_start_gap_m: f64,
    /// Unilateral normal force evaluated over the interval ending at `time_s` [N].
    pub interval_normal_force_n: f64,
    /// True when any accepted subinterval used the closed signed-gap branch.
    /// This is intentionally not inferred from the normal-force magnitude:
    /// a newly localized reimpact may begin at zero penetration.
    pub contact_active: bool,
    /// Every bounded, mechanically split contact transition retained within
    /// this interval, in chronological order. The locator probes four fixed
    /// subintervals before bisection. An empty list therefore does not prove
    /// that no open/reimpact excursion occurred entirely between adjacent scan
    /// nodes, and makes no continuum-scale chatter claim.
    pub contact_transitions: Vec<LocalizedContactTransition>,
    /// Bracket for the terminal-inclination crossing that ended this sample,
    /// when one occurred. The locator re-evolves the reduced mechanics within
    /// its bracket; it is a numerical event diagnostic, not a measurement or
    /// a claim of calibrated terminal-time accuracy.
    pub terminal_inclination_event: Option<LocalizedTerminalInclination>,
    /// Profile feature selected by the analytic support query at this
    /// sample's retained post-step pose. Cylinder-only compatibility runs do
    /// not expose a feature index.
    pub support_source_feature: Option<usize>,
    pub reimpact_count: u32,
    pub channels: ChannelOwnership,
    pub mechanical_energy_j: f64,
    pub energy_defect_j: f64,
}

/// Unilateral branch at one accepted runner state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoupledContactBranch {
    /// The disc support point is separated from the base.
    Open,
    /// The signed-gap contact branch is active at the endpoint. A localized
    /// root may have zero force.
    Closed,
}

/// Which unilateral-contact branch begins immediately after a localized root.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactTransitionKind {
    Opening,
    Reimpact,
}

/// A deterministic bracket for one contact transition inside an accepted macro
/// step. The bracket is a numerical diagnostic, not an experimental event
/// measurement or a claim that the reduced point-contact law resolves a finite
/// contact patch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalizedContactTransition {
    pub kind: ContactTransitionKind,
    pub time_s: f64,
    pub bracket_start_s: f64,
    pub bracket_end_s: f64,
}

/// A deterministic bracket for the inclination threshold that terminates an
/// accepted macro step. It records only the reduced-model root localization;
/// it does not certify the corresponding physical Euler-disc stopping time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalizedTerminalInclination {
    pub time_s: f64,
    pub bracket_start_s: f64,
    pub bracket_end_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoupledTerminal {
    TerminalInclination,
    HorizonReached,
    NumericalRefusal {
        reason: CoupledNumericalRefusalReason,
    },
}

/// Typed reason why the reduced numerical trajectory refused to continue.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoupledNumericalRefusalReason {
    /// The declared separation-to-contact transition budget was exceeded.
    ReimpactLimitExceeded,
    /// A bracketed contact event could not be resolved into a positive-length
    /// subinterval under the fixed deterministic localization budget.
    ContactEventLocalizationFailed,
    /// Energy accounting or the reduced base state became non-finite.
    NonFiniteEnergyOrBaseState,
}

/// Deterministic, restartable reduced state.  It contains only actual state,
/// never an encoded duration or target outcome.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledCheckpoint {
    pub state: RigidBodyState,
    pub time_s: f64,
    pub base_deflection_m: f64,
    pub base_velocity_m_per_s: f64,
    pub accumulated_channel_work_j: [f64; 5],
    pub accumulated_energy_defect_j: f64,
    /// Total body-plus-base-plus-contact energy at the start of this trajectory [J].
    pub initial_total_energy_j: f64,
    /// Stable identity of factors, restart-relevant controls, initial twist, and energy.
    pub configuration_fingerprint: u64,
    /// Deterministic binding of every restartable state and ledger field.
    /// Callers must treat this as opaque: changing a public checkpoint field
    /// without resealing it is refused on restart.
    checkpoint_fingerprint: u64,
    pub reimpact_count: u32,
    pub was_in_contact: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CoupledRun {
    pub samples: Vec<CoupledSample>,
    pub checkpoint: CoupledCheckpoint,
    pub terminal: CoupledTerminal,
    pub applicability: &'static str,
    pub model_disagreement: &'static str,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CoupledError {
    InvalidInput(&'static str),
    CheckpointMismatch,
    /// A checkpoint's state or accumulated ledger was modified after the
    /// runner sealed it, so it cannot be used as a physical restart state.
    CheckpointIntegrityMismatch,
    Dynamics(String),
}

impl fmt::Display for CoupledError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for CoupledError {}

enum CoupledGeometry<'profile, 'cx> {
    Cylinder(ContactDiscGeometry),
    Profile {
        profile: &'profile ResolvedDiscProfile,
        cx: &'profile Cx<'cx>,
    },
}

impl CoupledGeometry<'_, '_> {
    fn contact(&self, pose: Pose) -> Result<(ContactGeometry, Option<usize>), CoupledError> {
        match self {
            Self::Cylinder(geometry) => contact_geometry(*geometry, pose)
                .map(|contact| (contact, None))
                .map_err(|error| CoupledError::Dynamics(error.to_string())),
            Self::Profile { profile, cx } => {
                profile_contact_geometry(&profile.chart, profile.mass_properties, pose, cx)
                    .map(|geometry| (geometry.contact, Some(geometry.support_source_feature)))
                    .map_err(|error| CoupledError::Dynamics(error.to_string()))
            }
        }
    }

    fn state_at_ground_contact(
        &self,
        orientation: UnitQuaternion,
        linear_momentum_world: Vec3,
        angular_momentum_body: Vec3,
    ) -> Result<RigidBodyState, CoupledError> {
        match self {
            Self::Cylinder(geometry) => state_at_ground_contact(
                *geometry,
                orientation,
                linear_momentum_world,
                angular_momentum_body,
            )
            .map_err(|error| CoupledError::Dynamics(error.to_string())),
            Self::Profile { profile, cx } => profile_state_at_ground_contact(
                &profile.chart,
                profile.mass_properties,
                orientation,
                linear_momentum_world,
                angular_momentum_body,
                cx,
            )
            .map_err(|error| CoupledError::Dynamics(error.to_string())),
        }
    }

    fn identity_word(&self) -> u64 {
        match self {
            Self::Cylinder(_) => 0,
            Self::Profile { profile, .. } => profile.identity.0,
        }
    }
}

const CONTACT_EVENT_BISECTION_ITERATIONS: u32 = 48;
/// Fixed deterministic interior scan used before root bisection. Four equal
/// subintervals resolve the first crossing whose excursion persists to a scan
/// node without turning the reduced runner into an unbounded adaptive search.
const CONTACT_EVENT_SCAN_SUBDIVISIONS: u32 = 4;
const CONTACT_EVENT_GAP_TOLERANCE_M: f64 = 1.0e-12;
const CONTACT_EVENT_ROOT_TOLERANCE_M: f64 = 1.0e-15;
const TERMINAL_EVENT_BISECTION_ITERATIONS: u32 = 48;
const TERMINAL_EVENT_ROOT_TOLERANCE_RAD: f64 = 1.0e-12;
const CONTACT_EVENT_MIN_PROGRESS_FRACTION: f64 = 1.0e-12;
const MAX_CONTACT_EVENTS_PER_MACRO_STEP: u32 = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ContactBranch {
    ForceOpen,
    ForceClosed,
}

#[derive(Clone, Debug)]
struct SegmentAdvance {
    state_after: RigidBodyState,
    base_deflection_m: f64,
    base_velocity_m_per_s: f64,
    interval_start_gap_m: f64,
    interval_normal_force_n: f64,
    contact_active: bool,
    channels: ChannelOwnership,
    tangential_work_j: f64,
    contact_damping_work_j: f64,
    body_mechanical_energy_j: f64,
    post_support_gap_m: f64,
    post_contact_geometry: ContactGeometry,
    post_support_source_feature: Option<usize>,
}

/// All reduced channel quantities evaluated at one explicit state/base pair.
///
/// This is intentionally an internal stage value, not a retained physical
/// receipt: the macro-step exposes the midpoint approximation only through the
/// existing averaged channel and work fields.
#[derive(Clone, Copy, Debug)]
struct SegmentEvaluation {
    velocity_world_m_per_s: Vec3,
    omega_world_rad_per_s: Vec3,
    contact_arm_world_m: Vec3,
    interval_gap_m: f64,
    normal_force_n: f64,
    contact_active: bool,
    friction_force_world_n: Vec3,
    contact_force_world_n: Vec3,
    contact_torque_world_nm: Vec3,
    gas_force_world_n: Vec3,
    gas_torque_world_nm: Vec3,
    rolling_torque_world_nm: Vec3,
    rolling_power_w: f64,
    penetration_m: f64,
    relative_normal_speed_m_per_s: f64,
}

fn transition_kind(
    branch: ContactBranch,
    start_gap_m: f64,
    end_gap_m: f64,
) -> Option<ContactTransitionKind> {
    // The branch-specific tolerance is a deterministic hysteresis band. It
    // avoids re-detecting an event at the same root on the next subinterval,
    // and makes reimpact begin from a small but positive penetration rather
    // than an exact zero-gap force evaluation.
    match branch {
        ContactBranch::ForceClosed
            if start_gap_m <= CONTACT_EVENT_GAP_TOLERANCE_M
                && end_gap_m > CONTACT_EVENT_GAP_TOLERANCE_M =>
        {
            Some(ContactTransitionKind::Opening)
        }
        ContactBranch::ForceOpen
            if start_gap_m >= -CONTACT_EVENT_GAP_TOLERANCE_M
                && end_gap_m < -CONTACT_EVENT_GAP_TOLERANCE_M =>
        {
            Some(ContactTransitionKind::Reimpact)
        }
        _ => None,
    }
}

fn reimpact_budget_exceeded(base_count: u32, additional_events: u32, limit: u32) -> bool {
    base_count
        .checked_add(additional_events)
        .map_or(true, |count| count > limit)
}

/// Localizes the first transition found by a bounded deterministic interior
/// scan, then re-evaluates the actual reduced segment evolution to bisect that
/// bracket. An excursion entirely between adjacent scan nodes remains outside
/// this reduced runner's claim boundary; the event cap applies only to
/// transitions it resolves.
fn localize_endpoint_transition(
    interval_start_s: f64,
    duration_s: f64,
    branch: ContactBranch,
    start_gap_m: f64,
    end_gap_m: f64,
    mut gap_after: impl FnMut(f64) -> Result<f64, CoupledError>,
) -> Result<Option<LocalizedContactTransition>, CoupledError> {
    if !(interval_start_s.is_finite()
        && duration_s.is_finite()
        && duration_s > 0.0
        && start_gap_m.is_finite()
        && end_gap_m.is_finite())
    {
        return Ok(None);
    }
    let mut low_s = 0.0;
    let mut low_gap_m = start_gap_m;
    let mut bracket = None;
    for scan_index in 1..=CONTACT_EVENT_SCAN_SUBDIVISIONS {
        let high_s =
            duration_s * f64::from(scan_index) / f64::from(CONTACT_EVENT_SCAN_SUBDIVISIONS);
        let high_gap_m = if scan_index == CONTACT_EVENT_SCAN_SUBDIVISIONS {
            end_gap_m
        } else {
            gap_after(high_s)?
        };
        if !high_gap_m.is_finite() {
            return Ok(None);
        }
        if let Some(kind) = transition_kind(branch, low_gap_m, high_gap_m) {
            bracket = Some((kind, low_s, high_s));
            break;
        }
        low_s = high_s;
        low_gap_m = high_gap_m;
    }
    let Some((kind, mut low_s, mut high_s)) = bracket else {
        return Ok(None);
    };
    for _ in 0..CONTACT_EVENT_BISECTION_ITERATIONS {
        let midpoint_s = 0.5 * (low_s + high_s);
        let midpoint_gap_m = gap_after(midpoint_s)?;
        if !midpoint_gap_m.is_finite() {
            return Ok(None);
        }
        let event_threshold_m = match kind {
            ContactTransitionKind::Opening => CONTACT_EVENT_GAP_TOLERANCE_M,
            ContactTransitionKind::Reimpact => -CONTACT_EVENT_GAP_TOLERANCE_M,
        };
        if (midpoint_gap_m - event_threshold_m).abs() <= CONTACT_EVENT_ROOT_TOLERANCE_M {
            low_s = midpoint_s;
            high_s = midpoint_s;
            break;
        }
        let midpoint_is_start_side = match kind {
            ContactTransitionKind::Opening => midpoint_gap_m <= event_threshold_m,
            ContactTransitionKind::Reimpact => midpoint_gap_m >= event_threshold_m,
        };
        if midpoint_is_start_side {
            low_s = midpoint_s;
        } else {
            high_s = midpoint_s;
        }
    }
    let time_s = interval_start_s + 0.5 * (low_s + high_s);
    if !time_s.is_finite() {
        return Ok(None);
    }
    Ok(Some(LocalizedContactTransition {
        kind,
        time_s,
        bracket_start_s: interval_start_s + low_s,
        bracket_end_s: interval_start_s + high_s,
    }))
}

fn inclination_rad(state: RigidBodyState) -> f64 {
    state
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))
        .z
        .clamp(-1.0, 1.0)
        .acos()
}

/// Localizes a descending terminal-inclination crossing using actual reduced
/// segment evolution. It deliberately requires endpoint bracketing: just as
/// with contact transitions, this reduced runner makes no claim to find an
/// unbracketed sub-grid crossing that returns above the threshold.
fn localize_terminal_inclination(
    interval_start_s: f64,
    duration_s: f64,
    terminal_inclination_rad: f64,
    start_inclination_rad: f64,
    end_inclination_rad: f64,
    mut inclination_after: impl FnMut(f64) -> Result<f64, CoupledError>,
) -> Result<Option<LocalizedTerminalInclination>, CoupledError> {
    if !(interval_start_s.is_finite()
        && duration_s.is_finite()
        && duration_s > 0.0
        && terminal_inclination_rad.is_finite()
        && start_inclination_rad.is_finite()
        && end_inclination_rad.is_finite())
    {
        return Ok(None);
    }
    if !(start_inclination_rad > terminal_inclination_rad
        && end_inclination_rad <= terminal_inclination_rad)
    {
        return Ok(None);
    }
    let mut low_s = 0.0;
    let mut high_s = duration_s;
    for _ in 0..TERMINAL_EVENT_BISECTION_ITERATIONS {
        let midpoint_s = 0.5 * (low_s + high_s);
        let midpoint_inclination_rad = inclination_after(midpoint_s)?;
        if !midpoint_inclination_rad.is_finite() {
            return Ok(None);
        }
        if (midpoint_inclination_rad - terminal_inclination_rad).abs()
            <= TERMINAL_EVENT_ROOT_TOLERANCE_RAD
        {
            low_s = midpoint_s;
            high_s = midpoint_s;
            break;
        }
        if midpoint_inclination_rad > terminal_inclination_rad {
            low_s = midpoint_s;
        } else {
            high_s = midpoint_s;
        }
    }
    let time_s = interval_start_s + 0.5 * (low_s + high_s);
    if !time_s.is_finite() {
        return Ok(None);
    }
    Ok(Some(LocalizedTerminalInclination {
        time_s,
        bracket_start_s: interval_start_s + low_s,
        bracket_end_s: interval_start_s + high_s,
    }))
}

/// Advances one mechanically homogeneous contact branch. A caller that finds
/// an opening/reimpact root composes two or more of these segments, so neither
/// the contact force nor the base load is applied across the wrong side of a
/// localized transition.
fn advance_segment(
    state: RigidBodyState,
    base_deflection_m: f64,
    base_velocity_m_per_s: f64,
    duration_s: f64,
    branch: ContactBranch,
    mass: MassProperties,
    factors: CoupledFactors,
    gravity: Gravity,
    integrator: &RigidBodyIntegrator,
    geometry_model: &CoupledGeometry<'_, '_>,
) -> Result<SegmentAdvance, CoupledError> {
    let start = evaluate_segment_channels(
        state,
        base_deflection_m,
        base_velocity_m_per_s,
        branch,
        mass,
        factors,
        geometry_model,
    )?;
    let half_duration_s = 0.5 * duration_s;
    let predicted_state = integrator
        .step(
            state,
            mass,
            Wrench {
                force_world: start.contact_force_world_n.add(start.gas_force_world_n),
                torque_body: state.pose().orientation().rotate_world_to_body(
                    start
                        .contact_torque_world_nm
                        .add(start.rolling_torque_world_nm)
                        .add(start.gas_torque_world_nm),
                ),
            },
            half_duration_s,
        )
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?
        .state_after;
    let start_base_acceleration_m_per_s2 = (-start.normal_force_n
        - factors.base_stiffness_n_per_m * base_deflection_m
        - factors.base_damping_n_s_per_m * base_velocity_m_per_s)
        / factors.base_effective_mass_kg;
    let predicted_base_deflection_m = base_deflection_m + base_velocity_m_per_s * half_duration_s;
    let predicted_base_velocity_m_per_s =
        base_velocity_m_per_s + start_base_acceleration_m_per_s2 * half_duration_s;
    let midpoint = evaluate_segment_channels(
        predicted_state,
        predicted_base_deflection_m,
        predicted_base_velocity_m_per_s,
        branch,
        mass,
        factors,
        geometry_model,
    )?;
    let total_force = midpoint
        .contact_force_world_n
        .add(midpoint.gas_force_world_n);
    let total_torque_world = midpoint
        .contact_torque_world_nm
        .add(midpoint.rolling_torque_world_nm)
        .add(midpoint.gas_torque_world_nm);
    let receipt = integrator
        .step(
            state,
            mass,
            Wrench {
                force_world: total_force,
                // The force and torque are sampled at the predicted midpoint.
                // `Wrench::torque_body` is expressed in the body frame, so use
                // that same midpoint orientation rather than the start frame.
                torque_body: predicted_state
                    .pose()
                    .orientation()
                    .rotate_world_to_body(total_torque_world),
            },
            duration_s,
        )
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let work = |force: Vec3, torque: Vec3| {
        (force.dot(midpoint.velocity_world_m_per_s) + torque.dot(midpoint.omega_world_rad_per_s))
            * duration_s
    };
    let channels = ChannelOwnership {
        gravity: ChannelWrench {
            force_world_n: gravity.acceleration_world().scale(factors.mass_kg),
            torque_world_nm: Vec3::ZERO,
            work_j: work(
                gravity.acceleration_world().scale(factors.mass_kg),
                Vec3::ZERO,
            ),
        },
        contact: ChannelWrench {
            force_world_n: midpoint.contact_force_world_n,
            torque_world_nm: midpoint.contact_torque_world_nm,
            work_j: work(
                midpoint.contact_force_world_n,
                midpoint.contact_torque_world_nm,
            ),
        },
        rolling: ChannelWrench {
            force_world_n: Vec3::ZERO,
            torque_world_nm: midpoint.rolling_torque_world_nm,
            work_j: -midpoint.rolling_power_w * duration_s,
        },
        base: ChannelWrench {
            force_world_n: Vec3::ZERO,
            torque_world_nm: Vec3::ZERO,
            work_j: -factors.base_damping_n_s_per_m
                * predicted_base_velocity_m_per_s.powi(2)
                * duration_s,
        },
        gas: ChannelWrench {
            force_world_n: midpoint.gas_force_world_n,
            torque_world_nm: midpoint.gas_torque_world_nm,
            work_j: work(midpoint.gas_force_world_n, midpoint.gas_torque_world_nm),
        },
    };
    let midpoint_base_acceleration_m_per_s2 = (-midpoint.normal_force_n
        - factors.base_stiffness_n_per_m * predicted_base_deflection_m
        - factors.base_damping_n_s_per_m * predicted_base_velocity_m_per_s)
        / factors.base_effective_mass_kg;
    let base_velocity_m_per_s =
        base_velocity_m_per_s + midpoint_base_acceleration_m_per_s2 * duration_s;
    let base_deflection_m = base_deflection_m + predicted_base_velocity_m_per_s * duration_s;
    let (post_geometry, post_support_source_feature) =
        geometry_model.contact(receipt.state_after.pose())?;
    Ok(SegmentAdvance {
        state_after: receipt.state_after,
        base_deflection_m,
        base_velocity_m_per_s,
        interval_start_gap_m: start.interval_gap_m,
        interval_normal_force_n: midpoint.normal_force_n,
        contact_active: midpoint.contact_active,
        channels,
        tangential_work_j: work(
            midpoint.friction_force_world_n,
            midpoint
                .contact_arm_world_m
                .cross(midpoint.friction_force_world_n),
        ),
        contact_damping_work_j: if midpoint.penetration_m > 0.0 {
            -factors.contact_damping_n_s_per_m
                * midpoint.relative_normal_speed_m_per_s.min(0.0).powi(2)
                * duration_s
        } else {
            0.0
        },
        body_mechanical_energy_j: receipt.diagnostics_after.mechanical_energy,
        post_support_gap_m: post_geometry.gap_m - base_deflection_m,
        post_contact_geometry: post_geometry,
        post_support_source_feature,
    })
}

fn evaluate_segment_channels(
    state: RigidBodyState,
    base_deflection_m: f64,
    base_velocity_m_per_s: f64,
    branch: ContactBranch,
    mass: MassProperties,
    factors: CoupledFactors,
    geometry_model: &CoupledGeometry<'_, '_>,
) -> Result<SegmentEvaluation, CoupledError> {
    let pose = state.pose();
    let velocity_world_m_per_s = state
        .center_of_mass_velocity_world(mass)
        .map_err(|error| CoupledError::Dynamics(error.to_string()))?;
    let omega_world_rad_per_s = pose.orientation().rotate_body_to_world(
        mass.angular_velocity_body_checked(state.angular_momentum_body())
            .map_err(|error| CoupledError::Dynamics(error.to_string()))?,
    );
    let (geometry, _) = geometry_model.contact(pose)?;
    let contact_arm_world_m = geometry.radius_world_m;
    let interval_gap_m = geometry.gap_m - base_deflection_m;
    let contact_velocity_world_m_per_s =
        velocity_world_m_per_s.add(omega_world_rad_per_s.cross(contact_arm_world_m));
    let relative_normal_speed_m_per_s = contact_velocity_world_m_per_s.z - base_velocity_m_per_s;
    let penetration_m = (-interval_gap_m).max(0.0);
    let contact_active = matches!(branch, ContactBranch::ForceClosed);
    let normal_force_n = if contact_active {
        unilateral_normal_force(
            penetration_m,
            relative_normal_speed_m_per_s,
            factors.contact_stiffness_n_per_m,
            factors.contact_damping_n_s_per_m,
        )
    } else {
        0.0
    };
    let tangent_velocity_world_m_per_s = Vec3::new(
        contact_velocity_world_m_per_s.x,
        contact_velocity_world_m_per_s.y,
        0.0,
    );
    let tangent_speed_m_per_s = tangent_velocity_world_m_per_s.norm_squared().sqrt();
    let friction_force_world_n = if contact_active && tangent_speed_m_per_s > 1.0e-14 {
        tangent_velocity_world_m_per_s
            .scale(-factors.sliding_friction_coefficient * normal_force_n / tangent_speed_m_per_s)
    } else {
        Vec3::ZERO
    };
    let contact_force_world_n = Vec3::new(
        friction_force_world_n.x,
        friction_force_world_n.y,
        normal_force_n,
    );
    let contact_torque_world_nm = contact_arm_world_m.cross(contact_force_world_n);
    let gas_force_world_n =
        velocity_world_m_per_s.scale(-factors.gas_translation_damping_n_s_per_m);
    let gas_torque_world_nm = omega_world_rad_per_s.scale(-factors.gas_rotational_damping_n_m_s);
    let normal_world = pose
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let declared_rolling_power_w = factors.rolling_resistance_m
        * normal_force_n
        * (omega_world_rad_per_s.dot(normal_world).abs()
            + tangent_speed_m_per_s / factors.radius_m);
    let omega_squared = omega_world_rad_per_s.norm_squared();
    let rolling_power_w = if omega_squared > 1.0e-28 {
        declared_rolling_power_w
    } else {
        0.0
    };
    let rolling_torque_world_nm = if rolling_power_w > 0.0 {
        omega_world_rad_per_s.scale(-rolling_power_w / omega_squared)
    } else {
        Vec3::ZERO
    };
    Ok(SegmentEvaluation {
        velocity_world_m_per_s,
        omega_world_rad_per_s,
        contact_arm_world_m,
        interval_gap_m,
        normal_force_n,
        contact_active,
        friction_force_world_n,
        contact_force_world_n,
        contact_torque_world_nm,
        gas_force_world_n,
        gas_torque_world_nm,
        rolling_torque_world_nm,
        rolling_power_w,
        penetration_m,
        relative_normal_speed_m_per_s,
    })
}

fn accumulate_channel_wrench(
    total: &mut ChannelWrench,
    contribution: ChannelWrench,
    duration_s: f64,
) {
    total.force_world_n = total
        .force_world_n
        .add(contribution.force_world_n.scale(duration_s));
    total.torque_world_nm = total
        .torque_world_nm
        .add(contribution.torque_world_nm.scale(duration_s));
    total.work_j += contribution.work_j;
}

fn accumulate_channels(
    total: &mut ChannelOwnership,
    contribution: ChannelOwnership,
    duration_s: f64,
) {
    accumulate_channel_wrench(&mut total.gravity, contribution.gravity, duration_s);
    accumulate_channel_wrench(&mut total.contact, contribution.contact, duration_s);
    accumulate_channel_wrench(&mut total.rolling, contribution.rolling, duration_s);
    accumulate_channel_wrench(&mut total.base, contribution.base, duration_s);
    accumulate_channel_wrench(&mut total.gas, contribution.gas, duration_s);
}

fn average_channel_wrenches(channels: &mut ChannelOwnership, duration_s: f64) {
    for channel in [
        &mut channels.gravity,
        &mut channels.contact,
        &mut channels.rolling,
        &mut channels.base,
        &mut channels.gas,
    ] {
        channel.force_world_n = channel.force_world_n.scale(1.0 / duration_s);
        channel.torque_world_nm = channel.torque_world_nm.scale(1.0 / duration_s);
    }
}

/// Starts or resumes a closed reduced trajectory.  Contact is unilateral and
/// may separate/reimpact; rolling/base/gas are recomputed from each new state.
pub fn run_closed_reduced(
    factors: CoupledFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    restart: Option<CoupledCheckpoint>,
) -> Result<CoupledRun, CoupledError> {
    let geometry = CoupledGeometry::Cylinder(ContactDiscGeometry {
        radius_m: factors.radius_m,
        thickness_m: factors.thickness_m,
        mass_kg: factors.mass_kg,
    });
    run_closed_with_geometry(factors, controls, initial, restart, geometry)
}

/// Runs the same reduced channel model against a true resolved profile.
///
/// Geometry, support, mass, center of mass, and inertia all come from
/// `profile`. This removes the former cone-inertia/cylinder-contact surrogate;
/// it does not upgrade the reduced point-contact or loss laws to experimental
/// validation.
pub fn run_closed_profile_reduced(
    profile: &ResolvedDiscProfile,
    channels: CoupledChannelFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    restart: Option<CoupledCheckpoint>,
    cx: &Cx<'_>,
) -> Result<CoupledRun, CoupledError> {
    let factors = CoupledFactors {
        mass_kg: profile.mass_properties.mass,
        radius_m: profile.dimensions.outer_radius_m,
        thickness_m: profile.dimensions.thickness_m,
        transverse_inertia_kg_m2: profile.mass_properties.principal_inertia.transverse,
        axial_inertia_kg_m2: profile.mass_properties.principal_inertia.axial,
        gravity_m_per_s2: channels.gravity_m_per_s2,
        sliding_friction_coefficient: channels.sliding_friction_coefficient,
        rolling_resistance_m: channels.rolling_resistance_m,
        contact_stiffness_n_per_m: channels.contact_stiffness_n_per_m,
        contact_damping_n_s_per_m: channels.contact_damping_n_s_per_m,
        base_effective_mass_kg: channels.base_effective_mass_kg,
        base_stiffness_n_per_m: channels.base_stiffness_n_per_m,
        base_damping_n_s_per_m: channels.base_damping_n_s_per_m,
        gas_rotational_damping_n_m_s: channels.gas_rotational_damping_n_m_s,
        gas_translation_damping_n_s_per_m: channels.gas_translation_damping_n_s_per_m,
    };
    run_closed_with_geometry(
        factors,
        controls,
        initial,
        restart,
        CoupledGeometry::Profile { profile, cx },
    )
}

fn run_closed_with_geometry(
    factors: CoupledFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    restart: Option<CoupledCheckpoint>,
    geometry_model: CoupledGeometry<'_, '_>,
) -> Result<CoupledRun, CoupledError> {
    validate(factors, controls, initial)?;
    let mass = MassProperties::new(
        factors.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            factors.transverse_inertia_kg_m2,
            factors.transverse_inertia_kg_m2,
            factors.axial_inertia_kg_m2,
        ),
    )
    .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -factors.gravity_m_per_s2))
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let integrator = RigidBodyIntegrator::new(gravity);
    let mut configured_checkpoint = initial_checkpoint(factors, initial, mass, &geometry_model)?;
    let configured_initial_energy_j = total_energy(
        integrator
            .diagnostics(configured_checkpoint.state, mass)
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?
            .mechanical_energy,
        factors,
        configured_checkpoint.base_deflection_m,
        configured_checkpoint.base_velocity_m_per_s,
        initial_penetration_m(&geometry_model, &configured_checkpoint)?,
    );
    let configured_fingerprint = configuration_fingerprint(
        factors,
        controls,
        initial,
        configured_initial_energy_j,
        geometry_model.identity_word(),
    );
    let mut checkpoint = match restart {
        Some(checkpoint)
            if checkpoint.initial_total_energy_j.to_bits()
                == configured_initial_energy_j.to_bits()
                && checkpoint.configuration_fingerprint == configured_fingerprint =>
        {
            if checkpoint.checkpoint_fingerprint != checkpoint_fingerprint(&checkpoint) {
                return Err(CoupledError::CheckpointIntegrityMismatch);
            }
            checkpoint
        }
        Some(_) => return Err(CoupledError::CheckpointMismatch),
        None => {
            configured_checkpoint.initial_total_energy_j = configured_initial_energy_j;
            configured_checkpoint.configuration_fingerprint = configured_fingerprint;
            seal_checkpoint(&mut configured_checkpoint);
            configured_checkpoint
        }
    };
    let mut samples = Vec::with_capacity(controls.maximum_steps as usize);
    let mut previous_precession = qois(checkpoint.state, mass)?.1;

    for _ in 0..controls.maximum_steps {
        let state = checkpoint.state;
        let inclination = inclination_rad(state);
        if inclination <= controls.terminal_inclination_rad {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::TerminalInclination,
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
        let mut state = checkpoint.state;
        let mut base_deflection_m = checkpoint.base_deflection_m;
        let mut base_velocity_m_per_s = checkpoint.base_velocity_m_per_s;
        let mut remaining_s = controls.timestep_s;
        let mut elapsed_s = 0.0;
        // Restart carries the accepted branch state, so a resumed run cannot
        // reinterpret a tolerance-band contact boundary as the opposite mode.
        let mut branch = if checkpoint.was_in_contact {
            ContactBranch::ForceClosed
        } else {
            ContactBranch::ForceOpen
        };
        let mut localized_transitions =
            Vec::with_capacity(MAX_CONTACT_EVENTS_PER_MACRO_STEP as usize);
        let mut events_in_macro_step = 0_u32;
        let mut reimpact_events = 0_u32;
        let mut interval_start_gap_m = None;
        let mut interval_normal_force_n = None;
        let mut contact_active = false;
        let mut channels = ChannelOwnership::default();
        let mut channel_work_j = [0.0; 5];
        let mut body_mechanical_energy_j = None;
        let mut post_support_gap_m = None;
        let mut post_contact_geometry = None;
        let mut post_support_source_feature = None;
        let mut terminal_inclination_event = None;

        while remaining_s > 0.0 {
            let start_inclination_rad = inclination_rad(state);
            if start_inclination_rad <= controls.terminal_inclination_rad {
                let time_s = checkpoint.time_s + elapsed_s;
                terminal_inclination_event = Some(LocalizedTerminalInclination {
                    time_s,
                    bracket_start_s: time_s,
                    bracket_end_s: time_s,
                });
                break;
            }
            let trial = advance_segment(
                state,
                base_deflection_m,
                base_velocity_m_per_s,
                remaining_s,
                branch,
                mass,
                factors,
                gravity,
                &integrator,
                &geometry_model,
            )?;
            let terminal_event = localize_terminal_inclination(
                checkpoint.time_s + elapsed_s,
                remaining_s,
                controls.terminal_inclination_rad,
                start_inclination_rad,
                inclination_rad(trial.state_after),
                |candidate_duration_s| {
                    Ok(inclination_rad(
                        advance_segment(
                            state,
                            base_deflection_m,
                            base_velocity_m_per_s,
                            candidate_duration_s,
                            branch,
                            mass,
                            factors,
                            gravity,
                            &integrator,
                            &geometry_model,
                        )?
                        .state_after,
                    ))
                },
            )?;
            let event = localize_endpoint_transition(
                checkpoint.time_s + elapsed_s,
                remaining_s,
                branch,
                trial.interval_start_gap_m,
                trial.post_support_gap_m,
                |candidate_duration_s| {
                    Ok(advance_segment(
                        state,
                        base_deflection_m,
                        base_velocity_m_per_s,
                        candidate_duration_s,
                        branch,
                        mass,
                        factors,
                        gravity,
                        &integrator,
                        &geometry_model,
                    )?
                    .post_support_gap_m)
                },
            )?;
            let terminal_precedes_contact = terminal_event.is_some_and(|terminal| {
                event
                    .as_ref()
                    .is_none_or(|contact| terminal.time_s <= contact.time_s)
            });
            let (accepted, accepted_duration_s, next_branch) = if terminal_precedes_contact {
                let Some(terminal) = terminal_event else {
                    return Ok(CoupledRun {
                        samples,
                        checkpoint,
                        terminal: CoupledTerminal::NumericalRefusal {
                            reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                        },
                        applicability: applicability(),
                        model_disagreement: disagreement(),
                    });
                };
                let event_duration_s = terminal.time_s - (checkpoint.time_s + elapsed_s);
                if !(event_duration_s.is_finite()
                    && event_duration_s > controls.timestep_s * CONTACT_EVENT_MIN_PROGRESS_FRACTION
                    && event_duration_s <= remaining_s)
                {
                    return Ok(CoupledRun {
                        samples,
                        checkpoint,
                        terminal: CoupledTerminal::NumericalRefusal {
                            reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                        },
                        applicability: applicability(),
                        model_disagreement: disagreement(),
                    });
                }
                terminal_inclination_event = Some(terminal);
                (
                    advance_segment(
                        state,
                        base_deflection_m,
                        base_velocity_m_per_s,
                        event_duration_s,
                        branch,
                        mass,
                        factors,
                        gravity,
                        &integrator,
                        &geometry_model,
                    )?,
                    event_duration_s,
                    branch,
                )
            } else if let Some(event) = event {
                if events_in_macro_step >= MAX_CONTACT_EVENTS_PER_MACRO_STEP {
                    return Ok(CoupledRun {
                        samples,
                        checkpoint,
                        terminal: CoupledTerminal::NumericalRefusal {
                            reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                        },
                        applicability: applicability(),
                        model_disagreement: disagreement(),
                    });
                }
                let event_duration_s = event.time_s - (checkpoint.time_s + elapsed_s);
                if !(event_duration_s.is_finite()
                    && event_duration_s > controls.timestep_s * CONTACT_EVENT_MIN_PROGRESS_FRACTION
                    && event_duration_s <= remaining_s)
                {
                    return Ok(CoupledRun {
                        samples,
                        checkpoint,
                        terminal: CoupledTerminal::NumericalRefusal {
                            reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                        },
                        applicability: applicability(),
                        model_disagreement: disagreement(),
                    });
                }
                // Refuse at the next prohibited reimpact boundary. The
                // pre-event segment is retained in the restart checkpoint;
                // no post-event closed mechanics is committed past the budget.
                if event.kind == ContactTransitionKind::Reimpact
                    && reimpact_budget_exceeded(
                        checkpoint.reimpact_count,
                        reimpact_events + 1,
                        controls.reimpact_limit,
                    )
                {
                    debug_assert_eq!(event.kind, ContactTransitionKind::Reimpact);
                    let Some(committed_reimpact_count) =
                        checkpoint.reimpact_count.checked_add(reimpact_events)
                    else {
                        return Ok(CoupledRun {
                            samples,
                            checkpoint,
                            terminal: CoupledTerminal::NumericalRefusal {
                                reason: CoupledNumericalRefusalReason::ReimpactLimitExceeded,
                            },
                            applicability: applicability(),
                            model_disagreement: disagreement(),
                        });
                    };
                    let accepted = advance_segment(
                        state,
                        base_deflection_m,
                        base_velocity_m_per_s,
                        event_duration_s,
                        branch,
                        mass,
                        factors,
                        gravity,
                        &integrator,
                        &geometry_model,
                    )?;
                    channel_work_j[0] += accepted.channels.gravity.work_j;
                    channel_work_j[1] +=
                        accepted.tangential_work_j + accepted.contact_damping_work_j;
                    channel_work_j[2] += accepted.channels.rolling.work_j;
                    channel_work_j[3] += accepted.channels.base.work_j;
                    channel_work_j[4] += accepted.channels.gas.work_j;
                    checkpoint.state = accepted.state_after;
                    checkpoint.base_deflection_m = accepted.base_deflection_m;
                    checkpoint.base_velocity_m_per_s = accepted.base_velocity_m_per_s;
                    checkpoint.time_s += elapsed_s + event_duration_s;
                    checkpoint.reimpact_count = committed_reimpact_count;
                    checkpoint.was_in_contact = false;
                    checkpoint.accumulated_channel_work_j[0] += channel_work_j[0];
                    checkpoint.accumulated_channel_work_j[1] += channel_work_j[1];
                    checkpoint.accumulated_channel_work_j[2] += channel_work_j[2];
                    checkpoint.accumulated_channel_work_j[3] += channel_work_j[3];
                    checkpoint.accumulated_channel_work_j[4] += channel_work_j[4];
                    let total_energy = total_energy(
                        accepted.body_mechanical_energy_j,
                        factors,
                        checkpoint.base_deflection_m,
                        checkpoint.base_velocity_m_per_s,
                        (-accepted.post_support_gap_m).max(0.0),
                    );
                    checkpoint.accumulated_energy_defect_j = (total_energy
                        - checkpoint.initial_total_energy_j)
                        - (checkpoint.accumulated_channel_work_j[1]
                            + checkpoint.accumulated_channel_work_j[2]
                            + checkpoint.accumulated_channel_work_j[3]
                            + checkpoint.accumulated_channel_work_j[4]);
                    seal_checkpoint(&mut checkpoint);
                    return Ok(CoupledRun {
                        samples,
                        checkpoint,
                        terminal: CoupledTerminal::NumericalRefusal {
                            reason: CoupledNumericalRefusalReason::ReimpactLimitExceeded,
                        },
                        applicability: applicability(),
                        model_disagreement: disagreement(),
                    });
                }
                let accepted = advance_segment(
                    state,
                    base_deflection_m,
                    base_velocity_m_per_s,
                    event_duration_s,
                    branch,
                    mass,
                    factors,
                    gravity,
                    &integrator,
                    &geometry_model,
                )?;
                events_in_macro_step += 1;
                if event.kind == ContactTransitionKind::Reimpact {
                    reimpact_events += 1;
                }
                localized_transitions.push(event);
                let next_branch = match event.kind {
                    ContactTransitionKind::Opening => ContactBranch::ForceOpen,
                    ContactTransitionKind::Reimpact => ContactBranch::ForceClosed,
                };
                (accepted, event_duration_s, next_branch)
            } else {
                (trial, remaining_s, branch)
            };

            if interval_start_gap_m.is_none() {
                interval_start_gap_m = Some(accepted.interval_start_gap_m);
                interval_normal_force_n = Some(accepted.interval_normal_force_n);
            }
            contact_active |= accepted.contact_active;
            accumulate_channels(&mut channels, accepted.channels, accepted_duration_s);
            channel_work_j[0] += accepted.channels.gravity.work_j;
            channel_work_j[1] += accepted.tangential_work_j + accepted.contact_damping_work_j;
            channel_work_j[2] += accepted.channels.rolling.work_j;
            channel_work_j[3] += accepted.channels.base.work_j;
            channel_work_j[4] += accepted.channels.gas.work_j;
            state = accepted.state_after;
            base_deflection_m = accepted.base_deflection_m;
            base_velocity_m_per_s = accepted.base_velocity_m_per_s;
            body_mechanical_energy_j = Some(accepted.body_mechanical_energy_j);
            post_support_gap_m = Some(accepted.post_support_gap_m);
            post_contact_geometry = Some(accepted.post_contact_geometry);
            post_support_source_feature = accepted.post_support_source_feature;
            elapsed_s += accepted_duration_s;
            remaining_s -= accepted_duration_s;
            branch = next_branch;
            if terminal_inclination_event.is_some() {
                break;
            }
        }

        let (Some(interval_start_gap_m), Some(interval_normal_force_n)) =
            (interval_start_gap_m, interval_normal_force_n)
        else {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        };
        let (Some(body_mechanical_energy_j), Some(post_support_gap_m), Some(post_contact_geometry)) = (
            body_mechanical_energy_j,
            post_support_gap_m,
            post_contact_geometry,
        ) else {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        };
        if terminal_inclination_event.is_none()
            && (elapsed_s - controls.timestep_s).abs()
                > controls.timestep_s * CONTACT_EVENT_MIN_PROGRESS_FRACTION
        {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }

        if !(elapsed_s.is_finite() && elapsed_s > 0.0) {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::ContactEventLocalizationFailed,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
        average_channel_wrenches(&mut channels, elapsed_s);
        let Some(updated_reimpact_count) = checkpoint.reimpact_count.checked_add(reimpact_events)
        else {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::ReimpactLimitExceeded,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        };
        checkpoint.state = state;
        checkpoint.base_deflection_m = base_deflection_m;
        checkpoint.base_velocity_m_per_s = base_velocity_m_per_s;
        checkpoint.time_s += elapsed_s;
        checkpoint.was_in_contact = post_support_gap_m <= 0.0;
        checkpoint.reimpact_count = updated_reimpact_count;
        checkpoint.accumulated_channel_work_j[0] += channel_work_j[0];
        checkpoint.accumulated_channel_work_j[1] += channel_work_j[1];
        checkpoint.accumulated_channel_work_j[2] += channel_work_j[2];
        checkpoint.accumulated_channel_work_j[3] += channel_work_j[3];
        checkpoint.accumulated_channel_work_j[4] += channel_work_j[4];
        let total_energy = total_energy(
            body_mechanical_energy_j,
            factors,
            checkpoint.base_deflection_m,
            checkpoint.base_velocity_m_per_s,
            (-post_support_gap_m).max(0.0),
        );
        let defect = (total_energy - checkpoint.initial_total_energy_j)
            - (checkpoint.accumulated_channel_work_j[1]
                + checkpoint.accumulated_channel_work_j[2]
                + checkpoint.accumulated_channel_work_j[3]
                + checkpoint.accumulated_channel_work_j[4]);
        checkpoint.accumulated_energy_defect_j = defect;
        seal_checkpoint(&mut checkpoint);
        let (sample_inclination, precession, spin) = qois(checkpoint.state, mass)?;
        let precession_acceleration = (precession - previous_precession) / elapsed_s;
        previous_precession = precession;
        let center_of_mass_velocity_world_m_per_s = checkpoint
            .state
            .center_of_mass_velocity_world(mass)
            .map_err(|error| CoupledError::Dynamics(error.to_string()))?;
        if !defect.is_finite()
            || !total_energy.is_finite()
            || !checkpoint.base_deflection_m.is_finite()
            || !checkpoint.base_velocity_m_per_s.is_finite()
            || !center_of_mass_velocity_world_m_per_s.is_finite()
            || !post_support_gap_m.is_finite()
        {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::NonFiniteEnergyOrBaseState,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
        samples.push(CoupledSample {
            time_s: checkpoint.time_s,
            state: checkpoint.state,
            center_of_mass_velocity_world_m_per_s,
            base_deflection_m: checkpoint.base_deflection_m,
            base_velocity_m_per_s: checkpoint.base_velocity_m_per_s,
            contact_branch: if checkpoint.was_in_contact {
                CoupledContactBranch::Closed
            } else {
                CoupledContactBranch::Open
            },
            endpoint_contact_geometry: post_contact_geometry,
            endpoint_signed_gap_m: post_support_gap_m,
            inclination_rad: sample_inclination,
            precession_rad_per_s: precession,
            spin_rad_per_s: spin,
            precession_acceleration_rad_per_s2: precession_acceleration,
            interval_start_gap_m,
            interval_normal_force_n,
            contact_active,
            contact_transitions: localized_transitions,
            terminal_inclination_event,
            support_source_feature: post_support_source_feature,
            reimpact_count: checkpoint.reimpact_count,
            channels,
            mechanical_energy_j: total_energy,
            energy_defect_j: defect,
        });
        if checkpoint.reimpact_count > controls.reimpact_limit {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal {
                    reason: CoupledNumericalRefusalReason::ReimpactLimitExceeded,
                },
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
        if terminal_inclination_event.is_some()
            || sample_inclination <= controls.terminal_inclination_rad
        {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::TerminalInclination,
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
    }
    Ok(CoupledRun {
        samples,
        checkpoint,
        terminal: CoupledTerminal::HorizonReached,
        applicability: applicability(),
        model_disagreement: disagreement(),
    })
}

fn initial_checkpoint(
    factors: CoupledFactors,
    initial: CoupledInitialState,
    mass: MassProperties,
    geometry: &CoupledGeometry<'_, '_>,
) -> Result<CoupledCheckpoint, CoupledError> {
    let orientation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), initial.inclination_rad)
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let normal = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let omega_world =
        Vec3::new(0.0, 0.0, initial.precession_rad_per_s).add(normal.scale(initial.spin_rad_per_s));
    let omega_body = orientation.rotate_world_to_body(omega_world);
    let angular = Vec3::new(
        omega_body.x * mass.principal_inertia_body().x,
        omega_body.y * mass.principal_inertia_body().y,
        omega_body.z * mass.principal_inertia_body().z,
    );
    let ground_state = geometry.state_at_ground_contact(orientation, Vec3::ZERO, angular)?;
    let contact_arm = geometry.contact(ground_state.pose())?.0.radius_world_m;
    // This is the declared rolling initial condition: the current analytic
    // profile support point has zero world velocity before the preload offset.
    let center_of_mass_velocity = omega_world.cross(contact_arm).scale(-1.0);
    let preload_force_n = factors.mass_kg * factors.gravity_m_per_s2;
    let (base_deflection_m, contact_penetration_m) =
        if factors.contact_stiffness_n_per_m > 0.0 && factors.base_stiffness_n_per_m > 0.0 {
            (
                -preload_force_n / factors.base_stiffness_n_per_m,
                preload_force_n / factors.contact_stiffness_n_per_m,
            )
        } else {
            (0.0, 0.0)
        };
    let grounded_state = geometry.state_at_ground_contact(
        orientation,
        center_of_mass_velocity.scale(factors.mass_kg),
        angular,
    )?;
    let shifted_position = grounded_state.pose().position_world().add(Vec3::new(
        0.0,
        0.0,
        base_deflection_m - contact_penetration_m,
    ));
    let shifted_pose = Pose::new(shifted_position, orientation)
        .map_err(|error| CoupledError::Dynamics(error.to_string()))?;
    Ok(CoupledCheckpoint {
        state: RigidBodyState::new(
            shifted_pose,
            center_of_mass_velocity.scale(factors.mass_kg),
            angular,
        )
        .map_err(|error| CoupledError::Dynamics(error.to_string()))?,
        time_s: 0.0,
        base_deflection_m,
        base_velocity_m_per_s: 0.0,
        accumulated_channel_work_j: [0.0; 5],
        accumulated_energy_defect_j: 0.0,
        initial_total_energy_j: 0.0,
        configuration_fingerprint: 0,
        checkpoint_fingerprint: 0,
        reimpact_count: 0,
        was_in_contact: contact_penetration_m > 0.0,
    })
}

/// Decomposes the symmetry-axis motion without treating an arbitrary horizontal
/// angular-velocity projection as precession.
pub fn qois(state: RigidBodyState, mass: MassProperties) -> Result<(f64, f64, f64), CoupledError> {
    let normal = state
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let inclination = normal.z.clamp(-1.0, 1.0).acos();
    let omega_body = mass
        .angular_velocity_body_checked(state.angular_momentum_body())
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let omega = state.pose().orientation().rotate_body_to_world(omega_body);
    let normal_dot = omega.cross(normal);
    let horizontal_norm_squared = normal.x.mul_add(normal.x, normal.y * normal.y);
    if horizontal_norm_squared <= 1.0e-14 {
        return Err(CoupledError::InvalidInput(
            "precession undefined at symmetry-axis pole",
        ));
    }
    let precession = (normal.x * normal_dot.y - normal.y * normal_dot.x) / horizontal_norm_squared;
    // Intrinsic spin is the residual axial rate after removing the declared
    // world-vertical precession component.
    Ok((
        inclination,
        precession,
        omega.dot(normal) - precession * normal.z,
    ))
}

/// Returns the QoIs of the declared initial twist before any force is applied.
pub fn initial_qois(
    factors: CoupledFactors,
    initial: CoupledInitialState,
) -> Result<(f64, f64, f64), CoupledError> {
    validate(
        factors,
        CoupledControls {
            timestep_s: 1.0,
            maximum_steps: 1,
            terminal_inclination_rad: 0.001,
            reimpact_limit: 0,
        },
        initial,
    )?;
    let mass = MassProperties::new(
        factors.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            factors.transverse_inertia_kg_m2,
            factors.transverse_inertia_kg_m2,
            factors.axial_inertia_kg_m2,
        ),
    )
    .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let geometry = ContactDiscGeometry {
        radius_m: factors.radius_m,
        thickness_m: factors.thickness_m,
        mass_kg: factors.mass_kg,
    };
    qois(
        initial_checkpoint(factors, initial, mass, &CoupledGeometry::Cylinder(geometry))?.state,
        mass,
    )
}

/// Returns the finite-cylinder rim velocity at the initialized contact point.
/// The preload changes only position, so this verifies the rolling twist before
/// the preload offset as well.
pub fn initial_contact_point_velocity(
    factors: CoupledFactors,
    initial: CoupledInitialState,
) -> Result<Vec3, CoupledError> {
    validate(
        factors,
        CoupledControls {
            timestep_s: 1.0,
            maximum_steps: 1,
            terminal_inclination_rad: 0.001,
            reimpact_limit: 0,
        },
        initial,
    )?;
    let mass = MassProperties::new(
        factors.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            factors.transverse_inertia_kg_m2,
            factors.transverse_inertia_kg_m2,
            factors.axial_inertia_kg_m2,
        ),
    )
    .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let geometry = ContactDiscGeometry {
        radius_m: factors.radius_m,
        thickness_m: factors.thickness_m,
        mass_kg: factors.mass_kg,
    };
    let geometry = CoupledGeometry::Cylinder(geometry);
    let checkpoint = initial_checkpoint(factors, initial, mass, &geometry)?;
    let contact = geometry.contact(checkpoint.state.pose())?.0;
    let velocity = checkpoint
        .state
        .center_of_mass_velocity_world(mass)
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let omega_body = mass
        .angular_velocity_body_checked(checkpoint.state.angular_momentum_body())
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let omega_world = checkpoint
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(omega_body);
    Ok(velocity.add(omega_world.cross(contact.radius_world_m)))
}

fn initial_penetration_m(
    geometry: &CoupledGeometry<'_, '_>,
    checkpoint: &CoupledCheckpoint,
) -> Result<f64, CoupledError> {
    let contact = geometry.contact(checkpoint.state.pose())?.0;
    Ok((checkpoint.base_deflection_m - contact.gap_m).max(0.0))
}

fn total_energy(
    body_mechanical_energy_j: f64,
    factors: CoupledFactors,
    base_deflection_m: f64,
    base_velocity_m_per_s: f64,
    penetration_m: f64,
) -> f64 {
    body_mechanical_energy_j
        + 0.5 * factors.base_effective_mass_kg * base_velocity_m_per_s.powi(2)
        + 0.5 * factors.base_stiffness_n_per_m * base_deflection_m.powi(2)
        + 0.5 * factors.contact_stiffness_n_per_m * penetration_m.powi(2)
}

fn unilateral_normal_force(
    penetration_m: f64,
    relative_normal_speed_m_per_s: f64,
    stiffness_n_per_m: f64,
    damping_n_s_per_m: f64,
) -> f64 {
    // Kelvin-Voigt damping belongs to the closed contact branch. Applying the
    // dashpot while the gap is open would create a non-physical attractive
    // pre-contact force and could turn harmless approach/separation into
    // artificial reimpact chatter.
    if penetration_m <= 0.0 {
        return 0.0;
    }
    (stiffness_n_per_m * penetration_m
        + damping_n_s_per_m * (-relative_normal_speed_m_per_s).max(0.0))
    .max(0.0)
}

#[cfg(test)]
mod tests {
    use fs_mbd::{RigidBodyState, Vec3};

    use super::{
        CONTACT_EVENT_GAP_TOLERANCE_M, ContactBranch, ContactTransitionKind, CoupledControls,
        CoupledError, CoupledFactors, CoupledInitialState, CoupledTerminal,
        localize_endpoint_transition, reimpact_budget_exceeded, run_closed_reduced,
        unilateral_normal_force,
    };

    fn terminal_factors() -> CoupledFactors {
        let radius_m = 0.038;
        let thickness_m = 0.006;
        let mass_kg = std::f64::consts::PI * radius_m * radius_m * thickness_m * 2680.0;
        CoupledFactors {
            mass_kg,
            radius_m,
            thickness_m,
            transverse_inertia_kg_m2: mass_kg * (3.0 * radius_m * radius_m + thickness_m.powi(2))
                / 12.0,
            axial_inertia_kg_m2: 0.5 * mass_kg * radius_m * radius_m,
            gravity_m_per_s2: 9.806_65,
            sliding_friction_coefficient: 0.0,
            rolling_resistance_m: 0.0,
            contact_stiffness_n_per_m: 8.0e4,
            contact_damping_n_s_per_m: 3.0,
            base_effective_mass_kg: 0.25,
            base_stiffness_n_per_m: 4.0e4,
            base_damping_n_s_per_m: 4.0,
            gas_rotational_damping_n_m_s: 0.0,
            gas_translation_damping_n_s_per_m: 0.0,
        }
    }

    fn terminal_initial() -> CoupledInitialState {
        CoupledInitialState {
            inclination_rad: 0.08,
            precession_rad_per_s: 0.0,
            spin_rad_per_s: 0.0,
        }
    }

    #[test]
    fn kelvin_voigt_dashpot_is_inactive_across_an_open_gap() {
        assert_eq!(unilateral_normal_force(0.0, -3.0, 80_000.0, 3.0), 0.0);
        assert_eq!(unilateral_normal_force(-1.0e-6, -3.0, 80_000.0, 3.0), 0.0);
        assert!(unilateral_normal_force(1.0e-6, -3.0, 80_000.0, 3.0) > 0.0);
    }

    #[test]
    fn localized_opening_time_is_stable_across_equivalent_time_windows() {
        let gap_at = |time_s: f64| 1.0e-10 * (time_s - 0.3);
        let full = localize_endpoint_transition(
            0.0,
            1.0,
            ContactBranch::ForceClosed,
            gap_at(0.0),
            gap_at(1.0),
            |duration_s| Ok::<_, CoupledError>(gap_at(duration_s)),
        )
        .expect("finite gap")
        .expect("opening");
        let narrowed_start_s = 0.2;
        let narrowed = localize_endpoint_transition(
            narrowed_start_s,
            0.2,
            ContactBranch::ForceClosed,
            gap_at(narrowed_start_s),
            gap_at(narrowed_start_s + 0.2),
            |duration_s| Ok::<_, CoupledError>(gap_at(narrowed_start_s + duration_s)),
        )
        .expect("finite gap")
        .expect("opening");
        assert_eq!(full.kind, ContactTransitionKind::Opening);
        assert_eq!(narrowed.kind, ContactTransitionKind::Opening);
        // Each locator stops at the declared 1e-15 m root tolerance; with
        // this 1e-10 m/s synthetic gap slope, the two independent brackets
        // may differ by at most 2e-5 s plus floating-point roundoff.
        assert!((full.time_s - narrowed.time_s).abs() <= 2.1e-5);
        assert!(full.bracket_start_s <= full.time_s);
        assert!(full.time_s <= full.bracket_end_s);
    }

    #[test]
    fn bounded_interior_scan_finds_first_hidden_opening() {
        use core::cell::Cell;

        let tolerance = CONTACT_EVENT_GAP_TOLERANCE_M;
        let gap_at = |time_s: f64| {
            if time_s <= 0.25 {
                -2.0 * tolerance + 16.0 * tolerance * time_s
            } else if time_s <= 0.75 {
                2.0 * tolerance - 8.0 * tolerance * (time_s - 0.25)
            } else {
                -2.0 * tolerance
            }
        };
        let callback_count = Cell::new(0_u32);
        let event = localize_endpoint_transition(
            0.0,
            1.0,
            ContactBranch::ForceClosed,
            gap_at(0.0),
            gap_at(1.0),
            |duration_s| {
                callback_count.set(callback_count.get() + 1);
                Ok::<_, CoupledError>(gap_at(duration_s))
            },
        )
        .expect("finite gap");
        let event = event.expect("bounded scan must retain the first opening");
        assert_eq!(event.kind, ContactTransitionKind::Opening);
        assert!((event.time_s - 0.1875).abs() <= 1.0e-4);
        assert!(event.bracket_start_s <= event.time_s);
        assert!(event.time_s <= event.bracket_end_s);
        assert!(callback_count.get() > 0);
    }

    #[test]
    fn reimpact_budget_refuses_overflow_without_wrapping() {
        assert!(reimpact_budget_exceeded(u32::MAX, 1, u32::MAX));
        assert!(reimpact_budget_exceeded(7, 1, 7));
        assert!(!reimpact_budget_exceeded(6, 1, 7));
    }

    #[test]
    fn restart_refuses_modified_state_base_or_branch() {
        let factors = terminal_factors();
        let controls = CoupledControls {
            timestep_s: 1.0e-3,
            maximum_steps: 1,
            terminal_inclination_rad: 0.002,
            reimpact_limit: 8,
        };
        let initial = terminal_initial();
        let seed = run_closed_reduced(factors, controls, initial, None).expect("seed run");
        let mut state_mutation = seed.checkpoint.clone();
        state_mutation.state = RigidBodyState::new(
            state_mutation.state.pose(),
            state_mutation
                .state
                .linear_momentum_world()
                .add(Vec3::new(1.0e-9, 0.0, 0.0)),
            state_mutation.state.angular_momentum_body(),
        )
        .expect("finite mutation");
        let mut base_mutation = seed.checkpoint.clone();
        base_mutation.base_deflection_m += 1.0e-9;
        let mut branch_mutation = seed.checkpoint;
        branch_mutation.was_in_contact = !branch_mutation.was_in_contact;
        for mutated in [state_mutation, base_mutation, branch_mutation] {
            assert_eq!(
                run_closed_reduced(factors, controls, initial, Some(mutated)),
                Err(CoupledError::CheckpointIntegrityMismatch)
            );
        }
    }

    #[test]
    fn terminal_inclination_time_is_stable_across_coarse_and_fine_macrosteps() {
        let factors = terminal_factors();
        let initial = terminal_initial();
        let coarse_controls = CoupledControls {
            timestep_s: 1.0e-3,
            maximum_steps: 32,
            terminal_inclination_rad: 0.079,
            reimpact_limit: 8,
        };
        let fine_controls = CoupledControls {
            timestep_s: 5.0e-4,
            maximum_steps: 64,
            ..coarse_controls
        };
        let coarse = run_closed_reduced(factors, coarse_controls, initial, None)
            .expect("coarse terminal run");
        let fine =
            run_closed_reduced(factors, fine_controls, initial, None).expect("fine terminal run");
        assert_eq!(coarse.terminal, CoupledTerminal::TerminalInclination);
        assert_eq!(fine.terminal, CoupledTerminal::TerminalInclination);
        let coarse_event = coarse
            .samples
            .last()
            .and_then(|sample| sample.terminal_inclination_event)
            .expect("coarse terminal bracket");
        let fine_event = fine
            .samples
            .last()
            .and_then(|sample| sample.terminal_inclination_event)
            .expect("fine terminal bracket");
        assert!(
            coarse.checkpoint.time_s
                < coarse_controls.maximum_steps as f64 * coarse_controls.timestep_s
        );
        assert_eq!(coarse.checkpoint.time_s, coarse_event.time_s);
        assert_eq!(fine.checkpoint.time_s, fine_event.time_s);
        assert!(coarse_event.bracket_start_s <= coarse_event.time_s);
        assert!(coarse_event.time_s <= coarse_event.bracket_end_s);
        // The event locator is substep-accurate for each discrete trajectory,
        // but the trajectories themselves use a first-order rigid/base step.
        // Their observed terminal times must therefore agree within the
        // coarser 1 ms macrostep rather than a root-solver-only tolerance.
        assert!((coarse.checkpoint.time_s - fine.checkpoint.time_s).abs() <= 1.0e-3);
        assert!(
            (coarse
                .samples
                .last()
                .expect("coarse sample")
                .inclination_rad
                - coarse_controls.terminal_inclination_rad)
                .abs()
                <= 2.0e-12
        );
    }

    #[test]
    fn terminal_checkpoint_restarts_equivalently_after_a_nonterminal_prefix() {
        let factors = terminal_factors();
        let initial = terminal_initial();
        let controls = CoupledControls {
            timestep_s: 1.0e-3,
            maximum_steps: 32,
            terminal_inclination_rad: 0.079,
            reimpact_limit: 8,
        };
        let full = run_closed_reduced(factors, controls, initial, None).expect("full terminal run");
        let prefix = run_closed_reduced(
            factors,
            CoupledControls {
                maximum_steps: 1,
                ..controls
            },
            initial,
            None,
        )
        .expect("nonterminal prefix");
        assert_eq!(prefix.terminal, CoupledTerminal::HorizonReached);
        assert!(
            prefix
                .samples
                .iter()
                .all(|sample| sample.terminal_inclination_event.is_none())
        );
        let resumed = run_closed_reduced(
            factors,
            CoupledControls {
                maximum_steps: controls.maximum_steps - 1,
                ..controls
            },
            initial,
            Some(prefix.checkpoint.clone()),
        )
        .expect("restart through terminal event");
        assert_eq!(resumed.terminal, CoupledTerminal::TerminalInclination);
        let resumed_event = resumed
            .samples
            .last()
            .and_then(|sample| sample.terminal_inclination_event)
            .expect("resumed terminal bracket");
        assert!(resumed_event.time_s > prefix.checkpoint.time_s);
        assert_eq!(full.checkpoint, resumed.checkpoint);
    }
}

/// Public integration plan for replacing each reduced channel law with the
/// corresponding exported receipt adapter after focused compilation proves the
/// adapter boundary: `mechanics` supplies contact/base, `rolling_contact`
/// supplies rolling, and `air` supplies gas.  Until then, this runner makes no
/// adapter-consumption claim.
pub const ADAPTER_INTEGRATION_PLAN: &str = "after focused adapter compilation: mechanics=>contact/base;rolling_contact=>rolling;air=>gas;retain independent work ownership and energy closure";

fn applicability() -> &'static str {
    "reduced-rigid-Kelvin-Voigt-contact-and-one-mode-base;profile path uses one analytic axisymmetric chart for support and mass properties;gross-sliding-Coulomb-without-static-stick-resolution;not-finite-patch-or-resolved-gas-film"
}
fn disagreement() -> &'static str {
    "profile geometry is physical line/arc support but contact remains a reduced point Kelvin-Voigt law;rolling/base/gas are reduced channel laws and the runner does not yet consume the higher-fidelity transactional adapters"
}

fn configuration_fingerprint(
    factors: CoupledFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    initial_energy_j: f64,
    geometry_identity: u64,
) -> u64 {
    // A fixed FNV-1a mix of IEEE-754 encodings is portable and deterministic;
    // this is an identity binding, not a cryptographic digest.
    let mut fingerprint = 0xcbf2_9ce4_8422_2325_u64;
    for bits in [
        factors.mass_kg.to_bits(),
        factors.radius_m.to_bits(),
        factors.thickness_m.to_bits(),
        factors.transverse_inertia_kg_m2.to_bits(),
        factors.axial_inertia_kg_m2.to_bits(),
        factors.gravity_m_per_s2.to_bits(),
        factors.sliding_friction_coefficient.to_bits(),
        factors.rolling_resistance_m.to_bits(),
        factors.contact_stiffness_n_per_m.to_bits(),
        factors.contact_damping_n_s_per_m.to_bits(),
        factors.base_effective_mass_kg.to_bits(),
        factors.base_stiffness_n_per_m.to_bits(),
        factors.base_damping_n_s_per_m.to_bits(),
        factors.gas_rotational_damping_n_m_s.to_bits(),
        factors.gas_translation_damping_n_s_per_m.to_bits(),
        controls.timestep_s.to_bits(),
        controls.terminal_inclination_rad.to_bits(),
        u64::from(controls.reimpact_limit),
        initial.inclination_rad.to_bits(),
        initial.precession_rad_per_s.to_bits(),
        initial.spin_rad_per_s.to_bits(),
        initial_energy_j.to_bits(),
        geometry_identity,
    ] {
        fingerprint ^= bits;
        fingerprint = fingerprint.wrapping_mul(0x0000_0100_0000_01b3);
    }
    fingerprint
}

fn checkpoint_fingerprint(checkpoint: &CoupledCheckpoint) -> u64 {
    // This is an integrity binding for a restartable numerical state, not a
    // cryptographic authenticity mechanism. It prevents a caller from
    // silently grafting a different pose, base state, ledger, or branch onto
    // an otherwise valid configuration identity.
    let pose = checkpoint.state.pose();
    let position = pose.position_world();
    let orientation = pose.orientation().components();
    let linear_momentum = checkpoint.state.linear_momentum_world();
    let angular_momentum = checkpoint.state.angular_momentum_body();
    let mut fingerprint = 0xcbf2_9ce4_8422_2325_u64;
    for bits in [
        position.x.to_bits(),
        position.y.to_bits(),
        position.z.to_bits(),
        orientation[0].to_bits(),
        orientation[1].to_bits(),
        orientation[2].to_bits(),
        orientation[3].to_bits(),
        linear_momentum.x.to_bits(),
        linear_momentum.y.to_bits(),
        linear_momentum.z.to_bits(),
        angular_momentum.x.to_bits(),
        angular_momentum.y.to_bits(),
        angular_momentum.z.to_bits(),
        checkpoint.time_s.to_bits(),
        checkpoint.base_deflection_m.to_bits(),
        checkpoint.base_velocity_m_per_s.to_bits(),
        checkpoint.accumulated_channel_work_j[0].to_bits(),
        checkpoint.accumulated_channel_work_j[1].to_bits(),
        checkpoint.accumulated_channel_work_j[2].to_bits(),
        checkpoint.accumulated_channel_work_j[3].to_bits(),
        checkpoint.accumulated_channel_work_j[4].to_bits(),
        checkpoint.accumulated_energy_defect_j.to_bits(),
        checkpoint.initial_total_energy_j.to_bits(),
        checkpoint.configuration_fingerprint,
        u64::from(checkpoint.reimpact_count),
        if checkpoint.was_in_contact { 1 } else { 0 },
    ] {
        fingerprint ^= bits;
        fingerprint = fingerprint.wrapping_mul(0x0000_0100_0000_01b3);
    }
    fingerprint
}

fn seal_checkpoint(checkpoint: &mut CoupledCheckpoint) {
    checkpoint.checkpoint_fingerprint = checkpoint_fingerprint(checkpoint);
}

fn validate(
    f: CoupledFactors,
    c: CoupledControls,
    i: CoupledInitialState,
) -> Result<(), CoupledError> {
    for x in [
        f.mass_kg,
        f.radius_m,
        f.thickness_m,
        f.transverse_inertia_kg_m2,
        f.axial_inertia_kg_m2,
        f.gravity_m_per_s2,
        f.sliding_friction_coefficient,
        f.rolling_resistance_m,
        f.contact_stiffness_n_per_m,
        f.contact_damping_n_s_per_m,
        f.base_effective_mass_kg,
        f.base_stiffness_n_per_m,
        f.base_damping_n_s_per_m,
        f.gas_rotational_damping_n_m_s,
        f.gas_translation_damping_n_s_per_m,
        c.timestep_s,
        c.terminal_inclination_rad,
        i.inclination_rad,
        i.precession_rad_per_s,
        i.spin_rad_per_s,
    ] {
        if !x.is_finite() {
            return Err(CoupledError::InvalidInput("non-finite factor"));
        }
    }
    if f.mass_kg <= 0.0
        || f.radius_m <= 0.0
        || f.thickness_m <= 0.0
        || f.transverse_inertia_kg_m2 <= 0.0
        || f.axial_inertia_kg_m2 <= 0.0
        || f.gravity_m_per_s2 <= 0.0
        || f.sliding_friction_coefficient < 0.0
        || f.rolling_resistance_m < 0.0
        || f.contact_stiffness_n_per_m < 0.0
        || f.contact_damping_n_s_per_m < 0.0
        || f.base_effective_mass_kg <= 0.0
        || f.base_stiffness_n_per_m < 0.0
        || f.base_damping_n_s_per_m < 0.0
        || f.gas_rotational_damping_n_m_s < 0.0
        || f.gas_translation_damping_n_s_per_m < 0.0
        || c.timestep_s <= 0.0
        || c.maximum_steps == 0
        || c.terminal_inclination_rad <= 0.0
        || i.inclination_rad <= c.terminal_inclination_rad
    {
        return Err(CoupledError::InvalidInput("out-of-domain factor"));
    }
    Ok(())
}
