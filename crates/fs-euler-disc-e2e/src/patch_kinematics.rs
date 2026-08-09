//! Bounded pre-constitutive kinematics for one Euler-disc contact patch.
//!
//! This module turns a profile-support point and two free rigid-body point
//! velocities into an ordered relative kinematic record. It deliberately
//! contains neither a contact law nor a stick/slip decision.

use core::fmt;

use fs_couple::StableId;
use fs_mbd::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_rep_frep::AxisymmetricSupportAuthority;
use fs_tribo::InputAuthority;

use crate::contact_dynamics::ProfileContactGeometry;

/// Kinematic statuses only; no value is a constitutive contact regime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchContactStatus {
    /// Geometry/motion uncertainty leaves the bounded classification undecidable.
    Unknown,
    /// The support gap is certainly above the caller-declared separation bound.
    Separated,
    /// A near-contact pair closes more quickly than the approach threshold.
    Approaching,
    /// A near-contact pair opens faster than the stationary threshold while
    /// still inside the geometric contact envelope.
    Receding,
    /// A near-contact pair has stationary normal and tangential relative motion.
    Touching,
    /// A near-contact pair has stationary normal motion but non-stationary tangent motion.
    Grazing,
    /// A near-contact pair closes at or beyond the candidate-impact threshold.
    ImpactCandidate,
}

/// Which physical surface occupies the first slot of every relative quantity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceOrder {
    /// `relative = disc - base`; the normal points from base to disc.
    DiscThenBase,
    /// `relative = base - disc`; the normal points from disc to base.
    BaseThenDisc,
}

/// Stable identities and explicit ordering for the two contacting surfaces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrderedSurfacePair {
    first_surface: StableId,
    second_surface: StableId,
    order: SurfaceOrder,
}

impl OrderedSurfacePair {
    /// Binds an ordered pair of distinct surface identities.
    pub fn try_new(
        first_surface: StableId,
        second_surface: StableId,
        order: SurfaceOrder,
    ) -> Result<Self, PatchKinematicsError> {
        if first_surface == second_surface {
            return Err(PatchKinematicsError::IdenticalSurfaces);
        }
        Ok(Self {
            first_surface,
            second_surface,
            order,
        })
    }

    /// First surface in every `first - second` output.
    #[must_use]
    pub fn first_surface(&self) -> &StableId {
        &self.first_surface
    }

    /// Second surface in every `first - second` output.
    #[must_use]
    pub fn second_surface(&self) -> &StableId {
        &self.second_surface
    }

    /// Declared relation between the named surfaces and the disc/base inputs.
    #[must_use]
    pub const fn order(&self) -> SurfaceOrder {
        self.order
    }
}

/// Profile-support data retained at one disc material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProfileSupportKinematics {
    /// Selected disc contact arm from the true center of mass, in world metres.
    pub disc_arm_world_m: Vec3,
    /// Selected disc material point in world metres.
    pub disc_point_world_m: Vec3,
    /// Positive separation gap in metres, as supplied by profile support.
    pub gap_m: f64,
    /// Analytic meridian feature selected by the profile support query.
    pub source_feature: usize,
    /// Support-query authority. It is retained and never promoted here.
    pub support_authority: AxisymmetricSupportAuthority,
}

impl ProfileSupportKinematics {
    /// Retains exactly the support fields used by the profile-contact rung.
    #[must_use]
    pub fn from_profile_contact_geometry(geometry: ProfileContactGeometry) -> Self {
        Self {
            disc_arm_world_m: geometry.contact.radius_world_m,
            disc_point_world_m: geometry.contact.point_world_m,
            gap_m: geometry.contact.gap_m,
            source_feature: geometry.support_source_feature,
            support_authority: geometry.support_authority,
        }
    }
}

/// Curvature metadata supplied by a geometry owner, never inferred from speed.
#[derive(Clone, Debug, PartialEq)]
pub enum CurvatureMetadata {
    /// Two finite principal curvatures and an absolute uncertainty bound.
    Known {
        /// Stable identity of the curvature query/result.
        curvature_identity: StableId,
        /// Geometry owner's explicit authority; downstream laws may not promote it.
        authority: InputAuthority,
        /// First principal curvature in m⁻¹.
        first_principal_m_inverse: f64,
        /// Second principal curvature in m⁻¹.
        second_principal_m_inverse: f64,
        /// Non-negative absolute curvature uncertainty in m⁻¹.
        uncertainty_m_inverse: f64,
    },
    /// A geometry owner explicitly withheld curvature rather than guessing it.
    Unavailable {
        /// Stable identity of the unavailable curvature query.
        curvature_identity: StableId,
        /// Stable reason/capability identity for the refusal.
        reason_identity: StableId,
    },
}

/// Patch identity, selected feature, and geometry uncertainty.
#[derive(Clone, Debug, PartialEq)]
pub struct PatchGeometryMetadata {
    /// Stable patch identity.
    pub patch_identity: StableId,
    /// Support feature that must match the retained profile support result.
    pub source_feature: usize,
    /// Non-negative absolute support-gap uncertainty in metres.
    pub gap_uncertainty_m: f64,
    /// Retained curvature identity/value/refusal.
    pub curvature: CurvatureMetadata,
}

/// Caller-bound thresholds and equality tie-break identities.
#[derive(Clone, Debug, PartialEq)]
pub struct PatchKinematicThresholds {
    /// Stable identity of this complete threshold declaration.
    pub threshold_identity: StableId,
    /// Stable identity documenting boundary comparisons used below.
    pub tie_break_identity: StableId,
    /// A gap strictly above this value is certainly separated [m].
    pub separation_gap_m: f64,
    /// Largest near-contact absolute gap considered by velocity classification [m].
    pub touching_gap_m: f64,
    /// Maximum state/profile support-point reconstruction disagreement [m].
    pub support_point_coincidence_tolerance_m: f64,
    /// Maximum tangent separation between the declared disc/base counterpart points [m].
    pub tangent_counterpart_coincidence_tolerance_m: f64,
    /// Normal closure speed above this magnitude is approaching [m/s].
    pub approach_speed_m_per_s: f64,
    /// Normal closure speed at or above this magnitude is an impact candidate [m/s].
    pub impact_candidate_speed_m_per_s: f64,
    /// Tangential speed at or below this bound is stationary for `Touching` [m/s].
    pub stationary_tangent_speed_m_per_s: f64,
    /// Creepage is unavailable at or below this rolling-speed bound [m/s].
    pub minimum_reference_rolling_speed_m_per_s: f64,
    /// A projected gauge vector at or below this norm uses deterministic fallback.
    pub gauge_degeneracy_norm: f64,
}

impl PatchKinematicThresholds {
    fn validate(&self) -> Result<(), PatchKinematicsError> {
        let finite_nonnegative = [
            self.separation_gap_m,
            self.touching_gap_m,
            self.support_point_coincidence_tolerance_m,
            self.tangent_counterpart_coincidence_tolerance_m,
            self.approach_speed_m_per_s,
            self.stationary_tangent_speed_m_per_s,
        ];
        if finite_nonnegative
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || !self.impact_candidate_speed_m_per_s.is_finite()
            || self.impact_candidate_speed_m_per_s <= self.approach_speed_m_per_s
            || !self.minimum_reference_rolling_speed_m_per_s.is_finite()
            || self.minimum_reference_rolling_speed_m_per_s <= 0.0
            || !self.gauge_degeneracy_norm.is_finite()
            || self.gauge_degeneracy_norm <= 0.0
            || self.touching_gap_m > self.separation_gap_m
        {
            return Err(PatchKinematicsError::InvalidThresholds);
        }
        Ok(())
    }
}

/// Caller-selected tangent gauge with a deterministic in-plane rotation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentGaugeInput {
    /// Reference direction in world coordinates before projection into the tangent plane.
    pub reference_world: Vec3,
    /// Counterclockwise rotation around the ordered normal in radians.
    pub rotation_rad: f64,
}

/// Source of the tangent basis before the declared SO(2) rotation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TangentGaugeSource {
    /// The projected caller reference was non-degenerate.
    CallerReference,
    /// The caller reference was degenerate; the least aligned world axis was selected.
    DeterministicFallback,
}

/// A right-handed, ordered world tangent basis: `first cross second = normal`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentBasis {
    /// Ordered unit normal from second surface to first surface.
    pub normal_world: Vec3,
    /// First unit tangent.
    pub first_world: Vec3,
    /// Second unit tangent.
    pub second_world: Vec3,
    /// Whether caller reference or deterministic fallback formed the base gauge.
    pub source: TangentGaugeSource,
}

/// Coordinates of a world tangent vector in one [`TangentBasis`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentComponents {
    /// Component along `TangentBasis::first_world`.
    pub first: f64,
    /// Component along `TangentBasis::second_world`.
    pub second: f64,
}

impl TangentComponents {
    /// Reconstructs the world tangent vector in this gauge.
    #[must_use]
    pub fn reconstruct(self, basis: TangentBasis) -> Vec3 {
        basis
            .first_world
            .scale(self.first)
            .add(basis.second_world.scale(self.second))
    }

    /// Squared tangent magnitude in the declared physical units.
    #[must_use]
    pub fn squared_norm(self) -> f64 {
        self.first.mul_add(self.first, self.second * self.second)
    }
}

/// Dimensionless creepage relative to the declared rolling/entrainment convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Creepage {
    /// Relative tangent velocity divided componentwise by reference rolling speed.
    Available {
        /// First tangent creepage, dimensionless.
        longitudinal: f64,
        /// Second tangent creepage, dimensionless.
        lateral: f64,
        /// Tangential norm of `(disc_velocity + base_velocity) / 2` [m/s].
        reference_rolling_speed_m_per_s: f64,
    },
    /// Reference speed is too small for normalized creepage.
    Unavailable {
        /// Retained finite reference rolling speed [m/s].
        reference_rolling_speed_m_per_s: f64,
        /// Caller-declared lower bound [m/s].
        minimum_reference_rolling_speed_m_per_s: f64,
    },
}

/// All inputs needed for one ordered, pre-constitutive patch query.
#[derive(Clone, Debug, PartialEq)]
pub struct PatchKinematicsInput {
    /// Stable surface identities and the required subtraction order.
    pub surfaces: OrderedSurfacePair,
    /// Unit world normal directed from the second surface to the first surface.
    pub normal_world: Vec3,
    /// Profile-selected disc support point and support authority.
    pub profile_support: ProfileSupportKinematics,
    /// Patch and curvature metadata retained alongside this query.
    pub patch: PatchGeometryMetadata,
    /// Disc state and centroidal mass/inertia used to reconstruct its material velocity.
    pub disc_state: RigidBodyState,
    /// Disc centroidal mass/inertia.
    pub disc_mass_properties: MassProperties,
    /// Base state and centroidal mass/inertia used to reconstruct its material velocity.
    pub base_state: RigidBodyState,
    /// Base centroidal mass/inertia.
    pub base_mass_properties: MassProperties,
    /// Base contact arm from its center of mass in the base body frame [m].
    pub base_contact_arm_body_m: Vec3,
    /// Caller tangent gauge and its explicit rotation.
    pub tangent_gauge: TangentGaugeInput,
    /// Thresholds and equality tie-break identities.
    pub thresholds: PatchKinematicThresholds,
    /// Optional caller-supplied tangent effort probe [N] used only for the
    /// gauge-invariant scalar `tangential_power_w`; it is not a force law.
    pub tangent_effort_probe_world_n: Option<Vec3>,
}

/// Kinematics of one scalar, vertical flexible-base mode at a contact point.
///
/// This is intentionally not a rigid-body state: it carries no invented base
/// mass, orientation, angular velocity, or material arm.  It is only the
/// moving-one-mode base state that the reduced Euler runner actually owns.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingOneModeBaseState {
    /// Contact point in the undeformed base reference configuration [m].
    pub undeformed_contact_point_world_m: Vec3,
    /// Signed vertical displacement of that material point [m].
    pub vertical_displacement_m: f64,
    /// Signed vertical material-point velocity [m/s].
    pub vertical_velocity_m_per_s: f64,
}

/// Inputs for the geometry-to-patch bridge used by the moving-one-mode base.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingOneModePatchBridgeInput {
    /// Profile-selected disc support point and support authority.
    pub profile_support: ProfileSupportKinematics,
    /// Disc state used to reconstruct the selected material-point velocity.
    pub disc_state: RigidBodyState,
    /// Disc centroidal mass/inertia.
    pub disc_mass_properties: MassProperties,
    /// Explicit scalar base-mode state; never reinterpreted as a rigid body.
    pub base_mode: MovingOneModeBaseState,
    /// Unit world normal used by the downstream patch law.
    pub normal_world: Vec3,
    /// Caller-selected tangent gauge and deterministic in-plane rotation.
    pub tangent_gauge: TangentGaugeInput,
    /// Gauge and support-point coincidence tolerances.
    pub thresholds: PatchKinematicThresholds,
}

/// Inputs that complete a moving-one-mode bridge into a full patch record.
///
/// The base point is represented directly as a material-point kinematic
/// record with zero angular velocity; it is not promoted into an invented
/// rigid body.  This is the common bridge for a coupled solver that owns one
/// flexible base mode and one rigid Euler disc.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingOneModePatchKinematicsInput {
    /// Actual profile-support and moving-base bridge inputs.
    pub bridge: MovingOneModePatchBridgeInput,
    /// Explicit ordered disc/base surface identities.
    pub surfaces: OrderedSurfacePair,
    /// Resolved profile feature and relative-gap curvature metadata.
    pub patch: PatchGeometryMetadata,
    /// Optional diagnostic effort only; it is not a force law.
    pub tangent_effort_probe_world_n: Option<Vec3>,
}

/// Direct kinematics available before any normal/tangential constitutive law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingOneModePatchBridge {
    /// Disc material-point reconstruction from the actual rigid-body state.
    pub disc_point: PointKinematics,
    /// Deformed base contact position [m].
    pub base_contact_point_world_m: Vec3,
    /// Base material-point velocity [m/s], constrained to world vertical.
    pub base_contact_velocity_world_m_per_s: Vec3,
    /// Signed normal gap, `disc_point - base_point` dotted with the normal [m].
    /// It is intentionally unclassified; no contact decision is made here.
    pub normal_gap_m: f64,
    /// Tangential residual between the two declared counterpart points [m].
    pub tangent_counterpart_residual_m: f64,
    /// Disc-minus-base material-point velocity [m/s].
    pub relative_velocity_world_m_per_s: Vec3,
    /// Signed normal relative velocity [m/s].
    pub normal_relative_velocity_m_per_s: f64,
    /// Tangential disc-minus-base relative velocity [m/s].
    pub tangential_relative_velocity_world_m_per_s: Vec3,
    /// Deterministic unit normal retained for a downstream patch law.
    pub normal_world: Vec3,
    /// Deterministic right-handed tangent basis around `normal_world`.
    pub tangent_basis: TangentBasis,
}

/// One ordered pre-constitutive patch-kinematics record.
#[derive(Clone, Debug, PartialEq)]
pub struct PatchKinematics {
    /// Echoed stable surface ordering.
    pub surfaces: OrderedSurfacePair,
    /// Retained patch/curvature metadata and uncertainty.
    pub patch: PatchGeometryMetadata,
    /// Retained support-query authority.
    pub support_authority: AxisymmetricSupportAuthority,
    /// Disc material-point reconstruction.
    pub disc_point: PointKinematics,
    /// Base material-point reconstruction.
    pub base_point: PointKinematics,
    /// Ordered tangent gauge.
    pub tangent_basis: TangentBasis,
    /// `first_material_velocity - second_material_velocity` [m/s].
    pub relative_velocity_world_m_per_s: Vec3,
    /// Ordered normal relative speed, positive when the gap opens [m/s].
    pub normal_relative_velocity_m_per_s: f64,
    /// World tangent projection of relative velocity [m/s].
    pub tangential_relative_velocity_world_m_per_s: Vec3,
    /// Tangential relative velocity in the chosen gauge [m/s].
    pub tangential_relative_velocity: TangentComponents,
    /// `(disc_material_velocity + base_material_velocity) / 2` [m/s].
    pub rolling_entrainment_velocity_world_m_per_s: Vec3,
    /// Tangential rolling/entrainment velocity [m/s].
    pub rolling_entrainment_tangent_world_m_per_s: Vec3,
    /// Tangential norm of the rolling/entrainment velocity [m/s].
    pub reference_rolling_speed_m_per_s: f64,
    /// Ordered relative angular velocity projected on the normal [rad/s].
    pub normal_spin_rad_per_s: f64,
    /// Dimensionless creepage, or an explicit zero-speed refusal.
    pub creepage: Creepage,
    /// Optional caller-probed tangent effort dot relative tangent velocity [W].
    pub tangential_power_w: Option<f64>,
    /// Bounded kinematic state only; never a constitutive regime.
    pub status: PatchContactStatus,
}

/// Refusal from a bounded pre-constitutive kinematics query.
#[derive(Clone, Debug, PartialEq)]
pub enum PatchKinematicsError {
    /// The two surfaces had the same stable identity.
    IdenticalSurfaces,
    /// A named vector/scalar was non-finite or violated its declared range.
    InvalidInput { field: &'static str },
    /// A profile support field and patch feature did not describe one query.
    SupportFeatureMismatch {
        support_feature: usize,
        patch_feature: usize,
    },
    /// Reconstructed disc point and supplied support point disagreed beyond uncertainty.
    SupportPointMismatch { residual_m: f64, tolerance_m: f64 },
    /// Disc/base points were not counterparts of the declared profile patch.
    CounterpartPointMismatch {
        /// Absolute normal separation residual in metres.
        normal_residual_m: f64,
        /// Declared normal-separation tolerance in metres.
        normal_tolerance_m: f64,
        /// Tangential separation magnitude in metres.
        tangent_residual_m: f64,
        /// Declared tangent-coincidence tolerance in metres.
        tangent_tolerance_m: f64,
    },
    /// The moving-one-mode base point is not on the normal line through the
    /// selected disc support, so it cannot represent the same patch.
    MovingOneModeCounterpartMismatch {
        /// Tangential separation magnitude in metres.
        tangent_residual_m: f64,
        /// Declared tangent-coincidence tolerance in metres.
        tangent_tolerance_m: f64,
    },
    /// A lower-level rigid-body point kinematics query refused input.
    RigidBodyRefusal { detail: DynamicsError },
    /// A finite derived result could not be represented.
    NonFiniteDerived { field: &'static str },
    /// Threshold values were internally inconsistent.
    InvalidThresholds,
}

impl fmt::Display for PatchKinematicsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for PatchKinematicsError {}

/// Bridge the actual disc rigid-body point motion to the moving-one-mode base.
///
/// It intentionally does not classify contact, derive a normal force, or make
/// the scalar base mode into a fake rigid body.  The returned point velocities
/// and deterministic basis are the smallest common input for later patch laws.
pub fn bridge_moving_one_mode_patch_kinematics(
    input: MovingOneModePatchBridgeInput,
) -> Result<MovingOneModePatchBridge, PatchKinematicsError> {
    input.thresholds.validate()?;
    validate_profile_support(input.profile_support)?;
    finite_vec(
        input.base_mode.undeformed_contact_point_world_m,
        "base_mode.undeformed_contact_point_world_m",
    )?;
    finite_scalar(
        input.base_mode.vertical_displacement_m,
        "base_mode.vertical_displacement_m",
    )?;
    finite_scalar(
        input.base_mode.vertical_velocity_m_per_s,
        "base_mode.vertical_velocity_m_per_s",
    )?;
    if !input.tangent_gauge.reference_world.is_finite()
        || !input.tangent_gauge.rotation_rad.is_finite()
    {
        return Err(PatchKinematicsError::InvalidInput {
            field: "tangent gauge",
        });
    }
    let normal_world = unit(input.normal_world, "normal_world")?;
    let tangent_basis = tangent_basis(normal_world, input.tangent_gauge, &input.thresholds)?;
    let disc_arm_body_m = input
        .disc_state
        .pose()
        .orientation()
        .rotate_world_to_body(input.profile_support.disc_arm_world_m);
    finite_vec(disc_arm_body_m, "disc_contact_arm_body_m")?;
    let disc_point = input
        .disc_state
        .point_kinematics(input.disc_mass_properties, disc_arm_body_m)
        .map_err(rigid_refusal)?;
    let support_point_residual_m = norm(
        disc_point
            .point_world
            .sub(input.profile_support.disc_point_world_m),
        "support_point_residual_m",
    )?;
    if support_point_residual_m > input.thresholds.support_point_coincidence_tolerance_m {
        return Err(PatchKinematicsError::SupportPointMismatch {
            residual_m: support_point_residual_m,
            tolerance_m: input.thresholds.support_point_coincidence_tolerance_m,
        });
    }
    let base_contact_point_world_m = checked_add(
        input.base_mode.undeformed_contact_point_world_m,
        Vec3::new(0.0, 0.0, input.base_mode.vertical_displacement_m),
        "base_mode.deformed_contact_point_world_m",
    )?;
    let base_contact_velocity_world_m_per_s =
        Vec3::new(0.0, 0.0, input.base_mode.vertical_velocity_m_per_s);
    finite_vec(
        base_contact_velocity_world_m_per_s,
        "base_mode.vertical_contact_velocity_world_m_per_s",
    )?;
    // The scalar base mode has a single named point. A tangentially displaced
    // location is not the same contact patch, even if a caller could make its
    // vertical gap appear plausible; refuse it before any contact law sees it.
    let disc_minus_base_world_m = checked_sub(
        disc_point.point_world,
        base_contact_point_world_m,
        "disc_minus_base_counterpart_position",
    )?;
    let normal_gap_m = checked_dot(
        disc_minus_base_world_m,
        normal_world,
        "moving_one_mode_normal_gap",
    )?;
    let tangent_counterpart_world_m = checked_sub(
        disc_minus_base_world_m,
        normal_world.scale(normal_gap_m),
        "moving_one_mode_tangent_counterpart",
    )?;
    let tangent_counterpart_residual_m = norm(
        tangent_counterpart_world_m,
        "moving_one_mode_tangent_counterpart_residual",
    )?;
    if tangent_counterpart_residual_m > input.thresholds.tangent_counterpart_coincidence_tolerance_m
    {
        return Err(PatchKinematicsError::MovingOneModeCounterpartMismatch {
            tangent_residual_m: tangent_counterpart_residual_m,
            tangent_tolerance_m: input.thresholds.tangent_counterpart_coincidence_tolerance_m,
        });
    }
    let relative_velocity_world_m_per_s = checked_sub(
        disc_point.point_velocity_world,
        base_contact_velocity_world_m_per_s,
        "moving_one_mode_relative_velocity",
    )?;
    let normal_relative_velocity_m_per_s = checked_dot(
        relative_velocity_world_m_per_s,
        normal_world,
        "moving_one_mode_normal_relative_velocity",
    )?;
    let tangential_relative_velocity_world_m_per_s = checked_sub(
        relative_velocity_world_m_per_s,
        normal_world.scale(normal_relative_velocity_m_per_s),
        "moving_one_mode_tangential_relative_velocity",
    )?;
    Ok(MovingOneModePatchBridge {
        disc_point,
        base_contact_point_world_m,
        base_contact_velocity_world_m_per_s,
        normal_gap_m,
        tangent_counterpart_residual_m,
        relative_velocity_world_m_per_s,
        normal_relative_velocity_m_per_s,
        tangential_relative_velocity_world_m_per_s,
        normal_world,
        tangent_basis,
    })
}

/// Complete a full pre-constitutive patch record from the moving-one-mode bridge.
///
/// This shares the same threshold classification as the rigid-base path while
/// retaining the base as an explicitly translational material point.
pub fn compute_moving_one_mode_patch_kinematics(
    input: MovingOneModePatchKinematicsInput,
) -> Result<PatchKinematics, PatchKinematicsError> {
    input.bridge.thresholds.validate()?;
    validate_patch(&input.patch, input.bridge.profile_support.source_feature)?;
    let bridge = bridge_moving_one_mode_patch_kinematics(input.bridge.clone())?;
    let base_point = PointKinematics {
        arm_body: Vec3::ZERO,
        arm_world: Vec3::ZERO,
        point_world: bridge.base_contact_point_world_m,
        center_of_mass_velocity_world: bridge.base_contact_velocity_world_m_per_s,
        angular_velocity_body: Vec3::ZERO,
        angular_velocity_world: Vec3::ZERO,
        point_velocity_world: bridge.base_contact_velocity_world_m_per_s,
    };
    let (first_velocity, second_velocity, first_angular, second_angular) =
        match input.surfaces.order {
            SurfaceOrder::DiscThenBase => (
                bridge.disc_point.point_velocity_world,
                base_point.point_velocity_world,
                bridge.disc_point.angular_velocity_world,
                base_point.angular_velocity_world,
            ),
            SurfaceOrder::BaseThenDisc => (
                base_point.point_velocity_world,
                bridge.disc_point.point_velocity_world,
                base_point.angular_velocity_world,
                bridge.disc_point.angular_velocity_world,
            ),
        };
    let relative_velocity_world_m_per_s = checked_sub(
        first_velocity,
        second_velocity,
        "moving one-mode relative velocity",
    )?;
    let normal_relative_velocity_m_per_s = checked_dot(
        relative_velocity_world_m_per_s,
        bridge.normal_world,
        "moving one-mode normal relative velocity",
    )?;
    let tangential_relative_velocity_world_m_per_s = checked_sub(
        relative_velocity_world_m_per_s,
        bridge.normal_world.scale(normal_relative_velocity_m_per_s),
        "moving one-mode tangential relative velocity",
    )?;
    let tangential_relative_velocity = tangent_components(
        tangential_relative_velocity_world_m_per_s,
        bridge.tangent_basis,
        "moving one-mode tangent components",
    )?;
    let rolling_entrainment_velocity_world_m_per_s = checked_scale(
        bridge
            .disc_point
            .point_velocity_world
            .add(base_point.point_velocity_world),
        0.5,
        "moving one-mode rolling entrainment",
    )?;
    let rolling_normal = checked_dot(
        rolling_entrainment_velocity_world_m_per_s,
        bridge.normal_world,
        "moving one-mode rolling normal",
    )?;
    let rolling_entrainment_tangent_world_m_per_s = checked_sub(
        rolling_entrainment_velocity_world_m_per_s,
        bridge.normal_world.scale(rolling_normal),
        "moving one-mode rolling tangent",
    )?;
    let reference_rolling_speed_m_per_s = norm(
        rolling_entrainment_tangent_world_m_per_s,
        "moving one-mode reference rolling speed",
    )?;
    let creepage = if reference_rolling_speed_m_per_s
        > input
            .bridge
            .thresholds
            .minimum_reference_rolling_speed_m_per_s
    {
        Creepage::Available {
            longitudinal: finite_divide(
                tangential_relative_velocity.first,
                reference_rolling_speed_m_per_s,
                "moving one-mode longitudinal creepage",
            )?,
            lateral: finite_divide(
                tangential_relative_velocity.second,
                reference_rolling_speed_m_per_s,
                "moving one-mode lateral creepage",
            )?,
            reference_rolling_speed_m_per_s,
        }
    } else {
        Creepage::Unavailable {
            reference_rolling_speed_m_per_s,
            minimum_reference_rolling_speed_m_per_s: input
                .bridge
                .thresholds
                .minimum_reference_rolling_speed_m_per_s,
        }
    };
    let tangential_power_w = input
        .tangent_effort_probe_world_n
        .map(|effort| {
            finite_vec(effort, "moving one-mode tangent effort probe")?;
            checked_dot(
                effort,
                tangential_relative_velocity_world_m_per_s,
                "moving one-mode tangential power",
            )
        })
        .transpose()?;
    let status = classify_status(
        input.bridge.profile_support.gap_m,
        input.patch.gap_uncertainty_m,
        normal_relative_velocity_m_per_s,
        tangential_relative_velocity.squared_norm().sqrt(),
        &input.bridge.thresholds,
    )?;
    let normal_spin_rad_per_s = checked_dot(
        checked_sub(
            first_angular,
            second_angular,
            "moving one-mode relative angular velocity",
        )?,
        bridge.normal_world,
        "moving one-mode normal spin",
    )?;
    Ok(PatchKinematics {
        surfaces: input.surfaces,
        patch: input.patch,
        support_authority: input.bridge.profile_support.support_authority,
        disc_point: bridge.disc_point,
        base_point,
        tangent_basis: bridge.tangent_basis,
        relative_velocity_world_m_per_s,
        normal_relative_velocity_m_per_s,
        tangential_relative_velocity_world_m_per_s,
        tangential_relative_velocity,
        rolling_entrainment_velocity_world_m_per_s,
        rolling_entrainment_tangent_world_m_per_s,
        reference_rolling_speed_m_per_s,
        normal_spin_rad_per_s,
        creepage,
        tangential_power_w,
        status,
    })
}

/// Computes one bounded, pre-constitutive Euler contact-patch record.
pub fn compute_patch_kinematics(
    input: PatchKinematicsInput,
) -> Result<PatchKinematics, PatchKinematicsError> {
    input.thresholds.validate()?;
    validate_profile_support(input.profile_support)?;
    validate_patch(&input.patch, input.profile_support.source_feature)?;
    let normal_world = unit(input.normal_world, "normal_world")?;
    if !input.tangent_gauge.reference_world.is_finite()
        || !input.tangent_gauge.rotation_rad.is_finite()
        || !input.base_contact_arm_body_m.is_finite()
    {
        return Err(PatchKinematicsError::InvalidInput {
            field: "tangent gauge or base contact arm",
        });
    }
    let tangent_basis = tangent_basis(normal_world, input.tangent_gauge, &input.thresholds)?;
    let disc_arm_body_m = input
        .disc_state
        .pose()
        .orientation()
        .rotate_world_to_body(input.profile_support.disc_arm_world_m);
    finite_vec(disc_arm_body_m, "disc_contact_arm_body_m")?;
    let disc_point = input
        .disc_state
        .point_kinematics(input.disc_mass_properties, disc_arm_body_m)
        .map_err(rigid_refusal)?;
    let support_point_residual_m = norm(
        disc_point
            .point_world
            .sub(input.profile_support.disc_point_world_m),
        "support_point_residual_m",
    )?;
    let support_point_tolerance_m = input.thresholds.support_point_coincidence_tolerance_m;
    if support_point_residual_m > support_point_tolerance_m {
        return Err(PatchKinematicsError::SupportPointMismatch {
            residual_m: support_point_residual_m,
            tolerance_m: support_point_tolerance_m,
        });
    }
    let base_point = input
        .base_state
        .point_kinematics(input.base_mass_properties, input.base_contact_arm_body_m)
        .map_err(rigid_refusal)?;
    let (first_point_world_m, second_point_world_m) = match input.surfaces.order {
        SurfaceOrder::DiscThenBase => (disc_point.point_world, base_point.point_world),
        SurfaceOrder::BaseThenDisc => (base_point.point_world, disc_point.point_world),
    };
    let counterpart_separation_world_m = checked_sub(
        first_point_world_m,
        second_point_world_m,
        "counterpart_separation",
    )?;
    let counterpart_normal_separation_m = checked_dot(
        counterpart_separation_world_m,
        normal_world,
        "counterpart_normal_separation",
    )?;
    let normal_counterpart_residual_m =
        (counterpart_normal_separation_m - input.profile_support.gap_m).abs();
    if !normal_counterpart_residual_m.is_finite() {
        return Err(PatchKinematicsError::NonFiniteDerived {
            field: "normal_counterpart_residual_m",
        });
    }
    let normal_counterpart_tolerance_m =
        input.patch.gap_uncertainty_m + input.thresholds.support_point_coincidence_tolerance_m;
    finite_scalar(
        normal_counterpart_tolerance_m,
        "normal_counterpart_tolerance_m",
    )?;
    let tangent_counterpart_separation_world_m = checked_sub(
        counterpart_separation_world_m,
        normal_world.scale(counterpart_normal_separation_m),
        "tangent_counterpart_separation",
    )?;
    let tangent_counterpart_residual_m = norm(
        tangent_counterpart_separation_world_m,
        "tangent_counterpart_residual_m",
    )?;
    if normal_counterpart_residual_m > normal_counterpart_tolerance_m
        || tangent_counterpart_residual_m
            > input.thresholds.tangent_counterpart_coincidence_tolerance_m
    {
        return Err(PatchKinematicsError::CounterpartPointMismatch {
            normal_residual_m: normal_counterpart_residual_m,
            normal_tolerance_m: normal_counterpart_tolerance_m,
            tangent_residual_m: tangent_counterpart_residual_m,
            tangent_tolerance_m: input.thresholds.tangent_counterpart_coincidence_tolerance_m,
        });
    }
    let (first_velocity, second_velocity, first_angular, second_angular) =
        match input.surfaces.order {
            SurfaceOrder::DiscThenBase => (
                disc_point.point_velocity_world,
                base_point.point_velocity_world,
                disc_point.angular_velocity_world,
                base_point.angular_velocity_world,
            ),
            SurfaceOrder::BaseThenDisc => (
                base_point.point_velocity_world,
                disc_point.point_velocity_world,
                base_point.angular_velocity_world,
                disc_point.angular_velocity_world,
            ),
        };
    let relative_velocity_world_m_per_s =
        checked_sub(first_velocity, second_velocity, "relative_velocity")?;
    let normal_relative_velocity_m_per_s = checked_dot(
        relative_velocity_world_m_per_s,
        normal_world,
        "normal_relative_velocity",
    )?;
    let tangential_relative_velocity_world_m_per_s = checked_sub(
        relative_velocity_world_m_per_s,
        normal_world.scale(normal_relative_velocity_m_per_s),
        "tangential_relative_velocity",
    )?;
    let tangential_relative_velocity = tangent_components(
        tangential_relative_velocity_world_m_per_s,
        tangent_basis,
        "tangential_relative_velocity_components",
    )?;
    let rolling_entrainment_velocity_world_m_per_s = checked_scale(
        disc_point
            .point_velocity_world
            .add(base_point.point_velocity_world),
        0.5,
        "rolling_entrainment_velocity",
    )?;
    let rolling_normal = checked_dot(
        rolling_entrainment_velocity_world_m_per_s,
        normal_world,
        "rolling_entrainment_normal",
    )?;
    let rolling_entrainment_tangent_world_m_per_s = checked_sub(
        rolling_entrainment_velocity_world_m_per_s,
        normal_world.scale(rolling_normal),
        "rolling_entrainment_tangent",
    )?;
    let reference_rolling_speed_m_per_s = norm(
        rolling_entrainment_tangent_world_m_per_s,
        "reference_rolling_speed",
    )?;
    let normal_spin_rad_per_s = checked_dot(
        checked_sub(first_angular, second_angular, "relative_angular_velocity")?,
        normal_world,
        "normal_spin",
    )?;
    let creepage = if reference_rolling_speed_m_per_s
        > input.thresholds.minimum_reference_rolling_speed_m_per_s
    {
        Creepage::Available {
            longitudinal: finite_divide(
                tangential_relative_velocity.first,
                reference_rolling_speed_m_per_s,
                "longitudinal_creepage",
            )?,
            lateral: finite_divide(
                tangential_relative_velocity.second,
                reference_rolling_speed_m_per_s,
                "lateral_creepage",
            )?,
            reference_rolling_speed_m_per_s,
        }
    } else {
        Creepage::Unavailable {
            reference_rolling_speed_m_per_s,
            minimum_reference_rolling_speed_m_per_s: input
                .thresholds
                .minimum_reference_rolling_speed_m_per_s,
        }
    };
    let tangential_power_w = input
        .tangent_effort_probe_world_n
        .map(|effort| {
            finite_vec(effort, "tangent_effort_probe_world_n")?;
            checked_dot(
                effort,
                tangential_relative_velocity_world_m_per_s,
                "tangential_power",
            )
        })
        .transpose()?;
    let status = classify_status(
        input.profile_support.gap_m,
        input.patch.gap_uncertainty_m,
        normal_relative_velocity_m_per_s,
        tangential_relative_velocity.squared_norm().sqrt(),
        &input.thresholds,
    )?;
    Ok(PatchKinematics {
        surfaces: input.surfaces,
        patch: input.patch,
        support_authority: input.profile_support.support_authority,
        disc_point,
        base_point,
        tangent_basis,
        relative_velocity_world_m_per_s,
        normal_relative_velocity_m_per_s,
        tangential_relative_velocity_world_m_per_s,
        tangential_relative_velocity,
        rolling_entrainment_velocity_world_m_per_s,
        rolling_entrainment_tangent_world_m_per_s,
        reference_rolling_speed_m_per_s,
        normal_spin_rad_per_s,
        creepage,
        tangential_power_w,
        status,
    })
}

fn validate_profile_support(support: ProfileSupportKinematics) -> Result<(), PatchKinematicsError> {
    finite_vec(support.disc_arm_world_m, "profile.disc_arm_world_m")?;
    finite_vec(support.disc_point_world_m, "profile.disc_point_world_m")?;
    finite_scalar(support.gap_m, "profile.gap_m")
}

fn validate_patch(
    patch: &PatchGeometryMetadata,
    support_feature: usize,
) -> Result<(), PatchKinematicsError> {
    if patch.source_feature != support_feature {
        return Err(PatchKinematicsError::SupportFeatureMismatch {
            support_feature,
            patch_feature: patch.source_feature,
        });
    }
    if !patch.gap_uncertainty_m.is_finite() || patch.gap_uncertainty_m < 0.0 {
        return Err(PatchKinematicsError::InvalidInput {
            field: "patch.gap_uncertainty_m",
        });
    }
    if let CurvatureMetadata::Known {
        first_principal_m_inverse,
        second_principal_m_inverse,
        uncertainty_m_inverse,
        ..
    } = &patch.curvature
        && (!first_principal_m_inverse.is_finite()
            || !second_principal_m_inverse.is_finite()
            || !uncertainty_m_inverse.is_finite()
            || *uncertainty_m_inverse < 0.0)
    {
        return Err(PatchKinematicsError::InvalidInput {
            field: "patch.curvature",
        });
    }
    Ok(())
}

fn tangent_basis(
    normal_world: Vec3,
    gauge: TangentGaugeInput,
    thresholds: &PatchKinematicThresholds,
) -> Result<TangentBasis, PatchKinematicsError> {
    let projected = gauge.reference_world.sub(normal_world.scale(checked_dot(
        gauge.reference_world,
        normal_world,
        "gauge_normal_projection",
    )?));
    let (base_first, source) =
        if norm(projected, "projected_gauge_norm")? > thresholds.gauge_degeneracy_norm {
            (
                unit(projected, "projected_gauge")?,
                TangentGaugeSource::CallerReference,
            )
        } else {
            let candidates = [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ];
            let fallback = candidates
                .into_iter()
                .min_by(|left, right| {
                    normal_world
                        .dot(*left)
                        .abs()
                        .total_cmp(&normal_world.dot(*right).abs())
                })
                .ok_or(PatchKinematicsError::InvalidInput {
                    field: "fallback basis",
                })?;
            let projected = fallback.sub(normal_world.scale(checked_dot(
                fallback,
                normal_world,
                "fallback_normal_projection",
            )?));
            (
                unit(projected, "fallback_tangent")?,
                TangentGaugeSource::DeterministicFallback,
            )
        };
    let base_second = unit(normal_world.cross(base_first), "base_second_tangent")?;
    let (sin, cos) = gauge.rotation_rad.sin_cos();
    let first_world = unit(
        base_first.scale(cos).add(base_second.scale(sin)),
        "rotated_first_tangent",
    )?;
    let second_world = unit(normal_world.cross(first_world), "rotated_second_tangent")?;
    Ok(TangentBasis {
        normal_world,
        first_world,
        second_world,
        source,
    })
}

fn classify_status(
    gap_m: f64,
    gap_uncertainty_m: f64,
    normal_velocity_m_per_s: f64,
    tangential_speed_m_per_s: f64,
    thresholds: &PatchKinematicThresholds,
) -> Result<PatchContactStatus, PatchKinematicsError> {
    finite_scalar(normal_velocity_m_per_s, "normal_velocity_m_per_s")?;
    finite_scalar(tangential_speed_m_per_s, "tangential_speed_m_per_s")?;
    let lower_gap = gap_m - gap_uncertainty_m;
    let upper_gap = gap_m + gap_uncertainty_m;
    if !lower_gap.is_finite() || !upper_gap.is_finite() {
        return Err(PatchKinematicsError::NonFiniteDerived {
            field: "gap interval",
        });
    }
    if lower_gap > thresholds.separation_gap_m {
        return Ok(PatchContactStatus::Separated);
    }
    if upper_gap > thresholds.separation_gap_m
        || lower_gap < -thresholds.touching_gap_m
        || upper_gap < -thresholds.touching_gap_m
    {
        return Ok(PatchContactStatus::Unknown);
    }
    if normal_velocity_m_per_s <= -thresholds.impact_candidate_speed_m_per_s {
        return Ok(PatchContactStatus::ImpactCandidate);
    }
    if normal_velocity_m_per_s <= -thresholds.approach_speed_m_per_s {
        return Ok(PatchContactStatus::Approaching);
    }
    if normal_velocity_m_per_s.abs() <= thresholds.approach_speed_m_per_s {
        return if tangential_speed_m_per_s <= thresholds.stationary_tangent_speed_m_per_s {
            Ok(PatchContactStatus::Touching)
        } else {
            Ok(PatchContactStatus::Grazing)
        };
    }
    Ok(PatchContactStatus::Receding)
}

fn tangent_components(
    vector: Vec3,
    basis: TangentBasis,
    field: &'static str,
) -> Result<TangentComponents, PatchKinematicsError> {
    Ok(TangentComponents {
        first: checked_dot(vector, basis.first_world, field)?,
        second: checked_dot(vector, basis.second_world, field)?,
    })
}

fn rigid_refusal(detail: DynamicsError) -> PatchKinematicsError {
    PatchKinematicsError::RigidBodyRefusal { detail }
}

fn finite_scalar(value: f64, field: &'static str) -> Result<(), PatchKinematicsError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PatchKinematicsError::InvalidInput { field })
    }
}

fn finite_vec(value: Vec3, field: &'static str) -> Result<(), PatchKinematicsError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PatchKinematicsError::InvalidInput { field })
    }
}

fn norm(value: Vec3, field: &'static str) -> Result<f64, PatchKinematicsError> {
    finite_vec(value, field)?;
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    let scaled = value.scale(scale.recip());
    let result = scale * scaled.dot(scaled).sqrt();
    if result.is_finite() {
        Ok(result)
    } else {
        Err(PatchKinematicsError::NonFiniteDerived { field })
    }
}

fn unit(value: Vec3, field: &'static str) -> Result<Vec3, PatchKinematicsError> {
    let magnitude = norm(value, field)?;
    if magnitude == 0.0 {
        return Err(PatchKinematicsError::InvalidInput { field });
    }
    let result = value.scale(magnitude.recip());
    finite_vec(result, field)?;
    Ok(result)
}

fn checked_sub(left: Vec3, right: Vec3, field: &'static str) -> Result<Vec3, PatchKinematicsError> {
    let result = left.sub(right);
    finite_vec(result, field)?;
    Ok(result)
}

fn checked_add(left: Vec3, right: Vec3, field: &'static str) -> Result<Vec3, PatchKinematicsError> {
    let result = left.add(right);
    finite_vec(result, field)?;
    Ok(result)
}

fn checked_scale(
    value: Vec3,
    scalar: f64,
    field: &'static str,
) -> Result<Vec3, PatchKinematicsError> {
    finite_scalar(scalar, field)?;
    let result = value.scale(scalar);
    finite_vec(result, field)?;
    Ok(result)
}

fn checked_dot(left: Vec3, right: Vec3, field: &'static str) -> Result<f64, PatchKinematicsError> {
    let result = left.dot(right);
    if result.is_finite() {
        Ok(result)
    } else {
        Err(PatchKinematicsError::NonFiniteDerived { field })
    }
}

fn finite_divide(
    numerator: f64,
    denominator: f64,
    field: &'static str,
) -> Result<f64, PatchKinematicsError> {
    let result = numerator / denominator;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(PatchKinematicsError::NonFiniteDerived { field })
    }
}
