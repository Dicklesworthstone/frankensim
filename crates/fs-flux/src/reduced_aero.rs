//! Deterministic, reduced aerodynamic wrenches for moving rigid bodies.
//!
//! This module is deliberately a *screening* model, not CFD.  It accepts a
//! circular-disc geometry because a spinning, tilted disc is a useful consumer,
//! but its inputs and outputs are generic moving-body quantities in one declared
//! inertial/world frame.  All scalar fields use coherent SI units in their names.
//!
//! The returned wrench is the force and moment **on the body** about its supplied
//! reference point.  The passive, stationary-ambient identity is
//! `force dot relative_velocity + torque dot angular_velocity <= 0`.  With a
//! moving ambient, that relative-power identity remains available, but a caller
//! must close a total gas-plus-boundary accounting window before making an energy
//! claim.
//!
//! Thin-gap pressure and target-fitted/video-derived terms are rejected at model
//! admission.  They belong to their own spatial models and cannot be smuggled
//! into a total drag coefficient here.

use core::fmt;
use std::collections::BTreeSet;

const PI: f64 = std::f64::consts::PI;
const UNIT_TOLERANCE: f64 = 1.0e-12;

/// A three-vector expressed in the request's explicit world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// X component in the declared frame.
    pub x: f64,
    /// Y component in the declared frame.
    pub y: f64,
    /// Z component in the declared frame.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Construct a vector.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn norm(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }

    /// Scale all components by a dimensionless scalar.
    #[must_use]
    pub fn scaled(self, scalar: f64) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }

    fn minus(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    fn plus(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

/// Typed refusal from reduced-wrench admission or evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum ReducedAeroError {
    /// A named quantity is absent from an otherwise identified card.
    MissingGasProperty(&'static str),
    /// A named scalar or vector is non-finite or outside its physical domain.
    InvalidInput {
        /// Stable field name.
        field: &'static str,
    },
    /// A source, frame, or correlation identifier is empty or non-canonical.
    InvalidIdentity {
        /// Stable field name.
        field: &'static str,
    },
    /// The disc axis was not a unit direction in the declared frame.
    NonUnitDiscAxis {
        /// Measured norm.
        norm: f64,
    },
    /// A component was listed but not supplied, or supplied but not listed.
    ContributionDeclarationMismatch {
        /// Contribution family at issue.
        family: ContributionFamily,
    },
    /// A family is excluded from this generic exterior-flow model.
    ForbiddenContribution {
        /// Rejected family.
        family: ContributionFamily,
    },
    /// A configured correlation was used beyond its declared valid range.
    OutsideCorrelationDomain {
        /// Dimensionless quantity name.
        quantity: &'static str,
        /// Computed value.
        value: f64,
        /// Inclusive lower limit.
        minimum: f64,
        /// Inclusive upper limit.
        maximum: f64,
    },
    /// The same work-exchange key was recorded twice in one window.
    DuplicateWorkExchange {
        /// Caller-owned exact-once key.
        key: u64,
    },
    /// Work duration must be finite and non-negative.
    InvalidWorkDuration,
    /// At least one explicitly identified correlation candidate is required.
    EmptyAlternativeWrenchSet,
    /// A derived aerodynamic quantity overflowed or became non-finite.
    NonFiniteDerived {
        /// Stable derived quantity name.
        field: &'static str,
    },
    /// A public candidate wrench cannot be admitted to a work transaction.
    InvalidCandidateWrench {
        /// Stable malformed candidate field name.
        field: &'static str,
    },
    /// A passive-gas candidate reports positive relative mechanical power.
    PassivePowerViolation {
        /// Reported relative mechanical power [W].
        relative_power_w: f64,
    },
    /// Alternative candidates must retain distinct complete correlation identities.
    DuplicateCorrelationIdentity {
        /// Identity duplicated within one alternative set.
        correlation: CorrelationIdentity,
    },
}

impl fmt::Display for ReducedAeroError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for ReducedAeroError {}

fn checked_identity(value: &str, field: &'static str) -> Result<(), ReducedAeroError> {
    if value.is_empty()
        || !value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(ReducedAeroError::InvalidIdentity { field });
    }
    Ok(())
}

fn finite_nonnegative(value: f64, field: &'static str) -> Result<(), ReducedAeroError> {
    if !(value.is_finite() && value >= 0.0) {
        return Err(ReducedAeroError::InvalidInput { field });
    }
    Ok(())
}

fn finite_positive(value: f64, field: &'static str) -> Result<(), ReducedAeroError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(ReducedAeroError::InvalidInput { field });
    }
    Ok(())
}

fn finite_derived(value: f64, field: &'static str) -> Result<f64, ReducedAeroError> {
    if !value.is_finite() {
        return Err(ReducedAeroError::NonFiniteDerived { field });
    }
    Ok(value)
}

fn finite_vec3(value: Vec3, field: &'static str) -> Result<Vec3, ReducedAeroError> {
    if !value.finite() {
        return Err(ReducedAeroError::NonFiniteDerived { field });
    }
    Ok(value)
}

fn checked_dot(left: Vec3, right: Vec3, field: &'static str) -> Result<f64, ReducedAeroError> {
    finite_derived(left.dot(right), field)
}

fn checked_norm(value: Vec3, field: &'static str) -> Result<f64, ReducedAeroError> {
    finite_derived(value.norm(), field)
}

fn checked_product(left: f64, right: f64, field: &'static str) -> Result<f64, ReducedAeroError> {
    finite_derived(left * right, field)
}

fn checked_sum(left: f64, right: f64, field: &'static str) -> Result<f64, ReducedAeroError> {
    finite_derived(left + right, field)
}

/// Identifies the cited correlation and version used for one candidate.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CorrelationIdentity {
    /// Stable correlation name, not a target-outcome label.
    pub id: String,
    /// Correlation revision or source edition.
    pub version: String,
    /// Source-card identity retained in every estimate receipt.
    pub source_id: String,
}

impl CorrelationIdentity {
    /// Construct a correlation identity with canonical transport-safe fields.
    pub fn try_new(
        id: impl Into<String>,
        version: impl Into<String>,
        source_id: impl Into<String>,
    ) -> Result<Self, ReducedAeroError> {
        let identity = Self {
            id: id.into(),
            version: version.into(),
            source_id: source_id.into(),
        };
        identity.validate()?;
        Ok(identity)
    }

    fn validate(&self) -> Result<(), ReducedAeroError> {
        checked_identity(&self.id, "correlation.id")?;
        checked_identity(&self.version, "correlation.version")?;
        checked_identity(&self.source_id, "correlation.source_id")?;
        Ok(())
    }
}

/// A complete-or-refused gas-property card in coherent SI units.
#[derive(Debug, Clone, PartialEq)]
pub struct GasPropertyCard {
    /// Immutable input-card identity.
    pub source_id: String,
    /// Mass density [kg/m^3].  Zero is admitted as the vacuum limiting case.
    pub density_kg_per_m3: Option<f64>,
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: Option<f64>,
    /// Speed of sound [m/s], used for a tip-Mach validity gate.
    pub speed_of_sound_m_per_s: Option<f64>,
    /// Ambient gas velocity [m/s] in the explicit world frame.
    pub velocity_world_m_per_s: Vec3,
}

/// Validated gas properties used by the exterior-flow model.
#[derive(Debug, Clone, PartialEq)]
pub struct GasProperties {
    source_id: String,
    density_kg_per_m3: f64,
    dynamic_viscosity_pa_s: f64,
    speed_of_sound_m_per_s: f64,
    velocity_world_m_per_s: Vec3,
}

impl TryFrom<GasPropertyCard> for GasProperties {
    type Error = ReducedAeroError;

    fn try_from(card: GasPropertyCard) -> Result<Self, Self::Error> {
        checked_identity(&card.source_id, "gas.source_id")?;
        let density = card
            .density_kg_per_m3
            .ok_or(ReducedAeroError::MissingGasProperty("density_kg_per_m3"))?;
        finite_nonnegative(density, "gas.density_kg_per_m3")?;
        let viscosity = card
            .dynamic_viscosity_pa_s
            .ok_or(ReducedAeroError::MissingGasProperty(
                "dynamic_viscosity_pa_s",
            ))?;
        finite_positive(viscosity, "gas.dynamic_viscosity_pa_s")?;
        let sound_speed =
            card.speed_of_sound_m_per_s
                .ok_or(ReducedAeroError::MissingGasProperty(
                    "speed_of_sound_m_per_s",
                ))?;
        finite_positive(sound_speed, "gas.speed_of_sound_m_per_s")?;
        if !card.velocity_world_m_per_s.finite() {
            return Err(ReducedAeroError::InvalidInput {
                field: "gas.velocity_world_m_per_s",
            });
        }
        Ok(Self {
            source_id: card.source_id,
            density_kg_per_m3: density,
            dynamic_viscosity_pa_s: viscosity,
            speed_of_sound_m_per_s: sound_speed,
            velocity_world_m_per_s: card.velocity_world_m_per_s,
        })
    }
}

impl GasProperties {
    /// Source identity retained in the estimate receipt.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
}

/// Circular-disc geometry.  Thickness only supplies the exterior rim area;
/// it never denotes a fluid gap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscGeometry {
    /// Disc radius [m].
    pub radius_m: f64,
    /// Exterior rim thickness [m], zero for an ideal thin exterior rim.
    pub exterior_thickness_m: f64,
}

impl DiscGeometry {
    fn validate(self) -> Result<(), ReducedAeroError> {
        finite_positive(self.radius_m, "geometry.radius_m")?;
        finite_nonnegative(self.exterior_thickness_m, "geometry.exterior_thickness_m")
    }

    fn face_area_m2(self) -> f64 {
        PI * self.radius_m * self.radius_m
    }

    /// Exterior cylindrical rim area [m^2], retained for an edge-flow torque
    /// correlation.  It is not the rim's projected drag silhouette.
    fn rim_wetted_area_m2(self) -> f64 {
        2.0 * PI * self.radius_m * self.exterior_thickness_m
    }

    /// Edge-on projected rim silhouette [m^2].
    fn edge_on_rim_silhouette_m2(self) -> f64 {
        2.0 * self.radius_m * self.exterior_thickness_m
    }
}

/// Pose information needed by the reduced exterior-flow geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscPose {
    /// Nominal unit disc normal in the declared world frame. Accepted
    /// near-unit values are projected before force evaluation.
    pub normal_world: Vec3,
}

impl DiscPose {
    /// Admit a near-unit disc normal and store its unit projection. Values
    /// outside the explicit unit tolerance are refused rather than silently
    /// normalized.
    pub fn try_new(normal_world: Vec3) -> Result<Self, ReducedAeroError> {
        let pose = Self { normal_world };
        Ok(Self {
            normal_world: pose.unit_normal()?,
        })
    }

    fn validate(self) -> Result<(), ReducedAeroError> {
        if !self.normal_world.finite() {
            return Err(ReducedAeroError::InvalidInput {
                field: "pose.normal_world",
            });
        }
        let norm = self.normal_world.norm();
        if !(norm.is_finite() && (norm - 1.0).abs() <= UNIT_TOLERANCE) {
            return Err(ReducedAeroError::NonUnitDiscAxis { norm });
        }
        Ok(())
    }

    fn unit_normal(self) -> Result<Vec3, ReducedAeroError> {
        self.validate()?;
        let norm = checked_norm(self.normal_world, "pose.normal_world_norm")?;
        finite_vec3(
            self.normal_world.scaled(1.0 / norm),
            "pose.unit_normal_world",
        )
    }
}

/// Body translational and angular velocity in the explicit world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyKinematics {
    /// Disc-center position [m] in the declared world frame. All reported
    /// wrenches are about this center; this reduced model does not shift
    /// moments to an arbitrary off-center reference point.
    pub reference_point_world_m: Vec3,
    /// Disc-center velocity [m/s].
    pub linear_velocity_world_m_per_s: Vec3,
    /// Angular velocity [rad/s] (radian is dimensionless).
    pub angular_velocity_world_rad_per_s: Vec3,
}

impl BodyKinematics {
    fn validate(self) -> Result<(), ReducedAeroError> {
        if !self.reference_point_world_m.finite() {
            return Err(ReducedAeroError::InvalidInput {
                field: "kinematics.reference_point_world_m",
            });
        }
        if !self.linear_velocity_world_m_per_s.finite() {
            return Err(ReducedAeroError::InvalidInput {
                field: "kinematics.linear_velocity_world_m_per_s",
            });
        }
        if !self.angular_velocity_world_rad_per_s.finite() {
            return Err(ReducedAeroError::InvalidInput {
                field: "kinematics.angular_velocity_world_rad_per_s",
            });
        }
        Ok(())
    }
}

/// Roughness of the exterior body surface, with a retained input-card identity.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceRoughness {
    /// Roughness card/source identity.
    pub source_id: String,
    /// Equivalent exterior roughness height [m].
    pub height_m: f64,
}

impl SurfaceRoughness {
    fn validate(&self) -> Result<(), ReducedAeroError> {
        checked_identity(&self.source_id, "roughness.source_id")?;
        finite_nonnegative(self.height_m, "roughness.height_m")
    }
}

/// Inclusive validity interval for one dimensionless quantity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedRange {
    /// Inclusive lower boundary.
    pub minimum: f64,
    /// Inclusive upper boundary.
    pub maximum: f64,
}

impl ClosedRange {
    /// Construct a finite ordered inclusive range.
    pub fn try_new(minimum: f64, maximum: f64) -> Result<Self, ReducedAeroError> {
        if !(minimum.is_finite() && maximum.is_finite() && minimum >= 0.0 && maximum >= minimum) {
            return Err(ReducedAeroError::InvalidInput {
                field: "correlation.range",
            });
        }
        Ok(Self { minimum, maximum })
    }

    fn validate(self) -> Result<(), ReducedAeroError> {
        if !(self.minimum.is_finite()
            && self.maximum.is_finite()
            && self.minimum >= 0.0
            && self.maximum >= self.minimum)
        {
            return Err(ReducedAeroError::InvalidInput {
                field: "correlation.range",
            });
        }
        Ok(())
    }

    fn require(self, quantity: &'static str, value: f64) -> Result<(), ReducedAeroError> {
        self.validate()?;
        if !(self.minimum..=self.maximum).contains(&value) {
            return Err(ReducedAeroError::OutsideCorrelationDomain {
                quantity,
                value,
                minimum: self.minimum,
                maximum: self.maximum,
            });
        }
        Ok(())
    }
}

/// Explicit validity envelope; a model refuses rather than extrapolating.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApplicabilityEnvelope {
    /// Valid Reynolds range based on disc diameter and relative translation.
    pub translational_reynolds: ClosedRange,
    /// Valid Reynolds range based on `omega * radius^2`.
    pub rotational_reynolds: ClosedRange,
    /// Valid relative roughness `height / radius` range.
    pub relative_roughness: ClosedRange,
    /// Maximum rim-tip Mach number, inclusive.
    pub maximum_tip_mach: f64,
}

/// Caller-supplied coefficient uncertainty retained without upgrading the
/// correlation to validated physical authority.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationUncertainty {
    /// Source identity for the stated uncertainty convention.
    pub source_id: String,
    /// Non-negative relative half-width applied by downstream uncertainty
    /// propagation; this reduced model does not invent a propagation theorem.
    pub coefficient_relative_half_width: f64,
}

impl CorrelationUncertainty {
    fn validate(&self) -> Result<(), ReducedAeroError> {
        checked_identity(&self.source_id, "uncertainty.source_id")?;
        finite_nonnegative(
            self.coefficient_relative_half_width,
            "uncertainty.coefficient_relative_half_width",
        )
    }
}

impl ApplicabilityEnvelope {
    fn validate(self) -> Result<(), ReducedAeroError> {
        self.translational_reynolds.validate()?;
        self.rotational_reynolds.validate()?;
        self.relative_roughness.validate()?;
        finite_nonnegative(self.maximum_tip_mach, "correlation.maximum_tip_mach")
    }
}

/// Explicit component family.  The two forbidden variants exist so an
/// admission test can prove they are rejected rather than silently ignored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ContributionFamily {
    /// Quadratic translational/form drag on projected exterior area.
    TranslationalFormDrag,
    /// Face rotational skin-friction torque.
    RotationalSkinFriction,
    /// Exterior-rim rotational edge-flow torque.
    EdgeFlow,
    /// Torque damping disc-axis reorientation (not spin).
    OrientationRateDamping,
    /// Forbidden: spatial thin-gap pressure belongs to a gas-film model.
    ThinGapPressure,
    /// Forbidden: coefficients fitted to device/video target outcomes.
    TargetFitted,
}

impl ContributionFamily {
    fn admitted(self) -> bool {
        matches!(
            self,
            Self::TranslationalFormDrag
                | Self::RotationalSkinFriction
                | Self::EdgeFlow
                | Self::OrientationRateDamping
        )
    }
}

/// Nonnegative exterior translational/form-drag coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FormDrag {
    /// Dimensionless drag coefficient.
    pub coefficient: f64,
}

/// Nonnegative face rotational skin-friction coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotationalSkinFriction {
    /// Dimensionless torque coefficient.
    pub coefficient: f64,
}

/// Nonnegative exterior-rim edge-flow coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeFlow {
    /// Dimensionless torque coefficient.
    pub coefficient: f64,
}

/// Nonnegative damping coefficient for angular velocity perpendicular to the
/// disc normal.  It does not model an aerodynamic lifting force.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientationRateDamping {
    /// Dimensionless torque coefficient.
    pub coefficient: f64,
}

/// The optional, separately visible components of one correlation candidate.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ReducedAeroComponents {
    /// Optional translational/form term.
    pub form_drag: Option<FormDrag>,
    /// Optional face rotational-skin term.
    pub rotational_skin_friction: Option<RotationalSkinFriction>,
    /// Optional exterior-rim edge-flow term.
    pub edge_flow: Option<EdgeFlow>,
    /// Optional disc-axis-reorientation damping term.
    pub orientation_rate_damping: Option<OrientationRateDamping>,
}

impl ReducedAeroComponents {
    fn families(self) -> Vec<ContributionFamily> {
        let mut families = Vec::new();
        if self.form_drag.is_some() {
            families.push(ContributionFamily::TranslationalFormDrag);
        }
        if self.rotational_skin_friction.is_some() {
            families.push(ContributionFamily::RotationalSkinFriction);
        }
        if self.edge_flow.is_some() {
            families.push(ContributionFamily::EdgeFlow);
        }
        if self.orientation_rate_damping.is_some() {
            families.push(ContributionFamily::OrientationRateDamping);
        }
        families
    }

    fn validate(self) -> Result<(), ReducedAeroError> {
        for coefficient in [
            self.form_drag.map(|term| term.coefficient),
            self.rotational_skin_friction.map(|term| term.coefficient),
            self.edge_flow.map(|term| term.coefficient),
            self.orientation_rate_damping.map(|term| term.coefficient),
        ]
        .into_iter()
        .flatten()
        {
            finite_nonnegative(coefficient, "component.coefficient")?;
        }
        Ok(())
    }
}

/// A one-correlation, reduced exterior aerodynamic candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedAeroModel {
    correlation: CorrelationIdentity,
    envelope: ApplicabilityEnvelope,
    uncertainty: CorrelationUncertainty,
    components: ReducedAeroComponents,
}

impl ReducedAeroModel {
    /// Admit an explicit exterior-flow component set.
    ///
    /// `declared_families` is intentionally redundant with `components`: it
    /// prevents a metadata path from claiming a component that the executable
    /// model did not supply, and makes forbidden thin-gap/target-fit attempts
    /// fail closed.
    pub fn try_new(
        correlation: CorrelationIdentity,
        envelope: ApplicabilityEnvelope,
        uncertainty: CorrelationUncertainty,
        components: ReducedAeroComponents,
        declared_families: &[ContributionFamily],
    ) -> Result<Self, ReducedAeroError> {
        correlation.validate()?;
        envelope.validate()?;
        uncertainty.validate()?;
        components.validate()?;
        let declared: BTreeSet<_> = declared_families.iter().copied().collect();
        for family in &declared {
            if !family.admitted() {
                return Err(ReducedAeroError::ForbiddenContribution { family: *family });
            }
        }
        let supplied: BTreeSet<_> = components.families().into_iter().collect();
        for family in declared.symmetric_difference(&supplied) {
            return Err(ReducedAeroError::ContributionDeclarationMismatch { family: *family });
        }
        if supplied.is_empty() {
            return Err(ReducedAeroError::InvalidInput {
                field: "components",
            });
        }
        Ok(Self {
            correlation,
            envelope,
            uncertainty,
            components,
        })
    }

    /// The retained correlation identity.
    #[must_use]
    pub fn correlation(&self) -> &CorrelationIdentity {
        &self.correlation
    }

    /// Evaluate this candidate.  Output authority is always Estimate-only.
    pub fn evaluate(&self, input: &ReducedAeroInput) -> Result<CandidateWrench, ReducedAeroError> {
        self.validate()?;
        input.validate()?;
        let geometry = input.geometry;
        let density = input.gas.density_kg_per_m3;
        let relative_velocity = finite_vec3(
            input
                .kinematics
                .linear_velocity_world_m_per_s
                .minus(input.gas.velocity_world_m_per_s),
            "relative_velocity_world_m_per_s",
        )?;
        let speed = checked_norm(relative_velocity, "relative_speed_m_per_s")?;
        let angular = input.kinematics.angular_velocity_world_rad_per_s;
        let angular_speed = checked_norm(angular, "angular_speed_rad_per_s")?;
        let normal = input.pose.unit_normal()?;
        let relative_roughness = finite_derived(
            input.roughness.height_m / geometry.radius_m,
            "relative_roughness",
        )?;
        self.envelope
            .relative_roughness
            .require("relative_roughness", relative_roughness)?;

        if density > 0.0 && speed > 0.0 && self.components.form_drag.is_some() {
            let re = finite_derived(
                density * speed * (2.0 * geometry.radius_m) / input.gas.dynamic_viscosity_pa_s,
                "translational_reynolds",
            )?;
            self.envelope
                .translational_reynolds
                .require("translational_reynolds", re)?;
        }
        if density > 0.0 && angular_speed > 0.0 && self.has_rotational_component() {
            let re = finite_derived(
                density * angular_speed * geometry.radius_m * geometry.radius_m
                    / input.gas.dynamic_viscosity_pa_s,
                "rotational_reynolds",
            )?;
            self.envelope
                .rotational_reynolds
                .require("rotational_reynolds", re)?;
            let mach = finite_derived(
                angular_speed * geometry.radius_m / input.gas.speed_of_sound_m_per_s,
                "tip_mach",
            )?;
            if mach > self.envelope.maximum_tip_mach {
                return Err(ReducedAeroError::OutsideCorrelationDomain {
                    quantity: "tip_mach",
                    value: mach,
                    minimum: 0.0,
                    maximum: self.envelope.maximum_tip_mach,
                });
            }
        }

        let face_area = finite_derived(geometry.face_area_m2(), "face_area_m2")?;
        let rim_wetted_area = finite_derived(geometry.rim_wetted_area_m2(), "rim_wetted_area_m2")?;
        let edge_on_rim_silhouette = finite_derived(
            geometry.edge_on_rim_silhouette_m2(),
            "edge_on_rim_silhouette_m2",
        )?;
        let roughness_factor = finite_derived(1.0 + relative_roughness, "roughness_factor")?;
        let mut force = Vec3::ZERO;
        let mut torque = Vec3::ZERO;
        let mut form_force = Vec3::ZERO;
        let mut skin_torque = Vec3::ZERO;
        let mut edge_torque = Vec3::ZERO;
        let mut orientation_torque = Vec3::ZERO;

        if let Some(term) = self.components.form_drag {
            if speed > 0.0 && density > 0.0 {
                let direction = finite_vec3(
                    relative_velocity.scaled(1.0 / speed),
                    "relative_velocity_direction",
                )?;
                let face_projection = finite_derived(
                    checked_dot(normal, direction, "face_projection")?.abs(),
                    "face_projection",
                )?;
                let rim_projection = finite_derived(
                    (1.0 - face_projection * face_projection).max(0.0).sqrt(),
                    "rim_projection",
                )?;
                let projected_area = finite_derived(
                    face_area * face_projection + edge_on_rim_silhouette * rim_projection,
                    "projected_area_m2",
                )?;
                let force_scale = finite_derived(
                    -0.5 * density * term.coefficient * roughness_factor * projected_area * speed,
                    "form_force_scale",
                )?;
                form_force =
                    finite_vec3(relative_velocity.scaled(force_scale), "form_force_world_n")?;
                force = finite_vec3(force.plus(form_force), "force_world_n")?;
            }
        }

        let spin = finite_vec3(
            normal.scaled(checked_dot(angular, normal, "spin_projection")?),
            "spin_angular_velocity",
        )?;
        let reorientation = finite_vec3(angular.minus(spin), "reorientation_angular_velocity")?;
        let torque_scale = finite_derived(
            0.5 * density * roughness_factor * geometry.radius_m.powi(3),
            "torque_scale",
        )?;
        if let Some(term) = self.components.rotational_skin_friction {
            let scale = finite_derived(
                -torque_scale
                    * face_area
                    * term.coefficient
                    * checked_norm(spin, "spin_angular_speed")?,
                "rotational_skin_torque_scale",
            )?;
            skin_torque = finite_vec3(spin.scaled(scale), "rotational_skin_torque_world_n_m")?;
            torque = finite_vec3(torque.plus(skin_torque), "torque_world_n_m")?;
        }
        if let Some(term) = self.components.edge_flow {
            let scale = finite_derived(
                -torque_scale
                    * rim_wetted_area
                    * term.coefficient
                    * checked_norm(spin, "spin_angular_speed")?,
                "edge_flow_torque_scale",
            )?;
            edge_torque = finite_vec3(spin.scaled(scale), "edge_flow_torque_world_n_m")?;
            torque = finite_vec3(torque.plus(edge_torque), "torque_world_n_m")?;
        }
        if let Some(term) = self.components.orientation_rate_damping {
            let scale = finite_derived(
                -torque_scale
                    * face_area
                    * term.coefficient
                    * checked_norm(reorientation, "reorientation_angular_speed")?,
                "orientation_rate_torque_scale",
            )?;
            orientation_torque = finite_vec3(
                reorientation.scaled(scale),
                "orientation_rate_torque_world_n_m",
            )?;
            torque = finite_vec3(torque.plus(orientation_torque), "torque_world_n_m")?;
        }

        let torque_power_w = checked_dot(torque, angular, "torque_power_w")?;
        let relative_power_w = checked_sum(
            checked_dot(force, relative_velocity, "relative_force_power_w")?,
            torque_power_w,
            "relative_power_w",
        )?;
        if relative_power_w > 0.0 {
            return Err(ReducedAeroError::PassivePowerViolation { relative_power_w });
        }
        let body_power_w = checked_sum(
            checked_dot(
                force,
                input.kinematics.linear_velocity_world_m_per_s,
                "body_force_power_w",
            )?,
            torque_power_w,
            "body_power_w",
        )?;
        let ambient_boundary_power_w = checked_dot(
            force,
            input.gas.velocity_world_m_per_s,
            "ambient_boundary_power_w",
        )?;
        let candidate = CandidateWrench {
            correlation: self.correlation.clone(),
            force_world_n: force,
            torque_world_n_m: torque,
            components: ComponentWrenches {
                form_force_world_n: form_force,
                rotational_skin_torque_world_n_m: skin_torque,
                edge_flow_torque_world_n_m: edge_torque,
                orientation_rate_torque_world_n_m: orientation_torque,
            },
            receipt: ReducedAeroReceipt {
                authority: EstimateAuthority::EstimateOnly,
                gas_source_id: input.gas.source_id.clone(),
                roughness_source_id: input.roughness.source_id.clone(),
                correlation_uncertainty_source_id: self.uncertainty.source_id.clone(),
                coefficient_relative_half_width: self.uncertainty.coefficient_relative_half_width,
                relative_power_w,
                dissipated_relative_power_w: finite_derived(
                    -relative_power_w,
                    "dissipated_relative_power_w",
                )?,
                body_power_w,
                ambient_boundary_power_w,
                moving_ambient_requires_total_energy_accounting: input.gas.velocity_world_m_per_s
                    != Vec3::ZERO,
            },
        };
        candidate.validate()?;
        Ok(candidate)
    }

    fn has_rotational_component(&self) -> bool {
        self.components.rotational_skin_friction.is_some()
            || self.components.edge_flow.is_some()
            || self.components.orientation_rate_damping.is_some()
    }

    fn validate(&self) -> Result<(), ReducedAeroError> {
        self.correlation.validate()?;
        self.envelope.validate()?;
        self.uncertainty.validate()?;
        self.components.validate()
    }
}

/// Validated evaluation inputs, all in one named world frame.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedAeroInput {
    /// Frame identity shared by all vectors in this request.
    pub world_frame_id: String,
    /// Exterior circular-disc geometry.
    pub geometry: DiscGeometry,
    /// Disc orientation.
    pub pose: DiscPose,
    /// Body motion.
    pub kinematics: BodyKinematics,
    /// Ambient gas properties and far-field velocity.
    pub gas: GasProperties,
    /// Exterior surface roughness.
    pub roughness: SurfaceRoughness,
}

impl ReducedAeroInput {
    fn validate(&self) -> Result<(), ReducedAeroError> {
        checked_identity(&self.world_frame_id, "world_frame_id")?;
        self.geometry.validate()?;
        self.pose.validate()?;
        self.kinematics.validate()?;
        self.roughness.validate()
    }
}

/// Force/torque contributions retained separately so a total coefficient cannot
/// conceal the physical family that supplied it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComponentWrenches {
    /// Form-drag force [N].
    pub form_force_world_n: Vec3,
    /// Face rotational-skin torque [N m].
    pub rotational_skin_torque_world_n_m: Vec3,
    /// Exterior-rim edge-flow torque [N m].
    pub edge_flow_torque_world_n_m: Vec3,
    /// Disc-axis-rate damping torque [N m].
    pub orientation_rate_torque_world_n_m: Vec3,
}

/// The authority carried by every reduced model result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimateAuthority {
    /// A correlation calculation, not independently qualified physical truth.
    EstimateOnly,
}

/// Explicit energy/power information for a single wrench evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedAeroReceipt {
    /// Results remain estimate-only.
    pub authority: EstimateAuthority,
    /// Retained gas-card identity.
    pub gas_source_id: String,
    /// Retained roughness-card identity.
    pub roughness_source_id: String,
    /// Retained source for caller-supplied coefficient uncertainty.
    pub correlation_uncertainty_source_id: String,
    /// Caller-supplied coefficient relative half-width; not a generated or
    /// independently validated physical confidence interval.
    pub coefficient_relative_half_width: f64,
    /// `F dot (v_body-v_ambient) + M dot omega` [W]; passive terms make this non-positive.
    pub relative_power_w: f64,
    /// `-relative_power_w` [W], non-negative for this passive model.
    pub dissipated_relative_power_w: f64,
    /// `F dot v_body + M dot omega` [W], power into the body.
    pub body_power_w: f64,
    /// `F dot v_ambient` [W]; with moving ambient, use this in a total accounting window.
    pub ambient_boundary_power_w: f64,
    /// True exactly when a far-field boundary supplies/removes mechanical power.
    pub moving_ambient_requires_total_energy_accounting: bool,
}

/// One correlation candidate's wrench, expressed in its request's world frame.
#[derive(Debug, Clone, PartialEq)]
pub struct CandidateWrench {
    /// Correlation identity.
    pub correlation: CorrelationIdentity,
    /// Total exterior force on the body [N].
    pub force_world_n: Vec3,
    /// Total exterior torque on the body about the request reference point [N m].
    pub torque_world_n_m: Vec3,
    /// Separately retained component wrenches.
    pub components: ComponentWrenches,
    /// Estimate authority and energy-accounting information.
    pub receipt: ReducedAeroReceipt,
}

impl CandidateWrench {
    fn validate(&self) -> Result<(), ReducedAeroError> {
        self.correlation.validate()?;
        for (value, field) in [
            (self.force_world_n, "candidate.force_world_n"),
            (self.torque_world_n_m, "candidate.torque_world_n_m"),
            (
                self.components.form_force_world_n,
                "candidate.components.form_force_world_n",
            ),
            (
                self.components.rotational_skin_torque_world_n_m,
                "candidate.components.rotational_skin_torque_world_n_m",
            ),
            (
                self.components.edge_flow_torque_world_n_m,
                "candidate.components.edge_flow_torque_world_n_m",
            ),
            (
                self.components.orientation_rate_torque_world_n_m,
                "candidate.components.orientation_rate_torque_world_n_m",
            ),
        ] {
            if !value.finite() {
                return Err(ReducedAeroError::InvalidCandidateWrench { field });
            }
        }
        if !within_vec_roundoff(self.force_world_n, self.components.form_force_world_n) {
            return Err(ReducedAeroError::InvalidCandidateWrench {
                field: "candidate.force_world_n",
            });
        }
        let component_torque = finite_vec3(
            self.components
                .rotational_skin_torque_world_n_m
                .plus(self.components.edge_flow_torque_world_n_m)
                .plus(self.components.orientation_rate_torque_world_n_m),
            "candidate.components.torque_world_n_m",
        )?;
        if !within_vec_roundoff(self.torque_world_n_m, component_torque) {
            return Err(ReducedAeroError::InvalidCandidateWrench {
                field: "candidate.torque_world_n_m",
            });
        }
        checked_identity(
            &self.receipt.gas_source_id,
            "candidate.receipt.gas_source_id",
        )?;
        checked_identity(
            &self.receipt.roughness_source_id,
            "candidate.receipt.roughness_source_id",
        )?;
        checked_identity(
            &self.receipt.correlation_uncertainty_source_id,
            "candidate.receipt.correlation_uncertainty_source_id",
        )?;
        finite_nonnegative(
            self.receipt.coefficient_relative_half_width,
            "candidate.receipt.coefficient_relative_half_width",
        )?;
        for (value, field) in [
            (
                self.receipt.relative_power_w,
                "candidate.receipt.relative_power_w",
            ),
            (
                self.receipt.dissipated_relative_power_w,
                "candidate.receipt.dissipated_relative_power_w",
            ),
            (self.receipt.body_power_w, "candidate.receipt.body_power_w"),
            (
                self.receipt.ambient_boundary_power_w,
                "candidate.receipt.ambient_boundary_power_w",
            ),
        ] {
            if !value.is_finite() {
                return Err(ReducedAeroError::InvalidCandidateWrench { field });
            }
        }
        if self.receipt.relative_power_w > 0.0 {
            return Err(ReducedAeroError::PassivePowerViolation {
                relative_power_w: self.receipt.relative_power_w,
            });
        }
        if self.receipt.dissipated_relative_power_w < 0.0
            || !within_roundoff(
                self.receipt.dissipated_relative_power_w,
                -self.receipt.relative_power_w,
            )
        {
            return Err(ReducedAeroError::InvalidCandidateWrench {
                field: "candidate.receipt.dissipated_relative_power_w",
            });
        }
        let reconstructed_body_power_w = checked_sum(
            self.receipt.relative_power_w,
            self.receipt.ambient_boundary_power_w,
            "candidate.receipt.body_power_w",
        )?;
        if !within_roundoff(self.receipt.body_power_w, reconstructed_body_power_w) {
            return Err(ReducedAeroError::InvalidCandidateWrench {
                field: "candidate.receipt.body_power_w",
            });
        }
        Ok(())
    }
}

fn within_roundoff(left: f64, right: f64) -> bool {
    if left == right {
        return true;
    }
    let scale = left.abs().max(right.abs()).max(f64::MIN_POSITIVE);
    (left - right).abs() <= 64.0 * f64::EPSILON * scale
}

fn within_vec_roundoff(left: Vec3, right: Vec3) -> bool {
    within_roundoff(left.x, right.x)
        && within_roundoff(left.y, right.y)
        && within_roundoff(left.z, right.z)
}

/// Deterministic set of alternative correlation candidates.  Candidates are
/// sorted by canonical identity so input iteration order cannot affect output.
#[derive(Debug, Clone, PartialEq)]
pub struct AlternativeWrenchSet {
    /// Candidate results in deterministic identity order.
    pub candidates: Vec<CandidateWrench>,
}

impl AlternativeWrenchSet {
    /// Evaluate independently retained alternatives without averaging away their
    /// disagreement.
    pub fn evaluate(
        models: &[ReducedAeroModel],
        input: &ReducedAeroInput,
    ) -> Result<Self, ReducedAeroError> {
        input.validate()?;
        if models.is_empty() {
            return Err(ReducedAeroError::EmptyAlternativeWrenchSet);
        }
        let mut identities = BTreeSet::new();
        for model in models {
            model.validate()?;
            if !identities.insert(model.correlation.clone()) {
                return Err(ReducedAeroError::DuplicateCorrelationIdentity {
                    correlation: model.correlation.clone(),
                });
            }
        }
        let mut candidates = models
            .iter()
            .map(|model| model.evaluate(input))
            .collect::<Result<Vec<_>, _>>()?;
        candidates.sort_by(|left, right| left.correlation.cmp(&right.correlation));
        Ok(Self { candidates })
    }

    /// Whether the retained alternatives disagree in their total wrench.
    #[must_use]
    pub fn has_force_or_torque_disagreement(&self) -> bool {
        self.candidates.windows(2).any(|pair| match pair {
            [left, right] => {
                left.force_world_n != right.force_world_n
                    || left.torque_world_n_m != right.torque_world_n_m
            }
            _ => false,
        })
    }
}

/// Caller-owned, exact-once work accumulation for one accounting window.
///
/// This small receipt validates public candidate arithmetic and internal
/// self-consistency before preventing duplicate application at this domain
/// seam. It does not validate correlation authority or replace `fs-couple`'s
/// closed-window audit; a coupling driver must still submit the selected
/// gas/boundary chart there.
#[derive(Debug, Default)]
pub struct WorkWindow {
    recorded_keys: BTreeSet<u64>,
    body_work_j: f64,
    relative_dissipation_j: f64,
}

impl WorkWindow {
    /// Record one candidate exactly once for a finite, non-negative duration.
    pub fn record_once(
        &mut self,
        exchange_key: u64,
        duration_s: f64,
        wrench: &CandidateWrench,
    ) -> Result<WorkReceipt, ReducedAeroError> {
        if !(duration_s.is_finite() && duration_s >= 0.0) {
            return Err(ReducedAeroError::InvalidWorkDuration);
        }
        wrench.validate()?;
        finite_derived(self.body_work_j, "work_window.body_work_j")?;
        finite_derived(
            self.relative_dissipation_j,
            "work_window.relative_dissipation_j",
        )?;
        let body_work_j =
            checked_product(wrench.receipt.body_power_w, duration_s, "work.body_work_j")?;
        let relative_dissipation_j = checked_product(
            wrench.receipt.dissipated_relative_power_w,
            duration_s,
            "work.relative_dissipation_j",
        )?;
        let next_body_work_j =
            checked_sum(self.body_work_j, body_work_j, "work_window.body_work_j")?;
        let next_relative_dissipation_j = checked_sum(
            self.relative_dissipation_j,
            relative_dissipation_j,
            "work_window.relative_dissipation_j",
        )?;
        if self.recorded_keys.contains(&exchange_key) {
            return Err(ReducedAeroError::DuplicateWorkExchange { key: exchange_key });
        }
        self.recorded_keys.insert(exchange_key);
        self.body_work_j = next_body_work_j;
        self.relative_dissipation_j = next_relative_dissipation_j;
        Ok(WorkReceipt {
            exchange_key,
            duration_s,
            body_work_j,
            relative_dissipation_j,
        })
    }

    /// Accumulated body work [J].
    #[must_use]
    pub fn body_work_j(&self) -> f64 {
        self.body_work_j
    }

    /// Accumulated passive relative dissipation [J].
    #[must_use]
    pub fn relative_dissipation_j(&self) -> f64 {
        self.relative_dissipation_j
    }
}

/// One exactly-once work receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorkReceipt {
    /// Caller-supplied exact-once key.
    pub exchange_key: u64,
    /// Integration duration [s].
    pub duration_s: f64,
    /// Wrench work into body [J].
    pub body_work_j: f64,
    /// Passive relative dissipation [J].
    pub relative_dissipation_j: f64,
}
