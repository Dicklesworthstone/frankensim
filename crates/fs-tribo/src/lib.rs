#![forbid(unsafe_code)]
//! Typed dry-contact baseline with explicit caller provenance and refusal boundaries.
//!
//! This crate never queries a material database or upgrades caller input into an admitted claim.

/// Solver-independent finite-patch tangential partial-slip rung.
pub mod partial_slip;
/// Solver-independent rolling and contour-deformation loss candidates.
pub mod rolling_loss;

use core::fmt;

const EPSILON: f64 = 64.0 * f64::EPSILON;

/// The declared ceiling of the caller's numerical input. This is not a receipt or admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputAuthority {
    /// The caller supplied a named source, but this crate has not independently admitted it.
    CallerDeclared,
    /// An explicitly synthetic fixture, suitable for tests and analytic exercises only.
    SyntheticFixture,
    /// An explicitly estimated value with no validation claim.
    Estimated,
}

/// The medium at the dry-law boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterfaceMedium {
    /// Explicitly dry contact.
    Dry,
    /// A fluid or lubrication film is present.
    Lubricated,
    /// The medium was not declared.
    Undeclared,
}

/// Provenance retained on every constitutive response.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputProvenance {
    source_id: String,
    authority: InputAuthority,
}

impl InputProvenance {
    /// Named caller source; this crate does not verify or upgrade it.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
    /// Caller-declared authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }
}

/// A checked ordered interface and history identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InterfaceSystemRef {
    ordered_system_id: String,
    history_id: String,
    provenance: InputProvenance,
    medium: InterfaceMedium,
}

impl InterfaceSystemRef {
    /// Creates an ordered dry-interface identity. Reversing A/B requires a new `ordered_system_id`.
    pub fn new(
        ordered_system_id: impl Into<String>,
        history_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        medium: InterfaceMedium,
    ) -> Result<Self, TriboError> {
        let value = Self {
            ordered_system_id: ordered_system_id.into(),
            history_id: history_id.into(),
            provenance: InputProvenance {
                source_id: source_id.into(),
                authority,
            },
            medium,
        };
        value.validate()?;
        Ok(value)
    }

    /// Ordered interface-system identity.
    #[must_use]
    pub fn ordered_system_id(&self) -> &str {
        &self.ordered_system_id
    }
    /// Named runtime history identity.
    #[must_use]
    pub fn history_id(&self) -> &str {
        &self.history_id
    }
    /// Caller source and authority ceiling.
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
    /// Declared medium.
    #[must_use]
    pub const fn medium(&self) -> InterfaceMedium {
        self.medium
    }

    fn validate(&self) -> Result<(), TriboError> {
        nonblank(&self.ordered_system_id, "ordered_system_id")?;
        nonblank(&self.history_id, "history_id")?;
        nonblank(&self.provenance.source_id, "source_id")?;
        if self.medium != InterfaceMedium::Dry {
            return Err(TriboError::NotDryInterface {
                medium: self.medium,
            });
        }
        Ok(())
    }
}

/// A checked contact frame that identifies the normal excluded from tangential friction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactFrame {
    normal: [f64; 3],
}

impl ContactFrame {
    /// Builds a frame from a finite non-zero normal and stores its unit normal.
    pub fn new(normal: [f64; 3]) -> Result<Self, TriboError> {
        finite_vec(normal, "contact_normal")?;
        let length = stable_norm(normal, "contact_normal")?;
        positive_finite(length, "contact_normal_length")?;
        let unit = scale_checked(normal, 1.0 / length, "contact_normal")?;
        Ok(Self { normal: unit })
    }

    /// Unit normal pointing out of the support surface.
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }
}

/// A checked relative velocity constrained to the tangent plane of a `ContactFrame`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentialSlip {
    velocity_mps: [f64; 3],
}

impl TangentialSlip {
    /// Refuses a velocity with a material normal component instead of silently discarding it.
    pub fn new(frame: &ContactFrame, velocity_mps: [f64; 3]) -> Result<Self, TriboError> {
        finite_vec(velocity_mps, "slip_velocity")?;
        let speed = stable_norm(velocity_mps, "slip_velocity")?;
        let normal_component = stable_dot(frame.normal, velocity_mps, "slip_normal_component")?;
        if normal_component.abs() > EPSILON * speed.max(1.0) {
            return Err(TriboError::NormalSlipComponent { normal_component });
        }
        Ok(Self { velocity_mps })
    }

    /// Tangential relative velocity in m/s.
    #[must_use]
    pub const fn velocity_mps(&self) -> [f64; 3] {
        self.velocity_mps
    }
    fn speed(self) -> Result<f64, TriboError> {
        stable_norm(self.velocity_mps, "slip_velocity")
    }
}

/// Total refusal surface for this leaf.
#[derive(Debug, Clone, PartialEq)]
pub enum TriboError {
    /// A required identity is empty.
    MissingIdentity { field: &'static str },
    /// Dry law received a lubrication or unknown-medium interface.
    NotDryInterface { medium: InterfaceMedium },
    /// A finite scalar or a derived candidate was physically invalid or unrepresentable.
    InvalidInput { field: &'static str },
    /// A vector contains a non-finite component.
    NonFiniteVector { field: &'static str },
    /// Submitted velocity leaks out of the declared tangent plane.
    NormalSlipComponent { normal_component: f64 },
    /// A heat fraction does not sum to unity.
    InvalidHeatPartition { sum: f64 },
    /// A forged or derived dissipation increment breaks an invariant.
    InvalidDissipationStep { field: &'static str },
    /// A closed applicability interval is malformed.
    InvalidApplicabilityRange { field: &'static str },
    /// A query value lies outside the named card interval.
    OutsideApplicability {
        field: &'static str,
        value: f64,
        minimum: f64,
        maximum: f64,
    },
}

impl fmt::Display for TriboError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingIdentity { field } => write!(f, "nonblank identity required: {field}"),
            Self::NotDryInterface { medium } => write!(f, "dry law refuses {medium:?} interface"),
            Self::InvalidInput { field } => {
                write!(f, "invalid or non-finite physical input: {field}")
            }
            Self::NonFiniteVector { field } => write!(f, "non-finite vector: {field}"),
            Self::NormalSlipComponent { normal_component } => {
                write!(f, "slip has normal component {normal_component}")
            }
            Self::InvalidHeatPartition { sum } => {
                write!(f, "heat partition must sum to one, got {sum}")
            }
            Self::InvalidDissipationStep { field } => {
                write!(f, "invalid dissipation step: {field}")
            }
            Self::InvalidApplicabilityRange { field } => {
                write!(f, "invalid applicability range: {field}")
            }
            Self::OutsideApplicability {
                field,
                value,
                minimum,
                maximum,
            } => write!(
                f,
                "{field}={value} lies outside declared applicability [{minimum}, {maximum}]"
            ),
        }
    }
}
impl std::error::Error for TriboError {}

/// Dry-friction ladder. Parameters are caller-supplied and retain only caller authority.
#[derive(Debug, Clone, PartialEq)]
pub enum FrictionLaw {
    /// Coulomb static and kinetic coefficients.
    Coulomb { static_mu: f64, kinetic_mu: f64 },
    /// `mu_zero + slope_per_speed * |v|`, lower-clamped at zero.
    VelocityDependent {
        static_mu: f64,
        mu_zero: f64,
        slope_per_speed: f64,
    },
    /// Exponential Stribeck decay plus a non-negative viscous coefficient.
    Stribeck {
        static_mu: f64,
        kinetic_mu: f64,
        characteristic_speed: f64,
        viscous_per_speed: f64,
    },
}

impl FrictionLaw {
    /// Reports static capacity at rest; during sliding traction opposes the checked tangent slip.
    pub fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        normal_force: f64,
        slip: TangentialSlip,
    ) -> Result<FrictionResponse, TriboError> {
        interface.validate()?;
        nonnegative_finite(normal_force, "normal_force")?;
        self.validate()?;
        let static_limit = checked_mul(self.static_mu(), normal_force, "static_limit")?;
        let speed = slip.speed()?;
        if speed == 0.0 {
            return Ok(FrictionResponse {
                regime: FrictionRegime::Sticking,
                static_limit,
                kinetic_coefficient: None,
                traction_n: [0.0; 3],
                dissipated_power_w: 0.0,
                provenance: interface.provenance.clone(),
            });
        }
        let coefficient = self.kinetic_mu(speed)?;
        let magnitude = checked_mul(coefficient, normal_force, "friction_magnitude")?;
        let traction_n = scale_checked(slip.velocity_mps, -magnitude / speed, "traction")?;
        let dissipated_power_w = -stable_dot(traction_n, slip.velocity_mps, "friction_power")?;
        nonnegative_finite(dissipated_power_w, "friction_power")?;
        Ok(FrictionResponse {
            regime: FrictionRegime::Sliding,
            static_limit,
            kinetic_coefficient: Some(coefficient),
            traction_n,
            dissipated_power_w,
            provenance: interface.provenance.clone(),
        })
    }

    fn validate(&self) -> Result<(), TriboError> {
        let (static_mu, kinetic_mu) = match *self {
            Self::Coulomb {
                static_mu,
                kinetic_mu,
            }
            | Self::Stribeck {
                static_mu,
                kinetic_mu,
                ..
            } => (static_mu, kinetic_mu),
            Self::VelocityDependent {
                static_mu,
                mu_zero,
                slope_per_speed,
            } => {
                nonnegative_finite(slope_per_speed.abs(), "slope_per_speed")?;
                if !slope_per_speed.is_finite() {
                    return Err(TriboError::InvalidInput {
                        field: "slope_per_speed",
                    });
                }
                (static_mu, mu_zero)
            }
        };
        nonnegative_finite(static_mu, "static_mu")?;
        nonnegative_finite(kinetic_mu, "kinetic_mu")?;
        if kinetic_mu > static_mu {
            return Err(TriboError::InvalidInput {
                field: "kinetic_mu > static_mu",
            });
        }
        if let Self::Stribeck {
            characteristic_speed,
            viscous_per_speed,
            ..
        } = *self
        {
            positive_finite(characteristic_speed, "characteristic_speed")?;
            nonnegative_finite(viscous_per_speed, "viscous_per_speed")?;
        }
        Ok(())
    }

    fn static_mu(&self) -> f64 {
        match *self {
            Self::Coulomb { static_mu, .. }
            | Self::VelocityDependent { static_mu, .. }
            | Self::Stribeck { static_mu, .. } => static_mu,
        }
    }
    fn kinetic_mu(&self, speed: f64) -> Result<f64, TriboError> {
        let candidate = match *self {
            Self::Coulomb { kinetic_mu, .. } => kinetic_mu,
            Self::VelocityDependent {
                mu_zero,
                slope_per_speed,
                ..
            } => checked_add(
                mu_zero,
                checked_mul(slope_per_speed, speed, "velocity_mu")?,
                "velocity_mu",
            )?
            .max(0.0),
            Self::Stribeck {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                viscous_per_speed,
            } => {
                let ratio = speed / characteristic_speed;
                finite(ratio, "stribeck_speed_ratio")?;
                let decay = (-(ratio * ratio)).exp();
                finite(decay, "stribeck_decay")?;
                checked_add(
                    checked_add(
                        kinetic_mu,
                        checked_mul(static_mu - kinetic_mu, decay, "stribeck_mu")?,
                        "stribeck_mu",
                    )?,
                    checked_mul(viscous_per_speed, speed, "stribeck_mu")?,
                    "stribeck_mu",
                )?
            }
        };
        nonnegative_finite(candidate, "kinetic_coefficient")?;
        Ok(candidate)
    }
}

/// Contact solver branch represented by a friction response.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrictionRegime {
    Sticking,
    Sliding,
}

/// One friction query result; all force values are newtons and power is watts.
#[derive(Debug, Clone, PartialEq)]
pub struct FrictionResponse {
    pub regime: FrictionRegime,
    pub static_limit: f64,
    pub kinetic_coefficient: Option<f64>,
    traction_n: [f64; 3],
    dissipated_power_w: f64,
    provenance: InputProvenance,
}
impl FrictionResponse {
    #[must_use]
    pub const fn traction_n(&self) -> [f64; 3] {
        self.traction_n
    }
    #[must_use]
    pub const fn dissipated_power_w(&self) -> f64 {
        self.dissipated_power_w
    }
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
}

/// Closed non-negative applicability interval in the SI unit of its owning field.
///
/// The range itself carries no material admission or uncertainty authority.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApplicabilityRange {
    minimum: f64,
    maximum: f64,
}

impl ApplicabilityRange {
    /// Builds a finite closed interval with `0 <= minimum <= maximum`.
    pub fn new(minimum: f64, maximum: f64) -> Result<Self, TriboError> {
        nonnegative_finite(minimum, "applicability.minimum")?;
        nonnegative_finite(maximum, "applicability.maximum")?;
        if minimum > maximum {
            return Err(TriboError::InvalidApplicabilityRange {
                field: "minimum > maximum",
            });
        }
        Ok(Self { minimum, maximum })
    }

    /// Inclusive lower bound.
    #[must_use]
    pub const fn minimum(&self) -> f64 {
        self.minimum
    }

    /// Inclusive upper bound.
    #[must_use]
    pub const fn maximum(&self) -> f64 {
        self.maximum
    }

    fn require_contains(self, value: f64, field: &'static str) -> Result<(), TriboError> {
        nonnegative_finite(value, field)?;
        if value < self.minimum || value > self.maximum {
            return Err(TriboError::OutsideApplicability {
                field,
                value,
                minimum: self.minimum,
                maximum: self.maximum,
            });
        }
        Ok(())
    }
}

/// Validity domain declared by a dry-friction system card.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryFrictionApplicability {
    temperature_kelvin: ApplicabilityRange,
    contact_pressure_pa: ApplicabilityRange,
    slip_speed_mps: ApplicabilityRange,
}

impl DryFrictionApplicability {
    /// Builds an ordered dry-law validity domain. Temperature is absolute and
    /// therefore has a strictly positive lower bound.
    pub fn new(
        temperature_kelvin: ApplicabilityRange,
        contact_pressure_pa: ApplicabilityRange,
        slip_speed_mps: ApplicabilityRange,
    ) -> Result<Self, TriboError> {
        if temperature_kelvin.minimum == 0.0 {
            return Err(TriboError::InvalidApplicabilityRange {
                field: "temperature_kelvin.minimum",
            });
        }
        Ok(Self {
            temperature_kelvin,
            contact_pressure_pa,
            slip_speed_mps,
        })
    }

    /// Declared absolute-temperature domain [K].
    #[must_use]
    pub const fn temperature_kelvin(&self) -> ApplicabilityRange {
        self.temperature_kelvin
    }

    /// Declared nominal-contact-pressure domain [Pa].
    #[must_use]
    pub const fn contact_pressure_pa(&self) -> ApplicabilityRange {
        self.contact_pressure_pa
    }

    /// Declared tangential slip-speed domain [m/s].
    #[must_use]
    pub const fn slip_speed_mps(&self) -> ApplicabilityRange {
        self.slip_speed_mps
    }

    fn require_contains(self, state: DryFrictionState) -> Result<(), TriboError> {
        positive_finite(state.temperature_kelvin, "temperature_kelvin")?;
        self.temperature_kelvin
            .require_contains(state.temperature_kelvin, "temperature_kelvin")?;
        self.contact_pressure_pa
            .require_contains(state.contact_pressure_pa, "contact_pressure_pa")?;
        self.slip_speed_mps
            .require_contains(state.slip_speed_mps, "slip_speed_mps")
    }
}

/// Runtime state checked against a `DryFrictionApplicability` domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryFrictionState {
    temperature_kelvin: f64,
    contact_pressure_pa: f64,
    slip_speed_mps: f64,
}

impl DryFrictionState {
    /// Creates a finite SI state. Domain admission happens through the card.
    pub fn new(
        temperature_kelvin: f64,
        contact_pressure_pa: f64,
        slip_speed_mps: f64,
    ) -> Result<Self, TriboError> {
        positive_finite(temperature_kelvin, "temperature_kelvin")?;
        nonnegative_finite(contact_pressure_pa, "contact_pressure_pa")?;
        nonnegative_finite(slip_speed_mps, "slip_speed_mps")?;
        Ok(Self {
            temperature_kelvin,
            contact_pressure_pa,
            slip_speed_mps,
        })
    }

    /// Absolute temperature [K].
    #[must_use]
    pub const fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }

    /// Nominal contact pressure [Pa].
    #[must_use]
    pub const fn contact_pressure_pa(&self) -> f64 {
        self.contact_pressure_pa
    }

    /// Tangential slip-speed magnitude [m/s].
    #[must_use]
    pub const fn slip_speed_mps(&self) -> f64 {
        self.slip_speed_mps
    }
}

/// Dependency-independent, ordered dry-interface card.
///
/// This is deliberately card-*style*, rather than an alias of the pending
/// `fs-matdb` card: no seed-data admission is available in this leaf. It binds
/// the existing ordered surface/history identity, caller authority, law, and
/// validity domain so consumers cannot silently use a scalar coefficient beyond
/// the declared state regime.
#[derive(Debug, Clone, PartialEq)]
pub struct DryInterfaceSystemCard {
    interface: InterfaceSystemRef,
    friction_law: FrictionLaw,
    applicability: DryFrictionApplicability,
}

impl DryInterfaceSystemCard {
    /// Creates an immutable ordered interface card without upgrading caller authority.
    pub fn new(
        interface: InterfaceSystemRef,
        friction_law: FrictionLaw,
        applicability: DryFrictionApplicability,
    ) -> Result<Self, TriboError> {
        interface.validate()?;
        friction_law.validate()?;
        Ok(Self {
            interface,
            friction_law,
            applicability,
        })
    }

    /// The ordered surface/history identity retained by the card.
    #[must_use]
    pub const fn interface(&self) -> &InterfaceSystemRef {
        &self.interface
    }

    /// The caller-declared friction rung.
    #[must_use]
    pub const fn friction_law(&self) -> &FrictionLaw {
        &self.friction_law
    }

    /// The closed state domain required for queries.
    #[must_use]
    pub const fn applicability(&self) -> DryFrictionApplicability {
        self.applicability
    }

    /// Evaluates only after state-domain and tangential-speed consistency checks.
    pub fn query(
        &self,
        state: DryFrictionState,
        normal_force_n: f64,
        slip: TangentialSlip,
    ) -> Result<FrictionQueryReceipt, TriboError> {
        self.applicability.require_contains(state)?;
        let slip_speed_mps = slip.speed()?;
        let mismatch = (slip_speed_mps - state.slip_speed_mps).abs();
        if mismatch > EPSILON * slip_speed_mps.max(state.slip_speed_mps).max(1.0) {
            return Err(TriboError::InvalidInput {
                field: "state.slip_speed_mps",
            });
        }
        let response = self
            .friction_law
            .evaluate(&self.interface, normal_force_n, slip)?;
        Ok(FrictionQueryReceipt { state, response })
    }
}

/// A card-checked friction response with the exact admitted runtime state.
#[derive(Debug, Clone, PartialEq)]
pub struct FrictionQueryReceipt {
    state: DryFrictionState,
    response: FrictionResponse,
}

impl FrictionQueryReceipt {
    /// State checked against the immutable card domain.
    #[must_use]
    pub const fn state(&self) -> DryFrictionState {
        self.state
    }

    /// Friction response retaining the card interface provenance.
    #[must_use]
    pub const fn response(&self) -> &FrictionResponse {
        &self.response
    }
}

/// Hertz sphere/plane geometry and modulus, in metres and pascals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzSpherePlane {
    pub effective_radius: f64,
    pub reduced_modulus: f64,
}
/// Hertz sphere/plane response: newtons, metres, pascals respectively.
#[derive(Debug, Clone, PartialEq)]
pub struct HertzSphereResponse {
    pub normal_force_n: f64,
    pub contact_radius_m: f64,
    pub peak_pressure_pa: f64,
    provenance: InputProvenance,
}
impl HertzSphereResponse {
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
}

impl HertzSpherePlane {
    /// Elastic closed-form response for non-negative indentation in metres.
    pub fn response(
        self,
        interface: &InterfaceSystemRef,
        indentation_m: f64,
    ) -> Result<HertzSphereResponse, TriboError> {
        interface.validate()?;
        positive_finite(self.effective_radius, "effective_radius")?;
        positive_finite(self.reduced_modulus, "reduced_modulus")?;
        nonnegative_finite(indentation_m, "indentation_m")?;
        let contact_radius_m = checked_sqrt(
            checked_mul(self.effective_radius, indentation_m, "contact_radius")?,
            "contact_radius",
        )?;
        let normal_force_n = checked_mul(
            4.0 / 3.0,
            checked_mul(
                self.reduced_modulus,
                checked_mul(
                    self.effective_radius.sqrt(),
                    indentation_m.powf(1.5),
                    "normal_force",
                )?,
                "normal_force",
            )?,
            "normal_force",
        )?;
        let peak_pressure_pa = if contact_radius_m == 0.0 {
            0.0
        } else {
            let area = checked_mul(
                core::f64::consts::PI,
                checked_mul(contact_radius_m, contact_radius_m, "contact_area")?,
                "contact_area",
            )?;
            checked_mul(1.5, normal_force_n / area, "peak_pressure")?
        };
        nonnegative_finite(peak_pressure_pa, "peak_pressure")?;
        Ok(HertzSphereResponse {
            normal_force_n,
            contact_radius_m,
            peak_pressure_pa,
            provenance: interface.provenance.clone(),
        })
    }
}

/// Hertz cylinder/plane geometry and reduced modulus, in metres and pascals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzCylinderPlane {
    pub effective_radius: f64,
    pub reduced_modulus: f64,
}
/// Hertz cylinder/plane response: half-width in metres and pressure in pascals.
#[derive(Debug, Clone, PartialEq)]
pub struct HertzCylinderResponse {
    pub half_width_m: f64,
    pub peak_pressure_pa: f64,
    provenance: InputProvenance,
}
impl HertzCylinderResponse {
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
}

impl HertzCylinderPlane {
    /// Elastic line-contact response for normal line load in N/m.
    pub fn response(
        self,
        interface: &InterfaceSystemRef,
        normal_line_load_n_per_m: f64,
    ) -> Result<HertzCylinderResponse, TriboError> {
        interface.validate()?;
        positive_finite(self.effective_radius, "effective_radius")?;
        positive_finite(self.reduced_modulus, "reduced_modulus")?;
        nonnegative_finite(normal_line_load_n_per_m, "normal_line_load_n_per_m")?;
        let numerator = checked_mul(
            4.0,
            checked_mul(
                normal_line_load_n_per_m,
                self.effective_radius,
                "half_width",
            )?,
            "half_width",
        )?;
        let denominator = checked_mul(core::f64::consts::PI, self.reduced_modulus, "half_width")?;
        let half_width_m = checked_sqrt(numerator / denominator, "half_width")?;
        let peak_pressure_pa = if half_width_m == 0.0 {
            0.0
        } else {
            checked_mul(
                2.0,
                normal_line_load_n_per_m
                    / checked_mul(core::f64::consts::PI, half_width_m, "peak_pressure")?,
                "peak_pressure",
            )?
        };
        nonnegative_finite(peak_pressure_pa, "peak_pressure")?;
        Ok(HertzCylinderResponse {
            half_width_m,
            peak_pressure_pa,
            provenance: interface.provenance.clone(),
        })
    }
}

/// Generic rolling/contour loss protocol.
pub trait ResistanceLaw {
    /// Produces efforts opposing their matching signed rates and non-negative dissipation.
    fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        input: ResistanceInput,
    ) -> Result<ResistanceResponse, TriboError>;
}

/// Inputs shared by the supplied scalar resistance laws.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResistanceInput {
    pub normal_force_n: f64,
    pub angular_speed_rad_s: f64,
    pub contour_speed_mps: f64,
}
impl ResistanceInput {
    fn validate(self) -> Result<(), TriboError> {
        nonnegative_finite(self.normal_force_n, "normal_force_n")?;
        finite(self.angular_speed_rad_s, "angular_speed_rad_s")?;
        finite(self.contour_speed_mps, "contour_speed_mps")
    }
}
/// Response from one resistance mechanism.
#[derive(Debug, Clone, PartialEq)]
pub struct ResistanceResponse {
    pub rolling_moment_n_m: f64,
    pub contour_force_n: f64,
    dissipated_power_w: f64,
    provenance: InputProvenance,
}
impl ResistanceResponse {
    #[must_use]
    pub const fn dissipated_power_w(&self) -> f64 {
        self.dissipated_power_w
    }
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
}

/// A caller-supplied rolling moment arm in metres; it carries no calibration claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantRollingMoment {
    pub moment_arm_m: f64,
}
impl ResistanceLaw for ConstantRollingMoment {
    fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        input: ResistanceInput,
    ) -> Result<ResistanceResponse, TriboError> {
        interface.validate()?;
        input.validate()?;
        nonnegative_finite(self.moment_arm_m, "moment_arm_m")?;
        let rolling_moment_n_m = if input.angular_speed_rad_s == 0.0 {
            0.0
        } else {
            -input.angular_speed_rad_s.signum()
                * checked_mul(input.normal_force_n, self.moment_arm_m, "rolling_moment")?
        };
        let dissipated_power_w = -checked_mul(
            rolling_moment_n_m,
            input.angular_speed_rad_s,
            "rolling_power",
        )?;
        nonnegative_finite(dissipated_power_w, "rolling_power")?;
        Ok(ResistanceResponse {
            rolling_moment_n_m,
            contour_force_n: 0.0,
            dissipated_power_w,
            provenance: interface.provenance.clone(),
        })
    }
}

/// A caller-supplied contour force in newtons; it carries no geometry/calibration claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantContourForce {
    pub force_n: f64,
}
impl ResistanceLaw for ConstantContourForce {
    fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        input: ResistanceInput,
    ) -> Result<ResistanceResponse, TriboError> {
        interface.validate()?;
        input.validate()?;
        nonnegative_finite(self.force_n, "force_n")?;
        let contour_force_n = if input.contour_speed_mps == 0.0 {
            0.0
        } else {
            -input.contour_speed_mps.signum() * self.force_n
        };
        let dissipated_power_w =
            -checked_mul(contour_force_n, input.contour_speed_mps, "contour_power")?;
        nonnegative_finite(dissipated_power_w, "contour_power")?;
        Ok(ResistanceResponse {
            rolling_moment_n_m: 0.0,
            contour_force_n,
            dissipated_power_w,
            provenance: interface.provenance.clone(),
        })
    }
}

/// Validated fractional allocation of frictional work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatPartition {
    surface_a: f64,
    surface_b: f64,
    other: f64,
}
impl HeatPartition {
    pub fn new(surface_a: f64, surface_b: f64, other: f64) -> Result<Self, TriboError> {
        let value = Self {
            surface_a,
            surface_b,
            other,
        };
        value.validate()?;
        Ok(value)
    }
    #[must_use]
    pub const fn surface_a(&self) -> f64 {
        self.surface_a
    }
    #[must_use]
    pub const fn surface_b(&self) -> f64 {
        self.surface_b
    }
    #[must_use]
    pub const fn other(&self) -> f64 {
        self.other
    }
    fn validate(&self) -> Result<(), TriboError> {
        nonnegative_finite(self.surface_a, "surface_a")?;
        nonnegative_finite(self.surface_b, "surface_b")?;
        nonnegative_finite(self.other, "other")?;
        let sum = checked_add(
            checked_add(self.surface_a, self.surface_b, "heat_partition")?,
            self.other,
            "heat_partition",
        )?;
        if (sum - 1.0).abs() > EPSILON {
            return Err(TriboError::InvalidHeatPartition { sum });
        }
        Ok(())
    }
}

/// Thermal properties required by the bounded semi-infinite flash-temperature
/// candidate. They are caller data and carry no material-card admission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceThermalProperties {
    conductivity_w_per_m_k: f64,
    diffusivity_m2_per_s: f64,
}

impl SurfaceThermalProperties {
    /// Creates a finite, positive thermal-property pair in SI units.
    pub fn new(conductivity_w_per_m_k: f64, diffusivity_m2_per_s: f64) -> Result<Self, TriboError> {
        positive_finite(conductivity_w_per_m_k, "conductivity_w_per_m_k")?;
        positive_finite(diffusivity_m2_per_s, "diffusivity_m2_per_s")?;
        Ok(Self {
            conductivity_w_per_m_k,
            diffusivity_m2_per_s,
        })
    }

    /// Thermal conductivity [W/(m K)].
    #[must_use]
    pub const fn conductivity_w_per_m_k(&self) -> f64 {
        self.conductivity_w_per_m_k
    }

    /// Thermal diffusivity [m²/s].
    #[must_use]
    pub const fn diffusivity_m2_per_s(&self) -> f64 {
        self.diffusivity_m2_per_s
    }
}

/// Inputs for a bounded, uniform-flux, semi-infinite-body flash-temperature
/// candidate. `surface_*` may be absent, which yields a typed `Unknown` rather
/// than inventing a thermal material property.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlashTemperatureInput {
    dissipated_power_w: f64,
    partition: HeatPartition,
    contact_area_m2: f64,
    traverse_length_m: f64,
    slip_speed_mps: f64,
    surface_a: Option<SurfaceThermalProperties>,
    surface_b: Option<SurfaceThermalProperties>,
}

impl FlashTemperatureInput {
    /// Builds a flash candidate request. The dwell time is
    /// `traverse_length_m / slip_speed_mps`.
    pub fn new(
        dissipated_power_w: f64,
        partition: HeatPartition,
        contact_area_m2: f64,
        traverse_length_m: f64,
        slip_speed_mps: f64,
        surface_a: Option<SurfaceThermalProperties>,
        surface_b: Option<SurfaceThermalProperties>,
    ) -> Result<Self, TriboError> {
        nonnegative_finite(dissipated_power_w, "dissipated_power_w")?;
        partition.validate()?;
        positive_finite(contact_area_m2, "contact_area_m2")?;
        positive_finite(traverse_length_m, "traverse_length_m")?;
        positive_finite(slip_speed_mps, "slip_speed_mps")?;
        Ok(Self {
            dissipated_power_w,
            partition,
            contact_area_m2,
            traverse_length_m,
            slip_speed_mps,
            surface_a,
            surface_b,
        })
    }
}

/// Reason a flash-temperature request has no physical-property closure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlashTemperatureUnknown {
    /// Nonzero allocated heat reaches surface A without its thermal properties.
    MissingSurfaceAThermalProperties,
    /// Nonzero allocated heat reaches surface B without its thermal properties.
    MissingSurfaceBThermalProperties,
}

/// Explicit model-form result of a flash-temperature request.
#[derive(Debug, Clone, PartialEq)]
pub enum FlashTemperatureEstimate {
    /// A bounded-model candidate, not an error enclosure or temperature-port solution.
    Candidate(FlashTemperatureCandidate),
    /// Insufficient thermal data; no candidate was fabricated.
    Unknown(FlashTemperatureUnknown),
}

/// Uniform-flux, semi-infinite-body temperature-rise candidate.
///
/// For each receiving surface the reported rise is
/// `2 q'' sqrt(alpha * t / pi) / k`, where `q''` is that surface's explicit
/// partition of frictional power divided by contact area and
/// `t = traverse_length / slip_speed`. This is a model-form estimate under the
/// stated assumptions, never a flash-temperature bound or a thermal solve.
#[derive(Debug, Clone, PartialEq)]
pub struct FlashTemperatureCandidate {
    surface_a_rise_k: f64,
    surface_b_rise_k: f64,
    dwell_time_s: f64,
    surface_a_heat_flux_w_per_m2: f64,
    surface_b_heat_flux_w_per_m2: f64,
    provenance: InputProvenance,
}

impl FlashTemperatureCandidate {
    /// Surface-A temperature-rise candidate [K].
    #[must_use]
    pub const fn surface_a_rise_k(&self) -> f64 {
        self.surface_a_rise_k
    }

    /// Surface-B temperature-rise candidate [K].
    #[must_use]
    pub const fn surface_b_rise_k(&self) -> f64 {
        self.surface_b_rise_k
    }

    /// Uniform heat-flux dwell time [s].
    #[must_use]
    pub const fn dwell_time_s(&self) -> f64 {
        self.dwell_time_s
    }

    /// Explicit surface-A flux share [W/m²].
    #[must_use]
    pub const fn surface_a_heat_flux_w_per_m2(&self) -> f64 {
        self.surface_a_heat_flux_w_per_m2
    }

    /// Explicit surface-B flux share [W/m²].
    #[must_use]
    pub const fn surface_b_heat_flux_w_per_m2(&self) -> f64 {
        self.surface_b_heat_flux_w_per_m2
    }

    /// Caller provenance from the ordered dry interface.
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
}

/// Computes a flash-temperature candidate or reports exactly which receiving
/// surface lacks data. The ordered dry interface is validated first.
pub fn flash_temperature_candidate(
    interface: &InterfaceSystemRef,
    input: FlashTemperatureInput,
) -> Result<FlashTemperatureEstimate, TriboError> {
    interface.validate()?;
    let dwell_time_s = checked_div(
        input.traverse_length_m,
        input.slip_speed_mps,
        "flash_dwell_time_s",
    )?;
    let surface_a_heat_flux_w_per_m2 = checked_div(
        checked_mul(
            input.dissipated_power_w,
            input.partition.surface_a,
            "surface_a_heat_power",
        )?,
        input.contact_area_m2,
        "surface_a_heat_flux",
    )?;
    let surface_b_heat_flux_w_per_m2 = checked_div(
        checked_mul(
            input.dissipated_power_w,
            input.partition.surface_b,
            "surface_b_heat_power",
        )?,
        input.contact_area_m2,
        "surface_b_heat_flux",
    )?;
    let surface_a = match (surface_a_heat_flux_w_per_m2 == 0.0, input.surface_a) {
        (true, _) => 0.0,
        (false, Some(properties)) => flash_rise_k(
            surface_a_heat_flux_w_per_m2,
            dwell_time_s,
            properties,
            "surface_a_flash_temperature",
        )?,
        (false, None) => {
            return Ok(FlashTemperatureEstimate::Unknown(
                FlashTemperatureUnknown::MissingSurfaceAThermalProperties,
            ));
        }
    };
    let surface_b = match (surface_b_heat_flux_w_per_m2 == 0.0, input.surface_b) {
        (true, _) => 0.0,
        (false, Some(properties)) => flash_rise_k(
            surface_b_heat_flux_w_per_m2,
            dwell_time_s,
            properties,
            "surface_b_flash_temperature",
        )?,
        (false, None) => {
            return Ok(FlashTemperatureEstimate::Unknown(
                FlashTemperatureUnknown::MissingSurfaceBThermalProperties,
            ));
        }
    };
    Ok(FlashTemperatureEstimate::Candidate(
        FlashTemperatureCandidate {
            surface_a_rise_k: surface_a,
            surface_b_rise_k: surface_b,
            dwell_time_s,
            surface_a_heat_flux_w_per_m2,
            surface_b_heat_flux_w_per_m2,
            provenance: interface.provenance.clone(),
        },
    ))
}

fn flash_rise_k(
    heat_flux_w_per_m2: f64,
    dwell_time_s: f64,
    properties: SurfaceThermalProperties,
    field: &'static str,
) -> Result<f64, TriboError> {
    let diffusion_length_m = checked_sqrt(
        checked_div(
            checked_mul(properties.diffusivity_m2_per_s, dwell_time_s, field)?,
            core::f64::consts::PI,
            field,
        )?,
        field,
    )?;
    let rise_k = checked_div(
        checked_mul(
            2.0,
            checked_mul(heat_flux_w_per_m2, diffusion_length_m, field)?,
            field,
        )?,
        properties.conductivity_w_per_m_k,
        field,
    )?;
    nonnegative_finite(rise_k, field)?;
    Ok(rise_k)
}

/// A closed, validated dissipative-work increment in joules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DissipationStep {
    total_work_j: f64,
    surface_a_heat_j: f64,
    surface_b_heat_j: f64,
    other_work_j: f64,
}
impl DissipationStep {
    /// Multiplies a non-negative power (W) and duration (s), then applies a checked partition.
    pub fn from_power(
        dissipated_power_w: f64,
        duration_s: f64,
        partition: HeatPartition,
    ) -> Result<Self, TriboError> {
        nonnegative_finite(dissipated_power_w, "dissipated_power_w")?;
        nonnegative_finite(duration_s, "duration_s")?;
        partition.validate()?;
        let total_work_j = checked_mul(dissipated_power_w, duration_s, "total_work_j")?;
        let value = Self {
            total_work_j,
            surface_a_heat_j: checked_mul(total_work_j, partition.surface_a, "surface_a_heat_j")?,
            surface_b_heat_j: checked_mul(total_work_j, partition.surface_b, "surface_b_heat_j")?,
            other_work_j: checked_mul(total_work_j, partition.other, "other_work_j")?,
        };
        value.validate()?;
        Ok(value)
    }
    #[must_use]
    pub const fn total_work_j(&self) -> f64 {
        self.total_work_j
    }
    #[must_use]
    pub const fn surface_a_heat_j(&self) -> f64 {
        self.surface_a_heat_j
    }
    #[must_use]
    pub const fn surface_b_heat_j(&self) -> f64 {
        self.surface_b_heat_j
    }
    #[must_use]
    pub const fn other_work_j(&self) -> f64 {
        self.other_work_j
    }
    fn validate(&self) -> Result<(), TriboError> {
        nonnegative_finite(self.total_work_j, "total_work_j").map_err(|_| {
            TriboError::InvalidDissipationStep {
                field: "total_work_j",
            }
        })?;
        nonnegative_finite(self.surface_a_heat_j, "surface_a_heat_j").map_err(|_| {
            TriboError::InvalidDissipationStep {
                field: "surface_a_heat_j",
            }
        })?;
        nonnegative_finite(self.surface_b_heat_j, "surface_b_heat_j").map_err(|_| {
            TriboError::InvalidDissipationStep {
                field: "surface_b_heat_j",
            }
        })?;
        nonnegative_finite(self.other_work_j, "other_work_j").map_err(|_| {
            TriboError::InvalidDissipationStep {
                field: "other_work_j",
            }
        })?;
        let channels = checked_add(
            checked_add(
                self.surface_a_heat_j,
                self.surface_b_heat_j,
                "dissipation_channels",
            )
            .map_err(|_| TriboError::InvalidDissipationStep { field: "channels" })?,
            self.other_work_j,
            "dissipation_channels",
        )
        .map_err(|_| TriboError::InvalidDissipationStep { field: "channels" })?;
        if (channels - self.total_work_j).abs() > EPSILON * self.total_work_j.max(1.0) {
            return Err(TriboError::InvalidDissipationStep { field: "closure" });
        }
        Ok(())
    }
}

/// Caller-owned, checked cumulative work ledger.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WorkLedger {
    dissipated_work_j: f64,
    surface_a_heat_j: f64,
    surface_b_heat_j: f64,
    other_work_j: f64,
}
impl WorkLedger {
    #[must_use]
    pub const fn dissipated_work_j(&self) -> f64 {
        self.dissipated_work_j
    }
    /// Cumulative heat assigned to surface A [J].
    #[must_use]
    pub const fn surface_a_heat_j(&self) -> f64 {
        self.surface_a_heat_j
    }
    /// Cumulative heat assigned to surface B [J].
    #[must_use]
    pub const fn surface_b_heat_j(&self) -> f64 {
        self.surface_b_heat_j
    }
    /// Cumulative declared non-surface work share [J].
    #[must_use]
    pub const fn other_work_j(&self) -> f64 {
        self.other_work_j
    }
    /// Validates the complete increment and candidate total before changing state.
    pub fn record(&mut self, step: DissipationStep) -> Result<(), TriboError> {
        step.validate()?;
        let dissipated_work_j = checked_add(
            self.dissipated_work_j,
            step.total_work_j,
            "dissipated_work_j",
        )?;
        let surface_a_heat_j = checked_add(
            self.surface_a_heat_j,
            step.surface_a_heat_j,
            "surface_a_heat_j",
        )?;
        let surface_b_heat_j = checked_add(
            self.surface_b_heat_j,
            step.surface_b_heat_j,
            "surface_b_heat_j",
        )?;
        let other_work_j = checked_add(self.other_work_j, step.other_work_j, "other_work_j")?;
        nonnegative_finite(dissipated_work_j, "dissipated_work_j")?;
        nonnegative_finite(surface_a_heat_j, "surface_a_heat_j")?;
        nonnegative_finite(surface_b_heat_j, "surface_b_heat_j")?;
        nonnegative_finite(other_work_j, "other_work_j")?;
        let channels = checked_add(
            checked_add(surface_a_heat_j, surface_b_heat_j, "work_ledger_channels")?,
            other_work_j,
            "work_ledger_channels",
        )?;
        if (channels - dissipated_work_j).abs() > EPSILON * dissipated_work_j.max(1.0) {
            return Err(TriboError::InvalidDissipationStep {
                field: "ledger_closure",
            });
        }
        self.dissipated_work_j = dissipated_work_j;
        self.surface_a_heat_j = surface_a_heat_j;
        self.surface_b_heat_j = surface_b_heat_j;
        self.other_work_j = other_work_j;
        Ok(())
    }
}

/// Caller-supplied dimensionless Archard coefficient with no calibration claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArchardLaw {
    pub wear_coefficient: f64,
}
/// Caller-owned wear volume in cubic metres.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WearState {
    volume_m3: f64,
}
impl WearState {
    pub fn new(volume_m3: f64) -> Result<Self, TriboError> {
        nonnegative_finite(volume_m3, "wear_volume_m3")?;
        Ok(Self { volume_m3 })
    }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
}
/// Result of one Archard update, retaining caller provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct WearAdvance {
    pub increment_m3: f64,
    pub total_volume_m3: f64,
    provenance: InputProvenance,
}
impl WearAdvance {
    #[must_use]
    pub fn provenance(&self) -> &InputProvenance {
        &self.provenance
    }
}
impl ArchardLaw {
    /// Computes and validates `k * F * ds / H` before atomically committing the state candidate.
    pub fn advance(
        self,
        interface: &InterfaceSystemRef,
        state: &mut WearState,
        normal_force_n: f64,
        sliding_distance_m: f64,
        hardness_pa: f64,
    ) -> Result<WearAdvance, TriboError> {
        interface.validate()?;
        nonnegative_finite(self.wear_coefficient, "wear_coefficient")?;
        nonnegative_finite(normal_force_n, "normal_force_n")?;
        nonnegative_finite(sliding_distance_m, "sliding_distance_m")?;
        positive_finite(hardness_pa, "hardness_pa")?;
        nonnegative_finite(state.volume_m3, "wear_volume_m3")?;
        let increment_m3 = checked_mul(
            self.wear_coefficient,
            checked_mul(normal_force_n, sliding_distance_m, "wear_increment")?,
            "wear_increment",
        )? / hardness_pa;
        nonnegative_finite(increment_m3, "wear_increment")?;
        let total_volume_m3 = checked_add(state.volume_m3, increment_m3, "wear_volume_m3")?;
        nonnegative_finite(total_volume_m3, "wear_volume_m3")?;
        state.volume_m3 = total_volume_m3;
        Ok(WearAdvance {
            increment_m3,
            total_volume_m3,
            provenance: interface.provenance.clone(),
        })
    }
}

fn nonblank(value: &str, field: &'static str) -> Result<(), TriboError> {
    if value.trim().is_empty() {
        Err(TriboError::MissingIdentity { field })
    } else {
        Ok(())
    }
}
fn finite(value: f64, field: &'static str) -> Result<(), TriboError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(TriboError::InvalidInput { field })
    }
}
fn nonnegative_finite(value: f64, field: &'static str) -> Result<(), TriboError> {
    finite(value, field)?;
    if value >= 0.0 {
        Ok(())
    } else {
        Err(TriboError::InvalidInput { field })
    }
}
fn positive_finite(value: f64, field: &'static str) -> Result<(), TriboError> {
    finite(value, field)?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(TriboError::InvalidInput { field })
    }
}
fn finite_vec(values: [f64; 3], field: &'static str) -> Result<(), TriboError> {
    if values.into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(TriboError::NonFiniteVector { field })
    }
}
fn checked_mul(a: f64, b: f64, field: &'static str) -> Result<f64, TriboError> {
    let value = a * b;
    finite(value, field)?;
    Ok(value)
}
fn checked_add(a: f64, b: f64, field: &'static str) -> Result<f64, TriboError> {
    let value = a + b;
    finite(value, field)?;
    Ok(value)
}
fn checked_div(a: f64, b: f64, field: &'static str) -> Result<f64, TriboError> {
    positive_finite(b, field)?;
    let value = a / b;
    finite(value, field)?;
    Ok(value)
}
fn checked_sqrt(value: f64, field: &'static str) -> Result<f64, TriboError> {
    nonnegative_finite(value, field)?;
    let root = value.sqrt();
    finite(root, field)?;
    Ok(root)
}
fn stable_norm(values: [f64; 3], field: &'static str) -> Result<f64, TriboError> {
    let value = values[0].hypot(values[1]).hypot(values[2]);
    finite(value, field)?;
    Ok(value)
}
fn stable_dot(a: [f64; 3], b: [f64; 3], field: &'static str) -> Result<f64, TriboError> {
    let value = a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]));
    finite(value, field)?;
    Ok(value)
}
fn scale_checked(
    values: [f64; 3],
    scalar: f64,
    field: &'static str,
) -> Result<[f64; 3], TriboError> {
    let output = [
        checked_mul(values[0], scalar, field)?,
        checked_mul(values[1], scalar, field)?,
        checked_mul(values[2], scalar, field)?,
    ];
    finite_vec(output, field)?;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn declared_dry() -> InterfaceSystemRef {
        InterfaceSystemRef::new(
            "fixture/a->b",
            "fixture/history",
            "fixture/source",
            InputAuthority::CallerDeclared,
            InterfaceMedium::Dry,
        )
        .unwrap()
    }
    fn synthetic_dry() -> InterfaceSystemRef {
        InterfaceSystemRef::new(
            "fixture/a->b",
            "fixture/history",
            "fixture/synthetic",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap()
    }
    fn tangent(v: [f64; 3]) -> TangentialSlip {
        TangentialSlip::new(&ContactFrame::new([0.0, 0.0, 1.0]).unwrap(), v).unwrap()
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g0_refuses_bad_interface_and_normal_slip() {
        assert!(matches!(
            InterfaceSystemRef::new(
                "",
                "h",
                "s",
                InputAuthority::Estimated,
                InterfaceMedium::Dry
            ),
            Err(TriboError::MissingIdentity { .. })
        ));
        assert!(matches!(
            InterfaceSystemRef::new(
                "a",
                "h",
                "s",
                InputAuthority::Estimated,
                InterfaceMedium::Lubricated
            ),
            Err(TriboError::NotDryInterface { .. })
        ));
        let frame = ContactFrame::new([0.0, 0.0, 2.0]).unwrap();
        assert!(matches!(
            TangentialSlip::new(&frame, [1.0, 0.0, 1e-3]),
            Err(TriboError::NormalSlipComponent { .. })
        ));
    }

    #[test]
    fn g1_rigid_block_static_capacity_matches_incline_equilibrium() {
        let law = FrictionLaw::Coulomb {
            static_mu: 0.5,
            kinetic_mu: 0.4,
        };
        let response = law
            .evaluate(&synthetic_dry(), 20.0, tangent([0.0; 3]))
            .unwrap();
        assert_eq!(response.regime, FrictionRegime::Sticking);
        close(response.static_limit, 10.0);
        close((response.static_limit / 20.0).atan().tan(), 0.5);
    }

    #[test]
    fn g3_friction_reversal_scaling_and_provenance_hold() {
        let law = FrictionLaw::Stribeck {
            static_mu: 0.6,
            kinetic_mu: 0.3,
            characteristic_speed: 2.0,
            viscous_per_speed: 0.01,
        };
        let a = law
            .evaluate(&synthetic_dry(), 10.0, tangent([2.0, 0.0, 0.0]))
            .unwrap();
        let b = law
            .evaluate(&synthetic_dry(), 20.0, tangent([2.0, 0.0, 0.0]))
            .unwrap();
        let reversed = law
            .evaluate(&synthetic_dry(), 10.0, tangent([-2.0, 0.0, 0.0]))
            .unwrap();
        close(b.traction_n()[0], 2.0 * a.traction_n()[0]);
        close(reversed.traction_n()[0], -a.traction_n()[0]);
        close(reversed.dissipated_power_w(), a.dissipated_power_w());
        assert_eq!(a.provenance().authority(), InputAuthority::SyntheticFixture);
        assert_eq!(a.provenance().source_id(), "fixture/synthetic");
    }

    #[test]
    fn g1_hertz_independent_numeric_values_and_cross_relations() {
        let sphere = HertzSpherePlane {
            effective_radius: 0.02,
            reduced_modulus: 200e9,
        }
        .response(&declared_dry(), 2e-6)
        .unwrap();
        close(sphere.contact_radius_m, 2e-4);
        close(sphere.normal_force_n, 106.666_666_666_666_67);
        close(sphere.peak_pressure_pa, 1.273_239_544_735_162_8e9);
        close(
            sphere.normal_force_n,
            (2.0 / 3.0)
                * core::f64::consts::PI
                * sphere.contact_radius_m.powi(2)
                * sphere.peak_pressure_pa,
        );
        let cylinder = HertzCylinderPlane {
            effective_radius: 0.02,
            reduced_modulus: 200e9,
        }
        .response(&declared_dry(), 1_000.0)
        .unwrap();
        close(cylinder.half_width_m, 1.128_379_167_095_512_6e-5);
        close(cylinder.peak_pressure_pa, 56_418_958.354_775_63);
        close(
            1_000.0,
            0.5 * core::f64::consts::PI * cylinder.half_width_m * cylinder.peak_pressure_pa,
        );
    }

    #[test]
    fn g0_derived_overflow_refuses_without_state_mutation() {
        let partition = HeatPartition::new(0.5, 0.5, 0.0).unwrap();
        assert!(matches!(
            DissipationStep::from_power(f64::MAX, 2.0, partition),
            Err(TriboError::InvalidInput { .. })
        ));
        let mut ledger = WorkLedger {
            dissipated_work_j: f64::MAX,
            ..WorkLedger::default()
        };
        let step = DissipationStep::from_power(f64::MAX, 1.0, partition).unwrap();
        assert!(ledger.record(step).is_err());
        assert_eq!(ledger.dissipated_work_j(), f64::MAX);
        let mut wear = WearState::new(1.0).unwrap();
        let result = ArchardLaw {
            wear_coefficient: 2.0,
        }
        .advance(&declared_dry(), &mut wear, f64::MAX, 1.0, 1.0);
        assert!(result.is_err());
        assert_eq!(wear.volume_m3(), 1.0);
    }

    #[test]
    fn g0_forged_negative_or_nonfinite_work_is_rejected_without_ledger_mutation() {
        let mut ledger = WorkLedger::default();
        let negative = DissipationStep {
            total_work_j: 1.0,
            surface_a_heat_j: -1.0,
            surface_b_heat_j: 2.0,
            other_work_j: 0.0,
        };
        assert!(matches!(
            ledger.record(negative),
            Err(TriboError::InvalidDissipationStep { .. })
        ));
        assert_eq!(ledger.dissipated_work_j(), 0.0);
        let nonfinite = DissipationStep {
            total_work_j: 1.0,
            surface_a_heat_j: f64::INFINITY,
            surface_b_heat_j: 0.0,
            other_work_j: 0.0,
        };
        assert!(matches!(
            ledger.record(nonfinite),
            Err(TriboError::InvalidDissipationStep { .. })
        ));
        assert_eq!(ledger.dissipated_work_j(), 0.0);
    }

    #[test]
    fn g3_stable_norm_and_all_derived_extremes_refuse_or_remain_finite() {
        let law = FrictionLaw::Coulomb {
            static_mu: 0.5,
            kinetic_mu: 0.5,
        };
        let response = law
            .evaluate(&declared_dry(), 1.0, tangent([f64::MAX, 0.0, 0.0]))
            .unwrap();
        assert!(response.dissipated_power_w().is_finite());
        close(response.traction_n()[0], -0.5);
        let overflowing_law = FrictionLaw::Coulomb {
            static_mu: 2.0,
            kinetic_mu: 2.0,
        };
        assert!(matches!(
            overflowing_law.evaluate(&declared_dry(), f64::MAX, tangent([1.0, 0.0, 0.0])),
            Err(TriboError::InvalidInput { .. })
        ));
        assert!(
            HertzSpherePlane {
                effective_radius: f64::MAX,
                reduced_modulus: f64::MAX
            }
            .response(&declared_dry(), 1.0)
            .is_err()
        );
    }

    #[test]
    fn g0_resistance_signs_and_checked_partition_closure() {
        let input = ResistanceInput {
            normal_force_n: 100.0,
            angular_speed_rad_s: -3.0,
            contour_speed_mps: 4.0,
        };
        let rolling = ConstantRollingMoment {
            moment_arm_m: 0.002,
        }
        .evaluate(&synthetic_dry(), input)
        .unwrap();
        let contour = ConstantContourForce { force_n: 7.0 }
            .evaluate(&synthetic_dry(), input)
            .unwrap();
        assert!(rolling.rolling_moment_n_m > 0.0 && rolling.dissipated_power_w() > 0.0);
        assert!(contour.contour_force_n < 0.0 && contour.dissipated_power_w() > 0.0);
        let step =
            DissipationStep::from_power(80.0, 0.5, HeatPartition::new(0.25, 0.5, 0.25).unwrap())
                .unwrap();
        close(step.total_work_j(), 40.0);
        close(
            step.surface_a_heat_j() + step.surface_b_heat_j() + step.other_work_j(),
            40.0,
        );
    }

    #[test]
    fn g3_archard_scales_and_retains_caller_ceiling() {
        let law = ArchardLaw {
            wear_coefficient: 2e-6,
        };
        let mut state = WearState::default();
        let a = law
            .advance(&synthetic_dry(), &mut state, 1_000.0, 5.0, 2e9)
            .unwrap();
        let b = law
            .advance(&synthetic_dry(), &mut state, 2_000.0, 5.0, 2e9)
            .unwrap();
        close(b.increment_m3, 2.0 * a.increment_m3);
        close(state.volume_m3(), 3.0 * a.increment_m3);
        assert_eq!(b.provenance().authority(), InputAuthority::SyntheticFixture);
    }
}
