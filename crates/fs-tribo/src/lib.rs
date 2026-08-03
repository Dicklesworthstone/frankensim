#![forbid(unsafe_code)]
//! Typed dry-contact baseline with explicit authority and refusal boundaries.
//!
//! The crate deliberately does not provide a material-pair coefficient table. A caller must bind an
//! ordered interface system and declare whether its inputs are admitted or synthetic. Synthetic
//! inputs are useful for analytic tests, but never become material authority.

use core::fmt;

const EPSILON: f64 = 64.0 * f64::EPSILON;

/// A complete ordered interface identity and its authority source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InterfaceSystemRef {
    /// Ordered surface/system identity; reversing surfaces requires a new identity.
    pub ordered_system_id: String,
    /// Named history-state identity (wear, third body, and prior treatment are not implicit).
    pub history_id: String,
    /// Whether this input may support an admitted evaluation.
    pub authority: InputAuthority,
    /// Dry laws refuse any other medium declaration.
    pub medium: InterfaceMedium,
}

impl InterfaceSystemRef {
    /// Constructs a checked interface reference.
    pub fn new(
        ordered_system_id: impl Into<String>,
        history_id: impl Into<String>,
        authority: InputAuthority,
        medium: InterfaceMedium,
    ) -> Result<Self, TriboError> {
        let value = Self {
            ordered_system_id: ordered_system_id.into(),
            history_id: history_id.into(),
            authority,
            medium,
        };
        value.validate()?;
        Ok(value)
    }

    fn validate(&self) -> Result<(), TriboError> {
        if self.ordered_system_id.trim().is_empty() {
            return Err(TriboError::MissingInterfaceIdentity);
        }
        if self.history_id.trim().is_empty() {
            return Err(TriboError::MissingHistoryIdentity);
        }
        if self.medium != InterfaceMedium::Dry {
            return Err(TriboError::NotDryInterface {
                medium: self.medium,
            });
        }
        Ok(())
    }

    fn require_analytic_authority(&self) -> Result<(), TriboError> {
        self.validate()?;
        if self.authority != InputAuthority::Admitted {
            return Err(TriboError::UnsupportedAuthority {
                authority: self.authority,
            });
        }
        Ok(())
    }
}

/// Authority carried by a caller-provided card or fixture.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputAuthority {
    /// A caller has bound the input to an admitted ordered system card.
    Admitted,
    /// A declared synthetic fixture, suitable only for tests/oracles.
    SyntheticFixture,
    /// Data is present but has not been admitted for use.
    Estimated,
}

/// Medium classification at the dry-law boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterfaceMedium {
    /// Explicitly dry contact.
    Dry,
    /// A fluid or lubrication film is present; dry laws must refuse.
    Lubricated,
    /// The caller has not identified the medium.
    Undeclared,
}

/// Total refusals for this baseline.
#[derive(Debug, Clone, PartialEq)]
pub enum TriboError {
    /// Interface identity is absent.
    MissingInterfaceIdentity,
    /// History identity is absent.
    MissingHistoryIdentity,
    /// A dry law received an incompatible medium.
    NotDryInterface { medium: InterfaceMedium },
    /// The caller has not supplied the authority required for an analytic response.
    UnsupportedAuthority { authority: InputAuthority },
    /// A scalar is non-finite or outside its physical domain.
    InvalidInput { field: &'static str },
    /// A three-vector is not finite.
    NonFiniteVector { field: &'static str },
    /// A heat partition does not close.
    InvalidHeatPartition { sum: f64 },
    /// A zero-speed direction was requested where a signed loss is indeterminate.
    IndeterminateDirection,
}

impl fmt::Display for TriboError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingInterfaceIdentity => {
                f.write_str("ordered interface-system identity is required")
            }
            Self::MissingHistoryIdentity => {
                f.write_str("named interface history identity is required")
            }
            Self::NotDryInterface { medium } => write!(f, "dry law refuses {medium:?} interface"),
            Self::UnsupportedAuthority { authority } => {
                write!(
                    f,
                    "analytic response requires admitted authority, got {authority:?}"
                )
            }
            Self::InvalidInput { field } => write!(f, "invalid physical input: {field}"),
            Self::NonFiniteVector { field } => write!(f, "non-finite vector: {field}"),
            Self::InvalidHeatPartition { sum } => {
                write!(f, "heat partition must sum to one, got {sum}")
            }
            Self::IndeterminateDirection => {
                f.write_str("signed resistance requires non-zero speed")
            }
        }
    }
}

impl std::error::Error for TriboError {}

/// A friction constitutive ladder. All coefficients are dimensionless except documented speed terms.
#[derive(Debug, Clone, PartialEq)]
pub enum FrictionLaw {
    /// Coulomb static and kinetic limits.
    Coulomb { static_mu: f64, kinetic_mu: f64 },
    /// Kinetic coefficient `mu_zero + slope_per_speed * |v|`, clamped at zero.
    VelocityDependent {
        static_mu: f64,
        mu_zero: f64,
        slope_per_speed: f64,
    },
    /// Exponential Stribeck decay plus a non-negative viscous term.
    Stribeck {
        static_mu: f64,
        kinetic_mu: f64,
        characteristic_speed: f64,
        viscous_per_speed: f64,
    },
}

impl FrictionLaw {
    /// Evaluates the law. At zero slip only the stick capacity is reported; the contact solver owns
    /// the actual reaction traction.
    pub fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        normal_force: f64,
        slip_velocity: [f64; 3],
    ) -> Result<FrictionResponse, TriboError> {
        interface.validate()?;
        nonnegative_finite(normal_force, "normal_force")?;
        finite_vec(slip_velocity, "slip_velocity")?;
        self.validate()?;
        let speed = norm(slip_velocity);
        let static_limit = self.static_mu() * normal_force;
        if speed == 0.0 {
            return Ok(FrictionResponse {
                regime: FrictionRegime::Sticking,
                static_limit,
                kinetic_coefficient: None,
                traction: [0.0; 3],
                dissipated_power: 0.0,
            });
        }
        let coefficient = self.kinetic_mu(speed);
        let magnitude = coefficient * normal_force;
        let traction = scale(slip_velocity, -magnitude / speed);
        let dissipated_power = -dot(traction, slip_velocity);
        Ok(FrictionResponse {
            regime: FrictionRegime::Sliding,
            static_limit,
            kinetic_coefficient: Some(coefficient),
            traction,
            dissipated_power,
        })
    }

    fn validate(&self) -> Result<(), TriboError> {
        match *self {
            Self::Coulomb {
                static_mu,
                kinetic_mu,
            } => {
                nonnegative_finite(static_mu, "static_mu")?;
                nonnegative_finite(kinetic_mu, "kinetic_mu")?;
                if kinetic_mu > static_mu {
                    return Err(TriboError::InvalidInput {
                        field: "kinetic_mu > static_mu",
                    });
                }
            }
            Self::VelocityDependent {
                static_mu,
                mu_zero,
                slope_per_speed,
            } => {
                nonnegative_finite(static_mu, "static_mu")?;
                nonnegative_finite(mu_zero, "mu_zero")?;
                finite(slope_per_speed, "slope_per_speed")?;
                if mu_zero > static_mu {
                    return Err(TriboError::InvalidInput {
                        field: "mu_zero > static_mu",
                    });
                }
            }
            Self::Stribeck {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                viscous_per_speed,
            } => {
                nonnegative_finite(static_mu, "static_mu")?;
                nonnegative_finite(kinetic_mu, "kinetic_mu")?;
                positive_finite(characteristic_speed, "characteristic_speed")?;
                nonnegative_finite(viscous_per_speed, "viscous_per_speed")?;
                if kinetic_mu > static_mu {
                    return Err(TriboError::InvalidInput {
                        field: "kinetic_mu > static_mu",
                    });
                }
            }
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

    fn kinetic_mu(&self, speed: f64) -> f64 {
        match *self {
            Self::Coulomb { kinetic_mu, .. } => kinetic_mu,
            Self::VelocityDependent {
                mu_zero,
                slope_per_speed,
                ..
            } => (mu_zero + slope_per_speed * speed).max(0.0),
            Self::Stribeck {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                viscous_per_speed,
            } => {
                kinetic_mu
                    + (static_mu - kinetic_mu) * (-(speed / characteristic_speed).powi(2)).exp()
                    + viscous_per_speed * speed
            }
        }
    }
}

/// Whether the response leaves traction to the contact solve or supplies sliding traction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrictionRegime {
    /// Stick capacity only; no fictitious reaction is emitted.
    Sticking,
    /// Gross sliding traction is emitted.
    Sliding,
}

/// Typed result of one friction-law query.
#[derive(Debug, Clone, PartialEq)]
pub struct FrictionResponse {
    /// Stick threshold `mu_s * N` in newtons.
    pub static_limit: f64,
    /// Kinetic coefficient during sliding only.
    pub kinetic_coefficient: Option<f64>,
    /// Tangential traction force in newtons, opposing the supplied slip velocity.
    pub traction: [f64; 3],
    /// Non-negative dissipated power in watts.
    pub dissipated_power: f64,
    /// Response branch.
    pub regime: FrictionRegime,
}

/// Hertz sphere/plane parameters, in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzSpherePlane {
    /// Effective sphere radius in metres.
    pub effective_radius: f64,
    /// Reduced elastic modulus in pascals.
    pub reduced_modulus: f64,
}

impl HertzSpherePlane {
    /// Closed-form elastic response for non-negative indentation.
    pub fn response(
        self,
        interface: &InterfaceSystemRef,
        indentation: f64,
    ) -> Result<HertzSphereResponse, TriboError> {
        interface.require_analytic_authority()?;
        positive_finite(self.effective_radius, "effective_radius")?;
        positive_finite(self.reduced_modulus, "reduced_modulus")?;
        nonnegative_finite(indentation, "indentation")?;
        let contact_radius = (self.effective_radius * indentation).sqrt();
        let normal_force = (4.0 / 3.0)
            * self.reduced_modulus
            * self.effective_radius.sqrt()
            * indentation.powf(1.5);
        let peak_pressure = if contact_radius == 0.0 {
            0.0
        } else {
            3.0 * normal_force / (2.0 * core::f64::consts::PI * contact_radius.powi(2))
        };
        Ok(HertzSphereResponse {
            normal_force,
            contact_radius,
            peak_pressure,
        })
    }
}

/// Hertz sphere/plane response in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzSphereResponse {
    /// Normal force in newtons.
    pub normal_force: f64,
    /// Circular contact radius in metres.
    pub contact_radius: f64,
    /// Maximum Hertz pressure in pascals.
    pub peak_pressure: f64,
}

/// Hertz cylinder/plane parameters, expressed per unit axial length in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzCylinderPlane {
    /// Effective cylinder radius in metres.
    pub effective_radius: f64,
    /// Reduced elastic modulus in pascals.
    pub reduced_modulus: f64,
}

impl HertzCylinderPlane {
    /// Closed-form line-contact response for a non-negative normal line load in N/m.
    pub fn response(
        self,
        interface: &InterfaceSystemRef,
        normal_line_load: f64,
    ) -> Result<HertzCylinderResponse, TriboError> {
        interface.require_analytic_authority()?;
        positive_finite(self.effective_radius, "effective_radius")?;
        positive_finite(self.reduced_modulus, "reduced_modulus")?;
        nonnegative_finite(normal_line_load, "normal_line_load")?;
        let half_width = (4.0 * normal_line_load * self.effective_radius
            / (core::f64::consts::PI * self.reduced_modulus))
            .sqrt();
        let peak_pressure = if half_width == 0.0 {
            0.0
        } else {
            2.0 * normal_line_load / (core::f64::consts::PI * half_width)
        };
        Ok(HertzCylinderResponse {
            half_width,
            peak_pressure,
        })
    }
}

/// Hertz cylinder/plane response in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HertzCylinderResponse {
    /// Contact half-width in metres.
    pub half_width: f64,
    /// Peak line-contact pressure in pascals.
    pub peak_pressure: f64,
}

/// An interface for generic rolling, contour, or torsional loss models.
pub trait ResistanceLaw {
    /// Evaluates a signed resistance and its non-negative dissipation rate.
    fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        input: ResistanceInput,
    ) -> Result<ResistanceResponse, TriboError>;
}

/// Input shared by the included resistance rungs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResistanceInput {
    /// Compressive normal force in newtons.
    pub normal_force: f64,
    /// Signed angular speed in rad/s for rolling moment.
    pub angular_speed: f64,
    /// Signed contour speed in m/s for contour force.
    pub contour_speed: f64,
}

impl ResistanceInput {
    fn validate(self) -> Result<(), TriboError> {
        nonnegative_finite(self.normal_force, "normal_force")?;
        finite(self.angular_speed, "angular_speed")?;
        finite(self.contour_speed, "contour_speed")
    }
}

/// A resistance output. Exactly one generalized effort is non-zero for each supplied rung.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResistanceResponse {
    /// Rolling moment in N m, opposite angular speed.
    pub rolling_moment: f64,
    /// Contour force in N, opposite contour speed.
    pub contour_force: f64,
    /// Non-negative loss power in W.
    pub dissipated_power: f64,
}

/// Generic constant rolling-moment arm; the caller supplies the coefficient and authority.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantRollingMoment {
    /// Moment arm in metres; magnitude is `normal_force * moment_arm`.
    pub moment_arm: f64,
}

impl ResistanceLaw for ConstantRollingMoment {
    fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        input: ResistanceInput,
    ) -> Result<ResistanceResponse, TriboError> {
        interface.validate()?;
        positive_finite(self.moment_arm, "moment_arm")?;
        input.validate()?;
        if input.angular_speed == 0.0 {
            return Ok(ResistanceResponse {
                rolling_moment: 0.0,
                contour_force: 0.0,
                dissipated_power: 0.0,
            });
        }
        let rolling_moment = -input.angular_speed.signum() * input.normal_force * self.moment_arm;
        Ok(ResistanceResponse {
            rolling_moment,
            contour_force: 0.0,
            dissipated_power: -rolling_moment * input.angular_speed,
        })
    }
}

/// Generic constant contour-resistance force; it does not encode a geometry-specific calibration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantContourForce {
    /// Non-negative contour force magnitude in newtons.
    pub force: f64,
}

impl ResistanceLaw for ConstantContourForce {
    fn evaluate(
        &self,
        interface: &InterfaceSystemRef,
        input: ResistanceInput,
    ) -> Result<ResistanceResponse, TriboError> {
        interface.validate()?;
        nonnegative_finite(self.force, "force")?;
        input.validate()?;
        if input.contour_speed == 0.0 {
            return Ok(ResistanceResponse {
                rolling_moment: 0.0,
                contour_force: 0.0,
                dissipated_power: 0.0,
            });
        }
        let contour_force = -input.contour_speed.signum() * self.force;
        Ok(ResistanceResponse {
            rolling_moment: 0.0,
            contour_force,
            dissipated_power: -contour_force * input.contour_speed,
        })
    }
}

/// Declared closure of one frictional heat partition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatPartition {
    /// Fraction of frictional work assigned to surface A.
    pub surface_a: f64,
    /// Fraction of frictional work assigned to surface B.
    pub surface_b: f64,
    /// Explicit residual share (third body, radiation, or another declared port).
    pub other: f64,
}

impl HeatPartition {
    /// Makes a checked partition. The named `other` share prevents silent energy disappearance.
    pub fn new(surface_a: f64, surface_b: f64, other: f64) -> Result<Self, TriboError> {
        nonnegative_finite(surface_a, "surface_a")?;
        nonnegative_finite(surface_b, "surface_b")?;
        nonnegative_finite(other, "other")?;
        let sum = surface_a + surface_b + other;
        if (sum - 1.0).abs() > EPSILON {
            return Err(TriboError::InvalidHeatPartition { sum });
        }
        Ok(Self {
            surface_a,
            surface_b,
            other,
        })
    }
}

/// One declared dissipation increment, in joules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DissipationStep {
    /// Total non-negative friction work.
    pub total_work: f64,
    /// Heat assigned to surface A.
    pub surface_a_heat: f64,
    /// Heat assigned to surface B.
    pub surface_b_heat: f64,
    /// Explicit other loss/heat channel.
    pub other_work: f64,
}

impl DissipationStep {
    /// Converts non-negative dissipated power over a non-negative duration into a closed work split.
    pub fn from_power(
        dissipated_power: f64,
        duration: f64,
        partition: HeatPartition,
    ) -> Result<Self, TriboError> {
        nonnegative_finite(dissipated_power, "dissipated_power")?;
        nonnegative_finite(duration, "duration")?;
        let total_work = dissipated_power * duration;
        Ok(Self {
            total_work,
            surface_a_heat: total_work * partition.surface_a,
            surface_b_heat: total_work * partition.surface_b,
            other_work: total_work * partition.other,
        })
    }

    /// Checks that the declared channels reproduce total work within roundoff.
    pub fn closes(self) -> bool {
        (self.surface_a_heat + self.surface_b_heat + self.other_work - self.total_work).abs()
            <= EPSILON * self.total_work.max(1.0)
    }
}

/// Caller-owned monotonically accumulated dissipated work.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WorkLedger {
    /// Cumulative work in joules.
    pub dissipated_work: f64,
}

impl WorkLedger {
    /// Records a closed non-negative work increment.
    pub fn record(&mut self, step: DissipationStep) -> Result<(), TriboError> {
        nonnegative_finite(step.total_work, "total_work")?;
        if !step.closes() {
            return Err(TriboError::InvalidHeatPartition {
                sum: (step.surface_a_heat + step.surface_b_heat + step.other_work)
                    / step.total_work.max(1.0),
            });
        }
        self.dissipated_work += step.total_work;
        finite(self.dissipated_work, "dissipated_work")
    }
}

/// Archard wear coefficient and its explicit caller authority binding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArchardLaw {
    /// Dimensionless wear coefficient; no default is supplied.
    pub wear_coefficient: f64,
}

/// Caller-owned Archard wear state.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WearState {
    /// Cumulative wear volume in m^3.
    pub volume: f64,
}

impl ArchardLaw {
    /// Evolves the supplied state by `k * normal_force * sliding_distance / hardness`.
    pub fn advance(
        self,
        interface: &InterfaceSystemRef,
        state: &mut WearState,
        normal_force: f64,
        sliding_distance: f64,
        hardness: f64,
    ) -> Result<f64, TriboError> {
        interface.validate()?;
        nonnegative_finite(self.wear_coefficient, "wear_coefficient")?;
        nonnegative_finite(normal_force, "normal_force")?;
        nonnegative_finite(sliding_distance, "sliding_distance")?;
        positive_finite(hardness, "hardness")?;
        nonnegative_finite(state.volume, "wear_state.volume")?;
        let increment = self.wear_coefficient * normal_force * sliding_distance / hardness;
        state.volume += increment;
        finite(state.volume, "wear_state.volume")?;
        Ok(increment)
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

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
fn scale(v: [f64; 3], scalar: f64) -> [f64; 3] {
    [v[0] * scalar, v[1] * scalar, v[2] * scalar]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn admitted_dry() -> InterfaceSystemRef {
        InterfaceSystemRef::new(
            "fixture/steel-a->steel-b",
            "fixture/history-0",
            InputAuthority::Admitted,
            InterfaceMedium::Dry,
        )
        .unwrap()
    }

    fn synthetic_dry() -> InterfaceSystemRef {
        InterfaceSystemRef::new(
            "fixture/synthetic-a->b",
            "fixture/history-0",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g0_refuses_missing_identity_medium_and_unadmitted_hertz() {
        assert_eq!(
            InterfaceSystemRef::new("", "h", InputAuthority::Admitted, InterfaceMedium::Dry),
            Err(TriboError::MissingInterfaceIdentity)
        );
        assert!(matches!(
            InterfaceSystemRef::new(
                "a->b",
                "h",
                InputAuthority::Admitted,
                InterfaceMedium::Lubricated
            ),
            Err(TriboError::NotDryInterface { .. })
        ));
        let hertz = HertzSpherePlane {
            effective_radius: 0.01,
            reduced_modulus: 100e9,
        };
        assert!(matches!(
            hertz.response(&synthetic_dry(), 1e-6),
            Err(TriboError::UnsupportedAuthority { .. })
        ));
    }

    #[test]
    fn g0_coulomb_stick_slip_threshold_is_exact_for_rigid_block() {
        let law = FrictionLaw::Coulomb {
            static_mu: 0.5,
            kinetic_mu: 0.4,
        };
        let normal = 20.0;
        let stick = law.evaluate(&synthetic_dry(), normal, [0.0; 3]).unwrap();
        assert_eq!(stick.regime, FrictionRegime::Sticking);
        close(stick.static_limit, 10.0);
        // tan(theta)=mu_s is the analytic rigid-block threshold; the fixture is synthetic.
        close((stick.static_limit / normal).atan().tan(), 0.5);
        let slide = law
            .evaluate(&synthetic_dry(), normal, [3.0, 4.0, 0.0])
            .unwrap();
        assert_eq!(slide.regime, FrictionRegime::Sliding);
        close(slide.traction[0], -4.8);
        close(slide.traction[1], -6.4);
        close(slide.dissipated_power, 40.0);
    }

    #[test]
    fn g1_hertz_sphere_and_cylinder_match_closed_forms_and_limits() {
        let interface = admitted_dry();
        let sphere = HertzSpherePlane {
            effective_radius: 0.02,
            reduced_modulus: 200e9,
        };
        let r = sphere.response(&interface, 2e-6).unwrap();
        close(r.contact_radius, (0.02_f64 * 2e-6).sqrt());
        close(
            r.normal_force,
            (4.0 / 3.0) * 200e9 * 0.02_f64.sqrt() * (2e-6_f64).powf(1.5),
        );
        assert!(r.peak_pressure > 0.0);
        let cylinder = HertzCylinderPlane {
            effective_radius: 0.02,
            reduced_modulus: 200e9,
        };
        let c = cylinder.response(&interface, 1_000.0).unwrap();
        close(
            c.half_width,
            (4.0 * 1_000.0 * 0.02 / (core::f64::consts::PI * 200e9)).sqrt(),
        );
        assert!(c.peak_pressure > 0.0);
        assert_eq!(sphere.response(&interface, 0.0).unwrap().normal_force, 0.0);
        assert_eq!(cylinder.response(&interface, 0.0).unwrap().half_width, 0.0);
    }

    #[test]
    fn g3_friction_reversal_and_scaling_are_metamorphic() {
        let law = FrictionLaw::Stribeck {
            static_mu: 0.6,
            kinetic_mu: 0.3,
            characteristic_speed: 2.0,
            viscous_per_speed: 0.01,
        };
        let a = law
            .evaluate(&synthetic_dry(), 10.0, [2.0, 0.0, 0.0])
            .unwrap();
        let b = law
            .evaluate(&synthetic_dry(), 20.0, [2.0, 0.0, 0.0])
            .unwrap();
        let reversed = law
            .evaluate(&synthetic_dry(), 10.0, [-2.0, 0.0, 0.0])
            .unwrap();
        close(b.traction[0], 2.0 * a.traction[0]);
        close(reversed.traction[0], -a.traction[0]);
        close(reversed.dissipated_power, a.dissipated_power);
    }

    #[test]
    fn g0_resistance_opposes_motion_and_never_creates_power() {
        let input = ResistanceInput {
            normal_force: 100.0,
            angular_speed: -3.0,
            contour_speed: 4.0,
        };
        let rolling = ConstantRollingMoment { moment_arm: 0.002 }
            .evaluate(&synthetic_dry(), input)
            .unwrap();
        let contour = ConstantContourForce { force: 7.0 }
            .evaluate(&synthetic_dry(), input)
            .unwrap();
        assert!(rolling.rolling_moment > 0.0 && rolling.dissipated_power > 0.0);
        assert!(contour.contour_force < 0.0 && contour.dissipated_power > 0.0);
    }

    #[test]
    fn g0_heat_partition_closes_and_replay_is_deterministic() {
        let partition = HeatPartition::new(0.25, 0.5, 0.25).unwrap();
        let step = DissipationStep::from_power(80.0, 0.5, partition).unwrap();
        assert!(step.closes());
        close(step.total_work, 40.0);
        close(
            step.surface_a_heat + step.surface_b_heat + step.other_work,
            40.0,
        );
        let mut a = WorkLedger::default();
        let mut b = WorkLedger::default();
        for _ in 0..3 {
            a.record(step).unwrap();
        }
        for _ in 0..3 {
            b.record(step).unwrap();
        }
        assert_eq!(a, b);
        assert_eq!(
            HeatPartition::new(0.2, 0.2, 0.2),
            Err(TriboError::InvalidHeatPartition {
                sum: 0.6000000000000001
            })
        );
    }

    #[test]
    fn g3_archard_scales_and_accumulates_without_material_defaults() {
        let law = ArchardLaw {
            wear_coefficient: 2e-6,
        };
        let mut state = WearState::default();
        let first = law
            .advance(&synthetic_dry(), &mut state, 1_000.0, 5.0, 2e9)
            .unwrap();
        let second = law
            .advance(&synthetic_dry(), &mut state, 2_000.0, 5.0, 2e9)
            .unwrap();
        close(second, 2.0 * first);
        close(state.volume, 3.0 * first);
        assert!(matches!(
            law.advance(&synthetic_dry(), &mut state, 1.0, 1.0, 0.0),
            Err(TriboError::InvalidInput { field: "hardness" })
        ));
    }
}
