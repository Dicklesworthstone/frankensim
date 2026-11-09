//! Temperature-dependent parameters of the original dry-friction ladder.
//!
//! A source supplies parameters of an explicitly selected executable curve.
//! Local sliding speed stays an argument of that curve at every Hertz pressure
//! quadrature point; it is never replaced by a centre or effective patch speed.
//! This binding adds no constitutive history, calibration, thermal storage or
//! source authority. A parameter's source support does not by itself validate
//! the chosen curve over an unmeasured range of speed or contact geometry.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MatDbError, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::sliding_patch::rate_dependent::HertzFrictionPatch;
use fs_tribo::{FrictionLaw, InputAuthority, TriboError};

use super::{BoundDryInterfaceSystem, SourceThermalCoulombError, ThermalFrictionAxes};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const INVERSE_SPEED: Dims = Dims([-1, 0, 1, 0, 0, 0]);

/// Explicit parameter roles for the existing [`FrictionLaw`] formulas.
///
/// Coefficients and exponents are dimensionless; characteristic speeds are
/// m/s, and slopes are s/m. Equal requirements may deliberately share one
/// source. Different requirements cannot alias a property name. All values
/// are resolved together before the original law validates its physical domain.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalFrictionCurve {
    Coulomb {
        static_mu: ScalarPropertyRequirement,
        kinetic_mu: ScalarPropertyRequirement,
    },
    VelocityDependent {
        static_mu: ScalarPropertyRequirement,
        mu_zero: ScalarPropertyRequirement,
        slope_per_speed: ScalarPropertyRequirement,
    },
    Stribeck {
        static_mu: ScalarPropertyRequirement,
        kinetic_mu: ScalarPropertyRequirement,
        characteristic_speed: ScalarPropertyRequirement,
        viscous_per_speed: ScalarPropertyRequirement,
    },
    Hyperbolic {
        static_mu: ScalarPropertyRequirement,
        kinetic_mu: ScalarPropertyRequirement,
        characteristic_speed: ScalarPropertyRequirement,
    },
    ShiftedPower {
        static_mu: ScalarPropertyRequirement,
        kinetic_mu: ScalarPropertyRequirement,
        characteristic_speed: ScalarPropertyRequirement,
        exponent: ScalarPropertyRequirement,
    },
}

impl ThermalFrictionCurve {
    fn roles(&self) -> (u8, Vec<(&'static str, &ScalarPropertyRequirement, Dims)>) {
        use ThermalFrictionCurve::*;
        match self {
            Coulomb {
                static_mu,
                kinetic_mu,
            } => (
                0,
                vec![
                    ("static_mu", static_mu, Dims::NONE),
                    ("kinetic_mu", kinetic_mu, Dims::NONE),
                ],
            ),
            VelocityDependent {
                static_mu,
                mu_zero,
                slope_per_speed,
            } => (
                1,
                vec![
                    ("static_mu", static_mu, Dims::NONE),
                    ("mu_zero", mu_zero, Dims::NONE),
                    ("slope_per_speed", slope_per_speed, INVERSE_SPEED),
                ],
            ),
            Stribeck {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                viscous_per_speed,
            } => (
                2,
                vec![
                    ("static_mu", static_mu, Dims::NONE),
                    ("kinetic_mu", kinetic_mu, Dims::NONE),
                    ("characteristic_speed", characteristic_speed, SPEED),
                    ("viscous_per_speed", viscous_per_speed, INVERSE_SPEED),
                ],
            ),
            Hyperbolic {
                static_mu,
                kinetic_mu,
                characteristic_speed,
            } => (
                3,
                vec![
                    ("static_mu", static_mu, Dims::NONE),
                    ("kinetic_mu", kinetic_mu, Dims::NONE),
                    ("characteristic_speed", characteristic_speed, SPEED),
                ],
            ),
            ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                exponent,
            } => (
                4,
                vec![
                    ("static_mu", static_mu, Dims::NONE),
                    ("kinetic_mu", kinetic_mu, Dims::NONE),
                    ("characteristic_speed", characteristic_speed, SPEED),
                    ("exponent", exponent, Dims::NONE),
                ],
            ),
        }
    }

    fn law(
        &self,
        source: &ResolvedInterfaceStatePoint,
    ) -> Result<FrictionLaw, SourceHertzFrictionError> {
        let value = |requirement: &ScalarPropertyRequirement| {
            source
                .property(requirement.name())
                .map(|p| p.value_si())
                .ok_or(SourceThermalCoulombError::Invalid(
                    "selected friction parameter is missing or ambiguous",
                ))
        };
        use ThermalFrictionCurve::*;
        let law = match self {
            Coulomb {
                static_mu,
                kinetic_mu,
            } => FrictionLaw::Coulomb {
                static_mu: value(static_mu)?,
                kinetic_mu: value(kinetic_mu)?,
            },
            VelocityDependent {
                static_mu,
                mu_zero,
                slope_per_speed,
            } => FrictionLaw::VelocityDependent {
                static_mu: value(static_mu)?,
                mu_zero: value(mu_zero)?,
                slope_per_speed: value(slope_per_speed)?,
            },
            Stribeck {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                viscous_per_speed,
            } => FrictionLaw::Stribeck {
                static_mu: value(static_mu)?,
                kinetic_mu: value(kinetic_mu)?,
                characteristic_speed: value(characteristic_speed)?,
                viscous_per_speed: value(viscous_per_speed)?,
            },
            Hyperbolic {
                static_mu,
                kinetic_mu,
                characteristic_speed,
            } => FrictionLaw::Hyperbolic {
                static_mu: value(static_mu)?,
                kinetic_mu: value(kinetic_mu)?,
                characteristic_speed: value(characteristic_speed)?,
            },
            ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                exponent,
            } => FrictionLaw::ShiftedPower {
                static_mu: value(static_mu)?,
                kinetic_mu: value(kinetic_mu)?,
                characteristic_speed: value(characteristic_speed)?,
                exponent: value(exponent)?,
            },
        };
        law.validate()?;
        Ok(law)
    }
}

/// Immutable original interface, selected parameter claims and curve roles.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceThermalFriction {
    interface: BoundDryInterfaceSystem,
    initial_source: ResolvedInterfaceStatePoint,
    curve: ThermalFrictionCurve,
    axes: ThermalFrictionAxes,
    identity: ContentHash,
}

/// One complete thermal query and its original local-speed Hertz law.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzFriction {
    owner: ContentHash,
    source: ResolvedInterfaceStatePoint,
    law: HertzFrictionPatch,
}

impl SourceHertzFriction {
    #[must_use]
    pub const fn owner_identity(&self) -> ContentHash {
        self.owner
    }
    #[must_use]
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }
    #[must_use]
    pub const fn law(&self) -> &HertzFrictionPatch {
        &self.law
    }
}

/// Original source/query or constitutive/geometry refusal; no state is published.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceHertzFrictionError {
    Source(SourceThermalCoulombError),
    Tribo(TriboError),
}
impl From<SourceThermalCoulombError> for SourceHertzFrictionError {
    fn from(e: SourceThermalCoulombError) -> Self {
        Self::Source(e)
    }
}
impl From<MaterialStatePointError> for SourceHertzFrictionError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Source(e.into())
    }
}
impl From<MatDbError> for SourceHertzFrictionError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e.into())
    }
}
impl From<TriboError> for SourceHertzFrictionError {
    fn from(e: TriboError) -> Self {
        Self::Tribo(e)
    }
}
impl fmt::Display for SourceHertzFrictionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "source-backed thermal friction curve: {self:?}")
    }
}
impl std::error::Error for SourceHertzFrictionError {}

impl SourceThermalFriction {
    /// Select all parameter claims once on one ordered dry interface.
    ///
    /// Every response coordinate and applicability bound must be a static
    /// external condition. Velocity/history-dependent *parameter claims* refuse:
    /// local speed is owned by the chosen executable curve, not a frozen source
    /// coordinate. No claim is selected again during temperature feedback.
    pub fn resolve(
        interface: BoundDryInterfaceSystem,
        initial_query: QueryPoint,
        curve: ThermalFrictionCurve,
        selection: MaterialPropertySelection,
        axes: ThermalFrictionAxes,
    ) -> Result<Self, SourceHertzFrictionError> {
        super::validate_axes(&initial_query, &axes, &[])?;
        let (family, roles) = curve.roles();
        let mut requirements: Vec<ScalarPropertyRequirement> = Vec::with_capacity(roles.len());
        for (_, requirement, dimensions) in &roles {
            if requirement.quantity() != QuantitySpec::dimensional(*dimensions) {
                return Err(SourceThermalCoulombError::Invalid(
                    "friction parameter quantity does not match its curve role",
                )
                .into());
            }
            if let Some(previous) = requirements.iter().find(|r| r.name() == requirement.name()) {
                if previous != *requirement {
                    return Err(SourceThermalCoulombError::Invalid(
                        "distinct friction parameter requirements cannot alias a property name",
                    )
                    .into());
                }
            } else {
                requirements.push((*requirement).clone());
            }
        }
        let initial_source = resolve_interface_state_point(
            interface.card(),
            &initial_query,
            &requirements,
            selection,
        )?;
        if initial_source
            .properties()
            .iter()
            .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(SourceThermalCoulombError::Invalid(
                "friction curve parameters require static sources at fixed external conditions",
            )
            .into());
        }
        curve.law(&initial_source)?;
        let mut hash = DomainHasher::new("org.frankensim.source-thermal-friction-curve.v1");
        hash.update(initial_source.identity().as_bytes());
        hash.update(&[
            family,
            match interface.interface().provenance().authority() {
                InputAuthority::CallerDeclared => 0,
                InputAuthority::SyntheticFixture => 1,
                InputAuthority::Estimated => 2,
            },
        ]);
        hash.update(&[u8::from(axes.nominal_pressure.is_some())]);
        for text in [
            axes.temperature.as_str(),
            axes.nominal_pressure.as_deref().unwrap_or(""),
        ] {
            hash.update(&(text.len() as u64).to_le_bytes());
            hash.update(text.as_bytes());
        }
        for (role, requirement, _) in roles {
            for text in [role, requirement.name()] {
                hash.update(&(text.len() as u64).to_le_bytes());
                hash.update(text.as_bytes());
            }
        }
        Ok(Self {
            interface,
            initial_source,
            curve,
            axes,
            identity: hash.finalize(),
        })
    }

    #[must_use]
    pub const fn interface(&self) -> &BoundDryInterfaceSystem {
        &self.interface
    }
    #[must_use]
    pub const fn initial_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.initial_source
    }
    #[must_use]
    pub const fn curve(&self) -> &ThermalFrictionCurve {
        &self.curve
    }
    #[must_use]
    pub const fn axes(&self) -> &ThermalFrictionAxes {
        &self.axes
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Refresh the original complete parameter bundle at absolute T [K] and the
    /// explicitly owned nominal pressure [Pa]. All other coordinates stay fixed.
    pub fn coefficients_at(
        &self,
        temperature_k: f64,
        nominal_pressure_pa: Option<f64>,
    ) -> Result<(ResolvedInterfaceStatePoint, FrictionLaw), SourceHertzFrictionError> {
        let initial = self.initial_source.query();
        let mut point = initial
            .clone()
            .with(&self.axes.temperature, temperature_k)?;
        if let (Some(axis), Some(p)) = (&self.axes.nominal_pressure, nominal_pressure_pa) {
            point = point.with(axis, p)?;
        }
        super::query::validate_refresh(
            &initial,
            &point,
            &self.axes,
            &[],
            Some(temperature_k),
            nominal_pressure_pa,
        )?;
        let source = if point == initial {
            self.initial_source.clone()
        } else {
            self.initial_source
                .at_query_point(self.interface.card(), &point)?
        };
        let law = self.curve.law(&source)?;
        Ok((source, law))
    }

    /// Bind the complete original speed law to the caller's admitted Hertz
    /// normal force [N] and radius [m]. Geometry is not inferred from source
    /// pressure; its applicability stays with the original normal-contact owner.
    pub fn hertz_at(
        &self,
        temperature_k: f64,
        nominal_pressure_pa: Option<f64>,
        normal_force_n: f64,
        contact_radius_m: f64,
    ) -> Result<SourceHertzFriction, SourceHertzFrictionError> {
        let (source, law) = self.coefficients_at(temperature_k, nominal_pressure_pa)?;
        let law = HertzFrictionPatch::new(normal_force_n, contact_radius_m, law)?;
        Ok(SourceHertzFriction {
            owner: self.identity,
            source,
            law,
        })
    }
}
