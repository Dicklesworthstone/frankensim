#![forbid(unsafe_code)]

//! Typed L1 thermochemical law data and deterministic standard-state evaluation.
//!
//! These first slices deliberately reuse `fs-qty`'s exact chemistry artifacts
//! instead of creating a second species or conservation system. They add a
//! provenance-bound NASA-9 standard-state evaluator and a bounded,
//! frozen-composition ideal-gas mixture evaluator, the positive-state
//! mechanical ideal-gas `p`-`rho`-`T` closure, and a bounded Wilke dilute-gas
//! viscosity mixing closure, plus bounded ideal-gas standard-state reaction
//! equilibrium. Derived identities remain scoped to explicit ideal-gas
//! conventions. This crate owns no transport solver, evolving state,
//! composition-equilibrium solve, kinetics integrator, or L3 protocol.

pub mod eos;
pub mod equilibrium;
pub mod mixture;
pub mod transport;

use core::fmt;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError};
use fs_qty::{Dims, MolarMass, Pressure, Qty, Temperature};

pub use fs_qty::chemistry::{
    ChargeVector, ChemistryError, ConservationCertificate, ElementId, ElementalMatrix,
    MassAmountBasis, ReactionId, SpeciesId, StoichiometricMatrix, verify_conservation,
};
pub use fs_qty::semantic::{Composition, CompositionBasis, SemanticError};

/// Crate version retained in diagnostic and evaluation provenance.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Version of the fixed NASA-9 operation tree and receipt layout.
pub const NASA9_EVALUATOR_VERSION_V1: u32 = 1;
/// Exact `fs-matdb` law id required for one NASA-9 region card.
pub const NASA9_LAW_ID_V1: &str = "nasa9-standard-state";
/// Semantic version of the admitted NASA-9 card parameter layout.
pub const NASA9_LAW_VERSION_V1: u32 = 1;
/// Stateless cards use state-schema version zero.
pub const NASA9_STATE_SCHEMA_VERSION_V1: u32 = 0;
/// Hard bound on regions retained by one species model.
pub const MAX_NASA9_REGIONS_V1: usize = 16;
/// Universal molar gas constant in J mol^-1 K^-1.
///
/// The value follows the exact post-2019 SI definitions of the Boltzmann and
/// Avogadro constants.
pub const UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314_462_618_153_24;

/// Coherent-SI molar energy quantity, J/mol.
pub type MolarEnergyQuantityV1 = Qty<2, 1, -2, 0, 0, -1>;
/// Coherent-SI molar thermal quantity, J/(mol K).
pub type MolarThermalQuantityV1 = Qty<2, 1, -2, -1, 0, -1>;
/// Coherent-SI mass-specific thermal quantity, J/(kg K).
pub type MassSpecificThermalQuantityV1 = Qty<2, 0, -2, -1, 0, 0>;

macro_rules! molar_quantity {
    ($(#[$meta:meta])* $name:ident, $quantity:ty) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
        pub struct $name($quantity);

        impl $name {
            const fn new(value: f64) -> Self {
                Self(<$quantity>::new(value))
            }

            /// Raw coherent-SI scalar.
            #[must_use]
            pub const fn value(self) -> f64 {
                self.0.value()
            }

            /// Dimensioned coherent-SI quantity.
            #[must_use]
            pub const fn quantity(self) -> $quantity {
                self.0
            }
        }
    };
}

molar_quantity!(
    /// Standard-state molar heat capacity at constant pressure.
    MolarHeatCapacityV1,
    MolarThermalQuantityV1
);
molar_quantity!(
    /// Standard-state molar enthalpy.
    MolarEnthalpyV1,
    MolarEnergyQuantityV1
);
molar_quantity!(
    /// Standard-state molar entropy.
    MolarEntropyV1,
    MolarThermalQuantityV1
);
molar_quantity!(
    /// Standard-state molar internal energy derived under the retained EOS.
    MolarInternalEnergyV1,
    MolarEnergyQuantityV1
);
molar_quantity!(
    /// Standard-state molar Gibbs energy derived as `h - T s`.
    MolarGibbsEnergyV1,
    MolarEnergyQuantityV1
);

/// Canonical mass-specific gas constant `R / M` shared by mixtures and EOS
/// closures.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct MassSpecificGasConstantV1(MassSpecificThermalQuantityV1);

impl MassSpecificGasConstantV1 {
    pub(crate) const fn new(value: f64) -> Self {
        Self(MassSpecificThermalQuantityV1::new(value))
    }

    /// Raw coherent-SI scalar.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensioned coherent-SI quantity.
    #[must_use]
    pub const fn quantity(self) -> MassSpecificThermalQuantityV1 {
        self.0
    }
}

/// Explicit standard-state phase supported by the first evaluator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StandardStatePhaseV1 {
    /// Caller-declared gas-phase standard state.
    Gas,
}

/// Explicit reference EOS used for derived standard-state properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReferenceEquationOfStateV1 {
    /// `p v = R T` on a molar basis, used only to derive `u = h - R T`.
    IdealGas,
}

/// Validated opaque identity of the elemental reference convention.
///
/// The id names a convention; it does not authenticate the convention or
/// establish formation-property authority.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ElementalReferenceIdV1(String);

impl ElementalReferenceIdV1 {
    /// Admit a compact canonical ASCII convention id.
    ///
    /// # Errors
    /// Refuses empty, overlong, non-alphanumeric-leading, whitespace-bearing,
    /// or non-ASCII ids.
    pub fn new(value: impl Into<String>) -> Result<Self, ThermochemErrorV1> {
        let value = value.into();
        if value.is_empty() {
            return Err(ThermochemErrorV1::InvalidElementalReference {
                value,
                reason: "id must not be empty",
            });
        }
        if value.len() > 128 {
            return Err(ThermochemErrorV1::InvalidElementalReference {
                value,
                reason: "id exceeds the 128-byte limit",
            });
        }
        if !value.as_bytes()[0].is_ascii_alphanumeric() {
            return Err(ThermochemErrorV1::InvalidElementalReference {
                value,
                reason: "id must begin with an ASCII letter or digit",
            });
        }
        if !value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        }) {
            return Err(ThermochemErrorV1::InvalidElementalReference {
                value,
                reason: "id must use compact ASCII without whitespace",
            });
        }
        Ok(Self(value))
    }

    /// Canonical id text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for ElementalReferenceIdV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// Complete convention required before standard-state evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct StandardStateConventionV1 {
    phase: StandardStatePhaseV1,
    eos: ReferenceEquationOfStateV1,
    reference_pressure: Pressure,
    elemental_reference: ElementalReferenceIdV1,
}

impl StandardStateConventionV1 {
    /// Construct the explicit phase/EOS/reference convention.
    ///
    /// # Errors
    /// Refuses a non-positive or non-finite reference pressure.
    pub fn new(
        phase: StandardStatePhaseV1,
        eos: ReferenceEquationOfStateV1,
        reference_pressure: Pressure,
        elemental_reference: ElementalReferenceIdV1,
    ) -> Result<Self, ThermochemErrorV1> {
        let value = reference_pressure.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(ThermochemErrorV1::InvalidReferencePressure {
                bits: value.to_bits(),
            });
        }
        Ok(Self {
            phase,
            eos,
            reference_pressure,
            elemental_reference,
        })
    }

    /// Retained standard-state phase.
    #[must_use]
    pub const fn phase(&self) -> StandardStatePhaseV1 {
        self.phase
    }

    /// Retained reference EOS.
    #[must_use]
    pub const fn eos(&self) -> ReferenceEquationOfStateV1 {
        self.eos
    }

    /// Positive finite reference pressure.
    #[must_use]
    pub const fn reference_pressure(&self) -> Pressure {
        self.reference_pressure
    }

    /// Opaque elemental reference convention.
    #[must_use]
    pub const fn elemental_reference(&self) -> &ElementalReferenceIdV1 {
        &self.elemental_reference
    }
}

/// NASA-9 property named by a non-finite evaluation refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nasa9PropertyV1 {
    /// Molar heat capacity at constant pressure.
    HeatCapacity,
    /// Molar enthalpy.
    Enthalpy,
    /// Molar entropy.
    Entropy,
    /// Ideal-gas-derived molar internal energy.
    InternalEnergy,
    /// Derived molar Gibbs energy.
    GibbsEnergy,
}

/// Typed construction or evaluation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermochemErrorV1 {
    /// The supplied immutable material-data card failed its own gates.
    MaterialCard(MatDbError),
    /// Card law id is not the pinned NASA-9 law.
    CardLawMismatch {
        /// Offered law id.
        found: String,
    },
    /// Card law version is not supported.
    CardLawVersionMismatch {
        /// Offered version.
        found: u32,
        /// Supported version.
        supported: u32,
    },
    /// A NASA-9 card declared state despite the evaluator being stateless.
    CardStateConventionMismatch {
        /// Offered state-schema version.
        state_schema_version: u32,
        /// Offered initial-state policy.
        initial_state: InitialStatePolicy,
    },
    /// Required coefficient or reference-pressure parameter is absent.
    MissingCardParameter {
        /// Canonical parameter name.
        parameter: &'static str,
    },
    /// A foreign parameter would make the card layout ambiguous.
    UnexpectedCardParameter {
        /// Foreign parameter name.
        parameter: String,
    },
    /// A parameter has the wrong six-base dimension vector.
    CardParameterDimsMismatch {
        /// Parameter name.
        parameter: &'static str,
        /// Required dimensions.
        expected: Dims,
        /// Offered dimensions.
        found: Dims,
    },
    /// Card validity is not exactly one positive finite temperature range.
    InvalidTemperatureValidity {
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// Opaque elemental-reference id is not canonical.
    InvalidElementalReference {
        /// Rejected text.
        value: String,
        /// Stable guidance.
        reason: &'static str,
    },
    /// Molar mass must be finite and strictly positive.
    InvalidMolarMass {
        /// Exact rejected IEEE bits.
        bits: u64,
    },
    /// Reference pressure must be finite and strictly positive.
    InvalidReferencePressure {
        /// Exact rejected IEEE bits.
        bits: u64,
    },
    /// A model requires at least one temperature region.
    EmptyRegions,
    /// Region count exceeded the fixed metadata bound.
    TooManyRegions {
        /// Offered region count.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// Region cards overlap or arrive out of increasing-temperature order.
    OverlappingRegions {
        /// Earlier region index.
        lower: usize,
        /// Later region index.
        upper: usize,
    },
    /// Region card reference pressure disagrees with the model convention.
    RegionReferencePressureMismatch {
        /// Region index.
        region: usize,
        /// Exact convention bits.
        expected_bits: u64,
        /// Exact card bits.
        found_bits: u64,
    },
    /// Evaluation temperature is zero, negative, NaN, or infinite.
    InvalidEvaluationTemperature {
        /// Exact rejected IEEE bits.
        bits: u64,
    },
    /// Positive finite temperature lies outside every retained region.
    TemperatureOutOfRange {
        /// Exact requested IEEE bits.
        bits: u64,
    },
    /// Fixed-order arithmetic produced a non-finite output.
    NonFiniteEvaluation {
        /// First derived property that failed.
        property: Nasa9PropertyV1,
        /// Exact non-finite IEEE bits.
        bits: u64,
    },
}

impl fmt::Display for ThermochemErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaterialCard(error) => write!(f, "NASA-9 source card refused: {error}"),
            Self::CardLawMismatch { found } => write!(
                f,
                "NASA-9 region requires law id {NASA9_LAW_ID_V1:?}, found {found:?}"
            ),
            Self::CardLawVersionMismatch { found, supported } => write!(
                f,
                "NASA-9 law version {found} is unsupported; expected {supported}"
            ),
            Self::CardStateConventionMismatch {
                state_schema_version,
                initial_state,
            } => write!(
                f,
                "NASA-9 is stateless but card declared schema {state_schema_version} and {initial_state:?}"
            ),
            Self::MissingCardParameter { parameter } => {
                write!(f, "NASA-9 card is missing parameter {parameter:?}")
            }
            Self::UnexpectedCardParameter { parameter } => {
                write!(f, "NASA-9 card has unexpected parameter {parameter:?}")
            }
            Self::CardParameterDimsMismatch {
                parameter,
                expected,
                found,
            } => write!(
                f,
                "NASA-9 parameter {parameter:?} has dimensions {found:?}, expected {expected:?}"
            ),
            Self::InvalidTemperatureValidity { reason } => {
                write!(f, "NASA-9 temperature validity is invalid: {reason}")
            }
            Self::InvalidElementalReference { value, reason } => {
                write!(f, "elemental reference id {value:?} refused: {reason}")
            }
            Self::InvalidMolarMass { bits } => {
                write!(
                    f,
                    "molar mass must be positive and finite (bits {bits:#018x})"
                )
            }
            Self::InvalidReferencePressure { bits } => write!(
                f,
                "reference pressure must be positive and finite (bits {bits:#018x})"
            ),
            Self::EmptyRegions => f.write_str("NASA-9 model requires at least one region"),
            Self::TooManyRegions { offered, limit } => write!(
                f,
                "NASA-9 model offered {offered} regions, exceeding limit {limit}"
            ),
            Self::OverlappingRegions { lower, upper } => write!(
                f,
                "NASA-9 regions {lower} and {upper} overlap or are out of order"
            ),
            Self::RegionReferencePressureMismatch {
                region,
                expected_bits,
                found_bits,
            } => write!(
                f,
                "NASA-9 region {region} reference pressure bits {found_bits:#018x} disagree with convention {expected_bits:#018x}"
            ),
            Self::InvalidEvaluationTemperature { bits } => write!(
                f,
                "evaluation temperature must be positive and finite (bits {bits:#018x})"
            ),
            Self::TemperatureOutOfRange { bits } => write!(
                f,
                "temperature bits {bits:#018x} are outside every NASA-9 region"
            ),
            Self::NonFiniteEvaluation { property, bits } => write!(
                f,
                "NASA-9 {property:?} evaluation became non-finite (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for ThermochemErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::MaterialCard(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MatDbError> for ThermochemErrorV1 {
    fn from(value: MatDbError) -> Self {
        Self::MaterialCard(value)
    }
}

/// One validated, provenance-bound NASA-9 temperature region.
#[derive(Debug, Clone, PartialEq)]
pub struct Nasa9RegionV1 {
    temperature_min: Temperature,
    temperature_max: Temperature,
    coefficients: [f64; 9],
    reference_pressure: Pressure,
    card_identity: [u8; 32],
    card: ConstitutiveModelCard,
}

impl Nasa9RegionV1 {
    /// Validate and retain one immutable `fs-matdb` NASA-9 card.
    ///
    /// The exact parameter set is `a0` through `a8` plus
    /// `reference_pressure`. Coefficient dimensions follow the NASA-9
    /// temperature powers, and validity must contain exactly one `T` axis.
    ///
    /// # Errors
    /// Returns a typed card, law, state, parameter, dimension, or validity
    /// refusal. Nothing partial escapes.
    pub fn from_card(card: ConstitutiveModelCard) -> Result<Self, ThermochemErrorV1> {
        card.validate()?;
        if card.law.0 != NASA9_LAW_ID_V1 {
            return Err(ThermochemErrorV1::CardLawMismatch {
                found: card.law.0.clone(),
            });
        }
        if card.law_version != NASA9_LAW_VERSION_V1 {
            return Err(ThermochemErrorV1::CardLawVersionMismatch {
                found: card.law_version,
                supported: NASA9_LAW_VERSION_V1,
            });
        }
        if card.state_schema_version != NASA9_STATE_SCHEMA_VERSION_V1
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(ThermochemErrorV1::CardStateConventionMismatch {
                state_schema_version: card.state_schema_version,
                initial_state: card.initial_state,
            });
        }

        let required = required_nasa9_parameter_dims();
        for (name, expected) in required {
            let parameter = card
                .parameters
                .get(name)
                .ok_or(ThermochemErrorV1::MissingCardParameter { parameter: name })?;
            if parameter.dims != expected {
                return Err(ThermochemErrorV1::CardParameterDimsMismatch {
                    parameter: name,
                    expected,
                    found: parameter.dims,
                });
            }
        }
        if let Some(parameter) = card
            .parameters
            .keys()
            .find(|name| !is_required_nasa9_parameter(name))
        {
            return Err(ThermochemErrorV1::UnexpectedCardParameter {
                parameter: parameter.clone(),
            });
        }

        if card.validity.bounds().len() != 1 {
            return Err(ThermochemErrorV1::InvalidTemperatureValidity {
                reason: "card validity must contain exactly the T axis",
            });
        }
        let (temperature_min, temperature_max) =
            card.validity
                .bound("T")
                .ok_or(ThermochemErrorV1::InvalidTemperatureValidity {
                    reason: "card validity is missing the T axis",
                })?;
        if !temperature_min.is_finite()
            || !temperature_max.is_finite()
            || temperature_min <= 0.0
            || temperature_min >= temperature_max
        {
            return Err(ThermochemErrorV1::InvalidTemperatureValidity {
                reason: "T bounds must be positive, finite, and strictly increasing",
            });
        }

        let mut coefficients = [0.0; 9];
        for (index, coefficient) in coefficients.iter_mut().enumerate() {
            let name = NASA9_COEFFICIENT_NAMES[index];
            *coefficient = card
                .parameters
                .get(name)
                .ok_or(ThermochemErrorV1::MissingCardParameter { parameter: name })?
                .value;
        }
        let reference_pressure = Pressure::new(
            card.parameters
                .get("reference_pressure")
                .ok_or(ThermochemErrorV1::MissingCardParameter {
                    parameter: "reference_pressure",
                })?
                .value,
        );
        if reference_pressure.value() <= 0.0 {
            return Err(ThermochemErrorV1::InvalidReferencePressure {
                bits: reference_pressure.value().to_bits(),
            });
        }
        let card_identity = *card.content_hash().as_bytes();
        Ok(Self {
            temperature_min: Temperature::new(temperature_min),
            temperature_max: Temperature::new(temperature_max),
            coefficients,
            reference_pressure,
            card_identity,
            card,
        })
    }

    /// Inclusive lower temperature bound in kelvin.
    #[must_use]
    pub const fn temperature_min(&self) -> Temperature {
        self.temperature_min
    }

    /// Upper temperature bound in kelvin.
    #[must_use]
    pub const fn temperature_max(&self) -> Temperature {
        self.temperature_max
    }

    /// Exact nine coefficient scalars in `a0..a8` order.
    #[must_use]
    pub const fn coefficients(&self) -> &[f64; 9] {
        &self.coefficients
    }

    /// Reference pressure retained by this source card.
    #[must_use]
    pub const fn reference_pressure(&self) -> Pressure {
        self.reference_pressure
    }

    /// Immutable source-card content identity.
    #[must_use]
    pub const fn card_identity(&self) -> &[u8; 32] {
        &self.card_identity
    }

    /// Full immutable source card for provenance replay.
    #[must_use]
    pub const fn card(&self) -> &ConstitutiveModelCard {
        &self.card
    }
}

const NASA9_COEFFICIENT_NAMES: [&str; 9] = ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8"];

const fn required_nasa9_parameter_dims() -> [(&'static str, Dims); 10] {
    [
        ("a0", Dims([0, 0, 0, 2, 0, 0])),
        ("a1", Dims([0, 0, 0, 1, 0, 0])),
        ("a2", Dims::NONE),
        ("a3", Dims([0, 0, 0, -1, 0, 0])),
        ("a4", Dims([0, 0, 0, -2, 0, 0])),
        ("a5", Dims([0, 0, 0, -3, 0, 0])),
        ("a6", Dims([0, 0, 0, -4, 0, 0])),
        ("a7", Dims([0, 0, 0, 1, 0, 0])),
        ("a8", Dims::NONE),
        ("reference_pressure", Pressure::DIMS),
    ]
}

fn is_required_nasa9_parameter(name: &str) -> bool {
    NASA9_COEFFICIENT_NAMES.contains(&name) || name == "reference_pressure"
}

/// Validated multi-region NASA-9 model for one species.
#[derive(Debug, Clone, PartialEq)]
pub struct Nasa9StandardStateModelV1 {
    species: SpeciesId,
    molar_mass: MolarMass,
    convention: StandardStateConventionV1,
    regions: Vec<Nasa9RegionV1>,
}

impl Nasa9StandardStateModelV1 {
    /// Bind one species, molar mass, convention, and ordered source cards.
    ///
    /// The current `fs-matdb` law-card schema does not carry species, molar
    /// mass, phase, EOS, or elemental-reference identity. Those inputs are
    /// therefore caller-declared associations: the evaluation receipt binds
    /// them exactly but does not authenticate that they describe the retained
    /// source cards. In particular, callers must not label condensed-phase
    /// coefficients as gas data merely to obtain the ideal-gas `u` derivation.
    ///
    /// Regions may have gaps but may not overlap. At an exactly shared
    /// boundary the upper region is selected. Every card must carry the exact
    /// convention reference pressure.
    ///
    /// # Errors
    /// Refuses invalid molar mass, empty/excess/overlapping regions, or any
    /// reference-pressure drift.
    pub fn new(
        species: SpeciesId,
        molar_mass: MolarMass,
        convention: StandardStateConventionV1,
        regions: Vec<Nasa9RegionV1>,
    ) -> Result<Self, ThermochemErrorV1> {
        let molar_mass_value = molar_mass.value();
        if !molar_mass_value.is_finite() || molar_mass_value <= 0.0 {
            return Err(ThermochemErrorV1::InvalidMolarMass {
                bits: molar_mass_value.to_bits(),
            });
        }
        if regions.is_empty() {
            return Err(ThermochemErrorV1::EmptyRegions);
        }
        if regions.len() > MAX_NASA9_REGIONS_V1 {
            return Err(ThermochemErrorV1::TooManyRegions {
                offered: regions.len(),
                limit: MAX_NASA9_REGIONS_V1,
            });
        }
        for (index, region) in regions.iter().enumerate() {
            if region.reference_pressure().value().to_bits()
                != convention.reference_pressure().value().to_bits()
            {
                return Err(ThermochemErrorV1::RegionReferencePressureMismatch {
                    region: index,
                    expected_bits: convention.reference_pressure().value().to_bits(),
                    found_bits: region.reference_pressure().value().to_bits(),
                });
            }
        }
        for (lower, pair) in regions.windows(2).enumerate() {
            if pair[0].temperature_max().value() > pair[1].temperature_min().value() {
                return Err(ThermochemErrorV1::OverlappingRegions {
                    lower,
                    upper: lower + 1,
                });
            }
        }
        Ok(Self {
            species,
            molar_mass,
            convention,
            regions,
        })
    }

    /// Exact species id.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Positive finite molar mass.
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }

    /// Explicit standard-state convention.
    #[must_use]
    pub const fn convention(&self) -> &StandardStateConventionV1 {
        &self.convention
    }

    /// Caller-significant ordered temperature regions.
    #[must_use]
    pub fn regions(&self) -> &[Nasa9RegionV1] {
        &self.regions
    }

    /// Evaluate standard-state molar properties and retain exact provenance.
    ///
    /// The operation tree is the NASA/TP-2002-211556 nine-coefficient form.
    /// `fs_math::det::ln` avoids platform libm drift. `u` is derived only
    /// because this v1 convention is explicitly ideal gas; `g` is derived as
    /// `h - T s` under the same retained standard state.
    ///
    /// # Errors
    /// Refuses invalid/out-of-range temperature or the first non-finite
    /// computed property. No partial evaluation escapes.
    pub fn evaluate(
        &self,
        temperature: Temperature,
    ) -> Result<Nasa9EvaluationV1, ThermochemErrorV1> {
        let t = temperature.value();
        if !t.is_finite() || t <= 0.0 {
            return Err(ThermochemErrorV1::InvalidEvaluationTemperature { bits: t.to_bits() });
        }
        let (region_index, region) = self
            .region_for_temperature(t)
            .ok_or(ThermochemErrorV1::TemperatureOutOfRange { bits: t.to_bits() })?;
        let [a0, a1, a2, a3, a4, a5, a6, a7, a8] = region.coefficients;
        let inverse_t = 1.0 / t;
        let inverse_t2 = inverse_t * inverse_t;
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t2 * t2;
        // `t` is the coherent-SI kelvin scalar, so this is ln(T / 1 K).
        let ln_t = fs_math::det::ln(t);

        let cp_over_r =
            a0 * inverse_t2 + a1 * inverse_t + a2 + a3 * t + a4 * t2 + a5 * t3 + a6 * t4;
        let h_over_rt = -a0 * inverse_t2
            + a1 * ln_t * inverse_t
            + a2
            + (a3 / 2.0) * t
            + (a4 / 3.0) * t2
            + (a5 / 4.0) * t3
            + (a6 / 5.0) * t4
            + a7 * inverse_t;
        let s_over_r = -(a0 / 2.0) * inverse_t2 - a1 * inverse_t
            + a2 * ln_t
            + a3 * t
            + (a4 / 2.0) * t2
            + (a5 / 3.0) * t3
            + (a6 / 4.0) * t4
            + a8;

        let cp = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * cp_over_r;
        let h = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t * h_over_rt;
        let s = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * s_over_r;
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t;
        let u = h - rt;
        let g = h - t * s;
        for (property, value) in [
            (Nasa9PropertyV1::HeatCapacity, cp),
            (Nasa9PropertyV1::Enthalpy, h),
            (Nasa9PropertyV1::Entropy, s),
            (Nasa9PropertyV1::InternalEnergy, u),
            (Nasa9PropertyV1::GibbsEnergy, g),
        ] {
            if !value.is_finite() {
                return Err(ThermochemErrorV1::NonFiniteEvaluation {
                    property,
                    bits: value.to_bits(),
                });
            }
        }

        let properties = StandardStateMolarPropertiesV1 {
            cp: MolarHeatCapacityV1::new(cp),
            h: MolarEnthalpyV1::new(h),
            s: MolarEntropyV1::new(s),
            u: MolarInternalEnergyV1::new(u),
            g: MolarGibbsEnergyV1::new(g),
        };
        let receipt = Nasa9EvaluationReceiptV1 {
            evaluator_version: NASA9_EVALUATOR_VERSION_V1,
            fs_math_version: fs_math::VERSION,
            gas_constant_bits: UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits(),
            species: self.species.clone(),
            molar_mass_bits: self.molar_mass.value().to_bits(),
            temperature_bits: t.to_bits(),
            reference_pressure_bits: self.convention.reference_pressure().value().to_bits(),
            phase: self.convention.phase(),
            eos: self.convention.eos(),
            elemental_reference: self.convention.elemental_reference().clone(),
            region_index,
            region_min_bits: region.temperature_min().value().to_bits(),
            region_max_bits: region.temperature_max().value().to_bits(),
            coefficient_bits: region.coefficients.map(f64::to_bits),
            card_identity: region.card_identity,
        };
        Ok(Nasa9EvaluationV1 {
            temperature,
            properties,
            receipt,
        })
    }

    fn region_for_temperature(&self, temperature: f64) -> Option<(usize, &Nasa9RegionV1)> {
        for (index, region) in self.regions.iter().enumerate() {
            let min = region.temperature_min().value();
            let max = region.temperature_max().value();
            let shared_with_next = self
                .regions
                .get(index + 1)
                .is_some_and(|next| max.to_bits() == next.temperature_min().value().to_bits());
            if temperature >= min
                && (temperature < max || (!shared_with_next && temperature <= max))
            {
                return Some((index, region));
            }
        }
        None
    }
}

/// Typed standard-state molar properties from one evaluation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StandardStateMolarPropertiesV1 {
    cp: MolarHeatCapacityV1,
    h: MolarEnthalpyV1,
    s: MolarEntropyV1,
    u: MolarInternalEnergyV1,
    g: MolarGibbsEnergyV1,
}

impl StandardStateMolarPropertiesV1 {
    /// Standard-state molar heat capacity at constant pressure.
    #[must_use]
    pub const fn cp(self) -> MolarHeatCapacityV1 {
        self.cp
    }

    /// Standard-state molar enthalpy.
    #[must_use]
    pub const fn h(self) -> MolarEnthalpyV1 {
        self.h
    }

    /// Ideal-gas standard-state molar entropy.
    #[must_use]
    pub const fn s(self) -> MolarEntropyV1 {
        self.s
    }

    /// Molar internal energy derived as `h - R T` under the retained EOS.
    #[must_use]
    pub const fn u(self) -> MolarInternalEnergyV1 {
        self.u
    }

    /// Molar Gibbs energy derived as `h - T s`.
    #[must_use]
    pub const fn g(self) -> MolarGibbsEnergyV1 {
        self.g
    }
}

/// Immutable exact-field evaluation receipt.
///
/// The receipt binds the source card content identity and every convention or
/// scalar needed to reproduce region selection and arithmetic. It is not an
/// evidence color and does not authenticate the source coefficients.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Nasa9EvaluationReceiptV1 {
    evaluator_version: u32,
    fs_math_version: &'static str,
    gas_constant_bits: u64,
    species: SpeciesId,
    molar_mass_bits: u64,
    temperature_bits: u64,
    reference_pressure_bits: u64,
    phase: StandardStatePhaseV1,
    eos: ReferenceEquationOfStateV1,
    elemental_reference: ElementalReferenceIdV1,
    region_index: usize,
    region_min_bits: u64,
    region_max_bits: u64,
    coefficient_bits: [u64; 9],
    card_identity: [u8; 32],
}

impl Nasa9EvaluationReceiptV1 {
    /// Evaluator/operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Deterministic elementary-math crate version.
    #[must_use]
    pub const fn fs_math_version(&self) -> &'static str {
        self.fs_math_version
    }

    /// Exact universal-gas-constant bits used by the operation tree.
    #[must_use]
    pub const fn gas_constant_bits(&self) -> u64 {
        self.gas_constant_bits
    }

    /// Exact species id.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Exact retained molar-mass bits.
    #[must_use]
    pub const fn molar_mass_bits(&self) -> u64 {
        self.molar_mass_bits
    }

    /// Exact evaluation-temperature bits.
    #[must_use]
    pub const fn temperature_bits(&self) -> u64 {
        self.temperature_bits
    }

    /// Exact reference-pressure bits.
    #[must_use]
    pub const fn reference_pressure_bits(&self) -> u64 {
        self.reference_pressure_bits
    }

    /// Explicit phase convention.
    #[must_use]
    pub const fn phase(&self) -> StandardStatePhaseV1 {
        self.phase
    }

    /// Explicit reference EOS.
    #[must_use]
    pub const fn eos(&self) -> ReferenceEquationOfStateV1 {
        self.eos
    }

    /// Explicit elemental-reference convention id.
    #[must_use]
    pub const fn elemental_reference(&self) -> &ElementalReferenceIdV1 {
        &self.elemental_reference
    }

    /// Deterministically selected region index.
    #[must_use]
    pub const fn region_index(&self) -> usize {
        self.region_index
    }

    /// Exact selected lower-bound bits.
    #[must_use]
    pub const fn region_min_bits(&self) -> u64 {
        self.region_min_bits
    }

    /// Exact selected upper-bound bits.
    #[must_use]
    pub const fn region_max_bits(&self) -> u64 {
        self.region_max_bits
    }

    /// Exact selected coefficient bits in `a0..a8` order.
    #[must_use]
    pub const fn coefficient_bits(&self) -> &[u64; 9] {
        &self.coefficient_bits
    }

    /// Exact immutable `fs-matdb` source-card identity.
    #[must_use]
    pub const fn card_identity(&self) -> &[u8; 32] {
        &self.card_identity
    }
}

/// One successful NASA-9 evaluation and its exact-field receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct Nasa9EvaluationV1 {
    temperature: Temperature,
    properties: StandardStateMolarPropertiesV1,
    receipt: Nasa9EvaluationReceiptV1,
}

impl Nasa9EvaluationV1 {
    /// Evaluated absolute temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Derived typed molar properties.
    #[must_use]
    pub const fn properties(&self) -> StandardStateMolarPropertiesV1 {
        self.properties
    }

    /// Exact convention and source receipt.
    #[must_use]
    pub const fn receipt(&self) -> &Nasa9EvaluationReceiptV1 {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{LawId, LawParameter, Provenance};
    use std::collections::BTreeMap;

    fn card(
        temperature_min: f64,
        temperature_max: f64,
        reference_pressure: f64,
        coefficients: [f64; 9],
    ) -> ConstitutiveModelCard {
        let mut parameters = BTreeMap::new();
        for ((name, dims), value) in required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .zip(coefficients)
        {
            parameters.insert(name.to_string(), LawParameter { value, dims });
        }
        parameters.insert(
            "reference_pressure".to_string(),
            LawParameter {
                value: reference_pressure,
                dims: Pressure::DIMS,
            },
        );
        ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.to_string()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", temperature_min, temperature_max),
            sources: Vec::new(),
            provenance: Provenance {
                source: "NASA/TP-2002-211556 synthetic conformance fixture".to_string(),
                license: "test fixture".to_string(),
                artifact: None,
            },
        }
    }

    fn convention(reference_pressure: f64) -> StandardStateConventionV1 {
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(reference_pressure),
            ElementalReferenceIdV1::new("nasa-glenn-2002-reference-elements")
                .expect("valid reference id"),
        )
        .expect("valid convention")
    }

    fn model(coefficients: [f64; 9]) -> Nasa9StandardStateModelV1 {
        Nasa9StandardStateModelV1::new(
            SpeciesId::new("N2").expect("valid species"),
            MolarMass::new(0.028_013_4),
            convention(100_000.0),
            vec![
                Nasa9RegionV1::from_card(card(200.0, 6_000.0, 100_000.0, coefficients))
                    .expect("valid region"),
            ],
        )
        .expect("valid model")
    }

    fn assert_close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        assert!(
            (actual - expected).abs() <= 64.0 * f64::EPSILON * scale,
            "actual {actual:?} differs from expected {expected:?}",
        );
    }

    #[test]
    fn g0_nasa9_parameter_dimensions_are_independently_pinned() {
        assert_eq!(
            required_nasa9_parameter_dims(),
            [
                ("a0", Dims([0, 0, 0, 2, 0, 0])),
                ("a1", Dims([0, 0, 0, 1, 0, 0])),
                ("a2", Dims::NONE),
                ("a3", Dims([0, 0, 0, -1, 0, 0])),
                ("a4", Dims([0, 0, 0, -2, 0, 0])),
                ("a5", Dims([0, 0, 0, -3, 0, 0])),
                ("a6", Dims([0, 0, 0, -4, 0, 0])),
                ("a7", Dims([0, 0, 0, 1, 0, 0])),
                ("a8", Dims::NONE),
                ("reference_pressure", Dims([-1, 1, -2, 0, 0, 0])),
            ],
        );
    }

    #[test]
    fn g0_every_nasa9_coefficient_channel_matches_the_published_operation_tree() {
        let temperature = 2.0;
        let ln_temperature = fs_math::det::ln(temperature);
        let cp_over_r_weights = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 0.0, 0.0];
        let h_over_rt_weights = [
            -0.25,
            0.5 * ln_temperature,
            1.0,
            1.0,
            4.0 / 3.0,
            2.0,
            16.0 / 5.0,
            0.5,
            0.0,
        ];
        let s_over_r_weights = [
            -0.125,
            -0.5,
            ln_temperature,
            2.0,
            2.0,
            8.0 / 3.0,
            4.0,
            0.0,
            1.0,
        ];

        for coefficient_index in 0..9 {
            let mut coefficients = [0.0; 9];
            coefficients[coefficient_index] = 1.0;
            let source_card = card(1.0, 10.0, 100_000.0, coefficients);
            let expected_card_identity = *source_card.content_hash().as_bytes();
            let model = Nasa9StandardStateModelV1::new(
                SpeciesId::new("N2").expect("valid species"),
                MolarMass::new(0.028_013_4),
                convention(100_000.0),
                vec![Nasa9RegionV1::from_card(source_card).expect("valid region")],
            )
            .expect("valid model");
            let evaluation = model
                .evaluate(Temperature::new(temperature))
                .expect("in-range evaluation");
            let properties = evaluation.properties();
            assert_close(
                properties.cp().value() / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
                cp_over_r_weights[coefficient_index],
            );
            assert_close(
                properties.h().value() / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature),
                h_over_rt_weights[coefficient_index],
            );
            assert_close(
                properties.s().value() / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
                s_over_r_weights[coefficient_index],
            );
            assert_eq!(
                evaluation.receipt().coefficient_bits(),
                &coefficients.map(f64::to_bits),
            );
            assert_eq!(
                evaluation.receipt().card_identity(),
                &expected_card_identity
            );
        }
    }

    #[test]
    fn g0_constant_cp_fixture_matches_closed_form_and_derived_potentials() {
        let mut coefficients = [0.0; 9];
        coefficients[2] = 3.5;
        let model = model(coefficients);
        let evaluation = model
            .evaluate(Temperature::new(1_000.0))
            .expect("in-range evaluation");
        let properties = evaluation.properties();
        let expected_cp = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 3.5;
        let expected_h = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 1_000.0 * 3.5;
        let expected_s = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (3.5 * fs_math::det::ln(1_000.0));
        assert_eq!(properties.cp().value().to_bits(), expected_cp.to_bits());
        assert_eq!(properties.h().value().to_bits(), expected_h.to_bits());
        assert_eq!(properties.s().value().to_bits(), expected_s.to_bits());
        assert_eq!(
            properties.u().value().to_bits(),
            (expected_h - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 1_000.0).to_bits(),
        );
        assert_eq!(
            properties.g().value().to_bits(),
            (expected_h - 1_000.0 * expected_s).to_bits(),
        );
        assert_eq!(evaluation.receipt().species(), model.species());
        assert_eq!(
            evaluation.receipt().evaluator_version(),
            NASA9_EVALUATOR_VERSION_V1,
        );
        assert_eq!(evaluation.receipt().fs_math_version(), fs_math::VERSION);
        assert_eq!(
            evaluation.receipt().gas_constant_bits(),
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits(),
        );
        assert_eq!(
            evaluation.receipt().molar_mass_bits(),
            0.028_013_4_f64.to_bits(),
        );
        assert_eq!(
            evaluation.receipt().temperature_bits(),
            1_000.0_f64.to_bits(),
        );
        assert_eq!(
            evaluation.receipt().reference_pressure_bits(),
            100_000.0_f64.to_bits(),
        );
        assert_eq!(evaluation.receipt().phase(), StandardStatePhaseV1::Gas);
        assert_eq!(
            evaluation.receipt().eos(),
            ReferenceEquationOfStateV1::IdealGas,
        );
        assert_eq!(
            evaluation.receipt().elemental_reference().as_str(),
            "nasa-glenn-2002-reference-elements",
        );
        assert_eq!(evaluation.receipt().region_index(), 0);
        assert_eq!(evaluation.receipt().region_min_bits(), 200.0_f64.to_bits(),);
        assert_eq!(
            evaluation.receipt().region_max_bits(),
            6_000.0_f64.to_bits(),
        );
    }

    #[test]
    fn g5_replay_and_shared_boundary_selection_are_deterministic() {
        let mut lower_coefficients = [0.0; 9];
        lower_coefficients[2] = 3.0;
        let mut upper_coefficients = [0.0; 9];
        upper_coefficients[2] = 4.0;
        let model = Nasa9StandardStateModelV1::new(
            SpeciesId::new("N2").expect("valid species"),
            MolarMass::new(0.028_013_4),
            convention(100_000.0),
            vec![
                Nasa9RegionV1::from_card(card(200.0, 1_000.0, 100_000.0, lower_coefficients))
                    .expect("lower region"),
                Nasa9RegionV1::from_card(card(1_000.0, 6_000.0, 100_000.0, upper_coefficients))
                    .expect("upper region"),
            ],
        )
        .expect("valid two-region model");
        let first = model
            .evaluate(Temperature::new(1_000.0))
            .expect("shared boundary selects upper");
        let replay = model
            .evaluate(Temperature::new(1_000.0))
            .expect("deterministic replay");
        assert_eq!(first, replay);
        assert_eq!(first.receipt().region_index(), 1);
        assert_eq!(
            model
                .evaluate(Temperature::new(200.0))
                .expect("first lower endpoint")
                .receipt()
                .region_index(),
            0,
        );
        assert_eq!(
            model
                .evaluate(Temperature::new(6_000.0))
                .expect("last upper endpoint")
                .receipt()
                .region_index(),
            1,
        );
        assert_eq!(
            first.properties().cp().value().to_bits(),
            (4.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K).to_bits(),
        );
        assert_ne!(
            first.receipt(),
            model
                .evaluate(Temperature::new(1_001.0))
                .expect("nearby temperature")
                .receipt(),
        );
    }

    #[test]
    fn g3_card_schema_and_convention_boundaries_fail_closed() {
        let coefficients = [0.0; 9];

        let mut wrong_law = card(200.0, 1_000.0, 100_000.0, coefficients);
        wrong_law.law = LawId("foreign-law".to_string());
        assert!(matches!(
            Nasa9RegionV1::from_card(wrong_law),
            Err(ThermochemErrorV1::CardLawMismatch { .. })
        ));

        let mut wrong_version = card(200.0, 1_000.0, 100_000.0, coefficients);
        wrong_version.law_version += 1;
        assert!(matches!(
            Nasa9RegionV1::from_card(wrong_version),
            Err(ThermochemErrorV1::CardLawVersionMismatch { .. })
        ));

        let mut stateful = card(200.0, 1_000.0, 100_000.0, coefficients);
        stateful.state_schema_version = 1;
        stateful.initial_state = InitialStatePolicy::RequiresDeclaredState;
        assert!(matches!(
            Nasa9RegionV1::from_card(stateful),
            Err(ThermochemErrorV1::CardStateConventionMismatch { .. })
        ));

        let mut missing = card(200.0, 1_000.0, 100_000.0, coefficients);
        missing.parameters.remove("a8");
        assert!(matches!(
            Nasa9RegionV1::from_card(missing),
            Err(ThermochemErrorV1::MissingCardParameter { parameter: "a8" })
        ));

        let mut unexpected = card(200.0, 1_000.0, 100_000.0, coefficients);
        unexpected.parameters.insert(
            "foreign".to_string(),
            LawParameter {
                value: 1.0,
                dims: Dims::NONE,
            },
        );
        assert!(matches!(
            Nasa9RegionV1::from_card(unexpected),
            Err(ThermochemErrorV1::UnexpectedCardParameter { .. })
        ));

        let mut wrong_axis = card(200.0, 1_000.0, 100_000.0, coefficients);
        wrong_axis.validity = ValidityDomain::unconstrained().with("P", 1.0, 2.0);
        assert!(matches!(
            Nasa9RegionV1::from_card(wrong_axis),
            Err(ThermochemErrorV1::InvalidTemperatureValidity { .. })
        ));

        assert!(matches!(
            ElementalReferenceIdV1::new("-not-canonical"),
            Err(ThermochemErrorV1::InvalidElementalReference { .. })
        ));
    }

    #[test]
    fn g3_cards_refuse_wrong_dimensions_and_nonfinite_coefficients() {
        let coefficients = [0.0; 9];
        let mut wrong_dims = card(200.0, 1_000.0, 100_000.0, coefficients);
        wrong_dims.parameters.get_mut("a0").expect("a0").dims = Dims::NONE;
        assert!(matches!(
            Nasa9RegionV1::from_card(wrong_dims),
            Err(ThermochemErrorV1::CardParameterDimsMismatch {
                parameter: "a0",
                ..
            })
        ));

        let invalid_coefficient = {
            let mut values = coefficients;
            values[4] = f64::NAN;
            card(200.0, 1_000.0, 100_000.0, values)
        };
        assert!(matches!(
            Nasa9RegionV1::from_card(invalid_coefficient),
            Err(ThermochemErrorV1::MaterialCard(
                MatDbError::NonFiniteParameter { .. }
            ))
        ));
    }

    #[test]
    fn g3_model_construction_refuses_invalid_boundaries() {
        let coefficients = [0.0; 9];
        assert!(matches!(
            StandardStateConventionV1::new(
                StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas,
                Pressure::new(0.0),
                ElementalReferenceIdV1::new("reference").expect("reference id"),
            ),
            Err(ThermochemErrorV1::InvalidReferencePressure { .. })
        ));
        assert!(matches!(
            Nasa9StandardStateModelV1::new(
                SpeciesId::new("N2").expect("species"),
                MolarMass::new(0.028),
                convention(100_000.0),
                Vec::new(),
            ),
            Err(ThermochemErrorV1::EmptyRegions)
        ));
        assert!(matches!(
            Nasa9StandardStateModelV1::new(
                SpeciesId::new("N2").expect("species"),
                MolarMass::new(0.0),
                convention(100_000.0),
                Vec::new(),
            ),
            Err(ThermochemErrorV1::InvalidMolarMass { .. })
        ));

        let bounded_region =
            Nasa9RegionV1::from_card(card(200.0, 1_000.0, 100_000.0, coefficients))
                .expect("bounded region");
        assert!(matches!(
            Nasa9StandardStateModelV1::new(
                SpeciesId::new("N2").expect("species"),
                MolarMass::new(0.028),
                convention(100_000.0),
                vec![bounded_region.clone(); MAX_NASA9_REGIONS_V1 + 1],
            ),
            Err(ThermochemErrorV1::TooManyRegions { .. })
        ));
        assert!(matches!(
            Nasa9StandardStateModelV1::new(
                SpeciesId::new("N2").expect("species"),
                MolarMass::new(0.028),
                convention(101_325.0),
                vec![bounded_region],
            ),
            Err(ThermochemErrorV1::RegionReferencePressureMismatch { .. })
        ));

        let lower = Nasa9RegionV1::from_card(card(200.0, 1_100.0, 100_000.0, coefficients))
            .expect("lower region");
        let upper = Nasa9RegionV1::from_card(card(1_000.0, 2_000.0, 100_000.0, coefficients))
            .expect("upper region");
        assert!(matches!(
            Nasa9StandardStateModelV1::new(
                SpeciesId::new("N2").expect("species"),
                MolarMass::new(0.028),
                convention(100_000.0),
                vec![lower, upper],
            ),
            Err(ThermochemErrorV1::OverlappingRegions { .. })
        ));
    }

    #[test]
    fn g3_evaluation_refuses_invalid_gapped_and_overflowing_inputs() {
        let coefficients = [0.0; 9];
        let model = model(coefficients);
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                model.evaluate(Temperature::new(invalid)),
                Err(ThermochemErrorV1::InvalidEvaluationTemperature { .. })
            ));
        }
        assert!(matches!(
            model.evaluate(Temperature::new(199.0)),
            Err(ThermochemErrorV1::TemperatureOutOfRange { .. })
        ));

        let gapped = Nasa9StandardStateModelV1::new(
            SpeciesId::new("N2").expect("species"),
            MolarMass::new(0.028),
            convention(100_000.0),
            vec![
                Nasa9RegionV1::from_card(card(200.0, 1_000.0, 100_000.0, coefficients))
                    .expect("lower gapped region"),
                Nasa9RegionV1::from_card(card(1_100.0, 2_000.0, 100_000.0, coefficients))
                    .expect("upper gapped region"),
            ],
        )
        .expect("gapped models are admitted explicitly");
        assert_eq!(
            gapped
                .evaluate(Temperature::new(1_000.0))
                .expect("isolated lower upper endpoint")
                .receipt()
                .region_index(),
            0,
        );
        assert!(matches!(
            gapped.evaluate(Temperature::new(1_050.0)),
            Err(ThermochemErrorV1::TemperatureOutOfRange { .. })
        ));
        assert_eq!(
            gapped
                .evaluate(Temperature::new(1_100.0))
                .expect("isolated upper lower endpoint")
                .receipt()
                .region_index(),
            1,
        );

        let mut overflowing_coefficients = [0.0; 9];
        overflowing_coefficients[6] = f64::MAX;
        let overflowing = Nasa9StandardStateModelV1::new(
            SpeciesId::new("N2").expect("species"),
            MolarMass::new(0.028),
            convention(100_000.0),
            vec![
                Nasa9RegionV1::from_card(card(1.0, 10.0, 100_000.0, overflowing_coefficients))
                    .expect("finite source coefficients"),
            ],
        )
        .expect("structurally valid model");
        assert!(matches!(
            overflowing.evaluate(Temperature::new(10.0)),
            Err(ThermochemErrorV1::NonFiniteEvaluation {
                property: Nasa9PropertyV1::HeatCapacity,
                ..
            })
        ));
    }

    #[test]
    fn g0_reexported_chemistry_certificate_is_the_single_exact_authority() {
        let h2 = SpeciesId::new("H2").expect("H2");
        let o2 = SpeciesId::new("O2").expect("O2");
        let h2o = SpeciesId::new("H2O").expect("H2O");
        let elemental = ElementalMatrix::new(
            vec![
                ElementId::new("H").expect("H"),
                ElementId::new("O").expect("O"),
            ],
            vec![h2.clone(), o2.clone(), h2o.clone()],
            vec![vec![2, 0, 2], vec![0, 2, 1]],
        )
        .expect("elemental matrix");
        let stoichiometric = StoichiometricMatrix::new(
            vec![h2.clone(), o2.clone(), h2o.clone()],
            vec![ReactionId::new("water-formation").expect("reaction")],
            vec![vec![-2], vec![-1], vec![2]],
        )
        .expect("stoichiometric matrix");
        let charge = ChargeVector::new(vec![h2, o2, h2o], vec![0, 0, 0]).expect("charge vector");
        let certificate = verify_conservation(&elemental, &stoichiometric, &charge)
            .expect("exact balanced mechanism");
        assert!(certificate.matches(&elemental, &stoichiometric, &charge));
    }
}
