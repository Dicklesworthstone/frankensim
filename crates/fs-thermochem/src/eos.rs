//! Mechanical ideal-gas equation-of-state closure.
//!
//! This first EOS rung maps a positive finite `(T, p, M)` state to mass
//! density, molar volume, specific gas constant, and `Z = 1`. It grants only
//! the algebra of the declared ideal-gas model. It does not establish phase
//! identity, phase stability, coefficient authenticity, caloric properties,
//! or physical validity for a real fluid.

use core::fmt;

use fs_qty::{Density, Dimensionless, MolarMass, Pressure, Qty, Temperature};

use crate::{MassSpecificGasConstantV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

/// Version of the fixed mechanical ideal-gas EOS operation tree and receipt.
pub const IDEAL_GAS_EOS_EVALUATOR_VERSION_V1: u32 = 1;

/// Coherent-SI molar-volume quantity, m^3/mol.
pub type IdealGasMolarVolumeQuantityV1 = Qty<3, 0, 0, 0, 0, -1>;

macro_rules! eos_quantity {
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

eos_quantity!(
    /// Ideal-gas mass density.
    IdealGasMassDensityV1,
    Density
);
eos_quantity!(
    /// Ideal-gas molar volume.
    IdealGasMolarVolumeV1,
    IdealGasMolarVolumeQuantityV1
);
eos_quantity!(
    /// Ideal-gas compressibility factor, identically one in this rung.
    IdealGasCompressibilityFactorV1,
    Dimensionless
);

/// Mechanical equation-of-state rung named by an evaluation receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MechanicalEquationOfStateV1 {
    /// `p V_m = R T`, with compressibility factor fixed to one.
    IdealGas,
}

/// EOS arithmetic value named by a non-representable result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IdealGasEosArithmeticValueV1 {
    /// Specific gas constant `R / M`.
    SpecificGasConstant,
    /// Intermediate molar energy scale `R T` used as the volume numerator.
    GasConstantTemperatureProduct,
    /// Molar volume `R T / p`.
    MolarVolume,
    /// Mass density `M / V_m`.
    MassDensity,
}

/// Mechanical ideal-gas construction or evaluation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum IdealGasEosErrorV1 {
    /// Model molar mass must be finite and strictly positive.
    InvalidMolarMass {
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// Evaluation temperature must be finite and strictly positive.
    InvalidTemperature {
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// Evaluation pressure must be finite and strictly positive.
    InvalidPressure {
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// Positive finite inputs produced a zero, negative, or non-finite value.
    UnrepresentableArithmeticValue {
        /// First output or intermediate that failed in operation-tree order.
        value: IdealGasEosArithmeticValueV1,
        /// Exact rejected arithmetic-value bits.
        bits: u64,
    },
}

impl fmt::Display for IdealGasEosErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMolarMass { bits } => write!(
                formatter,
                "ideal-gas EOS molar mass must be finite and positive (bits {bits:#018x})"
            ),
            Self::InvalidTemperature { bits } => write!(
                formatter,
                "ideal-gas EOS temperature must be finite and positive (bits {bits:#018x})"
            ),
            Self::InvalidPressure { bits } => write!(
                formatter,
                "ideal-gas EOS pressure must be finite and positive (bits {bits:#018x})"
            ),
            Self::UnrepresentableArithmeticValue { value, bits } => write!(
                formatter,
                "ideal-gas EOS {value:?} is not positive finite (bits {bits:#018x})"
            ),
        }
    }
}

impl std::error::Error for IdealGasEosErrorV1 {}

/// Immutable mechanical ideal-gas EOS law for one caller-declared molar mass.
///
/// Molar mass is retained as a typed scalar but is not bound to a species or
/// authenticated material source by this primitive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdealGasEosV1 {
    molar_mass: MolarMass,
}

impl IdealGasEosV1 {
    /// Admit a mechanical ideal-gas law with positive finite molar mass.
    ///
    /// # Errors
    /// Refuses zero, negative, NaN, or infinite molar mass.
    pub fn new(molar_mass: MolarMass) -> Result<Self, IdealGasEosErrorV1> {
        let value = molar_mass.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(IdealGasEosErrorV1::InvalidMolarMass {
                bits: value.to_bits(),
            });
        }
        Ok(Self { molar_mass })
    }

    /// Caller-declared positive finite molar mass.
    #[must_use]
    pub const fn molar_mass(self) -> MolarMass {
        self.molar_mass
    }

    /// Evaluate the mechanical ideal-gas law at absolute temperature and
    /// pressure.
    ///
    /// Version 1 uses the fixed order `R_specific = R/M`, `RT = R*T`,
    /// `V_m = RT/p`, `rho = M/V_m`, and `Z = 1`. Every exposed derived scalar
    /// must remain strictly positive and finite.
    ///
    /// # Errors
    /// Refuses non-positive or non-finite state inputs and the first zero,
    /// negative, or non-finite derived value in operation-tree order.
    pub fn evaluate_pt(
        &self,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<IdealGasEosEvaluationV1, IdealGasEosErrorV1> {
        let temperature_value = temperature.value();
        if !temperature_value.is_finite() || temperature_value <= 0.0 {
            return Err(IdealGasEosErrorV1::InvalidTemperature {
                bits: temperature_value.to_bits(),
            });
        }
        let pressure_value = pressure.value();
        if !pressure_value.is_finite() || pressure_value <= 0.0 {
            return Err(IdealGasEosErrorV1::InvalidPressure {
                bits: pressure_value.to_bits(),
            });
        }

        let molar_mass_value = self.molar_mass.value();
        let specific_gas_constant = require_positive_finite(
            IdealGasEosArithmeticValueV1::SpecificGasConstant,
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass_value,
        )?;
        let gas_constant_temperature = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature_value;
        if !gas_constant_temperature.is_finite() || gas_constant_temperature <= 0.0 {
            return Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
                value: IdealGasEosArithmeticValueV1::GasConstantTemperatureProduct,
                bits: gas_constant_temperature.to_bits(),
            });
        }
        let molar_volume = require_positive_finite(
            IdealGasEosArithmeticValueV1::MolarVolume,
            gas_constant_temperature / pressure_value,
        )?;
        let mass_density = require_positive_finite(
            IdealGasEosArithmeticValueV1::MassDensity,
            molar_mass_value / molar_volume,
        )?;
        let compressibility_factor = 1.0f64;

        let receipt = IdealGasEosReceiptV1 {
            evaluator_version: IDEAL_GAS_EOS_EVALUATOR_VERSION_V1,
            fs_qty_version: fs_qty::VERSION,
            eos: MechanicalEquationOfStateV1::IdealGas,
            gas_constant_bits: UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits(),
            molar_mass_bits: molar_mass_value.to_bits(),
            temperature_bits: temperature_value.to_bits(),
            pressure_bits: pressure_value.to_bits(),
            mass_density_bits: mass_density.to_bits(),
            molar_volume_bits: molar_volume.to_bits(),
            specific_gas_constant_bits: specific_gas_constant.to_bits(),
            compressibility_factor_bits: compressibility_factor.to_bits(),
        };

        Ok(IdealGasEosEvaluationV1 {
            molar_mass: self.molar_mass,
            temperature,
            pressure,
            mass_density: IdealGasMassDensityV1::new(mass_density),
            molar_volume: IdealGasMolarVolumeV1::new(molar_volume),
            specific_gas_constant: MassSpecificGasConstantV1::new(specific_gas_constant),
            compressibility_factor: IdealGasCompressibilityFactorV1::new(compressibility_factor),
            receipt,
        })
    }
}

fn require_positive_finite(
    arithmetic_value: IdealGasEosArithmeticValueV1,
    value: f64,
) -> Result<f64, IdealGasEosErrorV1> {
    if !value.is_finite() || value <= 0.0 {
        return Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
            value: arithmetic_value,
            bits: value.to_bits(),
        });
    }
    Ok(value)
}

/// Immutable exact-field receipt for one mechanical ideal-gas evaluation.
///
/// The receipt binds the operation-tree version and exact scalar inputs and
/// outputs. It records provenance but does not authenticate molar mass, species,
/// phase, or material identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdealGasEosReceiptV1 {
    evaluator_version: u32,
    fs_qty_version: &'static str,
    eos: MechanicalEquationOfStateV1,
    gas_constant_bits: u64,
    molar_mass_bits: u64,
    temperature_bits: u64,
    pressure_bits: u64,
    mass_density_bits: u64,
    molar_volume_bits: u64,
    specific_gas_constant_bits: u64,
    compressibility_factor_bits: u64,
}

impl IdealGasEosReceiptV1 {
    /// EOS evaluator/operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Quantity-semantics crate version.
    #[must_use]
    pub const fn fs_qty_version(&self) -> &'static str {
        self.fs_qty_version
    }

    /// Explicit equation-of-state rung.
    #[must_use]
    pub const fn eos(&self) -> MechanicalEquationOfStateV1 {
        self.eos
    }

    /// Exact universal-gas-constant bits.
    #[must_use]
    pub const fn gas_constant_bits(&self) -> u64 {
        self.gas_constant_bits
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

    /// Exact evaluation-pressure bits.
    #[must_use]
    pub const fn pressure_bits(&self) -> u64 {
        self.pressure_bits
    }

    /// Exact mass-density output bits.
    #[must_use]
    pub const fn mass_density_bits(&self) -> u64 {
        self.mass_density_bits
    }

    /// Exact molar-volume output bits.
    #[must_use]
    pub const fn molar_volume_bits(&self) -> u64 {
        self.molar_volume_bits
    }

    /// Exact specific-gas-constant output bits.
    #[must_use]
    pub const fn specific_gas_constant_bits(&self) -> u64 {
        self.specific_gas_constant_bits
    }

    /// Exact compressibility-factor output bits.
    #[must_use]
    pub const fn compressibility_factor_bits(&self) -> u64 {
        self.compressibility_factor_bits
    }
}

/// One successful mechanical ideal-gas EOS evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasEosEvaluationV1 {
    molar_mass: MolarMass,
    temperature: Temperature,
    pressure: Pressure,
    mass_density: IdealGasMassDensityV1,
    molar_volume: IdealGasMolarVolumeV1,
    specific_gas_constant: MassSpecificGasConstantV1,
    compressibility_factor: IdealGasCompressibilityFactorV1,
    receipt: IdealGasEosReceiptV1,
}

impl IdealGasEosEvaluationV1 {
    /// Caller-declared molar mass used by the EOS.
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }

    /// Evaluated absolute temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Evaluated absolute pressure.
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }

    /// Positive finite ideal-gas mass density.
    #[must_use]
    pub const fn mass_density(&self) -> IdealGasMassDensityV1 {
        self.mass_density
    }

    /// Positive finite ideal-gas molar volume.
    #[must_use]
    pub const fn molar_volume(&self) -> IdealGasMolarVolumeV1 {
        self.molar_volume
    }

    /// Positive finite specific gas constant.
    #[must_use]
    pub const fn specific_gas_constant(&self) -> MassSpecificGasConstantV1 {
        self.specific_gas_constant
    }

    /// Compressibility factor, exactly one for this ideal-gas rung.
    #[must_use]
    pub const fn compressibility_factor(&self) -> IdealGasCompressibilityFactorV1 {
        self.compressibility_factor
    }

    /// Exact operation-tree and scalar receipt.
    #[must_use]
    pub const fn receipt(&self) -> &IdealGasEosReceiptV1 {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const AIR_MOLAR_MASS: f64 = 0.028_965_46;
    const STANDARD_PRESSURE: f64 = 101_325.0;

    fn assert_relative_eq(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        let tolerance = 64.0 * f64::EPSILON * scale;
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.17e}"
        );
    }

    fn air() -> IdealGasEosV1 {
        IdealGasEosV1::new(MolarMass::new(AIR_MOLAR_MASS)).expect("air molar mass")
    }

    // G0: all exposed fields satisfy the defining ideal-gas identities.
    #[test]
    fn g0_ideal_gas_pt_identities_hold_in_typed_outputs() {
        let temperature = Temperature::new(300.0);
        let pressure = Pressure::new(STANDARD_PRESSURE);
        let evaluation = air()
            .evaluate_pt(temperature, pressure)
            .expect("positive ideal-gas state");

        let molar_volume = evaluation.molar_volume().value();
        let density = evaluation.mass_density().value();
        let specific_gas_constant = evaluation.specific_gas_constant().value();
        assert_relative_eq(
            pressure.value() * molar_volume,
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value(),
        );
        assert_relative_eq(
            density * specific_gas_constant * temperature.value(),
            pressure.value(),
        );
        assert_relative_eq(density * molar_volume, AIR_MOLAR_MASS);
        assert_eq!(
            evaluation.compressibility_factor().value().to_bits(),
            1.0f64.to_bits()
        );

        assert_eq!(
            evaluation.mass_density().quantity().value().to_bits(),
            density.to_bits()
        );
        assert_eq!(
            evaluation.molar_volume().quantity().value().to_bits(),
            molar_volume.to_bits()
        );
        assert_eq!(
            evaluation
                .specific_gas_constant()
                .quantity()
                .value()
                .to_bits(),
            specific_gas_constant.to_bits()
        );
        assert_eq!(
            evaluation
                .compressibility_factor()
                .quantity()
                .value()
                .to_bits(),
            1.0f64.to_bits()
        );
    }

    // G3: equivalent pressure and temperature rescalings preserve the law.
    #[test]
    fn g3_pressure_and_temperature_scaling_is_metamorphic() {
        let base = air()
            .evaluate_pt(Temperature::new(300.0), Pressure::new(STANDARD_PRESSURE))
            .expect("base state");
        let doubled_pressure = air()
            .evaluate_pt(
                Temperature::new(300.0),
                Pressure::new(2.0 * STANDARD_PRESSURE),
            )
            .expect("pressure-scaled state");
        let doubled_temperature = air()
            .evaluate_pt(Temperature::new(600.0), Pressure::new(STANDARD_PRESSURE))
            .expect("temperature-scaled state");

        assert_relative_eq(
            doubled_pressure.mass_density().value(),
            2.0 * base.mass_density().value(),
        );
        assert_relative_eq(
            doubled_pressure.molar_volume().value(),
            0.5 * base.molar_volume().value(),
        );
        assert_relative_eq(
            doubled_temperature.mass_density().value(),
            0.5 * base.mass_density().value(),
        );
        assert_relative_eq(
            doubled_temperature.molar_volume().value(),
            2.0 * base.molar_volume().value(),
        );
        assert_eq!(
            doubled_temperature.specific_gas_constant(),
            base.specific_gas_constant()
        );
    }

    // G3: malformed states and finite arithmetic extremes fail closed.
    #[test]
    fn g3_invalid_and_unrepresentable_states_refuse_without_partial_output() {
        for value in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                IdealGasEosV1::new(MolarMass::new(value)),
                Err(IdealGasEosErrorV1::InvalidMolarMass { bits })
                    if bits == value.to_bits()
            ));
        }

        let model = air();
        for value in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                model.evaluate_pt(Temperature::new(value), Pressure::new(1.0)),
                Err(IdealGasEosErrorV1::InvalidTemperature { bits })
                    if bits == value.to_bits()
            ));
            assert!(matches!(
                model.evaluate_pt(Temperature::new(1.0), Pressure::new(value)),
                Err(IdealGasEosErrorV1::InvalidPressure { bits })
                    if bits == value.to_bits()
            ));
        }

        assert!(matches!(
            model.evaluate_pt(Temperature::new(f64::MAX), Pressure::new(1.0)),
            Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
                value: IdealGasEosArithmeticValueV1::GasConstantTemperatureProduct,
                bits,
            }) if bits == f64::INFINITY.to_bits()
        ));
        assert!(matches!(
            model.evaluate_pt(
                Temperature::new(f64::MIN_POSITIVE),
                Pressure::new(f64::MAX),
            ),
            Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
                value: IdealGasEosArithmeticValueV1::MolarVolume,
                bits,
            }) if bits == 0.0f64.to_bits()
        ));
        let tiny_mass = IdealGasEosV1::new(MolarMass::new(f64::MIN_POSITIVE))
            .expect("positive finite tiny molar mass is structurally valid");
        assert!(matches!(
            tiny_mass.evaluate_pt(Temperature::new(1.0), Pressure::new(1.0)),
            Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
                value: IdealGasEosArithmeticValueV1::SpecificGasConstant,
                bits,
            }) if bits == f64::INFINITY.to_bits()
        ));

        let huge_mass = IdealGasEosV1::new(MolarMass::new(f64::MAX))
            .expect("positive finite huge molar mass is structurally valid");
        assert!(matches!(
            huge_mass.evaluate_pt(
                Temperature::new(1.0),
                Pressure::new(f64::MAX),
            ),
            Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
                value: IdealGasEosArithmeticValueV1::MassDensity,
                bits,
            }) if bits == f64::INFINITY.to_bits()
        ));

        let small_mass = IdealGasEosV1::new(MolarMass::new(1.0e-307))
            .expect("positive finite small molar mass is structurally valid");
        assert!(matches!(
            small_mass.evaluate_pt(
                Temperature::new(1.0),
                Pressure::new(1.0e-307),
            ),
            Err(IdealGasEosErrorV1::UnrepresentableArithmeticValue {
                value: IdealGasEosArithmeticValueV1::MassDensity,
                bits,
            }) if bits == 0.0f64.to_bits()
        ));
    }

    // G5: exact operation inputs and outputs replay bit-identically.
    #[test]
    fn g5_receipt_replays_and_binds_every_scalar_field() {
        let model = air();
        let temperature = Temperature::new(350.0);
        let pressure = Pressure::new(250_000.0);
        let first = model
            .evaluate_pt(temperature, pressure)
            .expect("first replay state");
        let second = model
            .evaluate_pt(temperature, pressure)
            .expect("second replay state");
        assert_eq!(first, second);

        let receipt = first.receipt();
        assert_eq!(
            receipt.evaluator_version(),
            IDEAL_GAS_EOS_EVALUATOR_VERSION_V1
        );
        assert_eq!(receipt.fs_qty_version(), fs_qty::VERSION);
        assert_eq!(receipt.eos(), MechanicalEquationOfStateV1::IdealGas);
        assert_eq!(
            receipt.gas_constant_bits(),
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits()
        );
        assert_eq!(receipt.molar_mass_bits(), AIR_MOLAR_MASS.to_bits());
        assert_eq!(receipt.temperature_bits(), temperature.value().to_bits());
        assert_eq!(receipt.pressure_bits(), pressure.value().to_bits());
        assert_eq!(
            receipt.mass_density_bits(),
            first.mass_density().value().to_bits()
        );
        assert_eq!(
            receipt.molar_volume_bits(),
            first.molar_volume().value().to_bits()
        );
        assert_eq!(
            receipt.specific_gas_constant_bits(),
            first.specific_gas_constant().value().to_bits()
        );
        assert_eq!(
            receipt.compressibility_factor_bits(),
            first.compressibility_factor().value().to_bits()
        );

        let mut mutated = receipt.clone();
        mutated.temperature_bits = 300.0f64.to_bits();
        assert_ne!(mutated, *receipt);
        let pressure_changed = model
            .evaluate_pt(temperature, Pressure::new(300_000.0))
            .expect("changed pressure");
        assert_ne!(
            pressure_changed.receipt().pressure_bits(),
            receipt.pressure_bits()
        );
        assert_ne!(
            pressure_changed.receipt().mass_density_bits(),
            receipt.mass_density_bits()
        );
    }
}
