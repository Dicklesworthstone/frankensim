//! Semantic quantity kinds layered over the six-base dimensional algebra.
//!
//! Dimensions answer whether two values have compatible units. This module
//! answers the separate question of whether they mean the same thing. Its
//! carriers keep their fields private, validate dimensions and scalar domains
//! at construction, and require named conversions wherever a dimensionally
//! legal operation could hide an affine offset, `2*pi`, pole-pair, `sqrt(2)`,
//! or composition-basis error.
//!
//! The scope is deliberately small. It does not implement general signal
//! processing, complex arithmetic, density-based volume-fraction conversion,
//! or logarithms. In particular, acoustic levels store an explicit validated
//! reference but leave deterministic linear/log conversion to the owning
//! acoustics/math layer.

use crate::{Dims, MolarMass, QtyAny};
use core::fmt;
use core::num::NonZeroU32;

const TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const ANGULAR_VELOCITY_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);
const ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MASS_DIMS: Dims = Dims([0, 1, 0, 0, 0, 0]);
const AMOUNT_DIMS: Dims = Dims([0, 0, 0, 0, 0, 1]);
const MOLAR_MASS_DIMS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const MASS_CONCENTRATION_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const AMOUNT_CONCENTRATION_DIMS: Dims = Dims([-3, 0, 0, 0, 0, 1]);
const ENTROPY_DIMS: Dims = Dims([2, 1, -2, -1, 0, 0]);
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);
const COMPOSITION_SUM_TOLERANCE: f64 = 1.0e-12;
const RADIANS_PER_SECOND_PER_RPM: f64 = core::f64::consts::TAU / 60.0;
const RPM_PER_RADIAN_PER_SECOND: f64 = 60.0 / core::f64::consts::TAU;
const STATIC_FORM: u8 = 1 << 0;
const INSTANTANEOUS_FORM: u8 = 1 << 1;
const PEAK_FORM: u8 = 1 << 2;
const RMS_FORM: u8 = 1 << 3;
const PHASOR_FORM: u8 = 1 << 4;
const WAVEFORM_FORMS: u8 = STATIC_FORM | INSTANTANEOUS_FORM | PEAK_FORM | RMS_FORM | PHASOR_FORM;

/// Domain attached to an angle or angular velocity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AngleDomain {
    /// Mechanical shaft or geometry angle.
    Mechanical,
    /// Electrical phase angle.
    Electrical,
}

/// Convention used to store a strain component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StrainBasis {
    /// Tensor strain convention.
    Tensor,
    /// Engineering strain convention.
    Engineering,
}

/// Component role needed to convert strain conventions without doubling a
/// normal component accidentally.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StrainComponent {
    /// Normal strain; tensor and engineering values are numerically equal.
    Normal,
    /// Shear strain; engineering shear is twice tensor shear.
    Shear,
}

/// Basis used by a dimensionless composition vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CompositionBasis {
    /// Component mass fractions.
    MassFraction,
    /// Component amount-of-substance (mole) fractions.
    MoleFraction,
    /// Component volume fractions.
    VolumeFraction,
}

/// Physical meaning carried in addition to a dimension vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QuantityKind {
    /// Affine absolute thermodynamic temperature.
    AbsoluteTemperature,
    /// Linear temperature interval.
    TemperatureDifference,
    /// Angle in the named mechanical/electrical domain.
    Angle(AngleDomain),
    /// Angular velocity in the named mechanical/electrical domain.
    AngularVelocity(AngleDomain),
    /// Torque, distinct from energy despite sharing dimensions.
    Torque,
    /// Energy, distinct from torque despite sharing dimensions.
    Energy,
    /// Thermodynamic/fluid pressure.
    Pressure,
    /// Solid stress.
    Stress,
    /// A scalar strain component with an explicit convention.
    Strain {
        /// Tensor or engineering convention.
        basis: StrainBasis,
        /// Normal or shear component role.
        component: StrainComponent,
    },
    /// A component of a composition vector in the named basis.
    Composition(CompositionBasis),
    /// Mass.
    Mass,
    /// Amount of substance.
    Amount,
    /// Mass per amount of substance.
    MolarMass,
    /// Mass per volume.
    MassConcentration,
    /// Amount of substance per volume.
    AmountConcentration,
    /// Entropy.
    Entropy,
    /// Heat capacity, distinct from entropy despite sharing dimensions.
    HeatCapacity,
    /// Physical acoustic pressure.
    AcousticPressure,
    /// Physical acoustic power.
    AcousticPower,
}

impl QuantityKind {
    /// Required six-base dimension vector for this semantic kind.
    #[must_use]
    pub const fn expected_dims(self) -> Dims {
        match self {
            Self::AbsoluteTemperature | Self::TemperatureDifference => TEMPERATURE_DIMS,
            Self::Angle(_) | Self::Strain { .. } | Self::Composition(_) => Dims::NONE,
            Self::AngularVelocity(_) => ANGULAR_VELOCITY_DIMS,
            Self::Torque | Self::Energy => ENERGY_DIMS,
            Self::Pressure | Self::Stress | Self::AcousticPressure => PRESSURE_DIMS,
            Self::Mass => MASS_DIMS,
            Self::Amount => AMOUNT_DIMS,
            Self::MolarMass => MOLAR_MASS_DIMS,
            Self::MassConcentration => MASS_CONCENTRATION_DIMS,
            Self::AmountConcentration => AMOUNT_CONCENTRATION_DIMS,
            Self::Entropy | Self::HeatCapacity => ENTROPY_DIMS,
            Self::AcousticPower => POWER_DIMS,
        }
    }

    /// Whether the sealed kind/form matrix admits this real-scalar form.
    ///
    /// This exhaustive match is the authority boundary: adding a new kind
    /// requires choosing its form policy instead of inheriting a permissive
    /// default.
    #[must_use]
    pub const fn admits_scalar_form(self, form: ValueForm) -> bool {
        self.form_mask() & scalar_form_bit(form) != 0
    }

    /// Whether the sealed kind/form matrix admits a paired complex phasor.
    #[must_use]
    pub const fn admits_phasor(self) -> bool {
        self.form_mask() & PHASOR_FORM != 0
    }

    const fn form_mask(self) -> u8 {
        match self {
            Self::TemperatureDifference
            | Self::Angle(_)
            | Self::AngularVelocity(_)
            | Self::Torque
            | Self::Pressure
            | Self::Stress
            | Self::Strain { .. }
            | Self::AcousticPressure => WAVEFORM_FORMS,
            Self::AbsoluteTemperature
            | Self::Energy
            | Self::Composition(_)
            | Self::Mass
            | Self::Amount
            | Self::MolarMass
            | Self::MassConcentration
            | Self::AmountConcentration
            | Self::Entropy
            | Self::HeatCapacity
            | Self::AcousticPower => STATIC_FORM,
        }
    }
}

/// Scalar value form. A complex phasor is intentionally represented by the
/// paired [`PhasorQty`] carrier instead of two independently taggable scalars.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValueForm {
    /// State or aggregate scalar without waveform semantics.
    Static,
    /// Instantaneous sample.
    Instantaneous,
    /// Nonnegative peak amplitude.
    Peak,
    /// Nonnegative root-mean-square amplitude.
    Rms,
}

const fn scalar_form_bit(form: ValueForm) -> u8 {
    match form {
        ValueForm::Static => STATIC_FORM,
        ValueForm::Instantaneous => INSTANTANEOUS_FORM,
        ValueForm::Peak => PEAK_FORM,
        ValueForm::Rms => RMS_FORM,
    }
}

/// Requested semantic type of a real scalar.
///
/// This descriptor may be assembled for diagnostics. Value carriers enforce
/// the sealed [`QuantityKind`] form-admissibility matrix before retaining it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SemanticType {
    kind: QuantityKind,
    form: ValueForm,
}

impl SemanticType {
    /// Construct a semantic type.
    #[must_use]
    pub const fn new(kind: QuantityKind, form: ValueForm) -> Self {
        Self { kind, form }
    }

    /// Physical quantity kind.
    #[must_use]
    pub const fn kind(self) -> QuantityKind {
        self.kind
    }

    /// Scalar value form.
    #[must_use]
    pub const fn form(self) -> ValueForm {
        self.form
    }

    /// Dimension vector required by the quantity kind.
    #[must_use]
    pub const fn expected_dims(self) -> Dims {
        self.kind.expected_dims()
    }
}

/// Scalar-domain rule named by [`SemanticError::InvalidValue`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValueRequirement {
    /// A finite real value.
    Finite,
    /// A finite value greater than or equal to zero.
    FiniteNonnegative,
    /// A finite value strictly greater than zero.
    FinitePositive,
    /// A finite component fraction in the closed interval `[0, 1]`.
    UnitFraction,
}

/// Value-form policy named by [`SemanticError::UnsupportedForm`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FormRequirement {
    /// This kind admits only [`ValueForm::Static`].
    StaticOnly,
    /// A point value: [`ValueForm::Static`] or [`ValueForm::Instantaneous`].
    PointValue,
}

/// Structured semantic validation and conversion error.
#[derive(Debug, Clone, PartialEq)]
pub enum SemanticError {
    /// The runtime dimension vector does not match the semantic kind.
    DimensionMismatch {
        /// Named operation or boundary.
        operation: &'static str,
        /// Semantic type requested by the operation.
        semantic_type: SemanticType,
        /// Dimension vector supplied by the caller.
        actual: Dims,
        /// Required dimension vector.
        expected: Dims,
    },
    /// A value with one semantic type was supplied where another was required.
    KindMismatch {
        /// Named operation or boundary.
        operation: &'static str,
        /// Actual source type.
        source: SemanticType,
        /// Required or requested target type.
        target: SemanticType,
    },
    /// A scalar violated a finite/range requirement.
    InvalidValue {
        /// Named operation or boundary.
        operation: &'static str,
        /// Semantic type of the value.
        semantic_type: SemanticType,
        /// Rejected value.
        value: f64,
        /// Required scalar domain.
        requirement: ValueRequirement,
    },
    /// An aggregate amplitude form was supplied to a point-value operation.
    UnsupportedForm {
        /// Named operation or boundary.
        operation: &'static str,
        /// Exact semantic source type, including its rejected form.
        source: SemanticType,
        /// Form policy required by the operation.
        requirement: FormRequirement,
    },
    /// A pole-pair phase map was constructed with zero pole pairs.
    ZeroPolePairs,
    /// A composition vector was empty.
    EmptyComposition {
        /// Declared composition basis.
        basis: CompositionBasis,
    },
    /// A composition component was non-finite or outside `[0, 1]`.
    InvalidCompositionEntry {
        /// Declared composition basis.
        basis: CompositionBasis,
        /// Component index.
        index: usize,
        /// Rejected component value.
        value: f64,
    },
    /// Composition components did not sum to one within the fixed tolerance.
    InvalidCompositionSum {
        /// Declared composition basis.
        basis: CompositionBasis,
        /// Deterministic input-order sum.
        sum: f64,
        /// Absolute acceptance tolerance.
        tolerance: f64,
    },
    /// A whole-vector basis conversion received the wrong molar-mass count.
    CompositionLengthMismatch {
        /// Number of fractions.
        fractions: usize,
        /// Number of molar masses.
        molar_masses: usize,
    },
    /// A whole-vector basis conversion received a non-positive molar mass.
    InvalidMolarMass {
        /// Component index.
        index: usize,
        /// Rejected coherent-SI molar mass.
        value: f64,
    },
    /// Positive mass/amount basis conversion rounded to an unrepresentable
    /// zero result.
    UnrepresentablePositiveConversion {
        /// Named conversion boundary.
        operation: &'static str,
        /// Exact semantic source type.
        source: SemanticType,
        /// Exact semantic target type.
        target: SemanticType,
        /// Positive coherent-SI source value that would be lost.
        source_value: f64,
        /// Positive coherent-SI molar mass used by the conversion.
        molar_mass: f64,
    },
}

impl fmt::Display for SemanticError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionMismatch {
                operation,
                semantic_type,
                actual,
                expected,
            } => write!(
                f,
                "dimension mismatch in {operation} for {semantic_type:?}: actual [{actual:?}], expected [{expected:?}]"
            ),
            Self::KindMismatch {
                operation,
                source,
                target,
            } => write!(
                f,
                "semantic kind mismatch in {operation}: source {source:?}, target {target:?}"
            ),
            Self::InvalidValue {
                operation,
                semantic_type,
                value,
                requirement,
            } => write!(
                f,
                "invalid value in {operation} for {semantic_type:?}: {value} does not satisfy {requirement:?}"
            ),
            Self::UnsupportedForm {
                operation,
                source,
                requirement,
            } => write!(
                f,
                "unsupported value form in {operation}: source {source:?} does not satisfy {requirement:?}"
            ),
            Self::ZeroPolePairs => {
                f.write_str("pole-pair phase map requires at least one pole pair")
            }
            Self::EmptyComposition { basis } => {
                write!(
                    f,
                    "{basis:?} composition must contain at least one component"
                )
            }
            Self::InvalidCompositionEntry {
                basis,
                index,
                value,
            } => write!(
                f,
                "invalid {basis:?} component at index {index}: {value} is not a finite unit fraction"
            ),
            Self::InvalidCompositionSum {
                basis,
                sum,
                tolerance,
            } => write!(
                f,
                "invalid {basis:?} composition sum {sum}; expected 1 within {tolerance}"
            ),
            Self::CompositionLengthMismatch {
                fractions,
                molar_masses,
            } => write!(
                f,
                "composition basis conversion needs one molar mass per fraction: {fractions} fractions, {molar_masses} molar masses"
            ),
            Self::InvalidMolarMass { index, value } => write!(
                f,
                "invalid molar mass at component {index}: {value} must be finite and positive"
            ),
            Self::UnrepresentablePositiveConversion {
                operation,
                source,
                target,
                source_value,
                molar_mass,
            } => write!(
                f,
                "unrepresentable positive basis conversion in {operation}: {source:?} value {source_value} with molar mass {molar_mass} would round {target:?} to zero"
            ),
        }
    }
}

impl core::error::Error for SemanticError {}

/// Privately validated runtime quantity plus its semantic type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SemanticQty {
    quantity: QtyAny,
    semantic_type: SemanticType,
}

impl SemanticQty {
    /// Attach a semantic type after validating dimensions, finiteness, and the
    /// small set of universal scalar-domain rules.
    ///
    /// # Errors
    /// Returns a structured [`SemanticError`] when the value is inadmissible.
    pub fn new(quantity: QtyAny, semantic_type: SemanticType) -> Result<Self, SemanticError> {
        Self::new_for_operation(quantity, semantic_type, "construct semantic quantity")
    }

    fn new_for_operation(
        quantity: QtyAny,
        semantic_type: SemanticType,
        operation: &'static str,
    ) -> Result<Self, SemanticError> {
        let expected = semantic_type.expected_dims();
        if quantity.dims != expected {
            return Err(SemanticError::DimensionMismatch {
                operation,
                semantic_type,
                actual: quantity.dims,
                expected,
            });
        }
        if !semantic_type
            .kind()
            .admits_scalar_form(semantic_type.form())
        {
            return Err(SemanticError::UnsupportedForm {
                operation,
                source: semantic_type,
                requirement: FormRequirement::StaticOnly,
            });
        }
        validate_scalar(operation, semantic_type, quantity.value)?;
        Ok(Self {
            quantity,
            semantic_type,
        })
    }

    /// Runtime dimensional value in coherent SI units.
    #[must_use]
    pub const fn quantity(self) -> QtyAny {
        self.quantity
    }

    /// Raw coherent-SI scalar.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.quantity.value
    }

    /// Semantic type carried by the value.
    #[must_use]
    pub const fn semantic_type(self) -> SemanticType {
        self.semantic_type
    }

    /// Verify an exact semantic type without retagging the value.
    ///
    /// # Errors
    /// Returns [`SemanticError::KindMismatch`] with typed source and target.
    pub fn require_type(
        self,
        target: SemanticType,
        operation: &'static str,
    ) -> Result<Self, SemanticError> {
        require_type(self, target, operation)
    }
}

fn validate_scalar(
    operation: &'static str,
    semantic_type: SemanticType,
    value: f64,
) -> Result<(), SemanticError> {
    if !value.is_finite() {
        return Err(invalid_value(
            operation,
            semantic_type,
            value,
            ValueRequirement::Finite,
        ));
    }

    if matches!(semantic_type.form(), ValueForm::Peak | ValueForm::Rms) && value < 0.0 {
        return Err(invalid_value(
            operation,
            semantic_type,
            value,
            ValueRequirement::FiniteNonnegative,
        ));
    }

    match semantic_type.kind() {
        QuantityKind::AbsoluteTemperature
        | QuantityKind::Mass
        | QuantityKind::Amount
        | QuantityKind::MassConcentration
        | QuantityKind::AmountConcentration
            if value < 0.0 =>
        {
            Err(invalid_value(
                operation,
                semantic_type,
                value,
                ValueRequirement::FiniteNonnegative,
            ))
        }
        QuantityKind::MolarMass if value <= 0.0 => Err(invalid_value(
            operation,
            semantic_type,
            value,
            ValueRequirement::FinitePositive,
        )),
        QuantityKind::Composition(_) if !(0.0..=1.0).contains(&value) => Err(invalid_value(
            operation,
            semantic_type,
            value,
            ValueRequirement::UnitFraction,
        )),
        _ => Ok(()),
    }
}

const fn invalid_value(
    operation: &'static str,
    semantic_type: SemanticType,
    value: f64,
    requirement: ValueRequirement,
) -> SemanticError {
    SemanticError::InvalidValue {
        operation,
        semantic_type,
        value,
        requirement,
    }
}

fn require_type(
    source: SemanticQty,
    target: SemanticType,
    operation: &'static str,
) -> Result<SemanticQty, SemanticError> {
    if source.semantic_type == target {
        Ok(source)
    } else {
        Err(SemanticError::KindMismatch {
            operation,
            source: source.semantic_type,
            target,
        })
    }
}

fn require_point_form(
    source: SemanticQty,
    operation: &'static str,
) -> Result<SemanticQty, SemanticError> {
    if matches!(
        source.semantic_type.form(),
        ValueForm::Static | ValueForm::Instantaneous
    ) {
        Ok(source)
    } else {
        Err(SemanticError::UnsupportedForm {
            operation,
            source: source.semantic_type,
            requirement: FormRequirement::PointValue,
        })
    }
}

fn require_finite_result(
    value: f64,
    semantic_type: SemanticType,
    operation: &'static str,
) -> Result<f64, SemanticError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid_value(
            operation,
            semantic_type,
            value,
            ValueRequirement::Finite,
        ))
    }
}

fn semantic_value(
    value: f64,
    semantic_type: SemanticType,
    operation: &'static str,
) -> Result<SemanticQty, SemanticError> {
    SemanticQty::new_for_operation(
        QtyAny::new(value, semantic_type.expected_dims()),
        semantic_type,
        operation,
    )
}

const fn kind_type(kind: QuantityKind, form: ValueForm) -> SemanticType {
    SemanticType::new(kind, form)
}

/// Construct a static absolute temperature from kelvin.
pub fn absolute_temperature_kelvin(value: f64) -> Result<SemanticQty, SemanticError> {
    semantic_value(
        value,
        kind_type(QuantityKind::AbsoluteTemperature, ValueForm::Static),
        "construct absolute temperature from kelvin",
    )
}

/// Construct a static absolute temperature from degrees Celsius.
pub fn absolute_temperature_celsius(value: f64) -> Result<SemanticQty, SemanticError> {
    semantic_value(
        value + 273.15,
        kind_type(QuantityKind::AbsoluteTemperature, ValueForm::Static),
        "construct absolute temperature from Celsius",
    )
}

/// Construct a static temperature difference from kelvin.
pub fn temperature_difference_kelvin(value: f64) -> Result<SemanticQty, SemanticError> {
    semantic_value(
        value,
        kind_type(QuantityKind::TemperatureDifference, ValueForm::Static),
        "construct temperature difference from kelvin",
    )
}

/// Construct a static temperature difference from degrees Celsius. Celsius
/// intervals have the same scale as kelvin and therefore receive no offset.
pub fn temperature_difference_celsius(value: f64) -> Result<SemanticQty, SemanticError> {
    semantic_value(
        value,
        kind_type(QuantityKind::TemperatureDifference, ValueForm::Static),
        "construct temperature difference from Celsius",
    )
}

/// Subtract two absolute temperatures to obtain a temperature difference.
pub fn temperature_difference(
    upper: SemanticQty,
    lower: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "subtract absolute temperatures";
    let upper = require_point_form(upper, operation)?;
    let lower = require_point_form(lower, operation)?;
    let absolute = kind_type(
        QuantityKind::AbsoluteTemperature,
        upper.semantic_type.form(),
    );
    let upper = require_type(upper, absolute, operation)?;
    let lower = require_type(lower, absolute, operation)?;
    semantic_value(
        upper.value() - lower.value(),
        kind_type(QuantityKind::TemperatureDifference, absolute.form()),
        operation,
    )
}

/// Add a temperature difference to an absolute temperature.
pub fn add_temperature_difference(
    absolute: SemanticQty,
    difference: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "add temperature difference";
    let absolute = require_point_form(absolute, operation)?;
    let difference = require_point_form(difference, operation)?;
    let absolute_type = kind_type(
        QuantityKind::AbsoluteTemperature,
        absolute.semantic_type.form(),
    );
    let absolute = require_type(absolute, absolute_type, operation)?;
    let difference = require_type(
        difference,
        kind_type(QuantityKind::TemperatureDifference, absolute_type.form()),
        operation,
    )?;
    semantic_value(
        absolute.value() + difference.value(),
        absolute_type,
        operation,
    )
}

/// Subtract a temperature difference from an absolute temperature.
pub fn subtract_temperature_difference(
    absolute: SemanticQty,
    difference: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "subtract temperature difference";
    let absolute = require_point_form(absolute, operation)?;
    let difference = require_point_form(difference, operation)?;
    let absolute_type = kind_type(
        QuantityKind::AbsoluteTemperature,
        absolute.semantic_type.form(),
    );
    let absolute = require_type(absolute, absolute_type, operation)?;
    let difference = require_type(
        difference,
        kind_type(QuantityKind::TemperatureDifference, absolute_type.form()),
        operation,
    )?;
    semantic_value(
        absolute.value() - difference.value(),
        absolute_type,
        operation,
    )
}

/// Add two temperature differences.
pub fn add_temperature_differences(
    left: SemanticQty,
    right: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    combine_temperature_differences(left, right, "add temperature differences", |a, b| a + b)
}

/// Subtract two temperature differences.
pub fn subtract_temperature_differences(
    left: SemanticQty,
    right: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    combine_temperature_differences(left, right, "subtract temperature differences", |a, b| {
        a - b
    })
}

fn combine_temperature_differences(
    left: SemanticQty,
    right: SemanticQty,
    operation: &'static str,
    combine: fn(f64, f64) -> f64,
) -> Result<SemanticQty, SemanticError> {
    let left = require_point_form(left, operation)?;
    let right = require_point_form(right, operation)?;
    let difference_type = kind_type(
        QuantityKind::TemperatureDifference,
        left.semantic_type.form(),
    );
    let left = require_type(left, difference_type, operation)?;
    let right = require_type(right, difference_type, operation)?;
    semantic_value(
        combine(left.value(), right.value()),
        difference_type,
        operation,
    )
}

/// Construct an instantaneous angle from radians.
pub fn angle_from_radians(radians: f64, domain: AngleDomain) -> Result<SemanticQty, SemanticError> {
    semantic_value(
        radians,
        kind_type(QuantityKind::Angle(domain), ValueForm::Instantaneous),
        "construct angle from radians",
    )
}

/// Construct an instantaneous angle from revolutions.
pub fn angle_from_revolutions(
    revolutions: f64,
    domain: AngleDomain,
) -> Result<SemanticQty, SemanticError> {
    angle_from_radians(revolutions * core::f64::consts::TAU, domain)
}

/// Read an angle in radians after checking its domain.
pub fn angle_to_radians(angle: SemanticQty, domain: AngleDomain) -> Result<f64, SemanticError> {
    let expected = kind_type(QuantityKind::Angle(domain), angle.semantic_type.form());
    Ok(require_type(angle, expected, "read angle in radians")?.value())
}

/// Convert an angle to revolutions after checking its domain.
pub fn angle_to_revolutions(angle: SemanticQty, domain: AngleDomain) -> Result<f64, SemanticError> {
    Ok(angle_to_radians(angle, domain)? / core::f64::consts::TAU)
}

/// Construct an instantaneous angular velocity from radians per second.
pub fn angular_velocity_from_radians_per_second(
    value: f64,
    domain: AngleDomain,
) -> Result<SemanticQty, SemanticError> {
    semantic_value(
        value,
        kind_type(
            QuantityKind::AngularVelocity(domain),
            ValueForm::Instantaneous,
        ),
        "construct angular velocity from radians per second",
    )
}

/// Construct an instantaneous angular velocity from revolutions per minute.
pub fn angular_velocity_from_rpm(
    rpm: f64,
    domain: AngleDomain,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert RPM to radians per second";
    let semantic_type = kind_type(
        QuantityKind::AngularVelocity(domain),
        ValueForm::Instantaneous,
    );
    let value = require_finite_result(rpm * RADIANS_PER_SECOND_PER_RPM, semantic_type, operation)?;
    semantic_value(value, semantic_type, operation)
}

/// Read an angular velocity in radians per second after checking its domain.
pub fn angular_velocity_to_radians_per_second(
    velocity: SemanticQty,
    domain: AngleDomain,
) -> Result<f64, SemanticError> {
    let expected = kind_type(
        QuantityKind::AngularVelocity(domain),
        velocity.semantic_type.form(),
    );
    Ok(require_type(
        velocity,
        expected,
        "read angular velocity in radians per second",
    )?
    .value())
}

/// Convert an angular velocity to revolutions per minute after checking its
/// domain.
pub fn angular_velocity_to_rpm(
    velocity: SemanticQty,
    domain: AngleDomain,
) -> Result<f64, SemanticError> {
    let operation = "convert radians per second to RPM";
    let semantic_type = kind_type(
        QuantityKind::AngularVelocity(domain),
        velocity.semantic_type.form(),
    );
    let radians_per_second = angular_velocity_to_radians_per_second(velocity, domain)?;
    require_finite_result(
        radians_per_second * RPM_PER_RADIAN_PER_SECOND,
        semantic_type,
        operation,
    )
}

/// Validated mechanical/electrical phase map
/// `theta_e = p * theta_m + electrical_phase_offset`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolePairPhaseMap {
    pole_pairs: NonZeroU32,
    electrical_phase_offset: f64,
}

impl PolePairPhaseMap {
    /// Construct a zero-offset map with a strictly positive pole-pair count.
    ///
    /// # Errors
    /// Returns [`SemanticError::ZeroPolePairs`] for zero.
    pub fn new(pole_pairs: u32) -> Result<Self, SemanticError> {
        Self::with_electrical_phase_offset(pole_pairs, 0.0)
    }

    /// Construct a map with a finite electrical phase offset in radians.
    ///
    /// The offset changes angle maps but not angular-velocity maps.
    ///
    /// # Errors
    /// Returns [`SemanticError::ZeroPolePairs`] for zero or
    /// [`SemanticError::InvalidValue`] for a non-finite offset.
    pub fn with_electrical_phase_offset(
        pole_pairs: u32,
        electrical_phase_offset: f64,
    ) -> Result<Self, SemanticError> {
        let pole_pairs = NonZeroU32::new(pole_pairs).ok_or(SemanticError::ZeroPolePairs)?;
        if !electrical_phase_offset.is_finite() {
            return Err(invalid_value(
                "construct pole-pair phase map",
                kind_type(
                    QuantityKind::Angle(AngleDomain::Electrical),
                    ValueForm::Instantaneous,
                ),
                electrical_phase_offset,
                ValueRequirement::Finite,
            ));
        }
        Ok(Self {
            pole_pairs,
            electrical_phase_offset,
        })
    }

    /// Pole-pair count.
    #[must_use]
    pub const fn pole_pairs(self) -> u32 {
        self.pole_pairs.get()
    }

    /// Electrical phase offset in radians.
    #[must_use]
    pub const fn electrical_phase_offset(self) -> f64 {
        self.electrical_phase_offset
    }

    /// Map a mechanical angle to electrical phase.
    pub fn mechanical_to_electrical_angle(
        self,
        source: SemanticQty,
    ) -> Result<SemanticQty, SemanticError> {
        self.map_angle(source, AngleDomain::Mechanical, AngleDomain::Electrical)
    }

    /// Map electrical phase back to mechanical angle.
    pub fn electrical_to_mechanical_angle(
        self,
        source: SemanticQty,
    ) -> Result<SemanticQty, SemanticError> {
        self.map_angle(source, AngleDomain::Electrical, AngleDomain::Mechanical)
    }

    /// Map mechanical angular velocity to electrical angular velocity.
    pub fn mechanical_to_electrical_angular_velocity(
        self,
        source: SemanticQty,
    ) -> Result<SemanticQty, SemanticError> {
        self.map_angular_velocity(source, AngleDomain::Mechanical, AngleDomain::Electrical)
    }

    /// Map electrical angular velocity back to mechanical angular velocity.
    pub fn electrical_to_mechanical_angular_velocity(
        self,
        source: SemanticQty,
    ) -> Result<SemanticQty, SemanticError> {
        self.map_angular_velocity(source, AngleDomain::Electrical, AngleDomain::Mechanical)
    }

    fn map_angle(
        self,
        source: SemanticQty,
        source_domain: AngleDomain,
        target_domain: AngleDomain,
    ) -> Result<SemanticQty, SemanticError> {
        let source = require_point_form(source, "apply pole-pair phase map")?;
        self.map_domain_quantity(source, source_domain, target_domain, true)
    }

    fn map_angular_velocity(
        self,
        source: SemanticQty,
        source_domain: AngleDomain,
        target_domain: AngleDomain,
    ) -> Result<SemanticQty, SemanticError> {
        self.map_domain_quantity(source, source_domain, target_domain, false)
    }

    fn map_domain_quantity(
        self,
        source: SemanticQty,
        source_domain: AngleDomain,
        target_domain: AngleDomain,
        angle: bool,
    ) -> Result<SemanticQty, SemanticError> {
        let operation = "apply pole-pair phase map";
        let form = source.semantic_type.form();
        let source_kind = if angle {
            QuantityKind::Angle(source_domain)
        } else {
            QuantityKind::AngularVelocity(source_domain)
        };
        let target_kind = if angle {
            QuantityKind::Angle(target_domain)
        } else {
            QuantityKind::AngularVelocity(target_domain)
        };
        let source = require_type(source, kind_type(source_kind, form), operation)?;
        let factor = f64::from(self.pole_pairs.get());
        let value = match (angle, source_domain) {
            (true, AngleDomain::Mechanical) => {
                let scaled = source.value() * factor;
                if scaled.is_finite() {
                    scaled + self.electrical_phase_offset
                } else {
                    // The product may overflow even when the opposite-signed
                    // phase offset makes the final affine result representable.
                    (source.value() + self.electrical_phase_offset / factor) * factor
                }
            }
            (true, AngleDomain::Electrical) => {
                let shifted = source.value() - self.electrical_phase_offset;
                if shifted.is_finite() {
                    shifted / factor
                } else {
                    // Divide before subtracting when the finite operands have
                    // opposite signs and their intermediate difference overflows.
                    source.value() / factor - self.electrical_phase_offset / factor
                }
            }
            (false, AngleDomain::Mechanical) => source.value() * factor,
            (false, AngleDomain::Electrical) => source.value() / factor,
        };
        semantic_value(value, kind_type(target_kind, form), operation)
    }
}

/// Convert between tensor and engineering strain conventions. The target type
/// must name the same component role and value form as the source.
pub fn convert_strain(
    source: SemanticQty,
    target: SemanticType,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert strain convention";
    let QuantityKind::Strain {
        basis: source_basis,
        component: source_component,
    } = source.semantic_type.kind()
    else {
        return Err(SemanticError::KindMismatch {
            operation,
            source: source.semantic_type,
            target,
        });
    };
    let QuantityKind::Strain {
        basis: target_basis,
        component: target_component,
    } = target.kind()
    else {
        return Err(SemanticError::KindMismatch {
            operation,
            source: source.semantic_type,
            target,
        });
    };
    if source_component != target_component || source.semantic_type.form() != target.form() {
        return Err(SemanticError::KindMismatch {
            operation,
            source: source.semantic_type,
            target,
        });
    }

    let value = match (source_component, source_basis, target_basis) {
        (_, StrainBasis::Tensor, StrainBasis::Tensor)
        | (_, StrainBasis::Engineering, StrainBasis::Engineering)
        | (StrainComponent::Normal, _, _) => source.value(),
        (StrainComponent::Shear, StrainBasis::Tensor, StrainBasis::Engineering) => {
            source.value() * 2.0
        }
        (StrainComponent::Shear, StrainBasis::Engineering, StrainBasis::Tensor) => {
            source.value() / 2.0
        }
    };
    semantic_value(value, target, operation)
}

const fn expected_with_source_form(source: SemanticQty, kind: QuantityKind) -> SemanticType {
    kind_type(kind, source.semantic_type.form())
}

const fn static_type(kind: QuantityKind) -> SemanticType {
    kind_type(kind, ValueForm::Static)
}

fn basis_conversion_result(
    source: SemanticQty,
    target: SemanticType,
    molar_mass: f64,
    result: f64,
    operation: &'static str,
) -> Result<SemanticQty, SemanticError> {
    if source.value() > 0.0 && result == 0.0 {
        return Err(SemanticError::UnrepresentablePositiveConversion {
            operation,
            source: source.semantic_type(),
            target,
            source_value: source.value(),
            molar_mass,
        });
    }
    semantic_value(result, target, operation)
}

/// Convert mass to amount of substance using a positive typed molar mass.
pub fn mass_to_amount(
    mass: SemanticQty,
    molar_mass: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert mass to amount";
    let mass_type = expected_with_source_form(mass, QuantityKind::Mass);
    let mass = require_type(mass, mass_type, operation)?;
    let molar_mass = require_type(molar_mass, static_type(QuantityKind::MolarMass), operation)?;
    basis_conversion_result(
        mass,
        kind_type(QuantityKind::Amount, mass_type.form()),
        molar_mass.value(),
        mass.value() / molar_mass.value(),
        operation,
    )
}

/// Convert amount of substance to mass using a positive typed molar mass.
pub fn amount_to_mass(
    amount: SemanticQty,
    molar_mass: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert amount to mass";
    let amount_type = expected_with_source_form(amount, QuantityKind::Amount);
    let amount = require_type(amount, amount_type, operation)?;
    let molar_mass = require_type(molar_mass, static_type(QuantityKind::MolarMass), operation)?;
    basis_conversion_result(
        amount,
        kind_type(QuantityKind::Mass, amount_type.form()),
        molar_mass.value(),
        amount.value() * molar_mass.value(),
        operation,
    )
}

/// Convert mass concentration to amount concentration using a positive typed
/// molar mass.
pub fn mass_concentration_to_amount_concentration(
    concentration: SemanticQty,
    molar_mass: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert mass concentration to amount concentration";
    let source_type = expected_with_source_form(concentration, QuantityKind::MassConcentration);
    let concentration = require_type(concentration, source_type, operation)?;
    let molar_mass = require_type(molar_mass, static_type(QuantityKind::MolarMass), operation)?;
    basis_conversion_result(
        concentration,
        kind_type(QuantityKind::AmountConcentration, source_type.form()),
        molar_mass.value(),
        concentration.value() / molar_mass.value(),
        operation,
    )
}

/// Convert amount concentration to mass concentration using a positive typed
/// molar mass.
pub fn amount_concentration_to_mass_concentration(
    concentration: SemanticQty,
    molar_mass: SemanticQty,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert amount concentration to mass concentration";
    let source_type = expected_with_source_form(concentration, QuantityKind::AmountConcentration);
    let concentration = require_type(concentration, source_type, operation)?;
    let molar_mass = require_type(molar_mass, static_type(QuantityKind::MolarMass), operation)?;
    basis_conversion_result(
        concentration,
        kind_type(QuantityKind::MassConcentration, source_type.form()),
        molar_mass.value(),
        concentration.value() * molar_mass.value(),
        operation,
    )
}

/// Immutable, normalized whole-composition vector. Whole-vector storage is
/// required because a single mass or mole fraction cannot change basis in
/// isolation.
#[derive(Debug, Clone, PartialEq)]
pub struct Composition {
    basis: CompositionBasis,
    fractions: Vec<f64>,
}

impl Composition {
    /// Validate and own a non-empty composition whose components are finite
    /// unit fractions summing to one within a fixed absolute tolerance.
    ///
    /// # Errors
    /// Returns a structured composition error; values are never normalized
    /// silently.
    pub fn new(basis: CompositionBasis, fractions: Vec<f64>) -> Result<Self, SemanticError> {
        if fractions.is_empty() {
            return Err(SemanticError::EmptyComposition { basis });
        }
        for (index, &value) in fractions.iter().enumerate() {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(SemanticError::InvalidCompositionEntry {
                    basis,
                    index,
                    value,
                });
            }
        }
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > COMPOSITION_SUM_TOLERANCE {
            return Err(SemanticError::InvalidCompositionSum {
                basis,
                sum,
                tolerance: COMPOSITION_SUM_TOLERANCE,
            });
        }
        Ok(Self { basis, fractions })
    }

    /// Composition basis.
    #[must_use]
    pub const fn basis(&self) -> CompositionBasis {
        self.basis
    }

    /// Borrow the immutable component fractions.
    #[must_use]
    pub fn fractions(&self) -> &[f64] {
        &self.fractions
    }

    /// Convert mass fractions to mole fractions with one typed molar mass per
    /// component. Calling this on a volume-fraction composition refuses.
    pub fn to_mole_fractions(&self, molar_masses: &[MolarMass]) -> Result<Self, SemanticError> {
        match self.basis {
            CompositionBasis::MoleFraction => Ok(self.clone()),
            CompositionBasis::MassFraction => self.convert_mass_mole(molar_masses, true),
            CompositionBasis::VolumeFraction => Err(composition_kind_mismatch(
                self.basis,
                CompositionBasis::MoleFraction,
                "convert composition to mole fractions",
            )),
        }
    }

    /// Convert mole fractions to mass fractions with one typed molar mass per
    /// component. Calling this on a volume-fraction composition refuses.
    pub fn to_mass_fractions(&self, molar_masses: &[MolarMass]) -> Result<Self, SemanticError> {
        match self.basis {
            CompositionBasis::MassFraction => Ok(self.clone()),
            CompositionBasis::MoleFraction => self.convert_mass_mole(molar_masses, false),
            CompositionBasis::VolumeFraction => Err(composition_kind_mismatch(
                self.basis,
                CompositionBasis::MassFraction,
                "convert composition to mass fractions",
            )),
        }
    }

    fn convert_mass_mole(
        &self,
        molar_masses: &[MolarMass],
        mass_to_mole: bool,
    ) -> Result<Self, SemanticError> {
        if self.fractions.len() != molar_masses.len() {
            return Err(SemanticError::CompositionLengthMismatch {
                fractions: self.fractions.len(),
                molar_masses: molar_masses.len(),
            });
        }
        for (index, molar_mass) in molar_masses.iter().enumerate() {
            let value = molar_mass.value();
            if !value.is_finite() || value <= 0.0 {
                return Err(SemanticError::InvalidMolarMass { index, value });
            }
        }

        let target_basis = if mass_to_mole {
            CompositionBasis::MoleFraction
        } else {
            CompositionBasis::MassFraction
        };
        let mut scale = if mass_to_mole { f64::INFINITY } else { 0.0 };
        for (&fraction, molar_mass) in self.fractions.iter().zip(molar_masses) {
            if fraction > 0.0 {
                scale = if mass_to_mole {
                    scale.min(molar_mass.value())
                } else {
                    scale.max(molar_mass.value())
                };
            }
        }

        let weights: Vec<f64> = self
            .fractions
            .iter()
            .zip(molar_masses)
            .map(|(&fraction, molar_mass)| {
                if fraction <= 0.0 {
                    0.0
                } else if mass_to_mole {
                    fraction * (scale / molar_mass.value())
                } else {
                    fraction * (molar_mass.value() / scale)
                }
            })
            .collect();
        let total: f64 = weights.iter().sum();
        if !total.is_finite() || total <= 0.0 || !weights.iter().any(|&weight| weight > 0.0) {
            return Err(SemanticError::InvalidCompositionSum {
                basis: target_basis,
                sum: total,
                tolerance: COMPOSITION_SUM_TOLERANCE,
            });
        }
        let fractions = weights.into_iter().map(|weight| weight / total).collect();
        Self::new(target_basis, fractions)
    }
}

const fn composition_kind_mismatch(
    source: CompositionBasis,
    target: CompositionBasis,
    operation: &'static str,
) -> SemanticError {
    SemanticError::KindMismatch {
        operation,
        source: static_type(QuantityKind::Composition(source)),
        target: static_type(QuantityKind::Composition(target)),
    }
}

/// Convert a declared sinusoidal peak amplitude to RMS. No conversion is
/// offered for an arbitrary waveform or instantaneous sample.
pub fn sinusoidal_peak_to_rms(source: SemanticQty) -> Result<SemanticQty, SemanticError> {
    convert_sinusoidal_amplitude(source, ValueForm::Peak, ValueForm::Rms)
}

/// Convert a declared sinusoidal RMS amplitude to peak. No conversion is
/// offered for an arbitrary waveform or instantaneous sample.
pub fn sinusoidal_rms_to_peak(source: SemanticQty) -> Result<SemanticQty, SemanticError> {
    convert_sinusoidal_amplitude(source, ValueForm::Rms, ValueForm::Peak)
}

fn convert_sinusoidal_amplitude(
    source: SemanticQty,
    source_form: ValueForm,
    target_form: ValueForm,
) -> Result<SemanticQty, SemanticError> {
    let operation = "convert sinusoidal peak/RMS amplitude";
    let source_type = kind_type(source.semantic_type.kind(), source_form);
    let target_type = kind_type(source.semantic_type.kind(), target_form);
    let source = require_type(source, source_type, operation)?;
    let value = if source_form == ValueForm::Peak {
        source.value() / core::f64::consts::SQRT_2
    } else {
        source.value() * core::f64::consts::SQRT_2
    };
    semantic_value(value, target_type, operation)
}

/// Amplitude convention carried by a complete complex phasor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhasorAmplitude {
    /// Complex components use peak amplitude.
    Peak,
    /// Complex components use RMS amplitude.
    Rms,
}

impl PhasorAmplitude {
    const fn scalar_form(self) -> ValueForm {
        match self {
            Self::Peak => ValueForm::Peak,
            Self::Rms => ValueForm::Rms,
        }
    }
}

/// Complete Cartesian complex phasor. Keeping real and imaginary components
/// together prevents independently tagged components from drifting in kind,
/// dimension, or peak/RMS convention.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhasorQty {
    real: QtyAny,
    imaginary: QtyAny,
    kind: QuantityKind,
    amplitude: PhasorAmplitude,
}

impl PhasorQty {
    /// Construct a finite phasor whose two components have the dimension
    /// required by `kind`. Cartesian components may be signed.
    ///
    /// # Errors
    /// Returns a structured dimension or finite-value error.
    pub fn new(
        real: QtyAny,
        imaginary: QtyAny,
        kind: QuantityKind,
        amplitude: PhasorAmplitude,
    ) -> Result<Self, SemanticError> {
        let semantic_type = kind_type(kind, amplitude.scalar_form());
        if !kind.admits_phasor() {
            return Err(SemanticError::UnsupportedForm {
                operation: "construct phasor",
                source: semantic_type,
                requirement: FormRequirement::StaticOnly,
            });
        }
        validate_phasor_component(real, semantic_type, "construct phasor real component")?;
        validate_phasor_component(
            imaginary,
            semantic_type,
            "construct phasor imaginary component",
        )?;
        Ok(Self {
            real,
            imaginary,
            kind,
            amplitude,
        })
    }

    /// Real component in coherent SI units.
    #[must_use]
    pub const fn real(self) -> QtyAny {
        self.real
    }

    /// Imaginary component in coherent SI units.
    #[must_use]
    pub const fn imaginary(self) -> QtyAny {
        self.imaginary
    }

    /// Physical quantity kind.
    #[must_use]
    pub const fn kind(self) -> QuantityKind {
        self.kind
    }

    /// Peak or RMS component convention.
    #[must_use]
    pub const fn amplitude(self) -> PhasorAmplitude {
        self.amplitude
    }

    /// Convert both Cartesian components between peak and RMS conventions.
    /// This conversion is valid because a phasor declares a harmonic signal.
    pub fn to_amplitude(self, target: PhasorAmplitude) -> Result<Self, SemanticError> {
        if self.amplitude == target {
            return Ok(self);
        }
        let factor = if self.amplitude == PhasorAmplitude::Peak {
            1.0 / core::f64::consts::SQRT_2
        } else {
            core::f64::consts::SQRT_2
        };
        Self::new(
            QtyAny::new(self.real.value * factor, self.real.dims),
            QtyAny::new(self.imaginary.value * factor, self.imaginary.dims),
            self.kind,
            target,
        )
    }
}

fn validate_phasor_component(
    component: QtyAny,
    semantic_type: SemanticType,
    operation: &'static str,
) -> Result<(), SemanticError> {
    let expected = semantic_type.expected_dims();
    if component.dims != expected {
        return Err(SemanticError::DimensionMismatch {
            operation,
            semantic_type,
            actual: component.dims,
            expected,
        });
    }
    if !component.value.is_finite() {
        return Err(invalid_value(
            operation,
            semantic_type,
            component.value,
            ValueRequirement::Finite,
        ));
    }
    Ok(())
}

/// Logarithmic acoustic level family selected by its physical reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AcousticLevelKind {
    /// Sound pressure level with a positive RMS pressure reference.
    Pressure,
    /// Sound power level with a positive power reference.
    Power,
}

/// Positive, explicitly typed reference carried by an acoustic level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticReference {
    quantity: SemanticQty,
    kind: AcousticLevelKind,
}

impl AcousticReference {
    /// Construct an acoustic-pressure reference. The value must be a positive
    /// RMS [`QuantityKind::AcousticPressure`].
    pub fn pressure(quantity: SemanticQty) -> Result<Self, SemanticError> {
        Self::new(
            quantity,
            kind_type(QuantityKind::AcousticPressure, ValueForm::Rms),
            AcousticLevelKind::Pressure,
            "construct acoustic pressure reference",
        )
    }

    /// Construct an acoustic-power reference. The value must be a positive
    /// static [`QuantityKind::AcousticPower`].
    pub fn power(quantity: SemanticQty) -> Result<Self, SemanticError> {
        Self::new(
            quantity,
            kind_type(QuantityKind::AcousticPower, ValueForm::Static),
            AcousticLevelKind::Power,
            "construct acoustic power reference",
        )
    }

    fn new(
        quantity: SemanticQty,
        target: SemanticType,
        kind: AcousticLevelKind,
        operation: &'static str,
    ) -> Result<Self, SemanticError> {
        let quantity = require_type(quantity, target, operation)?;
        if quantity.value() <= 0.0 {
            return Err(invalid_value(
                operation,
                target,
                quantity.value(),
                ValueRequirement::FinitePositive,
            ));
        }
        Ok(Self { quantity, kind })
    }

    /// Referenced physical quantity.
    #[must_use]
    pub const fn quantity(self) -> SemanticQty {
        self.quantity
    }

    /// Pressure-level or power-level family.
    #[must_use]
    pub const fn kind(self) -> AcousticLevelKind {
        self.kind
    }
}

/// Decibel value with a mandatory physical pressure or power reference.
///
/// This carrier intentionally does not compute `log10` or exponentials. The
/// owning acoustics/math layer performs those deterministic conversions and
/// then constructs this semantic record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticLevel {
    decibels: f64,
    reference: AcousticReference,
}

impl AcousticLevel {
    /// Construct a finite decibel level with an explicit positive reference.
    pub fn new(decibels: f64, reference: AcousticReference) -> Result<Self, SemanticError> {
        if !decibels.is_finite() {
            return Err(invalid_value(
                "construct acoustic level",
                reference.quantity.semantic_type,
                decibels,
                ValueRequirement::Finite,
            ));
        }
        Ok(Self {
            decibels,
            reference,
        })
    }

    /// Decibel value; negative levels are valid relative levels.
    #[must_use]
    pub const fn decibels(self) -> f64 {
        self.decibels
    }

    /// Explicit physical reference.
    #[must_use]
    pub const fn reference(self) -> AcousticReference {
        self.reference
    }

    /// Pressure-level or power-level family.
    #[must_use]
    pub const fn kind(self) -> AcousticLevelKind {
        self.reference.kind
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        assert!(
            (actual - expected).abs() <= 64.0 * f64::EPSILON * scale,
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    fn quantity(value: f64, kind: QuantityKind, form: ValueForm) -> SemanticQty {
        let semantic_type = SemanticType::new(kind, form);
        SemanticQty::new(QtyAny::new(value, kind.expected_dims()), semantic_type)
            .expect("valid semantic test quantity")
    }

    #[test]
    fn g0_private_carrier_checks_dimensions_and_scalar_domains() {
        let target = SemanticType::new(QuantityKind::Pressure, ValueForm::Static);
        let error = SemanticQty::new(QtyAny::dimensionless(1.0), target)
            .expect_err("dimension mismatch must refuse");
        assert!(matches!(
            error,
            SemanticError::DimensionMismatch {
                semantic_type,
                actual: Dims::NONE,
                expected: PRESSURE_DIMS,
                ..
            } if semantic_type == target
        ));

        assert!(absolute_temperature_kelvin(-1.0).is_err());
        assert!(quantity_result(f64::NAN, QuantityKind::Energy, ValueForm::Static).is_err());
        assert!(matches!(
            quantity_result(1.0, QuantityKind::Energy, ValueForm::Rms),
            Err(SemanticError::UnsupportedForm {
                source,
                requirement: FormRequirement::StaticOnly,
                ..
            }) if source == SemanticType::new(QuantityKind::Energy, ValueForm::Rms)
        ));
        assert!(quantity_result(-1.0, QuantityKind::AcousticPressure, ValueForm::Rms).is_err());
    }

    fn quantity_result(
        value: f64,
        kind: QuantityKind,
        form: ValueForm,
    ) -> Result<SemanticQty, SemanticError> {
        SemanticQty::new(
            QtyAny::new(value, kind.expected_dims()),
            SemanticType::new(kind, form),
        )
    }

    fn assert_kind_mismatch<T: core::fmt::Debug>(
        result: Result<T, SemanticError>,
        expected_source: SemanticType,
        expected_target: SemanticType,
    ) {
        match result {
            Err(SemanticError::KindMismatch { source, target, .. }) => {
                assert_eq!(source, expected_source);
                assert_eq!(target, expected_target);
            }
            other => panic!("expected typed kind mismatch, got {other:?}"),
        }
    }

    #[test]
    fn g0_affine_temperature_algebra_keeps_absolute_and_delta_distinct() {
        let cold = absolute_temperature_celsius(20.0).expect("20 C");
        let hot = absolute_temperature_celsius(30.0).expect("30 C");
        assert_close(cold.value(), 293.15);

        let delta = temperature_difference(hot, cold).expect("absolute difference");
        assert_eq!(
            delta.semantic_type().kind(),
            QuantityKind::TemperatureDifference
        );
        assert_close(delta.value(), 10.0);
        assert_close(
            temperature_difference_celsius(10.0)
                .expect("Celsius interval")
                .value(),
            10.0,
        );
        assert_close(
            add_temperature_difference(cold, delta)
                .expect("absolute plus delta")
                .value(),
            hot.value(),
        );
        assert_close(
            subtract_temperature_difference(hot, delta)
                .expect("absolute minus delta")
                .value(),
            cold.value(),
        );

        assert!(matches!(
            quantity_result(
                300.0,
                QuantityKind::AbsoluteTemperature,
                ValueForm::Peak
            ),
            Err(SemanticError::UnsupportedForm {
                source,
                requirement: FormRequirement::StaticOnly,
                ..
            }) if source.kind() == QuantityKind::AbsoluteTemperature
                && source.form() == ValueForm::Peak
        ));
        let peak_delta = quantity(10.0, QuantityKind::TemperatureDifference, ValueForm::Peak);
        for result in [
            add_temperature_difference(cold, peak_delta),
            subtract_temperature_difference(hot, peak_delta),
            add_temperature_differences(peak_delta, peak_delta),
            subtract_temperature_differences(peak_delta, peak_delta),
        ] {
            assert!(matches!(
                result,
                Err(SemanticError::UnsupportedForm {
                    source,
                    requirement: FormRequirement::PointValue,
                    ..
                }) if source.form() == ValueForm::Peak
            ));
        }
    }

    #[test]
    fn g0_revolution_and_rpm_conversions_are_named_and_invertible() {
        let angle = angle_from_revolutions(1.25, AngleDomain::Mechanical).expect("angle");
        assert_close(
            angle_to_radians(angle, AngleDomain::Mechanical).expect("radians"),
            1.25 * core::f64::consts::TAU,
        );
        assert_close(
            angle_to_revolutions(angle, AngleDomain::Mechanical).expect("revolutions"),
            1.25,
        );

        let velocity =
            angular_velocity_from_rpm(120.0, AngleDomain::Mechanical).expect("angular velocity");
        assert_close(
            angular_velocity_to_radians_per_second(velocity, AngleDomain::Mechanical)
                .expect("rad/s"),
            2.0 * core::f64::consts::TAU,
        );
        assert_close(
            angular_velocity_to_rpm(velocity, AngleDomain::Mechanical).expect("rpm"),
            120.0,
        );

        let extreme_rpm = f64::MAX / 2.0;
        let extreme_velocity = angular_velocity_from_rpm(extreme_rpm, AngleDomain::Mechanical)
            .expect("precomputed RPM ratio avoids intermediate overflow");
        assert!(extreme_velocity.value().is_finite());
        assert_close(
            extreme_velocity.value(),
            extreme_rpm * RADIANS_PER_SECOND_PER_RPM,
        );

        let extreme_radians_per_second = f64::MAX / 32.0;
        let extreme_velocity = angular_velocity_from_radians_per_second(
            extreme_radians_per_second,
            AngleDomain::Mechanical,
        )
        .expect("finite extreme angular velocity");
        let extreme_rpm = angular_velocity_to_rpm(extreme_velocity, AngleDomain::Mechanical)
            .expect("precomputed inverse ratio avoids intermediate overflow");
        assert!(extreme_rpm.is_finite());
        assert_close(
            extreme_rpm,
            extreme_radians_per_second * RPM_PER_RADIAN_PER_SECOND,
        );

        let overflowing_rpm =
            angular_velocity_from_radians_per_second(f64::MAX, AngleDomain::Mechanical)
                .expect("maximum finite angular velocity");
        assert!(matches!(
            angular_velocity_to_rpm(overflowing_rpm, AngleDomain::Mechanical),
            Err(SemanticError::InvalidValue {
                requirement: ValueRequirement::Finite,
                ..
            })
        ));
        assert!(angular_velocity_from_rpm(f64::INFINITY, AngleDomain::Mechanical).is_err());
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One fail-closed pole-pair domain and refusal matrix.
    fn g0_pole_pair_map_changes_domain_and_rejects_wrong_source() {
        assert!(matches!(
            PolePairPhaseMap::new(0),
            Err(SemanticError::ZeroPolePairs)
        ));
        assert_close(
            PolePairPhaseMap::new(1)
                .expect("default phase map")
                .electrical_phase_offset(),
            0.0,
        );
        assert!(matches!(
            PolePairPhaseMap::with_electrical_phase_offset(3, f64::NAN),
            Err(SemanticError::InvalidValue {
                semantic_type,
                requirement: ValueRequirement::Finite,
                ..
            }) if semantic_type.kind()
                == QuantityKind::Angle(AngleDomain::Electrical)
        ));
        let map = PolePairPhaseMap::with_electrical_phase_offset(3, 0.25)
            .expect("three pole pairs with phase offset");
        assert_close(map.electrical_phase_offset(), 0.25);
        let mechanical = angle_from_radians(0.5, AngleDomain::Mechanical).expect("mechanical");
        let electrical = map
            .mechanical_to_electrical_angle(mechanical)
            .expect("electrical");
        assert_eq!(
            electrical.semantic_type().kind(),
            QuantityKind::Angle(AngleDomain::Electrical)
        );
        assert_close(electrical.value(), 1.75);
        assert_close(
            map.electrical_to_mechanical_angle(electrical)
                .expect("mechanical round trip")
                .value(),
            mechanical.value(),
        );

        let mechanical_velocity =
            angular_velocity_from_radians_per_second(2.0, AngleDomain::Mechanical)
                .expect("mechanical angular velocity");
        let electrical_velocity = map
            .mechanical_to_electrical_angular_velocity(mechanical_velocity)
            .expect("electrical angular velocity");
        assert_close(electrical_velocity.value(), 6.0);
        assert_close(
            map.electrical_to_mechanical_angular_velocity(electrical_velocity)
                .expect("mechanical angular-velocity round trip")
                .value(),
            mechanical_velocity.value(),
        );

        let peak_angle = quantity(
            0.5,
            QuantityKind::Angle(AngleDomain::Mechanical),
            ValueForm::Peak,
        );
        assert!(matches!(
            map.mechanical_to_electrical_angle(peak_angle),
            Err(SemanticError::UnsupportedForm {
                source,
                requirement: FormRequirement::PointValue,
                ..
            }) if source.form() == ValueForm::Peak
        ));
        let peak_velocity = quantity(
            2.0,
            QuantityKind::AngularVelocity(AngleDomain::Mechanical),
            ValueForm::Peak,
        );
        let electrical_peak_velocity = map
            .mechanical_to_electrical_angular_velocity(peak_velocity)
            .expect("phase offset does not affect angular-velocity amplitude");
        assert_eq!(
            electrical_peak_velocity.semantic_type().form(),
            ValueForm::Peak
        );
        assert_close(electrical_peak_velocity.value(), 6.0);

        let extreme_forward = PolePairPhaseMap::with_electrical_phase_offset(3, -f64::MAX)
            .expect("finite extreme phase offset")
            .mechanical_to_electrical_angle(
                angle_from_radians(f64::MAX / 2.0, AngleDomain::Mechanical)
                    .expect("finite extreme mechanical angle"),
            )
            .expect("scaled affine form avoids intermediate product overflow");
        assert!(extreme_forward.value().is_finite());
        assert_close(extreme_forward.value(), f64::MAX / 2.0);

        let extreme_inverse = PolePairPhaseMap::with_electrical_phase_offset(3, -f64::MAX)
            .expect("finite extreme phase offset")
            .electrical_to_mechanical_angle(
                angle_from_radians(f64::MAX, AngleDomain::Electrical)
                    .expect("finite extreme electrical angle"),
            )
            .expect("divide-first affine form avoids intermediate difference overflow");
        assert!(extreme_inverse.value().is_finite());
        assert_close(extreme_inverse.value(), 2.0 * (f64::MAX / 3.0));

        let wrong = quantity(1.0, QuantityKind::Energy, ValueForm::Static);
        assert!(matches!(
            map.mechanical_to_electrical_angle(wrong),
            Err(SemanticError::KindMismatch { source, target, .. })
                if source.kind() == QuantityKind::Energy
                    && target.kind() == QuantityKind::Angle(AngleDomain::Mechanical)
        ));
    }

    #[test]
    fn g0_same_dimension_semantics_do_not_interchange() {
        let energy = quantity(4.0, QuantityKind::Energy, ValueForm::Static);
        let torque_type = SemanticType::new(QuantityKind::Torque, ValueForm::Static);
        assert!(matches!(
            energy.require_type(torque_type, "torque boundary"),
            Err(SemanticError::KindMismatch { source, target, .. })
                if source.kind() == QuantityKind::Energy
                    && target.kind() == QuantityKind::Torque
        ));

        let stress = quantity(2.0, QuantityKind::Stress, ValueForm::Static);
        let pressure_type = SemanticType::new(QuantityKind::Pressure, ValueForm::Static);
        assert_kind_mismatch(
            stress.require_type(pressure_type, "pressure boundary"),
            stress.semantic_type(),
            pressure_type,
        );

        let entropy = quantity(1.0, QuantityKind::Entropy, ValueForm::Static);
        let heat_capacity = SemanticType::new(QuantityKind::HeatCapacity, ValueForm::Static);
        assert_kind_mismatch(
            entropy.require_type(heat_capacity, "heat-capacity boundary"),
            entropy.semantic_type(),
            heat_capacity,
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One crossed kind/form/domain Gauntlet battery.
    fn g0_crossed_kind_form_and_domain_battery_fails_closed() {
        let absolute = absolute_temperature_kelvin(300.0).expect("absolute temperature");
        let difference = temperature_difference_kelvin(5.0).expect("temperature difference");
        assert!(matches!(
            temperature_difference(difference, absolute),
            Err(SemanticError::KindMismatch { source, target, .. })
                if source.kind() == QuantityKind::TemperatureDifference
                    && target.kind() == QuantityKind::AbsoluteTemperature
        ));
        assert!(absolute_temperature_kelvin(f64::NAN).is_err());
        assert!(temperature_difference_kelvin(f64::INFINITY).is_err());

        let mechanical_angle =
            angle_from_radians(1.0, AngleDomain::Mechanical).expect("mechanical angle");
        assert!(matches!(
            angle_to_radians(mechanical_angle, AngleDomain::Electrical),
            Err(SemanticError::KindMismatch { source, target, .. })
                if source.kind() == QuantityKind::Angle(AngleDomain::Mechanical)
                    && target.kind() == QuantityKind::Angle(AngleDomain::Electrical)
        ));
        let electrical_velocity =
            angular_velocity_from_radians_per_second(1.0, AngleDomain::Electrical)
                .expect("electrical angular velocity");
        assert_kind_mismatch(
            angular_velocity_to_rpm(electrical_velocity, AngleDomain::Mechanical),
            electrical_velocity.semantic_type(),
            SemanticType::new(
                QuantityKind::AngularVelocity(AngleDomain::Mechanical),
                ValueForm::Instantaneous,
            ),
        );
        assert!(angle_from_radians(f64::NEG_INFINITY, AngleDomain::Mechanical).is_err());
        assert!(
            angular_velocity_from_radians_per_second(f64::NAN, AngleDomain::Electrical).is_err()
        );

        let tensor_shear = quantity(
            0.01,
            QuantityKind::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Shear,
            },
            ValueForm::Static,
        );
        for target in [
            SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Engineering,
                    component: StrainComponent::Normal,
                },
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Engineering,
                    component: StrainComponent::Shear,
                },
                ValueForm::Instantaneous,
            ),
        ] {
            assert!(matches!(
                convert_strain(tensor_shear, target),
                Err(SemanticError::KindMismatch {
                    source,
                    target: actual_target,
                    ..
                }) if source == tensor_shear.semantic_type() && actual_target == target
            ));
        }

        let molar_mass = quantity(0.018, QuantityKind::MolarMass, ValueForm::Static);
        let mass = quantity(1.0, QuantityKind::Mass, ValueForm::Static);
        let amount = quantity(1.0, QuantityKind::Amount, ValueForm::Static);
        assert_kind_mismatch(
            mass_to_amount(amount, molar_mass),
            amount.semantic_type(),
            SemanticType::new(QuantityKind::Mass, ValueForm::Static),
        );
        assert_kind_mismatch(
            mass_to_amount(mass, amount),
            amount.semantic_type(),
            molar_mass.semantic_type(),
        );
        let mass_concentration = quantity(1.0, QuantityKind::MassConcentration, ValueForm::Static);
        let amount_concentration =
            quantity(1.0, QuantityKind::AmountConcentration, ValueForm::Static);
        assert_kind_mismatch(
            mass_concentration_to_amount_concentration(amount_concentration, molar_mass),
            amount_concentration.semantic_type(),
            mass_concentration.semantic_type(),
        );
        assert_kind_mismatch(
            amount_concentration_to_mass_concentration(mass_concentration, molar_mass),
            mass_concentration.semantic_type(),
            amount_concentration.semantic_type(),
        );
        assert!(quantity_result(-1.0, QuantityKind::Mass, ValueForm::Static).is_err());
        assert!(quantity_result(0.0, QuantityKind::MolarMass, ValueForm::Static).is_err());
        assert!(quantity_result(-1.0, QuantityKind::MassConcentration, ValueForm::Static).is_err());
        assert!(
            quantity_result(
                f64::INFINITY,
                QuantityKind::AmountConcentration,
                ValueForm::Static
            )
            .is_err()
        );
    }

    #[test]
    fn g0_strain_conversion_doubles_only_shear() {
        let tensor_shear = quantity(
            0.01,
            QuantityKind::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Shear,
            },
            ValueForm::Static,
        );
        let engineering_shear = SemanticType::new(
            QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Shear,
            },
            ValueForm::Static,
        );
        assert_close(
            convert_strain(tensor_shear, engineering_shear)
                .expect("shear conversion")
                .value(),
            0.02,
        );

        let tensor_normal = quantity(
            0.01,
            QuantityKind::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Normal,
            },
            ValueForm::Static,
        );
        let engineering_normal = SemanticType::new(
            QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Normal,
            },
            ValueForm::Static,
        );
        assert_close(
            convert_strain(tensor_normal, engineering_normal)
                .expect("normal conversion")
                .value(),
            0.01,
        );
    }

    #[test]
    fn g0_mass_amount_and_concentration_conversions_are_typed() {
        let molar_mass = quantity(0.018, QuantityKind::MolarMass, ValueForm::Static);
        let mass = quantity(18.0, QuantityKind::Mass, ValueForm::Static);
        let amount = mass_to_amount(mass, molar_mass).expect("mass to amount");
        assert_eq!(amount.semantic_type().kind(), QuantityKind::Amount);
        assert_close(amount.value(), 1000.0);
        assert_close(
            amount_to_mass(amount, molar_mass)
                .expect("amount to mass")
                .value(),
            mass.value(),
        );

        let mass_concentration = quantity(18.0, QuantityKind::MassConcentration, ValueForm::Static);
        let amount_concentration =
            mass_concentration_to_amount_concentration(mass_concentration, molar_mass)
                .expect("concentration basis");
        assert_close(amount_concentration.value(), 1000.0);
        assert_close(
            amount_concentration_to_mass_concentration(amount_concentration, molar_mass)
                .expect("concentration round trip")
                .value(),
            mass_concentration.value(),
        );
    }

    #[test]
    fn g0_positive_basis_conversions_refuse_zero_underflow_at_f64_boundary() {
        let minimum = f64::from_bits(1);
        let next = f64::from_bits(2);
        let twice = quantity(2.0, QuantityKind::MolarMass, ValueForm::Static);
        let half = quantity(0.5, QuantityKind::MolarMass, ValueForm::Static);

        for result in [
            mass_to_amount(
                quantity(minimum, QuantityKind::Mass, ValueForm::Static),
                twice,
            ),
            amount_to_mass(
                quantity(minimum, QuantityKind::Amount, ValueForm::Static),
                half,
            ),
            mass_concentration_to_amount_concentration(
                quantity(minimum, QuantityKind::MassConcentration, ValueForm::Static),
                twice,
            ),
            amount_concentration_to_mass_concentration(
                quantity(
                    minimum,
                    QuantityKind::AmountConcentration,
                    ValueForm::Static,
                ),
                half,
            ),
        ] {
            assert!(matches!(
                result,
                Err(SemanticError::UnrepresentablePositiveConversion {
                    source_value,
                    molar_mass,
                    ..
                }) if source_value.to_bits() == minimum.to_bits()
                    && (molar_mass.to_bits() == 0.5f64.to_bits()
                        || molar_mass.to_bits() == 2.0f64.to_bits())
            ));
        }

        let amount = mass_to_amount(quantity(next, QuantityKind::Mass, ValueForm::Static), twice)
            .expect("the adjacent subnormal remains representable");
        assert_eq!(amount.value().to_bits(), minimum.to_bits());
        assert_eq!(
            amount_to_mass(amount, twice)
                .expect("boundary round trip")
                .value()
                .to_bits(),
            next.to_bits()
        );

        let mass = amount_to_mass(
            quantity(next, QuantityKind::Amount, ValueForm::Static),
            half,
        )
        .expect("the adjacent subnormal remains representable");
        assert_eq!(mass.value().to_bits(), minimum.to_bits());
        assert_eq!(
            mass_to_amount(mass, half)
                .expect("boundary round trip")
                .value()
                .to_bits(),
            next.to_bits()
        );

        let zero_mass = quantity(0.0, QuantityKind::Mass, ValueForm::Static);
        assert_eq!(
            mass_to_amount(zero_mass, twice)
                .expect("true zero remains representable")
                .value()
                .to_bits(),
            0.0f64.to_bits()
        );
        let maximum_amount = quantity(f64::MAX, QuantityKind::Amount, ValueForm::Static);
        assert!(matches!(
            amount_to_mass(maximum_amount, twice),
            Err(SemanticError::InvalidValue {
                requirement: ValueRequirement::Finite,
                ..
            })
        ));
    }

    #[test]
    fn g0_composition_conversion_is_whole_vector_and_refuses_volume_shortcut() {
        let mass = Composition::new(CompositionBasis::MassFraction, vec![0.5, 0.5])
            .expect("mass fractions");
        let molar_masses = [MolarMass::new(2.0), MolarMass::new(4.0)];
        let mole = mass
            .to_mole_fractions(&molar_masses)
            .expect("mole fractions");
        assert_close(mole.fractions()[0], 2.0 / 3.0);
        assert_close(mole.fractions()[1], 1.0 / 3.0);
        let round_trip = mole
            .to_mass_fractions(&molar_masses)
            .expect("mass-fraction round trip");
        assert_close(round_trip.fractions()[0], 0.5);
        assert_close(round_trip.fractions()[1], 0.5);

        assert!(Composition::new(CompositionBasis::MassFraction, vec![0.4, 0.4]).is_err());
        assert!(mass.to_mole_fractions(&molar_masses[..1]).is_err());
        let volume =
            Composition::new(CompositionBasis::VolumeFraction, vec![1.0]).expect("volume fraction");
        assert!(matches!(
            volume.to_mass_fractions(&[MolarMass::new(1.0)]),
            Err(SemanticError::KindMismatch { source, target, .. })
                if source.kind() == QuantityKind::Composition(CompositionBasis::VolumeFraction)
                    && target.kind()
                        == QuantityKind::Composition(CompositionBasis::MassFraction)
        ));

        let extreme_molar_masses = [MolarMass::new(f64::from_bits(1)), MolarMass::new(f64::MAX)];
        let equal_mass = Composition::new(CompositionBasis::MassFraction, vec![0.5, 0.5])
            .expect("equal extreme mass fractions");
        let extreme_mole = equal_mass
            .to_mole_fractions(&extreme_molar_masses)
            .expect("minimum-active-mass scaling avoids reciprocal underflow");
        assert_close(extreme_mole.fractions()[0], 1.0);
        assert_eq!(extreme_mole.fractions()[1].to_bits(), 0.0f64.to_bits());

        let equal_mole = Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5])
            .expect("equal extreme mole fractions");
        let extreme_mass = equal_mole
            .to_mass_fractions(&extreme_molar_masses)
            .expect("maximum-active-mass scaling avoids product overflow");
        assert_eq!(extreme_mass.fractions()[0].to_bits(), 0.0f64.to_bits());
        assert_close(extreme_mass.fractions()[1], 1.0);

        let subnormal_molar_masses = [
            MolarMass::new(f64::from_bits(1)),
            MolarMass::new(f64::from_bits(1)),
        ];
        let equal_mole_from_subnormal = equal_mass
            .to_mole_fractions(&subnormal_molar_masses)
            .expect("scaled reciprocal avoids all-infinite weights");
        assert_close(equal_mole_from_subnormal.fractions()[0], 0.5);
        assert_close(equal_mole_from_subnormal.fractions()[1], 0.5);
        let equal_mass_from_subnormal = equal_mole
            .to_mass_fractions(&subnormal_molar_masses)
            .expect("scaled product avoids all-zero weights");
        assert_close(equal_mass_from_subnormal.fractions()[0], 0.5);
        assert_close(equal_mass_from_subnormal.fractions()[1], 0.5);

        let zero_preserving = Composition::new(CompositionBasis::MassFraction, vec![0.0, 1.0])
            .expect("zero component");
        let converted = zero_preserving
            .to_mole_fractions(&extreme_molar_masses)
            .expect("inactive zero does not set the scaling anchor");
        assert_eq!(converted.fractions()[0].to_bits(), 0.0f64.to_bits());
        assert_close(converted.fractions()[1], 1.0);
    }

    #[test]
    fn g0_composition_negative_and_boundary_battery_fails_closed() {
        assert!(matches!(
            Composition::new(CompositionBasis::MassFraction, Vec::new()),
            Err(SemanticError::EmptyComposition {
                basis: CompositionBasis::MassFraction
            })
        ));
        for fractions in [
            vec![-0.1, 1.1],
            vec![1.1, -0.1],
            vec![f64::NAN, 1.0],
            vec![f64::INFINITY, 0.0],
            vec![0.25, 0.25],
        ] {
            assert!(Composition::new(CompositionBasis::MoleFraction, fractions).is_err());
        }

        let composition = Composition::new(CompositionBasis::MassFraction, vec![0.5, 0.5])
            .expect("valid boundary fixture");
        assert!(matches!(
            composition.to_mole_fractions(&[MolarMass::new(1.0)]),
            Err(SemanticError::CompositionLengthMismatch {
                fractions: 2,
                molar_masses: 1
            })
        ));
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                composition.to_mole_fractions(&[
                    MolarMass::new(1.0),
                    MolarMass::new(invalid),
                ]),
                Err(SemanticError::InvalidMolarMass { index: 1, value })
                    if value.to_bits() == invalid.to_bits()
            ));
        }
    }

    #[test]
    fn g0_sinusoidal_amplitude_conversion_requires_declared_form() {
        let peak = quantity(10.0, QuantityKind::AcousticPressure, ValueForm::Peak);
        let rms = sinusoidal_peak_to_rms(peak).expect("peak to RMS");
        assert_eq!(rms.semantic_type().form(), ValueForm::Rms);
        assert_close(rms.value(), 10.0 / core::f64::consts::SQRT_2);
        assert_close(
            sinusoidal_rms_to_peak(rms).expect("RMS to peak").value(),
            10.0,
        );

        let instantaneous = quantity(
            10.0,
            QuantityKind::AcousticPressure,
            ValueForm::Instantaneous,
        );
        assert_kind_mismatch(
            sinusoidal_peak_to_rms(instantaneous),
            instantaneous.semantic_type(),
            SemanticType::new(QuantityKind::AcousticPressure, ValueForm::Peak),
        );
        assert_kind_mismatch(
            sinusoidal_rms_to_peak(peak),
            peak.semantic_type(),
            SemanticType::new(QuantityKind::AcousticPressure, ValueForm::Rms),
        );
        assert!(
            quantity_result(
                -f64::MIN_POSITIVE,
                QuantityKind::AcousticPressure,
                ValueForm::Peak
            )
            .is_err()
        );
        assert!(matches!(
            quantity_result(1.0, QuantityKind::AcousticPower, ValueForm::Rms),
            Err(SemanticError::UnsupportedForm {
                source,
                requirement: FormRequirement::StaticOnly,
                ..
            }) if source.kind() == QuantityKind::AcousticPower
                && source.form() == ValueForm::Rms
        ));
    }

    #[test]
    fn g0_phasor_keeps_components_paired_and_converts_convention() {
        let dims = QuantityKind::AcousticPressure.expected_dims();
        let peak = PhasorQty::new(
            QtyAny::new(-2.0, dims),
            QtyAny::new(4.0, dims),
            QuantityKind::AcousticPressure,
            PhasorAmplitude::Peak,
        )
        .expect("signed Cartesian components are valid");
        let rms = peak.to_amplitude(PhasorAmplitude::Rms).expect("phasor RMS");
        assert_close(rms.real().value, -2.0 / core::f64::consts::SQRT_2);
        assert_close(rms.imaginary().value, 4.0 / core::f64::consts::SQRT_2);
        assert!(
            PhasorQty::new(
                QtyAny::new(1.0, dims),
                QtyAny::dimensionless(2.0),
                QuantityKind::AcousticPressure,
                PhasorAmplitude::Peak,
            )
            .is_err()
        );
        assert!(
            PhasorQty::new(
                QtyAny::new(f64::NAN, dims),
                QtyAny::new(0.0, dims),
                QuantityKind::AcousticPressure,
                PhasorAmplitude::Rms,
            )
            .is_err()
        );
        let overflowing = PhasorQty::new(
            QtyAny::new(f64::MAX, dims),
            QtyAny::new(0.0, dims),
            QuantityKind::AcousticPressure,
            PhasorAmplitude::Rms,
        )
        .expect("maximum finite RMS component");
        assert!(overflowing.to_amplitude(PhasorAmplitude::Peak).is_err());
        let mass_dims = QuantityKind::Mass.expected_dims();
        assert!(matches!(
            PhasorQty::new(
                QtyAny::new(1.0, mass_dims),
                QtyAny::new(0.0, mass_dims),
                QuantityKind::Mass,
                PhasorAmplitude::Peak,
            ),
            Err(SemanticError::UnsupportedForm {
                source,
                requirement: FormRequirement::StaticOnly,
                ..
            }) if source.kind() == QuantityKind::Mass
                && source.form() == ValueForm::Peak
        ));
    }

    #[test]
    fn g0_acoustic_level_requires_typed_positive_reference_without_raw_log_math() {
        let pressure_reference = quantity(20.0e-6, QuantityKind::AcousticPressure, ValueForm::Rms);
        let reference =
            AcousticReference::pressure(pressure_reference).expect("pressure reference");
        let level = AcousticLevel::new(-3.0, reference).expect("negative relative dB is valid");
        assert_eq!(level.kind(), AcousticLevelKind::Pressure);
        assert_close(level.decibels(), -3.0);
        assert_close(level.reference().quantity().value(), 20.0e-6);

        let zero_reference = quantity(0.0, QuantityKind::AcousticPressure, ValueForm::Rms);
        assert!(AcousticReference::pressure(zero_reference).is_err());
        let stress = quantity(1.0, QuantityKind::Stress, ValueForm::Rms);
        assert!(matches!(
            AcousticReference::pressure(stress),
            Err(SemanticError::KindMismatch { source, target, .. })
                if source.kind() == QuantityKind::Stress
                    && target.kind() == QuantityKind::AcousticPressure
        ));

        let power_reference = quantity(1.0e-12, QuantityKind::AcousticPower, ValueForm::Static);
        let power_reference =
            AcousticReference::power(power_reference).expect("sound-power reference");
        let power_level = AcousticLevel::new(0.0, power_reference).expect("power level");
        assert_eq!(power_level.kind(), AcousticLevelKind::Power);
        assert_close(power_level.reference().quantity().value(), 1.0e-12);

        assert!(matches!(
            quantity_result(1.0, QuantityKind::AcousticPower, ValueForm::Rms),
            Err(SemanticError::UnsupportedForm {
                source,
                requirement: FormRequirement::StaticOnly,
                ..
            }) if source.kind() == QuantityKind::AcousticPower
                && source.form() == ValueForm::Rms
        ));
        let zero_power = quantity(0.0, QuantityKind::AcousticPower, ValueForm::Static);
        assert!(AcousticReference::power(zero_power).is_err());
        for invalid_level in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                AcousticLevel::new(invalid_level, power_reference),
                Err(SemanticError::InvalidValue {
                    requirement: ValueRequirement::Finite,
                    ..
                })
            ));
        }
    }
}
