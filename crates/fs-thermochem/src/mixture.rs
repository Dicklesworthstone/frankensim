//! Frozen-composition ideal-gas mixture thermodynamics.
//!
//! This module combines admitted species standard states without owning a
//! reacting derivative, equilibrium solve, phase-stability decision, transport
//! rule, or evolving gas state. Composition basis is the shared `fs-qty`
//! authority; components are canonicalized by `SpeciesId`, and every successful
//! evaluation retains the complete component receipts.

use core::fmt;

use fs_qty::{MolarMass, Pressure, Temperature};

pub use crate::{MassSpecificGasConstantV1, MassSpecificThermalQuantityV1};

use crate::{
    Composition, CompositionBasis, ElementalReferenceIdV1, MolarEnthalpyV1, MolarEntropyV1,
    MolarGibbsEnergyV1, MolarHeatCapacityV1, MolarInternalEnergyV1, MolarThermalQuantityV1,
    Nasa9EvaluationReceiptV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SemanticError,
    SpeciesId, StandardStatePhaseV1, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};

/// Version of the frozen ideal-gas mixture operation tree and receipt layout.
pub const FROZEN_IDEAL_GAS_MIXTURE_EVALUATOR_VERSION_V1: u32 = 1;
/// Hard resource bound checked before any declared component is otherwise admitted.
pub const MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1: usize = 128;

macro_rules! thermal_quantity {
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

thermal_quantity!(
    /// Frozen-composition molar heat capacity at constant volume.
    MolarIsochoricHeatCapacityV1,
    MolarThermalQuantityV1
);
thermal_quantity!(
    /// Frozen-composition mass-specific heat capacity at constant pressure.
    MassSpecificIsobaricHeatCapacityV1,
    MassSpecificThermalQuantityV1
);
thermal_quantity!(
    /// Frozen-composition mass-specific heat capacity at constant volume.
    MassSpecificIsochoricHeatCapacityV1,
    MassSpecificThermalQuantityV1
);
/// Convention field named by a cross-component mismatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FrozenMixtureConventionFieldV1 {
    /// Standard-state phase.
    Phase,
    /// Reference equation of state.
    EquationOfState,
    /// Reference-pressure IEEE bits.
    ReferencePressure,
    /// Elemental-reference convention id.
    ElementalReference,
}

/// Exact convention value retained in an actionable mismatch diagnostic.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FrozenMixtureConventionValueV1 {
    /// Standard-state phase.
    Phase(StandardStatePhaseV1),
    /// Reference equation of state.
    EquationOfState(ReferenceEquationOfStateV1),
    /// Reference-pressure IEEE bits.
    ReferencePressureBits(u64),
    /// Elemental-reference convention id.
    ElementalReference(ElementalReferenceIdV1),
}

/// Derived mixture field named by a non-finite refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FrozenMixturePropertyV1 {
    /// Molar heat capacity at constant pressure.
    MolarCp,
    /// Molar heat capacity at constant volume.
    MolarCv,
    /// Molar enthalpy.
    MolarEnthalpy,
    /// Molar internal energy.
    MolarInternalEnergy,
    /// Molar entropy.
    MolarEntropy,
    /// Molar Gibbs energy.
    MolarGibbsEnergy,
    /// Mixture molar mass.
    MolarMass,
    /// Mass-specific heat capacity at constant pressure.
    MassSpecificCp,
    /// Mass-specific heat capacity at constant volume.
    MassSpecificCv,
    /// Mixture-specific gas constant.
    MassSpecificGasConstant,
}

/// Frozen ideal-gas mixture construction or evaluation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum FrozenIdealGasMixtureErrorV1 {
    /// At least one component model is required.
    EmptyComponents,
    /// The declared component count exceeded the fixed resource bound.
    TooManyComponents {
        /// Offered component count.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// Component models and composition fractions have different lengths.
    CompositionLengthMismatch {
        /// Model count.
        components: usize,
        /// Fraction count.
        fractions: usize,
    },
    /// The shared `fs-qty` composition boundary refused the vector or basis.
    Composition(SemanticError),
    /// A listed component has zero or negative-zero fraction and must be omitted.
    ZeroFraction {
        /// Canonical component index.
        component: usize,
        /// Affected species.
        species: SpeciesId,
        /// Exact zero IEEE bits.
        bits: u64,
    },
    /// Canonical fixed-order fractions do not sum to exactly one.
    FractionSumNotExact {
        /// Fraction basis whose canonical sum failed.
        basis: CompositionBasis,
        /// Exact canonical-sum IEEE bits.
        sum_bits: u64,
    },
    /// Two active entries name the same canonical species.
    DuplicateSpecies {
        /// Duplicated species id.
        species: SpeciesId,
    },
    /// A positive declared component vanished during mass-to-mole conversion.
    UnrepresentableMoleFraction {
        /// Canonical component index.
        component: usize,
        /// Affected species.
        species: SpeciesId,
    },
    /// Component standard-state conventions are not identical.
    ConventionMismatch {
        /// Canonical component index.
        component: usize,
        /// Affected species.
        species: SpeciesId,
        /// Convention field that disagrees with component zero.
        field: FrozenMixtureConventionFieldV1,
        /// Convention value required by canonical component zero.
        expected: FrozenMixtureConventionValueV1,
        /// Convention value declared by the affected component.
        found: FrozenMixtureConventionValueV1,
    },
    /// Mixture pressure must be finite and strictly positive.
    InvalidPressure {
        /// Exact rejected IEEE bits.
        bits: u64,
    },
    /// Canonical positive inputs produced a non-positive or non-finite mixture molar mass.
    InvalidMixtureMolarMass {
        /// Exact rejected derived-value bits.
        bits: u64,
    },
    /// One active species standard-state evaluation refused.
    ComponentEvaluation {
        /// Canonical component index.
        component: usize,
        /// Affected species.
        species: SpeciesId,
        /// Underlying typed refusal.
        source: Box<ThermochemErrorV1>,
    },
    /// Fixed-order mixture arithmetic produced a non-finite field.
    NonFiniteEvaluation {
        /// First field that failed.
        property: FrozenMixturePropertyV1,
        /// Exact non-finite IEEE bits.
        bits: u64,
    },
}

impl fmt::Display for FrozenIdealGasMixtureErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyComponents => f.write_str("frozen ideal-gas mixture requires a component"),
            Self::TooManyComponents { offered, limit } => write!(
                f,
                "frozen ideal-gas mixture offered {offered} components, exceeding limit {limit}"
            ),
            Self::CompositionLengthMismatch {
                components,
                fractions,
            } => write!(
                f,
                "frozen ideal-gas mixture has {components} models but {fractions} fractions"
            ),
            Self::Composition(error) => write!(f, "frozen mixture composition refused: {error}"),
            Self::ZeroFraction {
                component,
                species,
                bits,
            } => write!(
                f,
                "frozen mixture component {component} ({species}) has zero fraction bits {bits:#018x}; omit absent species"
            ),
            Self::FractionSumNotExact { basis, sum_bits } => write!(
                f,
                "frozen mixture {basis:?} fractions must sum to exact 1.0 in canonical order (sum bits {sum_bits:#018x})"
            ),
            Self::DuplicateSpecies { species } => {
                write!(f, "frozen mixture repeats active species {species}")
            }
            Self::UnrepresentableMoleFraction { component, species } => write!(
                f,
                "frozen mixture component {component} ({species}) lost its positive mole fraction during basis conversion"
            ),
            Self::ConventionMismatch {
                component,
                species,
                field,
                expected,
                found,
            } => write!(
                f,
                "frozen mixture component {component} ({species}) disagrees on {field:?}: expected {expected:?}, found {found:?}"
            ),
            Self::InvalidPressure { bits } => write!(
                f,
                "frozen mixture pressure must be positive and finite (bits {bits:#018x})"
            ),
            Self::InvalidMixtureMolarMass { bits } => write!(
                f,
                "frozen mixture molar mass must remain positive and finite (bits {bits:#018x})"
            ),
            Self::ComponentEvaluation {
                component,
                species,
                source,
            } => write!(
                f,
                "frozen mixture component {component} ({species}) refused: {source}"
            ),
            Self::NonFiniteEvaluation { property, bits } => write!(
                f,
                "frozen mixture {property:?} became non-finite (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for FrozenIdealGasMixtureErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Composition(error) => Some(error),
            Self::ComponentEvaluation { source, .. } => Some(source.as_ref()),
            _ => None,
        }
    }
}

impl From<SemanticError> for FrozenIdealGasMixtureErrorV1 {
    fn from(value: SemanticError) -> Self {
        Self::Composition(value)
    }
}

/// Canonical frozen-composition ideal-gas mixture model.
///
/// Component order is canonical `SpeciesId` order. Zero-fraction declarations
/// refuse and must be omitted. All active models must use exact-matching phase,
/// EOS, reference pressure, and elemental-reference conventions. Those
/// conventions remain caller assertions because the generic source-card schema
/// does not authenticate species or phase.
#[derive(Debug, Clone, PartialEq)]
pub struct FrozenIdealGasMixtureModelV1 {
    components: Vec<Nasa9StandardStateModelV1>,
    declared_composition: Composition,
    mole_composition: Composition,
}

impl FrozenIdealGasMixtureModelV1 {
    /// Canonicalize component order and convert the declared mass/mole basis.
    ///
    /// `composition.fractions()[i]` initially belongs to `components[i]`; each
    /// pair moves together when canonicalized. Zero fractions are semantically
    /// absent and must be omitted. Fractions must sum to bit-exact one in
    /// canonical species order, both before and after any mass-to-mole
    /// conversion. A volume-fraction composition refuses because no ideal-gas
    /// volume-basis conversion is implied by this L1 boundary.
    ///
    /// # Errors
    /// Refuses empty/excess/misaligned/duplicate components, a non-exact sum,
    /// unsupported or unrepresentable composition conversion, or convention
    /// drift.
    #[allow(clippy::too_many_lines)] // One ordered fail-closed canonicalization transaction.
    pub fn new(
        components: Vec<Nasa9StandardStateModelV1>,
        composition: Composition,
    ) -> Result<Self, FrozenIdealGasMixtureErrorV1> {
        if components.is_empty() {
            return Err(FrozenIdealGasMixtureErrorV1::EmptyComponents);
        }
        if components.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
            return Err(FrozenIdealGasMixtureErrorV1::TooManyComponents {
                offered: components.len(),
                limit: MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1,
            });
        }
        if components.len() != composition.fractions().len() {
            return Err(FrozenIdealGasMixtureErrorV1::CompositionLengthMismatch {
                components: components.len(),
                fractions: composition.fractions().len(),
            });
        }

        let basis = composition.basis();
        let mut paired: Vec<(Nasa9StandardStateModelV1, f64)> = components
            .into_iter()
            .zip(composition.fractions().iter().copied())
            .collect();
        paired.sort_by(|(left, _), (right, _)| left.species().cmp(right.species()));
        for pair in paired.windows(2) {
            if pair[0].0.species() == pair[1].0.species() {
                return Err(FrozenIdealGasMixtureErrorV1::DuplicateSpecies {
                    species: pair[0].0.species().clone(),
                });
            }
        }
        for (component, (model, fraction)) in paired.iter().enumerate() {
            if *fraction <= 0.0 {
                return Err(FrozenIdealGasMixtureErrorV1::ZeroFraction {
                    component,
                    species: model.species().clone(),
                    bits: fraction.to_bits(),
                });
            }
        }
        let declared_sum: f64 = paired.iter().map(|(_, fraction)| *fraction).sum();
        if declared_sum.to_bits() != 1.0_f64.to_bits() {
            return Err(FrozenIdealGasMixtureErrorV1::FractionSumNotExact {
                basis,
                sum_bits: declared_sum.to_bits(),
            });
        }

        let declared_composition = Composition::new(
            basis,
            paired.iter().map(|(_, fraction)| *fraction).collect(),
        )?;
        let molar_masses: Vec<MolarMass> =
            paired.iter().map(|(model, _)| model.molar_mass()).collect();
        let mole_composition = declared_composition.to_mole_fractions(&molar_masses)?;
        for (component, (&declared, &mole)) in declared_composition
            .fractions()
            .iter()
            .zip(mole_composition.fractions())
            .enumerate()
        {
            if declared > 0.0 && mole <= 0.0 {
                return Err(FrozenIdealGasMixtureErrorV1::UnrepresentableMoleFraction {
                    component,
                    species: paired[component].0.species().clone(),
                });
            }
        }
        let mole_sum: f64 = mole_composition.fractions().iter().sum();
        if mole_sum.to_bits() != 1.0_f64.to_bits() {
            return Err(FrozenIdealGasMixtureErrorV1::FractionSumNotExact {
                basis: CompositionBasis::MoleFraction,
                sum_bits: mole_sum.to_bits(),
            });
        }

        let reference = paired[0].0.convention();
        for (component, (model, _)) in paired.iter().enumerate().skip(1) {
            let convention = model.convention();
            let mismatch = if convention.phase() != reference.phase() {
                Some((
                    FrozenMixtureConventionFieldV1::Phase,
                    FrozenMixtureConventionValueV1::Phase(reference.phase()),
                    FrozenMixtureConventionValueV1::Phase(convention.phase()),
                ))
            } else if convention.eos() != reference.eos() {
                Some((
                    FrozenMixtureConventionFieldV1::EquationOfState,
                    FrozenMixtureConventionValueV1::EquationOfState(reference.eos()),
                    FrozenMixtureConventionValueV1::EquationOfState(convention.eos()),
                ))
            } else if convention.reference_pressure().value().to_bits()
                != reference.reference_pressure().value().to_bits()
            {
                Some((
                    FrozenMixtureConventionFieldV1::ReferencePressure,
                    FrozenMixtureConventionValueV1::ReferencePressureBits(
                        reference.reference_pressure().value().to_bits(),
                    ),
                    FrozenMixtureConventionValueV1::ReferencePressureBits(
                        convention.reference_pressure().value().to_bits(),
                    ),
                ))
            } else if convention.elemental_reference() != reference.elemental_reference() {
                Some((
                    FrozenMixtureConventionFieldV1::ElementalReference,
                    FrozenMixtureConventionValueV1::ElementalReference(
                        reference.elemental_reference().clone(),
                    ),
                    FrozenMixtureConventionValueV1::ElementalReference(
                        convention.elemental_reference().clone(),
                    ),
                ))
            } else {
                None
            };
            if let Some((field, expected, found)) = mismatch {
                return Err(FrozenIdealGasMixtureErrorV1::ConventionMismatch {
                    component,
                    species: model.species().clone(),
                    field,
                    expected,
                    found,
                });
            }
        }

        Ok(Self {
            components: paired.into_iter().map(|(model, _)| model).collect(),
            declared_composition,
            mole_composition,
        })
    }

    /// Active component models in canonical species order.
    #[must_use]
    pub fn components(&self) -> &[Nasa9StandardStateModelV1] {
        &self.components
    }

    /// Canonically ordered composition in the caller-declared basis.
    #[must_use]
    pub const fn declared_composition(&self) -> &Composition {
        &self.declared_composition
    }

    /// Canonically ordered mole fractions used by the evaluator.
    #[must_use]
    pub const fn mole_composition(&self) -> &Composition {
        &self.mole_composition
    }

    /// Evaluate frozen ideal-gas mixture properties at absolute `T` and `p`.
    ///
    /// Molar `cp`, `h`, and standard entropy are mole-fraction sums of species
    /// values. Entropy then applies ideal mixing and common-reference-pressure
    /// corrections. Frozen `cv = cp - R`, `u = h - R T`, and `g = h - T s`.
    /// Mass-specific values divide by the mole-fraction mixture molar mass.
    ///
    /// # Errors
    /// Refuses invalid pressure, any active species evaluation failure, or the
    /// first non-finite mixture field. Nothing partial escapes.
    #[allow(clippy::too_many_lines)] // One fixed operation tree and its exact receipt stay adjacent.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<FrozenIdealGasMixtureEvaluationV1, FrozenIdealGasMixtureErrorV1> {
        let pressure_value = pressure.value();
        if !pressure_value.is_finite() || pressure_value <= 0.0 {
            return Err(FrozenIdealGasMixtureErrorV1::InvalidPressure {
                bits: pressure_value.to_bits(),
            });
        }
        let t = temperature.value();
        let reference_pressure = self.components[0].convention().reference_pressure().value();
        let mole_fractions = self.mole_composition.fractions();
        let mut cp = 0.0;
        let mut h = 0.0;
        let mut standard_entropy = 0.0;
        let mut molar_mass = 0.0;
        let mut x_log_x = 0.0;
        let mut component_receipts = Vec::with_capacity(self.components.len());

        for (component, (model, &mole_fraction)) in
            self.components.iter().zip(mole_fractions).enumerate()
        {
            let evaluation = model.evaluate(temperature).map_err(|source| {
                FrozenIdealGasMixtureErrorV1::ComponentEvaluation {
                    component,
                    species: model.species().clone(),
                    source: Box::new(source),
                }
            })?;
            let properties = evaluation.properties();
            cp += mole_fraction * properties.cp().value();
            h += mole_fraction * properties.h().value();
            standard_entropy += mole_fraction * properties.s().value();
            molar_mass += mole_fraction * model.molar_mass().value();
            x_log_x += mole_fraction * fs_math::det::ln(mole_fraction);
            component_receipts.push(evaluation.receipt().clone());
        }

        let pressure_log = fs_math::det::ln(pressure_value) - fs_math::det::ln(reference_pressure);
        let entropy = standard_entropy
            - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * x_log_x
            - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * pressure_log;
        let cv = cp - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t;
        let internal_energy = h - rt;
        let gibbs_energy = h - t * entropy;
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err(FrozenIdealGasMixtureErrorV1::InvalidMixtureMolarMass {
                bits: molar_mass.to_bits(),
            });
        }
        let mass_specific_cp = cp / molar_mass;
        let mass_specific_gas_constant = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass;
        let mass_specific_cv = cv / molar_mass;

        for (property, value) in [
            (FrozenMixturePropertyV1::MolarCp, cp),
            (FrozenMixturePropertyV1::MolarCv, cv),
            (FrozenMixturePropertyV1::MolarEnthalpy, h),
            (
                FrozenMixturePropertyV1::MolarInternalEnergy,
                internal_energy,
            ),
            (FrozenMixturePropertyV1::MolarEntropy, entropy),
            (FrozenMixturePropertyV1::MolarGibbsEnergy, gibbs_energy),
            (FrozenMixturePropertyV1::MolarMass, molar_mass),
            (FrozenMixturePropertyV1::MassSpecificCp, mass_specific_cp),
            (FrozenMixturePropertyV1::MassSpecificCv, mass_specific_cv),
            (
                FrozenMixturePropertyV1::MassSpecificGasConstant,
                mass_specific_gas_constant,
            ),
        ] {
            if !value.is_finite() {
                return Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
                    property,
                    bits: value.to_bits(),
                });
            }
        }

        let properties = FrozenIdealGasMixturePropertiesV1 {
            cp: MolarHeatCapacityV1::new(cp),
            cv: MolarIsochoricHeatCapacityV1::new(cv),
            h: MolarEnthalpyV1::new(h),
            u: MolarInternalEnergyV1::new(internal_energy),
            s: MolarEntropyV1::new(entropy),
            g: MolarGibbsEnergyV1::new(gibbs_energy),
            molar_mass: MolarMass::new(molar_mass),
            mass_specific_cp: MassSpecificIsobaricHeatCapacityV1::new(mass_specific_cp),
            mass_specific_cv: MassSpecificIsochoricHeatCapacityV1::new(mass_specific_cv),
            mass_specific_gas_constant: MassSpecificGasConstantV1::new(mass_specific_gas_constant),
        };
        let receipt = FrozenIdealGasMixtureReceiptV1 {
            evaluator_version: FROZEN_IDEAL_GAS_MIXTURE_EVALUATOR_VERSION_V1,
            fs_math_version: fs_math::VERSION,
            fs_qty_version: fs_qty::VERSION,
            gas_constant_bits: UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits(),
            temperature_bits: t.to_bits(),
            pressure_bits: pressure_value.to_bits(),
            reference_pressure_bits: reference_pressure.to_bits(),
            phase: self.components[0].convention().phase(),
            eos: self.components[0].convention().eos(),
            elemental_reference: self.components[0]
                .convention()
                .elemental_reference()
                .clone(),
            declared_basis: self.declared_composition.basis(),
            declared_fraction_sum_bits: self
                .declared_composition
                .fractions()
                .iter()
                .sum::<f64>()
                .to_bits(),
            mole_fraction_sum_bits: mole_fractions.iter().sum::<f64>().to_bits(),
            x_log_x_bits: x_log_x.to_bits(),
            molar_mass_bits: molar_mass.to_bits(),
            species: self
                .components
                .iter()
                .map(|model| model.species().clone())
                .collect(),
            declared_fraction_bits: self
                .declared_composition
                .fractions()
                .iter()
                .map(|fraction| fraction.to_bits())
                .collect(),
            mole_fraction_bits: mole_fractions
                .iter()
                .map(|fraction| fraction.to_bits())
                .collect(),
            component_receipts,
        };
        Ok(FrozenIdealGasMixtureEvaluationV1 {
            temperature,
            pressure,
            properties,
            receipt,
        })
    }
}

/// Frozen ideal-gas mixture properties in molar and mass-specific bases.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrozenIdealGasMixturePropertiesV1 {
    cp: MolarHeatCapacityV1,
    cv: MolarIsochoricHeatCapacityV1,
    h: MolarEnthalpyV1,
    u: MolarInternalEnergyV1,
    s: MolarEntropyV1,
    g: MolarGibbsEnergyV1,
    molar_mass: MolarMass,
    mass_specific_cp: MassSpecificIsobaricHeatCapacityV1,
    mass_specific_cv: MassSpecificIsochoricHeatCapacityV1,
    mass_specific_gas_constant: MassSpecificGasConstantV1,
}

impl FrozenIdealGasMixturePropertiesV1 {
    /// Frozen molar heat capacity at constant pressure.
    #[must_use]
    pub const fn cp(self) -> MolarHeatCapacityV1 {
        self.cp
    }

    /// Frozen molar heat capacity at constant volume.
    #[must_use]
    pub const fn cv(self) -> MolarIsochoricHeatCapacityV1 {
        self.cv
    }

    /// Mixture molar enthalpy.
    #[must_use]
    pub const fn h(self) -> MolarEnthalpyV1 {
        self.h
    }

    /// Mixture molar internal energy.
    #[must_use]
    pub const fn u(self) -> MolarInternalEnergyV1 {
        self.u
    }

    /// Mixture molar entropy, including ideal mixing and pressure correction.
    #[must_use]
    pub const fn s(self) -> MolarEntropyV1 {
        self.s
    }

    /// Mixture molar Gibbs energy.
    #[must_use]
    pub const fn g(self) -> MolarGibbsEnergyV1 {
        self.g
    }

    /// Mole-fraction-weighted mixture molar mass.
    #[must_use]
    pub const fn molar_mass(self) -> MolarMass {
        self.molar_mass
    }

    /// Frozen mass-specific heat capacity at constant pressure.
    #[must_use]
    pub const fn mass_specific_cp(self) -> MassSpecificIsobaricHeatCapacityV1 {
        self.mass_specific_cp
    }

    /// Frozen mass-specific heat capacity at constant volume.
    #[must_use]
    pub const fn mass_specific_cv(self) -> MassSpecificIsochoricHeatCapacityV1 {
        self.mass_specific_cv
    }

    /// Mixture-specific gas constant `R / M_mix`.
    #[must_use]
    pub const fn mass_specific_gas_constant(self) -> MassSpecificGasConstantV1 {
        self.mass_specific_gas_constant
    }
}

/// Immutable exact-field frozen-mixture evaluation receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrozenIdealGasMixtureReceiptV1 {
    evaluator_version: u32,
    fs_math_version: &'static str,
    fs_qty_version: &'static str,
    gas_constant_bits: u64,
    temperature_bits: u64,
    pressure_bits: u64,
    reference_pressure_bits: u64,
    phase: StandardStatePhaseV1,
    eos: ReferenceEquationOfStateV1,
    elemental_reference: ElementalReferenceIdV1,
    declared_basis: CompositionBasis,
    declared_fraction_sum_bits: u64,
    mole_fraction_sum_bits: u64,
    x_log_x_bits: u64,
    molar_mass_bits: u64,
    species: Vec<SpeciesId>,
    declared_fraction_bits: Vec<u64>,
    mole_fraction_bits: Vec<u64>,
    component_receipts: Vec<Nasa9EvaluationReceiptV1>,
}

impl FrozenIdealGasMixtureReceiptV1 {
    /// Mixture evaluator/operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Deterministic elementary-math crate version.
    #[must_use]
    pub const fn fs_math_version(&self) -> &'static str {
        self.fs_math_version
    }

    /// Quantity-semantics crate version used for basis conversion.
    #[must_use]
    pub const fn fs_qty_version(&self) -> &'static str {
        self.fs_qty_version
    }

    /// Exact universal-gas-constant bits.
    #[must_use]
    pub const fn gas_constant_bits(&self) -> u64 {
        self.gas_constant_bits
    }

    /// Exact evaluation-temperature bits.
    #[must_use]
    pub const fn temperature_bits(&self) -> u64 {
        self.temperature_bits
    }

    /// Exact mixture-pressure bits.
    #[must_use]
    pub const fn pressure_bits(&self) -> u64 {
        self.pressure_bits
    }

    /// Exact common standard-state reference-pressure bits.
    #[must_use]
    pub const fn reference_pressure_bits(&self) -> u64 {
        self.reference_pressure_bits
    }

    /// Common caller-declared standard-state phase.
    #[must_use]
    pub const fn phase(&self) -> StandardStatePhaseV1 {
        self.phase
    }

    /// Common caller-declared reference EOS.
    #[must_use]
    pub const fn eos(&self) -> ReferenceEquationOfStateV1 {
        self.eos
    }

    /// Common caller-declared elemental-reference convention.
    #[must_use]
    pub const fn elemental_reference(&self) -> &ElementalReferenceIdV1 {
        &self.elemental_reference
    }

    /// Caller-declared composition basis.
    #[must_use]
    pub const fn declared_basis(&self) -> CompositionBasis {
        self.declared_basis
    }

    /// Exact canonical sum bits of the declared fractions.
    #[must_use]
    pub const fn declared_fraction_sum_bits(&self) -> u64 {
        self.declared_fraction_sum_bits
    }

    /// Exact canonical sum bits of the evaluated mole fractions.
    #[must_use]
    pub const fn mole_fraction_sum_bits(&self) -> u64 {
        self.mole_fraction_sum_bits
    }

    /// Exact `sum(x_i * ln(x_i))` accumulator bits.
    #[must_use]
    pub const fn x_log_x_bits(&self) -> u64 {
        self.x_log_x_bits
    }

    /// Exact evaluated mixture-molar-mass bits.
    #[must_use]
    pub const fn molar_mass_bits(&self) -> u64 {
        self.molar_mass_bits
    }

    /// Canonically ordered active species.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Exact declared-fraction bits in canonical species order.
    #[must_use]
    pub fn declared_fraction_bits(&self) -> &[u64] {
        &self.declared_fraction_bits
    }

    /// Exact mole-fraction bits used by the evaluator.
    #[must_use]
    pub fn mole_fraction_bits(&self) -> &[u64] {
        &self.mole_fraction_bits
    }

    /// Full active species evaluation receipts in canonical species order.
    #[must_use]
    pub fn component_receipts(&self) -> &[Nasa9EvaluationReceiptV1] {
        &self.component_receipts
    }
}

/// One successful frozen ideal-gas mixture evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct FrozenIdealGasMixtureEvaluationV1 {
    temperature: Temperature,
    pressure: Pressure,
    properties: FrozenIdealGasMixturePropertiesV1,
    receipt: FrozenIdealGasMixtureReceiptV1,
}

impl FrozenIdealGasMixtureEvaluationV1 {
    /// Evaluated absolute temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Evaluated mixture pressure.
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }

    /// Derived frozen ideal-gas properties.
    #[must_use]
    pub const fn properties(&self) -> FrozenIdealGasMixturePropertiesV1 {
        self.properties
    }

    /// Exact composition, convention, and component-provenance receipt.
    #[must_use]
    pub const fn receipt(&self) -> &FrozenIdealGasMixtureReceiptV1 {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
        NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, ReferenceEquationOfStateV1,
        StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
    use std::collections::BTreeMap;

    const REFERENCE_PRESSURE: f64 = 100_000.0;

    fn standard_state_model(
        species: &str,
        molar_mass: f64,
        cp_over_r: f64,
        reference_pressure: f64,
        elemental_reference: &str,
    ) -> Nasa9StandardStateModelV1 {
        let mut coefficients = [0.0; 9];
        coefficients[2] = cp_over_r;
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
        let card = ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.to_string()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 6_000.0),
            sources: Vec::new(),
            provenance: Provenance {
                source: format!("synthetic ideal-gas fixture for {species}"),
                license: "test fixture".to_string(),
                artifact: None,
            },
        };
        Nasa9StandardStateModelV1::new(
            SpeciesId::new(species).expect("canonical species"),
            MolarMass::new(molar_mass),
            StandardStateConventionV1::new(
                StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas,
                Pressure::new(reference_pressure),
                ElementalReferenceIdV1::new(elemental_reference)
                    .expect("canonical elemental reference"),
            )
            .expect("valid standard-state convention"),
            vec![Nasa9RegionV1::from_card(card).expect("valid synthetic NASA-9 card")],
        )
        .expect("valid synthetic species model")
    }

    fn composition(basis: CompositionBasis, fractions: Vec<f64>) -> Composition {
        Composition::new(basis, fractions).expect("valid test composition")
    }

    fn assert_close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        assert!(
            (actual - expected).abs() <= 128.0 * f64::EPSILON * scale,
            "actual {actual:?} differs from expected {expected:?}",
        );
    }

    #[test]
    fn g0_single_species_limit_recovers_standard_state_and_frozen_identities() {
        let species = standard_state_model(
            "N2",
            0.028_013_4,
            3.5,
            REFERENCE_PRESSURE,
            "reference-elements",
        );
        let standard = species
            .evaluate(Temperature::new(1_000.0))
            .expect("species standard state");
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![species],
            composition(CompositionBasis::MoleFraction, vec![1.0]),
        )
        .expect("single-species mixture");
        let evaluation = mixture
            .evaluate(Temperature::new(1_000.0), Pressure::new(REFERENCE_PRESSURE))
            .expect("single-species evaluation");
        let properties = evaluation.properties();
        let standard_properties = standard.properties();

        assert_eq!(
            properties.cp().value().to_bits(),
            standard_properties.cp().value().to_bits(),
        );
        assert_eq!(
            properties.h().value().to_bits(),
            standard_properties.h().value().to_bits(),
        );
        assert_eq!(
            properties.s().value().to_bits(),
            standard_properties.s().value().to_bits(),
        );
        assert_eq!(
            properties.u().value().to_bits(),
            standard_properties.u().value().to_bits(),
        );
        assert_eq!(
            properties.g().value().to_bits(),
            standard_properties.g().value().to_bits(),
        );
        assert_close(
            properties.cp().value() - properties.cv().value(),
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
        );
        assert_close(
            properties.mass_specific_cp().value() - properties.mass_specific_cv().value(),
            properties.mass_specific_gas_constant().value(),
        );
        assert_eq!(evaluation.receipt().species()[0].as_str(), "N2");
        assert_eq!(evaluation.receipt().component_receipts().len(), 1);
        assert_eq!(
            evaluation.receipt().declared_fraction_sum_bits(),
            1.0_f64.to_bits(),
        );
        assert_eq!(
            evaluation.receipt().mole_fraction_sum_bits(),
            1.0_f64.to_bits(),
        );
    }

    #[test]
    fn g0_binary_mixture_matches_weighted_and_ideal_mixing_closed_forms() {
        let component_a =
            standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements");
        let component_b =
            standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements");
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![component_a.clone(), component_b.clone()],
            composition(CompositionBasis::MoleFraction, vec![0.25, 0.75]),
        )
        .expect("binary mixture");
        let temperature = 1_000.0;
        let pressure = 2.0 * REFERENCE_PRESSURE;
        let evaluation = mixture
            .evaluate(Temperature::new(temperature), Pressure::new(pressure))
            .expect("binary evaluation");
        let properties = evaluation.properties();
        let expected_cp = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 4.5;
        let expected_h = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature * 4.5;
        let x_log_x = 0.25 * fs_math::det::ln(0.25) + 0.75 * fs_math::det::ln(0.75);
        let expected_entropy = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * (4.5 * fs_math::det::ln(temperature))
            - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * x_log_x
            - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * (fs_math::det::ln(pressure) - fs_math::det::ln(REFERENCE_PRESSURE));
        let expected_molar_mass = 0.25 * 0.002 + 0.75 * 0.004;

        assert_close(properties.cp().value(), expected_cp);
        assert_close(
            properties.cv().value(),
            expected_cp - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
        );
        assert_close(properties.h().value(), expected_h);
        assert_close(
            properties.u().value(),
            expected_h - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature,
        );
        assert_close(properties.s().value(), expected_entropy);
        assert_close(
            properties.g().value(),
            expected_h - temperature * expected_entropy,
        );
        assert_close(properties.molar_mass().value(), expected_molar_mass);
        assert_close(
            properties.mass_specific_cp().value(),
            expected_cp / expected_molar_mass,
        );
        assert_close(
            properties.mass_specific_cv().value(),
            properties.cv().value() / expected_molar_mass,
        );
        assert_close(
            properties.mass_specific_gas_constant().value(),
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / expected_molar_mass,
        );
        assert_close(
            properties.mass_specific_cp().value() - properties.mass_specific_cv().value(),
            properties.mass_specific_gas_constant().value(),
        );

        let standard_a = component_a
            .evaluate(Temperature::new(temperature))
            .expect("component A standard state")
            .properties();
        let standard_b = component_b
            .evaluate(Temperature::new(temperature))
            .expect("component B standard state")
            .properties();
        let alternative_gibbs = 0.25 * standard_a.g().value()
            + 0.75 * standard_b.g().value()
            + UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * temperature
                * (x_log_x + fs_math::det::ln(pressure) - fs_math::det::ln(REFERENCE_PRESSURE));
        assert_close(properties.g().value(), alternative_gibbs);
        assert_eq!(evaluation.receipt().x_log_x_bits(), x_log_x.to_bits(),);
        assert_eq!(
            evaluation.receipt().molar_mass_bits(),
            properties.molar_mass().value().to_bits(),
        );
    }

    #[test]
    fn g3_component_permutation_and_g5_replay_are_bit_identical() {
        let component_a =
            standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements");
        let component_b =
            standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements");
        let canonical = FrozenIdealGasMixtureModelV1::new(
            vec![component_a.clone(), component_b.clone()],
            composition(CompositionBasis::MoleFraction, vec![0.25, 0.75]),
        )
        .expect("canonical mixture");
        let permuted = FrozenIdealGasMixtureModelV1::new(
            vec![component_b, component_a],
            composition(CompositionBasis::MoleFraction, vec![0.75, 0.25]),
        )
        .expect("permuted mixture");
        assert_eq!(canonical, permuted);

        let temperature = Temperature::new(900.0);
        let pressure = Pressure::new(150_000.0);
        let first = canonical
            .evaluate(temperature, pressure)
            .expect("canonical evaluation");
        let replay = canonical
            .evaluate(temperature, pressure)
            .expect("replay evaluation");
        let permuted_evaluation = permuted
            .evaluate(temperature, pressure)
            .expect("permuted evaluation");
        assert_eq!(first, replay);
        assert_eq!(first, permuted_evaluation);
        assert_eq!(first.receipt().species()[0].as_str(), "A");
        assert_eq!(first.receipt().species()[1].as_str(), "B");
    }

    #[test]
    fn g5_receipt_binds_exact_mixture_and_nested_component_authority() {
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![
                standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements"),
                standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements"),
            ],
            composition(CompositionBasis::MoleFraction, vec![0.25, 0.75]),
        )
        .expect("binary mixture");
        let temperature = 850.0;
        let pressure = 125_000.0;
        let evaluation = mixture
            .evaluate(Temperature::new(temperature), Pressure::new(pressure))
            .expect("binary evaluation");
        let receipt = evaluation.receipt();

        assert_eq!(
            receipt.evaluator_version(),
            FROZEN_IDEAL_GAS_MIXTURE_EVALUATOR_VERSION_V1,
        );
        assert_eq!(receipt.fs_math_version(), fs_math::VERSION);
        assert_eq!(receipt.fs_qty_version(), fs_qty::VERSION);
        assert_eq!(
            receipt.gas_constant_bits(),
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits(),
        );
        assert_eq!(receipt.temperature_bits(), temperature.to_bits());
        assert_eq!(receipt.pressure_bits(), pressure.to_bits());
        assert_eq!(
            receipt.reference_pressure_bits(),
            REFERENCE_PRESSURE.to_bits(),
        );
        assert_eq!(receipt.phase(), StandardStatePhaseV1::Gas);
        assert_eq!(receipt.eos(), ReferenceEquationOfStateV1::IdealGas);
        assert_eq!(receipt.elemental_reference().as_str(), "reference-elements");
        assert_eq!(
            receipt.declared_fraction_bits(),
            &[0.25_f64.to_bits(), 0.75_f64.to_bits()],
        );
        assert_eq!(
            receipt.mole_fraction_bits(),
            receipt.declared_fraction_bits(),
        );
        for (species, component) in receipt.species().iter().zip(receipt.component_receipts()) {
            assert_eq!(species, component.species());
            assert_eq!(component.temperature_bits(), receipt.temperature_bits());
            assert_eq!(
                component.reference_pressure_bits(),
                receipt.reference_pressure_bits(),
            );
        }
    }

    #[test]
    fn g3_receipt_fields_change_with_pressure_temperature_and_composition() {
        let components = || {
            vec![
                standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements"),
                standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements"),
            ]
        };
        let mixture = FrozenIdealGasMixtureModelV1::new(
            components(),
            composition(CompositionBasis::MoleFraction, vec![0.25, 0.75]),
        )
        .expect("binary mixture");
        let base = mixture
            .evaluate(Temperature::new(800.0), Pressure::new(100_000.0))
            .expect("base evaluation");
        let pressure_changed = mixture
            .evaluate(Temperature::new(800.0), Pressure::new(200_000.0))
            .expect("pressure mutation");
        let temperature_changed = mixture
            .evaluate(Temperature::new(900.0), Pressure::new(100_000.0))
            .expect("temperature mutation");
        let composition_changed = FrozenIdealGasMixtureModelV1::new(
            components(),
            composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
        )
        .expect("composition mutation")
        .evaluate(Temperature::new(800.0), Pressure::new(100_000.0))
        .expect("composition-mutated evaluation");

        assert_ne!(
            base.receipt().pressure_bits(),
            pressure_changed.receipt().pressure_bits(),
        );
        assert_eq!(
            base.receipt().component_receipts(),
            pressure_changed.receipt().component_receipts(),
        );
        assert_ne!(
            base.receipt().temperature_bits(),
            temperature_changed.receipt().temperature_bits(),
        );
        assert_ne!(
            base.receipt().component_receipts(),
            temperature_changed.receipt().component_receipts(),
        );
        assert_ne!(
            base.receipt().mole_fraction_bits(),
            composition_changed.receipt().mole_fraction_bits(),
        );
        assert_ne!(
            base.receipt().x_log_x_bits(),
            composition_changed.receipt().x_log_x_bits(),
        );
    }

    #[test]
    fn g0_mass_basis_conversion_matches_equivalent_mole_basis() {
        let component_a =
            standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements");
        let component_b =
            standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements");
        let mole = FrozenIdealGasMixtureModelV1::new(
            vec![component_a.clone(), component_b.clone()],
            composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
        )
        .expect("mole-basis mixture");
        let mass = FrozenIdealGasMixtureModelV1::new(
            vec![component_a, component_b],
            composition(CompositionBasis::MassFraction, vec![1.0 / 3.0, 2.0 / 3.0]),
        )
        .expect("mass-basis mixture");
        assert_eq!(
            mass.mole_composition().fractions(),
            mole.mole_composition().fractions(),
        );

        let mole_evaluation = mole
            .evaluate(Temperature::new(700.0), Pressure::new(REFERENCE_PRESSURE))
            .expect("mole-basis evaluation");
        let mass_evaluation = mass
            .evaluate(Temperature::new(700.0), Pressure::new(REFERENCE_PRESSURE))
            .expect("mass-basis evaluation");
        assert_eq!(mole_evaluation.properties(), mass_evaluation.properties());
        assert_eq!(
            mass_evaluation.receipt().declared_basis(),
            CompositionBasis::MassFraction,
        );
        assert_eq!(
            mole_evaluation.receipt().declared_basis(),
            CompositionBasis::MoleFraction,
        );
    }

    #[test]
    fn g3_pressure_scaling_changes_only_entropy_and_gibbs() {
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![
                standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements"),
                standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements"),
            ],
            composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
        )
        .expect("binary mixture");
        let temperature = Temperature::new(800.0);
        let at_reference = mixture
            .evaluate(temperature, Pressure::new(REFERENCE_PRESSURE))
            .expect("reference-pressure evaluation");
        let doubled = mixture
            .evaluate(temperature, Pressure::new(2.0 * REFERENCE_PRESSURE))
            .expect("doubled-pressure evaluation");
        let reference = at_reference.properties();
        let shifted = doubled.properties();

        assert_eq!(
            reference.cp().value().to_bits(),
            shifted.cp().value().to_bits()
        );
        assert_eq!(
            reference.cv().value().to_bits(),
            shifted.cv().value().to_bits()
        );
        assert_eq!(
            reference.h().value().to_bits(),
            shifted.h().value().to_bits()
        );
        assert_eq!(
            reference.u().value().to_bits(),
            shifted.u().value().to_bits()
        );
        assert_eq!(
            reference.molar_mass().value().to_bits(),
            shifted.molar_mass().value().to_bits(),
        );
        let entropy_shift = -UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * fs_math::det::ln(2.0);
        assert_close(shifted.s().value() - reference.s().value(), entropy_shift);
        assert_close(
            shifted.g().value() - reference.g().value(),
            -temperature.value() * entropy_shift,
        );
    }

    #[test]
    fn g3_construction_refuses_noncanonical_compositions() {
        let model_a =
            || standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements");
        let model_b =
            || standard_state_model("B", 0.004, 5.0, REFERENCE_PRESSURE, "reference-elements");

        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                Vec::new(),
                composition(CompositionBasis::MoleFraction, vec![1.0]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::EmptyComponents)
        ));
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a()],
                composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::CompositionLengthMismatch { .. })
        ));
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a(), model_b()],
                composition(CompositionBasis::MoleFraction, vec![1.0, 0.0]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::ZeroFraction { .. })
        ));
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a(), model_b()],
                composition(
                    CompositionBasis::MoleFraction,
                    vec![0.5, 0.500_000_000_000_5],
                ),
            ),
            Err(FrozenIdealGasMixtureErrorV1::FractionSumNotExact { .. })
        ));
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a(), model_b()],
                composition(CompositionBasis::MoleFraction, vec![1.0, -0.0]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::ZeroFraction { .. })
        ));
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a(), model_a()],
                composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::DuplicateSpecies { .. })
        ));
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a()],
                composition(CompositionBasis::VolumeFraction, vec![1.0]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::Composition(_))
        ));

        let over_limit = vec![model_a(); MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 + 1];
        let mut over_limit_fractions = vec![0.0; MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1];
        over_limit_fractions.push(1.0);
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                over_limit,
                composition(CompositionBasis::MoleFraction, over_limit_fractions),
            ),
            Err(FrozenIdealGasMixtureErrorV1::TooManyComponents { .. })
        ));
    }

    #[test]
    fn g3_construction_refuses_incompatible_standard_state_conventions() {
        let model_a =
            || standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements");
        let pressure_mismatch =
            standard_state_model("B", 0.004, 5.0, 101_325.0, "reference-elements");
        let pressure_result = FrozenIdealGasMixtureModelV1::new(
            vec![model_a(), pressure_mismatch],
            composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
        );
        assert!(matches!(
            pressure_result,
            Err(FrozenIdealGasMixtureErrorV1::ConventionMismatch {
                field: FrozenMixtureConventionFieldV1::ReferencePressure,
                species,
                expected: FrozenMixtureConventionValueV1::ReferencePressureBits(expected),
                found: FrozenMixtureConventionValueV1::ReferencePressureBits(found),
                ..
            }) if species.as_str() == "B"
                && expected == REFERENCE_PRESSURE.to_bits()
                && found == 101_325.0_f64.to_bits()
        ));
        let reference_mismatch = standard_state_model(
            "B",
            0.004,
            5.0,
            REFERENCE_PRESSURE,
            "other-reference-elements",
        );
        assert!(matches!(
            FrozenIdealGasMixtureModelV1::new(
                vec![model_a(), reference_mismatch],
                composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
            ),
            Err(FrozenIdealGasMixtureErrorV1::ConventionMismatch {
                field: FrozenMixtureConventionFieldV1::ElementalReference,
                ..
            })
        ));
    }

    #[test]
    fn g3_evaluation_refuses_pressure_component_and_basis_underflow_failures() {
        let model_a =
            standard_state_model("A", 0.002, 3.0, REFERENCE_PRESSURE, "reference-elements");
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![model_a.clone()],
            composition(CompositionBasis::MoleFraction, vec![1.0]),
        )
        .expect("single-component mixture");
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                mixture.evaluate(Temperature::new(300.0), Pressure::new(invalid)),
                Err(FrozenIdealGasMixtureErrorV1::InvalidPressure { .. })
            ));
        }
        assert!(matches!(
            mixture.evaluate(Temperature::new(199.0), Pressure::new(REFERENCE_PRESSURE),),
            Err(FrozenIdealGasMixtureErrorV1::ComponentEvaluation { component: 0, .. })
        ));

        let extreme_mass = FrozenIdealGasMixtureModelV1::new(
            vec![
                standard_state_model(
                    "A",
                    f64::MIN_POSITIVE,
                    3.0,
                    REFERENCE_PRESSURE,
                    "reference-elements",
                ),
                standard_state_model("B", f64::MAX, 5.0, REFERENCE_PRESSURE, "reference-elements"),
            ],
            composition(CompositionBasis::MassFraction, vec![0.5, 0.5]),
        );
        assert!(matches!(
            extreme_mass,
            Err(FrozenIdealGasMixtureErrorV1::UnrepresentableMoleFraction { .. })
        ));

        let smallest_subnormal = f64::from_bits(1);
        let underflowing_molar_mass = FrozenIdealGasMixtureModelV1::new(
            vec![
                standard_state_model(
                    "A",
                    smallest_subnormal,
                    3.0,
                    REFERENCE_PRESSURE,
                    "reference-elements",
                ),
                standard_state_model(
                    "B",
                    smallest_subnormal,
                    5.0,
                    REFERENCE_PRESSURE,
                    "reference-elements",
                ),
            ],
            composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
        )
        .expect("positive component molar masses remain admissible")
        .evaluate(Temperature::new(300.0), Pressure::new(REFERENCE_PRESSURE));
        assert!(matches!(
            underflowing_molar_mass,
            Err(FrozenIdealGasMixtureErrorV1::InvalidMixtureMolarMass { bits })
                if bits == 0.0_f64.to_bits()
        ));
    }
}
