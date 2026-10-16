//! Bounded dilute-gas transport-property mixing closures.
//!
//! Version 1 implements only Wilke's gas-mixture viscosity rule. It combines
//! caller-supplied pure-species viscosities at one explicitly declared
//! temperature; it does not evaluate those viscosities, authenticate their
//! common state, establish dilute-gas validity, or own a transport solver.

use core::fmt;

use fs_qty::semantic::{Composition, CompositionBasis, SemanticError};
use fs_qty::{DynViscosity, MolarMass, Temperature};

use crate::SpeciesId;

/// Version of the fixed Wilke operation tree and receipt layout.
pub const WILKE_GAS_VISCOSITY_EVALUATOR_VERSION_V1: u32 = 1;
/// Hard bound checked before component scanning or quadratic mixing work.
pub const MAX_WILKE_GAS_COMPONENTS_V1: usize = 128;

/// Transport mixing rule named by an exact-field receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GasViscosityMixingRuleV1 {
    /// Wilke's low-density gas-mixture viscosity rule.
    Wilke,
}

/// Component field named by a non-positive or non-finite refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WilkeComponentFieldV1 {
    /// Caller-declared molar mass.
    MolarMass,
    /// Caller-supplied pure-species dynamic viscosity.
    DynamicViscosity,
}

/// Intermediate named by a fail-closed arithmetic refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WilkeArithmeticValueV1 {
    /// Pure-viscosity ratio `mu_i / mu_j`.
    ViscosityRatio,
    /// Molar-mass ratio `M_j / M_i`.
    ReverseMolarMassRatio,
    /// Nested square root `(M_j / M_i)^(1/4)`.
    MolarMassQuarterRoot,
    /// Product `sqrt(mu_i / mu_j) (M_j / M_i)^(1/4)`.
    InteractionRootProduct,
    /// Numerator base `1 + sqrt(mu_i / mu_j) (M_j / M_i)^(1/4)`.
    InteractionNumeratorBase,
    /// Squared interaction numerator.
    InteractionNumerator,
    /// Molar-mass ratio `M_i / M_j` in the denominator.
    ForwardMolarMassRatio,
    /// Denominator radicand `8 (1 + M_i / M_j)`.
    InteractionDenominatorRadicand,
    /// Interaction denominator.
    InteractionDenominator,
    /// Wilke interaction factor `Phi_ij`.
    InteractionFactor,
    /// Weighted interaction `x_j Phi_ij`.
    WeightedInteraction,
    /// Fixed-order sum `sum_j x_j Phi_ij`.
    ComponentDenominator,
    /// Scaled pure viscosity `mu_i / denominator_i`.
    ScaledComponentViscosity,
    /// Component contribution `x_i mu_i / denominator_i`.
    ComponentContribution,
    /// Fixed-order mixture-viscosity sum.
    MixtureViscosity,
}

/// Allocation stage named by a bounded fallible reservation refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WilkeAllocationStageV1 {
    /// Canonical component/fraction pairing.
    CanonicalComponents,
    /// Canonical declared-fraction storage.
    DeclaredFractions,
    /// Canonical molar-mass storage used for basis conversion.
    MolarMasses,
    /// Per-component denominator receipt storage.
    DenominatorReceipt,
    /// Exact-field component receipt storage.
    ComponentReceipt,
}

/// Wilke model construction or evaluation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum WilkeGasViscosityErrorV1 {
    /// At least one component is required.
    EmptyComponents,
    /// The component count exceeded the fixed admission bound.
    TooManyComponents {
        /// Offered component count.
        offered: usize,
        /// Hard component limit.
        limit: usize,
    },
    /// Component and composition lengths differ.
    CompositionLengthMismatch {
        /// Component count.
        components: usize,
        /// Fraction count.
        fractions: usize,
    },
    /// The shared quantity-semantics boundary refused a composition.
    Composition(SemanticError),
    /// One immutable component contains an invalid physical scalar.
    InvalidComponentValue {
        /// Affected species.
        species: SpeciesId,
        /// Invalid component field.
        field: WilkeComponentFieldV1,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// The common caller-declared temperature is not positive and finite.
    InvalidTemperature {
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// A listed component has zero or negative-zero fraction and must be omitted.
    ZeroFraction {
        /// Canonical component index.
        component: usize,
        /// Affected species.
        species: SpeciesId,
        /// Exact zero bits.
        bits: u64,
    },
    /// Canonical active fractions do not sum to exact one.
    FractionSumNotExact {
        /// Fraction basis whose sum failed.
        basis: CompositionBasis,
        /// Exact fixed-order sum bits.
        sum_bits: u64,
    },
    /// Two active components name the same species.
    DuplicateSpecies {
        /// Repeated canonical species id.
        species: SpeciesId,
    },
    /// A positive declared component vanished during mass-to-mole conversion.
    UnrepresentableMoleFraction {
        /// Canonical component index.
        component: usize,
        /// Affected species.
        species: SpeciesId,
    },
    /// A bounded vector reservation failed before publication.
    AllocationRefused {
        /// Allocation stage.
        stage: WilkeAllocationStageV1,
        /// Requested element capacity.
        requested: usize,
    },
    /// Positive finite inputs produced a zero, negative, or non-finite intermediate.
    UnrepresentableArithmeticValue {
        /// Outer component `i` in canonical order.
        component: usize,
        /// Inner component `j` for pairwise work, if applicable.
        paired_component: Option<usize>,
        /// First failed operation-tree value.
        value: WilkeArithmeticValueV1,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
}

impl fmt::Display for WilkeGasViscosityErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyComponents => {
                formatter.write_str("Wilke viscosity mixing requires a component")
            }
            Self::TooManyComponents { offered, limit } => write!(
                formatter,
                "Wilke viscosity mixing offered {offered} components, exceeding limit {limit}"
            ),
            Self::CompositionLengthMismatch {
                components,
                fractions,
            } => write!(
                formatter,
                "Wilke viscosity mixing has {components} components but {fractions} fractions"
            ),
            Self::Composition(error) => {
                write!(formatter, "Wilke viscosity composition refused: {error}")
            }
            Self::InvalidComponentValue {
                species,
                field,
                bits,
            } => write!(
                formatter,
                "Wilke component {species} has invalid {field:?} bits {bits:#018x}"
            ),
            Self::InvalidTemperature { bits } => write!(
                formatter,
                "Wilke mixture temperature must be positive and finite (bits {bits:#018x})"
            ),
            Self::ZeroFraction {
                component,
                species,
                bits,
            } => write!(
                formatter,
                "Wilke component {component} ({species}) has zero fraction bits {bits:#018x}; omit absent species"
            ),
            Self::FractionSumNotExact { basis, sum_bits } => write!(
                formatter,
                "Wilke {basis:?} fractions must sum to exact 1.0 in canonical order (sum bits {sum_bits:#018x})"
            ),
            Self::DuplicateSpecies { species } => {
                write!(formatter, "Wilke mixture repeats active species {species}")
            }
            Self::UnrepresentableMoleFraction { component, species } => write!(
                formatter,
                "Wilke component {component} ({species}) lost its positive mole fraction during basis conversion"
            ),
            Self::AllocationRefused { stage, requested } => write!(
                formatter,
                "Wilke {stage:?} allocation refused for {requested} elements"
            ),
            Self::UnrepresentableArithmeticValue {
                component,
                paired_component,
                value,
                bits,
            } => write!(
                formatter,
                "Wilke component {component}, pair {paired_component:?}, {value:?} is not positive finite (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for WilkeGasViscosityErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Composition(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SemanticError> for WilkeGasViscosityErrorV1 {
    fn from(value: SemanticError) -> Self {
        Self::Composition(value)
    }
}

/// One immutable pure-species viscosity declaration.
///
/// The enclosing model declares that every component value applies at its one
/// retained temperature. This type validates scalar representation only; it
/// does not authenticate the species law, temperature, pressure, phase, or
/// transport regime from which the viscosity was obtained.
#[derive(Debug, Clone, PartialEq)]
pub struct WilkeGasViscosityComponentV1 {
    species: SpeciesId,
    molar_mass: MolarMass,
    dynamic_viscosity: DynViscosity,
}

impl WilkeGasViscosityComponentV1 {
    /// Admit one caller-supplied pure-species viscosity.
    ///
    /// # Errors
    /// Refuses non-positive or non-finite molar mass and dynamic viscosity.
    pub fn new(
        species: SpeciesId,
        molar_mass: MolarMass,
        dynamic_viscosity: DynViscosity,
    ) -> Result<Self, WilkeGasViscosityErrorV1> {
        for (field, value) in [
            (WilkeComponentFieldV1::MolarMass, molar_mass.value()),
            (
                WilkeComponentFieldV1::DynamicViscosity,
                dynamic_viscosity.value(),
            ),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(WilkeGasViscosityErrorV1::InvalidComponentValue {
                    species,
                    field,
                    bits: value.to_bits(),
                });
            }
        }
        Ok(Self {
            species,
            molar_mass,
            dynamic_viscosity,
        })
    }

    /// Canonical chemical species identity.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Positive finite molar mass.
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }

    /// Positive finite caller-supplied pure-species viscosity.
    #[must_use]
    pub const fn dynamic_viscosity(&self) -> DynViscosity {
        self.dynamic_viscosity
    }
}

/// Canonical bounded Wilke dilute-gas viscosity mixture.
///
/// Components are sorted by `SpeciesId`. Fractions pair with the input order,
/// then move with their component during canonicalization. Mass fractions are
/// converted through typed molar masses; volume fractions refuse. Every listed
/// component must be active because zero fractions create ambiguous receipt
/// identity and needless quadratic work.
#[derive(Debug, Clone, PartialEq)]
pub struct WilkeGasViscosityMixtureV1 {
    temperature: Temperature,
    components: Vec<WilkeGasViscosityComponentV1>,
    declared_composition: Composition,
    mole_composition: Composition,
}

impl WilkeGasViscosityMixtureV1 {
    /// Admit and canonicalize one bounded Wilke mixture.
    ///
    /// # Errors
    /// Count and length gates run before component scanning or sorting.
    /// Thereafter construction refuses invalid temperature, allocation,
    /// duplicate/zero/non-exact composition, unsupported volume basis, or an
    /// unrepresentable mass-to-mole conversion. No partial model escapes.
    #[allow(clippy::too_many_lines)]
    pub fn new(
        temperature: Temperature,
        components: Vec<WilkeGasViscosityComponentV1>,
        composition: Composition,
    ) -> Result<Self, WilkeGasViscosityErrorV1> {
        if components.is_empty() {
            return Err(WilkeGasViscosityErrorV1::EmptyComponents);
        }
        if components.len() > MAX_WILKE_GAS_COMPONENTS_V1 {
            return Err(WilkeGasViscosityErrorV1::TooManyComponents {
                offered: components.len(),
                limit: MAX_WILKE_GAS_COMPONENTS_V1,
            });
        }
        if components.len() != composition.fractions().len() {
            return Err(WilkeGasViscosityErrorV1::CompositionLengthMismatch {
                components: components.len(),
                fractions: composition.fractions().len(),
            });
        }
        let temperature_value = temperature.value();
        if !temperature_value.is_finite() || temperature_value <= 0.0 {
            return Err(WilkeGasViscosityErrorV1::InvalidTemperature {
                bits: temperature_value.to_bits(),
            });
        }

        let basis = composition.basis();
        let component_count = components.len();
        let mut paired = reserve_vec(component_count, WilkeAllocationStageV1::CanonicalComponents)?;
        for (component, fraction) in components
            .into_iter()
            .zip(composition.fractions().iter().copied())
        {
            paired.push((component, fraction));
        }
        paired.sort_unstable_by(|(left, _), (right, _)| left.species.cmp(&right.species));

        for pair in paired.windows(2) {
            if pair[0].0.species == pair[1].0.species {
                return Err(WilkeGasViscosityErrorV1::DuplicateSpecies {
                    species: pair[0].0.species.clone(),
                });
            }
        }
        for (component, (input, fraction)) in paired.iter().enumerate() {
            if *fraction <= 0.0 {
                return Err(WilkeGasViscosityErrorV1::ZeroFraction {
                    component,
                    species: input.species.clone(),
                    bits: fraction.to_bits(),
                });
            }
        }
        let declared_sum: f64 = paired.iter().map(|(_, fraction)| *fraction).sum();
        if declared_sum.to_bits() != 1.0_f64.to_bits() {
            return Err(WilkeGasViscosityErrorV1::FractionSumNotExact {
                basis,
                sum_bits: declared_sum.to_bits(),
            });
        }

        let mut declared_fractions =
            reserve_vec(component_count, WilkeAllocationStageV1::DeclaredFractions)?;
        let mut molar_masses = reserve_vec(component_count, WilkeAllocationStageV1::MolarMasses)?;
        for (component, fraction) in &paired {
            declared_fractions.push(*fraction);
            molar_masses.push(component.molar_mass);
        }
        let declared_composition = Composition::new(basis, declared_fractions)?;
        let mole_composition = declared_composition.to_mole_fractions(&molar_masses)?;
        for (component, (&declared, &mole)) in declared_composition
            .fractions()
            .iter()
            .zip(mole_composition.fractions())
            .enumerate()
        {
            if declared > 0.0 && mole <= 0.0 {
                return Err(WilkeGasViscosityErrorV1::UnrepresentableMoleFraction {
                    component,
                    species: paired[component].0.species.clone(),
                });
            }
        }
        let mole_sum: f64 = mole_composition.fractions().iter().sum();
        if mole_sum.to_bits() != 1.0_f64.to_bits() {
            return Err(WilkeGasViscosityErrorV1::FractionSumNotExact {
                basis: CompositionBasis::MoleFraction,
                sum_bits: mole_sum.to_bits(),
            });
        }

        let mut canonical_components =
            reserve_vec(component_count, WilkeAllocationStageV1::CanonicalComponents)?;
        for (component, _) in paired {
            canonical_components.push(component);
        }

        Ok(Self {
            temperature,
            components: canonical_components,
            declared_composition,
            mole_composition,
        })
    }

    /// Common caller-declared evaluation temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Active components in canonical species order.
    #[must_use]
    pub fn components(&self) -> &[WilkeGasViscosityComponentV1] {
        &self.components
    }

    /// Canonical composition in the caller-declared basis.
    #[must_use]
    pub const fn declared_composition(&self) -> &Composition {
        &self.declared_composition
    }

    /// Canonical mole fractions used by Wilke's operation tree.
    #[must_use]
    pub const fn mole_composition(&self) -> &Composition {
        &self.mole_composition
    }

    /// Evaluate Wilke's fixed-order dilute-gas viscosity mixing rule.
    ///
    /// Every intermediate must remain positive and finite. Pairwise work is
    /// computed directly in canonical `(i, j)` order without retaining an
    /// `N x N` matrix; only one denominator per component enters the receipt.
    ///
    /// # Errors
    /// Refuses bounded receipt allocation or the first unrepresentable
    /// operation-tree value. Nothing partial escapes.
    pub fn evaluate(&self) -> Result<WilkeGasViscosityEvaluationV1, WilkeGasViscosityErrorV1> {
        let count = self.components.len();
        let fractions = self.mole_composition.fractions();
        let mut denominator_bits = reserve_vec(count, WilkeAllocationStageV1::DenominatorReceipt)?;
        let mut mixture_viscosity = 0.0;

        for (i, (component_i, &fraction_i)) in self.components.iter().zip(fractions).enumerate() {
            let mut denominator = 0.0;
            for (j, (component_j, &fraction_j)) in self.components.iter().zip(fractions).enumerate()
            {
                let phi = interaction_factor(component_i, component_j, i, j)?;
                let weighted = require_positive_finite(
                    i,
                    Some(j),
                    WilkeArithmeticValueV1::WeightedInteraction,
                    fraction_j * phi,
                )?;
                denominator = require_positive_finite(
                    i,
                    Some(j),
                    WilkeArithmeticValueV1::ComponentDenominator,
                    denominator + weighted,
                )?;
            }
            denominator_bits.push(denominator.to_bits());
            let scaled = require_positive_finite(
                i,
                None,
                WilkeArithmeticValueV1::ScaledComponentViscosity,
                component_i.dynamic_viscosity.value() / denominator,
            )?;
            let contribution = require_positive_finite(
                i,
                None,
                WilkeArithmeticValueV1::ComponentContribution,
                fraction_i * scaled,
            )?;
            mixture_viscosity = require_positive_finite(
                i,
                None,
                WilkeArithmeticValueV1::MixtureViscosity,
                mixture_viscosity + contribution,
            )?;
        }

        let mut species = reserve_vec(count, WilkeAllocationStageV1::ComponentReceipt)?;
        let mut molar_mass_bits = reserve_vec(count, WilkeAllocationStageV1::ComponentReceipt)?;
        let mut pure_viscosity_bits = reserve_vec(count, WilkeAllocationStageV1::ComponentReceipt)?;
        let mut mole_fraction_bits = reserve_vec(count, WilkeAllocationStageV1::ComponentReceipt)?;
        let mut declared_fraction_bits =
            reserve_vec(count, WilkeAllocationStageV1::ComponentReceipt)?;
        for (index, component) in self.components.iter().enumerate() {
            species.push(component.species.clone());
            molar_mass_bits.push(component.molar_mass.value().to_bits());
            pure_viscosity_bits.push(component.dynamic_viscosity.value().to_bits());
            mole_fraction_bits.push(fractions[index].to_bits());
            declared_fraction_bits.push(self.declared_composition.fractions()[index].to_bits());
        }

        Ok(WilkeGasViscosityEvaluationV1 {
            dynamic_viscosity: DynViscosity::new(mixture_viscosity),
            receipt: WilkeGasViscosityReceiptV1 {
                evaluator_version: WILKE_GAS_VISCOSITY_EVALUATOR_VERSION_V1,
                fs_math_version: fs_math::VERSION,
                fs_qty_version: fs_qty::VERSION,
                rule: GasViscosityMixingRuleV1::Wilke,
                temperature_bits: self.temperature.value().to_bits(),
                declared_basis: self.declared_composition.basis(),
                declared_fraction_sum_bits: self
                    .declared_composition
                    .fractions()
                    .iter()
                    .sum::<f64>()
                    .to_bits(),
                mole_fraction_sum_bits: fractions.iter().sum::<f64>().to_bits(),
                species,
                molar_mass_bits,
                pure_viscosity_bits,
                declared_fraction_bits,
                mole_fraction_bits,
                denominator_bits,
                mixture_viscosity_bits: mixture_viscosity.to_bits(),
            },
        })
    }
}

fn interaction_factor(
    component_i: &WilkeGasViscosityComponentV1,
    component_j: &WilkeGasViscosityComponentV1,
    i: usize,
    j: usize,
) -> Result<f64, WilkeGasViscosityErrorV1> {
    let viscosity_ratio = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::ViscosityRatio,
        component_i.dynamic_viscosity.value() / component_j.dynamic_viscosity.value(),
    )?;
    let reverse_mass_ratio = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::ReverseMolarMassRatio,
        component_j.molar_mass.value() / component_i.molar_mass.value(),
    )?;
    let mass_quarter_root = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::MolarMassQuarterRoot,
        fs_math::det::sqrt(fs_math::det::sqrt(reverse_mass_ratio)),
    )?;
    let root_product = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::InteractionRootProduct,
        fs_math::det::sqrt(viscosity_ratio) * mass_quarter_root,
    )?;
    let numerator_base = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::InteractionNumeratorBase,
        1.0 + root_product,
    )?;
    let numerator = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::InteractionNumerator,
        numerator_base * numerator_base,
    )?;
    let forward_mass_ratio = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::ForwardMolarMassRatio,
        component_i.molar_mass.value() / component_j.molar_mass.value(),
    )?;
    let denominator_radicand = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::InteractionDenominatorRadicand,
        8.0 * (1.0 + forward_mass_ratio),
    )?;
    let denominator = require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::InteractionDenominator,
        fs_math::det::sqrt(denominator_radicand),
    )?;
    require_positive_finite(
        i,
        Some(j),
        WilkeArithmeticValueV1::InteractionFactor,
        numerator / denominator,
    )
}

fn require_positive_finite(
    component: usize,
    paired_component: Option<usize>,
    value_name: WilkeArithmeticValueV1,
    value: f64,
) -> Result<f64, WilkeGasViscosityErrorV1> {
    if !value.is_finite() || value <= 0.0 {
        return Err(WilkeGasViscosityErrorV1::UnrepresentableArithmeticValue {
            component,
            paired_component,
            value: value_name,
            bits: value.to_bits(),
        });
    }
    Ok(value)
}

fn reserve_vec<T>(
    capacity: usize,
    stage: WilkeAllocationStageV1,
) -> Result<Vec<T>, WilkeGasViscosityErrorV1> {
    let mut values = Vec::new();
    values.try_reserve_exact(capacity).map_err(|_| {
        WilkeGasViscosityErrorV1::AllocationRefused {
            stage,
            requested: capacity,
        }
    })?;
    Ok(values)
}

/// Immutable exact-field receipt for one Wilke mixture evaluation.
///
/// The receipt binds the complete supplied scalar data and fixed operation
/// tree. It is provenance for replay, not evidence that component properties
/// share a physical state or that the Wilke regime is valid.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WilkeGasViscosityReceiptV1 {
    evaluator_version: u32,
    fs_math_version: &'static str,
    fs_qty_version: &'static str,
    rule: GasViscosityMixingRuleV1,
    temperature_bits: u64,
    declared_basis: CompositionBasis,
    declared_fraction_sum_bits: u64,
    mole_fraction_sum_bits: u64,
    species: Vec<SpeciesId>,
    molar_mass_bits: Vec<u64>,
    pure_viscosity_bits: Vec<u64>,
    declared_fraction_bits: Vec<u64>,
    mole_fraction_bits: Vec<u64>,
    denominator_bits: Vec<u64>,
    mixture_viscosity_bits: u64,
}

impl WilkeGasViscosityReceiptV1 {
    /// Evaluator and operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Deterministic elementary-math crate version.
    #[must_use]
    pub const fn fs_math_version(&self) -> &'static str {
        self.fs_math_version
    }

    /// Quantity-semantics crate version.
    #[must_use]
    pub const fn fs_qty_version(&self) -> &'static str {
        self.fs_qty_version
    }

    /// Explicit viscosity mixing rule.
    #[must_use]
    pub const fn rule(&self) -> GasViscosityMixingRuleV1 {
        self.rule
    }

    /// Exact common caller-declared temperature bits.
    #[must_use]
    pub const fn temperature_bits(&self) -> u64 {
        self.temperature_bits
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

    /// Canonically ordered active species.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Exact canonical molar-mass bits.
    #[must_use]
    pub fn molar_mass_bits(&self) -> &[u64] {
        &self.molar_mass_bits
    }

    /// Exact canonical pure-species viscosity bits.
    #[must_use]
    pub fn pure_viscosity_bits(&self) -> &[u64] {
        &self.pure_viscosity_bits
    }

    /// Exact declared-fraction bits in canonical species order.
    #[must_use]
    pub fn declared_fraction_bits(&self) -> &[u64] {
        &self.declared_fraction_bits
    }

    /// Exact evaluated mole-fraction bits in canonical species order.
    #[must_use]
    pub fn mole_fraction_bits(&self) -> &[u64] {
        &self.mole_fraction_bits
    }

    /// Exact denominator bits for every canonical outer component.
    #[must_use]
    pub fn denominator_bits(&self) -> &[u64] {
        &self.denominator_bits
    }

    /// Exact mixed dynamic-viscosity bits.
    #[must_use]
    pub const fn mixture_viscosity_bits(&self) -> u64 {
        self.mixture_viscosity_bits
    }
}

/// One successful Wilke gas-mixture viscosity evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct WilkeGasViscosityEvaluationV1 {
    dynamic_viscosity: DynViscosity,
    receipt: WilkeGasViscosityReceiptV1,
}

impl WilkeGasViscosityEvaluationV1 {
    /// Mixed dynamic viscosity in coherent SI (`Pa s`).
    #[must_use]
    pub const fn dynamic_viscosity(&self) -> DynViscosity {
        self.dynamic_viscosity
    }

    /// Complete exact-field evaluation receipt.
    #[must_use]
    pub const fn receipt(&self) -> &WilkeGasViscosityReceiptV1 {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn component(name: &str, molar_mass: f64, viscosity: f64) -> WilkeGasViscosityComponentV1 {
        WilkeGasViscosityComponentV1::new(
            SpeciesId::new(name).expect("canonical test species"),
            MolarMass::new(molar_mass),
            DynViscosity::new(viscosity),
        )
        .expect("valid test component")
    }

    fn composition(basis: CompositionBasis, fractions: Vec<f64>) -> Composition {
        Composition::new(basis, fractions).expect("valid test composition")
    }

    fn model(
        components: Vec<WilkeGasViscosityComponentV1>,
        basis: CompositionBasis,
        fractions: Vec<f64>,
    ) -> WilkeGasViscosityMixtureV1 {
        WilkeGasViscosityMixtureV1::new(
            Temperature::new(900.0),
            components,
            composition(basis, fractions),
        )
        .expect("valid test Wilke model")
    }

    fn assert_close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
        assert!(
            (actual - expected).abs() <= 128.0 * f64::EPSILON * scale,
            "actual {actual:?} differs from expected {expected:?}",
        );
    }

    fn oracle_phi(mu_i: f64, mu_j: f64, mass_i: f64, mass_j: f64) -> f64 {
        let root = fs_math::det::sqrt((mu_i / mu_j) * fs_math::det::sqrt(mass_j / mass_i));
        (1.0 + root) * (1.0 + root)
            / (fs_math::det::sqrt(8.0) * fs_math::det::sqrt(1.0 + mass_i / mass_j))
    }

    #[test]
    fn g0_single_species_and_equal_component_limits_are_exact() {
        let pure = model(
            vec![component("Ar", 0.039_948, 4.2e-5)],
            CompositionBasis::MoleFraction,
            vec![1.0],
        )
        .evaluate()
        .expect("single species evaluates");
        assert_eq!(
            pure.dynamic_viscosity().value().to_bits(),
            4.2e-5f64.to_bits()
        );
        assert_eq!(pure.receipt().denominator_bits(), &[1.0f64.to_bits()]);

        let equal = model(
            vec![component("A", 0.02, 3.0e-5), component("B", 0.02, 3.0e-5)],
            CompositionBasis::MoleFraction,
            vec![0.25, 0.75],
        )
        .evaluate()
        .expect("equal properties evaluate");
        assert_eq!(
            equal.dynamic_viscosity().value().to_bits(),
            3.0e-5f64.to_bits()
        );
        assert_eq!(
            equal.receipt().denominator_bits(),
            &[1.0f64.to_bits(), 1.0f64.to_bits()]
        );
    }

    #[test]
    fn g0_binary_result_matches_independent_wilke_expression() {
        let mu_a = 1.8e-5;
        let mu_b = 4.1e-5;
        let mass_a = 0.028;
        let mass_b = 0.044;
        let x_a = 0.375;
        let x_b = 0.625;
        let evaluation = model(
            vec![component("A", mass_a, mu_a), component("B", mass_b, mu_b)],
            CompositionBasis::MoleFraction,
            vec![x_a, x_b],
        )
        .evaluate()
        .expect("binary mixture evaluates");

        let phi_aa = oracle_phi(mu_a, mu_a, mass_a, mass_a);
        let phi_ab = oracle_phi(mu_a, mu_b, mass_a, mass_b);
        let phi_ba = oracle_phi(mu_b, mu_a, mass_b, mass_a);
        let phi_bb = oracle_phi(mu_b, mu_b, mass_b, mass_b);
        let expected =
            x_a * mu_a / (x_a * phi_aa + x_b * phi_ab) + x_b * mu_b / (x_a * phi_ba + x_b * phi_bb);
        assert_close(evaluation.dynamic_viscosity().value(), expected);
    }

    #[test]
    fn g3_canonical_permutation_and_common_scaling_are_stable() {
        let canonical = model(
            vec![
                component("A", 0.02, 2.0e-5),
                component("B", 0.03, 3.0e-5),
                component("C", 0.04, 5.0e-5),
            ],
            CompositionBasis::MoleFraction,
            vec![0.25, 0.25, 0.5],
        )
        .evaluate()
        .expect("canonical mixture evaluates");
        let permuted = model(
            vec![
                component("C", 0.04, 5.0e-5),
                component("A", 0.02, 2.0e-5),
                component("B", 0.03, 3.0e-5),
            ],
            CompositionBasis::MoleFraction,
            vec![0.5, 0.25, 0.25],
        )
        .evaluate()
        .expect("permuted mixture evaluates");
        assert_eq!(canonical, permuted);

        let scaled = model(
            vec![
                component("A", 0.02, 8.0e-5),
                component("B", 0.03, 12.0e-5),
                component("C", 0.04, 20.0e-5),
            ],
            CompositionBasis::MoleFraction,
            vec![0.25, 0.25, 0.5],
        )
        .evaluate()
        .expect("scaled mixture evaluates");
        assert_close(
            scaled.dynamic_viscosity().value(),
            4.0 * canonical.dynamic_viscosity().value(),
        );
    }

    #[test]
    fn g0_mass_basis_uses_typed_molar_mass_conversion() {
        let from_mass = model(
            vec![component("A", 0.02, 2.0e-5), component("B", 0.02, 4.0e-5)],
            CompositionBasis::MassFraction,
            vec![0.25, 0.75],
        );
        let from_mole = model(
            vec![component("A", 0.02, 2.0e-5), component("B", 0.02, 4.0e-5)],
            CompositionBasis::MoleFraction,
            vec![0.25, 0.75],
        );
        assert_eq!(
            from_mass
                .evaluate()
                .expect("mass basis")
                .dynamic_viscosity(),
            from_mole
                .evaluate()
                .expect("mole basis")
                .dynamic_viscosity()
        );
        assert_eq!(
            from_mass.declared_composition().basis(),
            CompositionBasis::MassFraction
        );
    }

    #[test]
    fn g3_admission_refuses_invalid_and_ambiguous_inputs_in_order() {
        let invalid_mass = WilkeGasViscosityComponentV1::new(
            SpeciesId::new("A").expect("species"),
            MolarMass::new(0.0),
            DynViscosity::new(1.0),
        );
        assert!(matches!(
            invalid_mass,
            Err(WilkeGasViscosityErrorV1::InvalidComponentValue {
                field: WilkeComponentFieldV1::MolarMass,
                ..
            })
        ));
        let invalid_viscosity = WilkeGasViscosityComponentV1::new(
            SpeciesId::new("A").expect("species"),
            MolarMass::new(1.0),
            DynViscosity::new(f64::NAN),
        );
        assert!(matches!(
            invalid_viscosity,
            Err(WilkeGasViscosityErrorV1::InvalidComponentValue {
                field: WilkeComponentFieldV1::DynamicViscosity,
                ..
            })
        ));

        let valid_composition = composition(CompositionBasis::MoleFraction, vec![1.0]);
        assert_eq!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                Vec::new(),
                valid_composition.clone(),
            ),
            Err(WilkeGasViscosityErrorV1::EmptyComponents)
        );
        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                vec![component("A", 1.0, 1.0), component("B", 1.0, 1.0)],
                valid_composition.clone(),
            ),
            Err(WilkeGasViscosityErrorV1::CompositionLengthMismatch { .. })
        ));
        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(0.0),
                vec![component("A", 1.0, 1.0)],
                valid_composition,
            ),
            Err(WilkeGasViscosityErrorV1::InvalidTemperature { .. })
        ));

        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                vec![component("A", 1.0, 1.0), component("A", 1.0, 2.0)],
                composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]),
            ),
            Err(WilkeGasViscosityErrorV1::DuplicateSpecies { .. })
        ));
        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                vec![component("A", 1.0, 1.0), component("B", 1.0, 2.0)],
                composition(CompositionBasis::MoleFraction, vec![1.0, 0.0]),
            ),
            Err(WilkeGasViscosityErrorV1::ZeroFraction { component: 1, .. })
        ));
        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                vec![component("A", 1.0, 1.0), component("B", 1.0, 2.0)],
                composition(CompositionBasis::MoleFraction, vec![0.5, 0.5 + 1.0e-14],),
            ),
            Err(WilkeGasViscosityErrorV1::FractionSumNotExact { .. })
        ));
        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                vec![component("A", 1.0, 1.0)],
                composition(CompositionBasis::VolumeFraction, vec![1.0]),
            ),
            Err(WilkeGasViscosityErrorV1::Composition(_))
        ));

        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                vec![
                    component("A", f64::MIN_POSITIVE, 1.0),
                    component("B", f64::MAX, 2.0),
                ],
                composition(CompositionBasis::MassFraction, vec![0.5, 0.5]),
            ),
            Err(WilkeGasViscosityErrorV1::UnrepresentableMoleFraction { component: 1, .. })
        ));
    }

    #[test]
    fn g3_count_and_arithmetic_extremes_refuse_without_partial_result() {
        let mut components = Vec::new();
        let mut fractions = Vec::new();
        for index in 0..=MAX_WILKE_GAS_COMPONENTS_V1 {
            components.push(component(&format!("S{index:03}"), 1.0, 1.0));
            fractions.push(1.0 / 129.0);
        }
        assert!(matches!(
            WilkeGasViscosityMixtureV1::new(
                Temperature::new(300.0),
                components,
                composition(CompositionBasis::MoleFraction, fractions),
            ),
            Err(WilkeGasViscosityErrorV1::TooManyComponents {
                offered: 129,
                limit: MAX_WILKE_GAS_COMPONENTS_V1,
            })
        ));

        let hostile = model(
            vec![
                component("A", 1.0, f64::MAX),
                component("B", 1.0, f64::MIN_POSITIVE),
            ],
            CompositionBasis::MoleFraction,
            vec![0.5, 0.5],
        );
        assert!(matches!(
            hostile.evaluate(),
            Err(WilkeGasViscosityErrorV1::UnrepresentableArithmeticValue {
                component: 0,
                paired_component: Some(1),
                value: WilkeArithmeticValueV1::ViscosityRatio,
                ..
            })
        ));
    }

    #[test]
    fn g5_receipt_is_complete_and_replay_stable() {
        let mixture = model(
            vec![component("A", 0.02, 2.0e-5), component("B", 0.04, 4.0e-5)],
            CompositionBasis::MoleFraction,
            vec![0.25, 0.75],
        );
        let first = mixture.evaluate().expect("first replay");
        let second = mixture.evaluate().expect("second replay");
        assert_eq!(first, second);

        let receipt = first.receipt();
        assert_eq!(
            receipt.evaluator_version(),
            WILKE_GAS_VISCOSITY_EVALUATOR_VERSION_V1
        );
        assert_eq!(receipt.fs_math_version(), fs_math::VERSION);
        assert_eq!(receipt.fs_qty_version(), fs_qty::VERSION);
        assert_eq!(receipt.rule(), GasViscosityMixingRuleV1::Wilke);
        assert_eq!(receipt.temperature_bits(), 900.0f64.to_bits());
        assert_eq!(receipt.declared_basis(), CompositionBasis::MoleFraction);
        assert_eq!(receipt.declared_fraction_sum_bits(), 1.0f64.to_bits());
        assert_eq!(receipt.mole_fraction_sum_bits(), 1.0f64.to_bits());
        assert_eq!(receipt.species()[0].as_str(), "A");
        assert_eq!(receipt.species()[1].as_str(), "B");
        assert_eq!(
            receipt.molar_mass_bits(),
            &[0.02f64.to_bits(), 0.04f64.to_bits()]
        );
        assert_eq!(
            receipt.pure_viscosity_bits(),
            &[2.0e-5f64.to_bits(), 4.0e-5f64.to_bits()]
        );
        assert_eq!(
            receipt.declared_fraction_bits(),
            &[0.25f64.to_bits(), 0.75f64.to_bits()]
        );
        assert_eq!(
            receipt.mole_fraction_bits(),
            receipt.declared_fraction_bits()
        );
        assert_eq!(receipt.denominator_bits().len(), 2);
        assert_eq!(
            receipt.mixture_viscosity_bits(),
            first.dynamic_viscosity().value().to_bits()
        );

        let changed = model(
            vec![component("A", 0.02, 2.0e-5), component("B", 0.04, 4.1e-5)],
            CompositionBasis::MoleFraction,
            vec![0.25, 0.75],
        )
        .evaluate()
        .expect("changed input evaluates");
        assert_ne!(first.receipt(), changed.receipt());

        let changed_temperature = WilkeGasViscosityMixtureV1::new(
            Temperature::new(901.0),
            vec![component("A", 0.02, 2.0e-5), component("B", 0.04, 4.0e-5)],
            composition(CompositionBasis::MoleFraction, vec![0.25, 0.75]),
        )
        .expect("changed temperature admits")
        .evaluate()
        .expect("changed temperature evaluates");
        assert_eq!(
            first.dynamic_viscosity(),
            changed_temperature.dynamic_viscosity()
        );
        assert_ne!(first.receipt(), changed_temperature.receipt());
    }

    #[test]
    fn g0_component_cap_boundary_is_admitted_before_quadratic_work() {
        let mut components = Vec::new();
        let mut fractions = Vec::new();
        for index in 0..MAX_WILKE_GAS_COMPONENTS_V1 {
            components.push(component(&format!("S{index:03}"), 0.02, 3.0e-5));
            fractions.push(1.0 / 128.0);
        }
        let evaluation = WilkeGasViscosityMixtureV1::new(
            Temperature::new(300.0),
            components,
            composition(CompositionBasis::MoleFraction, fractions),
        )
        .expect("cap boundary admits")
        .evaluate()
        .expect("cap boundary evaluates");
        assert_eq!(evaluation.receipt().species().len(), 128);
        assert_eq!(evaluation.receipt().denominator_bits().len(), 128);
        assert_close(evaluation.dynamic_viscosity().value(), 3.0e-5);
    }
}
