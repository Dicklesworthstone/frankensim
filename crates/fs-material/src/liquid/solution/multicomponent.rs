//! Homogeneous multicomponent liquids with their original material sources.
//!
//! A species-bound source chart uses N or N-1 explicit mole-fraction coordinates.
//! Any dependent species is declared, never inferred from a material label.
//! All-explicit charts preserve traces without subtracting from the carrier.
//! The same homogeneous composition supplies regular-solution thermodynamics
//! and every pinned transport claim, including in a spinodal region. This is
//! not a multicomponent flash or a mixing rule for pure-fluid properties.
//! The caloric adapter binds explicit pure-liquid sources for homogeneous
//! thermal mixing and the response of an accepted local phase candidate.
//! Pressure and all other source coordinates remain fixed.

use std::sync::Arc;

/// Original pure-liquid Cp sources for heating, mixing and local phase response.
pub mod caloric;
/// Original material properties resolved separately for accepted liquid phases.
pub mod equilibrium;

use crate::MaterialError;
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::liquid::ResolvedLiquidState;
use crate::state_point::MaterialStatePointError;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{MolarMass, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::multicomponent::{
    MAX_REGULAR_SOLUTION_SPECIES_V1, MulticomponentRegularSolutionV1,
    MulticomponentSolutionErrorV1, MulticomponentSolutionStateV1,
};

/// Explicit chemical identity, mass conversion and source-coordinate role.
/// At most one component has no axis and is the dependent mole fraction.
/// All N axes may instead be explicit, retaining trace source coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidComponent {
    pub species: SpeciesId,
    pub molar_mass: MolarMass,
    pub mole_fraction_axis: Option<String>,
}

/// Invalid source chart, unrepresentable inventory, or original-owner refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentLiquidError {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Source(MaterialStatePointError),
    Thermodynamics(MulticomponentSolutionErrorV1),
    Cancelled,
}
impl core::fmt::Display for MulticomponentLiquidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "homogeneous multicomponent liquid: {self:?}")
    }
}
impl std::error::Error for MulticomponentLiquidError {}
impl From<MaterialStatePointError> for MulticomponentLiquidError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
impl From<MulticomponentSolutionErrorV1> for MulticomponentLiquidError {
    fn from(error: MulticomponentSolutionErrorV1) -> Self {
        Self::Thermodynamics(error)
    }
}
type Result<T> = core::result::Result<T, MulticomponentLiquidError>;

#[derive(Clone, Debug, PartialEq)]
struct SourceOwner {
    thermodynamics: MulticomponentRegularSolutionV1,
    components: Box<[MulticomponentLiquidComponent]>,
    reference: ResolvedLiquidState,
    temperature_axis: String,
    dependent_component: Option<usize>,
}

/// Original mixture sources and a species-bound homogeneous free-energy law.
///
/// Bindings are canonicalized by species along with their molar masses and
/// coordinate names. Inventory inputs subsequently use `species()` order.
/// Chemical qualification of the mixture card and its pair coefficients is
/// the caller's responsibility; a label or matching dimensions cannot supply it.
/// With N-1 axes, prefer an abundant dependent component: its fraction must be
/// recoverable within bounded relative precision. Use N explicit coordinates
/// when no component stays abundant across all physical phases.
pub struct MulticomponentLiquidSolutionLaw<'card> {
    card: &'card MaterialCard,
    point: QueryPoint,
    owner: Arc<SourceOwner>,
}
impl<'card> MulticomponentLiquidSolutionLaw<'card> {
    /// Bind one original liquid-property bundle and N or N-1 distinct, explicitly
    /// typed mole-fraction axes. All species must occur exactly once. The
    /// initial source composition and temperature must be in the law's domain.
    pub fn try_new(
        thermodynamics: MulticomponentRegularSolutionV1,
        mut components: Vec<MulticomponentLiquidComponent>,
        card: &'card MaterialCard,
        reference: ResolvedLiquidState,
        temperature_axis: &str,
    ) -> Result<Self> {
        if reference.material().card_identity() != card.content_hash() {
            return Err(MulticomponentLiquidError::Invalid(
                "liquid solution requires the original material card",
            ));
        }
        components.sort_by(|a, b| a.species.cmp(&b.species));
        if !components
            .iter()
            .map(|c| &c.species)
            .eq(thermodynamics.species())
            || components.iter().any(|c| !positive(c.molar_mass.value()))
        {
            return Err(MulticomponentLiquidError::Invalid(
                "each original species requires its own positive finite molar mass",
            ));
        }
        let point = reference.material().query();
        crate::liquid_temperature::validate_temperature(&point, temperature_axis)?;
        let mole_fraction = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Composition(CompositionBasis::MoleFraction),
            ValueForm::Static,
        ));
        let mut dependent = None;
        let mut fractions = vec![0.0; components.len()];
        for (i, component) in components.iter().enumerate() {
            let Some(axis) = &component.mole_fraction_axis else {
                if dependent.replace(i).is_some() {
                    return Err(MulticomponentLiquidError::Invalid(
                        "at most one dependent mole-fraction component is permitted",
                    ));
                }
                continue;
            };
            if axis == temperature_axis
                || components[..i]
                    .iter()
                    .any(|previous| previous.mole_fraction_axis.as_ref() == Some(axis))
                || point.axis_quantities().get(axis) != Some(&mole_fraction)
                || !point
                    .axes()
                    .get(axis)
                    .is_some_and(|x| (0.0..=1.0).contains(x))
            {
                return Err(MulticomponentLiquidError::Invalid(
                    "distinct original mole-fraction source coordinates are required",
                ));
            }
            fractions[i] = point.axes()[axis];
        }
        if point.axis_quantities().iter().any(|(axis, quantity)| {
            *quantity == mole_fraction
                && !components
                    .iter()
                    .any(|component| component.mole_fraction_axis.as_ref() == Some(axis))
        }) {
            return Err(MulticomponentLiquidError::Invalid(
                "every mole-fraction source coordinate must name an independent component",
            ));
        }
        if let Some(dependent) = dependent {
            let independent_sum =
                sum_nonnegative(fractions.iter().copied(), "reference mole-fraction sum")?;
            if independent_sum > 1.0 {
                return Err(MulticomponentLiquidError::Invalid(
                    "reference mole fractions exceed the composition simplex",
                ));
            }
            fractions[dependent] = 1.0 - independent_sum;
        }
        thermodynamics.evaluate(
            Temperature::new(point.axes()[temperature_axis]),
            &fractions,
            || false,
        )?;
        Ok(Self {
            card,
            point,
            owner: Arc::new(SourceOwner {
                thermodynamics,
                components: components.into_boxed_slice(),
                reference,
                temperature_axis: temperature_axis.to_owned(),
                dependent_component: dependent,
            }),
        })
    }

    /// Resolve one finite homogeneous inventory, including metastable and
    /// spinodal states, without replacing it by a bulk phase equilibrium.
    /// Every component amount is in mol in the exact `species()` order.
    /// Source refusal, cancellation or loss of a positive trace is atomic.
    pub fn homogeneous(
        &self,
        temperature: Temperature,
        component_moles: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentLiquidCell> {
        poll(&mut cancelled)?;
        if component_moles.len() != self.species().len() {
            return Err(MulticomponentLiquidError::Invalid(
                "one component amount per original species is required",
            ));
        }
        let total_moles = positive_value(
            sum_nonnegative(component_moles.iter().copied(), "total component amount")?,
            "positive total component amount",
        )?;
        let mut fractions = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for (out, &amount) in fractions.iter_mut().zip(component_moles) {
            *out = amount / total_moles;
        }
        let thermodynamics = self.owner.thermodynamics.evaluate(
            temperature,
            &fractions[..component_moles.len()],
            &mut cancelled,
        )?;
        let fractions = thermodynamics.mole_fractions();
        if component_moles
            .iter()
            .zip(fractions)
            .any(|(n, x)| *n > 0.0 && *x <= 0.0)
        {
            return Err(MulticomponentLiquidError::Arithmetic(
                "positive component vanished in mole-fraction normalization",
            ));
        }
        let mut point = self
            .point
            .clone()
            .with(&self.owner.temperature_axis, temperature.value())
            .map_err(|_| {
                MulticomponentLiquidError::Invalid("invalid source temperature coordinate")
            })?;
        for (component, &fraction) in self.owner.components.iter().zip(fractions) {
            if let Some(axis) = &component.mole_fraction_axis {
                point = point.with(axis, fraction).map_err(|_| {
                    MulticomponentLiquidError::Invalid("invalid source mole-fraction coordinate")
                })?;
            }
        }
        if let Some(dependent) = self.owner.dependent_component {
            let source_remainder = 1.0
                - sum_nonnegative(
                    fractions
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != dependent)
                        .map(|(_, &x)| x),
                    "source mole-fraction sum",
                )?;
            let dependent_fraction = fractions[dependent];
            let source_matches = if dependent_fraction == 0.0 {
                source_remainder == 0.0
            } else {
                source_remainder > 0.0
                    && (source_remainder - dependent_fraction).abs()
                        <= 64.0 * f64::EPSILON * dependent_fraction
            };
            if !source_matches {
                return Err(MulticomponentLiquidError::Arithmetic(
                    "dependent component is not accurately representable in the source composition chart",
                ));
            }
        }
        poll(&mut cancelled)?;
        let liquid = self.owner.reference.at_query_point(self.card, &point)?;
        let mut component_mass: Vec<_> = component_moles
            .iter()
            .zip(&self.owner.components)
            .map(|(&amount, component)| {
                let mass = amount * component.molar_mass.value();
                if !mass.is_finite() || (amount > 0.0 && mass <= 0.0) {
                    Err(MulticomponentLiquidError::Arithmetic("component mass"))
                } else {
                    Ok(mass)
                }
            })
            .collect::<Result<_>>()?;
        let mass_kg = positive_value(
            sum_nonnegative(component_mass.iter().copied(), "total mass")?,
            "total mass",
        )?;
        let volume_m3 = positive_value(mass_kg / liquid.density_kg_m3(), "liquid volume")?;
        let molar_concentration_mol_m3 =
            positive_value(total_moles / volume_m3, "molar concentration")?;
        for mass in &mut component_mass {
            let fraction = *mass / mass_kg;
            if *mass > 0.0 && fraction == 0.0 {
                return Err(MulticomponentLiquidError::Arithmetic(
                    "positive component mass fraction",
                ));
            }
            *mass = fraction;
        }
        poll(&mut cancelled)?;
        Ok(ResolvedMulticomponentLiquidCell {
            owner: Arc::clone(&self.owner),
            thermodynamics,
            liquid,
            component_moles: component_moles.into(),
            mass_fractions: component_mass.into_boxed_slice(),
            mass_kg,
            volume_m3,
            molar_concentration_mol_m3,
        })
    }

    pub fn species(&self) -> &[SpeciesId] {
        self.owner.thermodynamics.species()
    }
    pub fn thermodynamics(&self) -> &MulticomponentRegularSolutionV1 {
        &self.owner.thermodynamics
    }
    pub fn components(&self) -> &[MulticomponentLiquidComponent] {
        &self.owner.components
    }
    pub fn reference_state(&self) -> &ResolvedLiquidState {
        &self.owner.reference
    }
}

/// Original homogeneous thermodynamics and actual source-resolved transport.
/// A negative simplex curvature remains available through `thermodynamics()`;
/// it is not silently discarded and does not assert equilibrium phase amounts.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMulticomponentLiquidCell {
    owner: Arc<SourceOwner>,
    thermodynamics: MulticomponentSolutionStateV1,
    liquid: ResolvedLiquidState,
    component_moles: Box<[f64]>,
    mass_fractions: Box<[f64]>,
    mass_kg: f64,
    volume_m3: f64,
    molar_concentration_mol_m3: f64,
}
impl ResolvedMulticomponentLiquidCell {
    pub const fn thermodynamics(&self) -> &MulticomponentSolutionStateV1 {
        &self.thermodynamics
    }
    pub const fn liquid(&self) -> &ResolvedLiquidState {
        &self.liquid
    }
    pub fn species(&self) -> &[SpeciesId] {
        self.owner.thermodynamics.species()
    }
    pub fn components(&self) -> &[MulticomponentLiquidComponent] {
        &self.owner.components
    }
    pub fn mole_fractions(&self) -> &[f64] {
        self.thermodynamics.mole_fractions()
    }
    pub fn component_moles(&self) -> &[f64] {
        &self.component_moles
    }
    pub fn mass_fractions(&self) -> &[f64] {
        &self.mass_fractions
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Source-derived volume; no pure-component additive-volume assumption.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Actual concentration for the original Maxwell-Stefan force path.
    pub const fn molar_concentration_mol_m3(&self) -> f64 {
        self.molar_concentration_mol_m3
    }
}
impl FluidTransport for ResolvedMulticomponentLiquidCell {
    fn transport_state(&self) -> core::result::Result<FluidTransportState, MaterialError> {
        self.liquid.transport_state()
    }
}
impl ThermalFluidTransport for ResolvedMulticomponentLiquidCell {
    fn thermal_transport_state(
        &self,
    ) -> core::result::Result<ThermalFluidTransportState, MaterialError> {
        ThermalFluidTransportState::try_new(
            self.transport_state()?,
            self.thermodynamics.temperature().value(),
            self.liquid.volumetric_expansion_per_k(),
        )
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn positive_value(value: f64, name: &'static str) -> Result<f64> {
    if positive(value) {
        Ok(value)
    } else {
        Err(MulticomponentLiquidError::Arithmetic(name))
    }
}
fn sum_nonnegative(values: impl IntoIterator<Item = f64>, name: &'static str) -> Result<f64> {
    // Neumaier accumulation keeps many small positive components from depending
    // on their position beside the largest one; logical order is always fixed.
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        if !value.is_finite() || value < 0.0 {
            return Err(MulticomponentLiquidError::Invalid(
                "finite nonnegative component values required",
            ));
        }
        let next = sum + value;
        correction += if sum >= value {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    let total = sum + correction;
    if total.is_finite() {
        Ok(total)
    } else {
        Err(MulticomponentLiquidError::Arithmetic(name))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentLiquidError::Cancelled)
    } else {
        Ok(())
    }
}
