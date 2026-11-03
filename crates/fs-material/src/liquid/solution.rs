//! Binary liquid equilibrium with the original composition-dependent sources.
//!
//! The L1 regular-solution law selects one homogeneous liquid or two coexisting
//! compositions. This adapter closes their component inventories and evaluates
//! the SAME selected material claims at each resulting composition and
//! temperature. Density, viscosity, heat capacity and conductivity are never
//! inferred by averaging pure-component properties. In particular, phase mole
//! fractions are not phase mass fractions or volume fractions.
//!
//! This is isothermal bulk liquid equilibrium at the source's fixed pressure.
//! The caller supplies compatible pair parameters, component molar masses and
//! mixture-property data. It does not infer chemical identity from a material
//! name, include vapor/solid phases, solve a heat balance, or evolve nucleation,
//! diffusion, entrainment or interfaces. Each phase can enter existing flow,
//! heat and buoyancy consumers through the shared transport traits. There is
//! deliberately no single-fluid transport implementation for a two-phase result.

use super::ResolvedLiquidState;
use crate::MaterialError;
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::state_point::MaterialStatePointError;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{MolarMass, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::{
    BinaryCoexistenceOptionsV1, BinaryPhaseEquilibriumV1, BinaryRegularSolutionV1,
    BinarySolutionErrorV1,
};

/// Invalid inventory, thermodynamic refusal, or loss of the original source.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryLiquidError {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Source(MaterialStatePointError),
    Thermodynamics(BinarySolutionErrorV1),
    Cancelled,
}

impl core::fmt::Display for BinaryLiquidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "binary liquid solution: {self:?}")
    }
}
impl std::error::Error for BinaryLiquidError {}
impl From<MaterialStatePointError> for BinaryLiquidError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
impl From<BinarySolutionErrorV1> for BinaryLiquidError {
    fn from(error: BinarySolutionErrorV1) -> Self {
        Self::Thermodynamics(error)
    }
}

/// One pair law bound to a mixture material card and its original claim choices.
///
/// Component order is the thermodynamic owner's order; the composition axis is
/// explicitly the MOLE fraction of its first species. Molar masses and source
/// chemistry must describe those same species. Every other source coordinate,
/// including pressure, remains frozen. A new evaluation cannot switch claims.
pub struct BinaryLiquidSolutionLaw<'a> {
    thermodynamics: BinaryRegularSolutionV1,
    molar_masses: [MolarMass; 2],
    card: &'a MaterialCard,
    reference: ResolvedLiquidState,
    point: QueryPoint,
    temperature_axis: String,
    composition_axis: String,
}

impl<'a> BinaryLiquidSolutionLaw<'a> {
    /// Bind a previously resolved mixture-property bundle, retaining all receipts.
    /// An untyped or mass/volume-fraction composition coordinate is refused.
    pub fn try_new(
        thermodynamics: BinaryRegularSolutionV1,
        molar_masses: [MolarMass; 2],
        card: &'a MaterialCard,
        reference: ResolvedLiquidState,
        temperature_axis: &str,
        composition_axis: &str,
    ) -> Result<Self, BinaryLiquidError> {
        if reference.material().card_identity() != card.content_hash() {
            return Err(BinaryLiquidError::Invalid(
                "mixture requires the original material card",
            ));
        }
        if molar_masses.iter().any(|m| !positive(m.value())) {
            return Err(BinaryLiquidError::Invalid(
                "component molar masses must be positive finite kg/mol",
            ));
        }
        let point = reference.material().query();
        crate::liquid_temperature::validate_temperature(&point, temperature_axis)?;
        let mole_fraction = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Composition(CompositionBasis::MoleFraction),
            ValueForm::Static,
        ));
        if temperature_axis == composition_axis
            || point.axis_quantities().get(composition_axis) != Some(&mole_fraction)
            || !point
                .axes()
                .get(composition_axis)
                .is_some_and(|x| (0.0..=1.0).contains(x))
        {
            return Err(BinaryLiquidError::Invalid(
                "composition axis must explicitly be first-species mole fraction",
            ));
        }
        thermodynamics.evaluate(
            Temperature::new(point.axes()[temperature_axis]),
            point.axes()[composition_axis],
        )?;
        Ok(Self {
            thermodynamics,
            molar_masses,
            card,
            reference,
            point,
            temperature_axis: temperature_axis.to_owned(),
            composition_axis: composition_axis.to_owned(),
        })
    }

    /// Equilibrate a closed, finite pair inventory at a prescribed temperature.
    ///
    /// Component amounts are in mol in the thermodynamic owner's species order.
    /// Both conserved amounts survive cooling/heating and phase-count changes.
    /// Every returned phase has its own source-resolved density and transport.
    /// Work is bounded by the coexistence controls and at most two source
    /// refreshes. Cancellation or a failed second phase publishes no result.
    pub fn equilibrate(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryLiquidEquilibrium, BinaryLiquidError> {
        poll(&mut cancelled)?;
        if component_moles.iter().any(|n| !n.is_finite() || *n < 0.0) {
            return Err(BinaryLiquidError::Invalid(
                "component amounts must be finite nonnegative mol",
            ));
        }
        let total_moles = positive_value(component_moles[0] + component_moles[1], "total amount")?;
        // Normalize the smaller component directly, then complement only the
        // larger one. Both fractions reach the thermodynamic owner: deriving
        // a trace as 1-x would destroy its relative accuracy and chemical
        // potential, even while its original inventory remained conserved.
        let fractions = if component_moles[0] <= component_moles[1] {
            let x = component_moles[0] / total_moles;
            [x, 1.0 - x]
        } else {
            let y = component_moles[1] / total_moles;
            [1.0 - y, y]
        };
        if (component_moles[0] > 0.0 && fractions[0] == 0.0)
            || (component_moles[1] > 0.0 && fractions[0] == 1.0)
        {
            return Err(BinaryLiquidError::Arithmetic(
                "nonzero component lost in the source's first-species coordinate",
            ));
        }
        let equilibrium = self.thermodynamics.equilibrate_composition(
            temperature,
            fractions,
            options,
            &mut cancelled,
        )?;
        let mut inventories = [[0.0; 2]; 2];
        let phases = equilibrium.phases();
        if phases.len() == 1 {
            inventories[0] = component_moles;
        } else {
            // Compute the smaller allocation directly and the larger as a
            // complement, separately for each species. This retains a dilute
            // phase without subtracting nearly equal inventories and closes
            // both balances with the original feed, to floating-point roundoff.
            let amounts = [
                positive_value(
                    total_moles * phases[0].fraction_of_total_moles(),
                    "first phase amount",
                )?,
                positive_value(
                    total_moles * phases[1].fraction_of_total_moles(),
                    "second phase amount",
                )?,
            ];
            for (species, &feed) in component_moles.iter().enumerate() {
                let a = positive_value(
                    amounts[0] * phases[0].mole_fractions()[species],
                    "first phase species amount",
                )?;
                let b = positive_value(
                    amounts[1] * phases[1].mole_fractions()[species],
                    "second phase species amount",
                )?;
                let (first, second) = if a <= b { (a, feed - a) } else { (feed - b, b) };
                inventories[0][species] = positive_value(first, "first conserved species amount")?;
                inventories[1][species] =
                    positive_value(second, "second conserved species amount")?;
            }
        }
        let mut resolved = Vec::with_capacity(phases.len());
        for (phase, amounts) in phases.iter().zip(inventories) {
            poll(&mut cancelled)?;
            let point = self
                .point
                .clone()
                .with(&self.temperature_axis, temperature.value())
                .and_then(|point| point.with(&self.composition_axis, phase.mole_fraction_a()))
                .map_err(|_| {
                    BinaryLiquidError::Invalid("equilibrium violates source coordinate conventions")
                })?;
            let liquid = self.reference.at_query_point(self.card, &point)?;
            let component_mass = self.component_masses(amounts)?;
            let mass = positive_value(component_mass[0] + component_mass[1], "phase mass")?;
            let volume = positive_value(mass / liquid.density_kg_m3(), "phase volume")?;
            let moles = positive_value(amounts[0] + amounts[1], "phase amount")?;
            let concentration = positive_value(moles / volume, "phase molar concentration")?;
            resolved.push(ResolvedBinaryLiquidPhase {
                liquid,
                temperature,
                species: self.thermodynamics.species().clone(),
                mole_fractions: phase.mole_fractions(),
                component_moles: amounts,
                mass_fractions: [component_mass[0] / mass, component_mass[1] / mass],
                mass_kg: mass,
                volume_m3: volume,
                molar_concentration_mol_m3: concentration,
            });
        }
        let mass = self.component_masses(component_moles)?;
        let total_mass = positive_value(mass[0] + mass[1], "total mass")?;
        let volume = positive_value(
            resolved.iter().map(|p| p.volume_m3).sum(),
            "total liquid volume",
        )?;
        poll(&mut cancelled)?;
        Ok(ResolvedBinaryLiquidEquilibrium {
            equilibrium,
            phases: resolved,
            temperature,
            component_moles,
            total_mass_kg: total_mass,
            total_volume_m3: volume,
        })
    }

    /// The explicit thermodynamic pair law; no material-name inference.
    pub const fn thermodynamics(&self) -> &BinaryRegularSolutionV1 {
        &self.thermodynamics
    }
    /// Original source state and selected property-use receipts.
    pub const fn reference_state(&self) -> &ResolvedLiquidState {
        &self.reference
    }
    /// Component molar masses in the SAME order as the thermodynamic species.
    pub const fn molar_masses(&self) -> [MolarMass; 2] {
        self.molar_masses
    }

    fn component_masses(&self, amounts: [f64; 2]) -> Result<[f64; 2], BinaryLiquidError> {
        let mut mass = [0.0; 2];
        for i in 0..2 {
            mass[i] = amounts[i] * self.molar_masses[i].value();
            if !mass[i].is_finite() || (amounts[i] > 0.0 && mass[i] <= 0.0) {
                return Err(BinaryLiquidError::Arithmetic("component mass"));
            }
        }
        Ok(mass)
    }
}

/// One homogeneous equilibrium liquid with a conserved component inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryLiquidPhase {
    liquid: ResolvedLiquidState,
    temperature: Temperature,
    species: [SpeciesId; 2],
    mole_fractions: [f64; 2],
    component_moles: [f64; 2],
    mass_fractions: [f64; 2],
    mass_kg: f64,
    volume_m3: f64,
    molar_concentration_mol_m3: f64,
}

impl ResolvedBinaryLiquidPhase {
    /// Original ordered chemical identities for all component amounts and fractions.
    pub const fn species(&self) -> &[SpeciesId; 2] {
        &self.species
    }
    /// The original source-resolved liquid at this phase's actual composition.
    pub const fn liquid(&self) -> &ResolvedLiquidState {
        &self.liquid
    }
    /// Equilibrium compositions in thermodynamic species order, on a molar basis.
    pub const fn mole_fractions(&self) -> [f64; 2] {
        self.mole_fractions
    }
    /// Conserved component amounts [mol], in the same species order.
    pub const fn component_moles(&self) -> [f64; 2] {
        self.component_moles
    }
    /// Mass fractions derived using the declared component molar masses.
    pub const fn mass_fractions(&self) -> [f64; 2] {
        self.mass_fractions
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Actual source-derived volume; no additive-volume assumption is made.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Concentration usable by the existing liquid Maxwell-Stefan force path.
    pub const fn molar_concentration_mol_m3(&self) -> f64 {
        self.molar_concentration_mol_m3
    }
}

impl FluidTransport for ResolvedBinaryLiquidPhase {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.liquid.transport_state()
    }
}
impl ThermalFluidTransport for ResolvedBinaryLiquidPhase {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        ThermalFluidTransportState::try_new(
            self.transport_state()?,
            self.temperature.value(),
            self.liquid.volumetric_expansion_per_k(),
        )
    }
}

/// A complete one- or two-phase result. Sources and inventories remain together.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryLiquidEquilibrium {
    equilibrium: BinaryPhaseEquilibriumV1,
    phases: Vec<ResolvedBinaryLiquidPhase>,
    temperature: Temperature,
    component_moles: [f64; 2],
    total_mass_kg: f64,
    total_volume_m3: f64,
}

impl ResolvedBinaryLiquidEquilibrium {
    /// Original ordered chemical identities, retained through splitting and remixing.
    pub fn species(&self) -> &[SpeciesId; 2] {
        self.phases[0].species()
    }
    /// One homogeneous phase, or two in increasing first-species mole fraction.
    /// This is composition ordering; hydrostatic consumers must order by density.
    pub fn phases(&self) -> &[ResolvedBinaryLiquidPhase] {
        &self.phases
    }
    pub const fn thermodynamic_equilibrium(&self) -> &BinaryPhaseEquilibriumV1 {
        &self.equilibrium
    }
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn total_component_moles(&self) -> [f64; 2] {
        self.component_moles
    }
    pub const fn total_mass_kg(&self) -> f64 {
        self.total_mass_kg
    }
    pub const fn total_volume_m3(&self) -> f64 {
        self.total_volume_m3
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn positive_value(value: f64, quantity: &'static str) -> Result<f64, BinaryLiquidError> {
    if positive(value) {
        Ok(value)
    } else {
        Err(BinaryLiquidError::Arithmetic(quantity))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BinaryLiquidError> {
    if cancelled() {
        Err(BinaryLiquidError::Cancelled)
    } else {
        Ok(())
    }
}
