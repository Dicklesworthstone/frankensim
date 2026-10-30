//! One species/temperature/composition evaluation for acoustics and fluid heat transport.
//!
//! NASA-9 thermodynamics, the ideal-gas EOS and Wilke viscosity remain owned by
//! `fs-thermochem`. Pure-species viscosity and conductivity are queried through
//! the existing material-card resolver at exactly the same state. Conductivity
//! uses Wassiljewa--Mason--Saxena with the already evaluated Wilke denominators.
//! No calorically-perfect gamma, air preset or second species registry is used.
//! WMS reference: Mason and Saxena, Phys. Fluids 1, 361 (1958),
//! doi:10.1063/1.1724352. The shared viscosity interaction is retained explicitly.
//!
//! This is a frozen-composition, dilute ideal gas with local-equilibrium heat
//! capacities. Sound speed is the frozen-composition adiabatic derivative;
//! chemical-equilibrium sound speed, molecular relaxation, condensation and
//! real-gas effects are NOT modeled. The legacy USSA gas evaluator is unchanged.

/// Conservative density/energy reconstruction and closed-volume gas heating.
pub mod energy;
/// Explicit nonideal mechanical/phase binding over the same canonical species.
pub mod nonideal;
/// Explicit binding of verified species, NASA-9 and transport packs.
pub mod packs;

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::{IdealGasEosEvaluationV1, IdealGasEosV1};
use fs_thermochem::mixture::{FrozenIdealGasMixtureEvaluationV1, FrozenIdealGasMixtureModelV1};
use fs_thermochem::transport::{
    WilkeGasViscosityComponentV1, WilkeGasViscosityEvaluationV1, WilkeGasViscosityMixtureV1,
};

use crate::MaterialError;
use crate::fluid_transport::{FluidTransport, FluidTransportState};
use crate::gas::GasState;
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};

/// Explicit assignment of pure-species transport data to a canonical species.
/// The assignment is caller authority; a material name does not identify a gas.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesTransportMaterial {
    /// Species identity shared with the NASA-9 mixture, not a new name system.
    pub species: SpeciesId,
    /// Source-backed pure viscosity [Pa s] and conductivity [W/(m K)].
    pub material: MaterialCard,
    /// Explicit source-selection policy for this species, including any pinned claims.
    pub selection: MaterialPropertySelection,
}

/// Re-evaluable physical medium. Changing composition or source data constructs
/// a new model; changing temperature/pressure evaluates the same source laws.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermochemicalGasModel {
    thermodynamics: FrozenIdealGasMixtureModelV1,
    transport: Vec<SpeciesTransportMaterial>,
    water_component: Option<usize>,
    source_packs: Vec<packs::ThermochemicalPackReceipt>,
}

fn invalid(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn failure(stage: &str, error: impl core::fmt::Display) -> MaterialError {
    invalid(format!("thermochemical gas {stage}: {error}"))
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

/// The same explicitly typed coordinate convention as source-backed GasState.
fn coordinates(point: &QueryPoint) -> Result<(Temperature, Pressure), MaterialError> {
    let coordinate = |axis: &str, kind| {
        if point.axis_quantities().get(axis) != Some(&typed(kind)) {
            return Err(invalid(format!("{axis} needs its explicit absolute physical quantity")));
        }
        point.axes().get(axis).copied().ok_or_else(|| invalid(format!("missing {axis}")))
    };
    let t = coordinate("temperature", QuantityKind::AbsoluteTemperature)?;
    let p = coordinate("pressure", QuantityKind::Pressure)?;
    // NASA-9 regions and every pure-transport claim own the supported domain.
    // The legacy calorically-perfect air fit's envelope is not this model's.
    // Source support still does not prove ideal-gas accuracy or phase stability.
    if !positive(t) || !positive(p) {
        return Err(invalid(
            "thermochemical gas requires finite positive absolute temperature and pressure",
        ));
    }
    Ok((Temperature::new(t), Pressure::new(p)))
}

impl ThermochemicalGasModel {
    /// Bind exactly one transport card to every active canonical species.
    /// Input order is immaterial; species/card pairs are sorted together.
    /// Zero species, duplicate assignments, unknown/missing species and an
    /// absent water label refuse before any material evaluation.
    pub fn new(
        thermodynamics: FrozenIdealGasMixtureModelV1,
        mut transport: Vec<SpeciesTransportMaterial>,
        water_species: Option<&SpeciesId>,
    ) -> Result<Self, MaterialError> {
        let components = thermodynamics.components();
        if transport.len() != components.len() {
            return Err(invalid("each active thermochemical species needs exactly one transport material"));
        }
        transport.sort_by(|a, b| a.species.cmp(&b.species));
        if transport.iter().zip(components).any(|(data, model)| &data.species != model.species()) {
            return Err(invalid("duplicate, missing or foreign species in transport material assignments"));
        }
        let water_component = water_species.map(|species| {
            components.iter().position(|model| model.species() == species)
                .ok_or_else(|| invalid("explicit water species is absent from the active mixture"))
        }).transpose()?;
        Ok(Self { thermodynamics, transport, water_component, source_packs: Vec::new() })
    }

    /// Retained NASA-9 species and canonical composition.
    #[must_use]
    pub const fn thermodynamics(&self) -> &FrozenIdealGasMixtureModelV1 { &self.thermodynamics }
    /// Source transport assignments in exactly the same canonical species order.
    #[must_use]
    pub fn transport_materials(&self) -> &[SpeciesTransportMaterial] {
        &self.transport
    }
    /// Pack bindings retained by a source-backed constructor; empty for models
    /// assembled directly from cards. Empty does not assert source authority.
    #[must_use]
    pub fn source_packs(&self) -> &[packs::ThermochemicalPackReceipt] {
        &self.source_packs
    }

    /// Explicitly identified water species, if any; never inferred from a name
    /// or molar mass. A zero-fraction composition replacement omits its label.
    #[must_use]
    pub fn water_species(&self) -> Option<&SpeciesId> {
        self.water_component.map(|index| self.thermodynamics.components()[index].species())
    }

    /// Replace composition while retaining each active species' NASA-9 law,
    /// pure transport source, phase label and original pack identity. Fractions refer
    /// to [`Self::thermodynamics`]' canonical component order, not the order in
    /// which the original model was constructed. This creates a new static
    /// constitutive model; it does not change the mass or energy of a gas parcel.
    /// Explicit zero fractions are omitted before calling the shared active-
    /// species owner, so a source-backed gas can enter a field from zero without
    /// an invented trace concentration. `self` retains the full source set.
    ///
    /// # Errors
    /// The existing mixture owner refuses an unsupported basis, missing
    /// fractions, or a composition whose canonical sum is not exact one.
    pub fn with_composition(
        &self,
        composition: fs_thermochem::Composition,
    ) -> Result<Self, MaterialError> {
        let components = self.thermodynamics.components();
        let fractions = composition.fractions();
        if fractions.len() != components.len() {
            return Err(invalid("replacement composition must cover every source species in canonical order"));
        }
        let active = fractions.iter().enumerate().filter(|(_, fraction)| **fraction > 0.0)
            .map(|(index, _)| index).collect::<Vec<_>>();
        let thermodynamics = FrozenIdealGasMixtureModelV1::new(
            active.iter().map(|&index| components[index].clone()).collect(),
            fs_thermochem::Composition::new(composition.basis(),
                active.iter().map(|&index| fractions[index]).collect())
                .map_err(|error| failure("active composition", error))?,
        )
        .map_err(|error| failure("composition replacement", error))?;
        let water = self.water_component.filter(|&index| fractions[index] > 0.0)
            .map(|index| components[index].species());
        let mut revised = Self::new(
            thermodynamics,
            active.iter().map(|&index| self.transport[index].clone()).collect(),
            water,
        )?;
        revised.source_packs = self.source_packs.iter().filter(|receipt|
            revised.thermodynamics.components().iter().any(|component| component.species() == &receipt.species))
            .cloned().collect();
        Ok(revised)
    }

    /// Resolve a generated composition from nonnegative species amounts [mol].
    ///
    /// Amounts use this bank's canonical species order. Only exact zeros are
    /// omitted; a positive amount that cannot be represented as a positive
    /// mole fraction refuses instead of silently losing a reaction product.
    /// The existing generated-simplex closure corrects floating-point roundoff
    /// before the original composition owner checks its exact unit sum. Every
    /// surviving NASA-9, transport, water and pack identity is retained.
    /// This changes an immutable constitutive model, not parcel mass or energy.
    pub fn with_mole_amounts(&self, amounts: &[fs_qty::Amount]) -> Result<Self, MaterialError> {
        if amounts.len() != self.thermodynamics.components().len()
            || amounts.iter().any(|amount| !amount.value().is_finite() || amount.value() < 0.0)
        {
            return Err(invalid("mole amounts must cover every source species and be finite/nonnegative"));
        }
        // Scaling first avoids overflowing a valid extensive inventory. It also
        // makes composition independent of the caller's total parcel size.
        let scale = amounts.iter().map(|amount| amount.value()).fold(0.0_f64, f64::max);
        if !positive(scale) {
            return Err(invalid("at least one positive species amount is required"));
        }
        let mut fractions: Vec<_> = amounts.iter().map(|amount| amount.value() / scale).collect();
        let total: f64 = fractions.iter().sum();
        for fraction in &mut fractions { *fraction /= total; }
        if fractions.iter().zip(amounts).any(|(&fraction, amount)| amount.value() > 0.0 && !positive(fraction)) {
            return Err(invalid("positive species amount has an unrepresentable mole fraction"));
        }
        fs_thermal::enthalpy::advection::species::close_generated_simplex(&mut fractions)
            .map_err(|error| failure("generated mole composition", error))?;
        let composition = fs_thermochem::Composition::new(
            fs_thermochem::CompositionBasis::MoleFraction, fractions,
        ).map_err(|error| failure("generated mole composition", error))?;
        self.with_composition(composition)
    }

    /// Resolve a declared vapor in equilibrium with liquid water. The caller
    /// chooses this phase condition explicitly; the ordinary frozen-gas
    /// evaluator makes no condensation claim. Uses the existing Buck owner
    /// within its stated liquid-water temperature window, without extrapolation.
    pub fn evaluate_with_liquid_water_limit(
        &self,
        point: &QueryPoint,
    ) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        let (temperature, pressure) = coordinates(point)?;
        let water = self
            .water_component
            .ok_or_else(|| invalid("liquid-water phase check needs an explicit water species"))?;
        let fraction = self.thermodynamics.mole_composition().fractions()[water];
        let saturation = crate::gas::saturation_pressure_water_pa(temperature.value())?;
        if fraction > saturation / pressure.value() {
            return Err(invalid(
                "supersaturated thermochemical gas needs condensation physics",
            ));
        }
        self.evaluate(point)
    }

    /// Resolve ALL physical coefficients at one common typed temperature/pressure.
    /// Other source conditions remain in `point` and are never invented here.
    /// Transport cards supply `dynamic_viscosity` and `thermal_conductivity`;
    /// heat capacity, molar mass and the gas constant belong solely to NASA-9.
    /// Pure transport is evaluated at total pressure, as in the dilute Wilke law.
    ///
    /// # Errors
    /// Refuses out-of-domain/ambiguous/incomplete source claims, invalid heat
    /// capacities, invalid thermodynamic states and unrepresentable arithmetic.
    /// No partially updated medium is returned.
    pub fn evaluate(
        &self, point: &QueryPoint,
    ) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        let (temperature, pressure) = coordinates(point)?;
        let thermodynamics = self.thermodynamics.evaluate(temperature, pressure)
            .map_err(|error| failure("thermodynamics", error))?;
        let properties = thermodynamics.properties();
        let cp = properties.mass_specific_cp().value();
        let cv = properties.mass_specific_cv().value();
        let gamma = cp / cv;
        if !positive(cp) || !positive(cv) || !gamma.is_finite()
            || gamma <= 1.0 || gamma > 5.0 / 3.0 + 1.0e-9 {
            return Err(invalid("NASA-9 gas needs positive cp/cv and physical frozen gamma in (1,5/3]"));
        }
        let eos = IdealGasEosV1::new(properties.molar_mass())
            .map_err(|error| failure("EOS construction", error))?
            .evaluate_pt(temperature, pressure).map_err(|error| failure("EOS", error))?;
        let requirements = [
            ScalarPropertyRequirement::try_new("dynamic_viscosity", DynViscosity::DIMS,
                ScalarAdmissibility::StrictlyPositive).map_err(|error| failure("viscosity requirement", error))?,
            ScalarPropertyRequirement::try_new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive).map_err(|error| failure("conductivity requirement", error))?,
        ];
        let mut resolved = Vec::with_capacity(self.transport.len());
        let mut viscosities = Vec::with_capacity(self.transport.len());
        for (data, species) in self.transport.iter().zip(self.thermodynamics.components()) {
            let material = resolve_material_state_point(&data.material, point, &requirements, data.selection.clone())
                .map_err(|error| failure("transport material", error))?;
            let mu = material.property("dynamic_viscosity")
                .ok_or_else(|| invalid("missing resolved species viscosity"))?.value_si();
            viscosities.push(WilkeGasViscosityComponentV1::new(
                data.species.clone(), species.molar_mass(), DynViscosity::new(mu),
            ).map_err(|error| failure("species viscosity", error))?);
            resolved.push(material);
        }
        let viscosity = WilkeGasViscosityMixtureV1::new(
            temperature, viscosities, self.thermodynamics.mole_composition().clone(),
        ).map_err(|error| failure("viscosity mixture", error))?
            .evaluate().map_err(|error| failure("Wilke evaluation", error))?;
        let fractions = self.thermodynamics.mole_composition().fractions();
        // These immutable evaluated denominators are the L1 Wilke operation's
        // actual outputs, NOT a second interaction formula or conductivity ratio.
        let mut conductivity = 0.0;
        for ((material, &fraction), &bits) in resolved.iter().zip(fractions)
            .zip(viscosity.receipt().denominator_bits()) {
            let pure = material.property(THERMAL_CONDUCTIVITY_PROPERTY)
                .ok_or_else(|| invalid("missing resolved species conductivity"))?.value_si();
            let contribution = fraction * (pure / f64::from_bits(bits));
            if !positive(contribution) { return Err(invalid("WMS conductivity contribution is unrepresentable")); }
            conductivity += contribution;
        }
        let density = eos.mass_density().value();
        let gas_constant = eos.specific_gas_constant().value();
        let sound_speed = fs_math::det::sqrt(gamma * gas_constant * temperature.value());
        let mu = viscosity.dynamic_viscosity().value();
        let transport = FluidTransportState::try_new(density, cp, conductivity, mu)?;
        let impedance = density * sound_speed;
        let molar_mass = properties.molar_mass().value();
        let internal_energy = properties.u().value() / molar_mass;
        let enthalpy = properties.h().value() / molar_mass;
        if !positive(sound_speed) || !positive(impedance)
            || !internal_energy.is_finite() || !enthalpy.is_finite() {
            return Err(invalid("thermochemical sound speed, impedance or specific energy is unrepresentable"));
        }
        let state = GasState {
            temperature: temperature.value(), pressure: pressure.value(), density, sound_speed,
            dynamic_viscosity: mu, thermal_conductivity: conductivity, gamma,
            specific_gas_constant: gas_constant, specific_heat_cp: cp,
            prandtl: transport.prandtl(), characteristic_impedance: impedance,
            water_mole_fraction: self.water_component.map_or(0.0, |index| fractions[index]),
        };
        Ok(ResolvedThermochemicalGasState {
            state, thermodynamics, eos, viscosity, components: resolved, transport,
            specific_internal_energy_j_kg: internal_energy, specific_enthalpy_j_kg: enthalpy,
            source_packs: self.source_packs.clone(),
        })
    }
}

/// Immutable common state retaining its original thermodynamic and material sources.
/// A borrowed GasState is a numerical view, not a replacement source authority.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermochemicalGasState {
    state: GasState,
    thermodynamics: FrozenIdealGasMixtureEvaluationV1,
    eos: IdealGasEosEvaluationV1,
    viscosity: WilkeGasViscosityEvaluationV1,
    components: Vec<ResolvedMaterialStatePoint>,
    transport: FluidTransportState,
    specific_internal_energy_j_kg: f64,
    specific_enthalpy_j_kg: f64,
    source_packs: Vec<packs::ThermochemicalPackReceipt>,
}
impl ResolvedThermochemicalGasState {
    /// Exact species/model/transport pack identities used to construct this
    /// state. Composition and evaluation receipts remain in `thermodynamics`.
    #[must_use]
    pub fn source_packs(&self) -> &[packs::ThermochemicalPackReceipt] {
        &self.source_packs
    }
    /// Complete medium consumed by the existing duct, radiation and gas solvers.
    #[must_use]
    pub const fn state(&self) -> &GasState { &self.state }
    /// Original species thermodynamics, including source enthalpy references.
    #[must_use]
    pub const fn thermodynamics(&self) -> &FrozenIdealGasMixtureEvaluationV1 { &self.thermodynamics }
    /// Mechanical EOS evaluated with that mixture's actual molar mass.
    #[must_use]
    pub const fn equation_of_state(&self) -> &IdealGasEosEvaluationV1 { &self.eos }
    /// Actual shared L1 viscosity evaluation.
    #[must_use]
    pub const fn viscosity(&self) -> &WilkeGasViscosityEvaluationV1 { &self.viscosity }
    /// Pure transport source bundles in canonical species order.
    #[must_use]
    pub fn components(&self) -> &[ResolvedMaterialStatePoint] { &self.components }
    /// Mass-specific internal energy [J/kg], including the declared reference.
    #[must_use]
    pub const fn specific_internal_energy_j_kg(&self) -> f64 { self.specific_internal_energy_j_kg }
    /// Mass-specific enthalpy [J/kg], including the declared reference.
    #[must_use]
    pub const fn specific_enthalpy_j_kg(&self) -> f64 { self.specific_enthalpy_j_kg }
}
impl FluidTransport for ResolvedThermochemicalGasState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> { Ok(self.transport) }
}

#[cfg(test)]
#[path = "thermochemical_gas_tests.rs"]
pub(crate) mod tests;
