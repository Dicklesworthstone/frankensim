//! Frozen-composition ambient media for every consumer of [`GasState`].
//!
//! Composition uses the shared `fs-qty` mass/mole-fraction authority. At one
//! common (T, p), heat capacities and molar mass use ideal-mixture sums;
//! viscosity uses Wilke and conductivity uses Wassiljewa-Mason-Saxena with
//! Wilke's interaction function. These are the same constitutive approximations
//! as the binary moist-air path, generalized to a bounded number of components.
//!
//! This is a low-pressure, nonreacting gas model, not chemical equilibrium,
//! combustion, condensation, real-gas thermodynamics, or an experimental error
//! bound. Component states must be pure-gas models or explicitly declared
//! pseudocomponents. Wilke mixing is not associative: do not silently replace
//! resolved species by previously mixed pseudocomponents. No species identity,
//! water content, transport law, or missing property is inferred from a name.
//! ISO 9613 is an air-specific law and must not be applied to arbitrary mixtures.
//!
//! References: Wilke, J. Chem. Phys. 18, 517 (1950), doi:10.1063/1.1747673;
//! Mason and Saxena, Phys. Fluids 1, 361 (1958), doi:10.1063/1.1724352.

use fs_qty::semantic::Composition;
use fs_qty::MolarMass;

use crate::gas::{GasSpec, GasState, R_USSA_1976, ResolvedGasState};
use crate::MaterialError;

/// Cold-path resource bound: mixing is quadratic in the number of components.
pub const MAX_GAS_MIXTURE_COMPONENTS: usize = 128;

fn invalid(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

fn check_count(count: usize, composition: &Composition) -> Result<(), MaterialError> {
    if count == 0 || count > MAX_GAS_MIXTURE_COMPONENTS {
        return Err(invalid("gas mixture needs 1..=128 components"));
    }
    if composition.fractions().len() != count {
        return Err(invalid("gas components and composition must have identical lengths"));
    }
    if composition.fractions().iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(invalid("omit zero-fraction gas components; all active fractions must be positive"));
    }
    Ok(())
}

// Coherence belongs to GasState. The additional pure-component requirement
// belongs only to this mixing law; it must not restrict downstream media.
fn validate_component(state: &GasState) -> Result<(), MaterialError> {
    state.validate_coherence()?;
    if state.water_mole_fraction != 0.0 {
        return Err(invalid("gas component is already humid; supply explicitly labeled pure components"));
    }
    Ok(())
}

/// Wilke's positive interaction coefficient. The same coefficient belongs to
/// the WMS conductivity model; it is not evaluated from conductivity ratios.
fn interaction(left: &GasState, right: &GasState, left_mass: f64, right_mass: f64) -> f64 {
    let numerator = 1.0 + (left.dynamic_viscosity / right.dynamic_viscosity).sqrt()
        * (right_mass / left_mass).sqrt().sqrt();
    numerator * numerator / (8.0 * (1.0 + left_mass / right_mass)).sqrt()
}

fn evaluate(
    components: &[GasState],
    composition: &Composition,
    water_component: Option<usize>,
) -> Result<(GasState, Composition), MaterialError> {
    check_count(components.len(), composition)?;
    if water_component.is_some_and(|index| index >= components.len()) {
        return Err(invalid("water component index is outside the mixture"));
    }
    let reference = components[0];
    let mut masses = Vec::with_capacity(components.len());
    for component in components {
        validate_component(component)?;
        if component.temperature != reference.temperature || component.pressure != reference.pressure {
            return Err(invalid("gas components must share the exact temperature and total pressure"));
        }
        let mass = R_USSA_1976 / component.specific_gas_constant;
        if !mass.is_finite() || mass <= 0.0 {
            return Err(invalid("gas component molar mass is unrepresentable"));
        }
        masses.push(MolarMass::new(mass));
    }
    let mole = composition.to_mole_fractions(&masses)
        .map_err(|error| invalid(format!("gas mixture composition: {error}")))?;
    let fractions = mole.fractions();
    let sum: f64 = fractions.iter().sum();
    if fractions.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || !sum.is_finite()
        || (sum - 1.0).abs() > 32.0 * f64::EPSILON * components.len() as f64
    {
        return Err(invalid("gas mole fractions must remain positive and sum to one within roundoff"));
    }
    // Preserve the existing pure-gas evaluator bit for bit.
    if components.len() == 1 {
        let mut state = reference;
        state.water_mole_fraction = if water_component.is_some() { 1.0 } else { 0.0 };
        return Ok((state, mole));
    }
    let mut molar_mass = 0.0;
    let mut cp_molar = 0.0;
    let mut viscosity = 0.0;
    let mut conductivity = 0.0;
    for (i, component) in components.iter().enumerate() {
        let xi = fractions[i];
        let mi = masses[i].value();
        molar_mass += xi * mi;
        cp_molar += xi * component.specific_heat_cp * mi;
        let mut denominator = 0.0;
        for (j, other) in components.iter().enumerate() {
            let phi = if i == j { 1.0 } else {
                interaction(component, other, mi, masses[j].value())
            };
            if !phi.is_finite() || phi <= 0.0 {
                return Err(invalid("gas transport interaction is unrepresentable"));
            }
            denominator += fractions[j] * phi;
        }
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(invalid("gas transport mixing denominator is unrepresentable"));
        }
        viscosity += xi * component.dynamic_viscosity / denominator;
        conductivity += xi * component.thermal_conductivity / denominator;
    }
    let cv_molar = cp_molar - R_USSA_1976;
    if !molar_mass.is_finite() || molar_mass <= 0.0 || !cv_molar.is_finite() || cv_molar <= 0.0 {
        return Err(invalid("gas mixture heat capacity or molar mass is nonphysical"));
    }
    let r = R_USSA_1976 / molar_mass;
    let cp = cp_molar / molar_mass;
    let gamma = cp_molar / cv_molar;
    let density = reference.pressure / (r * reference.temperature);
    let sound_speed = (gamma * r * reference.temperature).sqrt();
    let mut state = GasState {
        temperature: reference.temperature,
        pressure: reference.pressure,
        density,
        sound_speed,
        dynamic_viscosity: viscosity,
        thermal_conductivity: conductivity,
        gamma,
        specific_gas_constant: r,
        specific_heat_cp: cp,
        prandtl: viscosity * cp / conductivity,
        characteristic_impedance: density * sound_speed,
        water_mole_fraction: 0.0,
    };
    validate_component(&state)?;
    state.water_mole_fraction = water_component.map_or(0.0, |index| fractions[index]);
    Ok((state, mole))
}

/// Mix component states into the existing shared medium, not an audio-only model.
///
/// Components are evaluated at the same total mixture pressure, not their
/// partial pressures: the component transport laws are dilute-gas models.
/// The optional index explicitly labels water; it is not inferred from molar
/// mass. Water fraction alone does not establish saturation or ISO-air validity.
/// Input order is retained; repeated identical ordered evaluations are deterministic.
/// Evaluation allocates on the cold path and is bounded by 128 squared terms.
///
/// # Errors
/// Refuses invalid/incoherent components, unequal states, absent/zero components,
/// inadmissible mass/mole composition, and unrepresentable derived properties.
pub fn mix_gas_states(
    components: &[GasState],
    composition: &Composition,
    water_component: Option<usize>,
) -> Result<GasState, MaterialError> {
    evaluate(components, composition, water_component).map(|(state, _)| state)
}

impl GasState {
    /// Evaluate arbitrary declared gas specs and mix them through the shared
    /// ideal-gas/Wilke/WMS path. Species changes require data, not solver edits.
    ///
    /// # Errors
    /// Propagates component gas-model and [`mix_gas_states`] refusals.
    pub fn try_new_mixture(
        specs: &[GasSpec],
        composition: &Composition,
        temperature: f64,
        pressure: f64,
        water_component: Option<usize>,
    ) -> Result<Self, MaterialError> {
        check_count(specs.len(), composition)?;
        let components = specs.iter().map(|spec| Self::try_new(spec, temperature, pressure))
            .collect::<Result<Vec<_>, _>>()?;
        mix_gas_states(&components, composition, water_component)
    }
}

/// A source-backed mixture retaining every original property receipt and the
/// declared composition. Component bundle identities do not by themselves
/// identify this mixture: composition, ordering and water labeling are retained too.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedGasMixtureState {
    state: GasState,
    components: Vec<ResolvedGasState>,
    declared_composition: Composition,
    mole_composition: Composition,
    water_component: Option<usize>,
}

impl ResolvedGasMixtureState {
    /// Join existing `resolve_sutherland_gas_state` results without replacing
    /// any source-selected property with a representative gas preset.
    ///
    /// # Errors
    /// Propagates [`mix_gas_states`] refusals; publishes no partial mixture.
    pub fn try_new(
        components: Vec<ResolvedGasState>,
        composition: Composition,
        water_component: Option<usize>,
    ) -> Result<Self, MaterialError> {
        check_count(components.len(), &composition)?;
        let states: Vec<_> = components.iter().map(|component| *component.state()).collect();
        let (state, mole_composition) = evaluate(&states, &composition, water_component)?;
        Ok(Self { state, components, declared_composition: composition, mole_composition, water_component })
    }

    /// The same medium accepted by duct, radiation, flow and heat consumers.
    #[must_use]
    pub const fn state(&self) -> &GasState { &self.state }

    /// Source-resolved components, in the declared composition order.
    #[must_use]
    pub fn components(&self) -> &[ResolvedGasState] { &self.components }

    /// Caller-declared mass or mole fractions, without a second basis system.
    #[must_use]
    pub const fn declared_composition(&self) -> &Composition { &self.declared_composition }

    /// Mole fractions actually used by the thermodynamic and transport sums.
    #[must_use]
    pub const fn mole_composition(&self) -> &Composition { &self.mole_composition }

    /// Explicit water label; not a chemical-identity inference.
    #[must_use]
    pub const fn water_component(&self) -> Option<usize> { self.water_component }
}
