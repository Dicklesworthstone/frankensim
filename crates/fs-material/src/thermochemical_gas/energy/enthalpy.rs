//! Isobaric NASA-9 enthalpy closure and immutable caloric laws for transport.
//!
//! Enthalpy includes each species' declared formation-energy reference. Mixing
//! transports those joules; it never averages temperature or replaces cp(T) by
//! a guessed constant. Thermal relaxation remains owned by fs-thermal.

use super::{GasEnergyError, GasEnergySolve, map_root, poll, root, updated_point};
use crate::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel, coordinates};
use fs_matdb::QueryPoint;
use fs_qty::Temperature;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use fs_thermal::enthalpy::advection::species::{CaloricComposition, close_generated_simplex};
use fs_thermochem::{Composition, CompositionBasis, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

/// Accepted isobaric state with an absolute enthalpy residual.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEnthalpySolution {
    state: ResolvedThermochemicalGasState,
    residual_j_kg: f64,
    iterations: usize,
}
impl GasEnthalpySolution {
    /// Common thermodynamic, transport and acoustic state.
    pub fn state(&self) -> &ResolvedThermochemicalGasState { &self.state }
    /// Actual h(T,p) minus requested h [J/kg].
    pub fn residual_j_kg(&self) -> f64 { self.residual_j_kg }
    /// Interior iterations of the shared safeguarded root solver.
    pub fn iterations(&self) -> usize { self.iterations }
}

fn temperature(
    model: &ThermochemicalGasModel, context: &QueryPoint, enthalpy_j_kg: f64,
    controls: GasEnergySolve, cancelled: &mut impl FnMut() -> bool,
) -> Result<root::Solution, GasEnergyError> {
    controls.validate()?;
    let (_, pressure) = coordinates(context)?;
    if !enthalpy_j_kg.is_finite() {
        return Err(GasEnergyError::Invalid("finite specific enthalpy required"));
    }
    let thermal = |t| -> Result<root::Sample, GasEnergyError> {
        let evaluated = model.thermodynamics().evaluate(Temperature::new(t), pressure)
            .map_err(|error| crate::thermochemical_gas::failure("enthalpy thermodynamics", error))?;
        let properties = evaluated.properties();
        Ok(root::Sample {
            residual: properties.h().value() / properties.molar_mass().value() - enthalpy_j_kg,
            tangent: properties.mass_specific_cp().value(),
            pressure_pa: pressure.value(),
        })
    };
    root::solve(controls.minimum_temperature_k, controls.maximum_temperature_k,
        controls.energy_tolerance_j_kg, controls.maximum_iterations, thermal, cancelled)
        .map_err(map_root)
}

impl ThermochemicalGasModel {
    /// Reconstruct temperature at the explicitly typed pressure in context.
    /// All other source coordinates survive; the accepted medium is resolved
    /// once through the existing complete material/thermochemistry owner.
    /// The caller's bracket must describe one continuous increasing gas branch.
    pub fn reconstruct_pressure_enthalpy(
        &self, context: &QueryPoint, enthalpy_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnthalpySolution, GasEnergyError> {
        let solution = temperature(self, context, enthalpy_j_kg, controls, &mut cancelled)?;
        poll(&mut cancelled)?;
        let state = self.evaluate(&updated_point(context, solution.temperature_k,
            solution.sample.pressure_pa)?)?;
        let residual = state.specific_enthalpy_j_kg() - enthalpy_j_kg;
        if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j_kg {
            return Err(GasEnergyError::NotConverged {
                iterations: solution.iterations, residual_j_kg: residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(GasEnthalpySolution { state, residual_j_kg: residual, iterations: solution.iterations })
    }
}

/// Source-bound immutable h(T)/T(h) law for one stationary mixture.
/// Every cell law retains the FULL common species order in its descriptor,
/// even when its active thermochemical model omits zero-concentration species.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEnthalpyLaw {
    model: ThermochemicalGasModel,
    context: QueryPoint,
    controls: GasEnergySolve,
    species: Vec<String>,
    mass_fractions: Vec<f64>,
    mole_composition: Composition,
    enthalpy_reference: String,
    maximum_slope: f64,
}

fn conversion_error(error: impl core::fmt::Display) -> GasEnergyError {
    crate::thermochemical_gas::failure("transport composition", error).into()
}

/// Convert a complete source-order composition with the EXISTING typed basis
/// conversion. Only its generated floating-point simplex residual is closed.
pub(crate) fn converted(
    source: &ThermochemicalGasModel, composition: &Composition, basis: CompositionBasis,
) -> Result<Composition, GasEnergyError> {
    if composition.fractions().len() != source.thermodynamics().components().len()
        || composition.fractions().iter().sum::<f64>() != 1.0
    {
        return Err(GasEnergyError::Invalid("composition must cover the full canonical source with exact unit sum"));
    }
    let masses: Vec<_> = source.thermodynamics().components().iter()
        .map(|species| species.molar_mass()).collect();
    let result = match basis {
        CompositionBasis::MoleFraction => composition.to_mole_fractions(&masses),
        CompositionBasis::MassFraction => composition.to_mass_fractions(&masses),
        CompositionBasis::VolumeFraction => {
            return Err(GasEnergyError::Invalid("volume-fraction transport is not admitted"));
        }
    }.map_err(conversion_error)?;
    let mut fractions = result.fractions().to_vec();
    if composition.fractions().iter().zip(&fractions).any(|(&old, &new)| old > 0.0 && new == 0.0) {
        return Err(GasEnergyError::Invalid("composition basis conversion lost an active species"));
    }
    if composition.basis() != basis {
        close_generated_simplex(&mut fractions).map_err(conversion_error)?;
    }
    Composition::new(basis, fractions).map_err(conversion_error)
}

// Positive interval powers and outward-rounded term/sum bounds. This bounds
// source coefficients; actual cp and h are ALWAYS evaluated by NASA-9 itself.
#[derive(Clone, Copy)]
struct Interval { lo: f64, hi: f64 }
impl Interval {
    fn multiply_positive(self, rhs: Self) -> Self {
        Self { lo: (self.lo * rhs.lo).next_down(), hi: (self.hi * rhs.hi).next_up() }
    }
    fn scaled_lower(self, scale: f64) -> f64 {
        (scale * if scale >= 0.0 { self.lo } else { self.hi }).next_down()
    }
}

pub(super) fn minimum_cp(
    species: &Nasa9StandardStateModelV1, controls: GasEnergySolve,
) -> Result<f64, GasEnergyError> {
    let lo = controls.minimum_temperature_k;
    let hi = controls.maximum_temperature_k;
    // A single source region avoids asserting continuity across independently
    // fitted NASA-9 charts. Breath-temperature brackets fit within this scope.
    let regions = species.regions();
    let region = regions.iter().enumerate().find_map(|(index, region)| {
        let upper = region.temperature_max().value();
        if lo >= region.temperature_min().value()
            && (hi < upper || (index + 1 == regions.len() && hi <= upper))
        { Some(region) } else { None }
    }).ok_or(GasEnergyError::Invalid("caloric bracket must stay within one NASA-9 region per active species"))?;
    let t = Interval { lo, hi };
    let inverse = Interval { lo: (1.0 / hi).next_down(), hi: (1.0 / lo).next_up() };
    let square = t.multiply_positive(t);
    let cube = square.multiply_positive(t);
    let fourth = cube.multiply_positive(t);
    let basis = [inverse.multiply_positive(inverse), inverse,
        Interval { lo: 1.0, hi: 1.0 }, t, square, cube, fourth];
    let mut lower = 0.0;
    for (&coefficient, term) in region.coefficients()[..7].iter().zip(basis) {
        if coefficient != 0.0 {
            lower = (lower + term.scaled_lower(coefficient)).next_down();
        }
    }
    let cp = ((lower * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K).next_down()
        / species.molar_mass().value()).next_down();
    if !cp.is_finite() || cp <= 0.0 {
        return Err(GasEnergyError::Invalid("source coefficient bound cannot establish positive cp; narrow the caloric bracket"));
    }
    Ok(cp)
}

impl GasEnthalpyLaw {
    /// Construct a frozen mixture caloric law from one verified common source
    /// set. Its mass fractions are explicit, canonical and immutable.
    /// A conservative source-coefficient cp lower bound supplies the existing
    /// thermal CFL. No caller-invented heat capacity or temperature chart is used.
    pub fn new(
        source: &ThermochemicalGasModel, mass_composition: Composition,
        context: QueryPoint, controls: GasEnergySolve,
    ) -> Result<Self, GasEnergyError> {
        if mass_composition.basis() != CompositionBasis::MassFraction {
            return Err(GasEnergyError::Invalid("stationary caloric laws require mass fractions"));
        }
        let mole_composition = converted(source, &mass_composition, CompositionBasis::MoleFraction)?;
        Self::from_compositions(source, mass_composition, mole_composition, context, controls)
    }

    /// Retain authored inlet mole fractions exactly while deriving the mass
    /// fractions used by the common upwind transport owner.
    pub fn from_mole(
        source: &ThermochemicalGasModel, mole_composition: Composition,
        context: QueryPoint, controls: GasEnergySolve,
    ) -> Result<Self, GasEnergyError> {
        if mole_composition.basis() != CompositionBasis::MoleFraction {
            return Err(GasEnergyError::Invalid("inlet caloric laws require explicit mole fractions"));
        }
        let mass_composition = converted(source, &mole_composition, CompositionBasis::MassFraction)?;
        Self::from_compositions(source, mass_composition, mole_composition, context, controls)
    }

    fn from_compositions(
        source: &ThermochemicalGasModel, mass_composition: Composition,
        mole_composition: Composition, context: QueryPoint, controls: GasEnergySolve,
    ) -> Result<Self, GasEnergyError> {
        controls.validate()?;
        coordinates(&context)?;
        let model = source.with_composition(mole_composition.clone())?;
        let convention = source.thermodynamics().components()[0].convention();
        if convention.phase() != StandardStatePhaseV1::Gas
            || convention.eos() != ReferenceEquationOfStateV1::IdealGas
            || source.thermodynamics().components().iter().any(|component| component.convention() != convention)
        {
            return Err(GasEnergyError::Invalid("stationary gas enthalpy requires one common ideal-gas reference convention"));
        }
        // Every admitted mixture is a positive mass-weighted combination, so
        // the smallest active pure-species lower bound bounds the mixture.
        let mut lower: f64 = f64::INFINITY;
        for species in model.thermodynamics().components() {
            lower = lower.min(minimum_cp(species, controls)?);
        }
        let maximum_slope = (1.0 / lower).next_up();
        if !maximum_slope.is_finite() || maximum_slope <= 0.0 {
            return Err(GasEnergyError::Invalid("caloric CFL slope is not representable"));
        }
        let enthalpy_reference = format!("NASA9:{:?}:{:?}:{:016x}:{}",
            convention.phase(), convention.eos(), convention.reference_pressure().value().to_bits(),
            convention.elemental_reference().as_str());
        Ok(Self {
            model, context, controls,
            species: source.thermodynamics().components().iter()
                .map(|component| component.species().to_string()).collect(),
            mass_fractions: mass_composition.fractions().to_vec(),
            mole_composition, enthalpy_reference, maximum_slope,
        })
    }

    /// The active mixture still retains its original thermochemical/material receipts.
    pub fn model(&self) -> &ThermochemicalGasModel { &self.model }
    /// Full canonical mole composition for the shared axial-field owner.
    pub fn mole_composition(&self) -> &Composition { &self.mole_composition }
    /// Unchanged source conditions, including common pressure.
    pub fn context(&self) -> &QueryPoint { &self.context }

    /// Source NASA-9 specific enthalpy at an admitted temperature [J/kg].
    pub fn enthalpy_at_temperature(&self, temperature_k: f64) -> Result<f64, GasEnergyError> {
        if !temperature_k.is_finite() || temperature_k < self.controls.minimum_temperature_k
            || temperature_k > self.controls.maximum_temperature_k
        {
            return Err(GasEnergyError::Invalid("temperature outside the declared caloric bracket"));
        }
        let (_, pressure) = coordinates(&self.context)?;
        let state = self.model.thermodynamics().evaluate(Temperature::new(temperature_k), pressure)
            .map_err(conversion_error)?;
        let properties = state.properties();
        let h = properties.h().value() / properties.molar_mass().value();
        if !h.is_finite() { return Err(GasEnergyError::Invalid("source enthalpy is not finite")); }
        Ok(h)
    }

    /// Full common state at the specified h and retained p.
    pub fn resolve_enthalpy(
        &self, enthalpy_j_kg: f64, cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnthalpySolution, GasEnergyError> {
        self.model.reconstruct_pressure_enthalpy(&self.context, enthalpy_j_kg, self.controls, cancelled)
    }
}

impl ThermalLaw for GasEnthalpyLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        temperature(&self.model, &self.context, h, self.controls, &mut || false)
            .map(|solution| solution.temperature_k).map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 { self.maximum_slope }
    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        Some(CaloricComposition { species: &self.species, mass_fractions: &self.mass_fractions,
            enthalpy_reference: &self.enthalpy_reference })
    }
}

#[cfg(test)]
mod tests;
