//! Reacting ideal-gas enthalpy on the original finite thermal-storage owner.
//!
//! Thermal cells store source-relative h at fixed pressure and fixed elemental
//! ratios. The original Gibbs/HP law supplies T(h), changing composition,
//! reacting Cp, density and chemical entropy. Existing fs-thermal timelines and
//! links own finite heat transfer; this module supplies no second integrator.
//!
//! Only pV work is included in the enthalpy coordinate. Reference mesh geometry
//! does not become a rigid-vessel constraint: actual volume follows conserved
//! mass and the returned equilibrium density when a geometric consumer asks.
//! This is the infinitely fast, neutral ideal-gas equilibrium limit, with no
//! kinetics, condensed phases, plasma transport or finite pressure controller.

/// Same reacting calorics and an explicitly declared bulk conductivity curve.
pub mod conduction;

use super::{
    GasEquilibriumEnthalpyError as Error, GasEquilibriumEnthalpySolve,
    GasEquilibriumEnthalpyState, ThermochemicalGasModel, freeze, poll,
};
use crate::thermochemical_gas::{coordinates, positive};
use crate::thermochemical_gas::energy::{GasEnergySolve, enthalpy::minimum_cp};
use fs_matdb::QueryPoint;
use fs_math::dd::Dd;
use fs_qty::{Amount, Pressure, Temperature};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    IdealGasEquilibriumV1,
    enthalpy::{
        EquilibriumEnthalpyPressureSampleV1,
        isobaric::EquilibriumIsobaricCaloricV1,
    },
};
use std::cell::RefCell;

/// Immutable reacting caloric data with reusable private numerical scratch.
///
/// Storage is zero at the cold equilibrium, not at an arbitrary formation
/// energy datum. Any positive cell mass may use the same elemental ratios;
/// the retained feed mass only normalizes extensive chemistry into J/kg.
/// There is deliberately no frozen-species advection descriptor. Different
/// caloric anchors are not automatically compatible stream enthalpy data.
pub struct GasEquilibriumThermalStorage<'a> {
    bank: &'a ThermochemicalGasModel,
    context: QueryPoint,
    controls: GasEquilibriumEnthalpySolve,
    reference: GasEquilibriumEnthalpyState,
    caloric: RefCell<EquilibriumIsobaricCaloricV1<'a>>,
    mass_kg: f64,
    pressure: Pressure,
    temperature_domain: [Temperature; 2],
    enthalpy_domain_j_kg: [f64; 2],
    maximum_slope: f64,
}

impl<'a> GasEquilibriumThermalStorage<'a> {
    /// Prepare original equilibrium calorics and material support at both ends.
    ///
    /// Pressure and auxiliary material axes come from the typed context; the
    /// declared lower temperature defines the storage reference. Every eligible
    /// NASA source must support the complete continuous temperature interval.
    /// Original material/transport resolution remains a point-state operation:
    /// endpoint admission is not a promise about unsupported interior transport
    /// claims. Caloric trials use NASA data only; resolve accepted products
    /// before consuming their frozen transport or acoustic properties.
    pub fn new(
        bank: &'a ThermochemicalGasModel,
        chemistry: &'a IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        controls: GasEquilibriumEnthalpySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !positive(controls.enthalpy_tolerance_j) {
            return Err(Error::Invalid("positive finite thermal material enthalpy allowance required"));
        }
        bank.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)?;
        let (_, pressure) = coordinates(context)
            .map_err(super::super::GasEquilibriumChamberError::from)?;
        let mut caloric = EquilibriumIsobaricCaloricV1::new(
            chemistry, pressure, feed, controls.equilibrium, &mut cancelled,
        )?;
        let mass_kg = caloric.feed_mass_kg();
        let temperature_domain = caloric.temperature_domain();
        let enthalpy_domain_j_kg = caloric.specific_enthalpy_domain_j_kg();
        let cold_sample = caloric.sample_specific_enthalpy(0.0, &mut cancelled)?;
        let reference = freeze(
            bank, context, caloric.reference_state().clone(), cold_sample,
            mass_kg, None, controls, &mut cancelled,
        )?;

        // Reuse the original interval-coefficient Cp bound, not endpoint
        // samples. At stable ideal-gas equilibrium, minimizing Gibbs over
        // composition adds a nonnegative redistribution contribution to Cp:
        // Cp_eq >= sum(n_i cp_i). The pure mass-specific minimum therefore
        // bounds Cp_eq/product_mass. Include the original admitted discrepancy
        // between product mass and the conserved normalization feed mass.
        let source_controls = GasEnergySolve {
            minimum_temperature_k: temperature_domain[0].value(),
            maximum_temperature_k: temperature_domain[1].value(),
            energy_tolerance_j_kg: controls.equilibrium.enthalpy_tolerance_j / mass_kg,
            maximum_iterations: controls.equilibrium.maximum_temperature_iterations.min(128),
        };
        source_controls.validate().map_err(super::super::GasEquilibriumChamberError::from)?;
        let mut minimum: f64 = f64::INFINITY;
        for (source, amount) in chemistry.models().iter()
            .zip(reference.equilibrium().equilibrium().species()) {
            poll(&mut cancelled)?;
            if amount.log_mole_fraction.is_some() {
                minimum = minimum.min(minimum_cp(source, source_controls)
                    .map_err(super::super::GasEquilibriumChamberError::from)?);
            }
        }
        let retained_mass = (1.0 - controls.equilibrium.equilibrium.relative_balance_tolerance).next_down();
        let lower_capacity = (minimum * retained_mass).next_down();
        let maximum_slope = (1.0 / lower_capacity).next_up();
        if !positive(lower_capacity) || !positive(maximum_slope) {
            return Err(Error::Invalid("no representable source bound for reacting thermal storage"));
        }
        bound_sample(cold_sample, mass_kg, maximum_slope)?;
        let hot_sample = caloric.sample_temperature(temperature_domain[1], &mut cancelled)?;
        bound_sample(hot_sample, mass_kg, maximum_slope)?;
        let hot = caloric.publish_current(&mut cancelled)?;
        freeze(
            bank, context, hot, hot_sample, mass_kg, Some(reference.parcel()),
            controls, &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            bank, context: context.clone(), controls, reference,
            caloric: RefCell::new(caloric), mass_kg, pressure,
            temperature_domain, enthalpy_domain_j_kg, maximum_slope,
        })
    }

    /// Original cold material endpoint, source receipts and elemental reference.
    pub const fn reference_state(&self) -> &GasEquilibriumEnthalpyState { &self.reference }
    /// Source pressure held throughout this thermal enthalpy chart.
    pub const fn pressure(&self) -> Pressure { self.pressure }
    /// Original normalization mass [kg], not a restriction on thermal-cell mass.
    pub const fn feed_mass_kg(&self) -> f64 { self.mass_kg }
    /// Complete continuous original NASA temperature interval [K].
    pub const fn temperature_domain(&self) -> [Temperature; 2] { self.temperature_domain }
    /// Inclusive source-relative h interval [J/kg], with zero at the cold state.
    pub const fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] { self.enthalpy_domain_j_kg }

    /// One original HP inverse supplies T, density and reacting derivatives.
    /// Private scratch is reset to the identical anchor before every query;
    /// previous failures or query order cannot select a different warm start.
    pub fn caloric_state_at_specific_enthalpy(
        &self, h_j_kg: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnthalpyPressureSampleV1, Error> {
        poll(&mut cancelled)?;
        let mut caloric = self.caloric.try_borrow_mut()
            .map_err(|_| Error::Invalid("reacting caloric query is already in progress"))?;
        let sample = caloric.sample_specific_enthalpy(h_j_kg, &mut cancelled)?;
        bound_sample(sample, self.mass_kg, self.maximum_slope)?;
        poll(&mut cancelled)?;
        Ok(sample)
    }

    /// Source-integrated specific enthalpy and reacting Cp at one temperature.
    /// Returns (h [J/kg], Cp/feed_mass [J/(kg K)]), with formation terms once.
    pub fn enthalpy_and_heat_capacity_at_temperature(
        &self, temperature: Temperature, mut cancelled: impl FnMut() -> bool,
    ) -> Result<(f64, f64), Error> {
        poll(&mut cancelled)?;
        let mut caloric = self.caloric.try_borrow_mut()
            .map_err(|_| Error::Invalid("reacting caloric query is already in progress"))?;
        let sample = caloric.sample_temperature(temperature, &mut cancelled)?;
        bound_sample(sample, self.mass_kg, self.maximum_slope)?;
        let h = sample.enthalpy_change_j / self.mass_kg;
        let capacity = sample.response.heat_capacity_p_j_k / self.mass_kg;
        if !h.is_finite() || !positive(capacity)
            || !(self.enthalpy_domain_j_kg[0]..=self.enthalpy_domain_j_kg[1]).contains(&h) {
            return Err(Error::Invalid("reacting caloric forward state is not representable"));
        }
        poll(&mut cancelled)?;
        Ok((h, capacity))
    }

    /// Resolve a complete accepted material endpoint at this same stored h.
    ///
    /// Original source/transport conditions and the independent material H and
    /// mass balances are checked here, including source holes inside otherwise
    /// valid caloric support. The returned parcel's transport and acoustics are
    /// explicitly frozen; its separate native response is reacting.
    pub fn resolve_specific_enthalpy(
        &self, h_j_kg: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEquilibriumEnthalpyState, Error> {
        poll(&mut cancelled)?;
        let (sample, native) = {
            let mut caloric = self.caloric.try_borrow_mut()
                .map_err(|_| Error::Invalid("reacting caloric query is already in progress"))?;
            let sample = caloric.sample_specific_enthalpy(h_j_kg, &mut cancelled)?;
            bound_sample(sample, self.mass_kg, self.maximum_slope)?;
            let native = caloric.publish_current(&mut cancelled)?;
            (sample, native)
        };
        let state = freeze(
            self.bank, &self.context, native, sample, self.mass_kg,
            Some(self.reference.parcel()), self.controls, &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    /// Original reacting inverse-temperature mean times the actual stored dh.
    ///
    /// Direct local source entropy retains changing composition and mixing.
    /// Scaling that source mean by the stored increment avoids losing heat in
    /// two bounded inversions; the coincident-temperature limit is dh/T. Its
    /// accuracy still depends on the original absolute HP solve tolerance.
    pub fn thermal_entropy_change_j_kg_k(
        &self, initial_h: f64, final_h: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(&mut cancelled)?;
        let change = self.caloric.try_borrow_mut()
            .map_err(|_| Error::Invalid("reacting caloric query is already in progress"))?
            .change_between_specific_enthalpies(initial_h, final_h, &mut cancelled)?;
        bound_sample(change.initial, self.mass_kg, self.maximum_slope)?;
        bound_sample(change.final_state, self.mass_kg, self.maximum_slope)?;
        if initial_h == final_h {
            poll(&mut cancelled)?;
            return Ok(0.0);
        }
        let dh = final_h - initial_h;
        // A local source difference must agree with the requested stored heat
        // plus the two already reported endpoint inverse residuals. Do not
        // conceal chemical-inventory or formation-datum drift by rescaling an
        // unrelated positive entropy/enthalpy ratio.
        let local_residual = (Dd::from_f64(change.enthalpy_change_j)
            - Dd::from_f64(change.requested_enthalpy_change_j)
            - Dd::from_f64(change.final_state.enthalpy_residual_j)
            + Dd::from_f64(change.initial.enthalpy_residual_j)).to_f64();
        if !local_residual.is_finite()
            || local_residual.abs() > self.controls.equilibrium.enthalpy_tolerance_j {
            return Err(Error::Invalid("reacting local source heat does not close its stored interval"));
        }
        let (t0, t1) = (change.initial.temperature.value(), change.final_state.temperature.value());
        if (dh > 0.0 && t1 < t0) || (dh < 0.0 && t1 > t0) {
            return Err(Error::Invalid("reacting thermal inverse is not monotone at the stored interval"));
        }
        let inverse_temperature = if t0 == t1 {
            1.0 / t0
        } else {
            if change.enthalpy_change_j == 0.0
                || change.enthalpy_change_j.is_sign_positive() != dh.is_sign_positive()
                || change.entropy_change_j_k == 0.0
                || change.entropy_change_j_k.is_sign_positive() != dh.is_sign_positive() {
                return Err(Error::Invalid("reacting source entropy or enthalpy interval is unresolved"));
            }
            let mean = change.entropy_change_j_k / change.enthalpy_change_j;
            // Along this positive-Cp fixed-pressure branch, integral dH/T
            // divided by integral dH lies between the endpoint inverse
            // temperatures. Endpoint HP tolerances govern the heat gate above,
            // not this interval at the actual accepted temperatures. Allow only
            // explicit reciprocal/ratio roundoff; refuse rather than clamp an
            // unresolved chemical source mean.
            let lower = (1.0 / t0.max(t1)).next_down();
            let upper = (1.0 / t0.min(t1)).next_up();
            let allowance = (64.0 * f64::EPSILON * upper).next_up();
            let violation = if mean < lower { lower - mean }
                else if mean > upper { mean - upper } else { 0.0 };
            if !positive(mean) || !positive(lower) || !positive(upper)
                || !allowance.is_finite() || violation > allowance {
                return Err(Error::Invalid("reacting source entropy mean is outside its thermal interval"));
            }
            mean
        };
        let entropy = dh * inverse_temperature;
        if !positive(inverse_temperature) || !entropy.is_finite() || entropy == 0.0 {
            return Err(Error::Invalid("reacting thermal entropy increment is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(entropy)
    }
}

fn bound_sample(
    sample: EquilibriumEnthalpyPressureSampleV1, mass: f64, maximum_slope: f64,
) -> Result<(), Error> {
    let slope = mass / sample.response.heat_capacity_p_j_k;
    if !positive(slope) || slope > maximum_slope || !positive(sample.response.density_kg_m3) {
        return Err(Error::Invalid("reacting caloric response violates the original source capacity bound"));
    }
    Ok(())
}

impl ThermalLaw for GasEquilibriumThermalStorage<'_> {
    fn storage_basis(&self) -> ThermalStorageBasis { ThermalStorageBasis::Enthalpy }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h, || false)
            .map(|state| state.temperature.value()).map_err(|_| DepositionError::MaterialDomain)
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h, || false)
            .map(|state| Some(state.response.density_kg_m3)).map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 { self.maximum_slope }
    fn specific_entropy_change_j_kg_k(
        &self, initial_h: f64, final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.thermal_entropy_change_j_kg_k(initial_h, final_h, || false)
            .map(Some).map_err(|_| DepositionError::MaterialDomain)
    }
}

#[cfg(test)]
mod tests;
