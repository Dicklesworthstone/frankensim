//! Inventory-limited forward progress of one original, conserved reaction.
//!
//! The caller supplies a requested extent; available reactants, including O2
//! when present in the authored reaction, limit it. This is an instantaneous
//! mixing closure, not an Arrhenius integration, ignition/extinction model,
//! equilibrium calculation or calibrated fire law. Species are never inferred
//! from names. Every species remains an original ideal gas; condensed soot is
//! outside this owner. Reactions with a species on both sides are refused.
//!
//! Chemistry changes the sensible-energy gauge by `-sum(dm_i/dt * e_ref_i)`.
//! Physical internal energy is unchanged. The reported reaction enthalpy is a
//! diagnostic (or a requested-extent conversion), **not another heat source**.

use super::{Layer, LayerRates, MAX_SPECIES, MixtureGas, TwoZone, sum};
use crate::multicomponent::Error as GasError;
use fs_qty::{Amount, Energy, Power, Temperature, Time};
use fs_thermochem::kinetics::IdealGasKineticsV1;
use fs_thermochem::{MolarEnergyQuantityV1, ReactionId, SpeciesId};

#[cfg(test)]
mod tests;

/// A refusal returns no chemical rates and cannot alter the room or mechanism.
#[derive(Clone, Debug, PartialEq)]
pub enum BurningError {
    /// Invalid requested extent, duration, incoming rate or reference temperature.
    Invalid(&'static str),
    /// The complete original bank or energy gauge differs from the bound room.
    SourceMismatch,
    /// The selected reaction is absent or outside the single-pass closure.
    UnsupportedReaction(&'static str),
    /// A rate or common positive extent cannot be represented without depletion.
    Unrepresentable(&'static str),
    /// The original NASA source refused a thermodynamic evaluation.
    Gas(GasError),
    /// Cancellation was observed before a result was returned.
    Cancelled,
}
impl core::fmt::Display for BurningError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "inventory-limited reaction: {self:?}")
    }
}
impl std::error::Error for BurningError {}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BurningError> {
    if cancelled() {
        Err(BurningError::Cancelled)
    } else {
        Ok(())
    }
}

/// One time-averaged chemical source. No incoming mass/enthalpy is repeated here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurningReport {
    /// Chemical species rates and the corresponding sensible-energy correction.
    pub chemical_rates: LayerRates,
    /// Original caller request before reactant availability was applied.
    pub requested_extent: Amount,
    /// Admitted forward extent integrated over the supplied duration.
    pub actual_extent: Amount,
    /// Canonical species index setting the inventory/representability limit.
    /// Ties choose the first canonical species; a zero request has no limiter.
    pub limiting_species: Option<usize>,
    /// Explicit temperature at which the source reaction enthalpy was evaluated.
    pub reference_temperature: Temperature,
    /// Original `sum(nu_i h_i)` [J/mol of reaction extent], negative if exothermic.
    pub reaction_enthalpy: MolarEnergyQuantityV1,
    /// `-actual_extent * reaction_enthalpy`, for diagnostics only.
    pub released_enthalpy: Energy,
    /// The same source enthalpy release per supplied second, for diagnostics only.
    pub released_power: Power,
    /// Actual signed molecular-mass source discrepancy [kg/s], never repaired.
    pub mass_balance_rate_kg_s: f64,
}

/// Immutable adapter bound to an already conserved mechanism and original gas.
/// Evaluation allocates no workspace and changes neither owner nor room.
#[derive(Clone, Debug, PartialEq)]
pub struct InventoryLimitedReaction {
    gas: MixtureGas,
    mechanism: IdealGasKineticsV1,
    reaction_index: usize,
    coefficients: [f64; MAX_SPECIES],
    mass_per_extent: [f64; MAX_SPECIES],
}
impl InventoryLimitedReaction {
    /// Select one already element/charge/molecular-mass-checked source reaction.
    /// The mechanism's Arrhenius parameters are retained but never evaluated.
    pub fn new(
        gas: &MixtureGas,
        mechanism: IdealGasKineticsV1,
        reaction: &ReactionId,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BurningError> {
        poll(&mut cancelled)?;
        if gas.source() != mechanism.thermodynamics() {
            return Err(BurningError::SourceMismatch);
        }
        let reaction_index = mechanism
            .reactions()
            .binary_search_by(|r| r.reaction().cmp(reaction))
            .map_err(|_| BurningError::UnsupportedReaction("unknown reaction"))?;
        let selected = &mechanism.reactions()[reaction_index];
        if selected.reactants().iter().any(|(id, _)| {
            selected
                .products()
                .binary_search_by(|(product, _)| product.cmp(id))
                .is_ok()
        }) {
            return Err(BurningError::UnsupportedReaction(
                "overlapping reactant/product species",
            ));
        }
        let mut coefficients = [0.0; MAX_SPECIES];
        let mut mass_per_extent = [0.0; MAX_SPECIES];
        for (i, model) in gas.source().models().iter().enumerate() {
            poll(&mut cancelled)?;
            coefficients[i] = mechanism.stoichiometry().get(i, reaction_index).ok_or(
                BurningError::UnsupportedReaction("missing canonical stoichiometry"),
            )? as f64;
            mass_per_extent[i] = coefficients[i] * model.molar_mass().value();
            if !mass_per_extent[i].is_finite()
                || (coefficients[i] != 0.0 && mass_per_extent[i] == 0.0)
            {
                return Err(BurningError::Unrepresentable("stoichiometric mass"));
            }
        }
        let out = Self {
            gas: gas.clone(),
            mechanism,
            reaction_index,
            coefficients,
            mass_per_extent,
        };
        poll(&mut cancelled)?;
        Ok(out)
    }
    /// Complete original mechanism, including its exact conservation certificate.
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        &self.mechanism
    }
    /// Selected original reaction identity.
    pub fn reaction(&self) -> &ReactionId {
        self.mechanism.reactions()[self.reaction_index].reaction()
    }
    /// Canonical species identities used by the chemical rates and limiter index.
    pub fn species(&self) -> &[SpeciesId] {
        self.mechanism.stoichiometry().species()
    }

    /// Convert an authored enthalpy-release request through the original source
    /// at an explicit temperature. This requires a strictly exothermic reaction;
    /// actual available reactants are subsequently applied by [`Self::evaluate`].
    pub fn extent_for_power(
        &self,
        power: Power,
        duration: Time,
        reference_temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Amount, BurningError> {
        poll(&mut cancelled)?;
        if !power.value().is_finite()
            || power.value() < 0.0
            || !duration.value().is_finite()
            || duration.value() <= 0.0
        {
            return Err(BurningError::Invalid("power or duration"));
        }
        let h = self.enthalpy(reference_temperature, &mut cancelled)?;
        if h >= 0.0 {
            return Err(BurningError::Invalid(
                "power request requires an exothermic reaction",
            ));
        }
        let extent = (power.value() / -h) * duration.value();
        if !extent.is_finite() || (power.value() > 0.0 && extent == 0.0) {
            return Err(BurningError::Unrepresentable("requested enthalpy extent"));
        }
        poll(&mut cancelled)?;
        Ok(Amount::new(extent))
    }

    /// Limit one forward request by the selected layer's actual inventories.
    ///
    /// `incoming` supplies nonnegative species rates available during this same
    /// interval, so newly injected fuel may burn immediately. It is used only
    /// for availability: the caller must add that identical incoming rate and
    /// `chemical_rates` exactly once to the layer passed to `TwoZone::advance`.
    /// Other simultaneous withdrawals require their own combined admission.
    ///
    /// The final mass test uses the same `m + dt*(incoming + chemical)` arithmetic
    /// as the room. At most 16 predecessor rates admit a common representable
    /// extent; no individual species is clipped or artificially populated.
    pub fn evaluate(
        &self,
        room: &TwoZone,
        layer: Layer,
        requested_extent: Amount,
        duration: Time,
        incoming: Option<&LayerRates>,
        reference_temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BurningReport, BurningError> {
        poll(&mut cancelled)?;
        if room.gas() != &self.gas {
            return Err(BurningError::SourceMismatch);
        }
        let dt = duration.value();
        let mut extent = requested_extent.value();
        let n = self.species().len();
        if !dt.is_finite() || dt <= 0.0 || !extent.is_finite() || extent < 0.0 {
            return Err(BurningError::Invalid("extent or duration"));
        }
        let zero = LayerRates::default();
        let incoming = incoming.unwrap_or(&zero);
        if !incoming.referenced_energy_w.is_finite()
            || incoming
                .species_kg_s
                .iter()
                .any(|r| !r.is_finite() || *r < 0.0)
            || incoming.species_kg_s[n..].iter().any(|r| *r != 0.0)
        {
            return Err(BurningError::Invalid(
                "incoming species axis or nonnegative rates",
            ));
        }
        let masses = room.layer(layer).species_masses();
        let mut limiting_species = None;
        for (i, mass) in masses.iter().enumerate() {
            poll(&mut cancelled)?;
            let available = mass + dt * incoming.species_kg_s[i];
            if !available.is_finite() {
                return Err(BurningError::Unrepresentable("available species inventory"));
            }
            if self.mass_per_extent[i] < 0.0 {
                let cap = available / -self.mass_per_extent[i];
                if cap < extent || (cap == extent && extent > 0.0 && limiting_species.is_none()) {
                    extent = cap;
                    limiting_species = Some(i);
                }
            }
        }
        let reaction_enthalpy = self.enthalpy(reference_temperature, &mut cancelled)?;
        let mut progress_rate = extent / dt;
        if !progress_rate.is_finite() || (extent > 0.0 && progress_rate == 0.0) {
            return Err(BurningError::Unrepresentable("reaction progress rate"));
        }
        for _ in 0..16 {
            poll(&mut cancelled)?;
            let actual = progress_rate * dt;
            let mut rates = LayerRates::default();
            let mut reduce = actual > extent;
            for i in 0..n {
                let rate = self.mass_per_extent[i] * progress_rate;
                if !rate.is_finite()
                    || (self.mass_per_extent[i] != 0.0 && progress_rate > 0.0 && rate == 0.0)
                {
                    return Err(BurningError::Unrepresentable("chemical species rate"));
                }
                rates.species_kg_s[i] = rate;
                let final_mass = masses[i] + dt * (incoming.species_kg_s[i] + rate);
                if !final_mass.is_finite() {
                    return Err(BurningError::Unrepresentable("final species inventory"));
                }
                if final_mass < 0.0 {
                    reduce = true;
                    if limiting_species.is_none() {
                        limiting_species = Some(i);
                    }
                }
            }
            if reduce {
                progress_rate = progress_rate.next_down();
                if progress_rate <= 0.0 {
                    return Err(BurningError::Unrepresentable("positive feasible extent"));
                }
                continue;
            }
            rates.referenced_energy_w =
                -sum((0..n).map(|i| rates.species_kg_s[i] * self.gas.energy_references()[i]));
            let released_enthalpy = -actual * reaction_enthalpy;
            let released_power = -progress_rate * reaction_enthalpy;
            let mass_balance_rate_kg_s = sum(rates.species_kg_s[..n].iter().copied());
            if !rates.referenced_energy_w.is_finite()
                || !released_enthalpy.is_finite()
                || !released_power.is_finite()
                || !mass_balance_rate_kg_s.is_finite()
            {
                return Err(BurningError::Unrepresentable(
                    "chemical energy or mass balance",
                ));
            }
            poll(&mut cancelled)?;
            return Ok(BurningReport {
                chemical_rates: rates,
                requested_extent,
                actual_extent: Amount::new(actual),
                limiting_species,
                reference_temperature,
                reaction_enthalpy: MolarEnergyQuantityV1::new(reaction_enthalpy),
                released_enthalpy: Energy::new(released_enthalpy),
                released_power: Power::new(released_power),
                mass_balance_rate_kg_s,
            });
        }
        Err(BurningError::Unrepresentable(
            "bounded common-extent admission",
        ))
    }

    fn enthalpy(
        &self,
        temperature: Temperature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, BurningError> {
        let domain = self.gas.domain();
        if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&temperature.value())
        {
            return Err(BurningError::Invalid("source reference temperature"));
        }
        let mut terms = [0.0; MAX_SPECIES];
        for (i, model) in self.gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            if self.coefficients[i] != 0.0 {
                terms[i] = self.coefficients[i]
                    * model
                        .evaluate_properties(temperature)
                        .map_err(|e| BurningError::Gas(GasError::from(e)))?
                        .h()
                        .value();
            }
        }
        let value = sum(terms[..self.species().len()].iter().copied());
        if !value.is_finite() {
            return Err(BurningError::Unrepresentable("source reaction enthalpy"));
        }
        Ok(value)
    }
}
