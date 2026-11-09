//! Finite circular molecular channels with diffuse thermal wall scattering.
//!
//! The original aperture law supplies transmitted species and their source
//! molecular energies. The default fully accommodating wall re-emits a wall-hit
//! molecule with that species' molecular energy at the wall temperature. An
//! optional per-collision accommodation probability retains the full repeated-
//! collision thermal history. Returned molecules also exchange heat with their
//! original reservoir; a conductance multiplier alone misses both effects.
//!
//! The model assumes stationary Maxwellian reservoirs, collisionless flight,
//! velocity-independent diffuse reflection, equilibrated reservoir internal modes,
//! no adsorption/reaction, and negligible flight time and channel gas storage.
//! It supplies no outgoing angular/momentum law or moving-wall work. Santeler's
//! transmission approximation and the optional ring-resolved Clausing solution
//! retain explicit numerical estimates, not certified accuracy bounds.
//! Molecular energy accommodation follows Chiggiato, *Vacuum Technology for
//! Superconducting Devices*, section 2.2.2, <https://arxiv.org/abs/1501.07162>.

mod geometry;
pub use geometry::{
    AccommodationResolution, ClausingControls, ClausingResolution, CylindricalMolecularGeometry,
    MolecularAccommodation,
};
/// Molecular escape through the same diffuse channel into a true vacuum.
pub mod vacuum;

use super::{
    ClosedGasParcel, EffusiveFlux, EffusiveInterface, GasEnergyError, R, SpeciesId, invalid, poll,
    positive, positive_product_ratio, signed_product, source_error, valid_provenance,
};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::Nasa9StandardStateModelV1;

/// Original gas-qualified entrance mouths plus a finite diffuse-wall channel.
///
/// The entrance interface describes open circular mouths: every species has
/// unit transmission. Channel conductance comes exclusively from the retained
/// geometry. The mouth's thinness is not used as the channel's length.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseMolecularChannel {
    mouth: EffusiveInterface,
    geometry: CylindricalMolecularGeometry,
    accommodation: Option<(MolecularAccommodation, String)>,
}

impl DiffuseMolecularChannel {
    pub fn try_new(
        mouth: EffusiveInterface,
        geometry: CylindricalMolecularGeometry,
    ) -> Result<Self, GasEnergyError> {
        if mouth
            .transmissions()
            .iter()
            .any(|row| row.probability() != 1.0)
            || mouth.regime().aperture_diameter_m() != geometry.diameter_m()
            || geometry.minimum_mean_free_path_m() > mouth.regime().minimum_mean_free_path_m()
        {
            return Err(invalid(
                "diffuse channel needs unit-transmission matching mouths and no stronger collision qualification",
            ));
        }
        Ok(Self {
            mouth,
            geometry,
            accommodation: None,
        })
    }

    /// Resolve partial thermal accommodation on this SAME circular geometry.
    ///
    /// The authored coefficient is the probability that each wall collision
    /// equilibrates translational AND internal energy with the stationary wall.
    /// The complementary event changes only direction, diffusely, preserving
    /// molecular speed and internal state. This common coefficient must cover
    /// every original gas species and be independent of incident energy, angle
    /// and prior collisions throughout the caller's physical qualification.
    /// It is not a momentum-accommodation coefficient or an inferred bulk
    /// material property. A real surface needs independent constitutive data.
    ///
    /// Repeated collisions give each path weight `1-(1-alpha)^N`; multiplying
    /// the fully accommodating heat by alpha would miss that history. The
    /// Interior coefficients refine geometry and heat weights together; exact
    /// zero/unit coefficients retain the supplied geometry and its estimates.
    /// Construction performs the solve once; subsequent caloric/flux sampling
    /// reuses its immutable probabilities without further geometric solves.
    pub fn with_thermal_accommodation(
        mut self,
        coefficient: f64,
        controls: ClausingControls,
        provenance: impl Into<String>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        poll(&mut cancelled)?;
        let provenance = provenance.into();
        if !valid_provenance(&provenance) {
            return Err(invalid(
                "partial molecular accommodation needs bounded explicit surface qualification",
            ));
        }
        let (geometry, accommodation) =
            self.geometry
                .resolve_thermal_accommodation(coefficient, controls, &mut cancelled)?;
        self.geometry = geometry;
        self.accommodation = Some((accommodation, provenance));
        poll(&mut cancelled)?;
        Ok(self)
    }

    pub const fn mouth(&self) -> &EffusiveInterface {
        &self.mouth
    }
    pub const fn geometry(&self) -> &CylindricalMolecularGeometry {
        &self.geometry
    }

    /// Resolved partial-accommodation weights; None retains full accommodation.
    pub fn thermal_accommodation(&self) -> Option<&MolecularAccommodation> {
        self.accommodation.as_ref().map(|(model, _)| model)
    }
    pub fn thermal_accommodation_provenance(&self) -> Option<&str> {
        self.accommodation
            .as_ref()
            .map(|(_, source)| source.as_str())
    }

    fn heat_probabilities(&self) -> [f64; 2] {
        match self.thermal_accommodation() {
            Some(model) => [
                model.reflected_heat_probability(),
                model.transmitted_heat_probability(),
            ],
            None => [
                self.geometry.reflected_probability(),
                self.geometry.wall_transmitted_probability(),
            ],
        }
    }

    fn wall_heat_active(&self) -> bool {
        self.heat_probabilities()
            .into_iter()
            .any(|weight| weight > 0.0)
    }

    /// Preserve both the original mouth qualification and the full channel's
    /// collisional length at each actual gas state, including accepted endpoints.
    pub fn admit_parcel(
        &self,
        parcel: &ClosedGasParcel,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        self.mouth.admit_parcel(parcel, &mut cancelled)?;
        self.mouth.admit_collision_length(
            parcel,
            self.geometry.length_m().max(self.geometry.diameter_m()),
            &mut cancelled,
        )
    }

    /// Re-admit an exposed wall at the actual gas and wall endpoint without
    /// sampling another molecular beam or spending inventory.
    ///
    /// The original gas, collision and species sources remain authoritative.
    /// Only species present in this parcel need a caloric path from the gas
    /// temperature to the wall temperature. An absent catalog species does not
    /// narrow an active source's wall range. A zero-length or exactly elastic
    /// channel admits its gas but never queries an inactive thermal wall;
    /// callers skip geometrically closed channels entirely.
    pub fn admit_wall_for_parcel(
        &self,
        parcel: &ClosedGasParcel,
        wall_temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        self.admit_parcel(parcel, &mut cancelled)?;
        if !self.wall_heat_active() {
            return poll(&mut cancelled);
        }
        self.admit_wall_temperature(wall_temperature_k)?;
        let temperature = Temperature::new(parcel.state().state().temperature);
        let mixture = parcel.model().thermodynamics();
        for (species, &fraction) in mixture
            .components()
            .iter()
            .zip(mixture.mole_composition().fractions())
        {
            poll(&mut cancelled)?;
            if fraction == 0.0 {
                continue;
            }
            let _ = self
                .original_species(species.species())?
                .molar_enthalpy_change(temperature, Temperature::new(wall_temperature_k))
                .map_err(source_error)?;
        }
        poll(&mut cancelled)
    }

    /// Equivalent area for the original thin-aperture mass/energy sampler.
    /// Actual pore diameter and length remain independent of aggregate open area.
    pub fn transmitted_area_m2(&self, area_m2: f64) -> Result<f64, GasEnergyError> {
        if !area_m2.is_finite() || area_m2 < 0.0 {
            return Err(invalid(
                "molecular channel area must be finite and nonnegative",
            ));
        }
        if area_m2 == 0.0 {
            return Ok(0.0);
        }
        positive_product_ratio(&[area_m2, self.geometry.transmission_probability()], &[])
    }

    /// Sample rates without changing either gas. Reflected mass is in kg/s and
    /// wall heat is in W; a time integrator owns the corresponding inventories.
    pub fn sample(
        &self,
        left: &ClosedGasParcel,
        right: &ClosedGasParcel,
        area_m2: f64,
        wall_temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiffuseChannelFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        let area = self.transmitted_area_m2(area_m2)?;
        if area > 0.0 {
            self.admit_parcel(left, &mut cancelled)?;
            self.admit_parcel(right, &mut cancelled)?;
        }
        let transmitted = left.effusive_flux_to(right, &self.mouth, area, &mut cancelled)?;
        let wall = self.integrate_wall_exchange(
            [left, right],
            &transmitted,
            wall_temperature_k,
            1.0,
            &mut cancelled,
        )?;
        // The shared accounting kernel integrated a unit time interval. Expose
        // its numerical values as rates, distinct from finite-step exchanges.
        let wall = ChannelWallRates {
            reflected_species: wall
                .reflected_species
                .into_iter()
                .map(|row| ChannelReflectedSpeciesRate {
                    species: row.species,
                    mass_rate_kg_s: row.mass_kg,
                })
                .collect(),
            heat_into_gas_w: wall.heat_into_gas_j,
            heat_into_wall_w: wall.heat_into_wall_j,
        };
        Ok(DiffuseChannelFlux { transmitted, wall })
    }

    /// Use the ORIGINAL transmitted beam and the same rounded integrated mass
    /// as the native effusion transaction. Never resample a second gas state.
    pub(crate) fn integrate_wall_exchange(
        &self,
        gases: [&ClosedGasParcel; 2],
        transmitted: &EffusiveFlux,
        wall_temperature_k: f64,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelWallExchange, GasEnergyError> {
        poll(&mut cancelled)?;
        if !positive(duration_s) || transmitted.interface != self.mouth {
            return Err(invalid(
                "channel exchange lost its original interface or positive interval",
            ));
        }
        let mut result = ChannelWallExchange {
            reflected_species: Vec::with_capacity(transmitted.species.len()),
            heat_into_gas_j: [0.0; 2],
            heat_into_wall_j: 0.0,
        };
        // A closed mouth or zero-length channel never queries an inactive wall.
        if transmitted.species.is_empty() || self.geometry.length_m() == 0.0 {
            return Ok(result);
        }
        if self.wall_heat_active() {
            self.admit_wall_temperature(wall_temperature_k)?;
        }
        let temperatures = gases.map(|gas| gas.state().state().temperature);
        let mut heat = [Dd::ZERO; 2];
        for row in &transmitted.species {
            poll(&mut cancelled)?;
            let original = self.original_species(&row.species)?;
            let mut reflected = [0.0; 2];
            for (side, beam) in [row.left_to_right, row.right_to_left]
                .into_iter()
                .enumerate()
            {
                if beam.mass_rate_kg_s == 0.0 {
                    continue;
                }
                // Match the original network owner's rate*duration rounding.
                let mass = beam.mass_rate_kg_s * duration_s;
                let wall = self.integrate_beam_wall_exchange(
                    original,
                    mass,
                    temperatures[side],
                    wall_temperature_k,
                )?;
                reflected[side] = wall.reflected_mass_kg;
                heat[side] = heat[side] + Dd::from_f64(wall.heat_into_source_j);
                heat[1 - side] = heat[1 - side] + Dd::from_f64(wall.heat_into_transmitted_j);
            }
            if reflected.iter().any(|&mass| mass > 0.0) {
                result.reflected_species.push(ChannelReflectedSpecies {
                    species: row.species.clone(),
                    mass_kg: reflected,
                });
            }
        }
        result.heat_into_gas_j = heat.map(Dd::to_f64);
        result.heat_into_wall_j = (-(heat[0] + heat[1])).to_f64();
        if result
            .heat_into_gas_j
            .iter()
            .chain([&result.heat_into_wall_j])
            .any(|value| !value.is_finite())
        {
            return Err(invalid("channel wall heat is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn admit_wall_temperature(&self, wall_temperature_k: f64) -> Result<(), GasEnergyError> {
        let support = self.mouth.support().temperature_k;
        if !positive(wall_temperature_k)
            || wall_temperature_k < support.0
            || wall_temperature_k > support.1
        {
            return Err(invalid(
                "channel wall temperature exceeds the original caloric support",
            ));
        }
        Ok(())
    }

    fn original_species(
        &self,
        species: &SpeciesId,
    ) -> Result<&Nasa9StandardStateModelV1, GasEnergyError> {
        let components = self.mouth.gas_catalog().thermodynamics().components();
        let index = components
            .binary_search_by(|original| original.species().cmp(species))
            .map_err(|_| invalid("channel species is missing its original caloric source"))?;
        Ok(&components[index])
    }

    // A single caloric owner for finite receivers AND true vacuum. The caller
    // supplies its already-rounded transmitted mass. Return heat belongs to
    // the source; transmitted heat belongs to the receiving gas or escaped
    // beam. Never subtract two rounded absolute formation energies.
    fn integrate_beam_wall_exchange(
        &self,
        original: &Nasa9StandardStateModelV1,
        transmitted_mass_kg: f64,
        source_temperature_k: f64,
        wall_temperature_k: f64,
    ) -> Result<BeamWallExchange, GasEnergyError> {
        if !positive(transmitted_mass_kg) {
            return Err(invalid("integrated channel beam is unrepresentable"));
        }
        let tau = self.geometry.transmission_probability();
        let scaled = |probability| {
            if probability == 0.0 {
                Ok(0.0)
            } else {
                positive_product_ratio(&[transmitted_mass_kg, probability], &[tau])
            }
        };
        let reflected_mass_kg = scaled(self.geometry.reflected_probability())?;
        let [reflected_heat, transmitted_heat] = self.heat_probabilities();
        // Elastic diffuse scattering still returns a GROSS molecular inventory.
        // Only its thermal exchange vanishes; never discard the return budget.
        if reflected_heat == 0.0 && transmitted_heat == 0.0 {
            return Ok(BeamWallExchange {
                reflected_mass_kg,
                heat_into_source_j: 0.0,
                heat_into_transmitted_j: 0.0,
            });
        }
        let accommodated_return_mass_kg = scaled(reflected_heat)?;
        let accommodated_mass_kg = scaled(transmitted_heat)?;
        let molar_mass = original.molar_mass().value();
        let dh = original
            .molar_enthalpy_change(
                Temperature::new(source_temperature_k),
                Temperature::new(wall_temperature_k),
            )
            .map_err(source_error)?
            .value()
            / molar_mass;
        let selection = (R / molar_mass) * (0.5 * (wall_temperature_k - source_temperature_k));
        let de = (Dd::from_f64(dh) - Dd::from_f64(selection)).to_f64();
        if !de.is_finite() {
            return Err(invalid("channel accommodation energy is unrepresentable"));
        }
        Ok(BeamWallExchange {
            reflected_mass_kg,
            heat_into_source_j: signed_product(accommodated_return_mass_kg, de)?,
            heat_into_transmitted_j: signed_product(accommodated_mass_kg, de)?,
        })
    }
}

struct BeamWallExchange {
    reflected_mass_kg: f64,
    heat_into_source_j: f64,
    heat_into_transmitted_j: f64,
}

/// Molecules that leave and return to the SAME original reservoir. They spend
/// its gross inventory allowance despite giving no net species transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelReflectedSpecies {
    pub species: SpeciesId,
    pub mass_kg: [f64; 2],
}

/// Wall exchange additional to original transmitted molecular energy, once.
/// Positive heat enters the named gas or the thermostat. Opposite wall heat
/// belongs to the caller's thermal owner; this object does not invent storage.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelWallExchange {
    pub reflected_species: Vec<ChannelReflectedSpecies>,
    pub heat_into_gas_j: [f64; 2],
    pub heat_into_wall_j: f64,
}

/// Gross mass rates that return to their original left/right reservoirs.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelReflectedSpeciesRate {
    pub species: SpeciesId,
    pub mass_rate_kg_s: [f64; 2],
}

/// Wall power additional to the transmitted molecular energy rate, once.
/// Positive heat enters the named gas or thermostat. The caller's thermal
/// owner receives the opposite wall power; no wall storage is invented here.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelWallRates {
    pub reflected_species: Vec<ChannelReflectedSpeciesRate>,
    pub heat_into_gas_w: [f64; 2],
    pub heat_into_wall_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelFlux {
    pub transmitted: EffusiveFlux,
    pub wall: ChannelWallRates,
}

#[cfg(test)]
mod tests;
