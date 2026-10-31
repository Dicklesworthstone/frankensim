//! Dry onset and complete dry-out on an explicitly fully-wetting wall.
//!
//! Selecting this boundary admits a planar, continuously wettable interface
//! whose liquid is thermally equilibrated to the prescribed wall temperature.
//! The temperature comes from the ORIGINAL liquid chart at the supplied
//! interface enthalpy. This is an ideal wetting boundary, not a nucleation,
//! contact-angle, Kelvin-curvature, boiling, or moving-film model.
//!
//! The original Stefan-film rate uses the initial gas state. A dry wall cannot
//! evaporate; supersaturated vapor can create a real positive liquid inventory.
//! Evaporation stops exactly at m/(j*A) under that frozen-rate approximation,
//! after which the film is Dry for the rest of the REQUESTED interval. No seed,
//! epsilon film, liquid-fraction asymptote or renewed vapor donor budget is used.
//! Refine the full physical time step to converge the rate and event location.
//!
//! Wall heat is an explicit port: positive heat enters liquid plus gas. Existing
//! liquid is first equilibrated to the prescribed interface h, then phase
//! transfer carries upstream vapor enthalpy. NASA increments and the original
//! phase-law latent gap cancel formation references before evaluating that heat.
//! Sensible gas-wall convection belongs to a separate declared thermal owner.
//! A finite wall must debit the returned heat through its own caloric state.

use super::{
    ClosedGasParcel, EvaporationError, EvaporationLaw, GasEvaporationError, GasEvaporationSolve,
    GasFilm, GasLiquidBinding, GasSpeciesBoundaryCandidate, SpeciesBoundaryRelativeTransfer,
    ThermochemicalGasModel, admit_liquid_source, check_interface_endpoint, nonnegative, poll,
    positive, smaller_step, surface_flux, validate_controls,
};

/// No temperature, heat capacity or stored enthalpy is attached to a dry film.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GasWallFilmInventory {
    Dry,
    Wet {
        mass_kg: f64,
        /// Original liquid datum, separately from the explicit gas alignment.
        specific_enthalpy_j_kg: f64,
    },
}
impl GasWallFilmInventory {
    pub const fn mass_kg(self) -> f64 {
        match self {
            Self::Dry => 0.0,
            Self::Wet { mass_kg, .. } => mass_kg,
        }
    }
    pub const fn specific_enthalpy_j_kg(self) -> Option<f64> {
        match self {
            Self::Dry => None,
            Self::Wet {
                specific_enthalpy_j_kg,
                ..
            } => Some(specific_enthalpy_j_kg),
        }
    }
}

/// Explicit opt-in to the fully-wetting, thermally equilibrated wall closure.
/// Geometry remains on the wall when the liquid inventory becomes exactly Dry.
#[derive(Clone, Copy)]
pub struct GasWallFilmSurface<'a> {
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub film: GasWallFilmInventory,
    pub area_m2: f64,
    /// Liquid-chart h(T_wall), NOT the enthalpy of the wall material.
    pub interface_liquid_enthalpy_j_kg: f64,
    pub temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
    pub mass_transfer_coefficient_m_s: f64,
}

/// Full requested interval, with an exact event for its frozen Stefan rate.
#[derive(Clone, Debug, PartialEq)]
pub struct GasWallFilmStep {
    pub duration_s: f64,
    /// Time during which the nonzero signed phase flux acts.
    pub phase_active_duration_s: f64,
    pub dry_out_time_s: Option<f64>,
    pub film_formed: bool,
    pub liquid_mass_kg: [f64; 2],
    pub wall_temperature_k: f64,
    pub gas_temperature_k: [f64; 2],
    pub vapor_partial_pressure_pa: [f64; 2],
    /// Potential transfer if liquid were available throughout the full interval.
    pub requested_transferred_mass_kg: f64,
    /// Positive evaporation; negative condensation.
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    pub energy_into_gas_j: f64,
    /// Existing film's thermal equilibration to the prescribed wall [J].
    pub film_sensible_energy_j: f64,
    /// Actual represented thermal storage plus intended mass-reference change;
    /// the independently bounded mass-reference rounding follows the gas owner.
    pub film_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    /// Required wall heat into film plus gas; negative condensation heats wall.
    pub wall_heat_into_system_j: f64,
    pub energy_residual_j: f64,
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub liquid_mass_residual_kg: f64,
    /// Positive thickness loss; this report does not change wall geometry.
    pub recession_m: f64,
}

/// Joint finite gas, zero-capable film inventory and required wall heat port.
#[derive(Clone, Debug, PartialEq)]
pub struct GasWallFilmCandidate<'catalog> {
    film: GasWallFilmInventory,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    step: GasWallFilmStep,
}
impl<'catalog> GasWallFilmCandidate<'catalog> {
    pub const fn film(&self) -> GasWallFilmInventory {
        self.film
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    pub const fn step(&self) -> &GasWallFilmStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        GasWallFilmInventory,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasWallFilmStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.film, gas, catalog, self.step)
    }
}
fn invalid(message: &'static str) -> GasEvaporationError {
    GasEvaporationError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("wall-film arithmetic is unrepresentable"))
    }
}

/// Advance a fully-wetting prescribed wall through one frozen-rate interval.
///
/// A Wet film may transfer its ENTIRE mass at a represented dry-out time. The
/// existing inventory-fraction control applies to vapor withdrawal only here;
/// gas temperature and equilibrium overshoot guards still apply to evaporation.
/// A zero duration is strict identity, including no film thermal equilibration.
/// Failure or late cancellation publishes neither inventory nor wall heat.
pub fn advance_wall_film_with_gas<'catalog>(
    surface: GasWallFilmSurface<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasEvaporationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasWallFilmCandidate<'catalog>, GasEvaporationError> {
    poll(&mut cancelled)?;
    validate_controls(controls, duration_s, surface.mass_transfer_coefficient_m_s)?;
    if !nonnegative(surface.area_m2) || !surface.interface_liquid_enthalpy_j_kg.is_finite() {
        return Err(invalid(
            "finite wall area and original interface liquid enthalpy required",
        ));
    }
    let h = surface.interface_liquid_enthalpy_j_kg;
    let wall_t = (surface.temperature_from_enthalpy)(h)?;
    let source = admit_liquid_source(surface.law, surface.binding, catalog, wall_t, h, controls)?;
    let old_mass = surface.film.mass_kg();
    let delta_h = match surface.film {
        GasWallFilmInventory::Dry => 0.0,
        GasWallFilmInventory::Wet {
            mass_kg,
            specific_enthalpy_j_kg,
        } => {
            if !positive(mass_kg) || !specific_enthalpy_j_kg.is_finite() {
                return Err(invalid(
                    "Wet film needs positive finite mass and original enthalpy",
                ));
            }
            let initial_t = (surface.temperature_from_enthalpy)(specific_enthalpy_j_kg)?;
            admit_liquid_source(
                surface.law,
                surface.binding,
                catalog,
                initial_t,
                specific_enthalpy_j_kg,
                controls,
            )?;
            if duration_s > 0.0
                && (wall_t - initial_t).abs() > controls.maximum_temperature_increment_k
            {
                return Err(smaller_step(
                    duration_s,
                    "wall-film equilibration temperature departure",
                ));
            }
            if duration_s == 0.0 {
                0.0
            } else {
                finite(h - specific_enthalpy_j_kg)?
            }
        }
    };
    let gas_t = gas.state().state().temperature;
    let partial_pressure = gas.species_partial_pressure_pa(&surface.binding.species)?;
    let flux = surface_flux(
        surface.law,
        wall_t,
        GasFilm {
            total_pressure_pa: gas.state().state().pressure,
            vapor_partial_pressure_pa: partial_pressure,
            temperature_k: gas_t,
            mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
        },
    )?;
    let exposure = finite(surface.area_m2 * duration_s)?;
    let requested = finite(exposure * flux.mass_flux_kg_m2_s)?;
    if (duration_s > 0.0 && surface.area_m2 > 0.0 && exposure == 0.0)
        || (exposure > 0.0 && flux.mass_flux_kg_m2_s != 0.0 && requested == 0.0)
    {
        return Err(invalid("nonzero wall-film exposure or transfer underflows"));
    }
    let mut dm = requested;
    let mut active_duration = if dm == 0.0 { 0.0 } else { duration_s };
    let mut dry_out_time = None;
    if requested > 0.0 {
        if matches!(surface.film, GasWallFilmInventory::Dry) {
            dm = 0.0;
            active_duration = 0.0;
        } else {
            let rate = finite(surface.area_m2 * flux.mass_flux_kg_m2_s)?;
            if rate <= 0.0 {
                return Err(invalid("positive wall-film rate is unrepresentable"));
            }
            let event = old_mass / rate;
            if event <= duration_s {
                let event_mass_residual = finite(rate * event - old_mass)?;
                if !positive(event) || event_mass_residual.abs() > controls.gas.mass_tolerance_kg {
                    return Err(invalid(
                        "exact wall-film exhaustion event is unrepresentable",
                    ));
                }
                dm = old_mass;
                active_duration = event;
                dry_out_time = Some(event);
            } else if requested >= old_mass {
                return Err(invalid(
                    "wall-film time and inventory round to inconsistent exhaustion",
                ));
            }
        }
    }
    let old_vapor_mass = gas.species_mass_kg(&surface.binding.species)?;
    if dm < 0.0
        && (old_vapor_mass <= 0.0 || -dm / old_vapor_mass > controls.maximum_inventory_fraction)
    {
        let duration = duration_s * (controls.maximum_inventory_fraction * old_vapor_mass / -dm);
        return Err(GasEvaporationError::StepSizeLimit {
            maximum_duration_s: duration.next_down().max(0.0),
            reason: "initial wall-film vapor inventory withdrawal limit",
        });
    }
    let next_mass = if dry_out_time.is_some() {
        0.0
    } else {
        finite(old_mass - dm)?
    };
    let liquid_mass_residual = finite((next_mass - old_mass) + dm)?;
    if !nonnegative(next_mass)
        || (dm != 0.0 && next_mass == old_mass)
        || (next_mass == 0.0 && old_mass > 0.0 && dry_out_time.is_none())
        || liquid_mass_residual.abs() > controls.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "actual wall-film inventory cannot represent the phase event",
        ));
    }
    let vapor_delta_h = if dm == 0.0 {
        0.0
    } else {
        catalog.species_specific_enthalpy_change_j_kg(&surface.binding.species, gas_t, wall_t)?
    };
    let carried_gap = finite(
        source.phase.latent_heat_j_kg + source.alignment
            - if dm < 0.0 { vapor_delta_h } else { 0.0 },
    )?;
    let sensible = finite(old_mass * delta_h)?;
    let phase_wall_heat = finite(dm * carried_gap)?;
    if (old_mass > 0.0 && delta_h != 0.0 && sensible == 0.0)
        || (dm != 0.0 && carried_gap != 0.0 && phase_wall_heat == 0.0)
    {
        return Err(invalid("nonzero wall-film heat is unrepresentable"));
    }
    let (wall_heat, wall_heat_roundoff) = fs_math::eft::two_sum(sensible, phase_wall_heat);
    finite(wall_heat)?;
    finite(wall_heat_roundoff)?;
    let relative_energy = if dm > 0.0 {
        finite(dm * vapor_delta_h)?
    } else {
        0.0
    };
    poll(&mut cancelled)?;
    let endpoint = gas.exchange_species_boundary_relative(
        catalog,
        &[SpeciesBoundaryRelativeTransfer {
            species: surface.binding.species.clone(),
            mass_into_gas_kg: dm,
            relative_energy_into_gas_j: relative_energy,
        }],
        controls.gas,
        &mut cancelled,
    )?;
    let final_pressure = check_interface_endpoint(
        surface.law,
        surface.binding,
        wall_t,
        gas_t,
        dm,
        endpoint.gas(),
        surface.mass_transfer_coefficient_m_s,
        duration_s,
        controls,
    )?;
    let offset = surface.binding.liquid_enthalpy_offset_j_kg;
    // This form is exactly -m0*h0 at exhaustion; it does not subtract a
    // fictitious final film enthalpy or retain a zero-mass caloric state.
    let film_energy = match surface.film {
        GasWallFilmInventory::Dry => finite(-dm * (h + offset))?,
        GasWallFilmInventory::Wet {
            specific_enthalpy_j_kg,
            ..
        } => finite(next_mass * delta_h - dm * (specific_enthalpy_j_kg + offset))?,
    };
    let film_energy_residual = finite(liquid_mass_residual * delta_h)?;
    // Qwall is an actual port supplied to a different thermal owner. Any
    // addition residual lost while publishing that port is an energy defect,
    // even when opposing sensible/phase terms nearly cancel.
    let (partial_residual, partial_roundoff) =
        fs_math::eft::two_sum(film_energy_residual, wall_heat_roundoff);
    let (residual, residual_roundoff) =
        fs_math::eft::two_sum(partial_residual, endpoint.step().energy_residual_j);
    let residual = finite(residual + (partial_roundoff + residual_roundoff))?;
    let absolute_residual = finite(
        film_energy_residual.abs()
            + wall_heat_roundoff.abs()
            + endpoint.step().energy_residual_j.abs(),
    )?;
    let mass_residual = finite(
        (next_mass - old_mass)
            + (endpoint.gas().species_mass_kg(&surface.binding.species)? - old_vapor_mass),
    )?;
    let total_mass_residual =
        finite((next_mass - old_mass) + (endpoint.gas().mass_kg() - gas.mass_kg()))?;
    if mass_residual.abs() > controls.gas.mass_tolerance_kg
        || total_mass_residual.abs() > controls.gas.mass_tolerance_kg
        || liquid_mass_residual.abs() + endpoint.step().mass_residual_kg.abs()
            > controls.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "wall-film and actual gas inventory fail their shared mass allowance",
        ));
    }
    if absolute_residual > controls.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: absolute_residual,
            tolerance_j: controls.gas.energy_tolerance_j,
        });
    }
    let film = if next_mass == 0.0 {
        GasWallFilmInventory::Dry
    } else if duration_s == 0.0 {
        surface.film
    } else {
        GasWallFilmInventory::Wet {
            mass_kg: next_mass,
            specific_enthalpy_j_kg: h,
        }
    };
    let recession = if surface.area_m2 == 0.0 {
        0.0
    } else {
        finite(dm / source.phase.liquid_density_kg_m3 / surface.area_m2)?
    };
    let step = GasWallFilmStep {
        duration_s,
        phase_active_duration_s: active_duration,
        dry_out_time_s: dry_out_time,
        film_formed: matches!(surface.film, GasWallFilmInventory::Dry) && next_mass > 0.0,
        liquid_mass_kg: [old_mass, next_mass],
        wall_temperature_k: wall_t,
        gas_temperature_k: [gas_t, endpoint.gas().state().state().temperature],
        vapor_partial_pressure_pa: [partial_pressure, final_pressure],
        requested_transferred_mass_kg: requested,
        transferred_mass_kg: dm,
        latent_energy_j: finite(dm * source.phase.latent_heat_j_kg)?,
        energy_into_gas_j: endpoint.step().energy_into_gas_j,
        film_sensible_energy_j: sensible,
        film_enthalpy_change_j: film_energy,
        gas_internal_energy_change_j: endpoint.step().stored_energy_change_j,
        wall_heat_into_system_j: wall_heat,
        energy_residual_j: residual,
        sum_absolute_energy_residual_j: absolute_residual,
        mass_residual_kg: mass_residual,
        liquid_mass_residual_kg: liquid_mass_residual,
        recession_m: recession,
    };
    poll(&mut cancelled)?;
    Ok(GasWallFilmCandidate {
        film,
        gas: endpoint,
        step,
    })
}

#[cfg(test)]
mod tests;
