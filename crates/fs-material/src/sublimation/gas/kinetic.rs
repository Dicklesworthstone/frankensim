//! Local molecular exchange at a moving solid surface in a drifting ideal gas.
//!
//! The original gas species owns h(T), internal modes and the energy datum;
//! the original solid/vapor law owns saturation and sticking. The shared
//! shifted-Maxwellian owner supplies incident half-space moments, including
//! normal drift, tangential slip and the crossing energy. Explicit Maxwell
//! reflection supplies specular and thermalized diffuse branches. These are
//! complete molecular loads, including carrier pressure, not an extra stress
//! to add to another pressure boundary for the same molecules.
//!
//! This is an instantaneous constitutive rate, not a finite-reservoir update.
//! It does not evolve gas bulk momentum, body motion, volume or inventories.
//! An advancing owner must bound both gross phase streams, reconstruct the
//! original caloric stores, subtract its actual bulk kinetic-energy change
//! from laboratory energy, and handle geometry/volume work exactly once. No
//! continuum film, Knudsen-layer solution or non-Maxwellian gas is inferred.

use super::chamber::{
    GasChamberSurface,
    mixed::{CondensedInventory, SharpMeltingCaloric},
};
use super::momentum::{MaxwellReflection, MolecularSurfaceFrame};
use super::{
    ClosedGasParcel, GasSurfaceError, SpeciesId, ThermochemicalGasModel, VaporReservoir,
    compatibility, invalid, nonnegative, phase_at, poll, positive,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::admit_catalog;
use fs_qty::Temperature;
use fs_thermal::sublimation::molecular::{
    drift::{incident_maxwellian, stationary_crossing_energy_minus_enthalpy},
    molecular_mass_flux,
    momentum::flux_weighted_normal_speed,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use std::collections::BTreeMap;

/// Velocities in the same inertial SI frame as the surface centroid and origin.
/// The gas is Maxwellian about its declared bulk velocity. Surface rotation
/// can be represented locally by its velocity at this uniform patch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularInterfaceVelocity {
    pub gas_velocity_m_s: [f64; 3],
    pub surface_velocity_m_s: [f64; 3],
}

/// One original species' flux per unit exposed area. Positive transfers enter
/// the gas. Gross incident flux includes molecules that subsequently reflect.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticSpeciesFlux {
    pub species: SpeciesId,
    pub emitted_mass_flux_kg_m2_s: f64,
    pub incident_mass_flux_kg_m2_s: f64,
    pub sticking_mass_flux_kg_m2_s: f64,
    pub diffuse_reflected_mass_flux_kg_m2_s: f64,
    /// Evaluated without subtracting rounded gross beams near saturation.
    pub net_mass_flux_into_gas_kg_m2_s: f64,
    /// Momentum transfer using molecule velocities relative to the surface.
    pub relative_momentum_flux_into_gas_pa: [f64; 3],
    /// Relative momentum plus surface velocity times net transferred mass.
    pub laboratory_momentum_flux_into_gas_pa: [f64; 3],
    /// Surface-frame energy minus net mass times ORIGINAL h_i(T_gas).
    /// Formation constants cancel before evaluating this coordinate.
    pub relative_energy_flux_beyond_gas_enthalpy_w_m2: f64,
    /// Laboratory energy minus the SAME original partial enthalpy reference.
    pub laboratory_energy_flux_beyond_gas_enthalpy_w_m2: f64,
    /// Absolute laboratory energy, including the original formation datum.
    /// Use the source-relative coordinate for caloric reconstruction.
    pub laboratory_energy_flux_into_gas_w_m2: f64,
    /// Surface velocity dotted with the relative momentum transfer. This is
    /// mechanical wall power; advected bulk kinetic energy is reported apart.
    pub mechanical_power_into_gas_w_m2: f64,
    /// Net mass times half the squared surface velocity.
    pub transported_surface_kinetic_energy_flux_w_m2: f64,
}

/// Complete local rate, evaluated once from immutable original sources.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticInterfaceRate {
    pub surface_id: u64,
    pub frame: MolecularSurfaceFrame,
    pub origin_m: [f64; 3],
    pub velocity: MolecularInterfaceVelocity,
    pub surface_temperature_k: f64,
    pub gas_temperature_k: f64,
    /// Actual declared solid area; never inferred from a bulk phase fraction.
    pub exposed_solid_area_m2: f64,
    pub source_alignment_residual_j_kg: f64,
    /// Active gas species plus the emitting species, in canonical ID order.
    pub species: Vec<KineticSpeciesFlux>,
    pub mass_rate_into_gas_kg_s: f64,
    pub laboratory_momentum_rate_into_gas_n: [f64; 3],
    /// Force in the condensed body's variable-mass equation, excluding the
    /// transported bulk momentum. It is minus the relative molecular transfer.
    pub force_on_condensed_material_n: [f64; 3],
    pub torque_on_condensed_material_nm: [f64; 3],
    pub laboratory_energy_rate_into_gas_w: f64,
    pub mechanical_power_into_gas_w: f64,
    pub transported_surface_kinetic_energy_rate_w: f64,
}

#[derive(Default)]
struct Sum(f64, f64);
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), GasSurfaceError> {
        let (sum, error) = fs_math::eft::two_sum(self.0, finite(value)?);
        self.0 = finite(sum)?;
        self.1 = finite(self.1 + error)?;
        Ok(())
    }
    fn value(&self) -> Result<f64, GasSurfaceError> {
        finite(self.0 + self.1)
    }
}

fn finite(value: f64) -> Result<f64, GasSurfaceError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("molecular interface rate is unrepresentable"))
    }
}
fn product(a: f64, b: f64) -> Result<f64, GasSurfaceError> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(invalid("nonzero molecular interface rate underflowed"));
    }
    Ok(value)
}
fn sum(terms: impl IntoIterator<Item = f64>) -> Result<f64, GasSurfaceError> {
    let mut total = Sum::default();
    for term in terms {
        total.add(term)?;
    }
    total.value()
}
fn dot(a: [f64; 3], b: [f64; 3]) -> Result<f64, GasSurfaceError> {
    let mut total = Sum::default();
    for d in 0..3 {
        total.add(product(a[d], b[d])?)?;
    }
    total.value()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> Result<[f64; 3], GasSurfaceError> {
    let mut result = [0.0; 3];
    for d in 0..3 {
        let i = (d + 1) % 3;
        let j = (d + 2) % 3;
        result[d] = sum([product(a[i], b[j])?, -product(a[j], b[i])?])?;
    }
    Ok(result)
}

/// Evaluate one local molecular interface without changing either reservoir.
///
/// `caloric=None` uses the supplied pure-solid h/T callback. An explicit sharp
/// melting chart resolves bulk storage into its ORIGINAL solid constituent;
/// latent compatibility uses that solid h, never the mixture's bulk h. A body
/// with no remaining solid refuses. `scattering=None` means specular reflection
/// of nonsticking molecules. An explicit Maxwell map must cover the retained
/// species catalog; thermal accommodation alone does not infer angular physics.
/// All modes use one uniform surface velocity and the actual exposed solid area.
///
/// One distribution is evaluated per active species or emitting source; rate
/// storage is O(retained species), capped by the original gas model. Source
/// identity uses the existing boundary admission. Cancellation is checked
/// before work, per species and before publication.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_molecular_interface(
    surface: &GasChamberSurface<'_>,
    caloric: Option<SharpMeltingCaloric<'_>>,
    frame: MolecularSurfaceFrame,
    origin_m: [f64; 3],
    velocity: MolecularInterfaceVelocity,
    scattering: Option<&MaxwellReflection>,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    latent_compatibility_tolerance_j_kg: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<KineticInterfaceRate, GasSurfaceError> {
    poll(&mut cancelled)?;
    if frame.surface_id != surface.id
        || !positive(surface.solid.mass_kg)
        || !nonnegative(surface.solid.area_m2)
        || !surface.solid.specific_enthalpy_j_kg.is_finite()
        || !surface.binding.solid_enthalpy_offset_j_kg.is_finite()
        || !nonnegative(surface.binding.accommodation_coefficient)
        || surface.binding.accommodation_coefficient > 1.0
        || !nonnegative(latent_compatibility_tolerance_j_kg)
        || !frame
            .outward_normal
            .into_iter()
            .chain(frame.centroid_m)
            .chain(origin_m)
            .chain(velocity.gas_velocity_m_s)
            .chain(velocity.surface_velocity_m_s)
            .all(f64::is_finite)
    {
        return Err(invalid(
            "invalid moving molecular surface, binding or SI frame",
        ));
    }
    if (dot(frame.outward_normal, frame.outward_normal)? - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(invalid(
            "molecular surface normal must already be a unit vector",
        ));
    }
    admit_catalog(gas, catalog)?;
    let gas_temperature = gas.state().state().temperature;
    let surface_temperature =
        (surface.temperature_from_enthalpy)(surface.solid.specific_enthalpy_j_kg)?;
    let volatile = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|source| source.species() == &surface.binding.species)
        .ok_or(invalid(
            "moving phase interface requires its retained vapor source",
        ))?;
    let phase = phase_at(
        surface.law,
        surface_temperature,
        volatile.molar_mass().value(),
    )?;
    let solid_h = if let Some(source) = caloric {
        let inventory = CondensedInventory::new(source, surface.solid)?;
        if inventory.solid_mass_kg == 0.0 {
            return Err(invalid(
                "moving solid interface has no exposed solid inventory",
            ));
        }
        if (inventory.caloric.temperature_k() - surface_temperature).abs()
            > 64.0 * f64::EPSILON * surface_temperature.abs().max(1.0)
        {
            return Err(invalid(
                "surface callback disagrees with its original fusion chart",
            ));
        }
        source.admit_density(inventory.caloric, phase.solid_density_kg_m3)?;
        source.solid_enthalpy(inventory.caloric)
    } else {
        surface.solid.specific_enthalpy_j_kg
    };
    let vapor_h =
        catalog.species_specific_enthalpy_j_kg(&surface.binding.species, surface_temperature)?;
    let alignment = compatibility(
        vapor_h,
        solid_h,
        surface.binding.solid_enthalpy_offset_j_kg,
        phase.latent_heat_j_kg,
        surface_temperature,
        latent_compatibility_tolerance_j_kg,
    )?;
    if let Some(model) = scattering {
        model
            .thermal_accommodation()
            .admit_temperatures(surface_temperature, gas_temperature)?;
        let coefficients = model.thermal_accommodation().coefficients();
        if coefficients.len() != catalog.thermodynamics().components().len()
            || catalog
                .thermodynamics()
                .components()
                .iter()
                .any(|source| !coefficients.contains_key(source.species()))
        {
            return Err(invalid(
                "Maxwell scattering must cover the original retained species",
            ));
        }
    }
    let mut partial_pressures = BTreeMap::new();
    let active = gas.model().thermodynamics();
    for (source, &fraction) in active
        .components()
        .iter()
        .zip(active.mole_composition().fractions())
    {
        poll(&mut cancelled)?;
        if fraction != 0.0 {
            partial_pressures.insert(
                source.species(),
                product(fraction, gas.state().state().pressure)?,
            );
        }
    }
    let vapor_pressure = partial_pressures
        .get(&surface.binding.species)
        .copied()
        .unwrap_or(0.0);
    let stationary = molecular_mass_flux(
        surface.law,
        surface_temperature,
        surface.binding.accommodation_coefficient,
        VaporReservoir::Prescribed {
            partial_pressure_pa: vapor_pressure,
            temperature_k: gas_temperature,
        },
    )?;
    let relative_velocity =
        core::array::from_fn(|d| velocity.gas_velocity_m_s[d] - velocity.surface_velocity_m_s[d]);
    if !relative_velocity.into_iter().all(f64::is_finite) {
        return Err(invalid(
            "relative molecular bulk velocity is unrepresentable",
        ));
    }
    let normal_drift = dot(relative_velocity, frame.outward_normal)?;
    let mut rows = BTreeMap::new();
    for source in catalog.thermodynamics().components() {
        poll(&mut cancelled)?;
        let id = source.species();
        let is_volatile = id == &surface.binding.species;
        let pressure = partial_pressures.get(id).copied().unwrap_or(0.0);
        if pressure == 0.0 && !is_volatile {
            continue;
        }
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / source.molar_mass().value();
        let density = if pressure == 0.0 {
            0.0
        } else {
            // Divisions in two orders avoid an avoidable R*T overflow.
            let first = (pressure / r) / gas_temperature;
            let density = if positive(first) {
                first
            } else {
                (pressure / gas_temperature) / r
            };
            if !positive(density) {
                return Err(invalid("incident species density is unrepresentable"));
            }
            density
        };
        let incidence = incident_maxwellian(
            density,
            r,
            gas_temperature,
            relative_velocity,
            frame.outward_normal,
        )?;
        let incident = incidence.incoming_mass_flux_kg_m2_s;
        let beta = if is_volatile {
            surface.binding.accommodation_coefficient
        } else {
            0.0
        };
        let diffuse_fraction = scattering.map_or(0.0, |model| {
            model.thermal_accommodation().coefficients()[id]
        });
        let emitted = if is_volatile {
            stationary.outgoing_kg_m2_s
        } else {
            0.0
        };
        let sticking = product(beta, incident)?;
        let reflected = product(1.0 - beta, incident)?;
        let diffuse = product(diffuse_fraction, reflected)?;
        let specular = product(1.0 - diffuse_fraction, reflected)?;
        let net = if !is_volatile {
            0.0
        } else if normal_drift == 0.0 {
            stationary.net_kg_m2_s
        } else {
            let sigma = fs_math::det::sqrt(r) * fs_math::det::sqrt(gas_temperature);
            if normal_drift.abs() <= sigma {
                // Near equilibrium retain the small drift correction even when
                // incoming and outgoing gross beams round to the same value.
                sum([
                    stationary.net_kg_m2_s,
                    -product(
                        beta,
                        incidence.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
                    )?,
                ])?
            } else {
                sum([emitted, -sticking])?
            }
        };
        let emitted_speed = if emitted != 0.0 || diffuse != 0.0 {
            flux_weighted_normal_speed(r, surface_temperature)?
        } else {
            0.0
        };
        let (incident_velocity, incident_correction) =
            incidence.moments.map_or(([0.0; 3], 0.0), |moments| {
                (
                    moments.mean_relative_velocity_m_s,
                    moments.crossing_energy_minus_enthalpy_j_kg,
                )
            });
        let incident_normal = dot(incident_velocity, frame.outward_normal)?;
        let mut relative_momentum = [0.0; 3];
        let mut laboratory_momentum = [0.0; 3];
        for d in 0..3 {
            let outward = product(emitted_speed, frame.outward_normal[d])?;
            relative_momentum[d] = sum([
                product(emitted, outward)?,
                -product(sticking, incident_velocity[d])?,
                -product(
                    specular,
                    product(2.0 * incident_normal, frame.outward_normal[d])?,
                )?,
                product(diffuse, finite(outward - incident_velocity[d])?)?,
            ])?;
            laboratory_momentum[d] = sum([
                relative_momentum[d],
                product(net, velocity.surface_velocity_m_s[d])?,
            ])?;
        }
        let gas_h = finite(
            source
                .evaluate(Temperature::new(gas_temperature))
                .map_err(|_| {
                    invalid("original incident species caloric source refused temperature")
                })?
                .properties()
                .h()
                .value()
                / source.molar_mass().value(),
        )?;
        let thermalized = sum([sticking, diffuse])?;
        let relative_energy = if net == 0.0 && thermalized == 0.0 {
            0.0
        } else {
            let delta_h = finite(
                source
                    .molar_enthalpy_change(
                        Temperature::new(gas_temperature),
                        Temperature::new(surface_temperature),
                    )
                    .map_err(|_| invalid("original interface species caloric path refused"))?
                    .value()
                    / source.molar_mass().value(),
            )?;
            // The emitted beam's crossing energy below h(Ts), -R_i*Ts/2, is
            // fs-thermal's evaluation for a stationary incident beam, so at
            // Ts = Tg and zero drift the diffuse term below is exactly 0, as
            // accommodation.rs forms it.
            let outgoing_beyond_h = sum([
                delta_h,
                stationary_crossing_energy_minus_enthalpy(r, surface_temperature)?,
            ])?;
            let crossing_difference = sum([outgoing_beyond_h, -incident_correction])?;
            sum([
                product(net, outgoing_beyond_h)?,
                product(thermalized, crossing_difference)?,
            ])?
        };
        let mechanical_power = dot(velocity.surface_velocity_m_s, relative_momentum)?;
        let advected_kinetic = if net == 0.0 {
            0.0
        } else {
            // Multiply by net mass first so a representable flux does not
            // require a separately representable squared laboratory speed.
            let weighted = velocity.surface_velocity_m_s.map(|v| net * v);
            dot(weighted, velocity.surface_velocity_m_s.map(|v| 0.5 * v))?
        };
        let laboratory_relative_energy =
            sum([relative_energy, mechanical_power, advected_kinetic])?;
        let laboratory_energy = sum([product(net, gas_h)?, laboratory_relative_energy])?;
        rows.insert(
            id.clone(),
            KineticSpeciesFlux {
                species: id.clone(),
                emitted_mass_flux_kg_m2_s: emitted,
                incident_mass_flux_kg_m2_s: incident,
                sticking_mass_flux_kg_m2_s: sticking,
                diffuse_reflected_mass_flux_kg_m2_s: diffuse,
                net_mass_flux_into_gas_kg_m2_s: net,
                relative_momentum_flux_into_gas_pa: relative_momentum,
                laboratory_momentum_flux_into_gas_pa: laboratory_momentum,
                relative_energy_flux_beyond_gas_enthalpy_w_m2: relative_energy,
                laboratory_energy_flux_beyond_gas_enthalpy_w_m2: laboratory_relative_energy,
                laboratory_energy_flux_into_gas_w_m2: laboratory_energy,
                mechanical_power_into_gas_w_m2: mechanical_power,
                transported_surface_kinetic_energy_flux_w_m2: advected_kinetic,
            },
        );
    }
    let species: Vec<_> = rows.into_values().collect();
    let area = surface.solid.area_m2;
    let mut force = [0.0; 3];
    let mut lab_momentum = [0.0; 3];
    for d in 0..3 {
        force[d] = -product(
            area,
            sum(species
                .iter()
                .map(|row| row.relative_momentum_flux_into_gas_pa[d]))?,
        )?;
        lab_momentum[d] = product(
            area,
            sum(species
                .iter()
                .map(|row| row.laboratory_momentum_flux_into_gas_pa[d]))?,
        )?;
    }
    let arm = core::array::from_fn(|d| frame.centroid_m[d] - origin_m[d]);
    let result = KineticInterfaceRate {
        surface_id: surface.id,
        frame,
        origin_m,
        velocity,
        surface_temperature_k: surface_temperature,
        gas_temperature_k: gas_temperature,
        exposed_solid_area_m2: area,
        source_alignment_residual_j_kg: alignment,
        mass_rate_into_gas_kg_s: product(
            area,
            sum(species.iter().map(|row| row.net_mass_flux_into_gas_kg_m2_s))?,
        )?,
        laboratory_momentum_rate_into_gas_n: lab_momentum,
        force_on_condensed_material_n: force,
        torque_on_condensed_material_nm: cross(arm, force)?,
        laboratory_energy_rate_into_gas_w: product(
            area,
            sum(species
                .iter()
                .map(|row| row.laboratory_energy_flux_into_gas_w_m2))?,
        )?,
        mechanical_power_into_gas_w: product(
            area,
            sum(species.iter().map(|row| row.mechanical_power_into_gas_w_m2))?,
        )?,
        transported_surface_kinetic_energy_rate_w: product(
            area,
            sum(species
                .iter()
                .map(|row| row.transported_surface_kinetic_energy_flux_w_m2))?,
        )?,
        species,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
