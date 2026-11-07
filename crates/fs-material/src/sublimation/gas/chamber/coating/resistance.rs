//! Separate deposit/support temperatures against one finite, evolving gas.
//!
//! Molecular streams use the original deposit temperature and gas state. Their
//! opposite energy enters the deposit; an explicit finite conductance exchanges
//! heat with the nonvolatile support at the accepted endpoint temperatures.
//! The original two-body caloric solver and single gas reconstruction own the
//! conservation transaction. Neither an effective heat capacity nor a second
//! latent heater is supplied.
//!
//! The ordinary entry points require positive condensed inventories. The
//! explicit bare-birth entry point can form a positive film from incoming
//! molecules at the original support temperature, without a mass seed. Complete
//! extinction, localized birth time, spatial temperature profiles and
//! thickness-dependent conductance require separate physical models. These are
//! first-order frozen-rate intervals, not a complete coating lifecycle.

/// Accuracy-controlled histories carrying both condensed temperatures and one gas.
pub mod advance;
/// Original conductivity integrals and the actual inner coating contact.
pub mod conductive;
/// Explicit heat/power inputs with the original conservative history owner.
pub mod driven;
/// Closed support/coating/gas phase equilibrium on the original sources.
pub mod equilibrium;
/// Live gas-film transport in the original adaptive deposition history.
pub mod film;
/// Mass-derived exposed area and material conductivity in that same adaptive history.
pub mod growing;

use std::collections::{BTreeMap, BTreeSet};

use fs_thermal::sublimation::coating::resistance::{
    ResistiveHeatStep, ResistiveInterfaceLaw, ResistiveSurfaceState, heat_resistive_surface,
    heat_resistive_surface_with_interface,
};

pub use super::GasSupportedChamberSolve;
use super::{
    ChamberSurfaceOwner, ClosedGasParcel, GasChamberCandidate, GasChamberStep,
    GasChamberSurfaceStep, GasSpeciesBoundaryCandidate, GasSpeciesSnapshot, GasSurfaceBinding,
    GasSurfaceError, PreparedSurface, SolidSurfaceState, SublimationError, SublimationLaw,
    SurfaceKinetics, SurfaceSnapshot, ThermalAccommodation, ThermochemicalGasModel,
    advance_chamber_with_owner, compatibility, condensed_thermal_target, invalid,
    molecular_stream_energy, nonnegative, phase_at, poll, positive,
    snapshot_surface_at_temperature_with_kinetics, validate_storage, validate_surface,
};

/// One coating with two independent original caloric coordinates.
///
/// The gas contacts only the deposit. `conductance_w_k` is the complete frozen
/// support/deposit thermal conductance, not a per-area coefficient. Zero means
/// thermal insulation. Both inverse charts must be continuous, nondecreasing
/// and retain their original material domains. The fixed binding offset aligns
/// only deposit enthalpy with its gas species; support enthalpy is independent.
/// A zero-mass deposit is admitted only by the explicit bare-birth entry point.
#[derive(Clone, Copy)]
pub struct GasResistiveSurface<'a> {
    pub id: u64,
    pub law: &'a dyn SublimationLaw,
    pub binding: &'a GasSurfaceBinding,
    pub state: ResistiveSurfaceState,
    pub conductance_w_k: f64,
    pub substrate_temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    pub deposit_temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
}

/// Explicit heat delivered during one interval and an optional original contact law.
///
/// Positive heat enters the named condensed body; negative heat removes energy.
/// These are externally supplied joules, never a second molecular/latent source.
/// `interface: None` retains the surface's authored scalar conductance. A supplied
/// passive law replaces that conductance and is evaluated at both endpoint
/// temperatures inside the original bounded caloric solve.
#[derive(Clone, Copy, Default)]
pub struct GasResistiveDrive<'a> {
    pub substrate_heat_j: f64,
    pub deposit_heat_j: f64,
    pub interface: Option<&'a dyn ResistiveInterfaceLaw>,
}

impl ChamberSurfaceOwner for GasResistiveSurface<'_> {
    fn id(&self) -> u64 {
        self.id
    }
    fn binding(&self) -> &GasSurfaceBinding {
        self.binding
    }
    fn area_m2(&self) -> f64 {
        self.state.deposit.area_m2
    }
    fn condensed_mass_kg(&self) -> f64 {
        self.state.deposit.mass_kg
    }
}

/// Accepted original substrate/deposit storage and its allocated molecular stream.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasResistiveSurfaceStep {
    pub id: u64,
    pub state: ResistiveSurfaceState,
    /// Original substrate temperatures before/after this interval [K].
    /// Deposit temperatures are in `exchange.surface_temperature_k`.
    pub substrate_temperature_k: [f64; 2],
    /// Positive from substrate to deposit; internal, never added to gas heat.
    pub interface_heat_j: f64,
    /// External sources, kept separate even when equal/opposite or net zero.
    pub supplied_substrate_heat_j: f64,
    pub supplied_deposit_heat_j: f64,
    pub substrate_energy_change_j: f64,
    /// Includes actual deposit storage and the fixed gas-datum conversion.
    pub deposit_energy_change_j: f64,
    /// Separate heat balances keep opposing local defects from hiding.
    pub substrate_energy_residual_j: f64,
    pub deposit_energy_residual_j: f64,
    pub heat_evaluations: usize,
    /// Includes BOTH substrate and deposit in its solid storage account.
    pub exchange: GasChamberSurfaceStep,
    /// Already included in exchange.energy_into_gas_j.
    pub reflected_heat_into_gas_j: f64,
}

/// All two-temperature surfaces and the ONE original gas endpoint are accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasResistiveChamberCandidate<'catalog> {
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    surfaces: Vec<GasResistiveSurfaceStep>,
    step: GasChamberStep,
}
impl<'catalog> GasResistiveChamberCandidate<'catalog> {
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    #[must_use]
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    #[must_use]
    pub fn surfaces(&self) -> &[GasResistiveSurfaceStep] {
        &self.surfaces
    }
    #[must_use]
    pub const fn step(&self) -> &GasChamberStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        Vec<GasResistiveSurfaceStep>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasChamberStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.surfaces, gas, catalog, self.step)
    }
}

fn snapshot(
    surface: &GasResistiveSurface<'_>,
    gas: GasSpeciesSnapshot,
    catalog: &ThermochemicalGasModel,
    film: Option<f64>,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    admit_births: bool,
    drive: GasResistiveDrive<'_>,
) -> Result<SurfaceSnapshot, GasSurfaceError> {
    let state = surface.state;
    if !positive(state.substrate_mass_kg)
        || !state.substrate_specific_enthalpy_j_kg.is_finite()
        || !(positive(state.deposit.mass_kg) || (admit_births && state.deposit.mass_kg == 0.0))
        || (drive.interface.is_none() && !nonnegative(surface.conductance_w_k))
    {
        return Err(invalid(
            "resistive deposition requires positive existing support/film masses and finite conductance",
        ));
    }
    validate_surface(surface.binding, state.deposit)?;
    let support_t = (surface.substrate_temperature)(state.substrate_specific_enthalpy_j_kg)?;
    let deposit_t = (surface.deposit_temperature)(state.deposit.specific_enthalpy_j_kg)?;
    if !positive(support_t) {
        return Err(invalid("original support temperature is invalid"));
    }
    let temperature = if state.deposit.mass_kg == 0.0 {
        // This h is a caloric reference, not a mass/energy seed. The bare
        // substrate supplies the initial molecular boundary temperature.
        let tolerance = 16.0 * f64::EPSILON * support_t.abs().max(deposit_t.abs());
        if !positive(deposit_t) || (support_t - deposit_t).abs() > tolerance {
            return Err(invalid(
                "bare deposit caloric reference must match the original support temperature",
            ));
        }
        support_t
    } else {
        deposit_t
    };
    let mut snapshot = snapshot_surface_at_temperature_with_kinetics(
        surface.law,
        surface.binding,
        state.deposit,
        gas,
        catalog,
        temperature,
        SurfaceKinetics {
            accommodation_coefficient: surface.binding.accommodation_coefficient,
            mass_transfer_coefficient_m_s: film,
        },
        duration_s,
        controls.chamber.surface,
    )?;
    // A gas/deposit equilibrium is not stationary if support contact or EITHER
    // external source can change its state. Equal/opposite sources cannot
    // cancel this admission check merely because their combined heat is zero.
    let contact_stationary =
        support_t == deposit_t || (drive.interface.is_none() && surface.conductance_w_k == 0.0);
    snapshot.stationary_equilibrium &=
        contact_stationary && drive.substrate_heat_j == 0.0 && drive.deposit_heat_j == 0.0;
    Ok(snapshot)
}

struct Storage {
    heat: ResistiveHeatStep,
    initial_substrate_temperature_k: f64,
    deposit_energy_change_j: f64,
    supplied_substrate_heat_j: f64,
    supplied_deposit_heat_j: f64,
}

fn prepare(
    surface: &GasResistiveSurface<'_>,
    snapshot: SurfaceSnapshot,
    catalog: &ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    reflected_heat_into_gas_j: f64,
    drive: GasResistiveDrive<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(PreparedSurface, Storage), GasSurfaceError> {
    let energy = molecular_stream_energy(snapshot, reflected_heat_into_gas_j)?;
    let dm = snapshot.transferred_mass_kg;
    let next_mass = snapshot.solid.mass_kg - dm;
    if snapshot.solid.mass_kg == 0.0 && !positive(next_mass) {
        return Err(invalid(
            "bare resistive birth requires positive representable incoming deposit mass",
        ));
    }
    let mass_residual = (next_mass - snapshot.solid.mass_kg) + dm;
    if !positive(next_mass)
        || (dm != 0.0 && next_mass == snapshot.solid.mass_kg)
        || !mass_residual.is_finite()
        || mass_residual.abs() > controls.chamber.surface.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "resistive deposit exhausted or its inventory increment is unrepresentable",
        ));
    }
    // First account for transferred material at its original deposit enthalpy.
    // Remaining energy heats the ACTUAL final film mass; the independent support
    // retains its supplied heat and equal/opposite contact. Formation references cancel
    // before this sensible/latent increment, exactly as in the isothermal owner.
    let heat_j = condensed_thermal_target(
        snapshot,
        snapshot.solid.specific_enthalpy_j_kg,
        energy.difference_j_kg,
        reflected_heat_into_gas_j,
    );
    // Retain the error-free sum's low part in the physical residual. An
    // externally supplied heat must not disappear into a larger molecular
    // heat target while the reported source account still claims it arrived.
    let (deposit_heat_j, input_roundoff_j) = fs_math::eft::two_sum(heat_j, drive.deposit_heat_j);
    if !deposit_heat_j.is_finite() || !input_roundoff_j.is_finite() {
        return Err(invalid(
            "combined molecular and supplied coating heat is unrepresentable",
        ));
    }
    let thermal_state = ResistiveSurfaceState {
        deposit: SolidSurfaceState {
            mass_kg: next_mass,
            ..surface.state.deposit
        },
        ..surface.state
    };
    let checkpoint = || {
        if cancelled() {
            Err(SublimationError::Cancelled)
        } else {
            Ok(())
        }
    };
    let mut heat = match drive.interface {
        Some(interface) => heat_resistive_surface_with_interface(
            surface.law,
            thermal_state,
            duration_s,
            interface,
            drive.substrate_heat_j,
            deposit_heat_j,
            controls.storage,
            surface.substrate_temperature,
            surface.deposit_temperature,
            checkpoint,
        ),
        None => heat_resistive_surface(
            surface.law,
            thermal_state,
            duration_s,
            surface.conductance_w_k,
            drive.substrate_heat_j,
            deposit_heat_j,
            controls.storage,
            surface.substrate_temperature,
            surface.deposit_temperature,
            checkpoint,
        ),
    }?;
    heat.deposit_energy_residual_j -= input_roundoff_j;
    heat.energy_residual_j -= input_roundoff_j;
    let temperature = heat.deposit_temperature_k;
    let deposit_h = heat.state.deposit.specific_enthalpy_j_kg;
    let phase = phase_at(surface.law, temperature, snapshot.gas.molar_mass)?;
    let vapor_h = catalog.species_specific_enthalpy_j_kg(&surface.binding.species, temperature)?;
    compatibility(
        vapor_h,
        deposit_h,
        surface.binding.solid_enthalpy_offset_j_kg,
        phase.latent_heat_j_kg,
        temperature,
        controls.chamber.surface.latent_compatibility_tolerance_j_kg,
    )?;
    let deposit_energy_change_j = heat.deposit_energy_change_j
        - dm * (snapshot.solid.specific_enthalpy_j_kg + surface.binding.solid_enthalpy_offset_j_kg);
    let solid_energy_change_j = heat.substrate_energy_change_j + deposit_energy_change_j;
    let tolerance = controls.chamber.surface.gas.energy_tolerance_j;
    for residual in [
        heat.substrate_energy_residual_j,
        heat.deposit_energy_residual_j,
        heat.energy_residual_j,
    ] {
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(GasSurfaceError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
    }
    if !deposit_energy_change_j.is_finite() || !solid_energy_change_j.is_finite() {
        return Err(invalid("resistive deposit storage is unrepresentable"));
    }
    let recession_m = if surface.state.deposit.area_m2 == 0.0 {
        0.0
    } else {
        dm / phase.solid_density_kg_m3 / surface.state.deposit.area_m2
    };
    if !recession_m.is_finite() {
        return Err(invalid("resistive deposit recession is unrepresentable"));
    }
    let initial_substrate_temperature_k =
        (surface.substrate_temperature)(surface.state.substrate_specific_enthalpy_j_kg)?;
    Ok((
        PreparedSurface {
            solid: heat.state.deposit,
            surface_temperature_k: [snapshot.temperature_k, temperature],
            outgoing_mass_kg: snapshot.outgoing_mass_kg,
            incoming_mass_kg: snapshot.incoming_mass_kg,
            transferred_mass_kg: dm,
            outgoing_specific_energy_j_kg: snapshot.outgoing_specific_energy_j_kg,
            incoming_specific_energy_j_kg: snapshot.incoming_specific_energy_j_kg,
            energy_into_gas_j: energy.energy_into_gas_j,
            phase_relative_energy_into_gas_j: energy.phase_relative_energy_into_gas_j,
            relative_energy_into_gas_j: energy.relative_energy_into_gas_j,
            supplied_heat_j: drive.substrate_heat_j + drive.deposit_heat_j,
            solid_energy_change_j,
            energy_residual_j: heat.energy_residual_j,
            recession_m,
        },
        Storage {
            heat,
            initial_substrate_temperature_k,
            deposit_energy_change_j,
            supplied_substrate_heat_j: drive.substrate_heat_j,
            supplied_deposit_heat_j: drive.deposit_heat_j,
        },
    ))
}

/// Exchange existing two-temperature coatings with one finite original gas.
///
/// Gross streams share the original species-inventory limits; substrate mass
/// never enters a vapor allowance. Internal support heat disables the stationary
/// gas/deposit-equilibrium shortcut. Reflection is evaluated at the deposit
/// temperature and enters that deposit's heat balance exactly once.
/// Sources, inventories, local energy defects and cancellation are admitted
/// before any candidate is returned. An empty surface roster preserves gas.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveChamberCandidate<'catalog>, GasSurfaceError> {
    advance(
        surfaces,
        accommodation,
        None,
        None,
        gas,
        catalog,
        duration_s,
        controls,
        false,
        cancelled,
    )
}

/// The same conservative interval with explicitly supplied gas-film resistance.
///
/// A missing coefficient retains the original free-molecular limit; zero seals
/// molecular mass exchange while internal support/deposit heat can still flow.
/// Unknown IDs or invalid coefficients refuse. Thermal and mass-transfer
/// conductances are independent physical inputs and are never inferred from one another.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_with_films<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: &BTreeMap<u64, f64>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveChamberCandidate<'catalog>, GasSurfaceError> {
    advance(
        surfaces,
        accommodation,
        Some(films),
        None,
        gas,
        catalog,
        duration_s,
        controls,
        false,
        cancelled,
    )
}

/// Advance one interval with explicit external heating/cooling and contact laws.
///
/// Missing drive IDs retain zero supplied heat and the surface's scalar
/// conductance. Each heat is an interval amount in joules; zero duration may
/// apply a heat impulse but cannot transfer phase mass or contact heat. Both
/// external heat inputs are exposed separately in surface reports, and their
/// total in every joint source-aware energy receipt. The gas receives only its
/// original molecular/reflected energy transaction. Positive films are required.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_driven<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: Option<&BTreeMap<u64, f64>>,
    drives: &BTreeMap<u64, GasResistiveDrive<'_>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveChamberCandidate<'catalog>, GasSurfaceError> {
    advance(
        surfaces,
        accommodation,
        films,
        Some(drives),
        gas,
        catalog,
        duration_s,
        controls,
        false,
        cancelled,
    )
}

/// Form a thermally resistive film on explicitly admitted bare substrates.
///
/// A zero-mass row must carry the original deposit enthalpy at its substrate
/// temperature: its inverse chart is checked against that temperature. This
/// caloric reference stores no energy. The original bare-interface molecular
/// snapshot admits incoming molecules and emits no nonexistent deposit. Their
/// full energy heats the ACTUAL positive final film mass, then the original
/// finite-conductance solver shares heat with the support. Reflected molecular
/// heat participates in that same newborn-film balance exactly once.
///
/// Every bare row must receive a positive representable deposit in this
/// interval. Zero duration, zero area, sealed films or absent vapor therefore
/// refuse a bare birth; heating a still-bare substrate requires its own owner.
/// Positive-film rows retain ordinary behavior. Optional film coefficients and
/// every shared gas/inventory budget have the same meaning as above.
///
/// This explicitly freezes bare exposure over one caller-chosen interval. It
/// does NOT locate a nucleation/birth time or account for re-emission during
/// that first interval. The accepted positive film can then enter ordinary
/// resistive evolution. Do not use this operation as a generic equal-half-step
/// adaptive birth probe: a newborn's gross withdrawal fraction need not shrink
/// when both its age and the next trial interval are halved. No artificial mass
/// seed, evaporation of support, complete extinction or barrier law is added.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_with_births<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: Option<&BTreeMap<u64, f64>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveChamberCandidate<'catalog>, GasSurfaceError> {
    advance(
        surfaces,
        accommodation,
        films,
        None,
        gas,
        catalog,
        duration_s,
        controls,
        true,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: Option<&BTreeMap<u64, f64>>,
    drives: Option<&BTreeMap<u64, GasResistiveDrive<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    admit_births: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveChamberCandidate<'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    if surfaces.len() > controls.chamber.maximum_surfaces {
        return Err(invalid(
            "resistive surface roster exceeds its declared work bound",
        ));
    }
    if films.is_some() || drives.is_some() {
        if films.is_some_and(|map| map.len() > surfaces.len())
            || drives.is_some_and(|map| map.len() > surfaces.len())
        {
            return Err(invalid(
                "too many resistive gas-film coefficients or drive bindings",
            ));
        }
        let mut ids = BTreeSet::new();
        for surface in surfaces {
            poll(&mut cancelled)?;
            ids.insert(surface.id);
        }
        if let Some(films) = films {
            for (id, value) in films {
                poll(&mut cancelled)?;
                if !ids.contains(id) || !nonnegative(*value) {
                    return Err(invalid(
                        "gas-film coefficient references an unknown resistive surface or is invalid",
                    ));
                }
            }
        }
        if let Some(drives) = drives {
            for (id, drive) in drives {
                poll(&mut cancelled)?;
                if !ids.contains(id)
                    || !drive.substrate_heat_j.is_finite()
                    || !drive.deposit_heat_j.is_finite()
                    || !(drive.substrate_heat_j + drive.deposit_heat_j).is_finite()
                {
                    return Err(invalid(
                        "external heat references an unknown resistive surface or is unrepresentable",
                    ));
                }
            }
        }
    }
    let mut storage = Vec::with_capacity(surfaces.len());
    let (candidate, reflection) = advance_chamber_with_owner(
        surfaces,
        accommodation,
        gas,
        catalog,
        duration_s,
        controls.chamber,
        || validate_storage(controls.storage),
        |surface, gas| {
            snapshot(
                surface,
                gas,
                catalog,
                films.and_then(|coefficients| coefficients.get(&surface.id).copied()),
                duration_s,
                controls,
                admit_births,
                drives
                    .and_then(|map| map.get(&surface.id))
                    .copied()
                    .unwrap_or_default(),
            )
        },
        |surface, snapshot, reflected, cancelled| {
            let (prepared, retained) = prepare(
                surface,
                snapshot,
                catalog,
                duration_s,
                controls,
                reflected,
                drives
                    .and_then(|map| map.get(&surface.id))
                    .copied()
                    .unwrap_or_default(),
                cancelled,
            )?;
            storage.push(retained);
            Ok(prepared)
        },
        |_, _| Ok(()),
        &mut cancelled,
    )?;
    let GasChamberCandidate {
        gas,
        surfaces: exchanges,
        step,
    } = candidate;
    let mut surfaces = Vec::with_capacity(exchanges.len());
    for (exchange, stored) in exchanges.into_iter().zip(storage) {
        poll(&mut cancelled)?;
        let heat = stored.heat;
        let reflected_heat_into_gas_j = reflection
            .get(&exchange.id)
            .map_or(0.0, |account| account.heat_into_gas_j);
        surfaces.push(GasResistiveSurfaceStep {
            id: exchange.id,
            state: heat.state,
            substrate_temperature_k: [
                stored.initial_substrate_temperature_k,
                heat.substrate_temperature_k,
            ],
            interface_heat_j: heat.interface_heat_j,
            supplied_substrate_heat_j: stored.supplied_substrate_heat_j,
            supplied_deposit_heat_j: stored.supplied_deposit_heat_j,
            substrate_energy_change_j: heat.substrate_energy_change_j,
            deposit_energy_change_j: stored.deposit_energy_change_j,
            substrate_energy_residual_j: heat.substrate_energy_residual_j,
            deposit_energy_residual_j: heat.deposit_energy_residual_j,
            heat_evaluations: heat.evaluations,
            exchange,
            reflected_heat_into_gas_j,
        });
    }
    poll(&mut cancelled)?;
    Ok(GasResistiveChamberCandidate {
        gas,
        surfaces,
        step,
    })
}

#[cfg(test)]
#[path = "resistance/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "resistance/birth_tests.rs"]
mod birth_tests;
