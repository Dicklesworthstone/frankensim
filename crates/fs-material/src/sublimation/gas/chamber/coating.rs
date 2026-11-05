//! Supported deposition against one finite, nonisothermal original gas parcel.
//!
//! A nonvolatile substrate and its separately conserved deposit share one
//! temperature. A bare substrate can receive the explicitly admitted species;
//! it cannot emit a nonexistent coating. The substrate retains its original
//! mass and caloric coordinate, while only deposit mass enters the gas budget.
//! No numerical seed, nucleation barrier or substrate evaporation is inferred.
//!
//! Molecular streams, NASA relative energy, reflected heat, species inventory
//! limits and the single gas reconstruction belong to the original chamber.
//! The existing supported-storage solver distributes the opposite stream heat
//! between both original caloric charts. A deposit is never thermally massless.
//! Its fixed binding offset aligns the deposit chart with the gas species;
//! the substrate's energy reference remains independent.
//!
//! This is a thin, internally isothermal coating with fixed exposed area.
//! Rates use the OLD gas and surface state. Newly emitted vapor cannot fund
//! condensation in the same explicit interval; refine time to resolve transport.
//! An absent species can start in a finite carrier gas. Complete gas vacuum,
//! porous frost, evolving geometry and resolved coating resistance require
//! their own physical owners and are not silently supplied here.

use std::collections::BTreeMap;

use fs_thermal::sublimation::{
    SolverConfig, SurfaceKinetics,
    coating::{SupportedSurfaceState, heat_supported_surface},
};

use super::super::{
    compatibility, invalid, molecular_stream_energy, nonnegative, phase_at, poll, positive,
    snapshot_surface_at_temperature, snapshot_surface_at_temperature_with_kinetics,
    validate_surface,
};
use super::{
    ChamberSurfaceOwner, ClosedGasParcel, GasChamberCandidate, GasChamberSolve, GasChamberStep,
    GasChamberSurfaceStep, GasSpeciesBoundaryCandidate, GasSpeciesSnapshot, GasSurfaceBinding,
    GasSurfaceError, PreparedSurface, SolidSurfaceState, SublimationError, SublimationLaw,
    SurfaceSnapshot, ThermalAccommodation, ThermochemicalGasModel, advance_chamber_with_owner,
};

/// Original substrate/deposit storage and an explicitly admitted depositing species.
///
/// Both charts must be continuous and nondecreasing over their admitted solid
/// branch. `deposit_enthalpy(T)` retains its original datum; the binding supplies
/// its fixed conversion to the NASA species datum. `substrate_temperature(h)`
/// inverts the independent nonvolatile substrate chart. A zero deposit mass
/// admits incoming molecules without inventing a volatile substrate inventory.
#[derive(Clone, Copy)]
pub struct GasSupportedSurface<'a> {
    pub id: u64,
    pub law: &'a dyn SublimationLaw,
    pub binding: &'a GasSurfaceBinding,
    pub state: SupportedSurfaceState,
    pub substrate_temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    pub deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
}

impl ChamberSurfaceOwner for GasSupportedSurface<'_> {
    fn id(&self) -> u64 {
        self.id
    }
    fn binding(&self) -> &GasSurfaceBinding {
        self.binding
    }
    fn area_m2(&self) -> f64 {
        self.state.area_m2
    }
    fn condensed_mass_kg(&self) -> f64 {
        self.state.deposit_mass_kg
    }
}

/// Original gas/inventory controls and a bounded original supported-storage solve.
/// The chamber's absolute energy tolerance also bounds every final local defect.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSupportedChamberSolve {
    pub chamber: GasChamberSolve,
    pub storage: SolverConfig,
}

/// One accepted supported surface in canonical stable-id order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSupportedSurfaceStep {
    pub id: u64,
    pub state: SupportedSurfaceState,
    /// `solid` is the deposit alone. `solid_energy_change_j` and the caloric
    /// defect include BOTH substrate and deposit; stream energy includes latent
    /// and reflected heat exactly once.
    pub exchange: GasChamberSurfaceStep,
    /// Actual nonvolatile substrate energy change in its independent datum.
    pub substrate_energy_change_j: f64,
    /// Deposit storage change including the binding's fixed NASA datum offset.
    pub deposit_energy_change_j: f64,
    /// Reflected molecular heat INCLUDED in `exchange.energy_into_gas_j`.
    pub reflected_heat_into_gas_j: f64,
}

/// Publish every supported endpoint and the one original gas endpoint together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSupportedChamberCandidate<'catalog> {
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    surfaces: Vec<GasSupportedSurfaceStep>,
    step: GasChamberStep,
}
impl<'catalog> GasSupportedChamberCandidate<'catalog> {
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
    pub fn surfaces(&self) -> &[GasSupportedSurfaceStep] {
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
        Vec<GasSupportedSurfaceStep>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasChamberStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.surfaces, gas, catalog, self.step)
    }
}

fn validate_storage(config: SolverConfig) -> Result<(), GasSurfaceError> {
    // These are the existing thermal storage owner's controls, including its
    // hard iteration ceiling. Validate even an empty active selection.
    if !positive(config.relative_tolerance)
        || config.relative_tolerance >= 1.0
        || !nonnegative(config.absolute_energy_tolerance_j)
        || config.max_iterations == 0
        || config.max_iterations > 4096
    {
        return Err(SublimationError::Invalid("invalid sublimation nonlinear budget").into());
    }
    Ok(())
}

fn snapshot(
    surface: &GasSupportedSurface<'_>,
    gas: GasSpeciesSnapshot,
    catalog: &ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
) -> Result<SurfaceSnapshot, GasSurfaceError> {
    snapshot_with_film(surface, gas, catalog, None, duration_s, controls)
}

fn snapshot_with_film(
    surface: &GasSupportedSurface<'_>,
    gas: GasSpeciesSnapshot,
    catalog: &ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: Option<f64>,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
) -> Result<SurfaceSnapshot, GasSurfaceError> {
    if !positive(surface.state.substrate_mass_kg)
        || !surface.state.substrate_specific_enthalpy_j_kg.is_finite()
    {
        return Err(invalid(
            "finite positive original substrate storage required",
        ));
    }
    let temperature =
        (surface.substrate_temperature)(surface.state.substrate_specific_enthalpy_j_kg)?;
    let deposit = SolidSurfaceState {
        mass_kg: surface.state.deposit_mass_kg,
        specific_enthalpy_j_kg: (surface.deposit_enthalpy)(temperature)?,
        area_m2: surface.state.area_m2,
    };
    validate_surface(surface.binding, deposit)?;
    snapshot_surface_at_temperature_with_kinetics(
        surface.law,
        surface.binding,
        deposit,
        gas,
        catalog,
        temperature,
        SurfaceKinetics {
            accommodation_coefficient: surface.binding.accommodation_coefficient,
            mass_transfer_coefficient_m_s,
        },
        duration_s,
        controls.chamber.surface,
    )
}

struct SupportedStorage {
    state: SupportedSurfaceState,
    substrate_energy_change_j: f64,
    deposit_energy_change_j: f64,
}

fn prepare(
    surface: &GasSupportedSurface<'_>,
    snapshot: SurfaceSnapshot,
    catalog: &ThermochemicalGasModel,
    controls: GasSupportedChamberSolve,
    reflected_heat_into_gas_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(PreparedSurface, SupportedStorage), GasSurfaceError> {
    let energy = molecular_stream_energy(snapshot, reflected_heat_into_gas_j)?;
    let dm = snapshot.transferred_mass_kg;
    let next_mass = snapshot.solid.mass_kg - dm;
    let mass_residual = (next_mass - snapshot.solid.mass_kg) + dm;
    if !nonnegative(next_mass)
        || (dm != 0.0 && next_mass == snapshot.solid.mass_kg)
        || !mass_residual.is_finite()
        || mass_residual.abs() > controls.chamber.surface.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "supported deposit mass increment is unrepresentable",
        ));
    }
    // Remove/add the exchanged deposit at its INITIAL temperature, then share
    // the remaining molecular heat between the ACTUAL final deposit mass and
    // the fixed substrate. This is complete storage conservation, not a latent
    // source applied to the substrate after updating a thermally massless film.
    let heat_j = dm
        * (-snapshot.latent_heat_j_kg - snapshot.source_alignment_residual_j_kg
            + 0.5 * snapshot.gas.specific_gas_constant * snapshot.temperature_k)
        - snapshot.incoming_mass_kg * energy.difference_j_kg
        - reflected_heat_into_gas_j;
    let heat = heat_supported_surface(
        surface.law,
        SupportedSurfaceState {
            deposit_mass_kg: next_mass,
            ..surface.state
        },
        heat_j,
        controls.storage,
        surface.substrate_temperature,
        surface.deposit_enthalpy,
        || {
            if cancelled() {
                Err(SublimationError::Cancelled)
            } else {
                Ok(())
            }
        },
    )?;
    let final_temperature = heat.surface_temperature_k;
    let final_deposit_h = (surface.deposit_enthalpy)(final_temperature)?;
    let final_phase = phase_at(surface.law, final_temperature, snapshot.gas.molar_mass)?;
    let vapor_h =
        catalog.species_specific_enthalpy_j_kg(&surface.binding.species, final_temperature)?;
    compatibility(
        vapor_h,
        final_deposit_h,
        surface.binding.solid_enthalpy_offset_j_kg,
        final_phase.latent_heat_j_kg,
        final_temperature,
        controls.chamber.surface.latent_compatibility_tolerance_j_kg,
    )?;
    let deposit_energy_change_j = heat.deposit_energy_change_j
        - dm * (snapshot.solid.specific_enthalpy_j_kg + surface.binding.solid_enthalpy_offset_j_kg);
    let solid_energy_change_j = heat.substrate_energy_change_j + deposit_energy_change_j;
    let energy_residual_j = heat.energy_residual_j;
    if !deposit_energy_change_j.is_finite()
        || !solid_energy_change_j.is_finite()
        || !energy_residual_j.is_finite()
        || energy_residual_j.abs() > controls.chamber.surface.gas.energy_tolerance_j
    {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: energy_residual_j,
            tolerance_j: controls.chamber.surface.gas.energy_tolerance_j,
        });
    }
    let recession_m = if surface.state.area_m2 == 0.0 {
        0.0
    } else {
        dm / final_phase.solid_density_kg_m3 / surface.state.area_m2
    };
    if !recession_m.is_finite() {
        return Err(invalid("supported deposit recession is unrepresentable"));
    }
    Ok((
        PreparedSurface {
            solid: SolidSurfaceState {
                mass_kg: next_mass,
                specific_enthalpy_j_kg: final_deposit_h,
                area_m2: surface.state.area_m2,
            },
            surface_temperature_k: [snapshot.temperature_k, final_temperature],
            outgoing_mass_kg: snapshot.outgoing_mass_kg,
            incoming_mass_kg: snapshot.incoming_mass_kg,
            transferred_mass_kg: dm,
            outgoing_specific_energy_j_kg: snapshot.outgoing_specific_energy_j_kg,
            incoming_specific_energy_j_kg: snapshot.incoming_specific_energy_j_kg,
            energy_into_gas_j: energy.energy_into_gas_j,
            phase_relative_energy_into_gas_j: energy.phase_relative_energy_into_gas_j,
            relative_energy_into_gas_j: energy.relative_energy_into_gas_j,
            solid_energy_change_j,
            energy_residual_j,
            recession_m,
        },
        SupportedStorage {
            state: heat.state,
            substrate_energy_change_j: heat.substrate_energy_change_j,
            deposit_energy_change_j,
        },
    ))
}

/// Advance supported surfaces together using the original evolving gas owner.
///
/// Each surface retains an independent substrate chart and binding. The original
/// catalog admits its deposit species, molar mass and latent compatibility at
/// both endpoints. Only the old deposit bounds outgoing molecular turnover;
/// all incoming streams share their species' OLD gas allowance. An inventory
/// refusal supplies a conservative shorter duration for the entire interval.
/// In particular, a bare substrate at saturation is not a stationary coating.
///
/// Optional accommodation must cover exactly the active stable ids and retains
/// the existing explicit carrier/nonsticking coefficients. Its reflected heat
/// shares the same supported storage and gas endpoint, with no second heat step.
/// Source-domain, balance, inventory, work or cancellation failure publishes
/// no endpoint. An empty active selection preserves the original gas.
#[allow(clippy::too_many_arguments)]
pub fn advance_supported_chamber<'catalog>(
    surfaces: &[GasSupportedSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasSupportedChamberCandidate<'catalog>, GasSurfaceError> {
    advance_supported_chamber_impl(
        surfaces,
        accommodation,
        None,
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

/// Advance supported surfaces with optional prescribed gas-film resistances.
///
/// The coefficient map is sparse: an absent id retains the free-molecular
/// limit, while a present zero seals only that surface's phase stream. Entries
/// for unknown ids refuse the whole candidate. Finite positive values use the
/// same D/delta series resistance as the thermal sublimation owner, scaling
/// both opposing molecular streams before shared gas accounting.
#[allow(clippy::too_many_arguments)]
pub fn advance_supported_chamber_with_films<'catalog>(
    surfaces: &[GasSupportedSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    mass_transfer_coefficients_m_s: &BTreeMap<u64, f64>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasSupportedChamberCandidate<'catalog>, GasSurfaceError> {
    let ids: BTreeMap<_, _> = surfaces.iter().map(|surface| (surface.id, ())).collect();
    if mass_transfer_coefficients_m_s
        .keys()
        .any(|id| !ids.contains_key(id))
    {
        return Err(invalid("gas-film coefficient references an inactive supported surface"));
    }
    advance_supported_chamber_impl(
        surfaces,
        accommodation,
        Some(mass_transfer_coefficients_m_s),
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_supported_chamber_impl<'catalog>(
    surfaces: &[GasSupportedSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    mass_transfer_coefficients_m_s: Option<&BTreeMap<u64, f64>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedChamberSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasSupportedChamberCandidate<'catalog>, GasSurfaceError> {
    let mut storage = Vec::new();
    let (candidate, reflection) = advance_chamber_with_owner(
        surfaces,
        accommodation,
        gas,
        catalog,
        duration_s,
        controls.chamber,
        || validate_storage(controls.storage),
        |surface, gas| {
            let film = mass_transfer_coefficients_m_s
                .and_then(|coefficients| coefficients.get(&surface.id).copied());
            snapshot_with_film(surface, gas, catalog, film, duration_s, controls)
        },
        |surface, snapshot, reflected_heat, cancelled| {
            let (prepared, supported) = prepare(
                surface,
                snapshot,
                catalog,
                controls,
                reflected_heat,
                cancelled,
            )?;
            storage.push(supported);
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
        let reflected_heat_into_gas_j = reflection
            .get(&exchange.id)
            .map_or(0.0, |account| account.heat_into_gas_j);
        surfaces.push(GasSupportedSurfaceStep {
            id: exchange.id,
            state: stored.state,
            exchange,
            substrate_energy_change_j: stored.substrate_energy_change_j,
            deposit_energy_change_j: stored.deposit_energy_change_j,
            reflected_heat_into_gas_j,
        });
    }
    poll(&mut cancelled)?;
    Ok(GasSupportedChamberCandidate {
        gas,
        surfaces,
        step,
    })
}

#[cfg(test)]
mod tests;
