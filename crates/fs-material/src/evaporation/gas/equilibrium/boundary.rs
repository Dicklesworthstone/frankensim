//! Stop liquid-interface exchange at conservative, complete freezing.
//!
//! For the original frozen-rate film, `m1=m0-dm` and
//! `m1*(h1-h0)=-dm*g`, where g is the original carried vapor enthalpy minus
//! gas-aligned BULK h. At the original solid fusion knot h_s, the signed event
//! mass is `dm=m0*(h0-h_s)/(g+h0-h_s)`. This constant-work event calculation
//! also admits cooling by condensation when upstream vapor enthalpy is low.
//!
//! The shared mass/energy owner clips the active duration before donor checks,
//! selects the represented event mass within the original film mass budget,
//! and tests the actual extensive defect at the EXACT source knot. The original
//! gas is reconstructed once. No residual liquid, latent inventory, temperature
//! epsilon, or second fusion source is silently discarded.
//!
//! A completed step returns unconsumed time and an all-solid body for an
//! explicitly selected next owner. Stefan-film evaporation and free-molecular
//! sublimation have different transport assumptions; this operation never
//! switches between them. Exposed area, geometry and gas volume retain their
//! original owners. Event times belong to the frozen-rate discretization and
//! require physical timestep refinement for continuous-time accuracy.
//!
//! This is a phase-only interval. If heat is also split over the same interval,
//! a shorter event requires retrying the WHOLE heat/phase step. A full-duration
//! heat endpoint cannot simply be assigned the shorter event clock.

use super::{
    ClosedGasParcel, EquilibriumEvaporation, EvaporationLaw, GasEvaporationCandidate,
    GasEvaporationError, GasEvaporationSolve, GasLiquidBinding, LiquidSurfaceState,
    PhaseInventoryStep, ThermochemicalGasModel, finish_surface, invalid, nonnegative,
};
use crate::evaporation::gas::{admit_liquid_source, poll, prepare_surface_with_boundary};
use fs_thermal::evaporation::{GasFilm, surface_flux};

/// Admit a terminal interface with the original source/film owners. Exact
/// saturation can turn off vapor exchange even when its coefficient is nonzero.
/// This does not invent a liquid seed or admit a finite solid-interface step.
pub(in crate::evaporation::gas) fn zero_flux_at_freezing<L: EvaporationLaw + ?Sized>(
    inventory: EquilibriumEvaporation<'_>,
    law: &L,
    binding: &GasLiquidBinding,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    controls: GasEvaporationSolve,
) -> Result<bool, GasEvaporationError> {
    let caloric = inventory.caloric_with_completion(inventory.solid_at_melting_h_j_kg, true)?;
    let source = admit_liquid_source(
        law,
        binding,
        catalog,
        caloric.temperature_k,
        caloric.liquid_enthalpy_j_kg,
        controls,
    )?;
    caloric.admit_phase(source.phase, true)?;
    let flux = surface_flux(
        law,
        caloric.temperature_k,
        GasFilm {
            total_pressure_pa: gas.state().state().pressure,
            vapor_partial_pressure_pa: gas.species_partial_pressure_pa(&binding.species)?,
            temperature_k: gas.state().state().temperature,
            mass_transfer_coefficient_m_s,
        },
    )?;
    Ok(flux.mass_flux_kg_m2_s == 0.0)
}

/// Physical time consumed by the existing film and the terminal phase decision.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingBoundaryStep {
    /// Accepted interface duration [s]; also the reciprocal heat duration when
    /// this receipt is returned by the combined heat/phase owner.
    pub elapsed_duration_s: f64,
    /// Requested duration not consumed before completion [s].
    pub remaining_duration_s: f64,
    /// The body is exactly at the original all-solid fusion endpoint.
    pub freezing_completed: bool,
    /// Represented condensed storage minus the original phase/heat target [J],
    /// as documented by the consuming owner. Included in its energy residual,
    /// never an added heat source.
    pub condensed_energy_residual_j: f64,
}

/// Original gas, bulk enthalpy and phase inventory at the SAME accepted time.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumFreezingCandidate<'curve, 'catalog> {
    inventory: EquilibriumEvaporation<'curve>,
    exchange: GasEvaporationCandidate<'catalog>,
    phase: PhaseInventoryStep,
    boundary: FreezingBoundaryStep,
}
impl<'curve, 'catalog> EquilibriumFreezingCandidate<'curve, 'catalog> {
    /// Original gas and TOTAL condensed body; `liquid()` is the bulk store here.
    /// An all-solid result is terminal for the liquid-interface owner.
    #[must_use]
    pub const fn exchange(&self) -> &GasEvaporationCandidate<'catalog> {
        &self.exchange
    }
    /// Original-to-final constituent masses, including exact zero liquid at completion.
    #[must_use]
    pub const fn phase(&self) -> &PhaseInventoryStep {
        &self.phase
    }
    /// Accepted event time, unconsumed time and measured terminal storage defect.
    #[must_use]
    pub const fn boundary(&self) -> &FreezingBoundaryStep {
        &self.boundary
    }
    /// Unmodified original caloric source. Its default active-interface state
    /// method continues to refuse a completed all-solid endpoint.
    #[must_use]
    pub const fn inventory(&self) -> EquilibriumEvaporation<'curve> {
        self.inventory
    }
    /// Consume the mutually consistent endpoint and its handoff information.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasEvaporationCandidate<'catalog>,
        PhaseInventoryStep,
        FreezingBoundaryStep,
    ) {
        (self.exchange, self.phase, self.boundary)
    }
}

/// Advance at most `maximum_duration_s`, stopping at complete freezing if reached.
///
/// Starts with an existing liquid seed. If no event occurs, the result agrees
/// with the existing scalar step and consumes the entire duration. All original
/// donor, source, temperature, saturation and gas reconstruction gates remain
/// active. In particular, the caller's inventory fraction must admit the mass
/// required by the event; this function does not relax it to force completion.
/// Rounding of event mass/time must fit the same absolute mass/energy budgets.
/// No caller-owned state is modified on failure or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn advance_to_freezing<'curve, 'catalog, L: EvaporationLaw + ?Sized>(
    inventory: EquilibriumEvaporation<'curve>,
    law: &L,
    binding: &GasLiquidBinding,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    maximum_duration_s: f64,
    controls: GasEvaporationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumFreezingCandidate<'curve, 'catalog>, GasEvaporationError> {
    let prepared = prepare_surface_with_boundary(
        law,
        binding,
        bulk,
        gas,
        catalog,
        |h| inventory.caloric_with_completion(h, true),
        mass_transfer_coefficient_m_s,
        maximum_duration_s,
        controls,
        Ok,
        Some(inventory.solid_at_melting_h_j_kg),
        &mut cancelled,
    )?;
    let elapsed_duration_s = prepared.duration_s;
    let remaining_duration_s = maximum_duration_s - elapsed_duration_s;
    if !nonnegative(remaining_duration_s) {
        return Err(invalid(
            "terminal phase event cannot represent its remaining duration",
        ));
    }
    let freezing_completed = prepared.terminal_enthalpy_j_kg.is_some();
    let phase = inventory.partition_step_with_completion(
        bulk,
        prepared.liquid,
        prepared.dm,
        controls.gas.mass_tolerance_kg,
        freezing_completed,
    )?;
    let boundary = FreezingBoundaryStep {
        elapsed_duration_s,
        remaining_duration_s,
        freezing_completed,
        condensed_energy_residual_j: prepared.liquid_energy_residual_j,
    };
    let exchange = finish_surface(
        prepared,
        law,
        binding,
        gas,
        catalog,
        mass_transfer_coefficient_m_s,
        elapsed_duration_s,
        controls,
        &mut cancelled,
    )?;
    if freezing_completed {
        let absolute_residual = boundary.condensed_energy_residual_j.abs()
            + exchange.gas_exchange().step().energy_residual_j.abs();
        if !absolute_residual.is_finite() || absolute_residual > controls.gas.energy_tolerance_j {
            return Err(GasEvaporationError::EnergyBalance {
                residual_j: absolute_residual,
                tolerance_j: controls.gas.energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
    }
    Ok(EquilibriumFreezingCandidate {
        inventory,
        exchange,
        phase,
        boundary,
    })
}

#[cfg(test)]
mod tests;
