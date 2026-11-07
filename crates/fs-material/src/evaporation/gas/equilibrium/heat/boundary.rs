//! Complete freezing with reciprocal heat and vapor on one physical clock.
//!
//! Every trial restarts BOTH original owners from the caller's unchanged body
//! and gas. A vapor event after a full-duration heat solve is only a time
//! bracket; it is never published with that shorter vapor clock. Acceptance
//! requires the heat duration and the original film duration to be identical.
//! The original film's absolute mass allowance admits the final root rounding.
//!
//! With zero vapor flux, reciprocal heat may itself complete freezing.
//! Its represented endpoint can be placed at the exact original solid knot
//! only within the original extensive-energy AND heat-equation allowances.
//! With active vapor transport, the post-heat body must still contain liquid;
//! extinction during heat is a bracket, not permission to omit vapor transfer.
//! If that split has no admissible common-duration event, bounded work refuses
//! it. Refine the complete physical step or explicitly select the next owner.
//!
//! This locates an endpoint-bracketed event of the original heat-first Lie
//! split. It does not certify continuous-time event accuracy, monotonicity in
//! duration, or absence of multiple crossings; physical timestep refinement
//! remains necessary. No automatic solid-interface or transport-regime switch
//! occurs, and fusion is never another heat source.

use super::{
    ClosedGasParcel, EquilibriumEvaporation, EquilibriumGasHeatSolve, EvaporationError,
    EvaporationLaw, GasEvaporationError, GasEvaporationSolve, GasEvaporationStep, GasFilmHeatStep,
    GasLiquidBinding, LiquidSurfaceState, PhaseInventoryStep, ThermochemicalGasModel,
    advance_equilibrium_surface, combine_step, finite, invalid, nonnegative, poll, sensible_heat,
    smaller_step, validate_controls,
};
use crate::evaporation::gas::{
    check_interface_endpoint,
    equilibrium::boundary::{
        EquilibriumFreezingCandidate, FreezingBoundaryStep, advance_to_freezing,
        zero_flux_at_freezing,
    },
};

/// Bounded complete heat/phase trials used to locate one terminal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingHeatSolve {
    /// Maximum complete trials, in 1..=256. Source and heat inverses retain
    /// their own original bounds. Exhaustion returns a smaller-step request.
    pub maximum_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
enum Endpoint<'curve, 'catalog> {
    Transfer(EquilibriumFreezingCandidate<'curve, 'catalog>),
    HeatOnly {
        bulk: LiquidSurfaceState,
        gas: ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
    },
}

/// Original gas and total condensed body at one accepted heat/phase time.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumFreezingHeatCandidate<'curve, 'catalog> {
    inventory: EquilibriumEvaporation<'curve>,
    endpoint: Endpoint<'curve, 'catalog>,
    phase_initial_bulk: LiquidSurfaceState,
    heat: GasFilmHeatStep,
    step: GasEvaporationStep,
    phase: PhaseInventoryStep,
    boundary: FreezingBoundaryStep,
    sum_absolute_energy_residual_j: f64,
    trials: usize,
}
impl<'curve, 'catalog> EquilibriumFreezingHeatCandidate<'curve, 'catalog> {
    /// Final TOTAL condensed mass and exact original bulk enthalpy.
    #[must_use]
    pub const fn bulk(&self) -> LiquidSurfaceState {
        match &self.endpoint {
            Endpoint::Transfer(transfer) => transfer.exchange().liquid(),
            Endpoint::HeatOnly { bulk, .. } => *bulk,
        }
    }
    /// Final original finite gas at the same time as the body.
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        match &self.endpoint {
            Endpoint::Transfer(transfer) => transfer.exchange().gas(),
            Endpoint::HeatOnly { gas, .. } => gas,
        }
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        match &self.endpoint {
            Endpoint::Transfer(transfer) => transfer.exchange().catalog(),
            Endpoint::HeatOnly { catalog, .. } => catalog,
        }
    }
    #[must_use]
    pub const fn inventory(&self) -> EquilibriumEvaporation<'curve> {
        self.inventory
    }
    /// Original-to-final total storage and temperature accounts. Vapor-carried
    /// `energy_into_gas_j` excludes the separately reported reciprocal heat.
    #[must_use]
    pub const fn step(&self) -> &GasEvaporationStep {
        &self.step
    }
    #[must_use]
    pub const fn heat(&self) -> &GasFilmHeatStep {
        &self.heat
    }
    /// Original-to-final solid/liquid partition, with exact zero liquid at completion.
    #[must_use]
    pub const fn phase(&self) -> &PhaseInventoryStep {
        &self.phase
    }
    /// Phase-only owner after heat. Absent only when zero vapor flux
    /// permits heat itself to complete freezing. Its requested horizon is the
    /// accepted thermal trial, not this operation's original maximum horizon.
    #[must_use]
    pub const fn transfer(&self) -> Option<&EquilibriumFreezingCandidate<'curve, 'catalog>> {
        match &self.endpoint {
            Endpoint::Transfer(transfer) => Some(transfer),
            Endpoint::HeatOnly { .. } => None,
        }
    }
    #[must_use]
    pub const fn phase_initial_bulk(&self) -> LiquidSurfaceState {
        self.phase_initial_bulk
    }
    /// Common elapsed duration and original unconsumed duration. In heat-only
    /// completion the condensed defect is the measured heat endpoint projection.
    #[must_use]
    pub const fn boundary(&self) -> &FreezingBoundaryStep {
        &self.boundary
    }
    /// Absolute substep defects and any exact-knot correction, sharing ONE J budget.
    #[must_use]
    pub const fn sum_absolute_energy_residual_j(&self) -> f64 {
        self.sum_absolute_energy_residual_j
    }
    /// Complete thermal trials, including rejected brackets. Zero heat delegates
    /// the original analytic phase event directly and uses zero thermal trials.
    #[must_use]
    pub const fn trials(&self) -> usize {
        self.trials
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        LiquidSurfaceState,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasEvaporationStep,
        GasFilmHeatStep,
        PhaseInventoryStep,
        FreezingBoundaryStep,
    ) {
        let (bulk, gas, catalog) = match self.endpoint {
            Endpoint::Transfer(transfer) => {
                let (exchange, _, _) = transfer.into_parts();
                let (bulk, gas, catalog, _) = exchange.into_parts();
                (bulk, gas, catalog)
            }
            Endpoint::HeatOnly { bulk, gas, catalog } => (bulk, gas, catalog),
        };
        (
            bulk,
            gas,
            catalog,
            self.step,
            self.heat,
            self.phase,
            self.boundary,
        )
    }
}

fn phase_defect(transfer: &EquilibriumFreezingCandidate<'_, '_>) -> f64 {
    if transfer.boundary().freezing_completed {
        transfer.boundary().condensed_energy_residual_j.abs()
            + transfer
                .exchange()
                .gas_exchange()
                .step()
                .energy_residual_j
                .abs()
    } else {
        // Keep the unchanged scalar owner's admission on a nonterminal step.
        transfer.exchange().step().energy_residual_j.abs()
    }
}

/// One complete original heat/phase trial, without a nested event-time search.
/// A bracket carries no endpoint that a coupled caller could publish by mistake.
#[derive(Clone, Debug, PartialEq)]
pub(in crate::evaporation::gas) enum FreezingHeatTrial<'curve, 'catalog> {
    Accepted(EquilibriumFreezingHeatCandidate<'curve, 'catalog>),
    Bracket { residual_energy_j: f64 },
}

struct TrialStart {
    original_step: GasEvaporationStep,
    native_heat: super::GasFilmHeatSolve,
}

#[allow(clippy::too_many_arguments)]
fn prepare_trial_start<L: EvaporationLaw + ?Sized>(
    inventory: EquilibriumEvaporation<'_>,
    law: &L,
    binding: &GasLiquidBinding,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    phase_controls: GasEvaporationSolve,
    heat_controls: EquilibriumGasHeatSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<TrialStart, GasEvaporationError> {
    poll(&mut cancelled)?;
    validate_controls(phase_controls, duration_s, mass_transfer_coefficient_m_s)?;
    let mut native_heat =
        heat_controls.original_controls(inventory, bulk.specific_enthalpy_j_kg)?;
    let original = advance_equilibrium_surface(
        inventory,
        law,
        binding,
        bulk,
        gas,
        catalog,
        mass_transfer_coefficient_m_s,
        0.0,
        phase_controls,
        &mut cancelled,
    )?;
    // Preserve the exact source knot as the existing heat solver's far bound.
    // Its original bisection can then accept an exactly represented h_s.
    native_heat.pure_liquid_enthalpy_range_j_kg[0] =
        native_heat.pure_liquid_enthalpy_range_j_kg[0].max(inventory.solid_at_melting_h_j_kg);
    Ok(TrialStart {
        original_step: original.exchange().step().clone(),
        native_heat,
    })
}

fn zero_heat_candidate<'curve, 'catalog>(
    inventory: EquilibriumEvaporation<'curve>,
    bulk: LiquidSurfaceState,
    transfer: EquilibriumFreezingCandidate<'curve, 'catalog>,
) -> EquilibriumFreezingHeatCandidate<'curve, 'catalog> {
    EquilibriumFreezingHeatCandidate {
        inventory,
        phase_initial_bulk: bulk,
        heat: GasFilmHeatStep::default(),
        step: transfer.exchange().step().clone(),
        phase: *transfer.phase(),
        boundary: *transfer.boundary(),
        sum_absolute_energy_residual_j: phase_defect(&transfer),
        endpoint: Endpoint::Transfer(transfer),
        trials: 0,
    }
}

#[allow(clippy::too_many_arguments)]
fn prepared_trial<'curve, 'catalog, L: EvaporationLaw + ?Sized>(
    inventory: EquilibriumEvaporation<'curve>,
    law: &L,
    binding: &GasLiquidBinding,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    duration_s: f64,
    phase_controls: GasEvaporationSolve,
    start: &TrialStart,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FreezingHeatTrial<'curve, 'catalog>, GasEvaporationError> {
    poll(&mut cancelled)?;
    let original_step = &start.original_step;
    let native_heat = start.native_heat;
    let initial_t = original_step.surface_temperature_k[0];
    let boundary_h = inventory.solid_at_melting_h_j_kg;
    let temperature = |h| {
        inventory
            .curve()
            .temperature_at_specific_enthalpy(h)
            .map_err(|_| {
                EvaporationError::Invalid("original bulk fusion caloric chart refused heat trial")
            })
    };
    let heat_trial = sensible_heat(
        bulk,
        gas,
        initial_t,
        &temperature,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        phase_controls,
        native_heat,
        &mut cancelled,
    );
    let (heated_bulk, heated_gas, heat) = match heat_trial {
        Ok(candidate) => candidate,
        Err(GasEvaporationError::HeatOutsideBracket {
            residual_temperature_k,
        }) if native_heat.pure_liquid_enthalpy_range_j_kg[0] == boundary_h
            && gas.state().state().temperature < initial_t
            && phase_controls.gas.gas.maximum_temperature_k >= initial_t
            && residual_temperature_k[1] < 0.0 =>
        {
            // With the gas-to-initial-body-temperature end available, the
            // only cooling bracket cut off here is the exact solid knot.
            // Its original heat-equation residual locates a time bracket;
            // no gas or condensed endpoint from this refused trial is used.
            let residual = finite(
                residual_temperature_k[1]
                    * (heat_transfer_coefficient_w_m2_k * bulk.area_m2 * duration_s),
            )?;
            return Ok(FreezingHeatTrial::Bracket {
                residual_energy_j: residual,
            });
        }
        Err(error) => return Err(error),
    };
    let heat_gap = finite(heated_bulk.mass_kg * (heated_bulk.specific_enthalpy_j_kg - boundary_h))?;

    // Disabled transport or exact saturation permits the heat endpoint
    // itself to terminate the interface. Never suppress an active stream.
    if heat.heat_into_liquid_j < 0.0 {
        let exact_q = finite(bulk.mass_kg * (boundary_h - bulk.specific_enthalpy_j_kg))?;
        let correction = finite(exact_q - heat.liquid_enthalpy_change_j)?;
        let absolute_defect = finite(heat.energy_residual_j.abs() + correction.abs())?;
        if absolute_defect <= phase_controls.gas.energy_tolerance_j
            && (native_heat.pure_liquid_enthalpy_range_j_kg[0]
                ..=native_heat.pure_liquid_enthalpy_range_j_kg[1])
                .contains(&boundary_h)
            && zero_flux_at_freezing(
                inventory,
                law,
                binding,
                &heated_gas,
                catalog,
                mass_transfer_coefficient_m_s,
                phase_controls,
            )?
        {
            let caloric = inventory.caloric_with_completion(boundary_h, true)?;
            let exposure = finite(heat_transfer_coefficient_w_m2_k * bulk.area_m2 * duration_s)?;
            let residual_t = finite(
                (heated_gas.state().state().temperature - caloric.temperature_k)
                    - exact_q / exposure,
            )?;
            if exposure > 0.0
                && residual_t.abs() <= native_heat.temperature_tolerance_k
                && heated_gas.state().state().temperature <= caloric.temperature_k
            {
                let final_bulk = LiquidSurfaceState {
                    specific_enthalpy_j_kg: boundary_h,
                    ..bulk
                };
                let mut heat = heat;
                heat.heat_into_liquid_j = exact_q;
                heat.liquid_enthalpy_change_j = exact_q;
                heat.energy_residual_j = finite(exact_q + heat.gas_internal_energy_change_j)?;
                heat.temperature_equation_residual_k = residual_t;
                let mut empty_phase = original_step.clone();
                empty_phase.duration_s = duration_s;
                empty_phase.surface_temperature_k[1] = caloric.temperature_k;
                empty_phase.gas_temperature_k[1] = heated_gas.state().state().temperature;
                empty_phase.vapor_partial_pressure_pa[1] = check_interface_endpoint(
                    law,
                    binding,
                    caloric.temperature_k,
                    heated_gas.state().state().temperature,
                    0.0,
                    &heated_gas,
                    mass_transfer_coefficient_m_s,
                    duration_s,
                    phase_controls,
                )?;
                let absolute_defect = absolute_defect.max(heat.energy_residual_j.abs());
                let (step, absolute_defect) = combine_step(
                    original_step,
                    empty_phase,
                    heat,
                    absolute_defect,
                    phase_controls,
                )?;
                let phase = inventory.partition_step_with_completion(
                    bulk,
                    final_bulk,
                    0.0,
                    phase_controls.gas.mass_tolerance_kg,
                    true,
                )?;
                poll(&mut cancelled)?;
                return Ok(FreezingHeatTrial::Accepted(
                    EquilibriumFreezingHeatCandidate {
                        inventory,
                        endpoint: Endpoint::HeatOnly {
                            bulk: final_bulk,
                            gas: heated_gas,
                            catalog,
                        },
                        phase_initial_bulk: final_bulk,
                        heat,
                        step,
                        phase,
                        boundary: FreezingBoundaryStep {
                            elapsed_duration_s: duration_s,
                            remaining_duration_s: 0.0,
                            freezing_completed: true,
                            condensed_energy_residual_j: correction,
                        },
                        sum_absolute_energy_residual_j: absolute_defect,
                        trials: 1,
                    },
                ));
            }
        }
    }

    let residual = if heat_gap <= 0.0 {
        // The active phase owner cannot start from a completed solid. This
        // is only an upper time bracket; no vapor endpoint is fabricated.
        heat_gap
    } else {
        let transfer = advance_to_freezing(
            inventory,
            law,
            binding,
            heated_bulk,
            &heated_gas,
            catalog,
            mass_transfer_coefficient_m_s,
            duration_s,
            phase_controls,
            &mut cancelled,
        )?;
        let event = *transfer.boundary();
        if event.elapsed_duration_s == duration_s {
            let absolute_defect = finite(heat.energy_residual_j.abs() + phase_defect(&transfer))?;
            let (step, absolute_defect) = combine_step(
                original_step,
                transfer.exchange().step().clone(),
                heat,
                absolute_defect,
                phase_controls,
            )?;
            let phase = inventory.partition_step_with_completion(
                bulk,
                transfer.exchange().liquid(),
                step.transferred_mass_kg,
                phase_controls.gas.mass_tolerance_kg,
                event.freezing_completed,
            )?;
            poll(&mut cancelled)?;
            return Ok(FreezingHeatTrial::Accepted(
                EquilibriumFreezingHeatCandidate {
                    inventory,
                    endpoint: Endpoint::Transfer(transfer),
                    phase_initial_bulk: heated_bulk,
                    heat,
                    step,
                    phase,
                    boundary: FreezingBoundaryStep {
                        remaining_duration_s: 0.0,
                        ..event
                    },
                    sum_absolute_energy_residual_j: absolute_defect,
                    trials: 1,
                },
            ));
        } else {
            // Positive post-heat seed, but phase would finish too early.
            // The frozen-rate extensive residual has this sign and root.
            // Scaling by dt bounds it even for an enormous upper horizon.
            finite(heat_gap * ((event.elapsed_duration_s - duration_s) / duration_s))?
        }
    };
    Ok(FreezingHeatTrial::Bracket {
        residual_energy_j: residual,
    })
}

/// One common-duration trial for callers that already own a complete coupled
/// time search. Even zero heat returns a bracket when vapor would finish early:
/// preceding native heat must be replayed, never relabeled with that short clock.
#[allow(clippy::too_many_arguments)]
pub(in crate::evaporation::gas) fn single_trial_to_freezing_with_heat<
    'curve,
    'catalog,
    L: EvaporationLaw + ?Sized,
>(
    inventory: EquilibriumEvaporation<'curve>,
    law: &L,
    binding: &GasLiquidBinding,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    duration_s: f64,
    phase_controls: GasEvaporationSolve,
    heat_controls: EquilibriumGasHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<FreezingHeatTrial<'curve, 'catalog>, GasEvaporationError> {
    if !nonnegative(heat_transfer_coefficient_w_m2_k) {
        return Err(invalid(
            "finite nonnegative sensible heat coefficient required",
        ));
    }
    if heat_transfer_coefficient_w_m2_k == 0.0 {
        let transfer = advance_to_freezing(
            inventory,
            law,
            binding,
            bulk,
            gas,
            catalog,
            mass_transfer_coefficient_m_s,
            duration_s,
            phase_controls,
            cancelled,
        )?;
        if transfer.boundary().elapsed_duration_s != duration_s {
            let gap = finite(
                bulk.mass_kg * (bulk.specific_enthalpy_j_kg - inventory.solid_at_melting_h_j_kg),
            )?;
            return Ok(FreezingHeatTrial::Bracket {
                residual_energy_j: finite(
                    gap * ((transfer.boundary().elapsed_duration_s - duration_s) / duration_s),
                )?,
            });
        }
        return Ok(FreezingHeatTrial::Accepted(zero_heat_candidate(
            inventory, bulk, transfer,
        )));
    }
    let mut cancelled = cancelled;
    let start = prepare_trial_start(
        inventory,
        law,
        binding,
        bulk,
        gas,
        catalog,
        mass_transfer_coefficient_m_s,
        duration_s,
        phase_controls,
        heat_controls,
        &mut cancelled,
    )?;
    prepared_trial(
        inventory,
        law,
        binding,
        bulk,
        gas,
        catalog,
        mass_transfer_coefficient_m_s,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        phase_controls,
        &start,
        cancelled,
    )
}

/// Advance at most the requested duration, locating a shared heat/phase event.
///
/// Zero heat delegates exactly to the original analytic phase boundary, without
/// reading unused heat/event controls. Otherwise all original source, donor,
/// temperature, saturation and gas budgets remain active at every evaluated
/// endpoint. Invalid controls, refused original owners, bounded nonconvergence,
/// or cancellation publish neither caller-owned input.
#[allow(clippy::too_many_arguments)]
pub fn advance_to_freezing_with_heat<'curve, 'catalog, L: EvaporationLaw + ?Sized>(
    inventory: EquilibriumEvaporation<'curve>,
    law: &L,
    binding: &GasLiquidBinding,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    maximum_duration_s: f64,
    phase_controls: GasEvaporationSolve,
    heat_controls: EquilibriumGasHeatSolve,
    event_controls: FreezingHeatSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumFreezingHeatCandidate<'curve, 'catalog>, GasEvaporationError> {
    if !nonnegative(heat_transfer_coefficient_w_m2_k) {
        return Err(invalid(
            "finite nonnegative sensible heat coefficient required",
        ));
    }
    if heat_transfer_coefficient_w_m2_k == 0.0 {
        let transfer = advance_to_freezing(
            inventory,
            law,
            binding,
            bulk,
            gas,
            catalog,
            mass_transfer_coefficient_m_s,
            maximum_duration_s,
            phase_controls,
            cancelled,
        )?;
        return Ok(zero_heat_candidate(inventory, bulk, transfer));
    }
    if !(1..=256).contains(&event_controls.maximum_iterations) {
        return Err(invalid(
            "freezing event requires 1..=256 complete thermal trials",
        ));
    }
    let start = prepare_trial_start(
        inventory,
        law,
        binding,
        bulk,
        gas,
        catalog,
        mass_transfer_coefficient_m_s,
        maximum_duration_s,
        phase_controls,
        heat_controls,
        &mut cancelled,
    )?;
    let boundary_h = inventory.solid_at_melting_h_j_kg;
    let mut near_time = 0.0;
    let mut near_residual = finite(bulk.mass_kg * (bulk.specific_enthalpy_j_kg - boundary_h))?;
    let mut far_time = maximum_duration_s;
    let mut far_residual = 0.0;
    let mut next_time = |time: f64, residual: f64| {
        if residual > 0.0 {
            near_time = time;
            near_residual = residual;
        } else {
            far_time = time;
            far_residual = residual;
        }
        let width = far_time - near_time;
        let midpoint = near_time + 0.5 * width;
        if midpoint == near_time || midpoint == far_time {
            return None;
        }
        let denominator = near_residual - far_residual;
        let weight = near_residual / denominator;
        Some(
            if denominator.is_finite() && (0.125..=0.875).contains(&weight) {
                near_time + weight * width
            } else {
                midpoint
            },
        )
    };
    let mut duration_s = maximum_duration_s;
    for trials in 1..=event_controls.maximum_iterations {
        let trial = prepared_trial(
            inventory,
            law,
            binding,
            bulk,
            gas,
            catalog,
            mass_transfer_coefficient_m_s,
            heat_transfer_coefficient_w_m2_k,
            duration_s,
            phase_controls,
            &start,
            &mut cancelled,
        )?;
        let residual = match trial {
            FreezingHeatTrial::Accepted(mut candidate) => {
                if candidate.boundary.freezing_completed || duration_s == maximum_duration_s {
                    candidate.boundary.remaining_duration_s = maximum_duration_s - duration_s;
                    candidate.trials = trials;
                    return Ok(candidate);
                }
                let final_bulk = candidate.bulk();
                finite(final_bulk.mass_kg * (final_bulk.specific_enthalpy_j_kg - boundary_h))?
            }
            FreezingHeatTrial::Bracket { residual_energy_j } => residual_energy_j,
        };
        match next_time(duration_s, residual) {
            Some(next) => duration_s = next,
            None => break,
        }
    }
    Err(smaller_step(
        maximum_duration_s,
        "common-duration heat/phase freezing event did not converge",
    ))
}

#[cfg(test)]
mod tests;
