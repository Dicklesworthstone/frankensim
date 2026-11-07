//! Reciprocal gas-film heat and vapor transfer through a melting body.
//!
//! The existing implicit sensible-heat owner solves on the original BULK
//! enthalpy chart, including its isothermal fusion plateau. Heat acts first;
//! the original equilibrium surface then transfers vapor over the same physical
//! interval, using the actual liquid inventory after heat. This is a first-order
//! Lie split, not a simultaneous heat/phase solve. Refine the complete interval
//! to assess time error. Conservation alone does not certify splitting accuracy.
//!
//! The sensible coefficient must exclude vapor-carried sensible and latent
//! energy. Fusion remains in bulk storage and is never another heat source.
//! Original, post-heat and final accepted bodies all require a liquid seed.
//! Trial heat brackets may use the original solid branch without publishing a
//! solid endpoint as a liquid surface. Exposed area and geometry remain explicit.

/// Complete freezing with both original substeps on the same accepted clock.
pub mod boundary;

use super::{
    ClosedGasParcel, EquilibriumEvaporation, EquilibriumSurfaceCandidate, EvaporationLaw,
    GasEvaporationError, GasEvaporationSolve, GasLiquidBinding, LiquidSurfaceState,
    PhaseInventoryStep, ThermochemicalGasModel, advance_equilibrium_surface, invalid,
};
use crate::evaporation::gas::{
    GasEvaporationStep,
    heat::{GasFilmHeatSolve, GasFilmHeatStep, sensible_heat},
    nonnegative, poll, positive, smaller_step, validate_controls,
};
use fs_thermal::evaporation::EvaporationError;

/// Original bulk-enthalpy domain and bounded work for reciprocal sensible heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumGasHeatSolve {
    /// Admitted portion of the original body chart [J/kg], containing the
    /// initial bulk h. It may include fusion and solid trial brackets.
    pub bulk_enthalpy_range_j_kg: [f64; 2],
    /// Absolute backward-Euler heat-equation residual [K]. Conservation uses
    /// the single joule allowance in `GasEvaporationSolve::gas` instead.
    pub temperature_tolerance_k: f64,
    /// Maximum scalar bisection trials, in 1..=256.
    pub maximum_iterations: usize,
}
impl EquilibriumGasHeatSolve {
    fn original_controls(
        self,
        inventory: EquilibriumEvaporation<'_>,
        initial_h: f64,
    ) -> Result<GasFilmHeatSolve, GasEvaporationError> {
        let [lo, hi] = self.bulk_enthalpy_range_j_kg;
        let knots = inventory.curve().knots();
        if !lo.is_finite()
            || !hi.is_finite()
            || lo > hi
            || !(lo..=hi).contains(&initial_h)
            || lo < knots[0].specific_enthalpy_j_kg
            || hi > knots[knots.len() - 1].specific_enthalpy_j_kg
            || !positive(self.temperature_tolerance_k)
            || !(1..=256).contains(&self.maximum_iterations)
        {
            return Err(invalid(
                "heat needs an original bulk-enthalpy interval containing the initial state and positive bounded solve controls",
            ));
        }
        // Only the existing scalar heat solver consumes this legacy field
        // name. Its h-to-T callback below is the original BULK caloric chart.
        Ok(GasFilmHeatSolve {
            pure_liquid_enthalpy_range_j_kg: self.bulk_enthalpy_range_j_kg,
            temperature_tolerance_k: self.temperature_tolerance_k,
            maximum_iterations: self.maximum_iterations,
        })
    }
}

/// One mutually consistent heat, vapor and condensed-fusion endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasHeatCandidate<'curve, 'catalog> {
    transfer: EquilibriumSurfaceCandidate<'curve, 'catalog>,
    phase_initial_bulk: LiquidSurfaceState,
    heat: GasFilmHeatStep,
    step: GasEvaporationStep,
    phase: PhaseInventoryStep,
    sum_absolute_energy_residual_j: f64,
}
impl<'curve, 'catalog> EquilibriumGasHeatCandidate<'curve, 'catalog> {
    /// Final TOTAL condensed mass and original bulk enthalpy.
    #[must_use]
    pub const fn bulk(&self) -> LiquidSurfaceState {
        self.transfer.exchange().liquid()
    }
    /// Original finite gas after both substeps.
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.transfer.exchange().gas()
    }
    /// Retained original species and caloric catalog.
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.transfer.exchange().catalog()
    }
    /// Original-to-final temperatures and TOTAL storage changes. The legacy
    /// `liquid_enthalpy_change_j` field describes the whole condensed body.
    /// `energy_into_gas_j` is the vapor stream alone; subtract sensible heat
    /// into the body to obtain the net energy delivered to the gas.
    #[must_use]
    pub const fn step(&self) -> &GasEvaporationStep {
        &self.step
    }
    /// Existing heat receipt. Its `liquid` fields describe the TOTAL condensed
    /// body here; represented heat can change phase fraction with no change in T.
    #[must_use]
    pub const fn heat(&self) -> &GasFilmHeatStep {
        &self.heat
    }
    /// Fusion reconstructed from the original body to the complete endpoint.
    #[must_use]
    pub const fn phase(&self) -> &PhaseInventoryStep {
        &self.phase
    }
    /// Vapor-transfer substep and its own post-heat-to-final fusion inventory.
    #[must_use]
    pub const fn transfer(&self) -> &EquilibriumSurfaceCandidate<'curve, 'catalog> {
        &self.transfer
    }
    /// Exact bulk state passed from reciprocal heat to the vapor-transfer owner.
    #[must_use]
    pub const fn phase_initial_bulk(&self) -> LiquidSurfaceState {
        self.phase_initial_bulk
    }
    /// Heat and vapor substep defects share one absolute joule allowance;
    /// opposing numerical defects cannot cancel to pass that admission.
    #[must_use]
    pub const fn sum_absolute_energy_residual_j(&self) -> f64 {
        self.sum_absolute_energy_residual_j
    }
    /// Consume the original endpoints and complete heat/fusion accounts together.
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
    ) {
        let (exchange, _) = self.transfer.into_parts();
        let (bulk, gas, catalog, _) = exchange.into_parts();
        (bulk, gas, catalog, self.step, self.heat, self.phase)
    }
}

fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "combined bulk heat and vapor account is unrepresentable",
        ))
    }
}

// One receipt owner for the strict and terminal heat/phase facades. The caller
// supplies the absolute substep allowance, including any exact-knot projection.
fn combine_step(
    original: &GasEvaporationStep,
    mut step: GasEvaporationStep,
    heat: GasFilmHeatStep,
    sum_absolute_energy_residual_j: f64,
    controls: GasEvaporationSolve,
) -> Result<(GasEvaporationStep, f64), GasEvaporationError> {
    step.surface_temperature_k[0] = original.surface_temperature_k[0];
    step.gas_temperature_k[0] = original.gas_temperature_k[0];
    step.vapor_partial_pressure_pa[0] = original.vapor_partial_pressure_pa[0];
    step.liquid_enthalpy_change_j =
        finite(step.liquid_enthalpy_change_j + heat.liquid_enthalpy_change_j)?;
    step.gas_internal_energy_change_j =
        finite(step.gas_internal_energy_change_j + heat.gas_internal_energy_change_j)?;
    let sum_absolute_energy_residual_j = finite(sum_absolute_energy_residual_j)?;
    step.energy_residual_j = finite(step.energy_residual_j + heat.energy_residual_j)?;
    if sum_absolute_energy_residual_j > controls.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: sum_absolute_energy_residual_j,
            tolerance_j: controls.gas.energy_tolerance_j,
        });
    }
    if (step.surface_temperature_k[1] - step.surface_temperature_k[0]).abs()
        > controls.maximum_temperature_increment_k
        || (step.gas_temperature_k[1] - step.gas_temperature_k[0]).abs()
            > controls.maximum_temperature_increment_k
    {
        return Err(smaller_step(
            step.duration_s,
            "combined heat and phase temperature departure",
        ));
    }
    Ok((step, sum_absolute_energy_residual_j))
}

/// Implicit reciprocal gas heat, then original fusion-aware vapor transfer.
///
/// Zero sensible coefficient delegates exactly to `advance_equilibrium_surface`
/// without inspecting unused heat controls or invoking a heat solve. Otherwise,
/// the original phase/source admission precedes heat and both accepted substeps
/// preserve a liquid seed. All failures leave both caller-owned inputs intact.
#[allow(clippy::too_many_arguments)]
pub fn advance_equilibrium_surface_with_heat<'curve, 'catalog, L: EvaporationLaw + ?Sized>(
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
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumGasHeatCandidate<'curve, 'catalog>, GasEvaporationError> {
    if !nonnegative(heat_transfer_coefficient_w_m2_k) {
        return Err(invalid(
            "finite nonnegative sensible heat coefficient required",
        ));
    }
    if heat_transfer_coefficient_w_m2_k == 0.0 {
        let transfer = advance_equilibrium_surface(
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
        let step = transfer.exchange().step().clone();
        return Ok(EquilibriumGasHeatCandidate {
            phase: *transfer.phase(),
            transfer,
            phase_initial_bulk: bulk,
            heat: GasFilmHeatStep::default(),
            sum_absolute_energy_residual_j: step.energy_residual_j.abs(),
            step,
        });
    }
    poll(&mut cancelled)?;
    validate_controls(phase_controls, duration_s, mass_transfer_coefficient_m_s)?;
    let heat_controls = heat_controls.original_controls(inventory, bulk.specific_enthalpy_j_kg)?;
    // A zero-duration original phase call admits species, constituent calorics,
    // density, vapor law and the initial liquid seed before heat can change them.
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
    let initial_t = original.exchange().step().surface_temperature_k[0];
    let temperature = |h| {
        inventory
            .curve()
            .temperature_at_specific_enthalpy(h)
            .map_err(|_| {
                EvaporationError::Invalid("original bulk fusion caloric chart refused heat trial")
            })
    };
    let (heated_bulk, heated_gas, heat) = sensible_heat(
        bulk,
        gas,
        initial_t,
        &temperature,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        phase_controls,
        heat_controls,
        &mut cancelled,
    )?;
    let transfer = advance_equilibrium_surface(
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
    let (step, sum_absolute_energy_residual_j) = combine_step(
        original.exchange().step(),
        transfer.exchange().step().clone(),
        heat,
        transfer.exchange().step().energy_residual_j.abs() + heat.energy_residual_j.abs(),
        phase_controls,
    )?;
    let phase = inventory.partition_step(
        bulk,
        transfer.exchange().liquid(),
        step.transferred_mass_kg,
        phase_controls.gas.mass_tolerance_kg,
    )?;
    poll(&mut cancelled)?;
    Ok(EquilibriumGasHeatCandidate {
        transfer,
        phase_initial_bulk: heated_bulk,
        heat,
        step,
        phase,
        sum_absolute_energy_residual_j,
    })
}

#[cfg(test)]
mod tests;
