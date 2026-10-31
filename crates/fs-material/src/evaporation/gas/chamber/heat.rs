//! Simultaneous reciprocal sensible heat, followed by joint phase exchange.
//!
//! All active films see one finite gas temperature. At fixed masses the heat
//! stage solves backward Euler: Q_i = K_i (T_g1 - T_i1), Q_i = m_i Delta h_i,
//! and m_g Delta u_g + sum Q_i = 0, with K_i = htc_i A_i dt. Each original
//! liquid h-to-T chart and the original NASA gas increment remain authoritative.
//! Scalar trials never reconstruct gas transport or replace a caloric source.
//!
//! Heat and phase act over the SAME physical interval in a first-order Lie
//! split. Sensible coefficients exclude all vapor-carried sensible and latent
//! heat. Other films may move the gas through a film's initial temperature;
//! only endpoint flux direction, source domains and global extrema are imposed.
//! The final gas is reconstructed once for heat, then once for joint phase.

mod solve;

use super::super::{GasFilmHeatSolve, nonnegative, poll, smaller_step};
use super::{
    ClosedGasParcel, GasEvaporationError, GasLiquidChamberCandidate, GasLiquidChamberSolve,
    GasLiquidChamberStep, GasLiquidChamberSurface, GasLiquidChamberSurfaceStep, LiquidSurfaceState,
    ThermochemicalGasModel, advance_chamber_with_gas, advance_chamber_with_gas_mass_selection,
    finite, invalid,
};

/// Original physical surface plus its separately prescribed sensible heat law.
#[derive(Clone, Copy)]
pub struct GasLiquidChamberHeatSurface<'a> {
    pub surface: GasLiquidChamberSurface<'a>,
    pub heat_transfer_coefficient_w_m2_k: f64,
    pub heat_controls: GasFilmHeatSolve,
}

/// Bounded outer gas-temperature work. Each film's original heat controls bound
/// its local scalar response. Zero heat bypasses these unused controls; a single
/// active film uses its existing single-film solver and work ceiling instead.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasLiquidChamberHeatSolve {
    pub maximum_iterations: usize,
}

/// Exact intermediate caloric state for one retained film, in stable ID order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasLiquidChamberSurfaceHeatStep {
    pub id: u64,
    /// Actual state passed to the phase solver, without a Q/m round trip.
    pub phase_initial_liquid: LiquidSurfaceState,
    pub heat_into_liquid_j: f64,
    pub liquid_enthalpy_change_j: f64,
    pub temperature_k: [f64; 2],
    pub temperature_equation_residual_k: f64,
    /// Total nested liquid-response iterations, including rejected outer trials.
    /// The singleton owner has no nested response solve and reports zero here.
    pub iterations: usize,
}

/// One shared heat-stage account, separate from all vapor-carried energy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GasLiquidChamberHeatStep {
    pub heat_into_liquid_j: f64,
    pub liquid_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub sum_absolute_energy_residual_j: f64,
    pub maximum_temperature_equation_residual_k: f64,
    pub gas_temperature_k: [f64; 2],
    pub iterations: usize,
    pub liquid_iterations: usize,
}

/// Original heat THEN phase endpoints and their exact intermediate liquid states.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidChamberHeatCandidate<'catalog> {
    phase: GasLiquidChamberCandidate<'catalog>,
    surfaces: Vec<GasLiquidChamberSurfaceStep>,
    step: GasLiquidChamberStep,
    heat: GasLiquidChamberHeatStep,
    heat_surfaces: Vec<GasLiquidChamberSurfaceHeatStep>,
}
impl<'catalog> GasLiquidChamberHeatCandidate<'catalog> {
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.phase.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.phase.catalog()
    }
    /// Phase-only receipt, starting from the exact heated states.
    pub const fn phase(&self) -> &GasLiquidChamberCandidate<'catalog> {
        &self.phase
    }
    /// Combined original-to-final film reports. Carried energy remains vapor only.
    pub fn surfaces(&self) -> &[GasLiquidChamberSurfaceStep] {
        &self.surfaces
    }
    /// Combined original-to-final chamber report. Sensible heat is counted once
    /// in storage changes; energy_into_gas_j remains the vapor stream only.
    pub const fn step(&self) -> &GasLiquidChamberStep {
        &self.step
    }
    pub const fn heat(&self) -> &GasLiquidChamberHeatStep {
        &self.heat
    }
    pub fn heat_surfaces(&self) -> &[GasLiquidChamberSurfaceHeatStep] {
        &self.heat_surfaces
    }
    pub fn into_parts(
        self,
    ) -> (
        Vec<GasLiquidChamberSurfaceStep>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasLiquidChamberStep,
        GasLiquidChamberHeatStep,
        Vec<GasLiquidChamberSurfaceHeatStep>,
    ) {
        let (_, gas, catalog, _) = self.phase.into_parts();
        (
            self.surfaces,
            gas,
            catalog,
            self.step,
            self.heat,
            self.heat_surfaces,
        )
    }
}

struct HeatEndpoint {
    gas: ClosedGasParcel,
    heat: GasLiquidChamberHeatStep,
    surfaces: Vec<GasLiquidChamberSurfaceHeatStep>,
}

/// Reciprocal heat on one shared gas capacity, then the existing joint phase
/// owner. Every original source/phase binding is admitted before heating.
pub fn advance_chamber_with_gas_and_heat<'catalog>(
    surfaces: &[GasLiquidChamberHeatSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    phase_controls: GasLiquidChamberSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasLiquidChamberHeatCandidate<'catalog>, GasEvaporationError> {
    advance_chamber_with_gas_and_heat_mass_selection(
        surfaces,
        gas,
        catalog,
        duration_s,
        phase_controls,
        heat_controls,
        |_, mass| Ok(mass),
        cancelled,
    )
}

/// Select each native representable phase mass exactly once AFTER joint heat,
/// before the sole joint phase caloric solve. Initial zero-duration admission
/// uses identity mass selection and does not consume the caller's selector.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_with_gas_and_heat_mass_selection<'catalog>(
    surfaces: &[GasLiquidChamberHeatSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    phase_controls: GasLiquidChamberSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    mut select_mass: impl FnMut(u64, f64) -> Result<f64, GasEvaporationError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasLiquidChamberHeatCandidate<'catalog>, GasEvaporationError> {
    if surfaces.len() > phase_controls.maximum_surfaces {
        return Err(invalid("joint heat surface count exceeds its work ceiling"));
    }
    if surfaces
        .iter()
        .any(|s| !nonnegative(s.heat_transfer_coefficient_w_m2_k))
    {
        return Err(invalid(
            "finite nonnegative joint sensible heat coefficients required",
        ));
    }
    let mut ordered: Vec<_> = surfaces.iter().collect();
    ordered.sort_unstable_by_key(|surface| surface.surface.id);
    let original_surfaces: Vec<_> = ordered.iter().map(|surface| surface.surface).collect();
    if ordered
        .iter()
        .all(|surface| surface.heat_transfer_coefficient_w_m2_k == 0.0)
    {
        let phase = advance_chamber_with_gas_mass_selection(
            &original_surfaces,
            gas,
            catalog,
            duration_s,
            phase_controls,
            select_mass,
            &mut cancelled,
        )?;
        let heat_surfaces = ordered
            .iter()
            .zip(phase.surfaces())
            .map(|(surface, step)| GasLiquidChamberSurfaceHeatStep {
                id: surface.surface.id,
                phase_initial_liquid: surface.surface.liquid,
                heat_into_liquid_j: 0.0,
                liquid_enthalpy_change_j: 0.0,
                temperature_k: [step.surface_temperature_k[0]; 2],
                temperature_equation_residual_k: 0.0,
                iterations: 0,
            })
            .collect();
        poll(&mut cancelled)?;
        return Ok(GasLiquidChamberHeatCandidate {
            surfaces: phase.surfaces().to_vec(),
            step: phase.step().clone(),
            phase,
            heat: GasLiquidChamberHeatStep {
                gas_temperature_k: [gas.state().state().temperature; 2],
                ..GasLiquidChamberHeatStep::default()
            },
            heat_surfaces,
        });
    }
    poll(&mut cancelled)?;
    // Original zero-duration phase admission validates every caloric datum and
    // full gas catalog before a thermal change could conceal an invalid input.
    let original = advance_chamber_with_gas(
        &original_surfaces,
        gas,
        catalog,
        0.0,
        phase_controls,
        &mut cancelled,
    )?;
    super::super::validate_controls(phase_controls.surface, duration_s, 0.0)?;
    let heated = solve::sensible_heat(
        &ordered,
        original.surfaces(),
        gas,
        duration_s,
        phase_controls,
        heat_controls,
        &mut cancelled,
    )?;
    let heated_surfaces: Vec<_> = ordered
        .iter()
        .zip(&heated.surfaces)
        .map(|(surface, heat)| GasLiquidChamberSurface {
            liquid: heat.phase_initial_liquid,
            ..surface.surface
        })
        .collect();
    let phase = advance_chamber_with_gas_mass_selection(
        &heated_surfaces,
        &heated.gas,
        catalog,
        duration_s,
        phase_controls,
        &mut select_mass,
        &mut cancelled,
    )?;
    let mut combined = phase.surfaces().to_vec();
    for ((surface, heat), initial) in combined
        .iter_mut()
        .zip(&heated.surfaces)
        .zip(original.surfaces())
    {
        surface.surface_temperature_k[0] = initial.surface_temperature_k[0];
        surface.vapor_partial_pressure_pa[0] = initial.vapor_partial_pressure_pa[0];
        surface.liquid_enthalpy_change_j =
            finite(surface.liquid_enthalpy_change_j + heat.liquid_enthalpy_change_j)?;
        if (surface.surface_temperature_k[1] - surface.surface_temperature_k[0]).abs()
            > phase_controls.surface.maximum_temperature_increment_k
        {
            return Err(smaller_step(
                duration_s,
                "combined joint heat and phase temperature departure",
            ));
        }
    }
    let mut step = phase.step().clone();
    step.gas_temperature_k[0] = gas.state().state().temperature;
    step.liquid_enthalpy_change_j =
        finite(step.liquid_enthalpy_change_j + heated.heat.liquid_enthalpy_change_j)?;
    step.gas_internal_energy_change_j =
        finite(step.gas_internal_energy_change_j + heated.heat.gas_internal_energy_change_j)?;
    step.energy_residual_j = finite(step.energy_residual_j + heated.heat.energy_residual_j)?;
    step.sum_absolute_energy_residual_j =
        finite(step.sum_absolute_energy_residual_j + heated.heat.sum_absolute_energy_residual_j)?;
    if step.sum_absolute_energy_residual_j > phase_controls.surface.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: step.sum_absolute_energy_residual_j,
            tolerance_j: phase_controls.surface.gas.energy_tolerance_j,
        });
    }
    if (step.gas_temperature_k[1] - step.gas_temperature_k[0]).abs()
        > phase_controls.surface.maximum_temperature_increment_k
    {
        return Err(smaller_step(
            duration_s,
            "combined joint heat and phase temperature departure",
        ));
    }
    poll(&mut cancelled)?;
    Ok(GasLiquidChamberHeatCandidate {
        phase,
        surfaces: combined,
        step,
        heat: heated.heat,
        heat_surfaces: heated.surfaces,
    })
}

#[cfg(test)]
mod tests;
