//! Reciprocal sensible gas-film heat on the original liquid and gas calorics.
//!
//! Backward Euler solves `Q = h*A*dt*(T_g1-T_l1)` with
//! `m_l*(h_l1-h_l0)=Q` and `m_g*Delta u_g=-Q`. Positive Q heats the liquid.
//! Original NASA-9 increments cancel formation references before evaluating
//! energy. Scalar trials use the original gas Delta-u and liquid h-to-T inverse;
//! no frozen heat capacity, second EOS or per-trial source certification is used.
//!
//! Heat acts first, then the existing signed phase-transfer step acts over the
//! same interval. This Lie split requires full-timestep refinement. The supplied
//! sensible coefficient is additional to vapor-carried enthalpy: it must NOT
//! include that advected sensible heat or latent heat a second time. Gas volume,
//! liquid reference geometry and exposed area retain their original owners.

use super::{
    GasEvaporationCandidate, GasEvaporationError, GasEvaporationSolve, GasEvaporationStep,
    GasLiquidBinding, advance_surface_with_gas, advance_surface_with_gas_mass_selection,
    nonnegative, poll, positive, smaller_step,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_qty::Temperature;
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, LiquidSurfaceState};

/// Explicit original liquid branch and bounded implicit heat work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasFilmHeatSolve {
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
    /// Absolute residual of `(T_g1-T_l1)-Q/(h*A*dt)` [K]. This is distinct
    /// from conservation: the combined heat/phase energy defect must fit the
    /// ONE absolute joule allowance in `GasEvaporationSolve::gas`.
    pub temperature_tolerance_k: f64,
    /// Maximum scalar bisection trials, in 1..=256; original caloric solves
    /// retain their own bounded work. A failed solve publishes neither state.
    pub maximum_iterations: usize,
}

/// Only the reciprocal heat substep; vapor-carried energy is reported separately.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GasFilmHeatStep {
    /// Represented heat delivered to the liquid and debited from the gas [J].
    pub heat_into_liquid_j: f64,
    pub liquid_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub temperature_equation_residual_k: f64,
    pub iterations: usize,
}

/// Both original endpoints after heat THEN phase, accepted or discarded together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEvaporationHeatCandidate<'catalog> {
    phase: GasEvaporationCandidate<'catalog>,
    phase_initial_liquid: LiquidSurfaceState,
    heat: GasFilmHeatStep,
    step: GasEvaporationStep,
}
impl<'catalog> GasEvaporationHeatCandidate<'catalog> {
    pub const fn liquid(&self) -> LiquidSurfaceState {
        self.phase.liquid()
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.phase.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.phase.catalog()
    }
    /// Original-to-final temperatures and TOTAL storage/residual changes.
    /// `energy_into_gas_j` remains the advected vapor stream ONLY; net energy
    /// into the gas is that stream minus `heat().heat_into_liquid_j`.
    pub const fn step(&self) -> &GasEvaporationStep {
        &self.step
    }
    pub const fn heat(&self) -> &GasFilmHeatStep {
        &self.heat
    }
    /// Phase-only receipt, whose starting states already include sensible heat.
    pub const fn phase(&self) -> &GasEvaporationCandidate<'catalog> {
        &self.phase
    }
    /// Exact liquid state passed from sensible heat into the phase owner.
    /// Its enthalpy is retained directly, without a heat/mass round trip.
    pub const fn phase_initial_liquid(&self) -> LiquidSurfaceState {
        self.phase_initial_liquid
    }
    pub fn into_parts(
        self,
    ) -> (
        LiquidSurfaceState,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasEvaporationStep,
        GasFilmHeatStep,
    ) {
        let (liquid, gas, catalog, _) = self.phase.into_parts();
        (liquid, gas, catalog, self.step, self.heat)
    }
}

pub(crate) fn validate_heat_controls(
    coefficient: f64,
    controls: GasFilmHeatSolve,
) -> Result<(), GasEvaporationError> {
    let [lo, hi] = controls.pure_liquid_enthalpy_range_j_kg;
    if !nonnegative(coefficient)
        || !lo.is_finite()
        || !hi.is_finite()
        || lo > hi
        || !positive(controls.temperature_tolerance_k)
        || !(1..=256).contains(&controls.maximum_iterations)
    {
        return Err(GasEvaporationError::Invalid(
            "finite nonnegative heat coefficient, liquid branch and positive bounded heat controls required",
        ));
    }
    Ok(())
}
fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GasEvaporationError::Invalid(
            "sensible film heat is unrepresentable",
        ))
    }
}

#[derive(Clone, Copy)]
struct Sample {
    gas_t: f64,
    liquid_t: f64,
    h: f64,
    q: f64,
    residual: f64,
}

#[allow(clippy::too_many_arguments)]
fn sensible_heat(
    liquid: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    liquid_t: f64,
    temperature: &impl Fn(f64) -> Result<f64, EvaporationError>,
    coefficient: f64,
    duration_s: f64,
    phase: GasEvaporationSolve,
    controls: GasFilmHeatSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(LiquidSurfaceState, ClosedGasParcel, GasFilmHeatStep), GasEvaporationError> {
    use GasEvaporationError as Error;
    validate_heat_controls(coefficient, controls)?;
    let [lo, hi] = controls.pure_liquid_enthalpy_range_j_kg;
    if !(lo..=hi).contains(&liquid.specific_enthalpy_j_kg) {
        return Err(Error::Surface(EvaporationError::Invalid(
            "initial liquid leaves heat branch",
        )));
    }
    let gas_t = gas.state().state().temperature;
    let difference = finite(gas_t - liquid_t)?;
    if difference == 0.0 || coefficient == 0.0 || liquid.area_m2 == 0.0 || duration_s == 0.0 {
        return Ok((liquid, gas.clone(), GasFilmHeatStep::default()));
    }
    let conductance_time = finite(coefficient * liquid.area_m2 * duration_s)?;
    if conductance_time == 0.0 {
        return Err(Error::Invalid("nonzero heat exposure underflows"));
    }
    let sign = difference.signum();
    let mut gas_controls = phase.gas.gas;
    gas_controls.energy_tolerance_j_kg = gas_controls
        .energy_tolerance_j_kg
        .min(0.25 * phase.gas.energy_tolerance_j / gas.mass_kg());
    let branch = gas
        .model()
        .prepare_caloric_branch(gas_controls, &mut *cancelled)?;
    if !(gas_controls.minimum_temperature_k..=gas_controls.maximum_temperature_k).contains(&gas_t) {
        return Err(Error::Invalid(
            "heat gas branch must contain the original temperature",
        ));
    }
    let energy = |t: f64| -> Result<f64, Error> {
        let du = branch
            .model()
            .thermodynamics()
            .specific_internal_energy_change(Temperature::new(gas_t), Temperature::new(t))
            .map_err(|error| {
                super::GasEnergyError::Material(crate::MaterialError::Parameters {
                    what: format!("original gas film heat increment: {error}"),
                })
            })?
            .value();
        finite(-gas.mass_kg() * du)
    };
    // Restrict the search to BOTH original source domains. The gas may approach
    // the initial liquid temperature but cannot pass it during a heat-only step.
    let gas_limit = liquid_t
        .max(gas_controls.minimum_temperature_k)
        .min(gas_controls.maximum_temperature_k);
    let gas_q = energy(gas_limit)?;
    let liquid_limit = if sign > 0.0 { hi } else { lo };
    let liquid_q = finite(liquid.mass_kg * (liquid_limit - liquid.specific_enthalpy_j_kg))?;
    let (far_t, far_q, far_h) = if sign * gas_q <= sign * liquid_q {
        (gas_limit, gas_q, None)
    } else {
        let inverse = branch.reconstruct_density_internal_energy_increment(
            gas.point(),
            gas.mass_kg() / gas.volume_m3(),
            -liquid_q / gas.mass_kg(),
            &mut *cancelled,
        )?;
        (
            inverse.state().state().temperature,
            liquid_q,
            Some(liquid_limit),
        )
    };
    let evaluate = |t: f64, q: f64, boundary_h: Option<f64>| -> Result<Sample, Error> {
        let dh = finite(q / liquid.mass_kg)?;
        // Preserve an admitted source endpoint exactly: multiplying its h gap
        // by mass and dividing back can round one ulp outside that branch.
        let h = finite(boundary_h.unwrap_or(liquid.specific_enthalpy_j_kg + dh))?;
        if !(lo..=hi).contains(&h) {
            return Err(Error::Surface(EvaporationError::Invalid(
                "sensible heat leaves liquid branch",
            )));
        }
        if q != 0.0 && (dh == 0.0 || h == liquid.specific_enthalpy_j_kg) {
            return Err(Error::Invalid(
                "sensible liquid caloric increment rounds away",
            ));
        }
        let liquid_t = temperature(h)?;
        if !positive(liquid_t) {
            return Err(Error::Invalid("positive liquid heat temperature required"));
        }
        Ok(Sample {
            gas_t: t,
            liquid_t,
            h,
            q,
            residual: finite((t - liquid_t) - q / conductance_time)?,
        })
    };
    let mut near = Sample {
        gas_t,
        liquid_t,
        h: liquid.specific_enthalpy_j_kg,
        q: 0.0,
        residual: difference,
    };
    let mut far = evaluate(far_t, far_q, far_h)?;
    // Leave half the equation allowance for publishing the gas endpoint through
    // its original caloric inverse and for represented liquid enthalpy rounding.
    let equation_tolerance = 0.5 * controls.temperature_tolerance_k;
    if sign * far.residual > equation_tolerance || sign * far_q <= 0.0 {
        return Err(Error::HeatOutsideBracket {
            residual_temperature_k: [difference, far.residual],
        });
    }
    // Accept from the initial ordering side: an implicit heat step must not
    // reverse the gas/liquid temperature ordering even with a long duration.
    let converged = |sample: Sample| {
        sample.q != 0.0
            && sign * (sample.gas_t - sample.liquid_t) >= 0.0
            && sample.residual.abs() <= equation_tolerance
    };
    let mut accepted = if converged(far) { Some(far) } else { None };
    let mut iterations = 0;
    while accepted.is_none() && iterations < controls.maximum_iterations {
        poll(cancelled)?;
        let t = near.gas_t + 0.5 * (far.gas_t - near.gas_t);
        if t == near.gas_t || t == far.gas_t {
            break;
        }
        let sample = evaluate(t, energy(t)?, None)?;
        iterations += 1;
        if converged(sample) {
            accepted = Some(sample);
        } else if sign * sample.residual > 0.0 {
            near = sample;
        } else {
            far = sample;
        }
    }
    let endpoint = accepted.ok_or(Error::HeatNotConverged {
        iterations,
        residual_temperature_k: near.residual,
    })?;
    // Debit the SAME represented liquid heat, then let the original gas owner
    // reconstruct every gas property once. No model/card cloning occurs in trials.
    let q = finite(liquid.mass_kg * (endpoint.h - liquid.specific_enthalpy_j_kg))?;
    if sign * q <= 0.0 {
        return Err(Error::Invalid("nonzero reciprocal heat rounds away"));
    }
    gas_controls.energy_tolerance_j_kg = gas_controls
        .energy_tolerance_j_kg
        .min(0.25 * controls.temperature_tolerance_k * branch.minimum_cv_j_kg_k());
    let (heated_gas, stored) = gas.at_volume_internal_energy_increment(
        gas.volume_m3(),
        -q,
        gas_controls,
        &mut *cancelled,
    )?;
    let final_gas_t = heated_gas.state().state().temperature;
    let residual_t = finite((final_gas_t - endpoint.liquid_t) - q / conductance_time)?;
    if residual_t.abs() > controls.temperature_tolerance_k
        || sign * (final_gas_t - endpoint.liquid_t) < 0.0
    {
        return Err(Error::HeatNotConverged {
            iterations,
            residual_temperature_k: residual_t,
        });
    }
    let residual_j = finite(q + stored)?;
    if residual_j.abs() > phase.gas.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j,
            tolerance_j: phase.gas.energy_tolerance_j,
        });
    }
    if (final_gas_t - gas_t).abs() > phase.maximum_temperature_increment_k
        || (endpoint.liquid_t - liquid_t).abs() > phase.maximum_temperature_increment_k
    {
        return Err(smaller_step(
            duration_s,
            "sensible heat temperature departure",
        ));
    }
    poll(cancelled)?;
    Ok((
        LiquidSurfaceState {
            specific_enthalpy_j_kg: endpoint.h,
            ..liquid
        },
        heated_gas,
        GasFilmHeatStep {
            heat_into_liquid_j: q,
            liquid_enthalpy_change_j: q,
            gas_internal_energy_change_j: stored,
            energy_residual_j: residual_j,
            temperature_equation_residual_k: residual_t,
            iterations,
        },
    ))
}

/// Implicit reciprocal heat, then the original signed phase step over the SAME
/// duration. Refine the complete physical step to assess splitting error.
/// The caller admits one continuous increasing original liquid `T(h)` branch.
/// A zero heat coefficient delegates bit-for-bit to `advance_surface_with_gas`;
/// no heat-source branch or numerical control is inferred in that path.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface_with_gas_and_heat<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    liquid: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, EvaporationError>,
    mass_transfer_coefficient_m_s: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    duration_s: f64,
    phase_controls: GasEvaporationSolve,
    heat_controls: GasFilmHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasEvaporationHeatCandidate<'catalog>, GasEvaporationError> {
    advance_surface_with_gas_and_heat_mass_selection(
        law,
        binding,
        liquid,
        gas,
        catalog,
        temperature_from_enthalpy,
        mass_transfer_coefficient_m_s,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        phase_controls,
        heat_controls,
        Ok,
        cancelled,
    )
}

/// The same reciprocal heat and phase step with a native mass representation.
/// The selector receives the phase transfer requested AFTER sensible heating,
/// exactly once on a successful call, before either phase caloric update.
/// Initial zero-duration phase admission keeps its identity transfer and never
/// consumes the selector. All original selected-mass limits remain active.
/// A zero heat coefficient delegates directly to the selected phase owner.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface_with_gas_and_heat_mass_selection<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    liquid: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, EvaporationError>,
    mass_transfer_coefficient_m_s: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    duration_s: f64,
    phase_controls: GasEvaporationSolve,
    heat_controls: GasFilmHeatSolve,
    select_mass: impl FnOnce(f64) -> Result<f64, GasEvaporationError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasEvaporationHeatCandidate<'catalog>, GasEvaporationError> {
    if !nonnegative(heat_transfer_coefficient_w_m2_k) {
        return Err(GasEvaporationError::Invalid(
            "finite nonnegative sensible heat coefficient required",
        ));
    }
    if heat_transfer_coefficient_w_m2_k == 0.0 {
        let phase = advance_surface_with_gas_mass_selection(
            law,
            binding,
            liquid,
            gas,
            catalog,
            temperature_from_enthalpy,
            mass_transfer_coefficient_m_s,
            duration_s,
            phase_controls,
            select_mass,
            cancelled,
        )?;
        return Ok(GasEvaporationHeatCandidate {
            step: phase.step().clone(),
            phase,
            phase_initial_liquid: liquid,
            heat: GasFilmHeatStep::default(),
        });
    }
    poll(&mut cancelled)?;
    // Validate ORIGINAL phase/caloric compatibility before heat can shift it.
    super::validate_controls(phase_controls, duration_s, mass_transfer_coefficient_m_s)?;
    let original = advance_surface_with_gas(
        law,
        binding,
        liquid,
        gas,
        catalog,
        &temperature_from_enthalpy,
        mass_transfer_coefficient_m_s,
        0.0,
        phase_controls,
        &mut cancelled,
    )?;
    let initial_t = original.step().surface_temperature_k[0];
    let (heated_liquid, heated_gas, heat) = sensible_heat(
        liquid,
        gas,
        initial_t,
        &temperature_from_enthalpy,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        phase_controls,
        heat_controls,
        &mut cancelled,
    )?;
    let phase = advance_surface_with_gas_mass_selection(
        law,
        binding,
        heated_liquid,
        &heated_gas,
        catalog,
        temperature_from_enthalpy,
        mass_transfer_coefficient_m_s,
        duration_s,
        phase_controls,
        select_mass,
        &mut cancelled,
    )?;
    let mut step = phase.step().clone();
    step.surface_temperature_k[0] = initial_t;
    step.gas_temperature_k[0] = gas.state().state().temperature;
    step.vapor_partial_pressure_pa[0] = original.step().vapor_partial_pressure_pa[0];
    step.liquid_enthalpy_change_j =
        finite(step.liquid_enthalpy_change_j + heat.liquid_enthalpy_change_j)?;
    step.gas_internal_energy_change_j =
        finite(step.gas_internal_energy_change_j + heat.gas_internal_energy_change_j)?;
    step.energy_residual_j = finite(step.energy_residual_j + heat.energy_residual_j)?;
    if step.energy_residual_j.abs() > phase_controls.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: step.energy_residual_j,
            tolerance_j: phase_controls.gas.energy_tolerance_j,
        });
    }
    if (step.surface_temperature_k[1] - step.surface_temperature_k[0]).abs()
        > phase_controls.maximum_temperature_increment_k
        || (step.gas_temperature_k[1] - step.gas_temperature_k[0]).abs()
            > phase_controls.maximum_temperature_increment_k
    {
        return Err(smaller_step(
            duration_s,
            "combined heat and phase temperature departure",
        ));
    }
    poll(&mut cancelled)?;
    Ok(GasEvaporationHeatCandidate {
        phase,
        phase_initial_liquid: heated_liquid,
        heat,
        step,
    })
}

#[cfg(test)]
mod tests;
