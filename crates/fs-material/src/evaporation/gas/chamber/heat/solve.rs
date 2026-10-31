//! Monotone scalar reduction of the original nonlinear caloric star network.
use super::super::Total;
use super::{
    ClosedGasParcel, GasEvaporationError, GasLiquidChamberHeatSolve, GasLiquidChamberHeatStep,
    GasLiquidChamberHeatSurface, GasLiquidChamberSolve, GasLiquidChamberSurfaceHeatStep,
    GasLiquidChamberSurfaceStep, HeatEndpoint, finite, invalid, poll, smaller_step,
};
use crate::evaporation::gas::heat::{sensible_heat as single_heat, validate_heat_controls};
use crate::evaporation::gas::positive;
use fs_qty::Temperature;
use fs_thermal::evaporation::EvaporationError;

#[derive(Clone, Copy)]
struct LiquidSample {
    h: f64,
    t: f64,
    q: f64,
}
struct Response {
    index: usize,
    h0: f64,
    t0: f64,
    mass: f64,
    k: f64,
    lo: LiquidSample,
    hi: LiquidSample,
    temperature_tolerance: f64,
    maximum_iterations: usize,
}
impl Response {
    fn residual(&self, sample: LiquidSample, gas_t: f64) -> Result<f64, GasEvaporationError> {
        finite((sample.t - gas_t) + sample.q / self.k)
    }
    fn sample(
        &self,
        h: f64,
        temperature: &dyn Fn(f64) -> Result<f64, EvaporationError>,
    ) -> Result<LiquidSample, GasEvaporationError> {
        let q = finite(self.mass * (h - self.h0))?;
        if h != self.h0 && q == 0.0 {
            return Err(invalid("nonzero joint liquid heat rounds away"));
        }
        let t = temperature(h)?;
        if !positive(t) {
            return Err(invalid("positive original liquid temperature required"));
        }
        Ok(LiquidSample { h, t, q })
    }
}

/// Interior secant proposals are safeguarded by a fixed bracket fraction.
/// This uses no inferred heat capacity or numerical derivative.
fn interior(lo: f64, hi: f64, rlo: f64, rhi: f64) -> f64 {
    let fraction = -rlo / (rhi - rlo);
    let fraction = if fraction.is_finite() && (0.1..=0.9).contains(&fraction) {
        fraction
    } else {
        0.5
    };
    lo + fraction * (hi - lo)
}

#[allow(clippy::too_many_arguments)]
fn liquid_response(
    response: &Response,
    gas_t: f64,
    temperature: &dyn Fn(f64) -> Result<f64, EvaporationError>,
    energy_allowance: f64,
    warm: Option<(LiquidSample, LiquidSample)>,
    work: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(LiquidSample, f64), GasEvaporationError> {
    poll(cancelled)?;
    if gas_t == response.t0 {
        return Ok((
            LiquidSample {
                h: response.h0,
                t: response.t0,
                q: 0.0,
            },
            0.0,
        ));
    }
    let (mut lo, mut hi) = (response.lo, response.hi);
    if let Some((a, b)) = warm {
        if a.h <= b.h && response.residual(a, gas_t)? <= 0.0 && response.residual(b, gas_t)? >= 0.0
        {
            (lo, hi) = (a, b);
        }
    }
    let mut rlo = response.residual(lo, gas_t)?;
    let mut rhi = response.residual(hi, gas_t)?;
    if rlo > 0.0 || rhi < 0.0 {
        return Err(GasEvaporationError::HeatOutsideBracket {
            residual_temperature_k: [rlo, rhi],
        });
    }
    let converged = |sample: LiquidSample, residual: f64, width: f64| {
        // Since T(h) increases, dF/dQ >= 1/K. Either the retained Q bracket
        // or K*|F| bounds the local heat error without guessing a derivative.
        let uncertainty = width.abs().min(response.k * residual.abs());
        let direction = sample.q == 0.0 || sample.q.signum() * (gas_t - sample.t) >= 0.0;
        if sample.q != 0.0
            && direction
            && residual.abs() <= 0.25 * response.temperature_tolerance
            && uncertainty <= energy_allowance
        {
            Some(uncertainty)
        } else {
            None
        }
    };
    if let Some(error) = converged(lo, rlo, hi.q - lo.q) {
        return Ok((lo, error));
    }
    if let Some(error) = converged(hi, rhi, hi.q - lo.q) {
        return Ok((hi, error));
    }
    let mut count = 0usize;
    while count < response.maximum_iterations {
        poll(cancelled)?;
        let h = interior(lo.h, hi.h, rlo, rhi);
        if h == lo.h || h == hi.h {
            break;
        }
        let sample = response.sample(h, temperature)?;
        count += 1;
        *work = work
            .checked_add(1)
            .ok_or(invalid("joint heat work overflow"))?;
        let residual = response.residual(sample, gas_t)?;
        if let Some(error) = converged(sample, residual, hi.q - lo.q) {
            return Ok((sample, error));
        }
        if residual < 0.0 {
            lo = sample;
            rlo = residual;
        } else {
            hi = sample;
            rhi = residual;
        }
    }
    Err(GasEvaporationError::HeatNotConverged {
        iterations: count,
        residual_temperature_k: if rlo.abs() <= rhi.abs() { rlo } else { rhi },
    })
}

#[derive(Clone)]
struct Trial {
    gas_t: f64,
    liquids: Vec<LiquidSample>,
    q: f64,
    residual_j: f64,
    uncertainty_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn trial(
    gas_t: f64,
    responses: &[Response],
    surfaces: &[&GasLiquidChamberHeatSurface<'_>],
    energy: &impl Fn(f64) -> Result<f64, GasEvaporationError>,
    local_energy_allowance: f64,
    warm: Option<(&Trial, &Trial)>,
    work: &mut [usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial, GasEvaporationError> {
    poll(cancelled)?;
    let mut liquids = Vec::with_capacity(responses.len());
    let mut total = Total::default();
    let mut uncertainty = Total::default();
    for (index, response) in responses.iter().enumerate() {
        let hint = warm.map(|(lo, hi)| (lo.liquids[index], hi.liquids[index]));
        let (liquid, error) = liquid_response(
            response,
            gas_t,
            surfaces[response.index].surface.temperature_from_enthalpy,
            local_energy_allowance,
            hint,
            &mut work[index],
            cancelled,
        )?;
        total.add(liquid.q)?;
        uncertainty.add(error)?;
        liquids.push(liquid);
    }
    let q = total.value()?;
    Ok(Trial {
        gas_t,
        liquids,
        q,
        residual_j: finite(energy(gas_t)? + q)?,
        uncertainty_j: uncertainty.value()?,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn sensible_heat(
    surfaces: &[&GasLiquidChamberHeatSurface<'_>],
    original: &[GasLiquidChamberSurfaceStep],
    gas: &ClosedGasParcel,
    duration_s: f64,
    phase: GasLiquidChamberSolve,
    controls: GasLiquidChamberHeatSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HeatEndpoint, GasEvaporationError> {
    poll(cancelled)?;
    let gas_t = gas.state().state().temperature;
    let mut reports: Vec<_> = surfaces
        .iter()
        .zip(original)
        .map(|(surface, initial)| GasLiquidChamberSurfaceHeatStep {
            id: surface.surface.id,
            phase_initial_liquid: surface.surface.liquid,
            heat_into_liquid_j: 0.0,
            liquid_enthalpy_change_j: 0.0,
            temperature_k: [initial.surface_temperature_k[0]; 2],
            temperature_equation_residual_k: 0.0,
            iterations: 0,
        })
        .collect();
    let mut active = Vec::new();
    for (index, surface) in surfaces.iter().enumerate() {
        poll(cancelled)?;
        let coefficient = surface.heat_transfer_coefficient_w_m2_k;
        if coefficient == 0.0 {
            continue;
        }
        validate_heat_controls(coefficient, surface.heat_controls)?;
        let liquid = surface.surface.liquid;
        let [lo, hi] = surface.heat_controls.pure_liquid_enthalpy_range_j_kg;
        if !(lo..=hi).contains(&liquid.specific_enthalpy_j_kg) {
            return Err(invalid("initial liquid leaves joint heat branch"));
        }
        if liquid.area_m2 == 0.0 || duration_s == 0.0 {
            continue;
        }
        let k = finite(coefficient * liquid.area_m2 * duration_s)?;
        if k == 0.0 {
            return Err(invalid("nonzero joint heat exposure underflows"));
        }
        active.push((index, k));
    }
    if active.is_empty() {
        return Ok(HeatEndpoint {
            gas: gas.clone(),
            heat: GasLiquidChamberHeatStep {
                gas_temperature_k: [gas_t; 2],
                ..GasLiquidChamberHeatStep::default()
            },
            surfaces: reports,
        });
    }
    if active.len() == 1 {
        // Reuse the established singleton arithmetic, domain and iteration
        // semantics exactly; disconnected films consume no gas heat capacity.
        let index = active[0].0;
        let surface = surfaces[index];
        let (liquid, heated_gas, heat) = single_heat(
            surface.surface.liquid,
            gas,
            reports[index].temperature_k[0],
            &surface.surface.temperature_from_enthalpy,
            surface.heat_transfer_coefficient_w_m2_k,
            duration_s,
            phase.surface,
            surface.heat_controls,
            cancelled,
        )?;
        reports[index].phase_initial_liquid = liquid;
        reports[index].heat_into_liquid_j = heat.heat_into_liquid_j;
        reports[index].liquid_enthalpy_change_j = heat.liquid_enthalpy_change_j;
        reports[index].temperature_k[1] =
            (surface.surface.temperature_from_enthalpy)(liquid.specific_enthalpy_j_kg)?;
        reports[index].temperature_equation_residual_k = heat.temperature_equation_residual_k;
        return Ok(HeatEndpoint {
            heat: GasLiquidChamberHeatStep {
                heat_into_liquid_j: heat.heat_into_liquid_j,
                liquid_enthalpy_change_j: heat.liquid_enthalpy_change_j,
                gas_internal_energy_change_j: heat.gas_internal_energy_change_j,
                energy_residual_j: heat.energy_residual_j,
                sum_absolute_energy_residual_j: heat.energy_residual_j.abs(),
                maximum_temperature_equation_residual_k: heat.temperature_equation_residual_k.abs(),
                gas_temperature_k: [gas_t, heated_gas.state().state().temperature],
                iterations: heat.iterations,
                liquid_iterations: 0,
            },
            gas: heated_gas,
            surfaces: reports,
        });
    }
    if !(1..=256).contains(&controls.maximum_iterations) {
        return Err(invalid("joint heat needs bounded positive outer work"));
    }
    let mut gas_controls = phase.surface.gas.gas;
    gas_controls.energy_tolerance_j_kg = gas_controls
        .energy_tolerance_j_kg
        .min(0.125 * phase.surface.gas.energy_tolerance_j / gas.mass_kg());
    let branch = gas
        .model()
        .prepare_caloric_branch(gas_controls, &mut *cancelled)?;
    if !(gas_controls.minimum_temperature_k..=gas_controls.maximum_temperature_k).contains(&gas_t) {
        return Err(invalid(
            "joint heat gas branch must contain the initial temperature",
        ));
    }
    let energy = |temperature| -> Result<f64, GasEvaporationError> {
        let change = branch
            .model()
            .thermodynamics()
            .specific_internal_energy_change(Temperature::new(gas_t), Temperature::new(temperature))
            .map_err(|error| {
                crate::thermochemical_gas::energy::GasEnergyError::Material(
                    crate::MaterialError::Parameters {
                        what: format!("original joint gas heat increment: {error}"),
                    },
                )
            })?
            .value();
        finite(gas.mass_kg() * change)
    };
    let mut hull_lo = gas_t;
    let mut hull_hi = gas_t;
    let mut responses = Vec::with_capacity(active.len());
    let mut minimum_equation_tolerance = f64::INFINITY;
    for (index, k) in active {
        poll(cancelled)?;
        let surface = surfaces[index];
        let liquid = surface.surface.liquid;
        let t0 = reports[index].temperature_k[0];
        hull_lo = hull_lo.min(t0);
        hull_hi = hull_hi.max(t0);
        minimum_equation_tolerance =
            minimum_equation_tolerance.min(surface.heat_controls.temperature_tolerance_k);
        let [lo, hi] = surface.heat_controls.pure_liquid_enthalpy_range_j_kg;
        let temperature = surface.surface.temperature_from_enthalpy;
        let make = |h| -> Result<LiquidSample, GasEvaporationError> {
            let t = temperature(h)?;
            if !positive(t) {
                return Err(invalid("positive liquid branch temperature required"));
            }
            Ok(LiquidSample {
                h,
                t,
                q: finite(liquid.mass_kg * (h - liquid.specific_enthalpy_j_kg))?,
            })
        };
        let lo = make(lo)?;
        let hi = make(hi)?;
        if lo.t > t0 || hi.t < t0 || lo.t > hi.t {
            return Err(invalid(
                "joint heat needs an admitted increasing original liquid branch",
            ));
        }
        responses.push(Response {
            index,
            h0: liquid.specific_enthalpy_j_kg,
            t0,
            mass: liquid.mass_kg,
            k,
            lo,
            hi,
            temperature_tolerance: surface.heat_controls.temperature_tolerance_k,
            maximum_iterations: surface.heat_controls.maximum_iterations,
        });
    }
    let mut lo = hull_lo.max(gas_controls.minimum_temperature_k);
    let mut hi = hull_hi.min(gas_controls.maximum_temperature_k);
    for response in &responses {
        lo = lo.max(finite(response.lo.t + response.lo.q / response.k)?);
        hi = hi.min(finite(response.hi.t + response.hi.q / response.k)?);
    }
    if lo > hi {
        return Err(invalid(
            "original heat domains have no shared gas-temperature interval",
        ));
    }
    let energy_allowance = (0.125 * phase.surface.gas.energy_tolerance_j)
        .min(0.125 * minimum_equation_tolerance * gas.mass_kg() * branch.minimum_cv_j_kg_k());
    let local_energy_allowance = energy_allowance / (4.0 * responses.len() as f64);
    if !positive(energy_allowance) || !positive(local_energy_allowance) {
        return Err(invalid("joint heat allowances are unrepresentable"));
    }
    let mut work = vec![0usize; responses.len()];
    let mut left = trial(
        lo,
        &responses,
        surfaces,
        &energy,
        local_energy_allowance,
        None,
        &mut work,
        cancelled,
    )?;
    let mut right = if lo == hi {
        left.clone()
    } else {
        trial(
            hi,
            &responses,
            surfaces,
            &energy,
            local_energy_allowance,
            None,
            &mut work,
            cancelled,
        )?
    };
    let mut accepted = if left.residual_j.abs() <= energy_allowance {
        Some(left.clone())
    } else if right.residual_j.abs() <= energy_allowance {
        Some(right.clone())
    } else {
        None
    };
    if accepted.is_none()
        && (left.residual_j - left.uncertainty_j > 0.0
            || right.residual_j + right.uncertainty_j < 0.0
            || lo == hi)
    {
        return Err(GasEvaporationError::HeatEnergyOutsideBracket {
            residual_energy_j: [left.residual_j, right.residual_j],
        });
    }
    let mut iterations = 0usize;
    while accepted.is_none() && iterations < controls.maximum_iterations {
        poll(cancelled)?;
        let t = interior(left.gas_t, right.gas_t, left.residual_j, right.residual_j);
        if t == left.gas_t || t == right.gas_t {
            break;
        }
        let sample = trial(
            t,
            &responses,
            surfaces,
            &energy,
            local_energy_allowance,
            Some((&left, &right)),
            &mut work,
            cancelled,
        )?;
        iterations += 1;
        if sample.residual_j.abs() <= energy_allowance {
            accepted = Some(sample);
        } else if sample.residual_j + sample.uncertainty_j < 0.0 {
            left = sample;
        } else if sample.residual_j - sample.uncertainty_j > 0.0 {
            right = sample;
        } else {
            return Err(invalid(
                "joint heat outer sign is unresolved by local caloric bounds",
            ));
        }
    }
    let accepted = accepted.ok_or(GasEvaporationError::HeatEnergyNotConverged {
        iterations,
        residual_energy_j: if left.residual_j.abs() <= right.residual_j.abs() {
            left.residual_j
        } else {
            right.residual_j
        },
    })?;
    gas_controls.energy_tolerance_j_kg = gas_controls
        .energy_tolerance_j_kg
        .min(0.125 * minimum_equation_tolerance * branch.minimum_cv_j_kg_k());
    // Actual represented liquid increments define the debit, not a rounded
    // absolute enthalpy or an approximate outer target.
    let (heated_gas, stored) = if accepted.q == 0.0 {
        (gas.clone(), 0.0)
    } else {
        gas.at_volume_internal_energy_increment(
            gas.volume_m3(),
            -accepted.q,
            gas_controls,
            &mut *cancelled,
        )?
    };
    let final_gas_t = heated_gas.state().state().temperature;
    if !(hull_lo..=hull_hi).contains(&final_gas_t) {
        return Err(invalid(
            "joint heat gas endpoint violates original temperature extrema",
        ));
    }
    let mut maximum_residual = 0.0f64;
    let mut liquid_iterations = 0usize;
    for ((response, sample), work) in responses.iter().zip(accepted.liquids).zip(work) {
        poll(cancelled)?;
        let residual = response.residual(sample, final_gas_t)?;
        if residual.abs() > response.temperature_tolerance
            || (sample.q != 0.0 && sample.q.signum() * (final_gas_t - sample.t) < 0.0)
        {
            return Err(GasEvaporationError::HeatNotConverged {
                iterations,
                residual_temperature_k: residual,
            });
        }
        if !(hull_lo..=hull_hi).contains(&sample.t) {
            return Err(invalid(
                "joint liquid heat endpoint violates original temperature extrema",
            ));
        }
        if (sample.t - response.t0).abs() > phase.surface.maximum_temperature_increment_k {
            return Err(smaller_step(
                duration_s,
                "joint sensible heat temperature departure",
            ));
        }
        maximum_residual = maximum_residual.max(residual.abs());
        liquid_iterations = liquid_iterations
            .checked_add(work)
            .ok_or(invalid("joint liquid heat work overflow"))?;
        let report = &mut reports[response.index];
        report.phase_initial_liquid.specific_enthalpy_j_kg = sample.h;
        report.heat_into_liquid_j = sample.q;
        report.liquid_enthalpy_change_j = sample.q;
        report.temperature_k[1] = sample.t;
        report.temperature_equation_residual_k = -residual;
        report.iterations = work;
    }
    if (final_gas_t - gas_t).abs() > phase.surface.maximum_temperature_increment_k {
        return Err(smaller_step(
            duration_s,
            "joint sensible heat temperature departure",
        ));
    }
    let residual = finite(accepted.q + stored)?;
    if residual.abs() > phase.surface.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: residual,
            tolerance_j: phase.surface.gas.energy_tolerance_j,
        });
    }
    poll(cancelled)?;
    Ok(HeatEndpoint {
        gas: heated_gas,
        heat: GasLiquidChamberHeatStep {
            heat_into_liquid_j: accepted.q,
            liquid_enthalpy_change_j: accepted.q,
            gas_internal_energy_change_j: stored,
            energy_residual_j: residual,
            sum_absolute_energy_residual_j: residual.abs(),
            maximum_temperature_equation_residual_k: maximum_residual,
            gas_temperature_k: [gas_t, final_gas_t],
            iterations,
            liquid_iterations,
        },
        surfaces: reports,
    })
}
