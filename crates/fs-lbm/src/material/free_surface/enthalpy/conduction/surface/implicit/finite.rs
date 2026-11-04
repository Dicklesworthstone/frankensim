//! One original isochoric gas equation beside original liquid and solid BE rows.
//!
//! Coordinates are liquid/solid specific-enthalpy increments followed by one gas
//! specific-internal-energy increment. Each private gas candidate is recovered
//! by its original caloric/EOS owner. Surface rates use THAT gas and the actual
//! liquid/solid endpoints together; no accepted heat/history changes during a trial.

use super::*;
use fs_conduction::enthalpy::backward_euler::BackwardEulerTrial;

struct CoupledTrial<'t, 'a, 'm> {
    thermal: BackwardEulerTrial<'t, Conduction<'a, 'm>>,
    gas_residual_j: f64,
    gas_stored_j: f64,
    independent_heat_j: f64,
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn solve_interval<'t, 'a, 'm>(
    thermal: &DepositionSimulation<'t>,
    observer: &Conduction<'a, 'm>,
    target_s: f64,
    thermal_scales_kg: &[f64],
    solve: FreeSurfaceImplicitSolve,
    total_tolerance_j: f64,
    evaluations: &Cell<usize>,
) -> Result<
    (
        DepositionSimulation<'t>,
        Conduction<'a, 'm>,
        EnergyAudit,
        NewtonReport,
        f64,
        f64,
    ),
    Error,
> {
    let check = || checkpoint(observer.gate);
    check()?;
    let gas = observer.finite_atmosphere.as_ref().ok_or(Error::Invalid(
        "coupled implicit surface requires its original finite gas",
    ))?;
    let step = BackwardEulerStep::new(thermal, target_s, check)?;
    let count = step.cell_ids().len();
    let dimension = count
        .checked_add(1)
        .ok_or(Error::Invalid("implicit gas dimension overflows"))?;
    validate(solve, dimension)?;
    if thermal_scales_kg.len() != count
        || thermal_scales_kg.iter().any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(Error::Invalid(
            "implicit liquid/solid row scales must be complete and positive",
        ));
    }
    let dt = environment::finite(target_s - thermal.time_s())?;
    let masses: Vec<_> = step
        .cell_ids()
        .iter()
        .map(|id| thermal.cells()[id].mass_kg)
        .collect();
    let mut gas_scale_kg = gas.parcel.mass_kg();
    for surface in &observer.surface_last {
        check()?;
        gas_scale_kg = environment::finite(
            gas_scale_kg + dt * surface.gas_decay_conductance_w_k * gas.temperature_slope,
        )?;
    }
    for surface in &observer.body_atmosphere.last {
        check()?;
        gas_scale_kg = environment::finite(
            gas_scale_kg + dt * surface.conductance_w_k * gas.temperature_slope,
        )?;
    }
    if gas_scale_kg <= 0.0 {
        return Err(Error::Invalid(
            "implicit gas conditioning requires positive actual capacity",
        ));
    }
    // Independent authored/molecular heat. Radiation is appended from each
    // trial's original liquid/solid samples below. Gas heat is internal to the union
    // and is deliberately absent from this independent external-energy account.
    let mut authored_power = Sum::default();
    for &power in observer.powers.values() {
        check()?;
        authored_power.add(power)?;
    }
    let authored_power = authored_power.finish()?;
    let evaluate = |increments: &[f64]| -> Result<CoupledTrial<'t, 'a, 'm>, Error> {
        check()?;
        if increments.len() != dimension {
            return Err(Error::Invalid(
                "implicit liquid/solid/gas coordinates are incomplete",
            ));
        }
        if evaluations.get() >= solve.maximum_residual_evaluations {
            return Err(Error::ImplicitSurface(
                FreeSurfaceImplicitFailure::EvaluationBudget,
            ));
        }
        evaluations.set(evaluations.get() + 1);
        let energy = increments[..count]
            .iter()
            .zip(&masses)
            .map(|(&h, &mass)| {
                let joules = environment::finite(h * mass)?;
                if h != 0.0 && joules == 0.0 {
                    return Err(Error::Invalid(
                        "implicit specific-energy increment underflows actual cell energy",
                    ));
                }
                Ok(joules)
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let mut payload = (*observer).clone();
        let gas_stored_j = payload
            .finite_atmosphere
            .as_mut()
            .expect("original finite atmosphere")
            .trial_energy_increment(increments[count], &mut || check())?;
        let thermal = step.evaluate(
            &energy,
            |endpoint| {
                let sources_w = payload.sources(endpoint, &mut || check())?;
                let conductances_w_k = payload.conductances(endpoint, &mut || check())?;
                Ok::<_, Error>(StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                })
            },
            check,
        )?;
        // Match the original observer's surface order and dt multiplication.
        // Its accepted accounting consumes these exact same endpoint samples.
        let mut gas_heat_j = 0.0;
        let mut external = Sum::default();
        external.add(environment::finite(dt * authored_power)?)?;
        for surface in &thermal.payload().surface_last {
            check()?;
            gas_heat_j = environment::finite(
                gas_heat_j + environment::finite(dt * surface.heat_from_gas_w)?,
            )?;
            external.add(environment::finite(dt * surface.heat_from_radiation_w)?)?;
        }
        // The same atmosphere supplies EVERY original solid gas patch. Match
        // the accepted observer's liquid-surface then body-surface ordering;
        // omitting even one patch would invent an external gas energy supply.
        for surface in &thermal.payload().body_atmosphere.last {
            check()?;
            gas_heat_j = environment::finite(
                gas_heat_j + environment::finite(dt * surface.heat_into_solid_w)?,
            )?;
        }
        for surface in &thermal.payload().body_radiation.last {
            check()?;
            external.add(environment::finite(dt * surface.heat_into_solid_w)?)?;
        }
        Ok(CoupledTrial {
            thermal,
            gas_residual_j: environment::finite(gas_stored_j + gas_heat_j)?,
            gas_stored_j,
            independent_heat_j: external.finish()?,
        })
    };
    let (increments, nonlinear) = solve_nonlinear(
        dimension,
        solve,
        |x| {
            let trial = evaluate(x)?;
            let mut residual: Vec<_> = trial
                .thermal
                .residual_j()
                .iter()
                .zip(thermal_scales_kg)
                .map(|(&r, &scale)| r / scale)
                .collect();
            residual.push(trial.gas_residual_j / gas_scale_kg);
            Ok(residual)
        },
        || observer.gate.is_requested(),
    )?;
    let trial = evaluate(&increments)?;
    let maximum_gas = trial.gas_residual_j.abs();
    if maximum_gas > solve.cell_energy_tolerance_j {
        return Err(DepositionError::EnergyBalance(trial.gas_residual_j).into());
    }
    // This union gate reads original-law liquid/solid/gas storage and separately
    // accumulated external heat, not a difference of opposite gas receipts.
    let combined = environment::finite(
        trial.thermal.audit().stored_energy_change_j + trial.gas_stored_j - trial.independent_heat_j,
    )?;
    if combined.abs() > environment::finite(total_tolerance_j + solve.cell_energy_tolerance_j)? {
        return Err(DepositionError::EnergyBalance(combined).into());
    }
    let maximum_cell = trial.thermal.maximum_cell_residual_j();
    let (next, payload, audit) =
        trial
            .thermal
            .accept(solve.cell_energy_tolerance_j, total_tolerance_j, check)?;
    Ok((next, payload, audit, nonlinear, maximum_cell, maximum_gas))
}
