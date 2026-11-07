//! Original gas calorics against a prescribed wall temperature and an integrated
//! solid/film law. Shared finite walls use this same port equation, then debit
//! their one retained inventory by the accepted port heats.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(crate) fn exchange_cell_with_reservoir_conduction(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    conduction: WallEndpointConduction<'_>,
    control: WallControl,
    cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    exchange_caloric_with_reservoir_conduction(
        &MixtureCaloric {
            gas,
            partial_densities: state.partial_densities,
        },
        state.momentum,
        state.energy,
        volume_m3,
        duration,
        contact,
        conduction,
        control,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn exchange_caloric_with_reservoir_conduction(
    gas: &impl WallGasCaloric,
    initial_momentum: f64,
    initial_energy: f64,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    conduction: WallEndpointConduction<'_>,
    control: WallControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    contact.validate()?;
    let WallThermal::Reservoir {
        temperature_k: wall_temperature,
        conductance_w_per_k: film,
    } = contact.thermal
    else {
        return Err(WallExchangeError::Invalid("integrated reservoir contact"));
    };
    let sensitivity = conduction
        .gas_sensitivity_bound_w_per_k
        .filter(|value| nonnegative(*value))
        .ok_or(WallExchangeError::Invalid(
            "integrated reservoir gas sensitivity",
        ))?;
    let integrated = sensitivity * duration.value();
    if !nonnegative(integrated) {
        return Err(WallExchangeError::Invalid(
            "integrated reservoir conductance",
        ));
    }
    // Mechanical work, stationary drag and slip heating remain in their original
    // owner. This second stage adds only one conservative thermal energy change.
    let mut update = exchange_caloric(
        gas,
        initial_momentum,
        initial_energy,
        volume_m3,
        duration,
        WallContact {
            drag: contact.drag,
            thermal: WallThermal::Adiabatic,
        },
        control,
        &mut cancelled,
    )?;
    update.contact = contact;
    let original = gas.primitive(update.momentum, update.energy, &mut cancelled)?;
    if duration.value() == 0.0 || film == 0.0 {
        // An insulating film still requires the retained material temperature
        // to remain in its original source, even when another shared port moves it.
        let response = (conduction.resolve)(
            contact,
            [wall_temperature, original.temperature],
            &mut cancelled,
        )?;
        response.validate()?;
        update.contact.thermal = WallThermal::Reservoir {
            temperature_k: wall_temperature,
            conductance_w_per_k: response.conductance_w_per_k,
        };
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let kinetic = 0.5 * update.momentum * original.velocity;
    let target = update.energy - kinetic;
    let initial_temperature = refine_caloric_temperature_with_allowance(
        gas,
        original.temperature,
        target,
        volume_m3,
        integrated,
        control.tolerance(0.0),
        control,
        &mut cancelled,
    )?;
    let domain = gas.domain();
    // At an outer shared-wall trial, (Twall_trial,Tgas_initial) is not a physical
    // endpoint and can have an unsupported material face. Start instead at the
    // gas-domain temperature nearest Twall. Its face lies between any supported
    // passive root face and Twall; no gas calorics outside their domain are used.
    let anchor = wall_temperature.clamp(domain.minimum_temperature, domain.maximum_temperature);
    let mut equation = |temperature: f64, checkpoint: &mut dyn FnMut() -> bool| {
        if checkpoint() {
            return Err(WallExchangeError::Cancelled);
        }
        let (internal_energy, cv) = gas.caloric(temperature)?;
        let heat = (internal_energy - target) * volume_m3;
        let response = (conduction.resolve)(contact, [wall_temperature, temperature], checkpoint)?;
        response.validate()?;
        let flux = duration.value() * response.power_to_gas_w;
        let capacity = cv * volume_m3;
        let slope = response
            .temperature_derivative_w_per_k
            .map(|derivatives| capacity - duration.value() * derivatives[1]);
        let residual = heat - flux;
        if !positive(capacity)
            || !residual.is_finite()
            || !flux.is_finite()
            || slope.is_some_and(|value| !positive(value))
        {
            return Err(WallExchangeError::Invalid(
                "integrated reservoir heat equation",
            ));
        }
        Ok((residual, slope, heat.abs().max(flux.abs())))
    };
    // Reserve the other half for original gas inversion and the material/film
    // face solve; final acceptance rechecks their actual represented endpoint.
    let mut solve_control = WallControl {
        absolute_energy_tolerance_j: 0.5 * control.absolute_energy_tolerance_j,
        relative_energy_tolerance: 0.5 * control.relative_energy_tolerance,
        ..control
    };
    let anchor_value = equation(anchor, &mut cancelled)?;
    let (temperature, iterations) =
        if anchor_value.0.abs() <= solve_control.tolerance(anchor_value.2) {
            (anchor, 0)
        } else {
            let direction = (initial_temperature - anchor).signum();
            if initial_temperature == anchor
                || direction * anchor_value.0 > solve_control.tolerance(anchor_value.2)
            {
                return Err(FlowError::EnergyDomain.into());
            }
            let (far, contractions) = supported_profile_heat_bracket(
                anchor,
                initial_temperature,
                solve_control,
                &mut equation,
                &mut cancelled,
            )?;
            solve_control.maximum_iterations -= contractions;
            let (temperature, iterations) = solve_scalar(
                anchor.min(far),
                anchor.max(far),
                anchor,
                solve_control,
                &mut equation,
                &mut cancelled,
            )?;
            (temperature, iterations + contractions)
        };
    let heat = (gas.caloric(temperature)?.0 - target) * volume_m3;
    let energy = update.energy + heat / volume_m3;
    let actual_heat = (energy - update.energy) * volume_m3;
    if !actual_heat.is_finite() || (heat != 0.0 && actual_heat == 0.0) {
        return Err(WallExchangeError::Invalid(
            "integrated reservoir heat is unrepresentable",
        ));
    }
    let actual = gas.primitive(update.momentum, energy, &mut cancelled)?;
    let actual_temperature = refine_caloric_temperature_with_allowance(
        gas,
        actual.temperature,
        energy - 0.5 * update.momentum * actual.velocity,
        volume_m3,
        integrated,
        control.tolerance(actual_heat.abs()),
        control,
        &mut cancelled,
    )?;
    let response = (conduction.resolve)(
        contact,
        [wall_temperature, actual_temperature],
        &mut cancelled,
    )?;
    response.validate()?;
    let flux = duration.value() * response.power_to_gas_w;
    let residual = actual_heat - flux;
    if !residual.is_finite()
        || !flux.is_finite()
        || residual.abs() + duration.value() * response.flux_error_w
            > control.tolerance(actual_heat.abs().max(flux.abs()))
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: update.report.energy_residual_j,
            thermal_residual_j: residual,
        });
    }
    update.energy = energy;
    update.contact.thermal = WallThermal::Reservoir {
        temperature_k: wall_temperature,
        conductance_w_per_k: response.conductance_w_per_k,
    };
    update.report.heat_to_gas_j = actual_heat;
    update.report.external_heat_to_gas_j = actual_heat;
    update.report.thermal_residual_j = residual;
    update.report.iterations = iterations;
    poll(&mut cancelled)?;
    Ok(update)
}
