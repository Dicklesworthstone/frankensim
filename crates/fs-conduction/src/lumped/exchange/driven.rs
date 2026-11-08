//! Constant applied heat and held convection through the original phase chart.
//!
//! The original uniform body owns its caloric inventory and full-chart Biot
//! admission. Sensible heat and latent residence share one enthalpy account.

use super::{
    compensated_sum, product_ratio,
};
use super::super::{BiotGate, LumpedEnthalpyBody, lumped_error, phase_state};
use crate::ConductionError;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::phase::{
    EquilibriumPhaseState,
    exchange::{
        PhaseDrivenHeatExchange, PhaseDrivenHeatExchangeJacobian, PhaseHeatExchangeError,
    },
};
use fs_math::eft::two_prod;

/// One held-convection and constant-power interval on the original body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedDrivenExchangeConfig {
    pub initial_specific_enthalpy_j_kg: f64,
    pub reservoir_temperature_k: f64,
    /// Signed heat deposited in the body [W]; positive heats, negative cools.
    pub applied_power_w: f64,
    pub duration_s: f64,
    /// Absolute tolerance for represented storage and separate heat accounts [J].
    pub energy_tolerance_j: f64,
}

/// Complete proposal; no accepted body, source debit or timeline is mutated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedDrivenExchange {
    pub body_identity: ContentHash,
    pub duration_s: f64,
    pub maximum_biot: f64,
    pub caloric: PhaseDrivenHeatExchange,
    pub heat_into_body_j: f64,
    pub stored_enthalpy_change_j: f64,
    pub heat_residual_j: f64,
    /// Original applied power times duration, counted once.
    pub source_heat_j: f64,
    /// Signed heat from the original convective reservoir into the body.
    pub reservoir_heat_j: f64,
    /// Source plus reservoir heat minus modeled net heat, compensated [J].
    pub heat_partition_residual_j: f64,
    /// Source plus reservoir heat minus represented native storage [J].
    pub energy_balance_residual_j: f64,
}

/// Continuous-flow sensitivities of the same admitted physical proposal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedDrivenExchangeResponse {
    pub exchange: LumpedDrivenExchange,
    pub caloric_jacobian: PhaseDrivenHeatExchangeJacobian,
    pub heat_initial_enthalpy_derivative_kg: f64,
    pub heat_reservoir_derivative_j_k: f64,
    pub heat_source_power_derivative_s: f64,
    pub heat_duration_derivative_w: f64,
}

/// Evolve signed applied heat through sensible branches and latent plateaus.
///
/// The original phase curve solves dh/ds = Q0 - X*(T(h)-T_initial) on 0..1.
/// Q0 retains original applied power and convection; X=G*dt/m. Neither a
/// positive virtual bath nor a finite apparent capacity is required during
/// fusion. Zero convection remains exact native-energy calorimetry.
///
/// The body's original mass, geometry, convection coefficient and source
/// context are held. Radiation must be explicitly disabled. Its complete
/// conductivity support and conservative Biot gate remain authoritative.
/// This proposes the original chart's enthalpy/phase evolution; it does not
/// advance free-surface motion, mechanical work or any other physical owner.
pub fn solve_lumped_driven_exchange(
    cx: &Cx<'_>,
    body: &LumpedEnthalpyBody<'_>,
    gate: BiotGate,
    config: LumpedDrivenExchangeConfig,
) -> Result<LumpedDrivenExchange, ConductionError> {
    let input = prepare(cx, body, gate, config)?;
    let caloric = body.phase_curve
        .driven_heat_exchange_with_checkpoint(
            input.initial, input.drive_j_kg, input.exposure_j_kg_k,
            || cx.checkpoint().is_err(),
        )
        .map_err(|error| material_error(body, error))?;
    finish(cx, body, config, input, caloric)
}

/// Same proposal with analytic derivatives through sensible and latent residence.
///
/// At a zero-drive caloric corner with unequal one-sided slopes,
/// a unique two-sided response does not exist. The response entry refuses that
/// case; the value-only entry still returns its physically unchanged state.
pub fn solve_lumped_driven_exchange_response(
    cx: &Cx<'_>,
    body: &LumpedEnthalpyBody<'_>,
    gate: BiotGate,
    config: LumpedDrivenExchangeConfig,
) -> Result<LumpedDrivenExchangeResponse, ConductionError> {
    let input = prepare(cx, body, gate, config)?;
    let response = body.phase_curve
        .driven_heat_exchange_response_with_checkpoint(
            input.initial, input.drive_j_kg, input.exposure_j_kg_k,
            || cx.checkpoint().is_err(),
        )
        .map_err(|error| material_error(body, error))?;
    let caloric_jacobian = response.jacobian;
    let exchange = finish(cx, body, config, input, response.exchange)?;
    let inlet = -body.mass_kg * caloric_jacobian.inlet_enthalpy_loss;
    let reservoir = body.mass_kg
        * (input.exposure_j_kg_k * caloric_jacobian.source_heat_gain);
    let source = config.duration_s * caloric_jacobian.source_heat_gain;
    let initial_bath_power = two_prod(
        input.conductance_w_k,
        config.reservoir_temperature_k - input.initial.temperature_k(),
    );
    let original_power = compensated_sum(&[
        initial_bath_power.0, config.applied_power_w, initial_bath_power.1,
    ]);
    let duration = if let Some(power) = original_power {
        power * caloric_jacobian.inlet_enthalpy_gain
    } else if config.duration_s > 0.0 {
        product_ratio(
            caloric_jacobian.remaining_drive_j_kg, body.mass_kg, config.duration_s,
        ).ok_or_else(|| invalid(body, "endpoint heat rate is not representable"))?
    } else {
        return Err(invalid(body, "initial heat rate is not representable"));
    };
    if [inlet, reservoir, source, duration].iter().any(|v| !v.is_finite()) {
        return Err(invalid(body, "driven phase heat Jacobian is not representable"));
    }
    checkpoint(cx)?;
    Ok(LumpedDrivenExchangeResponse {
        exchange,
        caloric_jacobian,
        heat_initial_enthalpy_derivative_kg: inlet,
        heat_reservoir_derivative_j_k: reservoir,
        heat_source_power_derivative_s: source,
        heat_duration_derivative_w: duration,
    })
}

#[derive(Clone, Copy)]
struct Inputs {
    initial: EquilibriumPhaseState,
    maximum_biot: f64,
    conductance_w_k: f64,
    exposure_j_kg_k: f64,
    drive_j_kg: f64,
    source_j_kg: f64,
    bath_initial: (f64, f64),
}

fn prepare(
    cx: &Cx<'_>,
    body: &LumpedEnthalpyBody<'_>,
    gate: BiotGate,
    config: LumpedDrivenExchangeConfig,
) -> Result<Inputs, ConductionError> {
    checkpoint(cx)?;
    if !config.initial_specific_enthalpy_j_kg.is_finite()
        || !config.reservoir_temperature_k.is_finite()
        || config.reservoir_temperature_k <= 0.0
        || !config.applied_power_w.is_finite()
        || !config.duration_s.is_finite() || config.duration_s < 0.0
        || !config.energy_tolerance_j.is_finite() || config.energy_tolerance_j < 0.0
    {
        return Err(invalid(body, "invalid driven phase state, source, duration or energy tolerance"));
    }
    if body.transport.maximum_emissivity() != 0.0 {
        return Err(invalid(body, "exact driven convection requires radiation to be explicitly disabled"));
    }
    // Transport caches the original full-chart minimum; zero radiation needs
    // no repeated traversal of caloric temperatures.
    let maximum_biot = body.convection_w_per_m2_k * body.characteristic_length_m
        / body.transport.minimum_conductivity_w_per_m_k();
    if !maximum_biot.is_finite() || maximum_biot > gate.ceiling()
        || (body.convection_w_per_m2_k > 0.0 && maximum_biot == 0.0)
    {
        return Err(lumped_error(
            &body.name,
            format!("convection Biot {maximum_biot} exceeds ceiling {} or is not representable; use spatial enthalpy transport", gate.ceiling()),
        ));
    }
    let conductance_w_k = body.convection_w_per_m2_k * body.surface_area_m2;
    if !conductance_w_k.is_finite()
        || (body.convection_w_per_m2_k > 0.0 && conductance_w_k == 0.0)
    {
        return Err(invalid(body, "original convective conductance is not representable"));
    }
    let initial = phase_state(body, config.initial_specific_enthalpy_j_kg)?;
    let exposure_j_kg_k = product_ratio(conductance_w_k, config.duration_s, body.mass_kg)
        .ok_or_else(|| invalid(body, "convective exposure is not representable"))?;
    let source_j_kg = product_ratio(config.applied_power_w, config.duration_s, body.mass_kg)
        .ok_or_else(|| invalid(body, "specific applied heat is not representable"))?;
    let gap = config.reservoir_temperature_k - initial.temperature_k();
    let bath_initial = two_prod(exposure_j_kg_k, gap);
    if exposure_j_kg_k != 0.0 && gap != 0.0 && bath_initial.0 == 0.0 {
        return Err(invalid(body, "integrated convective drive is not representable"));
    }
    let drive_j_kg = compensated_sum(&[bath_initial.0, source_j_kg, bath_initial.1])
        .ok_or_else(|| invalid(body, "integrated source and bath forcing is not representable"))?;
    Ok(Inputs {
        initial, maximum_biot, conductance_w_k, exposure_j_kg_k, drive_j_kg,
        source_j_kg, bath_initial,
    })
}

fn finish(
    cx: &Cx<'_>,
    body: &LumpedEnthalpyBody<'_>,
    config: LumpedDrivenExchangeConfig,
    input: Inputs,
    caloric: PhaseDrivenHeatExchange,
) -> Result<LumpedDrivenExchange, ConductionError> {
    let heat = body.mass_kg * caloric.specific_heat_j_kg;
    let stored = body.mass_kg * caloric.represented_enthalpy_change_j_kg;
    let residual = body.mass_kg * caloric.heat_residual_j_kg;
    if !heat.is_finite() || !stored.is_finite() || !residual.is_finite()
        || (caloric.specific_heat_j_kg != 0.0 && heat == 0.0)
        || (caloric.represented_enthalpy_change_j_kg != 0.0 && stored == 0.0)
        || residual.abs() > config.energy_tolerance_j
        || (heat - stored).abs() > config.energy_tolerance_j
    {
        return Err(invalid(body, "driven phase heat exceeds the represented native-energy budget"));
    }
    // Relative mean travel retains weak convective loss during source heating.
    // Net heat minus source instead retains finite bath heat at saturation.
    let direct_bath = (-input.exposure_j_kg_k)
        .mul_add(caloric.mean_temperature_change_k, input.bath_initial.0)
        + input.bath_initial.1;
    let balanced_bath = caloric.specific_heat_j_kg - input.source_j_kg;
    let direct_scale = input.bath_initial.0.abs() + input.bath_initial.1.abs()
        + (input.exposure_j_kg_k * caloric.mean_temperature_change_k).abs();
    let balanced_scale = caloric.specific_heat_j_kg.abs() + input.source_j_kg.abs();
    let bath_specific = if input.exposure_j_kg_k > 0.0
        && balanced_bath.is_finite() && balanced_scale < direct_scale
    {
        balanced_bath
    } else {
        direct_bath
    };
    let reservoir = body.mass_kg * bath_specific;
    let source = config.applied_power_w * config.duration_s;
    if !reservoir.is_finite() || !source.is_finite()
        || (bath_specific != 0.0 && reservoir == 0.0)
        || (config.applied_power_w != 0.0 && config.duration_s > 0.0 && source == 0.0)
    {
        return Err(invalid(body, "individual source or reservoir heat is not representable"));
    }
    let partition = compensated_sum(&[source, reservoir, -heat])
        .ok_or_else(|| invalid(body, "source and bath heat summation is not representable"))?;
    let balance = compensated_sum(&[source, reservoir, -stored])
        .ok_or_else(|| invalid(body, "native enthalpy balance is not representable"))?;
    if partition.abs() > config.energy_tolerance_j || balance.abs() > config.energy_tolerance_j {
        return Err(lumped_error(
            &body.name,
            format!("source/bath partition {partition} J or native balance {balance} J exceeds {} J tolerance", config.energy_tolerance_j),
        ));
    }
    checkpoint(cx)?;
    Ok(LumpedDrivenExchange {
        body_identity: body.identity(),
        duration_s: config.duration_s,
        maximum_biot: input.maximum_biot,
        caloric,
        heat_into_body_j: heat,
        stored_enthalpy_change_j: stored,
        heat_residual_j: residual,
        source_heat_j: source,
        reservoir_heat_j: reservoir,
        heat_partition_residual_j: partition,
        energy_balance_residual_j: balance,
    })
}

fn material_error(body: &LumpedEnthalpyBody<'_>, error: PhaseHeatExchangeError) -> ConductionError {
    match error {
        PhaseHeatExchangeError::Cancelled => ConductionError::Cancelled {
            stage: "lumped-driven-exchange", at: 0,
        },
        error => lumped_error(&body.name, error.to_string()),
    }
}

fn invalid(body: &LumpedEnthalpyBody<'_>, what: &str) -> ConductionError {
    lumped_error(&body.name, what.to_owned())
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "lumped-driven-exchange", at: 0,
    })
}
