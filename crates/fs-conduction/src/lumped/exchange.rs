//! Exact fixed-bath convection through the original uniform-body phase curve.
//!
//! The existing body owns mass, area, conductivity and calorics. Its conservative
//! Biot gate remains mandatory. This entry requires radiation to be disabled;
//! it adds no internal source or second heat-transfer model.

use super::{BiotGate, LumpedEnthalpyBody, lumped_error, phase_state};
use crate::ConductionError;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::phase::exchange::{PhaseHeatExchange, PhaseHeatExchangeError};

/// One held-bath exchange. The body's original convection coefficient is held.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedLinearExchangeConfig {
    pub initial_specific_enthalpy_j_kg: f64,
    pub reservoir_temperature_k: f64,
    pub duration_s: f64,
    /// Absolute allowed difference between integrated heat and stored h [J].
    /// This admits floating-point state quantization, not a temporal error.
    pub energy_tolerance_j: f64,
}

/// Complete proposal, without modifying the body or any accepted history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedLinearExchange {
    pub body_identity: ContentHash,
    pub duration_s: f64,
    pub maximum_biot: f64,
    pub caloric: PhaseHeatExchange,
    pub heat_into_body_j: f64,
    pub stored_enthalpy_change_j: f64,
    pub heat_residual_j: f64,
}

/// Evolve a radiation-disabled body against one held-temperature fluid bath.
///
/// The original phase owner integrates each sensible branch and latent plateau
/// directly. No timestep, apparent heat capacity or endpoint-power quadrature
/// is introduced. The resulting mean temperature closes G*dt*(T_bath-T_mean).
/// With zero conductance or duration the state is preserved.
///
/// Existing source support, material identities and the conservative full-span
/// conductivity/Biot admission remain load bearing. The body's geometry and
/// conductance are held over the exchange; density and phase come from its
/// original chart. The caller retains that chart's pressure/work convention.
/// Nonzero emissivity refuses instead of silently omitting physical radiation.
pub fn solve_lumped_linear_exchange(
    cx: &Cx<'_>,
    body: &LumpedEnthalpyBody<'_>,
    gate: BiotGate,
    config: LumpedLinearExchangeConfig,
) -> Result<LumpedLinearExchange, ConductionError> {
    checkpoint(cx)?;
    if !config.initial_specific_enthalpy_j_kg.is_finite()
        || !config.reservoir_temperature_k.is_finite()
        || config.reservoir_temperature_k <= 0.0
        || !config.duration_s.is_finite()
        || config.duration_s < 0.0
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j < 0.0
    {
        return Err(lumped_error(
            &body.name,
            "invalid fixed-bath exchange state, duration or energy tolerance".into(),
        ));
    }
    if body.transport.maximum_emissivity != 0.0 {
        return Err(lumped_error(
            &body.name,
            "exact fixed-bath convection requires radiation to be explicitly disabled".into(),
        ));
    }
    let maximum_biot = body.maximum_biot(config.reservoir_temperature_k);
    if !maximum_biot.is_finite() || maximum_biot > gate.ceiling() {
        return Err(lumped_error(
            &body.name,
            format!(
                "convection Biot {maximum_biot} exceeds ceiling {}; escalate to spatial enthalpy transport",
                gate.ceiling(),
            ),
        ));
    }
    let conductance = body.convection_w_per_m2_k * body.surface_area_m2;
    let exposure = (conductance / body.mass_kg) * config.duration_s;
    if !conductance.is_finite()
        || !exposure.is_finite()
        || (body.convection_w_per_m2_k > 0.0 && conductance == 0.0)
        || (conductance > 0.0 && config.duration_s > 0.0 && exposure == 0.0)
    {
        return Err(lumped_error(
            &body.name,
            "fixed-bath exchange exposure is not representable".into(),
        ));
    }
    let initial = phase_state(body, config.initial_specific_enthalpy_j_kg)?;
    let caloric = body
        .phase_curve
        .linear_heat_exchange_with_checkpoint(
            initial,
            config.reservoir_temperature_k,
            exposure,
            || cx.checkpoint().is_err(),
        )
        .map_err(|error| match error {
            PhaseHeatExchangeError::Cancelled => ConductionError::Cancelled {
                stage: "lumped-linear-exchange",
                at: 0,
            },
            error => lumped_error(&body.name, error.to_string()),
        })?;
    let heat = body.mass_kg * caloric.specific_heat_j_kg;
    let stored = body.mass_kg * caloric.represented_enthalpy_change_j_kg;
    let residual = body.mass_kg * caloric.heat_residual_j_kg;
    if !heat.is_finite()
        || !stored.is_finite()
        || !residual.is_finite()
        || (caloric.specific_heat_j_kg != 0.0 && heat == 0.0)
        || (caloric.represented_enthalpy_change_j_kg != 0.0 && stored == 0.0)
        || residual.abs() > config.energy_tolerance_j
    {
        return Err(lumped_error(
            &body.name,
            format!(
                "fixed-bath stored-energy residual {residual} J exceeds the admitted representation or {} J tolerance",
                config.energy_tolerance_j,
            ),
        ));
    }
    checkpoint(cx)?;
    Ok(LumpedLinearExchange {
        body_identity: body.identity(),
        duration_s: config.duration_s,
        maximum_biot,
        caloric,
        heat_into_body_j: heat,
        stored_enthalpy_change_j: stored,
        heat_residual_j: residual,
    })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "lumped-linear-exchange",
        at: 0,
    })
}
