//! Direct fixed-resistance contact between two finite phase-changing bodies.
//!
//! This is one isolated reciprocal contact with fixed mass and geometry. It
//! introduces no applied heat, ambient bath or radiation channel. The original
//! series resistance supplies G=1/R, and each body's full-source conductivity
//! minimum retains the existing contact Biot admission.
//!
//! Heat conserves the pair's original declared enthalpy account. The caller
//! retains the charts' pressure/work convention: density/volume change is not
//! an additional pressure-work source or a constant-volume internal-energy claim.

use super::super::{LUMPED_BIOT_CEILING, LumpedEnthalpyBody, lumped_error, phase_state};
use crate::{ConductionError, interface::SeriesThermalResistance};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::phase::exchange::pair::{
    FinitePhaseBody, FinitePhaseExchange, FinitePhaseExchangeControls, FinitePhaseExchangeError,
    exchange_with_checkpoint,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedPairExchangeConfig {
    pub initial_specific_enthalpies_j_kg: [f64; 2],
    pub duration_s: f64,
    /// Absolute J gate for EACH represented stored-heat residual and their sum.
    /// Independent from the material trajectory heat/exposure error budget.
    pub energy_tolerance_j: f64,
    pub integration: FinitePhaseExchangeControls,
}

/// Both states and their heat/entropy accounts; no accepted history is modified.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedPairExchange {
    pub body_identities: [ContentHash; 2],
    pub duration_s: f64,
    pub resistance_k_per_w: f64,
    pub maximum_biot: f64,
    pub caloric: FinitePhaseExchange,
    pub total_energy_residual_j: f64,
}

/// Integrate the isolated fixed-R trajectory through both original phase charts.
///
/// Constant-Cp, authored affine T(h) and latent intervals use their analytic
/// trajectories. Varying original Cp uses bounded source-split quadrature and
/// inversion. There is no backward-Euler endpoint-power approximation or
/// prescribed phase-front timing. Equilibrium is never crossed.
///
/// Each declared surface area is the effective contact area for the original
/// Lc/(R*A*k_min) reduced Biot check; this does not solve spatial spreading.
/// Source exhaustion, lost absolute-h increments, either local energy residual,
/// the combined energy residual and cancellation can all refuse the proposal.
pub fn solve_lumped_pair_exchange(
    cx: &Cx<'_>,
    first: &LumpedEnthalpyBody<'_>,
    second: &LumpedEnthalpyBody<'_>,
    resistance: &SeriesThermalResistance,
    config: LumpedPairExchangeConfig,
) -> Result<LumpedPairExchange, ConductionError> {
    let fail = |what: &str| lumped_error("finite-body direct contact", what.to_owned());
    checkpoint(cx)?;
    if !config.duration_s.is_finite()
        || config.duration_s < 0.0
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j < 0.0
    {
        return Err(fail(
            "duration and absolute energy tolerance must be finite/nonnegative",
        ));
    }
    let r = resistance.budget().value_k_per_w;
    let exposure = config.duration_s / r;
    if !exposure.is_finite() || (config.duration_s > 0.0 && exposure == 0.0) {
        return Err(fail("fixed-resistance exposure is not representable"));
    }
    let mut maximum_biot = 0.0_f64;
    for body in [first, second] {
        let biot = body.characteristic_length_m
            / r
            / body.surface_area_m2
            / body.transport.minimum_conductivity_w_per_m_k;
        if !biot.is_finite() || biot > LUMPED_BIOT_CEILING {
            return Err(fail(
                "contact Biot exceeds the lumped limit; use spatial enthalpy transport",
            ));
        }
        maximum_biot = maximum_biot.max(biot);
    }
    let initial_first = phase_state(first, config.initial_specific_enthalpies_j_kg[0])?;
    let initial_second = phase_state(second, config.initial_specific_enthalpies_j_kg[1])?;
    let caloric = exchange_with_checkpoint(
        FinitePhaseBody {
            curve: first.phase_curve,
            initial: initial_first,
            mass_kg: first.mass_kg,
        },
        FinitePhaseBody {
            curve: second.phase_curve,
            initial: initial_second,
            mass_kg: second.mass_kg,
        },
        exposure,
        config.integration,
        || cx.checkpoint().is_err(),
    )
    .map_err(|error| match error {
        FinitePhaseExchangeError::Cancelled => ConductionError::Cancelled {
            stage: "lumped-pair-exchange",
            at: 0,
        },
        error => fail(&error.to_string()),
    })?;
    let total_residual = caloric.energy_residual_j[0] + caloric.energy_residual_j[1];
    if !total_residual.is_finite()
        || total_residual.abs() > config.energy_tolerance_j
        || caloric
            .energy_residual_j
            .iter()
            .any(|residual| !residual.is_finite() || residual.abs() > config.energy_tolerance_j)
    {
        return Err(fail(
            "represented finite-body stored heat exceeds the local or total energy budget",
        ));
    }
    checkpoint(cx)?;
    Ok(LumpedPairExchange {
        body_identities: [first.identity(), second.identity()],
        duration_s: config.duration_s,
        resistance_k_per_w: r,
        maximum_biot,
        caloric,
        total_energy_residual_j: total_residual,
    })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "lumped-pair-exchange",
        at: 0,
    })
}
