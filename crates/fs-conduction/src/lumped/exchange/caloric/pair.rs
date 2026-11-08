//! Reciprocal fixed-resistance exchange on two original native caloric cells.
//!
//! One signed heat transfer closes both declared storage accounts. Each body
//! retains its own Cp/h, Cv/u, reference-volume or amount-bound sensible law.
//! The original pressure/work conventions remain explicit; summing different
//! storage coordinates is not an assertion that both are internal energy.

use super::{LumpedCaloricExchangeConfig, prepare_caloric_exchange};
use super::super::compensated_sum;
use crate::{ConductionError, interface::SeriesThermalResistance};
use crate::lumped::{BiotGate, lumped_error};
use fs_exec::Cx;
use fs_material::phase::exchange::pair::FinitePhaseExchangeError;
use fs_material::thermal::transport::conduction::{
    CaloricConductionLayer,
    pair::{
        FiniteCaloricBody, FiniteCaloricExchange, FiniteCaloricExchangeControls,
        FiniteCaloricExchangeError, exchange_with_checkpoint,
    },
};
use fs_thermal::enthalpy::DepositedCell;

/// Held geometry, original resistance and independent numerical acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedCaloricPairExchangeConfig {
    pub duration_s: f64,
    /// Internal conduction distance for each body, in first/second order [m].
    pub characteristic_lengths_m: [f64; 2],
    /// Actual reduced exchange area on each body [m2], within its surface area.
    pub contact_areas_m2: [f64; 2],
    /// Absolute J gate for EACH stored-heat residual and their combined sum.
    pub energy_tolerance_j: f64,
    /// Original pair trajectory's separate heat/exposure and work budgets.
    pub integration: FiniteCaloricExchangeControls,
}

/// Atomic two-cell proposal, in supplied first/second order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedCaloricPairExchange {
    pub cell_ids: [u64; 2],
    pub caloric_identities: [[u8; 32]; 2],
    pub duration_s: f64,
    pub resistance_k_per_w: f64,
    pub biot_numbers: [f64; 2],
    pub caloric: FiniteCaloricExchange,
    pub total_energy_residual_j: f64,
}

/// Exchange heat through one isolated original fixed series resistance.
///
/// The shared finite-body driver integrates both original variable capacities
/// along the same reciprocal trajectory. Constant capacities use the analytic
/// relaxation; nonlinear source intervals retain bounded quadrature/inversion.
/// No held-temperature surrogate body or resampled caloric chart is introduced.
///
/// Original cell inventories, source limits and conductivity support apply at
/// both endpoints. The contact Biot number Lc/(R*A*k_min) uses each declared
/// contact area and the minimum over EVERY original conductivity knot. A
/// selected scalar conductivity does not establish unprovided isotropy or
/// spatial spreading. The caller owns the admitted geometry and work context.
///
/// Both cells must have zero active ambient/support conductance and zero
/// applied power for this isolated exchange. Other physical owners need their
/// own coupled accounting window. Neither input cell nor a source is mutated.
#[allow(clippy::too_many_arguments)]
pub fn solve_lumped_caloric_pair_exchange(
    cx: &Cx<'_>,
    first: &DepositedCell,
    first_layer: &CaloricConductionLayer<'_, '_>,
    second: &DepositedCell,
    second_layer: &CaloricConductionLayer<'_, '_>,
    resistance: &SeriesThermalResistance,
    gate: BiotGate,
    config: LumpedCaloricPairExchangeConfig,
) -> Result<LumpedCaloricPairExchange, ConductionError> {
    checkpoint(cx)?;
    if first.id == second.id {
        return Err(invalid("a finite contact requires two distinct cell identities"));
    }
    let r = resistance.budget().value_k_per_w;
    if !r.is_finite() || r <= 0.0 {
        return Err(invalid("fixed series resistance must be finite and positive"));
    }
    let (first_identity, first_biot) =
        admit_cell(cx, first, first_layer, gate, config, r, 0)?;
    let (second_identity, second_biot) =
        admit_cell(cx, second, second_layer, gate, config, r, 1)?;
    let exposure = config.duration_s / r;
    if !exposure.is_finite() || (config.duration_s > 0.0 && exposure == 0.0) {
        return Err(invalid("fixed-resistance conductance-time is not representable"));
    }
    let caloric = exchange_with_checkpoint(
        FiniteCaloricBody {
            layer: first_layer,
            initial_specific_energy_j_kg: first.specific_h_j_kg,
            mass_kg: first.mass_kg,
        },
        FiniteCaloricBody {
            layer: second_layer,
            initial_specific_energy_j_kg: second.specific_h_j_kg,
            mass_kg: second.mass_kg,
        },
        exposure,
        config.integration,
        || cx.checkpoint().is_err(),
    ).map_err(|error| match error {
        FiniteCaloricExchangeError::Integration(FinitePhaseExchangeError::Cancelled) =>
            ConductionError::Cancelled { stage: "lumped-caloric-pair", at: 0 },
        error => invalid(&error.to_string()),
    })?;

    admit_endpoint(first, first_layer, caloric.final_specific_energies_j_kg[0])?;
    admit_endpoint(second, second_layer, caloric.final_specific_energies_j_kg[1])?;
    let total_residual = compensated_sum(&caloric.energy_residual_j)
        .ok_or_else(|| invalid("combined native stored-heat residual is not representable"))?;
    if total_residual.abs() > config.energy_tolerance_j
        || caloric.energy_residual_j.iter().any(|residual| {
            !residual.is_finite() || residual.abs() > config.energy_tolerance_j
        })
    {
        return Err(invalid("finite-body stored heat exceeds a local or combined energy budget"));
    }
    checkpoint(cx)?;
    Ok(LumpedCaloricPairExchange {
        cell_ids: [first.id.0, second.id.0],
        caloric_identities: [first_identity, second_identity],
        duration_s: config.duration_s,
        resistance_k_per_w: r,
        biot_numbers: [first_biot, second_biot],
        caloric,
        total_energy_residual_j: total_residual,
    })
}

#[allow(clippy::too_many_arguments)]
fn admit_cell(
    cx: &Cx<'_>,
    cell: &DepositedCell,
    layer: &CaloricConductionLayer<'_, '_>,
    gate: BiotGate,
    config: LumpedCaloricPairExchangeConfig,
    resistance: f64,
    side: usize,
) -> Result<([u8; 32], f64), ConductionError> {
    let prepared = prepare_caloric_exchange(
        cx, cell, layer, gate,
        LumpedCaloricExchangeConfig {
            duration_s: config.duration_s,
            characteristic_length_m: config.characteristic_lengths_m[side],
            energy_tolerance_j: config.energy_tolerance_j,
        },
    )?;
    if prepared.total_conductance_w_k != 0.0 || cell.boundary.source_w != 0.0 {
        return Err(invalid("isolated finite contact requires zero ambient/support conductance and applied source"));
    }
    let area = config.contact_areas_m2[side];
    if !area.is_finite() || area <= 0.0 || area > cell.geometry.surface_area_m2 {
        return Err(invalid("contact area must be positive and lie within the original cell surface"));
    }
    let biot = config.characteristic_lengths_m[side] / resistance / area
        / prepared.minimum_conductivity_w_m_k;
    if !biot.is_finite() || biot <= 0.0 || biot > gate.ceiling() {
        return Err(lumped_error(
            "native caloric finite contact",
            format!("side {side} contact Biot {biot} exceeds ceiling {} or is not representable; use spatial caloric transport", gate.ceiling()),
        ));
    }
    Ok((prepared.caloric_identity, biot))
}

fn admit_endpoint(
    original: &DepositedCell,
    layer: &CaloricConductionLayer<'_, '_>,
    final_energy: f64,
) -> Result<(), ConductionError> {
    let mut endpoint = original.clone();
    endpoint.specific_h_j_kg = final_energy;
    layer.caloric_law().validate_cell(&endpoint).map_err(|error| {
        invalid(&format!("original inventory refuses the proposed finite-body state: {error:?}"))
    })?;
    layer.state_at_specific_energy(final_energy).map_err(|error| {
        invalid(&format!("proposed finite-body state exceeds caloric/conductivity support: {error:?}"))
    })?;
    Ok(())
}

fn invalid(what: &str) -> ConductionError {
    lumped_error("native caloric finite contact", what.to_owned())
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "lumped-caloric-pair", at: 0,
    })
}
