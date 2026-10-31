//! Closed original-to-final accounting for native finite gas transport.

use super::super::super::entropy;
use super::{ClosedGasParcel, DiffusionNetworkError, GasTransportError, GasTransportSolve};
use std::collections::BTreeMap;

/// Actual accepted endpoints, independent of the number of integration steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasTransportWindowAudit {
    pub total_mass_residual_kg: f64,
    /// Every source-qualified species has its own fixed absolute mass budget.
    pub maximum_species_mass_residual_kg: f64,
    /// Original NASA internal-energy difference with each conserved species'
    /// common formation reference removed before floating-point arithmetic.
    pub energy_change_j: f64,
    /// Source caloric entropy plus ideal mixing and actual mass redistribution.
    pub entropy_change_j_k: f64,
}

/// Audit one closed, rigid gas network directly from its original/final states.
///
/// Every parcel identity, volume and original thermochemical source is retained.
/// Chemical reactions, external mass/heat/work and moving boundaries need their
/// enclosing physical owner and are not admitted by this closed balance.
/// Species mass residuals remain separate; multiplying their rounding error by
/// an arbitrary standard-state energy or entropy datum would invent a defect.
/// All limits are the original absolute physical budgets, never step-count
/// scaled. This endpoint admission does not estimate temporal truncation error.
pub fn audit_closed_window(
    before: &BTreeMap<u64, ClosedGasParcel>,
    after: &BTreeMap<u64, ClosedGasParcel>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasTransportWindowAudit, GasTransportError> {
    if cancelled() {
        return Err(GasTransportError::Cancelled);
    }
    if before.is_empty()
        || !before.keys().eq(after.keys())
        || [
            controls.exchange.mass_tolerance_kg,
            controls.exchange.energy_tolerance_j,
            controls.entropy_tolerance_j_k,
        ]
        .into_iter()
        .any(|value| !value.is_finite() || value <= 0.0)
    {
        return Err(GasTransportError::Invalid(
            "closed gas window needs original nonempty parcel identities and positive finite physical budgets",
        ));
    }
    let change = entropy::closed_inventory_change(
        before.iter().map(|(&id, old)| (id, old, &after[&id])),
        controls.exchange.mass_tolerance_kg,
        &mut cancelled,
    )
    .map_err(|error| match error {
        entropy::Error::Invalid(reason) => {
            GasTransportError::Diffusion(DiffusionNetworkError::Invalid(reason))
        }
        entropy::Error::Parcel { id, error } => {
            GasTransportError::Diffusion(DiffusionNetworkError::Parcel { id, error })
        }
        entropy::Error::Cancelled => GasTransportError::Cancelled,
    })?;
    if change.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || change.energy_change_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(GasTransportError::Invalid(
            "original-to-final gas window exceeded its fixed mass or internal-energy budget",
        ));
    }
    if change.entropy_change_j_k < -controls.entropy_tolerance_j_k {
        return Err(GasTransportError::Diffusion(
            DiffusionNetworkError::EntropyDecrease {
                change_j_k: change.entropy_change_j_k,
                tolerance_j_k: controls.entropy_tolerance_j_k,
            },
        ));
    }
    if cancelled() {
        return Err(GasTransportError::Cancelled);
    }
    Ok(GasTransportWindowAudit {
        total_mass_residual_kg: change.total_mass_residual_kg,
        maximum_species_mass_residual_kg: change.maximum_species_mass_residual_kg,
        energy_change_j: change.energy_change_j,
        entropy_change_j_k: change.entropy_change_j_k,
    })
}
