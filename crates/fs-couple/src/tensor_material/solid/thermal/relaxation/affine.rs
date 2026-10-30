//! Affine stored tensor thermodynamics on the shared tensor/enthalpy driver.
use super::*;
use crate::tensor_material::TensorMaterialError;
use fs_material::visco_tensor::{TensorThermalState, source::thermal::TensorTemperatureDerivative};

/// Admission of source-operator changes internal to the solved thermal step.
#[derive(Clone, Copy, Debug)]
pub struct AffineThermalBudget {
    /// Maximum magnitude of native mechanical free-energy coefficient change
    /// per proposed interval [J]. This is not external supply or a heat tolerance.
    pub max_absolute_coefficient_work_j: f64,
    /// Native physical-memory reconstruction allowance, in (0, 1e-8].
    pub relative_memory_tolerance: f64,
}

/// Explicit per-element reversible exchange and material thermal potentials.
/// Temperatures and mechanical states are the SAME accepted physical endpoints.
#[derive(Clone, Debug)]
pub struct TensorThermalExchange {
    pub temperature_k: Vec<[f64; 2]>,
    pub before: Vec<TensorThermalState>,
    pub after: Vec<TensorThermalState>,
    /// psi(T_new, old strains) - psi(T_old, old strains), before motion [J].
    /// This is an INTERNAL thermal coefficient change, never external work.
    pub coefficient_work_j: Vec<f64>,
    /// T_new*S_new - T_old*S_old [J].
    pub entropy_energy_change_j: Vec<f64>,
    /// -coefficient_work - delta(T*S), supplied to the same caloric receiver.
    /// Signed reversible heat is distinct from positive dashpot dissipation.
    pub reversible_heat_j: Vec<f64>,
}

pub(super) fn states<E>(
    history: &TensorSolidHistory<'_, '_>,
    derivatives: &[TensorTemperatureDerivative],
    temperatures: &[f64],
    expansion_per_k: Option<&[[f64; 6]]>,
    tolerance: f64,
    gate: &CancelGate,
) -> Result<Vec<TensorThermalState>, Error<E>> {
    if derivatives.len() != history.solid.projection.element_count()
        || temperatures.len() != derivatives.len()
        || expansion_per_k.is_some_and(|a| a.len() != derivatives.len())
    {
        return Err(invalid(
            "complete thermal derivative and temperature per original element required",
        ));
    }
    let mut states = vec![TensorThermalState::default(); derivatives.len()];
    let mut at = 2 * history.inner.body.coordinates();
    for (point, &element) in history
        .solid
        .projection
        .material_points()
        .iter()
        .zip(history.solid.projection.point_elements())
    {
        check(gate)?;
        let count = 6 * point.branches().len();
        let state = point
            .thermal_state_with_expansion(
                &derivatives[element],
                history.inner.displacement_m(),
                &history.inner.state[at..at + count],
                temperatures[element],
                expansion_per_k.map_or([0.0; 6], |a| a[element]),
                tolerance,
                || gate.is_requested(),
            )
            .map_err(|e| solid(SolidError::Dynamics(TensorMaterialError::Material(e))))?;
        let total = &mut states[element];
        total.free_energy_j = finite(total.free_energy_j + state.free_energy_j)?;
        total.entropy_j_k = finite(total.entropy_j_k + state.entropy_j_k)?;
        total.internal_energy_j = finite(total.internal_energy_j + state.internal_energy_j)?;
        at += count;
    }
    check(gate)?;
    Ok(states)
}

pub(super) fn exchange<E>(
    before: &[TensorThermalState],
    frozen: &[TensorThermalState],
    after: &[TensorThermalState],
    temperatures_before: &[f64],
    temperatures_after: &[f64],
    native_coefficient_work_j: f64,
    energy_tolerance_j: f64,
    gate: &CancelGate,
) -> Result<TensorThermalExchange, Error<E>> {
    let mut result = TensorThermalExchange {
        temperature_k: Vec::with_capacity(before.len()),
        before: before.to_vec(),
        after: after.to_vec(),
        coefficient_work_j: Vec::with_capacity(before.len()),
        entropy_energy_change_j: Vec::with_capacity(before.len()),
        reversible_heat_j: Vec::with_capacity(before.len()),
    };
    let mut total = 0.0;
    for i in 0..before.len() {
        check(gate)?;
        let parameter = finite(frozen[i].free_energy_j - before[i].free_energy_j)?;
        let ts = finite(
            temperatures_after[i] * after[i].entropy_j_k
                - temperatures_before[i] * before[i].entropy_j_k,
        )?;
        let reversible = finite(-parameter - ts)?;
        result
            .temperature_k
            .push([temperatures_before[i], temperatures_after[i]]);
        result.coefficient_work_j.push(parameter);
        result.entropy_energy_change_j.push(ts);
        result.reversible_heat_j.push(reversible);
        total = finite(total + parameter)?;
    }
    if (total - native_coefficient_work_j).abs() > energy_tolerance_j {
        return Err(invalid(
            "material thermal coefficient energy disagrees with native source-operator transfer",
        ));
    }
    Ok(result)
}
