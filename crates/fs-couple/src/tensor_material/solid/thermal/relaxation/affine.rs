//! Source tensor thermodynamics on the shared tensor/enthalpy driver.
use super::*;
use crate::tensor_material::TensorMaterialError;
use fs_material::visco_tensor::{
    TensorThermalExpansion, TensorThermalState, source::thermal::TensorTemperatureDerivative,
};

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

/// Admit the complete caloric/stored-energy interval before any endpoint trial.
/// The source bound covers every original strain and retained-memory envelope,
/// including curvature with zero thermal expansion.
pub(super) fn admit_capacity<E>(
    history: &TensorSolidHistory<'_, '_>,
    thermal: &DepositionSimulation<'_>,
    cells: &[RelaxationThermalCell<'_>],
    derivatives: &[TensorTemperatureDerivative],
    temperatures_k: &[f64],
    expansion: Option<&[TensorThermalExpansion]>,
    gate: &CancelGate,
) -> Result<(), Error<E>> {
    if derivatives.len() != cells.len()
        || temperatures_k.len() != cells.len()
        || expansion.is_some_and(|a| a.len() != cells.len())
    {
        return Err(invalid(
            "complete thermal derivatives per original caloric receiver required",
        ));
    }
    let mut coefficients = vec![0.0; cells.len()];
    for (point, &element) in history
        .solid
        .projection
        .material_points()
        .iter()
        .zip(history.solid.projection.point_elements())
    {
        check(gate)?;
        let coefficient = if let Some(expansion) = expansion {
            point.thermal_capacity_bound_with_expansion_j_k2(
                &derivatives[element], temperatures_k[element], &expansion[element],
            )
        } else {
            point.thermal_capacity_bound_j_k2(
                &derivatives[element], temperatures_k[element], [0.0; 6],
            )
        }.map_err(|e| solid(SolidError::Dynamics(TensorMaterialError::Material(e))))?;
        coefficients[element] = finite(coefficients[element] + coefficient)?;
    }
    // The original caloric owner is piecewise linear in h(T). Bound
    // C_eff = mass*dh/dT - T*psi_TT on every intersecting segment, including
    // its hottest endpoint. Mechanical capacity stays in reversible exchange;
    // inserting it into the caloric curve would count the same energy twice.
    for (cell, coefficient) in cells.iter().zip(coefficients) {
        check(gate)?;
        let [low, high] = cell.temperature_interval_k;
        let mass = thermal
            .cells()
            .get(&cell.cell)
            .ok_or_else(|| invalid("tensor caloric receiver disappeared"))?
            .mass_kg;
        for segment in cell.caloric.curve().knots().windows(2) {
            let a = segment[0];
            let b = segment[1];
            if b.temperature_k < low || a.temperature_k > high {
                continue;
            }
            let capacity = mass * (b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg)
                / (b.temperature_k - a.temperature_k)
                - high.min(b.temperature_k) * coefficient;
            if !capacity.is_finite() || capacity <= 0.0 {
                return Err(invalid(
                    "tensor thermodynamics cannot establish positive total fixed-strain heat capacity over its admitted domain",
                ));
            }
        }
    }
    check(gate)
}

pub(super) fn states<E>(
    history: &TensorSolidHistory<'_, '_>,
    derivatives: &[TensorTemperatureDerivative],
    temperatures: &[f64],
    expansion: Option<&[TensorThermalExpansion]>,
    tolerance: f64,
    gate: &CancelGate,
) -> Result<Vec<TensorThermalState>, Error<E>> {
    if derivatives.len() != history.solid.projection.element_count()
        || temperatures.len() != derivatives.len()
        || expansion.is_some_and(|a| a.len() != derivatives.len())
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
        let state = if let Some(expansion) = expansion {
            point.thermal_state_with_expansion_derivatives(
                &derivatives[element], history.inner.displacement_m(),
                &history.inner.state[at..at + count], temperatures[element],
                &expansion[element], tolerance, || gate.is_requested(),
            )
        } else {
            point.thermal_state(
                &derivatives[element], history.inner.displacement_m(),
                &history.inner.state[at..at + count], temperatures[element],
                tolerance, || gate.is_requested(),
            )
        }.map_err(|e| solid(SolidError::Dynamics(TensorMaterialError::Material(e))))?;
        let total = &mut states[element];
        total.free_energy_j = finite(total.free_energy_j + state.free_energy_j)?;
        total.entropy_j_k = finite(total.entropy_j_k + state.entropy_j_k)?;
        total.internal_energy_j = finite(total.internal_energy_j + state.internal_energy_j)?;
        total.heat_capacity_j_k = finite(total.heat_capacity_j_k + state.heat_capacity_j_k)?;
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
