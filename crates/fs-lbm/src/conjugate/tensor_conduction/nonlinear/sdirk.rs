//! Second-order physical time integration with original nonlinear tensor laws.

use fs_exec::CancelGate;

use super::super::sdirk::{interval_record, stage_carry, stage_clock};
use super::{
    TensorConductivityLaw, TensorNonlinearConfig, TensorNonlinearSolution,
    advance_nonlinear_tensor_conduction, solve_inner,
};
use crate::conjugate::energy::{EnergyConfig, PseudoStep, ThermalSetup};
use crate::conjugate::tensor_conduction::invalid;
use crate::conjugate::transient::TransientRecord;
use crate::conjugate::{ChtError, SolidMaterial, VoxelDomain, poll};

/// One accepted SDIRK2 interval with both original-law stages retained.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorNonlinearSdirk2Step {
    /// Physical source/boundary sample times: start + gamma*dt, then start + dt.
    pub stage_time_s: [f64; 2],
    /// First stage and stiffly accurate endpoint, each with actual integrated
    /// tensor heat, both contact traces, and its nonlinear solve evidence.
    pub stages: [TensorNonlinearSolution; 2],
    /// Whole-interval storage and tableau-weighted external heat in joules.
    /// Iterations sum all Krylov work in both nonlinear stages.
    pub record: TransientRecord,
}

/// Advance original-law full-tensor conduction with physical SDIRK2 stages.
///
/// Uses the same two-stage tableau and physical heat ledger as
/// [`super::super::advance_tensor_conduction_sdirk2`], with
/// gamma = 1 - 1/sqrt(2). Every stage reevaluates the original conductivity
/// integrals and endpoint tensors throughout its nonlinear solve. Conductivity
/// is never frozen at the start or a stage-average temperature. The method is
/// second order for smooth forcing and laws with fixed constraints; its
/// tableau is L-stable. This does not establish unconditional nonlinear
/// convergence, positivity for arbitrary time steps, or a maximum principle
/// for full-tensor anisotropy.
///
/// `setups` are explicit physical source, boundary, contact and sink samples
/// at the two returned stage times. Supply both samples for time-dependent
/// forcing, and split intervals at discontinuities. Both stages use the same
/// immutable laws and strictly positive sensible cell capacities in J/K.
/// Capacity is not multiplied by volume again. This entry does not infer
/// latent heat or temperature-dependent storage.
///
/// Both storage equations retain the physical `previous_temperature_k`.
/// Stage 2 uses the accepted first-stage cell temperatures only as a nonlinear
/// seed, and carries its actual physical heat as an internal rate. Prescribed
/// cells include their reservoir rate and storage exactly once, including
/// when the prescribed mask changes. The carry is never an external source.
///
/// All three work limits in `config` are shared over the complete interval:
/// Krylov iterations, Newton corrections, and rejected line-search candidates.
/// Either stage's refusal or cancellation returns no partial interval and
/// changes no input.
///
/// # Errors
/// The original-law spatial and sensible-storage admission rules apply at
/// both stages. Unsupported original source intervals, invalid fields or
/// unresolved stage times, unrepresentable heat/storage, exhausted shared
/// budgets, failed nonlinear convergence, allocation refusal or cancellation
/// return an error.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_nonlinear_tensor_conduction_sdirk2(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setups: [&ThermalSetup; 2],
    laws: &[&dyn TensorConductivityLaw],
    heat_capacity_j_k: &[f64],
    previous_temperature_k: &[f64],
    start_time_s: f64,
    time_step_s: f64,
    config: &TensorNonlinearConfig,
    gate: &CancelGate,
) -> Result<TensorNonlinearSdirk2Step, ChtError> {
    poll(gate)?;
    let (stage_interval_s, stage_time_s) =
        stage_clock(start_time_s, time_step_s).ok_or_else(|| {
            invalid("nonlinear SDIRK2 needs finite, nonnegative and strictly ordered stage times")
        })?;
    let first = advance_nonlinear_tensor_conduction(
        domain,
        solids,
        setups[0],
        laws,
        heat_capacity_j_k,
        previous_temperature_k,
        start_time_s,
        stage_interval_s,
        config,
        gate,
    )?;
    let (coefficient_w_k, carried_power_w) = stage_carry(
        domain,
        setups[0],
        &first.solution.solution,
        heat_capacity_j_k,
        previous_temperature_k,
        stage_interval_s,
        gate,
    )?;
    let remaining = |limit: usize, used: usize| {
        limit
            .checked_sub(used)
            .ok_or_else(|| invalid("first nonlinear SDIRK2 stage exceeded its shared work budget"))
    };
    let second_config = TensorNonlinearConfig {
        linear: EnergyConfig {
            max_iterations: remaining(
                config.linear.max_iterations,
                first.solution.solution.report.iterations,
            )?,
            ..config.linear
        },
        max_newton_iterations: remaining(
            config.max_newton_iterations,
            first.solution.newton_iterations,
        )?,
        max_backtracks: remaining(config.max_backtracks, first.solution.backtracks)?,
    };
    let storage_step = PseudoStep {
        coefficient_w_k: &coefficient_w_k,
        previous: previous_temperature_k,
        junction_coefficient_w_k: &[],
    };
    let second = solve_inner(
        domain,
        solids,
        setups[1],
        laws,
        &first.solution.solution.temperature,
        &second_config,
        Some(&storage_step),
        Some(&carried_power_w),
        gate,
    )?;
    poll(gate)?;
    let iterations = first
        .solution
        .solution
        .report
        .iterations
        .checked_add(second.solution.report.iterations)
        .filter(|&count| count <= config.linear.max_iterations)
        .ok_or_else(|| invalid("nonlinear SDIRK2 stages exceeded the shared Krylov budget"))?;
    let record = interval_record(
        domain,
        setups[1],
        &first.record,
        &second.solution,
        heat_capacity_j_k,
        previous_temperature_k,
        &carried_power_w,
        stage_time_s[1],
        stage_interval_s,
        iterations,
        gate,
    )?;
    poll(gate)?;
    Ok(TensorNonlinearSdirk2Step {
        stage_time_s,
        stages: [first.solution, second],
        record,
    })
}
