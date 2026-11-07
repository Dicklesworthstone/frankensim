//! Physical sensible storage through the original nonlinear tensor face owner.

use fs_exec::CancelGate;

use super::super::transient::prepare_tensor_step;
use super::{TensorConductivityLaw, TensorNonlinearConfig, TensorNonlinearSolution, solve_inner};
use crate::conjugate::energy::ThermalSetup;
use crate::conjugate::transient::TransientRecord;
use crate::conjugate::{ChtError, SolidMaterial, VoxelDomain, poll};

/// One accepted physical interval with the original-law endpoint and heat ledger.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorNonlinearStep {
    /// Actual integrated face heat, cell/face temperatures, and bounded
    /// nonlinear/linear solve evidence at the physical endpoint.
    pub solution: TensorNonlinearSolution,
    /// Original-state storage and independently integrated source/reservoir
    /// heat in joules. Krylov iterations include every Newton correction.
    pub record: TransientRecord,
}

/// Advance original-law full-tensor conduction by physical backward Euler.
///
/// `laws` retains one immutable original tensor conductivity law per cell.
/// Every trial resolves its actual cell-to-face conductivity integrals and
/// endpoint tensors through the same nonlinear spatial owner. Prescribed
/// boundary traces and finite contacts keep their existing physical heat
/// balance; no conductivity is frozen at the old or mean temperature.
///
/// `heat_capacity_j_k` is complete, strictly positive sensible capacity C in
/// J/K per cell, fixed over this interval. It is not multiplied by volume
/// again. The common transient preparation and ledger owner supplies C/dt
/// and accounts C*(T_new-T_old) from the immutable physical start state.
/// `previous_temperature_k` is also the initial nonlinear seed, never replaced
/// by a Newton trial or an intermediate face temperature. Positive storage
/// admits insulated bodies without an artificial temperature anchor.
///
/// The supplied setup is the physical endpoint source/boundary sample.
/// A prescribed cell's reservoir supplies its actual conductive reaction and
/// storage jump exactly once. Faces and zero-thickness contacts store no heat.
/// The returned native face heat contains conduction alone, while the local
/// cell acceptance includes its original storage term.
///
/// Backward Euler remains first order in time. This entry does not infer a
/// caloric or latent-heat law from conductivity data, update C with temperature,
/// or advance material history. The nonlinear owner's documented consistency
/// and convergence limits for arbitrary K(T) apply unchanged. Native source,
/// Newton, backtracking, shared Krylov and cancellation gates remain active.
///
/// # Errors
/// Shared time/storage admission, original conductivity support and spatial
/// ownership rules apply. Invalid or unrepresentable fields, exhausted work
/// budgets, unsupported original-law intervals, failure to converge, or
/// cancellation return no partial interval and leave every input unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_nonlinear_tensor_conduction(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    heat_capacity_j_k: &[f64],
    previous_temperature_k: &[f64],
    start_time_s: f64,
    time_step_s: f64,
    config: &TensorNonlinearConfig,
    gate: &CancelGate,
) -> Result<TensorNonlinearStep, ChtError> {
    let prepared = prepare_tensor_step(
        domain,
        laws.len(),
        heat_capacity_j_k,
        previous_temperature_k,
        start_time_s,
        time_step_s,
        gate,
    )?;
    let step = prepared.storage();
    let solution = solve_inner(
        domain,
        solids,
        setup,
        laws,
        previous_temperature_k,
        config,
        Some(&step),
        gate,
    )?;
    poll(gate)?;
    let record = prepared.record(domain, setup, &solution.solution)?;
    poll(gate)?;
    Ok(TensorNonlinearStep { solution, record })
}
