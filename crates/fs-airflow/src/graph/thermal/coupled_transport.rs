//! Shared-solid fixed point with upstream heating and junction mixing.
//!
//! One callback solves all solid Robin regions per iteration. The graph then
//! transports the resulting heat downstream, including split/merge junctions.
//! Existing single-path probes enforce region ownership, applied references,
//! area, temperature and heat-rate consistency; they never rerun the solid.
//! Hydraulics remain frozen, and fixed-point convergence is not validation.
//!
//! The explicit IQN-ILS entry points retain vector secants across the ENTIRE
//! network, including feedback through mixing junctions and the common solid.
//! They use the same unrelaxed temperature and branch-local watt gates as the
//! fixed/Aitken entry points. Acceleration never supplies a stopping criterion.

use fs_couple::AitkenRelaxation;
use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig};
use fs_exec::Cx;

use crate::AirflowError;
use crate::conjugate::{ConjugateConfig, RegionBalance, Relaxation, SolidRegionState, solve_conjugate_from};
use super::transport::{TransportError, TransportMarch, TransportNetwork};
use super::{admit, checked, checked_sum, iqn_error};

/// A mutually consistent solid response and transported air state.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledTransportSolution {
    /// Exact reference vector used by the final solid callback, in
    /// `TransportNetwork::regions()` order. The transport's proposed vector
    /// is separate and differs by at most the configured temperature tolerance.
    pub reference_temperatures_k: Vec<f64>,
    /// Final solid response, evaluated against the exact vector above.
    pub solid: Vec<SolidRegionState>,
    /// Air state recomputed from those final wall temperatures.
    pub transport: TransportMarch,
    /// Admitted per-region solid/air heat rates, in the same region order.
    pub region_balances: Vec<RegionBalance>,
    /// Shared solid callbacks performed, not summed branch-probe iterations.
    pub iterations: usize,
    /// Maximum unrelaxed reference residual of the accepted exchange, K.
    pub max_reference_change_k: f64,
}

/// Solve wall temperatures and air mixing together under a frozen hydraulic
/// solution. The callback accepts every Robin reference in the network's region
/// order and returns every region in that same order. An all-adiabatic network
/// should use `TransportNetwork::march` instead and refuses here.
///
/// # Errors
/// Preserves typed channel/solid refusals, rejects invalid coupling budgets
/// before solid work, and returns no partially converged result. Cancellation
/// from any nested operation retains the complete current reference vector.
pub fn solve_coupled_transport<F>(
    cx: &Cx<'_>, network: &TransportNetwork<'_>, config: &ConjugateConfig, solid: F,
) -> Result<CoupledTransportSolution, TransportError>
where F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    let initial = network.initial_references(cx)?;
    solve_coupled_transport_from(cx, network, config, &initial, solid)
}

/// Continue from an explicit complete reference vector. With fixed relaxation
/// and the same deterministic callback, network and budgets, the numerical
/// tail is reproducible. Aitken's per-branch history restarts; its continuation
/// is not a bitwise replay. No hidden previous solid field is required.
///
/// # Errors
/// As for `solve_coupled_transport`, plus invalid resumed temperatures or arity.
pub fn solve_coupled_transport_from<F>(
    cx: &Cx<'_>, network: &TransportNetwork<'_>, config: &ConjugateConfig,
    initial_references_k: &[f64], solid: F,
) -> Result<CoupledTransportSolution, TransportError>
where F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    solve_driver(cx, network, config, initial_references_k, None, solid)
}

/// Accelerate the complete shared-solid/mixed-air fixed point with bounded
/// vector IQN-ILS. A single history spans all branches, not one scalar per
/// branch or per mixing node. Startup and rank-zero steps use the configured
/// relaxation. A nonpositive extrapolated absolute temperature discards the
/// secants and falls back to that same relaxation; invalid arithmetic refuses.
/// No rejected extrapolation invokes the solid callback.
///
/// # Errors
/// All refusals of [`solve_coupled_transport`], plus invalid acceleration policy
/// and attributed non-finite accelerator arithmetic. A real solid or transport
/// refusal remains terminal, never a request to retry with different physics.
pub fn solve_coupled_transport_iqn<F>(
    cx: &Cx<'_>, network: &TransportNetwork<'_>, config: &ConjugateConfig,
    acceleration: IqnIlsConfig, solid: F,
) -> Result<CoupledTransportSolution, TransportError>
where F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    let initial = network.initial_references(cx)?;
    solve_coupled_transport_iqn_from(cx, network, config, &initial, acceleration, solid)
}

/// Warm-start IQN-ILS from a complete reference vector with EMPTY secant and
/// Aitken histories. This is valid continuation, not bitwise replay of an
/// interrupted accelerated solve. Cancellation always retains the exact
/// references used by the current solid callback, never half an update.
///
/// # Errors
/// As for [`solve_coupled_transport_iqn`], plus malformed initial references.
pub fn solve_coupled_transport_iqn_from<F>(
    cx: &Cx<'_>, network: &TransportNetwork<'_>, config: &ConjugateConfig,
    initial_references_k: &[f64], acceleration: IqnIlsConfig, solid: F,
) -> Result<CoupledTransportSolution, TransportError>
where F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    solve_driver(cx, network, config, initial_references_k, Some(acceleration), solid)
}

fn solve_driver<F>(
    cx: &Cx<'_>, network: &TransportNetwork<'_>, config: &ConjugateConfig,
    initial_references_k: &[f64], acceleration: Option<IqnIlsConfig>, mut solid: F,
) -> Result<CoupledTransportSolution, TransportError>
where F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    checkpoint(cx, 0, initial_references_k)?;
    let seed = network.initial_references(cx)
        .map_err(|error| at_iteration(error, 0, initial_references_k))?;
    let mut active = Vec::new();
    let mut admission_paths = Vec::new();
    for branch in 0..network.hydraulics().branches.len() {
        checkpoint(cx, 0, initial_references_k)?;
        let range = network.row_range(branch);
        if !range.is_empty() {
            admission_paths.push(network.path(branch, seed[range.start])?);
            active.push((branch, range));
        }
    }
    admit(&admission_paths, config, initial_references_k)?;
    for &reference in initial_references_k {
        if reference <= 0.0 {
            return Err(TransportError::InvalidInput("resumed absolute reference temperatures must be positive"));
        }
    }
    let mut accelerator = acceleration.map(|policy| IqnIls::new(initial_references_k.len(), policy))
        .transpose().map_err(iqn_error)?;
    let probe_config = ConjugateConfig {
        max_iterations: 1, relaxation: Relaxation::Fixed { omega: 1.0 }, ..*config
    };
    let mut relaxers: Vec<Option<AitkenRelaxation>> = active.iter().map(|_| match config.relaxation {
        Relaxation::Fixed { .. } => None,
        Relaxation::Aitken { omega_init, omega_max } => Some(AitkenRelaxation::new(omega_init, omega_max)),
    }).collect();
    let mut reference = initial_references_k.to_vec();
    let mut last_change = 0.0_f64;
    for iteration in 0..config.max_iterations {
        checkpoint(cx, iteration, &reference)?;
        let response = solid(cx, &reference)
            .map_err(|error| at_iteration(TransportError::Airflow(error), iteration, &reference))?;
        checkpoint(cx, iteration, &reference)?;
        if response.len() != reference.len() {
            return Err(AirflowError::SolidResponseArity { expected: reference.len(), found: response.len() }.into());
        }
        let walls: Vec<f64> = response.iter().map(|state| state.mean_wall_temperature_k).collect();
        let transport = network.march(cx, &walls)
            .map_err(|error| at_iteration(error, iteration, &reference))?;
        let next = &transport.reference_temperatures_k;
        let mut accepted = true;
        let mut balances = Vec::with_capacity(reference.len());
        let mut omegas = Vec::with_capacity(active.len());
        last_change = 0.0;
        for (ordinal, (branch, range)) in active.iter().enumerate() {
            checkpoint(cx, iteration, &reference)?;
            let inlet = transport.branches[*branch].inlet_temperature_k
                .ok_or(TransportError::Branch { branch: *branch, reason: "heated branch has no transported inlet" })?;
            let path = network.path(*branch, inlet)?;
            let states = &response[range.clone()];
            match solve_conjugate_from(cx, &path, &probe_config, &reference[range.clone()], |_, _| Ok(states.to_vec())) {
                Ok(probe) => balances.extend(probe.balance.regions),
                Err(AirflowError::ConjugateNotConverged { .. }) => accepted = false,
                Err(error) => return Err(at_iteration(TransportError::Airflow(error), iteration, &reference)),
            }
            let total_area = checked_sum("transport branch wetted area", path.segments().iter().map(|segment| segment.area_m2()))?;
            let mut scalar = 0.0;
            for ((&old, &new), segment) in reference[range.clone()].iter().zip(&next[range.clone()]).zip(path.segments()) {
                let delta = checked("transport reference residual", new - old)?;
                last_change = last_change.max(delta.abs());
                scalar = checked("transport weighted reference residual", scalar + delta * (segment.area_m2() / total_area))?;
            }
            let omega = checked("transport relaxation omega", match config.relaxation {
                Relaxation::Fixed { omega } => omega,
                Relaxation::Aitken { omega_init, omega_max: _ } => relaxers[ordinal].as_mut()
                    .map_or(omega_init, |relaxer| relaxer.next_omega(scalar)),
            })?;
            omegas.push(omega);
        }
        checkpoint(cx, iteration, &reference)?;
        if accepted {
            // Publish the references actually used by this solid response,
            // not the next proposed vector. No unaudited final solve is needed.
            return Ok(CoupledTransportSolution {
                reference_temperatures_k: reference, solid: response, transport,
                region_balances: balances, iterations: iteration + 1,
                max_reference_change_k: last_change,
            });
        }
        if let Some(iqn) = accelerator.as_mut() {
            // Zero fallback makes the rank-zero proposal a no-op; branch-local
            // omegas below still own startup. History stores ACTUAL map samples.
            let step = iqn.step(&reference, next, 0.0).map_err(iqn_error)?;
            checkpoint(cx, iteration, &reference)?;
            if step.used_columns > 0 {
                if step.values.iter().all(|&value| value > 0.0) {
                    reference = step.values;
                    continue;
                }
                // Do not clip a vector extrapolation: clipping would mix
                // unrelated states and keep history that already left the
                // absolute-temperature domain. Restart with declared relaxation.
                iqn.reset();
            }
        }
        // Stage the entire relaxation. A cancelled update resumes the old
        // iteration, never a vector containing half old and half new entries.
        let mut relaxed = reference.clone();
        for ((_, range), omega) in active.iter().zip(omegas) {
            for slot in range.clone() {
                checkpoint(cx, iteration, &reference)?;
                let updated = checked("relaxed transport reference", reference[slot] + omega * (next[slot] - reference[slot]))?;
                if updated <= 0.0 {
                    return Err(TransportError::InvalidInput("relaxation produced a nonpositive absolute temperature"));
                }
                relaxed[slot] = updated;
            }
        }
        checkpoint(cx, iteration, &reference)?;
        reference = relaxed;
    }
    Err(AirflowError::ConjugateNotConverged {
        iterations: config.max_iterations, max_change_bits: last_change.to_bits(),
        tolerance_bits: config.temperature_tolerance_k.to_bits(),
    }.into())
}

fn checkpoint(cx: &Cx<'_>, iteration: usize, reference: &[f64]) -> Result<(), TransportError> {
    cx.checkpoint().map_err(|_| cancelled(iteration, reference))
}
fn cancelled(iteration: usize, reference: &[f64]) -> TransportError {
    AirflowError::Cancelled { iteration, references_k: reference.to_vec() }.into()
}
fn at_iteration(error: TransportError, iteration: usize, reference: &[f64]) -> TransportError {
    match error {
        TransportError::Interrupted | TransportError::Airflow(AirflowError::Cancelled { .. }) => cancelled(iteration, reference),
        other => other,
    }
}

#[cfg(test)]
mod tests;

/// Implicit FEM/air-network tangents and adjoints at a checked coupled state.
pub mod sensitivity;
