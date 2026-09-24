//! Goal comparison for one or more independent air paths on a shared solid.
//!
//! At fixed flow, h and geometry, the exponential air-reference law is affine
//! in wall means: r=R(WT). The physical solid residual is R_c=R_s-B(r-r_bound)
//! and its tangent is J_c=J_s-B*A*W. We solve the transpose interface equation
//! mu=B^T J_s^-T (goal+W^T A^T mu), using the actual material/contact tangent
//! and an analytic reverse air sweep. IQN history belongs to this linear
//! equation; no primal iteration trace or finite difference is differentiated.

use std::{collections::BTreeSet, fmt};

use fs_conduction::adjoint::{DiscreteGoalComparison, RobinGoalFeedback, RobinGoalLinearization};
use fs_conduction::{ConductionError, ConductionProblem, LinearConfig, ThermalInterfaces};
use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig, IqnIlsError};
use fs_exec::Cx;

use super::{
    AirPath, ConjugateConfig, Relaxation, SolidRegionState, admitted_exchange_terms,
    solve_conjugate_from,
};
use crate::AirflowError;
pub use crate::graph::thermal::coupled_transport::sensitivity::InterfaceSolveConfig;

/// Numerical budgets for a residual-checked conjugate goal comparison.
#[derive(Debug, Clone, Copy)]
pub struct CoupledGoalConfig {
    /// Full coupled primal/dual residual tolerance and per-pullback Krylov cap.
    /// The solid seam uses one tenth of this tolerance for inner solves.
    pub solid: LinearConfig,
    /// Existing temperature and per-branch heat-balance admission gates.
    pub primal: ConjugateConfig,
    /// Maximum transpose sweeps and unrelaxed interface residual tolerances.
    /// The goal is first normalized by its largest absolute nodal weight,
    /// making the absolute threshold independent of objective units.
    pub interface: InterfaceSolveConfig,
    /// Fresh bounded IQN history for this transpose equation only.
    pub acceleration: IqnIlsConfig,
}

/// Compared discrete goal and actual work/residual of the interface adjoint.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledGoalComparison {
    /// Full physical coupled residual contributions and explicit remainder.
    /// `dual_iterations` counts all inner Krylov work in the interface solve.
    pub goal: DiscreteGoalComparison,
    /// Completed interface sweeps, including the successful final sweep.
    pub interface_iterations: usize,
    /// Maximum unrelaxed residual for the normalized goal's multiplier.
    pub interface_residual: f64,
    /// Threshold applied to that normalized interface equation.
    pub interface_tolerance: f64,
}

/// A producer, binding, numerical-budget or interruption refusal.
#[derive(Debug)]
pub enum CoupledGoalError {
    /// Invalid paths, port binding, vectors or interface configuration.
    InvalidInput(&'static str),
    /// Solid primal, tangent, residual gate or field refusal.
    Solid(ConductionError),
    /// Air-path or per-branch physical balance refusal.
    Air(AirflowError),
    /// Invalid IQN policy or nonfinite acceleration arithmetic.
    Acceleration(IqnIlsError),
    /// Context refused more work; no partial comparison is returned.
    Interrupted,
    /// Unrelaxed transpose equation failed within its declared sweep budget.
    DidNotConverge {
        iterations: usize,
        residual: f64,
        tolerance: f64,
    },
}
impl fmt::Display for CoupledGoalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "coupled goal: {self:?}")
    }
}
impl std::error::Error for CoupledGoalError {}
impl From<ConductionError> for CoupledGoalError {
    fn from(error: ConductionError) -> Self {
        Self::Solid(error)
    }
}
impl From<AirflowError> for CoupledGoalError {
    fn from(error: AirflowError) -> Self {
        Self::Air(error)
    }
}
type Result<T> = std::result::Result<T, CoupledGoalError>;

/// Compare two full nodal fields using the same coupled solid/air equations.
/// The approximation may be a prolonged coarse field. References for BOTH
/// fields are recomputed by marching the actual paths; none are frozen from
/// the coarse solution. Fixed prescribed values remain unchanged. Ports must
/// share exact names/order/h and matching area with the reference mesh.
///
/// The reference solid is independently residual-checked, every branch passes
/// its own existing temperature/watt gate, and the final fully coupled primal
/// and transpose residuals pass `config.solid.tolerance`. Signed contributions
/// prioritize refinement; the measured linearization remainder must remain in
/// any estimate. This is not a continuum error certificate, a derivative of
/// hydraulics, a moving-interface derivative, or a thermal flow-field model.
///
/// # Errors
/// Refuses missing/duplicate/mismatched ports, invalid fields/budgets, failed
/// material or contact admission, a nonconverged reference, a nonsmooth
/// material tangent, exhausted derivative work, or cancellation. Neither
/// inputs nor a previously retained result are mutated.
#[allow(clippy::too_many_arguments)]
pub fn compare_discrete_goal(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath],
    config: CoupledGoalConfig,
    reference_temperature: &[f64],
    approximate_temperature: &[f64],
    full_nodal_weights: &[f64],
) -> Result<CoupledGoalComparison> {
    poll(cx)?;
    validate(config.interface)?;
    if paths.is_empty() {
        return Err(bad("at least one physical air path is required"));
    }
    let mut names = Vec::new();
    let mut seen = BTreeSet::new();
    for path in paths {
        for name in path.regions() {
            poll(cx)?;
            if !seen.insert(name) {
                return Err(bad("an air region may belong to only one path"));
            }
            names.push(name);
        }
    }
    let linear = RobinGoalLinearization::new(
        cx,
        problem,
        interfaces,
        config.solid,
        reference_temperature,
        approximate_temperature,
        &names,
    )?;
    let response = linear.response();
    let mut accelerator =
        IqnIls::new(names.len(), config.acceleration).map_err(CoupledGoalError::Acceleration)?;
    let approximate_walls = response.wall_means(cx, approximate_temperature)?;
    let mut reference_shift_k = Vec::with_capacity(names.len());
    let mut reference_difference_k = Vec::with_capacity(names.len());
    let mut offset = 0;
    for path in paths {
        poll(cx)?;
        let end = offset + path.segments().len();
        let mut states = Vec::with_capacity(end - offset);
        let mut references = Vec::with_capacity(end - offset);
        for (index, segment) in (offset..end).zip(path.segments()) {
            poll(cx)?;
            let port = &response.ports()[index];
            if port.htc_w_m2_k != segment.htc_w_per_m2_k()
                || (port.area_m2 - segment.area_m2()).abs()
                    > 128.0 * f64::EPSILON * port.area_m2.max(segment.area_m2())
            {
                return Err(bad(
                    "solid and air must share the same coefficient and wetted area",
                ));
            }
            let flux = response
                .robin_fluxes()
                .iter()
                .find(|flux| flux.region == port.name)
                .ok_or_else(|| bad("the reference solid is missing an air-region flux"))?;
            states.push(SolidRegionState::from_robin_flux(flux));
            references.push(port.reference_k);
        }
        let gate = ConjugateConfig {
            max_iterations: 1,
            relaxation: Relaxation::Fixed { omega: 1.0 },
            ..config.primal
        };
        let checked =
            solve_conjugate_from(cx, path, &gate, &references, |_, _| Ok(states.clone()))?;
        let approximate = path.march(&approximate_walls[offset..end])?;
        for (i, (reference, approximate)) in checked
            .march
            .segments
            .iter()
            .zip(&approximate.segments)
            .enumerate()
        {
            poll(cx)?;
            reference_shift_k.push(finite(reference.reference_temperature_k - references[i])?);
            reference_difference_k.push(finite(
                reference.reference_temperature_k - approximate.reference_temperature_k,
            )?);
        }
        offset = end;
    }
    if full_nodal_weights.len() != reference_temperature.len() {
        return Err(bad("one goal weight per full solid node is required"));
    }
    let mut goal_scale = 0.0_f64;
    for &weight in full_nodal_weights {
        poll(cx)?;
        goal_scale = goal_scale.max(finite(weight)?.abs());
    }
    if goal_scale == 0.0 {
        goal_scale = 1.0;
    }
    let weights: Vec<f64> = full_nodal_weights.iter().map(|w| w / goal_scale).collect();
    let mut current = vec![0.0; names.len()];
    let zero = vec![0.0; names.len()];
    let mut total_krylov = 0_usize;
    for iteration in 1..=config.interface.max_iterations {
        poll(cx)?;
        let wall_adjoint = reference_pullback(cx, paths, &current)?;
        let solid = response.pullback(cx, &weights, &wall_adjoint, &zero)?;
        total_krylov = total_krylov
            .checked_add(solid.iterations)
            .ok_or_else(|| bad("coupled derivative work counter overflow"))?;
        let (residual, tolerance) = equation(&current, &solid.references, config.interface)?;
        if residual <= tolerance {
            // Recompute the feedback from B^T lambda actually returned by the
            // solid, not from the interface iterate used on its right side.
            let wall_adjoint = reference_pullback(cx, paths, &solid.references)?
                .into_iter()
                .map(|w| finite(w * goal_scale))
                .collect::<Result<Vec<_>>>()?;
            let adjoint = solid
                .nodal_load
                .iter()
                .map(|w| finite(w * goal_scale))
                .collect::<Result<Vec<_>>>()?;
            let feedback = RobinGoalFeedback {
                reference_shift_k,
                reference_difference_k,
                wall_adjoint,
            };
            let goal = linear.compare(cx, full_nodal_weights, &adjoint, &feedback, total_krylov)?;
            poll(cx)?;
            return Ok(CoupledGoalComparison {
                goal,
                interface_iterations: iteration,
                interface_residual: residual,
                interface_tolerance: tolerance,
            });
        }
        if iteration == config.interface.max_iterations {
            return Err(CoupledGoalError::DidNotConverge {
                iterations: iteration,
                residual,
                tolerance,
            });
        }
        let next = accelerator
            .step(&current, &solid.references, config.interface.relaxation)
            .map_err(CoupledGoalError::Acceleration)?;
        poll(cx)?;
        current = next.values;
    }
    unreachable!("positive bounded interface solve returns")
}

/// A^T times reference weights. The reverse outlet multiplier propagates each
/// downstream reference's dependence on every upstream wall. Fixed inlets
/// terminate that propagation at each independent branch's entrance.
fn reference_pullback(cx: &Cx<'_>, paths: &[AirPath], references: &[f64]) -> Result<Vec<f64>> {
    let mut walls = vec![0.0; references.len()];
    let mut offset = 0;
    for path in paths {
        let mut outlet = 0.0;
        for (i, segment) in path.segments().iter().enumerate().rev() {
            poll(cx)?;
            let (ntu, eps, g) = admitted_exchange_terms(segment, path.capacity_rate_w_per_k())?;
            // Stable d(reference)/d(wall) at arbitrarily small admitted NTU.
            let one_minus_g = if ntu < 1e-3 {
                ntu * (0.5
                    + ntu
                        * (-1.0 / 6.0
                            + ntu
                                * (1.0 / 24.0
                                    + ntu * (-1.0 / 120.0 + ntu * (1.0 / 720.0 - ntu / 5040.0)))))
            } else {
                1.0 - g
            };
            let weight = references[offset + i];
            walls[offset + i] = finite(one_minus_g * weight + eps * outlet)?;
            outlet = finite(g * weight + fs_math::det::exp(-ntu) * outlet)?;
        }
        offset += path.segments().len();
    }
    Ok(walls)
}
fn equation(current: &[f64], next: &[f64], config: InterfaceSolveConfig) -> Result<(f64, f64)> {
    let mut residual = 0.0_f64;
    let mut scale = 0.0_f64;
    for (&a, &b) in current.iter().zip(next) {
        residual = residual.max(finite(b - a)?.abs());
        scale = scale.max(finite(a)?.abs()).max(finite(b)?.abs());
    }
    Ok((
        residual,
        finite(config.absolute_tolerance + config.relative_tolerance * scale)?,
    ))
}
fn validate(config: InterfaceSolveConfig) -> Result<()> {
    if config.max_iterations == 0
        || !config.absolute_tolerance.is_finite()
        || config.absolute_tolerance <= 0.0
        || !(0.0..1.0).contains(&config.relative_tolerance)
        || !config.relaxation.is_finite()
        || config.relaxation <= 0.0
        || config.relaxation > 1.0
    {
        return Err(bad("invalid coupled goal interface budget"));
    }
    Ok(())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("nonfinite coupled goal arithmetic"))
    }
}
fn bad(reason: &'static str) -> CoupledGoalError {
    CoupledGoalError::InvalidInput(reason)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| CoupledGoalError::Interrupted)
}
