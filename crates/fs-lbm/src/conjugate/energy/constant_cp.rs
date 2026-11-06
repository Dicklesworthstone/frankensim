//! Converged limited convection for the native constant-Cp energy equation.
//!
//! The existing upwind matrix carries a deferred minmod face correction.
//! Refreshing that correction is part of the solve, including at the accepted
//! endpoint; a single frozen correction is not a converged transport solution.
//! The same shared face sampler is used for source enthalpy. See Thakur & Shyy,
//! Numerical Heat Transfer B 24 (1993), for deferred correction and consistent
//! interior face treatment: <https://ntrs.nasa.gov/citations/19930071174>.
//!
//! This changes spatial convection only. Boundary stencils remain donor-cell;
//! no general multidimensional maximum principle or temporal-order claim follows.

use super::{
    ChtError, ConvectionScheme, EnergyConfig, EnergySolution, FlowField, FluidProperties,
    PseudoStep, SolidMaterial, ThermalSetup, VoxelDomain, solve_energy_field_once,
};
use fs_exec::CancelGate;

const MAX_RECONSTRUCTION_ITERATIONS: usize = 200;

/// Native temperature transport needs no inferred material chart: constant Cp
/// factors out of minmod. Molecular coefficients, old-time storage, sources,
/// compact components and contact laws retain their original owners.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>,
    nonlinear: Option<(&[f64], Option<&[f64]>)>,
    molecular_conductivity: Option<&[[f64; 3]]>,
    evaluate_only: bool,
    gate: &CancelGate,
) -> Result<EnergySolution, ChtError> {
    super::super::poll(gate)?;
    if evaluate_only {
        // An outer conductivity/material owner asks about this unchanged
        // candidate. Reconstruct from it and evaluate; never solve it away.
        return solve_energy_field_once(
            domain, fluid, solids, flow, setup, config, step, nonlinear,
            molecular_conductivity, None, true, gate,
        );
    }

    // A native upwind solve admits all original inputs and gives a robust
    // starting field without inventing a temperature or altering the history.
    let seed_config = EnergyConfig { scheme: ConvectionScheme::Upwind, ..*config };
    let mut solution = solve_energy_field_once(
        domain, fluid, solids, flow, setup, &seed_config, step, nonlinear,
        molecular_conductivity, None, false, gate,
    )?;
    let mut iterations = solution.report.iterations;
    let unit_scales = vec![1.0; domain.cell_count()];
    let scales = nonlinear.map_or(unit_scales.as_slice(), |(scales, _)| scales);
    let mut linear_config = *config;
    let mut endpoint_residual = f64::INFINITY;

    for _ in 0..MAX_RECONSTRUCTION_ITERATIONS {
        super::super::poll(gate)?;
        let next = solve_energy_field_once(
            domain, fluid, solids, flow, setup, &linear_config, step,
            Some((scales, Some(&solution.temperature))), molecular_conductivity,
            None, false, gate,
        )?;
        iterations = iterations.saturating_add(next.report.iterations);
        let (lo, hi) = next.temperature.iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(lo, hi), &t| (lo.min(t), hi.max(t)),
        );
        let change = next.temperature.iter().zip(&solution.temperature)
            .fold(0.0_f64, |maximum, (a, b)| maximum.max((a - b).abs()));
        solution = next;
        // At a nearly uniform Kelvin field the span can be smaller than an
        // ULP. This allowance only proposes an endpoint check; it does not
        // relax either of the original acceptance residuals below.
        let proposal_roundoff = 32.0 * f64::EPSILON * lo.abs().max(hi.abs());
        if change > 1e-10 * (hi - lo).max(1e-12) + proposal_roundoff {
            continue;
        }
        // This rebuilds the same conservative face law at the UNCHANGED new
        // temperature. Both algebraic and local physical heat residuals must
        // satisfy the caller's tolerance; small iterate changes alone do not.
        match solve_energy_field_once(
            domain, fluid, solids, flow, setup, config, step,
            Some((scales, Some(&solution.temperature))), molecular_conductivity,
            None, true, gate,
        ) {
            Ok(mut accepted) => {
                accepted.report.iterations = iterations;
                super::super::poll(gate)?;
                return Ok(accepted);
            }
            Err(ChtError::SolverNotConverged {
                system: "energy conductivity endpoint" | "energy conductivity heat",
                relative_residual,
                ..
            }) => {
                endpoint_residual = relative_residual;
                // Tighten only the inner solve; the physical acceptance
                // tolerance and the original accounting window stay fixed.
                let tighter = 0.1 * linear_config.tolerance;
                if tighter > 0.0 { linear_config.tolerance = tighter; }
            }
            Err(error) => return Err(error),
        }
    }
    Err(ChtError::SolverNotConverged {
        system: "energy limited convection",
        iterations: MAX_RECONSTRUCTION_ITERATIONS,
        relative_residual: endpoint_residual,
        tolerance: config.tolerance,
    })
}
