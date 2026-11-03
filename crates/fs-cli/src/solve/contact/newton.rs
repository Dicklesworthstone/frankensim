//! The native steady path uses the original complete bulk/contact Newton owner.

use fs_conduction::interface::temperature::TemperatureContacts;
use fs_conduction::interface::temperature::steady::newton::{
    ContactSteadyNewtonSolution, SteadyContactNewtonConfig,
};
use fs_conduction::{
    ConductionProblem, ConductionReport, ConductionSolution, LinearSolveEvidence, Nonlinearity,
    ResidualClaim, SolveConfig, StopReason,
};
use fs_exec::Cx;
use fs_solver::{Globalization, LineSearchConfig, NewtonKrylovConfig};

use super::super::{SolveRefusal, radiation::SolidSolution};
use super::{Prepared, Work, bad, json_string, lower, number, poll};

const SCOPE: &str = "Source-backed matching P1 contact and heterogeneous bulk conductivity share the original coupled damped Newton residual and analytic tangent. Continuous linear source curves use the native primal generalized tangent at knots; each physical trial still passes its pinned source query. The accepted field passes independent watt residual and heat-balance gates. One thermal response is one complete coupled solve; residual callbacks, tangent assemblies, outer iterations, inner iterations and accepted damping are retained separately. No contact adjoint, arbitrary-law convergence or uniqueness, continuum bound, moving interface, radiation or new material authority is inferred.";

const MAX_NEWTON_BACKTRACKS: usize = 4096;

/// Curves have a source-owned primal tangent even at their knots. Surface laws
/// keep the existing fixed point: a triangulated source need not have a unique
/// gradient at the caller's operating point.
pub(super) fn supported(prepared: &Prepared, config: &SolveConfig) -> bool {
    matches!(config.nonlinearity, Nonlinearity::Newton { .. })
        && prepared.supports_coupled_newton()
}

pub(super) fn controls(
    config: &SolveConfig,
    tolerance: f64,
    max_pairs: usize,
    vertices: usize,
) -> Result<SteadyContactNewtonConfig, SolveRefusal> {
    let Nonlinearity::Newton { line_search } = config.nonlinearity else {
        return Err(bad("the coupled Newton owner requires Newton controls"));
    };
    if !(config.linear.tolerance.is_finite() && config.linear.tolerance > 0.0) {
        return Err(bad(
            "contact Newton linear tolerance must be finite and positive",
        ));
    }
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err(bad(
            "contact Newton physical tolerance must be finite and positive",
        ));
    }
    if !(line_search.armijo_c.is_finite()
        && line_search.armijo_c > 0.0
        && line_search.armijo_c < 1.0
        && line_search.shrink.is_finite()
        && line_search.shrink > 0.0
        && line_search.shrink < 1.0)
    {
        return Err(bad(
            "contact Newton Armijo and contraction controls must lie in (0, 1)",
        ));
    }
    if line_search.max_backtracks > MAX_NEWTON_BACKTRACKS {
        return Err(bad("contact Newton backtracking allowance exceeds 4096"));
    }
    // Match the numerical owner's repeated multiplication exactly. `powi`
    // can round above the final allowed trial and silently remove it.
    let mut minimum_step = 1.0;
    for _ in 0..line_search.max_backtracks {
        minimum_step *= line_search.shrink;
        if minimum_step == 0.0 {
            return Err(bad("contact Newton minimum backtracking step underflowed"));
        }
    }
    if minimum_step * line_search.shrink >= minimum_step {
        return Err(bad(
            "contact Newton minimum backtracking step cannot shrink further",
        ));
    }
    // The heat defect sums free-node residuals. Leave room below its separate
    // watt gate even when a uniform rung contains many more free vertices.
    let absolute_tolerance = tolerance * 0.1 / fs_math::det::sqrt(vertices.max(1) as f64);
    if absolute_tolerance == 0.0 {
        return Err(bad("contact Newton numerical tolerance underflowed"));
    }
    let max_residual_evaluations = line_search
        .max_backtracks
        .checked_add(2)
        .and_then(|per_iteration| per_iteration.checked_mul(config.stop.max_iterations))
        .and_then(|evaluations| evaluations.checked_add(1))
        .ok_or_else(|| bad("contact Newton residual-work allowance overflowed"))?;
    // Keep the native project's linear tolerance. The original numerical owner
    // requires a forcing value strictly below one.
    let forcing = config.linear.tolerance.min(0.5);
    Ok(SteadyContactNewtonConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance,
            relative_tolerance: 0.0,
            linear_restart: config.linear.restart,
            max_linear_cycles: config.linear.max_iterations,
            forcing_minimum: forcing,
            forcing_maximum: forcing,
            globalization: Globalization::LineSearch(LineSearchConfig {
                armijo: line_search.armijo_c,
                contraction: line_search.shrink,
                minimum_step,
            }),
            ..NewtonKrylovConfig::default()
        },
        max_iterations: config.stop.max_iterations,
        max_residual_evaluations,
        max_face_pairs: max_pairs,
        residual_tolerance_w: tolerance,
        energy_tolerance_w: tolerance,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>,
    initial: &[f64],
    config: &SolveConfig,
    tolerance: f64,
    max_pairs: usize,
    work: &mut Work,
) -> Result<(SolidSolution, String), SolveRefusal> {
    poll(cx)?;
    let accepted = contacts
        .solve_steady_newton(
            cx,
            problem,
            work.last_temperature.as_deref().unwrap_or(initial),
            controls(config, tolerance, max_pairs, initial.len())?,
        )
        .map_err(lower)?;
    finish(cx, problem, accepted, tolerance, work, SCOPE, "")
}

/// Format a counted numerical response after its original physical owner has
/// accepted the unchanged field. A composition may supply a stronger endpoint
/// audit and separately declared work without inventing Newton history.
pub(super) fn finish(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    accepted: ContactSteadyNewtonSolution,
    tolerance: f64,
    work: &mut Work,
    scope: &str,
    composition: &str,
) -> Result<(SolidSolution, String), SolveRefusal> {
    let krylov = accepted
        .newton
        .history
        .iter()
        .try_fold(0usize, |total, row| {
            total
                .checked_add(row.linear_iterations)
                .ok_or_else(|| bad("counted contact Newton Krylov work overflowed"))
        })?;
    let mut history = Vec::with_capacity(accepted.newton.history.len());
    let mut linear = Vec::with_capacity(accepted.newton.history.len());
    for row in &accepted.newton.history {
        poll(cx)?;
        history.push(format!(
            concat!(
                "{{\"iteration\":{},\"residual_before_w\":{},\"residual_after_w\":{},",
                "\"forcing\":{},\"linear_relative_residual\":{},\"linear_iterations\":{},",
                "\"newton_step_norm_k\":{},\"accepted_step_length\":{}}}"
            ),
            row.iteration,
            number(row.residual_before)?,
            number(row.residual_after)?,
            number(row.forcing)?,
            number(row.linear_relative_residual)?,
            row.linear_iterations,
            number(row.newton_step_norm)?,
            number(row.step_length)?
        ));
        // These are the numerical owner's actual FGMRES observations, not a
        // replacement solve against a frozen-contact operator.
        linear.push(LinearSolveEvidence {
            nonlinear_iteration: row.iteration,
            method: "fgmres",
            iterations: row.linear_iterations,
            reported: ResidualClaim::TrueEuclidean(row.linear_relative_residual),
            true_relative_residual: row.linear_relative_residual,
            converged_true: row.linear_relative_residual <= row.forcing,
            stall: None,
        });
    }
    let mut states = Vec::with_capacity(accepted.physical_contacts.states().len());
    for state in accepted.physical_contacts.states() {
        poll(cx)?;
        states.push(format!(concat!(
            "{{\"interface\":{},\"side_a_face\":{},\"side_b_face\":{},",
            "\"side_a_temperature_k\":{},\"side_b_temperature_k\":{},\"query_temperature_k\":{},",
            "\"resistance_m2_k_w\":{},\"selected_claim\":{},\"usage_receipt\":{}}}"),
            json_string(&state.interface), state.pair.side_a, state.pair.side_b,
            number(state.side_a_temperature_k)?, number(state.side_b_temperature_k)?,
            number(state.interface_temperature_k)?, number(state.resistance.value_m2_k_per_w())?,
            json_string(&state.resistance.receipt().selected.0.to_hex()),
            json_string(&state.resistance.receipt().content_hash().to_hex())));
    }
    let interface_fluxes = accepted
        .physical_contacts
        .frozen_interfaces()
        .fluxes(&accepted.temperature)
        .map_err(lower)?;
    // One invocation solves the complete coupled physical system. Residual
    // callbacks and tangent assemblies are separate work, never extra solves.
    work.record(1, accepted.newton.iterations, krylov)?;
    let receipt = format!(
        concat!(
            "{{\"model\":\"source-backed-face-mean-temperature\",\"algorithm\":\"coupled-newton\",",
            "\"thermal_responses\":1,\"nonlinear_iterations\":{},\"krylov_iterations\":{},",
            "\"total_thermal_responses\":{},\"total_nonlinear_iterations\":{},\"total_krylov_iterations\":{},",
            "\"residual_evaluations\":{},\"jacobian_assemblies\":{},",
            "\"residual_norm_w\":{},\"residual_tolerance_w\":{},\"energy_closure_w\":{},",
            "\"energy_tolerance_w\":{},\"newton_history\":[{}],\"face_states\":[{}],\"scope\":{}{}}}"
        ),
        accepted.newton.iterations,
        krylov,
        work.responses,
        work.nonlinear,
        work.krylov,
        accepted.residual_evaluations,
        accepted.jacobian_assemblies,
        number(accepted.residual_norm_w)?,
        number(tolerance)?,
        number(accepted.energy.closure_w)?,
        number(tolerance)?,
        history.join(","),
        states.join(","),
        json_string(scope),
        composition,
    );
    let conduction = ConductionSolution {
        temperature: accepted.temperature,
        report: ConductionReport {
            iterations: accepted.newton.iterations,
            residual_history: accepted
                .newton
                .history
                .iter()
                .map(|row| row.residual_before)
                .collect(),
            final_residual: accepted.residual_norm_w,
            residual_threshold: tolerance,
            stop_reason: StopReason::ResidualTolerance,
            linear,
            energy: accepted.energy,
            material_provenance: problem.element_materials.map_or_else(
                || problem.material.provenance(),
                |assigned| assigned.provenance(),
            ),
            material_receipts: problem.element_materials.map_or_else(
                || problem.material.receipts().len(),
                |assigned| assigned.receipts().len(),
            ),
            element_material_identity: problem
                .element_materials
                .map(|assigned| assigned.identity()),
            interface_fluxes,
            robin_fluxes: accepted.robin_fluxes,
            free_dofs: accepted.free_residual_w.len(),
            elements: problem.mesh.element_count(),
        },
    };
    poll(cx)?;
    work.last_temperature = Some(conduction.temperature.clone());
    Ok((SolidSolution::from_conduction(conduction), receipt))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_conduction::LineSearch;

    #[test]
    fn controls_reject_invalid_tolerances_before_forcing_clamp() {
        for tolerance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -1.0] {
            let mut config = SolveConfig::default();
            config.linear.tolerance = tolerance;
            assert!(controls(&config, 1.0e-8, 1, 8).is_err());
            assert!(controls(&SolveConfig::default(), tolerance, 1, 8).is_err());
        }
    }

    #[test]
    fn controls_preserve_the_last_admitted_backtrack() {
        for shrink in [0.1, 0.3, 0.7, 0.9, 0.99] {
            for max_backtracks in [0, 1, 24, 100] {
                let config = SolveConfig {
                    nonlinearity: Nonlinearity::Newton {
                        line_search: LineSearch {
                            shrink,
                            max_backtracks,
                            ..LineSearch::default()
                        },
                    },
                    ..SolveConfig::default()
                };
                let controls = controls(&config, 1.0e-8, 1, 8).unwrap();
                let Globalization::LineSearch(line) = controls.newton.globalization else {
                    panic!("contact Newton requires line search");
                };
                let mut step = 1.0;
                let mut trials = 0;
                while step >= line.minimum_step {
                    trials += 1;
                    step *= line.contraction;
                }
                assert_eq!(trials, max_backtracks + 1);
            }
        }
    }

    #[test]
    fn controls_refuse_unrepresentable_backtracking() {
        for max_backtracks in [1075, MAX_NEWTON_BACKTRACKS + 1, usize::MAX] {
            let config = SolveConfig {
                nonlinearity: Nonlinearity::Newton {
                    line_search: LineSearch {
                        max_backtracks,
                        ..LineSearch::default()
                    },
                },
                ..SolveConfig::default()
            };
            assert!(controls(&config, 1.0e-8, 1, 8).is_err());
        }
    }
}
