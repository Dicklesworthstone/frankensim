//! Periodic physical driving and nonlinear port feedback in the same full ODE.

use super::{
    AutonomousPortFeedback, PeriodicPortInput, PhsFeedbackStateJacobian, PhsOrbitAudit,
    PhsOrbitError, PhsStorageHessianAction, audit_orbit_with_input_effort, bad,
    variational_feedback::{FeedbackVariationalProblem, feedback_tangent_problem},
};
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, ScaledHbOrbit,
    ShootBudget, VariationalShootControls, VariationalShootError, VariationalShootOrbit,
    solve_hb_adaptive_scaled_seeded, solve_hb_scaled_seeded, solve_shooting_variational,
};
use fs_phs::PortHamiltonian;

fn problem<'a>(
    system: &'a PortHamiltonian,
    input: &'a PeriodicPortInput,
    feedback: &'a dyn AutonomousPortFeedback,
    feedback_jacobian: &'a PhsFeedbackStateJacobian<'a>,
    budget: &HbBudget,
    hessian: Option<&'a PhsStorageHessianAction<'a>>,
) -> Result<FeedbackVariationalProblem<'a>, PhsOrbitError> {
    // Resolve the forcing bandwidth before preparing or evaluating any source.
    // Zero padding does not consume harmonic resolution; tiny nonzero lines do.
    if input.highest_nonzero_harmonic() > budget.harmonics {
        return Err(bad(
            "nonzero forcing harmonic exceeds the retained HB truncation",
        ));
    }
    if input.port_count() != system.port_dim() {
        return Err(bad(
            "periodic input must cover every original physical port",
        ));
    }
    let mut problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)?;
    if !problem.original.tangent.available() {
        return Err(bad(
            "complete original storage tangent required for analytic feedback HB",
        ));
    }
    problem.external = Some(input);
    Ok(problem)
}

/// Resolve a periodically driven nonlinear feedback loop in original pHS states.
///
/// Every residual uses `(J-R)*grad H + G*(u_periodic(t) + u_feedback(x,y))`,
/// where `y=G^T*grad H`. The original storage, native resistance, physical port
/// order and all material-memory coordinates remain present. The drive fixes
/// the frequency and phase; no autonomous phase constraint or unknown frequency
/// is inserted. The feedback law remains instantaneous and autonomous, although
/// the combined system is externally driven.
///
/// The original complete storage tangent is discovered automatically. The
/// authored TOTAL [`PhsFeedbackStateJacobian`] includes its flow/Hessian chain;
/// the state-independent drive contributes no additional state derivative.
/// This shared analytic AFT Jacobian is evaluated once per physical sample,
/// rather than perturbing every Fourier coefficient.
///
/// Pack `seed` using [`HbOrbit::packed_seed`] with a forced anchor. The initial
/// order must retain every nonzero drive harmonic. [`HbScales`] has the same
/// original-coordinate semantics as the ordinary forced pHS solver. This
/// returns a fixed-order sampled balance, not omitted-tail, continuous-orbit,
/// stability, passivity or excitation-law certification. Source callbacks must
/// be deterministic and read-only; their memory belongs in explicit states.
///
/// # Errors
/// Invalid drive/feedback ports, forcing bandwidth, seed or scales; unavailable
/// or refusing full tangents; original source refusal; or native HB failure.
/// Declared tangents never silently fall back to finite differences.
pub fn solve_forced_feedback_harmonic_balance_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    solve(
        system,
        input,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        scales,
        None,
    )
}

/// The same original forced feedback loop with a caller-authored COMPLETE
/// storage Hessian action. All memory and nonlinear cross derivatives must
/// follow [`PhsStorageHessianAction`]; the feedback Jacobian already includes
/// the flow chain through that same original storage.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_harmonic_balance_analytic_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
    hessian: &PhsStorageHessianAction<'_>,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    solve(
        system,
        input,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        scales,
        Some(hessian),
    )
}

#[allow(clippy::too_many_arguments)]
fn solve(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    let problem = problem(system, input, feedback, feedback_jacobian, budget, hessian)?;
    let result = solve_hb_scaled_seeded(
        &problem,
        HbAnchor::Forced {
            omega: input.omega_rad_s(),
        },
        seed,
        budget,
        scales,
    );
    if problem.original.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// Refine the same periodically driven nonlinear feedback response, retaining
/// the original external phase and every physical state at each harmonic order.
/// The shared owner checks sampling, retained and omitted equation harmonics,
/// and successive waveforms in the caller's physical tolerances. It does not
/// infer an infinite tail, global stability, or a physical constitutive model.
/// The initial seed must already resolve every nonzero external-drive harmonic.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_harmonic_balance_adaptive_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    solve_adaptive(
        system,
        input,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        refinement,
        scales,
        None,
    )
}

/// Adaptive forced feedback with an explicitly authored complete storage
/// Hessian. Refusal of either declared tangent stops the shared analytic solve;
/// source work, original port order and finite-band semantics are unchanged.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_harmonic_balance_analytic_adaptive_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
    hessian: &PhsStorageHessianAction<'_>,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    solve_adaptive(
        system,
        input,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        refinement,
        scales,
        Some(hessian),
    )
}

#[allow(clippy::too_many_arguments)]
fn solve_adaptive(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    let problem = problem(system, input, feedback, feedback_jacobian, budget, hessian)?;
    let result = solve_hb_adaptive_scaled_seeded(
        &problem,
        HbAnchor::Forced {
            omega: input.omega_rad_s(),
        },
        seed,
        budget,
        refinement,
        scales,
    );
    if problem.original.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// Shoot the same periodically driven feedback ODE with its full analytic tangent.
///
/// The original pHS and instantaneous feedback law supply every RK4 stage. The
/// existing total feedback tangent differentiates those SAME stages, including
/// all nonlinear storage and material-memory cross terms. Native complete
/// storage metadata/actions are discovered automatically. A merely quadratic
/// part of nonlinear storage returns TangentUnavailable; use the analytic entry
/// for an authored complete Hessian. Neither declared tangent is replaced by
/// finite differences. The feedback Jacobian must be total, including its flow
/// chain, as specified by [`PhsFeedbackStateJacobian`].
///
/// Period is exactly `2*pi/input.omega_rad_s()`, and every authored drive phase
/// is retained. The initial grid must resolve every nonzero external harmonic.
/// Physical seed and state scales keep all original coordinates. The returned
/// map and Floquet spectrum concern one discrete time grid: temporal convergence,
/// primitive period, global attraction and excitation physics need independent
/// evidence. The separate audit still reports drive, feedback and native R work.
/// Original sources are read-only; refusal or cancellation publishes no cycle.
///
/// Before preparation this adapter reserves `16*n²+8*n*m+64*n+48*m+256` bytes
/// and `32*n²+(4*maps*steps+1)*(4*n³+32*n²+16*n²*m+32*n*m+32*m*(h+2))`
/// work units from the supplied controls, with `h` the authored harmonic slots.
/// These cover original tangent metadata, scratch, contractions and periodic
/// sampling; arbitrary callback internals remain caller-owned. The remaining
/// caps admit the existing shared variational shooting owner. Cancellation is
/// checked before preparation and at the native stage/product boundaries.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_shooting_variational(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: &[f64],
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    shoot(
        system,
        input,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        controls,
        None,
        cancelled,
    )
}

/// The same forced feedback cycle with a complete original storage Hessian
/// action, following [`PhsStorageHessianAction`]. The authored feedback Jacobian
/// already includes this same flow/Hessian chain exactly once. Source, tangent,
/// bandwidth, resource and cancellation refusals have the automatic entry's
/// semantics, with no surrogate storage or fallback derivative.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_shooting_variational_analytic(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: &[f64],
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    shoot(
        system,
        input,
        feedback,
        feedback_jacobian,
        seed,
        budget,
        controls,
        Some(hessian),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn shoot(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    seed: &[f64],
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    if cancelled() {
        return Err(PhsOrbitError::Variational(VariationalShootError::Cancelled));
    }
    let (n, m) = (system.state_dim(), system.port_dim());
    if n == 0
        || n > 256
        || n > controls.maximum_dimension
        || input.port_count() != m
        || feedback.port_count() != m
        || seed.len() != n
        || seed.iter().any(|v| !v.is_finite())
        || controls.state_scales.len() != n
        || controls
            .state_scales
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        || !controls.period_scale.is_finite()
        || controls.period_scale <= 0.0
        || !controls.minimum_autonomous_excursion.is_finite()
        || controls.minimum_autonomous_excursion <= 0.0
        || !budget.tolerance.is_finite()
        || budget.tolerance <= 0.0
        || budget.tolerance >= 1.0
        || budget.steps_per_period == 0
        || !(1..=4096).contains(&budget.max_newton)
        || controls.maximum_period_maps == 0
        || input.highest_nonzero_harmonic() > (budget.steps_per_period - 1) / 2
    {
        return Err(bad(
            "complete finite forced feedback/state data, positive physical scales and a bounded shooting grid resolving every authored input harmonic required",
        ));
    }
    let period = core::f64::consts::TAU / input.omega_rad_s();
    if !period.is_finite() || period <= 0.0 {
        return Err(bad(
            "positive representable physical forcing period required",
        ));
    }
    let (n, m) = (n as u128, m as u128);
    let bytes = 16 * n * n + 8 * n * m + 64 * n + 48 * m + 256;
    let work = (controls.maximum_period_maps as u128)
        .checked_mul(budget.steps_per_period as u128)
        .and_then(|v| v.checked_mul(4))
        .and_then(|v| v.checked_add(1))
        .and_then(|v| {
            v.checked_mul(
                4 * n * n * n
                    + 32 * n * n
                    + 16 * n * n * m
                    + 32 * n * m
                    + 32 * m * (input.harmonic_count() as u128 + 2),
            )
        })
        .and_then(|v| v.checked_add(32 * n * n))
        .ok_or_else(|| bad("forced feedback variational work envelope overflows"))?;
    for (quantity, requested, cap) in [
        (
            "forced feedback variational payload bytes",
            bytes,
            controls.maximum_payload_bytes as u128,
        ),
        (
            "forced feedback variational work",
            work,
            controls.maximum_work as u128,
        ),
    ] {
        if requested > cap {
            return Err(PhsOrbitError::Variational(VariationalShootError::Budget {
                quantity,
                requested,
                cap,
            }));
        }
    }
    let mut admitted = controls.clone();
    admitted.maximum_payload_bytes -= bytes as usize;
    admitted.maximum_work -= work as u64;
    let mut problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)?;
    problem.external = Some(input);
    let result =
        solve_shooting_variational(&problem, seed, period, budget, &admitted, &mut cancelled);
    match result {
        Err(VariationalShootError::Cancelled) => {
            Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
        }
        _ if problem.original.source_refused() => Err(PhsOrbitError::SourceRefused),
        other => other.map_err(PhsOrbitError::Variational),
    }
}

/// Separate physical drive and feedback work on one sampled candidate cycle.
#[derive(Debug, Clone)]
pub struct PhsForcedFeedbackAudit {
    /// Original storage, TOTAL supplied work, native R loss and their balance.
    pub cycle: PhsOrbitAudit,
    /// Cycle mean of `grad(H)^T G u_periodic`, in watts; signed external supply.
    pub mean_external_drive_power_w: f64,
    /// Cycle mean of `grad(H)^T G u_feedback`, in watts. Negative means net
    /// absorption by the authored feedback law, separate from native R loss.
    /// This signed work does not by itself identify material heat or radiation.
    pub mean_feedback_power_w: f64,
}

/// Audit the actual original periodic drive, nonlinear feedback and pHS storage.
///
/// Each sample evaluates the original gradient once. Both signed component
/// powers use the original compensated physical map, not rounded port outputs.
/// The cycle residual uses their combined physical input in the same original
/// pHS balance; separately rounded large opposing component means need not sum
/// bit-for-bit to that combined account. Passive feedback absorption is never
/// silently folded into native resistance or labeled heat.
///
/// This is finite-grid quadrature. The requested 32..=65536 samples must resolve
/// the authored drive and retained state orders; nonlinear work can need finer
/// sampling. Compare successive grids before treating its balance as resolved.
/// No storage, law, input or candidate state is changed by success or refusal.
pub fn audit_forced_feedback_orbit(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    orbit: &HbOrbit,
    samples: usize,
) -> Result<PhsForcedFeedbackAudit, PhsOrbitError> {
    let state_order = orbit
        .coeffs
        .first()
        .map_or(0, |row| row.len().saturating_sub(1));
    if input.port_count() != system.port_dim()
        || feedback.port_count() != system.port_dim()
        || orbit.omega.to_bits() != input.omega_rad_s().to_bits()
        || !(32..=65_536).contains(&samples)
        || input.highest_nonzero_harmonic().max(state_order) > (samples - 1) / 2
    {
        return Err(bad(
            "forced feedback audit needs matching ports/frequency and a grid resolving drive and state harmonics",
        ));
    }
    let mut drive = vec![0.0; system.port_dim()];
    let mut sums = [0.0; 2];
    let mut corrections = [0.0; 2];
    let cycle = audit_orbit_with_input_effort(
        system,
        orbit,
        samples,
        |time, state, effort, flows, out| {
            if !feedback.input(state, flows, out) || !input.sample_into(time, &mut drive) {
                return false;
            }
            let Ok(external_work) = system.supplied_power_from_effort(effort, &drive) else {
                return false;
            };
            let Ok(feedback_work) = system.supplied_power_from_effort(effort, out) else {
                return false;
            };
            for (i, work) in [external_work, feedback_work].into_iter().enumerate() {
                let adjusted = work - corrections[i];
                let next = sums[i] + adjusted;
                corrections[i] = (next - sums[i]) - adjusted;
                sums[i] = next;
            }
            for (value, imposed) in out.iter_mut().zip(&drive) {
                *value += imposed;
            }
            out.iter().all(|v| v.is_finite())
        },
    )?;
    let [mean_external_drive_power_w, mean_feedback_power_w] = sums.map(|v| v / samples as f64);
    if !mean_external_drive_power_w.is_finite() || !mean_feedback_power_w.is_finite() {
        return Err(bad("forced feedback cycle work reduction overflow"));
    }
    Ok(PhsForcedFeedbackAudit {
        cycle,
        mean_external_drive_power_w,
        mean_feedback_power_w,
    })
}
