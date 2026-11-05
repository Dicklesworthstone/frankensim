//! Resolved forced-branch stability from the original exact state tangent.

use super::{PhsFloquetControls, PhsNumericalStability, admit_candidate, check_orbit, classify};
use crate::phs_orbit::{
    PeriodicPortInput, PhsOrbitError, PhsStorageHessianAction, solve_forced_shooting_variational,
    solve_forced_shooting_variational_analytic,
};
use fs_orbit::{
    HbOrbit, ShootBudget, VariationalShootControls, VariationalShootError, VariationalShootOrbit,
};
use fs_phs::PortHamiltonian;

/// One time resolution of the complete original forced response and tangent.
#[derive(Debug, Clone)]
pub struct PhsVariationalFloquetLevel {
    /// Fixed RK4 steps on the authored forcing period.
    pub steps_per_period: usize,
    /// Original physical state, complete scaled monodromy, all multipliers,
    /// normalized closure residual and actual source-evaluation counts.
    pub orbit: VariationalShootOrbit,
    /// Largest multiplier magnitude, without deleting any state direction.
    pub spectral_radius: f64,
}

/// Three-grid numerical stability assessment of the supplied forced HB branch.
#[derive(Debug, Clone)]
pub struct PhsVariationalForcedStability {
    /// Uses the same state/radius convergence and unit-circle guard as
    /// [`super::analyze_forced_orbit_stability`].
    pub classification: PhsNumericalStability,
    /// Exactly three levels, at N, 2N and 4N steps per forcing period.
    pub levels: Vec<PhsVariationalFloquetLevel>,
    /// ln(finest spectral radius)/period [1/s], when representable.
    pub dominant_growth_rate_per_s: Option<f64>,
}

/// Resolve forced resonance-branch stability using the complete declared
/// original storage tangent, including every explicit material-memory state.
/// Exact quadratic metadata or [`fs_phs::Storage::continuous_tangent`] is
/// discovered automatically. Other nonlinear storage can use the complete
/// Hessian action accepted by [`analyze_forced_orbit_stability_variational_analytic`].
///
/// The shared variational shooting owner differentiates its actual RK4 stages;
/// no trajectory perturbation or finite-difference amplitude is used. Each grid
/// starts from the SAME full physical HB seed and authored forcing phase. The
/// original physical seed-agreement limits, three-grid state/radius comparison
/// and unit-circle guard remain authoritative. No neutral multiplier is removed
/// from this externally forced system, even at zero input amplitude.
///
/// `variational` supplies physical coordinate scales and numerical limits for
/// EACH native time-grid solve, including its pHS adapter. The query runs exactly
/// three such solves and retains their results; the per-solve payload cap excludes
/// already returned levels. Retained numerical arrays are three monodromies,
/// three states/scales/spectra and three bounded Newton traces. Callback internals
/// remain caller-owned, as in the underlying shooting API. Cancellation is polled
/// before admission, between levels and at the native stage/product boundaries;
/// a refusal or cancellation returns no partial stability assessment.
///
/// The tangent is exact for the discrete map, not for the continuous orbit.
/// The measured temporal ladder does not certify harmonic truncation, Floquet
/// conditioning, global stability or hybrid contact events. Source and storage
/// callbacks must be deterministic and read-only.
pub fn analyze_forced_orbit_stability_variational(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    candidate: &HbOrbit,
    controls: &PhsFloquetControls,
    variational: &VariationalShootControls,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsVariationalForcedStability, PhsOrbitError> {
    analyze(
        system,
        input,
        candidate,
        controls,
        variational,
        None,
        cancelled,
    )
}

/// The same three-grid stability assessment with the COMPLETE continuous
/// Hessian action of the original nonlinear Hamiltonian. The supplied action
/// follows [`PhsStorageHessianAction`], including material-memory contributions.
/// Unavailable, incomplete, refusing or nonfinite tangents never silently fall
/// back to finite differences or a quadratic replacement of nonlinear storage.
#[allow(clippy::too_many_arguments)]
pub fn analyze_forced_orbit_stability_variational_analytic(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    candidate: &HbOrbit,
    controls: &PhsFloquetControls,
    variational: &VariationalShootControls,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsVariationalForcedStability, PhsOrbitError> {
    analyze(
        system,
        input,
        candidate,
        controls,
        variational,
        Some(hessian),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn analyze(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    candidate: &HbOrbit,
    controls: &PhsFloquetControls,
    variational: &VariationalShootControls,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhsVariationalForcedStability, PhsOrbitError> {
    let mut poll = || {
        if cancelled() {
            Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
        } else {
            Ok(())
        }
    };
    poll()?;
    let (seed, period, grids) = admit_candidate(system, input, candidate, controls)?;
    let mut levels = Vec::with_capacity(3);
    for steps_per_period in grids {
        let budget = ShootBudget {
            steps_per_period,
            ..controls.shooting
        };
        let orbit = match hessian {
            Some(action) => solve_forced_shooting_variational_analytic(
                system,
                input,
                &seed,
                &budget,
                variational,
                action,
                || poll().is_err(),
            ),
            None => solve_forced_shooting_variational(
                system,
                input,
                &seed,
                &budget,
                variational,
                || poll().is_err(),
            ),
        }?;
        poll()?;
        let spectral_radius = check_orbit(&orbit.orbit, &seed, period, controls)?;
        levels.push(PhsVariationalFloquetLevel {
            steps_per_period,
            orbit,
            spectral_radius,
        });
    }
    let views = [
        (&levels[0].orbit.orbit, levels[0].spectral_radius),
        (&levels[1].orbit.orbit, levels[1].spectral_radius),
        (&levels[2].orbit.orbit, levels[2].spectral_radius),
    ];
    let (classification, dominant_growth_rate_per_s) = classify(&views, period, controls);
    Ok(PhsVariationalForcedStability {
        classification,
        levels,
        dominant_growth_rate_per_s,
    })
}

#[cfg(test)]
mod tests;
