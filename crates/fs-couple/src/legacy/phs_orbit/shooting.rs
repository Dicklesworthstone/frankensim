//! Complete original pHS responses on the shared variational shooting owner.

use super::{
    PeriodicPortInput, PhsOrbitError, PhsStorageHessianAction, bad, forced_problem,
    tangent::StorageTangent,
};
use fs_orbit::{
    ShootBudget, VariationalShootControls, VariationalShootError, VariationalShootOrbit,
    solve_shooting_variational,
};
use fs_phs::PortHamiltonian;

/// Shoot the complete forced pHS ODE using its declared original storage
/// tangent, including explicit material-memory coordinates. Complete exact
/// quadratic metadata or [`fs_phs::Storage::continuous_tangent`] is discovered
/// automatically. A quadratic part of a nonlinear Hamiltonian is insufficient
/// and returns TangentUnavailable. Other laws can supply an action through
/// [`solve_forced_shooting_variational_analytic`].
///
/// The physical period is exactly 2*pi/input.omega_rad_s(), with the authored
/// source phases unchanged. `seed` and state scales retain original displacement,
/// momentum, material-memory or other physical coordinates; no mass-normalized
/// replacement model is inferred. The first grid must resolve every nonzero
/// forcing harmonic. Source/storage work still follows the original H/R/G laws;
/// this method adds no damping and manufactures no cycle-average heat account.
///
/// fs-orbit owns every RK4 stage, its exact discrete tangent, Newton correction
/// and complete Floquet spectrum. This returns one numerical time grid, without
/// inferring continuous-orbit accuracy, a stability classification or a temporal
/// convergence certificate. Refusal/cancellation publishes no partial response.
/// Original pHS sources and callbacks must be deterministic and read-only.
///
/// The supplied numerical caps cover BOTH the adapter and shared shooting owner.
/// Before source preparation, reserve 16*n*n+64*n+32*m+256 payload bytes for
/// complete storage/tangent metadata and scratch. The additional work score is
/// 32*n*n + 4*maps*steps*(4*n^3+32*n*n+16*n*m+32*m*(harmonics+1));
/// harmonics counts authored slots, including zero padding. This bounds owner
/// contractions/source sampling, not arbitrary callback internals.
pub fn solve_forced_shooting_variational(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: &[f64],
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    solve(system, input, seed, budget, controls, None, cancelled)
}

/// The same physical forced response with an authored COMPLETE continuous
/// Hessian action of the original nonlinear Hamiltonian. The shared source owner
/// forms `(J-R)*D(grad H)` at each exact primal stage. Every material contribution
/// and explicit memory state must be included; no sampled quadratic surrogate,
/// trajectory difference, island/port inversion or second RK4 implementation.
///
/// The action has [`PhsStorageHessianAction`]'s original-coordinate contract.
/// Its absence, refusal or incomplete/nonfinite result cannot silently switch
/// to finite differences. At most n actions per native Jacobian evaluation are
/// requested; native stage/map caps bound those requests. The action's internal
/// work and memory are caller-owned. Source availability, scaling, bandwidth and
/// numerical limits have the same meaning as the automatic storage-owned entry.
pub fn solve_forced_shooting_variational_analytic(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: &[f64],
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    solve(
        system,
        input,
        seed,
        budget,
        controls,
        Some(hessian),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn solve(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: &[f64],
    budget: &ShootBudget,
    controls: &VariationalShootControls,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VariationalShootOrbit, PhsOrbitError> {
    if cancelled() {
        return Err(PhsOrbitError::Variational(VariationalShootError::Cancelled));
    }
    let n = system.state_dim();
    let m = system.port_dim();
    if n == 0
        || n > 256
        || n > controls.maximum_dimension
        || input.port_count() != m
        || seed.len() != n
        || seed.iter().any(|v| !v.is_finite())
        || controls.state_scales.len() != n
        || controls
            .state_scales
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        || budget.steps_per_period == 0
        || input.highest_nonzero_harmonic() > (budget.steps_per_period - 1) / 2
    {
        return Err(bad(
            "complete finite pHS seed, positive physical state scales and a shooting grid resolving every authored input harmonic required",
        ));
    }
    let period = core::f64::consts::TAU / input.omega_rad_s();
    if !period.is_finite() || period <= 0.0 {
        return Err(bad(
            "positive representable physical forcing period required",
        ));
    }
    let (n, m) = (n as u128, m as u128);
    let bytes = 16 * n * n + 64 * n + 32 * m + 256;
    let work = (controls.maximum_period_maps as u128)
        .checked_mul(budget.steps_per_period as u128)
        .and_then(|v| v.checked_mul(4))
        .and_then(|v| {
            v.checked_mul(
                4 * n * n * n
                    + 32 * n * n
                    + 16 * n * m
                    + 32 * m * (input.harmonic_count() as u128 + 1),
            )
        })
        .and_then(|v| v.checked_add(32 * n * n))
        .ok_or_else(|| bad("pHS variational work envelope overflows"))?;
    for (quantity, requested, cap) in [
        (
            "pHS variational payload bytes",
            bytes,
            controls.maximum_payload_bytes as u128,
        ),
        ("pHS variational work", work, controls.maximum_work as u128),
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
    let mut problem = forced_problem(system, input)?;
    if let Some(hessian) = hessian {
        problem.tangent = StorageTangent::hessian(system.state_dim(), hessian);
    }
    let result =
        solve_shooting_variational(&problem, seed, period, budget, &admitted, &mut cancelled);
    match result {
        Err(VariationalShootError::Cancelled) => {
            Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
        }
        _ if problem.source_refused() => Err(PhsOrbitError::SourceRefused),
        other => other.map_err(PhsOrbitError::Variational),
    }
}
