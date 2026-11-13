//! Physical output and resolved numerical stability along a complete forced branch.

use super::{PhsScaledFrequencyResponsePoint, admit_response};
use crate::phs_orbit::{
    PeriodicPortInput, PhsFloquetControls, PhsOrbitError, PhsOrbitProblem, PhsStorageHessianAction,
    PhsVariationalForcedStability, analyze_forced_orbit_stability_variational,
    analyze_forced_orbit_stability_variational_analytic, audit_orbit_with_input, forced_problem,
    spectrum::{PortSpectrum, PortSpectrumControls, PortSpectrumError, observe_port_spectrum},
    tangent::StorageTangent,
};
use core::cell::{Cell, RefCell};
use fs_math::c64::C64;
use fs_orbit::{
    ContinuableProblem, ContinuationBudget, ContinuationScales, ForcedParameterJacobian,
    ForcedParameterPoint, HbScales, IslandStateJacobian, OrbitProblem, VariationalShootControls,
    continue_branch_scaled_seeded,
};
use fs_phs::PortHamiltonian;

/// Existing original-state Floquet owners, applied independently at every point.
#[derive(Clone, Debug)]
pub struct PhsForcedBranchStabilityControls {
    /// Three time grids, physical seed-agreement limits and unit-circle guard.
    pub floquet: PhsFloquetControls,
    /// Physical tangent scales and work/space bounds for EACH time-grid solve.
    pub variational: VariationalShootControls,
}

/// Physical continuation, observation and optional stability resolution.
#[derive(Clone, Debug)]
pub struct PhsForcedPortBranchControls {
    /// Original pseudo-arclength/Newton/AFT budgets. The frequency may reverse.
    pub continuation: ContinuationBudget,
    /// Complete original-state and equation scales; no memory row is dropped.
    pub harmonic_balance_scales: HbScales,
    /// Original Hamiltonian, signed port supply and resistance audit samples.
    pub audit_samples: usize,
    /// Refined output spectrum of the original nonlinear effort law. Its order
    /// is independent of the retained state order and includes every port.
    pub spectrum: PortSpectrumControls,
    /// None makes no stability assessment. Some retains the complete existing
    /// three-grid variational result, including unresolved classifications.
    pub stability: Option<PhsForcedBranchStabilityControls>,
    /// Conservative retained numerical-result payload, admitted for at most
    /// `continuation.max_steps + 1` points before source work. Includes all
    /// state/output coefficients, traces and optional three-grid Floquet maps.
    /// Excludes allocator metadata, caller inputs and the existing owners'
    /// transient workspaces. Their own work/space limits remain in force.
    pub maximum_retained_bytes: usize,
}

/// One accepted original forced branch point, in continuation traversal order.
#[derive(Clone, Debug)]
pub struct PhsForcedPortBranchPoint {
    /// Physical Fourier state and clock, normalized continuation defects and
    /// the independent original cycle-work account. Every memory is retained.
    pub response: PhsScaledFrequencyResponsePoint,
    /// `G^T gradient(H)` with DC and `exp(-i*k*omega*t)` peak phasors.
    /// Nonlinear output harmonics are evaluated from the original storage.
    pub ports: PortSpectrum,
    /// Numerical stability of this point at the SAME authored source phase.
    /// Time-grid corrections must satisfy the original per-state agreement
    /// limits; no neighboring point supplies a shooting seed or classification.
    pub stability: Option<PhsVariationalForcedStability>,
}

/// Original owner failures retain their stage and accepted-point index.
#[derive(Clone, Debug, PartialEq)]
pub enum PhsForcedPortBranchError {
    Invalid(&'static str),
    Budget {
        requested_bytes: u128,
        maximum_bytes: usize,
    },
    Continuation(PhsOrbitError),
    Audit {
        point: usize,
        error: PhsOrbitError,
    },
    Spectrum {
        point: usize,
        error: PortSpectrumError,
    },
    Stability {
        point: usize,
        error: PhsOrbitError,
    },
    Cancelled,
}
impl core::fmt::Display for PhsForcedPortBranchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original forced port branch refused: {self:?}")
    }
}
impl std::error::Error for PhsForcedPortBranchError {}

type Result<T> = core::result::Result<T, PhsForcedPortBranchError>;

/// Continue a complete forced response through frequency folds, observing the
/// original physical port spectrum and optionally resolving Floquet stability.
///
/// The shared pseudo-arclength owner solves all Fourier states and frequency
/// together. Returned points retain that traversal order; no sort, frequency
/// interpolation, frequency-graph inversion or automatic stable-branch choice
/// is made. The original DC/complex input coefficients and phase origin stay
/// fixed while the common fundamental changes. The complete continuous storage
/// tangent is required and discovered from the original storage; use the
/// analytic entry for a storage whose tangent is supplied by the caller.
///
/// Each point gets its original cycle audit and independently refined port
/// spectrum. Optional stability uses the existing three-grid variational
/// shooting owner, starting every grid from THAT point's full phase-zero state.
/// It includes every original multiplier; no forced neutral direction is
/// removed. The returned HB state is not replaced by a corrected shooting state.
///
/// Fixed state-order and measured output/time-grid evidence do not establish
/// infinite harmonic completeness, unique branch identity, exact fold locations,
/// a physical hysteresis trajectory or global stability. In particular a local
/// unresolved classification is preserved. The corrector's frequency column
/// is the exact `i*n*X_n` derivative for all original coordinates at fixed
/// authored phase amplitudes; its state columns use the complete declared
/// original Hessian, not a quadratic approximation to nonlinear H.
///
/// Work is bounded by the original continuation budget, plus at most
/// `max_steps + 1` cycle audits/output queries and, when requested, three native
/// variational solves per point. Continuation polls cancellation at original
/// AFT source and tangent evaluations; audits poll per sample, and spectrum and
/// shooting use their existing cancellation boundaries. Shared dense solves,
/// spectral extractions and single source calls remain atomic and bounded by
/// their admitted dimensions. Cancellation/refusal returns no partial branch.
pub fn continue_forced_port_response_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    omega_end_rad_s: f64,
    seed: Vec<f64>,
    controls: &PhsForcedPortBranchControls,
    cancelled: impl FnMut() -> bool,
) -> Result<Vec<PhsForcedPortBranchPoint>> {
    continue_response(
        system,
        input,
        omega_end_rad_s,
        seed,
        controls,
        None,
        cancelled,
    )
}

/// The same branch with an explicit COMPLETE continuous Hessian action.
///
/// The same action differentiates the original storage in both the harmonic
/// corrector and optional variational shooting. It must include all nonlinear
/// couplings and material-memory directions. Refusal never falls back to a
/// different tangent or assigns another orbit's stability to this branch.
pub fn continue_forced_port_response_analytic_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    omega_end_rad_s: f64,
    seed: Vec<f64>,
    controls: &PhsForcedPortBranchControls,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<Vec<PhsForcedPortBranchPoint>> {
    continue_response(
        system,
        input,
        omega_end_rad_s,
        seed,
        controls,
        Some(hessian),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn continue_response(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    omega_end_rad_s: f64,
    seed: Vec<f64>,
    controls: &PhsForcedPortBranchControls,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<PhsForcedPortBranchPoint>> {
    poll(&mut cancelled)?;
    if system.port_dim() == 0 {
        return Err(PhsForcedPortBranchError::Invalid(
            "at least one original physical port is required for a response branch",
        ));
    }
    admit_response(
        input,
        omega_end_rad_s,
        None,
        &controls.continuation,
        controls.audit_samples,
    )
    .map_err(PhsForcedPortBranchError::Continuation)?;
    if input.highest_nonzero_harmonic() > (controls.audit_samples - 1) / 2 {
        return Err(PhsForcedPortBranchError::Invalid(
            "cycle audit must resolve every authored forcing harmonic",
        ));
    }
    admit_retained_payload(system, controls)?;
    let mut original =
        forced_problem(system, input).map_err(PhsForcedPortBranchError::Continuation)?;
    if let Some(action) = hessian {
        original.tangent = StorageTangent::hessian(system.state_dim(), action);
    }
    if !original.tangent.available() {
        return Err(PhsForcedPortBranchError::Invalid(
            "complete continuous original storage tangent required for the branch",
        ));
    }
    let points = {
        let mut problem = CancellableBranch {
            original,
            cancelled: RefCell::new(&mut cancelled),
            stopped: Cell::new(false),
        };
        let result = continue_branch_scaled_seeded(
            &mut problem,
            |omega| omega,
            input.omega_rad_s(),
            omega_end_rad_s,
            seed,
            &controls.continuation,
            &ContinuationScales {
                hb: controls.harmonic_balance_scales.clone(),
                parameter: controls.harmonic_balance_scales.angular_frequency,
            },
        );
        if problem.stopped.get() {
            return Err(PhsForcedPortBranchError::Cancelled);
        }
        if problem.original.source_refused() {
            return Err(PhsForcedPortBranchError::Continuation(
                PhsOrbitError::SourceRefused,
            ));
        }
        result
            .map_err(PhsOrbitError::Orbit)
            .map_err(PhsForcedPortBranchError::Continuation)?
    };
    poll(&mut cancelled)?;
    let mut response = Vec::with_capacity(points.len());
    for (index, point) in points.into_iter().enumerate() {
        poll(&mut cancelled)?;
        if point.point.lambda.to_bits() != point.point.orbit.omega.to_bits() {
            return Err(PhsForcedPortBranchError::Invalid(
                "continuation frequency and original forcing clock disagree",
            ));
        }
        let drive = input.at_frequency(point.point.lambda).map_err(|error| {
            PhsForcedPortBranchError::Audit {
                point: index,
                error,
            }
        })?;
        let mut stopped = false;
        let audit = audit_orbit_with_input(
            system,
            &point.point.orbit,
            controls.audit_samples,
            |time, _, _, out| {
                stopped = cancelled();
                !stopped && drive.sample_into(time, out)
            },
        );
        if stopped {
            return Err(PhsForcedPortBranchError::Cancelled);
        }
        let audit = audit.map_err(|error| PhsForcedPortBranchError::Audit {
            point: index,
            error,
        })?;
        let ports = observe_port_spectrum(
            system,
            &point.point.orbit,
            &controls.spectrum,
            &mut cancelled,
        );
        let ports = match ports {
            Err(PortSpectrumError::Cancelled) => return Err(PhsForcedPortBranchError::Cancelled),
            other => other.map_err(|error| PhsForcedPortBranchError::Spectrum {
                point: index,
                error,
            })?,
        };
        let stability = if let Some(settings) = &controls.stability {
            let result = match hessian {
                Some(action) => analyze_forced_orbit_stability_variational_analytic(
                    system,
                    &drive,
                    &point.point.orbit,
                    &settings.floquet,
                    &settings.variational,
                    action,
                    &mut cancelled,
                ),
                None => analyze_forced_orbit_stability_variational(
                    system,
                    &drive,
                    &point.point.orbit,
                    &settings.floquet,
                    &settings.variational,
                    &mut cancelled,
                ),
            };
            Some(match result {
                Err(PhsOrbitError::Variational(fs_orbit::VariationalShootError::Cancelled)) => {
                    return Err(PhsForcedPortBranchError::Cancelled);
                }
                other => other.map_err(|error| PhsForcedPortBranchError::Stability {
                    point: index,
                    error,
                })?,
            })
        } else {
            None
        };
        poll(&mut cancelled)?;
        response.push(PhsForcedPortBranchPoint {
            response: PhsScaledFrequencyResponsePoint { point, audit },
            ports,
            stability,
        });
    }
    Ok(response)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PhsForcedPortBranchError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_retained_payload(system: &PortHamiltonian, c: &PhsForcedPortBranchControls) -> Result<()> {
    let n = system.state_dim() as u128;
    let m = system.port_dim() as u128;
    // Complete complex state rows (including their real DC slot), native HB
    // Newton trace, all port lines, and the observer's at-most-12-grid trace.
    let orbit = 16 * n * (c.continuation.hb.harmonics as u128 + 1)
        + 8 * (c.continuation.hb.max_newton as u128 + 2);
    let ports = 16 * m * (c.spectrum.harmonics as u128 + 1) + 8 * (32 * m + 16);
    let stability = c.stability.as_ref().map_or(0, |s| {
        // Per grid: full monodromy, physical seed, scales, every complex
        // multiplier and Newton trace. Forced shooting has no period column.
        3 * (8 * n * n + 32 * n + 8 * (s.floquet.shooting.max_newton as u128 + 1) + 256)
    });
    let per_point = orbit + ports + stability + 512;
    let requested = per_point
        .checked_mul(c.continuation.max_steps as u128 + 1)
        .ok_or(PhsForcedPortBranchError::Invalid(
            "retained branch payload extent overflow",
        ))?;
    if requested > c.maximum_retained_bytes as u128 {
        return Err(PhsForcedPortBranchError::Budget {
            requested_bytes: requested,
            maximum_bytes: c.maximum_retained_bytes,
        });
    }
    Ok(())
}

/// Only supplies cancellation boundaries to the existing original ODE adapter.
/// Its physical port and exact state tangent keep their original definitions.
struct CancellableBranch<'a, 'c> {
    original: PhsOrbitProblem<'a>,
    cancelled: RefCell<&'c mut dyn FnMut() -> bool>,
    stopped: Cell<bool>,
}
impl CancellableBranch<'_, '_> {
    fn stopped(&self) -> bool {
        let stopped = self.stopped.get() || (self.cancelled.borrow_mut())();
        self.stopped.set(stopped);
        stopped
    }
    fn state_jacobian(&self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        !self.stopped() && self.original.state_jacobian(time, state, out)
    }
    fn parameter_jacobian(&self, point: ForcedParameterPoint<'_>, out: &mut [f64]) -> bool {
        !self.stopped() && self.original.frequency_parameter_jacobian(point, out)
    }
}
impl OrbitProblem for CancellableBranch<'_, '_> {
    fn dim(&self) -> usize {
        self.original.dim()
    }
    fn autonomous(&self) -> bool {
        self.original.autonomous()
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        if self.stopped() {
            out.fill(f64::NAN);
        } else {
            self.original.island(time, state, out);
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        self.original.port(s)
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.original
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
}
impl ContinuableProblem for CancellableBranch<'_, '_> {
    fn set_parameter(&mut self, lambda: f64) {
        self.original.set_parameter(lambda);
    }
    fn forced_parameter_jacobian(&self) -> Option<ForcedParameterJacobian<Self>> {
        self.original
            .forced_parameter_jacobian()
            .is_some()
            .then_some(Self::parameter_jacobian)
    }
}
