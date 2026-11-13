//! Native autonomous shooting with an original-coordinate transverse section.

use super::{AutonomousPortFeedback, audit_autonomous_feedback_orbit, feedback_problem};
use crate::phs_orbit::{
    PhsFloquetControls, PhsNumericalStability, PhsOrbitAudit, PhsOrbitError, PhsOrbitProblem, bad,
};
use core::cell::RefCell;
use fs_math::{c64::C64, det};
use fs_orbit::{HbOrbit, OrbitProblem, ShootBudget, ShootOrbit, TimeDomainRhs, solve_shooting};
use fs_phs::PortHamiltonian;

/// Physical section, period and Floquet resolution for an autonomous orbit.
#[derive(Debug, Clone)]
pub struct PhsAutonomousControls {
    /// Three-grid shooting and original-coordinate state/radius controls.
    pub floquet: PhsFloquetControls,
    /// Positive relative period agreement with HB and between grids, below 1/2.
    pub relative_period_tolerance: f64,
    /// Distance from +1 admitted for the single phase multiplier; positive and
    /// smaller than floquet.unit_circle_guard, which also separates other modes.
    pub phase_multiplier_tolerance: f64,
    /// Original source-work quadrature samples on the supplied HB waveform.
    pub audit_samples: usize,
}

/// One native autonomous shooting result in ORIGINAL physical coordinates.
#[derive(Debug, Clone)]
pub struct PhsAutonomousLevel {
    /// Fixed RK4 steps used per solved period on this level.
    pub steps_per_period: usize,
    /// Includes every multiplier. A phase multiplier is never deleted.
    pub orbit: ShootOrbit,
    /// Present only for one isolated multiplier near +1 on a transverse orbit.
    /// Multiple neutral directions, including conservative families, give None.
    pub phase_multiplier_index: Option<usize>,
    /// Largest modulus after excluding that one identified phase direction.
    pub transverse_spectral_radius: Option<f64>,
}

/// Numerically resolved autonomous cycle and its original physical source work.
#[derive(Debug, Clone)]
pub struct PhsAutonomousStability {
    /// Local numerical transverse classification after the three-grid checks.
    pub classification: PhsNumericalStability,
    /// Original coordinate held at the candidate's phase-zero value. The native
    /// coordinate-0 chart is only a permutation; no phase/time shift is hidden.
    pub section_coordinate: usize,
    /// Candidate phase-zero value in that coordinate's original units.
    pub section_value: f64,
    /// Exactly three time grids; every x0 retains the original state ordering.
    pub levels: Vec<PhsAutonomousLevel>,
    /// Finest transverse log-modulus / period, when a phase mode is isolated.
    pub dominant_transverse_growth_rate_per_s: Option<f64>,
    /// Original H, feedback supply and R loss on the supplied HB waveform. This
    /// is not a time-domain quadrature of the shooting trajectories.
    pub candidate_audit: PhsOrbitAudit,
}

// Native shooting freezes coordinate 0. Put a genuinely flowing coordinate
// there, rather than freezing displacement at a turning point. Swapping both
// states and derivatives is a similarity change, not a changed physical model.
struct SectionProblem<'a> {
    original: PhsOrbitProblem<'a>,
    coordinate: usize,
    state: RefCell<Vec<f64>>,
    derivative: RefCell<Vec<f64>>,
}
impl OrbitProblem for SectionProblem<'_> {
    fn dim(&self) -> usize {
        self.original.dim()
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        if state.len() != self.dim() || out.len() != self.dim() {
            self.original.fail(out);
            return;
        }
        let mut physical = self.state.borrow_mut();
        physical.copy_from_slice(state);
        physical.swap(0, self.coordinate);
        let mut derivative = self.derivative.borrow_mut();
        self.original.island(time, &physical, &mut derivative);
        out.copy_from_slice(&derivative);
        out.swap(0, self.coordinate);
    }
}

/// Resolve a smooth self-excited orbit against the complete original ODE.
///
/// The strongest tolerance-scaled physical flow at HB phase zero supplies a
/// transverse coordinate section. Native shooting solves for the period; every
/// state is returned in original coordinates and must agree with that same HB
/// section. Orientation and resolvable section speed are rechecked after each
/// shot. Source work is independently audited on the original HB waveform.
///
/// An autonomous nonstationary flow has a phase-neutral tangent. Its numerical
/// multiplier is identified only when exactly one lies near +1 and all others
/// are separated from +1 by the declared guard. Conservative families, internal
/// resonance degeneracies and unresolved phase clusters remain Unresolved, or
/// retain native Newton/torus refusal. No amplitude constraint is invented for
/// a conservative family, and no stable-cycle label follows from convergence.
///
/// This reuses native RK4, forward-difference monodromy and eigensolve on three
/// grids. It is not a certified quotient, bifurcation classifier, global
/// stability theorem, hybrid-event map or fundamental-period proof. The callback
/// is an explicit external source law, not validated material/flow physics.
/// These bounded synchronous native solves have no cooperative Cx cancellation.
///
/// # Errors
/// Invalid controls or waveform, stationary/unresolved section, source refusal,
/// native solver failure, or state/period disagreement. Conservative families
/// can refuse at native Newton rather than producing a resolved quotient.
pub fn analyze_autonomous_feedback_stability(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    candidate: &HbOrbit,
    controls: &PhsAutonomousControls,
) -> Result<PhsAutonomousStability, PhsOrbitError> {
    let n = system.state_dim();
    let width = candidate.coeffs.first().map_or(0, Vec::len);
    let period = core::f64::consts::TAU / candidate.omega;
    let c = &controls.floquet;
    let steps = c.shooting.steps_per_period;
    let finest = steps
        .checked_mul(4)
        .ok_or_else(|| bad("shooting refinement step count overflow"))?;
    if n < 2
        || feedback.port_count() != system.port_dim()
        || candidate.coeffs.len() != n
        || width < 2
        || !candidate.omega.is_finite()
        || candidate.omega <= 0.0
        || !period.is_finite()
        || period <= 0.0
        || candidate.coeffs.iter().any(|row| {
            row.len() != width
                || row[0].im != 0.0
                || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        })
        || candidate
            .coeffs
            .iter()
            .all(|row| row.iter().skip(1).all(|z| z.re == 0.0 && z.im == 0.0))
    {
        return Err(bad(
            "autonomous shooting needs a finite nonconstant original-state HB waveform",
        ));
    }
    if steps == 0
        || (width - 1) > (steps - 1) / 2
        || c.shooting.max_newton == 0
        || !c.shooting.tolerance.is_finite()
        || c.shooting.tolerance <= 0.0
        || c.state_tolerances.len() != n
        || c.state_tolerances
            .iter()
            .any(|t| !t.is_finite() || *t <= 0.0)
        || !c.spectral_radius_tolerance.is_finite()
        || c.spectral_radius_tolerance <= 0.0
        || !c.unit_circle_guard.is_finite()
        || c.unit_circle_guard <= 0.0
        || c.unit_circle_guard >= 1.0
        || !controls.relative_period_tolerance.is_finite()
        || controls.relative_period_tolerance <= 0.0
        || controls.relative_period_tolerance >= 0.5
        || !controls.phase_multiplier_tolerance.is_finite()
        || controls.phase_multiplier_tolerance <= 0.0
        || controls.phase_multiplier_tolerance >= c.unit_circle_guard
        || !(32..=65_536).contains(&controls.audit_samples)
    {
        return Err(bad(
            "invalid autonomous shooting, physical agreement or phase-isolation controls",
        ));
    }
    let seed: Vec<f64> = (0..n).map(|i| candidate.sample(i, 0.0)).collect();
    if seed.iter().any(|x| !x.is_finite()) {
        return Err(bad("autonomous phase-zero state overflow"));
    }
    let original = feedback_problem(system, feedback)?;
    let mut tangent = vec![0.0; n];
    original.island(0.0, &seed, &mut tangent);
    if original.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    let mut coordinate = 0;
    for i in 1..n {
        if tangent[i].abs() / c.state_tolerances[i]
            > tangent[coordinate].abs() / c.state_tolerances[coordinate]
        {
            coordinate = i;
        }
    }
    if tangent[coordinate].abs() / c.state_tolerances[coordinate] * period <= 1.0 {
        return Err(bad(
            "HB section is stationary or unresolved in the original state tolerances",
        ));
    }
    let direction = tangent[coordinate].is_sign_positive();
    let problem = SectionProblem {
        original,
        coordinate,
        state: RefCell::new(vec![0.0; n]),
        derivative: RefCell::new(vec![0.0; n]),
    };
    let mut chart_seed = seed.clone();
    chart_seed.swap(0, coordinate);
    let mut levels = Vec::with_capacity(3);
    for steps_per_period in [steps, finest / 2, finest] {
        let budget = ShootBudget {
            steps_per_period,
            ..c.shooting
        };
        let result = solve_shooting(&problem, &chart_seed, period, &budget);
        if problem.original.source_refused() {
            return Err(PhsOrbitError::SourceRefused);
        }
        let mut orbit = result.map_err(PhsOrbitError::Orbit)?;
        if orbit.x0.len() != n
            || orbit.multipliers.len() != n
            || !orbit.period.is_finite()
            || orbit.period <= 0.0
            || orbit.x0.iter().any(|x| !x.is_finite())
            || orbit.residual_trace.is_empty()
            || orbit
                .residual_trace
                .iter()
                .any(|r| !r.is_finite() || *r < 0.0)
            || orbit
                .multipliers
                .iter()
                .any(|&(re, im)| !C64::new(re, im).abs().is_finite())
        {
            return Err(bad(
                "autonomous shooting returned an invalid orbit or spectrum",
            ));
        }
        orbit.x0.swap(0, coordinate);
        if (orbit.period / period - 1.0).abs() > controls.relative_period_tolerance {
            return Err(bad(
                "autonomous shooting period disagrees with the supplied HB cycle",
            ));
        }
        for i in 0..n {
            let change = (orbit.x0[i] - seed[i]).abs();
            if change > c.state_tolerances[i] {
                return Err(PhsOrbitError::SeedDisagreement {
                    coordinate: i,
                    absolute_change: change,
                    tolerance: c.state_tolerances[i],
                });
            }
        }
        problem.original.island(0.0, &orbit.x0, &mut tangent);
        if problem.original.source_refused() {
            return Err(PhsOrbitError::SourceRefused);
        }
        if orbit.x0[coordinate].to_bits() != seed[coordinate].to_bits()
            || tangent[coordinate].is_sign_positive() != direction
            || tangent[coordinate].abs() / c.state_tolerances[coordinate] * orbit.period <= 1.0
        {
            return Err(bad(
                "autonomous shooting left the original transverse oriented section",
            ));
        }
        let phase: Vec<usize> = orbit
            .multipliers
            .iter()
            .enumerate()
            .filter_map(|(i, &(re, im))| {
                (C64::new(re - 1.0, im).abs() <= controls.phase_multiplier_tolerance).then_some(i)
            })
            .collect();
        let phase_multiplier_index = if phase.len() == 1
            && orbit.multipliers.iter().enumerate().all(|(i, &(re, im))| {
                i == phase[0] || C64::new(re - 1.0, im).abs() > c.unit_circle_guard
            }) {
            Some(phase[0])
        } else {
            None
        };
        let transverse_spectral_radius = phase_multiplier_index.map(|phase| {
            orbit
                .multipliers
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != phase)
                .map(|(_, &(re, im))| C64::new(re, im).abs())
                .fold(0.0, f64::max)
        });
        levels.push(PhsAutonomousLevel {
            steps_per_period,
            orbit,
            phase_multiplier_index,
            transverse_spectral_radius,
        });
    }
    let resolved = levels.iter().all(|l| l.phase_multiplier_index.is_some())
        && levels.windows(2).all(|pair| {
            (pair[1].orbit.period / pair[0].orbit.period - 1.0).abs()
                <= controls.relative_period_tolerance
                && pair[1]
                    .orbit
                    .x0
                    .iter()
                    .zip(&pair[0].orbit.x0)
                    .zip(&c.state_tolerances)
                    .all(|((a, b), t)| (a - b).abs() <= *t)
                && match (
                    pair[1].transverse_spectral_radius,
                    pair[0].transverse_spectral_radius,
                ) {
                    (Some(a), Some(b)) => (a - b).abs() <= c.spectral_radius_tolerance,
                    _ => false,
                }
        });
    let classification = if resolved
        && levels.iter().all(|l| {
            l.transverse_spectral_radius
                .is_some_and(|r| r < 1.0 - c.unit_circle_guard)
        }) {
        PhsNumericalStability::Stable
    } else if resolved
        && levels.iter().all(|l| {
            l.transverse_spectral_radius
                .is_some_and(|r| r > 1.0 + c.unit_circle_guard)
        })
    {
        PhsNumericalStability::Unstable
    } else {
        PhsNumericalStability::Unresolved
    };
    let growth = levels[2].transverse_spectral_radius.and_then(|radius| {
        let value = det::ln(radius) / levels[2].orbit.period;
        value.is_finite().then_some(value)
    });
    let candidate_audit =
        audit_autonomous_feedback_orbit(system, feedback, candidate, controls.audit_samples)?;
    Ok(PhsAutonomousStability {
        classification,
        section_coordinate: coordinate,
        section_value: seed[coordinate],
        levels,
        dominant_transverse_growth_rate_per_s: growth,
        candidate_audit,
    })
}
