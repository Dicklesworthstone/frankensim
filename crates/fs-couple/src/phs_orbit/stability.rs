//! Shooting and numerical Floquet stability of original forced pHS branches.

use super::{PeriodicPortInput, PhsOrbitError, bad, forced_problem};
use fs_math::{c64::C64, det};
use fs_orbit::{HbOrbit, ShootBudget, ShootOrbit, solve_shooting};
use fs_phs::PortHamiltonian;

mod subharmonic;
pub(super) use subharmonic::{gcd, repeat_rms};
mod variational;
pub use variational::{
    PhsVariationalFloquetLevel, PhsVariationalForcedStability,
    analyze_forced_orbit_stability_variational,
    analyze_forced_orbit_stability_variational_analytic,
};

/// Explicit physical agreement and numerical resolution controls.
#[derive(Debug, Clone)]
pub struct PhsFloquetControls {
    /// Native shooting budget on the first grid. The next grids use twice and
    /// four times its steps per forcing period; Newton controls stay unchanged.
    pub shooting: ShootBudget,
    /// Positive finite absolute tolerances, one per ORIGINAL state coordinate.
    /// Every shot must agree with the HB phase-zero state within these limits.
    /// Successive time grids must also agree to classify the result as resolved.
    /// Displacement, momentum and material-memory coordinates keep their units.
    pub state_tolerances: Vec<f64>,
    /// Positive finite maximum change of spectral radius between adjacent grids.
    pub spectral_radius_tolerance: f64,
    /// Dimensionless exclusion band around |multiplier|=1, strictly in (0,1).
    /// This numerical margin is not an enclosure of integration or spectral error.
    pub unit_circle_guard: f64,
}

/// Local linear stability inferred from measured, resolved numerical maps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhsNumericalStability {
    /// All three spectral radii are below 1 minus the declared guard.
    Stable,
    /// All three spectral radii exceed 1 plus the declared guard.
    Unstable,
    /// A grid comparison failed or a radius lies in the guarded transition.
    Unresolved,
}

/// One native shooting result; no neutral multiplier is removed for forcing.
#[derive(Debug, Clone)]
pub struct PhsFloquetLevel {
    /// Fixed RK4 steps used on this level's forcing period.
    pub steps_per_period: usize,
    /// Native phase-zero state, period, complete multiplier list and residuals.
    pub orbit: ShootOrbit,
    /// Largest multiplier magnitude, including every original state direction.
    pub spectral_radius: f64,
}

/// Stability evidence for the supplied forced branch at its authored phase.
#[derive(Debug, Clone)]
pub struct PhsForcedStability {
    /// Numerical local classification only after both grid comparisons pass.
    pub classification: PhsNumericalStability,
    /// Exactly three levels, ordered by increasing steps per period.
    pub levels: Vec<PhsFloquetLevel>,
    /// ln(finest spectral radius)/period, in 1/s; positive means growth.
    /// None when a zero/underflowed radius or division prevents a finite value.
    pub dominant_growth_rate_per_s: Option<f64>,
}

/// Independently shoot a forced HB/continuation point and resolve its local
/// numerical Floquet stability on three successively finer time grids.
///
/// `input` and `candidate` must share the same fundamental and phase origin.
/// For a frequency branch, first use [`PeriodicPortInput::at_frequency`] with
/// the accepted point's frequency. Each solve starts at the candidate's full
/// phase-zero state, uses the ORIGINAL pHS ODE, and fixes T=2*pi/omega. No time
/// shift, autonomous phase quotient, deleted damping or memory reduction is
/// allowed. A correction outside the supplied physical state tolerances
/// refuses instead of assigning another orbit's stability to this point.
///
/// fs-orbit owns RK4, Newton, the finite-difference monodromy and its eigensolve.
/// The first time grid must resolve the authored retained harmonic count at
/// Nyquist; the comparison checks measured state and spectral-radius changes,
/// not unrepresented nonlinear harmonics or a rigorous integration-error bound.
/// The native monodromy perturbation remains fixed. This is numerical evidence,
/// not a certified orbit, eigenvalue enclosure, global/nonlinear stability
/// theorem, or identification of a fold/flip/torus bifurcation. Native suspected
/// torus and source-domain refusals remain errors. Smooth ODEs only: switching
/// contacts and hybrid events need an event-aware map outside this path's scope.
/// No Cx cancellation is offered
/// by the existing synchronous shooting owner.
///
/// # Errors
/// Invalid settings, incompatible forcing/seed, source refusal, native shooting
/// failure, nonfinite results, or [`PhsOrbitError::SeedDisagreement`]. A completed
/// but insufficiently resolved ladder returns `Unresolved`, with all its data.
pub fn analyze_forced_orbit_stability(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    candidate: &HbOrbit,
    controls: &PhsFloquetControls,
) -> Result<PhsForcedStability, PhsOrbitError> {
    let (seed, period, grids) = admit_candidate(system, input, candidate, controls)?;
    let problem = forced_problem(system, input)?;
    let mut levels = Vec::with_capacity(3);
    for steps_per_period in grids {
        let budget = ShootBudget {
            steps_per_period,
            ..controls.shooting
        };
        let result = solve_shooting(&problem, &seed, period, &budget);
        if problem.source_refused() {
            return Err(PhsOrbitError::SourceRefused);
        }
        let orbit = result.map_err(PhsOrbitError::Orbit)?;
        let spectral_radius = check_orbit(&orbit, &seed, period, controls)?;
        levels.push(PhsFloquetLevel {
            steps_per_period,
            orbit,
            spectral_radius,
        });
    }
    let views = [
        (&levels[0].orbit, levels[0].spectral_radius),
        (&levels[1].orbit, levels[1].spectral_radius),
        (&levels[2].orbit, levels[2].spectral_radius),
    ];
    let (classification, dominant_growth_rate_per_s) = classify(&views, period, controls);
    Ok(PhsForcedStability {
        classification,
        levels,
        dominant_growth_rate_per_s,
    })
}

fn admit_candidate(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    candidate: &HbOrbit,
    controls: &PhsFloquetControls,
) -> Result<(Vec<f64>, f64, [usize; 3]), PhsOrbitError> {
    let n = system.state_dim();
    let width = candidate.coeffs.first().map_or(0, Vec::len);
    let period = core::f64::consts::TAU / candidate.omega;
    let steps = controls.shooting.steps_per_period;
    let finest_steps = steps
        .checked_mul(4)
        .ok_or_else(|| bad("shooting refinement step count overflow"))?;
    if n == 0
        || system.port_dim() != input.port_count()
        || candidate.coeffs.len() != n
        || width < 2
        || candidate.omega.to_bits() != input.omega_rad_s().to_bits()
        || !period.is_finite()
        || period <= 0.0
        || candidate.coeffs.iter().any(|row| {
            row.len() != width
                || row[0].im != 0.0
                || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        })
    {
        return Err(bad(
            "shooting needs matching finite original-state Fourier rows and forcing clock",
        ));
    }
    if input.highest_nonzero_harmonic() >= width {
        return Err(bad(
            "nonzero forcing harmonic exceeds the candidate HB truncation",
        ));
    }
    if steps == 0
        || (width - 1) > (steps - 1) / 2
        || controls.shooting.max_newton == 0
        || !controls.shooting.tolerance.is_finite()
        || controls.shooting.tolerance <= 0.0
        || controls.state_tolerances.len() != n
        || controls
            .state_tolerances
            .iter()
            .any(|t| !t.is_finite() || *t <= 0.0)
        || !controls.spectral_radius_tolerance.is_finite()
        || controls.spectral_radius_tolerance <= 0.0
        || !controls.unit_circle_guard.is_finite()
        || controls.unit_circle_guard <= 0.0
        || controls.unit_circle_guard >= 1.0
    {
        return Err(bad(
            "finite positive shooting controls and a time grid resolving the candidate harmonics required",
        ));
    }
    let seed: Vec<f64> = (0..n).map(|i| candidate.sample(i, 0.0)).collect();
    if seed.iter().any(|value| !value.is_finite()) {
        return Err(bad("candidate phase-zero state overflow"));
    }
    Ok((seed, period, [steps, finest_steps / 2, finest_steps]))
}

fn check_orbit(
    orbit: &ShootOrbit,
    seed: &[f64],
    period: f64,
    controls: &PhsFloquetControls,
) -> Result<f64, PhsOrbitError> {
    let n = seed.len();
    if orbit.x0.len() != n
        || orbit.multipliers.len() != n
        || orbit.period.to_bits() != period.to_bits()
        || orbit.x0.iter().any(|x| !x.is_finite())
        || orbit.residual_trace.is_empty()
        || orbit
            .residual_trace
            .iter()
            .any(|r| !r.is_finite() || *r < 0.0)
        || orbit
            .multipliers
            .iter()
            .any(|(re, im)| !re.is_finite() || !im.is_finite())
    {
        return Err(bad(
            "shooting returned an invalid physical orbit or Floquet spectrum",
        ));
    }
    for (coordinate, ((x, original), &tolerance)) in orbit
        .x0
        .iter()
        .zip(seed)
        .zip(&controls.state_tolerances)
        .enumerate()
    {
        let absolute_change = (x - original).abs();
        if absolute_change > tolerance {
            return Err(PhsOrbitError::SeedDisagreement {
                coordinate,
                absolute_change,
                tolerance,
            });
        }
    }
    let spectral_radius = orbit
        .multipliers
        .iter()
        .map(|&(re, im)| C64::new(re, im).abs())
        .fold(0.0, f64::max);
    if !spectral_radius.is_finite() {
        return Err(bad("Floquet spectral radius overflow"));
    }
    Ok(spectral_radius)
}

fn classify(
    levels: &[(&ShootOrbit, f64); 3],
    period: f64,
    controls: &PhsFloquetControls,
) -> (PhsNumericalStability, Option<f64>) {
    let resolved = levels.windows(2).all(|pair| {
        (pair[1].1 - pair[0].1).abs() <= controls.spectral_radius_tolerance
            && pair[1]
                .0
                .x0
                .iter()
                .zip(&pair[0].0.x0)
                .zip(&controls.state_tolerances)
                .all(|((a, b), tolerance)| (a - b).abs() <= *tolerance)
    });
    let classification = if resolved
        && levels
            .iter()
            .all(|level| level.1 < 1.0 - controls.unit_circle_guard)
    {
        PhsNumericalStability::Stable
    } else if resolved
        && levels
            .iter()
            .all(|level| level.1 > 1.0 + controls.unit_circle_guard)
    {
        PhsNumericalStability::Unstable
    } else {
        PhsNumericalStability::Unresolved
    };
    let radius = levels[2].1;
    let growth = if radius > 0.0 {
        det::ln(radius) / period
    } else {
        f64::NAN
    };
    (classification, growth.is_finite().then_some(growth))
}

#[cfg(test)]
mod tests;
