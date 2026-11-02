//! Full-waveform seeds and unchanged forcing on a multiple-period clock.
//!
//! These are adapters to the existing fs-orbit solvers, not new integrators.
//! An input can be represented over several of its original periods without
//! changing any physical forcing frequency, phase, DC value or port amplitude.

use crate::phs_orbit::{
    MAX_INPUT_HARMONICS, PeriodicPortInput, PhsFloquetControls, PhsForcedStability,
    PhsOrbitAudit, PhsOrbitError, analyze_forced_orbit_stability, audit_forced_orbit,
    bad, forced_problem,
};
use fs_math::{c64::C64, det};
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbOrbit, HbRefinementBudget,
    solve_hb_adaptive_seeded, solve_hb_seeded,
};
use fs_phs::PortHamiltonian;

impl PeriodicPortInput {
    /// Express the SAME physical forcing over `periods` original drive periods.
    /// The new fundamental is omega/periods and coefficient n moves to n*periods.
    /// Unlike `at_frequency`, this does not retune the physical drive. The phase
    /// origin, complex amplitudes, DC and original port ordering are unchanged.
    ///
    /// This permits period-doubled and higher subharmonic response searches on
    /// the existing forced HB, audit and shooting paths. A longer chart alone
    /// does NOT establish a longer-period response: an ordinary drive-periodic
    /// solution also fits it. Use a full waveform seed to select a branch.
    ///
    /// Exact-zero trailing padding does not consume bandwidth. Every nonzero
    /// coefficient, including subnormals, must fit the original 256-harmonic
    /// input cap after remapping. The multiplier is explicitly bounded to 256;
    /// unsupported bandwidth or an unrepresentable period refuses before any
    /// source evaluation. Powers-of-two remaps preserve binary frequencies
    /// exactly when representable; other integer remaps have ordinary rounding.
    pub fn over_periods(&self, periods: usize) -> Result<Self, PhsOrbitError> {
        if !(1..=MAX_INPUT_HARMONICS).contains(&periods) {
            return Err(bad("forcing period multiplier must be in 1..=256"));
        }
        let highest = self.highest_nonzero_harmonic();
        let remapped = highest
            .checked_mul(periods)
            .filter(|&n| n <= MAX_INPUT_HARMONICS)
            .ok_or_else(|| bad("multiple-period forcing exceeds the input harmonic limit"))?;
        let omega = self.omega_rad_s / periods as f64;
        let duration = core::f64::consts::TAU / omega;
        if !omega.is_finite() || omega <= 0.0 || !duration.is_finite() || duration <= 0.0 {
            return Err(bad("multiple-period forcing clock is not representable"));
        }
        if periods == 1 {
            return Ok(self.clone());
        }
        let mut coefficients = Vec::with_capacity(self.port_count());
        for original in &self.coefficients {
            let mut row = vec![C64::ZERO; remapped + 1];
            row[0] = original[0];
            for n in 1..=highest {
                row[n * periods] = original[n];
            }
            coefficients.push(row);
        }
        Self::new(omega, coefficients)
    }

    /// Solve a forced branch from every original-state Fourier coefficient.
    ///
    /// The seed must have this input's fundamental and phase origin, one row per
    /// original state (including material memory), and exactly budget.harmonics
    /// retained harmonics. No coefficient is dropped, phase-aligned or inferred
    /// from a scalar displacement amplitude. Previous seed convergence metadata
    /// is not evidence: the native owner solves and measures the original
    /// equations again. Newton may reach another branch; no branch-identity or
    /// fundamental-period certificate is implied.
    ///
    /// Native work limits, harmonic truncation and synchronous/no-Cx cancellation
    /// behavior remain unchanged. Source-domain refusals remain SourceRefused.
    pub fn solve_harmonic_balance_seeded(
        &self,
        system: &PortHamiltonian,
        seed: &HbOrbit,
        budget: &HbBudget,
    ) -> Result<HbOrbit, PhsOrbitError> {
        admit_seed(self, system, seed, budget)?;
        let problem = forced_problem(system, self)?;
        let anchor = HbAnchor::Forced { omega: self.omega_rad_s() };
        let result = solve_hb_seeded(&problem, anchor, seed.packed_seed(&anchor), budget);
        if problem.source_refused() {
            return Err(PhsOrbitError::SourceRefused);
        }
        result.map_err(PhsOrbitError::Orbit)
    }

    /// Refine a full seeded forced waveform with the native harmonic ladder.
    /// Starting order must exactly match the seed and contain every nonzero
    /// forcing harmonic. Larger solved orders retain the original native
    /// waveform, equation-defect and frequency checks, including all memory
    /// coordinates. No unresolved input harmonic is silently truncated.
    pub fn solve_harmonic_balance_adaptive_seeded(
        &self,
        system: &PortHamiltonian,
        seed: &HbOrbit,
        budget: &HbBudget,
        refinement: &HbRefinementBudget,
    ) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
        admit_seed(self, system, seed, budget)?;
        let problem = forced_problem(system, self)?;
        let anchor = HbAnchor::Forced { omega: self.omega_rad_s() };
        let result = solve_hb_adaptive_seeded(
            &problem, anchor, seed.packed_seed(&anchor), budget, refinement,
        );
        if problem.source_refused() {
            return Err(PhsOrbitError::SourceRefused);
        }
        result.map_err(PhsOrbitError::Orbit)
    }


    /// Refuse a purported period-m response indistinguishable from a shorter
    /// forcing-synchronized cycle at the declared ORIGINAL-state RMS tolerances.
    ///
    /// For every proper divisor r of m, compare x(t+r*T_drive) with x(t) using
    /// Parseval on the complete retained Fourier waveform. At least one original
    /// state must distinguish each shorter repeat. All material-memory states
    /// count; a displacement node does not erase a genuine internal oscillation.
    /// DC/preload does not rescale the tolerances and tiny numerical tails below
    /// those tolerances do not establish a longer period. No phase alignment or
    /// extra sampling grid can hide a shorter cycle.
    ///
    /// Call this on the original drive, not its over_periods chart. Its nonzero
    /// AC harmonic indices must have gcd 1 so the authored drive period is not
    /// itself an artificially padded period. This is retained-waveform separation,
    /// NOT an existence proof, a bound on omitted harmonics, or a proof of the
    /// fundamental period among every possible nonsynchronized time shift.
    pub fn require_subharmonic_separation(
        &self,
        system: &PortHamiltonian,
        periods: usize,
        candidate: &HbOrbit,
        state_rms_tolerances: &[f64],
    ) -> Result<(), PhsOrbitError> {
        admit_subharmonic(self, system, periods, candidate, state_rms_tolerances).map(|_| ())
    }

    /// Audit physical work over the complete response cycle m*T_drive.
    /// Original source amplitudes/phases, storage and resistance are unchanged.
    /// Shorter synchronized repeats must first be resolved in the complete state.
    /// As with audit_forced_orbit, repeat the quadrature at a finer sample count;
    /// one finite quadrature is not a broadband or omitted-harmonic certificate.
    pub fn audit_subharmonic_orbit(
        &self,
        system: &PortHamiltonian,
        periods: usize,
        candidate: &HbOrbit,
        state_rms_tolerances: &[f64],
        samples: usize,
    ) -> Result<PhsOrbitAudit, PhsOrbitError> {
        let chart = admit_subharmonic(self, system, periods, candidate, state_rms_tolerances)?;
        audit_forced_orbit(system, &chart, candidate, samples)
    }

    /// Shoot the complete subharmonic cycle on the original physical ODE.
    /// All three grids and returned multipliers span m*T_drive; the controls'
    /// steps_per_period is therefore steps per RESPONSE cycle, and the reported
    /// growth rate divides by that full physical duration. The orbit is forced:
    /// no autonomous phase multiplier is removed. A single drive-period map is
    /// not exponentiated as a substitute for integrating the actual full orbit.
    ///
    /// Separate RMS tolerances admit the candidate's period; the existing
    /// Floquet controls independently limit shooting corrections, time-grid
    /// disagreement and the guarded unit circle. This establishes numerical
    /// stability only under the existing native shooting limitations. It does
    /// not certify that a corrected trajectory has a mathematically primitive
    /// period, identify a period-doubling bifurcation, or cover hybrid events.
    pub fn analyze_subharmonic_orbit_stability(
        &self,
        system: &PortHamiltonian,
        periods: usize,
        candidate: &HbOrbit,
        state_rms_tolerances: &[f64],
        controls: &PhsFloquetControls,
    ) -> Result<PhsForcedStability, PhsOrbitError> {
        let chart = admit_subharmonic(self, system, periods, candidate, state_rms_tolerances)?;
        analyze_forced_orbit_stability(system, &chart, candidate, controls)
    }
}

fn admit_waveform(
    input: &PeriodicPortInput,
    dimension: usize,
    seed: &HbOrbit,
) -> Result<usize, PhsOrbitError> {
    let width = seed.coeffs.first().map_or(0, Vec::len);
    let duration = core::f64::consts::TAU / seed.omega;
    if dimension == 0
        || seed.coeffs.len() != dimension
        || width < 2
        || seed.omega.to_bits() != input.omega_rad_s().to_bits()
        || !duration.is_finite()
        || duration <= 0.0
        || seed.coeffs.iter().any(|row| {
            row.len() != width || row[0].im != 0.0
                || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        })
    {
        return Err(bad("seed must contain finite complete original-state rows on the forcing clock"));
    }
    if input.highest_nonzero_harmonic() >= width {
        return Err(bad("nonzero forcing harmonic exceeds the seeded HB truncation"));
    }
    Ok(width - 1)
}

fn admit_seed(
    input: &PeriodicPortInput,
    system: &PortHamiltonian,
    seed: &HbOrbit,
    budget: &HbBudget,
) -> Result<(), PhsOrbitError> {
    if budget.harmonics == 0
        || budget.max_newton == 0
        || !budget.tolerance.is_finite()
        || budget.tolerance <= 0.0
        || budget.harmonics.checked_mul(2).and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(system.state_dim()))
            .and_then(|n| n.checked_mul(n)).is_none()
        || system.port_dim() != input.port_count()
    {
        return Err(bad("seeded forced HB needs valid budgets and every original physical port"));
    }
    let harmonics = admit_waveform(input, system.state_dim(), seed)?;
    if harmonics != budget.harmonics {
        return Err(bad("seed harmonic count must equal the declared starting order"));
    }
    Ok(())
}


fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 { let remainder = a % b; a = b; b = remainder; }
    a
}

fn admit_subharmonic(
    drive: &PeriodicPortInput,
    system: &PortHamiltonian,
    periods: usize,
    candidate: &HbOrbit,
    tolerances: &[f64],
) -> Result<PeriodicPortInput, PhsOrbitError> {
    let chart = drive.over_periods(periods)?;
    admit_waveform(&chart, system.state_dim(), candidate)?;
    if drive.port_count() != system.port_dim()
        || tolerances.len() != system.state_dim()
        || tolerances.iter().any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(bad("subharmonic separation needs every port and positive RMS tolerances per original state"));
    }
    let mut common = 0;
    for row in &drive.coefficients {
        for (n, value) in row.iter().enumerate().skip(1) {
            if value.re != 0.0 || value.im != 0.0 { common = gcd(common, n); }
        }
    }
    if common != 1 {
        return Err(bad("subharmonic order requires a nonconstant primitive authored drive clock"));
    }
    for repeat in 1..periods {
        if periods % repeat != 0 { continue; }
        let mut separated = false;
        for (row, &tolerance) in candidate.coeffs.iter().zip(tolerances) {
            if repeat_rms(row, periods, repeat)? > tolerance {
                separated = true;
                break;
            }
        }
        if !separated {
            return Err(bad("candidate is unresolved from a shorter forcing-synchronized cycle"));
        }
    }
    Ok(chart)
}

// ||x(theta+2*pi*r/m)-x(theta)||_RMS^2
// = 8 sum_{n>0} |X_n|^2 sin^2(pi*n*r/m). Exact integer remainders
// keep commensurate harmonics EXACTLY zero instead of sin(k*pi) roundoff.
// Stable hypotenuse accumulation avoids squaring physical state amplitudes.
fn repeat_rms(row: &[C64], periods: usize, repeat: usize) -> Result<f64, PhsOrbitError> {
    let mut norm = 0.0;
    for (n, coefficient) in row.iter().enumerate().skip(1) {
        let phase = ((n % periods) * repeat) % periods;
        if phase == 0 { continue; }
        let factor = 2.0 * core::f64::consts::SQRT_2
            * det::sin(core::f64::consts::PI * phase as f64 / periods as f64).abs();
        let value = factor * det::hypot(coefficient.re, coefficient.im);
        norm = det::hypot(norm, value);
        if !norm.is_finite() {
            return Err(bad("subharmonic waveform separation is not representable"));
        }
    }
    Ok(norm)
}

#[cfg(test)]
mod tests;
