//! Bounded simultaneous contact for ONE hammer, not another time integrator.
//!
//! The caller has already condensed strings, board and hammer into a reciprocal
//! displacement compliance, and adds each existing Prony compliance diagonally.
//! We solve e + A*f(e) = free using the same discrete felt force as the scalar
//! path. Separate local overlaps and histories are never averaged together.
use fs_la::LuWorkspace;
use fs_material::WoolFelt;
use super::{average, State};

/// Four sites on each of three unison strings. No contacts may be dropped.
pub const MAX_SITES: usize = 12;
const MAX_UPDATES: usize = 32;
const MAX_BACKTRACKS: usize = 24;
const POSITION_TOL: f64 = 1e-13; // same displacement tolerance as felt::solve
const FORCE_ABS: f64 = 1e-5; // original engine's independent force-equation gate
const FORCE_REL: f64 = 1e-8;

#[derive(Clone, Copy)]
pub struct Site<'a> {
    pub law: &'a WoolFelt,
    pub history: &'a State,
    pub start_m: f64,
    pub thickness_m: f64,
    pub area_m2: f64,
}
impl Site<'_> {
    fn force(self, end_m: f64) -> (f64, f64) {
        average(self.law, self.history, self.start_m, end_m, self.thickness_m, self.area_m2)
    }
    fn maximum(self) -> f64 { self.thickness_m * self.law.eps_densify }
}

/// Cold-allocated LU; all other scratch is fixed-size. Failure leaves caller
/// forces unchanged. Every call reinitializes scratch, so retry is independent
/// of failed iterates. Physical/material history remains read-only throughout.
pub struct Workspace {
    lu: LuWorkspace,
    end: [f64; MAX_SITES],
    trial: [f64; MAX_SITES],
    force: [f64; MAX_SITES],
    derivative: [f64; MAX_SITES],
    residual: [f64; MAX_SITES],
    delta: [f64; MAX_SITES],
    jacobian: [f64; MAX_SITES * MAX_SITES],
}
impl Workspace {
    pub fn new() -> Result<Self, String> {
        Ok(Self { lu: LuWorkspace::new(MAX_SITES).map_err(|e| e.to_string())?,
            end: [0.; MAX_SITES], trial: [0.; MAX_SITES], force: [0.; MAX_SITES],
            derivative: [0.; MAX_SITES], residual: [0.; MAX_SITES],
            delta: [0.; MAX_SITES], jacobian: [0.; MAX_SITES * MAX_SITES] })
    }

    /// A is the already-admitted, row-major n*n displacement compliance [m/N].
    /// `free` has the free Prony deformation removed. `initial_force` is only a
    /// numerical warm start; it contributes no impulse or physical work.
    /// Returns the number of Newton updates; no allocation occurs in this call.
    ///
    /// Newton acts on overlaps, not forces: separated sites can become active
    /// without forcing a negative-force trial or changing an active set. The
    /// exact force derivative differentiates the discrete average, not endpoint
    /// stress. Line search limits trial strain to the original densification
    /// domain. Neither the accepted overlap nor physical force is clamped.
    pub fn solve(&mut self, sites: &[Site<'_>], a: &[f64], free: &[f64],
        initial_force: &[f64], output: &mut [f64]) -> Result<usize, &'static str>
    {
        let n = sites.len();
        if n == 0 || n > MAX_SITES || a.len() != n*n || free.len() != n
            || initial_force.len() != n || output.len() != n
            || a.iter().chain(free).chain(initial_force).any(|v| !v.is_finite())
            || initial_force.iter().any(|f| *f < 0.)
            || (0..n).any(|i| a[i*n+i] <= 0.)
            || sites.iter().any(|s| !s.start_m.is_finite()
                || !s.thickness_m.is_finite() || s.thickness_m <= 0.
                || !s.area_m2.is_finite() || s.area_m2 <= 0.
                || !s.maximum().is_finite() || s.maximum() <= 0.) {
            return Err("invalid simultaneous hammer contact dimensions or physical inputs");
        }
        for i in 0..n {
            let end = free[i] - (0..n).map(|j| a[i*n+j]*initial_force[j]).sum::<f64>();
            if !end.is_finite() { return Err("hammer contact warm-start overflow"); }
            // Only the starting Newton guess is placed inside the trial domain.
            self.end[i] = end.min(sites[i].maximum());
        }
        for iteration in 0..=MAX_UPDATES {
            let rnorm = evaluate(sites, a, free, &self.end, &mut self.force,
                &mut self.derivative, &mut self.residual)?;
            let mut accepted = rnorm <= POSITION_TOL;
            for i in 0..n {
                // Check the force at the overlap that mechanics WILL produce,
                // not just at Newton's independent overlap variable.
                let actual = free[i] - (0..n).map(|j| a[i*n+j]*self.force[j]).sum::<f64>();
                if !actual.is_finite() || actual > sites[i].maximum() { accepted = false; continue; }
                let expected = sites[i].force(actual).0;
                if !expected.is_finite() || expected < 0.
                    || (self.force[i]-expected).abs() > FORCE_ABS + FORCE_REL*expected.abs() {
                    accepted = false;
                }
            }
            if accepted { output.copy_from_slice(&self.force[..n]); return Ok(iteration); }
            if iteration == MAX_UPDATES { break; }
            self.jacobian.fill(0.);
            // Pad the small block with decoupled identities so one prepared LU
            // handles 1..12 sites, including whole-string una-corda exclusions.
            for i in 0..MAX_SITES { self.jacobian[i*MAX_SITES+i] = 1.; }
            for i in 0..n { for j in 0..n {
                self.jacobian[i*MAX_SITES+j] += a[i*n+j]*self.derivative[j];
            }}
            self.lu.solve_into(&self.jacobian, &self.residual, &mut self.delta)
                .map_err(|_| "simultaneous hammer contact tangent solve refused")?;
            let mut fraction = 1.; let mut decreased = false;
            // Scratch-only trial evaluations: none can update felt or Prony.
            for _ in 0..MAX_BACKTRACKS {
                for i in 0..n { self.trial[i] = self.end[i] - fraction*self.delta[i]; }
                let admissible = (0..n).all(|i| self.trial[i].is_finite()
                    && self.trial[i] <= sites[i].maximum());
                if admissible {
                    let candidate = evaluate(sites, a, free, &self.trial, &mut self.force,
                        &mut self.derivative, &mut self.residual)?;
                    if candidate <= (1. - 1e-4*fraction)*rnorm || candidate <= POSITION_TOL {
                        self.end[..n].copy_from_slice(&self.trial[..n]); decreased = true; break;
                    }
                }
                fraction *= 0.5;
            }
            if !decreased { return Err("simultaneous hammer contact exhausted its bounded line search"); }
        }
        Err("simultaneous hammer contact exhausted 32 Newton updates")
    }
}

fn evaluate(sites: &[Site<'_>], a: &[f64], free: &[f64], end: &[f64; MAX_SITES],
    force: &mut [f64; MAX_SITES], derivative: &mut [f64; MAX_SITES],
    residual: &mut [f64; MAX_SITES]) -> Result<f64, &'static str>
{
    let n = sites.len(); residual.fill(0.);
    for i in 0..n {
        (force[i], derivative[i]) = sites[i].force(end[i]);
        if !force[i].is_finite() || force[i] < 0. || !derivative[i].is_finite() {
            return Err("simultaneous hammer contact has nonfinite or tensile felt force");
        }
    }
    let mut norm = 0.0_f64;
    for i in 0..n {
        residual[i] = end[i] + (0..n).map(|j| a[i*n+j]*force[j]).sum::<f64>() - free[i];
        if !residual[i].is_finite() { return Err("simultaneous hammer contact residual overflow"); }
        norm = norm.max(residual[i].abs());
    }
    Ok(norm)
}

#[cfg(test)]
#[path = "hammer_contact_tests.rs"]
mod tests;
