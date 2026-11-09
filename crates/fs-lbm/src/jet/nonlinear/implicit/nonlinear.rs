//! Source-retaining nonlinear solves for fixed-area IMEX viscosity stages.
//!
//! With M_i=(A_i+A_{i+1})/2, solve M*u-dt*D(3*A*mu(gamma)*Du/rho)=q_source.
//! gamma=sqrt(3)|Du| and d(mu*Du)/d(Du)=mu_t are the SAME original material
//! response. Its positive differential stress gives a symmetric positive
//! definite cyclic Jacobian: M plus positive weighted periodic springs.
//! Each Newton correction uses the native O(N) cyclic LDL/rank-one solver.
//!
//! Only algebraic line-search candidates may backtrack after source refusal.
//! Every accepted iterate and actual RK stage remains inside original support;
//! no rate, viscosity or tangent is clipped. A failed search containing no
//! supported candidate returns the original material diagnostic.

use super::*;

/// Explicit work and accuracy controls for nonlinear implicit viscosity.
/// All bounds apply to each stage; the outer step budget bounds stage count.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetImplicitControls {
    /// Maximum Newton corrections for EACH nonlinear implicit RK stage.
    pub maximum_iterations: usize,
    /// Maximum halvings per correction, in 0..=52.
    pub maximum_backtracks: usize,
    /// Relative velocity-scaled stage residual, strictly between zero and one.
    /// Reduce below the desired temporal error; this is not a time-step control.
    pub relative_tolerance: f64,
}

impl Default for SlenderJetImplicitControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 40,
            maximum_backtracks: 24,
            relative_tolerance: 1e-12,
        }
    }
}

impl SlenderJetImplicitControls {
    pub(super) fn validate(self) -> Result<(), NonlinearJetError> {
        if self.maximum_iterations == 0
            || self.maximum_backtracks > 52
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance <= 0.0
            || self.relative_tolerance >= 1.0
        {
            return Err(NonlinearJetError::InvalidInput(
                "nonlinear viscosity controls require iterations, 0..=52 backtracks and 0<tolerance<1",
            ));
        }
        Ok(())
    }
}

/// Five arrays, allocated only when nonlinear implicit viscosity is selected.
#[derive(Clone, Debug)]
pub(super) struct ConstitutiveSolve {
    velocity: Vec<f64>,
    trial: Vec<f64>,
    source: Vec<f64>,
    flux: Vec<f64>,
    tangent: Vec<f64>,
}

impl ConstitutiveSolve {
    pub(super) fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            velocity: zeros(n, gate)?,
            trial: zeros(n, gate)?,
            source: zeros(n, gate)?,
            flux: zeros(n, gate)?,
            tangent: zeros(n, gate)?,
        })
    }

    /// Return the number of rejected algebraic line-search candidates.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn apply(
        &mut self,
        linear: &mut CyclicSolve,
        state: &mut Fields,
        law: &CarreauSlenderJetParameters,
        dx: f64,
        dt: f64,
        controls: SlenderJetImplicitControls,
        gate: &CancelGate,
    ) -> Result<usize, NonlinearJetError> {
        poll(gate)?;
        let n = state.a.len();
        let mut momentum = Sum::default();
        let mut mass_sum = Sum::default();
        let first_velocity = finite(state.q[0] / average(state.a[0], state.a[1]))?;
        let mut uniform = true;
        for i in 0..n {
            poll_at(gate, i)?;
            positive(state.a[i])?;
            let mass = positive(average(state.a[i], state.a[next(i, n)]))?;
            self.source[i] = finite(state.q[i])?;
            momentum.add(self.source[i]);
            mass_sum.add(mass);
            uniform &= self.source[i] / mass == first_velocity;
            self.velocity[i] = 0.0;
        }
        // Preserve the exact viscous null mode. The enclosing stage evaluator
        // still checks the original material and geometric domain.
        if uniform {
            return Ok(0);
        }
        let mass_sum = positive(mass_sum.total())?;
        let mean = finite(momentum.total() / mass_sum)?;
        let mut scale = 0.0_f64;
        for i in 0..n {
            poll_at(gate, i)?;
            let mass = average(state.a[i], state.a[next(i, n)]);
            scale = scale.hypot((self.source[i] - mass * mean) / mass);
        }
        let scale = positive(scale)?;
        // Start from the source-weighted mean, which has zero strain and is
        // admitted by every Carreau jet. Centering keeps translation out of
        // the nonlinear residual scale and preserves the momentum null mode.
        let mut norm = evaluate(
            &state.a,
            &self.velocity,
            &self.source,
            mean,
            law,
            dx,
            dt,
            &mut self.flux,
            &mut self.tangent,
            &mut linear.rhs,
            gate,
        )?;
        let mut backtracks = 0;
        for iteration in 0..=controls.maximum_iterations {
            let relative = finite(norm / scale)?;
            if relative <= controls.relative_tolerance {
                self.finish(
                    state,
                    linear,
                    law,
                    dx,
                    dt,
                    mean,
                    mass_sum,
                    momentum.total(),
                    scale,
                    controls,
                    iteration,
                    backtracks,
                    gate,
                )?;
                return Ok(backtracks);
            }
            if iteration == controls.maximum_iterations {
                return Err(NonlinearJetError::ImplicitViscosity {
                    iterations: iteration,
                    backtracks,
                    relative_residual: relative,
                });
            }
            // evaluate left -R in rhs and retained the SAME iterate's tangent.
            linear.solve_edges(&state.a, &|i| self.tangent[i], 0.0, gate)?;
            let mut alpha = 1.0;
            let mut accepted = None;
            let mut supported = false;
            let mut refusal = None;
            for attempt in 0..=controls.maximum_backtracks {
                for i in 0..n {
                    poll_at(gate, i)?;
                    self.trial[i] = self.velocity[i] + alpha * linear.y[i];
                }
                match evaluate(
                    &state.a,
                    &self.trial,
                    &self.source,
                    mean,
                    law,
                    dx,
                    dt,
                    &mut self.flux,
                    &mut self.tangent,
                    &mut linear.rhs,
                    gate,
                ) {
                    Ok(trial_norm) => {
                        supported = true;
                        if trial_norm <= (1.0 - 1e-4 * alpha) * norm
                            || trial_norm / scale <= controls.relative_tolerance
                        {
                            accepted = Some(trial_norm);
                            break;
                        }
                    }
                    Err(
                        error @ (NonlinearJetError::Rheology { .. }
                        | NonlinearJetError::NonRepresentable),
                    ) => refusal = Some(error),
                    Err(error) => return Err(error),
                }
                backtracks += 1;
                if attempt < controls.maximum_backtracks {
                    alpha *= 0.5;
                }
            }
            let Some(trial_norm) = accepted else {
                if !supported {
                    if let Some(error) = refusal {
                        return Err(error);
                    }
                }
                return Err(NonlinearJetError::ImplicitViscosity {
                    iterations: iteration + 1,
                    backtracks,
                    relative_residual: relative,
                });
            };
            core::mem::swap(&mut self.velocity, &mut self.trial);
            // Accepted trial stress, residual and tangent are already current.
            norm = trial_norm;
        }
        unreachable!("bounded implicit viscosity solve accepts or refuses")
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &mut self,
        state: &mut Fields,
        linear: &mut CyclicSolve,
        law: &CarreauSlenderJetParameters,
        dx: f64,
        dt: f64,
        mean: f64,
        mass_sum: f64,
        source_momentum: f64,
        scale: f64,
        controls: SlenderJetImplicitControls,
        iterations: usize,
        backtracks: usize,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        let n = state.a.len();
        let mut solved_momentum = Sum::default();
        for i in 0..n {
            poll_at(gate, i)?;
            let mass = average(state.a[i], state.a[next(i, n)]);
            state.q[i] = finite(mass * (self.velocity[i] + mean))?;
            solved_momentum.add(state.q[i]);
        }
        let shift = finite((source_momentum - solved_momentum.total()) / mass_sum)?;
        for i in 0..n {
            poll_at(gate, i)?;
            let mass = average(state.a[i], state.a[next(i, n)]);
            state.q[i] = finite(state.q[i] + mass * shift)?;
            self.trial[i] = finite(state.q[i] / mass)?;
        }
        // Requery the velocities that actually survive q reconstruction and
        // mean restoration. This is an actual stage, so source errors escape.
        evaluate(
            &state.a,
            &self.trial,
            &self.source,
            0.0,
            law,
            dx,
            dt,
            &mut self.flux,
            &mut self.tangent,
            &mut linear.rhs,
            gate,
        )?;
        let mut norm = 0.0_f64;
        let mut roundoff = 0.0_f64;
        for i in 0..n {
            poll_at(gate, i)?;
            let r = next(i, n);
            let mass = average(state.a[i], state.a[r]);
            let defect = finite((state.q[i] - self.source[i]) - (self.flux[r] - self.flux[i]))?;
            norm = norm.hypot(defect / mass);
            roundoff = roundoff.hypot(finite(
                (state.q[i].abs() + self.source[i].abs() + self.flux[i].abs() + self.flux[r].abs())
                    / mass,
            )?);
        }
        // The requested centered residual tolerance plus the same arithmetic
        // backward-error allowance used by the native linear stage solver.
        if norm > controls.relative_tolerance * scale + 512.0 * f64::EPSILON * roundoff {
            return Err(NonlinearJetError::ImplicitViscosity {
                iterations,
                backtracks,
                relative_residual: finite(norm / scale)?,
            });
        }
        poll(gate)
    }
}

/// Return ||R/M||_2 and retain -R and the exact positive Jacobian springs.
#[allow(clippy::too_many_arguments)]
fn evaluate(
    area: &[f64],
    velocity: &[f64],
    source: &[f64],
    mean: f64,
    law: &CarreauSlenderJetParameters,
    dx: f64,
    dt: f64,
    flux: &mut [f64],
    tangent: &mut [f64],
    negative_residual: &mut [f64],
    gate: &CancelGate,
) -> Result<f64, NonlinearJetError> {
    let n = area.len();
    for i in 0..n {
        poll_at(gate, i)?;
        let strain = finite((velocity[i] - velocity[previous(i, n)]) / dx)?;
        let sample = law.sample(strain, Some(i))?;
        let nu = positive(sample.dynamic_viscosity_pa_s() / law.density_kg_m3)?;
        let differential = positive(sample.stress_tangent_pa_s() / law.density_kg_m3)?;
        flux[i] = finite((3.0 * nu * area[i] * strain) * (dt / dx))?;
        tangent[i] = positive(((3.0 * differential * dt) / dx / dx) * area[i])?;
    }
    let mut norm = 0.0_f64;
    for i in 0..n {
        poll_at(gate, i)?;
        let r = next(i, n);
        let mass = average(area[i], area[r]);
        let residual =
            finite((mass * velocity[i] - (source[i] - mass * mean)) - (flux[r] - flux[i]))?;
        negative_residual[i] = -residual;
        norm = norm.hypot(residual / mass);
    }
    finite(norm)
}

#[cfg(test)]
mod tests;
