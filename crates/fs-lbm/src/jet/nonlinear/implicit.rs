//! Optional third-order IMEX viscosity integration for the existing jet owner.
//!
//! Pareschi & Russo, "Implicit-explicit Runge-Kutta schemes and applications
//! to hyperbolic systems with relaxation", Table 6:
//! <https://arxiv.org/pdf/1009.2757>.
//! The IMEX-SSP3(4,3,3) explicit tableau is SSPRK3; the implicit tableau is
//! L-stable. Four cyclic tridiagonal solves remove the explicit viscosity
//! timestep restriction. Capillary/advection/drainage limits remain active.
//!
//! L-stability is not a claim of uniformly accurate unresolved viscous
//! transients: the scheme is not stiffly accurate. Its embedded second-order
//! difference is a diagnostic, not a global error certificate. Authored
//! maximum_step_s remains necessary to resolve dynamics of interest.
//!
//! Each stage has a known area field. Newtonian viscosity gives a linear
//! face-velocity solve; Carreau viscosity uses bounded damped Newton with its
//! original positive differential stress. Both use the same cyclic positive
//! definite Jacobian: a tridiagonal path plus one positive rank-one edge.
//! Direct LDL and Sherman-Morrison work is O(N) per correction, with bounded
//! cancellation polls. Nonlinear residual error must resolve temporal error.
//! The constant-velocity nullspace is restored to roundoff after the solve;
//! a componentwise residual check refuses a numerically failed solve.
//!
//! Final positive RK weights integrate physical and upwind powers separately.
//! The existing signed energy remainder includes time/space truncation; no
//! implicit numerical damping is relabelled as molecular heat.

use super::*;

mod nonlinear;
pub(super) use nonlinear::ConstitutiveSolve;
pub use nonlinear::SlenderJetImplicitControls;

/// Time integration policy; changing it does not change any physical coefficient.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SlenderJetIntegrationMethod {
    #[default]
    ExplicitSspRk3,
    /// Third-order IMEX-SSP3(4,3,3), with converged implicit viscous stages.
    ImplicitViscosity,
}

/// Largest accepted-step difference between the third-order state and an
/// embedded second-order estimate. These are local indicators in SI units,
/// not accumulated/global bounds and not an automatic tolerance controller.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SlenderJetTemporalError {
    pub maximum_radius_error_m: f64,
    pub maximum_velocity_error_m_s: f64,
}

const ALPHA: f64 = 0.241_694_260_788_21;
// Preserve the algebraic order identities at the stated coefficient precision.
const BETA: f64 = ALPHA / 4.0;
const ETA: f64 = (1.0 - 2.0 * ALPHA) / 4.0;
const GAMMA: f64 = 0.5 - BETA - ETA - ALPHA;

#[derive(Clone, Debug)]
pub(super) struct Workspace {
    f2: Fields,
    f3: Fields,
    g1: Vec<f64>,
    g2: Vec<f64>,
    g3: Vec<f64>,
    g4: Vec<f64>,
    solve: CyclicSolve,
    nonlinear: Option<ConstitutiveSolve>,
}
impl Workspace {
    fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            f2: Fields::zeros(n, gate)?,
            f3: Fields::zeros(n, gate)?,
            g1: zeros(n, gate)?,
            g2: zeros(n, gate)?,
            g3: zeros(n, gate)?,
            g4: zeros(n, gate)?,
            solve: CyclicSolve {
                diagonal: zeros(n, gate)?,
                y: zeros(n, gate)?,
                z: zeros(n, gate)?,
                rhs: zeros(n, gate)?,
            },
            nonlinear: None,
        })
    }

    pub(super) fn prepare_rheology(
        &mut self,
        law: &CarreauProfile,
        n: usize,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        if !law.is_uniform_newtonian() && self.nonlinear.is_none() {
            let nonlinear = ConstitutiveSolve::new(n, gate)?;
            poll(gate)?;
            self.nonlinear = Some(nonlinear);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn apply(
        &mut self,
        state: &mut Fields,
        p: SlenderJetParameters,
        carreau: Option<&CarreauProfile>,
        dx: f64,
        dt: f64,
        controls: SlenderJetImplicitControls,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        if let Some(law) = carreau.filter(|law| !law.is_uniform_newtonian()) {
            self.nonlinear
                .as_mut()
                .ok_or(NonlinearJetError::InvalidState(
                    "nonlinear viscosity workspace is absent",
                ))?
                .apply(&mut self.solve, state, law, dx, dt, controls, gate)
                .map(|_| ())
        } else {
            // Exact degenerate source laws share the original linear arithmetic.
            // The enclosing evaluator still queries their original rate support.
            self.solve
                .apply(state, p.kinematic_viscosity_m2_s, dx, dt, gate)
        }
    }
}

impl PeriodicSlenderJet {
    #[must_use]
    pub const fn integration_method(&self) -> SlenderJetIntegrationMethod {
        self.integration_method
    }

    #[must_use]
    pub const fn implicit_viscosity_controls(&self) -> SlenderJetImplicitControls {
        self.implicit_controls
    }

    /// Bound the original-law nonlinear work in each implicit RK stage.
    /// Newtonian/degenerate linear stages retain their native direct solve.
    /// This changes no accepted field, source, geometry, time or policy.
    pub fn set_implicit_viscosity_controls(
        &mut self,
        controls: SlenderJetImplicitControls,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        poll(gate)?;
        controls.validate()?;
        poll(gate)?;
        self.implicit_controls = controls;
        Ok(())
    }

    /// Select a numerical policy without modifying accepted state, source
    /// coefficients, geometry, time or diagnostics. Implicit storage is
    /// allocated once, retained across calls and reused after policy switches.
    /// Cancellation/allocation failure leaves the policy and state unchanged.
    pub fn set_integration_method(
        &mut self,
        method: SlenderJetIntegrationMethod,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        poll(gate)?;
        if method == SlenderJetIntegrationMethod::ImplicitViscosity
            && self.implicit_workspace.is_none()
        {
            let mut workspace = Workspace::new(self.accepted.a.len(), gate)?;
            if self
                .carreau
                .as_ref()
                .is_some_and(|law| !law.is_uniform_newtonian())
            {
                workspace.nonlinear = Some(ConstitutiveSolve::new(self.accepted.a.len(), gate)?);
            }
            poll(gate)?;
            self.implicit_workspace = Some(workspace);
        }
        poll(gate)?;
        self.integration_method = method;
        Ok(())
    }
}

type Attempt = Option<(Evaluation, f64, f64, Option<SlenderJetTemporalError>)>;

/// Reuse the parent transaction driver; only this substep integrator differs.
pub(super) fn step(
    base: &Fields,
    stage: &mut Fields,
    rhs: &mut Fields,
    spatial: &mut Spatial,
    work: &mut Workspace,
    p: SlenderJetParameters,
    carreau: Option<&CarreauProfile>,
    limits: SlenderJetLimits,
    dx: f64,
    dt: f64,
    cfl: f64,
    controls: SlenderJetImplicitControls,
    gate: &CancelGate,
) -> Result<Attempt, NonlinearJetError> {
    let n = base.a.len();
    // Explicit rows: [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,1/4,1/4,0].
    // Implicit rows: [a,0,0,0], [-a,a,0,0], [0,1-a,a,0],
    // [b,eta,1/2-b-eta-a,a]. Both final weights are [0,1/6,1/6,2/3].
    copy_fields(stage, base, gate)?;
    work.apply(stage, p, carreau, dx, ALPHA * dt, controls, gate)?;
    let first = evaluate_constitutive(stage, p, carreau, limits, dx, spatial, rhs, false, gate)?;
    if dt > stable_step(first.rate_per_s, cfl) {
        return Ok(None);
    }
    viscous_rhs(stage, spatial, p, carreau, dx, &mut work.g1, gate)?;

    for i in 0..n {
        poll_at(gate, i)?;
        stage.a[i] = base.a[i];
        stage.q[i] = base.q[i] - (ALPHA * dt) * work.g1[i];
    }
    work.apply(stage, p, carreau, dx, ALPHA * dt, controls, gate)?;
    let second = evaluate_constitutive(stage, p, carreau, limits, dx, spatial, rhs, false, gate)?;
    if dt > stable_step(second.rate_per_s, cfl) {
        return Ok(None);
    }
    spatial.start_viscous_heat(dt / 6.0, gate)?;
    copy_fields(&mut work.f2, rhs, gate)?;
    viscous_rhs(stage, spatial, p, carreau, dx, &mut work.g2, gate)?;

    for i in 0..n {
        poll_at(gate, i)?;
        stage.a[i] = base.a[i] + dt * work.f2.a[i];
        stage.q[i] = base.q[i] + dt * (work.f2.q[i] + (1.0 - ALPHA) * work.g2[i]);
    }
    work.apply(stage, p, carreau, dx, ALPHA * dt, controls, gate)?;
    let third = evaluate_constitutive(stage, p, carreau, limits, dx, spatial, rhs, false, gate)?;
    if dt > stable_step(third.rate_per_s, cfl) {
        return Ok(None);
    }
    spatial.add_viscous_heat(dt / 6.0, gate)?;
    copy_fields(&mut work.f3, rhs, gate)?;
    viscous_rhs(stage, spatial, p, carreau, dx, &mut work.g3, gate)?;

    for i in 0..n {
        poll_at(gate, i)?;
        stage.a[i] = base.a[i] + (0.25 * dt) * (work.f2.a[i] + work.f3.a[i]);
        stage.q[i] = base.q[i]
            + dt * (0.25 * (work.f2.q[i] + work.f3.q[i])
                + BETA * work.g1[i]
                + ETA * work.g2[i]
                + GAMMA * work.g3[i]);
    }
    work.apply(stage, p, carreau, dx, ALPHA * dt, controls, gate)?;
    let fourth = evaluate_constitutive(stage, p, carreau, limits, dx, spatial, rhs, false, gate)?;
    if dt > stable_step(fourth.rate_per_s, cfl) {
        return Ok(None);
    }
    spatial.add_viscous_heat((2.0 / 3.0) * dt, gate)?;
    viscous_rhs(stage, spatial, p, carreau, dx, &mut work.g4, gate)?;

    // Publish the RK weighted state, NOT the last implicit stage.
    for i in 0..n {
        poll_at(gate, i)?;
        stage.a[i] =
            base.a[i] + dt * ((work.f2.a[i] + work.f3.a[i]) / 6.0 + (2.0 / 3.0) * rhs.a[i]);
        stage.q[i] = base.q[i]
            + dt * ((work.f2.q[i] + work.g2[i] + work.f3.q[i] + work.g3[i]) / 6.0
                + (2.0 / 3.0) * (rhs.q[i] + work.g4[i]));
    }
    let end = evaluate_constitutive(stage, p, carreau, limits, dx, spatial, rhs, false, gate)?;
    if dt > stable_step(end.rate_per_s, cfl) {
        return Ok(None);
    }
    // b_hat=[0,1/2,1/2,0] satisfies both additive second-order conditions.
    // No new solve is needed. An inadmissible embedded area requests a retry.
    let mut estimate = SlenderJetTemporalError::default();
    for i in 0..n {
        poll_at(gate, i)?;
        let j = next(i, n);
        let area_hat = finite(base.a[i] + (0.5 * dt) * (work.f2.a[i] + work.f3.a[i]))?;
        let area_hat_next = finite(base.a[j] + (0.5 * dt) * (work.f2.a[j] + work.f3.a[j]))?;
        if area_hat <= 0.0 || area_hat_next <= 0.0 {
            return Ok(None);
        }
        let q_hat = finite(
            base.q[i] + (0.5 * dt) * (work.f2.q[i] + work.g2[i] + work.f3.q[i] + work.g3[i]),
        )?;
        let u_hat = finite(q_hat / average(area_hat, area_hat_next))?;
        estimate.maximum_radius_error_m = estimate
            .maximum_radius_error_m
            .max((spatial.h[i] - area_hat.sqrt()).abs());
        estimate.maximum_velocity_error_m_s = estimate
            .maximum_velocity_error_m_s
            .max(finite((spatial.u[i] - u_hat).abs())?);
    }
    let physical = finite(
        dt * (second.viscous_power_w / 6.0
            + third.viscous_power_w / 6.0
            + (2.0 / 3.0) * fourth.viscous_power_w),
    )?;
    let numerical = finite(
        dt * (second.upwind_power_w / 6.0
            + third.upwind_power_w / 6.0
            + (2.0 / 3.0) * fourth.upwind_power_w),
    )?;
    Ok(Some((end, physical, numerical, Some(estimate))))
}

fn copy_fields(to: &mut Fields, from: &Fields, gate: &CancelGate) -> Result<(), NonlinearJetError> {
    for i in 0..to.a.len() {
        poll_at(gate, i)?;
        to.a[i] = from.a[i];
        to.q[i] = from.q[i];
    }
    Ok(())
}

/// Same centered conservative constitutive stress used by the parent operator.
pub(super) fn viscous_rhs(
    state: &Fields,
    spatial: &Spatial,
    p: SlenderJetParameters,
    carreau: Option<&CarreauProfile>,
    dx: f64,
    result: &mut [f64],
    gate: &CancelGate,
) -> Result<(), NonlinearJetError> {
    let n = state.a.len();
    if let Some(law) = carreau.filter(|law| !law.is_uniform_newtonian()) {
        // Reuse the output for one source sample per cell, then differentiate
        // those SAME stresses in place. Retain the periodic wrap before writes.
        for i in 0..n {
            poll_at(gate, i)?;
            let strain = (spatial.u[i] - spatial.u[previous(i, n)]) / dx;
            let [viscosity, _] = law.coefficients(strain, i)?;
            let nu = positive(viscosity / p.density_kg_m3)?;
            result[i] = finite(3.0 * nu * state.a[i] * strain)?;
        }
        let first = result[0];
        for i in 0..n {
            poll_at(gate, i)?;
            let right = if i + 1 == n { first } else { result[i + 1] };
            result[i] = finite((right - result[i]) / dx)?;
        }
        return Ok(());
    }
    let nu = p.kinematic_viscosity_m2_s;
    for i in 0..n {
        poll_at(gate, i)?;
        let l = previous(i, n);
        let r = next(i, n);
        let left = 3.0 * nu * state.a[i] * ((spatial.u[i] - spatial.u[l]) / dx);
        let right = 3.0 * nu * state.a[r] * ((spatial.u[r] - spatial.u[i]) / dx);
        result[i] = finite((right - left) / dx)?;
    }
    Ok(())
}

#[derive(Clone, Debug)]
pub(super) struct CyclicSolve {
    diagonal: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    rhs: Vec<f64>,
}
impl CyclicSolve {
    #[cfg(feature = "implicit-enthalpy-coupling")]
    pub(super) fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            diagonal: zeros(n, gate)?,
            y: zeros(n, gate)?,
            z: zeros(n, gate)?,
            rhs: zeros(n, gate)?,
        })
    }

    /// Frozen positive node weights and periodic springs, reused by the
    /// thermoviscous block preconditioner. This is not a joint SPD claim.
    #[cfg(feature = "implicit-enthalpy-coupling")]
    pub(super) fn solve_weighted(
        &mut self,
        mass: &impl Fn(usize) -> f64,
        coupling: &impl Fn(usize) -> f64,
        rhs: &[f64],
        output: &mut [f64],
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        let n = self.rhs.len();
        if rhs.len() != n || output.len() != n {
            return Err(NonlinearJetError::InvalidInput("cyclic block dimensions"));
        }
        for i in 0..n {
            poll_at(gate, i)?;
            self.rhs[i] = finite(rhs[i])?;
        }
        self.solve_weights(n, mass, coupling, 0.0, gate)?;
        for i in 0..n {
            poll_at(gate, i)?;
            output[i] = self.y[i];
        }
        Ok(())
    }

    /// Solve q - dt*3nu*D(A*D(q/M)) = rhs with M_i=(A_i+A_i+1)/2.
    /// Areas remain fixed during this linear stage.
    pub(super) fn apply(
        &mut self,
        state: &mut Fields,
        nu: f64,
        dx: f64,
        dt: f64,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        poll(gate)?;
        if nu == 0.0 {
            return Ok(());
        }
        let n = state.a.len();
        let coefficient = positive((3.0 * nu * dt) / dx / dx)?;
        let mut total_rhs = Sum::default();
        let mut total_mass = Sum::default();
        let mut uniform = true;
        let first_velocity = finite(state.q[0] / average(state.a[0], state.a[1]))?;
        for i in 0..n {
            poll_at(gate, i)?;
            let mass = positive(average(state.a[i], state.a[next(i, n)]))?;
            if state.a[i] <= 0.0 || !state.a[i].is_finite() {
                return Err(NonlinearJetError::InvalidState(
                    "nonpositive implicit-stage area",
                ));
            }
            let source = finite(state.q[i])?;
            self.rhs[i] = source;
            total_rhs.add(source);
            total_mass.add(mass);
            uniform &= source / mass == first_velocity;
        }
        // A constant velocity is the exact viscous null mode, including on
        // a nonuniform area field; do not introduce solver roundoff into it.
        if uniform {
            return Ok(());
        }
        let mass_sum = positive(total_mass.total())?;
        let mean = finite(total_rhs.total() / mass_sum)?;
        self.solve_edges(&state.a, &|i| coefficient * state.a[i], mean, gate)?;
        let mut solved_momentum = Sum::default();
        for i in 0..n {
            poll_at(gate, i)?;
            let u = self.y[i];
            let mass = average(state.a[i], state.a[next(i, n)]);
            state.q[i] = finite(mass * u)?;
            solved_momentum.add(state.q[i]);
        }
        // Roundoff projection along the constant-velocity nullspace. It does
        // not change viscous gradients in exact arithmetic and enforces the
        // periodic momentum identity; the residual is checked afterwards.
        let shift = finite((total_rhs.total() - solved_momentum.total()) / mass_sum)?;
        for i in 0..n {
            poll_at(gate, i)?;
            let mass = average(state.a[i], state.a[next(i, n)]);
            state.q[i] = finite(state.q[i] + mass * shift)?;
            self.y[i] = finite(state.q[i] / mass)?;
        }
        for i in 0..n {
            poll_at(gate, i)?;
            let l = previous(i, n);
            let r = next(i, n);
            let left = finite(coefficient * state.a[i] * (self.y[i] - self.y[l]))?;
            let right = finite(coefficient * state.a[r] * (self.y[r] - self.y[i]))?;
            let residual = finite((state.q[i] - self.rhs[i]) - (right - left))?;
            let scale = finite(state.q[i].abs() + self.rhs[i].abs() + left.abs() + right.abs())?;
            if residual.abs() > 512.0 * f64::EPSILON * scale.max(f64::MIN_POSITIVE) {
                return Err(NonlinearJetError::InvalidState(
                    "implicit viscosity solve residual",
                ));
            }
        }
        poll(gate)
    }

    /// Solve (M + periodic weighted springs)*u = rhs. The same operation
    /// order is retained for constant-viscosity stages; nonlinear corrections
    /// supply the original differential stress at each edge.
    fn solve_edges(
        &mut self,
        area: &[f64],
        coupling: &impl Fn(usize) -> f64,
        mean: f64,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        self.solve_weights(
            area.len(),
            &|i| average(area[i], area[next(i, area.len())]),
            coupling,
            mean,
            gate,
        )
    }

    fn solve_weights(
        &mut self,
        n: usize,
        mass: &impl Fn(usize) -> f64,
        coupling: &impl Fn(usize) -> f64,
        mean: f64,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        let edge = finite(coupling(0))?;
        for i in 0..n {
            poll_at(gate, i)?;
            let mass = positive(mass(i))?;
            // Remove only the periodic spring: B=T+edge*v*v^T, v=e0-eN.
            let left = if i == 0 { 0.0 } else { coupling(i) };
            let right = if i + 1 == n { 0.0 } else { coupling(i + 1) };
            let diagonal = finite(mass + left + right)?;
            let source = finite(self.rhs[i] - mass * mean)?;
            let basis = if i == 0 {
                1.0
            } else if i + 1 == n {
                -1.0
            } else {
                0.0
            };
            if i == 0 {
                self.diagonal[i] = positive(diagonal)?;
                self.y[i] = source;
                self.z[i] = basis;
            } else {
                let off = -left;
                let factor = finite(off / self.diagonal[i - 1])?;
                self.diagonal[i] = positive(diagonal - factor * off)?;
                self.y[i] = finite(source - factor * self.y[i - 1])?;
                self.z[i] = finite(basis - factor * self.z[i - 1])?;
            }
        }
        for i in (0..n).rev() {
            poll_at(gate, i)?;
            let off = if i + 1 == n { 0.0 } else { -coupling(i + 1) };
            let next_y = if i + 1 == n { 0.0 } else { self.y[i + 1] };
            let next_z = if i + 1 == n { 0.0 } else { self.z[i + 1] };
            self.y[i] = finite((self.y[i] - off * next_y) / self.diagonal[i])?;
            self.z[i] = finite((self.z[i] - off * next_z) / self.diagonal[i])?;
        }
        let denominator = positive(1.0 + edge * (self.z[0] - self.z[n - 1]))?;
        let factor = finite(edge * (self.y[0] - self.y[n - 1]) / denominator)?;
        for i in 0..n {
            poll_at(gate, i)?;
            self.y[i] = finite(self.y[i] - factor * self.z[i] + mean)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
