//! Opt-in second-order Runge--Kutta--Legendre molecular super-time-stepping.
//!
//! The recurrence and negative-real stability interval follow Meyer, Balsara
//! and Aslam (JCP 2014, doi:10.1016/j.jcp.2013.08.021), also written explicitly
//! in Vaidya et al. (MNRAS 2017, doi:10.1093/mnras/stx2176), section 3.3 and
//! appendix A. Only time integration changes: every stage uses the original
//! finite-volume Maxwell--Stefan/enthalpy/heat/stress face owner.

use super::super::evolve::totals;
use super::super::{Mesh1d, MixtureEuler1d, MixtureGas, Totals, allocated, sum};
use super::entropy::EntropyTracker;
use super::evolve::{forward, spatial_with_wall_motion};
use super::reservoir::ReservoirEnds;
use super::{
    DissipationControl, DissipationError, DissipationOperator, DissipationReport, MAX_COMPONENTS,
    TransportLaw, poll, positive,
};
use fs_qty::Time;

/// Shared second-order Legendre recurrence for original molecular operators.
///
/// The numerical coefficients are independent of dimension and geometry.
/// Spatial divergence, physical/source admission and boundary quadrature remain
/// with each original owner. A negative-real linear stability interval does not
/// establish positivity, nonlinear stability or temporal accuracy for a gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rkl2Scheme {
    stages: usize,
}

impl Rkl2Scheme {
    /// Admit the same bounded stage count used by the original axial owner.
    pub fn new(stages: usize) -> Result<Self, DissipationError> {
        if !(2..=256).contains(&stages) {
            return Err(DissipationError::Invalid(
                "RKL2 stage count must be in [2,256]",
            ));
        }
        Ok(Self { stages })
    }

    pub const fn stage_count(self) -> usize {
        self.stages
    }

    /// Extension relative to the original forward-Euler parabolic interval.
    pub fn stability_extent(self) -> f64 {
        extent(self.stages)
    }

    /// Initial stage is `Y1 = Y0 + first_weight * dt * L(Y0)`.
    pub fn first_weight(self) -> f64 {
        b(1) * (1.0 / self.stability_extent())
    }

    /// Physical time fraction of state `Y_stage` in the native recurrence.
    ///
    /// The operator forming `Y_j` observes `Y_(j-1)` at `abscissa(j-1)`;
    /// its reused anchor observes `Y_0` at zero. These fractions solve
    /// `c_j = mu*c_(j-1) + nu*c_(j-2) + mu_tilde + gamma_tilde`.
    /// In particular `c_1` is the first forward weight, not `1 / stages`.
    /// The final admitted state is exactly at the macrostep endpoint.
    pub fn abscissa(self, stage: usize) -> Result<f64, DissipationError> {
        if stage > self.stages {
            return Err(DissipationError::Invalid(
                "RKL2 time stage outside prepared scheme",
            ));
        }
        Ok(match stage {
            0 => 0.0,
            1 => self.first_weight(),
            _ => {
                let j = stage as f64;
                (j * j + j - 2.0) / (4.0 * self.stability_extent())
            }
        })
    }

    /// Weights of the actual RHS evaluations L(Y_0), ..., L(Y_(s-1))
    /// in the complete macrostep boundary integral, divided by dt.
    ///
    /// Reverse accumulation of the native recurrence gives these weights in
    /// O(s) work and caller-owned storage. They include the negative anchor
    /// contribution (for s=2: [-1/2, 3/2]); absolute throughput must take each
    /// weighted flux's magnitude before summing. No weight normalization or
    /// replacement by positive SSP weights is applied.
    pub fn rhs_weights(self, out: &mut [f64]) -> Result<(), DissipationError> {
        if out.len() != self.stages {
            return Err(DissipationError::Invalid("RKL2 RHS weight dimension"));
        }
        out.fill(0.0);
        let mut adjoint = 1.0;
        let mut previous = 0.0;
        let mut anchor = 0.0;
        let mut correction = 0.0;
        for stage in (2..=self.stages).rev() {
            let [mu, nu, mu_tilde, gamma_tilde] = self.coefficients(stage)?;
            out[stage - 1] = adjoint * mu_tilde;
            accumulate(&mut anchor, &mut correction, adjoint * gamma_tilde);
            let older = adjoint * nu;
            adjoint = previous + adjoint * mu;
            previous = older;
        }
        accumulate(&mut anchor, &mut correction, adjoint * self.first_weight());
        out[0] = anchor;
        if !out.iter().all(|weight| weight.is_finite()) {
            return Err(DissipationError::Invalid("RKL2 RHS weight overflow"));
        }
        Ok(())
    }

    /// `[mu, nu, mu_tilde, gamma_tilde]` for each stage `j` in `2..=s`:
    /// `Yj = Y0 + mu*(Yj-1-Y0) + nu*(Yj-2-Y0)`
    /// `     + dt*(mu_tilde*L(Yj-1) + gamma_tilde*L(Y0))`.
    /// Every extensive boundary transfer uses these same signed weights.
    pub fn coefficients(self, stage: usize) -> Result<[f64; 4], DissipationError> {
        if !(2..=self.stages).contains(&stage) {
            return Err(DissipationError::Invalid(
                "RKL2 recurrence stage outside prepared scheme",
            ));
        }
        let j = stage as f64;
        let mu = ((2.0 * j - 1.0) / j) * b(stage) / b(stage - 1);
        let nu = -((j - 1.0) / j) * b(stage) / b(stage - 2);
        let mu_tilde = mu * (1.0 / self.stability_extent());
        let gamma_tilde = -(1.0 - b(stage - 1)) * mu_tilde;
        Ok([mu, nu, mu_tilde, gamma_tilde])
    }
}

/// Bounded RKL2 macrosteps on the existing stationary molecular operator.
///
/// `transport.maximum_step_s` caps the physical macrostep for accuracy, not
/// each internal stage. `transport.parabolic_cfl` retains the original
/// forward-step estimate; its RKL2 extension is `(s*s+s-2)/4` for `s` stages.
/// This is an Estimated engineering lane. The linear diffusion stability
/// interval does not prove nonlinear stability, positivity or entropy decay
/// for the coupled gas system. Every stage must pass the original caloric,
/// species, source-domain and constitutive checks; refused trials are halved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedDissipationControl {
    /// Original residual, macrostep, minimum-step, attempt and retry budgets.
    pub transport: DissipationControl,
    /// Largest number of RKL2 stages per attempted macrostep, in [2,256].
    pub maximum_stages: usize,
    /// Total RHS calls, including rejected trials and endpoint admission.
    /// In [1,1,000,000]; each accepted s-stage step needs s+1 calls.
    pub maximum_operator_evaluations: usize,
}
impl Default for StabilizedDissipationControl {
    fn default() -> Self {
        Self {
            transport: DissipationControl::default(),
            maximum_stages: 32,
            maximum_operator_evaluations: 100_000,
        }
    }
}
impl StabilizedDissipationControl {
    /// Check work and physical time-window controls before touching scratch.
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate()?;
        if !(2..=256).contains(&self.maximum_stages)
            || !(1..=1_000_000).contains(&self.maximum_operator_evaluations)
        {
            return Err(DissipationError::Invalid("stabilized molecular controls"));
        }
        Ok(())
    }
}

/// Complete-duration physical balances and actual stabilization work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedDissipationReport {
    /// Original mass, momentum, species-referenced and absolute-energy ledgers.
    /// Accepted steps are RKL2 macrosteps; the parabolic number is unscaled.
    pub transport: DissipationReport,
    /// Actual RHS calls, including rejected trials and endpoint admission.
    pub operator_evaluations: usize,
    /// Largest stage count among accepted macrosteps; zero for zero duration.
    pub maximum_stages_used: usize,
    /// Largest accepted `dt * rate / ((s*s+s-2)/4)`, bounded by the base CFL.
    pub maximum_scaled_parabolic_number: f64,
}

impl<L: TransportLaw> DissipationOperator<L> {
    /// Advance a complete stationary molecular interval using bounded RKL2.
    ///
    /// This opt-in path accelerates diffusion-limited intervals while retaining
    /// actual-state transport coefficients, Soret/Dufour coupling, diffusive
    /// enthalpy, viscous heating and conservative boundary transfers. It reuses
    /// the three original state arrays and lazily reserves one extra face
    /// array, with no per-stage allocation. Boundary transfers follow the same
    /// RKL2 recurrence as the conserved field, including its negative weights.
    ///
    /// The caller must choose an accuracy-appropriate maximum macrostep; there
    /// is no temporal error estimator. Rejected trials are never clipped.
    /// Publication is atomic for the entire requested duration, including
    /// late cancellation/work refusal; the flow clock is unchanged. This path
    /// does not integrate mesh motion or hyperbolic/chemical source terms.
    pub fn advance_stabilized_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedDissipationReport, DissipationError> {
        self.advance_stabilized_with_budget(flow, duration, controls, &mut 0, &mut cancelled)
    }

    // The composing owner retains this actual work counter through both halves,
    // coupled retries and complete physical windows. Failed calls consume work.
    pub(in crate::multicomponent) fn advance_stabilized_with_budget(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedDissipationReport, DissipationError> {
        advance(
            self,
            flow,
            duration,
            controls,
            [None; 2],
            evaluations,
            cancelled,
        )
    }
}

pub(in crate::multicomponent) fn extent(stages: usize) -> f64 {
    let s = stages as f64;
    (s * s + s - 2.0) / 4.0
}

pub(in crate::multicomponent) fn b(stage: usize) -> f64 {
    if stage < 2 {
        1.0 / 3.0
    } else {
        let j = stage as f64;
        (j * j + j - 2.0) / (2.0 * j * (j + 1.0))
    }
}

pub(in crate::multicomponent) fn reserve_evaluation(
    evaluations: &mut usize,
    maximum: usize,
    elapsed: f64,
) -> Result<(), DissipationError> {
    if *evaluations >= maximum {
        return Err(DissipationError::OperatorWorkLimit {
            elapsed_time_s: elapsed,
            operator_evaluations: *evaluations,
        });
    }
    *evaluations += 1;
    Ok(())
}

pub(in crate::multicomponent) fn check_rate(
    step: f64,
    rate: f64,
    maximum: f64,
) -> Result<(), DissipationError> {
    let number = step * rate;
    if !number.is_finite() || number > maximum * (1.0 + 8.0 * f64::EPSILON) {
        return Err(DissipationError::ParabolicNumber {
            observed: number,
            maximum,
        });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(in crate::multicomponent) fn recurrence(
    gas: &MixtureGas,
    mesh: &Mesh1d,
    anchor: &[f64],
    previous: &[f64],
    older_and_next: &mut [f64],
    anchor_faces: &[f64],
    faces: &[f64],
    step: f64,
    coefficients: [f64; 4],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let [mu, nu, mu_tilde, gamma_tilde] = coefficients;
    let m = gas.species_count() + 2;
    let uniform_area = mesh.uniform_cross_section();
    for (i, q) in older_and_next.chunks_exact_mut(m).enumerate() {
        poll(cancelled)?;
        for (k, value) in q.iter_mut().enumerate() {
            let index = i * m + k;
            let base = anchor[index];
            let change = if let Some(area) = uniform_area {
                let scale = step * (area / mesh.volumes()[i]);
                scale
                    * (mu_tilde * (faces[index] - faces[index + m])
                        + gamma_tilde * (anchor_faces[index] - anchor_faces[index + m]))
            } else {
                let left = mesh.face_areas()[i];
                let right = mesh.face_areas()[i + 1];
                (step / mesh.volumes()[i])
                    * (mu_tilde * (left * faces[index] - right * faces[index + m])
                        + gamma_tilde
                            * (left * anchor_faces[index] - right * anchor_faces[index + m]))
            };
            // Differences preserve an exactly stationary state without
            // multiplying a large constant caloric/partial-density offset.
            *value = base + mu * (previous[index] - base) + nu * (*value - base) + change;
        }
        if !gas.admissible(q) {
            return Err(DissipationError::Invalid(
                "stabilized molecular stage left species/caloric domain",
            ));
        }
    }
    Ok(())
}

fn accumulate(total: &mut f64, correction: &mut f64, value: f64) {
    let y = value - *correction;
    let next = *total + y;
    *correction = (next - *total) - y;
    *total = next;
}

pub(super) fn advance<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: StabilizedDissipationControl,
    reservoirs: ReservoirEnds<'_>,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StabilizedDissipationReport, DissipationError> {
    advance_inner::<false, L>(
        owner, flow, duration, controls, reservoirs, evaluations, None, cancelled,
    )
}

/// The same RKL2 recurrence with source-conjugate faces and measured endpoint
/// entropy. Intermediate Legendre states are not independent physical steps.
pub(super) fn advance_entropy_consistent<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: StabilizedDissipationControl,
    reservoirs: ReservoirEnds<'_>,
    account: &mut EntropyTracker,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StabilizedDissipationReport, DissipationError> {
    advance_entropy_consistent_budget(
        owner, flow, duration, controls, reservoirs, &mut 0, account, cancelled,
    )
}

/// Retain one actual RHS budget across both molecular halves and every
/// discarded coupled attempt. Entropy admission does not restart the counter.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance_entropy_consistent_budget<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: StabilizedDissipationControl,
    reservoirs: ReservoirEnds<'_>,
    evaluations: &mut usize,
    account: &mut EntropyTracker,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StabilizedDissipationReport, DissipationError> {
    advance_inner::<true, L>(
        owner, flow, duration, controls, reservoirs, evaluations, Some(account), cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_inner<const ENTROPY_FACES: bool, L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: StabilizedDissipationControl,
    reservoirs: ReservoirEnds<'_>,
    evaluations: &mut usize,
    mut entropy: Option<&mut EntropyTracker>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StabilizedDissipationReport, DissipationError> {
    poll(cancelled)?;
    controls.validate()?;
    if !duration.value().is_finite() || duration.value() < 0.0 {
        return Err(DissipationError::Invalid(
            "nonnegative finite molecular duration required",
        ));
    }
    if flow.gas != owner.gas
        || flow.mesh != owner.mesh
        || flow.boundaries != owner.boundaries
        || owner.law.source() != owner.gas.source()
    {
        return Err(DissipationError::Invalid(
            "stabilized molecular operator and flow ownership differ",
        ));
    }
    let n = owner.gas.species_count();
    let m = n + 2;
    let before = totals(&flow.states, flow.mesh.volumes(), n, cancelled)?;
    let mut boundary = [0.0; MAX_COMPONENTS];
    let mut boundary_correction = [0.0; MAX_COMPONENTS];
    let mut ends = [[0.0; MAX_COMPONENTS]; 2];
    let mut end_corrections = [[0.0; MAX_COMPONENTS]; 2];
    let mut accepted_steps = 0;
    let mut rejected_steps = 0;
    let initial_evaluations = *evaluations;
    let mut elapsed = 0.0;
    let mut maximum_stages_used = 0;
    let mut maximum_parabolic = 0.0_f64;
    let mut maximum_scaled_parabolic = 0.0_f64;
    let mut maximum_gradient_residual = 0.0_f64;
    let mut maximum_mass_frame_residual = 0.0_f64;
    let base = controls.transport;
    let DissipationOperator {
        gas,
        mesh,
        boundaries,
        law,
        maxwell_stefan,
        widths,
        mass_ratio,
        work,
        ..
    } = owner;
    work.current.copy_from_slice(&flow.states);
    if duration.value() > 0.0 && work.stabilized_faces.len() != work.faces.len() {
        work.stabilized_faces = allocated(work.faces.len(), 0.0)?;
    }
    let uniform_area = mesh.uniform_cross_section();
    let end_areas = [mesh.face_areas()[0], mesh.face_areas()[mesh.cell_count()]];
    // A previous aborted moving trial may have changed these scratch widths.
    for (i, width) in widths.iter_mut().enumerate() {
        poll(cancelled)?;
        *width = mesh.uniform_cell_width().map_or_else(
            || mesh.positions()[i + 1] - mesh.positions()[i],
            |w| w.value(),
        );
        if let Some(area) = uniform_area {
            let expected = area * *width;
            if !positive(expected)
                || (mesh.volumes()[i] - expected).abs()
                    > 64.0 * f64::EPSILON * expected.max(mesh.volumes()[i])
            {
                return Err(DissipationError::Invalid(
                    "molecular volume differs from area times physical width",
                ));
            }
        }
    }
    while elapsed < duration.value() {
        poll(cancelled)?;
        if accepted_steps + rejected_steps >= base.maximum_steps {
            return Err(DissipationError::WorkLimit {
                elapsed_time_s: elapsed,
                accepted_steps,
                rejected_steps,
            });
        }
        reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
        let a = spatial_with_wall_motion::<ENTROPY_FACES, _>(
            gas,
            mesh,
            *boundaries,
            reservoirs,
            widths,
            *mass_ratio,
            law,
            maxwell_stefan,
            &work.current,
            &mut work.primitive,
            &mut work.heat_capacity_v,
            &mut work.faces,
            &mut work.face_rates,
            &mut work.binary_diffusivities,
            &mut work.diffusion,
            base.diffusion,
            [0.0; 2],
            cancelled,
        )?;
        if let Some(account) = &mut entropy {
            account.begin_attempt(&work.primitive);
        }
        work.stabilized_faces.copy_from_slice(&work.faces);
        let remaining = duration.value() - elapsed;
        let mut step = remaining.min(base.maximum_step_s);
        if a.rate > 0.0 {
            step = step.min(base.parabolic_cfl * extent(controls.maximum_stages) / a.rate);
        }
        let mut accepted = false;
        let mut last = DissipationError::Invalid("no admissible stabilized molecular step");
        for _ in 0..=base.maximum_halvings {
            poll(cancelled)?;
            if accepted_steps + rejected_steps >= base.maximum_steps {
                return Err(DissipationError::WorkLimit {
                    elapsed_time_s: elapsed,
                    accepted_steps,
                    rejected_steps,
                });
            }
            if !positive(step)
                || elapsed + step == elapsed
                || (step < base.minimum_step_s && remaining > base.minimum_step_s)
            {
                return Err(DissipationError::StepFailure {
                    elapsed_time_s: elapsed,
                    step_s: step,
                    last: Box::new(last),
                });
            }
            let stages = (2..=controls.maximum_stages)
                .find(|s| step * a.rate <= base.parabolic_cfl * extent(*s))
                .unwrap_or(controls.maximum_stages);
            let scheme = Rkl2Scheme::new(stages)?;
            let extension = scheme.stability_extent();
            let limit = base.parabolic_cfl * extension;
            let first_weight = scheme.first_weight();
            let mut previous_ends = [[0.0; MAX_COMPONENTS]; 2];
            let mut older_ends = [[0.0; MAX_COMPONENTS]; 2];
            let mut peak = a;
            let mut entropy_weights = [0.0; 256];
            let attempt = (|| {
                check_rate(step, a.rate, limit)?;
                if let Some(account) = &mut entropy {
                    account.begin_stabilized_quadrature(
                        scheme, &mut entropy_weights, &a, step, end_areas, n, cancelled,
                    )?;
                }
                forward(
                    gas,
                    mesh,
                    &work.current,
                    &work.stabilized_faces,
                    step * first_weight,
                    &mut work.stage,
                    cancelled,
                )?;
                work.candidate.copy_from_slice(&work.current);
                for k in 0..m {
                    previous_ends[0][k] = step * first_weight * end_areas[0] * a.left[k];
                    previous_ends[1][k] = -step * first_weight * end_areas[1] * a.right[k];
                }
                // Evaluate every intermediate and the final stage with the
                // same original source, caloric and constitutive admission.
                for j in 2..=stages + 1 {
                    poll(cancelled)?;
                    reserve_evaluation(
                        evaluations,
                        controls.maximum_operator_evaluations,
                        elapsed,
                    )?;
                    let evaluated = spatial_with_wall_motion::<ENTROPY_FACES, _>(
                        gas,
                        mesh,
                        *boundaries,
                        reservoirs,
                        widths,
                        *mass_ratio,
                        law,
                        maxwell_stefan,
                        &work.stage,
                        &mut work.primitive,
                        &mut work.heat_capacity_v,
                        &mut work.faces,
                        &mut work.face_rates,
                        &mut work.binary_diffusivities,
                        &mut work.diffusion,
                        base.diffusion,
                        [0.0; 2],
                        cancelled,
                    )?;
                    check_rate(step, evaluated.rate, limit)?;
                    peak.rate = peak.rate.max(evaluated.rate);
                    peak.gradient_residual =
                        peak.gradient_residual.max(evaluated.gradient_residual);
                    peak.mass_frame_residual =
                        peak.mass_frame_residual.max(evaluated.mass_frame_residual);
                    if j > stages {
                        break;
                    }
                    if let Some(account) = &mut entropy {
                        account.observe_stabilized_quadrature(
                            &evaluated, step * entropy_weights[j - 1],
                            end_areas, n, cancelled,
                        )?;
                    }
                    let [mu, nu, mu_tilde, gamma_tilde] = scheme.coefficients(j)?;
                    recurrence(
                        gas,
                        mesh,
                        &work.current,
                        &work.stage,
                        &mut work.candidate,
                        &work.stabilized_faces,
                        &work.faces,
                        step,
                        [mu, nu, mu_tilde, gamma_tilde],
                        cancelled,
                    )?;
                    for k in 0..m {
                        for (end, (now, anchor)) in [
                            (evaluated.left[k], a.left[k]),
                            (-evaluated.right[k], -a.right[k]),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            let next = mu * previous_ends[end][k]
                                + nu * older_ends[end][k]
                                + step * end_areas[end] * (mu_tilde * now + gamma_tilde * anchor);
                            older_ends[end][k] = previous_ends[end][k];
                            previous_ends[end][k] = next;
                        }
                    }
                    core::mem::swap(&mut work.stage, &mut work.candidate);
                }
                if let Some(account) = &mut entropy {
                    // The final native RHS admission already inverted this
                    // exact endpoint. No extra constitutive or EOS query is
                    // needed, and refused endpoints cannot alter the anchor.
                    account.check_stabilized_candidate(
                        gas, mesh.volumes(), &work.current, &work.stage,
                        &work.primitive, &previous_ends, cancelled,
                    )?;
                }
                Ok::<_, DissipationError>(())
            })();
            match attempt {
                Ok(()) => {
                    if let Some(account) = &mut entropy {
                        account.accept_stabilized_candidate(n, &work.primitive)?;
                    }
                    for k in 0..m {
                        accumulate(
                            &mut boundary[k],
                            &mut boundary_correction[k],
                            previous_ends[0][k] + previous_ends[1][k],
                        );
                        for end in 0..2 {
                            accumulate(
                                &mut ends[end][k],
                                &mut end_corrections[end][k],
                                previous_ends[end][k],
                            );
                        }
                    }
                    maximum_stages_used = maximum_stages_used.max(stages);
                    maximum_parabolic = maximum_parabolic.max(step * peak.rate);
                    maximum_scaled_parabolic =
                        maximum_scaled_parabolic.max(step * peak.rate / extension);
                    maximum_gradient_residual =
                        maximum_gradient_residual.max(peak.gradient_residual);
                    maximum_mass_frame_residual =
                        maximum_mass_frame_residual.max(peak.mass_frame_residual);
                    core::mem::swap(&mut work.current, &mut work.stage);
                    elapsed = if step == remaining {
                        duration.value()
                    } else {
                        elapsed + step
                    };
                    accepted_steps += 1;
                    accepted = true;
                    break;
                }
                Err(
                    error @ (DissipationError::Cancelled
                    | DissipationError::OperatorWorkLimit { .. }),
                ) => {
                    return Err(error);
                }
                Err(error) => {
                    rejected_steps += 1;
                    last = error;
                    step *= 0.5;
                }
            }
        }
        if !accepted {
            return Err(DissipationError::StepFailure {
                elapsed_time_s: elapsed,
                step_s: step,
                last: Box::new(last),
            });
        }
    }
    let after = totals(&work.current, mesh.volumes(), n, cancelled)?;
    let mut residual = [0.0; MAX_COMPONENTS];
    for k in 0..m {
        residual[k] = after.get(k) - before.get(k) - boundary[k];
    }
    let absolute_energy_residual = sum(core::iter::once(residual[n + 1])
        .chain((0..n).map(|k| residual[k] * gas.energy_references()[k])));
    if !boundary[..m]
        .iter()
        .chain(&residual[..m])
        .chain(&ends[0][..m])
        .chain(&ends[1][..m])
        .all(|x| x.is_finite())
        || !absolute_energy_residual.is_finite()
    {
        return Err(DissipationError::Invalid(
            "non-finite stabilized molecular extensive balance",
        ));
    }
    if let Some(account) = &mut entropy {
        // One allowance applies to the complete requested interval as well
        // as each accepted macrostep; rejected work cannot publish a prefix.
        account.finish(gas, mesh.volumes(), &flow.states, &work.current, cancelled)?;
    }
    poll(cancelled)?;
    flow.states.copy_from_slice(&work.current);
    Ok(StabilizedDissipationReport {
        transport: DissipationReport {
            elapsed: duration,
            accepted_steps,
            rejected_steps,
            maximum_parabolic_number: maximum_parabolic,
            maximum_gradient_residual_per_m: maximum_gradient_residual,
            maximum_mass_frame_residual_kg_m2_s: maximum_mass_frame_residual,
            before,
            after,
            boundary_change: Totals::from_values(&boundary[..m], n),
            left_boundary_change: Totals::from_values(&ends[0][..m], n),
            right_boundary_change: Totals::from_values(&ends[1][..m], n),
            conservation_residual: Totals::from_values(&residual[..m], n),
            absolute_energy_residual_j: absolute_energy_residual,
        },
        operator_evaluations: *evaluations - initial_evaluations,
        maximum_stages_used,
        maximum_scaled_parabolic_number: maximum_scaled_parabolic,
    })
}
