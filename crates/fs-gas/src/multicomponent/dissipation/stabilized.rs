//! Opt-in second-order Runge--Kutta--Legendre molecular super-time-stepping.
//!
//! The recurrence and negative-real stability interval follow Meyer, Balsara
//! and Aslam (JCP 2014, doi:10.1016/j.jcp.2013.08.021), also written explicitly
//! in Vaidya et al. (MNRAS 2017, doi:10.1093/mnras/stx2176), section 3.3 and
//! appendix A. Only time integration changes: every stage uses the original
//! finite-volume Maxwell--Stefan/enthalpy/heat/stress face owner.

use super::super::evolve::totals;
use super::super::{Mesh1d, MixtureEuler1d, MixtureGas, Totals, allocated, sum};
use super::evolve::{forward, spatial};
use super::{
    DissipationControl, DissipationError, DissipationOperator, DissipationReport, MAX_COMPONENTS,
    TransportLaw, poll, positive,
};
use fs_qty::Time;

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
        advance(self, flow, duration, controls, &mut cancelled)
    }
}

fn extent(stages: usize) -> f64 {
    let s = stages as f64;
    (s * s + s - 2.0) / 4.0
}

fn b(stage: usize) -> f64 {
    if stage < 2 {
        1.0 / 3.0
    } else {
        let j = stage as f64;
        (j * j + j - 2.0) / (2.0 * j * (j + 1.0))
    }
}

fn reserve_evaluation(
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

fn check_rate(step: f64, rate: f64, maximum: f64) -> Result<(), DissipationError> {
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
fn recurrence(
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
    for (i, q) in older_and_next.chunks_exact_mut(m).enumerate() {
        poll(cancelled)?;
        let scale = step * (mesh.face_areas()[0] / mesh.volumes()[i]);
        for (k, value) in q.iter_mut().enumerate() {
            let index = i * m + k;
            let base = anchor[index];
            // Differences preserve an exactly stationary state without
            // multiplying a large constant caloric/partial-density offset.
            *value = base
                + mu * (previous[index] - base)
                + nu * (*value - base)
                + scale
                    * (mu_tilde * (faces[index] - faces[index + m])
                        + gamma_tilde * (anchor_faces[index] - anchor_faces[index + m]));
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

fn advance<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: StabilizedDissipationControl,
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
    let mut evaluations = 0;
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
    let area = mesh.face_areas()[0];
    // A previous aborted moving trial may have changed these scratch widths.
    for (i, width) in widths.iter_mut().enumerate() {
        poll(cancelled)?;
        *width = mesh.uniform_cell_width().map_or_else(
            || mesh.positions()[i + 1] - mesh.positions()[i],
            |w| w.value(),
        );
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
    while elapsed < duration.value() {
        poll(cancelled)?;
        if accepted_steps + rejected_steps >= base.maximum_steps {
            return Err(DissipationError::WorkLimit {
                elapsed_time_s: elapsed,
                accepted_steps,
                rejected_steps,
            });
        }
        reserve_evaluation(
            &mut evaluations,
            controls.maximum_operator_evaluations,
            elapsed,
        )?;
        let a = spatial(
            gas,
            mesh,
            *boundaries,
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
            cancelled,
        )?;
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
            let extension = extent(stages);
            let limit = base.parabolic_cfl * extension;
            let w1 = 1.0 / extension;
            let first_weight = b(1) * w1;
            let mut previous_ends = [[0.0; MAX_COMPONENTS]; 2];
            let mut older_ends = [[0.0; MAX_COMPONENTS]; 2];
            let mut peak = a;
            let attempt = (|| {
                check_rate(step, a.rate, limit)?;
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
                    previous_ends[0][k] = step * first_weight * area * a.left[k];
                    previous_ends[1][k] = -step * first_weight * area * a.right[k];
                }
                // Evaluate every intermediate and the final stage with the
                // same original source, caloric and constitutive admission.
                for j in 2..=stages + 1 {
                    poll(cancelled)?;
                    reserve_evaluation(
                        &mut evaluations,
                        controls.maximum_operator_evaluations,
                        elapsed,
                    )?;
                    let evaluated = spatial(
                        gas,
                        mesh,
                        *boundaries,
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
                    let jf = j as f64;
                    let mu = ((2.0 * jf - 1.0) / jf) * b(j) / b(j - 1);
                    let nu = -((jf - 1.0) / jf) * b(j) / b(j - 2);
                    let mu_tilde = mu * w1;
                    let gamma_tilde = -(1.0 - b(j - 1)) * mu_tilde;
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
                                + step * area * (mu_tilde * now + gamma_tilde * anchor);
                            older_ends[end][k] = previous_ends[end][k];
                            previous_ends[end][k] = next;
                        }
                    }
                    core::mem::swap(&mut work.stage, &mut work.candidate);
                }
                Ok::<_, DissipationError>(())
            })();
            match attempt {
                Ok(()) => {
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
        operator_evaluations: evaluations,
        maximum_stages_used,
        maximum_scaled_parabolic_number: maximum_scaled_parabolic,
    })
}
