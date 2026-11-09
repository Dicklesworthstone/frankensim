//! Continuous sinusoidal loads on the existing complete tensor-memory flow.
//!
//! Source records obey c'=omega*s, s'=-omega*c. They carry applied force, not
//! material storage. The original native quadratic integral observes F(t) dot v;
//! independent rotating velocity observers expose its conjugate work moments.

use super::*;
use fs_phs::linear_exact::exponential_with_budget;

/// Accepted physical motion and moments conjugate to the supplied coefficients.
/// Each outer vector follows angular_frequencies_rad_s(), each inner vector
/// follows the unchanged original displacement coordinates.
#[derive(Clone, Debug)]
pub struct PreparedTensorToneStep {
    pub motion: PreparedTensorStep,
    /// Integral_0^1 cos(omega*h*s) v(h*s) ds [m/s].
    pub cosine_velocity_moment_mps: Vec<Vec<f64>>,
    /// Integral_0^1 sin(omega*h*s) v(h*s) ds [m/s].
    pub sine_velocity_moment_mps: Vec<Vec<f64>>,
}

pub(super) struct TonePrepared {
    pub(super) frequencies: Vec<f64>,
    force_scale_n: f64,
    rotations: Vec<(f64, f64)>,
    // Original-order cosine/sine moment rows acting on conditioned full state.
    velocity_moments: Vec<f64>,
}

impl PreparedTensorBody {
    /// Prepare F(t0+t)=sum(A_j*cos(omega_j*t)+B_j*sin(omega_j*t)).
    ///
    /// Every original coordinate, tensor point and exact material time remains.
    /// Force coefficients are supplied relative to each interval's start; finite
    /// changes between intervals cause no impulse. Each authored frequency must
    /// lie in every original material law's admitted band. This is not a claim
    /// about transient content outside that band or strain inside an interval.
    ///
    /// A positive force scale conditions source records only. All source and
    /// observer coordinates count against the existing explicit native budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn from_mass_matrix_tones(
        coordinates: usize,
        mass_kg: Vec<f64>,
        points: Vec<ProjectedTensorProny>,
        heat_groups: Vec<usize>,
        dt_s: f64,
        angular_frequencies_rad_s: &[f64],
        force_scale_n: f64,
        controls: PreparedTensorControls,
        gate: &CancelGate,
    ) -> Result<Self> {
        Self::prepare_sources(
            coordinates,
            mass_kg,
            points,
            heat_groups,
            dt_s,
            Some(force_scale_n),
            Some(angular_frequencies_rad_s),
            controls,
            gate,
        )
    }

    pub fn angular_frequencies_rad_s(&self) -> Option<&[f64]> {
        self.tones
            .as_ref()
            .map(|tones| tones.frequencies.as_slice())
    }

    /// Rotate unchanged physical source amplitudes to the next interval origin.
    /// No mechanical history advances. Fresh coefficients can also be supplied
    /// directly, for example when an external source changes its amplitudes.
    pub fn advance_force_phase(
        &self,
        cosine_force_n: &[Vec<f64>],
        sine_force_n: &[Vec<f64>],
        gate: &CancelGate,
    ) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        self.rotate_force_phase(cosine_force_n, sine_force_n, self.coordinates(), gate)
    }

    pub(crate) fn rotate_force_phase(
        &self,
        cosine: &[Vec<f64>],
        sine: &[Vec<f64>],
        width: usize,
        gate: &CancelGate,
    ) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        poll(gate)?;
        let tones = self
            .tones
            .as_ref()
            .ok_or_else(|| invalid("force phase requires continuous-tone preparation"))?;
        check_coefficients(cosine, sine, tones.frequencies.len(), width)?;
        let mut next_cosine = Vec::with_capacity(cosine.len());
        let mut next_sine = Vec::with_capacity(sine.len());
        for (tone, &(s, c)) in tones.rotations.iter().enumerate() {
            poll(gate)?;
            let mut a = zeros(width)?;
            let mut b = zeros(width)?;
            for i in 0..width {
                a[i] = c.mul_add(cosine[tone][i], s * sine[tone][i]);
                b[i] = (-s).mul_add(cosine[tone][i], c * sine[tone][i]);
            }
            if !finite(&a) || !finite(&b) {
                return Err(invalid("continuous force phase overflow"));
            }
            next_cosine.push(a);
            next_sine.push(b);
        }
        Ok((next_cosine, next_sine))
    }
}

fn check_coefficients(
    cosine: &[Vec<f64>],
    sine: &[Vec<f64>],
    tones: usize,
    coordinates: usize,
) -> Result<()> {
    if cosine.len() != tones
        || sine.len() != tones
        || cosine
            .iter()
            .chain(sine)
            .any(|v| v.len() != coordinates || !finite(v))
    {
        return Err(invalid(
            "one complete finite cosine/sine force field per tone required",
        ));
    }
    Ok(())
}

/// A one-way oscillator per observed velocity gives both weighted integrals
/// using one native exponential. Its endpoint is rotated BACK by omega*h:
/// I_cos=c*r_cos-s*r_sin, I_sin=s*r_cos+c*r_sin. Observers have no feedback.
#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_observer(
    n: usize,
    generator: &[f64],
    scales: &[f64],
    inertia: &Inertia,
    frequencies: &[f64],
    dt_s: f64,
    force_scale_n: f64,
    budget: DenseExponentialBudget,
    gate: &CancelGate,
) -> Result<TonePrepared> {
    poll(gate)?;
    let states = scales.len();
    let observations = 2 * n * frequencies.len();
    let dimension = states + observations;
    let mut operator = zeros(dimension * dimension)?;
    for i in 0..states {
        poll(gate)?;
        for j in 0..states {
            operator[i * dimension + j] = dt_s * generator[i * states + j];
        }
    }
    let mut rotations = Vec::with_capacity(frequencies.len());
    for (tone, &omega) in frequencies.iter().enumerate() {
        poll(gate)?;
        rotations.push(fs_math::det::sin_cos(omega * dt_s));
        for i in 0..n {
            let cosine = states + (2 * tone) * n + i;
            let sine = cosine + n;
            operator[cosine * dimension + sine] = dt_s * omega;
            operator[sine * dimension + cosine] = -dt_s * omega;
            for j in 0..n {
                operator[cosine * dimension + n + j] =
                    dt_s * (inertia.inverse_entry(i, j) / scales[n + j]) * scales[i];
            }
        }
    }
    if !finite(&operator) {
        return Err(invalid("continuous velocity observer overflow"));
    }
    let exponential = exponential_with_budget(&operator, dimension, budget)?;
    poll(gate)?;
    let mut velocity_moments = zeros(observations * states)?;
    for (tone, &(s, c)) in rotations.iter().enumerate() {
        for i in 0..n {
            poll(gate)?;
            let cosine = (2 * tone) * n + i;
            let sine = cosine + n;
            for j in 0..states {
                let a = exponential[(states + cosine) * dimension + j];
                let b = exponential[(states + sine) * dimension + j];
                velocity_moments[cosine * states + j] = c.mul_add(a, -s * b) / scales[i] / dt_s;
                velocity_moments[sine * states + j] = s.mul_add(a, c * b) / scales[i] / dt_s;
            }
        }
    }
    if !finite(&velocity_moments) {
        return Err(invalid("continuous velocity moments overflow"));
    }
    Ok(TonePrepared {
        frequencies: frequencies.to_vec(),
        force_scale_n,
        rotations,
        velocity_moments,
    })
}

impl PreparedTensorHistory<'_> {
    /// Advance the whole continuous multi-tone force, including all interference
    /// terms in physical loss and work. No endpoint sample or independent-tone
    /// superposition is used. Refusal preserves motion, memory, clock and ledgers.
    pub fn advance_tones(
        &mut self,
        cosine_force_n: &[Vec<f64>],
        sine_force_n: &[Vec<f64>],
        gate: &CancelGate,
    ) -> Result<PreparedTensorToneStep> {
        poll(gate)?;
        let body = self.body;
        let tones = body
            .tones
            .as_ref()
            .ok_or_else(|| invalid("continuous tones require explicit preparation"))?;
        let n = body.coordinates();
        let physical = body.states;
        let propagated = body.propagation_state_dimension();
        check_coefficients(cosine_force_n, sine_force_n, tones.frequencies.len(), n)?;
        let step = self
            .steps
            .checked_add(1)
            .filter(|s| *s <= body.controls.maximum_steps)
            .ok_or_else(|| invalid("accepted step budget exhausted"))?;
        let time_s = step as f64 * body.dt_s;
        if !time_s.is_finite() || time_s <= self.time_s() {
            return Err(invalid("representable later fixed clock required"));
        }
        self.scratch.source_state[..physical].copy_from_slice(&self.state);
        for tone in 0..tones.frequencies.len() {
            let at = physical + 2 * tone * n;
            self.scratch.source_state[at..at + n].copy_from_slice(&cosine_force_n[tone]);
            self.scratch.source_state[at + n..at + 2 * n].copy_from_slice(&sine_force_n[tone]);
        }
        // The complete forces are rotating source states; held inputs are zero.
        self.scratch.source_increment.fill(0.0);
        body.flow.evaluate_into(
            &self.scratch.source_state,
            &self.scratch.source_increment,
            &mut self.scratch.augmented,
            &mut self.scratch.balanced,
            &mut self.scratch.endpoint,
            &mut self.scratch.integrals,
        )?;
        poll(gate)?;
        for (tone, &(s, c)) in tones.rotations.iter().enumerate() {
            for i in 0..n {
                let a = cosine_force_n[tone][i];
                let b = sine_force_n[tone][i];
                let expected = [c.mul_add(a, s * b), (-s).mul_add(a, c * b)];
                for (part, value) in expected.into_iter().enumerate() {
                    let actual = self.scratch.endpoint[physical + (2 * tone + part) * n + i];
                    let scale = tones
                        .force_scale_n
                        .max(a.abs())
                        .max(b.abs())
                        .max(value.abs());
                    let allowance = body.controls.maximum_relative_operator_error * scale;
                    let defect = (actual - value).abs();
                    if !defect.is_finite() || defect > allowance {
                        return Err(PreparedTensorError::Operator { defect, allowance });
                    }
                }
            }
        }
        self.scratch
            .physical_candidate
            .copy_from_slice(&self.scratch.endpoint[..physical]);
        body.admit_displacement(&self.scratch.physical_candidate[..n], gate)?;
        let velocity = body
            .inertia
            .velocity(&self.scratch.physical_candidate[n..2 * n]);
        if !finite(&velocity) {
            return Err(invalid("candidate velocity overflow"));
        }
        let after = body.energy(&self.scratch.physical_candidate, &velocity, gate)?;
        let change = body.energy_change(
            &self.state,
            &self.velocity,
            &self.scratch.physical_candidate,
            &velocity,
            gate,
        )?;
        for i in 0..propagated {
            self.scratch.heat_state[i] = body.heat_scales[i] * self.scratch.source_state[i];
        }
        let mut cosine_velocity_moment_mps = Vec::with_capacity(tones.frequencies.len());
        let mut sine_velocity_moment_mps = Vec::with_capacity(tones.frequencies.len());
        let (mut moment_work, mut moment_correction) = (0.0, 0.0);
        for tone in 0..tones.frequencies.len() {
            poll(gate)?;
            let mut cosine = zeros(n)?;
            let mut sine = zeros(n)?;
            for i in 0..n {
                for (part, value) in [&mut cosine[i], &mut sine[i]].into_iter().enumerate() {
                    let row = &tones.velocity_moments[((2 * tone + part) * n + i) * propagated
                        ..((2 * tone + part) * n + i + 1) * propagated];
                    let (mut sum, mut correction) = (0.0, 0.0);
                    for (&a, &b) in row.iter().zip(&self.scratch.heat_state) {
                        add(&mut sum, &mut correction, a * b);
                    }
                    *value = sum;
                }
                add(
                    &mut moment_work,
                    &mut moment_correction,
                    body.dt_s * cosine_force_n[tone][i] * cosine[i],
                );
                add(
                    &mut moment_work,
                    &mut moment_correction,
                    body.dt_s * sine_force_n[tone][i] * sine[i],
                );
            }
            if !finite(&cosine) || !finite(&sine) {
                return Err(invalid("continuous velocity moment is nonfinite"));
            }
            cosine_velocity_moment_mps.push(cosine);
            sine_velocity_moment_mps.push(sine);
        }
        let supplied_work_j = self.scratch.integrals[0];
        let defect_j = moment_work - supplied_work_j;
        if !defect_j.is_finite() || defect_j.abs() > body.controls.maximum_work_error_j {
            return Err(PreparedTensorError::Energy {
                defect_j,
                allowance_j: body.controls.maximum_work_error_j,
            });
        }
        let mut heat_by_group_j = zeros(body.heat.len())?;
        let (mut dissipated_j, mut heat_roundoff) = (0.0, 0.0);
        for (out, cache) in heat_by_group_j.iter_mut().zip(&body.heat) {
            poll(gate)?;
            *out = cache.evaluate_into_controlled(
                &self.scratch.heat_state,
                &self.scratch.source_increment,
                &mut self.scratch.heat_coordinates,
                || gate.is_requested(),
            )?;
            if !out.is_finite() || *out < 0.0 {
                return Err(invalid(
                    "original source heat must be finite and nonnegative",
                ));
            }
            add(&mut dissipated_j, &mut heat_roundoff, *out);
        }
        let (mut work_j, mut work_correction_j) = (self.work_j, self.work_correction_j);
        let (mut heat_j, mut heat_correction_j) = (self.heat_j, self.heat_correction_j);
        add(&mut work_j, &mut work_correction_j, supplied_work_j);
        add(&mut heat_j, &mut heat_correction_j, dissipated_j);
        let mut mean_velocity_mps = zeros(n)?;
        for (i, value) in mean_velocity_mps.iter_mut().enumerate() {
            *value = (self.scratch.physical_candidate[i] - self.state[i]) / body.dt_s;
        }
        if !finite(&[
            supplied_work_j,
            dissipated_j,
            work_j,
            work_correction_j,
            heat_j,
            heat_correction_j,
        ]) || !finite(&mean_velocity_mps)
            || heat_j < self.heat_j
        {
            return Err(invalid(
                "physical continuous-tone work/heat ledger overflow",
            ));
        }
        let energy = PreparedTensorEnergy {
            equilibrium_before_j: self.energy.equilibrium,
            equilibrium_after_j: after.equilibrium,
            kinetic_before_j: self.energy.kinetic,
            kinetic_after_j: after.kinetic,
            memory_before_j: self.energy.memory,
            memory_after_j: after.memory,
            supplied_work_j,
            dissipated_j,
            cumulative_work_j: work_j,
            cumulative_parameter_work_j: self.parameter_work_j,
            cumulative_dissipated_j: heat_j,
            balance_error_j: supplied_work_j - dissipated_j - change,
            cumulative_balance_error_j: work_j + self.parameter_work_j
                - heat_j
                - (after.total() - self.initial_energy_j),
        };
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > body.controls.maximum_work_error_j {
                return Err(PreparedTensorError::Energy {
                    defect_j,
                    allowance_j: body.controls.maximum_work_error_j,
                });
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.state, &mut self.scratch.physical_candidate);
        self.velocity = velocity;
        self.energy = after;
        self.steps = step;
        self.work_j = work_j;
        self.work_correction_j = work_correction_j;
        self.heat_j = heat_j;
        self.heat_correction_j = heat_correction_j;
        Ok(PreparedTensorToneStep {
            motion: PreparedTensorStep {
                step,
                time_s,
                mean_velocity_mps,
                heat_by_group_j,
                energy,
            },
            cosine_velocity_moment_mps,
            sine_velocity_moment_mps,
        })
    }
}
