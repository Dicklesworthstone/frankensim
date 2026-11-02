//! Time-domain rotor motion at one declared constant shaft speed.
//!
//! Implicit midpoint advances the ORIGINAL full equation
//! `M v' + (C + Omega G) v + K q = f`, with `q' = v`. Native real
//! factorizations are cached: preparation costs O(n^3), subsequent steps and
//! instantaneous impulse responses cost O(n^2), with n <= 32. The step solves
//! for a velocity increment, avoiding subtraction of a large unchanged M v.
//!
//! Skew gyroscopic force changes motion without physical dissipation. Bearing
//! C and K retain all cross-coupled and active entries; their signed work is
//! not called heat. For symmetric K and C=0, midpoint preserves the quadratic
//! mechanical energy up to solve/rounding error; its phase still needs timestep
//! refinement. A caller-supplied midpoint force is second-order only when it
//! approximates the actual smooth load there. Spin ramps, nonlinear bearings,
//! shaft-model validity and event timing remain outside this constant-speed lane.

use super::{MAX_ROTOR_DOFS, RotorBearingModel, RotorError, RotorOperatingPoint, finite, poll};
use fs_la::factor::{Cholesky, FactorError, Lu, cholesky, lu};
use fs_math::dd::Dd;

/// Original generalized displacement/angle and rate coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorTransientState {
    /// Original displacements/angles, using the model's coordinate units.
    pub displacement: Vec<f64>,
    /// Their time derivatives, conjugate to the supplied forces/torques.
    pub velocity: Vec<f64>,
}

/// A complete accepted finite step, including work of the original operators.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorTransientStep {
    /// Reusable physical state; the borrowed input is never modified.
    pub state_after: RotorTransientState,
    /// Relative residual of the ORIGINAL impulse balance after diagonal mass
    /// scaling; zero-force rows retain the complete vector's roundoff scale.
    pub relative_residual: f64,
    /// Signed generalized load work at the actual midpoint.
    pub input_work_j: f64,
    /// Signed work against C v. Active bearings can make this negative.
    pub velocity_force_work_j: f64,
    /// Signed work against K q. Only symmetric K defines elastic storage.
    pub stiffness_force_work_j: f64,
    /// Measured work against Omega G v; zero in exact arithmetic.
    pub gyroscopic_work_j: f64,
    /// Polarized kinetic change from the actual returned velocity increment.
    pub kinetic_energy_change_j: f64,
    /// Kinetic change plus original operator work minus supplied work.
    pub work_energy_residual_j: f64,
}

/// Instantaneous physical generalized impulse, with unchanged configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorImpulseResponse {
    /// Same displacements with the accepted full-inertia velocity response.
    pub state_after: RotorTransientState,
    /// Original mass-scaled momentum-increment residual, after endpoint rounding.
    pub relative_residual: f64,
    /// Physical impulse paired with the actual midpoint velocity.
    pub impulse_work_j: f64,
    /// Polarized kinetic change; includes the quadratic impulse energy.
    pub kinetic_energy_change_j: f64,
}

/// Reusable fixed-duration stepper borrowing one immutable physical model.
#[derive(Debug, Clone)]
pub struct RotorMidpoint<'a> {
    model: &'a RotorBearingModel,
    spin_rad_s: f64,
    duration_seconds: f64,
    relative_tolerance: f64,
    velocity_operator: Vec<f64>,
    mass_scale: Vec<f64>,
    mass_factor: Cholesky,
    step_row_scale: Vec<f64>,
    step_factor: Lu,
}

impl<'a> RotorOperatingPoint<'a> {
    /// Prepare a midpoint trajectory and physical-impulse operator once.
    ///
    /// M is scaled by its original diagonal, and the complete nonsymmetric
    /// step matrix is row-equilibrated before native partial-pivot LU. No
    /// matrix is symmetrized, diagonalized, inverted or regularized. A singular
    /// step refuses. `relative_tolerance` in (0,1) admits the actual original
    /// momentum equation after each solve, including endpoint rounding.
    /// Native dense factorizations/solves are indivisible but bounded to 32
    /// coordinates; cancellation is polled around them and per assembly row.
    pub fn midpoint_stepper(
        &self,
        duration_seconds: f64,
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMidpoint<'a>, RotorError> {
        poll(&mut cancelled)?;
        let half = 0.5 * duration_seconds;
        if !duration_seconds.is_finite()
            || half <= 0.0
            || !relative_tolerance.is_finite()
            || relative_tolerance <= 0.0
            || relative_tolerance >= 1.0
        {
            return Err(RotorError::Invalid(
                "positive representable duration and relative tolerance in (0,1) required",
            ));
        }
        let m = self.model;
        let mass_scale: Vec<_> = (0..m.n)
            .map(|i| fs_math::det::sqrt(m.mass[i * m.n + i]))
            .collect();
        let mut scaled_mass = vec![0.0; m.n * m.n];
        let mut matrix = vec![0.0; m.n * m.n];
        let mut step_row_scale = vec![0.0_f64; m.n];
        for i in 0..m.n {
            poll(&mut cancelled)?;
            for j in 0..m.n {
                let ij = i * m.n + j;
                scaled_mass[ij] = finite((m.mass[ij] / mass_scale[i]) / mass_scale[j])?;
                let entry = finite(half.mul_add(
                    half * m.stiffness[ij],
                    half.mul_add(self.velocity_operator[ij], m.mass[ij]),
                ))?;
                matrix[ij] = finite((entry / mass_scale[i]) / mass_scale[j])?;
                step_row_scale[i] = step_row_scale[i].max(matrix[ij].abs());
            }
            if step_row_scale[i] == 0.0 {
                return Err(RotorError::TransientLinear(FactorError::Singular {
                    index: i,
                }));
            }
            for entry in &mut matrix[i * m.n..(i + 1) * m.n] {
                *entry /= step_row_scale[i];
            }
        }
        let mass_factor = cholesky(&scaled_mass, m.n).map_err(RotorError::TransientLinear)?;
        poll(&mut cancelled)?;
        let step_factor = lu(&matrix, m.n).map_err(RotorError::TransientLinear)?;
        finite(step_factor.growth)?;
        poll(&mut cancelled)?;
        Ok(RotorMidpoint {
            model: m,
            spin_rad_s: self.spin_rad_s,
            duration_seconds,
            relative_tolerance,
            velocity_operator: self.velocity_operator.clone(),
            mass_scale,
            mass_factor,
            step_row_scale,
            step_factor,
        })
    }
}

impl RotorMidpoint<'_> {
    /// Fixed physical duration of each accepted trajectory step.
    pub const fn duration_seconds(&self) -> f64 {
        self.duration_seconds
    }

    /// Advance the actual rotor coordinates under a force/torque sampled at
    /// the temporal midpoint. All returned work uses actual endpoint midpoints.
    /// No input is changed if the solve, original balance or cancellation refuses.
    pub fn step(
        &self,
        state: &RotorTransientState,
        midpoint_force: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorTransientStep, RotorError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(midpoint_force)?;
        let n = self.model.n;
        let h = self.duration_seconds;
        let mut half_position = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in half_position.iter_mut().enumerate().take(n) {
            *value = (0.5 * h).mul_add(state.velocity[i], state.displacement[i]);
        }
        self.admit_vector(&half_position[..n])?;
        let velocity_force = action(&self.velocity_operator, &state.velocity, &mut cancelled)?;
        let spring_force = action(&self.model.stiffness, &half_position[..n], &mut cancelled)?;
        let mut increment = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            increment[i] = ((h * (midpoint_force[i] - velocity_force[i] - spring_force[i]))
                / self.mass_scale[i])
                / self.step_row_scale[i];
        }
        self.admit_vector(&increment[..n])?;
        self.step_factor.solve(&mut increment[..n]);
        poll(&mut cancelled)?;
        let velocity: Vec<_> = state
            .velocity
            .iter()
            .enumerate()
            .map(|(i, &v)| v + increment[i] / self.mass_scale[i])
            .collect();
        let midpoint_velocity = midpoint(&state.velocity, &velocity);
        let displacement = state
            .displacement
            .iter()
            .zip(&midpoint_velocity)
            .map(|(&q, &v)| h.mul_add(v, q))
            .collect();
        let state_after = RotorTransientState {
            displacement,
            velocity,
        };
        self.admit_state(&state_after)?;
        let midpoint_position = midpoint(&state.displacement, &state_after.displacement);
        let mut delta = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in delta.iter_mut().enumerate().take(n) {
            *value = state_after.velocity[i] - state.velocity[i];
        }
        let inertial = action(&self.model.mass, &delta[..n], &mut cancelled)?;
        let damping = scaled_action(
            &self.model.damping,
            &midpoint_velocity[..n],
            h,
            &mut cancelled,
        )?;
        let gyro = scaled_action(
            &self.model.gyro,
            &midpoint_velocity[..n],
            finite(h * self.spin_rad_s)?,
            &mut cancelled,
        )?;
        let stiffness = scaled_action(
            &self.model.stiffness,
            &midpoint_position[..n],
            h,
            &mut cancelled,
        )?;
        let mut input = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            input[i] = h * midpoint_force[i];
        }
        self.admit_vector(&input[..n])?;
        let relative_residual = self.admit_balance(
            &[&inertial[..n], &damping[..n], &gyro[..n], &stiffness[..n]],
            &input[..n],
        )?;
        let input_work_j = pairing(&input[..n], &midpoint_velocity[..n])?;
        let velocity_force_work_j = pairing(&damping[..n], &midpoint_velocity[..n])?;
        let stiffness_force_work_j = pairing(&stiffness[..n], &midpoint_velocity[..n])?;
        let gyroscopic_work_j = pairing(&gyro[..n], &midpoint_velocity[..n])?;
        let kinetic_energy_change_j = pairing(&inertial[..n], &midpoint_velocity[..n])?;
        let work_energy_residual_j = finite(
            (Dd::from_f64(kinetic_energy_change_j)
                + Dd::from_f64(velocity_force_work_j)
                + Dd::from_f64(stiffness_force_work_j)
                + Dd::from_f64(gyroscopic_work_j)
                - Dd::from_f64(input_work_j))
            .to_f64(),
        )?;
        poll(&mut cancelled)?;
        Ok(RotorTransientStep {
            state_after,
            relative_residual,
            input_work_j,
            velocity_force_work_j,
            stiffness_force_work_j,
            gyroscopic_work_j,
            kinetic_energy_change_j,
            work_energy_residual_j,
        })
    }

    /// Apply an instantaneous generalized impulse through the ORIGINAL full M.
    /// Bearings and gyroscopic forces are finite forces and supply no impulsive
    /// correction. Positions remain unchanged. An unrepresentable momentum
    /// increment refuses the original balance instead of silently clipping it.
    pub fn apply_impulse(
        &self,
        state: &RotorTransientState,
        impulse: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorImpulseResponse, RotorError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(impulse)?;
        let n = self.model.n;
        let mut increment = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            increment[i] = impulse[i] / self.mass_scale[i];
        }
        self.admit_vector(&increment[..n])?;
        self.mass_factor.solve(&mut increment[..n]);
        poll(&mut cancelled)?;
        let state_after = RotorTransientState {
            displacement: state.displacement.clone(),
            velocity: state
                .velocity
                .iter()
                .enumerate()
                .map(|(i, &v)| v + increment[i] / self.mass_scale[i])
                .collect(),
        };
        self.admit_state(&state_after)?;
        let mut delta = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in delta.iter_mut().enumerate().take(n) {
            *value = state_after.velocity[i] - state.velocity[i];
        }
        let momentum = action(&self.model.mass, &delta[..n], &mut cancelled)?;
        let relative_residual = self.admit_balance(&[&momentum[..n]], impulse)?;
        let midpoint_velocity = midpoint(&state.velocity, &state_after.velocity);
        let impulse_work_j = pairing(impulse, &midpoint_velocity[..n])?;
        let kinetic_energy_change_j = pairing(&momentum[..n], &midpoint_velocity[..n])?;
        poll(&mut cancelled)?;
        Ok(RotorImpulseResponse {
            state_after,
            relative_residual,
            impulse_work_j,
            kinetic_energy_change_j,
        })
    }

    fn admit_state(&self, state: &RotorTransientState) -> Result<(), RotorError> {
        self.admit_vector(&state.displacement)?;
        self.admit_vector(&state.velocity)
    }

    fn admit_vector(&self, vector: &[f64]) -> Result<(), RotorError> {
        if vector.len() != self.model.n || vector.iter().any(|x| !x.is_finite()) {
            Err(RotorError::Invalid(
                "matching finite original rotor coordinates required",
            ))
        } else {
            Ok(())
        }
    }

    fn admit_balance(&self, terms: &[&[f64]], rhs: &[f64]) -> Result<f64, RotorError> {
        let mut residual_norm = 0.0_f64;
        let mut rhs_norm = 0.0_f64;
        let mut norms = [0.0_f64; 4];
        let term_norms = &mut norms[..terms.len()];
        for i in 0..self.model.n {
            let mut residual = Dd::from_f64(-rhs[i]);
            rhs_norm = rhs_norm.max(finite(rhs[i].abs() / self.mass_scale[i])?);
            for (index, term) in terms.iter().enumerate() {
                residual = residual + Dd::from_f64(term[i]);
                term_norms[index] =
                    term_norms[index].max(finite(term[i].abs() / self.mass_scale[i])?);
            }
            residual_norm =
                residual_norm.max(finite(residual.to_f64().abs() / self.mass_scale[i])?);
        }
        let scale = term_norms.iter().copied().fold(rhs_norm, f64::max);
        let relative = if scale == 0.0 {
            residual_norm
        } else {
            finite(
                (residual_norm / scale)
                    / (rhs_norm / scale + term_norms.iter().map(|x| x / scale).sum::<f64>()),
            )?
        };
        if relative > self.relative_tolerance {
            Err(RotorError::Residual {
                relative,
                limit: self.relative_tolerance,
            })
        } else {
            Ok(relative)
        }
    }
}

fn midpoint(before: &[f64], after: &[f64]) -> [f64; MAX_ROTOR_DOFS] {
    let mut result = [0.0; MAX_ROTOR_DOFS];
    for (i, (&a, &b)) in before.iter().zip(after).enumerate() {
        result[i] = 0.5 * a + 0.5 * b;
    }
    result
}

fn action(
    matrix: &[f64],
    vector: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
    let n = vector.len();
    let mut result = [0.0_f64; MAX_ROTOR_DOFS];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            result[i] = matrix[i * n + j].mul_add(vector[j], result[i]);
        }
        finite(result[i])?;
    }
    Ok(result)
}

fn scaled_action(
    matrix: &[f64],
    vector: &[f64],
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
    let mut result = action(matrix, vector, cancelled)?;
    for value in &mut result[..vector.len()] {
        *value = finite(*value * scale)?;
    }
    Ok(result)
}

fn pairing(force: &[f64], velocity: &[f64]) -> Result<f64, RotorError> {
    let mut sum = Dd::ZERO;
    for (&f, &v) in force.iter().zip(velocity) {
        sum = sum + Dd::from_f64(f) * Dd::from_f64(v);
    }
    finite(sum.to_f64())
}
