//! Coupled relative-shaft drag and held drive effort on the native rotor bank.

use super::{
    HeldTorque, LieStepError, RotorBank, RotorImpulse, RotorState, RotorStep, components, dd,
    dot_dd, invalid, product, scalar, vector,
};
use fs_ga::{So3, Vec3};
use fs_la::factor::{Cholesky, cholesky};

/// Prepared midpoint shaft stages for `torque_i = drive_i - damping_i * relative_rate_i`.
///
/// A positive interval and nonnegative viscous coefficients are held fixed.
/// Diagonal rotor self-mobility is eliminated before ONE native 3-by-3 Cholesky
/// factorization. Each step is O(N), allocation-free, with both motor reaction
/// and full-tensor cross-axis recoil; shafts are not solved independently.
///
/// Midpoint is passive but not L-stable. Preparation refuses a half-stage
/// relaxation number above 2, before any unresolved Cayley sign reversal.
/// The bound is conservative for multiple rotors; smaller steps and refinement
/// are still necessary for accurate decay and gyroscopic motion. This fixed
/// affine law supplies neither electrical dynamics nor arbitrary feedback.
#[derive(Debug, Clone)]
pub struct ViscousRotorStepper<const N: usize> {
    pub(super) bank: RotorBank<N>,
    duration_s: f64,
    half_s: f64,
    damping_nms: [f64; N],
    drag_impulse_scale: [f64; N],
    drive_scale: [f64; N],
    recoil_weight: [f64; N],
    factor: Option<Cholesky>,
    factor_scale: f64,
    relaxation_rate_bound_per_s: f64,
    half_stage_relaxation_number: f64,
}

/// Actual trajectory and separate drive work / passive viscous heat.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousRotorStep<const N: usize> {
    /// Native trajectory, net shaft work, external work and separate drift defect.
    pub motion: RotorStep<N>,
    /// Signed drive work at actual midpoint relative rates, J. All-zero drag
    /// delegates to the native held-torque work of its actual rounded impulses.
    pub drive_work_j: [f64; N],
    /// Nonnegative physical viscous dissipation, sum(dt_stage*d*s_mid^2), J.
    pub dissipated_energy_j: [f64; N],
    /// Kinetic change - drive work + heat - external work - native drift defect.
    /// Includes measured finite-precision constitutive/work rounding; not heat.
    pub work_energy_residual_j: f64,
}

impl<const N: usize> RotorBank<N> {
    /// Upper estimate for the fastest viscous relaxation rate, in 1/s.
    ///
    /// The shared mobility is `M = diag(1/J) + A^T I^-1 A`. The trace of
    /// `sqrt(D) M sqrt(D)` bounds its largest eigenvalue and is exact for one
    /// shaft. A small floating-point margin is included; this is a numerical
    /// admission bound, not an interval spectral certificate. Zero drag gives 0.
    pub fn viscous_relaxation_rate_bound(
        &self,
        damping_nms: [f64; N],
    ) -> Result<f64, LieStepError> {
        let mut trace = dd(0.0);
        for (rotor, damping) in self.rotors.iter().zip(damping_nms) {
            scalar(damping, "relative-shaft viscous damping")?;
            if damping < 0.0 {
                return Err(invalid("nonnegative relative-shaft viscous damping"));
            }
            if damping != 0.0 {
                let inverse_axis = self.inertia.velocity(rotor.axis_body)?;
                trace = trace
                    + dd(damping / rotor.axial_inertia_kg_m2)
                    + dd(damping) * dot_dd(rotor.axis_body, inverse_axis);
            }
        }
        let bound = scalar(
            trace.to_f64() * (1.0 + 256.0 * f64::EPSILON),
            "viscous relaxation rate bound",
        )?;
        if bound < 0.0 {
            return Err(invalid("positive viscous relaxation mobility"));
        }
        Ok(bound)
    }

    /// Prepare passive relative-shaft drag for one positive held time interval.
    ///
    /// Refuses if `duration/2 * viscous_relaxation_rate_bound(damping) > 2`.
    /// Query that bound before selecting a timestep; this admission avoids
    /// unresolved midpoint sign flips and does not replace convergence checks.
    /// Reprepare when the interval, drag, rotor inventory or inertia changes.
    pub fn viscous_stepper(
        &self,
        damping_nms: [f64; N],
        duration_s: f64,
    ) -> Result<ViscousRotorStepper<N>, LieStepError> {
        scalar(duration_s, "viscous rotor duration")?;
        if duration_s <= 0.0 || 0.5 * duration_s == 0.0 {
            return Err(invalid("positive resolved viscous rotor duration"));
        }
        let half_s = 0.5 * duration_s;
        let rate = self.viscous_relaxation_rate_bound(damping_nms)?;
        let relaxation = scalar(half_s * rate, "half-stage viscous relaxation number")?;
        if relaxation > 2.0 {
            return Err(invalid(
                "viscous rotor interval does not resolve relaxation",
            ));
        }
        let mut drag_scale = [0.0; N];
        let mut drive_scale = [half_s; N];
        let mut weight = [0.0; N];
        let mut matrix = self.inertia.matrix().m.map(dd);
        let mut any_drag = false;
        for i in 0..N {
            let damping = damping_nms[i];
            drag_scale[i] = product(half_s, damping, "viscous shaft impulse scale")?;
            if damping == 0.0 {
                continue;
            }
            any_drag = true;
            let a = product(0.5, drag_scale[i], "viscous midpoint coefficient")?;
            let j = self.rotors[i].axial_inertia_kg_m2;
            // beta = a*J/(J+a), gamma = half_s*J/(J+a). No large
            // relaxation number or product a*J is formed by this elimination.
            weight[i] = if a <= j {
                a / (1.0 + a / j)
            } else {
                j / (1.0 + j / a)
            };
            drive_scale[i] = product(2.0, weight[i] / damping, "viscous drive scale")?;
            let axis = components(self.rotors[i].axis_body);
            for row in 0..3 {
                for col in row..3 {
                    matrix[3 * row + col] =
                        matrix[3 * row + col] + dd(weight[i]) * dd(axis[row]) * dd(axis[col]);
                    matrix[3 * col + row] = matrix[3 * row + col];
                }
            }
        }
        let mut values = matrix.map(|value| value.to_f64());
        let scale = values.iter().map(|value| value.abs()).fold(0.0, f64::max);
        scalar(scale, "viscous carrier factor scale")?;
        let factor = if any_drag {
            for value in &mut values {
                *value = scalar(*value / scale, "viscous carrier factor")?;
            }
            Some(cholesky(&values, 3).map_err(|_| invalid("viscous carrier factor not SPD"))?)
        } else {
            None
        };
        Ok(ViscousRotorStepper {
            bank: *self,
            duration_s,
            half_s,
            damping_nms,
            drag_impulse_scale: drag_scale,
            drive_scale,
            recoil_weight: weight,
            factor,
            factor_scale: scale,
            relaxation_rate_bound_per_s: rate,
            half_stage_relaxation_number: relaxation,
        })
    }
}

impl<const N: usize> ViscousRotorStepper<N> {
    /// Positive held full-step interval, s.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Prepared upper estimate for the fastest relative-shaft relaxation rate.
    pub const fn relaxation_rate_bound_per_s(&self) -> f64 {
        self.relaxation_rate_bound_per_s
    }
    /// Prepared half-step times the upper relaxation-rate estimate (at most 2).
    pub const fn half_stage_relaxation_number(&self) -> f64 {
        self.half_stage_relaxation_number
    }
    /// Nonnegative per-shaft viscous coefficients, N m s / rad.
    pub const fn damping_nms(&self) -> &[f64; N] {
        &self.damping_nms
    }

    /// Advance held drive torque and a declared body- or space-held external torque.
    ///
    /// The two implicit midpoint shaft stages surround one native gyrostat
    /// drift; the actual accepted endpoint rates also integrate relative phase.
    /// Every stage uses the existing simultaneous impulse/recoil admission.
    /// The constitutive impulse equation and full stage work are independently
    /// checked after native rounding, using local transfer/work scales without
    /// a background-energy allowance. Physical heat is `dt*d*s_mid^2`, never a
    /// drift or solve defect.
    /// With zero drag the returned `motion` equals the existing held-torque step
    /// bit for bit. All inputs are unchanged if either stage or phase refuses.
    pub fn step(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<ViscousRotorStep<N>, LieStepError> {
        if self.factor.is_none() {
            let motion =
                self.bank
                    .step(state, drive_torque_nm, external_torque, self.duration_s)?;
            return Ok(ViscousRotorStep {
                drive_work_j: motion.shaft_work_j,
                dissipated_energy_j: [0.0; N],
                work_energy_residual_j: motion.work_energy_residual_j,
                motion,
            });
        }
        let first = self.stage(
            state,
            drive_torque_nm,
            self.external_impulse(state.rotation, external_torque)?,
        )?;
        let mut drifted = first.impulse.state_after;
        (drifted.rotation, drifted.total_momentum_body) = self.bank.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.bank.rotor_momentum(&drifted.axial_momentum_nms)?,
            self.duration_s,
        )?;
        let second = self.stage(
            &drifted,
            drive_torque_nm,
            self.external_impulse(drifted.rotation, external_torque)?,
        )?;
        let motion = self.bank.finish_step(
            state,
            &first.impulse,
            &drifted,
            &second.impulse,
            self.duration_s,
            None,
        )?;
        let mut drive_work = [0.0; N];
        let mut heat = [0.0; N];
        let mut residual = dd(motion.kinetic_energy_change_j)
            - dd(motion.external_work_j)
            - dd(motion.drift_energy_defect_j);
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(first.drive_work[i]) + dd(second.drive_work[i])).to_f64(),
                "viscous rotor drive work",
            )?;
            heat[i] = scalar(
                (dd(first.heat[i]) + dd(second.heat[i])).to_f64(),
                "viscous rotor heat",
            )?;
            residual = residual - dd(drive_work[i]) + dd(heat[i]);
        }
        Ok(ViscousRotorStep {
            motion,
            drive_work_j: drive_work,
            dissipated_energy_j: heat,
            work_energy_residual_j: scalar(residual.to_f64(), "viscous rotor work balance")?,
        })
    }

    pub(super) fn external_impulse(
        &self,
        rotation: So3,
        torque: HeldTorque,
    ) -> Result<Vec3, LieStepError> {
        let body = match torque {
            HeldTorque::Body(value) => vector(value, "external body torque")?,
            HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
        };
        Ok(Vec3::new(
            product(body.x, self.half_s, "external torque impulse")?,
            product(body.y, self.half_s, "external torque impulse")?,
            product(body.z, self.half_s, "external torque impulse")?,
        ))
    }

    pub(super) fn stage(
        &self,
        state: &RotorState<N>,
        drive: [f64; N],
        external: Vec3,
    ) -> Result<ShaftStage<N>, LieStepError> {
        let mut impulses = [0.0; N];
        for i in 0..N {
            impulses[i] = product(drive[i], self.half_s, "held shaft drive impulse")?;
        }
        if let Some(factor) = &self.factor {
            let initial_rate = self.bank.relative_rates(state)?;
            // Use the actual rounded external transfer that apply_impulse will
            // accept, so drag and externally driven carrier motion are simultaneous.
            let next_total =
                super::super::super::torque_kick(state.total_momentum_body, external, 1.0)?;
            let mut rhs = components(next_total - state.total_momentum_body).map(dd);
            for i in 0..N {
                let rotor = self.bank.rotors[i];
                impulses[i] = scalar(
                    (dd(self.drive_scale[i]) * dd(drive[i])
                        - dd(2.0) * dd(self.recoil_weight[i]) * dd(initial_rate[i]))
                    .to_f64(),
                    "coupled viscous shaft impulse",
                )?;
                for (value, axis) in rhs.iter_mut().zip(components(rotor.axis_body)) {
                    *value = *value - dd(axis) * dd(impulses[i]);
                }
            }
            let mut solve = rhs.map(|value| value.to_f64() / self.factor_scale);
            factor.solve(&mut solve);
            let carrier_velocity_change = vector(
                Vec3::new(solve[0], solve[1], solve[2]),
                "coupled viscous carrier solve",
            )?;
            for i in 0..N {
                impulses[i] = scalar(
                    (dd(impulses[i])
                        + dd(self.recoil_weight[i])
                            * dot_dd(self.bank.rotors[i].axis_body, carrier_velocity_change))
                    .to_f64(),
                    "coupled viscous shaft correction",
                )?;
            }
        }
        let impulse = self.bank.apply_impulse(state, impulses, external)?;
        let mut drive_work = [0.0; N];
        let mut heat = [0.0; N];
        let mut work_defect = dd(impulse.kinetic_energy_change_j) - dd(impulse.external_work_j);
        let mut work_scale = impulse
            .kinetic_energy_change_j
            .abs()
            .max(impulse.external_work_j.abs());
        for i in 0..N {
            let rate = impulse.mean_relative_rate_rad_s[i];
            let drive_impulse = dd(self.half_s) * dd(drive[i]);
            let drag_impulse = dd(self.drag_impulse_scale[i]) * dd(rate);
            let expected = scalar(
                (drive_impulse - drag_impulse).to_f64(),
                "midpoint shaft law",
            )?;
            let actual = impulse.axial_momentum_change_nms[i];
            let residual = scalar(
                (dd(actual) - drive_impulse + drag_impulse).to_f64(),
                "midpoint shaft law residual",
            )?;
            let scale = actual
                .abs()
                .max(drive_impulse.to_f64().abs())
                .max(drag_impulse.to_f64().abs());
            if (expected != 0.0 && actual == 0.0) || residual.abs() > 4096.0 * f64::EPSILON * scale
            {
                return Err(invalid("unresolved implicit relative-shaft effort"));
            }
            drive_work[i] = scalar((drive_impulse * dd(rate)).to_f64(), "midpoint drive work")?;
            heat[i] = scalar(
                (drag_impulse * dd(rate)).to_f64(),
                "passive relative-shaft heat",
            )?;
            work_defect = work_defect - dd(drive_work[i]) + dd(heat[i]);
            work_scale = work_scale.max(drive_work[i].abs()).max(heat[i]);
        }
        // A huge preexisting carrier momentum can admit a partially rounded
        // recoil while leaving the shaft law accurate. Reject its actual work
        // defect instead of disguising it as dissipation or merely reporting it.
        let work_defect = scalar(work_defect.to_f64(), "implicit shaft stage work balance")?;
        let allowance = 4096.0 * f64::EPSILON * (2 * N + 2) as f64 * work_scale;
        if work_defect.abs() > allowance {
            return Err(invalid("unresolved implicit shaft stage work"));
        }
        Ok(ShaftStage {
            impulse,
            drive_work,
            heat,
        })
    }
}

pub(super) struct ShaftStage<const N: usize> {
    pub(super) impulse: RotorImpulse<N>,
    drive_work: [f64; N],
    heat: [f64; N],
}
