//! Finite torsional strain, internal spring reaction and physical shaft losses.

use super::{
    HeldTorque, LieStepError, RotorBank, RotorImpulse, RotorState, RotorStep, ViscousRotorStepper,
    dd, dot_dd, increment, invalid, product, scalar,
};
use fs_ga::Vec3;

/// Prepared rotor/carrier springs, held drive and passive relative-shaft drag.
///
/// Each potential is `K_i*(q_i-q_rest_i)^2/2` on the existing UNWRAPPED
/// relative phase. The spring acts between rotor and carrier, not the world;
/// its reaction preserves total angular momentum. Stiffness, rest phase,
/// viscous coefficient, geometry and inertia stay fixed over this owner.
///
/// A native gyrostat drift evolves `(R,L)` with `(q,p)` fixed. Two implicit
/// midpoint shaft half-stages evolve `(q,p)` at fixed `R` and total body `L`
/// (apart from the declared simultaneous external torque). Eliminating
/// `q_mid=q_before+delta*s_mid/2` gives the original coupled affine solve with
/// `drive_eff=drive-K*(q_before-q_rest)` and `d_eff=d+delta*K/2`.
/// Its ONE native 3-by-3 factor is reused; stepping is O(N), allocation-free.
///
/// This is a second-order split trajectory, not a global symplectic or exact
/// trajectory claim. The native gyrostat energy defect remains explicit.
/// Preparation inherits the affine half-stage resolution bound, including
/// elastic stiffness. Refine oscillation/decay and reject phase increments
/// that the actual stored unwrapped angles cannot adequately represent.
#[derive(Debug, Clone)]
pub struct TorsionalRotorStepper<const N: usize> {
    affine: ViscousRotorStepper<N>,
    stiffness_nm_per_rad: [f64; N],
    rest_angle_rad: [f64; N],
    damping_nms: [f64; N],
    angular_frequency_bound_rad_s: f64,
    has_spring: bool,
}

/// Accepted compliant-shaft motion and its actual mechanical energy account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionalRotorStep<const N: usize> {
    /// Native kinetic/momentum account; shaft work is the NET drive/spring/drag work.
    pub motion: RotorStep<N>,
    /// Signed authored drive work, from the actual shaft midpoint rates, J.
    pub drive_work_j: [f64; N],
    /// Nonnegative physical viscous heat; elastic storage is excluded, J.
    pub dissipated_energy_j: [f64; N],
    /// Polarized spring-potential change between actual stored phases, J.
    pub elastic_energy_change_j: [f64; N],
    /// Delta kinetic + Delta elastic - drive + heat - external - native drift defect.
    /// Measured rounding remains distinct from physical dissipation.
    pub work_energy_residual_j: f64,
}

impl<const N: usize> RotorBank<N> {
    /// Prepare independent torsional springs sharing the COMPLETE carrier recoil.
    ///
    /// Stiffness [N m/rad] and viscous damping [N m s/rad] must be nonnegative;
    /// rest phases are finite, unwrapped angles. Duration is positive. With
    /// `delta=duration/2`, preparation requires
    /// `delta*trace(D*M) + delta^2*trace(K*M)/2 <= 2` (including a small
    /// numerical margin), where `M=diag(1/J)+A^T I^-1 A`.
    /// This bounds the native affine stage; it does not establish adequate
    /// oscillation/decay accuracy. All-zero stiffness delegates exactly to the
    /// existing viscous owner, including its zero-drag held-impulse behavior.
    pub fn torsional_stepper(
        &self,
        stiffness_nm_per_rad: [f64; N],
        rest_angle_rad: [f64; N],
        damping_nms: [f64; N],
        duration_s: f64,
    ) -> Result<TorsionalRotorStepper<N>, LieStepError> {
        scalar(duration_s, "torsional rotor duration")?;
        if duration_s <= 0.0 || 0.5 * duration_s == 0.0 {
            return Err(invalid("positive resolved torsional rotor duration"));
        }
        let mut effective_damping = damping_nms;
        let mut frequency_squared = dd(0.0);
        let mut has_spring = false;
        for i in 0..N {
            let k = scalar(stiffness_nm_per_rad[i], "torsional stiffness")?;
            scalar(rest_angle_rad[i], "unwrapped spring rest phase")?;
            let d = scalar(damping_nms[i], "torsional viscous damping")?;
            if k < 0.0 || d < 0.0 {
                return Err(invalid("nonnegative torsional stiffness and damping"));
            }
            if k != 0.0 {
                has_spring = true;
                let spring_scale = product(
                    product(0.5, duration_s, "shaft half duration")?,
                    k,
                    "torsional midpoint scale",
                )?;
                effective_damping[i] = scalar(
                    d + product(0.5, spring_scale, "torsional midpoint coefficient")?,
                    "combined shaft midpoint coefficient",
                )?;
                let rotor = self.rotors[i];
                let inverse_axis = self.inertia.velocity(rotor.axis_body)?;
                frequency_squared = frequency_squared
                    + dd(k / rotor.axial_inertia_kg_m2)
                    + dd(k) * dot_dd(rotor.axis_body, inverse_axis);
            }
        }
        let frequency_squared = scalar(
            frequency_squared.to_f64() * (1.0 + 256.0 * f64::EPSILON),
            "torsional frequency bound",
        )?;
        let affine = self.viscous_stepper(effective_damping, duration_s)?;
        let angular_frequency_bound_rad_s = scalar(
            fs_math::det::sqrt(frequency_squared),
            "finite nonnegative torsional frequency bound",
        )?;
        Ok(TorsionalRotorStepper {
            affine,
            stiffness_nm_per_rad,
            rest_angle_rad,
            damping_nms,
            angular_frequency_bound_rad_s,
            has_spring,
        })
    }
}

impl<const N: usize> TorsionalRotorStepper<N> {
    /// Prepared positive full interval [s].
    pub const fn duration_s(&self) -> f64 {
        self.affine.duration_s()
    }
    /// Numerical PSD-trace upper estimate of the fastest undamped angular frequency.
    /// Exact for one shaft apart from the floating-point margin; not an interval proof.
    pub const fn angular_frequency_bound_rad_s(&self) -> f64 {
        self.angular_frequency_bound_rad_s
    }
    /// Combined viscous/elastic half-stage resolution number admitted at most 2.
    pub const fn half_stage_resolution_number(&self) -> f64 {
        self.affine.half_stage_relaxation_number()
    }
    /// Total stored spring energy of the actual unwrapped relative phases [J].
    pub fn potential_energy(&self, state: &RotorState<N>) -> Result<f64, LieStepError> {
        self.affine.bank.admit_state(state)?;
        let mut value = dd(0.0);
        for i in 0..N {
            if self.stiffness_nm_per_rad[i] != 0.0 {
                let strain = dd(state.relative_angle_rad[i]) - dd(self.rest_angle_rad[i]);
                value = value + dd(0.5) * dd(self.stiffness_nm_per_rad[i]) * strain * strain;
            }
        }
        scalar(value.to_f64(), "torsional stored potential")
    }

    /// Advance held drive torques, springs, drag and the declared external torque.
    ///
    /// Actual stored phases enter the SECOND shaft stage; the old full-step
    /// phase quadrature is bypassed. Each stage independently admits phase
    /// increments, the original spring/drag law and kinetic-plus-potential work
    /// at local transfer/work scales. Large initial phase/potential energy never
    /// supplies an allowance for an unresolved increment. Any refusal leaves
    /// caller state unchanged; no accepted half-step is published on its own.
    pub fn step(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<TorsionalRotorStep<N>, LieStepError> {
        if !self.has_spring {
            let step = self.affine.step(state, drive_torque_nm, external_torque)?;
            return Ok(TorsionalRotorStep {
                motion: step.motion,
                drive_work_j: step.drive_work_j,
                dissipated_energy_j: step.dissipated_energy_j,
                elastic_energy_change_j: [0.0; N],
                work_energy_residual_j: step.work_energy_residual_j,
            });
        }
        let first = self.stage(
            state,
            drive_torque_nm,
            self.affine
                .external_impulse(state.rotation, external_torque)?,
        )?;
        let mut drifted = first.impulse.state_after;
        drifted.relative_angle_rad = first.phase;
        (drifted.rotation, drifted.total_momentum_body) = self.affine.bank.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.affine
                .bank
                .rotor_momentum(&drifted.axial_momentum_nms)?,
            self.duration_s(),
        )?;
        let second = self.stage(
            &drifted,
            drive_torque_nm,
            self.affine
                .external_impulse(drifted.rotation, external_torque)?,
        )?;
        let motion = self.affine.bank.finish_step(
            state,
            &first.impulse,
            &drifted,
            &second.impulse,
            self.duration_s(),
            Some(second.phase),
        )?;
        let elastic_change = self.potential_changes(state.relative_angle_rad, second.phase)?;
        let mut drive_work = [0.0; N];
        let mut heat = [0.0; N];
        let mut balance = dd(motion.kinetic_energy_change_j)
            - dd(motion.external_work_j)
            - dd(motion.drift_energy_defect_j);
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(first.drive_work[i]) + dd(second.drive_work[i])).to_f64(),
                "torsional drive work",
            )?;
            heat[i] = scalar(
                (dd(first.heat[i]) + dd(second.heat[i])).to_f64(),
                "torsional viscous heat",
            )?;
            balance = balance + dd(elastic_change[i]) - dd(drive_work[i]) + dd(heat[i]);
        }
        Ok(TorsionalRotorStep {
            motion,
            drive_work_j: drive_work,
            dissipated_energy_j: heat,
            elastic_energy_change_j: elastic_change,
            work_energy_residual_j: scalar(balance.to_f64(), "torsional complete work balance")?,
        })
    }

    fn potential_changes(
        &self,
        before: [f64; N],
        after: [f64; N],
    ) -> Result<[f64; N], LieStepError> {
        let mut change = [0.0; N];
        for i in 0..N {
            if self.stiffness_nm_per_rad[i] != 0.0 {
                let mean_strain = (dd(before[i]) - dd(self.rest_angle_rad[i])) * dd(0.5)
                    + (dd(after[i]) - dd(self.rest_angle_rad[i])) * dd(0.5);
                change[i] = scalar(
                    (dd(self.stiffness_nm_per_rad[i])
                        * mean_strain
                        * (dd(after[i]) - dd(before[i])))
                    .to_f64(),
                    "actual spring potential increment",
                )?;
            }
        }
        Ok(change)
    }

    fn stage(
        &self,
        state: &RotorState<N>,
        drive: [f64; N],
        external: Vec3,
    ) -> Result<TorsionalStage<N>, LieStepError> {
        let delta = 0.5 * self.duration_s();
        self.affine.bank.admit_state(state)?;
        let mut effective_drive = drive;
        for i in 0..N {
            scalar(drive[i], "torsional drive torque")?;
            let strain = dd(state.relative_angle_rad[i]) - dd(self.rest_angle_rad[i]);
            scalar(strain.to_f64(), "actual incoming torsional strain")?;
            effective_drive[i] = scalar(
                (dd(drive[i]) - dd(self.stiffness_nm_per_rad[i]) * strain).to_f64(),
                "midpoint elastic drive",
            )?;
        }
        let impulse = self.affine.stage(state, effective_drive, external)?.impulse;
        let mut phase = state.relative_angle_rad;
        let relative_allowance = 4096.0 * f64::EPSILON * (2 * N + 2) as f64;
        for i in 0..N {
            let requested = product(
                delta,
                impulse.mean_relative_rate_rad_s[i],
                "torsional phase increment",
            )?;
            phase[i] = increment(phase[i], requested, "torsional relative phase")?;
            let actual = scalar(
                (dd(phase[i]) - dd(state.relative_angle_rad[i])).to_f64(),
                "actual torsional phase increment",
            )?;
            if (actual - requested).abs() > relative_allowance * actual.abs().max(requested.abs()) {
                return Err(invalid("unresolved torsional phase increment"));
            }
        }
        let elastic_change = self.potential_changes(state.relative_angle_rad, phase)?;
        let mut drive_work = [0.0; N];
        let mut heat = [0.0; N];
        let mut balance = dd(impulse.kinetic_energy_change_j) - dd(impulse.external_work_j);
        let mut scale = impulse
            .kinetic_energy_change_j
            .abs()
            .max(impulse.external_work_j.abs());
        for i in 0..N {
            let rate = impulse.mean_relative_rate_rad_s[i];
            let before_strain = dd(state.relative_angle_rad[i]) - dd(self.rest_angle_rad[i]);
            let after_strain = dd(phase[i]) - dd(self.rest_angle_rad[i]);
            scalar(after_strain.to_f64(), "actual returned torsional strain")?;
            let mean_strain = before_strain * dd(0.5) + after_strain * dd(0.5);
            let authored = dd(delta) * dd(drive[i]);
            let spring = dd(delta) * dd(self.stiffness_nm_per_rad[i]) * mean_strain;
            let drag = dd(delta) * dd(self.damping_nms[i]) * dd(rate);
            let expected = scalar(
                (authored - spring - drag).to_f64(),
                "original midpoint shaft law",
            )?;
            let actual = impulse.axial_momentum_change_nms[i];
            let residual = scalar(
                (dd(actual) - authored + spring + drag).to_f64(),
                "torsional shaft law residual",
            )?;
            let law_scale = actual
                .abs()
                .max(authored.to_f64().abs())
                .max(spring.to_f64().abs())
                .max(drag.to_f64().abs());
            if (expected != 0.0 && actual == 0.0) || residual.abs() > relative_allowance * law_scale
            {
                return Err(invalid("unresolved implicit torsional effort"));
            }
            drive_work[i] = scalar((authored * dd(rate)).to_f64(), "torsional shaft drive work")?;
            heat[i] = scalar(
                (drag * dd(rate)).to_f64(),
                "physical torsional viscous heat",
            )?;
            balance = balance + dd(elastic_change[i]) - dd(drive_work[i]) + dd(heat[i]);
            scale = scale
                .max(elastic_change[i].abs())
                .max(drive_work[i].abs())
                .max(heat[i]);
        }
        if scalar(balance.to_f64(), "torsional stage work balance")?.abs()
            > relative_allowance * scale
        {
            return Err(invalid("unresolved torsional stage potential work"));
        }
        Ok(TorsionalStage {
            impulse,
            phase,
            drive_work,
            heat,
        })
    }
}

struct TorsionalStage<const N: usize> {
    impulse: RotorImpulse<N>,
    phase: [f64; N],
    drive_work: [f64; N],
    heat: [f64; N],
}
