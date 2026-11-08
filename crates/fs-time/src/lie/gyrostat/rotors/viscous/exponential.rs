//! Exact frozen-attitude shaft exchange with the original carrier recoil.
//!
//! For the positively damped shaft subset, M=diag(1/J)+A^T I^-1 A is SPD.
//! Factoring `I_ref M=L L^T` and diagonalizing `L^T D L` gives independent
//! relaxation coordinates without inverting D. Undamped shafts keep their
//! exact drive impulses and force every damped channel through the carrier.
//! The existing native gyrostat still owns attitude and precession drift.

// Original shaft/mode indices make the dense transforms explicit. Every
// integer-to-float count is bounded by the admitted 64-rotor inventory.
#![allow(clippy::needless_range_loop, clippy::cast_precision_loss)]

use super::super::{MAX_AXIAL_ROTORS, increment};
use super::{
    HeldTorque, LieStepError, RotorBank, RotorImpulse, RotorState, components, dd, dot_dd, invalid,
    product, scalar, vector,
};
use crate::lie::{compensated, viscous::exponential::Response};
use fs_ga::Vec3;
use fs_la::{eigen::jacobi_eigh, factor::cholesky};
use fs_math::dd::Dd;

/// Prepared exponential relative-shaft damping on the original rotor bank.
///
/// All N original shafts and their full carrier recoil remain coupled. Only
/// K shafts with strictly positive damping enter the dense preparation; zero
/// damping is structural, never an approximate eigensolver nullspace. The
/// bank's explicit [`MAX_AXIAL_ROTORS`] bound (64) bounds K, preparation's
/// O(K^3) work/O(K^2) storage and each step's O(N+K^2) work. Cached shaft stages
/// allocate nothing and evaluate no factorization, eigensolve or exponential.
///
/// Exactness applies to the linear shaft stages with frozen attitude. Their
/// symmetric composition with native gyrostat drift remains second order;
/// refine the interval for tumbling and spatially held torque reorientation.
/// Aggregate physical slip heat is integrated independently; individual
/// shaft heat and continuous shaft-transfer work are not supplied.
#[derive(Debug, Clone)]
pub struct ExponentialViscousRotorStepper<const N: usize> {
    bank: RotorBank<N>,
    damping_nms: [f64; N],
    duration_s: f64,
    half_s: f64,
    reference_moment: f64,
    active: [usize; N],
    count: usize,
    relative_modes: Vec<f64>,
    dual_modes: Vec<f64>,
    response: Vec<Response>,
}

/// Actual motion and physical work of the two exponential shaft stages.
///
/// Unlike instantaneous kick work, drive and external work integrate their
/// actual continuous frozen-attitude responses. Heat is the total of the
/// original `sum_i damping_i * relative_rate_i^2` integral. No per-shaft
/// heat allocation or heat inferred from a kinetic-energy gap is published.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExponentialViscousRotorStep<const N: usize> {
    /// Complete accepted state for the next native rotor step.
    pub state_after: RotorState<N>,
    /// Actual accepted rotor momentum changes [N m s].
    pub axial_momentum_change_nms: [f64; N],
    /// Actual external angular impulse in inertial coordinates [N m s].
    pub external_angular_impulse_world: Vec3,
    /// Signed per-drive work, torque times integrated relative shaft angle [J].
    pub drive_work_j: [f64; N],
    /// Nonnegative total physical relative-shaft heat [J].
    pub dissipated_energy_j: f64,
    /// Total drive work minus the independently integrated total heat [J].
    pub net_shaft_work_j: f64,
    /// Integrated held external torque work on the carrier [J].
    pub external_work_j: f64,
    /// Actual full carrier-plus-rotor kinetic-energy increment [J].
    pub kinetic_energy_change_j: f64,
    /// Native gyrostat drift's numerical kinetic-energy defect [J], never heat.
    pub drift_energy_defect_j: f64,
    /// Kinetic change - net shaft work - external work - drift defect [J].
    pub work_energy_residual_j: f64,
    /// Actual spatial momentum change minus external impulse [N m s].
    pub angular_momentum_residual_world: Vec3,
}

impl<const N: usize> RotorBank<N> {
    /// Prepare exact linear shaft damping, including stiff finite relaxation.
    ///
    /// Nonnegative damping and a positive resolved duration are required.
    /// There is no midpoint relaxation-number cutoff. The original N<=64
    /// inventory bounds all work and memory; dense preparation acts only on
    /// the positive-damping subset. Original full-tensor mobility, source
    /// impulses, phase increments and physical work must remain numerically
    /// resolved. Zero damping delegates to the native held-torque step.
    ///
    /// # Errors
    /// Invalid parameters, an unresolved mobility/eigendecomposition, or
    /// unrepresentable response coefficients. No partial stepper is returned.
    pub fn viscous_exponential_stepper(
        &self,
        damping_nms: [f64; N],
        duration_s: f64,
    ) -> Result<ExponentialViscousRotorStepper<N>, LieStepError> {
        scalar(duration_s, "exponential rotor duration")?;
        let half_s = 0.5 * duration_s;
        if N > MAX_AXIAL_ROTORS || duration_s <= 0.0 || half_s == 0.0 {
            return Err(invalid("bounded positive exponential rotor duration"));
        }
        let mut active = [0; N];
        let mut count = 0;
        for (i, &damping) in damping_nms.iter().enumerate() {
            scalar(damping, "exponential shaft damping")?;
            if damping < 0.0 {
                return Err(invalid("nonnegative exponential shaft damping"));
            }
            if damping > 0.0 {
                active[count] = i;
                count += 1;
            }
        }
        let mut prepared = ExponentialViscousRotorStepper {
            bank: *self,
            damping_nms,
            duration_s,
            half_s,
            reference_moment: 1.0,
            active,
            count,
            relative_modes: Vec::new(),
            dual_modes: Vec::new(),
            response: Vec::new(),
        };
        if count != 0 {
            prepared.prepare()?;
        }
        Ok(prepared)
    }
}

impl<const N: usize> ExponentialViscousRotorStepper<N> {
    /// Fixed positive full-step duration [s].
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Original per-shaft viscous coefficients [N m s / rad].
    #[must_use]
    pub const fn damping_nms(&self) -> &[f64; N] {
        &self.damping_nms
    }

    /// Number K of coupled positive-damping relaxation modes prepared.
    #[must_use]
    pub const fn damped_rotor_count(&self) -> usize {
        self.count
    }

    #[allow(clippy::too_many_lines)] // One bounded mobility/spectrum admission with reciprocal basis checks.
    fn prepare(&mut self) -> Result<(), LieStepError> {
        let k = self.count;
        let moments = components(self.bank.inertia.principal_moments());
        self.reference_moment = self.active[..k]
            .iter()
            .map(|&i| self.bank.rotors[i].axial_inertia_kg_m2)
            .chain(moments)
            .fold(f64::INFINITY, f64::min);
        let root_reference = self.reference_moment.sqrt();
        let to_principal = self.bank.inertia.principal_to_body().inverse();
        let mut axes = [[0.0; 3]; N];
        for (row, &i) in self.active[..k].iter().enumerate() {
            let principal = components(to_principal.rotate(self.bank.rotors[i].axis_body)?);
            axes[row] =
                core::array::from_fn(|j| principal[j] * (root_reference / moments[j].sqrt()));
        }
        let mut mobility = vec![0.0; k * k];
        for i in 0..k {
            for j in i..k {
                let mut value = Dd::ZERO;
                for a in 0..3 {
                    value = value + dd(axes[i][a]) * dd(axes[j][a]);
                }
                if i == j {
                    value = value
                        + dd(self.reference_moment
                            / self.bank.rotors[self.active[i]].axial_inertia_kg_m2);
                }
                let value = scalar(value.to_f64(), "exponential rotor mobility")?;
                mobility[i * k + j] = value;
                mobility[j * k + i] = value;
            }
        }
        let factor =
            cholesky(&mobility, k).map_err(|_| invalid("resolved exponential rotor mobility"))?;
        let damping_scale = self.damping_nms.iter().copied().fold(0.0_f64, f64::max);
        let mut normalized = vec![0.0; k * k];
        for i in 0..k {
            for j in i..k {
                let mut value = Dd::ZERO;
                for row in 0..k {
                    value = value
                        + dd(self.damping_nms[self.active[row]] / damping_scale)
                            * dd(factor.l(row, i))
                            * dd(factor.l(row, j));
                }
                normalized[i * k + j] = value.to_f64();
                normalized[j * k + i] = value.to_f64();
            }
        }
        let scale = normalized.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if !scale.is_finite() || scale == 0.0 {
            return Err(invalid("resolved exponential rotor damping"));
        }
        for value in &mut normalized {
            *value /= scale;
        }
        let (_, modes) = jacobi_eigh(&normalized, k);
        self.relative_modes = vec![0.0; k * k];
        self.dual_modes = vec![0.0; k * k];
        self.response = vec![Response::default(); k];
        for mode in 0..k {
            let mut relative = [0.0; N];
            let mut normalized_rate = Dd::ZERO;
            let mut physical_rate = Dd::ZERO;
            for row in 0..k {
                let mut value = Dd::ZERO;
                for col in 0..k {
                    value = value + dd(factor.l(row, col)) * dd(modes[col * k + mode]);
                }
                relative[row] = scalar(value.to_f64(), "exponential rotor mode")?;
                let square = dd(relative[row]) * dd(relative[row]);
                normalized_rate = normalized_rate
                    + dd(self.damping_nms[self.active[row]] / damping_scale) * square;
                physical_rate = physical_rate + dd(self.damping_nms[self.active[row]]) * square;
            }
            let eigenvalue = normalized_rate.to_f64() / scale;
            if !eigenvalue.is_finite() || eigenvalue <= 0.0 {
                return Err(invalid("resolved positive exponential rotor rate"));
            }
            for row in 0..k {
                let mut applied = Dd::ZERO;
                let mut orthogonal = Dd::ZERO;
                for col in 0..k {
                    applied = applied + dd(normalized[row * k + col]) * dd(modes[col * k + mode]);
                    orthogonal = orthogonal + dd(modes[col * k + row]) * dd(modes[col * k + mode]);
                }
                let target = if row == mode { 1.0 } else { 0.0 };
                let allowance = 4096.0 * f64::EPSILON * k as f64;
                if (applied - dd(eigenvalue) * dd(modes[row * k + mode]))
                    .to_f64()
                    .abs()
                    > allowance
                    || (orthogonal.to_f64() - target).abs() > allowance
                {
                    return Err(invalid("resolved orthogonal exponential rotor spectrum"));
                }
            }
            let number =
                compensated::divide(dd(self.half_s) * physical_rate, dd(self.reference_moment))?
                    .to_f64();
            if !number.is_finite() || number <= 0.0 {
                return Err(invalid("resolved exponential rotor relaxation number"));
            }
            let response = Response::prepare(number, 1.0)?;
            self.response[mode] = response;
            let mut dual = relative;
            factor.solve(&mut dual[..k]);
            for row in 0..k {
                self.relative_modes[row * k + mode] = relative[row];
                self.dual_modes[row * k + mode] = scalar(dual[row], "exponential rotor dual mode")?;
            }
        }
        for i in 0..k {
            for j in 0..k {
                let mut value = Dd::ZERO;
                for row in 0..k {
                    value = value
                        + dd(self.dual_modes[row * k + i]) * dd(self.relative_modes[row * k + j]);
                }
                if (value.to_f64() - if i == j { 1.0 } else { 0.0 }).abs()
                    > 8192.0 * f64::EPSILON * k as f64
                {
                    return Err(invalid("resolved exponential rotor duality"));
                }
            }
        }
        Ok(())
    }

    /// Advance held shaft drives and body- or space-held external torque.
    ///
    /// The two exact affine shaft stages surround the existing native gyrostat
    /// drift. Their integrated rates update every original unwrapped phase.
    /// Actual accepted impulses/recoil and the full physical stage work must
    /// resolve at local transfer/work scales, without a background-energy
    /// allowance. Inputs remain unchanged on any refusal.
    ///
    /// # Errors
    /// Invalid/unresolved source, state, phase, recoil or physical work, or a
    /// refusal from the native gyrostat. No partial trajectory is returned.
    pub fn step(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<ExponentialViscousRotorStep<N>, LieStepError> {
        if self.count == 0 {
            let native =
                self.bank
                    .step(state, drive_torque_nm, external_torque, self.duration_s)?;
            let net = native
                .shaft_work_j
                .iter()
                .fold(Dd::ZERO, |sum, &v| sum + dd(v));
            return Ok(ExponentialViscousRotorStep {
                state_after: native.state_after,
                axial_momentum_change_nms: native.axial_momentum_change_nms,
                external_angular_impulse_world: native.external_angular_impulse_world,
                drive_work_j: native.shaft_work_j,
                dissipated_energy_j: 0.0,
                net_shaft_work_j: scalar(net.to_f64(), "exponential undamped shaft work")?,
                external_work_j: native.external_work_j,
                kinetic_energy_change_j: native.kinetic_energy_change_j,
                drift_energy_defect_j: native.drift_energy_defect_j,
                work_energy_residual_j: native.work_energy_residual_j,
                angular_momentum_residual_world: native.angular_momentum_residual_world,
            });
        }
        let first = self.stage(state, drive_torque_nm, external_torque)?;
        let mut drifted = first.impulse.state_after;
        (drifted.rotation, drifted.total_momentum_body) = self.bank.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.bank.rotor_momentum(&drifted.axial_momentum_nms)?,
            self.duration_s,
        )?;
        let second = self.stage(&drifted, drive_torque_nm, external_torque)?;
        let mut phase = [0.0; N];
        let mut drive_work = [0.0; N];
        let mut net = Dd::ZERO;
        for i in 0..N {
            let delta = scalar(
                (first.phase[i] + second.phase[i]).to_f64(),
                "exponential rotor phase increment",
            )?;
            phase[i] = increment(
                state.relative_angle_rad[i],
                delta,
                "exponential rotor phase",
            )?;
            drive_work[i] = scalar(
                (first.drive_work[i] + second.drive_work[i]).to_f64(),
                "exponential rotor drive work",
            )?;
            net = net + first.drive_work[i] + second.drive_work[i];
        }
        let native = self.bank.finish_step(
            state,
            &first.impulse,
            &drifted,
            &second.impulse,
            self.duration_s,
            Some(phase),
        )?;
        let heat = first.heat + second.heat;
        let external_work = first.external_work + second.external_work;
        net = net - heat;
        let residual = dd(native.kinetic_energy_change_j)
            - net
            - external_work
            - dd(native.drift_energy_defect_j);
        Ok(ExponentialViscousRotorStep {
            state_after: native.state_after,
            axial_momentum_change_nms: native.axial_momentum_change_nms,
            external_angular_impulse_world: native.external_angular_impulse_world,
            drive_work_j: drive_work,
            dissipated_energy_j: scalar(heat.to_f64(), "exponential rotor heat")?,
            net_shaft_work_j: scalar(net.to_f64(), "exponential net shaft work")?,
            external_work_j: scalar(external_work.to_f64(), "exponential external work")?,
            kinetic_energy_change_j: native.kinetic_energy_change_j,
            drift_energy_defect_j: native.drift_energy_defect_j,
            work_energy_residual_j: scalar(residual.to_f64(), "exponential rotor work balance")?,
            angular_momentum_residual_world: native.angular_momentum_residual_world,
        })
    }

    #[allow(clippy::too_many_lines)] // One atomic constitutive/impulse/work transaction on fixed stack storage.
    fn stage(
        &self,
        state: &RotorState<N>,
        drive: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<Stage<N>, LieStepError> {
        let k = self.count;
        let initial_rate = self.bank.relative_rates(state)?;
        let omega = self.bank.angular_velocity(state)?;
        let torque = match external_torque {
            HeldTorque::Body(value) => vector(value, "exponential external body torque")?,
            HeldTorque::Space(value) => state.rotation.inverse().rotate(value)?,
        };
        let requested_external = Vec3::new(
            product(torque.x, self.half_s, "exponential external impulse")?,
            product(torque.y, self.half_s, "exponential external impulse")?,
            product(torque.z, self.half_s, "exponential external impulse")?,
        );
        let next_total =
            crate::lie::torque_kick(state.total_momentum_body, requested_external, 1.0)?;
        let external = next_total - state.total_momentum_body;
        let mut impulses = [0.0; N];
        let mut mean_increment = [Dd::ZERO; N];
        let mut known_carrier = components(external).map(dd);
        for i in 0..N {
            scalar(drive[i], "exponential shaft drive")?;
            if self.damping_nms[i] == 0.0 {
                impulses[i] = product(self.half_s, drive[i], "undamped shaft impulse")?;
                mean_increment[i] = dd(0.5) * dd(impulses[i]);
                for (value, axis) in known_carrier
                    .iter_mut()
                    .zip(components(self.bank.rotors[i].axis_body))
                {
                    *value = *value - dd(axis) * dd(impulses[i]);
                }
            }
        }
        let known_velocity = self
            .bank
            .inertia
            .velocity(as_vector(known_carrier, "exponential carrier forcing")?)?;
        let mut known_rates = [Dd::ZERO; N];
        for row in 0..k {
            known_rates[row] =
                -dot_dd(self.bank.rotors[self.active[row]].axis_body, known_velocity);
        }
        let mut change = [Dd::ZERO; N];
        let mut mean_change = [Dd::ZERO; N];
        let mut mean = [Dd::ZERO; N];
        let mut heat = Dd::ZERO;
        let root_reference = dd(self.reference_moment.sqrt());
        for mode in 0..k {
            let mut initial = Dd::ZERO;
            let mut drive_mode = Dd::ZERO;
            let mut initial_effort = Dd::ZERO;
            let mut known = Dd::ZERO;
            for row in 0..k {
                initial = initial
                    + dd(self.dual_modes[row * k + mode]) * dd(initial_rate[self.active[row]]);
                drive_mode = drive_mode
                    + dd(self.relative_modes[row * k + mode]) * dd(drive[self.active[row]]);
                let i = self.active[row];
                initial_effort = initial_effort
                    + dd(self.relative_modes[row * k + mode])
                        * (dd(drive[i]) - dd(self.damping_nms[i]) * dd(initial_rate[i]));
                known = known + dd(self.dual_modes[row * k + mode]) * known_rates[row];
            }
            let source =
                compensated::divide(dd(self.half_s) * drive_mode, dd(self.reference_moment))?
                    + known;
            let r = self.response[mode];
            // Applying phi_1 to the original initial torque residual preserves
            // exact driven equilibrium, even after exponential saturation.
            let residual =
                compensated::divide(dd(self.half_s) * initial_effort, dd(self.reference_moment))?
                    + known;
            let delta = dd(r.response) * residual;
            change[mode] = delta - known;
            mean_change[mode] = dd(r.integrated_response) * residual - dd(0.5) * known;
            mean[mode] = dd(r.response) * initial + dd(r.integrated_response) * source;
            let heat_mean =
                root_reference * (dd(r.heat_mean[0]) * initial + dd(r.heat_mean[1]) * source);
            let deviation = root_reference
                * (dd(r.heat_deviation[0]) * initial + dd(r.heat_deviation[1]) * source);
            heat = heat + heat_mean * heat_mean + deviation * deviation;
        }
        for row in 0..k {
            let i = self.active[row];
            let mut delta = Dd::ZERO;
            for mode in 0..k {
                let dual = dd(self.reference_moment) * dd(self.dual_modes[row * k + mode]);
                delta = delta + dual * change[mode];
                mean_increment[i] = mean_increment[i] + dual * mean_change[mode];
            }
            impulses[i] = scalar(delta.to_f64(), "exponential shaft impulse")?;
        }
        let impulse = self
            .bank
            .apply_impulse(state, impulses, requested_external)?;
        let mut carrier_mean_change = components(external).map(|x| dd(0.5) * dd(x));
        for i in 0..N {
            for (value, axis) in carrier_mean_change
                .iter_mut()
                .zip(components(self.bank.rotors[i].axis_body))
            {
                *value = *value - dd(axis) * mean_increment[i];
            }
        }
        let mean_velocity_change = self.bank.inertia.velocity(as_vector(
            carrier_mean_change,
            "exponential mean carrier impulse",
        )?)?;
        let external_work = dot_dd(external, omega) + dot_dd(external, mean_velocity_change);
        let mut mean_rates = [Dd::ZERO; N];
        // Damped-channel integrals come directly from the exponential response.
        // Reconstructing them as initial rate plus a nearly opposite recoil
        // would erase the small, finite angle of a saturated fast relaxation.
        for row in 0..k {
            let mut predicted = Dd::ZERO;
            for mode in 0..k {
                predicted = predicted + dd(self.relative_modes[row * k + mode]) * mean[mode];
            }
            mean_rates[self.active[row]] = predicted;
        }
        let mut phase = [Dd::ZERO; N];
        let mut drive_work = [Dd::ZERO; N];
        let mut work_defect = dd(impulse.kinetic_energy_change_j) - external_work + heat;
        let mut work_scale = impulse
            .kinetic_energy_change_j
            .abs()
            .max(external_work.to_f64().abs())
            .max(heat.to_f64().abs());
        // Native kinetic increments contain opposing absolute-spin transfers
        // between carrier and rotors. Their rounding scales with each actual
        // local transfer, even when total slip work is much smaller. Include
        // those components, never an unchanged momentum/energy inventory.
        let carrier_change = self.bank.carrier_momentum(&impulse.state_after)?
            - self.bank.carrier_momentum(state)?;
        let mean_omega =
            omega.scale(0.5) + self.bank.angular_velocity(&impulse.state_after)?.scale(0.5);
        for (speed, change) in components(mean_omega)
            .into_iter()
            .zip(components(carrier_change))
        {
            work_scale = work_scale.max(
                scalar(
                    (dd(speed) * dd(change)).to_f64(),
                    "exponential carrier transfer work",
                )?
                .abs(),
            );
        }
        for i in 0..N {
            let mean_absolute_spin = (0.5 * state.axial_momentum_nms[i]
                + 0.5 * impulse.state_after.axial_momentum_nms[i])
                / self.bank.rotors[i].axial_inertia_kg_m2;
            work_scale = work_scale.max(
                scalar(
                    (dd(mean_absolute_spin) * dd(impulse.axial_momentum_change_nms[i])).to_f64(),
                    "exponential axial transfer work",
                )?
                .abs(),
            );
            let reconstructed = dd(initial_rate[i])
                + compensated::divide(
                    mean_increment[i],
                    dd(self.bank.rotors[i].axial_inertia_kg_m2),
                )?
                - dot_dd(self.bank.rotors[i].axis_body, mean_velocity_change);
            if self.damping_nms[i] == 0.0 {
                mean_rates[i] = reconstructed;
            } else {
                let scale = initial_rate[i]
                    .abs()
                    .max(reconstructed.to_f64().abs())
                    .max(mean_rates[i].to_f64().abs());
                if scalar(
                    (reconstructed - mean_rates[i]).to_f64(),
                    "exponential mean-rate reconstruction",
                )?
                .abs()
                    > 32768.0 * f64::EPSILON * k as f64 * scale
                {
                    return Err(invalid("resolved exponential shaft mean rate"));
                }
            }
            phase[i] = dd(self.half_s) * mean_rates[i];
            drive_work[i] = dd(drive[i]) * phase[i];
            let drive_impulse = dd(self.half_s) * dd(drive[i]);
            let drag_impulse = dd(self.damping_nms[i]) * phase[i];
            let residual = dd(impulse.axial_momentum_change_nms[i]) - drive_impulse + drag_impulse;
            let scale = impulse.axial_momentum_change_nms[i]
                .abs()
                .max(drive_impulse.to_f64().abs())
                .max(drag_impulse.to_f64().abs());
            if scalar(residual.to_f64(), "exponential shaft law residual")?.abs()
                > 32768.0 * f64::EPSILON * k as f64 * scale
            {
                return Err(invalid("resolved exponential shaft effort"));
            }
            work_defect = work_defect - drive_work[i];
            work_scale =
                work_scale.max(scalar(drive_work[i].to_f64(), "exponential shaft work")?.abs());
        }
        if scalar(heat.to_f64(), "exponential physical shaft heat")? < 0.0
            || scalar(work_defect.to_f64(), "exponential shaft stage work")?.abs()
                > 65536.0 * f64::EPSILON * (N + k + 1) as f64 * work_scale
        {
            return Err(invalid("resolved exponential physical shaft work"));
        }
        Ok(Stage {
            impulse,
            phase,
            drive_work,
            external_work,
            heat,
        })
    }
}

struct Stage<const N: usize> {
    impulse: RotorImpulse<N>,
    phase: [Dd; N],
    drive_work: [Dd; N],
    external_work: Dd,
    heat: Dd,
}

fn as_vector(value: [Dd; 3], context: &'static str) -> Result<Vec3, LieStepError> {
    vector(
        Vec3::new(value[0].to_f64(), value[1].to_f64(), value[2].to_f64()),
        context,
    )
}
