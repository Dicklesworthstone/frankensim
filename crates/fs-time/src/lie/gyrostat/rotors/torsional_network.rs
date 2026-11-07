//! Reciprocal elastic and viscous transmission between actual rotor shafts.
//!
//! For each row of B, x = q_first - ratio*q_second - rest, with
//! V = sum(k*x*x/2) and shaft torque = drive - B^T(k*x + d*B*q').
//! The existing rotor bank owns inertia, carrier recoil, impulses and drift.
//! At fixed attitude its relative-rate mobility is M = diag(1/J)+A^T I^-1 A.
//! Eliminating midpoint phase gives D_eff = D + delta*K/2 and the SPD solve
//! (1 + delta/2 sqrt(D_eff) B M B^T sqrt(D_eff)) z =
//! K*x_before/sqrt(D_eff) + sqrt(D_eff)*B*s_free_midpoint.
//! Coupling effort is sqrt(D_eff)*z; zero-coefficient rows contribute nothing.

use super::{
    HeldTorque, LieStepError, RotorBank, RotorImpulse, RotorState, RotorStep, components, dd,
    dot_dd, increment, invalid, product, scalar, vector,
};
use fs_ga::{So3, Vec3};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::dd::Dd;

/// Fixed preparation and stack-work bound for one transmission network.
pub const MAX_TORSIONAL_COUPLINGS: usize = 64;

/// One passive, finite-compliance transmission between two distinct shafts.
///
/// Strain is `q[first] - ratio*q[second] - rest_angle_rad`, using the bank's
/// unwrapped angles RELATIVE TO THE CARRIER. Positive effort acts as `-effort`
/// on the first shaft and `+ratio*effort` on the second. Thus negative ratios
/// describe opposite-sense transmission and carrier reaction is automatic.
/// A zero ratio is an explicit first-shaft/carrier spring/drag limit.
/// Coefficients and ratio remain fixed: this is not rigid gearing, tooth
/// contact, backlash, changing-ratio actuation or a thermal material model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionalCoupling {
    /// Original shaft receiving minus the coupling effort.
    pub first_rotor: usize,
    /// Distinct original shaft receiving ratio times the coupling effort.
    pub second_rotor: usize,
    /// Dimensionless signed coefficient of the second shaft's relative angle.
    pub ratio: f64,
    /// Nonnegative stiffness conjugate to the declared strain [N m/rad].
    pub stiffness_nm_per_rad: f64,
    /// Nonnegative drag conjugate to its rate [N m s/rad].
    pub damping_nms_per_rad: f64,
    /// Finite stress-free unwrapped relative phase [rad].
    pub rest_angle_rad: f64,
}

impl TorsionalCoupling {
    pub(super) fn ports(self) -> [(usize, f64); 2] {
        [(self.first_rotor, 1.0), (self.second_rotor, -self.ratio)]
    }
    pub(super) fn active(self) -> bool {
        self.stiffness_nm_per_rad != 0.0 || self.damping_nms_per_rad != 0.0
    }
    pub(super) fn mapped<const N: usize>(self, values: &[f64; N]) -> Dd {
        dd(values[self.first_rotor]) - dd(self.ratio) * dd(values[self.second_rotor])
    }
    pub(super) fn strain<const N: usize>(self, phase: &[f64; N]) -> Dd {
        self.mapped(phase) - dd(self.rest_angle_rad)
    }
}

/// Prepared simultaneous compliant transmission on a complete native rotor bank.
///
/// Preparation uses O(E^3) work and O(E^2) storage for E <= 64 couplings.
/// Each step uses O(N+E^2) bounded work, fixed stack arrays and no allocation.
/// Two midpoint shaft half-stages surround the original gyrostat drift. This
/// is a second-order split approximation; refine time for oscillations, decay
/// and gyroscopic motion. Parallel edges and closed coupling loops are allowed:
/// the prepared matrix remains SPD even when B has dependent rows.
#[derive(Debug, Clone)]
pub struct TorsionalNetworkStepper<const N: usize, const E: usize> {
    bank: RotorBank<N>,
    couplings: [TorsionalCoupling; E],
    duration_s: f64,
    effective_sqrt: [f64; E],
    factor: Option<Cholesky>,
    angular_frequency_bound_rad_s: f64,
    half_stage_resolution_number: f64,
}

/// Original mechanical trajectory, supplied work, coupling storage and heat.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionalNetworkStep<const N: usize, const E: usize> {
    /// Original bank receipt, including actual impulses and native drift defect.
    pub motion: RotorStep<N>,
    /// Signed authored work at each original relative-shaft drive port [J].
    pub drive_work_j: [f64; N],
    /// Nonnegative viscous heat in each authored coupling, in input order [J].
    pub dissipated_energy_j: [f64; E],
    /// Polarized change in each original spring potential [J].
    pub elastic_energy_change_j: [f64; E],
    /// Delta kinetic + Delta elastic + heat - drive - external - drift defect.
    /// Numerical error stays distinct from material dissipation.
    pub work_energy_residual_j: f64,
}

impl<const N: usize> RotorBank<N> {
    /// Prepare simultaneous shaft-to-shaft springs and dampers on fixed axes.
    ///
    /// Every pair names distinct original rotor indices; ratios may have either
    /// sign or be zero. Stiffness and damping must be nonnegative. One native
    /// Cholesky factor is cached; no diagonal carrier approximation, matrix
    /// inverse, sequential pair kick or new time integrator is substituted.
    /// With delta=duration/2, preparation requires
    /// delta*trace((D+delta*K/2)*B*M*B^T) <= 2, with a rounding margin.
    /// This conservative resolution check is not a trajectory-error bound.
    /// Changing topology, ratio, rest phase, coefficients or duration requires
    /// explicit repreparation; associated parameter-change work is not inferred.
    pub fn torsional_network_stepper<const E: usize>(
        &self,
        couplings: [TorsionalCoupling; E],
        duration_s: f64,
    ) -> Result<TorsionalNetworkStepper<N, E>, LieStepError> {
        scalar(duration_s, "torsional network duration")?;
        let delta = 0.5 * duration_s;
        if E > MAX_TORSIONAL_COUPLINGS || delta <= 0.0 {
            return Err(invalid(
                "bounded coupling count and positive resolved duration",
            ));
        }
        let mut axes = [Vec3::default(); E];
        let mut inverse_axes = [Vec3::default(); E];
        let mut effective = [0.0; E];
        let mut effective_sqrt = [0.0; E];
        let mut active = false;
        for (e, coupling) in couplings.iter().copied().enumerate() {
            if coupling.first_rotor >= N
                || coupling.second_rotor >= N
                || coupling.first_rotor == coupling.second_rotor
            {
                return Err(invalid("distinct original torsional coupling shafts"));
            }
            scalar(coupling.ratio, "torsional transmission ratio")?;
            scalar(coupling.rest_angle_rad, "torsional coupling rest phase")?;
            let k = scalar(coupling.stiffness_nm_per_rad, "coupling stiffness")?;
            let d = scalar(coupling.damping_nms_per_rad, "coupling damping")?;
            if k < 0.0 || d < 0.0 {
                return Err(invalid("nonnegative coupling stiffness and damping"));
            }
            if !coupling.active() {
                continue;
            }
            active = true;
            let spring = product(delta, k, "network midpoint spring scale")?;
            effective[e] = scalar(d + 0.5 * spring, "network effective drag")?;
            effective_sqrt[e] = fs_math::det::sqrt(effective[e]);
            if effective_sqrt[e] == 0.0 {
                return Err(invalid("resolved positive network effective drag"));
            }
            let a = components(self.rotors[coupling.first_rotor].axis_body);
            let b = components(self.rotors[coupling.second_rotor].axis_body);
            let axis = core::array::from_fn::<_, 3, _>(|i| {
                (dd(a[i]) - dd(coupling.ratio) * dd(b[i])).to_f64()
            });
            axes[e] = vector(Vec3::new(axis[0], axis[1], axis[2]), "coupling recoil axis")?;
            inverse_axes[e] = self.inertia.velocity(axes[e])?;
        }
        let mut matrix = vec![0.0; E * E];
        let mut frequency_squared = dd(0.0);
        let mut resolution = dd(0.0);
        for e in 0..E {
            matrix[e * E + e] = 1.0;
            if !couplings[e].active() {
                continue;
            }
            for f in 0..=e {
                if !couplings[f].active() {
                    continue;
                }
                let mut mobility = dot_dd(axes[e], inverse_axes[f]);
                for (i, a) in couplings[e].ports() {
                    for (j, b) in couplings[f].ports() {
                        if i == j {
                            mobility =
                                mobility + dd(a) * dd(b) / dd(self.rotors[i].axial_inertia_kg_m2);
                        }
                    }
                }
                let h = scalar(mobility.to_f64(), "network port mobility")?;
                if e == f {
                    if h <= 0.0 {
                        return Err(invalid("positive network self-mobility"));
                    }
                    frequency_squared =
                        frequency_squared + dd(couplings[e].stiffness_nm_per_rad) * dd(h);
                    resolution = resolution + dd(delta) * dd(effective[e]) * dd(h);
                }
                let value = scalar(
                    ((dd(delta) * dd(0.5)) * dd(effective_sqrt[e]) * dd(h) * dd(effective_sqrt[f])
                        + dd(if e == f { 1.0 } else { 0.0 }))
                    .to_f64(),
                    "network midpoint matrix",
                )?;
                matrix[e * E + f] = value;
                matrix[f * E + e] = value;
            }
        }
        let margin = 1.0 + 256.0 * f64::EPSILON;
        let resolution = scalar(
            resolution.to_f64() * margin,
            "network half-stage resolution",
        )?;
        if resolution > 2.0 {
            return Err(invalid("unresolved torsional network half-stage"));
        }
        let frequency = scalar(
            frequency_squared.to_f64() * margin,
            "network frequency estimate",
        )?;
        let factor = if active {
            Some(cholesky(&matrix, E).map_err(|_| invalid("network midpoint Cholesky factor"))?)
        } else {
            None
        };
        Ok(TorsionalNetworkStepper {
            bank: *self,
            couplings,
            duration_s,
            effective_sqrt,
            factor,
            angular_frequency_bound_rad_s: fs_math::det::sqrt(frequency),
            half_stage_resolution_number: resolution,
        })
    }
}

impl<const N: usize, const E: usize> TorsionalNetworkStepper<N, E> {
    /// Prepared positive full interval [s].
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Original authored topology and coefficients, in receipt order.
    pub const fn couplings(&self) -> &[TorsionalCoupling; E] {
        &self.couplings
    }
    /// Numerical PSD-trace upper estimate, not an interval eigenvalue proof.
    pub const fn angular_frequency_bound_rad_s(&self) -> f64 {
        self.angular_frequency_bound_rad_s
    }
    /// Combined elastic/viscous half-stage resolution number, admitted at most 2.
    pub const fn half_stage_resolution_number(&self) -> f64 {
        self.half_stage_resolution_number
    }

    /// Original total spring potential at the actual unwrapped phases [J].
    pub fn potential_energy(&self, state: &RotorState<N>) -> Result<f64, LieStepError> {
        self.bank.admit_state(state)?;
        let mut energy = dd(0.0);
        for coupling in self.couplings {
            if coupling.stiffness_nm_per_rad != 0.0 {
                let strain = coupling.strain(&state.relative_angle_rad);
                energy = energy + dd(0.5) * dd(coupling.stiffness_nm_per_rad) * strain * strain;
            }
        }
        scalar(energy.to_f64(), "network spring potential")
    }

    /// Advance all couplings and authored drives simultaneously.
    ///
    /// Actual endpoint phases, original port torques and kinetic/potential work
    /// must resolve at local increment scales. The existing carrier-recoil and
    /// impulse admission stays authoritative. A refusal publishes no half-step.
    /// Empty or entirely zero-law networks delegate exactly to `RotorBank::step`.
    pub fn step(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<TorsionalNetworkStep<N, E>, LieStepError> {
        if self.factor.is_none() {
            let motion =
                self.bank
                    .step(state, drive_torque_nm, external_torque, self.duration_s)?;
            return Ok(TorsionalNetworkStep {
                drive_work_j: motion.shaft_work_j,
                motion,
                dissipated_energy_j: [0.0; E],
                elastic_energy_change_j: [0.0; E],
                work_energy_residual_j: motion.work_energy_residual_j,
            });
        }
        let first = self.stage(
            state,
            drive_torque_nm,
            self.external_impulse(state.rotation, external_torque)?,
        )?;
        let mut drifted = first.impulse.state_after;
        drifted.relative_angle_rad = first.phase;
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
            Some(second.phase),
        )?;
        let elastic_energy_change_j =
            self.potential_changes(&state.relative_angle_rad, &second.phase)?;
        let mut drive_work_j = [0.0; N];
        let mut dissipated_energy_j = [0.0; E];
        let mut balance = dd(motion.kinetic_energy_change_j)
            - dd(motion.external_work_j)
            - dd(motion.drift_energy_defect_j);
        for i in 0..N {
            drive_work_j[i] = scalar(
                (dd(first.drive_work[i]) + dd(second.drive_work[i])).to_f64(),
                "network drive work",
            )?;
            balance = balance - dd(drive_work_j[i]);
        }
        for e in 0..E {
            dissipated_energy_j[e] = scalar(
                (dd(first.heat[e]) + dd(second.heat[e])).to_f64(),
                "network viscous heat",
            )?;
            balance = balance + dd(elastic_energy_change_j[e]) + dd(dissipated_energy_j[e]);
        }
        Ok(TorsionalNetworkStep {
            motion,
            drive_work_j,
            dissipated_energy_j,
            elastic_energy_change_j,
            work_energy_residual_j: scalar(balance.to_f64(), "network complete work balance")?,
        })
    }

    fn external_impulse(&self, rotation: So3, torque: HeldTorque) -> Result<Vec3, LieStepError> {
        let body = match torque {
            HeldTorque::Body(value) => vector(value, "network external body torque")?,
            HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
        };
        let delta = 0.5 * self.duration_s;
        Ok(Vec3::new(
            product(body.x, delta, "network external impulse")?,
            product(body.y, delta, "network external impulse")?,
            product(body.z, delta, "network external impulse")?,
        ))
    }

    fn potential_changes(
        &self,
        before: &[f64; N],
        after: &[f64; N],
    ) -> Result<[f64; E], LieStepError> {
        let mut changes = [0.0; E];
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            if coupling.stiffness_nm_per_rad == 0.0 {
                continue;
            }
            let mean = (coupling.strain(before) + coupling.strain(after)) * dd(0.5);
            let delta = (dd(after[coupling.first_rotor]) - dd(before[coupling.first_rotor]))
                - dd(coupling.ratio)
                    * (dd(after[coupling.second_rotor]) - dd(before[coupling.second_rotor]));
            changes[e] = scalar(
                (dd(coupling.stiffness_nm_per_rad) * mean * delta).to_f64(),
                "network elastic increment",
            )?;
        }
        Ok(changes)
    }

    fn stage(
        &self,
        state: &RotorState<N>,
        drive: [f64; N],
        external: Vec3,
    ) -> Result<NetworkStage<N, E>, LieStepError> {
        self.bank.admit_state(state)?;
        let delta = 0.5 * self.duration_s;
        let initial_rates = self.bank.relative_rates(state)?;
        let next_total =
            super::super::super::torque_kick(state.total_momentum_body, external, 1.0)?;
        let mut free_carrier = components(next_total - state.total_momentum_body).map(dd);
        let mut impulse_values = [dd(0.0); N];
        for i in 0..N {
            product(delta, drive[i], "network drive impulse")?;
            impulse_values[i] = dd(delta) * dd(drive[i]);
            for (value, axis) in free_carrier
                .iter_mut()
                .zip(components(self.bank.rotors[i].axis_body))
            {
                *value = *value - dd(axis) * impulse_values[i];
            }
        }
        let dv = self.bank.inertia.velocity(Vec3::new(
            free_carrier[0].to_f64(),
            free_carrier[1].to_f64(),
            free_carrier[2].to_f64(),
        ))?;
        let mut free_mid = [0.0; N];
        for i in 0..N {
            let rotor = self.bank.rotors[i];
            free_mid[i] = scalar(
                (dd(initial_rates[i])
                    + dd(0.5)
                        * (impulse_values[i] / dd(rotor.axial_inertia_kg_m2)
                            - dot_dd(rotor.axis_body, dv)))
                .to_f64(),
                "network free midpoint rate",
            )?;
        }
        let mut scaled_effort = [0.0; E];
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            if !coupling.active() {
                continue;
            }
            let spring = if coupling.stiffness_nm_per_rad == 0.0 {
                dd(0.0)
            } else {
                dd(coupling.stiffness_nm_per_rad) * coupling.strain(&state.relative_angle_rad)
            };
            scaled_effort[e] = scalar(
                (spring / dd(self.effective_sqrt[e])
                    + dd(self.effective_sqrt[e]) * coupling.mapped(&free_mid))
                .to_f64(),
                "network midpoint right hand side",
            )?;
        }
        self.factor
            .as_ref()
            .expect("active network factor")
            .solve(&mut scaled_effort);
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            scalar(scaled_effort[e], "network solved effort")?;
            let reaction = dd(delta) * dd(self.effective_sqrt[e]) * dd(scaled_effort[e]);
            for (i, coefficient) in coupling.ports() {
                impulse_values[i] = impulse_values[i] - dd(coefficient) * reaction;
            }
        }
        let impulse =
            self.bank
                .apply_impulse(state, impulse_values.map(|v| v.to_f64()), external)?;
        let allowance = 4096.0 * f64::EPSILON * (2 * N + 2 * E + 2) as f64;
        let mut phase = state.relative_angle_rad;
        for i in 0..N {
            let requested = product(
                delta,
                impulse.mean_relative_rate_rad_s[i],
                "network phase increment",
            )?;
            phase[i] = increment(phase[i], requested, "network relative phase")?;
            let actual = (dd(phase[i]) - dd(state.relative_angle_rad[i])).to_f64();
            if (actual - requested).abs() > allowance * actual.abs().max(requested.abs()) {
                return Err(invalid("unresolved network phase increment"));
            }
        }
        let elastic = self.potential_changes(&state.relative_angle_rad, &phase)?;
        let mut expected = drive.map(|v| dd(delta) * dd(v));
        let mut law_scale = expected.map(|v| v.to_f64().abs());
        let mut heat = [0.0; E];
        let mut balance = dd(impulse.kinetic_energy_change_j) - dd(impulse.external_work_j);
        let mut work_scale = impulse
            .kinetic_energy_change_j
            .abs()
            .max(impulse.external_work_j.abs());
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            if !coupling.active() {
                continue;
            }
            let rate = coupling.mapped(&impulse.mean_relative_rate_rad_s);
            let spring = if coupling.stiffness_nm_per_rad == 0.0 {
                dd(0.0)
            } else {
                dd(delta)
                    * dd(coupling.stiffness_nm_per_rad)
                    * dd(0.5)
                    * (coupling.strain(&state.relative_angle_rad) + coupling.strain(&phase))
            };
            let drag = dd(delta) * dd(coupling.damping_nms_per_rad) * rate;
            for (i, coefficient) in coupling.ports() {
                let reaction = dd(coefficient) * (spring + drag);
                expected[i] = expected[i] - reaction;
                law_scale[i] =
                    law_scale[i].max(scalar(reaction.to_f64(), "network original reaction")?.abs());
            }
            heat[e] = scalar((drag * rate).to_f64(), "network physical heat")?;
            if heat[e] < 0.0 {
                return Err(invalid("nonnegative network physical heat"));
            }
            balance = balance + dd(elastic[e]) + dd(heat[e]);
            work_scale = work_scale.max(elastic[e].abs()).max(heat[e]);
        }
        let mut drive_work = [0.0; N];
        for i in 0..N {
            let actual = impulse.axial_momentum_change_nms[i];
            let required = scalar(expected[i].to_f64(), "network original shaft impulse")?;
            let residual = scalar(
                (dd(actual) - expected[i]).to_f64(),
                "network shaft law residual",
            )?;
            if (required != 0.0 && actual == 0.0)
                || residual.abs() > allowance * law_scale[i].max(actual.abs())
            {
                return Err(invalid("unresolved reciprocal network torque"));
            }
            drive_work[i] = scalar(
                (dd(delta) * dd(drive[i]) * dd(impulse.mean_relative_rate_rad_s[i])).to_f64(),
                "network shaft drive work",
            )?;
            balance = balance - dd(drive_work[i]);
            work_scale = work_scale.max(drive_work[i].abs());
        }
        if scalar(balance.to_f64(), "network stage energy balance")?.abs() > allowance * work_scale
        {
            return Err(invalid("unresolved network potential work"));
        }
        Ok(NetworkStage {
            impulse,
            phase,
            drive_work,
            heat,
        })
    }
}

struct NetworkStage<const N: usize, const E: usize> {
    impulse: RotorImpulse<N>,
    phase: [f64; N],
    drive_work: [f64; N],
    heat: [f64; E],
}
