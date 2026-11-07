//! Simultaneous passive shaft loads inside an ideal multi-shaft transmission.
//!
//! With ideal rows A*s=0 and authored load rows C, the actual passive law is
//! torque = drive - C^T(K*(C*q-rest) + D*C*s). Each C row is the original
//! TorsionalCoupling; it is never applied as an unconstrained pair kick.
//!
//! P is the prepared projection of original shaft impulses into the retained
//! free generalized momenta. For a load row c, b=P*c and g=sum_i a_i*b_i.
//! The constrained port mobility is the positive semidefinite Gram form
//! H_ef=sum_i b_ei*b_fi/J_i + g_e.I_retained^-1.g_f. It admits dependent
//! load rows and immobile ports without discarding their static reactions.
//! Eliminating midpoint phase gives D_eff=D+delta*K/2 and the same native SPD
//! midpoint solve used by the original unconstrained passive network.

use super::super::{MAX_TORSIONAL_COUPLINGS, TorsionalCoupling};
use super::{
    HeldTorque, LieStepError, RotorState, RotorStep, RotorTransmission, So3, TransmissionImpulse,
    Vec3, dd, dot_dd, invalid, product, scalar, vector, work_residual,
};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::dd::Dd;

/// Prepared elastic/viscous shaft loads on an admitted ideal transmission.
///
/// All original free differential modes and the complete carrier recoil remain.
/// Preparation projects the load ports once and factors one E-by-E SPD matrix.
/// A step uses O(N*M+M*M+N+E*E) fixed stack work and no allocation or factorization.
/// The midpoint stages resolve finite intervals; they are not exact stiff
/// relaxation or an unconditional trajectory-accuracy guarantee.
#[derive(Debug, Clone)]
pub struct TransmissionNetworkStepper<const N: usize, const M: usize, const E: usize> {
    transmission: RotorTransmission<N, M>,
    couplings: [TorsionalCoupling; E],
    duration_s: f64,
    effective_sqrt: [f64; E],
    factor: Option<Cholesky>,
    angular_frequency_bound_rad_s: f64,
    half_stage_resolution_number: f64,
}

/// Native constrained motion with separate original loads, ideal reactions and work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransmissionNetworkStep<const N: usize, const M: usize, const E: usize> {
    /// Original carrier/rotor trajectory, including numerical gyrostat drift error.
    pub motion: RotorStep<N>,
    /// Ideal row reactions from both loaded stages and free precession [N m s].
    /// Their units and orientation remain those of the authored ShaftConstraints.
    pub constraint_impulse_nms: [f64; M],
    /// Positive coupling effort impulses; shaft impulse is minus C^T times this.
    /// Values follow the original TorsionalCoupling inventory [N m s].
    pub coupling_impulse_nms: [f64; E],
    /// Work at each original authored motor/brake port, excluding passive loads [J].
    pub drive_work_j: [f64; N],
    /// Nonnegative viscous heat from each actual relative coupling rate [J].
    pub dissipated_energy_j: [f64; E],
    /// Polarized change in each original spring potential [J].
    pub elastic_energy_change_j: [f64; E],
    /// Summed ideal reaction work at the two constitutive stages [J].
    pub reaction_work_j: f64,
    /// Delta kinetic + Delta elastic + heat - drives - external - drift defect [J].
    pub work_energy_residual_j: f64,
}

impl<const N: usize, const M: usize> RotorTransmission<N, M> {
    /// Prepare simultaneous original shaft springs/dampers in the retained modes.
    ///
    /// Reuses TorsionalCoupling's original relative phases, signed ratio, rest
    /// phase and passive coefficients. Ratios may be negative or zero; parallel
    /// or dependent load rows are allowed. They do not add ideal constraints.
    /// Positive duration and at most 64 load rows are required. With delta=dt/2,
    /// delta*trace((D+delta*K/2)*H) must be at most two, with a rounding margin.
    /// This numerical PSD-trace resolution gate is not an interval error bound.
    /// Fixed topology, coefficients, ideal rows, axes and inertia are required.
    pub fn torsional_network_stepper<const E: usize>(
        &self,
        couplings: [TorsionalCoupling; E],
        duration_s: f64,
    ) -> Result<TransmissionNetworkStepper<N, M, E>, LieStepError> {
        scalar(duration_s, "loaded transmission duration")?;
        let delta = 0.5 * duration_s;
        if E > MAX_TORSIONAL_COUPLINGS || delta <= 0.0 {
            return Err(invalid(
                "bounded loads and positive resolved transmission duration",
            ));
        }
        let mut projected = [[0.0; N]; E];
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
                return Err(invalid("distinct original loaded transmission shafts"));
            }
            scalar(coupling.ratio, "loaded shaft ratio")?;
            scalar(coupling.rest_angle_rad, "loaded shaft rest phase")?;
            let k = scalar(coupling.stiffness_nm_per_rad, "loaded shaft stiffness")?;
            let d = scalar(coupling.damping_nms_per_rad, "loaded shaft damping")?;
            if k < 0.0 || d < 0.0 {
                return Err(invalid("nonnegative loaded shaft stiffness and damping"));
            }
            if !coupling.active() {
                continue;
            }
            active = true;
            effective[e] = scalar(
                d + 0.5 * product(delta, k, "loaded midpoint spring scale")?,
                "loaded midpoint effective drag",
            )?;
            effective_sqrt[e] = fs_math::det::sqrt(effective[e]);
            if effective_sqrt[e] == 0.0 {
                return Err(invalid("resolved positive loaded effective drag"));
            }
            let mut port = [0.0; N];
            for (i, coefficient) in coupling.ports() {
                port[i] = coefficient;
            }
            let free = self.project_drive(&port)?;
            axes[e] = self.projected_rotor_momentum(&port, &free)?;
            inverse_axes[e] = self.inertia.velocity(axes[e])?;
            for i in 0..N {
                projected[e][i] = scalar(
                    free[i] / self.bank.rotors[i].axial_inertia_kg_m2.sqrt(),
                    "loaded free-mode port",
                )?;
            }
        }
        let mut matrix = vec![0.0; E * E];
        let mut frequency_squared = Dd::ZERO;
        let mut resolution = Dd::ZERO;
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
                for i in 0..N {
                    mobility = mobility + dd(projected[e][i]) * dd(projected[f][i]);
                }
                let h = scalar(mobility.to_f64(), "constrained load port mobility")?;
                if e == f {
                    if h < 0.0 {
                        return Err(invalid("nonnegative constrained load mobility"));
                    }
                    frequency_squared =
                        frequency_squared + dd(couplings[e].stiffness_nm_per_rad) * dd(h);
                    resolution = resolution + dd(delta) * dd(effective[e]) * dd(h);
                }
                let value = scalar(
                    (dd(0.5) * dd(delta) * dd(effective_sqrt[e]) * dd(h) * dd(effective_sqrt[f])
                        + dd(if e == f { 1.0 } else { 0.0 }))
                    .to_f64(),
                    "constrained load midpoint matrix",
                )?;
                matrix[e * E + f] = value;
                matrix[f * E + e] = value;
            }
        }
        let margin = 1.0 + 256.0 * f64::EPSILON;
        let resolution = scalar(resolution.to_f64() * margin, "loaded half-stage resolution")?;
        if resolution > 2.0 {
            return Err(invalid("unresolved loaded transmission half-stage"));
        }
        let frequency = scalar(
            frequency_squared.to_f64() * margin,
            "loaded frequency estimate",
        )?;
        let factor = if active {
            Some(cholesky(&matrix, E).map_err(|_| invalid("loaded transmission Cholesky factor"))?)
        } else {
            None
        };
        Ok(TransmissionNetworkStepper {
            transmission: *self,
            couplings,
            duration_s,
            effective_sqrt,
            factor,
            angular_frequency_bound_rad_s: fs_math::det::sqrt(frequency),
            half_stage_resolution_number: resolution,
        })
    }
}

impl<const N: usize, const M: usize, const E: usize> TransmissionNetworkStepper<N, M, E> {
    /// Prepared positive full interval [s].
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Original passive load laws and row ordering.
    pub const fn couplings(&self) -> &[TorsionalCoupling; E] {
        &self.couplings
    }
    /// Numerical PSD-trace upper estimate, not a certified eigenvalue bound.
    pub const fn angular_frequency_bound_rad_s(&self) -> f64 {
        self.angular_frequency_bound_rad_s
    }
    /// Combined elastic/viscous half-stage resolution number, admitted at most two.
    pub const fn half_stage_resolution_number(&self) -> f64 {
        self.half_stage_resolution_number
    }

    /// Original total spring potential at actual unwrapped shaft phases [J].
    pub fn potential_energy(&self, state: &RotorState<N>) -> Result<f64, LieStepError> {
        self.transmission.bank.admit_state(state)?;
        let mut value = Dd::ZERO;
        for coupling in self.couplings {
            if coupling.stiffness_nm_per_rad != 0.0 {
                let strain = coupling.strain(&state.relative_angle_rad);
                value = value + dd(0.5) * dd(coupling.stiffness_nm_per_rad) * strain * strain;
            }
        }
        scalar(value.to_f64(), "loaded transmission spring potential")
    }

    /// Advance all loads, original drives and external torque together.
    ///
    /// Two implicit midpoint shaft/phase half-stages surround the existing
    /// constrained gyrostat drift. Original ideal reactions remain active in
    /// every stage. Actual phases, constitutive impulses and work must resolve
    /// at local increment scales; any refusal returns no partial state.
    ///
    /// This is a second-order split approximation: refine for torsional
    /// oscillation, damping and gyroscopic motion. Midpoint is not exact stiff
    /// relaxation. Heat is mechanical viscous loss, without an inferred thermal
    /// destination. Changing coefficients/ratios requires explicit repreparation
    /// and separately accounted parameter work.
    pub fn step(
        &self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<TransmissionNetworkStep<N, M, E>, LieStepError> {
        if self.factor.is_none() {
            let result =
                self.transmission
                    .step(state, drive_torque_nm, external_torque, self.duration_s)?;
            return Ok(TransmissionNetworkStep {
                motion: result.motion,
                constraint_impulse_nms: result.constraint_impulse_nms,
                coupling_impulse_nms: [0.0; E],
                drive_work_j: result.drive_work_j,
                dissipated_energy_j: [0.0; E],
                elastic_energy_change_j: [0.0; E],
                reaction_work_j: result.reaction_work_j,
                work_energy_residual_j: result.work_energy_residual_j,
            });
        }
        let first = self.stage(
            state,
            drive_torque_nm,
            self.external_impulse(state.rotation, external_torque)?,
        )?;
        let mut carried = first.transfer.impulse.state_after;
        carried.relative_angle_rad = first.phase;
        let drift = self.transmission.drift_momenta(&carried, self.duration_s)?;
        let second = self.stage(
            &drift.state_after,
            drive_torque_nm,
            self.external_impulse(drift.state_after.rotation, external_torque)?,
        )?;
        let mut motion = self.transmission.bank.finish_step(
            state,
            &first.transfer.impulse,
            &drift.state_after,
            &second.transfer.impulse,
            self.duration_s,
            Some(second.phase),
        )?;
        motion.drift_energy_defect_j = drift.drift_energy_defect_j;
        motion.work_energy_residual_j = work_residual(
            motion.kinetic_energy_change_j,
            &motion.shaft_work_j,
            motion.external_work_j,
            motion.drift_energy_defect_j,
        )?;
        let mut constraints = [0.0; M];
        for e in 0..M {
            constraints[e] = scalar(
                (dd(first.transfer.constraint_impulse_nms[e])
                    + dd(drift.constraint_impulse_nms[e])
                    + dd(second.transfer.constraint_impulse_nms[e]))
                .to_f64(),
                "loaded integrated ideal reaction",
            )?;
        }
        let mut drive_work = [0.0; N];
        let mut coupling_impulse = [0.0; E];
        let mut heat = [0.0; E];
        let elastic = self.potential_changes(&state.relative_angle_rad, &second.phase)?;
        let mut balance = dd(motion.kinetic_energy_change_j)
            - dd(motion.external_work_j)
            - dd(motion.drift_energy_defect_j);
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(first.drive_work[i]) + dd(second.drive_work[i])).to_f64(),
                "loaded authored drive work",
            )?;
            balance = balance - dd(drive_work[i]);
        }
        for e in 0..E {
            coupling_impulse[e] = scalar(
                (dd(first.coupling_impulse[e]) + dd(second.coupling_impulse[e])).to_f64(),
                "loaded integrated coupling effort",
            )?;
            heat[e] = scalar(
                (dd(first.heat[e]) + dd(second.heat[e])).to_f64(),
                "loaded viscous heat",
            )?;
            balance = balance + dd(elastic[e]) + dd(heat[e]);
        }
        Ok(TransmissionNetworkStep {
            motion,
            constraint_impulse_nms: constraints,
            coupling_impulse_nms: coupling_impulse,
            drive_work_j: drive_work,
            dissipated_energy_j: heat,
            elastic_energy_change_j: elastic,
            reaction_work_j: scalar(
                (dd(first.transfer.reaction_work_j) + dd(second.transfer.reaction_work_j)).to_f64(),
                "loaded ideal reaction work",
            )?,
            work_energy_residual_j: scalar(balance.to_f64(), "loaded complete work balance")?,
        })
    }

    fn external_impulse(&self, rotation: So3, torque: HeldTorque) -> Result<Vec3, LieStepError> {
        let body = match torque {
            HeldTorque::Body(value) => vector(value, "loaded external body torque")?,
            HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
        };
        let delta = 0.5 * self.duration_s;
        Ok(Vec3::new(
            product(body.x, delta, "loaded external impulse")?,
            product(body.y, delta, "loaded external impulse")?,
            product(body.z, delta, "loaded external impulse")?,
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
                "loaded elastic increment",
            )?;
        }
        Ok(changes)
    }

    fn stage(
        &self,
        state: &RotorState<N>,
        drive: [f64; N],
        external: Vec3,
    ) -> Result<LoadedStage<N, M, E>, LieStepError> {
        self.transmission.admit_transmission(state)?;
        let bank = &self.transmission.bank;
        let delta = 0.5 * self.duration_s;
        let rates = bank.relative_rates(state)?;
        let delivered =
            super::super::super::super::torque_kick(state.total_momentum_body, external, 1.0)?
                - state.total_momentum_body;
        let mut drive_impulse = [0.0; N];
        for i in 0..N {
            drive_impulse[i] = product(delta, drive[i], "loaded shaft drive impulse")?;
        }
        let (free_transfer, delta_omega) = self
            .transmission
            .impulse_response(&drive_impulse, delivered)?;
        let mut free_mid = [0.0; N];
        for i in 0..N {
            if M == N || self.transmission.fixed_phase[i] {
                continue;
            }
            let rotor = bank.rotors[i];
            free_mid[i] = scalar(
                (dd(rates[i])
                    + dd(0.5)
                        * (dd(free_transfer[i]) / dd(rotor.axial_inertia_kg_m2)
                            - dot_dd(rotor.axis_body, delta_omega)))
                .to_f64(),
                "loaded free midpoint rate",
            )?;
        }
        let mut solved = [0.0; E];
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            if !coupling.active() {
                continue;
            }
            let spring = if coupling.stiffness_nm_per_rad == 0.0 {
                Dd::ZERO
            } else {
                dd(coupling.stiffness_nm_per_rad) * coupling.strain(&state.relative_angle_rad)
            };
            solved[e] = scalar(
                (spring / dd(self.effective_sqrt[e])
                    + dd(self.effective_sqrt[e]) * coupling.mapped(&free_mid))
                .to_f64(),
                "loaded midpoint right hand side",
            )?;
        }
        self.factor
            .as_ref()
            .expect("active loaded factor")
            .solve(&mut solved);
        let mut coupling_impulse = [0.0; E];
        let mut requested = drive.map(|value| dd(delta) * dd(value));
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            scalar(solved[e], "loaded solved effort")?;
            let value = dd(delta) * dd(self.effective_sqrt[e]) * dd(solved[e]);
            coupling_impulse[e] = scalar(value.to_f64(), "loaded original coupling impulse")?;
            for (i, coefficient) in coupling.ports() {
                requested[i] = requested[i] - dd(coefficient) * value;
            }
        }
        let transfer = self.transmission.apply_impulse(
            state,
            requested.map(|value| value.to_f64()),
            external,
        )?;
        let mut mean = transfer.impulse.mean_relative_rate_rad_s;
        for i in 0..N {
            if M == N || self.transmission.fixed_phase[i] {
                mean[i] = 0.0;
            }
        }
        let (phase, phase_increment) = self
            .transmission
            .advance_phase_with_increment(state, &mean, delta)?;
        let allowance = 4096.0 * f64::EPSILON * (2 * N + 2 * M + 2 * E + 2) as f64;
        for i in 0..N {
            let expected = product(delta, mean[i], "loaded phase quadrature")?;
            if (phase_increment[i] - expected).abs()
                > allowance * phase_increment[i].abs().max(expected.abs())
            {
                return Err(invalid("unresolved loaded phase quadrature"));
            }
        }
        let elastic = self.potential_changes(&state.relative_angle_rad, &phase)?;
        let mut heat = [0.0; E];
        let mut balance =
            dd(transfer.impulse.kinetic_energy_change_j) - dd(transfer.impulse.external_work_j);
        let mut work_scale =
            transfer.impulse.kinetic_energy_change_j.abs() + transfer.impulse.external_work_j.abs();
        for (e, coupling) in self.couplings.iter().copied().enumerate() {
            if !coupling.active() {
                continue;
            }
            let rate = coupling.mapped(&mean);
            let spring = if coupling.stiffness_nm_per_rad == 0.0 {
                Dd::ZERO
            } else {
                dd(delta)
                    * dd(coupling.stiffness_nm_per_rad)
                    * dd(0.5)
                    * (coupling.strain(&state.relative_angle_rad) + coupling.strain(&phase))
            };
            let drag = dd(delta) * dd(coupling.damping_nms_per_rad) * rate;
            let residual = scalar(
                (dd(coupling_impulse[e]) - spring - drag).to_f64(),
                "loaded original constitutive impulse residual",
            )?;
            let scale = scalar(
                (dd(coupling_impulse[e]).abs() + spring.abs() + drag.abs()).to_f64(),
                "loaded constitutive impulse scale",
            )?;
            if residual.abs() > allowance * scale {
                return Err(invalid("unresolved original loaded spring or damping law"));
            }
            heat[e] = scalar((drag * rate).to_f64(), "loaded physical viscous heat")?;
            if heat[e] < 0.0 {
                return Err(invalid("nonnegative loaded viscous heat"));
            }
            balance = balance + dd(elastic[e]) + dd(heat[e]);
            work_scale += elastic[e].abs() + heat[e];
        }
        let mut drive_work = [0.0; N];
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(delta) * dd(drive[i]) * dd(mean[i])).to_f64(),
                "loaded shaft drive work",
            )?;
            balance = balance - dd(drive_work[i]);
            work_scale += drive_work[i].abs();
        }
        if scalar(balance.to_f64(), "loaded stage energy balance")?.abs()
            > allowance * scalar(work_scale, "loaded stage work scale")?
        {
            return Err(invalid("unresolved loaded transmission potential work"));
        }
        Ok(LoadedStage {
            transfer,
            phase,
            coupling_impulse,
            drive_work,
            heat,
        })
    }
}

struct LoadedStage<const N: usize, const M: usize, const E: usize> {
    transfer: TransmissionImpulse<N, M>,
    phase: [f64; N],
    coupling_impulse: [f64; E],
    drive_work: [f64; N],
    heat: [f64; E],
}
