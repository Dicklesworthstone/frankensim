//! Ideal fixed-ratio transmission on the original finite rotor inventory.
//!
//! In each independent gear component, relative rates are s_i = c_i * nu.
//! With J_g = sum(J_i*c_i*c_i), b_g = sum(J_i*c_i*a_i)/J_g and
//! v_i = a_i - c_i*b_g, eliminating the constrained rates gives
//! I_retained = I + sum(J_i*v_i*v_i^T).
//! A change in carrier rate changes p_i by J_i*v_i.dot(delta_omega).
//! The remaining generalized momentum sum(c_i*p_i) is an actual free shaft
//! coordinate, not a prescribed speed. The original gyrostat advances this
//! retained tensor, and the original impulse owner admits every recoil.

use super::{
    HeldTorque, InertiaTensor3, LieStepError, MAX_AXIAL_ROTORS, Mat3, RotorBank, RotorImpulse,
    RotorState, RotorStep, So3, Vec3, components, dd, dot_dd, increment, invalid, product, scalar,
    subtract3, vector, work_residual,
};
use fs_math::dd::Dd;

const ROUNDING: f64 = 8192.0 * f64::EPSILON;
const NONE: usize = usize::MAX;

/// A constant ideal ratio between two distinct original relative shaft angles.
///
/// The retained relation is q_first' = ratio * q_second'. A positive constraint
/// impulse adds lambda to the first axial momentum and -ratio*lambda to the
/// second; carrier reaction follows from the original total momentum.
/// Ratios may have either sign but must be nonzero. A carrier brake belongs to
/// RotorBank::lock_mode. No tooth geometry, backlash or friction is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GearConstraint {
    /// First original shaft in the ordered relation.
    pub first_rotor: usize,
    /// Distinct second original shaft.
    pub second_rotor: usize,
    /// Signed multiplier of the second relative rate.
    pub ratio: f64,
}

/// Prepared independent fixed-ratio gear components on one native rotor bank.
///
/// Each connected component retains one free relative shaft coordinate.
/// Acyclic constraints give unique mesh reactions. Every cycle, including a
/// redundant consistent cycle or a contradictory ratio loop, refuses instead
/// of selecting arbitrary reactions or silently locking the entire train.
/// Preparation and each step use O(N+E) fixed-size work, apart from one native
/// 3x3 tensor preparation. No matrix factorization or allocation occurs in a step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GearTrain<const N: usize, const E: usize> {
    bank: RotorBank<N>,
    constraints: [GearConstraint; E],
    group: [usize; N],
    relative_scale: [f64; N],
    drive_share: [f64; N],
    recoil_axis: [Vec3; N],
    parent: [usize; N],
    parent_edge: [usize; N],
    order: [usize; N],
    inertia: InertiaTensor3,
}

/// Explicit perfectly inelastic engagement; relative phases are not reset.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GearEngagement<const N: usize, const E: usize> {
    /// Original accepted momentum transfer and its kinetic/work diagnostics.
    pub impulse: RotorImpulse<N>,
    /// Signed actual mesh impulses in constraint order, N m s.
    /// Reciprocal mesh impulses in original constraint order, N m s.
    pub constraint_impulse_nms: [f64; E],
    /// Nonnegative loss from the original carrier-plus-rotor kinetic form, J.
    pub dissipated_energy_j: f64,
    /// Remaining signed work/kinetic reconstruction defect, J.
    pub work_energy_residual_j: f64,
}

/// Simultaneous original drives, ideal mesh reactions and external impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GearedRotorImpulse<const N: usize, const E: usize> {
    /// Native total shaft transfers, including drive and mesh reaction.
    /// Original accepted momentum transfer and its kinetic/work diagnostics.
    pub impulse: RotorImpulse<N>,
    /// Reciprocal mesh impulses in original constraint order, N m s.
    pub constraint_impulse_nms: [f64; E],
    /// Signed work at each authored original shaft port, J.
    pub drive_work_j: [f64; N],
    /// Sum of the mesh reactions' signed works; admitted roundoff, never heat.
    /// Summed signed ideal mesh-reaction work, J.
    pub reaction_work_j: f64,
    /// Remaining signed work/kinetic reconstruction defect, J.
    pub work_energy_residual_j: f64,
}

/// Free constrained motion, retaining each component's generalized shaft spin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GearedRotorDrift<const N: usize, const E: usize> {
    /// Complete native state for the next interval.
    pub state_after: RotorState<N>,
    /// Actual absolute axial-momentum changes in inventory order, N m s.
    pub axial_momentum_change_nms: [f64; N],
    /// Reciprocal mesh impulses in original constraint order, N m s.
    pub constraint_impulse_nms: [f64; E],
    /// Change in the original carrier-plus-rotor kinetic energy, J.
    pub kinetic_energy_change_j: f64,
    /// Native reduced-gyrostat splitting error, never gear loss.
    pub drift_energy_defect_j: f64,
    /// Remaining signed work/kinetic reconstruction defect, J.
    pub work_energy_residual_j: f64,
    /// Spatial angular momentum change minus external impulse, N m s.
    pub angular_momentum_residual_world: Vec3,
}

/// Driven constrained attitude, shaft spin and shared relative-phase motion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GearedRotorStep<const N: usize, const E: usize> {
    /// Native trajectory; shaft_work_j includes ideal reaction transfers.
    pub motion: RotorStep<N>,
    /// Mesh impulses include both driven endpoints and constrained precession.
    /// Reciprocal mesh impulses in original constraint order, N m s.
    pub constraint_impulse_nms: [f64; E],
    /// Authored drive work, excluding the ideal constraint reactions.
    pub drive_work_j: [f64; N],
    /// Summed signed ideal mesh-reaction work, J.
    pub reaction_work_j: f64,
    /// Delta kinetic - drives - external work - native drift defect.
    /// Remaining signed work/kinetic reconstruction defect, J.
    pub work_energy_residual_j: f64,
}

impl<const N: usize> RotorBank<N> {
    /// Prepare fixed signed gear ratios without changing any mechanical state.
    ///
    /// Constraints must form a forest of distinct original shafts. A state
    /// must explicitly engage, or already satisfy the retained rate relations.
    /// Existing phase offsets are retained; no tooth-phase capture is invented.
    /// Unrepresentable composed ratios or reduced inertia refuse.
    pub fn gear_train<const E: usize>(
        &self,
        constraints: [GearConstraint; E],
    ) -> Result<GearTrain<N, E>, LieStepError> {
        if N > MAX_AXIAL_ROTORS || E > MAX_AXIAL_ROTORS || (E != 0 && E >= N) {
            return Err(invalid("independent bounded gear constraints"));
        }
        let mut head = [NONE; N];
        let mut next = [NONE; 2 * MAX_AXIAL_ROTORS];
        for (e, relation) in constraints.iter().copied().enumerate() {
            if relation.first_rotor >= N
                || relation.second_rotor >= N
                || relation.first_rotor == relation.second_rotor
            {
                return Err(invalid("distinct original gear shafts"));
            }
            if scalar(relation.ratio, "finite gear ratio")? == 0.0 {
                return Err(invalid("nonzero ideal gear ratio"));
            }
            next[2 * e] = head[relation.first_rotor];
            head[relation.first_rotor] = 2 * e;
            next[2 * e + 1] = head[relation.second_rotor];
            head[relation.second_rotor] = 2 * e + 1;
        }
        let mut group = [NONE; N];
        let mut parent = [NONE; N];
        let mut parent_edge = [NONE; N];
        let mut order = [0; N];
        let mut relative_scale = [1.0; N];
        let mut count = [0; N];
        let mut used = 0;
        for root in 0..N {
            if group[root] != NONE {
                continue;
            }
            let start = used;
            group[root] = root;
            parent[root] = root;
            order[used] = root;
            used += 1;
            let mut cursor = start;
            while cursor < used {
                let node = order[cursor];
                cursor += 1;
                let mut arc = head[node];
                while arc != NONE {
                    let e = arc / 2;
                    let relation = constraints[e];
                    let other = if node == relation.first_rotor {
                        relation.second_rotor
                    } else {
                        relation.first_rotor
                    };
                    if parent_edge[node] != e {
                        if group[other] != NONE {
                            return Err(invalid("cyclic or redundant ideal gear constraints"));
                        }
                        group[other] = root;
                        parent[other] = node;
                        parent_edge[other] = e;
                        relative_scale[other] = if node == relation.first_rotor {
                            scalar(relative_scale[node] / relation.ratio, "composed gear ratio")?
                        } else {
                            product(relative_scale[node], relation.ratio, "composed gear ratio")?
                        };
                        if relative_scale[other] == 0.0 {
                            return Err(invalid("unresolved composed gear ratio"));
                        }
                        order[used] = other;
                        used += 1;
                    }
                    arc = next[arc];
                }
            }
            count[root] = used - start;
            let scale = order[start..used]
                .iter()
                .map(|&i| relative_scale[i].abs())
                .fold(0.0, f64::max);
            for &i in &order[start..used] {
                relative_scale[i] /= scale;
                if relative_scale[i] == 0.0 {
                    return Err(invalid("unresolved normalized gear ratio"));
                }
            }
        }
        let mut polar = [Dd::ZERO; N];
        let mut weighted_axis = [[Dd::ZERO; 3]; N];
        for i in 0..N {
            let root = group[i];
            let weighted = dd(self.rotors[i].axial_inertia_kg_m2) * dd(relative_scale[i]);
            polar[root] = polar[root] + weighted * dd(relative_scale[i]);
            for (value, axis) in weighted_axis[root]
                .iter_mut()
                .zip(components(self.rotors[i].axis_body))
            {
                *value = *value + weighted * dd(axis);
            }
        }
        let mut drive_share = [0.0; N];
        let mut recoil_axis = [Vec3::default(); N];
        let mut projected = [[Dd::ZERO; 3]; N];
        for i in 0..N {
            if count[group[i]] == 1 {
                drive_share[i] = 1.0;
                continue;
            }
            let root = group[i];
            if scalar(polar[root].to_f64(), "gear component polar inertia")? <= 0.0 {
                return Err(invalid("positive gear component polar inertia"));
            }
            drive_share[i] = scalar(
                (dd(self.rotors[i].axial_inertia_kg_m2) * dd(relative_scale[i]) / polar[root])
                    .to_f64(),
                "geared source momentum share",
            )?;
            let axes = components(self.rotors[i].axis_body);
            let mut recoil = [0.0; 3];
            for d in 0..3 {
                projected[i][d] =
                    dd(axes[d]) - dd(relative_scale[i]) * weighted_axis[root][d] / polar[root];
                recoil[d] = scalar(
                    (dd(self.rotors[i].axial_inertia_kg_m2) * projected[i][d]).to_f64(),
                    "gear carrier recoil coefficient",
                )?;
            }
            recoil_axis[i] = Vec3::new(recoil[0], recoil[1], recoil[2]);
        }
        let inertia = if E == 0 {
            self.inertia
        } else {
            let mut matrix = self.inertia.matrix().m;
            for row in 0..3 {
                for col in row..3 {
                    let mut entry = dd(matrix[3 * row + col]);
                    for i in 0..N {
                        entry = entry
                            + dd(self.rotors[i].axial_inertia_kg_m2)
                                * projected[i][row]
                                * projected[i][col];
                    }
                    matrix[3 * row + col] = scalar(entry.to_f64(), "retained gear inertia")?;
                    matrix[3 * col + row] = matrix[3 * row + col];
                }
            }
            InertiaTensor3::new(Mat3 { m: matrix })?
        };
        Ok(GearTrain {
            bank: *self,
            constraints,
            group,
            relative_scale,
            drive_share,
            recoil_axis,
            parent,
            parent_edge,
            order,
            inertia,
        })
    }
}

impl<const N: usize, const E: usize> GearTrain<N, E> {
    /// Original ordered fixed-ratio relations.
    pub const fn constraints(&self) -> &[GearConstraint; E] {
        &self.constraints
    }

    /// Original reduced carrier inertia plus only the constrained axial modes.
    pub const fn retained_inertia(&self) -> InertiaTensor3 {
        self.inertia
    }

    /// Capture all declared ratios in one perfectly inelastic internal event.
    ///
    /// Total momentum, attitude and all original phases are unchanged.
    /// Loss is the positive original kinetic quadratic form of the actual
    /// momentum change. Already-admitted states are retained exactly.
    pub fn engage(&self, state: &RotorState<N>) -> Result<GearEngagement<N, E>, LieStepError> {
        self.bank.admit_state(state)?;
        let impulse = if self.admit_geared(state).is_ok() {
            self.bank.apply_impulse(state, [0.0; N], Vec3::default())?
        } else {
            let rates = self.bank.relative_rates(state)?;
            let grouped = self.group_rates(&rates)?;
            let mut transfer = [0.0; N];
            for i in 0..N {
                transfer[i] = scalar(
                    (dd(self.bank.rotors[i].axial_inertia_kg_m2)
                        * (dd(self.relative_scale[i]) * dd(grouped[self.group[i]]) - dd(rates[i])))
                    .to_f64(),
                    "gear engagement axial transfer",
                )?;
            }
            let recoil = self.bank.rotor_momentum(&transfer)?.scale(-1.0);
            let delta_omega = self.inertia.velocity(recoil)?;
            for (i, value) in transfer.iter_mut().enumerate() {
                *value = scalar(
                    (dd(*value) + dot_dd(self.recoil_axis[i], delta_omega)).to_f64(),
                    "coupled gear engagement transfer",
                )?;
            }
            self.bank.apply_impulse(state, transfer, Vec3::default())?
        };
        self.admit_geared(&impulse.state_after)?;
        let constraints =
            self.constraint_impulses(&impulse.axial_momentum_change_nms, &[0.0; N])?;
        let delta_carrier = self.bank.carrier_momentum(&impulse.state_after)?
            - self.bank.carrier_momentum(state)?;
        let mut loss = dot_dd(self.bank.inertia.velocity(delta_carrier)?, delta_carrier) * dd(0.5);
        for i in 0..N {
            let delta = impulse.axial_momentum_change_nms[i];
            loss = loss + dd(0.5 * delta / self.bank.rotors[i].axial_inertia_kg_m2) * dd(delta);
        }
        let loss = scalar(loss.to_f64(), "gear engagement loss")?;
        let residual = scalar(
            (dd(impulse.kinetic_energy_change_j) + dd(loss)).to_f64(),
            "gear engagement work",
        )?;
        if loss < 0.0 {
            return Err(invalid("nonnegative gear engagement loss"));
        }
        admit_work(residual.abs(), impulse.kinetic_energy_change_j.abs() + loss)?;
        Ok(GearEngagement {
            impulse,
            constraint_impulse_nms: constraints,
            dissipated_energy_j: loss,
            work_energy_residual_j: residual,
        })
    }

    /// Apply simultaneous original shaft drives and external body impulse.
    ///
    /// Gear reactions redistribute each component's generalized drive with the
    /// full original carrier recoil. They cannot change that generalized drive
    /// or create work. A drive is allowed on any original geared shaft.
    /// The input must already satisfy the ratios; engagement is never implicit.
    pub fn apply_impulse(
        &self,
        state: &RotorState<N>,
        drive_impulse_nms: [f64; N],
        external_impulse_body: Vec3,
    ) -> Result<GearedRotorImpulse<N, E>, LieStepError> {
        self.admit_geared(state)?;
        if E == 0 {
            let impulse =
                self.bank
                    .apply_impulse(state, drive_impulse_nms, external_impulse_body)?;
            return Ok(GearedRotorImpulse {
                constraint_impulse_nms: [0.0; E],
                drive_work_j: impulse.shaft_work_j,
                reaction_work_j: 0.0,
                work_energy_residual_j: impulse.work_energy_residual_j,
                impulse,
            });
        }
        let impulse = {
            vector(external_impulse_body, "geared external impulse")?;
            let delivered_external = super::super::super::torque_kick(
                state.total_momentum_body,
                external_impulse_body,
                1.0,
            )? - state.total_momentum_body;
            let mut generalized = [Dd::ZERO; N];
            for i in 0..N {
                scalar(drive_impulse_nms[i], "geared shaft drive impulse")?;
                generalized[self.group[i]] = generalized[self.group[i]]
                    + dd(self.relative_scale[i]) * dd(drive_impulse_nms[i]);
            }
            let mut transfer = [0.0; N];
            for i in 0..N {
                transfer[i] = scalar(
                    (dd(self.drive_share[i]) * generalized[self.group[i]]).to_f64(),
                    "geared generalized drive transfer",
                )?;
            }
            let delta_omega = self
                .inertia
                .velocity(delivered_external - self.bank.rotor_momentum(&transfer)?)?;
            for (i, value) in transfer.iter_mut().enumerate() {
                *value = scalar(
                    (dd(*value) + dot_dd(self.recoil_axis[i], delta_omega)).to_f64(),
                    "geared total shaft impulse",
                )?;
            }
            self.bank
                .apply_impulse(state, transfer, external_impulse_body)?
        };
        self.admit_geared(&impulse.state_after)?;
        let constraints =
            self.constraint_impulses(&impulse.axial_momentum_change_nms, &drive_impulse_nms)?;
        let mut drive_work = [0.0; N];
        let mut reaction_work = Dd::ZERO;
        let mut work_scale = impulse.kinetic_energy_change_j.abs() + impulse.external_work_j.abs();
        for i in 0..N {
            let rate = dd(impulse.mean_relative_rate_rad_s[i]);
            drive_work[i] = scalar(
                (dd(drive_impulse_nms[i]) * rate).to_f64(),
                "geared drive work",
            )?;
            reaction_work = reaction_work
                + (dd(impulse.axial_momentum_change_nms[i]) - dd(drive_impulse_nms[i])) * rate;
            work_scale += drive_work[i].abs() + impulse.shaft_work_j[i].abs();
        }
        let reaction_work = scalar(reaction_work.to_f64(), "ideal gear reaction work")?;
        let residual = work_residual(
            impulse.kinetic_energy_change_j,
            &drive_work,
            impulse.external_work_j,
            0.0,
        )?;
        admit_work(reaction_work.abs() + residual.abs(), work_scale)?;
        Ok(GearedRotorImpulse {
            impulse,
            constraint_impulse_nms: constraints,
            drive_work_j: drive_work,
            reaction_work_j: reaction_work,
            work_energy_residual_j: residual,
        })
    }

    /// Advance free precession in the retained gear mode, with shared phases.
    ///
    /// This is the original second-order gyrostat flow on the exactly reduced
    /// kinetic form. Reactions during precession change the original absolute
    /// rotor momenta while preserving their component generalized momenta.
    /// Signed nonzero duration permits reversal. No projection loss is added.
    pub fn drift(
        &self,
        state: &RotorState<N>,
        duration_s: f64,
    ) -> Result<GearedRotorDrift<N, E>, LieStepError> {
        let mut drift = self.drift_momenta(state, duration_s)?;
        let first = self.bank.relative_rates(state)?;
        let last = self.bank.relative_rates(&drift.state_after)?;
        let mean = core::array::from_fn(|i| 0.5 * first[i] + 0.5 * last[i]);
        drift.state_after.relative_angle_rad = self.advance_phase(state, &mean, duration_s)?;
        Ok(drift)
    }

    /// Advance held drives and a body- or space-held external torque.
    ///
    /// Two simultaneous constrained half-impulses surround one retained-mode
    /// drift. Each component advances one shared relative phase, with the same
    /// endpoint-rate quadrature as drive work. Spatial torques are transformed
    /// at both actual attitudes. The native drift defect remains numerical
    /// error, separate from drive work and the ideal reactions' zero work.
    /// This fixed-axis/fixed-ratio map is symmetric and second order; timestep
    /// refinement remains necessary. There is no backlash, finite capacity,
    /// tooth contact, parameter-change work or discrete adjoint in this owner.
    pub fn step(
        &self,
        state: &RotorState<N>,
        shaft_torque_nm: [f64; N],
        external_torque: HeldTorque,
        duration_s: f64,
    ) -> Result<GearedRotorStep<N, E>, LieStepError> {
        if E == 0 {
            let motion = self
                .bank
                .step(state, shaft_torque_nm, external_torque, duration_s)?;
            return Ok(GearedRotorStep {
                constraint_impulse_nms: [0.0; E],
                drive_work_j: motion.shaft_work_j,
                reaction_work_j: 0.0,
                work_energy_residual_j: motion.work_energy_residual_j,
                motion,
            });
        }
        scalar(duration_s, "geared step duration")?;
        let half = 0.5 * duration_s;
        if duration_s == 0.0 || half == 0.0 {
            return Err(invalid("nonzero resolved geared duration"));
        }
        let mut drive = [0.0; N];
        for i in 0..N {
            drive[i] = product(shaft_torque_nm[i], half, "geared shaft torque impulse")?;
        }
        let impulse_at = |rotation: So3| -> Result<Vec3, LieStepError> {
            let value = match external_torque {
                HeldTorque::Body(value) => vector(value, "geared external body torque")?,
                HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
            };
            Ok(Vec3::new(
                product(value.x, half, "geared external impulse")?,
                product(value.y, half, "geared external impulse")?,
                product(value.z, half, "geared external impulse")?,
            ))
        };
        let first = self.apply_impulse(state, drive, impulse_at(state.rotation)?)?;
        let drift = self.drift_momenta(&first.impulse.state_after, duration_s)?;
        let second = self.apply_impulse(
            &drift.state_after,
            drive,
            impulse_at(drift.state_after.rotation)?,
        )?;
        let mean = core::array::from_fn(|i| {
            0.5 * first.impulse.mean_relative_rate_rad_s[i]
                + 0.5 * second.impulse.mean_relative_rate_rad_s[i]
        });
        let phase = self.advance_phase(state, &mean, duration_s)?;
        let mut motion = self.bank.finish_step(
            state,
            &first.impulse,
            &drift.state_after,
            &second.impulse,
            duration_s,
            Some(phase),
        )?;
        motion.drift_energy_defect_j = drift.drift_energy_defect_j;
        motion.work_energy_residual_j = work_residual(
            motion.kinetic_energy_change_j,
            &motion.shaft_work_j,
            motion.external_work_j,
            motion.drift_energy_defect_j,
        )?;
        let mut drive_work = [0.0; N];
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(first.drive_work_j[i]) + dd(second.drive_work_j[i])).to_f64(),
                "geared step drive work",
            )?;
        }
        let mut constraints = [0.0; E];
        for e in 0..E {
            constraints[e] = scalar(
                (dd(first.constraint_impulse_nms[e])
                    + dd(drift.constraint_impulse_nms[e])
                    + dd(second.constraint_impulse_nms[e]))
                .to_f64(),
                "integrated gear reaction",
            )?;
        }
        let reaction_work = scalar(
            (dd(first.reaction_work_j) + dd(second.reaction_work_j)).to_f64(),
            "geared endpoint reaction work",
        )?;
        let residual = work_residual(
            motion.kinetic_energy_change_j,
            &drive_work,
            motion.external_work_j,
            motion.drift_energy_defect_j,
        )?;
        Ok(GearedRotorStep {
            motion,
            constraint_impulse_nms: constraints,
            drive_work_j: drive_work,
            reaction_work_j: reaction_work,
            work_energy_residual_j: residual,
        })
    }

    fn drift_momenta(
        &self,
        state: &RotorState<N>,
        duration_s: f64,
    ) -> Result<GearedRotorDrift<N, E>, LieStepError> {
        self.admit_geared(state)?;
        let omega = self.bank.angular_velocity(state)?;
        // Construct the retained gyrostatic vector from the actual carrier rate:
        // a carrier exactly at rest must not acquire a rounding-sized spin from
        // subtracting two independently accumulated generalized rotor sums.
        let stored = state.total_momentum_body - self.inertia.momentum(omega)?;
        let recovered = self.inertia.velocity(state.total_momentum_body - stored)?;
        admit_vector_error(recovered - omega, omega.norm().max(recovered.norm()))?;
        let mut after = *state;
        (after.rotation, after.total_momentum_body) = self.inertia.momentum_step(
            state.rotation,
            state.total_momentum_body,
            stored,
            duration_s,
        )?;
        let delta_total = after.total_momentum_body - state.total_momentum_body;
        let delta_omega = self.inertia.velocity(delta_total)?;
        let mut requested = [0.0; N];
        let mut change = [0.0; N];
        for i in 0..N {
            requested[i] = scalar(
                dot_dd(self.recoil_axis[i], delta_omega).to_f64(),
                "gear precession reaction",
            )?;
            after.axial_momentum_nms[i] = increment(
                state.axial_momentum_nms[i],
                requested[i],
                "geared transported axial momentum",
            )?;
            change[i] = after.axial_momentum_nms[i] - state.axial_momentum_nms[i];
        }
        self.admit_geared(&after)?;
        self.bank.admit_recoil(
            self.bank.carrier_momentum(state)?,
            self.bank.carrier_momentum(&after)?,
            delta_total,
            &change,
            delta_total,
            &requested,
        )?;
        admit_vector_error(
            self.bank.angular_velocity(&after)? - (omega + delta_omega),
            omega.norm().max((omega + delta_omega).norm()),
        )?;
        let constraints = self.constraint_impulses(&change, &[0.0; N])?;
        let mean = omega + delta_omega.scale(0.5);
        let native_change = scalar(
            dot_dd(mean, delta_total).to_f64(),
            "geared native drift energy",
        )?;
        let kinetic_change = self.bank.energy_change(state, &after)?;
        let residual = scalar(
            (dd(kinetic_change) - dd(native_change)).to_f64(),
            "geared kinetic reconstruction",
        )?;
        admit_work(residual.abs(), self.transport_work_scale(state, &after)?)?;
        let momentum_residual = subtract3(
            after.rotation.rotate(after.total_momentum_body)?,
            state.rotation.rotate(state.total_momentum_body)?,
            Vec3::default(),
        );
        vector(momentum_residual, "geared spatial momentum")?;
        Ok(GearedRotorDrift {
            state_after: after,
            axial_momentum_change_nms: change,
            constraint_impulse_nms: constraints,
            kinetic_energy_change_j: kinetic_change,
            drift_energy_defect_j: native_change,
            work_energy_residual_j: residual,
            angular_momentum_residual_world: momentum_residual,
        })
    }

    fn admit_geared(&self, state: &RotorState<N>) -> Result<(), LieStepError> {
        let omega = self.bank.angular_velocity(state)?;
        for relation in self.constraints {
            let i = relation.first_rotor;
            let j = relation.second_rotor;
            let spin_i = dd(state.axial_momentum_nms[i] / self.bank.rotors[i].axial_inertia_kg_m2);
            let spin_j = dd(state.axial_momentum_nms[j] / self.bank.rotors[j].axial_inertia_kg_m2);
            let carried_i = dot_dd(self.bank.rotors[i].axis_body, omega);
            let carried_j = dot_dd(self.bank.rotors[j].axis_body, omega);
            let residual = scalar(
                ((spin_i - carried_i) - dd(relation.ratio) * (spin_j - carried_j)).to_f64(),
                "geared relative-rate constraint",
            )?;
            let scale = scalar(
                (spin_i.abs()
                    + carried_i.abs()
                    + dd(relation.ratio.abs()) * (spin_j.abs() + carried_j.abs()))
                .to_f64(),
                "gear rate resolution scale",
            )?;
            if residual.abs() > ROUNDING * scale {
                return Err(invalid("state does not satisfy ideal gear ratios"));
            }
        }
        Ok(())
    }

    fn group_rates(&self, rates: &[f64; N]) -> Result<[f64; N], LieStepError> {
        let mut grouped = [Dd::ZERO; N];
        for i in 0..N {
            grouped[self.group[i]] =
                grouped[self.group[i]] + dd(self.drive_share[i]) * dd(rates[i]);
        }
        let mut result = [0.0; N];
        for i in 0..N {
            result[i] = scalar(grouped[i].to_f64(), "shared geared phase rate")?;
        }
        Ok(result)
    }

    fn advance_phase(
        &self,
        state: &RotorState<N>,
        rates: &[f64; N],
        duration: f64,
    ) -> Result<[f64; N], LieStepError> {
        let grouped = self.group_rates(rates)?;
        let mut phase = state.relative_angle_rad;
        let mut actual = [0.0; N];
        for i in 0..N {
            let delta = product(
                duration,
                product(
                    self.relative_scale[i],
                    grouped[self.group[i]],
                    "geared phase rate",
                )?,
                "geared phase increment",
            )?;
            phase[i] = increment(phase[i], delta, "geared relative phase")?;
            actual[i] = (dd(phase[i]) - dd(state.relative_angle_rad[i])).to_f64();
        }
        for relation in self.constraints {
            let a = dd(actual[relation.first_rotor]);
            let b = dd(relation.ratio) * dd(actual[relation.second_rotor]);
            if scalar((a - b).to_f64(), "geared phase relation")?.abs()
                > ROUNDING * scalar((a.abs() + b.abs()).to_f64(), "gear phase scale")?
            {
                return Err(invalid("unresolved geared phase increments"));
            }
        }
        Ok(phase)
    }

    fn constraint_impulses(
        &self,
        actual: &[f64; N],
        drives: &[f64; N],
    ) -> Result<[f64; E], LieStepError> {
        if E == 0 {
            return Ok([0.0; E]);
        }
        let mut reaction: [Dd; N] = core::array::from_fn(|i| dd(actual[i]) - dd(drives[i]));
        let mut result = [0.0; E];
        for &node in self.order.iter().rev() {
            if self.parent[node] == node {
                continue;
            }
            let e = self.parent_edge[node];
            let relation = self.constraints[e];
            let (child, parent) = if node == relation.first_rotor {
                (1.0, -relation.ratio)
            } else {
                (-relation.ratio, 1.0)
            };
            result[e] = scalar(
                (reaction[node] / dd(child)).to_f64(),
                "resolved gear mesh impulse",
            )?;
            reaction[self.parent[node]] = reaction[self.parent[node]] - dd(parent) * dd(result[e]);
        }
        // Reconstruct every original shaft, not merely a total that could hide
        // lost transfers. No stored-spin magnitude enters this local audit.
        let mut rebuilt = [Dd::ZERO; N];
        let mut scale = [Dd::ZERO; N];
        for (relation, value) in self.constraints.iter().zip(result) {
            let first = dd(value);
            let second = -dd(relation.ratio) * dd(value);
            rebuilt[relation.first_rotor] = rebuilt[relation.first_rotor] + first;
            rebuilt[relation.second_rotor] = rebuilt[relation.second_rotor] + second;
            scale[relation.first_rotor] = scale[relation.first_rotor] + first.abs();
            scale[relation.second_rotor] = scale[relation.second_rotor] + second.abs();
        }
        for i in 0..N {
            let residual = scalar(
                (dd(actual[i]) - dd(drives[i]) - rebuilt[i]).to_f64(),
                "gear mesh reaction reconstruction",
            )?;
            let scale = scalar(
                (scale[i] + dd(actual[i]).abs() + dd(drives[i]).abs()).to_f64(),
                "gear impulse resolution scale",
            )?;
            if residual.abs() > ROUNDING * scale {
                return Err(invalid("unresolved generalized gear drive or reaction"));
            }
        }
        Ok(result)
    }

    fn transport_work_scale(
        &self,
        before: &RotorState<N>,
        after: &RotorState<N>,
    ) -> Result<f64, LieStepError> {
        let c0 = self.bank.carrier_momentum(before)?;
        let c1 = self.bank.carrier_momentum(after)?;
        let mean =
            self.bank.inertia.velocity(c0)?.scale(0.5) + self.bank.inertia.velocity(c1)?.scale(0.5);
        let mut scale = components(mean)
            .into_iter()
            .zip(components(c1 - c0))
            .map(|(rate, impulse)| (rate * impulse).abs())
            .sum::<f64>();
        for i in 0..N {
            let mean = (0.5 * before.axial_momentum_nms[i] + 0.5 * after.axial_momentum_nms[i])
                / self.bank.rotors[i].axial_inertia_kg_m2;
            scale += (mean * (after.axial_momentum_nms[i] - before.axial_momentum_nms[i])).abs();
        }
        scalar(scale, "geared transport work scale")
    }
}

fn admit_work(residual: f64, scale: f64) -> Result<(), LieStepError> {
    scalar(residual, "geared work residual")?;
    scalar(scale, "geared work scale")?;
    if residual > ROUNDING * scale {
        return Err(invalid("unresolved ideal gear work"));
    }
    Ok(())
}

fn admit_vector_error(error: Vec3, scale: f64) -> Result<(), LieStepError> {
    vector(error, "geared carrier velocity reconstruction")?;
    scalar(scale, "geared carrier velocity scale")?;
    if error.norm() > ROUNDING * scale {
        return Err(invalid("unresolved retained gear carrier rate"));
    }
    Ok(())
}
