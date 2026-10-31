//! Impulsive loading at rotating material attachments on an existing rod.
//!
//! A batch acts at one frozen pose: forces, elastic strain and the clock do not
//! advance. The original central inertia and supported-DOF policies determine
//! the velocity jump. Point arms remain relative to the material node, even
//! when an eccentric attachment has moved the nodal mass centroid.

use super::{
    RodMotion, RodMotionError, RodMotionInertia, RodMotionObservables, RodSupportImpulse, add,
    checkpoint, cross, finite, rotate, sub,
};
use fs_exec::CancelGate;

/// An instantaneous force impulse at a point fixed in one material frame.
/// The optional angular impulse is a free couple, independent of the point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodPointImpulse {
    pub node: usize,
    /// Offset FROM THE MATERIAL NODE in its current material axes [m].
    pub point_body_m: [f64; 3],
    /// Linear impulse ON the rod, world axes [N s].
    pub linear_world_n_s: [f64; 3],
    /// Additional pure angular impulse, world axes [N m s].
    pub couple_world_nm_s: [f64; 3],
}

/// Absolute acceptance bounds for the ACTUALLY stored momentum increments.
/// All tolerances must be finite and strictly positive. Each bound applies
/// independently at every node and to the complete batch, preventing opposite
/// unresolved transfers from cancelling in a global account.
#[derive(Clone, Copy, Debug)]
pub struct RodImpulseControls {
    pub balance_origin_m: [f64; 3],
    /// Euclidean linear-impulse residual bound [N s].
    pub linear_tolerance_n_s: f64,
    /// Euclidean central-spin and world angular-impulse residual bound [N m s].
    pub angular_tolerance_nm_s: f64,
    /// Absolute kinetic-increment minus applied-work residual bound [J].
    pub work_tolerance_j: f64,
}

/// Actual kinetic and momentum changes at the unchanged configuration.
#[derive(Clone, Debug)]
pub struct RodImpulseReport {
    pub before: RodMotionObservables,
    pub after: RodMotionObservables,
    /// Impulse dotted with the mean pre/post point velocity, plus couple work.
    /// This is mechanical impulse work at a frozen mass matrix, not heat.
    pub applied_work_j: f64,
    /// Stable actual increment from momentum changes and mean velocities;
    /// subtracting the large `after` and `before` energies can lose this value.
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub actual_linear_momentum_change_kg_m_s: [f64; 3],
    pub actual_angular_momentum_change_kg_m2_s: [f64; 3],
    /// Caller-applied impulses, excluding the separately returned supports.
    pub applied_linear_impulse_n_s: [f64; 3],
    /// Caller-applied moment about the supplied world balance origin.
    pub applied_angular_impulse_nm_s: [f64; 3],
    pub linear_balance_residual_n_s: [f64; 3],
    pub angular_balance_residual_nm_s: [f64; 3],
    /// Stationary support impulses ON the rod, one entry per material node.
    pub support_impulses: Vec<RodSupportImpulse>,
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Apply simultaneous point impulses without advancing elastic motion.
    ///
    /// Each point force also changes central spin by `(r - c) x J`. Full
    /// eccentric section inertia is consumed through the same native angular
    /// response as time stepping. Fixed translations/frames receive the exact
    /// balancing support impulses and do no work. Previously unsupported
    /// eccentric moving-rotor position constraints remain inadmissible.
    ///
    /// Duplicate nodes are accumulated before updating velocity, so the work
    /// includes cross terms between simultaneous impulses. Cost is O(nodes +
    /// impulses), with O(nodes) scratch; no source evaluation, factorization or
    /// eigen solve is performed. Input order fixes the deterministic reduction.
    /// Per-node stored momentum changes and mean velocities close the accepted
    /// impulse and work; a sub-resolution kick against large existing momentum
    /// refuses when it exceeds the caller's absolute transfer tolerances.
    /// Failure or cancellation returns no candidate and never changes `self`.
    /// Pose, source inventory, cached response, time and step count are retained.
    pub fn apply_point_impulses(
        &self,
        impulses: &[RodPointImpulse],
        controls: RodImpulseControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodImpulseReport), RodMotionError> {
        checkpoint(gate)?;
        if [
            controls.linear_tolerance_n_s,
            controls.angular_tolerance_nm_s,
            controls.work_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(RodMotionError::Invalid(
                "finite positive rod impulse transfer tolerances required",
            ));
        }
        let balance_origin_m = controls.balance_origin_m;
        let before = self.observe(balance_origin_m)?;
        let count = self.kinetics.len();
        let mut linear = vec![[0.0; 3]; count];
        let mut couple_body = vec![[0.0; 3]; count];
        for impulse in impulses {
            checkpoint(gate)?;
            if impulse.node >= count {
                return Err(RodMotionError::Invalid(
                    "point impulse node is outside the rod",
                ));
            }
            finite(&impulse.point_body_m)?;
            finite(&impulse.linear_world_n_s)?;
            finite(&impulse.couple_world_nm_s)?;
            let node = impulse.node;
            let inverse = self.pose.frames[node].inverse();
            let point_moment = cross(
                impulse.point_body_m,
                rotate(inverse, impulse.linear_world_n_s)?,
            );
            linear[node] = add(linear[node], impulse.linear_world_n_s);
            couple_body[node] = add(
                couple_body[node],
                add(point_moment, rotate(inverse, impulse.couple_world_nm_s)?),
            );
            finite(&linear[node])?;
            finite(&couple_body[node])?;
        }

        let mut candidate = self.clone();
        let mut work = 0.0;
        let mut kinetic_change = 0.0;
        let mut energy_residual = 0.0;
        let mut applied_linear = [0.0; 3];
        let mut applied_angular = [0.0; 3];
        let mut actual_linear = [0.0; 3];
        let mut actual_angular = [0.0; 3];
        let mut linear_residual = [0.0; 3];
        let mut angular_residual = [0.0; 3];
        let mut linear_defect: f64 = 0.0;
        let mut angular_defect: f64 = 0.0;
        let mut work_defect: f64 = 0.0;
        let mut supports = Vec::with_capacity(count);
        for node in 0..count {
            checkpoint(gate)?;
            let k = self.kinetics[node];
            let frame = self.pose.frames[node];
            let central = sub(
                couple_body[node],
                cross(k.center, rotate(frame.inverse(), linear[node])?),
            );
            finite(&central)?;
            let old_omega = k.omega(self.spin[node])?;
            let mut support_f = [0.0; 3];
            let mut support_t = [0.0; 3];
            for c in 0..3 {
                if self.model.fixed_positions()[node][c] {
                    support_f[c] = -linear[node][c];
                } else {
                    candidate.linear[node][c] += linear[node][c];
                }
                if self.model.fixed_frames()[node] {
                    support_t[c] = -central[c];
                } else {
                    candidate.spin[node][c] += central[c];
                }
            }
            let delta_p = sub(candidate.linear[node], self.linear[node]);
            let delta_spin = sub(candidate.spin[node], self.spin[node]);
            let linear_error = sub(delta_p, add(linear[node], support_f));
            let spin_error = sub(delta_spin, add(central, support_t));
            finite(&delta_p)?;
            finite(&delta_spin)?;
            finite(&linear_error)?;
            finite(&spin_error)?;
            // A fixed material frame can translate. Shift the support couple
            // back to that node, including the constrained linear impulse.
            support_t = add(
                support_t,
                cross(k.center, rotate(frame.inverse(), support_f)?),
            );
            let omega = k.omega(candidate.spin[node])?;
            let mut node_work_residual = 0.0;
            for c in 0..3 {
                let mean_velocity = 0.5 * (self.linear[node][c] / k.mass)
                    + 0.5 * (candidate.linear[node][c] / k.mass);
                let mean_omega = 0.5 * old_omega[c] + 0.5 * omega[c];
                work += linear[node][c] * mean_velocity + central[c] * mean_omega;
                kinetic_change += delta_p[c] * mean_velocity + delta_spin[c] * mean_omega;
                // Subtract impulses before multiplying by velocity, so the
                // work gate never subtracts two large total kinetic energies.
                node_work_residual += (delta_p[c] - linear[node][c]) * mean_velocity
                    + (delta_spin[c] - central[c]) * mean_omega;
            }
            let support_t = rotate(frame, support_t)?;
            let couple_world = rotate(frame, couple_body[node])?;
            let arm = sub(self.pose.positions_m[node], balance_origin_m);
            let moment = add(cross(arm, linear[node]), couple_world);
            applied_linear = add(applied_linear, linear[node]);
            applied_angular = add(applied_angular, moment);
            let centroid_arm = add(arm, rotate(frame, k.center)?);
            let angular_change = add(cross(centroid_arm, delta_p), rotate(frame, delta_spin)?);
            let angular_error = add(
                cross(centroid_arm, linear_error),
                rotate(frame, spin_error)?,
            );
            finite(&angular_change)?;
            finite(&angular_error)?;
            finite(&[node_work_residual])?;
            actual_linear = add(actual_linear, delta_p);
            actual_angular = add(actual_angular, angular_change);
            linear_residual = add(linear_residual, linear_error);
            angular_residual = add(angular_residual, angular_error);
            energy_residual += node_work_residual;
            linear_defect = linear_defect.max(norm(linear_error));
            angular_defect = angular_defect
                .max(norm(angular_error))
                .max(norm(spin_error));
            work_defect = work_defect.max(node_work_residual.abs());
            finite(&support_f)?;
            finite(&support_t)?;
            supports.push(RodSupportImpulse {
                linear_n_s: support_f,
                angular_nm_s: support_t,
            });
        }
        let after = candidate.observe(balance_origin_m)?;
        finite(&[work, kinetic_change, energy_residual])?;
        finite(&applied_linear)?;
        finite(&applied_angular)?;
        finite(&actual_linear)?;
        finite(&actual_angular)?;
        finite(&linear_residual)?;
        finite(&angular_residual)?;
        linear_defect = linear_defect.max(norm(linear_residual));
        angular_defect = angular_defect.max(norm(angular_residual));
        work_defect = work_defect.max(energy_residual.abs());
        if linear_defect > controls.linear_tolerance_n_s
            || angular_defect > controls.angular_tolerance_nm_s
            || work_defect > controls.work_tolerance_j
        {
            return Err(RodMotionError::ImpulseBalance {
                linear_defect_n_s: linear_defect,
                angular_defect_nm_s: angular_defect,
                work_defect_j: work_defect,
            });
        }
        checkpoint(gate)?;
        Ok((
            candidate,
            RodImpulseReport {
                before,
                after,
                applied_work_j: work,
                kinetic_energy_change_j: kinetic_change,
                energy_residual_j: energy_residual,
                actual_linear_momentum_change_kg_m_s: actual_linear,
                actual_angular_momentum_change_kg_m2_s: actual_angular,
                applied_linear_impulse_n_s: applied_linear,
                applied_angular_impulse_nm_s: applied_angular,
                linear_balance_residual_n_s: linear_residual,
                angular_balance_residual_nm_s: angular_residual,
                support_impulses: supports,
            },
        ))
    }
}

fn norm(value: [f64; 3]) -> f64 {
    fs_math::det::hypot(fs_math::det::hypot(value[0], value[1]), value[2])
}
