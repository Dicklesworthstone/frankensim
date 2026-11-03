//! Impulsive loading at rotating material attachments on an existing rod.
//!
//! A batch acts at one frozen pose: forces, elastic strain and the clock do not
//! advance. The original central inertia and supported-DOF policies determine
//! the velocity jump. Point arms remain relative to the material node, even
//! when an eccentric attachment has moved the nodal mass centroid.

use super::{
    RodMotion, RodMotionError, RodMotionInertia, RodMotionObservables, RodSupportImpulse, add,
    checkpoint, finite, rotate, sub,
};
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_math::dd::Dd;

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
    /// Duplicate nodes are accumulated in double-double arithmetic before one
    /// final rounding into native momenta. Point moments and free couples keep
    /// their small remainder when large opposite impulses cancel; the work
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
        let mut linear = vec![[Dd::ZERO; 3]; count];
        let mut couple_body = vec![[Dd::ZERO; 3]; count];
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
            let point_moment = cross_dd(
                impulse.point_body_m.map(dd),
                rotate(inverse, impulse.linear_world_n_s)?.map(dd),
            );
            let free_couple = rotate(inverse, impulse.couple_world_nm_s)?;
            for c in 0..3 {
                linear[node][c] = linear[node][c] + dd(impulse.linear_world_n_s[c]);
                couple_body[node][c] = couple_body[node][c] + point_moment[c] + dd(free_couple[c]);
            }
            finite_dd(&linear[node])?;
            finite_dd(&couple_body[node])?;
        }

        let mut candidate = self.clone();
        let mut work = Dd::ZERO;
        let mut kinetic_change = Dd::ZERO;
        let mut energy_residual = Dd::ZERO;
        let mut applied_linear = [Dd::ZERO; 3];
        let mut applied_angular = [Dd::ZERO; 3];
        let mut actual_linear = [Dd::ZERO; 3];
        let mut actual_angular = [Dd::ZERO; 3];
        let mut linear_residual = [Dd::ZERO; 3];
        let mut angular_residual = [Dd::ZERO; 3];
        let mut linear_defect: f64 = 0.0;
        let mut angular_defect: f64 = 0.0;
        let mut work_defect: f64 = 0.0;
        let mut supports = Vec::with_capacity(count);
        for node in 0..count {
            checkpoint(gate)?;
            let k = self.kinetics[node];
            let frame = self.pose.frames[node];
            let centroid_moment =
                cross_dd(k.center.map(dd), rotate_dd(frame.inverse(), linear[node])?);
            let central: [Dd; 3] =
                std::array::from_fn(|c| couple_body[node][c] - centroid_moment[c]);
            finite_dd(&central)?;
            let old_omega = k.omega(self.spin[node])?;
            let mut support_f = [0.0; 3];
            let mut support_t = [0.0; 3];
            for c in 0..3 {
                if self.model.fixed_positions()[node][c] {
                    support_f[c] = (-linear[node][c]).to_f64();
                } else {
                    candidate.linear[node][c] =
                        (dd(self.linear[node][c]) + linear[node][c]).to_f64();
                }
                if self.model.fixed_frames()[node] {
                    support_t[c] = (-central[c]).to_f64();
                } else {
                    candidate.spin[node][c] = (dd(self.spin[node][c]) + central[c]).to_f64();
                }
            }
            let delta_p = difference(candidate.linear[node], self.linear[node]);
            let delta_spin = difference(candidate.spin[node], self.spin[node]);
            let linear_error: [Dd; 3] =
                std::array::from_fn(|c| delta_p[c] - linear[node][c] - dd(support_f[c]));
            let spin_error: [Dd; 3] =
                std::array::from_fn(|c| delta_spin[c] - central[c] - dd(support_t[c]));
            finite_dd(&delta_p)?;
            finite_dd(&delta_spin)?;
            finite_dd(&linear_error)?;
            finite_dd(&spin_error)?;
            // A fixed material frame can translate. Shift the support couple
            // back to that node, including the constrained linear impulse.
            let support_shift = cross_dd(
                k.center.map(dd),
                rotate_dd(frame.inverse(), support_f.map(dd))?,
            );
            let support_t = rotate_dd(
                frame,
                std::array::from_fn(|c| dd(support_t[c]) + support_shift[c]),
            )?
            .map(Dd::to_f64);
            let omega = k.omega(candidate.spin[node])?;
            let mut node_work_residual = Dd::ZERO;
            for c in 0..3 {
                let mean_velocity = (dd(self.linear[node][c]) * dd(0.5)
                    + dd(candidate.linear[node][c]) * dd(0.5))
                    / dd(k.mass);
                let mean_omega = dd(old_omega[c]) * dd(0.5) + dd(omega[c]) * dd(0.5);
                work = work + linear[node][c] * mean_velocity + central[c] * mean_omega;
                kinetic_change =
                    kinetic_change + delta_p[c] * mean_velocity + delta_spin[c] * mean_omega;
                // Subtract impulses before multiplying by velocity, so the
                // work gate never subtracts two large total kinetic energies.
                node_work_residual = node_work_residual
                    + (delta_p[c] - linear[node][c]) * mean_velocity
                    + (delta_spin[c] - central[c]) * mean_omega;
            }
            let couple_world = rotate_dd(frame, couple_body[node])?;
            let arm = sub(self.pose.positions_m[node], balance_origin_m);
            let orbital_moment = cross_dd(arm.map(dd), linear[node]);
            let moment: [Dd; 3] = std::array::from_fn(|c| orbital_moment[c] + couple_world[c]);
            let centroid_arm = add(arm, rotate(frame, k.center)?);
            let orbital_change = cross_dd(centroid_arm.map(dd), delta_p);
            let spin_change = rotate_dd(frame, delta_spin)?;
            let angular_change: [Dd; 3] =
                std::array::from_fn(|c| orbital_change[c] + spin_change[c]);
            let support_orbital = cross_dd(arm.map(dd), support_f.map(dd));
            let angular_error: [Dd; 3] = std::array::from_fn(|c| {
                angular_change[c] - moment[c] - support_orbital[c] - dd(support_t[c])
            });
            finite_dd(&angular_change)?;
            finite_dd(&angular_error)?;
            finite_dd(&[node_work_residual])?;
            for c in 0..3 {
                applied_linear[c] = applied_linear[c] + linear[node][c];
                applied_angular[c] = applied_angular[c] + moment[c];
                actual_linear[c] = actual_linear[c] + delta_p[c];
                actual_angular[c] = actual_angular[c] + angular_change[c];
                linear_residual[c] = linear_residual[c] + linear_error[c];
                angular_residual[c] = angular_residual[c] + angular_error[c];
            }
            energy_residual = energy_residual + node_work_residual;
            linear_defect = linear_defect.max(norm(linear_error.map(Dd::to_f64)));
            angular_defect = angular_defect
                .max(norm(angular_error.map(Dd::to_f64)))
                .max(norm(spin_error.map(Dd::to_f64)));
            work_defect = work_defect.max(node_work_residual.to_f64().abs());
            finite(&support_f)?;
            finite(&support_t)?;
            supports.push(RodSupportImpulse {
                linear_n_s: support_f,
                angular_nm_s: support_t,
            });
        }
        let after = candidate.observe(balance_origin_m)?;
        finite_dd(&[work, kinetic_change, energy_residual])?;
        finite_dd(&applied_linear)?;
        finite_dd(&applied_angular)?;
        finite_dd(&actual_linear)?;
        finite_dd(&actual_angular)?;
        finite_dd(&linear_residual)?;
        finite_dd(&angular_residual)?;
        let work = work.to_f64();
        let kinetic_change = kinetic_change.to_f64();
        let energy_residual = energy_residual.to_f64();
        let applied_linear = applied_linear.map(Dd::to_f64);
        let applied_angular = applied_angular.map(Dd::to_f64);
        let actual_linear = actual_linear.map(Dd::to_f64);
        let actual_angular = actual_angular.map(Dd::to_f64);
        let linear_residual = linear_residual.map(Dd::to_f64);
        let angular_residual = angular_residual.map(Dd::to_f64);
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

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

fn finite_dd(values: &[Dd]) -> Result<(), RodMotionError> {
    for value in values {
        finite(&[value.hi, value.lo])?;
    }
    Ok(())
}

fn difference(a: [f64; 3], b: [f64; 3]) -> [Dd; 3] {
    std::array::from_fn(|c| dd(a[c]) - dd(b[c]))
}

fn cross_dd(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Keep the accumulated low part through the ORIGINAL rigid rotation action.
/// Ordinary f64 vectors need only the original one rotation.
fn rotate_dd(rotation: So3, value: [Dd; 3]) -> Result<[Dd; 3], RodMotionError> {
    let high = rotate(rotation, value.map(|x| x.hi))?;
    if value.iter().all(|x| x.lo == 0.0) {
        return Ok(high.map(dd));
    }
    let low = rotate(rotation, value.map(|x| x.lo))?;
    let result = std::array::from_fn(|c| dd(high[c]) + dd(low[c]));
    finite_dd(&result)?;
    Ok(result)
}

fn norm(value: [f64; 3]) -> f64 {
    fs_math::det::hypot(fs_math::det::hypot(value[0], value[1]), value[2])
}
