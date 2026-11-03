//! Recoiling monopole-plus-quadrupole gravity on native rigid bodies.
//!
//! One symmetric pair potential owns both COM forces and both gravity-gradient
//! torques. The spin/orbital exchange therefore closes without fixing either
//! source, inventing a torque reaction, or projecting the returned momenta.
//! Terms through second mass moment are retained; quadrupole-quadrupole and
//! higher moments are omitted. Declared COM bounding spheres must be disjoint
//! and satisfy the caller's far-field ratio throughout the split drift.

use super::super::{quadrupole_correction, quadrupole_inertia_trace};
use super::{
    GravityError, MutualGravity, MutualGravityStepReceipt, Sum, VectorSum, binary_parts,
    check_drift_separations, diagnostics, force_component, pair_magnitudes, poll, same_frame,
    scaled_positive, storage,
};
use crate::{DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench};

/// An explicitly selected extended-body model on an existing mutual-gravity owner.
///
/// Radii follow exactly the original body/property order and are borrowed, not
/// inferred from inertia. Their sum divided by COM separation must not exceed
/// the declared ratio, which is strictly below one. This excludes overlapping
/// bounds; it is not a certified multipole truncation-error estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MutualQuadrupoleGravity<'a> {
    monopole: MutualGravity,
    body_radii_m: &'a [f64],
    maximum_radius_ratio: f64,
}

/// One simultaneous conservative field sample in the supplied body order.
#[derive(Clone, Debug, PartialEq)]
pub struct MutualQuadrupoleLoad {
    /// Equal/opposite pair forces in world coordinates and each body's own
    /// gravity-gradient torque about COM in its actual body coordinates.
    pub wrenches: Vec<Wrench>,
    /// Pair potential counted once, including both second-moment corrections [J].
    pub potential_energy_j: f64,
    /// Number of unordered pairs evaluated in this sample.
    pub pair_count: usize,
}

impl MutualGravity {
    /// Select finite-source quadrupole interactions using original native inertia.
    ///
    /// Radius values, body counts and mass-radius consistency are checked by
    /// each load/step under the original pair budget and cancellation policy.
    /// A ratio nearer zero selects a stricter far field. No softening, contact
    /// geometry, shape reconstruction or higher mass moments are inferred.
    pub fn quadrupole<'a>(
        self,
        body_radii_m: &'a [f64],
        maximum_radius_ratio: f64,
    ) -> Result<MutualQuadrupoleGravity<'a>, GravityError> {
        if !maximum_radius_ratio.is_finite()
            || maximum_radius_ratio <= 0.0
            || maximum_radius_ratio >= 1.0
        {
            return Err(GravityError::Invalid(
                "mutual quadrupole radius ratio must be in (0,1)",
            ));
        }
        Ok(MutualQuadrupoleGravity {
            monopole: self,
            body_radii_m,
            maximum_radius_ratio,
        })
    }
}

impl MutualQuadrupoleGravity<'_> {
    /// Evaluate one symmetric second-moment potential and all of its derivatives.
    ///
    /// For d = x_j - x_i, n = d/r and Q = m_j I_i + m_i I_j in world axes,
    /// U = -G m_i m_j/r + G (3 n.Q.n - tr Q)/(2 r^3).
    /// Both COM forces and the individual body torques differentiate this same U.
    /// The quadrupole correction is the original fixed-source implementation.
    ///
    /// Work is O(N^2), scratch is O(N), and no inertia is diagonalized per pair.
    /// Compensated reductions use the supplied fixed body order. Unrepresentable
    /// coefficients or results refuse; no multipole is silently discarded.
    pub fn load(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MutualQuadrupoleLoad, GravityError> {
        poll(&mut cancelled)?;
        same_frame(frame, self.monopole.frame)?;
        if states.len() != properties.len() || states.len() != self.body_radii_m.len() {
            return Err(GravityError::Invalid(
                "mutual quadrupole body/property/radius count mismatch",
            ));
        }
        let count = states.len();
        let pair_count = self.monopole.pair_count(count)?;
        let mut force_sums = storage(count)?;
        let mut torque_sums = storage(count)?;
        let mut traces = storage(count)?;
        for ((&state, &mass), &radius) in states.iter().zip(properties).zip(self.body_radii_m) {
            poll(&mut cancelled)?;
            state.validate()?;
            mass.validate()?;
            if !radius.is_finite() || radius <= 0.0 {
                return Err(GravityError::Invalid(
                    "mutual quadrupole body radius must be positive and finite",
                ));
            }
            traces.push(quadrupole_inertia_trace(mass, radius)?);
            force_sums.push(VectorSum::default());
            torque_sums.push(VectorSum::default());
        }
        let mut potential = Sum::default();
        for i in 0..count {
            for j in i + 1..count {
                poll(&mut cancelled)?;
                let displacement = states[j]
                    .pose()
                    .position_world()
                    .sub(states[i].pose().position_world());
                let radius = self.monopole.admit_separation(displacement)?;
                self.admit_pair_radius(i, j, radius)?;
                let n = Vec3::new(
                    displacement.x / radius,
                    displacement.y / radius,
                    displacement.z / radius,
                );
                let (magnitude, force_magnitude) = pair_magnitudes(
                    self.monopole.gravitational_constant,
                    properties[i].mass(),
                    properties[j].mass(),
                    radius,
                )?;
                let first = quadrupole_correction(
                    states[i].pose(),
                    properties[i],
                    n.scale(-1.0),
                    radius,
                    inverse_cube_parameter(
                        self.monopole.gravitational_constant,
                        properties[j].mass(),
                        radius,
                    )?,
                    traces[i],
                )?;
                let second = quadrupole_correction(
                    states[j].pose(),
                    properties[j],
                    n,
                    radius,
                    inverse_cube_parameter(
                        self.monopole.gravitational_constant,
                        properties[i].mass(),
                        radius,
                    )?,
                    traces[j],
                )?;
                // Each shape correction reacts on its source. Apply the SAME
                // pair resultant with opposite signs; do not update sources
                // from independently rounded evaluations of the whole force.
                let mut pair_force = VectorSum::default();
                pair_force.add(Vec3::new(
                    force_component(force_magnitude, displacement.x, radius)?,
                    force_component(force_magnitude, displacement.y, radius)?,
                    force_component(force_magnitude, displacement.z, radius)?,
                ))?;
                pair_force.add(first.wrench.force_world)?;
                pair_force.add(second.wrench.force_world.scale(-1.0))?;
                let force = pair_force.value()?;
                force_sums[i].add(force)?;
                force_sums[j].add(force.scale(-1.0))?;
                torque_sums[i].add(first.wrench.torque_body)?;
                torque_sums[j].add(second.wrench.torque_body)?;
                potential.add(-magnitude)?;
                potential.add(first.potential_energy_j)?;
                potential.add(second.potential_energy_j)?;
            }
        }
        let mut wrenches = storage(count)?;
        for (force, torque) in force_sums.into_iter().zip(torque_sums) {
            poll(&mut cancelled)?;
            wrenches.push(Wrench {
                force_world: force.value()?,
                torque_body: torque.value()?,
            });
        }
        let potential_energy_j = potential.value()?;
        poll(&mut cancelled)?;
        Ok(MutualQuadrupoleLoad {
            wrenches,
            potential_energy_j,
            pair_count,
        })
    }

    /// Advance all bodies through the original native Verlet stages.
    ///
    /// Initial loads are simultaneous. All original half-kicks and full free
    /// drifts precede the second simultaneous field sample and final half-kicks.
    /// Force and torque use the pose of their own endpoint. This conservative
    /// splitting is second order and symmetric at fixed resolved steps; its
    /// energy error is not zero. Linear and orbital-plus-spin angular momenta
    /// close to roundoff without an after-the-fact correction.
    ///
    /// Endpoints and the closest point of every straight relative split drift
    /// satisfy both the original separation cutoff and the combined-radius
    /// far-field bound. This is not an enclosure of the true curved trajectory.
    /// Close encounters require caller refinement. Cancellation, invalid
    /// arithmetic or any domain/budget refusal publishes no partial system.
    pub fn step(
        self,
        frame: u64,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        duration_seconds: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MutualGravityStepReceipt, GravityError> {
        if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration.into());
        }
        let before = self.load(frame, states, properties, &mut cancelled)?;
        let diagnostics_before = diagnostics(
            states,
            properties,
            before.potential_energy_j,
            &mut cancelled,
        )?;
        let native = RigidBodyIntegrator::new(Gravity::ZERO);
        let mut stages = storage(states.len())?;
        let mut candidates = storage(states.len())?;
        for ((&state, &mass), &wrench) in states.iter().zip(properties).zip(&before.wrenches) {
            poll(&mut cancelled)?;
            let stage = native.begin_verlet(state, mass, wrench, duration_seconds)?;
            candidates.push(stage.drifted_state());
            stages.push(stage);
        }
        check_drift_separations(
            states,
            &candidates,
            |i, j, radius| self.admit_pair_radius(i, j, radius),
            &mut cancelled,
        )?;
        let after = self.load(frame, &candidates, properties, &mut cancelled)?;
        for ((state, stage), &wrench) in candidates.iter_mut().zip(stages).zip(&after.wrenches) {
            poll(&mut cancelled)?;
            *state = stage.finish(wrench)?.state_after;
        }
        let diagnostics_after = diagnostics(
            &candidates,
            properties,
            after.potential_energy_j,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(MutualGravityStepReceipt {
            duration_seconds,
            states_after: candidates,
            diagnostics_before,
            diagnostics_after,
            pair_count: before.pair_count,
        })
    }

    fn admit_pair_radius(self, i: usize, j: usize, radius: f64) -> Result<(), GravityError> {
        self.monopole.admit_radius(radius)?;
        // Divide before summing so finite radii need not form a representable sum.
        let ratio = self.body_radii_m[i] / radius + self.body_radii_m[j] / radius;
        if !ratio.is_finite() || ratio > self.maximum_radius_ratio {
            Err(GravityError::OutsideExpansionDomain)
        } else {
            Ok(())
        }
    }
}

// Preserve the original pair coefficient's binary scaling rather than form
// G*m or r^3 first. A nonzero but unrepresentable coefficient refuses explicitly.
fn inverse_cube_parameter(g: f64, mass: f64, radius: f64) -> Result<f64, GravityError> {
    let (g, eg) = binary_parts(g);
    let (mass, em) = binary_parts(mass);
    let (r, er) = binary_parts(radius);
    scaled_positive(((g * mass) / r) / r / r, eg + em - 3 * er)
}
