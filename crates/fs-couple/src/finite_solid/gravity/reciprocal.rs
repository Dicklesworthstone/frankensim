//! A frozen-configuration Newtonian kick with finite-source recoil.
//!
//! Native CentralGravity and the original solid quadrature supply one discrete
//! pair potential and its forces. The full consistent mass changes the solid's
//! velocity; the native rigid momentum owner changes the source's translation.
//! The source is a monopole at its COM, so its spin is unchanged. This is a
//! Hamiltonian force subflow, not a complete time integrator or a claim that
//! kinetic plus potential energy is conserved during an isolated kick.

/// Moving solid/source trajectories through native symmetric subflows.
pub mod trajectory;

use super::{BoundGravity, StationaryGravity, checkpoint, dd};
use crate::finite_solid::impulse::{FiniteSolidImpulse, FiniteSolidImpulseControls};
use crate::finite_solid::{Cx, FiniteSolidBody, FiniteSolidError, inflate};
use fs_math::dd::Dd;
use fs_mbd::gravity::CentralGravity;
use fs_mbd::{MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};

/// Declared physical exposure and the original native impulse allowances.
#[derive(Clone, Copy, Debug)]
pub struct ReciprocalGravityKickControls {
    /// Both configurations use these same inertial Cartesian axes and meters.
    pub frame: u64,
    /// Positive Newtonian constant; no planet or unit-system preset is inferred.
    pub gravitational_constant_m3_kg_s2: f64,
    /// Positive kick duration [s]; neither original physical clock advances.
    pub duration_s: f64,
    pub balance_origin_m: [f64; 3],
    /// Also bounds each source residual and the complete pair independently.
    pub impulse: FiniteSolidImpulseControls,
}

/// Actual native responses at the unchanged geometry and physical time.
#[derive(Clone, Debug)]
pub struct ReciprocalGravityKick {
    pub frame: u64,
    pub duration_s: f64,
    pub balance_origin_m: [f64; 3],
    /// One pair potential from the native reference-mass quadrature, zero at
    /// infinity. It is unchanged by this momentum-only subflow; count it once.
    pub potential_energy_j: f64,
    pub solid: FiniteSolidImpulse,
    pub source: MomentumImpulseReceipt,
    /// Pair momentum change MINUS the solid's stationary support impulse.
    pub linear_residual_kg_m_s: [f64; 3],
    /// Includes source orbital momentum about the declared balance origin.
    pub angular_residual_kg_m2_s: [f64; 3],
    pub kinetic_energy_change_j: f64,
    pub impulse_work_j: f64,
    /// Actual native kinetic change minus mechanical impulse work, never heat.
    pub work_residual_j: f64,
}

impl<'a> FiniteSolidBody<'a> {
    /// Apply Newtonian tidal loading and the equal/opposite source recoil.
    ///
    /// The original cellwise or nodal reference density supplies gravity and
    /// inertia. Native quadrature evaluates the CURRENT deformed coordinates;
    /// no force is evaluated only at the solid COM. The source's full mass is
    /// its original `properties.mass()`, with a point-monopole field at its COM.
    /// There is no source spin torque, solid self-gravity or finite-source shape
    /// multipole. Near-source spatial quadrature still requires refinement.
    ///
    /// This additional interaction is not retained as a stationary field and
    /// must not also be supplied through the existing external-gravity binding.
    /// Other bound external fields and material forces are unchanged and do not
    /// act in this kick. Pose, elastic state, acceleration, time and native
    /// history are preserved. The caller composes and refines physical subflows.
    /// Every refusal/cancellation leaves both inputs unchanged.
    pub fn reciprocal_gravity_kick_candidate(
        &self,
        source: RigidBodyState,
        properties: MassProperties,
        controls: ReciprocalGravityKickControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, RigidBodyState, ReciprocalGravityKick), FiniteSolidError> {
        checkpoint(cx)?;
        if !controls.duration_s.is_finite()
            || controls.duration_s <= 0.0
            || !controls.gravitational_constant_m3_kg_s2.is_finite()
            || controls.gravitational_constant_m3_kg_s2 <= 0.0
            || controls.balance_origin_m.iter().any(|x| !x.is_finite())
            || [
                controls.impulse.linear_tolerance_kg_m_s,
                controls.impulse.angular_tolerance_kg_m2_s,
                controls.impulse.work_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(FiniteSolidError::Invalid(
                "finite positive gravity kick data and tolerances required",
            ));
        }
        let field = CentralGravity::new(
            controls.frame,
            source.pose().position_world(),
            controls.gravitational_constant_m3_kg_s2 * properties.mass(),
        )
        .map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        let gravity = BoundGravity::new(
            &self.problem,
            self.reference_density_field.as_deref(),
            StationaryGravity::new(&field, controls.frame),
            cx,
        )?;
        let displacement = inflate(&self.state.q, &self.free, &self.fixed);
        let (force, potential_energy_j) =
            gravity.force_and_potential(&self.prepared, &displacement, cx)?;
        let mut total = [Dd::ZERO; 3];
        let mut impulses = Vec::with_capacity(force.len() / 3);
        for (node, force) in force.chunks_exact(3).enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            let impulse = std::array::from_fn(|c| controls.duration_s * force[c]);
            if impulse.iter().any(|x| !x.is_finite()) {
                return Err(FiniteSolidError::Invalid(
                    "gravity impulse is not representable",
                ));
            }
            for c in 0..3 {
                total[c] = total[c] + dd(impulse[c]);
            }
            impulses.push((node, impulse));
        }
        let reaction = total.map(|x| (-x).to_f64());
        let source_receipt = source
            .apply_momentum_impulse(properties, vector(reaction), Vec3::ZERO)
            .map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        // Check the source separately: opposite unresolved increments cannot
        // disappear in a joint residual. The solid owner gates its own update.
        check(
            components(source_receipt.linear_residual_world),
            components(source_receipt.angular_residual_body),
            source_receipt.work.work_energy_residual_j,
            controls.impulse,
        )?;
        let (candidate, solid) = self.nodal_impulse_candidate(
            &impulses,
            controls.balance_origin_m,
            controls.impulse,
            cx,
        )?;
        let source_increment = components(source_receipt.linear_increment_world).map(dd);
        let arm = std::array::from_fn(|c| {
            dd(components(source.pose().position_world())[c]) - dd(controls.balance_origin_m[c])
        });
        let source_orbital = cross(arm, source_increment);
        // Native residuals were measured from M*(v_after-v_before), so these
        // identities retain actual solid increments without subtracting large
        // endpoint momenta or repeating a full mass-matrix multiplication.
        let linear = std::array::from_fn(|c| {
            (dd(solid.applied.linear_kg_m_s[c])
                + dd(solid.linear_residual_kg_m_s[c])
                + source_increment[c])
                .to_f64()
        });
        let angular = std::array::from_fn(|c| {
            (dd(solid.applied.angular_kg_m2_s[c])
                + dd(solid.angular_residual_kg_m2_s[c])
                + source_orbital[c])
                .to_f64()
        });
        let work_residual_j =
            (dd(solid.work_residual_j) + dd(source_receipt.work.work_energy_residual_j)).to_f64();
        check(linear, angular, work_residual_j, controls.impulse)?;
        let kinetic_energy_change_j =
            (dd(solid.kinetic_change_j) + dd(source_receipt.work.kinetic_energy_change_j)).to_f64();
        let impulse_work_j =
            (dd(solid.impulse_work_j) + dd(source_receipt.work.impulse_work_j)).to_f64();
        if !kinetic_energy_change_j.is_finite() || !impulse_work_j.is_finite() {
            return Err(FiniteSolidError::Invalid(
                "combined gravity kick work is not representable",
            ));
        }
        checkpoint(cx)?;
        Ok((
            candidate,
            source_receipt.state_after,
            ReciprocalGravityKick {
                frame: controls.frame,
                duration_s: controls.duration_s,
                balance_origin_m: controls.balance_origin_m,
                potential_energy_j,
                solid,
                source: source_receipt,
                linear_residual_kg_m_s: linear,
                angular_residual_kg_m2_s: angular,
                kinetic_energy_change_j,
                impulse_work_j,
                work_residual_j,
            },
        ))
    }
}

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn check(
    linear: [f64; 3],
    angular: [f64; 3],
    work: f64,
    controls: FiniteSolidImpulseControls,
) -> Result<(), FiniteSolidError> {
    let norm = |a: [f64; 3]| a[0].hypot(a[1]).hypot(a[2]);
    if linear.iter().chain(&angular).any(|x| !x.is_finite())
        || !work.is_finite()
        || norm(linear) > controls.linear_tolerance_kg_m_s
        || norm(angular) > controls.angular_tolerance_kg_m2_s
        || work.abs() > controls.work_tolerance_j
    {
        Err(FiniteSolidError::ImpulseBalance {
            linear_residual_kg_m_s: linear,
            angular_residual_kg_m2_s: angular,
            work_residual_j: work,
        })
    } else {
        Ok(())
    }
}
