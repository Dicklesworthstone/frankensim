//! Objective elastic history under a changed native contact projection.
//!
//! The existing tribology isometry transports shear; it does not apply a slip
//! increment or integrate the chord between old and new vector components.
//! Original normal approach, source, policy and chronology stay fixed. Native
//! geometry owns the new complete mobility and actual projected rates.
//!
//! K_parallel = v^T W^-1 v/2 is the contact-observable kinetic partition, not
//! total body kinetic energy. Changing the projection can change this number
//! without any physical impulse or energy supply. Its signed difference is a
//! diagnostic only; the native owner must independently authenticate body
//! geometry, momenta and total kinetic energy. No work or heat port is created.

use super::{
    JkrMixedModeMotion, JkrMixedModeMotionError, JkrMixedModeMotionState, JkrMixedModeResponse,
    finite, mul, norm, poll,
};
use crate::normal_patch::adhesive::mixed_mode::path::JkrMixedModePathError;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::tangential::transport::TangentialHistoryRotation;

/// A zero-mechanical-step transfer. Projected kinetic repartition is not work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeProjectionReport {
    pub before: JkrMixedModeResponse,
    pub after: JkrMixedModeResponse,
    pub rotation: TangentialHistoryRotation,
    /// Original source storage change retained as numerical roundoff [J].
    pub source_energy_residual_j: f64,
    pub source_energy_tolerance_j: f64,
    pub projected_kinetic_energy_before_j: f64,
    pub projected_kinetic_energy_after_j: f64,
    /// Signed change of the contact-observable kinetic partition [J].
    /// This is neither an energy source, a loss, heat nor a body-energy audit.
    pub projected_kinetic_repartition_j: f64,
    pub evaluations: usize,
}

/// Immutable successor proposal retaining the original material owners.
pub struct JkrMixedModeProjectionContinuation<'a> {
    previous: JkrMixedModeMotionState,
    motion: JkrMixedModeMotion<'a>,
    state: JkrMixedModeMotionState,
    report: JkrMixedModeProjectionReport,
}
impl<'a> JkrMixedModeProjectionContinuation<'a> {
    pub const fn previous_state(&self) -> JkrMixedModeMotionState {
        self.previous
    }
    pub const fn motion(&self) -> &JkrMixedModeMotion<'a> {
        &self.motion
    }
    pub const fn state(&self) -> JkrMixedModeMotionState {
        self.state
    }
    pub const fn report(&self) -> &JkrMixedModeProjectionReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrMixedModeMotion<'a>,
        JkrMixedModeMotionState,
        JkrMixedModeProjectionReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

impl<'a> JkrMixedModeMotion<'a> {
    /// Rebind the SAME accepted full-stick history to a new native projection.
    ///
    /// The existing checked tribology rotation preserves the elastic shear
    /// vector's length. A gauge change only re-expresses it; minimum-normal
    /// rotation and declared common material spin follow that owner's existing
    /// convention. Relative torsion, contact migration and reattachment are not
    /// supplied by this operation. No time interval or physical slip is added.
    ///
    /// The native caller must authenticate `mobility_per_kg` and `velocity_m_s`
    /// from actual compatible geometry and unchanged mechanical history. This
    /// source-level operation cannot establish that body-level correspondence.
    /// It preserves normal approach, source/shear claims, all budgets, time and
    /// step count; only shear coordinates and the native projection are rebound.
    /// A pure coordinate change additionally has W'=P W P^T and v'=P v, where
    /// P=diag(1,R). A physical change of projection need not preserve K_parallel.
    ///
    /// One original source query is used if shear coordinates are unchanged,
    /// otherwise two. Strict source/fold, static-capacity, force, speed and
    /// zero-storage-change gates remain active. No fictitious chord is passed
    /// through the contact law, and the physical-increment bound does not limit
    /// a coordinate re-expression. Failure changes no accepted owner or token.
    pub fn rebind_native_projection(
        &self,
        previous: JkrMixedModeMotionState,
        rotation: TangentialHistoryRotation,
        mobility_per_kg: [[f64; 3]; 3],
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<JkrMixedModeProjectionContinuation<'a>, JkrMixedModeMotionError> {
        poll(gate)?;
        if previous.identity != self.identity {
            return Err(JkrMixedModeMotionError::StateMismatch);
        }
        let shear = rotation
            .apply([previous.position_m[1], previous.position_m[2]])
            .map_err(|error| JkrMixedModePathError::Source(error.into()))?;
        let position_m = [previous.position_m[0], shear[0], shear[1]];
        let evaluations = if position_m == previous.position_m {
            1
        } else {
            2
        };
        if evaluations > self.budget.maximum_evaluations {
            return Err(JkrMixedModeMotionError::Budget);
        }
        let motion =
            JkrMixedModeMotion::new(self.source, self.shear, mobility_per_kg, self.budget)?;
        motion.admit_velocity(velocity_m_s)?;
        let before = self
            .source
            .discrete_path(previous.position_m, previous.position_m, gate)?
            .contact;
        let after = if evaluations == 1 {
            before
        } else {
            self.source
                .discrete_path(position_m, position_m, gate)?
                .contact
        };
        // Reuse the original constitutive capacity, including its real area
        // and repulsive-load semantics; no second junction law is introduced.
        let capacity = after
            .friction_contact_state()
            .and_then(|contact| self.shear.law().static_capacity_n(contact))
            .map_err(|error| JkrMixedModePathError::Source(error.into()))?;
        let shear_force = after.tangential_force_n[0].hypot(after.tangential_force_n[1]);
        if shear_force > capacity {
            return Err(JkrMixedModePathError::StaticCapacity {
                force_n: shear_force,
                capacity_n: capacity,
            }
            .into());
        }
        let force = [
            after.normal_force_n,
            after.tangential_force_n[0],
            after.tangential_force_n[1],
        ];
        if norm(force) > self.budget.maximum_force_n {
            return Err(JkrMixedModeMotionError::Range("transported source force"));
        }
        let source_energy_residual_j = finite(
            (Dd::from_f64(after.reversible_energy_j) - Dd::from_f64(before.reversible_energy_j))
                .to_f64(),
            "transported source storage difference",
        )?;
        let scale = finite(
            before.normal_elastic_energy_j + before.tangential_elastic_energy_j
                - before.surface_energy_j
                + after.normal_elastic_energy_j
                + after.tangential_elastic_energy_j
                - after.surface_energy_j,
            "transported source storage scale",
        )?;
        let source_energy_tolerance_j = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "transported source storage tolerance",
        )?;
        if source_energy_residual_j.abs() > source_energy_tolerance_j {
            return Err(JkrMixedModeMotionError::EnergyBalance {
                residual_j: source_energy_residual_j,
                tolerance_j: source_energy_tolerance_j,
            });
        }
        let projected_kinetic_energy_before_j = finite(
            0.5 * super::dot(previous.velocity_m_s, mul(self.mass, previous.velocity_m_s)),
            "old projected kinetic energy",
        )?;
        let projected_kinetic_energy_after_j = finite(
            0.5 * super::dot(velocity_m_s, mul(motion.mass, velocity_m_s)),
            "new projected kinetic energy",
        )?;
        let projected_kinetic_repartition_j = finite(
            (Dd::from_f64(projected_kinetic_energy_after_j)
                - Dd::from_f64(projected_kinetic_energy_before_j))
            .to_f64(),
            "projected kinetic repartition",
        )?;
        let state = if motion.identity == self.identity
            && position_m == previous.position_m
            && velocity_m_s == previous.velocity_m_s
        {
            previous
        } else {
            JkrMixedModeMotionState {
                identity: motion.identity,
                position_m,
                velocity_m_s,
                ..previous
            }
        };
        poll(gate)?;
        Ok(JkrMixedModeProjectionContinuation {
            previous,
            motion,
            state,
            report: JkrMixedModeProjectionReport {
                before,
                after,
                rotation,
                source_energy_residual_j,
                source_energy_tolerance_j,
                projected_kinetic_energy_before_j,
                projected_kinetic_energy_after_j,
                projected_kinetic_repartition_j,
                evaluations,
            },
        })
    }
}
