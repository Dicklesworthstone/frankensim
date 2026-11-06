//! Bounded held-world couples around the original native free-rotation map.

use super::{SphereContactError, poll};
use crate::spatial_loads::{PreparedSpatialLoads, SpatialLoadError};
use crate::{Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Independent rotation resolution for applied world couples.
///
/// Centre paths and impact times retain the original analytic queries. These
/// bounds resolve full-inertia spin and attitude between the detected events;
/// a small work or momentum residual is not an angular-phase error bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereTorqueControls {
    /// Positive finite maximum duration of each applied-torque rotation step [s].
    pub maximum_step_s: f64,
    /// Total native torque steps across all bodies, events and root-bracket
    /// trials. Must be in 1..=1_000_000. Zero couples retain the original drift
    /// and consume no steps; all active couples share this work budget.
    pub maximum_substeps: usize,
}

impl Default for SphereTorqueControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 0.01,
            maximum_substeps: 65_536,
        }
    }
}

pub(super) struct TorqueFlight {
    loads: Vec<Option<PreparedSpatialLoads>>,
    total_torque: [Dd; 3],
    controls: SphereTorqueControls,
    pub(super) substeps: usize,
}

impl TorqueFlight {
    pub(super) fn new(
        torques: &[Vec3],
        count: usize,
        controls: SphereTorqueControls,
        gate: &CancelGate,
    ) -> Result<Self, SphereContactError> {
        if torques.len() != count
            || torques.iter().any(|torque| !torque.is_finite())
            || !controls.maximum_step_s.is_finite()
            || controls.maximum_step_s <= 0.0
            || controls.maximum_substeps == 0
            || controls.maximum_substeps > 1_000_000
        {
            return Err(SphereContactError::Invalid(
                "one finite world couple per body and bounded torque controls required",
            ));
        }
        let mut loads = Vec::with_capacity(count);
        let mut total_torque = [Dd::ZERO; 3];
        for &torque in torques {
            poll(gate)?;
            for (sum, value) in total_torque.iter_mut().zip([torque.x, torque.y, torque.z]) {
                *sum = *sum + Dd::from_f64(value);
            }
            let prepared = if torque == Vec3::ZERO {
                None
            } else {
                Some(
                    PreparedSpatialLoads::new(Wrench::ZERO, torque, &[], || gate.is_requested())
                        .map_err(|error| match error {
                            SpatialLoadError::Dynamics(error) => {
                                SphereContactError::Dynamics(error)
                            }
                            SpatialLoadError::Cancelled => SphereContactError::Cancelled,
                            SpatialLoadError::TooManyPointForces => SphereContactError::Invalid(
                                "unexpected spatial-couple preparation refusal",
                            ),
                        })?,
                )
            };
            loads.push(prepared);
        }
        Ok(Self {
            loads,
            total_torque,
            controls,
            substeps: 0,
        })
    }

    pub(super) fn impulse(&self, duration: f64) -> Result<Vec3, SphereContactError> {
        Ok(crate::verlet::rounded(
            self.total_torque
                .map(|value| value * Dd::from_f64(duration)),
            "sphere_torque.angular_impulse",
        )?)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance(
        &mut self,
        index: usize,
        before: RigidBodyState,
        properties: MassProperties,
        gravity: Gravity,
        force: Vec3,
        duration: f64,
        gate: &CancelGate,
    ) -> Result<(RigidBodyState, f64), SphereContactError> {
        let Some(loads) = self.loads[index] else {
            return Ok((
                RigidBodyIntegrator::new(gravity)
                    .step(
                        before,
                        properties,
                        Wrench {
                            force_world: force,
                            torque_body: Vec3::ZERO,
                        },
                        duration,
                    )?
                    .state_after,
                0.0,
            ));
        };
        let count = (duration / self.controls.maximum_step_s).ceil().max(1.0);
        let remaining = self.controls.maximum_substeps - self.substeps;
        if !count.is_finite() || count > remaining as f64 {
            return Err(SphereContactError::TorqueBudget);
        }
        let count = count as usize;
        let step = duration / count as f64;
        if !step.is_finite() || step <= 0.0 {
            return Err(SphereContactError::Unresolved("torque substep duration"));
        }

        // Translation and rotation separate for a centre force and free couple.
        // Reuse the native exact centre drift ONCE, so rotational refinement
        // cannot move the analytic collision path through accumulated roundoff.
        let translation = RigidBodyIntegrator::new(gravity)
            .step(
                RigidBodyState::new(before.pose(), before.linear_momentum_world(), Vec3::ZERO)?,
                properties,
                Wrench {
                    force_world: force,
                    torque_body: Vec3::ZERO,
                },
                duration,
            )?
            .state_after;
        // The native spatial-load owner supplies full-inertia endpoint kicks,
        // free precession and couple work. A zero-COM translation factor avoids
        // obscuring its rotational work with unrelated translational energy.
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let mut rotation = RigidBodyState::new(
            Pose::new(Vec3::ZERO, before.pose().orientation())?,
            Vec3::ZERO,
            before.angular_momentum_body(),
        )?;
        let mut elapsed = 0.0;
        let mut work = Dd::ZERO;
        for index in 0..count {
            poll(gate)?;
            let end = if index + 1 == count {
                duration
            } else {
                step * (index + 1) as f64
            };
            let h = end - elapsed;
            if h <= 0.0 || !h.is_finite() {
                return Err(SphereContactError::Unresolved("torque substep clock"));
            }
            let result = integrator.step_spatial(rotation, properties, &loads, h)?;
            self.substeps += 1;
            rotation = result.step.state_after;
            work = work + Dd::from_f64(result.couple_work_j);
            elapsed = end;
        }
        let work = crate::finite_derived(work.to_f64(), "sphere_torque.couple_work")?;
        Ok((
            RigidBodyState::new(
                Pose::new(
                    translation.pose().position_world(),
                    rotation.pose().orientation(),
                )?,
                translation.linear_momentum_world(),
                rotation.angular_momentum_body(),
            )?,
            work,
        ))
    }
}
