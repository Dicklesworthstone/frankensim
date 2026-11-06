//! Complete physical intervals through native released-face impact events.
//!
//! Motion uses the original work-conjugate midpoint and immutable material
//! history; contact uses the original simultaneous consistent-mass impulse.
//! Each actual prefix samples the original affine loads at its own times.
//! Only a complete duration publishes a successor. No penetrating candidate,
//! artificial gap correction, reset history or unresolved zero-time event is
//! allowed to escape into the retained trajectory.

use fs_contact::multicontact::MultiImpactResult;
use fs_exec::Cx;

use super::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError, BrittleFaceContactLocation,
};
use super::{
    BrittleFaceImpactControls, BrittleFaceImpactError, BrittleFaceImpactRow, SolidImpulseReport,
};
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointEnergy;
use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::path::{Sum, interpolate};
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;

/// Bounds on the complete duration, with the original motion/contact controls.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceTimeControls {
    /// Largest free-flight interval proposed to the native prefix locator.
    /// This is a caller-selected temporal resolution, not an error estimate.
    pub maximum_interval_s: f64,
    pub maximum_intervals: usize,
    /// Includes instantaneous events at the original and final physical times.
    pub maximum_impacts: usize,
    /// All native motion solves, including discarded localization trials.
    pub maximum_trials: usize,
    pub approach: BrittleFaceApproachControls,
    pub impact: BrittleFaceImpactControls,
    /// Whole-duration signed balance AND sum of absolute local balance errors.
    pub total_energy_tolerance_j: f64,
}

/// One accepted physical midpoint interval, without a copy of its full fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceMotionRecord {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub energy: BrittleTetMidpointEnergy,
}

/// Original coupled-contact and body-impulse receipts at one fixed configuration.
#[derive(Debug)]
pub struct BrittleFaceTimedImpact {
    pub location: BrittleFaceContactLocation,
    pub contact: MultiImpactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub maximum_velocity_residual_mps: f64,
    pub relative_work_error_j: f64,
    pub energy_balance_error_j: f64,
    pub absolute_contact_couple_kg_m2_s: f64,
}

/// Whole-duration physical work. Impact loss remains separate from fracture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceTimeEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Compensated sum of native polarized kinetic increments, including events.
    pub kinetic_change_j: f64,
    pub material_work_j: f64,
    pub fracture_dissipation_j: f64,
    /// Mechanical loss from the contact owner; no thermal allocation is inferred.
    pub contact_dissipation_j: f64,
    pub dead_load_mean_work_j: f64,
    pub support_mean_work_j: f64,
    /// Signed reaction work at held, possibly moving velocity components.
    pub support_impulse_work_j: f64,
    pub material_balance_error_j: f64,
    pub mean_work_error_j: f64,
    pub absolute_mean_work_error_j: f64,
    /// External work minus elastic/kinetic changes and fracture/contact losses.
    pub total_balance_error_j: f64,
    /// Prevents opposite local numerical defects from hiding in the total.
    pub absolute_local_balance_error_j: f64,
}

impl BrittleFaceTimeEnergy {
    #[must_use]
    pub fn external_work_j(&self) -> f64 {
        (self.dead_load_mean_work_j + self.support_mean_work_j) + self.support_impulse_work_j
    }
}

#[derive(Debug)]
pub struct BrittleFaceTimeStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub intervals: Vec<BrittleFaceMotionRecord>,
    pub impacts: Vec<BrittleFaceTimedImpact>,
    pub trials: usize,
    /// Completed native trial telemetry; a refused native solve has no receipt.
    pub numerical_tangent_evaluations: usize,
    pub energy: BrittleFaceTimeEnergy,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleFaceTimeError {
    Invalid(&'static str),
    Allocation,
    Motion(BrittleTetDynamicError),
    Approach(BrittleFaceApproachError),
    Impact(BrittleFaceImpactError),
    WorkLimit {
        intervals: usize,
        impacts: usize,
        trials: usize,
    },
    /// No two events may occur within one localization time tolerance.
    /// An instantaneous impulse cannot replace a sustained contact-force solve.
    UnresolvedContact {
        time_s: f64,
    },
}

impl core::fmt::Display for BrittleFaceTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "released-face physical interval refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceTimeError {}
impl From<BrittleTetError> for BrittleFaceTimeError {
    fn from(error: BrittleTetError) -> Self {
        Self::Motion(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Advance the complete affine load interval, automatically resolving impacts.
    ///
    /// Supports follow an affine displacement path; their initial velocities
    /// must match its slope within the explicit velocity tolerance. The same
    /// support DOFs retain their velocities during each instantaneous impact.
    /// Applied-force pairs give values at this duration's original endpoints.
    /// Every shortened prefix and remainder recomputes its own mean forces.
    ///
    /// Each accepted midpoint consumes one native motion step. Impacts consume
    /// an event allowance and no physical timestep. Discarded localization
    /// trials consume only the solve budget. An exhausted budget, cancellation,
    /// unsupported geometry, physical-work refusal or unresolved contact returns
    /// no partially advanced motion; the original is available for retry.
    ///
    /// The native prefix locator admits nonpenetrating affine displacement
    /// traces on fixed paired faces. Its bracket is a local numerical event
    /// resolution, not a global continuous-collision certificate. Choose the
    /// maximum interval to resolve forcing and motion; split load reversals at
    /// their knots. Finite sliding, new fragment pairs, friction, persistent
    /// contact forces and thermal allocation of impact loss are outside this
    /// entry. Repeated events separated by no more than the localization time
    /// tolerance refuse instead of generating an unresolved impulse sequence.
    /// No temporal convergence order is asserted across impact events.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_time(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceTimeStep<'a>, BrittleFaceTimeError> {
        use BrittleFaceTimeError as E;
        checkpoint(cx)?;
        let target_time = self.time_s + duration_s;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !target_time.is_finite()
            || target_time <= self.time_s
            || !controls.maximum_interval_s.is_finite()
            || controls.maximum_interval_s <= 0.0
            || controls.maximum_intervals == 0
            || controls.maximum_trials == 0
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j < 0.0
            || controls.approach.activation_tolerance_m != controls.impact.activation_tolerance_m
            || controls.approach.velocity_tolerance_mps < controls.impact.velocity_tolerance_mps
            || restitution.len() != self.solid.released_faces.len()
            || restitution
                .iter()
                .any(|e| !e.is_finite() || !(0.0..=1.0).contains(e))
            || origin_m.iter().any(|x| !x.is_finite())
        {
            return Err(E::Invalid(
                "positive bounded duration, compatible event gates and complete restitution required",
            ));
        }
        let duration = target_time - self.time_s;
        let initial_u = self.solid.solution.displacement_m.as_flattened();
        let initial_v = self.velocity_mps.as_flattened();
        if prescribed_end_m.len() > initial_u.len() || force_intervals_n.len() > initial_u.len() {
            return Err(E::Invalid("load entries exceed the original DOF count"));
        }
        let mut supports = buffer(prescribed_end_m.len())?;
        let mut held = buffer(prescribed_end_m.len())?;
        for &(dof, end) in prescribed_end_m {
            checkpoint(cx)?;
            let Some(&start) = initial_u.get(dof) else {
                return Err(E::Invalid("prescribed DOF outside the retained mesh"));
            };
            let speed = (end - start) / duration;
            if !end.is_finite()
                || !speed.is_finite()
                || (speed - initial_v[dof]).abs() > controls.approach.velocity_tolerance_mps
            {
                return Err(E::Invalid(
                    "initial support velocity must match its affine path",
                ));
            }
            supports.push((dof, start));
            held.push(dof);
        }
        let mut forces = buffer(force_intervals_n.len())?;
        for &(dof, values) in force_intervals_n {
            checkpoint(cx)?;
            if dof >= initial_u.len() || values.iter().any(|x| !x.is_finite()) {
                return Err(E::Invalid(
                    "finite applied forces on existing DOFs required",
                ));
            }
            forces.push((dof, values));
        }
        let mut current: Option<BrittleTetMotion<'a>> = None;
        let mut intervals = buffer(controls.maximum_intervals.min(32))?;
        let mut impacts: Vec<BrittleFaceTimedImpact> = buffer(controls.maximum_impacts.min(8))?;
        let mut sums = [Sum::default(); 10];
        let mut trials = 0usize;
        let mut tangent_evaluations = 0usize;
        loop {
            checkpoint(cx)?;
            let accepted = current.as_ref().unwrap_or(self);
            if accepted.time_s == target_time {
                break;
            }
            if trials >= controls.maximum_trials || intervals.len() == controls.maximum_intervals {
                return Err(E::WorkLimit {
                    intervals: intervals.len(),
                    impacts: impacts.len(),
                    trials,
                });
            }
            let end = (accepted.time_s + controls.maximum_interval_s).min(target_time);
            if !end.is_finite() || end <= accepted.time_s {
                return Err(E::Invalid(
                    "increasing native physical time must be representable",
                ));
            }
            let start_fraction = (accepted.time_s - self.time_s) / duration;
            let end_fraction = if end == target_time {
                1.0
            } else {
                (end - self.time_s) / duration
            };
            for (i, (entry, &(dof, target))) in
                supports.iter_mut().zip(prescribed_end_m).enumerate()
            {
                if i % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = interpolate(initial_u[dof], target, end_fraction);
            }
            for (i, (entry, &(_, [start, target]))) in
                forces.iter_mut().zip(force_intervals_n).enumerate()
            {
                if i % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = [
                    interpolate(start, target, start_fraction),
                    interpolate(start, target, end_fraction),
                ];
            }
            let mut approach_controls = controls.approach;
            approach_controls.maximum_trials = approach_controls
                .maximum_trials
                .min(controls.maximum_trials - trials);
            let prefix = accepted
                .advance_toward_released_face_contact(
                    end - accepted.time_s,
                    &supports,
                    &forces,
                    settings,
                    approach_controls,
                    cx,
                )
                .map_err(E::Approach)?;
            trials = trials
                .checked_add(prefix.trials)
                .ok_or(E::Invalid("native trial count overflow"))?;
            tangent_evaluations = tangent_evaluations
                .checked_add(prefix.numerical_tangent_evaluations)
                .ok_or(E::Invalid("native tangent count overflow"))?;
            if prefix.advance.is_none() && prefix.contact.is_none() {
                return Err(E::Invalid(
                    "prefix must advance physical time or resolve a contact",
                ));
            }
            if let Some(step) = prefix.advance {
                if intervals.len() == controls.maximum_intervals {
                    return Err(E::WorkLimit {
                        intervals: intervals.len(),
                        impacts: impacts.len(),
                        trials,
                    });
                }
                if step.motion.time_s <= accepted.time_s || step.motion.time_s > end {
                    return Err(E::Invalid(
                        "native prefix escaped its requested physical interval",
                    ));
                }
                let energy = step.energy;
                accumulate(
                    &mut sums,
                    [
                        energy.material_work_j,
                        energy.fracture_dissipation_j,
                        energy.dead_load_mean_work_j,
                        energy.support_mean_work_j,
                        0.0,
                        energy.inertia.kinetic_change_j,
                        0.0,
                        energy.total_balance_error_j.abs(),
                        energy.mean_work_error_j.abs(),
                        energy.mean_internal_work_j,
                    ],
                )?;
                push(
                    &mut intervals,
                    BrittleFaceMotionRecord {
                        start_time_s: accepted.time_s,
                        end_time_s: step.motion.time_s,
                        energy,
                    },
                )?;
                current = Some(step.motion);
            }
            if let Some(location) = prefix.contact {
                let at_contact = current.as_ref().unwrap_or(self);
                if location.time_s != at_contact.time_s
                    || impacts.last().is_some_and(|last| {
                        at_contact.time_s - last.location.time_s
                            <= controls.approach.time_tolerance_s
                    })
                {
                    return Err(E::UnresolvedContact {
                        time_s: at_contact.time_s,
                    });
                }
                if impacts.len() == controls.maximum_impacts {
                    return Err(E::WorkLimit {
                        intervals: intervals.len(),
                        impacts: impacts.len(),
                        trials,
                    });
                }
                let hit = at_contact
                    .impact_released_faces(restitution, &held, origin_m, controls.impact, cx)
                    .map_err(E::Impact)?;
                if !hit.contact.normal_impulses.iter().any(|&p| p > 0.0) {
                    return Err(E::UnresolvedContact {
                        time_s: at_contact.time_s,
                    });
                }
                accumulate(
                    &mut sums,
                    [
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        hit.impulse.support_work_j,
                        hit.impulse.kinetic_change_j,
                        hit.contact.dissipated_energy,
                        hit.energy_balance_error_j.abs(),
                        0.0,
                        0.0,
                    ],
                )?;
                push(
                    &mut impacts,
                    BrittleFaceTimedImpact {
                        location,
                        contact: hit.contact,
                        impulse: hit.impulse,
                        rows: hit.rows,
                        maximum_velocity_residual_mps: hit.maximum_velocity_residual_mps,
                        relative_work_error_j: hit.relative_work_error_j,
                        energy_balance_error_j: hit.energy_balance_error_j,
                        absolute_contact_couple_kg_m2_s: hit.absolute_contact_couple_kg_m2_s,
                    },
                )?;
                current = Some(hit.motion);
            }
        }
        let motion = current.ok_or(E::Invalid("no physical interval accepted"))?;
        let [
            material,
            fracture,
            dead,
            support,
            impulse_support,
            kinetic,
            contact,
            absolute_balance,
            absolute_mean,
            internal,
        ] = sums.map(Sum::value);
        let storage_change = motion.solid.stored_energy_j - self.solid.stored_energy_j;
        let energy = BrittleFaceTimeEnergy {
            stored_before_j: self.solid.stored_energy_j,
            stored_after_j: motion.solid.stored_energy_j,
            kinetic_before_j: self.kinetic_energy_j,
            kinetic_after_j: motion.kinetic_energy_j,
            kinetic_change_j: kinetic,
            material_work_j: material,
            fracture_dissipation_j: fracture,
            contact_dissipation_j: contact,
            dead_load_mean_work_j: dead,
            support_mean_work_j: support,
            support_impulse_work_j: impulse_support,
            material_balance_error_j: (material - storage_change) - fracture,
            mean_work_error_j: internal - material,
            absolute_mean_work_error_j: absolute_mean,
            total_balance_error_j: (((((dead + support) + impulse_support) - storage_change)
                - fracture)
                - contact)
                - kinetic,
            absolute_local_balance_error_j: absolute_balance,
        };
        for (account, defect, tolerance) in [
            (
                "released-face time material work/storage/fracture",
                energy.material_balance_error_j,
                controls.approach.work.material_energy_tolerance_j,
            ),
            (
                "released-face time cumulative fracture loss",
                (motion.solid.fracture_dissipation_j - self.solid.fracture_dissipation_j)
                    - fracture,
                controls.approach.work.material_energy_tolerance_j,
            ),
            (
                "released-face time mean-stress work",
                energy.mean_work_error_j,
                controls.approach.work.material_energy_tolerance_j,
            ),
            (
                "released-face time absolute mean-stress work",
                absolute_mean,
                controls.approach.work.material_energy_tolerance_j,
            ),
            (
                "released-face time complete work balance",
                energy.total_balance_error_j,
                controls.total_energy_tolerance_j,
            ),
            (
                "released-face time absolute local work balance",
                absolute_balance,
                controls.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        checkpoint(cx)?;
        Ok(BrittleFaceTimeStep {
            motion,
            intervals,
            impacts,
            trials,
            numerical_tangent_evaluations: tangent_evaluations,
            energy,
        })
    }
}

fn buffer<T>(capacity: usize) -> Result<Vec<T>, BrittleFaceTimeError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| BrittleFaceTimeError::Allocation)?;
    Ok(result)
}
fn push<T>(values: &mut Vec<T>, value: T) -> Result<(), BrittleFaceTimeError> {
    if values.len() == values.capacity() {
        values
            .try_reserve(1)
            .map_err(|_| BrittleFaceTimeError::Allocation)?;
    }
    values.push(value);
    Ok(())
}
fn accumulate(sums: &mut [Sum; 10], values: [f64; 10]) -> Result<(), BrittleFaceTimeError> {
    for (sum, value) in sums.iter_mut().zip(values) {
        sum.add(value);
        if !value.is_finite() || !sum.value().is_finite() {
            return Err(BrittleFaceTimeError::Invalid(
                "complete physical work must be representable",
            ));
        }
    }
    Ok(())
}
