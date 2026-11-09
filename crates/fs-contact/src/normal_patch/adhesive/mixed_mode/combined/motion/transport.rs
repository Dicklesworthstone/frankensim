//! Objective full-stick history under a changed native contact projection.
//!
//! Only the two elastic shear components rotate. Normal approach and relative
//! twist retain their physical values on the same original common-radius
//! source. A common material spin or tangent gauge change is not a shear path,
//! a torsional increment, a material refresh or a new bond.
//!
//! The projected kinetic partition v^T W^-1 v/2 may change when the actual
//! native projection changes. Its difference is diagnostic, never supplied
//! work or heat. The native owner must separately admit the pose path, body
//! momenta, complete kinetic energy and integrated geometric travel.

use super::{
    JkrCombinedMotion, JkrCombinedMotionError, JkrCombinedMotionState, JkrCombinedPathError,
    JkrCombinedResponse, bounded, checked_array, finite, poll,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::tangential::transport::TangentialHistoryRotation;

/// Same source storage and chronology under an authenticated new projection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedProjectionReport {
    pub before: JkrCombinedResponse,
    pub after: JkrCombinedResponse,
    pub rotation: TangentialHistoryRotation,
    /// Numerical storage change only; this is not mechanical work or heat.
    pub source_energy_residual_j: f64,
    pub source_energy_tolerance_j: f64,
    /// [Normal elastic, vector-shear elastic, torsional elastic, surface] [J].
    pub storage_component_residual_j: [f64; 4],
    /// Sum of absolute component defects; opposite signs cannot hide a defect.
    pub absolute_storage_residual_j: f64,
    pub projected_kinetic_energy_before_j: f64,
    pub projected_kinetic_energy_after_j: f64,
    /// Diagnostic repartition, not actual body work, loss or heat [J].
    pub projected_kinetic_repartition_j: f64,
    /// Cached old endpoint; at most one new original source evaluation.
    pub evaluations: usize,
}

/// Original physical history and its new native projection, adopted together.
#[derive(Clone, Debug)]
pub struct JkrCombinedProjectionContinuation<'a> {
    previous: JkrCombinedMotionState,
    motion: JkrCombinedMotion<'a>,
    state: JkrCombinedMotionState,
    report: JkrCombinedProjectionReport,
}
impl<'a> JkrCombinedProjectionContinuation<'a> {
    pub const fn previous_state(&self) -> JkrCombinedMotionState {
        self.previous
    }
    pub const fn motion(&self) -> &JkrCombinedMotion<'a> {
        &self.motion
    }
    pub const fn state(&self) -> JkrCombinedMotionState {
        self.state
    }
    pub const fn report(&self) -> &JkrCombinedProjectionReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrCombinedMotion<'a>,
        JkrCombinedMotionState,
        JkrCombinedProjectionReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

impl<'a> JkrCombinedMotion<'a> {
    /// Retain attached I/II/III history under a new actual native projection.
    ///
    /// The checked tribology isometry transports only elastic shear. Normal
    /// approach, relative twist, original source, every budget, clock and step
    /// count stay fixed. A pure gauge change has W'=P W P^T and v'=P v with
    /// P=diag(1,R,1); a physical change of projection can repartition observable
    /// kinetic energy while leaving actual body energy unchanged.
    ///
    /// Native geometry owns the supplied mobility, rates and physical rotation
    /// path. This source operation cannot authenticate body kinematics. It does
    /// not integrate a chord between rotated shear components or apply the
    /// physical-increment limit to a coordinate re-expression. Original rate,
    /// effort, strict-attachment and source-applicability gates still apply.
    /// Both complete storage and the absolute defects of its four components
    /// must close at the original energy tolerance. No cancellation between
    /// storage terms can authorize physical work here.
    ///
    /// The old response is sealed in accepted history. Unchanged coordinates
    /// need no source query; changed shear coordinates need one original point
    /// admission. Failure or cancellation changes no accepted owner or state.
    pub fn rebind_projection(
        &self,
        previous: JkrCombinedMotionState,
        rotation: TangentialHistoryRotation,
        new_mobility: [[f64; 4]; 4],
        new_rates: [f64; 4],
        gate: &CancelGate,
    ) -> Result<JkrCombinedProjectionContinuation<'a>, JkrCombinedMotionError> {
        poll(gate)?;
        self.check(previous)?;
        let before = previous.contact;
        let shear = rotation
            .apply([before.coordinates[1], before.coordinates[2]])
            .map_err(|error| {
                JkrCombinedPathError::Source(super::super::JkrCombinedError::Source(error.into()))
            })?;
        let coordinates = [
            before.coordinates[0],
            shear[0],
            shear[1],
            before.coordinates[3],
        ];
        let motion = Self::new(self.source, new_mobility, self.budget)?;
        bounded(new_rates, self.budget.maximum_rate, "reprojected rates")?;
        let evaluations = usize::from(coordinates != before.coordinates);
        if evaluations > self.budget.maximum_evaluations {
            return Err(JkrCombinedMotionError::Budget);
        }
        let after = if evaluations == 0 {
            before
        } else {
            motion.initial_state(coordinates, new_rates, gate)?.contact
        };
        bounded(
            after.force,
            self.budget.maximum_effort,
            "reprojected source efforts",
        )?;
        let old_storage = storage(before);
        let new_storage = storage(after);
        let storage_component_residual_j = checked_array(
            core::array::from_fn(|i| {
                (Dd::from_f64(new_storage[i]) - Dd::from_f64(old_storage[i])).to_f64()
            }),
            "reprojected source storage components",
        )?;
        let absolute_storage_residual_j = finite(
            storage_component_residual_j
                .into_iter()
                .fold(Dd::from_f64(0.0), |sum, x| sum + Dd::from_f64(x.abs()))
                .to_f64(),
            "absolute reprojected source storage defect",
        )?;
        let source_energy_residual_j = finite(
            (Dd::from_f64(after.reversible_energy_j) - Dd::from_f64(before.reversible_energy_j))
                .to_f64(),
            "reprojected source storage difference",
        )?;
        let scale = finite(
            old_storage
                .into_iter()
                .chain(new_storage)
                .fold(Dd::from_f64(0.0), |sum, x| sum + Dd::from_f64(x.abs()))
                .to_f64(),
            "reprojected source storage scale",
        )?;
        let source_energy_tolerance_j = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "reprojected source storage tolerance",
        )?;
        if source_energy_residual_j
            .abs()
            .max(absolute_storage_residual_j)
            > source_energy_tolerance_j
        {
            return Err(JkrCombinedMotionError::EnergyBalance {
                residual_j: source_energy_residual_j
                    .abs()
                    .max(absolute_storage_residual_j),
                tolerance_j: source_energy_tolerance_j,
            });
        }
        let projected_kinetic_energy_before_j = self.kinetic_energy(previous.rates)?;
        let projected_kinetic_energy_after_j = motion.kinetic_energy(new_rates)?;
        let projected_kinetic_repartition_j = finite(
            (Dd::from_f64(projected_kinetic_energy_after_j)
                - Dd::from_f64(projected_kinetic_energy_before_j))
            .to_f64(),
            "projected kinetic repartition",
        )?;
        let state =
            if motion.identity == self.identity && after == before && new_rates == previous.rates {
                previous
            } else {
                JkrCombinedMotionState {
                    identity: motion.identity,
                    contact: after,
                    rates: new_rates,
                    ..previous
                }
            };
        poll(gate)?;
        Ok(JkrCombinedProjectionContinuation {
            previous,
            motion,
            state,
            report: JkrCombinedProjectionReport {
                before,
                after,
                rotation,
                source_energy_residual_j,
                source_energy_tolerance_j,
                storage_component_residual_j,
                absolute_storage_residual_j,
                projected_kinetic_energy_before_j,
                projected_kinetic_energy_after_j,
                projected_kinetic_repartition_j,
                evaluations,
            },
        })
    }
}

fn storage(response: JkrCombinedResponse) -> [f64; 4] {
    [
        response.normal_elastic_energy_j,
        response.tangential_elastic_energy_j,
        response.torsional_elastic_energy_j,
        response.surface_energy_j,
    ]
}
