//! Prescribed source temperature on the original finite-rotation rigid pair.
//!
//! Every point transfers once at its accepted geometry before Newton trials.
//! Reusable scratch carries those candidates; neither body nor accepted point
//! changes until the original native solve and complete energy/momentum audit
//! succeed. Mechanical iterations never repeat thermal damage or parameter work.

use super::*;
use fs_material::cohesive::source::ResolvedBilinearCohesive;

impl<'a> RigidCohesiveSystem<'a> {
    /// Advance the same two bodies with one temperature successor per original
    /// site, in original site order. Each source must retain its own card,
    /// selected claims, envelope, held coordinates and separation/rate bounds.
    /// All sites use the named original absolute-temperature axis for this call.
    ///
    /// Coefficients change first at fixed accepted geometry/damage; any further
    /// damage releases stored energy. Then the existing joint native body solve
    /// advances mechanics at these coefficients. Parameter work is a separate
    /// signed source input, and thermal damage is included once in total fracture
    /// loss. Subsequent ordinary `step` calls retain the accepted temperature.
    ///
    /// This first-order staggered model requires thermal/load step refinement;
    /// it does not supply heat transport, thermal expansion, entropy or a heat
    /// partition. Enabled damaged-area friction refuses because its separate
    /// source/history temperature transfer is not implemented.
    pub fn step_with_temperature(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        temperature_axis: &str,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_inner(loads, next_sources, temperature_axis, None)
    }

    /// Same native temperature/mechanics transaction with bounded cancellation.
    pub fn step_with_temperature_under_gate(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_inner(loads, next_sources, temperature_axis, Some(gate))
    }

    fn step_with_temperature_inner(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        temperature_axis: &str,
        gate: Option<&CancelGate>,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        poll(gate)?;
        if self.friction.is_some() {
            return Err(RigidCohesiveError::Invalid(
                "damaged-area friction requires its own source temperature/history transfer",
            ));
        }
        if next_sources.len() != self.sites.len() {
            return Err(RigidCohesiveError::Invalid(
                "one temperature successor per original cohesive site required",
            ));
        }
        // Inactive partial stages are scratch only and cannot affect ordinary
        // steps or retries. All entries are overwritten before activation.
        self.thermal_active = false;
        for (index, (site, source)) in self.sites.iter().zip(next_sources).enumerate() {
            poll(gate)?;
            if index >= self.config.maximum_material_evaluations {
                return Err(RigidCohesiveError::Budget(
                    "original cohesive temperature transfers",
                ));
            }
            self.staged_temperature[index] = Some(site.point.transfer_temperature(
                *source,
                temperature_axis,
                self.config.objective,
            )?);
        }
        poll(gate)?;
        self.thermal_active = true;
        // Drop the returned frame borrow before restoring the scratch selector.
        // Every Result path, including late physical/cancellation refusals, exits
        // with ordinary evaluation bound to the original accepted site owners.
        let result = self.step_inner(loads, gate).map(|_| ());
        self.thermal_active = false;
        result?;
        Ok(self
            .last_frame
            .as_ref()
            .expect("native step published frame"))
    }

    pub(super) fn point_for_step(
        &self,
        index: usize,
    ) -> Result<SourceBoundObjectiveCohesivePoint<'a>, RigidCohesiveError> {
        if self.thermal_active {
            self.staged_temperature[index]
                .as_ref()
                .map(SourceBoundObjectiveCohesiveTemperatureStep::point)
                .ok_or(RigidCohesiveError::Invalid("incomplete temperature stage"))
        } else {
            Ok(self.sites[index].point)
        }
    }

    pub(super) fn temperature_accounts(
        &self,
        gate: Option<&CancelGate>,
    ) -> Result<[f64; 3], RigidCohesiveError> {
        let mut accounts = [0.0; 3];
        if self.thermal_active {
            for staged in &self.staged_temperature {
                poll(gate)?;
                let step = staged
                    .as_ref()
                    .ok_or(RigidCohesiveError::Invalid("incomplete temperature stage"))?;
                accounts[0] += step.parameter_work_j;
                accounts[1] += step.stored_change_j;
                accounts[2] += step.damage_dissipated_j;
            }
        }
        for value in accounts {
            finite(value)?;
        }
        Ok(accounts)
    }
}
