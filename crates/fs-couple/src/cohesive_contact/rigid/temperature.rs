//! Prescribed source temperature on the original finite-rotation rigid pair.
//!
//! Every point transfers once at its accepted geometry before Newton trials.
//! Reusable scratch carries those candidates; neither body nor accepted point
//! changes until the original native solve and complete energy/momentum audit
//! succeed. Mechanical iterations never repeat thermal damage or parameter work.

use super::*;
use crate::cohesive_contact::reclosure::ResidualCoulombSource;
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
    /// partition. Enabled damaged-area friction requires the paired-source
    /// `step_with_temperature_and_friction` operation.
    pub fn step_with_temperature(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        temperature_axis: &str,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_inner(loads, next_sources, None, temperature_axis, None)
    }

    /// Same native temperature/mechanics transaction with bounded cancellation.
    pub fn step_with_temperature_under_gate(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_inner(loads, next_sources, None, temperature_axis, Some(gate))
    }

    /// Advance the original damaged-area friction and cohesive bond together at
    /// prescribed source temperatures. Both source arrays follow original site
    /// order. Each friction successor must retain its original card/claims and
    /// have exactly the same query as its paired cohesive successor.
    ///
    /// The existing rigid Coulomb specialization has no elastic friction spring
    /// or stored slip energy. Updating mu at fixed accepted jump therefore does
    /// zero friction work and creates zero heat. Original cohesive coefficient
    /// work and fixed-jump damage release are counted once; only the subsequent
    /// actual mechanical slip contributes friction loss. Static/kinetic branch
    /// history survives, including physical motion at an accepted endpoint.
    pub fn step_with_temperature_and_friction(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        next_friction_sources: &[&'a ResidualCoulombSource],
        temperature_axis: &str,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_inner(
            loads,
            next_sources,
            Some(next_friction_sources),
            temperature_axis,
            None,
        )
    }

    /// Same paired-source transaction with bounded cooperative cancellation.
    pub fn step_with_temperature_and_friction_under_gate(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        next_friction_sources: &[&'a ResidualCoulombSource],
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_inner(
            loads,
            next_sources,
            Some(next_friction_sources),
            temperature_axis,
            Some(gate),
        )
    }

    fn step_with_temperature_inner(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        next_friction_sources: Option<&[&'a ResidualCoulombSource]>,
        temperature_axis: &str,
        gate: Option<&CancelGate>,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_with_temperature_acceptance(
            loads,
            next_sources,
            next_friction_sources,
            temperature_axis,
            gate,
            |_, _| Ok::<(), RigidCohesiveError>(()),
        )?;
        Ok(self
            .last_frame
            .as_ref()
            .expect("native step published frame"))
    }

    /// Stage the same prescribed source transfer and native solve, then let a
    /// coupled owner audit its candidate before either accepted state changes.
    pub(crate) fn step_with_temperature_acceptance<P, E: From<RigidCohesiveError>>(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        next_sources: &[&'a ResolvedBilinearCohesive],
        next_friction_sources: Option<&[&'a ResidualCoulombSource]>,
        temperature_axis: &str,
        gate: Option<&CancelGate>,
        accept: impl FnOnce(&RigidCohesiveFrame, &[Option<RigidCohesiveSiteFrame>]) -> Result<P, E>,
    ) -> Result<P, E> {
        poll(gate)?;
        match (self.friction.is_some(), next_friction_sources) {
            (true, None) => {
                return Err(RigidCohesiveError::Invalid(
                    "damaged-area friction requires paired temperature successors",
                )
                .into());
            }
            (false, Some(_)) => {
                return Err(RigidCohesiveError::Invalid(
                    "temperature friction sources require the original enabled friction owner",
                )
                .into());
            }
            (_, Some(sources)) if sources.len() != self.sites.len() => {
                return Err(RigidCohesiveError::Invalid(
                    "one temperature successor per original friction site required",
                )
                .into());
            }
            _ => {}
        }
        if next_sources.len() != self.sites.len() {
            return Err(RigidCohesiveError::Invalid(
                "one temperature successor per original cohesive site required",
            )
            .into());
        }
        // Inactive partial stages are scratch only and cannot affect ordinary
        // steps or retries. All entries are overwritten before activation.
        self.thermal_active = false;
        for (index, (site, source)) in self.sites.iter().zip(next_sources).enumerate() {
            poll(gate)?;
            if index >= self.config.maximum_material_evaluations {
                return Err(
                    RigidCohesiveError::Budget("original cohesive temperature transfers").into(),
                );
            }
            self.staged_temperature[index] = Some(
                site.point
                    .transfer_temperature(*source, temperature_axis, self.config.objective)
                    .map_err(RigidCohesiveError::from)?,
            );
        }
        if let Some(next) = next_friction_sources {
            // Move only the reusable owner while validating its staged source
            // and branch history against the SAME canonical thermal points.
            // Restore it before propagating every refusal or cancellation.
            let mut friction = self.friction.take().expect("enabled friction checked");
            let staged = friction.stage_temperature(self, next, temperature_axis, gate);
            self.friction = Some(friction);
            staged?;
        }
        poll(gate)?;
        self.thermal_active = true;
        // Every Result path, including late physical/cancellation refusals, exits
        // with ordinary evaluation bound to the original accepted site owners.
        let result = self.step_with_acceptance(loads, gate, accept);
        self.thermal_active = false;
        result
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
