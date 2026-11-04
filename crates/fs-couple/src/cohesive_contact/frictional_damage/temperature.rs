//! Fixed-jump temperature continuation of the original progressively damaged bond and friction state.
use super::*;
use fs_material::cohesive::temperature::CohesiveTemperatureStep;

impl<'a> DamageFrictionLoadLine<'a> {
    /// Stage temperature successors for this SAME original bond and friction
    /// source without advancing mechanics. Both successors retain their own
    /// card, selected claims, held coordinates and original surface ordering,
    /// and must have exactly the same query as each other.
    ///
    /// The original material transfer runs once at the accepted jump. Its
    /// signed coefficient work and nonnegative damage release remain separate;
    /// changing Coulomb coefficients creates no friction work or heat. Original
    /// branch occupancy is retained exactly, including an occupied kinetic
    /// branch. No endpoint velocity, elapsed time or new bond is inferred.
    ///
    /// Area, full compliance and every original control/work allowance are
    /// retained. The new coefficients pass the original owner's admission
    /// checks before a candidate owner and its rebound state are returned.
    /// The caller can publish both atomically after its own coupled audit;
    /// this owner and the supplied state are unchanged on every refusal.
    pub fn transfer_temperature<'b>(
        &self,
        previous: DamageFrictionState<'a>,
        next_cohesive: &'b ResolvedBilinearCohesive,
        next_friction: &'b ResidualCoulombSource,
        temperature_axis: &str,
    ) -> Result<
        (
            DamageFrictionLoadLine<'b>,
            DamageFrictionState<'b>,
            CohesiveTemperatureStep<'b>,
        ),
        DamageFrictionError,
    > {
        self.require_state(previous)?;
        self.cohesive
            .require_temperature_successor(next_cohesive, temperature_axis)?;
        self.friction
            .require_temperature_successor(next_friction, temperature_axis)?;
        // Construction checks the paired query as well as physical admission
        // at the new coefficients. It cannot reset the accepted bond history.
        let candidate = DamageFrictionLoadLine::new(
            next_cohesive,
            next_friction,
            self.area,
            self.compliance,
            self.config,
        )?;
        let transfer = previous.point.transfer_temperature(
            next_cohesive.law(),
            self.config.contact.material_energy_tolerance_j_m2,
        )?;
        let parameter = finite(self.area * transfer.parameter_work_j_m2)?;
        let storage = finite(self.area * transfer.stored_change_j_m2)?;
        let damage = finite(self.area * transfer.damage_dissipated_j_m2)?;
        let residual = finite(parameter - storage - damage)?;
        let tolerance = finite(
            self.config.work_absolute_tolerance_j
                + self.config.work_relative_tolerance
                    * parameter.abs().max(storage.abs()).max(damage),
        )?;
        if residual.abs() > tolerance {
            return Err(DamageFrictionError::WorkMismatch {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        let state = DamageFrictionState {
            identity: candidate.identity,
            point: transfer.point,
            occupancy: previous.occupancy,
        };
        Ok((candidate, state, transfer))
    }
}
