//! Source continuation of the original joint damaged-area patch.
//! Only local source bindings and their numerical cone scales are replaced;
//! complete cross-compliance, prepared dense algebra and solver scratch survive.
use super::*;
use fs_material::cohesive::temperature::CohesiveTemperatureStep;

/// Private candidate for the same solver that prepared it. Its construction
/// performs every fallible source, material, work and scale check.
pub(crate) struct PreparedPatchTemperature<'a> {
    laws: Vec<LocalDamageFriction<'a>>,
    states: Vec<DamageFrictionState<'a>>,
    cone_scale: Vec<f64>,
    transfers: Vec<CohesiveTemperatureStep<'a>>,
}
impl<'a> PreparedPatchTemperature<'a> {
    pub(crate) fn transfers(&self) -> &[CohesiveTemperatureStep<'a>] {
        &self.transfers
    }
}

impl<'a> DamageFrictionPatchLoadLine<'a> {
    /// Stage each original bond once, preserving its friction occupancy and
    /// local identity. The outer owner supplies its material tolerance so both
    /// its original gate and this friction owner's original gate apply.
    pub(crate) fn prepare_temperature_successors(
        &self,
        previous: &[DamageFrictionState<'a>],
        next_cohesive: &[&'a ResolvedBilinearCohesive],
        next_friction: &[&'a ResidualCoulombSource],
        temperature_axis: &str,
        material_tolerance_j_m2: f64,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<PreparedPatchTemperature<'a>, DamageFrictionPatchError> {
        poll(&mut keep_going)?;
        let n = self.laws.len();
        if previous.len() != n || next_cohesive.len() != n || next_friction.len() != n {
            return Err(DamageFrictionPatchError::Invalid(
                "temperature continuation requires every original patch site",
            ));
        }
        positive(material_tolerance_j_m2)?;
        let tolerance = material_tolerance_j_m2
            .min(self.config.friction.contact.material_energy_tolerance_j_m2);
        // Check the complete original binding before any material transfer.
        for p in 0..n {
            poll(&mut keep_going)?;
            self.laws[p].require_state(previous[p]).map_err(|e| local_error(p, e))?;
            self.laws[p].cohesive
                .require_temperature_successor(next_cohesive[p], temperature_axis)
                .map_err(|e| local_error(p, e.into()))?;
            self.laws[p].friction
                .require_temperature_successor(next_friction[p], temperature_axis)
                .map_err(|e| local_error(p, e.into()))?;
            next_friction[p].require_cohesive_source(next_cohesive[p])
                .map_err(|e| local_error(p, e.into()))?;
        }
        let mut candidate = PreparedPatchTemperature {
            laws: reserved(n)?,
            states: reserved(n)?,
            cone_scale: reserved(n)?,
            transfers: reserved(n)?,
        };
        let dimension = 3 * n;
        for p in 0..n {
            poll(&mut keep_going)?;
            let transfer = previous[p].point().transfer_temperature(
                next_cohesive[p].law(), tolerance,
            ).map_err(|e| local_error(p, e.into()))?;
            let (law, state) = self.laws[p].adopt_temperature_successor(
                previous[p], next_cohesive[p], next_friction[p],
                transfer.point, temperature_axis,
            ).map_err(|e| local_error(p, e))?;
            let parameter = finite(law.area * transfer.parameter_work_j_m2)?;
            let storage = finite(law.area * transfer.stored_change_j_m2)?;
            let damage = finite(law.area * transfer.damage_dissipated_j_m2)?;
            let residual = finite(parameter - storage - damage)?;
            let work_tolerance = finite(
                law.config.work_absolute_tolerance_j
                    + law.config.work_relative_tolerance
                        * parameter.abs().max(storage.abs()).max(damage),
            )?;
            if residual.abs() > work_tolerance {
                return Err(local_error(p, DamageFrictionError::WorkMismatch {
                    residual_j: residual,
                    tolerance_j: work_tolerance,
                }));
            }
            // Reuse the existing complete response; no diagonal surrogate,
            // rebuilt solver, copied compliance or new friction spring.
            let response = (1..3)
                .map(|axis| {
                    self.compliance[(3 * p + axis) * dimension..(3 * p + axis + 1) * dimension]
                        .iter().map(|x| x.abs()).sum::<f64>()
                })
                .fold(0.0, f64::max);
            let material_length = next_cohesive[p].law().onset_separation_m();
            let material_force = law.area * next_cohesive[p].law().strength_pa();
            candidate.cone_scale.push(
                positive(1.0 / response.max(material_length / material_force))?,
            );
            candidate.laws.push(law);
            candidate.states.push(state);
            candidate.transfers.push(transfer);
        }
        poll(&mut keep_going)?;
        Ok(candidate)
    }

    /// Publish only after the enclosing owner's complete audit, receiver
    /// admission and final cancellation poll. No fallible operation remains.
    pub(crate) fn publish_temperature_successors(
        &mut self,
        candidate: PreparedPatchTemperature<'a>,
    ) -> Vec<DamageFrictionState<'a>> {
        self.laws.copy_from_slice(&candidate.laws);
        self.cone_scale.copy_from_slice(&candidate.cone_scale);
        candidate.states
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, DamageFrictionPatchError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count)
        .map_err(|_| DamageFrictionPatchError::Invalid("temperature candidate allocation"))?;
    Ok(values)
}
