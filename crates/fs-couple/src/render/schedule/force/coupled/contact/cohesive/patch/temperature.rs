//! Atomic temperature continuation of every original finite-patch site.
//!
//! Coefficients change at held local jumps, followed by the material owner's
//! irreversible relaxation. Per-site signed parameter work and new fracture
//! release remain separate. This does not advance time or assign release to heat.
//! Native mechanics, Maxwell memory, full cross-compliance and observations stay
//! in place; only prepared source bindings and canonical point histories publish.

use super::*;
use crate::cohesive_contact::patch::frictional_damage::temperature::PreparedPatchTemperature;
use crate::cohesive_contact::reclosure::ResidualCoulombSource;
use fs_math::dd::Dd;

/// One original quadrature site's fixed-jump source change; energies include
/// its own reference area exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesivePatchPointTemperatureReport {
    pub area_m2: f64,
    pub jump_m: [f64; 3],
    pub damage_before: f64,
    pub damage_after: f64,
    pub cohesive_energy_before_j: f64,
    pub cohesive_energy_after_j: f64,
    pub parameter_work_j: f64,
    pub damage_dissipation_j: f64,
    pub stored_change_j: f64,
    /// Storage change plus fracture release minus parameter work [J].
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

/// Complete source update at one unchanged native sample. The tolerance scales
/// with the interface accounts, independently of unrelated body energy.
#[derive(Clone, Debug, PartialEq)]
pub struct CohesivePatchModalTemperatureReport {
    pub sample: u64,
    /// Original site order, including heterogeneous areas and temperatures.
    pub points: Vec<CohesivePatchPointTemperatureReport>,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub cohesive_energy_before_j: f64,
    pub cohesive_energy_after_j: f64,
    /// Signed sum supplied by the explicit temperature driver [J].
    pub parameter_work_j: f64,
    /// Sum of new irreversible releases at held separations [J].
    pub damage_dissipation_j: f64,
    pub stored_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

struct PreparedTemperature<'a> {
    points: Vec<CohesivePoint<'a>>,
    friction: Option<PreparedPatchTemperature<'a>>,
}

impl<'a> CohesivePatchModalSystem<'a> {
    /// Continue each site's original card/query at an explicitly prescribed
    /// temperature. Supply exactly one successor per original site in order.
    /// Enabled damaged-area friction requires the paired method below.
    pub fn apply_temperature_successors(
        &mut self,
        next: &[&'a ResolvedBilinearCohesive],
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<CohesivePatchModalTemperatureReport, CohesivePatchModalError> {
        self.apply_temperature_successors_with_admission(
            next, None, temperature_axis, gate,
            |_, _| Ok::<(), CohesivePatchModalError>(()),
        ).map(|(report, ())| report)
    }

    /// Continue both original sources at every site. Each cohesive/friction pair
    /// must share its exact new query; different sites may have different
    /// temperatures. Friction occupancy is retained without inferring velocity,
    /// slip work, a new bond, or a new thermal reservoir.
    pub fn apply_temperature_and_friction_successors(
        &mut self,
        next: &[&'a ResolvedBilinearCohesive],
        next_friction: &[&'a ResidualCoulombSource],
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<CohesivePatchModalTemperatureReport, CohesivePatchModalError> {
        self.apply_temperature_successors_with_admission(
            next, Some(next_friction), temperature_axis, gate,
            |_, _| Ok::<(), CohesivePatchModalError>(()),
        ).map(|(report, ())| report)
    }

    /// Prepare a receiving owner's candidate account before any source or
    /// history publishes. The callback must return a candidate, not mutate an
    /// accepted receiver; installation after success must be infallible.
    /// Refusal or late cancellation preserves every accepted body, site,
    /// friction occupancy, clock and previous mechanical observation.
    ///
    /// Returns the complete report together with the admitted receiver
    /// candidate. Candidate storage is O(sites); no native network, dense solver or full
    /// compliance is cloned or reconstructed.
    pub fn apply_temperature_successors_with_admission<R, E>(
        &mut self,
        next: &[&'a ResolvedBilinearCohesive],
        next_friction: Option<&[&'a ResidualCoulombSource]>,
        temperature_axis: &str,
        gate: &CancelGate,
        admit: impl FnOnce(&CohesivePatchModalTemperatureReport, &[CohesivePoint<'a>]) -> Result<R, E>,
    ) -> Result<(CohesivePatchModalTemperatureReport, R), E>
    where
        E: From<CohesivePatchModalError>,
    {
        let (candidate, report) = self
            .prepare_temperature_successors(next, next_friction, temperature_axis, gate)
            .map_err(E::from)?;
        let receiver = admit(&report, &candidate.points)?;
        poll(Some(gate)).map_err(CohesivePatchModalError::from).map_err(E::from)?;
        // Everything below is infallible. The candidate was prepared against
        // these same immutable bindings and original site count.
        if let Some(friction) = candidate.friction {
            self.damage_friction.as_mut()
                .expect("prepared original patch friction")
                .publish_temperature_successors(friction);
        }
        self.points.copy_from_slice(&candidate.points);
        self.staged_points.copy_from_slice(&candidate.points);
        for (site, source) in self.sites.iter_mut().zip(next) {
            site.law = source.law();
            site.source = Some(*source);
        }
        Ok((report, receiver))
    }

    fn prepare_temperature_successors(
        &self,
        next: &[&'a ResolvedBilinearCohesive],
        next_friction: Option<&[&'a ResidualCoulombSource]>,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<(PreparedTemperature<'a>, CohesivePatchModalTemperatureReport), CohesivePatchModalError> {
        poll(Some(gate))?;
        let n = self.points.len();
        if next.len() != n || next_friction.is_some_and(|sources| sources.len() != n) {
            return Err(invalid("temperature continuation requires every original patch site").into());
        }
        for p in 0..n {
            poll(Some(gate))?;
            let previous = self.sites[p].source.ok_or_else(|| {
                invalid("temperature continuation requires every original cohesive source")
            })?;
            previous.require_temperature_successor(next[p], temperature_axis)?;
            if !core::ptr::eq(self.sites[p].law, previous.law())
                || !core::ptr::eq(self.points[p].law(), previous.law())
            {
                return Err(invalid("temperature continuation must retain every accepted source point").into());
            }
        }
        let material_tolerance = self.config.patch.contact.material_energy_tolerance_j_m2;
        let mut unpaired_transfers = Vec::new();
        let friction = match (&self.damage_friction, next_friction) {
            (None, None) => {
                unpaired_transfers.try_reserve_exact(n)
                    .map_err(|_| invalid("patch temperature material allocation refused"))?;
                for p in 0..n {
                    poll(Some(gate))?;
                    unpaired_transfers.push(
                        self.points[p].transfer_temperature(next[p].law(), material_tolerance)
                            .map_err(|source| CohesivePatchModalError::Material { point: p, source })?,
                    );
                }
                None
            }
            (Some(binding), Some(sources)) => Some(binding.prepare_temperature_successors(
                &self.points, next, sources, temperature_axis, material_tolerance, gate,
            )?),
            (Some(_), None) => return Err(invalid(
                "patch temperature continuation requires every paired original friction source",
            ).into()),
            (None, Some(_)) => return Err(invalid(
                "a patch temperature update cannot enable absent friction",
            ).into()),
        };
        let transfers = friction.as_ref()
            .map_or(unpaired_transfers.as_slice(), PreparedPatchTemperature::transfers);
        let mut points = reserved(n)?;
        let mut point_reports = reserved(n)?;
        let mut before = Dd::from_f64(0.0);
        let mut after = Dd::from_f64(0.0);
        let mut parameter = Dd::from_f64(0.0);
        let mut absolute_parameter = Dd::from_f64(0.0);
        let mut release = Dd::from_f64(0.0);
        for p in 0..n {
            poll(Some(gate))?;
            let transfer = transfers[p];
            if finite(transfer.defect_j_m2)?.abs() > material_tolerance {
                return Err(invalid("patch temperature transfer exceeds the original material energy tolerance").into());
            }
            if transfer.point.jump_m() != self.points[p].jump_m()
                || !core::ptr::eq(transfer.point.law(), next[p].law())
                || transfer.point.damage() < self.points[p].damage()
                || transfer.point.dissipated_j_m2() < self.points[p].dissipated_j_m2()
            {
                return Err(invalid("patch temperature transfer must retain the original jump and irreversible history").into());
            }
            let area = self.sites[p].area_m2;
            let local_before = finite(area * self.points[p].stored_energy_j_m2())?;
            let local_after = finite(area * transfer.point.stored_energy_j_m2())?;
            let work = finite(area * transfer.parameter_work_j_m2)?;
            let damage = finite(area * transfer.damage_dissipated_j_m2)?;
            if damage < 0.0 {
                return Err(invalid("patch temperature transfer cannot create fracture energy").into());
            }
            let change = Dd::from_f64(local_after) - Dd::from_f64(local_before);
            let stored_change = finite(change.to_f64())?;
            let residual = finite((change + Dd::from_f64(damage) - Dd::from_f64(work)).to_f64())?;
            let scale = local_before.abs().max(local_after.abs())
                .max(stored_change.abs()).max(work.abs()).max(damage);
            let tolerance = self.temperature_tolerance(scale)?;
            check_balance(residual, tolerance)?;
            point_reports.push(CohesivePatchPointTemperatureReport {
                area_m2: area,
                jump_m: self.points[p].jump_m(),
                damage_before: self.points[p].damage(),
                damage_after: transfer.point.damage(),
                cohesive_energy_before_j: local_before,
                cohesive_energy_after_j: local_after,
                parameter_work_j: work,
                damage_dissipation_j: damage,
                stored_change_j: stored_change,
                energy_residual_j: residual,
                energy_tolerance_j: tolerance,
            });
            points.push(transfer.point);
            before = before + Dd::from_f64(local_before);
            after = after + Dd::from_f64(local_after);
            parameter = parameter + Dd::from_f64(work);
            absolute_parameter = absolute_parameter + Dd::from_f64(work.abs());
            release = release + Dd::from_f64(damage);
        }
        let cohesive_before = finite(before.to_f64())?;
        let cohesive_after = finite(after.to_f64())?;
        let work = finite(parameter.to_f64())?;
        let damage = finite(release.to_f64())?;
        let change = after - before;
        let stored_change = finite(change.to_f64())?;
        let residual = finite((change + release - parameter).to_f64())?;
        let scale = cohesive_before.abs().max(cohesive_after.abs())
            .max(stored_change.abs()).max(finite(absolute_parameter.to_f64())?).max(damage);
        let tolerance = self.temperature_tolerance(scale)?;
        check_balance(residual, tolerance)?;
        let native = self.network.total_energy_j()?;
        let energy_before = finite(native + cohesive_before)?;
        let energy_after = finite(native + cohesive_after)?;
        limit("temperature-refreshed native patch energy", energy_after,
            self.network.config.maximum_total_energy_j)?;
        let report = CohesivePatchModalTemperatureReport {
            sample: self.samples_rendered(),
            points: point_reports,
            energy_before_j: energy_before,
            energy_after_j: energy_after,
            cohesive_energy_before_j: cohesive_before,
            cohesive_energy_after_j: cohesive_after,
            parameter_work_j: work,
            damage_dissipation_j: damage,
            stored_change_j: stored_change,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(Some(gate))?;
        Ok((PreparedTemperature { points, friction }, report))
    }

    fn temperature_tolerance(&self, scale: f64) -> Result<f64, CohesivePatchModalError> {
        Ok(finite(self.network.config.energy_absolute_tolerance_j
            + self.network.config.energy_relative_tolerance * scale)?)
    }
}

fn check_balance(residual: f64, tolerance: f64) -> Result<(), CohesivePatchModalError> {
    if residual.abs() > tolerance {
        Err(ModalCouplingError::EnergyBalance {
            residual_j: residual, tolerance_j: tolerance,
        }.into())
    } else {
        Ok(())
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, CohesivePatchModalError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count)
        .map_err(|_| invalid("patch temperature candidate allocation refused"))?;
    Ok(values)
}
