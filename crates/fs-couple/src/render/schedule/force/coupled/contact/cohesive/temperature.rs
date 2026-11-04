//! Explicit temperature continuation of the original native cohesive interface.
//!
//! The source driver changes coefficients at fixed accepted separation, then
//! the original material owner admits any irreversible damage release. Native
//! coordinates, material/Maxwell history, work columns, time and the last
//! mechanical observation remain untouched. Signed parameter work belongs to
//! that driver; this operation does not assign fracture release to heat.

use super::*;
use crate::cohesive_contact::reclosure::ResidualCoulombSource;
use fs_material::cohesive::temperature::CohesiveTemperatureStep;
use fs_math::dd::Dd;

/// One accepted source change at the existing native sample and separation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveModalTemperatureReport {
    pub sample: u64,
    pub jump_m: [f64; 3],
    pub damage_before: f64,
    pub damage_after: f64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub cohesive_energy_before_j: f64,
    pub cohesive_energy_after_j: f64,
    /// Signed coefficient-change work supplied by the temperature driver [J].
    pub parameter_work_j: f64,
    /// Additional irreversible release at the held separation [J].
    pub damage_dissipation_j: f64,
    /// Complete interface storage change, including the new damage [J].
    pub stored_change_j: f64,
    /// Storage change plus new fracture release minus parameter work [J].
    pub energy_residual_j: f64,
    /// Interface-sized tolerance, independent of unrelated native energy [J].
    pub energy_tolerance_j: f64,
}

struct PreparedTemperature<'a> {
    point: CohesivePoint<'a>,
    residual: Option<reclosure::ReclosureBinding<'a>>,
    damage_friction: Option<frictional_damage::DamageFrictionBinding<'a>>,
}

impl<'a> CohesiveModalSystem<'a> {
    /// Continue the original source at a prescribed temperature without a
    /// mechanical timestep. The returned signed parameter work must be supplied
    /// by the explicit driver. Cooling cannot heal the accepted bond.
    ///
    /// An enabled friction mode requires the paired method below, even if its
    /// coefficients happen to be unchanged: both sources must retain the same
    /// original card and exact external query.
    pub fn apply_temperature_successor(
        &mut self,
        next: &'a ResolvedBilinearCohesive,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<CohesiveModalTemperatureReport, CohesiveModalError> {
        self.apply_temperature_successors_with_admission(
            next,
            None,
            temperature_axis,
            gate,
            |report, _| Ok(*report),
        )
    }

    /// Continue the existing cohesive and residual/damaged-area friction
    /// sources together. Contact occupancy is retained exactly; a coefficient
    /// change supplies neither an inferred slip velocity nor frictional work.
    pub fn apply_temperature_and_friction_successors(
        &mut self,
        next: &'a ResolvedBilinearCohesive,
        next_friction: &'a ResidualCoulombSource,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<CohesiveModalTemperatureReport, CohesiveModalError> {
        self.apply_temperature_successors_with_admission(
            next,
            Some(next_friction),
            temperature_axis,
            gate,
            |report, _| Ok(*report),
        )
    }

    /// Stage a receiving owner's account before publishing any source or
    /// irreversible history. The callback returns a candidate, not an accepted
    /// receiver mutation; its installation after success must be infallible.
    /// Refusal or cancellation leaves every accepted native/interface state
    /// unchanged. No native body, propagator or observation is cloned.
    pub fn apply_temperature_successors_with_admission<R, E>(
        &mut self,
        next: &'a ResolvedBilinearCohesive,
        next_friction: Option<&'a ResidualCoulombSource>,
        temperature_axis: &str,
        gate: &CancelGate,
        admit: impl FnOnce(&CohesiveModalTemperatureReport, CohesivePoint<'a>) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<CohesiveModalError>,
    {
        let (candidate, report) = self
            .prepare_temperature_successors(next, next_friction, temperature_axis, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, candidate.point)?;
        poll(Some(gate))
            .map_err(CohesiveModalError::from)
            .map_err(E::from)?;
        self.point = candidate.point;
        self.source = Some(next);
        self.residual = candidate.residual;
        self.damage_friction = candidate.damage_friction;
        Ok(receiving)
    }

    fn prepare_temperature_successors(
        &self,
        next: &'a ResolvedBilinearCohesive,
        next_friction: Option<&'a ResidualCoulombSource>,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<(PreparedTemperature<'a>, CohesiveModalTemperatureReport), CohesiveModalError> {
        poll(Some(gate))?;
        let previous = self.source.ok_or_else(|| {
            invalid("temperature continuation requires the original cohesive source")
        })?;
        previous.require_temperature_successor(next, temperature_axis)?;
        if !core::ptr::eq(self.point.law(), previous.law()) {
            return Err(invalid("temperature continuation must retain the accepted source point").into());
        }
        let (transfer, next_residual, next_damage_friction): (
            CohesiveTemperatureStep<'a>,
            _,
            _,
        ) = match (&self.residual, &self.damage_friction, next_friction) {
            (None, None, None) => (
                self.point.transfer_temperature(
                    next.law(),
                    self.config.contact.material_energy_tolerance_j_m2,
                )?,
                None,
                None,
            ),
            (Some(binding), None, Some(friction)) => {
                let (candidate, transfer) =
                    binding.transfer_temperature(self.point, next, friction, temperature_axis)?;
                (transfer, Some(candidate), None)
            }
            (None, Some(binding), Some(friction)) => {
                let (candidate, transfer) =
                    binding.transfer_temperature(self.point, next, friction, temperature_axis)?;
                (transfer, None, Some(candidate))
            }
            (None, None, Some(_)) => {
                return Err(invalid("a temperature update cannot enable an absent friction mode").into());
            }
            _ => {
                return Err(invalid("temperature continuation requires the paired original friction source").into());
            }
        };
        if finite(transfer.defect_j_m2)?.abs()
            > self.config.contact.material_energy_tolerance_j_m2
        {
            return Err(invalid("temperature transfer exceeds the original material energy tolerance").into());
        }
        let area = self.area_m2();
        let cohesive_before = finite(area * self.point.stored_energy_j_m2())?;
        let cohesive_after = finite(area * transfer.point.stored_energy_j_m2())?;
        let work = finite(area * transfer.parameter_work_j_m2)?;
        let release = finite(area * transfer.damage_dissipated_j_m2)?;
        if release < 0.0 || transfer.point.damage() < self.point.damage() {
            return Err(invalid("temperature continuation cannot heal or create fracture energy").into());
        }
        let change = Dd::from_f64(cohesive_after) - Dd::from_f64(cohesive_before);
        let stored_change = finite(change.to_f64())?;
        let residual = finite(
            (change + Dd::from_f64(release) - Dd::from_f64(work)).to_f64(),
        )?;
        let scale = cohesive_before
            .abs()
            .max(cohesive_after.abs())
            .max(stored_change.abs())
            .max(work.abs())
            .max(release);
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let native = self.network.total_energy_j()?;
        let energy_before = finite(native + cohesive_before)?;
        let energy_after = finite(native + cohesive_after)?;
        limit(
            "temperature-refreshed native cohesive energy",
            energy_after,
            self.network.config.maximum_total_energy_j,
        )?;
        let report = CohesiveModalTemperatureReport {
            sample: self.samples_rendered(),
            jump_m: self.point.jump_m(),
            damage_before: self.point.damage(),
            damage_after: transfer.point.damage(),
            energy_before_j: energy_before,
            energy_after_j: energy_after,
            cohesive_energy_before_j: cohesive_before,
            cohesive_energy_after_j: cohesive_after,
            parameter_work_j: work,
            damage_dissipation_j: release,
            stored_change_j: stored_change,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(Some(gate))?;
        Ok((
            PreparedTemperature {
                point: transfer.point,
                residual: next_residual,
                damage_friction: next_damage_friction,
            },
            report,
        ))
    }
}
