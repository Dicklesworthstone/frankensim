//! Explicit source-coefficient changes at fixed accepted native motion.
//!
//! The original ordered interface must authorize the temperature successor.
//! Existing branch strains and elapsed mechanical time survive; signed
//! coefficient-change work belongs to the prescribed temperature driver.
//! Neither body motion nor an implicit heat source is introduced here.
use super::*;

/// An accepted original-source update without a mechanical timestep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticBondModalTemperatureReport {
    pub sample: u64,
    pub elapsed_s: f64,
    pub jump_m: [f64; 3],
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub bond_energy_before_j: f64,
    pub bond_energy_after_j: f64,
    /// Signed reversible coefficient-change work supplied by the driver [J].
    /// Negative work returns energy to that driver; it is not relaxation heat.
    pub parameter_work_j: f64,
    pub stored_change_j: f64,
    pub energy_residual_j: f64,
    /// Bond-sized tolerance independent of unrelated native stored energy [J].
    pub energy_tolerance_j: f64,
}

impl<'a> ViscoelasticBondModalSystem<'a> {
    /// Refresh only the original source's temperature coordinate. Its pinned
    /// coefficient claims, surfaces, geometry, branch order and runtime limits
    /// must remain the same. Rebuild the cached ramps for the existing native
    /// timestep while retaining accepted branch strains and time.
    ///
    /// The last mechanical frame remains the last mechanical observation.
    /// This operation has its own report and does not advance a body or sample.
    pub fn apply_temperature_successor(
        &mut self,
        next: &'a ResolvedViscoelasticBond,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<ViscoelasticBondModalTemperatureReport, CohesiveModalError> {
        self.apply_temperature_successor_with_admission(
            next,
            temperature_axis,
            gate,
            |report, _| Ok(*report),
        )
    }

    /// Stage the receiving driver's signed work account before any source,
    /// history or prepared coefficients are published. The callback returns a
    /// candidate to install infallibly only after this method succeeds.
    pub fn apply_temperature_successor_with_admission<R, E>(
        &mut self,
        next: &'a ResolvedViscoelasticBond,
        temperature_axis: &str,
        gate: &CancelGate,
        admit: impl FnOnce(
            &ViscoelasticBondModalTemperatureReport,
            &ViscoelasticBondPoint<'a>,
        ) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<CohesiveModalError>,
    {
        let (point, prepared, report) = self
            .prepare_temperature_successor(next, temperature_axis, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, &point)?;
        poll(Some(gate))
            .map_err(CohesiveModalError::from)
            .map_err(E::from)?;
        self.point = point;
        self.prepared = prepared;
        self.source = next;
        Ok(receiving)
    }

    fn prepare_temperature_successor(
        &self,
        next: &'a ResolvedViscoelasticBond,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<
        (
            ViscoelasticBondPoint<'a>,
            ViscoelasticBondRamp<'a>,
            ViscoelasticBondModalTemperatureReport,
        ),
        CohesiveModalError,
    > {
        poll(Some(gate))?;
        self.source
            .require_temperature_successor(next, temperature_axis)?;
        if !core::ptr::eq(self.point.law(), self.source.law()) {
            return Err(invalid(
                "temperature continuation must retain the original relaxing bond point",
            )
            .into());
        }
        let transfer = self.point.transfer_temperature(
            next.law(),
            self.config.contact.material_energy_tolerance_j_m2,
        )?;
        let prepared = next.law().prepare(self.network.dt)?;
        let area = self.area_m2();
        let before = finite(area * self.point.stored_energy_j_m2())?;
        let after = finite(area * transfer.point.stored_energy_j_m2())?;
        let parameter_work = finite(area * transfer.parameter_work_j_m2)?;
        // Retain the constitutive coefficient increments even when they are
        // small beside the original elastic preload; do not subtract stores.
        let stored_change = finite(area * transfer.stored_change_j_m2)?;
        let scale = stored_change.abs().max(parameter_work.abs());
        let tolerance = finite(
            area * self.config.contact.material_energy_tolerance_j_m2
                + 128.0 * f64::EPSILON * scale,
        )?;
        check_work(finite(area * transfer.defect_j_m2)?, tolerance)?;
        let residual =
            finite((Dd::from_f64(stored_change) - Dd::from_f64(parameter_work)).to_f64())?;
        check_work(residual, tolerance)?;
        let native = self.network.total_energy_j()?;
        let energy_before = finite(native + before)?;
        let energy_after = finite(native + after)?;
        limit(
            "temperature-refreshed native relaxing bond energy",
            energy_after,
            self.network.config.maximum_total_energy_j,
        )?;
        let report = ViscoelasticBondModalTemperatureReport {
            sample: self.samples_rendered(),
            elapsed_s: self.point.elapsed_s(),
            jump_m: self.point.jump_m(),
            energy_before_j: energy_before,
            energy_after_j: energy_after,
            bond_energy_before_j: before,
            bond_energy_after_j: after,
            parameter_work_j: parameter_work,
            stored_change_j: stored_change,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(Some(gate))?;
        Ok((transfer.point, prepared, report))
    }
}
