//! Original wall heat around the same accepted Cartesian pressure transaction.
use super::*;
use fs_qty::Time;

/// Thermal history and physical first-law receipt for one accepted interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirHeatReport {
    /// Original wall/bath state before the complete interval.
    pub before: WallThermal,
    /// Accepted finite wall, or unchanged prescribed-temperature reservoir.
    pub after: WallThermal,
    /// Sum of the two original implicit source receipts.
    pub source: WallCellReport,
    /// Tank plus finite-wall change minus port energy and external heat [J].
    pub gas_wall_energy_residual_j: f64,
}

/// A view is chosen by the caller, never by a mutable accepted/trial flag.
#[derive(Clone, Copy)]
pub(in crate::compressible) enum FiniteStage {
    Accepted,
    TransportStart,
    Trial,
}

impl<const D: usize, const S: usize> FiniteReservoirs<D, S> {
    pub(in crate::compressible) fn has_heat(&self) -> bool {
        self.heat_state.iter().any(Option::is_some)
    }

    pub(in crate::compressible) fn thermal(&self, slot: usize) -> Option<WallThermal> {
        self.heat_state[slot]
    }

    /// Dynamic storage in addition to the field's already counted metadata.
    pub(in crate::compressible) fn additional_heat_bytes(
        gas: &MixtureGas,
        count: usize,
    ) -> Result<usize, Error> {
        if count > S {
            return Err(Error::Invalid("too many finite mixture thermal sides"));
        }
        count
            .checked_mul(gas.species_count() * size_of::<f64>())
            .ok_or(Error::Invalid("finite mixture thermal workspace overflow"))
    }

    pub(in crate::compressible) fn bind_heat(
        &mut self,
        gas: &MixtureGas,
        slot: usize,
        thermal: WallThermal,
        control: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if slot >= S || self.heat_state[slot].is_some() {
            return Err(Error::Invalid(
                "unique finite mixture thermal side required",
            ));
        }
        let tank = self.accepted[slot].as_ref().ok_or(Error::Invalid(
            "finite mixture heat needs an already bound tank",
        ))?;
        let mut densities = Vec::new();
        densities
            .try_reserve_exact(gas.species_count())
            .map_err(|_| Error::Allocation)?;
        densities.extend_from_slice(tank.state().partial_densities);
        let bytes = self
            .allocated_bytes
            .checked_add(densities.capacity() * size_of::<f64>())
            .ok_or(Error::Invalid("finite mixture thermal workspace overflow"))?;
        let mut base = Plenum::new(
            tank.volume_m3(),
            Conserved {
                partial_densities: densities,
                momentum: 0.0,
                energy: tank.state().energy,
            },
        )?;
        // A zero-duration original-source call admits the same thermal law and
        // caloric data without changing either accepted inventory.
        tank.heat_candidate_into(
            gas,
            Time::new(0.0),
            thermal,
            control,
            &mut base,
            &mut *cancelled,
        )?;
        poll(cancelled)?;
        self.heat_base[slot] = Some(base);
        self.heat_state[slot] = Some(thermal);
        self.trial_heat_state[slot] = Some(thermal);
        self.heat_controls[slot] = Some(control);
        self.allocated_bytes = bytes;
        Ok(())
    }

    pub(in crate::compressible) fn state_at(
        &self,
        slot: usize,
        stage: FiniteStage,
    ) -> Option<Cell<'_>> {
        match stage {
            FiniteStage::Accepted => self.state(slot, false),
            FiniteStage::Trial => self.state(slot, true),
            FiniteStage::TransportStart => self.heat_base[slot]
                .as_ref()
                .or(self.accepted[slot].as_ref())
                .map(Plenum::state),
        }
    }

    /// Every complete transport attempt starts from accepted gas AND wall data.
    pub(in crate::compressible) fn begin_heat(
        &mut self,
        gas: &MixtureGas,
        duration: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if !self.has_heat() {
            return Ok(());
        }
        poll(cancelled)?;
        self.trial_heat_reports = [None; S];
        for slot in 0..S {
            let Some(thermal) = self.heat_state[slot] else {
                continue;
            };
            let update = self.accepted[slot]
                .as_ref()
                .expect("bound thermal tank")
                .heat_candidate_into(
                    gas,
                    Time::new(0.5 * duration),
                    thermal,
                    self.heat_controls[slot].expect("bound heat control"),
                    self.heat_base[slot]
                        .as_mut()
                        .expect("prepared thermal tank"),
                    &mut *cancelled,
                )?;
            self.trial_heat_state[slot] = Some(update.contact.thermal);
            self.trial_heat_reports[slot] = Some(update.report);
        }
        poll(cancelled)
    }

    /// The final port candidate supplies its actual mixed composition and Cp(T)
    /// to the second source. All outputs remain private until publish().
    pub(in crate::compressible) fn finish_heat(
        &mut self,
        gas: &MixtureGas,
        duration: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if !self.has_heat() {
            return Ok(());
        }
        poll(cancelled)?;
        for slot in 0..S {
            if self.heat_state[slot].is_none() {
                continue;
            }
            let mut report = *self.trial_reports[slot].first().ok_or(Error::Invalid(
                "finite tank heat needs the complete port candidate",
            ))?;
            let update = self.trial[slot]
                .as_ref()
                .expect("prepared port candidate")
                .heat_candidate_into(
                    gas,
                    Time::new(0.5 * duration),
                    self.trial_heat_state[slot].expect("prepared first heat half"),
                    self.heat_controls[slot].expect("bound heat control"),
                    self.heat_base[slot]
                        .as_mut()
                        .expect("prepared thermal output"),
                    &mut *cancelled,
                )?;
            let after = self.heat_base[slot]
                .as_ref()
                .expect("prepared thermal endpoint");
            let primitive = gas.primitive(after.state(), &mut *cancelled)?;
            let source = combine(
                self.trial_heat_reports[slot].expect("prepared first heat receipt"),
                update.report,
            )?;
            report.after = FiniteReservoirInventory::from_tank(after, primitive);
            report.energy_residual_j = report.after.energy
                - report.before.energy
                - report.transfer.energy
                - source.heat_to_gas_j;
            report.absolute_energy_residual_j = absolute_energy(
                gas,
                gas.species_count(),
                &report.species_residual_kg,
                report.energy_residual_j,
            )?;
            let coupled = report.energy_residual_j + source.energy_residual_j;
            if !report.energy_residual_j.is_finite() || !coupled.is_finite() {
                return Err(Error::Invalid(
                    "unrepresentable finite mixture gas-wall balance",
                ));
            }
            report.heat = Some(FiniteReservoirHeatReport {
                before: self.heat_state[slot].expect("accepted thermal state"),
                after: update.contact.thermal,
                source,
                gas_wall_energy_residual_j: coupled,
            });
            poll(cancelled)?;
            std::mem::swap(&mut self.trial[slot], &mut self.heat_base[slot]);
            self.trial_primitive[slot] = Some(primitive);
            self.trial_heat_state[slot] = Some(update.contact.thermal);
            self.trial_heat_reports[slot] = Some(source);
            self.trial_reports[slot].clear();
            self.trial_reports[slot].push(report);
        }
        poll(cancelled)
    }

    pub(super) fn heat_receipt(
        &self,
        slot: usize,
        energy_residual_j: f64,
    ) -> Option<FiniteReservoirHeatReport> {
        self.trial_heat_reports[slot].map(|source| FiniteReservoirHeatReport {
            before: self.heat_state[slot].expect("accepted thermal state"),
            after: self.trial_heat_state[slot].expect("trial thermal state"),
            source,
            gas_wall_energy_residual_j: energy_residual_j + source.energy_residual_j,
        })
    }
}

fn combine(a: WallCellReport, b: WallCellReport) -> Result<WallCellReport, Error> {
    let result = WallCellReport {
        elapsed: Time::new(a.elapsed.value() + b.elapsed.value()),
        wall_impulse_ns: a.wall_impulse_ns + b.wall_impulse_ns,
        mechanical_work_to_gas_j: a.mechanical_work_to_gas_j + b.mechanical_work_to_gas_j,
        drag_heating_j: a.drag_heating_j + b.drag_heating_j,
        heat_to_gas_j: a.heat_to_gas_j + b.heat_to_gas_j,
        wall_energy_change_j: a.wall_energy_change_j + b.wall_energy_change_j,
        external_heat_to_gas_j: a.external_heat_to_gas_j + b.external_heat_to_gas_j,
        energy_residual_j: a.energy_residual_j + b.energy_residual_j,
        thermal_residual_j: a.thermal_residual_j + b.thermal_residual_j,
        iterations: a
            .iterations
            .checked_add(b.iterations)
            .ok_or(Error::Invalid(
                "finite mixture thermal iteration count overflow",
            ))?,
    };
    if ![
        result.elapsed.value(),
        result.heat_to_gas_j,
        result.wall_energy_change_j,
        result.external_heat_to_gas_j,
        result.energy_residual_j,
        result.thermal_residual_j,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(Error::Invalid(
            "unrepresentable finite mixture thermal receipt",
        ));
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
