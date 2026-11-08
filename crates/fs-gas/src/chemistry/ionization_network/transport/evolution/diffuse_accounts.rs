//! Accepted diffuse traffic and direct-source/boundary audits.

use super::*;

impl Accounts {
    pub(super) fn add_diffuse(&mut self, increment: &Self) -> Result<(), Error> {
        add_rows(&mut self.diffuse_cells, &increment.diffuse_cells)?;
        add_rows(&mut self.faces, &increment.faces)
    }

    pub(super) fn integrate_diffuse(
        &mut self,
        field: &RadiationField<'_>,
        dt: f64,
    ) -> Result<(), Error> {
        for (c, target) in self.diffuse_cells.iter_mut().enumerate() {
            for (value, rate) in target.iter_mut().zip(field.diffuse_cell_rates(c)?) {
                *value = sum(*value, rate * Dd::from_f64(dt))?;
            }
        }
        for (target, rate) in self.faces.iter_mut().zip(field.interfaces()) {
            for (value, rate) in target.iter_mut().zip([
                rate.right_photons_per_s,
                rate.right_power_w,
                rate.left_photons_per_s,
                rate.left_power_w,
            ]) {
                *value = sum(*value, Dd::from_f64(product(rate, dt)?))?;
            }
        }
        Ok(())
    }

    pub(super) fn diffuse_error(
        &self,
        full: &Self,
        fine: &Self,
        controls: NetworkRadiationSolve,
    ) -> Result<f64, Error> {
        let mut error = 0.0_f64;
        for ((old, a), b) in self
            .diffuse_cells
            .iter()
            .zip(&full.diffuse_cells)
            .zip(&fine.diffuse_cells)
        {
            for j in 0..10 {
                let unit = if j < 4 { 1.0 / NA } else { 1.0 };
                let absolute = if j < 4 || j == 7 {
                    controls.network.local.absolute_amount_tolerance_mol
                } else {
                    controls.network.local.energy_tolerance_j
                };
                error = error.max(moment_error(old[j], a[j], b[j], unit, absolute, controls)?);
            }
        }
        for ((old, a), b) in self.faces.iter().zip(&full.faces).zip(&fine.faces) {
            for j in 0..4 {
                let (unit, absolute) = if j % 2 == 0 {
                    (
                        1.0 / NA,
                        controls.network.local.absolute_amount_tolerance_mol,
                    )
                } else {
                    (1.0, controls.network.local.energy_tolerance_j)
                };
                error = error.max(moment_error(old[j], a[j], b[j], unit, absolute, controls)?);
            }
        }
        Ok(error)
    }

    pub(super) fn publish_diffuse(
        &self,
        field: &RadiationField<'_>,
        gas: NetworkRadiationBalance,
        controls: NetworkRadiationSolve,
    ) -> Result<Option<NetworkDiffuseEvolutionSummary>, Error> {
        let Some(diffuse_controls) = field.diffuse_controls() else {
            return Ok(None);
        };
        let mut interfaces = reserved(self.faces.len())?;
        for face in &self.faces {
            let a = face.map(Dd::to_f64);
            if a.iter().any(|v| !nonnegative(*v)) {
                return Err(Error::Arithmetic("integrated diffuse interface"));
            }
            interfaces.push(DiffuseRecombinationExposure {
                right_photons: a[0],
                right_energy_j: a[1],
                left_photons: a[2],
                left_energy_j: a[3],
            });
        }
        let first = self
            .faces
            .first()
            .ok_or(Error::Arithmetic("missing left radiation boundary"))?;
        let last = self
            .faces
            .last()
            .ok_or(Error::Arithmetic("missing right radiation boundary"))?;
        let incoming = sum(first[0], last[2])?;
        let incoming_energy = sum(first[1], last[3])?;
        let outgoing = sum(first[2], last[0])?;
        let outgoing_energy = sum(first[3], last[1])?;
        let mut totals = [Dd::ZERO; 10];
        let mut net_extent_magnitude = Dd::ZERO;
        for cell in &self.diffuse_cells {
            for (total, value) in totals.iter_mut().zip(cell) {
                *total = sum(*total, *value)?;
            }
            net_extent_magnitude = sum(net_extent_magnitude, cell[7].abs())?;
        }
        let local = controls.network.local;
        let photon_residual = (totals[7] * Dd::from_f64(NA) + outgoing - incoming).to_f64();
        let photon_allowance = NA * local.absolute_amount_tolerance_mol
            + 128.0
                * f64::EPSILON
                * (NA * net_extent_magnitude.to_f64() + outgoing.to_f64() + incoming.to_f64());
        if !photon_residual.is_finite()
            || !nonnegative(photon_allowance)
            || photon_residual.abs() > photon_allowance
        {
            return Err(Error::Arithmetic("integrated diffuse photon balance"));
        }
        let boundary_roundoff = 128.0
            * f64::EPSILON
            * (incoming_energy + outgoing_energy + Dd::from_f64(gas.free_free_escaped_energy_j))
                .to_f64();
        let radiation_residual = (totals[9] + outgoing_energy - incoming_energy).to_f64();
        let gas_residual = (Dd::from_f64(gas.gas_energy_change_j)
            + outgoing_energy
            + Dd::from_f64(gas.free_free_escaped_energy_j)
            - incoming_energy
            - Dd::from_f64(gas.external_work_j))
        .to_f64();
        energy_check(
            radiation_residual,
            local.energy_tolerance_j + boundary_roundoff,
        )?;
        energy_check(gas_residual, local.energy_tolerance_j + boundary_roundoff)?;
        // Independent original chemical ledgers must describe the same captures
        // and direct source. Gross arithmetic allowances never relax gas storage.
        let emitted_delta = (totals[0] - Dd::from_f64(gas.emitted_photons)).to_f64();
        let emitted_allowance = NA * local.absolute_amount_tolerance_mol
            + 128.0 * f64::EPSILON * totals[0].to_f64().max(gas.emitted_photons);
        if !emitted_delta.is_finite()
            || !nonnegative(emitted_allowance)
            || emitted_delta.abs() > emitted_allowance
        {
            return Err(Error::Arithmetic(
                "network and diffuse capture counts disagree",
            ));
        }
        energy_check(
            (totals[5] - Dd::from_f64(gas.emitted_radiation_energy_j)).to_f64(),
            local.energy_tolerance_j
                + 128.0 * f64::EPSILON * totals[5].to_f64().max(gas.emitted_radiation_energy_j),
        )?;
        energy_check(
            (totals[9] - Dd::from_f64(gas.net_radiation_energy_j)).to_f64(),
            local.energy_tolerance_j,
        )?;
        let [spectral, angular, tail] = field.maxima();
        Ok(Some(NetworkDiffuseEvolutionSummary {
            interfaces,
            controls: diffuse_controls,
            balance: NetworkDiffuseEvolutionBalance {
                emitted_photons: totals[0].to_f64(),
                emitted_energy_j: totals[5].to_f64(),
                diffuse_absorbed_photons: (totals[1] + totals[2] + totals[3]).to_f64(),
                diffuse_absorbed_energy_j: totals[4].to_f64(),
                incident_photons: incoming.to_f64(),
                incident_energy_j: incoming_energy.to_f64(),
                escaped_photons: outgoing.to_f64(),
                escaped_energy_j: outgoing_energy.to_f64(),
                free_free_escaped_energy_j: gas.free_free_escaped_energy_j,
                net_ionization_extent_mol: totals[7].to_f64(),
                net_radiation_energy_j: totals[9].to_f64(),
                photon_residual,
                radiation_energy_residual_j: radiation_residual,
                energy_residual_j: gas_residual,
                boundary_energy_roundoff_j: boundary_roundoff,
                radiation_evaluations: field.evaluations(),
                radiation_work: field.work(),
                maximum_normalized_spectral_estimate: spectral,
                maximum_normalized_angular_estimate: angular,
                maximum_normalized_tail_bound: tail,
            },
        }))
    }
}

fn add_rows<const N: usize>(target: &mut [[Dd; N]], source: &[[Dd; N]]) -> Result<(), Error> {
    for (target, source) in target.iter_mut().zip(source) {
        for (value, increment) in target.iter_mut().zip(source) {
            *value = sum(*value, *increment)?;
        }
    }
    Ok(())
}

fn moment_error(
    old: Dd,
    full: Dd,
    fine: Dd,
    unit: f64,
    absolute: f64,
    controls: NetworkRadiationSolve,
) -> Result<f64, Error> {
    let difference = (full - fine).to_f64().abs() * unit;
    let scale = (old + full).to_f64().abs().max((old + fine).to_f64().abs()) * unit;
    let error = difference / (3.0 * (absolute + controls.network.local.relative_tolerance * scale));
    if !nonnegative(difference) || !nonnegative(scale) || !nonnegative(error) {
        Err(Error::Arithmetic("diffuse exposure local error"))
    } else {
        Ok(error)
    }
}
