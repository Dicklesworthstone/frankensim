use super::*;
use crate::chemistry::collisional_radiative::solve::{Point, Problem, Rates};
use crate::chemistry::photoionization::{depletion, transmitted_photons};
use fs_math::dd::Dd;
use fs_math::det::{exp, ln};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;

type Error = RadiationTransportError;

#[derive(Clone, Copy)]
struct Band {
    opacity_per_mol: f64,
    excess_j_mol: f64,
}
struct Column<'a> {
    problems: Vec<Problem<'a>>,
    groups: &'a [PhotonNumberFluxV1],
    bands: Vec<Band>,
    boundary: Vec<[f64; 2]>, // photon rate, complete-exposure photon count
    controls: RadiationTransportSolve,
}

pub(super) fn advance(
    cells: &[RadiationTransportCell<'_>],
    groups: &[PhotonNumberFluxV1],
    area: Area,
    duration: Time,
    controls: RadiationTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RadiationTransportResult, Error> {
    let column = Column::new(cells, groups, area, duration, controls, &mut cancelled)?;
    let mut points = buffer(cells.len(), Point::default())?;
    let mut statistics = buffer(cells.len(), CollisionalRadiativeBalance::default())?;
    let mut rates = Vec::new();
    rates
        .try_reserve_exact(cells.len())
        .map_err(|_| Error::Arithmetic("reaction rate allocation"))?;
    for (c, problem) in column.problems.iter().enumerate() {
        points[c] = problem.initial_point();
        rates.push(
            problem
                .rates(points[c], &mut statistics[c], &mut cancelled)
                .map_err(|error| cell_error(c, error))?,
        );
    }
    let mut full = buffer(cells.len(), Point::default())?;
    let mut fine = buffer(cells.len(), Point::default())?;
    let mut cell_errors = buffer(cells.len(), 0.0_f64)?;
    let pairs = column.bands.len();
    let mut accounts = buffer(pairs, [Dd::ZERO; 2])?;
    let mut full_accounts = buffer(pairs, [Dd::ZERO; 2])?;
    let mut fine_accounts = buffer(pairs, [Dd::ZERO; 2])?;
    let (mut time, mut step) = (0.0, controls.chemistry.initial_step_s);
    let (mut accepted, mut rejected) = (0, 0);
    while time < duration.value() {
        poll(&mut cancelled)?;
        if column.inactive(&points, &rates) {
            // No reaction or thermal work remains. Preserve the complete
            // remaining incident inventory, including subthreshold bands.
            for (g, boundary) in column.boundary.iter().enumerate() {
                let photons = boundary[0] * (duration.value() - time);
                if !nonnegative(photons) {
                    return Err(Error::Arithmetic("remaining photon exposure"));
                }
                for c in 0..cells.len() {
                    poll(&mut cancelled)?;
                    let account = &mut accounts[g * cells.len() + c];
                    account[1] = account[1] + Dd::from_f64(photons);
                }
            }
            break;
        }
        if accepted + rejected >= controls.chemistry.maximum_steps {
            return Err(Error::StepBudget);
        }
        let remaining = duration.value() - time;
        let dt = step.min(remaining);
        let end = if dt == remaining {
            duration.value()
        } else {
            time + dt
        };
        let dt = end - time;
        if end <= time || dt * 0.5 == 0.0 {
            return Err(Error::StepBudget);
        }
        full.copy_from_slice(&points);
        fine.copy_from_slice(&points);
        full_accounts.fill([Dd::ZERO; 2]);
        fine_accounts.fill([Dd::ZERO; 2]);
        cell_errors.fill(0.0);
        let trial = (|| {
            column.split_step(
                &mut full,
                Some(&rates),
                dt,
                &mut full_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            column.split_step(
                &mut fine,
                Some(&rates),
                dt * 0.5,
                &mut fine_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            column.split_step(
                &mut fine,
                None,
                dt * 0.5,
                &mut fine_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            let mut error = 0.0_f64;
            for (c, problem) in column.problems.iter().enumerate() {
                poll(&mut cancelled)?;
                cell_errors[c] = problem
                    .error(points[c], full[c], fine[c])
                    .map_err(|error| cell_error(c, error))?;
                error = error.max(cell_errors[c]);
            }
            for (k, ((old, coarse), next)) in accounts
                .iter()
                .zip(&full_accounts)
                .zip(&fine_accounts)
                .enumerate()
            {
                if k % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                for j in 0..2 {
                    // Photon integrals are cumulative ODE states. Subtract the
                    // local Dd increments before scaling to retain small errors.
                    let difference = (coarse[j] - next[j]).to_f64().abs() / NA;
                    let scale = (old[j] + coarse[j])
                        .to_f64()
                        .abs()
                        .max((old[j] + next[j]).to_f64().abs())
                        / NA;
                    let local = difference
                        / (3.0
                            * (controls.chemistry.absolute_amount_tolerance_mol
                                + controls.chemistry.relative_tolerance * scale));
                    if !nonnegative(local) {
                        return Err(Error::Arithmetic("band local error"));
                    }
                    error = error.max(local);
                    cell_errors[k % cells.len()] = cell_errors[k % cells.len()].max(local);
                }
            }
            Ok(error)
        })();
        match trial {
            Ok(error) if error <= 1.0 => {
                for (k, (total, local)) in accounts.iter_mut().zip(&fine_accounts).enumerate() {
                    if k % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    for j in 0..2 {
                        total[j] = total[j] + local[j];
                    }
                }
                core::mem::swap(&mut points, &mut fine);
                time = end;
                accepted += 1;
                for (c, problem) in column.problems.iter().enumerate() {
                    statistics[c].maximum_normalized_local_error = statistics[c]
                        .maximum_normalized_local_error
                        .max(cell_errors[c]);
                    rates[c] = problem
                        .rates(points[c], &mut statistics[c], &mut cancelled)
                        .map_err(|error| cell_error(c, error))?;
                }
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
                };
                step = (dt * factor).clamp(
                    controls.chemistry.minimum_step_s,
                    controls.chemistry.maximum_step_s,
                );
            }
            Ok(error) => {
                rejected += 1;
                step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
                if step < controls.chemistry.minimum_step_s {
                    return Err(Error::StepBudget);
                }
            }
            Err(
                error @ Error::Cell {
                    error: CollisionalRadiativeError::OutsideDomain,
                    ..
                },
            ) => {
                rejected += 1;
                step = 0.5 * dt;
                if step < controls.chemistry.minimum_step_s {
                    return Err(error);
                }
            }
            Err(error) => return Err(error),
        }
    }
    for statistics in &mut statistics {
        statistics.accepted_steps = accepted;
        statistics.rejected_steps = rejected;
    }
    column.publish(
        cells,
        &points,
        &accounts,
        statistics,
        area,
        duration,
        accepted,
        rejected,
        &mut cancelled,
    )
}

impl<'a> Column<'a> {
    fn new(
        cells: &[RadiationTransportCell<'a>],
        groups: &'a [PhotonNumberFluxV1],
        area: Area,
        duration: Time,
        controls: RadiationTransportSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        controls.chemistry.validate().map_err(Error::Chemistry)?;
        if !positive(area.value())
            || !nonnegative(duration.value())
            || !positive(controls.relative_photon_tolerance)
            || controls.relative_photon_tolerance >= 1.0
            || !(1..=MAX_PHOTON_GROUPS_V1).contains(&controls.maximum_groups)
        {
            return Err(Error::Invalid("geometry or photon controls"));
        }
        if cells.len() > controls.maximum_cells {
            return Err(Error::CellLimit {
                requested: cells.len(),
                maximum: controls.maximum_cells,
            });
        }
        if groups.len() > controls.maximum_groups {
            return Err(Error::GroupLimit {
                requested: groups.len(),
                maximum: controls.maximum_groups,
            });
        }
        let pairs = cells
            .len()
            .checked_mul(groups.len())
            .filter(|n| *n <= controls.maximum_cell_groups)
            .ok_or(Error::CellGroupLimit {
                cells: cells.len(),
                groups: groups.len(),
                maximum: controls.maximum_cell_groups,
            })?;
        let mut boundary = Vec::new();
        boundary
            .try_reserve_exact(groups.len())
            .map_err(|_| Error::Arithmetic("boundary allocation"))?;
        for (group, beam) in groups.iter().enumerate() {
            poll(cancelled)?;
            if !positive(beam.photon_energy.value()) || !nonnegative(beam.number_flux_per_m2_s) {
                return Err(Error::Band {
                    cell: None,
                    group,
                    error: PhotoionizationError::Invalid("boundary photon energy or flux"),
                });
            }
            let rate = beam.number_flux_per_m2_s * area.value();
            let photons = (beam.number_flux_per_m2_s * duration.value()) * area.value();
            if !nonnegative(rate)
                || !nonnegative(photons)
                || (beam.number_flux_per_m2_s > 0.0 && rate == 0.0)
            {
                return Err(Error::Arithmetic("boundary photon exposure"));
            }
            boundary.push([rate, photons]);
        }
        let mut problems = Vec::new();
        problems
            .try_reserve_exact(cells.len())
            .map_err(|_| Error::Arithmetic("cell preparation allocation"))?;
        for (c, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
            // Radiation is supplied by conservative column subflows. The local
            // reaction kernel receives no independent, replenished photon field.
            problems.push(
                Problem::new_with_thermal(
                    cell.initial,
                    cell.source,
                    &[],
                    cell.thermal,
                    controls.chemistry,
                    cancelled,
                )
                .map_err(|error| cell_error(c, error))?,
            );
        }
        let mut bands = Vec::new();
        bands
            .try_reserve_exact(pairs)
            .map_err(|_| Error::Arithmetic("column band allocation"))?;
        // Band-major layout keeps each incident-to-exit sweep contiguous.
        for (g, beam) in groups.iter().enumerate() {
            for (c, cell) in cells.iter().enumerate() {
                poll(cancelled)?;
                let source = cell.source.source();
                let sigma = source
                    .cross_section(beam.photon_energy)
                    .map_err(|error| Error::Band {
                        cell: Some(c),
                        group: g,
                        error: PhotoionizationError::Source(error),
                    })?
                    .value();
                let band = Band {
                    opacity_per_mol: (sigma * NA) / area.value(),
                    excess_j_mol: if sigma == 0.0 {
                        0.0
                    } else {
                        (beam.photon_energy.value()
                            - source.channel().ionization_energy_j_mol() / NA)
                            * NA
                    },
                };
                if !nonnegative(band.opacity_per_mol)
                    || !nonnegative(band.excess_j_mol)
                    || (sigma > 0.0 && band.opacity_per_mol == 0.0)
                {
                    return Err(Error::Arithmetic("prepared photon band"));
                }
                bands.push(band);
            }
        }
        Ok(Self {
            problems,
            groups,
            bands,
            boundary,
            controls,
        })
    }

    fn inactive(&self, points: &[Point], rates: &[Rates]) -> bool {
        rates
            .iter()
            .zip(points)
            .all(|(rate, point)| rate.inactive(*point))
            && self.boundary.iter().enumerate().all(|(g, input)| {
                input[0] == 0.0
                    || points.iter().enumerate().all(|(c, point)| {
                        point.neutral_amount() == 0.0
                            || self.bands[g * points.len() + c].opacity_per_mol == 0.0
                    })
            })
    }

    fn split_step(
        &self,
        points: &mut [Point],
        initial_rates: Option<&[Rates]>,
        dt: f64,
        accounts: &mut [[Dd; 2]],
        statistics: &mut [CollisionalRadiativeBalance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        // D(dt/2) P(dt) D(dt/2). Each D is the same second-order positive
        // exponential midpoint kernel as the local chemistry/thermal owner.
        for (c, problem) in self.problems.iter().enumerate() {
            poll(cancelled)?;
            let rate = match initial_rates {
                Some(rates) => rates[c],
                None => problem
                    .rates(points[c], &mut statistics[c], &mut *cancelled)
                    .map_err(|error| cell_error(c, error))?,
            };
            points[c] = problem
                .midpoint(
                    points[c],
                    rate,
                    dt * 0.5,
                    &mut statistics[c],
                    &mut *cancelled,
                )
                .map_err(|error| cell_error(c, error))?;
        }
        for g in 0..self.groups.len() {
            let h = if g + 1 == self.groups.len() {
                dt
            } else {
                dt * 0.5
            };
            self.photon_step(points, g, h, accounts, cancelled)?;
        }
        for g in (0..self.groups.len().saturating_sub(1)).rev() {
            self.photon_step(points, g, dt * 0.5, accounts, cancelled)?;
        }
        for (c, problem) in self.problems.iter().enumerate() {
            poll(cancelled)?;
            // Artificial individual-band birth temperatures need not lie inside
            // the interval; the completed physical spectrum must be supported.
            problem
                .admit_temperature(points[c].electron_temperature())
                .map_err(|error| cell_error(c, error))?;
            let rate = problem
                .rates(points[c], &mut statistics[c], &mut *cancelled)
                .map_err(|error| cell_error(c, error))?;
            points[c] = problem
                .midpoint(
                    points[c],
                    rate,
                    dt * 0.5,
                    &mut statistics[c],
                    &mut *cancelled,
                )
                .map_err(|error| cell_error(c, error))?;
        }
        Ok(())
    }

    fn photon_step(
        &self,
        points: &mut [Point],
        group: usize,
        dt: f64,
        accounts: &mut [[Dd; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let count = points.len();
        let mut photons = self.boundary[group][0] * dt;
        if !nonnegative(photons) {
            return Err(Error::Arithmetic("photon substep exposure"));
        }
        for (c, point) in points.iter_mut().enumerate() {
            poll(cancelled)?;
            let k = group * count + c;
            let band = self.bands[k];
            let neutral = point.neutral_amount();
            if photons == 0.0 || neutral == 0.0 || band.opacity_per_mol == 0.0 {
                accounts[k][1] = accounts[k][1] + Dd::from_f64(photons);
                continue;
            }
            let tau = band.opacity_per_mol * neutral;
            let dose = band.opacity_per_mol * (photons / NA);
            if !positive(tau) || !positive(dose) {
                return Err(Error::Arithmetic("positive band depth or dose"));
            }
            let (lost_fraction, remaining_fraction, _) =
                depletion(tau, dose).map_err(|error| Error::Band {
                    cell: Some(c),
                    group,
                    error,
                })?;
            let lost = neutral * lost_fraction;
            let remaining = neutral * remaining_fraction;
            let outgoing = transmitted_photons(photons, tau, dose);
            let absorbed = lost * NA;
            if !positive(lost)
                || !nonnegative(remaining)
                || !nonnegative(outgoing)
                || !nonnegative(absorbed)
            {
                return Err(Error::Arithmetic("photon depletion"));
            }
            *point = self.problems[c]
                .absorb_photons(
                    *point,
                    remaining,
                    lost,
                    lost * band.excess_j_mol,
                    absorbed * self.groups[group].photon_energy.value(),
                )
                .map_err(|error| cell_error(c, error))?;
            accounts[k][0] = accounts[k][0] + Dd::from_f64(absorbed);
            accounts[k][1] = accounts[k][1] + Dd::from_f64(outgoing);
            photons = outgoing;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn publish(
        &self,
        cells: &[RadiationTransportCell<'_>],
        points: &[Point],
        accounts: &[[Dd; 2]],
        mut statistics: Vec<CollisionalRadiativeBalance>,
        area: Area,
        duration: Time,
        accepted: usize,
        rejected: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RadiationTransportResult, Error> {
        let mut endpoints = Vec::new();
        endpoints
            .try_reserve_exact(cells.len())
            .map_err(|_| Error::Arithmetic("cell result allocation"))?;
        for (c, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
            let mut group_balances = Vec::new();
            group_balances
                .try_reserve_exact(self.groups.len())
                .map_err(|_| Error::Arithmetic("cell band result allocation"))?;
            let mut absorbed = [Dd::ZERO; 2];
            for g in 0..self.groups.len() {
                poll(cancelled)?;
                let k = g * cells.len() + c;
                let incoming = if c == 0 {
                    Dd::from_f64(self.boundary[g][1])
                } else {
                    accounts[k - 1][1]
                };
                let b = self.band_balance(Some(c), g, incoming, accounts[k][0], accounts[k][1])?;
                absorbed[0] = absorbed[0] + accounts[k][0];
                absorbed[1] = absorbed[1] + Dd::from_f64(b.absorbed_energy_j);
                group_balances.push(b);
            }
            let state = self.problems[c]
                .publish(points[c], &mut statistics[c], &mut *cancelled)
                .map_err(|error| cell_error(c, error))?;
            let b = statistics[c];
            let amount_residual = absorbed[0].to_f64() / NA - b.photoionization_extent_mol;
            let amount_allowance = self.controls.chemistry.absolute_amount_tolerance_mol
                + 64.0 * f64::EPSILON * b.photoionization_extent_mol;
            if !amount_residual.is_finite() || amount_residual.abs() > amount_allowance {
                return Err(Error::Arithmetic(
                    "cell absorbed photons and photoionization extent disagree",
                ));
            }
            let energy_residual =
                (absorbed[1] - Dd::from_f64(b.absorbed_radiation_energy_j)).to_f64();
            if !energy_residual.is_finite()
                || energy_residual.abs() > self.controls.chemistry.energy_tolerance_j
            {
                return Err(Error::EnergyBalance {
                    residual_j: energy_residual,
                });
            }
            endpoints.push(RadiationTransportCellResult {
                initial: cell.initial.clone(),
                state,
                source: cell.source.clone(),
                thermal: cell.thermal.cloned(),
                extent: Amount::new(points[c].reaction_extent_mol()),
                group_balances,
                balance: b,
            });
        }
        let mut group_balances = Vec::new();
        group_balances
            .try_reserve_exact(self.groups.len())
            .map_err(|_| Error::Arithmetic("column band result allocation"))?;
        let mut photons = Dd::ZERO;
        let mut photon_energy = Dd::ZERO;
        for (g, boundary) in self.boundary.iter().enumerate() {
            let mut absorbed = Dd::ZERO;
            for c in 0..cells.len() {
                if c % 256 == 0 {
                    poll(cancelled)?;
                }
                absorbed = absorbed + accounts[g * cells.len() + c][0];
            }
            let incoming = Dd::from_f64(boundary[1]);
            let outgoing = if cells.is_empty() {
                incoming
            } else {
                accounts[(g + 1) * cells.len() - 1][1]
            };
            let b = self.band_balance(None, g, incoming, absorbed, outgoing)?;
            photons = photons + absorbed;
            photon_energy = photon_energy + Dd::from_f64(b.absorbed_energy_j);
            group_balances.push(b);
        }
        let mut totals = [Dd::ZERO; 17];
        let mut balance = CollisionalRadiativeBalance::default();
        for (c, b) in statistics.iter().enumerate() {
            if c % 256 == 0 {
                poll(cancelled)?;
            }
            for (sum, value) in totals.iter_mut().zip([
                b.photoionization_extent_mol,
                b.collisional_extent_mol,
                b.radiative_recombination_extent_mol,
                b.absorbed_photons,
                b.emitted_photons,
                b.absorbed_radiation_energy_j,
                b.emitted_radiation_energy_j,
                b.chemical_storage_change_j,
                b.electron_sensible_energy_change_j,
                b.electron_deposited_energy_j,
                b.heavy_deposited_energy_j,
                b.heavy_sensible_energy_change_j,
                b.electron_to_heavy_energy_j,
                b.electron_to_heavy_neutral_energy_j,
                b.electron_to_heavy_coulomb_energy_j,
                b.electron_sensible_energy_residual_j,
                b.heavy_sensible_energy_residual_j,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            balance.maximum_normalized_local_error = balance
                .maximum_normalized_local_error
                .max(b.maximum_normalized_local_error);
            balance.maximum_relative_quadrature_error = balance
                .maximum_relative_quadrature_error
                .max(b.maximum_relative_quadrature_error);
            balance.maximum_relative_tail_bound = balance
                .maximum_relative_tail_bound
                .max(b.maximum_relative_tail_bound);
        }
        // Close actual endpoint gas energy against independently accumulated
        // finite beam absorption and escaping recombination radiation. Elastic
        // exchange is internal and must not be counted as external heating.
        let residual = (totals[7] + totals[8] + totals[11] + totals[6]
            - photon_energy
            - totals[9]
            - totals[10])
            .to_f64();
        if totals.iter().any(|x| !x.to_f64().is_finite()) || !residual.is_finite() {
            return Err(Error::Arithmetic("whole-column reaction/energy account"));
        }
        if residual.abs() > self.controls.chemistry.energy_tolerance_j {
            return Err(Error::EnergyBalance {
                residual_j: residual,
            });
        }
        balance.photoionization_extent_mol = totals[0].to_f64();
        balance.collisional_extent_mol = totals[1].to_f64();
        balance.radiative_recombination_extent_mol = totals[2].to_f64();
        balance.absorbed_photons = photons.to_f64();
        balance.emitted_photons = totals[4].to_f64();
        balance.absorbed_radiation_energy_j = photon_energy.to_f64();
        balance.emitted_radiation_energy_j = totals[6].to_f64();
        balance.chemical_storage_change_j = totals[7].to_f64();
        balance.electron_sensible_energy_change_j = totals[8].to_f64();
        balance.electron_deposited_energy_j = totals[9].to_f64();
        balance.heavy_deposited_energy_j = totals[10].to_f64();
        balance.heavy_sensible_energy_change_j = totals[11].to_f64();
        balance.electron_to_heavy_energy_j = totals[12].to_f64();
        balance.electron_to_heavy_neutral_energy_j = totals[13].to_f64();
        balance.electron_to_heavy_coulomb_energy_j = totals[14].to_f64();
        balance.electron_sensible_energy_residual_j = totals[15].to_f64();
        balance.heavy_sensible_energy_residual_j = totals[16].to_f64();
        balance.energy_residual_j = residual;
        balance.accepted_steps = accepted;
        balance.rejected_steps = rejected;
        poll(cancelled)?;
        Ok(RadiationTransportResult {
            cells: endpoints,
            groups: self.groups.to_vec(),
            group_balances,
            balance,
            area,
            elapsed: duration,
            controls: self.controls,
        })
    }

    fn band_balance(
        &self,
        cell: Option<usize>,
        group: usize,
        incoming: Dd,
        absorbed: Dd,
        outgoing: Dd,
    ) -> Result<RadiationTransportBandBalance, Error> {
        let residual = (absorbed + outgoing - incoming).to_f64();
        if !residual.is_finite()
            || residual.abs() > self.controls.relative_photon_tolerance * incoming.to_f64()
        {
            return Err(Error::PhotonBalance {
                cell,
                group,
                residual_photons: residual,
            });
        }
        let energy = absorbed.to_f64() * self.groups[group].photon_energy.value();
        if [
            incoming.to_f64(),
            absorbed.to_f64(),
            outgoing.to_f64(),
            energy,
        ]
        .iter()
        .any(|x| !nonnegative(*x))
        {
            return Err(Error::Arithmetic("band photon or energy balance"));
        }
        Ok(RadiationTransportBandBalance {
            incident_photons: incoming.to_f64(),
            absorbed_photons: absorbed.to_f64(),
            transmitted_photons: outgoing.to_f64(),
            absorbed_energy_j: energy,
            photon_residual: residual,
        })
    }
}

fn buffer<T: Clone>(count: usize, value: T) -> Result<Vec<T>, Error> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| Error::Arithmetic("bounded transport allocation"))?;
    result.resize(count, value);
    Ok(result)
}
fn cell_error(index: usize, error: CollisionalRadiativeError) -> Error {
    match error {
        CollisionalRadiativeError::Cancelled => Error::Cancelled,
        error => Error::Cell { index, error },
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
