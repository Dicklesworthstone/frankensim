use super::super::{
    DiffuseRecombinationCell,
    solve::{RadiationState, Workspace},
};
use super::*;
use crate::chemistry::collisional_radiative::solve::{DiffuseIncrement, Point, Problem};
use fs_math::dd::Dd;
use fs_math::det::{exp, ln};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;

type Error = DiffuseRecombinationEvolutionError;

struct Column<'a> {
    problems: Vec<Problem<'a>>,
    field: Workspace<'a>,
    raw: Vec<RadiationState>,
    controls: DiffuseRecombinationEvolutionSolve,
    balance: DiffuseRecombinationEvolutionBalance,
}

/// Per-cell self/nonlocal photons, net mol, net electron/gas J; per-face
/// right photons/J and left photons/J. Trial increments are distinct from the
/// cumulative accepted account so a rejected step cannot leak into publication.
struct Accounts {
    cells: Vec<[Dd; 5]>,
    faces: Vec<[Dd; 4]>,
}
impl Accounts {
    fn new(n: usize) -> Result<Self, Error> {
        Ok(Self {
            cells: buffer(n, [Dd::ZERO; 5])?,
            faces: buffer(n + 1, [Dd::ZERO; 4])?,
        })
    }
    fn clear(&mut self) {
        self.cells.fill([Dd::ZERO; 5]);
        self.faces.fill([Dd::ZERO; 4]);
    }
    fn add(&mut self, other: &Self) -> Result<(), Error> {
        for (total, increment) in self.cells.iter_mut().zip(&other.cells) {
            for j in 0..5 {
                total[j] = finite_sum(total[j], increment[j])?;
            }
        }
        for (total, increment) in self.faces.iter_mut().zip(&other.faces) {
            for j in 0..4 {
                total[j] = finite_sum(total[j], increment[j])?;
            }
        }
        Ok(())
    }
    fn integrate(&mut self, field: &Workspace<'_>, dt: f64) -> Result<(), Error> {
        for (account, cell) in self.cells.iter_mut().zip(&field.cells) {
            for (sum, rate) in account.iter_mut().zip([
                cell.self_absorbed_photons_per_s,
                cell.nonlocal_absorbed_photons_per_s,
                cell.net_ionization_rate_mol_per_s,
                cell.net_electron_power_w,
                cell.net_gas_power_w,
            ]) {
                *sum = finite_sum(*sum, Dd::from_f64(product(rate, dt)?))?;
            }
        }
        for (account, face) in self.faces.iter_mut().zip(&field.interfaces) {
            for (sum, rate) in account.iter_mut().zip([
                face.right_photons_per_s,
                face.right_power_w,
                face.left_photons_per_s,
                face.left_power_w,
            ]) {
                *sum = finite_sum(*sum, Dd::from_f64(product(rate, dt)?))?;
            }
        }
        Ok(())
    }
    fn error(
        &self,
        full: &Self,
        fine: &Self,
        controls: CollisionalRadiativeSolve,
    ) -> Result<f64, Error> {
        let mut error = 0.0_f64;
        for (c, old) in self.cells.iter().enumerate() {
            for j in 0..5 {
                let units = if j < 2 { NA } else { 1.0 };
                let absolute = if j < 3 {
                    controls.absolute_amount_tolerance_mol
                } else {
                    controls.energy_tolerance_j
                };
                error = error.max(local_error(
                    old[j],
                    full.cells[c][j],
                    fine.cells[c][j],
                    units,
                    absolute,
                    controls.relative_tolerance,
                )?);
            }
        }
        for (f, old) in self.faces.iter().enumerate() {
            for j in 0..4 {
                let units = if j % 2 == 0 { NA } else { 1.0 };
                let absolute = if j % 2 == 0 {
                    controls.absolute_amount_tolerance_mol
                } else {
                    controls.energy_tolerance_j
                };
                error = error.max(local_error(
                    old[j],
                    full.faces[f][j],
                    fine.faces[f][j],
                    units,
                    absolute,
                    controls.relative_tolerance,
                )?);
            }
        }
        Ok(error)
    }
}

pub(super) fn advance(
    cells: &[DiffuseRecombinationEvolutionCell<'_>],
    area: Area,
    duration: Time,
    controls: DiffuseRecombinationEvolutionSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DiffuseRecombinationEvolutionResult, Error> {
    let mut column = Column::new(cells, area, duration, controls, &mut cancelled)?;
    let mut points = buffer(cells.len(), Point::default())?;
    let mut statistics = buffer(cells.len(), CollisionalRadiativeBalance::default())?;
    for (c, problem) in column.problems.iter().enumerate() {
        points[c] = problem.initial_point();
    }
    let mut full = buffer(cells.len(), Point::default())?;
    let mut fine = buffer(cells.len(), Point::default())?;
    let mut predictor = buffer(cells.len(), Point::default())?;
    let mut accounts = Accounts::new(cells.len())?;
    let mut full_accounts = Accounts::new(cells.len())?;
    let mut fine_accounts = Accounts::new(cells.len())?;
    let (mut time, mut step) = (0.0, controls.chemistry.initial_step_s);
    while time < duration.value() {
        poll(&mut cancelled)?;
        if column.inactive(&points, &mut statistics, &mut cancelled)? {
            break;
        }
        if column.balance.accepted_steps + column.balance.rejected_steps
            >= controls.chemistry.maximum_steps
        {
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
        full_accounts.clear();
        fine_accounts.clear();
        let trial = (|| {
            column.split_step(
                &mut full,
                dt,
                &mut predictor,
                &mut full_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            column.split_step(
                &mut fine,
                0.5 * dt,
                &mut predictor,
                &mut fine_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            column.split_step(
                &mut fine,
                0.5 * dt,
                &mut predictor,
                &mut fine_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            let mut error = accounts.error(&full_accounts, &fine_accounts, controls.chemistry)?;
            for (c, problem) in column.problems.iter().enumerate() {
                poll(&mut cancelled)?;
                error = error.max(
                    problem
                        .error(points[c], full[c], fine[c])
                        .map_err(|e| cell_error(c, e))?,
                );
            }
            Ok(error)
        })();
        match trial {
            Ok(error) if error <= 1.0 => {
                accounts.add(&fine_accounts)?;
                core::mem::swap(&mut points, &mut fine);
                time = end;
                column.balance.accepted_steps += 1;
                column.balance.maximum_normalized_local_error =
                    column.balance.maximum_normalized_local_error.max(error);
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
                column.balance.rejected_steps += 1;
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
                column.balance.rejected_steps += 1;
                step = 0.5 * dt;
                if step < controls.chemistry.minimum_step_s {
                    return Err(error);
                }
            }
            Err(error) => return Err(error),
        }
    }
    column.publish(
        cells,
        &points,
        &accounts,
        statistics,
        area,
        duration,
        &mut cancelled,
    )
}

impl<'a> Column<'a> {
    fn new(
        cells: &[DiffuseRecombinationEvolutionCell<'a>],
        area: Area,
        duration: Time,
        controls: DiffuseRecombinationEvolutionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        controls
            .chemistry
            .validate()
            .map_err(|_| Error::Invalid("chemistry integration controls"))?;
        if !area.value().is_finite()
            || area.value() <= 0.0
            || !nonnegative(duration.value())
            || controls.maximum_radiation_evaluations == 0
            || controls.maximum_radiation_work == 0
        {
            return Err(Error::Invalid(
                "geometry, duration or radiation work controls",
            ));
        }
        if controls.radiation.maximum_cells > 4096 || cells.len() > controls.radiation.maximum_cells
        {
            return Err(Error::Invalid("radiation cell limit"));
        }
        let mut problems = Vec::new();
        let mut frozen = Vec::new();
        problems
            .try_reserve_exact(cells.len())
            .map_err(|_| Error::Arithmetic("prepared cell allocation"))?;
        frozen
            .try_reserve_exact(cells.len())
            .map_err(|_| Error::Arithmetic("source binding allocation"))?;
        for (c, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
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
            frozen.push(DiffuseRecombinationCell {
                state: cell.initial,
                source: cell.source,
            });
        }
        let mut radiation_controls = controls.radiation;
        radiation_controls.maximum_work = radiation_controls
            .maximum_work
            .min(controls.maximum_radiation_work);
        let field = Workspace::new(&frozen, area, radiation_controls, cancelled)
            .map_err(radiation_error)?;
        let work = field.preparation_work();
        if work > controls.maximum_radiation_work {
            return Err(Error::RadiationBudget);
        }
        Ok(Self {
            problems,
            field,
            raw: buffer(cells.len(), RadiationState::default())?,
            controls,
            balance: DiffuseRecombinationEvolutionBalance {
                radiation_work: work,
                ..Default::default()
            },
        })
    }

    fn inactive(
        &self,
        points: &[Point],
        statistics: &mut [CollisionalRadiativeBalance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Error> {
        // Only exact absence of every emitting pair permits the dark shortcut.
        // Gross self-recycling is still physical traffic even when its net is zero.
        for (c, problem) in self.problems.iter().enumerate() {
            let amounts = problem.active_amounts(points[c]);
            if amounts[1] > 0.0 && amounts[2] > 0.0 {
                return Ok(false);
            }
            let rate = problem
                .collisional_rates(points[c], &mut statistics[c], &mut *cancelled)
                .map_err(|error| cell_error(c, error))?;
            if !rate.inactive(points[c]) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn radiation(
        &mut self,
        points: &[Point],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if self.balance.radiation_evaluations >= self.controls.maximum_radiation_evaluations {
            return Err(Error::RadiationBudget);
        }
        for (c, problem) in self.problems.iter().enumerate() {
            self.raw[c] = RadiationState {
                amounts: problem.active_amounts(points[c]),
                electron_temperature_k: points[c].electron_temperature(),
            };
        }
        let remaining = self
            .controls
            .maximum_radiation_work
            .saturating_sub(self.balance.radiation_work);
        self.balance.radiation_evaluations += 1;
        self.field
            .evaluate(&self.raw, remaining, cancelled)
            .map_err(radiation_error)?;
        let b = self.field.balance;
        self.balance.radiation_work = self
            .balance
            .radiation_work
            .checked_add(b.work)
            .filter(|work| *work <= self.controls.maximum_radiation_work)
            .ok_or(Error::RadiationBudget)?;
        self.balance.maximum_normalized_spectral_estimate = self
            .balance
            .maximum_normalized_spectral_estimate
            .max(b.maximum_normalized_spectral_estimate);
        self.balance.maximum_normalized_angular_estimate = self
            .balance
            .maximum_normalized_angular_estimate
            .max(b.maximum_normalized_angular_estimate);
        self.balance.maximum_normalized_tail_bound = self
            .balance
            .maximum_normalized_tail_bound
            .max(b.maximum_normalized_tail_bound);
        Ok(())
    }

    fn split_step(
        &mut self,
        points: &mut [Point],
        dt: f64,
        predictor: &mut [Point],
        accounts: &mut Accounts,
        statistics: &mut [CollisionalRadiativeBalance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.local_half(points, 0.5 * dt, statistics, cancelled)?;
        self.radiation(points, cancelled)?;
        for (c, point) in points.iter().enumerate() {
            poll(cancelled)?;
            predictor[c] = self.apply(c, *point, 0.5 * dt)?;
        }
        self.radiation(predictor, cancelled)?;
        for (c, point) in points.iter_mut().enumerate() {
            poll(cancelled)?;
            *point = self.apply(c, *point, dt)?;
        }
        accounts.integrate(&self.field, dt)?;
        self.local_half(points, 0.5 * dt, statistics, cancelled)
    }

    fn local_half(
        &self,
        points: &mut [Point],
        dt: f64,
        statistics: &mut [CollisionalRadiativeBalance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for (c, problem) in self.problems.iter().enumerate() {
            poll(cancelled)?;
            let rates = problem
                .collisional_rates(points[c], &mut statistics[c], &mut *cancelled)
                .map_err(|error| cell_error(c, error))?;
            points[c] = problem
                .collisional_midpoint(points[c], rates, dt, &mut statistics[c], &mut *cancelled)
                .map_err(|error| cell_error(c, error))?;
        }
        Ok(())
    }

    fn apply(&self, c: usize, point: Point, dt: f64) -> Result<Point, Error> {
        let b = self.field.cells[c];
        self.problems[c]
            .apply_diffuse(
                point,
                DiffuseIncrement {
                    net_extent_mol: product(b.net_ionization_rate_mol_per_s, dt)?,
                    net_electron_energy_j: product(b.net_electron_power_w, dt)?,
                    net_gas_energy_j: product(b.net_gas_power_w, dt)?,
                    absorbed_mol: product(
                        b.self_absorbed_photons_per_s / NA + b.nonlocal_absorbed_photons_per_s / NA,
                        dt,
                    )?,
                    captured_mol: product(b.emitted_photons_per_s / NA, dt)?,
                    absorbed_energy_j: product(b.absorbed_power_w, dt)?,
                    captured_energy_j: product(b.captured_electron_power_w, dt)?,
                },
            )
            .map_err(|error| cell_error(c, error))
    }

    #[allow(clippy::too_many_arguments)]
    fn publish(
        mut self,
        cells: &[DiffuseRecombinationEvolutionCell<'_>],
        points: &[Point],
        accounts: &Accounts,
        mut statistics: Vec<CollisionalRadiativeBalance>,
        area: Area,
        duration: Time,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DiffuseRecombinationEvolutionResult, Error> {
        let mut endpoints = Vec::new();
        endpoints
            .try_reserve_exact(cells.len())
            .map_err(|_| Error::Arithmetic("endpoint allocation"))?;
        let mut gas = Dd::ZERO;
        let mut work = Dd::ZERO;
        let mut radiative_extent = Dd::ZERO;
        let mut net_radiation = Dd::ZERO;
        let mut net_traffic = 0.0;
        for (c, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
            statistics[c].accepted_steps = self.balance.accepted_steps;
            statistics[c].rejected_steps = self.balance.rejected_steps;
            statistics[c].maximum_normalized_local_error =
                self.balance.maximum_normalized_local_error;
            let state = self.problems[c]
                .publish(points[c], &mut statistics[c], &mut *cancelled)
                .map_err(|error| cell_error(c, error))?;
            let a = accounts.cells[c];
            let radiation = DiffuseRecombinationIntegratedRates {
                self_absorbed_photons: a[0].to_f64(),
                nonlocal_absorbed_photons: a[1].to_f64(),
                net_ionization_extent_mol: a[2].to_f64(),
                net_electron_energy_j: a[3].to_f64(),
                net_gas_energy_j: a[4].to_f64(),
            };
            let net_residual = (Dd::from_f64(points[c].net_radiation_energy_j()) - a[4]).to_f64();
            if !net_residual.is_finite()
                || net_residual.abs() > self.controls.chemistry.energy_tolerance_j
            {
                return Err(Error::EnergyBalance {
                    residual_j: net_residual,
                });
            }
            let b = statistics[c];
            for energy in [
                b.chemical_storage_change_j,
                b.electron_sensible_energy_change_j,
                b.heavy_sensible_energy_change_j,
            ] {
                gas = finite_sum(gas, Dd::from_f64(energy))?;
            }
            work = finite_sum(work, Dd::from_f64(b.electron_deposited_energy_j))?;
            work = finite_sum(work, Dd::from_f64(b.heavy_deposited_energy_j))?;
            radiative_extent = finite_sum(radiative_extent, a[2])?;
            net_radiation = finite_sum(net_radiation, a[4])?;
            net_traffic += radiation.net_ionization_extent_mol.abs();
            endpoints.push(DiffuseRecombinationEvolutionCellResult {
                initial: cell.initial.clone(),
                state,
                source: cell.source.clone(),
                thermal: cell.thermal.cloned(),
                extent: Amount::new(points[c].reaction_extent_mol()),
                balance: b,
                radiation,
            });
        }
        let mut interfaces = Vec::new();
        interfaces
            .try_reserve_exact(accounts.faces.len())
            .map_err(|_| Error::Arithmetic("exposure allocation"))?;
        for face in &accounts.faces {
            poll(cancelled)?;
            let v = face.map(Dd::to_f64);
            if v.iter().any(|x| !nonnegative(*x)) {
                return Err(Error::Arithmetic("integrated interface exposure"));
            }
            interfaces.push(DiffuseRecombinationExposure {
                right_photons: v[0],
                right_energy_j: v[1],
                left_photons: v[2],
                left_energy_j: v[3],
            });
        }
        let n = cells.len();
        let escaped = accounts.faces[0][2] + accounts.faces[n][0];
        let escaped_energy = accounts.faces[0][3] + accounts.faces[n][1];
        let photon_residual = (radiative_extent * Dd::from_f64(NA) + escaped).to_f64();
        let photon_allowance = NA * self.controls.chemistry.absolute_amount_tolerance_mol
            + 128.0 * f64::EPSILON * (NA * net_traffic + escaped.to_f64());
        if !photon_residual.is_finite()
            || !photon_allowance.is_finite()
            || photon_residual.abs() > photon_allowance
        {
            return Err(Error::PhotonBalance {
                residual_photons: photon_residual,
            });
        }
        let energy_residual = (gas + escaped_energy - work).to_f64();
        let radiation_residual = (net_radiation + escaped_energy).to_f64();
        for residual in [energy_residual, radiation_residual] {
            if !residual.is_finite() || residual.abs() > self.controls.chemistry.energy_tolerance_j
            {
                return Err(Error::EnergyBalance {
                    residual_j: residual,
                });
            }
        }
        self.balance.escaped_photons = escaped.to_f64();
        self.balance.escaped_energy_j = escaped_energy.to_f64();
        self.balance.net_radiative_extent_mol = radiative_extent.to_f64();
        self.balance.net_radiation_energy_j = net_radiation.to_f64();
        self.balance.gas_energy_change_j = gas.to_f64();
        self.balance.external_work_j = work.to_f64();
        self.balance.photon_residual = photon_residual;
        self.balance.energy_residual_j = energy_residual;
        poll(cancelled)?;
        Ok(DiffuseRecombinationEvolutionResult {
            cells: endpoints,
            interfaces,
            balance: self.balance,
            area,
            elapsed: duration,
            controls: self.controls,
        })
    }
}

fn local_error(
    old: Dd,
    full: Dd,
    fine: Dd,
    units: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, Error> {
    let difference = (full - fine).to_f64().abs() / units;
    let scale = (old + full).to_f64().abs().max((old + fine).to_f64().abs()) / units;
    let error = difference / (3.0 * (absolute + relative * scale));
    if !nonnegative(error) {
        return Err(Error::Arithmetic("integrated radiation local error"));
    }
    Ok(error)
}
fn product(rate: f64, dt: f64) -> Result<f64, Error> {
    let value = rate * dt;
    if !value.is_finite() || (rate != 0.0 && dt != 0.0 && value == 0.0) {
        return Err(Error::Arithmetic("radiation increment"));
    }
    Ok(value)
}
fn finite_sum(a: Dd, b: Dd) -> Result<Dd, Error> {
    let sum = a + b;
    if !sum.to_f64().is_finite() {
        return Err(Error::Arithmetic("cumulative radiation account"));
    }
    Ok(sum)
}
fn buffer<T: Clone>(n: usize, value: T) -> Result<Vec<T>, Error> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| Error::Arithmetic("bounded workspace allocation"))?;
    values.resize(n, value);
    Ok(values)
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn cell_error(index: usize, error: CollisionalRadiativeError) -> Error {
    if matches!(error, CollisionalRadiativeError::Cancelled) {
        Error::Cancelled
    } else {
        Error::Cell { index, error }
    }
}
fn radiation_error(error: DiffuseRecombinationError) -> Error {
    if matches!(error, DiffuseRecombinationError::Cancelled) {
        Error::Cancelled
    } else {
        Error::Radiation(error)
    }
}
