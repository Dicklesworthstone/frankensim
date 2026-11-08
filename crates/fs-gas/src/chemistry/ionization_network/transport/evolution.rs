use super::super::{IonizationChannelBalance, IonizationNetworkBalance};
use super::diffuse::{
    NetworkDiffuseCellExposure, NetworkDiffuseEvolutionBalance, NetworkDiffuseEvolutionSummary,
};
use super::*;
use crate::chemistry::collisional_radiative::solve::DiffuseIncrement;
use crate::chemistry::diffuse_recombination::DiffuseRecombinationSolve;
use crate::chemistry::diffuse_recombination::evolution::DiffuseRecombinationExposure;
use crate::chemistry::radiation_transport::RadiationTransportBandBalance;
use fs_math::dd::Dd;
use fs_math::det::{exp, ln};
use fs_qty::Time;
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;

mod diffuse_accounts;
mod radiation;
use radiation::RadiationField;

/// Original mixture, accepted endpoint and independent channel/beam accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkRadiationCellResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    sources: Vec<AtomicRadiativeRecombinationV1>,
    thermal: Option<CollisionalRadiativeThermal>,
    channels: Vec<IonizationChannelBalance>,
    beams: Vec<RadiationTransportBandBalance>,
    balance: IonizationNetworkBalance,
    diffuse: Option<NetworkDiffuseCellExposure>,
    free_free: Option<ElectronIonBremsstrahlungV1>,
}
impl NetworkRadiationCellResult {
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    pub fn sources(&self) -> &[AtomicRadiativeRecombinationV1] {
        &self.sources
    }
    pub const fn thermal_inputs(&self) -> Option<&CollisionalRadiativeThermal> {
        self.thermal.as_ref()
    }
    pub fn channel_balances(&self) -> &[IonizationChannelBalance] {
        &self.channels
    }
    pub fn beam_balances(&self) -> &[RadiationTransportBandBalance] {
        &self.beams
    }
    pub const fn balance(&self) -> IonizationNetworkBalance {
        self.balance
    }
    pub const fn diffuse_radiation(&self) -> Option<&NetworkDiffuseCellExposure> {
        self.diffuse.as_ref()
    }
    pub const fn free_free_cooling_source(&self) -> Option<&ElectronIonBremsstrahlungV1> {
        self.free_free.as_ref()
    }
}

/// Whole-column physical balances and all accepted/rejected numerical work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkRadiationBalance {
    pub absorbed_photons: f64,
    pub emitted_photons: f64,
    pub absorbed_radiation_energy_j: f64,
    /// Gross spontaneous emission; the optional diffuse summary records escape.
    pub emitted_radiation_energy_j: f64,
    /// Direct signed bound-free source, before gross-traffic subtraction.
    pub net_radiation_energy_j: f64,
    /// Escaping free-free energy under the separately selected thin closure.
    pub free_free_escaped_energy_j: f64,
    /// Direct total gas source: bound-free deposition minus free-free escape.
    /// This is retained before the two published component totals are rounded.
    pub net_gas_radiation_energy_j: f64,
    pub gas_energy_change_j: f64,
    pub external_work_j: f64,
    /// Actual gas storage minus the direct radiation source and held work.
    /// The original absolute chemistry tolerance applies once to the column.
    pub energy_residual_j: f64,
    /// Direct source plus emitted minus independently counted absorbed energy.
    pub photon_energy_residual_j: f64,
    /// Additional arithmetic allowance only for gross radiation-energy accounts.
    pub radiation_energy_roundoff_j: f64,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub maximum_normalized_local_error: f64,
    pub radiation_evaluations: usize,
    pub radiation_work: usize,
    /// Original rate queries across every cell and rejected trial, under the
    /// single `network.maximum_rate_evaluations` budget.
    pub rate_evaluations: usize,
    pub thermal_evaluations: usize,
    /// Included in rate_evaluations, across preparation and all trial work.
    pub free_free_evaluations: usize,
}

/// Complete source-retained evolution; refusal publishes no partial endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkRadiationResult {
    cells: Vec<NetworkRadiationCellResult>,
    beams: Vec<DiffuseRecombinationBeam>,
    beam_balances: Vec<RadiationTransportBandBalance>,
    balance: NetworkRadiationBalance,
    area: Area,
    elapsed: Time,
    controls: NetworkRadiationSolve,
    diffuse: Option<NetworkDiffuseEvolutionSummary>,
}
impl NetworkRadiationResult {
    pub fn cells(&self) -> &[NetworkRadiationCellResult] {
        &self.cells
    }
    pub fn beams(&self) -> &[DiffuseRecombinationBeam] {
        &self.beams
    }
    pub fn beam_balances(&self) -> &[RadiationTransportBandBalance] {
        &self.beam_balances
    }
    pub const fn balance(&self) -> NetworkRadiationBalance {
        self.balance
    }
    pub const fn area(&self) -> Area {
        self.area
    }
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    pub const fn controls(&self) -> NetworkRadiationSolve {
        self.controls
    }
    pub const fn diffuse_radiation(&self) -> Option<&NetworkDiffuseEvolutionSummary> {
        self.diffuse.as_ref()
    }
}

/// Evolve competing atomic ladders and shared thermal pools under finite beams.
///
/// Each radiation stage evaluates the original cross sections against current
/// populations. A common explicit-midpoint field updates every channel at once;
/// newly formed charge stages may absorb in subsequent stages of the same
/// interval. Symmetric composition surrounds it with the existing native dark
/// network kernel, including impact ionization, three-body and escaping
/// radiative recombination, optional elastic exchange and held pool powers.
/// Photons are never passed again as a maintained local intensity.
///
/// One full/two-half comparison controls the combined local second-order
/// estimate across populations, temperatures, reaction/energy ledgers and beam
/// exposures. Negative-population or out-of-domain trials reduce the step;
/// no clipping, population floor or artificial electron seed is introduced.
/// Source tables and work buffers persist across every trial. This is not an
/// unconditionally stiff method or a global-error certificate. Recombination
/// emission escapes in this entry. The diffuse network entry reuses the same
/// controller with spectral transport of spontaneous radiation.
pub fn advance_network_radiation_transport(
    cells: &[NetworkRadiationCell<'_>],
    area: Area,
    beams: &[DiffuseRecombinationBeam],
    duration: Time,
    controls: NetworkRadiationSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<NetworkRadiationResult, NetworkRadiationError> {
    advance_impl(cells, area, beams, duration, controls, None, cancelled)
}

pub(super) fn advance_impl(
    cells: &[NetworkRadiationCell<'_>],
    area: Area,
    beams: &[DiffuseRecombinationBeam],
    duration: Time,
    controls: NetworkRadiationSolve,
    diffuse: Option<DiffuseRecombinationSolve>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<NetworkRadiationResult, NetworkRadiationError> {
    controls.validate(cells, beams.len())?;
    poll(&mut cancelled)?;
    if !nonnegative(duration.value()) {
        return Err(Error::Invalid("photon exposure duration"));
    }
    let problems = prepare_problems(cells, controls, &mut cancelled)?;
    let field = RadiationField::new(
        cells,
        &problems,
        area,
        beams,
        controls,
        diffuse,
        &mut cancelled,
    )?;
    let mut column = Column {
        problems,
        field,
        controls,
        balance: NetworkRadiationBalance::default(),
    };
    let mut points: Vec<_> = column.problems.iter().map(Problem::initial_point).collect();
    let mut full = trial_points(&points)?;
    let mut fine = trial_points(&points)?;
    let mut predictor = trial_points(&points)?;
    let mut statistics = buffer(cells.len(), IonizationNetworkBalance::default())?;
    for (statistics, problem) in statistics.iter_mut().zip(&column.problems) {
        let preparation = problem.cooling_preparation_evaluations();
        statistics.rate_evaluations = preparation;
        statistics.free_free_evaluations = preparation;
        column.balance.rate_evaluations += preparation;
        column.balance.free_free_evaluations += preparation;
    }
    let mut increments = reserved(cells.len())?;
    for point in &points {
        increments.push(buffer(point.events.len(), DiffuseIncrement::default())?);
    }
    let mut accounts = Accounts::new(cells.len(), beams.len(), diffuse.is_some())?;
    let mut full_accounts = Accounts::new(cells.len(), beams.len(), diffuse.is_some())?;
    let mut fine_accounts = Accounts::new(cells.len(), beams.len(), diffuse.is_some())?;
    let local = controls.network.local;
    let mut time = 0.0;
    let mut step = local.initial_step_s;
    if cells.is_empty() && duration.value() > 0.0 {
        column.field.evaluate(&points, &mut cancelled)?;
        accounts.integrate(&column.field, duration.value())?;
        time = duration.value();
    }
    while time < duration.value() {
        poll(&mut cancelled)?;
        if !column.field.illuminated()
            && column
                .problems
                .iter()
                .zip(&points)
                .all(|(p, x)| p.inactive(x))
        {
            break;
        }
        if column.balance.accepted_steps + column.balance.rejected_steps >= local.maximum_steps {
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
        if end <= time || 0.5 * dt == 0.0 {
            return Err(Error::StepBudget);
        }
        copy_points(&mut full, &points);
        copy_points(&mut fine, &points);
        full_accounts.clear();
        fine_accounts.clear();
        let trial = (|| {
            column.split_step(
                &mut full,
                dt,
                &mut predictor,
                &mut increments,
                &mut full_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            column.split_step(
                &mut fine,
                0.5 * dt,
                &mut predictor,
                &mut increments,
                &mut fine_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            column.split_step(
                &mut fine,
                0.5 * dt,
                &mut predictor,
                &mut increments,
                &mut fine_accounts,
                &mut statistics,
                &mut cancelled,
            )?;
            let mut error = accounts.error(&full_accounts, &fine_accounts, beams, controls)?;
            for c in 0..cells.len() {
                poll(&mut cancelled)?;
                error = error.max(
                    column.problems[c]
                        .error(&points[c], &full[c], &fine[c])
                        .map_err(|error| cell_error(c, error))?,
                );
            }
            Ok(error)
        })();
        match trial {
            Ok(error) if error <= 1.0 => {
                core::mem::swap(&mut points, &mut fine);
                accounts.add(&fine_accounts)?;
                time = end;
                column.balance.accepted_steps += 1;
                column.balance.maximum_normalized_local_error =
                    column.balance.maximum_normalized_local_error.max(error);
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
                };
                step = (dt * factor)
                    .min(local.maximum_step_s)
                    .max(local.minimum_step_s);
            }
            Ok(error) => {
                column.balance.rejected_steps += 1;
                step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
                if step < local.minimum_step_s {
                    return Err(Error::StepBudget);
                }
            }
            Err(Error::Cell {
                error: IonizationNetworkError::OutsideDomain,
                ..
            }) => {
                column.balance.rejected_steps += 1;
                step = 0.5 * dt;
                if step < local.minimum_step_s {
                    return Err(Error::StepBudget);
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
        beams,
        duration,
        &mut cancelled,
    )
}

struct Column<'a> {
    problems: Vec<Problem<'a>>,
    field: RadiationField<'a>,
    controls: NetworkRadiationSolve,
    balance: NetworkRadiationBalance,
}
impl Column<'_> {
    fn dark_step(
        &mut self,
        points: &mut [Point],
        dt: f64,
        statistics: &mut [IonizationNetworkBalance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for c in 0..points.len() {
            poll(cancelled)?;
            if self.problems[c].inactive(&points[c]) {
                continue;
            }
            // Let the original kernel enforce one global rate budget, then
            // retain this cell's actual work count separately for publication.
            let before_total = self.balance.rate_evaluations;
            let before_cell = statistics[c].rate_evaluations;
            let thermal_before = statistics[c].thermal_evaluations;
            let free_free_before = statistics[c].free_free_evaluations;
            statistics[c].rate_evaluations = before_total;
            let result = if self.field.is_diffuse() {
                self.problems[c].collisional_compose(
                    &mut points[c],
                    dt,
                    &mut statistics[c],
                    cancelled,
                )
            } else {
                self.problems[c].compose(&mut points[c], dt, &mut statistics[c], cancelled)
            };
            let after_total = statistics[c].rate_evaluations;
            self.balance.rate_evaluations = after_total;
            statistics[c].rate_evaluations = before_cell + (after_total - before_total);
            self.balance.thermal_evaluations += statistics[c].thermal_evaluations - thermal_before;
            self.balance.free_free_evaluations +=
                statistics[c].free_free_evaluations - free_free_before;
            result.map_err(|error| cell_error(c, error))?;
        }
        Ok(())
    }

    fn split_step(
        &mut self,
        points: &mut [Point],
        dt: f64,
        predictor: &mut [Point],
        increments: &mut [Vec<DiffuseIncrement>],
        accounts: &mut Accounts,
        statistics: &mut [IonizationNetworkBalance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.dark_step(points, 0.5 * dt, statistics, cancelled)?;
        if self.field.radiation_enabled() {
            self.field.evaluate(points, cancelled)?;
            copy_points(predictor, points);
            self.absorb(predictor, 0.5 * dt, increments, cancelled)?;
            self.field.evaluate(predictor, cancelled)?;
            self.absorb(points, dt, increments, cancelled)?;
            accounts.integrate(&self.field, dt)?;
        }
        self.dark_step(points, 0.5 * dt, statistics, cancelled)
    }

    fn absorb(
        &self,
        points: &mut [Point],
        dt: f64,
        increments: &mut [Vec<DiffuseIncrement>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for c in 0..points.len() {
            poll(cancelled)?;
            let totals = self.field.increments(c, dt, &mut increments[c])?;
            let result = if let Some(totals) = totals {
                self.problems[c].apply_diffuse_with_totals(&mut points[c], &increments[c], totals)
            } else {
                self.problems[c].apply_diffuse(&mut points[c], &increments[c])
            };
            result.map_err(|error| cell_error(c, error))?;
        }
        Ok(())
    }

    fn publish(
        mut self,
        cells: &[NetworkRadiationCell<'_>],
        points: &[Point],
        accounts: &Accounts,
        mut statistics: Vec<IonizationNetworkBalance>,
        area: Area,
        beams: &[DiffuseRecombinationBeam],
        duration: Time,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkRadiationResult, Error> {
        let mut results = reserved(cells.len())?;
        // Gas storage, bound-free source, recombination, work, photons,
        // free-free escape and the direct combined radiation source.
        let mut totals = [Dd::ZERO; 7];
        for (c, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
            statistics[c].accepted_steps = self.balance.accepted_steps;
            statistics[c].rejected_steps = self.balance.rejected_steps;
            statistics[c].maximum_normalized_local_error =
                self.balance.maximum_normalized_local_error;
            let (state, channels) = self.problems[c]
                .publish(&points[c], duration.value(), &mut statistics[c], cancelled)
                .map_err(|error| cell_error(c, error))?;
            let b = statistics[c];
            let mut band_results = reserved(beams.len())?;
            let mut absorbed = Dd::ZERO;
            let mut absorbed_energy = Dd::ZERO;
            for (g, beam) in beams.iter().enumerate() {
                poll(cancelled)?;
                let a = accounts.cells[c * beams.len() + g];
                let result = band_balance(Some(c), g, a, *beam, self.controls)?;
                absorbed = sum(absorbed, a[1])?;
                absorbed_energy = sum(absorbed_energy, Dd::from_f64(result.absorbed_energy_j))?;
                band_results.push(result);
            }
            let diffuse = accounts.diffuse_cells.get(c).copied();
            if let Some(a) = diffuse {
                absorbed = sum(absorbed, a[1] + a[2] + a[3])?;
                absorbed_energy = sum(absorbed_energy, a[4])?;
            }
            let photon_residual = (absorbed - Dd::from_f64(b.absorbed_photons)).to_f64();
            let photon_allowance = NA * self.controls.network.local.absolute_amount_tolerance_mol
                + 128.0 * f64::EPSILON * absorbed.to_f64().max(b.absorbed_photons);
            if !photon_residual.is_finite() || photon_residual.abs() > photon_allowance {
                return Err(Error::Arithmetic(
                    "channel and transported absorption disagree",
                ));
            }
            let energy_residual =
                (absorbed_energy - Dd::from_f64(b.absorbed_radiation_energy_j)).to_f64();
            energy_check(
                energy_residual,
                self.controls.network.local.energy_tolerance_j
                    + 128.0
                        * f64::EPSILON
                        * absorbed_energy.to_f64().max(b.absorbed_radiation_energy_j),
            )?;
            for value in [
                b.chemical_storage_change_j,
                b.electron_sensible_energy_change_j,
                b.heavy_sensible_energy_change_j,
            ] {
                totals[0] = sum(totals[0], Dd::from_f64(value))?;
            }
            totals[1] = sum(totals[1], Dd::from_f64(b.net_radiation_energy_j))?;
            totals[2] = sum(totals[2], Dd::from_f64(b.emitted_radiation_energy_j))?;
            for value in [b.electron_deposited_energy_j, b.heavy_deposited_energy_j] {
                totals[3] = sum(totals[3], Dd::from_f64(value))?;
            }
            totals[4] = sum(totals[4], Dd::from_f64(b.emitted_photons))?;
            totals[5] = sum(totals[5], Dd::from_f64(b.free_free_escaped_energy_j))?;
            totals[6] = sum(totals[6], Dd::from_f64(b.net_gas_radiation_energy_j))?;
            results.push(NetworkRadiationCellResult {
                initial: cell.initial.clone(),
                state,
                sources: self.problems[c].sources().cloned().collect(),
                thermal: cell.thermal.cloned(),
                channels,
                beams: band_results,
                balance: b,
                diffuse: diffuse.map(radiation::exposure),
                free_free: cell.optically_thin_bremsstrahlung.cloned(),
            });
        }
        let mut band_results = reserved(beams.len())?;
        let mut absorbed = Dd::ZERO;
        let mut absorbed_energy = Dd::ZERO;
        for (g, beam) in beams.iter().enumerate() {
            poll(cancelled)?;
            let a = accounts.column[g];
            let b = band_balance(None, g, a, *beam, self.controls)?;
            absorbed = sum(absorbed, a[1])?;
            absorbed_energy = sum(absorbed_energy, Dd::from_f64(b.absorbed_energy_j))?;
            band_results.push(b);
        }
        for a in &accounts.diffuse_cells {
            absorbed = sum(absorbed, a[1] + a[2] + a[3])?;
            absorbed_energy = sum(absorbed_energy, a[4])?;
        }
        let [
            gas,
            radiation,
            emitted,
            work,
            emitted_photons,
            free_free,
            net_gas_radiation,
        ] = totals;
        let residual = (gas - net_gas_radiation - work).to_f64();
        energy_check(residual, self.controls.network.local.energy_tolerance_j)?;
        let photon_residual = (radiation + emitted - absorbed_energy).to_f64();
        let roundoff = 128.0
            * f64::EPSILON
            * emitted
                .to_f64()
                .max(absorbed_energy.to_f64())
                .max(radiation.to_f64().abs());
        energy_check(
            photon_residual,
            self.controls.network.local.energy_tolerance_j + roundoff,
        )?;
        self.balance.absorbed_photons = absorbed.to_f64();
        self.balance.emitted_photons = emitted_photons.to_f64();
        self.balance.absorbed_radiation_energy_j = absorbed_energy.to_f64();
        self.balance.emitted_radiation_energy_j = emitted.to_f64();
        self.balance.net_radiation_energy_j = radiation.to_f64();
        self.balance.free_free_escaped_energy_j = free_free.to_f64();
        self.balance.net_gas_radiation_energy_j = net_gas_radiation.to_f64();
        self.balance.gas_energy_change_j = gas.to_f64();
        self.balance.external_work_j = work.to_f64();
        self.balance.energy_residual_j = residual;
        self.balance.photon_energy_residual_j = photon_residual;
        self.balance.radiation_energy_roundoff_j = roundoff;
        self.balance.radiation_evaluations = self.field.evaluations();
        self.balance.radiation_work = self.field.work();
        let diffuse = accounts.publish_diffuse(&self.field, self.balance, self.controls)?;
        poll(cancelled)?;
        Ok(NetworkRadiationResult {
            cells: results,
            beams: beams.to_vec(),
            beam_balances: band_results,
            balance: self.balance,
            area,
            elapsed: duration,
            controls: self.controls,
            diffuse,
        })
    }
}

/// Direct nonnegative incoming/absorbed/outgoing counts, never reconstructed
/// by subtracting a bright transmitted beam from its incident exposure.
struct Accounts {
    cells: Vec<[Dd; 3]>,
    column: Vec<[Dd; 3]>,
    diffuse_cells: Vec<[Dd; 10]>,
    faces: Vec<[Dd; 4]>,
}
impl Accounts {
    fn new(cells: usize, beams: usize, diffuse: bool) -> Result<Self, Error> {
        let faces = if diffuse {
            cells.checked_add(1).ok_or(Error::WorkBudget)?
        } else {
            0
        };
        Ok(Self {
            cells: buffer(cells * beams, [Dd::ZERO; 3])?,
            column: buffer(beams, [Dd::ZERO; 3])?,
            diffuse_cells: buffer(if diffuse { cells } else { 0 }, [Dd::ZERO; 10])?,
            faces: buffer(faces, [Dd::ZERO; 4])?,
        })
    }
    fn clear(&mut self) {
        self.cells.fill([Dd::ZERO; 3]);
        self.column.fill([Dd::ZERO; 3]);
        self.diffuse_cells.fill([Dd::ZERO; 10]);
        self.faces.fill([Dd::ZERO; 4]);
    }
    fn add(&mut self, increment: &Self) -> Result<(), Error> {
        for (target, source) in self
            .cells
            .iter_mut()
            .chain(&mut self.column)
            .zip(increment.cells.iter().chain(&increment.column))
        {
            for j in 0..3 {
                target[j] = sum(target[j], source[j])?;
            }
        }
        self.add_diffuse(increment)?;
        Ok(())
    }
    fn integrate(&mut self, field: &RadiationField<'_>, dt: f64) -> Result<(), Error> {
        let beams = self.column.len();
        for (k, target) in self.cells.iter_mut().enumerate() {
            let rate = field.cell_beam_rates(k / beams)[k % beams];
            integrate_rate(target, rate, dt)?;
        }
        for (target, rate) in self.column.iter_mut().zip(field.beam_rates()) {
            integrate_rate(target, *rate, dt)?;
        }
        self.integrate_diffuse(field, dt)?;
        Ok(())
    }
    fn error(
        &self,
        full: &Self,
        fine: &Self,
        beams: &[DiffuseRecombinationBeam],
        controls: NetworkRadiationSolve,
    ) -> Result<f64, Error> {
        let mut error = 0.0_f64;
        let local = controls.network.local;
        for (k, ((old, a), b)) in self
            .cells
            .iter()
            .chain(&self.column)
            .zip(full.cells.iter().chain(&full.column))
            .zip(fine.cells.iter().chain(&fine.column))
            .enumerate()
        {
            let energy = beams[k % beams.len()].photons.photon_energy.value();
            for j in 0..3 {
                let delta = (a[j] - b[j]).to_f64().abs();
                let scale = (old[j] + a[j])
                    .to_f64()
                    .abs()
                    .max((old[j] + b[j]).to_f64().abs());
                let amount_error = (delta / NA)
                    / (3.0
                        * (local.absolute_amount_tolerance_mol
                            + local.relative_tolerance * (scale / NA)));
                let energy_error = (delta * energy)
                    / (3.0
                        * (local.energy_tolerance_j + local.relative_tolerance * scale * energy));
                if !nonnegative(amount_error) || !nonnegative(energy_error) {
                    return Err(Error::Arithmetic("beam exposure local error"));
                }
                error = error.max(amount_error).max(energy_error);
            }
        }
        Ok(error.max(self.diffuse_error(full, fine, controls)?))
    }
}
fn integrate_rate(
    target: &mut [Dd; 3],
    rate: NetworkRadiationBeamRate,
    dt: f64,
) -> Result<(), Error> {
    for (sum_value, value) in target.iter_mut().zip([
        rate.incident_photons_per_s,
        rate.absorbed_photons_per_s,
        rate.transmitted_photons_per_s,
    ]) {
        *sum_value = sum(*sum_value, Dd::from_f64(product(value, dt)?))?;
    }
    Ok(())
}
fn band_balance(
    cell: Option<usize>,
    beam: usize,
    a: [Dd; 3],
    source: DiffuseRecombinationBeam,
    controls: NetworkRadiationSolve,
) -> Result<RadiationTransportBandBalance, Error> {
    let values = a.map(Dd::to_f64);
    if values.iter().any(|v| !nonnegative(*v)) {
        return Err(Error::Arithmetic("integrated beam account"));
    }
    let residual = (a[1] + a[2] - a[0]).to_f64();
    if !residual.is_finite() || residual.abs() > controls.relative_photon_tolerance * values[0] {
        return Err(Error::PhotonBalance {
            cell,
            beam,
            residual_photons: residual,
        });
    }
    let energy = product(values[1], source.photons.photon_energy.value())?;
    Ok(RadiationTransportBandBalance {
        incident_photons: values[0],
        absorbed_photons: values[1],
        transmitted_photons: values[2],
        absorbed_energy_j: energy,
        photon_residual: residual,
    })
}
fn trial_points(points: &[Point]) -> Result<Vec<Point>, Error> {
    let mut trial = reserved(points.len())?;
    for point in points {
        trial.push(Point::new(point.amounts.len(), point.events.len()));
    }
    Ok(trial)
}
fn copy_points(target: &mut [Point], source: &[Point]) {
    for (target, source) in target.iter_mut().zip(source) {
        target.copy_from(source);
    }
}
fn product(a: f64, b: f64) -> Result<f64, Error> {
    let value = a * b;
    if !nonnegative(value) || (a != 0.0 && b != 0.0 && value == 0.0) {
        return Err(Error::Arithmetic("photon exposure increment"));
    }
    Ok(value)
}
fn sum(a: Dd, b: Dd) -> Result<Dd, Error> {
    let value = a + b;
    if !value.to_f64().is_finite() {
        return Err(Error::Arithmetic("cumulative radiation account"));
    }
    Ok(value)
}
fn energy_check(residual: f64, allowance: f64) -> Result<(), Error> {
    if !residual.is_finite() || !nonnegative(allowance) || residual.abs() > allowance {
        Err(Error::EnergyBalance {
            residual_j: residual,
        })
    } else {
        Ok(())
    }
}
