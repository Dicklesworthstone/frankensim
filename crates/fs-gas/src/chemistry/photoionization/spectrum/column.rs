//! Photon-conserving multigroup transport through ordered, mixed gas cells.
//!
//! Each band's exact subflow advances the entire column: a cell receives the
//! integrated photons escaping the evolving upstream gas, not the photons
//! transmitted by its initial opacity. Symmetric composition of these column
//! subflows resolves simultaneous competition between bands. Full/two-half
//! comparison controls second-order local errors in every cell's population,
//! electron temperature and band photon accounts. No average photon energy or
//! opacity replaces the spectrum, and trial steps allocate no storage.
//!
//! The beam area is fixed. Propagation is quasi-static; each cell is uniformly
//! mixed, with its original atomic source and instantaneous photoelectron
//! thermalization. Cell endpoints resolve shielding and spectral hardening
//! along the column, not structure within cells or light-travel delays. Other
//! absorbers, scattering, collisional reactions, recombination and motion are
//! absent during this exposure. Error estimates are local numerical controls,
//! not empirical, spectral-quadrature or spatial-discretization certificates.

use super::*;
use crate::chemistry::photoionization::column::PhotoionizationCell;

/// Source, time-integration and bounded spatial work controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralPhotoionizationColumnSolve {
    /// Common source-temperature interval and numerical allowances. The energy
    /// allowance applies to each cell AND once to the complete column.
    pub cell: SpectralPhotoionizationSolve,
    /// Maximum cells admitted before result or work allocation.
    pub maximum_cells: usize,
    /// Maximum product of cells and groups, bounding all band work buffers.
    pub maximum_cell_groups: usize,
}
impl SpectralPhotoionizationColumnSolve {
    /// Explicit electron source interval with bounded numerical defaults.
    pub fn new(minimum: Temperature, maximum: Temperature) -> Self {
        Self {
            cell: SpectralPhotoionizationSolve::new(minimum, maximum),
            maximum_cells: 4096,
            maximum_cell_groups: 1_048_576,
        }
    }
}

/// Atomic refusal; caller-owned cells and sources are never modified.
#[derive(Clone, Debug, PartialEq)]
pub enum SpectralPhotoionizationColumnError {
    /// Invalid common exposure or numerical controls.
    Invalid(SpectralPhotoionizationError),
    /// Explicit spatial work budget exceeded before allocation.
    CellLimit { requested: usize, maximum: usize },
    /// The cell/group product exceeds its budget or cannot be represented.
    CellGroupLimit {
        cells: usize,
        groups: usize,
        maximum: usize,
    },
    /// Original source, state or numerical update refused in this cell.
    Cell {
        index: usize,
        error: SpectralPhotoionizationError,
    },
    /// Accepted plus rejected attempts, minimum step or time resolution exhausted.
    StepBudget,
    /// A required physical value or bounded allocation cannot be represented.
    Arithmetic(&'static str),
    /// Summed cell absorption plus escaped photons failed a band's input account.
    PhotonBalance { group: usize, residual_photons: f64 },
    /// Whole-column material energy failed absorbed radiation closure.
    EnergyBalance { residual_j: f64 },
    /// Cancellation before publication of a complete column.
    Cancelled,
}
impl core::fmt::Display for SpectralPhotoionizationColumnError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral photoionization column: {self:?}")
    }
}
impl std::error::Error for SpectralPhotoionizationColumnError {}

/// Original-source endpoint and integrated radiation accounts for one cell.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPhotoionizationColumnCellResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicPhotoionizationV1,
    extent: Amount,
    group_balances: Vec<PhotoionizationBalance>,
    balance: PhotoionizationBalance,
}
impl SpectralPhotoionizationColumnCellResult {
    /// Actual endpoint for existing chemistry, electrical and optical consumers.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Immutable original gas inventory and temperatures.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original atomic cross section, levels and thermodynamic bank.
    pub const fn source(&self) -> &AtomicPhotoionizationV1 {
        &self.source
    }
    /// One ion/electron pair per absorbed photon [mol].
    pub const fn reaction_extent(&self) -> Amount {
        self.extent
    }
    /// Integrated incoming/absorbed/outgoing counts in the declared band order.
    /// Downstream illumination varies with time; these are not constant fluxes.
    /// Each band's sensible contribution is its deposited photoelectron energy.
    pub fn group_balances(&self) -> &[PhotoionizationBalance] {
        &self.group_balances
    }
    /// Total photons and energy change measured from the actual gas endpoint.
    pub const fn balance(&self) -> PhotoionizationBalance {
        self.balance
    }
}

/// Complete stratified gas endpoints and frequency-resolved column accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPhotoionizationColumnResult {
    cells: Vec<SpectralPhotoionizationColumnCellResult>,
    groups: Vec<PhotonNumberFluxV1>,
    group_balances: Vec<PhotoionizationBalance>,
    balance: PhotoionizationBalance,
    area: Area,
    elapsed: Time,
    controls: SpectralPhotoionizationColumnSolve,
    accepted_steps: usize,
    rejected_steps: usize,
    maximum_normalized_local_error: f64,
}
impl SpectralPhotoionizationColumnResult {
    /// Endpoints ordered from the illuminated face toward the exit.
    pub fn cells(&self) -> &[SpectralPhotoionizationColumnCellResult] {
        &self.cells
    }
    /// Original constant boundary illumination in deterministic band order.
    pub fn groups(&self) -> &[PhotonNumberFluxV1] {
        &self.groups
    }
    /// Boundary input, summed cell absorption and column exit, separately by band.
    pub fn group_balances(&self) -> &[PhotoionizationBalance] {
        &self.group_balances
    }
    /// Summed photons and actual original-source material-energy change.
    pub const fn balance(&self) -> PhotoionizationBalance {
        self.balance
    }
    /// Common area normal to every band [m²].
    pub const fn illuminated_area(&self) -> Area {
        self.area
    }
    /// Completed exposure [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Applied source, accuracy and work controls.
    pub const fn controls(&self) -> SpectralPhotoionizationColumnSolve {
        self.controls
    }
    /// Accepted complete steps, each using two symmetric half steps.
    pub const fn accepted_steps(&self) -> usize {
        self.accepted_steps
    }
    /// Rejected attempts charged to the same step budget.
    pub const fn rejected_steps(&self) -> usize {
        self.rejected_steps
    }
    /// Largest accepted local full/two-half difference divided by its allowance.
    pub const fn maximum_normalized_local_error(&self) -> f64 {
        self.maximum_normalized_local_error
    }
}

#[derive(Clone, Copy)]
struct Band {
    opacity_per_mol: f64,
    excess_j_mol: f64,
}

/// Advance a constant discrete spectrum through spatially ordered gas cells.
///
/// Cross sections prepare once per cell and band. Every positive column
/// subflow conserves photons and avoids a reaction-stiffness positivity limit;
/// temporal refinement controls the remaining competition between bands.
/// The one-band limit is the existing exact serial fluence solution. Empty
/// columns transmit all photons; dark/empty spectra preserve every gas cell.
/// Source bindings and caloric support are checked even when no photons arrive.
pub fn advance_spectral_column(
    cells: &[PhotoionizationCell<'_>],
    groups: &[PhotonNumberFluxV1],
    area: Area,
    duration: Time,
    controls: SpectralPhotoionizationColumnSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpectralPhotoionizationColumnResult, SpectralPhotoionizationColumnError> {
    use SpectralPhotoionizationColumnError as Error;
    column_poll(&mut cancelled)?;
    controls.cell.validate().map_err(Error::Invalid)?;
    if cells.len() > controls.maximum_cells {
        return Err(Error::CellLimit {
            requested: cells.len(),
            maximum: controls.maximum_cells,
        });
    }
    if groups.len() > controls.cell.maximum_groups {
        return Err(Error::Invalid(SpectralPhotoionizationError::GroupLimit {
            requested: groups.len(),
            maximum: controls.cell.maximum_groups,
        }));
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
    if !positive(area.value()) || !nonnegative(duration.value()) {
        return Err(Error::Invalid(SpectralPhotoionizationError::Invalid(
            "column area or duration",
        )));
    }
    let mut boundary = Vec::new();
    boundary
        .try_reserve_exact(groups.len())
        .map_err(|_| Error::Arithmetic("boundary allocation"))?;
    for (index, &beam) in groups.iter().enumerate() {
        column_poll(&mut cancelled)?;
        let (_, photons) = validate_exposure(
            PhotoionizationExposure {
                beam,
                illuminated_area: area,
                duration,
            },
            PhotoionizationControls {
                relative_charge_tolerance: controls.cell.relative_charge_tolerance,
                energy_tolerance_j: controls.cell.energy_tolerance_j,
            },
        )
        .map_err(|error| Error::Invalid(SpectralPhotoionizationError::Group { index, error }))?;
        let rate = beam.number_flux_per_m2_s * area.value();
        if !nonnegative(rate) || (beam.number_flux_per_m2_s > 0.0 && rate == 0.0) {
            return Err(Error::Arithmetic("boundary photon rate"));
        }
        boundary.push([rate, photons]);
    }
    let zero = Point {
        neutral: 0.0,
        extent: Dd::ZERO,
        deposited: Dd::ZERO,
    };
    let mut points = buffer(cells.len(), zero)?;
    for (index, cell) in cells.iter().enumerate() {
        column_poll(&mut cancelled)?;
        let channel = cell.source.channel();
        let e = channel.electron_index();
        if cell.initial.thermodynamics() != channel.thermodynamics()
            || cell.initial.electron_index() != Some(e)
            || cell.initial.relative_charge_residual() > controls.cell.relative_charge_tolerance
        {
            return Err(cell_error(
                index,
                SpectralPhotoionizationError::Invalid("original source bank"),
            ));
        }
        check_temperature(cell.initial.electron_temperature().value(), controls.cell)
            .map_err(|error| cell_error(index, error))?;
        admit_caloric_interval(
            &cell.initial.thermodynamics().models()[e],
            Temperature::new(controls.cell.minimum_electron_temperature_k),
            Temperature::new(controls.cell.maximum_electron_temperature_k),
            &mut cancelled,
        )
        .map_err(|error| match error {
            KineticsErrorV1::Cancelled => Error::Cancelled,
            error => cell_error(index, SpectralPhotoionizationError::Admission(error)),
        })?;
        points[index].neutral = cell.initial.amounts()[channel.neutral_index()].value();
    }
    // Band-major layout makes the hot incident-to-exit sweep contiguous.
    let mut bands = Vec::new();
    bands
        .try_reserve_exact(pairs)
        .map_err(|_| Error::Arithmetic("column band allocation"))?;
    let mut active = false;
    for (g, group) in groups.iter().enumerate() {
        for (c, cell) in cells.iter().enumerate() {
            column_poll(&mut cancelled)?;
            let sigma = cell
                .source
                .cross_section(group.photon_energy)
                .map_err(|error| {
                    cell_error(
                        c,
                        SpectralPhotoionizationError::Group {
                            index: g,
                            error: PhotoionizationError::Source(error),
                        },
                    )
                })?
                .value();
            let band = Band {
                opacity_per_mol: (sigma * NA) / area.value(),
                excess_j_mol: if sigma == 0.0 {
                    0.0
                } else {
                    (group.photon_energy.value()
                        - cell.source.channel().ionization_energy_j_mol() / NA)
                        * NA
                },
            };
            if !nonnegative(band.opacity_per_mol)
                || (sigma > 0.0 && band.opacity_per_mol == 0.0)
                || !nonnegative(band.excess_j_mol)
            {
                return Err(Error::Arithmetic("prepared column band"));
            }
            active |= sigma > 0.0 && boundary[g][0] > 0.0 && points[c].neutral > 0.0;
            bands.push(band);
        }
    }
    let mut full = buffer(cells.len(), zero)?;
    let mut half = buffer(cells.len(), zero)?;
    let mut middle_temperatures = buffer(cells.len(), 0.0)?;
    let mut full_accounts = buffer(pairs, [Dd::ZERO; 2])?;
    let mut half_accounts = buffer(pairs, [Dd::ZERO; 2])?;
    let mut accounts = buffer(pairs, [Dd::ZERO; 2])?;
    let (mut time, mut step) = (0.0, controls.cell.initial_step_s);
    let (mut accepted, mut rejected, mut maximum_error) = (0, 0, 0.0_f64);
    if !active {
        for (g, incoming) in boundary.iter().enumerate() {
            for c in 0..cells.len() {
                accounts[g * cells.len() + c][1] = Dd::from_f64(incoming[1]);
            }
        }
        time = duration.value();
    }
    while time < duration.value() {
        column_poll(&mut cancelled)?;
        if accepted + rejected >= controls.cell.maximum_steps {
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
        half.copy_from_slice(&points);
        full_accounts.fill([Dd::ZERO; 2]);
        half_accounts.fill([Dd::ZERO; 2]);
        column_split_step(
            &mut full,
            dt,
            &bands,
            &boundary,
            &mut full_accounts,
            &mut cancelled,
        )?;
        column_split_step(
            &mut half,
            dt * 0.5,
            &bands,
            &boundary,
            &mut half_accounts,
            &mut cancelled,
        )?;
        for (c, cell) in cells.iter().enumerate() {
            middle_temperatures[c] =
                half[c].temperature(cell.initial, cell.source.channel().electron_index());
        }
        column_split_step(
            &mut half,
            dt * 0.5,
            &bands,
            &boundary,
            &mut half_accounts,
            &mut cancelled,
        )?;
        let amount_error = |x: f64, y: f64| {
            (x - y).abs()
                / (3.0
                    * (controls.cell.absolute_amount_tolerance_mol
                        + controls.cell.relative_tolerance * x.abs().max(y.abs())))
        };
        let mut error = 0.0_f64;
        for (c, cell) in cells.iter().enumerate() {
            column_poll(&mut cancelled)?;
            let e = cell.source.channel().electron_index();
            let tf = full[c].temperature(cell.initial, e);
            let th = half[c].temperature(cell.initial, e);
            if !positive(tf) || !positive(th) {
                return Err(Error::Arithmetic("column trial temperature"));
            }
            error = error
                .max(amount_error(full[c].neutral, half[c].neutral))
                .max(amount_error(
                    full[c].extent.to_f64(),
                    half[c].extent.to_f64(),
                ))
                .max(
                    (tf - th).abs()
                        / (3.0
                            * (controls.cell.temperature_tolerance_k
                                + controls.cell.relative_tolerance * tf.abs().max(th.abs()))),
                );
        }
        for (index, (f, h)) in full_accounts.iter().zip(&half_accounts).enumerate() {
            if index % 256 == 0 {
                column_poll(&mut cancelled)?;
            }
            for k in 0..2 {
                error = error.max(amount_error(f[k].to_f64() / NA, h[k].to_f64() / NA));
            }
        }
        if !error.is_finite() {
            return Err(Error::Arithmetic("column local error estimate"));
        }
        if error <= 1.0 {
            for (c, cell) in cells.iter().enumerate() {
                check_temperature(middle_temperatures[c], controls.cell)
                    .map_err(|error| cell_error(c, error))?;
                check_temperature(
                    half[c].temperature(cell.initial, cell.source.channel().electron_index()),
                    controls.cell,
                )
                .map_err(|error| cell_error(c, error))?;
            }
            for (index, (sum, local)) in accounts.iter_mut().zip(&half_accounts).enumerate() {
                if index % 256 == 0 {
                    column_poll(&mut cancelled)?;
                }
                for k in 0..2 {
                    sum[k] = sum[k] + local[k];
                }
            }
            core::mem::swap(&mut points, &mut half);
            time = end;
            accepted += 1;
            maximum_error = maximum_error.max(error);
            let factor = if error == 0.0 {
                2.0
            } else {
                (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
            };
            step = (dt * factor).clamp(controls.cell.minimum_step_s, controls.cell.maximum_step_s);
        } else {
            rejected += 1;
            step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
            if step < controls.cell.minimum_step_s {
                return Err(Error::StepBudget);
            }
        }
    }
    let mut endpoints = Vec::new();
    endpoints
        .try_reserve_exact(cells.len())
        .map_err(|_| Error::Arithmetic("column result allocation"))?;
    for (c, &cell) in cells.iter().enumerate() {
        column_poll(&mut cancelled)?;
        endpoints.push(publish_cell(
            c,
            cell,
            points[c],
            groups,
            &bands,
            &boundary,
            &accounts,
            cells.len(),
            controls.cell,
            &mut cancelled,
        )?);
    }
    let mut group_balances = Vec::new();
    group_balances
        .try_reserve_exact(groups.len())
        .map_err(|_| Error::Arithmetic("column group result allocation"))?;
    let mut totals = [Dd::ZERO; 6];
    for (g, incoming) in boundary.iter().enumerate() {
        column_poll(&mut cancelled)?;
        let mut sums = [Dd::ZERO; 4];
        for (c, cell) in endpoints.iter().enumerate() {
            if c % 256 == 0 {
                column_poll(&mut cancelled)?;
            }
            let b = cell.group_balances[g];
            for (sum, value) in sums.iter_mut().zip([
                b.absorbed_photons,
                b.chemical_storage_change_j,
                b.electron_sensible_energy_change_j,
                b.absorbed_energy_j,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let out = endpoints.last().map_or(incoming[1], |cell| {
            cell.group_balances[g].transmitted_photons
        });
        let residual = (sums[0] + Dd::from_f64(out) - Dd::from_f64(incoming[1])).to_f64();
        if !residual.is_finite()
            || residual.abs() > controls.cell.relative_photon_tolerance * incoming[1]
        {
            return Err(Error::PhotonBalance {
                group: g,
                residual_photons: residual,
            });
        }
        let values = [
            incoming[1],
            sums[0].to_f64(),
            out,
            sums[1].to_f64(),
            sums[2].to_f64(),
            sums[3].to_f64(),
        ];
        for (sum, value) in totals.iter_mut().zip(values) {
            *sum = *sum + Dd::from_f64(value);
        }
        group_balances.push(balance(values, (sums[1] + sums[2] - sums[3]).to_f64()));
    }
    // Check actual endpoint sensible energy once over the whole column. Per-band
    // deposited excess energy alone would not verify the published temperatures.
    totals[4] = endpoints.iter().fold(Dd::ZERO, |sum, cell| {
        sum + Dd::from_f64(cell.balance.electron_sensible_energy_change_j)
    });
    let residual = (totals[3] + totals[4] - totals[5]).to_f64();
    if totals.iter().any(|x| !x.to_f64().is_finite()) || !residual.is_finite() {
        return Err(Error::Arithmetic("whole-column energy or photon accounts"));
    }
    if residual.abs() > controls.cell.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: residual,
        });
    }
    column_poll(&mut cancelled)?;
    Ok(SpectralPhotoionizationColumnResult {
        cells: endpoints,
        groups: groups.to_vec(),
        group_balances,
        balance: balance(totals.map(Dd::to_f64), residual),
        area,
        elapsed: duration,
        controls,
        accepted_steps: accepted,
        rejected_steps: rejected,
        maximum_normalized_local_error: maximum_error,
    })
}

fn column_split_step(
    points: &mut [Point],
    dt: f64,
    bands: &[Band],
    boundary: &[[f64; 2]],
    accounts: &mut [[Dd; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SpectralPhotoionizationColumnError> {
    for g in 0..boundary.len() {
        let h = if g + 1 == boundary.len() {
            dt
        } else {
            dt * 0.5
        };
        band_column_step(points, g, h, bands, boundary, accounts, cancelled)?;
    }
    for g in (0..boundary.len().saturating_sub(1)).rev() {
        band_column_step(points, g, dt * 0.5, bands, boundary, accounts, cancelled)?;
    }
    Ok(())
}
fn band_column_step(
    points: &mut [Point],
    group: usize,
    dt: f64,
    bands: &[Band],
    boundary: &[[f64; 2]],
    accounts: &mut [[Dd; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SpectralPhotoionizationColumnError> {
    use SpectralPhotoionizationColumnError as Error;
    let mut photons = boundary[group][0] * dt;
    if !nonnegative(photons) {
        return Err(Error::Arithmetic("column substep photon exposure"));
    }
    let count = points.len();
    for (c, point) in points.iter_mut().enumerate() {
        column_poll(cancelled)?;
        let k = group * count + c;
        let band = bands[k];
        let account = &mut accounts[k];
        if photons == 0.0 || point.neutral == 0.0 || band.opacity_per_mol == 0.0 {
            account[1] = account[1] + Dd::from_f64(photons);
            continue;
        }
        let tau = band.opacity_per_mol * point.neutral;
        let dose = band.opacity_per_mol * (photons / NA);
        if !positive(tau) || !positive(dose) {
            return Err(Error::Arithmetic("positive column band depth or dose"));
        }
        let (lost_fraction, retained_fraction, _) = depletion(tau, dose).map_err(|error| {
            cell_error(
                c,
                SpectralPhotoionizationError::Group {
                    index: group,
                    error,
                },
            )
        })?;
        let lost = point.neutral * lost_fraction;
        let retained = point.neutral * retained_fraction;
        let out = transmitted_photons(photons, tau, dose);
        let energy = lost * band.excess_j_mol;
        if !positive(lost) || !nonnegative(retained) || !nonnegative(out) || !nonnegative(energy) {
            return Err(Error::Arithmetic("column band photon depletion"));
        }
        point.neutral = retained;
        point.extent = point.extent + Dd::from_f64(lost);
        point.deposited = point.deposited + Dd::from_f64(energy);
        account[0] = account[0] + Dd::from_f64(lost * NA);
        account[1] = account[1] + Dd::from_f64(out);
        photons = out;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn publish_cell(
    c: usize,
    cell: PhotoionizationCell<'_>,
    point: Point,
    groups: &[PhotonNumberFluxV1],
    bands: &[Band],
    boundary: &[[f64; 2]],
    accounts: &[[Dd; 2]],
    count: usize,
    controls: SpectralPhotoionizationSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpectralPhotoionizationColumnCellResult, SpectralPhotoionizationColumnError> {
    use SpectralPhotoionizationColumnError as Error;
    let channel = cell.source.channel();
    let [a, i, e] = [
        channel.neutral_index(),
        channel.ion_index(),
        channel.electron_index(),
    ];
    let mut group_balances = Vec::new();
    group_balances
        .try_reserve_exact(groups.len())
        .map_err(|_| Error::Arithmetic("cell result allocation"))?;
    let mut totals = [Dd::ZERO; 6];
    for (g, group) in groups.iter().enumerate() {
        column_poll(cancelled)?;
        let k = g * count + c;
        let incoming = if c == 0 {
            Dd::from_f64(boundary[g][1])
        } else {
            accounts[k - 1][1]
        };
        let [absorbed, out] = accounts[k];
        let residual = (absorbed + out - incoming).to_f64();
        if !residual.is_finite()
            || residual.abs() > controls.relative_photon_tolerance * incoming.to_f64()
        {
            return Err(cell_error(
                c,
                SpectralPhotoionizationError::PhotonBalance {
                    group: g,
                    residual_photons: residual,
                },
            ));
        }
        let extent = absorbed.to_f64() / NA;
        let values = [
            incoming.to_f64(),
            absorbed.to_f64(),
            out.to_f64(),
            extent * channel.ionization_energy_j_mol(),
            extent * bands[k].excess_j_mol,
            absorbed.to_f64() * group.photon_energy.value(),
        ];
        if values.iter().any(|x| !nonnegative(*x)) {
            return Err(Error::Arithmetic("cell photon or energy accounts"));
        }
        for (sum, value) in totals.iter_mut().zip(values) {
            *sum = *sum + Dd::from_f64(value);
        }
        group_balances.push(balance(values, values[4] + (values[3] - values[5])));
    }
    let extent = point.extent.to_f64();
    let amount_residual = (Dd::from_f64(point.neutral) + point.extent
        - Dd::from_f64(cell.initial.amounts()[a].value()))
    .to_f64();
    if !amount_residual.is_finite()
        || amount_residual.abs() > controls.absolute_amount_tolerance_mol
    {
        return Err(cell_error(
            c,
            SpectralPhotoionizationError::AmountBalance {
                residual_mol: amount_residual,
            },
        ));
    }
    let te = point.temperature(cell.initial, e);
    let sensible = 1.5
        * R
        * (cell.initial.amounts()[e].value() * (te - cell.initial.electron_temperature().value())
            + extent * te);
    totals[4] = Dd::from_f64(sensible);
    let residual = (totals[3] + totals[4] - totals[5]).to_f64();
    if totals.iter().any(|x| !x.to_f64().is_finite()) || !residual.is_finite() {
        return Err(Error::Arithmetic("cell total photon or energy accounts"));
    }
    if residual.abs() > controls.energy_tolerance_j {
        return Err(cell_error(
            c,
            SpectralPhotoionizationError::EnergyBalance {
                residual_j: residual,
            },
        ));
    }
    let mut amounts = cell.initial.amounts().to_vec();
    amounts[a] = Amount::new(point.neutral);
    amounts[i] = Amount::new(cell.initial.amounts()[i].value() + extent);
    amounts[e] = Amount::new(cell.initial.amounts()[e].value() + extent);
    let state = FrozenIonizedGasStateV1::new(
        channel.thermodynamics(),
        &amounts,
        cell.initial.volume(),
        cell.initial.heavy_temperature(),
        Temperature::new(te),
        controls.relative_charge_tolerance,
        &mut *cancelled,
    )
    .map_err(|error| match error {
        FrozenIonizedGasStateErrorV1::Cancelled => Error::Cancelled,
        error => cell_error(c, SpectralPhotoionizationError::State(error)),
    })?;
    Ok(SpectralPhotoionizationColumnCellResult {
        initial: cell.initial.clone(),
        state,
        source: cell.source.clone(),
        extent: Amount::new(extent),
        group_balances,
        balance: balance(totals.map(Dd::to_f64), residual),
    })
}

fn balance(values: [f64; 6], residual: f64) -> PhotoionizationBalance {
    PhotoionizationBalance {
        incident_photons: values[0],
        absorbed_photons: values[1],
        transmitted_photons: values[2],
        chemical_storage_change_j: values[3],
        electron_sensible_energy_change_j: values[4],
        absorbed_energy_j: values[5],
        energy_residual_j: residual,
    }
}
fn buffer<T: Clone>(count: usize, value: T) -> Result<Vec<T>, SpectralPhotoionizationColumnError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| SpectralPhotoionizationColumnError::Arithmetic("column work allocation"))?;
    result.resize(count, value);
    Ok(result)
}
fn cell_error(
    index: usize,
    error: SpectralPhotoionizationError,
) -> SpectralPhotoionizationColumnError {
    match error {
        SpectralPhotoionizationError::Cancelled => SpectralPhotoionizationColumnError::Cancelled,
        error => SpectralPhotoionizationColumnError::Cell { index, error },
    }
}
fn column_poll(
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SpectralPhotoionizationColumnError> {
    if cancelled() {
        Err(SpectralPhotoionizationColumnError::Cancelled)
    } else {
        Ok(())
    }
}
