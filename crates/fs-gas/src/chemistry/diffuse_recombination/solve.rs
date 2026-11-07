use super::*;
use fs_math::dd::Dd;
use fs_math::det::{exp, expm1};
use fs_qty::Energy;
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::spectrum::RecombinationSpectrumV1;

type Error = DiffuseRecombinationError;
const CELL: usize = 10;
const FACE: usize = 4;

// Positive halves of the standard Gauss-Legendre 8- and 16-point rules.
const X8: [f64; 4] = [
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const W8: [f64; 4] = [
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];
const X16: [f64; 8] = [
    0.095_012_509_837_637_44,
    0.281_603_550_779_258_9,
    0.458_016_777_657_227_4,
    0.617_876_244_402_643_8,
    0.755_404_408_355_003,
    0.865_631_202_387_831_8,
    0.944_575_023_073_232_6,
    0.989_400_934_991_649_9,
];
const W16: [f64; 8] = [
    0.189_450_610_455_068_5,
    0.182_603_415_044_923_6,
    0.169_156_519_395_002_5,
    0.149_595_988_816_576_7,
    0.124_628_971_255_533_9,
    0.095_158_511_682_492_78,
    0.062_253_523_938_647_89,
    0.027_152_459_411_754_096,
];

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let increment = value - self.correction;
        let next = self.value + increment;
        self.correction = (next - self.value) - increment;
        self.value = next;
    }
}

/// The transfer kernel retains source/geometry bindings and raw populations;
/// it never clones a frozen state or thermodynamic bank during quadrature.
struct Binding<'a> {
    source: &'a AtomicRadiativeRecombinationV1,
    spectrum: RecombinationSpectrumV1<'a>,
    opacity_per_sigma: f64,
    pair_density_volume: f64,
    gap: f64,
    expected: [f64; 3], // emitted number, photon energy, captured energy
    cutoff: f64,
    tail: [f64; 3],
}
struct Prepared<'a> {
    cells: Vec<Binding<'a>>,
    shared_upper: f64,
    controls: DiffuseRecombinationSolve,
    panels: usize,
    work: usize,
}

struct Scratch {
    opacity: Vec<f64>,
    emission: Vec<f64>,
    fine_point: Vec<Sum>,
    angular_point: Vec<Sum>,
    fine_panel: Vec<Sum>,
    coarse_panel: Vec<Sum>,
    angular_panel: Vec<Sum>,
    values: Vec<Sum>,
    spectral_errors: Vec<Sum>,
    angular_errors: Vec<Sum>,
    tail: Vec<f64>,
}
impl Scratch {
    fn new(cells: usize) -> Result<Self, Error> {
        let n = CELL * cells + FACE * (cells + 1);
        Ok(Self {
            opacity: buffer(cells, 0.0)?,
            emission: buffer(cells, 0.0)?,
            fine_point: buffer(n, Sum::default())?,
            angular_point: buffer(n, Sum::default())?,
            fine_panel: buffer(n, Sum::default())?,
            coarse_panel: buffer(n, Sum::default())?,
            angular_panel: buffer(n, Sum::default())?,
            values: buffer(n, Sum::default())?,
            spectral_errors: buffer(n, Sum::default())?,
            angular_errors: buffer(n, Sum::default())?,
            tail: buffer(n, 0.0)?,
        })
    }
}

pub(super) fn evaluate(
    cells: &[DiffuseRecombinationCell<'_>],
    area: Area,
    controls: DiffuseRecombinationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DiffuseRecombinationResult, Error> {
    let mut prepared = Prepared::new(cells, area, controls, &mut cancelled)?;
    let mut scratch = Scratch::new(cells.len())?;
    let (mut extent, mut refinement, mut angular) = (32.0, 1_usize, 1_usize);
    let mut normalized = [0.0; 3];
    if prepared.cells.iter().any(|cell| cell.expected[0] > 0.0) {
        loop {
            poll(&mut cancelled)?;
            prepared.cutoffs(extent)?;
            let grid = prepared.grid(&mut cancelled)?;
            let fine_angles = angles(angular, &X16, &W16)?;
            let coarse_angles = angles(angular, &X8, &W8)?;
            prepared.integrate(
                &grid,
                refinement,
                &fine_angles,
                &coarse_angles,
                &mut scratch,
                &mut cancelled,
            )?;
            prepared.tail_vector(&mut scratch.tail);
            normalized = prepared.errors(&scratch)?;
            if normalized[2] > 0.25 {
                if extent >= 1024.0
                    || prepared
                        .cells
                        .iter()
                        .all(|c| c.pair_density_volume == 0.0 || c.cutoff == prepared.shared_upper)
                {
                    return Err(Error::UnresolvedTail {
                        maximum_normalized_bound: normalized[2],
                    });
                }
                extent *= 2.0;
                continue;
            }
            if normalized[0] > 0.25 {
                refinement = refinement.checked_mul(2).ok_or(Error::SpectralBudget)?;
                continue;
            }
            if normalized[1] > 0.25 {
                if angular >= controls.maximum_angular_subdivisions {
                    return Err(Error::AngularBudget);
                }
                angular = angular
                    .saturating_mul(2)
                    .min(controls.maximum_angular_subdivisions);
                continue;
            }
            break;
        }
    }
    prepared.publish(cells, area, &scratch, normalized, angular, &mut cancelled)
}

impl<'a> Prepared<'a> {
    fn new(
        cells: &[DiffuseRecombinationCell<'a>],
        area: Area,
        controls: DiffuseRecombinationSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if !positive(area.value())
            || !positive(controls.relative_tolerance)
            || controls.relative_tolerance >= 1.0
            || !positive(controls.absolute_photon_rate_tolerance_per_s)
            || !positive(controls.absolute_power_tolerance_w)
            || !positive(controls.relative_charge_tolerance)
            || controls.relative_charge_tolerance >= 1.0
            || controls.maximum_cells > 4096
            || cells.len() > controls.maximum_cells
            || controls.maximum_spectral_panels == 0
            || !(1..=1024).contains(&controls.maximum_angular_subdivisions)
            || controls.maximum_work == 0
            || !positive(controls.source_quadrature.relative_tolerance)
            || controls.source_quadrature.relative_tolerance > controls.relative_tolerance * 0.125
            || controls.source_quadrature.maximum_subdivisions == 0
        {
            return Err(Error::Invalid(
                "geometry, source tolerance, cell or quadrature controls",
            ));
        }
        let mut result = Self {
            cells: Vec::new(),
            shared_upper: f64::INFINITY,
            controls,
            panels: 0,
            work: 0,
        };
        result
            .cells
            .try_reserve_exact(cells.len())
            .map_err(|_| Error::Allocation)?;
        for (index, cell) in cells.iter().enumerate() {
            poll(cancelled)?;
            result.spend(1)?;
            let state = cell.state;
            let channel = cell.source.source().channel();
            if state.thermodynamics() != channel.thermodynamics()
                || state.electron_index() != Some(channel.electron_index())
                || state.relative_charge_residual() > controls.relative_charge_tolerance
            {
                return Err(Error::Invalid("state and original reciprocal source bank"));
            }
            let spectrum = cell
                .source
                .spectrum(state.electron_temperature())
                .map_err(|error| source_error(index, error))?;
            let lower = state.amounts()[channel.neutral_index()].value();
            let upper = state.amounts()[channel.ion_index()].value();
            let electron = state.amounts()[channel.electron_index()].value();
            let opacity_per_sigma = NA * (lower / area.value());
            let pair_density_volume = (NA * upper) * (NA * (electron / state.volume().value()));
            if !nonnegative(opacity_per_sigma)
                || !nonnegative(pair_density_volume)
                || (lower > 0.0 && opacity_per_sigma == 0.0)
                || (upper > 0.0 && electron > 0.0 && pair_density_volume == 0.0)
            {
                return Err(Error::Arithmetic(
                    "positive source density or column inventory",
                ));
            }
            let gap = channel.ionization_energy_j_mol() / NA;
            let mut expected = [0.0; 3];
            if pair_density_volume > 0.0 {
                let moments = cell
                    .source
                    .evaluate(
                        state.electron_temperature(),
                        controls.source_quadrature,
                        &mut *cancelled,
                    )
                    .map_err(|error| source_error(index, error))?;
                expected = [
                    moments.coefficient_m3_s,
                    moments.emitted_photon_energy_j_m3_s,
                    moments.captured_kinetic_energy_j_m3_s,
                ]
                .map(|value| value * pair_density_volume);
                if expected.iter().any(|x| !nonnegative(*x))
                    || (moments.coefficient_m3_s > 0.0 && expected[0] == 0.0)
                {
                    return Err(Error::Arithmetic("original emission moments"));
                }
            }
            if opacity_per_sigma > 0.0 || expected[0] > 0.0 {
                result.shared_upper = result
                    .shared_upper
                    .min(cell.source.source().maximum_energy().value());
            }
            result.cells.push(Binding {
                source: cell.source,
                spectrum,
                opacity_per_sigma,
                pair_density_volume,
                gap,
                expected,
                cutoff: 0.0,
                tail: [0.0; 3],
            });
        }
        Ok(result)
    }

    fn spend(&mut self, amount: usize) -> Result<(), Error> {
        self.work = self
            .work
            .checked_add(amount)
            .filter(|n| *n <= self.controls.maximum_work)
            .ok_or(Error::WorkBudget)?;
        Ok(())
    }

    fn cutoffs(&mut self, thermal_extent: f64) -> Result<(), Error> {
        for (index, cell) in self.cells.iter_mut().enumerate() {
            if cell.expected[0] == 0.0 {
                continue;
            }
            let end = cell.source.source().threshold().value()
                + thermal_extent * cell.spectrum.thermal_energy().value();
            if !positive(end) {
                return Err(Error::Arithmetic("emitter spectral cutoff"));
            }
            cell.cutoff = end.min(self.shared_upper);
            let tail = cell
                .spectrum
                .tail_bound(Energy::new(cell.cutoff))
                .map_err(|error| source_error(index, error))?;
            cell.tail = [
                tail.coefficient_m3_s,
                tail.emitted_photon_energy_j_m3_s,
                tail.captured_kinetic_energy_j_m3_s,
            ]
            .map(|v| (v * cell.pair_density_volume).max(f64::from_bits(1)));
            if cell.tail.iter().any(|x| !positive(*x)) {
                return Err(Error::Arithmetic("emission tail scaling"));
            }
        }
        Ok(())
    }

    fn grid(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>, Error> {
        let lower = self
            .cells
            .iter()
            .filter(|c| c.expected[0] > 0.0)
            .map(|c| c.source.source().threshold().value())
            .fold(f64::INFINITY, f64::min);
        let upper = self.cells.iter().map(|c| c.cutoff).fold(0.0, f64::max);
        if lower >= upper {
            return Ok(Vec::new());
        }
        let mut grid = Vec::new();
        let mut push = |x: f64| -> Result<(), Error> {
            if x < lower || x > upper {
                return Ok(());
            }
            if grid.len() > self.controls.maximum_spectral_panels {
                return Err(Error::SpectralBudget);
            }
            grid.try_reserve(1).map_err(|_| Error::Allocation)?;
            grid.push(x);
            Ok(())
        };
        push(lower)?;
        push(upper)?;
        for cell in &self.cells {
            poll(cancelled)?;
            let photo = cell.source.source();
            push(photo.threshold().value())?;
            push(cell.cutoff)?;
            if let Some(points) = photo.points() {
                for point in points {
                    poll(cancelled)?;
                    push(point.photon_energy.value())?;
                }
            }
            if cell.expected[0] > 0.0 {
                let mut x = 1.0;
                while x <= 1024.0 {
                    let energy =
                        photo.threshold().value() + x * cell.spectrum.thermal_energy().value();
                    if energy >= cell.cutoff {
                        break;
                    }
                    push(energy)?;
                    x *= 2.0;
                }
            }
        }
        grid.sort_by(f64::total_cmp);
        grid.dedup();
        Ok(grid)
    }

    fn integrate(
        &mut self,
        grid: &[f64],
        refinement: usize,
        fine_angles: &[(f64, f64)],
        coarse_angles: &[(f64, f64)],
        scratch: &mut Scratch,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        scratch.values.fill(Sum::default());
        scratch.spectral_errors.fill(Sum::default());
        scratch.angular_errors.fill(Sum::default());
        let needed = grid
            .len()
            .saturating_sub(1)
            .checked_mul(refinement)
            .ok_or(Error::SpectralBudget)?;
        if needed
            > self
                .controls
                .maximum_spectral_panels
                .saturating_sub(self.panels)
        {
            return Err(Error::SpectralBudget);
        }
        self.panels += needed;
        for interval in grid.windows(2) {
            let width = (interval[1] - interval[0]) / refinement as f64;
            for panel in 0..refinement {
                poll(cancelled)?;
                let lo = interval[0] + panel as f64 * width;
                let hi = if panel + 1 == refinement {
                    interval[1]
                } else {
                    interval[0] + (panel + 1) as f64 * width
                };
                let half = (hi - lo) * 0.5;
                let mid = lo + half;
                if !positive(half) {
                    return Err(Error::Arithmetic("positive spectral panel"));
                }
                scratch.fine_panel.fill(Sum::default());
                scratch.coarse_panel.fill(Sum::default());
                scratch.angular_panel.fill(Sum::default());
                for (&x, &w) in X16.iter().zip(&W16) {
                    for sign in [-1.0, 1.0] {
                        self.point(
                            mid + sign * half * x,
                            fine_angles,
                            Some(coarse_angles),
                            scratch,
                            cancelled,
                        )?;
                        for j in 0..scratch.values.len() {
                            scratch.fine_panel[j].add(half * w * scratch.fine_point[j].value);
                            scratch.angular_panel[j].add(half * w * scratch.angular_point[j].value);
                        }
                    }
                }
                for (&x, &w) in X8.iter().zip(&W8) {
                    for sign in [-1.0, 1.0] {
                        self.point(mid + sign * half * x, fine_angles, None, scratch, cancelled)?;
                        for j in 0..scratch.values.len() {
                            scratch.coarse_panel[j].add(half * w * scratch.fine_point[j].value);
                        }
                    }
                }
                for j in 0..scratch.values.len() {
                    scratch.values[j].add(scratch.fine_panel[j].value);
                    scratch.spectral_errors[j]
                        .add((scratch.fine_panel[j].value - scratch.coarse_panel[j].value).abs());
                    scratch.angular_errors[j]
                        .add((scratch.fine_panel[j].value - scratch.angular_panel[j].value).abs());
                }
            }
        }
        Ok(())
    }

    fn point(
        &mut self,
        energy: f64,
        fine_angles: &[(f64, f64)],
        coarse_angles: Option<&[(f64, f64)]>,
        scratch: &mut Scratch,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.spend(self.cells.len())?;
        scratch.fine_point.fill(Sum::default());
        scratch.angular_point.fill(Sum::default());
        for (index, cell) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            let mut sigma = 0.0;
            let (mut photons, mut captured) = (0.0, 0.0);
            if cell.expected[0] > 0.0 && energy <= cell.cutoff {
                let density = cell
                    .spectrum
                    .evaluate(Energy::new(energy))
                    .map_err(|error| source_error(index, error))?;
                sigma = density.photoionization_cross_section.value();
                photons = density.coefficient_m3_s_per_j * cell.pair_density_volume;
                captured = density.captured_kinetic_energy_m3_s * cell.pair_density_volume;
            } else if cell.opacity_per_sigma > 0.0 {
                // A cold emitter can have a rigorously omitted high-energy
                // tail while its original absorption remains fully active.
                sigma = cell
                    .source
                    .source()
                    .cross_section(Energy::new(energy))
                    .map_err(|error| Error::Absorption { index, error })?
                    .value();
            }
            let depth = sigma * cell.opacity_per_sigma;
            if !nonnegative(depth)
                || !nonnegative(photons)
                || !nonnegative(captured)
                || (sigma > 0.0 && cell.opacity_per_sigma > 0.0 && depth == 0.0)
            {
                return Err(Error::Arithmetic("spectral optical depth or emissivity"));
            }
            scratch.opacity[index] = depth;
            scratch.emission[index] = photons;
            for (slot, value) in [photons, energy * photons, captured]
                .into_iter()
                .enumerate()
            {
                if !nonnegative(value) {
                    return Err(Error::Arithmetic("spectral emission moment"));
                }
                scratch.fine_point[CELL * index + slot].value = value;
                scratch.angular_point[CELL * index + slot].value = value;
            }
        }
        self.sweep(
            energy,
            fine_angles,
            &scratch.opacity,
            &scratch.emission,
            &mut scratch.fine_point,
            cancelled,
        )?;
        if let Some(angles) = coarse_angles {
            self.sweep(
                energy,
                angles,
                &scratch.opacity,
                &scratch.emission,
                &mut scratch.angular_point,
                cancelled,
            )?;
        }
        Ok(())
    }

    fn sweep(
        &mut self,
        energy: f64,
        angles: &[(f64, f64)],
        opacity: &[f64],
        emission: &[f64],
        result: &mut [Sum],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let n = self.cells.len();
        self.spend(n.checked_mul(2 * angles.len()).ok_or(Error::WorkBudget)?)?;
        for &(mu, weight) in angles {
            for right in [true, false] {
                let mut incoming = 0.0;
                for step in 0..n {
                    poll(cancelled)?;
                    let c = if right { step } else { n - step - 1 };
                    let (transmission, absorption, escape, self_absorption) =
                        characteristic(opacity[c] / mu);
                    let local = emission[c] * (0.5 * weight);
                    let self_absorbed = local * self_absorption;
                    let external_absorbed = incoming * absorption;
                    let outgoing = incoming * transmission + local * escape;
                    let absorbed = self_absorbed + external_absorbed;
                    let heat = if absorbed == 0.0 {
                        0.0
                    } else {
                        (energy - self.cells[c].gap) * absorbed
                    };
                    if !nonnegative(outgoing) || !nonnegative(heat) {
                        return Err(Error::Arithmetic("formal spectral transport"));
                    }
                    result[CELL * c + 3].add(self_absorbed);
                    result[CELL * c + 4].add(external_absorbed);
                    result[CELL * c + 5].add(energy * absorbed);
                    result[CELL * c + 6].add(heat);
                    // Self-reabsorption returns the source photon and its
                    // captured electron energy to the same channel. Cancel
                    // that internal traffic before rounding: subtracting two
                    // integrated gross rates would erase physical escape in
                    // an optically thick slab.
                    let net = external_absorbed - local * escape;
                    result[CELL * c + 7].add(net);
                    result[CELL * c + 8].add((energy - self.cells[c].gap) * net);
                    result[CELL * c + 9].add(energy * net);
                    let face = if right { c + 1 } else { c };
                    let slot = CELL * n + FACE * face + if right { 0 } else { 2 };
                    result[slot].add(outgoing);
                    result[slot + 1].add(energy * outgoing);
                    incoming = outgoing;
                }
            }
        }
        Ok(())
    }

    fn tail_vector(&self, output: &mut [f64]) {
        output.fill(0.0);
        let mut total = [Dd::ZERO; 2];
        for cell in &self.cells {
            for j in 0..2 {
                total[j] = total[j] + Dd::from_f64(cell.tail[j]);
            }
        }
        for (c, cell) in self.cells.iter().enumerate() {
            output[CELL * c..CELL * c + 3].copy_from_slice(&cell.tail);
            if cell.opacity_per_sigma > 0.0 {
                output[CELL * c + 3] = cell.tail[0];
                output[CELL * c + 4] = (total[0] - Dd::from_f64(cell.tail[0])).to_f64().max(0.0);
                output[CELL * c + 5] = total[1].to_f64();
                output[CELL * c + 6] = total[1].to_f64();
            }
            output[CELL * c + 7] = total[0].to_f64();
            output[CELL * c + 8] = total[1].to_f64();
            output[CELL * c + 9] = total[1].to_f64();
        }
        // Only sources on the corresponding side can cross each interface.
        for right in [true, false] {
            let mut sum = [Dd::ZERO; 2];
            for step in 0..self.cells.len() {
                let c = if right {
                    step
                } else {
                    self.cells.len() - step - 1
                };
                for j in 0..2 {
                    sum[j] = sum[j] + Dd::from_f64(self.cells[c].tail[j]);
                }
                let face = if right { c + 1 } else { c };
                let slot = CELL * self.cells.len() + FACE * face + if right { 0 } else { 2 };
                for j in 0..2 {
                    output[slot + j] = 0.5 * sum[j].to_f64();
                }
            }
        }
    }

    fn allowance(&self, index: usize, value: f64) -> f64 {
        let number = if index < CELL * self.cells.len() {
            matches!(index % CELL, 0 | 3 | 4 | 7)
        } else {
            (index - CELL * self.cells.len()) % 2 == 0
        };
        (if number {
            self.controls.absolute_photon_rate_tolerance_per_s
        } else {
            self.controls.absolute_power_tolerance_w
        }) + self.controls.relative_tolerance * value.abs()
    }

    fn errors(&self, scratch: &Scratch) -> Result<[f64; 3], Error> {
        let mut errors = [0.0_f64; 3];
        for (j, value) in scratch.values.iter().enumerate() {
            let signed = j < CELL * self.cells.len() && j % CELL >= 7;
            if !value.value.is_finite() || (!signed && value.value < 0.0) {
                return Err(Error::Arithmetic("integrated radiation moment"));
            }
            let allowance = self.allowance(j, value.value);
            for (k, error) in [
                scratch.spectral_errors[j].value,
                scratch.angular_errors[j].value,
                scratch.tail[j],
            ]
            .into_iter()
            .enumerate()
            {
                let normalized = error / allowance;
                if !nonnegative(normalized) {
                    return Err(Error::Arithmetic("transport error estimate"));
                }
                errors[k] = errors[k].max(normalized);
            }
        }
        for (c, cell) in self.cells.iter().enumerate() {
            for j in 0..3 {
                // Independently integrated original source moments guard a
                // missed emission peak. Explicitly bounded omitted tails are
                // not mistaken for quadrature error.
                let mismatch = ((scratch.values[CELL * c + j].value - cell.expected[j]).abs()
                    - cell.tail[j])
                    .max(0.0);
                errors[0] =
                    errors[0].max(mismatch / self.allowance(CELL * c + j, cell.expected[j]));
            }
        }
        Ok(errors)
    }

    fn publish(
        &self,
        originals: &[DiffuseRecombinationCell<'_>],
        area: Area,
        scratch: &Scratch,
        normalized: [f64; 3],
        angular: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DiffuseRecombinationResult, Error> {
        let n = self.cells.len();
        let mut interfaces = Vec::new();
        interfaces
            .try_reserve_exact(n + 1)
            .map_err(|_| Error::Allocation)?;
        for face in 0..=n {
            let values = &scratch.values[CELL * n + FACE * face..CELL * n + FACE * (face + 1)];
            interfaces.push(DiffuseRecombinationFlux {
                right_photons_per_s: values[0].value,
                right_power_w: values[1].value,
                left_photons_per_s: values[2].value,
                left_power_w: values[3].value,
            });
        }
        let mut cells = Vec::new();
        cells.try_reserve_exact(n).map_err(|_| Error::Allocation)?;
        let mut sums = [Dd::ZERO; 6];
        for (c, binding) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            let v = &scratch.values[CELL * c..CELL * (c + 1)];
            let absorbed = Dd::from_f64(v[3].value) + Dd::from_f64(v[4].value);
            let net_events = Dd::from_f64(v[7].value);
            let net_electron = Dd::from_f64(v[8].value);
            let chemical = binding.gap * net_events.to_f64();
            let net_gas = Dd::from_f64(v[9].value);
            let gas_residual = (Dd::from_f64(chemical) + net_electron - net_gas).to_f64();
            let left = interfaces[c];
            let right = interfaces[c + 1];
            let net_photon_flux_residual = (net_events
                - Dd::from_f64(left.right_photons_per_s)
                - Dd::from_f64(right.left_photons_per_s)
                + Dd::from_f64(right.right_photons_per_s)
                + Dd::from_f64(left.left_photons_per_s))
            .to_f64();
            let net_power_flux_residual =
                (net_gas - Dd::from_f64(left.right_power_w) - Dd::from_f64(right.left_power_w)
                    + Dd::from_f64(right.right_power_w)
                    + Dd::from_f64(left.left_power_w))
                .to_f64();
            // Scale only the roundoff allowance to the actual interface
            // traffic. Enormous local self-recycling cannot hide a small,
            // incorrectly lost net reaction or thermal source here.
            let photon_roundoff = 128.0
                * f64::EPSILON
                * (left.right_photons_per_s
                    + right.left_photons_per_s
                    + right.right_photons_per_s
                    + left.left_photons_per_s);
            let power_roundoff = 128.0
                * f64::EPSILON
                * (left.right_power_w
                    + right.left_power_w
                    + right.right_power_w
                    + left.left_power_w);
            if !net_photon_flux_residual.is_finite()
                || !net_power_flux_residual.is_finite()
                || net_photon_flux_residual.abs()
                    > self.allowance(CELL * c + 7, net_events.to_f64()) + photon_roundoff
                || net_power_flux_residual.abs()
                    > self.allowance(CELL * c + 9, net_gas.to_f64()) + power_roundoff
            {
                return Err(Error::Balance {
                    photon_residual_per_s: net_photon_flux_residual,
                    energy_residual_w: net_power_flux_residual,
                });
            }
            let photon_residual = (Dd::from_f64(v[0].value)
                + Dd::from_f64(left.right_photons_per_s)
                + Dd::from_f64(right.left_photons_per_s)
                - absorbed
                - Dd::from_f64(right.right_photons_per_s)
                - Dd::from_f64(left.left_photons_per_s))
            .to_f64();
            let energy_residual = (Dd::from_f64(v[1].value)
                + Dd::from_f64(left.right_power_w)
                + Dd::from_f64(right.left_power_w)
                - Dd::from_f64(v[5].value)
                - Dd::from_f64(right.right_power_w)
                - Dd::from_f64(left.left_power_w))
            .to_f64();
            self.check_balance(
                photon_residual,
                energy_residual,
                v[0].value + left.right_photons_per_s + right.left_photons_per_s,
                v[1].value + left.right_power_w + right.left_power_w,
            )?;
            self.check_balance(0.0, gas_residual, 0.0, v[1].value + v[5].value)?;
            let balance = DiffuseRecombinationCellBalance {
                emitted_photons_per_s: v[0].value,
                emitted_power_w: v[1].value,
                captured_electron_power_w: v[2].value,
                self_absorbed_photons_per_s: v[3].value,
                nonlocal_absorbed_photons_per_s: v[4].value,
                absorbed_power_w: v[5].value,
                photoelectron_power_w: v[6].value,
                net_ionization_rate_mol_per_s: net_events.to_f64() / NA,
                net_electron_power_w: net_electron.to_f64(),
                chemical_storage_rate_w: chemical,
                net_gas_power_w: net_gas.to_f64(),
                photon_residual_per_s: photon_residual,
                energy_residual_w: energy_residual,
                gas_energy_residual_w: gas_residual,
                net_photon_flux_residual_per_s: net_photon_flux_residual,
                net_power_flux_residual_w: net_power_flux_residual,
            };
            for (sum, value) in sums.iter_mut().zip([
                v[0].value,
                v[1].value,
                absorbed.to_f64(),
                v[5].value,
                binding.tail[0],
                binding.tail[1],
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            cells.push(DiffuseRecombinationCellResult {
                state: originals[c].state.clone(),
                source: originals[c].source.clone(),
                balance,
                emission_cutoff_j: binding.cutoff,
            });
        }
        let escaped = Dd::from_f64(interfaces[0].left_photons_per_s)
            + Dd::from_f64(interfaces[n].right_photons_per_s);
        let escaped_power =
            Dd::from_f64(interfaces[0].left_power_w) + Dd::from_f64(interfaces[n].right_power_w);
        let photon_residual = (sums[0] - sums[2] - escaped).to_f64();
        let energy_residual = (sums[1] - sums[3] - escaped_power).to_f64();
        self.check_balance(
            photon_residual,
            energy_residual,
            sums[0].to_f64(),
            sums[1].to_f64(),
        )?;
        poll(cancelled)?;
        Ok(DiffuseRecombinationResult {
            cells,
            interfaces,
            area,
            controls: self.controls,
            balance: DiffuseRecombinationBalance {
                emitted_photons_per_s: sums[0].to_f64(),
                emitted_power_w: sums[1].to_f64(),
                absorbed_photons_per_s: sums[2].to_f64(),
                absorbed_power_w: sums[3].to_f64(),
                escaped_photons_per_s: escaped.to_f64(),
                escaped_power_w: escaped_power.to_f64(),
                photon_residual_per_s: photon_residual,
                energy_residual_w: energy_residual,
                omitted_photon_rate_bound_per_s: sums[4].to_f64(),
                omitted_power_bound_w: sums[5].to_f64(),
                maximum_normalized_spectral_estimate: normalized[0],
                maximum_normalized_angular_estimate: normalized[1],
                maximum_normalized_tail_bound: normalized[2],
                spectral_panels: self.panels,
                angular_subdivisions: angular,
                work: self.work,
            },
        })
    }

    fn check_balance(
        &self,
        photons: f64,
        energy: f64,
        photon_scale: f64,
        energy_scale: f64,
    ) -> Result<(), Error> {
        if !photons.is_finite()
            || !energy.is_finite()
            || photons.abs()
                > self.controls.absolute_photon_rate_tolerance_per_s
                    + self.controls.relative_tolerance * photon_scale
            || energy.abs()
                > self.controls.absolute_power_tolerance_w
                    + self.controls.relative_tolerance * energy_scale
        {
            return Err(Error::Balance {
                photon_residual_per_s: photons,
                energy_residual_w: energy,
            });
        }
        Ok(())
    }
}

/// Uniform volume source: emission escape phi(tau)=(1-exp(-tau))/tau.
/// Retain small absorption directly, rather than subtracting it from one.
fn characteristic(depth: f64) -> (f64, f64, f64, f64) {
    if depth == 0.0 {
        return (1.0, 0.0, 1.0, 0.0);
    }
    if depth.is_infinite() {
        return (0.0, 1.0, 0.0, 1.0);
    }
    let absorption = -expm1(-depth);
    let transmission = exp(-depth);
    let self_absorption = if depth < 1e-3 {
        depth
            * (0.5
                + depth
                    * (-1.0 / 6.0 + depth * (1.0 / 24.0 + depth * (-1.0 / 120.0 + depth / 720.0))))
    } else {
        1.0 - absorption / depth
    };
    (
        transmission,
        absorption,
        absorption / depth,
        self_absorption,
    )
}

fn angles<const N: usize>(
    panels: usize,
    nodes: &[f64; N],
    weights: &[f64; N],
) -> Result<Vec<(f64, f64)>, Error> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(2 * N * panels)
        .map_err(|_| Error::Allocation)?;
    let half = 0.5 / panels as f64;
    for panel in 0..panels {
        let mid = (panel as f64 + 0.5) / panels as f64;
        for (&x, &w) in nodes.iter().zip(weights) {
            for sign in [-1.0, 1.0] {
                let t = mid + sign * half * x;
                output.push((t * t * t * t, half * w * 4.0 * t * t * t));
            }
        }
    }
    Ok(output)
}
fn buffer<T: Clone>(n: usize, value: T) -> Result<Vec<T>, Error> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| Error::Allocation)?;
    v.resize(n, value);
    Ok(v)
}
fn source_error(index: usize, error: RadiativeRecombinationErrorV1) -> Error {
    if matches!(error, RadiativeRecombinationErrorV1::Cancelled) {
        Error::Cancelled
    } else {
        Error::Source { index, error }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
