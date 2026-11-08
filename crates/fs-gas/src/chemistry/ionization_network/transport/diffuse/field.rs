use super::super::field::PreparedField;
use super::*;
use fs_math::dd::Dd;
use fs_qty::{Energy, Temperature};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::{
    RadiativeRecombinationErrorV1, spectrum::RecombinationSpectrumV1,
};

mod spectral;
type Error = DiffuseRecombinationError;

// Per channel: emitted number/energy/captured energy, same-source/cross-source/
// nonlocal absorption, absorbed energy/photoelectron heat, signed net number/
// electron energy/gas energy. Cell totals are retained independently.
const CHANNEL: usize = 11;
const CELL: usize = 3;
const FACE: usize = 4;

struct Binding<'a> {
    source: &'a AtomicRadiativeRecombinationV1,
    spectrum: RecombinationSpectrumV1<'a>,
    lower: usize,
    upper: usize,
    electron: usize,
    opacity_per_sigma: f64,
    pair_density_volume: f64,
    gap: f64,
    expected: [f64; 3],
    cutoff: f64,
    tail: [f64; 3],
    volume: f64,
}

struct Prepared<'a> {
    channels: Vec<Binding<'a>>,
    offsets: Vec<usize>,
    area: f64,
    shared_upper: f64,
    controls: DiffuseRecombinationSolve,
    panels: usize,
    work: usize,
}

struct Scratch {
    opacity: Vec<f64>,
    fractions: Vec<Dd>,
    emission: Vec<f64>,
    cross: Vec<Dd>,
    cell_opacity: Vec<f64>,
    cell_opacity_scale: Vec<f64>,
    cell_emission: Vec<Dd>,
    cell_cross_heat: Vec<Dd>,
    fine_point: Vec<Dd>,
    angular_point: Vec<Dd>,
    fine_panel: Vec<Dd>,
    coarse_panel: Vec<Dd>,
    angular_panel: Vec<Dd>,
    values: Vec<Dd>,
    spectral_errors: Vec<Dd>,
    angular_errors: Vec<Dd>,
    tail: Vec<f64>,
}

impl Scratch {
    fn new(channels: usize, cells: usize) -> Result<Self, Error> {
        let count = CHANNEL
            .checked_mul(channels)
            .and_then(|n| CELL.checked_mul(cells).and_then(|m| n.checked_add(m)))
            .and_then(|n| FACE.checked_mul(cells + 1).and_then(|m| n.checked_add(m)))
            .ok_or(Error::WorkBudget)?;
        Ok(Self {
            opacity: buffer(channels, 0.0)?,
            fractions: buffer(channels, Dd::ZERO)?,
            emission: buffer(channels, 0.0)?,
            cross: buffer(channels, Dd::ZERO)?,
            cell_opacity: buffer(cells, 0.0)?,
            cell_opacity_scale: buffer(cells, 0.0)?,
            cell_emission: buffer(cells, Dd::ZERO)?,
            cell_cross_heat: buffer(cells, Dd::ZERO)?,
            fine_point: buffer(count, Dd::ZERO)?,
            angular_point: buffer(count, Dd::ZERO)?,
            fine_panel: buffer(count, Dd::ZERO)?,
            coarse_panel: buffer(count, Dd::ZERO)?,
            angular_panel: buffer(count, Dd::ZERO)?,
            values: buffer(count, Dd::ZERO)?,
            spectral_errors: buffer(count, Dd::ZERO)?,
            angular_errors: buffer(count, Dd::ZERO)?,
            tail: buffer(count, 0.0)?,
        })
    }
}

/// One common field for all original channels in each spatial slab. Beam and
/// diffuse quadratures remain separately controlled; persistent buffers avoid
/// state/bank cloning and allocation inside spectral or angular evaluations.
pub(in crate::chemistry::ionization_network::transport) struct PreparedDiffuseField<'a> {
    prepared: Prepared<'a>,
    beam: PreparedField,
    original_beams: Vec<DiffuseRecombinationBeam>,
    scratch: Scratch,
    grid: Vec<f64>,
    fine_angles: Vec<(f64, f64)>,
    coarse_angles: Vec<(f64, f64)>,
    angular: usize,
    rates: Vec<NetworkDiffuseChannelRate>,
    cell_balances: Vec<DiffuseRecombinationCellBalance>,
    net: Vec<[Dd; 3]>,
    preparation: usize,
    used: usize,
    beam_query_work: usize,
    pub beam_rates: Vec<NetworkRadiationBeamRate>,
    pub interfaces: Vec<DiffuseRecombinationFlux>,
    pub balance: DiffuseRecombinationBalance,
}

impl<'a> PreparedDiffuseField<'a> {
    pub fn new(
        cells: &[NetworkRadiationCell<'a>],
        problems: &[Problem<'a>],
        area: Area,
        beams: &[DiffuseRecombinationBeam],
        controls: NetworkDiffuseSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, NetworkRadiationError> {
        let prepared =
            Prepared::new(cells, area, controls.radiation, cancelled).map_err(diffuse_error)?;
        let beam = PreparedField::new(cells, problems, area, beams, controls.network, cancelled)?;
        let channels = prepared.channels.len();
        let preparation = prepared
            .work
            .checked_add(beam.work)
            .ok_or(NetworkRadiationError::WorkBudget)?;
        let beam_query_work = channels
            .checked_add(
                beams
                    .len()
                    .checked_mul(channels + cells.len() + 1)
                    .ok_or(NetworkRadiationError::WorkBudget)?,
            )
            .ok_or(NetworkRadiationError::WorkBudget)?;
        let mut original_beams = Vec::new();
        original_beams
            .try_reserve_exact(beams.len())
            .map_err(|_| diffuse_error(Error::Allocation))?;
        original_beams.extend_from_slice(beams);
        Ok(Self {
            scratch: Scratch::new(channels, cells.len()).map_err(diffuse_error)?,
            rates: buffer(channels, NetworkDiffuseChannelRate::default()).map_err(diffuse_error)?,
            cell_balances: buffer(cells.len(), DiffuseRecombinationCellBalance::default())
                .map_err(diffuse_error)?,
            net: buffer(cells.len(), [Dd::ZERO; 3]).map_err(diffuse_error)?,
            interfaces: buffer(cells.len() + 1, DiffuseRecombinationFlux::default())
                .map_err(diffuse_error)?,
            beam_rates: buffer(beams.len(), NetworkRadiationBeamRate::default())
                .map_err(diffuse_error)?,
            prepared,
            beam,
            original_beams,
            preparation,
            beam_query_work,
            grid: Vec::new(),
            fine_angles: Vec::new(),
            coarse_angles: Vec::new(),
            angular: 1,
            used: 0,
            balance: DiffuseRecombinationBalance::default(),
        })
    }

    pub fn preparation_work(&self) -> usize {
        self.preparation
    }
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn illuminated(&self) -> bool {
        self.beam.illuminated()
    }
    pub fn channel_rates(&self, cell: usize) -> &[NetworkDiffuseChannelRate] {
        &self.rates[self.prepared.offsets[cell]..self.prepared.offsets[cell + 1]]
    }
    pub fn cell_beam_rates(&self, cell: usize) -> &[NetworkRadiationBeamRate] {
        self.beam.cell_beam_rates(cell)
    }
    pub fn cell_net_rates(&self, cell: usize) -> [Dd; 3] {
        self.net[cell]
    }
    pub fn cell_balance(&self, cell: usize) -> DiffuseRecombinationCellBalance {
        self.cell_balances[cell]
    }

    pub fn evaluate(
        &mut self,
        points: &[Point],
        remaining_work: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), NetworkRadiationError> {
        self.used = 0;
        self.prepared.work = 0;
        self.prepared.panels = 0;
        if points.len() + 1 != self.prepared.offsets.len() {
            return Err(NetworkRadiationError::Invalid(
                "diffuse network state shape",
            ));
        }
        let per_query = self.prepared.controls.maximum_work;
        let limit = per_query.min(remaining_work);
        if self.beam_query_work > limit {
            return Err(NetworkRadiationError::WorkBudget);
        }
        let before = self.beam.work;
        let beam_result = self.beam.evaluate(points, cancelled);
        self.used = self.beam.work.saturating_sub(before);
        beam_result?;
        self.prepared.controls.maximum_work = limit.saturating_sub(self.used);
        let result = (|| {
            self.prepared.refresh(points, cancelled)?;
            let normalized = spectral::run(
                &mut self.prepared,
                &mut self.scratch,
                &mut self.grid,
                &mut self.angular,
                &mut self.fine_angles,
                &mut self.coarse_angles,
                cancelled,
            )?;
            self.publish(normalized, cancelled)
        })();
        self.used = self
            .used
            .checked_add(self.prepared.work)
            .ok_or(NetworkRadiationError::WorkBudget)?;
        self.prepared.controls.maximum_work = per_query;
        result.map_err(diffuse_error)?;
        self.balance.work = self.used;
        self.beam_rates.copy_from_slice(&self.beam.beam_rates);
        Ok(())
    }
}

impl<'a> Prepared<'a> {
    fn new(
        cells: &[NetworkRadiationCell<'a>],
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
                "network geometry, source or quadrature controls",
            ));
        }
        let count = cells
            .iter()
            .try_fold(0_usize, |n, cell| n.checked_add(cell.sources.len()))
            .ok_or(Error::WorkBudget)?;
        if count > controls.maximum_work {
            return Err(Error::WorkBudget);
        }
        let mut channels = Vec::new();
        channels
            .try_reserve_exact(count)
            .map_err(|_| Error::Allocation)?;
        let mut offsets = Vec::new();
        offsets
            .try_reserve_exact(cells.len() + 1)
            .map_err(|_| Error::Allocation)?;
        offsets.push(0);
        for cell in cells {
            poll(cancelled)?;
            let mut sources: Vec<_> = Vec::new();
            sources
                .try_reserve_exact(cell.sources.len())
                .map_err(|_| Error::Allocation)?;
            sources.extend(cell.sources.iter());
            sources.sort_by_key(|source| {
                let channel = source.source().channel();
                (channel.lower_index(), channel.upper_index())
            });
            for source in sources {
                poll(cancelled)?;
                let index = channels.len();
                let state = cell.initial;
                let channel = source.source().channel();
                if state.thermodynamics() != channel.thermodynamics()
                    || state.electron_index() != Some(channel.electron_index())
                    || state.relative_charge_residual() > controls.relative_charge_tolerance
                {
                    return Err(Error::Invalid(
                        "network state and original reciprocal source bank",
                    ));
                }
                let spectrum = source
                    .spectrum(state.electron_temperature())
                    .map_err(|error| source_error(index, error))?;
                channels.push(Binding {
                    source,
                    spectrum,
                    lower: channel.lower_index(),
                    upper: channel.upper_index(),
                    electron: channel.electron_index(),
                    opacity_per_sigma: 0.0,
                    pair_density_volume: 0.0,
                    gap: channel.ionization_energy_j_mol() / NA,
                    expected: [0.0; 3],
                    cutoff: 0.0,
                    tail: [0.0; 3],
                    volume: state.volume().value(),
                });
            }
            offsets.push(channels.len());
        }
        Ok(Self {
            channels,
            offsets,
            area: area.value(),
            shared_upper: f64::INFINITY,
            controls,
            panels: 0,
            work: count,
        })
    }

    fn refresh(
        &mut self,
        points: &[Point],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.shared_upper = f64::INFINITY;
        for (c, point) in points.iter().enumerate() {
            if !positive(point.electron_temperature())
                || point.amounts.iter().any(|v| !nonnegative(*v))
            {
                return Err(Error::Invalid(
                    "network radiation populations or electron temperature",
                ));
            }
            for index in self.offsets[c]..self.offsets[c + 1] {
                poll(cancelled)?;
                self.spend(1)?;
                let binding = &mut self.channels[index];
                let temperature = Temperature::new(point.electron_temperature());
                binding.spectrum = binding
                    .source
                    .spectrum(temperature)
                    .map_err(|error| source_error(index, error))?;
                let lower = point.amounts[binding.lower];
                let upper = point.amounts[binding.upper];
                let electron = point.amounts[binding.electron];
                binding.opacity_per_sigma = NA * (lower / self.area);
                binding.pair_density_volume = (NA * upper) * (NA * (electron / binding.volume));
                binding.cutoff = 0.0;
                binding.tail = [0.0; 3];
                binding.expected = [0.0; 3];
                if !nonnegative(binding.opacity_per_sigma)
                    || !nonnegative(binding.pair_density_volume)
                    || (lower > 0.0 && binding.opacity_per_sigma == 0.0)
                    || (upper > 0.0 && electron > 0.0 && binding.pair_density_volume == 0.0)
                {
                    return Err(Error::Arithmetic("network refreshed source inventory"));
                }
                if binding.pair_density_volume > 0.0 {
                    let moments = binding
                        .source
                        .evaluate(
                            temperature,
                            self.controls.source_quadrature,
                            &mut *cancelled,
                        )
                        .map_err(|error| source_error(index, error))?;
                    binding.expected = [
                        moments.coefficient_m3_s,
                        moments.emitted_photon_energy_j_m3_s,
                        moments.captured_kinetic_energy_j_m3_s,
                    ]
                    .map(|v| v * binding.pair_density_volume);
                    if binding.expected.iter().any(|v| !nonnegative(*v))
                        || (moments.coefficient_m3_s > 0.0 && binding.expected[0] == 0.0)
                    {
                        return Err(Error::Arithmetic("network original emission moments"));
                    }
                }
                if binding.opacity_per_sigma > 0.0 || binding.expected[0] > 0.0 {
                    self.shared_upper = self
                        .shared_upper
                        .min(binding.source.source().maximum_energy().value());
                }
            }
        }
        Ok(())
    }

    fn spend(&mut self, work: usize) -> Result<(), Error> {
        self.work = self
            .work
            .checked_add(work)
            .filter(|v| *v <= self.controls.maximum_work)
            .ok_or(Error::WorkBudget)?;
        Ok(())
    }
    fn cells(&self) -> usize {
        self.offsets.len() - 1
    }
    fn cell_slot(&self, c: usize) -> usize {
        CHANNEL * self.channels.len() + CELL * c
    }
    fn face_slot(&self, face: usize) -> usize {
        CHANNEL * self.channels.len() + CELL * self.cells() + FACE * face
    }
    fn allowance(&self, index: usize, value: f64) -> f64 {
        let channels = CHANNEL * self.channels.len();
        let faces = channels + CELL * self.cells();
        let number = if index < channels {
            matches!(index % CHANNEL, 0 | 3 | 4 | 5 | 8)
        } else if index < faces {
            (index - channels) % CELL == 0
        } else {
            (index - faces) % 2 == 0
        };
        (if number {
            self.controls.absolute_photon_rate_tolerance_per_s
        } else {
            self.controls.absolute_power_tolerance_w
        }) + self.controls.relative_tolerance * value.abs()
    }
    fn check_balance(
        &self,
        photons: f64,
        power: f64,
        photon_scale: f64,
        power_scale: f64,
    ) -> Result<(), Error> {
        if !photons.is_finite()
            || !power.is_finite()
            || photons.abs()
                > self.controls.absolute_photon_rate_tolerance_per_s
                    + self.controls.relative_tolerance * photon_scale
            || power.abs()
                > self.controls.absolute_power_tolerance_w
                    + self.controls.relative_tolerance * power_scale
        {
            Err(Error::Balance {
                photon_residual_per_s: photons,
                energy_residual_w: power,
            })
        } else {
            Ok(())
        }
    }
}

fn buffer<T: Clone>(n: usize, value: T) -> Result<Vec<T>, Error> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| Error::Allocation)?;
    values.resize(n, value);
    Ok(values)
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
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

impl PreparedDiffuseField<'_> {
    fn publish(
        &mut self,
        normalized: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let n = self.prepared.cells();
        self.prepared
            .spend(self.prepared.channels.len() + n + self.interfaces.len())?;
        // Beam characteristics are evaluated once by the existing shared-opacity
        // kernel. Insert their face traffic only AFTER diffuse error normalization.
        for (b, beam) in self.original_beams.iter().enumerate() {
            poll(cancelled)?;
            self.prepared.spend(n + 1)?;
            let right = beam.direction_cosine > 0.0;
            let direction = if right { 0 } else { 2 };
            let entry = if right { 0 } else { n };
            let energy = beam.photons.photon_energy.value();
            let entry_slot = self.prepared.face_slot(entry) + direction;
            let incoming = self.beam.beam_rates[b].incident_photons_per_s;
            self.scratch.values[entry_slot] =
                self.scratch.values[entry_slot] + Dd::from_f64(incoming);
            self.scratch.values[entry_slot + 1] =
                self.scratch.values[entry_slot + 1] + Dd::from_f64(incoming * energy);
            for c in 0..n {
                let face = if right { c + 1 } else { c };
                let slot = self.prepared.face_slot(face) + direction;
                let rate = self.beam.cell_beam_rates(c)[b].transmitted_photons_per_s;
                self.scratch.values[slot] = self.scratch.values[slot] + Dd::from_f64(rate);
                self.scratch.values[slot + 1] =
                    self.scratch.values[slot + 1] + Dd::from_f64(rate * energy);
            }
        }
        for face in 0..=n {
            let slot = self.prepared.face_slot(face);
            let v = &self.scratch.values[slot..slot + FACE];
            if v.iter().any(|v| !finite(*v) || v.to_f64() < 0.0) {
                return Err(Error::Arithmetic("network interface flux"));
            }
            self.interfaces[face] = DiffuseRecombinationFlux {
                right_photons_per_s: v[0].to_f64(),
                right_power_w: v[1].to_f64(),
                left_photons_per_s: v[2].to_f64(),
                left_power_w: v[3].to_f64(),
            };
        }
        // Emission number/power, total absorption number/power, omitted tails.
        let mut sums = [Dd::ZERO; 6];
        let mut direct_column = [Dd::ZERO; 2];
        for c in 0..n {
            poll(cancelled)?;
            // emitted N/E/captured E, same-cell N, nonlocal N, absorbed E,
            // photoelectron E, incident N/E, and original chemical source.
            let mut aggregate = [Dd::ZERO; 10];
            let cell_slot = self.prepared.cell_slot(c);
            let mut net = [
                self.scratch.values[cell_slot] / Dd::from_f64(NA),
                self.scratch.values[cell_slot + 1],
                self.scratch.values[cell_slot + 2],
            ];
            for (local, k) in (self.prepared.offsets[c]..self.prepared.offsets[c + 1]).enumerate() {
                let binding = &self.prepared.channels[k];
                let v = &self.scratch.values[CHANNEL * k..CHANNEL * (k + 1)];
                let beam = self.beam.channel_rates(c)[local];
                let incident = Dd::from_f64(beam.ionization_rate_mol_per_s) * Dd::from_f64(NA);
                let absorbed = v[6] + Dd::from_f64(beam.absorbed_power_w);
                let heat = v[7] + Dd::from_f64(beam.electron_heating_power_w);
                let extent = v[8] / Dd::from_f64(NA) + Dd::from_f64(beam.ionization_rate_mol_per_s);
                let electron = v[9] + Dd::from_f64(beam.electron_heating_power_w);
                let gas = v[10] + Dd::from_f64(beam.absorbed_power_w);
                let chemical = (v[8] + incident) * Dd::from_f64(binding.gap);
                let all_absorbed = v[3] + v[4] + v[5] + incident;
                let gross_number = (v[0] + extent * Dd::from_f64(NA) - all_absorbed).to_f64();
                let gross_energy = (v[1] + gas - absorbed).to_f64();
                self.prepared.check_balance(
                    gross_number,
                    gross_energy,
                    (v[0] + all_absorbed).to_f64(),
                    (v[1] + absorbed).to_f64(),
                )?;
                let closure = (chemical + electron - gas).to_f64();
                let allowance = self.prepared.controls.absolute_power_tolerance_w
                    + self.prepared.controls.relative_tolerance
                        * chemical
                            .abs()
                            .to_f64()
                            .max(electron.abs().to_f64())
                            .max(gas.abs().to_f64());
                if !closure.is_finite() || !nonnegative(allowance) || closure.abs() > allowance {
                    return Err(Error::Balance {
                        photon_residual_per_s: 0.0,
                        energy_residual_w: closure,
                    });
                }
                self.rates[k] = NetworkDiffuseChannelRate {
                    emitted_photons_per_s: v[0].to_f64(),
                    emitted_power_w: v[1].to_f64(),
                    captured_electron_power_w: v[2].to_f64(),
                    self_absorbed_photons_per_s: v[3].to_f64(),
                    cross_channel_absorbed_photons_per_s: v[4].to_f64(),
                    nonlocal_absorbed_photons_per_s: v[5].to_f64(),
                    incident_absorbed_photons_per_s: incident.to_f64(),
                    incident_absorbed_power_w: beam.absorbed_power_w,
                    diffuse_absorbed_power_w: v[6].to_f64(),
                    absorbed_power_w: absorbed.to_f64(),
                    photoelectron_power_w: heat.to_f64(),
                    net_ionization_rate_mol_per_s: extent.to_f64(),
                    net_electron_power_w: electron.to_f64(),
                    chemical_storage_rate_w: chemical.to_f64(),
                    net_gas_power_w: gas.to_f64(),
                    emission_cutoff_j: binding.cutoff,
                };
                for (sum, value) in aggregate.iter_mut().zip([
                    v[0],
                    v[1],
                    v[2],
                    v[3] + v[4],
                    v[5],
                    absorbed,
                    heat,
                    incident,
                    Dd::from_f64(beam.absorbed_power_w),
                    chemical,
                ]) {
                    *sum = *sum + value;
                }
                net[1] = net[1] + Dd::from_f64(beam.electron_heating_power_w);
                sums[4] = sums[4] + Dd::from_f64(binding.tail[0]);
                sums[5] = sums[5] + Dd::from_f64(binding.tail[1]);
            }
            // Use the beam kernel's independent cell photon/power totals, not a
            // rounded sum of channel products, for the shared direct source.
            for beam in self.beam.cell_beam_rates(c) {
                net[0] = net[0] + Dd::from_f64(beam.absorbed_photons_per_s) / Dd::from_f64(NA);
                net[2] = net[2] + Dd::from_f64(beam.absorbed_power_w);
            }
            self.net[c] = net;
            let left = self.interfaces[c];
            let right = self.interfaces[c + 1];
            let number = net[0] * Dd::from_f64(NA);
            let photon_flux_residual = (number
                - Dd::from_f64(left.right_photons_per_s)
                - Dd::from_f64(right.left_photons_per_s)
                + Dd::from_f64(right.right_photons_per_s)
                + Dd::from_f64(left.left_photons_per_s))
            .to_f64();
            let power_flux_residual =
                (net[2] - Dd::from_f64(left.right_power_w) - Dd::from_f64(right.left_power_w)
                    + Dd::from_f64(right.right_power_w)
                    + Dd::from_f64(left.left_power_w))
                .to_f64();
            let photon_traffic = left.right_photons_per_s
                + right.left_photons_per_s
                + right.right_photons_per_s
                + left.left_photons_per_s;
            let power_traffic =
                left.right_power_w + right.left_power_w + right.right_power_w + left.left_power_w;
            let number_allowance = self.prepared.allowance(cell_slot, number.to_f64())
                + 128.0 * f64::EPSILON * photon_traffic;
            let power_allowance = self.prepared.allowance(cell_slot + 2, net[2].to_f64())
                + 128.0 * f64::EPSILON * power_traffic;
            if !photon_flux_residual.is_finite()
                || !power_flux_residual.is_finite()
                || !nonnegative(number_allowance)
                || !nonnegative(power_allowance)
                || photon_flux_residual.abs() > number_allowance
                || power_flux_residual.abs() > power_allowance
            {
                return Err(Error::Balance {
                    photon_residual_per_s: photon_flux_residual,
                    energy_residual_w: power_flux_residual,
                });
            }
            let total_absorbed = aggregate[3] + aggregate[4] + aggregate[7];
            let gross_photon = (aggregate[0]
                + Dd::from_f64(left.right_photons_per_s)
                + Dd::from_f64(right.left_photons_per_s)
                - total_absorbed
                - Dd::from_f64(right.right_photons_per_s)
                - Dd::from_f64(left.left_photons_per_s))
            .to_f64();
            let gross_power = (aggregate[1]
                + Dd::from_f64(left.right_power_w)
                + Dd::from_f64(right.left_power_w)
                - aggregate[5]
                - Dd::from_f64(right.right_power_w)
                - Dd::from_f64(left.left_power_w))
            .to_f64();
            self.prepared.check_balance(
                gross_photon,
                gross_power,
                aggregate[0].to_f64() + photon_traffic,
                aggregate[1].to_f64() + power_traffic,
            )?;
            let gas_residual = (aggregate[9] + net[1] - net[2]).to_f64();
            let gas_allowance = self.prepared.controls.absolute_power_tolerance_w
                + 128.0
                    * f64::EPSILON
                    * (aggregate[9].abs().to_f64() + net[1].abs().to_f64() + net[2].abs().to_f64());
            if !gas_residual.is_finite()
                || !nonnegative(gas_allowance)
                || gas_residual.abs() > gas_allowance
            {
                return Err(Error::Balance {
                    photon_residual_per_s: 0.0,
                    energy_residual_w: gas_residual,
                });
            }
            self.cell_balances[c] = DiffuseRecombinationCellBalance {
                emitted_photons_per_s: aggregate[0].to_f64(),
                emitted_power_w: aggregate[1].to_f64(),
                captured_electron_power_w: aggregate[2].to_f64(),
                self_absorbed_photons_per_s: aggregate[3].to_f64(),
                nonlocal_absorbed_photons_per_s: aggregate[4].to_f64(),
                incident_absorbed_photons_per_s: aggregate[7].to_f64(),
                incident_absorbed_power_w: aggregate[8].to_f64(),
                absorbed_power_w: aggregate[5].to_f64(),
                photoelectron_power_w: aggregate[6].to_f64(),
                net_ionization_rate_mol_per_s: net[0].to_f64(),
                net_electron_power_w: net[1].to_f64(),
                chemical_storage_rate_w: aggregate[9].to_f64(),
                net_gas_power_w: net[2].to_f64(),
                photon_residual_per_s: gross_photon,
                energy_residual_w: gross_power,
                gas_energy_residual_w: gas_residual,
                net_photon_flux_residual_per_s: photon_flux_residual,
                net_power_flux_residual_w: power_flux_residual,
            };
            for (sum, value) in
                sums[..4]
                    .iter_mut()
                    .zip([aggregate[0], aggregate[1], total_absorbed, aggregate[5]])
            {
                *sum = *sum + value;
            }
            direct_column[0] = direct_column[0] + number;
            direct_column[1] = direct_column[1] + net[2];
        }
        let escaped = Dd::from_f64(self.interfaces[0].left_photons_per_s)
            + Dd::from_f64(self.interfaces[n].right_photons_per_s);
        let escaped_power = Dd::from_f64(self.interfaces[0].left_power_w)
            + Dd::from_f64(self.interfaces[n].right_power_w);
        let incident = Dd::from_f64(self.interfaces[0].right_photons_per_s)
            + Dd::from_f64(self.interfaces[n].left_photons_per_s);
        let incident_power = Dd::from_f64(self.interfaces[0].right_power_w)
            + Dd::from_f64(self.interfaces[n].left_power_w);
        let photon_residual = (sums[0] + incident - sums[2] - escaped).to_f64();
        let energy_residual = (sums[1] + incident_power - sums[3] - escaped_power).to_f64();
        self.prepared.check_balance(
            photon_residual,
            energy_residual,
            (sums[0] + incident).to_f64(),
            (sums[1] + incident_power).to_f64(),
        )?;
        let direct_photon = (direct_column[0] + escaped - incident).to_f64();
        let direct_power = (direct_column[1] + escaped_power - incident_power).to_f64();
        let photon_allowance = self.prepared.controls.absolute_photon_rate_tolerance_per_s
            + self.prepared.controls.relative_tolerance * direct_column[0].abs().to_f64()
            + 128.0 * f64::EPSILON * (escaped + incident).to_f64();
        let power_allowance = self.prepared.controls.absolute_power_tolerance_w
            + self.prepared.controls.relative_tolerance * direct_column[1].abs().to_f64()
            + 128.0 * f64::EPSILON * (escaped_power + incident_power).to_f64();
        if !direct_photon.is_finite()
            || !direct_power.is_finite()
            || !nonnegative(photon_allowance)
            || !nonnegative(power_allowance)
            || direct_photon.abs() > photon_allowance
            || direct_power.abs() > power_allowance
        {
            return Err(Error::Balance {
                photon_residual_per_s: direct_photon,
                energy_residual_w: direct_power,
            });
        }
        self.balance = DiffuseRecombinationBalance {
            incident_photons_per_s: incident.to_f64(),
            incident_power_w: incident_power.to_f64(),
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
            spectral_panels: self.prepared.panels,
            angular_subdivisions: self.angular,
            work: self.prepared.work,
        };
        Ok(())
    }
}
