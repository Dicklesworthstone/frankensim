use super::*;
use fs_math::dd::Dd;
use fs_math::det::{exp, expm1};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;

/// Immutable source table and reused rate buffers; no source queries or
/// population allocation occur in a radiation stage.
pub(super) struct PreparedField {
    offsets: Vec<usize>,
    lower: Vec<usize>,
    gaps: Vec<f64>,
    // Beam-major normal optical depth per mole of the corresponding lower stage.
    opacity: Vec<f64>,
    beams: Vec<DiffuseRecombinationBeam>,
    incoming: Vec<f64>,
    rates: Vec<NetworkPhotoionizationRate>,
    sums: Vec<[Dd; 3]>,
    cell_beams: Vec<NetworkRadiationBeamRate>,
    weights: Vec<f64>,
    pub beam_rates: Vec<NetworkRadiationBeamRate>,
    pub work: usize,
    pub evaluations: usize,
    controls: NetworkRadiationSolve,
}
impl PreparedField {
    pub fn new(
        cells: &[NetworkRadiationCell<'_>],
        problems: &[Problem<'_>],
        area: Area,
        beams: &[DiffuseRecombinationBeam],
        controls: NetworkRadiationSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if !positive(area.value()) {
            return Err(Error::Invalid("beam face area"));
        }
        let mut offsets = reserved(cells.len() + 1)?;
        let count: usize = cells.iter().map(|cell| cell.sources.len()).sum();
        let mut lower = reserved(count)?;
        let mut gaps = reserved(count)?;
        offsets.push(0);
        for problem in problems {
            for source in problem.sources() {
                poll(cancelled)?;
                lower.push(source.source().channel().lower_index());
                gaps.push(source.source().channel().ionization_energy_j_mol() / NA);
            }
            offsets.push(lower.len());
        }
        let count = lower.len();
        let mut result = Self {
            offsets,
            lower,
            gaps,
            opacity: reserved(count * beams.len())?,
            beams: beams.to_vec(),
            incoming: reserved(beams.len())?,
            rates: buffer(count, NetworkPhotoionizationRate::default())?,
            sums: buffer(count, [Dd::ZERO; 3])?,
            cell_beams: buffer(
                cells.len() * beams.len(),
                NetworkRadiationBeamRate::default(),
            )?,
            weights: buffer(count, 0.0)?,
            beam_rates: buffer(beams.len(), NetworkRadiationBeamRate::default())?,
            work: 0,
            evaluations: 0,
            controls,
        };
        for (b, beam) in beams.iter().enumerate() {
            poll(cancelled)?;
            let mu = beam.direction_cosine.abs();
            if !positive(mu)
                || mu > 1.0
                || !positive(beam.photons.photon_energy.value())
                || !nonnegative(beam.photons.number_flux_per_m2_s)
            {
                return Err(Error::Invalid("beam energy, number flux or direction"));
            }
            let incoming = (area.value() * mu) * beam.photons.number_flux_per_m2_s;
            if !nonnegative(incoming)
                || (beam.photons.number_flux_per_m2_s > 0.0 && incoming == 0.0)
                || !nonnegative(incoming * beam.photons.photon_energy.value())
            {
                return Err(Error::Arithmetic("boundary photon rate"));
            }
            result.incoming.push(incoming);
            for (c, problem) in problems.iter().enumerate() {
                for (channel, source) in problem.sources().enumerate() {
                    poll(cancelled)?;
                    result.spend(1)?;
                    let sigma = source
                        .source()
                        .cross_section(beam.photons.photon_energy)
                        .map_err(|error| Error::Source {
                            cell: c,
                            channel,
                            beam: b,
                            error,
                        })?
                        .value();
                    let opacity = (sigma * NA / area.value()) / mu;
                    if !nonnegative(opacity) || (sigma > 0.0 && opacity == 0.0) {
                        return Err(Error::Arithmetic("channel optical depth coefficient"));
                    }
                    result.opacity.push(opacity);
                }
            }
        }
        Ok(result)
    }
    pub fn channel_rates(&self, cell: usize) -> &[NetworkPhotoionizationRate] {
        &self.rates[self.offsets[cell]..self.offsets[cell + 1]]
    }
    pub fn cell_beam_rates(&self, cell: usize) -> &[NetworkRadiationBeamRate] {
        let n = self.beams.len();
        &self.cell_beams[cell * n..(cell + 1) * n]
    }
    pub fn illuminated(&self) -> bool {
        self.incoming.iter().any(|r| *r > 0.0)
    }
    pub fn evaluate(
        &mut self,
        points: &[Point],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let cells = self.offsets.len() - 1;
        if points.len() != cells {
            return Err(Error::Invalid("radiation state shape"));
        }
        if self.evaluations >= self.controls.maximum_radiation_evaluations {
            return Err(Error::WorkBudget);
        }
        self.evaluations += 1;
        self.spend(self.lower.len())?;
        self.sums.fill([Dd::ZERO; 3]);
        let channels = self.lower.len();
        for b in 0..self.beams.len() {
            poll(cancelled)?;
            self.spend(1)?;
            let beam = self.beams[b];
            let energy = beam.photons.photon_energy.value();
            let mut incoming = self.incoming[b];
            let mut absorbed_total = Dd::ZERO;
            for ordinal in 0..cells {
                poll(cancelled)?;
                let c = if beam.direction_cosine > 0.0 {
                    ordinal
                } else {
                    cells - 1 - ordinal
                };
                let range = self.offsets[c]..self.offsets[c + 1];
                self.spend(range.len() + 1)?;
                let mut maximum = 0.0_f64;
                for k in range.clone() {
                    let amount = points[c].amounts[self.lower[k]];
                    let depth = self.opacity[b * channels + k] * amount;
                    if !nonnegative(depth)
                        || (amount > 0.0 && self.opacity[b * channels + k] > 0.0 && depth == 0.0)
                    {
                        return Err(Error::Arithmetic("channel optical depth"));
                    }
                    self.weights[k] = depth;
                    maximum = maximum.max(depth);
                }
                let mut scaled = Dd::ZERO;
                if maximum > 0.0 {
                    for k in range.clone() {
                        scaled = scaled + Dd::from_f64(self.weights[k] / maximum);
                    }
                }
                let scale = scaled.to_f64();
                let depth = maximum * scale;
                let absorbed = if depth.is_infinite() {
                    incoming
                } else {
                    incoming * -expm1(-depth)
                };
                let outgoing = if depth.is_infinite() {
                    0.0
                } else {
                    incoming * exp(-depth)
                };
                if !nonnegative(absorbed) || !nonnegative(outgoing) {
                    return Err(Error::Arithmetic("beam characteristic"));
                }
                let mut allocated = Dd::ZERO;
                for k in range {
                    let fraction = if maximum == 0.0 {
                        0.0
                    } else {
                        (self.weights[k] / maximum) / scale
                    };
                    let photons = absorbed * fraction;
                    let heat = if photons == 0.0 {
                        0.0
                    } else {
                        photons * (energy - self.gaps[k])
                    };
                    let power = photons * energy;
                    let mol = photons / NA;
                    if !nonnegative(photons)
                        || !nonnegative(heat)
                        || !nonnegative(power)
                        || !nonnegative(mol)
                        || (photons > 0.0 && mol == 0.0)
                    {
                        return Err(Error::Arithmetic("allocated photoionization rate"));
                    }
                    allocated = allocated + Dd::from_f64(photons);
                    for (sum, value) in self.sums[k].iter_mut().zip([mol, power, heat]) {
                        *sum = *sum + Dd::from_f64(value);
                    }
                }
                let allocated = allocated.to_f64();
                let rate = self.audit(Some(c), b, incoming, allocated, outgoing)?;
                self.cell_beams[c * self.beams.len() + b] = rate;
                absorbed_total = absorbed_total + Dd::from_f64(allocated);
                incoming = outgoing;
            }
            self.beam_rates[b] =
                self.audit(None, b, self.incoming[b], absorbed_total.to_f64(), incoming)?;
        }
        for (rate, sums) in self.rates.iter_mut().zip(&self.sums) {
            let values = sums.map(Dd::to_f64);
            if values.iter().any(|x| !nonnegative(*x)) {
                return Err(Error::Arithmetic("summed photon rates"));
            }
            *rate = NetworkPhotoionizationRate {
                ionization_rate_mol_per_s: values[0],
                absorbed_power_w: values[1],
                electron_heating_power_w: values[2],
            };
        }
        Ok(())
    }
    fn audit(
        &self,
        cell: Option<usize>,
        beam: usize,
        incoming: f64,
        absorbed: f64,
        outgoing: f64,
    ) -> Result<NetworkRadiationBeamRate, Error> {
        let residual =
            (Dd::from_f64(absorbed) + Dd::from_f64(outgoing) - Dd::from_f64(incoming)).to_f64();
        if !residual.is_finite()
            || residual.abs() > self.controls.relative_photon_tolerance * incoming
        {
            return Err(Error::PhotonBalance {
                cell,
                beam,
                residual_photons: residual,
            });
        }
        let power = absorbed * self.beams[beam].photons.photon_energy.value();
        if !nonnegative(power) {
            return Err(Error::Arithmetic("absorbed photon power"));
        }
        Ok(NetworkRadiationBeamRate {
            incident_photons_per_s: incoming,
            absorbed_photons_per_s: absorbed,
            transmitted_photons_per_s: outgoing,
            absorbed_power_w: power,
            photon_residual_per_s: residual,
        })
    }
    fn spend(&mut self, work: usize) -> Result<(), Error> {
        self.work = self.work.checked_add(work).ok_or(Error::WorkBudget)?;
        if self.work > self.controls.maximum_radiation_work {
            return Err(Error::WorkBudget);
        }
        Ok(())
    }
}
