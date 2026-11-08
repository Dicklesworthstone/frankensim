use super::*;
use fs_math::det::{exp, expm1};

// Positive halves of the same Gauss-Legendre rules as the original transfer
// owner; angular integration uses mu=t^4 to resolve grazing characteristics.
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

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    prepared: &mut Prepared<'_>,
    scratch: &mut Scratch,
    grid: &mut Vec<f64>,
    angular: &mut usize,
    fine_angles: &mut Vec<(f64, f64)>,
    coarse_angles: &mut Vec<(f64, f64)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], Error> {
    scratch.values.fill(Dd::ZERO);
    scratch.spectral_errors.fill(Dd::ZERO);
    scratch.angular_errors.fill(Dd::ZERO);
    scratch.tail.fill(0.0);
    if !prepared.channels.iter().any(|c| c.expected[0] > 0.0) {
        return Ok([0.0; 3]);
    }
    let (mut extent, mut refinement) = (32.0, 1_usize);
    loop {
        poll(cancelled)?;
        prepared.cutoffs(extent)?;
        prepared.grid(grid, cancelled)?;
        angles(*angular, &X16, &W16, fine_angles)?;
        angles(*angular, &X8, &W8, coarse_angles)?;
        prepared.integrate(
            grid,
            refinement,
            fine_angles,
            coarse_angles,
            scratch,
            cancelled,
        )?;
        prepared.tail_vector(&mut scratch.tail);
        let normalized = prepared.errors(scratch)?;
        if normalized[2] > 0.25 {
            if extent >= 1024.0
                || prepared
                    .channels
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
            if *angular >= prepared.controls.maximum_angular_subdivisions {
                return Err(Error::AngularBudget);
            }
            *angular = angular
                .saturating_mul(2)
                .min(prepared.controls.maximum_angular_subdivisions);
            continue;
        }
        return Ok(normalized);
    }
}

impl Prepared<'_> {
    fn cutoffs(&mut self, thermal_extent: f64) -> Result<(), Error> {
        for (index, source) in self.channels.iter_mut().enumerate() {
            if source.expected[0] == 0.0 {
                continue;
            }
            let end = source.source.source().threshold().value()
                + thermal_extent * source.spectrum.thermal_energy().value();
            if !positive(end) {
                return Err(Error::Arithmetic("network spectral cutoff"));
            }
            source.cutoff = end.min(self.shared_upper);
            let tail = source
                .spectrum
                .tail_bound(Energy::new(source.cutoff))
                .map_err(|error| source_error(index, error))?;
            source.tail = [
                tail.coefficient_m3_s,
                tail.emitted_photon_energy_j_m3_s,
                tail.captured_kinetic_energy_j_m3_s,
            ]
            .map(|v| (v * source.pair_density_volume).max(f64::from_bits(1)));
            if source.tail.iter().any(|v| !positive(*v)) {
                return Err(Error::Arithmetic("network emission tail scaling"));
            }
        }
        Ok(())
    }

    fn grid(&self, grid: &mut Vec<f64>, cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
        grid.clear();
        let lower = self
            .channels
            .iter()
            .filter(|c| c.expected[0] > 0.0)
            .map(|c| c.source.source().threshold().value())
            .fold(f64::INFINITY, f64::min);
        let upper = self.channels.iter().map(|c| c.cutoff).fold(0.0, f64::max);
        if lower >= upper {
            return Ok(());
        }
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
        for source in &self.channels {
            poll(cancelled)?;
            let photo = source.source.source();
            push(photo.threshold().value())?;
            push(source.cutoff)?;
            if let Some(points) = photo.points() {
                for point in points {
                    poll(cancelled)?;
                    push(point.photon_energy.value())?;
                }
            }
            if source.expected[0] > 0.0 {
                let mut x = 1.0;
                while x <= 1024.0 {
                    let energy =
                        photo.threshold().value() + x * source.spectrum.thermal_energy().value();
                    if energy >= source.cutoff {
                        break;
                    }
                    push(energy)?;
                    x *= 2.0;
                }
            }
        }
        grid.sort_by(f64::total_cmp);
        grid.dedup();
        Ok(())
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
        scratch.values.fill(Dd::ZERO);
        scratch.spectral_errors.fill(Dd::ZERO);
        scratch.angular_errors.fill(Dd::ZERO);
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
                    return Err(Error::Arithmetic("network positive spectral panel"));
                }
                scratch.fine_panel.fill(Dd::ZERO);
                scratch.coarse_panel.fill(Dd::ZERO);
                scratch.angular_panel.fill(Dd::ZERO);
                for (&x, &w) in X16.iter().zip(&W16) {
                    for sign in [-1.0, 1.0] {
                        self.point(
                            mid + sign * half * x,
                            fine_angles,
                            Some(coarse_angles),
                            scratch,
                            cancelled,
                        )?;
                        let weight = Dd::from_f64(half * w);
                        for j in 0..scratch.values.len() {
                            scratch.fine_panel[j] =
                                scratch.fine_panel[j] + weight * scratch.fine_point[j];
                            scratch.angular_panel[j] =
                                scratch.angular_panel[j] + weight * scratch.angular_point[j];
                        }
                    }
                }
                for (&x, &w) in X8.iter().zip(&W8) {
                    for sign in [-1.0, 1.0] {
                        self.point(mid + sign * half * x, fine_angles, None, scratch, cancelled)?;
                        let weight = Dd::from_f64(half * w);
                        for j in 0..scratch.values.len() {
                            scratch.coarse_panel[j] =
                                scratch.coarse_panel[j] + weight * scratch.fine_point[j];
                        }
                    }
                }
                for j in 0..scratch.values.len() {
                    scratch.values[j] = scratch.values[j] + scratch.fine_panel[j];
                    scratch.spectral_errors[j] = scratch.spectral_errors[j]
                        + (scratch.fine_panel[j] - scratch.coarse_panel[j]).abs();
                    scratch.angular_errors[j] = scratch.angular_errors[j]
                        + (scratch.fine_panel[j] - scratch.angular_panel[j]).abs();
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
        self.spend(self.channels.len())?;
        scratch.fine_point.fill(Dd::ZERO);
        scratch.angular_point.fill(Dd::ZERO);
        scratch.cross.fill(Dd::ZERO);
        scratch.cell_cross_heat.fill(Dd::ZERO);
        for (index, source) in self.channels.iter().enumerate() {
            poll(cancelled)?;
            let mut sigma = 0.0;
            let (mut photons, mut captured) = (0.0, 0.0);
            if source.expected[0] > 0.0 && energy <= source.cutoff {
                let density = source
                    .spectrum
                    .evaluate(Energy::new(energy))
                    .map_err(|error| source_error(index, error))?;
                sigma = density.photoionization_cross_section.value();
                photons = density.coefficient_m3_s_per_j * source.pair_density_volume;
                captured = density.captured_kinetic_energy_m3_s * source.pair_density_volume;
            } else if source.opacity_per_sigma > 0.0 {
                // A bounded cold emission tail never switches off its absorber.
                sigma = source
                    .source
                    .source()
                    .cross_section(Energy::new(energy))
                    .map_err(|error| Error::Absorption { index, error })?
                    .value();
            }
            let depth = sigma * source.opacity_per_sigma;
            if !nonnegative(depth)
                || !nonnegative(photons)
                || !nonnegative(captured)
                || (sigma > 0.0 && source.opacity_per_sigma > 0.0 && depth == 0.0)
            {
                return Err(Error::Arithmetic(
                    "network spectral optical depth or emissivity",
                ));
            }
            scratch.opacity[index] = depth;
            scratch.emission[index] = photons;
            for (slot, value) in [photons, energy * photons, captured]
                .into_iter()
                .enumerate()
            {
                if !nonnegative(value) {
                    return Err(Error::Arithmetic("network spectral emission moment"));
                }
                scratch.fine_point[CHANNEL * index + slot] = Dd::from_f64(value);
                scratch.angular_point[CHANNEL * index + slot] = Dd::from_f64(value);
            }
        }
        for c in 0..self.cells() {
            poll(cancelled)?;
            let range = self.offsets[c]..self.offsets[c + 1];
            let maximum = range
                .clone()
                .map(|k| scratch.opacity[k])
                .fold(0.0, f64::max);
            let mut scaled = Dd::ZERO;
            let mut emission = Dd::ZERO;
            for k in range.clone() {
                emission = emission + Dd::from_f64(scratch.emission[k]);
                if maximum > 0.0 {
                    scaled = scaled + Dd::from_f64(scratch.opacity[k] / maximum);
                }
            }
            if !finite(emission) || !finite(scaled) {
                return Err(Error::Arithmetic("mixed source moment"));
            }
            scratch.cell_emission[c] = emission;
            scratch.cell_opacity[c] = maximum;
            scratch.cell_opacity_scale[c] = scaled.to_f64();
            for k in range.clone() {
                scratch.fractions[k] = if maximum == 0.0 {
                    Dd::ZERO
                } else {
                    Dd::from_f64(scratch.opacity[k] / maximum) / scaled
                };
            }
            // Compute each antisymmetric exchange ONCE per energy, outside all
            // angular sweeps. The same-channel term cancels analytically.
            // Sum_i cross_i=0, and the pair electron heat is (chi_j-chi_i)*q.
            let pairs = range
                .len()
                .checked_mul(range.len().saturating_sub(1))
                .ok_or(Error::WorkBudget)?
                / 2;
            self.spend(pairs)?;
            for i in range.clone() {
                for j in i + 1..range.end {
                    let exchange = scratch.fractions[i] * Dd::from_f64(scratch.emission[j])
                        - scratch.fractions[j] * Dd::from_f64(scratch.emission[i]);
                    scratch.cross[i] = scratch.cross[i] + exchange;
                    scratch.cross[j] = scratch.cross[j] - exchange;
                    scratch.cell_cross_heat[c] = scratch.cell_cross_heat[c]
                        + exchange
                            * (Dd::from_f64(self.channels[j].gap)
                                - Dd::from_f64(self.channels[i].gap));
                }
            }
        }
        self.sweep(energy, fine_angles, scratch, false, cancelled)?;
        if let Some(angles) = coarse_angles {
            self.sweep(energy, angles, scratch, true, cancelled)?;
        }
        Ok(())
    }

    fn sweep(
        &mut self,
        energy: f64,
        angles: &[(f64, f64)],
        scratch: &mut Scratch,
        coarse: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let n = self.cells();
        let needed = (n + self.channels.len())
            .checked_mul(2 * angles.len())
            .ok_or(Error::WorkBudget)?;
        self.spend(needed)?;
        let result = if coarse {
            &mut scratch.angular_point
        } else {
            &mut scratch.fine_point
        };
        for &(mu, w) in angles {
            let weight = Dd::from_f64(0.5 * w);
            for right in [true, false] {
                let mut incoming = Dd::ZERO;
                for ordinal in 0..n {
                    poll(cancelled)?;
                    let c = if right { ordinal } else { n - 1 - ordinal };
                    let (transmission, absorption, escape, self_absorption) =
                        characteristic(scratch.cell_opacity[c], scratch.cell_opacity_scale[c], mu);
                    let local = scratch.cell_emission[c] * weight;
                    let external = incoming * Dd::from_f64(absorption);
                    let escaped = local * Dd::from_f64(escape);
                    let outgoing = incoming * Dd::from_f64(transmission) + escaped;
                    let direct = external - escaped;
                    let mut direct_heat =
                        scratch.cell_cross_heat[c] * weight * Dd::from_f64(self_absorption);
                    for k in self.offsets[c]..self.offsets[c + 1] {
                        let source = &self.channels[k];
                        let own = Dd::from_f64(scratch.emission[k]) * weight;
                        let fraction = scratch.fractions[k];
                        let same = own * fraction * Dd::from_f64(self_absorption);
                        let cross_absorbed =
                            (local - own) * fraction * Dd::from_f64(self_absorption);
                        let external_absorbed = external * fraction;
                        let absorbed = same + cross_absorbed + external_absorbed;
                        let net = external_absorbed
                            + scratch.cross[k] * weight * Dd::from_f64(self_absorption)
                            - own * Dd::from_f64(escape);
                        let excess = Dd::from_f64(energy) - Dd::from_f64(source.gap);
                        let heat = if absorbed == Dd::ZERO {
                            Dd::ZERO
                        } else {
                            excess * absorbed
                        };
                        let net_heat = if net == Dd::ZERO {
                            Dd::ZERO
                        } else {
                            excess * net
                        };
                        direct_heat =
                            direct_heat + excess * (external_absorbed - own * Dd::from_f64(escape));
                        if [same, cross_absorbed, external_absorbed, absorbed, heat]
                            .iter()
                            .any(|v| !finite(*v) || v.to_f64() < 0.0)
                            || !finite(net)
                            || !finite(net_heat)
                        {
                            return Err(Error::Arithmetic("network mixed formal characteristic"));
                        }
                        for (slot, value) in [
                            (3, same),
                            (4, cross_absorbed),
                            (5, external_absorbed),
                            (6, Dd::from_f64(energy) * absorbed),
                            (7, heat),
                            (8, net),
                            (9, net_heat),
                            (10, Dd::from_f64(energy) * net),
                        ] {
                            result[CHANNEL * k + slot] = result[CHANNEL * k + slot] + value;
                        }
                    }
                    if !finite(outgoing) || outgoing.to_f64() < 0.0 || !finite(direct_heat) {
                        return Err(Error::Arithmetic("network diffuse outgoing flux"));
                    }
                    let slot = self.cell_slot(c);
                    result[slot] = result[slot] + direct;
                    result[slot + 1] = result[slot + 1] + direct_heat;
                    result[slot + 2] = result[slot + 2] + Dd::from_f64(energy) * direct;
                    let face = if right { c + 1 } else { c };
                    let slot = self.face_slot(face) + if right { 0 } else { 2 };
                    result[slot] = result[slot] + outgoing;
                    result[slot + 1] = result[slot + 1] + Dd::from_f64(energy) * outgoing;
                    incoming = outgoing;
                }
            }
        }
        Ok(())
    }

    fn tail_vector(&self, output: &mut [f64]) {
        output.fill(0.0);
        let mut total = [Dd::ZERO; 2];
        for source in &self.channels {
            for j in 0..2 {
                total[j] = total[j] + Dd::from_f64(source.tail[j]);
            }
        }
        for c in 0..self.cells() {
            let mut local = [Dd::ZERO; 2];
            for k in self.offsets[c]..self.offsets[c + 1] {
                for j in 0..2 {
                    local[j] = local[j] + Dd::from_f64(self.channels[k].tail[j]);
                }
            }
            for k in self.offsets[c]..self.offsets[c + 1] {
                let source = &self.channels[k];
                output[CHANNEL * k..CHANNEL * k + 3].copy_from_slice(&source.tail);
                if source.opacity_per_sigma > 0.0 {
                    output[CHANNEL * k + 3] = source.tail[0];
                    output[CHANNEL * k + 4] =
                        (local[0] - Dd::from_f64(source.tail[0])).to_f64().max(0.0);
                    output[CHANNEL * k + 5] = (total[0] - local[0]).to_f64().max(0.0);
                    output[CHANNEL * k + 6] = total[1].to_f64();
                    output[CHANNEL * k + 7] = total[1].to_f64();
                }
                output[CHANNEL * k + 8] = total[0].to_f64();
                output[CHANNEL * k + 9] = total[1].to_f64();
                output[CHANNEL * k + 10] = total[1].to_f64();
            }
            let slot = self.cell_slot(c);
            output[slot] = total[0].to_f64();
            output[slot + 1] = total[1].to_f64();
            output[slot + 2] = total[1].to_f64();
        }
        for right in [true, false] {
            let mut sum = [Dd::ZERO; 2];
            for ordinal in 0..self.cells() {
                let c = if right {
                    ordinal
                } else {
                    self.cells() - 1 - ordinal
                };
                for k in self.offsets[c]..self.offsets[c + 1] {
                    for j in 0..2 {
                        sum[j] = sum[j] + Dd::from_f64(self.channels[k].tail[j]);
                    }
                }
                let face = if right { c + 1 } else { c };
                let slot = self.face_slot(face) + if right { 0 } else { 2 };
                for j in 0..2 {
                    output[slot + j] = 0.5 * sum[j].to_f64();
                }
            }
        }
    }

    fn errors(&self, scratch: &Scratch) -> Result<[f64; 3], Error> {
        let mut errors = [0.0_f64; 3];
        let channel_end = CHANNEL * self.channels.len();
        let cell_end = channel_end + CELL * self.cells();
        for (j, value) in scratch.values.iter().enumerate() {
            let signed = if j < channel_end {
                matches!(j % CHANNEL, 8..=10)
            } else {
                j < cell_end
            };
            if !finite(*value) || (!signed && value.to_f64() < 0.0) {
                return Err(Error::Arithmetic("network integrated radiation moment"));
            }
            let allowance = self.allowance(j, value.to_f64());
            for (kind, error) in [
                scratch.spectral_errors[j].to_f64(),
                scratch.angular_errors[j].to_f64(),
                scratch.tail[j],
            ]
            .into_iter()
            .enumerate()
            {
                let normalized = error / allowance;
                if !nonnegative(normalized) {
                    return Err(Error::Arithmetic("network transport error estimate"));
                }
                errors[kind] = errors[kind].max(normalized);
            }
        }
        for (k, source) in self.channels.iter().enumerate() {
            for j in 0..3 {
                let mismatch = ((scratch.values[CHANNEL * k + j].to_f64() - source.expected[j])
                    .abs()
                    - source.tail[j])
                    .max(0.0);
                errors[0] =
                    errors[0].max(mismatch / self.allowance(CHANNEL * k + j, source.expected[j]));
            }
        }
        Ok(errors)
    }
}

/// Exact uniform-volume characteristic. Thin absorption uses expm1 and its
/// small-depth series; optically thick escape is retained as 1/tau.
fn characteristic(maximum: f64, scale: f64, mu: f64) -> (f64, f64, f64, f64) {
    if maximum == 0.0 {
        return (1.0, 0.0, 1.0, 0.0);
    }
    let depth = maximum * (scale / mu);
    if depth.is_infinite() {
        // The sum of individually finite channel depths can overflow even
        // though its uniform-source escape remains representable.
        return (0.0, 1.0, (mu / maximum) / scale, 1.0);
    }
    let absorption = -expm1(-depth);
    let transmission = exp(-depth);
    let escape = absorption / depth;
    let self_absorption = if depth < 1e-3 {
        depth
            * (0.5
                + depth
                    * (-1.0 / 6.0 + depth * (1.0 / 24.0 + depth * (-1.0 / 120.0 + depth / 720.0))))
    } else {
        1.0 - escape
    };
    (transmission, absorption, escape, self_absorption)
}

fn angles<const N: usize>(
    panels: usize,
    nodes: &[f64; N],
    weights: &[f64; N],
    output: &mut Vec<(f64, f64)>,
) -> Result<(), Error> {
    output.clear();
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
    Ok(())
}
