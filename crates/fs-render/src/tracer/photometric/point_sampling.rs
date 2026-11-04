//! Bounded point-source proposals on the existing calibrated path vertices.
//!
//! Each independent reservoir replaces its candidate with probability
//! `w_new / (w_previous_sum + w_new)`. With exact arithmetic and continuous
//! uniforms, its final marginal is therefore `w/W`.
//! Actual connections retain the native visibility, spectral transport and
//! receiver response; the returned packet averages `C(selected)/(w/W)`.
//! Numerical draws use the 53-bit grid in `[0,1)`. A comparison rounds its
//! Bernoulli probability upward by less than `2^-53`; this and floating-point
//! summation remain numerical approximations, not certified error bounds.
//!
//! The opt-in streaming CDF policy instead traverses the candidates twice and
//! draws only once per requested connection. Its positive integer intervals
//! partition that same 53-bit grid; each selected connection uses its actual
//! represented interval probability. Floating-point weights choose a proposal,
//! without changing the original physical source or connection evaluation.

use super::*;
use crate::photometry::path::MAX_POINT_CONNECTIONS_PER_VERTEX;

const VERTEX_DOMAIN: u32 = 0x706e_7476;
const CANDIDATE_DOMAIN: u32 = 0x706e_7463;
const CDF_VERTEX_DOMAIN: u32 = 0x7063_6476;
const CDF_DRAW_DOMAIN: u32 = 0x7063_6464;
const CAPACITY: usize = MAX_POINT_CONNECTIONS_PER_VERTEX as usize;
const UNIFORM_BINS: u64 = 1 << 53;
const UNIFORM_SPACING: f64 = 1.0 / 9_007_199_254_740_992.0;

fn bits53(high: u32, low: u32) -> u64 {
    (u64::from(high) << 21) | (u64::from(low) >> 11)
}

/// Two native counter words supply all 53 exactly represented fraction bits.
/// The greatest value is `1-2^-53`, retaining possible survival when a new
/// source has almost all of the accumulated importance.
fn unit53(high: u32, low: u32) -> f64 {
    bits53(high, low) as f64 * UNIFORM_SPACING
}

#[derive(Clone, Copy)]
struct Candidate {
    source: usize,
    direction: Vec3,
    distance: f64,
    incident_lux: f64,
    importance: f64,
}

/// All active source distances are admitted before beam/receiver rejection.
/// Both proposal policies use this same operation order and source owner.
fn candidate(
    source: &PhotometricSource,
    index: usize,
    point: Point3,
    options: &TransportOptions,
    proposal_response: &mut impl FnMut(Vec3) -> Result<f64, TracerError>,
) -> Result<Option<Candidate>, TracerError> {
    if source.flux_lm() == 0.0 {
        return Ok(None);
    }
    let displacement = source.position_m().delta_from(point);
    let distance =
        norm(displacement).map_err(|_| connection_error(index, "invalid source distance"))?;
    if distance <= options.ray_epsilon_m {
        return Err(connection_error(
            index,
            "receiver enters point-source exclusion radius",
        ));
    }
    if distance > options.max_distance_m {
        return Err(connection_error(index, "source exceeds distance budget"));
    }
    let direction = displacement.scale(1.0 / distance);
    let intensity = source
        .intensity_cd(direction.scale(-1.0))
        .map_err(TracerError::FixtureSource)?;
    if intensity == 0.0 {
        return Ok(None);
    }
    let geometric_response = proposal_response(direction)?;
    if !geometric_response.is_finite() || geometric_response < 0.0 {
        return Err(connection_error(index, "invalid fixture proposal response"));
    }
    if geometric_response == 0.0 {
        return Ok(None);
    }
    let incident_lux = intensity / distance / distance;
    Ok(Some(Candidate {
        source: index,
        direction,
        distance,
        incident_lux,
        importance: incident_lux * geometric_response,
    }))
}

/// Sum positive proposal weights in units of their largest member. No physical
/// intensity is clamped or changed. If distinct positive support cannot survive
/// this representation, refuse rather than silently removing a fixture.
struct WeightSum {
    maximum: f64,
    minimum: f64,
    scaled: f64,
}

impl WeightSum {
    const EMPTY: Self = Self {
        maximum: 0.0,
        minimum: f64::INFINITY,
        scaled: 0.0,
    };

    fn insert(&mut self, weight: f64, source: usize) -> Result<f64, TracerError> {
        if !weight.is_finite() || weight <= 0.0 {
            return Err(connection_error(
                source,
                "unrepresentable fixture importance",
            ));
        }
        if self.maximum == 0.0 {
            self.maximum = weight;
            self.minimum = weight;
            self.scaled = 1.0;
            return Ok(1.0);
        }
        let maximum = self.maximum.max(weight);
        let previous = self.scaled * (self.maximum / maximum);
        let incoming = weight / maximum;
        let sum = previous + incoming;
        if !(previous > 0.0
            && incoming > 0.0
            && sum.is_finite()
            && sum > previous
            && sum > incoming)
        {
            return Err(connection_error(
                source,
                "fixture importance loses positive proposal support",
            ));
        }
        self.maximum = maximum;
        self.minimum = self.minimum.min(weight);
        self.scaled = sum;
        let smallest = self.minimum / maximum;
        if smallest / sum < UNIFORM_SPACING || sum + smallest == sum {
            return Err(connection_error(
                source,
                "fixture importance loses positive proposal support",
            ));
        }
        Ok(incoming / sum)
    }

    fn probability(&self, candidate: Candidate) -> Result<f64, TracerError> {
        let probability = candidate.importance / self.maximum / self.scaled;
        if !probability.is_finite() || !(0.0 < probability && probability <= 1.0) {
            return Err(connection_error(
                candidate.source,
                "unrepresentable fixture selection probability",
            ));
        }
        Ok(probability)
    }
}

/// A streaming CDF on the exact integer draw grid. Repeating the original
/// prefix arithmetic makes the complete prefix normalize to exactly one.
/// Quantizing endpoints changes only the proposal; the returned PDF counts
/// its actual integer outcomes, including any endpoint-rounding adjustment.
struct CdfIntervals {
    prefix: WeightSum,
    end: u64,
}

impl CdfIntervals {
    const EMPTY: Self = Self {
        prefix: WeightSum::EMPTY,
        end: 0,
    };

    fn next(
        &mut self,
        weight: f64,
        source: usize,
        total: &WeightSum,
    ) -> Result<(u64, f64), TracerError> {
        self.prefix.insert(weight, source)?;
        let endpoint = (self.prefix.maximum / total.maximum) * (self.prefix.scaled / total.scaled);
        if !endpoint.is_finite() || !(0.0 < endpoint && endpoint <= 1.0) {
            return Err(connection_error(
                source,
                "unrepresentable fixture CDF endpoint",
            ));
        }
        // A raw integer draw d obeys d/2^53 < endpoint precisely when
        // d < ceil(endpoint*2^53). Both the threshold and its differences are
        // exactly represented integers, so the resulting PDF is dyadic.
        let end = (endpoint * UNIFORM_BINS as f64).ceil() as u64;
        if end <= self.end {
            return Err(connection_error(
                source,
                "fixture CDF loses positive proposal support",
            ));
        }
        let probability = (end - self.end) as f64 * UNIFORM_SPACING;
        self.end = end;
        Ok((end, probability))
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn contribution(
    scene: &Scene,
    fixtures: &FixtureLighting<'_>,
    cx: &Cx<'_>,
    medium: Option<PathMedium<'_, '_>>,
    point: Point3,
    wavelengths: &[f64; PACKET],
    rng: &PathRng,
    depth: u32,
    connections: u32,
    mut proposal_response: impl FnMut(Vec3) -> Result<f64, TracerError>,
    response: impl FnMut(Vec3, f64) -> Result<f64, TracerError>,
) -> Result<[f64; PACKET], TracerError> {
    if connections == 0 || connections > MAX_POINT_CONNECTIONS_PER_VERTEX {
        return Err(TracerError::InvalidInput);
    }
    let count = connections as usize;
    let mut selected = [None; CAPACITY];
    let mut weights = WeightSum::EMPTY;
    // A separate domain never advances or replaces a camera, BSDF, phase or
    // area-light draw. The source index retains the caller's stable ordering.
    let words = philox4x32_10([rng.pixel, rng.sample, depth, VERTEX_DOMAIN], rng.key);
    let key = [words[0] ^ words[2], words[1] ^ words[3]];
    for (index, source) in fixtures.sources.iter().enumerate() {
        cx.checkpoint()?;
        let Some(candidate) = candidate(
            source,
            index,
            point,
            &fixtures.options,
            &mut proposal_response,
        )?
        else {
            continue;
        };
        let replacement_probability = weights.insert(candidate.importance, index)?;
        let source_index = u64::try_from(index).map_err(|_| TracerError::InvalidInput)?;
        // Two independent 53-bit reservoirs share one native Philox block. The
        // inner loop has a fixed ceiling of 32; it performs no transport work.
        for block in 0..connections.div_ceil(2) {
            let draws = philox4x32_10(
                [
                    source_index as u32,
                    (source_index >> 32) as u32,
                    block,
                    CANDIDATE_DOMAIN,
                ],
                key,
            );
            for lane in 0..2 {
                let slot = block as usize * 2 + lane;
                let draw = unit53(draws[2 * lane], draws[2 * lane + 1]);
                if slot < count && draw < replacement_probability {
                    selected[slot] = Some(candidate);
                }
            }
        }
    }
    selected_contribution(
        scene,
        fixtures,
        cx,
        medium,
        wavelengths,
        connections,
        &selected,
        |_, candidate| weights.probability(candidate),
        response,
    )
}

/// Two complete source sweeps, bounded independent draws and no heap storage.
/// Candidate evaluation remains receiver-dependent, including every active
/// source's distance refusal. Only proposal construction differs from the
/// reservoir policy; selected connections use the same physical evaluator.
#[allow(clippy::too_many_arguments)]
pub(super) fn contribution_cdf(
    scene: &Scene,
    fixtures: &FixtureLighting<'_>,
    cx: &Cx<'_>,
    medium: Option<PathMedium<'_, '_>>,
    point: Point3,
    wavelengths: &[f64; PACKET],
    rng: &PathRng,
    depth: u32,
    connections: u32,
    mut proposal_response: impl FnMut(Vec3) -> Result<f64, TracerError>,
    response: impl FnMut(Vec3, f64) -> Result<f64, TracerError>,
) -> Result<[f64; PACKET], TracerError> {
    if connections == 0 || connections > MAX_POINT_CONNECTIONS_PER_VERTEX {
        return Err(TracerError::InvalidInput);
    }
    cx.checkpoint()?;
    let mut weights = WeightSum::EMPTY;
    for (index, source) in fixtures.sources.iter().enumerate() {
        cx.checkpoint()?;
        if let Some(candidate) = candidate(
            source,
            index,
            point,
            &fixtures.options,
            &mut proposal_response,
        )? {
            weights.insert(candidate.importance, index)?;
        }
    }
    cx.checkpoint()?;
    if weights.maximum == 0.0 {
        return Ok([0.0; PACKET]);
    }
    let count = connections as usize;
    let mut draws = [(0_u64, 0_usize); CAPACITY];
    let words = philox4x32_10([rng.pixel, rng.sample, depth, CDF_VERTEX_DOMAIN], rng.key);
    let key = [words[0] ^ words[2], words[1] ^ words[3]];
    for block in 0..connections.div_ceil(2) {
        let words = philox4x32_10([block, 0, 0, CDF_DRAW_DOMAIN], key);
        for lane in 0..2 {
            let slot = block as usize * 2 + lane;
            if slot < count {
                draws[slot] = (bits53(words[2 * lane], words[2 * lane + 1]), slot);
            }
        }
    }
    cx.checkpoint()?;
    draws[..count].sort_unstable();
    cx.checkpoint()?;
    let mut cdf = CdfIntervals::EMPTY;
    let mut selected = [None; CAPACITY];
    let mut probabilities = [0.0; CAPACITY];
    let mut next_draw = 0;
    for (index, source) in fixtures.sources.iter().enumerate() {
        cx.checkpoint()?;
        let Some(candidate) = candidate(
            source,
            index,
            point,
            &fixtures.options,
            &mut proposal_response,
        )?
        else {
            continue;
        };
        let (end, probability) = cdf.next(candidate.importance, index, &weights)?;
        while next_draw < count && draws[next_draw].0 < end {
            let slot = draws[next_draw].1;
            selected[slot] = Some(candidate);
            probabilities[slot] = probability;
            next_draw += 1;
        }
        // Even after all draws are selected, later positive sources must
        // retain represented support and all cancellation checks still run.
    }
    cx.checkpoint()?;
    if cdf.end != UNIFORM_BINS || next_draw != count {
        return Err(TracerError::InvalidInput);
    }
    selected_contribution(
        scene,
        fixtures,
        cx,
        medium,
        wavelengths,
        connections,
        &selected,
        |slot, _| Ok(probabilities[slot]),
        response,
    )
}

#[allow(clippy::too_many_arguments)]
fn selected_contribution(
    scene: &Scene,
    fixtures: &FixtureLighting<'_>,
    cx: &Cx<'_>,
    medium: Option<PathMedium<'_, '_>>,
    wavelengths: &[f64; PACKET],
    connections: u32,
    selected: &[Option<Candidate>; CAPACITY],
    mut probability: impl FnMut(usize, Candidate) -> Result<f64, TracerError>,
    mut response: impl FnMut(Vec3, f64) -> Result<f64, TracerError>,
) -> Result<[f64; PACKET], TracerError> {
    let count = connections as usize;
    let mut result = [0.0; PACKET];
    let mut correction = [0.0; PACKET];
    for slot in 0..count {
        cx.checkpoint()?;
        let Some(candidate) = selected[slot] else {
            continue;
        };
        if selected[..slot]
            .iter()
            .flatten()
            .any(|previous| previous.source == candidate.source)
        {
            continue;
        }
        // With-replacement draws can choose the same fixture repeatedly.
        // Coalescing its multiplicity keeps the estimator and sample stream,
        // while performing its shadow/Beer/BSDF evaluation exactly once.
        let multiplicity = selected[..count]
            .iter()
            .flatten()
            .filter(|entry| entry.source == candidate.source)
            .count();
        let probability = probability(slot, candidate)?;
        let factor = (multiplicity as f64 / f64::from(connections)) / probability;
        if !factor.is_finite() || factor <= 0.0 {
            return Err(connection_error(
                candidate.source,
                "unrepresentable fixture selection weight",
            ));
        }
        let source = &fixtures.sources[candidate.source];
        let ray = Ray {
            origin: source.position_m(),
            dir: candidate.direction.scale(-1.0),
        };
        if fixtures
            .geometry
            .intersect(
                &scene.primitives,
                cx,
                &ray,
                candidate.distance - fixtures.options.ray_epsilon_m,
            )?
            .is_some()
        {
            continue;
        }
        let transmittance = if let Some(medium) = medium {
            medium.transmittance(cx, ray.origin, ray.dir, candidate.distance, wavelengths)?
        } else {
            [1.0; PACKET]
        };
        for (lane, &lambda) in wavelengths.iter().enumerate() {
            let value = candidate.incident_lux
                * source.spectrum().value_w_per_lumen_nm(lambda)
                * response(candidate.direction, lambda)?
                * transmittance[lane]
                * factor;
            if !value.is_finite() || value < 0.0 {
                return Err(connection_error(
                    candidate.source,
                    "nonfinite spectral fixture contribution",
                ));
            }
            let sum = result[lane] + value;
            correction[lane] += if result[lane].abs() >= value.abs() {
                (result[lane] - sum) + value
            } else {
                (value - sum) + result[lane]
            };
            result[lane] = sum;
        }
    }
    for lane in 0..PACKET {
        result[lane] += correction[lane];
        if !result[lane].is_finite() || result[lane] < 0.0 {
            return Err(TracerError::InvalidInput);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::RigidTransform;
    use crate::photometry::source::{AngularDistribution, SourceSpectrum};

    fn summed(weights: &[f64]) -> WeightSum {
        let mut total = WeightSum::EMPTY;
        for (source, &weight) in weights.iter().enumerate() {
            total.insert(weight, source).unwrap();
        }
        total
    }

    fn intervals(weights: &[f64]) -> Vec<(u64, f64)> {
        let total = summed(weights);
        let mut cdf = CdfIntervals::EMPTY;
        weights
            .iter()
            .enumerate()
            .map(|(source, &weight)| cdf.next(weight, source, &total).unwrap())
            .collect()
    }

    #[test]
    fn fifty_three_bit_draws_retain_both_sides_of_extreme_reservoir_replacements() {
        let lowest = unit53(0, 0);
        let highest = unit53(u32::MAX, u32::MAX);
        assert_eq!(bits53(0, 0), 0);
        assert_eq!(bits53(0, 1 << 11), 1);
        assert_eq!(bits53(u32::MAX, u32::MAX), UNIFORM_BINS - 1);
        assert_eq!(lowest, 0.0);
        assert_eq!(unit53(0, 1 << 11), UNIFORM_SPACING);
        assert_eq!(highest, 1.0 - UNIFORM_SPACING);
        let mut weights = WeightSum::EMPTY;
        weights.insert(1.0e-12, 0).unwrap();
        let probability = weights.insert(1.0, 1).unwrap();
        // Every 32-bit [0,1) draw would replace the earlier source here.
        // The new grid admits both replacement and survival without changing
        // either physical source or clamping their importance ratio.
        assert!(probability > u32_unit(u32::MAX));
        assert!(lowest < probability && highest >= probability);
        let mut reverse = WeightSum::EMPTY;
        reverse.insert(1.0, 0).unwrap();
        let probability = reverse.insert(1.0e-12, 1).unwrap();
        assert!(lowest < probability && highest >= probability);

        // All positive candidates must retain at least one fraction interval.
        // A tiny early source cannot silently lose support as the sum grows.
        let mut unsupported = WeightSum::EMPTY;
        unsupported.insert(1.0e-14, 0).unwrap();
        let refusal = (1..=100).find(|&source| unsupported.insert(1.0, source).is_err());
        assert!(refusal.is_some());
    }

    #[test]
    fn cdf_boundaries_and_enumerated_dyadic_draws_preserve_inverse_probability_sum() {
        let bins = intervals(&[1.0, 2.0, 1.0]);
        assert_eq!(
            bins,
            [
                (UNIFORM_BINS / 4, 0.25),
                (3 * UNIFORM_BINS / 4, 0.5),
                (UNIFORM_BINS, 0.25),
            ]
        );
        for (draw, expected_source) in [
            (0, 0),
            (UNIFORM_BINS / 4 - 1, 0),
            (UNIFORM_BINS / 4, 1),
            (3 * UNIFORM_BINS / 4 - 1, 1),
            (3 * UNIFORM_BINS / 4, 2),
            (UNIFORM_BINS - 1, 2),
        ] {
            assert_eq!(
                bins.partition_point(|&(end, _)| draw >= end),
                expected_source
            );
        }
        // Eight equal strata exactly enumerate this dyadic distribution.
        // Use an independent binary lookup and unequal physical values so a
        // missing PDF correction cannot pass by counting selections alone.
        let values = [2.0, 7.0, 11.0];
        let mut counts = [0; 3];
        let mut estimate = 0.0;
        for stratum in 0..8 {
            let draw = stratum * (UNIFORM_BINS / 8);
            let source = bins.partition_point(|&(end, _)| draw >= end);
            counts[source] += 1;
            estimate += values[source] / bins[source].1 / 8.0;
        }
        assert_eq!(counts, [2, 4, 2]);
        assert_eq!(estimate, values.iter().sum::<f64>());
    }

    #[test]
    fn cdf_pdf_uses_represented_integer_width_instead_of_ideal_weight_ratio() {
        let bins = intervals(&[1.0; 3]);
        let third = UNIFORM_BINS / 3;
        assert_eq!(
            bins.iter().map(|&(end, _)| end).collect::<Vec<_>>(),
            [third + 1, 2 * third + 1, UNIFORM_BINS]
        );
        let expected_widths = [third + 1, third, third + 1];
        for ((_, pdf), width) in bins.iter().zip(expected_widths) {
            assert_eq!(pdf.to_bits(), (width as f64 * UNIFORM_SPACING).to_bits());
        }
        assert_ne!(bins[0].1.to_bits(), (1.0_f64 / 3.0).to_bits());
        assert_eq!(expected_widths.iter().sum::<u64>(), UNIFORM_BINS);
        assert_eq!(bins.iter().map(|&(_, pdf)| pdf).sum::<f64>(), 1.0);
    }

    #[test]
    fn cdf_rescaling_keeps_extreme_positive_support_and_overflowing_total_weights() {
        assert_eq!(
            intervals(&[f64::MAX, f64::MAX]),
            [(UNIFORM_BINS / 2, 0.5), (UNIFORM_BINS, 0.5)]
        );
        for weights in [
            [1.0e-12, 1.0, 1.0e-12],
            [1.0e-12, 1.0e-6, 1.0],
            [1.0, 1.0e-6, 1.0e-12],
        ] {
            let bins = intervals(&weights);
            let mut previous = 0;
            for &(end, pdf) in &bins {
                assert!(end > previous);
                assert_eq!(pdf, (end - previous) as f64 * UNIFORM_SPACING);
                previous = end;
            }
            assert_eq!(previous, UNIFORM_BINS);
        }
        let mut unsupported = WeightSum::EMPTY;
        unsupported.insert(1.0, 0).unwrap();
        assert!(matches!(
            unsupported.insert(1.0e-30, 1),
            Err(TracerError::FixtureConnection { source: 1, .. })
        ));
        let mut cdf = CdfIntervals::EMPTY;
        assert!(matches!(
            cdf.next(1.0, 7, &WeightSum::EMPTY),
            Err(TracerError::FixtureConnection { source: 7, .. })
        ));
    }

    #[test]
    fn cdf_refuses_quantized_support_loss_even_when_weight_sum_admits() {
        let weights = [1.0, 9.0 / 36_028_797_018_963_968.0, 1.1];
        let total = summed(&weights);
        let mut cdf = CdfIntervals::EMPTY;
        assert_eq!(
            cdf.next(weights[0], 0, &total).unwrap().0,
            4_289_142_502_257_616
        );
        // The second positive prefix rounds to the same integer endpoint.
        // Its nominal weight survives summation, but it has no selectable
        // integer outcome and must refuse even if earlier draws were filled.
        assert!(matches!(
            cdf.next(weights[1], 1, &total),
            Err(TracerError::FixtureConnection { source: 1, .. })
        ));
    }

    #[test]
    fn all_active_source_distances_refuse_before_zero_beam_or_receiver_response() {
        let options = TransportOptions {
            ray_epsilon_m: 1.0e-3,
            max_distance_m: 10.0,
            ..TransportOptions::default()
        };
        for (distance, reason) in [
            (
                options.ray_epsilon_m / 2.0,
                "receiver enters point-source exclusion radius",
            ),
            (
                options.max_distance_m * 2.0,
                "source exceeds distance budget",
            ),
        ] {
            // The +Z cone points away from the receiver. A proposal that
            // skips zero beam/response candidates before admission is wrong.
            let source = PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, distance]).unwrap(),
                100.0,
                1.0,
                AngularDistribution::uniform_cone(0.1).unwrap(),
                SourceSpectrum::equal_energy(),
            )
            .unwrap();
            let result = candidate(
                &source,
                17,
                Point3::new(0.0, 0.0, 0.0),
                &options,
                &mut |_| Ok(0.0),
            );
            assert!(matches!(
                result,
                Err(TracerError::FixtureConnection {
                    source: 17,
                    reason: observed,
                }) if observed == reason
            ));
        }
    }
}
