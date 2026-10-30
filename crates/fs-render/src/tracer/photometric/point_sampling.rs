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

use super::*;
use crate::photometry::path::MAX_POINT_CONNECTIONS_PER_VERTEX;

const VERTEX_DOMAIN: u32 = 0x706e_7476;
const CANDIDATE_DOMAIN: u32 = 0x706e_7463;
const CAPACITY: usize = MAX_POINT_CONNECTIONS_PER_VERTEX as usize;
const UNIFORM_SPACING: f64 = 1.0 / 9_007_199_254_740_992.0;

/// Two native counter words supply all 53 exactly represented fraction bits.
/// The greatest value is `1-2^-53`, retaining possible survival when a new
/// source has almost all of the accumulated importance.
fn unit53(high: u32, low: u32) -> f64 {
    let bits = (u64::from(high) << 21) | (u64::from(low) >> 11);
    bits as f64 * UNIFORM_SPACING
}

#[derive(Clone, Copy)]
struct Candidate {
    source: usize,
    direction: Vec3,
    distance: f64,
    incident_lux: f64,
    importance: f64,
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
    mut response: impl FnMut(Vec3, f64) -> Result<f64, TracerError>,
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
        if source.flux_lm() == 0.0 {
            continue;
        }
        let displacement = source.position_m().delta_from(point);
        let distance =
            norm(displacement).map_err(|_| connection_error(index, "invalid source distance"))?;
        if distance <= fixtures.options.ray_epsilon_m {
            return Err(connection_error(
                index,
                "receiver enters point-source exclusion radius",
            ));
        }
        if distance > fixtures.options.max_distance_m {
            return Err(connection_error(index, "source exceeds distance budget"));
        }
        let direction = displacement.scale(1.0 / distance);
        let intensity = source
            .intensity_cd(direction.scale(-1.0))
            .map_err(TracerError::FixtureSource)?;
        if intensity == 0.0 {
            continue;
        }
        let geometric_response = proposal_response(direction)?;
        if !geometric_response.is_finite() || geometric_response < 0.0 {
            return Err(connection_error(index, "invalid fixture proposal response"));
        }
        if geometric_response == 0.0 {
            continue;
        }
        let incident_lux = intensity / distance / distance;
        let candidate = Candidate {
            source: index,
            direction,
            distance,
            incident_lux,
            importance: incident_lux * geometric_response,
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
        let probability = weights.probability(candidate)?;
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
        if photometric_intersect(
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

    #[test]
    fn fifty_three_bit_draws_retain_both_sides_of_extreme_reservoir_replacements() {
        let lowest = unit53(0, 0);
        let highest = unit53(u32::MAX, u32::MAX);
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
}
