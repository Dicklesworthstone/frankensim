//! Continuous volumetric caloric data without invented mass-specific inputs.
//!
//! The storage primitive integrates the selected C_vol(T) at a fixed volume
//! convention. It is not automatically internal energy, enthalpy, or a moving
//! body's energy: the caller supplies the caloric constraint and geometry.
//! No mechanical work, latent heat, moisture transport or phase change is added.

use super::super::{positive, source_value};
use super::{
    ThermalTransportError, admit_keys, admit_temperature_point, source_line, transport_values,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};

/// Original volumetric capacity and scalar conductivity on one fixed context.
#[derive(Clone, Copy, Debug)]
pub struct VolumetricThermalPathQuery<'a> {
    pub volumetric_heat_capacity: &'a PropertyKey,
    pub conductivity: &'a PropertyKey,
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    /// Capacity, then conductivity; neither choice implies independent samples.
    pub selections: [ClaimSelection; 2],
}

/// Algebraic response of the retained source models, not a new measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumetricThermalSample {
    path_identity: ContentHash,
    temperature_k: f64,
    values: [f64; 4],
}

impl VolumetricThermalSample {
    #[must_use]
    pub const fn path_identity(&self) -> ContentHash {
        self.path_identity
    }
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    #[must_use]
    pub const fn volumetric_heat_capacity_j_m3_k(&self) -> f64 {
        self.values[0]
    }
    #[must_use]
    pub const fn conductivity_w_m_k(&self) -> f64 {
        self.values[1]
    }
    #[must_use]
    pub const fn diffusivity_m2_s(&self) -> f64 {
        self.values[2]
    }
    #[must_use]
    pub const fn effusivity_j_m2_k_sqrt_s(&self) -> f64 {
        self.values[3]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Knot {
    temperature: f64,
    values: [f64; 2],
}

/// Immutable source-line preparation for repeated transport and caloric queries.
///
/// Two canonical line budgets bound preparation and caloric traversal. Source
/// breaks and original receipts survive; no density or mass-specific heat is
/// fabricated. Point evaluation uses one binary search, then derives transport
/// from the interpolated independent C_vol and k. All successful runtime methods
/// allocate nothing, perform no source queries and leave the path unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedVolumetricThermalTransport {
    material_card: ContentHash,
    temperature_axis: String,
    endpoints: [QueryPoint; 2],
    sources: [PropertyLineSlice; 2],
    knots: Vec<Knot>,
    identity: ContentHash,
}

impl PreparedVolumetricThermalTransport {
    pub fn try_new(
        card: &MaterialCard,
        query: VolumetricThermalPathQuery<'_>,
    ) -> Result<Self, ThermalTransportError> {
        let keys = [query.volumetric_heat_capacity, query.conductivity];
        admit_keys(keys)?;
        for point in [query.lower, query.upper] {
            admit_temperature_point(point, query.temperature_axis)?;
        }
        if query.lower.axis_quantities() != query.upper.axis_quantities()
            || !query.lower.axes().keys().eq(query.upper.axes().keys())
            || query.lower.axes().iter().any(|(axis, value)| {
                axis != query.temperature_axis && query.upper.axes()[axis] != *value
            })
            || query.lower.axes()[query.temperature_axis]
                >= query.upper.axes()[query.temperature_axis]
        {
            return Err(invalid(
                "a positive-width temperature path must retain its complete fixed source context",
            ));
        }
        let resolve = |index| {
            source_line(
                card,
                keys[index],
                query.temperature_axis,
                query.lower,
                query.upper,
                query.selections[index],
            )
        };
        let sources = [resolve(0)?, resolve(1)?];
        let mut temperatures: Vec<_> = sources
            .iter()
            .flat_map(|source| source.knots.iter().map(|(t, _)| *t))
            .collect();
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        if temperatures.len() > 2 * fs_matdb::MAX_PROPERTY_LINE_KNOTS {
            return Err(invalid(
                "volumetric source union exceeds its breakpoint budget",
            ));
        }
        let material_card = card.content_hash();
        let mut hasher =
            DomainHasher::new("org.frankensim.fs-material.prepared-volumetric-thermal.v1");
        hasher.update(material_card.as_bytes());
        hasher.update(&(query.temperature_axis.len() as u64).to_le_bytes());
        hasher.update(query.temperature_axis.as_bytes());
        for source in &sources {
            hasher.update(&(source.knots.len() as u64).to_le_bytes());
            for (_, answer) in &source.knots {
                hasher.update(answer.receipt.content_hash().as_bytes());
            }
        }
        let identity = hasher.finalize();
        let mut knots = Vec::with_capacity(temperatures.len());
        for temperature in temperatures {
            let values = sources
                .each_ref()
                .map(|source| source_value(source, temperature));
            transport_values(values[0], values[1])?;
            knots.push(Knot {
                temperature,
                values,
            });
        }
        Ok(Self {
            material_card,
            temperature_axis: query.temperature_axis.to_owned(),
            endpoints: [query.lower.clone(), query.upper.clone()],
            sources,
            knots,
            identity,
        })
    }

    /// Evaluate independent source curves before forming k/C_vol and effusivity.
    /// No extrapolation occurs. Extreme interior arithmetic can refuse even when
    /// the endpoints are representable, without mutating the prepared object.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<VolumetricThermalSample, ThermalTransportError> {
        let values = self.values(temperature_k)?;
        let [diffusivity, effusivity] = transport_values(values[0], values[1])?;
        Ok(VolumetricThermalSample {
            path_identity: self.identity,
            temperature_k,
            values: [values[0], values[1], diffusivity, effusivity],
        })
    }

    /// Signed integral of C_vol(T) dT [J/m3] in the selected fixed volume basis.
    /// Exact for the piecewise-affine source model, subject to floating-point
    /// roundoff; it is not a claim about unsampled material or source uncertainty.
    /// Integrating crossed segments directly avoids subtracting large absolute
    /// primitives for a small temperature change. Reversal is antisymmetric.
    pub fn storage_change_j_m3(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.check_temperature(initial_k)?;
        self.check_temperature(final_k)?;
        if initial_k == final_k {
            return Ok(0.0);
        }
        let (low, high, sign) = if initial_k < final_k {
            (initial_k, final_k, 1.0)
        } else {
            (final_k, initial_k, -1.0)
        };
        let mut i = self.segment(low);
        let mut current = low;
        let mut sum = 0.0;
        let mut correction = 0.0;
        while current < high {
            let end = high.min(self.knots[i + 1].temperature);
            let amount = self.segment_storage(i, current, end)?;
            let corrected = amount - correction;
            let next = sum + corrected;
            if !positive(next) {
                return Err(invalid(
                    "volumetric storage is not finite and representable",
                ));
            }
            correction = (next - sum) - corrected;
            sum = next;
            current = end;
            i += 1;
        }
        Ok(sign * sum)
    }

    /// Invert the same fixed-volume storage change for heating or cooling.
    /// Positive capacity makes the supported inverse unique. Traversal is bounded
    /// by the retained knots; the final affine-capacity segment uses a scaled
    /// quadratic inverse without iteration. No clock, power, geometry or state
    /// is owned here. The caller converts its heat to this same volume basis.
    /// An unsupported target or unrepresentable nonzero temperature change refuses.
    pub fn temperature_after_storage_change(
        &self,
        initial_k: f64,
        change_j_m3: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.check_temperature(initial_k)?;
        if !change_j_m3.is_finite() {
            return Err(invalid("storage change must be finite"));
        }
        if change_j_m3 == 0.0 {
            return Ok(initial_k);
        }
        let heating = change_j_m3 > 0.0;
        let span = self.temperature_span_k();
        let boundary = if heating { span[1] } else { span[0] };
        let target = change_j_m3.abs();
        // Exact finite endpoint targets retain their original temperature.
        // A nonrepresentable full-interval integral does not preclude a small
        // representable change inside it; the bounded traversal still decides.
        if let Ok(available) = self.storage_change_j_m3(initial_k, boundary) {
            let available = available.abs();
            if target > available {
                return Err(invalid(
                    "storage target exceeds the original source temperature support",
                ));
            }
            if target == available {
                return Ok(boundary);
            }
        }
        let mut current = initial_k;
        let mut cumulative = 0.0;
        let mut correction = 0.0;
        let mut i = if heating {
            self.segment(initial_k)
        } else {
            self.knots
                .partition_point(|k| k.temperature < initial_k)
                .saturating_sub(1)
        };
        for _ in 0..self.knots.len() - 1 {
            let end = if heating {
                self.knots[i + 1].temperature
            } else {
                self.knots[i].temperature
            };
            let start_c = self.on_segment(i, current)[0];
            let end_c = self.on_segment(i, end)[0];
            let mean_c = 0.5 * start_c + 0.5 * end_c;
            let width = (end - current).abs();
            let amount = mean_c * width;
            if amount.is_nan() || amount <= 0.0 {
                return Err(invalid(
                    "storage inverse encountered an unrepresentable segment",
                ));
            }
            let corrected = amount - correction;
            let next = cumulative + corrected;
            if target <= next {
                if target == next {
                    return Ok(end);
                }
                let fraction_energy = if amount.is_finite() {
                    (target - cumulative) / amount
                } else {
                    ((target - cumulative) / mean_c) / width
                };
                let scale = start_c.max(end_c);
                let a = start_c / scale;
                let b = end_c / scale;
                // A convex combination of endpoint squares avoids cancellation
                // for decreasing capacity near its supported cold/hot endpoint.
                let root = ((1.0 - fraction_energy) * a * a + fraction_energy * b * b).sqrt();
                let fraction = fraction_energy * (a + b) / (a + root);
                let temperature = current + fraction * (end - current);
                if !fraction.is_finite()
                    || !(0.0..=1.0).contains(&fraction)
                    || !temperature.is_finite()
                    || temperature == initial_k
                    || temperature < current.min(end)
                    || temperature > current.max(end)
                {
                    return Err(invalid(
                        "storage inverse is not representable at this temperature scale",
                    ));
                }
                return Ok(temperature);
            }
            correction = (next - cumulative) - corrected;
            cumulative = next;
            current = end;
            if current == boundary {
                break;
            }
            if heating {
                i += 1;
            } else {
                i -= 1;
            }
        }
        Err(invalid(
            "storage inverse exceeded its bounded source traversal",
        ))
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    #[must_use]
    pub const fn query_endpoints(&self) -> &[QueryPoint; 2] {
        &self.endpoints
    }
    /// Capacity and conductivity, with unmodified source answers and uncertainty.
    #[must_use]
    pub const fn source_lines(&self) -> &[PropertyLineSlice; 2] {
        &self.sources
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [
            self.knots[0].temperature,
            self.knots[self.knots.len() - 1].temperature,
        ]
    }
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.knots.iter().map(|k| k.temperature)
    }

    fn check_temperature(&self, temperature: f64) -> Result<(), ThermalTransportError> {
        let [low, high] = self.temperature_span_k();
        if !temperature.is_finite() || temperature < low || temperature > high {
            return Err(invalid(
                "temperature is outside the prepared volumetric interval",
            ));
        }
        Ok(())
    }
    fn segment(&self, temperature: f64) -> usize {
        self.knots
            .partition_point(|k| k.temperature <= temperature)
            .clamp(1, self.knots.len() - 1)
            - 1
    }
    fn on_segment(&self, index: usize, temperature: f64) -> [f64; 2] {
        let a = self.knots[index];
        let b = self.knots[index + 1];
        if temperature == a.temperature {
            return a.values;
        }
        if temperature == b.temperature {
            return b.values;
        }
        let fraction = (temperature - a.temperature) / (b.temperature - a.temperature);
        std::array::from_fn(|i| (1.0 - fraction) * a.values[i] + fraction * b.values[i])
    }
    fn values(&self, temperature: f64) -> Result<[f64; 2], ThermalTransportError> {
        self.check_temperature(temperature)?;
        Ok(self.on_segment(self.segment(temperature), temperature))
    }
    fn segment_storage(
        &self,
        index: usize,
        low: f64,
        high: f64,
    ) -> Result<f64, ThermalTransportError> {
        let mean = 0.5 * self.on_segment(index, low)[0] + 0.5 * self.on_segment(index, high)[0];
        let amount = mean * (high - low);
        if !positive(amount) {
            return Err(invalid(
                "volumetric storage increment is not positive and representable",
            ));
        }
        Ok(amount)
    }
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
