//! Original material absorption along a prepared absolute-temperature path.

use super::{SpectralDielectricAbsorptionStatePoint, interpolate};
use crate::state_point::{
    MaterialStatePointError, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_OPTICAL_SAMPLE_COUNT, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
    SpectralDielectricAbsorptionSample,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
struct Prepared {
    source: SpectralDielectricAbsorptionStatePoint,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    // Exact source-line breakpoints and values, without the point-query receipt
    // payload in the hot search. The original source retains every requirement,
    // selected claim and immutable card identity, including auxiliary roles.
    lines: Vec<Vec<(f64, f64)>>,
    wavelengths: [usize; SPECTRAL_OPTICAL_SAMPLE_COUNT],
    absorption: [usize; SPECTRAL_OPTICAL_SAMPLE_COUNT],
    source_knot_count: usize,
    common_wavelength_bounds_m: Option<(f64, f64)>,
    identity: ContentHash,
}

/// Original alpha(T,lambda) and moving spectral support on one fixed condition.
///
/// Preparation uses the exact selected claims' continuous source lines, including
/// sampled-surface cell crossings, while pressure and every other coordinate
/// remain fixed. It retains all auxiliary requirements and refuses their source
/// gaps or constitutive-domain exits. No temperature fit, fixed wavelength grid,
/// extrapolation, or absorption/phase relation is introduced.
///
/// Clones share the prepared source. Temperature/wavelength queries perform
/// bounded searches and stack arithmetic only: no allocation or material query.
/// The nine native wavelengths may move independently with temperature.
/// Transparent spectra retain their finite support.
///
/// This is the original piecewise-affine source interpolation, not a measured
/// temperature law, an uncertainty enclosure, or a temperature-evolution model.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralDielectricAbsorptionTemperaturePath {
    prepared: Arc<Prepared>,
}

impl SpectralDielectricAbsorptionTemperaturePath {
    /// Admit the complete inclusive temperature interval before publishing it.
    ///
    /// The named axis must already mean static absolute temperature. The
    /// reference snapshot may lie outside the requested interval; every line
    /// still uses that snapshot's original claim and complete fixed context.
    /// A zero interval admits one temperature. All retained breakpoints across
    /// every requirement count toward max_source_knots, including duplicates
    /// shared by different properties. Each source line also retains matdb's
    /// own finite breakpoint cap.
    ///
    /// Cancellation is checked before and after each bounded line query and
    /// through breakpoint admission. Refusal leaves source and card unchanged.
    pub fn try_new(
        card: &MaterialCard,
        source: SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
        max_source_knots: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AbsorptionTemperaturePathError> {
        poll(&mut cancelled)?;
        if source.resolved().card_identity() != card.content_hash() {
            return Err(AbsorptionTemperaturePathError::SourceStateMismatch);
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if source.resolved().axis_quantities().get(temperature_axis) != Some(&absolute)
            || !source
                .resolved()
                .state_coordinate(temperature_axis)
                .is_some_and(|t| t.is_finite() && t > 0.0)
        {
            return Err(AbsorptionTemperaturePathError::TemperatureAxis);
        }
        let [lo, hi] = temperature_range_k;
        if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi < lo {
            return Err(AbsorptionTemperaturePathError::TemperatureRange);
        }
        let properties = source.resolved().properties();
        let minimum = properties.len() * if lo == hi { 1 } else { 2 };
        if minimum > max_source_knots {
            return Err(AbsorptionTemperaturePathError::Budget {
                required: minimum,
                maximum: max_source_knots,
            });
        }
        let context = source.resolved().query();
        let lower = context.clone().with(temperature_axis, lo).map_err(|error| {
            AbsorptionTemperaturePathError::Source(MaterialStatePointError::Query {
                property: format!("axis:{temperature_axis}"),
                source: error,
            })
        })?;
        let upper = context.with(temperature_axis, hi).map_err(|error| {
            AbsorptionTemperaturePathError::Source(MaterialStatePointError::Query {
                property: format!("axis:{temperature_axis}"),
                source: error,
            })
        })?;
        let mut lines = Vec::with_capacity(properties.len());
        let mut wavelengths = [usize::MAX; SPECTRAL_OPTICAL_SAMPLE_COUNT];
        let mut absorption = [usize::MAX; SPECTRAL_OPTICAL_SAMPLE_COUNT];
        let mut temperatures = Vec::with_capacity(minimum);
        let mut source_knot_count = 0usize;
        let mut hash = DomainHasher::new(
            "org.frankensim.material.spectral-absorption-temperature-path.v1",
        );
        hash.update(source.identity().as_bytes());
        hash.update(&(temperature_axis.len() as u64).to_le_bytes());
        hash.update(temperature_axis.as_bytes());
        for t in temperature_range_k {
            hash.update(&t.to_bits().to_le_bytes());
        }
        for (index, property) in properties.iter().enumerate() {
            poll(&mut cancelled)?;
            let requirement = property.requirement();
            let selected = property.answer().receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or(AbsorptionTemperaturePathError::SourceStateMismatch)?;
            let slice = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    temperature_axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|source| {
                    AbsorptionTemperaturePathError::Source(MaterialStatePointError::Support {
                        property: requirement.name().to_owned(),
                        source: Box::new(source),
                    })
                })?;
            poll(&mut cancelled)?;
            source_knot_count = source_knot_count
                .checked_add(slice.knots.len())
                .ok_or(AbsorptionTemperaturePathError::Budget {
                    required: usize::MAX,
                    maximum: max_source_knots,
                })?;
            if source_knot_count > max_source_knots {
                return Err(AbsorptionTemperaturePathError::Budget {
                    required: source_knot_count,
                    maximum: max_source_knots,
                });
            }
            let mut values = Vec::with_capacity(slice.knots.len());
            hash.update(&(slice.knots.len() as u64).to_le_bytes());
            for (temperature, answer) in slice.knots {
                poll(&mut cancelled)?;
                let value = answer.evidence.value.value;
                if !requirement.admissibility().admits(value) {
                    return Err(AbsorptionTemperaturePathError::Source(
                        MaterialStatePointError::OutsideConsumerDomain {
                            property: requirement.name().to_owned(),
                            value,
                            admissibility: requirement.admissibility(),
                        },
                    ));
                }
                hash.update(&temperature.to_bits().to_le_bytes());
                hash.update(&value.to_bits().to_le_bytes());
                temperatures.push(temperature);
                values.push((temperature, value));
            }
            if let Some(role) = SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES
                .iter()
                .position(|name| *name == requirement.name())
            {
                wavelengths[role] = index;
            }
            if let Some(role) = SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES
                .iter()
                .position(|name| *name == requirement.name())
            {
                absorption[role] = index;
            }
            lines.push(values);
        }
        if wavelengths.contains(&usize::MAX) || absorption.contains(&usize::MAX) {
            return Err(AbsorptionTemperaturePathError::SourceStateMismatch);
        }
        let mut prepared = Prepared {
            source,
            temperature_axis: temperature_axis.to_owned(),
            temperature_range_k,
            lines,
            wavelengths,
            absorption,
            source_knot_count,
            common_wavelength_bounds_m: None,
            identity: hash.finalize(),
        };
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        let (mut minimum, mut maximum) = (f64::NEG_INFINITY, f64::INFINITY);
        // Every wavelength, adjacent difference and alpha is affine between
        // union breakpoints. Their endpoint signs therefore admit the whole
        // continuous source interval, including interlaced property knots.
        for temperature_k in temperatures {
            poll(&mut cancelled)?;
            let samples = prepared.samples_at(temperature_k).map_err(|source| {
                AbsorptionTemperaturePathError::Spectrum {
                    temperature_k,
                    source,
                }
            })?;
            minimum = minimum.max(samples[0].wavelength_m);
            maximum = maximum.min(samples[SPECTRAL_OPTICAL_SAMPLE_COUNT - 1].wavelength_m);
        }
        if minimum <= maximum {
            prepared.common_wavelength_bounds_m = Some((minimum, maximum));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            prepared: Arc::new(prepared),
        })
    }

    pub fn source(&self) -> &SpectralDielectricAbsorptionStatePoint {
        &self.prepared.source
    }

    pub fn identity(&self) -> ContentHash {
        self.prepared.identity
    }

    pub fn temperature_axis(&self) -> &str {
        &self.prepared.temperature_axis
    }

    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.prepared.temperature_range_k
    }

    pub fn source_knot_count(&self) -> usize {
        self.prepared.source_knot_count
    }

    // The phase-temperature owner uses these SAME prepared source lines. It
    // neither re-queries claims nor copies a second interpolation/resolution
    // path; every auxiliary requirement remains in the original bundle.
    pub(in crate::state_point::spectral_optics) fn source_line(
        &self,
        name: &str,
    ) -> Option<&[(f64, f64)]> {
        self.prepared.source.resolved().properties().iter()
            .position(|property| property.requirement().name() == name)
            .map(|index| self.prepared.lines[index].as_slice())
    }

    /// Exact intersection of the original moving bands over the full path.
    /// None is permitted: each individual temperature may still have a valid
    /// spectrum. One shared endpoint is a usable monochromatic intersection.
    pub fn common_wavelength_bounds_m(&self) -> Option<(f64, f64)> {
        self.prepared.common_wavelength_bounds_m
    }

    pub fn wavelength_bounds_m_at(
        &self,
        temperature_k: f64,
    ) -> Result<(f64, f64), AbsorptionTemperatureSampleError> {
        let samples = self.prepared.samples_at(temperature_k)?;
        Ok((
            samples[0].wavelength_m,
            samples[SPECTRAL_OPTICAL_SAMPLE_COUNT - 1].wavelength_m,
        ))
    }

    /// Original internal intensity absorption [power Np/m] at (T,lambda).
    /// Both interpolation stages remain inside their original source support.
    pub fn absorption_per_m_at(
        &self,
        temperature_k: f64,
        wavelength_m: f64,
    ) -> Result<f64, AbsorptionTemperatureSampleError> {
        let samples = self.prepared.samples_at(temperature_k)?;
        if !wavelength_m.is_finite()
            || wavelength_m < samples[0].wavelength_m
            || wavelength_m > samples[SPECTRAL_OPTICAL_SAMPLE_COUNT - 1].wavelength_m
        {
            return Err(AbsorptionTemperatureSampleError::WavelengthOutsideRange);
        }
        interpolate(wavelength_m, |index| samples[index])
            .map_err(|_| AbsorptionTemperatureSampleError::NumericalRange)
    }
}

impl Prepared {
    fn samples_at(
        &self,
        temperature_k: f64,
    ) -> Result<[SpectralDielectricAbsorptionSample; SPECTRAL_OPTICAL_SAMPLE_COUNT],
                AbsorptionTemperatureSampleError> {
        let [lo, hi] = self.temperature_range_k;
        if !temperature_k.is_finite() || !(lo..=hi).contains(&temperature_k) {
            return Err(AbsorptionTemperatureSampleError::TemperatureOutsideRange);
        }
        let mut samples = [SpectralDielectricAbsorptionSample {
            wavelength_m: 0.0,
            absorption_per_m: 0.0,
        }; SPECTRAL_OPTICAL_SAMPLE_COUNT];
        for (i, sample) in samples.iter_mut().enumerate() {
            *sample = SpectralDielectricAbsorptionSample {
                wavelength_m: line_value(&self.lines[self.wavelengths[i]], temperature_k)?,
                absorption_per_m: line_value(&self.lines[self.absorption[i]], temperature_k)?,
            };
            if sample.wavelength_m <= 0.0 || sample.absorption_per_m < 0.0 {
                return Err(AbsorptionTemperatureSampleError::InvalidSpectrum);
            }
        }
        if samples.windows(2).any(|pair| pair[0].wavelength_m >= pair[1].wavelength_m) {
            return Err(AbsorptionTemperatureSampleError::InvalidSpectrum);
        }
        Ok(samples)
    }
}

pub(in crate::state_point::spectral_optics) fn line_value(
    knots: &[(f64, f64)],
    temperature_k: f64,
) -> Result<f64, AbsorptionTemperatureSampleError> {
    let right = knots.partition_point(|(t, _)| *t < temperature_k);
    let Some(&(hi, b)) = knots.get(right) else {
        return Err(AbsorptionTemperatureSampleError::TemperatureOutsideRange);
    };
    if hi == temperature_k {
        return Ok(b);
    }
    let Some(&(lo, a)) = right.checked_sub(1).and_then(|index| knots.get(index)) else {
        return Err(AbsorptionTemperatureSampleError::TemperatureOutsideRange);
    };
    let t = (temperature_k - lo) / (hi - lo);
    let value = if t <= 0.5 {
        a + t * (b - a)
    } else {
        b + (1.0 - t) * (a - b)
    };
    if !value.is_finite() {
        return Err(AbsorptionTemperatureSampleError::NumericalRange);
    }
    Ok(value)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AbsorptionTemperaturePathError> {
    if cancelled() {
        Err(AbsorptionTemperaturePathError::Cancelled)
    } else {
        Ok(())
    }
}

/// Construction refuses before any prepared path is published.
#[derive(Clone, Debug, PartialEq)]
pub enum AbsorptionTemperaturePathError {
    SourceStateMismatch,
    TemperatureAxis,
    TemperatureRange,
    Source(MaterialStatePointError),
    Budget { required: usize, maximum: usize },
    Cancelled,
    Spectrum { temperature_k: f64, source: AbsorptionTemperatureSampleError },
}

/// Allocation-free refusal from a prepared temperature/wavelength query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbsorptionTemperatureSampleError {
    TemperatureOutsideRange,
    WavelengthOutsideRange,
    NumericalRange,
    InvalidSpectrum,
}

impl core::fmt::Display for AbsorptionTemperaturePathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral absorption temperature path refused: {self:?}")
    }
}

impl core::error::Error for AbsorptionTemperaturePathError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Source(source) => Some(source),
            Self::Spectrum { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl core::fmt::Display for AbsorptionTemperatureSampleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral absorption temperature sample refused: {self:?}")
    }
}

impl core::error::Error for AbsorptionTemperatureSampleError {}

#[cfg(test)]
mod tests;
