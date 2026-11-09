//! Exact thermal refraction of the original tabulated material source.

use super::{
    AbsorptionTemperaturePathError, AbsorptionTemperatureSampleError, MaterialStatePointError,
    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
    SPECTRAL_OPTICAL_SAMPLE_COUNT, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
    SpectralDielectricAbsorptionStatePoint, SpectralDielectricAbsorptionTemperaturePath,
    SpectralTabulatedDielectricStatePoint,
};
use super::absorption::line_value;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialCard;
use std::sync::Arc;

/// Original wavelength, phase-index and intensity-loss claims along one typed
/// absolute-temperature line, with every other physical coordinate held fixed.
///
/// Preparation reuses the absorption owner's pinned source-line admission for
/// the COMPLETE original resolved bundle, including all phase and auxiliary
/// properties. No fit, resampling, index law or missing temperature dependence
/// is invented. Value queries retain piecewise source interpolation; smooth ray
/// jets additionally require an admitted interpolation cell at their wavelength.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralDielectricTemperaturePath {
    source: Arc<SpectralTabulatedDielectricStatePoint>,
    lines: SpectralDielectricAbsorptionTemperaturePath,
    maximum_wavelength_knot_span_m: f64,
    identity: ContentHash,
}

impl SpectralDielectricTemperaturePath {
    /// Admit original source support, every requirement and the native phase
    /// domain throughout the closed interval. Breakpoints of ALL properties
    /// count toward the source budget; cancellation publishes no partial path.
    pub fn try_new(
        card: &MaterialCard,
        source: SpectralTabulatedDielectricStatePoint,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
        max_source_knots: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AbsorptionTemperaturePathError> {
        let absorption = SpectralDielectricAbsorptionStatePoint::try_from_material(
            source.resolved().clone(),
        ).map_err(AbsorptionTemperaturePathError::Source)?;
        let lines = SpectralDielectricAbsorptionTemperaturePath::try_new(
            card, absorption, temperature_axis, temperature_range_k,
            max_source_knots, &mut cancelled,
        )?;
        for name in SPECTRAL_DIELECTRIC_IOR_PROPERTIES {
            let line = lines.source_line(name)
                .ok_or(AbsorptionTemperaturePathError::SourceStateMismatch)?;
            // Each phase value is affine between its original knots, so this
            // convex native domain is sufficient over the full interval, even
            // if the original requirement admitted any merely positive scalar.
            for &(_, value) in line {
                if cancelled() {
                    return Err(AbsorptionTemperaturePathError::Cancelled);
                }
                if !valid_phase(value) {
                    return Err(AbsorptionTemperaturePathError::Source(
                        MaterialStatePointError::InvalidDerived {
                            quantity: "spectral_dielectric_phase_index_temperature_path",
                        },
                    ));
                }
            }
        }
        // Each adjacent gap is affine between the union of its TWO original
        // wavelength-line breakpoint sets. These exact extrema let native
        // consumers retain their declared interpolation-gap ceiling as the
        // measured grid moves, without querying or resampling the source.
        let mut maximum_wavelength_knot_span_m: f64 = 0.0;
        for names in SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.windows(2) {
            let left = lines.source_line(names[0])
                .ok_or(AbsorptionTemperaturePathError::SourceStateMismatch)?;
            let right = lines.source_line(names[1])
                .ok_or(AbsorptionTemperaturePathError::SourceStateMismatch)?;
            for &(temperature_k, _) in left.iter().chain(right.iter()) {
                if cancelled() {
                    return Err(AbsorptionTemperaturePathError::Cancelled);
                }
                let gap = line_value(right, temperature_k)
                    .and_then(|hi| line_value(left, temperature_k).map(|lo| hi - lo))
                    .map_err(|source| AbsorptionTemperaturePathError::Spectrum {
                        temperature_k, source,
                    })?;
                if !gap.is_finite() || gap <= 0.0 {
                    return Err(AbsorptionTemperaturePathError::Spectrum {
                        temperature_k, source: AbsorptionTemperatureSampleError::InvalidSpectrum,
                    });
                }
                maximum_wavelength_knot_span_m = maximum_wavelength_knot_span_m.max(gap);
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.material.spectral-dielectric-temperature.v1");
        hash.update(source.identity().as_bytes());
        hash.update(lines.identity().as_bytes());
        if cancelled() {
            return Err(AbsorptionTemperaturePathError::Cancelled);
        }
        Ok(Self {
            source: Arc::new(source), lines, maximum_wavelength_knot_span_m,
            identity: hash.finalize(),
        })
    }

    pub fn source(&self) -> &SpectralTabulatedDielectricStatePoint { &self.source }
    pub fn identity(&self) -> ContentHash { self.identity }
    pub fn temperature_axis(&self) -> &str { self.lines.temperature_axis() }
    pub fn temperature_range_k(&self) -> [f64; 2] { self.lines.temperature_range_k() }
    pub fn source_knot_count(&self) -> usize { self.lines.source_knot_count() }

    /// Largest adjacent original wavelength gap anywhere in the full thermal
    /// interval, including interior source breakpoints. Consumers can preserve
    /// their interpolation-gap admission even when the measured grid moves.
    pub fn maximum_wavelength_knot_span_m(&self) -> f64 {
        self.maximum_wavelength_knot_span_m
    }

    /// Exact common support of the moving source grid. Disjoint global bands
    /// do not invalidate individual temperature-dependent value queries.
    pub fn common_wavelength_bounds_m(&self) -> Option<(f64, f64)> {
        self.lines.common_wavelength_bounds_m()
    }

    pub fn wavelength_bounds_m_at(
        &self, temperature_k: f64,
    ) -> Result<(f64, f64), SpectralDielectricTemperatureError> {
        self.lines.wavelength_bounds_m_at(temperature_k).map_err(Into::into)
    }

    pub fn absorption_per_m_at(
        &self, temperature_k: f64, wavelength_m: f64,
    ) -> Result<f64, SpectralDielectricTemperatureError> {
        self.lines.absorption_per_m_at(temperature_k, wavelength_m).map_err(Into::into)
    }

    /// Phase-only evaluation remains available at genuine source cusps. It
    /// does not claim a temperature or group derivative at those points.
    pub fn phase_index_at(
        &self, temperature_k: f64, wavelength_m: f64,
    ) -> Result<f64, SpectralDielectricTemperatureError> {
        let (lower, upper) = self.wavelength_bounds_m_at(temperature_k)?;
        if !wavelength_m.is_finite() || !(lower..=upper).contains(&wavelength_m) {
            return Err(SpectralDielectricTemperatureError::WavelengthOutsideRange);
        }
        let upper = self.bracket(temperature_k, wavelength_m)?;
        let left = self.value(SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[upper - 1], temperature_k)?;
        let right = self.value(SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[upper], temperature_k)?;
        let a = self.value(SPECTRAL_DIELECTRIC_IOR_PROPERTIES[upper - 1], temperature_k)?;
        let b = self.value(SPECTRAL_DIELECTRIC_IOR_PROPERTIES[upper], temperature_k)?;
        let phase = blend(a, b, (wavelength_m - left) / (right - left));
        if !valid_phase(phase) { return Err(SpectralDielectricTemperatureError::NumericalRange); }
        Ok(phase)
    }

    /// Prepare a monochromatic smooth thermal cell over the ENTIRE declared
    /// interval. The two active wavelength/index/loss lines must be affine
    /// throughout it, and the wavelength must stay in one source bracket.
    /// Exact equal-slope continuation is allowed; no tolerance smooths a cusp.
    /// A moving interior knot crossing or a genuine active temperature cusp
    /// refuses, even if the two extreme states happen to agree.
    ///
    /// At range endpoints derivatives are one-sided from inside the cell. A
    /// stationary interior spectral knot supports thermal jets, but group index
    /// is absent whenever its original left/right spectral slopes disagree.
    pub fn prepare(
        &self, wavelength_m: f64,
    ) -> Result<PreparedSpectralDielectricTemperatureCell<'_>, SpectralDielectricTemperatureError> {
        let (lo, hi) = self.common_wavelength_bounds_m()
            .ok_or(SpectralDielectricTemperatureError::WavelengthOutsideRange)?;
        if !wavelength_m.is_finite() || !(lo..=hi).contains(&wavelength_m) {
            return Err(SpectralDielectricTemperatureError::WavelengthOutsideRange);
        }
        let [minimum, maximum] = self.temperature_range_k();
        if minimum == maximum {
            return Err(SpectralDielectricTemperatureError::UndefinedDerivative);
        }
        let middle = minimum + 0.5 * (maximum - minimum);
        let upper = self.bracket(middle, wavelength_m)?;
        let affine = |names: &[&str; SPECTRAL_OPTICAL_SAMPLE_COUNT]| -> Result<_, SpectralDielectricTemperatureError> {
            Ok([
                Affine::try_new(self.line(names[upper - 1])?)?,
                Affine::try_new(self.line(names[upper])?)?,
            ])
        };
        let wavelength = affine(&SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES)?;
        let phase = affine(&SPECTRAL_DIELECTRIC_IOR_PROPERTIES)?;
        let absorption = affine(&SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES)?;
        for t in [minimum, maximum] {
            if wavelength_m < wavelength[0].value(t) || wavelength_m > wavelength[1].value(t) {
                return Err(SpectralDielectricTemperatureError::UndefinedDerivative);
            }
        }
        let neighbor = |index: usize| -> Result<_, SpectralDielectricTemperatureError> {
            Ok(Neighbor {
                wavelength: self.line(SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index])?,
                phase: self.line(SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index])?,
            })
        };
        Ok(PreparedSpectralDielectricTemperatureCell {
            path: self, wavelength_m, wavelength, phase, absorption,
            left: if upper > 1 { Some(neighbor(upper - 2)?) } else { None },
            right: if upper + 1 < SPECTRAL_OPTICAL_SAMPLE_COUNT {
                Some(neighbor(upper + 1)?)
            } else { None },
        })
    }

    fn line(&self, name: &str) -> Result<&[(f64, f64)], SpectralDielectricTemperatureError> {
        self.lines.source_line(name).ok_or(SpectralDielectricTemperatureError::InvalidSpectrum)
    }

    fn value(&self, name: &str, temperature_k: f64) -> Result<f64, SpectralDielectricTemperatureError> {
        line_value(self.line(name)?, temperature_k).map_err(Into::into)
    }

    fn bracket(&self, temperature_k: f64, wavelength_m: f64) -> Result<usize, SpectralDielectricTemperatureError> {
        for index in 1..SPECTRAL_OPTICAL_SAMPLE_COUNT {
            if self.value(SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index], temperature_k)? >= wavelength_m {
                return Ok(index);
            }
        }
        Err(SpectralDielectricTemperatureError::WavelengthOutsideRange)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Affine { range: [f64; 2], values: [f64; 2], slope: f64 }

impl Affine {
    fn try_new(line: &[(f64, f64)]) -> Result<Self, SpectralDielectricTemperatureError> {
        let (Some(&(lo, a)), Some(&(hi, b))) = (line.first(), line.last()) else {
            return Err(SpectralDielectricTemperatureError::InvalidSpectrum);
        };
        if lo >= hi { return Err(SpectralDielectricTemperatureError::UndefinedDerivative); }
        let slope = (b - a) / (hi - lo);
        if !slope.is_finite() { return Err(SpectralDielectricTemperatureError::NumericalRange); }
        if line.windows(2).any(|pair| (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0) != slope) {
            return Err(SpectralDielectricTemperatureError::UndefinedDerivative);
        }
        Ok(Self { range: [lo, hi], values: [a, b], slope })
    }

    fn value(self, temperature_k: f64) -> f64 {
        blend(self.values[0], self.values[1],
            (temperature_k - self.range[0]) / (self.range[1] - self.range[0]))
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Neighbor<'a> { wavelength: &'a [(f64, f64)], phase: &'a [(f64, f64)] }

/// Borrowed original interpolation cell. Hot sampling performs no card query,
/// allocation, fit or finite difference. Only exact interior spectral knots
/// additionally inspect the retained adjacent lines for a group derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedSpectralDielectricTemperatureCell<'a> {
    path: &'a SpectralDielectricTemperaturePath,
    wavelength_m: f64,
    wavelength: [Affine; 2],
    phase: [Affine; 2],
    absorption: [Affine; 2],
    left: Option<Neighbor<'a>>,
    right: Option<Neighbor<'a>>,
}

impl PreparedSpectralDielectricTemperatureCell<'_> {
    pub fn wavelength_m(&self) -> f64 { self.wavelength_m }
    pub fn temperature_range_k(&self) -> [f64; 2] { self.path.temperature_range_k() }

    /// Exact rational moving-grid derivatives of the declared source law.
    /// A finite nonpositive group index is retained as the source derivative;
    /// consumers separately decide whether it admits their causal clock.
    pub fn sample(&self, temperature_k: f64) -> Result<SpectralDielectricTemperatureSample, SpectralDielectricTemperatureError> {
        let [minimum, maximum] = self.temperature_range_k();
        if !temperature_k.is_finite() || !(minimum..=maximum).contains(&temperature_k) {
            return Err(SpectralDielectricTemperatureError::TemperatureOutsideRange);
        }
        let [left, right] = self.wavelength.map(|line| line.value(temperature_k));
        let [a, b] = self.phase.map(|line| line.value(temperature_k));
        let span = right - left;
        let q = (self.wavelength_m - left) / span;
        if !span.is_finite() || span <= 0.0 || !q.is_finite() || !(0.0..=1.0).contains(&q) {
            return Err(SpectralDielectricTemperatureError::NumericalRange);
        }
        let phase_index = blend(a, b, q);
        let [alpha_a, alpha_b] = self.absorption.map(|line| line.value(temperature_k));
        let absorption_per_m = blend(alpha_a, alpha_b, q);
        let spectral_slope = (b - a) / span;
        let width_derivative = self.wavelength[1].slope - self.wavelength[0].slope;
        let moving_coordinate = blend(self.wavelength[0].slope, self.wavelength[1].slope, q);
        let derivative_temperature_per_k = blend(self.phase[0].slope, self.phase[1].slope, q)
            - spectral_slope * moving_coordinate;
        let spectral_slope_derivative = ((self.phase[1].slope - self.phase[0].slope)
            - spectral_slope * width_derivative) / span;
        let second_derivative_temperature_per_k2 = -2.0 * moving_coordinate * spectral_slope_derivative;
        if !valid_phase(phase_index) || !absorption_per_m.is_finite() || absorption_per_m < 0.0
            || !derivative_temperature_per_k.is_finite()
            || !second_derivative_temperature_per_k2.is_finite() {
            return Err(SpectralDielectricTemperatureError::NumericalRange);
        }
        let mut group = phase_index - self.wavelength_m * spectral_slope;
        for (at_knot, neighbor, edge, index) in [
            (self.wavelength_m == left, self.left, left, a),
            (self.wavelength_m == right, self.right, right, b),
        ] {
            if at_knot && let Some(neighbor) = neighbor {
                let neighbor_wavelength = line_value(neighbor.wavelength, temperature_k)?;
                let neighbor_phase = line_value(neighbor.phase, temperature_k)?;
                let other_slope = (neighbor_phase - index) / (neighbor_wavelength - edge);
                if other_slope != spectral_slope { group = f64::NAN; }
            }
        }
        Ok(SpectralDielectricTemperatureSample {
            phase_index, group_index: group.is_finite().then_some(group), absorption_per_m,
            derivative_temperature_per_k, second_derivative_temperature_per_k2,
        })
    }
}

/// Optical value and exact temperature jet of one original smooth source cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDielectricTemperatureSample {
    pub phase_index: f64,
    pub group_index: Option<f64>,
    pub absorption_per_m: f64,
    pub derivative_temperature_per_k: f64,
    pub second_derivative_temperature_per_k2: f64,
}

/// Small allocation-free refusal used by native ray and thermal-field owners.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectralDielectricTemperatureError {
    TemperatureOutsideRange,
    WavelengthOutsideRange,
    UndefinedDerivative,
    NumericalRange,
    InvalidSpectrum,
}

impl From<AbsorptionTemperatureSampleError> for SpectralDielectricTemperatureError {
    fn from(error: AbsorptionTemperatureSampleError) -> Self {
        match error {
            AbsorptionTemperatureSampleError::TemperatureOutsideRange => Self::TemperatureOutsideRange,
            AbsorptionTemperatureSampleError::WavelengthOutsideRange => Self::WavelengthOutsideRange,
            AbsorptionTemperatureSampleError::NumericalRange => Self::NumericalRange,
            AbsorptionTemperatureSampleError::InvalidSpectrum => Self::InvalidSpectrum,
        }
    }
}

impl core::fmt::Display for SpectralDielectricTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral dielectric temperature cell refused: {self:?}")
    }
}
impl core::error::Error for SpectralDielectricTemperatureError {}

fn valid_phase(value: f64) -> bool {
    value.is_finite() && value >= 1.0 && value * value <= f64::MAX / 64.0
}

fn blend(a: f64, b: f64, fraction: f64) -> f64 {
    if fraction <= 0.5 { a + fraction * (b - a) }
    else { b + (1.0 - fraction) * (a - b) }
}

#[cfg(test)]
mod tests;
