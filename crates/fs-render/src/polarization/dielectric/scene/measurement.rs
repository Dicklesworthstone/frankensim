//! Calibrated spectral receiving power from the native polarized scene path.
//!
//! Record EVERY launched sample, including occluded, rejected, dark and limited
//! paths. Each wavelength bin has its own launch count: conditioning on detector
//! hits would bias power upward. Input weights are W/nm, including the emitter's
//! angular quadrature or inverse-PDF weight. No cosine/inverse-square factor is
//! applied a second time after geometry has selected the receiving surface.
//!
//! All samples must belong to one consistently weighted estimator. To combine
//! sources, weight by their source-selection PDF or sum separately normalized
//! source readings. The existing midpoint observer quadrature supplies nm and
//! lm/W factors. Missing spectral bins are refused, not silently made black.
//!
//! Limited paths contribute no *detected* power yet and remain explicitly
//! unresolved in the report; this is not an error bound. Ideal analyzers act
//! per ray in the detector's projected frame, not through a finite-NA lens or
//! a solved angle-dependent sensor/polarizer. Monte Carlo uncertainty, omitted
//! wavelengths and geometry/model errors are not certified here.

use fs_blake3::ContentHash;
use fs_geom::Vec3;
use fs_math::det;
use crate::photometry::source::{
    BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, SPECTRAL_BINS, photopic_weight_nm, wavelength_nm,
};
use crate::polarization::Stokes;
use super::{Interaction, PathResult, PathTermination, PolarizedScene};

/// Measurement admission/range errors leave the accumulator unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeterError {
    NotADetector,
    InvalidBin,
    WavelengthMismatch,
    DetectorGeometryMismatch,
    InvalidAnalyzer,
    MissingSamples { bin: usize },
    UnresolvedPaths { count: u64 },
    NumericalRange,
}
impl core::fmt::Display for MeterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized detector: {self:?}")
    }
}
impl core::error::Error for MeterError {}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), MeterError> {
        let next = self.value + value;
        let correction = self.correction + if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else { (value - next) + self.value };
        if !next.is_finite() || !correction.is_finite() || !(next + correction).is_finite() {
            return Err(MeterError::NumericalRange);
        }
        self.value = next; self.correction = correction;
        Ok(())
    }
    fn total(self) -> f64 { self.value + self.correction }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Bin {
    stokes: [Sum; 4],
    launched: u64,
    detected: u64,
    unresolved: u64,
}

/// Fixed-size, allocation-free per-detector spectral accumulator.
/// Accumulation order is explicit; arbitrary thread merges are not bit-identical.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPowerMeter {
    object_id: u64,
    frame_identity: ContentHash,
    reference_axis: Vec3,
    bins: [Bin; SPECTRAL_BINS],
}

impl PolarizedScene<'_> {
    /// Bind a receiving meter to an admitted detector's actual geometry and axis.
    pub fn power_meter(&self, object_id: u64) -> Result<SpectralPowerMeter, MeterError> {
        let index = self.bindings.binary_search_by_key(&object_id, |b| b.object_id)
            .map_err(|_| MeterError::NotADetector)?;
        let Interaction::Detector { reference_axis } = self.bindings[index].interaction else {
            return Err(MeterError::NotADetector);
        };
        Ok(SpectralPowerMeter {
            object_id, reference_axis,
            frame_identity: self.geometry.instances()[index].frame_identity(),
            bins: [Bin::default(); SPECTRAL_BINS],
        })
    }
}

impl SpectralPowerMeter {
    /// Record one original launch outcome, not only successful detector hits.
    /// Wavelength must match the chosen native midpoint exactly. Re-recording a
    /// result is another statistical sample; IDs/deduplication remain caller-owned.
    pub fn record(&mut self, bin: usize, path: &PathResult) -> Result<(), MeterError> {
        if bin >= SPECTRAL_BINS { return Err(MeterError::InvalidBin); }
        if path.beam().wavelength_nm() != wavelength_nm(bin) {
            return Err(MeterError::WavelengthMismatch);
        }
        let mut next = self.bins[bin];
        // Keep the normalizer exactly representable as a binary64 integer.
        if next.launched >= (1_u64 << 53) { return Err(MeterError::NumericalRange); }
        next.launched += 1;
        match path.termination() {
            PathTermination::Detected => {
                let hit = path.last_hit().ok_or(MeterError::DetectorGeometryMismatch)?;
                if hit.object_id == self.object_id {
                    if hit.frame_identity != self.frame_identity {
                        return Err(MeterError::DetectorGeometryMismatch);
                    }
                    let beam = path.beam().with_reference_axis(self.reference_axis)
                        .map_err(|_| MeterError::NumericalRange)?;
                    let s = beam.power();
                    for (sum, value) in next.stokes.iter_mut().zip([s.i, s.q, s.u, s.v]) {
                        sum.add(value)?;
                    }
                    next.detected += 1;
                }
            }
            PathTermination::SurfaceLimit | PathTermination::DistanceLimit => {
                if path.beam().power().i > 0.0 { next.unresolved += 1; }
            }
            PathTermination::Opaque | PathTermination::RejectedMicrofacet | PathTermination::Extinguished => {}
        }
        self.bins[bin] = next;
        Ok(())
    }

    /// Mean detected spectral Stokes flux [W/nm], divided by ALL launches in
    /// this bin. A valid zero-hit bin is zero; an unobserved bin is an error.
    pub fn spectral_mean_w_nm(&self, bin: usize) -> Result<Stokes, MeterError> {
        let state = self.bins.get(bin).ok_or(MeterError::InvalidBin)?;
        if state.launched == 0 { return Err(MeterError::MissingSamples { bin }); }
        let s = state.stokes.map(|sum| sum.total() / state.launched as f64);
        stokes(s)
    }

    /// Integrate the original 380--780 nm midpoint quadrature. The result is
    /// receiving flux, not irradiance; divide by an explicitly supplied area
    /// to obtain W/m² or lux. Call `require_complete` to reject truncated paths.
    pub fn reading(&self) -> Result<DetectorReading, MeterError> {
        let mut watts = [Sum::default(); 4];
        let mut lumens = [Sum::default(); 4];
        let (mut launched, mut detected, mut unresolved) = (0_u64, 0_u64, 0_u64);
        for bin in 0..SPECTRAL_BINS {
            let s = self.spectral_mean_w_nm(bin)?;
            let photopic = LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(bin);
            for (axis, value) in [s.i, s.q, s.u, s.v].into_iter().enumerate() {
                watts[axis].add(value * BIN_WIDTH_NM)?;
                lumens[axis].add(value * photopic)?;
            }
            launched += self.bins[bin].launched;
            detected += self.bins[bin].detected;
            unresolved += self.bins[bin].unresolved;
        }
        Ok(DetectorReading {
            radiant: stokes(watts.map(Sum::total))?,
            photopic: stokes(lumens.map(Sum::total))?,
            launched, detected, unresolved,
        })
    }
}

/// Integrated detected contributions with explicit finite-path completeness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DetectorReading {
    radiant: Stokes,
    photopic: Stokes,
    launched: u64,
    detected: u64,
    unresolved: u64,
}
impl DetectorReading {
    pub const fn radiant_stokes_w(self) -> Stokes { self.radiant }
    pub const fn photopic_stokes_lm(self) -> Stokes { self.photopic }
    pub const fn launched_paths(self) -> u64 { self.launched }
    pub const fn detected_paths(self) -> u64 { self.detected }
    pub const fn unresolved_paths(self) -> u64 { self.unresolved }

    /// Reject a result with surviving truncated paths. Success is not an MC
    /// confidence interval, geometry certificate or continuum error bound.
    pub fn require_complete(self) -> Result<Self, MeterError> {
        if self.unresolved != 0 {
            Err(MeterError::UnresolvedPaths { count: self.unresolved })
        } else { Ok(self) }
    }

    /// Ideal linear analyzer rotated from the receiver's positive-Q axis.
    pub fn linear_analyzer(self, angle_rad: f64) -> Result<AnalyzerReading, MeterError> {
        if !angle_rad.is_finite() { return Err(MeterError::InvalidAnalyzer); }
        let (s, c) = det::sin_cos(2.0 * angle_rad.rem_euclid(core::f64::consts::PI));
        self.analyze([c, s, 0.0])
    }

    /// Ideal circular analyzer selecting positive or negative Stokes V.
    pub fn circular_analyzer(self, positive_v: bool) -> Result<AnalyzerReading, MeterError> {
        self.analyze([0.0, 0.0, if positive_v { 1.0 } else { -1.0 }])
    }
    fn analyze(self, axis: [f64; 3]) -> Result<AnalyzerReading, MeterError> {
        Ok(AnalyzerReading {
            radiant_w: project(self.radiant, axis)?,
            luminous_lm: project(self.photopic, axis)?,
            unresolved_paths: self.unresolved,
        })
    }
}

/// Analyzer throughput of the detected contributions; retains truncation status.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnalyzerReading {
    pub radiant_w: f64,
    pub luminous_lm: f64,
    pub unresolved_paths: u64,
}
fn project(s: Stokes, axis: [f64; 3]) -> Result<f64, MeterError> {
    // Halve before adding: (I+Q)/2 would overflow for a large valid I=Q.
    let value = 0.5 * s.i + 0.5 * s.q * axis[0] + 0.5 * s.u * axis[1] + 0.5 * s.v * axis[2];
    let tolerance = 32.0 * f64::EPSILON * s.i;
    if !value.is_finite() || value < -tolerance || value - s.i > tolerance {
        return Err(MeterError::NumericalRange);
    }
    Ok(value.clamp(0.0, s.i))
}
fn stokes(s: [f64; 4]) -> Result<Stokes, MeterError> {
    Stokes::try_new(s[0], s[1], s[2], s[3]).ok_or(MeterError::NumericalRange)
}

#[cfg(test)]
mod tests;
