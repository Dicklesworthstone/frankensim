//! Conditional higher-order VOLUME scattering, independent of camera sampling.
//!
//! For each wavelength, external point fixtures give a bounded irradiance
//! J <= sum_i I_peak_i / distance(source_i, grid)^2. The original phase and
//! finite-grid chord envelopes give first-order radiance B1 <= J*P*h. If the
//! original positive volume operator has norm q<1, orders two and above are
//! bounded by B1*q/(1-q). The bound is global, not fitted to a sampled ray.
//! This is the Neumann-series remainder of the equation of transfer (PBRT 4e,
//! Light Transport II, The Equation of Transfer). Surface reflection and all
//! non-fixture emission are explicitly excluded; opaque meshes only absorb.

use super::*;
use crate::photometry::source::SPECTRAL_BINS;
use crate::volumes::transport::VolumeTransportError;

/// Version of the external-point-source and positive-volume-order tail proof.
pub const GRID_VOLUME_REMAINDER_SEMANTICS_VERSION: u32 = 1;

/// Global omitted-volume-order envelope for an admitted original grid and rig.
/// It excludes paths with a surface interaction, thermal emission, calibration,
/// construction and discretization errors. It is not total image accuracy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridVolumeRemainder {
    upper_cd_m2: f64,
    first_order_upper_cd_m2: f64,
    max_contraction: f64,
    min_source_distance_m: Option<f64>,
    grid_diameter_m: f64,
}
impl GridVolumeRemainder {
    /// Derivation version, separate from the unchanged single-order target.
    #[must_use]
    pub const fn semantics_version(self) -> u32 { GRID_VOLUME_REMAINDER_SEMANTICS_VERSION }
    /// Upper bound on volume-only scattering orders two and above, in cd/m².
    #[must_use]
    pub const fn upper_cd_m2(self) -> f64 { self.upper_cd_m2 }
    /// Global first-order envelope used in the proof, not a camera estimate.
    #[must_use]
    pub const fn first_order_upper_cd_m2(self) -> f64 { self.first_order_upper_cd_m2 }
    /// Largest retained contraction in an illuminated, scattering spectral lane.
    #[must_use]
    pub const fn max_contraction(self) -> f64 { self.max_contraction }
    /// Proved source-to-grid separation lower bound; None when no active source
    /// requires separation (for example an exactly non-scattering field).
    #[must_use]
    pub const fn min_source_distance_m(self) -> Option<f64> { self.min_source_distance_m }
    /// Original grid chord-length ceiling, in metres.
    #[must_use]
    pub const fn grid_diameter_m(self) -> f64 { self.grid_diameter_m }
}

impl SingleScatterReference<'_, '_> {
    /// Optionally prepare a global bound on omitted volume scattering once.
    /// `max_cells` caps scanned cells and retained phase-model slots; the scan
    /// reuses the original medium, optical maps and source spectral packets.
    /// It does not affect single-scattering samples, refinement or tolerances.
    ///
    /// Unknown is preserved for homogeneous infinite media, thermal sources,
    /// a point fixture inside/touching the support box, noncontracting lit
    /// wavelengths, or nonrepresentable envelopes. Cancellation/work exhaustion
    /// return an error, never a partial bound. Finite support may establish
    /// contraction even without absorption; dense conservative media may fail.
    ///
    /// The bound covers only paths with volume events and opaque absorption.
    /// It never includes surface reflection, additive emission or model error.
    pub fn with_volume_remainder(mut self, max_cells: usize, cx: &Cx<'_>)
        -> Result<Self, ScatterEnclosureError> {
        cx.checkpoint()?;
        if max_cells == 0 { return Err(ScatterEnclosureError::Budget); }
        self.volume_remainder = match &self.medium {
            ScatteringMedium::Grid(grid) => grid.volume_remainder(self.scene, max_cells, cx)?,
            ScatteringMedium::Gray { .. } => None,
        };
        cx.checkpoint()?;
        Ok(self)
    }

    /// The optional prepared envelope, independent of any accepted camera ray.
    #[must_use]
    pub const fn volume_remainder(&self) -> Option<GridVolumeRemainder> { self.volume_remainder }
}

impl grid::GridScattering<'_> {
    pub(super) fn volume_remainder(&self, scene: &PhotometricScene<'_>, max_cells: usize, cx: &Cx<'_>)
        -> Result<Option<GridVolumeRemainder>, ScatterEnclosureError> {
        let envelope = match self.packet.volume_contraction_envelope(cx, max_cells) {
            Ok(Some(value)) => value,
            Ok(None) => return Ok(None),
            Err(VolumeTransportError::Cancelled) => return Err(PhotometryError::Cancelled.into()),
            Err(VolumeTransportError::Budget) => return Err(ScatterEnclosureError::Budget),
            Err(_) => return Ok(None),
        };
        let mut result = GridVolumeRemainder { upper_cd_m2: 0.0, first_order_upper_cd_m2: 0.0,
            max_contraction: 0.0, min_source_distance_m: None, grid_diameter_m: envelope.diameter_m() };
        if !self.packet.has_scattering() { return Ok(Some(result)); }
        let support = envelope.bounds_m();
        let mut irradiance = [zero(); SPECTRAL_BINS];
        for (index, source) in scene.sources.iter().enumerate() {
            cx.checkpoint()?;
            if source.flux_lm() == 0.0 { continue; }
            let position = coordinates(source.position_m());
            let gap = std::array::from_fn(|i| {
                let low = (point(support[i].lo()) - position[i]).lo().max(0.0);
                let high = (position[i] - point(support[i].hi())).lo().max(0.0);
                point(low.max(high))
            });
            let Some(distance) = norm(gap) else { return Ok(None); };
            if distance.lo() <= 0.0 { return Ok(None); }
            let distance = distance.lo();
            result.min_source_distance_m = Some(result.min_source_distance_m.map_or(distance, |old| old.min(distance)));
            // The original admitted maximum bounds every represented angular
            // law, including arbitrary asymmetric tables; no new beam fit.
            let intensity = mul(point(source.flux_lm()), point(source.distribution.peak_per_steradian()));
            let incident = intensity / point(distance) / point(distance);
            if !finite(incident) { return Ok(None); }
            for lane in 0..SPECTRAL_BINS {
                if lane % 16 == 0 { cx.checkpoint()?; }
                irradiance[lane] = add(irradiance[lane], mul(incident, self.weights[index][lane]));
                if !finite(irradiance[lane]) { return Ok(None); }
            }
        }
        let mut first = zero();
        let mut tail = zero();
        for lane in 0..SPECTRAL_BINS {
            cx.checkpoint()?;
            let h = envelope.scattering_integral_upper()[lane];
            let peak = envelope.phase_peak_per_sr()[lane];
            let b1 = mul(mul(irradiance[lane], point(h)), point(peak));
            if !finite(b1) { return Ok(None); }
            if b1.hi() == 0.0 { continue; }
            let Some(upper) = envelope.tail_upper(lane, b1.hi()) else { return Ok(None); };
            first = add(first, b1);
            tail = add(tail, point(upper));
            if !finite(first) || !finite(tail) { return Ok(None); }
            result.max_contraction = result.max_contraction.max(envelope.contraction_upper()[lane]);
        }
        result.first_order_upper_cd_m2 = first.hi();
        result.upper_cd_m2 = tail.hi();
        cx.checkpoint()?;
        Ok(Some(result))
    }
}

impl EnclosedSingleScatter {
    /// Enclose ALL volume-scattering orders contributed along this same finite
    /// segment, when the optional global remainder is available. Only the upper
    /// endpoint grows: omitted orders are nonnegative. This is not surface
    /// radiance, direct source visibility, thermal light or total physical error.
    /// Unknown remainder and overflow return None, never the single-order bound.
    #[must_use]
    pub fn all_volume_orders_enclosure(self) -> Option<NumericalCertificate> {
        let tail = self.omitted_multiple_scattering_bound_cd_m2?;
        if !tail.is_finite() || tail < 0.0
            || !self.numerical.lo.is_finite() || !self.numerical.hi.is_finite()
            || self.numerical.lo < 0.0 || self.numerical.lo > self.numerical.hi
            || self.numerical.kind != fs_evidence::NumericalKind::Enclosure { return None; }
        let upper = if tail == 0.0 { self.numerical.hi } else {
            (point(self.numerical.hi)+point(tail)).hi()
        };
        upper.is_finite().then(|| NumericalCertificate::enclosure(self.numerical.lo,upper))
    }
}

#[cfg(test)]
mod tests;
