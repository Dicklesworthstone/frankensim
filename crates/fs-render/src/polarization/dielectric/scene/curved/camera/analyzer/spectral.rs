//! Source-owned spectral intensity rows for a passive polarization detector.
//!
//! Each row acts on the arriving Stokes radiance before spectral integration.
//! Interpolating absolute Mueller rows retains both chromatic throughput and
//! polarization sensitivity; interpolating angles or normalized rows would not.
//! The common reference axis is projected by the original launch/lens owner.
//! See NIST, *Accurate Spectral Characterization of Polarization-Dependent Loss*:
//! <https://www.nist.gov/publications/accurate-spectral-characterization-polarization-dependent-loss-0>.

use super::{CameraAnalyzer, CameraRadianceError, Stokes};
use crate::polarization::dielectric::scene::curved::camera::MAX_RADIANCE_KNOTS;
use fs_geom::Vec3;
use fs_math::det;

/// An immutable, finite-band passive detector law in one declared world frame.
///
/// Knots are positive vacuum wavelengths in nanometres and absolute Mueller
/// intensity rows `[a0,aQ,aU,aV]`, not normalized Stokes states. Their convex
/// interpolation is the declared model; measured accuracy, angular response,
/// finite polarizer thickness and group delay are not inferred from it.
/// All-zero rows retain the native camera's explicit refusal policy.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralCameraAnalyzer {
    reference_axis: Vec3,
    knots: Vec<(f64, Stokes)>,
}

impl SpectralCameraAnalyzer {
    /// Retain 2..=MAX_RADIANCE_KNOTS strictly increasing original knots. Every
    /// row must be passive for ALL physical inputs, including circular ones:
    /// `0 < a0 <= 1`, `|(aQ,aU,aV)| <= min(a0,1-a0)`.
    pub fn try_new(
        reference_axis: Vec3,
        knots: Vec<(f64, Stokes)>,
    ) -> Result<Self, CameraRadianceError> {
        if !(2..=MAX_RADIANCE_KNOTS).contains(&knots.len())
            || knots.iter().any(|(w, _)| !w.is_finite() || *w <= 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(CameraRadianceError::InvalidPolarization);
        }
        for &(_, row) in &knots {
            CameraAnalyzer::try_new(reference_axis, row)
                .map_err(|_| CameraRadianceError::InvalidPolarization)?;
        }
        Ok(Self {
            reference_axis,
            knots,
        })
    }

    pub const fn reference_axis(&self) -> Vec3 {
        self.reference_axis
    }

    pub fn knots(&self) -> &[(f64, Stokes)] {
        &self.knots
    }

    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }

    /// Resolve once at the actual sampled wavelength, before lens or scene
    /// traversal. Exact authored knots are retained bit for bit. Queries
    /// allocate nothing and never extrapolate or renormalize throughput.
    pub fn at(&self, wavelength_nm: f64) -> Result<CameraAnalyzer, CameraRadianceError> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(CameraRadianceError::InvalidPolarization);
        }
        let row = match self
            .knots
            .binary_search_by(|(w, _)| w.total_cmp(&wavelength_nm))
        {
            Ok(index) => self.knots[index].1,
            Err(index) if index > 0 && index < self.knots.len() => {
                let (lo, a) = self.knots[index - 1];
                let (hi, b) = self.knots[index];
                let fraction = (wavelength_nm - lo) / (hi - lo);
                // Anchor at the smaller component. Equal knots retain their
                // exact bits, including the smallest positive throughput.
                let blend = |a: f64, b: f64| {
                    if a == b {
                        a
                    } else if a < b {
                        a + fraction * (b - a)
                    } else {
                        b + (1.0 - fraction) * (a - b)
                    }
                };
                let mut row = Stokes {
                    i: blend(a.i, b.i),
                    q: blend(a.q, b.q),
                    u: blend(a.u, b.u),
                    v: blend(a.v, b.v),
                };
                // The passive-row set is convex. Only generated interpolation
                // roundoff can leave it; direct that roundoff inward without
                // changing a0 or repairing authored inputs/knots.
                let bound = row.i.min(1.0 - row.i);
                let degree = det::hypot(det::hypot(row.q, row.u), row.v);
                if degree > bound {
                    let scale = bound / degree;
                    row.q *= scale;
                    row.u *= scale;
                    row.v *= scale;
                    for _ in 0..4 {
                        if det::hypot(det::hypot(row.q, row.u), row.v) <= bound {
                            break;
                        }
                        let value = if row.q.abs() >= row.u.abs() && row.q.abs() >= row.v.abs() {
                            &mut row.q
                        } else if row.u.abs() >= row.v.abs() {
                            &mut row.u
                        } else {
                            &mut row.v
                        };
                        *value = value.abs().next_down().copysign(*value);
                    }
                }
                row
            }
            Err(_) => return Err(CameraRadianceError::InvalidPolarization),
        };
        CameraAnalyzer::try_new(self.reference_axis, row)
            .map_err(|_| CameraRadianceError::InvalidPolarization)
    }

    /// Require the whole integration band before drawing/tracing a pixel.
    /// A detector law cannot borrow missing spectral support from scene or QE.
    pub fn admit_band(&self, min_nm: f64, max_nm: f64) -> Result<(), CameraRadianceError> {
        let (lo, hi) = self.wavelength_bounds_nm();
        if !min_nm.is_finite()
            || !max_nm.is_finite()
            || min_nm > max_nm
            || min_nm < lo
            || max_nm > hi
        {
            return Err(CameraRadianceError::InvalidPolarization);
        }
        Ok(())
    }
}

// Prepared consumers borrow one immutable response. The original constant
// branch performs no spectral lookup and returns the same small Copy value.
#[derive(Clone, Copy)]
pub(crate) enum AnalyzerResponse<'a> {
    Constant(CameraAnalyzer),
    Spectral(&'a SpectralCameraAnalyzer),
}

impl AnalyzerResponse<'_> {
    pub(crate) fn at(self, wavelength_nm: f64) -> Result<CameraAnalyzer, CameraRadianceError> {
        match self {
            Self::Constant(analyzer) => Ok(analyzer),
            Self::Spectral(spectrum) => spectrum.at(wavelength_nm),
        }
    }

    pub(crate) fn admit_band(self, min_nm: f64, max_nm: f64) -> Result<(), CameraRadianceError> {
        match self {
            Self::Constant(_) => Ok(()),
            Self::Spectral(spectrum) => spectrum.admit_band(min_nm, max_nm),
        }
    }
}

#[cfg(test)]
mod tests;
