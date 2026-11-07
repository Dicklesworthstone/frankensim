//! Physical source Stokes fractions in an original emitting-source frame.
//!
//! The intensity spectrum stays with `SurfaceEmitter`. These coordinates are
//! its normalized (I=1) polarization at each vacuum wavelength, expressed in
//! the frame obtained by projecting the body-local reference axis transverse
//! to the actual outward light direction. Proper rigid placement precedes
//! that projection. Environment bindings instead use map-local coordinates
//! and the original map rotation. A parallel axis refuses for a polarized query.
//!
//! This is incoherent Stokes emission, not a coherent Jones amplitude or a
//! polarization law inferred from an RGB color or material name. Convex
//! interpolation of normalized Stokes states is an explicit source model;
//! it preserves the physical cone without fitting angles or extrapolating.
//! Frame convention: <https://mitsuba.readthedocs.io/en/latest/src/key_topics/polarization.html>.

use super::{CameraRadianceError, MAX_RADIANCE_KNOTS};
use crate::polarization::Stokes;
use fs_geom::Vec3;
use fs_math::det;

/// Authored constant or spectral normalized emission polarization. The
/// positive-Q axis is in the original geometry's local coordinates for a
/// surface source, or the original map coordinates for an environment. Both
/// emitting sides use the same declared Stokes coordinates in their own
/// right-handed outward travel frames; in particular V follows that frame.
/// Area-bundle preparation retains this declared frame even at unpolarized
/// spectral knots, and requires it to be transverse to each initial seed.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePolarization {
    reference_axis: Vec3,
    spectrum: PolarizationSpectrum,
}

#[derive(Clone, Debug, PartialEq)]
enum PolarizationSpectrum {
    Constant(Stokes),
    Tabulated(Vec<(f64, Stokes)>),
}

impl SourcePolarization {
    /// Declare wavelength-independent fractions. No intensity normalization
    /// is performed: I must be exactly one and |(Q,U,V)| must not exceed one.
    pub fn try_constant(
        reference_axis_local: Vec3,
        normalized_stokes: Stokes,
    ) -> Result<Self, CameraRadianceError> {
        validate_state(normalized_stokes)?;
        Ok(Self {
            reference_axis: admit_axis(reference_axis_local)?,
            spectrum: PolarizationSpectrum::Constant(normalized_stokes),
        })
    }

    /// Retain 2..=MAX_RADIANCE_KNOTS ordered positive vacuum wavelengths in
    /// nanometres. Fractional Q/U/V interpolate linearly between original
    /// knots, independently of the emitter's intensity spectrum. A zero or
    /// reversed interval, nonphysical state or nonfinite coordinate refuses.
    pub fn try_spectral(
        reference_axis_local: Vec3,
        knots: Vec<(f64, Stokes)>,
    ) -> Result<Self, CameraRadianceError> {
        if !(2..=MAX_RADIANCE_KNOTS).contains(&knots.len())
            || knots.iter().any(|(w, _)| !w.is_finite() || *w <= 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(CameraRadianceError::InvalidPolarization);
        }
        for &(_, state) in &knots {
            validate_state(state)?;
        }
        Ok(Self {
            reference_axis: admit_axis(reference_axis_local)?,
            spectrum: PolarizationSpectrum::Tabulated(knots),
        })
    }

    /// Unit reference axis before the native source's proper-rigid placement
    /// or environment-map rotation.
    pub const fn reference_axis(&self) -> Vec3 {
        self.reference_axis
    }

    /// Explicit spectral support; a constant declaration has no finite band.
    pub fn wavelength_bounds_nm(&self) -> Option<(f64, f64)> {
        match &self.spectrum {
            PolarizationSpectrum::Constant(_) => None,
            PolarizationSpectrum::Tabulated(knots) => Some((knots[0].0, knots[knots.len() - 1].0)),
        }
    }

    /// The normalized source state in its own projected outward frame. This
    /// bounded logarithmic lookup allocates nothing and retains exact knots.
    pub fn normalized_stokes(&self, wavelength_nm: f64) -> Result<Stokes, CameraRadianceError> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(CameraRadianceError::InvalidPolarization);
        }
        match &self.spectrum {
            PolarizationSpectrum::Constant(state) => Ok(*state),
            PolarizationSpectrum::Tabulated(knots) => {
                match knots.binary_search_by(|(w, _)| w.total_cmp(&wavelength_nm)) {
                    Ok(index) => Ok(knots[index].1),
                    Err(index) if index > 0 && index < knots.len() => {
                        let (lo, a) = knots[index - 1];
                        let (hi, b) = knots[index];
                        let fraction = (wavelength_nm - lo) / (hi - lo);
                        let blend = |a: f64, b: f64| (1.0 - fraction) * a + fraction * b;
                        Ok(Stokes {
                            i: 1.0,
                            q: blend(a.q, b.q),
                            u: blend(a.u, b.u),
                            v: blend(a.v, b.v),
                        })
                    }
                    Err(_) => Err(CameraRadianceError::InvalidPolarization),
                }
            }
        }
    }
}

fn validate_state(state: Stokes) -> Result<(), CameraRadianceError> {
    if state.i != 1.0
        || ![state.q, state.u, state.v].into_iter().all(f64::is_finite)
        || det::hypot(det::hypot(state.q, state.u), state.v) > 1.0
    {
        return Err(CameraRadianceError::InvalidPolarization);
    }
    Ok(())
}

fn admit_axis(axis: Vec3) -> Result<Vec3, CameraRadianceError> {
    if ![axis.x, axis.y, axis.z].into_iter().all(f64::is_finite) {
        return Err(CameraRadianceError::InvalidPolarization);
    }
    let scale = axis.x.abs().max(axis.y.abs()).max(axis.z.abs());
    if scale == 0.0 {
        return Err(CameraRadianceError::InvalidPolarization);
    }
    // Divide components, not a common reciprocal: a subnormal authored axis
    // is a valid direction, and its reciprocal need not be representable.
    let scaled = Vec3::new(axis.x / scale, axis.y / scale, axis.z / scale);
    let length = det::hypot(det::hypot(scaled.x, scaled.y), scaled.z);
    Ok(scaled.scale(1.0 / length))
}

#[cfg(test)]
mod tests;
