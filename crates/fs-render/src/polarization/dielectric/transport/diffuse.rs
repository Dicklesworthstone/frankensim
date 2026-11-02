//! Explicit ideal depolarizing Lambertian reflection in the native beam owner.
//!
//! The Mueller BRDF has only `M00 = reflectance / pi`. Cosine sampling cancels
//! its projected-solid-angle factor, leaving exactly `reflectance * I` power.
//! This reciprocal model also applies to an unpolarized camera estimator. It
//! does not infer polarization retention, subsurface transport, or a material
//! spectrum: the caller supplies the reflectance at the beam's wavelength.

use super::{BeamError, PolarizationFrame, PolarizedBeam, cross, unit};
use crate::polarization::Stokes;
use fs_geom::Vec3;
use fs_math::det;

/// A continuous reflection proposal with density in inverse steradians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseBeamSample {
    /// Reflected, fully depolarized beam in its original surrounding medium.
    pub beam: PolarizedBeam,
    /// Cosine-weighted hemisphere density relative to the geometric normal.
    pub pdf_solid_angle: f64,
}

impl PolarizedBeam {
    /// Sample a two-sided ideal depolarizing Lambertian reflector.
    ///
    /// `normal` is the geometric normal; its sign does not select a different
    /// material. Both random coordinates must lie in `[0, 1)`. The continuation
    /// retains the exact position, wavelength, and homogeneous or sampled
    /// gradient medium. No extra refractive-index ratio, roulette division,
    /// ray-origin offset, or polarization frame conversion is applied.
    ///
    /// Exact zero reflectance gives zero power. Positive power that underflows
    /// when multiplied by positive reflectance is an explicit range error.
    pub fn sample_depolarizing_diffuse(
        self,
        normal: Vec3,
        reflectance: f64,
        samples: [f64; 2],
    ) -> Result<DiffuseBeamSample, BeamError> {
        if !(0.0..=1.0).contains(&reflectance) || samples.iter().any(|u| !(0.0..1.0).contains(u)) {
            return Err(BeamError::InvalidDiffuseSample);
        }
        let mut normal = unit(normal)?;
        let incident_cosine = normal.dot(self.frame.direction);
        if !incident_cosine.is_finite() || incident_cosine == 0.0 {
            return Err(BeamError::InvalidFrame);
        }
        if incident_cosine > 0.0 {
            normal = normal.scale(-1.0);
        }
        // The least-aligned Cartesian axis bounds the cross-product away
        // from zero, including both poles, without a near-pole approximation.
        let tangent = unit(cross(reference_axis(normal), normal))?;
        let bitangent = cross(normal, tangent);
        let radius = det::sqrt(samples[0]);
        let cosine = det::sqrt(1.0 - samples[0]);
        let (sine_phi, cosine_phi) = det::sin_cos(core::f64::consts::TAU * samples[1]);
        let x = radius * cosine_phi;
        let y = radius * sine_phi;
        let direction = unit(Vec3::new(
            tangent.x * x + bitangent.x * y + normal.x * cosine,
            tangent.y * x + bitangent.y * y + normal.y * cosine,
            tangent.z * x + bitangent.z * y + normal.z * cosine,
        ))?;
        let outgoing_cosine = normal.dot(direction);
        if !outgoing_cosine.is_finite() || outgoing_cosine <= 0.0 {
            return Err(BeamError::NumericalRange);
        }
        let power = self.power.i * reflectance;
        if !power.is_finite() || (self.power.i > 0.0 && reflectance > 0.0 && power == 0.0) {
            return Err(BeamError::NumericalRange);
        }
        let frame = PolarizationFrame::try_new(direction, reference_axis(direction))?;
        Ok(DiffuseBeamSample {
            beam: Self {
                frame,
                power: Stokes::unpolarized(power),
                ..self
            },
            pdf_solid_angle: outgoing_cosine / core::f64::consts::PI,
        })
    }
}

fn reference_axis(direction: Vec3) -> Vec3 {
    let x = direction.x.abs();
    let y = direction.y.abs();
    let z = direction.z.abs();
    if x <= y && x <= z {
        Vec3::new(1.0, 0.0, 0.0)
    } else if y <= z {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

#[cfg(test)]
mod tests;
