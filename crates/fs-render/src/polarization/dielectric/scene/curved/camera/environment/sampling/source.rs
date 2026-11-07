//! One original texel supplies radiance, polarization and its angular PDF.

use super::super::super::{CameraRadianceError, LAMBDA_MAX, LAMBDA_MIN, finite_nonnegative};
use super::{EnvironmentCameraError, EnvironmentRadiance};
use crate::dielectric::DielectricError;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::CurvedSceneError;
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_geom::Vec3;

#[derive(Clone, Copy, Debug)]
pub(super) struct EnvironmentSourceSample {
    pub(super) direction: Vec3,
    pub(super) radiance: f64,
    pub(super) pdf: f64,
    pub(super) frame: PolarizationFrame,
    pub(super) stokes: Stokes,
}

impl EnvironmentRadiance<'_> {
    pub(super) fn sample_direction(
        &self,
        wavelength: f64,
        sample: [f64; 2],
    ) -> Result<Option<EnvironmentSourceSample>, EnvironmentCameraError> {
        if !wavelength.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength) {
            return Err(CameraRadianceError::Optics(DielectricError::InvalidWavelength).into());
        }
        if sample
            .iter()
            .any(|value| !value.is_finite() || !(0.0..1.0).contains(value))
        {
            return Err(EnvironmentCameraError::InvalidSample);
        }
        if self.scale_w_m2_sr_nm == 0.0 || self.map.is_black() {
            return Ok(None);
        }
        let (original, texel) = self
            .map
            .sample_with_texel(sample[0], sample[1])
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        let radiance = finite_nonnegative(
            original.emission.0.eval(wavelength) * original.emission.1 * self.scale_w_m2_sr_nm,
        )?;
        if radiance == 0.0 {
            return Ok(None);
        }
        let incoming = original.direction.scale(-1.0);
        let (stokes, axis) = if let Some(polarization) = self.polarization {
            let law = polarization
                .at(texel)
                .ok_or(CameraRadianceError::InvalidPolarization)?;
            let state = law.normalized_stokes(wavelength)?;
            let axis = if state.q == 0.0 && state.u == 0.0 {
                stable_axis(incoming)
            } else {
                polarization.world_axis(law.reference_axis())
            };
            (state, axis)
        } else {
            (Stokes::unpolarized(1.0), stable_axis(incoming))
        };
        let frame = PolarizationFrame::try_new(incoming, axis)
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?;
        Ok(Some(EnvironmentSourceSample {
            direction: frame.direction().scale(-1.0),
            radiance,
            pdf: original.pdf_solid_angle,
            frame,
            stokes,
        }))
    }
}

fn stable_axis(direction: Vec3) -> Vec3 {
    let [x, y, z] = [direction.x.abs(), direction.y.abs(), direction.z.abs()];
    if x <= y && x <= z {
        Vec3::new(1.0, 0.0, 0.0)
    } else if y <= z {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}
