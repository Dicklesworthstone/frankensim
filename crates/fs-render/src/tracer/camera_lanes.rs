//! Native per-wavelength camera geometry shared by cinematic transport.
//! Lens prescriptions own all optical admission and throughput; this adapter
//! only presents their accepted rays to the existing spectral path state.
use super::*;

#[derive(Clone, Copy)]
pub(super) struct SpectralCameraLaunch {
    pub(super) ray: Ray,
    pub(super) weight: f64,
    pub(super) scene_phase_index: f64,
}

pub(super) fn has_spectral_camera_geometry(camera: &PhysicalCamera) -> bool {
    camera.sequential_lens().is_some()
        || (camera.dispersive_lens().is_some() && !camera.has_diffraction())
}

/// Preserve a real lens's clipped samples and throughput, or the original
/// ideal quadratic phase lens's unit aperture-average radiance weight.
pub(super) fn spectral_camera_launches(
    cx: &Cx<'_>,
    camera: &PhysicalCamera,
    tangent_xy: [f64; 2],
    pupil: LensSample,
    wavelengths_nm: &[f64; PACKET],
) -> Result<[Option<SpectralCameraLaunch>; PACKET], CameraError> {
    if camera.sequential_lens().is_some() {
        return Ok(camera
            .generate_sequential_spectral_rays_from_tangent_offsets(
                cx,
                tangent_xy[0],
                tangent_xy[1],
                pupil,
                wavelengths_nm,
            )?
            .map(|sample| {
                sample.map(|sample| SpectralCameraLaunch {
                    ray: sample.ray,
                    weight: sample.weight,
                    scene_phase_index: sample.scene_phase_index,
                })
            }));
    }
    let lens = camera
        .dispersive_lens()
        .ok_or(CameraError::InvalidDispersiveLens)?;
    let rays = camera.generate_geometric_spectral_rays_from_tangent_offsets(
        cx,
        tangent_xy[0],
        tangent_xy[1],
        pupil,
        wavelengths_nm,
    )?;
    let mut launches = [None; PACKET];
    for (lane, ray) in rays.into_iter().enumerate() {
        launches[lane] = Some(SpectralCameraLaunch {
            ray,
            weight: 1.0,
            scene_phase_index: lens.exterior_phase_index(wavelengths_nm[lane])?,
        });
    }
    Ok(launches)
}
