//! Compare coherent, partially coherent and incoherent pupil illumination.
//! Run: cargo run -p fs-render --example coherent_pupil
//! Output is a CSV sensor cut, in watts per tangent-coordinate area.

use fs_geom::{Point3, Vec3};
use fs_math::c64::C64;
use fs_render::camera::{
    Aperture, CameraProjection, CoherentDiffractionLimits, CoherentPupilIllumination,
    CoherentPupilWave, FocusedDiffraction, PhysicalCamera, PupilCoherenceFactor,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let wavelength_nm = 550.0;
    let radius_m = 0.001;
    let scale = wavelength_nm * 1e-9 / (core::f64::consts::TAU * radius_m);
    let camera = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_vertical_fov(0.5)?,
        2.0,
        Aperture::try_circular(radius_m)?,
    )?
    .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0)?)?;
    let waves = [-1.5, 1.5].map(|q| CoherentPupilWave {
        tangent_direction: [q * scale, 0.0],
        amplitude_sqrt_w: C64::from_re(core::f64::consts::FRAC_1_SQRT_2),
    });
    let independent = [
        PupilCoherenceFactor { waves: &waves[..1] },
        PupilCoherenceFactor { waves: &waves[1..] },
    ];
    // Covariance diagonals stay 0.5 W; the correlation coefficient is 0.5.
    // Factor weights are AMPLITUDES, so the square roots are essential.
    let plus = waves.map(|mut wave| {
        wave.amplitude_sqrt_w = wave.amplitude_sqrt_w.scale(0.75_f64.sqrt());
        wave
    });
    let minus = [
        CoherentPupilWave {
            amplitude_sqrt_w: waves[0].amplitude_sqrt_w.scale(0.5),
            ..waves[0]
        },
        CoherentPupilWave {
            amplitude_sqrt_w: waves[1].amplitude_sqrt_w.scale(-0.5),
            ..waves[1]
        },
    ];
    let partial = [
        PupilCoherenceFactor { waves: &plus },
        PupilCoherenceFactor { waves: &minus },
    ];
    let limits = CoherentDiffractionLimits {
        maximum_waves: 4,
        maximum_factors: 2,
    };
    let coherent = camera.prepare_coherent_diffraction(
        wavelength_nm,
        CoherentPupilIllumination {
            mean: &waves,
            factors: &[],
        },
        limits,
        || false,
    )?;
    let incoherent = camera.prepare_coherent_diffraction(
        wavelength_nm,
        CoherentPupilIllumination {
            mean: &[],
            factors: &independent,
        },
        limits,
        || false,
    )?;
    let partially_coherent = camera.prepare_coherent_diffraction(
        wavelength_nm,
        CoherentPupilIllumination {
            mean: &[],
            factors: &partial,
        },
        limits,
        || false,
    )?;
    println!("tangent_x,coherent_w,partial_w,incoherent_w");
    for i in 0..161 {
        let tangent_x = (f64::from(i) - 80.0) * 0.1 * scale;
        let at = [tangent_x, 0.0];
        println!(
            "{tangent_x:.12e},{:.12e},{:.12e},{:.12e}",
            coherent.observe(at, || false)?.power_density_w,
            partially_coherent.observe(at, || false)?.power_density_w,
            incoherent.observe(at, || false)?.power_density_w
        );
    }
    Ok(())
}
