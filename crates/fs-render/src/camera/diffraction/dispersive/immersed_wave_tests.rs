use super::*;
use crate::camera::{Aperture, CameraProjection, DefocusedDiffraction, LayeredDiffraction};
use crate::dielectric::{DielectricIorSample, SellmeierIor, TabulatedIor};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

const PI: f64 = core::f64::consts::PI;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn constant(n: f64) -> DielectricIor {
    SellmeierIor::try_new_in_band([n * n - 1.0, 0.0, 0.0], [0.0; 3], (200.0, 1400.0))
        .unwrap()
        .into()
}

fn exterior() -> DielectricIor {
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            let wavelength = 350.0 + i as f64 * 50.0;
            DielectricIorSample::try_new_in_band(wavelength, 1.25 + (550.0 - wavelength) * 1e-4)
                .unwrap()
        }),
        50.0,
    )
    .unwrap()
    .into()
}

fn lens(exterior: DielectricIor, glass: f64) -> DispersiveThinLens {
    DispersiveThinLens::try_new_in_medium(
        550.0,
        &[ThinLensElement::try_new(constant(glass), 20.0, -20.0).unwrap()],
        exterior,
    )
    .unwrap()
}

fn camera(lens: &DispersiveThinLens, aperture: Aperture) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.002).unwrap(),
        2.0,
        aperture,
    )
    .unwrap()
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn pupils(radius: f64) -> [Aperture; 3] {
    [
        Aperture::try_circular(radius).unwrap(),
        Aperture::try_annular(radius, 0.6 * radius).unwrap(),
        Aperture::try_regular_polygon(radius, 5, 0.27).unwrap(),
    ]
}

// Direct two-dimensional complex pupil integration, independent of the native
// disk Bessel recurrence, annular disk subtraction and polygon edge integrals.
// The physical Fourier density is |F|^2 / ((2*pi)^2 * area * scale^2).
fn direct_pupil_density(aperture: &Aperture, q: [f64; 2], phase: f64, scale: f64) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    let mut add = |x: f64, y: f64, weight: f64| {
        let angle = phase * (x * x + y * y) - q[0] * x - q[1] * y;
        re += weight * angle.cos();
        im += weight * angle.sin();
    };
    let area = if let Some(blades) = aperture.blades() {
        let step = 2.0 * PI / f64::from(blades);
        let rotation = aperture.rotation_rad().unwrap();
        const N: usize = 160;
        for edge in 0..blades {
            let a = rotation + f64::from(edge) * step;
            let b = a + step;
            for radial in 0..N {
                let r = (radial as f64 + 0.5) / N as f64;
                for lateral in 0..N {
                    let t = (lateral as f64 + 0.5) / N as f64;
                    add(
                        r * ((1.0 - t) * a.cos() + t * b.cos()),
                        r * ((1.0 - t) * a.sin() + t * b.sin()),
                        r * step.sin() / (N * N) as f64,
                    );
                }
            }
        }
        0.5 * f64::from(blades) * step.sin()
    } else {
        let inner = aperture.inner_radius_m() / aperture.radius_m();
        let clear = 1.0 - inner * inner;
        const NR: usize = 256;
        const NA: usize = 256;
        for radial in 0..NR {
            let r = (inner * inner + clear * (radial as f64 + 0.5) / NR as f64).sqrt();
            for angular in 0..NA {
                let angle = 2.0 * PI * (angular as f64 + 0.5) / NA as f64;
                add(
                    r * angle.cos(),
                    r * angle.sin(),
                    PI * clear / (NR * NA) as f64,
                );
            }
        }
        PI * clear
    };
    (re * re + im * im) / (4.0 * PI * PI * area * scale * scale)
}

#[test]
fn immersed_fresnel_psfs_match_independent_complex_opening_integrals() {
    let lens = lens(exterior(), 2.0);
    let depth = 1.7;
    for aperture in pupils(0.0008) {
        let camera = camera(&lens, aperture.clone())
            .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(depth), 0.01).unwrap())
            .unwrap()
            .with_dispersive_lens(lens.clone())
            .unwrap();
        for wavelength in [400.0, 550.0, 700.0] {
            let ambient = 1.25 + (550.0 - wavelength) * 1e-4;
            let lambda = wavelength * 1e-9 / ambient;
            let radius = aperture.radius_m();
            let scale = lambda / (2.0 * PI * radius);
            let power = (2.0 / ambient - 1.0) * 40.0;
            let phase = PI * radius * radius / lambda * (1.0 / depth - 0.5 + 24.0 - power);
            for q in [[0.0, 0.0], [1.2, -0.7], [3.1, 2.4]] {
                let expected = direct_pupil_density(&aperture, q, phase, scale);
                let actual = camera
                    .defocused_diffraction_density(q.map(|v| v * scale), wavelength)
                    .unwrap();
                near(actual, expected, 1.5e-4 * expected.max(1.0));
            }
        }
    }
}

#[test]
fn common_index_scaling_changes_wave_resolution_and_retains_integrated_radiance() {
    with_cx(|cx| {
        let small = lens(constant(1.25), 2.0);
        let large = lens(constant(2.5), 4.0);
        assert_eq!(small.reference_power_per_m(), large.reference_power_per_m());
        for aperture in pupils(0.0008) {
            let cameras = [&small, &large].map(|lens| {
                camera(lens, aperture.clone())
                    .with_defocused_diffraction(
                        DefocusedDiffraction::try_new(Some(2.0), 0.0).unwrap(),
                    )
                    .unwrap()
                    .with_dispersive_lens(lens.clone())
                    .unwrap()
            });
            for angle in [[0.0, 0.0], [2e-5, -3e-5], [1e-4, 8e-5]] {
                let expected = 4.0
                    * cameras[0]
                        .defocused_diffraction_density(angle.map(|x| 2.0 * x), 600.0)
                        .unwrap();
                near(
                    cameras[1]
                        .defocused_diffraction_density(angle, 600.0)
                        .unwrap(),
                    expected,
                    1e-12 * expected.max(1.0),
                );
            }
        }
        // The actual camera sampler integrates every original vacuum lane.
        // A missing square Jacobian or use of vacuum-width ratios loses power.
        let dispersive = lens(exterior(), 2.0);
        let camera = camera(&dispersive, Aperture::try_circular(0.0008).unwrap())
            .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(2.0), 0.0).unwrap())
            .unwrap()
            .with_dispersive_lens(dispersive)
            .unwrap();
        let mut energy = [0.0; 3];
        const N: u32 = 8192;
        for sample in 0..N {
            let result = camera
                .generate_spectral_ray_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    LensSample::try_new((f64::from(sample) + 0.5) / f64::from(N), 0.37).unwrap(),
                    &[400.0, 550.0, 700.0],
                )
                .unwrap();
            for (sum, weight) in energy.iter_mut().zip(result.weights) {
                *sum += weight / f64::from(N);
            }
        }
        for sum in energy {
            near(sum, 1.0, 0.008);
        }
    });
}

#[test]
fn layered_immersion_uses_actual_source_depth_and_physical_lane_jacobians() {
    with_cx(|cx| {
        let lens = lens(exterior(), 2.0);
        let camera = camera(&lens, Aperture::try_annular(0.0008, 0.00048).unwrap())
            .with_layered_diffraction(LayeredDiffraction::try_new(&[1.7, 2.4], 0.001).unwrap())
            .unwrap()
            .with_dispersive_lens(lens)
            .unwrap();
        let waves = [400.0, 550.0, 700.0];
        let sample = LensSample::try_new(0.32, 0.27).unwrap();
        let mut resolved = [[0.0; 3]; 2];
        for (slot, depth) in resolved.iter_mut().zip([1.7, 2.4]) {
            let token = camera
                .sample_layered_diffraction_ray(cx, 0.0, 0.0, sample, &waves)
                .unwrap();
            let ray = token.ray();
            let angle = [
                ray.dir.dot(camera.right()) / ray.dir.dot(camera.forward()),
                ray.dir.dot(camera.up()) / ray.dir.dot(camera.forward()),
            ];
            let distance = depth / ray.dir.dot(camera.forward());
            let point = Point3::new(
                ray.dir.x * distance,
                ray.dir.y * distance,
                ray.dir.z * distance,
            );
            let result = token.resolve_source(cx, Some(point)).unwrap();
            let proposal = camera
                .layered_diffraction_density(angle, depth, waves[0])
                .unwrap()
                / result.weights[0];
            for (lane, &wave) in waves.iter().enumerate() {
                let density = camera
                    .layered_diffraction_density(angle, depth, wave)
                    .unwrap();
                near(
                    result.weights[lane] * proposal,
                    density,
                    2e-12 * density.max(1.0),
                );
            }
            *slot = result.weights;
        }
        assert_ne!(resolved[0], resolved[1]);
        let black = camera
            .sample_layered_diffraction_ray(cx, 0.0, 0.0, sample, &waves)
            .unwrap()
            .resolve_source(cx, None)
            .unwrap();
        assert_eq!(black.weights, [0.0; 3]);
    });
}

#[test]
fn depth_error_uses_shortest_medium_wavelength_not_shortest_vacuum_lane() {
    let exterior = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            let wave = 400.0 + i as f64 * 25.0;
            DielectricIorSample::try_new_in_band(wave, 1.2 + (wave - 400.0) * 0.009).unwrap()
        }),
        25.0,
    )
    .unwrap();
    let lens = lens(exterior.into(), 5.0);
    let radius = 0.0001;
    let plane_error = 0.08 / (PI * radius * radius / (400e-9 / 1.2));
    let depth = 1.0 / (0.5 + plane_error);
    let base = camera(&lens, Aperture::try_circular(radius).unwrap());
    let single = base
        .clone()
        .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(2.0), 0.1).unwrap())
        .unwrap()
        .with_dispersive_lens(lens.clone())
        .unwrap();
    let layered = base
        .with_layered_diffraction(LayeredDiffraction::try_new(&[2.0], 0.1).unwrap())
        .unwrap()
        .with_dispersive_lens(lens)
        .unwrap();
    let point = Some(Point3::new(0.0, 0.0, depth));
    for camera in [single, layered] {
        camera.admit_diffraction_depth(point, &[400.0]).unwrap();
        assert_eq!(
            camera.admit_diffraction_depth(point, &[400.0, 600.0]),
            Err(CameraError::DefocusOutsideDiffractionModel)
        );
    }
}

#[test]
fn proposal_coordinates_never_replace_original_vacuum_source_coordinates() {
    with_cx(|cx| {
        let waves = core::array::from_fn(|i| 400.0 + 37.5 * i as f64);
        let exterior = TabulatedIor::try_new_in_band(
            waves.map(|wave| DielectricIorSample::try_new_in_band(wave, 4.0).unwrap()),
            37.5,
        )
        .unwrap();
        // Every physical wavelength and their proposal midpoint lie outside
        // this source's VACUUM band. They are propagation coordinates only.
        assert!(exterior.eval(137.5).is_err());
        let lens = DispersiveThinLens::try_new_in_medium(
            402.0,
            &[ThinLensElement::try_new(constant(5.0), 20.0, -20.0).unwrap()],
            exterior.into(),
        )
        .unwrap();
        let base = camera(&lens, Aperture::try_circular(0.0008).unwrap());
        let single = base
            .clone()
            .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(2.0), 0.0).unwrap())
            .unwrap()
            .with_dispersive_lens(lens.clone())
            .unwrap();
        let layered = base
            .with_layered_diffraction(LayeredDiffraction::try_new(&[2.0], 0.0).unwrap())
            .unwrap()
            .with_dispersive_lens(lens)
            .unwrap();
        single
            .generate_spectral_ray_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                LensSample::CENTER,
                &[402.0, 698.0],
            )
            .unwrap();
        let pending = layered
            .sample_layered_diffraction_ray(cx, 0.0, 0.0, LensSample::CENTER, &[402.0, 698.0])
            .unwrap();
        let ray = pending.ray();
        let distance = 2.0 / ray.dir.z;
        pending
            .resolve_source(
                cx,
                Some(Point3::new(ray.dir.x * distance, ray.dir.y * distance, 2.0)),
            )
            .unwrap();
        assert!(
            single
                .defocused_diffraction_density([0.0; 2], 399.0)
                .is_err()
        );
        assert!(
            single
                .defocused_diffraction_density([0.0; 2], 701.0)
                .is_err()
        );
    });
}

#[cfg(feature = "material-state")]
#[test]
fn sourced_water_state_drives_immersed_wave_resolution_at_original_vacuum_wavelengths() {
    use crate::dielectric::IapwsWaterIor;
    with_cx(|cx| {
        let states = [(293.15, 998.2), (333.15, 983.2)];
        let mut peaks = [0.0; 2];
        for (slot, (temperature, density)) in peaks.iter_mut().zip(states) {
            let water = IapwsWaterIor::try_new(temperature, density).unwrap();
            let lens = lens(water.into(), 2.0);
            let aperture = Aperture::try_circular(0.0008).unwrap();
            let area = aperture.transmitting_area_m2().unwrap();
            let camera = camera(&lens, aperture)
                .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(2.0), 0.0).unwrap())
                .unwrap()
                .with_dispersive_lens(lens)
                .unwrap();
            let n = water.eval(550.0).unwrap();
            let expected = area * (n / 550e-9).powi(2);
            *slot = camera
                .defocused_diffraction_density([0.0; 2], 550.0)
                .unwrap();
            near(*slot, expected, expected * 1e-13);
            let draw = camera
                .generate_spectral_ray_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    LensSample::try_new(0.32, 0.61).unwrap(),
                    &[400.0, 550.0, 900.0],
                )
                .unwrap();
            assert!(draw.weights.into_iter().all(|x| x.is_finite() && x > 0.0));
            assert!(
                camera
                    .defocused_diffraction_density([0.0; 2], 1101.0)
                    .is_err()
            );
        }
        assert!(peaks[0] > peaks[1]);
    });
}
