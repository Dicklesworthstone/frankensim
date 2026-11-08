//! Native physical-camera opacity, Planck/Kirchhoff source and spectral-unit tests.

use super::*;
use crate::camera::{Aperture, CameraProjection, FocusedDiffraction};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{
    AbsorbingGridBudget, CellPalette, CellPaletteBudget, MassAbsorptionSpectrum,
    OpticalMixtureBudget, ThermalEmissionComponent, ThermalEmissionModel,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey as ExecStreamKey};
use std::sync::Arc;

const SOURCE: ContentHash = ContentHash([113; 32]);

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                ExecStreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn camera(eye_z: f64, target_z: f64, half_tan: f64, radius: f64) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, eye_z),
        Point3::new(0.0, 0.0, target_z),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_half_tangent(half_tan).unwrap(),
        (eye_z - target_z).abs(),
        Aperture::try_circular(radius).unwrap(),
    )
    .unwrap()
}

fn settings(resolution: [usize; 2], samples_per_pixel: u32) -> VolumeCameraSettings {
    VolumeCameraSettings {
        resolution,
        samples_per_pixel,
        seed: 71,
        max_pixels: 64,
        max_samples: 10_000,
        max_cell_visits: 1_000_000,
    }
}

fn budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 16,
        max_steps_per_segment: 64,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-12 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

// Independent SI Planck equation, using moderate physical temperatures and
// wavelengths that do not need the production owner's range stabilization.
fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength = wavelength_nm * 1e-9;
    2.0 * h * c * c / wavelength.powi(5) / (h * c / (wavelength * k * temperature)).exp_m1() * 1e-9
}

#[test]
fn g1_spectral_camera_matches_planck_slab_and_isothermal_background_equilibrium() {
    with_cx(|cx, _| {
        let density = [2.0; 2];
        let temperature = [1200.0; 2];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 0.5]),
            MassAbsorptionSpectrum::source_table(
                &[
                    (800.0, 0.0),
                    (1000.0, 0.0),
                    (3000.0, 0.4),
                    (8000.0, 1.2),
                    (10_000.0, 1.2),
                ],
                SOURCE,
            )
            .unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let wavelengths = [1000.0, 3000.0, 8000.0];
        let view = camera(3.0, 0.0, 0.2, 0.0);
        let image =
            render_absorbing_camera(cx, &medium, &view, wavelengths, settings([3, 2], 1)).unwrap();
        assert_eq!(image.resolution(), [3, 2]);
        assert_eq!(image.wavelengths_nm(), &wavelengths);
        assert_eq!(image.field_identity(), medium.content_identity());
        let blackbody = wavelengths.map(|wavelength| planck(1200.0, wavelength));
        for py in 0..2 {
            for px in 0..3 {
                let x = (2.0 * (px as f64 + 0.5) / 3.0 - 1.0) * 0.3;
                let y = (1.0 - 2.0 * (py as f64 + 0.5) / 2.0) * 0.2;
                let length = (1.0 + x * x + y * y).sqrt();
                let pixel = &image.pixels()[px + 3 * py];
                let composite = pixel.composite(blackbody).unwrap();
                let ray = view
                    .generate_ray_from_tangent_offsets(
                        cx,
                        x,
                        y,
                        LensSample::try_new(0.5, 0.5).unwrap(),
                    )
                    .unwrap();
                let native = medium
                    .transport(cx, ray.origin, ray.dir, f64::INFINITY, &wavelengths)
                    .unwrap();
                for lane in 0..3 {
                    let tau = 2.0 * [0.0, 0.4, 1.2][lane] * length;
                    let opacity = -(-tau).exp_m1();
                    close(pixel.transmittance()[lane], (-tau).exp());
                    close(pixel.opacity()[lane], opacity);
                    close(pixel.effective_optical_depth()[lane], tau);
                    close(
                        pixel.emitted_radiance_w_m2_sr_nm()[lane],
                        blackbody[lane] * opacity,
                    );
                    close(
                        pixel.emitted_radiance_w_m2_sr_nm()[lane],
                        native.emission_w_m2_sr_nm[lane],
                    );
                    close(composite[lane], blackbody[lane]);
                }
            }
        }
        assert_eq!(
            image.approximate_emission_xyz(cx, [1.0; 3]),
            Err(VolumeCameraError::ObserverWavelength { lane: 0 }),
        );
    });
}

#[test]
fn g1_original_constituent_temperatures_reach_the_camera_without_averaging_temperature() {
    with_cx(|cx, _| {
        let source = ThermalEmissionModel::try_new(
            cx,
            &[
                ThermalEmissionComponent {
                    mass_weight: 0.25,
                    absorption: MassAbsorptionSpectrum::table(&[(380.0, 0.4), (780.0, 1.2)])
                        .unwrap(),
                    temperature_k: Some(900.0),
                },
                ThermalEmissionComponent {
                    mass_weight: 0.75,
                    absorption: MassAbsorptionSpectrum::table(&[(380.0, 2.0), (780.0, 0.4)])
                        .unwrap(),
                    temperature_k: Some(1800.0),
                },
            ],
            OpticalMixtureBudget { max_components: 2 },
        )
        .unwrap();
        let sources = CellPalette::try_new(
            cx,
            vec![source],
            Arc::from(&[0_u32][..]),
            CellPaletteBudget {
                max_models: 1,
                max_cells: 1,
            },
        )
        .unwrap();
        let density = [2.0];
        let medium = AbsorbingGrid::try_new_thermal_sources(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.3]),
            sources,
            budget(),
        )
        .unwrap();
        let wavelengths = [380.0, 523.0, 780.0];
        let image = render_absorbing_camera(
            cx,
            &medium,
            &camera(-1.0, 0.0, 0.1, 0.0),
            wavelengths,
            settings([1, 1], 1),
        )
        .unwrap();
        for (lane, wavelength) in wavelengths.into_iter().enumerate() {
            let first = 0.25 * (0.4 + 0.002 * (wavelength - 380.0));
            let second = 0.75 * (2.0 - 0.004 * (wavelength - 380.0));
            let absorption = first + second;
            let source = (first * planck(900.0, wavelength) + second * planck(1800.0, wavelength))
                / absorption;
            let tau = 0.6 * absorption;
            let pixel = &image.pixels()[0];
            close(
                pixel.emitted_radiance_w_m2_sr_nm()[lane],
                source * -(-tau).exp_m1(),
            );
            close(pixel.transmittance()[lane], (-tau).exp());
        }
    });
}

#[test]
fn g1_layer_order_changes_thermal_radiance_with_the_actual_foreground_opacity() {
    with_cx(|cx, _| {
        let density = [1.0; 2];
        let temperature = [900.0, 1800.0];
        let a = core::f64::consts::LN_2;
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::source_table(&[(1000.0, a), (5000.0, a)], SOURCE).unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let wavelengths = [1500.0, 4000.0];
        for (eye, target, near, far) in [(3.0, 1.0, 1800.0, 900.0), (-1.0, 1.0, 900.0, 1800.0)] {
            let image = render_absorbing_camera(
                cx,
                &medium,
                &camera(eye, target, 0.1, 0.0),
                wavelengths,
                settings([1, 1], 1),
            )
            .unwrap();
            for (lane, wavelength) in wavelengths.into_iter().enumerate() {
                close(image.pixels()[0].transmittance()[lane], 0.25);
                close(
                    image.pixels()[0].emitted_radiance_w_m2_sr_nm()[lane],
                    0.5 * planck(near, wavelength) + 0.25 * planck(far, wavelength),
                );
            }
        }
    });
}

#[test]
fn g0_camera_keeps_thin_and_opaque_public_limits_without_weakening_strict_columns() {
    with_cx(|cx, _| {
        let density = [2.0];
        let temperature = [1500.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-0.5, -0.5, 0.0], [1.0; 3]),
            MassAbsorptionSpectrum::source_table(
                &[
                    (800.0, 1e-20),
                    (1000.0, 1e-20),
                    (3000.0, f64::MAX),
                    (10_000.0, f64::MAX),
                ],
                SOURCE,
            )
            .unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let wavelengths = [1000.0, 3000.0];
        let image = render_absorbing_camera(
            cx,
            &medium,
            &camera(-1.0, 0.0, 1e-10, 0.0),
            wavelengths,
            settings([1, 1], 16),
        )
        .unwrap();
        let pixel = &image.pixels()[0];
        assert_eq!(pixel.transmittance()[0], 1.0);
        close(pixel.opacity()[0], 2e-20);
        close(
            pixel.emitted_radiance_w_m2_sr_nm()[0],
            2e-20 * planck(1500.0, 1000.0),
        );
        assert_eq!(pixel.transmittance()[1], 0.0);
        assert_eq!(pixel.effective_optical_depth()[1], f64::INFINITY);
        close(
            pixel.emitted_radiance_w_m2_sr_nm()[1],
            planck(1500.0, 3000.0),
        );
        assert!(matches!(
            medium.transport_observed::<2, true>(
                cx,
                Point3::new(0.0, 0.0, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
                &wavelengths,
                1,
            ),
            Err(VolumeTransportError::Unrepresentable),
        ));

        // A finite but opaque foreground makes a farther positive source
        // unrepresentably faint. That source is not required for this image.
        let density = [1.0; 2];
        let temperature = [1500.0, 1600.0];
        let layered = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-0.5, -0.5, 0.0], [1.0; 3]),
            MassAbsorptionSpectrum::source_table(&[(1000.0, 800.0), (5000.0, 800.0)], SOURCE)
                .unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let image = render_absorbing_camera(
            cx,
            &layered,
            &camera(-1.0, 0.0, 0.1, 0.0),
            [3000.0],
            settings([1, 1], 1),
        )
        .unwrap();
        close(
            image.pixels()[0].emitted_radiance_w_m2_sr_nm()[0],
            planck(1500.0, 3000.0),
        );
        close(image.pixels()[0].effective_optical_depth()[0], 1600.0);
        assert!(matches!(
            layered.transport_observed::<1, true>(
                cx,
                Point3::new(0.0, 0.0, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
                &[3000.0],
                2,
            ),
            Err(VolumeTransportError::Unrepresentable),
        ));
    });
}

#[test]
fn g1_xyz_uses_explicit_nm_weights_for_spectral_density() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [2000.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-0.5, -0.5, 0.0], [1.0; 3]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let wavelengths = [500.0, 650.0];
        let image = render_absorbing_camera(
            cx,
            &medium,
            &camera(-1.0, 0.0, 0.1, 0.0),
            wavelengths,
            settings([1, 1], 1),
        )
        .unwrap();
        let weights = [20.0, 80.0];
        let xyz = image.approximate_emission_xyz(cx, weights).unwrap()[0];
        let doubled = image.approximate_emission_xyz(cx, [40.0, 160.0]).unwrap()[0];
        let first = [cie_x(500.0), cie_y(500.0), cie_z(500.0)];
        let second = [cie_x(650.0), cie_y(650.0), cie_z(650.0)];
        for channel in 0..3 {
            let expected = -(-1.0_f64).exp_m1()
                * (planck(2000.0, 500.0) * weights[0] * first[channel]
                    + planck(2000.0, 650.0) * weights[1] * second[channel])
                / y_integral();
            close(xyz[channel], expected);
            close(doubled[channel], 2.0 * xyz[channel]);
        }
        for bad in [[-1.0, 80.0], [f64::NAN, 80.0], [f64::INFINITY, 80.0]] {
            assert_eq!(
                image.approximate_emission_xyz(cx, bad),
                Err(VolumeCameraError::Transport(VolumeTransportError::Sample)),
            );
        }
    });
}

#[test]
fn g0_absorption_only_camera_charges_cells_once_per_packet_and_preserves_refusals() {
    with_cx(|cx, gate| {
        let density = [1.0; 2];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 0.5]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let view = camera(3.0, 0.0, 0.1, 0.0);
        let base = settings([2, 1], 1);
        let wavelengths = [410.0, 550.0, 700.0];
        let complete = render_absorbing_camera(
            cx,
            &medium,
            &view,
            wavelengths,
            VolumeCameraSettings {
                max_cell_visits: 4,
                ..base
            },
        )
        .unwrap();
        for pixel in complete.pixels() {
            assert_eq!(pixel.emitted_radiance_w_m2_sr_nm(), &[0.0; 3]);
            for value in pixel.composite([1.0; 3]).unwrap() {
                assert!(value > 0.0 && value < 1.0);
            }
        }
        assert_eq!(
            render_absorbing_camera(
                cx,
                &medium,
                &view,
                wavelengths,
                VolumeCameraSettings {
                    max_cell_visits: 3,
                    ..base
                },
            ),
            Err(VolumeCameraError::Transport(VolumeTransportError::Budget)),
        );
        assert_eq!(
            render_absorbing_camera(cx, &medium, &view, [350.0], base),
            Err(VolumeCameraError::Transport(
                VolumeTransportError::Wavelength
            )),
        );
        assert_eq!(
            render_absorbing_camera(cx, &medium, &view, [], base),
            Err(VolumeCameraError::Transport(VolumeTransportError::Sample)),
        );
        let diffraction = camera(-1.0, 0.0, 0.1, 0.1)
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
        assert_eq!(
            render_absorbing_camera(cx, &medium, &diffraction, wavelengths, base),
            Err(VolumeCameraError::Camera(
                CameraError::DiffractionRequiresSpectralPath
            )),
        );
        let zero_density = [0.0; 2];
        let vacuum = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &zero_density, [-2.0, -2.0, 0.0], [4.0, 4.0, 0.5]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let clear = render_absorbing_camera(
            cx,
            &vacuum,
            &view,
            wavelengths,
            VolumeCameraSettings {
                max_cell_visits: 0,
                ..base
            },
        )
        .unwrap();
        for pixel in clear.pixels() {
            assert_eq!(pixel.transmittance(), &[1.0; 3]);
            assert_eq!(pixel.opacity(), &[0.0; 3]);
        }
        gate.request();
        assert_eq!(
            render_absorbing_camera(cx, &medium, &view, wavelengths, base),
            Err(VolumeCameraError::Transport(
                VolumeTransportError::Cancelled
            )),
        );
        assert_eq!(
            complete.approximate_emission_xyz(cx, [1.0; 3]),
            Err(VolumeCameraError::Transport(
                VolumeTransportError::Cancelled
            )),
        );
    });
}
