//! Physical slab, spectral shielding, pupil coverage and bounded camera tests.

use super::*;
use crate::camera::{Aperture, CameraError, CameraProjection, FocusedDiffraction, LensSample};
use crate::spectral::{cie_x, cie_y, cie_z, y_integral};
use crate::volumes::{DvrError, VolumeGrid};
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, VolumeTransportError,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey as ExecStreamKey};
use fs_geom::{Point3, Vec3};

const SOURCE: ContentHash = ContentHash([109; 32]);
const FOUR_PI: f64 = 4.0 * core::f64::consts::PI;

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

fn camera(eye_z: f64, target_z: f64, half_tan: f64, focus: f64, radius: f64) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, eye_z),
        Point3::new(0.0, 0.0, target_z),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_half_tangent(half_tan).unwrap(),
        focus,
        Aperture::try_circular(radius).unwrap(),
    )
    .unwrap()
}

fn grid<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    dims: [usize; 3],
    origin: [f64; 3],
    cell: [f64; 3],
    absorption: MassAbsorptionSpectrum,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(dims, density, origin, cell),
        absorption,
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 64,
        },
    )
    .unwrap()
}

fn settings(resolution: [usize; 2], samples_per_pixel: u32) -> LuminescenceCameraSettings {
    LuminescenceCameraSettings {
        resolution,
        samples_per_pixel,
        seed: 71,
        max_pixels: 64,
        max_samples: 100_000,
        max_cell_visits: 1_000_000,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_camera_matches_spectral_slab_and_original_ray_integral_at_every_pixel() {
    with_cx(|cx, _| {
        let density = [1.0; 2];
        let medium = grid(
            cx,
            &density,
            [1, 1, 2],
            [-2.0, -2.0, 0.0],
            [4.0, 4.0, 0.5],
            MassAbsorptionSpectrum::table(&[
                (380.0, 0.25),
                (500.0, 0.25),
                (650.0, 1.5),
                (780.0, 1.5),
            ])
            .unwrap(),
        );
        // Each cell has volume 8 m³ and isotropic emissivity [2, 3].
        let powers = [(0, [16.0 * FOUR_PI, 24.0 * FOUR_PI]), (1, [16.0 * FOUR_PI, 24.0 * FOUR_PI])];
        let field = LuminescenceGrid::try_new(cx, &medium, SOURCE, [500.0, 650.0], &powers, 2)
            .unwrap();
        let view = camera(3.0, 0.0, 0.2, 3.0, 0.0);
        let image = render_luminescence_camera(cx, &field, &view, settings([3, 2], 1)).unwrap();
        assert_eq!(image.resolution(), [3, 2]);
        assert_eq!(image.wavelengths_nm(), &[500.0, 650.0]);
        assert_eq!(image.source_identity(), SOURCE);
        assert_eq!(image.field_identity(), field.content_identity());
        for py in 0..2 {
            for px in 0..3 {
                let x = (2.0 * (px as f64 + 0.5) / 3.0 - 1.0) * 0.3;
                let y = (1.0 - 2.0 * (py as f64 + 0.5) / 2.0) * 0.2;
                let length = (1.0 + x * x + y * y).sqrt();
                let pixel = &image.pixels()[px + 3 * py];
                for lane in 0..2 {
                    let alpha = [0.25, 1.5][lane];
                    let emissivity = [2.0, 3.0][lane];
                    let tau = alpha * length;
                    close(pixel.transmittance()[lane], (-tau).exp());
                    close(pixel.opacity()[lane], -(-tau).exp_m1());
                    close(pixel.effective_optical_depth()[lane], tau);
                    close(
                        pixel.emitted_radiance_w_m2_sr()[lane],
                        emissivity * -(-tau).exp_m1() / alpha,
                    );
                }
                let ray = view
                    .generate_ray_from_tangent_offsets(cx, x, y, LensSample::try_new(0.5, 0.5).unwrap())
                    .unwrap();
                let original = field.transport(cx, ray.origin, ray.dir, f64::INFINITY).unwrap();
                for lane in 0..2 {
                    close(pixel.transmittance()[lane], original.transmittance[lane]);
                    close(
                        pixel.emitted_radiance_w_m2_sr()[lane],
                        original.emitted_radiance_w_m2_sr[lane],
                    );
                }
            }
        }
    });
}

#[test]
fn g1_foreground_order_and_xyz_keep_each_original_integrated_spectral_power() {
    with_cx(|cx, _| {
        let density = [1.0; 2];
        let a = core::f64::consts::LN_2;
        let medium = grid(
            cx,
            &density,
            [1, 1, 2],
            [-1.0, -1.0, 0.0],
            [2.0, 2.0, 1.0],
            MassAbsorptionSpectrum::table(&[
                (380.0, a), (510.0, a), (650.0, 2.0 * a), (780.0, 2.0 * a),
            ])
            .unwrap(),
        );
        let field = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [510.0, 650.0],
            &[(1, [4.0 * FOUR_PI, 8.0 * FOUR_PI])], 1,
        ).unwrap();
        let near = render_luminescence_camera(
            cx, &field, &camera(-1.0, 1.0, 0.1, 2.0, 0.0), settings([1, 1], 1),
        ).unwrap();
        let far = render_luminescence_camera(
            cx, &field, &camera(3.0, 1.0, 0.1, 2.0, 0.0), settings([1, 1], 1),
        ).unwrap();
        let front = near.pixels()[0];
        let back = far.pixels()[0];
        let expected = [0.25 / a, 0.375 / (2.0 * a)];
        for lane in 0..2 {
            close(front.emitted_radiance_w_m2_sr()[lane], expected[lane]);
            close(front.transmittance()[lane], [0.25, 0.0625][lane]);
            close(back.transmittance()[lane], front.transmittance()[lane]);
            close(
                back.emitted_radiance_w_m2_sr()[lane],
                expected[lane] / [0.5, 0.25][lane],
            );
        }
        let composed = front.composite([3.0, 5.0]).unwrap();
        close(composed[0], expected[0] + 0.75);
        close(composed[1], expected[1] + 5.0 / 16.0);
        let xyz = near.approximate_emission_xyz(cx).unwrap()[0];
        let first = [cie_x(510.0), cie_y(510.0), cie_z(510.0)];
        let second = [cie_x(650.0), cie_y(650.0), cie_z(650.0)];
        for channel in 0..3 {
            // No inferred 140 nm interval, integration width, or photon factor.
            close(xyz[channel], (expected[0] * first[channel] + expected[1] * second[channel]) / y_integral());
        }
        assert_eq!(
            front.composite([-1.0, 0.0]),
            Err(LuminescenceCameraError::Image(DvrError::InvalidBackgroundRadiance)),
        );
    });
}

#[test]
fn g1_zero_density_emitters_remain_visible_and_rows_follow_camera_up() {
    with_cx(|cx, _| {
        let density = [0.0; 2];
        let medium = grid(
            cx, &density, [1, 2, 1], [-1.0, -1.0, 0.0], [2.0, 1.0, 1.0],
            MassAbsorptionSpectrum::constant(100.0).unwrap(),
        );
        assert!(medium.is_vacuum());
        let field = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [600.0],
            &[(0, [2.0 * FOUR_PI]), (1, [6.0 * FOUR_PI])], 2,
        ).unwrap();
        let image = render_luminescence_camera(
            cx, &field, &camera(3.0, 0.0, 0.2, 3.0, 0.0), settings([1, 2], 1),
        ).unwrap();
        let length = 1.01_f64.sqrt();
        close(image.pixels()[0].emitted_radiance_w_m2_sr()[0], 3.0 * length);
        close(image.pixels()[1].emitted_radiance_w_m2_sr()[0], length);
        for pixel in image.pixels() {
            assert_eq!(pixel.transmittance(), &[1.0]);
            assert_eq!(pixel.opacity(), &[0.0]);
        }
    });
}

#[test]
fn g1_defocused_spectral_strip_matches_circular_pupil_coverage_and_refocuses() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = grid(
            cx, &density, [1; 3], [-0.125, -2.0, 2.0], [0.25, 4.0, 1e-8],
            MassAbsorptionSpectrum::table(&[
                (380.0, 1e12), (500.0, 1e12), (650.0, 0.0), (780.0, 0.0),
            ]).unwrap(),
        );
        // Cell volume is 1e-8 m³; first-lane j/alpha is 2 W/(m² sr).
        let field = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [500.0, 650.0], &[(0, [2e4 * FOUR_PI, 0.0])], 1,
        ).unwrap();
        let defocused = camera(0.0, 2.0, 1e-10, 4.0, 1.0);
        let image = render_luminescence_camera(cx, &field, &defocused, settings([1, 1], 16_384))
            .unwrap();
        let pixel = image.pixels()[0];
        let q = 0.25_f64;
        let coverage = 2.0 / core::f64::consts::PI * (q.asin() + q * (1.0 - q * q).sqrt());
        assert!((pixel.opacity()[0] - coverage).abs() < 0.02);
        close(pixel.transmittance()[0], 1.0 - pixel.opacity()[0]);
        close(pixel.emitted_radiance_w_m2_sr()[0], 2.0 * pixel.opacity()[0]);
        close(pixel.effective_optical_depth()[0], -pixel.transmittance()[0].ln());
        assert_eq!(pixel.transmittance()[1], 1.0);
        assert_eq!(pixel.opacity()[1], 0.0);
        let focused = render_luminescence_camera(
            cx, &field, &camera(0.0, 2.0, 1e-10, 2.0, 1.0), settings([1, 1], 64),
        ).unwrap();
        assert_eq!(focused.pixels()[0].transmittance()[0], 0.0);
        close(focused.pixels()[0].emitted_radiance_w_m2_sr()[0], 2.0);
        let repeat = settings([1, 1], 128);
        assert_eq!(
            render_luminescence_camera(cx, &field, &defocused, repeat).unwrap(),
            render_luminescence_camera(cx, &field, &defocused, repeat).unwrap(),
        );
    });
}

#[test]
fn g0_thin_opacity_and_bright_backgrounds_survive_survival_rounding() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = grid(
            cx, &density, [1; 3], [-0.5, -0.5, 0.0], [1.0; 3],
            MassAbsorptionSpectrum::table(&[
                (380.0, 1e-20), (500.0, 1e-20), (650.0, 800.0), (780.0, 800.0),
            ]).unwrap(),
        );
        let field = LuminescenceGrid::try_new(cx, &medium, SOURCE, [500.0, 650.0], &[], 0)
            .unwrap();
        let image = render_luminescence_camera(
            cx, &field, &camera(-1.0, 0.0, 1e-8, 2.0, 0.0), settings([1, 1], 32),
        ).unwrap();
        let pixel = image.pixels()[0];
        assert_eq!(pixel.transmittance()[0], 1.0);
        close(pixel.opacity()[0], 1e-20);
        close(pixel.effective_optical_depth()[0], 1e-20);
        assert_eq!(pixel.transmittance()[1], 0.0);
        close(pixel.effective_optical_depth()[1], 800.0);
        let composed = pixel.composite([0.0, 1e300]).unwrap();
        close(composed[1], (1e300_f64.ln() - 800.0).exp());
        assert!(composed[1] > 0.0);
    });
}

#[test]
fn g0_finite_bright_camera_means_do_not_overflow_a_sample_sum() {
    with_cx(|cx, _| {
        let density = [0.0];
        let medium = grid(
            cx, &density, [1; 3], [-0.5, -0.5, 0.0], [1.0; 3],
            MassAbsorptionSpectrum::constant(0.0).unwrap(),
        );
        let field = LuminescenceGrid::try_new(cx, &medium, SOURCE, [550.0], &[(0, [1e308])], 1)
            .unwrap();
        let image = render_luminescence_camera(
            cx, &field, &camera(-1.0, 0.0, 1e-10, 2.0, 0.0), settings([1, 1], 64),
        ).unwrap();
        close(image.pixels()[0].emitted_radiance_w_m2_sr()[0], 1e308 / FOUR_PI);
    });
}

#[test]
fn g0_total_cell_limits_misses_camera_refusals_and_cancellation_are_preserved() {
    with_cx(|cx, gate| {
        let density = [0.0; 2];
        let medium = grid(
            cx, &density, [1, 1, 2], [-2.0, -2.0, 0.0], [4.0, 4.0, 0.5],
            MassAbsorptionSpectrum::constant(0.0).unwrap(),
        );
        let field = LuminescenceGrid::try_new(cx, &medium, SOURCE, [550.0], &[(1, [1.0])], 1)
            .unwrap();
        let view = camera(3.0, 0.0, 0.1, 3.0, 0.0);
        let base = settings([2, 1], 1);
        let complete = render_luminescence_camera(
            cx, &field, &view, LuminescenceCameraSettings { max_cell_visits: 4, ..base },
        ).unwrap();
        assert_eq!(complete.pixels().len(), 2);
        assert_eq!(
            render_luminescence_camera(
                cx, &field, &view, LuminescenceCameraSettings { max_cell_visits: 3, ..base },
            ),
            Err(LuminescenceCameraError::Transport(VolumeTransportError::Budget)),
        );
        assert_eq!(
            render_luminescence_camera(
                cx, &field, &view, LuminescenceCameraSettings { max_samples: 1, ..base },
            ),
            Err(LuminescenceCameraError::Image(DvrError::SampleBudgetExceeded {
                requested: 2, limit: 1,
            })),
        );
        let missed = render_luminescence_camera(
            cx, &field, &camera(3.0, 4.0, 0.1, 2.0, 0.0),
            LuminescenceCameraSettings { max_cell_visits: 0, ..base },
        ).unwrap();
        for pixel in missed.pixels() {
            assert_eq!(pixel.emitted_radiance_w_m2_sr(), &[0.0]);
            assert_eq!(pixel.transmittance(), &[1.0]);
        }
        let diffraction = camera(-1.0, 0.0, 0.1, 2.0, 0.1)
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
        assert_eq!(
            render_luminescence_camera(cx, &field, &diffraction, base),
            Err(LuminescenceCameraError::Camera(CameraError::DiffractionRequiresSpectralPath)),
        );
        let limited_medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 0.5]),
            MassAbsorptionSpectrum::constant(0.0).unwrap(),
            None,
            AbsorbingGridBudget { max_grid_cells: 2, max_steps_per_segment: 1 },
        ).unwrap();
        let limited = LuminescenceGrid::try_new(
            cx, &limited_medium, SOURCE, [550.0], &[(1, [1.0])], 1,
        ).unwrap();
        assert_eq!(
            render_luminescence_camera(cx, &limited, &view, base),
            Err(LuminescenceCameraError::Transport(VolumeTransportError::Budget)),
        );
        gate.request();
        assert_eq!(
            render_luminescence_camera(cx, &field, &view, base),
            Err(LuminescenceCameraError::Transport(VolumeTransportError::Cancelled)),
        );
        assert_eq!(
            complete.approximate_emission_xyz(cx),
            Err(LuminescenceCameraError::Transport(VolumeTransportError::Cancelled)),
        );
    });
}

#[test]
fn g1_infrared_nodes_remain_physical_and_do_not_gain_a_visible_observer() {
    with_cx(|cx, _| {
        let density = [0.0];
        let medium = grid(
            cx, &density, [1; 3], [-0.5, -0.5, 0.0], [1.0; 3],
            MassAbsorptionSpectrum::source_table(&[(800.0, 0.0), (1000.0, 0.0)], SOURCE)
                .unwrap(),
        );
        let field = LuminescenceGrid::try_new(cx, &medium, SOURCE, [900.0], &[(0, [FOUR_PI])], 1)
            .unwrap();
        let image = render_luminescence_camera(
            cx, &field, &camera(-1.0, 0.0, 0.1, 2.0, 0.0), settings([1, 1], 1),
        ).unwrap();
        close(image.pixels()[0].emitted_radiance_w_m2_sr()[0], 1.0);
        assert_eq!(image.wavelengths_nm(), &[900.0]);
        assert_eq!(
            image.approximate_emission_xyz(cx),
            Err(LuminescenceCameraError::ObserverWavelength { lane: 0 }),
        );
        close(image.pixels()[0].composite([2.0]).unwrap()[0], 3.0);
    });
}
