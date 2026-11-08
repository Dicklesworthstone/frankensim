//! End-to-end perspective, coverage and finite-aperture volume images.

use crate::camera::{Aperture, CameraError, CameraProjection, FocusedDiffraction, PhysicalCamera};
use crate::volumes::opacity::{
    TransferCameraError, TransferCameraSettings, TransferRayBudget, TransferRayError,
    TransferRayGrid, render_transfer_camera,
};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::{DvrError, TransferFunction, TransferPoint, VolumeGrid};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
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

fn field_budget() -> TransferRayBudget {
    TransferRayBudget {
        max_grid_cells: 128,
        max_steps_per_ray: 128,
    }
}

fn settings(resolution: [usize; 2], samples_per_pixel: u32) -> TransferCameraSettings {
    TransferCameraSettings {
        resolution,
        samples_per_pixel,
        seed: 71,
        max_pixels: 128,
        max_samples: 100_000,
        max_cell_visits: 1_000_000,
    }
}

fn constant(extinction: f64, source_rgb: [f64; 3]) -> TransferFunction {
    TransferFunction::new(vec![TransferPoint {
        scalar: 0.0,
        extinction,
        source_rgb,
    }])
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs() + f64::MIN_POSITIVE,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_perspective_image_matches_the_analytic_slab_at_every_pixel() {
    with_cx(|cx, _| {
        let values = [0.0; 16];
        let grid = VolumeGrid::new([2, 2, 4], &values, [-2.0, -1.0, 0.0], [2.0, 1.0, 0.25]);
        let source = [2.0, 4.0, 7.0];
        let transfer = constant(1.7, source);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, field_budget()).unwrap();
        let half_tan = 0.2;
        let image = render_transfer_camera(
            cx,
            &medium,
            &camera(3.0, 0.0, half_tan, 3.0, 0.0),
            settings([3, 2], 1),
        )
        .unwrap();
        assert_eq!(image.len(), 6);
        let background = [1.0, 3.0, 6.0];
        for py in 0..2 {
            for px in 0..3 {
                let dx = (2.0 * (px as f64 + 0.5) / 3.0 - 1.0) * 1.5 * half_tan;
                let dy = (1.0 - 2.0 * (py as f64 + 0.5) / 2.0) * half_tan;
                let depth = 1.7 * (1.0 + dx * dx + dy * dy).sqrt();
                let transmission = (-depth).exp();
                let opacity = -(-depth).exp_m1();
                let pixel = image[px + 3 * py];
                close(pixel.optical_depth(), depth);
                close(pixel.transmittance(), transmission);
                close(pixel.opacity(), opacity);
                let composite = pixel.composite(background).unwrap();
                for channel in 0..3 {
                    close(pixel.emission_rgb()[channel], source[channel] * opacity);
                    close(
                        composite[channel],
                        source[channel] * opacity + background[channel] * transmission,
                    );
                }
            }
        }
    });
}

#[test]
fn g1_camera_image_rows_follow_the_declared_up_direction() {
    with_cx(|cx, _| {
        let values = [0.0, 1.0];
        let grid = VolumeGrid::new([1, 2, 1], &values, [-1.0, -1.0, 0.0], [2.0, 1.0, 1.0]);
        let bottom = [3.0, 0.0, 1.0];
        let top = [0.0, 2.0, 4.0];
        let transfer = TransferFunction::new(vec![
            TransferPoint {
                scalar: 0.0,
                extinction: 800.0,
                source_rgb: bottom,
            },
            TransferPoint {
                scalar: 1.0,
                extinction: 800.0,
                source_rgb: top,
            },
        ])
        .unwrap();
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, field_budget()).unwrap();
        let image = render_transfer_camera(
            cx,
            &medium,
            &camera(3.0, 0.0, 0.2, 3.0, 0.0),
            settings([1, 2], 1),
        )
        .unwrap();
        assert_eq!(image[0].emission_rgb(), top);
        assert_eq!(image[1].emission_rgb(), bottom);
    });
}

#[test]
fn g1_defocused_opaque_strip_matches_circular_pupil_coverage_and_refocuses() {
    with_cx(|cx, _| {
        let values = [0.0];
        // A thin, opaque emitting strip halfway between pupil and focal plane.
        // Its projected half-width is one quarter of the circular pupil radius.
        let grid = VolumeGrid::new([1; 3], &values, [-0.125, -2.0, 2.0], [0.25, 4.0, 1e-8]);
        let source = [2.0, 4.0, 8.0];
        let transfer = constant(1e12, source);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, field_budget()).unwrap();
        let defocused = camera(0.0, 2.0, 1e-10, 4.0, 1.0);
        let sampling = settings([1, 1], 32_768);
        let pixel = render_transfer_camera(cx, &medium, &defocused, sampling).unwrap()[0];
        let q = 0.25_f64;
        let covered = 2.0 / core::f64::consts::PI * (q.asin() + q * (1.0 - q * q).sqrt());
        // Independent disk-strip area, not the renderer's sampled ray arithmetic.
        assert!((pixel.opacity() - covered).abs() < 0.015);
        assert!((pixel.transmittance() - (1.0 - covered)).abs() < 0.015);
        close(pixel.optical_depth(), -pixel.transmittance().ln());
        assert!(pixel.transmittance() > 0.5); // Mean opaque depth would wrongly give zero.
        let background = [10.0, 7.0, 4.0];
        let composite = pixel.composite(background).unwrap();
        for channel in 0..3 {
            close(
                pixel.emission_rgb()[channel],
                source[channel] * pixel.opacity(),
            );
            close(
                composite[channel],
                source[channel] * pixel.opacity() + background[channel] * pixel.transmittance(),
            );
        }
        let focused = camera(0.0, 2.0, 1e-10, 2.0, 1.0);
        let sharp = render_transfer_camera(cx, &medium, &focused, settings([1, 1], 64)).unwrap()[0];
        assert_eq!(sharp.transmittance(), 0.0);
        assert_eq!(sharp.opacity(), 1.0);
        assert_eq!(sharp.emission_rgb(), source);
        // A shorter repeated sample sequence proves deterministic complete-image replay.
        let repeat = settings([1, 1], 256);
        assert_eq!(
            render_transfer_camera(cx, &medium, &defocused, repeat).unwrap(),
            render_transfer_camera(cx, &medium, &defocused, repeat).unwrap(),
        );
    });
}

#[test]
fn g0_sampled_thin_opacity_survives_unit_transmission_rounding() {
    with_cx(|cx, _| {
        let values = [0.0];
        let grid = VolumeGrid::new([1; 3], &values, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]);
        let transfer = constant(1e-20, [1e12; 3]);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, field_budget()).unwrap();
        let pixel = render_transfer_camera(
            cx,
            &medium,
            &camera(-1.0, 0.0, 0.01, 2.0, 0.0),
            settings([1, 1], 64),
        )
        .unwrap()[0];
        assert_eq!(pixel.transmittance(), 1.0);
        assert!(pixel.opacity() >= 1e-20 && pixel.opacity() < 1.001e-20);
        close(pixel.optical_depth(), pixel.opacity());
        for emitted in pixel.emission_rgb() {
            close(emitted, 1e12 * pixel.opacity());
        }
    });
}

#[test]
fn g0_finite_bright_means_and_backgrounds_survive_extreme_opacity() {
    with_cx(|cx, _| {
        let values = [0.0];
        let grid = VolumeGrid::new([1; 3], &values, [-1.0, -1.0, 0.0], [2.0; 3]);
        let transfer = constant(f64::MAX, [1e308; 3]);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, field_budget()).unwrap();
        let view = camera(-1.0, 0.0, 0.001, 2.0, 0.0);
        let pixel = render_transfer_camera(cx, &medium, &view, settings([1, 1], 64)).unwrap()[0];
        // Summing these finite samples before dividing would overflow.
        assert_eq!(pixel.emission_rgb(), [1e308; 3]);
        assert_eq!(pixel.transmittance(), 0.0);
        assert_eq!(pixel.optical_depth(), f64::INFINITY);

        let dark_grid = VolumeGrid::new([1; 3], &values, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]);
        let dark_transfer = constant(800.0, [0.0; 3]);
        let dark =
            TransferRayGrid::try_new(cx, &dark_grid, &dark_transfer, field_budget()).unwrap();
        let pixel = render_transfer_camera(cx, &dark, &view, settings([1, 1], 64)).unwrap()[0];
        assert_eq!(pixel.transmittance(), 0.0);
        assert!(pixel.optical_depth().is_finite());
        let maximum_depth = 800.0 * (1.0 + 2e-6_f64).sqrt();
        assert!(
            pixel.optical_depth() >= 800.0 - 1e-10
                && pixel.optical_depth() <= maximum_depth + 1e-10
        );
        let lower = (1e300_f64.ln() - maximum_depth).exp();
        let upper = (1e300_f64.ln() - 800.0).exp();
        for light in pixel.composite([1e300; 3]).unwrap() {
            assert!(light >= lower * (1.0 - 1e-11) && light <= upper * (1.0 + 1e-11));
        }
    });
}

#[test]
fn g4_image_limits_misses_camera_refusals_and_cancellation_are_transactional() {
    with_cx(|cx, gate| {
        let values = [0.0];
        let grid = VolumeGrid::new([1; 3], &values, [-2.0, -2.0, 0.0], [4.0, 4.0, 1.0]);
        let transfer = constant(1.0, [3.0; 3]);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, field_budget()).unwrap();
        let view = camera(-1.0, 0.0, 0.1, 2.0, 0.0);
        let base = settings([2, 1], 1);
        assert_eq!(
            render_transfer_camera(
                cx,
                &medium,
                &view,
                TransferCameraSettings {
                    max_pixels: 1,
                    ..base
                }
            ),
            Err(TransferCameraError::Image(DvrError::PixelBudgetExceeded {
                requested: 2,
                limit: 1
            }))
        );
        assert_eq!(
            render_transfer_camera(
                cx,
                &medium,
                &view,
                TransferCameraSettings {
                    max_samples: 1,
                    ..base
                }
            ),
            Err(TransferCameraError::Image(DvrError::SampleBudgetExceeded {
                requested: 2,
                limit: 1
            }))
        );
        assert_eq!(
            render_transfer_camera(
                cx,
                &medium,
                &view,
                TransferCameraSettings {
                    samples_per_pixel: 0,
                    ..base
                }
            ),
            Err(TransferCameraError::Image(DvrError::InvalidSamplingBudget))
        );
        // The second pixel must not renew the whole-image traversal allowance.
        assert_eq!(
            render_transfer_camera(
                cx,
                &medium,
                &view,
                TransferCameraSettings {
                    max_cell_visits: 1,
                    ..base
                }
            ),
            Err(TransferCameraError::Ray(TransferRayError::Traversal(
                VolumeTransportError::Budget
            )))
        );
        assert_eq!(
            render_transfer_camera(
                cx,
                &medium,
                &view,
                TransferCameraSettings {
                    max_cell_visits: 2,
                    ..base
                }
            )
            .unwrap()
            .len(),
            2
        );
        let misses = render_transfer_camera(
            cx,
            &medium,
            &camera(-1.0, -2.0, 0.1, 2.0, 0.0),
            TransferCameraSettings {
                max_cell_visits: 0,
                ..base
            },
        )
        .unwrap();
        for clear in misses {
            assert_eq!(clear.emission_rgb(), [0.0; 3]);
            assert_eq!(clear.opacity(), 0.0);
            assert_eq!(clear.transmittance(), 1.0);
            assert_eq!(clear.composite([7.0; 3]).unwrap(), [7.0; 3]);
        }
        let diffraction = camera(-1.0, 0.0, 0.1, 2.0, 0.1)
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
        assert_eq!(
            render_transfer_camera(cx, &medium, &diffraction, base),
            Err(TransferCameraError::Camera(
                CameraError::DiffractionRequiresSpectralPath
            ))
        );
        gate.request();
        assert_eq!(
            render_transfer_camera(cx, &medium, &view, base),
            Err(TransferCameraError::Ray(TransferRayError::Traversal(
                VolumeTransportError::Cancelled
            )))
        );
    });
}
