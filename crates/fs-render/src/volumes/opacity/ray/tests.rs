//! Physical ray-length, source-ordering and camera-use regressions.

use crate::camera::{Aperture, CameraProjection, LensSample, PhysicalCamera};
use crate::volumes::opacity::{
    TransferRayBudget, TransferRayError, TransferRayGrid, render_transfer_transport,
};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::{
    AbsorptionImageSettings, DvrError, Ray, TransferFunction, TransferPoint, VolumeGrid,
};
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

fn budget() -> TransferRayBudget {
    TransferRayBudget {
        max_grid_cells: 4096,
        max_steps_per_ray: 128,
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
        (actual - expected).abs() <= 2e-12 * expected.abs() + f64::MIN_POSITIVE,
        "{actual:e} != {expected:e}",
    );
}

#[test]
fn g1_perspective_camera_rays_use_physical_slab_length_and_linear_background() {
    with_cx(|cx, _| {
        let values = [0.0; 16];
        let grid = VolumeGrid::new([2, 2, 4], &values, [-2.0, -1.0, 0.0], [2.0, 1.0, 0.25]);
        let source = [2.0, 4.0, 7.0];
        let transfer = constant(1.7, source);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 3.0),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.3).unwrap(),
            3.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let background = [1.0, 3.0, 6.0];
        for y in [-0.1, 0.1] {
            for x in [-0.2, 0.0, 0.2] {
                let camera_ray = camera
                    .generate_ray_from_tangent_offsets(
                        cx,
                        x,
                        y,
                        LensSample::try_new(0.5, 0.5).unwrap(),
                    )
                    .unwrap();
                let pixel = medium
                    .transport(
                        cx,
                        Ray {
                            origin: [
                                camera_ray.origin.x,
                                camera_ray.origin.y,
                                camera_ray.origin.z,
                            ],
                            dir: [camera_ray.dir.x, camera_ray.dir.y, camera_ray.dir.z],
                            t0: 0.0,
                            t1: f64::INFINITY,
                        },
                    )
                    .unwrap();
                // Every camera ray traverses the complete unit-thickness slab.
                let depth = 1.7 * (1.0 + x * x + y * y).sqrt();
                let transmission = (-depth).exp();
                let opacity = -(-depth).exp_m1();
                close(pixel.optical_depth(), depth);
                close(pixel.transmittance(), transmission);
                close(pixel.opacity(), opacity);
                let observed = pixel.composite(background).unwrap();
                for channel in 0..3 {
                    close(pixel.emission_rgb()[channel], source[channel] * opacity);
                    close(
                        observed[channel],
                        source[channel] * opacity + background[channel] * transmission,
                    );
                }
            }
        }
    });
}

#[test]
fn g1_oblique_anisotropic_layers_preserve_source_order_in_both_directions() {
    with_cx(|cx, _| {
        let values = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let grid = VolumeGrid::new([2, 2, 2], &values, [5.0, -4.0, 7.0], [2.0, 3.0, 0.5]);
        let sources = [[8.0, 1.0, 0.0], [0.0, 3.0, 5.0]];
        let transfer = TransferFunction::new(vec![
            TransferPoint {
                scalar: 0.0,
                extinction: 0.7,
                source_rgb: sources[0],
            },
            TransferPoint {
                scalar: 1.0,
                extinction: 1.3,
                source_rgb: sources[1],
            },
        ])
        .unwrap();
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        let layer_length = 0.5 * 1.25_f64.sqrt();
        let depths = [0.7 * layer_length, 1.3 * layer_length];
        for (origin, dir, order) in [
            ([5.5, -3.0, 6.0], [0.3, 0.4, 1.0], [0, 1]),
            ([6.4, -1.8, 9.0], [-0.3, -0.4, -1.0], [1, 0]),
        ] {
            let pixel = medium
                .transport(
                    cx,
                    Ray {
                        origin,
                        dir,
                        t0: 0.0,
                        t1: 4.0,
                    },
                )
                .unwrap();
            let [front, back] = order;
            let front_alpha = -(-depths[front]).exp_m1();
            let back_weight = (-depths[front]).exp() * -(-depths[back]).exp_m1();
            close(pixel.optical_depth(), depths[0] + depths[1]);
            for channel in 0..3 {
                close(
                    pixel.emission_rgb()[channel],
                    sources[front][channel] * front_alpha + sources[back][channel] * back_weight,
                );
            }
        }
    });
}

#[test]
fn g1_clipped_inside_segments_are_invariant_to_ray_parameter_scale() {
    with_cx(|cx, _| {
        let values = [-1.0, 1.0];
        let grid = VolumeGrid::new([1, 1, 2], &values, [-1.0, -1.0, 0.0], [2.0, 3.0, 1.0]);
        let transfer = TransferFunction::new(vec![
            TransferPoint {
                scalar: -1.0,
                extinction: 0.4,
                source_rgb: [3.0; 3],
            },
            TransferPoint {
                scalar: 1.0,
                extinction: 1.6,
                source_rgb: [7.0; 3],
            },
        ])
        .unwrap();
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        for scale in [1.0, 2.0, 1e-200, 1e200] {
            let pixel = medium
                .transport(
                    cx,
                    Ray {
                        origin: [0.0, 0.0, -1.0],
                        dir: [0.0, 0.0, scale],
                        t0: 1.5 / scale,
                        t1: 2.5 / scale,
                    },
                )
                .unwrap();
            // Clipped endpoints z=0.5 and z=1.5 retain half of each layer.
            close(pixel.optical_depth(), 1.0);
            close(pixel.transmittance(), (-1.0_f64).exp());
            let expected =
                3.0 * -(-0.2_f64).exp_m1() + 7.0 * (-0.2_f64).exp() * -(-0.8_f64).exp_m1();
            for emission in pixel.emission_rgb() {
                close(emission, expected);
            }
        }
    });
}

#[test]
fn g1_axis_aligned_queries_match_existing_orthographic_pixels() {
    with_cx(|cx, _| {
        let values = [
            -1.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0, -1.0, 0.5, -1.0, 1.0, 0.0,
        ];
        let grid = VolumeGrid::new([2, 2, 3], &values, [1e6, -3.0, 9.0], [0.5, 1.0, 0.25]);
        let transfer = TransferFunction::new(vec![
            TransferPoint {
                scalar: -1.0,
                extinction: 0.2,
                source_rgb: [0.0, 4.0, 2.0],
            },
            TransferPoint {
                scalar: 1.0,
                extinction: 2.0,
                source_rgb: [3.0, 0.0, 7.0],
            },
        ])
        .unwrap();
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        let image = render_transfer_transport(
            &grid,
            &transfer,
            AbsorptionImageSettings {
                resolution: [2, 2],
                max_pixels: 4,
                max_grid_cells: values.len(),
                max_cell_visits: 12,
            },
        )
        .unwrap();
        for y in 0..2 {
            for x in 0..2 {
                let pixel = medium
                    .transport(
                        cx,
                        Ray {
                            origin: [1e6 + (x as f64 + 0.5) * 0.5, -3.0 + y as f64 + 0.5, 12.0],
                            dir: [0.0, 0.0, -1.0],
                            t0: 0.0,
                            t1: f64::INFINITY,
                        },
                    )
                    .unwrap();
                assert_eq!(pixel, image[x + 2 * y]);
            }
        }
    });
}

#[test]
fn g0_tiny_luminous_cells_survive_underflowed_opacity_in_both_apis() {
    with_cx(|cx, _| {
        let values = [0.0];
        let grid = VolumeGrid::new([1; 3], &values, [0.0; 3], [1.0, 1.0, 1e-200]);
        let transfer = constant(1e-200, [1e300; 3]);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        let pixel = medium
            .transport(
                cx,
                Ray {
                    origin: [0.5, 0.5, 0.0],
                    dir: [0.0, 0.0, 1.0],
                    t0: 0.0,
                    t1: f64::INFINITY,
                },
            )
            .unwrap();
        assert_eq!(pixel.transmittance(), 1.0);
        assert_eq!(pixel.opacity(), 0.0); // True 1e-400 opacity cannot be represented.
        for emission in pixel.emission_rgb() {
            close(emission, 1e-100);
        }
        let column = render_transfer_transport(
            &grid,
            &transfer,
            AbsorptionImageSettings {
                resolution: [1, 1],
                max_pixels: 1,
                max_grid_cells: 1,
                max_cell_visits: 1,
            },
        )
        .unwrap();
        assert_eq!(pixel, column[0]);
    });
}

#[test]
fn g0_bright_source_behind_an_opaque_foreground_remains_representable() {
    with_cx(|cx, _| {
        let values = [0.0, 1.0];
        let grid = VolumeGrid::new([1, 1, 2], &values, [0.0; 3], [1.0; 3]);
        let transfer = TransferFunction::new(vec![
            TransferPoint {
                scalar: 0.0,
                extinction: 800.0,
                source_rgb: [0.0; 3],
            },
            TransferPoint {
                scalar: 1.0,
                extinction: 1.0,
                source_rgb: [1e300; 3],
            },
        ])
        .unwrap();
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        let pixel = medium
            .transport(
                cx,
                Ray {
                    origin: [0.5, 0.5, -1.0],
                    dir: [0.0, 0.0, 1.0],
                    t0: 0.0,
                    t1: f64::INFINITY,
                },
            )
            .unwrap();
        assert_eq!(pixel.transmittance(), 0.0);
        for emission in pixel.emission_rgb() {
            close(
                emission,
                (1e300_f64.ln() - 800.0).exp() * -(-1.0_f64).exp_m1(),
            );
        }
    });
}

#[test]
fn g0_zero_scalar_is_not_vacuum_but_misses_and_empty_segments_are() {
    with_cx(|cx, _| {
        let values = [0.0];
        let grid = VolumeGrid::new([1; 3], &values, [0.0; 3], [1.0; 3]);
        let transfer = constant(2.0, [4.0; 3]);
        let medium = TransferRayGrid::try_new(cx, &grid, &transfer, budget()).unwrap();
        let ray = Ray {
            origin: [0.5, 0.5, -1.0],
            dir: [0.0, 0.0, 1.0],
            t0: 0.0,
            t1: 3.0,
        };
        close(
            medium.transport(cx, ray).unwrap().transmittance(),
            (-2.0_f64).exp(),
        );
        for empty in [
            Ray {
                origin: [2.0, 0.5, -1.0],
                ..ray
            },
            Ray {
                t0: 1.5,
                t1: 1.5,
                ..ray
            },
        ] {
            let clear = medium.transport(cx, empty).unwrap();
            assert_eq!(clear.emission_rgb(), [0.0; 3]);
            assert_eq!(clear.opacity(), 0.0);
            assert_eq!(clear.composite([3.0; 3]).unwrap(), [3.0; 3]);
        }
    });
}

#[test]
fn g4_invalid_inputs_budgets_and_cancellation_return_no_partial_ray() {
    with_cx(|cx, gate| {
        let values = [0.0; 2];
        let grid = VolumeGrid::new([1, 1, 2], &values, [0.0; 3], [1.0; 3]);
        let transfer = constant(f64::MAX, [1.0; 3]);
        assert!(matches!(
            TransferRayGrid::try_new(
                cx,
                &grid,
                &transfer,
                TransferRayBudget {
                    max_grid_cells: 1,
                    ..budget()
                }
            ),
            Err(TransferRayError::Transfer(
                DvrError::GridCellBudgetExceeded { .. }
            ))
        ));
        let mut invalid_values = [0.0; 1026];
        invalid_values[1025] = f64::NAN;
        let invalid_grid = VolumeGrid::new([1, 1, 1026], &invalid_values, [0.0; 3], [1.0; 3]);
        assert!(matches!(
            TransferRayGrid::try_new(cx, &invalid_grid, &transfer, budget()),
            Err(TransferRayError::Transfer(DvrError::NonFiniteGridValue {
                index: 1025
            }))
        ));
        let collapsed = VolumeGrid::new([1, 1, 2], &values, [0.0, 0.0, 1e16], [1.0; 3]);
        assert!(matches!(
            TransferRayGrid::try_new(cx, &collapsed, &transfer, budget()),
            Err(TransferRayError::Transfer(DvrError::InvalidGrid))
        ));
        let medium = TransferRayGrid::try_new(
            cx,
            &grid,
            &transfer,
            TransferRayBudget {
                max_steps_per_ray: 1,
                ..budget()
            },
        )
        .unwrap();
        let ray = Ray {
            origin: [0.5, 0.5, -1.0],
            dir: [0.0, 0.0, 1.0],
            t0: 0.0,
            t1: f64::INFINITY,
        };
        // Even an opaque first cell cannot hide a budget refusal in the second.
        assert_eq!(
            medium.transport(cx, ray),
            Err(TransferRayError::Traversal(VolumeTransportError::Budget))
        );
        for invalid in [
            Ray {
                dir: [0.0; 3],
                ..ray
            },
            Ray {
                origin: [f64::NAN, 0.5, 0.0],
                ..ray
            },
            Ray {
                t0: 2.0,
                t1: 1.0,
                ..ray
            },
        ] {
            assert_eq!(
                medium.transport(cx, invalid),
                Err(TransferRayError::Traversal(VolumeTransportError::Ray))
            );
        }
        gate.request();
        assert_eq!(
            medium.transport(cx, ray),
            Err(TransferRayError::Traversal(VolumeTransportError::Cancelled))
        );
        assert!(matches!(
            TransferRayGrid::try_new(cx, &grid, &transfer, budget()),
            Err(TransferRayError::Traversal(VolumeTransportError::Cancelled))
        ));
    });
}
