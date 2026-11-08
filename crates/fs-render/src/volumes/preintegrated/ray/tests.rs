//! Independent physical ray references for continuous trilinear opacity.
use super::*;
use crate::volumes::TransferPoint;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn settings() -> PreintegratedRaySettings {
    PreintegratedRaySettings {
        optical_depth_tolerance: 1e-6,
        source_reconstruction_tolerance_rgb: [1e-5; 3],
        max_curve_subdivisions: 100_000,
        integration: PreintegrationSettings {
            source_tolerance: 1e-9,
            max_intervals: 500_000,
            max_quadrature_evaluations: 5_000_000,
        },
    }
}

fn budget() -> NodeRayBudget {
    NodeRayBudget {
        max_grid_nodes: 100_000,
        max_steps_per_ray: 128,
    }
}

fn transfer(knots: &[(f64, f64, [f64; 3])]) -> TransferFunction {
    TransferFunction::new(
        knots
            .iter()
            .map(|&(scalar, extinction, source_rgb)| TransferPoint {
                scalar,
                extinction,
                source_rgb,
            })
            .collect(),
    )
    .unwrap()
}

fn nodes(n: usize, scalar: impl Fn(f64, f64, f64) -> f64) -> Vec<f64> {
    let mut values = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                values.push(scalar(
                    x as f64 / (n - 1) as f64,
                    y as f64 / (n - 1) as f64,
                    z as f64 / (n - 1) as f64,
                ));
            }
        }
    }
    values
}

fn diagonal() -> Ray {
    Ray {
        origin: [0.0; 3],
        dir: [1.0; 3],
        t0: 0.0,
        t1: 1.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.17e}"
    );
}

fn within_bounds(actual: PreintegratedRayPixel, tau: f64, emission: [f64; 3]) {
    close(
        actual.transport.optical_depth(),
        tau,
        actual.optical_depth_error + 2e-10,
    );
    for channel in 0..3 {
        close(
            actual.transport.emission_rgb()[channel],
            emission[channel],
            actual.source_reconstruction_error_rgb[channel]
                + actual.transport.source_error_rgb()[channel]
                + 2e-10,
        );
    }
}

#[test]
fn cubic_depth_and_constant_source_preserve_the_trilinear_field_under_grid_refinement() {
    with_cx(|cx, _| {
        let source = [2.0, 4.0, 8.0];
        let tf = transfer(&[(0.0, 0.0, source), (1.0, 1.0, source)]);
        let tau = 3.0_f64.sqrt() / 4.0;
        let emission = source.map(|s| s * -(-tau).exp_m1());
        for n in [2, 5] {
            let values = nodes(n, |x, y, z| x * y * z);
            let grid = NodeVolumeGrid {
                dims: [n; 3],
                data: &values,
                origin: [0.0; 3],
                spacing: [1.0 / (n - 1) as f64; 3],
            };
            let owner = PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()).unwrap();
            let pixel = owner.transport(cx, diagonal(), settings()).unwrap();
            within_bounds(pixel, tau, emission);
            assert!(pixel.optical_depth_error <= settings().optical_depth_tolerance);
            for c in 0..3 {
                assert!(
                    pixel.source_reconstruction_error_rgb[c]
                        <= settings().source_reconstruction_tolerance_rgb[c]
                );
                close(
                    pixel.transport.composite(source).unwrap()[c],
                    source[c],
                    2e-12,
                );
            }
        }
    });
}

fn hidden_scalar(t: f64) -> f64 {
    t * (t - 0.5) * (t - 1.0)
}

fn scalar_root(level: f64, mut low: f64, mut high: f64, increasing: bool) -> f64 {
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if (hidden_scalar(middle) < level) == increasing {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

// Independent exact primitive on each piece of the triangular transfer peak.
// The six crossing locations are roots of the explicitly authored cubic,
// unrelated to the native Bézier controls, subdivision or chord integration.
fn hidden_peak_depth() -> f64 {
    let turn = (3.0 - 3.0_f64.sqrt()) / 6.0;
    let rising = [0.038, 0.040, 0.042].map(|s| scalar_root(s, 0.0, turn, true));
    let falling = [0.038, 0.040, 0.042].map(|s| scalar_root(s, turn, 0.5, false));
    let primitive = |t: f64| t.powi(4) / 4.0 - t.powi(3) / 2.0 + t * t / 4.0;
    let pieces = [
        (rising[0], rising[1], 10_000.0, -380.0),
        (rising[1], rising[2], -10_000.0, 420.0),
        (falling[2], falling[1], -10_000.0, 420.0),
        (falling[1], falling[0], 10_000.0, -380.0),
    ];
    3.0_f64.sqrt()
        * pieces
            .into_iter()
            .map(|(a, b, m, q)| m * (primitive(b) - primitive(a)) + q * (b - a))
            .sum::<f64>()
}

#[test]
fn interior_cubic_opacity_peak_survives_transparent_endpoints_and_midpoint() {
    with_cx(|cx, _| {
        let values = nodes(2, |x, y, z| x * (y - 0.5) * (z - 1.0));
        let grid = NodeVolumeGrid {
            dims: [2; 3],
            data: &values,
            origin: [0.0; 3],
            spacing: [1.0; 3],
        };
        let source = [2.0, 4.0, 8.0];
        let tf = transfer(&[
            (-0.5, 0.0, source),
            (0.038, 0.0, source),
            (0.040, 20.0, source),
            (0.042, 0.0, source),
            (0.5, 0.0, source),
        ]);
        for t in [0.0, 0.5, 1.0] {
            assert_eq!(hidden_scalar(t), 0.0);
        }
        let owner = PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()).unwrap();
        let mut controls = settings();
        controls.optical_depth_tolerance = 1e-5;
        controls.source_reconstruction_tolerance_rgb = [1e-3; 3];
        controls.integration.max_quadrature_evaluations = 0;
        let pixel = owner.transport(cx, diagonal(), controls).unwrap();
        let tau = hidden_peak_depth();
        assert!(tau > 0.1);
        assert!(pixel.transport.opacity() > 0.1);
        within_bounds(pixel, tau, source.map(|s| s * -(-tau).exp_m1()));
        assert!(pixel.optical_depth_error <= controls.optical_depth_tolerance);
    });
}

// Smooth distance-domain Simpson reference with an analytic foreground depth.
// No production transfer lookup, chord model or quadrature is used here.
fn cubic_source_reference(intervals: usize) -> [f64; 3] {
    let mut sum = [0.0; 3];
    let length = 3.0_f64.sqrt();
    for i in 0..=intervals {
        let t = i as f64 / intervals as f64;
        let scalar = t.powi(3);
        let source = [2.0 + 3.0 * scalar, 7.0 * scalar, 5.0];
        let attenuation = (-length * (t + t.powi(4) / 4.0)).exp();
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for c in 0..3 {
            sum[c] += weight * length * (1.0 + scalar) * source[c] * attenuation;
        }
    }
    sum.map(|value| value / (3.0 * intervals as f64))
}

#[test]
fn cubic_varying_source_and_foreground_attenuation_match_independent_distance_integral() {
    with_cx(|cx, _| {
        let values = nodes(2, |x, y, z| x * y * z);
        let grid = NodeVolumeGrid {
            dims: [2; 3],
            data: &values,
            origin: [0.0; 3],
            spacing: [1.0; 3],
        };
        let tf = transfer(&[(0.0, 1.0, [2.0, 0.0, 5.0]), (1.0, 2.0, [5.0, 7.0, 5.0])]);
        let owner = PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()).unwrap();
        let reference = cubic_source_reference(8192);
        for (a, b) in reference.into_iter().zip(cubic_source_reference(4096)) {
            close(a, b, 2e-11);
        }
        let mut coarse_controls = settings();
        coarse_controls.optical_depth_tolerance = 5e-3;
        coarse_controls.source_reconstruction_tolerance_rgb = [0.05; 3];
        let coarse = owner.transport(cx, diagonal(), coarse_controls).unwrap();
        let fine = owner.transport(cx, diagonal(), settings()).unwrap();
        let tau = 1.25 * 3.0_f64.sqrt();
        for pixel in [coarse, fine] {
            within_bounds(pixel, tau, reference);
        }
        assert!(fine.optical_depth_error < coarse.optical_depth_error);
        let largest_error = |pixel: PreintegratedRayPixel| {
            (0..3)
                .map(|c| (pixel.transport.emission_rgb()[c] - reference[c]).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(largest_error(fine) < largest_error(coarse));
        for c in 0..3 {
            assert!(
                fine.source_reconstruction_error_rgb[c]
                    <= settings().source_reconstruction_tolerance_rgb[c]
            );
        }
    });
}

#[test]
fn affine_field_clipping_nonunit_rays_and_unit_scaling_preserve_physical_depth() {
    with_cx(|cx, _| {
        let values = nodes(2, |x, y, z| 0.2 + 0.3 * x + 0.4 * y + 0.1 * z);
        let source = [2.0, 0.5, 5.0];
        // Along the anisotropic diagonal, ds=7*scale*du. The clipped interval
        // is u in [1/4,1], with sigma=(0.74+0.96*u)/scale.
        let tau = 7.0 * (0.74 * 0.75 + 0.48 * (1.0 - 0.25_f64.powi(2)));
        for scale in [0.125, 1.0, 16.0] {
            let origin = [7.0, -11.0, 5.0];
            let extent = [2.0 * scale, 3.0 * scale, 6.0 * scale];
            let grid = NodeVolumeGrid {
                dims: [2; 3],
                data: &values,
                origin,
                spacing: extent,
            };
            let tf = transfer(&[(0.0, 0.5 / scale, source), (1.0, 1.7 / scale, source)]);
            let owner = PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()).unwrap();
            let ray = Ray {
                origin: std::array::from_fn(|c| origin[c] - extent[c]),
                dir: extent.map(|value| 2.0 * value),
                t0: 0.625,
                t1: 1.5,
            };
            let mut controls = settings();
            controls.max_curve_subdivisions = 0;
            let pixel = owner.transport(cx, ray, controls).unwrap();
            close(pixel.transport.optical_depth(), tau, 5e-13);
            within_bounds(pixel, tau, source.map(|s| s * -(-tau).exp_m1()));
            let equivalent = Ray {
                origin: std::array::from_fn(|c| origin[c] + 0.25 * extent[c]),
                dir: extent,
                t0: 0.0,
                t1: 0.75,
            };
            let other = owner.transport(cx, equivalent, controls).unwrap();
            close(
                other.transport.optical_depth(),
                pixel.transport.optical_depth(),
                5e-13,
            );
            let miss = Ray {
                origin: [1e6; 3],
                dir: [1.0; 3],
                t0: 0.0,
                t1: 1.0,
            };
            let vacuum = owner.transport(cx, miss, controls).unwrap();
            assert_eq!(vacuum.transport.optical_depth(), 0.0);
            assert_eq!(vacuum.transport.emission_rgb(), [0.0; 3]);
            assert_eq!(vacuum.source_reconstruction_error_rgb, [0.0; 3]);
        }
    });
}

#[test]
fn axis_ray_matches_original_node_render_and_reversal_preserves_depth_not_color_order() {
    with_cx(|cx, _| {
        let values = nodes(3, |_, _, z| z);
        let grid = NodeVolumeGrid {
            dims: [3; 3],
            data: &values,
            origin: [0.0; 3],
            spacing: [0.5; 3],
        };
        let tf = transfer(&[(0.0, 0.1, [8.0, 0.0, 0.0]), (1.0, 4.0, [0.0, 8.0, 0.0])]);
        let owner = PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()).unwrap();
        let forward = owner
            .transport(
                cx,
                Ray {
                    origin: [0.5, 0.5, -1.0],
                    dir: [0.0, 0.0, 3.0],
                    t0: 0.0,
                    t1: 1.0,
                },
                settings(),
            )
            .unwrap();
        let reverse = owner
            .transport(
                cx,
                Ray {
                    origin: [0.5, 0.5, 2.0],
                    dir: [0.0, 0.0, -3.0],
                    t0: 0.0,
                    t1: 1.0,
                },
                settings(),
            )
            .unwrap();
        let orthographic = super::super::render_preintegrated_transfer(
            &grid,
            &tf,
            super::super::PreintegratedImageSettings {
                resolution: [1, 1],
                max_pixels: 1,
                max_grid_nodes: values.len(),
                max_segments: 2,
                integration: settings().integration,
            },
        )
        .unwrap()[0];
        close(reverse.transport.optical_depth(), 2.05, 1e-14);
        close(
            forward.transport.optical_depth(),
            reverse.transport.optical_depth(),
            1e-14,
        );
        for c in 0..3 {
            close(
                reverse.transport.emission_rgb()[c],
                orthographic.emission_rgb()[c],
                2e-9,
            );
        }
        assert!(forward.transport.emission_rgb()[0] > reverse.transport.emission_rgb()[0]);
        assert!(reverse.transport.emission_rgb()[1] > forward.transport.emission_rgb()[1]);
    });
}

#[test]
fn exhausted_curve_ray_and_integration_budgets_and_cancellation_refuse_partial_transport() {
    with_cx(|cx, gate| {
        let values = nodes(2, |x, y, z| x * y * z);
        let grid = NodeVolumeGrid {
            dims: [2; 3],
            data: &values,
            origin: [0.0; 3],
            spacing: [1.0; 3],
        };
        let tf = transfer(&[(0.0, 0.1, [0.0; 3]), (1.0, 4.0, [1.0; 3])]);
        let owner = PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()).unwrap();
        let mut controls = settings();
        controls.max_curve_subdivisions = 0;
        assert_eq!(
            owner.transport(cx, diagonal(), controls),
            Err(PreintegratedRayError::SubdivisionBudgetExceeded)
        );
        controls = settings();
        controls.integration.max_intervals = 0;
        assert_eq!(
            owner.transport(cx, diagonal(), controls),
            Err(PreintegratedRayError::Integration(
                PreintegrationError::TransferIntervalBudgetExceeded
            ))
        );
        controls = settings();
        controls.integration.max_quadrature_evaluations = 0;
        assert_eq!(
            owner.transport(cx, diagonal(), controls),
            Err(PreintegratedRayError::Integration(
                PreintegrationError::QuadratureBudgetExceeded
            ))
        );
        controls = settings();
        controls.optical_depth_tolerance = -1.0;
        assert_eq!(
            owner.transport(cx, diagonal(), controls),
            Err(PreintegratedRayError::InvalidReconstructionTolerance)
        );
        assert!(matches!(
            PreintegratedRayGrid::try_new(
                cx,
                &grid,
                &tf,
                NodeRayBudget {
                    max_grid_nodes: 7,
                    ..budget()
                }
            ),
            Err(PreintegratedRayError::Integration(
                PreintegrationError::Transfer(DvrError::GridCellBudgetExceeded { .. })
            ))
        ));

        let constant = [0.5; 27];
        let multi = NodeVolumeGrid {
            dims: [3; 3],
            data: &constant,
            origin: [0.0; 3],
            spacing: [0.5; 3],
        };
        let limited = PreintegratedRayGrid::try_new(
            cx,
            &multi,
            &tf,
            NodeRayBudget {
                max_steps_per_ray: 1,
                ..budget()
            },
        )
        .unwrap();
        let axis = Ray {
            origin: [0.25, 0.25, 0.0],
            dir: [0.0, 0.0, 1.0],
            t0: 0.0,
            t1: 1.0,
        };
        assert_eq!(
            limited.transport(cx, axis, settings()),
            Err(PreintegratedRayError::Traversal(
                VolumeTransportError::Budget
            ))
        );
        let unlimited = PreintegratedRayGrid::try_new(cx, &multi, &tf, budget()).unwrap();
        controls = settings();
        controls.integration.max_intervals = 1;
        assert_eq!(
            unlimited.transport(cx, axis, controls),
            Err(PreintegratedRayError::Integration(
                PreintegrationError::TransferIntervalBudgetExceeded
            ))
        );

        let mut invalid = [0.0; 8];
        invalid[7] = f64::NAN;
        let invalid = NodeVolumeGrid {
            dims: [2; 3],
            data: &invalid,
            origin: [0.0; 3],
            spacing: [1.0; 3],
        };
        assert!(matches!(
            PreintegratedRayGrid::try_new(cx, &invalid, &tf, budget()),
            Err(PreintegratedRayError::Integration(
                PreintegrationError::Transfer(DvrError::NonFiniteGridValue { index: 7 })
            ))
        ));
        // Earlier failures cannot consume the borrowed input or partially
        // publish an output; the same admitted field still transports normally.
        assert!(owner.transport(cx, diagonal(), settings()).is_ok());
        gate.request();
        assert_eq!(
            owner.transport(cx, diagonal(), settings()),
            Err(PreintegratedRayError::Traversal(
                VolumeTransportError::Cancelled
            ))
        );
        assert!(matches!(
            PreintegratedRayGrid::try_new(cx, &grid, &tf, budget()),
            Err(PreintegratedRayError::Traversal(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}
