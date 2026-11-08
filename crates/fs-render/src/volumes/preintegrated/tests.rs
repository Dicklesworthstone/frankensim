use super::*;
use crate::volumes::TransferPoint;

fn settings() -> PreintegrationSettings {
    PreintegrationSettings {
        source_tolerance: 1e-10,
        max_intervals: 10_000,
        max_quadrature_evaluations: 200_000,
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn image_settings(nodes: usize, pixels: [usize; 2], segments: u64) -> PreintegratedImageSettings {
    PreintegratedImageSettings {
        resolution: pixels,
        max_pixels: pixels[0] * pixels[1],
        max_grid_nodes: nodes,
        max_segments: segments,
        integration: settings(),
    }
}

#[test]
fn g1_homogeneous_and_linear_extinction_slabs_match_beer_lambert() {
    for (a, b, length) in [(0.7, 0.7, 2.0), (0.0, 4.0, 0.3), (4.0, 0.01, 3.0)] {
        let source = [2.0, 7.0, 0.25];
        let tf = transfer(&[(0.0, a, source), (1.0, b, source)]);
        let pixel = preintegrate_transfer_segment(
            &tf,
            0.0,
            1.0,
            length,
            PreintegrationSettings {
                max_quadrature_evaluations: 0,
                ..settings()
            },
        )
        .unwrap();
        let tau = 0.5 * (a + b) * length;
        close(pixel.optical_depth(), tau, 2e-15 * tau);
        close(pixel.transmittance(), (-tau).exp(), 2e-15);
        for (actual, expected) in pixel.emission_rgb().into_iter().zip(source) {
            close(actual, expected * -(-tau).exp_m1(), 2e-14);
        }
        // A source equal to the background conserves its radiance exactly in
        // the transfer equation, including a varying extinction coefficient.
        for (actual, expected) in pixel.composite(source).unwrap().into_iter().zip(source) {
            close(actual, expected, 2e-14);
        }
        assert_eq!(pixel.source_error_rgb(), [0.0; 3]);
    }
}

#[test]
fn g1_extreme_constant_profiles_match_cellwise_transport() {
    use crate::volumes::{AbsorptionImageSettings, VolumeGrid, render_transfer_transport};

    for (sigma, length, source) in [
        (1e-200, 1e-200, [1e200, 2e200, 0.0]),
        (1e-300, 1e-20, [1e300, 2e300, 0.0]),
        (1e308, 4.0, [1.0, 2.0, 0.0]),
    ] {
        let tf = transfer(&[(0.0, sigma, source)]);
        let pixel = preintegrate_transfer_segment(
            &tf,
            0.0,
            0.0,
            length,
            PreintegrationSettings {
                max_quadrature_evaluations: 0,
                ..settings()
            },
        )
        .unwrap();
        let grid = VolumeGrid {
            dims: [1; 3],
            data: &[0.0],
            origin: [0.0; 3],
            cell: [1.0, 1.0, length],
        };
        let cellwise = render_transfer_transport(
            &grid,
            &tf,
            AbsorptionImageSettings {
                resolution: [1; 2],
                max_pixels: 1,
                max_grid_cells: 1,
                max_cell_visits: 1,
            },
        )
        .unwrap()[0];
        assert_eq!(pixel.source_error_rgb(), [0.0; 3]);
        for channel in 0..3 {
            let expected = if length < 1.0 {
                // Thin Beer--Lambert limit, multiplied in an independent
                // order so the original positive optical factors survive.
                (source[channel] * sigma) * length
            } else {
                source[channel]
            };
            close(pixel.emission_rgb()[channel], expected, 8e-13 * expected);
            close(cellwise.emission_rgb()[channel], expected, 8e-13 * expected);
        }
        assert_eq!(pixel.opacity(), cellwise.opacity());
        assert_eq!(pixel.transmittance(), cellwise.transmittance());
        assert_eq!(pixel.optical_depth(), cellwise.optical_depth());
    }
}

#[test]
fn g1_subnormal_extinction_is_scaled_before_its_mean_is_rounded() {
    for bits in [1, 3, 7] {
        let sigma = f64::from_bits(bits);
        let length = 1e300;
        let tf = transfer(&[(0.0, 0.0, [1.0; 3]), (1.0, sigma, [1.0; 3])]);
        let pixel = preintegrate_transfer_segment(&tf, 0.0, 1.0, length, settings()).unwrap();
        // The exact triangular area can be represented after multiplication
        // by length, even though sigma/2 would round by up to 100 percent.
        let expected_depth = (sigma * length) / 2.0;
        close(
            pixel.optical_depth(),
            expected_depth,
            2e-15 * expected_depth,
        );
        close(
            pixel.emission_rgb()[0],
            expected_depth,
            2e-15 * expected_depth,
        );
    }
}

#[test]
fn g1_underflowed_linear_extinction_retains_the_source_moment() {
    let sigma = 1e-200;
    let length = 1e-200;
    let source = 1e200;
    let request = PreintegrationSettings {
        max_quadrature_evaluations: 0,
        ..settings()
    };
    for (front, back, moment) in [(0.0, sigma, 1.0 / 3.0), (sigma, 0.0, 1.0 / 6.0)] {
        let tf = transfer(&[(0.0, front, [0.0; 3]), (1.0, back, [source; 3])]);
        let pixel = preintegrate_transfer_segment(&tf, 0.0, 1.0, length, request).unwrap();
        // Integral of u*u or (1-u)*u over [0,1]; attenuation corrections
        // are smaller than binary64 relative precision by hundreds of orders.
        let expected = ((source * sigma) * length) * moment;
        close(pixel.emission_rgb()[0], expected, 8e-13 * expected);
        assert_eq!(pixel.optical_depth(), 0.0);
        assert_eq!(pixel.opacity(), 0.0);
        assert_eq!(pixel.source_error_rgb(), [0.0; 3]);

        let mut split = RayIntegral::default();
        let mut work = Work::new(request).unwrap();
        integrate_segment(&tf, 0.0, 0.5, 0.5 * length, &mut work, &mut split).unwrap();
        integrate_segment(&tf, 0.5, 1.0, 0.5 * length, &mut work, &mut split).unwrap();
        close(split.finish().emission_rgb()[0], expected, 8e-13 * expected);
    }
}

#[test]
fn g1_overflowed_homogeneous_depth_retains_the_back_source_tail() {
    let tf = transfer(&[
        (0.0, 1e308, [0.0, 1.0, 1e-300]),
        (1.0, 1e308, [1e308, 0.0, 8e307]),
    ]);
    let request = PreintegrationSettings {
        max_quadrature_evaluations: 0,
        ..settings()
    };
    for intervals in [1, 2, 4, 8] {
        let mut ray = RayIntegral::default();
        let mut work = Work::new(request).unwrap();
        for i in 0..intervals {
            integrate_segment(
                &tf,
                f64::from(i) / f64::from(intervals),
                f64::from(i + 1) / f64::from(intervals),
                4.0 / f64::from(intervals),
                &mut work,
                &mut ray,
            )
            .unwrap();
        }
        let pixel = ray.finish();
        // Homogeneous transport has source moment 1/(sigma*length) in
        // this opaque limit, not zero. The red and blue front sources alone
        // cannot represent the radiance that escapes the near-front layer.
        for (actual, expected) in pixel.emission_rgb().into_iter().zip([0.25, 1.0, 0.2]) {
            close(actual, expected, 8e-13 * expected);
        }
        assert_eq!(pixel.optical_depth(), f64::INFINITY);
        assert_eq!(pixel.opacity(), 1.0);
        assert_eq!(pixel.transmittance(), 0.0);
        assert_eq!(pixel.source_error_rgb(), [0.0; 3]);
        assert_eq!(pixel.composite([1e308; 3]).unwrap(), pixel.emission_rgb());
    }
}

#[test]
fn g0_overflowed_varying_profiles_keep_their_quadrature_admission() {
    let constant_source = transfer(&[(0.0, 0.0, [3.0; 3]), (1.0, 1e308, [3.0; 3])]);
    let opaque =
        preintegrate_transfer_segment(&constant_source, 0.0, 1.0, 4.0, settings()).unwrap();
    assert_eq!(opaque.emission_rgb(), [3.0; 3]);
    assert_eq!(opaque.opacity(), 1.0);

    let tf = transfer(&[(0.0, 0.0, [0.0; 3]), (1.0, 1e308, [1e308; 3])]);
    assert_eq!(
        preintegrate_transfer_segment(&tf, 0.0, 1.0, 4.0, settings()),
        Err(PreintegrationError::UnrepresentableOpticalDepth)
    );
}

#[test]
fn g1_thin_layers_and_constant_extinction_source_gradient() {
    for tau in [1e-20, 1e-7, 0.2, 2.0, 1e4] {
        let tf = transfer(&[(0.0, tau, [1.0, 0.0, 4.0]), (1.0, tau, [0.0, 1.0, 4.0])]);
        let pixel = preintegrate_transfer_segment(&tf, 0.0, 1.0, 1.0, settings()).unwrap();
        let alpha = -(-tau).exp_m1();
        // The total red+green source is constant; each colored source is a
        // different first spatial moment of the exponential distribution.
        close(
            pixel.emission_rgb()[0] + pixel.emission_rgb()[1],
            alpha,
            3e-15 * alpha,
        );
        close(pixel.emission_rgb()[2], 4.0 * alpha, 1e-14 * alpha);
        if tau < 1e-6 {
            close(
                pixel.emission_rgb()[1],
                tau * (0.5 - tau / 3.0),
                tau * 1e-13,
            );
        } else {
            let first_moment = (1.0 - (1.0 + tau) * (-tau).exp()) / tau;
            close(pixel.emission_rgb()[1], first_moment, 2e-14);
        }
        if tau == 1e-20 {
            assert_eq!(pixel.transmittance(), 1.0);
            assert!(pixel.opacity() > 0.0);
        }
    }
}

#[test]
fn g1_transparent_endpoints_do_not_hide_a_narrow_opacity_peak() {
    let width = 2_f64.powi(-20);
    let source = [2.0, 4.0, 8.0];
    let tf = transfer(&[
        (0.0, 0.0, source),
        (0.5 - width, 0.0, source),
        (0.5, 1.0 / width, source),
        (0.5 + width, 0.0, source),
        (1.0, 0.0, source),
    ]);
    for (front, back) in [(0.0_f64, 1.0), (1.0, 0.0), (-1.0, 2.0)] {
        let length = (back - front).abs();
        let pixel = preintegrate_transfer_segment(
            &tf,
            front,
            back,
            length,
            PreintegrationSettings {
                max_quadrature_evaluations: 0,
                ..settings()
            },
        )
        .unwrap();
        // Exact area of the triangular extinction peak is one.
        close(pixel.optical_depth(), 1.0, 2e-15);
        close(pixel.opacity(), -(-1.0_f64).exp_m1(), 2e-15);
    }
}

// Independent dense physical-distance Simpson integration. This oracle does
// not invert optical depth, classify cumulative probabilities or use Gauss
// quadrature. The extinction is a+b*x, so its attenuation integral is known.
fn distance_oracle(a: f64, b: f64, source_front: f64, source_back: f64, length: f64) -> f64 {
    let n = 65_536;
    let h = length / f64::from(n);
    let mut sum = 0.0;
    for i in 0..=n {
        let x = f64::from(i) * h;
        let u = x / length;
        let extinction = a + (b - a) * u;
        let source = source_front + (source_back - source_front) * u;
        let attenuation = (-a * x - 0.5 * (b - a) * x * u).exp();
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * extinction * source * attenuation;
    }
    sum * h / 3.0
}

#[test]
fn g1_varying_extinction_and_color_match_independent_distance_integral() {
    for (front_sigma, back_sigma, length) in [
        (0.0, 5.0, 1.0),
        (5.0, 0.0, 1.0),
        (0.2, 2.3, 3.0),
        (3.0, 1.0, 0.01),
        (1e-8, 2e-8, 1.0),
    ] {
        let source_front = [0.0, 7.0, 1.5];
        let source_back = [2.0, 0.0, 4.0];
        let tf = transfer(&[
            (0.0, front_sigma, source_front),
            (1.0, back_sigma, source_back),
        ]);
        let pixel = preintegrate_transfer_segment(&tf, 0.0, 1.0, length, settings()).unwrap();
        for channel in 0..3 {
            let expected = distance_oracle(
                front_sigma,
                back_sigma,
                source_front[channel],
                source_back[channel],
                length,
            );
            close(
                pixel.emission_rgb()[channel],
                expected,
                2e-10 * expected.abs().max(1e-10),
            );
            assert!(pixel.source_error_rgb()[channel] <= 1e-9);
        }
    }
}

#[test]
fn g1_optically_thick_ramp_resolves_the_foreground_without_small_distance_steps() {
    // sigma(x)=k*x and source(x)=x imply emission sqrt(pi/(2*k)) in
    // the optically thick limit; the omitted rear term is exp(-k/2).
    let k = 1e12;
    let tf = transfer(&[(0.0, 0.0, [0.0; 3]), (1.0, k, [1.0; 3])]);
    let pixel = preintegrate_transfer_segment(
        &tf,
        0.0,
        1.0,
        1.0,
        PreintegrationSettings {
            source_tolerance: 1e-12,
            max_quadrature_evaluations: 10_000,
            ..settings()
        },
    )
    .unwrap();
    let expected = (core::f64::consts::PI / (2.0 * k)).sqrt();
    close(pixel.emission_rgb()[0], expected, 2e-13);
    assert_eq!(pixel.opacity(), 1.0);
    assert_eq!(pixel.transmittance(), 0.0);
}

#[test]
fn g3_reversing_endpoints_changes_color_but_preserves_opacity() {
    let tf = transfer(&[(0.0, 0.1, [8.0, 0.0, 0.0]), (1.0, 4.0, [0.0, 8.0, 0.0])]);
    let forward = preintegrate_transfer_segment(&tf, 0.0, 1.0, 2.0, settings()).unwrap();
    let reverse = preintegrate_transfer_segment(&tf, 1.0, 0.0, 2.0, settings()).unwrap();
    close(forward.optical_depth(), reverse.optical_depth(), 1e-15);
    assert!(forward.emission_rgb()[0] > reverse.emission_rgb()[0]);
    assert!(reverse.emission_rgb()[1] > forward.emission_rgb()[1]);
}

#[test]
fn g3_node_grid_refinement_preserves_physical_thickness_and_source_transport() {
    let tf = transfer(&[
        (0.0, 0.0, [0.5, 3.0, 0.0]),
        (0.3, 2.0, [4.0, 1.0, 2.0]),
        (0.7, 0.1, [1.0, 2.0, 4.0]),
        (1.0, 3.0, [0.0, 5.0, 1.0]),
    ]);
    let expected = preintegrate_transfer_segment(&tf, 1.0, 0.0, 2.0, settings()).unwrap();
    for intervals in [1, 2, 7, 32] {
        let nodes = intervals + 1;
        let values: Vec<_> = (0..nodes)
            .flat_map(|z| [z as f64 / intervals as f64; 4])
            .collect();
        let grid = NodeVolumeGrid {
            dims: [2, 2, nodes],
            data: &values,
            origin: [0.0; 3],
            spacing: [1.0, 1.0, 2.0 / intervals as f64],
        };
        let pixel = render_preintegrated_transfer(
            &grid,
            &tf,
            image_settings(values.len(), [1, 1], intervals as u64),
        )
        .unwrap()[0];
        close(pixel.optical_depth(), expected.optical_depth(), 2e-14);
        for channel in 0..3 {
            close(
                pixel.emission_rgb()[channel],
                expected.emission_rgb()[channel],
                2e-9,
            );
        }
    }
}

#[test]
fn g1_node_reconstruction_interpolates_scalars_before_opacity_classification() {
    let tf = transfer(&[
        (0.0, 0.0, [1.0; 3]),
        (0.5, 2.0, [1.0; 3]),
        (1.0, 0.0, [1.0; 3]),
    ]);
    let values = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
    let mut grid = NodeVolumeGrid {
        dims: [2; 3],
        data: &values,
        origin: [0.0; 3],
        spacing: [2.0; 3],
    };
    let request = image_settings(8, [3, 1], 3);
    let first = render_preintegrated_transfer(&grid, &tf, request).unwrap();
    // At x=1/6,1/2,5/6, scalar classification produces sigma=2/3,2,2/3.
    for (pixel, tau) in first.iter().zip([4.0 / 3.0, 4.0, 4.0 / 3.0]) {
        close(pixel.optical_depth(), tau, 2e-15);
    }
    grid.origin = [1e15; 3];
    assert_eq!(
        first,
        render_preintegrated_transfer(&grid, &tf, request).unwrap()
    );
}

#[test]
fn g0_vacuum_zero_length_clamping_and_bright_background_limits() {
    let source = [2.0; 3];
    let tf = transfer(&[(0.0, 0.0, source), (1.0, 1.0, source)]);
    let vacuum = preintegrate_transfer_segment(&tf, -2.0, -1.0, 4.0, settings()).unwrap();
    assert_eq!(vacuum.emission_rgb(), [0.0; 3]);
    assert_eq!(vacuum.composite([3.0; 3]).unwrap(), [3.0; 3]);
    let zero = preintegrate_transfer_segment(
        &tf,
        0.0,
        1.0,
        0.0,
        PreintegrationSettings {
            max_intervals: 0,
            max_quadrature_evaluations: 0,
            ..settings()
        },
    )
    .unwrap();
    assert_eq!(zero, vacuum);
    let black = transfer(&[(0.0, 750.0, [0.0; 3])]);
    let thick = preintegrate_transfer_segment(&black, 0.0, 0.0, 1.0, settings()).unwrap();
    let bright = thick.composite([1e300; 3]).unwrap();
    let expected = (1e300_f64.ln() - 750.0).exp();
    close(bright[0], expected, expected * 2e-13);
    assert!(bright[0] > 0.0);
    assert_eq!(
        thick.composite([-1.0; 3]),
        Err(DvrError::InvalidBackgroundRadiance)
    );
}

#[test]
fn g0_admission_and_whole_image_work_budgets_refuse_without_partial_output() {
    let tf = transfer(&[(0.0, 0.1, [0.0; 3]), (1.0, 0.3, [1.0; 3])]);
    for endpoints in [
        (f64::NAN, 1.0, 1.0),
        (0.0, f64::INFINITY, 1.0),
        (0.0, 1.0, -1.0),
    ] {
        assert_eq!(
            preintegrate_transfer_segment(&tf, endpoints.0, endpoints.1, endpoints.2, settings()),
            Err(PreintegrationError::InvalidSegment)
        );
    }
    assert_eq!(
        preintegrate_transfer_segment(
            &tf,
            0.0,
            1.0,
            1.0,
            PreintegrationSettings {
                source_tolerance: 0.0,
                ..settings()
            }
        ),
        Err(PreintegrationError::InvalidSourceTolerance)
    );
    assert_eq!(
        preintegrate_transfer_segment(
            &tf,
            0.0,
            1.0,
            1.0,
            PreintegrationSettings {
                max_intervals: 0,
                ..settings()
            }
        ),
        Err(PreintegrationError::TransferIntervalBudgetExceeded)
    );
    assert_eq!(
        preintegrate_transfer_segment(
            &tf,
            0.0,
            1.0,
            1.0,
            PreintegrationSettings {
                max_quadrature_evaluations: 23,
                ..settings()
            }
        ),
        Err(PreintegrationError::QuadratureBudgetExceeded)
    );
    let values = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let grid = NodeVolumeGrid {
        dims: [2; 3],
        data: &values,
        origin: [0.0; 3],
        spacing: [1.0; 3],
    };
    let request = PreintegratedImageSettings {
        integration: PreintegrationSettings {
            source_tolerance: 0.1,
            max_quadrature_evaluations: 96,
            ..settings()
        },
        ..image_settings(8, [2, 2], 4)
    };
    assert!(render_preintegrated_transfer(&grid, &tf, request).is_ok());
    assert_eq!(
        render_preintegrated_transfer(
            &grid,
            &tf,
            PreintegratedImageSettings {
                integration: PreintegrationSettings {
                    max_quadrature_evaluations: 72,
                    ..request.integration
                },
                ..request
            }
        ),
        Err(PreintegrationError::QuadratureBudgetExceeded)
    );
    assert_eq!(
        render_preintegrated_transfer(
            &grid,
            &tf,
            PreintegratedImageSettings {
                integration: PreintegrationSettings {
                    max_intervals: 3,
                    ..request.integration
                },
                ..request
            }
        ),
        Err(PreintegrationError::TransferIntervalBudgetExceeded)
    );
    assert!(matches!(
        render_preintegrated_transfer(
            &grid,
            &tf,
            PreintegratedImageSettings {
                max_segments: 3,
                ..request
            }
        ),
        Err(PreintegrationError::Transfer(
            DvrError::CellVisitBudgetExceeded {
                requested: 4,
                limit: 3
            }
        ))
    ));
    assert!(matches!(
        render_preintegrated_transfer(
            &grid,
            &tf,
            PreintegratedImageSettings {
                max_grid_nodes: 7,
                ..request
            }
        ),
        Err(PreintegrationError::Transfer(
            DvrError::GridCellBudgetExceeded { .. }
        ))
    ));
    let invalid = NodeVolumeGrid {
        dims: [1, 2, 4],
        ..grid
    };
    assert_eq!(
        render_preintegrated_transfer(&invalid, &tf, request),
        Err(PreintegrationError::Transfer(DvrError::InvalidGrid))
    );
}
