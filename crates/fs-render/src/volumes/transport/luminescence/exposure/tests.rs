//! Accepted-time transport, native-camera parity and complete work admission.

use super::*;
use super::super::{render_luminescence_camera, render_scattered_luminescence_camera};
use crate::camera::{Aperture, CameraProjection};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey as ExecStreamKey};

const SOURCE: ContentHash = ContentHash([109; 32]);
const FOUR_PI: f64 = 4.0 * core::f64::consts::PI;
const EXPOSURE: LuminescenceExposureBudget = LuminescenceExposureBudget {
    max_nodes: 64,
    max_observations: 4096,
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                ExecStreamKey { seed: 1, kernel_id: 1, tile: 0, iteration: 0 },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn camera() -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 3.0),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_half_tangent(0.01).unwrap(),
        3.0,
        Aperture::try_circular(0.0).unwrap(),
    ).unwrap()
}

fn grid<'a>(cx: &Cx<'_>, density: &'a [f64], absorption: f64) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1, 1, density.len()], density, [-0.5, -0.5, 0.0], [1.0; 3]),
        MassAbsorptionSpectrum::constant(absorption).unwrap(),
        None,
        AbsorbingGridBudget { max_grid_cells: density.len(), max_steps_per_segment: 64 },
    ).unwrap()
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
        "{actual:e} != {expected:e}",
    );
}

#[test]
fn g1_changing_emission_and_foreground_opacity_converge_to_actual_time_integral() {
    with_cx(|cx, _| {
        let density = [0.0; 2];
        let medium = grid(cx, &density, 0.0);
        let view = camera();
        // Unit-length rear emitter j(t)=exp(-0.4t), nonemitting foreground
        // alpha(t)=1.3t. Each actual observation is exp(-1.7t).
        let exact = -(-1.7_f64).exp_m1() / 1.7;
        let mut previous_error = f64::INFINITY;
        for count in [2_usize, 4, 8, 16] {
            let fields: Vec<_> = (0..count).map(|index| {
                let time = (index as f64 + 0.5) / count as f64;
                LuminescenceGrid::try_new(
                    cx, &medium, SOURCE, [550.0],
                    &[(0, [FOUR_PI * (-0.4 * time).exp()])], 1,
                ).unwrap().with_additional_absorption(
                    cx, ContentHash([110; 32]), &[(1, [1.3 * time])], 1,
                ).unwrap()
            }).collect();
            let nodes: Vec<_> = fields.iter().enumerate().map(|(index, field)| {
                LuminescenceExposureNode {
                    time_s: (index as f64 + 0.5) / count as f64,
                    weight: 1.0 / count as f64,
                    field: LuminescenceExposureField::Direct(field),
                }
            }).collect();
            let image = render_luminescence_exposure(
                cx, &nodes, &view, settings([1, 1], 1), [0.0, 1.0], EXPOSURE,
            ).unwrap();
            let mean = image.mean_radiance_w_m2_sr()[0][0];
            let error = (mean - exact).abs();
            assert!(error < previous_error * 0.3);
            previous_error = error;
            assert_eq!(image.cell_visits(), 2 * count as u64);
            assert_eq!(image.integrated_radiance_j_m2_sr(cx).unwrap()[0][0], mean);
        }
        assert!(previous_error < 8e-4 * exact);
        let incorrectly_averaged_state = -(-0.4_f64).exp_m1() / 0.4 * (-0.65_f64).exp();
        assert!((incorrectly_averaged_state - exact).abs() > 20.0 * previous_error);
    });
}

#[test]
fn g1_constant_direct_and_scattered_exposures_retain_native_samples_and_charge_all_nodes() {
    with_cx(|cx, _| {
        for scattering in [0.0, 0.5] {
            let density = [1.0];
            let host = ScatteringGrid::try_new(
                cx, grid(cx, &density, 0.1),
                MassScatteringSpectrum::constant(scattering).unwrap(), PhaseFunction::Isotropic,
            ).unwrap();
            let field = LuminescenceGrid::try_new(
                cx, host.absorption_grid(), SOURCE, [550.0], &[(0, [FOUR_PI])], 1,
            ).unwrap().with_additional_absorption(
                cx, ContentHash([110; 32]), &[(0, [0.7])], 1,
            ).unwrap();
            let scattered = ScatteringLuminescenceGrid::try_new(cx, &field, &host).unwrap();
            let history = LuminescenceScatteringBudget {
                max_scattering_events: 64,
                max_cell_visits: 4096,
            };
            let nodes = [(0.0, 0.2), (0.5, 0.3), (1.0, 0.5)].map(|(time_s, weight)| {
                LuminescenceExposureNode {
                    time_s, weight,
                    field: LuminescenceExposureField::Scattered { field: &scattered, budget: history },
                }
            });
            let view = camera();
            let base = settings([2, 1], 2);
            let image = render_luminescence_exposure(
                cx, &nodes, &view, base, [0.0, 2.0], EXPOSURE,
            ).unwrap();
            let frozen = render_scattered_luminescence_camera(
                cx, &scattered, &view, base, history,
            ).unwrap();
            assert_eq!(image.mean_radiance_w_m2_sr(), frozen.radiance_w_m2_sr());
            assert_eq!(image.camera(), &view);
            assert_eq!(image.settings(), base);
            assert_eq!(image.source_identity(), SOURCE);
            assert_eq!(image.wavelengths_nm(), &[550.0]);
            assert_eq!(image.interval_s(), [0.0, 2.0]);
            for (integrated, mean) in image.integrated_radiance_j_m2_sr(cx).unwrap().iter()
                .zip(image.mean_radiance_w_m2_sr()) {
                close(integrated[0], 2.0 * mean[0]);
            }
            assert_eq!(image.approximate_mean_xyz(cx).unwrap(), frozen.approximate_xyz(cx).unwrap());
            assert!(image.cell_visits() >= 12);
            assert_eq!(
                render_luminescence_exposure(cx, &nodes, &view, base, [0.0, 2.0],
                    LuminescenceExposureBudget { max_observations: 11, ..EXPOSURE }),
                Err(LuminescenceExposureError::Budget),
            );
            if scattering == 0.0 {
                assert_eq!(image.cell_visits(), 12);
                assert_eq!(
                    render_luminescence_exposure(cx, &nodes, &view,
                        LuminescenceCameraSettings { max_cell_visits: 11, ..base },
                        [0.0, 2.0], EXPOSURE),
                    Err(VolumeTransportError::Budget.into()),
                );
                let direct_nodes = nodes.map(|node| LuminescenceExposureNode {
                    field: LuminescenceExposureField::Direct(&field), ..node
                });
                let direct = render_luminescence_exposure(
                    cx, &direct_nodes, &view, base, [0.0, 2.0], EXPOSURE,
                ).unwrap();
                let frozen = render_luminescence_camera(cx, &field, &view, base).unwrap();
                for (mean, pixel) in direct.mean_radiance_w_m2_sr().iter().zip(frozen.pixels()) {
                    assert_eq!(mean, pixel.emitted_radiance_w_m2_sr());
                }
                assert_eq!(direct.mean_radiance_w_m2_sr(), image.mean_radiance_w_m2_sr());
            }
        }
    });
}

#[test]
fn g1_a_short_bright_pulse_survives_both_weight_orders_and_large_finite_means() {
    with_cx(|cx, _| {
        let density = [0.0];
        let medium = grid(cx, &density, 0.0);
        let pulse = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [550.0], &[(0, [FOUR_PI * 1e18])], 1,
        ).unwrap();
        let dark = LuminescenceGrid::try_new(cx, &medium, SOURCE, [550.0], &[], 0).unwrap();
        for ordered in [[(&pulse, 1e-18), (&dark, 1.0)], [(&dark, 1.0), (&pulse, 1e-18)]] {
            let nodes = ordered.map(|(field, weight)| LuminescenceExposureNode {
                time_s: 0.0, weight, field: LuminescenceExposureField::Direct(field),
            });
            let nodes = [nodes[0], LuminescenceExposureNode { time_s: 1.0, ..nodes[1] }];
            let image = render_luminescence_exposure(
                cx, &nodes, &camera(), settings([1, 1], 1), [0.0, 1.0], EXPOSURE,
            ).unwrap();
            close(image.mean_radiance_w_m2_sr()[0][0], 1.0);
        }
        let bright = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [550.0], &[(0, [1e308])], 1,
        ).unwrap();
        let nodes = [0.0, 1.0].map(|time_s| LuminescenceExposureNode {
            time_s, weight: 0.5, field: LuminescenceExposureField::Direct(&bright),
        });
        let image = render_luminescence_exposure(
            cx, &nodes, &camera(), settings([1, 1], 1), [0.0, 100.0], EXPOSURE,
        ).unwrap();
        close(image.mean_radiance_w_m2_sr()[0][0], 1e308 / FOUR_PI);
        assert_eq!(image.integrated_radiance_j_m2_sr(cx), Err(VolumeTransportError::Unrepresentable.into()));
        let dim = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [550.0], &[(0, [0.25 * FOUR_PI])], 1,
        ).unwrap();
        let nodes = [LuminescenceExposureNode {
            time_s: 0.0, weight: 1.0, field: LuminescenceExposureField::Direct(&dim),
        }];
        let image = render_luminescence_exposure(
            cx, &nodes, &camera(), settings([1, 1], 1), [0.0, f64::from_bits(1)], EXPOSURE,
        ).unwrap();
        close(image.mean_radiance_w_m2_sr()[0][0], 0.25);
        assert_eq!(image.integrated_radiance_j_m2_sr(cx), Err(VolumeTransportError::Unrepresentable.into()));
    });
}

#[test]
fn g0_correspondence_time_quadrature_work_and_cancellation_refuse_complete_exposures() {
    with_cx(|cx, gate| {
        let density = [0.0];
        let medium = grid(cx, &density, 0.0);
        let field = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [550.0], &[(0, [FOUR_PI])], 1,
        ).unwrap();
        let other_source = LuminescenceGrid::try_new(
            cx, &medium, ContentHash([111; 32]), [550.0], &[(0, [FOUR_PI])], 1,
        ).unwrap();
        let other_spectrum = LuminescenceGrid::try_new(
            cx, &medium, SOURCE, [600.0], &[(0, [FOUR_PI])], 1,
        ).unwrap();
        let nodes = [0.0, 1.0].map(|time_s| LuminescenceExposureNode {
            time_s, weight: 0.5, field: LuminescenceExposureField::Direct(&field),
        });
        let view = camera();
        let base = settings([1, 1], 1);
        let render = |nodes: &[LuminescenceExposureNode<'_, '_, '_, '_, 1>]| {
            render_luminescence_exposure(cx, nodes, &view, base, [0.0, 1.0], EXPOSURE)
        };
        for field in [&other_source, &other_spectrum] {
            let invalid = [nodes[0], LuminescenceExposureNode {
                field: LuminescenceExposureField::Direct(field), ..nodes[1]
            }];
            assert_eq!(render(&invalid), Err(LuminescenceExposureError::Correspondence));
        }
        for time_s in [0.0, -1.0, 2.0, f64::NAN] {
            let invalid = [nodes[0], LuminescenceExposureNode { time_s, ..nodes[1] }];
            assert_eq!(render(&invalid), Err(LuminescenceExposureError::Time));
        }
        for weight in [0.0, -0.5, 0.4, 1.0, f64::INFINITY, f64::NAN] {
            let invalid = [nodes[0], LuminescenceExposureNode { weight, ..nodes[1] }];
            assert_eq!(render(&invalid), Err(LuminescenceExposureError::Quadrature));
        }
        assert_eq!(render(&[]), Err(LuminescenceExposureError::Quadrature));
        assert_eq!(
            render_luminescence_exposure(cx, &nodes, &view, base, [1.0, 1.0], EXPOSURE),
            Err(LuminescenceExposureError::Time),
        );
        assert_eq!(
            render_luminescence_exposure(cx, &nodes, &view, base, [0.0, 1.0],
                LuminescenceExposureBudget { max_nodes: 1, ..EXPOSURE }),
            Err(LuminescenceExposureError::Budget),
        );
        assert_eq!(
            render_luminescence_exposure(cx, &nodes, &view,
                LuminescenceCameraSettings { resolution: [usize::MAX, 2], ..base },
                [0.0, 1.0], EXPOSURE),
            Err(LuminescenceExposureError::Budget),
        );
        let complete = render(&nodes).unwrap();
        let changed = [nodes[0], LuminescenceExposureNode { time_s: 0.75, ..nodes[1] }];
        assert_ne!(complete.field_history_identity(), render(&changed).unwrap().field_history_identity());
        gate.request();
        assert_eq!(render(&nodes), Err(VolumeTransportError::Cancelled.into()));
        assert_eq!(complete.integrated_radiance_j_m2_sr(cx), Err(VolumeTransportError::Cancelled.into()));
        assert_eq!(complete.approximate_mean_xyz(cx), Err(VolumeTransportError::Cancelled.into()));
    });
}
