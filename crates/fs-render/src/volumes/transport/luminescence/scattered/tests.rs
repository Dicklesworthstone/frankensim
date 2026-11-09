//! G1/G3/G4 checks against independent radiative-transfer identities.

use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, CellPalette, CellPaletteBudget, MassAbsorptionSpectrum,
    luminescence::{
        LuminescenceCameraSettings, render_luminescence_camera,
        render_scattered_luminescence_camera,
    },
    scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringOptics},
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey as ExecStreamKey};
use std::sync::Arc;

const SOURCE: ContentHash = ContentHash([187; 32]);
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
                    seed: 71,
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

fn absorption<'a>(cx: &Cx<'_>, density: &'a [f64], alpha: f64) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(alpha).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 16,
        },
    )
    .unwrap()
}

fn medium<'a>(cx: &Cx<'_>, density: &'a [f64], a: f64, s: f64) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        absorption(cx, density, a),
        MassScatteringSpectrum::constant(s).unwrap(),
        PhaseFunction::Isotropic,
    )
    .unwrap()
}

fn budget() -> LuminescenceScatteringBudget {
    LuminescenceScatteringBudget {
        max_scattering_events: 256,
        max_cell_visits: 2048,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g0_large_path_weight_preserves_representable_source_radiance() {
    // A packet continuation can exceed one. Its exponential can overflow
    // although the weighted emission remains finite; retain the log fallback.
    let actual = segment_emission(1e-200, 1.0, 1.0, 1.0, 1.0, -720.0);
    let expected = (720.0 - 200.0 * 10.0_f64.ln()).exp() * -(-1.0_f64).exp_m1();
    assert!(actual.is_finite());
    close(actual, expected);
}

#[test]
fn g3_zero_scattering_keeps_exact_source_and_vacuum_cell_emission() {
    with_cx(|cx, _| {
        let density = [0.0, 2.0];
        let medium = medium(cx, &density, 0.7, 0.0);
        let field = LuminescenceGrid::try_new(
            cx,
            medium.absorption_grid(),
            SOURCE,
            [510.0, 650.0],
            &[(0, [FOUR_PI, 2.0 * FOUR_PI]), (1, [3.0 * FOUR_PI, FOUR_PI])],
            2,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &medium).unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let native = field
            .transport(cx, origin, direction, f64::INFINITY)
            .unwrap();
        let observed = scattering
            .sample_ray_with_samples(
                cx,
                origin,
                direction,
                LuminescenceScatteringBudget {
                    max_scattering_events: 0,
                    max_cell_visits: 2,
                },
                |_| panic!("zero scattering must consume no variates"),
            )
            .unwrap();
        assert_eq!(observed.radiance_w_m2_sr, native.emitted_radiance_w_m2_sr);
        assert_eq!(observed.cell_visits, 2);
        assert_eq!(observed.scattering_events, 0);
        close(
            observed.radiance_w_m2_sr[0],
            1.0 + 3.0 * -(-1.4_f64).exp_m1() / 1.4,
        );
        let different = self::medium(cx, &[0.0, 3.0], 0.7, 0.0);
        assert!(matches!(
            ScatteringLuminescenceGrid::try_new(cx, &field, &different),
            Err(VolumeTransportError::SourceIdentity)
        ));
    });
}

#[test]
fn g1_one_scattering_vertex_uses_full_source_and_only_one_continuation_weight() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [4_000.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::table(&[(380.0, 0.7), (780.0, 0.2)]).unwrap(),
                Some(&temperature),
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    max_steps_per_segment: 16,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::table(&[(380.0, 1.3), (780.0, 0.8)]).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        assert!(medium.has_thermal_emission());
        let field = LuminescenceGrid::try_new(
            cx,
            medium.absorption_grid(),
            SOURCE,
            [380.0, 780.0],
            &[(0, [1.25 * FOUR_PI, 2.0 * FOUR_PI])],
            1,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &medium).unwrap();
        let origin = Point3::new(0.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let samples = |event| {
            if event == 0 {
                [0.2, (-0.5_f64).exp(), 0.0, 0.0]
            } else {
                [0.2, 0.001, 0.5, 0.5]
            }
        };
        let result = scattering
            .sample_ray_with_samples(cx, origin, direction, budget(), samples)
            .unwrap();
        // The collision is 0.25 m into a 1 m primary ray, then travels back
        // toward x=0. Each lane uses the COMPLETE spectral mixture density,
        // not its own conditional PDF. Unrelated LTE emission is not added.
        let mixture_density = 0.5 * (2.0 * (-0.5_f64).exp() + (-0.25_f64).exp());
        for lane in 0..2 {
            let extinction = [2.0_f64, 1.0][lane];
            let scattering_coefficient = [1.3, 0.8][lane];
            let source = [1.25, 2.0][lane];
            let primary = source * -(-extinction).exp_m1() / extinction;
            let secondary = source * -(-0.25 * extinction).exp_m1() / extinction;
            let continuation =
                scattering_coefficient * (-0.25 * extinction).exp() / mixture_density;
            close(
                result.radiance_w_m2_sr[lane],
                primary + continuation * secondary,
            );
        }
        assert_eq!(
            (
                result.scattering_events,
                result.ray_segments,
                result.cell_visits
            ),
            (1, 2, 2)
        );
        assert_eq!(
            scattering.sample_ray_with_samples(
                cx,
                origin,
                direction,
                LuminescenceScatteringBudget {
                    max_scattering_events: 0,
                    ..budget()
                },
                samples
            ),
            Err(VolumeTransportError::Budget)
        );
        assert_eq!(
            scattering.sample_ray_with_samples(
                cx,
                origin,
                direction,
                LuminescenceScatteringBudget {
                    max_cell_visits: 1,
                    ..budget()
                },
                samples
            ),
            Err(VolumeTransportError::Budget)
        );
    });
}

#[test]
fn g1_scattering_uses_the_actual_collision_cell_phase() {
    with_cx(|cx, _| {
        let density = [1.0; 2];
        let s = MassScatteringSpectrum::constant(1.0).unwrap();
        let medium = ScatteringGrid::try_new_indexed(
            cx,
            absorption(cx, &density, 0.0),
            CellPalette::try_new(
                cx,
                vec![
                    ScatteringOptics {
                        scattering: s,
                        phase: PhaseFunction::Isotropic,
                    },
                    ScatteringOptics {
                        scattering: s,
                        phase: PhaseFunction::HenyeyGreenstein { g: -0.8 },
                    },
                ],
                Arc::from([0_u32, 1]),
                CellPaletteBudget {
                    max_models: 2,
                    max_cells: 2,
                },
            )
            .unwrap(),
        )
        .unwrap();
        let field = LuminescenceGrid::try_new(
            cx,
            medium.absorption_grid(),
            SOURCE,
            [550.0],
            &[(0, [FOUR_PI]), (1, [FOUR_PI])],
            2,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &medium).unwrap();
        let observed = scattering
            .sample_ray_with_samples(
                cx,
                Point3::new(1.5, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                budget(),
                |event| {
                    if event == 0 {
                        [0.5, (-0.1_f64).exp(), 0.5, 0.25]
                    } else {
                        [0.5, 0.001, 0.5, 0.5]
                    }
                },
            )
            .unwrap();
        // Independent HG inverse at u=1/2, rotated about +x. The secondary
        // ray exits y=1 after crossing from original cell 1 to cell 0.
        let g = -0.8_f64;
        let cosine = (1.0 + g * g - (1.0 - g * g).powi(2)) / (2.0 * g);
        let secondary_length = (1.6 / -cosine).min(0.5 / (1.0 - cosine * cosine).sqrt());
        let expected = -(-0.5_f64).exp_m1() - (-secondary_length).exp_m1();
        close(observed.radiance_w_m2_sr[0], expected);
        assert_eq!(observed.scattering_events, 1);
        assert_eq!(observed.cell_visits, 3);
    });
}

#[test]
fn g1_lossless_isotropic_scattering_preserves_total_emitted_power() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = medium(cx, &density, 0.0, 1.0);
        let field = LuminescenceGrid::try_new(
            cx,
            medium.absorption_grid(),
            SOURCE,
            [550.0],
            &[(0, [FOUR_PI])],
            1,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &medium).unwrap();
        let mut stream = StreamKey {
            seed: 1729,
            kernel: 0x666c_7578,
            tile: 0,
        }
        .stream();
        let (mut sum, mut squares) = (0.0, 0.0);
        let count = 24_000_u64;
        for id in 0..count {
            let origin = Point3::new(0.0, stream.next_f64(), stream.next_f64());
            let cosine = stream.next_f64().sqrt();
            let sine = (1.0 - cosine * cosine).sqrt();
            let phi = core::f64::consts::TAU * stream.next_f64();
            let direction = Vec3::new(cosine, sine * phi.cos(), sine * phi.sin());
            let value = scattering
                .sample_ray(cx, origin, direction, budget(), 97, id)
                .unwrap()
                .radiance_w_m2_sr[0];
            sum += value;
            squares += value * value;
        }
        // Every emitted joule leaves the lossless cube. Symmetry gives equal
        // face flux; cosine-weighted boundary rays give Q = 6*pi*mean(L).
        let mean = sum / count as f64;
        let error = ((squares / count as f64 - mean * mean) / (count - 1) as f64).sqrt();
        assert!(
            (mean - 2.0 / 3.0).abs() < 6.0 * error + 0.002,
            "boundary mean {mean}, reference 2/3, standard error {error}"
        );
        assert!(error < 0.01);
        let a = scattering
            .sample_ray(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                budget(),
                97,
                31,
            )
            .unwrap();
        let b = scattering
            .sample_ray(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                budget(),
                97,
                31,
            )
            .unwrap();
        assert_eq!(a, b);
    });
}

#[test]
fn g3_camera_keeps_zero_scattering_rays_and_whole_image_work_limit() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = medium(cx, &density, 0.5, 0.0);
        let field = LuminescenceGrid::try_new(
            cx,
            medium.absorption_grid(),
            SOURCE,
            [510.0, 650.0],
            &[(0, [FOUR_PI, 2.0 * FOUR_PI])],
            1,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &medium).unwrap();
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.5, 0.5, -1.0),
            Point3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.1).unwrap(),
            1.5,
            Aperture::try_circular(0.02).unwrap(),
        )
        .unwrap();
        let settings = LuminescenceCameraSettings {
            resolution: [2, 2],
            samples_per_pixel: 7,
            seed: 18,
            max_pixels: 4,
            max_samples: 28,
            max_cell_visits: 28,
        };
        let direct = render_luminescence_camera(cx, &field, &camera, settings).unwrap();
        let observed =
            render_scattered_luminescence_camera(cx, &scattering, &camera, settings, budget())
                .unwrap();
        assert_eq!(observed.resolution(), [2, 2]);
        assert_eq!(observed.source_identity(), SOURCE);
        assert_eq!(observed.field_identity(), scattering.content_identity());
        for (a, b) in direct.pixels().iter().zip(observed.radiance_w_m2_sr()) {
            assert_eq!(a.emitted_radiance_w_m2_sr(), b);
        }
        assert_eq!(
            direct.approximate_emission_xyz(cx).unwrap(),
            observed.approximate_xyz(cx).unwrap()
        );
        assert!(matches!(
            render_scattered_luminescence_camera(
                cx,
                &scattering,
                &camera,
                LuminescenceCameraSettings {
                    max_cell_visits: 27,
                    ..settings
                },
                budget()
            ),
            Err(super::super::LuminescenceCameraError::Transport(
                VolumeTransportError::Budget
            ))
        ));
    });
}

#[test]
fn g4_cancellation_and_bad_variates_publish_no_partial_radiance() {
    with_cx(|cx, gate| {
        let density = [1.0];
        let medium = medium(cx, &density, 0.2, 0.8);
        let field = LuminescenceGrid::try_new(
            cx,
            medium.absorption_grid(),
            SOURCE,
            [550.0],
            &[(0, [FOUR_PI])],
            1,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &medium).unwrap();
        let origin = Point3::new(0.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        assert_eq!(
            scattering.sample_ray_with_samples(cx, origin, direction, budget(), |_| [
                0.5,
                0.5,
                f64::NAN,
                0.5
            ]),
            Err(VolumeTransportError::Sample)
        );
        let cancelled =
            scattering.sample_ray_with_samples(cx, origin, direction, budget(), |event| {
                if event == 1 {
                    gate.request();
                }
                [0.5, (-0.1_f64).exp(), 0.0, 0.0]
            });
        assert!(matches!(cancelled, Err(VolumeTransportError::Cancelled)));
    });
}

#[test]
fn g1_active_molecular_opacity_attenuates_sources_and_background_in_original_cells() {
    with_cx(|cx, _| {
        let density = [0.0, 2.0];
        let host = absorption(cx, &density, 0.7);
        let field = || {
            LuminescenceGrid::try_new(
                cx,
                &host,
                SOURCE,
                [510.0, 650.0],
                &[(0, [FOUR_PI, 2.0 * FOUR_PI]), (1, [3.0 * FOUR_PI, FOUR_PI])],
                2,
            )
            .unwrap()
        };
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let passive = field()
            .transport(cx, origin, direction, f64::INFINITY)
            .unwrap();
        let active = field()
            .with_additional_absorption(
                cx,
                ContentHash([188; 32]),
                &[(1, [1.0, 0.0]), (0, [0.5, 0.0])],
                2,
            )
            .unwrap();
        let observed = active
            .transport(cx, origin, direction, f64::INFINITY)
            .unwrap();
        // Exact two-layer formal solution: j0=1, alpha0=.5 followed by
        // j1=3, alpha1=2*.7+1. The second spectral lane has no extra opacity.
        let expected =
            -(-0.5_f64).exp_m1() / 0.5 + (-0.5_f64).exp() * 3.0 * -(-2.4_f64).exp_m1() / 2.4;
        close(observed.emitted_radiance_w_m2_sr[0], expected);
        close(observed.transmittance[0], (-2.9_f64).exp());
        assert!(observed.emitted_radiance_w_m2_sr[0] < passive.emitted_radiance_w_m2_sr[0]);
        assert_eq!(
            observed.emitted_radiance_w_m2_sr[1],
            passive.emitted_radiance_w_m2_sr[1]
        );
        assert_eq!(observed.transmittance[1], passive.transmittance[1]);
        assert_eq!(
            active.content_identity(),
            field()
                .with_additional_absorption(
                    cx,
                    ContentHash([188; 32]),
                    &[(0, [0.5, 0.0]), (1, [1.0, 0.0])],
                    2,
                )
                .unwrap()
                .content_identity(),
        );
        assert_ne!(active.content_identity(), field().content_identity());
    });
}

#[test]
fn g0_active_opacity_preserves_representable_opaque_cell_radiance() {
    with_cx(|cx, _| {
        let density = [1e155];
        let host = absorption(cx, &density, 1e154);
        let field =
            LuminescenceGrid::try_new(cx, &host, SOURCE, [550.0], &[(0, [FOUR_PI * 1e300])], 1)
                .unwrap()
                .with_additional_absorption(cx, ContentHash([188; 32]), &[(0, [1e308])], 1)
                .unwrap();
        let observed = field
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        assert_eq!(observed.transmittance, [0.0]);
        // j/alpha is representable although alpha=1.1e309 and tau overflow.
        close(observed.emitted_radiance_w_m2_sr[0], 1e-9 / 1.1);
    });
}

#[test]
fn g4_active_opacity_refuses_invalid_bindings_and_mismatched_scattering() {
    with_cx(|cx, gate| {
        let density = [1.0, 1.0];
        let host = medium(cx, &density, 0.0, 0.3);
        let field = || {
            LuminescenceGrid::try_new(
                cx,
                host.absorption_grid(),
                SOURCE,
                [550.0],
                &[(1, [FOUR_PI])],
                1,
            )
            .unwrap()
        };
        let identity = ContentHash([188; 32]);
        assert!(matches!(
            field().with_additional_absorption(cx, identity, &[(0, [1.0])], 0),
            Err(VolumeTransportError::Budget),
        ));
        assert!(matches!(
            field().with_additional_absorption(cx, identity, &[(0, [1.0]), (0, [2.0])], 2),
            Err(VolumeTransportError::Palette),
        ));
        assert!(matches!(
            field().with_additional_absorption(cx, identity, &[(2, [1.0])], 1),
            Err(VolumeTransportError::Palette),
        ));
        assert!(matches!(
            field().with_additional_absorption(cx, identity, &[(0, [f64::NAN])], 1),
            Err(VolumeTransportError::AbsorptionSpectrum),
        ));
        // An absorbing foreground cell has no source; it must still attenuate.
        let active = field()
            .with_additional_absorption(cx, identity, &[(0, [1.0])], 1)
            .unwrap();
        let observed = active
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        close(observed.emitted_radiance_w_m2_sr[0], (-1.0_f64).exp());
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &active, &host).unwrap();
        let passive = field();
        assert_ne!(
            scattering.content_identity(),
            ScatteringLuminescenceGrid::try_new(cx, &passive, &host)
                .unwrap()
                .content_identity(),
        );
        let different = medium(cx, &[1.0, 2.0], 0.0, 0.3);
        assert!(matches!(
            ScatteringLuminescenceGrid::try_new(cx, &active, &different),
            Err(VolumeTransportError::SourceIdentity),
        ));
        let unbound = field();
        gate.request();
        assert!(matches!(
            unbound.with_additional_absorption(cx, identity, &[(0, [1.0])], 1),
            Err(VolumeTransportError::Cancelled),
        ));
    });
}

#[test]
fn g1_active_absorption_uses_total_linear_extinction_in_spectral_scattering() {
    with_cx(|cx, _| {
        let density = [2.0];
        let host = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::table(&[(380.0, 0.2), (780.0, 0.7)]).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    max_steps_per_segment: 16,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::table(&[(380.0, 0.3), (780.0, 0.1)]).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let field = LuminescenceGrid::try_new(
            cx,
            host.absorption_grid(),
            SOURCE,
            [380.0, 780.0],
            &[(0, [1.25 * FOUR_PI, 2.0 * FOUR_PI])],
            1,
        )
        .unwrap()
        .with_additional_absorption(cx, ContentHash([188; 32]), &[(0, [0.4, 0.0])], 1)
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &host).unwrap();
        let result = scattering
            .sample_ray_with_samples(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                budget(),
                |event| {
                    if event == 0 {
                        [0.2, (-0.35_f64).exp(), 0.0, 0.0]
                    } else {
                        [0.2, 0.001, 0.5, 0.5]
                    }
                },
            )
            .unwrap();
        // beta=rho*(kappa_a+kappa_s)+alpha; sigma=rho*kappa_s.
        // Hero lane 0 places the vertex at x=.25. The full primary source
        // integral and the reverse secondary ray use the same beta, while
        // continuation uses the mixture over BOTH physical collision rates.
        let beta = [1.4_f64, 1.6];
        let sigma = [0.6, 0.2];
        let proposal =
            0.5 * (beta[0] * (-0.25 * beta[0]).exp() + beta[1] * (-0.25 * beta[1]).exp());
        for lane in 0..2 {
            let source = [1.25, 2.0][lane];
            let primary = source * -(-beta[lane]).exp_m1() / beta[lane];
            let secondary = source * -(-0.25 * beta[lane]).exp_m1() / beta[lane];
            let continuation = sigma[lane] * (-0.25 * beta[lane]).exp() / proposal;
            close(
                result.radiance_w_m2_sr[lane],
                primary + continuation * secondary,
            );
        }
        assert_eq!(
            (
                result.scattering_events,
                result.ray_segments,
                result.cell_visits
            ),
            (1, 2, 2)
        );
    });
}

#[test]
fn g1_active_absorber_in_zero_density_host_preserves_the_complete_source_ray() {
    with_cx(|cx, _| {
        let density = [0.0, 1.0];
        let host = medium(cx, &density, 0.0, 1.0);
        let field = LuminescenceGrid::try_new(
            cx,
            host.absorption_grid(),
            SOURCE,
            [550.0],
            &[(1, [FOUR_PI])],
            1,
        )
        .unwrap()
        .with_additional_absorption(cx, ContentHash([188; 32]), &[(0, [0.7])], 1)
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &host).unwrap();
        let allowance = LuminescenceScatteringBudget {
            max_scattering_events: 0,
            max_cell_visits: 2,
        };
        let origin = Point3::new(0.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let result = scattering
            .sample_ray_with_samples(cx, origin, direction, allowance, |_| {
                [0.5, (-0.35_f64).exp(), 0.5, 0.5]
            })
            .unwrap();
        // The sampled event at x=.5 is pure absorption: no phase event follows.
        // The emitter lies beyond it, so stopping the full source ray there
        // would incorrectly erase this independently integrated contribution.
        close(
            result.radiance_w_m2_sr[0],
            (-0.7_f64).exp() * -(-1.0_f64).exp_m1(),
        );
        assert_eq!(
            (
                result.scattering_events,
                result.ray_segments,
                result.cell_visits
            ),
            (0, 1, 2)
        );
        assert_eq!(
            scattering.sample_ray_with_samples(
                cx,
                origin,
                direction,
                LuminescenceScatteringBudget {
                    max_cell_visits: 1,
                    ..allowance
                },
                |_| [0.5, (-0.35_f64).exp(), 0.5, 0.5],
            ),
            Err(VolumeTransportError::Budget)
        );
    });
}

#[test]
fn g3_active_opacity_with_zero_scattering_keeps_exact_source_without_variates() {
    with_cx(|cx, _| {
        let density = [0.0, 2.0];
        let host = medium(cx, &density, 0.7, 0.0);
        let field = LuminescenceGrid::try_new(
            cx,
            host.absorption_grid(),
            SOURCE,
            [510.0, 650.0],
            &[(0, [FOUR_PI, 2.0 * FOUR_PI]), (1, [3.0 * FOUR_PI, FOUR_PI])],
            2,
        )
        .unwrap()
        .with_additional_absorption(
            cx,
            ContentHash([188; 32]),
            &[(0, [0.5, 0.0]), (1, [1.0, 0.0])],
            2,
        )
        .unwrap();
        let scattering = ScatteringLuminescenceGrid::try_new(cx, &field, &host).unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let native = field
            .transport(cx, origin, direction, f64::INFINITY)
            .unwrap();
        let observed = scattering
            .sample_ray_with_samples(
                cx,
                origin,
                direction,
                LuminescenceScatteringBudget {
                    max_scattering_events: 0,
                    max_cell_visits: 2,
                },
                |_| panic!("zero scattering must consume no variates"),
            )
            .unwrap();
        assert_eq!(observed.radiance_w_m2_sr, native.emitted_radiance_w_m2_sr);
        close(
            observed.radiance_w_m2_sr[0],
            -(-0.5_f64).exp_m1() / 0.5 + (-0.5_f64).exp() * 3.0 * -(-2.4_f64).exp_m1() / 2.4,
        );
        assert_eq!(
            (
                observed.scattering_events,
                observed.ray_segments,
                observed.cell_visits
            ),
            (0, 1, 2)
        );
    });
}
