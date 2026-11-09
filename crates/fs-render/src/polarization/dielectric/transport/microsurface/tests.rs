use super::*;
use crate::conductor::{
    CoatingIorSample, CoherentAbsorbingCoating, CoherentDielectricCoating, ConductorDataStatus,
    ConductorIorSample, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::Stokes;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 41,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                asupersync::types::Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
fn glass(index: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}
fn metal() -> ConductorOptics {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([79; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        source,
    )
    .unwrap()
}
fn beam(index: f64, wavelength: f64) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        wavelength,
        glass(index),
        Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap(),
    )
    .unwrap()
}
const DRAWS: [[f64; 3]; 3] = [[0.5, 0.8, 0.125], [0.5, 0.1, 0.375], [0.5, 0.2, 0.3]];
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-12, "{a} != {b}");
}

#[test]
fn successive_native_mueller_responses_retain_phase_coatings_and_medium() {
    with_cx(|cx, _| {
        let bare = metal();
        let film = bare.with_coating(
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.45).unwrap(),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        let absorbing = bare.with_absorbing_coating(
            CoherentAbsorbingCoating::try_new(
                core::array::from_fn(|i| {
                    CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.7, 0.05).unwrap()
                }),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        for optics in [&bare, &film, &absorbing] {
            for wavelength in [400.0, 550.0, 750.0] {
                for index in [1.0, 1.33] {
                    let input = beam(index, wavelength);
                    let result = input
                        .sample_smith_conductor(
                            cx,
                            Vec3::new(0.0, 0.0, 1.0),
                            optics,
                            ConductorSurface::try_rough(1.0).unwrap(),
                            SmithWalkSettings::try_new(2).unwrap(),
                            |event| DRAWS[event as usize],
                        )
                        .unwrap();
                    assert_eq!(result.facet_events(), 2);
                    assert!(!result.is_delta());
                    assert_eq!(result.beam().medium(), input.medium());
                    assert_eq!(result.beam().ray().origin, input.ray().origin);
                    assert_eq!(result.beam().wavelength_nm(), wavelength);
                    assert!(
                        result.beam().power().i > 0.0 && result.beam().power().i <= input.power().i
                    );
                    assert!(result.beam().power().degree() <= 1.0 + 1e-13);
                    // Replay the two actual facet normals through the original
                    // smooth adapter. Replacing either step with scalar Fresnel,
                    // depolarization, or an extra G factor changes this result.
                    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
                    let mut expected = input;
                    let mut depth = 0.0;
                    let mut log_density = 0.0;
                    for samples in DRAWS {
                        let step = model
                            .sample(expected.frame().direction(), depth, samples)
                            .unwrap();
                        log_density += step.log_density;
                        match step.event {
                            SmithEvent::Reflection {
                                depth: next, facet, ..
                            } => {
                                expected = expected.reflect_conductor(facet, optics).unwrap().beam;
                                depth = next;
                            }
                            SmithEvent::Escape => break,
                        }
                    }
                    let actual = result
                        .beam()
                        .with_reference_axis(expected.frame().horizontal())
                        .unwrap();
                    for (a, b) in [
                        actual.power().i,
                        actual.power().q,
                        actual.power().u,
                        actual.power().v,
                    ]
                    .into_iter()
                    .zip([
                        expected.power().i,
                        expected.power().q,
                        expected.power().u,
                        expected.power().v,
                    ]) {
                        close(a, b);
                    }
                    close(actual.frame().direction().x, expected.frame().direction().x);
                    close(actual.frame().direction().y, expected.frame().direction().y);
                    close(actual.frame().direction().z, expected.frame().direction().z);
                    close(result.log_path_density(), log_density);
                }
            }
        }
    });
}

#[test]
fn mirror_limit_consumes_no_samples_and_facet_limit_never_returns_black() {
    with_cx(|cx, _| {
        let input = beam(1.33, 550.0);
        let optics = metal();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let result = input
            .sample_smith_conductor(
                cx,
                normal,
                &optics,
                ConductorSurface::SMOOTH,
                SmithWalkSettings::try_new(1).unwrap(),
                |_| panic!("delta mirror drew a sample"),
            )
            .unwrap();
        assert!(result.is_delta());
        assert_eq!(result.log_path_density(), 0.0);
        assert_eq!(
            result.beam(),
            input.reflect_conductor(normal, &optics).unwrap().beam
        );
        let original = input;
        assert_eq!(
            input.sample_smith_conductor(
                cx,
                normal,
                &optics,
                ConductorSurface::try_rough(1.0).unwrap(),
                SmithWalkSettings::try_new(1).unwrap(),
                |event| DRAWS[event as usize]
            ),
            Err(SmithBeamError::FacetLimit {
                reflected_facets: 1
            })
        );
        assert_eq!(input, original);
        assert_eq!(
            SmithWalkSettings::try_new(0),
            Err(SmithBeamError::InvalidBudget)
        );
        assert_eq!(
            SmithWalkSettings::try_new(MAX_SMITH_FACET_EVENTS + 1),
            Err(SmithBeamError::InvalidBudget)
        );
    });
}

#[test]
fn cancellation_is_observed_before_next_random_draw() {
    with_cx(|cx, gate| {
        let mut calls = 0;
        let result = beam(1.0, 550.0).sample_smith_conductor(
            cx,
            Vec3::new(0.0, 0.0, 1.0),
            &metal(),
            ConductorSurface::try_rough(1.0).unwrap(),
            SmithWalkSettings::try_new(8).unwrap(),
            |event| {
                calls += 1;
                gate.request();
                DRAWS[event as usize]
            },
        );
        assert_eq!(result, Err(SmithBeamError::Cancelled));
        assert_eq!(calls, 1);
    });
}

fn anisotropic_draw(event: u32) -> [f64; 3] {
    match event {
        0 => [0.5, 0.95, 0.125],
        1 => [0.5, 0.1, 0.375],
        _ => [0.01, 0.05, 0.25],
    }
}
fn same_power(a: Stokes, b: Stokes) {
    for (a, b) in [a.i, a.q, a.u, a.v].into_iter().zip([b.i, b.q, b.u, b.v]) {
        close(a, b);
    }
}

#[test]
fn g1_anisotropic_walk_replays_every_native_facet_in_the_authored_frame() {
    use crate::surface_detail::{BrushingPattern, ShadingFrame, SurfaceDetail};
    with_cx(|cx, _| {
        let bare = metal();
        let coated = bare.with_coating(
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.45).unwrap(),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        let absorbing = bare.with_absorbing_coating(
            CoherentAbsorbingCoating::try_new(
                core::array::from_fn(|i| {
                    CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.7, 0.05).unwrap()
                }),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        let detail = SurfaceDetail::try_new(0.25, 4.0, BrushingPattern::Radial, 0.0, None).unwrap();
        let mut multiple = false;
        for mark in [false, true] {
            let frame = ShadingFrame {
                geometric_normal: Vec3::new(0.0, 0.0, 1.0),
                shading_normal: Vec3::new(0.0, 0.0, 1.0),
                tangent: Vec3::new(1.0, 0.0, 0.0),
                bitangent: Vec3::new(0.0, 1.0, 0.0),
                mark_covered: mark,
            };
            let (x, y) = detail.alpha_pair_for_mark(mark);
            let model = SmithGgx::try_new(x, y).unwrap();
            for optics in [&bare, &coated, &absorbing] {
                for index in [1.0, 1.33] {
                    for wavelength in [450.0, 650.0] {
                        let input = beam(index, wavelength);
                        let actual = input
                            .sample_smith_brushed_conductor(
                                cx,
                                optics,
                                &detail,
                                &frame,
                                SmithWalkSettings::try_new(64).unwrap(),
                                anisotropic_draw,
                            )
                            .unwrap()
                            .unwrap();
                        let mut expected = input;
                        let mut depth = 0.0;
                        let mut log_density = 0.0;
                        let mut facets = 0;
                        for event in 0..=64 {
                            let step = model
                                .sample(expected.ray().dir, depth, anisotropic_draw(event))
                                .unwrap();
                            log_density += step.log_density;
                            match step.event {
                                SmithEvent::Escape => break,
                                SmithEvent::Reflection {
                                    depth: next, facet, ..
                                } => {
                                    facets += 1;
                                    expected =
                                        expected.reflect_conductor(facet, optics).unwrap().beam;
                                    depth = next;
                                }
                            }
                        }
                        multiple |= facets > 1;
                        assert_eq!(actual.facet_events(), facets);
                        close(actual.log_path_density(), log_density);
                        assert_eq!(actual.beam().medium(), input.medium());
                        assert_eq!(actual.beam().ray().origin, input.ray().origin);
                        assert_eq!(actual.beam().wavelength_nm(), wavelength);
                        same_power(actual.beam().power(), expected.power());
                        close(actual.beam().ray().dir.x, expected.ray().dir.x);
                        close(actual.beam().ray().dir.y, expected.ray().dir.y);
                        close(actual.beam().ray().dir.z, expected.ray().dir.z);
                        assert!(actual.beam().power().i <= input.power().i);
                    }
                }
            }
        }
        assert!(multiple, "fixture must traverse hidden facet returns");
    });
}

#[test]
fn g3_body_locked_smith_rotation_and_exact_isotropic_fallbacks() {
    use crate::instances::RigidTransform;
    use crate::surface_detail::{BrushingPattern, RadialMark, SurfaceDetail};
    with_cx(|cx, _| {
        let input = beam(1.33, 550.0);
        let optics = metal();
        let settings = SmithWalkSettings::try_new(64).unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let point = Point3::new(0.4, 0.0, 0.0);
        let angle = 0.7_f64;
        let transform = RigidTransform::try_new(
            [0.0, (angle / 2.0).sin(), 0.0, (angle / 2.0).cos()],
            [2.0, -3.0, 4.0],
        )
        .unwrap();
        let rotated_input = PolarizedBeam::try_new(
            transform.transform_point(input.ray().origin),
            PolarizationFrame::try_new(
                transform.transform_vector(input.ray().dir),
                transform.transform_vector(input.frame().horizontal()),
            )
            .unwrap(),
            input.wavelength_nm(),
            input.medium(),
            input.power(),
        )
        .unwrap();
        for pattern in [BrushingPattern::Radial, BrushingPattern::Circular] {
            let detail = SurfaceDetail::try_new(
                0.25,
                4.0,
                pattern,
                0.0,
                Some(RadialMark {
                    azimuth_rad: 0.0,
                    half_width_rad: 0.1,
                }),
            )
            .unwrap();
            let local = input
                .sample_instanced_smith_brushed_conductor(
                    cx,
                    &optics,
                    &detail,
                    RigidTransform::identity(),
                    point,
                    normal,
                    settings,
                    anisotropic_draw,
                )
                .unwrap()
                .unwrap();
            let world = rotated_input
                .sample_instanced_smith_brushed_conductor(
                    cx,
                    &optics,
                    &detail,
                    transform,
                    point,
                    transform.transform_vector(normal),
                    settings,
                    anisotropic_draw,
                )
                .unwrap()
                .unwrap();
            assert_eq!(world.facet_events(), local.facet_events());
            close(world.log_path_density(), local.log_path_density());
            let expected_direction = transform.transform_vector(local.beam().ray().dir);
            close(world.beam().ray().dir.x, expected_direction.x);
            close(world.beam().ray().dir.y, expected_direction.y);
            close(world.beam().ray().dir.z, expected_direction.z);
            let aligned = world
                .beam()
                .with_reference_axis(transform.transform_vector(local.beam().frame().horizontal()))
                .unwrap();
            same_power(aligned.power(), local.beam().power());
            assert_eq!(world.beam().ray().origin, rotated_input.ray().origin);

            let pole = input
                .sample_instanced_smith_brushed_conductor(
                    cx,
                    &optics,
                    &detail,
                    RigidTransform::identity(),
                    Point3::new(0.0, 0.0, 1.0),
                    normal,
                    settings,
                    anisotropic_draw,
                )
                .unwrap()
                .unwrap();
            let zero_detail = input
                .sample_instanced_smith_brushed_conductor(
                    cx,
                    &optics,
                    &SurfaceDetail::none(0.25),
                    RigidTransform::identity(),
                    point,
                    normal,
                    settings,
                    anisotropic_draw,
                )
                .unwrap()
                .unwrap();
            let isotropic = input
                .sample_smith_conductor(
                    cx,
                    normal,
                    &optics,
                    ConductorSurface::try_rough(0.25).unwrap(),
                    settings,
                    anisotropic_draw,
                )
                .unwrap();
            assert_eq!(pole, isotropic);
            assert_eq!(zero_detail, isotropic);
        }
    });
}

#[test]
fn g1_g4_tilted_smith_applies_one_endpoint_adjoint_and_never_retries_clipping() {
    use crate::surface_detail::{BrushingPattern, SurfaceDetail, shading_frame};
    with_cx(|cx, gate| {
        let input = beam(1.0, 550.0);
        let optics = metal();
        let detail = SurfaceDetail::try_new(0.25, 4.0, BrushingPattern::Radial, 0.2, None).unwrap();
        let frame = shading_frame(
            &detail,
            Point3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        )
        .unwrap();
        let settings = SmithWalkSettings::try_new(64).unwrap();
        let forward = input
            .sample_smith_brushed_conductor(
                cx,
                &optics,
                &detail,
                &frame,
                settings,
                anisotropic_draw,
            )
            .unwrap()
            .unwrap();
        let radiance = input
            .sample_smith_brushed_conductor_for_transport(
                cx,
                &optics,
                &detail,
                &frame,
                settings,
                anisotropic_draw,
                true,
            )
            .unwrap()
            .unwrap();
        assert!(forward.facet_events() > 1);
        assert_eq!(forward.facet_events(), radiance.facet_events());
        assert_eq!(forward.log_path_density(), radiance.log_path_density());
        assert_eq!(forward.beam().ray().dir, radiance.beam().ray().dir);
        let wo = input.ray().dir.scale(-1.0);
        let wi = forward.beam().ray().dir;
        let adjoint = (frame.shading_normal.dot(wo) / frame.geometric_normal.dot(wo))
            * (frame.geometric_normal.dot(wi) / frame.shading_normal.dot(wi));
        let p = radiance.beam().power();
        same_power(
            forward.beam().power(),
            Stokes {
                i: p.i * adjoint,
                q: p.q * adjoint,
                u: p.u * adjoint,
                v: p.v * adjoint,
            },
        );
        let mut rejected_escape = false;
        for i in 1..32 {
            let mut draws = 0;
            let result = input
                .sample_smith_brushed_conductor(cx, &optics, &detail, &frame, settings, |event| {
                    draws += 1;
                    if event == 0 {
                        [0.99, f64::from(i) / 32.0, 0.25]
                    } else {
                        [0.000_001, 0.05, 0.25]
                    }
                })
                .unwrap();
            if result.is_none() {
                assert_eq!(
                    draws, 2,
                    "geometric clipping must not retry the complete path"
                );
                rejected_escape = true;
            }
        }
        assert!(
            rejected_escape,
            "fixture must escape below the geometric horizon"
        );
        assert!(
            input
                .sample_smith_brushed_conductor(
                    cx,
                    &optics,
                    &SurfaceDetail::none(f64::NAN),
                    &frame,
                    settings,
                    |_| panic!("invalid detail drew"),
                )
                .is_err()
        );
        let mut malformed = frame;
        malformed.bitangent = malformed.tangent;
        assert!(
            input
                .sample_smith_brushed_conductor(
                    cx,
                    &optics,
                    &detail,
                    &malformed,
                    settings,
                    |_| panic!("invalid frame drew"),
                )
                .is_err()
        );
        assert_eq!(
            input.sample_smith_brushed_conductor(
                cx,
                &optics,
                &detail,
                &frame,
                SmithWalkSettings::try_new(1).unwrap(),
                anisotropic_draw,
            ),
            Err(SmithBeamError::FacetLimit {
                reflected_facets: 1
            }),
        );
        let mut calls = 0;
        let cancelled =
            input.sample_smith_brushed_conductor(cx, &optics, &detail, &frame, settings, |event| {
                calls += 1;
                gate.request();
                anisotropic_draw(event)
            });
        assert_eq!(cancelled, Err(SmithBeamError::Cancelled));
        assert_eq!(calls, 1);
    });
}
