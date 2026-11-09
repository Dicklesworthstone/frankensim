use super::*;
use crate::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use core::f64::consts::PI;
use fs_exec::{CancelGate, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
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
        );
        f(&cx, &gate);
    });
}
fn metal() -> ConductorOptics {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([39; 32]),
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
fn beam(index: f64, wavelength: f64, power: Stokes) -> PolarizedBeam {
    let medium = DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        wavelength,
        medium,
        power,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-12 * b.abs().max(1.0), "{a} != {b}");
}
const TWO: [[f64; 3]; 3] = [[0.5, 0.8, 0.125], [0.5, 0.1, 0.375], [0.5, 0.2, 0.3]];

#[test]
fn g1_analytic_first_and_internal_connection_keep_native_mueller_prefix() {
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
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let source_frame =
            PolarizationFrame::try_new(normal.scale(-1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let source = Stokes::try_new(1.0, 0.2, -0.4, 0.5).unwrap();
        for optics in [&bare, &coated] {
            for index in [1.0, 1.33] {
                for wavelength in [450.0, 650.0] {
                    let input = beam(
                        index,
                        wavelength,
                        Stokes::try_new(1.0, 0.3, -0.2, 0.6).unwrap(),
                    );
                    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
                    let result = input
                        .evaluate_smith_conductor_connection(
                            cx,
                            normal,
                            optics,
                            ConductorSurface::try_rough(1.0).unwrap(),
                            source_frame,
                            source,
                            SmithWalkSettings::try_new(2).unwrap(),
                            |i| TWO[i as usize],
                        )
                        .unwrap();
                    assert_eq!(result.facet_events(), 2);
                    // At alpha=1 and a normal source, G_escape=1 and D=1/pi.
                    // The second collision's geometric factor is exactly
                    // 1/[2*pi*(1-w.z)]. It acts on the full first-facet state.
                    let first = input.reflect_conductor(normal, optics).unwrap().beam;
                    let first_value =
                        reciprocal_contraction(first, source_frame, source).unwrap() / (4.0 * PI);
                    let step = model.sample(input.ray().dir, 0.0, TWO[0]).unwrap();
                    let SmithEvent::Reflection { facet, .. } = step.event else {
                        panic!("first facet")
                    };
                    let internal = input.reflect_conductor(facet, optics).unwrap().beam;
                    let w = internal.ray().dir;
                    let connect_facet = unit(Vec3::new(-w.x, -w.y, 1.0 - w.z)).unwrap();
                    let connected = internal
                        .reflect_conductor(connect_facet, optics)
                        .unwrap()
                        .beam;
                    let higher = reciprocal_contraction(connected, source_frame, source).unwrap()
                        / (2.0 * PI * (1.0 - w.z));
                    near(result.factor_per_sr(), first_value + higher);
                    assert!(result.factor_per_sr() > first_value);
                    let one = input
                        .evaluate_smith_conductor_connection(
                            cx,
                            normal,
                            optics,
                            ConductorSurface::try_rough(1.0).unwrap(),
                            source_frame,
                            source,
                            SmithWalkSettings::try_new(1).unwrap(),
                            |_| [0.5, 0.0, 0.0],
                        )
                        .unwrap();
                    near(one.factor_per_sr(), first_value);
                }
            }
        }
    });
}
#[test]
fn g3_connection_preserves_source_frame_camera_scale_and_native_fallback() {
    with_cx(|cx, _| {
        let optics = metal();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = beam(1.2, 550.0, Stokes::try_new(3.0, 0.6, -0.9, 1.2).unwrap());
        let source_frame =
            PolarizationFrame::try_new(Vec3::new(-0.2, 0.1, -1.0), Vec3::new(1.0, 0.0, 0.0))
                .unwrap();
        let source = Stokes::try_new(1.0, 0.2, -0.4, 0.3).unwrap();
        let settings = SmithWalkSettings::try_new(16).unwrap();
        let evaluate = |camera: PolarizedBeam, frame, state| {
            camera
                .evaluate_smith_conductor_connection(
                    cx,
                    normal,
                    &optics,
                    ConductorSurface::try_rough(0.25).unwrap(),
                    frame,
                    state,
                    settings,
                    |_| [0.5, 0.0, 0.0],
                )
                .unwrap()
        };
        let expected = evaluate(input, source_frame, source);
        let camera = input.with_reference_axis(Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let rotated =
            PolarizationFrame::try_new(source_frame.direction(), Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let rotated_source = source_frame.express_in(source, rotated);
        near(
            evaluate(camera, rotated, rotated_source).factor_per_sr(),
            expected.factor_per_sr(),
        );
        let normalized = beam(1.2, 550.0, Stokes::try_new(1.0, 0.2, -0.3, 0.4).unwrap());
        near(
            evaluate(normalized, source_frame, source).factor_per_sr(),
            expected.factor_per_sr(),
        );
        let fallback = input
            .evaluate_instanced_smith_brushed_conductor_connection(
                cx,
                &optics,
                &SurfaceDetail::none(0.25),
                RigidTransform::identity(),
                Point3::new(0.2, 0.0, 0.0),
                normal,
                source_frame,
                source,
                settings,
                |_| [0.5, 0.0, 0.0],
            )
            .unwrap();
        assert_eq!(fallback, expected);
    });
}
#[test]
fn g0_connection_limits_cancellation_and_admission_refuse_whole_estimates() {
    with_cx(|cx, gate| {
        let optics = metal();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = beam(1.0, 550.0, Stokes::unpolarized(1.0));
        let frame = input.frame();
        let surface = ConductorSurface::try_rough(1.0).unwrap();
        let limit = SmithWalkSettings::try_new(1).unwrap();
        assert!(matches!(
            input.evaluate_smith_conductor_connection(
                cx,
                normal,
                &optics,
                surface,
                frame,
                Stokes::unpolarized(1.0),
                limit,
                |i| TWO[i as usize]
            ),
            Err(SmithBeamError::FacetLimit {
                reflected_facets: 1
            })
        ));
        assert!(
            input
                .evaluate_smith_conductor_connection(
                    cx,
                    normal,
                    &optics,
                    surface,
                    frame,
                    Stokes::unpolarized(2.0),
                    limit,
                    |_| panic!("invalid source drew")
                )
                .is_err()
        );
        let mut calls = 0;
        let cancelled = input.evaluate_smith_conductor_connection(
            cx,
            normal,
            &optics,
            surface,
            frame,
            Stokes::unpolarized(1.0),
            limit,
            |_| {
                calls += 1;
                gate.request();
                [0.5, 0.0, 0.0]
            },
        );
        assert_eq!(cancelled, Err(SmithBeamError::Cancelled));
        assert_eq!(calls, 1);
    });
}
