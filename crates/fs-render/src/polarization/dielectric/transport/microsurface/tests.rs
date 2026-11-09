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
