use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorIorSample, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use fs_blake3::ContentHash;
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
                    seed: 43,
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

fn metal(eta: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, 3.0).unwrap()
        }),
        ConductorSource::try_new(ContentHash([83; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap()
}

fn camera(direction: Vec3, power: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        550.0,
        DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        ),
        power,
    )
    .unwrap()
}

fn emitter(
    optics: ConductorOptics,
    film_temperature: f64,
    substrate_temperature: f64,
) -> SmithCoatingThermalEmission {
    let source = if let Some(film) = optics.absorbing_coating() {
        CoatingThermalEmission::try_opaque(
            &[CoatingLayer::Absorbing(&film)],
            &[film_temperature],
            optics,
            substrate_temperature,
            ContentHash([84; 32]),
        )
    } else if let Some(film) = optics.coating() {
        CoatingThermalEmission::try_opaque(
            &[CoatingLayer::Dielectric(&film)],
            &[film_temperature],
            optics,
            substrate_temperature,
            ContentHash([84; 32]),
        )
    } else {
        CoatingThermalEmission::try_opaque(
            &[],
            &[],
            optics,
            substrate_temperature,
            ContentHash([84; 32]),
        )
    }
    .unwrap();
    SmithCoatingThermalEmission::try_new(source, optics).unwrap()
}

fn draw(path: u32, event: u32) -> [f64; 3] {
    let bits = fs_rand::philox::philox4x32_10([path, event, 0x7468_6572, 0], [37, 97]);
    [
        (f64::from(bits[0]) + 0.5) / 4_294_967_296.0,
        f64::from(bits[1]) / 4_294_967_296.0,
        f64::from(bits[2]) / 4_294_967_296.0,
    ]
}

fn planck(temperature: f64) -> f64 {
    let lambda = 550e-9_f64;
    2.0 * 6.626_070_15e-34 * 299_792_458.0_f64.powi(2)
        / (lambda.powi(5)
            * (6.626_070_15e-34 * 299_792_458.0 / (lambda * 1.380_649e-23 * temperature)).exp_m1())
        * 1e-9
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:e} expected={expected:e}"
    );
}

fn coated() -> ConductorOptics {
    let bare = metal(0.2);
    bare.with_absorbing_coating(
        CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.7, 0.15).unwrap()
            }),
            125e-9,
            bare.source(),
        )
        .unwrap(),
    )
}

#[test]
fn g1_each_isothermal_path_closes_planck_equilibrium_without_changing_reflection() {
    with_cx(|cx, _| {
        let bare = metal(0.2);
        let transparent = bare.with_coating(
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.45).unwrap(),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let surface = ConductorSurface::try_rough(0.85).unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let temperature = 4300.0;
        let b = planck(temperature);
        for optics in [bare, transparent, coated()] {
            let source = emitter(optics, temperature, temperature);
            for power in [
                Stokes::unpolarized(1.0),
                Stokes {
                    i: 1.0,
                    q: 0.3,
                    u: -0.4,
                    v: 0.5,
                },
                Stokes {
                    i: 1.0,
                    q: 1.0,
                    u: 0.0,
                    v: 0.0,
                },
            ] {
                let camera = camera(Vec3::new(0.8, 0.0, -0.6), power);
                for path in 0..512 {
                    let sample = source
                        .sample(cx, camera, normal, surface, settings, |event| {
                            draw(path, event)
                        })
                        .unwrap();
                    let original = camera
                        .sample_smith_conductor(cx, normal, &optics, surface, settings, |event| {
                            draw(path, event)
                        })
                        .unwrap();
                    assert_eq!(sample.reflected(), original);
                    close(
                        sample.radiance_w_m2_sr_nm() + b * original.beam().power().i,
                        b,
                        2e-11,
                    );
                }
            }
        }
    });
}

#[test]
fn g1_smooth_limit_uses_positive_original_absorption_and_no_draws() {
    with_cx(|cx, _| {
        let optics = metal(1e-18);
        let source = emitter(optics, 4300.0, 4300.0);
        let camera = camera(Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(1.0));
        let sample = source
            .sample(
                cx,
                camera,
                Vec3::new(0.0, 0.0, 1.0),
                ConductorSurface::SMOOTH,
                SmithWalkSettings::try_new(1).unwrap(),
                |_| panic!("smooth thermal source sampled roughness"),
            )
            .unwrap();
        assert!(sample.reflected().is_delta());
        assert_eq!(sample.reflected().beam().power().i, 1.0);
        assert!(sample.radiance_w_m2_sr_nm() > 0.0);
        let positive_absorption = camera
            .reflect_conductor_with_balance(Vec3::new(0.0, 0.0, 1.0), &optics)
            .unwrap()
            .absorbed_power_w();
        close(
            sample.radiance_w_m2_sr_nm(),
            positive_absorption * planck(4300.0),
            3e-11,
        );
    });
}

#[test]
fn g3_each_film_temperature_remains_a_separate_linear_source() {
    with_cx(|cx, _| {
        let optics = coated();
        let sources = [
            emitter(optics, 5200.0, 2200.0),
            emitter(optics, 2200.0, 5200.0),
            emitter(optics, 5200.0, 5200.0),
            emitter(optics, 2200.0, 2200.0),
            emitter(optics, 3700.0, 3700.0),
        ];
        let camera = camera(
            Vec3::new(0.0, 0.0, -1.0),
            Stokes {
                i: 1.0,
                q: 0.3,
                u: -0.4,
                v: 0.5,
            },
        );
        let samples = sources.map(|source| {
            source
                .sample(
                    cx,
                    camera,
                    Vec3::new(0.0, 0.0, 1.0),
                    ConductorSurface::try_rough(1.0).unwrap(),
                    SmithWalkSettings::try_new(2).unwrap(),
                    |event| DRAWS[event as usize],
                )
                .unwrap()
        });
        assert_eq!(samples[0].reflected().facet_events(), 2);
        for sample in &samples[1..] {
            assert_eq!(sample.reflected(), samples[0].reflected());
        }
        let r = samples.map(SmithThermalSample::radiance_w_m2_sr_nm);
        close(r[0] + r[1], r[2] + r[3], 2e-13);
        assert!(((r[0] + r[1]) / (2.0 * r[4]) - 1.0).abs() > 0.1);
        // Both physical temperatures matter; a cold film cannot be replaced
        // by its hot opaque substrate, nor the reverse.
        assert!(r[0] > r[3] && r[1] > r[3] && r[2] > r[0] && r[2] > r[1]);
    });
}

const DRAWS: [[f64; 3]; 3] = [[0.5, 0.8, 0.125], [0.5, 0.1, 0.375], [0.5, 0.2, 0.3]];

#[test]
fn g3_reciprocal_polarization_is_linear_with_the_original_phase_history() {
    with_cx(|cx, _| {
        let source = emitter(coated(), 5000.0, 3000.0);
        let direction = Vec3::new(0.8, 0.0, -0.6);
        let a = Stokes {
            i: 1.0,
            q: 0.0,
            u: 0.6,
            v: 0.8,
        };
        let b = Stokes {
            i: 1.0,
            q: 0.8,
            u: -0.6,
            v: 0.0,
        };
        let mix = Stokes {
            i: 1.0,
            q: 0.6,
            u: -0.3,
            v: 0.2,
        };
        for path in 0..128 {
            let [x, y, mixed] = [a, b, mix].map(|power| {
                source
                    .sample(
                        cx,
                        camera(direction, power),
                        Vec3::new(0.0, 0.0, 1.0),
                        ConductorSurface::try_rough(1.0).unwrap(),
                        SmithWalkSettings::try_new(1024).unwrap(),
                        |event| draw(path, event),
                    )
                    .unwrap()
                    .radiance_w_m2_sr_nm()
            });
            close(mixed, 0.25 * x + 0.75 * y, 2e-12);
        }
    });
}

#[test]
fn g3_anisotropic_physical_frame_rotates_covariantly() {
    with_cx(|cx, _| {
        let source = emitter(coated(), 4400.0, 2800.0);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = camera(
            Vec3::new(0.6, 0.0, -0.8),
            Stokes {
                i: 1.0,
                q: 0.3,
                u: -0.4,
                v: 0.5,
            },
        );
        let rotate = |v: Vec3| Vec3::new(-v.y, v.x, v.z);
        let turned = PolarizedBeam::try_new(
            input.ray().origin,
            PolarizationFrame::try_new(
                rotate(input.frame().direction()),
                rotate(input.frame().horizontal()),
            )
            .unwrap(),
            input.wavelength_nm(),
            input.medium(),
            input.power(),
        )
        .unwrap();
        for path in 0..128 {
            let a = source
                .sample_anisotropic(
                    cx,
                    input,
                    normal,
                    Vec3::new(1.0, 0.0, 0.0),
                    0.8,
                    0.2,
                    SmithWalkSettings::try_new(1024).unwrap(),
                    |event| draw(path, event),
                )
                .unwrap();
            let b = source
                .sample_anisotropic(
                    cx,
                    turned,
                    normal,
                    Vec3::new(0.0, 1.0, 0.0),
                    0.8,
                    0.2,
                    SmithWalkSettings::try_new(1024).unwrap(),
                    |event| draw(path, event),
                )
                .unwrap();
            close(a.radiance_w_m2_sr_nm(), b.radiance_w_m2_sr_nm(), 2e-12);
            assert_eq!(a.reflected().facet_events(), b.reflected().facet_events());
            close(
                a.reflected().beam().power().i,
                b.reflected().beam().power().i,
                2e-12,
            );
        }
    });
}

#[test]
fn g0_matching_vacuum_and_unit_camera_are_explicit_admission() {
    with_cx(|cx, _| {
        let optics = metal(0.2);
        let source = emitter(optics, 4000.0, 4000.0);
        assert!(matches!(
            SmithCoatingThermalEmission::try_new(source.source().clone(), metal(0.3)),
            Err(SmithThermalError::OpticalMismatch)
        ));
        let input = camera(Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(1.0));
        for (medium, power, expected) in [
            (1.0, 0.5, SmithThermalError::UnitCameraRequired),
            (1.33, 1.0, SmithThermalError::VacuumRequired),
        ] {
            let beam = PolarizedBeam::try_new(
                input.ray().origin,
                input.frame(),
                550.0,
                DielectricGlass::new(
                    CauchyIor::try_constant(medium).unwrap(),
                    BeerLambertAbsorption::CLEAR,
                    GlassProvenance::Custom,
                ),
                Stokes::unpolarized(power),
            )
            .unwrap();
            let result = source.sample(
                cx,
                beam,
                Vec3::new(0.0, 0.0, 1.0),
                ConductorSurface::try_rough(1.0).unwrap(),
                SmithWalkSettings::try_new(2).unwrap(),
                |_| panic!("invalid thermal mode consumed a draw"),
            );
            assert_eq!(result, Err(expected));
        }
    });
}

#[test]
fn g4_cancellation_and_budget_never_publish_partial_emission() {
    with_cx(|cx, gate| {
        let source = emitter(coated(), 4500.0, 3000.0);
        let input = camera(Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(1.0));
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let surface = ConductorSurface::try_rough(1.0).unwrap();
        assert_eq!(
            source.sample(
                cx,
                input,
                normal,
                surface,
                SmithWalkSettings::try_new(1).unwrap(),
                |event| DRAWS[event as usize]
            ),
            Err(SmithThermalError::Walk(SmithBeamError::FacetLimit {
                reflected_facets: 1
            }))
        );
        let mut draws = 0;
        let cancelled = source.sample(
            cx,
            input,
            normal,
            surface,
            SmithWalkSettings::try_new(2).unwrap(),
            |event| {
                draws += 1;
                gate.request();
                DRAWS[event as usize]
            },
        );
        assert_eq!(
            cancelled,
            Err(SmithThermalError::Walk(SmithBeamError::Cancelled))
        );
        assert_eq!(draws, 1);
    });
}
