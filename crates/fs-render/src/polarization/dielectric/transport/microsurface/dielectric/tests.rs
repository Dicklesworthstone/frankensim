use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricCoating, DielectricGlass, GlassProvenance,
};
use crate::polarization::Stokes;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_rand::philox::philox4x32_10;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
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
fn glass(n: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}
fn beam(index: f64, direction: Vec3, state: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        550.0,
        glass(index),
        state,
    )
    .unwrap()
}
fn draw(path: u32, event: u32, domain: u32) -> [f64; 4] {
    let block = philox4x32_10([path, event, domain, 0], [0x6173_676c, 0x6469_656c]);
    core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}
fn state_near(a: Stokes, b: Stokes) {
    for (a, b) in [a.i, a.q, a.u, a.v].into_iter().zip([b.i, b.q, b.u, b.v]) {
        near(a, b, 5e-12);
    }
}
#[derive(Default)]
struct Moments {
    n: u32,
    sum: f64,
    square: f64,
}
impl Moments {
    fn push(&mut self, x: f64) {
        self.n += 1;
        self.sum += x;
        self.square += x * x;
    }
    fn mean(&self) -> f64 {
        self.sum / f64::from(self.n)
    }
    fn variance_of_mean(&self) -> f64 {
        ((self.square / f64::from(self.n) - self.mean() * self.mean()).max(0.0)) / f64::from(self.n)
    }
}

#[test]
fn shared_height_flip_retains_small_probabilities_and_transmission_beta_law() {
    for d in [1e-100, 1e-20, 1e-8, 0.1, 0.7, 2.0, 20.0, 100.0, 700.0] {
        let opposite = opposite_depth(d).unwrap();
        assert!(opposite > 0.0 && opposite.is_finite());
        near(opposite_depth(opposite).unwrap() / d, 1.0, 6e-13);
    }
    assert!(opposite_depth(0.0).is_err());
    assert!(opposite_depth(f64::INFINITY).is_err());
    assert!(opposite_depth(1e6).is_err());
    with_cx(|cx, _| {
        let ci = 0.35_f64;
        let direction = Vec3::new((1.0 - ci * ci).sqrt(), 0.0, -ci);
        let input = beam(1.0, direction, Stokes::unpolarized(1.0));
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let model = SmithGgx::try_new(1.0, 1.0).unwrap();
        let step = model.sample(direction, 0.0, [0.5, 0.0, 0.0]).unwrap();
        let SmithEvent::Reflection { facet, .. } = step.event else {
            panic!("first collision")
        };
        let response = DielectricSurface::SMOOTH
            .interface(550.0, false)
            .unwrap()
            .sample_smooth_polarized(facet, direction.scale(-1.0), 1.0, 1.5, 0.5)
            .unwrap();
        assert_eq!(response.event, DielectricEvent::Transmission);
        let ct = -response.direction.z;
        assert!(ct > 0.0);
        let li = 0.5 * (1.0 / ci - 1.0);
        let lt = 0.5 * (1.0 / ct - 1.0);
        // Independent quadrature of Heitz's transmission G2:
        // Beta(1+Lambda_i,1+Lambda_t), not reflection's reciprocal sum.
        let bins = 100_000_u32;
        let beta = (0..bins)
            .map(|i| {
                let h = (f64::from(i) + 0.5) / f64::from(bins);
                h.powf(li) * (1.0 - h).powf(lt)
            })
            .sum::<f64>()
            / f64::from(bins);
        let expected = (1.0 + li) * beta;
        let wrong_reflection_law = (1.0 + li) / (1.0 + li + lt);
        assert!((expected - wrong_reflection_law).abs() > 0.05);
        let count = 32_768_u32;
        let mut escapes = 0_u32;
        for path in 0..count {
            match input.sample_smith_dielectric(
                cx,
                normal,
                glass(1.5),
                DielectricSurface::try_rough(1.0).unwrap(),
                SmithWalkSettings::try_new(1).unwrap(),
                |event| {
                    let mut u = draw(path, event, 1);
                    u[1] = 0.0;
                    u[2] = 0.0;
                    u[3] = 0.5;
                    u
                },
            ) {
                Ok(sample) => {
                    assert_eq!(sample.facet_events(), 1);
                    assert_eq!(sample.event(), DielectricEvent::Transmission);
                    escapes += 1;
                }
                Err(SmithDielectricError::FacetLimit { facet_events: 1 }) => {}
                other => panic!("unexpected first-event outcome {other:?}"),
            }
        }
        let actual = f64::from(escapes) / f64::from(count);
        let tolerance = 7.0 * (expected * (1.0 - expected) / f64::from(count)).sqrt() + 2e-6;
        near(actual, expected, tolerance);
        eprintln!(
            "shared-height first transmission: {actual:.6}, Beta {expected:.6}, reflection law {wrong_reflection_law:.6}"
        );
    });
}

#[test]
fn native_two_sided_walk_conserves_power_and_contains_refractive_returns() {
    with_cx(|cx, _| {
        let surface =
            DielectricSurface::try_anisotropic_rough(0.85, 0.35, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let mut repeated_transmissions = 0;
        let mut multiple_facets = 0;
        let mut elliptical = 0;
        for (inside, cosine, state) in [
            (false, 1.0, Stokes::unpolarized(1.0)),
            (false, 0.3, Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap()),
            (true, 0.45, Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap()),
        ] {
            let cosine: f64 = cosine;
            let sign = if inside { -1.0 } else { 1.0 };
            let initial = if inside { 1.5 } else { 1.0 };
            let target = if inside { 1.0 } else { 1.5 };
            let input = beam(
                initial,
                Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, -sign * cosine),
                state,
            );
            let mut power = Moments::default();
            for path in 0..32_768 {
                let sample = input
                    .sample_smith_dielectric(
                        cx,
                        Vec3::new(0.0, 0.0, 1.0),
                        glass(target),
                        surface,
                        settings,
                        |event| draw(path, event, 2 + u32::from(inside)),
                    )
                    .unwrap();
                let output = sample.beam();
                power.push(output.power().i);
                assert!(sample.log_path_density().is_finite());
                assert!(
                    Stokes::try_new(
                        output.power().i,
                        output.power().q,
                        output.power().u,
                        output.power().v
                    )
                    .is_some()
                );
                let transmitted = sample.event() == DielectricEvent::Transmission;
                assert_eq!(sample.transmissions() % 2 == 1, transmitted);
                assert_eq!(
                    output.medium(),
                    BeamMedium::from(glass(if transmitted { target } else { initial }))
                );
                assert!(
                    output.frame().direction().z * sign * (if transmitted { -1.0 } else { 1.0 })
                        > 0.0
                );
                repeated_transmissions += u32::from(sample.transmissions() >= 2);
                multiple_facets += u32::from(sample.facet_events() >= 3);
                elliptical += u32::from(inside && output.power().v.abs() > 0.1);
            }
            let tolerance = 8.0 * power.variance_of_mean().sqrt() + 2e-10;
            near(power.mean(), 1.0, tolerance);
            eprintln!(
                "white furnace inside={inside} cosine={cosine}: {} +/- {tolerance}",
                power.mean()
            );
        }
        assert!(
            repeated_transmissions > 50,
            "no return through the opposite half-space"
        );
        assert!(multiple_facets > 100, "walk collapsed to single scattering");
        assert!(elliptical > 100, "TIR phase/polarization was lost");
    });
}

#[test]
fn cosine_integrated_transmission_obeys_two_sided_etendue_reciprocity() {
    with_cx(|cx, _| {
        let surface =
            DielectricSurface::try_anisotropic_rough(0.8, 0.3, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let mut totals = [Moments::default(), Moments::default()];
        for (side, total) in totals.iter_mut().enumerate() {
            let (ni, nt, sign) = if side == 0 {
                (1.0, 1.5, 1.0)
            } else {
                (1.5, 1.0, -1.0)
            };
            for path in 0..65_536 {
                let u = draw(path, 0, 10 + side as u32);
                let radius = u[0].sqrt();
                let (s, c) = (2.0 * core::f64::consts::PI * u[1]).sin_cos();
                let direction = Vec3::new(radius * c, radius * s, -sign * (1.0 - u[0]).sqrt());
                let input = beam(ni, direction, Stokes::unpolarized(1.0));
                let sample = input
                    .sample_smith_dielectric(
                        cx,
                        Vec3::new(0.0, 0.0, 1.0),
                        glass(nt),
                        surface,
                        settings,
                        |event| draw(path, event, 12 + side as u32),
                    )
                    .unwrap();
                total.push(if sample.event() == DielectricEvent::Transmission {
                    sample.beam().power().i
                } else {
                    0.0
                });
            }
        }
        let outward = totals[0].mean();
        let inward = 2.25 * totals[1].mean();
        let tolerance = 8.0
            * (totals[0].variance_of_mean() + 2.25 * 2.25 * totals[1].variance_of_mean()).sqrt();
        near(outward, inward, tolerance);
        assert!(outward > 0.5);
        eprintln!(
            "n^2 transmission reciprocity: {outward:.6} versus {inward:.6}, tolerance {tolerance:.6}"
        );
    });
}

#[test]
fn smooth_tir_and_lossy_coated_facets_use_native_mueller_and_reverse_order() {
    with_cx(|cx, _| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let settings = SmithWalkSettings::try_new(8).unwrap();
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([42; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let cell = [
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.08, 85e-9)
                .unwrap(),
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.1).unwrap(), 0.02, 125e-9)
                .unwrap(),
        ];
        let coating =
            DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cell, 1, source).unwrap())
                .unwrap();
        let surface = DielectricSurface::try_rough(0.8)
            .unwrap()
            .with_coating(coating);
        let mut events = [false; 2];
        for inside in [false, true] {
            let (ni, nt, sign) = if inside {
                (1.5, 1.0, -1.0)
            } else {
                (1.0, 1.5, 1.0)
            };
            let input = beam(
                ni,
                normal.scale(-sign),
                Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap(),
            );
            for branch in [0.001, 0.999] {
                let native = surface
                    .interface(550.0, inside)
                    .unwrap()
                    .sample_smooth_polarized(
                        normal.scale(sign),
                        input.frame().direction().scale(-1.0),
                        ni,
                        nt,
                        branch,
                    )
                    .unwrap();
                let expected = input
                    .scatter_smooth_sample(normal.scale(sign), glass(nt).into(), ni, nt, native)
                    .unwrap();
                let sampled = input
                    .sample_smith_dielectric(cx, normal, glass(nt), surface, settings, |_| {
                        [0.5, 0.0, 0.0, branch]
                    })
                    .unwrap();
                assert_eq!(sampled.facet_events(), 1);
                assert_eq!(sampled.event(), expected.event);
                state_near(sampled.beam().power(), expected.beam.power());
                assert_eq!(sampled.beam().medium(), expected.beam.medium());
                near(
                    sampled
                        .beam()
                        .frame()
                        .direction()
                        .dot(expected.beam.frame().direction()),
                    1.0,
                    2e-14,
                );
                assert!(sampled.beam().power().i < 1.0);
                events[usize::from(sampled.event() == DielectricEvent::Transmission)] = true;
            }
        }
        assert_eq!(events, [true, true]);
        let ci = 0.5_f64.sqrt();
        let input = beam(
            1.5,
            Vec3::new(ci, 0.0, ci),
            Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(),
        );
        let sample = input
            .sample_smith_dielectric(
                cx,
                normal,
                glass(1.0),
                DielectricSurface::SMOOTH,
                settings,
                |_| [0.5, 0.2, 0.3, 0.99],
            )
            .unwrap();
        assert_eq!(sample.event(), DielectricEvent::Reflection);
        assert!(sample.is_delta());
        near(sample.beam().power().i, 1.0, 2e-14);
        near(sample.beam().power().v.abs(), 0.6, 2e-14);
    });
}

#[test]
fn null_media_and_coated_equal_index_atoms_keep_topology_without_fake_rough_pdf() {
    with_cx(|cx, _| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = beam(
            1.5,
            Vec3::new(0.8, 0.0, -0.6),
            Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap(),
        );
        let target = DielectricGlass::new(
            CauchyIor::try_constant(1.5).unwrap(),
            BeerLambertAbsorption::try_constant(12.0).unwrap(),
            GlassProvenance::Custom,
        );
        let settings = SmithWalkSettings::try_new(8).unwrap();
        let null = input
            .sample_smith_dielectric(
                cx,
                normal,
                target,
                DielectricSurface::try_rough(1.0).unwrap(),
                settings,
                |_| panic!("null boundary drew a facet"),
            )
            .unwrap();
        assert_eq!(null.event(), DielectricEvent::Transmission);
        assert_eq!(null.facet_events(), 0);
        assert!(null.is_delta());
        assert_eq!(null.beam().frame(), input.frame());
        assert_eq!(null.beam().power(), input.power());
        assert_eq!(null.beam().medium(), target.into());
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([43; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let cell =
            [
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.0, 100e-9)
                    .unwrap(),
            ];
        let surface = DielectricSurface::try_rough(0.8).unwrap().with_coating(
            DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cell, 1, source).unwrap())
                .unwrap(),
        );
        let input = beam(1.5, normal.scale(-1.0), Stokes::unpolarized(1.0));
        let sample = input
            .sample_smith_dielectric(cx, normal, target, surface, settings, |_| {
                [0.5, 0.0, 0.0, 0.999]
            })
            .unwrap();
        assert_eq!(sample.event(), DielectricEvent::Transmission);
        assert!(sample.is_delta());
        assert_eq!(sample.beam().frame().direction(), input.frame().direction());
        assert_eq!(sample.beam().medium(), target.into());
        near(sample.beam().power().i, 1.0, 1e-12);
    });
}

#[test]
fn budgets_cancellation_invalid_draws_and_replay_are_explicit() {
    with_cx(|cx, gate| {
        let input = beam(1.0, Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(1.0));
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let surface = DielectricSurface::try_rough(1.0).unwrap();
        let settings = SmithWalkSettings::try_new(512).unwrap();
        let mut chosen = None;
        for path in 0..256 {
            let sample = input
                .sample_smith_dielectric(cx, normal, glass(1.5), surface, settings, |event| {
                    draw(path, event, 21)
                })
                .unwrap();
            if sample.facet_events() > 1 {
                chosen = Some((path, sample));
                break;
            }
        }
        let (path, expected) = chosen.expect("a genuine multi-facet path");
        assert_eq!(
            input
                .sample_smith_dielectric(cx, normal, glass(1.5), surface, settings, |event| draw(
                    path, event, 21
                ))
                .unwrap(),
            expected
        );
        assert_eq!(
            input.sample_smith_dielectric(
                cx,
                normal,
                glass(1.5),
                surface,
                SmithWalkSettings::try_new(1).unwrap(),
                |event| draw(path, event, 21)
            ),
            Err(SmithDielectricError::FacetLimit { facet_events: 1 })
        );
        for bad in [
            [0.0, 0.2, 0.3, 0.4],
            [0.5, 0.2, 0.3, 1.0],
            [0.5, 0.2, 0.3, f64::NAN],
        ] {
            assert_eq!(
                input.sample_smith_dielectric(cx, normal, glass(1.5), surface, settings, |_| bad),
                Err(SmithDielectricError::Geometry(SmithError::InvalidSample))
            );
        }
        gate.request();
        assert_eq!(
            input.sample_smith_dielectric(cx, normal, glass(1.5), surface, settings, |_| panic!(
                "cancelled walk consumed a draw"
            )),
            Err(SmithDielectricError::Cancelled)
        );
    });
}

#[test]
fn arbitrary_coating_stack_preserves_reverse_order_and_rotated_anisotropy() {
    use crate::conductor::CoherentDielectricCoating;
    with_cx(|cx, _| {
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([44; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let indices = [
            CauchyIor::try_constant(1.3).unwrap(),
            CauchyIor::try_constant(2.1).unwrap(),
        ];
        let thicknesses = [85e-9, 125e-9];
        let films = core::array::from_fn::<_, 2, _>(|i| {
            CoherentDielectricCoating::try_new(indices[i], thicknesses[i], source).unwrap()
        });
        let stack = DielectricCoatingStack::try_new(&films).unwrap();
        let cell = core::array::from_fn::<_, 2, _>(|i| {
            CauchyCoatingLayer::try_new(indices[i], 0.0, thicknesses[i]).unwrap()
        });
        let compact =
            DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cell, 1, source).unwrap())
                .unwrap();
        let surface =
            DielectricSurface::try_anisotropic_rough(0.8, 0.3, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let rotate = |v: Vec3| Vec3::new(-v.y, v.x, v.z);
        let rotated = surface.with_tangent_axis(Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        for inside in [false, true] {
            let (ni, nt, sign) = if inside {
                (1.5, 1.0, -1.0)
            } else {
                (1.0, 1.5, 1.0)
            };
            let input = beam(
                ni,
                Vec3::new(0.8, 0.0, -sign * 0.6),
                Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap(),
            );
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
                let expected = input
                    .sample_smith_dielectric(
                        cx,
                        Vec3::new(0.0, 0.0, 1.0),
                        glass(nt),
                        surface.with_coating(compact),
                        settings,
                        |event| draw(path, event, 40),
                    )
                    .unwrap();
                let actual = turned
                    .sample_smith_stack_dielectric(
                        cx,
                        Vec3::new(0.0, 0.0, 1.0),
                        glass(nt),
                        rotated,
                        stack,
                        settings,
                        |event| draw(path, event, 40),
                    )
                    .unwrap();
                assert_eq!(actual.event(), expected.event());
                assert_eq!(actual.facet_events(), expected.facet_events());
                assert_eq!(actual.transmissions(), expected.transmissions());
                near(
                    actual
                        .beam()
                        .frame()
                        .direction()
                        .dot(rotate(expected.beam().frame().direction())),
                    1.0,
                    3e-13,
                );
                state_near(actual.beam().power(), expected.beam().power());
                near(
                    actual.log_path_density(),
                    expected.log_path_density(),
                    2e-11,
                );
            }
        }
    });
}
