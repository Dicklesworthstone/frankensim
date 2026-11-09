use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoherentDielectricCoating, CoherentPeriodicCoating, ConductorDataStatus,
    ConductorSource,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricCoating, DielectricGlass, GlassProvenance,
};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_rand::philox::philox4x32_10;

mod mixed;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 61,
                    kernel_id: 9,
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
        source_frame(direction),
        550.0,
        glass(index),
        state,
    )
    .unwrap()
}
fn source_frame(direction: Vec3) -> PolarizationFrame {
    PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap()
}
fn draw(path: u32, event: u32, domain: u32) -> [f64; 4] {
    philox4x32_10([path, event, domain, 0], [0x6173_676c, 0x6e65_6506])
        .map(|v| (f64::from(v) + 0.5) / 4_294_967_296.0)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
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
    fn variance(&self) -> f64 {
        (self.square / f64::from(self.n) - self.mean() * self.mean()).max(0.0) / f64::from(self.n)
    }
}
// A bounded source away from the geometric horizon. The separately tested
// all-underflow grazing path is an explicit numerical refusal, never averaged
// as a rejected history. The same source function weights BOTH estimators.
fn sky(v: Vec3) -> f64 {
    if v.z.abs() < 0.02 {
        0.0
    } else {
        0.8 + 0.3 * v.z + 0.2 * v.x * v.x
    }
}
fn sphere(path: u32, domain: u32) -> Vec3 {
    let u = draw(path, 0, domain);
    let z = 2.0 * u[0] - 1.0;
    let radius = (1.0 - z * z).sqrt();
    let (s, c) = det::sin_cos(core::f64::consts::TAU * u[1]);
    Vec3::new(radius * c, radius * s, z)
}

#[test]
fn first_transmission_connection_uses_beta_height_not_reflection_masking() {
    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
    let ci = 0.35_f64;
    let incoming = Vec3::new((1.0 - ci * ci).sqrt(), 0.0, -ci);
    let ni = 1.0;
    let nt = 1.5;
    let sample = DielectricSurface::SMOOTH
        .interface(550.0, false)
        .unwrap()
        .sample_smooth_polarized(Vec3::new(0.0, 0.0, 1.0), incoming.scale(-1.0), ni, nt, 0.5)
        .unwrap();
    let outgoing = sample.direction;
    let co = -outgoing.z;
    let li = 0.5 * (1.0 / ci - 1.0);
    let lo = 0.5 * (1.0 / co - 1.0);
    let bins = 131_072_u32;
    let beta = (0..bins)
        .map(|i| {
            let h = (f64::from(i) + 0.5) / f64::from(bins);
            h.powf(li) * (1.0 - h).powf(lo)
        })
        .sum::<f64>()
        / f64::from(bins);
    let facet = Vec3::new(0.0, 0.0, 1.0);
    let denominator = ni * ci - nt * co;
    let prefactor = model.ndf(facet) * ci * nt * nt * co
        / (denominator * denominator * model.projected_area(incoming.scale(-1.0)) * co);
    let expected = prefactor * (1.0 + li) * beta;
    let wrong = prefactor * (1.0 + li) / (1.0 + li + lo);
    assert!((wrong - expected).abs() > 0.05 * expected);
    let integrated = (0..bins)
        .map(|i| {
            let u = (f64::from(i) + 0.5) / f64::from(bins);
            let depth = -u.ln() / (1.0 + li);
            det::exp(
                connection_at_depth(model, incoming, Some(depth), outgoing, ni, nt, false)
                    .unwrap()
                    .unwrap()
                    .log_factor_per_sr,
            )
        })
        .sum::<f64>()
        / f64::from(bins);
    near(integrated, expected, 8e-6 * expected);
    let analytic = det::exp(
        connection_at_depth(model, incoming, None, outgoing, ni, nt, false)
            .unwrap()
            .unwrap()
            .log_factor_per_sr,
    );
    near(analytic, expected, 8e-6 * expected);
    eprintln!(
        "glass first-connection Beta: {integrated:.9}, expected {expected:.9}, wrong reflection {wrong:.9}"
    );
    let reflected = Vec3::new(incoming.x, 0.0, ci);
    let reflected_mean = (0..bins)
        .map(|i| {
            let depth = -((f64::from(i) + 0.5) / f64::from(bins)).ln() / (1.0 + li);
            det::exp(
                connection_at_depth(model, incoming, Some(depth), reflected, ni, nt, true)
                    .unwrap()
                    .unwrap()
                    .log_factor_per_sr,
            )
        })
        .sum::<f64>()
        / f64::from(bins);
    near(
        reflected_mean,
        model
            .single_scattering_connection(incoming, reflected)
            .unwrap()
            .unwrap()
            .factor_per_sr,
        1e-7,
    );
}

#[test]
fn full_polarized_directional_integral_matches_independent_walk_from_both_sides() {
    with_cx(|cx, _| {
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let source = Stokes::try_new(1.0, 0.1, 0.5, -0.2).unwrap();
        for (inside, ci) in [(false, 0.7_f64), (true, 0.45_f64)] {
            let (ni, nt, sign) = if inside {
                (1.5, 1.0, -1.0)
            } else {
                (1.0, 1.5, 1.0)
            };
            let input = beam(
                ni,
                Vec3::new((1.0 - ci * ci).sqrt(), 0.0, -sign * ci),
                Stokes::try_new(1.0, 0.3, -0.2, 0.4).unwrap(),
            );
            let surface =
                DielectricSurface::try_anisotropic_rough(0.85, 0.55, Vec3::new(1.0, 0.0, 0.0))
                    .unwrap();
            let mut traced = [Moments::default(), Moments::default()];
            let mut evaluated = [Moments::default(), Moments::default()];
            let mut multi = 0;
            let mut partition_difference = Moments::default();
            let p = 1.0 / (4.0 * core::f64::consts::PI);
            for path in 0..262_144_u32 {
                let walk = input
                    .sample_smith_dielectric(cx, normal, glass(nt), surface, settings, |event| {
                        draw(path, event, 10 + u32::from(inside))
                    })
                    .unwrap();
                let out = walk.beam().frame().direction();
                let transmitted = walk.event() == DielectricEvent::Transmission;
                let eta = if transmitted { (ni / nt).powi(2) } else { 1.0 };
                let response =
                    reciprocal_contraction(walk.beam(), source_frame(out.scale(-1.0)), source)
                        .unwrap();
                let path_value = eta * response * sky(out);
                for side in 0..2 {
                    traced[side].push(if side == usize::from(transmitted) {
                        path_value
                    } else {
                        0.0
                    });
                }
                let q = input
                    .evaluate_dielectric_connection(
                        normal,
                        glass(nt),
                        surface,
                        source_frame(out.scale(-1.0)),
                        Stokes::unpolarized(1.0),
                    )
                    .unwrap()
                    .pdf;
                let camera_weight = q * q / (p * p + q * q);
                let weighted_path = path_value * camera_weight;
                let out = sphere(path, 20 + u32::from(inside));
                if sky(out) == 0.0 {
                    for m in &mut evaluated {
                        m.push(0.0);
                    }
                    partition_difference.push(weighted_path - path_value);
                    continue;
                }
                let response = input
                    .evaluate_smith_dielectric_connection(
                        cx,
                        normal,
                        glass(nt),
                        surface,
                        source_frame(out.scale(-1.0)),
                        source,
                        settings,
                        |event| draw(path, event, 30 + u32::from(inside)),
                    )
                    .unwrap_or_else(|e| panic!("inside={inside}, path={path}, out={out:?}: {e:?}"));
                multi += usize::from(response.transmissions() >= 2);
                let side = usize::from(response.event() == DielectricEvent::Transmission);
                let value =
                    response.factor_per_sr() * out.z.abs() * 4.0 * core::f64::consts::PI * sky(out);
                for component in 0..2 {
                    evaluated[component].push(if side == component { value } else { 0.0 });
                }
                let q = input
                    .evaluate_dielectric_connection(
                        normal,
                        glass(nt),
                        surface,
                        source_frame(out.scale(-1.0)),
                        Stokes::unpolarized(1.0),
                    )
                    .unwrap()
                    .pdf;
                let light_weight = p * p / (p * p + q * q);
                partition_difference.push(weighted_path + value * light_weight - path_value);
            }
            assert!(multi > 0);
            let partition_tolerance = 7.0 * partition_difference.variance().sqrt() + 3e-5;
            near(partition_difference.mean(), 0.0, partition_tolerance);
            eprintln!(
                "glass matched angular partition inside={inside}: delta {:.8}, tolerance {:.8}",
                partition_difference.mean(),
                partition_tolerance
            );
            for side in 0..2 {
                let tolerance =
                    7.0 * (traced[side].variance() + evaluated[side].variance()).sqrt() + 3e-5;
                near(traced[side].mean(), evaluated[side].mean(), tolerance);
                eprintln!(
                    "glass connection side={inside}/{side}: path {:.8}, eval {:.8}, tolerance {:.8}",
                    traced[side].mean(),
                    evaluated[side].mean(),
                    tolerance
                );
            }
        }
    });
}

#[test]
fn reciprocal_transmission_retains_etendue_and_near_matched_contrast() {
    with_cx(|cx, _| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let a = unit(Vec3::new(0.4, 0.2, 1.0)).unwrap();
        let b = unit(Vec3::new(-0.3, 0.15, -1.0)).unwrap();
        let surface =
            DielectricSurface::try_anisotropic_rough(0.9, 0.6, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let mut forward = Moments::default();
        let mut backward = Moments::default();
        for path in 0..32_768_u32 {
            let f = beam(1.0, a.scale(-1.0), Stokes::unpolarized(1.0))
                .evaluate_smith_dielectric_connection(
                    cx,
                    normal,
                    glass(1.5),
                    surface,
                    source_frame(b.scale(-1.0)),
                    Stokes::unpolarized(1.0),
                    settings,
                    |e| draw(path, e, 40),
                )
                .unwrap();
            let r = beam(1.5, b.scale(-1.0), Stokes::unpolarized(1.0))
                .evaluate_smith_dielectric_connection(
                    cx,
                    normal,
                    glass(1.0),
                    surface,
                    source_frame(a.scale(-1.0)),
                    Stokes::unpolarized(1.0),
                    settings,
                    |e| draw(path, e, 41),
                )
                .unwrap();
            forward.push(f.factor_per_sr());
            backward.push(r.factor_per_sr() / 2.25);
        }
        let tolerance = 7.0 * (forward.variance() + backward.variance()).sqrt() + 1e-5;
        near(forward.mean(), backward.mean(), tolerance);
        eprintln!(
            "glass reciprocity f/n²: {:.8} vs {:.8}, tolerance {:.8}",
            forward.mean(),
            backward.mean(),
            tolerance
        );
        for nt in [1.0 + f64::EPSILON, 1.0 + 1e-10, 1.001] {
            let c = connection_at_depth(
                SmithGgx::try_new(0.7, 0.2).unwrap(),
                Vec3::new(0.0, 0.0, -1.0),
                Some(0.4),
                Vec3::new(0.0, 0.0, -1.0),
                1.0,
                nt,
                false,
            )
            .unwrap()
            .unwrap();
            let expected = det::ln(1.0 / (core::f64::consts::PI * 0.7 * 0.2)) + 2.0 * det::ln(nt)
                - 2.0 * det::ln(nt - 1.0);
            near(c.log_factor_per_sr, expected, 2e-13);
        }
    });
}

#[test]
fn delta_mixed_measure_replay_cancellation_and_facet_refusals_are_explicit() {
    with_cx(|cx, gate| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = beam(1.0, Vec3::new(0.0, 0.0, -1.0), Stokes::unpolarized(1.0));
        let source = source_frame(Vec3::new(0.0, 0.0, 1.0));
        let settings = SmithWalkSettings::try_new(512).unwrap();
        let rough = DielectricSurface::try_rough(0.9).unwrap();
        assert!(
            !input
                .smith_dielectric_connection_is_continuous(glass(1.0), rough)
                .unwrap()
        );
        assert!(
            input
                .smith_dielectric_connection_is_continuous(glass(1.5), rough)
                .unwrap()
        );
        for surface in [rough, DielectricSurface::SMOOTH] {
            let result = input
                .evaluate_smith_dielectric_connection(
                    cx,
                    normal,
                    glass(1.0),
                    surface,
                    source,
                    Stokes::unpolarized(1.0),
                    settings,
                    |_| panic!("delta drew a hidden facet"),
                )
                .unwrap();
            assert_eq!(result.factor_per_sr(), 0.0);
            assert_eq!(result.facet_events(), 0);
        }
        let empty: [CoherentDielectricCoating; 0] = [];
        let stack = DielectricCoatingStack::try_new(&empty).unwrap();
        let empty_response = input
            .evaluate_smith_stack_dielectric_connection(
                cx,
                normal,
                glass(1.0),
                rough,
                stack,
                source,
                Stokes::unpolarized(1.0),
                settings,
                |_| panic!("empty index-matched stack drew a facet"),
            )
            .unwrap();
        assert_eq!(empty_response.factor_per_sr(), 0.0);
        let run = || {
            input.evaluate_smith_dielectric_connection(
                cx,
                normal,
                glass(1.5),
                rough,
                source,
                Stokes::unpolarized(1.0),
                settings,
                |e| draw(97, e, 50),
            )
        };
        assert_eq!(run(), run());
        let bad = input.evaluate_smith_dielectric_connection(
            cx,
            normal,
            glass(1.5),
            rough,
            source,
            Stokes::unpolarized(1.0),
            settings,
            |_| [0.0, 0.5, 0.5, 0.5],
        );
        assert_eq!(
            bad,
            Err(SmithDielectricError::Geometry(SmithError::InvalidSample))
        );
        let mut refused = false;
        for path in 0..128 {
            let result = input.evaluate_smith_dielectric_connection(
                cx,
                normal,
                glass(1.5),
                rough,
                source,
                Stokes::unpolarized(1.0),
                SmithWalkSettings::try_new(1).unwrap(),
                |e| draw(path, e, 51),
            );
            if matches!(
                result,
                Err(SmithDielectricError::FacetLimit { facet_events: 1 })
            ) {
                refused = true;
                break;
            }
        }
        assert!(refused);
        let cancelled = input.evaluate_smith_dielectric_connection(
            cx,
            normal,
            glass(1.5),
            rough,
            source,
            Stokes::unpolarized(1.0),
            settings,
            |e| {
                gate.request();
                draw(0, e, 52)
            },
        );
        assert_eq!(cancelled, Err(SmithDielectricError::Cancelled));
        let already = input.evaluate_smith_dielectric_connection(
            cx,
            normal,
            glass(1.5),
            rough,
            source,
            Stokes::unpolarized(1.0),
            settings,
            |_| panic!("cancelled call drew"),
        );
        assert_eq!(already, Err(SmithDielectricError::Cancelled));
    });
}

#[test]
fn log_beta_oracles_and_grazing_range_boundary_are_explicit() {
    // 70-digit Decimal log-Gamma references (recurrence to 100 followed by
    // ten Bernoulli terms), generated outside the production evaluator.
    for (a, b, expected) in [
        (1.0, 1e9, -20.723265836946411),
        (1.125, 1.75, -0.725360279169711),
        (1.2, 20000.0, -11.96956515297667),
        (2.5, 1e12, -68.79286991935033),
        (8.0, 8.0, -10.848948661710063),
        (20.0, 60.0, -45.41738460237084),
        (1.000000001, 1e20, -46.051701906509835),
    ] {
        near(log_beta(a, b).unwrap(), expected, 2e-13);
        near(log_beta(b, a).unwrap(), expected, 2e-13);
    }
    with_cx(|cx, _| {
        let input = beam(
            1.5,
            Vec3::new((1.0 - 0.45_f64.powi(2)).sqrt(), 0.0, 0.45),
            Stokes::try_new(1.0, 0.3, -0.2, 0.4).unwrap(),
        );
        let out = sphere(53859, 21);
        near(out.z, 0.0005819711368530989, 0.0);
        let result = input.evaluate_smith_dielectric_connection(
            cx,
            Vec3::new(0.0, 0.0, 1.0),
            glass(1.0),
            DielectricSurface::try_anisotropic_rough(0.85, 0.55, Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            source_frame(out.scale(-1.0)),
            Stokes::try_new(1.0, 0.1, 0.5, -0.2).unwrap(),
            SmithWalkSettings::try_new(1024).unwrap(),
            |e| draw(53859, e, 31),
        );
        assert_eq!(
            result,
            Err(SmithDielectricError::Beam(BeamError::NumericalRange))
        );
    });
}

#[test]
fn lossy_films_preserve_polarized_linearity_reversal_and_source_integral() {
    with_cx(|cx, _| {
        let provenance = ConductorSource::try_new(
            fs_blake3::ContentHash([62; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let cell = [
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.08, 85e-9)
                .unwrap(),
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.1).unwrap(), 0.02, 125e-9)
                .unwrap(),
        ];
        let coating = DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(&cell, 1, provenance).unwrap(),
        )
        .unwrap();
        let surface =
            DielectricSurface::try_anisotropic_rough(0.85, 0.55, Vec3::new(1.0, 0.0, 0.0))
                .unwrap()
                .with_coating(coating);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let source = Stokes::try_new(1.0, 0.2, -0.3, 0.4).unwrap();
        let mut totals = [0.0; 2];
        for inside in [false, true] {
            let (ni, nt, sign) = if inside {
                (1.5, 1.0, -1.0)
            } else {
                (1.0, 1.5, 1.0)
            };
            let direction = Vec3::new(0.6, 0.0, -sign * 0.8);
            let input = beam(ni, direction, Stokes::try_new(1.0, 0.2, 0.3, -0.1).unwrap());
            let mut traced = Moments::default();
            let mut evaluated = Moments::default();
            for path in 0..65_536 {
                let walk = input
                    .sample_smith_dielectric(cx, normal, glass(nt), surface, settings, |e| {
                        draw(path, e, 70 + u32::from(inside))
                    })
                    .unwrap();
                let out = walk.beam().frame().direction();
                let eta = if walk.event() == DielectricEvent::Transmission {
                    (ni / nt).powi(2)
                } else {
                    1.0
                };
                traced.push(
                    eta * reciprocal_contraction(
                        walk.beam(),
                        source_frame(out.scale(-1.0)),
                        source,
                    )
                    .unwrap()
                        * sky(out),
                );
                let out = sphere(path, 80 + u32::from(inside));
                let value = if sky(out) == 0.0 {
                    0.0
                } else {
                    input
                        .evaluate_smith_dielectric_connection(
                            cx,
                            normal,
                            glass(nt),
                            surface,
                            source_frame(out.scale(-1.0)),
                            source,
                            settings,
                            |e| draw(path, e, 90 + u32::from(inside)),
                        )
                        .unwrap()
                        .factor_per_sr()
                        * out.z.abs()
                        * 4.0
                        * core::f64::consts::PI
                        * sky(out)
                };
                evaluated.push(value);
            }
            let tolerance = 7.0 * (traced.variance() + evaluated.variance()).sqrt() + 2e-5;
            near(traced.mean(), evaluated.mean(), tolerance);
            totals[usize::from(inside)] = evaluated.mean();
            eprintln!(
                "lossy film side={inside}: path {:.8}, eval {:.8}, tolerance {:.8}",
                traced.mean(),
                evaluated.mean(),
                tolerance
            );
            let source_frame = source_frame(unit(Vec3::new(-0.3, 0.2, sign)).unwrap());
            for path in 0..64 {
                let evaluate = |state| {
                    beam(ni, direction, state)
                        .evaluate_smith_dielectric_connection(
                            cx,
                            normal,
                            glass(nt),
                            surface,
                            source_frame,
                            source,
                            settings,
                            |e| draw(path, e, 100),
                        )
                        .unwrap()
                };
                let a = evaluate(Stokes::try_new(1.0, 0.6, -0.4, 0.2).unwrap());
                let b = evaluate(Stokes::try_new(1.0, -0.6, 0.4, -0.2).unwrap());
                let u = evaluate(Stokes::unpolarized(1.0));
                assert_eq!(
                    (a.facet_events(), a.transmissions()),
                    (b.facet_events(), b.transmissions())
                );
                near(
                    0.5 * (a.factor_per_sr() + b.factor_per_sr()),
                    u.factor_per_sr(),
                    2e-12 * (1.0 + u.factor_per_sr()),
                );
            }
        }
        assert!((totals[0] - totals[1]).abs() > 0.01);
    });
}

#[test]
fn arbitrary_stack_and_rotated_anisotropy_share_the_same_connection() {
    with_cx(|cx, _| {
        let provenance = ConductorSource::try_new(
            fs_blake3::ContentHash([63; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let indices = [
            CauchyIor::try_constant(1.3).unwrap(),
            CauchyIor::try_constant(2.1).unwrap(),
        ];
        let thickness = [85e-9, 125e-9];
        let films = core::array::from_fn::<_, 2, _>(|i| {
            CoherentDielectricCoating::try_new(indices[i], thickness[i], provenance).unwrap()
        });
        let stack = DielectricCoatingStack::try_new(&films).unwrap();
        let cells = core::array::from_fn::<_, 2, _>(|i| {
            CauchyCoatingLayer::try_new(indices[i], 0.0, thickness[i]).unwrap()
        });
        let compact = DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(&cells, 1, provenance).unwrap(),
        )
        .unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let surface =
            DielectricSurface::try_anisotropic_rough(0.8, 0.3, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let rotated = surface.with_tangent_axis(Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let settings = SmithWalkSettings::try_new(1024).unwrap();
        let rotate = |v: Vec3| Vec3::new(-v.y, v.x, v.z);
        for inside in [false, true] {
            let (ni, nt, sign) = if inside {
                (1.5, 1.0, -1.0)
            } else {
                (1.0, 1.5, 1.0)
            };
            let input = beam(
                ni,
                Vec3::new(0.6, 0.0, -sign * 0.8),
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
            for side in [-1.0, 1.0] {
                let light = source_frame(unit(Vec3::new(0.3, 0.1, side)).unwrap());
                let turned_light = PolarizationFrame::try_new(
                    rotate(light.direction()),
                    rotate(light.horizontal()),
                )
                .unwrap();
                for path in 0..128 {
                    let a = input
                        .evaluate_smith_dielectric_connection(
                            cx,
                            normal,
                            glass(nt),
                            surface.with_coating(compact),
                            light,
                            Stokes::try_new(1.0, -0.2, 0.6, 0.1).unwrap(),
                            settings,
                            |e| draw(path, e, 110),
                        )
                        .unwrap();
                    let b = turned
                        .evaluate_smith_stack_dielectric_connection(
                            cx,
                            normal,
                            glass(nt),
                            rotated,
                            stack,
                            turned_light,
                            Stokes::try_new(1.0, -0.2, 0.6, 0.1).unwrap(),
                            settings,
                            |e| draw(path, e, 110),
                        )
                        .unwrap();
                    assert_eq!(
                        (a.event(), a.facet_events(), a.transmissions()),
                        (b.event(), b.facet_events(), b.transmissions())
                    );
                    near(
                        a.factor_per_sr(),
                        b.factor_per_sr(),
                        3e-11 * (1.0 + a.factor_per_sr()),
                    );
                }
            }
        }
    });
}
