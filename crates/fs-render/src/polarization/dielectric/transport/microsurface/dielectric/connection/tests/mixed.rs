use super::*;

fn coating(loss: f64) -> DielectricCoating {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([71; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let cells = [
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.4).unwrap(), loss, 65e-9).unwrap(),
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.7).unwrap(), 0.4 * loss, 95e-9)
            .unwrap(),
    ];
    DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cells, 1, source).unwrap())
        .unwrap()
}

#[test]
fn matched_first_reflection_is_analytic_and_never_contains_the_null_atom() {
    with_cx(|cx, _| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = beam(
            1.0,
            normal.scale(-1.0),
            Stokes::try_new(1.0, 0.3, -0.4, 0.2).unwrap(),
        );
        let source = Stokes::try_new(1.0, -0.2, 0.5, -0.1).unwrap();
        let surface = DielectricSurface::try_rough(0.8)
            .unwrap()
            .with_coating(coating(0.03));
        let source_frame = source_frame(normal.scale(-1.0));
        let optics = surface
            .interface(550.0, false)
            .unwrap()
            .fresnel_polarized(1.0, 1.0, 1.0)
            .unwrap();
        let reflected = input
            .scatter(
                normal,
                normal,
                optics.reflection,
                DielectricEvent::Reflection,
                glass(1.0).into(),
            )
            .unwrap();
        let expected = reciprocal_contraction(reflected, source_frame, source).unwrap()
            / (4.0 * core::f64::consts::PI * 0.8 * 0.8);
        let mut ordinals = Vec::new();
        let response = input
            .evaluate_smith_dielectric_connection(
                cx,
                normal,
                glass(1.0),
                surface,
                source_frame,
                source,
                SmithWalkSettings::try_new(16).unwrap(),
                |ordinal| {
                    ordinals.push(ordinal);
                    // Normal facet, transmission continuation; both outward normal
                    // flights escape. The sole nonzero connection is analytic.
                    [0.5, 0.0, 0.0, 0.999]
                },
            )
            .unwrap();
        near(response.factor_per_sr(), expected, 2e-14);
        assert_eq!(ordinals, [0, 1, 2]);
        assert_eq!((response.facet_events(), response.transmissions()), (1, 1));
        assert!(
            !input
                .smith_dielectric_connection_is_continuous(glass(1.0), surface)
                .unwrap()
        );
    });
}

#[test]
fn matched_continuous_integral_and_separate_atom_match_full_polarized_walk() {
    with_cx(|cx, _| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let surface = DielectricSurface::try_anisotropic_rough(1.0, 0.65, Vec3::new(1.0, 0.0, 0.0))
            .unwrap()
            .with_coating(coating(0.015));
        let settings = SmithWalkSettings::try_new(4096).unwrap();
        let source = Stokes::try_new(1.0, 0.1, -0.4, 0.3).unwrap();
        for inside in [false, true] {
            let sign = if inside { -1.0 } else { 1.0 };
            let input = beam(
                1.0,
                Vec3::new(0.8, 0.0, -sign * 0.6),
                Stokes::try_new(1.0, -0.3, 0.2, 0.4).unwrap(),
            );
            let mut traced = [Moments::default(), Moments::default()];
            let mut evaluated = [Moments::default(), Moments::default()];
            let mut atom = Moments::default();
            let mut complete = Moments::default();
            for path in 0..131_072 {
                let walk = input
                    .sample_smith_dielectric(cx, normal, glass(1.0), surface, settings, |e| {
                        draw(path, e, 171 + u32::from(inside))
                    })
                    .unwrap();
                let out = walk.beam().frame().direction();
                let value =
                    reciprocal_contraction(walk.beam(), source_frame(out.scale(-1.0)), source)
                        .unwrap()
                        * sky(out);
                complete.push(value);
                atom.push(if walk.is_delta() { value } else { 0.0 });
                let side = usize::from(walk.event() == DielectricEvent::Transmission);
                for component in 0..2 {
                    traced[component].push(if !walk.is_delta() && side == component {
                        value
                    } else {
                        0.0
                    });
                }
                let out = sphere(path, 181 + u32::from(inside));
                let side = usize::from(sign * out.z < 0.0);
                let value = if sky(out) == 0.0 {
                    0.0
                } else {
                    input
                        .evaluate_smith_dielectric_connection(
                            cx,
                            normal,
                            glass(1.0),
                            surface,
                            source_frame(out.scale(-1.0)),
                            source,
                            settings,
                            |e| draw(path, e, 191 + u32::from(inside)),
                        )
                        .unwrap_or_else(|error| {
                            panic!("inside={inside} path={path} out={out:?}: {error:?}")
                        })
                        .factor_per_sr()
                        * out.z.abs()
                        * 4.0
                        * core::f64::consts::PI
                        * sky(out)
                };
                for component in 0..2 {
                    evaluated[component].push(if side == component { value } else { 0.0 });
                }
            }
            assert!(atom.mean() > 0.1);
            assert!(
                traced[1].mean() > 0.01,
                "continuous transmitted component disappeared"
            );
            for side in 0..2 {
                let tolerance =
                    7.0 * (traced[side].variance() + evaluated[side].variance()).sqrt() + 2e-5;
                near(traced[side].mean(), evaluated[side].mean(), tolerance);
                eprintln!(
                    "matched continuous inside={inside} side={side}: path {:.8}, eval {:.8}, tolerance {:.8}",
                    traced[side].mean(),
                    evaluated[side].mean(),
                    tolerance
                );
            }
            let recombined = atom.mean() + evaluated[0].mean() + evaluated[1].mean();
            let tolerance = 7.0
                * (complete.variance()
                    + atom.variance()
                    + evaluated.iter().map(Moments::variance).sum::<f64>())
                .sqrt()
                + 3e-5;
            near(recombined, complete.mean(), tolerance);
            eprintln!(
                "matched atom inside={inside}: {:.8}, full {:.8}, recombined {:.8}",
                atom.mean(),
                complete.mean(),
                recombined
            );
        }
    });
}

#[test]
fn matched_transmitted_response_is_reciprocal() {
    with_cx(|cx, _| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let a = unit(Vec3::new(0.7, -0.2, 1.0)).unwrap();
        let b = unit(Vec3::new(0.4, 0.3, -1.0)).unwrap();
        let surface = DielectricSurface::try_anisotropic_rough(0.95, 0.6, Vec3::new(1.0, 0.0, 0.0))
            .unwrap()
            .with_coating(coating(0.02));
        let settings = SmithWalkSettings::try_new(4096).unwrap();
        let mut forward = Moments::default();
        let mut backward = Moments::default();
        for path in 0..65_536 {
            for (direction, out, domain, result) in [
                (a.scale(-1.0), b, 201, &mut forward),
                (b.scale(-1.0), a, 202, &mut backward),
            ] {
                result.push(
                    beam(1.0, direction, Stokes::unpolarized(1.0))
                        .evaluate_smith_dielectric_connection(
                            cx,
                            normal,
                            glass(1.0),
                            surface,
                            source_frame(out.scale(-1.0)),
                            Stokes::unpolarized(1.0),
                            settings,
                            |e| draw(path, e, domain),
                        )
                        .unwrap()
                        .factor_per_sr(),
                );
            }
        }
        let tolerance = 7.0 * (forward.variance() + backward.variance()).sqrt() + 1e-5;
        near(forward.mean(), backward.mean(), tolerance);
        assert!(forward.mean() > 0.01);
        eprintln!(
            "matched reciprocal transmission: {:.8}, {:.8}, tolerance {:.8}",
            forward.mean(),
            backward.mean(),
            tolerance
        );
    });
}

#[test]
fn matched_stack_frame_covariance_and_polarized_linearity_share_native_history() {
    with_cx(|cx, _| {
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([72; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let index = CauchyIor::try_constant(2.1).unwrap();
        let films = [CoherentDielectricCoating::try_new(index, 80e-9, source).unwrap()];
        let stack = DielectricCoatingStack::try_new(&films).unwrap();
        let cell = [CauchyCoatingLayer::try_new(index, 0.0, 80e-9).unwrap()];
        let compact =
            DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cell, 1, source).unwrap())
                .unwrap();
        let surface =
            DielectricSurface::try_anisotropic_rough(0.9, 0.4, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let turned_surface = surface.with_tangent_axis(Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let rotate = |v: Vec3| Vec3::new(-v.y, v.x, v.z);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let direction = Vec3::new(0.6, 0.0, -0.8);
        let settings = SmithWalkSettings::try_new(2048).unwrap();
        let source_state = Stokes::try_new(1.0, 0.2, -0.4, 0.1).unwrap();
        for side in [-1.0, 1.0] {
            let light = source_frame(unit(Vec3::new(0.3, 0.2, side)).unwrap());
            let turned_light =
                PolarizationFrame::try_new(rotate(light.direction()), rotate(light.horizontal()))
                    .unwrap();
            for path in 0..64 {
                let evaluate = |state| {
                    beam(1.0, direction, state)
                        .evaluate_smith_dielectric_connection(
                            cx,
                            normal,
                            glass(1.0),
                            surface.with_coating(compact),
                            light,
                            source_state,
                            settings,
                            |e| draw(path, e, 211),
                        )
                        .unwrap()
                };
                let plus = evaluate(Stokes::try_new(1.0, 0.5, -0.3, 0.2).unwrap());
                let minus = evaluate(Stokes::try_new(1.0, -0.5, 0.3, -0.2).unwrap());
                let unpolarized = evaluate(Stokes::unpolarized(1.0));
                near(
                    0.5 * (plus.factor_per_sr() + minus.factor_per_sr()),
                    unpolarized.factor_per_sr(),
                    3e-12,
                );
                assert_eq!(
                    (plus.facet_events(), plus.transmissions()),
                    (minus.facet_events(), minus.transmissions())
                );
                let original = beam(1.0, direction, Stokes::unpolarized(1.0));
                let turned = PolarizedBeam::try_new(
                    original.ray().origin,
                    PolarizationFrame::try_new(
                        rotate(direction),
                        rotate(original.frame().horizontal()),
                    )
                    .unwrap(),
                    550.0,
                    glass(1.0),
                    original.power(),
                )
                .unwrap();
                let result = turned
                    .evaluate_smith_stack_dielectric_connection(
                        cx,
                        normal,
                        glass(1.0),
                        turned_surface,
                        stack,
                        turned_light,
                        source_state,
                        settings,
                        |e| draw(path, e, 211),
                    )
                    .unwrap();
                near(result.factor_per_sr(), unpolarized.factor_per_sr(), 3e-11);
                assert_eq!(
                    (result.facet_events(), result.transmissions()),
                    (unpolarized.facet_events(), unpolarized.transmissions())
                );
            }
        }
    });
}

#[test]
fn matched_tail_work_is_jointly_bounded_replayable_and_cancellable() {
    with_cx(|cx, gate| {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let input = beam(1.0, normal.scale(-1.0), Stokes::unpolarized(1.0));
        let light = source_frame(Vec3::new(-0.9, 0.0, (1.0_f64 - 0.9 * 0.9).sqrt()));
        let surface = DielectricSurface::try_rough(1.0)
            .unwrap()
            .with_coating(coating(0.02));
        let mut calls = 0;
        let refused = input.evaluate_smith_dielectric_connection(
            cx,
            normal,
            glass(1.0),
            surface,
            light,
            Stokes::unpolarized(1.0),
            SmithWalkSettings::try_new(1).unwrap(),
            |ordinal| {
                assert_eq!(ordinal, calls);
                calls += 1;
                [0.5, 0.5, 0.25, 0.99]
            },
        );
        assert_eq!(
            refused,
            Err(SmithDielectricError::FacetLimit { facet_events: 1 })
        );
        assert_eq!(
            calls, 2,
            "the null-tail collision must consume the shared facet budget"
        );
        let run = || {
            input.evaluate_smith_dielectric_connection(
                cx,
                normal,
                glass(1.0),
                surface,
                light,
                Stokes::unpolarized(1.0),
                SmithWalkSettings::try_new(2048).unwrap(),
                |e| draw(13, e, 221),
            )
        };
        assert_eq!(run(), run());
        let cancelled = input.evaluate_smith_dielectric_connection(
            cx,
            normal,
            glass(1.0),
            surface,
            light,
            Stokes::unpolarized(1.0),
            SmithWalkSettings::try_new(2048).unwrap(),
            |ordinal| {
                if ordinal == 1 {
                    gate.request();
                }
                [0.5, 0.5, 0.25, 0.99]
            },
        );
        assert_eq!(cancelled, Err(SmithDielectricError::Cancelled));
    });
}
