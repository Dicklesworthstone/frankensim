//! Independent two-stage equations, finite-bath refinement and original phases.
use super::*;
use fs_couple::moisture::layered::MoistureLayeredSdirk2Report;
use fs_material::osmosis::solubility::FiniteSolubilitySolution;

#[path = "sdirk/adaptive.rs"]
mod adaptive;

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

struct Stage {
    fields: Vec<[f64; 2]>,
    activity: f64,
    volume: f64,
}

struct Reference {
    first: Stage,
    last: Stage,
}

// Solve the independent two-store physical equations in activity coordinates.
// The right-hand side uses actual original water storage plus known kg loads.
fn linear_stage(
    p: &Physics,
    initial: [f64; 2],
    load: [f64; 2],
    dt: f64,
    activity: f64,
) -> [f64; 2] {
    let a = p.capacity[0] + dt * p.join;
    let b = dt * p.join;
    let d = p.capacity[1] + dt * (p.join + p.outer);
    let rhs = [
        p.capacity[0] * initial[0] + load[0],
        p.capacity[1] * initial[1] + load[1] + dt * p.outer * activity,
    ];
    let determinant = a * d - b * b;
    [
        (d * rhs[0] + b * rhs[1]) / determinant,
        (b * rhs[0] + a * rhs[1]) / determinant,
    ]
}

fn water(p: &Physics, initial: [f64; 2], next: [f64; 2]) -> f64 {
    (0..2).map(|i| p.capacity[i] * (next[i] - initial[i])).sum()
}

fn face_load(p: &Physics, a: [f64; 2], bath: f64, weighted_dt: f64) -> [f64; 2] {
    let internal = p.join * (a[1] - a[0]);
    let outer = p.outer * (bath - a[1]);
    [weighted_dt * internal, weighted_dt * (outer - internal)]
}

// A separate scalar activity equation eliminates solvent volume using the
// complete closed water inventory. No native field/root method is invoked.
fn finite_stage(
    physics: &[Physics],
    initial: &[[f64; 2]],
    known: &[[f64; 2]],
    dt: f64,
    volume: f64,
    solute: f64,
) -> Stage {
    let evaluate = |activity| {
        let fields: Vec<_> = physics
            .iter()
            .enumerate()
            .map(|(i, p)| linear_stage(p, initial[i], known[i], dt, activity))
            .collect();
        let stored: f64 = physics
            .iter()
            .enumerate()
            .map(|(i, p)| water(p, initial[i], fields[i]))
            .sum();
        let volume = volume - stored / WATER_DENSITY;
        assert!(volume > 0.0);
        let residual = activity - (-aqueous::WATER_VOLUME * solute / volume).exp();
        (
            Stage {
                fields,
                activity,
                volume,
            },
            residual,
        )
    };
    let (mut lo, mut hi) = (0.0, 1.0);
    assert!(evaluate(lo).1 < 0.0 && evaluate(hi).1 > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if evaluate(mid).1 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    evaluate(0.5 * (lo + hi)).0
}

fn reference(
    physics: &[Physics],
    initial: &[[f64; 2]],
    dt: f64,
    volume: f64,
    solute: f64,
) -> Reference {
    let first = finite_stage(
        physics,
        initial,
        &vec![[0.0; 2]; physics.len()],
        GAMMA * dt,
        volume,
        solute,
    );
    let known: Vec<_> = physics
        .iter()
        .enumerate()
        .map(|(i, p)| face_load(p, first.fields[i], first.activity, (1.0 - GAMMA) * dt))
        .collect();
    let last = finite_stage(physics, initial, &known, GAMMA * dt, volume, solute);
    Reference { first, last }
}

fn check_transfers(
    report: &MoistureLayeredSdirk2Report,
    p: &Physics,
    initial: [f64; 2],
    first: [f64; 2],
    last: [f64; 2],
    activity: [f64; 2],
    dt: f64,
) {
    let join = p.join * dt * ((1.0 - GAMMA) * (first[0] - first[1]) + GAMMA * (last[0] - last[1]));
    let outer =
        p.outer * dt * ((1.0 - GAMMA) * (activity[0] - first[1]) + GAMMA * (activity[1] - last[1]));
    near(report.endpoint.boundary_water_kg[1], outer, 2e-12);
    near(
        report.endpoint.interfaces[0].water_left_to_right_kg,
        join,
        2e-12,
    );
    near(p.capacity[0] * (last[0] - initial[0]), -join, 2e-12);
    near(p.capacity[1] * (last[1] - initial[1]), join + outer, 2e-12);
    near(report.endpoint.water_balance_residual_kg, 0.0, 2e-12);
    assert_eq!(report.endpoint.boundary_water_kg[0], 0.0);
    assert_eq!(report.endpoint.interval_s, [0.0, dt]);
    assert_eq!(report.first_stage.interval_s, [0.0, GAMMA * dt]);
    near(
        report.first_stage.boundary_water_kg[1],
        GAMMA * dt * p.outer * (activity[0] - first[1]),
        2e-12,
    );
    near(
        report.first_stage.interfaces[0].water_left_to_right_kg,
        GAMMA * dt * p.join * (first[0] - first[1]),
        2e-12,
    );
    for i in 0..2 {
        near(
            report.first_stage.transport_end[i].query().axes()["moisture_content"],
            SLOPE[i] * first[i],
            2e-11,
        );
        near(
            report.endpoint.transport_end[i].query().axes()["moisture_content"],
            SLOPE[i] * last[i],
            2e-11,
        );
    }
    // Interface/surface diagnostics describe final instantaneous faces, whereas
    // their integrated water accounts contain BOTH stage contributions.
    let face = &report.endpoint.interfaces[0];
    let rate = p.join * (last[0] - last[1]);
    near(face.relative_activity, last[0] - rate / p.half[0], 2e-10);
    near(
        face.right_relative_activity,
        last[1] + rate / p.half[1],
        2e-10,
    );
    near(
        report.endpoint.surface_activities[1].unwrap(),
        last[1] + p.outer * (activity[1] - last[1]) / p.outer_bulk,
        2e-10,
    );
}

#[test]
fn both_layered_geometries_use_two_distinct_finite_bath_stages_and_weighted_faces() {
    with_cx(|cx| {
        let volume = 1e-5;
        let solute = 0.02;
        let dt = 2.0;
        for symmetry in [
            None,
            Some(RadialMoistureSymmetry::Sphere),
            Some(RadialMoistureSymmetry::Cylinder { length_m: 0.1 }),
        ] {
            let p = coefficients(symmetry);
            for initial in [[0.1, 0.2], [0.98, 0.99]] {
                let expected = reference(std::slice::from_ref(&p), &[initial], dt, volume, solute);
                let bath = aqueous::reservoir(volume, vec![solute]);
                let saved = bath.solution().clone();
                let (next_bath, activity, time, mass, report) = if let Some(symmetry) = symmetry {
                    let original = radial_body(symmetry, initial, cx);
                    let (bath, next, report) = bath
                        .advance_layered_radial_sdirk2(
                            &original,
                            dt,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    for (old, new) in original.layers().iter().zip(next.layers()) {
                        assert_eq!(old.radial().dry_masses_kg(), new.radial().dry_masses_kg());
                        assert_eq!(old.sorption().identity(), new.sorption().identity());
                    }
                    assert_eq!(original.time_s(), 0.0);
                    (
                        bath,
                        next.activity().to_vec(),
                        next.time_s(),
                        next.water_mass_kg(),
                        report,
                    )
                } else {
                    let original = slab_body(initial, cx);
                    let (bath, next, report) = bath
                        .advance_layered_sdirk2(
                            &original,
                            dt,
                            [false, true],
                            config(),
                            nonlinear(),
                            aqueous::controls(volume),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    for (old, new) in original.layers().iter().zip(next.layers()) {
                        assert_eq!(old.slab().dry_mass_kg(), new.slab().dry_mass_kg());
                        assert_eq!(old.sorption().identity(), new.sorption().identity());
                    }
                    assert_eq!(original.time_s(), 0.0);
                    (
                        bath,
                        next.activity().to_vec(),
                        next.time_s(),
                        next.water_mass_kg(),
                        report,
                    )
                };
                for i in 0..2 {
                    near(activity[i], expected.last.fields[0][i], 2e-10);
                }
                check_transfers(
                    &report.transport,
                    &p,
                    initial,
                    expected.first.fields[0],
                    expected.last.fields[0],
                    [expected.first.activity, expected.last.activity],
                    dt,
                );
                assert!((expected.first.activity - expected.last.activity).abs() > 1e-7);
                near(report.water_activity[1], expected.last.activity, 2e-10);
                near(
                    next_bath.solution().solvent_volume_m3(),
                    expected.last.volume,
                    2e-16,
                );
                near(
                    report.water_to_solid_kg,
                    water(&p, initial, expected.last.fields[0]),
                    2e-12,
                );
                near(report.water_balance_residual_kg, 0.0, 1e-13);
                let initial_water: f64 = (0..2).map(|i| p.capacity[i] * initial[i]).sum();
                near(
                    next_bath.solvent_mass_kg() + mass,
                    bath.solvent_mass_kg() + initial_water,
                    2e-12,
                );
                assert_eq!(next_bath.solution().solute_moles(), &[solute]);
                assert_eq!(bath.solution(), &saved);
                assert_eq!(time, dt);
            }
        }
    });
}

// Independent continuous two-store ODE, with bath volume eliminated by water
// conservation. RK4 is only a test oracle, never a production adapter.
fn continuous(
    p: &Physics,
    initial: [f64; 2],
    duration: f64,
    volume: f64,
    solute: f64,
    steps: usize,
) -> [f64; 2] {
    let rhs = |state: [f64; 2]| {
        let v = volume - water(p, initial, state) / WATER_DENSITY;
        assert!(v > 0.0);
        let bath = (-aqueous::WATER_VOLUME * solute / v).exp();
        let load = face_load(p, state, bath, 1.0);
        std::array::from_fn::<_, 2, _>(|i| load[i] / p.capacity[i])
    };
    let h = duration / steps as f64;
    let mut state = initial;
    for _ in 0..steps {
        let a = rhs(state);
        let b = rhs(std::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| state[i] + h * c[i]));
        for i in 0..2 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn finite_bath_layered_sdirk_refines_at_second_order_for_slabs_spheres_and_cylinders() {
    with_cx(|cx| {
        let initial = [0.1, 0.2];
        let volume = 1e-4;
        let solute = 0.2;
        let duration = 8.0;
        for symmetry in [
            None,
            Some(RadialMoistureSymmetry::Sphere),
            Some(RadialMoistureSymmetry::Cylinder { length_m: 0.1 }),
        ] {
            let p = coefficients(symmetry);
            let exact = continuous(&p, initial, duration, volume, solute, 32768);
            let finer = continuous(&p, initial, duration, volume, solute, 65536);
            for i in 0..2 {
                near(exact[i], finer[i], 2e-12);
            }
            let mut errors = Vec::new();
            for intervals in [4, 8, 16] {
                let dt = duration / intervals as f64;
                let mut bath = aqueous::reservoir(volume, vec![solute]);
                let final_activity = if let Some(symmetry) = symmetry {
                    let mut field = radial_body(symmetry, initial, cx);
                    let total = field.water_mass_kg() + bath.solvent_mass_kg();
                    for _ in 0..intervals {
                        let physical = aqueous::controls(bath.solution().solvent_volume_m3());
                        let (next_bath, next, step) = bath
                            .advance_layered_radial_sdirk2(
                                &field,
                                dt,
                                [false, true],
                                config(),
                                nonlinear(),
                                physical,
                                cx,
                            )
                            .unwrap()
                            .into_parts();
                        near(step.water_balance_residual_kg, 0.0, 1e-13);
                        near(
                            next.water_mass_kg() + next_bath.solvent_mass_kg(),
                            total,
                            5e-12,
                        );
                        bath = next_bath;
                        field = next;
                    }
                    assert_eq!(field.time_s(), duration);
                    field.activity().to_vec()
                } else {
                    let mut field = slab_body(initial, cx);
                    let total = field.water_mass_kg() + bath.solvent_mass_kg();
                    for _ in 0..intervals {
                        let physical = aqueous::controls(bath.solution().solvent_volume_m3());
                        let (next_bath, next, step) = bath
                            .advance_layered_sdirk2(
                                &field,
                                dt,
                                [false, true],
                                config(),
                                nonlinear(),
                                physical,
                                cx,
                            )
                            .unwrap()
                            .into_parts();
                        near(step.water_balance_residual_kg, 0.0, 1e-13);
                        near(
                            next.water_mass_kg() + next_bath.solvent_mass_kg(),
                            total,
                            5e-12,
                        );
                        bath = next_bath;
                        field = next;
                    }
                    assert_eq!(field.time_s(), duration);
                    field.activity().to_vec()
                };
                errors.push(
                    (0..2)
                        .map(|i| (final_activity[i] - exact[i]).abs())
                        .fold(0.0_f64, f64::max),
                );
            }
            for pair in errors.windows(2) {
                assert!(
                    pair[0] / pair[1] > 3.7 && pair[0] / pair[1] < 4.4,
                    "{symmetry:?}: {errors:?}"
                );
            }
        }
    });
}

#[test]
fn shared_layered_sdirk_retains_opposing_fluxes_and_canonical_ids_at_both_stages() {
    with_cx(|cx| {
        let volume = 1e-5;
        let solute = 0.02;
        let dt = 2.0;
        let initial = [[0.1, 0.2], [0.98, 0.99]];
        let bath = aqueous::reservoir(volume, vec![solute]);
        let slabs = initial.map(|a| slab_body(a, cx));
        let p = [coefficients(None), coefficients(None)];
        let expected = reference(&p, &initial, dt, volume, solute);
        let run = |reverse| {
            let bodies = if reverse {
                [
                    slab_participant(23, &slabs[1]),
                    slab_participant(7, &slabs[0]),
                ]
            } else {
                [
                    slab_participant(7, &slabs[0]),
                    slab_participant(23, &slabs[1]),
                ]
            };
            bath.advance_shared_layered_sdirk2(&bodies, dt, 4, aqueous::controls(volume), cx)
                .unwrap()
        };
        let result = run(false);
        let reverse = run(true);
        assert_eq!(
            result.reservoir().solution(),
            reverse.reservoir().solution()
        );
        assert_eq!(
            result.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![7, 23]
        );
        for (i, (solid, reversed)) in result.solids().iter().zip(reverse.solids()).enumerate() {
            assert_eq!(solid.moisture.activity(), reversed.moisture.activity());
            for j in 0..2 {
                near(
                    solid.moisture.activity()[j],
                    expected.last.fields[i][j],
                    2e-10,
                );
            }
            check_transfers(
                &solid.transport,
                &p[i],
                initial[i],
                expected.first.fields[i],
                expected.last.fields[i],
                [expected.first.activity, expected.last.activity],
                dt,
            );
            assert_eq!(solid.transport.endpoint.boundary_water_kg[1] > 0.0, i == 0);
        }
        near(
            result.reservoir().solution().solvent_volume_m3(),
            expected.last.volume,
            2e-16,
        );
        near(
            result.step().water_to_solids_kg,
            result
                .solids()
                .iter()
                .map(|s| s.transport.endpoint.boundary_water_kg[1])
                .sum(),
            2e-12,
        );

        // Different curved geometry and storage may share the same physical bath.
        let symmetries = [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ];
        let radial = std::array::from_fn::<_, 2, _>(|i| radial_body(symmetries[i], initial[i], cx));
        let p = symmetries.map(|s| coefficients(Some(s)));
        let expected = reference(&p, &initial, dt, volume, solute);
        let run = |reverse| {
            let bodies = if reverse {
                [
                    radial_participant(23, &radial[1]),
                    radial_participant(7, &radial[0]),
                ]
            } else {
                [
                    radial_participant(7, &radial[0]),
                    radial_participant(23, &radial[1]),
                ]
            };
            bath.advance_shared_layered_radial_sdirk2(&bodies, dt, 4, aqueous::controls(volume), cx)
                .unwrap()
        };
        let result = run(false);
        let reverse = run(true);
        assert_eq!(
            result.reservoir().solution(),
            reverse.reservoir().solution()
        );
        assert_eq!(
            result.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![7, 23]
        );
        for (i, (solid, reversed)) in result.solids().iter().zip(reverse.solids()).enumerate() {
            assert_eq!(solid.moisture.activity(), reversed.moisture.activity());
            for j in 0..2 {
                near(
                    solid.moisture.activity()[j],
                    expected.last.fields[i][j],
                    2e-10,
                );
            }
            check_transfers(
                &solid.transport,
                &p[i],
                initial[i],
                expected.first.fields[i],
                expected.last.fields[i],
                [expected.first.activity, expected.last.activity],
                dt,
            );
            assert_eq!(solid.transport.endpoint.boundary_water_kg[1] > 0.0, i == 0);
        }
        near(
            result.reservoir().solution().solvent_volume_m3(),
            expected.last.volume,
            2e-16,
        );
        near(
            result.step().water_to_solids_kg,
            result
                .solids()
                .iter()
                .map(|s| s.transport.endpoint.boundary_water_kg[1])
                .sum(),
            2e-12,
        );
        for body in slabs {
            assert_eq!(body.time_s(), 0.0);
        }
        for body in radial {
            assert_eq!(body.time_s(), 0.0);
        }
        assert_eq!(bath.solution().solvent_volume_m3(), volume);
    });
}

#[test]
fn buffered_algebraic_phases_follow_each_sdirk_volume_without_losing_solute() {
    with_cx(|cx| {
        let volume = 1e-5;
        let saturation = 1000.0;
        let total = 0.1;
        let dt = 2.0;
        let law = IdealSolubilityV1::new(
            SpeciesId::new("water").unwrap(),
            SpeciesId::new("authored-layered-buffer").unwrap(),
            Temperature::new(293.15),
            AmountConcentration::new(saturation),
            0.0,
            [Temperature::new(293.15), Temperature::new(293.15)],
        )
        .unwrap();
        let bath = FiniteMoistureSolution::with_solubility(
            FiniteSolubilitySolution::new(
                Temperature::new(293.15),
                aqueous::WATER_VOLUME,
                Volume::new(volume),
                vec![(law, Amount::new(total))],
            )
            .unwrap(),
            0.018,
        )
        .unwrap();
        let initial = [0.1, 0.2];
        let original = radial_body(RadialMoistureSymmetry::Sphere, initial, cx);
        let p = coefficients(Some(RadialMoistureSymmetry::Sphere));
        let activity = (-aqueous::WATER_VOLUME * saturation).exp();
        let first = linear_stage(&p, initial, [0.0; 2], GAMMA * dt, activity);
        let last = linear_stage(
            &p,
            initial,
            face_load(&p, first, activity, (1.0 - GAMMA) * dt),
            GAMMA * dt,
            activity,
        );
        let final_volume = volume - water(&p, initial, last) / WATER_DENSITY;
        let result = bath
            .advance_layered_radial_sdirk2(
                &original,
                dt,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(volume),
                cx,
            )
            .unwrap();
        for i in 0..2 {
            near(result.moisture().activity()[i], last[i], 2e-10);
        }
        check_transfers(
            &result.step().transport,
            &p,
            initial,
            first,
            last,
            [activity; 2],
            dt,
        );
        for a in result.step().water_activity {
            near(a, activity, 2e-14);
        }
        near(
            result.reservoir().solution().solvent_volume_m3(),
            final_volume,
            2e-16,
        );
        let phase = result.reservoir().solubility().unwrap().phase_partitions()[0];
        near(phase.dissolved.value(), saturation * final_volume, 2e-13);
        near(
            phase.solid.value(),
            total - saturation * final_volume,
            2e-13,
        );
        near(phase.dissolved.value() + phase.solid.value(), total, 1e-16);
        assert!(
            phase.solid.value()
                > bath.solubility().unwrap().phase_partitions()[0]
                    .solid
                    .value()
        );
        near(bath.solution().solvent_volume_m3(), volume, 0.0);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn layered_sdirk_refuses_kinetic_stages_and_late_inadmissible_activity_atomically() {
    with_cx(|cx| {
        let kinetic = kinetic_bath();
        let slab = slab_body([0.1, 0.2], cx);
        let radial = radial_body(RadialMoistureSymmetry::Sphere, [0.8, 0.1], cx);
        assert!(matches!(
            kinetic.advance_layered_sdirk2(
                &slab,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(5e-5),
                cx,
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        assert!(matches!(
            kinetic.advance_layered_radial_sdirk2(
                &radial,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(5e-5),
                cx,
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        near(
            kinetic.kinetic_solubility().unwrap().last_steps()[0]
                .state()
                .dissolved()
                .value(),
            0.0,
            0.0,
        );
        let volume = 1e-4;
        let bath = aqueous::reservoir(volume, vec![0.2]);
        let saved = bath.solution().clone();
        let p = coefficients(Some(RadialMoistureSymmetry::Sphere));
        let rate = p.join * (1.0 / p.capacity[0] + 1.0 / p.capacity[1]);
        let dt = 10.0 / rate;
        // The first physical stage succeeds. The complete final stage would
        // drive core activity below zero and must refuse, without a BE fallback.
        bath.advance_layered_radial(
            &radial,
            GAMMA * dt,
            [false; 2],
            config(),
            nonlinear(),
            aqueous::controls(volume),
            cx,
        )
        .unwrap();
        assert!(matches!(
            bath.advance_layered_radial_sdirk2(
                &radial,
                dt,
                [false; 2],
                config(),
                nonlinear(),
                aqueous::controls(volume),
                cx,
            ),
            Err(SolutionMoistureError::Moisture(
                MoistureError::ActivityStageRange { .. }
            ))
        ));
        assert!(matches!(
            bath.advance_layered_radial_sdirk2(
                &radial,
                2.0,
                [true, false],
                config(),
                nonlinear(),
                aqueous::controls(volume),
                cx,
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        assert!(matches!(
            bath.advance_shared_layered_sdirk2(
                &[slab_participant(7, &slab), slab_participant(7, &slab)],
                2.0,
                4,
                aqueous::controls(volume),
                cx,
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        assert!(matches!(
            bath.advance_shared_layered_radial_sdirk2(
                &[
                    radial_participant(7, &radial),
                    radial_participant(9, &radial)
                ],
                2.0,
                3,
                aqueous::controls(volume),
                cx,
            ),
            Err(SolutionMoistureError::Invalid(_))
        ));
        let pressure = pressure_stack(cx);
        assert!(matches!(
            bath.advance_layered_sdirk2(
                &pressure,
                2.0,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(volume),
                cx,
            ),
            Err(SolutionMoistureError::Moisture(_))
        ));
        assert_eq!(bath.solution(), &saved);
        assert_eq!(slab.time_s(), 0.0);
        assert_eq!(radial.time_s(), 0.0);
        near(radial.activity()[0], 0.8, 1e-15);
        near(radial.activity()[1], 0.1, 1e-15);
    });
}
