//! Moving sources enter both native stages and the original contact/radiation balance.
use super::*;
use fs_conduction::COMPONENT_POWER_DIMS;
use fs_conduction::bc::history::march::enthalpy::scanning::{ScanRadiation, ScanningEnthalpyRun};
use fs_conduction::bc::history::march::moving::motion::ScanMotionLimit;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::interface::temperature::ContactIterationConfig;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Stage, SDIRK2_GAMMA as G};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn constant(value: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Constant(QtyAny::new(value, dims))
}
fn ramp(start: f64, speed: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Ramp {
        t_start: 0.0,
        t_end: 1.0,
        from: QtyAny::new(start, dims),
        to: QtyAny::new(start + speed, dims),
    }
}
fn moving(z: TimeSignal, power: TimeSignal) -> MovingHeatSource {
    MovingHeatSource::new(
        [constant(0.2, LENGTH), constant(0.2, LENGTH), z],
        [2.0; 3],
        power,
        0.0,
        MovingHeatBudget {
            max_vertices: 8,
            max_elements: 2,
            min_active_vertices: 8,
            power_relative_tolerance: 1e-12,
        },
    )
    .unwrap()
}

// Independent compact-footprint samples. Every vertex in these two unit
// tetrahedra has source-volume weight 1/24, regardless of its reference mass.
fn density(mesh: &ConductionMesh, heaters: &[(f64, f64)], background: f64) -> [f64; 8] {
    let mut q = [background; 8];
    for &(z, watts) in heaters {
        let weights: [f64; 8] = std::array::from_fn(|i| {
            let p = mesh.positions()[i];
            let r2 = ((p[0] - 0.2).powi(2) + (p[1] - 0.2).powi(2) + (p[2] - z).powi(2)) / 4.0;
            (1.0 - r2).max(0.0).powi(2)
        });
        let volume = weights.iter().sum::<f64>() / 24.0;
        for i in 0..8 {
            q[i] += watts * weights[i] / volume;
        }
    }
    q
}

// Exact consistent P1 source load, divided by the original nodal reference mass.
fn source_rates(q: &[f64; 8]) -> [f64; 8] {
    std::array::from_fn(|i| {
        let body = i / 4;
        (q[4 * body..4 * body + 4].iter().sum::<f64>() + q[i]) / if body == 0 { 5.0 } else { 10.0 }
    })
}

// R=1+2fa+fb and Q=50/R. At this latent plateau, source rates may be
// spatially nonuniform while the contact heat is uniform on the matching face.
fn contact_stage(old: &[f64; 8], dt: f64, q: &[f64; 8]) -> [f64; 8] {
    let rate = source_rates(q);
    let free: [f64; 8] = std::array::from_fn(|i| old[i] + dt * rate[i]);
    let a = 1.0 + free[..3].iter().sum::<f64>() / 150.0 + free[4..7].iter().sum::<f64>() / 300.0;
    let heat = 100.0 / (a + (a * a + 24.0 * dt).sqrt());
    std::array::from_fn(|i| {
        free[i]
            + dt * heat
                * if i % 4 == 3 {
                    0.0
                } else if i < 4 {
                    8.0
                } else {
                    -4.0
                }
    })
}

fn contact_stages(old: &[f64; 8], dt: f64, q: &[[f64; 8]; 2]) -> [[f64; 8]; 2] {
    let first = contact_stage(old, G * dt, &q[0]);
    let carry = std::array::from_fn(|i| old[i] + (1.0 - G) / G * (first[i] - old[i]));
    [first, contact_stage(&carry, G * dt, &q[1])]
}

fn admit(stage: &EnthalpySdirk2Stage, expected: &[f64; 8], watts: f64) {
    for i in 0..8 {
        near(stage.specific_enthalpy_j_kg[i], expected[i], 5e-8);
        near(stage.liquid_mass_fraction[i], expected[i] / 100.0, 5e-10);
    }
    assert_eq!(
        stage.temperature,
        [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]
    );
    near(stage.source_w, watts, 2e-12);
    assert!(stage.physical_contacts.is_some());
    assert!(stage.stage_residual_norm_j <= config().residual_tolerance_j);
    assert!(stage.stage_energy_residual_j.abs() <= config().step.energy_tolerance_j);
}

#[test]
fn moving_sdirk_samples_both_physical_stages_and_retains_original_phase_contact() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let heaters = vec![
            moving(
                ramp(-0.2, 4.0, LENGTH),
                ramp(1.0, 10.0, COMPONENT_POWER_DIMS),
            ),
            moving(constant(0.5, LENGTH), constant(0.5, COMPONENT_POWER_DIMS)),
        ];
        let mut run = ContactEnthalpyHistoryRun::new_newton(
            cx,
            base(
                cx,
                &boundary,
                &conductivity,
                &curves,
                &OLD,
                HeatSourceHistory::Fixed(ScalarField::Uniform(3.0)),
                2,
            ),
            contacts(
                cx,
                &mesh,
                &initial,
                &cards,
                &curves,
                &source,
                &context,
                sides(&mesh, 4, 2, 0.0),
                0,
                false,
            ),
            config(),
        )
        .unwrap()
        .into_sdirk2()
        .unwrap()
        .with_moving_sources(cx, heaters.clone(), 2)
        .unwrap();
        assert_eq!(run.moving_sources(), heaters);
        assert!(run.contacts().is_some());
        let dt = 0.04;
        let q = [G * dt, dt].map(|time| {
            density(
                &mesh,
                &[(-0.2 + 4.0 * time, 1.0 + 10.0 * time), (0.5, 0.5)],
                3.0,
            )
        });
        let expected = contact_stages(&OLD, dt, &q);
        let solved = run.advance_to(cx, dt, config().step).unwrap();
        admit(&solved.first_stage, &expected[0], 2.5 + 10.0 * G * dt);
        admit(&solved.second_stage, &expected[1], 2.5 + 10.0 * dt);
        let exact_energy = 2.5 * dt + 5.0 * dt * dt;
        near(solved.source_j, exact_energy, 2e-12);
        near(solved.stored_energy_change_j, exact_energy, 2e-9);
        near(
            solved.energy_residual_j,
            0.0,
            config().step.energy_tolerance_j,
        );
        assert_eq!(
            run.conduction().specific_enthalpy_j_kg(),
            solved.second_stage.specific_enthalpy_j_kg
        );
        assert_eq!(run.conduction().accepted_steps(), 1);
        assert_eq!(run.conduction().time_s(), dt);
    });
}

#[test]
fn moving_sdirk_resolves_held_path_and_power_events_with_left_endpoint_samples() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let heater = moving(
            TimeSignal::Table {
                times: vec![0.0, 0.02],
                values: vec![-0.2, 0.2],
                dims: LENGTH,
                interp: Interp::Hold,
            },
            TimeSignal::Table {
                times: vec![0.0, 0.01, 0.03],
                values: vec![0.0, 1.0, 0.0],
                dims: COMPONENT_POWER_DIMS,
                interp: Interp::Hold,
            },
        );
        let mut run = ContactEnthalpyHistoryRun::new_newton(
            cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 4),
            contacts(
                cx,
                &mesh,
                &initial,
                &cards,
                &curves,
                &source,
                &context,
                sides(&mesh, 4, 2, 0.0),
                0,
                false,
            ),
            config(),
        )
        .unwrap()
        .into_sdirk2()
        .unwrap()
        .with_moving_sources(cx, vec![heater], 1)
        .unwrap();
        let mut old = OLD;
        let mut from = 0.0;
        let mut energy = 0.0;
        for (to, z, watts) in [
            (0.01, -0.2, 0.0),
            (0.02, -0.2, 1.0),
            (0.03, 0.2, 1.0),
            (0.04, 0.2, 0.0),
        ] {
            let q = density(&mesh, &[(z, watts)], 0.0);
            let expected = contact_stages(&old, to - from, &[q; 2]);
            let accepted = run
                .advance_event_step(cx, 0.04, 1.0, config().step)
                .unwrap();
            assert_eq!(accepted.end_time_s, to);
            assert_eq!(accepted.start_time_s, from);
            admit(&accepted.report.first_stage, &expected[0], watts);
            admit(&accepted.report.second_stage, &expected[1], watts);
            energy += accepted.report.source_j;
            old = expected[1];
            from = to;
        }
        near(energy, 0.02, 2e-12);
        assert_eq!(run.conduction().accepted_steps(), 4);
    });
}

#[test]
fn moving_sdirk_motion_cap_controls_both_adaptive_halves_and_entropy_path() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let create = || {
            ContactEnthalpyHistoryRun::new_newton(
                cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
                contacts(
                    cx,
                    &mesh,
                    &initial,
                    &cards,
                    &curves,
                    &source,
                    &context,
                    sides(&mesh, 4, 2, 0.0),
                    0,
                    false,
                ),
                config(),
            )
            .unwrap()
            .into_sdirk2()
            .unwrap()
            .with_moving_sources(
                cx,
                vec![
                    moving(
                        ramp(-0.2, 4.0, LENGTH),
                        ramp(1.0, 10.0, COMPONENT_POWER_DIMS),
                    ),
                    moving(ramp(-0.2, 8.0, LENGTH), constant(0.0, COMPONENT_POWER_DIMS)),
                ],
                2,
            )
            .unwrap()
        };
        let policy = ScanMotionLimit {
            max_normalized_travel: 0.02,
            max_chebyshev_coefficients: 0,
        };
        let adaptive = AdaptiveHistoryConfig {
            absolute_temperature_k: 1e-9,
            relative_temperature: 0.0,
            min_step_s: 0.02,
            max_step_s: 0.1,
            max_attempts: 8,
        };
        let phase = PhaseErrorTolerance {
            absolute_enthalpy_j_kg: 100.0,
            relative_enthalpy: 0.0,
            absolute_liquid_fraction: 1.0,
        };
        for entropy in [false, true] {
            let mut run = create();
            run.set_motion_limit(cx, policy).unwrap();
            assert!(run.advance_to(cx, 0.04, config().step).is_err());
            assert_eq!(run.conduction().specific_enthalpy_j_kg(), &OLD);
            assert_eq!(run.conduction().time_s(), 0.0);
            let accepted = if entropy {
                run.advance_adaptive_step_with_entropy(
                    cx,
                    0.04,
                    0.1,
                    config().step,
                    adaptive,
                    phase,
                    1e-8,
                )
                .unwrap()
                .physical
            } else {
                run.advance_adaptive_step(cx, 0.04, 0.1, config().step, adaptive, phase)
                    .unwrap()
            };
            near(accepted.end_time_s, 0.005, 1e-14);
            assert_eq!(accepted.rejected_attempts, 0);
            let mut old = OLD;
            let mut from = 0.0;
            for (to, report) in [
                (accepted.midpoint_time_s, &accepted.first_half),
                (accepted.end_time_s, &accepted.second_half),
            ] {
                let dt = to - from;
                let times = [from + G * dt, to];
                let q = times
                    .map(|time| density(&mesh, &[(-0.2 + 4.0 * time, 1.0 + 10.0 * time)], 0.0));
                let expected = contact_stages(&old, dt, &q);
                admit(&report.first_stage, &expected[0], 1.0 + 10.0 * times[0]);
                admit(&report.second_stage, &expected[1], 1.0 + 10.0 * times[1]);
                old = expected[1];
                from = to;
            }
            assert_eq!(run.conduction().accepted_steps(), 2);
            assert!(run.set_motion_limit(cx, policy).is_err());
            assert_eq!(run.motion_limit(), Some(policy));
        }
    });
}

#[test]
fn moving_sdirk_future_source_and_cancellation_failures_preserve_resume() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let create = || {
            ContactEnthalpyHistoryRun::new_newton(
                cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
                contacts(
                    cx,
                    &mesh,
                    &initial,
                    &cards,
                    &curves,
                    &source,
                    &context,
                    sides(&mesh, 4, 2, 0.0),
                    0,
                    false,
                ),
                config(),
            )
            .unwrap()
            .into_sdirk2()
            .unwrap()
            .with_moving_sources(
                cx,
                vec![moving(
                    constant(0.0, LENGTH),
                    TimeSignal::Table {
                        times: vec![0.0, 0.01],
                        values: vec![1.0, -1.0],
                        dims: COMPONENT_POWER_DIMS,
                        interp: Interp::Hold,
                    },
                )],
                1,
            )
            .unwrap()
        };
        let mut run = create();
        assert!(run.advance_to(cx, 0.02, config().step).is_err());
        support::with_cancelled_cx(|cancelled| {
            assert!(run.advance_to(cancelled, 0.005, config().step).is_err())
        });
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), &OLD);
        assert_eq!(run.conduction().time_s(), 0.0);
        assert_eq!(run.conduction().accepted_steps(), 0);
        let accepted = run.advance_to(cx, 0.005, config().step).unwrap();
        let replay = create().advance_to(cx, 0.005, config().step).unwrap();
        assert_eq!(
            accepted.second_stage.specific_enthalpy_j_kg,
            replay.second_stage.specific_enthalpy_j_kg
        );
        assert!(run.advance_to(cx, 0.02, config().step).is_err());
        assert_eq!(
            run.conduction().specific_enthalpy_j_kg(),
            accepted.second_stage.specific_enthalpy_j_kg
        );
        assert_eq!(run.conduction().time_s(), 0.005);
        assert_eq!(run.conduction().accepted_steps(), 1);
    });
}

fn radiative_stage(old: &[f64; 8], dt: f64, q: &[f64; 8], bath: f64) -> [f64; 8] {
    let source = source_rates(q);
    let mut h = *old;
    for _ in 0..128 {
        let (outward, _) = flux(&h, 0.4, bath);
        h = std::array::from_fn(|i| {
            old[i] + dt * (source[i] - outward[i] * if i < 4 { 24.0 } else { 12.0 })
        });
    }
    assert!(h.iter().all(|h| *h > 0.0 && *h < 100.0));
    h
}

#[test]
fn moving_sdirk_retains_phase_contacts_and_material_reflection_at_each_stage() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let enclosure = open(cx, &mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new_open(
            cx,
            &mesh,
            &initial,
            &enclosure,
            350.0,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let dt = 0.0004;
        let bath = TimeSignal::Ramp {
            t_start: 0.0,
            t_end: dt,
            from: QtyAny::new(350.0, TEMPERATURE_DIMS),
            to: QtyAny::new(450.0, TEMPERATURE_DIMS),
        };
        let mut run = ContactEnthalpyHistoryRun::new_newton(
            cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
            contacts(
                cx,
                &mesh,
                &initial,
                &cards,
                &curves,
                &source,
                &context,
                sides(&mesh, 4, 2, 0.0),
                0,
                false,
            ),
            config(),
        )
        .unwrap()
        .with_material_enclosure(
            cx,
            &radiation,
            patches(&cards, &curves, &context),
            Some(bath),
        )
        .unwrap()
        .into_sdirk2()
        .unwrap()
        .with_moving_sources(
            cx,
            vec![moving(
                ramp(-0.2, 500.0, LENGTH),
                ramp(10.0, 5000.0, COMPONENT_POWER_DIMS),
            )],
            1,
        )
        .unwrap();
        let times = [G * dt, dt];
        let baths = [350.0 + 100.0 * G, 450.0];
        let q =
            times.map(|time| density(&mesh, &[(-0.2 + 500.0 * time, 10.0 + 5000.0 * time)], 0.0));
        let first = radiative_stage(&OLD, G * dt, &q[0], baths[0]);
        let carry = std::array::from_fn(|i| OLD[i] + (1.0 - G) / G * (first[i] - OLD[i]));
        let expected = [first, radiative_stage(&carry, G * dt, &q[1], baths[1])];
        let solved = run.advance_to(cx, dt, config().step).unwrap();
        for (index, stage) in [&solved.first_stage, &solved.second_stage]
            .into_iter()
            .enumerate()
        {
            admit(stage, &expected[index], 10.0 + 5000.0 * times[index]);
            let report = stage.pointwise_enclosure_radiation.as_ref().unwrap();
            near(
                report.ambient.as_ref().unwrap().temperature_k,
                baths[index],
                1e-12,
            );
            let (_, heat) = flux(&expected[index], 0.4, baths[index]);
            for (body, patch) in report.patches.iter().enumerate() {
                near(patch.outward_heat_w, heat[body], 3e-7);
                for state in &patch.material_states {
                    near(
                        state.emissivity.value(),
                        0.2 + 0.006 * expected[index][state.vertex],
                        4e-10,
                    );
                }
            }
        }
        let exact_source = 10.0 * dt + 2500.0 * dt * dt;
        near(solved.source_j, exact_source, 2e-12);
        near(
            solved.stored_energy_change_j + solved.radiation_out_j,
            exact_source,
            2e-9,
        );
        assert!(run.contacts().is_some());
        assert_eq!(run.material_enclosure_patches().unwrap().len(), 2);
        assert!(solved.radiation_out_j.abs() > 1e-4);
    });
}

#[test]
fn moving_sdirk_conversion_preserves_accepted_scan_state_sources_and_phase_contacts() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let heaters = vec![moving(
            ramp(-0.2, 4.0, LENGTH),
            ramp(1.0, 10.0, COMPONENT_POWER_DIMS),
        )];
        let create = |radiation| {
            ScanningEnthalpyRun::new(
                cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
                heaters.clone(),
                1,
                radiation,
            )
            .unwrap()
        };
        let mut scan = create(None)
            .with_contacts_newton(
                cx,
                contacts(
                    cx,
                    &mesh,
                    &initial,
                    &cards,
                    &curves,
                    &source,
                    &context,
                    sides(&mesh, 4, 2, 0.0),
                    0,
                    false,
                ),
                config(),
            )
            .unwrap();
        let policy = ScanMotionLimit {
            max_normalized_travel: 0.02,
            max_chebyshev_coefficients: 0,
        };
        scan.set_motion_limit(cx, policy).unwrap();
        let from = 0.005;
        let accepted = scan.advance_to(cx, from, config().step).unwrap();
        let expected_be = contact_stage(
            &OLD,
            from,
            &density(&mesh, &[(-0.2 + 4.0 * from, 1.0 + 10.0 * from)], 0.0),
        );
        let old: [f64; 8] = scan
            .conduction()
            .specific_enthalpy_j_kg()
            .try_into()
            .unwrap();
        for i in 0..8 {
            near(old[i], expected_be[i], 5e-8);
        }
        let mut continued = scan.into_sdirk2().unwrap();
        assert_eq!(continued.conduction().specific_enthalpy_j_kg(), old);
        assert_eq!(
            continued.conduction().temperature(),
            accepted.solution().temperature
        );
        assert_eq!(
            continued.conduction().liquid_mass_fraction(),
            accepted.solution().liquid_mass_fraction
        );
        assert_eq!(continued.conduction().time_s(), from);
        assert_eq!(continued.conduction().accepted_steps(), 1);
        assert_eq!(continued.conduction().max_steps(), 2);
        assert_eq!(continued.moving_sources(), heaters);
        assert_eq!(continued.motion_limit(), Some(policy));
        assert!(continued.contacts().is_some());
        assert!(std::ptr::eq(
            continued.conduction().storage().phase_curve_at(0).unwrap(),
            &curves[0]
        ));
        assert!(continued.advance_to(cx, 0.025, config().step).is_err());
        let next = continued
            .advance_event_step(cx, 0.025, 1.0, config().step)
            .unwrap();
        let dt = next.end_time_s - from;
        near(dt, 0.01, 1e-14);
        let times = [from + G * dt, next.end_time_s];
        let q = times.map(|time| density(&mesh, &[(-0.2 + 4.0 * time, 1.0 + 10.0 * time)], 0.0));
        let expected = contact_stages(&old, dt, &q);
        admit(
            &next.report.first_stage,
            &expected[0],
            1.0 + 10.0 * times[0],
        );
        admit(
            &next.report.second_stage,
            &expected[1],
            1.0 + 10.0 * times[1],
        );
        near(
            next.report.source_j,
            dt * (1.0 + 5.0 * (from + next.end_time_s)),
            2e-12,
        );
        assert_eq!(continued.conduction().accepted_steps(), 2);
        assert!(
            continued
                .advance_event_step(cx, 0.025, 1.0, config().step)
                .is_err()
        );

        let emissivity = SurfaceEmissivity::from_condition(
            REGIONS[0],
            material(&cards[0], &context),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0;
        let radiation = ScanRadiation {
            patches: vec![RadiationRegionHistory {
                region: REGIONS[0].into(),
                emissivity,
                ambient_temperature: bath(350.0),
            }],
            config: AmbientRadiationConfig::default(),
        };
        assert!(create(Some(radiation)).into_sdirk2().is_err());
        let law = TemperatureContact::new_two_sided(
            "joint",
            ThermalInterfaces::coincident_face_pairs(&mesh).unwrap(),
            &source,
            &context,
            ["aT", "bT"],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let legacy =
            TemperatureContacts::new(cx, &mesh, &initial, Vec::new(), vec![law], 1).unwrap();
        assert!(
            create(None)
                .with_contacts(
                    cx,
                    legacy,
                    ContactIterationConfig {
                        max_iterations: 64,
                        relaxation: 0.5,
                        residual_tolerance_j: 1e-9
                    },
                    1
                )
                .unwrap()
                .into_sdirk2()
                .is_err()
        );
    });
}
