//! Independent finite-resistance balances, nonlinear dry fronts and real vapor.
use super::super::atmosphere as air;
use super::*;
use MoistureActivityTransfer::{ConductanceKgM2S as Finite, IdealContact, Sealed};

const CAPACITY: [f64; 2] = [0.032, 0.012];
const HALF: [f64; 2] = [0.016, 0.048];

fn stack(initial: [f64; 2], degenerate: bool, cx: &Cx<'_>) -> MoistureLayeredSlab {
    MoistureLayeredSlab::new(
        vec![
            layer(
                material(
                    800.0,
                    vec![(0.0, 0.0), (1.0, 0.2)],
                    if degenerate { [0.0, 0.06] } else { [0.01; 2] },
                ),
                0.2,
                vec![0.2 * initial[0]],
                cx,
            ),
            layer(
                material(
                    300.0,
                    vec![(0.0, 0.0), (1.0, 0.4)],
                    if degenerate { [0.0, 0.12] } else { [0.02; 2] },
                ),
                0.1,
                vec![0.4 * initial[1]],
                cx,
            ),
        ],
        256,
        cx,
    )
    .unwrap()
}

#[test]
fn g1_finite_join_matches_independent_two_store_relaxation_in_both_directions() {
    with_cx(|cx| {
        for initial in [[0.8, 0.1], [0.1, 0.8]] {
            for h in [1e-12, 0.008, 1e60] {
                let old = stack(initial, false, cx)
                    .with_transfer([IdealContact; 2], vec![Finite(h / 0.001)], cx)
                    .unwrap();
                let dt = 0.7;
                let g = 1.0 / (1.0 / HALF[0] + 1.0 / h + 1.0 / HALF[1]);
                let difference = (initial[0] - initial[1])
                    / (1.0 + dt * g * (1.0 / CAPACITY[0] + 1.0 / CAPACITY[1]));
                let water = dt * g * difference;
                let expected = [
                    initial[0] - water / CAPACITY[0],
                    initial[1] + water / CAPACITY[1],
                ];
                let (next, report) = old
                    .advance(dt, &sealed(), config(), nonlinear(), cx)
                    .unwrap();
                for i in 0..2 {
                    near(next.activity()[i], expected[i], 3e-13);
                }
                let join = &report.interfaces[0];
                near(join.water_left_to_right_kg, water, 3e-14);
                near(
                    join.relative_activity,
                    expected[0] - water / (dt * HALF[0]),
                    3e-13,
                );
                near(
                    join.right_relative_activity,
                    expected[1] + water / (dt * HALF[1]),
                    3e-13,
                );
                near(next.water_mass_kg(), old.water_mass_kg(), 3e-14);
                if h == 0.008 {
                    near(
                        dt * h * (join.relative_activity - join.right_relative_activity),
                        water,
                        3e-14,
                    );
                    assert!((join.relative_activity - join.right_relative_activity).abs() > 0.1);
                }
                assert_eq!(report.nonlinear.iterations, 1);
            }
        }
    });
}

#[test]
fn g1_finite_surface_matches_implicit_reservoir_exchange_and_actual_surface() {
    with_cx(|cx| {
        for (initial, reservoir) in [(0.1, 0.9), (0.9, 0.1)] {
            let old = MoistureLayeredSlab::new(
                vec![layer(
                    material(800.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.01; 2]),
                    0.2,
                    vec![0.2 * initial],
                    cx,
                )],
                256,
                cx,
            )
            .unwrap()
            .with_transfer([Finite(4.0), Sealed], vec![], cx)
            .unwrap();
            let dt = 0.5;
            let h = 0.004;
            let g = 1.0 / (1.0 / HALF[0] + 1.0 / h);
            let expected = (CAPACITY[0] * initial + dt * g * reservoir) / (CAPACITY[0] + dt * g);
            let water = CAPACITY[0] * (expected - initial);
            let faces = [MoistureLayeredFace::PrescribedRelativeHumidity(reservoir); 2];
            let (next, report) = old.advance(dt, &faces, config(), nonlinear(), cx).unwrap();
            near(next.activity()[0], expected, 3e-13);
            near(report.boundary_water_kg[0], water, 3e-14);
            near(
                report.surface_activities[0].unwrap(),
                reservoir - water / (dt * h),
                3e-13,
            );
            assert_eq!(report.surface_activities[1], None);
            assert_eq!(report.boundary_water_kg[1], 0.0);
            near(next.water_mass_kg() - old.water_mass_kg(), water, 3e-14);
        }
    });
}

#[test]
fn g1_resistive_degenerate_paths_keep_nonlinear_kirchhoff_fluxes_and_wet_dry_cells() {
    with_cx(|cx| {
        for initial in [[0.8, 0.0], [0.0, 0.8]] {
            let old = stack(initial, true, cx)
                .with_transfer([Finite(5.0), Finite(2.0)], vec![Finite(8.0)], cx)
                .unwrap();
            let accepted_activity = old.activity().to_vec();
            let dt = 0.1;
            let faces = [
                MoistureLayeredFace::PrescribedRelativeHumidity(0.9),
                MoistureLayeredFace::PrescribedRelativeHumidity(0.1),
            ];
            let (next, report) = old.advance(dt, &faces, config(), nonlinear(), cx).unwrap();
            let a = next.activity();
            let face = &report.interfaces[0];
            let surfaces = [face.relative_activity, face.right_relative_activity];
            let half = [0.016, 0.096];
            let water = face.water_left_to_right_kg;
            near(
                dt * half[0] * (a[0] * a[0] - surfaces[0] * surfaces[0]),
                water,
                3e-13,
            );
            near(dt * 0.008 * (surfaces[0] - surfaces[1]), water, 3e-13);
            near(
                dt * half[1] * (surfaces[1] * surfaces[1] - a[1] * a[1]),
                water,
                3e-13,
            );
            for i in 0..2 {
                let surface = report.surface_activities[i].unwrap();
                let q = report.boundary_water_kg[i];
                near(dt * half[i] * (surface * surface - a[i] * a[i]), q, 3e-13);
                near(dt * [0.005, 0.002][i] * ([0.9, 0.1][i] - surface), q, 3e-13);
                assert!(a[i] > 0.0);
            }
            near(
                next.water_mass_kg() - old.water_mass_kg(),
                report.boundary_water_kg.iter().sum(),
                3e-13,
            );
            assert_eq!(old.time_s(), 0.0);
            assert_eq!(old.activity(), accepted_activity.as_slice());
            assert!(matches!(
                old.advance(
                    dt,
                    &faces,
                    config(),
                    MoistureNonlinearConfig {
                        maximum_iterations: 1,
                        ..nonlinear()
                    },
                    cx
                ),
                Err(MoistureError::NonlinearNotConverged(_))
            ));
        }
    });
}

#[test]
fn g1_finite_contact_crosses_a_zero_diffusivity_plateau_without_rejecting_dry_uptake() {
    with_cx(|cx| {
        let dry = layer(
            material_profile(
                800.0,
                vec![(0.0, 0.0), (1.0, 0.2)],
                vec![(0.0, 0.0), (0.1, 0.0), (0.6, 0.05)],
            ),
            0.2,
            vec![0.0],
            cx,
        );
        let wet = layer(
            material(300.0, vec![(0.0, 0.0), (1.0, 0.4)], [0.02; 2]),
            0.1,
            vec![0.32],
            cx,
        );
        let old = MoistureLayeredSlab::new(vec![dry, wet], 256, cx)
            .unwrap()
            .with_transfer([IdealContact; 2], vec![Finite(0.8)], cx)
            .unwrap();
        let dt = 0.1;
        // The initial full-path secant predicts a first left surface below
        // a=.5, entirely inside the true zero-D plateau. The physical solution
        // instead puts that surface above the threshold and wets the dry cell.
        let (next, report) = old
            .advance(dt, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        let join = &report.interfaces[0];
        let a = next.activity();
        assert!(a[0] > 0.0 && a[0] < 0.5);
        assert!(join.relative_activity > 0.5);
        assert!(join.relative_activity < join.right_relative_activity);
        let water = join.water_left_to_right_kg;
        assert!(water < 0.0);
        near(
            -dt * 0.016 * (join.relative_activity - 0.5).powi(2),
            water,
            3e-13,
        );
        near(
            dt * 0.0008 * (join.relative_activity - join.right_relative_activity),
            water,
            3e-13,
        );
        near(
            dt * HALF[1] * (join.right_relative_activity - a[1]),
            water,
            3e-13,
        );
        near(next.layers()[0].slab().water_mass_kg(), -water, 3e-13);
        near(next.water_mass_kg(), old.water_mass_kg(), 3e-13);
        assert_eq!(old.layers()[0].slab().moisture(), &[0.0]);
    });
}

#[test]
fn g3_finite_vapor_retains_series_resistance_and_actual_surface_source_receipts() {
    with_cx(|cx| {
        let initial = [0.2, 0.1];
        let moisture = stack(initial, false, cx)
            .with_transfer([Finite(4.0), Finite(6.0)], vec![Finite(8.0)], cx)
            .unwrap();
        let atmosphere = air::atmosphere(0.8, 0.1);
        let dt = 0.5;
        let vapor_capacity = atmosphere.vapor_mass_kg().unwrap() / 0.8;
        let outer = [
            1.0 / (1.0 / HALF[0] + 1.0 / 0.004),
            1.0 / (1.0 / HALF[1] + 1.0 / 0.006),
        ];
        let join = 1.0 / (1.0 / HALF[0] + 1.0 / 0.008 + 1.0 / HALF[1]);
        let d = [
            CAPACITY[0] + dt * (outer[0] + join),
            CAPACITY[1] + dt * (outer[1] + join),
        ];
        let b = dt * join;
        let determinant = d[0] * d[1] - b * b;
        let stored = [CAPACITY[0] * initial[0], CAPACITY[1] * initial[1]];
        let u = [
            (d[1] * stored[0] + b * stored[1]) / determinant,
            (b * stored[0] + d[0] * stored[1]) / determinant,
        ];
        let v = [
            dt * (d[1] * outer[0] + b * outer[1]) / determinant,
            dt * (b * outer[0] + d[0] * outer[1]) / determinant,
        ];
        let expected_rh = (vapor_capacity * 0.8
            + CAPACITY[0] * (initial[0] - u[0])
            + CAPACITY[1] * (initial[1] - u[1]))
            / (vapor_capacity + CAPACITY[0] * v[0] + CAPACITY[1] * v[1]);
        let candidate = atmosphere
            .advance_layered(
                &moisture,
                dt,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            expected_rh,
            3e-10,
        );
        for i in 0..2 {
            near(
                candidate.moisture().activity()[i],
                u[i] + v[i] * expected_rh,
                3e-10,
            );
            let surface = candidate.step().transport.surface_activities[i].unwrap();
            let source = candidate.step().sorption[i].as_ref().unwrap();
            near(
                source.property("sorption").unwrap().value_si(),
                [0.2, 0.4][i] * surface,
                3e-12,
            );
            assert!((surface - expected_rh).abs() > 0.001);
        }
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap() + candidate.moisture().water_mass_kg(),
            atmosphere.vapor_mass_kg().unwrap() + moisture.water_mass_kg(),
            3e-13,
        );
        near(
            air::carrier_mass(candidate.atmosphere()),
            air::carrier_mass(&atmosphere),
            3e-13,
        );
    });
}

#[test]
fn g0_ideal_and_sealed_limits_preserve_native_states_and_invalid_laws_refuse() {
    with_cx(|cx| {
        let old = stack([0.8, 0.1], false, cx);
        let ideal = old
            .clone()
            .with_transfer([IdealContact; 2], vec![IdealContact], cx)
            .unwrap();
        let (a, report_a) = old
            .advance(0.5, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        let (b, report_b) = ideal
            .advance(0.5, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(a.activity(), b.activity());
        assert_eq!(
            report_a.interfaces[0].water_left_to_right_kg,
            report_b.interfaces[0].water_left_to_right_kg
        );
        let sealed_stack = old
            .clone()
            .with_transfer([Sealed; 2], vec![Sealed], cx)
            .unwrap();
        let gas = air::atmosphere(0.5, 0.1);
        let candidate = gas
            .advance_layered(
                &sealed_stack,
                0.5,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.atmosphere(), &gas);
        assert_eq!(candidate.moisture().activity(), old.activity());
        assert_eq!(candidate.step().water_to_solid_kg, 0.0);
        assert!(candidate.step().sorption.iter().all(Option::is_none));
        for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                old.clone()
                    .with_transfer([Finite(h); 2], vec![IdealContact], cx)
                    .is_err()
            );
        }
        assert!(
            old.clone()
                .with_transfer([IdealContact; 2], vec![], cx)
                .is_err()
        );
        assert_eq!(old.time_s(), 0.0);
    });
}
