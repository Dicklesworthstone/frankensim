//! Continuous-time heterogeneous references and complete-window accounting.
use super::super::atmosphere as air;
use super::*;
use fs_couple::moisture::adaptive::{MoistureAdvanceError, MoistureTimeControls};

const CAPACITY: [f64; 2] = [0.032, 0.012];
const HALF: [f64; 2] = [0.016, 0.048];

fn stack(initial: [f64; 2], slopes: [f64; 2], cx: &Cx<'_>) -> MoistureLayeredSlab {
    MoistureLayeredSlab::new(
        vec![
            layer(
                material(800.0, vec![(0.0, 0.0), (1.0, slopes[0])], [0.01; 2]),
                0.2,
                vec![slopes[0] * initial[0]],
                cx,
            ),
            layer(
                material(300.0, vec![(0.0, 0.0), (1.0, slopes[1])], [0.02; 2]),
                0.1,
                vec![slopes[1] * initial[1]],
                cx,
            ),
        ],
        256,
        cx,
    )
    .unwrap()
}
fn controls(maximum: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: maximum,
        minimum_step_s: 1e-7,
        maximum_step_s: maximum,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: tolerance * 1e-3,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

#[test]
fn g1_layered_time_refinement_resolves_contact_limited_exchange_and_layer_accounts() {
    with_cx(|cx| {
        let initial = [0.9, 0.1];
        let contact = 0.008;
        let old = stack(initial, [0.2, 0.4], cx)
            .with_transfer(
                [MoistureActivityTransfer::IdealContact; 2],
                vec![MoistureActivityTransfer::ConductanceKgM2S(contact / 0.001)],
                cx,
            )
            .unwrap();
        let saved = old.activity().to_vec();
        let conductance = 1.0 / (1.0 / HALF[0] + 1.0 / contact + 1.0 / HALF[1]);
        let equilibrium =
            (CAPACITY[0] * initial[0] + CAPACITY[1] * initial[1]) / (CAPACITY[0] + CAPACITY[1]);
        let decay = (-conductance * (1.0 / CAPACITY[0] + 1.0 / CAPACITY[1])).exp();
        let exact = initial.map(|a| equilibrium + (a - equilibrium) * decay);
        let error = |state: &MoistureLayeredSlab| {
            state
                .activity()
                .iter()
                .zip(exact)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max)
        };
        let (loose, loose_report) = old
            .advance_adaptive_to(
                1.0,
                &sealed(),
                config(),
                nonlinear(),
                controls(1.0, 2e-3),
                cx,
            )
            .unwrap();
        let (fine, report) = old
            .advance_adaptive_to(
                1.0,
                &sealed(),
                config(),
                nonlinear(),
                controls(1.0, 2e-5),
                cx,
            )
            .unwrap();
        assert!(error(&fine) < 0.3 * error(&loose));
        assert!(report.summary.evaluations > loose_report.summary.evaluations);
        assert_eq!(report.summary.interval_s, [0.0, 1.0]);
        assert_eq!(report.summary.boundary_water_kg, [0.0; 2]);
        assert_eq!(report.layers.len(), 2);
        assert_eq!(report.layers[0].id, 0);
        assert_eq!(report.layers[1].id, 1);
        assert!(report.layers[0].boundary_water_kg[1] < 0.0);
        assert_eq!(
            report.layers[0].boundary_water_kg[1],
            -report.layers[1].boundary_water_kg[0]
        );
        for layer in &report.layers {
            near(
                layer.water_change_kg,
                layer.boundary_water_kg.iter().sum(),
                1e-12,
            );
            assert!(layer.maximum_accepted_error_ratio <= 1.0);
        }
        near(fine.water_mass_kg(), old.water_mass_kg(), 1e-12);
        assert_eq!(old.activity(), saved);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn g3_local_activity_accuracy_cannot_hide_behind_shallow_sorption_storage() {
    with_cx(|cx| {
        let old = stack([0.9, 0.1], [1e-6, 2e-6], cx);
        let mut loose = controls(1.0, 1.0);
        let mut tight = loose;
        // Every moisture difference is below2e-6, so only the independently
        // controlled chemical driving coordinate can force this refinement.
        loose.activity_absolute_tolerance = 0.1;
        tight.activity_absolute_tolerance = 1e-5;
        let (_, a) = old
            .advance_adaptive_to(1.0, &sealed(), config(), nonlinear(), loose, cx)
            .unwrap();
        let (_, b) = old
            .advance_adaptive_to(1.0, &sealed(), config(), nonlinear(), tight, cx)
            .unwrap();
        assert!(b.summary.evaluations > a.summary.evaluations);
        assert!(b.summary.rejected_intervals > 0);
        assert!(
            b.layers
                .iter()
                .all(|layer| layer.maximum_accepted_error_ratio <= 1.0)
        );
    });
}

/// Exact exponential of the independent two-dimensional activity-difference
/// generator. Conservation then reconstructs the actual finite reservoir RH.
fn finite_reference(initial: [f64; 3], gas_capacity: f64, duration: f64) -> [f64; 3] {
    let join = 0.012;
    let a = -(join + HALF[0]) / CAPACITY[0] - HALF[0] / gas_capacity;
    let b = join / CAPACITY[0] - HALF[1] / gas_capacity;
    let c = join / CAPACITY[1] - HALF[0] / gas_capacity;
    let d = -(join + HALF[1]) / CAPACITY[1] - HALF[1] / gas_capacity;
    let split = ((a - d).powi(2) + 4.0 * b * c).sqrt();
    let roots = [(a + d + split) * 0.5, (a + d - split) * 0.5];
    let e = roots.map(|root| (root * duration).exp());
    let x = [initial[0] - initial[2], initial[1] - initial[2]];
    let mx = [a * x[0] + b * x[1], c * x[0] + d * x[1]];
    let next = std::array::from_fn::<_, 2, _>(|i| {
        ((mx[i] - roots[1] * x[i]) * e[0] - (mx[i] - roots[0] * x[i]) * e[1]) / split
    });
    let total = CAPACITY[0] * initial[0] + CAPACITY[1] * initial[1] + gas_capacity * initial[2];
    let humidity = (total - CAPACITY[0] * next[0] - CAPACITY[1] * next[1])
        / (CAPACITY[0] + CAPACITY[1] + gas_capacity);
    [humidity + next[0], humidity + next[1], humidity]
}

#[test]
fn g1_finite_vapor_layered_time_refinement_matches_three_store_continuous_solution() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.8, 0.1);
        let old = stack([0.2, 0.1], [0.2, 0.4], cx);
        let humidity = atmosphere.relative_humidity().unwrap();
        let capacity = atmosphere.vapor_mass_kg().unwrap() / humidity;
        let exact = finite_reference([0.2, 0.1, humidity], capacity, 0.5);
        let run = |time| {
            atmosphere
                .advance_layered_adaptive_to(
                    &old,
                    0.5,
                    [true; 2],
                    config(),
                    nonlinear(),
                    air::controls(),
                    time,
                    cx,
                )
                .unwrap()
        };
        let (loose_air, loose, _) = run(controls(0.5, 2e-3));
        let (fine_air, fine, report) = run(controls(0.5, 2e-5));
        let error = |air: &fs_couple::moisture::atmosphere::FiniteMoistureAtmosphere,
                     field: &MoistureLayeredSlab| {
            [
                field.activity()[0],
                field.activity()[1],
                air.relative_humidity().unwrap(),
            ]
            .into_iter()
            .zip(exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine_air, &fine) < 0.3 * error(&loose_air, &loose));
        assert_eq!(report.summary.interval_s, [0.0, 0.5]);
        assert_eq!(fine.time_s(), 0.5);
        assert_eq!(
            fine_air.gas().state().state().temperature,
            atmosphere.gas().state().state().temperature
        );
        near(
            air::carrier_mass(&fine_air),
            air::carrier_mass(&atmosphere),
            1e-12,
        );
        near(
            fine_air.vapor_mass_kg().unwrap() + fine.water_mass_kg(),
            atmosphere.vapor_mass_kg().unwrap() + old.water_mass_kg(),
            1e-10,
        );
        near(
            report.summary.water_change_kg,
            report.summary.boundary_water_kg.iter().sum(),
            1e-10,
        );
        assert!(report.summary.maximum_carrier_mass_residual_kg.unwrap() < 1e-10);
        assert!(report.summary.total_mass_residual_kg.unwrap().abs() < 1e-10);
        for (i, layer_report) in report.layers.iter().enumerate() {
            near(
                layer_report.water_change_kg,
                fine.layers()[i].slab().water_mass_kg() - old.layers()[i].slab().water_mass_kg(),
                1e-12,
            );
        }
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(atmosphere.relative_humidity().unwrap(), humidity);
    });
}

#[test]
fn g3_late_window_work_failure_preserves_all_original_layers_and_gas() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.8, 0.1);
        let saved_air = atmosphere.clone();
        let old = stack([0.2, 0.1], [0.2, 0.4], cx);
        let saved = old.activity().to_vec();
        let mut time = controls(0.1, 1.0);
        time.maximum_evaluations = 3;
        assert!(matches!(
            atmosphere.advance_layered_adaptive_to(
                &old,
                1.0,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                time,
                cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                accepted_intervals: 1,
                evaluations: 3,
                ..
            })
        ));
        assert_eq!(atmosphere, saved_air);
        assert_eq!(old.activity(), saved);
        assert_eq!(old.time_s(), 0.0);
        assert!(
            old.layers()
                .iter()
                .all(|layer| layer.slab().time_s() == 0.0)
        );
    });
}
