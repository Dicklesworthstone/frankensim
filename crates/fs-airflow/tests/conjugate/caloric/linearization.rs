//! G1/G3/G4: source-Cp path actions against independent heat-exchange physics.

use super::*;
use fs_airflow::conjugate::caloric::{CaloricPathTangent, CaloricSegmentTangent};

// Cp=800+2T, mass flow=.02 kg/s, UA=40/N W/K. Solve the independent
// closed-form exchange primitive in relative travel z, then integrate Cp
// directly. This reference invokes no production material or marching API.
fn reference(inlet: f64, walls: &[f64]) -> (Vec<[f64; 3]>, [f64; 2]) {
    let conductance = 40.0 / walls.len() as f64;
    let exposure = conductance / 0.02;
    let mut previous = inlet;
    let mut heat = 0.0;
    let mut segments = Vec::with_capacity(walls.len());
    for &wall in walls {
        let (mut low, mut high) = (0.0, 1.0);
        for _ in 0..80 {
            let z: f64 = 0.5 * (low + high);
            let integral = -(800.0 + 2.0 * wall) * (-z).ln_1p() - 2.0 * (wall - previous) * z;
            if integral < exposure {
                low = z;
            } else {
                high = z;
            }
        }
        let outlet = previous + (wall - previous) * (0.5 * (low + high));
        let local_heat = 0.02 * (outlet - previous) * (800.0 + previous + outlet);
        segments.push([outlet, local_heat, wall - local_heat / conductance]);
        heat += local_heat;
        previous = outlet;
    }
    (segments, [previous, heat])
}

fn objective(
    inlet: f64,
    walls: &[f64],
    seeds: &[CaloricSegmentTangent],
    terminal: CaloricPathTangent,
) -> f64 {
    let (segments, final_state) = reference(inlet, walls);
    segments
        .iter()
        .zip(seeds)
        .map(|(state, seed)| {
            seed.outlet_temperature_k * state[0]
                + seed.heat_rate_w * state[1]
                + seed.reference_temperature_k * state[2]
        })
        .sum::<f64>()
        + terminal.outlet_temperature_k * final_state[0]
        + terminal.total_heat_rate_w * final_state[1]
}

#[test]
fn g1_variable_cp_heating_and_cooling_actions_match_independent_forward_and_objective_differences()
{
    let law = source(&[(250.0, 1300.0), (400.0, 1600.0), (650.0, 2100.0)], 300.0);
    let original = path(law, 300.0, 4);
    let walls = [620.0, 270.0, 590.0, 310.0];
    let direction = [0.2, -0.7, 0.5, 1.1];
    let inlet_direction = -0.4;
    let epsilon = 1e-3;
    with_cx(|cx| {
        let (march, linear) = original.linearize(cx, &walls).unwrap();
        let (expected, final_expected) = reference(300.0, &walls);
        for (index, state) in march.segments.iter().enumerate() {
            near(state.outlet_temperature_k, expected[index][0], 2e-9);
            near(state.heat_rate_w, expected[index][1], 2e-7);
            near(state.reference_temperature_k, expected[index][2], 2e-9);
            assert_eq!(state.heat_rate_w > 0.0, index % 2 == 0);
        }
        near(march.outlet_temperature_k, final_expected[0], 2e-9);
        near(march.total_heat_rate_w, final_expected[1], 2e-7);

        let plus_walls = core::array::from_fn::<_, 4, _>(|i| walls[i] + epsilon * direction[i]);
        let minus_walls = core::array::from_fn::<_, 4, _>(|i| walls[i] - epsilon * direction[i]);
        let (plus, final_plus) = reference(300.0 + epsilon * inlet_direction, &plus_walls);
        let (minus, final_minus) = reference(300.0 - epsilon * inlet_direction, &minus_walls);
        let mut output = [CaloricSegmentTangent::default(); 4];
        let total = linear
            .apply(cx, inlet_direction, &direction, &mut output)
            .unwrap();
        for (index, state) in output.iter().enumerate() {
            for (component, actual, tolerance) in [
                (0, state.outlet_temperature_k, 5e-8),
                (1, state.heat_rate_w, 2e-6),
                (2, state.reference_temperature_k, 5e-8),
            ] {
                near(
                    actual,
                    (plus[index][component] - minus[index][component]) / (2.0 * epsilon),
                    tolerance,
                );
            }
        }
        near(
            total.outlet_temperature_k,
            (final_plus[0] - final_minus[0]) / (2.0 * epsilon),
            5e-8,
        );
        near(
            total.total_heat_rate_w,
            (final_plus[1] - final_minus[1]) / (2.0 * epsilon),
            2e-6,
        );

        let seeds = core::array::from_fn::<_, 4, _>(|i| CaloricSegmentTangent {
            outlet_temperature_k: 0.2 * (i + 1) as f64,
            heat_rate_w: [0.01, -0.02, 0.03, -0.04][i],
            reference_temperature_k: [-0.7, 0.5, -0.3, 0.1][i],
        });
        let terminal = CaloricPathTangent {
            outlet_temperature_k: 0.6,
            total_heat_rate_w: 0.025,
        };
        let mut wall_gradient = [0.0; 4];
        let inlet_gradient = linear
            .apply_transpose(cx, &seeds, terminal, &mut wall_gradient)
            .unwrap();
        // Every returned gradient is checked against a separate scalar physical
        // objective, rather than only pairing forward and transpose operators.
        for input in 0..5 {
            let (mut plus_walls, mut minus_walls) = (walls, walls);
            let (mut plus_inlet, mut minus_inlet) = (300.0, 300.0);
            if input == 0 {
                plus_inlet += epsilon;
                minus_inlet -= epsilon;
            } else {
                plus_walls[input - 1] += epsilon;
                minus_walls[input - 1] -= epsilon;
            }
            let expected = (objective(plus_inlet, &plus_walls, &seeds, terminal)
                - objective(minus_inlet, &minus_walls, &seeds, terminal))
                / (2.0 * epsilon);
            near(
                if input == 0 {
                    inlet_gradient
                } else {
                    wall_gradient[input - 1]
                },
                expected,
                2e-7,
            );
        }
    });
}

fn relative(actual: f64, expected: f64) {
    // A vanishing weak-transfer or saturated-inlet derivative must not pass
    // simply because the physically correct derivative is much smaller than 1.
    near(
        actual,
        expected,
        if expected == 0.0 {
            2e-14
        } else {
            3e-11 * expected.abs()
        },
    );
}

#[test]
fn g1_constant_cp_actions_retain_weak_transfer_saturated_inlet_and_equilibrium_limits() {
    let law = source(&[(250.0, 1005.0), (650.0, 1005.0)], 300.0);
    let capacity_rate = 0.02 * 1005.0;
    for ntu in [1e-12_f64, 1e-4, 2.0, 80.0] {
        let gain = (-ntu).exp();
        let loss = -(-ntu).exp_m1();
        let mean_wall = if ntu < 1e-3 {
            ntu * (0.5 - ntu * (1.0 / 6.0 - ntu * (1.0 / 24.0 - ntu * (1.0 / 120.0 - ntu / 720.0))))
        } else {
            1.0 - loss / ntu
        };
        let mean_inlet = loss / ntu;
        for (inlet, wall) in [(300.0, 480.0), (560.0, 280.0), (400.0, 400.0)] {
            let original = CaloricAirPath::new(
                inlet,
                0.02,
                law.clone(),
                vec![AirSegment::new("analytic wall", 1.0, capacity_rate * ntu).unwrap()],
            )
            .unwrap();
            let walls = [wall];
            with_cx(|cx| {
                let (march, linear) = original.linearize(cx, &walls).unwrap();
                if inlet == wall {
                    assert_eq!(march.total_heat_rate_w, 0.0);
                }
                for (inlet_direction, wall_direction) in [(1.0, 0.0), (0.0, 1.0)] {
                    let mut output = [CaloricSegmentTangent::default()];
                    let total = linear
                        .apply(cx, inlet_direction, &[wall_direction], &mut output)
                        .unwrap();
                    let expected_temperature = gain * inlet_direction + loss * wall_direction;
                    let expected_heat = capacity_rate * loss * (wall_direction - inlet_direction);
                    let expected_mean = mean_inlet * inlet_direction + mean_wall * wall_direction;
                    relative(output[0].outlet_temperature_k, expected_temperature);
                    relative(output[0].heat_rate_w, expected_heat);
                    relative(output[0].reference_temperature_k, expected_mean);
                    relative(total.outlet_temperature_k, expected_temperature);
                    relative(total.total_heat_rate_w, expected_heat);
                }
                let seeds = [CaloricSegmentTangent {
                    outlet_temperature_k: 0.7,
                    heat_rate_w: 0.03,
                    reference_temperature_k: 0.4,
                }];
                let terminal = CaloricPathTangent {
                    outlet_temperature_k: -0.2,
                    total_heat_rate_w: 0.02,
                };
                let mut wall_gradient = [0.0];
                let inlet_gradient = linear
                    .apply_transpose(cx, &seeds, terminal, &mut wall_gradient)
                    .unwrap();
                relative(
                    inlet_gradient,
                    0.5 * gain - 0.05 * capacity_rate * loss + 0.4 * mean_inlet,
                );
                relative(
                    wall_gradient[0],
                    0.5 * loss + 0.05 * capacity_rate * loss + 0.4 * mean_wall,
                );
            });
        }
    }
}

#[test]
fn g4_source_endpoints_shapes_nonfinite_directions_and_cancellation_refuse_without_poisoning_reuse()
{
    let narrow = source(&[(300.0, 1000.0), (400.0, 1000.0)], 300.0);
    let edge = path(narrow.clone(), 300.0, 1);
    with_cx(|cx| {
        assert!(edge.march(&[350.0]).is_ok());
        assert!(matches!(
            edge.linearize(cx, &[350.0]),
            Err(AirflowError::InvalidConjugateInput { .. })
        ));
        let saturated = CaloricAirPath::new(
            350.0,
            1.0,
            narrow.clone(),
            vec![AirSegment::new("upper endpoint", 1.0, 80000.0).unwrap()],
        )
        .unwrap();
        assert_eq!(
            saturated.march(&[400.0]).unwrap().outlet_temperature_k,
            400.0
        );
        assert!(matches!(
            saturated.linearize(cx, &[400.0]),
            Err(AirflowError::InvalidConjugateInput { .. })
        ));
        // A wall is a prescribed reservoir, not an extra coolant state.
        let external = CaloricAirPath::new(
            350.0,
            1.0,
            narrow,
            vec![AirSegment::new("external wall", 1.0, 0.01).unwrap()],
        )
        .unwrap();
        assert!(external.linearize(cx, &[1000.0]).is_ok());
    });
    let original = path(source(&[(250.0, 1300.0), (650.0, 2100.0)], 300.0), 300.0, 2);
    let walls = [480.0, 320.0];
    let (_, linear) = with_cx(|cx| original.linearize(cx, &walls)).unwrap();
    let mut output = [CaloricSegmentTangent::default(); 2];
    let mut gradient = [0.0; 2];
    let seeds = [CaloricSegmentTangent::default(); 2];
    let terminal = CaloricPathTangent::default();
    with_cx(|cx| {
        assert!(original.linearize(cx, &[480.0]).is_err());
        assert!(original.linearize(cx, &[480.0, f64::NAN]).is_err());
        assert!(linear.apply(cx, 0.0, &[1.0], &mut output).is_err());
        assert!(linear.apply(cx, 0.0, &[1.0; 2], &mut output[..1]).is_err());
        assert!(linear.apply(cx, f64::NAN, &[1.0; 2], &mut output).is_err());
        assert!(
            linear
                .apply(cx, 0.0, &[1.0, f64::INFINITY], &mut output)
                .is_err()
        );
        assert!(
            linear
                .apply_transpose(cx, &seeds[..1], terminal, &mut gradient)
                .is_err()
        );
        assert!(
            linear
                .apply_transpose(cx, &seeds, terminal, &mut gradient[..1])
                .is_err()
        );
        for bad in [
            CaloricSegmentTangent {
                outlet_temperature_k: f64::NAN,
                ..CaloricSegmentTangent::default()
            },
            CaloricSegmentTangent {
                heat_rate_w: f64::INFINITY,
                ..CaloricSegmentTangent::default()
            },
            CaloricSegmentTangent {
                reference_temperature_k: f64::NAN,
                ..CaloricSegmentTangent::default()
            },
        ] {
            assert!(
                linear
                    .apply_transpose(cx, &[seeds[0], bad], terminal, &mut gradient)
                    .is_err()
            );
        }
        assert!(
            linear
                .apply_transpose(
                    cx,
                    &seeds,
                    CaloricPathTangent {
                        total_heat_rate_w: f64::NAN,
                        ..terminal
                    },
                    &mut gradient
                )
                .is_err()
        );
    });
    let expected = with_cx(|cx| linear.apply(cx, 0.3, &[0.2, -0.1], &mut output)).unwrap();
    let expected_segments = output;
    let gate = CancelGate::new();
    gate.request();
    with_gated_cx(&gate, |cx| {
        assert!(matches!(
            original.linearize(cx, &walls),
            Err(AirflowError::Cancelled { .. })
        ));
        assert!(matches!(
            linear.apply(cx, 0.3, &[0.2, -0.1], &mut output),
            Err(AirflowError::Cancelled { .. })
        ));
        assert!(matches!(
            linear.apply_transpose(cx, &seeds, terminal, &mut gradient),
            Err(AirflowError::Cancelled { .. })
        ));
    });
    let repeated = with_cx(|cx| linear.apply(cx, 0.3, &[0.2, -0.1], &mut output)).unwrap();
    assert_eq!(repeated, expected);
    assert_eq!(output, expected_segments);
}
