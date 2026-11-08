//! Three-dimensional native fields under adaptive passive source closure.
use super::*;
#[path = "../../common/adaptive_conductor.rs"]
mod adaptive_source;

fn adaptation() -> AdaptiveConductionRelaxation {
    AdaptiveConductionRelaxation {
        initial_weight: 0.5,
        minimum_weight: 0.02,
        maximum_weight: 1.0,
    }
}

#[test]
fn adaptive_resistance_closes_rotated_fields_and_equipotential_total_current() {
    // The common source has rho(J)=40J-35 near J=1. E=5 gives the
    // unique positive root J=1; the conductivity/E dual at J=5 gives E=1.
    // Fixed half relaxation is unstable at either root.
    for current_control in [false, true] {
        let (mut positions, tets) = slab(&[0.0, 1.0]);
        rotate(&mut positions);
        let original = SpatialConductionResponse::new(
            positions,
            tets,
            vec![adaptive_source::law(current_control); 6],
            limits(),
            &mut || false,
        )
        .unwrap();
        let fixed = original.clone().with_relaxation(0.5).unwrap();
        let adaptive = original
            .clone()
            .with_adaptive_relaxation(adaptation())
            .unwrap();
        let faces = face_sets(&original, &[0.0, 1.0]);
        let terminals = faces.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let matrix = SpatialConductanceControls {
            maximum_terminals: 2,
            maximum_field_solves: 1,
            maximum_response_values: 1000,
            maximum_dense_work: 1000,
            reference_voltage_v: 1.0,
            coefficient_tolerance_s: 1e-10,
        };
        let mixed = SpatialMixedConductanceControls {
            maximum_current_terminals: 1,
            maximum_response_values: 100,
            maximum_dense_work: 100,
        };
        for sign in [-1.0, 1.0] {
            let boundary = SpatialConductionBoundary {
                voltage: faces
                    .iter()
                    .enumerate()
                    .flat_map(|(terminal, faces)| {
                        faces.iter().map(move |&face| ConductionVoltageFace {
                            face,
                            voltage_v: [if terminal == 0 { 0.0 } else { sign * 5.0 }; 3],
                        })
                    })
                    .collect(),
                current: Vec::new(),
            };
            let voltages = [0.0];
            let currents = [sign * 5.0];
            let drive = if current_control {
                SpatialConductionResponseDrive::MixedTerminals {
                    terminals: &terminals,
                    current_terminals: &[1],
                    voltages_v: &voltages,
                    currents_a: &currents,
                    response: matrix,
                    mixed,
                }
            } else {
                SpatialConductionResponseDrive::Boundary(&boundary)
            };
            assert!(matches!(
                fixed.solve_candidate(
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    &mut || false,
                ),
                Err(SpatialConductionResponseError::NotConverged { .. }),
            ));
            let candidate = adaptive
                .solve_candidate(
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    &mut || false,
                )
                .unwrap();
            let report = candidate.report();
            assert!(report.iterations < 32);
            assert_eq!(report.field_solves, report.iterations);
            let (e, j) = if current_control {
                (1.0, 5.0)
            } else {
                (5.0, 1.0)
            };
            for cell in 0..6 {
                near(report.field.electric_field_v_m[cell][0], 0.0);
                near(report.field.electric_field_v_m[cell][1], 0.0);
                near(report.field.electric_field_v_m[cell][2], -sign * e);
                near(report.field.current_density_a_m2[cell][2], -sign * j);
                near(report.field.cell_power_w[cell], 5.0 / 6.0);
                near(
                    report.evaluated[cell]
                        .conductor
                        .material()
                        .state_coordinate("operating")
                        .unwrap(),
                    1.0,
                );
                assert_eq!(
                    report.evaluated[cell]
                        .conductor
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                    original.conductors()[cell]
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                );
            }
            near(report.field.power_w, 5.0);
            near(report.field.supplied_power_w, 5.0);
            near(report.field.current_residual_a, 0.0);
            assert_eq!(candidate.problem().adaptive_relaxation(), Some(adaptation()));
            let budget = SpatialConductionResponseControls {
                maximum_field_solves: 2,
                ..response_controls()
            };
            assert!(matches!(
                adaptive.solve_candidate(drive, limits(), controls(), budget, &mut || false),
                Err(SpatialConductionResponseError::FieldSolveBudget { completed: 2, .. }),
            ));
            assert!(matches!(
                adaptive.solve_candidate(
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    &mut || true,
                ),
                Err(SpatialConductionResponseError::Cancelled),
            ));
        }
        assert_eq!(adaptive.conductors(), original.conductors());
        assert_eq!(
            adaptive.with_relaxation(0.5).unwrap().adaptive_relaxation(),
            None,
        );
    }
}
