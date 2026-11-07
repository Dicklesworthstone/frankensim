use super::*;
use fs_em::conduction::ConductionPowerSource;

#[test]
fn spatial_power_closes_nonlinear_sources_at_actual_amplitude_and_polarity() {
    for (basis, coordinate, intercept, slope, e) in [
        (
            ElectricalPropertyBasis::Conductivity,
            ConductorMagnitude::ElectricField,
            1.0,
            0.5,
            2.0,
        ),
        (
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            0.5,
            0.125,
            4.0,
        ),
    ] {
        let law = law_with_basis(Some(coordinate), basis, intercept, slope, 20.0);
        let (mut positions, tets) = slab(&[0.0, 1.0]);
        rotate(&mut positions);
        let problem =
            SpatialConductionResponse::new(positions, tets, vec![law; 6], limits(), &mut || false)
                .unwrap()
                .with_relaxation(0.5)
                .unwrap();
        let faces = face_sets(&problem, &[0.0, 1.0]);
        let source = ConductionPowerSource {
            power_w: 4.0 * e,
            maximum_voltage_v: 1.01 * e,
            maximum_current_a: 4.04,
        };
        for reference_voltage_v in [0.5, 3.0] {
            for reversed in [false, true] {
                let (positive, negative, sign) = if reversed { (0, 1, 1.0) } else { (1, 0, -1.0) };
                let candidate = problem
                    .solve_candidate(
                        SpatialConductionResponseDrive::Power {
                            positive_faces: &faces[positive],
                            negative_faces: &faces[negative],
                            source,
                            reference_voltage_v,
                        },
                        limits(),
                        controls(),
                        response_controls(),
                        &mut || false,
                    )
                    .unwrap();
                let report = candidate.report();
                near(report.field.power_w, source.power_w);
                near(report.field.supplied_power_w, source.power_w);
                for cell in 0..6 {
                    near(report.field.electric_field_v_m[cell][2], sign * e);
                    near(report.field.current_density_a_m2[cell][2], sign * 4.0);
                    near(report.field.cell_power_w[cell], source.power_w / 6.0);
                    near(
                        report.evaluated[cell].conductor.material().query().axes()["operating"],
                        if coordinate == ConductorMagnitude::ElectricField {
                            e
                        } else {
                            4.0
                        },
                    );
                }
                assert!(report.iterations > 1);
                assert_eq!(report.field_solves, report.iterations);
            }
        }
        let drive = SpatialConductionResponseDrive::Power {
            positive_faces: &faces[1],
            negative_faces: &faces[0],
            source,
            reference_voltage_v: 1.0,
        };
        assert!(matches!(
            problem.solve_candidate(
                drive,
                limits(),
                controls(),
                SpatialConductionResponseControls {
                    maximum_field_solves: 1,
                    ..response_controls()
                },
                &mut || false
            ),
            Err(SpatialConductionResponseError::FieldSolveBudget {
                completed: 1,
                required: 1,
                limit: 1
            })
        ));
        assert!(matches!(
            problem.solve_candidate(
                drive,
                limits(),
                controls(),
                response_controls(),
                &mut || true
            ),
            Err(SpatialConductionResponseError::Cancelled)
        ));
        for refused in [
            ConductionPowerSource {
                maximum_voltage_v: 0.99 * e,
                ..source
            },
            ConductionPowerSource {
                maximum_current_a: 3.99,
                ..source
            },
        ] {
            assert!(
                problem
                    .solve_candidate(
                        SpatialConductionResponseDrive::Power {
                            positive_faces: &faces[1],
                            negative_faces: &faces[0],
                            source: refused,
                            reference_voltage_v: 1.0,
                        },
                        limits(),
                        controls(),
                        response_controls(),
                        &mut || false
                    )
                    .is_err()
            );
        }
    }
}
