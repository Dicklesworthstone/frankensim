use super::*;

#[test]
fn held_power_closes_actual_e_and_j_before_enforcing_supply_compliance() {
    // Unit length/width, depth 0.2 m. Two independently known operating
    // points: sigma(E)=1+0.5E at E=2,J=4; rho(J)=0.5+0.125J at E=4,J=4.
    // Initial coefficient iterates exceed V compliance in the first case and
    // I compliance in the second; both physical fixed points are compliant.
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
        let law = response(basis, coordinate, intercept, slope, 20.0);
        let problem = specimen(law.clone(), law, 0.2)
            .with_relaxation(0.5)
            .unwrap();
        let source = ConductionPowerSource {
            power_w: 0.8 * e,
            maximum_voltage_v: 1.01 * e,
            maximum_current_a: 0.81,
        };
        for reversed in [false, true] {
            let (positive, negative, sign) = if reversed {
                (&[0, 3][..], &[2, 5][..], 1.0)
            } else {
                (&[2, 5][..], &[0, 3][..], -1.0)
            };
            let candidate = problem
                .solve_power_candidate(
                    positive,
                    negative,
                    source,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let field = candidate.report();
            near(field.power_w, source.power_w);
            near(field.power_balance_residual_w, 0.0);
            for cell in 0..4 {
                near(field.electric_field_v_m[cell][0], sign * e);
                near(field.current_density_a_m2[cell][0], sign * 4.0);
                near(field.cell_power_w[cell], source.power_w / 4.0);
                near(
                    field.evaluated[cell].conductor.material().query().axes()["operating"],
                    if coordinate == ConductorMagnitude::ElectricField {
                        e
                    } else {
                        4.0
                    },
                );
            }
            assert!(field.iterations > 1);
            assert_eq!(field.field_solves, field.iterations);
            assert_eq!(candidate.problem().relaxation(), 0.5);
        }
        // The original, unadvanced model remains available after refusals.
        for refused in [
            ConductionPowerSource {
                maximum_voltage_v: 0.99 * e,
                ..source
            },
            ConductionPowerSource {
                maximum_current_a: 0.79,
                ..source
            },
            ConductionPowerSource {
                power_w: f64::NAN,
                ..source
            },
        ] {
            assert!(
                problem
                    .solve_power_candidate(
                        &[2, 5],
                        &[0, 3],
                        refused,
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .is_err()
            );
        }
        assert!(matches!(
            problem.solve_power_candidate(
                &[2, 5],
                &[0, 3],
                source,
                limits(),
                controls(),
                policy(),
                &mut || true,
            ),
            Err(ConductionResponseError::Cancelled)
        ));
        let zero = problem
            .solve_power_candidate(
                &[2, 5],
                &[0, 3],
                ConductionPowerSource {
                    power_w: 0.0,
                    ..source
                },
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .unwrap();
        near(zero.report().power_w, 0.0);
        assert!(
            zero.report()
                .current_density_a_m2
                .iter()
                .all(|j| *j == [0.0; 2])
        );
    }
}
