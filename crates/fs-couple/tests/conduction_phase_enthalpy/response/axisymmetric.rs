//! G1/G3: the original latent chart and source-owned meridional E/J close together.
use super::*;
use fs_em::conduction::axisymmetric::response::AxisymmetricConductionResponse;

fn radial_mesh() -> TriComplex2 {
    // Radius 1, length 1, angle 2: actual swept volume and both disk areas = 1.
    TriComplex2::from_indexed_triangles(
        "confined nonlinear phase cylinder",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::axisymmetric(2.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn axisymmetric_phase_field_response_continues_through_latent_and_sensible_heat() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let field_law = coordinate == ConductorMagnitude::ElectricField;
        for slope in [1.0, 3.0] {
            let card = source_card(slope, 1.0, coordinate);
            let caloric = curve(&card, 1.0);
            let mut thermal = initial(&caloric, 10.0);
            let context = field_point(305.0, 0.0);
            // Reuse the established source-law fixture, never its planar field.
            let sources = nonlinear(&card, &context).sources().to_vec();
            let mut problem =
                AxisymmetricConductionResponse::new(radial_mesh(), sources, limits(), &mut || {
                    false
                })
                .unwrap();
            let cells = [binding(&card, &caloric, &context)];
            let layout = network(&cells);
            let mut old_h = 10.0;
            for interval in 1..=2 {
                // I=2, A=L=1. E-dependent rho gives E=2*(1+s*f)/0.8;
                // J-dependent rho gives E=2*(1+s*f+0.2).
                let (base_heat, slope_heat) = if field_law {
                    (25.0, 25.0 * slope)
                } else {
                    (24.0, 20.0 * slope)
                };
                let trial_fraction = (old_h + base_heat - 20.0) / (100.0 - slope_heat);
                let (h, t, fraction) = if trial_fraction <= 1.0 {
                    (20.0 + 100.0 * trial_fraction, 310.0, trial_fraction)
                } else {
                    let h = old_h + base_heat + slope_heat;
                    (h, 310.0 + (h - 120.0) / 2.0, 1.0)
                };
                let rho = if field_law {
                    (1.0 + slope * fraction) / 0.8
                } else {
                    1.2 + slope * fraction
                };
                let next = layout
                    .step_axisymmetric_response_candidate(
                        &thermal,
                        5.0 * f64::from(interval),
                        &problem,
                        ConductionDrive::Current {
                            positive: &[2, 3],
                            negative: &[0, 1],
                            current_a: 2.0,
                        },
                        limits(),
                        controls(),
                        response_controls(),
                        config(),
                        &[h - old_h],
                        || false,
                    )
                    .unwrap();
                near(next.thermal().energy_j(), h);
                near(next.thermal().temperature_k(CellId(7)).unwrap(), t);
                near(
                    next.report().source.cells[0].phase.liquid_mass_fraction(),
                    fraction,
                );
                near(next.report().electrical.field.voltage_v[2], 2.0 * rho);
                near(next.report().electrical.heat_j, h - old_h);
                near(next.report().electrical.supplied_work_j, h - old_h);
                near(next.report().deposited_joule_j, h - old_h);
                for (triangle, conductor) in next.problem().conductors().iter().enumerate() {
                    near(conductor.material().query().axes()["T"], t);
                    near(conductor.material().query().axes()[FRACTION], fraction);
                    near(conductor.material().query().axes()["field"], 2.0 * rho);
                    near(conductor.material().query().axes()["current"], 2.0);
                    assert_eq!(conductor.liquid_mass_fraction_axis(), Some(FRACTION));
                    assert_eq!(
                        conductor.source_property().answer().receipt.selected,
                        problem.conductors()[triangle]
                            .source_property()
                            .answer()
                            .receipt
                            .selected
                    );
                    near(
                        next.report().electrical.field.current_density_a_m2[triangle][0],
                        0.0,
                    );
                    near(
                        next.report().electrical.field.current_density_a_m2[triangle][1],
                        -2.0,
                    );
                    near(
                        next.report().electrical.field.cell_power_w[triangle],
                        4.0 * rho * problem.mesh().face_measure(triangle).unwrap(),
                    );
                }
                assert!(
                    next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
                );
                assert!(
                    layout
                        .step_axisymmetric_response_candidate(
                            next.thermal(),
                            5.0 * f64::from(interval + 1),
                            &problem,
                            ConductionDrive::Current {
                                positive: &[2, 3],
                                negative: &[0, 1],
                                current_a: 2.0,
                            },
                            limits(),
                            controls(),
                            response_controls(),
                            config(),
                            &[1.0],
                            || false,
                        )
                        .is_err()
                );
                old_h = h;
                (thermal, problem, _) = next.into_parts();
            }
            near(context.axes()[FRACTION], 0.0);
            near(context.axes()["field"], 0.0);
        }
    }
}
