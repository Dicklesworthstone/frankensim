//! The original liquid-fraction coordinate also drives meridional conduction.
use super::*;
use fs_em::conduction::axisymmetric::AxisymmetricConduction;

#[test]
fn axisymmetric_current_and_voltage_refresh_the_original_latent_phase_source() {
    // Unit-radius, unit-height cylinder sector with angle=2 has volume=area=1.
    // Its lumped thermal receiver is the same original mass/enthalpy owner.
    for slope in [1.0, 3.0] {
        let card = card(slope, 1.0);
        let curve = curve(&card, 1.0);
        let thermal = initial(&curve, 10.0);
        let context = point(305.0, 0.0);
        let mesh = TriComplex2::from_indexed_triangles(
            "confined axisymmetric phase conductor",
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
            vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
            Metric2::axisymmetric(2.0).unwrap(),
        )
        .unwrap();
        let old = AxisymmetricConduction::new(
            mesh,
            vec![source(&card, &context, false); 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let cells = [binding(&card, &curve, &context)];
        let layout = network(&cells);
        // At fixed I=2 A: h1=10+20*(1+slope*f), h1=20+100*f.
        let fraction = 10.0 / (100.0 - 20.0 * slope);
        let enthalpy = 20.0 + 100.0 * fraction;
        let resistance = 1.0 + slope * fraction;
        let voltage = [
            (0, 0.0),
            (1, 0.0),
            (2, 2.0 * resistance),
            (3, 2.0 * resistance),
        ];
        for drive in [
            ConductionDrive::Current {
                positive: &[2, 3],
                negative: &[0, 1],
                current_a: 2.0,
            },
            ConductionDrive::Voltage(&voltage),
        ] {
            let next = layout
                .step_axisymmetric_candidate(
                    &thermal,
                    5.0,
                    &old,
                    drive,
                    limits(),
                    controls(),
                    config(),
                    &[20.0],
                    || false,
                )
                .unwrap();
            near(next.thermal().energy_j(), enthalpy);
            near(next.thermal().temperature_k(CellId(7)).unwrap(), 310.0);
            near(
                next.report().source.cells[0].phase.liquid_mass_fraction(),
                fraction,
            );
            near(next.report().electrical.heat_j, enthalpy - 10.0);
            near(next.report().electrical.supplied_work_j, enthalpy - 10.0);
            near(next.report().electrical.voltage_v[2], 2.0 * resistance);
            for (triangle, conductor) in next.problem().conductors().iter().enumerate() {
                assert_eq!(conductor.liquid_mass_fraction_axis(), Some(FRACTION));
                near(conductor.material().query().axes()["T"], 310.0);
                near(conductor.material().query().axes()[FRACTION], fraction);
                assert_eq!(conductor.frame(), old.conductors()[triangle].frame());
                for (a, b) in conductor
                    .source_properties()
                    .zip(old.conductors()[triangle].source_properties())
                {
                    assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
                }
                near(
                    next.report().electrical.current_density_a_m2[triangle][0],
                    0.0,
                );
                near(
                    next.report().electrical.current_density_a_m2[triangle][1],
                    -2.0,
                );
                near(
                    next.report().electrical.cell_power_w[triangle],
                    4.0 * resistance * old.mesh().face_measure(triangle).unwrap(),
                );
            }
            assert!(
                next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
            );
        }
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.energy_j(), 10.0);
    }
}
