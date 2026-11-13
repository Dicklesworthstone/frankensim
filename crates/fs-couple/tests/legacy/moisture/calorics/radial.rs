//! Actual shell inventories and independent isothermal source integrals.
use super::*;
use fs_couple::moisture::radial::{MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry};

fn field(material: MaterialCard, initial: Vec<f64>, cx: &Cx<'_>) -> MoistureRadial {
    MoistureRadial::new(
        material,
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureRadialGeometry {
            symmetry: RadialMoistureSymmetry::Sphere,
            faces_m: vec![0.0, 0.005, 0.02],
            maximum_moisture: 0.3,
            maximum_cells: 2,
        },
        initial,
        cx,
    )
    .unwrap()
}

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn expected_heat(old: &MoistureRadial, new: &MoistureRadial, rt: f64, q: [f64; 2]) -> Vec<f64> {
    old.moisture()
        .iter()
        .zip(new.moisture())
        .zip(old.dry_masses_kg())
        .map(|((&a, &b), &dry)| dry * (b - a) * (rt - q[0] - (q[1] - q[0]) * (a + b) / 0.6))
        .collect()
}

#[test]
fn radial_sorption_energy_uses_each_unequal_shell_for_both_time_methods() {
    with_cx(|cx| {
        let q = [2e6, 3e6];
        let material = source(q);
        let energy = law(&material, [0.0, 0.3]);
        let original = field(material, vec![0.08, 0.13], cx);
        assert!(original.dry_masses_kg()[1] > 60.0 * original.dry_masses_kg()[0]);
        let faces = [None, Some(air::sorption(0.25))];
        for rh in [0.2, 0.8] {
            let atmosphere = air::atmosphere(rh, 1.0);
            for second_order in [false, true] {
                let (native_field, native_air, gas_change, calorics) = if second_order {
                    let expected = atmosphere
                        .advance_radial_sdirk2(
                            &original,
                            0.005,
                            &faces,
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap();
                    let result = atmosphere
                        .advance_radial_sdirk2_with_sorption_energy(
                            &original,
                            0.005,
                            &faces,
                            config(),
                            nonlinear(),
                            air::controls(),
                            &energy,
                            energy_controls(),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(
                        result.native().moisture().moisture(),
                        expected.moisture().moisture()
                    );
                    assert_eq!(result.native().atmosphere(), expected.atmosphere());
                    let (native, calorics) = result.into_parts();
                    let (air, field, step) = native.into_parts();
                    (field, air, step.gas_internal_energy_change_j, calorics)
                } else {
                    let expected = atmosphere
                        .advance_radial(
                            &original,
                            0.005,
                            &faces,
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap();
                    let result = atmosphere
                        .advance_radial_with_sorption_energy(
                            &original,
                            0.005,
                            &faces,
                            config(),
                            nonlinear(),
                            air::controls(),
                            &energy,
                            energy_controls(),
                            cx,
                        )
                        .unwrap();
                    assert_eq!(
                        result.native().moisture().moisture(),
                        expected.moisture().moisture()
                    );
                    assert_eq!(result.native().atmosphere(), expected.atmosphere());
                    let (native, calorics) = result.into_parts();
                    let (air, field, step) = native.into_parts();
                    (field, air, step.gas_internal_energy_change_j, calorics)
                };
                let reference = expected_heat(&original, &native_field, rt(&atmosphere), q);
                for (&actual, expected) in
                    calorics.cell_binding_energy_change_j.iter().zip(&reference)
                {
                    near(actual, *expected, 3e-9);
                }
                near(calorics.thermostat_heat_j, reference.iter().sum(), 3e-9);
                near(
                    calorics.mean_thermostat_power_w,
                    calorics.thermostat_heat_j / 0.005,
                    2e-8,
                );
                assert_eq!(calorics.gas_internal_energy_change_j, gas_change);
                near(
                    calorics.bound_internal_energy_change_j + gas_change
                        - calorics.thermostat_heat_j,
                    0.0,
                    energy_controls().energy_balance_tolerance_j,
                );
                near(
                    native_air.gas().mass_kg() + native_field.wet_mass_kg(),
                    atmosphere.gas().mass_kg() + original.wet_mass_kg(),
                    3e-12,
                );
                assert_eq!(calorics.interval_s, [0.0, 0.005]);
                assert_eq!(calorics.sorption_energy_source, energy.identity());
            }
        }
        assert_eq!(original.moisture(), &[0.08, 0.13]);
    });
}

#[test]
fn sealed_radial_redistribution_has_heat_duty_without_reservoir_water_exchange() {
    with_cx(|cx| {
        let q = [2e6, 3e6];
        let material = source(q);
        let energy = law(&material, [0.0, 0.3]);
        let original = field(material, vec![0.25, 0.05], cx);
        let atmosphere = air::atmosphere(0.5, 1.0);
        let result = atmosphere
            .advance_radial_with_sorption_energy(
                &original,
                0.01,
                &[None, None],
                config(),
                nonlinear(),
                air::controls(),
                &energy,
                energy_controls(),
                cx,
            )
            .unwrap();
        let next = result.native().moisture();
        assert!(next.moisture()[0] < original.moisture()[0]);
        assert!(next.moisture()[1] > original.moisture()[1]);
        near(next.water_mass_kg(), original.water_mass_kg(), 1e-13);
        assert_eq!(result.native().atmosphere(), &atmosphere);
        let expected: f64 = expected_heat(&original, next, rt(&atmosphere), q)
            .iter()
            .sum();
        assert!(expected.abs() > 1.0);
        near(result.energy().thermostat_heat_j, expected, 3e-9);
        near(result.energy().gas_internal_energy_change_j, 0.0, 1e-10);
    });
}

#[test]
fn late_radial_sorption_source_budget_failure_publishes_no_water_endpoint() {
    with_cx(|cx| {
        let material = source([2.5e6; 2]);
        let energy = law(&material, [0.0, 0.3]);
        let original = field(material, vec![0.1; 2], cx);
        let atmosphere = air::atmosphere(0.8, 1.0);
        let saved = atmosphere.clone();
        let mut budget = energy_controls();
        budget.maximum_source_segments = 0;
        let faces = [None, Some(air::sorption(0.25))];
        let failure = atmosphere
            .advance_radial_sdirk2_with_sorption_energy(
                &original,
                0.005,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                &energy,
                budget,
                cx,
            )
            .unwrap_err();
        assert!(matches!(failure, SorptionEnergyError::SourceBudget { .. }));
        assert_eq!(atmosphere, saved);
        assert_eq!(original.moisture(), &[0.1; 2]);
        assert_eq!(original.time_s(), 0.0);
    });
}
