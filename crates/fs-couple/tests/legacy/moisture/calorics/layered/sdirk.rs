//! Independent heat primitives for native two-stage finite-air transport.
use super::*;

#[test]
fn sealed_sdirk_heat_uses_one_physical_endpoint_not_first_stage_diagnostic() {
    with_cx(|cx| {
        let capacity = [
            RHO * AREA * LENGTH[0] * SLOPE[0],
            RHO * AREA * LENGTH[1] * SLOPE[1],
        ];
        let conductance = 1.0
            / (0.5 * LENGTH[0] / (RHO * DIFFUSION * SLOPE[0] * AREA)
                + 1.0 / CONTACT
                + 0.5 * LENGTH[1] / (RHO * DIFFUSION * SLOPE[1] * AREA));
        let dt = 0.5;
        let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
        let rate = conductance * (1.0 / capacity[0] + 1.0 / capacity[1]);
        // The independent two-store activity difference decays with the
        // SDIRK stability function; conserved capacity fixes each endpoint.
        let z = -rate * dt;
        let decay = (1.0 + (1.0 - 2.0 * gamma) * z) / (1.0 - gamma * z).powi(2);
        for initial in [[0.8, 0.15], [0.15, 0.8]] {
            let (old, laws) = stack(CONSTANT_Q, initial.map(|a| vec![a]), cx);
            let air = air::atmosphere(0.5, 1.0);
            let result = air
                .advance_layered_sdirk2_with_sorption_energy(
                    &old,
                    dt,
                    [false; 2],
                    config(),
                    nonlinear(),
                    air::controls(),
                    &laws,
                    energy_controls(),
                    cx,
                )
                .unwrap();
            let transfer =
                (initial[0] - initial[1]) * (1.0 - decay) / (1.0 / capacity[0] + 1.0 / capacity[1]);
            let expected = [
                initial[0] - transfer / capacity[0],
                initial[1] + transfer / capacity[1],
            ];
            for (&actual, target) in result.native().moisture().activity().iter().zip(expected) {
                near(actual, target, 2e-12);
            }
            near(result.energy().layers[0].water_change_kg, -transfer, 3e-14);
            near(result.energy().layers[1].water_change_kg, transfer, 3e-14);
            let transport = &result.native().step().transport;
            near(
                transport.endpoint.interfaces[0].water_left_to_right_kg,
                transfer,
                3e-14,
            );
            let first = transport.first_stage.interfaces[0].water_left_to_right_kg;
            assert!(first.abs() > 1e-4);
            let expected_heat = (CONSTANT_Q[0][0] - CONSTANT_Q[1][0]) * transfer;
            check_whole(
                result.energy(),
                expected_heat,
                result.native().step().gas_internal_energy_change_j,
                dt,
            );
            assert_eq!(result.native().atmosphere(), &air);
            assert_eq!(result.energy().source_segments, 2);
            assert!(result.energy().thermostat_heat_j.abs() > 500.0);
            assert_eq!(
                result.energy().thermostat_heat_j.signum(),
                -transfer.signum()
            );
            assert_eq!(old.time_s(), 0.0);
        }
    });
}

#[test]
fn sdirk_finite_air_heat_integrates_loading_on_actual_planar_and_curved_cells() {
    with_cx(|cx| {
        for rh in [0.85, 0.05] {
            let air = air::atmosphere(rh, 1.0);
            let (old, laws) = stack(LOADING_Q, [vec![0.2, 0.5], vec![0.1, 0.3, 0.6]], cx);
            let dt = 0.05;
            let result = air
                .advance_layered_sdirk2_with_sorption_energy(
                    &old,
                    dt,
                    [true; 2],
                    config(),
                    nonlinear(),
                    air::controls(),
                    &laws,
                    energy_controls(),
                    cx,
                )
                .unwrap();
            let mut heat = 0.0;
            for i in 0..2 {
                let before = old.layers()[i].slab();
                let after = result.native().moisture().layers()[i].slab();
                let dry = vec![
                    before.dry_mass_kg() / before.moisture().len() as f64;
                    before.moisture().len()
                ];
                heat += check_layer(
                    &result.energy().layers[i],
                    before.moisture(),
                    after.moisture(),
                    &dry,
                    LOADING_Q[i],
                    rt(&air),
                );
                assert_eq!(
                    result.energy().layers[i].sorption_energy_source,
                    laws[i].identity()
                );
            }
            check_whole(
                result.energy(),
                heat,
                result.native().step().gas_internal_energy_change_j,
                dt,
            );
            assert_eq!(result.energy().source_segments, 5);
            near(
                result.native().atmosphere().vapor_mass_kg().unwrap()
                    + result.native().moisture().water_mass_kg(),
                air.vapor_mass_kg().unwrap() + old.water_mass_kg(),
                2e-12,
            );
            assert_eq!(
                result.native().step().water_to_solid_kg.signum(),
                (rh - 0.5_f64).signum()
            );
            for symmetry in [
                RadialMoistureSymmetry::Sphere,
                RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
            ] {
                let (old, laws) = coated(symmetry, cx);
                let dt = 1e-5;
                let result = air
                    .advance_layered_radial_sdirk2_with_sorption_energy(
                        &old,
                        dt,
                        [false, true],
                        config(),
                        nonlinear(),
                        air::controls(),
                        &laws,
                        energy_controls(),
                        cx,
                    )
                    .unwrap();
                let mut heat = 0.0;
                for i in 0..2 {
                    let before = old.layers()[i].radial();
                    let after = result.native().moisture().layers()[i].radial();
                    heat += check_layer(
                        &result.energy().layers[i],
                        before.moisture(),
                        after.moisture(),
                        before.dry_masses_kg(),
                        LOADING_Q[i],
                        rt(&air),
                    );
                    assert_eq!(
                        result.energy().layers[i].sorption_energy_source,
                        laws[i].identity()
                    );
                }
                check_whole(
                    result.energy(),
                    heat,
                    result.native().step().gas_internal_energy_change_j,
                    dt,
                );
                assert_eq!(result.energy().source_segments, 4);
                near(
                    result.native().atmosphere().vapor_mass_kg().unwrap()
                        + result.native().moisture().water_mass_kg(),
                    air.vapor_mass_kg().unwrap() + old.water_mass_kg(),
                    2e-12,
                );
                assert_eq!(old.time_s(), 0.0);
            }
        }
    });
}

#[test]
fn caloric_refusal_after_two_native_stages_discards_all_layered_endpoints() {
    with_cx(|cx| {
        let air = air::atmosphere(0.5, 1.0);
        let saved_air = air.clone();
        let (old, laws) = stack(CONSTANT_Q, [vec![0.8], vec![0.15]], cx);
        let saved_activity = old.activity().to_vec();
        // Both native stages fit; the common caloric budget is exhausted only
        // after integrating the first material's complete endpoint increment.
        let native = air
            .advance_layered_sdirk2(
                &old,
                0.5,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(native.moisture().time_s(), 0.5);
        assert!(matches!(
            air.advance_layered_sdirk2_with_sorption_energy(
                &old,
                0.5,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                &laws,
                SorptionEnergyControls {
                    maximum_source_segments: 1,
                    ..energy_controls()
                },
                cx,
            ),
            Err(SorptionEnergyError::SourceBudget {
                required: 2,
                maximum: 1
            })
        ));
        let narrow = [
            laws[0].clone(),
            law(&material(SLOPE[1], CONSTANT_Q[1]), [0.0, 0.045000000001]),
        ];
        assert!(matches!(
            air.advance_layered_sdirk2_with_sorption_energy(
                &old,
                0.5,
                [false; 2],
                config(),
                nonlinear(),
                air::controls(),
                &narrow,
                energy_controls(),
                cx,
            ),
            Err(SorptionEnergyError::Moisture(MoistureError::Material(_)))
        ));
        assert_eq!(old.activity(), saved_activity);
        assert_eq!(old.time_s(), 0.0);
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let (old, laws) = coated(symmetry, cx);
            let saved = old.activity().to_vec();
            let native = air
                .advance_layered_radial_sdirk2(
                    &old,
                    1e-5,
                    [false, true],
                    config(),
                    nonlinear(),
                    air::controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(native.moisture().time_s(), 1e-5);
            assert!(matches!(
                air.advance_layered_radial_sdirk2_with_sorption_energy(
                    &old,
                    1e-5,
                    [false, true],
                    config(),
                    nonlinear(),
                    air::controls(),
                    &laws,
                    SorptionEnergyControls {
                        maximum_source_segments: 2,
                        ..energy_controls()
                    },
                    cx,
                ),
                Err(SorptionEnergyError::SourceBudget {
                    required: 3,
                    maximum: 2
                })
            ));
            assert_eq!(old.activity(), saved);
            assert_eq!(old.time_s(), 0.0);
        }
        assert_eq!(air, saved_air);
    });
}
