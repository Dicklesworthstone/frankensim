//! G1/G3 finite-vapor radial stages against independent inventories and ODEs.

use super::*;

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

fn stage_rh(stage: &MoistureStepReport) -> f64 {
    stage.surface_sorption[1].as_ref().unwrap().query().axes()["relative_humidity"]
}

#[test]
fn finite_radial_sdirk_matches_closed_linear_decay_and_keeps_first_stage_source_bounds() {
    with_cx(|cx| {
        let radius: f64 = 0.05;
        let diffusion = 1e-4;
        let slope = 0.25;
        let initial = 0.04;
        let initial_rh = 0.8;
        let air = gas(initial_rh, 1.0);
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![initial],
            [diffusion; 2],
            cx,
        );
        let faces = [None, Some(sorption_law(slope))];
        let dry = 500.0 * 4.0 * core::f64::consts::PI * radius.powi(3) / 3.0;
        let capacity = air.vapor_mass_kg().unwrap() / initial_rh;
        near(field.dry_mass_kg(), dry, 1e-15);
        // Sphere center R/2 gives dM/dt=k*(s*RH-M). Eliminating RH
        // with the independent closed water inventory leaves one exponential.
        let k = 3.0 * diffusion / radius.powi(2);
        let water = dry * initial + capacity * initial_rh;
        let equilibrium = slope * water / (capacity + slope * dry);
        let rate = k * (1.0 + slope * dry / capacity);
        let dt = 0.25;
        let z = rate * dt;
        let first = equilibrium + (initial - equilibrium) / (1.0 + GAMMA * z);
        let first_rh = (water - dry * first) / capacity;
        let factor = (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2);
        let expected = equilibrium + (initial - equilibrium) * factor;
        let expected_rh = (water - dry * expected) / capacity;

        // The scalar gas solver visits RH=0 after its nonzero initial residual.
        // Its second-stage field still receives the accepted first-stage load.
        // This legitimate probe lies ABOVE old M despite its current dry face;
        // a source chart prepared only from old/current boundaries rejects it.
        let known = (1.0 - GAMMA) * dt * k * (slope * first_rh - first);
        let dry_probe = (initial + known) / (1.0 + GAMMA * dt * k);
        assert!(dry_probe > initial && dry_probe < slope * first_rh);
        assert_eq!(finite_controls().relative_humidity_bracket[0], 0.0);
        let candidate = air
            .advance_radial_sdirk2(
                &field,
                dt,
                &faces,
                config(),
                nonlinear(),
                finite_controls(),
                cx,
            )
            .unwrap();
        near(candidate.moisture().moisture()[0], expected, 3e-12);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            expected_rh,
            3e-11,
        );
        let report = &candidate.step().transport;
        near(report.first_stage.stored_water_kg[1] / dry, first, 3e-12);
        near(stage_rh(&report.first_stage), first_rh, 3e-11);
        near(stage_rh(&report.endpoint), expected_rh, 3e-11);
        assert!(first_rh < initial_rh && first_rh > expected_rh);
        near(
            report.first_stage.boundary_water_kg[1],
            GAMMA * dt * dry * k * (slope * first_rh - first),
            2e-12,
        );
        let exchange = dry * (expected - initial);
        near(report.endpoint.boundary_water_kg[1], exchange, 2e-12);
        near(candidate.step().water_to_solid_kg, exchange, 2e-12);
        near(
            air.vapor_mass_kg().unwrap() - candidate.atmosphere().vapor_mass_kg().unwrap(),
            exchange,
            2e-12,
        );
        assert_eq!(report.first_stage.interval_s, [0.0, GAMMA * dt]);
        assert_eq!(report.endpoint.interval_s, [0.0, dt]);

        let duration = 2.0;
        let exact = equilibrium + (initial - equilibrium) * (-rate * duration).exp();
        let errors = [8, 16, 32].map(|steps| {
            let mut current_air = air.clone();
            let mut current_field = field.clone();
            for _ in 0..steps {
                let (next_air, next_field, _) = current_air
                    .advance_radial_sdirk2(
                        &current_field,
                        duration / f64::from(steps),
                        &faces,
                        config(),
                        nonlinear(),
                        finite_controls(),
                        cx,
                    )
                    .unwrap()
                    .into_parts();
                near(
                    next_air.vapor_mass_kg().unwrap() + next_field.water_mass_kg(),
                    water,
                    5e-12,
                );
                current_air = next_air;
                current_field = next_field;
            }
            assert_eq!(current_field.time_s(), duration);
            near(
                current_air.relative_humidity().unwrap(),
                (water - dry * current_field.moisture()[0]) / capacity,
                4e-10,
            );
            (current_field.moisture()[0] - exact).abs()
        });
        assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
        assert!(errors[1] / errors[2] < 4.5, "{errors:?}");
        assert!(errors[2] < 1e-6, "{errors:?}");
        assert_eq!(field.moisture(), &[initial]);
        assert_eq!(field.time_s(), 0.0);
    });
}

#[test]
fn finite_radial_sdirk_nonlinear_curved_surface_uses_each_stages_actual_humidity() {
    with_cx(|cx| {
        let radius: f64 = 0.05;
        let radii = vec![0.0, 0.01, 0.03, radius];
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            radii,
            vec![0.04; 3],
            [0.0001, 0.0031],
            cx,
        );
        let air = gas(0.8, 0.4);
        let saved = air.clone();
        let dt = 0.1;
        let faces = [None, Some(crate::surface::law([0.001, 0.005]))];
        let candidate = air
            .advance_radial_sdirk2(
                &field,
                dt,
                &faces,
                config(),
                nonlinear(),
                finite_controls(),
                cx,
            )
            .unwrap();
        let report = &candidate.step().transport;
        let stages = [&report.first_stage, &report.endpoint];
        let area = 4.0 * core::f64::consts::PI * radius.powi(2);
        let center = 0.04;
        let radial_resistance = (1.0 / center - 1.0 / radius) / (4.0 * core::f64::consts::PI);
        // Independent K(M) for D(M)=0.0001+0.01*M, evaluated at the
        // physical surface, gives the same rate as its external film law.
        let potential = |m: f64| 0.0001 * m + 0.005 * m * m;
        let rates = stages.map(|stage| {
            let rh = stage_rh(stage);
            let m = stage.transport_end.last().unwrap().query().axes()["moisture_content"];
            let surface = stage.surface_moisture[1].unwrap();
            let h = 0.001 + 0.004 * rh;
            let source = stage.surface_transfer[1].as_ref().unwrap();
            near(source.query().axes()["relative_humidity"], rh, 1e-15);
            near(
                source
                    .property(fs_couple::moisture::surface::DRY_BASIS_SURFACE_EMISSION_PROPERTY)
                    .unwrap()
                    .value_si(),
                h,
                2e-15,
            );
            let bulk = 500.0 * (potential(surface) - potential(m)) / radial_resistance;
            let film = 500.0 * area * h * (0.25 * rh - surface);
            assert!(bulk > 0.0 && surface > m && surface < 0.25 * rh);
            near(bulk, film, 2e-12);
            film
        });
        let rh_first = stage_rh(stages[0]);
        let rh_final = stage_rh(stages[1]);
        assert!(rh_first - rh_final > 1e-4);
        near(
            report.first_stage.boundary_water_kg[1],
            GAMMA * dt * rates[0],
            2e-12,
        );
        let weighted = dt * ((1.0 - GAMMA) * rates[0] + GAMMA * rates[1]);
        near(report.endpoint.boundary_water_kg[1], weighted, 2e-12);
        near(
            local_water_change(&field, candidate.moisture()),
            weighted,
            2e-12,
        );
        near(candidate.step().water_to_solid_kg, weighted, 2e-12);
        near(
            air.vapor_mass_kg().unwrap() - candidate.atmosphere().vapor_mass_kg().unwrap(),
            weighted,
            2e-12,
        );
        let first_change: f64 = field
            .dry_masses_kg()
            .iter()
            .zip(&report.first_stage.transport_end)
            .zip(field.moisture())
            .map(|((&dry, source), &old)| dry * (source.query().axes()["moisture_content"] - old))
            .sum();
        let capacity = air.vapor_mass_kg().unwrap() / 0.8;
        near(capacity * (rh_first - 0.8) + first_change, 0.0, 2e-12);
        near(
            carrier_mass(candidate.atmosphere()),
            carrier_mass(&air),
            3e-13,
        );
        near(
            candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
            air.gas().mass_kg() + field.wet_mass_kg(),
            3e-12,
        );
        assert_eq!(air, saved);
        assert_eq!(field.moisture(), &[0.04; 3]);
        assert_eq!(field.time_s(), 0.0);
    });
}

#[test]
fn shared_radial_sdirk_orders_distinct_geometries_and_debits_vapor_once() {
    with_cx(|cx| {
        let dry = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.003, 0.009, 0.02],
            vec![0.01, 0.03, 0.06],
            [1e-5; 2],
            cx,
        );
        let wet = specimen(
            RadialMoistureSymmetry::Cylinder { length_m: 0.04 },
            vec![0.0, 0.005, 0.015, 0.025],
            vec![0.24, 0.22, 0.20],
            [2e-5; 2],
            cx,
        );
        let air = gas(0.5, 0.5);
        let dt = 0.01;
        let faces = [None, Some(sorption_law(0.25))];
        let forward = air
            .advance_shared_radial_sdirk2(
                &[participant(9, &dry, &faces), participant(2, &wet, &faces)],
                dt,
                6,
                finite_controls(),
                cx,
            )
            .unwrap();
        let reverse = air
            .advance_shared_radial_sdirk2(
                &[participant(2, &wet, &faces), participant(9, &dry, &faces)],
                dt,
                6,
                finite_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(forward.atmosphere(), reverse.atmosphere());
        assert_eq!(
            forward.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![2, 9]
        );
        let mut change = 0.0;
        let mut flux = 0.0;
        let first_rh = stage_rh(&forward.solids()[0].transport.first_stage);
        let end_rh = forward.atmosphere().relative_humidity().unwrap();
        for (a, b) in forward.solids().iter().zip(reverse.solids()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.moisture.moisture(), b.moisture.moisture());
            assert_eq!(a.moisture.time_s(), dt);
            let before = if a.id == 2 { &wet } else { &dry };
            let local = local_water_change(before, &a.moisture);
            assert_eq!(local.is_sign_positive(), a.id == 9);
            change += local;
            flux += a.transport.endpoint.boundary_water_kg.iter().sum::<f64>();
            near(stage_rh(&a.transport.first_stage), first_rh, 1e-15);
            near(stage_rh(&a.transport.endpoint), end_rh, 1e-14);
            assert_eq!(a.transport.first_stage.interval_s, [0.0, GAMMA * dt]);
            assert_eq!(a.transport.endpoint.interval_s, [0.0, dt]);
        }
        near(change, flux, 2e-12);
        near(forward.step().water_to_solids_kg, change, 2e-12);
        near(
            air.vapor_mass_kg().unwrap() - forward.atmosphere().vapor_mass_kg().unwrap(),
            change,
            2e-12,
        );
        near(forward.step().water_balance_residual_kg, 0.0, 1e-13);
        near(forward.step().total_mass_residual_kg, 0.0, 1e-13);
        near(
            carrier_mass(forward.atmosphere()),
            carrier_mass(&air),
            3e-13,
        );
        assert_eq!(dry.time_s(), 0.0);
        assert_eq!(wet.time_s(), 0.0);
    });
}

#[test]
fn finite_radial_sdirk_refuses_nonphysical_second_stage_atomically() {
    with_cx(|cx| {
        let radius: f64 = 0.05;
        let diffusion = 0.001;
        let initial = 0.2;
        let air = gas(0.01, 100.0);
        let saved = air.clone();
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![initial],
            [diffusion; 2],
            cx,
        );
        let dt = 10.0 * radius.powi(2) / (3.0 * diffusion);
        let faces = [None, Some(sorption_law(0.25))];
        let mut controls = finite_controls();
        controls.mass_tolerance_kg = 1e-11;
        let first = air
            .advance_radial(
                &field,
                GAMMA * dt,
                &faces,
                config(),
                nonlinear(),
                controls,
                cx,
            )
            .unwrap();
        assert!(first.moisture().moisture()[0] > 0.0);
        assert!(first.moisture().moisture()[0] < initial);
        assert!(first.atmosphere().relative_humidity().unwrap() > 0.01);
        // The exact dissipative ODE stays positive; a stiff SDIRK stage has a
        // negative amplification here. Refuse it instead of clipping, mutating
        // the accepted gas, or silently substituting the positive BE endpoint.
        let capacity = air.vapor_mass_kg().unwrap() / 0.01;
        let dry = field.dry_mass_kg();
        let equilibrium = 0.25 * (dry * initial + capacity * 0.01) / (capacity + 0.25 * dry);
        let z = 10.0 * (1.0 + 0.25 * dry / capacity);
        let factor = (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2);
        assert!(equilibrium + (initial - equilibrium) * factor < 0.0);
        assert!(matches!(
            air.advance_radial_sdirk2(&field, dt, &faces, config(), nonlinear(), controls, cx),
            Err(AtmosphereMoistureError::Moisture(
                MoistureError::StageRange { .. }
            ))
        ));
        assert_eq!(air, saved);
        assert_eq!(field.moisture(), &[initial]);
        assert_eq!(field.time_s(), 0.0);
    });
}
