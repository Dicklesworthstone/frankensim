//! Analytic finite-vapor SDIRK limits and conservative physical stage accounts.
//! Gas, sorption and solid coefficients are authored verification data.
use super::*;
use fs_couple::moisture::atmosphere::AtmosphereMoistureError;

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

fn decay(z: f64) -> f64 {
    (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2)
}

fn vapor_capacity(volume_m3: f64) -> f64 {
    let temperature = 293.15;
    fs_material::gas::saturation_pressure_water_pa(temperature).unwrap() * 0.018_015_28 * volume_m3
        / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature)
}

#[test]
fn g1_kirchhoff_finite_vapor_retains_the_first_stage_source_envelope() {
    with_cx(|cx| {
        let diffusion = 1e-4;
        let geometry = grid_geometry([0.1; 3], [1; 3]);
        let initial_moisture = 0.04;
        let slope = 0.25;
        let initial_rh = 0.8;
        let initial = grid(
            card([diffusion; 2], 0.3),
            geometry,
            vec![initial_moisture],
            cx,
        )
        .with_kirchhoff_diffusivity();
        let atmosphere = air::atmosphere(initial_rh, 1.0);
        let faces = std::array::from_fn(|_| Some(air::sorption(slope)));
        let dry = initial.dry_mass_kg();
        let capacity = vapor_capacity(1.0);
        let rate = 12.0 * diffusion / geometry.lengths_m[0].powi(2);
        let total_water = dry * initial_moisture + capacity * initial_rh;
        let equilibrium = slope * total_water / (capacity + slope * dry);
        let lambda = rate * (1.0 + slope * dry / capacity);
        let dt = 0.25;
        let first = equilibrium
            + (initial_moisture - equilibrium) / (1.0 + GAMMA * lambda * dt);
        let first_rh = (total_water - dry * first) / capacity;
        let known = (1.0 - GAMMA) * dt * rate * (slope * first_rh - first);
        // The reservoir root visits RH=0 in stage two. The physical first-stage
        // influx persists at that dry probe, so its field exceeds the old/current
        // boundary envelope while remaining inside the admitted stage envelope.
        let dry_probe = (initial_moisture + known) / (1.0 + GAMMA * dt * rate);
        assert!(dry_probe > initial_moisture && dry_probe < slope * first_rh);
        let expected = equilibrium + (initial_moisture - equilibrium) * decay(lambda * dt);
        let candidate = atmosphere
            .advance_cartesian_sdirk2(
                &initial,
                dt,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        let report = &candidate.step().transport;
        near(candidate.moisture().moisture()[0], expected, 3e-12);
        near(report.first_stage.stored_water_kg[1] / dry, first, 3e-12);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            (total_water - dry * expected) / capacity,
            4e-11,
        );
        assert!(report.first_stage.kirchhoff_sources.is_some());
        assert!(report.endpoint.kirchhoff_sources.is_some());
        let exchange = dry * (expected - initial_moisture);
        near(report.endpoint.boundary_water_kg.iter().sum(), exchange, 2e-12);
        near(candidate.step().water_to_solid_kg, exchange, 2e-12);
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap() + candidate.moisture().water_mass_kg(),
            total_water,
            2e-12,
        );
        assert_eq!(initial.moisture(), &[initial_moisture]);
        assert_eq!(initial.time_s(), 0.0);
        near(atmosphere.relative_humidity().unwrap(), initial_rh, 2e-15);
    });
}

#[test]
fn g1_finite_cartesian_sdirk_has_second_order_adsorption_and_desorption() {
    with_cx(|cx| {
        let diffusion = 0.01;
        let geometry = grid_geometry([0.1; 3], [1; 3]);
        let initial = grid(card([diffusion; 2], 0.3), geometry, vec![0.1], cx);
        let slope = 0.25;
        let faces = std::array::from_fn(|_| Some(air::sorption(slope)));
        let dry = initial.dry_mass_kg();
        let capacity = vapor_capacity(1.0);
        let rate = 12.0 * diffusion / geometry.lengths_m[0].powi(2);
        let lambda = rate * (1.0 + slope * dry / capacity);
        let duration = 0.01;

        for rh in [0.8, 0.2] {
            let original_air = air::atmosphere(rh, 1.0);
            let r0 = original_air.relative_humidity().unwrap();
            near(original_air.vapor_mass_kg().unwrap(), capacity * r0, 2e-15);
            // Eliminate vapor from the two physical inventories:
            // M' = -lambda*(M-M_eq), RH = RH0-dry*(M-M0)/capacity.
            let equilibrium = slope * (dry * 0.1 + capacity * r0) / (capacity + slope * dry);
            let exact = equilibrium + (0.1 - equilibrium) * (-lambda * duration).exp();
            let errors = [4usize, 8, 16].map(|steps| {
                let dt = duration / steps as f64;
                let mut atmosphere = original_air.clone();
                let mut field = initial.clone();
                for _ in 0..steps {
                    let (next_air, next_field, report) = atmosphere
                        .advance_cartesian_sdirk2(
                            &field,
                            dt,
                            &faces,
                            config(),
                            nonlinear(),
                            air::controls(),
                            cx,
                        )
                        .unwrap()
                        .into_parts();
                    assert_eq!(
                        report.transport.endpoint.interval_s,
                        [field.time_s(), field.time_s() + dt]
                    );
                    atmosphere = next_air;
                    field = next_field;
                }
                let discrete =
                    equilibrium + (0.1 - equilibrium) * decay(lambda * dt).powi(steps as i32);
                near(field.moisture()[0], discrete, 2e-11);
                near(
                    atmosphere.relative_humidity().unwrap(),
                    r0 - dry * (discrete - 0.1) / capacity,
                    5e-10,
                );
                near(field.time_s(), duration, 1e-15);
                assert_eq!(field.moisture()[0] > 0.1, rh > 0.5);
                (field.moisture()[0] - exact).abs()
            });
            assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
            assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
            assert!(errors[0] / errors[1] < 4.5, "{errors:?}");
            assert!(errors[1] / errors[2] < 4.5, "{errors:?}");
            near(original_air.relative_humidity().unwrap(), r0, 0.0);
        }
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1]);
    });
}

#[test]
fn g1_finite_cartesian_sdirk_weights_six_nonlinear_face_rates_once() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.6, 2.0);
        let before_air = atmosphere.clone();
        let geometry = grid_geometry([0.1, 0.2, 0.05], [1; 3]);
        let initial = grid(card([0.01, 0.04], 0.3), geometry, vec![0.1], cx);
        let slopes = [0.2, 0.25, 0.15, 0.3, 0.18, 0.28];
        let faces = slopes.map(|slope| Some(air::sorption(slope)));
        let geometric: [f64; 6] =
            std::array::from_fn(|side| 2.0 / geometry.lengths_m[side / 2].powi(2));
        let dry = initial.dry_mass_kg();
        let capacity = vapor_capacity(2.0);
        let r0 = atmosphere.relative_humidity().unwrap();
        let dt = 0.003;
        let weighted_slope: f64 = geometric.iter().zip(slopes).map(|(g, s)| g * s).sum();
        // Eliminate finite vapor before solving either stage:
        // F(M)=(.01+.1*M)*(A-B*M). Each implicit stage is one quadratic.
        let a = weighted_slope * (r0 + dry * 0.1 / capacity);
        let b = geometric.iter().sum::<f64>() + weighted_slope * dry / capacity;
        let solve = |base: f64| {
            let h = GAMMA * dt;
            let qa = h * 0.1 * b;
            let qb = 1.0 - h * (0.1 * a - 0.01 * b);
            let positive = base + h * 0.01 * a;
            2.0 * positive / (qb + (qb * qb + 4.0 * qa * positive).sqrt())
        };
        let first = solve(0.1);
        let first_rh = r0 - dry * (first - 0.1) / capacity;
        let first_d = 0.01 + 0.1 * first;
        let endpoint = solve(0.1 + (1.0 - GAMMA) * dt * first_d * (a - b * first));
        let endpoint_rh = r0 - dry * (endpoint - 0.1) / capacity;
        let endpoint_d = 0.01 + 0.1 * endpoint;
        let candidate = atmosphere
            .advance_cartesian_sdirk2(
                &initial,
                dt,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        let report = candidate.step();
        let transport = &report.transport;
        near(candidate.moisture().moisture()[0], endpoint, 2e-11);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            endpoint_rh,
            5e-10,
        );
        assert_eq!(transport.endpoint.interval_s, [0.0, dt]);
        assert_eq!(transport.first_stage.interval_s, [0.0, GAMMA * dt]);
        assert!(transport.first_stage.nonlinear_iterations > 1);
        assert!(transport.endpoint.nonlinear_iterations > 1);
        near(
            transport.first_stage.transport_end[0]
                .property("moisture_diffusivity")
                .unwrap()
                .value_si(),
            first_d,
            3e-12,
        );
        near(
            transport.endpoint.transport_end[0]
                .property("moisture_diffusivity")
                .unwrap()
                .value_si(),
            endpoint_d,
            3e-12,
        );
        near(
            transport.endpoint.transport_start[0].query().axes()["moisture_content"],
            0.1,
            0.0,
        );
        for side in 0..6 {
            let first_rate = dry * geometric[side] * first_d * (slopes[side] * first_rh - first);
            let final_rate =
                dry * geometric[side] * endpoint_d * (slopes[side] * endpoint_rh - endpoint);
            near(
                transport.first_stage.boundary_water_kg[side],
                GAMMA * dt * first_rate,
                2e-12,
            );
            near(
                transport.endpoint.boundary_water_kg[side],
                dt * ((1.0 - GAMMA) * first_rate + GAMMA * final_rate),
                2e-12,
            );
            near(
                transport.endpoint.reservoir_moisture[side].unwrap(),
                slopes[side] * endpoint_rh,
                2e-10,
            );
        }
        assert!(transport.endpoint.boundary_water_kg[2] < 0.0);
        assert!(transport.endpoint.boundary_water_kg[5] > 0.0);
        let water = dry * (endpoint - 0.1);
        near(report.water_to_solid_kg, water, 2e-12);
        near(
            transport.endpoint.boundary_water_kg.iter().sum(),
            water,
            2e-12,
        );
        near(
            report.vapor_mass_kg[0] - report.vapor_mass_kg[1],
            water,
            2e-12,
        );
        near(report.water_balance_residual_kg, 0.0, 1e-13);
        near(transport.endpoint.water_balance_residual_kg, 0.0, 1e-11);
        near(report.maximum_carrier_mass_residual_kg, 0.0, 1e-13);
        near(
            air::carrier_mass(candidate.atmosphere()),
            air::carrier_mass(&atmosphere),
            1e-13,
        );
        near(
            candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
            atmosphere.gas().mass_kg() + initial.wet_mass_kg(),
            2e-13,
        );
        assert_eq!(atmosphere, before_air);
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g1_shared_cartesian_sdirk_uses_one_stage_humidity_in_stable_id_order() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 1.0);
        let ga = grid_geometry([0.1, 0.05, 0.02], [1; 3]);
        let gb = grid_geometry([0.05, 0.02, 0.02], [1; 3]);
        let a = grid(card([0.01; 2], 0.3), ga, vec![0.02], cx);
        let b = grid(card([0.003; 2], 0.3), gb, vec![0.2], cx);
        let mut fa = std::array::from_fn(|_| Some(air::sorption(0.25)));
        fa[0] = Some(emission::law([0.0; 2]));
        let h = 0.05;
        let fb = std::array::from_fn(|_| Some(emission::law([h; 2])));
        let dt = 0.003;
        let rates: [f64; 2] = [
            (1..6)
                .map(|side| 2.0 * 0.01 / ga.lengths_m[side / 2].powi(2))
                .sum(),
            (0..6)
                .map(|side| {
                    let dx = gb.lengths_m[side / 2];
                    1.0 / dx / (1.0 / h + dx / (2.0 * 0.003))
                })
                .sum(),
        ];
        let dry = [a.dry_mass_kg(), b.dry_mass_kg()];
        let old = [0.02, 0.2];
        let r0 = atmosphere.relative_humidity().unwrap();
        let capacity = vapor_capacity(1.0);
        // Eliminate both one-cell fields from the common finite-vapor balance.
        // Stage two retains the ORIGINAL inventories plus its known first rate.
        let solve = |base: [f64; 2]| {
            let alpha = rates.map(|rate| GAMMA * dt * rate);
            let numerator = capacity * r0
                - (0..2)
                    .map(|i| dry[i] * (base[i] / (1.0 + alpha[i]) - old[i]))
                    .sum::<f64>();
            let denominator = capacity
                + (0..2)
                    .map(|i| dry[i] * alpha[i] * 0.25 / (1.0 + alpha[i]))
                    .sum::<f64>();
            let rh = numerator / denominator;
            let content: [f64; 2] =
                std::array::from_fn(|i| (base[i] + alpha[i] * 0.25 * rh) / (1.0 + alpha[i]));
            (content, rh)
        };
        let (first, first_rh) = solve(old);
        let base = std::array::from_fn(|i| {
            old[i] + (1.0 - GAMMA) * dt * rates[i] * (0.25 * first_rh - first[i])
        });
        let (expected, expected_rh) = solve(base);
        let run = |reverse| {
            let participants = if reverse {
                [participant(3, &b, &fb), participant(7, &a, &fa)]
            } else {
                [participant(7, &a, &fa), participant(3, &b, &fb)]
            };
            atmosphere
                .advance_shared_cartesian_sdirk2(&participants, dt, 2, air::controls(), cx)
                .unwrap()
        };
        let accepted = run(false);
        let replay = run(true);
        assert_eq!(
            accepted
                .solids()
                .iter()
                .map(|body| body.id)
                .collect::<Vec<_>>(),
            vec![3, 7]
        );
        assert_eq!(accepted.atmosphere(), replay.atmosphere());
        near(
            accepted.atmosphere().relative_humidity().unwrap(),
            expected_rh,
            2e-10,
        );
        for (body, again) in accepted.solids().iter().zip(replay.solids()) {
            let i = if body.id == 7 { 0 } else { 1 };
            near(body.moisture.moisture()[0], expected[i], 2e-11);
            assert_eq!(body.moisture.moisture(), again.moisture.moisture());
            assert_eq!(body.transport.endpoint.interval_s, [0.0, dt]);
            assert_eq!(body.transport.first_stage.interval_s, [0.0, GAMMA * dt]);
            near(
                body.transport.endpoint.boundary_water_kg.iter().sum(),
                dry[i] * (expected[i] - old[i]),
                2e-12,
            );
        }
        assert!(accepted.solids()[0].moisture.moisture()[0] < 0.2);
        assert!(accepted.solids()[1].moisture.moisture()[0] > 0.02);
        assert_eq!(
            accepted.solids()[1].transport.endpoint.boundary_water_kg[0],
            0.0
        );
        assert_eq!(
            accepted.solids()[1].transport.first_stage.boundary_water_kg[0],
            0.0
        );
        let water: f64 = (0..2).map(|i| dry[i] * (expected[i] - old[i])).sum();
        near(accepted.step().water_to_solids_kg, water, 2e-12);
        near(accepted.step().water_balance_residual_kg, 0.0, 1e-13);
        near(accepted.step().total_mass_residual_kg, 0.0, 1e-13);
        near(
            air::carrier_mass(accepted.atmosphere()),
            air::carrier_mass(&atmosphere),
            1e-13,
        );
        assert_eq!(a.moisture(), &[0.02]);
        assert_eq!(b.moisture(), &[0.2]);
    });
}

#[test]
fn g4_finite_cartesian_sdirk_dormant_faces_cannot_admit_stiff_overshoot() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.8, 1.0);
        let before_air = atmosphere.clone();
        let initial = MoistureCartesian::new(
            directional_card([0.01, 0.0, 0.0]),
            point(0.1),
            ["diffusion_x", "diffusion_y", "diffusion_z"],
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry([0.01; 3], [1; 3]),
            vec![0.1],
            cx,
        )
        .unwrap();
        let mut faces = std::array::from_fn(|_| None);
        faces[0] = Some(air::sorption(0.2));
        faces[1] = Some(emission::law([0.0; 2]));
        faces[2] = Some(air::sorption(0.3));
        // The active face has z approximately 3 and equilibrium near .16.
        // Dormant surface transfer (.20) and impermeable y (.24) must not
        // authorize the numerical overshoot or an extrapolated material query.
        let error = atmosphere
            .advance_cartesian_sdirk2(
                &initial,
                0.015,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap_err();
        assert!(
            matches!(
                &error,
                AtmosphereMoistureError::Moisture(MoistureError::StageRange {
                    moisture, bounds, ..
                }) if *moisture > bounds[1] && bounds[1] < 0.17
            ),
            "{error:?}"
        );
        assert_eq!(atmosphere, before_air);
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
        let smaller = atmosphere
            .advance_cartesian_sdirk2(
                &initial,
                0.0025,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert!(smaller.moisture().moisture()[0] > 0.1);
        assert!(smaller.moisture().moisture()[0] < 0.16);
        assert!(smaller.atmosphere().relative_humidity().unwrap() < 0.8);
    });
}

#[test]
fn g0_finite_cartesian_sdirk_preserves_stiff_equilibrium_and_cancellation() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 1.0);
        let before_air = atmosphere.clone();
        let equilibrium = 0.25 * atmosphere.relative_humidity().unwrap();
        let initial = grid(
            card([0.01; 2], 0.3),
            grid_geometry([0.1; 3], [2; 3]),
            vec![equilibrium; 8],
            cx,
        );
        let faces = std::array::from_fn(|_| Some(air::sorption(0.25)));
        let accepted = atmosphere
            .advance_cartesian_sdirk2(
                &initial,
                1e8,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert_eq!(accepted.moisture().moisture(), initial.moisture());
        assert_eq!(accepted.atmosphere(), &atmosphere);
        assert_eq!(accepted.step().water_to_solid_kg, 0.0);
        assert_eq!(
            accepted.step().transport.endpoint.boundary_water_kg,
            [0.0; 6]
        );
        assert_eq!(
            accepted.step().transport.first_stage.boundary_water_kg,
            [0.0; 6]
        );
        assert_eq!(accepted.step().transport.endpoint.linear_iterations, 0);
        assert_eq!(accepted.step().transport.first_stage.linear_iterations, 0);
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                atmosphere.advance_cartesian_sdirk2(
                    &initial,
                    1.0,
                    &faces,
                    config(),
                    nonlinear(),
                    air::controls(),
                    &cancelled,
                ),
                Err(AtmosphereMoistureError::Moisture(MoistureError::Cancelled))
            ));
        });
        assert_eq!(atmosphere, before_air);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[equilibrium; 8]);
    });
}
