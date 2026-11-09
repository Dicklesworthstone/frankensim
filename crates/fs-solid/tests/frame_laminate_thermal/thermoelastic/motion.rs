//! Native motion and temperature derivatives against compatible frame states.

use super::*;
use fs_solid::Rayleigh;
use fs_solid::frame2d::laminate_thermal::thermoelastic::motion::FrameLaminateMotionHeatRate;

const UNDAMPED: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };
const STEP: f64 = 0.001;

fn shifted(temperatures: &[Vec<f64>], rates: &[Vec<f64>], step: f64) -> Vec<Vec<f64>> {
    temperatures
        .iter()
        .zip(rates)
        .map(|(temperatures, rates)| {
            temperatures
                .iter()
                .zip(rates)
                .map(|(&t, &rate)| t + step * rate)
                .collect()
        })
        .collect()
}

fn close_rate(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-6 * expected.abs().max(scale.abs()).max(1e-30),
        "native derivative {actual:e} != finite-difference derivative {expected:e}"
    );
}

fn velocity(frame: &mut Frame2d, target: &[f64], load: &[f64]) {
    // The original linear mass action maps a velocity to generalized momentum
    // just as it maps acceleration to force. Its actual impulse owner installs
    // the rate and resolves consistent supported acceleration; no state setter.
    let momentum = frame.inertial_forces(target).unwrap();
    frame
        .apply_impulse(&momentum, load, UNDAMPED, 0.01)
        .unwrap();
    for (&actual, &expected) in frame.velocities().iter().zip(target) {
        near(
            actual,
            expected,
            target.iter().fold(0.0_f64, |m, x| m.max(x.abs())),
        );
    }
}

fn tip_load(angle: f64, axial: f64, transverse: f64, moment: f64) -> [f64; 6] {
    [
        0.0,
        0.0,
        0.0,
        axial * angle.cos() - transverse * angle.sin(),
        axial * angle.sin() + transverse * angle.cos(),
        moment,
    ]
}

#[test]
fn constrained_spatial_heating_has_compatible_internal_redistribution_and_native_force_rates() {
    with_cx(|cx| {
        let temperatures = vec![vec![310.0, 330.0], vec![330.0, 360.0], vec![350.0, 380.0]];
        let rates = vec![vec![1.2, -0.7], vec![-2.0, 1.0], vec![0.3, 2.2]];
        let mut frame = member(
            typed_source(&bimetal(), 0.0, false),
            0.0005,
            0.1,
            0.63,
            &vec![vec![300.0; 2]; 3],
            cx,
        );
        frame.fix(1, [true; 3]);
        let frame = frame
            .with_laminate_thermoelasticity(0, &[[260.0, 400.0]; 2], cx)
            .unwrap()
            .with_laminate_temperatures(0, &temperatures, cx)
            .unwrap();
        let response = frame
            .laminate_thermoelastic_motion_heat_rate(0, &rates, cx)
            .unwrap();
        assert_eq!(response.basic_deformation_rate, [0.0; 3]);
        assert_eq!(frame.velocities(), &[0.0; 6]);
        let plus = frame
            .with_laminate_temperatures(0, &shifted(&temperatures, &rates, STEP), cx)
            .unwrap();
        let minus = frame
            .with_laminate_temperatures(0, &shifted(&temperatures, &rates, -STEP), cx)
            .unwrap();
        assert_eq!(plus.displacements(), frame.displacements());
        assert_eq!(minus.displacements(), frame.displacements());
        let mut compatible = [0.0; 3];
        let mut internal_rate = 0.0_f64;
        for station in &response.stations {
            let i = station.station;
            let before = minus.section_deformations(0)[i];
            let after = plus.section_deformations(0)[i];
            let expected = [
                (after.0 - before.0) / (2.0 * STEP),
                (after.1 - before.1) / (2.0 * STEP),
            ];
            for axis in 0..2 {
                close_rate(
                    station.generalized_strain_rate[axis],
                    expected[axis],
                    expected[axis],
                );
            }
            let [eps, kappa] = station.generalized_strain_rate;
            let weight = station.reference_length_weight_m;
            let xi = station.station_coordinate;
            compatible[0] += weight * eps;
            compatible[1] += weight * (xi - 1.0) * kappa;
            compatible[2] += weight * xi * kappa;
            internal_rate = internal_rate.max(eps.abs()).max(kappa.abs());
        }
        assert!(internal_rate > 1e-4);
        for residual in compatible {
            near(residual, 0.0, 0.1 * internal_rate);
        }
        for axis in 0..3 {
            let expected =
                (plus.element_forces(0)[axis] - minus.element_forces(0)[axis]) / (2.0 * STEP);
            close_rate(response.basic_force_rate[axis], expected, expected);
        }
    });
}

#[test]
fn actual_native_motion_and_source_temperature_rates_follow_the_same_equilibrium_derivative() {
    with_cx(|cx| {
        let plies = bimetal();
        let width = 0.0005;
        let temperatures = vec![vec![310.0, 330.0], vec![330.0, 350.0], vec![350.0, 370.0]];
        let rates = vec![vec![1.0, -0.5], vec![2.0, 0.7], vec![-0.3, 1.4]];
        let mut baseline: Option<FrameLaminateMotionHeatRate> = None;
        for angle in [0.0, 0.73] {
            let load = tip_load(angle, 0.5, 0.02, 0.001);
            let load_rate = tip_load(angle, 0.7, -0.01, 0.0008);
            let mut frame = member(
                typed_source(&plies, 0.0, false),
                width,
                0.1,
                angle,
                &temperatures,
                cx,
            )
            .with_laminate_thermoelasticity(0, &[[260.0, 400.0]; 2], cx)
            .unwrap();
            frame.static_load(&load, 1).unwrap();
            let mut plus = frame
                .with_laminate_temperatures(0, &shifted(&temperatures, &rates, STEP), cx)
                .unwrap();
            let mut minus = frame
                .with_laminate_temperatures(0, &shifted(&temperatures, &rates, -STEP), cx)
                .unwrap();
            let plus_load: Vec<_> = load
                .iter()
                .zip(load_rate)
                .map(|(&p, rate)| p + STEP * rate)
                .collect();
            let minus_load: Vec<_> = load
                .iter()
                .zip(load_rate)
                .map(|(&p, rate)| p - STEP * rate)
                .collect();
            plus.static_load(&plus_load, 1).unwrap();
            minus.static_load(&minus_load, 1).unwrap();
            let target: Vec<_> = plus
                .displacements()
                .iter()
                .zip(minus.displacements())
                .map(|(&after, &before)| (after - before) / (2.0 * STEP))
                .collect();
            assert!(target.iter().any(|v| v.abs() > 1e-5));
            let pose = frame.displacements().to_vec();
            velocity(&mut frame, &target, &load);
            let response = frame
                .laminate_thermoelastic_motion_heat_rate(0, &rates, cx)
                .unwrap();
            assert_eq!(frame.displacements(), pose);
            for station in &response.stations {
                let i = station.station;
                let before = minus.section_deformations(0)[i];
                let after = plus.section_deformations(0)[i];
                let expected = [
                    (after.0 - before.0) / (2.0 * STEP),
                    (after.1 - before.1) / (2.0 * STEP),
                ];
                for axis in 0..2 {
                    close_rate(
                        station.generalized_strain_rate[axis],
                        expected[axis],
                        expected[axis],
                    );
                }
                let (eps, kappa) = frame.section_deformations(0)[i];
                let (_, _, faces) = oracle(&plies, &temperatures[i], width);
                for layer in 0..2 {
                    let mid = 0.5 * (faces[layer][0] + faces[layer][1]);
                    let source =
                        potential(plies[layer], temperatures[i][layer], width, mid, eps, kappa);
                    let deformation = temperatures[i][layer]
                        * (source[6] * expected[0] + source[7] * expected[1]);
                    let thermal = -source[5] * rates[i][layer];
                    close_rate(
                        station.layers[layer].deformation_heat_w_m,
                        deformation,
                        deformation,
                    );
                    near(station.layers[layer].temperature_heat_w_m, thermal, thermal);
                }
            }
            for axis in 0..3 {
                let expected =
                    (plus.element_forces(0)[axis] - minus.element_forces(0)[axis]) / (2.0 * STEP);
                close_rate(response.basic_force_rate[axis], expected, expected);
            }
            if let Some(old) = &baseline {
                for (a, b) in response.stations.iter().zip(&old.stations) {
                    for axis in 0..2 {
                        close_rate(
                            a.generalized_strain_rate[axis],
                            b.generalized_strain_rate[axis],
                            b.generalized_strain_rate[axis],
                        );
                    }
                    for (a, b) in a.layers.iter().zip(&b.layers) {
                        close_rate(
                            a.reversible_heat_w_m,
                            b.reversible_heat_w_m,
                            b.reversible_heat_w_m,
                        );
                    }
                }
            }
            baseline = Some(response);
        }
    });
}

#[test]
fn actual_free_expansion_velocity_cancels_reversible_heat_for_both_cte_signs() {
    with_cx(|cx| {
        for sign in [1.0, -1.0] {
            let p = Ply {
                h: 0.005,
                rho: 4300.0,
                young: 2e9,
                poisson: 0.31,
                alpha: sign * 1e-5,
                softening: 0.4,
            };
            let mut frame = member(
                typed_source(&[p], 0.0, false),
                0.0005,
                0.1,
                0.0,
                &vec![vec![350.0]; 3],
                cx,
            )
            .with_laminate_thermoelasticity(0, &[[260.0, 400.0]], cx)
            .unwrap();
            frame.static_load(&[0.0; 6], 1).unwrap();
            let strain_rate = 1.5 * p.alpha * 4.0;
            velocity(
                &mut frame,
                &[0.0, 0.0, 0.0, 0.1 * strain_rate, 0.0, 0.0],
                &[0.0; 6],
            );
            let response = frame
                .laminate_thermoelastic_motion_heat_rate(0, &vec![vec![4.0]; 3], cx)
                .unwrap();
            for station in response.stations {
                near(station.generalized_strain_rate[0], strain_rate, strain_rate);
                near(station.generalized_strain_rate[1], 0.0, strain_rate / p.h);
                let heat = station.layers[0];
                assert!(heat.deformation_heat_w_m.abs() > 0.0);
                assert!(heat.elastic_heat_capacity_j_m_k < 0.0);
                near(heat.reversible_heat_w_m, 0.0, heat.deformation_heat_w_m);
            }
        }
    });
}

#[test]
fn native_rigid_motion_has_no_strain_heat_and_unsupported_context_refuses() {
    with_cx(|cx| {
        let length: f64 = 0.1;
        let angle: f64 = 0.63;
        let end = [length * angle.cos(), length * angle.sin()];
        for geometry in [Geometry::Linear, Geometry::Corotational] {
            let mut frame =
                Frame2d::with_mass_formulation(geometry, FrameMassFormulation::ConsistentReference);
            frame.add_node(0.0, 0.0);
            frame.add_node(end[0], end[1]);
            frame
                .add_laminate_thermal_element(
                    0,
                    1,
                    typed_source(&bimetal(), 0.0, false),
                    0.0005,
                    &vec![vec![300.0; 2]; 3],
                    0.02,
                    cx,
                )
                .unwrap();
            let mut frame = frame
                .with_laminate_thermoelasticity(0, &[[260.0, 400.0]; 2], cx)
                .unwrap();
            let zero_rates = vec![vec![0.0; 2]; 3];
            if geometry == Geometry::Corotational {
                assert!(
                    frame
                        .laminate_thermoelastic_motion_heat_rate(0, &zero_rates, cx)
                        .is_err()
                );
                continue;
            }
            let omega = 0.04;
            let target = [
                0.2,
                -0.3,
                omega,
                0.2 - omega * end[1],
                -0.3 + omega * end[0],
                omega,
            ];
            velocity(&mut frame, &target, &[0.0; 6]);
            let before = frame.velocities().to_vec();
            let result = frame.laminate_thermoelastic_motion_heat_rate(0, &zero_rates, cx);
            let linear_speed = before.chunks_exact(3).fold(0.0_f64, |speed, node| {
                speed.max(node[0].abs()).max(node[1].abs())
            });
            let angular_speed = before
                .chunks_exact(3)
                .fold(0.0_f64, |speed, node| speed.max(node[2].abs()));
            // Allow roundoff from the actual six-coordinate mass solve and
            // geometric/tangent actions. Curvature carries one further 1/L;
            // neither dimensional rate receives an arbitrary absolute floor.
            let roundoff = 8192.0 * f64::EPSILON;
            let strain_allowance = roundoff * (linear_speed / length + angular_speed);
            let curvature_allowance =
                roundoff * (linear_speed / (length * length) + angular_speed / length);
            for station in result.unwrap().stations {
                assert!(station.generalized_strain_rate[0].abs() <= strain_allowance);
                assert!(station.generalized_strain_rate[1].abs() <= curvature_allowance);
                let source = frame
                    .laminate_thermoelastic_response(0, station.station, cx)
                    .unwrap();
                for layer in station.layers {
                    let [n_t, m_t] = source.layers[layer.layer]
                        .source
                        .generalized_force_temperature_derivative;
                    let heat_allowance = layer.temperature_k
                        * (n_t.abs() * strain_allowance + m_t.abs() * curvature_allowance)
                        * (1.0 + 16.0 * f64::EPSILON);
                    assert!(layer.reversible_heat_w_m.abs() <= heat_allowance);
                    assert_eq!(layer.temperature_heat_w_m, 0.0);
                }
            }
            assert_eq!(frame.velocities(), before);
            assert!(
                frame
                    .laminate_thermoelastic_motion_heat_rate(0, &zero_rates[..2], cx)
                    .is_err()
            );
            assert!(
                frame
                    .laminate_thermoelastic_motion_heat_rate(0, &vec![vec![f64::NAN; 2]; 3], cx)
                    .is_err()
            );
            let gate = CancelGate::new_clock_free();
            gate.request();
            ArenaPool::new(ArenaConfig::default()).scope(|arena| {
                let cancelled = Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 17,
                        kernel_id: 21,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                );
                assert!(matches!(
                    frame.laminate_thermoelastic_motion_heat_rate(0, &zero_rates, &cancelled),
                    Err(SolidError::Cancelled)
                ));
            });
            assert_eq!(frame.velocities(), before);
        }
    });
}
