use super::*;
use fs_modal::harmonic::HarmonicFrequencyControls;
use fs_modal::harmonic::resonance::{
    ResonanceControls, TransferExtremumKind, refine_transfer_extremum,
};

fn frequency_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
        maximum_additional_payload_bytes: 2_000_000,
        maximum_work: 100_000_000,
    }
}

fn frequency_budget() -> AssemblyPortFrequencyBudget {
    AssemblyPortFrequencyBudget {
        maximum_payload_bytes: 2_000_000,
        maximum_observation_terms: 100_000_000,
    }
}

#[test]
fn maxwell_loaded_shell_tangent_and_physical_receiver_resonances_match_analytic_response() {
    let mut part = Part::rectangle(0, 2, 2, 0.0, true);
    for (node, force) in [(2, 25.0), (5, 50.0), (8, 25.0)] {
        part.loads[6 * node] = force;
    }
    let tip = 6 * 8 + 2;
    part.loads[tip] = -0.01;
    let body = solve(&[part.input()]);
    let before = body.displacements.clone();
    let reference = dense_reference(&body);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration = eigen_windows(
            &body,
            &WindowPlan {
                windows: vec![(0.0, 0.5 * (reference[0].lambda + reference[1].lambda))],
                threads: 1,
                maximum_windows: 1,
                maximum_bytes: 64 << 20,
                maximum_relative_halfwidth: 1e-8,
                coordinates: None,
            },
            &modal_options(),
            output_budget(),
            cx,
        )
        .unwrap();
        let model = AssemblyHarmonicModel::new(&vibration, model_budget(), cx).unwrap();
        let shape_tip = vibration.expand_mode(0, cx).unwrap()[0][tip];
        let mut weights = vec![0.0; part.loads.len()];
        weights[0] = 2.0;
        weights[tip] = 1.3;
        let ports = AssemblyHarmonicPorts::new(
            &model,
            vec![AssemblyHarmonicPort {
                component_force_shapes: vec![weights],
            }],
            map_budget(),
            cx,
        )
        .unwrap();
        let (m, k, b) = (
            model.mass_matrix()[0],
            model.stiffness_matrix()[0],
            ports.force_map()[0],
        );
        let natural = (k / m).sqrt();
        let resistance = 0.03 * m * natural / (b * b);
        let spring = 0.4 * k / (b * b);
        let inverse_time = 0.8 * natural;
        let impedance = |omega: f64| {
            let denominator = C64::new(inverse_time, -omega);
            (
                C64::from_re(resistance) + C64::from_re(spring) / denominator,
                C64::new(0.0, spring) / (denominator * denominator),
            )
        };
        let source = C64::new(0.02, -0.003);
        let source_slope = C64::new(-0.004, 0.007).scale(1.0 / natural);
        let mut force = vec![C64::ZERO; part.loads.len()];
        let mut slope = vec![C64::ZERO; part.loads.len()];
        slope[tip] = source_slope;
        let mut receiver = vec![0.0; part.loads.len()];
        receiver[tip] = 1.0;
        for ratio in [0.0, 0.8, 1.0, 1.2, 1.5] {
            let omega = ratio * natural;
            let (z, dz) = impedance(omega);
            force[tip] = source + source_slope.scale(omega);
            let prepared = ports
                .prepare(omega, &[z], controls(), port_controls(), cx)
                .unwrap();
            let actual = prepared
                .solve_with_frequency_derivative(
                    &[&force],
                    Some(&[&slope]),
                    &[dz],
                    frequency_controls(),
                    frequency_budget(),
                    cx,
                )
                .unwrap();
            let s = C64::new(0.0, -omega);
            let a = C64::from_re(k - omega * omega * m) + (s * z).scale(b * b);
            let da = C64::from_re(-2.0 * omega * m)
                + (C64::new(0.0, -1.0) * (z + dz.scale(omega))).scale(b * b);
            let q = force[tip].scale(shape_tip) / a;
            let dq = (source_slope.scale(shape_tip) - da * q) / a;
            let velocity = s * q;
            let dv = C64::new(0.0, -1.0) * q + s * dq;
            let flow = velocity.scale(b);
            let dflow = dv.scale(b);
            let effort = z * flow;
            let de = dz * flow + z * dflow;
            near(actual.assembly.modal.response.displacement[0], q);
            near(actual.assembly.modal.displacement_derivative[0], dq);
            near(actual.physical_port_flow_derivatives[0], dflow);
            near(actual.reaction_force_derivatives[0][tip], de.scale(-1.3));
            near(actual.reaction_force_derivatives[0][0], de.scale(-2.0));
            assert_eq!(actual.assembly.velocity_derivatives[0][0], C64::ZERO);
            let displacement = actual
                .assembly
                .component_displacement_transfer(0, &receiver)
                .unwrap();
            let speed = actual
                .assembly
                .component_velocity_transfer(0, &receiver)
                .unwrap();
            near(displacement.amplitude, q.scale(shape_tip));
            near(displacement.angular_derivative, dq.scale(shape_tip));
            near(speed.amplitude, velocity.scale(shape_tip));
            near(speed.angular_derivative, dv.scale(shape_tip));
            near(
                actual.port_flow_transfer(0).unwrap().angular_derivative,
                dflow,
            );
            near(
                actual.port_effort_transfer(0).unwrap().angular_derivative,
                de,
            );
            let load_slope = 0.5 * (dflow.conj() * effort + flow.conj() * de).re;
            if omega == 0.0 {
                assert_eq!(actual.average_physical_load_power_derivative, 0.0);
            } else {
                near(
                    C64::from_re(actual.average_physical_load_power_derivative),
                    C64::from_re(load_slope),
                );
            }
            assert!(actual.port_flow_transfer(1).is_err());
            assert!(
                actual
                    .assembly
                    .component_displacement_transfer(1, &receiver)
                    .is_err()
            );
            assert!(
                actual
                    .assembly
                    .component_velocity_transfer(0, &receiver[..tip])
                    .is_err()
            );
            assert!(
                prepared
                    .solve_with_frequency_derivative(
                        &[&force],
                        None,
                        &[],
                        frequency_controls(),
                        frequency_budget(),
                        cx,
                    )
                    .is_err()
            );
            for budget in [
                AssemblyPortFrequencyBudget {
                    maximum_payload_bytes: 1,
                    ..frequency_budget()
                },
                AssemblyPortFrequencyBudget {
                    maximum_observation_terms: 0,
                    ..frequency_budget()
                },
            ] {
                assert!(
                    prepared
                        .solve_with_frequency_derivative(
                            &[&force],
                            None,
                            &[dz],
                            frequency_controls(),
                            budget,
                            cx,
                        )
                        .is_err()
                );
            }
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            with_cx(&cancelled, |other| {
                assert!(
                    prepared
                        .solve_with_frequency_derivative(
                            &[&force],
                            None,
                            &[dz],
                            frequency_controls(),
                            frequency_budget(),
                            other,
                        )
                        .is_err()
                )
            });
            assert_eq!(
                prepared
                    .solve(&[&force], cx)
                    .unwrap()
                    .assembly
                    .modal
                    .displacement,
                actual.assembly.modal.response.displacement
            );
        }
        force[tip] = source;
        let refinement = ResonanceControls {
            absolute_angular_tolerance_rad_s: 1e-12 * natural,
            relative_angular_tolerance: 1e-12,
            maximum_evaluations: 100,
        };
        for (observe_flow, expected_ratio) in [
            (false, 1.126_298_086_627_302_5),
            (true, 1.136_193_386_021_996_6),
        ] {
            let peak =
                refine_transfer_extremum((0.9 * natural, 1.4 * natural), refinement, cx, |omega| {
                    let (z, dz) = impedance(omega);
                    let response = ports
                        .prepare(omega, &[z], controls(), port_controls(), cx)?
                        .solve_with_frequency_derivative(
                            &[&force],
                            None,
                            &[dz],
                            frequency_controls(),
                            frequency_budget(),
                            cx,
                        )?;
                    if observe_flow {
                        response.port_flow_transfer(0)
                    } else {
                        response
                            .assembly
                            .component_displacement_transfer(0, &receiver)
                    }
                })
                .unwrap();
            assert_eq!(peak.kind, TransferExtremumKind::Maximum);
            assert!((peak.angular_frequency_rad_s / natural - expected_ratio).abs() < 1e-9);
        }
        assert_eq!(body.displacements, before);
    });
}

#[test]
fn mutual_port_and_phased_source_tangents_preserve_nonmatching_slave_work() {
    let left = Part::rectangle(0, 2, 2, 0.0, true);
    let mut right = Part::rectangle(2, 2, 3, 0.0, false);
    let last = right.loads.len() - 4;
    right.loads[last] = -0.1;
    let body = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[seam(3, TieOffset::Midsurface)],
        &[],
        tie_budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration = eigen_windows(
            &body,
            &plan(&body, 4),
            &modal_options(),
            output_budget(),
            cx,
        )
        .unwrap();
        let model = AssemblyHarmonicModel::new(&vibration, model_budget(), cx).unwrap();
        let mut shapes = vec![
            AssemblyHarmonicPort {
                component_force_shapes: vec![
                    vec![0.0; left.loads.len()],
                    vec![0.0; right.loads.len()],
                ]
            };
            2
        ];
        shapes[0].component_force_shapes[1][6 * 3 + 2] = 1.0;
        shapes[1].component_force_shapes[0][6 * 8 + 2] = -0.4;
        shapes[1].component_force_shapes[1][last] = 1.0;
        let ports = AssemblyHarmonicPorts::new(&model, shapes, map_budget(), cx).unwrap();
        let omega = 0.94 * vibration.modes().next().unwrap().lambda.sqrt();
        let z = [
            C64::new(1.2, 0.4),
            C64::new(0.3, 0.15),
            C64::new(0.1, 0.05),
            C64::new(0.9, -0.2),
        ];
        let dz = z.map(|entry| entry.scale(0.75 / omega));
        let mut forces = [
            vec![C64::ZERO; left.loads.len()],
            vec![C64::ZERO; right.loads.len()],
        ];
        let mut slopes = forces.clone();
        forces[0][6 * 8 + 2] = C64::new(0.1, -0.03);
        forces[1][last] = C64::new(0.04, 0.02);
        slopes[0][6 * 8 + 2] = C64::new(-0.05, 0.08).scale(1.0 / omega);
        slopes[1][6 * 3 + 2] = C64::new(0.02, -0.03).scale(1.0 / omega);
        let loads = [&forces[0][..], &forces[1][..]];
        let source_slopes = [&slopes[0][..], &slopes[1][..]];
        let prepared = ports
            .prepare(omega, &z, controls(), port_controls(), cx)
            .unwrap();
        let actual = prepared
            .solve_with_frequency_derivative(
                &loads,
                Some(&source_slopes),
                &dz,
                frequency_controls(),
                frequency_budget(),
                cx,
            )
            .unwrap();
        let step = 1e-4 * omega;
        let cases: Vec<_> = [-2.0, -1.0, 1.0, 2.0]
            .into_iter()
            .map(|offset| {
                let delta = offset * step;
                let shifted_z: Vec<_> = z
                    .iter()
                    .zip(dz)
                    .map(|(&z, dz)| z + dz.scale(delta))
                    .collect();
                let shifted_forces: Vec<Vec<_>> = forces
                    .iter()
                    .zip(&slopes)
                    .map(|(f, df)| {
                        f.iter()
                            .zip(df)
                            .map(|(&f, &df)| f + df.scale(delta))
                            .collect()
                    })
                    .collect();
                let shifted_loads: Vec<_> = shifted_forces.iter().map(Vec::as_slice).collect();
                ports
                    .prepare(omega + delta, &shifted_z, controls(), port_controls(), cx)
                    .unwrap()
                    .solve(&shifted_loads, cx)
                    .unwrap()
            })
            .collect();
        let difference = |values: [C64; 4]| {
            (values[0] - values[1].scale(8.0) + values[2].scale(8.0) - values[3])
                .scale(1.0 / (12.0 * step))
        };
        let near_slope = |actual: C64, expected: C64| {
            assert!(
                (actual - expected).abs() <= 3e-6 * actual.abs().max(expected.abs()).max(1e-14),
                "{actual:?} != {expected:?}"
            )
        };
        for c in 0..forces.len() {
            for j in 0..forces[c].len() {
                near_slope(
                    actual.assembly.displacement_derivatives[c][j],
                    difference(core::array::from_fn(|i| {
                        cases[i].assembly.perturbation_displacements[c][j]
                    })),
                );
                near_slope(
                    actual.assembly.velocity_derivatives[c][j],
                    difference(core::array::from_fn(|i| cases[i].assembly.velocities[c][j])),
                );
                near_slope(
                    actual.reaction_force_derivatives[c][j],
                    difference(core::array::from_fn(|i| cases[i].reaction_forces[c][j])),
                );
            }
        }
        for port in 0..2 {
            near_slope(
                actual.physical_port_flow_derivatives[port],
                difference(core::array::from_fn(|i| cases[i].physical_port_flows[port])),
            );
        }
        near_slope(
            C64::from_re(actual.average_physical_load_power_derivative),
            difference(core::array::from_fn(|i| {
                C64::from_re(cases[i].average_physical_load_power_w)
            })),
        );
        near_slope(
            C64::from_re(actual.assembly.average_component_power_derivative),
            difference(core::array::from_fn(|i| {
                C64::from_re(cases[i].assembly.average_component_power_w)
            })),
        );
        let fixed_load = prepared
            .solve_with_frequency_derivative(
                &loads,
                Some(&source_slopes),
                &[C64::ZERO; 4],
                frequency_controls(),
                frequency_budget(),
                cx,
            )
            .unwrap();
        assert!(
            actual
                .assembly
                .modal
                .displacement_derivative
                .iter()
                .zip(&fixed_load.assembly.modal.displacement_derivative)
                .any(|(&a, &b)| (a - b).abs() > 1e-5 * a.abs())
        );
    });
}
