//! G1 same accepted shell, actual attachment and independently differentiated scalar equation.
use super::*;
use fs_modal::harmonic::HarmonicFrequencyCurvatureControls;

#[test]
fn loaded_shell_curvature_retains_source_phase_support_reactions_and_work() {
    let mut part = Part::rectangle(0, 2, 2, 0.0, true);
    for (node, force) in [(2, 25.0), (5, 50.0), (8, 25.0)] {
        part.loads[6 * node] = force;
    }
    part.loads[6 * 8 + 2] = -0.01;
    let body = solve(&[part.input()]);
    let reference = dense_reference(&body);
    let windows = WindowPlan {
        windows: vec![(0.0, 0.5 * (reference[0].lambda + reference[1].lambda))],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 64 << 20,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let model = AssemblyHarmonicModel::new(&vibration, model_budget(), cx).unwrap();
        let shape = vibration.expand_mode(0, cx).unwrap();
        let tip = 6 * 8 + 2;
        let mut attachment = vec![0.0; part.loads.len()];
        attachment[0] = 2.0;
        attachment[tip] = 1.3;
        let ports = AssemblyHarmonicPorts::new(
            &model,
            vec![AssemblyHarmonicPort {
                component_force_shapes: vec![attachment.clone()],
            }],
            map_budget(),
            cx,
        )
        .unwrap();
        let b = ports.force_map()[0];
        let (m, k) = (model.mass_matrix()[0], model.stiffness_matrix()[0]);
        let natural = (k / m).sqrt();
        let resistance = 0.03 * m * natural / (b * b);
        let spring = 0.4 * k / (b * b);
        let rate = 0.8 * natural;
        let curvature_controls = HarmonicFrequencyCurvatureControls {
            maximum_additional_payload_bytes: 4 << 20,
            maximum_work: 100_000_000,
        };
        let budget = AssemblyPortCurvatureBudget {
            maximum_payload_bytes: 4 << 20,
            maximum_observation_terms: 100_000_000,
        };
        let mean = body.displacements.clone();
        for ratio in [0.8, 1.0, 1.2, 1.5] {
            let w = ratio * natural;
            let s = C64::new(0.0, -w);
            let distance = C64::new(rate, -w);
            let z = C64::from_re(resistance) + C64::from_re(spring) / distance;
            let zw = C64::new(0.0, spring) / distance / distance;
            let zww = C64::from_re(-2.0 * spring) / distance / distance / distance;
            let prepared = ports
                .prepare(w, &[z], controls(), port_controls(), cx)
                .unwrap();
            let delay = 0.2 / natural;
            let mut f = vec![C64::ZERO; part.loads.len()];
            let mut fw = f.clone();
            let mut fww = f.clone();
            f[tip] = C64::new(0.02, -0.003) * C64::new((w * delay).cos(), (w * delay).sin());
            fw[tip] = C64::new(0.0, delay) * f[tip];
            fww[tip] = f[tip].scale(-delay * delay);
            let actual = prepared
                .solve_with_frequency_curvature(
                    &[&f],
                    Some(&[&fw]),
                    Some(&[&fww]),
                    &[zw],
                    &[zww],
                    curvature_controls,
                    budget,
                    cx,
                )
                .unwrap();
            let force = f[tip].scale(shape[0][tip]);
            let force_w = fw[tip].scale(shape[0][tip]);
            let force_ww = fww[tip].scale(shape[0][tip]);
            let d = C64::from_re(k - w * w * m) + s * z.scale(b * b);
            let dw = C64::new(-2.0 * w * m, 0.0) + C64::new(0.0, -b * b) * z + s * zw.scale(b * b);
            let dww =
                C64::from_re(-2.0 * m) + C64::new(0.0, -2.0 * b * b) * zw + s * zww.scale(b * b);
            let q = force / d;
            let qw = (force_w - dw * q) / d;
            let qww = (force_ww - dww * q - (dw * qw).scale(2.0)) / d;
            let v = s * q;
            let vw = s * qw + C64::new(q.im, -q.re);
            let vww = s * qww + C64::new(2.0 * qw.im, -2.0 * qw.re);
            let flow = v.scale(b);
            let flow_w = vw.scale(b);
            let flow_ww = vww.scale(b);
            let effort_ww = zww * flow + (zw * flow_w).scale(2.0) + z * flow_ww;
            near(actual.modal_displacement_second_derivative[0], qww);
            near(actual.modal_velocity_second_derivative[0], vww);
            near(actual.physical_port_flow_second_derivatives[0], flow_ww);
            for (j, &weight) in attachment.iter().enumerate() {
                if weight != 0.0 {
                    near(
                        actual.reaction_force_second_derivatives[0][j],
                        effort_ww.scale(-weight),
                    );
                }
                near_or_zero(
                    actual.displacement_second_derivatives[0][j],
                    qww.scale(shape[0][j]),
                );
                near_or_zero(
                    actual.velocity_second_derivatives[0][j],
                    vww.scale(shape[0][j]),
                );
            }
            assert_eq!(actual.velocity_second_derivatives[0][0], C64::ZERO);
            assert_ne!(actual.reaction_force_second_derivatives[0][0], C64::ZERO);
            let expected_power = 0.5
                * ((force_ww.conj() * v).re
                    + 2.0 * (force_w.conj() * vw).re
                    + (force.conj() * vww).re);
            let scale = expected_power.abs().max(1e-100);
            assert!(
                (actual.average_component_power_second_derivative - expected_power).abs()
                    < 1e-9 * scale
            );
            assert!(
                (actual.average_physical_load_power_second_derivative - expected_power).abs()
                    < 1e-9 * scale
            );
            assert_eq!(actual.average_viscous_power_second_derivative, 0.0);
            assert_eq!(
                actual
                    .port_power_curvature_sample()
                    .angular_second_derivative,
                actual.average_physical_load_power_second_derivative
            );
            assert!(
                actual
                    .physical_power_balance_residual_second_derivative
                    .abs()
                    < 1e-9 * scale
            );
            assert!(
                prepared
                    .solve_with_frequency_curvature(
                        &[&f],
                        None,
                        None,
                        &[zw],
                        &[zww],
                        curvature_controls,
                        AssemblyPortCurvatureBudget {
                            maximum_payload_bytes: 0,
                            ..budget
                        },
                        cx,
                    )
                    .is_err()
            );
            assert!(
                prepared
                    .solve_with_frequency_curvature(
                        &[&f],
                        None,
                        None,
                        &[zw],
                        &[],
                        curvature_controls,
                        budget,
                        cx,
                    )
                    .is_err()
            );
            assert_eq!(
                prepared
                    .solve(&[&f], cx)
                    .unwrap()
                    .assembly
                    .modal
                    .displacement,
                actual
                    .frequency_response
                    .assembly
                    .modal
                    .response
                    .displacement
            );
            assert_eq!(body.displacements, mean);
        }
    });
}
fn near_or_zero(actual: C64, expected: C64) {
    if expected == C64::ZERO {
        assert_eq!(actual, expected);
    } else {
        near(actual, expected);
    }
}
