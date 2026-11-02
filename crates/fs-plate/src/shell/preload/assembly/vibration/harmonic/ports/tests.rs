use super::super::AssemblyHarmonicBudget;
use super::*;
use crate::shell::junction::TieOffset;
use crate::shell::preload::assembly::equilibrate_shell_assembly;
use crate::shell::preload::assembly::tests::{Part, budget as tie_budget, options, seam};
use crate::shell::preload::assembly::vibration::tests::{
    dense_reference, modal_options, output_budget, plan, solve, with_cx,
};
use crate::shell::preload::assembly::vibration::{WindowPlan, eigen_windows};
use crate::shell::preload::tests::reference_newton;
use fs_exec::CancelGate;

fn model_budget() -> AssemblyHarmonicBudget {
    AssemblyHarmonicBudget {
        maximum_modes: 16,
        maximum_projection_terms: 100_000_000,
        maximum_projection_bytes: 2 << 20,
    }
}
fn map_budget() -> AssemblyPortProjectionBudget {
    AssemblyPortProjectionBudget {
        maximum_ports: 8,
        maximum_projection_terms: 100_000_000,
        maximum_payload_bytes: 2 << 20,
    }
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 16,
        maximum_payload_bytes: 2 << 20,
        maximum_relative_residual: 1e-10,
    }
}
fn port_controls() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 8,
        maximum_projection_terms: 100_000_000,
        maximum_payload_bytes: 2 << 20,
    }
}
fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 5e-10 * expected.abs().max(1e-100),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_loaded_shell_with_maxwell_port_matches_analytic_resonance_and_original_work() {
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
        attachment[0] = 2.0; // Supported force still has a physical reaction, but does no work.
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
        let b = 1.3 * shape[0][tip];
        assert_eq!(ports.force_map(), &[b]);
        let (m, k) = (model.mass_matrix()[0], model.stiffness_matrix()[0]);
        let natural = (k / m).sqrt();
        let resistance = 0.03 * m * natural / (b * b);
        let spring = 0.4 * k / (b * b);
        let inverse_time = 0.8 * natural;
        let mut loads = vec![C64::ZERO; part.loads.len()];
        loads[tip] = C64::new(0.02, -0.003);
        let force = loads[tip].scale(shape[0][tip]);
        let mean = body.displacements.clone();
        for ratio in [0.8, 1.0, 1.2, 1.5] {
            let omega = ratio * natural;
            let s = C64::new(0.0, -omega);
            // A parallel resistor and Maxwell arm at the declared attachment.
            // This is an analytic fixture law, not an inferred shell material.
            let z =
                C64::from_re(resistance) + C64::from_re(spring) / C64::new(inverse_time, -omega);
            let prepared = ports
                .prepare(omega, &[z], controls(), port_controls(), cx)
                .unwrap();
            let actual = prepared.solve(&[&loads], cx).unwrap();
            let expected = force / (C64::from_re(k - omega * omega * m) + s * z.scale(b * b));
            near(actual.assembly.modal.displacement[0], expected);
            let flow = (s * expected).scale(b);
            let effort = z * flow;
            near(actual.physical_port_flows[0], flow);
            near(
                actual.assembly.modal.port_load.as_ref().unwrap().effort[0],
                effort,
            );
            for (j, &weight) in attachment.iter().enumerate() {
                if weight != 0.0 {
                    near(actual.reaction_forces[0][j], effort.scale(-weight));
                } else {
                    assert_eq!(actual.reaction_forces[0][j], C64::ZERO);
                }
            }
            assert_ne!(actual.reaction_forces[0][0], C64::ZERO);
            assert_eq!(actual.assembly.velocities[0][0], C64::ZERO);
            let heat = 0.5 * z.re * flow.norm_sq();
            near(
                C64::from_re(actual.average_physical_load_power_w),
                C64::from_re(heat),
            );
            near(
                C64::from_re(actual.assembly.average_component_power_w),
                C64::from_re(heat),
            );
            assert_eq!(actual.assembly.modal.average_viscous_power_w, 0.0);
            assert_eq!(body.displacements, mean);
            assert_eq!(
                prepared
                    .solve(&[&loads], cx)
                    .unwrap()
                    .assembly
                    .modal
                    .displacement,
                actual.assembly.modal.displacement
            );
            let wrong_axis =
                force / (C64::from_re(k - omega * omega * m) + s * z.conj().scale(b * b));
            assert!((expected - wrong_axis).abs() > 1e-3 * expected.abs());
        }
    });
}

#[test]
fn g1_mutual_ports_preserve_slave_work_and_g0_g4_refusals_leave_body_unchanged() {
    let left = Part::rectangle(0, 2, 2, 0.0, true);
    let mut right = Part::rectangle(2, 2, 3, 0.0, false);
    right.loads[6 * 8 + 2] = -0.1;
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
                    vec![0.0; right.loads.len()]
                ],
            };
            2
        ];
        shapes[0].component_force_shapes[1][6 * 3 + 2] = 1.0;
        shapes[1].component_force_shapes[0][6 * 8 + 2] = -0.4;
        shapes[1].component_force_shapes[1][right.loads.len() - 4] = 1.0;
        for limited in [
            AssemblyPortProjectionBudget {
                maximum_ports: 1,
                ..map_budget()
            },
            AssemblyPortProjectionBudget {
                maximum_projection_terms: 0,
                ..map_budget()
            },
            AssemblyPortProjectionBudget {
                maximum_payload_bytes: 1,
                ..map_budget()
            },
        ] {
            assert!(AssemblyHarmonicPorts::new(&model, shapes.clone(), limited, cx).is_err());
        }
        let mut malformed = shapes.clone();
        malformed[0].component_force_shapes.pop();
        assert!(AssemblyHarmonicPorts::new(&model, malformed, map_budget(), cx).is_err());
        let mut invalid = shapes.clone();
        invalid[0].component_force_shapes[0][0] = f64::NAN;
        assert!(AssemblyHarmonicPorts::new(&model, invalid, map_budget(), cx).is_err());
        let ports = AssemblyHarmonicPorts::new(&model, shapes, map_budget(), cx).unwrap();
        let omega = 0.94 * vibration.modes().next().unwrap().lambda.sqrt();
        // General passive nonreciprocal load: preserve every complex entry.
        let z = [
            C64::new(1.2, 0.4),
            C64::new(0.3, 0.15),
            C64::new(0.1, 0.05),
            C64::new(0.9, -0.2),
        ];
        let prepared = ports
            .prepare(omega, &z, controls(), port_controls(), cx)
            .unwrap();
        let mut forces = [
            vec![C64::ZERO; left.loads.len()],
            vec![C64::ZERO; right.loads.len()],
        ];
        forces[0][6 * 8 + 2] = C64::new(0.1, -0.03);
        forces[1][right.loads.len() - 4] = C64::new(0.04, 0.02);
        let loads = [&forces[0][..], &forces[1][..]];
        let before = body.displacements.clone();
        let actual = prepared.solve(&loads, cx).unwrap();
        let flow = &actual.physical_port_flows;
        let efforts = [
            z[0] * flow[0] + z[1] * flow[1],
            z[2] * flow[0] + z[3] * flow[1],
        ];
        let work = 0.5 * (flow[0].conj() * efforts[0] + flow[1].conj() * efforts[1]).re;
        assert!(work > 0.0);
        near(
            C64::from_re(actual.average_physical_load_power_w),
            C64::from_re(work),
        );
        near(
            C64::from_re(actual.assembly.average_component_power_w),
            C64::from_re(work),
        );
        for c in 0..forces.len() {
            for j in 0..forces[c].len() {
                let reaction = -efforts[0]
                    .scale(ports.physical_ports()[0].component_force_shapes[c][j])
                    - efforts[1].scale(ports.physical_ports()[1].component_force_shapes[c][j]);
                if reaction == C64::ZERO {
                    assert_eq!(actual.reaction_forces[c][j], C64::ZERO);
                } else {
                    near(actual.reaction_forces[c][j], reaction);
                }
            }
        }
        let diagonal = [z[0], C64::ZERO, C64::ZERO, z[3]];
        let wrong = ports
            .prepare(omega, &diagonal, controls(), port_controls(), cx)
            .unwrap()
            .solve(&loads, cx)
            .unwrap();
        let magnitude = actual
            .assembly
            .modal
            .displacement
            .iter()
            .map(|z| z.abs())
            .fold(0.0, f64::max);
        let difference = actual
            .assembly
            .modal
            .displacement
            .iter()
            .zip(&wrong.assembly.modal.displacement)
            .map(|(&a, &b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(difference > 1e-5 * magnitude);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, |other| {
            assert!(prepared.solve(&loads, other).is_err());
            assert!(
                ports
                    .prepare(omega, &z, controls(), port_controls(), other)
                    .is_err()
            );
        });
        assert_eq!(body.displacements, before);
        assert_eq!(
            prepared
                .solve(&loads, cx)
                .unwrap()
                .assembly
                .modal
                .displacement,
            actual.assembly.modal.displacement
        );
    });
}
