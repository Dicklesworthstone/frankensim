//! A material-derived string tail retains the physical fiber and bridge work.
use super::*;

fn omitted_sum(card: &fs_couple::bowed_string::BowedStringCard, x: f64, y: f64) -> f64 {
    // Independent convergent sine series, not the production static Green
    // function or its Schur-complement helper. EI>0 makes this tail O(n^-4).
    (card.mode_count + 1..=20_000)
        .map(|mode| {
            let angle = mode as f64 * core::f64::consts::PI;
            let wave = angle / card.length_m;
            2.0 * (angle * x).sin() * (angle * y).sin()
                / (card.length_m
                    * (card.tension_n * wave * wave + card.bending_stiffness_n_m2 * wave.powi(4)))
        })
        .sum()
}

#[test]
fn receiving_residual_uses_source_string_tail_without_replacing_fiber_or_bridge_work() {
    with_cx(|cx, gate| {
        let make = || {
            base(cx, gate, 7800.0)
                .with_residual_flexibility(gate)
                .unwrap()
                .with_tensioned_fiber_ribbon(
                    fibers(3, 2, 1300.0, 3e9, 50.0, 40.0, 48000, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap()
        };
        let mut reference = make();
        let original_ports: Vec<_> = (0..3)
            .map(|i| reference.system().port(i).unwrap().clone())
            .collect();
        let original_modes: Vec<_> = reference
            .system()
            .components()
            .iter()
            .map(|m| m.modes().to_vec())
            .collect();
        let mut actual = make().with_bow_residual_flexibility(gate).unwrap();
        let card = &actual.cards()[0];
        let string_count = actual.cards().len();
        let member = &actual.source().bank().strings()[0];
        let connection = &member.connections[0];
        let norm = (0.5 * card.linear_density_kg_m * card.length_m).sqrt();
        let bridge: Vec<_> = (1..=card.mode_count)
            .map(|mode| {
                (mode as f64 * core::f64::consts::PI * connection.station_fraction).sin() / norm
            })
            .collect();
        let bridge_component = string_count + connection.body.component;
        let body_shapes = connection.body.shapes.clone();
        let gbb = omitted_sum(
            card,
            connection.station_fraction,
            connection.station_fraction,
        );
        let stiffness = connection.stiffness_n_m / (1.0 + connection.stiffness_n_m * gbb);
        let stations = actual.source().ribbon().unwrap().stations().to_vec();
        let cross: Vec<_> = stations
            .iter()
            .map(|&x| omitted_sum(card, x, connection.station_fraction))
            .collect();
        let tail = actual.system().friction_residual_compliance().unwrap();
        for i in 0..3 {
            assert_eq!(
                actual.system().port(i).unwrap().right,
                original_ports[i].right
            );
            for j in 0..3 {
                let independent =
                    omitted_sum(card, stations[i], stations[j]) - stiffness * cross[i] * cross[j];
                assert!(
                    (tail[3 * i + j] - independent).abs() < 1e-13,
                    "source tail {i},{j}"
                );
            }
            for (&mapped, (&original, &bridge_shape)) in actual
                .system()
                .port(i)
                .unwrap()
                .left
                .shapes
                .iter()
                .zip(original_ports[i].left.shapes.iter().zip(&bridge))
            {
                assert!((mapped - original + stiffness * cross[i] * bridge_shape).abs() < 1e-9);
            }
        }
        for (model, original) in actual.system().components().iter().zip(original_modes) {
            assert_eq!(
                model.modes(),
                original,
                "condensation must not retune retained string, body or fiber modes"
            );
        }
        let zero = vec![0.0; actual.system().mode_count()];
        let frame = actual
            .step(
                &zero,
                &[drive(0.001, 0.1), drive(-0.0005, 0.1), drive(0.0015, 0.1)],
            )
            .unwrap()
            .clone();
        let mut applied = zero;
        let mut offsets = Vec::new();
        let mut offset = 0;
        for model in reference.system().components() {
            offsets.push(offset);
            offset += model.modes().len();
        }
        for i in 0..3 {
            let force = frame.ports[i].traction_n;
            assert!(force != 0.0);
            let factor = stiffness * cross[i];
            let left = &original_ports[i].left;
            for (k, (&original, &bridge_shape)) in left.shapes.iter().zip(&bridge).enumerate() {
                applied[offsets[left.component] + k] += force * (original - factor * bridge_shape);
            }
            let right = original_ports[i].right.as_ref().unwrap();
            for (k, &shape) in right.shapes.iter().enumerate() {
                applied[offsets[right.component] + k] -= force * shape;
            }
            for (k, &shape) in body_shapes.iter().enumerate() {
                applied[offsets[bridge_component] + k] += force * factor * shape;
            }
        }
        // Apply the independently reconstructed virtual-work forces to the
        // identical original network with every contact released. This checks
        // the actual retained body/fiber states, not only a scalar energy sum.
        let expected = reference
            .step(&applied, &[drive(0.0, 0.0); 3])
            .unwrap()
            .clone();
        for (a, b) in states(&actual).iter().zip(states(&reference)) {
            assert!((a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg).abs() < 1e-13);
            assert!((a.velocity_m_sqrt_kg_per_s - b.velocity_m_sqrt_kg_per_s).abs() < 1e-9);
        }
        assert!(frame.observer_pressure_pa != 0.0);
        assert!((frame.observer_pressure_pa - expected.observer_pressure_pa).abs() < 1e-8);
        assert!(frame.friction_residual_energy_j > 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(
            make()
                .with_bow_residual_flexibility(gate)
                .unwrap()
                .with_fiber_support_motion(0.0, support_limits(), gate)
                .is_err(),
            "a fixed omitted tail cannot silently become a moving contact"
        );
        let nonlinear = fibers(3, 2, 1300.0, 3e9, 50.0, 40.0, 48000, gate)
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .unwrap();
        let nonlinear = base(cx, gate, 7800.0)
            .with_residual_flexibility(gate)
            .unwrap()
            .with_tensioned_fiber_ribbon(nonlinear, 0.4, ribbon_limits(), gate)
            .unwrap();
        assert!(
            nonlinear.with_bow_residual_flexibility(gate).is_err(),
            "the linear omitted tail cannot silently ignore nonlinear fiber strain"
        );
    });
}
