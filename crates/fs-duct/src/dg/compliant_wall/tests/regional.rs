use super::*;

fn base(spatial: bool) -> DgDuct {
    if !spatial {
        return duct(DgBoundary::Closed);
    }
    let field = crate::profile::gas_field::tests::field([300.0, 450.0], [0.95, 0.15], 101325.0);
    DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: field.length_m(),
            }],
        },
        &field,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
}

/// G1: reconstruct the local physical equations and independently integrate
/// storage/loss. Inlet, mean or one-material wall coefficients cannot pass.
#[test]
fn original_element_laws_drive_local_motion_and_physical_energy() {
    let pins = [
        wall(300.0),
        WallPin {
            surface_density: 0.7,
            stiffness_per_area: 4e7,
            resistance: 80.0,
        },
    ];
    for spatial in [false, true] {
        let dg = base(spatial).with_compliant_wall_regions(&pins).unwrap();
        assert!(dg.compliant_wall().is_none());
        assert_eq!(dg.compliant_wall_at_element(0), Some(pins[0]));
        assert_eq!(dg.compliant_wall_at_element(1), Some(pins[1]));
        assert_eq!(dg.compliant_wall_at_element(2), None);
        let m = dg.field_dofs();
        let n = dg.total_dofs();
        let mut state = dg.project_initial_state(|x| 10.0 + x, |_| 0.0);
        for (i, value) in state.tilde[2 * m..].iter_mut().enumerate() {
            *value = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
        }
        let mut derivative = state.clone();
        for (i, row) in dg.l_tilde.chunks_exact(n).enumerate() {
            derivative.tilde[i] = row.iter().zip(&state.tilde).map(|(a, x)| a * x).sum();
        }
        let motion = dg.compliant_wall_state(&state).unwrap();
        let rate = dg.compliant_wall_state(&derivative).unwrap();
        let (pressure, _) = dg.from_energy_coords(&state.tilde);
        let (pressure_rate, _) = dg.from_energy_coords(&derivative.tilde);
        let local = dg.mesh.order + 1;
        for i in 0..m {
            let pin = pins[i / local];
            close(rate.displacement_m[i], motion.velocity_m_s[i], 1.0);
            let reaction = pressure[i]
                - pin.stiffness_per_area * motion.displacement_m[i]
                - pin.resistance * motion.velocity_m_s[i];
            close(
                pin.surface_density * rate.velocity_m_s[i],
                reaction,
                reaction.abs().max(1.0),
            );
        }
        let (points, weights) = gauss_legendre_rule(32);
        let perimeter = core::f64::consts::TAU * RADIUS;
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        let (mut energy, mut loss) = (0.0, 0.0);
        for (e, element) in dg.elements.iter().enumerate() {
            let pin = pins[e];
            let mut residual = vec![0.0; local];
            let mut scale = vec![0.0; local];
            for (&point, &weight) in points.iter().zip(&weights) {
                let shape = dg.lgl.eval_basis(point);
                let interp = |values: &[f64]| -> f64 {
                    (0..local).map(|i| shape[i] * values[e * local + i]).sum()
                };
                let x = element.x_left + 0.5 * (point + 1.0) * element.length;
                let gas = dg
                    .gas_field()
                    .map(|field| *field.state_at(x).unwrap().state())
                    .unwrap_or(*dg.gas());
                let bulk = gas.density * gas.sound_speed.powi(2);
                let w = interp(&motion.displacement_m);
                let v = interp(&motion.velocity_m_s);
                let p_dot = interp(&pressure_rate);
                let dx = 0.5 * element.length * weight;
                energy += 0.5
                    * dx
                    * perimeter
                    * (pin.surface_density * v * v + pin.stiffness_per_area * w * w);
                loss += dx * perimeter * pin.resistance * v * v;
                for i in 0..local {
                    residual[i] += dx * shape[i] * (area / bulk * p_dot + perimeter * v);
                    scale[i] +=
                        dx * shape[i].abs() * ((area / bulk * p_dot).abs() + (perimeter * v).abs());
                }
            }
            for i in 0..local {
                close(residual[i], 0.0, scale[i]);
            }
        }
        close(motion.stored_energy_j, energy, energy);
        close(motion.dissipated_power_w, loss, loss);
        let storage_rate: f64 = state
            .tilde
            .iter()
            .zip(&derivative.tilde)
            .map(|(x, d)| x * d)
            .sum();
        close(storage_rate, -loss, loss);
    }
}

#[test]
fn regional_wall_keeps_gas_memories_and_conjugate_ports_in_the_same_step() {
    let gas = duct(DgBoundary::Port(PortId(0)))
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap();
    let original = gas.total_dofs();
    let dg = gas
        .clone()
        .with_compliant_wall_regions(&[
            wall(300.0),
            WallPin {
                surface_density: 0.7,
                stiffness_per_area: 4e7,
                resistance: 80.0,
            },
        ])
        .unwrap();
    let n = dg.total_dofs();
    for i in 0..original {
        assert_eq!(
            &dg.l_tilde[i * n..i * n + original],
            &gas.l_tilde[i * original..(i + 1) * original]
        );
    }
    assert_eq!(dg.ports[0].0, gas.ports[0].0);
    assert_eq!(&dg.ports[0].1[..original], gas.ports[0].1.as_slice());
    let mut current = dg.zero_state();
    let propagator = dg.propagator(2e-6, 1e-13).unwrap();
    let mut work = 0.0;
    for _ in 0..60 {
        let frame = propagator
            .step(
                &mut current,
                &[PortInput {
                    id: PortId(0),
                    effort: 1e-7,
                }],
            )
            .unwrap();
        work += frame.ledger.port_work_step_j;
    }
    assert!(current.stored_energy_j() > 0.0 && current.stored_energy_j() < work);
    assert!(dg.compliant_wall_state(&current).unwrap().stored_energy_j > 0.0);
    assert!(dg.compliant_wall_state(&gas.zero_state()).is_err());
}

#[test]
fn uniform_regions_keep_original_identity_and_invalid_inputs_do_not_bind() {
    for spatial in [false, true] {
        let gas = base(spatial);
        let original = gas.operator_identity();
        let direct = gas.clone().with_compliant_wall(wall(300.0)).unwrap();
        let same = gas
            .clone()
            .with_compliant_wall_regions(&[wall(300.0); 2])
            .unwrap();
        assert_eq!(same.operator_identity(), direct.operator_identity());
        assert_eq!(same.l_tilde, direct.l_tilde);
        for pins in [
            vec![],
            vec![wall(300.0)],
            vec![wall(300.0); 3],
            vec![wall(300.0), wall(-1.0)],
        ] {
            assert!(gas.clone().with_compliant_wall_regions(&pins).is_err());
        }
        for stop in [1, 4, 9, 15] {
            let mut polls = 0;
            let result = gas.clone().with_compliant_wall_regions_controlled(
                &[
                    wall(300.0),
                    WallPin {
                        surface_density: 0.7,
                        ..wall(50.0)
                    },
                ],
                &mut || {
                    polls += 1;
                    polls == stop
                },
            );
            assert!(matches!(result, Err(DuctError::Cancelled { .. })));
            assert_eq!(gas.operator_identity(), original);
            assert!(gas.compliant_wall_at_element(0).is_none());
        }
    }
}
