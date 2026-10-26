use super::*;
use crate::profile::gas_field::tests::field;
use fs_material::visco::GeneralizedMaxwell;

fn build(field: &AxialGasField, radius: f64, inlet: DgBoundary) -> DgDuct {
    DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: field.length_m(),
            }],
        },
        field,
        DgMesh {
            elements: 2,
            order: 2,
        },
        inlet,
        DgBoundary::Closed,
    )
    .unwrap()
}

/// Independent physical quadrature checks every pressure, displacement,
/// velocity and Maxwell coordinate against the original local equations.
/// A heterogeneous gas mass cannot pass by using inlet B or a mean B.
#[test]
fn pointwise_gas_and_geometric_wall_masses_close_the_physical_equations() {
    let field = field([300.0, 450.0], [0.95, 0.15], 101325.0);
    let pin = wall(300.0);
    let law = GeneralizedMaxwell::new(1e6, vec![(8e5, 0.003), (4e5, 0.0005)]).unwrap();
    let dg = build(&field, RADIUS, DgBoundary::Closed)
        .with_compliant_wall_relaxation(pin, &law)
        .unwrap();
    let m = dg.field_dofs();
    let n = dg.total_dofs();
    let mut state = dg.project_initial_state(|x| 10.0 + x, |_| 0.0);
    for (i, x) in state.tilde[2 * m..].iter_mut().enumerate() {
        *x = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
    }
    let mut derivative = state.clone();
    for i in 0..n {
        derivative.tilde[i] = dg.l_tilde[i * n..(i + 1) * n]
            .iter()
            .zip(&state.tilde)
            .map(|(a, x)| a * x)
            .sum();
    }
    let motion = dg.compliant_wall_state(&state).unwrap();
    let rate = dg.compliant_wall_state(&derivative).unwrap();
    let (pressure, _) = dg.from_energy_coords(&state.tilde);
    let (pressure_rate, _) = dg.from_energy_coords(&derivative.tilde);
    for i in 0..m {
        close(rate.displacement_m[i], motion.velocity_m_s[i], 1.0);
        let reaction: f64 = motion.relaxation_pressure_pa.iter().map(|z| z[i]).sum();
        let force = pressure[i]
            - pin.stiffness_per_area * motion.displacement_m[i]
            - pin.resistance * motion.velocity_m_s[i]
            - reaction;
        close(
            pin.surface_density * rate.velocity_m_s[i],
            force,
            force.abs().max(1.0),
        );
        for (j, (modulus, tau)) in law.terms.iter().enumerate() {
            let stiffness = pin.stiffness_per_area * modulus / law.e_inf;
            let expected =
                stiffness * motion.velocity_m_s[i] - motion.relaxation_pressure_pa[j][i] / tau;
            close(
                rate.relaxation_pressure_pa[j][i],
                expected,
                expected.abs().max(1.0),
            );
        }
    }
    let (points, weights) = gauss_legendre_rule(32);
    let perimeter = core::f64::consts::TAU * RADIUS;
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let local = dg.mesh.order + 1;
    let (mut energy, mut loss) = (0.0, 0.0);
    for (e, element) in dg.elements.iter().enumerate() {
        let base = e * local;
        let mut residual = vec![0.0; local];
        let mut scale = vec![0.0; local];
        for (&xi, &weight) in points.iter().zip(&weights) {
            let basis = dg.lgl.eval_basis(xi);
            let interpolate =
                |values: &[f64]| (0..local).map(|i| basis[i] * values[base + i]).sum::<f64>();
            let x = element.x_left + 0.5 * (xi + 1.0) * element.length;
            let gas = field.state_at(x).unwrap();
            let gas = gas.state();
            let bulk = gas.density * gas.sound_speed * gas.sound_speed;
            let v = interpolate(&motion.velocity_m_s);
            let w = interpolate(&motion.displacement_m);
            let p_dot = interpolate(&pressure_rate);
            let dx = 0.5 * element.length * weight;
            for i in 0..local {
                residual[i] += dx * basis[i] * (area / bulk * p_dot + perimeter * v);
                scale[i] +=
                    dx * basis[i].abs() * ((area / bulk * p_dot).abs() + (perimeter * v).abs());
            }
            energy += 0.5
                * dx
                * perimeter
                * (pin.surface_density * v * v + pin.stiffness_per_area * w * w);
            loss += dx * perimeter * pin.resistance * v * v;
            for (j, (modulus, tau)) in law.terms.iter().enumerate() {
                let stiffness = pin.stiffness_per_area * modulus / law.e_inf;
                let z = interpolate(&motion.relaxation_pressure_pa[j]);
                energy += 0.5 * dx * perimeter * z * z / stiffness;
                loss += dx * perimeter * z * z / (stiffness * tau);
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

#[test]
fn spatial_gas_memory_and_ports_are_retained_and_only_physical_losses_remain() {
    let field = field([300.0, 450.0], [0.9, 0.2], 101325.0);
    let gas = build(&field, 0.001, DgBoundary::Port(PortId(7)))
        .with_side_ports(&[DgSidePort {
            id: PortId(8),
            downstream_element: 1,
        }])
        .unwrap()
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 8,
            quadrature_points: 3,
            band_rad_s: [10.0, 5000.0],
            maximum_relative_error: 0.01,
        })
        .unwrap();
    let old = gas.total_dofs();
    let law = GeneralizedMaxwell::new(1e6, vec![(8e5, 0.003)]).unwrap();
    let dg = gas
        .clone()
        .with_compliant_wall_relaxation(wall(300.0), &law)
        .unwrap();
    let n = dg.total_dofs();
    assert_eq!(n, old + 3 * dg.field_dofs());
    assert_eq!(dg.gas_field(), Some(&field));
    assert_eq!(dg.side_ports(), gas.side_ports());
    assert_eq!(
        dg.side_port_source(PortId(8)).unwrap(),
        gas.side_port_source(PortId(8)).unwrap()
    );
    for i in 0..old {
        assert_eq!(
            &dg.l_tilde[i * n..i * n + old],
            &gas.l_tilde[i * old..(i + 1) * old]
        );
    }
    for ((id, port), (old_id, old_port)) in dg.ports.iter().zip(&gas.ports) {
        assert_eq!(id, old_id);
        assert_eq!(&port[..old], old_port.as_slice());
        assert!(port[old..].iter().all(|x| *x == 0.0));
    }
    let mut state = dg.zero_state();
    for (i, x) in state.tilde.iter_mut().enumerate() {
        *x = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
    }
    let storage_rate: f64 = dg
        .l_tilde
        .chunks_exact(n)
        .zip(&state.tilde)
        .map(|(row, x)| {
            x * row
                .iter()
                .zip(&state.tilde)
                .map(|(a, y)| a * y)
                .sum::<f64>()
        })
        .sum();
    let loss = dg.wall_dissipated_power_w(&state).unwrap();
    close(storage_rate, -loss, loss);
    assert!(loss > dg.compliant_wall_state(&state).unwrap().dissipated_power_w);
    assert!(dg.compliant_wall_state(&gas.zero_state()).is_err());

    // The original side pressure remains conjugate to its injected volume,
    // now with actual material wall storage in the same exact time evolution.
    let lossless = build(&field, RADIUS, DgBoundary::Closed)
        .with_side_ports(&[DgSidePort {
            id: PortId(8),
            downstream_element: 1,
        }])
        .unwrap()
        .with_compliant_wall(wall(0.0))
        .unwrap();
    let propagation = lossless.propagator(2e-6, 1e-13).unwrap();
    let mut current = lossless.zero_state();
    let mut work = 0.0;
    for _ in 0..30 {
        work += propagation
            .step(
                &mut current,
                &[PortInput {
                    id: PortId(8),
                    effort: 1e-7,
                }],
            )
            .unwrap()
            .ledger
            .port_work_step_j;
    }
    close(current.stored_energy_j(), work, work);
    assert!(
        lossless
            .compliant_wall_state(&current)
            .unwrap()
            .stored_energy_j
            > 0.0
    );
}

#[test]
fn uniform_route_is_identical_and_cancelled_spatial_candidates_are_not_published() {
    let uniform = field([300.0; 2], [0.8; 2], 101325.0);
    let bound = build(&uniform, RADIUS, DgBoundary::Closed);
    let direct = DgDuct::new(
        &bound.geometry,
        bound.gas(),
        bound.mesh,
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    let law = GeneralizedMaxwell::new(1e6, vec![(8e5, 0.003)]).unwrap();
    let a = bound
        .with_compliant_wall_relaxation(wall(0.0), &law)
        .unwrap();
    let b = direct
        .with_compliant_wall_relaxation(wall(0.0), &law)
        .unwrap();
    assert_eq!(a.l_tilde, b.l_tilde);
    assert_eq!(a.operator_identity, b.operator_identity);
    assert_eq!(
        a.compliant_wall_state(&a.zero_state()).unwrap(),
        b.compliant_wall_state(&b.zero_state()).unwrap()
    );
    let spatial = field([300.0, 450.0], [0.9, 0.2], 101325.0);
    let original = build(&spatial, RADIUS, DgBoundary::Closed);
    let identity = original.operator_identity;
    for stop in [1, 20, 40] {
        let mut polls = 0;
        assert!(matches!(
            original.clone().with_compliant_wall_relaxation_controlled(
                wall(0.0),
                &law,
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(DuctError::Cancelled { completed_steps: 0 })
        ));
        assert_eq!(original.operator_identity, identity);
        assert!(original.compliant_wall().is_none());
    }
    let elastic = original.clone().with_compliant_wall(wall(0.0)).unwrap();
    let empty = GeneralizedMaxwell::new(1e6, vec![(0.0, 0.003)]).unwrap();
    assert_eq!(
        elastic.operator_identity,
        original
            .clone()
            .with_compliant_wall_relaxation(wall(0.0), &empty)
            .unwrap()
            .operator_identity
    );
    let over_budget = GeneralizedMaxwell::new(1e6, vec![(8e5, 0.003); 400]).unwrap();
    assert!(
        original
            .with_compliant_wall_relaxation(wall(0.0), &over_budget)
            .is_err()
    );
}
