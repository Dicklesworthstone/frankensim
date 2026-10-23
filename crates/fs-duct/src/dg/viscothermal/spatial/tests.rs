use super::*;
use crate::{LossModel, Termination, input_impedance};
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasSpec;
use fs_math::c64::C64;

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

fn spatial(dg: &DgDuct) -> &SpatialViscothermal {
    match dg.viscothermal.as_ref().unwrap() {
        ViscothermalMemory::Spatial(wall) => wall,
        _ => panic!("expected spatial wall"),
    }
}

fn geometry() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.57,
        }],
    }
}

fn resolution() -> DgWallResolution {
    DgWallResolution {
        poles_per_field: 8,
        quadrature_points: 3,
        band_rad_s: [10.0, 5000.0],
        maximum_relative_error: 0.001,
    }
}

fn build(gas: &GasState, mesh: DgMesh, inlet: DgBoundary) -> DgDuct {
    DgDuct::new(&geometry(), gas, mesh, inlet, DgBoundary::Closed)
        .unwrap()
        .with_spatial_viscothermal(DgWallResolution {
            quadrature_points: mesh.order + 1,
            ..resolution()
        })
        .unwrap()
}

// Eliminate retained wall coordinates from THIS assembled time operator.
// The independent reference below is the exact-Bessel cylinder transfer,
// so copying a loss coefficient into an unused diagnostic cannot pass.
fn harmonic_input(dg: &DgDuct, omega: f64) -> C64 {
    let n = dg.total_dofs();
    let m = 2 * dg.field_dofs();
    let s = C64::new(0.0, omega);
    let mut matrix = vec![C64::ZERO; m * m];
    for i in 0..m {
        for j in 0..m {
            matrix[i * m + j] = C64::from_re(-dg.l_tilde[i * n + j]);
            if i == j {
                matrix[i * m + j] = matrix[i * m + j] + s;
            }
            for k in m..n {
                let value = dg.l_tilde[i * n + k] * dg.l_tilde[k * n + j];
                if value != 0.0 {
                    matrix[i * m + j] = matrix[i * m + j]
                        - C64::from_re(value) / (s - C64::from_re(dg.l_tilde[k * n + k]));
                }
            }
        }
    }
    let mut state: Vec<_> = dg.ports[0].1[..m]
        .iter()
        .map(|x| C64::from_re(*x))
        .collect();
    lu_complex(&matrix, m).unwrap().solve(&mut state);
    state
        .iter()
        .zip(&dg.ports[0].1)
        .fold(C64::ZERO, |sum, (x, b)| sum + x.scale(*b))
}

#[test]
fn viscothermal_dg_converges_to_exact_bessel_cylinder_impedance() {
    let gas = gas(101325.0);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let dg = build(
            &gas,
            DgMesh { elements, order: 3 },
            DgBoundary::Port(PortId(7)),
        );
        let mut error = 0.0_f64;
        for omega in [400.0, 1000.0, 1700.0] {
            // Frequency-domain duct owner uses exp(-i omega t).
            let exact = input_impedance(
                &geometry(),
                &gas,
                omega,
                LossModel::Bessel,
                Termination::Closed,
            )
            .unwrap()
            .impedance
            .conj();
            let actual = harmonic_input(&dg, omega);
            error = error.max((actual - exact).abs() / exact.abs());
            assert!(actual.re > 0.0);
        }
        errors.push(error);
    }
    assert!(errors[1] < errors[0], "{errors:?}");
    assert!(errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 3e-4, "{errors:?}");
}

#[test]
fn viscothermal_wall_power_is_the_independent_negative_energy_derivative() {
    let dg = build(
        &gas(101325.0),
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
    );
    let mut state = dg.zero_state();
    for (i, value) in state.tilde.iter_mut().enumerate() {
        *value = det::sin(0.37 * (i + 1) as f64) * 0.01;
    }
    state.initial_energy_j = state.stored_energy_j();
    let n = dg.total_dofs();
    let derivative: f64 = (0..n)
        .map(|i| {
            state.tilde[i]
                * (0..n)
                    .map(|j| dg.l_tilde[i * n + j] * state.tilde[j])
                    .sum::<f64>()
        })
        .sum();
    let power = dg.wall_dissipated_power_w(&state).unwrap();
    assert!(power > 0.0);
    assert!(
        (derivative + power).abs() < 2e-11 * power,
        "{derivative} {power}"
    );
    assert!(state.wall_stored_energy_j() > 0.0);
    let prop = dg.propagator(1e-5, 1e-9).unwrap();
    let before = state.stored_energy_j();
    prop.step(&mut state, &[]).unwrap();
    assert!(state.stored_energy_j() < before);
}

#[test]
fn retained_thermal_memory_recovers_isothermal_closed_volume_and_poiseuille_drag() {
    let gas = gas(101325.0);
    let mesh = DgMesh {
        elements: 1,
        order: 1,
    };
    let dg = DgDuct::new(
        &geometry(),
        &gas,
        mesh,
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_spatial_viscothermal(DgWallResolution {
        poles_per_field: 2,
        quadrature_points: 2,
        band_rad_s: [0.1, 1.0],
        maximum_relative_error: 0.01,
    })
    .unwrap();
    let wall = spatial(&dg);
    for channel in &wall.channels {
        let storage: f64 = channel
            .terms
            .iter()
            .map(|t| t.gain_per_s / t.rate_per_s)
            .sum();
        let expected = if channel.field_start < dg.field_dofs {
            gas.gamma - 1.0
        } else {
            1.0 / 3.0
        };
        assert!((storage - expected).abs() < 2e-12);
    }
    // Start with a uniform 1-Pa adiabatic disturbance and no previous wall
    // heat history. The same sealed gas mass cools to p0/gamma, not zero.
    let mut state = dg.project_initial_state(|_| 1.0, |_| 0.0);
    let initial = state.stored_energy_j();
    let prop = dg.propagator(0.02, 1e-8).unwrap();
    for _ in 0..50 {
        let row = prop.step(&mut state, &[]).unwrap().ledger;
        assert!(!row.balance_checked); // cumulative loss is still a deficit
        assert!(row.cumulative_dissipated_j >= -1e-12 * initial);
    }
    prop.sync_nodal(&mut state).unwrap();
    for pressure in state.pressure() {
        assert!((pressure - 1.0 / gas.gamma).abs() < 1e-7);
    }
    assert!((state.stored_energy_j() / initial - 1.0 / gas.gamma).abs() < 1e-7);
    assert!(state.wall_stored_energy_j() > 0.0);

    // Uniform flow has R = 8 mu L/(pi a^4), independently of gas density.
    // Equilibrate only the viscous arm coordinates algebraically to isolate
    // the exact Poiseuille term; no steady-flow solution is claimed here.
    let flow = 1e-9;
    let mut state = dg.zero_state();
    state.tilde = dg.to_energy_coords(&vec![0.0; dg.field_dofs], &vec![flow; dg.field_dofs]);
    for channel in &wall.channels {
        let value: f64 = channel
            .projection
            .iter()
            .enumerate()
            .map(|(i, w)| w * state.tilde[channel.field_start + i])
            .sum();
        for (arm, term) in channel.terms.iter().enumerate() {
            state.tilde[channel.memory_start + arm] =
                det::sqrt(term.gain_per_s / term.rate_per_s) * value;
        }
    }
    let expected = 8.0 * gas.dynamic_viscosity * 0.57 / (core::f64::consts::PI * 0.001_f64.powi(4))
        * flow
        * flow;
    assert!((dg.wall_dissipated_power_w(&state).unwrap() / expected - 1.0).abs() < 2e-12);
}

#[test]
fn pressure_and_source_gas_field_change_actual_wall_history_and_port_response() {
    let mesh = DgMesh {
        elements: 1,
        order: 2,
    };
    let a = build(&gas(101325.0), mesh, DgBoundary::Port(PortId(1)));
    let b = build(&gas(202650.0), mesh, DgBoundary::Port(PortId(1)));
    assert_eq!(a.gas().sound_speed, b.gas().sound_speed);
    assert_ne!(a.operator_identity(), b.operator_identity());
    for (a, b) in spatial(&a).channels.iter().zip(&spatial(&b).channels) {
        for (a, b) in a.terms.iter().zip(&b.terms) {
            assert_eq!(a.rate_per_s, 2.0 * b.rate_per_s);
            assert_eq!(a.gain_per_s, 2.0 * b.gain_per_s);
        }
    }
    let za = harmonic_input(&a, 1000.0);
    let zb = harmonic_input(&b, 1000.0);
    assert!((zb - za.scale(2.0)).abs() > 1e-3 * za.abs());

    let uniform = crate::profile::gas_field::tests::field([300.0; 2], [0.8; 2], 101325.0);
    let changed = crate::profile::gas_field::tests::field([300.0, 330.0], [0.8, 0.2], 101325.0);
    let constant = build(
        uniform.uniform_state().unwrap().state(),
        mesh,
        DgBoundary::Port(PortId(1)),
    );
    let from_field = |field: &AxialGasField| {
        DgDuct::new_with_gas_field(
            &geometry(),
            field,
            mesh,
            DgBoundary::Port(PortId(1)),
            DgBoundary::Closed,
        )
        .unwrap()
        .with_spatial_viscothermal(resolution())
        .unwrap()
    };
    let uniform_dg = from_field(&uniform);
    assert_eq!(constant.l_tilde, uniform_dg.l_tilde);
    assert_eq!(constant.operator_identity(), uniform_dg.operator_identity());
    let changed_dg = from_field(&changed);
    assert_eq!(changed_dg.gas_field(), Some(&changed));
    assert_ne!(changed_dg.l_tilde, uniform_dg.l_tilde);
    assert!(
        (harmonic_input(&changed_dg, 1000.0) - harmonic_input(&uniform_dg, 1000.0)).abs()
            > 1e-3 * harmonic_input(&uniform_dg, 1000.0).abs()
    );
}

#[test]
fn wall_state_binding_band_and_total_work_limits_refuse_without_fallback() {
    let gas = gas(101325.0);
    let mesh = DgMesh {
        elements: 1,
        order: 2,
    };
    let dg = build(&gas, mesh, DgBoundary::Port(PortId(1)));
    let other = DgDuct::new(
        &geometry(),
        &gas,
        mesh,
        DgBoundary::Port(PortId(1)),
        DgBoundary::Closed,
    )
    .unwrap();
    let prop = dg.propagator(1e-5, 1e-9).unwrap();
    let mut foreign = other.zero_state();
    let before = foreign.clone();
    assert!(prop.step(&mut foreign, &[]).is_err());
    assert_eq!(foreign, before);
    assert!(dg.wall_dissipated_power_w(&foreign).is_err());
    for resolution in [
        DgWallResolution {
            poles_per_field: 0,
            ..resolution()
        },
        DgWallResolution {
            quadrature_points: 2,
            ..resolution()
        },
        DgWallResolution {
            band_rad_s: [0.0, 10.0],
            ..resolution()
        },
        DgWallResolution {
            poles_per_field: 2,
            maximum_relative_error: 1e-8,
            ..resolution()
        },
        DgWallResolution {
            poles_per_field: 64,
            quadrature_points: 64,
            ..resolution()
        },
    ] {
        assert!(
            DgDuct::new(
                &geometry(),
                &gas,
                mesh,
                DgBoundary::Closed,
                DgBoundary::Closed
            )
            .unwrap()
            .with_spatial_viscothermal(resolution)
            .is_err()
        );
    }
    let input = [PortInput {
        id: PortId(1),
        effort: 1e-9,
    }];
    let mut a = dg.zero_state();
    let mut b = dg.zero_state();
    let fine = dg.propagator(5e-6, 1e-9).unwrap();
    for _ in 0..4 {
        prop.step(&mut a, &input).unwrap();
        fine.step(&mut b, &input).unwrap();
        fine.step(&mut b, &input).unwrap();
    }
    let scale = det::sqrt(a.tilde.iter().map(|x| x * x).sum::<f64>());
    let error = det::sqrt(
        a.tilde
            .iter()
            .zip(&b.tilde)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>(),
    );
    assert!(error / scale < 1e-8);
    assert!(a.wall_stored_energy_j() > 0.0);
}

#[test]
fn uniform_exact_mass_and_spatial_wall_routes_have_the_same_pressure_and_energy() {
    let gas = gas(101325.0);
    let base = DgDuct::new(
        &geometry(),
        &gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(1)),
        DgBoundary::Closed,
    )
    .unwrap();
    let r = resolution();
    let uniform = base
        .clone()
        .with_viscothermal(
            r.poles_per_field,
            r.band_rad_s.map(|w| w / core::f64::consts::TAU),
            r.maximum_relative_error,
        )
        .unwrap();
    let spatial = base.with_spatial_viscothermal(r).unwrap();
    assert!(uniform.has_viscothermal_memory() && spatial.has_viscothermal_memory());
    for w in [200.0, 700.0, 1700.0] {
        let a = harmonic_input(&uniform, w);
        assert!((harmonic_input(&spatial, w) - a).abs() / a.abs() < 2e-10);
    }
    let pu = uniform.propagator(1e-5, 1e-9).unwrap();
    let ps = spatial.propagator(1e-5, 1e-9).unwrap();
    let mut a = uniform.zero_state();
    let mut b = spatial.zero_state();
    for _ in 0..12 {
        let input = [PortInput {
            id: PortId(1),
            effort: 1e-9,
        }];
        let response = ps.port_response(&b, PortId(1)).unwrap();
        pu.step(&mut a, &input).unwrap();
        let output = ps.step(&mut b, &input).unwrap();
        let expected = response.history_pressure_pa + response.impedance_pa_s_m3 * 1e-9;
        assert!(
            (output.port_outputs[0].conjugate_rate - expected).abs()
                < 2e-11 * expected.abs().max(1e-30)
        );
    }
    pu.sync_nodal(&mut a).unwrap();
    ps.sync_nodal(&mut b).unwrap();
    let scale = a.pressure().iter().map(|p| p.abs()).fold(0.0_f64, f64::max);
    for (a, b) in a.pressure().iter().zip(b.pressure()) {
        assert!((a - b).abs() < 2e-9 * scale);
    }
    assert!((a.stored_energy_j() / b.stored_energy_j() - 1.0).abs() < 2e-9);
    assert!(
        (uniform.wall_dissipated_power_w(&a).unwrap()
            / spatial.wall_dissipated_power_w(&b).unwrap()
            - 1.0)
            .abs()
            < 2e-9
    );
}

#[test]
fn tapered_source_field_wall_quadrature_is_a_live_independent_axis() {
    let field = crate::profile::gas_field::tests::field([280.0, 360.0], [0.9, 0.1], 101325.0);
    let shape = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.0007,
            outlet_radius: 0.002,
            length: field.length_m(),
        }],
    };
    let base = DgDuct::new_with_gas_field(
        &shape,
        &field,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Port(PortId(3)),
        DgBoundary::Closed,
    )
    .unwrap();
    let response = |points| {
        let dg = base
            .clone()
            .with_spatial_viscothermal(DgWallResolution {
                quadrature_points: points,
                maximum_relative_error: 0.02,
                ..resolution()
            })
            .unwrap();
        assert!(dg.viscothermal_model().is_none());
        assert!(dg.spatial_viscothermal_admission().is_some());
        harmonic_input(&dg, 700.0)
    };
    let reference = response(16);
    let mut previous = f64::INFINITY;
    for points in [2, 4, 8] {
        let error = (response(points) - reference).abs() / reference.abs();
        assert!(
            error < previous,
            "quadrature={points} error={error} previous={previous}"
        );
        previous = error;
    }
    assert!(previous < 1e-5);
    let dg = base
        .with_spatial_viscothermal(DgWallResolution {
            quadrature_points: 4,
            maximum_relative_error: 0.02,
            ..resolution()
        })
        .unwrap();
    let (nodes, _) = gauss_legendre_rule(4);
    let channels = &spatial(&dg).channels;
    for (point, xi) in nodes.iter().enumerate() {
        let x = 0.5 * field.length_m() * (xi + 1.0);
        let local = *field.state_at(x).unwrap().state();
        let radius = 0.0007 + (0.002 - 0.0007) * x / field.length_m();
        let expected = 8.0 * local.dynamic_viscosity / (local.density * radius * radius);
        assert!((channels[2 * point + 1].drag_per_s / expected - 1.0).abs() < 2e-12);
    }
    assert_ne!(channels[1].drag_per_s, channels[7].drag_per_s);
}
