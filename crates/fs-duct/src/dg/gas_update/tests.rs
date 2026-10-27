use super::*;
use fs_material::{gas::GasSpec, visco::GeneralizedMaxwell};
use fs_phs::WallPin;

fn gas(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}
fn cylinder(length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.01,
            length,
        }],
    }
}
fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 3e-10 * scale.max(1e-30), "{a:e} != {b:e}");
}
fn bulk(gas: GasState) -> f64 {
    gas.density * gas.sound_speed.powi(2)
}

#[test]
fn live_uniform_gas_continues_the_independent_linear_element_wave_with_one_work_account() {
    let length = 0.5;
    let area = core::f64::consts::PI * 0.01_f64.powi(2);
    let old = gas(300.0, 100_000.0);
    let new = gas(600.0, 150_000.0);
    let mut duct = DgDuct::new(
        &cylinder(length),
        &old,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    let dt = 1e-5;
    let previous_propagation = duct.propagator(dt, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|x| 10.0 * (2.0 * x / length - 1.0), |_| 1e-6);
    for _ in 0..10 {
        previous_propagation.step(&mut state, &[]).unwrap();
    }
    previous_propagation.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let pressure = state.pressure()[1];
    let flow = state.volume_flow()[0];
    let account = duct
        .rebind_uniform_gas(&mut state, new, 1.0, &mut || false)
        .unwrap();
    let p0 = pressure * bulk(new) / bulk(old);
    let q0 = flow * old.density / new.density;
    near(state.pressure()[0], -p0, 10.0);
    near(state.pressure()[1], p0, 10.0);
    for &q in state.volume_flow() {
        near(q, q0, 1e-5);
    }
    // Integral of the original linear compression squared is L*P²/(3*B²);
    // uniform axial momentum gives L*(rho*Q/A)². Geometry is unchanged.
    let work = 0.5
        * area
        * length
        * ((bulk(new) - bulk(old)) * (pressure / bulk(old)).powi(2) / 3.0
            + (1.0 / new.density - 1.0 / old.density) * (old.density * flow / area).powi(2));
    near(account.parameter_work_j, work, before.stored_energy_j());
    near(
        account.energy_after_j - account.energy_before_j,
        work,
        before.stored_energy_j(),
    );
    assert!(work > 0.0);
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    assert_eq!(state.step_count(), 10);
    let checkpoint = state.clone();
    assert!(previous_propagation.step(&mut state, &[]).is_err());
    assert_eq!(state, checkpoint);

    // Independent Galerkin equations for p=P(2x/L-1), q=Q with closed fluxes:
    // P'=6B*Q/(A*L), Q'=-2A*P/(rho*L), omega=sqrt(12)*c/L.
    let a = 6.0 * bulk(new) / (area * length);
    let b = 2.0 * area / (new.density * length);
    let omega = 12.0_f64.sqrt() * new.sound_speed / length;
    let propagation = duct.propagator(dt, 1e-13).unwrap();
    for step in 1..=30 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        let (sin, cos) = (omega * dt * step as f64).sin_cos();
        let p = p0 * cos + a * q0 / omega * sin;
        let q = q0 * cos - b * p0 / omega * sin;
        near(state.pressure()[0], -p, 10.0);
        near(state.pressure()[1], p, 10.0);
        for &actual in state.volume_flow() {
            near(actual, q, 1e-5);
        }
        assert_eq!(
            frame.ledger.cumulative_parameter_work_j,
            account.parameter_work_j
        );
        near(
            frame.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
        near(
            state.stored_energy_j(),
            account.energy_after_j,
            before.initial_energy_j(),
        );
    }
    // Returning to the old temperature/pressure coefficients extracts energy.
    // A second accepted intervention must ADD its signed work to the original
    // account rather than resetting the initial energy or repaying old work.
    let returning_pressure = state.pressure()[1];
    let returning_flow = state.volume_flow()[0];
    let expected_return = 0.5
        * area
        * length
        * ((bulk(old) - bulk(new)) * (returning_pressure / bulk(new)).powi(2) / 3.0
            + (1.0 / old.density - 1.0 / new.density)
                * (new.density * returning_flow / area).powi(2));
    let returning = duct
        .rebind_uniform_gas(&mut state, old, 1.0, &mut || false)
        .unwrap();
    assert!(returning.parameter_work_j < 0.0);
    near(
        returning.parameter_work_j,
        expected_return,
        before.initial_energy_j(),
    );
    assert_eq!(state.step_count(), 40);
    let frame = duct
        .propagator(dt, 1e-13)
        .unwrap()
        .step(&mut state, &[])
        .unwrap();
    assert_eq!(
        frame.ledger.cumulative_parameter_work_j,
        account.parameter_work_j + returning.parameter_work_j
    );
    near(
        frame.ledger.energy_balance_residual_j,
        0.0,
        before.initial_energy_j(),
    );
}

#[test]
fn regional_maxwell_history_ports_and_matched_loads_use_the_same_updated_gas() {
    let old = gas(293.15, 101_325.0);
    let new = gas(390.0, 140_000.0);
    let pins = [
        WallPin {
            surface_density: 0.2,
            stiffness_per_area: 1e7,
            resistance: 30.0,
        },
        WallPin {
            surface_density: 0.7,
            stiffness_per_area: 4e7,
            resistance: 80.0,
        },
    ];
    let laws = [
        GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001)]).unwrap(),
        GeneralizedMaxwell::new(4e6, vec![(5e5, 0.002), (9e5, 0.0004)]).unwrap(),
    ];
    let side = DgSidePort {
        id: PortId(7),
        downstream_element: 1,
    };
    let build = |gas| {
        DgDuct::new(
            &cylinder(0.3),
            &gas,
            DgMesh {
                elements: 2,
                order: 2,
            },
            DgBoundary::Matched,
            DgBoundary::Port(PortId(1)),
        )
        .unwrap()
        .with_side_ports(&[side])
        .unwrap()
        .with_compliant_wall_relaxation_regions(&pins, &[Some(&laws[0]), Some(&laws[1])])
        .unwrap()
    };
    let mut duct = build(old);
    let propagator = duct.propagator(1e-6, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|x| 10.0 + x, |_| 0.0);
    for _ in 0..20 {
        propagator.step(&mut state, &[]).unwrap();
    }
    propagator.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let motion = duct.compliant_wall_state(&state).unwrap();
    assert!(motion.relaxation_stored_energy_j > 0.0);
    let original_ports = duct.ports.clone();
    let account = duct
        .rebind_uniform_gas(&mut state, new, 1.0, &mut || false)
        .unwrap();
    let after = duct.compliant_wall_state(&state).unwrap();
    let m = duct.field_dofs();
    assert_eq!(&state.tilde[2 * m..], &before.tilde[2 * m..]);
    assert_eq!(after.stored_energy_j, motion.stored_energy_j);
    assert_eq!(after.dissipated_power_w, motion.dissipated_power_w);
    for node in 0..m {
        near(
            after.displacement_m[node],
            motion.displacement_m[node],
            1e-6,
        );
        near(after.velocity_m_s[node], motion.velocity_m_s[node], 1e-2);
        for (a, b) in after
            .relaxation_pressure_pa
            .iter()
            .zip(&motion.relaxation_pressure_pa)
        {
            near(a[node], b[node], 10.0);
        }
        near(
            state.pressure()[node] / bulk(new),
            before.pressure()[node] / bulk(old),
            1e-4,
        );
        near(
            new.density * state.volume_flow()[node],
            old.density * before.volume_flow()[node],
            1e-5,
        );
    }
    assert_eq!(duct.side_port_source(side.id).unwrap().gas, new);
    assert_eq!(duct.side_ports(), &[side]);
    for ((id, column), (old_id, previous)) in duct.ports.iter().zip(&original_ports) {
        assert_eq!(id, old_id);
        let effort: f64 = column.iter().zip(&state.tilde).map(|(b, x)| b * x).sum();
        let old_effort: f64 = previous.iter().zip(&before.tilde).map(|(b, x)| b * x).sum();
        near(effort, old_effort * bulk(new) / bulk(old), 10.0);
    }
    // The independent cold native assembly must give the same physical operator
    // and port maps, including the material-dependent matched termination.
    let cold = build(new);
    for (&a, &b) in duct.l_tilde.iter().zip(&cold.l_tilde) {
        near(a, b, b.abs().max(1.0));
    }
    for ((_, a), (_, b)) in duct.ports.iter().zip(&cold.ports) {
        for (&a, &b) in a.iter().zip(b) {
            near(a, b, b.abs().max(1.0));
        }
    }
    let mut derivative = state.clone();
    let n = duct.total_dofs();
    for (i, row) in duct.l_tilde.chunks_exact(n).enumerate() {
        derivative.tilde[i] = row.iter().zip(&state.tilde).map(|(a, x)| a * x).sum();
    }
    let rate = duct.compliant_wall_state(&derivative).unwrap();
    for node in 0..m {
        let pin = pins[node / (duct.mesh.order + 1)];
        let reaction: f64 = after
            .relaxation_pressure_pa
            .iter()
            .map(|arm| arm[node])
            .sum();
        let force = state.pressure()[node]
            - pin.stiffness_per_area * after.displacement_m[node]
            - pin.resistance * after.velocity_m_s[node]
            - reaction;
        near(pin.surface_density * rate.velocity_m_s[node], force, 10.0);
    }
    let propagation = duct.propagator(1e-6, 1e-13).unwrap();
    for _ in 0..10 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        assert_eq!(
            frame.ledger.cumulative_parameter_work_j,
            account.parameter_work_j
        );
        near(
            frame.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
    }
}

#[test]
fn unresolved_gas_history_budget_stale_state_and_late_cancellation_preserve_the_original() {
    let old = gas(300.0, 100_000.0);
    let new = gas(600.0, 150_000.0);
    let original = DgDuct::new(
        &cylinder(0.57),
        &old,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    let initial = original.project_initial_state(|_| 10.0, |_| 1e-6);
    let mut duct = original.clone();
    let mut state = initial.clone();
    assert!(
        duct.rebind_uniform_gas(&mut state, new, 0.0, &mut || false)
            .is_err()
    );
    assert_eq!(state, initial);
    assert_eq!(duct.operator_identity(), original.operator_identity());
    let mut polls = 0;
    original
        .clone()
        .rebind_uniform_gas(&mut initial.clone(), new, 1.0, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        duct.rebind_uniform_gas(&mut state, new, 1.0, &mut || {
            count += 1;
            count == polls
        }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(state, initial);
    assert_eq!(duct.operator_identity(), original.operator_identity());
    let mut foreign = DgDuct::new(
        &cylinder(0.57),
        &new,
        original.mesh,
        original.inlet,
        original.outlet,
    )
    .unwrap()
    .zero_state();
    let saved = foreign.clone();
    assert!(
        duct.rebind_uniform_gas(&mut foreign, new, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(foreign, saved);

    let spatial = crate::profile::gas_field::tests::field([300.0, 450.0], [0.9, 0.1], 100_000.0);
    let uniform = crate::profile::gas_field::tests::field([400.0; 2], [0.5; 2], 130_000.0);
    for mut unsupported in [
        original
            .clone()
            .with_viscothermal(8, [50.0, 2000.0], 0.01)
            .unwrap(),
        DgDuct::new_with_gas_field(
            &cylinder(0.57),
            &spatial,
            original.mesh,
            original.inlet,
            original.outlet,
        )
        .unwrap(),
    ] {
        let mut current = unsupported.project_initial_state(|x| 10.0 + x, |_| 1e-6);
        current.tilde[2 * unsupported.field_dofs()..].fill(1e-6);
        let saved = current.clone();
        let identity = unsupported.operator_identity();
        assert!(
            unsupported
                .rebind_uniform_gas_field(&mut current, &uniform, 1.0, &mut || false)
                .is_err()
        );
        assert_eq!(current, saved);
        assert_eq!(unsupported.operator_identity(), identity);
    }
    let mut sourced = DgDuct::new_with_gas_field(
        &cylinder(0.57),
        &uniform,
        original.mesh,
        original.inlet,
        original.outlet,
    )
    .unwrap();
    let mut current = sourced.zero_state();
    let saved = current.clone();
    assert!(
        sourced
            .rebind_uniform_gas(&mut current, new, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(current, saved);
    assert_eq!(sourced.gas_field(), Some(&uniform));
}
