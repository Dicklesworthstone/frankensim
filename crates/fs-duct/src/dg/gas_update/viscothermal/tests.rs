use super::*;
use fs_material::{gas::GasSpec, visco::GeneralizedMaxwell};
use fs_phs::WallPin;

fn gas(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}
fn geometry(radius: f64, length: f64) -> Duct {
    Duct { segments: vec![Segment::Cylinder { radius, length }] }
}
fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 5e-9 * scale.max(1e-30), "{a:e} != {b:e}");
}

// Closed-form TWO-state uniform-flow diffusion, independent of the DG matrix
// exponential. The one-section Bessel fraction gives drag=8/tau, g=16/tau,
// p=48/tau. Distinct negative real roots avoid a numerical matrix integrator.
fn viscous_pair(x: f64, y: f64, tau: f64, time: f64) -> (f64, f64) {
    let a = -24.0 / tau;
    let d = -48.0 / tau;
    let b = (16.0_f64 * 48.0).sqrt() / tau;
    let center = 0.5 * (a + d);
    let delta = (0.25 * (a - d).powi(2) + b * b).sqrt();
    let plus = ((center + delta) * time).exp();
    let minus = ((center - delta) * time).exp();
    let cosh = 0.5 * (plus + minus);
    let sinh = 0.5 * (plus - minus) / delta;
    (cosh * x + sinh * ((a - center) * x + b * y),
     cosh * y + sinh * (b * x + (d - center) * y))
}

#[test]
fn continuing_viscous_section_matches_independent_diffusion_and_transport_only_work_is_zero() {
    let old = gas(300.0, 100_000.0);
    let new = gas(600.0, 150_000.0);
    let radius = 1e-4;
    let mut duct = DgDuct::new(
        &geometry(radius, 0.05), &old, DgMesh { elements: 1, order: 1 },
        DgBoundary::IdealOpen, DgBoundary::IdealOpen,
    ).unwrap().with_viscothermal(1, [20.0, 60.0], 0.01).unwrap();
    let dt = 1e-6;
    let previous = duct.propagator(dt, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|_| 0.0, |_| 1e-6);
    for _ in 0..10 { previous.step(&mut state, &[]).unwrap(); }
    previous.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let m = duct.field_dofs();
    assert!(energy(&before.tilde[2 * m..3 * m]) > 0.0);
    let expected_work = (old.density / new.density - 1.0)
        * energy(&before.tilde[m..3 * m]);
    let report = duct.rebind_uniform_viscothermal_gas(
        &mut state, new, 0.01, 1.0, &mut || false,
    ).unwrap();
    near(report.parameter_work_j, expected_work, before.stored_energy_j());
    near(report.energy_after_j - report.energy_before_j, expected_work, before.stored_energy_j());
    assert_eq!(state.step_count(), 10);
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(state.cumulative_dissipated_j(), before.cumulative_dissipated_j());
    let checkpoint = state.clone();
    assert!(previous.step(&mut state, &[]).is_err());
    assert_eq!(state, checkpoint);
    let propagation = duct.propagator(dt, 1e-13).unwrap();
    let tau = new.density * radius * radius / new.dynamic_viscosity;
    for step in 1..=30 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        for i in 0..m {
            let (x, y) = viscous_pair(
                checkpoint.tilde[m + i], checkpoint.tilde[2 * m + i], tau, dt * step as f64,
            );
            near(state.tilde[m + i], x, checkpoint.tilde[m + i].abs());
            near(state.tilde[2 * m + i], y, checkpoint.tilde[m + i].abs());
        }
        assert_eq!(frame.ledger.cumulative_parameter_work_j, report.parameter_work_j);
        near(frame.ledger.energy_balance_residual_j, 0.0, before.initial_energy_j());
    }
    propagation.sync_nodal(&mut state).unwrap();
    let unchanged = state.clone();
    let mut control = state.clone();
    let old_model = duct.viscothermal_model().unwrap().clone();
    let mut transport = new;
    transport.dynamic_viscosity *= 2.0;
    transport.thermal_conductivity *= 3.0;
    transport.prandtl = transport.dynamic_viscosity * transport.specific_heat_cp
        / transport.thermal_conductivity;
    let transport_work = duct.rebind_uniform_viscothermal_gas(
        &mut state, transport, 0.01, 0.0, &mut || false,
    ).unwrap();
    assert_eq!(transport_work.parameter_work_j, 0.0);
    assert_eq!(state.tilde, unchanged.tilde);
    let new_model = duct.viscothermal_model().unwrap();
    near(new_model.series_terms()[0].rate_per_s, 2.0 * old_model.series_terms()[0].rate_per_s,
         old_model.series_terms()[0].rate_per_s);
    near(new_model.shunt_terms()[0].rate_per_s, 3.0 * old_model.shunt_terms()[0].rate_per_s,
         old_model.shunt_terms()[0].rate_per_s);
    let refreshed = duct.propagator(dt, 1e-13).unwrap();
    for _ in 0..20 {
        refreshed.step(&mut state, &[]).unwrap();
        propagation.step(&mut control, &[]).unwrap();
    }
    assert!((state.tilde[m] - control.tilde[m]).abs() > 1e-3 * control.tilde[m].abs());
}

#[test]
fn species_change_preserves_native_arm_momentum_compression_and_regional_wall_histories() {
    let old_field = crate::profile::gas_field::tests::field([300.0; 2], [1.0; 2], 100_000.0);
    let new_field = crate::profile::gas_field::tests::field([420.0; 2], [0.2; 2], 150_000.0);
    let old = *old_field.uniform_state().unwrap().state();
    let new = *new_field.uniform_state().unwrap().state();
    assert_ne!(old.gamma, new.gamma);
    let shape = geometry(0.001, old_field.length_m());
    let pins = [
        WallPin { surface_density: 0.2, stiffness_per_area: 1e7, resistance: 30.0 },
        WallPin { surface_density: 0.7, stiffness_per_area: 4e7, resistance: 80.0 },
    ];
    let laws = [
        GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001)]).unwrap(),
        GeneralizedMaxwell::new(4e6, vec![(5e5, 0.002), (9e5, 0.0004)]).unwrap(),
    ];
    let side = DgSidePort { id: PortId(7), downstream_element: 1 };
    let build = |field: &AxialGasField| {
        DgDuct::new_with_gas_field(
            &shape, field, DgMesh { elements: 2, order: 1 },
            DgBoundary::Matched, DgBoundary::Port(PortId(1)),
        ).unwrap().with_side_ports(&[side]).unwrap()
            .with_viscothermal(16, [50.0, 2000.0], 0.001).unwrap()
            .with_compliant_wall_relaxation_regions(&pins, &[Some(&laws[0]), Some(&laws[1])])
            .unwrap()
    };
    let mut duct = build(&old_field);
    let propagation = duct.propagator(1e-6, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|x| 10.0 + x, |_| 1e-6);
    for _ in 0..20 { propagation.step(&mut state, &[]).unwrap(); }
    propagation.sync_nodal(&mut state).unwrap();
    let original = duct.clone();
    let before = state.clone();
    let old_wall = duct.compliant_wall_state(&state).unwrap();
    assert!(old_wall.relaxation_stored_energy_j > 0.0);
    let m = duct.field_dofs();
    let model = duct.viscothermal_model().unwrap().clone();
    let ns = model.series_terms().len();
    let nt = model.shunt_terms().len();
    let prefix = m * (2 + ns + nt);
    assert!(energy(&before.tilde[2 * m..(2 + ns) * m]) > 0.0);
    assert!(energy(&before.tilde[(2 + ns) * m..prefix]) > 0.0);
    let old_bulk = old.density * old.sound_speed.powi(2);
    let new_bulk = new.density * new.sound_speed.powi(2);
    let expected_work = (new_bulk / old_bulk - 1.0) * energy(&before.tilde[..m])
        + (old.density / new.density - 1.0) * energy(&before.tilde[m..(2 + ns) * m])
        + (new_bulk / old_bulk * (old.gamma - 1.0) / (new.gamma - 1.0) - 1.0)
            * energy(&before.tilde[(2 + ns) * m..prefix]);
    let report = duct.rebind_uniform_viscothermal_gas_field(
        &mut state, &new_field, 0.001, 1.0, &mut || false,
    ).unwrap();
    near(report.parameter_work_j, expected_work, before.stored_energy_j());
    assert_eq!(&state.tilde[prefix..], &before.tilde[prefix..]);
    assert_eq!(duct.gas_field(), Some(&new_field));
    assert_eq!(duct.side_port_source(side.id).unwrap().gas, new);
    assert_eq!(duct.viscothermal_band_rad_s(), original.viscothermal_band_rad_s());
    let next_model = duct.viscothermal_model().unwrap();
    // Independently reconstruct the WEAK physical arm momentum/compression:
    // sqrt(g/p)*M^(1/2)*y = M_arm * physical_arm_flow_or_pressure.
    // This checks every original local basis function, not only total energy.
    for (thermal, old_terms, new_terms, start) in [
        (false, model.series_terms(), next_model.series_terms(), 2 * m),
        (true, model.shunt_terms(), next_model.shunt_terms(), (2 + ns) * m),
    ] {
        for (arm, (a, b)) in old_terms.iter().zip(new_terms).enumerate() {
            for (element, (ea, eb)) in original.elements.iter().zip(&duct.elements).enumerate() {
                let ma = if thermal { &ea.mp_sqrt } else { &ea.mq_sqrt };
                let mb = if thermal { &eb.mp_sqrt } else { &eb.mq_sqrt };
                for i in 0..2 {
                    let mut pa = 0.0;
                    let mut pb = 0.0;
                    for j in 0..2 {
                        let k = start + arm * m + 2 * element + j;
                        pa += ma[2 * i + j] * before.tilde[k];
                        pb += mb[2 * i + j] * state.tilde[k];
                    }
                    pa *= (a.gain_per_s / a.rate_per_s).sqrt();
                    pb *= (b.gain_per_s / b.rate_per_s).sqrt();
                    near(pa, pb, pa.abs().max(pb.abs()));
                }
            }
        }
    }
    let wall = duct.compliant_wall_state(&state).unwrap();
    near(wall.stored_energy_j, old_wall.stored_energy_j, old_wall.stored_energy_j);
    near(wall.dissipated_power_w, old_wall.dissipated_power_w, old_wall.dissipated_power_w);
    for i in 0..m {
        near(wall.displacement_m[i], old_wall.displacement_m[i], 1e-6);
        near(wall.velocity_m_s[i], old_wall.velocity_m_s[i], 1e-2);
        for (a, b) in wall.relaxation_pressure_pa.iter().zip(&old_wall.relaxation_pressure_pa) {
            near(a[i], b[i], 10.0);
        }
    }
    let cold = build(&new_field);
    for (&a, &b) in duct.l_tilde.iter().zip(&cold.l_tilde) {
        near(a, b, b.abs().max(1.0));
    }
    for ((id, a), (other_id, b)) in duct.ports.iter().zip(&cold.ports) {
        assert_eq!(id, other_id);
        for (&a, &b) in a.iter().zip(b) { near(a, b, b.abs().max(1.0)); }
    }
    let n = duct.total_dofs();
    let actual_loss = -(0..n).map(|i| state.tilde[i] * (0..n)
        .map(|j| duct.l_tilde[i * n + j] * state.tilde[j]).sum::<f64>()).sum::<f64>();
    let mut physical_loss = wall.dissipated_power_w
        + next_model.series_drag_per_s() * state.tilde[m..2 * m].iter().map(|v| v * v).sum::<f64>();
    for (branch, (offset, term)) in next_model.series_terms().iter().map(|t| (m, t))
        .chain(next_model.shunt_terms().iter().map(|t| (0, t))).enumerate()
    {
        for i in 0..m {
            physical_loss += (term.gain_per_s.sqrt() * state.tilde[offset + i]
                - term.rate_per_s.sqrt() * state.tilde[(2 + branch) * m + i]).powi(2);
        }
    }
    let z = new.density * new.sound_speed / duct.elements[0].areas[0];
    physical_loss += 0.5 * (state.pressure()[0].powi(2) / z + z * state.volume_flow()[0].powi(2));
    near(actual_loss, physical_loss, physical_loss);
    let continuation = duct.propagator(1e-6, 1e-13).unwrap();
    for _ in 0..10 {
        let frame = continuation.step(&mut state, &[]).unwrap();
        assert_eq!(frame.ledger.cumulative_parameter_work_j, report.parameter_work_j);
        near(frame.ledger.energy_balance_residual_j, 0.0, before.initial_energy_j());
    }
}

#[test]
fn unresolved_band_support_budget_stale_state_and_late_cancel_retain_nonzero_memory() {
    let old = gas(300.0, 100_000.0);
    let new = gas(600.0, 150_000.0);
    let bare = DgDuct::new(
        &geometry(1e-4, 0.05), &old, DgMesh { elements: 1, order: 1 },
        DgBoundary::Closed, DgBoundary::Closed,
    ).unwrap();
    let original = bare.clone().with_viscothermal(1, [20.0, 60.0], 0.01).unwrap();
    let stepping = original.propagator(1e-6, 1e-13).unwrap();
    let mut initial = original.project_initial_state(|_| 10.0, |_| 1e-6);
    for _ in 0..10 { stepping.step(&mut initial, &[]).unwrap(); }
    let mut no_thermal = new;
    no_thermal.gamma = 1.0;
    for (gas, error, budget) in [(new, 1e-16, 1.0), (new, 0.01, 0.0), (no_thermal, 0.01, 1.0)] {
        let mut duct = original.clone();
        let mut state = initial.clone();
        assert!(duct.rebind_uniform_viscothermal_gas(&mut state, gas, error, budget, &mut || false).is_err());
        assert_eq!(state, initial);
        assert_eq!(duct.operator_identity(), original.operator_identity());
        assert_eq!(duct.viscothermal_model(), original.viscothermal_model());
    }
    let mut polls = 0;
    original.clone().rebind_uniform_viscothermal_gas(
        &mut initial.clone(), new, 0.01, 1.0, &mut || { polls += 1; false },
    ).unwrap();
    let mut count = 0;
    let mut duct = original.clone();
    let mut state = initial.clone();
    assert!(matches!(duct.rebind_uniform_viscothermal_gas(
        &mut state, new, 0.01, 1.0, &mut || { count += 1; count == polls },
    ), Err(DuctError::Cancelled { .. })));
    assert_eq!(state, initial);
    assert_eq!(duct.operator_identity(), original.operator_identity());
    let mut foreign = bare.zero_state();
    let foreign_before = foreign.clone();
    assert!(duct.rebind_uniform_viscothermal_gas(&mut foreign, new, 0.01, 1.0, &mut || false).is_err());
    assert_eq!(foreign, foreign_before);
    let mut control = initial.clone();
    assert_eq!(stepping.step(&mut state, &[]).unwrap(), stepping.step(&mut control, &[]).unwrap());
}
