//! G1/G3/G4: physical gas heat follows the original solved moving source.
use super::*;

fn same_scalar(a: &ApertureDg, b: &ApertureDg) {
    same_load(a.load(), b.load());
    assert_eq!(a.aperture().state(), b.aperture().state());
    assert_eq!(a.aperture().spec(), b.aperture().spec());
    assert_eq!(a.aperture().viscous_slit(), b.aperture().viscous_slit());
    assert_eq!(a.aperture().accepted_steps(), b.aperture().accepted_steps());
    assert_eq!(a.stored_energy_j(), b.stored_energy_j());
}

#[test]
fn moving_scalar_returns_resolved_channel_shear_heat_without_depositing_aperture_losses() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let mut live = super::super::scalar(&cold);
    let mut gas = gas(&cold);
    for sample in 0..24 {
        live.step(drive(sample)).unwrap();
    }
    let (initial, old_gas) = (live.clone(), gas.clone());
    let mut frozen = initial.clone();
    let native = frozen.step_with_force(drive(24), 2e-5).unwrap();
    let velocity = frozen.load().terminal_channels()[0]
        .channel().unwrap().velocities_m_s();
    let dy = 4e-5 / velocity.len() as f64;
    // Independent mu*(du/dy)^2 integral, including both half-cell wall faces.
    let gradients = velocity.windows(2).map(|v| (v[1] - v[0]).powi(2)).sum::<f64>()
        + 2.0 * (velocity[0].powi(2) + velocity[velocity.len() - 1].powi(2));
    let shear = old_gas.state().state().dynamic_viscosity * 0.001 * 0.002 / dy * gradients * DT;
    let passage = frozen.load().last_terminal_channel_frames()[0].step.unwrap();
    assert!(shear > 0.0);
    assert!((passage.viscous_dissipation_j - shear).abs() < 2e-12 * shear);
    assert!(native.aperture.dissipated_energy_j > 0.0);
    let external = 2e-6;
    let result = live.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut gas, drive(24), 2e-5, external, energy_solve(), controls(), no_budgets(),
        &mut || false,
    ).unwrap();
    let (mut reference, mut reference_gas) = (initial.clone(), old_gas.clone());
    let expected = reference.step_with_finite_gas_heat_controlled(
        &mut reference_gas, drive(24), 2e-5, external + passage.viscous_dissipation_j,
        energy_solve(), controls(), &mut || false,
    ).unwrap();
    same_scalar(&live, &reference);
    assert_eq!(gas, reference_gas);
    assert_eq!(result.coupled, expected.coupled);
    assert_eq!(result.coupled.aperture, native.aperture);
    assert_eq!(result.gas.molecular_heat_j, 0.0);
    assert_eq!(result.gas.viscous_heat_j, 0.0);
    assert_eq!(result.gas.terminal_viscous_heat_j, passage.viscous_dissipation_j);
    assert_eq!(result.gas.supplied_heat_j, external);
    assert_eq!(result.gas.exported_acoustic_dissipation_j,
        result.gas.acoustic_dissipation_j - passage.viscous_dissipation_j);
    assert_eq!(live.aperture().accepted_steps(), 25);
    assert_eq!(live.duct_state().step_count(), 25);
    let delta_u = caloric_change(&old_gas, &gas);
    assert!((delta_u + result.gas.refresh.parameter_work_j - external - shear).abs() < 2e-13);
    let frame = &result.coupled;
    let balance = delta_u + frame.storage_change_j + frame.inferred_total_dissipation_j
        - shear + frame.terminals.iter().map(|f| f.load_work_j).sum::<f64>()
        - frame.upstream_work_j - frame.body_work_j - frame.aperture.mechanical_work_j - external;
    assert!(balance.abs() < controls().energy_tolerance_j);
    assert!(result.combined_energy_residual_j.abs() < controls().energy_tolerance_j);
}

#[test]
fn moving_pair_bore_heat_uses_actual_solved_flow_and_retains_both_mechanical_histories() {
    let (load, mut gas, solve) = unadvanced_molecular_load();
    let mut live = LoadedApertureDg::new_loaded(
        paired_source(&load, load.impedance_pa_s_m3()), load,
    ).unwrap();
    for sample in 0..12 {
        live.step(drive(sample)).unwrap();
    }
    let (initial, old_gas) = (live.clone(), gas.clone());
    let forces = [2e-5, -3e-5];
    let mut frozen = initial.clone();
    let native = frozen.step_with_forces(drive(12), forces).unwrap();
    let ports = [
        PortInput { id: PortId(0), effort: native.aperture.bore_flow_m3_s },
        PortInput { id: PortId(1), effort: -frozen.load().last_terminal_frame().unwrap().outward_flow_m3_s },
    ];
    assert_ne!(ports[0].effort, drive(12).body_flow_m3_s);
    assert!(ports[1].effort.abs() > 0.0);
    let molecular = initial.load().duct()
        .prepare_molecular_heat(initial.load().propagator(), budget(), &mut || false).unwrap()
        .integrated_heat_j(initial.duct_state(), &ports, &mut || false).unwrap();
    let viscous = initial.load().duct()
        .prepare_viscous_heat(initial.load().propagator(), budget(), &mut || false).unwrap()
        .integrated_heat_j(initial.duct_state(), &ports, &mut || false).unwrap();
    assert!(molecular > 0.0 && viscous > 0.0);
    let internal = molecular + viscous;
    let external = 2e-6;
    let result = live.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut gas, drive(12), forces, [0.0; 2], external, energy_solve(), solve,
        budgets(), &mut || false,
    ).unwrap();
    let (mut reference, mut reference_gas) = (initial.clone(), old_gas.clone());
    let expected = reference.step_with_finite_gas_heat_controlled(
        &mut reference_gas, drive(12), forces, [0.0; 2], external + internal,
        energy_solve(), solve, &mut || false,
    ).unwrap();
    same_load(live.load(), reference.load());
    assert_eq!(gas, reference_gas);
    assert_eq!(result.coupled, expected.coupled);
    assert_eq!(result.coupled.aperture, native.aperture);
    assert_eq!(pair(&live).blade_states(), pair(&frozen).blade_states());
    assert_ne!(pair(&live).blade_states()[0], pair(&live).blade_states()[1]);
    assert_eq!(pair(&live).gas(), Some(gas.state().state()));
    assert_eq!(live.aperture().accepted_steps(), 13);
    assert_eq!(live.duct_state().step_count(), 13);
    assert_eq!(result.gas.molecular_heat_j, molecular);
    assert_eq!(result.gas.viscous_heat_j, viscous);
    assert_eq!(result.gas.terminal_viscous_heat_j, 0.0);
    assert_eq!(result.gas.supplied_heat_j, external);
    assert_eq!(result.gas.exported_acoustic_dissipation_j,
        result.gas.acoustic_dissipation_j - internal);
    let delta_u = caloric_change(&old_gas, &gas);
    assert!((delta_u + result.gas.refresh.parameter_work_j - external - internal).abs() < 2e-13);
    let frame = &result.coupled;
    let balance = delta_u + frame.storage_change_j + frame.inferred_total_dissipation_j
        - internal + frame.terminals.iter().map(|f| f.load_work_j).sum::<f64>()
        - frame.upstream_work_j - frame.body_work_j - frame.aperture.mechanical_work_j - external;
    assert!(balance.abs() < solve.energy_tolerance_j);
    assert!(result.combined_energy_residual_j.abs() < solve.energy_tolerance_j);
}

#[test]
fn independent_branch_and_aperture_losses_preserve_zero_heat_identity() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let branch = field(370.0, 0.0, 101_325.0, 0.015);
    let load = super::super::super::regional::branched(&cold, &branch);
    let mut live = ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap();
    let mut gas = inventory(&cold, main_volume());
    for sample in 0..24 {
        live.step(drive(sample)).unwrap();
    }
    let (mut reference, mut reference_gas) = (live.clone(), gas.clone());
    let expected = reference.step_with_finite_gas_heat_controlled(
        &mut reference_gas, drive(24), 2e-5, 2e-6, energy_solve(), controls(), &mut || false,
    ).unwrap();
    let branch_heat = reference.load().last_terminal_channel_frames()[0].step.unwrap()
        .viscous_dissipation_j;
    assert!(branch_heat > 0.0 && expected.coupled.aperture.dissipated_energy_j > 0.0);
    let result = live.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut gas, drive(24), 2e-5, 2e-6, energy_solve(), controls(), no_budgets(),
        &mut || false,
    ).unwrap();
    same_scalar(&live, &reference);
    assert_eq!(gas, reference_gas);
    assert_eq!(result.coupled, expected.coupled);
    assert_eq!(result.gas.molecular_heat_j, 0.0);
    assert_eq!(result.gas.viscous_heat_j, 0.0);
    assert_eq!(result.gas.terminal_viscous_heat_j, 0.0);
    assert_eq!(result.gas.exported_acoustic_dissipation_j, result.gas.acoustic_dissipation_j);
    assert_eq!(result.combined_energy_residual_j, expected.combined_energy_residual_j);
    assert_eq!(live.load().branches()[0].duct().gas_field(), Some(&branch));
}

#[test]
fn moving_source_budget_failure_and_final_cancellation_keep_all_accepted_owners() {
    let (load, mut gas, solve) = unadvanced_molecular_load();
    let mut attempted = ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap();
    for sample in 0..12 {
        attempted.step(drive(sample)).unwrap();
    }
    let (original, old_gas) = (attempted.clone(), gas.clone());
    let denied = FiniteGasDgHeatBudget { molecular: budget(), viscous: no_budget() };
    assert!(attempted.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut gas, drive(12), 2e-5, 2e-6, energy_solve(), solve, denied, &mut || false,
    ).is_err());
    same_scalar(&attempted, &original);
    assert_eq!(gas, old_gas);
    let (mut reference, mut reference_gas) = (original.clone(), old_gas.clone());
    let mut polls = 0;
    let expected = reference.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut reference_gas, drive(12), 2e-5, 2e-6, energy_solve(), solve, budgets(),
        &mut || { polls += 1; false },
    ).unwrap();
    let mut count = 0;
    assert!(matches!(attempted.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut gas, drive(12), 2e-5, 2e-6, energy_solve(), solve, budgets(),
        &mut || { count += 1; count == polls },
    ), Err(FiniteGasDgError::Cancelled)));
    assert_eq!(count, polls);
    same_scalar(&attempted, &original);
    assert_eq!(gas, old_gas);
    let retry = attempted.step_with_finite_gas_viscous_molecular_heat_controlled(
        &mut gas, drive(12), 2e-5, 2e-6, energy_solve(), solve, budgets(), &mut || false,
    ).unwrap();
    same_scalar(&attempted, &reference);
    assert_eq!(gas, reference_gas);
    assert_eq!(retry.coupled, expected.coupled);
    assert_eq!(retry.gas.molecular_heat_j, expected.gas.molecular_heat_j);
    assert_eq!(retry.gas.viscous_heat_j, expected.gas.viscous_heat_j);
}
