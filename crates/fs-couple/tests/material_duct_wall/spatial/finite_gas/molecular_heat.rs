//! Original molecular heat returns to the same finite gas and continuing wall.
use super::*;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_thermochem::SpeciesId;

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 100_000_000_000,
    }
}
fn solve() -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        loss: FiniteGasDgLoss::UniformMolecularDiffusion {
            maximum_relative_error: 0.02,
        },
        ..controls(false)
    }
}
fn declaration(active: bool) -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["A", "B"]
            .into_iter()
            .enumerate()
            .map(|(i, name)| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: if active {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 900.0 + 600.0 * i as f64,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 3.0 + 7.0 * i as f64,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 500.0,
                        },
                    }]
                } else {
                    Vec::new()
                },
            })
            .collect(),
    )
    .unwrap()
}
fn molecular_duct(diffusion: bool, active: bool) -> DgDuct {
    let duct = if diffusion {
        DgDuct::new_with_gas_field(
            &fs_duct::Duct {
                segments: vec![fs_duct::Segment::Cylinder {
                    radius: RADIUS,
                    length: LENGTH,
                }],
            },
            &initial_field(),
            DgMesh {
                elements: 1,
                order: 1,
            },
            DgBoundary::Port(PortId(7)),
            DgBoundary::Closed,
        )
        .unwrap()
        .with_viscothermal(8, [50.0, 2000.0], 0.02)
        .unwrap()
    } else {
        acoustic(&initial_field(), false)
    };
    duct.with_molecular_relaxation_field(declaration(active))
        .unwrap()
}
fn regional_molecular() -> (
    RegionalMaterialWallDg,
    DgState,
    DgPropagator,
    ClosedGasParcel,
) {
    let source = source_wall(1800.0, 4e8);
    let hereditary = LocalMaterialWall::new(
        source.material().clone(),
        2e-5,
        source.axial_constraint(),
        source.loss().clone(),
        source.band_hz(),
    )
    .unwrap();
    let bound = RegionalMaterialWallDg::new(
        vec![kelvin_wall(1000.0, 1e8, 2.0, 2e-5), hereditary],
        molecular_duct(false, true),
    )
    .unwrap();
    let propagator = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound
        .duct()
        .project_initial_state(|x| 3000.0 + 1000.0 * x, |_| 2e-6);
    for _ in 0..4 {
        propagator.step(&mut state, &[]).unwrap();
    }
    (bound, state, propagator, parcel(1.0, 300.0, 101325.0))
}
fn same_propagation(a: &DgPropagator, b: &DgPropagator, initial: &DgState) {
    assert_eq!(a.dt(), b.dt());
    assert_eq!(a.tolerance(), b.tolerance());
    let (mut left, mut right) = (initial.clone(), initial.clone());
    assert_eq!(
        a.step(&mut left, &[]).unwrap(),
        b.step(&mut right, &[]).unwrap()
    );
    assert_eq!(left, right);
}

#[test]
fn molecular_heat_returns_once_to_original_finite_gas_and_changes_continuing_regional_motion() {
    let (mut bound, mut state, mut propagator, mut gas) = regional_molecular();
    let (old_bound, initial, old_prop, old_gas) = (
        bound.clone(),
        state.clone(),
        propagator.clone(),
        gas.clone(),
    );
    let ports = [PortInput {
        id: PortId(7),
        effort: 1e-7,
    }];
    let channel = bound
        .duct()
        .prepare_molecular_heat(&propagator, budget(), &mut || false)
        .unwrap();
    let expected_heat = channel
        .integrated_heat_j(&initial, &ports, &mut || false)
        .unwrap();
    let mut marched = initial.clone();
    let native = propagator.step(&mut marched, &ports).unwrap();
    let wall_motion = bound.duct().compliant_wall_state(&marched).unwrap();
    assert!(wall_motion.relaxation_stored_energy_j > 0.0);
    let step = bound
        .step_with_finite_gas_molecular_heat_controlled(
            &mut state,
            &mut propagator,
            &mut gas,
            &ports,
            0.0,
            energy_solve(),
            solve(),
            budget(),
            &mut || false,
        )
        .unwrap();
    assert!(step.molecular_heat_j > 0.0);
    assert_eq!(step.molecular_heat_j, expected_heat);
    assert_eq!(step.supplied_heat_j, 0.0);
    assert_eq!(step.transported_enthalpy_j, 0.0);
    assert_eq!(
        step.exported_acoustic_dissipation_j,
        step.acoustic_dissipation_j - step.molecular_heat_j
    );
    assert!(step.exported_acoustic_dissipation_j > 0.0); // Existing wall viscosity.
    assert!(gas.state().state().temperature > old_gas.state().state().temperature);
    assert!(gas.state().state().pressure > old_gas.state().state().pressure);
    assert_eq!(gas.mass_kg(), old_gas.mass_kg());
    assert_eq!(gas.model().source_packs(), old_gas.model().source_packs());
    assert_eq!(
        gas.point().axis_quantities().get("material-lot"),
        old_gas.point().axis_quantities().get("material-lot")
    );
    assert_eq!(bound.walls(), old_bound.walls());
    same_motion(
        &bound.duct().compliant_wall_state(&state).unwrap(),
        &wall_motion,
    );
    assert_eq!(state.step_count(), initial.step_count() + 1);
    assert_eq!(state.initial_energy_j(), initial.initial_energy_j());
    assert_eq!(propagator.dt(), old_prop.dt());
    assert_eq!(
        step.acoustic.ledger.port_work_step_j,
        native.ledger.port_work_step_j
    );
    assert_eq!(
        step.acoustic.port_outputs[0].delta_x,
        native.port_outputs[0].delta_x
    );
    assert_eq!(
        bound.duct().gas_field().unwrap().uniform_state(),
        Some(gas.state())
    );
    assert_ne!(
        bound.duct().operator_identity(),
        old_bound.duct().operator_identity()
    );
    assert_ne!(
        bound.duct().molecular_relaxation().unwrap(),
        old_bound.duct().molecular_relaxation().unwrap()
    );
    let independently_closed = energy_change(&old_gas, &gas) + state.stored_energy_j()
        - initial.stored_energy_j()
        + step.exported_acoustic_dissipation_j
        - native.ledger.port_work_step_j;
    assert!(independently_closed.abs() <= solve().energy_tolerance_j);
    assert!(step.gas_energy_residual_j.abs() <= solve().energy_tolerance_j);
    assert!(step.work_residual_j.abs() <= solve().work_tolerance_j);
    assert!(step.combined_energy_residual_j.abs() <= solve().energy_tolerance_j);

    // The explicitly exported generic path remains available and heats no gas.
    let (mut exported_bound, mut exported_state, mut exported_prop, mut exported_gas) =
        (old_bound, initial, old_prop, old_gas);
    let exported = exported_bound
        .step_with_finite_gas_heat_controlled(
            &mut exported_state,
            &mut exported_prop,
            &mut exported_gas,
            &ports,
            0.0,
            energy_solve(),
            solve(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(exported.molecular_heat_j, 0.0);
    assert_eq!(
        exported.exported_acoustic_dissipation_j,
        exported.acoustic_dissipation_j
    );
    assert!(gas.state().state().temperature > exported_gas.state().state().temperature);
    let mut difference = 0.0;
    for _ in 0..8 {
        propagator.step(&mut state, &[]).unwrap();
        exported_prop.step(&mut exported_state, &[]).unwrap();
        propagator.sync_nodal(&mut state).unwrap();
        exported_prop.sync_nodal(&mut exported_state).unwrap();
        difference += state
            .pressure()
            .iter()
            .zip(exported_state.pressure())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>();
    }
    assert!(difference > 0.0);
}

#[test]
fn molecular_heat_excludes_viscothermal_wall_transport_and_retains_all_radial_memory() {
    let mut bound = kelvin_wall(1500.0, 1e9, 2.0, 2e-5)
        .bind_dg(molecular_duct(true, true))
        .unwrap();
    let mut propagator = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound
        .duct()
        .project_initial_state(|x| 2000.0 + 1000.0 * x, |_| 2e-6);
    for _ in 0..3 {
        propagator.step(&mut state, &[]).unwrap();
    }
    let mut gas = parcel(1.0, 300.0, 101325.0);
    let (old_gas, initial, original) = (gas.clone(), state.clone(), bound.duct().clone());
    let step = bound
        .step_with_finite_gas_molecular_heat_controlled(
            &mut state,
            &mut propagator,
            &mut gas,
            &[],
            1e-5,
            energy_solve(),
            solve(),
            budget(),
            &mut || false,
        )
        .unwrap();
    assert!(step.molecular_heat_j > 0.0);
    assert!(step.exported_acoustic_dissipation_j > 0.0);
    assert!(step.acoustic_dissipation_j > step.molecular_heat_j);
    assert_eq!(step.supplied_heat_j, 1e-5);
    assert_eq!(bound.duct().total_dofs(), original.total_dofs());
    assert_eq!(
        bound
            .duct()
            .viscothermal_model()
            .unwrap()
            .series_terms()
            .len(),
        original.viscothermal_model().unwrap().series_terms().len()
    );
    assert_ne!(
        bound.duct().viscothermal_model().unwrap().series_terms()[0].rate_per_s,
        original.viscothermal_model().unwrap().series_terms()[0].rate_per_s
    );
    let residual = energy_change(&old_gas, &gas) + state.stored_energy_j()
        - initial.stored_energy_j()
        + step.exported_acoustic_dissipation_j
        - 1e-5;
    assert!(residual.abs() <= solve().energy_tolerance_j);
}

#[test]
fn molecular_heat_failure_and_late_cancellation_leave_every_accepted_owner_retryable() {
    let (mut bound, mut state, mut propagator, mut gas) = regional_molecular();
    let (old_bound, old_state, old_prop, old_gas) = (
        bound.clone(),
        state.clone(),
        propagator.clone(),
        gas.clone(),
    );
    for (admission, allowance) in [
        (controls(false), budget()),
        (
            solve(),
            DenseExponentialBudget {
                max_multiply_terms: 0,
                ..budget()
            },
        ),
        (
            FiniteGasDgSolve {
                maximum_temperature_increment_k: 1e-14,
                ..solve()
            },
            budget(),
        ),
    ] {
        assert!(
            bound
                .step_with_finite_gas_molecular_heat_controlled(
                    &mut state,
                    &mut propagator,
                    &mut gas,
                    &[],
                    0.0,
                    energy_solve(),
                    admission,
                    allowance,
                    &mut || false,
                )
                .is_err()
        );
        assert_eq!(state, old_state);
        assert_eq!(gas, old_gas);
        assert_eq!(
            bound.duct().operator_identity(),
            old_bound.duct().operator_identity()
        );
        same_propagation(&propagator, &old_prop, &old_state);
    }
    let (mut expected_bound, mut expected_state, mut expected_prop, mut expected_gas) = (
        old_bound.clone(),
        old_state.clone(),
        old_prop.clone(),
        old_gas.clone(),
    );
    let mut successful_polls = 0;
    let expected = expected_bound
        .step_with_finite_gas_molecular_heat_controlled(
            &mut expected_state,
            &mut expected_prop,
            &mut expected_gas,
            &[],
            0.0,
            energy_solve(),
            solve(),
            budget(),
            &mut || {
                successful_polls += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        bound.step_with_finite_gas_molecular_heat_controlled(
            &mut state,
            &mut propagator,
            &mut gas,
            &[],
            0.0,
            energy_solve(),
            solve(),
            budget(),
            &mut || {
                polls += 1;
                polls == successful_polls
            },
        ),
        Err(FiniteGasDgError::Cancelled)
    ));
    assert_eq!(polls, successful_polls);
    assert_eq!(state, old_state);
    assert_eq!(gas, old_gas);
    assert_eq!(
        bound.duct().operator_identity(),
        old_bound.duct().operator_identity()
    );
    same_propagation(&propagator, &old_prop, &old_state);
    assert_eq!(bound.walls(), old_bound.walls());
    let retry = bound
        .step_with_finite_gas_molecular_heat_controlled(
            &mut state,
            &mut propagator,
            &mut gas,
            &[],
            0.0,
            energy_solve(),
            solve(),
            budget(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(state, expected_state);
    assert_eq!(gas, expected_gas);
    assert_eq!(retry.molecular_heat_j, expected.molecular_heat_j);
    assert_eq!(
        bound.duct().operator_identity(),
        expected_bound.duct().operator_identity()
    );
    same_propagation(&propagator, &expected_prop, &expected_state);
}

#[test]
fn declared_empty_molecular_source_keeps_existing_caloric_path_bitwise() {
    let bound = kelvin_wall(1500.0, 1e9, 2.0, 2e-5)
        .bind_dg(molecular_duct(false, false))
        .unwrap();
    let prop = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let state = bound.duct().project_initial_state(|_| 1000.0, |_| 0.0);
    let gas = parcel(1.0, 300.0, 101325.0);
    let (mut heated_bound, mut heated_prop, mut heated_state, mut heated_gas) =
        (bound.clone(), prop.clone(), state.clone(), gas.clone());
    let (mut export_bound, mut export_prop, mut export_state, mut export_gas) =
        (bound, prop, state, gas);
    let integrated = heated_bound
        .step_with_finite_gas_molecular_heat_controlled(
            &mut heated_state,
            &mut heated_prop,
            &mut heated_gas,
            &[],
            1e-5,
            energy_solve(),
            solve(),
            budget(),
            &mut || false,
        )
        .unwrap();
    let exported = export_bound
        .step_with_finite_gas_heat_controlled(
            &mut export_state,
            &mut export_prop,
            &mut export_gas,
            &[],
            1e-5,
            energy_solve(),
            solve(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(integrated.molecular_heat_j, 0.0);
    assert_eq!(
        integrated.acoustic_dissipation_j,
        exported.acoustic_dissipation_j
    );
    assert_eq!(heated_state, export_state);
    assert_eq!(heated_gas, export_gas);
    assert_eq!(
        heated_bound.duct().operator_identity(),
        export_bound.duct().operator_identity()
    );
    same_propagation(&heated_prop, &export_prop, &export_state);
}
