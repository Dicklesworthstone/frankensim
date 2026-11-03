//! G1/G3: both original gas-loss channels pay one caloric and native-wave account.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::FiniteGasDgHeatBudget;
use fs_duct::dg::DgWallResolution;
use fs_material::relaxation::RelaxationAbsorption;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_thermochem::SpeciesId;

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn budgets() -> FiniteGasDgHeatBudget {
    FiniteGasDgHeatBudget {
        molecular: budget(),
        viscous: budget(),
    }
}

fn no_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 0,
        max_dense_bytes: 0,
        max_multiply_terms: 0,
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

fn regional_both() -> (
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
    // This manufactured source declares one reservoir on A; catalog species B
    // has no relaxation. The declaration precedes assembly: no native state or
    // spatial quadrature is removed to meet the complete-state integration cap.
    let mut species = molecular_heat::declaration(true).species().to_vec();
    for row in &mut species {
        if row.species == SpeciesId::new("B").unwrap() {
            row.processes.clear();
        }
    }
    let relaxation = RelaxationAbsorption::new(species).unwrap();
    let duct = acoustic(&initial_field(), false)
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 8,
            quadrature_points: 2,
            band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
            maximum_relative_error: 0.02,
        })
        .unwrap()
        .with_molecular_relaxation_field(relaxation)
        .unwrap();
    let bound =
        RegionalMaterialWallDg::new(vec![kelvin_wall(1000.0, 1e8, 2.0, 2e-5), hereditary], duct)
            .unwrap();
    assert_eq!(
        bound
            .duct()
            .molecular_relaxation()
            .unwrap()
            .processes()
            .len(),
        1
    );
    assert!(
        bound.duct().total_dofs() <= 128,
        "complete native state must fit the source owner"
    );
    let prop = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound
        .duct()
        .project_initial_state(|x| 3000.0 + 1000.0 * x, |_| 2e-5);
    for _ in 0..4 {
        prop.step(&mut state, &[]).unwrap();
    }
    (bound, state, prop, parcel(1.0, 300.0, 101325.0))
}

#[test]
fn joint_viscous_molecular_heat_uses_one_native_step_and_closes_the_actual_caloric_balance() {
    let (mut bound, mut state, mut prop, mut gas) = regional_both();
    let (mut reference, mut reference_state, mut reference_prop, mut reference_gas) =
        (bound.clone(), state.clone(), prop.clone(), gas.clone());
    let (initial, old_gas) = (state.clone(), gas.clone());
    let ports = [PortInput {
        id: PortId(7),
        effort: 1e-7,
    }];
    let molecular = bound
        .duct()
        .prepare_molecular_heat(&prop, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&state, &ports, &mut || false)
        .unwrap();
    let viscous = bound
        .duct()
        .prepare_viscous_heat(&prop, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&state, &ports, &mut || false)
        .unwrap();
    assert!(molecular > 0.0 && viscous > 0.0);
    let mut marched = state.clone();
    let native = prop.step(&mut marched, &ports).unwrap();
    let motion = bound.duct().compliant_wall_state(&marched).unwrap();
    assert!(motion.relaxation_stored_energy_j > 0.0);
    let external = 2e-6;
    let internal = molecular + viscous;
    let report = bound
        .step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut state,
            &mut prop,
            &mut gas,
            &ports,
            external,
            energy_solve(),
            solve(),
            budgets(),
            &mut || false,
        )
        .unwrap();
    // Independent composition through the pre-existing caloric owner: supply
    // the measured physical source heat externally and compare its endpoint.
    reference
        .step_with_finite_gas_heat_controlled(
            &mut reference_state,
            &mut reference_prop,
            &mut reference_gas,
            &ports,
            external + internal,
            energy_solve(),
            solve(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(state, reference_state);
    assert_eq!(gas, reference_gas);
    assert_eq!(
        bound.duct().operator_identity(),
        reference.duct().operator_identity()
    );
    assert_eq!(report.molecular_heat_j, molecular);
    assert_eq!(report.viscous_heat_j, viscous);
    assert_eq!(report.supplied_heat_j, external);
    assert_eq!(
        report.exported_acoustic_dissipation_j,
        report.acoustic_dissipation_j - internal
    );
    assert!(report.exported_acoustic_dissipation_j > 0.0); // Thermal and solid loss remain.
    assert_eq!(state.step_count(), initial.step_count() + 1);
    assert_eq!(
        report.acoustic.ledger.port_work_step_j,
        native.ledger.port_work_step_j
    );
    assert_eq!(
        report.acoustic.port_outputs[0].delta_x,
        native.port_outputs[0].delta_x
    );
    same_motion(&bound.duct().compliant_wall_state(&state).unwrap(), &motion);
    assert_eq!(bound.walls(), reference.walls());
    assert_eq!(
        bound.duct().gas_field().unwrap().uniform_state(),
        Some(gas.state())
    );
    assert!(gas.state().state().temperature > old_gas.state().state().temperature);
    let caloric = energy_change(&old_gas, &gas) + report.work_by_gas_j - internal - external;
    assert!(caloric.abs() <= solve().energy_tolerance_j);
    let complete = energy_change(&old_gas, &gas) + state.stored_energy_j()
        - initial.stored_energy_j()
        + report.exported_acoustic_dissipation_j
        - native.ledger.port_work_step_j
        - external;
    assert!(complete.abs() <= solve().energy_tolerance_j);
}

#[test]
fn joint_heat_absent_and_single_sources_preserve_the_existing_native_owner() {
    for source in 0..3 {
        let duct = match source {
            0 => acoustic(&initial_field(), false),
            1 => acoustic(&initial_field(), false)
                .with_viscothermal(8, [50.0, 2000.0], 0.02)
                .unwrap(),
            _ => molecular_heat::molecular_duct(false, true),
        };
        let admission = match source {
            0 => controls(false),
            1 => FiniteGasDgSolve {
                loss: FiniteGasDgLoss::UniformViscothermal {
                    maximum_relative_error: 0.02,
                },
                ..controls(false)
            },
            _ => solve(),
        };
        let mut bound = kelvin_wall(1500.0, 1e9, 2.0, 2e-5).bind_dg(duct).unwrap();
        let mut prop = bound.duct().propagator(1e-6, 1e-13).unwrap();
        let mut state = bound
            .duct()
            .project_initial_state(|x| 2000.0 + 1000.0 * x, |_| 2e-5);
        for _ in 0..2 {
            prop.step(&mut state, &[]).unwrap();
        }
        let mut gas = parcel(1.0, 300.0, 101325.0);
        let (mut original, mut old_state, mut old_prop, mut old_gas) =
            (bound.clone(), state.clone(), prop.clone(), gas.clone());
        // Unused sources need no preparation allowance; active sources keep
        // precisely their original complete-state budget and heat value.
        let allowances = FiniteGasDgHeatBudget {
            molecular: if source == 2 { budget() } else { no_budget() },
            viscous: if source == 1 { budget() } else { no_budget() },
        };
        let combined = bound
            .step_with_finite_gas_viscous_molecular_heat_controlled(
                &mut state,
                &mut prop,
                &mut gas,
                &[],
                1e-5,
                energy_solve(),
                admission,
                allowances,
                &mut || false,
            )
            .unwrap();
        let reference = match source {
            0 => original.step_with_finite_gas_heat_controlled(
                &mut old_state,
                &mut old_prop,
                &mut old_gas,
                &[],
                1e-5,
                energy_solve(),
                admission,
                &mut || false,
            ),
            1 => original.step_with_finite_gas_viscous_heat_controlled(
                &mut old_state,
                &mut old_prop,
                &mut old_gas,
                &[],
                1e-5,
                energy_solve(),
                admission,
                budget(),
                &mut || false,
            ),
            _ => original.step_with_finite_gas_molecular_heat_controlled(
                &mut old_state,
                &mut old_prop,
                &mut old_gas,
                &[],
                1e-5,
                energy_solve(),
                admission,
                budget(),
                &mut || false,
            ),
        }
        .unwrap();
        assert_eq!(state, old_state);
        assert_eq!(gas, old_gas);
        assert_eq!(
            bound.duct().operator_identity(),
            original.duct().operator_identity()
        );
        assert_eq!(combined.acoustic, reference.acoustic);
        assert_eq!(combined.viscous_heat_j, reference.viscous_heat_j);
        assert_eq!(combined.molecular_heat_j, reference.molecular_heat_j);
        assert_eq!(
            combined.exported_acoustic_dissipation_j,
            reference.exported_acoustic_dissipation_j
        );
    }
}

#[test]
fn joint_heat_has_independent_source_budgets_and_atomic_late_cancellation() {
    let (mut bound, mut state, mut prop, mut gas) = regional_both();
    let (old_bound, old_state, old_prop, old_gas) =
        (bound.clone(), state.clone(), prop.clone(), gas.clone());
    for allowances in [
        FiniteGasDgHeatBudget {
            molecular: no_budget(),
            ..budgets()
        },
        FiniteGasDgHeatBudget {
            viscous: no_budget(),
            ..budgets()
        },
    ] {
        assert!(
            bound
                .step_with_finite_gas_viscous_molecular_heat_controlled(
                    &mut state,
                    &mut prop,
                    &mut gas,
                    &[],
                    0.0,
                    energy_solve(),
                    solve(),
                    allowances,
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
    }
    let mut successful_polls = 0;
    old_bound
        .clone()
        .step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut old_state.clone(),
            &mut old_prop.clone(),
            &mut old_gas.clone(),
            &[],
            0.0,
            energy_solve(),
            solve(),
            budgets(),
            &mut || {
                successful_polls += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        bound.step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut state,
            &mut prop,
            &mut gas,
            &[],
            0.0,
            energy_solve(),
            solve(),
            budgets(),
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
    assert_eq!(bound.walls(), old_bound.walls());
    assert_eq!(
        bound.duct().operator_identity(),
        old_bound.duct().operator_identity()
    );
    let mut reference = old_state;
    assert_eq!(
        prop.step(&mut state, &[]).unwrap(),
        old_prop.step(&mut reference, &[]).unwrap()
    );
    assert_eq!(state, reference);
}
