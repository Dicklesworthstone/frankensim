//! Fluid viscous heat, actual source thermodynamics and continuing wall memory.
use super::*;
use fs_duct::dg::DgWallResolution;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_thermochem::SpeciesId;

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn solve() -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        loss: FiniteGasDgLoss::SpatialViscothermal {
            maximum_relative_error: 0.02,
        },
        ..controls(false)
    }
}

fn regional_viscous() -> (
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
    let duct = acoustic(&initial_field(), false)
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 8,
            quadrature_points: 2,
            band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
            maximum_relative_error: 0.02,
        })
        .unwrap();
    let bound =
        RegionalMaterialWallDg::new(vec![kelvin_wall(1000.0, 1e8, 2.0, 2e-5), hereditary], duct)
            .unwrap();
    let propagator = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound
        .duct()
        .project_initial_state(|x| 3000.0 + 1000.0 * x, |_| 2e-5);
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
fn spatial_viscous_heat_returns_once_to_source_gas_and_retains_regional_histories() {
    let (mut bound, mut state, mut propagator, mut gas) = regional_viscous();
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
    let expected_heat = bound
        .duct()
        .prepare_viscous_heat(&propagator, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&initial, &ports, &mut || false)
        .unwrap();
    assert_eq!(state, initial);
    let mut marched = initial.clone();
    let native = propagator.step(&mut marched, &ports).unwrap();
    let motion = bound.duct().compliant_wall_state(&marched).unwrap();
    assert!(motion.relaxation_stored_energy_j > 0.0);
    let step = bound
        .step_with_finite_gas_viscous_heat_controlled(
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
    assert!(step.viscous_heat_j > 0.0);
    assert_eq!(step.viscous_heat_j, expected_heat);
    assert_eq!(step.molecular_heat_j, 0.0);
    assert_eq!(step.supplied_heat_j, 0.0);
    assert_eq!(step.transported_enthalpy_j, 0.0);
    assert_eq!(
        step.exported_acoustic_dissipation_j,
        step.acoustic_dissipation_j - step.viscous_heat_j
    );
    assert!(step.exported_acoustic_dissipation_j > 0.0);
    assert!(gas.state().state().temperature > old_gas.state().state().temperature);
    // These original synthetic transport cards explicitly declare constant mu.
    // Re-evaluation must honor that source instead of inventing a thermal fit.
    assert_eq!(
        gas.state().state().dynamic_viscosity,
        old_gas.state().state().dynamic_viscosity
    );
    assert_eq!(gas.mass_kg(), old_gas.mass_kg());
    assert_eq!(gas.model().source_packs(), old_gas.model().source_packs());
    assert_eq!(
        gas.point().axis_quantities().get("material-lot"),
        old_gas.point().axis_quantities().get("material-lot")
    );
    assert_eq!(bound.walls(), old_bound.walls());
    same_motion(&bound.duct().compliant_wall_state(&state).unwrap(), &motion);
    assert_eq!(bound.duct().total_dofs(), old_bound.duct().total_dofs());
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
    assert_eq!(
        bound.duct().side_port_source(PortId(7)).unwrap().gas,
        *gas.state().state()
    );
    let residual = energy_change(&old_gas, &gas) + state.stored_energy_j()
        - initial.stored_energy_j()
        + step.exported_acoustic_dissipation_j
        - native.ledger.port_work_step_j;
    assert!(residual.abs() <= solve().energy_tolerance_j);
    assert!(step.gas_energy_residual_j.abs() <= solve().energy_tolerance_j);
    assert!(step.work_residual_j.abs() <= solve().work_tolerance_j);
    assert!(step.combined_energy_residual_j.abs() <= solve().energy_tolerance_j);

    // The existing exported path receives no internal heat. The new source
    // deposit must alter the actual thermodynamic state and later sound.
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
    assert_eq!(exported.viscous_heat_j, 0.0);
    assert!(gas.state().state().temperature > exported_gas.state().state().temperature);
    assert_eq!(
        gas.state().state().dynamic_viscosity,
        exported_gas.state().state().dynamic_viscosity
    );
    let mut difference = 0.0;
    for _ in 0..4 {
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
fn uniform_viscous_heat_keeps_molecular_thermal_and_solid_loss_separate() {
    let declaration = RelaxationAbsorption::new(
        ["A", "B"]
            .into_iter()
            .enumerate()
            .map(|(i, name)| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: vec![RelaxationProcess {
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
                }],
            })
            .collect(),
    )
    .unwrap();
    let duct = DgDuct::new_with_gas_field(
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
    .with_molecular_relaxation_field(declaration)
    .unwrap();
    let mut bound = kelvin_wall(1500.0, 1e9, 2.0, 2e-5).bind_dg(duct).unwrap();
    let mut propagator = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound
        .duct()
        .project_initial_state(|x| 2000.0 + 1000.0 * x, |_| 2e-5);
    for _ in 0..3 {
        propagator.step(&mut state, &[]).unwrap();
    }
    let mut gas = parcel(1.0, 300.0, 101325.0);
    let old_gas = gas.clone();
    let ports = [PortInput {
        id: PortId(7),
        effort: 1e-7,
    }];
    let molecular = bound
        .duct()
        .prepare_molecular_heat(&propagator, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&state, &ports, &mut || false)
        .unwrap();
    let viscous = bound
        .duct()
        .prepare_viscous_heat(&propagator, budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&state, &ports, &mut || false)
        .unwrap();
    assert!(molecular > 0.0);
    assert!(viscous > 0.0);
    let admission = FiniteGasDgSolve {
        loss: FiniteGasDgLoss::UniformMolecularDiffusion {
            maximum_relative_error: 0.02,
        },
        ..controls(false)
    };
    let step = bound
        .step_with_finite_gas_viscous_heat_controlled(
            &mut state,
            &mut propagator,
            &mut gas,
            &ports,
            0.0,
            energy_solve(),
            admission,
            budget(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(step.viscous_heat_j, viscous);
    assert_eq!(step.molecular_heat_j, 0.0);
    assert!(step.exported_acoustic_dissipation_j >= molecular - admission.energy_tolerance_j);
    assert_eq!(
        step.exported_acoustic_dissipation_j,
        step.acoustic_dissipation_j - viscous
    );
    let gas_residual = energy_change(&old_gas, &gas) + step.work_by_gas_j - viscous;
    assert!(gas_residual.abs() <= admission.energy_tolerance_j);
    assert!(step.combined_energy_residual_j.abs() <= admission.energy_tolerance_j);
}

#[test]
fn absent_viscous_source_keeps_the_existing_caloric_path_bitwise() {
    let bound = kelvin_wall(1500.0, 1e9, 2.0, 2e-5)
        .bind_dg(acoustic(&initial_field(), false))
        .unwrap();
    let prop = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let state = bound.duct().project_initial_state(|_| 1000.0, |_| 0.0);
    let gas = parcel(1.0, 300.0, 101325.0);
    let (mut heated_bound, mut heated_prop, mut heated_state, mut heated_gas) =
        (bound.clone(), prop.clone(), state.clone(), gas.clone());
    let (mut exported_bound, mut exported_prop, mut exported_state, mut exported_gas) =
        (bound, prop, state, gas);
    let integrated = heated_bound
        .step_with_finite_gas_viscous_heat_controlled(
            &mut heated_state,
            &mut heated_prop,
            &mut heated_gas,
            &[],
            1e-5,
            energy_solve(),
            controls(false),
            budget(),
            &mut || false,
        )
        .unwrap();
    let exported = exported_bound
        .step_with_finite_gas_heat_controlled(
            &mut exported_state,
            &mut exported_prop,
            &mut exported_gas,
            &[],
            1e-5,
            energy_solve(),
            controls(false),
            &mut || false,
        )
        .unwrap();
    assert_eq!(integrated.viscous_heat_j, 0.0);
    assert_eq!(integrated.molecular_heat_j, 0.0);
    assert_eq!(
        integrated.acoustic_dissipation_j,
        exported.acoustic_dissipation_j
    );
    assert_eq!(heated_state, exported_state);
    assert_eq!(heated_gas, exported_gas);
    assert_eq!(
        heated_bound.duct().operator_identity(),
        exported_bound.duct().operator_identity()
    );
    same_propagation(&heated_prop, &exported_prop, &exported_state);
}

#[test]
fn viscous_heat_refusal_and_late_cancellation_preserve_all_accepted_owners() {
    let (mut bound, mut state, mut prop, mut gas) = regional_viscous();
    let (old_bound, old_state, old_prop, old_gas) =
        (bound.clone(), state.clone(), prop.clone(), gas.clone());
    assert!(
        bound
            .step_with_finite_gas_viscous_heat_controlled(
                &mut state,
                &mut prop,
                &mut gas,
                &[],
                0.0,
                energy_solve(),
                solve(),
                DenseExponentialBudget {
                    max_dimension: 1,
                    ..budget()
                },
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
    same_propagation(&prop, &old_prop, &old_state);
    let mut successful_polls = 0;
    old_bound
        .clone()
        .step_with_finite_gas_viscous_heat_controlled(
            &mut old_state.clone(),
            &mut old_prop.clone(),
            &mut old_gas.clone(),
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
        bound.step_with_finite_gas_viscous_heat_controlled(
            &mut state,
            &mut prop,
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
    assert_eq!(bound.walls(), old_bound.walls());
    assert_eq!(
        bound.duct().operator_identity(),
        old_bound.duct().operator_identity()
    );
    assert_eq!(bound.duct().gas_field(), old_bound.duct().gas_field());
    same_propagation(&prop, &old_prop, &old_state);
}
