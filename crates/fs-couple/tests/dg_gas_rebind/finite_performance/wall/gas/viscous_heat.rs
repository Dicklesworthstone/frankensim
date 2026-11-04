//! Distinct gas and solid losses on the same original moving-source interval.
use super::*;
use fs_couple::dg_port::finite_gas::{FiniteGasDgHeatBudget, FiniteGasDgLoss};
use fs_couple::dg_port::{DgTerminalChannelInput, DgTerminalChannelKind};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn budgets() -> FiniteGasDgHeatBudget {
    FiniteGasDgHeatBudget {
        molecular: heat_budget(),
        viscous: heat_budget(),
    }
}

fn no_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 0,
        max_dense_bytes: 0,
        max_multiply_terms: 0,
    }
}

fn molecular_controls() -> FiniteGasWallThermalSolve {
    let mut solve = joint_controls();
    solve.gas.loss = FiniteGasDgLoss::UniformMolecularDiffusion {
        maximum_relative_error: 0.02,
    };
    solve
}

fn molecular_source(card: &MaterialCard) -> ApertureDg {
    let cold = field(280.0, 1.0, 101_325.0, LENGTH);
    let relaxation = RelaxationAbsorption::new(
        ["A", "B"]
            .into_iter()
            .map(|name| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: if name == "A" {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 900.0,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 3.0,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 500.0,
                        },
                    }]
                } else {
                    vec![]
                },
            })
            .collect(),
    )
    .unwrap();
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &cold,
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap()
    .with_viscothermal(8, [50.0, 2000.0], 0.02)
    .unwrap()
    .with_molecular_relaxation_field(relaxation)
    .unwrap();
    let original = wall_load(card, true, true);
    let load = DgPortLoad::new_with_regional_material_wall_and_terminal(
        RegionalMaterialWallDg::new(
            vec![
                wall(card, 7.0, THICKNESS, true),
                wall(card, 8.0, 2.0 * THICKNESS, true),
            ],
            duct,
        )
        .unwrap(),
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        original.terminal().unwrap().clone(),
    )
    .unwrap();
    assert!(load.duct().total_dofs() <= 128);
    ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
}

fn assert_caloric_balances(
    result: &fs_couple::bernoulli_aperture::dg::thermal_wall::gas::FiniteGasWallApertureDgStep,
    original_gas: &ClosedGasParcel,
    gas: &ClosedGasParcel,
    input: UniformEnthalpyStepInput<'_>,
    state: EquilibriumPhaseState,
) {
    let du = gas_energy_change(
        original_gas,
        original_gas.state().state().temperature,
        gas.state().state().temperature,
    );
    let dh = input.mass_kg * CP * (state.temperature_k() - input.initial.temperature_k());
    let internal = result.gas.molecular_heat_j
        + (result.gas.viscous_heat_j + result.gas.terminal_viscous_heat_j);
    assert!(
        (du + result.gas.refresh.parameter_work_j
            - result.external_gas_heat_j
            - result.exchanged_heat_j
            - internal)
            .abs()
            < 1e-11,
    );
    assert!(
        (dh + result.wall.parameter_work_j + result.exchanged_heat_j
            - input.internal_heat_j
            - result.wall_heat_j)
            .abs()
            < 1e-11,
    );
    assert_eq!(
        result.gas.supplied_heat_j,
        result.external_gas_heat_j + result.exchanged_heat_j,
    );
    assert_eq!(
        result.exported_dissipation_j,
        (result.coupled.inferred_total_dissipation_j - result.wall_heat_j) - internal,
    );
    let frame = &result.coupled;
    let balance = du
        + dh
        + frame.storage_change_j
        + result.exported_dissipation_j
        + frame
            .terminals
            .iter()
            .map(|t| t.port_work_j - t.channel_work_j)
            .sum::<f64>()
        - frame.upstream_work_j
        - frame.body_work_j
        - frame.aperture.mechanical_work_j
        - result.external_gas_heat_j
        - input.internal_heat_j;
    assert!(balance.abs() < 1e-11);
    assert!(result.combined_energy_residual_j.abs() < 1e-11);
    assert_eq!(gas.mass_kg(), original_gas.mass_kg());
    assert_eq!(gas.volume_m3(), original_gas.volume_m3());
}

#[test]
fn joint_bore_heat_uses_original_solved_flows_and_separate_gas_and_solid_inventories() {
    let (card, caloric) = material(0.0, 800.0);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let mut live = molecular_source(&card);
    for _ in 0..12 {
        live.step(input_drive()).unwrap();
    }
    let before = live.clone();
    let mut frozen = before.clone();
    let native = frozen.step_with_force(input_drive(), 2e-5).unwrap();
    let ports = [
        PortInput {
            id: PortId(0),
            effort: native.aperture.bore_flow_m3_s,
        },
        PortInput {
            id: PortId(1),
            effort: -frozen
                .load()
                .last_terminal_frame()
                .unwrap()
                .outward_flow_m3_s,
        },
    ];
    assert_ne!(ports[0].effort, input_drive().body_flow_m3_s);
    assert!(ports[1].effort.abs() > 0.0);
    let duct = before.load().duct();
    let propagator = before.load().propagator();
    let molecular = duct
        .prepare_molecular_heat(propagator, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(before.duct_state(), &ports, &mut || false)
        .unwrap();
    let viscous = duct
        .prepare_viscous_heat(propagator, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(before.duct_state(), &ports, &mut || false)
        .unwrap();
    let solid = duct
        .prepare_compliant_wall_heat(propagator, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(before.duct_state(), &ports, &mut || false)
        .unwrap();
    assert!(molecular > 0.0 && viscous > 0.0 && solid > 0.0);
    let mut gas = parcel(280.0);
    let original_gas = gas.clone();
    let mut state = initial(caloric.curve());
    let input = body(caloric.curve(), state, true, 0.001);
    let result = live
        .step_with_finite_gas_wall_viscous_molecular_heat_controlled(
            &mut gas,
            &mut state,
            input,
            &interface,
            input_drive(),
            2e-5,
            0.0002,
            molecular_controls(),
            budgets(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.gas.molecular_heat_j, molecular);
    assert_eq!(result.gas.viscous_heat_j, viscous);
    assert_eq!(result.gas.terminal_viscous_heat_j, 0.0);
    assert_eq!(result.wall_heat_j, solid);
    assert!(result.exchanged_heat_j > 0.0);
    assert_eq!(result.coupled.aperture, native.aperture);
    assert_eq!(live.aperture().accepted_steps(), 13);
    assert_eq!(live.duct_state().step_count(), 13);
    assert_eq!(
        live.load().terminal().unwrap().pole_states(),
        frozen.load().terminal().unwrap().pole_states()
    );
    assert_eq!(
        live.load().duct().gas_field().unwrap().uniform_state(),
        Some(gas.state())
    );
    for wall in sources(live.load()) {
        assert_eq!(wall.material().query().axes()["T"], state.temperature_k());
    }
    assert_caloric_balances(&result, &original_gas, &gas, input, state);
}

#[test]
fn paired_joint_heat_returns_terminal_shear_without_depositing_numerical_or_aperture_loss() {
    let (card, caloric) = material(0.0, 800.0);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let load = wall_load(&card, true, true)
        .with_terminal_channels(
            vec![DgTerminalChannelInput {
                port: PortId(1),
                kind: DgTerminalChannelKind::Slit {
                    width_m: 0.002,
                    height_m: 0.0001,
                    length_m: 0.001,
                    transverse_cells: 8,
                    maximum_mach: 0.2,
                },
            }],
            8,
        )
        .unwrap();
    let cold = parcel(300.0);
    let mut gas = ClosedGasParcel::new(
        cold.model().clone(),
        cold.point().clone(),
        cold.volume_m3() + 0.002 * 0.0001 * 0.001,
    )
    .unwrap();
    let original_gas = gas.clone();
    let mut live =
        LoadedApertureDg::new_loaded(paired_source(&load, load.impedance_pa_s_m3()), load).unwrap();
    for _ in 0..12 {
        live.step(input_drive()).unwrap();
    }
    let mut frozen = live.clone();
    let forces = [2e-5, -3e-5];
    let native = frozen.step_with_forces(input_drive(), forces).unwrap();
    let channel = frozen.load().terminal_channels()[0].channel().unwrap();
    let velocity = channel.velocities_m_s();
    let dy = 0.0001 / velocity.len() as f64;
    let gradients = velocity
        .windows(2)
        .map(|v| (v[1] - v[0]).powi(2))
        .sum::<f64>()
        + 2.0 * (velocity[0].powi(2) + velocity[velocity.len() - 1].powi(2));
    let shear = gas.state().state().dynamic_viscosity * 0.001 * 0.002 / dy * gradients * DT;
    let passage = frozen.load().last_terminal_channel_frames()[0]
        .step
        .unwrap();
    assert!(shear > 0.0 && passage.numerical_dissipation_j > 0.0);
    assert!((passage.viscous_dissipation_j - shear).abs() < 2e-12 * shear);
    let mut state = initial(caloric.curve());
    let input = body(caloric.curve(), state, true, 0.001);
    let empty = FiniteGasDgHeatBudget {
        molecular: no_budget(),
        viscous: no_budget(),
    };
    let result = live
        .step_with_finite_gas_wall_viscous_molecular_heat_controlled(
            &mut gas,
            &mut state,
            input,
            &interface,
            input_drive(),
            forces,
            [0.0; 2],
            0.0002,
            joint_controls(),
            empty,
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.gas.molecular_heat_j, 0.0);
    assert_eq!(result.gas.viscous_heat_j, 0.0);
    assert_eq!(
        result.gas.terminal_viscous_heat_j,
        passage.viscous_dissipation_j
    );
    assert!(result.wall_heat_j > 0.0);
    assert!(result.exported_dissipation_j >= passage.numerical_dissipation_j - 1e-11);
    assert_eq!(result.coupled.aperture, native.aperture);
    assert_eq!(pair(&live).blade_states(), pair(&frozen).blade_states());
    assert_ne!(pair(&live).blade_states()[0], pair(&live).blade_states()[1]);
    assert_eq!(pair(&live).gas(), Some(gas.state().state()));
    assert_eq!(live.aperture().accepted_steps(), 13);
    assert_eq!(live.duct_state().step_count(), 13);
    assert_caloric_balances(&result, &original_gas, &gas, input, state);
}

#[test]
fn absent_gas_sources_preserve_existing_joint_transaction_without_integral_budgets() {
    let (card, caloric) = material(0.0, 800.0);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let mut original = joint_source(&card, 280.0);
    for _ in 0..12 {
        original.step(input_drive()).unwrap();
    }
    let mut heated = original.clone();
    let mut gas = parcel(280.0);
    let mut reference_gas = gas.clone();
    let mut state = initial(caloric.curve());
    let mut reference_state = state;
    let input = body(caloric.curve(), state, true, 0.001);
    let empty = FiniteGasDgHeatBudget {
        molecular: no_budget(),
        viscous: no_budget(),
    };
    let actual = heated
        .step_with_finite_gas_wall_viscous_molecular_heat_controlled(
            &mut gas,
            &mut state,
            input,
            &interface,
            input_drive(),
            0.0,
            0.0002,
            joint_controls(),
            empty,
            &mut || false,
        )
        .unwrap();
    let expected = original
        .step_with_finite_gas_wall_heat_controlled(
            &mut reference_gas,
            &mut reference_state,
            input,
            &interface,
            input_drive(),
            0.0,
            0.0002,
            joint_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.coupled, expected.coupled);
    assert_eq!(actual.gas.refresh, expected.gas.refresh);
    assert_eq!(actual.wall, expected.wall);
    assert_eq!(
        actual.exported_dissipation_j,
        expected.exported_dissipation_j
    );
    assert_eq!(actual.exchanged_heat_j, expected.exchanged_heat_j);
    assert_eq!(gas, reference_gas);
    assert_eq!(state, reference_state);
    equal_source(&heated, &original);
}

#[test]
fn gas_integral_refusal_and_final_cancellation_preserve_every_joint_owner() {
    let (card, caloric) = material(0.0, 800.0);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let mut before = molecular_source(&card);
    for _ in 0..12 {
        before.step(input_drive()).unwrap();
    }
    let original_gas = parcel(280.0);
    let original_state = initial(caloric.curve());
    let input = body(caloric.curve(), original_state, true, 0.001);
    let mut clean = before.clone();
    let mut clean_gas = original_gas.clone();
    let mut clean_state = original_state;
    let mut polls = 0;
    let expected = clean
        .step_with_finite_gas_wall_viscous_molecular_heat_controlled(
            &mut clean_gas,
            &mut clean_state,
            input,
            &interface,
            input_drive(),
            0.0,
            0.0002,
            molecular_controls(),
            budgets(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for failure in 0..3 {
        let mut attempted = before.clone();
        let mut gas = original_gas.clone();
        let mut state = original_state;
        let mut budget = budgets();
        match failure {
            0 => budget.molecular = no_budget(),
            1 => budget.viscous = no_budget(),
            _ => {}
        }
        let mut count = 0;
        let result = attempted.step_with_finite_gas_wall_viscous_molecular_heat_controlled(
            &mut gas,
            &mut state,
            input,
            &interface,
            input_drive(),
            0.0,
            0.0002,
            molecular_controls(),
            budget,
            &mut || {
                count += 1;
                failure == 2 && count == polls
            },
        );
        assert!(result.is_err());
        if failure == 2 {
            assert!(matches!(result, Err(AcousticRealizeError::Cancelled)));
            assert_eq!(count, polls);
        }
        assert_eq!(gas, original_gas);
        assert_eq!(state, original_state);
        equal_source(&attempted, &before);
        let retry = attempted
            .step_with_finite_gas_wall_viscous_molecular_heat_controlled(
                &mut gas,
                &mut state,
                input,
                &interface,
                input_drive(),
                0.0,
                0.0002,
                molecular_controls(),
                budgets(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(retry.coupled, expected.coupled);
        assert_eq!(retry.wall, expected.wall);
        assert_eq!(retry.exchanged_heat_j, expected.exchanged_heat_j);
        assert_eq!(gas, clean_gas);
        assert_eq!(state, clean_state);
        equal_source(&attempted, &clean);
    }
}
