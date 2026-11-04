//! Physical source heat on one complete held-flow trajectory and finite gas.
use super::*;
use fs_couple::dg_port::finite_gas::FiniteGasDgHeatBudget;
use fs_duct::dg::PortInput;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_thermochem::SpeciesId;

#[path = "heat/moving_aperture.rs"]
mod moving_aperture;

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

fn no_budgets() -> FiniteGasDgHeatBudget {
    FiniteGasDgHeatBudget {
        molecular: no_budget(),
        viscous: no_budget(),
    }
}

fn inventory(cold: &AxialGasField, volume: f64) -> ClosedGasParcel {
    let source = gas(cold);
    ClosedGasParcel::new(source.model().clone(), source.point().clone(), volume).unwrap()
}

fn main_volume() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS * LENGTH
}

fn warm(load: &mut DgPortLoad, gas: &mut ClosedGasParcel, solve: FiniteGasDgSolve) {
    for _ in 0..12 {
        load.step_with_finite_gas_heat_controlled(
            gas,
            8e-7,
            0.0,
            energy_solve(),
            solve,
            &mut || false,
        )
        .unwrap();
    }
}

// Integral of the independently authored NASA Cp/R = 3.5 + .001 T source.
// Formation references cancel before subtracting nearly equal temperatures.
fn caloric_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    let initial = old.state().state();
    let t = next.state().state().temperature;
    old.mass_kg()
        * initial.specific_gas_constant
        * (t - initial.temperature)
        * (2.5 + 0.0005 * (t + initial.temperature))
}

fn same_load(a: &DgPortLoad, b: &DgPortLoad) {
    assert_eq!(a.state(), b.state());
    assert_eq!(a.duct().operator_identity(), b.duct().operator_identity());
    assert_eq!(a.duct().gas_field(), b.duct().gas_field());
    assert_eq!(a.stored_energy_j(), b.stored_energy_j());
    assert_eq!(a.cumulative_dissipated_j(), b.cumulative_dissipated_j());
    assert_eq!(a.last_terminal_frames(), b.last_terminal_frames());
    assert_eq!(
        a.last_terminal_channel_frames(),
        b.last_terminal_channel_frames()
    );
    assert_eq!(a.last_branch_frames(), b.last_branch_frames());
    assert_eq!(
        a.terminal().map(|v| v.pole_states()),
        b.terminal().map(|v| v.pole_states())
    );
    assert_eq!(a.terminal_channels().len(), b.terminal_channels().len());
    for (x, y) in a.terminal_channels().iter().zip(b.terminal_channels()) {
        assert_eq!(x.gas(), y.gas());
        assert_eq!(
            x.channel().map(|v| v.velocities_m_s()),
            y.channel().map(|v| v.velocities_m_s())
        );
    }
    assert_eq!(a.branches().len(), b.branches().len());
    for (x, y) in a.branches().iter().zip(b.branches()) {
        assert_eq!(x.state(), y.state());
        assert_eq!(x.duct().gas_field(), y.duct().gas_field());
    }
}

#[test]
fn complete_passage_heat_matches_resolved_shear_and_independent_caloric_balance() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    for kind in [
        kind(),
        DgTerminalChannelKind::Circular {
            radius_m: 4e-5,
            length_m: 0.001,
            radial_cells: 8,
            maximum_mach: 0.2,
        },
    ] {
        let passage_volume = match kind {
            DgTerminalChannelKind::Slit {
                width_m,
                height_m,
                length_m,
                ..
            } => width_m * height_m * length_m,
            DgTerminalChannelKind::Circular {
                radius_m, length_m, ..
            } => core::f64::consts::PI * radius_m * radius_m * length_m,
            _ => unreachable!(),
        };
        let mut load = load(&cold, false)
            .with_terminal_channels(
                vec![DgTerminalChannelInput {
                    port: PortId(1),
                    kind,
                }],
                8,
            )
            .unwrap();
        let mut gas = inventory(&cold, main_volume() + passage_volume);
        warm(&mut load, &mut gas, controls());
        let (initial, old_gas) = (load.clone(), gas.clone());
        let (mut frozen, mut frozen_gas) = (load.clone(), gas.clone());
        let original = frozen
            .step_with_finite_gas_heat_controlled(
                &mut frozen_gas,
                -1e-5,
                0.0,
                energy_solve(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let channel = frozen.terminal_channels()[0].channel().unwrap();
        let velocity = channel.velocities_m_s();
        let previous = initial.terminal_channels()[0]
            .channel()
            .unwrap()
            .velocities_m_s();
        let material = old_gas.state().state();
        // Independently integrate mu*(du/dn)^2 over the actual planar or annular
        // dual faces, including half-cell no-slip distances at the outer walls.
        let (shear, numerical) = match kind {
            DgTerminalChannelKind::Slit {
                width_m,
                height_m,
                length_m,
                ..
            } => {
                let dy = height_m / velocity.len() as f64;
                let gradients = velocity
                    .windows(2)
                    .map(|v| (v[1] - v[0]).powi(2))
                    .sum::<f64>()
                    + 2.0 * (velocity[0].powi(2) + velocity[velocity.len() - 1].powi(2));
                let shear = material.dynamic_viscosity * length_m * width_m / dy * gradients * DT;
                let numerical = 0.5
                    * material.density
                    * length_m
                    * width_m
                    * dy
                    * velocity
                        .iter()
                        .zip(previous)
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f64>();
                (shear, numerical)
            }
            DgTerminalChannelKind::Circular {
                radius_m, length_m, ..
            } => {
                let n = velocity.len();
                let dr = radius_m / n as f64;
                let gradients = (1..n)
                    .map(|i| i as f64 * (velocity[i] - velocity[i - 1]).powi(2))
                    .sum::<f64>()
                    + 2.0 * n as f64 * velocity[n - 1].powi(2);
                let shear = 2.0
                    * core::f64::consts::PI
                    * material.dynamic_viscosity
                    * length_m
                    * gradients
                    * DT;
                let numerical = 0.5
                    * material.density
                    * length_m
                    * core::f64::consts::PI
                    * dr
                    * dr
                    * velocity
                        .iter()
                        .zip(previous)
                        .enumerate()
                        .map(|(i, (a, b))| (2 * i + 1) as f64 * (a - b).powi(2))
                        .sum::<f64>();
                (shear, numerical)
            }
            _ => unreachable!(),
        };
        assert!(shear > 0.0 && numerical > 2e-15);
        let result = load
            .step_with_finite_gas_viscous_molecular_heat_controlled(
                &mut gas,
                -1e-5,
                0.0,
                energy_solve(),
                controls(),
                no_budgets(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(result.molecular_heat_j, 0.0);
        assert_eq!(result.viscous_heat_j, 0.0);
        assert!((result.terminal_viscous_heat_j - shear).abs() < 2e-12 * shear);
        assert!(
            (result.exported_acoustic_dissipation_j - numerical).abs() < 1e-3 * numerical + 2e-18
        );
        assert_eq!(
            result.acoustic_dissipation_j,
            original.acoustic_dissipation_j
        );
        assert_eq!(result.terminal_work_j, original.terminal_work_j);
        assert_eq!(result.supplied_heat_j, 0.0);
        assert_eq!(load.last_terminal_frames(), frozen.last_terminal_frames());
        assert_eq!(
            load.last_terminal_channel_frames(),
            frozen.last_terminal_channel_frames()
        );
        assert_eq!(load.state().step_count(), initial.state().step_count() + 1);
        assert!(gas.state().state().temperature > old_gas.state().state().temperature);
        let delta_u = caloric_change(&old_gas, &gas);
        assert!((delta_u + result.refresh.parameter_work_j - shear).abs() < 2e-14);
        let compression = compression_energy(frozen.state(), material);
        let expected_work = compression * (bulk(gas.state().state()) / bulk(material) - 1.0);
        assert!((result.refresh.parameter_work_j - expected_work).abs() < 2e-13);
        let balance = delta_u + load.stored_energy_j() - initial.stored_energy_j() + numerical
            - result.acoustic.ledger.port_work_step_j;
        assert!(balance.abs() < controls().energy_tolerance_j);
        assert!(result.combined_energy_residual_j.abs() < controls().energy_tolerance_j);
    }
}

fn unadvanced_molecular_load() -> (DgPortLoad, ClosedGasParcel, FiniteGasDgSolve) {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let source = load(&cold, false);
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
            elements: 1,
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
    assert!(
        duct.total_dofs() <= 128,
        "retain the complete native source state"
    );
    let load = DgPortLoad::new_with_terminal(
        duct,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        source.terminal().unwrap().clone(),
    )
    .unwrap();
    let gas = inventory(&cold, main_volume());
    let solve = FiniteGasDgSolve {
        loss: FiniteGasDgLoss::UniformMolecularDiffusion {
            maximum_relative_error: 0.02,
        },
        ..controls()
    };
    (load, gas, solve)
}

fn molecular_load() -> (DgPortLoad, ClosedGasParcel, FiniteGasDgSolve) {
    let (mut load, mut gas, solve) = unadvanced_molecular_load();
    warm(&mut load, &mut gas, solve);
    (load, gas, solve)
}

#[test]
fn complete_bore_heat_uses_actual_implicit_terminal_flow_and_one_native_step() {
    let (mut load, mut gas, solve) = molecular_load();
    let (initial, old_gas) = (load.clone(), gas.clone());
    let (mut frozen, mut frozen_gas) = (load.clone(), gas.clone());
    let native = frozen
        .step_with_finite_gas_heat_controlled(
            &mut frozen_gas,
            -4e-7,
            0.0,
            energy_solve(),
            solve,
            &mut || false,
        )
        .unwrap();
    let ports = [
        PortInput {
            id: PortId(0),
            effort: -4e-7,
        },
        PortInput {
            id: PortId(1),
            effort: -frozen.last_terminal_frame().unwrap().outward_flow_m3_s,
        },
    ];
    assert!(ports[1].effort.abs() > 0.0);
    let molecular_observer = initial
        .duct()
        .prepare_molecular_heat(initial.propagator(), budget(), &mut || false)
        .unwrap();
    let molecular = molecular_observer
        .integrated_heat_j(initial.state(), &ports, &mut || false)
        .unwrap();
    let viscous = initial
        .duct()
        .prepare_viscous_heat(initial.propagator(), budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(initial.state(), &ports, &mut || false)
        .unwrap();
    assert!(molecular > 0.0 && viscous > 0.0);
    let missing_outlet = molecular_observer
        .integrated_heat_j(initial.state(), &ports[..1], &mut || false)
        .unwrap();
    assert!((missing_outlet - molecular).abs() > 1e-6 * molecular);
    let external = 2e-6;
    let result = load
        .step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut gas,
            -4e-7,
            external,
            energy_solve(),
            solve,
            budgets(),
            &mut || false,
        )
        .unwrap();
    let (mut reference, mut reference_gas) = (initial.clone(), old_gas.clone());
    let expected = reference
        .step_with_finite_gas_heat_controlled(
            &mut reference_gas,
            -4e-7,
            external + (molecular + viscous),
            energy_solve(),
            solve,
            &mut || false,
        )
        .unwrap();
    same_load(&load, &reference);
    assert_eq!(gas, reference_gas);
    assert_eq!(result.acoustic, expected.acoustic);
    assert_eq!(result.molecular_heat_j, molecular);
    assert_eq!(result.viscous_heat_j, viscous);
    assert_eq!(result.terminal_viscous_heat_j, 0.0);
    assert_eq!(result.supplied_heat_j, external);
    assert_eq!(result.acoustic_dissipation_j, native.acoustic_dissipation_j);
    assert_eq!(
        result.exported_acoustic_dissipation_j,
        native.acoustic_dissipation_j - (molecular + viscous)
    );
    assert_eq!(
        result.acoustic.ledger.port_work_step_j,
        native.acoustic.ledger.port_work_step_j
    );
    assert_eq!(load.state().step_count(), initial.state().step_count() + 1);
    assert_eq!(
        load.duct().molecular_relaxation_field(),
        initial.duct().molecular_relaxation_field()
    );
    let delta_u = caloric_change(&old_gas, &gas);
    assert!(
        (delta_u + result.refresh.parameter_work_j - external - molecular - viscous).abs() < 2e-13
    );
    let balance = delta_u + load.stored_energy_j() - initial.stored_energy_j()
        + result.exported_acoustic_dissipation_j
        - result.acoustic.ledger.port_work_step_j
        - external;
    assert!(balance.abs() < solve.energy_tolerance_j);
}

#[test]
fn independent_branch_passage_heat_stays_outside_main_inventory_with_zero_source_budgets() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let branch = field(370.0, 0.0, 101_325.0, 0.015);
    let mut load = super::super::regional::branched(&cold, &branch);
    let mut gas = inventory(&cold, main_volume());
    warm(&mut load, &mut gas, controls());
    let (mut reference, mut reference_gas) = (load.clone(), gas.clone());
    let expected = reference
        .step_with_finite_gas_heat_controlled(
            &mut reference_gas,
            -4e-7,
            0.0,
            energy_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    let branch_heat = reference.last_terminal_channel_frames()[0]
        .step
        .unwrap()
        .viscous_dissipation_j;
    assert!(branch_heat > 0.0);
    let result = load
        .step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut gas,
            -4e-7,
            0.0,
            energy_solve(),
            controls(),
            no_budgets(),
            &mut || false,
        )
        .unwrap();
    same_load(&load, &reference);
    assert_eq!(gas, reference_gas);
    assert_eq!(result.acoustic, expected.acoustic);
    assert_eq!(result.molecular_heat_j, 0.0);
    assert_eq!(result.viscous_heat_j, 0.0);
    assert_eq!(result.terminal_viscous_heat_j, 0.0);
    assert_eq!(
        result.acoustic_dissipation_j,
        expected.acoustic_dissipation_j
    );
    assert_eq!(
        result.exported_acoustic_dissipation_j,
        expected.acoustic_dissipation_j
    );
    assert!(result.exported_acoustic_dissipation_j >= branch_heat - 1e-14);
    assert_eq!(load.branches()[0].duct().gas_field(), Some(&branch));
}

#[test]
fn complete_heat_budget_failure_and_last_poll_cancellation_preserve_every_owner() {
    let (mut attempted, mut gas, solve) = molecular_load();
    let (original, old_gas) = (attempted.clone(), gas.clone());
    let denied = FiniteGasDgHeatBudget {
        molecular: budget(),
        viscous: no_budget(),
    };
    assert!(
        attempted
            .step_with_finite_gas_viscous_molecular_heat_controlled(
                &mut gas,
                -4e-7,
                2e-6,
                energy_solve(),
                solve,
                denied,
                &mut || false,
            )
            .is_err()
    );
    same_load(&attempted, &original);
    assert_eq!(gas, old_gas);
    let (mut reference, mut reference_gas) = (original.clone(), old_gas.clone());
    let mut polls = 0;
    let expected = reference
        .step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut reference_gas,
            -4e-7,
            2e-6,
            energy_solve(),
            solve,
            budgets(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        attempted.step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut gas,
            -4e-7,
            2e-6,
            energy_solve(),
            solve,
            budgets(),
            &mut || {
                count += 1;
                count == polls
            },
        ),
        Err(FiniteGasDgError::Cancelled)
    ));
    assert_eq!(count, polls);
    same_load(&attempted, &original);
    assert_eq!(gas, old_gas);
    let retry = attempted
        .step_with_finite_gas_viscous_molecular_heat_controlled(
            &mut gas,
            -4e-7,
            2e-6,
            energy_solve(),
            solve,
            budgets(),
            &mut || false,
        )
        .unwrap();
    same_load(&attempted, &reference);
    assert_eq!(gas, reference_gas);
    assert_eq!(retry.acoustic, expected.acoustic);
    assert_eq!(retry.molecular_heat_j, expected.molecular_heat_j);
    assert_eq!(retry.viscous_heat_j, expected.viscous_heat_j);
    assert_eq!(
        retry.terminal_viscous_heat_j,
        expected.terminal_viscous_heat_j
    );
}
