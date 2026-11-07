use super::*;
use fs_material::osmosis::nonideal::VirialSolutionModel;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const RT: f64 = R * 300.0;
const SATURATION: f64 = 1e-5;
const VOLUME: f64 = 1e-6;
const REFERENCE: f64 = 2.0;

fn law(saturation: f64) -> SurfactantLaw {
    SurfactantLaw::Langmuir {
        clean_tension_n_m: 0.2,
        thermal_energy_j_mol: RT,
        saturation_concentration_mol_m2: saturation,
    }
}
fn bath() -> AqueousSolution {
    AqueousSolution::nonideal(
        Arc::new(
            VirialSolutionModel::new(300.0, 18e-6, REFERENCE, vec![20.0; 2], vec![0.0; 4]).unwrap(),
        ),
        VOLUME,
        vec![4e-6, 1e-6],
    )
    .unwrap()
}
fn physical(step: f64, faces: usize, edges: usize) -> SolubleSurfaceControls {
    SolubleSurfaceControls {
        maximum_step_s: step,
        maximum_substeps: 6000,
        maximum_root_evaluations: 1_000_000,
        maximum_surface_pair_exchanges: 50_000,
        maximum_faces: faces,
        maximum_edges: edges,
        root_tolerance_mol: 1e-20,
        amount_tolerance_mol: 1e-15,
        energy_tolerance_j: 1e-10,
    }
}
fn mechanical() -> SurfaceDynamicsControls {
    let mut solve = ImplicitSolveConfig::default();
    solve.max_newton_iterations = 32;
    solve.newton.absolute_tolerance = 1e-13;
    solve.newton.relative_tolerance = 1e-12;
    SurfaceDynamicsControls {
        solve,
        maximum_krylov_restart: 32,
        impulse_tolerance_n_s: 1e-11,
        kinematic_tolerance_m: 1e-12,
        energy_tolerance_j: 1e-10,
    }
}
fn accuracy(step: f64, tolerance: f64) -> AdaptiveSolubleSurfaceDynamicsControls {
    AdaptiveSolubleSurfaceDynamicsControls {
        transport: AdaptiveSolubleSurfaceControls {
            initial_step_s: step,
            minimum_step_s: 1e-7,
            absolute_surface_amount_tolerance_mol: 1e-5 * tolerance,
            absolute_bulk_amount_tolerance_mol: 1e-5 * tolerance,
            absolute_bulk_activity_tolerance: tolerance,
            relative_tolerance: 0.0,
            maximum_trials: 1000,
        },
        absolute_position_tolerance_m: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
        maximum_mechanical_steps: 3000,
    }
}
fn loose(step: f64) -> AdaptiveSolubleSurfaceDynamicsControls {
    let mut a = accuracy(step, 1.0);
    a.transport.absolute_surface_amount_tolerance_mol = 1.0;
    a.transport.absolute_bulk_amount_tolerance_mol = 1.0;
    a
}
fn inviscid() -> SurfaceViscosityField<'static> {
    SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; allowance {tolerance:e}"
    );
}
fn triangle() -> [Vector; 3] {
    [
        [1.0, 0.0, 0.0],
        [-0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [-0.5, -3.0_f64.sqrt() / 2.0, 0.0],
    ]
}
fn ring() -> CapillarySurfaceDynamics<'static> {
    CapillarySurfaceDynamics::try_new(&[0.8; 3], &[[0, 1, 2]], 3, 1).unwrap()
}
fn ring_state(time: f64, speed: f64, amount: f64) -> SecondOrderState {
    let q = triangle();
    ring()
        .initial_state_surfactant(
            time,
            &q,
            &q.map(|p| p.map(|x| speed * x)),
            SurfactantField::Uniform {
                law: law(SATURATION),
                amounts_mol: &[amount],
            },
            inviscid(),
            &[[0.0; 3]; 3],
        )
        .unwrap()
}
fn ring_transport<'a>(
    laws: &'a [SurfactantLaw],
    rates: &'a [f64],
    step: f64,
) -> SolubleSurfaceTransport<'a> {
    SolubleSurfaceTransport {
        laws,
        kinetics: SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: RT * 2.0_f64.ln(),
            desorption_rates_s_inverse: rates,
        },
        edges: &[],
        controls: physical(step, 1, 0),
    }
}

// Independent continuous radial force/site equations for an equilateral
// material triangle. The reference never calls any native transport or
// mechanical stepper; its RK4 mesh is doubled before comparing production.
fn reference(speed: f64, rate: f64, duration: f64, steps: usize) -> [f64; 3] {
    let area = 3.0 * 3.0_f64.sqrt() / 4.0;
    let total = 5e-6;
    let f = |y: [f64; 3]| {
        let capacity = area * y[0] * y[0] * SATURATION;
        let tension = 0.2 + RT * SATURATION * (-y[2] / capacity).ln_1p();
        [
            y[1],
            -3.0_f64.sqrt() * tension * y[0] / (2.0 * 0.8),
            rate * (0.5 * (total - y[2]) / (VOLUME * REFERENCE) * (capacity - y[2]) - y[2]),
        ]
    };
    let h = duration / (steps as f64);
    let mut y = [1.0, speed, 1e-6];
    for _ in 0..steps {
        let k1 = f(y);
        let k2 = f(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = f(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = f(core::array::from_fn(|i| y[i] + h * k3[i]));
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}
fn physical_accounts(report: &AdaptiveSolubleSurfaceDynamicsReport) {
    assert_eq!(
        report.accepted_mechanical_steps,
        2 * report.accepted_intervals
    );
    assert_eq!(report.evaluated_mechanical_steps, 3 * report.trials);
    assert_eq!(
        report.accepted_chemical_substeps,
        4 * report.accepted_intervals
    );
    assert_eq!(report.evaluated_chemical_substeps, 6 * report.trials);
    assert!(report.maximum_accepted_normalized_error <= 1.0);
    assert!(report.evaluated_root_evaluations >= report.accepted_root_evaluations);
    close(
        report.total_amount_mol[0],
        report.total_amount_mol[1],
        1e-15,
    );
    close(report.amount_defect_mol, 0.0, 1e-15);
    close(report.energy_account_residual_j, 0.0, 1e-10);
    close(report.energy_residual_j, 0.0, 1e-10);
    close(
        report.kinetic_energy_j[1] - report.kinetic_energy_j[0],
        report.kinetic_change_j,
        2e-14,
    );
}

#[test]
fn moving_soluble_accuracy_resolves_continuous_breathing_and_finite_depletion() {
    let model = ring();
    let laws = [law(SATURATION)];
    let duration = 0.3;
    for speed in [-0.3, 0.3] {
        let expected = reference(speed, 3.0, duration, 4096);
        let refined = reference(speed, 3.0, duration, 8192);
        for i in 0..3 {
            close(expected[i], refined[i], if i == 2 { 1e-19 } else { 2e-13 });
        }
        let mut errors = Vec::new();
        let mut work = Vec::new();
        for tolerance in [2e-3, 2e-5] {
            let mut state = ring_state(0.0, speed, 1e-6);
            let mut bulk = bath();
            let mut amounts = [1e-6];
            let report = model
                .advance_soluble_surfactant_adaptive_to(
                    &mut state,
                    &mut bulk,
                    &mut amounts,
                    ring_transport(&laws, &[3.0], duration),
                    inviscid(),
                    &[[0.0; 3]; 3],
                    duration,
                    mechanical(),
                    accuracy(duration, tolerance),
                    || false,
                )
                .unwrap();
            physical_accounts(&report);
            assert_eq!(state.t, duration);
            assert_eq!(state.steps, report.accepted_mechanical_steps);
            assert_eq!(state.history.len(), state.steps);
            for pair in state.history.windows(2) {
                assert_eq!(pair[0].t_start + pair[0].h, pair[1].t_start);
                assert_eq!(pair[0].step + 1, pair[1].step);
            }
            close(bulk.solute_moles()[0] + amounts[0], 5e-6, 1e-18);
            assert_eq!(bulk.solute_moles()[1], 1e-6);
            assert_eq!(bulk.solvent_volume_m3(), VOLUME);
            assert!(report.exchange_free_energy_dissipation_j > 0.0);
            errors.push(
                (state.q[0] - expected[0])
                    .abs()
                    .max((state.v[0] - expected[1]).abs())
                    .max((amounts[0] - expected[2]).abs() / SATURATION),
            );
            work.push(report.evaluated_mechanical_steps);
        }
        assert!(errors[1] < 0.3 * errors[0], "speed={speed}: {errors:?}");
        assert!(work[1] > work[0], "{work:?}");
    }
}

fn patch() -> ([Vector; 4], [[usize; 3]; 2], [f64; 4], [Vector; 4]) {
    let q = [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.2, 1.0, 0.2],
    ];
    let v = q.map(|p| [0.2 * p[0], 0.1 * p[1], 0.05 * p[2]]);
    (q, [[0, 1, 2], [1, 3, 2]], [0.7, 0.8, 0.9, 1.0], v)
}

#[test]
fn accepted_moving_fine_path_preserves_history_and_excludes_probe_losses() {
    let (q, faces, masses, v) = patch();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 2).unwrap();
    let laws = [law(SATURATION), law(2.0 * SATURATION)];
    let edges = [SurfacePartitionDiffusionEdge {
        faces: [0, 1],
        conductance_mol_s: 3e-5,
    }];
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity {
        shear_n_s_m: 0.01,
        dilatational_n_s_m: 0.02,
    });
    let forces = [[0.01, 0.0, 0.0]; 4];
    let transport = SolubleSurfaceTransport {
        laws: &laws,
        kinetics: SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: RT * 0.7,
            desorption_rates_s_inverse: &[4.0, 1.0],
        },
        edges: &edges,
        controls: physical(0.125, 2, 1),
    };
    let mut amounts = [2e-7, 2e-6];
    let coating = core::array::from_fn::<_, 2, _>(|i| SurfactantFace {
        law: laws[i],
        amount_mol: amounts[i],
    });
    let mut state = model
        .initial_state_surfactant(
            0.5,
            &q,
            &v,
            SurfactantField::PerFace(&coating),
            viscosity,
            &forces,
        )
        .unwrap();
    let mut bulk = bath();
    model
        .advance_soluble_surfactant(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            viscosity,
            &forces,
            0.0625,
            mechanical(),
            || false,
        )
        .unwrap();
    let original = state.clone();
    let original_bulk = bulk.clone();
    let original_amounts = amounts;
    let mut expected = state.clone();
    let mut expected_bulk = bulk.clone();
    let mut expected_amounts = amounts;
    let first = model
        .advance_soluble_surfactant(
            &mut expected,
            &mut expected_bulk,
            &mut expected_amounts,
            transport,
            viscosity,
            &forces,
            0.0625,
            mechanical(),
            || false,
        )
        .unwrap();
    let second = model
        .advance_soluble_surfactant(
            &mut expected,
            &mut expected_bulk,
            &mut expected_amounts,
            transport,
            viscosity,
            &forces,
            0.0625,
            mechanical(),
            || false,
        )
        .unwrap();
    let mut coarse = original.clone();
    let mut coarse_bulk = original_bulk.clone();
    let mut coarse_amounts = original_amounts;
    let full = model
        .advance_soluble_surfactant(
            &mut coarse,
            &mut coarse_bulk,
            &mut coarse_amounts,
            transport,
            viscosity,
            &forces,
            0.125,
            mechanical(),
            || false,
        )
        .unwrap();
    let target = state.t + 0.125;
    let report = model
        .advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            viscosity,
            &forces,
            target,
            mechanical(),
            loose(0.125),
            || false,
        )
        .unwrap();
    assert_eq!(state, expected);
    assert_eq!(bulk, expected_bulk);
    assert_eq!(amounts, expected_amounts);
    assert_eq!(
        &state.history[..original.history.len()],
        original.history.as_slice()
    );
    assert_eq!(state.history.len(), original.history.len() + 2);
    assert_eq!(report.accepted_intervals, 1);
    assert_eq!(report.accepted_surface_pair_exchanges, 8);
    assert_eq!(report.evaluated_surface_pair_exchanges, 12);
    assert_eq!(
        report.accepted_root_evaluations,
        first.total_root_evaluations + second.total_root_evaluations
    );
    assert_eq!(
        report.evaluated_root_evaluations,
        report.accepted_root_evaluations + full.total_root_evaluations
    );
    close(
        report.exchange_free_energy_dissipation_j,
        first.exchange_free_energy_dissipation_j + second.exchange_free_energy_dissipation_j,
        1e-15,
    );
    close(
        report.surface_diffusion_free_energy_dissipation_j,
        first.surface_diffusion_free_energy_dissipation_j
            + second.surface_diffusion_free_energy_dissipation_j,
        1e-15,
    );
    close(
        report.surface_viscous_dissipation_j,
        first.mechanics.surface_viscous_dissipation_j
            + second.mechanics.surface_viscous_dissipation_j,
        1e-15,
    );
    assert!(report.surface_viscous_dissipation_j > 0.0);
    assert!(report.surface_diffusion_free_energy_dissipation_j > 0.0);
    let displacement: f64 = state
        .q
        .chunks_exact(3)
        .zip(original.q.chunks_exact(3))
        .map(|(b, a)| 0.01 * (b[0] - a[0]))
        .sum();
    close(report.external_work_j, displacement, 1e-15);
    close(
        report.free_energy_change_j,
        first.free_energy_change_j + second.free_energy_change_j,
        1e-13,
    );
    physical_accounts(&report);
}

#[test]
fn moving_first_order_rejection_spends_global_work_and_publishes_nothing() {
    let model = ring();
    let laws = [law(SATURATION)];
    let h = 0.25;
    let transport = ring_transport(&laws, &[3.0], h);
    let original = ring_state(0.0, 0.3, 1e-6);
    let original_bulk = bath();
    let original_amounts = [1e-6];
    let mut full = original.clone();
    let mut full_bulk = original_bulk.clone();
    let mut full_amounts = original_amounts;
    model
        .advance_soluble_surfactant(
            &mut full,
            &mut full_bulk,
            &mut full_amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 3],
            h,
            mechanical(),
            || false,
        )
        .unwrap();
    let mut fine = original.clone();
    let mut fine_bulk = original_bulk.clone();
    let mut fine_amounts = original_amounts;
    for _ in 0..2 {
        model
            .advance_soluble_surfactant(
                &mut fine,
                &mut fine_bulk,
                &mut fine_amounts,
                transport,
                inviscid(),
                &[[0.0; 3]; 3],
                0.5 * h,
                mechanical(),
                || false,
            )
            .unwrap();
    }
    let difference = (full_amounts[0] - fine_amounts[0]).abs();
    assert!(difference > 1e-10);
    let mut limits = loose(h);
    limits.transport.absolute_surface_amount_tolerance_mol = 0.5 * difference;
    limits.transport.minimum_step_s = h;
    let mut state = original.clone();
    let mut bulk = original_bulk.clone();
    let mut amounts = original_amounts;
    let error = model
        .advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 3],
            h,
            mechanical(),
            limits,
            || false,
        )
        .unwrap_err();
    match error {
        AdaptiveSolubleSurfaceDynamicsError::TemporalAccuracy {
            step_s,
            normalized_error,
        } => {
            assert_eq!(step_s, h);
            close(normalized_error, 2.0, 1e-11);
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(state, original);
    assert_eq!(bulk, original_bulk);
    assert_eq!(amounts, original_amounts);
    limits.transport.minimum_step_s = 1e-7;
    limits.maximum_mechanical_steps = 3;
    assert!(matches!(
        model.advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 3],
            h,
            mechanical(),
            limits,
            || false,
        ),
        Err(AdaptiveSolubleSurfaceDynamicsError::MechanicalStepBudget {
            completed: 3,
            required: 3
        })
    ));
    assert_eq!(state, original);
    assert_eq!(bulk, original_bulk);
    assert_eq!(amounts, original_amounts);
}

#[test]
fn moving_duration_rolls_back_private_prefix_on_each_shared_budget_and_late_cancel() {
    let (q, faces, masses, v) = patch();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 2).unwrap();
    let laws = [law(SATURATION); 2];
    let edges = [SurfacePartitionDiffusionEdge {
        faces: [0, 1],
        conductance_mol_s: 2e-5,
    }];
    let transport = SolubleSurfaceTransport {
        laws: &laws,
        kinetics: SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: RT * 0.5,
            desorption_rates_s_inverse: &[2.0, 1.0],
        },
        edges: &edges,
        controls: physical(0.125, 2, 1),
    };
    let original_amounts = [1e-6, 2e-6];
    let original_bulk = bath();
    let original = model
        .initial_state_surfactant(
            0.0,
            &q,
            &v,
            SurfactantField::Uniform {
                law: laws[0],
                amounts_mol: &original_amounts,
            },
            inviscid(),
            &[[0.0; 3]; 4],
        )
        .unwrap();
    let mut state = original.clone();
    let mut bulk = original_bulk.clone();
    let mut amounts = original_amounts;
    let first = model
        .advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 4],
            0.125,
            mechanical(),
            loose(0.125),
            || false,
        )
        .unwrap();
    assert_eq!(first.accepted_intervals, 1);
    for kind in 0..5 {
        let mut limits = loose(0.125);
        let mut selected = transport;
        match kind {
            0 => limits.maximum_mechanical_steps = 3,
            1 => selected.controls.maximum_substeps = 6,
            2 => selected.controls.maximum_surface_pair_exchanges = 12,
            3 => selected.controls.maximum_root_evaluations = first.evaluated_root_evaluations,
            _ => limits.transport.maximum_trials = 1,
        }
        let mut state = original.clone();
        let mut bulk = original_bulk.clone();
        let mut amounts = original_amounts;
        let error = model
            .advance_soluble_surfactant_adaptive_to(
                &mut state,
                &mut bulk,
                &mut amounts,
                selected,
                inviscid(),
                &[[0.0; 3]; 4],
                0.25,
                mechanical(),
                limits,
                || false,
            )
            .unwrap_err();
        match (kind, error) {
            (
                0,
                AdaptiveSolubleSurfaceDynamicsError::MechanicalStepBudget {
                    completed: 3,
                    required: 3,
                },
            )
            | (
                1,
                AdaptiveSolubleSurfaceDynamicsError::ChemicalSubstepBudget {
                    completed: 6,
                    required: 6,
                },
            )
            | (
                2,
                AdaptiveSolubleSurfaceDynamicsError::PairExchangeBudget {
                    completed: 12,
                    required: 12,
                },
            )
            | (
                4,
                AdaptiveSolubleSurfaceDynamicsError::TrialBudget {
                    completed_trials: 1,
                },
            ) => {}
            (
                3,
                AdaptiveSolubleSurfaceDynamicsError::Dynamics(
                    SolubleSurfaceDynamicsError::Transport {
                        error:
                            SolubleSurfaceError::RootBudget {
                                completed_evaluations,
                            },
                        ..
                    },
                ),
            ) => assert_eq!(completed_evaluations, first.evaluated_root_evaluations),
            (_, other) => panic!("budget {kind}: {other:?}"),
        }
        assert_eq!(state, original);
        assert_eq!(bulk, original_bulk);
        assert_eq!(amounts, original_amounts);
    }
    let mut count = 0usize;
    let mut state = original.clone();
    let mut bulk = original_bulk.clone();
    let mut amounts = original_amounts;
    let report = model
        .advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 4],
            0.25,
            mechanical(),
            loose(0.125),
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(report.accepted_intervals, 2);
    let mut seen = 0usize;
    let mut state = original.clone();
    let mut bulk = original_bulk.clone();
    let mut amounts = original_amounts;
    assert!(matches!(
        model.advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 4],
            0.25,
            mechanical(),
            loose(0.125),
            || {
                seen += 1;
                seen == count
            },
        ),
        Err(AdaptiveSolubleSurfaceDynamicsError::Dynamics(
            SolubleSurfaceDynamicsError::Mechanics(SurfaceDynamicsError::Cancelled)
        ))
    ));
    assert_eq!(seen, count);
    assert_eq!(state, original);
    assert_eq!(bulk, original_bulk);
    assert_eq!(amounts, original_amounts);
}

#[test]
fn moving_accuracy_controls_motion_when_all_surface_and_bulk_amounts_are_constant() {
    let model = ring();
    let laws = [law(SATURATION)];
    let duration = 0.6;
    let expected = reference(0.3, 0.0, duration, 8192);
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1.0, 1e-6] {
        let mut controls = loose(duration);
        controls.absolute_position_tolerance_m = tolerance;
        controls.absolute_velocity_tolerance_m_s = tolerance;
        let mut state = ring_state(0.0, 0.3, 1e-6);
        let mut bulk = bath();
        let original_bulk = bulk.clone();
        let mut amounts = [1e-6];
        let report = model
            .advance_soluble_surfactant_adaptive_to(
                &mut state,
                &mut bulk,
                &mut amounts,
                ring_transport(&laws, &[0.0], duration),
                inviscid(),
                &[[0.0; 3]; 3],
                duration,
                mechanical(),
                controls,
                || false,
            )
            .unwrap();
        assert_eq!(bulk, original_bulk);
        assert_eq!(amounts, [1e-6]);
        assert_eq!(report.gross_adsorption_mol, 0.0);
        assert_eq!(report.gross_desorption_mol, 0.0);
        physical_accounts(&report);
        errors.push(
            (state.q[0] - expected[0])
                .abs()
                .max((state.v[0] - expected[1]).abs()),
        );
        work.push(report.evaluated_mechanical_steps);
    }
    assert_eq!(work[0], 3);
    assert!(work[1] > work[0]);
    assert!(errors[1] < 0.02 * errors[0], "{errors:?}");
}

#[test]
fn moving_accuracy_respects_absolute_clock_and_actual_unequal_half_intervals() {
    let model = ring();
    let laws = [law(SATURATION)];
    let start = 2.0_f64.powi(40);
    let ulp = 2.0_f64.powi(-12);
    let h = 3.0 * ulp;
    let end = start + h;
    let middle = start + 0.5 * h;
    assert_eq!(middle - start, 2.0 * ulp);
    assert_eq!(end - middle, ulp);
    let transport = ring_transport(&laws, &[3.0], h);
    let original = ring_state(start, 0.3, 1e-6);
    let original_bulk = bath();
    let original_amounts = [1e-6];
    let mut coarse = original.clone();
    let mut coarse_bulk = original_bulk.clone();
    let mut coarse_amounts = original_amounts;
    model
        .advance_soluble_surfactant(
            &mut coarse,
            &mut coarse_bulk,
            &mut coarse_amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 3],
            h,
            mechanical(),
            || false,
        )
        .unwrap();
    let mut fine = original.clone();
    let mut fine_bulk = original_bulk.clone();
    let mut fine_amounts = original_amounts;
    for step in [middle - start, end - middle] {
        model
            .advance_soluble_surfactant(
                &mut fine,
                &mut fine_bulk,
                &mut fine_amounts,
                transport,
                inviscid(),
                &[[0.0; 3]; 3],
                step,
                mechanical(),
                || false,
            )
            .unwrap();
    }
    let discrepancy = coarse
        .v
        .iter()
        .zip(&fine.v)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(discrepancy > 1e-16);
    let mut controls = loose(h);
    // The actual 2:1 halves give first-order fine correction 1.25, not 1.
    controls.absolute_velocity_tolerance_m_s = 2.0 * 1.25 * discrepancy;
    let mut state = original.clone();
    let mut bulk = original_bulk.clone();
    let mut amounts = original_amounts;
    let report = model
        .advance_soluble_surfactant_adaptive_to(
            &mut state,
            &mut bulk,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 3],
            end,
            mechanical(),
            controls,
            || false,
        )
        .unwrap();
    assert_eq!(state, fine);
    assert_eq!(bulk, fine_bulk);
    assert_eq!(amounts, fine_amounts);
    assert_eq!(state.t, end);
    assert_eq!(state.history[0].h, 2.0 * ulp);
    assert_eq!(state.history[1].h, ulp);
    close(report.maximum_accepted_normalized_error, 0.5, 1e-10);
    assert_eq!(report.accepted_intervals, 1);
    let mut state = original.clone();
    let mut bulk = original_bulk.clone();
    let mut amounts = original_amounts;
    assert!(matches!(model.advance_soluble_surfactant_adaptive_to(
        &mut state,&mut bulk,&mut amounts,transport,inviscid(),&[[0.0;3];3],
        start+ulp,mechanical(),loose(h),||false,
    ),Err(AdaptiveSolubleSurfaceDynamicsError::TimeResolution {time_s}) if time_s==start));
    assert_eq!(state, original);
    assert_eq!(bulk, original_bulk);
    assert_eq!(amounts, original_amounts);
}
