//! G0/G1: total-current excitation remains an equipotential field problem.
use super::*;

#[test]
fn current_electrodes_solve_inductive_voltage_and_the_independent_rl_transient() {
    // Independent six-tet cube reduction: int W_x=1/6, int |W|^2=1/5.
    // Eliminating the UNKNOWN electrode voltage gives
    // sigma*(1/5-1/36)*a' + K*a = I/6, U=I/sigma+a'/6.
    // The 1/36 term specifically detects missing grouped-node cross terms.
    let (problem, sources) = cube(1, None);
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let mut state = voltage
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let (sigma, mass, k, incoming, dt): (f64, f64, f64, f64, f64) =
        (5.0, 31.0 / 36.0, 2.0, 1.7, 0.125);
    let mut exact = 0.0;
    for _ in 0..96 {
        let old = exact;
        exact = ((mass - 0.5 * dt * k) * old + dt * incoming / 6.0) / (mass + 0.5 * dt * k);
        let rate = (exact - old) / dt;
        let u = incoming / sigma + rate / 6.0;
        let candidate = state
            .step_current_candidate(
                &[incoming, -incoming],
                &zero,
                &velocity,
                dt,
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        near(
            candidate.state().potential_wb()[diagonal(&problem)],
            exact,
            3e-9,
        );
        near(step.terminals[0].incoming_current_a, incoming, 3e-9);
        near(step.terminals[1].incoming_current_a, -incoming, 3e-9);
        near(step.terminals[0].voltage_to_reference_v, 0.0, 1e-15);
        near(step.terminals[1].voltage_to_reference_v, -u, 3e-9);
        near(step.terminal_work_j, dt * incoming * u, 3e-9);
        near(
            step.heat_j,
            dt * sigma * (u * u - u * rate / 3.0 + rate * rate / 5.0),
            3e-9,
        );
        near(step.stored_after_j, 0.5 * k * exact * exact, 3e-9);
        near(step.energy_residual_j, 0.0, 3e-9);
        state = candidate.into_parts().0;
    }
    near(
        state.potential_wb()[diagonal(&problem)],
        incoming / (6.0 * k),
        1e-9,
    );
    let dc = state
        .step_current_candidate(
            &[incoming, -incoming],
            &zero,
            &velocity,
            dt,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(
        -dc.step().terminals[1].voltage_to_reference_v / incoming,
        1.0 / sigma,
        1e-9,
    );
    let old_a = state.potential_wb().to_vec();
    let old_time = state.time_s();
    let voltage = state
        .into_voltage_driven(preparation(), &mut || false)
        .unwrap();
    assert_eq!(voltage.potential_wb(), old_a);
    assert_eq!(voltage.time_s(), old_time);
}

#[test]
fn voltage_and_current_control_recover_the_same_refined_mesh_state_and_electrode_profile() {
    let (problem, sources) = cube(2, None);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], true),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let expected = voltage
        .step_voltage_candidate(
            &[1.3, -0.2],
            &zero,
            &velocity,
            0.08,
            controls(),
            &mut || false,
        )
        .unwrap();
    let prescribed: Vec<_> = expected
        .step()
        .terminals
        .iter()
        .map(|t| t.incoming_current_a)
        .collect();
    let current = voltage
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let actual = current
        .step_current_candidate(&prescribed, &zero, &velocity, 0.08, controls(), &mut || {
            false
        })
        .unwrap();
    for (&a, &b) in actual
        .state()
        .potential_wb()
        .iter()
        .zip(expected.state().potential_wb())
    {
        near(a, b, 2e-8);
    }
    near(
        actual.step().terminals[1].voltage_to_reference_v,
        -1.5,
        2e-8,
    );
    near(actual.step().heat_j, expected.step().heat_j, 2e-8);
    near(
        actual.step().terminal_work_j,
        expected.step().terminal_work_j,
        2e-8,
    );
    for (a, b) in actual.step().cells.iter().zip(&expected.step().cells) {
        for (u, v) in a
            .current_density_a_m2
            .iter()
            .flatten()
            .zip(b.current_density_a_m2.iter().flatten())
        {
            near(*u, *v, 2e-7);
        }
        for (local, &vertex) in problem.mesh().tets[a.cell].iter().enumerate() {
            let x = problem.positions_m()[vertex as usize][0];
            if x == 0.0 {
                near(a.potential_v[local], 0.0, 1e-15);
            }
            if x == 1.0 {
                near(a.potential_v[local], -1.5, 2e-8);
            }
        }
    }
}

#[test]
fn open_current_electrodes_recover_motional_emf_and_generator_work_remains_separate() {
    let (problem, sources) = cube(2, None);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0, 0.7, 0.0]; problem.mesh().vertex_count];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0, 0.0, 1.0], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let state =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
    let open = state
        .step_current_candidate(&[0.0, 0.0], &zero, &velocity, 0.03, controls(), &mut || {
            false
        })
        .unwrap();
    near(open.step().terminals[1].voltage_to_reference_v, 0.7, 2e-8);
    near(open.step().heat_j, 0.0, 1e-14);
    near(open.step().terminal_work_j, 0.0, 2e-9);
    let generation = state
        .step_current_candidate(
            &[1.0, -1.0],
            &zero,
            &velocity,
            0.03,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(generation.step().terminal_work_j < 0.0);
    assert!(generation.step().mechanical_work_j < 0.0);
    assert!(generation.step().heat_j > 0.0);
    let cells: Vec<_> = (0..problem.mesh().tets.len()).collect();
    let loads = generation
        .lorentz_loads(
            &cells,
            SpatialLorentzLoadLimits {
                max_cells: cells.len(),
                max_vertices: problem.mesh().vertex_count,
            },
            &mut || false,
        )
        .unwrap();
    near(
        0.03 * loads.power_w(&velocity, &mut || false).unwrap(),
        generation.step().mechanical_work_j,
        2e-9,
    );
    near(generation.step().energy_residual_j, 0.0, 2e-8);
}

#[test]
fn global_current_balance_cannot_join_insulator_separated_conductors() {
    let (positions, tets) = grid(3);
    let (metal, conductor) = source(2.0, Some(5.0), None);
    let (air, _) = source(2.0, None, None);
    let mut materials = vec![];
    let mut sources = vec![];
    for tet in &tets {
        let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
        let conducting = x < 1.0 / 3.0 || x > 2.0 / 3.0;
        materials.push(if conducting {
            metal.clone()
        } else {
            air.clone()
        });
        sources.push(if conducting { conductor.clone() } else { None });
    }
    let problem =
        SpatialVectorMagnetostatics::new(positions, tets, materials, assembly(), &mut || false)
            .unwrap();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let state =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
    // These two electrodes touch distinct conductors. Their globally zero sum
    // cannot authorize a current through the intervening insulating material.
    assert!(
        state
            .step_current_candidate(&[1.0, -1.0], &zero, &velocity, 0.1, controls(), &mut || {
                false
            })
            .is_err()
    );
    let open = state
        .step_current_candidate(&[0.0, 0.0], &zero, &velocity, 0.1, controls(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(open.step().heat_j, 0.0);
}

#[test]
fn current_preparation_drives_and_cancellation_fail_without_changing_history() {
    let (problem, sources) = cube(1, None);
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let mut limits = preparation();
    limits.maximum_unknowns = 0;
    assert!(
        voltage
            .clone()
            .into_current_driven(limits, &mut || false)
            .is_err()
    );
    limits = preparation();
    limits.maximum_assembly_terms = 0;
    assert!(
        voltage
            .clone()
            .into_current_driven(limits, &mut || false)
            .is_err()
    );
    assert!(matches!(
        voltage
            .clone()
            .into_current_driven(preparation(), &mut || true),
        Err(SpatialMagneticError::Cancelled)
    ));
    let state = voltage
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let original = state.potential_wb().to_vec();
    for inputs in [vec![1.0], vec![1.0, 0.0], vec![f64::NAN, 0.0]] {
        assert!(
            state
                .step_current_candidate(&inputs, &zero, &velocity, 0.1, controls(), &mut || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        state.step_current_candidate(&[1.0, -1.0], &zero, &velocity, 0.1, controls(), &mut || {
            polls += 1;
            polls > 20
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    assert_eq!(state.potential_wb(), original);
    assert_eq!(state.time_s(), 0.0);
}
