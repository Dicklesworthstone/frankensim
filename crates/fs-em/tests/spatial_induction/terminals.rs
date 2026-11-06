//! Physical terminal currents and work in the original coupled A–phi field.
use super::*;
#[path = "terminals/circuit.rs"]
mod circuit;
#[path = "terminals/current.rs"]
mod current;
#[cfg(feature = "magnetic-hysteresis")]
#[path = "terminals/hysteretic.rs"]
mod hysteretic;
#[path = "terminals/nonlinear.rs"]
mod nonlinear;
use fs_em::spatial::vector_potential::induction::{
    mechanical::SpatialLorentzLoadLimits, terminals::SpatialTerminalInduction,
};

fn electrodes(problem: &SpatialVectorMagnetostatics) -> Vec<Vec<usize>> {
    [0.0, 1.0]
        .into_iter()
        .map(|x| {
            problem
                .exterior_faces()
                .into_iter()
                .filter(|&f| {
                    problem.mesh().faces[f]
                        .iter()
                        .all(|&v| problem.positions_m()[v as usize][0] == x)
                })
                .collect()
        })
        .collect()
}

#[test]
fn voltage_drives_real_induction_and_exact_terminal_work_then_recovers_dc_resistance() {
    // Six tetrahedra, one interior body-diagonal Whitney edge W. Independent
    // integration gives int W_x=1/6, int |W|²=1/5 and K=4/mu. All electric
    // nodes lie on x=0 or x=1 electrodes, hence phi=-U*x. The actual equation
    // is M a'+K a=sigma U/6, M=sigma/5, not a static resistor approximation.
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
    let mut state =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let (sigma, k, mass, u, dt): (f64, f64, f64, f64, f64) = (5.0, 2.0, 1.0, 1.2, 0.1);
    let mut exact = 0.0;
    for _ in 0..120 {
        let old = exact;
        exact = ((mass - 0.5 * k * dt) * old + dt * sigma * u / 6.0) / (mass + 0.5 * k * dt);
        let rate = (exact - old) / dt;
        let incoming = sigma * (u - rate / 6.0);
        let expected_heat = dt * sigma * (u * u - u * rate / 3.0 + rate * rate / 5.0);
        let candidate = state
            .step_voltage_candidate(&[u, 0.0], &zero, &velocity, dt, controls(), &mut || false)
            .unwrap();
        near(
            candidate.state().potential_wb()[diagonal(&problem)],
            exact,
            3e-9,
        );
        near(
            candidate.step().terminals[0].incoming_current_a,
            incoming,
            3e-9,
        );
        near(
            candidate.step().terminals[1].incoming_current_a,
            -incoming,
            3e-9,
        );
        near(candidate.step().terminal_work_j, dt * u * incoming, 3e-9);
        near(candidate.step().heat_j, expected_heat, 3e-9);
        near(
            candidate.step().stored_after_j,
            0.5 * k * exact * exact,
            3e-9,
        );
        near(candidate.step().energy_residual_j, 0.0, 3e-9);
        near(candidate.step().impressed_work_j, 0.0, 1e-15);
        let accepted_time = state.time_s();
        let (next, _) = candidate.into_parts();
        near(state.time_s(), accepted_time, 1e-15);
        state = next;
    }
    near(
        state.potential_wb()[diagonal(&problem)],
        sigma * u / (6.0 * k),
        1e-9,
    );
    let dc = state
        .step_voltage_candidate(&[u, 0.0], &zero, &velocity, dt, controls(), &mut || false)
        .unwrap();
    near(
        u / dc.step().terminals[0].incoming_current_a,
        1.0 / sigma,
        1e-9,
    );
    near(dc.step().heat_j / dt, sigma * u * u, 1e-8);
}

#[test]
fn voltage_transient_has_second_order_time_refinement_against_independent_rl_solution() {
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
    let total: f64 = 0.4;
    let exact = 5.0 / 12.0 * (1.0 - (-2.0 * total).exp());
    let mut errors = vec![];
    for count in [4, 8, 16] {
        let native =
            SpatialInduction::new(&initial, sources.clone(), preparation(), 2e-9, &mut || {
                false
            })
            .unwrap();
        let mut state = SpatialTerminalInduction::new(
            native,
            &electrodes(&problem),
            preparation(),
            &mut || false,
        )
        .unwrap();
        for _ in 0..count {
            state = state
                .step_voltage_candidate(
                    &[1.0, 0.0],
                    &zero,
                    &velocity,
                    total / count as f64,
                    controls(),
                    &mut || false,
                )
                .unwrap()
                .into_parts()
                .0;
        }
        errors.push((state.potential_wb()[diagonal(&problem)] - exact).abs());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (order - 2.0).abs() < 0.05,
            "terminal transient orders: {errors:?}"
        );
    }
}

#[test]
fn free_charge_nodes_preserve_terminal_balance_and_exact_common_mode_invariance() {
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
    let state =
        SpatialTerminalInduction::new(native, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let run = |v: [f64; 2]| {
        state
            .step_voltage_candidate(&v, &zero, &velocity, 0.07, controls(), &mut || false)
            .unwrap()
    };
    let a = run([1.0, 0.0]);
    let b = run([1073741825.0, 1073741824.0]);
    assert_eq!(a.state().potential_wb(), b.state().potential_wb());
    assert_eq!(a.step().terminal_work_j, b.step().terminal_work_j);
    assert_eq!(a.step().terminals, b.step().terminals);
    near(
        a.step()
            .terminals
            .iter()
            .map(|t| t.incoming_current_a)
            .sum(),
        0.0,
        2e-8,
    );
    assert!(a.step().terminals[0].incoming_current_a > 0.0);
    assert!(a.step().stored_after_j > 0.0 && a.step().heat_j > 0.0);
    assert!(a.step().maximum_ampere_residual_a < controls().current_tolerance_a);
    assert!(a.step().maximum_charge_residual_a < controls().current_tolerance_a);
}

#[test]
fn motion_generates_native_short_circuit_current_and_work_and_open_voltage_cancels_it() {
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
            .unwrap();
    let open = state
        .step_voltage_candidate(&[0.0, 0.7], &zero, &velocity, 0.03, controls(), &mut || {
            false
        })
        .unwrap();
    for t in &open.step().terminals {
        near(t.incoming_current_a, 0.0, 2e-8);
    }
    near(open.step().heat_j, 0.0, 1e-15);
    let short = state
        .step_voltage_candidate(&[0.0, 0.0], &zero, &velocity, 0.03, controls(), &mut || {
            false
        })
        .unwrap();
    assert!(short.step().terminals[0].incoming_current_a > 0.0);
    assert!(short.step().heat_j > 0.0 && short.step().mechanical_work_j < 0.0);
    near(short.step().terminal_work_j, 0.0, 1e-15);
    let cells: Vec<_> = (0..problem.mesh().tets.len()).collect();
    let loads = short
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
        short.step().mechanical_work_j,
        1e-10,
    );
    near(short.step().energy_residual_j, 0.0, 2e-8);
}

#[test]
fn terminal_geometry_budgets_bad_drives_and_cancellation_refuse_without_advancing_history() {
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
    let faces = electrodes(&problem);
    assert!(
        SpatialTerminalInduction::new(native.clone(), &[], preparation(), &mut || false).is_err()
    );
    assert!(
        SpatialTerminalInduction::new(
            native.clone(),
            &[faces[0].clone(), faces[0].clone()],
            preparation(),
            &mut || false
        )
        .is_err()
    );
    let mut limits = preparation();
    limits.maximum_unknowns = 0;
    assert!(SpatialTerminalInduction::new(native.clone(), &faces, limits, &mut || false).is_err());
    // Distinct fixed-A trace components leave a relative gauge that changes
    // their scalar voltage difference. Refuse an arbitrary tree-cut voltage.
    let mut isolated_trace = trace(&problem, [0.0; 3], false);
    isolated_trace
        .potential
        .retain(|t| faces.iter().any(|set| set.contains(&t.face)));
    let isolated_initial = problem
        .solve(&zero, &isolated_trace, static_controls(), &mut || false)
        .unwrap();
    let isolated = SpatialInduction::new(
        &isolated_initial,
        native.conductors().to_vec(),
        preparation(),
        2e-9,
        &mut || false,
    )
    .unwrap();
    assert!(SpatialTerminalInduction::new(isolated, &faces, preparation(), &mut || false).is_err());
    let state =
        SpatialTerminalInduction::new(native, &faces, preparation(), &mut || false).unwrap();
    let original = state.potential_wb().to_vec();
    for volts in [vec![1.0], vec![f64::NAN, 0.0], vec![f64::MAX, -f64::MAX]] {
        assert!(
            state
                .step_voltage_candidate(&volts, &zero, &velocity, 0.1, controls(), &mut || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        state.step_voltage_candidate(&[1.0, 0.0], &zero, &velocity, 0.1, controls(), &mut || {
            polls += 1;
            polls > 20
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    assert_eq!(state.potential_wb(), original);
    assert_eq!(state.time_s(), 0.0);
}
