//! Coupled end flux must participate in the original nonlinear field/circuit
//! history, including terminal polarity and physical versus numerical heat.
use super::*;
use fs_em::winding::{end::EndWindingMode, nonlinear::resistors::GraphResistorMaterial};

#[test]
fn mutual_end_flux_matches_independent_two_loop_evolution_and_survives_graph_rebind() {
    let (magnetic, conductor) = material(false, 18.0);
    let problem = cube(magnetic);
    let profiles = profiles(&problem);
    let e0 = solve(&problem, &profiles[0]).stored_energy_j();
    let e1 = solve(&problem, &profiles[1]).stored_energy_j();
    let es = solve(&problem, &current(&profiles, &[1.0, 1.0])).stored_energy_j();
    // Independent expanded matrix, used only as this modest-scale oracle.
    // End modes are 0.04*[1,-2]*[1,-2]^T + 0.015*[2,1]*[2,1]^T.
    let ends = [0.12, -0.05, -0.05, 0.205];
    let l = [
        2.0 * e0 + ends[0],
        es - e0 - e1 + ends[1],
        es - e0 - e1 + ends[2],
        2.0 * e1 + ends[3],
    ];
    let modes = vec![
        EndWindingMode {
            inductance_h: 0.04,
            turns: vec![1.0, -2.0],
        },
        EndWindingMode {
            inductance_h: 0.015,
            turns: vec![2.0, 1.0],
        },
    ];
    let mut old = [0.15, -0.1];
    let original = SpatialNonlinearWindingCircuit::new_with_end_modes(
        history(&problem, &current(&profiles, &old), vec![None; 48]),
        CircuitGraph {
            node_count: 5,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: 0.7 }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
                (3, 0, Branch::VoltageSource { port: 1 }),
                (3, 4, Branch::Resistor { ohms: 1.1 }),
                (4, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            // Winding order deliberately differs from native branch order.
            branch_indices: vec![5, 2],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: vec![0.02, 0.03],
        },
        vec![NonlinearWindingPort::CurrentInput; 2],
        modes.clone(),
        sources(&profiles, &conductor, 0.5),
        initial(old.to_vec(), vec![0.0, 0.0]),
        config(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(original.end_inductance_h(), [0.02, 0.03]);
    assert_eq!(original.end_inductance().modes(), modes);
    for (&actual, expected) in original.end_inductance_matrix_h().iter().zip(ends) {
        near(actual, expected, 1e-15);
    }
    for w in 0..2 {
        near(
            original.flux_linkages_wb()[w],
            l[2 * w] * old[0] + l[2 * w + 1] * old[1],
            3e-7,
        );
    }

    // This takes the real graph reassembly path: the physical resistor sources
    // are newly bound, so simply retaining the original Layout cannot pass.
    let graph_sources = [(4, 1.1), (1, 0.7)]
        .into_iter()
        .map(|(branch, ohms)| GraphResistorMaterial {
            branch,
            conductor: conductor
                .clone()
                .into_uniform_conductor(18.0 * ohms, 1.0)
                .unwrap(),
        })
        .collect();
    let mut accepted = original
        .with_graph_resistors(graph_sources, &mut || false)
        .unwrap();
    assert_eq!(accepted.end_inductance(), original.end_inductance());
    assert_eq!(
        accepted.field().potential_wb(),
        original.field().potential_wb()
    );
    assert_eq!(accepted.flux_linkages_wb(), original.flux_linkages_wb());
    let rw = [
        resistance(&profiles[0], 18.0, 0.5),
        resistance(&profiles[1], 18.0, 0.5),
    ];
    let r = [rw[0] + 1.1, rw[1] + 0.7];
    let voltage = [-0.6, 1.3];
    let dt = 0.04;
    let energy =
        |i: [f64; 2]| 0.5 * (l[0] * i[0] * i[0] + 2.0 * l[1] * i[0] * i[1] + l[3] * i[1] * i[1]);
    for index in 1..=3 {
        let (a, d) = (l[0] + dt * r[0], l[3] + dt * r[1]);
        let rhs = [
            l[0] * old[0] + l[1] * old[1] + dt * voltage[0],
            l[2] * old[0] + l[3] * old[1] + dt * voltage[1],
        ];
        let determinant = a * d - l[1] * l[2];
        let expected = [
            (d * rhs[0] - l[1] * rhs[1]) / determinant,
            (a * rhs[1] - l[2] * rhs[0]) / determinant,
        ];
        let candidate = accepted
            .step_to_candidate(
                f64::from(index) * dt,
                &[[0.0; 3]; 27],
                &[1.3, -0.6],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        let next = candidate.circuit();
        for w in 0..2 {
            near(next.currents_a()[w], expected[w], 3e-7);
            near(
                next.flux_linkages_wb()[w],
                l[2 * w] * expected[0] + l[2 * w + 1] * expected[1],
                3e-7,
            );
            near(
                step.terminal_voltage_v[w],
                voltage[w] - [1.1, 0.7][w] * expected[w],
                3e-7,
            );
        }
        near(step.stored_before_j, energy(old), 3e-7);
        near(step.stored_after_j, energy(expected), 3e-7);
        near(step.stored_change_j, energy(expected) - energy(old), 3e-7);
        near(
            step.numerical_dissipation_j,
            energy([expected[0] - old[0], expected[1] - old[1]]),
            3e-7,
        );
        near(
            step.heat_j,
            dt * (r[0] * expected[0].powi(2) + r[1] * expected[1].powi(2)),
            3e-7,
        );
        near(
            step.supplied_work_j,
            dt * (voltage[0] * expected[0] + voltage[1] * expected[1]),
            3e-7,
        );
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert_eq!(next.end_inductance(), original.end_inductance());
        assert_eq!(step.graph_materials.len(), 2);
        old = expected;
        accepted = candidate.into_parts().0;
    }
    assert_eq!(original.currents_a(), [0.15, -0.1]);
    assert_eq!(original.time_s(), 0.0);
}

fn current_driven<'m>(
    problem: &'m SpatialNonlinearVectorMagnetostatics,
    profiles: Vec<Vec<[f64; 3]>>,
    conductor: &ResolvedConductorState,
    orientation: f64,
) -> SpatialNonlinearWindingCircuit<'m> {
    SpatialNonlinearWindingCircuit::new_with_end_modes(
        history(problem, &current(&profiles, &[0.0, 0.0]), vec![None; 48]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: f64::NAN }),
                (2, 0, Branch::CurrentSource { port: 1 }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![1, 3],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: vec![0.0, 0.0],
        },
        vec![
            NonlinearWindingPort::VoltageInput {
                initial_voltage_v: 0.0,
                residual_scale_ohm: 1.0,
                current_tolerance_a: 2e-7,
            };
            2
        ],
        vec![EndWindingMode {
            inductance_h: 0.2,
            turns: vec![1.0, -2.0 * orientation],
        }],
        sources(&profiles, conductor, 0.5),
        initial(vec![0.0, 0.0], vec![0.0, 0.0]),
        config(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn saturated_open_secondary_sees_end_mutual_voltage_without_fictitious_heat() {
    let (magnetic, conductor) = material(true, 18.0);
    let problem = cube(magnetic);
    let profiles = profiles(&problem);
    let original = current_driven(&problem, profiles.clone(), &conductor, 1.0);
    let mut reversed_profiles = profiles.clone();
    for density in &mut reversed_profiles[1] {
        *density = density.map(|j| -j);
    }
    let mut reversed = current_driven(&problem, reversed_profiles, &conductor, -1.0);
    let mut accepted = original.clone();
    let rw = [
        resistance(&profiles[0], 18.0, 0.5),
        resistance(&profiles[1], 18.0, 0.5),
    ];
    let mut old_linkage = [0.0, 0.0];
    let mut old_mode = 0.0;
    for (index, (target, imposed)) in [(0.1, [1.0, 0.0]), (0.18, [0.4, -0.2])]
        .into_iter()
        .enumerate()
    {
        let dt = target - accepted.time_s();
        let direct = solve(&problem, &current(&profiles, &imposed));
        let mode = imposed[0] - 2.0 * imposed[1];
        let linkage = [
            pair(&problem, direct.potential_wb(), &profiles[0]) + 0.2 * mode,
            pair(&problem, direct.potential_wb(), &profiles[1]) - 0.4 * mode,
        ];
        // These are the same original-source refresh entry points used by
        // thermal compositions; neither may lose the retained end modes.
        let candidate = if index == 0 {
            accepted.step_with_conductors_candidate(
                target,
                &[[0.0; 3]; 27],
                vec![None; 48],
                accepted.conductors().to_vec(),
                &imposed,
                controls(),
                &mut || false,
            )
        } else {
            accepted.step_with_materials_candidate(
                target,
                &[[0.0; 3]; 27],
                problem.materials().to_vec(),
                vec![None; 48],
                accepted.conductors().to_vec(),
                &imposed,
                controls(),
                &mut || false,
            )
        }
        .unwrap();
        let reverse = reversed
            .step_to_candidate(
                target,
                &[[0.0; 3]; 27],
                &[imposed[0], -imposed[1]],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        for w in 0..2 {
            near(candidate.circuit().currents_a()[w], imposed[w], 2e-7);
            near(candidate.circuit().flux_linkages_wb()[w], linkage[w], 3e-7);
            let expected_voltage = rw[w] * imposed[w] + (linkage[w] - old_linkage[w]) / dt;
            near(step.terminal_voltage_v[w], expected_voltage, 3e-7);
            near(
                reverse.step().terminal_voltage_v[w],
                [1.0, -1.0][w] * expected_voltage,
                3e-7,
            );
        }
        if index == 0 {
            // A zero-current secondary still receives the shared end-flux EMF.
            assert!(step.terminal_voltage_v[1].abs() > 1.0);
            assert!(
                direct
                    .flux_density_t()
                    .iter()
                    .any(|&b| dot(b, b).sqrt() > 0.01)
            );
        }
        near(
            step.stored_after_j,
            direct.stored_energy_j() + 0.1 * mode * mode,
            3e-7,
        );
        near(
            step.numerical_dissipation_j,
            step.field.numerical_dissipation_j + 0.1 * (mode - old_mode).powi(2),
            3e-7,
        );
        near(
            step.heat_j,
            dt * (rw[0] * imposed[0].powi(2) + rw[1] * imposed[1].powi(2)),
            3e-7,
        );
        near(step.heat_j, reverse.step().heat_j, 3e-7);
        near(step.stored_after_j, reverse.step().stored_after_j, 3e-7);
        near(
            step.numerical_dissipation_j,
            reverse.step().numerical_dissipation_j,
            3e-7,
        );
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert_eq!(
            candidate.circuit().end_inductance(),
            original.end_inductance()
        );
        old_linkage = linkage;
        old_mode = mode;
        accepted = candidate.into_parts().0;
        reversed = reverse.into_parts().0;
    }
    assert_eq!(original.currents_a(), [0.0, 0.0]);
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn late_cancelled_shared_end_candidate_preserves_original_field_and_circuit() {
    let (magnetic, conductor) = material(false, 18.0);
    let problem = cube(magnetic);
    let original = current_driven(&problem, profiles(&problem), &conductor, 1.0);
    let original_field = original.field().potential_wb().to_vec();
    let original_state = original.state().to_vec();
    let original_linkage = original.flux_linkages_wb().to_vec();
    let mut polls = 0;
    let success = original
        .step_to_candidate(0.1, &[[0.0; 3]; 27], &[0.7, -0.3], controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 20);
    let mut calls = 0;
    assert!(matches!(
        original.step_to_candidate(0.1, &[[0.0; 3]; 27], &[0.7, -0.3], controls(), &mut || {
            calls += 1;
            calls == polls - 1
        },),
        Err(SpatialNonlinearWindingError::Cancelled)
    ));
    assert_eq!(original.field().potential_wb(), original_field);
    assert_eq!(original.state(), original_state);
    assert_eq!(original.flux_linkages_wb(), original_linkage);
    assert_eq!(original.currents_a(), [0.0, 0.0]);
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(
        original.end_inductance(),
        success.circuit().end_inductance()
    );
}
