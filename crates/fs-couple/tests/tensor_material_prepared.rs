//! Authored constitutive fixtures and independent physical oracles, not
//! experimental qualification of a specimen or a material-frequency band.
use fs_couple::tensor_material::{
    InitialMemory, TensorMaterialBody, TensorMaterialBudget,
    prepared::{PreparedTensorBody, PreparedTensorControls, PreparedTensorHistory},
};
use fs_exec::CancelGate;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::visco_tensor::{ProjectedTensorProny, TensorProjectionBudget};

fn controls() -> PreparedTensorControls {
    PreparedTensorControls {
        maximum_coordinates: 32,
        maximum_material_points: 128,
        maximum_state_coordinates: 192,
        maximum_heat_groups: 16,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_preparation_multiply_terms: 10_000_000_000,
        maximum_step_multiply_terms: 1_000_000_000,
        maximum_relative_operator_error: 1e-10,
        maximum_work_error_j: 1e-10,
        maximum_steps: 1000,
    }
}

fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}

fn project(law: &TensorProny, columns: &[[f64; 6]], volume: f64) -> ProjectedTensorProny {
    law.project_strain(
        columns,
        volume,
        TensorProjectionBudget {
            max_coordinates: columns.len(),
            max_scalar_terms: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance * (1.0 + expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}, allowance={allowance:.3e}"
    );
}

fn close_slice(actual: &[f64], expected: &[f64], allowance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        close(a, b, allowance);
    }
}

fn snapshot(history: &PreparedTensorHistory<'_>) -> Vec<u64> {
    let mut result = Vec::new();
    result.extend(history.displacement_m().iter().map(|value| value.to_bits()));
    result.extend(history.velocity_mps().iter().map(|value| value.to_bits()));
    result.extend(
        history
            .memory_displacement_m()
            .iter()
            .map(|value| value.to_bits()),
    );
    result.extend(
        [
            history.time_s(),
            history.stored_energy_j(),
            history.supplied_work_j(),
            history.dissipated_heat_j(),
            history.balance_residual_j(),
        ]
        .map(f64::to_bits),
    );
    result.push(history.steps());
    result
}

fn scalar_point(tau: f64, limit: f64) -> ProjectedTensorProny {
    let law = TensorProny::new_mandel(
        diagonal([1.0; 6]),
        vec![PronyBranch::new(diagonal([10.0; 6]), tau).unwrap()],
        limit,
    )
    .unwrap();
    project(&law, &[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]], 1.0)
}

/// m=1, K_eq=1, K_Maxwell=10, tau=1/6. Eliminating the branch
/// force gives u''' + 6u'' + 11u' + 6u = 6F, with roots -1,-2,-3.
/// Heat below integrates the actual branch-force square, independently of H.
fn zener_oracle(
    time: f64,
    initial_u: f64,
    initial_v: f64,
    held_force: f64,
    memory: InitialMemory,
) -> (f64, f64, f64, f64) {
    let initial_branch = if memory == InitialMemory::Relaxed {
        0.0
    } else {
        10.0 * initial_u
    };
    let initial_a = held_force - initial_u - initial_branch;
    let shift = initial_u - held_force;
    let amplitudes = [
        3.0 * shift + 2.5 * initial_v + 0.5 * initial_a,
        -3.0 * shift - 4.0 * initial_v - initial_a,
        shift + 1.5 * initial_v + 0.5 * initial_a,
    ];
    let rates = [1.0, 2.0, 3.0];
    let mut displacement = held_force;
    let mut velocity = 0.0;
    let mut branch_force = 0.0;
    let branch_amplitudes =
        core::array::from_fn::<_, 3, _>(|i| -(1.0 + rates[i] * rates[i]) * amplitudes[i]);
    for i in 0..3 {
        let decay = (-rates[i] * time).exp();
        displacement += amplitudes[i] * decay;
        velocity -= rates[i] * amplitudes[i] * decay;
        branch_force += branch_amplitudes[i] * decay;
    }
    let mut heat = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let rate = rates[i] + rates[j];
            heat +=
                0.6 * branch_amplitudes[i] * branch_amplitudes[j] * -(-rate * time).exp_m1() / rate;
        }
    }
    (displacement, velocity, branch_force / 10.0, heat)
}

#[test]
fn g1_prepared_standard_linear_solid_matches_exact_motion_and_dashpot_heat() {
    let gate = CancelGate::new_clock_free();
    let dt = 0.025;
    let body = PreparedTensorBody::from_mass_matrix(
        1,
        vec![1.0],
        vec![scalar_point(1.0 / 6.0, 1.0)],
        vec![0],
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(body.state_dimension(), 3);
    assert_eq!(body.relaxation_times_s(), &[1.0 / 6.0]);
    let (u0, v0, force) = (0.025, -0.04, 0.012);
    for memory in [InitialMemory::Relaxed, InitialMemory::Instantaneous] {
        let mut history = body.history(&[u0], &[v0], memory, &gate).unwrap();
        let memory0 = if memory == InitialMemory::Relaxed {
            0.0
        } else {
            u0
        };
        let initial_energy = 0.5 * (u0 * u0 + v0 * v0 + 10.0 * memory0 * memory0);
        close(history.stored_energy_j(), initial_energy, 1e-13);
        let mut previous_heat = 0.0;
        let mut previous_u = u0;
        for step_number in 1..=40 {
            let step = history.advance(&[force], &gate).unwrap();
            let time = step_number as f64 * dt;
            let (u, v, spring_memory, heat) = zener_oracle(time, u0, v0, force, memory);
            close(history.displacement_m()[0], u, 2e-12);
            close(history.velocity_mps()[0], v, 2e-12);
            close(history.memory_displacement_m()[0], spring_memory, 2e-12);
            close(history.dissipated_heat_j(), heat, 3e-12);
            close(step.energy.dissipated_j, heat - previous_heat, 3e-12);
            close(step.heat_by_group_j[0], heat - previous_heat, 3e-12);
            close(step.energy.supplied_work_j, force * (u - previous_u), 2e-13);
            close(step.mean_velocity_mps[0], (u - previous_u) / dt, 2e-11);
            close(step.energy.equilibrium_after_j, 0.5 * u * u, 2e-13);
            close(step.energy.kinetic_after_j, 0.5 * v * v, 2e-13);
            close(
                step.energy.memory_after_j,
                5.0 * spring_memory * spring_memory,
                2e-13,
            );
            close(
                history.stored_energy_j() + heat,
                initial_energy + force * (u - u0),
                4e-12,
            );
            close(history.supplied_work_j(), force * (u - u0), 2e-13);
            assert!(step.energy.dissipated_j >= 0.0);
            assert!(step.heat_by_group_j[0] >= 0.0);
            assert_eq!(step.step, step_number);
            close(step.time_s, time, 1e-14);
            previous_heat = heat;
            previous_u = u;
        }
    }
}

#[test]
fn g1_independent_zener_elements_keep_their_own_analytical_heat() {
    let gate = CancelGate::new_clock_free();
    let original = scalar_point(1.0 / 6.0, 1.0);
    let axial = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let points = vec![
        project(original.material(), &[axial, [0.0; 6]], 1.0),
        project(original.material(), &[[0.0; 6], axial], 2.0),
    ];
    // Both original elements share the exact source relaxation time; their
    // volumes and physical inertia differ. Group labels deliberately reverse
    // original point order, so a global-loss redistribution cannot pass.
    let dt = 0.075;
    let body = PreparedTensorBody::from_mass_matrix(
        2,
        vec![1.0, 0.0, 0.0, 2.0],
        points,
        vec![1, 0],
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(body.relaxation_times_s().len(), 1);
    assert_eq!(body.state_dimension(), 6);
    let u0 = [0.025, -0.013];
    let v0 = [-0.04, 0.027];
    let force = [0.012, -0.018];
    let mut history = body
        .history(&u0, &v0, InitialMemory::Instantaneous, &gate)
        .unwrap();
    let mut previous_heat = [0.0; 2];
    for step_number in 1..=8 {
        let step = history.advance(&force, &gate).unwrap();
        let mut expected_heat = [0.0; 2];
        for coordinate in 0..2 {
            let volume = (coordinate + 1) as f64;
            let (u, v, memory, heat) = zener_oracle(
                step_number as f64 * dt,
                u0[coordinate],
                v0[coordinate],
                force[coordinate] / volume,
                InitialMemory::Instantaneous,
            );
            close(history.displacement_m()[coordinate], u, 3e-12);
            close(history.velocity_mps()[coordinate], v, 3e-12);
            close(history.memory_displacement_m()[coordinate], memory, 3e-12);
            expected_heat[coordinate] = volume * heat;
            close(
                step.heat_by_group_j[1 - coordinate],
                volume * heat - previous_heat[coordinate],
                4e-12,
            );
        }
        close(
            history.dissipated_heat_j(),
            expected_heat.iter().sum(),
            5e-12,
        );
        previous_heat = expected_heat;
    }
}

fn coupled_points(times: &[f64]) -> Vec<ProjectedTensorProny> {
    times
        .iter()
        .enumerate()
        .map(|(point, &tau)| {
            let scale = 1.0 + 0.3 * point as f64;
            let mut equilibrium = diagonal([4.0, 6.0, 9.0, 3.0, 5.0, 7.0]);
            equilibrium[0][1] = 0.4;
            equilibrium[1][0] = 0.4;
            let mut branch = diagonal([5.0, 3.0, 7.0, 2.0, 4.0, 6.0]);
            branch[0][1] = -0.7;
            branch[1][0] = -0.7;
            branch = branch.map(|row| row.map(|entry| entry * scale));
            let mut kelvin = diagonal([0.08, 0.19, 0.15, 0.12, 0.06, 0.09]);
            kelvin[0][1] = 0.015;
            kelvin[1][0] = 0.015;
            let law = TensorProny::new_mandel(
                equilibrium,
                vec![PronyBranch::new(branch, tau).unwrap()],
                0.25,
            )
            .unwrap()
            .with_kelvin_viscosity_mandel(kelvin)
            .unwrap();
            project(
                &law,
                &[
                    [1.0, 0.15, 0.0, 0.1, 0.0, 0.0],
                    [-0.2, 0.8, 0.1, 0.0, 0.2, 0.0],
                ],
                0.6 + 0.2 * point as f64,
            )
        })
        .collect()
}

#[test]
fn g1_shared_memory_preserves_consistent_mass_coupling_and_original_local_heat() {
    let gate = CancelGate::new_clock_free();
    let points = coupled_points(&[0.12, 0.37, 0.12]);
    let mass = vec![2.0, 0.35, 0.35, 1.4];
    let native = TensorMaterialBody::from_mass_matrix(
        2,
        mass.clone(),
        points.clone(),
        TensorMaterialBudget {
            max_state_coordinates: 64,
            max_material_points: 4,
        },
        &gate,
    )
    .unwrap();
    let groups = [1, 0, 1];
    let dt = 0.017;
    let body = PreparedTensorBody::from_mass_matrix(
        2,
        mass.clone(),
        points.clone(),
        groups.to_vec(),
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(body.coordinates(), 2);
    assert_eq!(body.state_dimension(), 8);
    assert_eq!(native.state_dimension(), 22);
    assert_eq!(body.mass_matrix_kg(), mass.as_slice());
    assert_eq!(body.relaxation_times_s(), &[0.12, 0.37]);
    assert_eq!(body.heat_group_count(), 2);
    for memory in [InitialMemory::Relaxed, InitialMemory::Instantaneous] {
        let mut reference = native
            .history(&[0.003, -0.002], &[0.011, 0.007], memory, &gate)
            .unwrap();
        let mut history = body
            .history(&[0.003, -0.002], &[0.011, 0.007], memory, &gate)
            .unwrap();
        for step_number in 0..6 {
            let force = [
                0.03 + 0.007 * step_number as f64,
                -0.02 + 0.003 * step_number as f64,
            ];
            let old = reference
                .advance_with_material_heat(&force, dt, 1e-10, &gate)
                .unwrap();
            let step = history.advance(&force, &gate).unwrap();
            close_slice(history.displacement_m(), reference.displacement_m(), 3e-11);
            close_slice(&history.velocity_mps(), &reference.velocity_m_s(), 3e-11);
            close(
                history.stored_energy_j(),
                reference.stored_energy_j(),
                3e-11,
            );
            close(step.energy.supplied_work_j, old.dynamics.supplied, 3e-11);
            close(step.energy.dissipated_j, old.dynamics.dissipated, 3e-11);
            let expected = [
                old.point_heat_j[1],
                old.point_heat_j[0] + old.point_heat_j[2],
            ];
            close_slice(&step.heat_by_group_j, &expected, 3e-11);
            close(
                step.heat_by_group_j.iter().sum(),
                step.energy.dissipated_j,
                3e-11,
            );
            // Every original quadrature spring coordinate is reconstructed,
            // including those with no direct independent mechanical direction.
            let mut offset = 4;
            for point in &points {
                for branch in point.branches() {
                    let family = body
                        .relaxation_times_s()
                        .iter()
                        .position(|tau| tau.to_bits() == branch.relaxation_time_s().to_bits())
                        .unwrap();
                    let y = &history.memory_displacement_m()[2 * family..2 * family + 2];
                    for row in branch.rows() {
                        let z = row.iter().zip(y).map(|(a, y)| a * y).sum::<f64>();
                        close(z, reference.state()[offset], 3e-11);
                        offset += 1;
                    }
                }
            }
            assert_eq!(offset, reference.state().len());
        }
    }
}

#[test]
fn g0_exact_relaxation_times_and_heat_grouping_never_merge_nearby_materials() {
    let gate = CancelGate::new_clock_free();
    let tau = 0.125_f64;
    let next = f64::from_bits(tau.to_bits() + 1);
    for (times, expected) in [([tau, tau], 1), ([tau, next], 2)] {
        let body = PreparedTensorBody::from_mass_matrix(
            2,
            vec![2.0, 0.35, 0.35, 1.4],
            coupled_points(&times),
            vec![0, 1],
            0.013,
            controls(),
            &gate,
        )
        .unwrap();
        assert_eq!(body.relaxation_times_s().len(), expected);
        assert_eq!(body.state_dimension(), 4 + 2 * expected);
        assert_eq!(body.relaxation_times_s()[0].to_bits(), tau.to_bits());
        if expected == 2 {
            assert_eq!(body.relaxation_times_s()[1].to_bits(), next.to_bits());
        }
    }
    let points = coupled_points(&[tau, tau]);
    let separated = PreparedTensorBody::from_mass_matrix(
        2,
        vec![2.0, 0.35, 0.35, 1.4],
        points.clone(),
        vec![0, 1],
        0.013,
        controls(),
        &gate,
    )
    .unwrap();
    let grouped = PreparedTensorBody::from_mass_matrix(
        2,
        vec![2.0, 0.35, 0.35, 1.4],
        points,
        vec![0, 0],
        0.013,
        controls(),
        &gate,
    )
    .unwrap();
    let mut a = separated
        .history(
            &[0.004, -0.002],
            &[0.01, 0.02],
            InitialMemory::Instantaneous,
            &gate,
        )
        .unwrap();
    let mut b = grouped
        .history(
            &[0.004, -0.002],
            &[0.01, 0.02],
            InitialMemory::Instantaneous,
            &gate,
        )
        .unwrap();
    let astep = a.advance(&[0.01, -0.03], &gate).unwrap();
    let bstep = b.advance(&[0.01, -0.03], &gate).unwrap();
    close_slice(a.displacement_m(), b.displacement_m(), 2e-12);
    close_slice(&a.velocity_mps(), &b.velocity_mps(), 2e-12);
    assert!(astep.heat_by_group_j.iter().all(|value| *value > 0.0));
    assert!((astep.heat_by_group_j[0] - astep.heat_by_group_j[1]).abs() > 1e-10);
    close(
        astep.heat_by_group_j.iter().sum(),
        bstep.heat_by_group_j[0],
        3e-12,
    );
}

#[test]
fn g0_prepared_source_capacity_and_memory_scope_refuse_without_substitution() {
    let gate = CancelGate::new_clock_free();
    let point = scalar_point(1.0 / 6.0, 1.0);
    let base = controls();
    let exact = PreparedTensorControls {
        maximum_coordinates: 1,
        maximum_material_points: 1,
        maximum_state_coordinates: 3,
        maximum_heat_groups: 1,
        ..base
    };
    assert!(
        PreparedTensorBody::from_mass_matrix(
            1,
            vec![1.0],
            vec![point.clone()],
            vec![0],
            0.025,
            exact,
            &gate
        )
        .is_ok()
    );
    for insufficient in [
        PreparedTensorControls {
            maximum_coordinates: 0,
            ..exact
        },
        PreparedTensorControls {
            maximum_material_points: 0,
            ..exact
        },
        PreparedTensorControls {
            maximum_state_coordinates: 2,
            ..exact
        },
        PreparedTensorControls {
            maximum_heat_groups: 0,
            ..exact
        },
        PreparedTensorControls {
            maximum_preparation_bytes: 0,
            ..exact
        },
        PreparedTensorControls {
            maximum_preparation_multiply_terms: 0,
            ..exact
        },
        PreparedTensorControls {
            maximum_step_multiply_terms: 0,
            ..exact
        },
    ] {
        assert!(
            PreparedTensorBody::from_mass_matrix(
                1,
                vec![1.0],
                vec![point.clone()],
                vec![0],
                0.025,
                insufficient,
                &gate
            )
            .is_err()
        );
    }
    assert!(
        PreparedTensorBody::from_mass_matrix(
            1,
            vec![-1.0],
            vec![point.clone()],
            vec![0],
            0.025,
            base,
            &gate
        )
        .is_err()
    );
    assert!(
        PreparedTensorBody::from_mass_matrix(
            1,
            vec![1.0],
            vec![point.clone()],
            vec![1],
            0.025,
            base,
            &gate
        )
        .is_err()
    );
    let eigenstrained = point
        .with_eigenstrain_mandel([0.001, 0.0, 0.0, 0.0, 0.0, 0.0])
        .unwrap();
    assert!(
        PreparedTensorBody::from_mass_matrix(
            1,
            vec![1.0],
            vec![eigenstrained],
            vec![0],
            0.025,
            base,
            &gate
        )
        .is_err()
    );
}

#[test]
fn g0_prepared_refusals_preserve_motion_memory_work_heat_and_clock_bitwise() {
    let gate = CancelGate::new_clock_free();
    let body = PreparedTensorBody::from_mass_matrix(
        1,
        vec![1.0],
        vec![scalar_point(1.0 / 6.0, 0.01)],
        vec![0],
        0.025,
        controls(),
        &gate,
    )
    .unwrap();
    assert!(
        body.history(&[0.011], &[0.0], InitialMemory::Relaxed, &gate)
            .is_err()
    );
    let mut history = body
        .history(&[0.002], &[0.004], InitialMemory::Instantaneous, &gate)
        .unwrap();
    history.advance(&[0.01], &gate).unwrap();
    let before = snapshot(&history);
    assert!(history.advance(&[f64::NAN], &gate).is_err());
    assert_eq!(snapshot(&history), before);
    assert!(history.advance(&[], &gate).is_err());
    assert_eq!(snapshot(&history), before);
    // This finite force drives the original material outside its strain domain.
    assert!(history.advance(&[1000.0], &gate).is_err());
    assert_eq!(snapshot(&history), before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(history.advance(&[0.01], &cancelled).is_err());
    assert_eq!(snapshot(&history), before);
    let strict = PreparedTensorControls {
        maximum_work_error_j: 0.0,
        ..controls()
    };
    let strict_body = PreparedTensorBody::from_mass_matrix(
        1,
        vec![1.0],
        vec![scalar_point(1.0 / 6.0, 1.0)],
        vec![0],
        0.13,
        strict,
        &gate,
    )
    .unwrap();
    let mut strict_history = strict_body
        .history(&[0.023], &[0.007], InitialMemory::Instantaneous, &gate)
        .unwrap();
    let strict_before = snapshot(&strict_history);
    assert!(strict_history.advance(&[0.31], &gate).is_err());
    assert_eq!(snapshot(&strict_history), strict_before);
    let once = PreparedTensorBody::from_mass_matrix(
        1,
        vec![1.0],
        vec![scalar_point(1.0 / 6.0, 1.0)],
        vec![0],
        0.025,
        PreparedTensorControls {
            maximum_steps: 1,
            ..controls()
        },
        &gate,
    )
    .unwrap();
    let mut once_history = once
        .history(&[0.002], &[0.004], InitialMemory::Relaxed, &gate)
        .unwrap();
    once_history.advance(&[0.01], &gate).unwrap();
    let accepted = snapshot(&once_history);
    assert!(once_history.advance(&[0.01], &gate).is_err());
    assert_eq!(snapshot(&once_history), accepted);
}
