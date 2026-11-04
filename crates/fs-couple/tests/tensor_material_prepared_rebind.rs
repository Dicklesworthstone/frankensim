//! Prescribed material changes with independent motion, heat and work oracles.
//! Authored constitutive fixtures, not a closed thermomechanical model.
use fs_couple::tensor_material::{
    InitialMemory,
    prepared::{
        PreparedTensorBody, PreparedTensorControls, PreparedTensorHistory,
        PreparedTensorMaterialChangeBudget,
    },
};
use fs_exec::CancelGate;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::visco_tensor::{ProjectedTensorProny, TensorProjectionBudget};

const DT: f64 = 0.025;

fn controls() -> PreparedTensorControls {
    PreparedTensorControls {
        maximum_coordinates: 1,
        maximum_material_points: 1,
        maximum_state_coordinates: 4,
        maximum_heat_groups: 1,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_preparation_multiply_terms: 10_000_000_000,
        maximum_step_multiply_terms: 100_000,
        maximum_relative_operator_error: 1e-10,
        maximum_work_error_j: 1e-10,
        maximum_steps: 100,
    }
}

fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}

fn project(law: &TensorProny, strain_scale: f64) -> ProjectedTensorProny {
    law.project_strain(
        &[[strain_scale, 0.0, 0.0, 0.0, 0.0, 0.0]],
        1.0,
        TensorProjectionBudget {
            max_coordinates: 1,
            max_scalar_terms: 10_000,
        },
        || false,
    )
    .unwrap()
}

fn point(equilibrium: f64, branches: &[(f64, f64)], viscosity: f64) -> ProjectedTensorProny {
    let law = TensorProny::new_mandel(
        diagonal([equilibrium; 6]),
        branches
            .iter()
            .map(|&(stiffness, tau)| PronyBranch::new(diagonal([stiffness; 6]), tau).unwrap())
            .collect(),
        1.0,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(diagonal([viscosity; 6]))
    .unwrap();
    project(&law, 1.0)
}

fn body(point: ProjectedTensorProny, mass: f64, dt: f64, gate: &CancelGate) -> PreparedTensorBody {
    PreparedTensorBody::from_mass_matrix(1, vec![mass], vec![point], vec![0], dt, controls(), gate)
        .unwrap()
}

fn allowance(work_j: f64) -> PreparedTensorMaterialChangeBudget {
    PreparedTensorMaterialChangeBudget {
        maximum_absolute_parameter_work_j: work_j,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}, allowance={tolerance:e}"
    );
}

fn snapshot(history: &PreparedTensorHistory<'_>) -> Vec<u64> {
    let mut values: Vec<_> = history
        .state()
        .iter()
        .chain(history.velocity_mps())
        .map(|x| x.to_bits())
        .collect();
    values.extend(
        [
            history.time_s(),
            history.stored_energy_j(),
            history.supplied_work_j(),
            history.parameter_work_j(),
            history.dissipated_heat_j(),
            history.balance_residual_j(),
        ]
        .map(f64::to_bits),
    );
    values.push(history.steps());
    values
}

fn same_instant(previous: &PreparedTensorHistory<'_>, next: &PreparedTensorHistory<'_>) {
    assert_eq!(previous.displacement_m(), next.displacement_m());
    assert_eq!(previous.velocity_mps(), next.velocity_mps());
    assert_eq!(previous.state()[1].to_bits(), next.state()[1].to_bits()); // momentum
    assert_eq!(previous.time_s().to_bits(), next.time_s().to_bits());
    assert_eq!(previous.steps(), next.steps());
    assert_eq!(
        previous.supplied_work_j().to_bits(),
        next.supplied_work_j().to_bits()
    );
    assert_eq!(
        previous.dissipated_heat_j().to_bits(),
        next.dissipated_heat_j().to_bits()
    );
}

/// After a change to m=1, K_eq=1, K_b=10, tau=1/6 and zero Kelvin/force,
/// u'''+6u''+11u'+6u=0 has roots -1,-2,-3. The initial spring displacement
/// is arbitrary accepted history, not an assumed relaxed/instantaneous reset.
/// Heat integrates K_b*y^2/tau independently of the energy balance.
fn analytical(t: f64, u0: f64, v0: f64, y0: f64) -> [f64; 4] {
    let a0 = -u0 - 10.0 * y0;
    let amplitudes = [
        3.0 * u0 + 2.5 * v0 + 0.5 * a0,
        -3.0 * u0 - 4.0 * v0 - a0,
        u0 + 1.5 * v0 + 0.5 * a0,
    ];
    let rates: [f64; 3] = [1.0, 2.0, 3.0];
    let spring: [f64; 3] =
        core::array::from_fn(|i| -(1.0 + rates[i] * rates[i]) * amplitudes[i] / 10.0);
    let (mut u, mut v, mut y, mut heat) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..3 {
        let decay = (-rates[i] * t).exp();
        u += amplitudes[i] * decay;
        v -= rates[i] * amplitudes[i] * decay;
        y += spring[i] * decay;
        for j in 0..3 {
            let rate = rates[i] + rates[j];
            heat += 60.0 * spring[i] * spring[j] * -(-rate * t).exp_m1() / rate;
        }
    }
    [u, v, y, heat]
}

#[test]
fn g1_material_changes_preserve_memory_and_exact_motion_with_signed_parameter_work() {
    let gate = CancelGate::new_clock_free();
    let original = body(point(0.25, &[(2.0, 0.4)], 0.05), 1.0, DT, &gate);
    let stiffer = body(point(4.0, &[(20.0, 0.3)], 0.1), 1.0, DT, &gate);
    // The new branch is rank one, unlike the original full-rank factors.
    // Complete compatible displacement memory supplies its physical strain.
    let target_law = TensorProny::new_mandel(
        diagonal([1.0; 6]),
        vec![PronyBranch::new(diagonal([10.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 1.0 / 6.0).unwrap()],
        1.0,
    )
    .unwrap();
    let target = PreparedTensorBody::from_mass_matrix_affine(
        1,
        vec![1.0],
        vec![project(&target_law, 1.0)],
        vec![0],
        DT,
        0.1,
        controls(),
        &gate,
    )
    .unwrap();
    let (u_initial, v_initial) = (0.025, -0.04);
    let initial_energy = 0.5 * (v_initial * v_initial + 2.25 * u_initial * u_initial);
    let mut history = original
        .history(
            &[u_initial],
            &[v_initial],
            InitialMemory::Instantaneous,
            &gate,
        )
        .unwrap();
    for _ in 0..3 {
        history.advance(&[0.012], &gate).unwrap();
    }
    assert!(history.dissipated_heat_j() > 0.0);
    assert!(history.supplied_work_j().abs() > 0.0);
    let before = snapshot(&history);
    let (u, y) = (
        history.displacement_m()[0],
        history.memory_displacement_m()[0],
    );
    let equilibrium_up = 0.5 * (4.0 - 0.25) * u * u;
    let memory_up = 0.5 * (20.0 - 2.0) * y * y;
    let up = history
        .rebind_materials(&stiffer, allowance(1.0), &gate)
        .unwrap();
    assert!(up.parameter_work_j > 0.0);
    close(up.equilibrium_change_j, equilibrium_up, 2e-14);
    close(up.memory_change_j, memory_up, 2e-14);
    close(up.parameter_work_j, equilibrium_up + memory_up, 2e-14);
    same_instant(&history, &up.history);
    assert_eq!(
        history.memory_displacement_m(),
        up.history.memory_displacement_m()
    );
    assert_eq!(snapshot(&history), before);
    let mut changed = up.history;
    let held = changed.advance(&[-0.009], &gate).unwrap();
    close(
        held.energy.cumulative_parameter_work_j,
        equilibrium_up + memory_up,
        2e-14,
    );
    close(changed.balance_residual_j(), 0.0, 2e-12);

    let (u, v, y) = (
        changed.displacement_m()[0],
        changed.velocity_mps()[0],
        changed.memory_displacement_m()[0],
    );
    let equilibrium_down = 0.5 * (1.0 - 4.0) * u * u;
    let memory_down = 0.5 * (10.0 - 20.0) * y * y;
    let down = changed
        .rebind_materials(&target, allowance(1.0), &gate)
        .unwrap();
    assert!(down.parameter_work_j < 0.0);
    close(down.equilibrium_change_j, equilibrium_down, 2e-14);
    close(down.memory_change_j, memory_down, 2e-14);
    close(down.parameter_work_j, equilibrium_down + memory_down, 2e-14);
    same_instant(&changed, &down.history);
    assert_eq!(
        changed.memory_displacement_m(),
        down.history.memory_displacement_m()
    );
    let parameter_work = equilibrium_up + memory_up + equilibrium_down + memory_down;
    let (heat_before, work_before, time_before, steps_before) = (
        changed.dissipated_heat_j(),
        changed.supplied_work_j(),
        changed.time_s(),
        changed.steps(),
    );
    let mut rebound = down.history;
    let mut previous_heat = 0.0;
    for step in 1..=20 {
        let report = rebound.advance_affine(&[0.0], &[0.0], &gate).unwrap();
        let expected = analytical(step as f64 * DT, u, v, y);
        close(rebound.displacement_m()[0], expected[0], 3e-12);
        close(rebound.velocity_mps()[0], expected[1], 3e-12);
        close(rebound.memory_displacement_m()[0], expected[2], 3e-12);
        close(
            rebound.dissipated_heat_j(),
            heat_before + expected[3],
            3e-12,
        );
        close(
            report.motion.heat_by_group_j[0],
            expected[3] - previous_heat,
            3e-12,
        );
        close(rebound.parameter_work_j(), parameter_work, 2e-14);
        close(
            report.motion.energy.cumulative_parameter_work_j,
            parameter_work,
            2e-14,
        );
        close(rebound.supplied_work_j(), work_before, 2e-14);
        close(
            rebound.stored_energy_j() + rebound.dissipated_heat_j()
                - rebound.supplied_work_j()
                - rebound.parameter_work_j(),
            initial_energy,
            5e-12,
        );
        close(rebound.balance_residual_j(), 0.0, 5e-12);
        close(rebound.time_s(), time_before + step as f64 * DT, 1e-14);
        assert_eq!(rebound.steps(), steps_before + step);
        assert!(report.motion.heat_by_group_j[0] >= 0.0);
        previous_heat = expected[3];
    }
}

#[test]
fn g1_relaxation_groups_split_without_reset_and_refuse_unequal_history_merge() {
    let gate = CancelGate::new_clock_free();
    let shared = body(point(1.0, &[(2.0, 0.2), (3.0, 0.2)], 0.0), 1.0, DT, &gate);
    let split = body(point(1.0, &[(2.0, 0.05), (3.0, 0.5)], 0.0), 1.0, DT, &gate);
    let mut history = shared
        .history(&[0.03], &[0.015], InitialMemory::Instantaneous, &gate)
        .unwrap();
    history.advance(&[0.01], &gate).unwrap();
    let before = snapshot(&history);
    let y = history.memory_displacement_m()[0];
    let changed = history
        .rebind_materials(&split, allowance(0.0), &gate)
        .unwrap();
    same_instant(&history, &changed.history);
    assert_eq!(changed.history.memory_displacement_m(), &[y, y]);
    assert_eq!(changed.parameter_work_j, 0.0);
    assert_eq!(snapshot(&history), before);
    // A same-instant reverse merge is exact because both histories coincide.
    let reverse = changed
        .history
        .rebind_materials(&shared, allowance(0.0), &gate)
        .unwrap();
    assert_eq!(reverse.history.state(), history.state());
    same_instant(&history, &reverse.history);
    let mut evolving = changed.history;
    evolving.advance(&[0.0], &gate).unwrap();
    assert!(
        (evolving.memory_displacement_m()[0] - evolving.memory_displacement_m()[1]).abs() > 1e-5
    );
    let accepted = snapshot(&evolving);
    assert!(
        evolving
            .rebind_materials(&shared, allowance(1.0), &gate)
            .is_err()
    );
    assert_eq!(snapshot(&evolving), accepted);
}

#[test]
fn g4_kinetic_changes_have_zero_instantaneous_work_and_refusals_are_atomic() {
    let gate = CancelGate::new_clock_free();
    let original = body(point(1.0, &[(10.0, 0.2)], 0.1), 1.0, DT, &gate);
    let different_tau = body(point(1.0, &[(10.0, 0.4)], 0.1), 1.0, DT, &gate);
    let different_kelvin = body(point(1.0, &[(10.0, 0.4)], 0.7), 1.0, DT, &gate);
    let mut history = original
        .history(&[0.02], &[-0.03], InitialMemory::Instantaneous, &gate)
        .unwrap();
    history.advance(&[0.01], &gate).unwrap();
    let tau_change = history
        .rebind_materials(&different_tau, allowance(0.0), &gate)
        .unwrap();
    assert_eq!(tau_change.parameter_work_j, 0.0);
    same_instant(&history, &tau_change.history);
    assert_eq!(history.state(), tau_change.history.state());
    let damping_change = tau_change
        .history
        .rebind_materials(&different_kelvin, allowance(0.0), &gate)
        .unwrap();
    assert_eq!(damping_change.parameter_work_j, 0.0);
    same_instant(&tau_change.history, &damping_change.history);
    assert_eq!(tau_change.history.state(), damping_change.history.state());
    let retained = damping_change.history;
    let accepted = snapshot(&retained);
    let next_point = point(1.0, &[(10.0, 0.4)], 0.7);
    let changed_mass = body(next_point.clone(), 2.0, DT, &gate);
    let changed_map = body(project(next_point.material(), 2.0), 1.0, DT, &gate);
    let changed_compatible_map = body(
        next_point
            .with_compatible_strain_columns(&[[2.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
            .unwrap(),
        1.0,
        DT,
        &gate,
    );
    let changed_clock = body(next_point, 1.0, DT / 2.0, &gate);
    for next in [
        &changed_mass,
        &changed_map,
        &changed_compatible_map,
        &changed_clock,
    ] {
        assert!(
            retained
                .rebind_materials(next, allowance(1.0), &gate)
                .is_err()
        );
        assert_eq!(snapshot(&retained), accepted);
    }
    let guarded = changed_compatible_map
        .history(&[0.01], &[0.0], InitialMemory::Relaxed, &gate)
        .unwrap();
    assert!(
        guarded
            .rebind_materials(&different_kelvin, allowance(1.0), &gate)
            .is_err()
    );
    let stiffer = body(point(2.0, &[(20.0, 0.4)], 0.7), 1.0, DT, &gate);
    assert!(
        retained
            .rebind_materials(&stiffer, allowance(0.0), &gate)
            .is_err()
    );
    assert_eq!(snapshot(&retained), accepted);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        retained
            .rebind_materials(&different_kelvin, allowance(0.0), &cancelled)
            .is_err()
    );
    assert_eq!(snapshot(&retained), accepted);
}
