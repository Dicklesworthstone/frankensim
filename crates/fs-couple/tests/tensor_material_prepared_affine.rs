//! Independent analytical affine-force motion, work and material heat.
//! Authored constitutive fixtures; no measured-material qualification claim.
use fs_couple::tensor_material::{
    InitialMemory,
    prepared::{PreparedTensorBody, PreparedTensorControls, PreparedTensorHistory},
};
use fs_exec::CancelGate;
use fs_material::visco::tensor::{PronyBranch, TensorProny};
use fs_material::visco_tensor::{ProjectedTensorProny, TensorProjectionBudget};

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

fn point(limit: f64) -> ProjectedTensorProny {
    let diagonal = |value| {
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0.0 }))
    };
    TensorProny::new_mandel(
        diagonal(1.0),
        vec![PronyBranch::new(diagonal(10.0), 1.0 / 6.0).unwrap()],
        limit,
    )
    .unwrap()
    .project_strain(
        &[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
        1.0,
        TensorProjectionBudget {
            max_coordinates: 1,
            max_scalar_terms: 10_000,
        },
        || false,
    )
    .unwrap()
}

fn body(
    dt: f64,
    force_scale: f64,
    controls: PreparedTensorControls,
    gate: &CancelGate,
) -> PreparedTensorBody {
    PreparedTensorBody::from_mass_matrix_affine(
        1,
        vec![1.0],
        vec![point(1.0)],
        vec![0],
        dt,
        force_scale,
        controls,
        gate,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}, allowance={tolerance:e}"
    );
}

/// m=1, K_eq=1, K_Maxwell=10, tau=1/6, F(t)=F0+b*t.
/// Eliminating the branch gives u'''+6u''+11u'+6u=6F+b. Its homogeneous
/// roots are -1,-2,-3 and its particular solution is F0+b*t-(5/3)b.
/// Work integrates F*v; heat integrates 0.6*branch_force^2 independently.
/// Returns [u, v, spring displacement, heat, work, integral(t*v dt)].
fn analytical(t: f64, u0: f64, v0: f64, f0: f64, b: f64, memory: InitialMemory) -> [f64; 6] {
    let branch0 = if memory == InitialMemory::Instantaneous {
        10.0 * u0
    } else {
        0.0
    };
    let a0 = f0 - u0 - branch0;
    let shifted_u = u0 - f0 + (5.0 / 3.0) * b;
    let shifted_v = v0 - b;
    let amplitudes = [
        3.0 * shifted_u + 2.5 * shifted_v + 0.5 * a0,
        -3.0 * shifted_u - 4.0 * shifted_v - a0,
        shifted_u + 1.5 * shifted_v + 0.5 * a0,
    ];
    let rates: [f64; 3] = [1.0, 2.0, 3.0];
    let branch_amplitudes: [f64; 3] =
        core::array::from_fn(|i| -(1.0 + rates[i] * rates[i]) * amplitudes[i]);
    let branch_constant = (5.0 / 3.0) * b;
    let mut u = f0 + b * t - branch_constant;
    let mut v = b;
    let mut branch = branch_constant;
    let mut time_velocity = 0.5 * b * t * t;
    let mut heat = 0.6 * branch_constant * branch_constant * t;
    for i in 0..3 {
        let decay = (-rates[i] * t).exp();
        let integral_decay = -(-rates[i] * t).exp_m1() / rates[i];
        u += amplitudes[i] * decay;
        v -= rates[i] * amplitudes[i] * decay;
        branch += branch_amplitudes[i] * decay;
        time_velocity -= amplitudes[i] * (integral_decay - t * decay);
        heat += 1.2 * branch_constant * branch_amplitudes[i] * integral_decay;
        for j in 0..3 {
            let rate = rates[i] + rates[j];
            heat +=
                0.6 * branch_amplitudes[i] * branch_amplitudes[j] * -(-rate * t).exp_m1() / rate;
        }
    }
    [
        u,
        v,
        branch / 10.0,
        heat,
        f0 * (u - u0) + b * time_velocity,
        time_velocity,
    ]
}

#[test]
fn g1_affine_zener_matches_analytical_motion_work_heat_and_velocity_moments() {
    let gate = CancelGate::new_clock_free();
    let dt = 0.025;
    let (u0, v0, f0, slope) = (0.025, -0.04, 0.012, 0.08);
    let mut largest_average_force_work_error = 0.0_f64;
    for force_scale in [0.01, 1.0, 100.0] {
        let body = body(dt, force_scale, controls(), &gate);
        assert_eq!(body.state_dimension(), 3);
        assert_eq!(body.propagation_state_dimension(), 4);
        for memory in [InitialMemory::Relaxed, InitialMemory::Instantaneous] {
            let mut history = body.history(&[u0], &[v0], memory, &gate).unwrap();
            let mut previous = analytical(0.0, u0, v0, f0, slope, memory);
            for step in 1..=40 {
                let t0 = (step - 1) as f64 * dt;
                let t1 = step as f64 * dt;
                let start = f0 + slope * t0;
                let end = f0 + slope * t1;
                let report = history.advance_affine(&[start], &[end], &gate).unwrap();
                let expected = analytical(t1, u0, v0, f0, slope, memory);
                let du = expected[0] - previous[0];
                let work = expected[4] - previous[4];
                let first_velocity = (expected[5] - previous[5] - t0 * du) / (dt * dt);
                close(history.displacement_m()[0], expected[0], 2e-12);
                close(history.velocity_mps()[0], expected[1], 2e-12);
                close(history.memory_displacement_m()[0], expected[2], 2e-12);
                close(history.dissipated_heat_j(), expected[3], 3e-12);
                close(history.supplied_work_j(), expected[4], 3e-12);
                close(report.motion.energy.supplied_work_j, work, 3e-12);
                close(
                    report.motion.heat_by_group_j[0],
                    expected[3] - previous[3],
                    3e-12,
                );
                close(report.motion.mean_velocity_mps[0], du / dt, 2e-11);
                close(report.first_velocity_moment_mps[0], first_velocity, 2e-10);
                assert!(report.motion.heat_by_group_j[0] >= 0.0);
                assert_eq!(history.state().len(), 3);
                largest_average_force_work_error =
                    largest_average_force_work_error.max((work - 0.5 * (start + end) * du).abs());
                previous = expected;
            }
        }
    }
    // A mean force paired only with mean velocity cannot reproduce this work.
    assert!(largest_average_force_work_error > 1e-9);
    println!(
        "affine Zener: two prehistories, three force scales, average-force work error={largest_average_force_work_error:e} J"
    );
}

#[test]
fn g3_affine_equal_endpoints_reproduce_original_held_material_dynamics() {
    let gate = CancelGate::new_clock_free();
    let dt = 0.05;
    let affine = body(dt, 0.1, controls(), &gate);
    let held = PreparedTensorBody::from_mass_matrix(
        1,
        vec![1.0],
        vec![point(1.0)],
        vec![0],
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    let mut a = affine
        .history(&[0.025], &[-0.04], InitialMemory::Instantaneous, &gate)
        .unwrap();
    let mut convenience = a.clone();
    let mut h = held
        .history(&[0.025], &[-0.04], InitialMemory::Instantaneous, &gate)
        .unwrap();
    for step in 0..20 {
        let force = [if step < 10 { 0.012 } else { -0.009 }];
        let ar = a.advance_affine(&force, &force, &gate).unwrap();
        let cr = convenience.advance(&force, &gate).unwrap();
        let hr = h.advance(&force, &gate).unwrap();
        for (x, y) in a.state().iter().zip(h.state()) {
            close(*x, *y, 2e-12);
        }
        close(
            ar.motion.energy.supplied_work_j,
            hr.energy.supplied_work_j,
            2e-12,
        );
        close(ar.motion.energy.dissipated_j, hr.energy.dissipated_j, 2e-12);
        assert_eq!(a.state(), convenience.state());
        assert_eq!(ar.motion.energy, cr.energy);
    }
}

fn snapshot(h: &PreparedTensorHistory<'_>) -> Vec<u64> {
    let mut result: Vec<_> = h
        .state()
        .iter()
        .chain(h.velocity_mps())
        .map(|x| x.to_bits())
        .collect();
    result.extend(
        [
            h.time_s(),
            h.stored_energy_j(),
            h.supplied_work_j(),
            h.dissipated_heat_j(),
        ]
        .map(f64::to_bits),
    );
    result.push(h.steps());
    result
}

#[test]
fn g4_affine_source_caps_and_failed_intervals_preserve_physical_history() {
    let gate = CancelGate::new_clock_free();
    let mut cap = controls();
    cap.maximum_state_coordinates = 3;
    assert!(
        PreparedTensorBody::from_mass_matrix_affine(
            1,
            vec![1.0],
            vec![point(1.0)],
            vec![0],
            0.025,
            1.0,
            cap,
            &gate,
        )
        .is_err()
    );
    let mut budget = controls();
    budget.maximum_steps = 1;
    let body = body(0.025, 1.0, budget, &gate);
    let mut history = body
        .history(&[0.025], &[-0.04], InitialMemory::Instantaneous, &gate)
        .unwrap();
    let before = snapshot(&history);
    for (start, end) in [
        (vec![0.01], vec![]),
        (vec![0.01], vec![f64::NAN]),
        (vec![f64::MAX], vec![-f64::MAX]),
    ] {
        assert!(history.advance_affine(&start, &end, &gate).is_err());
        assert_eq!(snapshot(&history), before);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        history
            .advance_affine(&[0.01], &[0.02], &cancelled)
            .is_err()
    );
    assert_eq!(snapshot(&history), before);
    history.advance_affine(&[0.01], &[0.02], &gate).unwrap();
    let accepted = snapshot(&history);
    assert!(history.advance_affine(&[0.02], &[0.03], &gate).is_err());
    assert_eq!(snapshot(&history), accepted);
    let held = PreparedTensorBody::from_mass_matrix(
        1,
        vec![1.0],
        vec![point(1.0)],
        vec![0],
        0.025,
        controls(),
        &gate,
    )
    .unwrap();
    let mut held_history = held
        .history(&[0.025], &[-0.04], InitialMemory::Instantaneous, &gate)
        .unwrap();
    let held_before = snapshot(&held_history);
    assert!(
        held_history
            .advance_affine(&[0.01], &[0.02], &gate)
            .is_err()
    );
    assert_eq!(snapshot(&held_history), held_before);
}
