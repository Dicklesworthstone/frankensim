//! Full consistent inertia, original Kelvin/Maxwell tensors, continuous forcing.
//! Independent original-coordinate RK4 oracle; all coefficients are authored fixtures.
use fs_couple::tensor_material::{
    InitialMemory,
    prepared::{PreparedTensorBody, PreparedTensorControls, PreparedTensorHistory},
};
use fs_exec::CancelGate;
use fs_material::visco::tensor::{PronyBranch, TensorProny};
use fs_material::visco_tensor::{ProjectedTensorProny, TensorProjectionBudget};

const MASS: [f64; 4] = [2.0, 0.35, 0.35, 1.4];
const B: [[f64; 2]; 2] = [[1.0, -0.4], [0.3, 1.0]];
const VOLUME: [f64; 2] = [0.5, 0.7];
const ELASTIC: [f64; 2] = [3.0, 4.0];
const KELVIN: [f64; 2] = [0.2, 0.07];
const BRANCH: [f64; 2] = [10.0, 3.0];
const TAU: [f64; 2] = [0.2, 0.7];
const OMEGA: [f64; 2] = [0.7, 3.0];
const COSINE: [[f64; 2]; 2] = [[0.12, -0.09], [-0.07, 0.11]];
const SINE: [[f64; 2]; 2] = [[0.03, 0.06], [0.09, -0.04]];

fn controls() -> PreparedTensorControls {
    PreparedTensorControls {
        maximum_coordinates: 2,
        maximum_material_points: 2,
        maximum_state_coordinates: 16,
        maximum_heat_groups: 2,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_preparation_multiply_terms: 10_000_000_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-10,
        maximum_work_error_j: 1e-10,
        maximum_steps: 100,
    }
}
fn points() -> Vec<ProjectedTensorProny> {
    let diagonal = |value| {
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0.0 }))
    };
    (0..2)
        .map(|p| {
            TensorProny::new(
                diagonal(ELASTIC[p]),
                Some(diagonal(KELVIN[p])),
                vec![PronyBranch::new(diagonal(BRANCH[p]), TAU[p]).unwrap()],
                (0.0, 10.0),
            )
            .unwrap()
            .project_strain(
                &B[p].map(|b| [b, 0.0, 0.0, 0.0, 0.0, 0.0]),
                VOLUME[p],
                TensorProjectionBudget {
                    max_coordinates: 2,
                    max_scalar_terms: 10_000,
                },
                || false,
            )
            .unwrap()
        })
        .collect()
}
fn body(dt: f64, frequencies: &[f64], scale: f64, gate: &CancelGate) -> PreparedTensorBody {
    PreparedTensorBody::from_mass_matrix_tones(
        2,
        MASS.to_vec(),
        points(),
        vec![0, 1],
        dt,
        frequencies,
        scale,
        controls(),
        gate,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:e} != {b:e}, allowance={tolerance:e}"
    );
}
fn snapshot(history: &PreparedTensorHistory<'_>) -> Vec<u64> {
    history
        .state()
        .iter()
        .chain(history.velocity_mps())
        .copied()
        .chain([
            history.time_s(),
            history.supplied_work_j(),
            history.dissipated_heat_j(),
        ])
        .map(f64::to_bits)
        .collect()
}

// q[2], v[2], spring-displacement histories[4], heat[2], work,
// then cosine/sine velocity integrals in each tone's local interval phase.
fn derivative(x: &[f64; 19], t: f64, t0: f64) -> [f64; 19] {
    let mut out = [0.0; 19];
    let mut force = [0.0; 2];
    for tone in 0..2 {
        let (s, c) = (OMEGA[tone] * t).sin_cos();
        let (sl, cl) = (OMEGA[tone] * (t - t0)).sin_cos();
        for i in 0..2 {
            force[i] += COSINE[tone][i] * c + SINE[tone][i] * s;
            out[11 + 4 * tone + i] = cl * x[2 + i];
            out[11 + 4 * tone + 2 + i] = sl * x[2 + i];
        }
    }
    out[10] = force[0] * x[2] + force[1] * x[3];
    for p in 0..2 {
        let strain = B[p][0] * x[0] + B[p][1] * x[1];
        let rate = B[p][0] * x[2] + B[p][1] * x[3];
        let spring = B[p][0] * x[4 + 2 * p] + B[p][1] * x[5 + 2 * p];
        let effort = VOLUME[p] * (ELASTIC[p] * strain + KELVIN[p] * rate + BRANCH[p] * spring);
        for i in 0..2 {
            force[i] -= B[p][i] * effort;
            out[4 + 2 * p + i] = x[2 + i] - x[4 + 2 * p + i] / TAU[p];
        }
        out[8 + p] = VOLUME[p] * (KELVIN[p] * rate * rate + BRANCH[p] * spring * spring / TAU[p]);
    }
    out[..2].copy_from_slice(&x[2..4]);
    let determinant = MASS[0] * MASS[3] - MASS[1] * MASS[2];
    out[2] = (MASS[3] * force[0] - MASS[1] * force[1]) / determinant;
    out[3] = (MASS[0] * force[1] - MASS[2] * force[0]) / determinant;
    out
}
fn reference_step(x: &mut [f64; 19], t0: f64, dt: f64) {
    x[11..].fill(0.0);
    let h = dt / 512.0;
    for j in 0..512 {
        let t = t0 + j as f64 * h;
        let a = derivative(x, t, t0);
        let b = derivative(
            &core::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
            t + 0.5 * h,
            t0,
        );
        let c = derivative(
            &core::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
            t + 0.5 * h,
            t0,
        );
        let d = derivative(&core::array::from_fn(|i| x[i] + h * c[i]), t + h, t0);
        for i in 0..19 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
}

#[test]
fn g1_multitone_full_mass_memory_heat_and_work_match_independent_equations() {
    let gate = CancelGate::new_clock_free();
    let dt = 0.04;
    for scale in [0.1, 10.0] {
        let body = body(dt, &OMEGA, scale, &gate);
        assert_eq!(body.state_dimension(), 8);
        assert_eq!(body.propagation_state_dimension(), 16);
        let mut history = body
            .history(
                &[0.025, -0.03],
                &[-0.04, 0.06],
                InitialMemory::Instantaneous,
                &gate,
            )
            .unwrap();
        let mut x = [0.0; 19];
        x[..8].copy_from_slice(&[0.025, -0.03, -0.04, 0.06, 0.025, -0.03, 0.025, -0.03]);
        let mut cosine: Vec<_> = COSINE.iter().map(|v| v.to_vec()).collect();
        let mut sine: Vec<_> = SINE.iter().map(|v| v.to_vec()).collect();
        let mut largest_mean_only_error = 0.0_f64;
        for step in 0..20 {
            let old_heat = [x[8], x[9]];
            let old_work = x[10];
            reference_step(&mut x, step as f64 * dt, dt);
            let result = history.advance_tones(&cosine, &sine, &gate).unwrap();
            for i in 0..2 {
                close(history.displacement_m()[i], x[i], 3e-11);
                close(history.velocity_mps()[i], x[2 + i], 3e-11);
                close(
                    result.motion.heat_by_group_j[i],
                    x[8 + i] - old_heat[i],
                    3e-11,
                );
                assert!(result.motion.heat_by_group_j[i] >= 0.0);
            }
            for i in 0..4 {
                close(history.memory_displacement_m()[i], x[4 + i], 3e-11);
            }
            for tone in 0..2 {
                for i in 0..2 {
                    close(
                        result.cosine_velocity_moment_mps[tone][i],
                        x[11 + 4 * tone + i] / dt,
                        5e-10,
                    );
                    close(
                        result.sine_velocity_moment_mps[tone][i],
                        x[13 + 4 * tone + i] / dt,
                        5e-10,
                    );
                }
            }
            close(
                result.motion.energy.supplied_work_j,
                x[10] - old_work,
                3e-11,
            );
            close(history.dissipated_heat_j(), x[8] + x[9], 3e-11);
            close(history.supplied_work_j(), x[10], 3e-11);
            assert!(result.motion.energy.balance_error_j.abs() < 1e-10);
            let sampled = (0..2)
                .map(|i| {
                    dt * result.motion.mean_velocity_mps[i]
                        * (0..2).map(|tone| cosine[tone][i]).sum::<f64>()
                })
                .sum::<f64>();
            largest_mean_only_error =
                largest_mean_only_error.max((sampled - (x[10] - old_work)).abs());
            (cosine, sine) = body.advance_force_phase(&cosine, &sine, &gate).unwrap();
        }
        assert!(largest_mean_only_error > 1e-7);
    }
}

#[test]
fn g3_zero_frequency_reduces_to_the_original_held_owner_without_sine_work() {
    let gate = CancelGate::new_clock_free();
    let tones = body(0.04, &[0.0], 1.0, &gate);
    let held = PreparedTensorBody::from_mass_matrix(
        2,
        MASS.to_vec(),
        points(),
        vec![0, 1],
        0.04,
        controls(),
        &gate,
    )
    .unwrap();
    let mut a = tones
        .history(
            &[0.025, -0.03],
            &[-0.04, 0.06],
            InitialMemory::Relaxed,
            &gate,
        )
        .unwrap();
    let mut b = held
        .history(
            &[0.025, -0.03],
            &[-0.04, 0.06],
            InitialMemory::Relaxed,
            &gate,
        )
        .unwrap();
    let ar = a
        .advance_tones(&[vec![0.1, -0.2]], &[vec![0.7, -0.3]], &gate)
        .unwrap();
    let br = b.advance(&[0.1, -0.2], &gate).unwrap();
    for (x, y) in a.state().iter().zip(b.state()) {
        close(*x, *y, 3e-11);
    }
    close(
        ar.motion.energy.supplied_work_j,
        br.energy.supplied_work_j,
        3e-11,
    );
    close(ar.motion.energy.dissipated_j, br.energy.dissipated_j, 3e-11);
    for value in &ar.sine_velocity_moment_mps[0] {
        close(*value, 0.0, 1e-15);
    }
}

#[test]
fn g0_frequency_capacity_and_atomic_refusals_keep_every_history() {
    let gate = CancelGate::new_clock_free();
    for frequencies in [vec![], vec![-1.0], vec![f64::NAN], vec![11.0]] {
        assert!(
            PreparedTensorBody::from_mass_matrix_tones(
                2,
                MASS.to_vec(),
                points(),
                vec![0, 1],
                0.04,
                &frequencies,
                1.0,
                controls(),
                &gate,
            )
            .is_err()
        );
    }
    let mut short = controls();
    short.maximum_state_coordinates = 15;
    assert!(
        PreparedTensorBody::from_mass_matrix_tones(
            2,
            MASS.to_vec(),
            points(),
            vec![0, 1],
            0.04,
            &OMEGA,
            1.0,
            short,
            &gate,
        )
        .is_err()
    );
    let body = body(0.04, &OMEGA, 1.0, &gate);
    let mut history = body
        .history(
            &[0.025, -0.03],
            &[-0.04, 0.06],
            InitialMemory::Instantaneous,
            &gate,
        )
        .unwrap();
    let before = snapshot(&history);
    let cosine: Vec<_> = COSINE.iter().map(|v| v.to_vec()).collect();
    let sine: Vec<_> = SINE.iter().map(|v| v.to_vec()).collect();
    assert!(history.advance(&[0.1, 0.1], &gate).is_err());
    assert!(history.advance_tones(&cosine[..1], &sine, &gate).is_err());
    let mut invalid = cosine.clone();
    invalid[1][0] = f64::NAN;
    assert!(history.advance_tones(&invalid, &sine, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(history.advance_tones(&cosine, &sine, &cancelled).is_err());
    assert_eq!(snapshot(&history), before);
    let mut replay = history.clone();
    history.advance_tones(&cosine, &sine, &gate).unwrap();
    replay.advance_tones(&cosine, &sine, &gate).unwrap();
    assert_eq!(snapshot(&history), snapshot(&replay));
}
