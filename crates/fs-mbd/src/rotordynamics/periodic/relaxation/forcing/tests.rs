//! Independent original-coordinate held-load and hereditary-orbit oracles.
use super::*;
use crate::rotordynamics::RotorBearingModel;
use crate::rotordynamics::periodic::prepare_maxwell_cycle;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}");
}

fn arm(strain: &[f64], stiffness: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch { strain: strain.to_vec(), stiffness, relaxation_time_s: tau }
}

fn oscillator(stiffness: f64, damping: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[1.0], &[damping], &[stiffness], &[0.0], || false).unwrap()
}

// Independent exact physical diagonalization of M=1, C=2, K_relaxed=1.5,
// k=1.5, tau=1/4, b=1. The characteristic polynomial is
// (s+1)(s+2)(s+3), with e_r=[1,r,-(r²+2r+1.5)].
// A unit mechanical force has modal coefficients beta=[1.5,-2,0.5].
// No production exponential, factorization or eigenvalue is used here.
const RATES: [f64; 3] = [-1.0, -2.0, -3.0];
const BETA: [f64; 3] = [1.5, -2.0, 0.5];

fn exact_step(before: [f64; 3], force: f64, h: f64) -> [f64; 3] {
    std::array::from_fn(|i| {
        let rate = RATES[i];
        (rate * h).exp() * before[i] + BETA[i] * force * (rate * h).exp_m1() / rate
    })
}

fn physical(modes: [f64; 3]) -> RotorMaxwellState {
    let mut x = [0.0; 3];
    for i in 0..3 {
        let rate = RATES[i];
        x[0] += modes[i];
        x[1] += rate * modes[i];
        x[2] -= (rate * rate + 2.0 * rate + 1.5) * modes[i];
    }
    RotorMaxwellState {
        displacement: vec![x[0]], velocity: vec![x[1]], branch_stress: vec![x[2]],
    }
}

fn state_close(actual: &RotorMaxwellState, expected: &RotorMaxwellState, tolerance: f64) {
    assert_eq!(actual.displacement.len(), expected.displacement.len());
    assert_eq!(actual.velocity.len(), expected.velocity.len());
    assert_eq!(actual.branch_stress.len(), expected.branch_stress.len());
    for (&a, &b) in actual.displacement.iter().chain(&actual.velocity).chain(&actual.branch_stress)
        .zip(expected.displacement.iter().chain(&expected.velocity).chain(&expected.branch_stress)) {
        close(a, b, tolerance);
    }
}

#[test]
fn held_force_induces_original_memory_and_relaxes_to_original_static_stiffness() {
    let model = oscillator(1.5, 2.0);
    let branches = [arm(&[1.0], 1.5, 0.25)];
    let segments = [
        RotorPeriodicSegment { model: &model, duration_seconds: 0.2 },
        RotorPeriodicSegment { model: &model, duration_seconds: 0.5 },
    ];
    let force = [3.0];
    let prepared = prepare_maxwell_forced_cycle(7.0, &segments, &branches,
        &[&force, &force], 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    let zero = physical([0.0; 3]);
    let trajectory = prepared.advance(&zero, 1e-10, || false).unwrap();
    for (&time, state) in trajectory.boundary_times_seconds.iter().zip(&trajectory.boundary_states) {
        state_close(state, &physical(exact_step([0.0; 3], 3.0, time)), 3e-12);
    }
    assert!(prepared.forcing_offset()[2].abs() > 0.01);
    for (&a, &b) in prepared.forcing_offset().iter().zip(
        state_values(trajectory.boundary_states.last().unwrap(), 1, 1).unwrap().iter()) {
        close(a, b, 3e-12);
    }
    let equilibrium = RotorMaxwellState {
        displacement: vec![2.0], velocity: vec![0.0], branch_stress: vec![0.0],
    };
    let orbit = prepared.periodic_orbit(1e-10, 1e10, || false).unwrap();
    for state in &orbit.trajectory.boundary_states {
        state_close(state, &equilibrium, 3e-11);
    }
    assert!(orbit.relative_fixed_point_residual < 1e-10);
    assert!(orbit.relative_closure_residual < 1e-10);
    // A static generalized load sees K_relaxed, not K_relaxed + k*b*b^T.
    close(orbit.trajectory.boundary_states[0].displacement[0], 3.0 / 1.5, 3e-11);
    assert_eq!(orbit, prepared.periodic_orbit(1e-10, 1e10, || false).unwrap());
    let homogeneous = prepare_maxwell_cycle(7.0, &segments, &branches, 1e-12,
        RotorFloquetBudget::default(), || false).unwrap();
    assert_eq!(prepared.floquet_report().memory_states, 1);
    assert_eq!(prepared.floquet_report().multipliers.len(), 3);
    for (&a, &b) in prepared.floquet_report().monodromy.iter()
        .zip(&homogeneous.floquet_report().monodromy) {
        close(a, b, 3e-12);
    }
}

#[test]
fn asymmetric_load_cycle_matches_exact_modal_periodic_orbit_and_keeps_switch_memory() {
    let model = oscillator(1.5, 2.0);
    let branches = [arm(&[1.0], 1.5, 0.25)];
    let (a, b) = (0.3_f64, 0.4_f64);
    let (f1, f2) = ([0.7], [-0.2]);
    let segments = [
        RotorPeriodicSegment { model: &model, duration_seconds: a },
        RotorPeriodicSegment { model: &model, duration_seconds: b },
    ];
    let initial_modes = std::array::from_fn(|i| {
        let rate = RATES[i];
        let first = BETA[i] * f1[0] * (rate * a).exp_m1() / rate;
        let second = BETA[i] * f2[0] * (rate * b).exp_m1() / rate;
        ((rate * b).exp() * first + second) / -(rate * (a + b)).exp_m1()
    });
    let first_modes = exact_step(initial_modes, f1[0], a);
    let final_modes = exact_step(first_modes, f2[0], b);
    let expected = [physical(initial_modes), physical(first_modes), physical(final_modes)];
    assert!(expected[0].branch_stress[0].abs() > 1e-4);
    assert!(expected[1].branch_stress[0].abs() > 1e-4);
    let prepared = prepare_maxwell_forced_cycle(0.0, &segments, &branches,
        &[&f1, &f2], 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    let orbit = prepared.periodic_orbit(1e-10, 1e10, || false).unwrap();
    assert_eq!(prepared.boundary_times_seconds(), &[0.0, a, a + b]);
    for (state, target) in orbit.trajectory.boundary_states.iter().zip(&expected) {
        state_close(state, target, 5e-11);
    }
    let zero_modes = exact_step(exact_step([0.0; 3], f1[0], a), f2[0], b);
    let exact_offset = state_values(&physical(zero_modes), 1, 1).unwrap();
    for (&actual, &target) in prepared.forcing_offset().iter().zip(&exact_offset) {
        close(actual, target, 5e-12);
    }
    // Splitting a held interval changes no physical input or constitutive
    // history. A memory reset at either artificial boundary fails this oracle.
    let split = [
        RotorPeriodicSegment { model: &model, duration_seconds: 0.1 },
        RotorPeriodicSegment { model: &model, duration_seconds: 0.2 },
        segments[1],
    ];
    let refined = prepare_maxwell_forced_cycle(0.0, &split, &branches,
        &[&f1, &f1, &f2], 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    let refined_orbit = refined.periodic_orbit(1e-10, 1e10, || false).unwrap();
    state_close(&refined_orbit.trajectory.boundary_states[0], &expected[0], 5e-11);
    state_close(&refined_orbit.trajectory.boundary_states[1],
        &physical(exact_step(initial_modes, f1[0], 0.1)), 5e-11);
    state_close(&refined_orbit.trajectory.boundary_states[2], &expected[1], 5e-11);
    state_close(&refined_orbit.trajectory.boundary_states[3], &expected[2], 5e-11);
}

#[test]
fn full_mass_gyro_and_conjugate_force_units_preserve_the_complete_affine_orbit() {
    let mass = [2.0, 0.4, 0.4, 1.5];
    let damping = [0.12, -0.02, 0.07, 0.18];
    let gyro = [0.0, 0.8, -0.8, 0.0];
    let stiffness = [[4.0, 0.5, 0.5, 3.0], [8.0, -1.0, -1.0, 6.0]];
    let target = [0.7, -0.2];
    let forces: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|row| stiffness[i][2 * row] * target[0]
            + stiffness[i][2 * row + 1] * target[1])
    });
    let branches = [
        arm(&[1.2, -0.4], 3.5, 0.7), arm(&[-0.3, 0.8], 2.0, 0.2),
        arm(&[0.0, 0.0], 1.0, 0.4),
    ];
    let models: [RotorBearingModel; 2] = std::array::from_fn(|i| {
        RotorBearingModel::new(2, &mass, &damping, &stiffness[i], &gyro, || false).unwrap()
    });
    let segments = [
        RotorPeriodicSegment { model: &models[0], duration_seconds: 0.17 },
        RotorPeriodicSegment { model: &models[1], duration_seconds: 0.29 },
    ];
    let units = [1e3, 1e-2];
    let congruence = |matrix: &[f64; 4]| -> [f64; 4] {
        std::array::from_fn(|i| matrix[i] * units[i / 2] * units[i % 2])
    };
    let transformed: [RotorBearingModel; 2] = std::array::from_fn(|i| {
        RotorBearingModel::new(2, &congruence(&mass), &congruence(&damping),
            &congruence(&stiffness[i]), &congruence(&gyro), || false).unwrap()
    });
    let new_segments = [
        RotorPeriodicSegment { model: &transformed[0], duration_seconds: 0.17 },
        RotorPeriodicSegment { model: &transformed[1], duration_seconds: 0.29 },
    ];
    let new_branches: Vec<_> = branches.iter().map(|branch| arm(
        &[units[0] * branch.strain[0], units[1] * branch.strain[1]],
        branch.stiffness, branch.relaxation_time_s,
    )).collect();
    let new_forces: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| units[j] * forces[i][j])
    });
    let scale = [units[0], units[1], units[0], units[1], 1.0, 1.0, 1.0];
    for spin in [-2.7, 2.7] {
        let original = prepare_maxwell_forced_cycle(spin, &segments, &branches,
            &[&forces[0], &forces[1]], 1e-11, RotorFloquetBudget::default(), || false).unwrap();
        let transformed = prepare_maxwell_forced_cycle(spin, &new_segments, &new_branches,
            &[&new_forces[0], &new_forces[1]], 1e-11, RotorFloquetBudget::default(), || false).unwrap();
        for i in 0..7 {
            close(transformed.forcing_offset()[i], original.forcing_offset()[i] / scale[i], 2e-10);
        }
        let orbit = transformed.periodic_orbit(2e-9, 1e10, || false).unwrap();
        assert_eq!(transformed.floquet_report().multipliers.len(), 7);
        for state in &orbit.trajectory.boundary_states {
            close(state.displacement[0], target[0] / units[0], 2e-9);
            close(state.displacement[1], target[1] / units[1], 2e-9);
            for &value in state.velocity.iter().chain(&state.branch_stress) {
                close(value, 0.0, 2e-9);
            }
        }
    }
}

#[test]
fn unique_unstable_orbit_survives_but_free_and_unresolved_material_directions_refuse() {
    let model = oscillator(-1.0, 0.1);
    let branches = [arm(&[1.0], 0.2, 0.5)];
    let segments = [RotorPeriodicSegment { model: &model, duration_seconds: 0.3 }];
    let force = [2.0];
    let prepared = prepare_maxwell_forced_cycle(0.0, &segments, &branches,
        &[&force], 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    assert!(prepared.floquet_report().spectral_radius > 1.0);
    let orbit = prepared.periodic_orbit(1e-10, 1e10, || false).unwrap();
    for state in &orbit.trajectory.boundary_states {
        close(state.displacement[0], -2.0, 3e-11);
        close(state.velocity[0], 0.0, 3e-11);
        close(state.branch_stress[0], 0.0, 3e-11);
    }
    let free = oscillator(0.0, 2.0);
    let free_segments = [RotorPeriodicSegment { model: &free, duration_seconds: 0.3 }];
    for force in [[0.0], [1.0]] {
        let plan = prepare_maxwell_forced_cycle(0.0, &free_segments, &branches,
            &[&force], 1e-12, RotorFloquetBudget::default(), || false).unwrap();
        assert!(matches!(plan.periodic_orbit(1e-10, 1e10, || false),
            Err(RotorError::PeriodicLinear(_))));
    }
    // Even an unforced and mechanically dark original material direction
    // contributes to the complete periodic resolvent; it cannot be removed.
    let slow = [arm(&[0.0], 1.0, 1e8)];
    let stable = oscillator(4.0, 0.2);
    let slow_segments = [RotorPeriodicSegment { model: &stable, duration_seconds: 0.1 }];
    let plan = prepare_maxwell_forced_cycle(0.0, &slow_segments, &slow,
        &[&force], 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    assert!(matches!(plan.periodic_orbit(1e-10, 1e6, || false),
        Err(RotorError::PeriodicOrbitUnresolved { .. })));
}

#[test]
fn forcing_admission_cancellation_and_dense_budgets_return_no_partial_orbit() {
    let model = oscillator(1.5, 2.0);
    let branches = [arm(&[1.0], 1.5, 0.25)];
    let segments = [RotorPeriodicSegment { model: &model, duration_seconds: 0.2 }];
    let force = [1.0];
    let budget = RotorFloquetBudget::default();
    for rows in [vec![], vec![&[][..]], vec![&[f64::NAN][..]], vec![&force[..], &force[..]]] {
        assert!(prepare_maxwell_forced_cycle(0.0, &segments, &branches, &rows,
            1e-12, budget, || false).is_err());
    }
    assert!(matches!(prepare_maxwell_forced_cycle(0.0, &segments, &branches, &[&force],
        1e-12, budget, || true), Err(RotorError::Cancelled)));
    let mut polls = 0;
    assert!(matches!(prepare_maxwell_forced_cycle(0.0, &segments, &branches, &[&force],
        1e-12, budget, || { polls += 1; polls >= 20 }), Err(RotorError::Cancelled)));
    let tiny = RotorFloquetBudget { max_dense_bytes: 1, ..budget };
    assert!(matches!(prepare_maxwell_forced_cycle(0.0, &segments, &branches, &[&force],
        1e-12, tiny, || false), Err(RotorError::PeriodicBudget(_))));
    let plan = prepare_maxwell_forced_cycle(0.0, &segments, &branches, &[&force],
        1e-12, budget, || false).unwrap();
    assert!(matches!(plan.periodic_orbit(1e-10, 1e10, || true), Err(RotorError::Cancelled)));
    assert!(plan.periodic_orbit(0.0, 1e10, || false).is_err());
    assert!(plan.periodic_orbit(1e-10, f64::NAN, || false).is_err());
    assert!(plan.periodic_orbit(1e-10, 0.5, || false).is_err());
}
