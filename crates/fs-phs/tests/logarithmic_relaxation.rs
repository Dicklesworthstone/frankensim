//! G0/G1/G3: nonpolynomial original strain, memory and coherent physical work.
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget, PreparedNonlinearTones};
use fs_phs::mechanical::LinearMechanicalSystem;
use fs_phs::{
    LogarithmicStrainMap, MappedRelaxationBranch, PortHamiltonian, QuadraticStorage,
    RelaxationStrainDiscreteGradient, RelaxationStrainMap, ResistanceWitness,
};
use std::rc::Rc;

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tol * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 64,
        max_omega_dt: 0.04,
    }
}
fn dense() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}
const M: f64 = 1.3;
const K: f64 = 1.7;
const C: f64 = 0.15;
const B: f64 = 0.65;
const ARM: f64 = 0.8;
const TAU: f64 = 0.3;
fn system() -> PortHamiltonian {
    LinearMechanicalSystem::new(&[M], &[C], &[K], &[1.0], 1, 1, vec![])
        .unwrap()
        .into_system()
        .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
            strain: Rc::new(LogarithmicStrainMap::new(1.0, vec![B, 0.0]).unwrap()),
            stiffness: ARM,
            relaxation_time_s: TAU,
        }])
        .unwrap()
}

#[test]
fn original_logarithmic_component_work_and_both_complete_tangents() {
    let sys = system();
    let (a, b) = ([0.2, 4.0, 0.1], [0.8, 9.0, -0.2]);
    let mid = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(a[i], b[i]));
    let rule = sys.storage().discrete_gradient_rule().unwrap();
    let mut effort = [0.0; 3];
    rule.evaluate(&a, &b, &mid, &mut effort);
    let stretch = |q: f64| 1.0 + B * q;
    let strain = |q: f64| stretch(q).ln();
    let stress = |x: [f64; 3]| ARM.sqrt() * strain(x[0]) - x[2];
    let mean = 0.5 * (stress(a) + stress(b));
    let dg = (strain(b[0]) - strain(a[0])) / (b[0] - a[0]);
    close(effort[0], K * mid[0] + mean * ARM.sqrt() * dg, 3e-15);
    close(effort[1], mid[1] / M, 3e-15);
    close(effort[2], -mean, 3e-15);
    let energy =
        |x: [f64; 3]| 0.5 * K * x[0] * x[0] + 0.5 * x[1] * x[1] / M + 0.5 * stress(x).powi(2);
    close(
        (0..3).map(|i| effort[i] * (b[i] - a[i])).sum(),
        energy(b) - energy(a),
        3e-15,
    );
    let mut swapped = [0.0; 3];
    rule.evaluate(&b, &a, &mid, &mut swapped);
    for i in 0..3 {
        close(swapped[i], effort[i], 3e-15);
    }

    let v = [0.3, -0.2, 0.4];
    let mut action = [0.0; 3];
    assert!(rule.end_jacobian_action(&a, &b, &v, &mut action));
    let ds = 0.5 * (ARM.sqrt() * B / stretch(b[0]) * v[0] - v[2]);
    let ddg = (B / stretch(b[0]) - dg) / (b[0] - a[0]) * v[0];
    close(
        action[0],
        0.5 * K * v[0] + ARM.sqrt() * (ds * dg + mean * ddg),
        2e-14,
    );
    close(action[1], 0.5 * v[1] / M, 3e-15);
    close(action[2], -ds, 3e-15);
    let mut hessian = [0.0; 3];
    assert!(
        sys.storage()
            .continuous_tangent()
            .unwrap()
            .apply(&a, &v, &mut hessian)
    );
    assert!(rule.end_jacobian_action(&a, &a, &v, &mut action));
    for i in 0..3 {
        close(action[i], 0.5 * hessian[i], 3e-15);
    }

    // Geometry is nonlinear even when Q propagates part of the model exactly.
    let part = sys.storage().quadratic_part().unwrap();
    assert!(!part.exact);
    assert_eq!(part.matrix[2 * 3 + 2], 1.0);
    assert_eq!(part.matrix[2], 0.0);
    assert_eq!(part.matrix[6], 0.0);
    let high_a = [0.2, 1e12, 1e8];
    let high_b = [0.8, 1e12 + 1.0, 1e8];
    let high_mid = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(high_a[i], high_b[i]));
    rule.evaluate(&high_a, &high_b, &high_mid, &mut effort);
    close(effort[1], high_mid[1] / M, 3e-15);
    // The independent strain remainder still resolves cubic geometric work
    // beside large prestress and independent kinetic energy.
    let rem = sys
        .storage()
        .midpoint_energy_remainder(&high_a, &high_b, &high_mid)
        .unwrap();
    let ea = stress(high_a);
    let eb = stress(high_b);
    let em = stress(high_mid);
    let change = ARM.sqrt() * (strain(high_b[0]) - strain(high_a[0]));
    let expected = 0.5 * (ea + eb) * change
        - em * ARM.sqrt() * B / stretch(high_mid[0]) * (high_b[0] - high_a[0]);
    close(rem, expected, 2e-12);
}

#[test]
fn logarithmic_divided_differences_keep_near_coincidence_and_geometric_units() {
    for scale in [1e-100, 1.0, 1e100] {
        let map = LogarithmicStrainMap::new(scale, vec![scale, 0.0]).unwrap();
        for delta in [0.0, 1e-12, 0.1, 0.8, -0.6] {
            let a = [0.2, 4.0];
            let b = [0.2 + delta, 9.0];
            let mid = [f64::midpoint(a[0], b[0]), 6.5];
            let mut out = [f64::NAN; 2];
            assert!(RelaxationStrainDiscreteGradient::evaluate(
                &map, &a, &b, &mid, &mut out
            ));
            let reference = if delta.abs() < 1e-8 {
                1.0 / (1.0 + mid[0])
            } else {
                ((1.0 + b[0]) / (1.0 + a[0])).ln() / (b[0] - a[0])
            };
            close(out[0], reference, 2e-14);
            assert_eq!(out[1], 0.0);
            let mut action = [0.0; 2];
            assert!(RelaxationStrainDiscreteGradient::end_jacobian_action(
                &map,
                &a,
                &b,
                &[1.0, 5.0],
                &mut action
            ));
            let expected = if delta.abs() < 1e-8 {
                -0.5 / (1.0 + mid[0]).powi(2)
            } else {
                (1.0 / (1.0 + b[0]) - reference) / (b[0] - a[0])
            };
            close(action[0], expected, 2e-12);
            assert_eq!(action[1], 0.0);
        }
    }
}

const FREQUENCIES: [f64; 3] = [0.0, 1.1, 2.7];
const C0: [f64; 3] = [0.1, 0.8, -0.25];
const S0: [f64; 3] = [0.0, -0.2, 0.3];
const C1: [f64; 3] = [-0.1, 0.35, 0.2];
const S1: [f64; 3] = [0.0, 0.3, -0.1];
const DURATION: f64 = 1.2;
fn rhs(t: f64, x: [f64; 6], affine: bool) -> [f64; 6] {
    let r = if affine { t / DURATION } else { 0.0 };
    let u: f64 = (0..3)
        .map(|j| {
            let (s, c) = (FREQUENCIES[j] * t).sin_cos();
            (C0[j] + r * (C1[j] - C0[j])) * c + (S0[j] + r * (S1[j] - S0[j])) * s
        })
        .sum();
    let v = x[1] / M;
    let stretch = 1.0 + B * x[0];
    let stress = ARM.sqrt() * stretch.ln() - x[2];
    [
        v,
        u - K * x[0] - C * v - ARM.sqrt() * B / stretch * stress,
        stress / TAU,
        u * v,
        C * v * v,
        stress * stress / TAU,
    ]
}
fn reference(affine: bool, steps: usize) -> [f64; 6] {
    let mut x = [0.3, -0.1, -0.05, 0.0, 0.0, 0.0];
    let h = DURATION / steps as f64;
    let add = |x: [f64; 6], v: [f64; 6], w: f64| core::array::from_fn(|i| x[i] + w * v[i]);
    for k in 0..steps {
        let t = k as f64 * h;
        let a = rhs(t, x, affine);
        let b = rhs(t + 0.5 * h, add(x, a, 0.5 * h), affine);
        let c = rhs(t + 0.5 * h, add(x, b, 0.5 * h), affine);
        let d = rhs(t + h, add(x, c, h), affine);
        for i in 0..6 {
            x[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    x
}
fn simulate(affine: bool, steps: usize) -> [f64; 6] {
    let h = DURATION / steps as f64;
    let channels = core::array::from_fn::<_, 2, _>(|i| {
        let mut matrix = vec![0.0; 9];
        let mut row = vec![0.0; 3];
        let coordinate = i + 1;
        let rate = if i == 0 { C } else { 1.0 / TAU };
        matrix[coordinate * 3 + coordinate] = rate;
        row[coordinate] = rate.sqrt();
        ResistanceWitness::try_from_source_rows(3, matrix, row).unwrap()
    });
    let constructor = if affine {
        PreparedNonlinearTones::new_owned_affine_with_resistance_channels
    } else {
        PreparedNonlinearTones::new_owned_with_resistance_channels
    };
    let mut prepared = constructor(system(), &FREQUENCIES, h, split(), dense(), &channels).unwrap();
    let mut x = [0.3, -0.1, -0.05];
    let (mut next, mut mean, mut heats) = ([0.0; 3], [0.0; 1], [0.0; 2]);
    let mut totals = [0.0; 3];
    for step in 0..steps {
        let t = step as f64 * h;
        let mut knots = [[0.0; 3]; 4];
        for j in 0..3 {
            let (sn, cs) = (FREQUENCIES[j] * t).sin_cos();
            for (offset, time) in [(0, t), (2, t + h)] {
                let r = if affine { time / DURATION } else { 0.0 };
                let c = C0[j] + r * (C1[j] - C0[j]);
                let s = S0[j] + r * (S1[j] - S0[j]);
                knots[offset][j] = c * cs + s * sn;
                knots[offset + 1][j] = -c * sn + s * cs;
            }
        }
        let record = if affine {
            prepared.step_affine_observed_into(
                &x, &knots[0], &knots[1], &knots[2], &knots[3], &mut next, &mut mean, &mut heats,
            )
        } else {
            prepared.step_observed_into(&x, &knots[0], &knots[1], &mut next, &mut mean, &mut heats)
        }
        .unwrap();
        close(record.balance_residual(), 0.0, 2e-12);
        close(record.dissipated, heats[0] + heats[1], 2e-12);
        assert!(heats.iter().all(|v| *v >= 0.0));
        totals[0] += record.supplied;
        totals[1] += heats[0];
        totals[2] += heats[1];
        x = next;
    }
    [x[0], x[1], x[2], totals[0], totals[1], totals[2]]
}
#[test]
fn nonlinear_logarithmic_memory_and_separate_heat_converge_under_coherent_affine_tones() {
    for affine in [false, true] {
        let coarse = reference(affine, 8192);
        let reference = reference(affine, 16384);
        for i in 0..6 {
            close(coarse[i], reference[i], 2e-12);
        }
        let mut errors = [0.0_f64; 3];
        for (rung, steps) in [24, 48, 96].into_iter().enumerate() {
            let actual = simulate(affine, steps);
            errors[rung] = (0..6)
                .map(|i| (actual[i] - reference[i]).abs())
                .fold(0.0, f64::max);
        }
        assert!(
            errors[0] > 1e-8 && errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1e-5, "{errors:?}");
    }
}

#[test]
fn exact_fixed_extension_memory_decay_and_domain_cancellation_are_original() {
    let fixed = PortHamiltonian::new(
        2,
        1,
        vec![0.0; 4],
        vec![0.0; 4],
        vec![0.0; 2],
        Box::new(QuadraticStorage::new(vec![0.0; 4], 2).unwrap()),
    )
    .unwrap()
    .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
        strain: Rc::new(LogarithmicStrainMap::new(1.0, vec![B, 0.0]).unwrap()),
        stiffness: ARM,
        relaxation_time_s: 0.02,
    }])
    .unwrap();
    let mut prepared =
        PreparedNonlinearTones::new_owned(fixed, &[2.1], 0.4, split(), dense()).unwrap();
    let initial = [0.3, 0.0, -0.1];
    let target = ARM.sqrt() * (1.0 + B * initial[0]).ln();
    let (mut next, mut mean) = ([7.0; 3], [9.0; 1]);
    assert!(
        prepared
            .step_into_controlled(&initial, &[0.7], &[-0.2], &mut next, &mut mean, || true)
            .is_err()
    );
    assert_eq!(next, [7.0; 3]);
    assert_eq!(mean, [9.0; 1]);
    assert!(
        prepared
            .step_into(&[-2.0, 0.0, 0.0], &[0.7], &[-0.2], &mut next, &mut mean)
            .is_err()
    );
    assert_eq!(next, [7.0; 3]);
    assert_eq!(mean, [9.0; 1]);
    let record = prepared
        .step_into(&initial, &[0.7], &[-0.2], &mut next, &mut mean)
        .unwrap();
    assert_eq!(next[0], initial[0]);
    assert_eq!(next[1], initial[1]);
    close(
        next[2],
        target + (initial[2] - target) * (-20.0_f64).exp(),
        3e-13,
    );
    close(
        record.dissipated,
        0.5 * (target - initial[2]).powi(2) * (1.0 - (-40.0_f64).exp()),
        3e-13,
    );
    assert_eq!(record.supplied, 0.0);
    close(record.balance_residual(), 0.0, 3e-13);
    assert!(LogarithmicStrainMap::new(0.0, vec![1.0]).is_err());
}
