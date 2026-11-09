use super::*;
use std::cell::Cell;

fn matrix(n: usize, m: usize, values: &[f64]) -> Csr {
    let mut a = Coo::new(n, m);
    for i in 0..n {
        for j in 0..m {
            if values[i * m + j] != 0.0 {
                a.push(i, j, values[i * m + j]);
            }
        }
    }
    a.assemble()
}
fn zero(n: usize, m: usize) -> Csr {
    Coo::new(n, m).assemble()
}
fn scalar(value: f64) -> Csr {
    matrix(1, 1, &[value])
}
fn oscillator(h: f64) -> SparseMidpoint {
    SparseMidpoint::new(
        SparseForms {
            mass: [scalar(2.0), scalar(0.3)],
            stiffness: [scalar(4.0), scalar(0.0)],
            damping: [scalar(0.0), scalar(0.0)],
            coupling: scalar(0.7),
        },
        h,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap()
}
fn close(x: f64, y: f64, tol: f64) {
    assert!(
        (x - y).abs() <= tol * (1.0 + x.abs() + y.abs()),
        "{x:e} != {y:e}"
    );
}

#[test]
fn closed_piston_compression_matches_analytic_wave_and_second_order_time_refinement() {
    // Physical sealed uniform cavity: Ma=V/(rho*c^2), B=piston area.
    // p=-B*u/Ma adds B^2/Ma to the original structural spring.
    let omega = (4.0_f64 / 2.0 + 0.7 * 0.7 / (2.0 * 0.3)).sqrt();
    let total = 1.1;
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let model = oscillator(total / steps as f64);
        let mut state = model
            .state(0.0, &[1.0], &[0.0], &[0.0], &[-0.7 / 0.3], &|| false)
            .unwrap();
        let initial = model.energy(&state, &|| false).unwrap().total_j();
        let mut exchanged = 0.0;
        for _ in 0..steps {
            let step = model.step(&state, &[0.0], &[0.0], &|| false).unwrap();
            exchanged += step.work.fluid_to_structure_work_j;
            close(step.work.structural_balance_residual_j(), 0.0, 1e-12);
            close(step.work.acoustic_balance_residual_j(), 0.0, 1e-12);
            close(step.work.interface_balance_residual_j(), 0.0, 1e-13);
            close(
                step.work.structural_dissipation_j + step.work.wall_dissipation_j,
                0.0,
                0.0,
            );
            state = step.state;
        }
        close(
            model.energy(&state, &|| false).unwrap().total_j(),
            initial,
            1e-12,
        );
        close(
            state.pressure_pa()[0],
            -0.7 / 0.3 * state.structural_displacement()[0],
            1e-12,
        );
        assert!(
            exchanged.abs() > 0.1,
            "nonzero reversible cavity/structure work must remain visible"
        );
        let exact = (omega * total).cos();
        let error = (state.structural_displacement()[0] - exact).abs();
        errors.push(error);
        close(
            state.pressure_potential_pa_s()[0],
            -(0.7 / 0.3) * (omega * total).sin() / omega,
            0.004,
        );
        close(
            state.structural_velocity()[0],
            -omega * (omega * total).sin(),
            0.004,
        );
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "midpoint phase must refine quadratically: {errors:?}"
    );
    assert!(errors[2] < 2e-4);
}

fn chain(n: usize, mass: f64, edge: f64) -> Csr {
    let mut a = Coo::new(n, n);
    for i in 0..n {
        a.push(i, i, mass);
    }
    for i in 1..n {
        a.push(i, i, edge);
        a.push(i - 1, i - 1, edge);
        a.push(i, i - 1, -edge);
        a.push(i - 1, i, -edge);
    }
    a.assemble()
}
fn large_forms() -> SparseForms {
    let ns = 48;
    let na = 64;
    let mut b = Coo::new(na, ns);
    for a in 0..na {
        b.push(a, a % ns, 0.07 + 0.0001 * a as f64);
        b.push(a, (a + 3) % ns, -0.02);
    }
    let mut wall = Coo::new(na, na);
    for (i, j) in [(0, 1), (na - 2, na - 1)] {
        wall.push(i, i, 0.09);
        wall.push(j, j, 0.09);
        wall.push(i, j, 0.03);
        wall.push(j, i, 0.03);
    }
    SparseForms {
        mass: [chain(ns, 1.7, 0.13), chain(na, 0.4, 0.02)],
        stiffness: [chain(ns, 0.0, 7.0), chain(na, 0.0, 1.3)],
        damping: [chain(ns, 0.03, 0.013), wall.assemble()],
        coupling: b.assemble(),
    }
}
#[test]
fn more_than_forty_two_full_coordinates_keep_consistent_mass_and_both_domain_work() {
    let model = SparseMidpoint::new(
        large_forms(),
        0.01,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap();
    let ns = model.structural_dofs();
    let na = model.pressure_dofs();
    assert_eq!(ns + na, 112);
    let u = (0..ns)
        .map(|i| 0.2 * (0.13 * i as f64).sin())
        .collect::<Vec<_>>();
    let v = (0..ns)
        .map(|i| 0.11 * (0.2 * i as f64).cos())
        .collect::<Vec<_>>();
    let psi = (0..na)
        .map(|i| 0.02 * (0.3 * i as f64).sin() + 0.7)
        .collect::<Vec<_>>();
    let p = (0..na)
        .map(|i| 0.14 * (0.12 * i as f64).cos() + 0.2)
        .collect::<Vec<_>>();
    let mut state = model.state(0.0, &u, &v, &psi, &p, &|| false).unwrap();
    let initial = model.energy(&state, &|| false).unwrap();
    let mut supply = 0.0;
    let mut heat = 0.0;
    let mut structure_work = 0.0;
    let mut fluid_work = 0.0;
    let mut exchange = 0.0;
    let mut hs = 0.0;
    let mut ha = 0.0;
    for i in 0..40 {
        let force = (0..ns)
            .map(|j| 0.01 * ((i + j) as f64 * 0.07).sin())
            .collect::<Vec<_>>();
        let source = (0..na)
            .map(|j| 0.008 * ((i + 3 * j) as f64 * 0.12).cos())
            .collect::<Vec<_>>();
        let step = model.step(&state, &force, &source, &|| false).unwrap();
        let w = step.work;
        assert!(w.structural_dissipation_j > 0.0 && w.wall_dissipation_j > 0.0);
        assert!(w.maximum_relative_residual < 2e-13);
        close(w.structural_balance_residual_j(), 0.0, 2e-13);
        close(w.acoustic_balance_residual_j(), 0.0, 2e-13);
        close(w.interface_balance_residual_j(), 0.0, 1e-14);
        // Independent original structural and acoustic endpoint equations.
        let vm = state
            .structural_velocity()
            .iter()
            .zip(step.state.structural_velocity())
            .map(|(&x, &y)| 0.5 * (x + y))
            .collect::<Vec<_>>();
        let pm = state
            .pressure_pa()
            .iter()
            .zip(step.state.pressure_pa())
            .map(|(&x, &y)| 0.5 * (x + y))
            .collect::<Vec<_>>();
        for a in 0..na {
            let (cols, vals) = model.forms.coupling.row(a);
            let outward: f64 = cols.iter().zip(vals).map(|(&s, &b)| b * vm[s]).sum();
            close(step.displaced_flow_m3_s[a], outward, 1e-14);
        }
        for s in 0..ns {
            let force: f64 = (0..na)
                .map(|a| model.forms.coupling.get(a, s) * pm[a])
                .sum();
            close(step.fluid_force[s], force, 1e-14);
        }
        supply += w.structural_input_work_j + w.acoustic_input_work_j;
        heat += w.structural_dissipation_j + w.wall_dissipation_j;
        structure_work += w.structural_input_work_j;
        fluid_work += w.acoustic_input_work_j;
        exchange += w.fluid_to_structure_work_j;
        hs += w.structural_dissipation_j;
        ha += w.wall_dissipation_j;
        state = step.state;
    }
    let final_energy = model.energy(&state, &|| false).unwrap();
    close(
        final_energy.total_j() - initial.total_j(),
        supply - heat,
        2e-12,
    );
    close(
        final_energy.structural_j() - initial.structural_j(),
        structure_work + exchange - hs,
        2e-12,
    );
    close(
        final_energy.acoustic_j() - initial.acoustic_j(),
        fluid_work - exchange - ha,
        2e-12,
    );
    assert!(exchange.abs() > 1e-5);
    assert_eq!(state.structural_displacement().len(), 48);
    assert_eq!(state.pressure_pa().len(), 64);
}

#[test]
fn structural_rigid_and_acoustic_uniform_modes_are_not_removed_or_regularized() {
    let ns = 45;
    let na = 2;
    let model = SparseMidpoint::new(
        SparseForms {
            mass: [Csr::identity(ns), matrix(na, na, &[0.5, 0.1, 0.1, 0.5])],
            stiffness: [zero(ns, ns), zero(na, na)],
            damping: [zero(ns, ns), zero(na, na)],
            coupling: zero(na, ns),
        },
        0.125,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap();
    let mut state = model
        .state(
            0.0,
            &vec![0.4; ns],
            &vec![0.2; ns],
            &[7.0, 7.0],
            &[0.3, 0.3],
            &|| false,
        )
        .unwrap();
    let e = model.energy(&state, &|| false).unwrap();
    for _ in 0..8 {
        let step = model
            .step(&state, &vec![0.0; ns], &[0.0, 0.0], &|| false)
            .unwrap();
        close(step.work.balance_residual_j(), 0.0, 1e-13);
        state = step.state;
    }
    for &u in state.structural_displacement() {
        close(u, 0.6, 2e-14);
    }
    for &v in state.structural_velocity() {
        close(v, 0.2, 2e-14);
    }
    for &psi in state.pressure_potential_pa_s() {
        close(psi, 7.3, 2e-14);
    }
    for &p in state.pressure_pa() {
        close(p, 0.3, 2e-14);
    }
    close(
        model.energy(&state, &|| false).unwrap().total_j(),
        e.total_j(),
        2e-14,
    );
}

#[test]
fn constant_force_and_volume_injection_do_their_own_physical_work() {
    let model = SparseMidpoint::new(
        SparseForms {
            mass: [scalar(2.0), scalar(0.4)],
            stiffness: [scalar(0.0), scalar(0.0)],
            damping: [scalar(0.0), scalar(0.0)],
            coupling: scalar(0.0),
        },
        0.2,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap();
    let step = model
        .step(&model.zero_state(), &[3.0], &[0.5], &|| false)
        .unwrap();
    close(step.state.structural_velocity()[0], 0.3, 1e-15);
    close(step.state.structural_displacement()[0], 0.03, 1e-15);
    close(step.state.pressure_pa()[0], 0.25, 1e-15);
    close(step.state.pressure_potential_pa_s()[0], 0.025, 1e-15);
    close(step.work.structural_input_work_j, 0.09, 1e-15);
    close(step.work.acoustic_input_work_j, 0.0125, 1e-15);
    close(step.work.fluid_to_structure_work_j, 0.0, 0.0);
}

fn small_forms() -> SparseForms {
    SparseForms {
        mass: [
            matrix(2, 2, &[2.0, 0.4, 0.4, 1.2]),
            matrix(3, 3, &[0.3, 0.04, 0.01, 0.04, 0.2, 0.02, 0.01, 0.02, 0.25]),
        ],
        stiffness: [matrix(2, 2, &[12.0, -3.0, -3.0, 7.0]), chain(3, 0.0, 2.0)],
        damping: [
            matrix(2, 2, &[0.2, -0.03, -0.03, 0.08]),
            matrix(3, 3, &[0.1, 0.02, 0.0, 0.02, 0.1, 0.0, 0.0, 0.0, 0.0]),
        ],
        coupling: matrix(3, 2, &[0.2, -0.1, -0.07, 0.2, 0.11, 0.05]),
    }
}
fn gaussian(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        for j in k..n {
            a.swap(k * n + j, pivot * n + j);
        }
        b.swap(k, pivot);
        for i in k + 1..n {
            let multiplier = a[i * n + k] / a[k * n + k];
            for j in k..n {
                a[i * n + j] -= multiplier * a[k * n + j];
            }
            b[i] -= multiplier * b[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            b[i] -= a[i * n + j] * b[j];
        }
        b[i] /= a[i * n + i];
    }
    b
}
#[test]
fn full_consistent_blocks_match_independent_nonsymmetric_original_equations() {
    let h = 0.07;
    let model = SparseMidpoint::new(
        small_forms(),
        h,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap();
    let q = [0.2, -0.1, 0.01, 0.04, -0.02];
    let w = [0.1, 0.2, 0.3, -0.1, 0.07];
    let f = [1.2, -0.2, 0.4, 0.02, -0.03];
    let state = model
        .state(0.0, &q[..2], &w[..2], &q[2..], &w[2..], &|| false)
        .unwrap();
    let result = model.step(&state, &f[..2], &f[2..], &|| false).unwrap();
    // Solve directly for NEXT velocities, WITHOUT the production acoustic
    // sign change, equilibration, sparse LDLT or midpoint-velocity unknowns.
    let n = 5;
    let mut a = vec![0.0; n * n];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        let d = usize::from(i >= 2);
        let off = if d == 0 { 0 } else { 2 };
        let r = i - off;
        rhs[i] = f[i];
        for j in off..if d == 0 { 2 } else { 5 } {
            let c = j - off;
            let m = model.forms.mass[d].get(r, c);
            let c0 = model.forms.damping[d].get(r, c);
            let k = model.forms.stiffness[d].get(r, c);
            a[i * n + j] = m / h + 0.5 * c0 + 0.25 * h * k;
            rhs[i] += (m / h - 0.5 * c0 - 0.25 * h * k) * w[j] - k * q[j];
        }
    }
    for i in 0..3 {
        for j in 0..2 {
            let b = model.forms.coupling.get(i, j);
            a[j * n + 2 + i] = -0.5 * b;
            rhs[j] += 0.5 * b * w[2 + i];
            a[(2 + i) * n + j] = 0.5 * b;
            rhs[2 + i] -= 0.5 * b * w[j];
        }
    }
    let expected = gaussian(a, rhs);
    for i in 0..n {
        close(result.state.velocity[i], expected[i], 2e-13);
        close(
            result.state.position[i],
            q[i] + 0.5 * h * (w[i] + expected[i]),
            2e-13,
        );
    }
    assert!(result.work.structural_dissipation_j > 0.0);
    assert!(result.work.wall_dissipation_j > 0.0);
}

#[test]
fn diagonal_coordinate_units_preserve_full_response_and_physical_accounts() {
    let reference = SparseMidpoint::new(
        small_forms(),
        0.04,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap();
    let mut changed = small_forms();
    let scales = [1e8, 0.003, 1e-5, 200.0, 0.7];
    for d in 0..2 {
        let off = if d == 0 { 0 } else { 2 };
        for a in [
            &mut changed.mass[d],
            &mut changed.stiffness[d],
            &mut changed.damping[d],
        ] {
            for i in 0..a.nrows() {
                let columns = a.row(i).0.to_vec();
                for (&j, value) in columns.iter().zip(a.row_values_mut(i)) {
                    let (x, y) = if i < j { (i, j) } else { (j, i) };
                    *value = (*value * scales[off + x]) * scales[off + y];
                }
            }
        }
    }
    for i in 0..3 {
        let columns = changed.coupling.row(i).0.to_vec();
        for (&j, value) in columns.iter().zip(changed.coupling.row_values_mut(i)) {
            *value *= scales[2 + i] * scales[j];
        }
    }
    let transformed = SparseMidpoint::new(
        changed,
        0.04,
        SparseStructuralTimeLimits::default(),
        &|| false,
    )
    .unwrap();
    let q = [0.2, -0.1, 0.01, 0.04, -0.02];
    let v = [0.1, 0.2, 0.3, -0.1, 0.07];
    let f = [1.2, -0.2, 0.4, 0.02, -0.03];
    let q2 = (0..5).map(|i| q[i] / scales[i]).collect::<Vec<_>>();
    let v2 = (0..5).map(|i| v[i] / scales[i]).collect::<Vec<_>>();
    let f2 = (0..5).map(|i| f[i] * scales[i]).collect::<Vec<_>>();
    let a = reference
        .state(0.0, &q[..2], &v[..2], &q[2..], &v[2..], &|| false)
        .unwrap();
    let b = transformed
        .state(0.0, &q2[..2], &v2[..2], &q2[2..], &v2[2..], &|| false)
        .unwrap();
    let x = reference.step(&a, &f[..2], &f[2..], &|| false).unwrap();
    let y = transformed.step(&b, &f2[..2], &f2[2..], &|| false).unwrap();
    for i in 0..5 {
        close(x.state.position[i], y.state.position[i] * scales[i], 1e-12);
        close(x.state.velocity[i], y.state.velocity[i] * scales[i], 1e-12);
    }
    for (a, b) in [
        (x.work.final_energy.total_j(), y.work.final_energy.total_j()),
        (
            x.work.structural_input_work_j,
            y.work.structural_input_work_j,
        ),
        (x.work.acoustic_input_work_j, y.work.acoustic_input_work_j),
        (
            x.work.fluid_to_structure_work_j,
            y.work.fluid_to_structure_work_j,
        ),
        (
            x.work.structural_dissipation_j,
            y.work.structural_dissipation_j,
        ),
        (x.work.wall_dissipation_j, y.work.wall_dissipation_j),
    ] {
        close(a, b, 1e-12);
    }
}

#[test]
fn nonpassive_nonsymmetric_and_oversized_original_models_refuse_without_repair() {
    let mut forms = small_forms();
    forms.stiffness[0] = matrix(2, 2, &[1.0, 2.0, 2.0, 1.0]);
    assert!(
        SparseMidpoint::new(forms, 0.1, SparseStructuralTimeLimits::default(), &|| false).is_err()
    );
    let mut forms = small_forms();
    forms.damping[0] = matrix(2, 2, &[1.0, 0.3, -0.3, 1.0]);
    assert!(matches!(
        SparseMidpoint::new(forms, 0.1, SparseStructuralTimeLimits::default(), &|| false),
        Err(SparseStructuralTimeError::Invalid(_))
    ));
    let mut forms = small_forms();
    forms.mass[0] = matrix(2, 2, &[0.0, 0.0, 0.0, 1.0]);
    assert!(
        SparseMidpoint::new(forms, 0.1, SparseStructuralTimeLimits::default(), &|| false).is_err()
    );
    let mut limits = SparseStructuralTimeLimits::default();
    limits.maximum_coordinates = 4;
    assert!(matches!(
        SparseMidpoint::new(small_forms(), 0.1, limits, &|| false),
        Err(SparseStructuralTimeError::Limit { .. })
    ));
    let mut limits = SparseStructuralTimeLimits::default();
    limits.maximum_factor_bytes = 1;
    assert!(matches!(
        SparseMidpoint::new(small_forms(), 0.1, limits, &|| false),
        Err(SparseStructuralTimeError::Factor(_))
    ));
    let mut limits = SparseStructuralTimeLimits::default();
    limits.factor.max_numeric_work = 1;
    assert!(SparseMidpoint::new(small_forms(), 0.1, limits, &|| false).is_err());
}

#[test]
fn cancellation_foreign_state_and_nonfinite_source_leave_exact_retry_available() {
    let polls = Cell::new(0usize);
    let cancelled = || {
        polls.set(polls.get() + 1);
        polls.get() > 30
    };
    assert!(matches!(
        SparseMidpoint::new(
            large_forms(),
            0.01,
            SparseStructuralTimeLimits::default(),
            &cancelled
        ),
        Err(SparseStructuralTimeError::Cancelled)
    ));
    let model = oscillator(0.1);
    let other = oscillator(0.1);
    let state = model
        .state(0.0, &[0.7], &[0.2], &[0.3], &[0.1], &|| false)
        .unwrap();
    let position = state.position.clone();
    let velocity = state.velocity.clone();
    let first = model.step(&state, &[0.2], &[0.03], &|| false).unwrap();
    let polls = Cell::new(0usize);
    let cancelled = || {
        polls.set(polls.get() + 1);
        polls.get() > 8
    };
    assert!(matches!(
        model.step(&state, &[0.2], &[0.03], &cancelled),
        Err(SparseStructuralTimeError::Cancelled)
    ));
    assert!(other.step(&state, &[0.2], &[0.03], &|| false).is_err());
    assert!(model.step(&state, &[f64::NAN], &[0.03], &|| false).is_err());
    assert_eq!(state.position, position);
    assert_eq!(state.velocity, velocity);
    assert_eq!(state.time_s, 0.0);
    let retry = model.step(&state, &[0.2], &[0.03], &|| false).unwrap();
    assert_eq!(first.state.position, retry.state.position);
    assert_eq!(first.state.velocity, retry.state.velocity);
}
