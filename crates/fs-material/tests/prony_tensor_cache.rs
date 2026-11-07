//! A public branch coefficient edit must reach every constitutive consumer.
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::visco_tensor::TensorProjectionBudget;

fn diagonal(value: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0.0 }))
}

fn changed_tensor() -> Mandel6 {
    let mut c = [[0.0; 6]; 6];
    c[0][0] = 9.0;
    c[1][1] = 4.0;
    c[0][1] = -2.0;
    c[1][0] = -2.0;
    c
}

fn edited_branch(c: Mandel6, from_rows: bool) -> PronyBranch {
    let mut branch = if from_rows {
        PronyBranch::from_gram_rows(diagonal(1.0), 0.3).unwrap()
    } else {
        PronyBranch::new(diagonal(1.0), 0.3).unwrap()
    };
    branch.c_pa = c;
    branch
}

fn admit(branch: PronyBranch, mandel_constructor: bool) -> TensorProny {
    if mandel_constructor {
        TensorProny::new_mandel(diagonal(2.0), vec![branch], 1.0).unwrap()
    } else {
        TensorProny::new(diagonal(2.0), None, vec![branch], (0.0, 1e6)).unwrap()
    }
}

fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn apply(c: &Mandel6, x: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| dot(&c[i], x))
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * actual.abs().max(expected.abs()).max(1e-250),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn budget() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 2,
        max_scalar_terms: 10_000,
    }
}

// Analytic single-ramp result from the CURRENT raw source matrix. The integral
// of (1-exp(-t/tau))^2 is evaluated independently of admitted Gram factors and
// the implementation's PronyRamp heat routine.
fn check_first_ramp(law: &TensorProny, source: &Mandel6) {
    let strain = [0.02, -0.01, 0.0, 0.0, 0.0, 0.0];
    let dt = 0.4;
    let r = dt / 0.3_f64;
    let b = -(-r).exp_m1() / r;
    let mean = (1.0 - b) / r;
    let branch_stress = apply(source, &strain);
    let branch_work = dot(&strain, &branch_stress);
    let expected_heat =
        branch_work * (r + 2.0 * (-r).exp_m1() - 0.5 * (-2.0 * r).exp_m1()) / (r * r);
    let expected_stored = dot(&strain, &strain) + 0.5 * b * b * branch_work;
    let mut direct = law.state();
    let a = law.step_with_mean(&mut direct, strain, dt).unwrap();
    // Only normal strains are used, so physical tensor-shear and Mandel
    // coordinates coincide without any conversion in this oracle.
    let mut history = law.history();
    let h = history.advance(strain, dt, || false).unwrap();
    for i in 0..6 {
        let endpoint = 2.0 * strain[i] + b * branch_stress[i];
        let average = strain[i] + mean * branch_stress[i];
        close(a.stress[i], endpoint);
        close(h.stress[i], endpoint);
        close(a.mean_stress[i], average);
        close(h.mean_stress[i], average);
    }
    close(direct.stored, expected_stored);
    close(h.stored_energy, expected_stored);
    close(a.step_dissipated, expected_heat);
    close(h.dissipated, expected_heat);
    close(a.step_work, expected_stored + expected_heat);
    close(h.work, expected_stored + expected_heat);
}

#[test]
fn edited_cached_coefficients_drive_both_physical_history_apis() {
    let c = changed_tensor();
    for from_rows in [false, true] {
        for mandel_constructor in [false, true] {
            let actual = admit(edited_branch(c, from_rows), mandel_constructor);
            let fresh = admit(PronyBranch::unvalidated(c, 0.3), mandel_constructor);
            check_first_ramp(&actual, &c);
            check_first_ramp(&fresh, &c);
            assert_eq!(actual.branches(), fresh.branches());
            assert_eq!(actual.factors(), fresh.factors());
        }
    }
}

#[test]
fn edited_source_survives_projection_rotation_and_harmonic_queries() {
    let c = changed_tensor();
    let columns = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let quarter_turn = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    for mandel_constructor in [false, true] {
        let law = admit(edited_branch(c, true), mandel_constructor);
        for rotated in [false, true] {
            let (law, source) = if rotated {
                let mut rotated_c = c;
                rotated_c[0][0] = c[1][1];
                rotated_c[1][1] = c[0][0];
                (law.rotated(quarter_turn).unwrap(), rotated_c)
            } else {
                (law.clone(), c)
            };
            check_first_ramp(&law, &source);
            let projection = law
                .project_strain(&columns, 0.25, budget(), || false)
                .unwrap();
            for i in 0..2 {
                for j in 0..2 {
                    let actual: f64 = projection.branches()[0]
                        .rows()
                        .iter()
                        .map(|row| row[i] * row[j])
                        .sum();
                    close(actual, 0.25 * source[i][j]);
                }
            }
            for omega in [0.0, 0.1, 2.0, 1e4] {
                let (storage, loss) = law.complex_modulus(omega).unwrap();
                let z = omega * 0.3;
                for i in 0..6 {
                    for j in 0..6 {
                        let equilibrium = if i == j { 2.0 } else { 0.0 };
                        close(
                            storage[i][j],
                            equilibrium + source[i][j] * z * z / (1.0 + z * z),
                        );
                        close(loss[i][j], source[i][j] * z / (1.0 + z * z));
                    }
                }
            }
        }
    }
}

#[test]
fn removing_a_cached_branch_stiffness_removes_its_loss_without_removing_history() {
    for mandel_constructor in [false, true] {
        let law = admit(edited_branch([[0.0; 6]; 6], true), mandel_constructor);
        check_first_ramp(&law, &[[0.0; 6]; 6]);
        assert_eq!(law.branches().len(), 1);
        assert_eq!(law.branches()[0].tau_s, 0.3);
        assert_eq!(law.complex_modulus(2.0).unwrap().1, [[0.0; 6]; 6]);
    }
}

#[test]
fn a_cached_factor_cannot_admit_an_active_or_nonfinite_edited_source() {
    let mut asymmetric = changed_tensor();
    asymmetric[0][1] = 0.0;
    for source in [
        diagonal(-1.0),
        diagonal(f64::NAN),
        diagonal(f64::INFINITY),
        asymmetric,
    ] {
        for from_rows in [false, true] {
            let branch = edited_branch(source, from_rows);
            assert!(TensorProny::new(diagonal(2.0), None, vec![branch], (0.0, 1e6)).is_err());
            assert!(TensorProny::new_mandel(diagonal(2.0), vec![branch], 1.0).is_err());
        }
    }
}

#[test]
fn unchanged_factor_defined_weak_channels_keep_their_original_coordinates() {
    let a = 2.0_f64.powi(-80);
    let mut rows = [[0.0; 6]; 6];
    rows[5][0] = 2.0;
    rows[5][1] = -1.0;
    rows[3][2] = a;
    let branch = PronyBranch::from_gram_rows(rows, 0.3).unwrap();
    let columns = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    ];
    for mandel_constructor in [false, true] {
        let law = admit(branch, mandel_constructor);
        let expected_factor = core::array::from_fn(|i| core::array::from_fn(|j| rows[j][i]));
        assert_eq!(law.factors()[0], expected_factor);
        let p = law
            .project_strain(&columns, 1.0, budget(), || false)
            .unwrap();
        for i in 0..6 {
            assert_eq!(p.branches()[0].rows()[i], [rows[i][0], rows[i][2]]);
        }
        let mut direct = law.state();
        let mut history = law.history();
        let mut expected_q = 0.0;
        let mut before = 0.0;
        for (next, dt) in [(0.02, 0.4), (-0.01, 0.1), (0.0, 0.7)] {
            let strain = [0.0, 0.0, next, 0.0, 0.0, 0.0];
            let ratio = dt / 0.3_f64;
            expected_q =
                (-ratio).exp() * expected_q - (-ratio).exp_m1() / ratio * a * (next - before);
            let plain = law.step_with_mean(&mut direct, strain, dt).unwrap();
            let factored = history.advance(strain, dt, || false).unwrap();
            close(plain.step_dissipated, factored.dissipated);
            close(
                history.branch_stored_energy(),
                0.5 * expected_q * expected_q,
            );
            assert!(plain.step_dissipated > 0.0 && factored.dissipated > 0.0);
            before = next;
        }
        close(direct.stored, history.branch_stored_energy());
    }
}
