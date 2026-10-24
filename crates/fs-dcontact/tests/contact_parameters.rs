//! G0/G3 parameter derivatives against the original storage and loss owners.
use fs_dcontact::{ContactStorage, DContactError, Obstacle};
use fs_phs::Storage;

struct Zero;
impl Storage for Zero {
    fn hamiltonian(&self, _: &[f64]) -> f64 {
        0.0
    }
    fn gradient(&self, _: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
}

fn point(k: f64, chi: f64, gap: f64, weight: f64, alpha: f64) -> Obstacle {
    Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![gap],
        vec![weight],
        k,
        alpha,
        "synthetic contact calibration fixture".into(),
    )
    .unwrap()
    .with_internal_loss(chi)
    .unwrap()
}

// Exercise the pre-existing conservative gradient and negative loss force;
// do not reconstruct the new parameter derivative's formula as its oracle.
fn objective(law: &Obstacle, n: usize, x: &[f64], v: &[f64], seed: &[f64]) -> f64 {
    let storage = ContactStorage::new(Box::new(Zero), n, vec![law.clone()]).unwrap();
    let mut gradient = vec![0.0; x.len()];
    storage.gradient(x, &mut gradient);
    let loss = law.dissipative_modal_forces(n, x, v).unwrap();
    (0..n).map(|k| seed[k] * (gradient[2 * k] - loss[k])).sum()
}

fn perturbed(law: &Obstacle, n: usize, parameter: usize, increment: f64) -> Obstacle {
    let mut k = law.stiffness();
    let mut chi = law.internal_loss();
    let mut gaps = law.gaps().to_vec();
    match parameter {
        0 => k += increment,
        1 => chi += increment,
        _ => gaps[parameter - 2] += increment,
    }
    Obstacle::new(
        law.collocation().to_vec(),
        law.n_points(),
        n,
        gaps,
        law.weights().to_vec(),
        k,
        law.alpha(),
        law.provenance().into(),
    )
    .unwrap()
    .with_internal_loss(chi)
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-7 * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn check_central(law: &Obstacle, n: usize, x: &[f64], v: &[f64], seed: &[f64]) {
    let mut derivative = vec![f64::NAN; 2 + law.n_points()];
    law.resistance_parameter_vjp_into(n, x, v, seed, &mut derivative)
        .unwrap();
    for (i, actual) in derivative.into_iter().enumerate() {
        let value = match i {
            0 => law.stiffness(),
            1 => law.internal_loss(),
            _ => law.gaps()[i - 2],
        };
        let h = 1e-6 * (1.0 + value.abs());
        let plus = objective(&perturbed(law, n, i, h), n, x, v, seed);
        let minus = objective(&perturbed(law, n, i, -h), n, x, v, seed);
        close(actual, (plus - minus) / (2.0 * h));
    }
}

#[test]
fn smooth_loading_and_unloading_match_force_finite_differences() {
    for alpha in [1.0, 1.5, 2.4] {
        for rate in [-0.4, 0.2] {
            check_central(
                &point(12.0, 1.5, 0.1, 2.3, alpha),
                1,
                &[0.4, 19.0],
                &[rate],
                &[-0.7],
            );
        }
    }
}

#[test]
fn multipoint_vjp_preserves_physical_measures_and_mixed_branches() {
    let law = Obstacle::new(
        vec![
            1.0, 0.2, -0.1, -0.5, 1.1, 0.4, 0.0, 1.0, 0.0, 0.7, -0.4, 0.8,
        ],
        4,
        3,
        vec![0.05, 0.2, -0.5, 0.1],
        vec![0.7, 2.3, 1.1, 0.3],
        17.0,
        1.5,
        "synthetic mixed active/open/capped sites".into(),
    )
    .unwrap()
    .with_internal_loss(1.5)
    .unwrap();
    // Two active sites, one open site, and one engaged but unloading-capped
    // site. Momentum values and the final internal coordinate are irrelevant.
    let x = [0.6, 15.0, -0.2, -23.0, 0.3, 4.0, 7.0];
    let velocity = [0.8, -1.6, 0.4];
    let seed = [-0.3, 0.6, 0.8];
    check_central(&law, 3, &x, &velocity, &seed);
    let mut derivative = [0.0; 6];
    law.resistance_parameter_vjp_into(3, &x, &velocity, &seed, &mut derivative)
        .unwrap();
    assert_eq!(derivative[3], 0.0); // Open site's gap.
    assert_eq!(derivative[4], 0.0); // Strictly capped site's gap.
}

#[test]
fn zero_coefficients_keep_the_admitted_side_derivatives() {
    let x = [0.4, 0.0];
    let v = [0.6];
    let seed = [1.7];
    for (law, parameter) in [
        (point(0.0, 0.3, 0.1, 2.0, 1.5), 0),
        (point(12.0, 0.0, 0.1, 2.0, 1.5), 1),
    ] {
        let mut derivative = [0.0; 3];
        law.resistance_parameter_vjp_into(1, &x, &v, &seed, &mut derivative)
            .unwrap();
        let h = 1e-6;
        let fd = (objective(&perturbed(&law, 1, parameter, h), 1, &x, &v, &seed)
            - objective(&law, 1, &x, &v, &seed))
            / h;
        close(derivative[parameter], fd);
        assert!(derivative[parameter] > 0.0);
        if parameter == 0 {
            assert_eq!(&derivative[1..], &[0.0, 0.0]);
        }
    }
    let mut zero_weight = [9.0; 3];
    point(12.0, 0.3, 0.1, 0.0, 1.5)
        .resistance_parameter_vjp_into(1, &x, &v, &seed, &mut zero_weight)
        .unwrap();
    assert_eq!(zero_weight, [0.0; 3]);
}

#[test]
fn separation_and_cap_corners_follow_the_existing_one_sided_policy() {
    let mut derivative = [9.0; 3];
    for q in [-0.2, 0.0] {
        point(12.0, 0.3, 0.0, 2.0, 1.0)
            .resistance_parameter_vjp_into(1, &[q, 0.0], &[0.4], &[0.7], &mut derivative)
            .unwrap();
        assert_eq!(derivative, [0.0; 3]);
    }
    let law = point(12.0, 2.0, 0.0, 2.0, 1.5);
    let (x, v, seed) = ([0.25, 0.0], [-0.5], [0.4]);
    law.resistance_parameter_vjp_into(1, &x, &v, &seed, &mut derivative)
        .unwrap();
    assert_eq!(derivative[0], 0.0);
    assert_eq!(derivative[2], 0.0);
    assert!(derivative[1] < 0.0);
    let h = 1e-6;
    close(
        derivative[1],
        (objective(&law, 1, &x, &v, &seed)
            - objective(&perturbed(&law, 1, 1, -h), 1, &x, &v, &seed))
            / h,
    );
    let capped = point(12.0, 3.0, 0.0, 2.0, 1.5);
    capped
        .resistance_parameter_vjp_into(1, &x, &v, &seed, &mut derivative)
        .unwrap();
    assert_eq!(derivative, [0.0; 3]);
    check_central(&capped, 1, &x, &v, &seed);
}

#[test]
fn shape_and_nonfinite_input_refusals_leave_output_untouched() {
    let law = point(12.0, 0.3, 0.0, 2.0, 1.5);
    for (n, x, v, seed) in [
        (1, vec![0.2], vec![0.1], vec![1.0]),
        (1, vec![0.2, 0.0], vec![], vec![1.0]),
        (1, vec![0.2, 0.0], vec![0.1], vec![1.0, 2.0]),
        (usize::MAX, vec![0.2, 0.0], vec![0.1], vec![1.0]),
    ] {
        let mut out = [19.0; 3];
        assert!(matches!(
            law.resistance_parameter_vjp_into(n, &x, &v, &seed, &mut out),
            Err(DContactError::Shape { .. })
        ));
        assert_eq!(out, [19.0; 3]);
    }
    let mut wrong_output = [19.0; 4];
    assert!(
        law.resistance_parameter_vjp_into(1, &[0.2, 0.0], &[0.1], &[1.0], &mut wrong_output)
            .is_err()
    );
    assert_eq!(wrong_output, [19.0; 4]);
    for (x, v, seed) in [
        (vec![0.2, 0.0, f64::NAN], [0.1], [1.0]),
        (vec![0.2, 0.0], [f64::INFINITY], [1.0]),
        (vec![0.2, 0.0], [0.1], [f64::NAN]),
    ] {
        let mut out = [19.0; 3];
        assert!(matches!(
            law.resistance_parameter_vjp_into(1, &x, &v, &seed, &mut out),
            Err(DContactError::Parameter { .. })
        ));
        assert_eq!(out, [19.0; 3]);
    }
}

#[test]
fn malformed_raw_obstacles_refuse_before_indexing_or_output_changes() {
    for (rows, gaps, weights, k, alpha) in [
        (vec![], vec![0.0], vec![1.0], 1.0, 1.5),
        (vec![1.0], vec![], vec![1.0], 1.0, 1.5),
        (vec![1.0], vec![0.0], vec![], 1.0, 1.5),
        (vec![f64::NAN], vec![0.0], vec![1.0], 1.0, 1.5),
        (vec![1.0], vec![f64::INFINITY], vec![1.0], 1.0, 1.5),
        (vec![1.0], vec![0.0], vec![-1.0], 1.0, 1.5),
        (vec![1.0], vec![0.0], vec![1.0], -1.0, 1.5),
        (vec![1.0], vec![0.0], vec![1.0], 1.0, 0.5),
    ] {
        let law = Obstacle::from_raw_parts(rows, 1, gaps, weights, k, alpha, "raw mutation".into());
        let mut out = [19.0; 3];
        assert!(
            law.resistance_parameter_vjp_into(1, &[0.2, 0.0], &[0.1], &[1.0], &mut out)
                .is_err()
        );
        assert_eq!(out, [19.0; 3]);
    }
}

#[test]
fn arithmetic_overflow_is_a_refusal_and_empty_owner_has_zero_derivatives() {
    let mut out = [0.0; 3];
    assert!(
        point(f64::MAX / 2.0, 0.3, 0.0, 1.0, 2.0)
            .resistance_parameter_vjp_into(1, &[2.0, 0.0], &[0.1], &[1.0], &mut out)
            .is_err()
    );
    let empty = Obstacle::new(
        vec![],
        0,
        0,
        vec![],
        vec![],
        0.0,
        1.0,
        "empty contact".into(),
    )
    .unwrap();
    let mut no_parameters = [19.0; 2];
    empty
        .resistance_parameter_vjp_into(0, &[], &[], &[], &mut no_parameters)
        .unwrap();
    assert_eq!(no_parameters, [0.0; 2]);
}
