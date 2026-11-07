//! Exact material-frequency slopes of complex, constrained 3-D stress laws.
//! Synthetic anisotropic fixtures; no measured-material or assembled-shell claim.

use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco::ViscoError;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;

const PLANES: [MaterialPlane; 4] = [
    MaterialPlane::Axes12,
    MaterialPlane::Axes21,
    MaterialPlane::Axes13,
    MaterialPlane::Axes23,
];

fn outer(diagonal: f64, v: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| v[i] * v[j] + if i == j { diagonal } else { 0.0 })
    })
}

fn fixture() -> (TensorProny, Mandel6, Mandel6, [(Mandel6, f64); 2]) {
    let equilibrium = outer(5.0, [0.4, -0.5, 0.7, 0.2, -0.3, 0.6]);
    let kelvin = outer(0.02, [0.02, -0.03, 0.04, -0.05, 0.06, -0.07]);
    let branches = [
        (outer(0.6, [1.0, -2.0, 0.5, 0.25, -0.75, 1.5]), 0.125),
        (outer(0.9, [-0.5, 0.75, 1.0, 1.25, 0.5, -1.0]), 2.0),
    ];
    let law = TensorProny::new(
        equilibrium,
        Some(kelvin),
        branches
            .iter()
            .map(|&(c, tau)| PronyBranch::new(c, tau).unwrap())
            .collect(),
        (0.0, 100.0),
    )
    .unwrap();
    (law, equilibrium, kelvin, branches)
}

fn indices(plane: MaterialPlane) -> ([usize; 3], [usize; 3]) {
    match plane {
        MaterialPlane::Axes12 => ([0, 1, 3], [2, 4, 5]),
        MaterialPlane::Axes21 => ([1, 0, 3], [2, 4, 5]),
        MaterialPlane::Axes13 => ([0, 2, 5], [1, 3, 4]),
        MaterialPlane::Axes23 => ([1, 2, 4], [0, 3, 5]),
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(actual.re.is_finite() && actual.im.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

fn branch(omega: f64, tau: f64) -> (C64, C64) {
    let x = omega * tau;
    let d = 1.0 + x * x;
    (
        C64::new(x * x / d, x / d),
        C64::new(2.0 * tau * x / (d * d), tau * (1.0 - x * x) / (d * d)),
    )
}

// Independent moderate-scale 6x6 constraint solve for the manufactured test.
// In-plane strain rows are prescribed identities; out-of-plane rows enforce
// zero physical stress. This neither forms nor differentiates a 3x3 Schur map.
fn constrained_solve(mut a: [[C64; 6]; 6], mut b: [C64; 6]) -> [C64; 6] {
    for k in 0..6 {
        let mut pivot = k;
        for i in k + 1..6 {
            if a[i][k].abs() > a[pivot][k].abs() {
                pivot = i;
            }
        }
        assert!(a[pivot][k].abs() > 1e-10);
        a.swap(k, pivot);
        b.swap(k, pivot);
        let diagonal = a[k][k];
        for j in k..6 {
            a[k][j] = a[k][j] / diagonal;
        }
        b[k] = b[k] / diagonal;
        for i in 0..6 {
            if i == k {
                continue;
            }
            let factor = a[i][k];
            for j in k..6 {
                a[i][j] = a[i][j] - factor * a[k][j];
            }
            b[i] = b[i] - factor * b[k];
        }
    }
    b
}

#[test]
fn anisotropic_slopes_match_independent_full_stress_constraints_in_every_plane() {
    let (law, equilibrium, kelvin, branches) = fixture();
    let commutator: f64 = (0..6)
        .map(|k| {
            branches[0].0[0][k] * branches[1].0[k][1] - branches[1].0[0][k] * branches[0].0[k][1]
        })
        .sum();
    assert!(
        commutator.abs() > 0.1,
        "branches must not share principal axes"
    );
    for omega in [0.0, 0.2, 0.5, 1.0, 4.0, 8.0, 13.0] {
        let mut c = [[C64::ZERO; 6]; 6];
        let mut dc = [[C64::ZERO; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                c[i][j] = C64::new(equilibrium[i][j], omega * kelvin[i][j]);
                dc[i][j] = C64::new(0.0, kelvin[i][j]);
                for &(coefficient, tau) in &branches {
                    let (value, slope) = branch(omega, tau);
                    c[i][j] = c[i][j] + value.scale(coefficient[i][j]);
                    dc[i][j] = dc[i][j] + slope.scale(coefficient[i][j]);
                }
            }
        }
        for plane in PLANES {
            let (inside, outside) = indices(plane);
            let engineering = [1.0, 1.0, core::f64::consts::FRAC_1_SQRT_2];
            let (re, im) = law
                .plane_stress_complex_frequency_derivative(plane, omega)
                .unwrap();
            let mut constraint = c;
            for row in inside {
                constraint[row] = [C64::ZERO; 6];
                constraint[row][row] = C64::from_re(1.0);
            }
            for column in 0..3 {
                let mut imposed = [C64::ZERO; 6];
                imposed[inside[column]] = C64::from_re(engineering[column]);
                let strain = constrained_solve(constraint, imposed);
                let mut derivative_force = [C64::ZERO; 6];
                for row in outside {
                    derivative_force[row] =
                        (0..6).fold(C64::ZERO, |sum, k| sum - dc[row][k] * strain[k]);
                }
                let strain_rate = constrained_solve(constraint, derivative_force);
                let stress_rate: [C64; 6] = core::array::from_fn(|row| {
                    (0..6).fold(C64::ZERO, |sum, k| {
                        sum + dc[row][k] * strain[k] + c[row][k] * strain_rate[k]
                    })
                });
                for row in outside {
                    close(stress_rate[row], C64::ZERO, 2e-13);
                }
                for row in 0..3 {
                    let at = 3 * row + column;
                    close(
                        C64::new(re[at], im[at]),
                        stress_rate[inside[row]].scale(engineering[row]),
                        3e-13,
                    );
                    assert_eq!(re[at], re[3 * column + row]);
                    assert_eq!(im[at], im[3 * column + row]);
                }
            }
        }
    }
}

#[test]
fn all_reduced_entries_match_native_frequency_differences_around_relaxation_peaks() {
    let (law, _, _, _) = fixture();
    for plane in PLANES {
        for omega in [0.499_f64, 0.5, 0.501, 7.99, 8.0, 8.01] {
            let h = 1e-5 * omega.max(1.0);
            let (pr, pi) = law.plane_stress_complex(plane, omega + h).unwrap();
            let (mr, mi) = law.plane_stress_complex(plane, omega - h).unwrap();
            let (dr, di) = law
                .plane_stress_complex_frequency_derivative(plane, omega)
                .unwrap();
            for i in 0..9 {
                close(
                    C64::new(dr[i], di[i]),
                    C64::new((pr[i] - mr[i]) / (2.0 * h), (pi[i] - mi[i]) / (2.0 * h)),
                    3e-8,
                );
            }
        }
    }
}

fn isotropic(lambda: f64, mu: f64) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                lambda + if i == j { 2.0 * mu } else { 0.0 }
            } else if i == j {
                2.0 * mu
            } else {
                0.0
            }
        })
    })
}

#[test]
fn isotropic_complex_lame_formula_retains_signed_loss_slopes_and_engineering_shear() {
    let law = TensorProny::new(
        isotropic(4.0, 3.0),
        None,
        vec![
            PronyBranch::new(isotropic(2.0, 0.7), 0.2).unwrap(),
            PronyBranch::new(isotropic(0.8, 1.4), 1.5).unwrap(),
        ],
        (0.0, 100.0),
    )
    .unwrap();
    for omega in [0.0, 0.2, 1.0, 5.0, 20.0] {
        let (mut lambda, mut mu) = (C64::from_re(4.0), C64::from_re(3.0));
        let (mut dl, mut dm) = (C64::ZERO, C64::ZERO);
        for (l, g, tau) in [(2.0, 0.7, 0.2), (0.8, 1.4, 1.5)] {
            let (value, slope) = branch(omega, tau);
            lambda = lambda + value.scale(l);
            mu = mu + value.scale(g);
            dl = dl + slope.scale(l);
            dm = dm + slope.scale(g);
        }
        let a = lambda + mu.scale(2.0);
        let da = dl + dm.scale(2.0);
        let eliminated = (lambda * dl).scale(2.0) / a - lambda * lambda * da / (a * a);
        let expected = [
            da - eliminated,
            dl - eliminated,
            C64::ZERO,
            dl - eliminated,
            da - eliminated,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            dm,
        ];
        for plane in PLANES {
            let (re, im) = law
                .plane_stress_complex_frequency_derivative(plane, omega)
                .unwrap();
            for i in 0..9 {
                close(C64::new(re[i], im[i]), expected[i], 3e-13);
            }
            if omega == 20.0 {
                assert!(im[0] < 0.0 && im[8] < 0.0, "loss slopes must remain signed");
            }
        }
    }
}

#[test]
fn plane_derivative_does_not_require_representable_effective_viscosity_slope() {
    let equilibrium = outer(5.0, [0.4, -0.5, 0.7, 0.2, -0.3, 0.6]);
    let coefficient = outer(0.6, [1.0, -2.0, 0.5, 0.25, -0.75, 1.5]);
    let make = |tau| {
        TensorProny::new(
            equilibrium,
            None,
            vec![PronyBranch::new(coefficient, tau).unwrap()],
            (0.0, 2.0),
        )
        .unwrap()
    };
    let tau = 2.0_f64.powi(600);
    let slow = make(tau);
    let reference = make(1.0);
    assert!(slow.complex_modulus(1.0 / tau).is_ok());
    assert!(
        slow.complex_modulus_frequency_derivative(1.0 / tau)
            .is_err()
    );
    for plane in PLANES {
        let (re, im) = slow
            .plane_stress_complex_frequency_derivative(plane, 1.0 / tau)
            .unwrap();
        let (rr, ri) = reference
            .plane_stress_complex_frequency_derivative(plane, 1.0)
            .unwrap();
        for i in 0..9 {
            close(
                C64::new(re[i] / tau, im[i] / tau),
                C64::new(rr[i], ri[i]),
                3e-13,
            );
        }
    }
}

#[test]
fn huge_eliminated_normal_slope_does_not_erase_weak_in_plane_response() {
    let mut slow = [[0.0; 6]; 6];
    slow[2][2] = 1.0;
    let mut weak = [[0.0; 6]; 6];
    weak[0][0] = 1e-100;
    let equilibrium =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    let law = TensorProny::new(
        equilibrium,
        None,
        vec![
            PronyBranch::new(slow, 1e300).unwrap(),
            PronyBranch::new(weak, 1.0).unwrap(),
        ],
        (0.0, 1.0),
    )
    .unwrap();
    let full = law.complex_modulus_frequency_derivative(0.0).unwrap();
    assert_eq!(full.loss_pa_s[2][2], 1e300);
    assert_eq!(full.loss_pa_s[0][0], 1e-100);
    let (storage, loss) = law
        .plane_stress_complex_frequency_derivative(MaterialPlane::Axes12, 0.0)
        .unwrap();
    assert_eq!(storage, [0.0; 9]);
    assert!((loss[0] / 1e-100 - 1.0).abs() < 3e-14);
    assert!(loss[1..].iter().all(|value| *value == 0.0));
}

#[test]
fn constant_law_and_typed_band_admission_keep_the_original_material_unchanged() {
    let constant = TensorProny::new(isotropic(4.0, 3.0), None, vec![], (0.0, 100.0)).unwrap();
    for plane in PLANES {
        for omega in [0.0, 1.0, 100.0] {
            assert_eq!(
                constant
                    .plane_stress_complex_frequency_derivative(plane, omega)
                    .unwrap(),
                ([0.0; 9], [0.0; 9])
            );
        }
    }
    let law = TensorProny::new(
        isotropic(4.0, 3.0),
        None,
        vec![PronyBranch::new(isotropic(2.0, 0.7), 0.2).unwrap()],
        (1.0, 8.0),
    )
    .unwrap();
    let before = law.clone();
    for plane in PLANES {
        for omega in [1.0, 8.0] {
            assert!(
                law.plane_stress_complex_frequency_derivative(plane, omega)
                    .is_ok()
            );
        }
        for omega in [-1.0, 0.0, 9.0] {
            assert!(matches!(
                law.plane_stress_complex_frequency_derivative(plane, omega),
                Err(ViscoError::OutOfBand { .. })
            ));
        }
        for omega in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                law.plane_stress_complex_frequency_derivative(plane, omega),
                Err(ViscoError::Parameters { .. })
            ));
        }
    }
    assert_eq!(law, before);
}
