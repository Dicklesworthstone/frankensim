use super::*;
use crate::SpeciesId;
use crate::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_qty::Temperature;

fn model(n: usize, pairs: &[(usize, usize, f64)]) -> MulticomponentRegularSolutionV1 {
    let mut chi = vec![0.0; n * n];
    for &(i, j, value) in pairs {
        chi[i * n + j] = value;
        chi[j * n + i] = value;
    }
    MulticomponentRegularSolutionV1::new(
        (0..n)
            .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
            .collect(),
        chi,
        vec![0.0; n * n],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}
fn state(n: usize, pairs: &[(usize, usize, f64)], x: &[f64]) -> MulticomponentSolutionStateV1 {
    model(n, pairs)
        .evaluate(Temperature::new(300.0), x, || false)
        .unwrap()
}
fn assess(state: &MulticomponentSolutionStateV1) -> MulticomponentLocalStabilityV1 {
    state
        .local_stability(MulticomponentStabilityControlV1::default(), || false)
        .unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * (1.0 + b.abs()),
        "actual={a:e}, expected={b:e}"
    );
}
fn check_direction(state: &MulticomponentSolutionStateV1, result: &MulticomponentLocalStabilityV1) {
    near(result.direction().iter().sum(), 0.0, 2e-14);
    near(
        sum(result
            .direction()
            .iter()
            .zip(state.mole_fractions())
            .filter(|(_, x)| **x > 0.0)
            .map(|(v, x)| {
                let q = v / sqrt(*x);
                q * q
            })),
        1.0,
        2e-13,
    );
    near(
        state.reduced_curvature(result.direction()).unwrap(),
        result.smallest_eigenvalue().unwrap(),
        2e-11 + result.numerical_uncertainty(),
    );
}

#[test]
fn g1_binary_entropy_metric_matches_the_independent_spinodal_formula() {
    for chi in [-3.0, 0.0, 1.5, 2.0, 3.0, 8.0] {
        for x in [1e-20, 0.1, 0.2, 0.5, 0.8] {
            let s = state(2, &[(0, 1, chi)], &[x, 1.0 - x]);
            let result = assess(&s);
            let exact = 1.0 - 2.0 * chi * x * (1.0 - x);
            near(result.smallest_eigenvalue().unwrap(), exact, 2e-13);
            check_direction(&s, &result);
            let expected = if exact < -1e-9 {
                MulticomponentLocalStabilityKindV1::Spinodal
            } else if exact > 1e-9 {
                MulticomponentLocalStabilityKindV1::LocallyStable
            } else {
                MulticomponentLocalStabilityKindV1::Unresolved
            };
            assert_eq!(result.kind(), expected);
        }
    }
}

#[test]
fn g1_coupled_spinodal_is_detected_when_every_pairwise_exchange_is_positive() {
    let s = state(
        3,
        &[(0, 1, 2.9), (0, 2, 2.9), (1, 2, -4.0)],
        &[1.0 / 3.0; 3],
    );
    for d in [[1.0, -1.0, 0.0], [1.0, 0.0, -1.0], [0.0, 1.0, -1.0]] {
        assert!(s.reduced_curvature(&d).unwrap() > 0.0);
    }
    let result = assess(&s);
    assert_eq!(result.kind(), MulticomponentLocalStabilityKindV1::Spinodal);
    assert_eq!(result.negative_modes(), 1);
    near(result.eigenvalues()[0], -13.2 / 18.0, 2e-13);
    near(result.eigenvalues()[1], 14.0 / 6.0, 2e-13);
    check_direction(&s, &result);
    let d = result.direction();
    near(d[1], d[2], 2e-13);
    near(d[0], -2.0 * d[1], 2e-13);
}

#[test]
fn g1_metastable_ternary_has_positive_local_curvature_and_negative_tangent_plane() {
    let s = state(3, &[(0, 1, 2.9), (0, 2, 2.9), (1, 2, 2.9)], &[1.0 / 3.0; 3]);
    let result = assess(&s);
    assert_eq!(
        result.kind(),
        MulticomponentLocalStabilityKindV1::LocallyStable
    );
    for &value in result.eigenvalues() {
        near(value, 1.0 - 2.9 / 3.0, 3e-13);
    }
    // Independent scalar common-tangent equation ln(a/b)+chi*(b-a)=0.
    let (mut low, mut high): (f64, f64) = (2.0 / 3.0, 1.0 - 1e-12);
    for _ in 0..80 {
        let a = 0.5 * (low + high);
        let b = (1.0 - a) / 2.0;
        if (a / b).ln() + 2.9 * (b - a) > 0.0 {
            high = a;
        } else {
            low = a;
        }
    }
    let a = 0.5 * (low + high);
    let b = (1.0 - a) / 2.0;
    let tpd = s.tangent_plane_distance(&[a, b, b], || false).unwrap();
    near(tpd, -0.014813560085012667, 2e-13);
    assert!(tpd < 0.0);
    for chi in [3.0, 4.0] {
        let s = state(3, &[(0, 1, chi), (0, 2, chi), (1, 2, chi)], &[1.0 / 3.0; 3]);
        let result = assess(&s);
        assert_eq!(result.negative_modes(), if chi > 3.0 { 2 } else { 0 });
        assert_eq!(
            result.kind(),
            if chi > 3.0 {
                MulticomponentLocalStabilityKindV1::Spinodal
            } else {
                MulticomponentLocalStabilityKindV1::Unresolved
            }
        );
    }
}

#[test]
fn g0_ideal_maximum_species_and_subnormal_traces_need_no_inverse_fraction() {
    let mut x = vec![1.0 / 127.0; 128];
    x[0] = 1e-320;
    let s = state(128, &[], &x);
    let result = assess(&s);
    assert_eq!(result.eigenvalues().len(), 127);
    assert_eq!(
        result.kind(),
        MulticomponentLocalStabilityKindV1::LocallyStable
    );
    assert_eq!(result.completed_sweeps(), 0);
    for &v in result.eigenvalues() {
        near(v, 1.0, 1e-15);
    }
    assert!(result.direction()[0] > 0.0);
    check_direction(&s, &result);
}

#[test]
fn g0_closed_support_excludes_absent_species_and_pure_states_have_no_mode() {
    let pairs = [(0, 1, 3.0), (0, 2, 8.0), (1, 2, 9.0)];
    let s = state(3, &pairs, &[0.4, 0.6, 0.0]);
    let result = assess(&s);
    assert_eq!(result.eigenvalues().len(), 1);
    assert_eq!(result.direction()[2], 0.0);
    near(result.eigenvalues()[0], -0.44, 2e-13);
    check_direction(&s, &result);
    assert!(matches!(
        s.tangent_plane_distance(&[0.3, 0.3, 0.4], || false),
        Err(MulticomponentSolutionErrorV1::Invalid(_))
    ));
    let pure = state(3, &pairs, &[0.0, 1.0, 0.0]);
    let result = assess(&pure);
    assert_eq!(
        result.kind(),
        MulticomponentLocalStabilityKindV1::NoCompositionDirections
    );
    assert!(result.eigenvalues().is_empty());
    assert_eq!(result.direction(), [0.0; 3]);
    assert_eq!(
        pure.tangent_plane_distance(&[0.0, 1.0, 0.0], || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn g1_ideal_tangent_plane_matches_relative_entropy_and_self_distance_is_zero() {
    let x = [0.1, 0.3, 0.6];
    let y = [0.4, 0.0, 0.6];
    let s = state(3, &[], &x);
    let exact = sum((0..3)
        .filter(|&i| y[i] > 0.0)
        .map(|i| y[i] * (y[i] / x[i]).ln()));
    near(
        s.tangent_plane_distance(&y, || false).unwrap(),
        exact,
        2e-13,
    );
    assert_eq!(s.tangent_plane_distance(&x, || false).unwrap(), 0.0);
    assert!(matches!(
        s.tangent_plane_distance(&y, || true),
        Err(MulticomponentSolutionErrorV1::Cancelled)
    ));
}

#[test]
fn g3_species_relabeling_preserves_the_spectrum_and_curvature() {
    let x = [0.1, 0.2, 0.3, 0.4];
    let pairs = [
        (0, 1, 2.0),
        (0, 2, -1.0),
        (0, 3, 0.5),
        (1, 2, 3.0),
        (1, 3, 4.0),
        (2, 3, 1.0),
    ];
    let s = state(4, &pairs, &x);
    let original = assess(&s);
    let p = [3, 0, 2, 1];
    let law = model(4, &pairs);
    let mut permuted_pairs = Vec::new();
    for i in 0..4 {
        for j in i + 1..4 {
            permuted_pairs.push((i, j, law.chi_constants()[p[i] * 4 + p[j]]));
        }
    }
    let permuted = state(4, &permuted_pairs, &p.map(|i| x[i]));
    let result = assess(&permuted);
    for (a, b) in original.eigenvalues().iter().zip(result.eigenvalues()) {
        near(*a, *b, 2e-12);
    }
    assert_eq!(original.kind(), result.kind());
    check_direction(&permuted, &result);
}

#[test]
fn g4_bounded_work_memory_and_cancellation_return_no_partial_assessment() {
    let s = state(
        4,
        &[
            (0, 1, 2.0),
            (0, 2, -1.0),
            (0, 3, 0.5),
            (1, 2, 3.0),
            (1, 3, 4.0),
            (2, 3, 1.0),
        ],
        &[0.1, 0.2, 0.3, 0.4],
    );
    let unchanged = s.clone();
    let control = MulticomponentStabilityControlV1::default();
    assert!(matches!(
        s.local_stability(
            MulticomponentStabilityControlV1 {
                max_workspace_bytes: 0,
                ..control
            },
            || false
        ),
        Err(MulticomponentStabilityErrorV1::WorkspaceBudget { .. })
    ));
    assert!(matches!(
        s.local_stability(
            MulticomponentStabilityControlV1 {
                max_sweeps: 1,
                relative_eigensolve_tolerance: 1e-16,
                ..control
            },
            || false
        ),
        Err(MulticomponentStabilityErrorV1::WorkBudget { .. })
    ));
    assert!(matches!(
        s.local_stability(
            MulticomponentStabilityControlV1 {
                max_sweeps: 0,
                ..control
            },
            || false
        ),
        Err(MulticomponentStabilityErrorV1::InvalidControl)
    ));
    let mut polls = 0;
    s.local_stability(control, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in 0..polls {
        let mut current = 0;
        let result = s.local_stability(control, || {
            let cancel = current == stop;
            current += 1;
            cancel
        });
        assert!(matches!(
            result,
            Err(MulticomponentStabilityErrorV1::Cancelled)
        ));
        assert_eq!(s, unchanged);
    }
}
