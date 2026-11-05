use fs_la::eigen_complex::lu_complex;
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::model::PoleTerm;

#[test]
fn conjugate_mutual_model_derivative_matches_independent_state_resolvent() {
    let model = MimoModel {
        ports: 2,
        poles: vec![
            PoleTerm::Real {
                pole: -2.0,
                residue: 0.0,
            },
            PoleTerm::Pair {
                pole: C64::new(-0.04, 7.0),
                residue: C64::ZERO,
            },
        ],
        residues: vec![
            vec![
                C64::from_re(0.8),
                C64::from_re(-0.2),
                C64::from_re(0.3),
                C64::from_re(0.6),
            ],
            vec![
                C64::new(0.5, 0.1),
                C64::new(-0.2, 0.3),
                C64::new(0.1, -0.4),
                C64::new(0.6, 0.2),
            ],
        ],
        d: vec![1.0, 0.2, -0.1, 0.8],
        e: vec![0.1, -0.02, 0.03, 0.05],
        reciprocal: false,
    };
    let state = model.state_space();
    let (n, p) = (state.n, state.ports);
    for s in [
        C64::ZERO,
        C64::new(0.0, 7.0),
        C64::new(0.0, -7.0),
        C64::new(0.13, 5.0),
        C64::new(1.2, -0.8),
    ] {
        let actual = model.eval_derivative(s);
        let conjugate = model.eval_derivative(s.conj());
        let mut resolvent = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                resolvent[i * n + j] = C64::from_re(-state.a[i * n + j]);
            }
            resolvent[i * n + i] = resolvent[i * n + i] + s;
        }
        let factor = lu_complex(&resolvent, n).unwrap();
        for column in 0..p {
            let mut x: Vec<_> = (0..n)
                .map(|i| C64::from_re(state.b[i * p + column]))
                .collect();
            factor.solve(&mut x);
            factor.solve(&mut x);
            for row in 0..p {
                let mut expected = C64::from_re(state.e[row * p + column]);
                for (i, &value) in x.iter().enumerate() {
                    expected = expected - value.scale(state.c[row * n + i]);
                }
                let j = row * p + column;
                assert!(
                    (actual[j] - expected).abs() <= 3e-11 * expected.abs().max(1.0),
                    "s={s:?} entry={j}: {:?} != {expected:?}",
                    actual[j]
                );
                assert!(
                    (actual[j].conj() - conjugate[j]).abs() <= 1e-14 * actual[j].abs().max(1.0)
                );
            }
        }
    }
}

#[test]
fn original_pole_distances_are_not_squared_or_regularized() {
    for rate in [1e-200, 1e200] {
        let model = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -rate,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(rate)]],
            d: vec![0.0],
            e: vec![0.0],
            reciprocal: true,
        };
        let derivative = model.eval_derivative(C64::ZERO)[0];
        let expected = -1.0 / rate;
        assert!(derivative.re.is_finite());
        assert!((derivative.re - expected).abs() <= 4e-15 * expected.abs());
        assert_eq!(derivative.im, 0.0);
        let singular = model.eval_derivative(C64::from_re(-rate))[0];
        assert!(!singular.re.is_finite() || !singular.im.is_finite());
    }
    let direct = MimoModel {
        ports: 1,
        poles: vec![],
        residues: vec![],
        d: vec![8.0],
        e: vec![0.25],
        reciprocal: true,
    };
    assert_eq!(
        direct.eval_derivative(C64::new(0.7, -4.0)),
        vec![C64::from_re(0.25)]
    );
}
