//! G1/G3: exact source jets, original real-state resolvents, and physical scales.
use fs_la::eigen_complex::lu_complex;
use fs_math::c64::C64;
use fs_vfit::mimo::{MimoLaplaceParameterJet, MimoModel, MimoParameterDirection};
use fs_vfit::model::PoleTerm;

fn model() -> MimoModel {
    MimoModel {
        ports: 2,
        poles: vec![
            PoleTerm::Real {
                pole: -3.0,
                residue: 0.0,
            },
            PoleTerm::Pair {
                pole: C64::new(-0.25, 2.0),
                residue: C64::ZERO,
            },
        ],
        residues: vec![
            vec![
                C64::from_re(2.0),
                C64::from_re(-0.5),
                C64::from_re(0.7),
                C64::from_re(3.0),
            ],
            vec![
                C64::new(1.0, 0.2),
                C64::new(0.4, -0.1),
                C64::new(-0.3, 0.15),
                C64::new(0.8, -0.05),
            ],
        ],
        d: vec![1.0, -0.2, 0.3, 1.5],
        e: vec![0.02, -0.01, 0.03, 0.04],
        reciprocal: false,
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale.max(f64::MIN_POSITIVE),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn arrays(jet: &MimoLaplaceParameterJet) -> [&[C64]; 5] {
    [
        &jet.frequency.impedance,
        &jet.frequency.laplace_derivative,
        &jet.frequency.laplace_second_derivative,
        &jet.parameter_derivative,
        &jet.laplace_parameter_derivative,
    ]
}

// Differentiate the original REAL state-space equations. The oracle does not
// use a pole quotient or differentiate samples/solver iterations.
fn resolvent(model: &MimoModel, s: C64, direction: MimoParameterDirection<'_>) -> [Vec<C64>; 5] {
    let source = model.state_space();
    let (n, ports) = (source.n, source.ports);
    let mut matrix = vec![C64::ZERO; n * n];
    for i in 0..n {
        for j in 0..n {
            matrix[i * n + j] = C64::from_re(-source.a[i * n + j]);
        }
        matrix[i * n + i] = matrix[i * n + i] + s;
    }
    let factor = lu_complex(&matrix, n).unwrap();
    let mut da = vec![0.0; n * n];
    let mut dc = vec![0.0; ports * n];
    let mut offset = 0;
    for (group, term) in model.poles.iter().enumerate() {
        let dp = direction.poles.map_or(C64::ZERO, |v| v[group]);
        for input in 0..ports {
            let at = offset + input;
            da[at * n + at] = dp.re;
            if matches!(term, PoleTerm::Pair { .. }) {
                da[at * n + at + ports] = dp.im;
                da[(at + ports) * n + at] = -dp.im;
                da[(at + ports) * n + at + ports] = dp.re;
            }
            for output in 0..ports {
                let dr = direction
                    .residues
                    .map_or(C64::ZERO, |v| v[group][output * ports + input]);
                match term {
                    PoleTerm::Real { .. } => dc[output * n + at] = dr.re,
                    PoleTerm::Pair { .. } => {
                        dc[output * n + at] = 2.0 * dr.re;
                        dc[output * n + at + ports] = 2.0 * dr.im;
                    }
                }
            }
        }
        offset += term.state_dim() * ports;
    }
    let action = |v: &[C64]| -> Vec<C64> {
        (0..n)
            .map(|i| (0..n).fold(C64::ZERO, |a, j| a + v[j].scale(da[i * n + j])))
            .collect()
    };
    let mut result = std::array::from_fn(|_| vec![C64::ZERO; ports * ports]);
    for input in 0..ports {
        let mut x: Vec<_> = (0..n)
            .map(|i| C64::from_re(source.b[i * ports + input]))
            .collect();
        factor.solve(&mut x);
        let mut xs: Vec<_> = x.iter().map(|&z| -z).collect();
        factor.solve(&mut xs);
        let mut xss: Vec<_> = xs.iter().map(|&z| z.scale(-2.0)).collect();
        factor.solve(&mut xss);
        let mut xt = action(&x);
        factor.solve(&mut xt);
        let mut xst = action(&xs);
        for i in 0..n {
            xst[i] = xst[i] - xt[i];
        }
        factor.solve(&mut xst);
        for output in 0..ports {
            let at = output * ports + input;
            let dd = direction.direct.map_or(0.0, |v| v[at]);
            let de = direction.improper.map_or(0.0, |v| v[at]);
            result[0][at] = C64::from_re(source.d[at]) + s.scale(source.e[at]);
            result[1][at] = C64::from_re(source.e[at]);
            result[3][at] = C64::from_re(dd) + s.scale(de);
            result[4][at] = C64::from_re(de);
            for i in 0..n {
                let (c, d) = (source.c[output * n + i], dc[output * n + i]);
                result[0][at] = result[0][at] + x[i].scale(c);
                result[1][at] = result[1][at] + xs[i].scale(c);
                result[2][at] = result[2][at] + xss[i].scale(c);
                result[3][at] = result[3][at] + xt[i].scale(c) + x[i].scale(d);
                result[4][at] = result[4][at] + xst[i].scale(c) + xs[i].scale(d);
            }
        }
    }
    result
}

#[test]
fn full_mutual_jet_matches_original_real_state_resolvent() {
    let source = model();
    let poles = [C64::from_re(0.1), C64::new(0.03, -0.08)];
    let residues = vec![
        vec![
            C64::from_re(0.3),
            C64::from_re(-0.1),
            C64::from_re(0.2),
            C64::from_re(-0.4),
        ],
        vec![
            C64::new(-0.2, 0.03),
            C64::new(0.05, -0.07),
            C64::new(-0.09, -0.01),
            C64::new(0.1, 0.08),
        ],
    ];
    let direct = [0.1, -0.03, 0.02, -0.2];
    let improper = [-0.01, 0.02, -0.03, 0.04];
    let direction = MimoParameterDirection {
        poles: Some(&poles),
        residues: Some(&residues),
        direct: Some(&direct),
        improper: Some(&improper),
    };
    for s in [
        C64::new(0.0, -2.0),
        C64::new(0.0, 2.1),
        C64::new(1.0, 0.7),
        C64::new(0.0, -0.01),
    ] {
        let actual = source.eval_parameter_jet(s, direction).unwrap();
        let expected = resolvent(&source, s, direction);
        for (actual, expected) in arrays(&actual).iter().zip(&expected) {
            for (&a, &e) in actual.iter().zip(expected) {
                close(a, e, 2.0e-12);
            }
        }
        assert_eq!(actual.frequency, source.eval_jet(s).unwrap());
        for (&a, &e) in actual.frequency.impedance.iter().zip(&source.eval(s)) {
            close(a, e, 3.0e-14);
        }
        for (&a, &e) in actual
            .frequency
            .laplace_derivative
            .iter()
            .zip(&source.eval_derivative(s))
        {
            close(a, e, 3.0e-14);
        }
    }
}

#[test]
fn conjugacy_and_physical_time_dilation_keep_all_chain_rules() {
    let source = model();
    let poles: Vec<_> = source
        .poles
        .iter()
        .map(|p| match p {
            PoleTerm::Real { pole, .. } => C64::from_re(*pole),
            PoleTerm::Pair { pole, .. } => *pole,
        })
        .collect();
    let improper: Vec<_> = source.e.iter().map(|&e| -e).collect();
    let direction = MimoParameterDirection {
        poles: Some(&poles),
        residues: Some(&source.residues),
        improper: Some(&improper),
        ..Default::default()
    };
    // Choosing the other stored conjugate changes neither the physical source
    // nor the authored direction. Do not reject an equivalent original chart.
    let mut reversed = source.clone();
    if let PoleTerm::Pair { pole, .. } = &mut reversed.poles[1] {
        *pole = pole.conj();
    }
    for residue in &mut reversed.residues[1] {
        *residue = residue.conj();
    }
    let reversed_poles = [poles[0], poles[1].conj()];
    let reversed_direction = MimoParameterDirection {
        poles: Some(&reversed_poles),
        residues: Some(&reversed.residues),
        improper: Some(&improper),
        ..Default::default()
    };
    for omega in [0.01, 0.7, 2.0, 3.5] {
        let s = C64::new(0.0, -omega);
        let negative = source.eval_parameter_jet(s, direction).unwrap();
        let positive = source.eval_parameter_jet(s.conj(), direction).unwrap();
        assert_eq!(
            negative,
            reversed.eval_parameter_jet(s, reversed_direction).unwrap()
        );
        for (a, b) in arrays(&negative).iter().zip(arrays(&positive)) {
            for (&a, &b) in a.iter().zip(b) {
                assert_eq!(a, b.conj());
            }
        }
        // Authored pθ=p, Rθ=R, Eθ=-E realizes Z(s/(1+θ)).
        for i in 0..4 {
            let zs = negative.frequency.laplace_derivative[i];
            let zss = negative.frequency.laplace_second_derivative[i];
            close(negative.parameter_derivative[i], -(s * zs), 3.0e-14);
            close(
                negative.laplace_parameter_derivative[i],
                -zs - s * zss,
                3.0e-14,
            );
        }
    }
}

#[test]
fn maxwell_strength_and_relaxation_time_differentiate_the_physical_load() {
    let (strength, tau) = (1.375, 0.8);
    let source = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Real {
            pole: -1.0 / tau,
            residue: 0.0,
        }],
        residues: vec![vec![C64::from_re(strength)]],
        d: vec![0.0],
        e: vec![0.0],
        reciprocal: true,
    };
    let dp = [C64::from_re(1.0 / (tau * tau))];
    let dr = vec![vec![C64::from_re(1.0)]];
    for omega in [0.03, 0.5, 2.0, 8.0] {
        let s = C64::new(0.0, -omega);
        let denominator = C64::from_re(1.0) + s.scale(tau);
        let kt = source
            .eval_parameter_jet(
                s,
                MimoParameterDirection {
                    residues: Some(&dr),
                    ..Default::default()
                },
            )
            .unwrap();
        let tt = source
            .eval_parameter_jet(
                s,
                MimoParameterDirection {
                    poles: Some(&dp),
                    ..Default::default()
                },
            )
            .unwrap();
        close(
            tt.frequency.impedance[0],
            C64::from_re(strength * tau) / denominator,
            2.0e-14,
        );
        close(
            kt.parameter_derivative[0],
            C64::from_re(tau) / denominator,
            2.0e-14,
        );
        close(
            kt.laplace_parameter_derivative[0],
            -C64::from_re(tau * tau) / denominator / denominator,
            2.0e-14,
        );
        close(
            tt.parameter_derivative[0],
            C64::from_re(strength) / denominator / denominator,
            2.0e-14,
        );
        close(
            tt.laplace_parameter_derivative[0],
            -C64::from_re(2.0 * strength * tau) / denominator / denominator / denominator,
            2.0e-14,
        );
    }
}

#[test]
fn extreme_pole_scales_and_weak_quadrature_survive_before_rounding() {
    for (distance, residue) in [
        (2.0f64.powi(-600), f64::from_bits(1)),
        (2.0f64.powi(400), 2.0f64.powi(1000)),
    ] {
        let source = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -distance,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(residue)]],
            d: vec![0.0],
            e: vec![0.0],
            reciprocal: true,
        };
        let dp = [C64::from_re(-distance)];
        let jet = source
            .eval_parameter_jet(
                C64::ZERO,
                MimoParameterDirection {
                    poles: Some(&dp),
                    ..Default::default()
                },
            )
            .unwrap();
        let value = residue / distance;
        let first = -value / distance;
        let second = -2.0 * (first / distance);
        assert_eq!(jet.frequency.impedance[0], C64::from_re(value));
        assert_eq!(jet.frequency.laplace_derivative[0], C64::from_re(first));
        assert_eq!(
            jet.frequency.laplace_second_derivative[0],
            C64::from_re(second)
        );
        assert_eq!(jet.parameter_derivative[0], C64::from_re(-value));
        assert_eq!(
            jet.laplace_parameter_derivative[0],
            C64::from_re(-2.0 * first)
        );
    }
    // The original curvature rounds to zero, but multiplying its unrounded
    // weak REAL component by the physical pole tangent gives a finite mixed jet.
    let scale = 2.0f64.powi(1000);
    let source = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Pair {
            pole: C64::new(-1.0, 2.0f64.powi(600)),
            residue: C64::ZERO,
        }],
        residues: vec![vec![C64::from_re(scale)]],
        d: vec![0.0],
        e: vec![0.0],
        reciprocal: true,
    };
    let dp = [C64::from_re(scale)];
    let jet = source
        .eval_parameter_jet(
            C64::ZERO,
            MimoParameterDirection {
                poles: Some(&dp),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(jet.frequency.laplace_second_derivative[0], C64::ZERO);
    close(
        jet.frequency.impedance[0],
        C64::from_re(2.0f64.powi(-199)),
        4.0e-15,
    );
    close(
        jet.parameter_derivative[0],
        C64::from_re(-2.0f64.powi(801)),
        4.0e-15,
    );
    close(
        jet.laplace_parameter_derivative[0],
        C64::from_re(12.0 * 2.0f64.powi(-400)),
        4.0e-15,
    );
}

#[test]
fn signed_mutual_entries_keep_physical_unit_rescaling() {
    let source = model();
    let s = C64::new(0.0, -1.7);
    let dp = [C64::from_re(-0.2), C64::new(0.1, -0.3)];
    let direction = MimoParameterDirection {
        poles: Some(&dp),
        ..Default::default()
    };
    let original = source.eval_parameter_jet(s, direction).unwrap();
    for scale in [2.0f64.powi(-600), 2.0f64.powi(600)] {
        let mut scaled = source.clone();
        for group in &mut scaled.residues {
            for r in group {
                *r = r.scale(scale);
            }
        }
        for c in scaled.d.iter_mut().chain(&mut scaled.e) {
            *c *= scale;
        }
        let actual = scaled.eval_parameter_jet(s, direction).unwrap();
        for (actual, expected) in arrays(&actual).iter().zip(arrays(&original)) {
            for (&a, &e) in actual.iter().zip(expected) {
                close(a, e.scale(scale), 3.0e-14);
            }
        }
    }
}

#[test]
fn incomplete_directions_changed_grouping_and_true_poles_refuse() {
    let source = model();
    let s = C64::new(0.0, -1.0);
    assert!(
        source
            .eval_parameter_jet(
                s,
                MimoParameterDirection {
                    poles: Some(&[C64::ZERO]),
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        source
            .eval_parameter_jet(
                s,
                MimoParameterDirection {
                    poles: Some(&[C64::new(0.0, 1.0), C64::ZERO]),
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        source
            .eval_parameter_jet(
                s,
                MimoParameterDirection {
                    direct: Some(&[f64::NAN; 4]),
                    ..Default::default()
                }
            )
            .is_err()
    );
    let bad = vec![vec![C64::ZERO; 4], vec![C64::ZERO; 3]];
    assert!(
        source
            .eval_parameter_jet(
                s,
                MimoParameterDirection {
                    residues: Some(&bad),
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(source.eval_jet(C64::from_re(-3.0)).is_err());
    assert!(source.eval_jet(C64::new(-0.25, -2.0)).is_err());
    assert!(source.eval_jet(C64::new(f64::NAN, 0.0)).is_err());
    let mut malformed = source.clone();
    malformed.residues[0][0].im = 1.0;
    assert!(malformed.eval_jet(s).is_err());
    let mut divergent = source.clone();
    divergent.poles[0] = PoleTerm::Real {
        pole: -f64::from_bits(1),
        residue: 0.0,
    };
    assert!(divergent.eval_jet(C64::ZERO).is_err());
    let zero = source
        .eval_parameter_jet(s, MimoParameterDirection::default())
        .unwrap();
    assert!(
        zero.parameter_derivative
            .iter()
            .chain(&zero.laplace_parameter_derivative)
            .all(|&v| v == C64::ZERO)
    );
}
