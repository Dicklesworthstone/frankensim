use super::*;

fn close(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-12 * b.abs().max(1.0), "{a:?} != {b:?}");
}

fn recovered(action: FrequencyAction) -> Vec<C64> {
    assert!(action.envelope.iter().all(|v| v.is_finite()));
    action
        .value
        .into_iter()
        .map(|v| scale_complex(v, action.exponent).unwrap())
        .collect()
}

#[test]
fn diagonal_and_dense_actions_preserve_original_coefficients() {
    let state = [C64::new(0.7, -0.2), C64::new(-0.4, 1.1)];
    for frequency in [0.0, 0.125, 3.0] {
        let diagonal = FrequencyPencil {
            n: 2,
            layout: MatrixLayout::Diagonal,
            scaled_frequency: frequency,
            time_exponent: 0,
            coordinate_exponents: &[0, 0],
            scaled_mass: &[2.0, 3.0],
            scaled_damping: &[0.2, 0.4],
            port: None,
        };
        let dense = FrequencyPencil {
            n: 2,
            layout: MatrixLayout::Dense,
            scaled_frequency: frequency,
            time_exponent: 0,
            coordinate_exponents: &[0, 0],
            scaled_mass: &[2.0, 0.0, 0.0, 3.0],
            scaled_damping: &[0.2, 0.0, 0.0, 0.4],
            port: None,
        };
        for order in [
            FrequencyOrder::First(None),
            FrequencyOrder::Second {
                slope: None,
                curvature: None,
            },
        ] {
            let a = recovered(diagonal.apply(&state, order, || Ok(())).unwrap());
            let b = recovered(dense.apply(&state, order, || Ok(())).unwrap());
            assert_eq!(a, b);
            for i in 0..2 {
                let coefficient = match order {
                    FrequencyOrder::First(_) => C64::new(
                        2.0 * frequency * diagonal.scaled_mass[i],
                        diagonal.scaled_damping[i],
                    ),
                    FrequencyOrder::Second { .. } => C64::from_re(2.0 * diagonal.scaled_mass[i]),
                };
                close(a[i], coefficient * state[i]);
            }
        }
    }
}

#[test]
fn mutual_port_actions_respect_coordinate_and_clock_units() {
    let n = 2;
    let state = [C64::new(1.3, -0.7), C64::new(-0.2, 1.1)];
    let mass = [2.0, 0.25, 0.25, 3.0];
    let damping = [0.2, 0.75, -0.5, 0.3];
    let map = [0.5, -0.25, 0.75, 1.0];
    let impedance = [
        C64::new(2.0, 0.5),
        C64::new(0.4, -0.2),
        C64::new(0.4, -0.2),
        C64::new(1.0, 0.7),
    ];
    let slope = [
        C64::new(0.1, -0.3),
        C64::new(-0.2, 0.4),
        C64::new(0.5, 0.1),
        C64::new(-0.4, 0.2),
    ];
    let curvature = [
        C64::new(0.3, 0.2),
        C64::new(-0.1, 0.6),
        C64::new(0.2, -0.3),
        C64::new(0.5, -0.1),
    ];
    let frequency = 1.75;
    // Independent complete original polynomial contractions at moderate units.
    let mut expected = [vec![C64::ZERO; n], vec![C64::ZERO; n]];
    for i in 0..n {
        for j in 0..n {
            let mut first = C64::new(2.0 * frequency * mass[i * n + j], damping[i * n + j]);
            let mut second = C64::from_re(2.0 * mass[i * n + j]);
            for a in 0..n {
                for b in 0..n {
                    let weight = map[i * n + a] * map[j * n + b];
                    let z = impedance[a * n + b] + slope[a * n + b].scale(frequency);
                    let dz = slope[a * n + b].scale(2.0) + curvature[a * n + b].scale(frequency);
                    first = first + C64::new(-z.im, z.re).scale(weight);
                    second = second + C64::new(-dz.im, dz.re).scale(weight);
                }
            }
            expected[0][i] = expected[0][i] + first * state[j];
            expected[1][i] = expected[1][i] + second * state[j];
        }
    }
    let coordinates = [-80, 90];
    let ports = [30, -20];
    for clock in [-300, 0, 300] {
        let physical_map: Vec<_> = (0..n)
            .flat_map(|i| (0..n).map(move |a| (i, a)))
            .map(|(i, a)| scale(map[i * n + a], coordinates[i] + ports[a]).unwrap())
            .collect();
        let transform = |values: &[C64], time: i32| -> Vec<C64> {
            (0..n)
                .flat_map(|a| (0..n).map(move |b| (a, b)))
                .map(|(a, b)| scale_complex(values[a * n + b], time - ports[a] - ports[b]).unwrap())
                .collect()
        };
        let z = transform(&impedance, clock);
        let dz = transform(&slope, 0);
        let ddz = transform(&curvature, -clock);
        let pencil = FrequencyPencil {
            n,
            layout: MatrixLayout::Dense,
            scaled_frequency: frequency,
            time_exponent: clock,
            coordinate_exponents: &coordinates,
            scaled_mass: &mass,
            scaled_damping: &damping,
            port: Some(FrequencyPort {
                ports: n,
                force_map: &physical_map,
                impedance: &z,
                exponents: &ports,
                minimum_time_exponent: clock + 8,
            }),
        };
        for (order, wanted) in [
            (FrequencyOrder::First(Some(&dz)), &expected[0]),
            (
                FrequencyOrder::Second {
                    slope: Some(&dz),
                    curvature: Some(&ddz),
                },
                &expected[1],
            ),
        ] {
            let actual = recovered(pencil.apply(&state, order, || Ok(())).unwrap());
            for (&a, &b) in actual.iter().zip(wanted) {
                close(a, b);
            }
        }
    }
}

#[test]
fn dc_curvature_ignores_a_huge_coefficient_multiplied_by_exact_zero() {
    let tiny = 2.0_f64.powi(-100);
    let mass = [1.0, tiny, tiny, 1.0];
    let pencil = FrequencyPencil {
        n: 2,
        layout: MatrixLayout::Dense,
        scaled_frequency: 0.0,
        time_exponent: 1,
        coordinate_exponents: &[0, 0],
        scaled_mass: &mass,
        scaled_damping: &[0.0; 4],
        port: Some(FrequencyPort {
            ports: 1,
            force_map: &[1.0, 1.0],
            impedance: &[C64::ONE],
            exponents: &[0],
            minimum_time_exponent: 0,
        }),
    };
    let values = recovered(
        pencil
            .apply(
                &[C64::ZERO, C64::ONE],
                FrequencyOrder::Second {
                    slope: Some(&[C64::ONE]),
                    curvature: Some(&[C64::from_re(2.0_f64.powi(1000))]),
                },
                || Ok(()),
            )
            .unwrap(),
    );
    assert_eq!(values[0].re, 2.0 * tiny);
    close(values[0], C64::new(2.0 * tiny, 2.0));
    close(values[1], C64::new(2.0, 2.0));
}

#[test]
fn reactive_spring_frequency_terms_cancel_before_range_recovery() {
    let stiffness = 2.0_f64.powi(500);
    let impedance = [C64::new(0.0, 2.0 * stiffness)];
    let slope = [C64::new(0.0, -4.0 * stiffness)];
    let curvature = [C64::new(0.0, 16.0 * stiffness)];
    let state = [C64::new(0.75, -0.5)];
    let pencil = FrequencyPencil {
        n: 1,
        layout: MatrixLayout::Diagonal,
        scaled_frequency: 0.5,
        time_exponent: 0,
        coordinate_exponents: &[0],
        scaled_mass: &[3.0],
        scaled_damping: &[0.0],
        port: Some(FrequencyPort {
            ports: 1,
            force_map: &[1.0],
            impedance: &impedance,
            exponents: &[0],
            minimum_time_exponent: 501,
        }),
    };
    let first = recovered(
        pencil
            .apply(&state, FrequencyOrder::First(Some(&slope)), || Ok(()))
            .unwrap(),
    );
    let second = recovered(
        pencil
            .apply(
                &state,
                FrequencyOrder::Second {
                    slope: Some(&slope),
                    curvature: Some(&curvature),
                },
                || Ok(()),
            )
            .unwrap(),
    );
    // -iw Z is constant stiffness; neither derivative has a load contribution.
    close(first[0], state[0].scale(3.0));
    close(second[0], state[0].scale(6.0));
}

#[test]
fn maxwell_curvature_matches_the_independent_original_cubic_polynomial() {
    let (mass, drag, relaxed, branch, tau) = (2.0, 0.3, 4.0, 1.7, 0.4);
    let force = C64::new(0.7, -0.2);
    for frequency in [0.0, 0.3, 2.0, 5.0] {
        let numerator = C64::new(1.0, -frequency * tau);
        let impedance = [C64::from_re(branch * tau) / numerator];
        let slope = [C64::new(0.0, branch * tau * tau) / (numerator * numerator)];
        let curvature =
            [C64::from_re(-2.0 * branch * tau * tau * tau) / (numerator * numerator * numerator)];
        let pencil = FrequencyPencil {
            n: 1,
            layout: MatrixLayout::Diagonal,
            scaled_frequency: frequency,
            time_exponent: 0,
            coordinate_exponents: &[0],
            scaled_mass: &[mass],
            scaled_damping: &[drag],
            port: Some(FrequencyPort {
                ports: 1,
                force_map: &[1.0],
                impedance: &impedance,
                exponents: &[0],
                minimum_time_exponent: 4,
            }),
        };
        let structure = C64::new(relaxed - mass * frequency * frequency, -drag * frequency);
        let operator = structure + C64::new(0.0, -frequency) * impedance[0];
        let q = force / operator;
        let dq = recovered(
            pencil
                .apply(&[q], FrequencyOrder::First(Some(&slope)), || Ok(()))
                .unwrap(),
        )[0] / operator;
        let a2 = recovered(
            pencil
                .apply(
                    &[q],
                    FrequencyOrder::Second {
                        slope: Some(&slope),
                        curvature: Some(&curvature),
                    },
                    || Ok(()),
                )
                .unwrap(),
        )[0];
        let a1 = recovered(
            pencil
                .apply(&[dq], FrequencyOrder::First(Some(&slope)), || Ok(()))
                .unwrap(),
        )[0];
        let ddq = (a2 + a1.scale(2.0)) / operator;
        // q=F*N/P: differentiate the independent original cubic denominator.
        let dn = C64::new(0.0, -tau);
        let ds = C64::new(-2.0 * mass * frequency, -drag);
        let p = structure * numerator + C64::new(0.0, -frequency * branch * tau);
        let dp = ds * numerator + structure * dn + C64::new(0.0, -branch * tau);
        let ddp = numerator.scale(-2.0 * mass) + (ds * dn).scale(2.0);
        let expected = force
            * (-(dn * dp).scale(2.0) / (p * p) - numerator * ddp / (p * p)
                + (numerator * dp * dp).scale(2.0) / (p * p * p));
        close(ddq, expected);
    }
}

#[test]
fn gyroscopic_dc_curvature_retains_cross_motion_without_heat() {
    let mass = [2.0, 0.0, 0.0, 3.0];
    let damping = [0.0, 0.8, -0.8, 0.0];
    let stiffness = [5.0, 7.0];
    let force = [C64::new(1.0, 0.2), C64::new(-0.4, 0.3)];
    let q = [
        force[0].scale(1.0 / stiffness[0]),
        force[1].scale(1.0 / stiffness[1]),
    ];
    let pencil = FrequencyPencil {
        n: 2,
        layout: MatrixLayout::Dense,
        scaled_frequency: 0.0,
        time_exponent: 0,
        coordinate_exponents: &[0, 0],
        scaled_mass: &mass,
        scaled_damping: &damping,
        port: None,
    };
    let first = recovered(
        pencil
            .apply(&q, FrequencyOrder::First(None), || Ok(()))
            .unwrap(),
    );
    let dq = [
        first[0].scale(1.0 / stiffness[0]),
        first[1].scale(1.0 / stiffness[1]),
    ];
    let cross = recovered(
        pencil
            .apply(&dq, FrequencyOrder::First(None), || Ok(()))
            .unwrap(),
    );
    let second = recovered(
        pencil
            .apply(
                &q,
                FrequencyOrder::Second {
                    slope: None,
                    curvature: None,
                },
                || Ok(()),
            )
            .unwrap(),
    );
    for i in 0..2 {
        let ddq = (second[i] + cross[i].scale(2.0)).scale(1.0 / stiffness[i]);
        close(
            ddq,
            q[i].scale(2.0 * (mass[i * 2 + i] + 0.64 / stiffness[1 - i]) / stiffness[i]),
        );
    }
    let work_curvature: f64 = force
        .iter()
        .zip(&dq)
        .map(|(&f, &d)| (f.conj() * C64::new(d.im, -d.re)).re)
        .sum();
    assert!(work_curvature.abs() < 1e-15);
    assert!(dq.iter().any(|v| v.abs() > 0.01));
}

#[test]
fn cancellation_returns_no_partial_action() {
    let pencil = FrequencyPencil {
        n: 2,
        layout: MatrixLayout::Dense,
        scaled_frequency: 1.0,
        time_exponent: 0,
        coordinate_exponents: &[0, 0],
        scaled_mass: &[1.0, 0.0, 0.0, 1.0],
        scaled_damping: &[0.0; 4],
        port: None,
    };
    let mut calls = 0;
    let result = pencil.apply(
        &[C64::ONE; 2],
        FrequencyOrder::Second {
            slope: None,
            curvature: None,
        },
        || {
            calls += 1;
            if calls == 2 {
                Err(HarmonicResponseError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(result, Err(HarmonicResponseError::Cancelled)));
    assert_eq!(calls, 2);
}
