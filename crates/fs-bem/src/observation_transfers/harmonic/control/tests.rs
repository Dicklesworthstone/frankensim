use super::*;

const MEDIUM: Medium = Medium {
    density: 2.0,
    sound_speed: 3.0,
};
const LIMITS: HarmonicControlLimits = HarmonicControlLimits {
    maximum_quadrature_points: 20,
    maximum_products: 1_000_000,
    maximum_workspace_bytes: 1_000_000,
};
const WORK: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 1_000_000,
};

fn owner(columns: &[Vec<[C64; 4]>]) -> PreparedHarmonicObservation {
    let omega = 7.0;
    let sources: Vec<_> = (0..columns.len())
        .map(|i| SourceId::StructuralMode(format!("drive-{i}")))
        .collect();
    let receivers: Vec<_> = (0..columns[0].len())
        .map(|i| ReceiverSpec {
            id: format!("receiver-{i}"),
            position_m: [i as f64, 0.2, -0.3],
        })
        .collect();
    let mut samples = Vec::new();
    for column in columns {
        for field in column {
            samples.push(field[0]);
            samples.extend(
                field[1..]
                    .iter()
                    .map(|&u| C64::new(0.0, omega * MEDIUM.density) * u),
            );
        }
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        MEDIUM,
        &sources,
        &receivers,
        ContentHash([23; 32]),
        samples,
    )
    .unwrap()
}

fn general() -> PreparedHarmonicObservation {
    let columns: Vec<Vec<[C64; 4]>> = (0..3)
        .map(|s| {
            (0..2)
                .map(|r| {
                    core::array::from_fn(|c| {
                        C64::new(
                            0.13 * (1 + s + 2 * r + c) as f64,
                            0.17 * (s as f64 - r as f64 - c as f64) + 0.09,
                        )
                    })
                })
                .collect()
        })
        .collect();
    owner(&columns)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:e}, expected {expected:e}"
    );
}
fn complex_close(actual: C64, expected: C64, tolerance: f64) {
    close(actual.re, expected.re, tolerance);
    close(actual.im, expected.im, tolerance);
}

#[test]
fn hermitian_source_pullback_matches_complex_duality_without_reciprocity() {
    let field = general();
    let coefficients = [C64::new(0.3, -0.2), C64::new(-0.7, 0.4), C64::new(0.1, 0.9)];
    let dual = [
        HarmonicFieldCovector {
            pressure: C64::new(0.4, -0.2),
            particle_velocity: [C64::new(0.7, 0.3), C64::new(-0.2, 0.6), C64::new(0.9, -0.1)],
        },
        HarmonicFieldCovector {
            pressure: C64::new(-0.5, 0.8),
            particle_velocity: [C64::new(-0.3, 0.2), C64::new(0.1, -0.8), C64::new(0.2, 0.5)],
        },
    ];
    let mut observed = [HarmonicField::default(); 2];
    field
        .evaluate_coherent(&coefficients, &mut observed, || false)
        .unwrap();
    let mut gradient = [C64::ZERO; 3];
    field
        .source_pullback(&dual, WORK, &mut gradient, || false)
        .unwrap();
    let mut left = C64::ZERO;
    for (actual, dual) in observed.iter().zip(dual) {
        left = left + actual.pressure_pa.conj() * dual.pressure;
        for c in 0..3 {
            left = left + actual.particle_velocity_m_s[c].conj() * dual.particle_velocity[c];
        }
    }
    let right = coefficients
        .iter()
        .zip(gradient)
        .fold(C64::ZERO, |sum, (&a, g)| sum + a.conj() * g);
    complex_close(left, right, 3e-14);
    // A transpose instead of a Hermitian adjoint is demonstrably different
    // on these complex, asymmetric retained transfers.
    let wrong = (0..2).fold(C64::ZERO, |sum, r| {
        sum + field.transfers[r][0] * dual[r].pressure
    });
    assert!((wrong - gradient[0]).abs() > 0.1);
}

fn time_energy(
    field: &PreparedHarmonicObservation,
    coefficients: &[C64],
    targets: &[HarmonicField],
    q: &[HarmonicVolumeQuadrature],
) -> f64 {
    let mut value = 0.0;
    for sample in q {
        let r = sample.receiver;
        let mut residual = [C64::ZERO; 4];
        for (s, &a) in coefficients.iter().enumerate() {
            for c in 0..4 {
                residual[c] = residual[c] + field.transfers[s * targets.len() + r][c] * a;
            }
        }
        let target = targets[r];
        let t = [
            target.pressure_pa,
            target.particle_velocity_m_s[0],
            target.particle_velocity_m_s[1],
            target.particle_velocity_m_s[2],
        ];
        for c in 0..4 {
            residual[c] = residual[c] - t[c];
        }
        for n in 0..512 {
            let angle = core::f64::consts::TAU * (n as f64 + 0.5) / 512.0;
            let v = residual.map(|z| z.re * angle.cos() + z.im * angle.sin());
            value += sample.volume_m3
                * 0.5
                * (v[0] * v[0] / (MEDIUM.density * MEDIUM.sound_speed * MEDIUM.sound_speed)
                    + MEDIUM.density * (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]))
                / 512.0;
        }
    }
    value
}

#[test]
fn energy_gradient_and_hessian_match_physical_time_energy_and_source_perturbations() {
    let field = general();
    let coefficients = [C64::new(0.3, -0.2), C64::new(-0.7, 0.4), C64::new(0.1, 0.9)];
    let targets = [
        HarmonicField {
            pressure_pa: C64::new(0.2, -0.1),
            particle_velocity_m_s: [C64::new(0.03, -0.02); 3],
        },
        HarmonicField {
            pressure_pa: C64::new(-0.3, 0.4),
            particle_velocity_m_s: [C64::new(-0.01, 0.07); 3],
        },
    ];
    let q = [
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: 0.4,
        },
        HarmonicVolumeQuadrature {
            receiver: 1,
            volume_m3: 0.9,
        },
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: 0.7,
        },
    ];
    let mut objective =
        PreparedHarmonicEnergyControl::new(&field, &q, Some(&targets), LIMITS, || false).unwrap();
    let result = objective.value_gradient(&coefficients, || false).unwrap();
    close(
        result.mean_residual_energy_j,
        time_energy(&field, &coefficients, &targets, &q),
        2e-13,
    );
    let gradient = result.source_gradient.to_vec();
    // An independent centered perturbation of the real-time energy oracle.
    // The objective is exactly quadratic, so this has no truncation term.
    let step = 0.001;
    for source in 0..coefficients.len() {
        for imaginary in [false, true] {
            let delta = if imaginary {
                C64::new(0.0, step)
            } else {
                C64::from_re(step)
            };
            let mut plus = coefficients;
            let mut minus = coefficients;
            plus[source] = plus[source] + delta;
            minus[source] = minus[source] - delta;
            let slope = (time_energy(&field, &plus, &targets, &q)
                - time_energy(&field, &minus, &targets, &q))
                / (2.0 * step);
            close(
                if imaginary {
                    gradient[source].im
                } else {
                    gradient[source].re
                },
                slope,
                5e-10,
            );
        }
    }
    let direction = [C64::new(-0.2, 0.6), C64::new(0.7, 0.1), C64::new(0.3, -0.4)];
    let hessian = objective
        .hessian_vector(&direction, || false)
        .unwrap()
        .to_vec();
    let plus = core::array::from_fn::<_, 3, _>(|i| coefficients[i] + direction[i].scale(step));
    let minus = core::array::from_fn::<_, 3, _>(|i| coefficients[i] - direction[i].scale(step));
    let gp = objective
        .value_gradient(&plus, || false)
        .unwrap()
        .source_gradient
        .to_vec();
    let gm = objective
        .value_gradient(&minus, || false)
        .unwrap()
        .source_gradient
        .to_vec();
    for i in 0..3 {
        complex_close(hessian[i], (gp[i] - gm[i]).scale(0.5 / step), 3e-11);
    }
    let curvature = direction
        .iter()
        .zip(&hessian)
        .map(|(&d, &h)| (d.conj() * h).re)
        .sum::<f64>();
    assert!(curvature >= 0.0);
    let mut quiet = PreparedHarmonicEnergyControl::new(&field, &q, None, LIMITS, || false).unwrap();
    assert_eq!(hessian, quiet.hessian_vector(&direction, || false).unwrap());
}

#[test]
fn standing_wave_velocity_energy_and_phase_invariant_control_are_retained() {
    let right = [C64::from_re(6.0), C64::ONE, C64::ZERO, C64::ZERO];
    let left = [C64::from_re(6.0), -C64::ONE, C64::ZERO, C64::ZERO];
    let field = owner(&[vec![right], vec![left]]);
    let q = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 1.0,
    }];
    let mut objective =
        PreparedHarmonicEnergyControl::new(&field, &q, None, LIMITS, || false).unwrap();
    for coefficients in [[C64::ONE, C64::ONE], [C64::ONE, -C64::ONE]] {
        let result = objective.value_gradient(&coefficients, || false).unwrap();
        close(result.mean_residual_energy_j, 2.0, 1e-14);
        for i in 0..2 {
            complex_close(result.source_gradient[i], coefficients[i].scale(2.0), 1e-14);
        }
    }
    let coefficients = [C64::new(0.3, -0.4), C64::new(-0.2, 0.7)];
    let before = objective.value_gradient(&coefficients, || false).unwrap();
    let value = before.mean_residual_energy_j;
    let gradient = before.source_gradient.to_vec();
    let phase = C64::new(0.6, 0.8);
    let after = objective
        .value_gradient(&coefficients.map(|a| a * phase), || false)
        .unwrap();
    close(after.mean_residual_energy_j, value, 2e-14);
    for i in 0..2 {
        complex_close(after.source_gradient[i], gradient[i] * phase, 2e-14);
    }
    // One native gradient step with the exact Hessian removes this field.
    let cancelled = core::array::from_fn::<_, 2, _>(|i| coefficients[i] - gradient[i].scale(0.5));
    assert!(
        objective
            .value_gradient(&cancelled, || false)
            .unwrap()
            .mean_residual_energy_j
            < 1e-28
    );
}

#[test]
fn weak_target_mismatch_is_formed_before_rounding_the_large_forward_field() {
    let field = owner(&[
        vec![[C64::from_re(1e16), C64::ZERO, C64::ZERO, C64::ZERO]],
        vec![[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]],
    ]);
    let targets = [HarmonicField {
        pressure_pa: C64::from_re(1e16),
        particle_velocity_m_s: [C64::ZERO; 3],
    }];
    let q = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 1.0,
    }];
    let mut objective =
        PreparedHarmonicEnergyControl::new(&field, &q, Some(&targets), LIMITS, || false).unwrap();
    let result = objective.value_gradient(&[C64::ONE; 2], || false).unwrap();
    close(result.mean_residual_energy_j, 1.0 / 72.0, 1e-16);
    complex_close(result.source_gradient[0], C64::from_re(1e16 / 36.0), 1e-15);
    complex_close(result.source_gradient[1], C64::from_re(1.0 / 36.0), 1e-16);
}

#[test]
fn physical_admission_and_cancellation_preserve_reusable_control() {
    let mut field = general();
    let q = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 0.4,
    }];
    let coefficients = [C64::ONE; 3];
    let mut control =
        PreparedHarmonicEnergyControl::new(&field, &q, None, LIMITS, || false).unwrap();
    assert_eq!(control.work_products(), 2 * (8 * 3 + 16));
    let bytes = control.workspace_bytes();
    assert!(
        PreparedHarmonicEnergyControl::new(
            &field,
            &q,
            None,
            HarmonicControlLimits {
                maximum_workspace_bytes: bytes - 1,
                ..LIMITS
            },
            || false
        )
        .is_err()
    );
    assert!(
        PreparedHarmonicEnergyControl::new(
            &field,
            &q,
            None,
            HarmonicControlLimits {
                maximum_products: control.work_products() - 1,
                ..LIMITS
            },
            || false
        )
        .is_err()
    );
    let reference = control.value_gradient(&coefficients, || false).unwrap();
    let value = reference.mean_residual_energy_j;
    let gradient = reference.source_gradient.to_vec();
    let mut polls = 0;
    assert!(matches!(
        control.value_gradient(&coefficients, || {
            polls += 1;
            polls == 5
        }),
        Err(TransferError::Cancelled)
    ));
    let replay = control.value_gradient(&coefficients, || false).unwrap();
    assert_eq!(replay.mean_residual_energy_j, value);
    assert_eq!(replay.source_gradient, gradient);
    assert!(control.value_gradient(&[C64::ZERO; 2], || false).is_err());
    let bad_q = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: -1.0,
    }];
    assert!(PreparedHarmonicEnergyControl::new(&field, &bad_q, None, LIMITS, || false).is_err());
    let mut dual = [HarmonicFieldCovector::default(); 2];
    dual[1].pressure = C64::new(f64::NAN, 0.0);
    let sentinel = C64::new(17.0, 9.0);
    let mut output = [sentinel; 3];
    assert!(
        field
            .source_pullback(&dual, WORK, &mut output, || false)
            .is_err()
    );
    assert_eq!(output, [sentinel; 3]);
    field.radiation_momentum = false;
    assert!(PreparedHarmonicEnergyControl::new(&field, &q, None, LIMITS, || false).is_err());
    // A general field pullback needs no nondispersive constitutive energy law.
    field
        .source_pullback(
            &[HarmonicFieldCovector::default(); 2],
            WORK,
            &mut output,
            || false,
        )
        .unwrap();
    assert_eq!(output, [C64::ZERO; 3]);
}
