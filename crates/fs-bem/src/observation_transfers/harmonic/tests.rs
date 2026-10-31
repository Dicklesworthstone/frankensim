use super::*;
use crate::helmholtz::bc::{
    Admittance,
    incident::{self, IncidentField, Monopole},
};
use crate::near_field::Options;

const BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 1_000_000,
};

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-9),
        "actual {actual:?}, expected {expected:?}"
    );
}

// Construct independently declared pressure/velocity transfers. The production
// bank accepts p/grad and must perform the Euler conversion with the right sign.
fn analytic(fields: &[[C64; 4]], sources: usize, receivers: usize) -> PreparedHarmonicObservation {
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.0,
    };
    let omega = 1200.0;
    let ids: Vec<_> = (0..sources)
        .map(|i| SourceId::StructuralMode(format!("coordinate-{i}")))
        .collect();
    let positions: Vec<_> = (0..receivers)
        .map(|i| ReceiverSpec {
            id: format!("receiver-{i}"),
            position_m: [i as f64, 0.0, 0.0],
        })
        .collect();
    let mut samples = Vec::new();
    for f in fields {
        samples.push(f[0]);
        samples.extend(
            f[1..]
                .iter()
                .map(|&v| C64::new(0.0, omega * medium.density) * v),
        );
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        medium,
        &ids,
        &positions,
        ContentHash([3; 32]),
        samples,
    )
    .unwrap()
}

/// G1: two counterpropagating waves at a quarter standing-wave phase retain
/// reactive flow while their mean transported energy cancels. A pressure-only
/// plane-wave power approximation would incorrectly report positive flux.
#[test]
fn standing_wave_statistics_retain_reactive_velocity_and_microphone_projection() {
    let z = 1.2 * 343.0;
    let right = phase(core::f64::consts::FRAC_PI_4);
    let left = right.conj();
    let owner = analytic(
        &[
            [right, right.scale(1.0 / z), C64::ZERO, C64::ZERO],
            [left, left.scale(-1.0 / z), C64::ZERO, C64::ZERO],
        ],
        2,
        1,
    );
    let mut output = [HarmonicFieldStatistics::default()];
    owner
        .evaluate_coherence_modes(&[C64::ONE; 2], 1, BUDGET, &mut output, || false)
        .unwrap();
    close(C64::from_re(output[0].mean_square_pressure_pa2), C64::ONE);
    close(output[0].complex_intensity_w_m2[0], C64::new(0.0, -1.0 / z));
    close(
        output[0].velocity_cross_m2_s2[0][0],
        C64::from_re(1.0 / (z * z)),
    );
    let front = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    let mut microphone = owner.microphones(&[front]).unwrap();
    close(
        microphone
            .evaluate_coherence_modes(&[C64::ONE; 2], 1, BUDGET, || false)
            .unwrap()[0],
        C64::from_re(0.5),
    );
    owner
        .evaluate_independent_sources(&[1.0, 1.0], BUDGET, &mut output, || false)
        .unwrap();
    close(C64::from_re(output[0].mean_square_pressure_pa2), C64::ONE);
    assert!(output[0].complex_intensity_w_m2[0].abs() < 1e-18);
}

/// G0/G3: factor propagation agrees with an independently assembled dense
/// H C H* contraction, including phase-bearing cross-receiver terms. Scaling
/// every source by a common phase leaves the quadratic observation unchanged.
#[test]
fn microphone_cross_spectra_match_dense_covariance_and_beamformer_power() {
    let fields = [
        [
            C64::new(2.0, -1.0),
            C64::new(0.01, 0.03),
            C64::ZERO,
            C64::ZERO,
        ],
        [
            C64::new(-0.5, 3.0),
            C64::new(-0.02, 0.01),
            C64::ZERO,
            C64::ZERO,
        ],
        [
            C64::new(1.0, 2.0),
            C64::new(0.02, -0.01),
            C64::ZERO,
            C64::ZERO,
        ],
        [
            C64::new(4.0, 0.5),
            C64::new(0.03, 0.02),
            C64::ZERO,
            C64::ZERO,
        ],
    ];
    let owner = analytic(&fields, 2, 2);
    let factors = [
        C64::new(0.4, 0.2),
        C64::new(-0.1, 0.3),
        C64::new(0.8, -0.1),
        C64::new(0.2, -0.4),
        C64::new(-0.3, 0.7),
        C64::new(0.5, 0.6),
    ];
    let patterns = [
        FirstOrder::new(1.0, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.25, [1.0, 0.0, 0.0]).unwrap(),
    ];
    let mut microphones = owner.microphones(&patterns).unwrap();
    let actual = microphones
        .evaluate_coherence_modes(&factors, 3, BUDGET, || false)
        .unwrap()
        .to_vec();
    let mut covariance = [[C64::ZERO; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for mode in 0..3 {
                covariance[i][j] =
                    covariance[i][j] + factors[2 * mode + i] * factors[2 * mode + j].conj();
            }
        }
    }
    let h = [
        [fields[0][0], fields[2][0]],
        [
            fields[1][0].scale(0.25) - fields[1][1].scale(0.75 * 1.2 * 343.0),
            fields[3][0].scale(0.25) - fields[3][1].scale(0.75 * 1.2 * 343.0),
        ],
    ];
    for i in 0..2 {
        for j in 0..2 {
            let mut expected = C64::ZERO;
            for s in 0..2 {
                for t in 0..2 {
                    expected = expected + h[i][s] * covariance[s][t] * h[j][t].conj();
                }
            }
            close(actual[2 * i + j], expected.scale(0.5));
            assert_eq!(actual[2 * i + j], actual[2 * j + i].conj());
        }
    }
    assert!(
        actual[1].im.abs() > 0.1,
        "spatial phase must not be discarded"
    );
    let weights = [C64::new(0.2, 0.6), C64::new(-0.7, 0.1)];
    let mut matrix_power = C64::ZERO;
    for i in 0..2 {
        for j in 0..2 {
            matrix_power = matrix_power + weights[i].conj() * actual[2 * i + j] * weights[j];
        }
    }
    let mut direct_power = 0.0;
    for mode in 0..3 {
        let mut output = C64::ZERO;
        for r in 0..2 {
            for s in 0..2 {
                output = output + weights[r].conj() * h[r][s] * factors[2 * mode + s];
            }
        }
        direct_power += 0.5 * output.norm_sq();
    }
    close(matrix_power, C64::from_re(direct_power));
    let rotated: Vec<_> = factors.iter().map(|&z| z * phase(0.73)).collect();
    for (&a, &b) in microphones
        .evaluate_coherence_modes(&rotated, 3, BUDGET, || false)
        .unwrap()
        .iter()
        .zip(&actual)
    {
        close(a, b);
    }
}

/// G0/G3: incoherent sources are the diagonal factor limit but retain spatial
/// microphone coherence; truly coherent opposite sources can cancel exactly.
#[test]
fn independent_fast_paths_and_rank_one_cancellation() {
    let p = C64::new(1.3, -0.7);
    let field = [
        p,
        C64::new(0.003, 0.004),
        C64::new(-0.001, 0.002),
        C64::ZERO,
    ];
    let owner = analytic(&[field, field], 2, 1);
    let mut coherent = [HarmonicFieldStatistics::default()];
    owner
        .evaluate_coherence_modes(&[C64::ONE, -C64::ONE], 1, BUDGET, &mut coherent, || false)
        .unwrap();
    assert_eq!(coherent[0], HarmonicFieldStatistics::default());
    let factors = [C64::from_re(2.0), C64::ZERO, C64::ZERO, C64::from_re(3.0)];
    owner
        .evaluate_coherence_modes(&factors, 2, BUDGET, &mut coherent, || false)
        .unwrap();
    let mut independent = [HarmonicFieldStatistics::default()];
    owner
        .evaluate_independent_sources(&[4.0, 9.0], BUDGET, &mut independent, || false)
        .unwrap();
    close(
        C64::from_re(coherent[0].mean_square_pressure_pa2),
        C64::from_re(independent[0].mean_square_pressure_pa2),
    );
    for i in 0..3 {
        close(
            coherent[0].complex_intensity_w_m2[i],
            independent[0].complex_intensity_w_m2[i],
        );
        for j in 0..3 {
            close(
                coherent[0].velocity_cross_m2_s2[i][j],
                independent[0].velocity_cross_m2_s2[i][j],
            );
        }
    }
    close(
        C64::from_re(independent[0].mean_square_pressure_pa2),
        C64::from_re(6.5 * p.norm_sq()),
    );
    let mut microphones = owner
        .microphones(&[FirstOrder::new(0.3, [0.0, 1.0, 0.0]).unwrap()])
        .unwrap();
    let factor_power = microphones
        .evaluate_coherence_modes(&factors, 2, BUDGET, || false)
        .unwrap()[0];
    close(
        microphones
            .evaluate_independent_sources(&[4.0, 9.0], BUDGET, || false)
            .unwrap()[0],
        factor_power,
    );
}

/// G0: signed compensated amplitude reduction preserves a weak real signal
/// between two large opposite emitters; statistics are taken after reduction.
#[test]
fn quiet_residual_is_retained_before_squaring() {
    let owner = analytic(&[[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]; 3], 3, 1);
    let drives = [C64::from_re(1e16), C64::ONE, C64::from_re(-1e16)];
    let mut fields = [HarmonicField::default()];
    owner
        .evaluate_coherent(&drives, &mut fields, || false)
        .unwrap();
    assert_eq!(fields[0].pressure_pa, C64::ONE);
    let mut stats = [HarmonicFieldStatistics::default()];
    owner
        .evaluate_coherence_modes(&drives, 1, BUDGET, &mut stats, || false)
        .unwrap();
    close(
        C64::from_re(stats[0].mean_square_pressure_pa2),
        C64::from_re(0.5),
    );
}

/// G1: the bank is connected to the actual rigid-obstacle solution. Compare
/// independently solving a simultaneously driven source array, including its
/// incident field, against combining the prepared unit-source total fields.
#[test]
fn prepared_diffraction_uses_coherent_total_field_before_statistics() {
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.0,
    };
    let omega = core::f64::consts::TAU * 700.0;
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let sources = [
        PointSourceSpec {
            id: "a".into(),
            position_m: [0.02, 0.01, 0.09],
        },
        PointSourceSpec {
            id: "b".into(),
            position_m: [-0.04, 0.03, 0.07],
        },
    ];
    let receivers = [
        ReceiverSpec {
            id: "shadow".into(),
            position_m: [-0.03, 0.02, -0.08],
        },
        ReceiverSpec {
            id: "side".into(),
            position_m: [0.07, -0.03, 0.04],
        },
    ];
    let owner = ScatteringObservation::new(
        &surface,
        medium,
        &sources,
        &receivers,
        omega,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap();
    let prepared = owner.prepare_harmonic(omega, || false).unwrap();
    assert_eq!(prepared.sources(), owner.sources());
    assert_eq!(prepared.producer_hash(), owner.operator_hash());
    assert!(matches!(
        prepared.producer_uncertainty(),
        ProducerUncertainty::Unresolved { .. }
    ));
    let amplitudes = [C64::new(0.002, 0.001), C64::new(-0.003, 0.002)];
    let monopoles = std::array::from_fn::<_, 2, _>(|i| Monopole {
        position_m: sources[i].position_m,
        volume_velocity_m3_s: amplitudes[i],
    });
    let incident = IncidentField::Monopoles(&monopoles);
    let solved = incident::solve_scattering(
        &surface,
        omega / medium.sound_speed,
        medium,
        Admittance::Rigid,
        incident,
        1e-3,
        helmholtz::Formulation::BurtonMiller,
    )
    .unwrap();
    let points = receivers.map(|r| r.position_m);
    let total = solved
        .total_pressure_velocity_at(&points, 1e-3, Options::default())
        .unwrap();
    let mut fields = [HarmonicField::default(); 2];
    prepared
        .evaluate_coherent(&amplitudes, &mut fields, || false)
        .unwrap();
    let mut statistics = [HarmonicFieldStatistics::default(); 2];
    prepared
        .evaluate_coherence_modes(&amplitudes, 1, BUDGET, &mut statistics, || false)
        .unwrap();
    for r in 0..2 {
        close(fields[r].pressure_pa, total.scalar.pressure[r]);
        for axis in 0..3 {
            close(
                fields[r].particle_velocity_m_s[axis],
                total.particle_velocity_m_s[r][axis],
            );
            close(
                statistics[r].complex_intensity_w_m2[axis],
                (total.scalar.pressure[r] * total.particle_velocity_m_s[r][axis].conj()).scale(0.5),
            );
        }
    }
    assert!(owner.prepare_harmonic(0.0, || false).is_err());
    assert!(owner.prepare_harmonic(omega, || true).is_err());
}

/// G5: fail before unbounded factor work, refuse invalid statistics, preserve
/// cancellation, and explicitly admit the zero-rank silent process.
#[test]
fn shapes_budget_nonfinite_and_cancellation_are_explicit() {
    let owner = analytic(&[[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]], 1, 1);
    let mut out = [HarmonicFieldStatistics::default()];
    assert!(matches!(
        owner.evaluate_coherence_modes(
            &[C64::ONE],
            1,
            HarmonicWorkBudget {
                maximum_products: 1
            },
            &mut out,
            || false
        ),
        Err(TransferError::Cap { .. })
    ));
    assert!(
        owner
            .evaluate_coherence_modes(&[C64::ONE], 2, BUDGET, &mut out, || false)
            .is_err()
    );
    assert!(
        owner
            .evaluate_coherence_modes(&[C64::new(f64::NAN, 0.0)], 1, BUDGET, &mut out, || false)
            .is_err()
    );
    assert!(
        owner
            .evaluate_independent_sources(&[-1.0], BUDGET, &mut out, || false)
            .is_err()
    );
    assert!(matches!(
        owner.evaluate_coherence_modes(&[C64::ONE], 1, BUDGET, &mut out, || true),
        Err(TransferError::Cancelled)
    ));
    owner
        .evaluate_coherence_modes(&[], 0, BUDGET, &mut out, || false)
        .unwrap();
    assert_eq!(out[0], HarmonicFieldStatistics::default());
    let mut microphones = owner
        .microphones(&[FirstOrder::new(1.0, [1.0, 0.0, 0.0]).unwrap()])
        .unwrap();
    assert!(
        microphones
            .evaluate_independent_sources(&[f64::INFINITY], BUDGET, || false)
            .is_err()
    );
    assert_eq!(
        microphones
            .evaluate_coherence_modes(&[], 0, BUDGET, || false)
            .unwrap(),
        &[C64::ZERO]
    );
    assert!(
        admit_frequency(
            1.0,
            Medium {
                density: 0.0,
                sound_speed: 343.0
            }
        )
        .is_err()
    );
}
