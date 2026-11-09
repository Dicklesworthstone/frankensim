//! Manufactured spectral laws; no empirical isotope calibration is implied.

use super::*;

fn moments(
    shape: AnalyticFissionSpectrum,
    theta: f64,
    cutoff: f64,
    edges: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let mut counts = vec![PositiveIntegral::new(); edges.len() - 1];
    let mut energies = counts.clone();
    integrate(
        shape,
        theta,
        cutoff,
        edges,
        &mut counts,
        &mut energies,
        &mut 0,
        DecayEmissionLimits::default(),
        &mut || Ok(()),
    )
    .unwrap();
    (
        counts.into_iter().map(|v| v.value().unwrap()).collect(),
        energies.into_iter().map(|v| v.value().unwrap()).collect(),
    )
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_finite_maxwell_and_evaporation_match_independent_gamma_references() {
    // 100-digit mpmath incomplete-gamma references using the exact binary64
    // input endpoints. Neither reference uses the production midpoint series.
    let cases = [
        (
            AnalyticFissionSpectrum::Maxwellian,
            [
                0.022414506931348845,
                0.40525138164730152,
                0.31099536363223778,
                0.26133874778911185,
            ],
            [
                0.0013294338445573454,
                0.22499142828288213,
                0.44966995517112978,
                0.82238892548380658,
            ],
        ),
        (
            AnalyticFissionSpectrum::Evaporation,
            [
                0.0046811779370895152,
                0.25969196746458140,
                0.32991779332413397,
                0.40570906127419512,
            ],
            [
                0.00030946068495631730,
                0.16037357844342704,
                0.48628722495551580,
                1.3484874745380296,
            ],
        ),
    ];
    for (shape, expected_counts, expected_energies) in cases {
        let (counts, energies) = moments(shape, 1.0, 10.0, &[0.0, 0.1, 1.0, 2.0, 10.0]);
        for (&actual, expected) in counts.iter().zip(expected_counts) {
            close(actual, expected, 3e-14);
        }
        for (&actual, expected) in energies.iter().zip(expected_energies) {
            close(actual, expected, 3e-14);
        }
        close(counts.iter().sum(), 1.0, 3e-14);
    }
}

#[test]
fn g1_tiny_cutoff_keeps_the_normalized_power_law_and_energy_scale() {
    // R/theta=1e-300 erases a raw gamma normalizer. The normalized limit on
    // E/R is s*u^(s-1), with a correction far below binary64 precision.
    let cutoff = 1e-200;
    for (shape, s) in [
        (AnalyticFissionSpectrum::Maxwellian, 1.5_f64),
        (AnalyticFissionSpectrum::Evaporation, 2.0_f64),
    ] {
        let fractions = [0.0_f64, 0.1, 0.5, 1.0];
        let edges: Vec<_> = fractions.iter().map(|value| value * cutoff).collect();
        let (counts, energies) = moments(shape, 1e100, cutoff, &edges);
        for (index, pair) in fractions.windows(2).enumerate() {
            close(counts[index], pair[1].powf(s) - pair[0].powf(s), 3e-14);
            close(
                energies[index],
                cutoff * s / (s + 1.0) * (pair[1].powf(s + 1.0) - pair[0].powf(s + 1.0)),
                3e-14,
            );
        }
    }
}

#[test]
fn g1_g3_narrow_bins_and_far_tails_retain_their_positive_moments() {
    // Tail references are positive 100-digit anchored integrals. Differences
    // of ordinary CDFs in this region would round to zero.
    for (shape, narrow_count, narrow_energy, tail_count, tail_energy) in [
        (
            AnalyticFissionSpectrum::Maxwellian,
            3.7760216468126753e-13,
            3.7760216468143925e-13,
            6.8224705935711757e-269,
            4.0934823561458824e-266,
        ),
        (
            AnalyticFissionSpectrum::Evaporation,
            3.3475157731569058e-13,
            3.3475157731584281e-13,
            1.4810244842090722e-267,
            8.8861469052613301e-265,
        ),
    ] {
        let (counts, energies) =
            moments(shape, 1.0, 10.0, &[0.0, 1.0, 1.0 + 2.0_f64.powi(-40), 10.0]);
        close(counts[1], narrow_count, 4e-14);
        close(energies[1], narrow_energy, 4e-14);
        let (counts, energies) = moments(
            shape,
            1.0,
            650.0,
            &[0.0, 600.0, 600.0 + 2.0_f64.powi(-30), 650.0],
        );
        close(counts[1], tail_count, 5e-13);
        close(energies[1], tail_energy, 5e-13);
        assert!(counts[2] > 0.0 && energies[2] > 0.0);
    }
    // The exact midpoint of adjacent floats is not itself representable.
    // Preserve its exposure before exp, including a non-exact theta quotient.
    // These references integrate the exact binary64 endpoints at 110 digits.
    for (theta, left, cutoff, expected) in [
        (
            1.0,
            600.0_f64,
            650.0,
            [
                [8.3282111776673599e-273, 4.9969267066004164e-270],
                [1.8078912169158869e-271, 1.0847347301495322e-268],
            ],
        ),
        (
            core::f64::consts::PI,
            1884.9555921538758_f64,
            2042.0352248333654,
            [
                [5.3019039041558987e-273, 9.9938534132011289e-270],
                [1.1509392949783757e-271, 2.1694694602991287e-268],
            ],
        ),
    ] {
        let right = f64::from_bits(left.to_bits() + 1);
        for (shape, expected) in [
            AnalyticFissionSpectrum::Maxwellian,
            AnalyticFissionSpectrum::Evaporation,
        ]
        .into_iter()
        .zip(expected)
        {
            let (counts, energies) = moments(shape, theta, cutoff, &[0.0, left, right, cutoff]);
            close(counts[1], expected[0], 3e-15);
            close(energies[1], expected[1], 3e-15);
        }
    }
}

#[test]
fn g3_group_refinement_and_energy_unit_rescaling_preserve_both_moments() {
    for shape in [
        AnalyticFissionSpectrum::Maxwellian,
        AnalyticFissionSpectrum::Evaporation,
    ] {
        let coarse = moments(shape, 1.0, 4.0, &[0.0, 4.0]);
        let edges = [0.0, 0.001, 0.1, 0.5, 1.5, 2.5, 4.0, 5.0];
        let fine = moments(shape, 1.0, 4.0, &edges);
        close(fine.0.iter().sum(), coarse.0[0], 3e-14);
        close(fine.1.iter().sum(), coarse.1[0], 3e-14);
        assert_eq!(fine.0[6], 0.0);
        assert_eq!(fine.1[6], 0.0);
        for factor in [1e-120, 1e120] {
            let scaled_edges: Vec<_> = edges.iter().map(|edge| edge * factor).collect();
            let scaled = moments(shape, factor, 4.0 * factor, &scaled_edges);
            for group in 0..edges.len() - 1 {
                close(scaled.0[group], fine.0[group], 5e-14);
                close(scaled.1[group], fine.1[group] * factor, 5e-14);
            }
        }
    }
}

#[test]
fn g0_g4_support_arithmetic_work_and_final_cancellation_refuse() {
    let shape = AnalyticFissionSpectrum::Maxwellian;
    let run = |theta,
               cutoff,
               edges: &[f64],
               limits,
               poll: &mut dyn FnMut() -> Result<(), FissionError>| {
        let mut counts = vec![PositiveIntegral::new(); edges.len() - 1];
        let mut energies = counts.clone();
        integrate(
            shape,
            theta,
            cutoff,
            edges,
            &mut counts,
            &mut energies,
            &mut 0,
            limits,
            &mut || poll(),
        )
    };
    for (theta, cutoff, edges) in [
        (0.0, 10.0, vec![0.0, 10.0]),
        (1.0, 10.0, vec![0.1, 10.0]),
        (1.0, 10.0, vec![0.0, 9.0]),
        (1.0, 10.0, vec![0.0, 5.0, 5.0, 10.0]),
        (1e300, 1e-300, vec![0.0, 1e-300]),
        (1.0, 1000.0, vec![0.0, 800.0, 1000.0]),
    ] {
        assert!(
            run(
                theta,
                cutoff,
                &edges,
                DecayEmissionLimits::default(),
                &mut || Ok(())
            )
            .is_err()
        );
    }
    for limits in [
        DecayEmissionLimits {
            max_groups: 1,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_group_values: 1,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_work: 1,
            ..DecayEmissionLimits::default()
        },
    ] {
        assert!(run(1.0, 4.0, &[0.0, 0.5, 1.5, 4.0], limits, &mut || Ok(())).is_err());
    }
    let mut polls = 0;
    run(
        1.0,
        4.0,
        &[0.0, 0.5, 1.5, 4.0],
        DecayEmissionLimits::default(),
        &mut || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut remaining = polls - 1;
    assert!(matches!(
        run(
            1.0,
            4.0,
            &[0.0, 0.5, 1.5, 4.0],
            DecayEmissionLimits::default(),
            &mut || {
                if remaining == 0 {
                    Err(FissionError::Cancelled)
                } else {
                    remaining -= 1;
                    Ok(())
                }
            }
        ),
        Err(EndfError::Physics(FissionError::Cancelled))
    ));
}
