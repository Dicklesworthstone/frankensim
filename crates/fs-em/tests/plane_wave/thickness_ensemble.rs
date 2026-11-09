//! Discrete thickness ensembles preserve coherent Maxwell realizations and
//! separate reflected field fluctuations from physical material absorption.
use super::*;
use fs_em::plane_wave::layers::{
    DielectricLayer, LayerThicknessRealization, PlanarDielectricLayers, ThicknessEnsembleBudget,
};

fn budget(rows: usize, layers: usize) -> ThicknessEnsembleBudget {
    ThicknessEnsembleBudget {
        maximum_realizations: rows,
        maximum_layer_realizations: rows * layers.max(1),
    }
}

#[test]
fn one_realization_and_identical_rows_retain_the_original_native_response() {
    let outside = material(1.0, 1.0);
    let film = material(9.0, 1.0);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: 0.7,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &films,
    )
    .unwrap();
    let geometry = [0.7];
    let single = [LayerThicknessRealization {
        thicknesses_m: &geometry,
        probability: 1.0,
    }];
    let duplicate = [
        LayerThicknessRealization {
            thicknesses_m: &geometry,
            probability: 0.2,
        },
        LayerThicknessRealization {
            thicknesses_m: &geometry,
            probability: 0.8,
        },
    ];
    for rows in [&single[..], &duplicate[..]] {
        let ensemble = stack
            .prepare_thickness_ensemble(rows, budget(rows.len(), 1), || false)
            .unwrap();
        assert!(std::ptr::eq(ensemble.source(), &stack));
        assert!(std::ptr::eq(ensemble.realizations(), rows));
        assert_eq!(ensemble.layer_realization_work(), rows.len());
        for polarization in [Polarization::Te, Polarization::Tm] {
            let original = stack.response(0.6, polarization).unwrap();
            let result = ensemble.evaluate(0.6, polarization, || false).unwrap();
            assert_eq!(
                result.mean_reflected_tangential_e,
                original.reflected_tangential_e
            );
            assert_eq!(
                result.mean_reflected_power_fraction,
                original.reflected_power_fraction
            );
            assert_eq!(
                result.mean_transmitted_power_fraction,
                original.transmitted_power_fraction
            );
            assert_eq!(
                result.mean_absorbed_power_fraction,
                original.absorbed_power_fraction
            );
            near(
                result.coherent_reflected_power_fraction,
                original.reflected_power_fraction,
                2e-16,
            );
            assert_eq!(result.fluctuating_reflected_power_fraction, 0.0);
            assert_eq!(result.reflected_power_variance, 0.0);
            assert_eq!(result.evaluated_realizations, rows.len());
        }
    }
    assert_eq!(stack.layer(0).unwrap().thickness_m, 0.7);
}

#[test]
fn uniform_phase_refinement_recovers_incoherent_power_and_retains_coherent_field() {
    let outside = material(1.0, 1.0);
    let film = material(9.0, 1.0);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: 0.0,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &films,
    )
    .unwrap();
    // A symmetric lossless film has single-interface R0=1/4. Summing
    // independent intensity round trips gives R=2R0/(1+R0)=0.4 and T=0.6.
    // Averaging the coherent Airy amplitude leaves r01=-1/2, so the
    // fluctuating reflected component is 0.4-0.25=0.15, with zero heat.
    let mut previous_error = f64::INFINITY;
    for count in [4, 8, 16, 64] {
        let geometry: Vec<_> = (0..count)
            .map(|i| [std::f64::consts::PI * (i as f64 + 0.5) / (3.0 * count as f64)])
            .collect();
        let rows: Vec<_> = geometry
            .iter()
            .map(|d| LayerThicknessRealization {
                thicknesses_m: d,
                probability: 1.0 / count as f64,
            })
            .collect();
        let ensemble = stack
            .prepare_thickness_ensemble(&rows, budget(count, 1), || false)
            .unwrap();
        let result = ensemble.evaluate(1.0, Polarization::Te, || false).unwrap();
        let error = (result.mean_reflected_power_fraction - 0.4).abs();
        assert!(error <= previous_error + 2e-15);
        previous_error = error;
        near(
            result.mean_reflected_power_fraction + result.mean_transmitted_power_fraction,
            1.0,
            2e-14,
        );
        assert_eq!(result.mean_absorbed_power_fraction, 0.0);
        assert!(result.maximum_absolute_power_balance_error < 2e-14);
        if count == 64 {
            near(result.mean_reflected_power_fraction, 0.4, 2e-14);
            near(result.mean_transmitted_power_fraction, 0.6, 2e-14);
            complex_near(
                result.mean_reflected_tangential_e,
                C64::from_re(-0.5),
                2e-14,
            );
            near(result.coherent_reflected_power_fraction, 0.25, 2e-14);
            near(result.fluctuating_reflected_power_fraction, 0.15, 2e-14);
            assert!(result.reflected_power_variance > 0.01);
        }
    }
}

#[test]
fn joint_rows_preserve_cross_layer_thickness_correlations() {
    let outside = material(1.0, 1.0);
    let film = material(9.0, 1.0);
    let films = [
        DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.2,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.8,
        },
    ];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &films,
    )
    .unwrap();
    let geometry = [[0.2, 0.2], [0.2, 0.8], [0.8, 0.2], [0.8, 0.8]];
    let correlated = [
        LayerThicknessRealization {
            thicknesses_m: &geometry[1],
            probability: 0.5,
        },
        LayerThicknessRealization {
            thicknesses_m: &geometry[2],
            probability: 0.5,
        },
    ];
    let independent: Vec<_> = geometry
        .iter()
        .map(|d| LayerThicknessRealization {
            thicknesses_m: d,
            probability: 0.25,
        })
        .collect();
    let fixed_total = stack
        .prepare_thickness_ensemble(&correlated, budget(2, 2), || false)
        .unwrap()
        .evaluate(1.0, Polarization::Te, || false)
        .unwrap();
    let independent = stack
        .prepare_thickness_ensemble(&independent, budget(4, 2), || false)
        .unwrap()
        .evaluate(1.0, Polarization::Te, || false)
        .unwrap();
    // Identical materials with complementary thicknesses always total 1 m;
    // removing the correlation gives exactly the same marginal distributions
    // but introduces real phase/fringe fluctuations.
    let nominal = stack.response(1.0, Polarization::Te).unwrap();
    near(
        fixed_total.mean_reflected_power_fraction,
        nominal.reflected_power_fraction,
        2e-14,
    );
    complex_near(
        fixed_total.mean_reflected_tangential_e,
        nominal.reflected_tangential_e,
        2e-14,
    );
    assert!(fixed_total.fluctuating_reflected_power_fraction < 1e-28);
    assert!(fixed_total.reflected_power_variance < 1e-28);
    assert!(
        (independent.mean_reflected_power_fraction - fixed_total.mean_reflected_power_fraction)
            .abs()
            > 0.01
    );
    assert!(independent.fluctuating_reflected_power_fraction > 0.01);
    assert!(independent.reflected_power_variance > 0.001);
}

// Independent one-film Airy sum in the fixed tangential-E convention.
// This is intentionally an unscaled analytic oracle for moderate fixtures;
// production retains its bounded E/H solve and positive absorption modes.
fn airy(
    epsilon: C64,
    thickness: f64,
    cosine: f64,
    polarization: Polarization,
) -> (C64, f64, f64, f64) {
    let q = (epsilon - C64::from_re((1.0 - cosine) * (1.0 + cosine))).sqrt();
    let incident_y = match polarization {
        Polarization::Te => cosine,
        Polarization::Tm => 1.0 / cosine,
    };
    let film_y = match polarization {
        Polarization::Te => q,
        Polarization::Tm => epsilon / q,
    };
    let incident_y = C64::from_re(incident_y);
    let r01 = (incident_y - film_y) / (incident_y + film_y);
    let r12 = -r01;
    let phase = q.re * thickness;
    let p = C64::new(phase.cos(), -phase.sin()).scale((q.im * thickness).exp());
    let denominator = C64::ONE + r01 * r12 * p * p;
    let r = (r01 + r12 * p * p) / denominator;
    let t = (C64::ONE + r01) * (C64::ONE + r12) * p / denominator;
    let reflected = r.norm_sq();
    let transmitted = t.norm_sq(); // the original exterior media are equal
    (r, reflected, transmitted, 1.0 - reflected - transmitted)
}

#[test]
fn lossy_realizations_match_airy_and_keep_heat_distinct_from_reflected_fluctuations() {
    let outside = material(1.0, 1.0);
    let state = source(
        &[
            (EPS_NAME, EPS, 2.0),
            (MU_NAME, MU, 1.0),
            ("delta", EPS, 2.0),
            ("tau", Time::DIMS, 1.0),
        ],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    let poles = [DebyePoleProperties {
        permittivity_increment: "delta".into(),
        relaxation_time: "tau".into(),
    }];
    let dielectric =
        ResolvedDebyeCapacitor::try_from_material(state, EPS_NAME, &poles, 1.0, 1.0).unwrap();
    let films = [DielectricLayer {
        response: DielectricResponse::Debye(&dielectric),
        thickness_m: 0.7,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &films,
    )
    .unwrap();
    let epsilon = C64::new(3.0, -1.0);
    complex_near(stack.layer_permittivity_f_m(0).unwrap(), epsilon, 2e-15);
    let geometry = [[0.0], [0.2], [0.7], [2.0]];
    let probabilities = [0.1, 0.2, 0.3, 0.4];
    let rows: Vec<_> = geometry
        .iter()
        .zip(probabilities)
        .map(|(d, probability)| LayerThicknessRealization {
            thicknesses_m: d,
            probability,
        })
        .collect();
    let ensemble = stack
        .prepare_thickness_ensemble(&rows, budget(4, 1), || false)
        .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let result = ensemble.evaluate(0.6, polarization, || false).unwrap();
        let samples: Vec<_> = rows
            .iter()
            .map(|row| airy(epsilon, row.thicknesses_m[0], 0.6, polarization))
            .collect();
        let mut mean_r = C64::ZERO;
        let (mut reflected, mut transmitted, mut absorbed) = (0.0, 0.0, 0.0);
        for (sample, probability) in samples.iter().zip(probabilities) {
            mean_r = mean_r + sample.0.scale(probability);
            reflected += probability * sample.1;
            transmitted += probability * sample.2;
            absorbed += probability * sample.3;
        }
        let mut field_variance = 0.0;
        let mut power_variance = 0.0;
        for (sample, probability) in samples.iter().zip(probabilities) {
            field_variance += probability * (sample.0 - mean_r).norm_sq();
            power_variance += probability * (sample.1 - reflected).powi(2);
        }
        complex_near(result.mean_reflected_tangential_e, mean_r, 2e-14);
        near(result.mean_reflected_power_fraction, reflected, 2e-14);
        near(result.mean_transmitted_power_fraction, transmitted, 2e-14);
        near(result.mean_absorbed_power_fraction, absorbed, 2e-14);
        near(
            result.fluctuating_reflected_power_fraction,
            field_variance,
            2e-14,
        );
        near(result.reflected_power_variance, power_variance, 2e-14);
        near(
            result.coherent_reflected_power_fraction + result.fluctuating_reflected_power_fraction,
            reflected,
            2e-14,
        );
        assert!(result.mean_absorbed_power_fraction > 0.2);
        assert!(result.fluctuating_reflected_power_fraction > 1e-4);
        assert!(result.maximum_absolute_power_balance_error < 2e-14);
    }
    assert_eq!(stack.layer(0).unwrap().thickness_m, 0.7);
}

#[test]
fn admission_budget_cancellation_and_replay_are_bounded_and_do_not_change_the_source() {
    let outside = material(1.0, 1.0);
    let film = material(4.0, 1.0);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: 0.5,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &films,
    )
    .unwrap();
    let geometry = [[0.0], [0.2], [0.7]];
    let rows = [
        LayerThicknessRealization {
            thicknesses_m: &geometry[0],
            probability: 0.0,
        },
        LayerThicknessRealization {
            thicknesses_m: &geometry[1],
            probability: 0.25,
        },
        LayerThicknessRealization {
            thicknesses_m: &geometry[2],
            probability: 0.75,
        },
    ];
    assert!(
        stack
            .prepare_thickness_ensemble(&[], budget(3, 1), || false)
            .is_err()
    );
    assert!(
        stack
            .prepare_thickness_ensemble(&rows, budget(2, 1), || false)
            .is_err()
    );
    assert!(
        stack
            .prepare_thickness_ensemble(
                &rows,
                ThicknessEnsembleBudget {
                    maximum_realizations: 3,
                    maximum_layer_realizations: 2,
                },
                || false
            )
            .is_err()
    );
    assert!(
        stack
            .prepare_thickness_ensemble(&rows, budget(3, 1), || true)
            .is_err()
    );
    for invalid in [f64::NAN, f64::INFINITY, -0.1] {
        let thickness = [invalid];
        let bad_row = [LayerThicknessRealization {
            thicknesses_m: &thickness,
            probability: 1.0,
        }];
        assert!(
            stack
                .prepare_thickness_ensemble(&bad_row, budget(1, 1), || false)
                .is_err()
        );
    }
    for probability in [f64::NAN, f64::INFINITY, -0.1, 0.0, 0.9, 1.1] {
        let bad_row = [LayerThicknessRealization {
            thicknesses_m: &geometry[0],
            probability,
        }];
        assert!(
            stack
                .prepare_thickness_ensemble(&bad_row, budget(1, 1), || false)
                .is_err()
        );
    }
    let missing = [LayerThicknessRealization {
        thicknesses_m: &[],
        probability: 1.0,
    }];
    assert!(
        stack
            .prepare_thickness_ensemble(&missing, budget(1, 1), || false)
            .is_err()
    );
    // Even a dark row's authored geometry must be valid.
    let invalid = [f64::NAN];
    let dark = [
        LayerThicknessRealization {
            thicknesses_m: &invalid,
            probability: 0.0,
        },
        LayerThicknessRealization {
            thicknesses_m: &geometry[1],
            probability: 1.0,
        },
    ];
    assert!(
        stack
            .prepare_thickness_ensemble(&dark, budget(2, 1), || false)
            .is_err()
    );
    let ensemble = stack
        .prepare_thickness_ensemble(&rows, budget(3, 1), || false)
        .unwrap();
    let expected = ensemble.evaluate(0.8, Polarization::Tm, || false).unwrap();
    assert_eq!(expected.evaluated_realizations, 2);
    let mut calls = 0;
    assert!(
        ensemble
            .evaluate(0.8, Polarization::Tm, || {
                calls += 1;
                calls == 4 // one illuminated realization has already been evaluated
            })
            .is_err()
    );
    assert_eq!(
        ensemble.evaluate(0.8, Polarization::Tm, || false).unwrap(),
        expected
    );
    for cosine in [0.0, -0.1, 1.1, f64::NAN] {
        assert!(
            ensemble
                .evaluate(cosine, Polarization::Te, || false)
                .is_err()
        );
    }
    assert_eq!(stack.layer(0).unwrap().thickness_m, 0.5);
}
