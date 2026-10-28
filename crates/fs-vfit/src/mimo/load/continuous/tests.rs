use super::*;
use crate::mimo::MimoModel;
use crate::model::PoleTerm;
use crate::passivity::CertificateClass;
use fs_math::c64::C64;

const DT: f64 = 1.0 / 48_000.0;

fn source() -> MimoModel {
    MimoModel {
        ports: 2,
        poles: vec![
            PoleTerm::Real {
                pole: -400.0,
                residue: 0.0,
            },
            PoleTerm::Pair {
                pole: C64::new(-80.0, 1200.0),
                residue: C64::ZERO,
            },
        ],
        residues: vec![
            vec![40.0, 10.0, 10.0, 30.0]
                .into_iter()
                .map(C64::from_re)
                .collect(),
            vec![300.0, 75.0, 75.0, 200.0]
                .into_iter()
                .map(|value| C64::new(value, value * 80.0 / 1200.0))
                .collect(),
        ],
        d: vec![2.0, 0.2, 0.2, 1.5],
        e: vec![0.0; 4],
        reciprocal: true,
    }
}

fn prepare(source: &MimoModel) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        source,
        DT,
        super::super::MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 33,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}

/// G0/G3: original partial fractions independently check the full real
/// realization over frequencies, phases and mutual channels. Nonzero accepted
/// history does not replace the source by a sampled or reduced operator.
#[test]
fn continuous_original_multiport_realization_matches_complex_partial_fractions() {
    let model = source();
    let mut load = prepare(&model);
    load.step_mean_into(&[0.7, -0.3], &mut [0.0; 2]).unwrap();
    let accepted = load.pole_states().to_vec();
    let view = load.continuous_operator();
    assert_eq!(view.dimension(), 6);
    assert_eq!(view.ports(), 2);
    let original = view.state_space(1 << 20, &mut || false).unwrap();
    assert_eq!(original, model.state_space());
    assert_eq!(original.d, model.d);
    assert_ne!(original.d, load.mean_direct_term());
    let input = [C64::new(0.7, 0.2), C64::new(-0.4, 0.9)];
    for omega in [1.0, 80.0, 600.0, 1200.0, 12000.0, 30000.0] {
        view.admit_frequency(omega).unwrap();
        let s = C64::new(0.0, -omega);
        let rational = model.eval(s);
        let realized = original.eval(s).unwrap();
        for row in 0..2 {
            let expected = (0..2).fold(C64::ZERO, |sum, j| sum + rational[row * 2 + j] * input[j]);
            let actual = (0..2).fold(C64::ZERO, |sum, j| sum + realized[row * 2 + j] * input[j]);
            assert!(
                (actual - expected).abs() < 2e-11 * expected.abs() + 1e-13,
                "omega={omega}, row={row}: {actual:?} != {expected:?}"
            );
        }
    }
    // The original real chart uses NEGATIVE pair imaginary components.
    assert!(accepted[2].im != 0.0 && accepted[3].im != 0.0);
    assert_eq!(
        view.real_coordinates(),
        vec![
            accepted[0].re,
            accepted[1].re,
            accepted[2].re,
            accepted[3].re,
            -accepted[2].im,
            -accepted[3].im,
        ]
    );
    // Adopt an endpoint produced by the actual original ZOH owner, then
    // continue the original held-input path without resetting any pole.
    let endpoint_input = [-0.2, 0.9];
    let mut evolved = load.clone();
    evolved
        .step_mean_into(&endpoint_input, &mut [0.0; 2])
        .unwrap();
    let endpoint = evolved.continuous_operator().real_coordinates();
    let mut candidate = view
        .candidate_from_coordinates(&endpoint, &endpoint_input)
        .unwrap();
    assert_eq!(candidate.pole_states(), evolved.pole_states());
    assert_eq!(candidate.source_model(), load.source_model());
    assert_eq!(candidate.limits(), load.limits());
    assert_eq!(candidate.time_step_s(), load.time_step_s());
    assert_eq!(candidate.reports().0.class, load.reports().0.class);
    assert_eq!(candidate.reports().1.class, load.reports().1.class);
    let mut actual_next = [0.0; 2];
    let mut expected_next = [0.0; 2];
    candidate
        .step_mean_into(&[0.3, 0.4], &mut actual_next)
        .unwrap();
    evolved
        .step_mean_into(&[0.3, 0.4], &mut expected_next)
        .unwrap();
    assert_eq!(actual_next, expected_next);
    assert_eq!(candidate.pole_states(), evolved.pole_states());
    assert_eq!(load.pole_states(), accepted);
    assert_eq!(load.source_model(), &model);
    assert_eq!(load.reports().1.class, CertificateClass::GridOnly);

    // A zero-state, zero-D terminal is a valid original harmonic constraint;
    // an ODE consumer may separately require an invertible D.
    let static_model = MimoModel {
        ports: 2,
        poles: vec![],
        residues: vec![],
        d: vec![0.0; 4],
        e: vec![0.0; 4],
        reciprocal: true,
    };
    let static_load = prepare(&static_model);
    let static_view = static_load.continuous_operator();
    assert_eq!(static_view.dimension(), 0);
    let static_space = static_view.state_space(64, &mut || false).unwrap();
    assert!(static_space.a.is_empty() && static_space.b.is_empty() && static_space.c.is_empty());
    assert_eq!(static_space.d, static_model.d);
    assert!(static_view.state_space(63, &mut || false).is_err());
}

/// G0/G4: caller byte limits, validity band and cancellation are checked
/// without changing the admitted source, accepted history or report strength.
#[test]
fn continuous_budget_band_and_cancellation_refusals_preserve_the_admitted_owner() {
    let model = source();
    let mut load = prepare(&model);
    load.step_mean_into(&[-0.2, 0.6], &mut [0.0; 2]).unwrap();
    let before = load.pole_states().to_vec();
    let limits = *load.limits();
    let view = load.continuous_operator();
    let n = view.dimension();
    let p = view.ports();
    let exact_bytes = 8 * (n * n + 2 * n * p + 2 * p * p);
    assert!(view.state_space(exact_bytes - 1, &mut || false).is_err());
    assert!(view.state_space(exact_bytes, &mut || false).is_ok());
    for omega in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        limits.band_rad_s[0] * 0.5,
        limits.band_rad_s[1] * 1.01,
    ] {
        assert!(view.admit_frequency(omega).is_err());
    }
    for omega in limits.band_rad_s {
        view.admit_frequency(omega).unwrap();
    }
    assert!(view.state_space(exact_bytes, &mut || true).is_err());
    let mut polls = 0;
    assert!(
        view.state_space(exact_bytes, &mut || {
            polls += 1;
            polls == 2
        })
        .is_err()
    );
    assert_eq!(polls, 2);
    let coordinates = view.real_coordinates();
    assert!(
        view.candidate_from_coordinates(&coordinates[..n - 1], &[0.0; 2])
            .is_err()
    );
    assert!(
        view.candidate_from_coordinates(&coordinates, &[0.0])
            .is_err()
    );
    assert!(
        view.candidate_from_coordinates(&coordinates, &[f64::INFINITY, 0.0])
            .is_err()
    );
    let mut invalid_coordinates = coordinates;
    invalid_coordinates[3] = f64::NAN;
    assert!(
        view.candidate_from_coordinates(&invalid_coordinates, &[0.0; 2])
            .is_err()
    );
    assert_eq!(load.pole_states(), before);
    assert_eq!(load.source_model(), &model);
    assert_eq!(load.limits(), &limits);
    assert_eq!(load.time_step_s(), DT);
    assert_eq!(load.reports().1.class, CertificateClass::GridOnly);
}
