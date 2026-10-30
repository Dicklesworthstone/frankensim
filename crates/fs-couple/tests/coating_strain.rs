//! Complete source-film coverage and coherent harmonic strain observation.
#![deny(unsafe_code)]

#[path = "coating_shell.rs"]
mod fixture;

use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::shear_shell::{
    CoatingSide, MaterialShellImage, ShellProjectionBudget, ShellReductionOptions,
    coating::strain::ShellCoatingStrainBudget,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_plate::shell::mitc::ShellMass;

fn image(gate: &CancelGate) -> MaterialShellImage {
    let (shell, directors, constraints) = fixture::coupon(gate);
    let shell = shell
        .with_coatings(
            vec![
                fixture::layer(CoatingSide::Positive, 20e-6),
                fixture::layer(CoatingSide::Positive, 15e-6),
                fixture::layer(CoatingSide::Negative, 25e-6),
            ],
            gate,
        )
        .unwrap();
    let mut port = vec![0.0; shell.mesh().nodes.len() * 6];
    port[6 * 5 + 2] = 1.0;
    shell
        .reduce(
            &directors,
            &ShellReductionOptions {
                constraints,
                mass: ShellMass::Consistent,
                drilling_alpha: 1e-3,
                eigenvalue_window: (0.0, 1e9),
                n_modes: 3,
                assembly_budget: fixture::budget(&shell),
                projection_budget: ShellProjectionBudget {
                    max_modes: 3,
                    max_projection_terms: 100_000_000,
                    reference_resistance_n_s_m: 1.0,
                },
            },
            &[port],
            gate,
        )
        .unwrap()
}

fn budget(evaluations: usize) -> ShellCoatingStrainBudget {
    ShellCoatingStrainBudget {
        maximum_bytes: 8_000_000,
        maximum_work: 100_000_000,
        evaluations,
    }
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 5e-11 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn all_original_films_faces_and_elements_share_one_coherent_body_motion() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let saved_modes = image.physical_modes().to_vec();
    let mut observer = image.prepare_coating_strain(budget(3), &gate).unwrap();
    assert_eq!(observer.film_count(), 3);
    assert_eq!(observer.resources().coated_elements, 18);
    assert_eq!(observer.resources().film_element_pairs, 54);
    assert_eq!(observer.film_sources().count(), 3);
    let q = [
        C64 {
            re: 2e-6,
            im: -3e-7,
        },
        C64 {
            re: -0.6e-6,
            im: 1e-6,
        },
        C64 {
            re: 0.4e-6,
            im: -0.2e-6,
        },
    ];
    let original = observer.evaluate(&q, &gate).unwrap();
    let rotated = observer
        .evaluate(
            &q.map(|z| C64 {
                re: -2.0 * z.im,
                im: 2.0 * z.re,
            }),
            &gate,
        )
        .unwrap();
    for (i, (a, b)) in original.films.iter().zip(&rotated.films).enumerate() {
        assert_eq!(a.film_index, i);
        assert_eq!(a.component_index, None);
        assert!(std::ptr::eq(a.source, &image.source().coatings()[i]));
        assert_eq!(a.covered_elements, 18);
        assert!(a.ruling_element < 18 && a.ruling_control_point < 6 && a.ruling_z_face < 2);
        assert!(a.frobenius_upper_bound > 0.0);
        near(b.frobenius_upper_bound, 2.0 * a.frobenius_upper_bound);
        near(b.limit_ratio, 2.0 * a.limit_ratio);
        near(
            b.sufficient_amplitude_scale.unwrap(),
            0.5 * a.sufficient_amplitude_scale.unwrap(),
        );
    }
    let faces: Vec<_> = original.films.iter().map(|f| f.source.faces_m()).collect();
    near(faces[0].0, 0.0015);
    near(faces[0].1, 0.00152);
    near(faces[1].0, faces[0].1);
    near(faces[1].1, 0.001535);
    near(faces[2].0, -0.001525);
    near(faces[2].1, -0.0015);
    near(
        rotated.sufficient_common_amplitude_scale.unwrap(),
        0.5 * original.sufficient_common_amplitude_scale.unwrap(),
    );
    let zero = observer.evaluate(&[C64::ZERO; 3], &gate).unwrap();
    assert!(zero.within_authored_envelope);
    assert_eq!(zero.sufficient_common_amplitude_scale, None);
    assert!(
        zero.films
            .iter()
            .all(|f| f.frobenius_upper_bound == 0.0 && f.limit_ratio == 0.0)
    );
    assert_eq!(image.physical_modes(), saved_modes);
}

#[test]
fn complete_setup_and_query_allowances_refuse_before_partial_reports() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let resources = *image
        .prepare_coating_strain(budget(1), &gate)
        .unwrap()
        .resources();
    let exact = ShellCoatingStrainBudget {
        maximum_bytes: resources.retained_bytes + resources.evaluation_bytes,
        maximum_work: resources.total_work,
        evaluations: 1,
    };
    assert!(
        image
            .prepare_coating_strain(
                ShellCoatingStrainBudget {
                    maximum_bytes: exact.maximum_bytes - 1,
                    ..exact
                },
                &gate
            )
            .is_err()
    );
    assert!(
        image
            .prepare_coating_strain(
                ShellCoatingStrainBudget {
                    maximum_work: exact.maximum_work - 1,
                    ..exact
                },
                &gate
            )
            .is_err()
    );
    assert!(
        image
            .prepare_coating_strain(
                ShellCoatingStrainBudget {
                    maximum_work: usize::MAX,
                    evaluations: usize::MAX,
                    ..exact
                },
                &gate
            )
            .is_err()
    );
    let mut observer = image.prepare_coating_strain(exact, &gate).unwrap();
    assert!(observer.evaluate(&[C64::ZERO; 2], &gate).is_err());
    assert!(
        observer
            .evaluate(
                &[C64 {
                    re: f64::NAN,
                    im: 0.0
                }; 3],
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        observer.evaluate(&[C64::ZERO; 3], &cancelled),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(observer.completed_evaluations(), 0);
    observer.evaluate(&[C64::ZERO; 3], &gate).unwrap();
    assert_eq!(observer.completed_evaluations(), 1);
    assert!(observer.evaluate(&[C64::ZERO; 3], &gate).is_err());
}
