//! End-to-end battery: a diverse fabricable program family is illuminated, the
//! best matches the target, and every simplification's certificate is re-verified.

use fs_evidence::Color;
use fs_grammar_e2e::{build_program, run_campaign, target};
use fs_shapeprog::max_sdf_discrepancy;

#[test]
fn a_fabricable_program_family_is_illuminated_and_simplified_soundly() {
    let report = run_campaign(0.2, 0.03);
    // ILLUMINATION: a diverse archive, not one model.
    assert!(
        report.num_elites >= 5,
        "too few niches: {}",
        report.num_elites
    );
    assert!(report.coverage > 0.0 && report.qd_score > 0.0); // fitness = 1/(1+discrepancy)
    // the best program genuinely matches the peanut target.
    assert!(
        report.best_discrepancy < 0.2,
        "best discrepancy {}",
        report.best_discrepancy
    );
    // CERTIFICATE-PRESERVING SIMPLIFICATION: some programs shrank, and EVERY
    // simplification's re-measured error stayed within its certified bound.
    assert!(report.simplified_count > 0, "nothing simplified");
    assert!(report.size_after < report.size_before, "no size reduction");
    assert!(
        report.simplification_sound,
        "a rewrite certificate was violated"
    );
    assert!(
        report.max_certified_error <= 0.03 + 1e-9,
        "certified error {}",
        report.max_certified_error
    );
    // FABRICABILITY: the constraint genuinely discriminates — some elites pass
    // the minimum feature size and some (the thin-sphere ones) do not.
    assert!(report.fab_satisfied > 0);
    assert!(
        report.fab_satisfied < report.num_elites,
        "fab did not discriminate: {}/{}",
        report.fab_satisfied,
        report.num_elites
    );
    // the headline claim is Verified (matches + fabricable + sound).
    assert!(matches!(report.headline_color, Color::Verified { .. }));
    println!(
        "{{\"campaign\":\"grammarforge\",\"niches\":{},\"coverage\":{:.3},\"best_disc\":{:.4},\
         \"best_params\":{:?},\"simplified\":{},\"size\":{}->{},\"max_cert_err\":{:.4},\
         \"sound\":{},\"fab_ok\":{}}}",
        report.num_elites,
        report.coverage,
        report.best_discrepancy,
        report.best_params,
        report.simplified_count,
        report.size_before,
        report.size_after,
        report.max_certified_error,
        report.simplification_sound,
        report.fab_satisfied,
    );
}

#[test]
fn the_certified_simplification_bound_actually_holds() {
    // a program with a droppable tiny offset simplifies with a certified bound
    // that its re-measured SDF change respects.
    let prog = build_program(1.0, 1.0, 0.8, 0.02);
    let simp = fs_shapeprog::simplify(&prog, 0.03);
    let samples: Vec<[f64; 3]> = (0..5)
        .flat_map(|i| (0..5).map(move |j| [-2.0 + f64::from(i), (-1.0 + f64::from(j) * 0.5), 0.0]))
        .collect();
    let actual = max_sdf_discrepancy(&prog, &simp.program, &samples);
    assert!(
        actual <= simp.max_error + 1e-9,
        "actual {} > bound {}",
        actual,
        simp.max_error
    );
}

#[test]
fn the_target_matches_itself_exactly() {
    let t = target();
    let samples: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [1.5, 0.0, 0.0]];
    assert!(max_sdf_discrepancy(&t, &t, &samples).abs() < 1e-12);
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(0.2, 0.03);
    let b = run_campaign(0.2, 0.03);
    assert_eq!(a.num_elites, b.num_elites);
    assert_eq!(a.best_discrepancy.to_bits(), b.best_discrepancy.to_bits());
    assert_eq!(a.size_after, b.size_after);
}
