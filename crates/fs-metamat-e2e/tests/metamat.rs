//! End-to-end battery: a stiffness-density frontier whose every point is
//! certified stable (PSD) and admissible (≤ Voigt), proving the solid optimal.

use fs_evidence::Color;
use fs_metamat_e2e::{default_radii, run_campaign};

#[test]
fn the_frontier_is_stable_admissible_and_solid_optimal() {
    let report = run_campaign(10, &default_radii());
    // a real frontier with positive solid stiffness.
    assert!(report.frontier.len() >= 4);
    assert!(report.c_solid > 0.0, "c_solid {}", report.c_solid);
    // STABILITY, PROVEN: every effective tensor is positive-definite.
    assert!(report.all_stable, "a cell was not PSD-stable");
    // ADMISSIBILITY, PROVEN: every C11 respects the Voigt upper bound.
    assert!(report.all_admissible, "a cell violated the Voigt bound");
    // stiffness falls as porosity grows, and density falls too.
    assert!(
        report.stiffness_monotone,
        "stiffness not monotone in porosity"
    );
    let dens: Vec<f64> = report.frontier.iter().map(|p| p.density).collect();
    assert!(
        dens.windows(2).all(|w| w[1] <= w[0] + 1e-12),
        "density not monotone"
    );
    // the Voigt bound PROVES no porous cell beats solid on specific stiffness.
    assert!(report.solid_is_specific_optimal);
    // the frontier is Verified.
    assert!(matches!(report.stability_color, Color::Verified { .. }));
    println!(
        "{{\"campaign\":\"metamatcert\",\"points\":{},\"c_solid\":{:.4},\"all_stable\":{},\
         \"all_admissible\":{},\"monotone\":{},\"solid_optimal\":{},\"frontier\":{:?}}}",
        report.frontier.len(),
        report.c_solid,
        report.all_stable,
        report.all_admissible,
        report.stiffness_monotone,
        report.solid_is_specific_optimal,
        report
            .frontier
            .iter()
            .map(|p| (p.density, p.c11))
            .collect::<Vec<_>>(),
    );
}

#[test]
fn porosity_actually_reduces_stiffness() {
    let report = run_campaign(10, &default_radii());
    let solid = report.frontier.first().unwrap();
    let porous = report.frontier.last().unwrap();
    assert!(
        porous.c11 < solid.c11,
        "porous {} !< solid {}",
        porous.c11,
        solid.c11
    );
    assert!(porous.density < solid.density);
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(10, &default_radii());
    let b = run_campaign(10, &default_radii());
    assert_eq!(a.c_solid.to_bits(), b.c_solid.to_bits());
    assert_eq!(
        a.frontier.last().unwrap().c11.to_bits(),
        b.frontier.last().unwrap().c11.to_bits()
    );
}
