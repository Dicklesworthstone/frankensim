//! Integration test for deterministic simulation report generation
//! (bead `frankensim-euler-disc-emergent-flagship-t6314.8.7`).

#![allow(missing_docs)]

use fs_euler_disc_e2e::report::{EulerSimulationReport, PhysicalComparisonSection};

#[test]
fn test_simulation_report_generation_with_no_data_placeholder() {
    let report = EulerSimulationReport::new_simulation_only(
        "campaign_test_01",
        "specimen_45mm_chrome",
        12.45,
        15.0,
        0.5,
        1.250,
        0.005,
    );

    assert_eq!(report.campaign_id, "campaign_test_01");
    assert_eq!(report.specimen_id, "specimen_45mm_chrome");
    assert_eq!(report.duration_s, 12.45);
    assert!(report.energy_defect_j > 0.0);
    assert!(matches!(
        report.physical_section,
        PhysicalComparisonSection::NoData { .. }
    ));

    let md = report.render_markdown();
    assert!(md.contains("# Euler Disc Simulation Report"));
    assert!(md.contains("campaign_test_01"));
    assert!(md.contains("NO DATA"));
    assert!(md.contains("No-Claim Disclosure"));
}

#[test]
fn test_simulation_report_with_admitted_physical_scorecard() {
    let mut report = EulerSimulationReport::new_simulation_only(
        "campaign_test_02",
        "specimen_45mm_chrome",
        15.0,
        15.0,
        0.5,
        1.250,
        0.005,
    );

    report.physical_section = PhysicalComparisonSection::AdmittedPhysicalScorecard {
        scorecard_digest: "1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef".to_string(),
        disposition_authority: "independent-disposition-auth-v1".to_string(),
        spin_time_relative_error: 0.035,
    };

    let md = report.render_markdown();
    assert!(md.contains("Scorecard Digest"));
    assert!(md.contains("3.50%"));
}

#[test]
fn test_display_trait_formatting() {
    let report = EulerSimulationReport::new_simulation_only(
        "campaign_test_display",
        "specimen_45mm_chrome",
        10.0,
        15.0,
        0.5,
        1.250,
        0.005,
    );

    let display_str = format!("{}", report);
    assert_eq!(display_str, report.render_markdown());

    let no_data = PhysicalComparisonSection::no_data();
    assert!(format!("{}", no_data).starts_with("NO DATA"));
}
