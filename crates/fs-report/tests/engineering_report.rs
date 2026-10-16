//! Tests for deterministic HTML engineering report and JSON twin.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.9`

use fs_evidence::Color;
use fs_ladder::{ConvergenceEvaluator, ConvergencePlan, MeshRung};
use fs_report::{EngineeringReport, MaterialReportItem, NoClaimItem, QoiReportItem};
use fs_uq::{ParameterUncertainty, PropagationMethod, UqPlan, UqPropagator};

#[test]
fn test_engineering_report_html_rendering() {
    let mut report = EngineeringReport::new("run_cooling_001", "Dual-Sided Cold Plate Simulation");

    // Add QoI
    report = report.with_qoi(QoiReportItem {
        name: "junction_maximum".to_string(),
        description: "Maximum junction temperature across silicon dies".to_string(),
        nominal_value: 342.15,
        unit: "K".to_string(),
        color: Color::Verified {
            lo: 341.0,
            hi: 343.3,
        },
        discretization_error: 0.25,
        parameter_uncertainty: 0.85,
        surrogate_error: 0.0,
        total_uncertainty_budget: 1.10,
        source_root: "qoi://cooling/thermal_qoi/junction_maximum".to_string(),
    });

    // Convergence
    let conv_plan = ConvergencePlan::new("junction_maximum", 2.0)
        .with_rung(MeshRung::new(
            0,
            "mesh_coarse",
            0.04,
            "m",
            1000,
            "junction_maximum",
            350.0,
            "K",
        ))
        .with_rung(MeshRung::new(
            1,
            "mesh_medium",
            0.02,
            "m",
            8000,
            "junction_maximum",
            342.5,
            "K",
        ))
        .with_rung(MeshRung::new(
            2,
            "mesh_fine",
            0.01,
            "m",
            64000,
            "junction_maximum",
            340.625,
            "K",
        ));
    let conv_res = ConvergenceEvaluator::evaluate(&conv_plan);
    report = report.with_convergence(conv_res);

    // Uncertainty
    let uq_plan = UqPlan::new("junction_maximum", PropagationMethod::MonteCarlo, 200)
        .with_parameter(ParameterUncertainty::gaussian(
            "ambient_temp",
            300.0,
            5.0,
            "K",
        ))
        .with_compliance_threshold(355.0);
    let uq_res = UqPropagator::run(&uq_plan, |p| p[0] + 42.15);
    report = report.with_uncertainty(uq_res);

    // Material
    report.materials.push(MaterialReportItem {
        region_name: "cold_plate_base".to_string(),
        material_card_id: "mat_al6061_t6".to_string(),
        thermal_conductivity: "167.0 W/(m·K)".to_string(),
        specific_heat: "896.0 J/(kg·K)".to_string(),
        density: "2700.0 kg/m^3".to_string(),
        source_pack: "frankensim-materials-2026.1".to_string(),
    });

    // No-claims
    report.no_claims.push(NoClaimItem {
        component: "radiation_coupling".to_string(),
        status: "Unmodeled".to_string(),
        statement:
            "Surface-to-surface radiative heat transfer is neglected (<1% contribution at <400K)."
                .to_string(),
    });

    report.replay_command = "frankensim solve --project cooling --seed 0x0517".to_string();

    let html = report.render_html();
    assert!(html.contains("<!DOCTYPE html>"));
    assert!(html.contains("Dual-Sided Cold Plate Simulation"));
    assert!(html.contains("junction_maximum"));
    assert!(html.contains("342.15"));
    assert!(html.contains("Verified"));
    assert!(html.contains("mat_al6061_t6"));
    assert!(html.contains("radiation_coupling"));
    assert!(html.contains("frankensim solve"));

    let json = report.render_json();
    assert!(json.contains("\"schema\": \"frankensim.report.engineering.v1\""));
    assert!(json.contains("\"name\": \"junction_maximum\""));
    assert!(json.contains("\"value\": 342.150000"));
}

#[test]
fn test_engineering_report_determinism() {
    let mut r1 = EngineeringReport::new("run_001", "Deterministic Test");
    r1 = r1.with_qoi(QoiReportItem {
        name: "qoi_test".to_string(),
        description: "Test QoI".to_string(),
        nominal_value: 100.0,
        unit: "W".to_string(),
        color: Color::Verified {
            lo: 99.0,
            hi: 101.0,
        },
        discretization_error: 0.1,
        parameter_uncertainty: 0.2,
        surrogate_error: 0.0,
        total_uncertainty_budget: 0.3,
        source_root: "qoi://test/1".to_string(),
    });

    let r2 = r1.clone();

    assert_eq!(
        r1.render_html(),
        r2.render_html(),
        "HTML renders bit-identically"
    );
    assert_eq!(
        r1.render_json(),
        r2.render_json(),
        "JSON renders bit-identically"
    );
    assert_eq!(r1.content_hash(), r2.content_hash(), "Content hashes match");
}

#[test]
fn test_engineering_report_escaping() {
    let mut r = EngineeringReport::new("run_injection", "<script>alert('xss')</script>");
    r.no_claims.push(NoClaimItem {
        component: "<b>bold</b>".to_string(),
        status: "Unmodeled & Risky".to_string(),
        statement: "\"Quotes\" & <tags>".to_string(),
    });

    let html = r.render_html();
    assert!(
        !html.contains("<script>alert('xss')</script>"),
        "script tag must be escaped"
    );
    assert!(html.contains("&lt;script&gt;alert(&#39;xss&#39;)&lt;/script&gt;"));
    assert!(html.contains("&lt;b&gt;bold&lt;/b&gt;"));
    assert!(html.contains("&amp;"));
}

#[test]
fn g0_uq_report_preserves_refusals_and_empirical_authority() {
    let plan = UqPlan::new("constant", PropagationMethod::MonteCarlo, 2)
        .with_parameter(ParameterUncertainty::gaussian("a", 0.0, 1.0, "1"))
        .with_compliance_threshold(3.0);
    let result = UqPropagator::run(&plan, |_| 2.0);
    let report = EngineeringReport::new("uq-complete", "UQ").with_uncertainty(result);
    let html = report.render_html();
    assert!(html.contains("Empirical compliance frequency"));
    assert!(html.contains("Observed range") && html.contains("Estimated"));
    assert!(!html.contains("Probability of Compliance"));
    let json = report.render_json();
    assert!(json.contains("\"mean\": 2,") && json.contains("\"p_compliance\": 1,"));
    assert!(!json.contains("Some(") && !json.contains("None"));
    // Retain the actual rendered twin for an independent JSON-reader check.
    println!("uq-report-json-begin\n{json}\nuq-report-json-end");

    let mut unsupported = plan;
    unsupported.method = PropagationMethod::QuasiMonteCarlo;
    let mut result = UqPropagator::run(&unsupported, |_| panic!("unsupported model call"));
    result.rejection_reason = Some("missing <joint> & \"measure\"\nrequest declaration".into());
    let report = EngineeringReport::new("uq-refused", "UQ").with_uncertainty(result);
    let html = report.render_html();
    assert!(html.contains("refused") && html.contains("missing &lt;joint&gt; &amp;"));
    assert!(!html.contains("Empirical compliance frequency"));
    assert!(!html.contains("Sample statistics:"));
    let json = report.render_json();
    assert!(json.contains("\"mean\": null") && json.contains("\"p_compliance\": null"));
    assert!(json.contains("\"observed_range\": null"));
    assert!(!json.contains("Some(") && !json.contains("None"));
    println!("uq-report-json-begin\n{json}\nuq-report-json-end");
}
