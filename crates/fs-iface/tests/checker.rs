//! Battery for the coupling-graph static checker (addendum Proposal 13).
//! Covers legal/illegal continuity, certified/unstable/unknown saddle
//! pairings (the conservative-by-default soundness rule), malformed graphs
//! (missing field, self-loop), the empty graph, multi-fault localization,
//! determinism, and the FEEC exact-sequence type lattice.

use fs_iface::{
    CheckReport, CouplingGraph, CouplingRole, PairingRegistry, PairingVerdict, Severity, SpaceType,
    check,
};

fn rejects_on(report: &CheckReport, check_name: &str) -> usize {
    report
        .findings
        .iter()
        .filter(|f| f.severity == Severity::Reject && f.check == check_name)
        .count()
}

#[test]
fn feec_periodic_table_is_the_exact_sequence() {
    // H(grad) --grad--> H(curl) --curl--> H(div) --div--> L²  (form degrees 0..3)
    assert_eq!(SpaceType::HGrad.form_degree(), 0);
    assert_eq!(SpaceType::HCurl.form_degree(), 1);
    assert_eq!(SpaceType::HDiv.form_degree(), 2);
    assert_eq!(SpaceType::L2.form_degree(), 3);
    assert_eq!(SpaceType::HDiv.name(), "H(div)");
}

#[test]
fn empty_graph_is_vacuously_legal() {
    let report = check(&CouplingGraph::new(), &PairingRegistry::standard());
    assert!(report.admitted, "{}", report.diagnosis());
    assert!(report.findings.is_empty());
}

#[test]
fn same_space_continuity_is_legal() {
    // H(div) normal-flux continuity across an interface: both sides H(div).
    let g = CouplingGraph::new()
        .field("flux_left", SpaceType::HDiv)
        .field("flux_right", SpaceType::HDiv)
        .couple(
            "interface",
            "flux_left",
            "flux_right",
            CouplingRole::Continuity,
        );
    assert!(check(&g, &PairingRegistry::standard()).admitted);
}

#[test]
fn cross_space_continuity_is_illegal() {
    // The addendum's example: coupling a pressure (L²) to a displacement
    // trace (H(grad)) as a continuity coupling is nonsense — incompatible
    // trace spaces.
    let g = CouplingGraph::new()
        .field("pressure", SpaceType::L2)
        .field("displacement", SpaceType::HGrad)
        .couple(
            "bad-iface",
            "pressure",
            "displacement",
            CouplingRole::Continuity,
        );
    let report = check(&g, &PairingRegistry::standard());
    assert!(!report.admitted);
    assert_eq!(rejects_on(&report, "coupling.continuity"), 1);
    // localized to the offending coupling.
    assert_eq!(report.findings[0].coupling, "bad-iface");
}

#[test]
fn certified_saddle_pairing_is_legal() {
    // Mixed Poisson / Darcy: H(div) flux with L² pressure — LBB-stable.
    let g = CouplingGraph::new()
        .field("flux", SpaceType::HDiv)
        .field("pressure", SpaceType::L2)
        .couple("darcy", "flux", "pressure", CouplingRole::Saddle);
    assert!(check(&g, &PairingRegistry::standard()).admitted);
    // Stokes Taylor-Hood: H¹ velocity with L² pressure — LBB-stable.
    let g2 = CouplingGraph::new()
        .field("velocity", SpaceType::HGrad)
        .field("pressure", SpaceType::L2)
        .couple("stokes", "velocity", "pressure", CouplingRole::Saddle);
    assert!(check(&g2, &PairingRegistry::standard()).admitted);
}

#[test]
fn known_unstable_saddle_pairing_is_rejected() {
    // Equal-order H¹ velocity–pressure (P1–P1) violates LBB.
    let g = CouplingGraph::new()
        .field("velocity", SpaceType::HGrad)
        .field("pressure", SpaceType::HGrad)
        .couple("p1p1", "velocity", "pressure", CouplingRole::Saddle);
    let report = check(&g, &PairingRegistry::standard());
    assert!(!report.admitted);
    assert_eq!(rejects_on(&report, "coupling.infsup"), 1);
    // the teaching fix names a certified partner (L2).
    assert!(
        report.findings[0].fix.contains("L2"),
        "fix: {}",
        report.findings[0].fix
    );
}

#[test]
fn unknown_saddle_pairing_is_rejected_conservatively() {
    // The load-bearing soundness rule: a pairing NOT in the registry is
    // illegal-until-certified, never silently admitted.
    let g = CouplingGraph::new()
        .field("a", SpaceType::HCurl)
        .field("b", SpaceType::HDiv)
        .couple("exotic", "a", "b", CouplingRole::Saddle);
    let report = check(&g, &PairingRegistry::standard());
    assert!(
        !report.admitted,
        "an unrecognized pairing must NOT be admitted by default"
    );
    assert_eq!(rejects_on(&report, "coupling.infsup"), 1);
    assert!(report.findings[0].what.contains("illegal until certified"));
}

#[test]
fn missing_field_reference_is_a_localized_error() {
    let g = CouplingGraph::new().field("flux", SpaceType::HDiv).couple(
        "dangling",
        "flux",
        "ghost",
        CouplingRole::Saddle,
    );
    let report = check(&g, &PairingRegistry::standard());
    assert!(!report.admitted);
    assert_eq!(rejects_on(&report, "graph.field"), 1);
    assert!(
        report.findings[0].what.contains("ghost"),
        "names the missing field"
    );
}

#[test]
fn self_loop_coupling_is_rejected() {
    let g = CouplingGraph::new().field("u", SpaceType::HGrad).couple(
        "selfie",
        "u",
        "u",
        CouplingRole::Continuity,
    );
    let report = check(&g, &PairingRegistry::standard());
    assert!(!report.admitted);
    assert_eq!(rejects_on(&report, "graph.self-loop"), 1);
}

#[test]
fn multiple_faults_are_each_localized() {
    let g = CouplingGraph::new()
        .field("v", SpaceType::HGrad)
        .field("p", SpaceType::HGrad)
        .field("flux", SpaceType::HDiv)
        .field("mult", SpaceType::L2)
        .couple("good", "flux", "mult", CouplingRole::Saddle) // certified → no finding
        .couple("p1p1", "v", "p", CouplingRole::Saddle) // unstable
        .couple("dangling", "v", "nope", CouplingRole::Saddle); // missing field
    let report = check(&g, &PairingRegistry::standard());
    assert!(!report.admitted);
    // two rejects, each naming its own coupling; the certified one is silent.
    assert_eq!(report.findings.len(), 2);
    let ids: Vec<&str> = report
        .findings
        .iter()
        .map(|f| f.coupling.as_str())
        .collect();
    assert!(ids.contains(&"p1p1") && ids.contains(&"dangling") && !ids.contains(&"good"));
}

#[test]
fn check_is_deterministic() {
    let g = CouplingGraph::new()
        .field("v", SpaceType::HGrad)
        .field("p", SpaceType::HGrad)
        .couple("p1p1", "v", "p", CouplingRole::Saddle);
    let r = PairingRegistry::standard();
    assert_eq!(check(&g, &r), check(&g, &r));
}

#[test]
fn empty_registry_rejects_everything_conservatively() {
    // With no certified pairings, even a normally-legal saddle is unknown.
    let g = CouplingGraph::new()
        .field("flux", SpaceType::HDiv)
        .field("p", SpaceType::L2)
        .couple("darcy", "flux", "p", CouplingRole::Saddle);
    let report = check(&g, &PairingRegistry::empty());
    assert!(!report.admitted, "an empty registry certifies nothing");
    assert!(matches!(
        PairingRegistry::empty().classify_saddle(SpaceType::HDiv, SpaceType::L2),
        PairingVerdict::Unknown
    ));
}

#[test]
fn registry_certifies_the_standard_pairings() {
    let r = PairingRegistry::standard();
    assert!(matches!(
        r.classify_saddle(SpaceType::HDiv, SpaceType::L2),
        PairingVerdict::Certified { .. }
    ));
    assert!(matches!(
        r.classify_saddle(SpaceType::L2, SpaceType::L2),
        PairingVerdict::Unstable { .. }
    ));
    assert_eq!(r.certified_partners(SpaceType::HDiv), vec![SpaceType::L2]);
}
