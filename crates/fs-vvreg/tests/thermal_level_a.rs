//! G0/G1/G3 checks for the reference-only Level-A thermal corpus.

use std::collections::BTreeSet;

use fs_evidence::{ColorRank, NumericalKind};
use fs_qty::QtyAny;
use fs_vvreg::corpus::{ContextValue, CorpusEnvelope, EvidenceLevel, PayloadRetention, corpus};
use fs_vvreg::partition::{DatasetPurpose, PartitionLedger};
use fs_vvreg::thermal_level_a::{
    ThermalLevelAAcceptance, ThermalLevelAFamily, ThermalLevelAKind, thermal_level_a_cases,
};

const MANIFEST: &[u8] =
    include_bytes!("../../../data/vv-corpus/thermal-level-a/thermal-level-a-v1.tsv");

fn case(id: &str) -> &'static fs_vvreg::thermal_level_a::ThermalLevelACase {
    thermal_level_a_cases()
        .iter()
        .find(|case| case.id == id)
        .unwrap_or_else(|| panic!("missing Level-A case {id}"))
}

fn assert_reference(id: &str, expected: f64) {
    let observed = case(id).reference_value_si;
    let scale = expected.abs().max(1.0);
    assert!(
        (observed - expected).abs() <= 2.0e-14 * scale,
        "{id}: retained {observed:.17e}, independent formula {expected:.17e}"
    );
}

#[test]
fn catalog_is_unique_and_covers_every_requested_family() {
    let cases = thermal_level_a_cases();
    assert_eq!(cases.len(), 19);

    let ids = cases.iter().map(|case| case.id).collect::<BTreeSet<_>>();
    assert_eq!(ids.len(), cases.len());
    let families = cases
        .iter()
        .map(|case| case.family)
        .collect::<BTreeSet<_>>();
    assert_eq!(
        families,
        BTreeSet::from([
            ThermalLevelAFamily::SteadyConduction,
            ThermalLevelAFamily::Fin,
            ThermalLevelAFamily::LumpedTransient,
            ThermalLevelAFamily::ConvectionLimit,
            ThermalLevelAFamily::Radiation,
            ThermalLevelAFamily::Contact,
            ThermalLevelAFamily::ManufacturedPrimal,
            ThermalLevelAFamily::ManufacturedAdjoint,
        ])
    );
    assert_eq!(
        cases
            .iter()
            .filter(|case| case.kind == ThermalLevelAKind::AnalyticReference)
            .count(),
        12
    );
    assert_eq!(
        cases
            .iter()
            .filter(|case| case.kind == ThermalLevelAKind::ManufacturedTarget)
            .count(),
        7
    );
}

#[test]
fn retained_closed_form_values_reproduce_from_independent_formulas() {
    assert_reference("thermal-a-slab-dirichlet", 20.0 * 40.0 / 0.2);
    assert_reference("thermal-a-slab-robin", 100.0 / (0.1 / 10.0 + 1.0 / 100.0));
    assert_reference(
        "thermal-a-slab-uniform-source",
        100_000.0 * 0.1_f64.powi(2) / (8.0 * 10.0),
    );
    assert_reference(
        "thermal-a-rectangle-linear",
        300.0 + 20.0 * 0.5 + 40.0 * 0.25,
    );
    assert_reference(
        "thermal-a-cylinder-shell",
        2.0 * std::f64::consts::PI * 15.0 / 2.0_f64.ln(),
    );
    assert_reference(
        "thermal-a-sphere-shell",
        4.0 * std::f64::consts::PI * 15.0 / (1.0 / 0.05 - 1.0 / 0.1),
    );
    assert_reference("thermal-a-fin-efficiency", 1.0_f64.tanh());
    assert_reference("thermal-a-lumped-transient", (-1.0_f64).exp());
    assert_reference("thermal-a-duct-nu-cwt", 3.66);
    assert_reference("thermal-a-duct-nu-chf", 4.36);
    assert_reference("thermal-a-parallel-plate-view-factor", 1.0);
    assert_reference(
        "thermal-a-contact-series",
        0.01 / (10.0 * 0.01) + 0.1 + 0.02 / (20.0 * 0.01),
    );
}

#[test]
fn retained_manifest_rows_match_the_typed_catalog_exactly() {
    let text = std::str::from_utf8(MANIFEST).expect("tracked manifest is UTF-8");
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some(
            "schema_version\tcase_id\tfamily\tkind\tmetric\treference_value_si\tformula\tacceptance\tcontext\tstatus"
        )
    );
    let rows = lines.collect::<Vec<_>>();
    assert_eq!(rows.len(), thermal_level_a_cases().len());
    for (row, case) in rows.iter().zip(thermal_level_a_cases()) {
        let fields = row.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), 10, "malformed manifest row: {row}");
        assert_eq!(fields[0], "1");
        assert_eq!(fields[1], case.id);
        assert_eq!(fields[2], case.family.name());
        assert_eq!(fields[3], case.kind.name());
        assert_eq!(fields[4], case.metric);
        assert_eq!(
            fields[5]
                .parse::<f64>()
                .expect("reference parses")
                .to_bits(),
            case.reference_value_si.to_bits(),
            "{} reference drifted",
            case.id
        );
        assert_eq!(fields[6], case.formula);
        assert!(!fields[7].is_empty());
        assert!(!fields[8].is_empty());
        assert_eq!(
            fields[9],
            if case.kind == ThermalLevelAKind::AnalyticReference {
                "reference-only"
            } else {
                "target-only"
            }
        );
    }
}

#[test]
fn seeded_registry_binds_every_reference_but_returns_no_solver_claim() {
    let manifest_hash = fs_blake3::hash_bytes(MANIFEST);
    let partitions = PartitionLedger::capture(corpus());
    for case in thermal_level_a_cases() {
        let dataset = corpus().dataset(case.id).expect("seeded Level-A row");
        assert_eq!(dataset.evidence_level(), EvidenceLevel::Analytic);
        assert_eq!(dataset.physical_claim_cap(), ColorRank::Estimated);
        assert!(matches!(
            dataset.raw_payload(),
            PayloadRetention::DerivedOnly { .. }
        ));
        assert_eq!(dataset.raw_payload().artifact().digest, manifest_hash);
        assert_eq!(
            dataset.raw_payload().artifact().byte_len,
            MANIFEST.len() as u64
        );
        assert_eq!(dataset.acceptance_envelopes().len(), 1);

        let context = case
            .context
            .iter()
            .map(|axis| ContextValue {
                name: axis.name.to_string(),
                value: QtyAny::new(axis.lo + (axis.hi - axis.lo) * 0.5, axis.dims),
            })
            .collect::<Vec<_>>();
        let evidence = corpus()
            .query(&partitions, case.id, DatasetPurpose::Validation, &context)
            .expect("exact Level-A context admits");
        assert_eq!(evidence.numerical.kind, NumericalKind::NoClaim);
        assert!(evidence.statistical.rel_width(1.0).is_infinite());
        assert!(evidence.model.discrepancy_rel.is_infinite());
        for axis in case.context {
            assert_eq!(
                evidence.model.validity.bound(axis.name),
                Some((axis.lo, axis.hi))
            );
        }
    }
}

#[test]
fn manufactured_targets_cover_two_orders_boundary_types_and_adjoint() {
    let targets = thermal_level_a_cases()
        .iter()
        .filter(|case| case.kind == ThermalLevelAKind::ManufacturedTarget)
        .collect::<Vec<_>>();
    let degrees = targets
        .iter()
        .map(|case| {
            case.context
                .iter()
                .find(|axis| axis.name == "element-degree")
                .expect("MMS target declares element degree")
                .lo
                .to_bits()
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(
        degrees,
        BTreeSet::from([1.0_f64.to_bits(), 2.0_f64.to_bits()])
    );
    assert!(targets.iter().any(|case| case.id.contains("neumann")));
    assert!(targets.iter().any(|case| case.id.contains("robin")));
    assert!(targets.iter().any(|case| {
        case.family == ThermalLevelAFamily::ManufacturedAdjoint && case.id.contains("p1-adjoint")
    }));
    assert!(targets.iter().any(|case| {
        case.family == ThermalLevelAFamily::ManufacturedAdjoint && case.id.contains("p2-adjoint")
    }));

    for case in targets {
        let ThermalLevelAAcceptance::OrderGate {
            theoretical,
            tolerance,
        } = case.acceptance
        else {
            panic!("{} must use an order gate", case.id);
        };
        assert_eq!(tolerance, 0.2);
        let dataset = corpus().dataset(case.id).expect("seeded MMS target");
        assert_eq!(
            dataset.acceptance_envelopes()[0].envelope,
            CorpusEnvelope::Interval {
                lo: theoretical - 0.2,
                hi: theoretical + 0.2,
            }
        );
    }
}

#[test]
fn audit_surfaces_reference_only_status_as_warn_not_green_solver_evidence() {
    let report = corpus().audit();
    assert!(report.is_clean());
    for case in thermal_level_a_cases() {
        let row = report
            .rows()
            .iter()
            .find(|row| row.dataset_id() == case.id)
            .expect("Level-A audit row");
        assert_eq!(row.status(), "WARN");
        assert!(report.warnings().iter().any(|warning| {
            warning.contains(&format!(
                "dataset={} claim_gap=raw_payload.original",
                case.id
            ))
        }));
    }
}
