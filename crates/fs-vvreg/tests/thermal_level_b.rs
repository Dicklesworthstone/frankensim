//! G0 checks for the Level-B thermal cross-code corpus: catalog sanity,
//! fail-closed manifest verification, bit-exact spec-echo/mesh-parity
//! binding, deck hash binding, and corpus registration boundaries.

use fs_blake3::hash_bytes;
use fs_evidence::ColorRank;
use fs_vvreg::corpus::{CorpusEnvelope, EvidenceLevel, PayloadRetention, corpus};
use fs_vvreg::thermal_level_b::{
    ThermalLevelBError, ThermalLevelBMaterial, ThermalLevelBSource, parse_thermal_level_b_manifest,
    thermal_level_b_case, thermal_level_b_cases, thermal_level_b_deck_bytes,
    thermal_level_b_references, verified_thermal_level_b_references, verify_probe_grid,
    verify_spec_echo,
};

const MANIFEST: &str =
    include_str!("../../../data/vv-corpus/thermal-level-b/thermal-level-b-references-v1.tsv");

#[test]
fn catalog_is_unique_and_well_formed() {
    let cases = thermal_level_b_cases();
    assert_eq!(cases.len(), 4);
    for (index, case) in cases.iter().enumerate() {
        assert!(
            cases[index + 1..].iter().all(|other| other.id != case.id),
            "duplicate id {}",
            case.id
        );
        for axis in 0..3 {
            assert!(case.mesh_counts[axis] > 0, "{}: zero cell count", case.id);
            assert!(
                case.mesh_extent[axis] > 0.0,
                "{}: non-positive extent",
                case.id
            );
        }
        assert!(case.acceptance_atol_k > 0.0);
        assert!(!case.probes.is_empty());
        for probe in case.probes {
            for axis in 0..3 {
                assert!(
                    probe[axis] <= case.mesh_counts[axis],
                    "{}: probe {probe:?} outside the vertex grid",
                    case.id
                );
            }
        }
        match case.material {
            ThermalLevelBMaterial::Isotropic { k } => assert!(k > 0.0),
            ThermalLevelBMaterial::Tensor { k } => {
                // Symmetric positive definite via leading principal minors.
                for i in 0..3 {
                    for j in 0..3 {
                        assert_eq!(
                            k[i][j].to_bits(),
                            k[j][i].to_bits(),
                            "{}: asymmetric",
                            case.id
                        );
                    }
                }
                let m1 = k[0][0];
                let m2 = k[0][0] * k[1][1] - k[0][1] * k[1][0];
                let m3 = k[0][0] * (k[1][1] * k[2][2] - k[1][2] * k[2][1])
                    - k[0][1] * (k[1][0] * k[2][2] - k[1][2] * k[2][0])
                    + k[0][2] * (k[1][0] * k[2][1] - k[1][1] * k[2][0]);
                assert!(m1 > 0.0 && m2 > 0.0 && m3 > 0.0, "{}: not PD", case.id);
            }
            ThermalLevelBMaterial::LinearKt { knots } => {
                assert!(knots.len() >= 2);
                for pair in knots.windows(2) {
                    assert!(pair[0].0 < pair[1].0, "{}: knots not increasing", case.id);
                }
                assert!(knots.iter().all(|(_, k)| *k > 0.0));
            }
        }
        let wants_picard = matches!(case.material, ThermalLevelBMaterial::LinearKt { .. });
        assert_eq!(
            case.picard.is_some(),
            wants_picard,
            "{}: picard controls must accompany exactly the nonlinear material",
            case.id
        );
        if let ThermalLevelBSource::PolyXy { q0 } = case.source {
            assert!(q0.is_finite());
        }
    }
}

#[test]
fn committed_manifest_passes_fail_closed_verification() {
    let references = verified_thermal_level_b_references().expect("committed manifest must verify");
    assert_eq!(references.len(), thermal_level_b_cases().len());
    for reference in &references {
        let case = thermal_level_b_case(&reference.case_id).expect("verified case exists");
        verify_spec_echo(case, reference).expect("spec echo bound");
        verify_probe_grid(case, reference).expect("probe grid bound");
        assert!(reference.external_code.starts_with("scikit-fem "));
        assert!(reference.linear_solver.contains("SuperLU"));
        assert!(reference.picard_iterations >= 1);
        assert!(reference.t_min_k <= reference.t_max_k);
        assert_eq!(reference.mesh_blake3.len(), 64);
    }
    // The shared accessor exposes the same verified slice.
    let shared = thermal_level_b_references().expect("shared accessor agrees");
    assert_eq!(shared.len(), references.len());
}

#[test]
fn deck_bytes_are_hash_bound_to_the_manifest() {
    for reference in verified_thermal_level_b_references().expect("verified") {
        let deck = thermal_level_b_deck_bytes(&reference.case_id)
            .expect("every catalog case commits its deck");
        assert_eq!(
            hash_bytes(deck).to_hex(),
            reference.deck_blake3,
            "{}: committed deck bytes are not the deck the external run consumed",
            reference.case_id
        );
    }
}

#[test]
fn nonlinear_case_records_a_real_iteration_count() {
    let references = verified_thermal_level_b_references().expect("verified");
    for reference in &references {
        let case = thermal_level_b_case(&reference.case_id).expect("case");
        if matches!(case.material, ThermalLevelBMaterial::LinearKt { .. }) {
            assert!(
                reference.picard_iterations > 1,
                "{}: a nonlinear solve converging in one Picard step is a red flag",
                reference.case_id
            );
        } else {
            assert_eq!(reference.picard_iterations, 1, "{}", reference.case_id);
        }
    }
}

#[test]
fn parser_refuses_malformed_manifests() {
    // Unknown row kind.
    let err = parse_thermal_level_b_manifest(b"bogus\tcase\tk\tv\n").unwrap_err();
    assert!(matches!(
        err,
        ThermalLevelBError::UnknownKind { line: 1, .. }
    ));

    // Wrong column count.
    let err = parse_thermal_level_b_manifest(b"case_meta\tonly-three\tk\n").unwrap_err();
    assert!(matches!(err, ThermalLevelBError::Columns { line: 1, .. }));

    // Unparsable probe number.
    let err =
        parse_thermal_level_b_manifest(b"probe\tc\t0\t0\t0\t0\t0.0\t0.0\t0.0\tnot-a-number\n")
            .unwrap_err();
    assert!(matches!(err, ThermalLevelBError::Number { line: 1, .. }));

    // Non-finite probe value.
    let err =
        parse_thermal_level_b_manifest(b"probe\tc\t0\t0\t0\t0\t0.0\t0.0\t0.0\tinf\n").unwrap_err();
    assert!(matches!(err, ThermalLevelBError::Number { line: 1, .. }));

    // Duplicate metadata key.
    let err = parse_thermal_level_b_manifest(
        b"case_meta\tc\tdeck_blake3\taa\ncase_meta\tc\tdeck_blake3\tbb\n",
    )
    .unwrap_err();
    assert!(matches!(err, ThermalLevelBError::Duplicate { .. }));

    // Probe indices must be dense from zero.
    let err = parse_thermal_level_b_manifest(b"probe\tc\t1\t0\t0\t0\t0.0\t0.0\t0.0\t300.0\n")
        .unwrap_err();
    assert!(matches!(
        err,
        ThermalLevelBError::ProbeOrder {
            expected: 0,
            observed: 1,
            ..
        }
    ));

    // Missing mandatory metadata.
    let err = parse_thermal_level_b_manifest(b"case_meta\tc\tself_check\tpass\n").unwrap_err();
    assert!(matches!(err, ThermalLevelBError::MissingMeta { .. }));
}

#[test]
fn failed_external_self_check_is_refused() {
    let tampered = MANIFEST.replace("\tself_check\tpass", "\tself_check\tfail");
    let references = parse_thermal_level_b_manifest(tampered.as_bytes());
    assert!(matches!(
        references.unwrap_err(),
        ThermalLevelBError::SelfCheck { .. }
    ));
}

#[test]
fn tampered_spec_echo_is_refused_bit_exactly() {
    let case = thermal_level_b_case("thermal-b-orthotropic-rotated-v1").expect("case");
    // Flip the last digit of the declared off-diagonal tensor entry: one
    // ULP-scale drift must already refuse.
    let tampered = MANIFEST.replace(
        "material.tensor.k01\t10.825317547305483",
        "material.tensor.k01\t10.825317547305484",
    );
    assert_ne!(tampered, MANIFEST, "tamper target must exist");
    let references = parse_thermal_level_b_manifest(tampered.as_bytes()).expect("parses");
    let reference = references
        .iter()
        .find(|r| r.case_id == case.id)
        .expect("case block");
    let err = verify_spec_echo(case, reference).unwrap_err();
    assert!(matches!(
        err,
        ThermalLevelBError::EchoValue { ref key, .. } if key == "material.tensor.k01"
    ));
}

#[test]
fn tampered_probe_position_is_refused_bit_exactly() {
    let case = thermal_level_b_case("thermal-b-kt-nonlinear-slab-v1").expect("case");
    // One-ULP-scale drift on the recorded x position must already refuse
    // (0.02500000000000001 parses to a different f64 than 0.025).
    let tampered = MANIFEST.replace(
        "\t0.025\t0.01\t0.01\t",
        "\t0.02500000000000001\t0.01\t0.01\t",
    );
    assert_ne!(tampered, MANIFEST, "position tamper target must exist");
    let references = parse_thermal_level_b_manifest(tampered.as_bytes()).expect("parses");
    let reference = references
        .iter()
        .find(|r| r.case_id == case.id)
        .expect("case block");
    let err = verify_probe_grid(case, reference).unwrap_err();
    assert!(matches!(err, ThermalLevelBError::ProbeGrid { .. }));
}

#[test]
fn corpus_registers_every_level_b_case_with_honest_boundaries() {
    for case in thermal_level_b_cases() {
        let dataset = corpus()
            .dataset(case.id)
            .unwrap_or_else(|| panic!("{} missing from the seeded corpus", case.id));
        assert_eq!(dataset.evidence_level(), EvidenceLevel::CrossCode);
        // Cross-code agreement can never exceed Estimated: derived-only
        // payload, no instrument, no acquisition date.
        assert_eq!(dataset.physical_claim_cap(), ColorRank::Estimated);
        let PayloadRetention::DerivedOnly { retained, .. } = dataset.raw_payload() else {
            panic!("{}: Level-B payload must be derived-only", case.id);
        };
        assert_eq!(
            retained.locator,
            "data/vv-corpus/thermal-level-b/thermal-level-b-references-v1.tsv"
        );
        let envelopes = dataset.acceptance_envelopes();
        assert_eq!(envelopes.len(), 1);
        assert_eq!(envelopes[0].metric, "probe-temperature-k");
        match envelopes[0].envelope {
            CorpusEnvelope::Tolerance { atol, rtol } => {
                assert_eq!(atol.to_bits(), case.acceptance_atol_k.to_bits());
                assert_eq!(rtol.to_bits(), 0.0_f64.to_bits());
            }
            ref other => panic!("{}: envelope must stay pinned, got {other:?}", case.id),
        }
    }
}
