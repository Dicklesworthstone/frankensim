//! G0/G1/G3: USGS ionic references close a source-backed calcite acid balance.
//! Reference-state heat and an inventory ceiling do not predict dissolution rate.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_evidence::NumericalKind;
use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{FormationEnthalpySource, minerals::MineralReaction},
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const AQUEOUS: &str = "source-aqueous-unit-activity-infinite-dilution";
const IONIC: &str = "source-ionic-reference-hydrogen-zero";
const SLUGS: [&str; 3] = [
    "hydrogen-ion-aqueous-usgs-b2131",
    "calcium-ion-aqueous-usgs-b2131",
    "bicarbonate-ion-aqueous-usgs-b2131",
];
const HF: &str = "standard-molar-formation-enthalpy";
const S: &str = "conventional-standard-molar-entropy";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-mining-aqueous-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .expect("run the actual material compiler");
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes, "{slug}: canonical round trip");
    assert_eq!(pack.pack_id(), slug);
    pack
}

fn point(calcite: bool) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![
        ("source-standard-state", 1.0),
        (REFERENCE, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if calcite {
        flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-calcite", 1.0),
            ("source-pure-species", 1.0),
        ]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            (AQUEOUS, 1.0),
            (IONIC, 1.0),
            ("source-individual-ion", 1.0),
        ]);
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}",
    );
}

#[test]
fn g0_g3_ionic_source_values_preserve_convention_and_reported_uncertainty() {
    // USGS B2131 printed pp.7 (H+), 6 (Ca2+), 15 (HCO3-).
    // Negative Ca2+ entropy is conventional, not an absolute ionic entropy.
    for (slug, mass, hf, entropy, reported) in [
        (SLUGS[0], 0.001008, 0.0, 0.0, None),
        (SLUGS[1], 0.040078, -543_000.0, -56.2, Some((1000.0, 1.0))),
        (SLUGS[2], 0.061017, -689_900.0, 98.4, Some((200.0, 0.5))),
    ] {
        let pack = compile(slug);
        assert_eq!(
            pack.claims().claim_count(),
            if reported.is_some() { 5 } else { 3 }
        );
        let query = point(false);
        let mut values = vec![
            ("molar-mass", MolarMass::DIMS, mass),
            (HF, MOLAR_ENERGY, hf),
            (S, MOLAR_ENTROPY, entropy),
        ];
        if let Some((h, s)) = reported {
            values.extend([
                (
                    "reported-formation-enthalpy-uncertainty-magnitude",
                    MOLAR_ENERGY,
                    h,
                ),
                (
                    "reported-conventional-entropy-uncertainty-magnitude",
                    MOLAR_ENTROPY,
                    s,
                ),
            ]);
        }
        for (property, dims, expected) in values {
            let key = PropertyKey::new(property, dims);
            let answer = pack
                .claims()
                .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert_eq!(answer.evidence.numerical.kind, NumericalKind::NoClaim);
            assert!(answer.receipt.observation_backed);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            for (axis, coordinate) in query.axes() {
                assert_eq!(claim.validity.bound(axis), Some((*coordinate, *coordinate)));
                assert_eq!(
                    claim.validity.axis_quantities().get(axis),
                    query.axis_quantities().get(axis),
                );
            }
            assert_eq!(claim.validity.bound("source-pure-species"), None);
            for wrong in [
                QueryPoint::new(),
                point(true),
                query.clone().with("temperature", 300.0).unwrap(),
                query.clone().with("pressure", 101_325.0).unwrap(),
                query.clone().with(REFERENCE, 0.0).unwrap(),
                query.clone().with(AQUEOUS, 0.0).unwrap(),
                query.clone().with(IONIC, 0.0).unwrap(),
                query.clone().with("source-individual-ion", 0.0).unwrap(),
                query
                    .clone()
                    .with("source-specimen-qualified", 1.0)
                    .unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{slug} {property}: admitted wrong state {wrong:?}",
                );
            }
        }
        // A client cannot accidentally consume conventional entropy as absolute.
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new("standard-molar-entropy", MOLAR_ENTROPY),
                    &query,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn g1_real_packs_close_calcite_bicarbonate_heat_and_acid_limited_inventory() {
    let packs = [
        "calcite-solid-usgs-b2131-thermochemistry",
        SLUGS[0],
        SLUGS[1],
        SLUGS[2],
    ]
    .map(compile);
    let ids = ["CaCO3(s)", "H+(aq)", "Ca[2+](aq)", "HCO3[-](aq)"]
        .map(|name| SpeciesId::new(name).unwrap());
    let points = [point(true), point(false), point(false), point(false)];
    let source_queries = || {
        [0, 1, 2, 3].map(|i| FormationEnthalpySource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
    };
    let sources = source_queries();
    let reaction = MineralReaction::CalciteBicarbonateNeutralization
        .build()
        .unwrap();
    let heat = reaction.standard_enthalpy_change(&sources).unwrap();
    // Independent literal USGS rows in J/mol: -543000-689900-(-1207400)-0.
    assert_eq!(heat.molar().value(), -25_500.0);
    assert_eq!(heat.terms().len(), 4);
    for term in heat.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        assert!(!term.answer.receipt.source_hashes.is_empty());
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    let masses: [_; 4] = std::array::from_fn(|i| {
        let answer = packs[i]
            .claims()
            .query_typed(
                &PropertyKey::new("molar-mass", MolarMass::DIMS),
                &points[i],
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        packs[i].claims().verify_receipt(&answer.receipt).unwrap();
        answer.evidence.value.value
    });
    let extent = reaction
        .extent_from_reacted_mass(&ids[0], Mass::new(100.087), MolarMass::new(masses[0]))
        .unwrap();
    close(extent.value(), 1000.0);
    close(heat.at_extent(extent).unwrap().value(), -25_500_000.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (i, expected_mass) in [(0, 100.087), (1, 1.008), (2, 40.078), (3, 61.017)] {
        let amount = amounts
            .iter()
            .find(|row| row.term.species == ids[i])
            .unwrap();
        close(amount.amount.value() * masses[i], expected_mass);
    }
    close(masses[0] + masses[1], masses[2] + masses[3]);

    // This is a stoichiometric ceiling from supplied amounts, not kinetics or
    // a claim that the acid-limited extent will be reached by a particular time.
    let inventory = [
        (ids[0].clone(), Amount::new(10.0)),
        (ids[1].clone(), Amount::new(6.0)),
        (ids[2].clone(), Amount::new(0.0)),
        (ids[3].clone(), Amount::new(0.0)),
    ];
    let limit = reaction.limiting_extent(&inventory[..2]).unwrap();
    assert_eq!(limit.value(), 6.0);
    assert_eq!(heat.at_extent(limit).unwrap().value(), -153_000.0);
    let after = reaction.advance_inventory(&inventory, limit).unwrap();
    for (i, expected) in [(0, 4.0), (1, 0.0), (2, 6.0), (3, 6.0)] {
        assert_eq!(
            after
                .iter()
                .find(|(id, _)| *id == ids[i])
                .unwrap()
                .1
                .value(),
            expected
        );
    }
    assert!(
        reaction
            .advance_inventory(&inventory, Amount::new(7.0))
            .is_err()
    );

    // Zero-convention H+ remains a required source, and its phase/reference
    // cannot be replaced by a JANAF or finite-concentration caller label.
    assert!(reaction.standard_enthalpy_change(&sources[..3]).is_err());
    let mut other_reference = source_queries();
    other_reference[1].reference_axis = "source-elemental-reference-janaf";
    assert!(reaction.standard_enthalpy_change(&other_reference).is_err());
    let wrong_ionic = point(false).with(IONIC, 0.0).unwrap();
    let mut incompatible = source_queries();
    incompatible[1].at = &wrong_ionic;
    assert!(reaction.standard_enthalpy_change(&incompatible).is_err());
    // Bicarbonate data do not silently become the different CO2(aq) endpoint.
    assert!(
        MineralReaction::CalciteAcidNeutralization
            .build()
            .unwrap()
            .standard_enthalpy_change(&sources)
            .is_err()
    );
}
