//! G0/G1/G3: real USGS packs close dolomite acid demand and reference-state heat.
//! Supplied reacted mass and an inventory ceiling do not predict dissolution rate.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_evidence::NumericalKind;
use fs_matdb::{
    ClaimSelection, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
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
const SLUGS: [&str; 5] = [
    "dolomite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
    "calcium-ion-aqueous-usgs-b2131",
    "magnesium-ion-aqueous-usgs-b2131",
    "bicarbonate-ion-aqueous-usgs-b2131",
];
const HF: &str = "standard-molar-formation-enthalpy";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);

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
        "frankensim-mining-dolomite-{}-{nonce}-{}",
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

fn point(dolomite: bool, temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
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
    if dolomite {
        flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-dolomite", 1.0),
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
fn g0_g3_dolomite_and_magnesium_preserve_source_values_and_reference_states() {
    // Independent literal values: USGS B2131 printed pp.25 (dolomite), 8 (Mg2+).
    for (slug, solid, mass, hf, entropy, hf_uncertainty, s_uncertainty) in [
        (SLUGS[0], true, 0.184401, -2_324_500.0, 155.2, 1500.0, 0.3),
        (SLUGS[3], false, 0.024305, -467_000.0, -137.0, 600.0, 4.0),
    ] {
        let pack = compile(slug);
        assert_eq!(pack.claims().claim_count(), if solid { 6 } else { 5 });
        let query = point(solid, 298.15);
        let (entropy_key, uncertainty_key) = if solid {
            (
                "standard-molar-entropy",
                "reported-standard-entropy-uncertainty-magnitude",
            )
        } else {
            (
                "conventional-standard-molar-entropy",
                "reported-conventional-entropy-uncertainty-magnitude",
            )
        };
        for (property, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, mass),
            (HF, MOLAR_ENERGY, hf),
            (entropy_key, MOLAR_ENTROPY, entropy),
            (
                "reported-formation-enthalpy-uncertainty-magnitude",
                MOLAR_ENERGY,
                hf_uncertainty,
            ),
            (uncertainty_key, MOLAR_ENTROPY, s_uncertainty),
        ] {
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
            assert!(!answer.receipt.source_hashes.is_empty());
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            for (axis, coordinate) in query.axes() {
                let upper = if solid && property == "molar-mass" && axis == "temperature" {
                    500.0
                } else {
                    *coordinate
                };
                assert_eq!(claim.validity.bound(axis), Some((*coordinate, upper)));
                assert_eq!(
                    claim.validity.axis_quantities().get(axis),
                    query.axis_quantities().get(axis)
                );
            }
            for wrong in [
                QueryPoint::new(),
                point(!solid, 298.15),
                query.clone().with("pressure", 101_325.0).unwrap(),
                query.clone().with(REFERENCE, 0.0).unwrap(),
                query
                    .clone()
                    .with("source-specimen-qualified", 1.0)
                    .unwrap(),
                query
                    .clone()
                    .with(
                        if solid {
                            "source-phase-dolomite"
                        } else {
                            IONIC
                        },
                        0.0,
                    )
                    .unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{slug} {property}: admitted wrong state {wrong:?}"
                );
            }
            if property != "molar-mass" {
                assert!(
                    pack.claims()
                        .query_typed(&key, &point(solid, 300.0), SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
            }
        }
        if !solid {
            // Negative individual-ion entropy cannot be consumed as absolute.
            assert!(
                pack.claims()
                    .query_typed(
                        &PropertyKey::new("standard-molar-entropy", MOLAR_ENTROPY),
                        &query,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
            for axis in [AQUEOUS, "source-individual-ion"] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &PropertyKey::new(HF, MOLAR_ENERGY),
                            &query.clone().with(axis, 0.0).unwrap(),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
        }
    }
}

#[test]
fn g1_dolomite_heat_capacity_matches_points_and_independent_heat_content() {
    let pack = compile(SLUGS[0]);
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let line = pack
        .claims()
        .query_line_typed(
            &cp,
            "temperature",
            &point(true, 298.15),
            &point(true, 500.0),
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )
        .unwrap();
    assert_eq!(line.knots.len(), 4);
    for ((temperature, answer), (source_t, source_cp)) in line.knots.iter().zip([
        (298.15, 157.51),
        (300.0, 158.06),
        (400.0, 183.54),
        (500.0, 201.89),
    ]) {
        assert_eq!(*temperature, source_t);
        close(answer.evidence.value.value * 0.184401, source_cp);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        pack.claims().verify_receipt(&answer.receipt).unwrap();
    }
    let heating = pack
        .claims()
        .query_line_typed(
            &cp,
            "temperature",
            &point(true, 300.0),
            &point(true, 400.0),
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )
        .unwrap();
    let heat_j: f64 = 2.0
        * heating
            .knots
            .windows(2)
            .map(|w| {
                0.5 * (w[1].0 - w[0].0)
                    * (w[0].1.evidence.value.value + w[1].1.evidence.value.value)
            })
            .sum::<f64>();
    // USGS p.269 independently prints (H(T)-H(298.15))/T, not delta H.
    // The 0.5% allowance compares the authored coarse linearization with this
    // separate source column; it is not a source uncertainty certificate.
    let source_heat_j = 2.0 / 0.184401 * (400.0 * 43.61 - 300.0 * 0.97);
    assert!(heat_j > 0.0);
    assert!((heat_j - source_heat_j).abs() / source_heat_j < 0.005);
    for temperature in [298.0, 501.0] {
        assert!(
            pack.claims()
                .query_typed(
                    &cp,
                    &point(true, temperature),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn g1_real_packs_close_dolomite_heat_mass_and_acid_limited_inventory() {
    let packs = SLUGS.map(compile);
    let ids = [
        "CaMg(CO3)2(s)",
        "H+(aq)",
        "Ca[2+](aq)",
        "Mg[2+](aq)",
        "HCO3[-](aq)",
    ]
    .map(|name| SpeciesId::new(name).unwrap());
    let points: [_; 5] = std::array::from_fn(|i| point(i == 0, 298.15));
    let source_queries = || {
        [0, 1, 2, 3, 4].map(|i| FormationEnthalpySource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
    };
    let sources = source_queries();
    let reaction = MineralReaction::DolomiteBicarbonateNeutralization
        .build()
        .unwrap();
    let heat = reaction.standard_enthalpy_change(&sources).unwrap();
    // Independent literal rows: -543000-467000+2*(-689900)-(-2324500)-2*0.
    assert_eq!(heat.molar().value(), -65_300.0);
    assert_eq!(heat.terms().len(), 5);
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
    let masses: [_; 5] = std::array::from_fn(|i| {
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
        .extent_from_reacted_mass(&ids[0], Mass::new(184.401), MolarMass::new(masses[0]))
        .unwrap();
    close(extent.value(), 1000.0);
    close(heat.at_extent(extent).unwrap().value(), -65_300_000.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (i, expected_mass) in [
        (0, 184.401),
        (1, 2.016),
        (2, 40.078),
        (3, 24.305),
        (4, 122.034),
    ] {
        let amount = amounts
            .iter()
            .find(|row| row.term.species == ids[i])
            .unwrap();
        close(amount.amount.value() * masses[i], expected_mass);
    }
    close(
        masses[0] + 2.0 * masses[1],
        masses[2] + masses[3] + 2.0 * masses[4],
    );

    // Two acid equivalents per reacted dolomite, with product inventories zero.
    let inventory: [_; 5] =
        std::array::from_fn(|i| (ids[i].clone(), Amount::new([10.0, 6.0, 0.0, 0.0, 0.0][i])));
    let limit = reaction.limiting_extent(&inventory[..2]).unwrap();
    assert_eq!(limit.value(), 3.0);
    assert_eq!(heat.at_extent(limit).unwrap().value(), -195_900.0);
    let after = reaction.advance_inventory(&inventory, limit).unwrap();
    for (i, expected) in [(0, 7.0), (1, 0.0), (2, 3.0), (3, 3.0), (4, 6.0)] {
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
            .advance_inventory(&inventory, Amount::new(4.0))
            .is_err()
    );

    // Mg2+ is required; its elemental and ionic conventions cannot be relabeled.
    let without_magnesium = [0, 1, 2, 4].map(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    assert!(
        reaction
            .standard_enthalpy_change(&without_magnesium)
            .is_err()
    );
    let mut other_reference = source_queries();
    other_reference[3].reference_axis = "source-elemental-reference-janaf";
    assert!(reaction.standard_enthalpy_change(&other_reference).is_err());
    let wrong_ionic = point(false, 298.15).with(IONIC, 0.0).unwrap();
    let mut incompatible = source_queries();
    incompatible[3].at = &wrong_ionic;
    assert!(reaction.standard_enthalpy_change(&incompatible).is_err());
}
