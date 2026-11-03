//! G0/G1/G3: compiled USGS references close copper-oxide acid demand and heat.
//! A reacted amount or inventory ceiling does not predict recovery versus time.
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
    stoichiometry::{
        FormationEnthalpySource, minerals::MineralReaction, thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const AQUEOUS: &str = "source-aqueous-unit-activity-infinite-dilution";
const IONIC: &str = "source-ionic-reference-hydrogen-zero";
const TENORITE: &str = "source-phase-tenorite";
const SLUGS: [&str; 4] = [
    "copper-ii-oxide-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
    "copper-ii-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
];
const HF: &str = "standard-molar-formation-enthalpy";
const GF: &str = "standard-molar-formation-gibbs-energy";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);

#[derive(Clone, Copy)]
enum Phase {
    Tenorite,
    Ion,
    Water,
}

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
        "frankensim-mining-copper-leach-{}-{nonce}-{}",
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

fn point(phase: Phase) -> QueryPoint {
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
    match phase {
        Phase::Tenorite => flags.extend([
            ("source-phase-solid", 1.0),
            (TENORITE, 1.0),
            ("source-pure-species", 1.0),
        ]),
        Phase::Ion => flags.extend([
            ("source-phase-aqueous", 1.0),
            (AQUEOUS, 1.0),
            (IONIC, 1.0),
            ("source-individual-ion", 1.0),
        ]),
        Phase::Water => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
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
fn g0_g3_copper_sources_preserve_literal_values_and_reference_states() {
    // USGS B2131 printed pp.16 (tenorite), 6 (Cu2+). Retain this compilation's
    // +64.9 kJ/mol and -98.0 J/(mol K) for Cu2+, including the signs.
    for (slug, phase, mass, hf, gf, entropy, h_unc, s_unc) in [
        (
            SLUGS[0],
            Phase::Tenorite,
            0.079545,
            -156_100.0,
            -128_300.0,
            42.6,
            2000.0,
            0.2,
        ),
        (
            SLUGS[2],
            Phase::Ion,
            0.063546,
            64_900.0,
            65_100.0,
            -98.0,
            1000.0,
            4.0,
        ),
    ] {
        let pack = compile(slug);
        assert_eq!(pack.claims().claim_count(), 6);
        let query = point(phase);
        let ionic = matches!(phase, Phase::Ion);
        let (entropy_key, uncertainty_key) = if ionic {
            (
                "conventional-standard-molar-entropy",
                "reported-conventional-entropy-uncertainty-magnitude",
            )
        } else {
            (
                "standard-molar-entropy",
                "reported-standard-entropy-uncertainty-magnitude",
            )
        };
        for (property, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, mass),
            (HF, MOLAR_ENERGY, hf),
            (GF, MOLAR_ENERGY, gf),
            (entropy_key, MOLAR_ENTROPY, entropy),
            (
                "reported-formation-enthalpy-uncertainty-magnitude",
                MOLAR_ENERGY,
                h_unc,
            ),
            (uncertainty_key, MOLAR_ENTROPY, s_unc),
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
                assert_eq!(claim.validity.bound(axis), Some((*coordinate, *coordinate)));
                assert_eq!(
                    claim.validity.axis_quantities().get(axis),
                    query.axis_quantities().get(axis)
                );
            }
            for wrong in [
                QueryPoint::new(),
                point(Phase::Water),
                query.clone().with("temperature", 300.0).unwrap(),
                query.clone().with("pressure", 101_325.0).unwrap(),
                query.clone().with(REFERENCE, 0.0).unwrap(),
                query
                    .clone()
                    .with("source-specimen-qualified", 1.0)
                    .unwrap(),
                query
                    .clone()
                    .with(if ionic { IONIC } else { TENORITE }, 0.0)
                    .unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{slug} {property}: admitted wrong state {wrong:?}"
                );
            }
        }
        if ionic {
            // Conventional individual-ion entropy is a distinct property.
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
fn g1_real_packs_close_copper_leach_heat_gibbs_mass_and_acid_inventory() {
    let packs = SLUGS.map(compile);
    let ids =
        ["CuO(s)", "H+(aq)", "Cu[2+](aq)", "H2O(l)"].map(|name| SpeciesId::new(name).unwrap());
    let points = [
        point(Phase::Tenorite),
        point(Phase::Ion),
        point(Phase::Ion),
        point(Phase::Water),
    ];
    let source_queries = || {
        [0, 1, 2, 3].map(|i| FormationEnthalpySource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
    };
    let sources = source_queries();
    let reaction = MineralReaction::CupricOxideAcidLeach.build().unwrap();
    let heat = reaction.standard_enthalpy_change(&sources).unwrap();
    // Literal USGS rows, J/mol: 64900-285800-(-156100)-2*0.
    close(heat.molar().value(), -64_800.0);
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
    let gibbs_sources = [0, 1, 2, 3].map(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let gibbs = reaction.standard_gibbs_change(&gibbs_sources).unwrap();
    // Independently tabulated formation Gibbs values, not Hf minus T*S_species.
    close(gibbs.molar().value(), -43_700.0);
    assert_eq!(gibbs.terms().len(), 4);
    for term in gibbs.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
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
        .extent_from_reacted_mass(&ids[0], Mass::new(79.545), MolarMass::new(masses[0]))
        .unwrap();
    close(extent.value(), 1000.0);
    close(heat.at_extent(extent).unwrap().value(), -64_800_000.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (i, expected_mass) in [(0, 79.545), (1, 2.016), (2, 63.546), (3, 18.015)] {
        let amount = amounts
            .iter()
            .find(|row| row.term.species == ids[i])
            .unwrap();
        close(amount.amount.value() * masses[i], expected_mass);
    }
    close(masses[0] + 2.0 * masses[1], masses[2] + masses[3]);

    // Two acid equivalents per reacted CuO, with product inventories zero.
    let inventory: [_; 4] =
        std::array::from_fn(|i| (ids[i].clone(), Amount::new([10.0, 6.0, 0.0, 0.0][i])));
    let limit = reaction.limiting_extent(&inventory[..2]).unwrap();
    assert_eq!(limit.value(), 3.0);
    close(heat.at_extent(limit).unwrap().value(), -194_400.0);
    let after = reaction.advance_inventory(&inventory, limit).unwrap();
    for (i, expected) in [(0, 7.0), (1, 0.0), (2, 3.0), (3, 3.0)] {
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

    assert!(reaction.standard_enthalpy_change(&sources[..3]).is_err());
    assert!(reaction.standard_gibbs_change(&gibbs_sources[..3]).is_err());
    let mut other_reference = source_queries();
    other_reference[0].reference_axis = "source-elemental-reference-janaf";
    assert!(reaction.standard_enthalpy_change(&other_reference).is_err());
    let wrong_ionic = point(Phase::Ion).with(IONIC, 0.0).unwrap();
    let mut incompatible = source_queries();
    incompatible[2].at = &wrong_ionic;
    assert!(reaction.standard_enthalpy_change(&incompatible).is_err());
    let mut incompatible_gibbs = gibbs_sources;
    incompatible_gibbs[2].at = &wrong_ionic;
    assert!(reaction.standard_gibbs_change(&incompatible_gibbs).is_err());
}
