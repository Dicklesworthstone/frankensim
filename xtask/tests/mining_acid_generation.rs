//! G0/G1/G3: compatible USGS sources close aqueous pyrite acid generation.
//! A declared ferrous endpoint and inventory ceiling do not predict oxidation rate.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

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
const HF: &str = "standard-molar-formation-enthalpy";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 6] = [
    "pyrite-solid-usgs-b2131-thermochemistry",
    "oxygen-gas-usgs-b2131",
    "water-liquid-usgs-b2131",
    "iron-ii-ion-aqueous-usgs-b2131",
    "sulfate-ion-aqueous-usgs-b2131",
    "hydrogen-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 6] = [
    "FeS2(s)",
    "O2(g)",
    "H2O(l)",
    "Fe[2+](aq)",
    "SO4[2-](aq)",
    "H+(aq)",
];

#[derive(Clone, Copy)]
enum State {
    Pyrite,
    Gas,
    Liquid,
    Ion,
}
const STATES: [State; 6] = [
    State::Pyrite,
    State::Gas,
    State::Liquid,
    State::Ion,
    State::Ion,
    State::Ion,
];

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
        "frankensim-acid-generation-{}-{nonce}-{}",
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

fn point(state: State) -> QueryPoint {
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
    match state {
        State::Ion => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
            ("source-individual-ion", 1.0),
        ]),
        State::Pyrite => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-pyrite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        State::Gas => flags.extend([("source-phase-gas", 1.0), ("source-pure-species", 1.0)]),
        State::Liquid => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
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
fn g0_g3_actual_usgs_rows_preserve_values_phases_and_ionic_conventions() {
    // Independent printed USGS B2131 rows: pp.12,8,16,7,18, respectively.
    for (i, mass, hf, entropy, h_uncertainty, s_uncertainty) in [
        (0, 0.119979, -171_500.0, 52.9, Some(1700.0), 0.1),
        (1, 0.031999, 0.0, 205.15, None, 0.02),
        (2, 0.018015, -285_800.0, 70.0, Some(100.0), 0.1),
        (3, 0.055847, -91_100.0, -107.1, Some(3000.0), 2.0),
        (4, 0.096064, -909_300.0, 18.5, Some(400.0), 0.4),
    ] {
        let pack = compile(SLUGS[i]);
        let at = point(STATES[i]);
        let ionic = matches!(STATES[i], State::Ion);
        let entropy_name = if ionic {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        let entropy_uncertainty_name = if ionic {
            "reported-conventional-entropy-uncertainty-magnitude"
        } else {
            "reported-standard-entropy-uncertainty-magnitude"
        };
        let mut values = vec![
            ("molar-mass", MolarMass::DIMS, mass),
            (HF, MOLAR_ENERGY, hf),
            (entropy_name, MOLAR_ENTROPY, entropy),
            (entropy_uncertainty_name, MOLAR_ENTROPY, s_uncertainty),
        ];
        if let Some(value) = h_uncertainty {
            values.push((
                "reported-formation-enthalpy-uncertainty-magnitude",
                MOLAR_ENERGY,
                value,
            ));
        }
        // Source profiles also carry a separate standard formation Gibbs claim.
        assert_eq!(pack.claims().claim_count(), values.len() + 1);
        for (property, dims, expected) in values {
            let key = PropertyKey::new(property, dims);
            let answer = pack
                .claims()
                .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            assert_eq!(claim.validity.bound(REFERENCE), Some((1.0, 1.0)));
            if ionic {
                assert_eq!(claim.validity.bound("source-pure-species"), None);
            }
            for wrong in [
                QueryPoint::new(),
                at.clone().with("temperature", 300.0).unwrap(),
                at.clone().with("pressure", 101_325.0).unwrap(),
                at.clone().with(REFERENCE, 0.0).unwrap(),
                at.clone().with("source-specimen-qualified", 1.0).unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
            }
        }
        let incompatible_entropy = if ionic {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(incompatible_entropy, MOLAR_ENTROPY),
                    &at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn g1_real_packs_close_pyrite_oxygen_water_heat_and_acid_inventory() {
    let packs = SLUGS.map(compile);
    let ids = SPECIES.map(|name| SpeciesId::new(name).unwrap());
    let points = STATES.map(point);
    let sources: [_; 6] = std::array::from_fn(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let reaction = MineralReaction::PyriteOxidationFerrous.build().unwrap();
    let heat = reaction.standard_enthalpy_change(&sources).unwrap();
    // 2*(-91100)+4*(-909300)-2*(-171500)-2*(-285800); oxygen/H+ remain sources.
    assert_eq!(heat.molar().value(), -2_904_800.0);
    assert_eq!(heat.terms().len(), 6);
    for term in heat.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        assert!(!term.answer.receipt.source_hashes.is_empty());
    }
    let masses: [_; 6] = std::array::from_fn(|i| {
        packs[i]
            .claims()
            .query_typed(
                &PropertyKey::new("molar-mass", MolarMass::DIMS),
                &points[i],
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap()
            .evidence
            .value
            .value
    });
    let extent = reaction
        .extent_from_reacted_mass(&ids[0], Mass::new(119.979), MolarMass::new(masses[0]))
        .unwrap();
    close(extent.value(), 500.0);
    close(heat.at_extent(extent).unwrap().value(), -1_452_400_000.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (i, expected_mass) in [119.979, 111.9965, 18.015, 55.847, 192.128, 2.016]
        .into_iter()
        .enumerate()
    {
        let amount = amounts
            .iter()
            .find(|row| row.term.species == ids[i])
            .unwrap();
        close(amount.amount.value() * masses[i], expected_mass);
    }
    // The original rounded formula weights give +0.0005 kg, not an exact zero.
    close(
        1000.0
            * (masses[3] + 2.0 * masses[4] + 2.0 * masses[5]
                - masses[0]
                - 3.5 * masses[1]
                - masses[2]),
        0.0005,
    );

    let initial = [10.0, 14.0, 20.0, 0.0, 0.0, 1.0];
    let inventory: [_; 6] = std::array::from_fn(|i| (ids[i].clone(), Amount::new(initial[i])));
    let limit = reaction.limiting_extent(&inventory[..3]).unwrap();
    assert_eq!(limit.value(), 2.0, "oxygen limits the supplied inventory");
    assert_eq!(heat.at_extent(limit).unwrap().value(), -5_809_600.0);
    let after = reaction.advance_inventory(&inventory, limit).unwrap();
    for (i, expected) in [6.0, 0.0, 16.0, 4.0, 8.0, 9.0].into_iter().enumerate() {
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
            .advance_inventory(&inventory, Amount::new(3.0))
            .is_err()
    );
}

#[test]
fn g3_reference_or_phase_substitution_cannot_supply_a_complete_reaction() {
    let packs = SLUGS.map(compile);
    let ids = SPECIES.map(|name| SpeciesId::new(name).unwrap());
    let points = STATES.map(point);
    let make_sources = || {
        std::array::from_fn::<_, 6, _>(|i| FormationEnthalpySource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
    };
    let reaction = MineralReaction::PyriteOxidationFerrous.build().unwrap();
    assert!(
        reaction
            .standard_enthalpy_change(&make_sources()[..5])
            .is_err(),
        "zero-reference H+ is still required"
    );
    let janaf = compile("pyrite-solid-janaf-thermochemistry");
    let both_labels = points[0]
        .clone()
        .with_quantity(
            "source-elemental-reference-janaf",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    let mut mislabeled = make_sources();
    mislabeled[0].pack = &janaf;
    mislabeled[0].at = &both_labels;
    assert!(
        reaction.standard_enthalpy_change(&mislabeled).is_err(),
        "a caller's USGS label cannot relabel the selected JANAF source"
    );
    for (i, axis) in [
        (0, "source-phase-pyrite"),
        (2, "source-phase-liquid"),
        (3, "source-ionic-reference-hydrogen-zero"),
        (4, "source-individual-ion"),
    ] {
        let wrong = points[i].clone().with(axis, 0.0).unwrap();
        let mut sources = make_sources();
        sources[i].at = &wrong;
        assert!(reaction.standard_enthalpy_change(&sources).is_err());
    }
    let mut reversed = make_sources();
    reversed.reverse();
    assert_eq!(
        reaction
            .standard_enthalpy_change(&reversed)
            .unwrap()
            .molar()
            .value()
            .to_bits(),
        reaction
            .standard_enthalpy_change(&make_sources())
            .unwrap()
            .molar()
            .value()
            .to_bits()
    );
}
