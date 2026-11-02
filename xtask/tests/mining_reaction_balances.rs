//! G0/G1/G3: real mineral/reagent packs feed the shared reaction-heat closure.
//! A caller's labels cannot supply missing source reference-state authority.
#![deny(unsafe_code)]

use std::{path::Path, process::Command};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, Dims, Mass, MolarMass, QuantitySpec};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{FormationEnthalpySource, ReactionBalanceError, minerals::MineralReaction},
};

const REFERENCE: &str = "source-elemental-reference-janaf";
const SLUGS: [&str; 3] = [
    "calcium-oxide-solid-nist-mining",
    "water-liquid-janaf-mining",
    "calcium-hydroxide-solid-nist-mining",
];

fn compile(slug: &str, omit_formation_axis: Option<&str>) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!(
        "frankensim-reaction-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&out).unwrap();
    let source = root.join(format!("data/matdb/seed-v1/{slug}"));
    let manifest = if let Some(axis) = omit_formation_axis {
        // Intentional mutation: the offered datum keeps its numeric value but
        // loses one source constraint. Extra query labels cannot repair that.
        std::fs::copy(source.join("manifest.tsv"), out.join("manifest.tsv")).unwrap();
        let text = std::fs::read_to_string(source.join("properties.tsv")).unwrap();
        let text = text
            .lines()
            .filter(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                !(fields.first() == Some(&"validity")
                    && fields.get(1) == Some(&"formation")
                    && fields.get(2) == Some(&axis))
            })
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        std::fs::write(out.join("properties.tsv"), text).unwrap();
        out.join("manifest.tsv")
    } else {
        source.join("manifest.tsv")
    };
    let output = out.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(manifest)
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap()
}

fn point(liquid: bool) -> QueryPoint {
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
    for (axis, value) in [
        ("source-standard-state", 1.0),
        (REFERENCE, 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            if liquid {
                "source-phase-liquid"
            } else {
                "source-phase-solid"
            },
            1.0,
        ),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn ids() -> [SpeciesId; 3] {
    ["CaO(s)", "H2O(l)", "Ca(OH)2(s)"].map(|name| SpeciesId::new(name).unwrap())
}

fn sources<'a>(
    packs: &'a [NormalizedPack; 3],
    ids: &'a [SpeciesId; 3],
    points: &'a [QueryPoint; 3],
) -> [FormationEnthalpySource<'a>; 3] {
    std::array::from_fn(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    })
}

#[test]
fn g1_g3_sourced_slaking_heat_reaches_conserved_reaction_and_mass_basis() {
    let packs = SLUGS.map(|slug| compile(slug, None));
    let ids = ids();
    let points = [point(false), point(true), point(false)];
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let heat = reaction
        .standard_enthalpy_change(&sources(&packs, &ids, &points))
        .unwrap();
    // Independently read source values: -986090 - (-635090) - (-285830).
    assert_eq!(heat.molar().value(), -65_170.0);
    assert_eq!(heat.terms().len(), 3);
    for term in heat.terms() {
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert!(!term.answer.receipt.source_hashes.is_empty());
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
    }
    let mass = packs[0]
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(
                "molar-mass",
                QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
            ),
            &points[0],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    let extent = reaction
        .extent_from_reacted_mass(
            &ids[0],
            Mass::new(56.077),
            MolarMass::new(mass.evidence.value.value),
        )
        .unwrap();
    assert!((extent.value() - 1000.0).abs() < 1e-10);
    assert!((heat.at_extent(extent).unwrap().value() + 65_170_000.0).abs() < 1e-7);
    let transfers = reaction.amounts(extent).unwrap();
    let water = transfers
        .iter()
        .find(|term| term.term.species == ids[1])
        .unwrap();
    assert_eq!(water.term.coefficient, -1);
    assert!((water.amount.value() - 1000.0).abs() < 1e-10);

    // Caller ordering never changes the canonical summation/receipt order.
    let mut reversed = sources(&packs, &ids, &points);
    reversed.reverse();
    let other = reaction.standard_enthalpy_change(&reversed).unwrap();
    assert_eq!(
        other.molar().value().to_bits(),
        heat.molar().value().to_bits()
    );
    assert_eq!(
        other.terms().iter().map(|t| &t.term).collect::<Vec<_>>(),
        heat.terms().iter().map(|t| &t.term).collect::<Vec<_>>()
    );
    assert_eq!(heat.at_extent(Amount::new(0.0)).unwrap().value(), 0.0);
}

#[test]
fn g0_source_must_declare_reference_axes_not_only_accept_caller_labels() {
    let ids = ids();
    let points = [point(false), point(true), point(false)];
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    for omitted in [
        "temperature",
        "pressure",
        "reference_convention",
        "standard_state",
        "phase",
    ] {
        let packs = [
            compile(SLUGS[0], Some(omitted)),
            compile(SLUGS[1], None),
            compile(SLUGS[2], None),
        ];
        // The generic material owner may answer an underconstrained datum.
        // Reaction thermochemistry must additionally prove its common reference.
        packs[0]
            .claims()
            .query_typed(
                &PropertyKey::with_quantity(
                    "standard-molar-formation-enthalpy",
                    QuantitySpec::dimensional(Dims([2, 1, -2, 0, 0, -1])),
                ),
                &points[0],
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert!(
            matches!(
                reaction.standard_enthalpy_change(&sources(&packs, &ids, &points)),
                Err(ReactionBalanceError::Invalid(_))
            ),
            "accepted source lacking {omitted}"
        );
    }
}

#[test]
fn g3_reaction_refuses_incomplete_duplicate_or_incompatible_source_contexts() {
    let packs = SLUGS.map(|slug| compile(slug, None));
    let ids = ids();
    let points = [point(false), point(true), point(false)];
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let mut offered = sources(&packs, &ids, &points);
    assert!(reaction.standard_enthalpy_change(&offered[..2]).is_err());
    offered[2].species = &ids[0];
    assert!(matches!(
        reaction.standard_enthalpy_change(&offered),
        Err(ReactionBalanceError::Species { .. })
    ));
    let mut offered = sources(&packs, &ids, &points);
    offered[0].reference_axis = "source-phase-solid";
    assert!(reaction.standard_enthalpy_change(&offered).is_err());
    for changed in [
        point(true),
        point(false).with("temperature", 300.0).unwrap(),
        point(false).with("pressure", 101_325.0).unwrap(),
        point(false).with(REFERENCE, 0.0).unwrap(),
        point(false).with("source-pure-species", 0.0).unwrap(),
    ] {
        let changed_points = [changed, point(true), point(false)];
        assert!(
            reaction
                .standard_enthalpy_change(&sources(&packs, &ids, &changed_points))
                .is_err()
        );
    }
}
