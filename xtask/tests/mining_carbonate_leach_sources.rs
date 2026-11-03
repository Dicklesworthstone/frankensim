//! Actual source packs close copper-carbonate leaching and acidic gangue demand.
//! The source-defined H2CO3 component retains its formal CO2(aq)+water basis.
#![deny(unsafe_code)]

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
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const CARBONIC: &str = "source-carbonic-reference-usgs-b2131";
const MALACHITE: &str = "malachite-solid-usgs-b2131-thermochemistry";
const HPLUS: &str = "hydrogen-ion-aqueous-usgs-b2131";
const CUPRIC: &str = "copper-ii-ion-aqueous-usgs-b2131";
const H2CO3: &str = "carbonic-acid-aqueous-usgs-b2131";
const WATER: &str = "water-liquid-usgs-b2131";
const CALCITE: &str = "calcite-solid-usgs-b2131-thermochemistry";
const CALCIUM: &str = "calcium-ion-aqueous-usgs-b2131";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);

#[derive(Clone, Copy)]
enum Phase {
    Malachite,
    Calcite,
    Carbonic,
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
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-carbonate-leach-{}-{stamp}-{}",
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
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    pack
}

fn point(phase: Phase) -> QueryPoint {
    let mut p = QueryPoint::new()
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
        Phase::Malachite | Phase::Calcite => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-pure-species", 1.0),
            (
                if matches!(phase, Phase::Malachite) {
                    "source-phase-malachite"
                } else {
                    "source-phase-calcite"
                },
                1.0,
            ),
        ]),
        Phase::Carbonic => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            (CARBONIC, 1.0),
        ]),
        Phase::Ion => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-individual-ion", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
        Phase::Water => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
    }
    for (axis, value) in flags {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    p
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn value(pack: &NormalizedPack, at: &QueryPoint, name: &str, dims: Dims) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    answer.evidence.value.value
}

#[test]
fn source_rows_keep_neutral_carbonic_identity_and_malachite_reference() {
    // USGS B2131 printed pp.25 and15; retain its numbers rather than substitute
    // NBS values. NBS identifies the formal H2CO3 basis as CO2(ao)+H2O(l).
    for (slug, phase, expected) in [
        (
            MALACHITE,
            Phase::Malachite,
            [0.221116, -1_054_000.0, -890_200.0, 166.3, 2100.0, 2.5],
        ),
        (
            H2CO3,
            Phase::Carbonic,
            [0.062025, -699_700.0, -623_200.0, 184.7, 100.0, 0.9],
        ),
    ] {
        let pack = compile(slug);
        assert_eq!(pack.claims().claim_count(), 6);
        let at = point(phase);
        for ((property, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", ENERGY),
            ("standard-molar-formation-gibbs-energy", ENERGY),
            ("standard-molar-entropy", ENTROPY),
            ("reported-formation-enthalpy-uncertainty-magnitude", ENERGY),
            ("reported-standard-entropy-uncertainty-magnitude", ENTROPY),
        ]
        .into_iter()
        .zip(expected)
        {
            close(value(&pack, &at, property, dims), expected);
            for wrong in [
                QueryPoint::new(),
                point(Phase::Water),
                point(Phase::Ion),
                at.clone().with("temperature", 300.0).unwrap(),
                at.clone().with("pressure", 101_325.0).unwrap(),
                at.clone().with(REFERENCE, 0.0).unwrap(),
                at.clone()
                    .with(
                        if matches!(phase, Phase::Carbonic) {
                            CARBONIC
                        } else {
                            "source-phase-malachite"
                        },
                        0.0,
                    )
                    .unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(
                            &PropertyKey::new(property, dims),
                            &wrong,
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
        }
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new("conventional-standard-molar-entropy", ENTROPY),
                    &at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn malachite_leach_closes_copper_acid_water_heat_and_direct_gibbs() {
    let packs = [MALACHITE, HPLUS, CUPRIC, H2CO3, WATER].map(compile);
    let ids = [
        "Cu2CO3(OH)2(s)",
        "H+(aq)",
        "Cu[2+](aq)",
        "H2CO3(aq)",
        "H2O(l)",
    ]
    .map(|s| SpeciesId::new(s).unwrap());
    let points = [
        Phase::Malachite,
        Phase::Ion,
        Phase::Ion,
        Phase::Carbonic,
        Phase::Water,
    ]
    .map(point);
    let sources = std::array::from_fn::<_, 5, _>(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let gibbs_sources = std::array::from_fn::<_, 5, _>(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let reaction = MineralReaction::MalachiteCarbonicAcidLeach.build().unwrap();
    let heat = reaction.standard_enthalpy_change(&sources).unwrap();
    let gibbs = reaction.standard_gibbs_change(&gibbs_sources).unwrap();
    close(heat.molar().value(), -87_500.0);
    close(gibbs.molar().value(), -77_000.0);
    let caloric = reaction.standard_equilibrium(&sources).unwrap();
    close(caloric.molar_entropy().value(), -37.6);
    close(caloric.molar_gibbs().value(), -76_289.56);
    for term in heat.terms() {
        packs
            .iter()
            .find(|p| p.pack_id() == term.pack_id)
            .unwrap()
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
    for term in gibbs.terms() {
        packs
            .iter()
            .find(|p| p.pack_id() == term.pack_id)
            .unwrap()
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
    let masses = std::array::from_fn::<_, 5, _>(|i| {
        value(&packs[i], &points[i], "molar-mass", MolarMass::DIMS)
    });
    let extent = reaction
        .extent_from_reacted_mass(&ids[0], Mass::new(221.116), MolarMass::new(masses[0]))
        .unwrap();
    close(extent.value(), 1000.0);
    close(heat.at_extent(extent).unwrap().value(), -87_500_000.0);
    let amounts = reaction.amounts(extent).unwrap();
    for (i, expected) in [221.116, 4.032, 127.092, 62.025, 36.030]
        .into_iter()
        .enumerate()
    {
        let row = amounts.iter().find(|a| a.term.species == ids[i]).unwrap();
        close(row.amount.value() * masses[i], expected);
    }
    // Do not alter printed formula weights to conceal the 1g/kmol residual.
    close(
        1000.0 * (2.0 * masses[2] + masses[3] + 2.0 * masses[4] - masses[0] - 4.0 * masses[1]),
        -0.001,
    );
    let inventory = std::array::from_fn::<_, 5, _>(|i| {
        (ids[i].clone(), Amount::new([10.0, 6.0, 0.0, 0.0, 0.0][i]))
    });
    close(
        reaction.limiting_extent(&inventory[..2]).unwrap().value(),
        1.5,
    );
    let final_inventory = reaction
        .advance_inventory(&inventory, Amount::new(1.5))
        .unwrap();
    for (id, expected) in ids.iter().zip([8.5, 0.0, 3.0, 1.5, 3.0]) {
        close(
            final_inventory
                .iter()
                .find(|(s, _)| s == id)
                .unwrap()
                .1
                .value(),
            expected,
        );
    }
    assert!(
        reaction
            .advance_inventory(&inventory, Amount::new(2.0))
            .is_err()
    );
    assert!(reaction.standard_enthalpy_change(&sources[..4]).is_err());
    let wrong = points[3].clone().with(CARBONIC, 0.0).unwrap();
    let mut incompatible = sources;
    incompatible[3].at = &wrong;
    assert!(reaction.standard_equilibrium(&incompatible).is_err());
    let mut incompatible = gibbs_sources;
    incompatible[3].at = &wrong;
    assert!(reaction.standard_gibbs_change(&incompatible).is_err());
}

#[test]
fn acidic_calcite_endpoint_uses_two_protons_and_preserves_source_heat() {
    let packs = [CALCITE, HPLUS, CALCIUM, H2CO3].map(compile);
    let ids = ["CaCO3(s)", "H+(aq)", "Ca[2+](aq)", "H2CO3(aq)"].map(|s| SpeciesId::new(s).unwrap());
    let points = [Phase::Calcite, Phase::Ion, Phase::Ion, Phase::Carbonic].map(point);
    let sources = std::array::from_fn::<_, 4, _>(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let gibbs_sources = std::array::from_fn::<_, 4, _>(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let reaction = MineralReaction::CalciteCarbonicAcidNeutralization
        .build()
        .unwrap();
    close(
        reaction
            .standard_enthalpy_change(&sources)
            .unwrap()
            .molar()
            .value(),
        -35_300.0,
    );
    close(
        reaction
            .standard_gibbs_change(&gibbs_sources)
            .unwrap()
            .molar()
            .value(),
        -48_300.0,
    );
    let caloric = reaction.standard_equilibrium(&sources).unwrap();
    close(caloric.molar_entropy().value(), 36.79);
    close(caloric.molar_gibbs().value(), -46_268.9385);
    let inventory = [
        (ids[0].clone(), Amount::new(10.0)),
        (ids[1].clone(), Amount::new(3.0)),
    ];
    close(reaction.limiting_extent(&inventory).unwrap().value(), 1.5);
    close(
        MineralReaction::CalciteBicarbonateNeutralization
            .build()
            .unwrap()
            .limiting_extent(&inventory)
            .unwrap()
            .value(),
        3.0,
    );
}
