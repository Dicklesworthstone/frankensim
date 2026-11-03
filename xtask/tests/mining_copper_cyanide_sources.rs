//! G0/G3 original corrected NBS rows through conserved Cu(I)-cyanide chemistry.
//! Formation constants and ligand demand are source-backed standard-state
//! calculations; these tests do not assert a coupled solution or ore recovery.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource, ReactionBalanceError,
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const RT: f64 = 8.314_462_618_153_24 * 298.15;

struct Row {
    id: &'static str,
    pack: &'static str,
    // Explicit Cu/C/N atom counts and charge, not a label-parsing shortcut.
    atoms: [i128; 3],
    charge: i128,
    mass: f64,             // original g/mol
    enthalpy: Option<f64>, // original kJ/mol; None means unavailable, not zero
    gibbs: f64,            // original kJ/mol
    entropy: Option<f64>,  // original J/(mol K)
}
// NIST corrected CSV (SHA256 verified), Table36:Cu, prepared1968; CN is Table23:C.
const ROWS: [Row; 7] = [
    Row {
        id: "Cu[+](aq)",
        pack: "copper-i-ion-aqueous-nbs1989",
        atoms: [1, 0, 0],
        charge: 1,
        mass: 63.5400,
        enthalpy: Some(71.67),
        gibbs: 49.98,
        entropy: Some(40.6),
    },
    Row {
        id: "Cu[2+](aq)",
        pack: "copper-ii-ion-aqueous-nbs1989",
        atoms: [1, 0, 0],
        charge: 2,
        mass: 63.5400,
        enthalpy: Some(64.77),
        gibbs: 65.49,
        entropy: Some(-99.6),
    },
    Row {
        id: "CuCN(s)",
        pack: "copper-i-cyanide-solid-nbs1989",
        atoms: [1, 1, 1],
        charge: 0,
        mass: 89.5579,
        enthalpy: Some(96.2),
        gibbs: 111.3,
        entropy: Some(84.5),
    },
    Row {
        id: "Cu(CN)2[-](aq)",
        pack: "dicyanocuprate-i-ion-aqueous-nbs1989",
        atoms: [1, 2, 2],
        charge: -1,
        mass: 115.5758,
        enthalpy: None,
        gibbs: 257.8,
        entropy: None,
    },
    Row {
        id: "Cu(CN)3[2-](aq)",
        pack: "tricyanocuprate-i-ion-aqueous-nbs1989",
        atoms: [1, 3, 3],
        charge: -2,
        mass: 141.5937,
        enthalpy: None,
        gibbs: 403.8,
        entropy: None,
    },
    Row {
        id: "Cu(CN)4[3-](aq)",
        pack: "tetracyanocuprate-i-ion-aqueous-nbs1989",
        atoms: [1, 4, 4],
        charge: -3,
        mass: 167.6116,
        enthalpy: None,
        gibbs: 566.6,
        entropy: None,
    },
    Row {
        id: "CN[-](aq)",
        pack: "cyanide-ion-aqueous-nbs1989",
        atoms: [0, 1, 1],
        charge: -1,
        mass: 26.0179,
        enthalpy: Some(150.6),
        gibbs: 172.4,
        entropy: Some(94.1),
    },
];

struct Source {
    id: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
}
fn workspace() -> PathBuf {
    std::env::var_os("CARGO_WORKSPACE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .into()
        })
}
fn point(row: &Row) -> QueryPoint {
    let mut at = QueryPoint::new()
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
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if row.charge == 0 {
        flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
            ("source-individual-ion", 1.0),
        ]);
    }
    for (name, value) in flags {
        at = at
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}
fn sources() -> &'static [Source; 7] {
    static SOURCES: OnceLock<[Source; 7]> = OnceLock::new();
    SOURCES.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-copper-cyanide-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        std::array::from_fn(|i| {
            let row = &ROWS[i];
            let target = out.join(format!("{}.fsmatpk", row.pack));
            let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                .args(["matdb-pack", "--manifest"])
                .arg(workspace().join(format!("data/matdb/seed-v1/{}/manifest.tsv", row.pack)))
                .arg("--out")
                .arg(&target)
                .env("CARGO_WORKSPACE_DIR", workspace())
                .current_dir(workspace())
                .output()
                .unwrap();
            assert!(
                run.status.success(),
                "{}: {}",
                row.pack,
                String::from_utf8_lossy(&run.stderr)
            );
            let bytes = std::fs::read(target).unwrap();
            let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
            assert_eq!(pack.pack_id(), row.pack);
            assert_eq!(pack.to_bytes(), bytes);
            Source {
                id: SpeciesId::new(row.id).unwrap(),
                pack,
                at: point(row),
            }
        })
    })
}
fn reaction(
    name: &str,
    terms: &[(usize, i128)],
) -> Result<ConservedReaction, ReactionBalanceError> {
    let ids = terms
        .iter()
        .map(|&(i, _)| sources()[i].id.clone())
        .collect::<Vec<_>>();
    let matrix = ElementalMatrix::new(
        ["Cu", "C", "N"]
            .into_iter()
            .map(ElementId::new)
            .collect::<Result<Vec<_>, _>>()?,
        ids.clone(),
        (0..3)
            .map(|axis| terms.iter().map(|&(i, _)| ROWS[i].atoms[axis]).collect())
            .collect(),
    )?;
    let id = ReactionId::new(name)?;
    let coefficients = StoichiometricMatrix::new(
        ids.clone(),
        vec![id.clone()],
        terms.iter().map(|&(_, n)| vec![n]).collect(),
    )?;
    let charge = ChargeVector::new(ids, terms.iter().map(|&(i, _)| ROWS[i].charge).collect())?;
    ConservedReaction::new(&matrix, &coefficients, &charge, id)
}
fn formation(n: usize) -> ConservedReaction {
    reaction(
        &format!("copper-i-{n}-cyanide-formation"),
        &[(0, -1), (6, -(n as i128)), (n + 1, 1)],
    )
    .unwrap()
}
fn gibbs(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let s = sources().iter().find(|s| s.id == term.species).unwrap();
            FormationGibbsSource {
                species: &s.id,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_literal_copper_rows_keep_absent_complex_enthalpy_and_entropy_unavailable() {
    for (row, source) in ROWS.iter().zip(sources()) {
        assert_eq!(
            source.pack.claims().claim_count(),
            if row.enthalpy.is_some() { 4 } else { 2 }
        );
        let entropy = if row.charge == 0 {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        for (name, dims, value) in [
            ("molar-mass", MASS, Some(row.mass * 0.001)),
            (
                "standard-molar-formation-gibbs-energy",
                ENERGY,
                Some(row.gibbs * 1000.0),
            ),
            (
                "standard-molar-formation-enthalpy",
                ENERGY,
                row.enthalpy.map(|v| v * 1000.0),
            ),
            (entropy, ENTROPY, row.entropy),
        ] {
            let query = source.pack.claims().query_typed(
                &PropertyKey::new(name, dims),
                &source.at,
                SelectionPolicy::SingleClaimOnly,
            );
            if let Some(value) = value {
                let answer = query.unwrap();
                close(answer.evidence.value.value, value);
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
                assert!(
                    answer.receipt.observation_backed && !answer.receipt.source_hashes.is_empty()
                );
                source
                    .pack
                    .claims()
                    .verify_receipt(&answer.receipt)
                    .unwrap();
            } else {
                assert!(query.is_err(), "{} invents absent {name}", row.pack);
            }
        }
    }
}

#[test]
fn g0_cumulative_complex_formation_retains_free_ligand_demand_and_exact_charge() {
    // Independent 55-digit evaluation from the literal standard Gibbs rows.
    let expected = [
        (2, -136_980.0, 23.997_857_736_424_443),
        (3, -163_380.0, 28.622_937_633_063_407),
        (4, -172_980.0, 30.304_784_868_204_848),
    ];
    for (n, delta_g, log_beta) in expected {
        let reaction = formation(n);
        let standard = reaction.standard_gibbs_change(&gibbs(&reaction)).unwrap();
        close(standard.molar().value(), delta_g);
        close(
            standard.log_equilibrium_constant().value() / std::f64::consts::LN_10,
            log_beta,
        );
        for term in standard.terms() {
            let source = sources()
                .iter()
                .find(|s| s.id == term.term.species)
                .unwrap();
            source
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        let initial = [
            (sources()[0].id.clone(), Amount::new(0.1)),
            (sources()[6].id.clone(), Amount::new(0.15)),
            (sources()[n + 1].id.clone(), Amount::new(0.0)),
        ];
        let extent = reaction.limiting_extent(&initial[..2]).unwrap();
        close(extent.value(), 0.15 / n as f64);
        let amounts = reaction.advance_inventory(&initial, extent).unwrap();
        let amount = |index: usize| {
            amounts
                .iter()
                .find(|(id, _)| id == &sources()[index].id)
                .unwrap()
                .1
                .value()
        };
        close(amount(0) + amount(n + 1), 0.1);
        close(amount(6) + n as f64 * amount(n + 1), 0.15);
        close(
            amount(0) - amount(6) + (1.0 - n as f64) * amount(n + 1),
            -0.05,
        );
        assert!(
            reaction
                .advance_inventory(&initial, Amount::new(extent.value() + 1e-5))
                .is_err()
        );
        // A standard heat cannot be returned from the complex's absent H row.
        let heat = gibbs(&reaction)
            .into_iter()
            .map(|s| FormationEnthalpySource {
                species: s.species,
                pack: s.pack,
                at: s.at,
                reference_axis: s.reference_axis,
            })
            .collect::<Vec<_>>();
        assert!(reaction.standard_enthalpy_change(&heat).is_err());
    }
    // Replacing Cu(I) by Cu(II) without a redox participant violates exact charge.
    assert!(
        reaction(
            "unsupported-cu-ii-complexation",
            &[(1, -1), (6, -2), (3, 1)]
        )
        .is_err()
    );
}

#[test]
fn g3_stepwise_binding_and_solid_dissolution_obey_the_same_gibbs_cycle() {
    let dissolution = reaction("copper-i-cyanide-dissolution", &[(2, -1), (0, 1), (6, 1)]).unwrap();
    let dissolved = dissolution
        .standard_gibbs_change(&gibbs(&dissolution))
        .unwrap();
    close(dissolved.molar().value(), 111_080.0);
    close(
        dissolved.log_equilibrium_constant().value() / std::f64::consts::LN_10,
        -19.460_374_049_949_095,
    );
    let assisted = reaction(
        "cyanide-assisted-cucn-dissolution",
        &[(2, -1), (6, -1), (3, 1)],
    )
    .unwrap();
    let assisted = assisted.standard_gibbs_change(&gibbs(&assisted)).unwrap();
    close(assisted.molar().value(), -25_900.0);
    let dicyano = formation(2);
    close(
        assisted.molar().value(),
        dissolved.molar().value()
            + dicyano
                .standard_gibbs_change(&gibbs(&dicyano))
                .unwrap()
                .molar()
                .value(),
    );
    for (from, to, expected) in [(3, 4, -26_400.0), (4, 5, -9600.0)] {
        let step = reaction("one-additional-cyanide", &[(from, -1), (6, -1), (to, 1)]).unwrap();
        close(
            step.standard_gibbs_change(&gibbs(&step))
                .unwrap()
                .molar()
                .value(),
            expected,
        );
    }
    for n in 2..=4 {
        let reaction = formation(n);
        let standard = reaction.standard_gibbs_change(&gibbs(&reaction)).unwrap();
        let mut activities = vec![
            (sources()[0].id.clone(), Dimensionless::new(1e-24)),
            (sources()[6].id.clone(), Dimensionless::new(0.001)),
            (sources()[n + 1].id.clone(), Dimensionless::new(0.00001)),
        ];
        let before = standard.at_activities(&activities).unwrap().molar.value();
        activities[1].1 = Dimensionless::new(0.002);
        let after = standard.at_activities(&activities).unwrap().molar.value();
        close(after - before, -(n as f64) * RT * std::f64::consts::LN_2);
    }
}

#[test]
fn g0_missing_participants_and_unsupported_original_source_conditions_refuse() {
    let reaction = formation(3);
    let selected = gibbs(&reaction);
    assert!(reaction.standard_gibbs_change(&selected[..2]).is_err());
    let original = &sources()[4];
    for at in [
        original
            .at
            .clone()
            .with("source-phase-aqueous", 0.0)
            .unwrap(),
        original
            .at
            .clone()
            .with("source-ionic-reference-hydrogen-zero", 0.0)
            .unwrap(),
        original
            .at
            .clone()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                310.0,
            )
            .unwrap(),
    ] {
        let mut altered = gibbs(&reaction);
        altered
            .iter_mut()
            .find(|s| s.species == &original.id)
            .unwrap()
            .at = &at;
        assert!(reaction.standard_gibbs_change(&altered).is_err());
    }
}
