//! Corrected NBS data for finite cyanide oxidation and retained effluent nitrogen.
//! Source thermodynamics does not assert a treatment rate or removal efficiency.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource, ReactionBalanceError,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const CP: &str = "conventional-standard-molar-isobaric-heat-capacity";
const PEROXIDE: usize = 0;
const HYDROPEROXIDE: usize = 1;
const CYANATE: usize = 2;
const AMMONIA: usize = 3;
const AMMONIUM: usize = 4;
const CYANIDE: usize = 5;
const WATER: usize = 6;
const PROTON: usize = 7;
const BICARBONATE: usize = 8;
const CARBONATE: usize = 9;
const HYDROXIDE: usize = 10;

struct Row {
    id: &'static str,
    pack: &'static str,
    // Explicit C/H/N/O counts and charge; names do not authenticate formulas.
    atoms: [i128; 4],
    charge: i128,
    // Original g/mol, kJ/mol, kJ/mol, J/(mol K), without reconciliation.
    values: [f64; 4],
    heat_capacity: Option<f64>,
}
// NIST corrected CSV DOI10.18434/M32124, SHA256
// cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44.
// Tables2:H(prepared1963),18:N(prepared1964),23:C(prepared1966); every solute is ao.
const ROWS: &[Row] = &[
    Row {
        id: "H2O2(aq)",
        pack: "hydrogen-peroxide-aqueous-nbs1989",
        atoms: [0, 2, 0, 2],
        charge: 0,
        values: [34.0148, -191.17, -134.03, 143.9],
        heat_capacity: None,
    },
    Row {
        id: "HO2[-](aq)",
        pack: "hydroperoxide-ion-aqueous-nbs1989",
        atoms: [0, 1, 0, 2],
        charge: -1,
        values: [33.0068, -160.33, -67.3, 23.8],
        heat_capacity: None,
    },
    Row {
        id: "OCN[-](aq)",
        pack: "cyanate-ion-aqueous-nbs1989",
        atoms: [1, 0, 1, 1],
        charge: -1,
        values: [42.0173, -146.0, -97.4, 106.7],
        heat_capacity: None,
    },
    Row {
        id: "NH3(aq)",
        pack: "ammonia-aqueous-nbs1989",
        atoms: [0, 3, 1, 0],
        charge: 0,
        values: [17.0307, -80.29, -26.50, 111.3],
        heat_capacity: None,
    },
    Row {
        id: "NH4[+](aq)",
        pack: "ammonium-ion-aqueous-nbs1989",
        atoms: [0, 4, 1, 0],
        charge: 1,
        values: [18.0387, -132.51, -79.31, 113.4],
        heat_capacity: Some(79.9),
    },
    Row {
        id: "CN[-](aq)",
        pack: "cyanide-ion-aqueous-nbs1989",
        atoms: [1, 0, 1, 0],
        charge: -1,
        values: [26.0179, 150.6, 172.4, 94.1],
        heat_capacity: None,
    },
    Row {
        id: "H2O(l)",
        pack: "water-liquid-nbs1989-thermochemistry",
        atoms: [0, 2, 0, 1],
        charge: 0,
        values: [18.0154, -285.830, -237.129, 69.91],
        heat_capacity: None,
    },
    Row {
        id: "H+(aq)",
        pack: "hydrogen-ion-aqueous-nbs1989",
        atoms: [0, 1, 0, 0],
        charge: 1,
        values: [1.0080, 0.0, 0.0, 0.0],
        heat_capacity: None,
    },
    Row {
        id: "HCO3[-](aq)",
        pack: "bicarbonate-ion-aqueous-nbs1989",
        atoms: [1, 1, 0, 3],
        charge: -1,
        values: [61.0174, -691.99, -586.77, 91.2],
        heat_capacity: None,
    },
    Row {
        id: "CO3[2-](aq)",
        pack: "carbonate-ion-aqueous-nbs1989",
        atoms: [1, 0, 0, 3],
        charge: -2,
        values: [60.0094, -677.14, -527.81, -56.9],
        heat_capacity: None,
    },
    Row {
        id: "OH[-](aq)",
        pack: "hydroxide-ion-aqueous-nbs1989",
        atoms: [0, 1, 0, 1],
        charge: -1,
        values: [17.0074, -229.994, -157.244, -10.75],
        heat_capacity: None,
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
fn point(index: usize) -> QueryPoint {
    let row = &ROWS[index];
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
    if index == WATER {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if row.charge != 0 {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    for (name, value) in flags {
        at = at
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}
fn sources() -> &'static [Source] {
    static SOURCES: OnceLock<Vec<Source>> = OnceLock::new();
    SOURCES.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-cyanide-effluent-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        ROWS.iter()
            .enumerate()
            .map(|(index, row)| {
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
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), row.pack);
                assert_eq!(pack.to_bytes(), bytes);
                Source {
                    id: SpeciesId::new(row.id).unwrap(),
                    pack,
                    at: point(index),
                }
            })
            .collect()
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
    let elements = ElementalMatrix::new(
        ["C", "H", "N", "O"]
            .into_iter()
            .map(ElementId::new)
            .collect::<Result<Vec<_>, _>>()?,
        ids.clone(),
        (0..4)
            .map(|axis| terms.iter().map(|&(i, _)| ROWS[i].atoms[axis]).collect())
            .collect(),
    )?;
    let id = ReactionId::new(name)?;
    let stoichiometry = StoichiometricMatrix::new(
        ids.clone(),
        vec![id.clone()],
        terms.iter().map(|&(_, n)| vec![n]).collect(),
    )?;
    let charges = ChargeVector::new(ids, terms.iter().map(|&(i, _)| ROWS[i].charge).collect())?;
    ConservedReaction::new(&elements, &stoichiometry, &charges, id)
}
fn gibbs(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let source = sources().iter().find(|s| s.id == term.species).unwrap();
            FormationGibbsSource {
                species: &source.id,
                pack: &source.pack,
                at: &source.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn heat(reaction: &ConservedReaction) -> Vec<FormationEnthalpySource<'static>> {
    gibbs(reaction)
        .into_iter()
        .map(|s| FormationEnthalpySource {
            species: s.species,
            pack: s.pack,
            at: s.at,
            reference_axis: s.reference_axis,
        })
        .collect()
}
fn inventory(values: &[(usize, f64)]) -> Vec<(SpeciesId, Amount)> {
    values
        .iter()
        .map(|&(index, n)| (sources()[index].id.clone(), Amount::new(n)))
        .collect()
}
fn amount(values: &[(SpeciesId, Amount)], index: usize) -> f64 {
    values
        .iter()
        .find(|(id, _)| *id == sources()[index].id)
        .unwrap()
        .1
        .value()
}
fn close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn oxidation() -> ConservedReaction {
    reaction(
        "cyanide-peroxide-oxidation",
        &[(CYANIDE, -1), (PEROXIDE, -1), (CYANATE, 1), (WATER, 1)],
    )
    .unwrap()
}
fn hydrolysis() -> ConservedReaction {
    reaction(
        "cyanate-hydrolysis-to-ammonia-bicarbonate",
        &[(CYANATE, -1), (WATER, -2), (AMMONIA, 1), (BICARBONATE, 1)],
    )
    .unwrap()
}

#[test]
fn literal_nbs_rows_keep_phase_charge_precision_and_unstated_uncertainty() {
    for (row, source) in ROWS[..5].iter().zip(&sources()[..5]) {
        assert_eq!(
            source.pack.claims().claim_count(),
            if row.heat_capacity.is_some() { 5 } else { 4 }
        );
        let entropy = if row.charge == 0 {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        for (name, dims, expected) in [
            ("molar-mass", MASS, Some(row.values[0] * 0.001)),
            (
                "standard-molar-formation-enthalpy",
                ENERGY,
                Some(row.values[1] * 1000.0),
            ),
            (
                "standard-molar-formation-gibbs-energy",
                ENERGY,
                Some(row.values[2] * 1000.0),
            ),
            (entropy, ENTROPY, Some(row.values[3])),
            (CP, ENTROPY, row.heat_capacity),
        ] {
            let answer = source.pack.claims().query_typed(
                &PropertyKey::new(name, dims),
                &source.at,
                SelectionPolicy::SingleClaimOnly,
            );
            if let Some(value) = expected {
                let answer = answer.unwrap();
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
                assert!(
                    answer.is_err(),
                    "{} invents unavailable heat capacity",
                    row.pack
                );
            }
        }
    }
}

#[test]
fn source_backed_oxidation_and_both_hydrolysis_forms_close_the_gibbs_cycle() {
    let ammonium_route = reaction(
        "cyanate-hydrolysis-to-ammonium-carbonate",
        &[(CYANATE, -1), (WATER, -2), (AMMONIUM, 1), (CARBONATE, 1)],
    )
    .unwrap();
    for (reaction, expected_g, expected_h) in [
        (oxidation(), -372_899.0, -391_260.0),
        (hydrolysis(), -41_612.0, -54_620.0),
        (ammonium_route, -35_462.0, -91_990.0),
    ] {
        let standard = reaction.standard_gibbs_change(&gibbs(&reaction)).unwrap();
        close(standard.molar().value(), expected_g);
        close(
            reaction
                .standard_enthalpy_change(&heat(&reaction))
                .unwrap()
                .molar()
                .value(),
            expected_h,
        );
        for term in standard.terms() {
            sources()
                .iter()
                .find(|s| s.id == term.term.species)
                .unwrap()
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
    }
    // The two nitrogen/carbon product choices differ by proton transfer.
    let transfer = reaction(
        "ammonium-carbonate-proton-transfer",
        &[
            (AMMONIUM, -1),
            (CARBONATE, -1),
            (AMMONIA, 1),
            (BICARBONATE, 1),
        ],
    )
    .unwrap();
    close(
        transfer
            .standard_gibbs_change(&gibbs(&transfer))
            .unwrap()
            .molar()
            .value(),
        -6_150.0,
    );
    // Removing water, or silently exchanging NH3 for NH4+, cannot conserve the equation.
    assert!(
        reaction(
            "missing-water",
            &[(CYANIDE, -1), (PEROXIDE, -1), (CYANATE, 1)]
        )
        .is_err()
    );
    assert!(
        reaction(
            "invented-charge",
            &[(CYANATE, -1), (WATER, -2), (AMMONIUM, 1), (BICARBONATE, 1)]
        )
        .is_err()
    );
}

#[test]
fn finite_peroxide_cannot_destroy_more_cyanide_or_erase_effluent_nitrogen() {
    let oxidation = oxidation();
    let initial = inventory(&[
        (CYANIDE, 1.0),
        (PEROXIDE, 0.25),
        (CYANATE, 0.0),
        (WATER, 55.5),
    ]);
    let extent = oxidation
        .limiting_extent(&inventory(&[(CYANIDE, 1.0), (PEROXIDE, 0.25)]))
        .unwrap();
    close(extent.value(), 0.25);
    let oxidized = oxidation.advance_inventory(&initial, extent).unwrap();
    close(amount(&oxidized, CYANIDE), 0.75);
    close(amount(&oxidized, PEROXIDE), 0.0);
    close(amount(&oxidized, CYANATE), 0.25);
    assert!(
        oxidation
            .advance_inventory(&initial, Amount::new(0.250001))
            .is_err()
    );

    let hydrolysis = hydrolysis();
    let hydrolysis_initial = inventory(&[
        (CYANATE, amount(&oxidized, CYANATE)),
        (WATER, amount(&oxidized, WATER)),
        (AMMONIA, 0.0),
        (BICARBONATE, 0.0),
    ]);
    let hydrolyzed = hydrolysis
        .advance_inventory(&hydrolysis_initial, extent)
        .unwrap();
    let final_inventory = inventory(&[
        (CYANIDE, amount(&oxidized, CYANIDE)),
        (PEROXIDE, amount(&oxidized, PEROXIDE)),
        (CYANATE, amount(&hydrolyzed, CYANATE)),
        (WATER, amount(&hydrolyzed, WATER)),
        (AMMONIA, amount(&hydrolyzed, AMMONIA)),
        (BICARBONATE, amount(&hydrolyzed, BICARBONATE)),
    ]);
    close(amount(&final_inventory, AMMONIA), 0.25);
    close(amount(&final_inventory, BICARBONATE), 0.25);
    close(amount(&final_inventory, WATER), 55.25);
    let totals = |values: &[(SpeciesId, Amount)]| {
        let mut totals = [0.0; 6];
        for (id, n) in values {
            let row = ROWS.iter().find(|row| row.id == id.as_str()).unwrap();
            for (total, atoms) in totals[..4].iter_mut().zip(row.atoms) {
                *total += atoms as f64 * n.value();
            }
            totals[4] += row.charge as f64 * n.value();
            totals[5] += row.values[0] * 0.001 * n.value();
        }
        totals
    };
    for (before, after) in totals(&initial).into_iter().zip(totals(&final_inventory)) {
        close(after, before);
    }
    let energy = oxidation
        .standard_enthalpy_change(&heat(&oxidation))
        .unwrap()
        .at_extent(extent)
        .unwrap()
        .value()
        + hydrolysis
            .standard_enthalpy_change(&heat(&hydrolysis))
            .unwrap()
            .at_extent(extent)
            .unwrap()
            .value();
    close(energy, -111_470.0);
    // Water is finite in this reaction-account API, unlike the acid/base solvent basis.
    close(
        hydrolysis
            .limiting_extent(&inventory(&[(CYANATE, 0.25), (WATER, 0.2)]))
            .unwrap()
            .value(),
        0.1,
    );
}

#[test]
fn source_bound_ammonium_and_peroxide_families_share_charge_balanced_ph() {
    let mut families = Vec::new();
    for (name, acid, base, expected_pka) in [
        ("peroxide", PEROXIDE, HYDROPEROXIDE, 11.690_590_208_436_29),
        ("ammonia", AMMONIUM, AMMONIA, 9.251_911_717_481_2),
    ] {
        let reaction = reaction(
            &format!("{name}-dissociation"),
            &[(acid, -1), (base, 1), (PROTON, 1)],
        )
        .unwrap();
        let selected = gibbs(&reaction);
        let acid = AcidDissociation::from_sources(
            reaction,
            sources()[acid].id.clone(),
            sources()[base].id.clone(),
            sources()[PROTON].id.clone(),
            &selected,
        )
        .unwrap();
        close(
            -acid.standard().log_equilibrium_constant().value() / std::f64::consts::LN_10,
            expected_pka,
        );
        families.push(AcidFamily::new(name, vec![acid]).unwrap());
    }
    let water = reaction(
        "water-autoprotolysis",
        &[(WATER, -1), (PROTON, 1), (HYDROXIDE, 1)],
    )
    .unwrap();
    let selected = gibbs(&water);
    let water = WaterAutoprotolysis::from_sources(
        water,
        sources()[WATER].id.clone(),
        sources()[PROTON].id.clone(),
        sources()[HYDROXIDE].id.clone(),
        &selected,
    )
    .unwrap();
    let system = IdealAcidBaseSystem::new(water, families).unwrap();
    let totals = [
        ("ammonia".to_owned(), Amount::new(0.001)),
        ("peroxide".to_owned(), Amount::new(0.002)),
    ];
    let neutral = system
        .speciate(
            &totals,
            &[],
            Amount::new(1.0),
            SpeciationOptions::default(),
            || false,
        )
        .unwrap();
    let alkaline = system
        .speciate(
            &totals,
            &[SpectatorIon {
                species: SpeciesId::new("Na[+](aq)").unwrap(),
                charge: 1,
                amount: Amount::new(0.01),
            }],
            Amount::new(1.0),
            SpeciationOptions::default(),
            || false,
        )
        .unwrap();
    assert!(alkaline.ph.value() > neutral.ph.value());
    let species_amount =
        |answer: &fs_thermochem::stoichiometry::thermodynamics::acid_base::AcidBaseSpeciation,
         index: usize| {
            answer
                .families
                .iter()
                .flat_map(|family| &family.forms)
                .find(|form| form.species == sources()[index].id)
                .unwrap()
                .amount
                .value()
        };
    assert!(species_amount(&alkaline, AMMONIA) > species_amount(&neutral, AMMONIA));
    assert!(species_amount(&alkaline, HYDROPEROXIDE) > species_amount(&neutral, HYDROPEROXIDE));
    for answer in [&neutral, &alkaline] {
        assert!(answer.charge_balance_residual.value().abs() <= 1e-10);
        assert!(answer.max_log_mass_action_residual.value() < 1e-9);
        for family in &answer.families {
            close(
                family.forms.iter().map(|form| form.amount.value()).sum(),
                family.total.value(),
            );
        }
    }
}

#[test]
fn missing_sources_or_unsupported_reference_conditions_cannot_supply_thermodynamics() {
    let reaction = oxidation();
    let selected = gibbs(&reaction);
    assert!(
        reaction
            .standard_gibbs_change(&selected[..selected.len() - 1])
            .is_err()
    );
    let source = &sources()[CYANATE];
    for at in [
        source.at.clone().with("source-phase-aqueous", 0.0).unwrap(),
        source
            .at
            .clone()
            .with("source-ionic-reference-hydrogen-zero", 0.0)
            .unwrap(),
        source
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
        let mut selected = gibbs(&reaction);
        selected
            .iter_mut()
            .find(|s| *s.species == source.id)
            .unwrap()
            .at = &at;
        assert!(reaction.standard_gibbs_change(&selected).is_err());
    }
}
