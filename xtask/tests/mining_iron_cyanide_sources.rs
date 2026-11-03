//! Source-backed legacy reference chemistry, finite inventories and exact balances.
//! These are arithmetic/closure tests, not experimental validation of the NBS beta.
#![deny(unsafe_code)]

use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

use fs_matdb::{PropertyKey, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        thermodynamics::acid_base::complexation::LigandSpeciation,
    },
};

#[allow(dead_code)]
#[path = "../../crates/fs-thermochem/examples/iron_cyanide_reference.rs"]
mod consumer;

const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);

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

fn compiled() -> &'static Path {
    static DIRECTORY: OnceLock<PathBuf> = OnceLock::new();
    DIRECTORY
        .get_or_init(|| {
            let nonce = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let out = std::env::temp_dir().join(format!(
                "frankensim-iron-cyanide-{}-{nonce}",
                std::process::id()
            ));
            std::fs::create_dir(&out).unwrap();
            for name in consumer::PACKS {
                let target = out.join(format!("{name}.fsmatpk"));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{name}/manifest.tsv")))
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{name}: {}",
                    String::from_utf8_lossy(&run.stderr)
                );
            }
            out
        })
        .as_path()
}

fn sources() -> &'static [consumer::Source] {
    static SOURCES: OnceLock<Vec<consumer::Source>> = OnceLock::new();
    SOURCES.get_or_init(|| consumer::load_sources(compiled()).unwrap())
}

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn relative(actual: f64, expected: f64) {
    assert!(actual.is_finite() && actual > 0.0 && expected > 0.0);
    assert!(
        (actual / expected - 1.0).abs() < 3e-8,
        "{actual} != {expected}"
    );
}

fn member(solution: &LigandSpeciation, id: &str) -> f64 {
    solution
        .metals
        .iter()
        .flat_map(|family| &family.members)
        .find(|entry| entry.species.as_str() == id)
        .unwrap()
        .amount
        .value()
}

fn balances(solution: &LigandSpeciation, amounts: [f64; 6]) {
    let [cyanide, fe_ii, fe_iii, sodium, chloride, _] = amounts;
    for (id, expected) in [(consumer::FE_II, fe_ii), (consumer::FE_III, fe_iii)] {
        let family = solution
            .metals
            .iter()
            .find(|family| family.id == id)
            .unwrap();
        close(family.total.value(), expected);
        close(
            family
                .members
                .iter()
                .map(|entry| entry.amount.value())
                .sum(),
            expected,
        );
        assert_eq!(
            family
                .members
                .iter()
                .map(|entry| entry.ligand_count)
                .collect::<Vec<_>>(),
            [0, 6]
        );
        assert!(
            family
                .members
                .iter()
                .all(|entry| entry.amount.value().is_finite() && entry.amount.value() >= 0.0)
        );
    }
    // Independent C and N accounts: both HCN and CN carry one, both complexes six.
    let aqueous_cn: f64 = solution
        .aqueous
        .families
        .iter()
        .find(|family| family.id == consumer::LIGAND)
        .unwrap()
        .forms
        .iter()
        .map(|entry| entry.amount.value())
        .sum();
    let bound_cn =
        6.0 * (member(solution, "Fe(CN)6[4-](aq)") + member(solution, "Fe(CN)6[3-](aq)"));
    close(aqueous_cn + bound_cn, cyanide);
    close(bound_cn, solution.bound_ligand.value());
    close(aqueous_cn, solution.uncomplexed_ligand.value());
    assert!(solution.bound_ligand.value() <= cyanide + 1e-15);
    let charge = solution.aqueous.proton_amount.value()
        - solution.aqueous.hydroxide_amount.value()
        - solution.free_ligand.value()
        + sodium
        - chloride
        + 2.0 * member(solution, "Fe[2+](aq)")
        + 3.0 * member(solution, "Fe[3+](aq)")
        - 4.0 * member(solution, "Fe(CN)6[4-](aq)")
        - 3.0 * member(solution, "Fe(CN)6[3-](aq)");
    assert!(charge.abs() < 3e-12, "charge residual {charge}");
    assert!(solution.ligand_balance_residual.value().abs() < 3e-12);
    assert!(solution.aqueous.charge_balance_residual.value().abs() < 3e-12);
    assert!(solution.max_log_mass_action_residual.value() < 1e-8);
}

#[test]
fn g0_literal_corrected_nbs_iron_rows_keep_conventions_receipts_and_absent_cp() {
    // Table 41:Fe, prepared 1966; all four exact ao rows, not salts labelled ai.
    for (index, mass, enthalpy, gibbs, entropy) in [
        (5, 55.8470, -89.1, -78.90, -137.7),
        (6, 55.8470, -48.5, -4.7, -315.9),
        (7, 211.9544, 455.6, 695.08, 95.0),
        (8, 211.9544, 561.9, 729.4, 270.3),
    ] {
        let source = &sources()[index];
        assert_eq!(source.pack.claims().claim_count(), 4);
        for (name, dims, expected) in [
            ("molar-mass", MASS, mass * 0.001),
            (
                "standard-molar-formation-enthalpy",
                ENERGY,
                enthalpy * 1000.0,
            ),
            (
                "standard-molar-formation-gibbs-energy",
                ENERGY,
                gibbs * 1000.0,
            ),
            ("conventional-standard-molar-entropy", ENTROPY, entropy),
        ] {
            let answer = source
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &source.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed && !answer.receipt.source_hashes.is_empty());
            source
                .pack
                .claims()
                .verify_receipt(&answer.receipt)
                .unwrap();
        }
        assert!(
            source
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new("molar-isobaric-heat-capacity", ENTROPY),
                    &source.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn g0_exact_six_ligand_inventory_and_legacy_gibbs_are_not_redox_conversion() {
    for (oxidation, parent, complex, delta_g, delta_h, log_beta) in [
        (2, 5, 7, -260_420.0, -358_900.0, 45.623_610_101_618_14),
        (3, 6, 8, -300_300.0, -293_200.0, 52.610_283_824_268_22),
    ] {
        let reaction = consumer::association(oxidation).unwrap();
        let selected = consumer::selected(&reaction, sources());
        let standard = reaction.standard_gibbs_change(&selected).unwrap();
        close(standard.molar().value(), delta_g);
        close(
            standard.log_equilibrium_constant().value() / std::f64::consts::LN_10,
            log_beta,
        );
        let heat = selected
            .iter()
            .map(|source| FormationEnthalpySource {
                species: source.species,
                pack: source.pack,
                at: source.at,
                reference_axis: source.reference_axis,
            })
            .collect::<Vec<_>>();
        close(
            reaction
                .standard_enthalpy_change(&heat)
                .unwrap()
                .molar()
                .value(),
            delta_h,
        );
        let initial = [
            (sources()[parent].id.clone(), Amount::new(0.1)),
            (sources()[2].id.clone(), Amount::new(0.3)),
            (sources()[complex].id.clone(), Amount::new(0.0)),
        ];
        let limit = reaction.limiting_extent(&initial[..2]).unwrap();
        close(limit.value(), 0.05);
        let after = reaction.advance_inventory(&initial, limit).unwrap();
        let amount = |i: usize| {
            after
                .iter()
                .find(|(id, _)| id == &sources()[i].id)
                .unwrap()
                .1
                .value()
        };
        close(amount(parent) + amount(complex), 0.1);
        close(amount(2) + 6.0 * amount(complex), 0.3);
        close(
            oxidation as f64 * amount(parent) - amount(2)
                + (oxidation as f64 - 6.0) * amount(complex),
            oxidation as f64 * 0.1 - 0.3,
        );
        assert!(
            reaction
                .advance_inventory(&initial, Amount::new(0.05001))
                .is_err()
        );
    }
    // Fe(III) + 6CN- -> ferrocyanide loses one unit of charge without an
    // electron donor. A label swap cannot create an implicit redox pathway.
    let ids = [
        consumer::SPECIES[6],
        consumer::SPECIES[2],
        consumer::SPECIES[7],
    ]
    .into_iter()
    .map(SpeciesId::new)
    .collect::<Result<Vec<_>, _>>()
    .unwrap();
    let matrix = ElementalMatrix::new(
        ["Fe", "C", "N"]
            .into_iter()
            .map(ElementId::new)
            .collect::<Result<Vec<_>, _>>()
            .unwrap(),
        ids.clone(),
        vec![vec![1, 0, 1], vec![0, 1, 6], vec![0, 1, 6]],
    )
    .unwrap();
    let id = ReactionId::new("unsupported-iron-redox").unwrap();
    let stoich = StoichiometricMatrix::new(
        ids.clone(),
        vec![id.clone()],
        vec![vec![-1], vec![-6], vec![1]],
    )
    .unwrap();
    let charge = ChargeVector::new(ids, vec![3, -1, -4]).unwrap();
    assert!(ConservedReaction::new(&matrix, &stoich, &charge, id).is_err());
}

#[test]
fn g3_two_iron_pools_share_a_finite_cn_inventory_with_independent_reference() {
    let model = consumer::model_from_sources(sources()).unwrap();
    for (amounts, free_cn, hcn, fe2, fe3, ph) in [
        (
            [0.001, 0.0001, 0.00005, 0.001, 0.00035, 1.0],
            6.618_713_951_384_4e-5,
            3.381_286_048_615_599_5e-5,
            2.829_751_024_716_783e-25,
            1.458_963_770_742_021e-32,
            9.524_332_250_750_164,
        ),
        (
            [0.0006, 0.0001, 0.00005, 0.0006, 0.00035, 1.0],
            2.487_369_727_798_223_8e-8,
            6.479_133_718_510_471e-7,
            5.011_212_599_921_015e-5,
            5.178_978_022_807_479e-12,
            7.816_863_952_595_053,
        ),
    ] {
        // Independent 80-digit Decimal scalar oracle from literal source Gs:
        // Ka=exp(-52700/RT), Kw=exp(-79885/RT), beta2=exp(260420/RT),
        // beta3=exp(300300/RT); w_i=beta_i*c^6; complex_i=F_i*w_i/(1+w_i).
        // These counterions give h=sqrt(Kw/(1+c/Ka)); solve
        // c+c*h/Ka+6*sum(complex_i)=total_CN/n_ref.
        let solution = consumer::solve(&model, amounts).unwrap();
        relative(solution.free_ligand.value(), free_cn);
        let hcn_amount = solution.aqueous.families[0]
            .forms
            .iter()
            .find(|entry| entry.species.as_str() == "HCN(aq)")
            .unwrap()
            .amount
            .value();
        relative(hcn_amount, hcn);
        relative(member(&solution, "Fe[2+](aq)"), fe2);
        relative(member(&solution, "Fe[3+](aq)"), fe3);
        assert!((solution.aqueous.ph.value() - ph).abs() < 3e-8);
        balances(&solution, amounts);
        for family in solution.model().families() {
            assert_eq!(family.steps()[0].ligand_count(), 6);
            for term in family.steps()[0].standard().terms() {
                sources()
                    .iter()
                    .find(|source| source.id == term.term.species)
                    .unwrap()
                    .pack
                    .claims()
                    .verify_receipt(&term.answer.receipt)
                    .unwrap();
            }
        }
    }
}

#[test]
fn g0_zero_stocks_scaling_and_competition_keep_separate_oxidation_state_accounts() {
    let model = consumer::model_from_sources(sources()).unwrap();
    let scarce = [0.0006, 0.0001, 0.00005, 0.0006, 0.00035, 1.0];
    let both = consumer::solve(&model, scarce).unwrap();
    let fe_ii_only = [0.0006, 0.0001, 0.0, 0.0006, 0.0002, 1.0];
    let single = consumer::solve(&model, fe_ii_only).unwrap();
    assert!(member(&both, "Fe(CN)6[4-](aq)") < member(&single, "Fe(CN)6[4-](aq)"));
    assert!(both.free_ligand.value() < single.free_ligand.value());
    let scaled = consumer::solve(&model, scarce.map(|amount| 7.0 * amount)).unwrap();
    relative(scaled.free_ligand.value(), 7.0 * both.free_ligand.value());
    assert!((scaled.aqueous.ph.value() - both.aqueous.ph.value()).abs() < 3e-9);
    for family in &both.metals {
        for entry in &family.members {
            relative(
                member(&scaled, entry.species.as_str()),
                7.0 * entry.amount.value(),
            );
        }
    }
    for amounts in [
        [0.001, 0.0, 0.0, 0.001, 0.0, 1.0],
        [0.0, 0.0001, 0.00005, 0.0, 0.00035, 1.0],
    ] {
        let solution = consumer::solve(&model, amounts).unwrap();
        assert_eq!(solution.bound_ligand.value(), 0.0);
        close(member(&solution, "Fe[2+](aq)"), amounts[1]);
        close(member(&solution, "Fe[3+](aq)"), amounts[2]);
        balances(&solution, amounts);
    }
    balances(&scaled, scarce.map(|amount| 7.0 * amount));
}

#[test]
fn g3_independent_counterion_stocks_change_ph_without_adding_iron_or_cyanide() {
    let model = consumer::model_from_sources(sources()).unwrap();
    let amounts = [0.001, 0.0001, 0.00005, 0.003, 0.00035, 1.0];
    let solution = consumer::solve(&model, amounts).unwrap();
    relative(solution.free_ligand.value(), 9.914_405_685_830_678e-5);
    assert!((solution.aqueous.ph.value() - 11.296_462_320_532_133).abs() < 3e-8);
    balances(&solution, amounts);
}

#[test]
fn g0_missing_sources_and_outside_original_source_conditions_refuse() {
    for index in 5..9 {
        for (axis, value) in [
            ("source-phase-aqueous", 0.0),
            ("source-aqueous-unit-activity-infinite-dilution", 0.0),
            ("source-ionic-reference-hydrogen-zero", 0.0),
            (consumer::REFERENCE, 0.0),
        ] {
            let mut altered = consumer::load_sources(compiled()).unwrap();
            altered[index].at = altered[index].at.clone().with(axis, value).unwrap();
            assert!(
                consumer::model_from_sources(&altered).is_err(),
                "accepted altered {index}/{axis}"
            );
        }
        let mut altered = consumer::load_sources(compiled()).unwrap();
        altered[index].at = altered[index]
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
            .unwrap();
        assert!(consumer::model_from_sources(&altered).is_err());
        let mut missing = consumer::load_sources(compiled()).unwrap();
        missing.remove(index);
        assert!(consumer::model_from_sources(&missing).is_err());
    }
}

#[test]
fn g0_runnable_consumer_requires_reference_and_explicit_finite_stocks() {
    let valid = [
        "--reference",
        "nbs1989-legacy",
        "--materials",
        "/example",
        "--amounts",
        ".001",
        ".0001",
        ".00005",
        ".003",
        ".00035",
        "1",
    ]
    .map(str::to_string)
    .to_vec();
    assert!(consumer::parse_amounts(&valid).is_ok());
    assert!(consumer::parse_amounts(&valid[2..]).is_err());
    for (index, replacement) in [
        (1, "best-available"),
        (5, "NaN"),
        (6, "-1"),
        (7, "inf"),
        (8, "-.1"),
        (9, "NaN"),
        (10, "0"),
    ] {
        let mut invalid = valid.clone();
        invalid[index] = replacement.into();
        assert!(consumer::parse_amounts(&invalid).is_err());
    }
    let model = consumer::model_from_sources(sources()).unwrap();
    for invalid in [[0.0; 6], [1.0, -1.0, 0.0, 1.0, 0.0, 1.0]] {
        assert!(consumer::solve(&model, invalid).is_err());
    }
}
