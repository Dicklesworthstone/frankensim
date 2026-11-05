//! Original material compiler and standalone consumer, with independent Ca/C/S
//! balances and a two-phase known answer manufactured from literal NBS energies.
#![deny(unsafe_code)]

use fs_matdb::{PropertyKey, SelectionPolicy, UncertaintyModel};
use fs_qty::{Amount, Dims, MolarMass};
use fs_thermochem::stoichiometry::{
    FormationEnthalpySource,
    minerals::MineralReaction,
    thermodynamics::acid_base::precipitation::{
        PrecipitationStatus, salt::assemblage::SaltAssemblageSpeciation,
    },
};
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

#[allow(dead_code)]
#[path = "../../crates/fs-thermochem/examples/limestone_sulfate_equilibrium.rs"]
mod consumer;

const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const RT: f64 = 8.314_462_618_153_24 * 298.15;

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .into()
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
                "frankensim-limestone-sulfate-{}-{nonce}",
                std::process::id()
            ));
            std::fs::create_dir(&out).unwrap();
            for (_, name, _) in consumer::ROWS {
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{name}/manifest.tsv")))
                    .arg("--out")
                    .arg(out.join(format!("{name}.fsmatpk")))
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
            eprintln!("compiled calcite/gypsum source packs: {}", out.display());
            out
        })
        .as_path()
}
fn sources() -> &'static [consumer::Source] {
    static SOURCES: OnceLock<Vec<consumer::Source>> = OnceLock::new();
    SOURCES.get_or_init(|| consumer::load_sources(compiled()).unwrap())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; tolerance {tolerance:e}"
    );
}
fn value(s: &consumer::Source, key: &str, dims: Dims) -> f64 {
    let a = s
        .pack
        .claims()
        .query_typed(
            &PropertyKey::new(key, dims),
            &s.at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
    assert!(a.receipt.observation_backed && !a.receipt.source_hashes.is_empty());
    s.pack.claims().verify_receipt(&a.receipt).unwrap();
    a.evidence.value.value
}

#[test]
fn five_literal_calcite_facts_retain_the_original_nbs_phase_and_receipts() {
    let s = &sources()[0];
    assert_eq!(s.pack.claims().claim_count(), 5);
    for (key, dims, expected) in [
        ("molar-mass", MolarMass::DIMS, 0.100_089_4),
        ("standard-molar-formation-enthalpy", ENERGY, -1_206_920.0),
        (
            "standard-molar-formation-gibbs-energy",
            ENERGY,
            -1_128_790.0,
        ),
        ("standard-molar-entropy", ENTROPY, 92.9),
        ("standard-molar-isobaric-heat-capacity", ENTROPY, 81.88),
    ] {
        close(value(s, key, dims), expected, 2e-15 * expected.abs());
    }
}

#[test]
fn calcite_source_refuses_temperature_reference_phase_and_ore_rate_substitution() {
    let s = &sources()[0];
    let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
    for at in [
        s.at.clone().with("temperature", 300.0).unwrap(),
        s.at.clone().with("pressure", 101_325.0).unwrap(),
        s.at.clone().with(consumer::REFERENCE, 0.0).unwrap(),
        s.at.clone().with("source-phase-calcite", 0.0).unwrap(),
        s.at.clone().with("source-pure-species", 0.0).unwrap(),
        s.at.clone().with("source-specimen-qualified", 1.0).unwrap(),
    ] {
        assert!(
            s.pack
                .claims()
                .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
    for key in [
        "ore-neutralization-rate",
        "limestone-purity",
        "gypsum-passivation-rate",
    ] {
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new(key, Dims::NONE),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn original_conserved_equations_use_direct_nbs_gibbs_without_reconciling_columns() {
    let model = consumer::model_from_sources(sources()).unwrap();
    assert_eq!(model.reaction_network().reactions().len(), 6);
    for (recipe, expected) in [
        (MineralReaction::CalcitePrecipitation, -47_400.0),
        (MineralReaction::GypsumPrecipitation, -24_912.0),
        (
            MineralReaction::CalciteCarbonicAcidNeutralization,
            -47_870.0,
        ),
        (MineralReaction::CalciteBicarbonateNeutralization, -11_560.0),
        (MineralReaction::CarbonicComponentDissociation, 36_310.0),
        (MineralReaction::BicarbonateDissociation, 58_960.0),
        (MineralReaction::HydrogenSulfateDissociation, 11_380.0),
        (MineralReaction::WaterDissociation, 79_885.0),
    ] {
        let r = recipe.build().unwrap();
        let g = r
            .standard_gibbs_change(&consumer::selected(&r, sources()))
            .unwrap();
        close(g.molar().value(), expected, 1e-8);
        close(g.temperature().value(), 298.15, 0.0);
        close(g.reference_pressure().value(), 100_000.0, 0.0);
    }
    let r = MineralReaction::CalcitePrecipitation.build().unwrap();
    let heats: Vec<_> = r
        .terms()
        .iter()
        .map(|term| {
            let s = sources().iter().find(|s| s.id == term.species).unwrap();
            FormationEnthalpySource {
                species: &s.id,
                pack: &s.pack,
                at: &s.at,
                reference_axis: consumer::REFERENCE,
            }
        })
        .collect();
    close(
        r.standard_enthalpy_change(&heats).unwrap().molar().value(),
        13_050.0,
        1e-8,
    );
    // Independent entropy column gives 202.9 J/(mol K), hence H-T*S=-47444.635.
    // The directly tabulated formation-G route remains -47400, not a fitted mean.
    let entropy = value(&sources()[0], "standard-molar-entropy", ENTROPY)
        - value(
            &sources()[2],
            "conventional-standard-molar-entropy",
            ENTROPY,
        )
        - value(
            &sources()[5],
            "conventional-standard-molar-entropy",
            ENTROPY,
        );
    close(entropy, 202.9, 1e-12);
    close(13_050.0 - 298.15 * entropy, -47_444.635, 1e-8);
}

fn solid(a: &SaltAssemblageSpeciation, id: &str) -> f64 {
    a.phases
        .iter()
        .find(|p| p.solid.as_str() == id)
        .unwrap()
        .amount
        .value()
}
fn independent_balances(a: &SaltAssemblageSpeciation, charge: consumer::Charge) {
    let carbon = &a
        .aqueous
        .families
        .iter()
        .find(|f| f.id == consumer::CARBON)
        .unwrap()
        .forms;
    let sulfate = &a
        .aqueous
        .families
        .iter()
        .find(|f| f.id == consumer::SULFATE)
        .unwrap()
        .forms;
    let calcite = solid(a, consumer::CALCITE);
    let gypsum = solid(a, consumer::GYPSUM);
    let ca = a.dissolved_cation.value();
    close(ca + calcite + gypsum, charge.calcium.value(), 1e-13);
    close(
        carbon.iter().map(|f| f.amount.value()).sum::<f64>() + calcite,
        charge.carbon.value(),
        1e-13,
    );
    close(
        sulfate.iter().map(|f| f.amount.value()).sum::<f64>() + gypsum,
        charge.sulfur.value(),
        1e-13,
    );
    let q = 2.0 * ca + charge.sodium.value() - charge.chloride.value()
        + a.aqueous.proton_amount.value()
        - a.aqueous.hydroxide_amount.value()
        - carbon[1].amount.value()
        - 2.0 * carbon[2].amount.value()
        - sulfate[0].amount.value()
        - 2.0 * sulfate[1].amount.value();
    close(q, 0.0, 3e-13);
    for (phase, anion, energy) in [
        (consumer::CALCITE, carbon[2].activity.value(), 47_400.0),
        (consumer::GYPSUM, sulfate[1].activity.value(), 24_912.0),
    ] {
        let p = a.phases.iter().find(|p| p.solid.as_str() == phase).unwrap();
        if p.status != PrecipitationStatus::EmptyComponent {
            let residual = (ca / charge.basis.value()).ln() + anion.ln() + energy / RT;
            if p.status == PrecipitationStatus::SolidPresent {
                assert!(p.amount.value() > 0.0);
                close(residual, 0.0, 1e-9);
            } else {
                assert_eq!(p.amount.value(), 0.0);
                assert!(residual <= 0.0);
            }
        }
    }
    assert!(a.aqueous_solves <= a.options.max_aqueous_solves);
}

#[test]
fn independent_source_polynomials_recover_a_manufactured_two_phase_state() {
    // Choose pH=8, aCa=.005, calcite=.02 and gypsum=.03 mol independently.
    // Literal source energies determine terminal-ion activities and the complete
    // protonation polynomials; electroneutrality determines the sodium input.
    // No numerical result from the owner is used to manufacture these stocks.
    let h = 1e-8;
    let ca = 0.005;
    let co3 = (-47_400.0 / RT).exp() / ca;
    let hco3 = co3 * h / (-58_960.0 / RT).exp();
    let h2co3 = hco3 * h / (-36_310.0 / RT).exp();
    let so4 = (-24_912.0 / RT).exp() / ca;
    let hso4 = so4 * h / (-11_380.0 / RT).exp();
    let oh = (-79_885.0 / RT).exp() / h;
    let sodium = hco3 + 2.0 * co3 + hso4 + 2.0 * so4 + oh - h - 2.0 * ca;
    assert!(sodium > 0.0);
    let charge = consumer::Charge {
        calcium: Amount::new(ca + 0.02 + 0.03),
        carbon: Amount::new(h2co3 + hco3 + co3 + 0.02),
        sulfur: Amount::new(hso4 + so4 + 0.03),
        sodium: Amount::new(sodium),
        chloride: Amount::new(0.0),
        basis: Amount::new(1.0),
    };
    let model = consumer::model_from_sources(sources()).unwrap();
    let answer = consumer::solve(&model, charge).unwrap();
    close(answer.aqueous.ph.value(), 8.0, 2e-8);
    close(answer.dissolved_cation.value(), ca, 2e-11);
    close(solid(&answer, consumer::CALCITE), 0.02, 2e-11);
    close(solid(&answer, consumer::GYPSUM), 0.03, 2e-11);
    independent_balances(&answer, charge);
}

#[test]
fn finite_acid_addition_dissolves_calcite_and_forms_gypsum_without_reusing_calcium() {
    let model = consumer::model_from_sources(sources()).unwrap();
    let charges = consumer::demonstration_charges();
    let before = consumer::solve(&model, charges[0]).unwrap();
    let after = consumer::solve(&model, charges[1]).unwrap();
    for (a, charge) in [(&before, charges[0]), (&after, charges[1])] {
        independent_balances(a, charge);
        assert!(
            a.phases
                .iter()
                .all(|p| p.status == PrecipitationStatus::SolidPresent)
        );
    }
    assert!(after.aqueous.ph.value() < before.aqueous.ph.value());
    assert!(solid(&after, consumer::CALCITE) < solid(&before, consumer::CALCITE));
    assert!(solid(&after, consumer::GYPSUM) > solid(&before, consumer::GYPSUM));
    assert!(after.dissolved_cation.value() > before.dissolved_cation.value());
    // The released calcite calcium is independently partitioned between the
    // additional gypsum and the increase in dissolved calcium, with no makeup.
    close(
        solid(&before, consumer::CALCITE) - solid(&after, consumer::CALCITE),
        solid(&after, consumer::GYPSUM) - solid(&before, consumer::GYPSUM)
            + after.dissolved_cation.value()
            - before.dissolved_cation.value(),
        1e-13,
    );
}
