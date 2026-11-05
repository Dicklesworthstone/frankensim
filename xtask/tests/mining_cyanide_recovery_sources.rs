//! G0/G3: actual material compiler, NBS identities and finite Cu/S/CN recovery.
//! The working standalone consumer is used directly, with no alternate parser.
#![deny(unsafe_code)]
use fs_matdb::{PropertyKey, SelectionPolicy, UncertaintyModel};
use fs_qty::{Dims, MolarMass};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        FormationEnthalpySource,
        minerals::MineralReaction,
        network::ReactionNetwork,
        thermodynamics::acid_base::{
            complexation::precipitation::CopperSulfideSpeciation,
            precipitation::PrecipitationStatus,
        },
    },
};
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

#[allow(dead_code)]
#[path = "../../crates/fs-thermochem/examples/copper_cyanide_recovery.rs"]
mod consumer;
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);

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
                "frankensim-cyanide-recovery-{}-{nonce}",
                std::process::id()
            ));
            std::fs::create_dir(&out).unwrap();
            for (_, name, _) in consumer::ROWS {
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
            eprintln!("compiled Cu2S/cyanide recovery packs: {}", out.display());
            out
        })
        .as_path()
}
fn sources() -> &'static [consumer::Source] {
    static SOURCES: OnceLock<Vec<consumer::Source>> = OnceLock::new();
    SOURCES.get_or_init(|| consumer::load_sources(compiled()).unwrap())
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
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-13 * b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn thirteen_literal_nbs_facts_preserve_molecular_ionic_and_alpha_phase_identity() {
    for (i, row) in [
        [0.0340800, -39_700.0, -27_830.0, 121.0],
        [0.0330720, -17_600.0, 12_080.0, 62.8],
        [0.1591440, -79_500.0, -86_200.0, 120.9],
    ]
    .into_iter()
    .enumerate()
    {
        let s = &sources()[i];
        assert_eq!(s.pack.claims().claim_count(), if i == 2 { 5 } else { 4 });
        for ((key, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", ENERGY),
            ("standard-molar-formation-gibbs-energy", ENERGY),
            (
                if i == 1 {
                    "conventional-standard-molar-entropy"
                } else {
                    "standard-molar-entropy"
                },
                ENTROPY,
            ),
        ]
        .into_iter()
        .zip(row)
        {
            close(value(s, key, dims), expected);
        }
    }
    close(
        value(
            &sources()[2],
            "standard-molar-isobaric-heat-capacity",
            ENTROPY,
        ),
        76.32,
    );
    for s in &sources()[..2] {
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new("standard-molar-isobaric-heat-capacity", ENTROPY),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new(
                        "conventional-standard-molar-isobaric-heat-capacity",
                        ENTROPY
                    ),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
    assert!(
        sources()[1]
            .pack
            .claims()
            .query_typed(
                &PropertyKey::new("standard-molar-entropy", ENTROPY),
                &sources()[1].at,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn source_queries_refuse_extrapolation_wrong_phase_and_unprovided_rate() {
    let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
    for (s, qualifier) in sources()[..3].iter().zip([
        "source-phase-aqueous",
        "source-individual-ion",
        "source-phase-copper-i-sulfide-nbs-alpha",
    ]) {
        for p in [
            s.at.clone().with("temperature", 300.0).unwrap(),
            s.at.clone().with("pressure", 101_325.0).unwrap(),
            s.at.clone().with(qualifier, 0.0).unwrap(),
            s.at.clone().with(consumer::REFERENCE, 0.0).unwrap(),
            s.at.clone().with("source-specimen-qualified", 1.0).unwrap(),
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &p, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new("sart-recovery-rate", Dims::NONE),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn source_derived_phase_and_sart_gibbs_use_the_original_compatible_species() {
    let m = consumer::model_from_sources(sources()).unwrap();
    close(m.standard().molar().value(), -198_240.0);
    close(m.standard().temperature().value(), 298.15);
    close(m.standard().reference_pressure().value(), 100_000.0);
    let acid = MineralReaction::HydrogenSulfideDissociation
        .build()
        .unwrap();
    let g = acid
        .standard_gibbs_change(&consumer::selected(&acid, sources()))
        .unwrap();
    close(g.molar().value(), 39_910.0);
    let net = MineralReaction::CopperCyanideHydrosulfideRecovery
        .build()
        .unwrap();
    close(
        net.standard_gibbs_change(&consumer::selected(&net, sources()))
            .unwrap()
            .molar()
            .value(),
        -187_680.0,
    );
    let mut reactions = m.reaction_network().reactions().to_vec();
    reactions.push(net.clone());
    ReactionNetwork::new(reactions).unwrap();
    let heat = |r: &fs_thermochem::stoichiometry::ConservedReaction| {
        r.terms()
            .iter()
            .map(|t| {
                let s = sources().iter().find(|s| s.id == t.species).unwrap();
                FormationEnthalpySource {
                    species: &s.id,
                    pack: &s.pack,
                    at: &s.at,
                    reference_axis: consumer::REFERENCE,
                }
            })
            .collect::<Vec<_>>()
    };
    close(
        m.reaction()
            .standard_enthalpy_change(&heat(m.reaction()))
            .unwrap()
            .molar()
            .value(),
        -205_240.0,
    );
    assert!(
        net.standard_enthalpy_change(&heat(&net)).is_err(),
        "Cu complex enthalpies are absent; no SART heat can be invented"
    );
}

fn independent_balances(a: &CopperSulfideSpeciation, cn: f64, cu: f64, s: f64, na: f64, cl: f64) {
    let family = &a.ligand.metals[0];
    let actual_cu: f64 = family.members.iter().map(|m| m.amount.value()).sum();
    let bound: f64 = family
        .members
        .iter()
        .map(|m| m.ligand_count as f64 * m.amount.value())
        .sum();
    let acids = &a.ligand.aqueous.families;
    let uncomplexed: f64 = acids
        .iter()
        .find(|f| f.id == consumer::CYANIDE)
        .unwrap()
        .forms
        .iter()
        .map(|f| f.amount.value())
        .sum();
    let sulfur = acids.iter().find(|f| f.id == consumer::SULFIDE).unwrap();
    let actual_s: f64 = sulfur.forms.iter().map(|f| f.amount.value()).sum();
    close(actual_cu + 2.0 * a.solid_amount.value(), cu);
    close(actual_s + a.solid_amount.value(), s);
    close(uncomplexed + bound, cn);
    let q = na - cl + a.ligand.aqueous.proton_amount.value()
        - a.ligand.aqueous.hydroxide_amount.value()
        - a.ligand.free_ligand.value()
        - sulfur.forms[1].amount.value()
        + family
            .members
            .iter()
            .map(|m| f64::from(m.charge) * m.amount.value())
            .sum::<f64>();
    assert!(q.abs() < 3e-13, "independent charge {q:e}");
    if a.status == PrecipitationStatus::SolidPresent {
        let h = a.ligand.aqueous.proton_activity.value();
        let free_cu = family.members[0].activity.value();
        let hs = sulfur.forms[1].activity.value();
        let r = 2.0 * free_cu.ln() + hs.ln() - h.ln() + 198_240.0 / (8.31446261815324 * 298.15);
        assert!(r.abs() < 1e-7, "independent saturation {r:e}");
    }
}

#[test]
fn actual_sources_release_bound_cyanide_after_explicit_copper_solid_separation() {
    let model = consumer::model_from_sources(sources()).unwrap();
    let a = consumer::run_stages(
        &model,
        [0.010, 0.002, 0.009, 0.002, 0.0011, 0.0082, 0.008, 10.0],
    )
    .unwrap();
    independent_balances(&a.feed, 0.010, 0.002, 0.0, 0.009, 0.002);
    independent_balances(&a.acidified, 0.010, 0.002, 0.0011, 0.0101, 0.0102);
    independent_balances(
        &a.recycled,
        0.010,
        a.acidified.dissolved_copper.value(),
        a.acidified.dissolved_sulfide.value(),
        0.0181,
        0.0102,
    );
    assert_eq!(a.acidified.status, PrecipitationStatus::SolidPresent);
    assert!(
        a.acidified.dissolved_copper.value() > 0.0
            && a.acidified.dissolved_copper.value() < 0.00001
    );
    assert!(a.acidified.ligand.bound_ligand.value() < a.feed.ligand.bound_ligand.value() / 100.0);
    assert!(a.recycled.ligand.free_ligand.value() > 2.0 * a.feed.ligand.free_ligand.value());
    assert!((4.9..5.1).contains(&a.acidified.ligand.aqueous.ph.value()));
    assert!((9.0..10.0).contains(&a.recycled.ligand.aqueous.ph.value()));
    close(
        a.recycled.dissolved_copper.value()
            + 2.0 * (a.separated_solid + a.recycled.solid_amount.value()),
        0.002,
    );
    close(
        a.recycled.dissolved_sulfide.value() + a.separated_solid + a.recycled.solid_amount.value(),
        0.0011,
    );
    let repeated = consumer::run_stages(
        &model,
        [0.010, 0.002, 0.009, 0.002, 0.0011, 0.0082, 0.008, 10.0],
    )
    .unwrap();
    assert_eq!(
        a.acidified.solid_amount.value().to_bits(),
        repeated.acidified.solid_amount.value().to_bits()
    );
    assert_eq!(
        a.recycled.ligand.free_ligand.value().to_bits(),
        repeated.recycled.ligand.free_ligand.value().to_bits()
    );
}

#[test]
fn real_sources_do_not_ignore_sulfide_shortage_or_high_ph_recycle() {
    let model = consumer::model_from_sources(sources()).unwrap();
    let limited = consumer::solve(&model, [0.010, 0.002, 0.0001, 0.0091, 0.0092, 10.0]).unwrap();
    independent_balances(&limited, 0.010, 0.002, 0.0001, 0.0091, 0.0092);
    assert!(limited.solid_amount.value() <= 0.0001);
    assert!(limited.dissolved_copper.value() >= 0.0018);
    assert!(
        consumer::run_stages(
            &model,
            [0.010, 0.002, 0.009, 0.002, 0.0011, 0.0082, 0.1, 10.0]
        )
        .is_err()
    );
    assert_eq!(model.solid(), &SpeciesId::new(consumer::SOLID).unwrap());
}
