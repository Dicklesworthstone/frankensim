//! G0/G1/G3: actual NBS sulfate/bisulfate packs, coupled water/charge closure.
//! The analytical pool is HSO4- + SO4(2-); metal pairs and H2SO4 are excluded.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        FormationEnthalpySource,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidBaseSpeciation, AcidDissociation, AcidFamily, IdealAcidBaseSystem,
                SpeciationOptions, SpectatorIon, WaterAutoprotolysis,
            },
        },
    },
};
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-sulfate";
const PACKS: [&str; 5] = [
    "hydrogen-sulfate-ion-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "sulfate-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
];
const SPECIES: [&str; 5] = [
    "HSO4[-](aq)",
    "H+(aq)",
    "SO4[2-](aq)",
    "OH[-](aq)",
    "H2O(l)",
];

struct Materials {
    packs: Vec<NormalizedPack>,
    points: Vec<QueryPoint>,
    ids: Vec<SpeciesId>,
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
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if index == 4 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-individual-ion", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn compile(name: &str) -> NormalizedPack {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out =
        std::env::temp_dir().join(format!("frankensim-sulfate-{}-{nonce}", std::process::id()));
    std::fs::create_dir(&out).unwrap();
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
    let bytes = std::fs::read(target).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.pack_id(), name);
    pack
}

fn materials() -> &'static Materials {
    static MATERIALS: OnceLock<Materials> = OnceLock::new();
    MATERIALS.get_or_init(|| Materials {
        packs: PACKS.iter().map(|name| compile(name)).collect(),
        points: (0..5).map(point).collect(),
        ids: SPECIES
            .iter()
            .map(|s| SpeciesId::new(*s).unwrap())
            .collect(),
    })
}

fn sources(indices: &[usize]) -> Vec<FormationGibbsSource<'static>> {
    let m = materials();
    indices
        .iter()
        .map(|&i| FormationGibbsSource {
            species: &m.ids[i],
            pack: &m.packs[i],
            at: &m.points[i],
            reference_axis: REFERENCE,
        })
        .collect()
}

fn system() -> IdealAcidBaseSystem {
    let ids = &materials().ids;
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build().unwrap(),
        ids[4].clone(),
        ids[1].clone(),
        ids[3].clone(),
        &sources(&[4, 1, 3]),
    )
    .unwrap();
    let acid = AcidDissociation::from_sources(
        MineralReaction::HydrogenSulfateDissociation
            .build()
            .unwrap(),
        ids[0].clone(),
        ids[2].clone(),
        ids[1].clone(),
        &sources(&[0, 1, 2]),
    )
    .unwrap();
    IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, vec![acid]).unwrap()]).unwrap()
}

fn solve(
    system: &IdealAcidBaseSystem,
    total: f64,
    sodium: f64,
    chloride: f64,
    basis: f64,
) -> AcidBaseSpeciation {
    let spectators = [("Na[+](aq)", 1, sodium), ("Cl[-](aq)", -1, chloride)]
        .into_iter()
        .map(|(name, charge, amount)| SpectatorIon {
            species: SpeciesId::new(name).unwrap(),
            charge,
            amount: Amount::new(amount),
        })
        .collect::<Vec<_>>();
    system
        .speciate(
            &[(FAMILY.to_owned(), Amount::new(total))],
            &spectators,
            Amount::new(basis),
            SpeciationOptions {
                ph_tolerance: Dimensionless::new(1e-11),
                absolute_charge_tolerance: Amount::new(1e-14),
                relative_charge_tolerance: Dimensionless::new(1e-12),
                ..SpeciationOptions::default()
            },
            || false,
        )
        .unwrap()
}

fn member(answer: &AcidBaseSpeciation, species: &str) -> f64 {
    answer.families[0]
        .forms
        .iter()
        .find(|s| s.species.as_str() == species)
        .unwrap()
        .amount
        .value()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}
fn constants() -> (f64, f64) {
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    ((-11_380.0 / rt).exp(), (-79_885.0 / rt).exp())
}
fn check_closure(answer: &AcidBaseSpeciation, total: f64) {
    let bisulfate = member(answer, SPECIES[0]);
    let sulfate = member(answer, SPECIES[2]);
    close(bisulfate + sulfate, total, 2e-15 * total.max(1.0));
    let charge =
        answer.proton_amount.value() - answer.hydroxide_amount.value() - bisulfate - 2.0 * sulfate
            + answer
                .spectators
                .iter()
                .map(|s| s.amount.value() * f64::from(s.charge))
                .sum::<f64>();
    close(charge, 0.0, 1e-12);
    close(answer.charge_balance_residual.value(), charge, 1e-16);
    let (ka, kw) = constants();
    close(
        answer.proton_activity.value() * answer.hydroxide_activity.value(),
        kw,
        2e-10 * kw,
    );
    if total > 0.0 {
        close(
            answer.proton_activity.value() * sulfate / bisulfate,
            ka,
            2e-10 * ka,
        );
    }
    assert!(answer.max_log_mass_action_residual.value() < 1e-10);
}

#[test]
fn g0_g1_original_nbs_values_support_distinct_direct_gibbs_and_caloric_sums() {
    let m = materials();
    for (i, expected) in [
        (0, [0.0970696, -887_340.0, -755_910.0, 131.8]),
        (2, [0.0960616, -909_270.0, -744_530.0, 20.1]),
    ] {
        assert_eq!(m.packs[i].claims().claim_count(), 4);
        for ((name, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            (
                "standard-molar-formation-enthalpy",
                Dims([2, 1, -2, 0, 0, -1]),
            ),
            (
                "standard-molar-formation-gibbs-energy",
                Dims([2, 1, -2, 0, 0, -1]),
            ),
            (
                "conventional-standard-molar-entropy",
                Dims([2, 1, -2, -1, 0, -1]),
            ),
        ]
        .into_iter()
        .zip(expected)
        {
            let answer = m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected, 1e-9);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            m.packs[i].claims().verify_receipt(&answer.receipt).unwrap();
            assert!(answer.receipt.observation_backed);
        }
    }
    let system = system();
    let acid = &system.families()[0].steps()[0];
    close(acid.standard().molar().value(), 11_380.0, 1e-8);
    close(
        -acid.standard().log_equilibrium_constant().value() / core::f64::consts::LN_10,
        1.993_689_743_323_917,
        1e-12,
    );
    let sources = sources(&[0, 1, 2]);
    let heat = sources
        .iter()
        .map(|s| FormationEnthalpySource {
            species: s.species,
            pack: s.pack,
            at: s.at,
            reference_axis: s.reference_axis,
        })
        .collect::<Vec<_>>();
    let caloric = acid.reaction().standard_equilibrium(&heat).unwrap();
    close(caloric.enthalpy().molar().value(), -21_930.0, 1e-8);
    close(caloric.molar_entropy().value(), -111.7, 1e-10);
    close(caloric.molar_gibbs().value(), 11_373.355, 1e-8);
    assert_ne!(
        caloric.molar_gibbs().value(),
        acid.standard().molar().value()
    );
    for standard in [acid.standard(), system.water().standard()] {
        for term in standard.terms() {
            let pack = m
                .packs
                .iter()
                .find(|p| p.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        }
    }
}

#[test]
fn g1_sulfuric_acid_pool_and_bisulfate_salt_obey_independent_quadratic_bounds() {
    let system = system();
    let (ka, kw) = constants();
    // Ideal analytical pool at T=.01 on a fixed molal basis. Without water,
    // h²+(Ka+C-T)h+Ka(C-2T)=0, C being fixed positive spectator equivalents.
    // The full charge function has derivative >=1, so omitted water shifts
    // h upward by at most Kw/h0. It is bounded here, not silently discarded.
    for (sodium, expected_h0) in [
        (0.0, 0.014_172_259_756_707_81),
        (0.1, 0.006_205_151_187_189_782),
    ] {
        let answer = solve(&system, 0.1, sodium, 0.0, 10.0);
        let t = 0.01;
        let c = sodium / 10.0;
        let b = ka + c - t;
        let h0 = 2.0 * ka * (2.0 * t - c) / (b + (b * b + 4.0 * ka * (2.0 * t - c)).sqrt());
        close(h0, expected_h0, 1e-16);
        let h = answer.proton_activity.value();
        assert!(h >= h0 - 1e-12 * h0 && h <= h0 + kw / h0 + 1e-12 * h0);
        check_closure(&answer, 0.1);
    }
}

#[test]
fn g1_g3_neutralization_changes_ph_and_free_sulfate_with_conserved_total_and_charge() {
    let system = system();
    let states: Vec<_> = [0.0, 0.1, 0.2, 0.22]
        .into_iter()
        .map(|sodium| solve(&system, 0.1, sodium, 0.0, 10.0))
        .collect();
    for state in &states {
        check_closure(state, 0.1);
    }
    for pair in states.windows(2) {
        assert!(pair[0].ph.value() < pair[1].ph.value());
        assert!(member(&pair[0], SPECIES[2]) < member(&pair[1], SPECIES[2]));
    }
    let acidified = solve(&system, 0.1, 0.2, 0.05, 10.0);
    check_closure(&acidified, 0.1);
    assert!(acidified.ph.value() < states[2].ph.value());
    assert!(member(&acidified, SPECIES[0]) > member(&states[2], SPECIES[0]));
    let scaled = solve(&system, 0.3, 0.3, 0.0, 30.0);
    close(scaled.ph.value(), states[1].ph.value(), 1e-11);
    close(
        member(&scaled, SPECIES[2]),
        3.0 * member(&states[1], SPECIES[2]),
        1e-12,
    );
    let pure_water = solve(&system, 0.0, 0.0, 0.0, 10.0);
    check_closure(&pure_water, 0.0);
    close(pure_water.ph.value(), 6.997_623_248_920_524, 1e-10);
}

#[test]
fn g3_reference_mixing_missing_sources_and_out_of_reference_temperature_refuse() {
    let reaction = MineralReaction::HydrogenSulfateDissociation
        .build()
        .unwrap();
    let mut sources = sources(&[0, 1, 2]);
    assert!(reaction.standard_gibbs_change(&sources[..2]).is_err());
    let wrong_temperature = materials().points[0]
        .clone()
        .with("temperature", 299.15)
        .unwrap();
    sources[0].at = &wrong_temperature;
    assert!(reaction.standard_gibbs_change(&sources).is_err());
    sources[0].at = &materials().points[0];
    let original_usgs = compile("sulfate-ion-aqueous-usgs-b2131");
    let both_flags = materials().points[2]
        .clone()
        .with_quantity(
            "source-elemental-reference-usgs-b2131",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    sources[2].pack = &original_usgs;
    sources[2].at = &both_flags;
    // A caller NBS query flag cannot relabel a source-owned USGS evaluation.
    assert!(reaction.standard_gibbs_change(&sources).is_err());
}
