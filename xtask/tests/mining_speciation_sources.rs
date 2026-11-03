//! Real NBS packs through the coupled ideal aqueous acid/base owner.
//! Five original source profiles, explicit spectator ions, fixed molal basis.
//! Quadratic comparisons below bound the water correction rather than omit it.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
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

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-cyanide";
const PACKS: [&str; 5] = [
    "hydrogen-cyanide-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "cyanide-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
];
const SPECIES: [&str; 5] = ["HCN(aq)", "H+(aq)", "CN[-](aq)", "OH[-](aq)", "H2O(l)"];

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
    if index == 4 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if index != 0 {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

fn materials() -> &'static Materials {
    static MATERIALS: OnceLock<Materials> = OnceLock::new();
    MATERIALS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-speciation-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        let packs = PACKS
            .iter()
            .map(|name| {
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
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), *name);
                pack
            })
            .collect();
        Materials {
            packs,
            points: (0..5).map(point).collect(),
            ids: SPECIES
                .iter()
                .map(|s| SpeciesId::new(*s).unwrap())
                .collect(),
        }
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
        MineralReaction::HydrogenCyanideDissociation
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
        .filter(|(_, _, amount)| *amount > 0.0)
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

fn member(answer: &AcidBaseSpeciation, name: &str) -> f64 {
    answer.families[0]
        .forms
        .iter()
        .find(|s| s.species.as_str() == name)
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

fn relative(actual: f64, expected: f64) {
    close(actual, expected, 2e-9 * expected.abs());
}

fn constants() -> (f64, f64) {
    // Independently transcribed primary formation-G reaction sums [J/mol].
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    ((-52_700.0 / rt).exp(), (-79_885.0 / rt).exp())
}

fn check_closure(answer: &AcidBaseSpeciation, total: f64) {
    let hcn = member(answer, "HCN(aq)");
    let cn = member(answer, "CN[-](aq)");
    close(hcn + cn, total, 2e-15 * total.max(1.0));
    let charge = answer.proton_amount.value() - answer.hydroxide_amount.value() - cn
        + answer
            .spectators
            .iter()
            .map(|s| s.amount.value() * f64::from(s.charge))
            .sum::<f64>();
    close(charge, 0.0, 1e-12);
    close(answer.charge_balance_residual.value(), charge, 1e-16);
    assert!(answer.max_log_mass_action_residual.value() < 1e-10);
    let (ka, kw) = constants();
    relative(
        answer.proton_activity.value() * answer.hydroxide_activity.value(),
        kw,
    );
    if total > 0.0 {
        relative(answer.proton_activity.value() * cn / hcn, ka);
    }
    assert!(answer.iterations <= 128);
}

#[test]
fn g1_original_source_constants_and_pure_water_close_simultaneously() {
    let system = system();
    let acid = system.families()[0].steps()[0].standard();
    let water = system.water().standard();
    close(acid.molar().value(), 52_700.0, 1e-8);
    close(water.molar().value(), 79_885.0, 1e-8);
    close(
        -acid.log_equilibrium_constant().value() / core::f64::consts::LN_10,
        9.232640551245206,
        1e-12,
    );
    close(
        -water.log_equilibrium_constant().value() / core::f64::consts::LN_10,
        13.995246497841048,
        1e-12,
    );
    for standard in [acid, water] {
        assert_eq!(standard.reference_axis(), REFERENCE);
        close(standard.temperature().value(), 298.15, 0.0);
        close(standard.reference_pressure().value(), 100_000.0, 0.0);
        for term in standard.terms() {
            let pack = materials()
                .packs
                .iter()
                .find(|p| p.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
            assert!(term.answer.receipt.observation_backed);
        }
    }
    let answer = solve(&system, 0.0, 0.0, 0.0, 10.0);
    let (_, kw) = constants();
    close(answer.ph.value(), 6.997623248920524, 1e-10);
    relative(answer.proton_activity.value(), kw.sqrt());
    relative(answer.hydroxide_activity.value(), kw.sqrt());
    assert_eq!(member(&answer, "HCN(aq)"), 0.0);
    assert_eq!(member(&answer, "CN[-](aq)"), 0.0);
    check_closure(&answer, 0.0);
}

// With no water contribution, the positive weak-acid/base root obeys
// x^2+k*x-k*C=0. In the complete system, F(x)=x-C*k/(k+x)=Kw/x.
// F'(x)>=1 implies x0 <= x <= x0+Kw/x0. This is an independent analytic
// bound including water, not an assertion that a water-free quadratic is exact.
fn quadratic_bounds(k: f64, total_activity: f64, kw: f64) -> (f64, f64) {
    let root = 2.0 * k * total_activity / (k + (k * k + 4.0 * k * total_activity).sqrt());
    (root, root + kw / root)
}

#[test]
fn g1_neutral_hcn_and_matching_cyanide_salt_obey_independent_quadratic_bounds() {
    let system = system();
    let (ka, kw) = constants();
    let basis = 10.0;
    let total = 0.1;
    let acid = solve(&system, total, 0.0, 0.0, basis);
    let salt = solve(&system, total, total, 0.0, basis);
    for (answer, k, actual) in [
        (&acid, ka, acid.proton_activity.value()),
        (&salt, kw / ka, salt.hydroxide_activity.value()),
    ] {
        let (low, high) = quadratic_bounds(k, total / basis, kw);
        let allowance = 1e-9 * low;
        assert!(
            actual >= low - allowance && actual <= high + allowance,
            "complete water-coupled root {actual} outside [{low}, {high}]"
        );
        check_closure(answer, total);
    }
    assert!(acid.ph.value() < 7.0);
    assert!(salt.ph.value() > 10.0);
    assert!(member(&acid, "CN[-](aq)") < member(&salt, "CN[-](aq)"));
    assert!(
        salt.spectators
            .iter()
            .any(|s| s.species.as_str() == "Na[+](aq)"
                && s.charge == 1
                && s.amount.value() == total)
    );
}

#[test]
fn g3_counterion_acid_base_changes_and_solvent_scaling_preserve_coupled_balances() {
    let system = system();
    let salt = solve(&system, 0.1, 0.1, 0.0, 10.0);
    let base = solve(&system, 0.1, 0.12, 0.0, 10.0);
    let acid = solve(&system, 0.1, 0.1, 0.05, 10.0);
    for answer in [&salt, &base, &acid] {
        check_closure(answer, 0.1);
    }
    assert!(base.ph.value() > salt.ph.value() && salt.ph.value() > acid.ph.value());
    assert!(member(&base, "CN[-](aq)") > member(&salt, "CN[-](aq)"));
    assert!(member(&salt, "CN[-](aq)") > member(&acid, "CN[-](aq)"));
    assert_eq!(acid.spectators.len(), 2);
    let scaled = solve(&system, 0.2, 0.2, 0.1, 20.0);
    close(scaled.ph.value(), acid.ph.value(), 1e-10);
    relative(member(&scaled, "HCN(aq)"), 2.0 * member(&acid, "HCN(aq)"));
    relative(
        member(&scaled, "CN[-](aq)"),
        2.0 * member(&acid, "CN[-](aq)"),
    );
    relative(
        scaled.proton_amount.value(),
        2.0 * acid.proton_amount.value(),
    );
    check_closure(&scaled, 0.2);

    // Independently choose spectator charge so h=Ka, including the water ions:
    // Na = C/2 + n_ref*(Kw/Ka - Ka). Both family forms must then be equal.
    let (ka, kw) = constants();
    let sodium = 0.05 + 10.0 * (kw / ka - ka);
    let halfway = solve(&system, 0.1, sodium, 0.0, 10.0);
    close(halfway.ph.value(), 9.232640551245206, 1e-10);
    relative(member(&halfway, "HCN(aq)"), 0.05);
    relative(member(&halfway, "CN[-](aq)"), 0.05);
    check_closure(&halfway, 0.1);
}
