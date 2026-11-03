//! Matched NBS carbonate components consumed with the cyanide acid/base family.
//! Formal H2CO3 means the source's CO2(aq)+H2O(l) component, not hydrated H2CO3.
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

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-cyanide";
const CARBON: &str = "inorganic-carbon";
const COMPONENT: &str = "source-carbonic-component-co2-plus-water";
const PACKS: [&str; 8] = [
    "hydrogen-cyanide-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "cyanide-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
    "carbonic-component-aqueous-nbs1989",
    "bicarbonate-ion-aqueous-nbs1989",
    "carbonate-ion-aqueous-nbs1989",
];
const SPECIES: [&str; 8] = [
    "HCN(aq)",
    "H+(aq)",
    "CN[-](aq)",
    "OH[-](aq)",
    "H2O(l)",
    "H2CO3(aq)",
    "HCO3[-](aq)",
    "CO3[2-](aq)",
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
        if index != 0 && index != 5 {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    if index == 5 {
        flags.push((COMPONENT, 1.0));
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
            "frankensim-carbonate-speciation-{}-{nonce}",
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
            points: (0..8).map(point).collect(),
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
    let step = |recipe: MineralReaction, acid: usize, base: usize| {
        AcidDissociation::from_sources(
            recipe.build().unwrap(),
            ids[acid].clone(),
            ids[base].clone(),
            ids[1].clone(),
            &sources(&[acid, 1, base]),
        )
        .unwrap()
    };
    let cyanide = AcidFamily::new(
        FAMILY,
        vec![step(MineralReaction::HydrogenCyanideDissociation, 0, 2)],
    )
    .unwrap();
    let carbonate = AcidFamily::new(
        CARBON,
        vec![
            step(MineralReaction::CarbonicComponentDissociation, 5, 6),
            step(MineralReaction::BicarbonateDissociation, 6, 7),
        ],
    )
    .unwrap();
    IdealAcidBaseSystem::new(water, vec![cyanide, carbonate]).unwrap()
}

fn solve(
    model: &IdealAcidBaseSystem,
    carbon: f64,
    cyanide: f64,
    sodium: f64,
    chloride: f64,
    basis: f64,
) -> AcidBaseSpeciation {
    let spectators = [("Na[+](aq)", 1, sodium), ("Cl[-](aq)", -1, chloride)]
        .into_iter()
        .filter(|(_, _, n)| *n > 0.0)
        .map(|(id, charge, amount)| SpectatorIon {
            species: SpeciesId::new(id).unwrap(),
            charge,
            amount: Amount::new(amount),
        })
        .collect::<Vec<_>>();
    model
        .speciate(
            &[
                (FAMILY.to_owned(), Amount::new(cyanide)),
                (CARBON.to_owned(), Amount::new(carbon)),
            ],
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

fn member(answer: &AcidBaseSpeciation, id: &str) -> f64 {
    answer
        .families
        .iter()
        .flat_map(|f| &f.forms)
        .find(|s| s.species.as_str() == id)
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
fn constants() -> [f64; 4] {
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    [36_310.0, 58_960.0, 52_700.0, 79_885.0].map(|g| (-g / rt).exp())
}

// Independent two-step polynomial weights, not the solver's log recurrence.
fn fractions(ph: f64) -> ([f64; 3], f64) {
    let [k1, k2, kcn, _] = constants();
    let h = 10.0_f64.powf(-ph);
    let denominator = h * h + k1 * h + k1 * k2;
    (
        [
            h * h / denominator,
            k1 * h / denominator,
            k1 * k2 / denominator,
        ],
        kcn / (kcn + h),
    )
}
fn sodium_at_ph(ph: f64, carbon: f64, cyanide: f64, basis: f64) -> f64 {
    let (_, _, _, kw) = {
        let k = constants();
        (k[0], k[1], k[2], k[3])
    };
    let h = 10.0_f64.powf(-ph);
    let (alpha, beta) = fractions(ph);
    carbon * (alpha[1] + 2.0 * alpha[2]) + cyanide * beta + basis * (kw / h - h)
}
fn closure(answer: &AcidBaseSpeciation, carbon: f64, cyanide: f64) {
    let c0 = member(answer, "H2CO3(aq)");
    let c1 = member(answer, "HCO3[-](aq)");
    let c2 = member(answer, "CO3[2-](aq)");
    let hcn = member(answer, "HCN(aq)");
    let cn = member(answer, "CN[-](aq)");
    close(c0 + c1 + c2, carbon, 2e-15 * carbon.max(1.0));
    close(hcn + cn, cyanide, 2e-15 * cyanide.max(1.0));
    let charge =
        answer.proton_amount.value() - answer.hydroxide_amount.value() - c1 - 2.0 * c2 - cn
            + answer
                .spectators
                .iter()
                .map(|s| s.amount.value() * f64::from(s.charge))
                .sum::<f64>();
    close(charge, 0.0, 1e-12);
    close(charge, answer.charge_balance_residual.value(), 1e-16);
    assert!(answer.max_log_mass_action_residual.value() < 1e-10);
    let [k1, k2, kcn, kw] = constants();
    let h = answer.proton_activity.value();
    if carbon > 0.0 {
        relative(h * c1 / c0, k1);
        relative(h * c2 / c1, k2);
    }
    if cyanide > 0.0 {
        relative(h * cn / hcn, kcn);
    }
    relative(h * answer.hydroxide_activity.value(), kw);
}

#[test]
fn g0_g1_literal_nbs_components_and_two_dissociation_energies_retain_their_basis() {
    let m = materials();
    for (i, mass, h, g, s) in [
        (5, 62.0254, -699.65, -623.08, 187.4),
        (6, 61.0174, -691.99, -586.77, 91.2),
        (7, 60.0094, -677.14, -527.81, -56.9),
    ] {
        let entropy = if i == 5 {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        for (name, dims, value) in [
            ("molar-mass", Dims([0, 1, 0, 0, 0, -1]), mass * 0.001),
            (
                "standard-molar-formation-enthalpy",
                Dims([2, 1, -2, 0, 0, -1]),
                h * 1000.0,
            ),
            (
                "standard-molar-formation-gibbs-energy",
                Dims([2, 1, -2, 0, 0, -1]),
                g * 1000.0,
            ),
            (entropy, Dims([2, 1, -2, -1, 0, -1]), s),
        ] {
            let answer = m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(
                answer.evidence.value.value,
                value,
                2e-12 * value.abs().max(1.0),
            );
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            m.packs[i].claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    let key = PropertyKey::new(
        "standard-molar-formation-gibbs-energy",
        Dims([2, 1, -2, 0, 0, -1]),
    );
    // Generic neutral-aqueous admission cannot erase the formal CO2+water component.
    assert!(
        m.packs[5]
            .claims()
            .query_typed(&key, &point(0), SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    assert!(
        m.packs[5]
            .claims()
            .query_typed(
                &key,
                &m.points[5].clone().with(COMPONENT, 0.0).unwrap(),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
    assert!(
        m.packs[5]
            .claims()
            .query_typed(
                &key,
                &m.points[5]
                    .clone()
                    .with(REFERENCE, 0.0)
                    .unwrap()
                    .with("source-elemental-reference-usgs-b2131", 1.0)
                    .unwrap(),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );

    for (recipe, indices, h, s, g, caloric, pk) in [
        (
            MineralReaction::CarbonicComponentDissociation,
            [5, 1, 6],
            7_660.0,
            -96.2,
            36_310.0,
            36_342.03,
            6.361236782081848,
        ),
        (
            MineralReaction::BicarbonateDissociation,
            [6, 1, 7],
            14_850.0,
            -148.1,
            58_960.0,
            59_006.015,
            10.329345102493686,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let gs = sources(&indices);
        let hs = gs
            .iter()
            .map(|source| FormationEnthalpySource {
                species: source.species,
                pack: source.pack,
                at: source.at,
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>();
        let direct = reaction.standard_gibbs_change(&gs).unwrap();
        let thermal = reaction.standard_equilibrium(&hs).unwrap();
        close(direct.molar().value(), g, 1e-8);
        close(thermal.enthalpy().molar().value(), h, 1e-8);
        close(thermal.molar_entropy().value(), s, 1e-10);
        close(thermal.molar_gibbs().value(), caloric, 1e-8);
        close(
            -direct.log_equilibrium_constant().value() / core::f64::consts::LN_10,
            pk,
            1e-11,
        );
    }
}

#[test]
fn g1_two_source_families_share_one_ph_and_match_independent_polynomial_fractions() {
    let model = system();
    let carbon = 0.1;
    let cyanide = 0.01;
    let basis = 10.0; // 10 kg solvent times the 1 mol/kg reference, not a volume.
    for ph in [6.5, 8.5, 10.5] {
        let sodium = sodium_at_ph(ph, carbon, cyanide, basis);
        let answer = solve(&model, carbon, cyanide, sodium, 0.0, basis);
        close(answer.ph.value(), ph, 1e-10);
        let (alpha, beta) = fractions(ph);
        for (id, fraction) in ["H2CO3(aq)", "HCO3[-](aq)", "CO3[2-](aq)"]
            .into_iter()
            .zip(alpha)
        {
            relative(member(&answer, id), carbon * fraction);
        }
        relative(member(&answer, "CN[-](aq)"), cyanide * beta);
        relative(member(&answer, "HCN(aq)"), cyanide * (1.0 - beta));
        closure(&answer, carbon, cyanide);
    }
}

#[test]
fn g3_carbonate_buffer_reduces_the_same_acid_pulse_and_preserves_free_cyanide() {
    let model = system();
    let initial_ph = 10.0;
    let cyanide = 0.01;
    let basis = 10.0;
    let acid = 0.005; // Explicit mol of Cl- spectator from a strong-acid addition.
    let mut after = Vec::new();
    for carbon in [0.0, 0.1] {
        // Equal initial pH, cyanide inventory and solvent basis. Carbonate's
        // matching Na+ is explicit, rather than an unbalanced extra anion pool.
        let sodium = sodium_at_ph(initial_ph, carbon, cyanide, basis);
        let initial = solve(&model, carbon, cyanide, sodium, 0.0, basis);
        let final_state = solve(&model, carbon, cyanide, sodium, acid, basis);
        close(initial.ph.value(), initial_ph, 1e-10);
        assert!(final_state.ph.value() < initial.ph.value());
        closure(&final_state, carbon, cyanide);
        println!(
            "carbon={carbon} mol: pH {} -> {}, free CN- {} -> {} mol",
            initial.ph.value(),
            final_state.ph.value(),
            member(&initial, "CN[-](aq)"),
            member(&final_state, "CN[-](aq)")
        );
        after.push(final_state);
    }
    assert!(after[1].ph.value() > after[0].ph.value());
    assert!(member(&after[1], "CN[-](aq)") > member(&after[0], "CN[-](aq)"));
    let sodium = sodium_at_ph(initial_ph, 0.1, cyanide, basis);
    let scaled = solve(
        &model,
        0.2,
        2.0 * cyanide,
        2.0 * sodium,
        2.0 * acid,
        2.0 * basis,
    );
    close(scaled.ph.value(), after[1].ph.value(), 1e-10);
    for id in [
        "H2CO3(aq)",
        "HCO3[-](aq)",
        "CO3[2-](aq)",
        "HCN(aq)",
        "CN[-](aq)",
    ] {
        relative(member(&scaled, id), 2.0 * member(&after[1], id));
    }
    closure(&scaled, 0.2, 2.0 * cyanide);
}
