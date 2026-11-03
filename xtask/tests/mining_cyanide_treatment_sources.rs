//! Real NBS sources, finite peroxide stages and a separately declared ideal
//! acid/base solution. Requested extents do not infer treatment rates or time.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource, ReactionBalanceError,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
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
const COMPONENT: &str = "source-carbonic-component-co2-plus-water";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const ROWS: &[(&str, &str, &str, bool)] = &[
    (
        "H2O2(aq)",
        "hydrogen-peroxide-aqueous-nbs1989",
        "source-phase-aqueous",
        false,
    ),
    (
        "OCN[-](aq)",
        "cyanate-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "NH3(aq)",
        "ammonia-aqueous-nbs1989",
        "source-phase-aqueous",
        false,
    ),
    (
        "NH4[+](aq)",
        "ammonium-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "CN[-](aq)",
        "cyanide-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "H2O(l)",
        "water-liquid-nbs1989-thermochemistry",
        "source-phase-liquid",
        false,
    ),
    (
        "CO3[2-](aq)",
        "carbonate-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "O2(g)",
        "oxygen-gas-nbs1989-thermochemistry",
        "source-phase-gas",
        false,
    ),
    (
        "O2(aq)",
        "oxygen-aqueous-nbs1989-thermochemistry",
        "source-phase-aqueous",
        false,
    ),
    (
        "H+(aq)",
        "hydrogen-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "OH[-](aq)",
        "hydroxide-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "HCN(aq)",
        "hydrogen-cyanide-aqueous-nbs1989",
        "source-phase-aqueous",
        false,
    ),
    (
        "H2CO3(aq)",
        "carbonic-component-aqueous-nbs1989",
        "source-phase-aqueous",
        false,
    ),
    (
        "HCO3[-](aq)",
        "bicarbonate-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
];

struct Source {
    species: SpeciesId,
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

fn point(species: &str, phase: &str, ion: bool) -> QueryPoint {
    let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            298.15,
        )
        .unwrap()
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        (phase, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if phase == "source-phase-aqueous" {
        flags.push(("source-aqueous-unit-activity-infinite-dilution", 1.0));
        if ion {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    } else {
        flags.push(("source-pure-species", 1.0));
    }
    if species == "H2CO3(aq)" {
        // This NBS row is the formal CO2(aq)+H2O(l) component.
        flags.push((COMPONENT, 1.0));
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
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
            "frankensim-cyanide-treatment-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        ROWS.iter()
            .map(|(species, slug, phase, ion)| {
                let target = out.join(format!("{slug}.fsmatpk"));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{slug}: {}",
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), *slug);
                assert_eq!(pack.to_bytes(), bytes);
                // Ammonium also retains its independently published point Cp.
                let expected_claims = if *species == "NH4[+](aq)" { 5 } else { 4 };
                assert_eq!(pack.claims().claim_count(), expected_claims);
                Source {
                    species: SpeciesId::new(*species).unwrap(),
                    pack,
                    at: point(species, phase, *ion),
                }
            })
            .collect()
    })
}

fn source(species: &SpeciesId) -> &'static Source {
    sources()
        .iter()
        .find(|source| source.species == *species)
        .unwrap()
}
fn heat_sources(
    species: impl IntoIterator<Item = SpeciesId>,
) -> Vec<FormationEnthalpySource<'static>> {
    species
        .into_iter()
        .map(|species| {
            let s = source(&species);
            FormationEnthalpySource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn gibbs_sources(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let s = source(&term.species);
            FormationGibbsSource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn inventory(network: &ReactionNetwork, stocks: &[(&str, f64)]) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|species| {
            let value = stocks
                .iter()
                .find(|(name, _)| *name == species.as_str())
                .map_or(0.0, |(_, n)| *n);
            (species.clone(), Amount::new(value))
        })
        .collect()
}
fn amount(stocks: &[(SpeciesId, Amount)], species: &str) -> f64 {
    stocks
        .iter()
        .find(|(name, _)| name.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn mass(stocks: &[(SpeciesId, Amount)]) -> f64 {
    stocks
        .iter()
        .map(|(species, n)| {
            let s = source(species);
            let m = s
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new("molar-mass", MolarMass::DIMS),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            n.value() * m.evidence.value.value
        })
        .sum()
}
fn treatment() -> ReactionNetwork {
    ReactionNetwork::new(
        [
            MineralReaction::CyanidePeroxideOxidation,
            MineralReaction::HydrogenPeroxideDecomposition,
            MineralReaction::CyanateHydrolysisAmmoniumCarbonate,
        ]
        .into_iter()
        .map(|recipe| recipe.build().unwrap())
        .collect(),
    )
    .unwrap()
}
fn extents(oxidation: f64, decomposition: f64, hydrolysis: f64) -> Vec<(ReactionId, Amount)> {
    [
        ("cyanide-peroxide-oxidation", oxidation),
        ("hydrogen-peroxide-decomposition", decomposition),
        ("cyanate-hydrolysis-ammonium-carbonate", hydrolysis),
    ]
    .into_iter()
    .map(|(id, n)| (ReactionId::new(id).unwrap(), Amount::new(n)))
    .collect()
}

#[test]
fn sixteen_literal_nbs_facts_retain_units_uncertainty_and_replayable_receipts() {
    // NBS corrected 1989 rows: g/mol, kJ/mol, kJ/mol, J/(mol K).
    for (index, (s, values)) in sources()[..4]
        .iter()
        .zip([
            [34.0148, -191.17, -134.03, 143.9],
            [42.0173, -146.0, -97.4, 106.7],
            [17.0307, -80.29, -26.50, 111.3],
            [18.0387, -132.51, -79.31, 113.4],
        ])
        .enumerate()
    {
        let entropy = if ROWS[index].3 {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        for (name, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, values[0] * 0.001),
            (
                "standard-molar-formation-enthalpy",
                ENERGY,
                values[1] * 1000.0,
            ),
            (
                "standard-molar-formation-gibbs-energy",
                ENERGY,
                values[2] * 1000.0,
            ),
            (entropy, ENTROPY, values[3]),
        ] {
            let answer = s
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            s.pack.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    // Preserve the additional point Cp already on main; it is not a T-law.
    let ammonium = &sources()[3];
    let cp = ammonium
        .pack
        .claims()
        .query_typed(
            &PropertyKey::new(
                "conventional-standard-molar-isobaric-heat-capacity",
                ENTROPY,
            ),
            &ammonium.at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    close(cp.evidence.value.value, 79.9);
    assert_eq!(cp.evidence.value.uncertainty, UncertaintyModel::Unstated);
    ammonium.pack.claims().verify_receipt(&cp.receipt).unwrap();
}

#[test]
fn temperature_phase_and_ionic_convention_cannot_be_relabelled() {
    for (index, s) in sources()[..4].iter().enumerate() {
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        for at in [
            QueryPoint::new(),
            s.at.clone()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    300.0,
                )
                .unwrap(),
            s.at.clone()
                .with(REFERENCE, 0.0)
                .unwrap()
                .with("source-elemental-reference-usgs-b2131", 1.0)
                .unwrap(),
            s.at.clone()
                .with("source-phase-aqueous", 0.0)
                .unwrap()
                .with("source-phase-gas", 1.0)
                .unwrap(),
            s.at.clone()
                .with("source-aqueous-unit-activity-infinite-dilution", 0.0)
                .unwrap(),
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        let forbidden_entropy = if ROWS[index].3 {
            let changed =
                s.at.clone()
                    .with("source-ionic-reference-hydrogen-zero", 0.0)
                    .unwrap();
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &changed, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new(forbidden_entropy, ENTROPY),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn five_source_reactions_retain_independent_enthalpy_entropy_and_direct_gibbs() {
    // Independent literal NBS sums per mole of each written integer equation.
    for (recipe, h, g, entropy) in [
        (
            MineralReaction::CyanidePeroxideOxidation,
            -391_260.0,
            -372_899.0,
            -61.39,
        ),
        (
            MineralReaction::CyanateHydrolysisAmmoniumCarbonate,
            -91_990.0,
            -35_462.0,
            -190.02,
        ),
        (
            MineralReaction::AmmoniumDissociation,
            52_220.0,
            52_810.0,
            -2.1,
        ),
        (
            MineralReaction::HydrogenPeroxideDecomposition,
            -189_320.0,
            -206_198.0,
            57.158,
        ),
        (
            MineralReaction::HydrogenPeroxideDecompositionDissolvedOxygen,
            -201_020.0,
            -189_798.0,
            -37.08,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let h_sources = heat_sources(reaction.terms().iter().map(|term| term.species.clone()));
        let caloric = reaction.standard_equilibrium(&h_sources).unwrap();
        let direct = reaction
            .standard_gibbs_change(&gibbs_sources(&reaction))
            .unwrap();
        close(caloric.enthalpy().molar().value(), h);
        close(caloric.molar_entropy().value(), entropy);
        close(caloric.molar_gibbs().value(), h - 298.15 * entropy);
        close(direct.molar().value(), g);
        // Tabulated G and caloric H-TS are not silently reconciled.
        assert!((caloric.molar_gibbs().value() - direct.molar().value()).abs() > 1.0);
        for term in caloric.enthalpy().terms() {
            source(&term.term.species)
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        for term in direct.terms() {
            source(&term.term.species)
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
    }
}

#[test]
fn finite_peroxide_and_cyanate_stages_close_historical_source_mass_and_heat() {
    let network = treatment();
    let initial = inventory(
        &network,
        &[("CN[-](aq)", 10.0), ("H2O2(aq)", 8.0), ("H2O(l)", 100.0)],
    );
    let h = heat_sources(network.species().iter().cloned());
    let oxidized = network
        .advance_with_standard_enthalpy(&initial, &extents(6.0, 1.0, 0.0), &h)
        .unwrap();
    for (species, expected) in [
        ("CN[-](aq)", 4.0),
        ("H2O2(aq)", 0.0),
        ("OCN[-](aq)", 6.0),
        ("H2O(l)", 108.0),
        ("O2(g)", 1.0),
        ("NH4[+](aq)", 0.0),
        ("CO3[2-](aq)", 0.0),
    ] {
        close(
            amount(&oxidized.inventory.final_inventory, species),
            expected,
        );
    }
    close(amount(&oxidized.inventory.consumed, "H2O2(aq)"), 8.0);
    close(oxidized.enthalpy.value(), -2_536_880.0);
    let hydrolyzed = network
        .advance_with_standard_enthalpy(
            &oxidized.inventory.final_inventory,
            &extents(0.0, 0.0, 6.0),
            &h,
        )
        .unwrap();
    for (species, expected) in [
        ("CN[-](aq)", 4.0),
        ("H2O2(aq)", 0.0),
        ("OCN[-](aq)", 0.0),
        ("H2O(l)", 96.0),
        ("O2(g)", 1.0),
        ("NH4[+](aq)", 6.0),
        ("CO3[2-](aq)", 6.0),
    ] {
        close(
            amount(&hydrolyzed.inventory.final_inventory, species),
            expected,
        );
    }
    close(hydrolyzed.enthalpy.value(), -551_940.0);
    close(
        oxidized.enthalpy.value() + hydrolyzed.enthalpy.value(),
        -3_088_820.0,
    );
    // Historical source masses include 18.0154 g/mol water; no row is replaced
    // with modern atomic weights to repair the independently closing ledger.
    close(mass(&initial), 2.333_837_4);
    close(mass(&initial), mass(&oxidized.inventory.final_inventory));
    close(mass(&initial), mass(&hydrolyzed.inventory.final_inventory));
    let mut reversed = initial.clone();
    reversed.reverse();
    let mut requests = extents(6.0, 1.0, 0.0);
    requests.reverse();
    let replay = network
        .advance_with_standard_enthalpy(&reversed, &requests, &h)
        .unwrap();
    assert_eq!(
        oxidized.inventory.final_inventory,
        replay.inventory.final_inventory
    );
    assert_eq!(
        oxidized.enthalpy.value().to_bits(),
        replay.enthalpy.value().to_bits()
    );
}

#[test]
fn shared_reagent_overdraw_same_stage_borrowing_water_and_bad_extents_are_refused() {
    let network = treatment();
    let initial = inventory(
        &network,
        &[("CN[-](aq)", 10.0), ("H2O2(aq)", 8.0), ("H2O(l)", 100.0)],
    );
    let original = initial.clone();
    let h = heat_sources(network.species().iter().cloned());
    assert!(
        matches!(network.advance_with_standard_enthalpy(&initial, &extents(7.0, 1.0, 0.0), &h),
        Err(ReactionNetworkError::InsufficientInventory { species, available, required })
        if species.as_str() == "H2O2(aq)" && available.value() == 8.0 && required.value() == 9.0)
    );
    assert!(
        matches!(network.advance_with_standard_enthalpy(&initial, &extents(6.0, 1.0, 6.0), &h),
        Err(ReactionNetworkError::InsufficientInventory { species, available, required })
        if species.as_str() == "OCN[-](aq)" && available.value() == 0.0 && required.value() == 6.0)
    );
    let dry = inventory(&network, &[("OCN[-](aq)", 6.0), ("H2O(l)", 11.0)]);
    assert!(
        matches!(network.advance_with_standard_enthalpy(&dry, &extents(0.0, 0.0, 6.0), &h),
        Err(ReactionNetworkError::InsufficientInventory { species, available, required })
        if species.as_str() == "H2O(l)" && available.value() == 11.0 && required.value() == 12.0)
    );
    for extent in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            network.advance_with_standard_enthalpy(&initial, &extents(extent, 0.0, 0.0), &h),
            Err(ReactionNetworkError::Balance(
                ReactionBalanceError::Invalid(_)
            ))
        ));
    }
    assert_eq!(initial, original);
    assert_eq!(amount(&dry, "H2O(l)"), 11.0);
}

fn acid_system() -> IdealAcidBaseSystem {
    let species = |id: &str| SpeciesId::new(id).unwrap();
    let water_reaction = MineralReaction::WaterDissociation.build().unwrap();
    let water = WaterAutoprotolysis::from_sources(
        water_reaction.clone(),
        species("H2O(l)"),
        species("H+(aq)"),
        species("OH[-](aq)"),
        &gibbs_sources(&water_reaction),
    )
    .unwrap();
    let step = |recipe: MineralReaction, acid: &str, base: &str| {
        let reaction = recipe.build().unwrap();
        AcidDissociation::from_sources(
            reaction.clone(),
            species(acid),
            species(base),
            species("H+(aq)"),
            &gibbs_sources(&reaction),
        )
        .unwrap()
    };
    let ammonium = step(
        MineralReaction::AmmoniumDissociation,
        "NH4[+](aq)",
        "NH3(aq)",
    );
    let pka = -ammonium.standard().log_equilibrium_constant().value() / 10.0_f64.ln();
    close(
        pka,
        52_810.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15 * 10.0_f64.ln()),
    );
    close(pka, 9.251_911_717_481_200);
    IdealAcidBaseSystem::new(
        water,
        vec![
            AcidFamily::new("total-ammonia", vec![ammonium]).unwrap(),
            AcidFamily::new(
                "inorganic-carbon",
                vec![
                    step(
                        MineralReaction::CarbonicComponentDissociation,
                        "H2CO3(aq)",
                        "HCO3[-](aq)",
                    ),
                    step(
                        MineralReaction::BicarbonateDissociation,
                        "HCO3[-](aq)",
                        "CO3[2-](aq)",
                    ),
                ],
            )
            .unwrap(),
            AcidFamily::new(
                "uncomplexed-cyanide",
                vec![step(
                    MineralReaction::HydrogenCyanideDissociation,
                    "HCN(aq)",
                    "CN[-](aq)",
                )],
            )
            .unwrap(),
        ],
    )
    .unwrap()
}
fn member(answer: &AcidBaseSpeciation, id: &str) -> f64 {
    answer
        .families
        .iter()
        .flat_map(|family| &family.forms)
        .find(|form| form.species.as_str() == id)
        .unwrap()
        .amount
        .value()
}

#[test]
fn coupled_ammonia_carbonate_and_cyanide_share_a_source_bound_charge_balance() {
    let model = acid_system();
    let spectators = [("Na[+](aq)", 1, 11.0), ("Cl[-](aq)", -1, 1.0)]
        .into_iter()
        .map(|(id, charge, n)| SpectatorIon {
            species: SpeciesId::new(id).unwrap(),
            charge,
            amount: Amount::new(n),
        })
        .collect::<Vec<_>>();
    let answer = model
        .speciate(
            &[
                ("total-ammonia".to_owned(), Amount::new(6.0)),
                ("inorganic-carbon".to_owned(), Amount::new(6.0)),
                ("uncomplexed-cyanide".to_owned(), Amount::new(4.0)),
            ],
            &spectators,
            Amount::new(1_000.0),
            SpeciationOptions {
                ph_tolerance: Dimensionless::new(1e-11),
                absolute_charge_tolerance: Amount::new(1e-14),
                relative_charge_tolerance: Dimensionless::new(1e-12),
                ..SpeciationOptions::default()
            },
            || false,
        )
        .unwrap();

    // Independent 60-digit decimal root of
    // 1000(h-Kw/h)+6h/(h+Kn)-6(K1*h+2*K1*K2)/(h*h+K1*h+K1*K2)
    // -4Kcn/(h+Kcn)+11-1=0, h=10^-pH; source literal delta G values below.
    assert!((answer.ph.value() - 9.916_420_261_204_192).abs() < 2e-10);
    assert!((member(&answer, "NH3(aq)") - 4.932_114_667_865_892).abs() < 2e-10);
    assert!((member(&answer, "NH4[+](aq)") - 1.067_885_332_134_108).abs() < 2e-10);
    assert!((member(&answer, "CN[-](aq)") - 3.313_674_499_171_000).abs() < 2e-10);
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    let [kn, k1, k2, kcn, kw] =
        [52_810.0, 36_310.0, 58_960.0, 52_700.0, 79_885.0].map(|g| (-g / rt).exp());
    let h = answer.proton_activity.value();
    let d = h * h + k1 * h + k1 * k2;
    // Independent polynomial weights at the SAME solved pH, not separately
    // buffered single-family calculations or the solver's log recurrence.
    for (id, expected) in [
        ("NH4[+](aq)", 6.0 * h / (h + kn)),
        ("NH3(aq)", 6.0 * kn / (h + kn)),
        ("H2CO3(aq)", 6.0 * h * h / d),
        ("HCO3[-](aq)", 6.0 * k1 * h / d),
        ("CO3[2-](aq)", 6.0 * k1 * k2 / d),
        ("HCN(aq)", 4.0 * h / (h + kcn)),
        ("CN[-](aq)", 4.0 * kcn / (h + kcn)),
    ] {
        close(member(&answer, id), expected);
    }
    for family in &answer.families {
        let expected = if family.id == "uncomplexed-cyanide" {
            4.0
        } else {
            6.0
        };
        close(family.total.value(), expected);
        close(
            family.forms.iter().map(|s| s.amount.value()).sum(),
            expected,
        );
        close(family.forms.iter().map(|s| s.fraction.value()).sum(), 1.0);
    }
    // C and N family totals close; ideal acid/base has implicit solvent/proton
    // exchange, so this is not finite-water H/O bookkeeping.
    close(
        member(&answer, "H2CO3(aq)")
            + member(&answer, "HCO3[-](aq)")
            + member(&answer, "CO3[2-](aq)")
            + member(&answer, "HCN(aq)")
            + member(&answer, "CN[-](aq)"),
        10.0,
    );
    close(
        member(&answer, "NH4[+](aq)")
            + member(&answer, "NH3(aq)")
            + member(&answer, "HCN(aq)")
            + member(&answer, "CN[-](aq)"),
        10.0,
    );
    for spectator in &spectators {
        let retained = answer
            .spectators
            .iter()
            .find(|s| s.species == spectator.species)
            .unwrap();
        assert_eq!(retained.charge, spectator.charge);
        assert_eq!(retained.amount, spectator.amount);
    }
    let charge = answer.proton_amount.value() - answer.hydroxide_amount.value()
        + member(&answer, "NH4[+](aq)")
        - member(&answer, "HCO3[-](aq)")
        - 2.0 * member(&answer, "CO3[2-](aq)")
        - member(&answer, "CN[-](aq)")
        + 11.0
        - 1.0;
    assert!(charge.abs() < 5e-11);
    close(charge, answer.charge_balance_residual.value());
    assert!(answer.max_log_mass_action_residual.value() < 1e-10);
    assert!((h * answer.hydroxide_activity.value() / kw - 1.0).abs() < 2e-9);
    for (ratio, k) in [
        (
            h * member(&answer, "NH3(aq)") / member(&answer, "NH4[+](aq)"),
            kn,
        ),
        (
            h * member(&answer, "HCO3[-](aq)") / member(&answer, "H2CO3(aq)"),
            k1,
        ),
        (
            h * member(&answer, "CO3[2-](aq)") / member(&answer, "HCO3[-](aq)"),
            k2,
        ),
        (
            h * member(&answer, "CN[-](aq)") / member(&answer, "HCN(aq)"),
            kcn,
        ),
    ] {
        assert!((ratio / k - 1.0).abs() < 2e-9);
    }
}
