//! Complete USGS-reference copper-sulfide leach and iron-regeneration balances.
//! Input H/G/S are literal B2131 pp.6, 8 and 11 values, not candidate-generated
//! oracles. Source compilation, canonical decoding and original physics owners
//! are exercised together; no claim of leach kinetics or plant calibration.
#![deny(unsafe_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, QueryPoint, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::{FormationGibbsSource, StandardReactionGibbs},
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
type Sources = BTreeMap<SpeciesId, (NormalizedPack, QueryPoint)>;

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn point(phase: &str, mineral: &str) -> QueryPoint {
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
        (phase, 1.0),
    ];
    if phase == "source-phase-aqueous" {
        flags.extend([
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    } else {
        flags.push(("source-pure-species", 1.0));
    }
    if !mineral.is_empty() {
        flags.push((mineral, 1.0));
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn sources() -> &'static Sources {
    static SOURCES: OnceLock<Sources> = OnceLock::new();
    SOURCES.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "frankensim-copper-sulfide-reactions-{}-{nonce}",
            std::process::id(),
        ));
        std::fs::create_dir(&directory).unwrap();
        let mut sources = Sources::new();
        for (species, slug, phase, mineral) in [
            (
                "CuFeS2(s)",
                "chalcopyrite-solid-usgs-b2131-thermochemistry",
                "source-phase-solid",
                "source-phase-chalcopyrite",
            ),
            (
                "Cu2S(s)",
                "chalcocite-solid-usgs-b2131-thermochemistry",
                "source-phase-solid",
                "source-phase-chalcocite",
            ),
            (
                "CuS(s)",
                "covellite-solid-usgs-b2131-thermochemistry",
                "source-phase-solid",
                "source-phase-covellite",
            ),
            (
                "Cu5FeS4(s)",
                "bornite-solid-usgs-b2131-thermochemistry",
                "source-phase-solid",
                "source-phase-bornite",
            ),
            (
                "S(s)",
                "sulfur-orthorhombic-solid-usgs-b2131",
                "source-phase-solid",
                "source-phase-sulfur-orthorhombic",
            ),
            (
                "Cu[2+](aq)",
                "copper-ii-ion-aqueous-usgs-b2131",
                "source-phase-aqueous",
                "",
            ),
            (
                "Fe[3+](aq)",
                "iron-iii-ion-aqueous-usgs-b2131",
                "source-phase-aqueous",
                "",
            ),
            (
                "Fe[2+](aq)",
                "iron-ii-ion-aqueous-usgs-b2131",
                "source-phase-aqueous",
                "",
            ),
            ("O2(g)", "oxygen-gas-usgs-b2131", "source-phase-gas", ""),
            (
                "H2O(l)",
                "water-liquid-usgs-b2131",
                "source-phase-liquid",
                "",
            ),
            (
                "H+(aq)",
                "hydrogen-ion-aqueous-usgs-b2131",
                "source-phase-aqueous",
                "",
            ),
        ] {
            let output = directory.join(format!("{slug}.fsmatpk"));
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
            sources.insert(
                SpeciesId::new(species).unwrap(),
                (pack, point(phase, mineral)),
            );
        }
        sources
    })
}

fn source<'a>(species: &'a SpeciesId, sources: &'a Sources) -> FormationEnthalpySource<'a> {
    let (pack, at) = sources.get(species).unwrap();
    FormationEnthalpySource {
        species,
        pack,
        at,
        reference_axis: REFERENCE,
    }
}

fn gibbs(reaction: &ConservedReaction, sources: &Sources) -> StandardReactionGibbs {
    let inputs = reaction
        .terms()
        .iter()
        .map(|term| {
            let input = source(&term.species, sources);
            FormationGibbsSource {
                species: input.species,
                pack: input.pack,
                at: input.at,
                reference_axis: input.reference_axis,
            }
        })
        .collect::<Vec<_>>();
    reaction.standard_gibbs_change(&inputs).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g1_four_sulfide_leaches_retain_source_heat_gibbs_entropy_and_activities() {
    use MineralReaction::*;
    for (recipe, expected_h, expected_g, expected_s) in [
        (FerricChalcopyriteLeach, 3_900.0, -123_000.0, 425.7),
        (FerricChalcociteLeach, 48_900.0, -73_800.0, 411.45),
        (FerricCovelliteLeach, 37_100.0, -26_200.0, 212.45),
        (FerricBorniteLeach, 110_600.0, -247_400.0, 1207.4),
    ] {
        let reaction = recipe.build().unwrap();
        let inputs = reaction
            .terms()
            .iter()
            .map(|term| source(&term.species, sources()))
            .collect::<Vec<_>>();
        let heat = reaction.standard_enthalpy_change(&inputs).unwrap();
        let standard = gibbs(&reaction, sources());
        let caloric = reaction.standard_equilibrium(&inputs).unwrap();
        close(heat.molar().value(), expected_h);
        close(
            heat.at_extent(Amount::new(1000.0)).unwrap().value(),
            expected_h * 1000.0,
        );
        close(standard.molar().value(), expected_g);
        close(caloric.molar_entropy().value(), expected_s);
        // H/G/S are separately rounded source observations. Preserve both routes;
        // especially bornite, where the caloric-derived G differs by 1.98631 kJ/mol.
        close(
            caloric.molar_gibbs().value(),
            expected_h - 298.15 * expected_s,
        );
        close(
            standard.log_equilibrium_constant().value(),
            -expected_g / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15),
        );
        for term in heat.terms() {
            let (pack, _) = sources().get(&term.term.species).unwrap();
            assert_eq!(
                term.answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        }
        for term in standard.terms() {
            let (pack, _) = sources().get(&term.term.species).unwrap();
            assert_eq!(
                term.answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        }
        let unit_activities = reaction
            .terms()
            .iter()
            .map(|term| (term.species.clone(), Dimensionless::new(1.0)))
            .collect::<Vec<_>>();
        close(
            standard
                .at_activities(&unit_activities)
                .unwrap()
                .molar
                .value(),
            expected_g,
        );
        let ferric = SpeciesId::new("Fe[3+](aq)").unwrap();
        let coefficient = reaction
            .terms()
            .iter()
            .find(|term| term.species == ferric)
            .unwrap()
            .coefficient;
        let mut dilute_ferric = unit_activities;
        dilute_ferric
            .iter_mut()
            .find(|(id, _)| *id == ferric)
            .unwrap()
            .1 = Dimensionless::new(0.1);
        let driving = standard.at_activities(&dilute_ferric).unwrap();
        // Lower ferric activity makes the forward driving force less favorable;
        // source reference conditions still remain exactly 298.15 K and 1 bar.
        close(
            driving.molar.value() - expected_g,
            coefficient as f64 * 0.1_f64.ln() * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15,
        );
        assert!(driving.molar.value() > expected_g);
    }
}

#[test]
fn g1_g3_leach_then_regeneration_closes_copper_iron_acid_oxygen_and_source_energy() {
    use MineralReaction::*;
    let regeneration = FerrousOxidationFerric.build().unwrap();
    let regeneration_g = gibbs(&regeneration, sources()).molar().value();
    close(regeneration_g, -181_000.0);
    for (recipe, mineral, ferric, copper, mineral_iron, sulfur, net_h, net_g) in [
        (
            FerricChalcopyriteLeach,
            "CuFeS2(s)",
            4.0,
            1.0,
            1.0,
            2.0,
            -402_900.0,
            -304_000.0,
        ),
        (
            FerricChalcociteLeach,
            "Cu2S(s)",
            4.0,
            2.0,
            0.0,
            1.0,
            -357_900.0,
            -254_800.0,
        ),
        (
            FerricCovelliteLeach,
            "CuS(s)",
            2.0,
            1.0,
            0.0,
            1.0,
            -166_300.0,
            -116_700.0,
        ),
        (
            FerricBorniteLeach,
            "Cu5FeS4(s)",
            12.0,
            5.0,
            1.0,
            4.0,
            -1_109_800.0,
            -790_400.0,
        ),
    ] {
        let leach = recipe.build().unwrap();
        let network = ReactionNetwork::new(vec![leach.clone(), regeneration.clone()]).unwrap();
        let inputs = network
            .species()
            .iter()
            .map(|id| source(id, sources()))
            .collect::<Vec<_>>();
        let initial = network
            .species()
            .iter()
            .map(|id| {
                let amount = match id.as_str() {
                    name if name == mineral => 1.0,
                    "Fe[3+](aq)" | "H+(aq)" => ferric,
                    "O2(g)" => ferric / 4.0,
                    _ => 0.0,
                };
                (id.clone(), Amount::new(amount))
            })
            .collect::<Vec<_>>();
        let requested = |leach_extent, regeneration_extent| {
            [
                (leach.id().clone(), Amount::new(leach_extent)),
                (regeneration.id().clone(), Amount::new(regeneration_extent)),
            ]
        };
        // The net cycle could conceal absent Fe(II). A simultaneous stage must
        // refuse borrowing the first reaction's products before they exist.
        assert!(matches!(
            network.advance_with_standard_enthalpy(
                &initial,
                &requested(1.0, ferric / 4.0),
                &inputs
            ),
            Err(ReactionNetworkError::InsufficientInventory { .. })
        ));
        let leached = network
            .advance_with_standard_enthalpy(&initial, &requested(1.0, 0.0), &inputs)
            .unwrap();
        let regenerated = network
            .advance_with_standard_enthalpy(
                &leached.inventory.final_inventory,
                &requested(0.0, ferric / 4.0),
                &inputs,
            )
            .unwrap();
        close(
            leached.enthalpy.value() + regenerated.enthalpy.value(),
            net_h,
        );
        close(
            gibbs(&leach, sources()).molar().value() + ferric / 4.0 * regeneration_g,
            net_g,
        );
        for (id, amount) in &regenerated.inventory.final_inventory {
            let expected = match id.as_str() {
                "Fe[3+](aq)" => ferric,
                "Fe[2+](aq)" => mineral_iron,
                "Cu[2+](aq)" => copper,
                "S(s)" => sulfur,
                "H2O(l)" => ferric / 2.0,
                _ => 0.0,
            };
            close(amount.value(), expected);
        }
        // Stoichiometric cycle closure does not erase either stage's receipts.
        assert_eq!(leached.reactions.len(), 2);
        assert_eq!(regenerated.reactions.len(), 2);
        for stage in [&leached, &regenerated] {
            for reaction in &stage.reactions {
                assert_eq!(reaction.standard.terms().len(), 5);
            }
        }
    }
}
