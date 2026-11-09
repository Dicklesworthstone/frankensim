use super::*;
use fs_qty::AmountConcentration;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, StoichiometricMatrix,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    stoichiometry::{
        ConservedReaction,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, WaterAutoprotolysis,
                complexation::{IdealLigandSystem, LigandAssociation, MetalLigandFamily},
            },
        },
    },
};

// Reuse the original test-only source selectors and property-pack construction.
#[allow(dead_code)]
#[path = "../tests/fixture.rs"]
mod source_fixture;
use source_fixture::{CN, H, HC, OH, REFERENCE, WATER, id, pack, point};

const T: f64 = 298.15;
const V: f64 = 1e-3;
const CU: &str = "Cu[+](aq)";
const MONO: &str = "CuCN(aq)";
const BIS: &str = "Cu(CN)2[-](aq)";
const HS: &str = "HS[-](aq)";
const H2S: &str = "H2S(aq)";
const SOLID: &str = "Cu2S(s-nbs-alpha)";

fn model(free_copper: f64) -> CopperSulfidePrecipitation {
    // Authored, independently known state: pH6, free CN=HS=1e-4;
    // Ka(HCN)=1e-6, Ka(H2S)=1e-7 and successive binding constants 1e3.
    // These are constructed test energies, not measured copper/cyanide data.
    let rt = R * T;
    let k_phase = 1e-6 / (free_copper * free_copper * 1e-4);
    let records: Vec<_> = [
        (H, 1, 0.0),
        (OH, -1, -rt * 1e-14_f64.ln()),
        (WATER, 0, 0.0),
        (CN, -1, 0.0),
        (HC, 0, rt * 1e-6_f64.ln()),
        (HS, -1, 0.0),
        (H2S, 0, rt * 1e-7_f64.ln()),
        (CU, 1, 0.0),
        (MONO, 0, -rt * 1e3_f64.ln()),
        (BIS, -1, -rt * 1e6_f64.ln()),
        (SOLID, 0, -rt * k_phase.ln()),
    ]
    .into_iter()
    .map(|(name, charge, gibbs)| {
        let phase = if name == WATER {
            "source-phase-liquid"
        } else if name == SOLID {
            "source-phase-solid"
        } else {
            "source-phase-aqueous"
        };
        let query = point(T, phase, charge != 0);
        (id(name), query.clone(), pack(name, gibbs, &query))
    })
    .collect();
    let sources = |reaction: &ConservedReaction| {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let record = records
                    .iter()
                    .find(|record| record.0 == term.species)
                    .unwrap();
                FormationGibbsSource {
                    species: &record.0,
                    at: &record.1,
                    pack: &record.2,
                    reference_axis: REFERENCE,
                }
            })
            .collect::<Vec<_>>()
    };
    let reaction = MineralReaction::WaterDissociation.build().unwrap();
    let water = WaterAutoprotolysis::from_sources(
        reaction.clone(),
        id(WATER),
        id(H),
        id(OH),
        &sources(&reaction),
    )
    .unwrap();
    let families = [
        (
            "cyanide",
            MineralReaction::HydrogenCyanideDissociation,
            HC,
            CN,
        ),
        (
            "sulfide",
            MineralReaction::HydrogenSulfideDissociation,
            H2S,
            HS,
        ),
    ]
    .into_iter()
    .map(|(name, recipe, acid, base)| {
        let reaction = recipe.build().unwrap();
        let step = AcidDissociation::from_sources(
            reaction.clone(),
            id(acid),
            id(base),
            id(H),
            &sources(&reaction),
        )
        .unwrap();
        AcidFamily::new(name, vec![step]).unwrap()
    })
    .collect();
    let aqueous = IdealAcidBaseSystem::new(water, families).unwrap();
    let steps = [(CU, MONO, 0_i128), (MONO, BIS, 1_i128)]
        .into_iter()
        .map(|(parent, complex, count)| {
            let species = vec![id(parent), id(CN), id(complex)];
            let reaction =
                ReactionId::new(format!("authored-material-copper-step-{count}")).unwrap();
            let elements = ElementalMatrix::new(
                ["Cu", "C", "N"]
                    .into_iter()
                    .map(|name| ElementId::new(name).unwrap())
                    .collect(),
                species.clone(),
                vec![
                    vec![1, 0, 1],
                    vec![count, 1, count + 1],
                    vec![count, 1, count + 1],
                ],
            )
            .unwrap();
            let matrix = StoichiometricMatrix::new(
                species.clone(),
                vec![reaction.clone()],
                vec![vec![-1], vec![-1], vec![1]],
            )
            .unwrap();
            let charges = ChargeVector::new(species, vec![1 - count, -1, -count]).unwrap();
            let reaction = ConservedReaction::new(&elements, &matrix, &charges, reaction).unwrap();
            LigandAssociation::from_sources(
                reaction.clone(),
                id(parent),
                id(CN),
                id(complex),
                &sources(&reaction),
            )
            .unwrap()
        })
        .collect();
    let copper = MetalLigandFamily::new("copper", steps).unwrap();
    let ligand = IdealLigandSystem::new(aqueous, "cyanide", vec![copper]).unwrap();
    let reaction = MineralReaction::CopperSulfideHydrosulfidePrecipitation
        .build()
        .unwrap();
    CopperSulfidePrecipitation::from_sources(
        ligand,
        reaction.clone(),
        "copper",
        "sulfide",
        id(SOLID),
        &sources(&reaction),
    )
    .unwrap()
}

fn material(free_copper: f64, solid: f64) -> ComplexedElectrolyteSolubilitySolution {
    let mono = 1e3 * 1e-4 * free_copper;
    let bis = 1e6 * 1e-4 * 1e-4 * free_copper;
    let sodium = 2e-4 + bis + 1e-8 - 1e-6 - free_copper;
    ComplexedElectrolyteSolubilitySolution::from_copper_sulfide(
        model(free_copper),
        SolventActivityBasis::new(18e-6, AmountConcentration::new(1000.0)).unwrap(),
        Volume::new(V),
        Amount::new(2e-4 + mono + 2.0 * bis),
        vec![(
            "copper".into(),
            Amount::new(free_copper + mono + bis + 2.0 * solid),
        )],
        vec![("sulfide".into(), Amount::new(11e-4 + solid))],
        vec![
            SpectatorIon {
                species: id("Na[+](aq)"),
                charge: 1,
                amount: Amount::new(sodium),
            },
            SpectatorIon {
                species: id("Cl[-](aq)"),
                charge: -1,
                amount: Amount::new(0.0),
            },
        ],
        CopperSulfideOptions::default(),
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-290),
        "{actual:e} != {expected:e}"
    );
}
fn particle(material: &ComplexedElectrolyteSolubilitySolution, name: &str) -> f64 {
    let index = material
        .species()
        .iter()
        .position(|species| species == &id(name))
        .unwrap();
    material.solution().solute_moles()[index]
}
fn check_balances(material: &ComplexedElectrolyteSolubilitySolution) {
    near(
        material.dissolved_copper().value() + 2.0 * material.solid_amount().value(),
        material.metal_totals()[0].1.value(),
        1e-12,
    );
    near(
        material.dissolved_sulfide().value() + material.solid_amount().value(),
        material.other_acid_totals()[0].1.value(),
        1e-12,
    );
    near(
        material.ligand().uncomplexed_ligand.value() + material.ligand().bound_ligand.value(),
        material.ligand_total().value(),
        2e-9,
    );
    near(
        particle(material, CU) + particle(material, MONO) + particle(material, BIS),
        material.dissolved_copper().value(),
        1e-14,
    );
}
fn snapshot(material: &ComplexedElectrolyteSolubilitySolution) -> Vec<u64> {
    let mut values = vec![
        material.solution().solvent_volume_m3().to_bits(),
        material.solid_amount().value().to_bits(),
        material.aqueous().ph.value().to_bits(),
    ];
    values.extend(
        material
            .solution()
            .solute_moles()
            .iter()
            .map(|amount| amount.to_bits()),
    );
    values
}

#[test]
fn actual_complexes_and_solid_match_an_independent_known_state() {
    let material = material(1e-4, 0.001);
    near(material.aqueous().ph.value(), 6.0, 2e-8);
    near(material.solid_amount().value(), 0.001, 3e-7);
    near(particle(&material, CU), 1e-4, 3e-7);
    near(particle(&material, MONO), 1e-5, 3e-7);
    near(particle(&material, BIS), 1e-6, 3e-7);
    near(material.ligand().free_ligand.value(), 1e-4, 3e-7);
    assert_eq!(material.species().len(), 11);
    assert_eq!(material.species().last(), Some(&id(MONO)));
    assert!(!material.species().contains(&id(SOLID)));
    assert!(!material.species().contains(&id(WATER)));
    let unique: std::collections::BTreeSet<_> = material.species().iter().collect();
    assert_eq!(unique.len(), material.species().len());
    // Independent particle total includes neutral CuCN once, never its bound
    // ligand as another particle; Cu2S contributes no osmotic pressure.
    near(
        material.solution().solute_moles().iter().sum(),
        0.00151202,
        3e-7,
    );
    near(
        material.solution().osmotic_pressure_pa().unwrap(),
        R * T * 0.00151202 / V,
        3e-7,
    );
    check_balances(&material);
}

#[test]
fn dilution_and_explicit_acid_dosing_reuse_sources_and_complete_stocks() {
    let original = material(1e-4, 0.001);
    let before = snapshot(&original);
    let diluted = original.with_solvent_volume_m3(2.0 * V).unwrap();
    assert!(diluted.solid_amount().value() < original.solid_amount().value());
    assert!(diluted.dissolved_copper().value() > original.dissolved_copper().value());
    assert_eq!(diluted.metal_totals(), original.metal_totals());
    assert_eq!(diluted.other_acid_totals(), original.other_acid_totals());
    assert_eq!(diluted.ligand_total(), original.ligand_total());
    assert_eq!(diluted.spectators(), original.spectators());
    assert_eq!(diluted.aqueous().amount_at_unit_activity, Amount::new(2.0));
    assert!(std::ptr::eq(diluted.species(), original.species()));
    assert_eq!(
        format!("{:?}", diluted.model().standard()),
        format!("{:?}", original.model().standard())
    );
    check_balances(&diluted);

    let spectators: Vec<_> = original
        .spectators()
        .iter()
        .map(|ion| Amount::new(ion.amount.value() + if ion.charge == -1 { 1e-4 } else { 0.0 }))
        .collect();
    let acidified = original
        .with_analytical_totals(
            original.ligand_total(),
            &[original.metal_totals()[0].1],
            &[original.other_acid_totals()[0].1],
            &spectators,
            || false,
        )
        .unwrap();
    assert!(acidified.aqueous().ph.value() < original.aqueous().ph.value());
    assert!(acidified.solid_amount().value() < original.solid_amount().value());
    assert!(std::ptr::eq(acidified.species(), original.species()));
    assert_eq!(
        format!("{:?}", acidified.model().standard()),
        format!("{:?}", original.model().standard())
    );
    check_balances(&acidified);
    assert_eq!(snapshot(&original), before);
}

#[test]
fn independent_trace_complexes_survive_nearly_complete_precipitation() {
    let material = material(1e-25, 1.0);
    assert_eq!(
        material.equilibrium().copper_total.value() - 2.0 * material.solid_amount().value(),
        0.0
    );
    assert!(material.dissolved_copper().value() > 0.0);
    near(particle(&material, CU), 1e-25, 3e-7);
    near(particle(&material, MONO), 1e-26, 3e-7);
    near(particle(&material, BIS), 1e-27, 3e-7);
    check_balances(&material);
}

#[test]
fn domain_layout_and_late_cancellation_refuse_without_partial_material() {
    let original = material(1e-4, 0.001);
    let before = snapshot(&original);
    assert!(original.with_solvent_volume_m3(0.0).is_err());
    assert!(
        original
            .with_analytical_totals(original.ligand_total(), &[], &[], &[], || false)
            .is_err()
    );
    let strong_acid: Vec<_> = original
        .spectators()
        .iter()
        .map(|ion| Amount::new(ion.amount.value() + if ion.charge == -1 { 0.1 } else { 0.0 }))
        .collect();
    assert!(
        original
            .with_analytical_totals(
                original.ligand_total(),
                &[original.metal_totals()[0].1],
                &[original.other_acid_totals()[0].1],
                &strong_acid,
                || false,
            )
            .is_err()
    ); // the source-selected final pH domain remains [2,10]
    let mut calls = 0;
    original
        .with_solvent_volume_m3_cancelled(2.0 * V, || {
            calls += 1;
            false
        })
        .unwrap();
    let mut observed = 0;
    assert!(
        original
            .with_solvent_volume_m3_cancelled(2.0 * V, || {
                observed += 1;
                observed == calls
            })
            .is_err()
    );
    assert_eq!(observed, calls);
    assert_eq!(snapshot(&original), before);
}
