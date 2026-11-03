//! Actual matched NBS equilibrium sources with an explicitly authored rate law.
//!
//! Run `copper_competitive_gold_leach --materials PACK_DIR`; --help lists packs.
//! Compare identical finite gold/cyanide/oxygen charges with zero and finite
//! dissolved Cu(I). Every kinetic stage re-solves Cu complexes, free CN and pH.
//! Inventories and the kinetic coefficient are numerical illustrations, not an
//! ore assay, measured rate, accessible-copper estimate or recovery forecast.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::cyanide_leach::{
    self, CyanidationAdvance, CyanidationCharge, CyanidationNumerics, CyanidationState,
    GoldCyanidation, GoldCyanidationBatch, competition::CompetingMetals,
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, QuantitySpec, Temperature, Time, Volume,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
    stoichiometry::{
        ConservedReaction,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                complexation::{LigandAssociation, LigandSpeciationOptions, MetalLigandFamily},
            },
        },
    },
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const CYANIDE: &str = "uncomplexed-cyanide";
const COPPER: &str = "dissolved-copper-i";
const TEMPERATURE: f64 = 298.15;
const LIQUID_VOLUME: f64 = 0.01;
const ACTIVITY_AMOUNT: f64 = 10.0;
const INITIAL_GOLD: f64 = 0.001;
const INITIAL_CYANIDE: f64 = 0.01;
const INITIAL_OXYGEN: f64 = 0.0005;
const DISSOLVED_COPPER: f64 = 0.002;
const DURATION: f64 = 120.0;
const DECLARATION: &str = "SYNTHETIC r_Au=k*n_Au*c_free_CN*c_DO, k=0.5 m6/(mol2 s); \
    fixed T=298.15 K, V_liquid=0.01 m3, n_ref=10 mol; pH8..13, freeCN0..1 mol/m3, \
    DO0..0.05 mol/m3, accessibleAu0..0.001 mol, elapsed0..300 s. No measured ore calibration.";
const PACKS: [(&str, &str, bool); 9] = [
    ("hydrogen-cyanide-aqueous-nbs1989", "HCN(aq)", false),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)", true),
    ("cyanide-ion-aqueous-nbs1989", "CN[-](aq)", true),
    ("hydroxide-ion-aqueous-nbs1989", "OH[-](aq)", true),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)", false),
    ("copper-i-ion-aqueous-nbs1989", "Cu[+](aq)", true),
    (
        "dicyanocuprate-i-ion-aqueous-nbs1989",
        "Cu(CN)2[-](aq)",
        true,
    ),
    (
        "tricyanocuprate-i-ion-aqueous-nbs1989",
        "Cu(CN)3[2-](aq)",
        true,
    ),
    (
        "tetracyanocuprate-i-ion-aqueous-nbs1989",
        "Cu(CN)4[3-](aq)",
        true,
    ),
];

fn point(index: usize) -> Result<QueryPoint> {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TEMPERATURE,
        )?
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )?;
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if PACKS[index].1 == "H2O(l)" {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if PACKS[index].2 {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn copper_step(
    parent: SpeciesId,
    complex: SpeciesId,
    before: i128,
    after: i128,
) -> Result<ConservedReaction> {
    let species = vec![parent, SpeciesId::new("CN[-](aq)")?, complex];
    let elements = ElementalMatrix::new(
        vec![
            ElementId::new("Cu")?,
            ElementId::new("C")?,
            ElementId::new("N")?,
        ],
        species.clone(),
        vec![
            vec![1, 0, 1],
            vec![before, 1, after],
            vec![before, 1, after],
        ],
    )?;
    let id = ReactionId::new(format!("copper-i-{after}-cyanide-association"))?;
    let coefficients = StoichiometricMatrix::new(
        species.clone(),
        vec![id.clone()],
        vec![vec![-1], vec![before - after], vec![1]],
    )?;
    let charge = ChargeVector::new(species, vec![1 - before, -1, 1 - after])?;
    Ok(ConservedReaction::new(
        &elements,
        &coefficients,
        &charge,
        id,
    )?)
}

fn load(directory: &Path) -> Result<(IdealAcidBaseSystem, MetalLigandFamily)> {
    let packs = PACKS
        .iter()
        .map(|(name, _, _)| {
            let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
            let decoded = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
            if pack.pack_id() != *name || pack.to_bytes() != bytes {
                return Err(format!("expected canonical source profile {name}").into());
            }
            println!("material={} hash={}", pack.pack_id(), pack.content_hash());
            Ok(pack)
        })
        .collect::<Result<Vec<_>>>()?;
    let ids = PACKS
        .iter()
        .map(|(_, id, _)| SpeciesId::new(*id))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let points = (0..PACKS.len()).map(point).collect::<Result<Vec<_>>>()?;
    let sources = |indices: &[usize]| {
        indices
            .iter()
            .map(|&i| FormationGibbsSource {
                species: &ids[i],
                pack: &packs[i],
                at: &points[i],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>()
    };
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build()?,
        ids[4].clone(),
        ids[1].clone(),
        ids[3].clone(),
        &sources(&[4, 1, 3]),
    )?;
    let cyanide = AcidDissociation::from_sources(
        MineralReaction::HydrogenCyanideDissociation.build()?,
        ids[0].clone(),
        ids[2].clone(),
        ids[1].clone(),
        &sources(&[0, 1, 2]),
    )?;
    let mut associations = Vec::new();
    for (parent, complex, before, after) in [(5, 6, 0, 2), (6, 7, 2, 3), (7, 8, 3, 4)] {
        let reaction = copper_step(ids[parent].clone(), ids[complex].clone(), before, after)?;
        let step = LigandAssociation::from_sources(
            reaction,
            ids[parent].clone(),
            ids[2].clone(),
            ids[complex].clone(),
            &sources(&[parent, 2, complex]),
        )?;
        println!(
            "association={} source_deltaG_J_per_mol={:.6} source_lnK={:.12}",
            step.reaction().id().as_str(),
            step.standard().molar().value(),
            step.standard().log_equilibrium_constant().value()
        );
        associations.push(step);
    }
    Ok((
        IdealAcidBaseSystem::new(water, vec![AcidFamily::new(CYANIDE, vec![cyanide])?])?,
        MetalLigandFamily::new(COPPER, associations)?,
    ))
}

fn kinetic_model(chemistry: IdealAcidBaseSystem) -> Result<GoldCyanidation> {
    use cyanide_leach::*;
    let artifact = hash_domain("copper-competitive-gold-example-v1", DECLARATION.as_bytes());
    let axes = [
        (PH_AXIS, Dims::NONE, 8.0, 13.0),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS, 0.0, 1.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 0.05),
        (GOLD_AXIS, Amount::DIMS, 0.0, INITIAL_GOLD),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 300.0),
        (
            TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (VOLUME_AXIS, Volume::DIMS, LIQUID_VOLUME, LIQUID_VOLUME),
        (
            ACTIVITY_REFERENCE_AXIS,
            Amount::DIMS,
            ACTIVITY_AMOUNT,
            ACTIVITY_AMOUNT,
        ),
        (TOTAL_CYANIDE_AXIS, Amount::DIMS, 0.0, INITIAL_CYANIDE),
    ];
    let validity = axes.iter().fold(
        ValidityDomain::unconstrained(),
        |domain, &(axis, dims, lo, hi)| {
            domain.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi)
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(CYANIDATION_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(
            RATE_PARAMETER.into(),
            LawParameter {
                value: 0.5,
                dims: CyanidationRateConstant::DIMS,
            },
        )]
        .into_iter()
        .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: DECLARATION.into(),
            license: "project-authored numerical example".into(),
            artifact: Some(artifact),
        },
    };
    Ok(GoldCyanidation::from_card(chemistry, CYANIDE, card)?)
}

fn charge(copper: f64) -> Result<CyanidationCharge> {
    Ok(CyanidationCharge {
        accessible_gold: Amount::new(INITIAL_GOLD),
        uncomplexed_cyanide: Amount::new(INITIAL_CYANIDE),
        dissolved_oxygen: Amount::new(INITIAL_OXYGEN),
        gold_complex: Amount::new(0.0),
        amount_at_unit_activity: Amount::new(ACTIVITY_AMOUNT),
        solution_volume: Volume::new(LIQUID_VOLUME),
        temperature: Temperature::new(TEMPERATURE),
        other_families: vec![],
        spectators: vec![
            SpectatorIon {
                species: SpeciesId::new("Na[+](aq)")?,
                charge: 1,
                amount: Amount::new(INITIAL_CYANIDE),
            },
            SpectatorIon {
                species: SpeciesId::new("Cl[-](aq)")?,
                charge: -1,
                amount: Amount::new(copper),
            },
        ],
    })
}

fn print_state(scenario: &str, stage: &str, state: &CyanidationState, initial_copper: f64) {
    let competition = state
        .competition
        .as_ref()
        .expect("explicit competing-metal model");
    let copper = &competition.metals[0];
    let copper_sum: f64 = copper.members.iter().map(|m| m.amount.value()).sum();
    println!(
        "scenario={scenario} stage={stage} time_s={:.3} calculated_Au_recovery_fraction={:.9} \
        Au_remaining_mol={:.12e} Au_complex_mol={:.12e} rate_Au_mol_per_s={:.12e}",
        state.elapsed.value(),
        state.gold_dissolved.value() / INITIAL_GOLD,
        state.accessible_gold.value(),
        state.gold_complex.value(),
        state.rate.value(),
    );
    println!(
        "  pH={:.9} freeCN_mol={:.12e} HCN_plus_CN_mol={:.12e} Cu_bound_CN_mol={:.12e} \
        Au_bound_CN_mol={:.12e} dissolved_O2_mol={:.12e} conserved_Cu_mol={:.12e}",
        state.speciation.ph.value(),
        state.free_cyanide.value(),
        state.uncomplexed_cyanide.value(),
        competition.bound_ligand.value(),
        2.0 * state.gold_complex.value(),
        state.dissolved_oxygen.value(),
        copper_sum,
    );
    for member in &copper.members {
        println!(
            "  Cu_member={} amount_mol={:.12e} ligand_count={} activity={:.12e}",
            member.species.as_str(),
            member.amount.value(),
            member.ligand_count,
            member.activity.value()
        );
    }
    println!(
        "  residual_Au_mol={:.3e} residual_CN_mol={:.3e} residual_Cu_mol={:.3e} residual_O2_mol={:.3e} \
        charge_residual_mol_equiv={:.3e} equilibrium_max_abs_lnQ_over_K={:.3e}",
        state.accessible_gold.value() + state.gold_complex.value() - INITIAL_GOLD,
        state.uncomplexed_cyanide.value()
            + competition.bound_ligand.value()
            + 2.0 * state.gold_complex.value()
            - INITIAL_CYANIDE,
        copper_sum - initial_copper,
        state.dissolved_oxygen.value() + 0.25 * state.gold_dissolved.value() - INITIAL_OXYGEN,
        state.speciation.charge_balance_residual.value(),
        competition.max_log_mass_action_residual.value(),
    );
}

fn run(
    label: &str,
    copper: f64,
    model: &GoldCyanidation,
    family: &MetalLigandFamily,
) -> Result<CyanidationAdvance> {
    let mut bath = GoldCyanidationBatch::new_with_competing_metals(
        model.clone(),
        charge(copper)?,
        QueryPoint::new(),
        CompetingMetals {
            families: vec![family.clone()],
            totals: vec![(COPPER.into(), Amount::new(copper))],
            numerics: LigandSpeciationOptions::default(),
        },
    )?;
    print_state(label, "initial", bath.state(), copper);
    let advance = bath.advance(
        Time::new(DURATION),
        CyanidationNumerics {
            absolute_gold_error: Amount::new(1e-10),
            relative_gold_error: Dimensionless::new(1e-6),
            maximum_work: 10_000_000,
            ..CyanidationNumerics::default()
        },
        || false,
    )?;
    print_state(label, "final", &advance.state, copper);
    println!(
        "  accepted={} rejected={} work={} accumulated_local_Au_error_estimate_mol={:.3e}",
        advance.accepted_steps,
        advance.rejected_steps,
        advance.work,
        advance.accumulated_error_estimate.value()
    );
    Ok(advance)
}

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("usage: copper_competitive_gold_leach --materials PACK_DIR");
        println!("Compile these profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for (name, _, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "Use: cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!(
            "Runs zero/finite dissolved-Cu(I) comparisons with actual NBS equilibria and an authored kinetic law."
        );
        return Ok(());
    }
    if args.len() != 2 || args[0] != "--materials" {
        return Err("usage: copper_competitive_gold_leach --materials PACK_DIR".into());
    }
    let (chemistry, copper) = load(Path::new(&args[1]))?;
    let model = kinetic_model(chemistry)?;
    println!("{DECLARATION}");
    println!(
        "SYNTHETIC charge: Au={INITIAL_GOLD} mol, total ligand={INITIAL_CYANIDE} mol, \
        Na+={INITIAL_CYANIDE} mol, dissolved_O2={INITIAL_OXYGEN} mol, Cu(I)=0 or {DISSOLVED_COPPER} mol \
        with equal chloride counterion. Fixed activity amount={ACTIVITY_AMOUNT} mol and \
        liquid volume={LIQUID_VOLUME} m3 are independent declarations."
    );
    let no_copper = run("no_copper", 0.0, &model, &copper)?;
    let competing = run("finite_copper", DISSOLVED_COPPER, &model, &copper)?;
    if competing.state.gold_dissolved.value() >= no_copper.state.gold_dissolved.value() {
        return Err(
            "authored comparison did not demonstrate the resolved free-cyanide competition".into(),
        );
    }
    println!(
        "The equilibrium re-solves at every rate stage as gold binds cyanide; copper binding is \
        reversible and its dissolved component total remains fixed. Oxygen is a finite dissolved \
        stock. Ore copper accessibility, Cu(II) reduction, CuCN precipitation, cyanide destruction, \
        gas exchange, nonideal activities and measured ore recovery require additional supported inputs."
    );
    Ok(())
}
