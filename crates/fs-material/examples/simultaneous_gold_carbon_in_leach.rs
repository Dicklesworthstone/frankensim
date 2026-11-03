//! Simultaneous gold dissolution and competing carbon capture in one finite bath.
//!
//! Run with --materials PACK_DIR; --help lists the nine original NBS profiles.
//! Source-backed HCN/water/Cu equilibria are solved while gold dissolves and
//! separately declared product carbon and natural carbon adsorb/desorb neutral
//! Na+/Au(CN)2- pairs. Dissolved oxygen is finite. All kinetic/isotherm inputs
//! below are explicit synthetic illustrations, not an ore or carbon calibration.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::cyanide_leach::{
    self, CyanidationCharge, GoldCyanidation, GoldCyanidationBatch, competition::CompetingMetals,
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

use fs_material::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
    CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER,
    flow::{LOCAL_LDF_LAW, LocalLangmuirTransfer, RetainedSorbent},
};
use fs_material::cyanide_leach::cil::{
    CarbonInLeachBatch, CarbonPopulation, CarbonRole, CilAdvance, CilNumerics, CilState,
    IonPairTransfer,
};
use fs_qty::{Frequency, Mass, chemistry::ElementId as AnalyticalElementId};

fn carbon_card(
    law: &str,
    declaration: &str,
    parameters: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "simultaneous-carbon-in-leach-example-v1",
        declaration.as_bytes(),
    );
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: parameters
            .iter()
            .map(|&(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: declaration.into(),
            license: "project-authored numerical example".into(),
            artifact: Some(artifact),
        },
    }
}

fn carbon(
    name: &str,
    role: CarbonRole,
    mass: f64,
    capacity: f64,
    half: f64,
    rate: f64,
) -> Result<CarbonPopulation> {
    let declaration = format!(
        "SYNTHETIC {name}: dry mass={mass} kg, saturation={capacity} mol-Au/kg, \
        half concentration={half} mol-Au/m3, local LDF relaxation={rate}/s; \
        aqueous Au0..0.2 mol/m3, pH8..13, T298.15 K, V0.01 m3. \
        No measured carbon, preg-robbing index, or ore calibration."
    );
    println!("{declaration}");
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            CONCENTRATION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            0.2,
        )
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TEMPERATURE,
            TEMPERATURE,
        )
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 8.0, 13.0);
    let equilibrium = AqueousLangmuir::from_card(
        AnalyticalElementId::new("Au")?,
        carbon_card(
            AQUEOUS_LANGMUIR_LAW,
            &declaration,
            &[
                (CAPACITY_PARAMETER, capacity, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
            ],
            domain.clone(),
        ),
    )?;
    let transfer = LocalLangmuirTransfer::from_card(
        equilibrium,
        carbon_card(
            LOCAL_LDF_LAW,
            &declaration,
            &[(
                fs_material::aqueous_adsorption::flow::RATE_PARAMETER,
                rate,
                Frequency::DIMS,
            )],
            domain
                .with_quantity(
                    fs_material::aqueous_adsorption::flow::LOADING_AXIS,
                    QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                    0.0,
                    capacity,
                )
                .with_quantity(
                    fs_material::aqueous_adsorption::flow::VOLUME_AXIS,
                    QuantitySpec::dimensional(Volume::DIMS),
                    LIQUID_VOLUME,
                    LIQUID_VOLUME,
                ),
        ),
    )?;
    Ok(CarbonPopulation::new(
        role,
        RetainedSorbent::new(
            name,
            transfer,
            Mass::new(mass),
            Amount::new(0.0),
            QueryPoint::new(),
        )?,
    )?)
}

fn print_state(scenario: &str, state: &CilState) -> Result<()> {
    let solution = &state.solution;
    let copper = solution
        .competition
        .as_ref()
        .ok_or("missing live copper equilibrium")?;
    let adsorbed: f64 = state.populations.iter().map(|p| p.adsorbed.value()).sum();
    let copper_amount: f64 = copper.metals[0]
        .members
        .iter()
        .map(|m| m.amount.value())
        .sum();
    let gold_residual =
        solution.accessible_gold.value() + solution.gold_complex.value() + adsorbed - INITIAL_GOLD;
    let cyanide_residual = solution.uncomplexed_cyanide.value()
        + copper.bound_ligand.value()
        + 2.0 * (solution.gold_complex.value() + adsorbed)
        - INITIAL_CYANIDE;
    let sodium_residual = state.dissolved_counterion.value() + adsorbed - INITIAL_CYANIDE;
    let oxygen_residual =
        solution.dissolved_oxygen.value() + 0.25 * solution.gold_dissolved.value() - INITIAL_OXYGEN;
    println!(
        "scenario={scenario} time_s={:.3} dissolved_from_ore_fraction={:.9} \
        solid_Au_mol={:.12e} aqueous_Au_mol={:.12e} all_carbon_Au_mol={:.12e}",
        state.elapsed.value(),
        solution.gold_dissolved.value() / INITIAL_GOLD,
        solution.accessible_gold.value(),
        solution.gold_complex.value(),
        adsorbed,
    );
    println!(
        "  pH={:.9} free_CN_mol={:.12e} HCN_plus_CN_mol={:.12e} \
        Cu_bound_CN_mol={:.12e} dissolved_Na_mol={:.12e} dissolved_O2_mol={:.12e}",
        solution.speciation.ph.value(),
        solution.free_cyanide.value(),
        solution.uncomplexed_cyanide.value(),
        copper.bound_ligand.value(),
        state.dissolved_counterion.value(),
        solution.dissolved_oxygen.value(),
    );
    for p in &state.populations {
        println!(
            "  carbon={} role={:?} Au_mol={:.12e} fraction_of_original_Au={:.9} \
            vacant_capacity_mol={:.12e} loading_mol_per_kg={:.12e}",
            p.id,
            p.role,
            p.adsorbed.value(),
            p.adsorbed.value() / INITIAL_GOLD,
            p.vacant_capacity.value(),
            p.loading.value(),
        );
    }
    println!(
        "  whole_bath_residuals_mol: Au={gold_residual:.3e} CN={cyanide_residual:.3e} \
        Na={sodium_residual:.3e} O2={oxygen_residual:.3e} Cu={:.3e}; \
        charge_residual_mol_equiv={:.3e}",
        copper_amount - DISSOLVED_COPPER,
        solution.speciation.charge_balance_residual.value(),
    );
    if gold_residual.abs() > 1e-11
        || cyanide_residual.abs() > 1e-9
        || sodium_residual.abs() > 1e-11
        || oxygen_residual.abs() > 1e-11
        || (copper_amount - DISSOLVED_COPPER).abs() > 1e-10
    {
        return Err("independent whole-bath CIL component account failed".into());
    }
    Ok(())
}

fn run(
    scenario: &str,
    with_product_carbon: bool,
    model: &GoldCyanidation,
    copper: &MetalLigandFamily,
) -> Result<CilAdvance> {
    let fresh = GoldCyanidationBatch::new_with_competing_metals(
        model.clone(),
        charge(DISSOLVED_COPPER)?,
        QueryPoint::new(),
        CompetingMetals {
            families: vec![copper.clone()],
            totals: vec![(COPPER.into(), Amount::new(DISSOLVED_COPPER))],
            numerics: LigandSpeciationOptions::default(),
        },
    )?;
    let mut populations = vec![carbon(
        "natural-carbon",
        CarbonRole::PregRobbing,
        0.1,
        0.005,
        0.005,
        0.10,
    )?];
    if with_product_carbon {
        populations.push(carbon(
            "activated-carbon",
            CarbonRole::ProductCarbon,
            0.2,
            0.01,
            0.01,
            0.04,
        )?);
    }
    let mut bath = CarbonInLeachBatch::from_batch(fresh, IonPairTransfer::sodium()?, populations)?;
    print_state(scenario, bath.state())?;
    let mut last = None;
    for duration in [10.0, 30.0, DURATION - 40.0] {
        let report = bath.advance(
            Time::new(duration),
            CilNumerics {
                maximum_step: Time::new(2.0),
                absolute_amount_error: Amount::new(2e-10),
                relative_amount_error: Dimensionless::new(1e-6),
                maximum_attempts: 8192,
                maximum_work: 10_000_000,
            },
            || false,
        )?;
        print_state(scenario, &report.state)?;
        println!(
            "  interval: dissolved_Au_mol={:.12e} product_capture_mol={:.12e} \
            preg_robbing_capture_mol={:.12e}; accepted={} rejected={} work={} \
            accumulated_local_amount_error_estimate_mol={:.3e}",
            report.interval_gold_dissolved.value(),
            report.interval_product_capture.value(),
            report.interval_preg_robbing_capture.value(),
            report.accepted_steps,
            report.rejected_steps,
            report.work,
            report.accumulated_error_estimate.value(),
        );
        last = Some(report);
    }
    last.ok_or_else(|| "no CIL intervals requested".into())
}

fn captured(report: &CilAdvance, name: &str) -> f64 {
    report
        .state
        .populations
        .iter()
        .find(|p| p.id == name)
        .map_or(0.0, |p| p.adsorbed.value())
}

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("usage: simultaneous_gold_carbon_in_leach --materials PACK_DIR");
        println!("Compile these original NBS profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for (name, _, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "Use: cargo run -p xtask -- matdb-pack --manifest \
            data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!(
            "Compare natural-carbon capture with simultaneous product-carbon capture \
            using finite DO, live Cu competition, and explicit Na+/Au(CN)2- pairs."
        );
        return Ok(());
    }
    if args.len() != 2 || args[0] != "--materials" {
        return Err("usage: simultaneous_gold_carbon_in_leach --materials PACK_DIR".into());
    }
    let (chemistry, copper) = load(Path::new(&args[1]))?;
    let model = kinetic_model(chemistry)?;
    println!("{DECLARATION}");
    println!(
        "SYNTHETIC charge: accessible Au={INITIAL_GOLD} mol, initial cyanide={INITIAL_CYANIDE} mol, \
        Na+={INITIAL_CYANIDE} mol, Cu(I)={DISSOLVED_COPPER} mol with equal Cl-, \
        dissolved O2={INITIAL_OXYGEN} mol. Sorbed pair=one Au(CN)2- plus one Na+; \
        uptake retains two cyanide ligands and desorption returns the complete pair."
    );
    let without = run("natural_carbon_only", false, &model, &copper)?;
    let with = run("simultaneous_cil", true, &model, &copper)?;
    let preg_without = captured(&without, "natural-carbon");
    let preg_with = captured(&with, "natural-carbon");
    let product = captured(&with, "activated-carbon");
    if !(product > 0.0 && preg_with < preg_without) {
        return Err("authored populations did not demonstrate competing CIL capture".into());
    }
    if (with.state.solution.gold_dissolved.value() - without.state.solution.gold_dissolved.value())
        .abs()
        > 5e-9
    {
        return Err("neutral pair capture spuriously changed the declared dissolution law".into());
    }
    println!(
        "SYNTHETIC comparison at {DURATION}s: product_capture_fraction={:.9}; \
        natural_carbon_capture_fraction_without_product={:.9}; with_product={:.9}.",
        product / INITIAL_GOLD,
        preg_without / INITIAL_GOLD,
        preg_with / INITIAL_GOLD,
    );
    println!(
        "The stated rate has no gold-product inhibition, so carbon changes where released gold \
        resides, while neutral pair transfer preserves free-CN/pH/Cu chemistry at fixed leach \
        progress. These synthetic captures are not a recovery forecast. This example excludes \
        oxygen replenishment, through-flow, carbon transfer between tanks, reagent dosing, \
        nonideal activities, cyanide destruction, and adsorption of copper complexes."
    );
    Ok(())
}
