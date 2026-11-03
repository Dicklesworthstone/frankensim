//! Continuous aerated gold leaching with explicitly retained carbon populations.
//!
//! Run --materials PACK_DIR; --help lists eleven original NBS profiles.
//! Feed composition, slurry flow and oxygen partial pressure change across four
//! operating phases, advanced in one-second windows. Carbon stays in the tank. The competing
//! carbon is not an entrained-ore or carbon-washout model. Local leach, kLa and
//! adsorption coefficients are explicit synthetic choices, not measured recovery.
#![deny(unsafe_code)]

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::{
    aqueous_adsorption::{
        AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
        CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER,
        flow::{LOCAL_LDF_LAW, LocalLangmuirTransfer, RetainedSorbent},
    },
    cyanide_leach::{
        self, CyanidationCharge,
        aeration::{self, OxygenTransfer},
        cil::{CarbonPopulation, CarbonRole, IonPairTransfer},
        competition::CompetingMetals,
        flow::{
            CyanidationCstr, CyanidationFeed, CyanidationInventory, FlowNumerics,
            LOCAL_CYANIDATION_LAW, LocalGoldCyanidation, SpectatorFeed,
            cil::{AeratedCilCstr, CilFlowAdvance},
        },
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Frequency, Mass, Pressure, QuantitySpec,
    Temperature, Time, Volume, VolumetricFlowRate,
    chemistry::ElementId as AnalyticalElementId,
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
            gas_liquid::IdealGasLiquidPartition,
        },
    },
};
use std::{error::Error, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const CYANIDE: &str = "uncomplexed-cyanide";
const COPPER: &str = "dissolved-copper-i";
const TEMPERATURE: f64 = 298.15;
const LIQUID_VOLUME: f64 = 0.01;
const ACTIVITY_AMOUNT: f64 = 10.0;
const INITIAL_GOLD: f64 = 0.001;
const INITIAL_CYANIDE: f64 = 0.01;
const INITIAL_OXYGEN: f64 = 0.0;
const DISSOLVED_COPPER: f64 = 0.002;
const INITIAL_PRESSURE: f64 = 21_000.0;
const PACKS: [(&str, &str, bool); 11] = [
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
    ("oxygen-aqueous-nbs1989-thermochemistry", "O2(aq)", false),
    ("oxygen-gas-nbs1989-thermochemistry", "O2(g)", false),
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
    } else if PACKS[index].1 == "O2(g)" {
        flags.extend([("source-phase-gas", 1.0), ("source-pure-species", 1.0)]);
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

fn load(
    directory: &Path,
) -> Result<(
    IdealAcidBaseSystem,
    MetalLigandFamily,
    IdealGasLiquidPartition,
)> {
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
    let oxygen = IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build()?,
        ids[10].clone(),
        ids[9].clone(),
        &sources(&[10, 9]),
    )?;
    println!(
        "oxygen_source_lnK={:.12}",
        oxygen.standard().log_equilibrium_constant().value()
    );
    Ok((
        IdealAcidBaseSystem::new(water, vec![AcidFamily::new(CYANIDE, vec![cyanide])?])?,
        MetalLigandFamily::new(COPPER, associations)?,
        oxygen,
    ))
}

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

fn local_model(chemistry: IdealAcidBaseSystem) -> Result<LocalGoldCyanidation> {
    let declaration = "SYNTHETIC local r_Au=k*n_Au*c_free_CN*c_DO, k=0.5 m6/(mol2 s); \
        T298.15K, V0.01m3, n_ref10mol, pH8..13, freeCN0..2mol/m3, DO0..1mol/m3, \
        accessibleAu0..0.002mol. Separate local law; no batch-age or ore calibration.";
    println!("{declaration}");
    let mut domain = ValidityDomain::unconstrained();
    for (axis, dims, lo, hi) in [
        (cyanide_leach::PH_AXIS, Dims::NONE, 8.0, 13.0),
        (
            cyanide_leach::FREE_CYANIDE_AXIS,
            AmountConcentration::DIMS,
            0.0,
            2.0,
        ),
        (
            cyanide_leach::OXYGEN_AXIS,
            AmountConcentration::DIMS,
            0.0,
            1.0,
        ),
        (cyanide_leach::GOLD_AXIS, Amount::DIMS, 0.0, 0.002),
        (
            cyanide_leach::TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (
            cyanide_leach::VOLUME_AXIS,
            Volume::DIMS,
            LIQUID_VOLUME,
            LIQUID_VOLUME,
        ),
        (
            cyanide_leach::ACTIVITY_REFERENCE_AXIS,
            Amount::DIMS,
            ACTIVITY_AMOUNT,
            ACTIVITY_AMOUNT,
        ),
        (cyanide_leach::TOTAL_CYANIDE_AXIS, Amount::DIMS, 0.0, 0.02),
    ] {
        domain = domain.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    Ok(LocalGoldCyanidation::from_card(
        chemistry,
        CYANIDE,
        carbon_card(
            LOCAL_CYANIDATION_LAW,
            declaration,
            &[(
                cyanide_leach::RATE_PARAMETER,
                0.5,
                cyanide_leach::CyanidationRateConstant::DIMS,
            )],
            domain,
        ),
    )?)
}

fn oxygen_transfer(solubility: IdealGasLiquidPartition) -> Result<OxygenTransfer> {
    let declaration = "SYNTHETIC maintained-reservoir kLa=0.04/s; \
        actual NBS O2(g)/O2(aq) equilibrium, no measured slurry or aerator calibration. \
        pO2 0..100000Pa, DO/saturation0..1mol/m3, elapsed0..200s, \
        T298.15K, V0.01m3, n_ref10mol.";
    println!("{declaration}");
    let mut domain = ValidityDomain::unconstrained();
    for (axis, dims, lo, hi) in [
        (
            aeration::OXYGEN_PARTIAL_PRESSURE_AXIS,
            Pressure::DIMS,
            0.0,
            100_000.0,
        ),
        (
            cyanide_leach::OXYGEN_AXIS,
            AmountConcentration::DIMS,
            0.0,
            1.0,
        ),
        (cyanide_leach::ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 200.0),
        (
            cyanide_leach::TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (
            cyanide_leach::VOLUME_AXIS,
            Volume::DIMS,
            LIQUID_VOLUME,
            LIQUID_VOLUME,
        ),
        (
            aeration::SATURATION_CONCENTRATION_AXIS,
            AmountConcentration::DIMS,
            0.0,
            1.0,
        ),
        (
            cyanide_leach::ACTIVITY_REFERENCE_AXIS,
            Amount::DIMS,
            ACTIVITY_AMOUNT,
            ACTIVITY_AMOUNT,
        ),
    ] {
        domain = domain.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    Ok(OxygenTransfer::from_card(
        solubility,
        carbon_card(
            aeration::OXYGEN_TRANSFER_LAW,
            declaration,
            &[(aeration::KLA_PARAMETER, 0.04, Frequency::DIMS)],
            domain,
        ),
    )?)
}

fn charge() -> Result<CyanidationCharge> {
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
                amount: Amount::new(DISSOLVED_COPPER),
            },
        ],
    })
}

fn feed(cyanide: f64, sodium: f64) -> Result<CyanidationFeed> {
    Ok(CyanidationFeed {
        accessible_gold: AmountConcentration::new(0.1),
        cyanide_ligand: AmountConcentration::new(cyanide),
        dissolved_oxygen: AmountConcentration::new(0.0),
        gold_complex: AmountConcentration::new(0.0),
        other_families: vec![],
        spectators: vec![
            SpectatorFeed {
                species: SpeciesId::new("Na[+](aq)")?,
                charge: 1,
                concentration: AmountConcentration::new(sodium),
            },
            SpectatorFeed {
                species: SpeciesId::new("Cl[-](aq)")?,
                charge: -1,
                concentration: AmountConcentration::new(0.2),
            },
        ],
        competing_metals: vec![(COPPER.into(), AmountConcentration::new(0.2))],
    })
}
fn sodium(inventory: &CyanidationInventory) -> f64 {
    inventory
        .spectators
        .iter()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .map_or(0.0, |s| s.amount.value())
}
fn copper(inventory: &CyanidationInventory) -> f64 {
    inventory
        .competing_metals
        .iter()
        .map(|(_, n)| n.value())
        .sum()
}
fn retained(report: &CilFlowAdvance, id: &str) -> f64 {
    report
        .state
        .populations
        .iter()
        .find(|p| p.id == id)
        .map_or(0.0, |p| p.adsorbed.value())
}
fn gold_supplied(report: &CilFlowAdvance) -> f64 {
    INITIAL_GOLD
        + report.flow.cumulative_inlet.accessible_gold.value()
        + report.flow.cumulative_inlet.gold_complex.value()
}

fn print_report(scenario: &str, phase: &str, report: &CilFlowAdvance, show: bool) -> Result<()> {
    let state = &report.state;
    let current = &state.solution.inventory;
    let incoming = &report.flow.cumulative_inlet;
    let outgoing = &report.flow.cumulative_outlet;
    let adsorbed: f64 = state.populations.iter().map(|p| p.adsorbed.value()).sum();
    let gold = current.accessible_gold.value()
        + current.gold_complex.value()
        + adsorbed
        + outgoing.accessible_gold.value()
        + outgoing.gold_complex.value()
        - gold_supplied(report);
    let cyanide = current.cyanide_ligand.value()
        + 2.0 * (current.gold_complex.value() + adsorbed)
        + outgoing.cyanide_ligand.value()
        + 2.0 * outgoing.gold_complex.value()
        - INITIAL_CYANIDE
        - incoming.cyanide_ligand.value()
        - 2.0 * incoming.gold_complex.value();
    let sodium = sodium(current) + adsorbed + sodium(outgoing) - INITIAL_CYANIDE - sodium(incoming);
    let oxygen = current.dissolved_oxygen.value()
        + outgoing.dissolved_oxygen.value()
        + 0.25 * report.flow.cumulative_gold_dissolved.value()
        + report.cumulative_oxygen_out.value()
        - INITIAL_OXYGEN
        - incoming.dissolved_oxygen.value()
        - report.cumulative_oxygen_in.value();
    let copper = copper(current) + copper(outgoing) - DISSOLVED_COPPER - copper(incoming);
    if show {
        println!(
            "scenario={scenario} phase={phase} window_start_s={:.3} time_s={:.3} flow_m3_per_s={:.6e} pO2_Pa={:.3} \
        supplied_Au_mol={:.12e} cumulative_Au_dissolved_mol={:.12e}",
            report.flow.previous_elapsed.value(),
            state.elapsed.value(),
            report.flow.flow.value(),
            report.oxygen_partial_pressure.value(),
            gold_supplied(report),
            report.flow.cumulative_gold_dissolved.value(),
        );
        println!(
            "  tank: accessible_Au_mol={:.12e} aqueous_Au_mol={:.12e} \
        pH={:.9} freeCN_mol={:.12e} DO_mol={:.12e}; \
        cumulative_liquid_exports: accessible_Au_mol={:.12e} aqueous_Au_mol={:.12e}",
            current.accessible_gold.value(),
            current.gold_complex.value(),
            state.solution.speciation.ph.value(),
            state.solution.free_cyanide.value(),
            current.dissolved_oxygen.value(),
            outgoing.accessible_gold.value(),
            outgoing.gold_complex.value(),
        );
        for p in &state.populations {
            println!(
                "  retained_population={} role={:?} loaded_Au_mol={:.12e} vacant_mol={:.12e}",
                p.id,
                p.role,
                p.adsorbed.value(),
                p.vacant_capacity.value()
            );
        }
        println!(
            "  gas: interval_O2_in_mol={:.12e} interval_O2_out_mol={:.12e} \
        cumulative_O2_in_mol={:.12e} cumulative_O2_out_mol={:.12e}; \
        whole_history_residuals_mol: Au={gold:.3e} CN={cyanide:.3e} Na={sodium:.3e} O2={oxygen:.3e} Cu={copper:.3e}",
            report.oxygen_in.value(),
            report.oxygen_out.value(),
            report.cumulative_oxygen_in.value(),
            report.cumulative_oxygen_out.value(),
        );
        println!(
            "  accepted={} rejected={} work={} accumulated_local_error_estimate_mol={:.3e}",
            report.flow.accepted_steps,
            report.flow.rejected_steps,
            report.flow.work,
            report.flow.accumulated_error_estimate.value()
        );
    }
    for residual in [gold, cyanide, sodium, oxygen, copper] {
        if residual.abs() > 1e-9 {
            return Err("independent continuous-CIL whole-history component balance failed".into());
        }
    }
    if state
        .solution
        .speciation
        .charge_balance_residual
        .value()
        .abs()
        > 1e-9
    {
        return Err("continuous-CIL charge balance failed".into());
    }
    Ok(())
}

fn run(
    name: &str,
    product_carbon: bool,
    model: &LocalGoldCyanidation,
    metal: &MetalLigandFamily,
    oxygen: &OxygenTransfer,
) -> Result<CilFlowAdvance> {
    let liquid = CyanidationCstr::new_with_competing_metals(
        model.clone(),
        charge()?,
        QueryPoint::new(),
        CompetingMetals {
            families: vec![metal.clone()],
            totals: vec![(COPPER.into(), Amount::new(DISSOLVED_COPPER))],
            numerics: LigandSpeciationOptions::default(),
        },
    )?;
    let mut populations = vec![carbon(
        "retained-competing-carbon",
        CarbonRole::PregRobbing,
        0.1,
        0.005,
        0.005,
        0.10,
    )?];
    if product_carbon {
        populations.push(carbon(
            "retained-product-carbon",
            CarbonRole::ProductCarbon,
            0.2,
            0.01,
            0.01,
            0.04,
        )?);
    }
    let mut tank = AeratedCilCstr::from_tank(
        liquid,
        oxygen.clone(),
        IonPairTransfer::sodium()?,
        populations,
        Pressure::new(INITIAL_PRESSURE),
        &QueryPoint::new(),
    )?;
    let mut last = None;
    for (phase, duration, q, cn, na, pressure) in [
        ("startup", 30.0, 0.0001, 1.0, 1.0, INITIAL_PRESSURE),
        (
            "higher_cyanide_feed",
            40.0,
            0.0001,
            1.2,
            1.2,
            INITIAL_PRESSURE,
        ),
        (
            "slurry_stop_with_oxygen_stripping",
            20.0,
            0.0,
            1.2,
            1.2,
            0.0,
        ),
        (
            "restart_at_higher_flow",
            30.0,
            0.0002,
            1.0,
            1.0,
            INITIAL_PRESSURE,
        ),
    ] {
        let windows = duration as usize;
        let inlet = feed(cn, na)?;
        for window in 0..windows {
            // Keep the original phase's physical inputs and accuracy allocation,
            // while respecting the nested chemistry owner's per-call work cap.
            let report = tank.advance(
                Time::new(duration / windows as f64),
                VolumetricFlowRate::new(q),
                &inlet,
                Pressure::new(pressure),
                &QueryPoint::new(),
                FlowNumerics {
                    maximum_step: Time::new(2.0),
                    absolute_amount_error: Amount::new(2e-10 / windows as f64),
                    relative_amount_error: Dimensionless::new(1e-6 / windows as f64),
                    maximum_attempts: 8192,
                    maximum_work: 10_000_000,
                },
                || false,
            )?;
            print_report(name, phase, &report, window + 1 == windows)?;
            if q == 0.0
                && (report.flow.outlet.accessible_gold.value() != 0.0
                    || report.flow.outlet.gold_complex.value() != 0.0)
            {
                return Err("stopped slurry flow produced a liquid export".into());
            }
            last = Some(report);
        }
    }
    last.ok_or_else(|| "no continuous-CIL intervals".into())
}

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("usage: continuous_aerated_carbon_in_leach --materials PACK_DIR");
        println!("Compile these canonical profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for (name, _, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "cargo run -p xtask -- matdb-pack --manifest \
            data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!("All carbon is explicitly retained. Initial liquid and feed contain zero oxygen.");
        return Ok(());
    }
    if args.len() != 2 || args[0] != "--materials" {
        return Err("usage: continuous_aerated_carbon_in_leach --materials PACK_DIR".into());
    }
    let (chemistry, copper, solubility) = load(Path::new(&args[1]))?;
    let model = local_model(chemistry)?;
    let transfer = oxygen_transfer(solubility)?;
    println!(
        "SYNTHETIC operating schedule: initial Au0.001mol/CN0.01mol/Cu(I)0.002mol, \
        zero initial/feed O2, explicit liquid-volume basis. Every carbon population \
        remains in the tank; entrained preg-robbing carbon does not enter or wash out. \
        pO2=0 is an oxygen-free maintained reservoir that strips DO, not a gas-off switch."
    );
    let competing_only = run(
        "retained_competing_carbon_only",
        false,
        &model,
        &copper,
        &transfer,
    )?;
    let cil = run(
        "retained_product_and_competing_carbon",
        true,
        &model,
        &copper,
        &transfer,
    )?;
    let product = retained(&cil, "retained-product-carbon");
    let competing_with = retained(&cil, "retained-competing-carbon");
    let competing_without = retained(&competing_only, "retained-competing-carbon");
    if !(product > 0.0 && competing_with < competing_without) {
        return Err("authored retained populations did not demonstrate competing capture".into());
    }
    if cil.flow.cumulative_gold_dissolved.value() <= 0.0
        || cil.cumulative_oxygen_in.value() <= 0.0
        || cil.cumulative_oxygen_out.value() <= 0.0
    {
        return Err("continuous-CIL oxygen supply/stripping demonstration failed".into());
    }
    if (cil.flow.cumulative_gold_dissolved.value()
        - competing_only.flow.cumulative_gold_dissolved.value())
    .abs()
        > 1e-7
    {
        return Err("neutral pair capture changed the declared dissolution-only rate law".into());
    }
    println!(
        "SYNTHETIC final comparison: supplied_Au_mol={:.12e}, product_carbon_Au_mol={product:.12e}, \
        competing_carbon_Au_with_product_mol={competing_with:.12e}, without_product_mol={competing_without:.12e}.",
        gold_supplied(&cil),
    );
    println!(
        "Loaded carbon is retained in-circuit inventory; this example exports no recovered product. \
        All reported captures depend on authored coefficients and are not an ore-recovery forecast. \
        No carbon harvesting, entrained-carbon advection, tank cascade, gas-feed capacity, \
        solids slip, variable volume, cyanide destruction, nonideal activities or heat feedback is inferred."
    );
    Ok(())
}
