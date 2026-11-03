//! Source-backed chemistry -> finite gold leach -> analytical-Au carbon contact.
//!
//! Run `gold_leach_to_carbon --materials PACK_DIR [--aerated]`; --help lists packs.
//! NBS supplies equilibrium constants, gold molar mass and standard reaction G.
//! Inventories, kinetic rates, carbon capacities and operating boxes are authored
//! illustrations, not measured ore data or calibrated CIP/CIL.
//! Carbon receives analytical Au only: no ligand/counterion adsorption chemistry
//! or complete speciated post-carbon liquor is inferred.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, PropertyKey,
    Provenance, QueryPoint, SelectionPolicy,
};
use fs_material::{
    aqueous_adsorption::{
        AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
        CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER, TEMPERATURE_AXIS,
        kinetics::{
            AdsorptionBatch, AdsorptionNumerics, ELAPSED_TIME_AXIS, LDF_LAW, LOADING_AXIS,
            LangmuirKinetics, RATE_PARAMETER, SorbentPopulation, VOLUME_AXIS,
        },
    },
    cyanide_leach,
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Frequency, Mass, MolarMass, Pressure, QuantitySpec,
    Temperature, Time, Volume,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            gas_liquid::IdealGasLiquidPartition,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
            },
        },
    },
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const CYANIDE: &str = "uncomplexed-cyanide";
const CARBONATE: &str = "inorganic-carbon";
const VOLUME: f64 = 0.01;
const TEMPERATURE: f64 = 298.15;
const INITIAL_GOLD: f64 = 0.001;
const INITIAL_CYANIDE: f64 = 0.01;
const OXYGEN_GAS_PACK: &str = "oxygen-gas-nbs1989-thermochemistry";
const PACKS: [(&str, &str); 11] = [
    ("hydrogen-cyanide-aqueous-nbs1989", "HCN(aq)"),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)"),
    ("cyanide-ion-aqueous-nbs1989", "CN[-](aq)"),
    ("hydroxide-ion-aqueous-nbs1989", "OH[-](aq)"),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)"),
    ("gold-solid-nbs1989-thermochemistry", "Au(s)"),
    ("dicyanoaurate-ion-aqueous-nbs1989", "Au(CN)2[-](aq)"),
    ("oxygen-aqueous-nbs1989-thermochemistry", "O2(aq)"),
    ("carbonic-component-aqueous-nbs1989", "H2CO3(aq)"),
    ("bicarbonate-ion-aqueous-nbs1989", "HCO3[-](aq)"),
    ("carbonate-ion-aqueous-nbs1989", "CO3[2-](aq)"),
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
    match index {
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        5 => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
        11 => flags.extend([("source-phase-gas", 1.0), ("source-pure-species", 1.0)]),
        _ => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if ![0, 7, 8].contains(&index) {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            }
            if index == 8 {
                flags.push(("source-carbonic-component-co2-plus-water", 1.0));
            }
        }
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn sources<'a>(
    indices: &[usize],
    ids: &'a [SpeciesId],
    packs: &'a [NormalizedPack],
    points: &'a [QueryPoint],
) -> Vec<FormationGibbsSource<'a>> {
    indices
        .iter()
        .map(|&i| FormationGibbsSource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
        .collect()
}

fn material(directory: &Path, name: &str) -> Result<NormalizedPack> {
    let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.pack_id() != name || pack.to_bytes() != bytes {
        return Err(format!("expected canonical source profile {name}").into());
    }
    println!("material={} hash={}", pack.pack_id(), pack.content_hash());
    Ok(pack)
}

fn oxygen_solubility(directory: &Path) -> Result<IdealGasLiquidPartition> {
    let gas = material(directory, OXYGEN_GAS_PACK)?;
    let aqueous = material(directory, PACKS[7].0)?;
    let gas_id = SpeciesId::new("O2(g)")?;
    let aqueous_id = SpeciesId::new("O2(aq)")?;
    let gas_at = point(11)?;
    let aqueous_at = point(7)?;
    let model = IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build()?,
        gas_id.clone(),
        aqueous_id.clone(),
        &[
            FormationGibbsSource {
                species: &gas_id, pack: &gas, at: &gas_at, reference_axis: REFERENCE,
            },
            FormationGibbsSource {
                species: &aqueous_id, pack: &aqueous, at: &aqueous_at, reference_axis: REFERENCE,
            },
        ],
    )?;
    let saturation = model.equilibrium_aqueous_amount(Pressure::new(21_000.0), Amount::new(10.0))?;
    println!(
        "oxygen_dissolution_standard_G_J_mol={:.3} lnK={:.12} \
        prescribed_pO2_Pa=21000 equilibrium_DO_mol={:.12} equilibrium_DO_mol_m3={:.12}",
        model.standard().molar().value(),
        model.standard().log_equilibrium_constant().value(),
        saturation.value(),
        saturation.value() / VOLUME,
    );
    Ok(model)
}

fn chemistry(directory: &Path) -> Result<(IdealAcidBaseSystem, MolarMass)> {
    let packs = PACKS
        .iter()
        .map(|(name, _)| material(directory, name))
        .collect::<Result<Vec<_>>>()?;
    let points = (0..PACKS.len()).map(point).collect::<Result<Vec<_>>>()?;
    let ids = PACKS
        .iter()
        .map(|(_, s)| SpeciesId::new(*s))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build()?,
        ids[4].clone(),
        ids[1].clone(),
        ids[3].clone(),
        &sources(&[4, 1, 3], &ids, &packs, &points),
    )?;
    let mut families = Vec::new();
    for (id, equations) in [
        (
            CYANIDE,
            vec![(MineralReaction::HydrogenCyanideDissociation, 0, 2)],
        ),
        (
            CARBONATE,
            vec![
                (MineralReaction::CarbonicComponentDissociation, 8, 9),
                (MineralReaction::BicarbonateDissociation, 9, 10),
            ],
        ),
    ] {
        let steps = equations
            .into_iter()
            .map(|(reaction, acid, base)| {
                AcidDissociation::from_sources(
                    reaction.build()?,
                    ids[acid].clone(),
                    ids[base].clone(),
                    ids[1].clone(),
                    &sources(&[acid, 1, base], &ids, &packs, &points),
                )
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        for step in &steps {
            println!(
                "family={id} step={} pKa={:.9}",
                step.reaction().id().as_str(),
                -step.standard().log_equilibrium_constant().value() / std::f64::consts::LN_10
            );
        }
        families.push(AcidFamily::new(id, steps)?);
    }
    let standard = MineralReaction::GoldCyanidationDissolvedOxygen
        .build()?
        .standard_gibbs_change(&sources(&[5, 2, 7, 4, 6, 3], &ids, &packs, &points))?;
    println!(
        "Elsner_standard_G_J_per_written_extent={:.3}; this is not a kinetic rate",
        standard.molar().value()
    );
    let mass = packs[5].claims().query_typed(
        &PropertyKey::with_quantity("molar-mass", QuantitySpec::dimensional(MolarMass::DIMS)),
        &points[5],
        SelectionPolicy::SingleClaimOnly,
    )?;
    packs[5].claims().verify_receipt(&mass.receipt)?;
    Ok((
        IdealAcidBaseSystem::new(water, families)?,
        MolarMass::new(mass.evidence.value.value),
    ))
}

fn card(
    law: &str,
    declaration: &str,
    parameters: impl IntoIterator<Item = (&'static str, f64, Dims)>,
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain("gold-leach-to-carbon-example-v1", declaration.as_bytes());
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: parameters
            .into_iter()
            .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: declaration.into(),
            license: "project-authored example".into(),
            artifact: Some(artifact),
        },
    }
}

fn carbon(
    id: &str,
    mass: f64,
    capacity: f64,
    half: f64,
    rate: f64,
    ph: f64,
) -> Result<SorbentPopulation> {
    let declaration = format!(
        "SYNTHETIC {id}: dry mass {mass} kg, Langmuir capacity {capacity} mol-Au/kg, \
        half concentration {half} mol-Au/m3, independent LDF rate {rate}/s, \
        T=298.15 K, V=0.01 m3, externally held pH={ph}; no measured ore/carbon calibration."
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
            TEMPERATURE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TEMPERATURE,
            TEMPERATURE,
        )
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), ph, ph);
    let equilibrium = AqueousLangmuir::from_card(
        ElementId::new("Au")?,
        card(
            AQUEOUS_LANGMUIR_LAW,
            &declaration,
            [
                (CAPACITY_PARAMETER, capacity, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
            ],
            domain.clone(),
        ),
    )?;
    let model = LangmuirKinetics::from_card(
        equilibrium,
        card(
            LDF_LAW,
            &declaration,
            [(RATE_PARAMETER, rate, Frequency::DIMS)],
            domain
                .with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    0.0,
                    500.0,
                )
                .with_quantity(
                    LOADING_AXIS,
                    QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                    0.0,
                    capacity,
                )
                .with_quantity(
                    VOLUME_AXIS,
                    QuantitySpec::dimensional(Volume::DIMS),
                    VOLUME,
                    VOLUME,
                ),
        ),
    )?;
    Ok(SorbentPopulation::new(
        id,
        model,
        Mass::new(mass),
        Amount::new(0.0),
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), ph)?,
    )?)
}

fn carbon_contact(
    dissolved: Amount,
    ph: f64,
    remaining_gold: Amount,
    molar_mass: MolarMass,
) -> Result<()> {
    let mut bath = AdsorptionBatch::new(
        vec![
            carbon("activated-carbon", 0.2, 0.01, 0.01, 0.04, ph)?,
            carbon("natural-carbon", 0.1, 0.005, 0.005, 0.10, ph)?,
        ],
        dissolved,
        Volume::new(VOLUME),
        Temperature::new(TEMPERATURE),
    )?;
    let result = bath.advance(
        Time::new(100.0),
        AdsorptionNumerics {
            maximum_step: Time::new(2.0),
            absolute_amount_error: Amount::new(1e-11),
            relative_amount_error: 1e-6,
            max_work: 2_000_000,
            ..AdsorptionNumerics::default()
        },
        || false,
    )?;
    let adsorbed: f64 = result.populations.iter().map(|s| s.adsorbed.value()).sum();
    let residual = remaining_gold.value() + result.dissolved.value() + adsorbed - INITIAL_GOLD;
    for state in &result.populations {
        println!(
            "carbon={} Au_mol={:.10} Au_mg={:.6} fraction_of_original_gold={:.8}",
            state.id,
            state.adsorbed.value(),
            state.adsorbed.value() * molar_mass.value() * 1e6,
            state.adsorbed.value() / INITIAL_GOLD
        );
    }
    println!(
        "after_carbon: unleached_solid_Au_mol={:.10} solution_Au_mol={:.10} \
        all_stage_Au_residual_mol={:.3e}",
        remaining_gold.value(),
        result.dissolved.value(),
        residual
    );
    if residual.abs() > 1e-10 * INITIAL_GOLD {
        return Err("leach-to-carbon analytical Au balance failed".into());
    }
    Ok(())
}

fn leach_model(
    chemistry: IdealAcidBaseSystem,
    oxygen_ceiling: f64,
) -> Result<cyanide_leach::GoldCyanidation> {
    let declaration = format!(
        "SYNTHETIC empirical gold law r=k*n_Au*c_free_CN*c_DO with k=1 \
        m6/(mol2 s), 298.15 K, V=0.01 m3, pH 8..13, time 0..200 s, \
        accessible Au 0..0.001 mol, free CN 0..1 mol/m3, DO 0..{oxygen_ceiling} mol/m3. \
        No measured ore, particle-size, liberation or reagent-consumption calibration."
    );
    println!("{declaration}");
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (cyanide_leach::PH_AXIS, Dims::NONE, 8.0, 13.0),
        (
            cyanide_leach::FREE_CYANIDE_AXIS,
            AmountConcentration::DIMS,
            0.0,
            1.0,
        ),
        (
            cyanide_leach::OXYGEN_AXIS,
            AmountConcentration::DIMS,
            0.0,
            oxygen_ceiling,
        ),
        (cyanide_leach::GOLD_AXIS, Amount::DIMS, 0.0, INITIAL_GOLD),
        (cyanide_leach::ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 200.0),
        (
            cyanide_leach::TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (cyanide_leach::VOLUME_AXIS, Volume::DIMS, VOLUME, VOLUME),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    Ok(cyanide_leach::GoldCyanidation::from_card(
        chemistry,
        CYANIDE,
        card(
            cyanide_leach::CYANIDATION_LAW,
            &declaration,
            [(
                cyanide_leach::RATE_PARAMETER,
                1.0,
                cyanide_leach::CyanidationRateConstant::DIMS,
            )],
            validity,
        ),
    )?)
}

fn print_leach(state: &cyanide_leach::CyanidationState) {
    let hcn = state
        .speciation
        .families
        .iter()
        .find(|family| family.id == CYANIDE)
        .and_then(|family| {
            family
                .forms
                .iter()
                .find(|form| form.species.as_str() == "HCN(aq)")
        })
        .expect("original cyanide family contains HCN");
    println!(
        "leach_s={:.1} pH={:.8} solid_Au_mol={:.10} dissolved_Au_mol={:.10} \
        uncomplexed_CN_mol={:.10} free_CN_mol={:.10} HCN_mol={:.10} \
        complex_bound_CN_mol={:.10} DO_mol={:.10} rate_Au_mol_s={:.6e}",
        state.elapsed.value(),
        state.speciation.ph.value(),
        state.accessible_gold.value(),
        state.gold_complex.value(),
        state.uncomplexed_cyanide.value(),
        state.free_cyanide.value(),
        hcn.amount.value(),
        2.0 * state.gold_complex.value(),
        state.dissolved_oxygen.value(),
        state.rate.value(),
    );
}

fn leach_charge(sodium: f64, inorganic_carbon: f64, oxygen: f64) -> Result<cyanide_leach::CyanidationCharge> {
    Ok(cyanide_leach::CyanidationCharge {
        accessible_gold: Amount::new(INITIAL_GOLD),
        uncomplexed_cyanide: Amount::new(INITIAL_CYANIDE),
        dissolved_oxygen: Amount::new(oxygen),
        gold_complex: Amount::new(0.0),
        amount_at_unit_activity: Amount::new(10.0),
        solution_volume: Volume::new(VOLUME),
        temperature: Temperature::new(TEMPERATURE),
        spectators: vec![SpectatorIon {
            species: SpeciesId::new("Na[+](aq)")?,
            charge: 1,
            amount: Amount::new(sodium),
        }],
        other_families: vec![(CARBONATE.into(), Amount::new(inorganic_carbon))],
    })
}

fn run(
    label: &str,
    model: cyanide_leach::GoldCyanidation,
    sodium: f64,
    inorganic_carbon: f64,
    oxygen: f64,
    molar_mass: MolarMass,
) -> Result<()> {
    println!(
        "\nscenario={label}; sodium_mol={sodium}; carbonate_family_mol={inorganic_carbon}; \
        initial_DO_mol={oxygen}; accessible_Au_mg={:.6}; finite_DO_gold_ceiling_mol={:.10}",
        INITIAL_GOLD * molar_mass.value() * 1e6,
        (4.0 * oxygen).min(INITIAL_GOLD)
    );
    let charge = leach_charge(sodium, inorganic_carbon, oxygen)?;
    let mut batch = cyanide_leach::GoldCyanidationBatch::new(
        model,
        charge,
        QueryPoint::new(),
        SpeciationOptions::default(),
    )?;
    print_leach(batch.state());
    for duration in [10.0, 40.0, 150.0] {
        let result = batch.advance(
            Time::new(duration),
            cyanide_leach::CyanidationNumerics {
                maximum_step: Time::new(2.0),
                absolute_gold_error: Amount::new(1e-11),
                relative_gold_error: Dimensionless::new(1e-6),
                ..cyanide_leach::CyanidationNumerics::default()
            },
            || false,
        )?;
        print_leach(&result.state);
        println!(
            "reaction_extent_mol={:.10} interval_dissolved_Au_mol={:.10} \
            Au_balance_mol={:.3e} CN_balance_mol={:.3e} DO_balance_mol={:.3e} \
            accepted={} rejected={} work={} accumulated_local_error_estimate_mol={:.3e}",
            result.reaction_extent.value(),
            result.interval_gold_dissolved.value(),
            result.conservation.gold.value(),
            result.conservation.cyanide.value(),
            result.conservation.dissolved_oxygen.value(),
            result.accepted_steps,
            result.rejected_steps,
            result.work,
            result.accumulated_error_estimate.value()
        );
    }
    let final_state = batch.state();
    println!(
        "Carbon handoff: only analytical dissolved Au enters the finite carbon bath; \
        pH is an externally held condition. No ligand/counterion transfer or post-carbon \
        cyanide speciation is asserted. Natural-carbon uptake is separate from activated-carbon recovery."
    );
    carbon_contact(
        final_state.gold_complex,
        final_state.speciation.ph.value(),
        final_state.accessible_gold,
        molar_mass,
    )
}

fn oxygen_transfer(
    solubility: IdealGasLiquidPartition,
    coefficient: f64,
) -> Result<cyanide_leach::aeration::OxygenTransfer> {
    use cyanide_leach::aeration;
    let declaration = format!(
        "SYNTHETIC oxygen-transfer card: kLa={coefficient}/s; fixed pO2=21000 Pa, \
        T=298.15 K, V=0.01 m3, actual DO and equilibrium target 0..0.4 mol/m3, \
        time 0..200 s. Maintained gas reservoir; no measured ore/slurry/aerator calibration."
    );
    println!("{declaration}");
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (aeration::OXYGEN_PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 21_000.0, 21_000.0),
        (cyanide_leach::OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 0.4),
        (cyanide_leach::ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 200.0),
        (cyanide_leach::TEMPERATURE_AXIS, Temperature::DIMS, TEMPERATURE, TEMPERATURE),
        (cyanide_leach::VOLUME_AXIS, Volume::DIMS, VOLUME, VOLUME),
        (aeration::SATURATION_CONCENTRATION_AXIS, AmountConcentration::DIMS, 0.0, 0.4),
        (cyanide_leach::ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 10.0, 10.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    Ok(aeration::OxygenTransfer::from_card(
        solubility,
        card(
            aeration::OXYGEN_TRANSFER_LAW,
            &declaration,
            [(aeration::KLA_PARAMETER, coefficient, Frequency::DIMS)],
            validity,
        ),
    )?)
}

fn run_aerated(
    label: &str,
    model: cyanide_leach::GoldCyanidation,
    solubility: IdealGasLiquidPartition,
    coefficient: f64,
    oxygen: f64,
    molar_mass: MolarMass,
) -> Result<()> {
    use cyanide_leach::aeration;
    println!(
        "\nscenario={label}; sodium_mol=0.02; carbonate_family_mol=0; \
        initial_DO_mol={oxygen}; kLa_per_s={coefficient}; prescribed_pO2_Pa=21000"
    );
    let mut batch = aeration::AeratedGoldCyanidationBatch::new(
        model,
        oxygen_transfer(solubility, coefficient)?,
        leach_charge(0.020, 0.0, oxygen)?,
        QueryPoint::new(),
        SpeciationOptions::default(),
    )?;
    println!(
        "finite_Au_CN_gold_ceiling_mol={:.10}; oxygen availability evolves through the gas boundary",
        batch.maximum_gold_dissolution().value()
    );
    print_leach(batch.state());
    for duration in [10.0, 40.0, 150.0] {
        let result = batch.advance(
            Time::new(duration),
            Pressure::new(21_000.0),
            &QueryPoint::new(),
            aeration::AerationNumerics {
                maximum_step: Time::new(2.0),
                absolute_gold_error: Amount::new(1e-11),
                absolute_oxygen_error: Amount::new(1e-11),
                absolute_transfer_error: Amount::new(1e-11),
                relative_error: Dimensionless::new(1e-6),
                ..aeration::AerationNumerics::default()
            },
            || false,
        )?;
        print_leach(&result.state);
        println!(
            "equilibrium_DO_mol={:.12} cumulative_O2_in_mol={:.12} \
            cumulative_O2_out_mol={:.12} net_O2_transfer_mol={:.12} \
            interval_O2_in_mol={:.12} interval_O2_out_mol={:.12} \
            reaction_O2_consumed_mol={:.12} Au_balance_mol={:.3e} CN_balance_mol={:.3e} \
            open_O2_balance_mol={:.3e} O2_roundoff_allowance_mol={:.3e} \
            accepted={} rejected={} work={} accumulated_local_estimates_mol={:?}",
            result.oxygen_saturation.value(),
            result.oxygen_in.value(),
            result.oxygen_out.value(),
            result.net_oxygen_transfer.value(),
            result.interval_oxygen_in.value(),
            result.interval_oxygen_out.value(),
            result.reaction_extent.value(),
            result.conservation.gold.value(),
            result.conservation.cyanide.value(),
            result.conservation.dissolved_oxygen.value(),
            result.oxygen_balance_allowance.value(),
            result.accepted_steps,
            result.rejected_steps,
            result.work,
            [
                result.accumulated_error_estimates.gold.value(),
                result.accumulated_error_estimates.dissolved_oxygen.value(),
                result.accumulated_error_estimates.oxygen_in.value(),
                result.accumulated_error_estimates.oxygen_out.value(),
            ],
        );
    }
    let state = batch.state();
    println!(
        "Carbon handoff: analytical dissolved Au, with pH externally held; \
        no ligand/counterion adsorption or post-carbon cyanide speciation is inferred."
    );
    carbon_contact(
        state.gold_complex,
        state.speciation.ph.value(),
        state.accessible_gold,
        molar_mass,
    )
}

fn main() -> Result<()> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if arguments.is_empty() || arguments.as_slice() == ["--help"] {
        println!("usage: gold_leach_to_carbon --materials PACK_DIR [--aerated]");
        println!("Compile each source with xtask matdb-pack into PACK_DIR/PACK_ID.fsmatpk:");
        for (pack, _) in PACKS {
            println!("  {pack}");
        }
        println!("Runs lower/higher alkalinity, finite oxygen limitation, and carbonate buffering.");
        println!("With --aerated, additionally compile {OXYGEN_GAS_PACK}.");
        println!(
            "Aerated mode compares a sealed control and two explicit kLa values, including zero-DO startup. \
            Gas pressure is prescribed; no finite gas-feed supply or aerator calibration is inferred."
        );
        println!(
            "Rates, inventories and carbon cards are authored illustrations, not plant predictions."
        );
        return Ok(());
    }
    if !(arguments.len() == 2 || (arguments.len() == 3 && arguments[2] == "--aerated"))
        || arguments[0] != "--materials"
    {
        return Err("usage: gold_leach_to_carbon --materials PACK_DIR [--aerated]".into());
    }
    let aerated = arguments.len() == 3;
    println!(
        "Fixed 10 kg solvent gives n_ref=10 mol at unit molal activity; \
        V=0.01 m3 is supplied separately for the empirical concentration rate law."
    );
    println!(
        "Water retains unit activity; finite Au and cyanide-family stocks are accounted. \
        The formal NBS H2CO3 component means CO2(aq)+H2O(l), not a separate hydration model."
    );
    let (chemistry, molar_mass) = chemistry(Path::new(&arguments[1]))?;
    if aerated {
        let solubility = oxygen_solubility(Path::new(&arguments[1]))?;
        let model = leach_model(chemistry, 0.4)?;
        for (label, coefficient, oxygen) in [
            ("sealed-control", 0.0, 0.00005),
            ("weak-aeration-zero-initial-DO", 0.0001, 0.0),
            ("stronger-aeration-zero-initial-DO", 0.01, 0.0),
        ] {
            run_aerated(label, model.clone(), solubility.clone(), coefficient, oxygen, molar_mass)?;
        }
        return Ok(());
    }
    let model = leach_model(chemistry, 0.1)?;
    for (label, sodium, carbon, oxygen) in [
        ("lower-alkalinity", 0.006, 0.0, 0.001),
        ("higher-alkalinity", 0.020, 0.0, 0.001),
        ("oxygen-limited", 0.020, 0.0, 0.00005),
        ("carbonate-buffered", 0.030, 0.010, 0.001),
    ] {
        run(label, model.clone(), sodium, carbon, oxygen, molar_mass)?;
    }
    Ok(())
}
