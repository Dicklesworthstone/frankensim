//! Actual NBS equilibrium sources with explicitly synthetic leach/transfer laws.
//!
//! Run `aerated_gold_leach --materials PACK_DIR`; --help lists all seven packs.
//! Compare a closed dissolved-oxygen inventory with finite headspace exchange,
//! then stop contact, meter a gas-phase dose and resume contact. Every oxygen
//! molecule is assigned to gas, liquid or original-equation consumption.
//! Rates, inventories and source boxes are authored numerical illustrations,
//! not measured ore data, air-flow settings or a calibrated recovery forecast.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::cyanide_leach::{
    self, GoldCyanidationBatch,
    gas_transfer::{
        AeratedCyanidationState, AeratedGoldCyanidation, AerationAdvance, AerationNumerics,
        GAS_VOLUME_AXIS, KLA_PARAMETER, OXYGEN_TRANSFER_LAW, OxygenTransfer, PARTIAL_PRESSURE_AXIS,
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dims, Frequency, Pressure, QuantitySpec, Temperature, Time,
    Volume,
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
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
            },
            gas_liquid::IdealGasLiquidPartition,
        },
    },
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const CYANIDE: &str = "uncomplexed-cyanide";
const TEMPERATURE: f64 = 298.15;
const LIQUID_VOLUME: f64 = 0.01;
const HEADSPACE_VOLUME: f64 = 0.0001;
const ACTIVITY_AMOUNT: f64 = 10.0;
const INITIAL_GOLD: f64 = 0.001;
const INITIAL_DO: f64 = 0.000001;
const INITIAL_GAS: f64 = 0.0001;
const OXYGEN_DOSE: f64 = 0.0002;
const PACKS: [(&str, &str); 7] = [
    ("hydrogen-cyanide-aqueous-nbs1989", "HCN(aq)"),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)"),
    ("cyanide-ion-aqueous-nbs1989", "CN[-](aq)"),
    ("hydroxide-ion-aqueous-nbs1989", "OH[-](aq)"),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)"),
    ("oxygen-gas-nbs1989-thermochemistry", "O2(g)"),
    ("oxygen-aqueous-nbs1989-thermochemistry", "O2(aq)"),
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
        5 => flags.extend([("source-phase-gas", 1.0), ("source-pure-species", 1.0)]),
        _ => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if [1, 2, 3].contains(&index) {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            }
        }
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn load(directory: &Path) -> Result<(IdealAcidBaseSystem, IdealGasLiquidPartition)> {
    let packs = PACKS
        .iter()
        .map(|(name, _)| {
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
        .map(|(_, id)| SpeciesId::new(*id))
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
    let phases = vec![ids[5].clone(), ids[6].clone()];
    let elements =
        ElementalMatrix::new(vec![ElementId::new("O")?], phases.clone(), vec![vec![2, 2]])?;
    let transfer_id = ReactionId::new("oxygen-gas-aqueous-partition")?;
    let coefficients = StoichiometricMatrix::new(
        phases.clone(),
        vec![transfer_id.clone()],
        vec![vec![-1], vec![1]],
    )?;
    let charge = ChargeVector::new(phases, vec![0, 0])?;
    let partition = IdealGasLiquidPartition::from_sources(
        ConservedReaction::new(&elements, &coefficients, &charge, transfer_id)?,
        ids[5].clone(),
        ids[6].clone(),
        &sources(&[5, 6]),
    )?;
    println!(
        "source_T_K={} source_p_ref_Pa={} oxygen_partition_lnK={:.12}",
        TEMPERATURE,
        partition.standard().reference_pressure().value(),
        partition.standard().log_equilibrium_constant().value(),
    );
    Ok((
        IdealAcidBaseSystem::new(water, vec![AcidFamily::new(CYANIDE, vec![cyanide])?])?,
        partition,
    ))
}

fn card(
    law: &str,
    declaration: &str,
    parameter: &str,
    value: f64,
    dims: Dims,
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain("aerated-gold-leach-example-v1", declaration.as_bytes());
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(parameter.into(), LawParameter { value, dims })]
            .into_iter()
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

fn domain(axes: &[(&str, Dims, f64, f64)]) -> ValidityDomain {
    axes.iter().fold(
        ValidityDomain::unconstrained(),
        |result, &(name, dims, lo, hi)| {
            result.with_quantity(name, QuantitySpec::dimensional(dims), lo, hi)
        },
    )
}

fn bath(chemistry: IdealAcidBaseSystem) -> Result<GoldCyanidationBatch> {
    use cyanide_leach::*;
    let declaration = "SYNTHETIC empirical r_Au=k*n_Au*c_free_CN*c_DO, k=0.5 m6/(mol2 s); \
        fixed T=298.15 K, V_liquid=0.01 m3, n_ref=10 mol; pH8..13, freeCN0..2 mol/m3, \
        DO0..0.1 mol/m3, accessibleAu0..0.001 mol, elapsed0..2000 s. No measured ore calibration.";
    println!("{declaration}");
    let validity = domain(&[
        (PH_AXIS, Dims::NONE, 8.0, 13.0),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS, 0.0, 2.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 0.1),
        (GOLD_AXIS, Amount::DIMS, 0.0, INITIAL_GOLD),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 2000.0),
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
    ]);
    let model = GoldCyanidation::from_card(
        chemistry,
        CYANIDE,
        card(
            CYANIDATION_LAW,
            declaration,
            RATE_PARAMETER,
            0.5,
            CyanidationRateConstant::DIMS,
            validity,
        ),
    )?;
    Ok(GoldCyanidationBatch::new(
        model,
        CyanidationCharge {
            accessible_gold: Amount::new(INITIAL_GOLD),
            uncomplexed_cyanide: Amount::new(0.01),
            dissolved_oxygen: Amount::new(INITIAL_DO),
            gold_complex: Amount::new(0.0),
            amount_at_unit_activity: Amount::new(ACTIVITY_AMOUNT),
            solution_volume: Volume::new(LIQUID_VOLUME),
            temperature: Temperature::new(TEMPERATURE),
            spectators: vec![SpectatorIon {
                species: SpeciesId::new("Na[+](aq)")?,
                charge: 1,
                amount: Amount::new(0.01),
            }],
            other_families: vec![],
        },
        QueryPoint::new(),
        SpeciationOptions::default(),
    )?)
}

fn transfer(partition: IdealGasLiquidPartition) -> Result<OxygenTransfer> {
    use cyanide_leach::*;
    let declaration = "SYNTHETIC liquid-film exchange J=kLa*(D_equilibrium(G)-D), kLa=0.01/s; \
        fixed T=298.15 K, V_liquid=0.01 m3, V_gas=0.0001 m3, n_ref=10 mol; \
        oxygen_partial_pressure0..100000 Pa, DO0..0.1 mol/m3. No measured slurry coefficient.";
    println!("{declaration}");
    let validity = domain(&[
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 0.1),
        (PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 0.0, 100_000.0),
        (
            TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (VOLUME_AXIS, Volume::DIMS, LIQUID_VOLUME, LIQUID_VOLUME),
        (
            GAS_VOLUME_AXIS,
            Volume::DIMS,
            HEADSPACE_VOLUME,
            HEADSPACE_VOLUME,
        ),
        (
            ACTIVITY_REFERENCE_AXIS,
            Amount::DIMS,
            ACTIVITY_AMOUNT,
            ACTIVITY_AMOUNT,
        ),
    ]);
    Ok(OxygenTransfer::from_card(
        partition,
        card(
            OXYGEN_TRANSFER_LAW,
            declaration,
            KLA_PARAMETER,
            0.01,
            Frequency::DIMS,
            validity,
        ),
    )?)
}

fn print_state(label: &str, state: &AeratedCyanidationState) {
    let s = &state.solution;
    println!(
        "stage={label} time_s={:.3} contact={} Au_dissolved_mol={:.12e} Au_remaining_mol={:.12e} \
        freeCN_mol={:.12e} uncomplexedCN_mol={:.12e} pH={:.9}",
        s.elapsed.value(),
        state.gas_contact_active,
        s.gold_dissolved.value(),
        s.accessible_gold.value(),
        s.free_cyanide.value(),
        s.uncomplexed_cyanide.value(),
        s.speciation.ph.value(),
    );
    println!(
        "  O2_gas_mol={:.12e} O2_liquid_mol={:.12e} oxygen_partial_pressure_Pa={:.9} \
        equilibrium_liquid_O2_mol={:.12e} cumulative_gas_to_liquid_O2_mol={:.12e}",
        state.gas_oxygen.value(),
        s.dissolved_oxygen.value(),
        state.oxygen_partial_pressure.value(),
        state.equilibrium_dissolved_oxygen.value(),
        state.cumulative_oxygen_transfer.value(),
    );
}

fn print_advance(label: &str, advance: &AerationAdvance) {
    print_state(label, &advance.state);
    let r = advance.conservation;
    println!(
        "  interval_Au_mol={:.12e} interval_transfer_O2_mol={:.12e} cumulative_dose_O2_mol={:.12e} \
        residual_Au_mol={:.3e} residual_CN_mol={:.3e} residual_liquid_O2_mol={:.3e} \
        residual_gas_O2_mol={:.3e} residual_total_O2_mol={:.3e} charge_residual_mol_equiv={:.3e}",
        advance.interval_gold_dissolved.value(),
        advance.interval_oxygen_transfer.value(),
        advance.cumulative_oxygen_dose.value(),
        r.gold.value(),
        r.cyanide.value(),
        r.liquid_oxygen.value(),
        r.gas_oxygen.value(),
        r.total_oxygen.value(),
        advance
            .state
            .solution
            .speciation
            .charge_balance_residual
            .value(),
    );
    println!(
        "  accepted={} rejected={} work={} accumulated_local_error_Au_equivalent_mol={:.3e}",
        advance.accepted_steps,
        advance.rejected_steps,
        advance.work,
        advance.accumulated_error_estimate.value(),
    );
}

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("usage: aerated_gold_leach --materials PACK_DIR");
        println!("Compile these profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for (name, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "Use: cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!(
            "Runs synthetic finite-headspace leach/stop/dose/restart scenarios using original NBS equilibrium data."
        );
        return Ok(());
    }
    if args.len() != 2 || args[0] != "--materials" {
        return Err("usage: aerated_gold_leach --materials PACK_DIR".into());
    }
    let (chemistry, partition) = load(Path::new(&args[1]))?;
    let bath = bath(chemistry)?;
    let initial_do_ceiling = bath.maximum_gold_dissolution().value();
    let transfer = transfer(partition)?;
    println!(
        "SYNTHETIC charge: Au={INITIAL_GOLD} mol, HCN+CN=0.01 mol, Na+=0.01 mol, \
        initial_liquid_O2={INITIAL_DO} mol, initial_gas_O2={INITIAL_GAS} mol. \
        Activity amount and liquid volume are independently declared."
    );
    let mut isolated = AeratedGoldCyanidation::from_batch(
        bath.clone(),
        transfer.clone(),
        Amount::new(INITIAL_GAS),
        Volume::new(HEADSPACE_VOLUME),
        QueryPoint::new(),
    )?;
    let no_contact =
        isolated.advance(Time::new(120.0), false, AerationNumerics::default(), || {
            false
        })?;
    print_advance("comparison_contact_off", &no_contact);
    let mut coupled = AeratedGoldCyanidation::from_batch(
        bath,
        transfer,
        Amount::new(INITIAL_GAS),
        Volume::new(HEADSPACE_VOLUME),
        QueryPoint::new(),
    )?;
    print_state("coupled_initial", coupled.state());
    let active = coupled.advance(Time::new(120.0), true, AerationNumerics::default(), || {
        false
    })?;
    print_advance("finite_headspace_exchange", &active);
    println!("Initial-DO-only stoichiometric Au ceiling [mol]: {initial_do_ceiling:.12e}");
    if no_contact.state.solution.gold_dissolved.value() > initial_do_ceiling
        || active.state.solution.gold_dissolved.value() <= initial_do_ceiling
    {
        return Err(
            "authored comparison did not demonstrate the finite-headspace oxygen contribution"
                .into(),
        );
    }
    let stopped = coupled.advance(Time::new(60.0), false, AerationNumerics::default(), || {
        false
    })?;
    print_advance("gas_contact_stopped", &stopped);
    coupled.dose_oxygen(Amount::new(OXYGEN_DOSE), || false)?;
    print_state("explicit_headspace_dose", coupled.state());
    println!(
        "Metered O2 dose [mol]: {OXYGEN_DOSE}; the dose does not directly alter liquid stocks or elapsed time."
    );
    let restarted = coupled.advance(Time::new(240.0), true, AerationNumerics::default(), || {
        false
    })?;
    print_advance("contact_restarted", &restarted);
    println!(
        "Final oxygen accounting [mol O2]: initial+metered={:.12e}; gas+liquid+reaction_consumed={:.12e}",
        INITIAL_DO + INITIAL_GAS + coupled.cumulative_oxygen_dose().value(),
        restarted.state.gas_oxygen.value()
            + restarted.state.solution.dissolved_oxygen.value()
            + 0.25 * restarted.state.solution.gold_dissolved.value(),
    );
    println!(
        "Only oxygen partial pressure is modeled; no air composition, unmetered gas reservoir, heat balance, cyanide-complex competition or carbon adsorption is inferred."
    );
    Ok(())
}
