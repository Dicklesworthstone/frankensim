//! Actual NBS aqueous/gas equilibrium with explicitly authored local leach and
//! transfer coefficients. Liquid feed contains no dissolved oxygen: the maintained
//! gas reservoir makes oxygen available while slurry and gold leave the tank.
#![deny(unsafe_code)]

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::cyanide_leach::{
    self, CyanidationCharge,
    aeration::{
        KLA_PARAMETER, OXYGEN_PARTIAL_PRESSURE_AXIS, OXYGEN_TRANSFER_LAW, OxygenTransfer,
        SATURATION_CONCENTRATION_AXIS,
    },
    flow::{
        CyanidationCstr, CyanidationFeed, FlowNumerics, LOCAL_CYANIDATION_LAW,
        LocalGoldCyanidation, SpectatorFeed,
        aeration::{AeratedCyanidationCstr, AeratedFlowAdvance},
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Frequency, Pressure, QuantitySpec,
    Temperature, Time, Volume, VolumetricFlowRate,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
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
use std::{error::Error, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-cyanide";
const T: f64 = 298.15;
const V: f64 = 0.01;
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
            T,
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
        .map(|(name, _)| -> Result<NormalizedPack> {
            let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
            let decoded = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
            if pack.pack_id() != *name || pack.to_bytes() != bytes {
                return Err(format!("expected canonical source pack {name}").into());
            }
            println!(
                "material={} pack_hash={}",
                pack.pack_id(),
                pack.content_hash()
            );
            Ok(pack)
        })
        .collect::<Result<Vec<_>>>()?;
    let ids = PACKS
        .iter()
        .map(|(_, name)| SpeciesId::new(*name))
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
    let oxygen = IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build()?,
        ids[5].clone(),
        ids[6].clone(),
        &sources(&[5, 6]),
    )?;
    for standard in [water.standard(), cyanide.standard(), oxygen.standard()] {
        for term in standard.terms() {
            println!(
                "equilibrium_source={} selection_receipt={}",
                term.pack_id,
                term.answer.receipt.content_hash()
            );
        }
    }
    Ok((
        IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, vec![cyanide])?])?,
        oxygen,
    ))
}
fn card(
    law: &str,
    parameter: &str,
    value: f64,
    dims: Dims,
    validity: ValidityDomain,
    note: &str,
) -> ConstitutiveModelCard {
    let artifact = hash_domain("continuous-aerated-gold-leach-example-v1", note.as_bytes());
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(parameter.into(), LawParameter { value, dims })].into(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: note.into(),
            license: "project example".into(),
            artifact: Some(artifact),
        },
    }
}
fn local_source() -> ConstitutiveModelCard {
    use cyanide_leach::*;
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, lo, hi) in [
        (PH_AXIS, Dims::NONE, 4.0, 14.0),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS, 0.0, 100.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (GOLD_AXIS, Amount::DIMS, 0.0, 1.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, T, T),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 1.0),
        (TOTAL_CYANIDE_AXIS, Amount::DIMS, 0.0, 1.0),
        (GOLD_COMPLEX_AXIS, Amount::DIMS, 0.0, 1.0),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 1.0, 100.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    card(
        LOCAL_CYANIDATION_LAW,
        RATE_PARAMETER,
        0.005,
        CyanidationRateConstant::DIMS,
        validity,
        "authored k=0.005 m6/(mol2 s) and local support; no ore calibration",
    )
}
fn transfer_source() -> ConstitutiveModelCard {
    use cyanide_leach::*;
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, lo, hi) in [
        (OXYGEN_PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 0.0, 100_000.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 10_000.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, T, T),
        (VOLUME_AXIS, Volume::DIMS, V, V),
        (
            SATURATION_CONCENTRATION_AXIS,
            AmountConcentration::DIMS,
            0.0,
            10.0,
        ),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 10.0, 10.0),
        (PH_AXIS, Dims::NONE, 4.0, 14.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    card(
        OXYGEN_TRANSFER_LAW,
        KLA_PARAMETER,
        0.02,
        Frequency::DIMS,
        validity,
        "authored kLa=0.02 s^-1 and transfer support; no measured slurry/aerator calibration",
    )
}
fn report(label: &str, answer: &AeratedFlowAdvance) {
    let f = &answer.flow;
    println!(
        "{label}: elapsed_s={:.3} Q_m3_per_s={:.6e} pO2_Pa={:.1} pH={:.6}",
        f.elapsed.value(),
        f.flow.value(),
        answer.oxygen_partial_pressure.value(),
        f.state.speciation.ph.value()
    );
    println!(
        "  DO_mol={:.9e} saturation_mol={:.9e} free_CN_mol={:.9e}",
        f.state.inventory.dissolved_oxygen.value(),
        answer.oxygen_saturation.value(),
        f.state.free_cyanide.value()
    );
    println!(
        "  interval_gas_O2_in_mol={:.9e} interval_gas_O2_out_mol={:.9e} liquid_O2_export_mol={:.9e}",
        answer.oxygen_in.value(),
        answer.oxygen_out.value(),
        f.outlet.dissolved_oxygen.value()
    );
    println!(
        "  cumulative_gas_O2_in_mol={:.9e} cumulative_gas_O2_out_mol={:.9e}",
        answer.cumulative_oxygen_in.value(),
        answer.cumulative_oxygen_out.value()
    );
    println!(
        "  cumulative_Au_dissolved_mol={:.9e} retained_Au_complex_mol={:.9e} cumulative_Au_complex_export_mol={:.9e}",
        f.cumulative_gold_dissolved.value(),
        f.state.inventory.gold_complex.value(),
        f.cumulative_outlet.gold_complex.value()
    );
    let b = &f.cumulative_conservation;
    println!(
        "  balance_residuals_mol Au={:.3e} CN={:.3e} O2={:.3e} Na={:.3e}; local_error_max_mol={:.3e}",
        b.gold.value(),
        b.cyanide.value(),
        b.dissolved_oxygen.value(),
        b.sodium.value(),
        f.maximum_local_error.value()
    );
    for residual in [b.gold, b.cyanide, b.dissolved_oxygen, b.sodium] {
        assert!(
            residual.value().abs() < 1e-9,
            "independent material balance failed"
        );
    }
}
fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("usage: continuous_aerated_gold_leach --materials PACK_DIR");
        println!("Compile these original source profiles with xtask matdb-pack:");
        for (name, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "NBS equilibrium; synthetic local k, kLa, inventories, flows and pressure schedule."
        );
        return Ok(());
    }
    if args.len() != 2 || args[0] != "--materials" {
        return Err("usage: continuous_aerated_gold_leach --materials PACK_DIR".into());
    }
    let (chemistry, gas) = load(Path::new(&args[1]))?;
    let model = LocalGoldCyanidation::from_card(chemistry, FAMILY, local_source())?;
    let transfer = OxygenTransfer::from_card(gas, transfer_source())?;
    println!(
        "SYNTHETIC rates: local_k={} m6/(mol2 s), kLa={} s^-1",
        model.rate_constant().value(),
        transfer.coefficient().value()
    );
    println!(
        "Fixed solution V={V} m3, ideal solvent basis=10 mol; initial/feed DO=0. Maintained gas is not a finite oxygen-supply or power model."
    );
    let sodium = SpeciesId::new("Na[+](aq)")?;
    let original = CyanidationCharge {
        accessible_gold: Amount::new(0.001),
        uncomplexed_cyanide: Amount::new(0.1),
        dissolved_oxygen: Amount::new(0.0),
        gold_complex: Amount::new(0.0),
        amount_at_unit_activity: Amount::new(10.0),
        solution_volume: Volume::new(V),
        temperature: Temperature::new(T),
        spectators: vec![SpectatorIon {
            species: sodium.clone(),
            charge: 1,
            amount: Amount::new(0.095),
        }],
        other_families: vec![],
    };
    let mut inlet = CyanidationFeed {
        accessible_gold: AmountConcentration::new(0.1),
        cyanide_ligand: AmountConcentration::new(10.0),
        dissolved_oxygen: AmountConcentration::new(0.0),
        gold_complex: AmountConcentration::new(0.0),
        spectators: vec![SpectatorFeed {
            species: sodium,
            charge: 1,
            concentration: AmountConcentration::new(9.5),
        }],
        other_families: vec![],
        competing_metals: vec![],
    };
    let liquid = CyanidationCstr::new(
        model,
        original.clone(),
        QueryPoint::new(),
        SpeciationOptions::default(),
    )?;
    let mut tank = AeratedCyanidationCstr::from_tank(liquid, transfer)?;
    let numerics = FlowNumerics {
        maximum_step: Time::new(5.0),
        absolute_amount_error: Amount::new(1e-10),
        relative_amount_error: Dimensionless::new(1e-6),
        maximum_work: 10_000_000,
        ..FlowNumerics::default()
    };
    for (label, seconds, q, pressure) in [
        ("gas-supported initial flow", 120.0, 0.00005, 21_000.0),
        ("higher flow and enriched oxygen", 120.0, 0.0001, 80_000.0),
        ("flow stop and oxygen stripping", 30.0, 0.0, 0.0),
    ] {
        let answer = tank.advance(
            Time::new(seconds),
            VolumetricFlowRate::new(q),
            &inlet,
            Pressure::new(pressure),
            &QueryPoint::new(),
            numerics,
            || false,
        )?;
        if q == 0.0 {
            assert_eq!(answer.flow.outlet.dissolved_oxygen.value(), 0.0);
            assert_eq!(answer.flow.outlet.gold_complex.value(), 0.0);
            assert!(answer.oxygen_out.value() > 0.0);
        }
        assert!(answer.flow.cumulative_gold_dissolved.value() > 0.0);
        report(label, &answer);
    }
    inlet.cyanide_ligand = AmountConcentration::new(12.0);
    inlet.spectators[0].concentration = AmountConcentration::new(11.5);
    let answer = tank.advance(
        Time::new(120.0),
        VolumetricFlowRate::new(0.000025),
        &inlet,
        Pressure::new(21_000.0),
        &QueryPoint::new(),
        numerics,
        || false,
    )?;
    report("restart with changed cyanide feed", &answer);
    assert!(answer.cumulative_oxygen_in.value() > 0.0);
    assert!(answer.cumulative_oxygen_out.value() > 0.0);
    assert_eq!(tank.tank().initial_charge(), &original);
    Ok(())
}
