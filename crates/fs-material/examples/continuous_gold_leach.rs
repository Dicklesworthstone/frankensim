//! Actual NBS acid/base data with a separately authored local cyanidation rate.
//! Demonstrates a throughput change, flow stop, feed change and restart in one
//! persistent mixed slurry tank. Dissolved gold leaving the tank is an exported
//! solution inventory; it is not inferred carbon recovery or refined metal.
#![deny(unsafe_code)]

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::cyanide_leach::{
    self, CyanidationCharge,
    flow::{
        CyanidationCstr, CyanidationFeed, CyanidationFlowAdvance, FlowNumerics,
        LOCAL_CYANIDATION_LAW, LocalGoldCyanidation, SpectatorFeed,
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, QuantitySpec, Temperature, Time, Volume,
    VolumetricFlowRate,
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
        },
    },
};
use std::{error::Error, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-cyanide";
const T: f64 = 298.15;
const V: f64 = 0.01;
const PACKS: [(&str, &str); 5] = [
    ("hydrogen-cyanide-aqueous-nbs1989", "HCN(aq)"),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)"),
    ("cyanide-ion-aqueous-nbs1989", "CN[-](aq)"),
    ("hydroxide-ion-aqueous-nbs1989", "OH[-](aq)"),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)"),
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
    if index == 4 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if index != 0 {
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
fn load(directory: &Path) -> Result<IdealAcidBaseSystem> {
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
    for standard in [water.standard(), cyanide.standard()] {
        for term in standard.terms() {
            println!(
                "equilibrium_source={} selection_receipt={}",
                term.pack_id,
                term.answer.receipt.content_hash()
            );
        }
    }
    Ok(IdealAcidBaseSystem::new(
        water,
        vec![AcidFamily::new(FAMILY, vec![cyanide])?],
    )?)
}
fn local_source() -> ConstitutiveModelCard {
    use cyanide_leach::*;
    let artifact = hash_domain(
        "continuous-gold-leach-example-v1",
        b"authored local rate 0.005 m6/(mol2 s); no measured ore kinetics or plant calibration",
    );
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, lo, hi) in [
        (PH_AXIS, Dims::NONE, 4.0, 14.0),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS, 0.0, 100.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 1.0),
        (GOLD_AXIS, Amount::DIMS, 0.0, 1.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, T, T),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 1.0),
        (TOTAL_CYANIDE_AXIS, Amount::DIMS, 0.0, 1.0),
        (GOLD_COMPLEX_AXIS, Amount::DIMS, 0.0, 1.0),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 1.0, 100.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    ConstitutiveModelCard { law: LawId(LOCAL_CYANIDATION_LAW.into()), law_version: 1,
        state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(RATE_PARAMETER.into(), LawParameter { value: 0.005, dims: CyanidationRateConstant::DIMS })].into(),
        validity, sources: vec![artifact], provenance: Provenance {
            source: "authored local kinetics and support for a numerical illustration; no measured ore calibration".into(),
            license: "project example".into(), artifact: Some(artifact),
        } }
}
fn report(label: &str, result: &CyanidationFlowAdvance) {
    println!(
        "{label}: time_s={:.3} flow_m3_per_s={:.6e} pH={:.6} free_CN_mol={:.9e}",
        result.elapsed.value(),
        result.flow.value(),
        result.state.speciation.ph.value(),
        result.state.free_cyanide.value()
    );
    println!(
        "  Au_solid_in_tank_mol={:.9e} Au_complex_in_tank_mol={:.9e} Au_reacted_this_interval_mol={:.9e}",
        result.state.inventory.accessible_gold.value(),
        result.state.inventory.gold_complex.value(),
        result.interval_gold_dissolved.value()
    );
    println!(
        "  interval_solid_Au_export_mol={:.9e} interval_dissolved_Au_export_mol={:.9e} cumulative_dissolved_Au_export_mol={:.9e}",
        result.outlet.accessible_gold.value(),
        result.outlet.gold_complex.value(),
        result.cumulative_outlet.gold_complex.value()
    );
    let balance = &result.cumulative_conservation;
    println!(
        "  cumulative_residuals_mol: Au={:.3e} CN={:.3e} O2={:.3e} Na={:.3e}; local_error_max_mol={:.3e}",
        balance.gold.value(),
        balance.cyanide.value(),
        balance.dissolved_oxygen.value(),
        balance.sodium.value(),
        result.maximum_local_error.value()
    );
    for residual in [
        balance.gold,
        balance.cyanide,
        balance.dissolved_oxygen,
        balance.sodium,
    ] {
        assert!(
            residual.value().abs() <= 1e-9,
            "independent material balance failed"
        );
    }
}
fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("usage: continuous_gold_leach --materials PACK_DIR");
        println!("Compile each profile to PACK_DIR/PACK_ID.fsmatpk with xtask matdb-pack:");
        for (name, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "Uses original NBS acid/base data and separately authored local kinetics; all flows, compositions and support boxes are numerical examples."
        );
        return Ok(());
    }
    if args.len() != 2 || args[0] != "--materials" {
        return Err("usage: continuous_gold_leach --materials PACK_DIR".into());
    }
    let chemistry = load(Path::new(&args[1]))?;
    let model = LocalGoldCyanidation::from_card(chemistry, FAMILY, local_source())?;
    println!(
        "LOCAL KINETICS: k={} m6/(mol2 s), law={}, independently authored; no ore or throughput calibration",
        model.rate_constant().value(),
        model.source().law.0
    );
    let sodium = SpeciesId::new("Na[+](aq)")?;
    let original = CyanidationCharge {
        accessible_gold: Amount::new(0.001),
        uncomplexed_cyanide: Amount::new(0.1),
        dissolved_oxygen: Amount::new(0.001),
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
        dissolved_oxygen: AmountConcentration::new(0.1),
        gold_complex: AmountConcentration::new(0.0),
        spectators: vec![SpectatorFeed {
            species: sodium,
            charge: 1,
            concentration: AmountConcentration::new(9.5),
        }],
        other_families: vec![],
        competing_metals: vec![],
    };
    let numerics = FlowNumerics {
        maximum_step: Time::new(5.0),
        absolute_amount_error: Amount::new(1e-10),
        relative_amount_error: Dimensionless::new(1e-6),
        maximum_work: 10_000_000,
        ..FlowNumerics::default()
    };
    let mut tank = CyanidationCstr::new(
        model,
        original.clone(),
        QueryPoint::new(),
        SpeciationOptions::default(),
    )?;
    println!(
        "SYNTHETIC operation: fixed V={V} m3, common ideal slurry residence; O2 is finite liquid/feed only. Exported dissolved Au is not carbon recovery."
    );
    for (label, seconds, q) in [
        ("initial flow", 120.0, 0.00005),
        ("higher throughput", 120.0, 0.0001),
    ] {
        report(
            label,
            &tank.advance(
                Time::new(seconds),
                VolumetricFlowRate::new(q),
                &inlet,
                numerics,
                || false,
            )?,
        );
    }
    let export_before = tank.cumulative_outlet();
    let progress_before = tank.cumulative_gold_dissolved();
    let stopped = tank.advance(
        Time::new(30.0),
        VolumetricFlowRate::new(0.0),
        &inlet,
        numerics,
        || false,
    )?;
    assert_eq!(tank.cumulative_outlet(), export_before);
    assert!(tank.cumulative_gold_dissolved().value() > progress_before.value());
    report("flow stopped; chemistry continues", &stopped);
    inlet.cyanide_ligand = AmountConcentration::new(12.0);
    inlet.spectators[0].concentration = AmountConcentration::new(11.5);
    report(
        "restart with changed dissolved NaCN feed",
        &tank.advance(
            Time::new(120.0),
            VolumetricFlowRate::new(0.000025),
            &inlet,
            numerics,
            || false,
        )?,
    );
    assert_eq!(tank.initial_charge(), &original);
    Ok(())
}
