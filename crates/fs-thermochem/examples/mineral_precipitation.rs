//! Simultaneous iron/goethite partition and aqueous pH from original USGS data.
//!
//! Compile the five profiles listed by --help into PACK_DIR, then run:
//! ```text
//! cargo run -p fs-thermochem --example mineral_precipitation -- \
//!   --materials PACK_DIR 0.001 0.003 0 10
//! ```
//! The four quantities are total Fe component, fixed chloride, fixed sodium
//! and amount at unit aqueous activity, all in mol. The example's explicit
//! inputs resemble a dilute ferric-chloride analytical inventory; they are
//! authored model scenarios, not a reagent recipe or measured plant solution.
//! Raising chloride at fixed Fe models separately declared acid equivalents;
//! raising sodium models separately declared base equivalents. The default mode
//! reports equilibrium. Commercial reagent purity and dosing are not calculated.
//!
//! Add `--kinetics DISSOLVED_FE_MOL GOETHITE_MOL CHLORIDE_MOL SODIUM_MOL
//! BASIS_MOL SURFACE_RATE_MOL_S SECONDS` after `--materials PACK_DIR` to
//! evolve explicit actual phases with fast aqueous speciation and finite
//! mineral exchange. The supplied surface rate is an independent constitutive
//! input; the USGS equilibrium sources do not establish that kinetic rate.
//!
//! Only free Fe3+, water ions, fixed spectator ions and the selected goethite
//! phase participate. Solvent basis and water/solid unit activities stay fixed.
//! Other iron complexes, competing solids and nonideal activities are omitted.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                IdealAcidBaseSystem, SpectatorIon, WaterAutoprotolysis,
                precipitation::{
                    HydrolysisPrecipitation, PrecipitationOptions,
                    kinetics::{
                        MineralGrowthPolicyV1, MineralKineticControlsV1,
                        MineralKineticStateV1, MineralKineticsV1, MineralSurfaceRateV1,
                    },
                },
            },
        },
    },
};
use std::{error::Error, path::Path};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const PACKS: [&str; 5] = [
    "iron-iii-ion-aqueous-usgs-b2131",
    "goethite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
    "hydroxide-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
];
const SPECIES: [&str; 5] = ["Fe[3+](aq)", "FeO(OH)(s)", "H+(aq)", "OH[-](aq)", "H2O(l)"];
const USAGE: &str = "usage: mineral_precipitation --materials PACK_DIR TOTAL_FE_MOL CHLORIDE_MOL SODIUM_MOL AQUEOUS_UNIT_ACTIVITY_MOL";
const KINETIC_USAGE: &str = "       mineral_precipitation --materials PACK_DIR --kinetics DISSOLVED_FE_MOL GOETHITE_MOL CHLORIDE_MOL SODIUM_MOL AQUEOUS_UNIT_ACTIVITY_MOL SURFACE_RATE_MOL_S SECONDS";

fn point(index: usize) -> Result<QueryPoint, Box<dyn Error>> {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
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
        1 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-goethite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-individual-ion", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn load(directory: &Path) -> Result<HydrolysisPrecipitation, Box<dyn Error>> {
    let packs = PACKS
        .iter()
        .map(|name| {
            let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
            let decoded = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
            if pack.pack_id() != *name {
                return Err(format!(
                    "expected source profile {name}, received {}",
                    pack.pack_id()
                )
                .into());
            }
            Ok(pack)
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let points = (0..5).map(point).collect::<Result<Vec<_>, _>>()?;
    let ids = SPECIES
        .iter()
        .map(|s| SpeciesId::new(*s))
        .collect::<Result<Vec<_>, _>>()?;
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
        ids[2].clone(),
        ids[3].clone(),
        &sources(&[4, 2, 3]),
    )?;
    let aqueous = IdealAcidBaseSystem::new(water, vec![])?;
    Ok(HydrolysisPrecipitation::from_sources(
        aqueous,
        MineralReaction::FerricHydrolysisGoethite.build()?,
        ids[0].clone(),
        ids[1].clone(),
        &sources(&[0, 4, 1, 2]),
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("{KINETIC_USAGE}");
        println!("Compile these source profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for pack in PACKS {
            println!("  {pack}");
        }
        println!(
            "Use: cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!(
            "All stock arguments are finite nonnegative moles; the unit-activity amount is positive."
        );
        println!("TOTAL_FE_MOL is free Fe3+ plus the selected goethite on a one-Fe formula basis.");
        println!(
            "The fixed solvent basis gives a=n/n_ref, gamma=1; pure water and goethite have unit activity."
        );
        println!(
            "Chloride and sodium are explicit fixed spectator stocks, not inferred counterions or doses."
        );
        println!(
            "Kinetic mode preserves declared dissolved/solid stocks and solves aqueous pH at every phase trial."
        );
        println!(
            "Its independent fixed-surface rate is K_A*(1-IAP/K), in mol formula units/s; existing solid is required for growth."
        );
        return Ok(());
    }
    if args.len() == 10 && args[0] == "--materials" && args[2] == "--kinetics" {
        return run_kinetic(load(Path::new(&args[1]))?, &args[3..]);
    }
    if args.len() != 6 || args[0] != "--materials" {
        return Err(format!("{USAGE}\n{KINETIC_USAGE}").into());
    }
    let model = load(Path::new(&args[1]))?;
    let numbers = args[2..]
        .iter()
        .map(|s| s.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let total = Amount::new(numbers[0]);
    let spectators = [("Cl[-](aq)", -1, numbers[1]), ("Na[+](aq)", 1, numbers[2])]
        .into_iter()
        .map(|(name, charge, amount)| {
            Ok(SpectatorIon {
                species: SpeciesId::new(name)?,
                charge,
                amount: Amount::new(amount),
            })
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let basis = Amount::new(numbers[3]);
    let answer = model.speciate(
        total,
        &[],
        &spectators,
        basis,
        PrecipitationOptions::default(),
        || false,
    )?;
    println!("Coupled Fe3+(aq) + 2 H2O(l) <=> goethite + 3 H+(aq) and water autoprotolysis.");
    println!(
        "Original USGS reference: 298.15 K, 100000 Pa; fixed ideal aqueous basis {} mol.",
        basis.value()
    );
    println!("All aqueous activities use a=n/n_ref and gamma=1; water and goethite activity=1.");
    println!("Total Fe component [mol]: {}", total.value());
    println!(
        "Fixed chloride/sodium [mol]: {} / {}",
        numbers[1], numbers[2]
    );
    println!("Status: {:?}", answer.status);
    println!("pH=-log10(aH): {:.12}", answer.aqueous.ph.value());
    println!("Free Fe3+ [mol]: {:.16e}", answer.dissolved_amount.value());
    println!("Goethite [mol]: {:.16e}", answer.solid_amount.value());
    println!(
        "Free H+ / OH- [mol]: {:.12e} / {:.12e}",
        answer.aqueous.proton_amount.value(),
        answer.aqueous.hydroxide_amount.value()
    );
    println!(
        "Component balance residual [mol]: {:.3e}",
        answer.component_balance_residual.value()
    );
    println!(
        "Aqueous charge residual [mol equivalents]: {:.3e}",
        answer.aqueous.charge_balance_residual.value()
    );
    match answer.log_saturation_ratio {
        Some(value) => println!("Natural log of dissolution IAP/K: {:.12e}", value.value()),
        None => println!("Phase saturation is undefined for an empty Fe component."),
    }
    if let Some(spread) = answer.log_saturation_numerical_spread {
        println!(
            "Remaining log-saturation variation over the inner pH bracket: {:.12e}",
            spread.value()
        );
    }
    println!(
        "Numerical dissolved bracket [mol]: [{:.16e}, {:.16e}]",
        answer.dissolved_bracket[0].value(),
        answer.dissolved_bracket[1].value()
    );
    println!(
        "Outer bisections / aqueous solves / inner pH bisections: {} / {} / {}",
        answer.iterations, answer.aqueous_solves, answer.aqueous_iterations
    );
    println!(
        "An absent phase retains its undersaturation inequality; no mineral mass is invented."
    );
    println!(
        "No other Fe complexes, competing minerals, nonideal activities, kinetics or finite-solvent H/O balance are included."
    );
    for pack in PACKS {
        println!("Original material source: {pack}");
    }
    Ok(())
}

fn kinetic_row(seconds: f64, state: &MineralKineticStateV1) {
    let saturation = state
        .log_saturation_ratio()
        .map_or_else(|| "zero IAP".to_owned(), |value| format!("{:.8e}", value.value()));
    println!(
        "{seconds:.8e}\t{:.12e}\t{:.12e}\t{:.10}\t{saturation}\t{:.3e}\t{:.3e}",
        state.dissolved_cation().value(),
        state.solid_amount().value(),
        state.aqueous().ph.value(),
        state.cation_balance_residual().value(),
        state.aqueous().charge_balance_residual.value(),
    );
}

fn run_kinetic(model: HydrolysisPrecipitation, args: &[String]) -> Result<(), Box<dyn Error>> {
    let values = args
        .iter()
        .map(|value| value.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let seconds = values[6];
    if !seconds.is_finite() || seconds <= 0.0 || seconds / 50.0 == 0.0 {
        return Err("kinetic observation duration must be finite, positive and representable".into());
    }
    let spectators = [("Cl[-](aq)", -1, values[2]), ("Na[+](aq)", 1, values[3])]
        .into_iter()
        .map(|(name, charge, amount)| {
            Ok(SpectatorIon {
                species: SpeciesId::new(name)?,
                charge,
                amount: Amount::new(amount),
            })
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let kinetics = MineralKineticsV1::from_hydrolysis(
        model.clone(),
        MineralSurfaceRateV1::new(values[5], MineralGrowthPolicyV1::ExistingSolidOnly)?,
    )?;
    let options = PrecipitationOptions::default().aqueous;
    let mut state = kinetics.prepare(
        Amount::new(values[0]),
        Amount::new(values[1]),
        vec![],
        spectators,
        Amount::new(values[4]),
        options,
        || false,
    )?;
    let controls = MineralKineticControlsV1 {
        aqueous: options,
        initial_step: Time::new(seconds / 50.0),
        minimum_step: Time::new(seconds * 1e-10),
        maximum_step: Time::new(seconds / 10.0),
        relative_tolerance: Dimensionless::new(1e-4),
        amount_absolute_tolerance: Amount::new(1e-12),
        extent_absolute_tolerance: Amount::new(1e-14),
        ..MineralKineticControlsV1::default()
    };
    println!("Finite-rate Fe3+/goethite exchange with the original USGS equilibrium sources.");
    println!("Fixed 298.15 K and 100000 Pa; ideal aqueous a=n/n_ref; unit water/solid activity.");
    println!(
        "Explicit independent surface rate: {:.8e} mol/s; existing-solid growth; no nucleation.",
        values[5]
    );
    println!("time_s\tdissolved_Fe_mol\tgoethite_mol\tpH\tln_IAP_over_K\tFe_residual_mol\tcharge_residual_mol");
    kinetic_row(0.0, &state);
    let mut accepted = 0;
    let mut rejected = 0;
    let mut aqueous_solves = 0;
    let mut blocked = false;
    let mut exhausted = false;
    for window in 1..=10 {
        let result = kinetics.advance(&state, Time::new(seconds / 10.0), controls, || false)?;
        accepted += result.accepted_steps();
        rejected += result.rejected_steps();
        aqueous_solves += result.aqueous_solves();
        blocked |= result.precipitation_blocked();
        exhausted |= result.solid_exhausted();
        state = result.into_state();
        kinetic_row(seconds * (window as f64 / 10.0), &state);
    }
    // Explicitly compare the accepted finite-rate state with the distinct
    // infinite-rate limit; this diagnostic never feeds back into the trajectory.
    let equilibrium = model.speciate(
        state.component_total(),
        state.family_totals(),
        state.spectators(),
        state.amount_at_unit_activity(),
        PrecipitationOptions::default(),
        || false,
    )?;
    println!(
        "Equilibrium reference: free Fe={:.12e} mol; goethite={:.12e} mol; pH={:.10}",
        equilibrium.dissolved_amount.value(),
        equilibrium.solid_amount.value(),
        equilibrium.aqueous.ph.value(),
    );
    println!(
        "Accepted/rejected intervals: {accepted}/{rejected}; aqueous solves: {aqueous_solves}; growth blocked: {blocked}; solid exhausted: {exhausted}"
    );
    println!("The implicit endpoint does not supply an exact continuous solid-exhaustion time.");
    println!("No other Fe complexes, competing minerals, nonideal activities or reaction heat are included.");
    Ok(())
}
