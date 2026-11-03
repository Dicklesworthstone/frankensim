//! Solve finite-inventory goethite precipitation or dissolution at 298.15 K.
//!
//! Compile the four real source packs with `cargo run -p xtask -- matdb-pack
//! --manifest data/matdb/seed-v1/PROFILE/manifest.tsv --out PROFILE.fsmatpk`.
//! The required profile names and order are printed by `--help`.
//!
//! Example after compiling those packs:
//! ```text
//! cargo run -p fs-thermochem --example mineral_equilibrium -- \
//!   iron-iii-ion-aqueous-usgs-b2131.fsmatpk \
//!   water-liquid-usgs-b2131.fsmatpk \
//!   goethite-solid-usgs-b2131-thermochemistry.fsmatpk \
//!   hydrogen-ion-aqueous-usgs-b2131.fsmatpk \
//!   0.001 55.5 0 0.001 1
//! ```
//! The last five arguments are initial Fe3+, water, goethite and H+ amounts
//! [mol], followed by the aqueous amount [mol] at unit relative activity.
//! All amounts and that fixed solvent basis are caller inputs, not a plant dose.
//! Aqueous activity coefficients are one; the pure liquid and solid have unit
//! activity while their stocks remain finite. This one-reaction ideal closure
//! omits other aqueous species, electroneutrality closure, activity corrections,
//! precipitate competition and kinetics. It does not calculate full speciation.
#![deny(unsafe_code)]

use std::error::Error;

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            equilibrium::{EquilibriumInventory, EquilibriumSolveOptions, IdealActivityModel},
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const SLUGS: [&str; 4] = [
    "iron-iii-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
    "goethite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 4] = ["Fe[3+](aq)", "H2O(l)", "FeO(OH)(s)", "H+(aq)"];
const USAGE: &str = "usage: mineral_equilibrium FE3_PACK WATER_PACK GOETHITE_PACK HPLUS_PACK FE3_MOL WATER_MOL GOETHITE_MOL HPLUS_MOL AQUEOUS_UNIT_ACTIVITY_MOL";

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
        1 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        2 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-goethite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if arguments.is_empty() || arguments.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("Required source profiles, in order:");
        for slug in SLUGS {
            println!("  {slug}");
        }
        println!(
            "Initial amounts are finite nonnegative moles, including zero for an absent phase."
        );
        println!(
            "AQUEOUS_UNIT_ACTIVITY_MOL must be positive and defines a = n/n_ref for both ions."
        );
        println!(
            "The fixed solvent basis, gamma=1 and unit pure-phase activities are explicit ideal approximations."
        );
        println!(
            "Only Fe3+ + 2 H2O <=> goethite + 3 H+ is equilibrated; no full speciation or kinetic model."
        );
        return Ok(());
    }
    if arguments.len() != 9 {
        return Err(USAGE.into());
    }
    let mut packs = Vec::with_capacity(4);
    for (path, expected) in arguments[..4].iter().zip(SLUGS) {
        let bytes = std::fs::read(path)?;
        let decoded = NormalizedPack::from_bytes(&bytes)?;
        let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
        if pack.pack_id() != expected {
            return Err(format!(
                "expected source profile {expected}, received {}",
                pack.pack_id()
            )
            .into());
        }
        packs.push(pack);
    }
    let points = (0..4).map(point).collect::<Result<Vec<_>, _>>()?;
    let ids = SPECIES
        .iter()
        .map(|name| SpeciesId::new(*name))
        .collect::<Result<Vec<_>, _>>()?;
    let amounts = arguments[4..8]
        .iter()
        .map(|amount| amount.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let scale = Amount::new(arguments[8].parse::<f64>()?);
    let sources = std::array::from_fn::<_, 4, _>(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let inventory = std::array::from_fn::<_, 4, _>(|i| EquilibriumInventory {
        species: ids[i].clone(),
        amount: Amount::new(amounts[i]),
        activity: if i == 0 || i == 3 {
            IdealActivityModel::IdealDilute {
                amount_at_unit_activity: scale,
            }
        } else {
            IdealActivityModel::PureCondensed
        },
    });
    let options = EquilibriumSolveOptions {
        extent_tolerance: Amount::new(1e-12),
        log_residual_tolerance: Dimensionless::new(1e-9),
        max_iterations: 200,
    };
    let result = MineralReaction::FerricHydrolysisGoethite
        .build()?
        .equilibrate_ideal(&sources, &inventory, options)?;
    println!("Fe3+(aq) + 2 H2O(l) <=> FeO(OH)(goethite) + 3 H+(aq)");
    println!(
        "Reference state: {} K, {} Pa",
        result.standard.temperature().value(),
        result.standard.reference_pressure().value()
    );
    println!(
        "Declared aqueous model: fixed n_ref = {} mol, a = n/n_ref, gamma = 1",
        scale.value()
    );
    println!("Declared pure-phase model: water and goethite activity = 1, finite stocks.");
    println!(
        "Standard reaction Gibbs energy [J/mol]: {}",
        result.standard.molar().value()
    );
    println!("Status: {:?}", result.status);
    println!("Signed extent [mol]: {}", result.extent.value());
    println!(
        "Numerical extent bracket [mol]: [{}, {}]",
        result.extent_bracket[0].value(),
        result.extent_bracket[1].value()
    );
    println!(
        "Numerical tolerances: {} mol extent, {} ln(Q/K); iterations: {} / {}",
        options.extent_tolerance.value(),
        options.log_residual_tolerance.value(),
        result.iterations,
        options.max_iterations
    );
    for (initial, (species, final_amount)) in result.initial.iter().zip(&result.final_inventory) {
        println!(
            "{} [mol]: initial {}, final {}",
            species.as_str(),
            initial.amount.value(),
            final_amount.value()
        );
    }
    match result.log_quotient_over_equilibrium {
        Some(residual) => println!("ln(Q/K) [dimensionless]: {}", residual.value()),
        None => println!("ln(Q/K): undefined for the inventory-blocked state"),
    }
    println!("At a pure-phase inventory limit, the one-sided residual need not be zero.");
    println!(
        "This is one-reaction ideal equilibrium with a fixed solvent basis, not full aqueous speciation."
    );
    for term in result.standard.terms() {
        println!(
            "Formation source for {}: {}",
            term.term.species.as_str(),
            term.pack_id
        );
    }
    Ok(())
}
