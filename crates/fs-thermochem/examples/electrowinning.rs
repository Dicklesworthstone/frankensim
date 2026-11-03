//! Finite-inventory copper electrowinning for an explicitly authored scenario.
//!
//! Run with `--help` for constant-current and passed-charge inputs. All inputs
//! are supplied operating assumptions, not a measured plant calibration.
//! Copper's 0.063546 kg/mol matches the existing copper-solid-janaf-mining pack
//! and NIST WebBook: <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440508>.
//! No source-pack selection or evidence receipt is implied by this example.

use std::error::Error;

use fs_qty::{Amount, Current, Dimensionless, ElectricCharge, MolarMass, Time};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{electrolysis::ElectrolysisReaction, minerals::MineralReaction},
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!(
            "Usage: electrowinning [CURRENT_A SECONDS EFFICIENCY COPPER_ION_MOL WATER_MOL]\n\
                   electrowinning --charge COULOMBS EFFICIENCY COPPER_ION_MOL WATER_MOL\n\n\
Defaults: 1000 A, 3600 s, efficiency 0.95, 1000 mol Cu(II), 10000 mol water.\n\
All initial product amounts are zero. Efficiency is a fraction in [0,1].\n\
Reports ONLY the productive copper/water/oxygen/proton pathway. Unallocated\n\
charge has unmodeled side chemistry; this is not a complete operating cell balance."
        );
        return Ok(());
    }
    let passed_charge = args.first().is_some_and(|arg| arg == "--charge");
    if !args.is_empty() && args.len() != 5 {
        return Err("expected five arguments or no arguments; use --help".into());
    }
    let number = |index: usize, default: f64| -> Result<f64, Box<dyn Error>> {
        args.get(index).map_or(Ok(default), |arg| Ok(arg.parse()?))
    };
    let efficiency = Dimensionless::new(number(2, 0.95)?);
    let inventory = [
        ("Cu[2+](aq)", number(3, 1000.0)?),
        ("H2O(l)", number(4, 10000.0)?),
        ("Cu(s)", 0.0),
        ("O2(g)", 0.0),
        ("H+(aq)", 0.0),
    ]
    .into_iter()
    .map(|(name, amount)| Ok((SpeciesId::new(name)?, Amount::new(amount))))
    .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let model = ElectrolysisReaction::new(
        MineralReaction::CopperElectrowinning.build()?,
        SpeciesId::new("Cu(s)")?,
        MolarMass::new(0.063546),
        2,
    )?;
    let result = if passed_charge {
        model.advance_charge(&inventory, ElectricCharge::new(number(1, 0.0)?), efficiency)?
    } else {
        model.advance_current(
            &inventory,
            Current::new(number(0, 1000.0)?),
            Time::new(number(1, 3600.0)?),
            efficiency,
        )?
    };
    println!("Authored copper EW scenario; caller-supplied current efficiency.");
    println!("2 Cu[2+](aq) + 2 H2O(l) -> 2 Cu(s) + O2(g) + 4 H+(aq)");
    println!(
        "Electrons per written equation: {}",
        model.electrons_per_equation()
    );
    println!("Passed charge [C]: {}", result.passed_charge.value());
    println!(
        "Productive charge [C]: {}",
        result.productive_charge.value()
    );
    println!(
        "Unallocated charge [C]: {}",
        result.unallocated_charge.value()
    );
    println!("Written-equation extent [mol]: {}", result.extent.value());
    println!(
        "Copper deposited [mol]: {}",
        result.product_yield.amount().value()
    );
    println!(
        "Copper deposited [kg]: {}",
        result.product_yield.mass().value()
    );
    println!("species,initial_mol,modeled_final_mol");
    for (species, after) in &result.modeled_final_inventory {
        let before = inventory
            .iter()
            .find(|(name, _)| name == species)
            .ok_or("missing input species")?
            .1;
        println!("{},{},{}", species.as_str(), before.value(), after.value());
    }
    println!(
        "Inventory changes cover only the productive reaction; unallocated-charge side reactions, cell voltage, energy and kinetics remain unspecified."
    );
    Ok(())
}
