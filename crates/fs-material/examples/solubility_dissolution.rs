//! Finite-rate dissolution, concentration-induced crystallization, and heat.
//!
//! Run: cargo run -p fs-material --example solubility_dissolution
//! All coefficients below are authored demonstration values for one neutral
//! solute, not measured data or a calibration for a named real substance.

use fs_material::osmosis::solubility::{
    FiniteSolubilitySolution, kinetics::KineticSolubilitySolution,
};
use fs_qty::{Amount, AmountConcentration, Temperature, Time, Volume, VolumetricFlowRate};
use fs_thermochem::{
    MolarThermalQuantityV1, SpeciesId,
    solution::{
        caloric::{MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1},
        solubility::{
            IdealSolubilityV1,
            caloric::{
                SolubilityCaloricComponentV1, SolubilityCaloricModelV1, SolubilityCaloricSolveV1,
            },
            kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
        },
    },
};

fn print_state(
    label: &str,
    bath: &KineticSolubilitySolution,
) -> Result<(), Box<dyn std::error::Error>> {
    let phase = bath.last_steps()[0].state();
    println!(
        "{label}: dissolved={:.9} mol, solid={:.9} mol, equilibrium dissolved={:.9} mol, activity={:.9}",
        phase.dissolved().value(),
        phase.solid().value(),
        bath.last_steps()[0].equilibrium().dissolved.value(),
        bath.solution().water_activity()?,
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let solvent = SpeciesId::new("demonstration-solvent")?;
    let solute = SpeciesId::new("demonstration-neutral-solute")?;
    let temperature = Temperature::new(300.0);
    let molar_volume = 18e-6;
    let volume = Volume::new(1e-3);
    let knots = |capacity| {
        [280.0, 350.0]
            .into_iter()
            .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(capacity)))
            .collect()
    };
    let law = IdealSolubilityV1::new(
        solvent.clone(),
        solute.clone(),
        temperature,
        AmountConcentration::new(100.0),
        6000.0,
        [Temperature::new(280.0), Temperature::new(350.0)],
    )?
    .with_phase_heat_capacities(
        MolarHeatCapacityCurveV1::new(solute.clone(), temperature, knots(40.0))?,
        MolarHeatCapacityCurveV1::new(solute, temperature, knots(90.0))?,
    )?;

    // A well-mixed isothermal bath has a 100 s relaxation time (V/kA).
    // The existing solid supplies the growth surface; no nuclei are invented.
    let mut kinetic = KineticSolubilitySolution::new(
        temperature,
        molar_volume,
        volume,
        vec![(
            SolubilityKineticsV1::new(
                law.clone(),
                VolumetricFlowRate::new(1e-5),
                PrecipitationPolicyV1::ExistingSolidOnly,
            )?,
            SolubilityKineticStateV1::new(Amount::new(0.2), Amount::new(0.0))?,
        )],
    )?;
    print_state("Initial finite-rate bath", &kinetic)?;
    for seconds in [10.0, 20.0, 30.0, 120.0] {
        kinetic.advance(Time::new(seconds), || false)?;
        print_state(&format!("After another {seconds:.0} s"), &kinetic)?;
    }
    // Prescribed solvent removal retains both actual phase stocks immediately.
    // Crystallization then takes time under the same explicit fixed-area law.
    let mut concentrated = kinetic.with_solvent_volume_m3(2.5e-4)?;
    print_state("Immediately after solvent removal", &concentrated)?;
    concentrated.advance(Time::new(60.0), || false)?;
    print_state("After 60 s of crystal growth", &concentrated)?;

    // Separate equilibrium bath: a closed heat input uses the SAME distinct
    // phase Cp curves for its saturation and dissolution enthalpy. This is not
    // a claim of simultaneous finite-rate thermal coupling.
    let thermal = SolubilityCaloricModelV1::new(
        molar_volume,
        PureLiquidHeatCapacityV1::new(solvent, temperature, knots(75.0))?,
        vec![SolubilityCaloricComponentV1::from_phase_heat_capacities(
            law.clone(),
        )?],
    )?;
    let equilibrium = FiniteSolubilitySolution::new(
        temperature,
        molar_volume,
        volume,
        vec![(law, Amount::new(0.2))],
    )?;
    let heated = equilibrium.add_heat(
        &thermal,
        20_000.0,
        SolubilityCaloricSolveV1 {
            temperature_bracket: [Temperature::new(300.0), Temperature::new(330.0)],
            temperature_tolerance_k: 1e-9,
            energy_tolerance_j: 1e-5,
            maximum_iterations: 128,
        },
        || false,
    )?;
    println!(
        "Separate equilibrium bath after 20000 J: T={:.9} K, dissolved={:.9} mol, energy residual={:.3e} J",
        heated.solution().solution().temperature_k(),
        heated.solution().phase_partitions()[0].dissolved.value(),
        heated.balance().energy_residual_j(),
    );
    Ok(())
}
