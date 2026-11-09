//! Fixed-area versus shrinking-particle dissolution, crystal growth and coupled heat.
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
                kinetics::{SolubilityThermalKineticControlsV1, SolubilityThermalKineticsV1},
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
        "{label}: dissolved={:.9} mol, solid={:.9} mol, equilibrium dissolved={:.9} mol, activity={:.9}, kA={:.9e} m3/s",
        phase.dissolved().value(),
        phase.solid().value(),
        bath.last_steps()[0].equilibrium().dissolved.value(),
        bath.solution().water_activity()?,
        bath.models()[0].conductance_at(phase.solid())?.value(),
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
    // The same film coefficient and initial area, now carried by a fixed number
    // of self-similar crystals. Their kA changes as (solid/reference solid)^(2/3).
    // Solvent-volume changes preserve that original reference stock.
    let mut particles = KineticSolubilitySolution::new(
        temperature,
        molar_volume,
        volume,
        vec![(
            kinetic.models()[0]
                .clone()
                .with_particle_surface(Amount::new(0.2))?,
            *kinetic.last_steps()[0].state(),
        )],
    )?;
    print_state("Initial finite-rate bath", &kinetic)?;
    for seconds in [10.0, 20.0, 30.0, 120.0] {
        kinetic.advance(Time::new(seconds), || false)?;
        particles.advance(Time::new(seconds), || false)?;
        print_state(
            &format!("Fixed area after another {seconds:.0} s"),
            &kinetic,
        )?;
        print_state(
            &format!("Shrinking particles after another {seconds:.0} s"),
            &particles,
        )?;
    }
    // Prescribed solvent removal retains both actual phase stocks immediately.
    // The same crystals then grow, increasing area under the retained model.
    let mut concentrated = particles.with_solvent_volume_m3(2.5e-4)?;
    print_state("Immediately after solvent removal", &concentrated)?;
    concentrated.advance(Time::new(60.0), || false)?;
    print_state("After 60 s of crystal growth", &concentrated)?;

    // The SAME distinct phase Cp curves control saturation and the dissolution
    // enthalpy. Finite-rate adiabatic dissolution now cools the actual bath;
    // the temperature feedback and evolving phases share one energy account.
    let thermal = SolubilityCaloricModelV1::new(
        molar_volume,
        PureLiquidHeatCapacityV1::new(solvent, temperature, knots(75.0))?,
        vec![SolubilityCaloricComponentV1::from_phase_heat_capacities(
            law.clone(),
        )?],
    )?;
    let coupled = SolubilityThermalKineticsV1::new(thermal.clone(), particles.models().to_vec())?;
    let adiabatic = particles.advance_with_heat(
        &coupled,
        Time::new(60.0),
        0.0,
        SolubilityThermalKineticControlsV1 {
            root: SolubilityCaloricSolveV1 {
                temperature_bracket: [Temperature::new(290.0), Temperature::new(310.0)],
                temperature_tolerance_k: 1e-9,
                energy_tolerance_j: 1e-5,
                maximum_iterations: 128,
            },
            initial_step: Time::new(10.0),
            minimum_step: Time::new(1e-4),
            maximum_step: Time::new(10.0),
            relative_tolerance: 1e-4,
            amount_absolute_tolerance_mol: 1e-10,
            temperature_absolute_tolerance_k: 1e-7,
            maximum_attempts: 4096,
        },
        || false,
    )?;
    print_state(
        "After 60 s of adiabatic finite-rate dissolution",
        adiabatic.solution(),
    )?;
    println!(
        "Coupled bath: T={:.9} K, energy residual={:.3e} J, accepted steps={}",
        adiabatic.solution().solution().temperature_k(),
        adiabatic.balance().energy_residual_j(),
        adiabatic.balance().accepted_steps(),
    );

    // An independent equilibrium bath demonstrates the infinite-rate limit
    // under a closed heat input, using the same constitutive sources.
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
