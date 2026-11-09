//! A neutral cosolute changes saturation, then dilution redissolves the crystals.
//!
//! Run: cargo run -p fs-material --example cosolute_crystallization
//! Add `-- --heat` for an explicitly temperature-qualified interaction family
//! with adiabatic crystal growth and prescribed heating/cooling.
//! Add `-- --jacket` for the same bath cooled/heated through a finite thermal
//! conductance to a prescribed jacket temperature. Heat then follows the
//! evolving bath temperature and is solved inside the original time integrator.
//! Add `-- --finite-jacket` to retain the coolant's finite sensible heat store:
//! coolant and solution temperatures evolve with equal/opposite transferred heat.
//! Values are authored demonstration coefficients, not measurements for a named
//! substance. The second species has zero phase-transfer conductance, so its
//! dissolved stock changes neither during crystallization nor during imposed
//! dilution. Its actual concentration still changes the first species' activity.
//! The default fixed-temperature model supplies no heat balance; the thermal
//! mode separately declares the interaction temperature law and reference Cp.

use std::sync::Arc;

use fs_material::osmosis::{
    nonideal::VirialSolutionModel,
    solubility::kinetics::{KineticSolubilitySolution, VirialSolubilityThermalModel},
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
                kinetics::{
                    SensibleHeatReservoirV1, SolubilityThermalKineticControlsV1,
                    SolubilityThermalKineticsV1,
                },
            },
            kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
        },
    },
};

fn initial() -> Result<KineticSolubilitySolution, Box<dyn std::error::Error>> {
    let temperature = Temperature::new(300.0);
    let mut components = Vec::new();
    for (name, dissolved, solid, conductance) in [
        ("crystallizing-neutral-solute", 0.1, 0.05, 1e-5),
        ("retained-neutral-cosolute", 0.05, 0.0, 0.0),
    ] {
        let law = IdealSolubilityV1::new(
            SpeciesId::new("demonstration-solvent")?,
            SpeciesId::new(name)?,
            temperature,
            AmountConcentration::new(100.0),
            6000.0,
            [Temperature::new(280.0), Temperature::new(330.0)],
        )?;
        let kinetics = SolubilityKineticsV1::new(
            law,
            VolumetricFlowRate::new(conductance),
            PrecipitationPolicyV1::ExistingSolidOnly,
        )?
        .with_particle_surface(Amount::new(0.05))?;
        components.push((
            kinetics,
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))?,
        ));
    }
    Ok(KineticSolubilitySolution::new(
        temperature,
        18e-6,
        Volume::new(1e-3),
        components,
    )?)
}

fn print_state(
    label: &str,
    solution: &KineticSolubilitySolution,
    phase_work_j: Option<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    let aqueous = solution.solution();
    let concentration: Vec<_> = aqueous
        .solute_moles()
        .iter()
        .map(|amount| amount / aqueous.solvent_volume_m3())
        .collect();
    let state = solution.last_steps()[0].state();
    let target = solution.conditional_saturation(0, &concentration)?;
    let work = phase_work_j.map_or_else(String::new, |value| format!("{value:.12e}"));
    println!(
        "{label},{:.9e},{:.9},{:.9},{:.12e},{:.12e},{:.9},{:.9},{:.12e},{work}",
        aqueous.solvent_volume_m3(),
        concentration[0],
        concentration[1],
        state.dissolved().value(),
        state.solid().value(),
        target.value(),
        aqueous.osmotic_pressure_pa()?,
        solution.models()[0].conductance_at(state.solid())?.value(),
    );
    Ok(())
}

fn advance(
    solution: &mut KineticSolubilitySolution,
    seconds: f64,
) -> Result<f64, Box<dyn std::error::Error>> {
    let before = solution.solution().solute_moles().to_vec();
    let volume = solution.solution().solvent_volume_m3();
    solution.advance(Time::new(seconds), || false)?;
    Ok(solution.solid_relative_inventory_free_energy_change_j(
        volume,
        &before,
        volume,
        solution.solution().solute_moles(),
    )?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().any(|argument| argument == "--finite-jacket") {
        return thermal_demo(ThermalMode::FiniteJacket);
    }
    if std::env::args().any(|argument| argument == "--jacket") {
        return thermal_demo(ThermalMode::Jacket);
    }
    if std::env::args().any(|argument| argument == "--heat") {
        return thermal_demo(ThermalMode::SuppliedHeat);
    }
    let mut ideal = initial()?;
    let original_virial = Arc::new(VirialSolutionModel::new(
        300.0,
        18e-6,
        1.0,
        vec![150.0, 150.0],
        vec![0.002, 0.004, 0.004, 0.002],
    )?);
    let mut interacting = initial()?.with_nonideal_thermodynamics(original_virial)?;
    println!(
        "stage,solvent_volume_m3,solute_mol_m3,cosolute_mol_m3,dissolved_mol,solid_mol,conditional_saturation_mol_m3,osmotic_pressure_Pa,particle_kA_m3_s,phase_work_J"
    );
    print_state("initial_ideal", &ideal, None)?;
    print_state("initial_interacting", &interacting, None)?;
    for seconds in [60.0, 120.0, 300.0] {
        let ideal_work = advance(&mut ideal, seconds)?;
        let interacting_work = advance(&mut interacting, seconds)?;
        print_state(
            &format!("ideal_after_another_{seconds}_s"),
            &ideal,
            Some(ideal_work),
        )?;
        print_state(
            &format!("interacting_after_another_{seconds}_s"),
            &interacting,
            Some(interacting_work),
        )?;
    }
    // External pure-solvent addition is prescribed, not claimed as a closed
    // water/mechanical-energy balance. Both component inventories are retained.
    interacting = interacting.with_solvent_volume_m3(2e-3)?;
    print_state("immediately_after_dilution", &interacting, None)?;
    let work = advance(&mut interacting, 1000.0)?;
    print_state("after_dilution_and_1000_s", &interacting, Some(work))?;
    Ok(())
}

enum ThermalDrive {
    SuppliedHeat(f64),
    Jacket {
        temperature: Temperature,
        conductance_w_k: f64,
    },
}

enum ThermalMode {
    SuppliedHeat,
    Jacket,
    FiniteJacket,
}

fn thermal_demo(mode: ThermalMode) -> Result<(), Box<dyn std::error::Error>> {
    let original = initial()?;
    let reference_t = Temperature::new(300.0);
    let knots = |capacity| {
        [280.0, 330.0]
            .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(capacity)))
            .to_vec()
    };
    // These are the explicitly declared REFERENCE phase capacities. The
    // interaction potential adds its consistent mixture/reference calorics;
    // they are not independent measured effective mixture heat capacities.
    let components = original
        .models()
        .iter()
        .map(|model| {
            SolubilityCaloricComponentV1::new(
                model.law().clone(),
                MolarHeatCapacityCurveV1::new(
                    model.law().solute().clone(),
                    reference_t,
                    knots(80.0),
                )?,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let caloric = SolubilityCaloricModelV1::new(
        18e-6,
        PureLiquidHeatCapacityV1::new(original.solvent().clone(), reference_t, knots(75.0))?,
        components,
    )?;
    let reference = Arc::new(VirialSolutionModel::new(
        300.0,
        18e-6,
        1.0,
        vec![150.0, 150.0],
        vec![0.002, 0.004, 0.004, 0.002],
    )?);
    let family = Arc::new(VirialSolubilityThermalModel::new(
        caloric.clone(),
        reference,
        // C [K m³/mol] in B(T)=B_ref+C*(1/T-1/T_ref).
        vec![0.2, 0.6, 0.6, 0.1],
        [Temperature::new(290.0), Temperature::new(315.0)],
    )?);
    let thermal = SolubilityThermalKineticsV1::new(caloric, original.models().to_vec())?;
    let mut bath = original.with_thermal_nonideal_thermodynamics(family)?;
    let controls = SolubilityThermalKineticControlsV1 {
        root: SolubilityCaloricSolveV1 {
            temperature_bracket: [Temperature::new(290.0), Temperature::new(315.0)],
            temperature_tolerance_k: 1e-9,
            energy_tolerance_j: 1e-6,
            maximum_iterations: 128,
        },
        initial_step: Time::new(5.0),
        minimum_step: Time::new(1e-5),
        maximum_step: Time::new(10.0),
        relative_tolerance: 1e-5,
        amount_absolute_tolerance_mol: 1e-11,
        temperature_absolute_tolerance_k: 1e-7,
        maximum_attempts: 8192,
    };
    if matches!(mode, ThermalMode::FiniteJacket) {
        return finite_thermal_demo(bath, &thermal, controls);
    }
    let mut elapsed = 0.0;
    let schedule = if matches!(mode, ThermalMode::Jacket) {
        // The prescribed reservoir represents a controlled jacket. Its heat
        // capacity, coolant flow and control-system power are not inferred.
        [
            ("jacket_cooling", 600.0, 292.0),
            ("cold_hold", 300.0, 292.0),
            ("jacket_heating", 900.0, 310.0),
            ("jacket_recooling", 900.0, 295.0),
        ]
        .map(|(label, seconds, temperature)| {
            (
                label,
                seconds,
                ThermalDrive::Jacket {
                    temperature: Temperature::new(temperature),
                    conductance_w_k: 25.0,
                },
            )
        })
    } else {
        [
            ("adiabatic_growth", 60.0, 0.0),
            ("adiabatic_relaxation", 120.0, 0.0),
            ("prescribed_heating", 120.0, 2500.0),
            ("prescribed_cooling", 240.0, -2500.0),
        ]
        .map(|(label, seconds, heat)| (label, seconds, ThermalDrive::SuppliedHeat(heat)))
    };
    println!(
        "stage,time_s,temperature_K,dissolved_mol,solid_mol,cosolute_mol,conditional_saturation_mol_m3,osmotic_pressure_Pa,reservoir_temperature_K,conductance_W_K,window_heat_J,window_enthalpy_change_J,energy_residual_J,accepted_steps"
    );
    for (label, seconds, drive) in schedule {
        let (step, reservoir, conductance) = match drive {
            ThermalDrive::SuppliedHeat(heat) => (
                bath.advance_with_heat(&thermal, Time::new(seconds), heat, controls, || false)?,
                String::new(),
                String::new(),
            ),
            ThermalDrive::Jacket {
                temperature,
                conductance_w_k,
            } => (
                bath.advance_with_heat_transfer(
                    &thermal,
                    Time::new(seconds),
                    temperature,
                    conductance_w_k,
                    controls,
                    || false,
                )?,
                format!("{:.3}", temperature.value()),
                format!("{conductance_w_k:.6}"),
            ),
        };
        let (next, balance) = step.into_parts();
        bath = next;
        elapsed += seconds;
        let state = bath.last_steps()[0].state();
        let aqueous = bath.solution();
        let concentration: Vec<_> = aqueous
            .solute_moles()
            .iter()
            .map(|amount| amount / aqueous.solvent_volume_m3())
            .collect();
        println!(
            "{label},{elapsed:.3},{:.9},{:.12e},{:.12e},{:.12e},{:.9},{:.9},{reservoir},{conductance},{:.9},{:.12e},{:.12e},{}",
            aqueous.temperature_k(),
            state.dissolved().value(),
            state.solid().value(),
            aqueous.solute_moles()[1],
            bath.conditional_saturation(0, &concentration)?.value(),
            aqueous.osmotic_pressure_pa()?,
            balance.heat_input_j(),
            balance.enthalpy_change_j(),
            balance.energy_residual_j(),
            balance.accepted_steps(),
        );
    }
    Ok(())
}

fn finite_thermal_demo(
    mut bath: KineticSolubilitySolution,
    thermal: &SolubilityThermalKineticsV1,
    controls: SolubilityThermalKineticControlsV1,
) -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "stage,contact_time_s,solution_temperature_K,coolant_temperature_K,dissolved_mol,solid_mol,cosolute_mol,window_heat_into_solution_J,solution_enthalpy_change_J,coolant_enthalpy_change_J,closed_energy_residual_J,coolant_energy_since_preparation_J,accepted_steps"
    );
    // Each block connects a newly prepared finite coolant charge. Within a
    // block, all four windows retain that SAME reservoir's original reference
    // and accumulated energy. The solution and its particle stocks continue
    // across both blocks; replacing coolant is an explicit apparatus action.
    for (label, initial_coolant_k) in [("finite_cooling", 292.0), ("finite_heating", 310.0)] {
        let mut coolant = SensibleHeatReservoirV1::new(
            Temperature::new(initial_coolant_k),
            4000.0, // Complete explicitly authored coolant heat capacity [J/K].
        )?;
        let mut elapsed = 0.0;
        for seconds in [30.0, 120.0, 300.0, 900.0] {
            let step = bath.advance_with_finite_heat_reservoir(
                thermal,
                Time::new(seconds),
                &coolant,
                25.0, // Complete fixed contact conductance [W/K].
                controls,
                || false,
            )?;
            let (next_bath, outcome) = step.into_parts();
            let coolant_change = outcome.reservoir_enthalpy_change_j();
            let closed_residual = outcome.closed_system_energy_residual_j();
            let (balance, next_coolant) = outcome.into_parts();
            bath = next_bath;
            coolant = next_coolant;
            elapsed += seconds;
            let phase = bath.last_steps()[0].state();
            println!(
                "{label},{elapsed:.3},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{coolant_change:.12e},{closed_residual:.12e},{:.12e},{}",
                bath.solution().temperature_k(),
                coolant.temperature().value(),
                phase.dissolved().value(),
                phase.solid().value(),
                bath.solution().solute_moles()[1],
                balance.heat_input_j(),
                balance.enthalpy_change_j(),
                coolant.energy_change_j(),
                balance.accepted_steps(),
            );
        }
    }
    Ok(())
}
