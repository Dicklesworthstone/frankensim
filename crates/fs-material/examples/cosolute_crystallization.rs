//! A neutral cosolute changes saturation, then dilution redissolves the crystals.
//!
//! Run: cargo run -p fs-material --example cosolute_crystallization
//! Values are authored demonstration coefficients, not measurements for a named
//! substance. The second species has zero phase-transfer conductance, so its
//! dissolved stock changes neither during crystallization nor during imposed
//! dilution. Its actual concentration still changes the first species' activity.
//! The fixed-temperature model supplies no excess calorics or heat balance.

use std::sync::Arc;

use fs_material::osmosis::{
    nonideal::VirialSolutionModel, solubility::kinetics::KineticSolubilitySolution,
};
use fs_qty::{Amount, AmountConcentration, Temperature, Time, Volume, VolumetricFlowRate};
use fs_thermochem::{
    SpeciesId,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
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
