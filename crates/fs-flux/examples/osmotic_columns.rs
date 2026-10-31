//! Two finite solution columns build hydrostatic counterpressure by osmosis.
//! Run: cargo run -p fs-flux --example osmotic_columns
//! Coefficients are explicit illustrative inputs, not a calibrated membrane.

use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticMembrane, OsmoticNetwork, OsmoticPressureVolume, OsmoticStepControls,
};
use fs_material::osmosis::{
    HydraulicPermeance, IdealDiluteSolution, MembraneTransport, ReflectionCoefficient,
    SolutePermeance,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let temperature_k = 298.15;
    let solvent_molar_volume_m3_mol = 18e-6;
    let initial_volume_m3 = 1e-6;
    let column_area_m2 = 1e-5;
    let solvent_density_kg_m3 = 1000.0;
    let gravity_m_s2 = 9.81;
    let pressure = OsmoticPressureVolume::new(
        initial_volume_m3,
        101_325.0 + solvent_density_kg_m3 * gravity_m_s2 * initial_volume_m3 / column_area_m2,
        solvent_density_kg_m3 * gravity_m_s2 / column_area_m2,
    )?;
    // One neutral independent solute: 50 and 100 mol/m³ equal 50 and 100 mM.
    // Already enumerated ions would occupy their own conserved species slots.
    let cells = [50.0, 100.0]
        .into_iter()
        .map(|concentration_mol_m3| {
            OsmoticCompartment::new(
                IdealDiluteSolution::new(
                    temperature_k,
                    solvent_molar_volume_m3_mol,
                    initial_volume_m3,
                    vec![concentration_mol_m3 * initial_volume_m3],
                )?,
                pressure,
                1e-9,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let membrane = OsmoticMembrane::new(
        0,
        1,
        1e-4,
        MembraneTransport::new(
            HydraulicPermeance::new(2e-12)?,
            vec![(ReflectionCoefficient::new(1.0)?, SolutePermeance::new(0.0)?)],
        )?,
    )?;
    let mut network = OsmoticNetwork::new(
        vec!["retained-neutral-solute".into()],
        cells,
        vec![membrane],
        1 << 20,
    )?;
    let controls = OsmoticStepControls {
        maximum_step_s: 3600.0,
        relative_tolerance: 1e-6,
        maximum_steps: 50_000,
        ..OsmoticStepControls::default()
    };
    println!(
        "time_s,left_volume_ml,right_volume_ml,left_concentration_mol_m3,right_concentration_mol_m3,hydrostatic_drop_pa,osmotic_drop_pa"
    );
    for time in [0.0, 3600.0, 21_600.0, 86_400.0, 604_800.0] {
        if time > 0.0 {
            let report = network.advance(time - network.time_s(), controls, || false)?;
            eprintln!(
                "t={time:.0}s accepted={} rejected={} water_balance_m3={:.3e} solute_balance_mol={:.3e}",
                report.accepted_steps,
                report.rejected_steps,
                report.water_balance_residual_m3,
                report.solute_balance_residual_mol[0]
            );
        }
        let cells = network.compartments();
        let left = cells[0].solution();
        let right = cells[1].solution();
        println!(
            "{time:.0},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
            left.solvent_volume_m3() * 1e6,
            right.solvent_volume_m3() * 1e6,
            left.solute_moles()[0] / left.solvent_volume_m3(),
            right.solute_moles()[0] / right.solvent_volume_m3(),
            cells[0].pressure_pa()? - cells[1].pressure_pa()?,
            left.osmotic_pressure_pa()? - right.osmotic_pressure_pa()?
        );
    }
    Ok(())
}
