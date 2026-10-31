//! Finite forward/reverse osmosis with and without concentration films.
//! Run: cargo run -p fs-flux --example osmotic_filtration
//! All coefficients are illustrative inputs, not a calibrated filtration device.
//!
//! The films have negligible storage and a prescribed mass-transfer coefficient.
//! They change the transient rate; a perfectly retained solute still reaches the
//! same equilibrium. Equal pressures give final left/right volumes of 2/3 and
//! 4/3 mL. Raising the right pressure by 0.4 MPa reverses initial water transport.

use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticMembrane, OsmoticNetwork, OsmoticPressureVolume, OsmoticStepControls,
};
use fs_material::osmosis::{
    HydraulicPermeance, IdealDiluteSolution, MembraneTransport, ReflectionCoefficient,
    SolutePermeance,
    polarization::{FilmPermeance, PolarizedMembraneTransport},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let initial_volume_m3 = 1e-6;
    let membrane_area_m2 = 1e-4;
    let controls = OsmoticStepControls {
        maximum_step_s: 3600.0,
        relative_tolerance: 1e-6,
        maximum_steps: 50_000,
        ..OsmoticStepControls::default()
    };
    println!(
        "case,boundary,time_s,left_volume_ml,right_volume_ml,left_bulk_mM,right_bulk_mM,interval_mean_left_to_right_flux_m_s,left_solute_umol,right_solute_umol"
    );
    for (case, right_overpressure_pa) in [("forward", 0.0), ("reverse", 400_000.0)] {
        for with_films in [false, true] {
            let boundary = if with_films { "films" } else { "well_mixed" };
            // Explicit constant mechanical pressures; the ideal dilute solution
            // separately supplies its osmotic pressure. No pressure is counted twice.
            let cells = [
                (50.0, 101_325.0),
                (100.0, 101_325.0 + right_overpressure_pa),
            ]
            .into_iter()
            .map(|(concentration, pressure)| {
                OsmoticCompartment::new(
                    IdealDiluteSolution::new(
                        298.15,
                        18e-6,
                        initial_volume_m3,
                        vec![concentration * initial_volume_m3],
                    )?,
                    OsmoticPressureVolume::new(initial_volume_m3, pressure, 0.0)?,
                    1e-9,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
            let material = MembraneTransport::new(
                HydraulicPermeance::new(2e-11)?,
                vec![(ReflectionCoefficient::new(1.0)?, SolutePermeance::new(0.0)?)],
            )?;
            let membrane = if with_films {
                OsmoticMembrane::new_polarized(
                    0,
                    1,
                    membrane_area_m2,
                    PolarizedMembraneTransport::new(
                        material,
                        vec![Some(FilmPermeance::new(1e-6)?)],
                        vec![Some(FilmPermeance::new(1e-6)?)],
                    )?,
                )?
            } else {
                OsmoticMembrane::new(0, 1, membrane_area_m2, material)?
            };
            let mut network = OsmoticNetwork::new(
                vec!["retained-neutral-solute".into()],
                cells,
                vec![membrane],
                1 << 20,
            )?;
            let mut maximum_water_residual: f64 = 0.0;
            for time in [0.0, 600.0, 3600.0, 21_600.0, 86_400.0, 604_800.0] {
                let interval = time - network.time_s();
                // This is a conserved interval-mean flux, not an instantaneous
                // membrane-face diagnostic. The initial row has no interval.
                let mean_flux = if interval > 0.0 {
                    let before = network.compartments()[1].solution().solvent_volume_m3();
                    let report = network.advance(interval, controls, || false)?;
                    maximum_water_residual =
                        maximum_water_residual.max(report.water_balance_residual_m3.abs());
                    let after = network.compartments()[1].solution().solvent_volume_m3();
                    format!("{:.9e}", (after - before) / (membrane_area_m2 * interval))
                } else {
                    String::new()
                };
                let left = network.compartments()[0].solution();
                let right = network.compartments()[1].solution();
                println!(
                    "{case},{boundary},{time:.0},{:.9},{:.9},{:.9},{:.9},{mean_flux},{:.9},{:.9}",
                    left.solvent_volume_m3() * 1e6,
                    right.solvent_volume_m3() * 1e6,
                    left.solute_moles()[0] / left.solvent_volume_m3(),
                    right.solute_moles()[0] / right.solvent_volume_m3(),
                    left.solute_moles()[0] * 1e6,
                    right.solute_moles()[0] * 1e6,
                );
            }
            eprintln!(
                "{case}/{boundary}: maximum water balance residual = {maximum_water_residual:.3e} m³"
            );
        }
    }
    Ok(())
}
