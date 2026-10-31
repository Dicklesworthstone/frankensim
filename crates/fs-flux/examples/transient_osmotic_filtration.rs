//! Finite film storage and startup compared with the steady-film approximation.
//! Run: cargo run -p fs-flux --example transient_osmotic_filtration
//! Geometry and coefficients are illustrative, not a calibrated filtration unit.
//!
//! Both cases start with the SAME total water and solute. The steady comparison
//! puts each film's initial inventory in its adjacent mixed bulk; the transient
//! model keeps that water and solute in 32 actual film cells. Compare transferred
//! water and film profiles. The transient step is first-order implicit with an
//! explicit one-second step; refine both step size and mesh for accuracy studies.

use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticMembrane, OsmoticNetwork, OsmoticNetworkError,
    OsmoticPressureVolume, OsmoticStepControls,
    transient::{TransientFilm, TransientOsmoticPair, TransientStepControls},
};
use fs_material::osmosis::{
    HydraulicPermeance, IdealDiluteSolution, MembraneTransport, ReflectionCoefficient,
    SolutePermeance,
    polarization::{FilmPermeance, PolarizedMembraneTransport},
};

fn bulk(volume: f64, concentration: f64) -> Result<OsmoticCompartment, OsmoticNetworkError> {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(298.15, 18e-6, volume, vec![volume * concentration])?,
        OsmoticPressureVolume::new(volume, 101_325.0, 0.0)?,
        1e-9,
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let area_m2 = 1e-4;
    let thickness_m = 1e-3;
    let diffusivity_m2_s = 1e-9;
    let bulk_volume_m3 = 1e-6;
    let film_cells = 32;
    let material = MembraneTransport::new(
        HydraulicPermeance::new(2e-11)?,
        vec![(ReflectionCoefficient::new(1.0)?, SolutePermeance::new(0.0)?)],
    )?;
    let mut transient = TransientOsmoticPair::new(
        [bulk(bulk_volume_m3, 50.0)?, bulk(bulk_volume_m3, 100.0)?],
        area_m2,
        material.clone(),
        Some(TransientFilm::new(
            thickness_m,
            film_cells,
            vec![diffusivity_m2_s],
            vec![50.0; film_cells],
        )?),
        Some(TransientFilm::new(
            thickness_m,
            film_cells,
            vec![diffusivity_m2_s],
            vec![100.0; film_cells],
        )?),
        1 << 20,
    )?;
    let mixed_volume_m3 = bulk_volume_m3 + area_m2 * thickness_m;
    let film = FilmPermeance::from_diffusion_layer(diffusivity_m2_s, thickness_m)?;
    let mut steady = OsmoticNetwork::new(
        vec!["retained-neutral-solute".into()],
        vec![bulk(mixed_volume_m3, 50.0)?, bulk(mixed_volume_m3, 100.0)?],
        vec![OsmoticMembrane::new_polarized(
            0,
            1,
            area_m2,
            PolarizedMembraneTransport::new(material, vec![Some(film)], vec![Some(film)])?,
        )?],
        1 << 20,
    )?;
    let steady_controls = OsmoticStepControls {
        maximum_step_s: 60.0,
        relative_tolerance: 1e-6,
        ..OsmoticStepControls::default()
    };
    let initial_water = transient.total_solvent_volume_m3();
    let mut initial_solute = [0.0];
    transient.total_solute_moles_into(&mut initial_solute)?;
    let mut last_faces = None;
    println!(
        "time_s,transient_water_transfer_ul,steady_water_transfer_ul,left_film_solute_umol,right_film_solute_umol,left_membrane_face_mM,right_membrane_face_mM"
    );
    for time in [0.0, 1.0, 10.0, 100.0, 1000.0, 6000.0] {
        while transient.time_s() < time {
            let step = (time - transient.time_s()).min(1.0);
            let report = transient.step(step, TransientStepControls::default(), || false)?;
            last_faces = Some((
                report.membrane_left_concentration_mol_m3[0],
                report.membrane_right_concentration_mol_m3[0],
            ));
        }
        if steady.time_s() < time {
            steady.advance(time - steady.time_s(), steady_controls, || false)?;
        }
        let stored = |side| {
            (0..film_cells)
                .map(|cell| transient.film_solute_moles(side, cell).unwrap()[0])
                .sum::<f64>()
                * 1e6
        };
        let faces = last_faces.map_or_else(
            || ",".to_string(),
            |(left, right)| format!("{left:.9},{right:.9}"),
        );
        println!(
            "{time:.0},{:.9},{:.9},{:.9},{:.9},{faces}",
            (transient.bulks()[1].solution().solvent_volume_m3() - bulk_volume_m3) * 1e9,
            (steady.compartments()[1].solution().solvent_volume_m3() - mixed_volume_m3) * 1e9,
            stored(0),
            stored(1),
        );
    }
    let mut final_solute = [0.0];
    transient.total_solute_moles_into(&mut final_solute)?;
    eprintln!(
        "bulk-plus-film water balance = {:.3e} m³; solute balance = {:.3e} mol",
        transient.total_solvent_volume_m3() - initial_water,
        final_solute[0] - initial_solute[0]
    );
    Ok(())
}
