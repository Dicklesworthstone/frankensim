//! Finite film storage and startup compared with the steady-film approximation.
//! Run: cargo run -p fs-flux --example transient_osmotic_filtration
//! Geometry and coefficients are illustrative, not a calibrated filtration unit.
//!
//! Both cases start with the SAME total water and solute. The steady comparison
//! puts each film's initial inventory in its adjacent mixed bulk; the transient
//! model keeps that water and solute in 32 nonuniform cells per film. The three
//! layers have different diffusivities; the steady comparison uses their true
//! series resistance, sum(dx/D), rather than an averaged diffusivity.
//! `advance` adapts first-order implicit steps to local amount, water-volume and
//! membrane-face errors. Spatial refinement remains a separate accuracy study.
//! CSV work counts include coarse and rejected trials. Equal prescribed bulk
//! pressures make mechanical pressure work zero here; it is not an account of
//! osmotic free-energy loss or heat production.

use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticMembrane, OsmoticNetwork, OsmoticNetworkError,
    OsmoticPressureVolume, OsmoticStepControls,
    transient::{
        TransientAdvanceControls, TransientAdvanceReport, TransientFilm, TransientOsmoticPair,
        TransientStepControls,
    },
};
use fs_material::osmosis::{
    HydraulicPermeance, IdealDiluteSolution, MembraneTransport, ReflectionCoefficient,
    SolutePermeance,
    polarization::{FilmPermeance, PolarizedMembraneTransport},
};

fn bulk(volume: f64, solute_moles: f64) -> Result<OsmoticCompartment, OsmoticNetworkError> {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(298.15, 18e-6, volume, vec![solute_moles])?,
        OsmoticPressureVolume::new(volume, 101_325.0, 0.0)?,
        1e-9,
    )
}

fn layered_film(
    concentration: f64,
    right_side: bool,
) -> Result<(TransientFilm, FilmPermeance), Box<dyn std::error::Error>> {
    // Illustrative layers from bulk to membrane: (thickness [m], D [m²/s], cells).
    // The small, slow cells nearest the membrane resolve concentration buildup.
    let layers = [(0.7e-3, 1e-9, 8), (0.2e-3, 4e-10, 8), (0.1e-3, 1e-10, 16)];
    let mut widths = Vec::new();
    let mut diffusivities = Vec::new();
    for (thickness, diffusivity, cells) in layers {
        widths.extend(std::iter::repeat_n(thickness / cells as f64, cells));
        diffusivities.extend(std::iter::repeat_n(diffusivity, cells));
    }
    // Both TransientFilm sides use spatial left-to-right order.
    if right_side {
        widths.reverse();
        diffusivities.reverse();
    }
    let resistance: f64 = widths
        .iter()
        .zip(&diffusivities)
        .map(|(dx, d)| dx / d)
        .sum();
    let profile = vec![concentration; widths.len()];
    Ok((
        TransientFilm::from_cells(widths, 1, diffusivities, profile)?,
        FilmPermeance::new(1.0 / resistance)?,
    ))
}

fn mixed_initial_bulk(
    pair: &TransientOsmoticPair,
    side: usize,
) -> Result<OsmoticCompartment, OsmoticNetworkError> {
    let cells = pair.film_cell_count(side).unwrap();
    let film_water: f64 = (0..cells)
        .map(|cell| pair.film_cell_volume_m3(side, cell).unwrap())
        .sum();
    let film_solute: f64 = (0..cells)
        .map(|cell| pair.film_solute_moles(side, cell).unwrap()[0])
        .sum();
    let initial = pair.bulks()[side].solution();
    bulk(
        initial.solvent_volume_m3() + film_water,
        initial.solute_moles()[0] + film_solute,
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let area_m2 = 1e-4;
    let bulk_volume_m3 = 1e-6;
    let (left_film, left_steady_permeance) = layered_film(50.0, false)?;
    let (right_film, right_steady_permeance) = layered_film(100.0, true)?;
    let material = MembraneTransport::new(
        HydraulicPermeance::new(2e-11)?,
        vec![(ReflectionCoefficient::new(1.0)?, SolutePermeance::new(0.0)?)],
    )?;
    let mut transient = TransientOsmoticPair::new(
        [
            bulk(bulk_volume_m3, 50.0 * bulk_volume_m3)?,
            bulk(bulk_volume_m3, 100.0 * bulk_volume_m3)?,
        ],
        area_m2,
        material.clone(),
        Some(left_film),
        Some(right_film),
        1 << 20,
    )?;
    let mut workspace = transient.adaptive_workspace(1 << 20)?;
    // Transfer actual initial cell inventories into the steady comparison;
    // film storage must not silently disappear from either closed system.
    let mixed = vec![
        mixed_initial_bulk(&transient, 0)?,
        mixed_initial_bulk(&transient, 1)?,
    ];
    let initial_steady_right_volume = mixed[1].solution().solvent_volume_m3();
    let mut steady = OsmoticNetwork::new(
        vec!["retained-neutral-solute".into()],
        mixed,
        vec![OsmoticMembrane::new_polarized(
            0,
            1,
            area_m2,
            PolarizedMembraneTransport::new(
                material,
                vec![Some(left_steady_permeance)],
                vec![Some(right_steady_permeance)],
            )?,
        )?],
        1 << 20,
    )?;
    let transient_controls = TransientAdvanceControls {
        maximum_step_s: 60.0,
        minimum_step_s: 1e-8,
        relative_tolerance: 1e-4,
        absolute_amount_tolerance_mol: 1e-14,
        absolute_water_tolerance_m3: 1e-15,
        absolute_face_concentration_tolerance_mol_m3: 1e-6,
        maximum_steps: 10_000,
        maximum_rejections: 128,
        // These bound algebraic root/linear error, independently of time error.
        step: TransientStepControls::default(),
    };
    let steady_controls = OsmoticStepControls {
        maximum_step_s: 60.0,
        relative_tolerance: 1e-6,
        ..OsmoticStepControls::default()
    };
    let initial_water = transient.total_solvent_volume_m3();
    let mut initial_solute = [0.0];
    transient.total_solute_moles_into(&mut initial_solute)?;
    let mut last_report: Option<TransientAdvanceReport> = None;
    let mut cumulative_steps = 0;
    let mut cumulative_roots = 0;
    let mut cumulative_pressure_work_j = 0.0;
    eprintln!(
        "series-matched film permeances: left={:.9e}, right={:.9e} m/s; initial water={initial_water:.9e} m³, solute={:.9e} mol",
        left_steady_permeance.value(),
        right_steady_permeance.value(),
        initial_solute[0],
    );
    println!(
        "time_s,transient_water_transfer_ul,steady_water_transfer_ul,left_film_solute_umol,right_film_solute_umol,left_membrane_face_mM,right_membrane_face_mM,window_accepted_steps,window_rejected_steps,window_implicit_steps,window_root_evaluations,window_maximum_scaled_error,cumulative_implicit_steps,cumulative_root_evaluations,window_pressure_work_J,cumulative_pressure_work_J"
    );
    for time in [0.0, 1.0, 10.0, 100.0, 1000.0, 6000.0] {
        if transient.time_s() < time {
            let report = transient.advance(
                time - transient.time_s(),
                transient_controls,
                &mut workspace,
                || false,
            )?;
            cumulative_steps += report.step_evaluations;
            cumulative_roots += report.root_evaluations;
            cumulative_pressure_work_j += report.pressure_work_into_solution_j;
            last_report = Some(report);
        }
        if steady.time_s() < time {
            steady.advance(time - steady.time_s(), steady_controls, || false)?;
        }
        let stored = |side| {
            (0..transient.film_cell_count(side).unwrap())
                .map(|cell| transient.film_solute_moles(side, cell).unwrap()[0])
                .sum::<f64>()
                * 1e6
        };
        let faces = last_report.as_ref().map_or_else(
            || ",".to_string(),
            |report| {
                format!(
                    "{:.9},{:.9}",
                    report.endpoint.membrane_left_concentration_mol_m3[0],
                    report.endpoint.membrane_right_concentration_mol_m3[0]
                )
            },
        );
        let (accepted, rejected, steps, roots, error, work) =
            last_report
                .as_ref()
                .map_or((0, 0, 0, 0, 0.0, 0.0), |report| {
                    (
                        report.accepted_steps,
                        report.rejected_steps,
                        report.step_evaluations,
                        report.root_evaluations,
                        report.maximum_scaled_error,
                        report.pressure_work_into_solution_j,
                    )
                });
        println!(
            "{time:.0},{:.9},{:.9},{:.9},{:.9},{faces},{accepted},{rejected},{steps},{roots},{error:.3e},{cumulative_steps},{cumulative_roots},{work:.9e},{cumulative_pressure_work_j:.9e}",
            (transient.bulks()[1].solution().solvent_volume_m3() - bulk_volume_m3) * 1e9,
            (steady.compartments()[1].solution().solvent_volume_m3() - initial_steady_right_volume)
                * 1e9,
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
