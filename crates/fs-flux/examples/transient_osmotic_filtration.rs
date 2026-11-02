//! Sorbing finite-film startup compared with a nonsorbing steady-film reference.
//! Run: cargo run -p fs-flux --example transient_osmotic_filtration
//! Geometry and coefficients are illustrative, not a calibrated filtration unit.
//!
//! Both cases start with the SAME total water and solute. The steady comparison
//! dissolves each film's initial TOTAL inventory into its adjacent mixed bulk;
//! the transient model keeps mobile and sorbed solute in 32 nonuniform cells per
//! film. The three layers have different diffusivities and linear sorption
//! ratios R=1,2,6. Initial mobile concentration is explicitly c_total/R.
//! The steady reference uses the true mobile diffusion resistance, sum(dx/D),
//! but has no immobile storage: it is a different storage model, not the steady
//! limit of the sorbing system. Equal totals do not imply equal mobile profiles.
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
    initial_total_concentration: f64,
    right_side: bool,
) -> Result<(TransientFilm, FilmPermeance), Box<dyn std::error::Error>> {
    // Illustrative bulk-to-membrane layers: (thickness [m], D [m²/s], cells, R).
    // The small, slow cells nearest the membrane resolve concentration buildup.
    let layers = [
        (0.7e-3, 1e-9, 8, 1.0),
        (0.2e-3, 4e-10, 8, 2.0),
        (0.1e-3, 1e-10, 16, 6.0),
    ];
    let mut widths = Vec::new();
    let mut diffusivities = Vec::new();
    let mut ratios = Vec::new();
    for (thickness, diffusivity, cells, ratio) in layers {
        widths.extend(std::iter::repeat_n(thickness / cells as f64, cells));
        diffusivities.extend(std::iter::repeat_n(diffusivity, cells));
        ratios.extend(std::iter::repeat_n(ratio, cells));
    }
    // Both TransientFilm sides use spatial left-to-right order.
    if right_side {
        widths.reverse();
        diffusivities.reverse();
        ratios.reverse();
    }
    let resistance: f64 = widths
        .iter()
        .zip(&diffusivities)
        .map(|(dx, d)| dx / d)
        .sum();
    let profile = ratios
        .iter()
        .map(|ratio| initial_total_concentration / ratio)
        .collect();
    Ok((
        TransientFilm::from_cells(widths, 1, diffusivities, profile)?
            .with_linear_sorption(ratios)?,
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
    // The nonsorbing reference dissolves the same total starting inventory.
    // It does not retain the transient film's mobile/bound equilibrium split.
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
        "nonsorbing reference has matched mobile diffusion resistance: left={:.9e}, right={:.9e} m/s; initial water={initial_water:.9e} m³, total solute={:.9e} mol",
        left_steady_permeance.value(),
        right_steady_permeance.value(),
        initial_solute[0],
    );
    println!(
        "time_s,transient_water_transfer_ul,steady_water_transfer_ul,left_film_solute_umol,right_film_solute_umol,left_membrane_face_mM,right_membrane_face_mM,window_accepted_steps,window_rejected_steps,window_implicit_steps,window_root_evaluations,window_maximum_scaled_error,cumulative_implicit_steps,cumulative_root_evaluations,window_pressure_work_J,cumulative_pressure_work_J,left_film_sorbed_umol,right_film_sorbed_umol"
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
        let sorbed = |side| {
            (0..transient.film_cell_count(side).unwrap())
                .map(|cell| transient.film_sorbed_solute_moles(side, cell, 0).unwrap())
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
            "{time:.0},{:.9},{:.9},{:.9},{:.9},{faces},{accepted},{rejected},{steps},{roots},{error:.3e},{cumulative_steps},{cumulative_roots},{work:.9e},{cumulative_pressure_work_j:.9e},{:.9},{:.9}",
            (transient.bulks()[1].solution().solvent_volume_m3() - bulk_volume_m3) * 1e9,
            (steady.compartments()[1].solution().solvent_volume_m3() - initial_steady_right_volume)
                * 1e9,
            stored(0),
            stored(1),
            sorbed(0),
            sorbed(1),
        );
    }
    let mut final_solute = [0.0];
    transient.total_solute_moles_into(&mut final_solute)?;
    eprintln!(
        "bulk-plus-film water balance = {:.3e} m³; mobile-plus-sorbed solute balance = {:.3e} mol",
        transient.total_solvent_volume_m3() - initial_water,
        final_solute[0] - initial_solute[0]
    );
    Ok(())
}
