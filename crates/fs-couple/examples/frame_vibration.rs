//! Native frame -> consistent inertia -> shared exact modal history -> SI CSV.
//! `cargo run -p fs-couple --example frame_vibration > frame.csv`
//! Manufactured beam/forcing data, not a calibrated material or building.
use fs_couple::frame_vibration::{FrameResponseLimits, FrameVibrationBasis,
    FrameVibrationBudget, FrameVibrationRuntime};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, Geometry, SectionModel};
use std::io::{BufWriter, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let elements = 8;
    for i in 0..=elements { frame.add_node(i as f64/elements as f64, 0.0); }
    frame.fix(0, [true; 3]);
    for i in 0..elements {
        frame.add_element(i, i+1, 3, 2.0, &|| SectionModel::Elastic { ea: 2.0e6, ei: 80.0 })?;
    }
    let basis = FrameVibrationBasis::from_frame(&frame, 6,
        FrameVibrationBudget { max_frame_dofs: 27, max_dense_entries: 729 }, || false)?;
    eprintln!("retained {} of {} modes; vertical effective mass {} / full influence mass {} kg",
        basis.modes().len(), basis.available_mode_count(),
        basis.retained_effective_mass_kg()[1], basis.base_influence_mass_kg()[1]);
    let limits = FrameResponseLimits {
        maximum_incremental_translation_m: 0.05, maximum_incremental_rotation_rad: 0.1,
        maximum_translation_speed_m_s: 10.0, maximum_angular_speed_rad_s: 100.0,
        maximum_absolute_translation_acceleration_m_s2: 1.0e5,
        maximum_angular_acceleration_rad_s2: 1.0e6,
    };
    let mut sim = FrameVibrationRuntime::new(basis, 4000, &[0.01; 6],
        ModalAcousticTimeBudget::audible_reference(), limits,
        &[ModalAcousticState::default(); 6], || false)?;
    let mut load = vec![0.0; frame.ndof()];
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    writeln!(out, "time_s,tip_y_m,tip_vy_m_s,tip_ay_absolute_m_s2,energy_j,nodal_work_j,base_work_j,viscous_loss_j")?;
    for sample in 0..8000 {
        let time = sample as f64/4000.0;
        // Explicit sample-and-hold actuator and support records, no hidden smoothing.
        load[3*elements+1] = if sample < 80 { 1.0 } else { 0.0 };
        let ground = if sample < 4000 { 0.1*(std::f64::consts::TAU*3.0*time).sin() } else { 0.0 };
        let r = sim.step(&load, [0.0,ground], || false)?;
        writeln!(out, "{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            r.time_s, r.displacement[3*elements+1], r.relative_velocity[3*elements+1],
            r.absolute_acceleration[3*elements+1], r.retained_energy_j,
            r.nodal_work_j, r.base_inertial_work_j, r.viscous_dissipation_j)?;
    }
    out.flush()?;
    Ok(())
}
