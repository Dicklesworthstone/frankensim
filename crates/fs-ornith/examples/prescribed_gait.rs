//! Run with: cargo run -p fs-ornith --example prescribed_gait --release
use fs_bem::wake2d::PanelWakeTransportOptions;
use fs_ornith::OrnithCandidate;
use fs_ornith::screen::{GaitScreenControls, GaitTransport, SectionGait, screen_section_gait};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let candidate = OrnithCandidate::from_genes(&[0.5, 0.4, 0.5, 0.5, 0.5]);
    let gait = SectionGait {
        heave_amplitude_chord: 0.02,
        pitch_amplitude_rad: 0.015,
        reduced_frequency: 1.0,
        heave_phase_rad: 0.0,
        pitch_phase_rad: std::f64::consts::FRAC_PI_2,
        pivot_chord: 0.25,
    };
    let result = screen_section_gait(
        &candidate,
        gait,
        GaitScreenControls {
            panels: 32,
            intervals: 80,
            interval_nondimensional: 0.04,
            core_radius_chord: 0.05,
            transport: GaitTransport::Fixed(PanelWakeTransportOptions {
                substeps: 8,
                max_coupled_solves: 17,
                max_displacement_fraction: Some(0.8),
            }),
            maximum_coupled_solves: 1360,
        },
        || false,
    )?;
    let last = result
        .wake
        .history()
        .last()
        .expect("positive interval count");
    println!("2-D prescribed heave/pitch; fixed c, U_ref and c/U_ref units");
    println!(
        "time={:.3}, bound circulation={:.6}, physical wake packets={}",
        last.t,
        last.bound,
        result.wake.wake().len()
    );
    println!(
        "circulation-weighted streamwise travel speed={:.6} U_ref",
        result.circulation_weighted_streamwise_speed
    );
    println!(
        "maximum Kelvin residual={:.3e}, coupled solves={}",
        result.maximum_kelvin_residual, result.coupled_solves
    );
    println!(
        "This kinematic wake metric does not predict thrust, efficiency or unsteady pressure loads."
    );
    Ok(())
}
