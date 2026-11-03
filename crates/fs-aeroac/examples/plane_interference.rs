//! Direct/reflected interference from a prescribed volume source.
//!
//! cargo run -p fs-aeroac --example plane_interference -- rigid
//! cargo run -p fs-aeroac --example plane_interference -- release
//!
//! All numbers are declared illustrative SI inputs, not calibrated source data.
use fs_aeroac::curle3d::{
    CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary, PreparedCompactSourceObserver,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let boundary = match std::env::args().nth(1).as_deref() {
        Some("rigid") => IdealPlaneBoundary::Rigid,
        Some("release") => IdealPlaneBoundary::PressureRelease,
        _ => return Err("usage: plane_interference rigid|release".into()),
    };
    let medium = FluidAcousticState::try_new(1.25, 340.0)?;
    let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundary)?;
    let source = CompactSourcePhasor::monopole(C64::from_re(1e-4));
    println!("# Authored point source: Q_peak=1e-4 m^3/s, z_source=1 m, z_receiver=3 m");
    println!("# Medium: rho=1.25 kg/m^3, c=340 m/s; infinite {boundary:?} plane at z=0");
    println!("# Paths: direct=2 m, reflected=4 m. Source phases refer to one emission clock.");
    println!(
        "frequency_hz,free_pressure_rms_pa,half_space_pressure_rms_pa,pressure_power_ratio,active_z_w_m2,reactive_z_w_m2"
    );
    for frequency in [42.5, 85.0, 127.5, 170.0, 255.0, 340.0, 425.0, 510.0] {
        let observer = PreparedCompactSourceObserver::new(
            medium,
            frequency,
            [0.0, 0.0, 3.0],
            &[[0.0, 0.0, 1.0]],
            1,
        )?
        .with_plane(plane, 2)?;
        let free = observer.direct().evaluate(&[source])?;
        let reflected = observer.evaluate_independent_sources(&[source])?;
        let free_square = 0.5 * free.pressure_pa.norm_sq();
        println!(
            "{frequency:.1},{:.9e},{:.9e},{:.9},{:.9e},{:.9e}",
            free_square.sqrt(),
            reflected.mean_square_pressure_pa2.sqrt(),
            reflected.mean_square_pressure_pa2 / free_square,
            reflected.active_intensity_w_m2[2],
            reflected.reactive_intensity_w_m2[2]
        );
    }
    Ok(())
}
