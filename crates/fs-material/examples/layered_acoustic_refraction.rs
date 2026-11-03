//! Refraction and coherent transmission through a finite warm-air layer.
//!
//! Run: `cargo run -p fs-material --example layered_acoustic_refraction -- 353.15 0.05 1000`
//! Arguments: layer temperature [K], nominal thickness [m], frequency [Hz].
//! Both exterior half-spaces contain dry air at 293.15 K and 101325 Pa. The
//! finite layer uses the same explicit USSA-1976 calorically perfect gas model
//! at the requested temperature and pressure. The maintained planar states
//! are prescribed; no temperature evolution or interface stability is solved.
//!
//! The calculation retains coherent pressure interference and conserved
//! tangential wavenumber. A finite layer can transmit above the single-interface
//! critical angle through an evanescent field. This lossless, infinite-plane
//! example does not include absorption, finite beam diffraction or flow.

use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::{Layer, TransmittingLayers};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 3 {
        return Err(
            "usage: layered_acoustic_refraction [layer_temperature_k] [thickness_m] [frequency_hz]"
                .into(),
        );
    }
    let temperature: f64 = args.first().map_or(Ok(353.15), |value| value.parse())?;
    let thickness: f64 = args.get(1).map_or(Ok(0.05), |value| value.parse())?;
    let frequency: f64 = args.get(2).map_or(Ok(1000.0), |value| value.parse())?;
    if !thickness.is_finite()
        || thickness <= 0.0
        || !(10.0 * thickness).is_finite()
        || !frequency.is_finite()
        || frequency <= 0.0
    {
        return Err("positive finite thickness and frequency are required".into());
    }
    let model = GasSpec::dry_air_ussa1976();
    let exterior = GasState::try_new(&model, 293.15, 101_325.0)?;
    let layer = GasState::try_new(&model, temperature, 101_325.0)?;
    let exterior_acoustic = exterior.acoustic_state()?;
    let layer_acoustic = layer.acoustic_state()?;
    let speed_ratio = exterior_acoustic.sound_speed_m_s() / layer_acoustic.sound_speed_m_s();
    eprintln!(
        "USSA-1976 dry-air model: exterior 293.15 K, layer {temperature} K, pressure 101325 Pa"
    );
    eprintln!(
        "exterior: density={:.9} kg/m3, sound_speed={:.9} m/s; layer: density={:.9} kg/m3, sound_speed={:.9} m/s",
        exterior_acoustic.density_kg_m3(),
        exterior_acoustic.sound_speed_m_s(),
        layer_acoustic.density_kg_m3(),
        layer_acoustic.sound_speed_m_s(),
    );
    let mut angles = vec![0.0_f64, 30.0, 60.0, 80.0];
    if speed_ratio < 1.0 {
        let critical = speed_ratio.asin().to_degrees();
        eprintln!("single-interface critical incidence: {critical:.9} degrees");
        angles.push(critical);
        if critical > 0.5 {
            angles.push(critical - 0.5);
        }
        if critical < 89.5 {
            angles.push(critical + 0.5);
        }
    }
    angles.sort_by(f64::total_cmp);
    angles.dedup();
    println!(
        "layer_temperature_k,thickness_m,frequency_hz,incidence_deg,reflected_power_fraction,transmitted_power_fraction,absorbed_power_fraction,power_balance_error,transmitted_pressure_real,transmitted_pressure_imag"
    );
    for distance in [0.0, thickness * 0.1, thickness, thickness * 10.0] {
        let stack = TransmittingLayers::new(
            &exterior,
            &exterior,
            core::f64::consts::TAU * frequency,
            &[Layer::FluidLayer {
                fluid: &layer,
                thickness_m: distance,
            }],
        )?;
        for &angle in &angles {
            let value = stack.scattering(angle.to_radians().cos())?;
            let balance =
                value.reflected_fraction + value.transmitted_fraction + value.absorbed_fraction
                    - 1.0;
            println!(
                "{temperature:.9},{distance:.9e},{frequency:.9},{angle:.9},{:.12e},{:.12e},{:.12e},{balance:.12e},{:.12e},{:.12e}",
                value.reflected_fraction,
                value.transmitted_fraction,
                value.absorbed_fraction,
                value.pressure_transmission.re,
                value.pressure_transmission.im,
            );
        }
    }
    Ok(())
}
