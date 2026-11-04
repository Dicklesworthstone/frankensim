//! Refraction and coherent transmission through a finite warm-air layer.
//!
//! Run: `cargo run -p fs-material --example layered_acoustic_refraction -- 353.15 0.05 1000`
//! Arguments: layer temperature [K], nominal thickness [m], frequency [Hz],
//! and optional `field` or `interior` to emit a coherent pressure/velocity map
//! in the receiving half-space or inside the finite layer, respectively.
//! Both exterior half-spaces contain dry air at 293.15 K and 101325 Pa. The
//! finite layer uses the same explicit USSA-1976 calorically perfect gas model
//! at the requested temperature and pressure. The maintained planar states
//! are prescribed; no temperature evolution or interface stability is solved.
//!
//! The calculation retains coherent pressure interference and conserved
//! tangential wavenumber. A finite layer can transmit above the single-interface
//! critical angle through an evanescent field. This lossless, infinite-plane
//! example does not include absorption or flow. The `field` mode propagates
//! nine explicitly weighted plane waves together, including their interference
//! and angular filtering by the layer. This finite angular sum is not a
//! compact-aperture source or a converged continuous angular-spectrum integral.

use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::angular_spectrum::{AngularSpectrumBudget, IncidentPlaneWave};
use fs_material::porous::layers::{Layer, TransmittingLayers};
use fs_math::c64::C64;
use std::fmt::Write as _;
use std::io::Write as _;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 4
        || args
            .get(3)
            .is_some_and(|mode| mode != "field" && mode != "interior")
    {
        return Err(
            "usage: layered_acoustic_refraction [layer_temperature_k] [thickness_m] [frequency_hz] [field|interior]"
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
    if let Some(mode) = args.get(3) {
        return coherent_field(&exterior, &layer, thickness, frequency, mode == "interior");
    }
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

fn coherent_field(
    exterior: &GasState,
    layer: &GasState,
    thickness_m: f64,
    frequency_hz: f64,
    inside_layer: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let stack = TransmittingLayers::new(
        exterior,
        exterior,
        core::f64::consts::TAU * frequency_hz,
        &[Layer::FluidLayer {
            fluid: layer,
            thickness_m,
        }],
    )?;
    // An explicitly authored angular distribution. Positive amplitudes sum to
    // one pascal at the entrance origin; none is a power probability or a
    // second copy of the layer's transmission coefficient.
    let weights = [1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0];
    let components = core::array::from_fn::<_, 9, _>(|i| {
        let angle = (48.0 + 3.0 * i as f64).to_radians();
        IncidentPlaneWave {
            direction: [angle.sin(), 0.0, angle.cos()],
            pressure_pa: C64::new(weights[i] / 256.0, 0.0),
        }
    });
    let spectrum = stack.prepare_angular_spectrum(
        &components,
        AngularSpectrumBudget { max_components: 9 },
        &mut || false,
    )?;
    let interior = if inside_layer {
        Some(spectrum.prepare_layer_field(0, &mut || false)?)
    } else {
        None
    };
    let wavelength_m = exterior.acoustic_state()?.sound_speed_m_s() / frequency_hz;
    if !wavelength_m.is_finite() || wavelength_m <= 0.0 {
        return Err("the exterior wavelength must be representable".into());
    }
    eprintln!(
        "coherent field: 9 plane waves at 48..72 degrees in 3-degree steps; binomial pressure weights /256, total incident peak pressure 1 Pa at the entrance origin"
    );
    let depth_name = if inside_layer {
        eprintln!(
            "interior depth is measured in metres from the layer entrance, including both faces"
        );
        "layer_depth_m"
    } else {
        eprintln!("receiving depth is measured in metres from the layer exit");
        "receiving_depth_m"
    };
    eprintln!(
        "x is metres; time convention exp(+i omega t); intensity is computed after summing complex fields"
    );
    let mut csv = format!(
        "x_m,{depth_name},pressure_real_pa,pressure_imag_pa,velocity_x_real_m_s,velocity_x_imag_m_s,velocity_z_real_m_s,velocity_z_imag_m_s,active_x_w_m2,active_z_w_m2,reactive_x_w_m2,reactive_z_w_m2\n",
    );
    for iz in 0..=10 {
        let z = if inside_layer {
            thickness_m * (iz as f64 / 10.0)
        } else {
            0.2 * iz as f64 * wavelength_m
        };
        for ix in 0..=40 {
            let x = (-2.0 + 0.15 * ix as f64) * wavelength_m;
            let sample = if let Some(interior) = &interior {
                interior.field([x, 0.0, z], &mut || false)?
            } else {
                spectrum.transmitted_field([x, 0.0, z], &mut || false)?
            };
            let p = sample.pressure_pa;
            let vx = sample.particle_velocity_m_s[0];
            let vz = sample.particle_velocity_m_s[2];
            let active = sample.active_intensity_w_m2;
            let reactive = sample.reactive_intensity_w_m2;
            writeln!(
                csv,
                "{x:.12e},{z:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                p.re,
                p.im,
                vx.re,
                vx.im,
                vz.re,
                vz.im,
                active[0],
                active[2],
                reactive[0],
                reactive[2],
            )?;
        }
    }
    // Keep a failed field query from publishing a partial map as a complete run.
    std::io::stdout().lock().write_all(csv.as_bytes())?;
    Ok(())
}
