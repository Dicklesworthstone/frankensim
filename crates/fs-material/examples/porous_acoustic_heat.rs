//! Coherent pressure, velocity and local acoustic heating in a finite foam layer.
//!
//! cargo run -p fs-material --example porous_acoustic_heat -- 1000 0.06 1 0
//! Arguments: frequency [Hz], thickness [m], authored peak-pressure scale [Pa],
//! and relative phase [rad]. Defaults are 1000, 0.06, 1 and 0, respectively.
//! Two coherent plane waves have directions (+/-0.6, 0, 0.8), each with half
//! the pressure scale. Changing their relative phase moves the interference
//! pattern; it does not change either wave's separately supplied amplitude.
//!
//! The existing Basotect Foam-02 card supplies the rigid porous skeleton. Its
//! characterization viscosity is explicitly assumed to be 1.81e-5 Pa s here.
//! The pore fluid and both half-spaces use the same USSA-1976 dry-air model at
//! 293.15 K and 101325 Pa. Source frequency limits remain enforced by the law.
//! These are authored linear-model inputs, not a new material calibration.
//!
//! CSV covers one transverse interference period and the complete layer depth.
//! Heat uses the original complex density/bulk modulus AFTER coherent field
//! superposition, retaining transverse viscous loss even at pressure nodes.
//! The original finite-layer propagation includes both interfaces, standing
//! waves and attenuation. No finite-aperture diffraction, frame motion,
//! temperature rise or thermal feedback is modeled. Spatial sample counts are
//! visualization choices and do not certify a resolved heating integral.

use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::{
    AngularSpectrumBudget, IncidentPlaneWave, Layer, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::{c64::C64, det};
use std::fmt::Write as _;
use std::io::Write as _;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 4 {
        return Err("usage: porous_acoustic_heat [frequency_hz] [thickness_m] [pressure_pa] [relative_phase_rad]".into());
    }
    let frequency: f64 = args.first().map_or(Ok(1000.0), |v| v.parse())?;
    let thickness: f64 = args.get(1).map_or(Ok(0.06), |v| v.parse())?;
    let pressure: f64 = args.get(2).map_or(Ok(1.0), |v| v.parse())?;
    let relative_phase: f64 = args.get(3).map_or(Ok(0.0), |v| v.parse())?;
    if !frequency.is_finite()
        || frequency <= 0.0
        || !thickness.is_finite()
        || thickness <= 0.0
        || !pressure.is_finite()
        || pressure < 0.0
        || !relative_phase.is_finite()
    {
        return Err("require finite positive frequency/thickness, nonnegative pressure and finite relative phase".into());
    }
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)?;
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5)?;
    let stack = TransmittingLayers::new(
        &air,
        &air,
        core::f64::consts::TAU * frequency,
        &[Layer::Porous {
            law: &law,
            thickness_m: thickness,
        }],
    )?;
    let (sin, cos) = det::sin_cos(relative_phase);
    let waves = [
        IncidentPlaneWave {
            direction: [0.6, 0.0, 0.8],
            pressure_pa: C64::from_re(0.5 * pressure),
        },
        IncidentPlaneWave {
            direction: [-0.6, 0.0, 0.8],
            pressure_pa: C64::new(cos, sin).scale(0.5 * pressure),
        },
    ];
    let spectrum = stack.prepare_angular_spectrum(
        &waves,
        AngularSpectrumBudget { max_components: 2 },
        &mut || false,
    )?;
    let interior = spectrum.prepare_layer_field(0, &mut || false)?;
    let period = air.acoustic_state()?.sound_speed_m_s() / (1.2 * frequency);
    if !period.is_finite() || period <= 0.0 {
        return Err("the transverse interference period must be representable".into());
    }
    eprintln!(
        "# Basotect Foam-02; reference viscosity=1.81e-5 Pa s; USSA-1976 dry air at 293.15 K, 101325 Pa"
    );
    eprintln!(
        "# frequency_hz={frequency},thickness_m={thickness},pressure_scale_pa={pressure},relative_phase_rad={relative_phase},transverse_period_m={period}"
    );
    // Buffer this bounded map so a refused observation emits no partial CSV.
    let mut csv = String::from(
        "x_m,depth_m,pressure_real_pa,pressure_imag_pa,velocity_x_real_m_s,velocity_x_imag_m_s,velocity_y_real_m_s,velocity_y_imag_m_s,velocity_z_real_m_s,velocity_z_imag_m_s,active_x_w_m2,active_y_w_m2,active_z_w_m2,reactive_x_w_m2,reactive_y_w_m2,reactive_z_w_m2,viscous_heat_w_m3,thermal_heat_w_m3\n",
    );
    for iz in 0..=20 {
        let depth = thickness * (f64::from(iz) / 20.0);
        for ix in 0..=40 {
            let x = period * (f64::from(ix) / 40.0);
            let sample = interior.field_with_heat([x, 0.0, depth], &mut || false)?;
            write!(csv, "{x:.12e},{depth:.12e}")?;
            for value in
                core::iter::once(sample.field.pressure_pa).chain(sample.field.particle_velocity_m_s)
            {
                write!(csv, ",{:.12e},{:.12e}", value.re, value.im)?;
            }
            for value in sample
                .field
                .active_intensity_w_m2
                .into_iter()
                .chain(sample.field.reactive_intensity_w_m2)
                .chain([sample.viscous_heating_w_m3, sample.thermal_heating_w_m3])
            {
                write!(csv, ",{value:.12e}")?;
            }
            writeln!(csv)?;
        }
    }
    std::io::stdout().lock().write_all(csv.as_bytes())?;
    Ok(())
}
