//! A file-independent material-to-duct impedance and power sweep.
//!
//! Example with explicit terrestrial-air composition and the existing ISO rates:
//! `cargo run -p fs-duct --example molecular_duct -- \
//!   'mixture 293.15 101325 mole 4 N2 0.781 O2 0.209 Ar 0.0096 CO2 0.0004' \
//!   'iso9613-terrestrial-air'`
//!
//! Both arguments use the existing material input grammar. Change the actual
//! gas state, composition, reservoir strength or rate declaration to change the
//! solve. A `declared` molecular law needs every species' relaxation entry, including explicit zero
//! reservoirs for a nonrelaxing species. No gas or rate is inferred. The example
//! composes existing engineering models; it is not experimental validation.
//!
//! Geometry is a 1 mm-radius, 0.2 m rigid cylinder with coupled circular wall
//! diffusion. Its outlet is the explicitly declared equilibrium characteristic
//! RESISTANCE, not an inferred radiation load or a perfectly matched dispersive
//! termination. Reported duct dissipation is net input minus output power.

use fs_duct::compose::ComposeOptions;
use fs_duct::molecular::{MolecularDuct, MolecularWallResolution};
use fs_duct::{ConeSlices, Duct, Segment};
use fs_material::gas_species::input::read_gas_record;
use fs_material::relaxation::input::{medium_from_record, read_relaxation};
use fs_math::c64::C64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 {
        return Err("usage: molecular_duct '<gas declaration>' '<relaxation declaration>'".into());
    }
    let mut gas_fields = args[0].split_ascii_whitespace();
    let gas = read_gas_record(&mut gas_fields, None, None)?;
    if gas_fields.next().is_some() {
        return Err("extra gas fields".into());
    }
    let mut law_fields = args[1].split_ascii_whitespace();
    let law = read_relaxation(&mut law_fields)?;
    if law_fields.next().is_some() {
        return Err("extra relaxation fields".into());
    }
    let state = law.evaluate(&medium_from_record(&gas)?)?;
    let radius = 0.001;
    let area = core::f64::consts::PI * radius * radius;
    let load = gas.state().characteristic_impedance / area;
    let model = MolecularDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 0.2,
            }],
        },
        state,
        Some(MolecularWallResolution {
            poles_per_field: 12,
            band_rad_s: [
                core::f64::consts::TAU * 100.0,
                core::f64::consts::TAU * 5000.0,
            ],
            maximum_relative_error: 0.01,
        }),
        ConeSlices::Fixed(1),
        1,
    )?;
    let inlet_flow_amplitude = 1e-6; // peak m^3/s, retained at every frequency
    println!(
        "frequency_hz,zin_real_pa_s_m3,zin_imag_pa_s_m3,outlet_flow_peak_m3_s,input_power_w,output_power_w,duct_dissipation_w"
    );
    for frequency in [100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0] {
        let response = model.response(
            core::f64::consts::TAU * frequency,
            Some(C64::from_re(load)),
            ComposeOptions::default(),
        )?;
        let outlet_flow = inlet_flow_amplitude * response.u_out_over_u_in.abs();
        let input_power = 0.5 * response.z_in.re * inlet_flow_amplitude.powi(2);
        let output_power = 0.5 * load * outlet_flow.powi(2);
        println!(
            "{frequency:.1},{:.12e},{:.12e},{outlet_flow:.12e},{input_power:.12e},{output_power:.12e},{:.12e}",
            response.z_in.re,
            response.z_in.im,
            input_power - output_power
        );
    }
    Ok(())
}
