//! Material-scaled throughflow with an incoming acoustic pulse and signed
//! boundary exchange on the original compressible NSF grid.
//!
//! cargo run -p fs-lbm --example nsf_open_channel -- 100000 300 5 240
//! Arguments: reference pressure [Pa], temperature [K], throughflow [m/s],
//! steps. Pressure and temperature re-evaluate the shared dry-air GasSpec;
//! no density, sound speed, viscosity or conductivity is authored separately.
//! The same grid API accepts a retained ThermochemicalNsfGas (NASA9 mixture
//! and source-selected local transport); see nsf3_open_boundaries.rs.
//! This is the local linear characteristic boundary rung, not a resolved
//! instrument or a qualified -40 dB flow/DG interface.

use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfGrid3, NsfSpec, WallThermal,
};
use fs_material::gas::{GasSpec, GasState};
use std::{error::Error, io};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let pressure: f64 = args.next().map_or(Ok(100000.0), |s| s.parse())?;
    let temperature: f64 = args.next().map_or(Ok(300.0), |s| s.parse())?;
    let speed: f64 = args.next().map_or(Ok(5.0), |s| s.parse())?;
    let steps: usize = args.next().map_or(Ok(240), |s| s.parse())?;
    if args.next().is_some() || !(1..=2000).contains(&steps) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "usage: nsf_open_channel [pressure_pa temperature_k speed_m_s steps(1..2000)]").into());
    }
    let gas = GasSpec::dry_air_ussa1976();
    let state = GasState::try_new(&gas, temperature, pressure)?;
    let dx = 5e-5;
    let velocity_unit = (3.0 * state.specific_gas_constant * temperature).sqrt();
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: pressure,
        reference_temperature_k: temperature,
        reference_velocity_m_s: [speed, 0.0, 0.0],
        incoming_characteristic_m_s: 0.0,
    };
    let mut grid = NsfGrid3::<D3q27>::new(NsfSpec {
        dims: [96, 4, 4], periodic: [false, true, true], gas,
        temperature_k: temperature, pressure_pa: pressure, dx_m: dx,
        temperature_window_k: [0.95 * temperature, 1.05 * temperature],
        maximum_mach: 0.2, maximum_temperature_ratio: 0.05,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 near 293 K; held constant for this example",
        maximum_cells: 1536, wall: WallThermal::Adiabatic,
        workers: 1, transport_points: 17,
    }, |_| false, |_| (1.0, [speed / velocity_unit, 0.0, 0.0], temperature))?
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.02,
        })?;
    let initial_mass = grid.total_mass_kg();
    let initial_energy = grid.total_energy_j();
    let (mut mass_out, mut energy_out) = (0.0, 0.0);
    println!("step,time_s,incoming_low_m_s,low_p_pa,low_q_out_m3_s,high_p_pa,high_q_out_m3_s,mass_residual_kg,energy_residual_j");
    for step in 0..steps {
        let incoming = if step < 64 {
            -0.1 * (1.0 - fs_math::det::cos(2.0 * core::f64::consts::PI * step as f64 / 64.0))
        } else { 0.0 };
        grid.set_incoming_characteristics([incoming, 0.0])?;
        grid.step()?;
        let exchange = grid.last_boundary_exchange().expect("one accepted open step");
        mass_out += exchange.iter().map(|f| f.mass_out_kg).sum::<f64>();
        energy_out += exchange.iter().map(|f| f.total_energy_out_j).sum::<f64>();
        if step % 16 == 0 || step + 1 == steps {
            println!("{},{:.9e},{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.6e},{:.6e}",
                step + 1, (step + 1) as f64 * grid.law().dt_s, incoming,
                exchange[0].mean_pressure_pa, exchange[0].normal_volume_flow_m3_s,
                exchange[1].mean_pressure_pa, exchange[1].normal_volume_flow_m3_s,
                grid.total_mass_kg() - initial_mass + mass_out,
                grid.total_energy_j() - initial_energy + energy_out);
        }
    }
    Ok(())
}
