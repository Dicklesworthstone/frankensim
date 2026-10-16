//! Run the same heated, force-driven fluid with data-only material substitutions.
//!
//! cargo run -p fs-lbm --example material_heated_flow -- \
//!   crates/fs-lbm/examples/material_heated_flow/high_capacity.txt 200
//! Substitute low_capacity.txt without changing any physics code.
//! Synthetic frozen-property BGK/Boussinesq slab, not a calibrated substance,
//! resolved combustion, melting, electromagnetic source model or speed claim.

use std::{error::Error, io};
use fs_lbm::material::{MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

fn usage(message: &str) -> io::Error { io::Error::new(io::ErrorKind::InvalidInput, message) }

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let path = args.next().ok_or_else(|| usage("usage: material_heated_flow MATERIAL.txt [STEPS=200]"))?;
    let steps: usize = args.next().map_or(Ok(200), |s| s.parse())?;
    if args.next().is_some() || steps < 2 || steps > 5000 {
        return Err(usage("expected one material file and a step count in 2..=5000").into());
    }
    let text = std::fs::read_to_string(path)?;
    let values = text.lines().flat_map(|line| line.split('#').next().unwrap_or("").split_whitespace())
        .map(str::parse::<f64>).collect::<Result<Vec<_>, _>>()?;
    if values.len() != 4 { return Err(usage("material file requires exactly rho, cp, conductivity, viscosity in SI units").into()); }
    let material = FluidTransportState::try_new(values[0], values[1], values[2], values[3])?;
    let dx = 0.001; let dt = 0.001;
    let plan = plan_scaling_from_material(&material, dx, dt, 0.01)?;
    let mut slab = plan.thermal_slab3(MaterialThermalSlab3 {
        nx: 4, ny: 4, nz: 12,
        bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 9.81, expansion_per_k: 2e-4, max_cells: 224,
    })?;
    let cells = slab.lattice.gpop.len(); let plane = 16; let fluid = plane..cells-plane;
    let force = vec![[2.0, 0.0, 0.0]; cells];
    let mut power = vec![0.0; cells];
    let limits = HeatSourceLimits { max_cells: 224, temperature_range_k: [290.0, 330.0], balance_tolerance_j: 1e-9 };
    let mut supplied = 0.0; let mut stored = 0.0; let mut boundary = 0.0;
    println!("step,time_s,mean_temperature_k,max_speed_m_s,source_j,stored_change_j,bottom_wall_j,top_wall_j,defect_j");
    for step in 0..steps {
        power[fluid.clone()].fill(if step < steps/2 { 2e7 } else { 0.0 });
        let report = slab.step_with_power_density(&power, &force, limits)?;
        supplied += report.source_j; stored += report.stored_change_j;
        boundary += report.bottom_wall_j + report.top_wall_j;
        let mut temperature = 0.0; let mut max_speed = 0.0_f64;
        for cell in fluid.clone() {
            temperature += slab.temperature_si(slab.lattice.gpop[cell].iter().sum())?;
            let (_, u) = slab.lattice.grid.moments3(cell);
            max_speed = max_speed.max(plan.velocity_si(u.iter().map(|x| x*x).sum::<f64>().sqrt())?);
        }
        if max_speed > 0.01 { return Err(usage("trajectory exceeded the declared physical speed bound").into()); }
        if step % 20 == 0 || step+1 == steps {
            println!("{},{:.6},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.3e}",
                step+1, (step+1) as f64*dt, temperature/fluid.len() as f64, max_speed,
                report.source_j, report.stored_change_j, report.bottom_wall_j, report.top_wall_j, report.defect_j);
        }
    }
    eprintln!("source={supplied:.12e} J; storage={stored:.12e} J; walls={boundary:.12e} J; residual={:.3e} J", stored-supplied-boundary);
    Ok(())
}
