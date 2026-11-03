//! A broadband acoustic packet crosses an original gas/gas or gas/water interface.
//!
//! cargo run -p fs-couple --example fluid_interface_pulse -- air-water 12 1600 1e-6 8
//! Arguments: CASE GRID STEPS DT_S ANGLE_DEG. CASE is air-water or cold-hot-air.
//! GRID must be even so the interface x=0.5 m conforms to tetrahedron faces.
//!
//! The one-metre cube starts with a compact pressure-potential packet in its
//! left fluid. Its propagation direction is tilted by ANGLE_DEG in the xy plane.
//! Dry air uses the existing USSA1976 gas specification at 100 kPa and 300 K;
//! the right side uses either the same gas at 600 K or original IF97 region-1
//! water at 100 kPa and 300 K. No density or sound speed is independently fitted.
//! All exterior faces have local real admittance 1/(rho*c) from their own fluid.
//! That outlet is exact for normally incident plane waves only, so oblique waves
//! may reflect there. The internal fluid interface is resolved by the full FE
//! pressure space. This example does not impose a ray/Snell approximation.
//!
//! Stdout contains the final original nodal pressure/potential field as CSV;
//! stderr contains time, physical energy and independently evaluated step work.
//! Refine GRID and DT_S independently: implicit stability is not phase accuracy,
//! and the compact pulse needs several elements across its spatial width. The
//! materials and interface are frozen; this does not evolve heat or mass transfer.

use std::error::Error;

use fs_couple::cavity_modal::interior::heterogeneous::regions::FluidRegionsInterior;
use fs_feec::acoustics::helmholtz::heterogeneous::time::{
    HeterogeneousTimeLimits, ResistiveBoundary,
};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97State};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        return Err(
            "usage: fluid_interface_pulse air-water|cold-hot-air GRID STEPS DT_S ANGLE_DEG".into(),
        );
    }
    let grid: usize = args[2].parse()?;
    let steps: usize = args[3].parse()?;
    let dt: f64 = args[4].parse()?;
    let angle_deg: f64 = args[5].parse()?;
    if !(4..=24).contains(&grid)
        || grid % 2 != 0
        || !(1..=200_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || !(dt * steps as f64).is_finite()
        || !angle_deg.is_finite()
        || angle_deg.abs() > 10.0
    {
        return Err("require even GRID in 4..=24, STEPS in 1..=200000, finite positive duration, and |ANGLE_DEG| <= 10".into());
    }
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0)?;
    let right: Box<dyn FluidAcoustics> = match args[1].as_str() {
        "cold-hot-air" => Box::new(GasState::try_new(
            &GasSpec::dry_air_ussa1976(),
            600.0,
            100_000.0,
        )?),
        "air-water" => Box::new(If97State::at_pressure_temperature(
            Pressure::new(100_000.0),
            Temperature::new(300.0),
            If97Region::Region1,
        )?),
        _ => return Err("CASE must be air-water or cold-hot-air".into()),
    };
    let (complex, positions) = fs_feec::fixtures::kuhn_cube(grid);
    let left = gas.acoustic_state()?;
    let geometry = HelmholtzFe::new(
        &complex,
        &positions,
        1,
        Medium {
            density_kg_m3: left.density_kg_m3(),
            sound_speed_m_s: left.sound_speed_m_s(),
        },
    )?;
    let element_regions: Vec<_> = geometry
        .tetrahedra()
        .iter()
        .map(|tet| {
            let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
            usize::from(x > 0.5)
        })
        .collect();
    let acoustic = FluidRegionsInterior::new(
        &geometry,
        &[&gas, right.as_ref()],
        &element_regions,
        3,
        &mut || false,
    )?;
    let mut boundaries = Vec::with_capacity(2);
    for region in 0..2 {
        let state = acoustic
            .region_state(region)
            .ok_or("unused example fluid")?;
        let faces = geometry
            .boundary_faces()
            .into_iter()
            .filter(|face| {
                let x = face.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 3.0;
                usize::from(x > 0.5) == region
            })
            .collect();
        boundaries.push(ResistiveBoundary {
            faces,
            admittance_m_pa_s: 1.0 / state.characteristic_impedance_pa_s_m(),
        });
        eprintln!(
            "# region={region},rho_kg_m3={},c_m_s={}",
            state.density_kg_m3(),
            state.sound_speed_m_s()
        );
    }
    let prepared = acoustic.operator().prepare_time(
        dt,
        &boundaries,
        HeterogeneousTimeLimits {
            max_dofs: 16_000,
            max_factor_bytes: 512 * 1024 * 1024,
            ..HeterogeneousTimeLimits::default()
        },
        &mut || false,
    )?;
    let angle = angle_deg * core::f64::consts::PI / 180.0;
    let direction = [det::cos(angle), det::sin(angle)];
    let width = 0.15;
    let amplitude_pa = 0.25;
    let mut potential = Vec::with_capacity(geometry.ndof());
    let mut pressure = Vec::with_capacity(geometry.ndof());
    // psi=A*w/c*f(q), p=-A*f'(q), f(q)=(1-q²)^4 on |q|<1.
    // This is the right-going homogeneous acoustic solution before the packet
    // meets an interface/exterior face. Its complete initial support lies in air.
    for &[x, y, _] in &positions {
        let q = (direction[0] * (x - 0.25) + direction[1] * (y - 0.5)) / width;
        if q.abs() < 1.0 {
            let envelope = 1.0 - q * q;
            potential.push(amplitude_pa * width / left.sound_speed_m_s() * envelope.powi(4));
            pressure.push(8.0 * amplitude_pa * q * envelope.powi(3));
        } else {
            potential.push(0.0);
            pressure.push(0.0);
        }
    }
    let mut state = prepared.state(&potential, &pressure)?;
    eprintln!(
        "time_s,kinetic_j,compression_j,input_work_j,boundary_dissipation_j,balance_residual_j"
    );
    let energy = prepared.energy(&state)?;
    eprintln!("0,{},{},0,0,0", energy.kinetic_j, energy.compression_j);
    for index in 1..=steps {
        let step = prepared.step(&state, &[], &mut || false)?;
        let energy = prepared.energy(&step.state)?;
        eprintln!(
            "{},{},{},{},{},{}",
            index as f64 * dt,
            energy.kinetic_j,
            energy.compression_j,
            step.work.input_work_j,
            step.work.boundary_dissipation_j,
            step.work.balance_residual_j()
        );
        state = step.state;
    }
    println!("x_m,y_m,z_m,pressure_pa,pressure_potential_pa_s");
    for (vertex, point) in positions.iter().enumerate() {
        println!(
            "{},{},{},{},{}",
            point[0],
            point[1],
            point[2],
            state.pressure_pa()[vertex],
            state.pressure_potential_pa_s()[vertex]
        );
    }
    Ok(())
}
