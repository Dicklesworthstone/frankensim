//! Native LVEL flow over a heated conducting slab, followed by the original
//! joint fluid/solid energy solve with explicit equilibrium thermal wall laws.
//! This small mesh demonstrates the complete path; it is not a converged
//! turbulent channel benchmark or a universal heat-transfer correlation.
//!
//! cargo run --release -p fs-lbm --example lvel_thermal_wall

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, EnergyConfig, FluidProperties, FvBoundary, SimpleConfig, SolidMaterial,
    ThermalFace, ThermalSetup, Turbulence, Voxel, VoxelDomain, simple_flow,
};
use fs_lbm::conjugate::turbulence::thermal_wall::{
    WallHeatConfig, WallHeatFace, solve_energy_with_wall_heat,
};
use fs_lbm::d3q19::Face3;

fn run() -> Result<(), ChtError> {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.004; 16], vec![0.004; 6], vec![0.004; 2]],
        |point| if point[1] < 0.004 { Voxel::Solid(0) } else { Voxel::Fluid },
    )?;
    let fluid = FluidProperties::dry_air_300k();
    let solids = [SolidMaterial::new("declared conducting slab", 16.0)];
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet { velocity: [4.0, 0.0, 0.0] };
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::wall();
    let mut flow_config = SimpleConfig::new(faces);
    flow_config.turbulence = Turbulence::Lvel;
    flow_config.max_iterations = 5000;
    let flow = simple_flow(&domain, &fluid, &flow_config, &gate)?;
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    thermal.faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    thermal.add_uniform_power(&domain, 3.0, |point| point[1] < 0.004);
    let walls: Vec<_> = (0..2).flat_map(|z| (0..16).map(move |x| (x, z)))
        .map(|(x, z)| WallHeatFace { fluid_cell: domain.index(x, 1, z), face: Face3::ALL[2] })
        .collect();
    let config = WallHeatConfig {
        turbulent_prandtl: 0.9, relative_tolerance: 1e-8, max_refinements: 14,
    };
    let result = solve_energy_with_wall_heat(&domain, &fluid, &solids, &flow, &thermal,
        &walls, config, &EnergyConfig::default(), &gate)?;
    let max_solid = result.energy.max_where(|cell| !domain.is_fluid(cell))
        .expect("the explicitly constructed slab contains solid cells");
    println!("maximum_solid_temperature_k,{}", max_solid.1);
    println!("solid_to_fluid_heat_w,{}", result.energy.solid_to_fluid_heat_w(&domain));
    println!("boundary_heat_out_w,{}", result.energy.report.balance.boundary_outflow_w);
    println!("energy_residual_w,{}", result.energy.report.balance.residual_w);
    println!("cell,y_plus,friction_velocity_m_s,wall_h_w_m2_k,added_resistance_m2_k_w");
    for face in &result.faces {
        println!("{},{},{},{},{}", face.wall.fluid_cell, face.transfer.y_plus,
            face.transfer.friction_velocity_m_s, face.transfer.heat_transfer_w_m2_k,
            face.correction_m2_k_w);
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("thermal wall simulation refused: {error:?}");
        std::process::exit(1);
    }
}
