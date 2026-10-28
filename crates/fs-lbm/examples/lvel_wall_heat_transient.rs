//! Native LVEL flow/heat ramp over a conducting slab with internal and exterior
//! equilibrium thermal walls. Writes one CSV row per accepted native step.
//! The short coarse run demonstrates the coupled application path; it is not
//! a resolved turbulent flow or an empirical heat-transfer validation.
//!
//! cargo run --release -p fs-lbm --example lvel_wall_heat_transient

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, EnergyConfig, FluidProperties, FvBoundary, PressureSolver,
    SimpleConfig, SolidMaterial, ThermalFace, ThermalSetup, Turbulence,
    UnsteadyConfig, Voxel, VoxelDomain,
};
use fs_lbm::conjugate::turbulence::thermal_wall::{
    BoundaryWallHeat, WallHeatConfig, WallHeatFace, WallHeatMarchConfig,
    march_conjugate_with_wall_heat,
};
use fs_lbm::d3q19::Face3;

fn run() -> Result<(), ChtError> {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::from_fn(8, 6, 2, 0.004, |point| {
        if point[1] < 0.008 { Voxel::Solid(0) } else { Voxel::Fluid }
    })?;
    let fluid = FluidProperties::dry_air_300k();
    let mut slab = SolidMaterial::new("declared conducting slab", 16.0);
    slab.volumetric_heat_capacity_j_m3_k = Some(2.43e6);
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet { velocity: [2.0, 0.0, 0.0] };
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut flow = SimpleConfig::new(faces);
    flow.turbulence = Turbulence::Lvel;
    flow.pressure_solver = PressureSolver::IluBicgstab;
    let mut clock = UnsteadyConfig::new(0.002, 50);
    clock.inner_iterations = 500;
    clock.inner_tolerance = 1e-7;
    clock.probe = Some(domain.index(4, 3, 0));
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    thermal.faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    thermal.add_uniform_power(&domain, 3.0, |point| point[1] < 0.008);
    let mut internal_walls = Vec::new();
    let mut boundary_walls = Vec::new();
    for z in 0..2 {
        for x in 0..8 {
            internal_walls.push(WallHeatFace {
                fluid_cell: domain.index(x, 2, z), face: Face3::ALL[2],
            });
            boundary_walls.push(BoundaryWallHeat {
                wall: WallHeatFace { fluid_cell: domain.index(x, 5, z), face: Face3::ALL[3] },
                condition: ThermalFace::Convective { h: 25.0, ambient: 300.0 },
                contact_m2_k_w: 0.001,
            });
        }
    }
    let walls = WallHeatMarchConfig {
        internal_walls, boundary_walls,
        transport: WallHeatConfig {
            turbulent_prandtl: 0.9, relative_tolerance: 1e-8, max_refinements: 14,
        },
    };
    let result = march_conjugate_with_wall_heat(
        &domain, &fluid, &[slab], &thermal, &flow, &clock, &walls, None,
        &vec![300.0; domain.cell_count()],
        |time| (time / 0.02).min(1.0),
        |time| if time <= 0.06 { 1.0 } else { 0.0 },
        &EnergyConfig::default(), None, &gate,
    )?;
    println!("time_s,kinetic_energy_j,max_solid_k,wall_outflow_w,max_y_plus,max_wall_h_w_m2_k,energy_closure_j,wall_quadrature_evaluations");
    for ((momentum, energy), wall) in result.march.flow.records.iter()
        .zip(&result.march.energy_records).zip(&result.wall_records)
    {
        println!("{},{},{},{},{},{},{},{}", energy.time_s,
            momentum.kinetic_energy_j, energy.max_solid_temperature_k,
            wall.boundary_outflow_w, wall.maximum_y_plus,
            wall.maximum_heat_transfer_w_m2_k, energy.closure_j,
            wall.quadrature_evaluations);
    }
    eprintln!("integrated_exterior_wall_outflow_j={}", result.boundary_outflow_j);
    eprintln!("Wall watts are already included in the native step balances; do not add them again.");
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("transient thermal wall simulation refused: {error:?}");
        std::process::exit(1);
    }
}
