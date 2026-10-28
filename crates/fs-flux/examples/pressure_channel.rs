//! Run with: cargo run -p fs-flux --example pressure_channel
//!
//! Fully developed Stokes flow on a 2 m by 1 m channel, per unit depth.
//! Kinematic viscosity is 0.4 m²/s. The kinematic pressure falls from
//! 5.3 to 3.7 m²/s², giving u_x=y(1-y) m/s and flow 1/6 m²/s.
//! Open ends prescribe TOTAL traction, including the known tangential
//! viscous stress of the developed profile, rather than prescribing velocity.

use fs_flux::{FluxParams, FluxSystem, StokesBoundary, TriMesh};
use fs_solid::Mesh2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 12, 6));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams { nu: 0.4, ..FluxParams::default() };
    let wall = |_: [f64; 2]| [0.0; 2];
    let total_traction = |point: [f64; 2], normal: [f64; 2]| {
        let pressure = 5.3 - 0.8 * point[0];
        let shear = params.nu * (1.0 - 2.0 * point[1]);
        [
            -pressure * normal[0] + shear * normal[1],
            shear * normal[0] - pressure * normal[1],
        ]
    };
    let boundaries: Vec<_> = mesh.edges.iter().enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[0].abs() > 0.5 {
                StokesBoundary::TotalTraction { edge, value: &total_traction }
            } else {
                StokesBoundary::Velocity { edge, value: &wall }
            }
        })
        .collect();
    let solution = system
        .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)?
        .solve()?;
    let outlet_flux: f64 = mesh.edges.iter().enumerate()
        .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[0] > 0.5)
        .map(|(edge, geometry)| geometry.len * solution.field.x[2 * edge])
        .sum();
    let (velocity_error, divergence) = system.velocity_error(
        &solution.field.x, &|point| [point[1] * (1.0 - point[1]), 0.0],
    );
    let pressure_mean = solution.pressure.iter().zip(&mesh.areas)
        .map(|(pressure, area)| pressure * area).sum::<f64>()
        / mesh.areas.iter().sum::<f64>();
    println!("pressure_reference={:?}", solution.pressure_reference);
    println!("outlet_flow_m2_per_s={outlet_flux:.12e}");
    println!("analytic_flow_m2_per_s={:.12e}", 1.0 / 6.0);
    println!("mean_kinematic_pressure_m2_per_s2={pressure_mean:.12e}");
    println!("velocity_l2_error={velocity_error:.12e}");
    println!("divergence_l2={divergence:.12e}");
    println!("relative_residual={:.12e}", solution.field.rel_residual);
    Ok(())
}
