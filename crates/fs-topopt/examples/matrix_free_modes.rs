//! Real SI steel-cantilever modes and optional per-node CSV export.
//! cargo run -p fs-topopt --example matrix_free_modes -- 4 modes.csv
use std::error::Error;
use std::io::{BufWriter, Write};
use std::ops::ControlFlow;
use fs_topopt::{DensityElasticity, SolveBudget, SolveControl};
use fs_topopt::modal::{MatrixFreeEigenOptions, controlled_matrix_free_eigenpairs};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let subdivisions = args.first().map(|value| value.parse::<usize>()).transpose()?.unwrap_or(4);
    if args.len() > 2 || !(1..=12).contains(&subdivisions) {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput,
            "usage: matrix_free_modes [subdivisions: 1..12] [new-mode-field.csv]").into());
    }
    let (mesh, positions) = fs_feec::kuhn_cube(subdivisions);
    // Unit cube in metres, E = 210 GPa, nu = 0.3, rho = 7800 kg/m^3.
    let elasticity = DensityElasticity::new(&mesh, &positions, 210e9, 0.3, &|p| p[0] == 0.0);
    let mass = vec![7800.0; elasticity.cells()];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {
        per_solve_iterations: 10_000, total_iterations: 500_000,
    }, &mut callback);
    let options = MatrixFreeEigenOptions { count: 6, oversampling: 6, ..Default::default() };
    let report = controlled_matrix_free_eigenpairs(&elasticity, &mass, options, &mut control)?;
    println!("mode,omega_rad_per_s,frequency_hz,relative_residual");
    for (mode, (&lambda, &residual)) in report.values.iter().zip(&report.relative_residuals).enumerate() {
        let omega = lambda.sqrt();
        println!("{},{omega:.17e},{:.17e},{residual:.17e}", mode + 1, omega / std::f64::consts::TAU);
    }
    if let Some(path) = args.get(1) {
        // Refuse to overwrite an existing result artifact.
        let file = std::fs::OpenOptions::new().write(true).create_new(true).open(path)?;
        let mut csv = BufWriter::new(file);
        writeln!(csv, "mode,node,x_m,y_m,z_m,ux_mass_normalized,uy_mass_normalized,uz_mass_normalized")?;
        for (mode, phi) in report.modes.iter().enumerate() {
            let mode = mode + 1;
            for (node, p) in positions.iter().enumerate() {
                writeln!(csv, "{mode},{node},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                    p[0], p[1], p[2], phi[3 * node], phi[3 * node + 1], phi[3 * node + 2])?;
            }
        }
        csv.flush()?;
    }
    eprintln!("{} tetrahedra, {} free DOFs, {} outer steps, {} CG solves / {} CG iterations",
        elasticity.cells(), elasticity.free().iter().filter(|free| **free).count(),
        report.iterations, report.work.linear_solves, report.work.linear_iterations);
    Ok(())
}
