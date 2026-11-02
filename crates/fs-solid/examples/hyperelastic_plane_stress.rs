//! A thin hyperelastic coupon under a dead nominal traction.
//!
//! Run: cargo run -p fs-solid --example hyperelastic_plane_stress
//! Exact volume constraint: append `-- --incompressible`.
//! Coefficients are declared model inputs, not a fit to a named material.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let incompressible = match arguments.as_slice() {
        [] => false,
        [argument] if argument == "--incompressible" => true,
        _ => return Err("usage: hyperelastic_plane_stress [--incompressible]".into()),
    };
    let mesh = Mesh2::quads(0.02, 0.01, 4, 2);
    let reference_thickness_m = 0.001;
    let material = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e6,
            lambda: 40e6,
        },
        4.0,
    )?;
    let clamp = |_: f64, _: f64| [0.0, 0.0];
    let pull = |_: f64, _: f64| [1e6, 0.0];
    let problem = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![(Patch::Left, &clamp)],
        traction: vec![(Patch::Right, &pull)],
        settings: NewtonSettings {
            tol: 1e-6,
            max_iters: 40,
            load_steps: 5,
            max_backtracks: 30,
        },
    };
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<(), SolidError> {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 49,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let local = PlaneStressSettings::default();
        let (displacement, report) = if incompressible {
            problem.solve_incompressible_plane_stress(local, &cx)?
        } else {
            problem.solve_plane_stress(local, &cx)?
        };
        let flat: Vec<f64> = displacement.iter().flatten().copied().collect();
        let thickness_stretches: Vec<f64> = if incompressible {
            problem
                .incompressible_plane_stress_states(&flat, local, &cx)?
                .iter()
                .map(|state| state.deformation_gradient[8])
                .collect()
        } else {
            problem
                .plane_stress_states(&flat, local, &cx)?
                .iter()
                .map(|state| state.deformation_gradient[8])
                .collect()
        };
        let min_thickness = thickness_stretches
            .iter()
            .map(|stretch| stretch * reference_thickness_m)
            .fold(f64::INFINITY, f64::min);
        let max_thickness = thickness_stretches
            .iter()
            .map(|stretch| stretch * reference_thickness_m)
            .fold(0.0, f64::max);
        let (residual, _) = if incompressible {
            problem.incompressible_plane_stress_residual_and_tangent(&flat, 1.0, local, &cx)?
        } else {
            problem.plane_stress_residual_and_tangent(&flat, 1.0, local, &cx)?
        };
        // The native 2-D solve uses unit REFERENCE thickness. Convert the
        // integrated reaction to this coupon's force only after the solve.
        let reaction_n = reference_thickness_m
            * mesh
                .patch_nodes(Patch::Left)
                .iter()
                .map(|&i| residual[2 * i])
                .sum::<f64>();
        let max_extension = displacement.iter().map(|u| u[0]).fold(0.0, f64::max);
        println!("Plane-stress hyperelastic coupon: 20 x 10 x 1 mm reference geometry");
        println!("Exact incompressibility: {incompressible}");
        println!("Maximum axial displacement: {max_extension:.9e} m");
        println!("Relaxed thickness: {min_thickness:.9e} .. {max_thickness:.9e} m");
        println!("Axial support reaction: {reaction_n:.9e} N");
        println!("Converged load steps: {}", report.histories.len());
        Ok(())
    })?;
    Ok(())
}
