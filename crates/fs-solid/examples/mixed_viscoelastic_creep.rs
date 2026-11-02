//! Real traction-controlled incompressible creep with retained P2/P1 histories.
//! Run `cargo run -p fs-solid --example mixed_viscoelastic_creep`.
//! Add `-- --dynamic` to release the loaded solid with consistent P2 inertia:
//! density 1000 kg/m³, 40 backward-Euler intervals of 0.01 s. Physical relaxation
//! and numerical kinetic losses are printed separately, with explicit gates.
//! The unused bulk modulus is deliberately huge; pressure enforces the volume
//! constraint and is determined by lateral equilibrium. CSV is written to stdout.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::hyper3::mixed::MixedUpSettings;
use fs_solid::hyper3::mixed::visco::dynamics::{FiniteViscoMotionControls, MixedViscoTetMotion};
use fs_solid::hyper3::mixed::visco::{MixedViscoLawField, MixedViscoTetSolid};
use fs_solid::hyper3::visco::FiniteViscoStepControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cancel = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&cancel, arena,
            StreamKey { seed: 61, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let corners = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let (nodes, elements) = promote_p1_to_p2(&corners, &[[0, 1, 2, 3]], None, &cx)?;
        let material = ViscoHyperelastic::new(Hyperelastic::new(HyperelasticModel::MooneyRivlin {
            c10: 1e3, c01: 0.0, kappa: 1e300,
        }, 5.0)?, vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 4e3 }, tau_s: 0.2,
        }])?;
        let rollers: Vec<_> = nodes.iter().enumerate().flat_map(|(node, point)| {
            (0..3).filter_map(move |axis| (point[axis] == 0.0).then_some((3 * node + axis, 0.0)))
        }).collect();
        // On x+y+z=1: the nominal traction is (600/sqrt(3),0,0) Pa.
        // Exact straight quadratic face loads are 100 N per midside, zero at corners.
        let forces: Vec<_> = TET10_EDGES.iter().enumerate().filter(|(_, edge)| !edge.contains(&0))
            .map(|(edge, _)| (3 * elements[0][4 + edge], 100.0)).collect();
        let settings = MixedUpSettings {
            load_steps: 1, maximum_newton_iterations: 30,
            newton: NewtonKrylovConfig { absolute_tolerance: 2e-11, relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default() },
            length_scale_m: 1.0, stress_scale_pa: 6e3,
            maximum_assembly_bytes: 2 * 1024 * 1024,
            maximum_factor_bytes: 2 * 1024 * 1024,
        };
        let controls = FiniteViscoStepControls {
            equilibrium_work_tolerance_j: 2e-6, maximum_endpoint_work_defect_j: 5.0,
        };
        let mut solid = MixedViscoTetSolid::try_new(&nodes, &elements,
            MixedViscoLawField::Uniform(&material), VolumetricTreatment::ExactIncompressible,
            TetAssemblyBudget::standard(), 4, 64, &cx)?;
        if std::env::args().any(|argument| argument == "--dynamic") {
            let loaded = solid.advance(0.0, &rollers, &forces, [0.; 3], None,
                settings, controls, &cx)?;
            let mut motion = MixedViscoTetMotion::try_new(loaded.solid, &[1000.0],
                &vec![[0.; 3]; nodes.len()], 2 * 1024 * 1024, &cx)?;
            let motion_controls = FiniteViscoMotionControls {
                material: FiniteViscoStepControls {
                    maximum_endpoint_work_defect_j: 0.1, ..controls
                },
                maximum_numerical_kinetic_loss_j: 1.0,
            };
            println!("time_s,corner_x_displacement_m,corner_x_velocity_mps,stored_j,kinetic_j,interval_viscous_loss_j,interval_be_kinetic_loss_j,signed_endpoint_defect_j,min_det_f");
            for _ in 0..40 {
                // Coordinate-plane rollers remain; the loading face is released.
                let accepted = motion.advance(0.01, &rollers, &[], [0.; 3], None,
                    settings, motion_controls, &cx)?;
                println!("{:.3},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                    accepted.motion.solid().time_s(), accepted.motion.solid().displacement_m()[1][0],
                    accepted.motion.velocity_mps()[1][0], accepted.material_energy.stored_after_j,
                    accepted.inertia.kinetic_after_j, accepted.material_energy.dissipated_j,
                    accepted.inertia.numerical_loss_j, accepted.material_energy.endpoint_work_defect_j,
                    accepted.equilibrium.min_det_f);
                motion = accepted.motion;
            }
            return Ok(());
        }
        println!("time_s,axial_stretch,pressure_pa,stored_j,interval_viscous_loss_j,signed_endpoint_defect_j");
        for step in 0..=8 {
            let accepted = solid.advance(if step == 0 { 0.0 } else { 0.05 }, &rollers,
                &forces, [0.; 3], None, settings, controls, &cx)?;
            println!("{:.3},{:.9},{:.9},{:.9},{:.9},{:.9}", accepted.solid.time_s(),
                1.0 + accepted.solid.displacement_m()[1][0], accepted.solid.pressure_pa()[0],
                accepted.energy.stored_after_j, accepted.energy.dissipated_j,
                accepted.energy.endpoint_work_defect_j);
            solid = accepted.solid;
        }
        Ok(())
    })
}
