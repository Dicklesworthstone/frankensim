//! Maintain an actual nodal load on a finite-strain viscoelastic tetrahedron.
//! The lateral supports constrain contraction. Material numbers are an explicit
//! synthetic example, not a specimen calibration or an incompressibility claim.
//! With `--dynamic`, release an initially stretched body and resolve its actual
//! damped motion with fixed reference density and separate numerical damping.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::hyper3::HyperTetSettings;
use fs_solid::hyper3::visco::dynamics::{FiniteViscoMotionControls, FiniteViscoTetMotion};
use fs_solid::hyper3::visco::{FiniteViscoStepControls, FiniteViscoTetSolid, ViscoelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dynamic = match std::env::args().nth(1).as_deref() {
        None => false,
        Some("--dynamic") => true,
        _ => return Err("usage: viscoelastic_creep [--dynamic]".into()),
    };
    let law = ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 10_000.0,
                lambda: 20_000.0,
            },
            4.0,
        )?,
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 30_000.0 },
            tau_s: 0.15,
        }],
    )?;
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tetrahedra = [[0, 1, 2, 3]];
    let prescribed: Vec<_> = (0..12)
        .filter(|&dof| dof != 3)
        .map(|dof| (dof, 0.0))
        .collect();
    let settings = HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        linear_tolerance: 1e-12,
        ..HyperTetSettings::default()
    };
    let controls = FiniteViscoStepControls {
        equilibrium_work_tolerance_j: 2e-7,
        maximum_endpoint_work_defect_j: 100.0,
    };
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(
        |arena| -> Result<(), Box<dyn std::error::Error>> {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let mut solid = FiniteViscoTetSolid::try_new(
                &nodes,
                &tetrahedra,
                ViscoelasticLawField::Uniform(&law),
                TetAssemblyBudget::default(),
                &cx,
            )?;
            if dynamic {
                let stretched: Vec<_> = (0..12).map(|dof| (dof, if dof == 3 { 0.05 } else { 0.0 })).collect();
                solid = solid.advance(0.0, &stretched, &[], [0.0; 3], settings, controls, &cx)?.solid;
                let mut motion = FiniteViscoTetMotion::try_new(solid, &[1200.0], &[[0.0; 3]; 4], &cx)?;
                let motion_controls = FiniteViscoMotionControls { material: controls,
                    maximum_numerical_kinetic_loss_j: 2.0 };
                println!("time_s,extension_m,velocity_mps,stored_j,kinetic_j,viscous_loss_j,numerical_kinetic_loss_j,endpoint_defect_j");
                for _ in 0..40 {
                    let next = motion.advance(0.005, &prescribed, &[], [0.0; 3], settings, motion_controls, &cx)?;
                    println!("{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                        next.motion.solid().time_s(), next.motion.solid().displacement_m()[1][0],
                        next.motion.velocity_mps()[1][0], next.material_energy.stored_after_j,
                        next.inertia.kinetic_after_j, next.material_energy.dissipated_j,
                        next.inertia.numerical_loss_j, next.material_energy.endpoint_work_defect_j);
                    motion = next.motion;
                }
                return Ok(());
            }
            println!("time_s,extension_m,axial_force_n,stored_j,interval_loss_j,endpoint_defect_j");
            for step in 0..=20 {
                let dt = if step == 0 { 0.0 } else { 0.03 };
                let next = solid.advance(
                    dt,
                    &prescribed,
                    &[(3, 900.0)],
                    [0.0; 3],
                    settings,
                    controls,
                    &cx,
                )?;
                println!(
                    "{:.6},{:.9},{:.6},{:.9},{:.9},{:.9}",
                    next.solid.time_s(),
                    next.solid.displacement_m()[1][0],
                    next.solid.element_first_piola_pa()[0][0] / 6.0,
                    next.energy.stored_after_j,
                    next.energy.dissipated_j,
                    next.energy.endpoint_work_defect_j
                );
                solid = next.solid;
            }
            Ok(())
        },
    )
}
