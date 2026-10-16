//! Two declared cohesive planes in series: support bond versus a notched feature.
//! Illustrative parameters only; no calibration or bulk crack-path prediction.
use std::error::Error;
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};

struct Actuator { stiffness_n_m: f64 }
impl LinearOp for Actuator {
    fn n(&self) -> usize { 2 }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = 0.0;
        y[1] = self.stiffness_n_m*x[1];
    }
}
fn main() -> Result<(),Box<dyn Error>> {
    let weak = BilinearCohesive::new(1e9,2e5,100.0,2e9,0.01)?;
    let strong = BilinearCohesive::new(1e9,1e6,2000.0,2e9,0.01)?;
    let area = 1e-4;
    let actuator = Actuator { stiffness_n_m: 1e6 };
    for (case,bond,feature,failed_index) in [
        ("support-interface-first",&weak,&strong,0),
        ("declared-feature-plane-first",&strong,&weak,1),
    ] {
        let patches = vec![
            CohesivePatch::new(1,area,bond,vec![JumpTerm { dof: 0, local: [1.0,0.0,0.0] }])?,
            CohesivePatch::new(2,area,feature,vec![
                JumpTerm { dof: 0, local: [-1.0,0.0,0.0] },
                JumpTerm { dof: 1, local: [1.0,0.0,0.0] },
            ])?,
        ];
        let system = CohesiveSystem::new(&actuator,patches,Budget {
            max_dofs: 2, max_patches: 2, max_jump_terms: 3, max_krylov_basis_entries: 256,
        })?;
        let config = Config {
            newton: NewtonKrylovConfig { linear_restart: 2, absolute_tolerance: 1e-9,
                relative_tolerance: 1e-13, forcing_maximum: 0.01, ..NewtonKrylovConfig::default() },
            max_newton_iterations: 80, max_operator_evaluations: 2000,
            force_tolerance_n: 1e-7, local_energy_tolerance_j_m2: 1e-7,
            total_cohesive_energy_tolerance_j: 1e-9, maximum_displacement_m: 0.005,
        };
        let mut state = system.bonded_state();
        let mut work = 0.0;
        println!("case={case}; area_m2={area}; columns=actuator_m,reaction_n,bond_jump_m,feature_jump_m,bond_damage,feature_damage");
        for step in 0..=120 {
            let imposed = f64::from(step)*1e-5;
            let result = state.solve_controlled(&[0.0,actuator.stiffness_n_m*imposed],
                &[[0.0;3];2],config,&mut || false)?;
            work += result.cohesive_energy_j.work;
            state = result.state;
            let force = actuator.stiffness_n_m*(imposed-state.displacements_m()[1]);
            if step % 10 == 0 {
                println!("{imposed:.8e},{force:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
                    state.points()[0].jump_m()[0],state.points()[1].jump_m()[0],
                    state.points()[0].damage(),state.points()[1].damage());
            }
            // Both joints carry the same actual force while in series.
            for point in state.points() {
                assert!((area*point.traction_pa()[0]-force).abs() < config.force_tolerance_n);
            }
        }
        assert!(state.points()[failed_index].failed());
        assert_eq!(state.points()[1-failed_index].damage(),0.0);
        assert!((work-area*weak.fracture_energy_j_m2()).abs() < 1e-8);
        println!("fractured_patch={}; cohesive_path_work_j={work:.8e}; no calibrated strength claim",failed_index+1);
    }
    Ok(())
}
