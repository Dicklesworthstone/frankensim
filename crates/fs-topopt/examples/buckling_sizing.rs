//! Fixed-geometry modulus allocation with a freshly recomputed buckling floor.
//! Run: cargo run --release -p fs-topopt --example buckling_sizing
use fs_ascent::projected_al::ProjectedAlOptions;
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::stability::sensitivity::element_scale_sensitivities;
use fs_solid::stability::{buckling_loads, reduced_pencil_with_element_scales};
use fs_topopt::buckling::{BucklingDesignOptions, BucklingDesignStudy};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mesh = Mesh2::quads(3.0, 0.7, 3, 1);
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: 1000.0,
        poisson: 0.28,
        plane: PlaneKind::Stress,
        formulation: Formulation::BBar,
        body_force: Some(&|_, _| [0.02, -0.03]),
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, y| [-1.2 - 0.6 * y, 0.08])],
        symmetry: vec![],
    };
    let initial = [1.3; 3];
    let baseline = element_scale_sensitivities(&problem, &initial, 1, 600)?;
    let floor = 0.8 * baseline.buckling().loads[0];
    let options = BucklingDesignOptions {
        minimum_load_factor: floor,
        spectral_steps: 600,
        optimizer: ProjectedAlOptions {
            max_evaluations: 1024,
            ..Default::default()
        },
    };
    let mut study = BucklingDesignStudy::new(
        &problem,
        &initial,
        &[0.3; 3],
        &[2.0; 3],
        &[1.0 / 3.0; 3],
        options,
    )?;
    let baseline_cost = study.accepted().sample().objective;
    let report = study.run(160)?;
    println!(
        "Discrete B-bar stiffness allocation; modulus cost has no mass or area interpretation."
    );
    println!("initial_cost={baseline_cost:.9e} required_load_factor={floor:.9e}");
    println!(
        "stop={:?} accepted_cost={:.9e} load_factor={:.9e}",
        report.stop,
        report.objective,
        study.accepted().physics().buckling().loads[0]
    );
    println!(
        "normalized_constraint={:.9e} kkt={:?} optimizer_work={:?}",
        report.constraint, report.kkt, report.work
    );
    if let Some(best) = study.best_feasible() {
        let scales = best.physics().element_scales();
        // Independent fresh equilibrium/pencil evaluation of the retained
        // physically feasible design; no original prestress is frozen.
        let (k, g, map, _) = reduced_pencil_with_element_scales(&problem, scales)?;
        let fresh = buckling_loads(&k, &g, &map, 2, 600)?;
        println!(
            "best_feasible_cost={:.9e} rechecked_factor={:.9e} meets_floor={}",
            best.sample().objective,
            fresh.loads[0],
            fresh.loads[0] >= floor
        );
        println!("element_modulus_multipliers={scales:?}");
    }
    println!(
        "Residual-qualified discrete first-mode screen; no continuum or global-optimum claim."
    );
    Ok(())
}
