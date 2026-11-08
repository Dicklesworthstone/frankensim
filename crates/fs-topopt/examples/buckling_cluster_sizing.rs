//! Smooth stiffness allocation through repeated continuum buckling modes.
//! Run: cargo run --release -p fs-topopt --example buckling_cluster_sizing
use fs_ascent::projected_al::ProjectedAlOptions;
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::stability::{buckling_loads, reduced_pencil_with_element_scales};
use fs_topopt::buckling::{
    AggregateBucklingDesignOptions, AggregateBucklingDesignStudy, BucklingDesignOptions,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let second = mesh.clone();
    let offset = mesh.nodes.len();
    mesh.nodes
        .extend(second.nodes.iter().map(|p| [p[0], p[1] + 2.0]));
    mesh.elems.extend(
        second
            .elems
            .iter()
            .map(|e| e.iter().map(|n| n + offset).collect()),
    );
    for ((_, edges), (_, extra)) in mesh.patches.iter_mut().zip(&second.patches) {
        edges.extend(extra.iter().map(|(a, b)| (a + offset, b + offset)));
    }
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: 1000.0,
        poisson: 0.28,
        plane: PlaneKind::Stress,
        formulation: Formulation::BBar,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-1.0, 0.0])],
        symmetry: vec![],
    };
    let (k, g, map, _) = reduced_pencil_with_element_scales(&problem, &[1.0; 2])?;
    let unit = buckling_loads(&k, &g, &map, 3, 600)?.loads[0];
    let floor = 0.8 * unit;
    let sharpness = 8.0 / unit;
    let options = AggregateBucklingDesignOptions {
        design: BucklingDesignOptions {
            minimum_load_factor: floor,
            spectral_steps: 600,
            optimizer: ProjectedAlOptions {
                max_evaluations: 2048,
                tolerance: 1e-6,
                ..Default::default()
            },
        },
        modes: 2,
        sharpness,
    };
    let mut study = AggregateBucklingDesignStudy::new(
        &problem, &[1.3; 2], &[0.2; 2], &[2.0; 2], &[0.5; 2], options,
    )?;
    let initial_cost = study.accepted().sample().objective;
    let report = study.run(240)?;
    println!("Two identical force-loaded B-bar columns; complete repeated-mode packet.");
    println!("Element multipliers scale Young's modulus; cost is not a mass or section-area law.");
    println!("initial_cost={initial_cost:.9e} required_load_factor={floor:.9e}");
    println!(
        "stop={:?} cost={:.9e} ks={:.9e} first_factor={:.9e} kkt={:?}",
        report.stop,
        report.objective,
        study.accepted().physics().value(),
        study.accepted().physics().buckling().loads[0],
        report.kkt,
    );
    println!("optimizer_work={:?}", report.work);
    println!(
        "analytic_equal_scale_optimum={:.9e}",
        0.8 + 2.0_f64.ln() / 8.0
    );
    if let Some(best) = study.best_feasible() {
        let scales = best.physics().element_scales();
        // Fresh equilibrium, geometric stiffness and physical eigenvalues;
        // independently reconstruct the two-mode KS envelope.
        let (k, g, map, _) = reduced_pencil_with_element_scales(&problem, scales)?;
        let fresh = buckling_loads(&k, &g, &map, 3, 600)?;
        let first = fresh.loads[0];
        let ks = first - (1.0 + (-sharpness * (fresh.loads[1] - first)).exp()).ln() / sharpness;
        println!(
            "best_feasible_cost={:.9e} rechecked_ks={ks:.9e} rechecked_first={first:.9e} meets_floor={}",
            best.sample().objective,
            ks >= floor,
        );
        println!("element_modulus_multipliers={scales:?}");
        println!(
            "cost_reduction_percent={:.6}",
            100.0 * (1.0 - best.sample().objective / initial_cost)
        );
    }
    Ok(())
}
