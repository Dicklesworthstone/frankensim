//! Existing example's minimum-volume/stress-cap mode; same implicit geometry.
use super::*;
use fs_topopt::sdf3::design::{StressDesignOptions3, StressDesignStudy3};

pub fn run(updates: usize, iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll)?;
    let op = CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3])?,
        [4, 2, 2],
        &CurvedCantilever,
        &IsotropicElastic::new(1.0, 0.3, 1.0)?,
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut quadrature,
    )?;
    let y = op.body_load(&|_| [0.0, -1.0, 0.0], || ControlFlow::Continue(()))?;
    let z = op.body_load(&|_| [0.0, 0.0, -1.0], || ControlFlow::Continue(()))?;
    let mut study = CutDensityStudy3::new(op, 0.15, SimpParams::default());
    let loads = [
        LoadCase {
            force: &y,
            weight: 0.3,
        },
        LoadCase {
            force: &z,
            weight: 0.7,
        },
    ];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: iterations,
            ..Default::default()
        },
        &mut poll,
    );
    // The named fixture chooses the cap from a fully solved uniform reference;
    // that same cap and fixed loads remain unchanged throughout optimization.
    let reference = study.evaluate_stress(
        &vec![0.55; study.cells()],
        &loads,
        Default::default(),
        &mut control,
    )?;
    let options = StressDesignOptions3 {
        stress_limit: reference.aggregate,
        density_floor: 0.05,
        ..Default::default()
    };
    let rho = vec![0.75; study.cells()];
    let mut design = StressDesignStudy3::new(&mut study, &loads, &rho, options, &mut control)?;
    let result = design.run(updates);
    println!(
        "iteration,volume_fraction,stress_aggregate,sampled_relaxed_max,sampled_physical_max,constraint_violation,feasible"
    );
    for row in design.history() {
        println!(
            "{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}",
            row.iteration,
            row.volume_fraction,
            row.stress_aggregate,
            row.sampled_relaxed_max,
            row.sampled_physical_max,
            row.constraint_violation,
            row.feasible
        );
    }
    eprintln!(
        "stress_limit={:.17e}; measure=normalized-volume-and-load-weighted-qp-von-mises; q={}; p={}; sampled_maximum_cap_claimed=false; reference_volume={:.17e}",
        options.stress_limit,
        options.stress.relaxation_power,
        options.stress.aggregation_power,
        reference.volume_fraction
    );
    match &result {
        Ok(report) => eprintln!(
            "stress_stop={:?}; kkt={:?}; multiplier={:.17e}; penalty={:.17e}; accepted_feasible={}; optimizer_work={:?}; solve_work={:?}",
            report.stop,
            report.kkt,
            report.multiplier,
            report.penalty,
            design.feasible(),
            report.work,
            design.work()
        ),
        Err(error) => eprintln!(
            "stress_stop={error}; accepted_feasible={}; optimizer_work={:?}; solve_work={:?}",
            design.feasible(),
            design.optimizer_work(),
            design.work()
        ),
    }
    let retained = design.best_feasible().unwrap_or_else(|| design.accepted());
    eprintln!(
        "exported_design={}; exported_feasible={}; volume_fraction={:.17e}; aggregate={:.17e}; sampled_relaxed_max={:.17e}; sampled_physical_max={:.17e}; geometry_rebuilt=false; continuum_stress_certified=false; global_optimum_certified=false",
        if design.best_feasible().is_some() {
            "best-feasible"
        } else {
            "last-accepted"
        },
        retained.aggregate / options.stress_limit - 1.0 <= options.optimizer.tolerance,
        retained.volume_fraction,
        retained.aggregate,
        retained.sampled_relaxed_max,
        retained.sampled_physical_max
    );
    println!(
        "cell_x,cell_y,cell_z,raw_density,projected_density,relaxed_max_case_y,relaxed_max_case_z"
    );
    for (i, key) in design.study().operator().cell_keys().iter().enumerate() {
        println!(
            "{},{},{},{:.17e},{:.17e},{:.17e},{:.17e}",
            key[0],
            key[1],
            key[2],
            retained.rho[i],
            retained.projected_rho[i],
            retained.cell_relaxed_max[0][i],
            retained.cell_relaxed_max[1][i]
        );
    }
    println!("node_x,node_y,node_z,uy_load_x,uy_load_y,uy_load_z,uz_load_x,uz_load_y,uz_load_z");
    for (i, p) in design.study().operator().nodes().iter().enumerate() {
        println!(
            "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0],
            p[1],
            p[2],
            retained.displacements[0][3 * i],
            retained.displacements[0][3 * i + 1],
            retained.displacements[0][3 * i + 2],
            retained.displacements[1][3 * i],
            retained.displacements[1][3 * i + 1],
            retained.displacements[1][3 * i + 2]
        );
    }
    result?;
    Ok(())
}
