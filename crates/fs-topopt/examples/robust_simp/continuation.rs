//! Three immutable physical stages with one lifetime work account.
use fs_ascent::newton_continuation::{NewtonContinuation, NewtonStage};
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions};
use fs_ascent::trust::KrylovOptions;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::PhysicsBinding;
use fs_opt::reverse::physics::box_map::LogisticBox;
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_newton::{SimpNewtonModel, SimpNewtonResources};
use fs_topopt::{DensityElasticity, DensityFilter, RobustPipeline, SimpParams};

pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions { p[0] *= p[0]; }
    let elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let n = elasticity.cells();
    let volumes: Vec<_> = fs_feec::element_geometry(&mesh, &positions).vol_signed.iter().map(|v| v.abs()).collect();
    let mut a = vec![0.0; elasticity.n()]; let mut b = a.clone();
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 { a[3 * i + 2] = -0.2; b[3 * i + 1] = 0.13 * (1.0 + p[2]); }
    }
    let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }];
    // Explicit continuation policy: penalization, sharpness, threshold offset,
    // and the stage's dimensionless numerical KKT tolerance.
    let schedule = [(1.0, 0.0, 0.0, 1e-4), (2.0, 1.0, 0.05, 3e-5), (3.0, 2.0, 0.1, 1e-5)];
    let pipelines: Vec<_> = schedule.iter().map(|&(penal, beta, eta_offset, _)| RobustPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams { e_min: 0.01, penal, beta, eta: 0.5 }, eta_offset,
    }).collect();
    let resources = SimpNewtonResources::new(ComplianceNewtonBudget {
        solve_reservations: 200_000, linear_iterations: 20_000_000,
    });
    let models: Vec<_> = pipelines.iter().map(|p| SimpNewtonModel::new_robust_with_resources(
        p, &elasticity, &loads, &volumes, ComplianceSolveOptions { tolerance: 1e-11, max_iterations: 4000 },
        resources.clone(),
    )).collect::<Result<_, _>>()?;
    let raw_c: Vec<_> = models.iter().map(|m| m.compliance("eroded-compliance")).collect::<Result<_, _>>()?;
    let raw_v: Vec<_> = models.iter().map(|m| m.volume("dilated-material")).collect::<Result<_, _>>()?;
    // The SAME coordinate chart is essential: the runner checks variable schema,
    // but cannot infer that different callbacks secretly changed coordinates.
    let chart = LogisticBox::new(&vec![0.05; n], &vec![0.95; n], Dims::NONE)?;
    let compliance: Vec<_> = raw_c.iter().map(|m| chart.bind(m)).collect::<Result<_, _>>()?;
    let volume: Vec<_> = raw_v.iter().map(|m| chart.bind(m)).collect::<Result<_, _>>()?;
    let cap = 0.6;
    let compliance_scale = 100.0; // Explicit joule normalization, fixed across stages.
    let mut definitions = Vec::new();
    for _ in &schedule {
        let mut builder = ProblemBuilder::new();
        let x = builder.var("bounded density", Manifold::Rn { dim: n as u32 }, Dims::NONE)?;
        let cn = builder.pde_residual("eroded-compliance", x, true, COMPLIANCE_DIMS)?;
        let vn = builder.pde_residual("dilated-material", x, true, Dims::NONE)?;
        let scale = builder.konst(compliance_scale, COMPLIANCE_DIMS)?;
        let objective = builder.div(cn, scale)?;
        builder.objective(objective, Sense::Minimize, 1.0)?;
        let bound = builder.konst(cap, Dims::NONE)?; let residual = builder.sub(vn, bound)?;
        builder.constraint(residual, ConstraintKind::LeZero, "dilated volume fraction")?;
        definitions.push((builder.finish(), cn, vn));
    }
    let oracles: Vec<_> = definitions.iter().enumerate().map(|(i, (p, cn, vn))|
        ReverseProblem::new_with_physics(p, ReverseLimits { max_nodes: 256, max_scalar_slots: 8192 },
            &[PhysicsBinding { node: *cn, model: &compliance[i] }, PhysicsBinding { node: *vn, model: &volume[i] }], None)
    ).collect::<Result<_, _>>()?;
    let stages: Vec<_> = oracles.iter().zip(&schedule).map(|(oracle, &(_, _, _, tolerance))| NewtonStage {
        oracle, tolerance, options: AugmentedNewtonOptions::default(),
        krylov: KrylovOptions::new(0.05, 8).expect("fixed valid Krylov policy"),
    }).collect();
    let pool = TilePool::for_host(1, 74); let gate = CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 74, kernel_id: 96, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let start = chart.encode(&vec![0.5; n], Some(&cx))?;
        let mut continuation = NewtonContinuation::new(&stages, &start, Some(&cx))?;
        let report = match continuation.run(1000, AugmentedNewtonBudget {
            objective_evals: 10_000, hessian_products: 40_000, multiplier_updates: 300,
        }, Some(&cx)) {
            Ok(report) => report,
            Err(error) => {
                if let Some(active) = continuation.active() {
                    eprintln!("retained_stage={}; accepted_point={:?}; optimizer_work={:?}; physics_work={:?}",
                        continuation.active_stage(), active.inner().snapshot().x, continuation.work(), resources.work());
                }
                return Err(error.into());
            }
        };
        // Reassemble and re-solve every published stage through the older
        // first-order path. Verification work is separate from optimization.
        let total: f64 = volumes.iter().sum();
        for (i, result) in continuation.completed_stages().iter().chain(report.solution.iter()).enumerate() {
            let raw = chart.decode(&result.x, Some(&cx))?;
            let pipeline = &pipelines[i];
            let mut reference = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
            let mut fresh_c = 0.0;
            for load in &loads {
                fresh_c += load.weight * pipeline.eroded_compliance_and_gradient(&mut reference, &raw, load.force).0;
            }
            let fields = pipeline.three_fields(&raw);
            let fraction = |r: &[f64]| r.iter().zip(&volumes).map(|(r, v)| r * v).sum::<f64>() / total;
            let (ve, vn, vd) = (fraction(&fields.eroded), fraction(&fields.nominal), fraction(&fields.dilated));
            if !fresh_c.is_finite() || !vd.is_finite()
                || (fresh_c / compliance_scale - result.objective).abs() > 1e-7 * (1.0 + result.objective.abs())
                || (vd - cap - result.constraints[0]).abs() > 1e-10 {
                return Err(std::io::Error::other("fresh stage physics disagrees with continuation").into());
            }
            println!("stage={i}; policy={:?}; normalized_compliance={:.9}; volumes=({ve:.9},{vn:.9},{vd:.9}); chart_kkt={:?}",
                schedule[i], result.objective, result.kkt);
            println!("stage={i}; raw_densities={raw:?}");
        }
        println!("stop={:?}; active_stage={}; optimizer_work={:?}; shared_physics_work={:?}",
            report.stop, report.stage, report.work, resources.work());
        println!("Only Complete means all declared stage KKT tolerances passed. No binary, boundary-KKT or continuum certificate.");
        Ok(())
    })
}
