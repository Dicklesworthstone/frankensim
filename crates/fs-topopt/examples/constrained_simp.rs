//! 48 independent raw-density variables, two structural loads, physical volume.
//!
//! cargo run --release -p fs-topopt --features newton-design --example constrained_simp
//!
//! Synthetic fixed-mesh linear elasticity; no calibrated-material, continuum,
//! global-optimality, closed-box KKT or binary-topology claim.
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStudy};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsModel};
use fs_opt::reverse::physics::box_map::LogisticBox;
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_newton::SimpNewtonModel;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, SimpParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions { p[0] *= p[0]; }
    let elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let n = elasticity.cells();
    let volumes: Vec<_> = fs_feec::element_geometry(&mesh, &positions).vol_signed.iter().map(|v| v.abs()).collect();
    let mut a = vec![0.0; elasticity.n()]; let mut b = a.clone();
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 { a[3*i+2] = -0.2; b[3*i+1] = 0.13*(1.0+p[2]); }
    }
    let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }];
    let pipeline = DesignPipeline { filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams { e_min: 0.01, penal: 3.0, beta: 2.0, eta: 0.5 } };
    let physics = SimpNewtonModel::new(&pipeline, &elasticity, &loads, &volumes,
        ComplianceSolveOptions { tolerance: 1e-11, max_iterations: 4000 },
        ComplianceNewtonBudget { solve_reservations: 100_000, linear_iterations: 10_000_000 })?;
    let raw_compliance = physics.compliance("filtered-compliance")?;
    let raw_volume = physics.volume("physical-material")?;
    // Both bounded views still share exactly one native state and work owner.
    let coordinates = LogisticBox::new(&vec![0.05;n], &vec![0.95;n], Dims::NONE)?;
    let compliance = coordinates.bind(&raw_compliance)?;
    let volume = coordinates.bind(&raw_volume)?;
    let pool = TilePool::for_host(1, 73); let gate = CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 73, kernel_id: 94, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let start = coordinates.encode(&vec![0.5;n], Some(&cx))?;
        let baseline = compliance.value_gradient(&start, Some(&cx))?.value;
        let mut builder = ProblemBuilder::new();
        let var = builder.var("bounded raw-density coordinates", Manifold::Rn { dim: n as u32 }, Dims::NONE)?;
        let cn = builder.pde_residual("filtered-compliance", var, true, COMPLIANCE_DIMS)?;
        let vn = builder.pde_residual("physical-material", var, true, Dims::NONE)?;
        let scale = builder.konst(baseline, COMPLIANCE_DIMS)?;
        let normalized = builder.div(cn, scale)?;
        builder.objective(normalized, Sense::Minimize, 1.0)?;
        let allowance = builder.konst(0.5, Dims::NONE)?;
        let constraint = builder.sub(vn, allowance)?;
        builder.constraint(constraint, ConstraintKind::LeZero, "projected physical volume fraction")?;
        let problem = builder.finish();
        let oracle = ReverseProblem::new_with_physics(&problem,
            ReverseLimits { max_nodes: 256, max_scalar_slots: 8192 },
            &[PhysicsBinding { node: cn, model: &compliance }, PhysicsBinding { node: vn, model: &volume }], None)?;
        let mut study = AugmentedNewtonStudy::new(&oracle, &start, AugmentedNewtonOptions::default(), Some(&cx))?;
        let report = study.run(1e-5, 300, AugmentedNewtonBudget {
            objective_evals: 2500, hessian_products: 10_000, multiplier_updates: 100,
        }, Some(&cx))?;
        let raw = coordinates.decode(&report.x, Some(&cx))?;
        // Independently reassemble/re-solve at the candidate with the retained
        // first-order path. Do not validate only by rereading our shared cache.
        let mut reference = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
        let fresh = pipeline.multi_load_compliance_and_gradient(&mut reference, &raw, &loads);
        let (fresh_volume, _) = pipeline.volume_and_gradient(&raw, &volumes);
        let difference = (fresh.compliance/baseline-report.objective).abs();
        if !fresh.compliance.is_finite() || !fresh_volume.is_finite()
            || difference > 1e-7*(1.0+report.objective.abs())
            || (fresh_volume-0.5-report.constraints[0]).abs() > 1e-10 {
            return Err(std::io::Error::other("independent physical check disagreed with Newton report").into());
        }
        println!("stop={:?}; chart_kkt={:?}", report.stop, report.kkt);
        println!("normalized_compliance={:.12}; projected_volume={:.12}; raw_range=[{:.9},{:.9}]",
            fresh.compliance/baseline, fresh_volume,
            raw.iter().copied().fold(f64::INFINITY, f64::min),
            raw.iter().copied().fold(f64::NEG_INFINITY, f64::max));
        println!("physics_work={:?}; newton_iterations={}; hessian_products={}",
            physics.work(), report.newton_iterations, report.hessian_products);
        println!("raw_densities={raw:?}");
        println!("Numerical candidate on this fixed mesh. A budget/iteration/inner stop is not convergence.");
        Ok(())
    })
}
