//! Minimum material under an actuator-force window and displacement ceilings.
//! cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example minimum_material_sdf3 -- 80 250000 0.0012 0.003
//! Arguments: accepted updates, cumulative Krylov iterations, min/max reaction.
//! Two independent extension/compression experiments on a dimensionless slab.
//! E=1, nu=0, imposed axial strain +/-0.02. Bounds are declared design limits,
//! not calibrated specimen data. No continuum, unique-recovery or actuator-work claim.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams, SolveBudget, SolveControl};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::limits::*;
use fs_topopt::sdf3::response::limits::design::*;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x-Interval::new(0.17, 0.17))*(x-Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X { Interval::new(2.0, 2.0)*Interval::new(lo[0], hi[0])-Interval::new(1.0, 1.0) }
        else { Interval::new(0.0, 0.0) }
    }
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] { [0.02*(p[0]-0.17), 0.0, 0.0] }
fn compression(p: [f64; 3], n: [f64; 3]) -> [f64; 3] { motion(p, n).map(|v| -v) }
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] { if n[0] > 0.0 { [1.0, 0.0, 0.0] } else { [0.0; 3] } }
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 4 { return Err("usage: minimum_material_sdf3 [UPDATES [KRYLOV [MIN_REACTION [MAX_REACTION]]]]".into()); }
    let steps = args.first().map_or(Ok(80), |s| s.parse::<usize>())?;
    let budget = args.get(1).map_or(Ok(250000), |s| s.parse::<usize>())?;
    let lo = args.get(2).map_or(Ok(0.0012), |s| s.parse::<f64>())?;
    let hi = args.get(3).map_or(Ok(0.003), |s| s.parse::<f64>())?;
    if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || lo >= hi { return Err("require finite 0 < MIN_REACTION < MAX_REACTION".into()); }
    let mut gp = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(), &mut gp)?;
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0; 3], [1.0; 3])?,
        &Octree3::uniform(1, 4, 4096)?, &Slab, &IsotropicElastic::new(1.0, 0.0, 1.0)?, &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut geometry)?;
    let force = op.body_load(&|_| [0.0; 3], || ControlFlow::Continue(()))?;
    let q = op.body_load(&|_| [1.0, 0.0, 0.0], || ControlFlow::Continue(()))?;
    let mut study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15, SimpParams::default());
    let tension = [
        ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtLeast, bound: lo, scale: 0.002 },
        ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtMost, bound: hi, scale: 0.002 },
        ResponseLimit3 { response: BoundedResponse3::Displacement(&q), sense: ResponseLimitSense3::AtMost, bound: 0.006, scale: 0.005 },
    ];
    let compression_limits = [ResponseLimit3 { sense: ResponseLimitSense3::AtMost, bound: -lo, ..tension[0] },
        ResponseLimit3 { sense: ResponseLimitSense3::AtLeast, bound: -hi, ..tension[1] },
        ResponseLimit3 { sense: ResponseLimitSense3::AtLeast, bound: -0.006, ..tension[2] }];
    let cases = [ResponseLimitCase3 { force: &force, prescribed: Some(&motion), limits: &tension },
        ResponseLimitCase3 { force: &force, prescribed: Some(&compression), limits: &compression_limits }];
    let rho = vec![0.5; study.cells()]; let options = MinimumMaterialOptions3::default();
    let mut p = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut p);
    let mut session = MinimumMaterialStudy3::new(&mut study, &cases, &rho, options, &mut control)?;
    let initial_volume = session.accepted().volume_fraction;
    let outcome = session.run(steps);
    println!("iteration,volume_fraction,smooth_constraint,worst_original_residual");
    for r in session.history() { println!("{},{:.17e},{:.17e},{:.17e}", r.iteration, r.volume_fraction, r.constraint, r.worst_residual); }
    // A failed fit must still expose its retained complete evidence. The latest
    // AL report can be infeasible; an earlier feasible incumbent is distinct.
    let selected = session.best_feasible().unwrap_or(session.accepted());
    let feasible = selected.feasible(options.optimizer.tolerance);
    println!("case,limit,sense,value,bound,scale,original_residual");
    for (i, case) in cases.iter().enumerate() { for (j, limit) in case.limits.iter().enumerate() {
        println!("{i},{j},{:?},{:.17e},{:.17e},{:.17e},{:.17e}", limit.sense,
            selected.measurements[i][j], limit.bound, limit.scale, selected.residuals[i][j]);
    } }
    println!("cell,raw_density,projected_density");
    for (i, (raw, physical)) in selected.rho.iter().zip(&selected.projected_rho).enumerate() {
        println!("{i},{raw:.17e},{physical:.17e}");
    }
    eprintln!("initial_volume={initial_volume:.17e}; retained_volume={:.17e}; feasible={feasible}; selected_latest={}; beta={}; work={:?}; continuum_certified=false",
        selected.volume_fraction, selected.rho == session.accepted().rho, options.limits.sharpness, session.work());
    let report = outcome?;
    eprintln!("latest_al_stop={:?}; latest_al_kkt={:?}; selected_endpoint_stationarity_not_claimed=true", report.stop, report.kkt);
    let _selected = session.into_best_feasible()?;
    Ok(())
}
