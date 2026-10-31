//! Run: cargo run -p fs-topopt --features cutfem-marquee --release \
//!   --example material_mass_sdf3 -- 24 250000 4 1.5 250 1 0.0002
//! Args: accepted updates, total Krylov iterations, B/A density ratio,
//! B/A stiffness ratio, minimum right-support reaction [N], downward acceleration
//! [m/s^2], maximum reference-volume-average downward displacement [m].
//! All geometry is m. Default acceleration is 1 m/s^2, not standard Earth gravity.
//! Synthetic illustrative properties/IDs, not calibrated catalog materials:
//! A has diagonal Mandel stiffness 1 MPa and density 500 kg/m^3; B scales these
//! independently. Use stiffness ratio 1 to isolate density's design influence.
//! The 0.66 m long slab is extended by 0.2%; reaction and self-weight deflection
//! limits constrain minimum projected mass. The same native density supplies
//! kg inventory and body force. A completed update allowance does not promise
//! convergence; these synthetic small-strain inputs are not a continuum certificate.
use std::ops::ControlFlow;
use fs_blake3::ContentHash;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticMaterialField3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::{MaterialError, elastic::acoustic::ElasticAcousticState};
use fs_topopt::{SimpParams, SolveBudget, SolveControl};
use fs_topopt::sdf3::{CutDensityStudy3, Sdf3Elasticity};
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
fn material(modulus_pa: f64, density_kg_m3: f64, illustrative_id: u8)
    -> Result<ElasticAcousticState, MaterialError> {
    let c = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { modulus_pa } else { 0.0 }));
    ElasticAcousticState::try_new_mandel(c, density_kg_m3, ContentHash([illustrative_id; 32]), 1.0)
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] { [0.002*(p[0]-0.17), 0.0, 0.0] }
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    if n[0] > 0.0 { [1.0, 0.0, 0.0] } else { [0.0; 3] }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 7 {
        return Err("usage: material_mass_sdf3 [UPDATES [TOTAL_KRYLOV [DENSITY_RATIO [STIFFNESS_RATIO [MIN_REACTION_N [ACCELERATION_M_S2 [MAX_DOWNWARD_M]]]]]]]".into());
    }
    let steps = args.first().map_or(Ok(24), |s| s.parse::<usize>())?;
    let budget = args.get(1).map_or(Ok(250000), |s| s.parse::<usize>())?;
    let density_ratio = args.get(2).map_or(Ok(4.0), |s| s.parse::<f64>())?;
    let stiffness_ratio = args.get(3).map_or(Ok(1.5), |s| s.parse::<f64>())?;
    let minimum_n = args.get(4).map_or(Ok(250.0), |s| s.parse::<f64>())?;
    let acceleration = args.get(5).map_or(Ok(1.0), |s| s.parse::<f64>())?;
    let maximum_downward = args.get(6).map_or(Ok(0.0002), |s| s.parse::<f64>())?;
    if steps > 200 || !(1..=2_000_000).contains(&budget)
        || [density_ratio, stiffness_ratio].iter().any(|v| !v.is_finite() || !(0.1..=10.0).contains(v))
        || !minimum_n.is_finite() || minimum_n <= 0.0
        || !acceleration.is_finite() || !(0.0..=10.0).contains(&acceleration)
        || !maximum_downward.is_finite() || maximum_downward <= 0.0 || maximum_downward > 0.01 {
        return Err("require updates <= 200, Krylov 1..=2000000, ratios 0.1..=10, positive reaction bound, acceleration 0..=10 m/s^2, deflection bound (0,0.01] m".into());
    }
    let a = material(1e6, 500.0, 0xa1)?;
    let b = material(1e6*stiffness_ratio, 500.0*density_ratio, 0xb2)?;
    // The source callback assigns immutable native states to reference cells;
    // stiffness, support traction and physical inventory use these SAME states.
    let field = |cell: HexCell| Ok(if cell.lo()[1] < 0.5 { a.clone() } else { b.clone() });
    let mut gp = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3 {
        max_boxes: 100_000, max_points: 500_000, ..Default::default()
    }, &mut gp)?;
    let tree = Octree3::uniform(1, 2, 64)?;
    let op = AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3])?, &tree, &Slab, ElasticMaterialField3::ByCell(&field),
        &|_| false, &|_, _| true, Default::default(), Default::default(), Default::default(), &mut geometry)?;
    let force = op.body_load(&|_| [0.0; 3], || ControlFlow::Continue(()))?;
    let reference_volume: f64 = op.volumes().iter().sum();
    // Fixed observation weights measure mean downward displacement over the
    // reference domain; the observation itself does not depend on occupancy.
    let downward = op.body_load(&|_| [0.0, 0.0, -1.0/reference_volume], || ControlFlow::Continue(()))?;
    let body = op.prepare_material_body_load(&|_| [0.0, 0.0, -acceleration], 100_000,
        || ControlFlow::Continue(()))?;
    let mut study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 1_000_000), 0.15,
        SimpParams { beta: 2.0, penal: 3.0, e_min: 0.02, eta: 0.5 }).with_material_body_load(body)?;
    let reference_mass: f64 = study.operator().reference_cell_masses_kg().ok_or("missing native mass")?.iter().sum();
    let limits = [ResponseLimit3 { response: BoundedResponse3::Reaction(&right),
        sense: ResponseLimitSense3::AtLeast, bound: minimum_n, scale: minimum_n },
        ResponseLimit3 { response: BoundedResponse3::Displacement(&downward),
            sense: ResponseLimitSense3::AtMost, bound: maximum_downward, scale: maximum_downward }];
    let cases = [ResponseLimitCase3 { force: &force, prescribed: Some(&motion), limits: &limits }];
    let mut options = MinimumMaterialOptions3 {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: reference_mass }, ..Default::default()
    };
    options.optimizer.max_dimension = 64;
    options.optimizer.max_evaluations = 2048;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { per_solve_iterations: 10_000, total_iterations: budget }, &mut poll);
    let rho = vec![0.7; study.cells()];
    let mut session = MinimumMaterialStudy3::new(&mut study, &cases, &rho, options, &mut control)?;
    let outcome = session.run(steps);
    eprintln!("illustrative material A: C_M=1e6*I Pa, rho=500 kg/m^3, id=a1 repeated; B: C_M={}*I Pa, rho={} kg/m^3, id=b2 repeated",
        1e6*stiffness_ratio, 500.0*density_ratio);
    println!("iteration,mass_kg,volume_fraction,volume_m3,scaled_objective,smooth_constraint,worst_original_residual");
    for r in session.history() {
        println!("{},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}", r.iteration,
            r.mass_kg.ok_or("mass objective lost its inventory")?, r.volume_fraction,
            r.volume_fraction*reference_volume, r.objective, r.constraint, r.worst_residual);
    }
    eprintln!("latest_optimizer={outcome:?}; cumulative_work={:?}", session.work());
    let best = session.best_feasible().ok_or("no accepted design satisfies both reaction and deflection limits")?;
    eprintln!("retained_mass_kg={:.9e}; retained_volume_m3={:.9e}; right_reaction_N={:.9e}; minimum_N={minimum_n:.9e}; original_residual={:.9e}; selected_latest={}",
        best.mass_kg.ok_or("missing mass")?, best.volume_fraction*reference_volume,
        best.measurements[0][0], best.residuals[0][0], best.rho == session.accepted().rho);
    eprintln!("mean_downward_m={:.9e}; maximum_m={maximum_downward:.9e}; original_residual={:.9e}; downward_acceleration_m_s2={acceleration}",
        best.measurements[0][1], best.residuals[0][1]);
    // Install the feasible incumbent even if a later trial exhausted its budget.
    let retained = session.into_best_feasible()?;
    let body_force = study.material_body_load().ok_or("missing prepared material body load")?
        .load(&retained.projected_rho, || ControlFlow::Continue(()))?;
    let body_external_work: f64 = body_force.iter().zip(&retained.displacements[0]).map(|(f,u)| f*u).sum();
    let downward_force: f64 = body_force.chunks_exact(3).map(|f| -f[2]).sum();
    // This is actual body-force work f_body^T u [J], excluding support lifting
    // and actuator work. Reporting it uses the retained stencil without a solve.
    eprintln!("body_external_work_J={body_external_work:.9e}; downward_body_force_N={downward_force:.9e}");
    println!("cell,raw_density,projected_density");
    for (i, (raw, projected)) in retained.rho.iter().zip(&retained.projected_rho).enumerate() {
        println!("{i},{raw:.9e},{projected:.9e}");
    }
    outcome?;
    Ok(())
}
