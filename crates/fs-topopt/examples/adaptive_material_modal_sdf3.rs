//! Minimum material mass with an adaptive, native 3-D frequency constraint.
//!
//! Run with:
//! cargo run -p fs-topopt --example adaptive_material_modal_sdf3 --features cutfem-marquee
//!
//! The same source density supplies inertia and the kg objective. Each probe
//! copies the accepted physical material field; each accepted refit rebuilds
//! its own filter, constrained stiffness, material mass and spectral problem.
//! Frequency and subspace changes are numerical discretization diagnostics,
//! not continuum error certificates. The vibration model has zero prestress.
use std::ops::ControlFlow;

use fs_ascent::projected_al::ProjectedAlOptions;
use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::{ElasticMaterialField3, ElasticityError3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3, QuadratureWork3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_topopt::modal::MatrixFreeEigenOptions;
use fs_topopt::sdf3::modal::Sdf3EigenfrequencyOptions3;
use fs_topopt::sdf3::modal::campaign::{
    Sdf3EigenfrequencyCampaign3, Sdf3EigenfrequencyCampaignOptions3,
    Sdf3EigenfrequencyMarking3,
};
use fs_topopt::sdf3::modal::design::{
    Sdf3EigenfrequencyDesignOptions3, Sdf3EigenfrequencyGoal3,
};
use fs_topopt::sdf3::modal::refit::Sdf3EigenfrequencyRefitOptions3;
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;
use fs_topopt::sdf3::{CutDensityStudy3, Sdf3Elasticity};
use fs_topopt::sdf3_goal::GoalRefinementError3;
use fs_topopt::{SimpParams, SolveBudget, SolveControl};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.53 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.53, 0.53)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}

fn build(
    tree: &Octree3,
    limits: QuadratureOptions3,
    spent: &mut QuadratureWork3,
    checkpoint: &mut dyn FnMut() -> ControlFlow<()>,
) -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
    // Illustrative anisotropic source in SI units, shared by every rebuild.
    let diagonal = [20.0, 17.0, 29.0, 8.0, 13.0, 11.0];
    let tensor = std::array::from_fn(|i| {
        std::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 })
    });
    let material = ElasticAcousticState::try_new_mandel(
        tensor, 3.0, ContentHash([10; 32]), 1e-3,
    ).expect("positive definite illustrative source");
    // Every initial, probe and refit build spends the same geometry quota.
    let remaining = QuadratureOptions3 {
        max_boxes: limits.max_boxes.saturating_sub(spent.boxes),
        max_points: limits.max_points.saturating_sub(spent.points),
        ..limits
    };
    let mut callback = |_| checkpoint();
    let mut quadrature = QuadratureControl3::new(remaining, &mut callback)
        .map_err(ElasticityError3::from)?;
    let result = AdaptiveElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.6, 0.9, 0.7]).expect("positive box"),
        tree,
        &Top,
        ElasticMaterialField3::Uniform(&material),
        &|p| p[0] == 0.0 || p[1] == 0.0 || p[2] == 0.0,
        Default::default(),
        &mut quadrature,
    );
    // Charge failed and cancelled builds before propagating their error.
    let work = quadrature.work();
    spent.boxes += work.boxes;
    spent.points += work.points;
    spent.field_evaluations += work.field_evaluations;
    Ok(AdaptiveSolveSpace3::jacobi(result?, 10_000_000))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let tree = Octree3::uniform(0, 3, 512)?;
    let geometry_limits = QuadratureOptions3 { depth: 1, ..Default::default() };
    let mut geometry = QuadratureWork3::default();
    let operator = build(
        &tree, geometry_limits, &mut geometry, &mut || ControlFlow::Continue(()),
    )?;
    let mass = operator.prepare_material_mass(100_000, || ControlFlow::Continue(()))?;
    let mut study = CutDensityStudy3::new(
        operator,
        0.2,
        SimpParams { e_min: 0.02, penal: 3.0, beta: 2.0, eta: 0.5 },
    );
    let frequency = Sdf3EigenfrequencyOptions3 {
        solver: MatrixFreeEigenOptions {
            count: 2,
            oversampling: 6,
            max_iterations: 192,
            relative_tolerance: 1e-8,
            linear_tolerance: 1e-11,
            seed: 17,
        },
        aggregation_beta: 30.0,
    };
    let rho = vec![0.72; study.cells()];
    let baseline = study.evaluate_eigenfrequency(&mass, &rho, frequency, &mut control)?;
    // A deliberately loose floor keeps the small example feasible through
    // the substantial discretization change from one cell to eight cells.
    let eigenvalue_floor = 0.05 * baseline.objective.aggregate;
    let goal = Sdf3EigenfrequencyGoal3::MinimumMaterial {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: baseline.physical_mass_kg },
        volume_scale: 1.0,
    };
    let options = Sdf3EigenfrequencyRefitOptions3 {
        design: Sdf3EigenfrequencyDesignOptions3 {
            frequency,
            eigenvalue_reference: eigenvalue_floor,
            density_floor: 0.3,
            optimizer: ProjectedAlOptions {
                max_spectral_step: 0.01,
                ..Default::default()
            },
            ..Default::default()
        },
        accepted_steps: 3,
        ..Default::default()
    };
    let initial = study.solve_eigenfrequency_design(&rho, goal, options, &mut control)?;
    if !initial.completed_and_rechecked() {
        return Err(format!("initial design stopped before a feasible recheck: {:?}", initial.outcome).into());
    }
    let mut campaign = Sdf3EigenfrequencyCampaign3::new(
        initial,
        tree,
        Sdf3EigenfrequencyCampaignOptions3 {
            max_refits: 1,
            max_candidate_attempts: 2,
            marking: Sdf3EigenfrequencyMarking3::ModalResidual { theta: 0.5, max_marks: 4 },
            ..Default::default()
        },
        &mut control,
        |tree, _purpose, checkpoint| build(tree, geometry_limits, &mut geometry, checkpoint),
    )?;
    let stop = campaign.run(1)?;
    let result = campaign.current();
    let before = &result.initial.as_ref().expect("solved current-grid baseline").evaluation;
    let after = &result.recheck.as_ref().expect("fresh recheck").as_ref()
        .expect("retained feasible result").evaluation;
    println!("campaign stop: {stop:?}");
    println!("accepted grids: {}; independent probes: {}", campaign.grids().len(), campaign.probes().len());
    println!("current active cells: {}", result.study.cells());
    println!(
        "current-grid initial mass: {:.8} kg; rechecked mass: {:.8} kg",
        before.physical_mass_kg, after.physical_mass_kg,
    );
    println!(
        "lowest rechecked frequency: {:.8} Hz; required floor: {:.8} Hz",
        after.objective.modal.values[0].sqrt() / std::f64::consts::TAU,
        eigenvalue_floor.sqrt() / std::f64::consts::TAU,
    );
    println!("rechecked physical volume fraction: {:.8}", after.volume_fraction);
    println!("material and spectral feasibility rechecked: {}", result.completed_and_rechecked());
    if let Some(probe) = campaign.current_probe() {
        println!(
            "final probe: eigenvalue change {:.3e}; subspace change {:.3e}; both grids feasible {}",
            probe.eigenvalue_change, probe.subspace_change, probe.physically_feasible(),
        );
    }
    drop(campaign);
    println!(
        "cumulative geometry: {} boxes; {} points; {} field evaluations",
        geometry.boxes, geometry.points, geometry.field_evaluations,
    );
    Ok(())
}
