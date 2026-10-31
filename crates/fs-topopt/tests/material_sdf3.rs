//! G0/G1/G3/G4: native material inventory through CutFEM design and constraints.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::response::limits::design::{
    MinimumMaterialObjective3, MinimumMaterialOptions3, MinimumMaterialStudy3,
};
use fs_topopt::sdf3::response::limits::{
    BoundedResponse3, ResponseLimit3, ResponseLimitCase3, ResponseLimitSense3,
};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3, Sdf3Elasticity};
use fs_topopt::{
    EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress, heaviside,
    heaviside_derivative,
};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        (p[0] - 0.17) * (p[0] - 0.83)
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x - Interval::new(0.17, 0.17)) * (x - Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X {
            Interval::new(2.0, 2.0) * Interval::new(lo[0], hi[0]) - Interval::new(1.0, 1.0)
        } else {
            Interval::new(0.0, 0.0)
        }
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn params() -> SimpParams {
    SimpParams {
        beta: 2.0,
        penal: 3.0,
        e_min: 0.02,
        eta: 0.5,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}
fn material(density: f64) -> ElasticAcousticState {
    let c = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    ElasticAcousticState::try_new_mandel(c, density, ContentHash([(density as u8).max(1); 32]), 1.0)
        .unwrap()
}
fn cartesian(densities: Option<[f64; 2]>) -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    if let Some(densities) = densities {
        let field = |cell: HexCell| Ok(material(densities[usize::from(cell.lo()[0] >= 0.5)]));
        CutElasticity3::build_with_material_field(
            domain(),
            [2; 3],
            &Top,
            ElasticMaterialField3::ByCell(&field),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut q,
        )
        .unwrap()
    } else {
        CutElasticity3::build(
            domain(),
            [2; 3],
            &Top,
            &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut q,
        )
        .unwrap()
    }
}
fn protected(op: CutElasticity3) -> CutDensityStudy3 {
    let mut regions = vec![PhysicalRegion3::Design; op.cells()];
    regions[0] = PhysicalRegion3::Solid;
    let last = regions.len() - 1;
    regions[last] = PhysicalRegion3::Void;
    CutDensityStudy3::new(op, 0.15, params())
        .with_physical_regions(regions)
        .unwrap()
}

#[test]
fn heterogeneous_mass_uses_cut_volumes_and_the_same_protected_design_pullback() {
    let mut study = protected(cartesian(Some([2.0, 8.0])));
    let other = protected(cartesian(Some([8.0, 2.0])));
    let volumes = study.operator().volumes();
    let weights: Vec<_> = study
        .operator()
        .cell_keys()
        .iter()
        .zip(&volumes)
        .map(|(key, v)| v * if key[0] == 0 { 2.0 } else { 8.0 })
        .collect();
    close(volumes.iter().sum(), 0.73, 1e-9);
    close(weights.iter().sum(), 3.65, 1e-9);
    assert_eq!(
        study.operator().reference_cell_masses_kg().unwrap(),
        weights
    );
    let scales = study.operator().scales().to_vec();
    let rho: Vec<_> = (0..study.cells())
        .map(|i| 0.28 + 0.055 * i as f64)
        .collect();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (mass, gradient) = study.mass_and_gradient(&rho, &mut control).unwrap();
    assert_eq!(control.work().linear_solves, 2);
    assert_eq!(control.work().preconditioner_operator_applications, 0);
    let volume = study.volume_and_gradient(&rho, &mut control).unwrap();
    // Swapping physical densities must not change the geometric graph filter.
    assert_eq!(
        volume,
        other.volume_and_gradient(&rho, &mut control).unwrap()
    );
    assert!((mass - other.mass_and_gradient(&rho, &mut control).unwrap().0).abs() > 1e-3);
    for i in 0..rho.len() {
        let mut plus = rho.clone();
        let mut minus = rho.clone();
        plus[i] += 1e-4;
        minus[i] -= 1e-4;
        let fd = (study.mass_and_gradient(&plus, &mut control).unwrap().0
            - study.mass_and_gradient(&minus, &mut control).unwrap().0)
            / 2e-4;
        close(fd, gradient[i], 2e-6);
    }
    // Raw controls in fixed physical cells still influence neighboring free cells.
    assert!(gradient[0] > 0.0 && gradient[rho.len() - 1] > 0.0);
    assert_eq!(study.operator().scales(), scales);
    let zero = vec![0.0; study.operator().n()];
    let evaluation = study
        .evaluate(
            &rho,
            &[LoadCase {
                force: &zero,
                weight: 1.0,
            }],
            &mut control,
        )
        .unwrap();
    assert_eq!(evaluation.projected_rho[0], 1.0);
    assert_eq!(evaluation.projected_rho[rho.len() - 1], 0.0);
    close(
        mass,
        evaluation
            .projected_rho
            .iter()
            .zip(&weights)
            .map(|(r, w)| r * w)
            .sum(),
        1e-12,
    );
    close(
        volume.0,
        evaluation
            .projected_rho
            .iter()
            .zip(&volumes)
            .map(|(r, v)| r * v)
            .sum::<f64>()
            / volumes.iter().sum::<f64>(),
        1e-12,
    );
    let installed = study.operator().scales().to_vec();
    assert_eq!(
        study.mass_and_gradient(&rho, &mut control).unwrap(),
        (mass, gradient)
    );
    assert_eq!(study.operator().scales(), installed);
}

#[test]
fn absent_density_and_cancellation_publish_no_inventory_or_operator_change() {
    let legacy = CutDensityStudy3::new(cartesian(None), 0.15, params());
    let rho = vec![0.45; legacy.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: 0,
            ..Default::default()
        },
        &mut callback,
    );
    assert!(matches!(
        legacy.mass_and_gradient(&rho, &mut control),
        Err(EvaluationStop::Breakdown {
            stage: "sdf3-material-mass-unavailable"
        })
    ));
    assert_eq!(control.work().linear_solves, 0);
    assert_eq!(control.work().linear_iterations, 0);
    let study = protected(cartesian(Some([2.0, 8.0])));
    let scales = study.operator().scales().to_vec();
    let mut stop = |p: SolveProgress| {
        if p.stage == "sdf3-filter-transpose" {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut stop);
    assert_eq!(
        study.mass_and_gradient(&rho, &mut control),
        Err(EvaluationStop::Cancelled)
    );
    assert!(control.work().linear_iterations > 0);
    assert_eq!(control.work().preconditioner_operator_applications, 0);
    assert_eq!(study.operator().scales(), scales);
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let accepted = study.mass_and_gradient(&rho, &mut control).unwrap();
    assert_eq!(
        accepted,
        study.mass_and_gradient(&rho, &mut control).unwrap()
    );
}

fn adaptive(mixed: bool, densities: [f64; 2]) -> CutDensityStudy3<AdaptiveSolveSpace3> {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = if mixed {
        tree.refined(&[*tree.leaves().iter().next().unwrap()], poll)
            .unwrap()
    } else {
        tree
    };
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    let field = |cell: HexCell| Ok(material(densities[usize::from(cell.lo()[1] >= 0.5)]));
    let op = AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
        domain(),
        &tree,
        &Slab,
        ElasticMaterialField3::ByCell(&field),
        &|_| false,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut q,
    )
    .unwrap();
    CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15, params())
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] {
    [0.02 * (p[0] - 0.17), 0.0, 0.0]
}
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    if n[0] > 0.0 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0; 3]
    }
}

#[test]
fn adaptive_inventory_survives_refinement_and_drives_real_minimum_mass_design() {
    let mut study = adaptive(false, [2.0, 8.0]);
    let fine = adaptive(true, [2.0, 8.0]);
    let weights = study
        .operator()
        .reference_cell_masses_kg()
        .unwrap()
        .to_vec();
    let volumes = study.operator().volumes();
    let total: f64 = weights.iter().sum();
    close(total, 3.3, 1e-9);
    close(
        fine.operator()
            .reference_cell_masses_kg()
            .unwrap()
            .iter()
            .sum(),
        total,
        1e-9,
    );
    let rho = vec![0.7; study.cells()];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (initial_mass, initial_gradient) = study.mass_and_gradient(&rho, &mut control).unwrap();
    let fine_mass = fine
        .mass_and_gradient(&vec![0.7; fine.cells()], &mut control)
        .unwrap()
        .0;
    close(
        initial_mass,
        total * heaviside(0.7, params().beta, params().eta),
        1e-9,
    );
    close(
        initial_gradient.iter().sum(),
        total * heaviside_derivative(0.7, params().beta, params().eta),
        1e-9,
    );
    close(fine_mass, initial_mass, 1e-9);
    assert_eq!(control.work().preconditioner_operator_applications, 0);
    let force = vec![0.0; study.operator().n()];
    let rows = [ResponseLimit3 {
        response: BoundedResponse3::Reaction(&right),
        sense: ResponseLimitSense3::AtLeast,
        bound: 0.002,
        scale: 0.002,
    }];
    let cases = [ResponseLimitCase3 {
        force: &force,
        prescribed: Some(&motion),
        limits: &rows,
    }];
    let options = MinimumMaterialOptions3 {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: total },
        ..Default::default()
    };
    let mut session =
        MinimumMaterialStudy3::new(&mut study, &cases, &rho, options, &mut control).unwrap();
    close(session.accepted().mass_kg.unwrap(), initial_mass, 1e-12);
    for (a, b) in session
        .accepted()
        .mass_gradient
        .as_ref()
        .unwrap()
        .iter()
        .zip(&initial_gradient)
    {
        close(*a, *b, 1e-11);
    }
    session.run(12).unwrap();
    let best = session.best_feasible().unwrap();
    assert!(best.feasible(options.optimizer.tolerance));
    assert!(best.mass_kg.unwrap() < initial_mass - 1e-4);
    assert!(best.measurements[0][0] >= rows[0].bound - options.optimizer.tolerance * rows[0].scale);
    // Identical stiffness and symmetric loading: a volume objective would not
    // distinguish the sides. The mass objective preferentially removes density 8.
    let mut occupancy = [0.0; 2];
    let mut measure = [0.0; 2];
    for ((r, v), w) in best.projected_rho.iter().zip(&volumes).zip(&weights) {
        let side = usize::from(w / v > 5.0);
        occupancy[side] += r * v;
        measure[side] += v;
    }
    assert!(occupancy[1] / measure[1] + 1e-4 < occupancy[0] / measure[0]);
    for row in session.history() {
        close(row.objective, row.mass_kg.unwrap() / total, 1e-13);
    }
    let installed = session.into_best_feasible().unwrap();
    close(
        study
            .mass_and_gradient(&installed.rho, &mut control)
            .unwrap()
            .0,
        installed.mass_kg.unwrap(),
        1e-12,
    );
}

#[test]
fn positive_material_mass_cannot_normalize_to_a_zero_objective() {
    use fs_ascent::projected_al::ProjectedAlError;
    use fs_topopt::sdf3::response::ResponseError3;

    let mut study = adaptive(false, [1e-16; 2]);
    let rho = vec![0.7; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (mass, _) = study.mass_and_gradient(&rho, &mut control).unwrap();
    assert!(mass > 0.0);
    assert_eq!(mass / 1e308, 0.0);

    // A zero-load observation keeps the fixture cheap while exercising the
    // public optimizer constructor and its actual material-inventory sample.
    let zero = vec![0.0; study.operator().n()];
    let rows = [ResponseLimit3 {
        response: BoundedResponse3::Displacement(&zero),
        sense: ResponseLimitSense3::AtMost,
        bound: 1.0,
        scale: 1.0,
    }];
    let cases = [ResponseLimitCase3 {
        force: &zero,
        prescribed: None,
        limits: &rows,
    }];
    let options = MinimumMaterialOptions3 {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: 1e308 },
        ..Default::default()
    };
    assert!(matches!(
        MinimumMaterialStudy3::new(&mut study, &cases, &rho, options, &mut control),
        Err(ProjectedAlError::Evaluation(ResponseError3::Invalid(
            "missing or unrepresentable minimum-material objective"
        )))
    ));
    assert_eq!(study.operator().scales(), scales);
}
