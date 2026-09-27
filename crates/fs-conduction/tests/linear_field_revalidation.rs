//! Exercise physical publication with ordinary tetrahedral solves, not an
//! injected residual or a precomputed ConductionReport.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    InitialGuess, ScalarField, SolveConfig, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 7306, kernel_id: 73, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic,
    )));
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel, source: ScalarField,
}
impl Fixture {
    fn new() -> Self {
        let (mesh, positions) = fs_conduction::fixtures::unit_cube(1);
        let mesh = ConductionMesh::new(mesh, positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap()).unwrap().finish().unwrap();
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(6.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            element_materials: None, source: &self.source }
    }
}
fn config() -> SolveConfig {
    SolveConfig { initial: InitialGuess::Uniform(300.0), ..SolveConfig::default() }
}

#[test]
fn refreshed_field_and_boundary_match_a_fresh_native_solve() {
    let f = Fixture::new();
    with_cx(|_, cx| {
        let original = fs_conduction::solve(cx, f.problem(), config()).unwrap();
        let saved = original.clone();
        let candidate: Vec<_> = original.temperature.iter().map(|t| t + 30.0).collect();
        let (updated, boundary) = original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &candidate, &[("air", 330.0)], 1e-7, 1e-7,
        ).unwrap();
        let independent = fs_conduction::solve(cx,
            ConductionProblem { boundary: &boundary, ..f.problem() }, config()).unwrap();
        for (&a, &b) in updated.temperature.iter().zip(&independent.temperature) {
            assert!((a - b).abs() < 1e-7);
        }
        let flux = &updated.report.robin_fluxes[0];
        assert!((flux.mean_reference_temperature_k - 330.0).abs() < 1e-10);
        assert!((flux.heat_rate_w - 6.0).abs() < 1e-7);
        assert!(updated.report.final_residual <= 1e-7);
        assert!(updated.report.energy.relative_closure() <= 1e-7);
        assert_eq!(updated.report.linear, original.report.linear);
        assert_eq!(updated.report.residual_history, original.report.residual_history);
        assert_eq!(original, saved);
        assert_ne!(boundary, f.boundary);
        let again = original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &candidate, &[("air", 330.0)], 1e-7, 1e-7,
        ).unwrap();
        assert_eq!((updated, boundary), again);
    });
}

#[test]
fn physical_residual_and_energy_are_both_required() {
    let f = Fixture::new();
    with_cx(|_, cx| {
        let original = fs_conduction::solve(cx, f.problem(), config()).unwrap();
        let saved = original.clone();
        // Stale temperatures cannot inherit a newly marched reference.
        assert!(matches!(original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &original.temperature, &[("air", 330.0)], 1e-7, 1.0,
        ), Err(ConductionError::NotConverged { .. })));
        // A deliberately loose residual gate must not bypass energy closure.
        let error = original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &original.temperature, &[("air", 330.0)], 1e6, 1e-7,
        ).unwrap_err();
        assert!(matches!(error, ConductionError::Config { parameter: "corrected field energy", .. }));
        assert_eq!(original, saved);
    });
}

#[test]
fn bad_references_fields_metadata_and_cancel_do_not_change_the_original() {
    let f = Fixture::new();
    with_cx(|gate, cx| {
        let original = fs_conduction::solve(cx, f.problem(), config()).unwrap();
        for references in [vec![("missing", 300.0)], vec![("air", f64::NAN)],
            vec![("air", 300.0), ("air", 301.0)]]
        {
            assert!(original.revalidate_linear_robin_temperature(
                cx, f.problem(), None, &original.temperature, &references, 1e-7, 1e-7,
            ).is_err());
        }
        for (residual, energy) in [(-1.0, 0.0), (f64::NAN, 0.0), (1.0, -1.0), (1.0, f64::INFINITY)] {
            assert!(original.revalidate_linear_robin_temperature(
                cx, f.problem(), None, &original.temperature, &[], residual, energy,
            ).is_err());
        }
        let mut bad = original.clone(); bad.report.free_dofs += 1;
        assert!(bad.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &original.temperature, &[], 1e-7, 1e-7,
        ).is_err());
        let mut field = original.temperature.clone(); field[0] = f64::NAN;
        assert!(original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &field, &[], 1e-7, 1e-7,
        ).is_err());
        gate.request();
        assert!(matches!(original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &original.temperature, &[], 1e-7, 1e-7,
        ), Err(ConductionError::Cancelled { .. })));
    });
}

#[test]
fn fixed_nodes_and_nonlinear_material_cannot_be_relabelled_as_linear_publication() {
    let mut f = Fixture::new();
    f.boundary = ThermalBoundaryBuilder::new(&f.mesh)
        .region("anchor", |face| face.centroid[0] < 1e-9, ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .region("air", |face| face.centroid[0] > 1.0 - 1e-9, ThermalBc::robin(2.0, 300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    with_cx(|_, cx| {
        let original = fs_conduction::solve(cx, f.problem(), config()).unwrap();
        let mut changed = original.temperature.clone();
        changed[f.boundary.dirichlet()[0].0] += 1.0;
        assert!(original.revalidate_linear_robin_temperature(
            cx, f.problem(), None, &changed, &[], 1e6, 1.0,
        ).is_err());
        let material = ConductivityModel::isotropic(
            fs_conduction::material::ConductivityTable::declared_curve(vec![(290.0, 9.0), (350.0, 11.0)]).unwrap());
        assert!(original.revalidate_linear_robin_temperature(
            cx, ConductionProblem { material: &material, ..f.problem() }, None,
            &original.temperature, &[], 1e6, 1.0,
        ).is_err());
    });
}
