//! Actual regional-mean producer with independent continuum material admission.
#![cfg(feature = "thermal-verification")]
mod support;
use fs_conduction::*;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::verification::{ConductionBoundError, FluxBudget, TetError};
use fs_conduction::verification::region::{
    GoalResidualLimits, RegionMeanConfig, bound_temperature_region_mean,
};
use fs_conduction::verification::region::material_range::
    bound_temperature_region_mean_with_slab_validity as bound;
use support::{with_cx, with_cancelled_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([4, 2, 2], [0.1, 0.02, 0.03]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn boundary(mesh: &ConductionMesh, reverse: bool) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("cooling", |face| on_box_face(face.centroid[0], if reverse { 0.0 } else { 0.1 }),
            ThermalBc::robin(20.0, 300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}
fn material(low: f64, high: f64, k: f64) -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(low, k), (high, k)]).unwrap())
}
fn config(mesh: &ConductionMesh) -> RegionMeanConfig {
    RegionMeanConfig {
        dual: LinearConfig { tolerance: 1e-11, ..LinearConfig::default() },
        residual_limits: GoalResidualLimits {
            max_rows: mesh.vertex_count(), max_nonzeros: mesh.vertex_count().pow(2),
        },
        flux: FluxBudget { max_cells: mesh.element_count(), max_iterations: 128 },
    }
}
fn solve_config() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(300.0);
    config.stop.step_atol = 0.0;
    config
}

#[test]
fn finite_material_range_encloses_heating_cooling_and_reversed_slab_means() {
    with_cx(|cx| {
        let mesh = mesh();
        let material = material(290.0, 310.0, 10.0);
        let original = material.clone();
        let cells = (0..mesh.element_count()).collect::<Vec<_>>();
        for reverse in [false, true] { for density in [-1000.0, 0.0, 1000.0] {
            let boundary = boundary(&mesh, reverse);
            let source = ScalarField::Uniform(density);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            let solved = fs_conduction::solve(cx, problem, solve_config()).unwrap();
            assert!(bound_temperature_region_mean(cx, problem, &[], &solved.temperature,
                &cells, config(&mesh)).is_err(), "legacy admission must not silently broaden");
            let checked = bound(cx, problem, &solved.temperature, &cells, 0, reverse, config(&mesh)).unwrap();
            let exact_mean = 300.0 + density * (0.1/20.0 + 0.1*0.1/(3.0*10.0));
            let enclosure = checked.region.bound.enclosure;
            assert!(enclosure.lo <= exact_mean && exact_mean <= enclosure.hi,
                "exact slab mean {exact_mean} outside {enclosure:?}");
            let extrema = 300.0 + density * (0.1/20.0 + 0.1*0.1/(2.0*10.0));
            let [low, high] = checked.continuum_temperature_range_k;
            assert!(low <= extrema.min(300.0) && high >= extrema.max(300.0));
            assert!(low >= 290.0 && high <= 310.0);
            assert_eq!(checked.axis, 0); assert_eq!(checked.reverse, reverse);
        }}
        assert_eq!(material, original);
        assert!(material.tensor_at(320.0).is_err(), "proof must not widen the original table");
    });
}

#[test]
fn an_in_span_nodal_field_cannot_replace_the_continuum_validity_proof() {
    with_cx(|cx| {
        let mesh = mesh(); let boundary = boundary(&mesh, false);
        let source = ScalarField::Uniform(1000.0);
        let temperature = vec![300.0; mesh.vertex_count()];
        let narrow = material(299.0, 305.4, 10.0);
        let cells = (0..mesh.element_count()).collect::<Vec<_>>();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &narrow, element_materials: None, source: &source };
        assert!(matches!(bound(cx, problem, &temperature, &cells, 0, false, config(&mesh)),
            Err(ConductionBoundError::Verification(TetError::Unsupported(_)))));
        let adequate = material(299.0, 306.0, 10.0);
        let checked = bound(cx, ConductionProblem { material: &adequate, ..problem },
            &temperature, &cells, 0, false, config(&mesh)).unwrap();
        assert!(checked.continuum_temperature_range_k[1] >= 305.5);
        // The inaccurate all-300 field is allowed as an approximation; the
        // original residual correction/majorants must still enclose the PDE mean.
        let exact = 300.0 + 1000.0*(0.1/20.0 + 0.01/30.0);
        assert!(checked.region.bound.enclosure.lo <= exact);
        assert!(checked.region.bound.enclosure.hi >= exact);
    });
}

#[test]
fn assigned_tensors_keep_tangential_heterogeneity_and_ignore_unused_fallbacks() {
    with_cx(|cx| {
        let mesh = mesh(); let boundary = boundary(&mesh, false);
        let axes = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let orthotropic = |ky, kz| ConductivityModel::orthotropic(axes,
            [10.0, ky, kz].map(|k| ConductivityTable::declared_curve(
                vec![(295.0,k),(308.0,k)]).unwrap())).unwrap();
        let labels = mesh.complex().tets.iter().map(|tet|
            if mesh.positions()[tet[0] as usize][1] < 0.01 { MaterialId(7) } else { MaterialId(9) }
        ).collect::<Vec<_>>();
        let assigned = ElementMaterials::new(MaterialTable::new([
            (MaterialId(7), orthotropic(2.0,30.0)), (MaterialId(9),orthotropic(40.0,3.0)),
            (MaterialId(99), material(0.0,1.0,1.0)), // unused table entry
        ]).unwrap(), labels).unwrap();
        let saved = assigned.clone();
        let unused = material(0.0,1.0,100.0);
        let source = ScalarField::Uniform(1000.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &unused, element_materials: Some(&assigned), source: &source };
        let solved = fs_conduction::solve(cx, problem, solve_config()).unwrap();
        let cells = (0..mesh.element_count()).collect::<Vec<_>>();
        let checked = bound(cx, problem, &solved.temperature, &cells, 0, false, config(&mesh)).unwrap();
        assert_eq!(checked.material_temperature_span, TemperatureSpan::Sampled { low:295.0, high:308.0 });
        assert!(checked.region.bound.enclosure.lo <= 305.0 + 1.0/3.0);
        assert!(checked.region.bound.enclosure.hi >= 305.0 + 1.0/3.0);
        assert_eq!(assigned, saved);
    });
}

#[test]
fn unsupported_physics_and_original_budgets_do_not_produce_a_bound() {
    let mesh = mesh(); let boundary = boundary(&mesh, false);
    let material = material(290.0,310.0,10.0);
    let source = ScalarField::Uniform(1000.0);
    let temperature = vec![300.0; mesh.vertex_count()];
    let cells = (0..mesh.element_count()).collect::<Vec<_>>();
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    with_cx(|cx| {
        assert!(bound(cx, problem, &temperature, &cells, 3, false, config(&mesh)).is_err());
        assert!(bound(cx, problem, &temperature, &cells, 1, false, config(&mesh)).is_err());
        let nonlinear = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(290.0,9.0),(310.0,11.0)]).unwrap());
        assert!(bound(cx, ConductionProblem { material:&nonlinear, ..problem },
            &temperature, &cells, 0, false, config(&mesh)).is_err());
        let mut limits = config(&mesh); limits.flux.max_cells -= 1;
        assert!(matches!(bound(cx, problem, &temperature, &cells, 0, false, limits),
            Err(ConductionBoundError::Verification(TetError::Budget))));
        let flux_boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("left", |f| on_box_face(f.centroid[0],0.0),
                ThermalBc::Neumann { outward_flux: ScalarField::Uniform(1.0) }).unwrap()
            .region("right", |f| on_box_face(f.centroid[0],0.1), ThermalBc::robin(20.0,300.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        assert!(bound(cx, ConductionProblem { boundary:&flux_boundary, ..problem },
            &temperature, &cells, 0, false, config(&mesh)).is_err());
    });
    assert!(matches!(with_cancelled_cx(|cx|
        bound(cx, problem, &temperature, &cells, 0, false, config(&mesh))),
        Err(ConductionBoundError::Verification(TetError::Cancelled))));
}
