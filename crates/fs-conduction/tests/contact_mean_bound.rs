#![cfg(feature = "thermal-verification")]
//! Actual native primal/dual solves and the same-contact continuum bound.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::verification::{ConductionBoundError, FluxBudget, MeanSolveConfig, TetError,
    bound_temperature_mean_with_contacts, solve_with_contact_mean_bound};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, ElementMaterials,
    InitialGuess, InterfaceResistance, InterfaceSurface, MaterialId, MaterialTable, Nonlinearity,
    ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder, ThermalInterfaces,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel};
use fs_rep_mesh::TetComplex;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 73, kernel_id: 11, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx, &gate)
    })
}
fn config() -> MeanSolveConfig {
    let mut primal = fs_conduction::SolveConfig::default();
    primal.nonlinearity = Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 8 };
    primal.stop.residual_rtol = 1e-11;
    primal.stop.residual_atol = 1e-12;
    primal.linear.tolerance = 1e-13;
    let mut dual = primal.clone(); dual.initial = InitialGuess::Uniform(0.0);
    MeanSolveConfig { primal, dual, flux: FluxBudget::default() }
}
fn mesh(n: usize) -> (ConductionMesh, usize) {
    let mut positions = Vec::new(); let mut tets = Vec::new(); let mut per = 0;
    for side in 0..2 {
        let (complex, points) = box_grid([n; 3], [1.0; 3]); per = points.len();
        let offset = positions.len() as u32;
        tets.extend(complex.tets.into_iter().map(|tet| tet.map(|v| v+offset)));
        positions.extend(points.into_iter().map(|[x,y,z]| [x+side as f64,y,z]));
    }
    (ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap(), per)
}
fn boundary(mesh: &ConductionMesh, left: Option<f64>) -> ThermalBoundary {
    let mut builder = ThermalBoundaryBuilder::new(mesh);
    if let Some(t) = left {
        builder = builder.region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(t).unwrap()).unwrap();
    }
    builder.region("cold", |f| on_box_face(f.centroid[0], 2.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}
fn materials(mesh: &ConductionMesh, per: usize) -> ElementMaterials {
    let table = MaterialTable::new([
        (MaterialId(1), ConductivityModel::constant_tensor([[2.0,0.0,0.0],[0.0,3.0,0.5],[0.0,0.5,2.0]]).unwrap()),
        (MaterialId(2), ConductivityModel::constant_tensor([[1.0,0.0,0.0],[0.0,3.0,0.5],[0.0,0.5,2.0]]).unwrap()),
    ]).unwrap();
    let ids = mesh.complex().tets.iter().map(|tet| if (tet[0] as usize) < per { MaterialId(1) } else { MaterialId(2) }).collect();
    ElementMaterials::new(table, ids).unwrap()
}
fn resistance(r: f64) -> InterfaceResistance {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),
        value: PropertyValue::Scalar { value: r, dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS },
        validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: Vec::new(),
        provenance: Provenance { source: "declared continuum-contact fixture".into(), license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId { chemistry: name.into(), phase: "solid".into(), process: "as-fixtured".into(), revision: 0 },
        texture_frame: "contact-normal".into(),
    };
    let card = InterfaceSystemCard::assemble(surface("a"), surface("b"),
        SystemContext { medium: "dry".into(), third_body: None, environment: "vacuum".into(), history: "unaged".into() },
        claims, Vec::new()).unwrap();
    InterfaceResistance::from_card("joint", &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly).unwrap()
}
fn surfaces(mesh: &ConductionMesh, r: f64) -> Vec<InterfaceSurface> {
    vec![InterfaceSurface::new("joint", ThermalInterfaces::coincident_face_pairs(mesh).unwrap(), resistance(r)).unwrap()]
}
fn contains(bound: &fs_conduction::verification::MeanTemperatureSolution, truth: f64) {
    assert!(bound.bound.enclosure.lo <= truth && truth <= bound.bound.enclosure.hi,
        "{:?} excludes {truth}", bound.bound.enclosure);
}

#[test]
fn original_native_contact_field_reports_and_analytic_mean_are_preserved() {
    with_cx(|cx, _| {
        let (mesh, per) = mesh(2); let boundary = boundary(&mesh, Some(400.0));
        let assigned = materials(&mesh, per);
        let fallback = ConductivityModel::isotropic_declared(99.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &fallback,
            element_materials: Some(&assigned), source: &source };
        for r in [0.5, 2.5] {
            let surfaces = surfaces(&mesh, r); let before = surfaces.clone();
            let interfaces = ThermalInterfaces::new(&mesh, &boundary, surfaces.clone()).unwrap();
            let reference = fs_conduction::solve_with_interfaces(cx, problem, &interfaces, config().primal).unwrap();
            let result = solve_with_contact_mean_bound(cx, problem, &surfaces, config()).unwrap();
            assert_eq!(reference.temperature, result.primal.temperature);
            assert_eq!(interfaces.fluxes(&reference.temperature).unwrap(), interfaces.fluxes(&result.primal.temperature).unwrap());
            let q = 100.0/(1.5+r); contains(&result, 350.0+q/8.0);
            assert!(result.bound.enclosure.hi-result.bound.enclosure.lo < 1e-5);
            for (i,p) in mesh.positions().iter().enumerate() {
                let expected = if i < per { 400.0-q*p[0]/2.0 } else { 300.0+q*(2.0-p[0]) };
                assert!((result.primal.temperature[i]-expected).abs() < 1e-7);
            }
            assert_eq!(surfaces, before);
        }
    });
}

#[test]
fn affine_heated_body_drains_only_through_contact_and_inexact_fields_stay_bounded() {
    with_cx(|cx, _| {
        let (mesh, per) = mesh(2); let boundary = boundary(&mesh, None);
        let assigned = materials(&mesh, per);
        let fallback = ConductivityModel::isotropic_declared(99.0).unwrap();
        let source = ScalarField::Nodal(mesh.positions().iter().enumerate()
            .map(|(i,p)| if i < per { 6.0*p[0] } else { 0.0 }).collect());
        let surfaces = surfaces(&mesh, 0.5);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &fallback,
            element_materials: Some(&assigned), source: &source };
        let result = solve_with_contact_mean_bound(cx, problem, &surfaces, config()).unwrap();
        contains(&result, 303.1875);
        let reused = bound_temperature_mean_with_contacts(cx, problem, &surfaces,
            &result.primal.temperature, config().dual, config().flux).unwrap();
        assert_eq!(reused.bound.enclosure, result.bound.enclosure);
        assert_eq!(reused.dual.temperature, result.dual.temperature);
        let unfinished = vec![300.0; mesh.vertex_count()];
        let bounded = bound_temperature_mean_with_contacts(cx, problem, &surfaces,
            &unfinished, config().dual, config().flux).unwrap();
        assert!(bounded.bound.enclosure.lo <= 303.1875 && bounded.bound.enclosure.hi >= 303.1875);
        assert!(bounded.bound.integral.residual_correction.hi > 0.0);
        assert_eq!(unfinished, vec![300.0; mesh.vertex_count()]);
    });
}

#[test]
fn mapped_face_resistances_reorder_with_their_pairs_in_both_owners() {
    with_cx(|cx, _| {
        let (mesh, per) = mesh(2); let boundary = boundary(&mesh, Some(400.0));
        let assigned = materials(&mesh, per); let fallback = ConductivityModel::isotropic_declared(99.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
        let base = resistance(0.5);
        let mut zipped: Vec<_> = pairs.into_iter().enumerate().map(|(i,pair)| {
            (pair, base.with_measured_value(if i % 2 == 0 { 0.25 } else { 2.0 },
                UncertaintyModel::Unstated, "declared heterogeneous fixture").unwrap())
        }).collect();
        let declared = |items: &Vec<(fs_conduction::InterfaceFacePair, InterfaceResistance)>|
            vec![InterfaceSurface::new_mapped("joint", items.iter().map(|p| p.0).collect(), items.iter().map(|p| p.1.clone()).collect()).unwrap()];
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &fallback,
            element_materials: Some(&assigned), source: &source };
        let a = declared(&zipped); let before = a.clone();
        let first = solve_with_contact_mean_bound(cx, problem, &a, config()).unwrap();
        zipped.reverse(); let b = declared(&zipped);
        let second = solve_with_contact_mean_bound(cx, problem, &b, config()).unwrap();
        assert_eq!(first.primal.temperature, second.primal.temperature);
        assert_eq!(first.dual.temperature, second.dual.temperature);
        assert_eq!(first.bound.enclosure, second.bound.enclosure);
        assert_eq!(a, before);
        let original = ThermalInterfaces::new(&mesh, &boundary, a.clone()).unwrap();
        assert!(original.surface_is_mapped("joint").unwrap());
        let native = fs_conduction::solve_with_interfaces(cx, problem, &original, config().primal).unwrap();
        assert_eq!(native.temperature, first.primal.temperature);
    });
}

#[test]
fn missing_contact_invalid_field_budget_and_cancellation_do_not_return_a_bound() {
    with_cx(|cx, gate| {
        let (mesh, _) = mesh(1); let boundary = boundary(&mesh, Some(400.0));
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0); let surfaces = surfaces(&mesh, 0.5);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source };
        assert!(solve_with_contact_mean_bound(cx, problem, &[], config()).is_err());
        let mut low = config(); low.flux.max_cells = 1;
        assert!(matches!(solve_with_contact_mean_bound(cx, problem, &surfaces, low),
            Err(ConductionBoundError::Verification(TetError::Budget))));
        for bad in [vec![], vec![f64::NAN;mesh.vertex_count()], vec![300.0;mesh.vertex_count()]] {
            assert!(bound_temperature_mean_with_contacts(cx, problem, &surfaces, &bad, config().dual, config().flux).is_err());
        }
        gate.request();
        assert!(matches!(solve_with_contact_mean_bound(cx, problem, &surfaces, config()),
            Err(ConductionBoundError::Verification(TetError::Cancelled))));
    });
}
