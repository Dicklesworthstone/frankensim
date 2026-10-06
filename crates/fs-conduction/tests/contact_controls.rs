//! Real matching-contact FEM solves; no reconstructed contact operator in the gradient.
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, InterfaceFacePair, InterfaceResistance, InterfaceSurface,
    ScalarField, SolveConfig, InitialGuess, ThermalBc, ThermalBoundaryBuilder,
    ThermalInterfaces, ConductionError};
use fs_conduction::adjoint::robin::RobinLinearization;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::interface::{NonmatchingOptions, NonmatchingSurface};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 59, kernel_id: 841, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn resistance() -> InterfaceResistance {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Scalar { value: 0.1, dims },
        validity: fs_evidence::ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(), provenance: Provenance {
            source: "synthetic contact-control test, not experimental data".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec { material: MaterialStateId {
        chemistry: name.into(), phase: "solid".into(), process: "fixture".into(), revision: 0,
    }, texture_frame: "declared contact normal".into() };
    let card = InterfaceSystemCard::assemble(side("a"), side("b"), SystemContext {
        medium: "dry".into(), third_body: Some("synthetic bondline".into()),
        environment: "vacuum".into(), history: "unaged".into(),
    }, claims, Vec::new()).unwrap();
    InterfaceResistance::from_card("bondline", &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly).unwrap()
}
fn mesh() -> ConductionMesh {
    let mut tets = Vec::new(); let mut positions = Vec::new();
    for block in 0..3 {
        let (cells, points) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
        let offset = positions.len() as u32;
        tets.extend(cells.tets.into_iter().map(|t| t.map(|v| v + offset)));
        positions.extend(points.into_iter().map(|[x,y,z]| [x + block as f64,y,z]));
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn config() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(310.0);
    config.linear.tolerance = 1e-12;
    config.stop.residual_rtol = 1e-13;
    config.stop.step_atol = 0.0;
    config.stop.max_iterations = 200;
    config
}

fn physical(cx: &Cx<'_>, scales: [f64; 2], nonlinear: bool, mapped: bool,
    reverse: bool, differentiated: bool) -> (f64, Vec<f64>) {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("fixed", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(330.0).unwrap()).unwrap()
        .region("cold", |f| on_box_face(f.centroid[0], 3.0), ThermalBc::robin(17.0,290.0).unwrap()).unwrap()
        .region("side", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::robin(7.0,295.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let base = resistance();
    let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
    let mut surfaces = Vec::new();
    for (i, &scale) in scales.iter().enumerate() {
        let mut faces: Vec<_> = pairs.iter().copied().filter(|p|
            on_box_face(mesh.boundary()[p.side_a].centroid[0], (i + 1) as f64)).collect();
        let mut values: Vec<_> = faces.iter().map(|p| {
            let c = mesh.boundary()[p.side_a].centroid;
            let spatial = if mapped { 0.5 + c[1] + 0.25*c[2] } else { 1.0 };
            base.with_measured_value(0.1*scale*spatial, UncertaintyModel::Unstated,
                "explicit numerical resistance map, not metrology").unwrap()
        }).collect();
        if reverse { faces.reverse(); values.reverse(); }
        let name = format!("bondline-{i}");
        surfaces.push(if mapped { InterfaceSurface::new_mapped(name, faces, values).unwrap() }
            else { InterfaceSurface::new(name, faces, values[0].clone()).unwrap() });
    }
    if reverse { surfaces.reverse(); }
    let interfaces = ThermalInterfaces::new(&mesh, &boundary, surfaces).unwrap();
    let knots = if nonlinear { vec![(250.0,3.0),(450.0,15.0)] }
        else { vec![(250.0,10.0),(450.0,10.0)] };
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(knots).unwrap());
    let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 4.0 + 8.0*p[0] + 3.0*p[1]).collect());
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source };
    let response = RobinLinearization::new_with_interfaces(cx, problem, &interfaces, config(), &["cold", "side"]).unwrap();
    let t = &response.primal().temperature;
    let value = t[mesh.vertex_count()-1];
    if !differentiated { return (value, Vec::new()); }
    let before = response.primal().clone();
    let mut weights = vec![0.0; t.len()]; weights[t.len()-1] = 1.0;
    let gradient = response.pullback(cx, &weights, &[0.0,0.0], &[0.0,0.0]).unwrap();
    for &(v, _) in boundary.dirichlet() { assert_eq!(gradient.nodal_load[v], 0.0); }
    let controls = interfaces.matching_resistance_scale_pullback(cx, t, &gradient.nodal_load, pairs.len()).unwrap();
    assert_eq!(controls.len(), 2);
    for (i, row) in controls.iter().enumerate() {
        assert_eq!(row.interface, format!("bondline-{i}"));
        assert_eq!(row.mapped, mapped); assert_eq!(row.card_identity, base.card_identity());
        assert_eq!(row.face_pairs, pairs.len()/2);
    }
    assert_eq!(&before, response.primal(), "contraction does not replace the accepted physical solution");
    (value, controls.into_iter().map(|r| r.derivative).collect())
}

#[test]
fn contact_controls_match_independent_linear_nonlinear_and_mapped_fem_resolves() {
    context(&CancelGate::new_clock_free(), |cx| {
        for nonlinear in [false,true] { for mapped in [false,true] {
            let scales = [0.8,1.3];
            let (value, actual) = physical(cx, scales, nonlinear, mapped, false, true);
            let other = physical(cx, scales, nonlinear, mapped, true, true);
            assert_eq!((value, actual.clone()), other, "declaration order cannot detach face resistances");
            let epsilon = 2e-4_f64;
            for i in 0..2 {
                let mut plus = scales; plus[i] *= epsilon.exp();
                let mut minus = scales; minus[i] *= (-epsilon).exp();
                let expected = (physical(cx,plus,nonlinear,mapped,false,false).0
                    - physical(cx,minus,nonlinear,mapped,false,false).0)/(2.0*epsilon);
                assert!(expected.abs() > 1e-5, "each named contact must influence the chosen node");
                assert!((actual[i]-expected).abs() < 3e-5*expected.abs().max(0.01),
                    "nonlinear={nonlinear} mapped={mapped} contact={i}: {} vs {expected}", actual[i]);
            }
        } }
    });
}

fn two_tets() -> ConductionMesh {
    let points = vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0],
        [0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
    ConductionMesh::new(TetComplex::from_tets(8,vec![[0,1,2,3],[4,5,6,7]]),points).unwrap()
}
#[test]
fn contact_controls_retain_nonuniform_jump_modes_and_refuse_invalid_work() {
    let mesh = two_tets();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let pair = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap()[0];
    let make = |pair| ThermalInterfaces::new(&mesh,&boundary,vec![
        InterfaceSurface::new("bondline",vec![pair],resistance()).unwrap()]).unwrap();
    let interfaces = make(pair);
    let t = [301.0,299.0,300.0,300.0,300.0,300.0,300.0,300.0];
    let z = [1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0];
    context(&CancelGate::new_clock_free(), |cx| {
        let call = |t: &[f64], z: &[f64], cap| interfaces.matching_resistance_scale_pullback(cx,t,z,cap);
        let expected = 5.0/6.0; // area 1/2, R=1/10; both mean jumps are zero.
        let actual = call(&t,&z,1).unwrap()[0].derivative;
        assert!((actual-expected).abs() < 1e-14);
        let reversed = make(InterfaceFacePair {side_a:pair.side_b,side_b:pair.side_a});
        assert!((reversed.matching_resistance_scale_pullback(cx,&t,&z,1).unwrap()[0].derivative-actual).abs()<1e-14);
        assert_eq!(call(&t.map(|v|v+1000.0),&z.map(|v|v+17.0),1).unwrap()[0].derivative,actual);
        assert_eq!(call(&t,&[0.0;8],1).unwrap()[0].derivative,0.0);
        assert!(call(&t,&z,0).is_err());
        assert!(call(&t,&z[..7],1).is_err());
        assert!(call(&t[..1],&z[..1],1).is_err());
        let mut wrong = t; wrong[0]=f64::NAN; assert!(call(&wrong,&z,1).is_err());
        let mut wrong = z; wrong[1]=f64::INFINITY; assert!(call(&t,&wrong,1).is_err());
        let surface = NonmatchingSurface::new("overlap",vec![pair.side_a],vec![pair.side_b],resistance(),
            NonmatchingOptions {plane_tolerance_m:0.0,coverage_relative_tolerance:1e-10,
                max_pair_tests:16,max_overlap_triangles:16}).unwrap();
        let delegated = ThermalInterfaces::with_nonmatching(cx,&mesh,&boundary,Vec::new(),vec![surface]).unwrap();
        assert!(delegated.matching_resistance_scale_pullback(cx,&t,&z,1).is_err(),
            "never relabel the delegated matching portion as a complete nonmatching derivative");
    });
    let gate = CancelGate::new_clock_free(); gate.request();
    context(&gate, |cx| assert!(matches!(interfaces.matching_resistance_scale_pullback(cx,&t,&z,1),
        Err(ConductionError::Cancelled {..}))));
}
