use super::*;
use crate::adjoint::robin::{RobinLinearization, RobinResponse};
use crate::assemble::assemble_operator;
use crate::{ConductionProblem, ConductivityModel, ConductivityTable, InitialGuess,
    InterfaceResistance, InterfaceSurface, SolveConfig, ThermalBoundaryBuilder,
    ThermalInterfaces, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use crate::fixtures::{box_grid, on_box_face};
use crate::transient::{VolumetricHeatCapacity, backward_euler::{BackwardEuler,
    NonlinearStepConfig, StepConfig}};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel};
use fs_rep_mesh::TetComplex;

fn context<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 61, kernel_id: 849, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() < tolerance * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn skew_tet() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0,1,2,3]]),
        vec![[0.0,0.0,0.0],[1.2,0.1,0.0],[0.2,1.1,0.0],[0.1,0.2,0.9]]).unwrap()
}
fn boundary(mesh: &ConductionMesh, q: [f64;3], normal: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("fixed", |f| f.centroid[2] == 0.0, ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .region("vector", |f| f.centroid[2] != 0.0 && f.vertices.contains(&0),
            ThermalBc::neumann_vector(q).unwrap()).unwrap()
        .remainder("normal", ThermalBc::neumann(normal).unwrap()).unwrap().finish().unwrap()
}

#[test]
fn scalar_and_vector_controls_follow_original_skew_faces_and_native_loads() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = skew_tet();
        let original = boundary(&mesh, [3.0,-4.0,2.0], 0.7);
        let controls = PrescribedFluxControls::new(cx, &mesh, &original,
            &["vector", "normal"], mesh.boundary().len()).unwrap();
        let direction = [0.3,-0.7,0.2];
        let load = controls.load_direction(cx, &[
            ThermalBc::neumann_vector(direction).unwrap(), ThermalBc::neumann(-0.4).unwrap()]).unwrap();
        let material = ConductivityModel::isotropic_declared(7.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let eps = 1e-3;
        let mut systems = Vec::new();
        for sign in [-1.0,1.0] {
            let row = boundary(&mesh, std::array::from_fn(|i| [3.0,-4.0,2.0][i] + sign*eps*direction[i]),
                0.7 - sign*eps*0.4);
            systems.push(assemble_operator(cx, &mesh, &row, &material, &source, &[300.0;4]).unwrap());
        }
        for (v, value) in load.iter().enumerate() {
            let expected = if original.dirichlet().iter().any(|&(fixed,_)| fixed == v) { 0.0 }
                else { (systems[1].load[v]-systems[0].load[v])/(2.0*eps) };
            close(*value, expected, 1e-12);
        }
        let lambda = [0.0,0.0,0.0,1.7];
        let gradient = controls.pullback(cx, &lambda).unwrap();
        let PrescribedFluxGradient::WorldVector(vector) = gradient[0] else { panic!("vector coordinate") };
        let PrescribedFluxGradient::Normal(normal) = gradient[1] else { panic!("normal coordinate") };
        close(dot(&lambda, &load), dot(&vector, &direction) - 0.4*normal, 1e-14);

        // Independently authored scalar face controls retain each normal.
        // Pulling them back through q·n must reproduce the vector derivative.
        let mut builder = ThermalBoundaryBuilder::new(&mesh)
            .region("fixed", |f| f.centroid[2] == 0.0, ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        let mut selected = Vec::new();
        let mut normals = Vec::new();
        for face in mesh.boundary().iter().filter(|f| f.centroid[2] != 0.0 && f.vertices.contains(&0)) {
            let name = format!("face-{}", face.face);
            builder = builder.region(&name, |f| f.face == face.face,
                ThermalBc::neumann(dot(&[3.0,-4.0,2.0], &face.outward_normal)).unwrap()).unwrap();
            selected.push(name); normals.push(face.outward_normal);
        }
        let scalar = builder.remainder("normal", ThermalBc::neumann(0.7).unwrap()).unwrap().finish().unwrap();
        let scalar_controls = PrescribedFluxControls::new(cx, &mesh, &scalar,
            &selected.iter().map(String::as_str).collect::<Vec<_>>(), mesh.boundary().len()).unwrap();
        let names: Vec<_> = selected.iter().map(String::as_str).collect();
        let previous = crate::adjoint::MeanFluxFeedback::new(cx, &mesh, &scalar, &names,
            &vec![0.0; names.len()*names.len()], 2*mesh.vertex_count()*names.len()).unwrap();
        let previous = previous.flux_pullback(cx, &lambda).unwrap();
        let scalar_gradient = scalar_controls.pullback(cx, &lambda).unwrap();
        for (old, gradient) in previous.iter().zip(&scalar_gradient) {
            let PrescribedFluxGradient::Normal(value) = gradient else { panic!("normal coordinate") };
            close(*value, -*old, 1e-14);
        }
        assert_ne!(normals[0], normals[1]);
        let mut expected = [0.0;3];
        for (gradient, normal) in scalar_gradient.iter().zip(normals) {
            let PrescribedFluxGradient::Normal(value) = gradient else { panic!("scalar coordinate") };
            for i in 0..3 { expected[i] += value * normal[i]; }
        }
        for i in 0..3 { close(vector[i], expected[i], 1e-14); }
    });
}

fn with_contact<R>(q: [f64;3], run: impl FnOnce(ConductionProblem<'_>, &ThermalInterfaces) -> R) -> R {
    let mut reference = Vec::new(); let mut tets = Vec::new();
    for side in 0..2 {
        let (cells, points) = box_grid([2,1,1], [1.0;3]);
        let offset = reference.len() as u32;
        tets.extend(cells.tets.into_iter().map(|t| t.map(|v| v+offset)));
        reference.extend(points.into_iter().map(|[x,y,z]| [x+side as f64,y,z]));
    }
    let positions = reference.iter().map(|&[x,y,z]| [x+0.2*y+0.1*z,0.1*x+y+0.15*z,0.05*x+0.1*y+z]).collect();
    let mesh = ConductionMesh::new(TetComplex::from_tets(reference.len(),tets),positions).unwrap();
    let coord = |f: &crate::BoundaryFace, axis| f.vertices.iter().map(|&v| reference[v as usize][axis]).sum::<f64>()/3.0;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("flux", |f| on_box_face(coord(f,0),0.0)
            || (on_box_face(coord(f,1),0.0) && coord(f,0)<1.0), ThermalBc::neumann_vector(q).unwrap()).unwrap()
        .region("cold", |f| on_box_face(coord(f,0),2.0), ThermalBc::robin(30.0,300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP,RD),
        value: PropertyValue::Scalar { value: 0.1, dims: RD }, validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic vector-flux adjoint contact".into(),
            license: "CC0-1.0".into(), artifact: None } }).unwrap();
    let side = |name: &str| SurfaceSpec { material: MaterialStateId { chemistry: name.into(),
        phase: "solid".into(), process: "declared".into(), revision: 0 }, texture_frame: "world".into() };
    let card = InterfaceSystemCard::assemble(side("a"),side("b"),SystemContext {
        medium: "declared".into(), third_body: None, environment: "declared".into(), history: "declared".into(),
    },claims,vec![]).unwrap();
    let resistance = InterfaceResistance::from_card("bond",&card,&QueryPoint::new(),SelectionPolicy::SingleClaimOnly).unwrap();
    let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
    assert!(!pairs.is_empty());
    let interfaces = ThermalInterfaces::new(&mesh,&boundary,
        vec![InterfaceSurface::new("bond",pairs,resistance).unwrap()]).unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0,2.0),(450.0,12.0)]).unwrap());
    let source = ScalarField::Uniform(2.0);
    run(ConductionProblem { mesh:&mesh,boundary:&boundary,material:&material,
        element_materials:None,source:&source },&interfaces)
}
fn config() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(300.0);
    config.linear.tolerance = 1e-11;
    config.stop.residual_rtol = 1e-12;
    config.stop.step_atol = 0.0;
    config
}
fn objective(cx: &Cx<'_>, q: [f64;3], transient: bool, differentiated: bool) -> (f64,[f64;3]) {
    with_contact(q, |problem, interfaces| {
        let controls = PrescribedFluxControls::new(cx,problem.mesh,problem.boundary,
            &["flux"],problem.mesh.boundary().len()).unwrap();
        let weights: Vec<_> = (0..problem.mesh.vertex_count()).map(|i| (1.0+0.1*(i%3) as f64)/problem.mesh.vertex_count() as f64).collect();
        let evaluate = |response: &RobinResponse| {
            let field = response.temperature().to_vec();
            let value = dot(&weights,&field);
            if !differentiated { return (value,[0.0;3]); }
            assert!(response.uses_nonlinear_jacobian());
            let dual = response.pullback(cx,&weights,&[0.0],&[0.0]).unwrap();
            let PrescribedFluxGradient::WorldVector(gradient) = controls.pullback(cx,&dual.nodal_load).unwrap()[0]
                else { panic!("native vector control") };
            let mut d = response.zero_direction();
            let direction = [0.2,-0.3,0.4];
            d.nodal_load_w = controls.load_direction(cx,&[ThermalBc::neumann_vector(direction).unwrap()]).unwrap();
            let tangent = response.apply(cx,&d).unwrap();
            close(dot(&weights,&tangent.temperature_k),dot(&gradient,&direction),2e-9);
            assert_eq!(field,response.temperature(),"derivatives retain the original accepted field");
            (value,gradient)
        };
        if transient {
            let owner = BackwardEuler::uniform(cx,problem.mesh,VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
            let old = vec![300.0;problem.mesh.vertex_count()];
            let response = owner.linearize_step(cx,problem,Some(interfaces),&old,0.03,
                StepConfig { linear:config().linear,energy_tolerance_j:1e-8 },
                Some(NonlinearStepConfig { residual_rtol:1e-12,residual_atol_j:1e-12,
                    ..NonlinearStepConfig::default() }),&["cold"]).unwrap();
            assert!(response.primal().energy_residual_j.abs()<1e-8);
            assert_eq!(old,vec![300.0;problem.mesh.vertex_count()]);
            evaluate(&response)
        } else {
            let response = RobinLinearization::new_with_interfaces(cx,problem,interfaces,config(),&["cold"]).unwrap();
            evaluate(&response)
        }
    })
}

#[test]
fn vector_controls_use_complete_nonlinear_contact_steady_and_transient_adjoints() {
    context(&CancelGate::new_clock_free(), |cx| {
        let q = [4.0,-1.0,0.8];
        for transient in [false,true] {
            let (_,gradient) = objective(cx,q,transient,true);
            for axis in 0..3 {
                let eps = 0.001;
                let mut low=q; low[axis]-=eps;
                let mut high=q; high[axis]+=eps;
                let difference=(objective(cx,high,transient,false).0-objective(cx,low,transient,false).0)/(2.0*eps);
                close(gradient[axis],difference,3e-6);
                assert!(gradient[axis].abs()>1e-5,"skew traces must exercise every world component");
            }
        }
    });
}

#[test]
fn controls_refuse_changed_coordinates_unmasked_duals_budgets_and_cancellation() {
    let mesh=skew_tet(); let boundary=boundary(&mesh,[1.0,2.0,3.0],0.0);
    let controls=context(&CancelGate::new_clock_free(),|cx| {
        assert!(PrescribedFluxControls::new(cx,&mesh,&boundary,&["vector"],0).is_err());
        assert!(PrescribedFluxControls::new(cx,&mesh,&boundary,&["vector","vector"],4).is_err());
        assert!(PrescribedFluxControls::new(cx,&mesh,&boundary,&["fixed"],4).is_err());
        let controls=PrescribedFluxControls::new(cx,&mesh,&boundary,&["vector"],4).unwrap();
        assert!(controls.pullback(cx,&[1.0;4]).is_err());
        assert!(controls.pullback(cx,&[0.0,0.0,0.0,f64::NAN]).is_err());
        assert!(controls.load_direction(cx,&[ThermalBc::neumann(1.0).unwrap()]).is_err());
        controls
    });
    let gate=CancelGate::new_clock_free(); gate.request();
    context(&gate,|cx| {
        assert!(matches!(controls.pullback(cx,&[0.0;4]),Err(ConductionError::Cancelled{..})));
        assert!(matches!(controls.load_direction(cx,&[ThermalBc::neumann_vector([0.0;3]).unwrap()]),Err(ConductionError::Cancelled{..})));
    });
}
