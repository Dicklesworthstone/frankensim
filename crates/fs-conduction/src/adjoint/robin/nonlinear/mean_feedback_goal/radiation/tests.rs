//! Physical nonlinear FEM references; the comparisons never run a primal.
use super::*;
use crate::{ConductionMesh, ConductivityModel, ConductivityTable, ScalarField,
    ThermalBoundary, ThermalBoundaryBuilder, ThermalBc, SurfaceEmissivity,
    AmbientRadiationConfig, InitialGuess, SolveConfig, solve_with_ambient_radiation,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY};
use crate::fixtures::{box_grid, on_box_face};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, InterpolationPolicy, UncertaintyModel};

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 51, kernel_id: 821, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn patch(reservoir: f64) -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: 0.8, dims: EMISSIVITY_DIMS },
        validity: fs_evidence::ValidityDomain::unconstrained().with("T",250.0,500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(), provenance: Provenance {
            source: "synthetic radiative goal fixture; not experimental".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "gray".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, Vec::new()).unwrap();
    AmbientRadiationPatch::new("right", SurfaceEmissivity::from_card(
        "right", &card, 350.0, SelectionPolicy::SingleClaimOnly).unwrap(), reservoir).unwrap()
}
fn config() -> SolveConfig {
    let mut c = SolveConfig::default();
    c.initial = InitialGuess::Uniform(350.0);
    c.stop.step_atol = 0.0;
    c.stop.residual_atol = 0.0;
    c.stop.residual_rtol = 1e-11;
    c.linear.tolerance = 1e-12;
    c
}
fn radiation_config() -> AmbientRadiationConfig {
    AmbientRadiationConfig { max_iterations: 200, temperature_tolerance_k: 1e-10,
        balance_tolerance_w: 1e-9, balance_relative_tolerance: 0.0, relaxation: 0.7 }
}
fn mesh() -> ConductionMesh {
    let (t, p) = box_grid([3,2,2], [1.0;3]);
    ConductionMesh::new(t,p).unwrap()
}
fn material() -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0,3.0),(500.0,15.5)]).unwrap())
}
fn mean(mesh: &ConductionMesh, t: &[f64], x: f64) -> f64 {
    let mut area = 0.0; let mut sum = 0.0;
    for face in mesh.boundary().iter().filter(|f| on_box_face(f.centroid[0],x)) {
        area += face.area;
        for &v in &face.vertices { sum += face.area/3.0*t[v as usize]; }
    }
    sum/area
}
fn comparison_config() -> LinearConfig { LinearConfig { tolerance: 1e-9, ..config().linear } }

#[test]
fn radiation_goal_uses_both_physical_states_and_has_quadratic_remainder() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let material = material();
        let source = ScalarField::Nodal(mesh.positions().iter().map(|x| 10.0+20.0*x[1]).collect());
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("left", |f| on_box_face(f.centroid[0],0.0), ThermalBc::dirichlet(340.0).unwrap()).unwrap()
            .region("right", |f| on_box_face(f.centroid[0],1.0), ThermalBc::robin(5.0,300.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let n = mesh.vertex_count();
        let mut weights = vec![0.0;n]; weights[n-1] = 1.0;
        for reservoir in [275.0,400.0] {
            let patches = [patch(reservoir)];
            let physical = solve_with_ambient_radiation(cx, problem, None, &patches,
                config(), radiation_config()).unwrap();
            let reference = &physical.conduction.temperature;
            let original = reference.to_vec();
            let mut previous = None;
            for step in [1.0,0.5,0.25] {
                let approximate: Vec<_> = reference.iter().zip(mesh.positions())
                    .map(|(t,x)| t-step*x[0]*(1.0+x[1])).collect();
                let run = |weights: &[f64]| AmbientRadiationGradient::compare_goal_at(cx, problem,
                    None, comparison_config(), reference, &approximate, &["right"],
                    &[[5.0,300.0,0.0,0.0]], &[[5.0,300.0]], &[0.0], &patches, weights, 4*n);
                let report = run(&weights).unwrap();
                let expected = reference[n-1]-approximate[n-1];
                assert!((report.signed_goal_change-expected).abs() < 1e-12);
                assert!((report.signed_residual_change+report.linearization_remainder-expected).abs() < 1e-12);
                assert!((report.nodal_contributions.iter().sum::<f64>()-report.signed_residual_change).abs() < 1e-12);
                assert!(report.primal_relative_residual < comparison_config().tolerance);
                assert!(report.dual_relative_residual < comparison_config().tolerance);
                assert!(report.uses_nonlinear_jacobian);
                let remainder = report.linearization_remainder.abs();
                if let Some(old) = previous { assert!(remainder < 0.35*old && remainder > 0.15*old,
                    "halving field difference must quarter the remainder: {old:e} -> {remainder:e}"); }
                previous = Some(remainder);
                let frozen = RobinResponse::compare_mean_robin_goal_at(cx,
                    ConductionProblem { boundary: &physical.combined_boundary, ..problem }, None,
                    comparison_config(), reference, &approximate, &[], &[], &[], &weights, 0).unwrap();
                assert!((report.signed_residual_change-frozen.signed_residual_change).abs() > 1e-3*step,
                    "the fixture must detect freezing the radiation law");
                let scaled: Vec<_> = weights.iter().map(|v| -3.0*v).collect();
                let scaled = run(&scaled).unwrap();
                assert!((scaled.signed_residual_change+3.0*report.signed_residual_change).abs() < 1e-7);
                let zero = run(&vec![0.0;n]).unwrap();
                assert_eq!(zero.signed_residual_change,0.0);
                assert_eq!(zero.dual_iterations,0);
            }
            assert_eq!(reference, &original);
            assert_eq!(physical.radiation.nonlinear_radiation_out_w.is_sign_negative(), reservoir == 400.0);
        }
    });
}

fn convection(means: [f64;2]) -> [[f64;4];2] {
    let [a,b] = means.map(|t| t-300.0);
    [[5.0+0.02*a,300.0+0.10*a+0.25*b,0.02,0.0],
     [8.0+0.04*b,285.0+0.05*a+0.12*b,0.04,0.0]]
}
fn coupled_boundary(mesh: &ConductionMesh, p: [[f64;4];2]) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("left", |f| on_box_face(f.centroid[0],0.0), ThermalBc::robin(p[0][0],p[0][1]).unwrap()).unwrap()
        .region("right", |f| on_box_face(f.centroid[0],1.0), ThermalBc::robin(p[1][0],p[1][1]).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

#[test]
fn radiation_goal_keeps_cross_wall_reference_feedback_and_its_radiative_weight() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let material = material();
        let source = ScalarField::Nodal(mesh.positions().iter().map(|x| 30.0+20.0*x[1]).collect());
        let patches = [patch(290.0)];
        let mut means = [330.0,325.0];
        let mut accepted = None;
        for _ in 0..160 {
            let boundary = coupled_boundary(&mesh, convection(means));
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            let solution = solve_with_ambient_radiation(cx, problem, None, &patches,
                config(), radiation_config()).unwrap();
            let t = solution.conduction.temperature;
            let next = [mean(&mesh,&t,0.0),mean(&mesh,&t,1.0)];
            if next.iter().zip(means).all(|(a,b)| (a-b).abs() < 1e-10) {
                accepted = Some((boundary,t,next)); break;
            }
            means = next;
        }
        let (boundary,reference,means) = accepted.expect("physical two-wall fixed point converges");
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let n = mesh.vertex_count(); let mut weights = vec![0.0;n]; weights[n-1] = 1.0;
        let mut remainders = Vec::new();
        for step in [1.0,0.5,0.25] {
            let approximate: Vec<_> = reference.iter().zip(mesh.positions())
                .map(|(t,x)| t-step*(1.0+0.3*x[0]+0.7*x[1])).collect();
            let other = convection([mean(&mesh,&approximate,0.0),mean(&mesh,&approximate,1.0)])
                .map(|p| [p[0],p[1]]);
            let run = |feedback: &[f64]| AmbientRadiationGradient::compare_goal_at(cx,problem,None,
                comparison_config(),&reference,&approximate,&["left","right"],&convection(means),
                &other,feedback,&patches,&weights,8*n).unwrap();
            let full = run(&[0.10,0.25,0.05,0.12]);
            let frozen = run(&[0.0;4]);
            assert!(full.linearization_remainder.abs() < 0.1*frozen.linearization_remainder.abs());
            remainders.push(full.linearization_remainder.abs());
        }
        for pair in remainders.windows(2) { assert!(pair[1] < 0.35*pair[0] && pair[1] > 0.15*pair[0]); }
        let p = convection(means);
        let other = p.map(|p| [p[0],p[1]]);
        let run = |cx: &Cx<'_>, t: &[f64], entries| AmbientRadiationGradient::compare_goal_at(
            cx,problem,None,comparison_config(),t,&reference,&["left","right"],&p,&other,
            &[0.10,0.25,0.05,0.12],&patches,&vec![0.0;n],entries);
        assert!(run(cx,&reference,8*n-1).is_err());
        assert!(run(cx,&vec![400.0;n],8*n).is_err(), "zero goal cannot admit an unconverged primal");
        assert!(run(cx,&vec![501.0;n],8*n).is_err(), "the original emissivity domain is not extended");
        let gate = CancelGate::new_clock_free(); gate.request();
        context(&gate,|cancelled| assert!(matches!(run(cancelled,&reference,8*n),
            Err(ConductionError::Cancelled { .. }))));
    });
}
