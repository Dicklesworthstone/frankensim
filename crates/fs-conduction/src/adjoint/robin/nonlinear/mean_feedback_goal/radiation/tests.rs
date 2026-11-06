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
                // Keep the physical residual at BOTH states and freeze only
                // the radiation feedback in the reference tangent. Freezing
                // both laws instead compares two self-consistent models whose
                // estimates agree to first order, masking a missing derivative.
                let combined = |temperature: &[f64]| {
                    let wall = mean(&mesh, temperature, 1.0);
                    let hr = patches[0].secant_coefficient_w_m2_k(wall).unwrap();
                    let h = 5.0 + hr;
                    [h, (5.0 / h) * 300.0 + (hr / h) * reservoir]
                };
                let [h, r] = combined(reference);
                let other = combined(&approximate);
                let frozen = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                    comparison_config(), reference, &approximate, &["right"],
                    &[[h, r, 0.0, 0.0]], &[other], &weights, 2*n).unwrap();
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
            // Moving local reference slopes out of the matrix diagonal must
            // preserve the same derivative, including its radiative weight.
            let mut local = convection(means);
            local[0][3] = 0.10;
            local[1][3] = 0.12;
            let stale = coupled_boundary(&mesh, [[2.0,270.0,0.0,0.0], [3.0,310.0,0.0,0.0]]);
            let rebound = AmbientRadiationGradient::compare_goal_at(cx,
                ConductionProblem { boundary: &stale, ..problem }, None, comparison_config(),
                &reference, &approximate, &["left","right"], &local, &other,
                &[0.0,0.25,0.05,0.0], &patches, &weights, 8*n).unwrap();
            assert!((rebound.signed_residual_change-full.signed_residual_change).abs() < 1e-10,
                "the supplied physical convection rows and each reference derivative must enter once");
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

#[test]
fn radiation_only_goal_matches_scalar_heat_balance_and_preserves_source_support() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let boundary_at = |hot| ThermalBoundaryBuilder::new(&mesh)
            .region("left", |f| on_box_face(f.centroid[0],0.0), ThermalBc::dirichlet(hot).unwrap()).unwrap()
            .region("right", |f| on_box_face(f.centroid[0],1.0), ThermalBc::neumann(0.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let boundary = boundary_at(340.0);
        let original = boundary.clone();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let patches = [patch(300.0)];
        let n = mesh.vertex_count();
        let linear = config().linear;
        let mut weights = vec![0.0;n]; weights[n-1] = 1.0;
        // Unit-length, unit-area conduction to a radiation-only planar wall:
        // 10 (T_wall - 340) + epsilon sigma (T_wall^4 - 300^4) = 0.
        // Its exact linear P1 field and scalar tangent form an independent
        // heat-balance oracle, without another radiation or thermal solve.
        let emission = 0.8 * crate::STEFAN_BOLTZMANN_W_M2_K4;
        let residual = |wall: f64| 10.0*(wall-340.0)
            + emission*(wall.powi(4)-300.0_f64.powi(4));
        let (mut lower, mut upper) = (300.0_f64,340.0_f64);
        for _ in 0..80 {
            let wall = 0.5*(lower+upper);
            if residual(wall) > 0.0 { upper = wall; } else { lower = wall; }
        }
        let wall = 0.5*(lower+upper);
        let field = |wall| mesh.positions().iter()
            .map(|x| 340.0+x[0]*(wall-340.0)).collect::<Vec<_>>();
        let reference = field(wall);
        let compare = |approximate: &[f64], point: [f64;4], other: [f64;2], feedback: f64| {
            AmbientRadiationGradient::compare_goal_at(cx, problem, None, linear,
                &reference, approximate, &["right"], &[point], &[other], &[feedback],
                &patches, &weights, 4*n)
        };
        let mut previous = None;
        for step in [1.0,0.5,0.25] {
            let approximate = field(wall-step);
            // The placeholder temperatures deliberately differ from the bath.
            // Zero convection means they cannot affect either physical state.
            let full = compare(&approximate, [0.0,280.0,0.0,0.0], [0.0,370.0], 0.0).unwrap();
            let expected = (residual(wall)-residual(wall-step))
                /(10.0+4.0*emission*wall.powi(3));
            assert!((full.signed_residual_change-expected).abs() < 2e-10);
            assert!((full.signed_goal_change-step).abs() < 1e-12);
            assert!((full.nodal_contributions.iter().sum::<f64>()-expected).abs() < 2e-10);
            assert!(full.primal_relative_residual < linear.tolerance);
            assert!(full.dual_relative_residual < linear.tolerance);
            let remainder = full.linearization_remainder.abs();
            assert!(remainder > 1e-6);
            if let Some(old) = previous { assert!((3.8..4.2).contains(&(old/remainder))); }
            previous = Some(remainder);
        }
        for endpoint in [250.0,500.0] {
            let at_endpoint = compare(&field(endpoint), [0.0,300.0,0.0,0.0], [0.0,300.0], 0.0);
            assert!(at_endpoint.is_ok(), "an approximate source endpoint needs only a physical value: {at_endpoint:?}");
        }
        for outside in [249.0,501.0] {
            assert!(compare(&field(outside), [0.0,300.0,0.0,0.0], [0.0,300.0], 0.0).is_err());
        }
        for point in [[1.0,300.0,0.0,0.0], [0.0,300.0,1.0,0.0], [0.0,300.0,0.0,1.0]] {
            assert!(compare(&reference, point, [0.0,300.0], 0.0).is_err(),
                "an original radiation-only trace cannot acquire convection or its slopes");
        }
        assert!(compare(&reference, [0.0,300.0,0.0,0.0], [1.0,300.0], 0.0).is_err());
        assert!(compare(&reference, [0.0,300.0,0.0,0.0], [0.0,300.0], 0.1).is_err());
        let index = boundary.region_names().iter().position(|name| name == "right").unwrap();
        let prescribed_flux = boundary.with_uniform_neumann_additions(&[(index,1.0)]).unwrap();
        assert!(AmbientRadiationGradient::compare_goal_at(cx,
            ConductionProblem { boundary: &prescribed_flux, ..problem }, None, linear,
            &reference, &reference, &["right"], &[[0.0,300.0,0.0,0.0]], &[[0.0,300.0]],
            &[0.0], &patches, &weights, 4*n).is_err(), "nonzero authored heat must not be erased");
        // This is an exact isothermal equilibrium and a zero goal. Its
        // reference lies at the emissivity support edge, so no two-sided
        // tangent is supplied even though the physical secant is defined.
        let edge_boundary = boundary_at(250.0);
        assert!(AmbientRadiationGradient::compare_goal_at(cx,
            ConductionProblem { boundary: &edge_boundary, ..problem }, None, linear,
            &vec![250.0;n], &vec![250.0;n], &["right"], &[[0.0,300.0,0.0,0.0]],
            &[[0.0,300.0]], &[0.0], &[patch(250.0)], &vec![0.0;n], 4*n).is_err());
        assert_eq!(boundary, original);
    });
}
