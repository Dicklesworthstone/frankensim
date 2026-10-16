//! The CLI adapter against the actual fan, correlation, and FEM producers.
//! The slab oracle uses closed-form resistance addition, not a second call to
//! the production conjugate driver. Its two inlets differ and one branch gives
//! heat TO the solid, so accidentally chaining the air paths changes the answer.

use super::*;
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, InitialGuess, SolveConfig, solve};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_qty::{Pressure, VolumetricFlowRate};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena, StreamKey {
            seed: 11, kernel_id: 73, tile: 0, iteration: 0,
        }, Budget::INFINITE, ExecMode::Deterministic))
    })
}

fn card() -> CorrelationId {
    CorrelationId::ALL.iter().copied()
        .find(|id| id.name() == "convection.circular-duct-hausen-developing")
        .expect("existing developing laminar card")
}

fn law(target: &str, branch: &str, order: u32, inlet: f64) -> AirflowLaw {
    AirflowLaw {
        target: target.to_string(), branch: branch.to_string(), order,
        inlet_temperature_k: inlet, hydraulic_diameter_m: 0.01,
        flow_area_m2: 0.001, channel_length_m: 0.1, correlation: card(),
    }
}

fn operating() -> OperatingPoint {
    let source = SourceProvenance::new("analytic adapter test", "adapter-test-v1");
    let curve = FanCurve::new("test", vec![
        FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(120.0)),
        FanPoint::new(VolumetricFlowRate::new(0.002), Pressure::new(0.0)),
    ], source.clone(), 0.0, ToleranceBasis::Analytic,
        VolumetricFlowRate::new(1.0e-8), (0.5, 2.0)).expect("curve");
    let fan = FanBank::new(curve, 1, FanArrangement::Series, 1.0).expect("bank");
    let loss = |name, resistance| LossElement::new(name, LossResistance::new(resistance),
        0.0, source.clone(), ToleranceBasis::Analytic).expect("loss");
    let primary = LossNetwork::parallel(vec![
        LossNetwork::Element(loss("vent:cold", 4.0e8)),
        LossNetwork::Element(loss("vent:hot", 1.0e8)),
    ]).expect("parallel branches");
    solve_operating_point(&fan, &EnclosureNetwork::new(primary,
        LeakageElement::new(loss("leakage", 1.0e10)))).expect("certified nominal operating point")
}

fn paths() -> ConjugatePath {
    derive_air_path(&[
        law("hot-face", "hot", 0, 330.0),
        law("cold-face", "cold", 0, 290.0),
    ], &operating(), 1.2, |_| Some(0.01), 1.0).expect("two derived paths")
}

fn fixed_walls(path: &ConjugatePath, refs: &BTreeMap<String, f64>) -> Vec<SolidRegionState> {
    path.segments.iter().map(|segment| SolidRegionState {
        region: segment.target.clone(), area_m2: segment.wetted_area_m2,
        mean_wall_temperature_k: 310.0,
        heat_rate_w: segment.htc_w_m2_k * segment.wetted_area_m2 * (310.0 - refs[&segment.target]),
        mean_reference_temperature_k: Some(refs[&segment.target]),
    }).collect()
}

#[test]
fn ordering_and_inlet_consistency_are_per_branch() {
    let mut laws = vec![
        law("hot-face", "hot", 0, 330.0),
        law("cold-downstream", "cold", 5, 290.0),
        law("cold-upstream", "cold", 0, 290.0),
    ];
    order_laws(&mut laws).expect("independent branch inlets and order zero are legal");
    assert_eq!(laws.iter().map(|l| l.target.as_str()).collect::<Vec<_>>(),
        vec!["cold-upstream", "cold-downstream", "hot-face"]);
    laws[1].inlet_temperature_k = 291.0;
    assert_eq!(order_laws(&mut laws).unwrap_err().code, "cli-solve-conduction-airflow-inlet");
    laws[1].inlet_temperature_k = 290.0;
    laws[1].order = 0;
    assert_eq!(order_laws(&mut laws).unwrap_err().code, "cli-solve-conduction-airflow-order");
}

#[test]
fn duplicate_solid_ownership_is_not_silently_overwritten_in_reference_map() {
    let mut laws = vec![law("shared", "cold", 0, 290.0), law("shared", "hot", 0, 330.0)];
    assert_eq!(order_laws(&mut laws).unwrap_err().code, "cli-solve-conduction-airflow-target");
}

#[test]
fn each_branch_derives_from_its_own_operating_flow_and_card_inputs() {
    let path = paths();
    assert_eq!(path.branches.len(), 2);
    let cold = &path.branches[0];
    let hot = &path.branches[1];
    assert_eq!(path.segments[0].target, "cold-face");
    assert_eq!(cold.path_name, "vent:cold");
    assert_eq!(hot.path_name, "vent:hot");
    assert!((hot.flow_mid_m3_s / cold.flow_mid_m3_s - 2.0).abs() < 1.0e-10);
    assert!((hot.mass_flow_kg_s / cold.mass_flow_kg_s - 2.0).abs() < 1.0e-10);
    assert!((hot.segments[0].reynolds / cold.segments[0].reynolds - 2.0).abs() < 1.0e-10);
    assert!(hot.segments[0].htc_w_m2_k > cold.segments[0].htc_w_m2_k);
    assert_eq!(cold.air_path.inlet_temperature_k(), 290.0);
    assert_eq!(hot.air_path.inlet_temperature_k(), 330.0);
    assert_eq!(derived_coefficients(&path).len(), 2);
}

#[test]
fn missing_second_hydraulic_branch_refuses_instead_of_reusing_the_first() {
    let error = derive_air_path(&[
        law("cold-face", "cold", 0, 290.0),
        law("other-face", "missing", 0, 330.0),
    ], &operating(), 1.2, |_| Some(0.01), 1.0).unwrap_err();
    assert_eq!(error.code, "cli-solve-conduction-airflow-handoff");
}

#[test]
fn shared_fem_slab_matches_two_independent_air_resistances_in_series_with_the_solid() {
    let path = paths();
    const LENGTH: f64 = 0.02;
    const AREA: f64 = 0.01;
    const K: f64 = 2.0;
    let (complex, positions) = box_grid([2, 1, 1], [LENGTH, 0.1, 0.1]);
    let mesh = ConductionMesh::new(complex, positions).expect("slab mesh");
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let coefficients = derived_coefficients(&path);
    with_cx(|cx| {
        let solve_once = |refs: &BTreeMap<String, f64>| {
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region("cold-face", |f| on_box_face(f.centroid[0], 0.0),
                    ThermalBc::robin(coefficients["cold-face"], refs["cold-face"]).expect("cold Robin"))
                .expect("cold region")
                .region("hot-face", |f| on_box_face(f.centroid[0], LENGTH),
                    ThermalBc::robin(coefficients["hot-face"], refs["hot-face"]).expect("hot Robin"))
                .expect("hot region")
                .adiabatic_remainder().finish().expect("partition");
            let mut config = SolveConfig::default();
            config.initial = InitialGuess::Uniform(310.0);
            config.linear.tolerance = 1.0e-13;
            solve(cx, ConductionProblem {
                mesh: &mesh, boundary: &boundary, material: &material,
                element_materials: None, source: &source,
            }, config).expect("actual shared-solid FEM solve")
        };
        let mut calls = 0;
        let outcome = run_exchange(cx, &path, false, |_, refs| {
            calls += 1;
            let field = solve_once(refs);
            Ok(path.segments.iter().map(|segment| {
                SolidRegionState::from_robin_flux(field.report.robin_fluxes.iter()
                    .find(|flux| flux.region == segment.target).expect("retained target flux"))
            }).collect())
        }).expect("shared conjugate solve");
        assert_eq!(outcome.solution.iterations, calls);
        let refs = path.segments.iter().zip(&outcome.solution.reference_temperatures_k)
            .map(|(s, &r)| (s.target.clone(), r)).collect();
        let field = solve_once(&refs);
        let outcome = cross_check_decomposition(outcome, field.report.energy.robin_out_w, 0.0)
            .expect("final field decomposition agrees");
        let effective = |branch: &ConjugateBranch| {
            let capacity = branch.air_path.capacity_rate_w_per_k();
            let conductance = branch.segments[0].htc_w_m2_k * AREA;
            capacity * (1.0 - (-conductance / capacity).exp())
        };
        let cold = effective(&path.branches[0]);
        let hot = effective(&path.branches[1]);
        let q = 40.0 / (1.0 / cold + LENGTH / (K * AREA) + 1.0 / hot);
        let left = 290.0 + q / cold;
        let right = 330.0 - q / hot;
        for (&temperature, position) in field.temperature.iter().zip(mesh.positions()) {
            let expected = left + (right - left) * position[0] / LENGTH;
            assert!((temperature - expected).abs() < 1.0e-7,
                "FEM temperature {temperature}, closed form {expected}");
        }
        assert!((outcome.solution.branches[0].march.total_heat_rate_w - q).abs() < 1.0e-7);
        assert!((outcome.solution.branches[1].march.total_heat_rate_w + q).abs() < 1.0e-7);
        let receipt = receipt_fragment(&path, &outcome).expect("receipt");
        assert!(receipt.starts_with("{\"schema\":\"independent-branches-shared-solid-v1\""));
        assert!(receipt.contains("\"branch_count\":2"));
        assert!(receipt.contains("\"branch\":\"cold\""));
        assert!(receipt.contains("\"branch\":\"hot\""));
        assert_eq!(receipt.matches("\"decomposition_residual_w\":null").count(), 2);
        assert_eq!(receipt.matches("\"method\":\"iqn-ils\"").count(), 3);
        assert!(receipt.contains("aggregate only"));
    });
}

#[test]
fn aggregate_publication_gate_detects_dropped_boundary_heat() {
    let path = paths();
    let outcome = with_cx(|cx| run_exchange(cx, &path, false,
        |_, refs| Ok(fixed_walls(&path, refs)))).expect("exchange");
    assert!(receipt_fragment(&path, &outcome).is_err(), "unverified final solve must not publish");
    let total: f64 = outcome.solution.branches.iter().map(|b| b.balance.solid_total_w).sum();
    let error = cross_check_decomposition(outcome.clone(), total + 1.0, 0.0).unwrap_err();
    assert_eq!(error.code, "cli-solve-conduction-airflow-decomposition");
    let checked = cross_check_decomposition(outcome, total + 1.0, 1.0).expect("known off-path heat excluded");
    assert!(receipt_fragment(&path, &checked).is_ok());
}

#[test]
fn single_branch_keeps_its_branch_fields_and_records_the_acceleration() {
    let path = derive_air_path(&[law("cold-face", "cold", 0, 290.0)],
        &operating(), 1.2, |_| Some(0.01), 1.0).expect("single path");
    let outcome = with_cx(|cx| run_exchange(cx, &path, false,
        |_, refs| Ok(fixed_walls(&path, refs)))).expect("exchange");
    let total = outcome.solution.branches[0].balance.solid_total_w;
    let checked = cross_check_decomposition(outcome, total, 0.0).expect("checked");
    let receipt = receipt_fragment(&path, &checked).expect("single receipt");
    assert!(receipt.starts_with("{\"branch\":\"cold\",\"path\":\"vent:cold\""));
    assert!(!receipt.contains("\"branches\":"));
    assert!(!receipt.contains("\"schema\":"));
    assert!(!receipt.contains("\"decomposition_residual_w\":null"));
    assert_eq!(receipt.matches("\"method\":\"iqn-ils\"").count(), 1);
    assert!(receipt.contains(CONJUGATE_NO_CLAIM));
}

#[test]
fn shared_solid_refusals_keep_the_original_code_and_message() {
    let path = paths();
    let error = with_cx(|cx| run_exchange(cx, &path, false, |_, _| {
        Err(conduction_error("test-solid-refusal", "sentinel solid diagnosis", "sentinel fix"))
    })).unwrap_err();
    assert_eq!(error.code, "test-solid-refusal");
    assert_eq!(error.what, "sentinel solid diagnosis");
}

#[test]
fn production_exchange_resolves_stiff_card_derived_paths_without_increasing_the_budget() {
    // Real fan and correlation producers, with an analytic lumped solid.
    // The deliberately large wetted area makes plain staggering nearly
    // stationary. This is a numerical coupling fixture, not a geometry claim.
    let path = derive_air_path(&[
        law("cold-face", "cold", 0, 290.0),
        law("hot-face", "hot", 0, 330.0),
    ], &operating(), 1.2, |_| Some(10.0), 1.0).expect("stiff derived paths");
    let solid = |refs: &BTreeMap<String, f64>| -> Vec<SolidRegionState> {
        path.segments.iter().map(|segment| {
            let power = if segment.target == "cold-face" { 5.0 } else { 3.0 };
            let conductance = segment.htc_w_m2_k * segment.wetted_area_m2;
            let reference = refs[&segment.target];
            SolidRegionState {
                region: segment.target.clone(),
                area_m2: segment.wetted_area_m2,
                mean_wall_temperature_k: reference + power / conductance,
                heat_rate_w: power,
                mean_reference_temperature_k: Some(reference),
            }
        }).collect()
    };
    with_cx(|cx| {
        let air_paths: Vec<_> = path.branches.iter().map(|b| b.air_path.clone()).collect();
        let plain = fs_airflow::graph::thermal::solve_conjugate_branches(
            cx, &air_paths, &ConjugateConfig::default(), |_, references| {
                let refs = path.segments.iter().zip(references)
                    .map(|(segment, &value)| (segment.target.clone(), value)).collect();
                Ok(solid(&refs))
            },
        );
        assert!(matches!(plain, Err(AirflowError::ConjugateNotConverged { .. })));
        let mut calls = 0;
        let outcome = run_exchange(cx, &path, false, |_, refs| {
            calls += 1;
            Ok(solid(refs))
        }).expect("production IQN exchange");
        assert_eq!(outcome.solution.iterations, calls);
        assert!(calls < 20, "bounded vector acceleration took {calls} solid solves");
        let heat: f64 = outcome.solution.branches.iter().map(|b| b.balance.air_total_w).sum();
        assert!((heat - 8.0).abs() < 1.0e-5);
        let outcome = cross_check_decomposition(outcome, 8.0, 0.0).expect("lumped solid balance");
        let receipt = receipt_fragment(&path, &outcome).expect("accelerated receipt");
        assert!(receipt.contains("\"scope\":\"all-branch-interfaces\""));
        assert!(receipt.contains(&format!("\"max_history\":{}", CONJUGATE_IQN_CONFIG.max_history)));
        assert!(receipt.contains("\"fallback\":\"fixed\""));
    });
}

#[test]
fn model_form_scale_reaches_the_air_path_the_exchange_uses() {
    // The scale must act before the air path is built: scaling only the
    // reported derivation would leave the fixed point on the nominal h.
    let nominal = derive_air_path(&[law("cold-face", "cold", 0, 290.0)],
        &operating(), 1.2, |_| Some(0.01), 1.0).expect("nominal path");
    let weaker = derive_air_path(&[law("cold-face", "cold", 0, 290.0)],
        &operating(), 1.2, |_| Some(0.01), 0.85).expect("scaled path");
    let h = nominal.segments[0].htc_w_m2_k;
    assert_eq!(weaker.segments[0].htc_w_m2_k.to_bits(), (h * 0.85).to_bits());
    assert_eq!(derived_coefficients(&weaker)["cold-face"].to_bits(), (h * 0.85).to_bits());
    let exchange = |path: &ConjugatePath| with_cx(|cx| run_exchange(cx, path, false,
        |_, refs| Ok(fixed_walls(path, refs)))).expect("exchange");
    let (a, b) = (exchange(&nominal), exchange(&weaker));
    assert_ne!(
        a.solution.reference_temperatures_k[0].to_bits(),
        b.solution.reference_temperatures_k[0].to_bits(),
        "a weaker coefficient must change the converged air reference"
    );
}
