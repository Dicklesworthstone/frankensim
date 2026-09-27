//! G0/G3/G4: real cooling correction, publication pairing and refusal atomicity.
use super::*;
use crate::json_read::JsonValue as Json;
use fs_airflow::conjugate::AirSegment;
use fs_conduction::{ConductionMesh, ConductivityModel, InitialGuess, ScalarField,
    SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_solver::goal::GoalResidualLimits;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 7311, kernel_id: 73,
            tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx);
    });
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel,
    source: ScalarField, paths: Vec<AirPath>,
}
impl Fixture {
    fn new(multiple: bool) -> Self {
        let (complex, points) = fs_conduction::fixtures::unit_cube(1);
        let mesh = ConductionMesh::new(complex, points).unwrap();
        let boundary = if multiple {
            ThermalBoundaryBuilder::new(&mesh)
                .region("hot", |f| f.centroid[0] < 1e-9, ThermalBc::robin(2.0, 300.0).unwrap()).unwrap()
                .region("cold", |f| f.centroid[0] > 1.0-1e-9, ThermalBc::robin(3.0, 300.0).unwrap()).unwrap()
                .region("static", |f| f.centroid[1] > 1.0-1e-9, ThermalBc::robin(1.0, 310.0).unwrap()).unwrap()
                .adiabatic_remainder().finish().unwrap()
        } else {
            ThermalBoundaryBuilder::new(&mesh).region("air", |_| true,
                ThermalBc::robin(2.0, 300.0).unwrap()).unwrap().finish().unwrap()
        };
        let segment = |name: &str, h| {
            let area = mesh.boundary().iter().enumerate().filter_map(|(i, face)| {
                boundary.region_for(i).filter(|&r| boundary.region_names()[r] == name).map(|_| face.area)
            }).sum();
            AirSegment::new(name, area, h).unwrap()
        };
        let paths = if multiple { vec![
            AirPath::new(330.0, 1.0, 12.0, vec![segment("hot", 2.0)]).unwrap(),
            AirPath::new(290.0, 1.0, 8.0, vec![segment("cold", 3.0)]).unwrap(),
        ] } else { vec![AirPath::new(330.0, 1.0, 24.0, vec![segment("air", 2.0)]).unwrap()] };
        Self { mesh, boundary, paths, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(6.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            source: &self.source, element_materials: None }
    }
    fn original(&self, cx: &Cx<'_>) -> ConductionSolution {
        fs_conduction::solve(cx, self.problem(), SolveConfig {
            initial: InitialGuess::Uniform(300.0), ..SolveConfig::default()
        }).unwrap()
    }
    fn history(&self) -> String {
        let branches: Vec<_> = self.paths.iter().enumerate().map(|(i, path)| {
            let rows: Vec<_> = path.segments().iter().map(|s| format!(
                "{{\"target\":{},\"card\":\"synthetic-test-card\",\"order\":0,\"air_in_k\":-999,\"air_out_k\":-999,\"reference_k\":-999,\"solid_heat_rate_w\":-999,\"air_heat_rate_w\":-999}}",
                json_string(s.region()))).collect();
            format!("{{\"branch\":\"branch-{i}\",\"path\":\"vent:branch-{i}\",\"inlet_k\":{},\"outlet_k\":-999,\"flow_m3_s\":{{\"lo\":1e-6,\"mid\":2e-6,\"hi\":3e-6}},\"iterations\":17,\"segments\":[{}]}}",
                path.inlet_temperature_k(), rows.join(","))
        }).collect();
        if branches.len() == 1 { branches[0].clone() }
        else { format!("{{\"schema\":\"independent-branches-shared-solid-v1\",\"branch_count\":{},\"branches\":[{}],\"iterations\":17}}", branches.len(), branches.join(",")) }
    }
}
fn configs() -> (LinearConfig, LinearGoalAnalysisConfig) {
    (LinearConfig { tolerance: 1e-12, max_iterations: 500, restart: 20 },
        LinearGoalAnalysisConfig { residual_limits: GoalResidualLimits {
            max_rows: 100, max_nonzeros: 100_000 }, max_stability_iterations: 500 })
}
fn prepared(f: &Fixture, cx: &Cx<'_>, original: &ConductionSolution, history: &str) -> Publication {
    let (linear, solid) = configs();
    prepare(cx, f.problem(), None, &f.paths, linear, original, history,
        &(0..f.mesh.vertex_count()).collect::<Vec<_>>(), 64*1024*1024, solid, 1e-7).unwrap()
}

#[test]
fn one_field_drives_the_qoi_bound_physical_report_and_live_air_receipt() {
    let f = Fixture::new(false);
    with_cx(|_, cx| {
        let original = f.original(cx); let saved = original.clone(); let history = f.history();
        let publication = prepared(&f, cx, &original, &history);
        let replacement = publication.replacement.unwrap();
        let control = Json::parse(publication.evidence.control_json.as_ref().unwrap()).unwrap();
        assert_eq!(control.str_field("schema"), Some(SCHEMA));
        assert_eq!(control.get("physical_accepted"), Some(&Json::Bool(true)));
        assert_eq!(control.get("candidate_accepted"), Some(&Json::Bool(true)));
        let gain = control.f64_field("feedback_gain_infinity_upper").unwrap();
        assert!(gain > 0.0 && gain < 1.0);
        assert_eq!(control.str_field("inverse_method"), Some("state-contraction"));
        assert_eq!(control.get("schur_inverse_infinity_upper"), Some(&Json::Null));
        assert!(control.f64_field("maximum_response_residual_upper").unwrap() >= 0.0);
        assert!(publication.evidence.primal_iterations > 0);
        let receipt = Json::parse(&replacement.conjugate).unwrap();
        assert_eq!(receipt.get("initial_exchange"), Some(&Json::parse(&history).unwrap()));
        assert_eq!(receipt.path(&["flow_m3_s", "lo"]).unwrap().number_raw(), Some("1e-6"));
        assert!((receipt.f64_field("outlet_k").unwrap() - 330.25).abs() < 1e-6);
        assert!((receipt.f64_field("air_total_w").unwrap() - 6.0).abs() < 1e-6);
        let row = &receipt.get("segments").unwrap().as_array().unwrap()[0];
        let flux = &replacement.solution.report.robin_fluxes[0];
        assert_eq!(row.f64_field("solid_heat_rate_w"), Some(flux.heat_rate_w));
        assert_eq!(row.str_field("card"), Some("synthetic-test-card"));
        assert_eq!(row.f64_field("wall_temperature_k"), Some(flux.mean_wall_temperature_k));
        let fresh = fs_conduction::solve(cx, ConductionProblem { boundary: &replacement.boundary, ..f.problem() },
            SolveConfig { initial: InitialGuess::Uniform(330.0), ..SolveConfig::default() }).unwrap();
        for (a, b) in fresh.temperature.iter().zip(&replacement.solution.temperature) { assert!((a-b).abs() < 1e-6); }
        assert_eq!(original, saved);
        let replay = prepared(&f, cx, &original, &history);
        assert_eq!(replay.replacement.unwrap().conjugate, replacement.conjugate);
        assert_eq!(replay.evidence.control_json, publication.evidence.control_json);
    });
}

#[test]
fn independent_inlets_keep_metadata_and_static_robin_fluxes() {
    let f = Fixture::new(true);
    with_cx(|_, cx| {
        let original = f.original(cx);
        let replacement = prepared(&f, cx, &original, &f.history()).replacement.unwrap();
        let receipt = Json::parse(&replacement.conjugate).unwrap();
        let branches = receipt.get("branches").unwrap().as_array().unwrap();
        assert_eq!(branches[0].f64_field("inlet_k"), Some(330.0));
        assert_eq!(branches[1].f64_field("inlet_k"), Some(290.0));
        for (i, branch) in branches.iter().enumerate() {
            assert_eq!(branch.str_field("branch"), Some(format!("branch-{i}").as_str()));
            assert_eq!(branch.get("decomposition_residual_w"), Some(&Json::Null));
            assert_eq!(branch.path(&["acceleration", "method"]).unwrap().as_str(), Some("fgmres-physical-goal-polish"));
        }
        assert!(receipt.f64_field("decomposition_residual_w").unwrap().abs() < 1e-6);
        let flux = replacement.solution.report.robin_fluxes.iter().find(|row| row.region == "static").unwrap();
        assert!((flux.mean_reference_temperature_k-310.0).abs() < 1e-10);
    });
}

#[test]
fn a_rejected_changed_candidate_never_supplies_the_unchanged_fields_bound() {
    let f = Fixture::new(false);
    with_cx(|_, cx| {
        let original = f.original(cx); let (linear, solid) = configs();
        let feedback = policy(64*1024*1024, solid, linear);
        let limits = spectral_policy(64*1024*1024, feedback);
        let gates = PhysicalCoolingGates { energy_relative_tolerance: 1e-6, reference_tolerance_k: 1e-8,
            balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-7 };
        let mut result = polish_linear_maximum_with_spectral(cx, f.problem(), None, &f.paths,
            linear, &original, &(0..f.mesh.vertex_count()).collect::<Vec<_>>(), solid, feedback,
            LinearGoalSolveConfig { absolute_tolerance: 1e-7, max_primal_iterations: 500,
                check_every: 8, max_defect_corrections: 2 },
            SpectralMaximumControl { initial_shift: 0.001, limits }, gates).unwrap();
        assert!(result.accepted.is_some());
        let initial = result.correction.initial_analysis.algebraic_half_width_k();
        let original_coupled = result.correction.initial_analysis.coupled();
        let original_gain = original_coupled.gain_infinity_upper();
        let original_method = original_coupled.inverse_method().map(|method| method.tag());
        let original_schur = original_coupled.schur_inverse_infinity_upper();
        let original_response = original_coupled.response_residual_infinity_upper()
            .iter().copied().reduce(f64::max);
        assert_ne!(initial, result.correction.solution.solid.analysis.algebraic_half_width_k());
        // Exercise the projector's refusal contract independently of which gate failed.
        result.accepted = None;
        result.physical_refusal = Some(fs_airflow::conjugate::goal::maximum::physical::PhysicalCoolingRefusal::Decomposition {
            imbalance_w: 1.0, limit_w: 1e-8 });
        let output = project(cx, result, &f.history(), &f.paths, 1e-7, 500, gates, limits).unwrap();
        assert!(output.replacement.is_none());
        let control = Json::parse(output.evidence.control_json.as_ref().unwrap()).unwrap();
        assert_eq!(control.f64_field("final_bound_k"), initial);
        for key in ["feedback_gain_infinity_upper", "inverse_method",
            "schur_inverse_infinity_upper", "maximum_response_residual_upper"]
        { assert!(control.get(key).is_some(), "missing coupled witness {key}"); }
        assert_eq!(control.f64_field("feedback_gain_infinity_upper"), original_gain);
        assert_eq!(control.str_field("inverse_method"), original_method);
        assert_eq!(control.f64_field("schur_inverse_infinity_upper"), original_schur);
        assert_eq!(control.f64_field("maximum_response_residual_upper"), original_response);
        assert!(output.evidence.primal_iterations > 0);
        assert_eq!(control.get("candidate_accepted"), Some(&Json::Bool(false)));
    });
}

#[test]
fn malformed_history_or_cancellation_cannot_partially_replace_caller_state() {
    let f = Fixture::new(false);
    with_cx(|gate, cx| {
        let original = f.original(cx); let saved = original.clone(); let (linear, solid) = configs();
        let vertices: Vec<_> = (0..f.mesh.vertex_count()).collect();
        for history in ["not-json".into(), f.history().replace("\"target\":\"air\"", "\"target\":\"wrong\"")] {
            assert!(prepare(cx, f.problem(), None, &f.paths, linear, &original, &history, &vertices,
                64*1024*1024, solid, 1e-7).is_err());
            assert_eq!(original, saved);
        }
        gate.request();
        assert!(prepare(cx, f.problem(), None, &f.paths, linear, &original, &f.history(), &vertices,
            64*1024*1024, solid, 1e-7).is_err());
        assert_eq!(original, saved);
    });
}
