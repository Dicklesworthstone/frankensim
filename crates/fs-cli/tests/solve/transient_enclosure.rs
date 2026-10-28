//! Native capacity and total-enthalpy enclosure histories, controls and reports.
//! Equal-area tetrahedral patches are an algebraic exchange fixture, not a
//! claim that these two finite faces geometrically close a real enclosure.

use super::*;

#[path = "transient_enclosure/adjoint.rs"]
mod adjoint;
#[path = "transient_enclosure/capacity.rs"]
mod capacity;
#[path = "transient_enclosure/capacity_adjoint.rs"]
mod capacity_adjoint;
use fs_project::{
    ConductionEnclosure, ConductionTransient, EnclosureRadiatingSurface,
    EnclosureViewFactorEvidence, TransientEnthalpy, TransientEnthalpyKnot,
    TransientPowerStep, TransientRegionCapacity, TransientRegionEnthalpy,
    TransientRegionPower, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const INITIAL_K: f64 = 293.15;
const DENSITY: f64 = 6.0;
const CP: f64 = 10.0;
const EPSILON: f64 = 0.85;
const SIGMA: f64 = 5.670_374_419e-8;
const MATERIAL_SOURCE: &str = "synthetic total-enthalpy enclosure oracle; not measured properties";
const POWER_SOURCE: &str = "synthetic delivered pulse: 10 W to 0.3 s, zero to 0.7 s, 20 W to 1 s";

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap_or_else(|| panic!("missing numeric {key}: {row:?}"))
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}");
}

fn values(row: &J, key: &str) -> Vec<f64> {
    row.get(key).unwrap().as_array().unwrap().iter()
        .map(|value| value.as_f64().unwrap()).collect()
}

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap()
}

fn knots(rows: &[(f64, f64, f64)]) -> Vec<TransientEnthalpyKnot> {
    rows.iter().map(|&(enthalpy, temperature, fraction)| TransientEnthalpyKnot {
        specific_enthalpy: QtyAny::new(enthalpy, dims::SPECIFIC_ENTHALPY),
        temperature: QtyAny::new(temperature, dims::TEMPERATURE),
        liquid_mass_fraction: fraction,
    }).collect()
}

fn project(card: &str) -> ProjectSpec {
    let mut spec = conduction_fixture_project(7, &tetra_stl());
    transient_product::declaration(&mut spec, DENSITY * CP, 1.0);
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    for (name, normal) in [("wall-x", [1.0, 0.0, 0.0]), ("wall-z", [0.0, 0.0, 1.0])] {
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".into(), name: name.into(),
            display: "Manufactured exchanging patch".into(), expect_id: None,
        });
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".into(), target: name.into(), length_unit: "m".into(),
            selector: MeshSelector::HalfSpace {
                normal, offset: 0.0, side: HalfSpaceSide::AtMost, tolerance: 0.0,
            },
            allow_overlap: true,
        });
    }
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries = [("wall-x", 290.0), ("wall-z", 310.0)].into_iter()
        .map(|(target, reference)| ThermalBoundary {
            target: target.into(),
            condition: ThermalBoundaryCondition::Convection {
                coefficient: QtyAny::new(10.0, dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(reference, dims::TEMPERATURE),
            },
        }).collect();
    setup.adiabatic_remainder = true;
    setup.enclosure_radiation = Some(ConductionEnclosure {
        surfaces: ["wall-x", "wall-z"].into_iter().enumerate()
            .map(|(index, target)| EnclosureRadiatingSurface {
                name: format!("radiator-{index}"), target: target.into(),
                card: card.into(), claim: None,
                query_temperature: QtyAny::new(300.0, dims::TEMPERATURE),
            }).collect(),
        view_factors: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        evidence: EnclosureViewFactorEvidence::Analytic {
            geometry: "manufactured equal-area exchange; not finite-geometry view factors".into(),
        },
        closure_tolerance: 1e-12,
        reciprocity_relative_tolerance: 1e-12,
        max_iterations: 128,
        temperature_tolerance: QtyAny::new(1e-9, dims::TEMPERATURE),
        heat_tolerance: QtyAny::new(1e-6, dims::POWER),
        relaxation: 0.5,
    });
    let policy = time(&mut spec);
    policy.max_step = QtyAny::new(0.4, dims::TIME);
    policy.max_steps = 27;
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            phase_conductivity: None,
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0).unwrap(),
            source: MATERIAL_SOURCE.into(),
            reference_density: QtyAny::new(DENSITY, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            knots: knots(&[(-600.0, 233.15, 0.0), (1050.0, 398.15, 0.0)]),
        }],
    });
    policy.power_schedules = vec![TransientRegionPower {
        region: "air".into(), source: POWER_SOURCE.into(),
        steps: [(0.3, 10.0), (0.7, 0.0), (1.0, 20.0)].into_iter()
            .map(|(until, watts)| TransientPowerStep {
                until: QtyAny::new(until, dims::TIME),
                watts: QtyAny::new(watts, dims::POWER),
            }).collect(),
    }];
    // The explicit schedule replaces, rather than scales or adds, this static load.
    spec.power.as_mut().unwrap()[0].duty = 0.25;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    spec.envelope.as_mut().unwrap().ambient_lo.value = 300.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 300.0;
    spec
}

#[derive(Debug, PartialEq)]
struct Solved {
    run: String,
    receipts: Vec<String>,
    qoi: String,
    conduction: String,
    field: String,
    html: String,
}

fn run(spec: &ProjectSpec, cards: &CardPackSet) -> Solved {
    let text = print_sexpr(spec).unwrap();
    assert!(text.contains(":enclosure-radiation (enclosure-radiation"));
    let decoded = fs_project::parse_sexpr(&text).unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, tetra_stl());
    let (run, qoi, conduction, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    let retained = J::parse(&conduction).unwrap();
    let history = retained.get("transient").unwrap();
    let last = history.get("fine").unwrap().as_array().unwrap().last().unwrap();
    let energy = retained.get("energy").unwrap();
    // The common endpoint report and trajectory must describe the same fresh
    // physical heat law, including dynamic convection and prescribed storage.
    let net_w = n(energy, "source_w") + n(energy, "dirichlet_in_w")
        - n(energy, "neumann_out_w") - n(energy, "robin_out_w");
    close(net_w, n(last, "net_input_w"), 2e-8);
    close(n(energy, "storage_w"), n(last, "stored_energy_change_j") / n(last, "dt_s"), 2e-8);
    close(n(energy, "closure_w"), n(last, "energy_residual_j") / n(last, "dt_s"), 2e-8);
    let receipts = stage_receipt_hashes(&ledger, &run);
    let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
    let html = String::from_utf8(artifact_bytes(&ledger,
        &receipt_str_field(&report, "report_html"))).unwrap();
    let before = ledger.visible_op_ids(fs_ledger::MAIN_BRANCH, None).unwrap();
    let mut progress = Vec::new();
    let cached = run_solve(&ledger, &CancelGate::new_clock_free(), &mut benign_clock(),
        &decoded, cards, &mut progress).unwrap();
    assert!(matches!(cached.status, SolveRunStatus::Completed));
    assert_eq!(cached.run, run);
    assert_eq!(cached.prior_stages, 7);
    assert!(cached.stages.is_empty());
    assert!(progress.is_empty());
    assert_eq!(before, ledger.visible_op_ids(fs_ledger::MAIN_BRANCH, None).unwrap());
    assert_eq!(receipts, stage_receipt_hashes(&ledger, &run));
    Solved { run, receipts, qoi, conduction, field, html }
}

fn check_history(solved: &Solved, card: &str, latent: bool) {
    let receipt = J::parse(&solved.conduction).unwrap();
    let field = J::parse(&solved.field).unwrap();
    let qoi = J::parse(&solved.qoi).unwrap();
    let history = receipt.get("transient").unwrap();
    assert_eq!(history.str_field("schema"), Some("fs-cli-transient-enthalpy-v1"));
    assert_eq!(history.str_field("authority"), Some("Estimated"));
    assert_eq!(history.get("radiation").unwrap().str_field("method"),
        Some("implicit-closed-gray-diffuse-enclosure"));
    close(n(history, "coarse_steps"), 3.0, 0.0);
    close(n(history, "fine_steps"), 6.0, 0.0);
    close(n(history, "total_steps"), 9.0, 0.0);
    close(n(history, "final_time_s"), 1.0, 0.0);
    close(n(history.get("workload").unwrap(), "integrated_scheduled_input_j"), 9.0, 1e-12);
    let mut fine_stored = 0.0;
    for grid in ["coarse", "fine"] {
        let rows = history.get(grid).unwrap().as_array().unwrap();
        for switch in [0.3, 0.7, 1.0] {
            assert!(rows.iter().any(|row| n(row, "time_s") == switch),
                "{grid} must land on the authored event {switch}");
        }
        let (mut start, mut input, mut stored, mut external) = (0.0, 0.0, 0.0, 0.0);
        for row in rows {
            let end = n(row, "time_s");
            let dt = n(row, "dt_s");
            close(dt, end - start, 1e-15);
            let watts = if start < 0.3 { 10.0 } else if start < 0.7 { 0.0 } else { 20.0 };
            close(n(row, "source_w"), watts, 1e-10);
            let radiation = row.get("radiation").unwrap();
            assert_eq!(radiation.str_field("model"), Some("closed-gray-diffuse-enclosure"));
            close(n(radiation, "applied_radiation_out_w"), 0.0, 1e-6);
            close(n(radiation, "nonlinear_radiation_out_w"), 0.0, 1e-6);
            close(n(radiation, "prescribed_neumann_out_w"), 0.0, 1e-12);
            close(n(radiation, "physical_dirichlet_in_w"), 0.0, 1e-12);
            assert!(n(radiation, "solid_solves") >= 1.0);
            assert!(n(radiation, "max_heat_mismatch_w") <= 1e-6);
            assert!(n(radiation, "physical_residual_j")
                <= n(radiation, "physical_residual_threshold_j"));
            assert!(n(radiation, "physical_energy_residual_j").abs()
                <= n(history, "energy_tolerance_j"));
            let convection = n(radiation, "convective_out_w");
            close(n(row, "net_input_w"), watts - convection, 2e-6);
            close(n(row, "stored_energy_change_j"), dt * (watts - convection), 2e-6);
            close(n(row, "energy_residual_j"), n(row, "stored_energy_change_j")
                - dt * (watts - convection), 2e-6);
            if latent {
                close(n(row, "final_region_max_k"), INITIAL_K, 2e-9);
                close(convection, -68.5, 2e-8);
            }
            input += dt * watts;
            external += dt * convection;
            stored += n(row, "stored_energy_change_j");
            start = end;
        }
        close(start, 1.0, 0.0);
        close(input, 9.0, 1e-10);
        close(stored, input - external, 2e-5);
        if latent { close(stored, 77.5, 2e-6); }
        if grid == "fine" { fine_stored = stored; }
    }
    let endpoint = history.get("endpoint").unwrap();
    let masses = values(endpoint, "reference_nodal_mass_kg");
    let enthalpy = values(endpoint, "specific_enthalpy_j_kg");
    let temperature = values(&field, "temperature");
    let fraction = values(endpoint, "liquid_mass_fraction");
    assert_eq!(masses.len(), enthalpy.len());
    assert_eq!(enthalpy.len(), temperature.len());
    close(masses.iter().sum(), DENSITY / 6.0, 2e-12);
    let retained_storage: f64 = masses.iter().zip(&enthalpy).map(|(mass, h)| mass * h).sum();
    close(retained_storage, fine_stored, 2e-6);
    close(n(history.get("energy").unwrap(), "stored_change_j"), retained_storage, 2e-6);
    close(n(history.get("energy").unwrap(), "integrated_net_input_j"), retained_storage, 2e-5);
    close(n(history.get("energy").unwrap(), "window_residual_j"), 0.0, 2e-5);
    for ((&h, &t), &liquid) in enthalpy.iter().zip(&temperature).zip(&fraction) {
        if latent {
            assert!((-20.0..980.0).contains(&h), "fixture must stay on its latent plateau");
            close(t, INITIAL_K, 2e-9);
            close(liquid, (h + 20.0) / 1000.0, 2e-10);
        } else {
            close(h, CP * (t - INITIAL_K), 2e-7);
            close(liquid, 0.0, 0.0);
        }
    }
    let maximum = temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    close(n(&qoi.get("qoi").unwrap().as_array().unwrap()[0], "value"), maximum, 0.0);
    let radiation = receipt.get("radiation").unwrap();
    assert_eq!(radiation.str_field("model"), Some("closed-gray-diffuse-enclosure"));
    assert_eq!(radiation.str_field("time_integration"), Some("backward-euler-total-enthalpy"));
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 2);
    for (index, target) in ["wall-x", "wall-z"].into_iter().enumerate() {
        assert_eq!(surfaces[index].str_field("target"), Some(target));
        assert_eq!(surfaces[index].str_field("card"), Some(card));
        close(n(&surfaces[index], "area_m2"), 0.5, 1e-12);
        close(n(&surfaces[index], "emissivity"), EPSILON, 0.0);
    }
    let a = n(&surfaces[0], "mean_temperature_k");
    let b = n(&surfaces[1], "mean_temperature_k");
    let exchange = SIGMA * 0.5 * (a - b) * (a + b) * (a * a + b * b)
        / (2.0 / EPSILON - 1.0);
    if !latent { assert!(exchange.abs() > 1e-5, "fixture must transfer radiation"); }
    close(n(&surfaces[0], "nonlinear_heat_w"), exchange, 1e-6);
    close(n(&surfaces[1], "nonlinear_heat_w"), -exchange, 1e-6);
    close(n(radiation, "nonlinear_radiative_out_w"), 0.0, 1e-6);
    close(n(radiation, "convective_out_w"), 5.0 * (a - 290.0) + 5.0 * (b - 310.0), 2e-8);
    assert!(n(radiation, "nonlinear_residual_j")
        <= n(radiation, "nonlinear_residual_threshold_j"));
    assert!(n(radiation, "physical_energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
    let energy = receipt.get("energy").unwrap();
    close(n(energy, "source_w") + n(energy, "dirichlet_in_w")
        - n(energy, "neumann_out_w") - n(energy, "robin_out_w") - n(energy, "storage_w"),
        n(energy, "closure_w"), 2e-10);
    close(n(energy, "closure_w"), 0.0, 2e-5);
    assert!(solved.html.contains("closed-gray-diffuse-enclosure"));
    assert!(solved.html.contains("total-enthalpy"));
    assert!(solved.html.contains("Estimated"));
}

#[test]
fn g1_g3_g5_native_enthalpy_enclosure_preserves_exchange_workload_energy_and_cached_replay() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EPSILON);
    let spec = project(&card);
    let first = run(&spec, &cards);
    check_history(&first, &card, false);
    let second = run(&spec, &cards);
    assert_eq!(first, second, "same project, cards and mesh must replay every retained artifact");
}

#[test]
fn g1_native_enclosure_latent_storage_accumulates_once_while_temperature_is_constant() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EPSILON);
    let mut spec = project(&card);
    time(&mut spec).enthalpy.as_mut().unwrap().materials[0].knots = knots(&[
        (-500.0, 243.15, 0.0), (-20.0, INITIAL_K, 0.0),
        (980.0, INITIAL_K, 1.0), (1480.0, 343.15, 1.0),
    ]);
    check_history(&run(&spec, &cards), &card, true);
}

#[test]
fn g3_native_enthalpy_enclosure_ladder_rebinds_trace_and_shares_trajectory_allowance() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EPSILON);
    let mut spec = project(&card);
    spec.solver.as_mut().unwrap().fidelity = "ladder".into();
    time(&mut spec).max_steps = 18;
    let solved = run(&spec, &cards);
    check_history(&solved, &card, false);
    let receipt = J::parse(&solved.conduction).unwrap();
    let ladder = receipt.get("ladder").unwrap();
    let rungs = ladder.get("rungs").unwrap().as_array().unwrap();
    assert_eq!(rungs.len(), 2);
    close(n(&rungs[1], "tets"), 8.0 * n(&rungs[0], "tets"), 0.0);
    close(n(receipt.get("mesh").unwrap(), "elements"), n(&rungs[1], "tets"), 0.0);
    assert_eq!(ladder.str_field("stop"), Some("time-step-budget"));
    let work = ladder.get("space_time").unwrap();
    close(n(work, "forward_steps_completed"), 18.0, 0.0);
    close(n(work, "max_forward_steps"), 18.0, 0.0);
    assert_eq!(work.get("rungs").unwrap().as_array().unwrap().len(), 2);
}

#[test]
fn g0_native_transient_enclosure_refuses_incompatible_goals_and_adaptive_claims() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EPSILON);
    for variant in ["steady-adjoint", "adaptive", "mean-bound"] {
        let mut spec = project(&card);
        match variant {
            "adaptive" => spec.solver.as_mut().unwrap().fidelity = "adaptive".into(),
            _ => spec.outputs.as_mut().unwrap().push(OutputRequest {
                name: if variant == "steady-adjoint" {
                    "temperature-max-adjoint"
                } else {
                    "temperature-volume-mean-bound"
                }.into(),
                kind: "report".into(), region: None,
            }),
        }
        let decoded = fs_project::parse_sexpr(&print_sexpr(&spec).unwrap()).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, tetra_stl());
        let error = run_solve(&ledger, &CancelGate::new_clock_free(), &mut benign_clock(),
            &decoded, &cards, &mut Vec::new()).unwrap_err();
        assert!(error.what.contains("steady") || error.what.contains("adaptive")
            || error.what.contains("bound") || error.what.contains("trajectory"),
            "{variant}: {error:?}");
        if let Some(run) = &error.run {
            assert!(stage_receipt_hashes(&ledger, run).len() <= 4,
                "{variant}: an unsupported objective must publish no conduction field");
        }
    }
}

#[path = "transient_enclosure/vacuum.rs"]
mod vacuum;
