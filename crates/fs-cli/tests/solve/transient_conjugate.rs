//! Native imported-mesh fan cooling with independent transient balances.
//! The regular tetra has one uniform thermal mode; the driven unit tetra has
//! one free temperature. Both oracles use the exponential air law directly.

use super::*;
use fs_project::spec::dims;

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

#[path = "transient_conjugate/material_controls.rs"]
mod material_controls;
#[path = "transient_conjugate/prescribed_controls.rs"]
mod prescribed_controls;

const INLET_K: f64 = 293.15;
const INITIAL_K: f64 = 330.0;
const DRIVEN_INLET_K: f64 = 315.0;
const CAPACITY_J_K: f64 = 100.0 * 8.0 / 3.0;
const AIR_CP: f64 = 1007.0;
const EMISSIVITY: f64 = 0.85;
const SIGMA: f64 = 5.670_374_419e-8;

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn bisect(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) <= 0.0 && residual(high) >= 0.0);
    for _ in 0..100 {
        let value = 0.5 * (low + high);
        if residual(value) > 0.0 {
            high = value;
        } else {
            low = value;
        }
    }
    0.5 * (low + high)
}

fn area_m2() -> f64 {
    8.0 * 3.0_f64.sqrt()
}

fn regular_project() -> (ProjectSpec, Vec<u8>) {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    transient_product::declaration(&mut spec, 100.0, 1.0);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient.as_mut().unwrap().initial_temperature =
        fs_project::TransientInitialTemperature::Uniform(QtyAny::new(INITIAL_K, dims::TEMPERATURE));
    setup.regions[0].seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    (spec, bytes)
}

struct ThermalRun {
    qoi: J,
    receipt: J,
    temperature: Vec<f64>,
}

impl ThermalRun {
    fn transient(&self) -> &J {
        self.receipt.get("transient").unwrap()
    }

    fn rows(&self, grid: &str) -> &[J] {
        self.transient().get(grid).unwrap().as_array().unwrap()
    }
}

fn run(spec: &ProjectSpec, cards: &CardPackSet, bytes: &[u8]) -> ThermalRun {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, bytes.to_vec());
    let (_, qoi, receipt, field) =
        run_conjugate_to_completion(&ledger, &decode(spec), cards);
    let field = J::parse(&field).unwrap();
    ThermalRun {
        qoi: J::parse(&qoi).unwrap(),
        receipt: J::parse(&receipt).unwrap(),
        temperature: field
            .get("temperature")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_f64().unwrap())
            .collect(),
    }
}

fn checked_coupling(row: &J) -> &J {
    let coupled = row.get("conjugate").expect("retain each endpoint's air exchange");
    let exchange = coupled.get("exchange").unwrap();
    let iterations = n(exchange, "iterations");
    assert!(iterations > 1.0 && iterations <= 100.0);
    assert!(n(coupled, "solid_solves") >= iterations);
    close(n(row, "krylov_iterations"), n(coupled, "krylov_iterations"), 0.0);
    assert!(n(coupled, "coupled_energy_residual_j").abs() <= 1e-6);
    let physical_target = n(coupled, "physical_residual_tolerance_j");
    let physical_residual = n(coupled, "physical_residual_norm_j");
    assert!(physical_target.is_finite() && physical_target > 0.0);
    assert!(physical_residual.is_finite() && physical_residual >= 0.0
        && physical_residual <= physical_target);
    assert!(n(coupled, "physical_energy_residual_j").abs() <= 1e-6);
    // The actual predictor and accepted field are audited directly;
    // intervening reference changes use the analytic load bound.
    close(n(coupled, "endpoint_evaluations"), 2.0, 0.0);
    let tolerance = n(exchange, "balance_tolerance_w");
    assert!(n(exchange, "max_region_imbalance_w").abs() <= tolerance);
    assert!(n(exchange, "decomposition_residual_w").abs() <= tolerance);
    close(n(coupled, "air_heat_gain_w"), n(exchange, "air_total_w"), 1e-12);
    close(n(coupled, "off_path_convective_out_w"), 0.0, 1e-12);
    let applied = coupled.get("applied_references_k").unwrap();
    for segment in exchange.get("segments").unwrap().as_array().unwrap() {
        close(
            n(applied, segment.str_field("target").unwrap()),
            n(segment, "reference_k"),
            n(coupled, "reference_tolerance_k") + 2e-12,
        );
    }
    coupled
}

/// Re-derive the path conductance from its fan/vent handoff and card evidence.
/// No production air-march or conduction solver supplies the oracle.
fn path_conductance(exchange: &J, area: f64, count: usize, inlet_k: f64) -> (f64, f64) {
    assert_eq!(exchange.str_field("branch"), Some("air"));
    assert_eq!(exchange.str_field("path"), Some("vent:air"));
    close(n(exchange, "inlet_k"), inlet_k, 0.0);
    close(
        n(exchange.get("air_properties").unwrap(), "specific_heat_j_kg_k"),
        AIR_CP,
        0.0,
    );
    close(
        n(exchange, "mass_flow_kg_s"),
        n(exchange, "air_density_kg_m3")
            * n(exchange.get("flow_m3_s").unwrap(), "mid"),
        1e-14,
    );
    let capacity_rate = n(exchange, "mass_flow_kg_s") * AIR_CP;
    assert!(capacity_rate > 0.0);
    let segments = exchange.get("segments").unwrap().as_array().unwrap();
    assert_eq!(segments.len(), count);
    let h = n(&segments[0], "htc_w_m2_k");
    assert!(h > 0.0);
    for (index, segment) in segments.iter().enumerate() {
        assert_eq!(segment.str_field("card"), Some("convection.gnielinski"));
        assert!(matches!(segment.get("in_domain"), Some(J::Bool(true))));
        assert!(n(segment, "reynolds") > 3000.0);
        close(n(segment, "order"), index as f64, 0.0);
        close(n(segment, "wetted_area_m2"), area / count as f64, 2e-12);
        close(n(segment, "htc_w_m2_k"), h, 1e-12);
        close(h, n(segment, "nusselt") * 0.0263 / 0.02, 1e-12);
        let ntu = h * area / count as f64 / capacity_rate;
        close(n(segment, "ntu"), ntu, 2e-12);
        close(n(segment, "effectiveness"), -(-ntu).exp_m1(), 2e-12);
    }
    (h, capacity_rate * -(-h * area / capacity_rate).exp_m1())
}

fn check_air_heat(exchange: &J, wall: f64, conductance: f64, inlet_k: f64) {
    let capacity_rate = n(exchange, "mass_flow_kg_s") * AIR_CP;
    let air_w = conductance * (wall - inlet_k);
    close(n(exchange, "air_total_w"), air_w, 2e-6);
    close(n(exchange, "outlet_k"), inlet_k + air_w / capacity_rate, 2e-8);
    let mut inlet = inlet_k;
    for segment in exchange.get("segments").unwrap().as_array().unwrap() {
        close(n(segment, "air_in_k"), inlet, 2e-8);
        let outlet = inlet + (wall - inlet) * -(-n(segment, "ntu")).exp_m1();
        close(n(segment, "air_out_k"), outlet, 2e-8);
        close(n(segment, "air_heat_rate_w"), capacity_rate * (outlet - inlet), 2e-6);
        inlet = outlet;
    }
    close(inlet, n(exchange, "outlet_k"), 2e-8);
}

fn check_final_uniform(solved: &ThermalRun) {
    assert_eq!(solved.temperature.len(), 4, "oracle requires one imported tetra");
    let last = solved.rows("fine").last().unwrap();
    let final_k = n(last, "final_region_max_k");
    for &temperature in &solved.temperature {
        close(temperature, final_k, 2e-8);
    }
    let qoi = &solved.qoi.get("qoi").unwrap().as_array().unwrap()[0];
    close(n(qoi, "value"), final_k, 0.0);
    let scope = solved.qoi.get("thermal_time_scope").unwrap();
    assert_eq!(scope.str_field("kind"), Some("final"));
    close(n(scope, "time_s"), 1.0, 0.0);
    close(n(solved.transient(), "total_steps"), 6.0, 0.0);
    close(
        n(solved.transient().get("energy").unwrap(), "stored_change_j"),
        CAPACITY_J_K * (final_k - INITIAL_K),
        2e-6,
    );
    close(
        n(solved.receipt.get("energy").unwrap(), "storage_w"),
        n(last, "stored_energy_change_j") / n(last, "dt_s"),
        2e-9,
    );
    close(
        n(solved.receipt.get("conjugate").unwrap(), "air_total_w"),
        n(last.get("conjugate").unwrap(), "air_heat_gain_w"),
        1e-12,
    );
}

#[test]
fn g1_native_transient_airflow_matches_every_endpoint_with_immutable_history() {
    let (spec, bytes) = regular_project();
    let solved = run(&spec, &fixture_cards(), &bytes);
    for grid in ["coarse", "fine"] {
        let mut old_expected = INITIAL_K;
        let mut old_actual = INITIAL_K;
        let mut previous_time = 0.0;
        for row in solved.rows(grid) {
            let coupled = checked_coupling(row);
            let exchange = coupled.get("exchange").unwrap();
            let (_, conductance) = path_conductance(exchange, area_m2(), 1, INLET_K);
            let dt = n(row, "dt_s");
            close(n(row, "time_s") - previous_time, dt, 1e-15);
            // C(T-old)/dt = P - m*cp*(1-exp(-hA/(m*cp)))*(T-Tin).
            let expected = (CAPACITY_J_K * old_expected / dt
                + 5.0 + conductance * INLET_K)
                / (CAPACITY_J_K / dt + conductance);
            let actual = n(row, "final_region_max_k");
            close(actual, expected, 2e-8);
            check_air_heat(exchange, actual, conductance, INLET_K);
            let stored = CAPACITY_J_K * (actual - old_actual);
            close(n(row, "stored_energy_change_j"), stored, 2e-6);
            close(n(row, "source_w"), 5.0, 1e-12);
            close(n(coupled, "physical_dirichlet_in_w"), 0.0, 1e-12);
            close(n(coupled, "radiative_out_w"), 0.0, 1e-12);
            close(n(coupled, "radiation_trials"), 0.0, 0.0);
            close(n(coupled, "nonlinear_iterations"), 0.0, 0.0);
            close(n(coupled, "nonlinear_backtracks"), 0.0, 0.0);
            let residual = stored - dt * (5.0 - conductance * (actual - INLET_K));
            close(residual, 0.0, 2e-6);
            close(n(coupled, "coupled_energy_residual_j"), residual, 2e-6);
            old_expected = expected;
            old_actual = actual;
            previous_time = n(row, "time_s");
        }
        close(previous_time, 1.0, 0.0);
    }
    check_final_uniform(&solved);
    assert_eq!(
        solved.transient().get("temporal_error").unwrap().str_field("status"),
        Some("estimated-order-one")
    );
}

#[test]
fn g1_native_transient_airflow_keeps_radiation_out_of_air_enthalpy() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    // The warmer reservoir sends heat into the solid while the air removes it.
    for reservoir in [270.0, 350.0] {
        let (mut spec, bytes) = regular_project();
        radiation_product::declare(&mut spec, &card, "air", reservoir);
        let solved = run(&spec, &cards, &bytes);
        for grid in ["coarse", "fine"] {
            let mut old_expected = INITIAL_K;
            let mut old_actual = INITIAL_K;
            for row in solved.rows(grid) {
                let coupled = checked_coupling(row);
                let exchange = coupled.get("exchange").unwrap();
                let (_, conductance) = path_conductance(exchange, area_m2(), 1, INLET_K);
                let dt = n(row, "dt_s");
                let radiative = |t: f64| {
                    EMISSIVITY * SIGMA * area_m2() * (t.powi(4) - reservoir.powi(4))
                };
                let expected = bisect(200.0, 450.0, |t| {
                    CAPACITY_J_K * (t - old_expected) / dt
                        + conductance * (t - INLET_K) + radiative(t) - 5.0
                });
                let actual = n(row, "final_region_max_k");
                close(actual, expected, 3e-8);
                check_air_heat(exchange, actual, conductance, INLET_K);
                let radiation_w = radiative(actual);
                assert!(radiation_w.abs() > 1.0);
                assert_eq!(radiation_w.is_sign_positive(), reservoir < actual);
                close(n(coupled, "radiative_out_w"), radiation_w, 2e-6);
                let radiation = row.get("radiation").unwrap();
                close(n(radiation, "nonlinear_radiative_out_w"), radiation_w, 2e-6);
                close(
                    n(radiation, "convective_out_w"),
                    n(coupled, "air_heat_gain_w"),
                    2e-6,
                );
                assert!(n(coupled, "radiation_trials") > n(exchange, "iterations"));
                assert!(n(coupled, "radiation_trials") >= n(radiation, "iterations"));
                assert!(n(coupled, "solid_solves") >= n(coupled, "radiation_trials"));
                let stored = CAPACITY_J_K * (actual - old_actual);
                close(n(row, "stored_energy_change_j"), stored, 2e-6);
                let residual = stored
                    - dt * (5.0 - conductance * (actual - INLET_K) - radiation_w);
                close(residual, 0.0, 2e-6);
                close(n(coupled, "coupled_energy_residual_j"), residual, 2e-6);
                close(n(coupled, "physical_dirichlet_in_w"), 0.0, 1e-12);
                old_expected = expected;
                old_actual = actual;
            }
        }
        check_final_uniform(&solved);
    }
}

fn driven_air_project(bytes: &[u8]) -> ProjectSpec {
    let mut spec = transient_product::driven_face_project(bytes);
    let template = conjugate_fixture_project(7, bytes, "convection.gnielinski")
        .cooling.unwrap().conduction.unwrap().boundaries[0].condition.clone();
    for axis in 0..3 {
        let target = format!("air-face-{axis}");
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".to_string(),
            name: target.clone(),
            display: "Air-cooled coordinate face".to_string(),
            expect_id: None,
        });
        let mut normal = [0.0; 3];
        normal[axis] = 1.0;
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".to_string(),
            target: target.clone(),
            length_unit: "m".to_string(),
            selector: MeshSelector::HalfSpace {
                normal,
                offset: 0.0,
                side: HalfSpaceSide::AtMost,
                tolerance: 0.0,
            },
            allow_overlap: true,
        });
        let mut condition = template.clone();
        if let ThermalBoundaryCondition::AirflowConvection {
            order, inlet_temperature, ..
        } = &mut condition
        {
            *order = axis as u32;
            inlet_temperature.value = DRIVEN_INLET_K;
        }
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .boundaries.push(ThermalBoundary { target, condition });
    }
    // The 315 K inlet keeps even the first unaccelerated trial inside the
    // conductivity card: its Robin free-row reservoir is 2*315-330=300 K.
    spec.envelope.as_mut().unwrap().ambient_lo.value = DRIVEN_INLET_K;
    spec.envelope.as_mut().unwrap().ambient_hi.value = DRIVEN_INLET_K;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    spec
}

#[test]
fn g1_native_transient_airflow_and_nonlinear_conductivity_share_one_endpoint() {
    let mut free_temperatures = Vec::new();
    for kind in ["nonlinear", "constant-scalar"] {
        let bytes = tetra_stl();
        let mut spec = driven_air_project(&bytes);
        let cards = transient_product::transient_conductivity_cards(kind);
        rebind_to(&mut spec, &cards);
        let solved = run(&spec, &cards, &bytes);
        let mut field = solved.temperature.clone();
        field.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(field.len(), 4);
        for &fixed in &field[1..] {
            close(fixed, 330.0, 0.0);
        }
        for grid in ["coarse", "fine"] {
            let mut old_expected = 300.0;
            let mut old_actual = 300.0;
            for (index, row) in solved.rows(grid).iter().enumerate() {
                let coupled = checked_coupling(row);
                let exchange = coupled.get("exchange").unwrap();
                let (h, conductance) = path_conductance(exchange, 1.5, 3, DRIVEN_INLET_K);
                let dt = n(row, "dt_s");
                // The free mass row is 15/24 J/K. The P1 Robin row on the
                // three coordinate faces is Qair/3 + h*(T-330)/12; using
                // the whole patch heat as the free-node load is incorrect.
                let expected = bisect(280.0, 330.0, |t| {
                    let k = if kind == "nonlinear" {
                        0.02 * ((t + 3.0 * 330.0) / 4.0) - 4.0
                    } else {
                        2.0
                    };
                    let wall = (t + 2.0 * 330.0) / 3.0;
                    0.625 * (t - old_expected) / dt
                        + (0.5 * k + h / 12.0) * (t - 330.0)
                        + conductance * (wall - DRIVEN_INLET_K) / 3.0
                });
                // The initial clamp moves three prescribed nodes from 300
                // to 330 exactly once on each grid, not once per air trial.
                let clamp_change = if index == 0 { 90.0 } else { 0.0 };
                let stored = n(row, "stored_energy_change_j");
                let actual = old_actual + stored / 0.625 - clamp_change;
                close(actual, expected, 3e-7);
                let wall = (actual + 2.0 * 330.0) / 3.0;
                check_air_heat(exchange, wall, conductance, DRIVEN_INLET_K);
                close(n(row, "source_w"), 0.0, 0.0);
                close(
                    n(coupled, "physical_dirichlet_in_w"),
                    stored / dt + conductance * (wall - DRIVEN_INLET_K),
                    1e-5,
                );
                if kind == "nonlinear" {
                    assert!(n(coupled, "nonlinear_iterations") > 0.0);
                    let nonlinear = row.get("nonlinear").unwrap();
                    assert!(n(nonlinear, "residual_j") <= n(nonlinear, "threshold_j"));
                } else {
                    close(n(coupled, "nonlinear_iterations"), 0.0, 0.0);
                    assert!(matches!(row.get("nonlinear"), Some(J::Null)));
                }
                old_expected = expected;
                old_actual = actual;
            }
            if grid == "fine" {
                close(field[0], old_expected, 3e-7);
            }
        }
        free_temperatures.push(field[0]);
    }
    assert!(
        (free_temperatures[0] - free_temperatures[1]).abs() > 0.1,
        "freezing k(T) changes the nonuniform fan-cooled transient"
    );
}

#[test]
fn g0_native_transient_airflow_refusals_preserve_the_stage_prefix_without_an_endpoint() {
    for failure in ["radiation-work", "card-domain"] {
        let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
        let (mut spec, bytes) = regular_project();
        let expected_code = if failure == "radiation-work" {
            radiation_product::declare(&mut spec, &card, "air", INLET_K);
            spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .radiation.as_mut().unwrap().max_iterations = 1;
            "cli-solve-conduction-transient"
        } else {
            if let ThermalBoundaryCondition::AirflowConvection { channel_length, .. } =
                &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .boundaries[0].condition
            {
                channel_length.value = 0.1; // L/Dh=5, outside Gnielinski's floor.
            }
            "cli-solve-conduction-airflow-card-domain"
        };
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let error = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decode(&spec),
            &cards,
            &mut Vec::new(),
        ).unwrap_err();
        assert_eq!(error.code, expected_code, "{failure}: {error:?}");
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4,
            "{failure} must retain only the completed input/geometry/material/flow stages"
        );
        assert_eq!(
            artifacts_of_kind(&ledger, "solve-conduction-solution"),
            0,
            "{failure} must not seal a partial thermal field"
        );
    }
}


fn request_air_adjoint(spec: &mut ProjectSpec, name: &str) {
    spec.outputs.as_mut().unwrap().push(OutputRequest {
        name: name.into(), kind: "report".into(), region: None,
    });
}

fn final_k(solved: &ThermalRun) -> f64 {
    solved.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn air_control<'a>(adjoint: &'a J, list: &str, key: &str, name: &str) -> &'a J {
    adjoint.get("convection").unwrap().get(list).unwrap().as_array().unwrap()
        .iter().find(|row| row.str_field(key) == Some(name)).unwrap()
}

fn set_air_inlet(spec: &mut ProjectSpec, value: f64) {
    for row in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        if let ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } =
            &mut row.condition
        {
            inlet_temperature.value = value;
        }
    }
}

fn set_air_speed(spec: &mut ProjectSpec, speed: f64) {
    let banks = &mut spec.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks;
    assert_eq!(banks.len(), 1);
    assert_eq!(banks[0].bank_id, "fixture-bank");
    banks[0].speed_ratio = speed;
}

fn air_gradients(adjoint: &J, radiative: bool) -> [f64; 7] {
    assert!(adjoint.get("robin").unwrap().as_array().unwrap().is_empty(),
        "derived air references and coefficients are not prescribed Robin controls");
    let power = &adjoint.get("regional_power").unwrap().as_array().unwrap()[0];
    let inlet = air_control(adjoint, "air_inlets", "branch", "air");
    let fan = air_control(adjoint, "fan_banks", "bank", "fixture-bank");
    assert_eq!(fan.str_field("status"), Some("available"));
    let mut values = [
        n(adjoint, "dtemperature_duniform_initial_k_per_k"),
        n(power, "dtemperature_duniform_added_watts_k_per_w"),
        n(adjoint, "dtemperature_dcapacity_multiplier_k"),
        n(inlet, "dtemperature_dinlet_temperature_k_per_k"),
        n(fan, "dtemperature_dspeed_ratio_k"), 0.0, 0.0,
    ];
    if radiative {
        let surface = &adjoint.get("radiation").unwrap()
            .get("surfaces").unwrap().as_array().unwrap()[0];
        values[5] = n(surface, "dtemperature_demissivity_k");
        values[6] = n(surface, "dtemperature_dreservoir_temperature_k_per_k");
    }
    values
}

/// Independent uniform BE recurrence with the exact exponential air law.
/// Single-bank affinity and fixed quadratic losses give Q proportional to
/// absolute speed. The Gnielinski formula supplies its analytic h elasticity.
fn air_trajectory_derivatives(rows: &[J], speed: f64, reservoir: Option<f64>) -> [f64; 7] {
    let mut old = INITIAL_K;
    let mut derivative = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for row in rows {
        let exchange = row.get("conjugate").unwrap().get("exchange").unwrap();
        let segment = &exchange.get("segments").unwrap().as_array().unwrap()[0];
        let re = n(segment, "reynolds");
        let z = 0.79 * re.ln() - 1.64;
        let a = 1.0 / (8.0 * z * z);
        let correction = 12.7 * a.sqrt() * (0.707_f64.powf(2.0 / 3.0) - 1.0);
        let h = a * (re - 1000.0) * 0.707 / (1.0 + correction) * 0.0263 / 0.02;
        close(h, n(segment, "htc_w_m2_k"), 1e-10 * h);
        let capacity_rate = n(exchange, "mass_flow_kg_s") * AIR_CP;
        let ntu = h * area_m2() / capacity_rate;
        let conductance = capacity_rate * -(-ntu).exp_m1();
        let friction_elasticity = -1.58 / z;
        let h_elasticity = friction_elasticity + re / (re - 1000.0)
            - 0.5 * correction / (1.0 + correction) * friction_elasticity;
        let conductance_speed = (conductance
            + capacity_rate * (-ntu).exp() * ntu * (h_elasticity - 1.0)) / speed;
        let dt = n(row, "dt_s");
        let storage = CAPACITY_J_K / dt;
        let radiation = |temperature: f64| reservoir.map_or(0.0, |ambient|
            EMISSIVITY * SIGMA * area_m2() * (temperature.powi(4) - ambient.powi(4)));
        let temperature = bisect(200.0, 450.0, |temperature|
            storage * (temperature - old) + conductance * (temperature - INLET_K)
                + radiation(temperature) - 5.0);
        close(temperature, n(row, "final_region_max_k"), 4e-8);
        let jacobian = storage + conductance + reservoir.map_or(0.0, |_|
            4.0 * EMISSIVITY * SIGMA * area_m2() * temperature.powi(3));
        let direct = [
            0.0, 1.0, -storage * (temperature - old), conductance,
            -conductance_speed * (temperature - INLET_K),
            reservoir.map_or(0.0, |ambient|
                -SIGMA * area_m2() * (temperature.powi(4) - ambient.powi(4))),
            reservoir.map_or(0.0, |ambient|
                4.0 * EMISSIVITY * SIGMA * area_m2() * ambient.powi(3)),
        ];
        for (gradient, rhs) in derivative.iter_mut().zip(direct) {
            *gradient = (storage * *gradient + rhs) / jacobian;
        }
        old = temperature;
    }
    derivative
}

#[test]
fn g3_native_capacity_air_adjoint_matches_implicit_trajectory_and_physical_controls() {
    let speed = 1.2; // The project control is absolute speed, not log(speed).
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for reservoir in [None, Some(350.0)] {
        let (mut spec, bytes) = regular_project();
        set_air_speed(&mut spec, speed);
        if let Some(ambient) = reservoir {
            radiation_product::declare(&mut spec, &card, "air", ambient);
        }
        let plain = run(&spec, &cards, &bytes);
        let mut requested = spec.clone();
        request_air_adjoint(&mut requested, "temperature-final-adjoint");
        let solved = run(&requested, &cards, &bytes);
        assert_eq!(plain.temperature, solved.temperature);
        let adjoint = solved.transient().get("adjoint").unwrap();
        close(n(adjoint, "value_k"), final_k(&plain), 0.0);
        close(n(adjoint, "reconstructed_solid_endpoints"), 0.0, 0.0);
        close(n(adjoint, "reconstruction_krylov_iterations"), 0.0, 0.0);
        close(n(adjoint, "revalidated_solid_endpoints"), solved.rows("fine").len() as f64, 0.0);
        let actual = air_gradients(adjoint, reservoir.is_some());
        let expected = air_trajectory_derivatives(solved.rows("fine"), speed, reservoir);
        for (&actual, &expected) in actual.iter().zip(&expected) {
            close(actual, expected, 4e-7 * expected.abs().max(1.0));
        }
        assert!(actual[0] > 0.0 && actual[0] < 1.0);
        assert!(actual[1] > 0.0 && actual[3] > 0.0 && actual[4] < 0.0);
        let count = if reservoir.is_some() { 7 } else { 5 };
        for (control, &derivative) in actual[..count].iter().enumerate() {
            let epsilon = if control == 4 { 0.0005 } else { 0.001 };
            let evaluate = |delta: f64| {
                let mut probe = spec.clone();
                let mut emissivity = EMISSIVITY;
                match control {
                    0 => probe.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .transient.as_mut().unwrap().initial_temperature =
                        fs_project::TransientInitialTemperature::Uniform(
                            QtyAny::new(INITIAL_K + delta, dims::TEMPERATURE)),
                    1 => probe.power.as_mut().unwrap()[0].watts.value += delta,
                    2 => {
                        for row in &mut probe.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                            .transient.as_mut().unwrap().capacities {
                            row.volumetric_heat_capacity.value *= 1.0 + delta;
                        }
                    },
                    3 => set_air_inlet(&mut probe, INLET_K + delta),
                    4 => set_air_speed(&mut probe, speed + delta),
                    5 => emissivity += delta,
                    6 => probe.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += delta,
                    _ => unreachable!(),
                }
                let (probe_cards, probe_card) =
                    radiation_product::emissivity_cards(&fixture_cards(), emissivity);
                if reservoir.is_some() {
                    probe.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .radiation.as_mut().unwrap().surfaces[0].card = probe_card;
                }
                assert_eq!(probe.envelope, spec.envelope);
                final_k(&run(&probe, &probe_cards, &bytes))
            };
            close(derivative, (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon),
                3e-5 * derivative.abs().max(1.0));
        }
    }
}

#[test]
fn g1_native_capacity_air_initial_peak_has_zero_future_controls() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let (mut spec, bytes) = regular_project();
    radiation_product::declare(&mut spec, &card, "air", 270.0);
    request_air_adjoint(&mut spec, "temperature-sampled-peak-adjoint");
    let solved = run(&spec, &cards, &bytes);
    let adjoint = solved.transient().get("adjoint").unwrap();
    close(n(adjoint, "state_index"), 0.0, 0.0);
    close(n(adjoint, "value_k"), INITIAL_K, 0.0);
    close(n(adjoint, "adjoint_krylov_iterations"), 0.0, 0.0);
    assert_eq!(air_gradients(adjoint, true), [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    for interval in adjoint.get("intervals").unwrap().as_array().unwrap() {
        close(n(interval, "dtemperature_dsource_multiplier_k"), 0.0, 0.0);
        for region in interval.get("regional_power").unwrap().as_array().unwrap() {
            close(n(region, "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
        }
    }
}


fn mixed_adjoint_project(card: &str) -> (ProjectSpec, Vec<u8>) {
    let bytes = tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    transient_product::declaration(&mut spec, 100.0, 0.5);
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap().initial_temperature =
        fs_project::TransientInitialTemperature::Uniform(
            QtyAny::new(INITIAL_K, dims::TEMPERATURE));
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    for (name, normal) in [
        ("upstream-wall", [1.0, 0.0, 0.0]),
        ("natural-wall", [0.0, 0.0, 1.0]),
    ] {
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".into(), name: name.into(), display: name.into(), expect_id: None,
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
    let first = setup.boundaries[0].condition.clone();
    if let ThermalBoundaryCondition::AirflowConvection { order, .. } =
        &mut setup.boundaries[0].condition { *order = 1; }
    setup.boundaries.push(ThermalBoundary { target: "upstream-wall".into(), condition: first });
    setup.boundaries.push(ThermalBoundary {
        target: "natural-wall".into(),
        condition: ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: QtyAny::new(0.1, dims::LENGTH),
            ambient_temperature: QtyAny::new(350.0, dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        },
    });
    // Radiation on the downstream air port exercises its combined Robin
    // reference's off-diagonal response to the upstream wall.
    radiation_product::declare(&mut spec, card, "air", 270.0);
    let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap();
    radiation.max_iterations = 24;
    radiation.relaxation = 1.0;
    set_air_speed(&mut spec, 1.2);
    (spec, bytes)
}

#[test]
fn g3_native_capacity_mixed_adjoint_keeps_upstream_air_and_independent_reservoir_controls() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let (spec, bytes) = mixed_adjoint_project(&card);
    let plain = run(&spec, &cards, &bytes);
    let mut requested = spec.clone();
    request_air_adjoint(&mut requested, "temperature-final-adjoint");
    let solved = run(&requested, &cards, &bytes);
    assert_eq!(plain.temperature, solved.temperature);
    let minimum = solved.temperature.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(final_k(&solved) - minimum > 0.01);
    let air = solved.receipt.get("conjugate").unwrap();
    let segments = air.get("segments").unwrap().as_array().unwrap();
    assert_eq!(segments.len(), 2);
    assert_eq!(segments[0].str_field("target"), Some("upstream-wall"));
    assert_eq!(segments[1].str_field("target"), Some("air"));
    close(n(&segments[0], "air_out_k"), n(&segments[1], "air_in_k"), 1e-12);
    assert!(n(&segments[1], "air_in_k") > INLET_K + 0.1);
    assert!(n(solved.receipt.get("radiation").unwrap(), "nonlinear_radiative_out_w") > 0.0);
    let natural = &solved.receipt.get("natural_convection").unwrap()
        .get("laws").unwrap().as_array().unwrap()[0];
    assert!(n(natural, "heat_rate_w") < 0.0);
    let adjoint = solved.transient().get("adjoint").unwrap();
    let gradients = air_gradients(adjoint, true);
    let ambient = air_control(adjoint, "natural_reservoirs", "target", "natural-wall");
    let direction = 2.0 * gradients[0] + 3.0 * gradients[1] + 0.04 * gradients[2]
        + 0.3 * gradients[3] + 0.02 * gradients[4] + 0.015 * gradients[5]
        + 0.4 * gradients[6] + 0.5 * n(ambient, "dtemperature_dambient_temperature_k_per_k");
    let evaluate = |delta: f64| {
        let (probe_cards, probe_card) =
            radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY + 0.015 * delta);
        let mut probe = mixed_adjoint_project(&probe_card).0;
        let setup = probe.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        let time = setup.transient.as_mut().unwrap();
        time.initial_temperature = fs_project::TransientInitialTemperature::Uniform(
            QtyAny::new(INITIAL_K + 2.0 * delta, dims::TEMPERATURE));
        for row in &mut time.capacities { row.volumetric_heat_capacity.value *= 1.0 + 0.04 * delta; }
        let boundary = setup.boundaries.iter_mut().find(|row| row.target == "natural-wall").unwrap();
        let ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } =
            &mut boundary.condition else { panic!("expected natural wall") };
        ambient_temperature.value += 0.5 * delta;
        setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += 0.4 * delta;
        probe.power.as_mut().unwrap()[0].watts.value += 3.0 * delta;
        set_air_inlet(&mut probe, INLET_K + 0.3 * delta);
        set_air_speed(&mut probe, 1.2 + 0.02 * delta);
        assert_eq!(probe.envelope, spec.envelope);
        final_k(&run(&probe, &probe_cards, &bytes))
    };
    let epsilon = 0.01;
    close(direction, (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon),
        3e-5 * direction.abs().max(1.0));
}
