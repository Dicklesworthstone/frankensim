//! G1/G3 mixed natural/forced endpoint physics through the native pipeline.
//! A unit tetra has three prescribed vertices and one free temperature.
//! The scalar oracle integrates its P1 Robin rows and independent air law;
//! every grid must apply the initial prescribed-temperature storage only once.

use super::*;
use fs_project::spec::dims;

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const INITIAL_K: f64 = 300.0;
const FIXED_K: f64 = 330.0;
const INLET_K: f64 = 315.0;
const NODE_CAPACITY_J_K: f64 = 0.625;
const RADIATION_K: f64 = 295.0;
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

// Independent Churchill-Chu evaluation with the declared native air properties.
fn natural_h(wall: f64, ambient: f64) -> f64 {
    let film = 0.5 * (wall + ambient);
    let density = 101325.0 / (287.05 * film);
    let nu = 1.846e-5 / density;
    let ra = 9.80665 * (wall - ambient).abs() * 0.1_f64.powi(3)
        / (film * nu * (nu / 0.707));
    let nu_number = (0.825 + 0.387 * ra.powf(1.0 / 6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0)).powi(2);
    nu_number * 0.0263 / 0.1
}

fn bisect(residual: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (280.0, 350.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let trial = 0.5 * (low + high);
        if residual(trial) > 0.0 { high = trial; } else { low = trial; }
    }
    0.5 * (low + high)
}

fn project(bytes: &[u8], ambient: f64) -> ProjectSpec {
    let mut spec = transient_product::driven_face_project(bytes);
    let template = conjugate_fixture_project(7, bytes, "convection.gnielinski")
        .cooling.unwrap().conduction.unwrap().boundaries[0].condition.clone();
    for axis in 0..3 {
        let target = if axis == 0 { "natural-wall".to_string() }
            else { format!("air-face-{axis}") };
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".to_string(), name: target.clone(),
            display: "Manufactured mixed-convection coordinate face".to_string(),
            expect_id: None,
        });
        let mut normal = [0.0; 3];
        normal[axis] = 1.0;
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".to_string(), target: target.clone(),
            length_unit: "m".to_string(),
            selector: MeshSelector::HalfSpace {
                normal, offset: 0.0, side: HalfSpaceSide::AtMost, tolerance: 0.0,
            },
            allow_overlap: true,
        });
        let condition = if axis == 0 {
            ThermalBoundaryCondition::NaturalConvection {
                characteristic_length: QtyAny::new(0.1, dims::LENGTH),
                ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
                correlation: "convection.churchill-chu-vertical-plate".to_string(),
            }
        } else {
            let mut condition = template.clone();
            if let ThermalBoundaryCondition::AirflowConvection {
                order, inlet_temperature, ..
            } = &mut condition {
                *order = (axis - 1) as u32;
                inlet_temperature.value = INLET_K;
            }
            condition
        };
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .boundaries.push(ThermalBoundary { target, condition });
    }
    spec.envelope.as_mut().unwrap().ambient_lo.value = INLET_K;
    spec.envelope.as_mut().unwrap().ambient_hi.value = INLET_K;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    spec
}

fn run(spec: &ProjectSpec, cards: &CardPackSet, bytes: &[u8]) -> (String, J, J) {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, bytes.to_vec());
    let (id, _, receipt, field) = run_conjugate_to_completion(&ledger, &decode(spec), cards);
    (id, J::parse(&receipt).unwrap(), J::parse(&field).unwrap())
}

fn check_run(receipt: &J, field: &J, ambient: f64, nonlinear: bool, radiation: bool) -> f64 {
    let transient = receipt.get("transient").unwrap();
    let mut final_free = f64::NAN;
    let mut all_air_iterations = 0.0;
    let mut all_solid_solves = 0.0;
    for grid in ["coarse", "fine"] {
        let mut old_expected = INITIAL_K;
        let mut old_actual = INITIAL_K;
        for (index, row) in transient.get(grid).unwrap().as_array().unwrap().iter().enumerate() {
            let natural = row.get("natural").unwrap();
            let law = &natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap()[0];
            let air = row.get("conjugate").unwrap();
            let exchange = air.get("exchange").unwrap();
            let segments = exchange.get("segments").unwrap().as_array().unwrap();
            assert_eq!(segments.len(), 2);
            let h_air = n(&segments[0], "htc_w_m2_k");
            let mcp = n(exchange, "mass_flow_kg_s") * 1007.0;
            let conductance = mcp * -(-h_air / mcp).exp_m1();
            for (order, segment) in segments.iter().enumerate() {
                close(n(segment, "wetted_area_m2"), 0.5, 1e-12);
                close(n(segment, "htc_w_m2_k"), h_air, 1e-12);
                close(n(segment, "order"), order as f64, 0.0);
                assert!(matches!(segment.get("in_domain"), Some(J::Bool(true))));
            }
            let dt = n(row, "dt_s");
            let expected = bisect(|t| {
                let wall = (t + 2.0 * FIXED_K) / 3.0;
                let h_natural = natural_h(wall, ambient);
                let h_rad = if radiation {
                    EMISSIVITY * SIGMA * (wall + RADIATION_K)
                        * (wall * wall + RADIATION_K * RADIATION_K)
                } else { 0.0 };
                let k = if nonlinear { 0.02 * ((t + 3.0 * FIXED_K) / 4.0) - 4.0 }
                    else { 2.0 };
                let natural_w = 0.5 * h_natural * (wall - ambient);
                let air_w = conductance * (wall - INLET_K);
                let radiation_w = 0.5 * h_rad * (wall - RADIATION_K);
                // Each area-1/2 face contributes Q/3 + h*(T-Tfixed)/36
                // to the free P1 row. Two air faces share the same wall mean
                // but have different marched air reference temperatures.
                NODE_CAPACITY_J_K * (t - old_expected) / dt
                    + 0.5 * k * (t - FIXED_K)
                    + (h_natural + 2.0 * h_air + h_rad) * (t - FIXED_K) / 36.0
                    + (natural_w + air_w + radiation_w) / 3.0
            });
            let clamp_change = if index == 0 { 3.0 * (FIXED_K - INITIAL_K) } else { 0.0 };
            let stored = n(row, "stored_energy_change_j");
            let actual = old_actual + stored / NODE_CAPACITY_J_K - clamp_change;
            close(actual, expected, 4e-7);
            let wall = (actual + 2.0 * FIXED_K) / 3.0;
            let h_natural = natural_h(wall, ambient);
            let natural_w = 0.5 * h_natural * (wall - ambient);
            let air_w = conductance * (wall - INLET_K);
            let radiation_w = if radiation {
                0.5 * EMISSIVITY * SIGMA * (wall.powi(4) - RADIATION_K.powi(4))
            } else { 0.0 };
            close(n(law, "mean_wall_k"), wall, 2e-7);
            close(n(law, "htc_w_m2_k"), h_natural, 2e-8);
            close(n(law, "heat_rate_w"), natural_w, 2e-6);
            assert_eq!(natural_w.is_sign_positive(), ambient < wall);
            assert!(air_w > 0.0);
            close(n(air, "air_heat_gain_w"), air_w, 3e-6);
            close(n(air, "off_path_convective_out_w"), natural_w, 3e-6);
            close(n(exchange, "outlet_k"), INLET_K + air_w / mcp, 2e-7);
            close(n(air, "radiative_out_w"), radiation_w, 3e-6);
            close(n(natural, "physical_radiative_out_w"), radiation_w, 3e-6);
            close(n(natural, "physical_convective_out_w"), natural_w + air_w, 3e-6);
            close(n(natural, "physical_dirichlet_in_w"),
                stored / dt + natural_w + air_w + radiation_w, 1e-5);
            close(n(natural, "coupled_net_input_w"), stored / dt, 4e-6);
            assert!(n(natural, "physical_residual_norm_j")
                <= n(natural, "physical_residual_tolerance_j"));
            assert!(n(natural, "physical_energy_residual_j").abs() <= 1e-6);
            assert!(n(natural, "coupled_energy_residual_j").abs() <= 1e-6);
            assert!(n(natural, "iterations") > 1.0 && n(natural, "iterations") <= 80.0);
            assert!(n(natural, "airflow_iterations") >= n(exchange, "iterations"));
            assert!(n(natural, "solid_solves") >= n(air, "solid_solves"));
            close(n(row, "krylov_iterations"), n(natural, "krylov_iterations"), 0.0);
            close(n(natural, "endpoint_evaluations"), 2.0 * n(natural, "iterations") - 1.0, 0.0);
            if nonlinear { assert!(n(natural, "nonlinear_iterations") > 0.0); }
            if radiation { assert!(n(natural, "radiation_trials") > 0.0); }
            all_air_iterations += n(natural, "airflow_iterations");
            all_solid_solves += n(natural, "solid_solves");
            old_expected = expected;
            old_actual = actual;
        }
        if grid == "fine" { final_free = old_actual; }
    }
    let mut temperatures = field.get("temperature").unwrap().as_array().unwrap()
        .iter().map(|value| value.as_f64().unwrap()).collect::<Vec<_>>();
    temperatures.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(temperatures.len(), 4);
    close(temperatures[0], final_free, 4e-7);
    for &fixed in &temperatures[1..] { close(fixed, FIXED_K, 0.0); }
    let natural = transient.get("natural").unwrap();
    let conjugate = transient.get("conjugate").unwrap();
    close(n(conjugate, "total_iterations"), all_air_iterations, 0.0);
    close(n(conjugate, "total_solid_solves"), all_solid_solves, 0.0);
    close(n(natural, "total_solid_solves"), all_solid_solves, 0.0);
    assert!(n(conjugate, "maximum_coupled_energy_residual_j") <= 1e-6);
    close(n(natural.get("energy").unwrap(), "stored_change_j"),
        NODE_CAPACITY_J_K * (final_free - INITIAL_K + 3.0 * (FIXED_K - INITIAL_K)), 2e-6);
    assert!(n(natural.get("energy").unwrap(), "window_residual_j").abs() <= 3e-6);
    final_free
}

#[test]
fn g1_native_transient_mixed_convection_preserves_signed_reservoirs_and_history() {
    let bytes = tetra_stl();
    let cards = transient_product::transient_conductivity_cards("constant-scalar");
    for ambient in [300.0, 340.0] {
        let mut spec = project(&bytes, ambient);
        rebind_to(&mut spec, &cards);
        let solved = run(&spec, &cards, &bytes);
        check_run(&solved.1, &solved.2, ambient, false, false);
        let replayed = run(&spec, &cards, &bytes);
        assert_eq!(solved, replayed, "the same mixed transient replays its fields and work");
    }
}

#[test]
fn g1_native_transient_mixed_convection_composes_endpoint_conductivity_and_radiation() {
    let bytes = tetra_stl();
    let mut final_temperatures = Vec::new();
    for kind in ["nonlinear", "constant-scalar"] {
        let base = transient_product::transient_conductivity_cards(kind);
        let (cards, card) = radiation_product::emissivity_cards(&base, EMISSIVITY);
        let mut spec = project(&bytes, 340.0);
        rebind_to(&mut spec, &base);
        radiation_product::declare(&mut spec, &card, "natural-wall", RADIATION_K);
        let (_, receipt, field) = run(&spec, &cards, &bytes);
        final_temperatures.push(check_run(&receipt, &field, 340.0, kind == "nonlinear", true));
    }
    assert!((final_temperatures[0] - final_temperatures[1]).abs() > 0.01,
        "freezing k(T) must change this nonuniform coupled endpoint");
}

#[test]
fn g0_native_transient_mixed_failures_publish_no_partial_endpoint() {
    let bytes = tetra_stl();
    for failure in ["air-card-domain", "radiation-work", "transient-adjoint"] {
        let base = transient_product::transient_conductivity_cards("constant-scalar");
        let (cards, card) = radiation_product::emissivity_cards(&base, EMISSIVITY);
        let mut spec = project(&bytes, 340.0);
        rebind_to(&mut spec, &base);
        let expected = match failure {
            "air-card-domain" => {
                for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
                    if let ThermalBoundaryCondition::AirflowConvection { channel_length, .. } = &mut boundary.condition {
                        channel_length.value = 0.1;
                    }
                }
                "cli-solve-conduction-airflow-card-domain"
            }
            "radiation-work" => {
                radiation_product::declare(&mut spec, &card, "natural-wall", RADIATION_K);
                spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .radiation.as_mut().unwrap().max_iterations = 1;
                "cli-solve-conduction-transient"
            }
            _ => {
                spec.outputs.as_mut().unwrap().push(OutputRequest {
                    name: "temperature-final-adjoint".to_string(),
                    kind: "report".to_string(), region: None,
                });
                "cli-solve-conduction-transient"
            }
        };
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes.clone());
        let error = run_solve(
            &ledger, &CancelGate::new_clock_free(), &mut benign_clock(),
            &decode(&spec), &cards, &mut Vec::new(),
        ).unwrap_err();
        assert_eq!(error.code, expected, "{failure}: {error:?}");
        assert_eq!(artifacts_of_kind(&ledger, "solve-conduction-solution"), 0,
            "{failure} must preserve only the completed predecessor stages");
    }
}
