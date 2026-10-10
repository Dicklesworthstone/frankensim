//! Native imported-mesh transient radiation with independent discrete balances.
//! The regular tetra admits an exact uniform P1 mode. The driven unit tetra
//! reduces to one free temperature, including the actual Robin trace weights.

use super::*;
use fs_project::{ConductionTransient, RadiatingSurface, spec::dims};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITY: f64 = 0.85;
const AMBIENT_K: f64 = 293.15;
const INITIAL_K: f64 = 330.0;
const CAPACITY_J_K: f64 = 100.0 * 8.0 / 3.0;

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

fn regular_tetra_stl() -> Vec<u8> {
    // Integer coordinates survive the STL importer's f32 conversion exactly.
    // Edge length is sqrt(8), total area 8*sqrt(3), volume 8/3.
    let p = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let mut stl = String::from("solid enclosure\n");
    for [a, b, c] in [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]] {
        stl.push_str(&facet(p[a], p[b], p[c]));
    }
    stl.push_str("endsolid enclosure\n");
    stl.into_bytes()
}

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
}

fn regular_project(card: &str, reservoir: f64) -> (ProjectSpec, Vec<u8>) {
    let bytes = regular_tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    transient_product::declaration(&mut spec, 100.0, 1.0);
    time(&mut spec).initial_temperature.value = INITIAL_K;
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .regions[0]
        .seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    radiation_product::declare(&mut spec, card, "air", reservoir);
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

    fn final_k(&self) -> f64 {
        n(&self.qoi.get("qoi").unwrap().as_array().unwrap()[0], "value")
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

fn radiative_w(temperature: f64, reservoir: f64, area: f64) -> f64 {
    EMISSIVITY * SIGMA * area * (temperature.powi(4) - reservoir.powi(4))
}

fn uniform_step(old: f64, dt: f64, reservoir: f64, radiation: bool) -> f64 {
    bisect(200.0, 450.0, |t| {
        CAPACITY_J_K * (t - old) / dt
            + 10.0 * area_m2() * (t - AMBIENT_K)
            + if radiation { radiative_w(t, reservoir, area_m2()) } else { 0.0 }
            - 5.0
    })
}

fn checked_radiation(row: &J) -> &J {
    let radiation = row.get("radiation").expect("retain per-step radiation evidence");
    assert!(n(radiation, "iterations") > 1.0, "fixture must require nonlinear trials");
    assert!(n(radiation, "iterations") <= 128.0);
    assert!(
        n(radiation, "physical_residual_norm_j")
            <= n(radiation, "physical_residual_tolerance_j")
    );
    assert!(
        n(radiation, "physical_energy_residual_j").abs() <= 1e-6,
        "the declared energy gate is not widened to a radiation watt tolerance"
    );
    radiation
}

#[test]
fn g1_native_transient_radiation_matches_each_endpoint_without_advancing_history_per_trial() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for reservoir in [270.0, 350.0] {
        let (spec, bytes) = regular_project(&card, reservoir);
        let solved = run(&spec, &cards, &bytes);
        assert_eq!(solved.temperature.len(), 4, "oracle requires one imported tetra");
        close(n(solved.transient(), "total_steps"), 6.0, 0.0);
        assert!(!matches!(solved.transient().get("radiation"), None | Some(J::Null)));
        for grid in ["coarse", "fine"] {
            let mut old_expected = INITIAL_K;
            let mut old_actual = INITIAL_K;
            let mut previous_time = 0.0;
            for row in solved.rows(grid) {
                let dt = n(row, "dt_s");
                close(n(row, "time_s") - previous_time, dt, 1e-15);
                let expected = uniform_step(old_expected, dt, reservoir, true);
                let actual = n(row, "final_region_max_k");
                close(actual, expected, 2e-8);
                let radiation = checked_radiation(row);
                let convective = 10.0 * area_m2() * (actual - AMBIENT_K);
                let nonlinear = radiative_w(actual, reservoir, area_m2());
                let storage = CAPACITY_J_K * (actual - old_actual);
                close(n(row, "source_w"), 5.0, 1e-12);
                close(n(row, "stored_energy_change_j"), storage, 2e-6);
                close(n(radiation, "convective_out_w"), convective, 2e-6);
                close(n(radiation, "nonlinear_radiative_out_w"), nonlinear, 2e-6);
                close(n(radiation, "physical_dirichlet_in_w"), 0.0, 1e-12);
                close(
                    n(radiation, "radiative_out_w"),
                    nonlinear,
                    n(radiation, "max_heat_mismatch_w") + 2e-6,
                );
                let physical_residual = storage - dt * (5.0 - convective - nonlinear);
                close(physical_residual, 0.0, 2e-6);
                close(
                    n(radiation, "physical_energy_residual_j"),
                    physical_residual,
                    2e-6,
                );
                old_expected = expected;
                old_actual = actual;
                previous_time = n(row, "time_s");
            }
            close(previous_time, 1.0, 0.0);
        }
        let last = solved.rows("fine").last().unwrap();
        let final_k = n(last, "final_region_max_k");
        for &temperature in &solved.temperature {
            close(temperature, final_k, 2e-8);
        }
        close(solved.final_k(), final_k, 0.0);
        let scope = solved.qoi.get("thermal_time_scope").unwrap();
        assert_eq!(scope.str_field("kind"), Some("final"));
        close(n(scope, "time_s"), 1.0, 0.0);
        close(
            n(solved.transient().get("energy").unwrap(), "stored_change_j"),
            CAPACITY_J_K * (final_k - INITIAL_K),
            2e-6,
        );
        let endpoint_radiation = solved.receipt.get("radiation").unwrap();
        let last_radiation = last.get("radiation").unwrap();
        for key in ["convective_out_w", "radiative_out_w", "nonlinear_radiative_out_w"] {
            close(n(endpoint_radiation, key), n(last_radiation, key), 1e-12);
        }
        close(
            n(endpoint_radiation, "solid_solves"),
            n(last_radiation, "iterations"),
            0.0,
        );
        let surfaces = endpoint_radiation.get("surfaces").unwrap().as_array().unwrap();
        assert_eq!(surfaces.len(), 1);
        assert_eq!(surfaces[0].str_field("card"), Some(card.as_str()));
        close(n(&surfaces[0], "area_m2"), area_m2(), 2e-12);
        close(n(&surfaces[0], "mean_temperature_k"), final_k, 2e-8);
        assert_eq!(
            n(endpoint_radiation, "nonlinear_radiative_out_w").is_sign_positive(),
            reservoir < final_k,
            "outward and inward exchange retain their physical signs"
        );
        close(
            n(solved.receipt.get("energy").unwrap(), "storage_w"),
            n(last, "stored_energy_change_j") / n(last, "dt_s"),
            2e-9,
        );
    }
}

#[test]
fn g1_native_radiation_changes_the_transient_and_preserves_the_steady_model() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let (spec, bytes) = regular_project(&card, AMBIENT_K);
    let on = run(&spec, &cards, &bytes);
    let mut off_spec = spec.clone();
    off_spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().radiation = None;
    let off = run(&off_spec, &cards, &bytes);
    let mut expected_off = INITIAL_K;
    for row in off.rows("fine") {
        expected_off = uniform_step(expected_off, n(row, "dt_s"), AMBIENT_K, false);
    }
    close(off.final_k(), expected_off, 2e-8);
    assert!(on.final_k() < off.final_k() - 0.1);
    assert!(matches!(off.receipt.get("radiation"), None | Some(J::Null)));

    let mut steady_spec = spec;
    steady_spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient = None;
    let steady = run(&steady_spec, &cards, &bytes);
    let expected_steady = bisect(AMBIENT_K, INITIAL_K, |t| {
        10.0 * area_m2() * (t - AMBIENT_K)
            + radiative_w(t, AMBIENT_K, area_m2())
            - 5.0
    });
    close(steady.final_k(), expected_steady, 2e-8);
    assert!(steady.final_k() < on.final_k());
    assert!(matches!(steady.receipt.get("transient"), None | Some(J::Null)));
    let radiation = steady.receipt.get("radiation").unwrap();
    close(
        n(radiation, "convective_out_w") + n(radiation, "nonlinear_radiative_out_w"),
        5.0,
        2e-6,
    );
}

fn driven_project(card: &str, bytes: &[u8]) -> ProjectSpec {
    let mut spec = transient_product::driven_face_project(bytes);
    let mut targets = Vec::new();
    for axis in 0..3 {
        let target = format!("radiating-face-{axis}");
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".to_string(),
            name: target.clone(),
            display: "Radiating coordinate face".to_string(),
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
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries.push(
            ThermalBoundary {
                target: target.clone(),
                condition: ThermalBoundaryCondition::Convection {
                    coefficient: QtyAny::new(1.0, dims::HEAT_TRANSFER_COEFFICIENT),
                    reference_temperature: QtyAny::new(330.0, dims::TEMPERATURE),
                },
            },
        );
        targets.push(target);
    }
    radiation_product::declare(&mut spec, card, &targets[0], 330.0);
    let policy = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap();
    let template = policy.surfaces[0].clone();
    policy.surfaces = targets.into_iter().map(|target| RadiatingSurface {
        name: target.clone(),
        target,
        ..template.clone()
    }).collect();
    spec
}

fn driven_step(old: f64, dt: f64, nonlinear: bool) -> f64 {
    // Each row-sum capacity is 15/24 J/K. The free conduction row is
    // k((T+3G)/4)*(T-G)/2. Each coordinate face has area 1/2, and the sum
    // of its P1 Robin trace rows is (h+h_rad)*(T-G)/4 for G=Tair=Trad=330.
    // The three equal patch means are (T+2G)/3, not the free-node temperature.
    bisect(300.0, 330.0, |t| {
        let mean = (t + 2.0 * 330.0) / 3.0;
        let h_rad = EMISSIVITY * SIGMA * (mean + 330.0) * (mean * mean + 330.0 * 330.0);
        let k = if nonlinear { 0.02 * ((t + 3.0 * 330.0) / 4.0) - 4.0 } else { 2.0 };
        0.625 * (t - old) + dt * (0.5 * k + 0.25 * (1.0 + h_rad)) * (t - 330.0)
    })
}

#[test]
fn g1_native_transient_radiation_and_nonlinear_conductivity_share_the_physical_endpoint() {
    let mut free_temperatures = Vec::new();
    for kind in ["nonlinear", "constant-scalar"] {
        let base = transient_product::transient_conductivity_cards(kind);
        let (cards, card) = radiation_product::emissivity_cards(&base, EMISSIVITY);
        let bytes = tetra_stl();
        let mut spec = driven_project(&card, &bytes);
        rebind_to(&mut spec, &base);
        let solved = run(&spec, &cards, &bytes);
        let mut temperature = solved.temperature.clone();
        temperature.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(temperature.len(), 4);
        for &fixed in &temperature[1..] {
            close(fixed, 330.0, 0.0);
        }
        for grid in ["coarse", "fine"] {
            let mut old = 300.0;
            for (index, row) in solved.rows(grid).iter().enumerate() {
                let dt = n(row, "dt_s");
                let expected = driven_step(old, dt, kind == "nonlinear");
                let stored = 0.625 * (expected - old + if index == 0 { 90.0 } else { 0.0 });
                let radiation = checked_radiation(row);
                let mean = (expected + 2.0 * 330.0) / 3.0;
                let convective = 1.5 * (mean - 330.0);
                let radiative = radiative_w(mean, 330.0, 1.5);
                close(n(row, "stored_energy_change_j"), stored, 2e-6);
                close(n(radiation, "convective_out_w"), convective, 2e-6);
                close(n(radiation, "nonlinear_radiative_out_w"), radiative, 2e-6);
                close(
                    n(radiation, "physical_dirichlet_in_w"),
                    stored / dt + convective + radiative,
                    1e-5,
                );
                old = expected;
            }
            if grid == "fine" {
                close(temperature[0], old, 2e-7);
            }
        }
        if kind == "nonlinear" {
            let nonlinear = solved.transient().get("nonlinear").unwrap();
            assert_eq!(nonlinear.str_field("method"), Some("newton-fgmres"));
            assert!(n(nonlinear, "total_updates") > 0.0);
        }
        free_temperatures.push(temperature[0]);
    }
    assert!(
        (free_temperatures[0] - free_temperatures[1]).abs() > 0.5,
        "freezing k at 300 K must change the radiating, nonuniform transient"
    );
}

#[test]
fn g0_native_transient_radiation_refuses_exhausted_trials_and_surface_domain_escape() {
    for failure in ["outer-iteration-cap", "surface-temperature-domain"] {
        let (cards, card) = if failure == "surface-temperature-domain" {
            // Query and reservoir are valid at 300 K. The first physical
            // endpoint remains above 315 K, so only its actual trace is invalid.
            radiation_product::emissivity_cards_with_domain(
                &fixture_cards(), EMISSIVITY, 295.0, 315.0,
            )
        } else {
            radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY)
        };
        let (mut spec, bytes) = regular_project(&card, 300.0);
        if failure == "outer-iteration-cap" {
            spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .radiation.as_mut().unwrap().max_iterations = 1;
        }
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
        assert_eq!(error.code, "cli-solve-conduction-transient", "{failure}: {error:?}");
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4,
            "{failure} cannot publish a partial field or final QoI"
        );
    }
}
