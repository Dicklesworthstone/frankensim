//! Native transient natural convection with independent physical balances.
//! Uniform regular-tetra modes and a one-free-node driven tetra exercise the
//! real imported mesh, card, storage, nonlinear material and radiation paths.

use super::*;
use fs_project::{RadiatingSurface, spec::dims};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const AMBIENT_K: f64 = 293.15;
const LENGTH_M: f64 = 0.1;
const PRESSURE_PA: f64 = 101325.0;
const CAPACITY_J_K: f64 = 100.0 * 8.0 / 3.0;
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

/// Independent Churchill-Chu formula with the native model's stated frozen
/// transport properties and ideal-gas film density. No production correlation
/// evaluator supplies this oracle. Its zero-difference limit is used only for
/// bisection brackets; every accepted state is checked against the Ra domain.
fn coefficient(wall: f64, ambient: f64, length: f64) -> (f64, f64, f64) {
    let film = 0.5 * (wall + ambient);
    let density = PRESSURE_PA / (287.05 * film);
    let kinematic = 1.846e-5 / density;
    let diffusivity = kinematic / 0.707;
    let rayleigh = 9.80665 * (wall - ambient).abs() * length.powi(3)
        / (film * kinematic * diffusivity);
    let factor = 0.825
        + 0.387 * rayleigh.powf(1.0 / 6.0)
            / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let nusselt = factor * factor;
    (nusselt * 0.0263 / length, rayleigh, nusselt)
}

fn natural_law(ambient: f64) -> ThermalBoundaryCondition {
    ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(LENGTH_M, dims::LENGTH),
        ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".to_string(),
    }
}

fn regular_project(initial: f64, ambient: f64, watts: f64) -> (ProjectSpec, Vec<u8>) {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    transient_product::declaration(&mut spec, 100.0, 1.0);
    spec.power.as_mut().unwrap()[0].watts.value = watts;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient.as_mut().unwrap().initial_temperature.value = initial;
    setup.regions[0].seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    setup.boundaries[0].condition = natural_law(ambient);
    spec.envelope.as_mut().unwrap().ambient_lo.value = ambient;
    spec.envelope.as_mut().unwrap().ambient_hi.value = ambient;
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
        temperature: field.get("temperature").unwrap().as_array().unwrap()
            .iter().map(|value| value.as_f64().unwrap()).collect(),
    }
}

fn checked_natural(row: &J, wall: f64, ambient: f64, area: f64, count: usize) -> f64 {
    let natural = row.get("natural").expect("retain each accepted natural endpoint");
    let exchange = natural.get("exchange").unwrap();
    let iterations = n(natural, "iterations");
    assert!(iterations > 1.0 && iterations <= 80.0);
    close(n(exchange, "iterations"), iterations, 0.0);
    assert!(n(natural, "solid_solves") >= iterations);
    close(n(row, "krylov_iterations"), n(natural, "krylov_iterations"), 0.0);
    assert!(
        n(natural, "physical_residual_norm_j")
            <= n(natural, "physical_residual_tolerance_j")
    );
    assert!(n(natural, "physical_energy_residual_j").abs() <= 1e-6);
    let (h, rayleigh, nusselt) = coefficient(wall, ambient, LENGTH_M);
    assert!((0.1..=1e12).contains(&rayleigh));
    let laws = exchange.get("laws").unwrap().as_array().unwrap();
    assert_eq!(laws.len(), count);
    let applied = natural.get("applied_coefficients_w_m2_k").unwrap();
    let physical = natural.get("physical_coefficients_w_m2_k").unwrap();
    let mut sum_w = 0.0;
    for law in laws {
        assert_eq!(law.str_field("card"), Some("convection.churchill-chu-vertical-plate"));
        assert!(matches!(law.get("in_domain"), Some(J::Bool(true))));
        close(n(law, "characteristic_length_m"), LENGTH_M, 0.0);
        close(n(law, "ambient_k"), ambient, 0.0);
        close(n(law, "mean_wall_k"), wall, 2e-8);
        close(n(law, "delta_t_k"), wall - ambient, 2e-8);
        close(n(law, "rayleigh"), rayleigh, 1e-8 * rayleigh.max(1.0));
        close(n(law, "nusselt"), nusselt, 1e-8);
        close(n(law, "htc_w_m2_k"), h, 1e-8);
        let target = law.str_field("target").unwrap();
        close(n(physical, target), h, 1e-8);
        close(
            n(applied, target),
            n(physical, target),
            1e-10 * n(applied, target).abs() + 1e-12,
        );
        close(n(law, "heat_rate_w"), h * area / count as f64 * (wall - ambient), 2e-6);
        sum_w += n(law, "heat_rate_w");
    }
    close(n(natural, "physical_convective_out_w"), h * area * (wall - ambient), 2e-6);
    close(sum_w, n(natural, "physical_convective_out_w"), 2e-6);
    h
}

fn checked_uniform(initial: f64, watts: f64) -> ThermalRun {
    let (spec, bytes) = regular_project(initial, AMBIENT_K, watts);
    let solved = run(&spec, &fixture_cards(), &bytes);
    assert_eq!(solved.temperature.len(), 4, "oracle requires one imported tetra");
    close(n(solved.transient(), "total_steps"), 6.0, 0.0);
    for grid in ["coarse", "fine"] {
        let mut old_expected = initial;
        let mut old_actual = initial;
        let mut previous_time = 0.0;
        for row in solved.rows(grid) {
            let dt = n(row, "dt_s");
            close(n(row, "time_s") - previous_time, dt, 1e-15);
            let expected = bisect(250.0, 370.0, |t| {
                let h = coefficient(t, AMBIENT_K, LENGTH_M).0;
                CAPACITY_J_K * (t - old_expected) / dt
                    + h * area_m2() * (t - AMBIENT_K) - watts
            });
            let actual = n(row, "final_region_max_k");
            close(actual, expected, 2e-8);
            let h = checked_natural(row, actual, AMBIENT_K, area_m2(), 1);
            let natural = row.get("natural").unwrap();
            let stored = CAPACITY_J_K * (actual - old_actual);
            let convective = h * area_m2() * (actual - AMBIENT_K);
            close(n(row, "source_w"), watts, 1e-12);
            close(n(row, "stored_energy_change_j"), stored, 2e-6);
            close(n(natural, "physical_dirichlet_in_w"), 0.0, 1e-12);
            close(n(natural, "physical_radiative_out_w"), 0.0, 1e-12);
            close(n(natural, "radiation_trials"), 0.0, 0.0);
            close(n(natural, "nonlinear_iterations"), 0.0, 0.0);
            close(n(natural, "nonlinear_backtracks"), 0.0, 0.0);
            let residual = stored - dt * (watts - convective);
            close(residual, 0.0, 2e-6);
            close(n(natural, "physical_energy_residual_j"), residual, 2e-6);
            if watts == 0.0 {
                assert_eq!(convective.is_sign_positive(), initial > AMBIENT_K);
                assert!((actual - AMBIENT_K).abs() < (old_actual - AMBIENT_K).abs());
            }
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
    let qoi = &solved.qoi.get("qoi").unwrap().as_array().unwrap()[0];
    close(n(qoi, "value"), final_k, 0.0);
    let scope = solved.qoi.get("thermal_time_scope").unwrap();
    assert_eq!(scope.str_field("kind"), Some("final"));
    close(n(scope, "time_s"), 1.0, 0.0);
    close(
        n(solved.transient().get("energy").unwrap(), "stored_change_j"),
        CAPACITY_J_K * (final_k - initial),
        2e-6,
    );
    close(
        n(solved.receipt.get("energy").unwrap(), "storage_w"),
        n(last, "stored_energy_change_j") / n(last, "dt_s"),
        2e-9,
    );
    let endpoint = solved.receipt.get("natural").unwrap().get("laws").unwrap()
        .as_array().unwrap();
    close(n(&endpoint[0], "mean_wall_k"), final_k, 2e-8);
    solved
}

#[test]
fn g1_native_transient_natural_powered_startup_evaluates_the_actual_endpoint() {
    let solved = checked_uniform(AMBIENT_K, 5.0);
    let initial_guess_h = coefficient(AMBIENT_K + 10.0, AMBIENT_K, LENGTH_M).0;
    for grid in ["coarse", "fine"] {
        for row in solved.rows(grid) {
            let actual = n(row, "final_region_max_k");
            assert!(actual > AMBIENT_K);
            let h = n(
                row.get("natural").unwrap().get("physical_coefficients_w_m2_k").unwrap(),
                "air",
            );
            assert!(
                h < 0.5 * initial_guess_h,
                "startup must update the +10 K seed coefficient to the actual small temperature rise"
            );
        }
    }
}

#[test]
fn g1_native_transient_natural_hot_and_cold_walls_preserve_sign_and_history() {
    for initial in [330.0, 280.0] {
        checked_uniform(initial, 0.0);
    }
}

fn driven_project(bytes: &[u8], card: &str) -> ProjectSpec {
    let mut spec = transient_product::driven_face_project(bytes);
    let mut targets = Vec::new();
    for axis in 0..3 {
        let target = format!("natural-face-{axis}");
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".to_string(),
            name: target.clone(),
            display: "Natural-convection coordinate face".to_string(),
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
            ThermalBoundary { target: target.clone(), condition: natural_law(330.0) },
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
    spec.envelope.as_mut().unwrap().ambient_lo.value = 330.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 330.0;
    spec
}

#[test]
fn g1_native_transient_natural_radiation_and_conductivity_share_the_physical_endpoint() {
    let mut free_temperatures = Vec::new();
    for kind in ["nonlinear", "constant-scalar"] {
        let base = transient_product::transient_conductivity_cards(kind);
        let (cards, card) = radiation_product::emissivity_cards(&base, EMISSIVITY);
        let bytes = tetra_stl();
        let mut spec = driven_project(&bytes, &card);
        rebind_to(&mut spec, &base);
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
                let dt = n(row, "dt_s");
                let expected = bisect(300.0, 330.0, |t| {
                    let wall = (t + 2.0 * 330.0) / 3.0;
                    let h = coefficient(wall, 330.0, LENGTH_M).0;
                    let h_rad = EMISSIVITY * SIGMA * (wall + 330.0)
                        * (wall * wall + 330.0 * 330.0);
                    let k = if kind == "nonlinear" {
                        0.02 * ((t + 3.0 * 330.0) / 4.0) - 4.0
                    } else {
                        2.0
                    };
                    // Three area-1/2 P1 trace rows sum to
                    // (h+h_rad)*(T-330)/4; Cfree=15/24 J/K.
                    0.625 * (t - old_expected) / dt
                        + (0.5 * k + 0.25 * (h + h_rad)) * (t - 330.0)
                });
                let clamp_change = if index == 0 { 90.0 } else { 0.0 };
                let stored = n(row, "stored_energy_change_j");
                let actual = old_actual + stored / 0.625 - clamp_change;
                close(actual, expected, 3e-7);
                let wall = (actual + 2.0 * 330.0) / 3.0;
                let h = checked_natural(row, wall, 330.0, 1.5, 3);
                let natural = row.get("natural").unwrap();
                let convective = 1.5 * h * (wall - 330.0);
                let radiative = EMISSIVITY * SIGMA * 1.5 * (wall.powi(4) - 330.0_f64.powi(4));
                assert!(convective < 0.0 && radiative < 0.0);
                close(n(natural, "physical_radiative_out_w"), radiative, 2e-6);
                close(
                    n(natural, "physical_dirichlet_in_w"),
                    stored / dt + convective + radiative,
                    1e-5,
                );
                assert!(n(natural, "radiation_trials") > n(natural, "iterations"));
                assert!(n(natural, "solid_solves") >= n(natural, "radiation_trials"));
                let radiation = row.get("radiation").unwrap();
                assert!(n(natural, "radiation_trials") >= n(radiation, "iterations"));
                if kind == "nonlinear" {
                    assert!(n(natural, "nonlinear_iterations") > 0.0);
                    let nonlinear = row.get("nonlinear").unwrap();
                    assert!(n(nonlinear, "residual_j") <= n(nonlinear, "threshold_j"));
                } else {
                    close(n(natural, "nonlinear_iterations"), 0.0, 0.0);
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
        "freezing conductivity changes this nonuniform coupled transient"
    );
}

#[test]
fn g0_native_transient_natural_refuses_equilibrium_actual_domain_escape_and_exhausted_work() {
    for failure in ["equilibrium", "actual-rayleigh-domain", "radiation-work"] {
        let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
        let (mut spec, bytes) = if failure == "equilibrium" {
            // Exactly representable temperatures define the intended equilibrium;
            // finite-element row sums or the surface mean can still carry roundoff.
            regular_project(256.0, 256.0, 0.0)
        } else {
            regular_project(350.0, AMBIENT_K, 0.0)
        };
        let expected_code = match failure {
            "equilibrium" => "cli-solve-conduction-natural-unheated",
            "actual-rayleigh-domain" => {
                let length = 7.0;
                let (seed_h, seed_ra, _) = coefficient(AMBIENT_K + 10.0, AMBIENT_K, length);
                assert!((0.1..=1e12).contains(&seed_ra));
                let candidate = (CAPACITY_J_K * 350.0 / 0.5
                    + seed_h * area_m2() * AMBIENT_K)
                    / (CAPACITY_J_K / 0.5 + seed_h * area_m2());
                assert!(coefficient(candidate, AMBIENT_K, length).1 > 1e12);
                if let ThermalBoundaryCondition::NaturalConvection { characteristic_length, .. } =
                    &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .boundaries[0].condition
                {
                    characteristic_length.value = length;
                }
                "cli-solve-conduction-natural-card"
            }
            _ => {
                radiation_product::declare(&mut spec, &card, "air", AMBIENT_K);
                spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .radiation.as_mut().unwrap().max_iterations = 1;
                "cli-solve-conduction-transient"
            }
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
        if failure == "equilibrium" && error.code != expected_code {
            assert!(
                error.code == "cli-solve-conduction-natural-card"
                    || error.code == "cli-solve-conduction-natural-card-domain",
                "{failure}: {error:?}"
            );
            if error.code == "cli-solve-conduction-natural-card" {
                assert!(error.what.contains("outside its validity domain"), "{error:?}");
            }
            let rayleigh = error.what
                .split_once("Ra ")
                .and_then(|(_, suffix)| suffix.split_whitespace().next())
                .map(|value| value.trim_end_matches(')'))
                .and_then(|value| value.parse::<f64>().ok())
                .expect("an equilibrium domain refusal must report its Rayleigh number");
            assert!(
                rayleigh.is_finite() && rayleigh > 0.0 && rayleigh < 0.1,
                "equilibrium roundoff must remain below the card's minimum Ra: {error:?}"
            );
        } else {
            assert_eq!(error.code, expected_code, "{failure}: {error:?}");
        }
        assert_eq!(error.stage, Some("conduction"));
        if failure == "actual-rayleigh-domain" {
            assert!(error.what.contains("Ra"), "{error:?}");
        }
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4,
            "{failure} must retain only the completed predecessor stages"
        );
        assert_eq!(
            artifacts_of_kind(&ledger, "solve-conduction-solution"),
            0,
            "{failure} must not publish a partial endpoint"
        );
    }
}
