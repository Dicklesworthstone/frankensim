//! Native total-enthalpy storage under card-derived natural convection.
//!
//! A regular tetrahedron has a uniform nodal mode, so an independent scalar
//! energy root checks the entire history without using a production solver.
//! Constant-temperature latent plateaus also expose accidental advancement of
//! the physical old enthalpy inside a boundary-coefficient or radiation trial.

use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot, TransientRegionEnthalpy,
    spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const AMBIENT_K: f64 = 293.15;
const LENGTH_M: f64 = 0.1;
const PRESSURE_PA: f64 = 101325.0;
const DENSITY: f64 = 10.0;
const CP: f64 = 10.0;
const EMISSIVITY: f64 = 0.85;
const SIGMA: f64 = 5.670_374_419e-8;
const CARD: &str = "convection.churchill-chu-vertical-plate";
const SOURCE: &str =
    "synthetic equilibrium enthalpy curve for an independent natural-convection energy oracle";

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn values<'a>(row: &'a J, key: &str) -> impl Iterator<Item = f64> + 'a {
    row.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
}

fn area_m2() -> f64 {
    8.0 * 3.0_f64.sqrt()
}

fn mass_kg() -> f64 {
    DENSITY * 8.0 / 3.0
}

/// Churchill-Chu independently evaluated with the stated native air model.
/// The zero-difference limit is only a root bracket, never a claimed card
/// evaluation. Actual accepted endpoints must remain inside the Ra domain.
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
    (0.0263 * nusselt / length, rayleigh, nusselt)
}

fn convection_w(wall: f64, ambient: f64) -> f64 {
    coefficient(wall, ambient, LENGTH_M).0 * area_m2() * (wall - ambient)
}

fn radiation_w(wall: f64, reservoir: f64) -> f64 {
    EMISSIVITY * SIGMA * area_m2() * (wall.powi(4) - reservoir.powi(4))
}

fn uniform_step(old: f64, dt: f64, ambient: f64, source: f64) -> f64 {
    let residual = |temperature: f64| {
        mass_kg() * CP * (temperature - old) / dt
            + convection_w(temperature, ambient)
            - source
    };
    let (mut low, mut high) = (200.0, 450.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (low + high);
        if residual(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    0.5 * (low + high)
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

fn project(initial: f64, ambient: f64, watts: f64) -> ProjectSpec {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    transient_product::declaration(&mut spec, DENSITY * CP, 1.0);
    spec.power.as_mut().unwrap()[0].watts.value = watts;
    let setup = spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    setup.regions[0].seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    setup.boundaries[0].condition = ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(LENGTH_M, dims::LENGTH),
        ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
        correlation: CARD.into(),
    };
    spec.envelope.as_mut().unwrap().ambient_lo.value = ambient;
    spec.envelope.as_mut().unwrap().ambient_hi.value = ambient;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    let policy = time(&mut spec);
    policy.initial_temperature =
        fs_project::TransientInitialTemperature::Uniform(QtyAny::new(initial, dims::TEMPERATURE));
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0)
                .unwrap(),
            source: SOURCE.into(),
            reference_density: QtyAny::new(DENSITY, dims::DENSITY),
            phase_conductivity: None,
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            knots: [200.0, 450.0]
                .into_iter()
                .map(|temperature| TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(
                        CP * (temperature - initial),
                        dims::SPECIFIC_ENTHALPY,
                    ),
                    temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                    liquid_mass_fraction: 0.0,
                })
                .collect(),
        }],
    });
    spec
}

struct ThermalRun {
    receipt: J,
    qoi: J,
    field: J,
}

impl ThermalRun {
    fn transient(&self) -> &J {
        self.receipt.get("transient").unwrap()
    }

    fn rows(&self, grid: &str) -> &[J] {
        self.transient().get(grid).unwrap().as_array().unwrap()
    }
}

fn run(spec: &ProjectSpec, cards: &CardPackSet) -> ThermalRun {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, transient_radiation::regular_tetra_stl());
    let decoded = decode(spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let (_, qoi, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        qoi: J::parse(&qoi).unwrap(),
        field: J::parse(&field).unwrap(),
    }
}

fn check_natural(row: &J, wall: f64, ambient: f64, energy_tolerance: f64) {
    let natural = row
        .get("natural")
        .expect("retain the accepted natural endpoint");
    let exchange = natural.get("exchange").unwrap();
    let iterations = n(natural, "iterations");
    assert!(iterations > 1.0 && iterations <= 80.0);
    close(n(exchange, "iterations"), iterations, 0.0);
    assert!(
        n(natural, "physical_residual_norm_j")
            <= n(natural, "physical_residual_tolerance_j")
    );
    assert!(n(natural, "physical_energy_residual_j").abs() <= energy_tolerance);
    close(
        n(row, "residual_j"),
        n(natural, "physical_residual_norm_j"),
        0.0,
    );
    close(
        n(row, "threshold_j"),
        n(natural, "physical_residual_tolerance_j"),
        0.0,
    );
    let (h, rayleigh, nusselt) = coefficient(wall, ambient, LENGTH_M);
    assert!((0.1..=1e12).contains(&rayleigh));
    let laws = exchange.get("laws").unwrap().as_array().unwrap();
    assert_eq!(laws.len(), 1);
    let law = &laws[0];
    assert_eq!(law.str_field("target"), Some("air"));
    assert_eq!(law.str_field("card"), Some(CARD));
    assert_eq!(law.get("in_domain"), Some(&J::Bool(true)));
    close(n(law, "characteristic_length_m"), LENGTH_M, 0.0);
    close(n(law, "ambient_k"), ambient, 0.0);
    close(n(law, "mean_wall_k"), wall, 2e-8);
    close(n(law, "delta_t_k"), wall - ambient, 2e-8);
    close(n(law, "rayleigh"), rayleigh, 1e-8 * rayleigh.max(1.0));
    close(n(law, "nusselt"), nusselt, 1e-8);
    close(n(law, "htc_w_m2_k"), h, 1e-8);
    close(n(law, "heat_rate_w"), convection_w(wall, ambient), 2e-6);
    close(
        n(natural, "physical_convective_out_w"),
        convection_w(wall, ambient),
        2e-6,
    );
    let applied = n(
        natural.get("applied_coefficients_w_m2_k").unwrap(),
        "air",
    );
    let physical = n(
        natural.get("physical_coefficients_w_m2_k").unwrap(),
        "air",
    );
    close(physical, h, 1e-8);
    close(applied, physical, 1e-10 * applied.abs() + 1e-12);
}

fn check_work(solved: &ThermalRun, with_radiation: bool) {
    let history = solved.transient();
    let work = history.get("natural").unwrap();
    assert_eq!(
        work.str_field("method"),
        Some("implicit-area-mean-natural-convection")
    );
    close(n(work, "max_iterations_per_step"), 80.0, 0.0);
    let counters = [
        ("iterations", "iterations", "max_iterations_per_step"),
        ("solid_solves", "solid_solves", "max_solid_solves_per_step"),
        (
            "nonlinear_iterations",
            "newton_updates",
            "max_newton_updates_per_step",
        ),
        (
            "krylov_iterations",
            "krylov_iterations",
            "max_krylov_iterations_per_step",
        ),
    ];
    let mut totals = [0.0; 4];
    let mut endpoint_evaluations = 0.0;
    let mut maximum_energy_residual = 0.0_f64;
    for grid in ["coarse", "fine"] {
        let mut grid_totals = [0.0; 4];
        for row in solved.rows(grid) {
            let natural = row.get("natural").unwrap();
            for (index, &(key, _, bound)) in counters.iter().enumerate() {
                assert!(n(natural, key) <= n(work, bound));
                grid_totals[index] += n(natural, key);
            }
            assert!(n(natural, "solid_solves") >= n(natural, "iterations"));
            assert!(n(natural, "nonlinear_iterations") >= n(row, "newton_iterations"));
            assert!(n(natural, "krylov_iterations") >= n(row, "krylov_iterations"));
            let audits = n(natural, "endpoint_evaluations");
            assert!(audits >= n(natural, "iterations"));
            assert!(audits <= n(work, "max_endpoint_evaluations_per_step"));
            endpoint_evaluations += audits;
            maximum_energy_residual = maximum_energy_residual
                .max(n(natural, "physical_energy_residual_j").abs());
            if with_radiation {
                assert!(n(natural, "radiation_trials") >= n(natural, "iterations"));
                assert!(
                    n(natural, "radiation_trials")
                        >= n(row.get("radiation").unwrap(), "iterations")
                );
            } else {
                close(n(natural, "radiation_trials"), 0.0, 0.0);
            }
        }
        for (index, &(_, suffix, _)) in counters.iter().enumerate() {
            close(n(work, &format!("{grid}_{suffix}")), grid_totals[index], 0.0);
            totals[index] += grid_totals[index];
        }
    }
    for (index, &(_, suffix, _)) in counters.iter().enumerate() {
        close(n(work, &format!("total_{suffix}")), totals[index], 0.0);
    }
    close(
        n(work, "total_endpoint_evaluations"),
        endpoint_evaluations,
        0.0,
    );
    close(
        n(work, "maximum_physical_energy_residual_j"),
        maximum_energy_residual,
        0.0,
    );
    assert!(
        totals[0] > n(history, "total_steps"),
        "boundary trials are counted independently from advances of the physical history"
    );
}

fn check_endpoint(solved: &ThermalRun, temperature: f64, horizon: f64) {
    let history = solved.transient();
    assert_eq!(
        history.str_field("schema"),
        Some("fs-cli-transient-enthalpy-v1")
    );
    assert_eq!(history.str_field("authority"), Some("Estimated"));
    let materials = history.get("materials").unwrap().as_array().unwrap();
    assert_eq!(materials.len(), 1);
    assert_eq!(materials[0].str_field("source"), Some(SOURCE));
    let identity = fixture_card_identity().0;
    assert_eq!(
        materials[0].str_field("material_card_identity"),
        Some(identity.as_str())
    );
    assert_eq!(values(&solved.field, "temperature").count(), 4);
    for value in values(&solved.field, "temperature") {
        close(value, temperature, 2e-8);
    }
    let qoi = &solved.qoi.get("qoi").unwrap().as_array().unwrap()[0];
    close(n(qoi, "value"), temperature, 0.0);
    let scope = solved.qoi.get("thermal_time_scope").unwrap();
    assert_eq!(scope.str_field("kind"), Some("final"));
    close(n(scope, "time_s"), horizon, 0.0);
    let laws = solved
        .receipt
        .get("natural_convection")
        .unwrap()
        .get("laws")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(laws.len(), 1);
    close(n(&laws[0], "mean_wall_k"), temperature, 2e-8);
    assert_eq!(laws[0].str_field("card"), Some(CARD));
}

#[test]
fn g1_native_enthalpy_natural_cooling_heating_and_ambient_startup_match_energy_roots() {
    for (initial, watts) in [(330.0, 0.0), (270.0, 0.0), (AMBIENT_K, 5.0)] {
        let solved = run(&project(initial, AMBIENT_K, watts), &fixture_cards());
        let history = solved.transient();
        close(n(history, "total_steps"), 6.0, 0.0);
        check_work(&solved, false);
        for (grid, count) in [("coarse", 2), ("fine", 4)] {
            let rows = solved.rows(grid);
            assert_eq!(rows.len(), count);
            let (mut expected_old, mut actual_old, mut elapsed) = (initial, initial, 0.0);
            for row in rows {
                let dt = n(row, "dt_s");
                let expected = uniform_step(expected_old, dt, AMBIENT_K, watts);
                let actual = n(row, "final_region_max_k");
                close(actual, expected, 2e-8);
                close(n(row, "time_s") - elapsed, dt, 1e-15);
                check_natural(row, actual, AMBIENT_K, n(history, "energy_tolerance_j"));
                let natural = row.get("natural").unwrap();
                let boundary_out = convection_w(actual, AMBIENT_K);
                let stored = mass_kg() * CP * (actual - actual_old);
                let residual = stored - dt * (watts - boundary_out);
                close(n(row, "source_w"), watts, 2e-12);
                close(n(row, "stored_energy_change_j"), stored, 2e-6);
                close(n(row, "net_input_w"), watts - boundary_out, 2e-6);
                close(n(row, "energy_residual_j"), residual, 2e-6);
                close(n(natural, "physical_energy_residual_j"), residual, 2e-6);
                close(n(natural, "physical_radiative_out_w"), 0.0, 0.0);
                close(residual, 0.0, 2e-6);
                if watts == 0.0 {
                    assert_eq!(boundary_out > 0.0, initial > AMBIENT_K);
                    assert!(
                        (actual - AMBIENT_K).abs() < (actual_old - AMBIENT_K).abs()
                    );
                } else {
                    assert!(actual > actual_old);
                    let seed_h = coefficient(AMBIENT_K + 10.0, AMBIENT_K, LENGTH_M).0;
                    assert!(
                        coefficient(actual, AMBIENT_K, LENGTH_M).0 < 0.5 * seed_h,
                        "powered startup must replace the initial +10 K guess with the actual wall state"
                    );
                }
                expected_old = expected;
                actual_old = actual;
                elapsed = n(row, "time_s");
            }
            close(elapsed, 1.0, 0.0);
        }
        let final_k = n(
            solved.rows("fine").last().unwrap(),
            "final_region_max_k",
        );
        check_endpoint(&solved, final_k, 1.0);
        let endpoint = history.get("endpoint").unwrap();
        close(
            values(endpoint, "reference_nodal_mass_kg").sum(),
            mass_kg(),
            2e-12,
        );
        for h in values(endpoint, "specific_enthalpy_j_kg") {
            close(h, CP * (final_k - initial), 5e-7);
        }
        for fraction in values(endpoint, "liquid_mass_fraction") {
            close(fraction, 0.0, 0.0);
        }
        let energy = history.get("energy").unwrap();
        close(
            n(energy, "stored_change_j"),
            mass_kg() * CP * (final_k - initial),
            3e-6,
        );
        close(
            n(energy, "stored_change_j"),
            n(energy, "integrated_net_input_j"),
            3e-6,
        );
    }
}

#[test]
fn g1_native_enthalpy_natural_and_radiation_change_latent_fraction_only_once_per_step() {
    const WALL: f64 = 330.0;
    const RESERVOIR: f64 = 270.0;
    const INITIAL_H: f64 = 5000.0;
    const LATENT_H: f64 = 10000.0;
    const HORIZON: f64 = 0.2;
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for with_radiation in [false, true] {
        let mut spec = project(WALL, AMBIENT_K, 5.0);
        let policy = time(&mut spec);
        policy.horizon.value = HORIZON;
        policy.max_step.value = 0.1;
        let material = &mut policy.enthalpy.as_mut().unwrap().materials[0];
        material.reference_density.value = 3.0 / 8.0; // V=8/3 m³ contains exactly 1 kg.
        material.initial_specific_enthalpy.value = INITIAL_H;
        material.knots = [
            (-1000.0, 200.0, 0.0),
            (0.0, WALL, 0.0),
            (LATENT_H, WALL, 1.0),
            (11000.0, 450.0, 1.0),
        ]
        .into_iter()
        .map(|(h, temperature, fraction)| TransientEnthalpyKnot {
            specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
            temperature: QtyAny::new(temperature, dims::TEMPERATURE),
            liquid_mass_fraction: fraction,
        })
        .collect();
        if with_radiation {
            radiation_product::declare(&mut spec, &card, "air", RESERVOIR);
            // Seed the radiation driving state at the exact plateau temperature.
            // The small cap then bounds the declared nested numerical allowance.
            let radiation = spec
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .radiation
                .as_mut()
                .unwrap();
            radiation.surfaces[0].query_temperature.value = WALL;
            radiation.max_iterations = 4;
        }
        let solved = run(&spec, &cards);
        let history = solved.transient();
        check_work(&solved, with_radiation);
        let radiative = if with_radiation {
            radiation_w(WALL, RESERVOIR)
        } else {
            0.0
        };
        let convective = convection_w(WALL, AMBIENT_K);
        let net_input = 5.0 - convective - radiative;
        let final_h = INITIAL_H + HORIZON * net_input;
        assert!(final_h > 0.0 && final_h < LATENT_H);
        for grid in ["coarse", "fine"] {
            let (mut old_h, mut elapsed) = (INITIAL_H, 0.0);
            for row in solved.rows(grid) {
                let dt = n(row, "dt_s");
                let next_h = old_h + dt * net_input;
                assert!(next_h > 0.0 && next_h < LATENT_H);
                close(n(row, "time_s") - elapsed, dt, 1e-15);
                close(n(row, "final_region_max_k"), WALL, 0.0);
                check_natural(row, WALL, AMBIENT_K, n(history, "energy_tolerance_j"));
                let natural = row.get("natural").unwrap();
                close(n(natural, "physical_radiative_out_w"), radiative, 2e-6);
                close(n(row, "source_w"), 5.0, 2e-12);
                close(n(row, "net_input_w"), net_input, 2e-6);
                close(n(row, "stored_energy_change_j"), next_h - old_h, 2e-6);
                close(n(row, "energy_residual_j"), 0.0, 2e-6);
                close(n(row, "liquid_mass_kg"), next_h / LATENT_H, 2e-9);
                close(
                    n(row, "mean_liquid_mass_fraction"),
                    next_h / LATENT_H,
                    2e-9,
                );
                if with_radiation {
                    let radiation = row.get("radiation").unwrap();
                    close(
                        n(radiation, "nonlinear_radiative_out_w"),
                        radiative,
                        2e-6,
                    );
                }
                old_h = next_h;
                elapsed = n(row, "time_s");
            }
            close(old_h, final_h, 2e-9);
            close(elapsed, HORIZON, 0.0);
        }
        check_endpoint(&solved, WALL, HORIZON);
        let endpoint = history.get("endpoint").unwrap();
        close(
            values(endpoint, "reference_nodal_mass_kg").sum(),
            1.0,
            2e-12,
        );
        for h in values(endpoint, "specific_enthalpy_j_kg") {
            close(h, final_h, 2e-6);
        }
        for fraction in values(endpoint, "liquid_mass_fraction") {
            close(fraction, final_h / LATENT_H, 2e-9);
        }
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), HORIZON * net_input, 2e-6);
        close(
            n(energy, "integrated_net_input_j"),
            HORIZON * net_input,
            2e-6,
        );
        assert_eq!(
            history
                .get("temporal_error")
                .unwrap()
                .get("estimated_half_width_k"),
            Some(&J::Null)
        );
        if with_radiation {
            let radiation = solved.receipt.get("radiation").unwrap();
            close(n(radiation, "nonlinear_radiative_out_w"), radiative, 2e-6);
            close(n(radiation, "convective_out_w"), convective, 2e-6);
        }
    }
}

#[test]
fn g0_native_enthalpy_natural_refuses_equilibrium_domains_and_nested_work_without_a_field() {
    for failure in [
        "equilibrium",
        "actual-rayleigh-domain",
        "enthalpy-domain",
        "radiation-work",
    ] {
        let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
        let mut spec = if failure == "equilibrium" {
            project(256.0, 256.0, 0.0)
        } else {
            project(350.0, AMBIENT_K, 0.0)
        };
        if failure == "actual-rayleigh-domain" {
            let length = 7.0;
            let (seed_h, seed_ra, _) = coefficient(AMBIENT_K + 10.0, AMBIENT_K, length);
            assert!((0.1..=1e12).contains(&seed_ra));
            let candidate = (mass_kg() * CP * 350.0 / 0.5
                + seed_h * area_m2() * AMBIENT_K)
                / (mass_kg() * CP / 0.5 + seed_h * area_m2());
            assert!(coefficient(candidate, AMBIENT_K, length).1 > 1e12);
            if let ThermalBoundaryCondition::NaturalConvection {
                characteristic_length,
                ..
            } = &mut spec
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .boundaries[0]
                .condition
            {
                characteristic_length.value = length;
            }
        } else if failure == "enthalpy-domain" {
            time(&mut spec).enthalpy.as_mut().unwrap().materials[0].knots =
                [-0.1, 0.1]
                    .into_iter()
                    .map(|delta| TransientEnthalpyKnot {
                        specific_enthalpy: QtyAny::new(CP * delta, dims::SPECIFIC_ENTHALPY),
                        temperature: QtyAny::new(350.0 + delta, dims::TEMPERATURE),
                        liquid_mass_fraction: 0.0,
                    })
                    .collect();
        } else if failure == "radiation-work" {
            radiation_product::declare(&mut spec, &card, "air", AMBIENT_K);
            spec.cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .radiation
                .as_mut()
                .unwrap()
                .max_iterations = 1;
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, transient_radiation::regular_tetra_stl());
        let decoded = decode(&spec);
        assert!(
            decoded.findings().is_empty(),
            "{failure}: {:?}",
            decoded.findings()
        );
        let error = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decoded,
            &cards,
            &mut Vec::new(),
        )
        .unwrap_err();
        if failure == "equilibrium" && error.code != "cli-solve-conduction-natural-unheated" {
            assert!(
                error.code == "cli-solve-conduction-natural-card"
                    || error.code == "cli-solve-conduction-natural-card-domain",
                "{error:?}"
            );
            let rayleigh = error
                .what
                .split_once("Ra ")
                .and_then(|(_, suffix)| suffix.split_whitespace().next())
                .map(|value| value.trim_end_matches(')'))
                .and_then(|value| value.parse::<f64>().ok())
                .expect("an equilibrium roundoff refusal must retain the actual Rayleigh number");
            assert!(
                rayleigh.is_finite() && rayleigh > 0.0 && rayleigh < 0.1,
                "{error:?}"
            );
        } else if failure == "actual-rayleigh-domain" {
            assert!(
                error.code == "cli-solve-conduction-natural-card"
                    || error.code == "cli-solve-conduction-natural-card-domain",
                "{error:?}"
            );
            assert!(error.what.contains("Ra"), "{error:?}");
        } else if failure != "equilibrium" {
            assert_eq!(
                error.code,
                "cli-solve-conduction-transient",
                "{failure}: {error:?}"
            );
        }
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4,
            "{failure} must not publish a conduction result, QoI or report"
        );
        assert_eq!(
            artifacts_of_kind(&ledger, "solve-conduction-solution"),
            0,
            "{failure} must not publish a partial temperature field"
        );
    }
}
