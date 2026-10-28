//! G1 native total-enthalpy radiation against an independent uniform P1 mode.
//!
//! The regular tetrahedron makes every nodal mass and boundary trace equal.
//! Its accepted backward-Euler state therefore obeys one scalar energy law;
//! the oracle never calls the conduction, phase-curve or radiation solvers.

use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot, TransientPowerStep,
    TransientRegionEnthalpy, TransientRegionPower, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITY: f64 = 0.85;
const AMBIENT_K: f64 = 293.15;
const INITIAL_K: f64 = 330.0;
const DENSITY: f64 = 10.0;
const CP: f64 = 10.0;
const PHASE_SOURCE: &str =
    "synthetic equilibrium enthalpy chart for the uniform radiation oracle; not measured data";
const PULSE_SOURCE: &str = "synthetic delivered pulse: 20 W for 0.075 s, then zero";

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn values(row: &J, key: &str) -> Vec<f64> {
    row.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
        .collect()
}

fn area_m2() -> f64 {
    8.0 * 3.0_f64.sqrt()
}

fn mass_kg(density: f64) -> f64 {
    density * (8.0 / 3.0)
}

fn convective_w(temperature: f64) -> f64 {
    10.0 * area_m2() * (temperature - AMBIENT_K)
}

fn radiative_w(temperature: f64, reservoir: f64) -> f64 {
    EMISSIVITY * SIGMA * area_m2() * (temperature.powi(4) - reservoir.powi(4))
}

fn uniform_step(old: f64, dt: f64, reservoir: f64, source: f64) -> f64 {
    let residual = |t: f64| {
        mass_kg(DENSITY) * CP * (t - old) / dt + convective_w(t) + radiative_w(t, reservoir)
            - source
    };
    let (mut low, mut high) = (233.15, 398.15);
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

fn project(card: &str, reservoir: f64) -> ProjectSpec {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    transient_product::declaration(&mut spec, DENSITY * CP, 1.0);
    radiation_product::declare(&mut spec, card, "air", reservoir);
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .regions[0]
        .seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    let policy = time(&mut spec);
    policy.initial_temperature =
        fs_project::TransientInitialTemperature::Uniform(QtyAny::new(INITIAL_K, dims::TEMPERATURE));
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            phase_conductivity: None,
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0)
                .unwrap(),
            source: PHASE_SOURCE.into(),
            reference_density: QtyAny::new(DENSITY, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            knots: [233.15, 398.15]
                .into_iter()
                .map(|temperature| TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(
                        CP * (temperature - INITIAL_K),
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
    temperature: Vec<f64>,
    html: String,
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
    let (run, qoi, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    let receipts = stage_receipt_hashes(&ledger, &run);
    let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
    let html = String::from_utf8(artifact_bytes(
        &ledger,
        &receipt_str_field(&report, "report_html"),
    ))
    .unwrap();
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        qoi: J::parse(&qoi).unwrap(),
        temperature: values(&J::parse(&field).unwrap(), "temperature"),
        html,
    }
}

fn check_sources_and_report(solved: &ThermalRun, card: &str, reservoir: f64) {
    let history = solved.transient();
    assert_eq!(
        history.str_field("schema"),
        Some("fs-cli-transient-enthalpy-v1")
    );
    assert_eq!(history.str_field("authority"), Some("Estimated"));
    let materials = history.get("materials").unwrap().as_array().unwrap();
    assert_eq!(materials.len(), 1);
    assert_eq!(materials[0].str_field("source"), Some(PHASE_SOURCE));
    let material_card = fixture_card_identity().0;
    assert_eq!(
        materials[0].str_field("material_card_identity"),
        Some(material_card.as_str())
    );
    let radiation = solved.receipt.get("radiation").unwrap();
    assert_eq!(radiation.str_field("authority"), Some("Estimated"));
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 1);
    let surface = &surfaces[0];
    assert_eq!(surface.str_field("card"), Some(card));
    assert_eq!(surface.str_field("target"), Some("air"));
    assert!(
        fs_blake3::ContentHash::from_hex(surface.str_field("emissivity_receipt").unwrap()).is_some()
    );
    close(n(surface, "emissivity"), EMISSIVITY, 0.0);
    close(n(surface, "query_temperature_k"), 300.0, 0.0);
    close(n(surface, "ambient_temperature_k"), reservoir, 0.0);
    close(n(surface, "area_m2"), area_m2(), 2e-12);
    assert!(solved.html.contains("total-enthalpy"));
    assert!(solved.html.contains("nonlinear_radiative_out_w"));
    assert!(solved.html.contains("Estimated"));
}

#[test]
fn g1_native_enthalpy_prescribed_radiating_bath_counts_physical_reaction_once() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    const BATH_K: f64 = 340.0;
    const RESERVOIR_K: f64 = 270.0;
    let mut spec = project(&card, RESERVOIR_K);
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.clear();
    setup.adiabatic_remainder = true;
    setup.radiation.as_mut().unwrap().surfaces[0].target = "radiating-face".into();
    // Two separate fixed faces cover all four vertices. Radiation acts on a
    // third, explicitly zero-flux face; its incident vertices are constrained
    // by those other faces, without replacing a fixed boundary's physical law.
    for (target, normal, fixed) in [
        ("bath-a", [1.0, -1.0, -1.0], true),
        ("bath-b", [-1.0, 1.0, -1.0], true),
        ("radiating-face", [1.0, 1.0, 1.0], false),
    ] {
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".into(), name: target.into(),
            display: target.into(), expect_id: None,
        });
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".into(), target: target.into(), length_unit: "m".into(),
            selector: MeshSelector::HalfSpace {
                normal, offset: -1.0, side: HalfSpaceSide::AtMost, tolerance: 0.0,
            },
            allow_overlap: true,
        });
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries.push(
            ThermalBoundary {
                target: target.into(),
                condition: if fixed {
                    ThermalBoundaryCondition::FixedTemperature {
                        temperature: QtyAny::new(BATH_K, dims::TEMPERATURE),
                    }
                } else {
                    ThermalBoundaryCondition::HeatFlux {
                        outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
                    }
                },
            },
        );
    }
    let solved = run(&spec, &cards);
    let stored_j = mass_kg(DENSITY) * CP * (BATH_K - INITIAL_K);
    let radiation_w = radiative_w(BATH_K, RESERVOIR_K) / 4.0;
    close(n(solved.receipt.get("solver").unwrap(), "free_dofs"), 0.0, 0.0);
    for &temperature in &solved.temperature {
        close(temperature, BATH_K, 2e-12);
    }
    for grid in ["coarse", "fine"] {
        let mut reservoir_j = 0.0;
        for (index, row) in solved.rows(grid).iter().enumerate() {
            let dt = n(row, "dt_s");
            let step_storage = if index == 0 { stored_j } else { 0.0 };
            let reaction_w = step_storage / dt - 5.0 + radiation_w;
            let radiative = row.get("radiation").unwrap();
            close(n(radiative, "convective_out_w"), 0.0, 2e-9);
            close(n(radiative, "nonlinear_radiative_out_w"), radiation_w, 2e-8);
            close(n(radiative, "physical_dirichlet_in_w"), reaction_w, 2e-8);
            close(n(row, "physical_dirichlet_in_w"), reaction_w, 2e-8);
            close(n(row, "stored_energy_change_j"), step_storage, 2e-8);
            close(n(row, "net_input_w"), step_storage / dt, 2e-8);
            close(n(row, "energy_residual_j"), 0.0, 2e-8);
            assert!(n(radiative, "physical_residual_norm_j")
                <= n(radiative, "physical_residual_tolerance_j"));
            reservoir_j += dt * n(row, "physical_dirichlet_in_w");
        }
        close(reservoir_j, stored_j + radiation_w - 5.0, 2e-8);
    }
    let energy = solved.transient().get("energy").unwrap();
    close(n(energy, "stored_change_j"), stored_j, 2e-8);
    close(n(energy, "integrated_net_input_j"), stored_j, 2e-8);
    close(n(energy, "window_residual_j"), 0.0, 2e-8);
}

#[test]
fn g1_native_enthalpy_radiation_matches_cooling_and_heating_energy_roots() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for reservoir in [270.0, 400.0] {
        let solved = run(&project(&card, reservoir), &cards);
        check_sources_and_report(&solved, &card, reservoir);
        assert_eq!(
            solved.temperature.len(),
            4,
            "oracle requires one regular tetra"
        );
        let history = solved.transient();
        let work = history.get("radiation").unwrap();
        let mut totals = [0.0; 3];
        let mut fine_final = INITIAL_K;
        for (grid, step_count) in [("coarse", 2), ("fine", 4)] {
            let rows = solved.rows(grid);
            assert_eq!(rows.len(), step_count);
            let (mut expected_old, mut actual_old, mut previous_time) = (INITIAL_K, INITIAL_K, 0.0);
            let mut grid_work = [0.0; 3];
            for row in rows {
                let dt = n(row, "dt_s");
                close(n(row, "time_s") - previous_time, dt, 1e-15);
                let expected = uniform_step(expected_old, dt, reservoir, 5.0);
                let actual = n(row, "final_region_max_k");
                close(actual, expected, 2e-8);
                let radiation = row.get("radiation").unwrap();
                let convection = convective_w(actual);
                let nonlinear = radiative_w(actual, reservoir);
                let storage = mass_kg(DENSITY) * CP * (actual - actual_old);
                let physical_residual = storage - dt * (5.0 - convection - nonlinear);
                close(n(row, "source_w"), 5.0, 1e-12);
                close(n(row, "stored_energy_change_j"), storage, 2e-6);
                close(n(row, "net_input_w"), 5.0 - convection - nonlinear, 2e-6);
                close(n(radiation, "convective_out_w"), convection, 2e-6);
                close(n(radiation, "nonlinear_radiative_out_w"), nonlinear, 2e-6);
                close(physical_residual, 0.0, 2e-6);
                close(n(row, "energy_residual_j"), physical_residual, 2e-6);
                close(
                    n(radiation, "physical_energy_residual_j"),
                    physical_residual,
                    2e-6,
                );
                assert!(
                    n(radiation, "physical_energy_residual_j").abs()
                        <= n(history, "energy_tolerance_j")
                );
                assert!(
                    n(radiation, "physical_residual_norm_j")
                        <= n(radiation, "physical_residual_tolerance_j")
                );
                assert!(
                    n(radiation, "iterations") > 1.0,
                    "the fixture must exercise coupling trials"
                );
                for (index, (key, bound)) in [
                    ("iterations", "max_trials_per_step"),
                    ("solid_iterations", "max_newton_updates_per_step"),
                    ("krylov_iterations", "max_krylov_iterations_per_step"),
                ]
                .into_iter()
                .enumerate()
                {
                    assert!(n(radiation, key) <= n(work, bound));
                    grid_work[index] += n(radiation, key);
                }
                assert!(n(radiation, "solid_iterations") >= n(row, "newton_iterations"));
                assert!(n(radiation, "krylov_iterations") >= n(row, "krylov_iterations"));
                expected_old = expected;
                actual_old = actual;
                previous_time = n(row, "time_s");
            }
            close(previous_time, 1.0, 0.0);
            for (index, suffix) in ["trials", "newton_updates", "krylov_iterations"]
                .into_iter()
                .enumerate()
            {
                close(n(work, &format!("{grid}_{suffix}")), grid_work[index], 0.0);
                totals[index] += grid_work[index];
            }
            if grid == "fine" {
                fine_final = actual_old;
            }
        }
        for (index, key) in [
            "total_trials",
            "total_newton_updates",
            "total_krylov_iterations",
        ]
        .into_iter()
        .enumerate()
        {
            close(n(work, key), totals[index], 0.0);
        }
        close(n(history, "total_steps"), 6.0, 0.0);
        assert!(totals[0] > n(history, "total_steps"));
        assert_eq!(fine_final > INITIAL_K, reservoir > INITIAL_K);
        let endpoint = history.get("endpoint").unwrap();
        close(
            values(endpoint, "reference_nodal_mass_kg").iter().sum(),
            mass_kg(DENSITY),
            2e-12,
        );
        for enthalpy in values(endpoint, "specific_enthalpy_j_kg") {
            close(enthalpy, CP * (fine_final - INITIAL_K), 5e-7);
        }
        for fraction in values(endpoint, "liquid_mass_fraction") {
            close(fraction, 0.0, 0.0);
        }
        for &temperature in &solved.temperature {
            close(temperature, fine_final, 2e-8);
        }
        let energy = history.get("energy").unwrap();
        close(
            n(energy, "stored_change_j"),
            mass_kg(DENSITY) * CP * (fine_final - INITIAL_K),
            3e-6,
        );
        close(
            n(energy, "stored_change_j"),
            n(energy, "integrated_net_input_j"),
            3e-6,
        );
        close(
            n(
                &solved.qoi.get("qoi").unwrap().as_array().unwrap()[0],
                "value",
            ),
            fine_final,
            0.0,
        );
        let scope = solved.qoi.get("thermal_time_scope").unwrap();
        assert_eq!(scope.str_field("kind"), Some("final"));
        close(n(scope, "time_s"), 1.0, 0.0);
        let final_radiation = solved.receipt.get("radiation").unwrap();
        assert_eq!(
            n(final_radiation, "nonlinear_radiative_out_w") < 0.0,
            reservoir > fine_final
        );
        let last_radiation = solved
            .rows("fine")
            .last()
            .unwrap()
            .get("radiation")
            .unwrap();
        close(
            n(final_radiation, "solid_solves"),
            n(last_radiation, "iterations"),
            0.0,
        );
        close(
            n(final_radiation, "nonlinear_radiative_out_w"),
            radiative_w(fine_final, reservoir),
            2e-6,
        );
    }
}

#[test]
fn g1_native_radiating_latent_plateau_accounts_for_pulse_and_boundary_energy_once() {
    const HORIZON: f64 = 0.2;
    const SWITCH: f64 = 0.075;
    const INITIAL_H: f64 = 5000.0;
    const LATENT_H: f64 = 10_000.0;
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let mut spec = project(&card, 270.0);
    let policy = time(&mut spec);
    policy.horizon.value = HORIZON;
    policy.max_step.value = 0.1;
    let material = &mut policy.enthalpy.as_mut().unwrap().materials[0];
    material.reference_density.value = 3.0 / 8.0; // The regular tetra contains exactly 1 kg.
    material.initial_specific_enthalpy.value = INITIAL_H;
    material.knots = [
        (-1000.0, 233.15, 0.0),
        (0.0, INITIAL_K, 0.0),
        (LATENT_H, INITIAL_K, 1.0),
        (11_000.0, 398.15, 1.0),
    ]
    .into_iter()
    .map(|(h, temperature, fraction)| TransientEnthalpyKnot {
        specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
        temperature: QtyAny::new(temperature, dims::TEMPERATURE),
        liquid_mass_fraction: fraction,
    })
    .collect();
    policy.power_schedules = vec![TransientRegionPower {
        region: "air".into(),
        source: PULSE_SOURCE.into(),
        steps: [(SWITCH, 20.0), (HORIZON, 0.0)]
            .into_iter()
            .map(|(until, watts)| TransientPowerStep {
                until: QtyAny::new(until, dims::TIME),
                watts: QtyAny::new(watts, dims::POWER),
            })
            .collect(),
    }];
    let solved = run(&spec, &cards);
    check_sources_and_report(&solved, &card, 270.0);
    let boundary_out = convective_w(INITIAL_K) + radiative_w(INITIAL_K, 270.0);
    let expected_change = 20.0 * SWITCH - HORIZON * boundary_out;
    let expected_final_h = INITIAL_H + expected_change;
    assert!(expected_final_h > 0.0 && expected_final_h < LATENT_H);
    for grid in ["coarse", "fine"] {
        let rows = solved.rows(grid);
        assert!(rows.iter().any(|row| n(row, "time_s") == SWITCH));
        let (mut old_h, mut old_time, mut source_j) = (INITIAL_H, 0.0, 0.0);
        for row in rows {
            let dt = n(row, "dt_s");
            let source = if old_time < SWITCH { 20.0 } else { 0.0 };
            let next_h = old_h + dt * (source - boundary_out);
            assert!(next_h > 0.0 && next_h < LATENT_H);
            close(n(row, "source_w"), source, 2e-12);
            close(n(row, "final_region_max_k"), INITIAL_K, 0.0);
            close(n(row, "stored_energy_change_j"), next_h - old_h, 2e-6);
            close(n(row, "mean_liquid_mass_fraction"), next_h / LATENT_H, 2e-9);
            close(n(row, "liquid_mass_kg"), next_h / LATENT_H, 2e-9);
            close(n(row, "net_input_w"), source - boundary_out, 2e-6);
            source_j += dt * source;
            old_h = next_h;
            old_time = n(row, "time_s");
        }
        close(source_j, 20.0 * SWITCH, 2e-12);
        close(old_time, HORIZON, 0.0);
        close(old_h, expected_final_h, 2e-9);
    }
    let history = solved.transient();
    let endpoint = history.get("endpoint").unwrap();
    close(
        values(endpoint, "reference_nodal_mass_kg").iter().sum(),
        1.0,
        2e-12,
    );
    for h in values(endpoint, "specific_enthalpy_j_kg") {
        close(h, expected_final_h, 2e-6);
    }
    for fraction in values(endpoint, "liquid_mass_fraction") {
        close(fraction, expected_final_h / LATENT_H, 2e-9);
    }
    for &temperature in &solved.temperature {
        close(temperature, INITIAL_K, 0.0);
    }
    let energy = history.get("energy").unwrap();
    close(n(energy, "stored_change_j"), expected_change, 2e-6);
    close(n(energy, "integrated_net_input_j"), expected_change, 2e-6);
    assert_eq!(
        history
            .get("temporal_error")
            .unwrap()
            .get("estimated_half_width_k"),
        Some(&J::Null)
    );
    let workload = history.get("workload").unwrap();
    assert_eq!(
        workload.str_field("watts_basis"),
        Some("absolute-delivered-replaces-static-region-power")
    );
    close(
        n(workload, "integrated_scheduled_input_j"),
        20.0 * SWITCH,
        2e-12,
    );
    let schedules = workload.get("schedules").unwrap().as_array().unwrap();
    assert_eq!(schedules.len(), 1);
    assert_eq!(schedules[0].str_field("source"), Some(PULSE_SOURCE));
    assert_eq!(schedules[0].str_field("region"), Some("air"));
}

#[test]
fn g1_pointwise_ambient_uses_one_newton_solve_and_matches_independent_heating_cooling_roots() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for reservoir in [270.0, 400.0] {
        let mut spec = project(&card, reservoir);
        let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .radiation.as_mut().unwrap();
        radiation.spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
        // A single radiation solve suffices because emission is in Newton.
        radiation.max_iterations = 1;
        let solved = run(&spec, &cards);
        assert_eq!(solved.temperature.len(), 4);
        let history = solved.transient();
        let work = history.get("radiation").unwrap();
        assert_eq!(work.str_field("method"), Some("implicit-pointwise-p1-gray-radiation"));
        close(n(work, "total_trials"), n(history, "total_steps"), 0.0);
        close(n(work, "max_trials_per_boundary_iteration"), 1.0, 0.0);
        let final_report = solved.receipt.get("radiation").unwrap();
        assert_eq!(final_report.str_field("model"),
            Some("pointwise-p1-gray-surface-to-black-reservoir"));
        assert_eq!(final_report.str_field("spatial_model"), Some("pointwise"));
        let surface = &final_report.get("surfaces").unwrap().as_array().unwrap()[0];
        assert_eq!(surface.str_field("card"), Some(card.as_str()));
        close(n(surface, "area_m2"), area_m2(), 2e-12);
        for grid in ["coarse", "fine"] {
            let mut old = INITIAL_K;
            for row in solved.rows(grid) {
                let dt = n(row, "dt_s");
                let expected = uniform_step(old, dt, reservoir, 5.0);
                close(n(row, "final_region_max_k"), expected, 2e-8);
                let radiation = row.get("radiation").unwrap();
                close(n(radiation, "solid_solves"), 1.0, 0.0);
                close(n(radiation, "convective_out_w"), convective_w(expected), 3e-6);
                close(n(radiation, "nonlinear_radiative_out_w"),
                    radiative_w(expected, reservoir), 3e-6);
                close(n(row, "stored_energy_change_j"),
                    mass_kg(DENSITY) * CP * (expected - old), 3e-6);
                assert!(n(radiation, "physical_residual_j")
                    <= n(radiation, "physical_residual_threshold_j"));
                assert!(n(radiation, "physical_energy_residual_j").abs()
                    <= n(history, "energy_tolerance_j"));
                old = expected;
            }
            if grid == "fine" {
                for &temperature in &solved.temperature {
                    close(temperature, old, 2e-8);
                }
                close(n(final_report, "nonlinear_radiative_out_w"),
                    radiative_w(old, reservoir), 3e-6);
            }
        }
    }
}

#[test]
fn g0_pointwise_radiation_refuses_unsupported_storage_and_derivatives_before_publishing_a_field() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for variant in ["steady", "capacity", "adjoint"] {
        let mut spec = project(&card, 270.0);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .radiation.as_mut().unwrap().spatial_model =
                Some(fs_project::RadiationSpatialModel::Pointwise);
        match variant {
            "steady" => spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient = None,
            "capacity" => {
                let time = time(&mut spec);
                time.enthalpy = None;
                time.capacities = vec![fs_project::TransientRegionCapacity {
                    region: "air".into(),
                    volumetric_heat_capacity: QtyAny::new(DENSITY * CP, dims::VOLUMETRIC_HEAT_CAPACITY),
                    source: "synthetic unchanged reference capacity".into(),
                }];
            }
            _ => spec.outputs.as_mut().unwrap().push(OutputRequest {
                name: "temperature-final-adjoint".into(), kind: "report".into(), region: None,
            }),
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, transient_radiation::regular_tetra_stl());
        let decoded = decode(&spec);
        assert!(decoded.findings().is_empty(), "{variant}: {:?}", decoded.findings());
        let error = run_solve(&ledger, &CancelGate::new_clock_free(), &mut benign_clock(),
            &decoded, &cards, &mut Vec::new()).unwrap_err();
        assert_eq!(error.code, "cli-solve-conduction-pointwise-radiation", "{variant}: {error:?}");
        if let Some(run) = &error.run {
            assert!(stage_receipt_hashes(&ledger, run).len() <= 4,
                "{variant}: no unsupported pointwise conduction field may be published");
        }
    }
}

#[test]
fn g0_native_enthalpy_radiation_refuses_work_and_source_domain_failures_without_a_field() {
    for failure in [
        "outer-iteration-cap",
        "emissivity-domain",
        "enthalpy-domain",
    ] {
        let (cards, card) = if failure == "emissivity-domain" {
            radiation_product::emissivity_cards_with_domain(
                &fixture_cards(),
                EMISSIVITY,
                295.0,
                315.0,
            )
        } else {
            radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY)
        };
        let mut spec = project(
            &card,
            if failure == "emissivity-domain" {
                300.0
            } else {
                400.0
            },
        );
        if failure == "outer-iteration-cap" {
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
        } else if failure == "enthalpy-domain" {
            time(&mut spec).enthalpy.as_mut().unwrap().materials[0].knots = [-0.1, 0.1]
                .into_iter()
                .map(|delta| TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(CP * delta, dims::SPECIFIC_ENTHALPY),
                    temperature: QtyAny::new(INITIAL_K + delta, dims::TEMPERATURE),
                    liquid_mass_fraction: 0.0,
                })
                .collect();
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
        assert_eq!(
            error.code, "cli-solve-conduction-transient",
            "{failure}: {error:?}"
        );
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4,
            "{failure} must publish no partial conduction field, QoI or report"
        );
    }
}
