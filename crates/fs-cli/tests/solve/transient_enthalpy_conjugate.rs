//! Native total-enthalpy storage coupled to a real quasi-steady air branch.
//! The imported regular tetrahedron has one uniform thermal mode. Independent
//! scalar energy roots and exponential air heating therefore check every step
//! without calling the production enthalpy, radiation or air-march operators.
//! Mixed-boundary plateaus additionally check signed natural heat, radiation
//! and fan transport while allowing a nonuniform spatial phase fraction.

use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot,
    TransientInitialTemperature, TransientPowerStep, TransientRegionEnthalpy,
    TransientRegionPower, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const INITIAL_K: f64 = 330.0;
const INLET_K: f64 = 293.15;
const DENSITY: f64 = 10.0;
const CP: f64 = 10.0;
const AIR_CP: f64 = 1007.0;
const POWER_W: f64 = 5.0;
const EMISSIVITY: f64 = 0.85;
const SIGMA: f64 = 5.670_374_419e-8;
const LATENT_H: f64 = 200.0;
const INITIAL_LATENT_H: f64 = 150.0;

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
    row.get(key).unwrap().as_array().unwrap().iter()
        .map(|value| value.as_f64().unwrap()).collect()
}

fn area_m2() -> f64 {
    8.0 * 3.0_f64.sqrt()
}

fn mass_kg() -> f64 {
    DENSITY * 8.0 / 3.0
}

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap()
}

fn capacity_project(horizon: f64) -> (ProjectSpec, Vec<u8>) {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    transient_product::declaration(&mut spec, DENSITY * CP, horizon);
    time(&mut spec).initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(INITIAL_K, dims::TEMPERATURE));
    time(&mut spec).max_steps = 12;
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().regions[0].seed =
        [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    (spec, bytes)
}

fn enthalpy(spec: &mut ProjectSpec, initial_h: f64, knots: &[(f64, f64, f64)]) {
    let policy = time(spec);
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 16,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0)
                .unwrap(),
            source: "synthetic enthalpy and density for an independent air-cooling energy oracle; not measured AA6061 data".into(),
            reference_density: QtyAny::new(DENSITY, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(initial_h, dims::SPECIFIC_ENTHALPY),
            knots: knots.iter().map(|&(h, temperature, fraction)| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: fraction,
            }).collect(),
            phase_conductivity: None,
        }],
    });
}

fn sensible_enthalpy(spec: &mut ProjectSpec) {
    enthalpy(spec, 0.0, &[
        (CP * (233.15 - INITIAL_K), 233.15, 0.0),
        (CP * (398.15 - INITIAL_K), 398.15, 0.0),
    ]);
}

struct ThermalRun {
    receipt: J,
    qoi: J,
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
    let decoded = decode(spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let (_, qoi, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        qoi: J::parse(&qoi).unwrap(),
        temperature: values(&J::parse(&field).unwrap(), "temperature"),
    }
}

/// Re-derive G = m_dot Cp (1 - exp(-h A / (m_dot Cp))) from the independently
/// lowered fan/vent flow and the admitted correlation's coefficient.
fn conductance(exchange: &J) -> f64 {
    conductance_for_area(exchange, area_m2(), "air")
}

fn conductance_for_area(exchange: &J, area: f64, branch: &str) -> f64 {
    assert_eq!(exchange.str_field("branch"), Some(branch));
    let path = format!("vent:{branch}");
    assert_eq!(exchange.str_field("path"), Some(path.as_str()));
    close(n(exchange, "inlet_k"), INLET_K, 0.0);
    close(n(exchange.get("air_properties").unwrap(), "specific_heat_j_kg_k"), AIR_CP, 0.0);
    close(
        n(exchange, "mass_flow_kg_s"),
        n(exchange, "air_density_kg_m3") * n(exchange.get("flow_m3_s").unwrap(), "mid"),
        2e-14,
    );
    let capacity_rate = n(exchange, "mass_flow_kg_s") * AIR_CP;
    assert!(capacity_rate > 0.0);
    let segments = exchange.get("segments").unwrap().as_array().unwrap();
    assert_eq!(segments.len(), 1);
    let segment = &segments[0];
    assert_eq!(segment.str_field("card"), Some("convection.gnielinski"));
    assert_eq!(segment.get("in_domain"), Some(&J::Bool(true)));
    assert!(n(segment, "reynolds") > 3000.0);
    close(n(segment, "wetted_area_m2"), area, 2e-12);
    let h = n(segment, "htc_w_m2_k");
    close(h, n(segment, "nusselt") * 0.0263 / 0.02, 1e-12);
    assert!(h > 0.0);
    let ntu = h * area / capacity_rate;
    close(n(segment, "ntu"), ntu, 2e-12);
    close(n(segment, "effectiveness"), -(-ntu).exp_m1(), 2e-12);
    capacity_rate * -(-ntu).exp_m1()
}

fn check_exchange(row: &J, temperature: f64, radiative_w: f64, enthalpy_mode: bool) -> f64 {
    let coupled = row.get("conjugate").unwrap();
    let exchange = coupled.get("exchange").unwrap();
    let air_w = conductance(exchange) * (temperature - INLET_K);
    assert!(air_w > 0.0, "the hot solid must warm the admitted air stream");
    close(n(exchange, "air_total_w"), air_w, 3e-6);
    close(n(exchange, "solid_total_w"), air_w, 3e-6);
    close(n(coupled, "air_heat_gain_w"), air_w, 3e-6);
    close(n(coupled, "off_path_convective_out_w"), 0.0, 1e-12);
    close(n(coupled, "radiative_out_w"), radiative_w, 3e-6);
    assert!(n(coupled, "coupled_energy_residual_j").abs() <= 1e-6);
    if enthalpy_mode {
        assert!(n(coupled, "physical_residual_norm_j")
            <= n(coupled, "physical_residual_tolerance_j"));
    }
    let tolerance = n(exchange, "balance_tolerance_w");
    assert!(n(exchange, "max_region_imbalance_w").abs() <= tolerance);
    assert!(n(exchange, "decomposition_residual_w").abs() <= tolerance);
    let outlet = INLET_K + air_w / (n(exchange, "mass_flow_kg_s") * AIR_CP);
    close(n(exchange, "outlet_k"), outlet, 3e-8);
    let segment = &exchange.get("segments").unwrap().as_array().unwrap()[0];
    close(n(segment, "air_in_k"), INLET_K, 0.0);
    close(n(segment, "air_out_k"), outlet, 3e-8);
    close(n(segment, "air_heat_rate_w"), air_w, 3e-6);
    air_w
}

fn check_final(solved: &ThermalRun, expected_k: f64, horizon: f64) {
    assert_eq!(solved.temperature.len(), 4, "oracle requires one imported regular tetra");
    for &temperature in &solved.temperature {
        close(temperature, expected_k, 5e-8);
    }
    let qoi = &solved.qoi.get("qoi").unwrap().as_array().unwrap()[0];
    close(n(qoi, "value"), solved.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max), 0.0);
    let scope = solved.qoi.get("thermal_time_scope").unwrap();
    assert_eq!(scope.str_field("kind"), Some("final"));
    close(n(scope, "time_s"), horizon, 0.0);
    let last = solved.rows("fine").last().unwrap().get("conjugate").unwrap();
    close(
        n(solved.receipt.get("conjugate").unwrap(), "air_total_w"),
        n(last, "air_heat_gain_w"),
        1e-12,
    );
}

#[test]
fn g1_native_enthalpy_airflow_recovers_constant_capacity_and_exponential_air_heating() {
    let (capacity_spec, bytes) = capacity_project(1.0);
    let mut enthalpy_spec = capacity_spec.clone();
    sensible_enthalpy(&mut enthalpy_spec);
    let cards = fixture_cards();
    let capacity = run(&capacity_spec, &cards, &bytes);
    let solved = run(&enthalpy_spec, &cards, &bytes);
    assert_eq!(solved.transient().str_field("schema"), Some("fs-cli-transient-enthalpy-v1"));
    let mut final_expected = INITIAL_K;
    for grid in ["coarse", "fine"] {
        let (a, b) = (capacity.rows(grid), solved.rows(grid));
        assert_eq!(a.len(), b.len());
        let mut previous = [INITIAL_K; 2];
        let mut expected_old = INITIAL_K;
        for (capacity_row, enthalpy_row) in a.iter().zip(b) {
            let dt = n(enthalpy_row, "dt_s");
            close(n(capacity_row, "dt_s"), dt, 0.0);
            close(n(capacity_row, "time_s"), n(enthalpy_row, "time_s"), 0.0);
            let g = conductance(enthalpy_row.get("conjugate").unwrap().get("exchange").unwrap());
            let expected = (mass_kg() * CP * expected_old / dt + POWER_W + g * INLET_K)
                / (mass_kg() * CP / dt + g);
            for (index, row) in [capacity_row, enthalpy_row].into_iter().enumerate() {
                let actual = n(row, "final_region_max_k");
                close(actual, expected, 5e-8);
                let air_w = check_exchange(row, actual, 0.0, index == 1);
                let stored = mass_kg() * CP * (actual - previous[index]);
                close(n(row, "stored_energy_change_j"), stored, 3e-6);
                close(n(row, "source_w"), POWER_W, 1e-12);
                close(stored - dt * (POWER_W - air_w), 0.0, 3e-6);
                previous[index] = actual;
            }
            close(n(capacity_row, "final_region_max_k"),
                n(enthalpy_row, "final_region_max_k"), 5e-8);
            expected_old = expected;
        }
        if grid == "fine" {
            final_expected = expected_old;
        }
    }
    check_final(&capacity, final_expected, 1.0);
    check_final(&solved, final_expected, 1.0);
    let endpoint = solved.transient().get("endpoint").unwrap();
    let h = values(endpoint, "specific_enthalpy_j_kg");
    let masses = values(endpoint, "reference_nodal_mass_kg");
    assert_eq!(h.len(), 4);
    assert_eq!(masses.len(), 4);
    close(masses.iter().sum(), mass_kg(), 2e-12);
    for (index, &temperature) in solved.temperature.iter().enumerate() {
        close(temperature, capacity.temperature[index], 5e-8);
        close(h[index], CP * (temperature - INITIAL_K), 2e-9);
        assert!(masses[index] > 0.0);
    }
    for fraction in values(endpoint, "liquid_mass_fraction") {
        close(fraction, 0.0, 0.0);
    }
    let stored: f64 = masses.iter().zip(&h).map(|(m, h)| m * h).sum();
    assert!(stored < 0.0, "air cooling removes initial sensible energy");
    for history in [capacity.transient(), solved.transient()] {
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), stored, 1e-5);
        close(n(energy, "integrated_net_input_j"), stored, 1e-5);
    }
}

fn temperature_at(h: f64) -> f64 {
    if h < 0.0 {
        INITIAL_K + h / CP
    } else if h <= LATENT_H {
        INITIAL_K
    } else {
        INITIAL_K + (h - LATENT_H) / CP
    }
}

fn radiation_w(temperature: f64, enabled: bool) -> f64 {
    if enabled {
        EMISSIVITY * SIGMA * area_m2() * (temperature.powi(4) - INLET_K.powi(4))
    } else {
        0.0
    }
}

fn phase_step(old_h: f64, dt: f64, g: f64, radiation: bool) -> f64 {
    let residual = |h: f64| {
        let t = temperature_at(h);
        mass_kg() * (h - old_h) / dt + g * (t - INLET_K)
            + radiation_w(t, radiation) - POWER_W
    };
    let mut low = CP * (233.15 - INITIAL_K);
    let mut high = LATENT_H + CP * (398.15 - INITIAL_K);
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

#[test]
fn g1_native_enthalpy_airflow_removes_latent_heat_and_counts_radiation_separately() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for with_radiation in [false, true] {
        let (mut spec, bytes) = capacity_project(2.0);
        enthalpy(&mut spec, INITIAL_LATENT_H, &[
            (CP * (233.15 - INITIAL_K), 233.15, 0.0),
            (0.0, INITIAL_K, 0.0),
            (LATENT_H, INITIAL_K, 1.0),
            (LATENT_H + CP * (398.15 - INITIAL_K), 398.15, 1.0),
        ]);
        if with_radiation {
            radiation_product::declare(&mut spec, &card, "air", INLET_K);
            let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
            let radiation = setup.radiation.as_mut().unwrap();
            radiation.surfaces[0].query_temperature.value = INITIAL_K;
            radiation.relaxation = 1.0;
            radiation.max_iterations = 24;
            // One uniform piecewise-linear storage mode needs few Newton
            // updates; this also bounds the nested air/radiation allowance.
            setup.transient.as_mut().unwrap().enthalpy.as_mut().unwrap()
                .max_newton_iterations = 8;
        }
        let solved = run(&spec, &cards, &bytes);
        let mut final_h = INITIAL_LATENT_H;
        for grid in ["coarse", "fine"] {
            let mut old_h = INITIAL_LATENT_H;
            let mut elapsed = 0.0;
            let mut air_j = 0.0;
            let mut radiation_j = 0.0;
            let mut saw_plateau = false;
            for row in solved.rows(grid) {
                let dt = n(row, "dt_s");
                close(n(row, "time_s") - elapsed, dt, 1e-15);
                let g = conductance(row.get("conjugate").unwrap().get("exchange").unwrap());
                let expected_h = phase_step(old_h, dt, g, with_radiation);
                let expected_k = temperature_at(expected_h);
                let expected_fraction = (expected_h / LATENT_H).clamp(0.0, 1.0);
                saw_plateau |= expected_fraction > 0.0 && expected_fraction < 1.0;
                close(n(row, "final_region_max_k"), expected_k, 5e-8);
                close(n(row, "mean_liquid_mass_fraction"), expected_fraction, 2e-8);
                close(n(row, "liquid_mass_kg"), mass_kg() * expected_fraction, 1e-6);
                let radiative = radiation_w(expected_k, with_radiation);
                let air = check_exchange(row, expected_k, radiative, true);
                close(n(row, "source_w"), POWER_W, 1e-12);
                close(n(row, "stored_energy_change_j"), mass_kg() * (expected_h - old_h), 3e-6);
                close(n(row, "net_input_w"), POWER_W - air - radiative, 3e-6);
                assert!(n(row, "energy_residual_j").abs() <= 1e-6);
                if with_radiation {
                    close(n(row.get("radiation").unwrap(), "nonlinear_radiative_out_w"),
                        radiative, 3e-6);
                }
                air_j += dt * air;
                radiation_j += dt * radiative;
                old_h = expected_h;
                elapsed = n(row, "time_s");
            }
            assert!(saw_plateau, "the trajectory must retain a genuine latent state");
            close(elapsed, 2.0, 0.0);
            let stored = mass_kg() * (old_h - INITIAL_LATENT_H);
            close(stored + air_j + radiation_j, POWER_W * elapsed, 1e-5);
            assert!(air_j > 0.0);
            if with_radiation {
                assert!(radiation_j > 0.0);
                assert!(old_h < 0.0, "combined cooling must leave the plateau and cool the solid");
            }
            if grid == "fine" {
                final_h = old_h;
            }
        }
        check_final(&solved, temperature_at(final_h), 2.0);
        let endpoint = solved.transient().get("endpoint").unwrap();
        let masses = values(endpoint, "reference_nodal_mass_kg");
        let enthalpies = values(endpoint, "specific_enthalpy_j_kg");
        close(masses.iter().sum(), mass_kg(), 2e-12);
        assert_eq!(masses.len(), 4);
        assert_eq!(enthalpies.len(), 4);
        for h in &enthalpies {
            close(*h, final_h, 1e-6);
        }
        for fraction in values(endpoint, "liquid_mass_fraction") {
            close(fraction, (final_h / LATENT_H).clamp(0.0, 1.0), 2e-8);
        }
        let stored: f64 = masses.iter().zip(&enthalpies)
            .map(|(mass, h)| mass * (h - INITIAL_LATENT_H)).sum();
        let energy = solved.transient().get("energy").unwrap();
        close(n(energy, "stored_change_j"), stored, 3e-6);
        close(stored, mass_kg() * (final_h - INITIAL_LATENT_H), 1e-5);
        close(n(energy, "integrated_net_input_j"), stored, 1e-5);
    }
}

#[test]
fn g0_native_enthalpy_airflow_refuses_source_domains_and_exhausted_nested_work() {
    for failure in ["air-card-domain", "enthalpy-domain", "radiation-work", "adjoint"] {
        let (mut spec, bytes) = capacity_project(1.0);
        sensible_enthalpy(&mut spec);
        let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
        let expected_code = if failure == "air-card-domain" {
            if let ThermalBoundaryCondition::AirflowConvection { channel_length, .. } =
                &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .boundaries[0].condition
            {
                channel_length.value = 0.1; // L/Dh=5 violates the admitted correlation.
            }
            "cli-solve-conduction-airflow-card-domain"
        } else {
            if failure == "enthalpy-domain" {
                enthalpy(&mut spec, 0.0, &[
                    (-0.1, INITIAL_K - 0.01, 0.0),
                    (0.1, INITIAL_K + 0.01, 0.0),
                ]);
            } else if failure == "radiation-work" {
                radiation_product::declare(&mut spec, &card, "air", INLET_K);
                spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .radiation.as_mut().unwrap().max_iterations = 1;
            } else {
                spec.outputs.as_mut().unwrap().push(OutputRequest {
                    name: "temperature-final-adjoint".into(),
                    kind: "report".into(),
                    region: None,
                });
            }
            "cli-solve-conduction-transient"
        };
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let decoded = decode(&spec);
        assert!(decoded.findings().is_empty(), "{failure}: {:?}", decoded.findings());
        let error = run_solve(
            &ledger, &CancelGate::new_clock_free(), &mut benign_clock(),
            &decoded, &cards, &mut Vec::new(),
        ).unwrap_err();
        assert_eq!(error.code, expected_code, "{failure}: {error:?}");
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(), 4);
        assert_eq!(artifacts_of_kind(&ledger, "solve-conduction-solution"), 0,
            "{failure} must not publish an unaccepted enthalpy field");
    }
}

const MIXED_INITIAL_H: f64 = 5000.0;
const MIXED_LATENT_H: f64 = 10_000.0;
const MIXED_HORIZON: f64 = 0.2;
const MIXED_SWITCH: f64 = 0.075;
const NATURAL_K: f64 = 350.0;
const RADIATION_K: f64 = 270.0;

fn mixed_plateau_project(card: &str) -> (ProjectSpec, Vec<u8>) {
    let (mut spec, bytes) = capacity_project(MIXED_HORIZON);
    enthalpy(&mut spec, MIXED_INITIAL_H, &[
        (-1000.0, 233.15, 0.0),
        (0.0, INITIAL_K, 0.0),
        (MIXED_LATENT_H, INITIAL_K, 1.0),
        (11_000.0, 398.15, 1.0),
    ]);
    let policy = time(&mut spec);
    policy.max_step.value = 0.1;
    policy.max_steps = 16;
    let model = policy.enthalpy.as_mut().unwrap();
    model.max_newton_iterations = 8;
    model.materials[0].reference_density.value = 3.0 / 8.0; // V=8/3: 1 kg.
    policy.power_schedules = vec![TransientRegionPower {
        region: "air".into(),
        source: "declared 20 W pulse ending at 0.075 s; replaces static power".into(),
        steps: [(MIXED_SWITCH, 20.0), (MIXED_HORIZON, 0.0)]
            .into_iter().map(|(until, watts)| TransientPowerStep {
                until: QtyAny::new(until, dims::TIME),
                watts: QtyAny::new(watts, dims::POWER),
            }).collect(),
    }];
    spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
        parent: "enclosure".into(),
        name: "natural-wall".into(),
        display: "Separate warm natural-convection reservoir".into(),
        expect_id: None,
    });
    let assignments = spec.assignments.as_mut().unwrap();
    assignments[0].allow_overlap = true;
    assignments.push(GeometryAssignment {
        artifact: "enclosure".into(),
        target: "natural-wall".into(),
        length_unit: "m".into(),
        // The face opposite (1,1,1) lies on x+y+z=-1 and has area 2*sqrt(3).
        selector: MeshSelector::HalfSpace {
            normal: [1.0, 1.0, 1.0],
            offset: -1.0,
            side: HalfSpaceSide::AtMost,
            tolerance: 0.0,
        },
        allow_overlap: true,
    });
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .boundaries.push(ThermalBoundary {
            target: "natural-wall".into(),
            condition: ThermalBoundaryCondition::NaturalConvection {
                characteristic_length: QtyAny::new(0.1, dims::LENGTH),
                ambient_temperature: QtyAny::new(NATURAL_K, dims::TEMPERATURE),
                correlation: "convection.churchill-chu-vertical-plate".into(),
            },
        });
    radiation_product::declare(&mut spec, card, "natural-wall", RADIATION_K);
    let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap();
    radiation.relaxation = 1.0;
    radiation.max_iterations = 24;
    (spec, bytes)
}

/// Direct Churchill-Chu law with the stated native air-property constants.
fn mixed_natural_heat(area: f64) -> f64 {
    let film = 0.5 * (INITIAL_K + NATURAL_K);
    let density = 101325.0 / (287.05 * film);
    let nu = 1.846e-5 / density;
    let ra = 9.80665 * (INITIAL_K - NATURAL_K).abs() * 0.1_f64.powi(3)
        / (film * nu * (nu / 0.707));
    let nusselt = (0.825 + 0.387 * ra.powf(1.0 / 6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0)).powi(2);
    (0.0263 * nusselt / 0.1) * area * (INITIAL_K - NATURAL_K)
}

fn check_mixed_plateau(solved: &ThermalRun, air_area: f64, natural_area: f64, branch: &str) {
    let natural = mixed_natural_heat(natural_area);
    let radiation = EMISSIVITY * SIGMA * natural_area
        * (INITIAL_K.powi(4) - RADIATION_K.powi(4));
    assert!(natural < 0.0 && radiation > 0.0);
    let mut final_change = 0.0;
    for grid in ["coarse", "fine"] {
        let (mut old_time, mut cumulative, mut source_j) = (0.0, 0.0, 0.0);
        assert!(solved.rows(grid).iter().any(|row| n(row, "time_s") == MIXED_SWITCH));
        for row in solved.rows(grid) {
            let coupled = row.get("conjugate").unwrap();
            let exchange = coupled.get("exchange").unwrap();
            let g = conductance_for_area(exchange, air_area, branch);
            let air = g * (INITIAL_K - INLET_K);
            close(n(exchange, "air_total_w"), air, 3e-6);
            close(
                AIR_CP * n(exchange, "mass_flow_kg_s")
                    * (n(exchange, "outlet_k") - n(exchange, "inlet_k")),
                air, 3e-6,
            );
            close(n(coupled, "air_heat_gain_w"), air, 3e-6);
            let tolerance = n(exchange, "balance_tolerance_w");
            assert!(n(exchange, "max_region_imbalance_w").abs() <= tolerance);
            assert!(n(exchange, "decomposition_residual_w").abs() <= tolerance);
            assert!(n(coupled, "physical_residual_norm_j")
                <= n(coupled, "physical_residual_tolerance_j"));
            assert!(n(coupled, "coupled_energy_residual_j").abs() <= 1e-6);
            close(n(coupled, "physical_dirichlet_in_w"), 0.0, 0.0);
            let dt = n(row, "dt_s");
            close(n(row, "time_s") - old_time, dt, 1e-15);
            let source = if old_time < MIXED_SWITCH { 20.0 } else { 0.0 };
            let stored = dt * (source - air - natural - radiation);
            cumulative += stored;
            source_j += dt * source;
            let mean_h = MIXED_INITIAL_H + cumulative; // Each fixture contains 1 kg.
            assert!(mean_h > 0.0 && mean_h < MIXED_LATENT_H);
            close(n(row, "source_w"), source, 2e-12);
            close(n(row, "final_region_max_k"), INITIAL_K, 2e-10);
            close(n(row, "stored_energy_change_j"), stored, 3e-6);
            close(n(row, "net_input_w"), source - air - natural - radiation, 3e-6);
            close(n(row, "mean_liquid_mass_fraction"), mean_h / MIXED_LATENT_H, 2e-9);
            close(n(row, "liquid_mass_kg"), mean_h / MIXED_LATENT_H, 2e-9);
            close(n(coupled, "off_path_convective_out_w"), natural, 3e-6);
            close(n(coupled, "radiative_out_w"), radiation, 3e-6);
            let laws = row.get("natural").unwrap().get("exchange").unwrap()
                .get("laws").unwrap().as_array().unwrap();
            assert_eq!(laws.len(), 1);
            assert_eq!(laws[0].str_field("target"), Some("natural-wall"));
            close(n(&laws[0], "mean_wall_k"), INITIAL_K, 2e-10);
            close(n(&laws[0], "heat_rate_w"), natural, 3e-6);
            close(n(row.get("radiation").unwrap(), "nonlinear_radiative_out_w"),
                radiation, 3e-6);
            old_time = n(row, "time_s");
        }
        close(source_j, 20.0 * MIXED_SWITCH, 2e-12);
        close(old_time, MIXED_HORIZON, 0.0);
        if grid == "coarse" {
            final_change = cumulative;
        } else {
            close(cumulative, final_change, 3e-6);
        }
    }
    assert!(final_change < 0.0, "the fan-cooled material must partially solidify");
    let history = solved.transient();
    let endpoint = history.get("endpoint").unwrap();
    let masses = values(endpoint, "reference_nodal_mass_kg");
    let enthalpies = values(endpoint, "specific_enthalpy_j_kg");
    let fractions = values(endpoint, "liquid_mass_fraction");
    assert_eq!(masses.len(), 4);
    assert_eq!(enthalpies.len(), masses.len());
    assert_eq!(fractions.len(), masses.len());
    close(masses.iter().sum(), 1.0, 2e-12);
    for (&h, &fraction) in enthalpies.iter().zip(&fractions) {
        assert!(h > 0.0 && h < MIXED_LATENT_H);
        close(fraction, h / MIXED_LATENT_H, 2e-12);
    }
    close(
        masses.iter().zip(&enthalpies)
            .map(|(mass, h)| mass * (h - MIXED_INITIAL_H)).sum(),
        final_change, 3e-6,
    );
    let energy = history.get("energy").unwrap();
    close(n(energy, "stored_change_j"), final_change, 3e-6);
    close(n(energy, "integrated_net_input_j"), final_change, 3e-6);
    close(n(history.get("workload").unwrap(), "integrated_scheduled_input_j"),
        20.0 * MIXED_SWITCH, 2e-12);
    assert_eq!(history.get("temporal_error").unwrap()
        .get("estimated_half_width_k"), Some(&J::Null));
    check_final(solved, INITIAL_K, MIXED_HORIZON);
}

#[test]
fn g1_pointwise_radiation_composes_with_natural_and_forced_air_without_double_counting_latent_heat() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let (mut spec, bytes) = mixed_plateau_project(&card);
    let declaration = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap();
    declaration.spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
    declaration.max_iterations = 1;
    let solved = run(&spec, &cards, &bytes);
    // This existing oracle independently closes Churchill-Chu reservoir heat,
    // exponential branch air heating, radiation, and the latent mass account.
    check_mixed_plateau(&solved, 0.75 * area_m2(), 0.25 * area_m2(), "air");
    assert_eq!(solved.receipt.get("radiation").unwrap().str_field("spatial_model"), Some("pointwise"));
    for grid in ["coarse", "fine"] {
        for row in solved.rows(grid) {
            let radiation = row.get("radiation").unwrap();
            assert_eq!(radiation.str_field("model"), Some("pointwise-p1"));
            close(n(radiation, "solid_solves"), 1.0, 0.0);
            assert!(n(row.get("conjugate").unwrap(), "physical_residual_norm_j")
                <= n(row.get("conjugate").unwrap(), "physical_residual_tolerance_j"));
        }
    }
}

#[test]
fn g1_native_mixed_enthalpy_keeps_signed_natural_heat_and_radiation_out_of_air() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let (spec, bytes) = mixed_plateau_project(&card);
    let solved = run(&spec, &cards, &bytes);
    check_mixed_plateau(&solved, 0.75 * area_m2(), 0.25 * area_m2(), "air");
}

#[test]
fn g1_native_mixed_enthalpy_example_runs_with_committed_material_and_surface_cards() {
    let decoded = fs_project::parse_sexpr(include_str!(
        "../../../../data/reference-project/cooling-enthalpy-airflow.fsim"
    )).unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let cards = CardPackSet::admit(vec![
        raw_pack(
            CardPackKind::Material, "aa6061.fsmcdpk",
            include_bytes!("../../../../data/reference-project/aa6061.fsmcdpk").to_vec(),
        ),
        raw_pack(
            CardPackKind::Material, "gray-surface.fsmcdpk",
            include_bytes!("../../../../data/reference-project/gray-surface.fsmcdpk").to_vec(),
        ),
    ]).unwrap();
    let solved = run(
        &decoded.spec, &cards,
        include_bytes!("../../../../data/reference-project/plate.stl"),
    );
    assert_eq!(solved.transient().str_field("schema"), Some("fs-cli-transient-enthalpy-v1"));
    // Unit right tetra: z=0 natural face has area 1/2; the other three faces
    // provide forced-air area 1+sqrt(3)/2. Density 6 kg/m3 makes its mass 1 kg.
    check_mixed_plateau(&solved, 1.0 + 0.5 * 3.0_f64.sqrt(), 0.5, "solid");
    let radiation = solved.receipt.get("radiation").unwrap();
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 1);
    assert_eq!(surfaces[0].str_field("card"),
        Some("63485429663ba24d53d67a7d3b03ab0611f17f1e3e445e6a1ef4f636093a6e4f"));
    close(n(&surfaces[0], "emissivity"), EMISSIVITY, 0.0);
}
