//! Regional capacity, live contacts and area-mean radiation share one history.

use super::*;
use fs_project::{TransientPowerStep, TransientRegionPower};

const START_K: f64 = 300.0;
const BATH_K: f64 = 290.0;
const RESERVOIR_K: f64 = 297.0;
const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITY: f64 = 0.85; // Existing gray-surface.fsmcdpk fixture declaration.

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Cooling { Prescribed, Natural, Air, Mixed }

fn n(row: &J, key: &str) -> f64 { row.f64_field(key).unwrap() }

pub(super) fn radiating_fixture(cooling: Cooling, bath: f64) -> Fixture {
    let mut f = fixture(false, false);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let reference = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
        root.join("data/reference-project/cooling-radiation.fsim")).unwrap()).unwrap().decoded.spec;
    let mut radiation = reference.cooling.unwrap().conduction.unwrap().radiation.unwrap();
    radiation.surfaces[0].name = "cold-radiator".into();
    radiation.surfaces[0].target = "cold".into();
    radiation.surfaces[0].reservoir_temperature.value = bath;
    radiation.temperature_tolerance.value = 1e-11;
    radiation.heat_tolerance.value = 1e-10;
    radiation.max_iterations = 256;
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.radiation = Some(radiation);
    for boundary in &mut setup.boundaries {
        if boundary.target == "cold" {
            boundary.condition = match cooling {
                Cooling::Prescribed => ThermalBoundaryCondition::Convection {
                    coefficient: QtyAny::new(2.0, dims::HEAT_TRANSFER_COEFFICIENT),
                    reference_temperature: QtyAny::new(RESERVOIR_K, dims::TEMPERATURE),
                },
                Cooling::Natural => natural(),
                Cooling::Air | Cooling::Mixed => ThermalBoundaryCondition::AirflowConvection {
                    branch: "cold".into(), order: 0,
                    inlet_temperature: QtyAny::new(INITIAL_K, dims::TEMPERATURE),
                    hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH),
                    flow_area: QtyAny::new(0.004, dims::AREA),
                    channel_length: QtyAny::new(0.3, dims::LENGTH),
                    correlation: "convection.gnielinski".into(),
                },
            };
        } else if cooling == Cooling::Mixed { boundary.condition = natural(); }
    }
    let policy = setup.transient.as_mut().unwrap();
    policy.initial_temperature = TransientInitialTemperature::Uniform(
        QtyAny::new(START_K, dims::TEMPERATURE));
    policy.power_schedules = vec![TransientRegionPower {
        region: "hot".into(), source: "synthetic held power with exactly 10 J input".into(),
        steps: [(0.25, 20.0), (0.75, 5.0), (1.0, 10.0)].into_iter().map(|(until, watts)|
            TransientPowerStep { until: QtyAny::new(until, dims::TIME),
                watts: QtyAny::new(watts, dims::POWER) }).collect(),
    }];
    f
}

fn natural() -> ThermalBoundaryCondition {
    ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(RESERVOIR_K, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    }
}

fn run_radiating(f: &Fixture, source: &Path, contact: &Path, ledger: &Path) -> fs_cli::CommandOutput {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    fs_cli::run(vec!["--json".into(), "solve".into(), source.display().to_string(),
        ledger.display().to_string(), "--materials".into(),
        f.dir.join("solid.fsmcdpk").display().to_string(), "--materials".into(),
        root.join("data/reference-project/gray-surface.fsmcdpk").display().to_string(),
        "--interfaces".into(), contact.display().to_string()])
}

pub(super) fn solve_radiating(f: &Fixture, ordinal: usize, pack: (Vec<u8>, String, String)) -> (J, J) {
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, pack);
    let output = run_radiating(f, &source, &contact, &ledger);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", output.stdout, output.stderr);
    let output = J::parse(&output.stdout).unwrap();
    let run = artifact(&ledger, output.str_field("run_receipt").unwrap());
    let stage = run.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn radiation_heat(row: &J, bath: f64) -> (f64, f64) {
    let radiation = row.get("radiation").unwrap();
    let patches = radiation.get("patches").unwrap().as_array().unwrap();
    assert_eq!(patches.len(), 1);
    let patch = &patches[0];
    assert_eq!(patch.str_field("region"), Some("cold"));
    // The contact face is internal; the cold unit cube has five exposed m2.
    close(n(patch, "area_m2"), 5.0, 1e-12);
    let wall = n(patch, "mean_temperature_k");
    let heat = SIGMA * EMISSIVITY * 5.0 * (wall.powi(4) - bath.powi(4));
    close(n(patch, "nonlinear_heat_w"), heat, 2e-8);
    close(n(radiation, "nonlinear_radiative_out_w"), heat, 2e-8);
    assert!(n(patch, "heat_mismatch_w").abs() <= n(patch, "heat_tolerance_w"));
    assert!(n(radiation, "physical_residual_norm_j") <= n(radiation, "physical_residual_tolerance_j"));
    assert!(n(radiation, "physical_energy_residual_j").abs() <= ENERGY_TOLERANCE_J);
    close(n(radiation, "physical_dirichlet_in_w"), 0.0, 0.0);
    (heat, wall)
}

pub(super) fn reservoir_heat(row: &J, target: &str) -> f64 {
    let natural = row.get("natural").unwrap();
    let laws = natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
    assert_eq!(laws.len(), 1);
    let law = &laws[0];
    assert_eq!(law.str_field("target"), Some(target));
    assert_eq!(law.get("in_domain"), Some(&J::Bool(true)));
    close(n(law, "ambient_k"), RESERVOIR_K, 0.0);
    let wall = n(law, "mean_wall_k");
    let film = 0.5 * (wall + RESERVOIR_K);
    let nu = 1.846e-5 / (101325.0 / (287.05 * film));
    let ra = 9.80665 * (wall - RESERVOIR_K).abs() * 0.1_f64.powi(3)
        / (film * nu * (nu / 0.707));
    let factor = 0.825 + 0.387 * ra.powf(1.0 / 6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let h = 0.0263 * factor * factor / 0.1;
    let heat = 5.0 * h * (wall - RESERVOIR_K);
    close(n(law, "htc_w_m2_k"), h, 1e-8);
    close(n(law, "heat_rate_w"), heat, 3e-6);
    assert!(heat > 0.1);
    heat
}

fn air_heat(row: &J, radiative: f64, reservoir: f64) -> f64 {
    let coupled = row.get("conjugate").unwrap();
    let exchange = coupled.get("exchange").unwrap();
    assert_eq!(exchange.str_field("path"), Some("vent:cold"));
    close(n(exchange, "inlet_k"), INITIAL_K, 0.0);
    close(n(exchange.get("air_properties").unwrap(), "specific_heat_j_kg_k"), 1007.0, 0.0);
    let heat = 1007.0 * n(exchange, "mass_flow_kg_s") * (n(exchange, "outlet_k") - INITIAL_K);
    assert!(heat > 0.1);
    close(n(exchange, "air_total_w"), heat, 3e-6);
    close(n(coupled, "air_heat_gain_w"), heat, 3e-6);
    close(n(coupled, "radiative_out_w"), radiative, 3e-6);
    close(n(coupled, "off_path_convective_out_w"), reservoir, 3e-6);
    assert!((n(exchange, "air_total_w") - (heat + radiative + reservoir)).abs() > 0.1,
        "radiation and off-path reservoirs must not heat transported air");
    heat
}

#[test]
fn native_capacity_contact_radiation_preserves_heat_accounts_across_all_cooling_paths() {
    // Reuse only the existing enthalpy owner's mesh mass weights. Its pointwise
    // radiation trajectory is NOT a reference for this area-mean emission law.
    let mut mass_fixture = fixture(true, false);
    time(&mut mass_fixture).max_step.value = 1.0;
    let (mass_receipt, _) = solve(&mass_fixture, 180);
    let masses = values(mass_receipt.get("transient").unwrap().get("endpoint").unwrap(), "reference_nodal_mass_kg");
    close(masses.iter().sum(), 20.0, 1e-10);
    for (case, cooling) in [Cooling::Prescribed, Cooling::Natural, Cooling::Air, Cooling::Mixed].into_iter().enumerate() {
        let f = radiating_fixture(cooling, BATH_K);
        let mut curve_endpoint = Vec::new();
        for surface in [false, true] {
            let pack = if surface { temperature_surface_pack() } else { temperature_pack() };
            let (card, claim) = (pack.1.clone(), pack.2.clone());
            let (receipt, field) = solve_radiating(&f, 181 + 2 * case + usize::from(surface), pack);
            let trajectory = receipt.get("transient").unwrap();
            let work = trajectory.get("temperature_contact").unwrap();
            assert_eq!(trajectory.get("radiation").unwrap().str_field("method"), Some("implicit-area-mean-gray-radiation"));
            let method = if surface { "source-backed-contact-fixed-point" } else { "source-backed-contact-newton" };
            assert_eq!(work.str_field("method"), Some(method));
            let (mut responses, mut fine_energy, mut min_query, mut max_query) = (0.0, 0.0, f64::INFINITY, f64::NEG_INFINITY);
            for (grid, ends) in [("coarse", vec![0.25, 0.75, 1.0]),
                ("fine", vec![0.125, 0.25, 0.5, 0.75, 0.875, 1.0])] {
                let rows = trajectory.get(grid).unwrap().as_array().unwrap();
                assert_eq!(rows.len(), ends.len());
                let (mut start, mut input, mut integrated, mut grid_responses) = (0.0, 0.0, 0.0, 0.0);
                for (row, end) in rows.iter().zip(ends) {
                    let dt = end - start;
                    close(n(row, "time_s"), end, 0.0); close(n(row, "dt_s"), dt, 0.0);
                    let query = check_contact(row, &claim);
                    min_query = min_query.min(query); max_query = max_query.max(query);
                    let contact = row.get("temperature_contact").unwrap();
                    assert_eq!(contact.str_field("method"), Some(method));
                    for face in contact.get("face_states").unwrap().as_array().unwrap() {
                        assert_eq!(face.str_field("interface_card"), Some(card.as_str()));
                    }
                    let (radiative, wall) = radiation_heat(row, BATH_K);
                    assert!(radiative > 0.1);
                    let reservoir = match cooling {
                        Cooling::Prescribed => 5.0 * 2.0 * (wall - RESERVOIR_K),
                        Cooling::Natural => reservoir_heat(row, "cold"),
                        Cooling::Mixed => reservoir_heat(row, "hot"),
                        Cooling::Air => 0.0,
                    };
                    let air = if matches!(cooling, Cooling::Air | Cooling::Mixed) { air_heat(row, radiative, reservoir) } else { 0.0 };
                    close(n(row.get("radiation").unwrap(), "convective_out_w"), air + reservoir, 3e-6);
                    let owner = row.get(match cooling {
                        Cooling::Prescribed => "radiation", Cooling::Air => "conjugate",
                        Cooling::Natural | Cooling::Mixed => "natural",
                    }).unwrap();
                    assert!(n(owner, "physical_residual_norm_j") <= n(owner, "physical_residual_tolerance_j"));
                    for key in ["physical_residual_norm_j", "physical_energy_residual_j"] {
                        close(n(contact, key), n(owner, key), 0.0);
                    }
                    if cooling != Cooling::Prescribed {
                        close(n(contact, "thermal_responses"), n(owner, "solid_solves"), 0.0);
                        assert!(n(contact, "thermal_responses") >= n(owner, "radiation_trials"));
                        if surface { assert!(n(contact, "thermal_responses") > n(owner, "radiation_trials")); }
                    }
                    let power = if start < 0.25 { 20.0 } else if start < 0.75 { 5.0 } else { 10.0 };
                    let net = power - radiative - air - reservoir;
                    close(n(row, "source_w"), power, 1e-10);
                    close(n(row, "net_input_w"), net, 2.0 * ENERGY_TOLERANCE_J / dt + 1e-8);
                    close(n(row, "stored_energy_change_j"), dt * net, 2.0 * ENERGY_TOLERANCE_J + 1e-8);
                    grid_responses += n(contact, "thermal_responses");
                    input += dt * power; integrated += dt * net; start = end;
                }
                close(input, 10.0, 0.0);
                close(n(work, &format!("{grid}_thermal_responses")), grid_responses, 0.0);
                responses += grid_responses;
                if grid == "fine" { fine_energy = integrated; }
            }
            close(n(work, "total_thermal_responses"), responses, 0.0);
            assert!(max_query - min_query > 1e-4, "{cooling:?}: source law must follow cooling");
            let temperatures = values(&field, "temperature");
            assert_eq!(temperatures.len(), masses.len());
            let stored: f64 = temperatures.iter().zip(&masses).map(|(t, mass)| mass * CP * (t - START_K)).sum();
            close(stored, fine_energy, 8.0 * ENERGY_TOLERANCE_J);
            let energy = trajectory.get("energy").unwrap();
            close(n(energy, "stored_change_j"), stored, 1e-8);
            close(n(energy, "integrated_net_input_j"), fine_energy, 4e-6);
            assert!(stored < 0.0, "the baths remove stored heat despite 10 J supplied power");
            // Final receipt work belongs to the last inner radiation solve;
            // natural/air contact rows also count discarded outer trials.
            let last = trajectory.get("fine").unwrap().as_array().unwrap().last().unwrap();
            let final_radiation = receipt.get("radiation").unwrap();
            close(n(final_radiation, "radiation_trials"), n(last.get("radiation").unwrap(), "iterations"), 0.0);
            assert!(n(final_radiation, "solid_solves") >= n(final_radiation, "radiation_trials"));
            if surface { assert!(n(final_radiation, "solid_solves") > n(final_radiation, "radiation_trials")); }
            if cooling == Cooling::Prescribed {
                close(n(final_radiation, "solid_solves"), n(last.get("temperature_contact").unwrap(), "thermal_responses"), 0.0);
            }
            if surface {
                for (actual, expected) in temperatures.iter().zip(&curve_endpoint) { close(*actual, *expected, 3e-6); }
            } else { curve_endpoint = temperatures; }
        }
    }
}

#[test]
fn native_capacity_radiation_contact_refuses_late_source_exit_and_retries_original_history() {
    let mut f = radiating_fixture(Cooling::Prescribed, 293.0);
    for boundary in &mut f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        boundary.condition = ThermalBoundaryCondition::HeatFlux { outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX) };
    }
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 293.0; envelope.ambient_hi.value = 293.25;
    let requirement = &mut f.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 293.25; requirement.margin.value = 0.0;
    let policy = time(&mut f);
    policy.initial_temperature = TransientInitialTemperature::Uniform(QtyAny::new(INITIAL_K, dims::TEMPERATURE));
    policy.max_step.value = 0.25;
    policy.power_schedules[0].steps = [(0.5, 0.0), (1.0, 200.0)].into_iter().map(|(until, watts)|
        TransientPowerStep { until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER) }).collect();
    let narrow = temperature_pack_from_knots(vec![(293.0, 0.293), (293.3, 0.2933)]);
    let (source, contact, ledger) = prepare_with_pack(&f, 190, None, narrow);
    let output = run_radiating(&f, &source, &contact, &ledger);
    assert_ne!(output.exit_code, fs_cli::exit::SUCCESS);
    let refusal = format!("{}{}", output.stdout, output.stderr);
    assert!(refusal.contains("cli-solve-conduction-transient"), "{refusal}");
    assert!(!refusal.contains("run_receipt"), "unsupported source state published history: {refusal}");
    let pack = temperature_pack(); let claim = pack.2.clone();
    let (receipt, _) = solve_radiating(&f, 191, pack);
    for grid in ["coarse", "fine"] {
        let (mut supplied, mut crossed) = (0.0, false);
        for row in receipt.get("transient").unwrap().get(grid).unwrap().as_array().unwrap() {
            let query = check_contact(row, &claim);
            let (radiated, _) = radiation_heat(row, 293.0);
            let power = if n(row, "time_s") <= 0.5 { 0.0 } else { 200.0 };
            close(n(row.get("radiation").unwrap(), "convective_out_w"), 0.0, 0.0);
            close(n(row, "source_w"), power, 1e-10);
            close(n(row, "stored_energy_change_j"), n(row, "dt_s") * (power - radiated), 2.0 * ENERGY_TOLERANCE_J);
            if power == 0.0 {
                assert!(n(row, "stored_energy_change_j") < 0.0);
                assert!(query < INITIAL_K, "retry must restart the original cooling history");
            }
            supplied += n(row, "dt_s") * power;
            crossed |= query > 293.3;
        }
        close(supplied, 100.0, 1e-10);
        assert!(crossed, "the wider source must admit the previously unsupported endpoint query");
    }
}
