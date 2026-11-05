//! Derived boundary controls retain the same complete live-contact history.
use super::*;

const START: f64 = 300.0;
const INLET: f64 = 293.15;
const NATURAL: &str = "convection.churchill-chu-vertical-plate";

fn cooling_fixture(mixed: bool, radiating: bool) -> Fixture {
    let mut f = fixture(true, 0.6);
    f.project.power.as_mut().unwrap()[0].watts.value = 10.0;
    f.project.power.as_mut().unwrap()[0].duty = 1.0;
    // Every physical perturbation retains these air-property/hydraulic inputs.
    // The independent inlet and natural reservoir remain inside this envelope.
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 292.0;
    envelope.ambient_hi.value = 298.0;
    let natural = |ambient| ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
        correlation: NATURAL.into(),
    };
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = natural(INLET);
    if mixed {
        setup.boundaries[0].condition = ThermalBoundaryCondition::AirflowConvection {
            branch: "cold".into(), order: 0,
            inlet_temperature: QtyAny::new(INLET, dims::TEMPERATURE),
            hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH),
            flow_area: QtyAny::new(0.004, dims::AREA),
            channel_length: QtyAny::new(0.3, dims::LENGTH),
            correlation: "convection.gnielinski".into(),
        };
        setup.boundaries[1].condition = natural(297.0);
    }
    if radiating {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let reference = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(
            root.join("data/reference-project/cooling-radiation.fsim")).unwrap()).unwrap().decoded.spec;
        let mut radiation = reference.cooling.unwrap().conduction.unwrap().radiation.unwrap();
        radiation.surfaces[0].name = "cold-radiator".into();
        radiation.surfaces[0].target = "cold".into();
        radiation.surfaces[0].reservoir_temperature.value = 290.0;
        radiation.temperature_tolerance.value = 1e-11;
        radiation.heat_tolerance.value = 1e-10;
        radiation.max_iterations = 256;
        setup.radiation = Some(radiation);
    }
    let policy = time(&mut f);
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(START, dims::TEMPERATURE));
    policy.horizon.value = 1.0;
    for capacity in &mut policy.capacities { capacity.volumetric_heat_capacity.value = 100.0; }
    policy.power_schedules[0].steps = [(0.25, 20.0), (0.75, 5.0), (1.0, 10.0)]
        .into_iter().map(|(until, watts)| TransientPowerStep {
            until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER),
        }).collect();
    f
}

fn run_cooling(f: &Fixture, ordinal: usize, output: Option<&str>, surface: bool)
    -> (fs_cli::CommandOutput, PathBuf)
{
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, output, contact_pack(surface));
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let result = fs_cli::run(vec!["--json".into(), "solve".into(),
        source.display().to_string(), ledger.display().to_string(),
        "--materials".into(), f.dir.join("solid.fsmcdpk").display().to_string(),
        "--materials".into(), root.join("data/reference-project/gray-surface.fsmcdpk").display().to_string(),
        "--interfaces".into(), contact.display().to_string()]);
    (result, ledger)
}

fn solve_cooling(f: &Fixture, ordinal: usize, output: Option<&str>, surface: bool) -> (J, J) {
    let (result, ledger) = run_cooling(f, ordinal, output, surface);
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", result.stdout, result.stderr);
    let result = J::parse(&result.stdout).unwrap();
    let completed = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "physical cooling history {actual:e} versus {expected:e}; tolerance {tolerance:e}");
}

fn control<'a>(derivative: &'a J, list: &str, key: &str, name: &str) -> &'a J {
    let rows: Vec<_> = derivative.get("convection").unwrap().get(list).unwrap()
        .as_array().unwrap().iter().filter(|row| row.str_field(key) == Some(name)).collect();
    assert_eq!(rows.len(), 1, "one physical {list}/{name} control");
    rows[0]
}

fn check_physical_history(receipt: &J, mixed: bool, radiating: bool) {
    let transient = receipt.get("transient").unwrap();
    let tolerance = n(transient, "energy_tolerance_j");
    for grid in ["coarse", "fine"] {
        for row in transient.get(grid).unwrap().as_array().unwrap() {
            let contact = row.get("temperature_contact").unwrap();
            assert!(n(contact, "physical_residual_norm_j") <= n(contact, "physical_residual_tolerance_j"));
            assert!(n(contact, "physical_energy_residual_j").abs() <= tolerance);
            let natural = row.get("natural").unwrap();
            let laws = natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
            assert_eq!(laws.len(), 1);
            let law = &laws[0];
            assert_eq!(law.str_field("target"), Some(if mixed { "hot" } else { "cold" }));
            assert_eq!(law.str_field("card"), Some(NATURAL));
            assert_eq!(law.get("in_domain"), Some(&J::Bool(true)));
            assert!(n(law, "mean_wall_k") - n(law, "ambient_k") > 0.1,
                "these physical probes must stay away from the equilibrium refusal");
            let natural_w = n(law, "heat_rate_w");
            let mut outward_w = natural_w;
            if mixed {
                let radiation_w = if radiating {
                    let heat = n(row.get("radiation").unwrap(), "nonlinear_radiative_out_w");
                    assert!(heat > 0.1);
                    heat
                } else {
                    assert_eq!(row.get("radiation"), Some(&J::Null));
                    0.0
                };
                let coupled = row.get("conjugate").unwrap();
                let exchange = coupled.get("exchange").unwrap();
                assert_eq!(exchange.str_field("branch"), Some("cold"));
                let cp = n(exchange.get("air_properties").unwrap(), "specific_heat_j_kg_k");
                let air_w = n(exchange, "mass_flow_kg_s") * cp
                    * (n(exchange, "outlet_k") - n(exchange, "inlet_k"));
                close(n(exchange, "air_total_w"), air_w, 3e-6);
                close(n(coupled, "air_heat_gain_w"), air_w, 3e-6);
                close(n(coupled, "off_path_convective_out_w"), natural_w, 3e-6);
                close(n(coupled, "radiative_out_w"), radiation_w, 3e-6);
                assert!(n(exchange, "max_region_imbalance_w").abs() <= n(exchange, "balance_tolerance_w"));
                assert!(n(coupled, "coupled_energy_residual_j").abs() <= tolerance);
                if !radiating { assert_eq!(n(coupled, "radiative_out_w"), 0.0); }
                assert!(natural_w > 0.1);
                assert!((n(exchange, "air_total_w") - (air_w + radiation_w + natural_w)).abs() > 0.1,
                    "off-path natural and radiative heat must remain outside the air stream");
                outward_w += radiation_w + air_w;
            }
            close(n(row, "net_input_w"), n(row, "source_w") - outward_w, 4e-6);
            close(n(row, "stored_energy_change_j"), n(row, "dt_s") * n(row, "net_input_w"),
                2.0 * tolerance);
        }
    }
}

#[test]
fn live_contact_capacity_natural_air_and_fan_controls_match_full_cooling_histories() {
    for (mixed, radiating) in [(false, false), (true, false), (true, true)] {
        let f = cooling_fixture(mixed, radiating);
        let surface = mixed;
        let (plain, plain_field) = solve_cooling(&f, 200, None, surface);
        let (receipt, field) = solve_cooling(&f, 201, Some(FINAL), surface);
        unchanged(&plain, &plain_field, &receipt, &field);
        assert_eq!(receipt.get("transient").unwrap().get("temperature_contact").unwrap()
            .str_field("method"), Some(if surface { "source-backed-contact-fixed-point" }
                else { "source-backed-contact-newton" }));
        for key in ["natural", "conjugate", "radiation"] {
            assert_eq!(plain.get("transient").unwrap().get(key),
                receipt.get("transient").unwrap().get(key), "forward cooling {key} changed");
        }
        check_physical_history(&receipt, mixed, radiating);
        let derivative = report(&receipt);
        assert_eq!(n(derivative, "value_k"), qoi(&plain, false));
        let convection = derivative.get("convection").unwrap();
        assert_eq!(convection.str_field("method"), Some(if mixed {
            "complete-capacity-cooling-tangent"
        } else { "complete-natural-capacity-tangent" }));
        assert_eq!(n(convection, "control_steps"), n(derivative, "state_index"));
        assert_eq!(n(convection, "natural_control_steps"), n(derivative, "state_index"));
        assert!(derivative.get("robin").unwrap().as_array().unwrap().is_empty(),
            "derived cooling must retain its physical control owners");
        let natural = control(derivative, "natural_reservoirs", "target", if mixed { "hot" } else { "cold" });
        let mut controls = vec![("natural", n(natural, "dtemperature_dambient_temperature_k_per_k"), 0.02)];
        if mixed {
            assert_eq!(n(convection, "air_control_steps"), n(derivative, "state_index"));
            let inlet = control(derivative, "air_inlets", "branch", "cold");
            let fan = control(derivative, "fan_banks", "bank", "fixture-bank");
            assert_eq!(fan.str_field("status"), Some("available"));
            controls.push(("inlet", n(inlet, "dtemperature_dinlet_temperature_k_per_k"), 0.02));
            controls.push(("fan", n(fan, "dtemperature_dspeed_ratio_k"), 0.002));
        }
        let selected = n(derivative, "active_vertex") as usize;
        for (coordinate, (name, actual, step)) in controls.into_iter().enumerate() {
            assert!(actual.abs() > 1e-7, "the {name} control must reach the hot body through this history");
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                let cooling = changed.project.cooling.as_mut().unwrap();
                if name == "fan" { cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio += sign * step; }
                for boundary in &mut cooling.conduction.as_mut().unwrap().boundaries {
                    match &mut boundary.condition {
                        ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } if name == "natural" =>
                            ambient_temperature.value += sign * step,
                        ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } if name == "inlet" =>
                            inlet_temperature.value += sign * step,
                        _ => {},
                    }
                }
                assert_eq!(changed.project.envelope, f.project.envelope);
                let (_, physical_field) = solve_cooling(&changed, 202 + 2 * coordinate + side, None, surface);
                values.push(physical_field.get("temperature").unwrap().as_array().unwrap()[selected]
                    .as_f64().unwrap());
            }
            let expected = (values[1] - values[0]) / (2.0 * step);
            close(actual, expected, 8e-5 * expected.abs().max(0.01));
        }
    }
}

#[test]
fn live_contact_capacity_cooling_keeps_natural_equilibrium_and_card_domain_refusals() {
    for equilibrium in [true, false] {
        let mut f = cooling_fixture(false, false);
        let ThermalBoundaryCondition::NaturalConvection { characteristic_length, ambient_temperature, .. } =
            &mut f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
            else { panic!("original natural wall"); };
        if equilibrium { ambient_temperature.value = START; }
        else { characteristic_length.value = 100.0; }
        if equilibrium {
            let envelope = f.project.envelope.as_mut().unwrap();
            envelope.ambient_lo.value = 299.0;
            envelope.ambient_hi.value = 301.0;
            f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
            time(&mut f).power_schedules.clear();
        }
        let (result, _) = run_cooling(&f, 210 + usize::from(equilibrium), Some(FINAL), false);
        assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
        let text = format!("{}{}", result.stdout, result.stderr);
        assert!(text.contains("natural"), "{text}");
        assert!(if equilibrium { text.contains("equilibrium") || text.contains("nonzero") }
            else { text.contains("domain") || text.contains("Rayleigh") }, "{text}");
        assert!(!text.contains("run_receipt"), "an inadmissible cooling law published a complete history: {text}");
    }
}
