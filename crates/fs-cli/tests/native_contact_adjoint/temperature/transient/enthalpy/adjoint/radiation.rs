//! Local emission and reflected heat share the original live-contact h history.
use super::*;

const EPSILON: f64 = 0.73;

fn radiating_fixture(open: bool) -> Fixture {
    let mut f = fixture(0.4, 0.1);
    super::super::radiation(&mut f, if open { 3 } else { 1 });
    f.project.power.as_mut().unwrap()[0].duty = 1.0;
    if open {
        // Enough input to rise above the initial hot state despite local
        // emission, then a source-free cooling interval before the horizon.
        let schedule = &mut policy(&mut f).power_schedules[0];
        schedule.steps[0].watts.value = 200.0;
        schedule.steps[1].watts.value = 0.0;
    }
    f
}

fn gray_source(emissivity: f64) -> (Vec<u8>, String) {
    let name = "pointwise-contact-gray";
    let bytes = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: name.into(), phase: "solid".into(), process: "synthetic".into(), revision: 0,
    }, pack(name, fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_qty::Dims::NONE,
        PropertyValue::Scalar { value: emissivity, dims: fs_qty::Dims::NONE },
        InterpolationPolicy::ConstantWithinValidity)).unwrap().to_bytes();
    let admitted = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Material, source: name.into(), bytes: bytes.clone(), expect: None,
    }]).unwrap();
    (bytes, admitted.materials()[0].card().to_hex())
}

fn solve_radiating(f: &Fixture, ordinal: usize, output: Option<&str>, open: bool,
    emissivities: [f64; 2], resistance: f64) -> (J, J)
{
    let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
    let mut materials = vec![f.dir.join("solid.fsmcdpk")];
    let mut cards: Vec<String> = Vec::new();
    for index in 0..if open { 2 } else { 1 } {
        if index == 1 && emissivities[1] == emissivities[0] {
            cards.push(cards[0].clone());
            continue;
        }
        let (bytes, card) = gray_source(emissivities[index]);
        let path = f.dir.join(format!("pointwise-gray-{ordinal}-{index}.fsmcdpk"));
        std::fs::write(&path, bytes).unwrap();
        materials.push(path);
        cards.push(card);
    }
    let setup = changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    if open {
        for (surface, card) in setup.enclosure_radiation.as_mut().unwrap().surfaces.iter_mut().zip(cards) {
            surface.card = card;
            surface.claim = None;
        }
    } else {
        setup.radiation.as_mut().unwrap().surfaces[0].card = cards[0].clone();
        setup.radiation.as_mut().unwrap().surfaces[0].claim = None;
    }
    let (source, contact, ledger) = prepare_fixed_pressure_surface(&changed, ordinal, output, source_pack(open, resistance));
    let mut args = vec!["--json".to_string(), "solve".into(),
        source.display().to_string(), ledger.display().to_string(),
        "--interfaces".into(), contact.display().to_string()];
    for path in materials { args.extend(["--materials".into(), path.display().to_string()]); }
    let result = fs_cli::run(args);
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", result.stdout, result.stderr);
    let result = J::parse(&result.stdout).unwrap();
    let completed = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn radiative_controls(receipt: &J, open: bool) -> &J {
    let gradient = report(receipt);
    let radiation = gradient.get("radiation").unwrap();
    assert_eq!(radiation.str_field("method"), Some(if open {
        "complete-pointwise-p1-open-gray-diffuse-enthalpy-tangent"
    } else { "complete-pointwise-p1-gray-enthalpy-tangent" }));
    assert_eq!(radiation.str_field("spatial_model"), Some("pointwise"));
    assert_eq!(radiation.str_field("residual_target"), Some("original-accepted-forward-physical-tolerance"));
    assert_eq!(n(radiation, "control_steps"), n(gradient, "state_index"));
    assert_eq!(n(radiation, "parameter_linear_solves"), 0.0);
    assert_eq!(n(radiation, "parameter_krylov_iterations"), 0.0);
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), if open { 2 } else { 1 });
    for (index, surface) in surfaces.iter().enumerate() {
        assert_eq!(surface.str_field("target"), Some(if index == 0 { "cold" } else { "hot" }));
        assert_eq!(n(surface, "emissivity"), EPSILON);
        if open {
            assert!(surface.get("reservoir_temperature_k").is_none());
            assert!(surface.get("dtemperature_dreservoir_temperature_k_per_k").is_none());
        }
    }
    if open {
        assert_eq!(surfaces[0].get("card"), surfaces[1].get("card"),
            "independent physical surfaces start from one shared source");
        assert_ne!(surfaces[0].get("name"), surfaces[1].get("name"));
        assert_eq!(n(radiation.get("black_environment").unwrap(), "temperature_k"), 300.0);
    } else {
        assert!(radiation.get("black_environment").is_none());
        assert_eq!(n(&surfaces[0], "reservoir_temperature_k"), 300.0);
    }
    radiation
}

fn bath_derivative(radiation: &J, open: bool) -> f64 {
    if open {
        n(radiation.get("black_environment").unwrap(), "dtemperature_dambient_temperature_k_per_k")
    } else {
        n(row(radiation, "surfaces", "target", "cold"), "dtemperature_dreservoir_temperature_k_per_k")
    }
}

fn check_heat(receipt: &J, open: bool) {
    let time = receipt.get("transient").unwrap();
    let tolerance = n(time, "energy_tolerance_j");
    let mut fine_net = 0.0;
    let mut exchanged = false;
    for grid in ["coarse", "fine"] {
        for step in time.get(grid).unwrap().as_array().unwrap() {
            let radiation = step.get("radiation").unwrap();
            assert_eq!(radiation.str_field("model"), Some("pointwise-p1"));
            assert!(n(radiation, "physical_residual_j") <= n(radiation, "physical_residual_threshold_j"));
            assert!(n(radiation, "physical_energy_residual_j").abs() <= tolerance);
            assert_eq!(n(radiation, "convective_out_w"), 0.0);
            assert_eq!(n(radiation, "prescribed_neumann_out_w"), 0.0);
            assert_eq!(n(radiation, "physical_dirichlet_in_w"), 0.0);
            let radiative = n(radiation, "nonlinear_radiative_out_w");
            exchanged |= radiative.abs() > 1e-3;
            let net = n(step, "source_w") - radiative;
            close(n(step, "net_input_w"), net, 2e-6);
            close(n(step, "stored_energy_change_j"), n(step, "dt_s") * net, tolerance);
            if grid == "fine" { fine_net += n(step, "dt_s") * net; }
        }
    }
    assert!(exchanged, "the radiative reservoir must exchange actual stored energy");
    close(stored(time), INITIAL_ENERGY_J + fine_net,
        time.get("fine").unwrap().as_array().unwrap().len() as f64 * tolerance);
    let radiation = receipt.get("radiation").unwrap();
    assert_eq!(radiation.str_field("spatial_model"), Some("pointwise"));
    if open {
        close(n(radiation, "ambient_absorbed_w"), n(radiation, "radiative_out_w"), 1e-6);
    }
}

#[test]
fn live_enthalpy_pointwise_ambient_and_open_peak_controls_match_full_physical_histories() {
    for open in [false, true] {
        let f = radiating_fixture(open);
        let output = if open { PEAK } else { FINAL };
        let (plain, original) = solve_radiating(&f, 540, None, open, [EPSILON; 2], 1.0);
        let (receipt, field) = solve_radiating(&f, 541, Some(output), open, [EPSILON; 2], 1.0);
        unchanged(&plain, &original, &receipt, &field);
        assert_eq!(plain.get("radiation"), receipt.get("radiation"));
        assert_eq!(plain.get("transient").unwrap().get("radiation"),
            receipt.get("transient").unwrap().get("radiation"));
        check_source(&receipt, open);
        check_heat(&receipt, open);
        let time = receipt.get("transient").unwrap();
        assert!(values(time.get("endpoint").unwrap(), "liquid_mass_fraction").iter()
            .any(|&fraction| fraction > 1e-5 && fraction < 1.0 - 1e-5));
        let gradient = report(&receipt);
        let radiation = radiative_controls(&receipt, open);
        let selected = n(gradient, "active_vertex") as usize;
        if open {
            assert!(n(gradient, "time_s") > 0.0 && n(gradient, "time_s") < 0.4);
            assert_eq!(n(gradient, "value_k"), n(time.get("sampled_peak").unwrap(), "value_k"));
            assert!(n(gradient, "value_k") > n(time.get("temporal_error").unwrap(), "fine_final_k") + 1e-3);
            future_is_zero(gradient);
        } else {
            assert_eq!(n(gradient, "value_k"), values(&field, "temperature")[selected]);
        }
        let joint = row(gradient.get("material_controls").unwrap(),
            "contact_interfaces", "interface", "cold-hot-joint");
        let mut controls = vec![("cold-emissivity", n(row(radiation, "surfaces", "target", "cold"),
            "dtemperature_demissivity_k"), 0.0005),
            ("bath", bath_derivative(radiation, open), 0.01),
            ("source", source_multiplier(gradient), 0.002), ("joint", n(joint, DR), 0.002),
            ("hot-h", n(row(gradient, "regional_storage", "region", "hot"), DH), 0.02)];
        if open {
            controls.push(("hot-emissivity", n(row(radiation, "surfaces", "target", "hot"),
                "dtemperature_demissivity_k"), 0.0005));
        }
        for (coordinate, (name, actual, step)) in controls.into_iter().enumerate() {
            assert!(actual.abs() > 1e-7, "active pointwise {name} control, open={open}");
            let mut probes = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                let delta = sign * step;
                let mut emissivities = [EPSILON; 2];
                match name {
                    "cold-emissivity" => emissivities[0] += delta,
                    "hot-emissivity" => emissivities[1] += delta,
                    "bath" => {
                        let setup = changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
                        if open {
                            setup.enclosure_radiation.as_mut().unwrap().ambient_temperature.as_mut().unwrap().value += delta;
                        } else { setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += delta; }
                    }
                    "source" => {
                        for source in changed.project.power.as_mut().unwrap() { source.watts.value *= 1.0 + delta; }
                        for schedule in &mut policy(&mut changed).power_schedules {
                            for step in &mut schedule.steps { step.watts.value *= 1.0 + delta; }
                        }
                    }
                    "hot-h" => policy(&mut changed).enthalpy.as_mut().unwrap().materials.iter_mut()
                        .find(|row| row.region == "hot").unwrap().initial_specific_enthalpy.value += delta,
                    "joint" => {},
                    _ => unreachable!(),
                }
                let (physical, field) = solve_radiating(&changed, 542 + 2 * coordinate + side, None,
                    open, emissivities, 1.0 + if name == "joint" { delta } else { 0.0 });
                if open {
                    let sample = physical.get("transient").unwrap().get("sampled_peak").unwrap();
                    assert_eq!(n(sample, "state"), n(gradient, "state_index"), "physical peak branch changed");
                    assert_eq!(n(sample, "time_s"), n(gradient, "time_s"));
                    probes.push(n(sample, "value_k"));
                } else { probes.push(values(&field, "temperature")[selected]); }
            }
            let expected = (probes[1] - probes[0]) / (2.0 * step);
            assert!((actual - expected).abs() <= 8e-5 * expected.abs().max(0.01),
                "open={open}, {name}: complete pointwise history {actual:e} versus physical resolve {expected:e}");
        }
    }
}

#[test]
fn live_enthalpy_initial_open_peak_has_zero_future_emissivity_and_shared_bath_controls() {
    let mut f = radiating_fixture(true);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    policy(&mut f).power_schedules.clear();
    let (plain, original) = solve_radiating(&f, 560, None, true, [EPSILON; 2], 1.0);
    let (receipt, field) = solve_radiating(&f, 561, Some(PEAK), true, [EPSILON; 2], 1.0);
    unchanged(&plain, &original, &receipt, &field);
    assert_eq!(plain.get("radiation"), receipt.get("radiation"));
    check_source(&receipt, true);
    let gradient = report(&receipt);
    assert_eq!(n(gradient, "state_index"), 0.0);
    close(n(gradient, "value_k"), 305.0, 1e-12);
    assert_eq!(n(gradient, "adjoint_krylov_iterations"), 0.0);
    let radiation = radiative_controls(&receipt, true);
    assert_eq!(bath_derivative(radiation, true), 0.0);
    for surface in radiation.get("surfaces").unwrap().as_array().unwrap() {
        assert_eq!(n(surface, "dtemperature_demissivity_k"), 0.0);
    }
    assert_eq!(n(row(gradient.get("material_controls").unwrap(), "contact_interfaces",
        "interface", "cold-hot-joint"), DR), 0.0);
    future_is_zero(gradient);
}
