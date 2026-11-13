//! Full physical enthalpy histories retain latent storage and live source laws.
use super::*;
#[path = "adjoint/radiation.rs"]
mod radiation;

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";
const DH: &str = "dtemperature_dinitial_specific_enthalpy_k_kg_per_j";
const DRHO: &str = "dtemperature_dreference_density_k_m3_per_kg";
const DK: &str = "dtemperature_dconductivity_multiplier_k";
const DR: &str = "dtemperature_dresistance_multiplier_k";

fn policy(f: &mut Fixture) -> &mut ConductionTransient {
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap()
}

fn source_pack(surface: bool, scale: f64) -> (Vec<u8>, String, String) {
    if !surface {
        return temperature_pack_from_knots(vec![(200.0, 0.2 * scale), (450.0, 0.45 * scale)]);
    }
    // Same R(T) slice at the original 1 MPa, with an independent pressure
    // axis. All accepted temperatures stay below the triangle seam at 325 K.
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: [[200.0, 0.5e6], [450.0, 0.5e6], [200.0, 1.5e6], [450.0, 1.5e6]]
            .into_iter().map(|point| (point,
                scale * (0.001 * point[0] + 1e-8 * (point[1] - 1e6)))).collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

fn scaled_material(scale: f64) -> (Vec<u8>, String, String) {
    let (bytes, _, _) = material_pack(true);
    let original = NormalizedMaterialCardPack::from_bytes(&bytes).unwrap();
    let claim = original.card().claims().claims_ordered().next().unwrap().1;
    let PropertyValue::Curve { abscissa, abscissa_dims, knots, dims } = &claim.value
        else { panic!("original nonlinear conductivity source"); };
    let claims = pack(original.pack_id(), "thermal-conductivity", *dims,
        PropertyValue::Curve { abscissa: abscissa.clone(), abscissa_dims: *abscissa_dims,
            knots: knots.iter().map(|&(temperature, value)| (temperature, scale * value)).collect(),
            dims: *dims }, claim.interpolation);
    let scaled = NormalizedMaterialCardPack::new(original.card().id().clone(), claims).unwrap();
    (scaled.to_bytes(), scaled.card().content_hash().to_hex(), scaled.card().id().to_string())
}

/// Re-author the complete R and k curves for each physical perturbation. The
/// enthalpy declaration retains the matching conductivity-card identity; its
/// chart and density do not change under a conductivity control.
fn solve_history(f: &Fixture, ordinal: usize, output: Option<&str>, surface: bool,
    resistance: f64, conductivity: f64) -> (J, J)
{
    let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
    let mut materials = vec![f.dir.join("solid.fsmcdpk")];
    if conductivity != 1.0 {
        let (bytes, card, state) = scaled_material(conductivity);
        let path = f.dir.join(format!("enthalpy-hot-{ordinal}.fsmcdpk"));
        std::fs::write(&path, bytes).unwrap();
        materials.push(path);
        let binding = changed.project.materials.as_mut().unwrap().iter_mut()
            .find(|row| row.region == "hot").unwrap();
        binding.card = card.clone();
        binding.state = state;
        binding.claim = None;
        policy(&mut changed).enthalpy.as_mut().unwrap().materials.iter_mut()
            .find(|row| row.region == "hot").unwrap().material_card_identity =
            fs_blake3::ContentHash::from_hex(&card).unwrap();
    }
    let (source, contact, ledger) = prepare_fixed_pressure_surface(&changed, ordinal, output, source_pack(surface, resistance));
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

fn report(receipt: &J) -> &J {
    let result = receipt.get("transient").unwrap().get("adjoint").unwrap();
    assert_eq!(result.str_field("schema"), Some("fs-cli-native-enthalpy-adjoint-v1"));
    result
}

fn row<'a>(report: &'a J, list: &str, key: &str, identity: &str) -> &'a J {
    let matches: Vec<_> = report.get(list).unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field(key) == Some(identity)).collect();
    assert_eq!(matches.len(), 1, "one {list}/{identity} physical control");
    matches[0]
}

fn source_multiplier(report: &J) -> f64 {
    report.get("intervals").unwrap().as_array().unwrap().iter()
        .map(|interval| n(interval, "dtemperature_dsource_multiplier_k")).sum()
}

fn unchanged(plain: &J, plain_field: &J, differentiated: &J, field: &J) {
    assert_eq!(plain_field.get("temperature"), field.get("temperature"));
    for key in ["energy", "solver", "interfaces"] {
        assert_eq!(
            forward_section(plain, key),
            forward_section(differentiated, key),
            "forward {key} changed"
        );
    }
    let plain = plain.get("transient").unwrap();
    let differentiated = differentiated.get("transient").unwrap();
    for key in ["coarse", "fine", "energy", "temperature_contact",
        "materials", "endpoint", "sampled_peak", "temporal_error"] {
        assert_eq!(plain.get(key), differentiated.get(key), "accepted enthalpy history {key} changed");
    }
    let plain_nonlinear = plain.get("nonlinear").unwrap();
    let differentiated_nonlinear = differentiated.get("nonlinear").unwrap();
    let forward_fields = |nonlinear: &J| {
        let J::Object(fields) = nonlinear else { panic!("original nonlinear report"); };
        fields.iter().filter(|(key, _)| key != "total_forward_and_adjoint_krylov_limit")
            .cloned().collect::<Vec<_>>()
    };
    assert_eq!(forward_fields(plain_nonlinear), forward_fields(differentiated_nonlinear),
        "accepted forward nonlinear work or controls changed");
    let gradient = differentiated.get("adjoint").unwrap();
    assert_eq!(n(plain_nonlinear, "total_forward_and_adjoint_krylov_limit"),
        n(plain_nonlinear, "total_krylov_limit"));
    assert_eq!(n(differentiated_nonlinear, "total_forward_and_adjoint_krylov_limit"),
        n(differentiated_nonlinear, "total_krylov_limit") + n(gradient, "max_adjoint_krylov_iterations"));
    assert_eq!(n(gradient, "revalidated_solid_endpoints"), n(gradient, "state_index"));
    assert!(n(gradient, "maximum_dual_relative_residual") < 1e-10);
    assert!(n(gradient, "adjoint_krylov_iterations") <= n(gradient, "max_adjoint_krylov_iterations"));
}

fn check_source(receipt: &J, surface: bool) {
    let (_, card, claim) = source_pack(surface, 1.0);
    let time = receipt.get("transient").unwrap();
    let method = if surface { "source-backed-contact-enthalpy-fixed-point" }
        else { "source-backed-contact-enthalpy-newton" };
    assert_eq!(time.get("temperature_contact").unwrap().str_field("method"), Some(method));
    for grid in ["coarse", "fine"] {
        for step in time.get(grid).unwrap().as_array().unwrap() {
            let contact = step.get("temperature_contact").unwrap();
            assert_eq!(contact.str_field("method"), Some(method));
            assert!(n(contact, "physical_residual_norm_j") <= n(contact, "physical_residual_tolerance_j"));
            assert!(n(contact, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j"));
            for state in contact.get("face_states").unwrap().as_array().unwrap() {
                assert_eq!(state.str_field("interface_card"), Some(card.as_str()));
                assert_eq!(state.str_field("selected_claim"), Some(claim.as_str()));
                close(n(state, "resistance_m2_k_w"), 0.001 * n(state, "query_temperature_k"), 1e-13);
            }
        }
    }
    let material = report(receipt).get("material_controls").unwrap();
    assert_eq!(n(material, "control_steps"), n(report(receipt), "state_index"));
    let contact = row(material, "contact_interfaces", "interface", "cold-hot-joint");
    assert_eq!(contact.str_field("interface_card"), Some(card.as_str()));
    assert_eq!(contact.get("mapped"), Some(&J::Bool(true)));
    assert!(n(contact, "face_pairs") > 0.0);
    assert_eq!(n(contact, "reference_value"), 1.0);
    assert_eq!(n(contact, "log_reference_value"), 0.0);
}

fn future_is_zero(gradient: &J) {
    let selected = n(gradient, "state_index") as usize;
    let intervals = gradient.get("intervals").unwrap().as_array().unwrap();
    assert!(selected < intervals.len());
    for interval in &intervals[selected..] {
        assert_eq!(n(interval, "dtemperature_dsource_multiplier_k"), 0.0);
        for source in interval.get("regional_power").unwrap().as_array().unwrap() {
            assert_eq!(n(source, "dtemperature_ddelivered_watts_k_per_w"), 0.0);
        }
    }
}

#[test]
fn live_enthalpy_contact_final_and_peak_controls_match_complete_latent_histories() {
    for peak in [false, true] {
        let mut f = fixture(2.0, 0.25);
        f.project.power.as_mut().unwrap()[0].duty = 1.0;
        if peak {
            let schedule = &mut policy(&mut f).power_schedules[0];
            schedule.steps[0].watts.value = 80.0;
            schedule.steps[1].watts.value = 0.0;
        }
        let surface = peak;
        let output = if peak { PEAK } else { FINAL };
        let (plain, original) = solve_history(&f, 500, None, surface, 1.0, 1.0);
        let (receipt, field) = solve_history(&f, 501, Some(output), surface, 1.0, 1.0);
        unchanged(&plain, &original, &receipt, &field);
        check_source(&receipt, surface);
        let time = receipt.get("transient").unwrap();
        let endpoint = time.get("endpoint").unwrap();
        close(values(endpoint, "reference_nodal_mass_kg").iter().sum(), 3.0, 1e-12);
        assert!(values(endpoint, "liquid_mass_fraction").iter().any(|&fraction|
            fraction > 1e-5 && fraction < 1.0 - 1e-5), "the history must retain a real latent state");
        let gradient = report(&receipt);
        assert_eq!(gradient.str_field("output"), Some(output));
        let selected = n(gradient, "active_vertex") as usize;
        if peak {
            let sample = time.get("sampled_peak").unwrap();
            assert_eq!(n(gradient, "value_k"), n(sample, "value_k"));
            assert_eq!(n(gradient, "time_s"), n(sample, "time_s"));
            assert!(n(gradient, "time_s") > 0.0 && n(gradient, "time_s") < 2.0);
            assert!(n(sample, "value_k") > n(time.get("temporal_error").unwrap(), "fine_final_k") + 1e-3);
            future_is_zero(gradient);
        } else {
            assert_eq!(n(gradient, "state_index"), time.get("fine").unwrap().as_array().unwrap().len() as f64);
            assert_eq!(n(gradient, "value_k"), values(&field, "temperature")[selected]);
        }
        let cold = row(gradient, "regional_storage", "region", "cold");
        let hot = row(gradient, "regional_storage", "region", "hot");
        assert_eq!(n(cold, "reference_density_kg_m3"), 1.0);
        assert_eq!(n(hot, "reference_density_kg_m3"), 2.0);
        close(n(cold, DH) + n(hot, DH),
            n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"), 1e-12);
        close(n(cold, DRHO) + 2.0 * n(hot, DRHO),
            n(gradient, "dtemperature_dreference_density_multiplier_k"), 1e-12);
        let material = gradient.get("material_controls").unwrap();
        let conductivity = row(material, "conductivity_regions", "region", "hot");
        let joint = row(material, "contact_interfaces", "interface", "cold-hot-joint");
        let controls = [("source", source_multiplier(gradient), 0.002),
            ("cold-h", n(cold, DH), 0.02), ("hot-h", n(hot, DH), 0.02),
            ("hot-density", n(hot, DRHO), 0.002),
            ("hot-k", n(conductivity, DK), 0.002), ("joint", n(joint, DR), 0.002)];
        for (coordinate, (name, actual, step)) in controls.into_iter().enumerate() {
            assert!(actual.abs() > 1e-7, "active enthalpy {name} control, peak={peak}");
            let mut probes = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                let delta = sign * step;
                match name {
                    "source" => {
                        for source in changed.project.power.as_mut().unwrap() { source.watts.value *= 1.0 + delta; }
                        for schedule in &mut policy(&mut changed).power_schedules {
                            for step in &mut schedule.steps { step.watts.value *= 1.0 + delta; }
                        }
                    }
                    "cold-h" | "hot-h" | "hot-density" => {
                        let region = if name == "cold-h" { "cold" } else { "hot" };
                        let storage = policy(&mut changed).enthalpy.as_mut().unwrap().materials.iter_mut()
                            .find(|row| row.region == region).unwrap();
                        if name == "hot-density" { storage.reference_density.value += delta; }
                        else { storage.initial_specific_enthalpy.value += delta; }
                    }
                    "hot-k" | "joint" => {},
                    _ => unreachable!(),
                }
                let (physical, field) = solve_history(&changed, 502 + 2 * coordinate + side, None, surface,
                    1.0 + if name == "joint" { delta } else { 0.0 },
                    1.0 + if name == "hot-k" { delta } else { 0.0 });
                if peak {
                    let sample = physical.get("transient").unwrap().get("sampled_peak").unwrap();
                    assert_eq!(n(sample, "state"), n(gradient, "state_index"), "physical peak branch changed");
                    assert_eq!(n(sample, "time_s"), n(gradient, "time_s"));
                    probes.push(n(sample, "value_k"));
                } else {
                    probes.push(values(&field, "temperature")[selected]);
                }
            }
            let expected = (probes[1] - probes[0]) / (2.0 * step);
            close(actual, expected, 8e-5 * expected.abs().max(0.01));
        }
    }
}

#[test]
fn live_enthalpy_initial_peak_retains_material_identities_with_zero_future_controls() {
    let mut f = fixture(2.0, 0.25);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    policy(&mut f).power_schedules.clear();
    let (plain, original) = solve_history(&f, 520, None, true, 1.0, 1.0);
    let (receipt, field) = solve_history(&f, 521, Some(PEAK), true, 1.0, 1.0);
    unchanged(&plain, &original, &receipt, &field);
    check_source(&receipt, true);
    let gradient = report(&receipt);
    assert_eq!(n(gradient, "state_index"), 0.0);
    assert_eq!(n(gradient, "time_s"), 0.0);
    close(n(gradient, "value_k"), 305.0, 1e-12);
    assert_eq!(gradient.get("unique_sampled_maximum"), Some(&J::Bool(false)));
    assert_eq!(n(gradient, "revalidated_solid_endpoints"), 0.0);
    assert_eq!(n(gradient, "adjoint_krylov_iterations"), 0.0);
    assert_eq!(n(gradient, "dtemperature_dreference_density_multiplier_k"), 0.0);
    let active = n(gradient, "active_vertex") as usize;
    for (vertex, derivative) in values(gradient, "dtemperature_dinitial_specific_enthalpy_k_kg_per_j")
        .into_iter().enumerate() {
        close(derivative, if vertex == active { 0.1 } else { 0.0 }, 1e-14);
    }
    close(n(row(gradient, "regional_storage", "region", "hot"), DH), 0.1, 1e-14);
    assert_eq!(n(row(gradient, "regional_storage", "region", "cold"), DH), 0.0);
    for storage in gradient.get("regional_storage").unwrap().as_array().unwrap() {
        assert_eq!(n(storage, DRHO), 0.0);
    }
    let material = gradient.get("material_controls").unwrap();
    for region in material.get("conductivity_regions").unwrap().as_array().unwrap() {
        assert_eq!(n(region, DK), 0.0);
    }
    for joint in material.get("contact_interfaces").unwrap().as_array().unwrap() {
        assert_eq!(n(joint, DR), 0.0);
    }
    future_is_zero(gradient);
}
