//! Imported bodies compose live source contacts, reflected heat and air transport.
use super::*;
use fs_project::{ConductionEnclosure, EnclosureRadiatingSurface, EnclosureViewFactorEvidence};
use fs_project::spec::dims;
use fs_qty::QtyAny;

fn fixture() -> Fixture {
    let mut f = super::fixture(false);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let ambient = setup.radiation.take().unwrap();
    setup.enclosure_radiation = Some(ConductionEnclosure {
        spatial_model: None, ambient_temperature: None,
        surfaces: ["cold", "hot"].into_iter().map(|target| EnclosureRadiatingSurface {
            name: format!("{target}-radiator"), target: target.into(),
            card: ambient.surfaces[0].card.clone(), claim: ambient.surfaces[0].claim.clone(),
            query_temperature: QtyAny::new(300.0, dims::TEMPERATURE),
        }).collect(),
        view_factors: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        evidence: EnclosureViewFactorEvidence::Analytic {
            geometry: "synthetic equal-area exchange; not a visibility reconstruction".into(),
        },
        closure_tolerance: 1e-12, reciprocity_relative_tolerance: 1e-12,
        max_iterations: 256, temperature_tolerance: QtyAny::new(1e-10, dims::TEMPERATURE),
        heat_tolerance: QtyAny::new(1e-7, dims::POWER), relaxation: 0.05,
    });
    f
}

fn surface_pack() -> (Vec<u8>, String, String) {
    scaled_surface_pack(1.0)
}

fn scaled_surface_pack(scale: f64) -> (Vec<u8>, String, String) {
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: [([200.0, 0.8e6], 0.18), ([450.0, 0.8e6], 0.43),
            ([200.0, 1.2e6], 0.22), ([450.0, 1.2e6], 0.47)]
            .into_iter().map(|(point, value)| (point, scale * value)).collect(),
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

fn open_fixture() -> Fixture {
    let mut f = fixture();
    let radiation = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap();
    radiation.ambient_temperature = Some(QtyAny::new(300.0, dims::TEMPERATURE));
    radiation.view_factors = vec![vec![0.0, 0.25], vec![0.25, 0.0]];
    f
}

#[test]
fn source_curve_and_surface_contact_share_closed_enclosure_and_original_air_heat() {
    let f = fixture();
    let curve = receipt(&f, 190);
    let surface_pack = surface_pack();
    let selected = surface_pack.2.clone();
    let surface = receipt_with_pack(&f, 191, surface_pack);
    let field = |r: &J| artifact(&f.dir.join("temperature.db"),
        r.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone();
    for (a, b) in field(&curve).as_array().unwrap().iter().zip(field(&surface).as_array().unwrap()) {
        assert!((a.as_f64().unwrap() - b.as_f64().unwrap()).abs() < 2e-7);
    }
    for (result, picard) in [(&curve, false), (&surface, true)] {
        let contact = result.get("temperature_contact").unwrap();
        let composition = contact.get("enclosure_radiation").unwrap();
        let radiation = result.get("radiation").unwrap();
        let air = result.get("conjugate").unwrap();
        assert_eq!(radiation.str_field("model"), Some("closed-gray-diffuse-enclosure"));
        assert_eq!(contact.str_field("algorithm"), Some(if picard {
            "bounded-contact-picard" } else { "coupled-newton" }));
        assert_eq!(composition.get("physical_endpoint_reassembled"), Some(&J::Bool(true)));
        assert_eq!(number(composition, "thermal_responses"), number(radiation, "solid_solves"));
        assert_eq!(number(composition, "radiation_responses"), number(radiation, "surface_iterations") + 1.0);
        if picard {
            assert!(number(composition, "thermal_responses") > number(composition, "radiation_responses"));
        } else {
            assert_eq!(number(composition, "thermal_responses"), number(composition, "radiation_responses"));
        }
        assert!(number(contact, "total_thermal_responses") >= number(composition, "thermal_responses"));
        assert!(number(contact, "total_nonlinear_iterations") >= number(radiation, "solid_iterations"));
        assert!(number(contact, "total_krylov_iterations") >= number(radiation, "krylov_iterations"));
        assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
        assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
        assert_eq!(number(result.get("energy").unwrap(), "closure_w"), number(contact, "energy_closure_w"));
        for state in contact.get("face_states").unwrap().as_array().unwrap() {
            let mean = 0.5 * number(state, "side_a_temperature_k") + 0.5 * number(state, "side_b_temperature_k");
            assert_eq!(number(state, "query_temperature_k"), mean);
            assert!((number(state, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
            if picard { assert_eq!(state.str_field("selected_claim"), Some(selected.as_str())); }
        }
        let patches = radiation.get("surfaces").unwrap().as_array().unwrap();
        assert!(patches.iter().any(|p| number(p, "nonlinear_heat_w") > 0.1));
        assert!(patches.iter().any(|p| number(p, "nonlinear_heat_w") < -0.1));
        assert!(number(radiation, "nonlinear_radiative_out_w").abs() < 1e-7);
        let convection = number(radiation, "convective_out_w");
        assert!((convection - 4.44).abs() < 1e-7, "closed exchange remains internal");
        assert!((number(air, "solid_total_w") - convection).abs() < 1e-7);
        assert!((number(air, "air_total_w") - convection).abs() <= number(air, "balance_tolerance_w"));
        assert!(number(radiation, "decomposition_residual_w").abs() < 1e-7);
    }
}

#[test]
fn source_curve_and_surface_contact_share_open_cavity_without_heating_air_with_radiation() {
    let f = open_fixture();
    let curve = receipt(&f, 192);
    let surface_pack = surface_pack();
    let selected = surface_pack.2.clone();
    let surface = receipt_with_pack(&f, 193, surface_pack);
    let field = |r: &J| artifact(&f.dir.join("temperature.db"),
        r.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone();
    for (a, b) in field(&curve).as_array().unwrap().iter().zip(field(&surface).as_array().unwrap()) {
        assert!((a.as_f64().unwrap() - b.as_f64().unwrap()).abs() < 2e-7);
    }
    for (result, picard) in [(&curve, false), (&surface, true)] {
        let contact = result.get("temperature_contact").unwrap();
        let composition = contact.get("enclosure_radiation").unwrap();
        let radiation = result.get("radiation").unwrap();
        let air = result.get("conjugate").unwrap();
        assert_eq!(radiation.str_field("model"), Some("open-gray-diffuse-enclosure"));
        assert_eq!(contact.str_field("algorithm"), Some(if picard {
            "bounded-contact-picard" } else { "coupled-newton" }));
        assert_eq!(composition.get("physical_endpoint_reassembled"), Some(&J::Bool(true)));
        assert_eq!(number(composition, "thermal_responses"), number(radiation, "solid_solves"));
        assert_eq!(number(composition, "radiation_responses"), number(radiation, "surface_iterations") + 1.0);
        assert_eq!(number(composition, "radiation_responses"), number(radiation, "radiation_responses"));
        if picard {
            assert!(number(composition, "thermal_responses") > number(composition, "radiation_responses"));
        } else {
            assert_eq!(number(composition, "thermal_responses"), number(composition, "radiation_responses"));
        }
        assert!(number(contact, "total_thermal_responses") >= number(composition, "thermal_responses"));
        assert!(number(contact, "total_nonlinear_iterations") >= number(radiation, "solid_iterations"));
        assert!(number(contact, "total_krylov_iterations") >= number(radiation, "krylov_iterations"));
        assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
        assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
        let energy = result.get("energy").unwrap();
        assert_eq!(number(energy, "closure_w"), number(contact, "energy_closure_w"));
        for state in contact.get("face_states").unwrap().as_array().unwrap() {
            let mean = 0.5 * number(state, "side_a_temperature_k") + 0.5 * number(state, "side_b_temperature_k");
            assert_eq!(number(state, "query_temperature_k"), mean);
            assert!((number(state, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
            if picard { assert_eq!(state.str_field("selected_claim"), Some(selected.as_str())); }
        }
        let patches = radiation.get("surfaces").unwrap().as_array().unwrap();
        assert_eq!(patches.len(), 2, "the black bath must not become a third solid surface");
        let physical_outward: f64 = patches.iter().map(|p| number(p, "nonlinear_heat_w")).sum();
        assert!(physical_outward.abs() > 0.1);
        assert!((physical_outward - number(radiation, "ambient_absorbed_w")).abs() < 1e-7);
        assert!((physical_outward - number(radiation, "nonlinear_radiative_out_w")).abs() < 1e-7);
        assert!(number(radiation, "enclosure_energy_closure_w").abs() < 1e-7);
        let convection = number(radiation, "convective_out_w");
        assert!((convection + physical_outward - 4.44).abs() < 1e-7);
        assert!((number(energy, "robin_out_w") - convection - physical_outward).abs() < 1e-7);
        assert!((number(air, "solid_total_w") - convection).abs() < 1e-7);
        assert!((number(air, "air_total_w") - convection).abs() <= number(air, "balance_tolerance_w"));
        assert!(number(radiation, "decomposition_residual_w").abs() < 1e-7);
    }
}

fn run_output(f: &Fixture, ordinal: usize, output: Option<&str>,
    pack: (Vec<u8>, String, String)) -> (fs_cli::CommandOutput, PathBuf) {
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, output, pack);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = fs_cli::run(vec!["--json".into(), "solve".into(),
        source.display().to_string(), ledger.display().to_string(), "--materials".into(),
        f.dir.join("solid.fsmcdpk").display().to_string(), "--materials".into(),
        root.join("data/reference-project/gray-surface.fsmcdpk").display().to_string(),
        "--interfaces".into(), contact.display().to_string()]);
    (output, ledger)
}

fn solve_output(f: &Fixture, ordinal: usize, output: Option<&str>,
    pack: (Vec<u8>, String, String)) -> (J, J) {
    let (output, ledger) = run_output(f, ordinal, output, pack);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", output.stdout, output.stderr);
    let result = J::parse(&output.stdout).unwrap();
    let run = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = run.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn selected_temperature(field: &J, vertex: usize) -> f64 {
    field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap()
}

fn control_derivative(report: &J, target: &str, entity: &str) -> f64 {
    let rows: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field("target") == Some(target)
            && row.str_field("entity") == Some(entity)).collect();
    assert_eq!(rows.len(), 1, "one original physical control: {target}/{entity}");
    number(rows[0], "derivative")
}

#[test]
fn complete_live_contact_enclosure_adjoints_match_independent_physical_resolves() {
    use fs_project::ThermalBoundaryCondition as B;
    // Two complementary compositions, not a product of unrelated solver modes:
    // the closed curve uses Newton; the open pressure/T source uses Picard.
    for open in [false, true] {
        let mut f = if open { open_fixture() } else { fixture() };
        // Keep the gas/hydraulic operating point fixed while perturbing the
        // independently authored inlet inside the same admitted envelope.
        let envelope = f.project.envelope.as_mut().unwrap();
        envelope.ambient_lo.value = 296.0;
        envelope.ambient_hi.value = 298.0;
        let pack = |scale| if open { scaled_surface_pack(scale) } else {
            temperature_pack_from_knots(vec![(200.0, 0.2 * scale), (450.0, 0.45 * scale)])
        };
        let (plain, plain_field) = solve_output(&f, 200, None, pack(1.0));
        let (nominal, field) = solve_output(&f, 201, Some(OUTPUT), pack(1.0));
        assert_eq!(plain_field.get("temperature"), field.get("temperature"));
        for key in ["energy", "temperature_contact", "radiation", "conjugate"] {
            assert_eq!(plain.get(key), nominal.get(key), "adjoint changed the primal {key}");
        }
        let report = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("mode"), Some(if open {
            "nonlinear-contact-open-enclosure-full-radiosity-feedback"
        } else { "nonlinear-contact-enclosure-full-radiosity-feedback" }));
        let selected = number(report, "selected_vertex") as usize;
        assert_eq!(number(report, "value_k"), selected_temperature(&field, selected));
        for (case, (target, entity, step)) in [
            (TARGET, "cold-hot-joint", 0.002), ("power", "hot", 0.01),
            ("radiation-reservoir-temperature", "enclosure-environment", 0.01),
            ("air-inlet-temperature", "cold", 0.01),
        ].into_iter().enumerate() {
            if !open && case >= 2 { continue; }
            let derivative = control_derivative(report, target, entity);
            assert!(derivative.abs() > 1e-6, "{target} must exercise the complete dual");
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
                let scale = if case == 0 { 1.0 + sign * step } else { 1.0 };
                if case == 1 {
                    changed.project.power.as_mut().unwrap()[0].watts.value += sign * step;
                } else if case == 2 {
                    changed.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .enclosure_radiation.as_mut().unwrap().ambient_temperature.as_mut().unwrap()
                        .value += sign * step;
                } else if case == 3 {
                    let B::AirflowConvection { inlet_temperature, .. } = &mut changed.project.cooling
                        .as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
                        else { panic!("original transported-air inlet"); };
                    inlet_temperature.value += sign * step;
                }
                let (_, perturbed) = solve_output(&changed, 202 + 2 * case + side, None, pack(scale));
                // Differentiate the SAME nominal selected node; a max-node
                // switch in a perturbed solve must not change the objective.
                values.push(selected_temperature(&perturbed, selected));
            }
            let finite_difference = (values[1] - values[0]) / (2.0 * step);
            assert!((derivative - finite_difference).abs()
                < 5e-5 * finite_difference.abs().max(1e-3),
                "open={open} {target}: complete adjoint {derivative:e} vs resolve {finite_difference:e}");
        }
        if open {
            let (ordinary, ordinary_field) = solve_output(&f, 210, Some("temperature-max-adjoint"), pack(1.0));
            assert_eq!(ordinary_field.get("temperature"), field.get("temperature"));
            let ordinary = ordinary.get("nominal_adjoint").unwrap();
            assert_eq!(ordinary.get("selected_vertex"), report.get("selected_vertex"));
            assert_eq!(control_derivative(ordinary, "power", "hot"),
                control_derivative(report, "power", "hot"));
        }
    }
}

#[test]
fn live_contact_enclosure_prescribed_temperature_controls_remain_unsupported() {
    for f in [fixture(), open_fixture()] {
        let (output, _) = run_output(&f, 211,
            Some("temperature-max-contact-boundary-adjoint"), temperature_pack());
        assert_ne!(output.exit_code, fs_cli::exit::SUCCESS);
        assert!(format!("{}{}", output.stdout, output.stderr)
            .contains("cli-solve-conduction-contact-composition"));
    }
}
