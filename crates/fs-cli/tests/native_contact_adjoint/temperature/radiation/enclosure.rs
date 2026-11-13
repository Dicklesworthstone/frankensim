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
        emissivity_response: None,
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
    let (source, contact, ledger) = prepare_fixed_pressure_surface(f, ordinal, output, pack);
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

fn adaptive_output(f: &Fixture, ordinal: usize) -> (J, J, J) {
    let (output, ledger) = run_output(f, ordinal, None, temperature_pack());
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", output.stdout, output.stderr);
    let output = J::parse(&output.stdout).unwrap();
    let run = artifact(&ledger, output.str_field("run_receipt").unwrap());
    let stages = run.get("stages").unwrap().as_array().unwrap();
    let stage = |name| {
        let row = stages.iter().find(|row| row.str_field("stage") == Some(name)).unwrap();
        artifact(&ledger, row.str_field("receipt").unwrap())
    };
    let receipt = stage("conduction");
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field, stage("qoi"))
}

fn independent_enclosure_heat(radiation: &J, open: bool) -> Vec<f64> {
    const SIGMA: f64 = 5.670_374_419e-8;
    let patches = radiation.get("surfaces").unwrap().as_array().unwrap();
    let epsilon = number(&patches[0], "emissivity");
    assert!(epsilon > 0.0 && epsilon < 1.0, "reflection must remain active");
    for patch in patches {
        assert!((number(patch, "area_m2") - 5.0).abs() < 1e-12);
        assert_eq!(number(patch, "emissivity"), epsilon);
    }
    let t: Vec<_> = patches.iter().map(|patch| number(patch, "mean_temperature_k")).collect();
    if !open {
        assert_eq!(t.len(), 2);
        let heat = 5.0 * SIGMA * (t[0] - t[1]) * (t[0] + t[1])
            * (t[0]*t[0] + t[1]*t[1]) / (2.0 / epsilon - 1.0);
        return vec![heat, -heat];
    }
    let bath = number(radiation, "ambient_temperature_k");
    assert_eq!(bath, 300.0);
    // Solve the authored one- or two-patch radiosity equations directly,
    // relative to the black bath, without invoking the production solver.
    let emission: Vec<_> = t.iter().map(|&t|
        SIGMA * (t - bath) * (t + bath) * (t*t + bath*bath)).collect();
    let view = 0.25;
    let reflected = (1.0 - epsilon) * view;
    if patches.len() == 1 {
        return vec![5.0 * (1.0 - view) * epsilon * emission[0] / (1.0 - reflected)];
    }
    assert_eq!(patches.len(), 2);
    let delta_j = [
        epsilon * (emission[0] + reflected * emission[1]) / (1.0 - reflected*reflected),
        epsilon * (emission[1] + reflected * emission[0]) / (1.0 - reflected*reflected),
    ];
    vec![5.0 * (delta_j[0] - view*delta_j[1]), 5.0 * (delta_j[1] - view*delta_j[0])]
}

#[test]
fn live_contact_adaptive_enclosures_preserve_reflection_air_heat_and_physical_gates() {
    use fs_project::ThermalBoundaryCondition as B;
    for mode in ["closed-air", "open-air", "open-vacuum"] {
        let open = mode != "closed-air";
        let vacuum = mode == "open-vacuum";
        let mut f = if open { open_fixture() } else { fixture() };
        let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        if vacuum {
            setup.boundaries[0].condition = B::HeatFlux {
                outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
            };
            setup.boundaries[1].condition = B::FixedTemperature {
                temperature: QtyAny::new(310.0, dims::TEMPERATURE),
            };
        }
        let declaration = setup.enclosure_radiation.as_mut().unwrap();
        declaration.temperature_tolerance.value = 1e-11;
        declaration.heat_tolerance.value = 1e-10;
        if vacuum {
            // No selected Robin rows: reflection must use the original cold
            // flux trace while the separate hot reservoir anchors conduction.
            declaration.surfaces.retain(|surface| surface.target == "cold");
            declaration.view_factors = vec![vec![0.25]];
            f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
        }
        f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
        f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
        f.project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
        let (receipt, field, qoi) = adaptive_output(&f, 280);
        let adaptive = receipt.get("adaptive").unwrap();
        assert_eq!(adaptive.str_field("status"), Some("observed-tolerance-met"), "{mode}");
        assert!(number(adaptive, "solved_meshes") > 1.0);
        assert_eq!(adaptive.get("continuum_error_bound"), Some(&J::Bool(false)));
        assert_eq!(receipt.get("solver_control"), Some(&J::Null));
        let history = adaptive.get("history").unwrap().as_array().unwrap();
        assert!(!history.is_empty());
        for row in history {
            let balance = row.get("algebraic_balance").unwrap();
            assert_eq!(balance.str_field("comparison_model"), Some("source-contact-enclosure-full-transpose"));
            assert_eq!(balance.str_field("status"), Some("unsupported-model"));
            for key in ["contact_temperature_feedback", "radiation_feedback", "enclosure_feedback"] {
                assert_eq!(balance.get(key), Some(&J::Bool(true)), "{mode}: {key}");
            }
            assert_eq!(balance.get("open_enclosure"), Some(&J::Bool(open)));
            assert_eq!(balance.get("air_reference_feedback"), Some(&J::Bool(!vacuum)));
            assert_eq!(row.get("uses_nonlinear_jacobian"), Some(&J::Bool(true)));
            assert_eq!(row.get("natural_coefficient_feedback"), Some(&J::Bool(false)));
            assert_eq!(row.get("coupled_adjoint"), Some(&J::Null));
            assert_eq!(number(row, "probe_tets"), 8.0 * number(row, "tets"));
            assert!(number(row, "interface_pairs") > 0.0);
            assert_eq!(number(row, "probe_interface_pairs"), 4.0 * number(row, "interface_pairs"));
            assert!(number(row, "absolute_contribution_sum_k") > 0.0);
            let measured = number(row, "measured_change_k");
            assert!((number(row, "signed_linear_change_k") + number(row, "linearization_remainder_k")
                + number(row, "maximum_remainder_k") - measured).abs() < 1e-9);
            assert!(number(row, "estimated_change_k") + 1e-12 >= measured.abs());
            assert!(number(row, "dual_residual") <= 1e-9);
            if vacuum {
                // The 300 K bath, below the 310 K fixed anchor, defines the
                // physical rise even though radiation has no Robin owner.
                let tolerance = 0.2 * (number(row, "t_max_k") - 300.0).abs();
                assert_eq!(number(row, "tolerance_k").to_bits(), tolerance.to_bits());
            }
        }
        let last = history.last().unwrap();
        let goal = &qoi.get("qoi").unwrap().as_array().unwrap()[0];
        assert_eq!(goal.str_field("region"), Some(if vacuum { "cold" } else { "hot" }));
        assert_eq!(number(goal, "value"), selected_temperature(&field, number(goal, "witness_vertex") as usize));
        assert!((number(goal, "value") - number(last, "t_max_k")).abs() < 1e-9);
        assert_eq!(goal.str_field("color"), Some("estimated"));
        let contact = receipt.get("temperature_contact").unwrap();
        assert_eq!(contact.str_field("algorithm"), Some("coupled-newton"));
        assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
        assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
        assert_eq!(contact.get("enclosure_radiation").unwrap().get("physical_endpoint_reassembled"),
            Some(&J::Bool(true)));
        let (_, _, selected) = temperature_pack();
        let faces = contact.get("face_states").unwrap().as_array().unwrap();
        assert!(!faces.is_empty());
        for face in faces {
            let mean = 0.5 * number(face, "side_a_temperature_k") + 0.5 * number(face, "side_b_temperature_k");
            assert_eq!(number(face, "query_temperature_k"), mean);
            assert!((number(face, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
            assert_eq!(face.str_field("selected_claim"), Some(selected.as_str()));
            assert_eq!(face.str_field("usage_receipt").unwrap().len(), 64);
        }
        let radiation = receipt.get("radiation").unwrap();
        assert_eq!(radiation.str_field("model"), Some(if open {
            "open-gray-diffuse-enclosure" } else { "closed-gray-diffuse-enclosure" }));
        let patches = radiation.get("surfaces").unwrap().as_array().unwrap();
        assert_eq!(patches.len(), if vacuum { 1 } else { 2 }, "the bath has no solid trace");
        let expected = independent_enclosure_heat(radiation, open);
        assert!(expected.iter().any(|heat| heat.abs() > 0.1), "{mode}: nontrivial reflected exchange");
        for (patch, expected) in patches.iter().zip(expected) {
            assert!((number(patch, "nonlinear_heat_w") - expected).abs() < 1e-7);
        }
        let radiative: f64 = patches.iter().map(|patch| number(patch, "nonlinear_heat_w")).sum();
        assert!((radiative - number(radiation, "nonlinear_radiative_out_w")).abs() < 1e-7);
        if open {
            assert!(radiative.abs() > 0.1);
            assert!((radiative - number(radiation, "ambient_absorbed_w")).abs() < 1e-7);
            assert!(number(radiation, "enclosure_energy_closure_w").abs() < 1e-7);
        } else {
            assert!(radiative.abs() < 1e-7, "closed reflected exchange is internal");
        }
        let convection = number(radiation, "convective_out_w");
        let energy = receipt.get("energy").unwrap();
        assert!((number(energy, "source_w") - 4.44).abs() < 1e-9);
        assert_eq!(number(energy, "closure_w"), number(contact, "energy_closure_w"));
        assert!((number(energy, "source_w") + number(energy, "dirichlet_in_w")
            - convection - radiative).abs() <= number(contact, "energy_tolerance_w"));
        assert!(number(radiation, "decomposition_residual_w").abs() < 1e-7);
        assert_eq!(receipt.get("conjugate").unwrap() != &J::Null, !vacuum);
        if vacuum {
            assert_eq!(convection, 0.0, "radiation-only traces must not invent convection");
        } else {
            let air = receipt.get("conjugate").unwrap();
            assert!((number(air, "solid_total_w") - convection).abs() < 1e-9);
            assert!((number(air, "air_total_w") - convection).abs() <= number(contact, "residual_tolerance_w"));
            if open {
                assert!((number(air, "air_total_w") - convection - radiative).abs() > 0.1,
                    "the independent black bath must not heat transported air");
            }
        }
        let term = qoi.get("budget").unwrap().as_array().unwrap()[0].get("terms").unwrap()
            .as_array().unwrap().iter().find(|term| term.str_field("kind") == Some("discretization")).unwrap();
        assert_eq!(qoi.str_field("authority"), Some("estimated-candidate"));
        assert_eq!(term.str_field("state"), Some("interval"));
        assert_eq!(term.get("derivation").unwrap().str_field("method"), Some("adaptive-enrichment-order-one"));
        let half_width = 2.0 * number(last, "estimated_change_k").abs().max(number(last, "measured_change_k").abs());
        assert_eq!(number(term, "upper_kelvin").to_bits(), half_width.to_bits());

        if mode == "closed-air" {
            let (replay, replay_field, _) = adaptive_output(&f, 281);
            assert_eq!(receipt.get("adaptive"), replay.get("adaptive"));
            assert_eq!(field.get("temperature"), replay_field.get("temperature"));
            // One independent first enrichment, bounded to two uniform rungs;
            // this checks contact and radiosity rebinding on the identical mesh.
            let first = &history[0];
            f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
            f.project.budgets.as_mut().unwrap().memory_bytes = (number(first, "probe_tets") as u64)
                .checked_mul(256).unwrap();
            let (uniform, _, _) = adaptive_output(&f, 282);
            let rungs = uniform.path(&["ladder", "rungs"]).unwrap().as_array().unwrap();
            assert_eq!(rungs.len(), 2);
            assert_eq!(number(first, "probe_tets"), number(&rungs[1], "tets"));
            assert!((number(first, "probe_t_max_k") - number(&rungs[1], "t_max_k")).abs() < 1e-6);
        }
    }
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
fn live_contact_fixed_temperature_controls_include_ambient_and_radiosity_feedback() {
    use fs_project::ThermalBoundaryCondition as B;
    for enclosure in [false, true] {
        let mut f = if enclosure { open_fixture() } else { super::fixture(true) };
        let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.boundaries[0].condition = B::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
        setup.boundaries[1].condition = B::FixedTemperature {
            temperature: QtyAny::new(310.0, dims::TEMPERATURE),
        };
        if enclosure {
            // Keep the prescribed hot patch separate from radiation. The cold
            // concave patch reflects to itself and escapes to the black bath;
            // its original zero-flux row exercises the mean-flux lift as well.
            let radiation = setup.enclosure_radiation.as_mut().unwrap();
            radiation.surfaces.retain(|surface| surface.target == "cold");
            radiation.view_factors = vec![vec![0.25]];
            radiation.temperature_tolerance.value = 1e-11;
            radiation.heat_tolerance.value = 1e-10;
        }
        f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
        let (plain, original) = solve_output(&f, 250, None, surface_pack());
        let (combined, field) = solve_output(&f, 251,
            Some("temperature-max-contact-boundary-adjoint"), surface_pack());
        let (boundary, boundary_field) = solve_output(&f, 252,
            Some("temperature-max-boundary-adjoint"), surface_pack());
        assert_eq!(original.get("temperature"), field.get("temperature"));
        assert_eq!(boundary_field.get("temperature"), field.get("temperature"));
        for receipt in [&combined, &boundary] {
            for key in ["energy", "interfaces", "temperature_contact", "radiation", "conjugate"] {
                assert_eq!(
                    forward_section(&plain, key),
                    forward_section(receipt, key),
                    "adjoint changed primal {key}"
                );
            }
        }
        let report = combined.get("nominal_adjoint").unwrap();
        let boundary = boundary.get("nominal_adjoint").unwrap();
        let shared: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target") != Some(TARGET)).collect();
        assert_eq!(shared, boundary.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
        assert_eq!(report.get("selected_vertex"), boundary.get("selected_vertex"));
        let selected = number(report, "selected_vertex") as usize;
        let actual = control_derivative(report, "fixed-temperature", "hot");
        assert!(actual > 1e-6);
        assert!(control_derivative(report, TARGET, "cold-hot-joint").abs() > 1e-6);
        let step = 0.01;
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
            let B::FixedTemperature { temperature } = &mut changed.project.cooling.as_mut().unwrap()
                .conduction.as_mut().unwrap().boundaries[1].condition else { panic!("original hot reservoir"); };
            temperature.value += sign * step;
            let (_, field) = solve_output(&changed, 253 + side, None, surface_pack());
            values.push(selected_temperature(&field, selected));
        }
        let expected = (values[1] - values[0]) / (2.0 * step);
        assert!((actual - expected).abs() < 5e-5 * expected.abs().max(1e-3),
            "enclosure={enclosure}: complete prescribed lift {actual:e} vs resolve {expected:e}");
    }
}
