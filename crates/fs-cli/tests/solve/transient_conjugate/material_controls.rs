//! Capacity-storage trajectory material controls through the native pipeline.
//! Fresh probes rebuild immutable conductivity or contact cards, request no
//! adjoint, and never reuse the nominal thermal response.

use super::*;

const FINAL: &str = "temperature-final-material-adjoint";
const PEAK: &str = "temperature-sampled-peak-material-adjoint";
const DK: &str = "dtemperature_dconductivity_multiplier_k";
const DR: &str = "dtemperature_dresistance_multiplier_k";

fn rows<'a>(value: &'a J, name: &str) -> &'a [J] {
    value.get(name).unwrap().as_array().unwrap()
}

fn controls(solved: &ThermalRun) -> &J {
    let adjoint = solved.transient().get("adjoint").unwrap();
    assert_eq!(adjoint.str_field("schema"), Some("fs-cli-native-transient-adjoint-v1"));
    assert!([FINAL, PEAK].contains(&adjoint.str_field("output").unwrap()));
    let scope = adjoint.str_field("scope").unwrap();
    assert!(scope.contains("Independent regional conductivity multipliers"));
    assert!(!scope.contains("prescribed-temperature, conductivity, contact, geometry"));
    let controls = adjoint.get("material_controls").unwrap();
    close(n(controls, "control_steps"), n(adjoint, "state_index"), 0.0);
    assert_eq!(controls.str_field("units"), Some("K"));
    controls
}

fn assert_bindings(controls: &J, spec: &ProjectSpec) {
    let actual = rows(controls, "conductivity_regions");
    let bindings = spec.materials.as_ref().unwrap();
    assert_eq!(actual.len(), bindings.len());
    for (row, binding) in actual.iter().zip(bindings) {
        assert_eq!(row.str_field("region"), Some(binding.region.as_str()));
        assert_eq!(row.str_field("material_card"), Some(binding.card.as_str()));
        assert_eq!(row.str_field("material_state"), Some(binding.state.as_str()));
        assert!(n(row, "elements") > 0.0);
        close(n(row, "reference_value"), 1.0, 0.0);
        close(n(row, "log_reference_value"), 0.0, 0.0);
        assert!(n(row, DK).is_finite());
    }
}

fn mixed_fixture(scale: f64) -> (ProjectSpec, CardPackSet, Vec<u8>) {
    let base = CardPackSet::admit(vec![raw_pack(
        CardPackKind::Material,
        "fixtures/capacity-conductivity-control.fsmcdpk",
        material_pack_bytes_with_property(
            "AA6061", &format!("capacity-conductivity-control-{scale}"), 280.0, 360.0,
            fs_matdb::PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
                // Scale all ordinates, not temperature, source power or C.
                // There is a nonzero smooth K-prime throughout the trajectory.
                knots: vec![(280.0, 1.6 * scale), (360.0, 3.2 * scale)],
                dims: CONDUCTIVITY_DIMS,
            },
            fs_matdb::InterpolationPolicy::LinearInside,
        ),
    )]).unwrap();
    let (cards, radiation_card) = radiation_product::emissivity_cards(&base, EMISSIVITY);
    let (mut spec, bytes) = mixed_adjoint_project(&radiation_card);
    rebind_to(&mut spec, &cards);
    (spec, cards, bytes)
}

#[test]
fn g3_native_capacity_material_nonlinear_mixed_derivative_matches_fresh_cards() {
    let (spec, cards, bytes) = mixed_fixture(1.0);
    let plain = run(&spec, &cards, &bytes);
    let mut requested = spec.clone();
    request_air_adjoint(&mut requested, FINAL);
    let solved = run(&requested, &cards, &bytes);
    assert_eq!(plain.temperature, solved.temperature);
    assert_eq!(plain.rows("coarse"), solved.rows("coarse"));
    assert_eq!(plain.rows("fine"), solved.rows("fine"));
    let result = controls(&solved);
    assert_bindings(result, &spec);
    assert!(rows(result, "contact_interfaces").is_empty());
    let material = &rows(result, "conductivity_regions")[0];
    let derivative = n(material, DK);
    assert!(derivative.abs() > 1e-6, "the nonuniform fixture must observe conductivity");
    let adjoint = solved.transient().get("adjoint").unwrap();
    close(n(adjoint, "value_k"), final_k(&plain), 0.0);
    close(n(adjoint, "reconstructed_solid_endpoints"), 0.0, 0.0);
    close(n(adjoint, "revalidated_solid_endpoints"), solved.rows("fine").len() as f64, 0.0);
    assert!(air_gradients(adjoint, true).iter().all(|value| value.is_finite()));
    assert!(n(air_control(adjoint, "natural_reservoirs", "target", "natural-wall"),
        "dtemperature_dambient_temperature_k_per_k").is_finite());

    // The explicit material request adds contractions, not a different state
    // or another transpose. Existing report names retain their old scope.
    let mut legacy_spec = spec.clone();
    request_air_adjoint(&mut legacy_spec, "temperature-final-adjoint");
    let legacy = run(&legacy_spec, &cards, &bytes);
    let legacy_adjoint = legacy.transient().get("adjoint").unwrap();
    assert!(legacy_adjoint.get("material_controls").is_none());
    assert_eq!(legacy.temperature, solved.temperature);
    for key in [
        "value_k", "time_s", "state_index", "active_vertex", "active_ties",
        "dtemperature_duniform_initial_k_per_k", "dtemperature_dinitial_temperatures",
        "dtemperature_dcapacity_multiplier_k", "regional_power", "intervals", "robin",
        "convection", "radiation", "adjoint_krylov_iterations",
        "maximum_dual_relative_residual", "revalidated_solid_endpoints",
    ] {
        assert_eq!(legacy_adjoint.get(key), adjoint.get(key), "{key}");
    }
    assert!(n(adjoint, "additional_storage_estimate_bytes")
        > n(legacy_adjoint, "additional_storage_estimate_bytes"));

    let evaluate = |scale: f64| {
        let (probe, probe_cards, probe_bytes) = mixed_fixture(scale);
        assert_eq!(bytes, probe_bytes);
        assert_eq!(spec.envelope, probe.envelope);
        assert_ne!(spec.materials.as_ref().unwrap()[0].card,
            probe.materials.as_ref().unwrap()[0].card);
        final_k(&run(&probe, &probe_cards, &probe_bytes))
    };
    // Nonunit fine timesteps expose an accidental second dt multiplication.
    assert!(solved.rows("fine").iter().all(|row| n(row, "dt_s") < 0.5));
    for epsilon in [0.002, 0.001] {
        let finite_difference = (evaluate(1.0 + epsilon) - evaluate(1.0 - epsilon))
            / (2.0 * epsilon);
        close(derivative, finite_difference, 5e-5 * derivative.abs().max(1.0));
    }
    let repeated = run(&requested, &cards, &bytes);
    assert_eq!(result, controls(&repeated));
}

fn contact_fixture(scale: f64) -> (ProjectSpec, CardPackSet) {
    let cards = contact_cards_with_interface(
        interface_pack_bytes_with_resistance_and_peak(0.1 * scale, None),
    );
    let mut spec = multi_region_contact_project();
    spec.interface_cards.as_mut().unwrap()[0].card = cards.interfaces()[0].card().to_hex();
    transient_product::declaration(&mut spec, 100.0, 2.0);
    spec.requirements.as_mut().unwrap()[0].region = "hot".into();
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
    }
    (spec, cards)
}

fn run_contact(spec: &ProjectSpec, cards: &CardPackSet, output: Option<&str>) -> ThermalRun {
    let mut spec = spec.clone();
    if let Some(output) = output { request_air_adjoint(&mut spec, output); }
    let decoded = decode(&spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&ledger, &spec);
    let (_, qoi, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    let field = J::parse(&field).unwrap();
    ThermalRun {
        qoi: J::parse(&qoi).unwrap(),
        receipt: J::parse(&receipt).unwrap(),
        temperature: rows(&field, "temperature").iter().map(|value| value.as_f64().unwrap()).collect(),
    }
}

fn contact_row<'a>(solved: &'a ThermalRun, spec: &ProjectSpec) -> &'a J {
    let result = controls(solved);
    assert_bindings(result, spec);
    let material_rows = rows(result, "conductivity_regions");
    assert_eq!(material_rows.len(), 2);
    assert_eq!(material_rows[0].str_field("material_card"), material_rows[1].str_field("material_card"));
    assert_ne!(material_rows[0].str_field("region"), material_rows[1].str_field("region"));
    let contacts = rows(result, "contact_interfaces");
    assert_eq!(contacts.len(), 1);
    let contact = &contacts[0];
    let binding = &spec.interface_cards.as_ref().unwrap()[0];
    assert_eq!(contact.str_field("interface"), Some(binding.interface.as_str()));
    assert_eq!(contact.str_field("interface_card"), Some(binding.card.as_str()));
    assert_eq!(contact.get("mapped"), Some(&J::Bool(false)));
    assert!(n(contact, "face_pairs") > 0.0);
    close(n(contact, "reference_value"), 1.0, 0.0);
    close(n(contact, "log_reference_value"), 0.0, 0.0);
    contact
}

#[test]
fn g3_native_capacity_material_contact_derivative_matches_fresh_interface_cards() {
    let (spec, cards) = contact_fixture(1.0);
    let plain = run_contact(&spec, &cards, None);
    let solved = run_contact(&spec, &cards, Some(FINAL));
    assert_eq!(plain.temperature, solved.temperature);
    assert_eq!(plain.rows("fine"), solved.rows("fine"));
    let derivative = n(contact_row(&solved, &spec), DR);
    assert!(derivative > 1e-6, "joint resistance must retain heat in the powered body");
    let adjoint = solved.transient().get("adjoint").unwrap();
    close(n(adjoint, "value_k"), final_k(&plain), 0.0);
    close(n(adjoint, "reconstructed_solid_endpoints"), solved.rows("fine").len() as f64, 0.0);
    let evaluate = |scale: f64| {
        let (probe, probe_cards) = contact_fixture(scale);
        assert_ne!(spec.interface_cards.as_ref().unwrap()[0].card,
            probe.interface_cards.as_ref().unwrap()[0].card);
        final_k(&run_contact(&probe, &probe_cards, None))
    };
    for epsilon in [0.002, 0.001] {
        close(derivative, (evaluate(1.0 + epsilon) - evaluate(1.0 - epsilon)) / (2.0 * epsilon),
            5e-5 * derivative.abs().max(1.0));
    }
}

#[test]
fn g1_native_capacity_material_initial_peak_keeps_zero_material_and_contact_controls() {
    let (mut spec, cards) = contact_fixture(1.0);
    for power in spec.power.iter_mut().flatten() { power.watts.value = 0.0; }
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient.as_mut().unwrap().initial_temperature =
        fs_project::TransientInitialTemperature::Uniform(QtyAny::new(330.0, dims::TEMPERATURE));
    for boundary in &mut setup.boundaries {
        // Both bodies cool from the same initial temperature. Future contact
        // and conduction work cannot influence a selected initial peak.
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.5, dims::HEAT_FLUX),
        };
    }
    let solved = run_contact(&spec, &cards, Some(PEAK));
    assert!(final_k(&solved) < 330.0);
    let adjoint = solved.transient().get("adjoint").unwrap();
    close(n(adjoint, "state_index"), 0.0, 0.0);
    close(n(adjoint, "value_k"), 330.0, 0.0);
    close(n(adjoint, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(n(adjoint, "reconstructed_solid_endpoints"), 0.0, 0.0);
    close(n(contact_row(&solved, &spec), DR), 0.0, 0.0);
    for material in rows(controls(&solved), "conductivity_regions") {
        close(n(material, DK), 0.0, 0.0);
    }
}

#[test]
fn g0_native_capacity_material_report_requests_reject_ambiguous_selection() {
    for invalid in ["duplicate", "scalar", "region"] {
        let (mut spec, bytes) = regular_project();
        request_air_adjoint(&mut spec, FINAL);
        if invalid == "duplicate" {
            request_air_adjoint(&mut spec, "temperature-final-adjoint");
        } else {
            let output = spec.outputs.as_mut().unwrap().last_mut().unwrap();
            if invalid == "scalar" { output.kind = "scalar".into(); }
            else { output.region = Some("air".into()); }
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let result = run_solve(
            &ledger, &CancelGate::new_clock_free(), &mut benign_clock(),
            &decode(&spec), &fixture_cards(), &mut Vec::new(),
        );
        assert!(result.is_err(), "{invalid} must not silently select an adjoint");
        assert_eq!(artifacts_of_kind(&ledger, "solve-conduction-solution"), 0);
    }
}
