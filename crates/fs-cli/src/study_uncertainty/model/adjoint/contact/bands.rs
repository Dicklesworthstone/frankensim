//! Manufactured bands reach actual native contact operators and thermal budgets.
use super::*;

fn receipt(ledger: &Ledger, sample: &Sample) -> J {
    let sealed = crate::solve::load_completed_run(ledger, &sample.run).unwrap();
    let hash = ContentHash::from_hex(&sealed.stages.iter()
        .find(|row| row.0 == "conduction").unwrap().2).unwrap();
    let bytes = ledger.get_artifact(&hash).unwrap().unwrap();
    J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}

#[test]
fn declared_contact_band_is_measured_by_actual_point_state_resolves() {
    assert!(crate::SOLVE_DRIVER_VERSION >= 46, "requires the complete budget integration");
    let f = Fixture::new();
    let project_path = f.inputs.join("contact.fsim");
    let mut spec = fs_project::parse_sexpr_migrating(&std::fs::read_to_string(&project_path).unwrap())
        .unwrap().decoded.spec;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    // Isolate the manufactured parameter. No material tolerance is invented.
    for material in spec.materials.as_mut().unwrap() { material.conductivity_tolerance = None; }
    let binding = &mut spec.interface_cards.as_mut().unwrap()[0];
    let fs_project::InterfaceState::DryContact { pressure, pressure_half_width, .. } = &mut binding.state
        else { panic!("dry-contact fixture") };
    pressure.value = 900000.0;
    pressure_half_width.value = 300000.0;
    std::fs::write(&project_path, fs_project::print_sexpr(&spec).unwrap()).unwrap();
    let source = SOURCE.replace("700000Pa", "600000Pa").replace("1100000Pa", "1200000Pa");
    let source_path = f.inputs.join("band.fsim");
    std::fs::write(&source_path, source).unwrap();
    let model = Model::load(&source_path).unwrap();
    let ledger = Ledger::open(f.root.join("band.db").to_str().unwrap()).unwrap();
    let gate = CancelGate::new_clock_free();
    let sample = model.sample(&ledger, &gate, &[900000.0, 5.0], 120.0).unwrap().unwrap();
    let nominal = receipt(&ledger, &sample);
    let propagation = nominal.get("propagation").unwrap();
    let term = propagation.get("parameters").unwrap();
    assert_eq!(term.str_field("state"), Some("measured"));
    assert_eq!(term.str_field("method"), Some("joint-state-resistance-corner-resolve"));
    let vertices = term.get("vertices").unwrap().as_array().unwrap();
    assert_eq!(vertices.len(), 2);
    assert!(term.str_field("detail").unwrap().contains("unprovided material uncertainty is not bounded"));
    // Bind a second real model with point states. Reapplying the original
    // half-width here would query a DIFFERENT manufacturing envelope.
    if let fs_project::InterfaceState::DryContact { pressure_half_width, .. }
        = &mut spec.interface_cards.as_mut().unwrap()[0].state { pressure_half_width.value = 0.0; }
    std::fs::write(&project_path, fs_project::print_sexpr(&spec).unwrap()).unwrap();
    let point_model = Model::load(&source_path).unwrap();
    let mut values = Vec::new();
    for pressure in [1200000.0, 600000.0] { // min resistance, max resistance
        let point = point_model.sample(&ledger, &gate, &[pressure, 5.0], 120.0).unwrap().unwrap();
        let result = receipt(&ledger, &point);
        values.push(result.get("propagation").unwrap().f64_field("nominal_base_k").unwrap());
        assert_ne!(result.get("propagation").unwrap().get("parameters").unwrap().str_field("method"),
            Some("joint-state-resistance-corner-resolve"), "zero-band models keep the legacy path");
    }
    let t0 = propagation.f64_field("nominal_base_k").unwrap();
    let width = values.iter().map(|t| (t - t0).abs()).fold(0.0, f64::max);
    assert!(width > 1e-6);
    assert!((term.f64_field("half_width_k").unwrap() - width).abs() < 1e-7);
    for (vertex, value) in vertices.iter().zip(values) {
        assert!((vertex.f64_field("t_max_k").unwrap() - value).abs() < 1e-7);
    }
    // The first model retained the original band despite replacing its file.
    assert!(matches!(&model.bound.base().interface_cards.as_ref().unwrap()[0].state,
        fs_project::InterfaceState::DryContact { pressure_half_width, .. } if pressure_half_width.value == 300000.0));
}
