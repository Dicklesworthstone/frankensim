//! Real imported bodies retain contact and radiation in the same steady solve.

use super::*;
#[path = "radiation/enclosure.rs"]
mod enclosure;

fn receipt(f: &Fixture, ordinal: usize) -> J {
    receipt_with_pack(f, ordinal, temperature_pack())
}

fn receipt_with_pack(f: &Fixture, ordinal: usize, pack: (Vec<u8>, String, String)) -> J {
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, pack);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let output = fs_cli::run(vec![
        "--json".into(),
        "solve".into(),
        source.display().to_string(),
        ledger.display().to_string(),
        "--materials".into(),
        f.dir.join("solid.fsmcdpk").display().to_string(),
        "--materials".into(),
        root.join("data/reference-project/gray-surface.fsmcdpk")
            .display()
            .to_string(),
        "--interfaces".into(),
        contact.display().to_string(),
    ]);
    assert_eq!(
        output.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        output.stdout,
        output.stderr
    );
    let run = J::parse(&output.stdout).unwrap();
    let run = artifact(&ledger, run.str_field("run_receipt").unwrap());
    let stage = run
        .get("stages")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("stage") == Some("conduction"))
        .unwrap();
    artifact(&ledger, stage.str_field("receipt").unwrap())
}

fn number(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn fixture(vacuum: bool) -> Fixture {
    let mut f = Fixture::new(true, true);
    // This is a nominal physical-composition fixture, with one fixed ambient.
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 300.0;
    envelope.ambient_hi.value = 300.0;
    if vacuum {
        f.project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .boundaries[0]
            .condition = fs_project::ThermalBoundaryCondition::HeatFlux {
            outward_flux: fs_qty::QtyAny::new(0.0, fs_project::spec::dims::HEAT_FLUX),
        };
    }
    f
}

fn check_contact(receipt: &J) {
    let contact = receipt.get("temperature_contact").unwrap();
    let radiation = receipt.get("radiation").unwrap();
    let coupling = contact.get("ambient_radiation").unwrap();
    assert_eq!(
        coupling.str_field("method"),
        Some("area-mean-secant-with-complete-contact-newton")
    );
    assert_eq!(contact.str_field("algorithm"), Some("coupled-newton"));
    assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
    assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
    assert_eq!(
        number(receipt.get("energy").unwrap(), "closure_w"),
        number(contact, "energy_closure_w")
    );
    assert_eq!(
        number(receipt.get("solver").unwrap(), "final_residual"),
        number(contact, "residual_norm_w")
    );
    assert_eq!(
        number(coupling, "thermal_responses"),
        number(radiation, "solid_solves")
    );
    assert!(number(contact, "total_thermal_responses") >= number(coupling, "thermal_responses"));
    assert!(number(contact, "total_nonlinear_iterations") >= number(radiation, "solid_iterations"));
    assert!(number(contact, "total_krylov_iterations") >= number(radiation, "krylov_iterations"));
    let states = contact.get("face_states").unwrap().as_array().unwrap();
    assert!(!states.is_empty());
    let (_, _, selected) = temperature_pack();
    for state in states {
        let mean = 0.5 * number(state, "side_a_temperature_k")
            + 0.5 * number(state, "side_b_temperature_k");
        assert_eq!(number(state, "query_temperature_k"), mean);
        assert!((number(state, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
        assert_eq!(state.str_field("selected_claim"), Some(selected.as_str()));
    }
    let interface_w: f64 = receipt
        .get("interfaces")
        .unwrap()
        .get("fluxes")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|row| number(row, "heat_rate_a_to_b_w"))
        .sum();
    assert!((interface_w.abs() - 4.44).abs() < 1e-7);
    let external =
        number(radiation, "convective_out_w") + number(radiation, "nonlinear_radiative_out_w");
    assert!(
        (external - 4.44).abs() < 1e-7,
        "contact is internal transfer, not extra external heat"
    );
    assert_eq!(receipt.get("solver_control"), Some(&J::Null));
}

#[test]
fn native_radiative_contact_cools_in_vacuum_without_a_fictitious_convection_bath() {
    let f = fixture(true);
    let receipt = receipt(&f, 80);
    check_contact(&receipt);
    let radiation = receipt.get("radiation").unwrap();
    assert_eq!(number(radiation, "convective_out_w"), 0.0);
    assert!((number(radiation, "nonlinear_radiative_out_w") - 4.44).abs() < 1e-7);
    assert_eq!(receipt.get("conjugate"), Some(&J::Null));
}

#[test]
fn native_radiative_contact_sends_only_original_convection_into_transported_air() {
    let f = fixture(false);
    let receipt = receipt(&f, 81);
    check_contact(&receipt);
    let radiation = receipt.get("radiation").unwrap();
    let air = receipt.get("conjugate").unwrap();
    let convective = number(radiation, "convective_out_w");
    let radiative = number(radiation, "nonlinear_radiative_out_w");
    assert!(
        radiative.abs() > 0.1,
        "the radiative reservoir must change this coupled heat flow"
    );
    assert!((number(air, "solid_total_w") - convective).abs() < 1e-7);
    assert!((number(air, "air_total_w") - convective).abs() <= number(air, "balance_tolerance_w"));
    assert!(
        (number(air, "air_total_w") - (convective + radiative)).abs() > 0.1,
        "radiation must exchange with its reservoir rather than heating the air stream"
    );
}

#[test]
fn source_surface_contact_reuses_ambient_radiation_and_preserves_the_air_heat_split() {
    for vacuum in [false, true] {
        let f = fixture(vacuum);
        // The full pressure/T source has the same nominal R(T) slice as the
        // independently selected curve. Its primal route remains Picard.
        let pack = temperature_pack_from_value(PropertyValue::Surface2D {
            axes: ["T".into(), "normal_pressure".into()],
            axis_dims: [fs_project::spec::dims::TEMPERATURE, fs_project::spec::dims::PRESSURE],
            samples: vec![([200.0, 0.8e6], 0.18), ([450.0, 0.8e6], 0.43),
                ([200.0, 1.2e6], 0.22), ([450.0, 1.2e6], 0.47)],
            triangles: vec![[0, 1, 3], [0, 3, 2]],
            dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        });
        let selected = pack.2.clone();
        let curve = receipt(&f, 90);
        let surface = receipt_with_pack(&f, 91, pack);
        let temperatures = |r: &J| artifact(&f.dir.join("temperature.db"),
            r.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone();
        let curve_t = temperatures(&curve);
        let surface_t = temperatures(&surface);
        for (a, b) in surface_t.as_array().unwrap().iter().zip(curve_t.as_array().unwrap()) {
            assert!((a.as_f64().unwrap() - b.as_f64().unwrap()).abs() < 2e-7);
        }
        let contact = surface.get("temperature_contact").unwrap();
        let composition = contact.get("ambient_radiation").unwrap();
        let radiation = surface.get("radiation").unwrap();
        assert_eq!(contact.str_field("algorithm"), Some("bounded-contact-picard"));
        assert_eq!(composition.str_field("method"), Some("area-mean-secant-with-bounded-contact-picard"));
        assert_eq!(composition.get("physical_endpoint_reassembled"), Some(&J::Bool(true)));
        assert_eq!(number(composition, "thermal_responses"), number(radiation, "solid_solves"));
        assert_eq!(number(composition, "radiation_responses"), number(radiation, "radiation_trials"));
        assert!(number(composition, "thermal_responses") > number(composition, "radiation_responses"));
        assert!(number(contact, "total_thermal_responses") >= number(composition, "thermal_responses"));
        assert!(number(contact, "total_nonlinear_iterations") >= number(radiation, "solid_iterations"));
        assert!(number(contact, "total_krylov_iterations") >= number(radiation, "krylov_iterations"));
        assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
        assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
        for state in contact.get("face_states").unwrap().as_array().unwrap() {
            let mean = 0.5 * number(state, "side_a_temperature_k") + 0.5 * number(state, "side_b_temperature_k");
            assert_eq!(number(state, "query_temperature_k"), mean);
            assert!((number(state, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
            assert_eq!(state.str_field("selected_claim"), Some(selected.as_str()));
        }
        let convection = number(radiation, "convective_out_w");
        let emission = number(radiation, "nonlinear_radiative_out_w");
        assert!((convection + emission - 4.44).abs() < 1e-7);
        if vacuum {
            assert_eq!(convection, 0.0);
            assert_eq!(surface.get("conjugate"), Some(&J::Null));
        } else {
            let air = surface.get("conjugate").unwrap();
            assert!(emission.abs() > 0.1);
            assert!((number(air, "solid_total_w") - convection).abs() < 1e-7);
            assert!((number(air, "air_total_w") - convection).abs() <= number(air, "balance_tolerance_w"));
        }
    }
}
