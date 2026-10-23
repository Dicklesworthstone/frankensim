//! Manufactured upstream-facing geometry through real material/duct playback.
use super::*;

fn records(side: &str) -> String {
    let faces = [6,7,8,9,10,11,18,19,20,21,22,23].map(|triangle|
        format!("wetted_facet {triangle} {side}\n")).concat();
    format!("wetted_surfaces 1\nwetted_surface 0 1000000 manufactured-exposed-tip\nwetted_facets 12\n{faces}")
}
fn source(side: &str) -> String {
    input().replace("audio 768000 48000 512 0.005", "audio 768000 48000 256 0.005")
        .replacen("compile_limits", &format!("{}compile_limits", records(side)), 1)
}

#[test]
fn original_wetted_faces_and_material_changes_reach_the_actual_microphone() {
    let text = source("positive");
    let performance = build(&text, 37);
    let original = build(&input(), 37);
    let model = valve(&performance);
    let port = model.material_plates()[0].wetted_surface().unwrap();
    assert_eq!(port.surface().provenance(), "manufactured-exposed-tip");
    assert_eq!(port.surface().upstream_facets().len(), 12);
    assert!((port.area_m2()-0.000125).abs() < 1e-18);
    assert_eq!(model.viscosity_matrix_n_s_m(), valve(&original).viscosity_matrix_n_s_m());
    assert_eq!(model.modal_states(), valve(&original).modal_states());
    for (a,b) in model.spec().modes.iter().zip(&valve(&original).spec().modes) {
        assert_eq!(a.mass_kg, b.mass_kg);
        assert_eq!(a.stiffness_n_m, b.stiffness_n_m);
        assert_eq!(a.gap_projection, b.gap_projection);
        assert_eq!(a.force_projection, b.force_projection);
    }
    let (pressure, renderer) = render(performance);
    let (partitioned, _) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert!(renderer.source().last_frame().unwrap().balance_residual_j().abs() < 1e-14);
    let (flipped, _) = render(build(&source("negative"), 37));
    let (material, _) = render(build(
        &text.replacen("section blade.matpack", "section altered.matpack", 1), 37));
    assert!(pressure.iter().zip(flipped).any(|(a,b)| (a-b).abs() > 1e-15));
    assert!(pressure.iter().zip(material).any(|(a,b)| (a-b).abs() > 1e-15));
}

#[test]
fn wetted_geometry_enters_full_viscosity_and_both_material_memory_operators() {
    let tensor = "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000";
    for (propagator, loss) in [
        ("midpoint", tensor),
        ("exact-zoh", "loss proportional-maxwell relaxed 1 boundary_modulus boundary_time"),
        ("exact-zoh", "loss regional-maxwell relaxed 1 boundary_modulus boundary_time"),
    ] {
        let text = source("negative").replace("propagator exact-zoh", &format!("propagator {propagator}"))
            .replace(tensor, loss).replacen("state 0 0", "state 0.00000001 0", 1);
        let performance = build(&text, 37);
        if loss != tensor {
            let baseline = build(&text.replace(&records("negative"), ""), 37);
            let memory = valve(&performance).material_relaxation().unwrap().memory_sqrt_j();
            assert_eq!(memory, valve(&baseline).material_relaxation().unwrap().memory_sqrt_j());
            assert!(memory.iter().flatten().any(|value| *value != 0.0));
        }
        let (pressure, renderer) = render(performance);
        assert!(pressure.iter().any(|p| p.abs() > 1e-15));
        let frame = renderer.source().last_frame().unwrap();
        assert!(frame.balance_residual_j().abs() < 1e-14);
        let LoadedAperture::Modal(model) = renderer.source().system().valve() else { panic!("original modal owner"); };
        assert!(model.material_plates()[0].wetted_surface().is_some());
        if loss != tensor {
            assert!(model.material_relaxation().unwrap().memory_sqrt_j().iter().flatten().any(|v| *v != 0.0));
        }
    }
}

#[test]
fn wetted_faces_coexist_with_distributed_contact_and_independent_pressure_ports() {
    let text = include_str!("../../../../../../examples/material_modal_valve/paired-pressure.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
        .replacen("pressure_shapes", &format!("{}pressure_shapes", records("negative")), 1);
    let (_, p) = render(build(&text, 37));
    let LoadedAperture::Modal(model) = p.source().system().valve() else { panic!("modal source"); };
    assert!(model.material_pressure_ports().iter().all(Option::is_some));
    assert!(model.material_plates()[0].wetted_surface().is_some());
    let text = include_str!("../../../../../../examples/material_modal_valve/paired-spatial.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
        .replacen("closure spatial", &format!("{}closure spatial", records("positive")), 1);
    let (_, p) = render(build(&text, 37));
    let LoadedAperture::Modal(model) = p.source().system().valve() else { panic!("modal source"); };
    assert!(model.plate_closure().is_some());
    assert!(model.material_plates()[0].wetted_surface().is_some());
}

#[test]
fn omitted_and_zero_wetted_bindings_preserve_original_playback_bitwise() {
    let text = input().replace("audio 768000 48000 512 0.005", "audio 768000 48000 128 0.005");
    let (a, _) = render(build(&text, 37));
    let (b, _) = render(build(&text.replacen("compile_limits", "wetted_surfaces 0\ncompile_limits", 1), 37));
    assert_eq!(a, b);
}

#[test]
fn acoustic_clock_admission_uses_actual_exposed_geometry_during_material_construction() {
    // A deliberately narrow lossless upstream passage raises Z A A^T for a
    // whole exposed blade. It does not load a physically sheltered pressure face.
    // Material damping and structural angular admission remain active in both.
    let text = input().replace("propagator exact-zoh", "propagator midpoint")
        .replace("duct_section 1 0 0.12 0.007 0.002", "duct_section 1 0 0.12 0.00005 0.002");
    let refused = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false))).err().unwrap();
    assert!(matches!(refused, PlateValveInputError::Physics(
        crate::acoustic_realize::AcousticRealizeError::TimeResolution(_))));
    let sheltered = "wetted_surfaces 2\nwetted_surface 0 1000000 manufactured-sheltered\nwetted_facets 0\nwetted_surface 1 1000000 manufactured-sheltered\nwetted_facets 0\n";
    let loaded = build(&text.replacen("compile_limits", &format!("{sheltered}compile_limits"), 1), 37);
    assert!(valve(&loaded).spec().modes.iter().all(|mode| mode.pressure_area_m2 == 0.0));
    assert!(valve(&loaded).viscosity_matrix_n_s_m().is_some());
    assert_eq!(loaded.info().propagator, BladePropagator::Midpoint);
}

#[test]
fn invalid_face_geometry_refuses_before_audio() {
    let text = source("positive");
    for changed in [
        text.replace("wetted_surface 0", "wetted_surface 2"),
        text.replace("wetted_facet 6 positive", "wetted_facet 24 positive"),
        text.replace("wetted_facet 7 positive", "wetted_facet 6 negative"),
        text.replace("wetted_facet 6 positive", "wetted_facet 6 channel"),
        text.replace("wetted_surface 0 1000000", "wetted_surface 0 1"),
    ] {
        assert!(MaterialModalValvePerformance::from_bytes_with_sources(changed.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false))).is_err());
    }
}
