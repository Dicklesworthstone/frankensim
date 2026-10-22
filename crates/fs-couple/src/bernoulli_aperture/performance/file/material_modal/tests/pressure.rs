use super::*;

const PRESSURE: &str = include_str!("../../../../../../examples/material_modal_valve/paired-pressure.performance");
fn source() -> String {
    PRESSURE.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
        .replace("pressure-pa\t5e0", "pressure-pa\t0e0")
}

#[test]
fn prescribed_asymmetric_pressure_drives_actual_sound_and_material_response() {
    let text = source();
    let a = build(&text, 37);
    let v = valve(&a);
    let first = v.material_pressure_ports()[0].as_ref().unwrap();
    assert_eq!(first.shape().nodal_shape().len(), 21);
    assert_eq!(first.shape().provenance(), "manufactured-balanced-pressure-gradient");
    assert!(first.shape().nodal_shape().iter().any(|p| *p < 0.0));
    assert!(first.shape().nodal_shape().iter().any(|p| *p > 0.0));
    assert!(first.modal_area_m2().iter().any(|a| a.abs() > 1e-8));
    let (pressure, source) = render(a);
    let (partitioned, _) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(encode_pcm16_wav(&pressure, 48000, 0.005).unwrap(),
        encode_pcm16_wav(&partitioned, 48000, 0.005).unwrap());
    let (changed, _) = render(build(&text.replacen("section blade.matpack", "section altered.matpack", 1), 37));
    assert!(pressure.iter().zip(changed).any(|(a,b)| (a-b).abs() > 1e-15));
    assert!(source.source().system().spatial_pressure_work_j().iter().all(|w| w.is_finite() && *w != 0.0));
    assert_eq!(source.held_material_pressures_pa(), [0.0; 2]);
    assert_eq!(source.applied_material_pressure_controls(0).len(), 4);
    assert_eq!(source.applied_material_pressure_controls(1).len(), 2);
    assert!(source.pending_material_pressure_controls(0).is_empty());
    assert!(source.source().last_frame().unwrap().balance_residual_j().abs() < 1e-14);
}

#[test]
fn pressure_controls_run_inside_audio_frames_and_keep_each_source_work_distinct() {
    let mut p = build(&source(), 37).into_renderer();
    assert!(p.set_spatial_pressures([1.0, 0.0]).is_err());
    p.block(&mut [0.0]).unwrap();
    assert_eq!(p.source().mechanical_samples_rendered(), 16);
    assert_eq!(p.held_material_pressures_pa(), [10.0, 3.0]);
    assert_eq!(p.applied_material_pressure_controls(0).len(), 2);
    assert_eq!(p.applied_material_pressure_controls(1).len(), 1);
    let frame = p.source().last_frame().unwrap();
    let pressure_work: f64 = frame.spatial_pressure.iter().map(|p| p.work_j).sum();
    assert!((pressure_work-frame.aperture.mechanical_work_j).abs() < 1e-22);
    assert_eq!(frame.source_work_j, 0.0); // no copied blowing-pressure source
    assert!(p.source().system().spatial_pressure_work_j().iter().all(|w| *w != 0.0));
}

#[test]
fn spatial_pressure_coexists_with_regional_memory_and_distributed_closure() {
    let text = source().replace("loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000",
        "loss regional-maxwell relaxed 1 boundary_modulus boundary_time");
    let (trace, source) = render(build(&text, 37));
    assert!(trace.iter().any(|p| p.abs() > 1e-15));
    let LoadedAperture::Modal(v) = source.source().system().valve() else { panic!("modal source"); };
    assert!(v.material_relaxation().unwrap().memory_sqrt_j().iter().flatten().any(|z| *z != 0.0));
    // The existing spatial-contact file provides the original geometry. Insert
    // the same independent pressure records before its closure and compiler.
    let spatial = include_str!("../../../../../../examples/material_modal_valve/paired-spatial.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005");
    let records = PRESSURE.split("pressure_shapes 2\n").nth(1).unwrap();
    let shapes = records.split("source_pressures ").next().unwrap();
    let events = records.split("source_pressures ").nth(1).unwrap().split("compile_limits").next().unwrap();
    let combined = spatial.replacen("closure spatial", &format!("pressure_shapes 2\n{shapes}closure spatial"), 1)
        .replacen("compile_limits", &format!("source_pressures {events}compile_limits"), 1);
    let (_, p) = render(build(&combined, 37));
    let LoadedAperture::Modal(v) = p.source().system().valve() else { panic!("modal source"); };
    assert!(v.plate_closure().is_some());
    assert!(v.material_pressure_ports().iter().all(Option::is_some));
}

#[test]
fn pressure_frontdoor_refuses_missing_duplicate_nonfinite_and_out_of_window_data() {
    let text = source();
    for altered in [
        text.replace("pressure_nodes 21", "pressure_nodes 20"),
        text.replace("pressure_shape 1 ", "pressure_shape 0 "),
        text.replace("pressure_shape 1 ", "pressure_shape 2 "),
        text.replace("pressure_node -1", "pressure_node NaN"),
        text.replace("source_pressure 1 0 8", "source_pressure 4096 0 8"),
        text.replace("source_pressure 5 1 3", "source_pressure 5 2 3"),
        text.replace("pressure_shape 0 1000000", "pressure_shape 0 1"),
        text.replace("compile_limits 100000 8192", "compile_limits 1 8192"),
    ] {
        assert!(MaterialModalValvePerformance::from_bytes_with_sources(altered.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false))).is_err());
    }
}
