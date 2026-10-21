use super::*;

const REACTION: &str = include_str!("../../../../examples/flue-mouth-reaction.performance");
fn source() -> String {
    REACTION.replace("output 24000 5", "output 1024 1")
}
fn render(text: &str, block: usize) -> (Vec<f64>, f64) {
    let file =
        FlueFilePerformance::from_bytes(text.as_bytes(), block, &CancelGate::new_clock_free())
            .unwrap();
    let mut pressure = vec![0.0; file.info().samples as usize];
    let mut renderer = file.into_renderer();
    for part in pressure.chunks_mut(block) {
        renderer.block(part).unwrap();
    }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-10));
    for d in renderer.voice().diagnostics() {
        let scale = d.source_work_j.abs()
            + d.pipe_work_j.abs()
            + d.radiation_loss_j
            + d.separation_loss_j
            + d.stored_delta_j.abs();
        assert!(d.ledger_defect_j.abs() <= 1e-10 * scale.max(1e-30));
        assert!(d.radiation_loss_j >= 0.0);
    }
    (pressure, renderer.voice().mouth_flow_m3_s())
}

#[test]
fn source_selected_rectangle_reaction_reaches_the_normal_pcm_owner() {
    let text = source();
    let file = FlueFilePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new_clock_free())
        .unwrap();
    let renderer = file.renderer();
    let load = renderer.voice().mouth_radiation_load().unwrap();
    assert_eq!(
        load.exterior(),
        &renderer.radiation_input().unwrap().exterior
    );
    assert_ne!(
        load.exterior().temperature,
        renderer.voice().gas().temperature,
        "reaction uses the actual exterior, not the inlet region"
    );
    assert_eq!(
        load.size_m(),
        [file.geometry().flue_width_m, file.geometry().cut_up_m]
    );
    assert_eq!(load.cells(), renderer.mouth_radiation().unwrap().cells);
    assert_eq!(
        load.maximum_frequency_hz(),
        renderer.config().receiver.maximum_frequency_hz
    );
    let expected = render(&text, 37);
    for block in [1, 137] {
        assert_eq!(render(&text, block), expected);
    }
    let wav = encode_pcm16_wav(&expected.0, 48_000, 1.0).unwrap().0;
    let dir = directory();
    let input = dir.join("rectangular-reaction.performance");
    std::fs::write(&input, &text).unwrap();
    let output = dir.join("rectangular-reaction.wav");
    success(run(&input, &output, 137, false));
    assert_eq!(std::fs::read(&output).unwrap(), wav);
    let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in [
        "\"model\":\"compact-rectangular-baffled-v1\"",
        "\"legacy_load_replaced\":true",
        "\"reaction_added\":true",
        "\"radiation_feedback_added\":true",
        "\"mutual_mouth_terminal_reaction\":false",
        "same-sample compact rectangular mouth reaction",
    ] {
        assert!(sidecar.contains(token), "{token}: {sidecar}");
    }
}

#[test]
fn exterior_gas_data_changes_the_source_motion_and_pcm_not_only_propagation() {
    let text = source();
    let alternate = text.replace(
        "exterior same",
        "exterior mixture 293.15 101325 mole 2 N2 0.5 CO2 0.5",
    );
    let base = render(&text, 37);
    let changed = render(&alternate, 37);
    assert_ne!(
        base.1.to_bits(),
        changed.1.to_bits(),
        "the exterior gas must feed back into the actual mouth motion"
    );
    assert_ne!(base.0, changed.0);
    assert_ne!(
        encode_pcm16_wav(&base.0, 48_000, 1.0).unwrap().0,
        encode_pcm16_wav(&changed.0, 48_000, 1.0).unwrap().0
    );
    // An observer-only exterior change cannot change the source mechanics.
    let legacy = text.replace("mouth_reaction compact-baffled\n", "");
    let legacy_alternate = alternate.replace("mouth_reaction compact-baffled\n", "");
    assert_eq!(
        render(&legacy, 37).1.to_bits(),
        render(&legacy_alternate, 37).1.to_bits()
    );
    assert_ne!(base.1.to_bits(), render(&legacy, 37).1.to_bits());
    let none = legacy.replace("controls 4096", "mouth_reaction none\ncontrols 4096");
    assert_eq!(render(&legacy, 37), render(&none, 37));
}

#[test]
fn mismatched_or_double_counted_mouth_declarations_refuse_before_output() {
    let text = source();
    for invalid in [
        text.replace(
            "flue 0.010 0.001 0.004 0 0 0",
            "flue 0.010 0.001 0.004 0 0 0.02",
        ),
        text.replace("mouth_radiation -0.3 0 8 4\n", ""),
        text.replace("mouth_radiation -0.3 0 8 4", "mouth_radiation -0.3 0 33 33"),
        text.replace("exterior same", "exterior 293.15 80000 0"),
        text.replace("mouth_reaction compact-baffled", "mouth_reaction guessed"),
    ] {
        assert!(
            FlueFilePerformance::from_bytes(invalid.as_bytes(), 37, &CancelGate::new_clock_free())
                .is_err()
        );
    }
    let dir = directory();
    let input = dir.join("duplicate-mouth-loss.performance");
    std::fs::write(
        &input,
        text.replace(
            "flue 0.010 0.001 0.004 0 0 0",
            "flue 0.010 0.001 0.004 0 0 0.02",
        ),
    )
    .unwrap();
    let output = dir.join("refused.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
