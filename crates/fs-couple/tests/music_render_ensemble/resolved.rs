//! Original resolved-fluid pressure alongside independently rendered mechanics.
use super::{directory, modal, run, succeeded};
use fs_couple::nsf_dg::performance::file::ResolvedDgFilePerformance;
use fs_couple::pcm_wav::{encode_pcm16_wav, observation::{DecimatedRenderer, PressureRenderer}};
use fs_couple::render::schedule::force::file::ModalPerformance;
use fs_exec::CancelGate;
use std::sync::Arc;

fn input() -> String {
    include_str!("../../examples/resolved-dg-pressure.performance")
        .replace("output 2048 20", "output 512 20")
        .replace("hold 2048 0", "hold 512 0")
}

fn trace(mut renderer: impl PressureRenderer, samples: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; samples];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    pressure
}

#[test]
fn resolved_pressure_keeps_its_physical_clock_gas_and_receiver_in_a_mixed_phrase() {
    let dir = directory();
    let source_path = dir.join("resolved.performance");
    let modal_path = dir.join("modal.performance");
    let modal_text = modal(48000, 256, 200.0);
    std::fs::write(&modal_path, &modal_text).unwrap();
    let dry = trace(ModalPerformance::from_bytes(modal_text.as_bytes(), 37).unwrap().into_renderer(), 256);
    let mut waves = Vec::new();
    // Only the resolved source gas changes. Its original mesh determines a
    // different physical clock; its declared exterior and observation clock stay fixed.
    for (case, text) in [input(), input().replace("gas dry-air 300 100000", "gas dry-air 310 100000")]
        .into_iter().enumerate()
    {
        let p = ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 37,
            Arc::new(CancelGate::new_clock_free())).unwrap();
        let info = p.info();
        let renderer = p.renderer();
        let source_delay = renderer.baffled_receiver().unwrap().filter_delay_samples() / 2.0;
        let dt = renderer.system().zone().law().dt_s;
        let physical_steps = info.physical_intervals * renderer.system().limits().substeps as u64;
        let exterior = *renderer.system().load().radiation_record().unwrap().gas_record().state();
        let observed = DecimatedRenderer::new(p.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
        let delay = source_delay + observed.info().delay_output_samples;
        assert!(delay.fract() == 0.0 && delay < 256.0);
        let mut pressure = trace(observed, 256);
        assert!(pressure.iter().any(|p| p.abs() > 1e-6));
        // Independent physical traces, filtering and indexed pure delay form
        // the reference; the ensemble implementation is not its own oracle.
        for (sum, value) in pressure[delay as usize..].iter_mut().zip(&dry) { *sum += value; }
        let expected = encode_pcm16_wav(&pressure, 48000, 1.0).unwrap().0;
        std::fs::write(&source_path, text).unwrap();
        for block in [1, 37] {
            let output = dir.join(format!("resolved-{case}-{block}.wav"));
            succeeded(run(&output, &[("--resolved", &source_path), ("--modal", &modal_path)], block, true));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["\"kind\":\"resolved\"", "frankensim-resolved-dg-performance-v1",
                "original-NsfGrid3-D3Q27", "\"pressure_point\":\"baffled-outlet\"",
                "\"pressure_sample_rate_hz\":96000", "\"pressure_samples\":512",
                "\"compiled_controls\":3", "\"state_authority\":\"complete_gas_state\"",
                &format!("\"mechanics_sample_rate_hz\":{}", dt.recip()),
                &format!("\"mechanics_samples\":{physical_steps}"),
                &format!("\"density_kg_m3\":{:e}", exterior.density),
                &format!("\"source_filter_delay_output_samples\":{source_delay}"),
                &format!("\"common_delay_output_samples\":{delay}")]
            {
                assert!(metadata.contains(token), "{token}: {metadata}");
            }
        }
        waves.push(expected);
    }
    assert_ne!(waves[0][44..], waves[1][44..], "source gas must alter the mixed physical pressure");
}

#[test]
fn resolved_admission_refuses_incompatible_medium_clock_or_source_before_output() {
    let dir = directory();
    let a = dir.join("a.performance");
    let b = dir.join("b.performance");
    std::fs::write(&a, input()).unwrap();
    for (case, text, decimate, expected) in [
        (0, input(), false, "needs explicit --decimate"),
        (1, input().replace("dg_exterior 300 100000 0", "dg_exterior 310 100000 0"), true, "one ensemble has one listener gas"),
        (2, input().replace("output 512 20", "output 514 20").replace("hold 512 0", "hold 514 0"), true, "same physical duration"),
        (3, input().replace("gas dry-air 300 100000", "gas thermochemical 300 100000 mole 1 N2 1"), true, "original thermochemical source packs"),
        (4, "frankensim-flue-dg-performance-v1\n".into(), true, "unsupported resolved source schema"),
    ] {
        std::fs::write(&b, text).unwrap();
        let output = dir.join(format!("refused-{case}.wav"));
        let result = run(&output, &[("--resolved", &a), ("--resolved", &b)], 37, decimate);
        assert!(!result.status.success());
        let message = format!("{}{}", String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
        assert!(message.contains(expected), "{expected}: {message}");
        assert!(!output.exists() && !output.with_extension("provenance.json").exists());
    }
}

#[test]
fn resolved_internal_trace_stays_explicit_and_uses_the_original_pressure_samples() {
    let dir = directory();
    let text = input()
        .replace("observation baffled 0 0 0.2 8 32 2000", "observation inlet-pressure")
        .replace("observer_gas stokes-kirchhoff 4096 1000000 50000000\n", "");
    let source = dir.join("internal.performance");
    std::fs::write(&source, &text).unwrap();
    let p = ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 37,
        Arc::new(CancelGate::new_clock_free())).unwrap();
    let observed = DecimatedRenderer::new(p.into_renderer(), 96000, 48000, 37).unwrap();
    let expected = encode_pcm16_wav(&trace(observed, 256), 48000, 1.0).unwrap().0;
    let output = dir.join("internal.wav");
    succeeded(run(&output, &[("--resolved", &source)], 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(metadata.contains("internal resolved DG inlet trace; no exterior microphone"));
    assert!(metadata.contains("\"listener_medium\":{\"declared\":\"none\""));
    assert!(metadata.contains("\"pressure_point\":\"mean-inlet\""));
}
