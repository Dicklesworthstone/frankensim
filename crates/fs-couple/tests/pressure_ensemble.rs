//! G3: real modal trajectories, existing decimators and one physical output clock.
use std::io::Cursor;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::pcm_wav::observation::ensemble::{PressureEnsemble, PressureEnsembleConfig};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_couple::render::{ControlDelta, ModalStringVoice, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_couple::render::schedule::force::ensemble::{EnsembleConfig, EnsembleRenderer};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> PressureEnsembleConfig {
    PressureEnsembleConfig { sample_rate_hz: 48_000, samples: 257, max_block: 257, max_parts: 4 }
}
fn source(rate: u32, capacity: usize, fail_at: Option<u64>) -> ScheduledRenderer {
    let mut model = ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: 1000.0, damping_ratio: 0.01,
        pressure_per_modal_velocity: C64::new(1.0, 0.1),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: 0.0001, velocity_m_sqrt_kg_per_s: 0.0,
    }]).unwrap();
    let voice = ModalStringVoice::new(model, vec![0.0]).unwrap();
    let mut controls: Vec<_> = [(0, 3.0), (5, -2.0), (8, 0.0), (149, 1.25), (1028, 0.0)]
        .into_iter().map(|(sample, force)| ScheduledControl { sample, delta: ControlDelta::SetModalForce {
            voice: 0, mode: 0, force_n_per_sqrt_kg: force,
        }}).collect();
    if let Some(sample) = fail_at {
        controls.push(ScheduledControl { sample, delta: ControlDelta::SetModalForce {
            voice: 0, mode: 0, force_n_per_sqrt_kg: 1e100,
        }});
    }
    ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::ModalString(voice)], capacity), controls, 6).unwrap()
}
fn observed(ratio: u32) -> DecimatedRenderer {
    DecimatedRenderer::new(source(48_000 * ratio, 257, None), 48_000 * ratio, 48_000, 257).unwrap()
}
fn mixed(ratios: &[u32]) -> PressureEnsemble {
    PressureEnsemble::new(ratios.iter().copied().map(observed).collect(), config()).unwrap()
}

// Independent whole-history filtering and indexed shifts, not the mixer's
// ring-buffer implementation or DecimatedRenderer's grouped callback loop.
fn expected(ratios: &[u32], samples: usize) -> Vec<f64> {
    let common = ratios.iter().map(|&r| Decimator::new(r as usize, 1).unwrap().delay_output_frames() as usize)
        .max().unwrap();
    let mut sum = vec![0.0; samples];
    for &ratio in ratios {
        let mut physical = source(48_000 * ratio, samples * ratio as usize, None);
        let mut input = vec![0.0; samples * ratio as usize];
        physical.block(&mut input).unwrap();
        let mut filter = Decimator::new(ratio as usize, 1).unwrap();
        let delay = common - filter.delay_output_frames() as usize;
        for (i, chunk) in input.chunks_exact(ratio as usize).enumerate() {
            let p = filter.preview(chunk).unwrap()[0]; filter.commit();
            if i + delay < samples { sum[i + delay] += p; }
        }
    }
    sum
}
fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|v| v.to_bits()).collect() }

#[test]
fn mixed_solver_clocks_preserve_controls_and_align_observation_latency() {
    let ratios = [1, 2, 3, 4];
    let reference = expected(&ratios, 257);
    assert!(reference.iter().any(|p| p.abs() > 1e-5));
    for partition in [1, 7, 37, 257] {
        let mut mix = mixed(&ratios);
        assert_eq!(mix.delay_output_samples(), 80);
        assert_eq!(mix.part_info().iter().map(|p| p.alignment_delay_samples).collect::<Vec<_>>(),
            vec![80, 40, 0, 36]);
        let mut out = vec![0.0; 257];
        for block in out.chunks_mut(partition) { mix.block(block).unwrap(); }
        assert_eq!(bits(&out), bits(&reference));
        assert_eq!(mix.remaining_samples(), 0);
        for (part, ratio) in mix.parts().iter().zip(ratios) {
            assert_eq!(part.samples_rendered(), 257);
            assert_eq!(part.source().samples_rendered(), 257 * u64::from(ratio));
            assert_eq!(part.source().applied_controls().iter().map(|e| e.sample).collect::<Vec<_>>(),
                vec![0, 5, 8, 149]);
            assert_eq!(part.source().pending_controls()[0].sample, 1028);
        }
    }
}

#[test]
fn identity_paths_introduce_neither_gain_nor_delay() {
    let mut a = source(48_000, 257, None);
    let mut out = [0.0; 257]; a.block(&mut out).unwrap();
    let reference: Vec<_> = out.iter().map(|v| 0.0 + *v + *v).collect();
    let mut mix = mixed(&[1, 1]);
    assert_eq!(mix.delay_output_samples(), 0);
    for block in out.chunks_mut(37) { mix.block(block).unwrap(); }
    assert_eq!(bits(&out), bits(&reference));
}

#[test]
fn streaming_pcm_cancellation_retains_alignment_and_all_filter_histories() {
    let (expected_wav, clips) = encode_pcm16_wav(&expected(&[1, 4], 257), 48_000, 0.05).unwrap();
    let mut mix = mixed(&[1, 4]);
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 0.05, 37).unwrap();
    let mut scratch = [0.0; 37];
    render_pressure_pcm16(&mut mix, &mut stream, &CancelGate::new(), &mut scratch, 37).unwrap();
    let pending = mix.parts()[1].source().pending_controls().to_vec();
    assert_eq!(pending[0].sample, 149); // next high-rate group, not an output boundary
    let before = stream.sink().get_ref().clone();
    let gate = CancelGate::new(); gate.request(); scratch.fill(12345.0);
    assert_eq!(render_pressure_pcm16(&mut mix, &mut stream, &gate, &mut scratch, 220).unwrap(),
        ScheduledWavProgress::Cancelled { samples: 0 });
    assert_eq!(scratch, [12345.0; 37]);
    assert_eq!(stream.sink().get_ref(), &before);
    assert_eq!(mix.parts()[1].source().samples_rendered(), 148);
    assert_eq!(mix.parts()[1].source().pending_controls(), pending.as_slice());
    render_pressure_pcm16(&mut mix, &mut stream, &CancelGate::new(), &mut scratch, 220).unwrap();
    let (output, summary) = stream.finish().unwrap();
    assert_eq!(output.into_inner(), expected_wav);
    assert_eq!(summary.clipped_samples, clips as u64);
    assert_eq!(mix.samples_rendered(), 257);
}

#[test]
fn full_window_and_shape_admission_precede_output_or_any_source_motion() {
    let mut mix = mixed(&[1, 4]);
    let mut out = [12345.0; 258];
    assert!(mix.block(&mut out).is_err());
    assert!(mix.validate_sample_count(258).is_err());
    assert!(mix.block(&mut []).is_err());
    assert_eq!(out, [12345.0; 258]);
    assert_eq!(mix.samples_rendered(), 0);
    assert!(mix.parts().iter().all(|p| p.source().applied_controls().is_empty()));
    mix.block(&mut out[..257]).unwrap();
    assert_eq!(bits(&out[..257]), bits(&expected(&[1, 4], 257)));
    out[0] = 12345.0;
    assert!(mix.block(&mut out[..1]).is_err());
    assert_eq!(out[0], 12345.0);
}

#[test]
fn construction_refuses_clock_mismatch_missing_history_and_storage_limits() {
    for cfg in [
        PressureEnsembleConfig { sample_rate_hz: 44_100, ..config() },
        PressureEnsembleConfig { max_parts: 1, ..config() },
        PressureEnsembleConfig { max_block: 258, ..config() },
        PressureEnsembleConfig { samples: u64::MAX, ..config() },
    ] { assert!(PressureEnsemble::new(vec![observed(1), observed(4)], cfg).is_err()); }
    assert!(PressureEnsemble::<ScheduledRenderer>::new(Vec::new(), config()).is_err());
    let mut running = observed(4); running.block(&mut [0.0]).unwrap();
    assert!(PressureEnsemble::new(vec![observed(1), running], config()).is_err());
}

#[test]
fn mixed_producer_types_keep_finite_source_horizons_through_rate_conversion() {
    let finite = |samples| EnsembleRenderer::from_parts(vec![source(96_000, 257, None)], EnsembleConfig {
        sample_rate_hz: 96_000, samples, max_block: 257, max_voices: 1, max_events: 6,
    }).unwrap();
    let build = |horizon| {
        let sources: Vec<(Box<dyn PressureRenderer>, u32)> = vec![
            (Box::new(source(48_000, 257, None)), 48_000), (Box::new(finite(horizon)), 96_000),
        ];
        let parts = sources.into_iter().map(|(s, rate)| DecimatedRenderer::new(s, rate, 48_000, 257).unwrap()).collect();
        PressureEnsemble::new(parts, config())
    };
    assert!(build(513).is_err(), "all 514 high-rate steps must be admitted at construction");
    let mut mix = build(514).unwrap();
    let mut out = vec![0.0; 257]; mix.block(&mut out).unwrap();
    assert_eq!(bits(&out), bits(&expected(&[1, 2], 257)));
}

#[test]
fn late_physical_refusal_poisons_the_mix_without_publishing_a_partial_pcm_callback() {
    let failing = DecimatedRenderer::new(source(192_000, 257, Some(12)), 192_000, 48_000, 257).unwrap();
    let mut mix = PressureEnsemble::new(vec![observed(1), failing], config()).unwrap();
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 1.0, 2).unwrap();
    assert!(render_pressure_pcm16(&mut mix, &mut stream, &CancelGate::new(), &mut [0.0; 2], 4).is_err());
    assert_eq!(stream.samples_written(), 2);
    assert_eq!(mix.samples_rendered(), 2);
    assert_eq!(mix.parts()[0].samples_rendered(), 4, "earlier part has moved; no partial retry");
    let mut sentinel = [12345.0; 2];
    assert!(matches!(mix.block(&mut sentinel), Err(RenderError::Poisoned)));
    assert!(matches!(mix.validate_sample_count(1), Err(RenderError::Poisoned)));
    assert_eq!(sentinel, [12345.0; 2]);
}
