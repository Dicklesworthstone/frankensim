use std::io::Cursor;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_couple::render::{ControlDelta, ModalStringVoice, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn source(rate: u32, capacity: usize) -> ScheduledRenderer {
    let mut model = ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: 1000.0, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::new(1.0,0.1),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: 0.0001, velocity_m_sqrt_kg_per_s: 0.0,
    }]).unwrap();
    let voice = ModalStringVoice::new(model, vec![0.0]).unwrap();
    let events = [(0,3.0),(5,-2.0),(8,0.0),(149,1.25),(1028,0.0)].into_iter()
        .map(|(sample,force)| ScheduledControl { sample, delta: ControlDelta::SetModalForce {
            voice:0,mode:0,force_n_per_sqrt_kg:force,
        }}).collect();
    ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::ModalString(voice)],capacity),events,5).unwrap()
}
fn reference(ratio: usize, samples: usize) -> Vec<f64> {
    let mut source = source(48_000 * ratio as u32, ratio*samples);
    let mut input = vec![0.0;ratio*samples]; source.block(&mut input).unwrap();
    let mut filter = Decimator::new(ratio,1).unwrap();
    input.chunks_exact(ratio).map(|chunk| {
        let value=filter.preview(chunk).unwrap()[0]; filter.commit(); value
    }).collect()
}

#[test]
fn source_sample_events_and_pressure_match_the_original_filter_across_partitions() {
    let expected=reference(4,257);
    for capacity in [1,3,4,17] {
        for block in [1,7,37,257] {
            let mut r=DecimatedRenderer::new(source(192_000,capacity),192_000,48_000,257).unwrap();
            let mut actual=vec![0.0;257];
            for chunk in actual.chunks_mut(block) { r.block(chunk).unwrap(); }
            assert_eq!(actual.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),
                expected.iter().map(|v|v.to_bits()).collect::<Vec<_>>());
            assert_eq!(r.samples_rendered(),257);assert_eq!(r.source().samples_rendered(),1028);
            assert_eq!(r.source().applied_controls().iter().map(|e|e.sample).collect::<Vec<_>>(),vec![0,5,8,149]);
            assert_eq!(r.source().pending_controls()[0].sample,1028);
            assert_eq!(r.info().delay_output_samples,44.0);
            assert_eq!(r.info().first_output_source_index,3);
        }
    }
}

#[test]
fn bypass_preserves_the_original_callback_output_and_control_logs_exactly() {
    let mut direct=source(48_000,512);
    let mut bypass=DecimatedRenderer::new(source(48_000,512),48_000,48_000,512).unwrap();
    for count in [1,37,64,155] {
        let mut expected=vec![0.0;count];let mut actual=expected.clone();
        direct.block(&mut expected).unwrap();bypass.block(&mut actual).unwrap();
        assert_eq!(actual,expected);
        assert_eq!(direct.context().control_log(),bypass.source().context().control_log());
    }
    assert_eq!(bypass.info().filter_profile,"identity");
    assert_eq!(bypass.info().delay_output_samples,0.0);
}

#[test]
fn bad_clocks_fractional_intervals_and_oversized_callbacks_refuse_before_motion() {
    for (input,output) in [(192_000,0),(192_000,44_100),(192_000,384_000),(192_000,1000),(96_000,48_000)] {
        assert!(DecimatedRenderer::new(source(192_000,1),input,output,17).is_err());
    }
    let mut advanced=source(192_000,1);advanced.block(&mut [0.0]).unwrap();
    assert!(DecimatedRenderer::new(advanced,192_000,48_000,17).is_err());
    let mut r=DecimatedRenderer::new(source(192_000,1),192_000,48_000,17).unwrap();
    assert!(r.output_samples_for(69).is_err());
    assert_eq!(r.output_samples_for(68).unwrap(),17);
    assert!(r.validate_sample_count(u64::MAX).is_err());
    assert!(r.validate_sample_rate(192_000).is_err());
    assert!(r.block(&mut []).is_err());let mut untouched=[7.0;18];
    assert!(r.block(&mut untouched).is_err());assert_eq!(untouched,[7.0;18]);
    assert_eq!(r.samples_rendered(),0);assert_eq!(r.source().samples_rendered(),0);
    assert!(r.source().applied_controls().is_empty());
    r.block(&mut [0.0;17]).unwrap();assert_eq!(r.source().samples_rendered(),68);
}

#[test]
fn stream_cancellation_preserves_filter_history_and_pending_source_events() {
    let expected=reference(4,257);
    let (wav,clips)=encode_pcm16_wav(&expected,48_000,1.0).unwrap();
    let mut r=DecimatedRenderer::new(source(192_000,3),192_000,48_000,37).unwrap();
    let mut stream=Pcm16WavStream::new(Cursor::new(Vec::new()),48_000,1.0,37).unwrap();
    let mut scratch=[0.0;37];
    render_pressure_pcm16(&mut r,&mut stream,&CancelGate::new(),&mut scratch,37).unwrap();
    assert_eq!(r.source().samples_rendered(),148);
    let pending=r.source().pending_controls().to_vec();
    let bytes=stream.sink().get_ref().clone();
    let gate=CancelGate::new();gate.request();scratch.fill(7.0);
    assert_eq!(render_pressure_pcm16(&mut r,&mut stream,&gate,&mut scratch,220).unwrap(),
        ScheduledWavProgress::Cancelled {samples:0});
    assert_eq!(scratch,[7.0;37]);assert_eq!(stream.sink().get_ref(),&bytes);
    assert_eq!(r.source().pending_controls(),pending.as_slice());
    assert_eq!(r.source().samples_rendered(),148);
    render_pressure_pcm16(&mut r,&mut stream,&CancelGate::new(),&mut scratch,220).unwrap();
    let (output,summary)=stream.finish().unwrap();
    assert_eq!(output.into_inner(),wav);assert_eq!(summary.samples,257);assert_eq!(summary.clipped_samples,clips as u64);
}

#[test]
fn late_physics_failure_poisoning_keeps_only_complete_output_callbacks() {
    let mut model=ModalAcousticTimeModel::try_new(192_000,vec![ModalAcousticMode {
        angular_frequency_rad_s:1000.0,damping_ratio:0.0,pressure_per_modal_velocity:C64::new(1.0,0.0),
    }],ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState {displacement_m_sqrt_kg:0.0001,velocity_m_sqrt_kg_per_s:0.0}]).unwrap();
    let voice=ModalStringVoice::new(model,vec![0.0]).unwrap();
    let source=ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::ModalString(voice)],4),vec![
        ScheduledControl {sample:12,delta:ControlDelta::SetModalForce {voice:0,mode:0,force_n_per_sqrt_kg:1e100}}
    ],1).unwrap();
    let mut r=DecimatedRenderer::new(source,192_000,48_000,2).unwrap();
    let mut stream=Pcm16WavStream::new(Cursor::new(Vec::new()),48_000,1.0,2).unwrap();
    assert!(render_pressure_pcm16(&mut r,&mut stream,&CancelGate::new(),&mut [0.0;2],4).is_err());
    assert_eq!(stream.samples_written(),2);assert_eq!(r.samples_rendered(),2);
    assert_eq!(r.source().samples_rendered(),12);
    let live = [ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 0.0 }];
    assert!(matches!(r.validate_controls(&live), Err(RenderError::Poisoned)));
    assert!(matches!(r.apply_controls(&live), Err(RenderError::Poisoned)));
    let mut untouched=[7.0;2];assert!(matches!(r.block(&mut untouched),Err(RenderError::Poisoned)));
    assert_eq!(untouched,[7.0;2]);
    let (output,summary)=stream.finish().unwrap();assert_eq!(summary.samples,2);assert_eq!(output.into_inner().len(),48);
}

#[test]
fn actual_high_rate_modal_tone_is_filtered_before_it_can_fold_into_audio() {
    let tone=|frequency:f64| {
        let mut model=ModalAcousticTimeModel::try_new(96_000,vec![ModalAcousticMode {
            angular_frequency_rad_s:core::f64::consts::TAU*frequency,damping_ratio:0.0,
            pressure_per_modal_velocity:C64::new(1.0,0.0),
        }],ModalAcousticTimeBudget::audible_reference()).unwrap();
        model.restore_states(&[ModalAcousticState {displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:1.0}]).unwrap();
        let source=ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::ModalString(
            ModalStringVoice::new(model,vec![0.0]).unwrap())],2),vec![],0).unwrap();
        let mut r=DecimatedRenderer::new(source,96_000,48_000,1536).unwrap();
        let mut out=vec![0.0;1536];r.block(&mut out).unwrap();
        (out[512..].iter().map(|v|v*v).sum::<f64>()/1024.0).sqrt()
    };
    let passband=tone(6000.0);let folded=tone(36_000.0);
    assert!((passband-std::f64::consts::FRAC_1_SQRT_2).abs()<2e-4,"{passband}");
    assert!(folded<1e-5,"{folded}");
    // Directly selecting every second 36-kHz sample yields a 12-kHz alias
    // with RMS sqrt(1/2), so this test would fail for unfiltered subsampling.
    let naive=(0..1024).map(|i| (core::f64::consts::TAU*36_000.0*(i+1) as f64/48_000.0).cos().powi(2))
        .sum::<f64>()/1024.0;
    assert!(naive.sqrt()>0.7);
}

#[test]
fn live_inputs_reach_retained_mechanics_through_borrow_erasure_and_decimation() {
    use fs_couple::render::schedule::force::ensemble::{EnsembleConfig, EnsembleRenderer};

    let live = ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 7.0 };
    let oracle = source(96_000, 1);
    let events = oracle.pending_controls().to_vec();
    let mut direct = oracle.into_context();
    let mut mechanics = vec![0.0; 256];
    for (sample, value) in mechanics.iter_mut().enumerate() {
        if sample == 64 { direct.apply_controls(&[live]).unwrap(); }
        for event in events.iter().filter(|event| event.sample == sample as u64) {
            direct.apply_controls(&[event.delta]).unwrap();
        }
        direct.block(std::slice::from_mut(value)).unwrap();
    }
    let mut filter = Decimator::new(2, 1).unwrap();
    let expected: Vec<_> = mechanics.chunks_exact(2).map(|chunk| {
        let value = filter.preview(chunk).unwrap()[0];
        filter.commit();
        value
    }).collect();

    let mut retained = EnsembleRenderer::from_parts(vec![source(96_000, 1)], EnsembleConfig {
        sample_rate_hz: 96_000, max_block: 1, samples: 256, max_voices: 1, max_events: 5,
    }).unwrap();
    let mut actual = vec![0.0; 128];
    {
        let erased: Box<dyn PressureRenderer + '_> = Box::new(&mut retained);
        let mut decimated = DecimatedRenderer::new(erased, 96_000, 48_000, 32).unwrap();
        let renderer: &mut dyn PressureRenderer = &mut decimated;
        renderer.block(&mut actual[..32]).unwrap();
        renderer.validate_controls(&[live]).unwrap();
        renderer.apply_controls(&[live]).unwrap();
        assert_eq!(renderer.samples_rendered(), 32);
        for chunk in actual[32..].chunks_mut(13) { renderer.block(chunk).unwrap(); }
        assert_eq!(renderer.samples_rendered(), 128);
        assert_eq!(decimated.source().samples_rendered(), 256);
    }
    assert_eq!(actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
    assert_eq!(retained.renderer().context().control_log(), direct.control_log());
    assert_eq!(retained.renderer().pending_controls()[0].sample, 1028);
    assert_eq!(retained.remaining_samples(), 0);
}

#[test]
fn invalid_live_batches_refuse_atomically_without_losing_source_events_or_history() {
    let mut renderer = DecimatedRenderer::new(source(96_000, 1), 96_000, 48_000, 16).unwrap();
    let mut reference = DecimatedRenderer::new(source(96_000, 1), 96_000, 48_000, 16).unwrap();
    renderer.block(&mut [0.0; 6]).unwrap();
    reference.block(&mut [0.0; 6]).unwrap();
    let log = renderer.source().context().control_log().to_vec();
    let pending = renderer.source().pending_controls().to_vec();
    let applied = renderer.source().applied_controls().to_vec();
    let deltas = [
        ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 7.0 },
        ControlDelta::SetModalForce { voice: 1, mode: 0, force_n_per_sqrt_kg: 2.0 },
    ];
    assert!(matches!(renderer.validate_controls(&deltas), Err(RenderError::UnknownVoice { index: 1 })));
    assert!(matches!(renderer.apply_controls(&deltas), Err(RenderError::UnknownVoice { index: 1 })));
    assert_eq!(renderer.samples_rendered(), 6);
    assert_eq!(renderer.source().samples_rendered(), 12);
    assert_eq!(renderer.source().context().control_log(), log);
    assert_eq!(renderer.source().pending_controls(), pending);
    assert_eq!(renderer.source().applied_controls(), applied);
    let mut actual = [0.0; 16];
    let mut expected = [0.0; 16];
    renderer.block(&mut actual).unwrap();
    reference.block(&mut expected).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(renderer.source().context().control_log(), reference.source().context().control_log());
}

#[test]
fn pressure_only_producers_refuse_inputs_and_keep_their_stream_usable() {
    use fs_couple::pcm_wav::observation::ensemble::{PressureEnsemble, PressureEnsembleConfig};

    let part = DecimatedRenderer::new(source(96_000, 1), 96_000, 48_000, 8).unwrap();
    let mut renderer = PressureEnsemble::new(vec![part], PressureEnsembleConfig {
        sample_rate_hz: 48_000, samples: 8, max_block: 8, max_parts: 1,
    }).unwrap();
    let deltas = [ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 7.0 }];
    assert!(matches!(renderer.validate_controls(&deltas), Err(RenderError::Control { .. })));
    assert!(matches!(renderer.apply_controls(&deltas), Err(RenderError::Control { .. })));
    renderer.validate_controls(&[]).unwrap();
    renderer.apply_controls(&[]).unwrap();
    assert_eq!(renderer.samples_rendered(), 0);
    let mut actual = [0.0; 8];
    renderer.block(&mut actual).unwrap();
    assert_eq!(actual.as_slice(), reference(2, 8));
}
