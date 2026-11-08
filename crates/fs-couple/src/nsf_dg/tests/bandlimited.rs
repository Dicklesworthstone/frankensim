//! The original resolved/DG owners remain on their physical clock.
use super::*;
use crate::pcm_wav::observation::resample::{BandlimitedObservationPolicy, BandlimitedResampler};

fn policy() -> BandlimitedObservationPolicy {
    BandlimitedObservationPolicy {
        passband_hz: 12000.0,
        stopband_hz: 24000.0,
        phase_count: 128,
        max_preparation_terms: 1 << 20,
        max_history_values: 1 << 16,
        max_output_terms: 1 << 16,
    }
}

fn filtered(
    rate: u32,
    samples: u64,
    block: usize,
    gate: Arc<CancelGate>,
) -> NsfDgPerformance<D3q27> {
    let original = system(100000.0, true, 1, DgBoundary::Matched);
    let count = source_intervals_for(original.load().propagator().dt(), rate, samples).unwrap();
    NsfDgPerformance::new_resampled(
        original,
        ApertureObservation::Inlet,
        rate,
        vec![0.0; count],
        samples,
        block,
        gate,
    )
    .unwrap()
    .with_bandlimited_observation(policy())
    .unwrap()
}

#[test]
fn bandlimited_callbacks_read_the_same_original_physical_trajectory() {
    let gate = Arc::new(CancelGate::new_clock_free());
    for (rate, samples) in [(72000, 48), (144000, 96)] {
        let mut direct = system(100000.0, true, 1, DgBoundary::Matched);
        let mut filter =
            BandlimitedResampler::new(direct.load().propagator().dt(), rate, 1, samples, policy())
                .unwrap();
        let mut expected = vec![0.0; samples as usize];
        for output in &mut expected {
            let required = filter.required_source_intervals().unwrap();
            while (direct.accepted_steps() as usize) < required {
                let frame = direct.step(0.0, &gate).unwrap();
                filter.push_frame(&[frame.inlet_pressure_pa]).unwrap();
            }
            filter.render_next(core::slice::from_mut(output)).unwrap();
        }
        let mut whole = filtered(rate, samples, samples as usize, Arc::clone(&gate));
        let mut parts = filtered(rate, samples, 13, Arc::clone(&gate));
        let mut a = vec![0.0; samples as usize];
        let mut b = a.clone();
        whole.block(&mut a).unwrap();
        for chunk in b.chunks_mut(13) {
            parts.block(chunk).unwrap();
        }
        assert_eq!(a, expected);
        assert_eq!(a, b);
        assert!(
            a.iter().any(|p| p.abs() > 1.0e-9),
            "retain actual nonzero source pressure"
        );
        for renderer in [&whole, &parts] {
            assert_eq!(
                renderer.system().zone().state_hash(),
                direct.zone().state_hash()
            );
            assert_eq!(renderer.system().load().state(), direct.load().state());
            assert_eq!(renderer.system().accepted_steps(), 64);
            assert_eq!(renderer.samples_rendered(), samples);
            assert_eq!(renderer.remaining_samples(), 0);
            assert_eq!(
                renderer.bandlimited_observation().unwrap().report(),
                filter.report()
            );
        }
    }
}

#[test]
fn bandlimited_cached_output_honors_cancellation_without_touching_the_prefix() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let mut renderer = filtered(RATE * 3, 6, 3, Arc::clone(&gate));
    let mut prefix = [123.0];
    renderer.block(&mut prefix).unwrap();
    assert_eq!(renderer.system().accepted_steps(), 1);
    let before = (
        renderer.system().zone().state_hash(),
        renderer.system().load().state().clone(),
    );
    gate.request();
    let mut output = [456.0, 789.0];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(output, [456.0, 789.0]);
    assert_eq!(
        before,
        (
            renderer.system().zone().state_hash(),
            renderer.system().load().state().clone()
        )
    );
    assert_eq!(renderer.samples_rendered(), 1);
    assert_eq!(
        renderer
            .bandlimited_observation()
            .unwrap()
            .rendered_samples(),
        1
    );
    assert!(matches!(
        renderer.block(&mut output),
        Err(crate::render::RenderError::Poisoned)
    ));
}
