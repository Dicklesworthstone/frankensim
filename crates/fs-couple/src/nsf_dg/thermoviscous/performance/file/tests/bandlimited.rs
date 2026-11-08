//! The authored filter moves with its original cold source into the complete field.
use super::*;
use crate::pcm_wav::observation::resample::BandlimitedObservationPolicy;

#[test]
fn authored_filter_keeps_its_admitted_policy_channels_and_window_through_field_handoff() {
    let text = input().replace(
        "output 384 1\n",
        "output 384 1\nobservation_filter bandlimited 2000 12000 128 1048576 65536 65536\n",
    );
    let source = read(&text, 37);
    let original = source.source().renderer();
    let filter = original.bandlimited_observation().unwrap();
    let policy = filter.policy();
    assert_eq!(
        policy,
        BandlimitedObservationPolicy {
            passband_hz: 2000.0,
            stopband_hz: 12000.0,
            phase_count: 128,
            max_preparation_terms: 1 << 20,
            max_history_values: 1 << 16,
            max_output_terms: 1 << 16,
        }
    );
    let report = filter.report();
    let physical_intervals = filter.source_intervals();
    let samples = filter.output_samples();
    let zone_hash = original.system().zone().state_hash();
    let duct_state = original.system().load().state().clone();
    let inventory = original.system().finite_reservoir_inventory().unwrap();
    assert_eq!(filter.accepted_source_intervals(), 0);
    assert_eq!(filter.rendered_samples(), 0);
    assert!(
        report.channels >= 3,
        "pressure, terminal and original escape tiles must all survive"
    );
    let result = source
        .with_physical_renderer(|renderer, info| {
            let transferred = renderer.bandlimited_observation().unwrap();
            assert_eq!(transferred.policy(), policy);
            assert_eq!(transferred.report(), report);
            assert_eq!(transferred.source_intervals(), physical_intervals);
            assert_eq!(transferred.output_samples(), samples);
            assert_eq!(transferred.accepted_source_intervals(), 0);
            assert_eq!(transferred.rendered_samples(), 0);
            assert_eq!(info.samples, samples);
            assert_eq!(renderer.remaining_samples(), samples);
            assert_eq!(renderer.system().zone().state_hash(), zone_hash);
            assert_eq!(
                renderer.system().load().system().load().state(),
                &duct_state
            );
            assert_eq!(
                renderer.system().finite_reservoir_inventory(),
                Some(inventory)
            );
            let mut pressure = vec![0.0; samples as usize];
            for chunk in pressure.chunks_mut(37) {
                renderer.block(chunk).unwrap();
            }
            assert!(pressure.iter().all(|p| p.is_finite()));
            assert!(pressure.iter().any(|p| p.abs() > 0.0));
            assert_eq!(renderer.remaining_samples(), 0);
            assert_eq!(
                renderer.system().accepted_steps(),
                physical_intervals as u64
            );
            let transferred = renderer.bandlimited_observation().unwrap();
            assert_eq!(transferred.accepted_source_intervals(), physical_intervals);
            assert_eq!(transferred.rendered_samples(), samples);
            assert_eq!(transferred.report(), report);
            Ok((
                pressure,
                renderer.system().zone().state_hash(),
                renderer.system().finite_reservoir_inventory().unwrap(),
            ))
        })
        .unwrap();
    assert_eq!(result, render(read(&text, 23), 23));
}
