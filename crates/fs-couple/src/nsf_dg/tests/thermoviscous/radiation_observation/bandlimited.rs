//! All original source bases pass through the filter before the live receivers.
use super::*;
use crate::pcm_wav::observation::resample::{BandlimitedObservationPolicy, BandlimitedResampler};

#[test]
fn every_terminal_and_reservoir_tile_is_filtered_before_rayleigh_superposition() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = || attach(system(&field, medium.clone(), &primary, 256, &gate), &gate);
    let rate = ACOUSTIC_RATE * 3 / 2;
    let samples = 256;
    let policy = BandlimitedObservationPolicy {
        passband_hz: 300.0,
        stopband_hz: 1000.0,
        phase_count: 128,
        max_preparation_terms: 1 << 24,
        max_history_values: 1 << 18,
        max_output_terms: 1 << 18,
    };
    let mut whole = observed(make(), rate, samples, samples as usize, Arc::clone(&gate))
        .with_bandlimited_observation(policy)
        .unwrap();
    let mut parts = observed(make(), rate, samples, 7, Arc::clone(&gate))
        .with_bandlimited_observation(policy)
        .unwrap();
    let mut direct = make();
    let mut bore = whole.baffled_receiver().unwrap().clone();
    let mut mouths: Vec<_> = [1, 2]
        .map(|aperture| {
            (
                aperture,
                whole.reservoir_receiver_at(aperture).unwrap().clone(),
            )
        })
        .to_vec();
    let terminal_count = primary.terminal_ports().len();
    let channels = 1
        + terminal_count
        + mouths
            .iter()
            .map(|(_, receiver)| receiver.input_count())
            .sum::<usize>();
    let mut filter =
        BandlimitedResampler::new(direct.load().duration_s(), rate, channels, samples, policy)
            .unwrap();
    assert_eq!(
        whole.bandlimited_observation().unwrap().report(),
        filter.report()
    );
    let controls = whole.drive_history().to_vec();
    let mut filtered = vec![0.0; channels];
    let mut expected = Vec::new();
    let mut different_surface_flows = false;
    for _ in 0..samples {
        let required = filter.required_source_intervals().unwrap();
        while (direct.accepted_steps() as usize) < required {
            let frame = direct
                .step_with_drive(controls[direct.accepted_steps() as usize], &gate)
                .unwrap();
            let mut values = vec![
                direct
                    .load()
                    .system()
                    .last_frame()
                    .unwrap()
                    .mean_pressure_pa,
            ];
            values.extend(
                direct
                    .load()
                    .system()
                    .load()
                    .last_terminal_frames()
                    .iter()
                    .map(|terminal| terminal.outward_flow_m3_s),
            );
            let mut first_flows = Vec::new();
            for (aperture, receiver) in &mouths {
                let flows: Vec<_> = frame
                    .reservoir_flows
                    .iter()
                    .filter(|flow| flow.aperture == *aperture)
                    .map(|flow| flow.outward_volume_flow_m3_s)
                    .collect();
                assert_eq!(flows.len(), receiver.input_count());
                first_flows.push(flows[0]);
                values.extend(flows);
            }
            different_surface_flows |= first_flows[0] != first_flows[1];
            filter.push_frame(&values).unwrap();
        }
        filter.render_next(&mut filtered).unwrap();
        let mut next = 1 + terminal_count;
        let mut pressure = bore.step_trace(&filtered[1..next]).unwrap();
        for (_, receiver) in &mut mouths {
            let end = next + receiver.input_count();
            pressure += receiver.step_trace(&filtered[next..end]).unwrap();
            next = end;
        }
        expected.push(pressure);
    }
    let mut a = vec![0.0; samples as usize];
    let mut b = a.clone();
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(7) {
        parts.block(chunk).unwrap();
    }
    assert_eq!(a, expected);
    assert_eq!(a, b);
    assert!(
        a.iter().any(|pressure| pressure.abs() > 0.0),
        "cover filter and receiver latency"
    );
    assert!(
        different_surface_flows,
        "independent mouths must keep distinct accepted flow"
    );
    assert_eq!(whole.system().accepted_steps() as usize, controls.len());
    assert_eq!(whole.samples_rendered(), samples);
    assert_eq!(whole.remaining_samples(), 0);
    same_physical_owners(whole.system(), &direct);
    same_physical_owners(whole.system(), parts.system());
}
