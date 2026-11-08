//! Selected physical duct pressure, without retiming NSF or either acoustic field.
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

#[test]
fn complete_field_filter_keeps_selected_duct_and_original_source_clock() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let source = playback_field(&mesh, &unit);
    let rate = RATE * 3 / 4;
    let samples = 48;
    let controls = vec![NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0)); 64];
    let make = |block| {
        NsfThermoviscousDgPerformance::new_resampled_with_drives(
            playback_system(&source, &gate),
            0,
            ApertureObservation::Inlet,
            rate,
            controls.clone(),
            samples,
            block,
            Arc::clone(&gate),
        )
        .unwrap()
        .with_bandlimited_observation(policy())
        .unwrap()
    };
    let mut whole = make(48);
    let mut parts = make(7);
    let mut direct = playback_system(&source, &gate);
    let mut filter =
        BandlimitedResampler::new(direct.load().duration_s(), rate, 1, samples, policy()).unwrap();
    let mut expected = [0.0; 48];
    for output in &mut expected {
        let required = filter.required_source_intervals().unwrap();
        while (direct.accepted_steps() as usize) < required {
            direct
                .step_with_drive(controls[direct.accepted_steps() as usize], &gate)
                .unwrap();
            let frame = direct.load().system().last_frame().unwrap();
            let selected = frame
                .additional_ducts
                .iter()
                .find(|duct| duct.port == 0)
                .unwrap();
            filter.push_frame(&[selected.mean_pressure_pa]).unwrap();
        }
        filter.render_next(core::slice::from_mut(output)).unwrap();
    }
    let mut a = [0.0; 48];
    let mut b = [0.0; 48];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(7) {
        parts.block(chunk).unwrap();
    }
    assert_eq!(a, expected);
    assert_eq!(a, b);
    assert!(a.iter().any(|p| p.abs() > 1.0e-9));
    assert_eq!(snapshot!(whole.system()), snapshot!(&direct));
    assert_eq!(snapshot!(parts.system()), snapshot!(&direct));
    assert_eq!(whole.system().accepted_steps(), 64);
    assert_eq!(whole.samples_rendered(), samples);
    assert_eq!(whole.remaining_samples(), 0);
    assert_eq!(
        whole.bandlimited_observation().unwrap().report(),
        filter.report()
    );
}

#[test]
fn complete_field_filter_refuses_cached_outputs_after_cancellation() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let source = playback_field(&mesh, &unit);
    let mut renderer = NsfThermoviscousDgPerformance::new_resampled_with_drives(
        playback_system(&source, &gate),
        0,
        ApertureObservation::Inlet,
        RATE * 3,
        vec![NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0)); 2],
        6,
        3,
        Arc::clone(&gate),
    )
    .unwrap()
    .with_bandlimited_observation(policy())
    .unwrap();
    renderer.block(&mut [0.0]).unwrap();
    assert_eq!(renderer.system().accepted_steps(), 1);
    let before = snapshot!(renderer.system());
    gate.request();
    let mut suffix = [123.0, 456.0];
    assert!(renderer.block(&mut suffix).is_err());
    assert_eq!(suffix, [123.0, 456.0]);
    assert_eq!(snapshot!(renderer.system()), before);
    assert_eq!(renderer.samples_rendered(), 1);
    assert_eq!(
        renderer
            .bandlimited_observation()
            .unwrap()
            .rendered_samples(),
        1
    );
    assert!(matches!(
        renderer.block(&mut suffix),
        Err(crate::render::RenderError::Poisoned)
    ));
}
