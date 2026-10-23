//! Real acoustic solves, fitted pressure observers and optional passive loads.
//! Motions are prescribed numerical fixtures, not measured instruments.
use super::*;

fn source(image: radiation_spec::Image) -> Boundary {
    let sphere = SpherePanels::icosphere(0.02, 0).unwrap();
    let mut triangles = sphere.triangles().unwrap().to_vec();
    let weights = if image == radiation_spec::Image::Reciprocal {
        // Same successful thin-source fixture as reciprocal_playback_tests.
        for p in triangles.iter_mut().flatten() { p[2] *= 0.05; }
        let surface = SpherePanels::from_triangles(triangles.clone()).unwrap();
        surface.normals().iter().map(|n| n[2]).collect()
    } else {
        // Same breathing source as radiation_feedback_tests, including scale.
        vec![100.0; triangles.len()]
    };
    Boundary { triangles, weights: vec![weights], state_modes: vec![0] }
}

fn specification(image: radiation_spec::Image) -> radiation_spec::Spec {
    radiation_spec::Spec {
        image, band_hz: [40.0, 400.0], training_intervals: 16, max_order: 8,
        subdivisions: 0, max_panels: 80, max_dense_work: 2_000_000_000,
    }
}

fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|v| v.to_bits()).collect() }

fn same_bakes(a: &[Bake], b: &[Bake]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(bits(&[a.range_m, a.propagation_delay_s, a.pressure_gain,
            a.medium.density, a.medium.sound_speed]),
            bits(&[b.range_m, b.propagation_delay_s, b.pressure_gain,
                b.medium.density, b.medium.sound_speed]));
        assert_eq!(a.filters.len(), b.filters.len());
        for (a, b) in a.filters.iter().zip(&b.filters) {
            assert_eq!(a.n, b.n);
            assert_eq!(bits(&a.a), bits(&b.a));
            assert_eq!(bits(&a.b), bits(&b.b));
            assert_eq!(bits(&a.c), bits(&b.c));
            assert_eq!(bits(&[a.d, a.e_leftover, a.t_s]), bits(&[b.d, b.e_leftover, b.t_s]));
        }
    }
}

fn impulse_response(baked: &[Bake]) -> (Vec<u64>, Vec<u8>) {
    let mut observers: Vec<_> = baked.iter().map(|b| b.runtime().unwrap()).collect();
    let mut pressure = Vec::new();
    for tick in 0..128 {
        let acceleration = [if tick == 0 { 1.0 } else { 0.0 }];
        for observer in &mut observers { pressure.push(observer.step(&acceleration).unwrap()); }
    }
    assert!(pressure.iter().any(|p| *p != 0.0), "non-vacuous pressure impulse");
    let wav = encode_pcm16_wav_interleaved(&pressure, OUTPUT_RATE, baked.len() as u16, 20.0).unwrap().0;
    (bits(&pressure), wav)
}

#[test]
fn parallel_real_bem_preserves_stereo_filters_pcm_and_passive_loads() {
    let receivers = [Receiver::FinitePoint([0.0, 0.0, 0.12]), Receiver::FinitePoint([0.1, 0.0, 0.08])];
    for image in [radiation_spec::Image::Collocation, radiation_spec::Image::Reciprocal] {
        let boundary = source(image);
        let spec = specification(image);
        let gate = CancelGate::new_clock_free();
        // The collocation fixture exercises actual feedback fitting. The thin
        // reciprocal fixture exercises the one-way branch with NO load bank.
        let loaded = image == radiation_spec::Image::Collocation;
        let (serial, serial_load) = bake_scene_image_with_threads(
            &boundary, &receivers, spec, Medium::air(), &gate, loaded, false, 1).unwrap();
        let (parallel, parallel_load) = bake_scene_image_with_threads(
            &boundary, &receivers, spec, Medium::air(), &gate, loaded, false, 4).unwrap();
        same_bakes(&serial, &parallel);
        assert_eq!(impulse_response(&serial), impulse_response(&parallel));
        assert_eq!(serial_load.is_some(), loaded);
        assert_eq!(parallel_load.is_some(), loaded);
        if let (Some(a), Some(b)) = (serial_load, parallel_load) {
            assert_eq!(a.ports, b.ports);
            assert_eq!(a.poles.len(), b.poles.len());
            for (a, b) in a.poles.iter().zip(&b.poles) {
                assert_eq!(bits(&[a.omega, a.zeta]), bits(&[b.omega, b.zeta]));
                assert_eq!(bits(&a.coupling), bits(&b.coupling));
            }
        }
    }
}

#[test]
fn parallel_real_bem_keeps_budget_band_and_cancellation_refusals() {
    let receivers = [Receiver::FinitePoint([0.0, 0.0, 0.12])];
    for image in [radiation_spec::Image::Collocation, radiation_spec::Image::Reciprocal] {
        let boundary = source(image);
        let spec = specification(image);
        for cfg in [radiation_spec::Spec { max_dense_work: 1, ..spec },
            radiation_spec::Spec { band_hz: [40.0, 20_000.0], ..spec }] {
            let gate = CancelGate::new_clock_free();
            let refusal = |threads| bake_scene_image_with_threads(
                &boundary, &receivers, cfg, Medium::air(), &gate, false, false, threads)
                .err().expect("bounded work or unresolved band must refuse").to_string();
            assert_eq!(refusal(1), refusal(4));
        }
        let gate = CancelGate::new_clock_free();
        gate.request();
        for threads in [1, 4] {
            let error = bake_scene_image_with_threads(
                &boundary, &receivers, spec, Medium::air(), &gate, false, false, threads)
                .err().expect("cancelled preparation must refuse");
            assert_eq!(error.to_string(), "radiation preparation cancelled");
        }
    }
}
