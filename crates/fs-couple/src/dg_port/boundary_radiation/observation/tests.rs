use super::*;
use fs_vfit::observation::FirKernel;

#[test]
fn accepted_aperture_flow_stream_preserves_original_fir_delays_and_transactions() {
    let mut load = super::super::tests::load();
    let source = load.boundary_radiation_record().unwrap();
    // Deliberately hand-authored routing oracle, not a physical BEM bake or
    // qualification. The mechanical load itself is the actual native BEM bake.
    let fir = ObservationFir::new(
        2,
        4,
        vec![
            FirKernel {
                input: 0,
                output: 0,
                integer_delay_samples: 1,
                taps: vec![2.0, -0.5],
            },
            FirKernel {
                input: 1,
                output: 0,
                integer_delay_samples: 0,
                taps: vec![-1.0, 0.25],
            },
            FirKernel {
                input: 0,
                output: 1,
                integer_delay_samples: 0,
                taps: vec![3.0],
            },
        ],
    )
    .unwrap();
    let mut observer = DgBoundaryRadiationObserver {
        fir,
        ports: load.terminal_ports().to_vec(),
        declaration_identity: source.declaration_identity(),
        realization_identity: source.realization_identity(),
        transfer_declaration_hash: ContentHash([1; 32]),
        time_step_s: load.propagator().dt(),
        observed_steps: 0,
        receivers: vec![ReceiverSpec {
            id: "routing-oracle".into(),
            position_m: [0.1, 0.1, 0.1],
        }],
        lower_window: None,
        pass_band_upper_hz: 1000.0,
        filter_latency_samples: 1,
    };
    assert!(observer.step(&load).is_err());
    {
        let abandoned = load.preview_step(1e-9).unwrap();
        assert_eq!(abandoned.state().step_count(), 1);
    }
    assert!(observer.step(&load).is_err());
    let mut previous = [0.0; 2];
    let mut preceding_left = 0.0;
    for sample in 0..32 {
        let mut replay = observer.clone();
        load.preview_step(if sample < 7 { 1e-9 } else { 0.0 })
            .unwrap()
            .accept();
        let frames = load.last_terminal_frames();
        let flows = [frames[0].outward_flow_m3_s, frames[1].outward_flow_m3_s];
        let actual = observer.step(&load).unwrap();
        assert_eq!(actual, replay.step(&load).unwrap());
        let pressure = ((2.0 * previous[0] - 0.5 * preceding_left) - flows[1]) + 0.25 * previous[1];
        let scale = 2.0 * previous[0].abs()
            + 0.5 * preceding_left.abs()
            + flows[1].abs()
            + 0.25 * previous[1].abs();
        assert!((actual[0] - pressure).abs() <= 4e-15 * scale.max(1e-300));
        assert_eq!(actual[1], 3.0 * flows[0]);
        assert_eq!(&actual[2..], &[0.0, 0.0]);
        assert!(
            observer.step(&load).is_err(),
            "a duplicate step cannot replay radiation twice"
        );
        assert_eq!(observer.observed_steps(), sample + 1);
        preceding_left = previous[0];
        previous = flows;
    }
    load.preview_step(0.0).unwrap().accept();
    load.preview_step(0.0).unwrap().accept();
    assert!(
        observer.step(&load).is_err(),
        "discarded mechanical history cannot be fabricated"
    );
    assert_eq!(observer.observed_steps(), 32);
}
