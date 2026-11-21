//! Point-sampled references retain their timing and identity through replay.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent, Proposal};
use fs_system::SystemBuilder;
use fs_system::elements::{ElementError, SignalClock, SignalSource, Waveform};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 9,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

#[test]
fn vector_references_are_point_samples_independent_of_trial_duration() {
    let wave = Waveform::sine(2.0, 0.5, 0.1, 3.0).unwrap();
    let source = SignalSource::new(vec![Waveform::constant(-4.0).unwrap(), wave.clone()]).unwrap();
    let old = source.decode_state(&[0.25_f64.to_bits(), 0]).unwrap();
    let stale = SignalClock::default();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for dt_s in [0.01, 0.1, 0.25] {
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s,
                inputs: &[],
            };
            let Observed::Stepped {
                outputs, receipt, ..
            } = check_component(&source, &case, cx).unwrap()
            else {
                panic!("smooth reference does not shorten");
            };
            assert_eq!(outputs, [-4.0, wave.value_at(0.25)]);
            assert_eq!(receipt.stored_change_j, 0.0);
            assert_eq!(receipt.dissipated_j, 0.0);
            let mut state = stale;
            source
                .propose(&old, dt_s, &[], &mut state, &mut [0.0; 2], cx)
                .unwrap();
            assert_eq!(state.time_s(), 0.25 + dt_s);
            check_port_tangent(&source, &case, 1e-6, 1e-9, cx).unwrap();
        }
    });
}

#[test]
fn table_events_and_compensated_clock_replay_with_the_system() {
    let reference = SignalSource::new(vec![
        Waveform::table(vec![(0.0, 2.0), (0.25, 2.0), (0.5, -1.0)]).unwrap(),
    ])
    .unwrap();
    let mut builder = SystemBuilder::new();
    let node = builder
        .add_with_tangent("reference", reference, SignalClock::default())
        .unwrap();
    let output = builder.port(node, "output").unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let first = system.step(1.0, cx).unwrap();
        assert_eq!(first.time_s, 0.25);
        assert_eq!(system.last_output(output), [2.0]);
        let checkpoint = system.checkpoint().unwrap();
        assert_eq!(system.step(1.0, cx).unwrap().time_s, 0.5);
        assert_eq!(system.last_output(output), [2.0]);
        system.step(0.125, cx).unwrap();
        assert_eq!(system.last_output(output), [-1.0]);
        let expected = *system.state::<SignalClock>(node).unwrap();
        system.restore(&checkpoint).unwrap();
        system.step(1.0, cx).unwrap();
        system.step(0.125, cx).unwrap();
        assert_eq!(system.state::<SignalClock>(node), Some(&expected));
        assert_eq!(system.ledger().energy_residual_j(), 0.0);
    });
}

#[test]
fn reference_identity_and_codec_reject_changed_or_invalid_state() {
    assert!(SignalSource::new(vec![]).is_err());
    let source = SignalSource::new(vec![Waveform::constant(2.0).unwrap()]).unwrap();
    let other = SignalSource::new(vec![Waveform::constant(-2.0).unwrap()]).unwrap();
    let (mut a, mut b) = (Vec::new(), Vec::new());
    source.identity_words(&mut a);
    other.identity_words(&mut b);
    assert_ne!(a, b);
    for words in [
        vec![],
        vec![0],
        vec![0, 0, 0],
        vec![f64::NAN.to_bits(), 0],
        vec![(-1.0_f64).to_bits(), 0],
        vec![0, f64::INFINITY.to_bits()],
    ] {
        assert!(source.decode_state(&words).is_none());
    }
    let clock = source
        .decode_state(&[0.5_f64.to_bits(), 2.0_f64.powi(-56).to_bits()])
        .unwrap();
    let mut words = Vec::new();
    assert!(source.encode_state(&clock, &mut words));
    assert_eq!(source.decode_state(&words), Some(clock));
}

#[test]
fn invalid_calls_and_cancellation_publish_no_clock_or_output() {
    let source = SignalSource::new(vec![Waveform::constant(2.0).unwrap()]).unwrap();
    let old = SignalClock::default();
    let stale = source.decode_state(&[0.5_f64.to_bits(), 0]).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for (dt, inputs) in [(0.0, &[][..]), (f64::NAN, &[][..]), (0.1, &[1.0][..])] {
            let mut new = stale;
            let mut out = [91.0];
            assert!(
                source
                    .propose(&old, dt, inputs, &mut new, &mut out, cx)
                    .is_err()
            );
            assert_eq!(new, stale);
            assert_eq!(out, [91.0]);
        }
        gate.request();
        let mut new = stale;
        let mut out = [91.0];
        assert_eq!(
            source.propose(&old, 0.1, &[], &mut new, &mut out, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            source.jvp(&old, 0.1, &[], &[], &mut out, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(new, stale);
        assert_eq!(out, [91.0]);
    });
}

#[test]
fn tangents_refuse_intervals_whose_primal_cannot_advance_or_sample() {
    let constant = SignalSource::new(vec![Waveform::constant(2.0).unwrap()]).unwrap();
    let stuck = constant.decode_state(&[1.0_f64.to_bits(), 0]).unwrap();
    let table = SignalSource::new(vec![
        Waveform::table(vec![(0.0, 1.0), (0.25, 2.0)]).unwrap(),
    ])
    .unwrap();
    let sine = SignalSource::new(vec![Waveform::sine(1.0, 1.0e149, 0.0, 0.0).unwrap()]).unwrap();
    let late = sine.decode_state(&[1.0e200_f64.to_bits(), 0]).unwrap();
    // The first clock cannot advance at this dt; the second interval crosses
    // a scheduled knot; the third clock advances but its sine phase overflows.
    let cases = [
        (constant, stuck, f64::MIN_POSITIVE, false),
        (table, SignalClock::default(), 0.5, true),
        (sine, late, 1.0e190, false),
    ];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (source, old, dt_s, shortens) in cases {
            let mut new = SignalClock::default();
            let mut output = [91.0];
            let primal = source.propose(&old, dt_s, &[], &mut new, &mut output, cx);
            if shortens {
                assert_eq!(primal, Ok(Proposal::Shorten { dt_s: 0.25 }));
            } else {
                assert!(primal.is_err());
            }
            assert_eq!(new, SignalClock::default());
            assert_eq!(output, [91.0]);
            assert!(source.jvp(&old, dt_s, &[], &[], &mut output, cx).is_err());
            assert_eq!(output, [91.0]);
            assert!(source.vjp(&old, dt_s, &[], &[1.0], &mut [], cx).is_err());
        }
    });
}
