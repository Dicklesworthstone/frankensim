use super::*;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};

const N: usize = 32;
const BIN: usize = 3;
const LIMITS: NetworkCoherenceBudget = NetworkCoherenceBudget {
    maximum_rank: 64,
    maximum_complex_products: 10_000,
};

fn recorded(drives: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<_> = drives
        .iter()
        .flat_map(|a| {
            (0..N).map(move |j| {
                let phase = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                a.map(|z| z.re * phase.cos() + z.im * phase.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        1.0 / 4096.0,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 4096,
            maximum_transforms: 128,
            maximum_workspace_bytes: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-12 * b.abs().max(1.0), "{a:e} != {b:e}");
}

impl PreparedNetworkSources {
    // Exact retained maps of a common pressure node: an open-pressure source
    // through Rs=2, an ideal flow source, and an exterior resistance Rl=3.
    // This is an analytical transfer fixture, not an executed graph assembly.
    fn recorded_test_mixer(omega: f64) -> Self {
        let internal = vec![
            C64::from_re(0.4),
            C64::from_re(-1.2),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ];
        Self {
            omega,
            ports: vec!["a_pressure".into(), "b_flow".into(), "z_load".into()],
            coordinates: vec![
                NetworkSourceCoordinate {
                    label: "a_pressure".into(),
                    quantity: NetworkSourceQuantity::OpenCircuitPressure,
                },
                NetworkSourceCoordinate {
                    label: "b_flow".into(),
                    quantity: NetworkSourceQuantity::VolumeFlow,
                },
            ],
            defaults: vec![C64::ZERO; 2],
            source_ports: vec![0, 1],
            source_index: vec![Some(0), Some(1), None],
            canonical: vec![0, 1, 2],
            pressure: [0.6, 1.2, 0.6, 1.2, 0.6, 1.2].map(C64::from_re).to_vec(),
            flow: [0.2, -0.6, 0.0, 1.0, -0.2, -0.4].map(C64::from_re).to_vec(),
            dissipative: internal.clone(),
            internal,
        }
    }

    // Lossless uniform tube with p(0)=x and inward U(L)=q. These original
    // source coordinates also bind the independently analytical field tests.
    pub(in crate::network) fn recorded_test_tube(omega: f64, angle: f64, z: f64) -> Self {
        let sec = 1.0 / angle.cos();
        let tan = angle.tan();
        Self {
            omega,
            ports: vec!["left_pressure".into(), "right_flow".into()],
            coordinates: vec![
                NetworkSourceCoordinate {
                    label: "left_pressure".into(),
                    quantity: NetworkSourceQuantity::OpenCircuitPressure,
                },
                NetworkSourceCoordinate {
                    label: "right_flow".into(),
                    quantity: NetworkSourceQuantity::VolumeFlow,
                },
            ],
            defaults: vec![C64::ZERO; 2],
            source_ports: vec![0, 1],
            source_index: vec![Some(0), Some(1)],
            canonical: vec![0, 1],
            pressure: vec![
                C64::ONE,
                C64::ZERO,
                C64::from_re(sec),
                C64::new(0.0, -z * tan),
            ],
            flow: vec![
                C64::new(0.0, -tan / z),
                C64::from_re(-sec),
                C64::ZERO,
                C64::ONE,
            ],
            internal: vec![C64::ZERO; 4],
            dissipative: vec![C64::ZERO; 4],
        }
    }
}

#[test]
fn g1_recorded_mixed_sources_keep_all_cross_port_moments_and_work() {
    let drives = [
        [C64::new(1.2, -0.3), C64::new(0.2, 0.1)],
        [C64::new(-0.4, 0.9), C64::new(-0.1, 0.15)],
        [C64::new(0.7, 0.2), C64::new(0.05, -0.25)],
    ];
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0),
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = recorded(&drives, window);
        let owner = PreparedNetworkSources::recorded_test_mixer(
            core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap(),
        );
        let mut output = owner.statistics_buffer();
        let storage = output.pressure.values.as_ptr();
        owner
            .observe_welch_bin_into(
                &spectrum,
                BIN,
                owner.source_coordinates(),
                LIMITS,
                &mut output,
                || false,
            )
            .unwrap();
        let (mut pp, mut qq, mut pq) = ([C64::ZERO; 9], [C64::ZERO; 9], [C64::ZERO; 9]);
        let (mut load, mut loss, mut generator) = (0.0, 0.0, 0.0);
        for [x, q] in drives {
            // Independent resistor equations, not the stored transfer matrix.
            let p = (x + q.scale(2.0)).scale(3.0 / 5.0);
            let u = [(x - p).scale(0.5), q, -p.scale(1.0 / 3.0)];
            let weight = fraction / 6.0;
            for i in 0..3 {
                for j in 0..3 {
                    pp[3 * i + j] = pp[3 * i + j] + (p * p.conj()).scale(weight);
                    qq[3 * i + j] = qq[3 * i + j] + (u[i] * u[j].conj()).scale(weight);
                    pq[3 * i + j] = pq[3 * i + j] + (p * u[j].conj()).scale(weight);
                }
            }
            load += weight * p.norm_sq() / 3.0;
            loss += weight * 2.0 * u[0].norm_sq();
            generator += weight * (x * u[0].conj() + p * q.conj()).re;
        }
        for (actual, expected) in [
            (output.pressure_second_moments().unwrap(), pp),
            (output.flow_second_moments().unwrap(), qq),
            (output.pressure_flow_moments().unwrap(), pq),
        ] {
            for (a, b) in actual.iter().zip(expected) {
                close(a.re, b.re);
                close(a.im, b.im);
            }
        }
        let work = output.work().unwrap();
        close(work.source_power_w, load);
        close(work.external_load_power_w, load);
        close(work.network_power_w, 0.0);
        close(work.generator_power_w, generator);
        close(work.source_dissipation_w, loss);
        close(generator, load + loss);
        assert_eq!(storage, output.pressure.values.as_ptr());
    }
}

#[test]
fn g3_recorded_dark_pressure_keeps_source_heat_and_accepts_more_than_32_segments() {
    let x = C64::new(1.0, -0.5);
    let spectrum = recorded(&[[x, -x.scale(0.5)]; 40], WelchWindow::Rectangular);
    let owner = PreparedNetworkSources::recorded_test_mixer(
        core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap(),
    );
    let mut output = owner.statistics_buffer();
    owner
        .observe_welch_bin_into(
            &spectrum,
            BIN,
            owner.source_coordinates(),
            LIMITS,
            &mut output,
            || false,
        )
        .unwrap();
    for value in output.pressure_second_moments().unwrap() {
        assert!(value.abs() < 1e-28);
    }
    let work = output.work().unwrap();
    close(work.external_load_power_w, 0.0);
    close(work.source_power_w, 0.0);
    close(work.source_dissipation_w, x.norm_sq() / 4.0);
    close(work.generator_power_w, work.source_dissipation_w);
    assert!(
        owner
            .coherence_into(&vec![C64::ZERO; 80], 40, LIMITS, &mut output)
            .is_err()
    );
}

#[test]
fn g0_recorded_port_admission_and_cancellation_invalidate_previous_output() {
    let spectrum = recorded(
        &[[C64::ONE, C64::new(0.0, 0.2)]; 3],
        WelchWindow::Rectangular,
    );
    let owner = PreparedNetworkSources::recorded_test_mixer(
        core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap(),
    );
    let mut output = owner.statistics_buffer();
    owner
        .observe_welch_bin_into(
            &spectrum,
            BIN,
            owner.source_coordinates(),
            LIMITS,
            &mut output,
            || false,
        )
        .unwrap();
    let expected = output.work().unwrap();
    for bin in [0, N / 2, BIN + 1] {
        assert!(
            owner
                .observe_welch_bin_into(
                    &spectrum,
                    bin,
                    owner.source_coordinates(),
                    LIMITS,
                    &mut output,
                    || false
                )
                .is_err()
        );
        assert!(output.work().is_none());
    }
    let mut wrong = owner.source_coordinates().to_vec();
    wrong[0].quantity = NetworkSourceQuantity::VolumeFlow;
    assert!(
        owner
            .observe_welch_bin_into(&spectrum, BIN, &wrong, LIMITS, &mut output, || false)
            .is_err()
    );
    wrong = owner.source_coordinates().iter().rev().cloned().collect();
    assert!(
        owner
            .observe_welch_bin_into(&spectrum, BIN, &wrong, LIMITS, &mut output, || false)
            .is_err()
    );
    for budget in [
        NetworkCoherenceBudget {
            maximum_rank: 2,
            ..LIMITS
        },
        NetworkCoherenceBudget {
            maximum_complex_products: 218,
            ..LIMITS
        },
    ] {
        assert!(
            owner
                .observe_welch_bin_into(
                    &spectrum,
                    BIN,
                    owner.source_coordinates(),
                    budget,
                    &mut output,
                    || false
                )
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        owner.observe_welch_bin_into(
            &spectrum,
            BIN,
            owner.source_coordinates(),
            LIMITS,
            &mut output,
            || {
                polls += 1;
                polls == 3
            }
        ),
        Err(DuctError::Cancelled { .. })
    ));
    assert!(!output.is_valid());
    assert!(output.pressure_second_moments().is_none());
    owner
        .observe_welch_bin_into(
            &spectrum,
            BIN,
            owner.source_coordinates(),
            LIMITS,
            &mut output,
            || false,
        )
        .unwrap();
    assert_eq!(
        output.work().unwrap().generator_power_w.to_bits(),
        expected.generator_power_w.to_bits()
    );
}
