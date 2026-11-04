use super::*;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};

// Original gas, segment propagation and loaded graph assembly: this is an
// integration gate, distinct from the retained-map native operator tests.
#[test]
fn recorded_counterpropagating_waves_retain_native_loss_and_boundary_work() {
    const N: usize = 32;
    const BIN: usize = 3;
    let gas = gas();
    let length = 0.37;
    let net = cylinder_net(length);
    let omega = core::f64::consts::TAU * 384.0;
    let points = probes(length, &[0.0, 0.17, 0.6, 1.0]);
    let waves = [
        [C64::new(1.2, -0.3), C64::new(-0.4, 0.7)],
        [C64::new(0.2, 0.5), C64::new(0.8, -0.1)],
        [C64::new(-0.7, 0.1), C64::new(0.3, 0.4)],
    ];
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let wave = crate::segment_wave(&gas, 0.008, omega, loss).unwrap();
        let z = wave.characteristic_impedance.abs();
        let field = |x, pair: [C64; 2]| {
            analytic(
                wave.wavenumber,
                wave.characteristic_impedance,
                length,
                x,
                pair[0],
                pair[1],
            )
        };
        let sources = [
            pressure_source("right", C64::ZERO),
            pressure_source("left", C64::ZERO),
        ];
        let query = request(omega, loss, &sources, &[], &points);
        let owner = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
        let ports = owner.port_response();
        for (window, fraction) in [
            (WelchWindow::Rectangular, 1.0),
            (WelchWindow::PeriodicHann, 2.0 / 3.0),
        ] {
            let mut rows = Vec::new();
            for pair in waves {
                let endpoints = [field(0.0, pair)[0], field(length, pair)[0]];
                for j in 0..N {
                    let angle = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                    rows.push(endpoints.map(|p| p.re * angle.cos() + p.im * angle.sin()));
                }
            }
            let spectrum = WelchCrossSpectrum::estimate_interleaved(
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
                    maximum_samples_per_channel: 96,
                    maximum_transforms: 6,
                    maximum_workspace_bytes: 1_000_000,
                },
                || false,
            )
            .unwrap();
            let values = owner
                .observe_welch_bin(
                    &spectrum,
                    BIN,
                    ports.source_coordinates(),
                    NetworkCoherenceBudget::default(),
                    || false,
                )
                .unwrap();
            let mut boundary = ports.statistics_buffer();
            ports
                .observe_welch_bin_into(
                    &spectrum,
                    BIN,
                    ports.source_coordinates(),
                    NetworkCoherenceBudget::default(),
                    &mut boundary,
                    || false,
                )
                .unwrap();
            for (value, location) in values.iter().zip(owner.locations()) {
                let mut expected = [0.0; 2];
                let mut pq = C64::ZERO;
                for pair in waves {
                    let [p, q] = field(location.axial_position_m, pair);
                    expected[0] += p.norm_sq() * fraction / 6.0;
                    expected[1] += q.norm_sq() * fraction / 6.0;
                    pq = pq + (p * q.conj()).scale(fraction / 6.0);
                }
                assert!((value.mean_square_pressure_pa2 - expected[0]).abs() < 1e-9);
                assert!((value.mean_square_volume_flow_m6_s2 - expected[1]).abs() < 1e-9 / (z * z));
                assert_close(value.pressure_flow_moment_w, pq, 10.0 / z);
                assert!((value.mean_power_w - pq.re).abs() < 1e-9 / z);
                assert!(
                    (value.mean_intensity_w_m2 - pq.re / location.area_m2).abs()
                        < 1e-9 / (z * location.area_m2)
                );
            }
            let flux = boundary.pressure_flow_moments().unwrap();
            assert_close(flux[0], values[0].pressure_flow_moment_w, 10.0 / z);
            assert_close(flux[3], -values[3].pressure_flow_moment_w, 10.0 / z);
            let graph_work = values[0].mean_power_w - values[3].mean_power_w;
            assert!((boundary.work().unwrap().network_power_w - graph_work).abs() < 1e-9 / z);
            if loss == LossModel::Lossless {
                assert!(graph_work.abs() < 1e-9 / z);
            } else {
                assert!(graph_work > 0.0);
            }
        }
    }
}
