use super::*;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};

const N: usize = 32;
const BIN: usize = 3;
const ANGLE: f64 = 0.63;
const AREA: f64 = 0.01;
const Z: f64 = 1.2 * 340.0 / AREA;
const BUDGET: NetworkCoherenceBudget = NetworkCoherenceBudget {
    maximum_rank: 64,
    maximum_complex_products: 10_000,
};

fn spectrum(drives: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<_> = drives
        .iter()
        .flat_map(|a| {
            (0..N).map(move |j| {
                let t = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                a.map(|v| v.re * t.cos() + v.im * t.sin())
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

fn fields(omega: f64) -> PreparedNetworkFields {
    let mut pressure = Vec::new();
    let mut flow = Vec::new();
    let mut locations = Vec::new();
    for phase in [0.0, ANGLE * 0.5, ANGLE] {
        // Exact analytical uniform-tube boundary fields: p(0)=x, U(L)=-q.
        // Retained-map fixture only; graph/material assembly is not substituted.
        pressure.extend([
            C64::from_re((ANGLE - phase).cos() / ANGLE.cos()),
            C64::new(0.0, -Z * phase.sin() / ANGLE.cos()),
        ]);
        flow.extend([
            C64::new(0.0, -(ANGLE - phase).sin() / (Z * ANGLE.cos())),
            C64::from_re(-phase.cos() / ANGLE.cos()),
        ]);
        let x = phase * 340.0 / omega;
        locations.push(NetworkFieldLocation {
            probe: NetworkFieldProbe {
                edge_id: "tube".into(),
                section: NetworkFieldSection::Cylinder {
                    segment: 0,
                    offset_m: x,
                },
            },
            axial_position_m: x,
            area_m2: AREA,
            density_kg_m3: 1.2,
            sound_speed_m_s: 340.0,
        });
    }
    PreparedNetworkFields {
        ports: PreparedNetworkSources::recorded_test_tube(omega, ANGLE, Z),
        source_clock: "recorded source clock".into(),
        locations,
        pressure,
        flow,
    }
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 4e-12 * scale, "{a:e} != {b:e}");
}

#[test]
fn g1_recorded_native_fields_keep_travel_direction_partial_coherence_and_hann_units() {
    let x = C64::new(1.0, 0.3);
    for directions in [[1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]] {
        let drives = directions.map(|direction| {
            let phase = direction * ANGLE;
            [
                x,
                (x * C64::new(phase.cos(), phase.sin())).scale(-direction / Z),
            ]
        });
        for (window, fraction) in [
            (WelchWindow::Rectangular, 1.0),
            (WelchWindow::PeriodicHann, 2.0 / 3.0),
        ] {
            let s = spectrum(&drives, window);
            let owner = fields(core::f64::consts::TAU * s.frequency_hz(BIN).unwrap());
            let values = owner
                .observe_welch_bin(&s, BIN, owner.ports.source_coordinates(), BUDGET, || false)
                .unwrap();
            let pp = 0.5 * x.norm_sq() * fraction;
            let power = pp * (directions[0] + directions[1]) / (2.0 * Z);
            for value in values {
                close(value.mean_square_pressure_pa2, pp, pp);
                close(
                    value.mean_square_volume_flow_m6_s2,
                    pp / (Z * Z),
                    pp / (Z * Z),
                );
                close(
                    value.mean_square_particle_velocity_m2_s2,
                    pp / (Z * Z * AREA * AREA),
                    pp / (Z * Z * AREA * AREA),
                );
                close(value.mean_power_w, power, pp / Z);
                close(value.mean_intensity_w_m2, power / AREA, pp / (Z * AREA));
                close(value.pressure_flow_moment_w.im, 0.0, pp / Z);
            }
        }
    }
}

#[test]
fn g1_recorded_closed_tube_keeps_reactive_standing_wave_and_zero_antinode_flux() {
    let x = C64::new(1.0, 0.2);
    let s = spectrum(&[[x, C64::ZERO]; 2], WelchWindow::Rectangular);
    let owner = fields(core::f64::consts::TAU * s.frequency_hz(BIN).unwrap());
    let values = owner
        .observe_welch_bin(&s, BIN, owner.ports.source_coordinates(), BUDGET, || false)
        .unwrap();
    for (value, phase) in values.iter().zip([0.0, ANGLE * 0.5, ANGLE]) {
        let p = (ANGLE - phase).cos() / ANGLE.cos();
        let u = (ANGLE - phase).sin() / (Z * ANGLE.cos());
        let square = 0.5 * x.norm_sq();
        close(value.mean_square_pressure_pa2, square * p * p, square);
        close(
            value.mean_square_volume_flow_m6_s2,
            square * u * u,
            square / (Z * Z),
        );
        close(value.mean_power_w, 0.0, square / Z);
        close(value.pressure_flow_moment_w.im, square * p * u, square / Z);
    }
    assert!(values[2].mean_square_pressure_pa2 > 0.5);
    assert_eq!(values[2].mean_square_volume_flow_m6_s2, 0.0);
    assert_eq!(values[2].mean_intensity_w_m2, 0.0);
    let mut ports = owner.ports.statistics_buffer();
    owner
        .ports
        .observe_welch_bin_into(
            &s,
            BIN,
            owner.ports.source_coordinates(),
            BUDGET,
            &mut ports,
            || false,
        )
        .unwrap();
    close(ports.work().unwrap().network_power_w, 0.0, 1.0 / Z);
}

#[test]
fn g0_recorded_spatial_admission_cancellation_and_replay_are_bounded() {
    let s = spectrum(
        &[[C64::ONE, C64::from_re(1.0 / Z)]; 2],
        WelchWindow::Rectangular,
    );
    let owner = fields(core::f64::consts::TAU * s.frequency_hz(BIN).unwrap());
    let order = owner.ports.source_coordinates();
    let expected = owner
        .observe_welch_bin(&s, BIN, order, BUDGET, || false)
        .unwrap();
    for bin in [0, N / 2, BIN + 1] {
        assert!(
            owner
                .observe_welch_bin(&s, bin, order, BUDGET, || false)
                .is_err()
        );
    }
    let mut wrong = order.to_vec();
    wrong[0].quantity = NetworkSourceQuantity::VolumeFlow;
    assert!(
        owner
            .observe_welch_bin(&s, BIN, &wrong, BUDGET, || false)
            .is_err()
    );
    // 2 factors*(2*3 probes*2 sources+3*3+2*2)=50 products.
    assert!(
        owner
            .observe_welch_bin(
                &s,
                BIN,
                order,
                NetworkCoherenceBudget {
                    maximum_complex_products: 49,
                    ..BUDGET
                },
                || false
            )
            .is_err()
    );
    assert!(
        owner
            .observe_welch_bin(
                &s,
                BIN,
                order,
                NetworkCoherenceBudget {
                    maximum_rank: 1,
                    ..BUDGET
                },
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.observe_welch_bin(&s, BIN, order, BUDGET, || {
            polls += 1;
            polls == 5
        }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(
        owner
            .observe_welch_bin(&s, BIN, order, BUDGET, || false)
            .unwrap(),
        expected
    );
}
