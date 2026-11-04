use super::*;
use fs_bem::observation_transfers::PreparedRadiationOptions;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};

const N: usize = 32;
const BIN: usize = 3;
const OMEGA: f64 = core::f64::consts::TAU * 1200.0;
const LIMITS: RadiatingWelchLimits = RadiatingWelchLimits {
    maximum_segments: 64,
    maximum_mapping_products: 16_384,
    maximum_factor_bytes: 32_768,
};
const NETWORK: NetworkCoherenceBudget = NetworkCoherenceBudget {
    maximum_rank: 64,
    maximum_complex_products: 262_144,
};
const ANGULAR: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 1_000_000,
};

fn recording(drives: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<_> = drives
        .iter()
        .flat_map(|drive| {
            (0..N).map(move |j| {
                let phase = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                drive.map(|z| z.re * phase.cos() + z.im * phase.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        1.0 / 12_800.0,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: rows.len(),
            maximum_transforms: 2 * drives.len(),
            maximum_workspace_bytes: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn recorded_mixed_generators_reach_original_loaded_fields_microphones_and_radiation() {
    let mut model = two_source_model();
    // Distinct source attachments produce different opening and angular phases.
    model.driving.network.edges.last_mut().unwrap().to = 2;
    let sources = [
        PortSource {
            label: "inlet".into(),
            excitation: PortExcitation::VolumeFlow(C64::ZERO),
        },
        PortSource {
            label: "aux".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::ZERO,
                impedance_pa_s_m3: C64::new(2e6, -1e6),
            },
        },
    ];
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources_with_radiation(
            OMEGA,
            LossModel::Bessel,
            &sources,
            None,
            &[],
            &[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            PreparedRadiationOptions {
                angular: helmholtz::DirectivityProjectionOptions::default(),
                initial_latitudes: 8,
                refinements: 3,
                maximum_kernel_evaluations: 4_000_000,
                maximum_source_entries: 4096,
            },
            &gate,
        )
        .unwrap();
    assert_eq!(
        prepared.source_coordinates()[0].quantity,
        NetworkSourceQuantity::OpenCircuitPressure
    );
    assert_eq!(
        prepared.source_coordinates()[1].quantity,
        NetworkSourceQuantity::VolumeFlow
    );
    let drives = [
        [C64::new(0.8, -0.2), C64::new(0.4e-6, 0.2e-6)],
        [C64::new(-0.3, 0.7), C64::new(0.1e-6, -0.5e-6)],
        [C64::new(0.1, 0.4), C64::new(-0.3e-6, 0.1e-6)],
    ];
    let patterns = [
        FirstOrder::new(0.4, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.2, [0.0, 1.0, 0.0]).unwrap(),
    ];
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0_f64),
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = recording(&drives, window);
        let measured = prepared
            .prepare_welch_bin(&spectrum, BIN, prepared.source_coordinates(), LIMITS, &gate)
            .unwrap();
        let factors: Vec<_> = drives
            .iter()
            .flatten()
            .map(|z| z.scale((fraction / 3.0).sqrt()))
            .collect();
        let mut actual = prepared.statistics_buffer();
        let mut reference = prepared.statistics_buffer();
        measured
            .evaluate_into(NETWORK, FIELD_BUDGET, &mut actual, &gate)
            .unwrap();
        prepared
            .coherence_into(&factors, 3, NETWORK, FIELD_BUDGET, &mut reference, &gate)
            .unwrap();
        for (a, b) in actual
            .receiver_fields()
            .unwrap()
            .iter()
            .zip(reference.receiver_fields().unwrap())
        {
            real_close(a.mean_square_pressure_pa2, b.mean_square_pressure_pa2);
            for i in 0..3 {
                close(a.complex_intensity_w_m2[i], b.complex_intensity_w_m2[i]);
                for j in 0..3 {
                    close(a.velocity_cross_m2_s2[i][j], b.velocity_cross_m2_s2[i][j]);
                }
            }
        }
        for (a, b) in actual
            .network()
            .unwrap()
            .pressure_flow_moments()
            .unwrap()
            .iter()
            .zip(
                reference
                    .network()
                    .unwrap()
                    .pressure_flow_moments()
                    .unwrap(),
            )
        {
            close(*a, *b);
        }
        real_close(
            actual.exterior_surface_power_w().unwrap(),
            reference.exterior_surface_power_w().unwrap(),
        );
        let work = actual.network().unwrap().work().unwrap();
        real_close(
            work.generator_power_w,
            work.source_dissipation_w + work.network_power_w + work.external_load_power_w,
        );
        real_close(
            actual.exterior_surface_power_w().unwrap(),
            work.external_load_power_w,
        );
        let mut microphones = prepared.microphones(&patterns).unwrap();
        let expected = microphones
            .coherence(&factors, 3, FIELD_BUDGET, &gate)
            .unwrap()
            .to_vec();
        for (&a, b) in microphones
            .welch_bin(&measured, FIELD_BUDGET, &gate)
            .unwrap()
            .iter()
            .zip(expected)
        {
            close(a, b);
        }
        let far = prepared.far_field().unwrap();
        let mut density = [0.0; 2];
        let mut expected_density = [0.0; 2];
        let power = far
            .evaluate_welch_bin(&measured, ANGULAR, &mut density, &gate)
            .unwrap();
        let expected_power = far
            .evaluate_coherence_modes(&factors, 3, ANGULAR, &mut expected_density, &gate)
            .unwrap();
        real_close(power.total_power_w, expected_power.total_power_w);
        for (a, b) in density.into_iter().zip(expected_density) {
            real_close(a, b);
        }
        for i in 0..3 {
            close(
                C64::from_re(power.first_directional_moment_w[i]),
                C64::from_re(expected_power.first_directional_moment_w[i]),
            );
        }
    }
}

#[test]
fn forty_recorded_antisymmetric_generators_keep_heat_and_transactional_refusals() {
    let model = two_source_model();
    let bank = bank(OMEGA);
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources(OMEGA, LossModel::Lossless, &[], Some(&bank), &[], &gate)
        .unwrap();
    let spectrum = recording(&[[C64::ONE, -C64::ONE]; 40], WelchWindow::Rectangular);
    let measured = prepared
        .prepare_welch_bin(&spectrum, BIN, prepared.source_coordinates(), LIMITS, &gate)
        .unwrap();
    assert_eq!(measured.segment_count(), 40);
    let mut response = prepared.response_buffer();
    prepared
        .drive_into(&[C64::ONE, -C64::ONE], &mut response, &gate)
        .unwrap();
    let mut statistics = prepared.statistics_buffer();
    measured
        .evaluate_into(NETWORK, FIELD_BUDGET, &mut statistics, &gate)
        .unwrap();
    let work = statistics.network().unwrap().work().unwrap();
    real_close(
        work.source_dissipation_w,
        response.network().unwrap().source_dissipation_w,
    );
    assert!(work.source_dissipation_w > 0.0);
    for field in statistics.receiver_fields().unwrap() {
        assert!(field.mean_square_pressure_pa2 < 1e-20);
    }
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(
        measured
            .evaluate_into(NETWORK, FIELD_BUDGET, &mut statistics, &stopped)
            .is_err()
    );
    assert!(!statistics.is_valid());
    assert!(statistics.network().is_none() && statistics.receiver_fields().is_none());
    measured
        .evaluate_into(NETWORK, FIELD_BUDGET, &mut statistics, &gate)
        .unwrap();
    assert!(
        measured
            .evaluate_into(
                NetworkCoherenceBudget::default(),
                FIELD_BUDGET,
                &mut statistics,
                &gate
            )
            .is_err()
    );
    assert!(!statistics.is_valid());
    let mut wrong_units = prepared.source_coordinates().to_vec();
    wrong_units[0].quantity = NetworkSourceQuantity::VolumeFlow;
    assert!(
        prepared
            .prepare_welch_bin(&spectrum, BIN, &wrong_units, LIMITS, &gate)
            .is_err()
    );
    assert!(
        prepared
            .prepare_welch_bin(&spectrum, 0, prepared.source_coordinates(), LIMITS, &gate)
            .is_err()
    );
    assert!(
        prepared
            .prepare_welch_bin(
                &spectrum,
                BIN,
                prepared.source_coordinates(),
                RadiatingWelchLimits {
                    maximum_segments: 39,
                    ..LIMITS
                },
                &gate
            )
            .is_err()
    );
}
