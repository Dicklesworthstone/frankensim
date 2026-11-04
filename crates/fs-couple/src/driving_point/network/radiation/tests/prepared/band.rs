//! Public original graph/BEM preparation and stationary selected-bin consumers.
use super::*;
use fs_bem::observation_transfers::PreparedRadiationOptions;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};

#[test]
fn sparse_loaded_band_reaches_all_original_observers_and_retains_propagation_identity() {
    const INDICES: [usize; 3] = [2, 5, 9];
    const DIRECTIONS: [[f64; 3]; 2] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let model = two_source_model();
    let gate = CancelGate::new_clock_free();
    let frequencies = INDICES.map(|b| core::f64::consts::TAU * 128.0 * b as f64);
    let drives = [
        [C64::new(1.0, 0.2), C64::new(0.3, -0.4)],
        [C64::new(-0.8, 0.1), C64::new(0.2, 0.5)],
        [C64::new(0.3, -0.7), C64::new(-0.4, 0.2)],
    ];
    let banks: Vec<_> = frequencies
        .iter()
        .map(|&omega| {
            let mut source = bank(omega);
            // The same reciprocal sources retain their authored resistance and
            // frequency-dependent inertance; no common scalar load is substituted.
            for z in &mut source.impedance.z {
                z.im *= omega / frequencies[0];
            }
            source
        })
        .collect();
    let owners: Vec<_> = frequencies
        .iter()
        .zip(&banks)
        .map(|(&omega, bank)| {
            model
                .prepare_sources_with_radiation(
                    omega,
                    LossModel::Bessel,
                    &[],
                    Some(bank),
                    &[],
                    &DIRECTIONS,
                    PreparedRadiationOptions {
                        angular: helmholtz::DirectivityProjectionOptions::default(),
                        initial_latitudes: 8,
                        refinements: 3,
                        maximum_kernel_evaluations: 4_000_000,
                        maximum_source_entries: 4096,
                    },
                    &gate,
                )
                .unwrap()
        })
        .collect();
    let rows: Vec<[f64; 2]> = (0..32)
        .map(|sample| {
            core::array::from_fn(|channel| {
                (0..3)
                    .map(|b| {
                        let phase =
                            core::f64::consts::TAU * INDICES[b] as f64 * sample as f64 / 32.0;
                        drives[b][channel].re * phase.cos() + drives[b][channel].im * phase.sin()
                    })
                    .sum()
            })
        })
        .collect();
    let spectrum = WelchCrossSpectrum::estimate_interleaved(
        &rows,
        1.0 / 4096.0,
        WelchConfig {
            segment_samples: 32,
            hop_samples: 32,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 32,
            maximum_transforms: 2,
            maximum_workspace_bytes: 100_000,
        },
        || false,
    )
    .unwrap();
    let bin_limits = RadiatingWelchLimits {
        maximum_segments: 1,
        maximum_mapping_products: 1000,
        maximum_factor_bytes: 1024,
    };
    let bins: Vec<_> = owners
        .iter()
        .zip(INDICES)
        .map(|(owner, index)| {
            owner
                .prepare_welch_bin(
                    &spectrum,
                    index,
                    owner.source_coordinates(),
                    bin_limits,
                    &gate,
                )
                .unwrap()
        })
        .collect();
    let limits = RadiatingWelchBandLimits {
        maximum_bins: 3,
        maximum_total_mapping_products: 3000,
        maximum_total_factor_bytes: 3072,
    };
    let band = PreparedRadiatingWelchBand::stationary(&bins, limits, &gate).unwrap();
    assert_eq!(band.coverage().selected_bins, 3);
    assert_eq!(band.coverage().omitted_interior_bins(), 12);
    assert!(!band.coverage().covers_all_interior_bins());
    let whole = HarmonicWorkBudget {
        maximum_products: 10_000_000,
    };
    let angular = HarmonicWorkBudget {
        maximum_products: 1_000_000,
    };
    let mut statistics = band.statistics(1_000_000, &gate).unwrap();
    statistics
        .evaluate(
            NetworkCoherenceBudget::default(),
            FIELD_BUDGET,
            whole,
            &gate,
        )
        .unwrap();
    let patterns = [
        FirstOrder::new(0.4, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.2, [0.0, 1.0, 0.0]).unwrap(),
    ];
    let mut microphones = band.microphones(&patterns, 1_000_000, &gate).unwrap();
    let microphone_moments = microphones.evaluate(FIELD_BUDGET, whole, &gate).unwrap();
    let mut radiation = band.far_field(1_000_000, &gate).unwrap();
    let power = radiation.evaluate(angular, whole, &gate).unwrap();
    let mut expected_pressure = [0.0; 2];
    let mut expected_microphones = [C64::ZERO; 4];
    let mut expected_work = [0.0; 3];
    let mut expected_power = 0.0;
    for (index, owner) in owners.iter().enumerate() {
        let mut response = owner.response_buffer();
        owner
            .drive_into(&drives[index], &mut response, &gate)
            .unwrap();
        let fields = response.receiver_fields().unwrap();
        let m: [C64; 2] = core::array::from_fn(|r| {
            patterns[r]
                .observe(
                    fields[r].pressure_pa,
                    fields[r].particle_velocity_m_s,
                    Medium::from(&model.gas()),
                )
                .unwrap()
        });
        for i in 0..2 {
            expected_pressure[i] += 0.5 * fields[i].pressure_pa.norm_sq();
            for j in 0..2 {
                expected_microphones[i * 2 + j] =
                    expected_microphones[i * 2 + j] + (m[i] * m[j].conj()).scale(0.5);
            }
        }
        let work = response.network().unwrap();
        expected_work[0] += work.generator_power_w;
        expected_work[1] += work.source_dissipation_w;
        expected_work[2] += response.exterior_surface_power_w().unwrap();
        expected_power += owner
            .far_field()
            .unwrap()
            .evaluate_coherent(&drives[index], angular, &mut [C64::ZERO; 2], &gate)
            .unwrap()
            .total_power_w;
    }
    for (field, expected) in statistics
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(expected_pressure)
    {
        real_close(field.mean_square_pressure_pa2, expected);
    }
    for (&actual, expected) in microphone_moments.iter().zip(expected_microphones) {
        close(actual, expected);
    }
    real_close(
        statistics.network_work().unwrap().generator_power_w,
        expected_work[0],
    );
    real_close(
        statistics.network_work().unwrap().source_dissipation_w,
        expected_work[1],
    );
    real_close(
        statistics.exterior_surface_power_w().unwrap(),
        expected_work[2],
    );
    real_close(power.total_power_w, expected_power);
    assert_eq!(radiation.bin_power_moments().unwrap().len(), 3);

    // Geometry and labels alone cannot certify the same propagation physics.
    let wrong = model
        .prepare_sources(
            frequencies[1],
            LossModel::Lossless,
            &[],
            Some(&banks[1]),
            &[],
            &gate,
        )
        .unwrap();
    let mismatched = [
        owners[0]
            .prepare_welch_bin(
                &spectrum,
                INDICES[0],
                owners[0].source_coordinates(),
                bin_limits,
                &gate,
            )
            .unwrap(),
        wrong
            .prepare_welch_bin(
                &spectrum,
                INDICES[1],
                wrong.source_coordinates(),
                bin_limits,
                &gate,
            )
            .unwrap(),
    ];
    assert!(PreparedRadiatingWelchBand::stationary(&mismatched, limits, &gate).is_err());
}
