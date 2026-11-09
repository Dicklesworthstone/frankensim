use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_READOUT_QUERY_SCRATCH_BYTES, HarmonicReadoutWorkspace, HarmonicWindowProfile,
};

fn owner() -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega: 128.0,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![SourceId::StructuralMode("source".into())],
        receivers: vec![ReceiverSpec {
            id: "physical".into(),
            position_m: [0.2, 0.3, 0.4],
        }],
        producer_hash: ContentHash([19; 32]),
        radiation_momentum: true,
        transfers: vec![[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]],
    }
}

fn array(
    reference: &PreparedHarmonicObservation,
    statistics: HarmonicExposureFactorStatistics,
) -> PreparedHarmonicPressureArray<'_> {
    let lines = [1.0, -2.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(i, a)| PressureLine {
            omega: f64::from_bits(128.0_f64.to_bits() + i as u64),
            amplitudes: vec![
                [Dd::from_f64(a), Dd::ZERO],
                [Dd::from_f64(-0.5 * a), Dd::ZERO],
            ],
        })
        .collect();
    PreparedHarmonicPressureArray {
        reference,
        source_clock: "original-clock",
        channels: ["first", "second"]
            .map(|id| HarmonicPressureChannel {
                id: id.into(),
                taps: Vec::new(),
            })
            .into(),
        microphone_patterns: None,
        microphone_calibration: None,
        banks: vec![
            PressureBank {
                lines: Vec::new(),
                statistics: HarmonicExposureFactorStatistics::Real,
            },
            PressureBank { lines, statistics },
        ],
        workspace_bytes: 0,
        query_products: 0,
    }
}

#[test]
fn array_averaged_covariance_preserves_the_original_carrier_dark_readout() {
    let reference = owner();
    let array = array(&reference, HarmonicExposureFactorStatistics::ProperComplex);
    for (profile, mean) in [
        (
            HarmonicWindowProfile::Uniform,
            -2.940312874907294897526168082634475199e-30_f64,
        ),
        (
            HarmonicWindowProfile::Triangular,
            3.305785098798042449504765811355513971e-32_f64,
        ),
    ] {
        let window = HarmonicWindow::with_profile(-0.5, 1.0, profile).unwrap();
        let actual = array.observe_window(window).unwrap();
        for (&covariance, gain) in actual
            .averaged_covariance_pa2()
            .iter()
            .zip([1.0, -0.5, -0.5, 0.25])
        {
            let expected = 0.5 * gain * mean * mean;
            assert!(
                (covariance - expected).abs() < expected.abs() * 5e-13,
                "actual={covariance:e}, expected={expected:e}"
            );
        }
        assert!(actual.moments()[0].mean_square_pressure_pa2 > actual.averaged_covariance_pa2()[0]);
        assert_eq!(actual.moments()[0].mean_pressure_pa, 0.0);
        assert_eq!(actual.source_clock(), "original-clock");
        let mut polls = 0;
        array
            .observe_window_cancellable(window, || {
                polls += 1;
                false
            })
            .unwrap();
        for stop in [1, polls / 2, polls] {
            let mut calls = 0;
            assert!(matches!(
                array.observe_window_cancellable(window, || {
                    calls += 1;
                    calls >= stop
                }),
                Err(TransferError::Cancelled)
            ));
            assert_eq!(
                array
                    .observe_window(window)
                    .unwrap()
                    .averaged_covariance_pa2(),
                actual.averaged_covariance_pa2()
            );
        }
    }
}

#[test]
fn linear_array_projection_retains_both_components_independent_of_factor_statistics() {
    let reference = owner();
    for statistics in [
        HarmonicExposureFactorStatistics::Real,
        HarmonicExposureFactorStatistics::ProperComplex,
    ] {
        let mut array = array(&reference, statistics);
        for line in &mut array.banks[1].lines {
            for signal in &mut line.amplitudes {
                signal[1] = signal[0] * Dd::from_f64(0.25);
            }
        }
        let window = HarmonicWindow::new(-0.5, 1.0).unwrap();
        let mut means = [[Dd::ZERO; 2]; 2];
        assert!(readout::accumulate(&array.banks[1], window, &mut means, &mut || false).unwrap());
        for (value, gain) in means.into_iter().zip([1.0, -0.5]) {
            let real = gain * -2.940312874907294897526168082634475199e-30;
            assert!((value[0].to_f64() - real).abs() < real.abs() * 2e-13);
            assert!((value[1].to_f64() - 0.25 * real).abs() < real.abs() * 2e-13);
        }
    }
}

#[test]
fn readout_work_and_fixed_scratch_are_admitted_separately_from_power() {
    let (products, scratch) = readout::admission(3, 2).unwrap();
    assert!(products > 2 * 3 * 2);
    assert_eq!(
        scratch,
        core::mem::size_of::<HarmonicReadoutWorkspace>()
            + core::mem::size_of::<[f64; 64]>()
            + core::mem::size_of::<readout::Fields>()
            + HARMONIC_READOUT_QUERY_SCRATCH_BYTES
    );
    assert_eq!(readout::admission(1, 2).unwrap(), (4, 0));
    assert_eq!(readout::admission(3, 0).unwrap(), (0, 0));
    assert!(readout::admission(3, usize::MAX).is_err());
    let reference = owner();
    let array = array(&reference, HarmonicExposureFactorStatistics::ProperComplex);
    let mut means = [[Dd::ONE; 2]; 2];
    let original = means;
    assert!(
        !readout::accumulate(
            &array.banks[1],
            HarmonicWindow::new(0.0, 1e-3).unwrap(),
            &mut means,
            &mut || false
        )
        .unwrap()
    );
    assert_eq!(means, original);
}
