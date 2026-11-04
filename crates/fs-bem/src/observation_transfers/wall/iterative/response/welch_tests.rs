//! Original retained-impedance work under measured coherent source factors.
use super::*;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};

const N: usize = 32;
const BIN: usize = 3;
const DT: f64 = 1.0 / 4096.0;
const WORK: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 30,
};

fn spectrum(segments: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<_> = segments
        .iter()
        .flat_map(|values| {
            (0..N).map(move |j| {
                let theta = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                values.map(|a| a.re * theta.cos() + a.im * theta.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 3 * N,
            maximum_transforms: 6,
            maximum_workspace_bytes: 100_000,
        },
        || false,
    )
    .unwrap()
}

fn segments() -> [[C64; 2]; 3] {
    [
        [C64::new(0.003, -0.002), C64::new(2e-5, 1e-5)],
        [C64::new(-0.001, 0.004), C64::new(-1e-5, 2e-5)],
        [C64::new(0.002, 0.001), C64::new(3e-5, -1e-5)],
    ]
}

fn owner(spectrum: &WelchCrossSpectrum) -> BoundaryPortImpedance {
    // A raw retained-matrix fixture exercises the existing quadratic owner;
    // it is not a material law, boundary solve or passivity certificate.
    BoundaryPortImpedance {
        omega: core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap(),
        operator_hash: ContentHash([39; 32]),
        sources: vec![
            SourceId::StructuralMode("surface".into()),
            SourceId::ApertureVolumeFlow("opening".into()),
        ],
        impedance: vec![
            C64::new(800.0, 300.0),
            C64::new(10_000.0, 2_000.0),
            C64::new(8_000.0, -1_000.0),
            C64::new(4e6, 2e6),
        ],
        reciprocity_defect: 0.2,
        passivity_margin: None,
        projection_panel_products: 0,
    }
}

fn oracle(matrix: &[C64], rates: [C64; 2]) -> f64 {
    let mut power = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            power += 0.5 * (rates[i].conj() * matrix[2 * i + j] * rates[j]).re;
        }
    }
    power
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-15),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn measured_port_work_retains_mixed_coordinates_phase_and_hann_bin_units() {
    let inputs = segments();
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0),
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = spectrum(&inputs, window);
        let owner = owner(&spectrum);
        let original_matrix = owner.impedance().to_vec();
        let expected = inputs
            .iter()
            .map(|&q| oracle(owner.impedance(), q))
            .sum::<f64>()
            * fraction
            / 3.0;
        let mut scratch = [C64::ZERO; 2];
        let actual = owner
            .input_power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
                false
            })
            .unwrap();
        close(actual, expected);
        for rates in inputs {
            close(
                owner.input_power_w(&rates).unwrap(),
                oracle(owner.impedance(), rates),
            );
        }
        assert_eq!(owner.impedance(), original_matrix);
        assert_eq!(owner.operator_hash(), ContentHash([39; 32]));
        assert_eq!(owner.passivity_margin(), None);
    }
}

#[test]
fn measured_port_work_keeps_dark_interference_and_signed_raw_loading() {
    let q = C64::new(0.002, -0.001);
    let quiet_spectrum = spectrum(&[[q, -q]; 2], WelchWindow::Rectangular);
    let mut owner = owner(&quiet_spectrum);
    owner.impedance.fill(C64::from_re(17.0));
    let mut scratch = [C64::ZERO; 2];
    let quiet = owner
        .input_power_welch_bin(
            &quiet_spectrum,
            BIN,
            owner.sources(),
            WORK,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(quiet, 0.0);
    let partial = spectrum(&[[q, -q], [q, q]], WelchWindow::Rectangular);
    let actual = owner
        .input_power_welch_bin(&partial, BIN, owner.sources(), WORK, &mut scratch, || false)
        .unwrap();
    close(actual, 17.0 * q.norm_sq());

    // Raw negative work remains diagnosable; no positive clipping or repair.
    owner.impedance = vec![C64::from_re(-2.0), C64::ZERO, C64::ZERO, C64::ZERO];
    let actual = owner
        .input_power_welch_bin(&partial, BIN, owner.sources(), WORK, &mut scratch, || false)
        .unwrap();
    close(actual, -q.norm_sq());
}

#[test]
fn measured_port_work_admits_order_scratch_budget_and_cancellation_before_result() {
    let spectrum = spectrum(&segments(), WelchWindow::Rectangular);
    let owner = owner(&spectrum);
    let mut scratch = [C64::ZERO; 2];
    for bin in [0, N / 2, BIN + 1, usize::MAX] {
        assert!(
            owner
                .input_power_welch_bin(&spectrum, bin, owner.sources(), WORK, &mut scratch, || {
                    false
                })
                .is_err()
        );
    }
    let wrong_kind = [
        SourceId::ApertureVolumeFlow("surface".into()),
        owner.sources()[1].clone(),
    ];
    let reversed = [owner.sources()[1].clone(), owner.sources()[0].clone()];
    for order in [&wrong_kind[..], &reversed[..]] {
        assert!(
            owner
                .input_power_welch_bin(&spectrum, BIN, order, WORK, &mut scratch, || false)
                .is_err()
        );
    }
    assert!(
        owner
            .input_power_welch_bin(
                &spectrum,
                BIN,
                owner.sources(),
                WORK,
                &mut scratch[..1],
                || false
            )
            .is_err()
    );
    assert!(
        owner
            .input_power_welch_bin(
                &spectrum,
                BIN,
                owner.sources(),
                HarmonicWorkBudget {
                    maximum_products: 29
                },
                &mut scratch,
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.input_power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
            polls += 1;
            polls == 6
        }),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    let actual = owner
        .input_power_welch_bin(&spectrum, BIN, owner.sources(), WORK, &mut scratch, || {
            false
        })
        .unwrap();
    let expected = segments()
        .iter()
        .map(|&q| oracle(owner.impedance(), q))
        .sum::<f64>()
        / 3.0;
    close(actual, expected);
    assert!(owner.input_power_w(&[C64::ONE]).is_err());
    assert!(
        owner
            .input_power_w(&[C64::new(f64::NAN, 0.0), C64::ZERO])
            .is_err()
    );
}
