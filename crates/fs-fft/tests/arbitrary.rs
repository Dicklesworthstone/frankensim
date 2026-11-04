//! G0/G1/G3: exact physical record lengths, independent DFT and replay.

use fs_fft::{C64, Fft, GeneralFft, GeneralFftBudget, GeneralFftError};

fn budget() -> GeneralFftBudget {
    GeneralFftBudget {
        max_samples: 65_536,
        max_workspace_bytes: 32 * 1024 * 1024,
    }
}

fn record(n: usize) -> Vec<C64> {
    let mut state = 0xa1b2_c3d4_e5f6_0123_u64;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    };
    (0..n).map(|_| C64::new(next(), next())).collect()
}

// Independent direct summation with platform trigonometry, not the chirp
// factorization or production strict twiddle table.
fn direct(input: &[C64], inverse: bool) -> Vec<C64> {
    let n = input.len();
    let sign = if inverse { 1.0 } else { -1.0 };
    (0..n)
        .map(|k| {
            let mut re = 0.0;
            let mut im = 0.0;
            for (j, z) in input.iter().enumerate() {
                let angle = sign * core::f64::consts::TAU * (j * k) as f64 / n as f64;
                let (sin, cos) = angle.sin_cos();
                re += z.re * cos - z.im * sin;
                im += z.re * sin + z.im * cos;
            }
            let scale = if inverse { 1.0 / n as f64 } else { 1.0 };
            C64::new(re * scale, im * scale)
        })
        .collect()
}

fn close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (k, (a, b)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a.re - b.re).abs().max((a.im - b.im).abs()) <= tolerance,
            "bin {k}: {a:?}, expected {b:?}, tolerance {tolerance:e}"
        );
    }
}

fn bits(values: &[C64]) -> Vec<[u64; 2]> {
    values
        .iter()
        .map(|z| [z.re.to_bits(), z.im.to_bits()])
        .collect()
}

#[test]
fn every_small_length_matches_independent_forward_and_inverse_dft() {
    for n in (1..=65).chain([97, 127, 257]) {
        let plan = GeneralFft::new(n, budget()).unwrap();
        assert_eq!(plan.len(), n);
        assert!(!plan.is_empty());
        let input = record(n);
        let sentinel = C64::new(137.0, -29.0);
        let mut scratch = vec![sentinel; plan.scratch_len() + 3];
        let mut values = input.clone();
        plan.forward(&mut values, &mut scratch).unwrap();
        close(&values, &direct(&input, false), 2e-12 * n as f64);
        assert_eq!(&scratch[plan.scratch_len()..], &[sentinel; 3]);
        plan.inverse(&mut values, &mut scratch).unwrap();
        close(&values, &input, 2e-13);
        values.copy_from_slice(&input);
        plan.inverse(&mut values, &mut scratch).unwrap();
        close(&values, &direct(&input, true), 2e-12);
    }
}

#[test]
fn odd_record_retains_its_original_grid_and_highest_positive_frequency() {
    let n = 1001;
    let timestep_s = 0.0007;
    let plan = GeneralFft::new(n, budget()).unwrap();
    let mut scratch = vec![C64::default(); plan.scratch_len()];
    for bin in [29, n / 2] {
        let physical_frequency_hz = bin as f64 / (n as f64 * timestep_s);
        assert!(physical_frequency_hz < 0.5 / timestep_s);
        let mut values: Vec<_> = (0..n)
            .map(|j| {
                let phase = core::f64::consts::TAU * physical_frequency_hz * j as f64 * timestep_s;
                C64::new(1.2 * phase.cos() - 0.7 * phase.sin(), 0.0)
            })
            .collect();
        plan.forward(&mut values, &mut scratch).unwrap();
        close(
            &values[bin..=bin],
            &[C64::new(0.6 * n as f64, 0.35 * n as f64)],
            2e-10,
        );
        close(&values[n - bin..=n - bin], &[values[bin].conj()], 2e-10);
        for (k, z) in values.iter().enumerate() {
            if k != bin && k != n - bin {
                assert!(
                    z.re.abs().max(z.im.abs()) < 2e-10,
                    "spurious bin {k}: {z:?}"
                );
            }
        }
    }
}

#[test]
fn prime_length_parseval_roundtrip_and_replay() {
    let n = 4093;
    let input = record(n);
    let plan = GeneralFft::new(n, budget()).unwrap();
    let mut values = input.clone();
    let mut scratch = vec![C64::default(); plan.scratch_len()];
    plan.forward(&mut values, &mut scratch).unwrap();
    let time_energy: f64 = input.iter().map(|z| z.norm_sq()).sum();
    let spectral_energy = values.iter().map(|z| z.norm_sq()).sum::<f64>() / n as f64;
    assert!((spectral_energy / time_energy - 1.0).abs() < 3e-14);
    let mut repeated = input.clone();
    plan.forward(&mut repeated, &mut scratch).unwrap();
    assert_eq!(bits(&repeated), bits(&values));
    plan.inverse(&mut values, &mut scratch).unwrap();
    close(&values, &input, 3e-14);
}

#[test]
fn power_of_two_route_preserves_original_transform_bits() {
    let sizes = [1, 2, 4, 8, 16, 64, 256, 2048];
    for n in sizes
        .into_iter()
        .chain(if cfg!(feature = "frontier-sixstep") {
            Some(65_536)
        } else {
            None
        })
    {
        let original = Fft::new(n);
        let plan = GeneralFft::new(n, budget()).unwrap();
        assert_eq!(plan.scratch_len(), n);
        let mut first = record(n);
        let mut second = first.clone();
        let mut a = vec![C64::default(); n];
        let mut b = a.clone();
        original.forward(&mut first, &mut a);
        plan.forward(&mut second, &mut b).unwrap();
        assert_eq!(bits(&first), bits(&second));
        original.inverse(&mut first, &mut a);
        plan.inverse(&mut second, &mut b).unwrap();
        assert_eq!(bits(&first), bits(&second));
    }
}

#[test]
fn physical_amplitude_scaling_survives_convolution_intermediates() {
    let n = 31;
    let plan = GeneralFft::new(n, budget()).unwrap();
    let mut scratch = vec![C64::default(); plan.scratch_len()];
    // An impulse has the same amplitude in every frequency bin. Its large
    // version has a finite DFT even when unscaled convolution intermediates
    // overflow. The small version exercises subnormal physical output.
    for scale in [1e-310, 1e-200, 1.0, 1e200, f64::MAX / 4.0] {
        let mut values = vec![C64::default(); n];
        values[0] = C64::new(scale, 0.0);
        plan.forward(&mut values, &mut scratch).unwrap();
        for value in &values {
            assert!((value.re / scale - 1.0).abs() < 2e-13);
            assert!((value.im / scale).abs() < 2e-13);
        }
        plan.inverse(&mut values, &mut scratch).unwrap();
        assert!((values[0].re / scale - 1.0).abs() < 2e-12);
        for value in values.iter().skip(1) {
            assert!((value.re / scale).abs().max((value.im / scale).abs()) < 2e-12);
        }
    }
}

#[test]
fn explicit_size_and_payload_limits_refuse_before_preparation() {
    assert!(matches!(
        GeneralFft::new(0, budget()),
        Err(GeneralFftError::InvalidSize { .. })
    ));
    assert!(GeneralFft::required_workspace_bytes(usize::MAX).is_err());
    assert!(matches!(
        GeneralFft::new(
            13,
            GeneralFftBudget {
                max_samples: 12,
                ..budget()
            }
        ),
        Err(GeneralFftError::SampleLimit {
            samples: 13,
            maximum: 12
        })
    ));
    let bytes = GeneralFft::required_workspace_bytes(13).unwrap();
    assert_eq!(bytes, (13 + 32 * 3 + 32 / 2) * core::mem::size_of::<C64>());
    assert!(matches!(
        GeneralFft::new(13, GeneralFftBudget { max_workspace_bytes: bytes - 1, ..budget() }),
        Err(GeneralFftError::WorkspaceBudget { required_bytes, .. }) if required_bytes == bytes
    ));
    let plan = GeneralFft::new(
        13,
        GeneralFftBudget {
            max_workspace_bytes: bytes,
            ..budget()
        },
    )
    .unwrap();
    assert_eq!(plan.workspace_bytes(), bytes);
    assert_eq!(plan.scratch_len(), 64);
}

#[test]
fn malformed_buffers_nonfinite_inputs_and_overflow_are_reported() {
    for n in [3, 4] {
        let plan = GeneralFft::new(n, budget()).unwrap();
        let input = record(n);
        let mut values = input.clone();
        let mut short = vec![C64::new(7.0, 9.0); plan.scratch_len() - 1];
        let original_scratch = short.clone();
        assert!(matches!(
            plan.forward(&mut values, &mut short),
            Err(GeneralFftError::BufferLength { .. })
        ));
        assert_eq!(values, input);
        assert_eq!(short, original_scratch);
        let mut scratch = vec![C64::new(13.0, 11.0); plan.scratch_len()];
        assert!(matches!(
            plan.forward(&mut values[..n - 1], &mut scratch),
            Err(GeneralFftError::BufferLength { .. })
        ));
        assert_eq!(values, input);
        values[1].im = f64::NAN;
        let before = bits(&values);
        let scratch_before = scratch.clone();
        assert_eq!(
            plan.inverse(&mut values, &mut scratch),
            Err(GeneralFftError::NonFinite { what: "input" })
        );
        assert_eq!(bits(&values), before);
        assert_eq!(scratch, scratch_before);
        values.fill(C64::new(f64::MAX, 0.0));
        assert_eq!(
            plan.forward(&mut values, &mut scratch),
            Err(GeneralFftError::NonFinite { what: "output" })
        );
        values.fill(C64::default());
        plan.forward(&mut values, &mut scratch).unwrap();
        assert!(values.iter().all(|z| z.re == 0.0 && z.im == 0.0));
    }
}
