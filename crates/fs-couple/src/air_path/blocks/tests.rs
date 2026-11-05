use super::*;
use crate::air_path::{ObserverFilter, process_fir};
use fs_material::iso9613::TerrestrialAir;

fn input(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| (i as f64 * 0.147).sin() + 0.2 * (i as f64 * 0.713).cos())
        .collect()
}

#[test]
fn g1_overlap_save_matches_independent_linear_convolution_and_retains_tail() {
    let taps: Vec<_> = (0..1537)
        .map(|i| (-0.004 * i as f64).exp() * (i as f64 * 0.13).cos() / 40.0)
        .collect();
    let mut signal = input(1703);
    signal.resize(4000, 0.0);
    let expected: Vec<f64> = (0..signal.len())
        .map(|i| {
            (0..taps.len().min(i + 1))
                .map(|j| taps[j] * signal[i - j])
                .sum()
        })
        .collect();
    let mut blocks = FirBlocks::prepare(&taps, 512).unwrap();
    let mut history = vec![0.0; taps.len()];
    let mut head = 0;
    let mut actual = signal.clone();
    // Exercise FFT blocks, small direct blocks, non-power-of-two boundaries,
    // scalar interleaving and a final partial block on one unbroken history.
    let sizes = [512, 1, 31, 509, 128];
    let mut start = 0;
    let mut step = 0;
    while start < actual.len() {
        let count = sizes[step % sizes.len()].min(actual.len() - start);
        if count == 1 {
            actual[start] = process_fir(&taps, &mut history, &mut head, actual[start]);
        } else {
            blocks
                .process(
                    &taps,
                    &mut history,
                    &mut head,
                    &mut actual[start..start + count],
                )
                .unwrap();
        }
        start += count;
        step += 1;
    }
    for (i, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
        assert!((a - b).abs() < 2e-12, "sample {i}: {a} != {b}");
    }
    assert!(actual[1703..].iter().any(|v| v.abs() > 1e-4));
    assert!(actual[3240..].iter().all(|v| v.abs() < 2e-12));
}

#[test]
fn g1_delayed_impulse_has_no_block_delay_or_circular_wrap() {
    let mut taps = vec![0.0; 2049];
    taps[2048] = 1.0;
    let mut blocks = FirBlocks::prepare(&taps, 512).unwrap();
    let mut history = vec![0.0; taps.len()];
    let mut head = 0;
    let mut pressure = vec![0.0; 3584];
    pressure[0] = 1.0;
    pressure[1023] = -0.5;
    for block in pressure.chunks_mut(512) {
        blocks
            .process(&taps, &mut history, &mut head, block)
            .unwrap();
    }
    for (i, &value) in pressure.iter().enumerate() {
        let expected = match i {
            2048 => 1.0,
            3071 => -0.5,
            _ => 0.0,
        };
        assert!((value - expected).abs() < 2e-12, "sample {i}: {value}");
    }
}

#[test]
fn g0_refused_blocks_preserve_input_and_history_then_recover() {
    let mut taps = vec![0.0; 513];
    taps[0] = 2.0;
    let mut blocks = FirBlocks::prepare(&taps, 512).unwrap();
    let mut history = vec![0.25; taps.len()];
    let before = history.clone();
    let mut head = 75;
    for mut bad in [
        vec![f64::NAN],
        vec![f64::INFINITY],
        vec![1.0; 513],
        vec![f64::MAX; 512],
    ] {
        let original: Vec<_> = bad.iter().map(|x| x.to_bits()).collect();
        assert!(
            blocks
                .process(&taps, &mut history, &mut head, &mut bad)
                .is_err()
        );
        assert_eq!(
            original,
            bad.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(history, before);
        assert_eq!(head, 75);
    }
    let mut good = vec![0.125; 512];
    blocks
        .process(&taps, &mut history, &mut head, &mut good)
        .unwrap();
    assert!(good.iter().all(|x| (*x - 0.25).abs() < 1e-12));
}

#[test]
fn g3_original_air_response_and_receipt_survive_block_preparation_and_clone() {
    let (gas, air) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
    let filters = [
        ObserverFilter::iso9613(&air, &gas, 0.5, 30.0, 48_000).unwrap(),
        ObserverFilter::stokes_kirchhoff(&gas, 30.0, 48_000).unwrap(),
        ObserverFilter::lossless(&gas, 30.0, 48_000).unwrap(),
    ];
    for mut actual in filters {
        let receipt = format!("{:?}", actual.receipt());
        let taps = actual.taps().to_vec();
        let mut reference = actual.clone();
        assert!(actual.process_block(&mut [1.0]).is_err());
        for x in input(75) {
            assert_eq!(actual.process(x), reference.process(x));
        }
        actual.prepare_blocks(512).unwrap();
        assert_eq!(actual.block_capacity(), Some(512));
        assert!(actual.prepare_blocks(0).is_err());
        assert!(
            actual
                .prepare_blocks(MAX_OBSERVER_BLOCK_SAMPLES + 1)
                .is_err()
        );
        assert_eq!(actual.block_capacity(), Some(512));
        for n in [512, 113, 1, 511] {
            let mut pressure = input(n);
            let expected: Vec<_> = pressure.iter().map(|&x| reference.process(x)).collect();
            actual.process_block(&mut pressure).unwrap();
            for (&a, &b) in pressure.iter().zip(&expected) {
                assert!((a - b).abs() < 3e-12);
            }
        }
        let mut resumed = actual.clone();
        let mut a = input(512);
        let mut b = a.clone();
        actual.process_block(&mut a).unwrap();
        resumed.process_block(&mut b).unwrap();
        assert_eq!(a, b);
        actual.prepare_blocks(1024).unwrap();
        let mut a = input(512);
        let expected: Vec<_> = a.iter().map(|&x| resumed.process(x)).collect();
        actual.process_block(&mut a).unwrap();
        assert!(a.iter().zip(expected).all(|(&x, y)| (x - y).abs() < 3e-12));
        assert_eq!(format!("{:?}", actual.receipt()), receipt);
        assert_eq!(actual.taps(), taps);
    }
}
