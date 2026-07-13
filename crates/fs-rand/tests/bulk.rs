//! Bulk-generation bitwise-equivalence gates (bead frankensim-1za9, item 2).
//! The batched fill path MUST produce the exact same stream as sequential
//! draws (the fs-simd twin doctrine) and advance the index identically.

use fs_rand::{Stream, StreamCheckpoint, StreamKey};

fn stream() -> Stream {
    StreamKey {
        seed: 0xDEAD_BEEF_1234,
        kernel: 9,
        tile: 5,
    }
    .stream()
}

#[test]
fn fill_f64_is_bitwise_sequential() {
    // A non-multiple-of-8 length exercises the batched body AND the scalar tail.
    for len in [0usize, 1, 7, 8, 9, 16, 1000, 1003] {
        let mut bulk = vec![0.0f64; len];
        let mut a = stream();
        a.fill_f64(&mut bulk);

        let mut b = stream();
        for (i, want) in (0..len).map(|_| b.next_f64()).enumerate() {
            assert_eq!(bulk[i].to_bits(), want.to_bits(), "len {len} idx {i}");
        }
        // The index advanced by exactly `len`, so the streams stay in lockstep.
        assert_eq!(a.index(), b.index(), "index mismatch at len {len}");
    }
}

#[test]
fn fill_u64_is_bitwise_sequential() {
    for len in [0usize, 3, 8, 15, 500] {
        let mut bulk = vec![0u64; len];
        let mut a = stream();
        a.fill_u64(&mut bulk);

        let mut b = stream();
        for (i, item) in bulk.iter().enumerate() {
            assert_eq!(*item, b.next_u64(), "len {len} idx {i}");
        }
        assert_eq!(a.index(), b.index(), "index mismatch at len {len}");
    }
}

#[test]
fn bulk_is_deterministic() {
    let mut a = vec![0.0f64; 4096];
    let mut b = vec![0.0f64; 4096];
    stream().fill_f64(&mut a);
    stream().fill_f64(&mut b);
    assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
}

#[test]
fn counter_wraps_match_sequential_draws() {
    let key = StreamKey {
        seed: 0xDEAD_BEEF_1234,
        kernel: 9,
        tile: 5,
    };
    let checkpoint = StreamCheckpoint::current(key, u64::MAX - 3);
    let mut scalar = Stream::resume(checkpoint).expect("current checkpoint is supported");
    let want: Vec<u64> = (0..10).map(|_| scalar.next_u64()).collect();
    assert_eq!(scalar.index(), 6);

    let mut bulk = Stream::resume(checkpoint).expect("current checkpoint is supported");
    let mut got = vec![0u64; 10];
    bulk.fill_u64(&mut got);
    assert_eq!(got, want);
    assert_eq!(bulk.index(), scalar.index());
}
