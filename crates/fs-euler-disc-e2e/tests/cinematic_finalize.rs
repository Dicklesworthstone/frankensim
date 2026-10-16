//! Public-boundary tests for independent cinematic bundle finalization.
//!
//! The complete valid/hostile bundle matrix lives beside the verifier, where
//! it can reuse the deliberately small image/audio fixture.  This target
//! instead checks that an external producer can persist the public receipt
//! wire format deterministically and receives stable refusals at explicit
//! budget and cancellation boundaries.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_euler_disc_e2e::cinematic_finalize::{
    CinematicReceiptError, encode_audio_video_alignment_receipt,
};
use fs_euler_disc_e2e::{AudioVideoAlignment, AudioVideoSyncMarker};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5055_424c_4943_5f46,
                kernel_id: 0x494e_414c,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn take<const N: usize>(bytes: &[u8], cursor: &mut usize) -> [u8; N] {
    let end = cursor.checked_add(N).expect("wire offset must fit usize");
    let value = bytes
        .get(*cursor..end)
        .expect("canonical receipt must contain the declared field");
    *cursor = end;
    value
        .try_into()
        .expect("slice has the requested exact width")
}

/// Decode from the documented canonical wire independently of the producer.
///
/// This deliberately tiny test oracle catches field-order, width, endianness,
/// count, and trailing-byte drift at the public persistence boundary.
fn decode_alignment_wire(bytes: &[u8]) -> AudioVideoAlignment {
    let mut cursor = 0;
    assert_eq!(&take::<8>(bytes, &mut cursor), b"FSAVSYN1");
    assert_eq!(u16::from_le_bytes(take(bytes, &mut cursor)), 1);
    let audio_frames_per_video_frame = u32::from_le_bytes(take(bytes, &mut cursor));
    let endpoint_drift_audio_frames = i64::from_le_bytes(take(bytes, &mut cursor));
    let marker_count = u32::from_le_bytes(take(bytes, &mut cursor));
    let markers = (0..marker_count)
        .map(|_| AudioVideoSyncMarker {
            video_tick: i64::from_le_bytes(take(bytes, &mut cursor)),
            audio_tick: i64::from_le_bytes(take(bytes, &mut cursor)),
            audio_frame_offset: u64::from_le_bytes(take(bytes, &mut cursor)),
        })
        .collect();
    assert_eq!(
        cursor,
        bytes.len(),
        "canonical receipt has no trailing bytes"
    );
    AudioVideoAlignment {
        audio_frames_per_video_frame,
        markers,
        endpoint_drift_audio_frames,
    }
}

fn alignment_fixture() -> AudioVideoAlignment {
    AudioVideoAlignment {
        audio_frames_per_video_frame: 2_000,
        markers: vec![
            AudioVideoSyncMarker {
                video_tick: -2,
                audio_tick: -4_000,
                audio_frame_offset: 0,
            },
            AudioVideoSyncMarker {
                video_tick: -1,
                audio_tick: -2_000,
                audio_frame_offset: 2_000,
            },
            AudioVideoSyncMarker {
                video_tick: 0,
                audio_tick: 0,
                audio_frame_offset: 4_000,
            },
        ],
        endpoint_drift_audio_frames: -1,
    }
}

#[test]
fn g0_public_alignment_receipt_is_canonical_and_independently_decodable() {
    let alignment = alignment_fixture();
    let first = with_cx(false, |cx| {
        encode_audio_video_alignment_receipt(&alignment, 3, 98, cx)
            .expect("fixture fits exact receipt ceilings")
    });
    let second = with_cx(false, |cx| {
        encode_audio_video_alignment_receipt(&alignment, 3, 98, cx)
            .expect("same input remains admissible")
    });

    assert_eq!(
        first, second,
        "canonical bytes and identity are deterministic"
    );
    assert_eq!(first.bytes().len(), 98);
    assert_eq!(decode_alignment_wire(first.bytes()), alignment);
}

#[test]
fn g4_public_alignment_receipt_refuses_budget_and_cancellation() {
    let alignment = alignment_fixture();
    assert_eq!(
        with_cx(false, |cx| {
            encode_audio_video_alignment_receipt(&alignment, 2, 98, cx)
        }),
        Err(CinematicReceiptError::BudgetExceeded),
    );
    assert_eq!(
        with_cx(false, |cx| {
            encode_audio_video_alignment_receipt(&alignment, 3, 97, cx)
        }),
        Err(CinematicReceiptError::BudgetExceeded),
    );
    assert_eq!(
        with_cx(true, |cx| {
            encode_audio_video_alignment_receipt(&alignment, 3, 98, cx)
        }),
        Err(CinematicReceiptError::Cancelled),
    );
}
