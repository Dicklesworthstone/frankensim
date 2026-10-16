//! The fs-exec → fs-rand stream-key bridge contract (bead wf9.7.1):
//! checked field widths, refusal-not-truncation at every boundary, and
//! the iteration no-slot rule — plus the collision demonstration that
//! motivates refusing (truncation WOULD collide distinct identities).

use fs_rand::{EXEC_KEY_BRIDGE_VERSION, ExecKeyBridgeError, StreamKey};

#[test]
fn bridge_accepts_in_range_fields_losslessly() {
    let k = StreamKey::from_exec_parts(u64::MAX, u64::from(u32::MAX), u64::from(u32::MAX), 0)
        .expect("boundary values fit");
    assert_eq!(k.seed, u64::MAX);
    assert_eq!(k.kernel, u32::MAX);
    assert_eq!(k.tile, u32::MAX);
    assert_eq!(EXEC_KEY_BRIDGE_VERSION, 1, "bump only with justification");
}

#[test]
fn bridge_refuses_at_every_truncation_boundary() {
    // One past the u32 boundary refuses — and names the field.
    assert_eq!(
        StreamKey::from_exec_parts(7, 1u64 << 32, 0, 0),
        Err(ExecKeyBridgeError::KernelOverflow {
            kernel_id: 1u64 << 32
        })
    );
    assert_eq!(
        StreamKey::from_exec_parts(7, 0, 1u64 << 32, 0),
        Err(ExecKeyBridgeError::TileOverflow { tile: 1u64 << 32 })
    );
    assert_eq!(
        StreamKey::from_exec_parts(7, 0, 0, 3),
        Err(ExecKeyBridgeError::IterationUnrepresentable { iteration: 3 })
    );
    // The COLLISION the refusal prevents: truncating kernel_id would
    // alias (1 << 32) | 5 onto 5 — two distinct logical streams, one
    // physical stream. The bridge refuses instead.
    let truncated_twin = StreamKey::from_exec_parts(7, 5, 9, 0).expect("in range");
    let would_collide = (1u64 << 32) | 5;
    assert!(
        StreamKey::from_exec_parts(7, would_collide, 9, 0).is_err(),
        "the identity that would alias {truncated_twin:?} is refused"
    );
    // Errors are teaching, not bare bools.
    let msg = ExecKeyBridgeError::IterationUnrepresentable { iteration: 3 }.to_string();
    assert!(msg.contains("draw index") && msg.contains("v1"), "{msg}");
}
