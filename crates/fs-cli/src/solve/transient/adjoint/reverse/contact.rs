//! Compose the original live-contact tangent with native capacity controls.

use super::{SolveRefusal, bad};

/// Additional sparse contact operators, complete capacity tangent and bounded
/// FGMRES workspace at one accepted endpoint. The native contact owner admits
/// up to 432 staged entries of 24 bytes per face pair/overlap triangle and eight
/// nodal work vectors. The retained receipt allowance is separate from source
/// packs already admitted by the forward path. This is not a peak-RSS claim.
pub(in crate::solve::transient) fn storage_bytes(
    nodes: usize,
    elements: usize,
    max_pairs: usize,
    restart: usize,
) -> Result<usize, SolveRefusal> {
    if nodes == 0 || elements == 0 || max_pairs == 0 || restart == 0 {
        return Err(bad(
            "capacity contact adjoints require a nonempty mesh, paired-face allowance and positive FGMRES restart",
        ));
    }
    // Match the original Robin tangent owner's bounded Arnoldi dimension.
    let restart = restart.min(nodes).min(64);
    max_pairs
        .checked_mul(432 * 24 + 2048)
        .and_then(|bytes| {
            elements
                .checked_mul(4096)
                .and_then(|v| bytes.checked_add(v))
        })
        .and_then(|bytes| {
            nodes
                .checked_mul(8 * (2 * restart + 24))
                .and_then(|v| bytes.checked_add(v))
        })
        .and_then(|bytes| {
            restart
                .checked_mul(restart + 8)
                .and_then(|v| v.checked_mul(8))
                .and_then(|v| bytes.checked_add(v))
        })
        .and_then(|bytes| bytes.checked_add(2048))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| bad("capacity contact adjoint workspace overflows"))
}

pub(super) fn fragment(steps: usize, max_pairs: usize, workspace: usize) -> String {
    format!(
        concat!(
            ",\"temperature_contact\":{{\"method\":\"complete-capacity-contact-tangent\",",
            "\"control_steps\":{steps},\"max_face_pairs\":{max_pairs},",
            "\"workspace_estimate_bytes\":{workspace},",
            "\"residual_target\":\"original-accepted-forward-physical-tolerance\",",
            "\"source_context\":\"original-pinned-law-at-independent-endpoint-traces\"}}"
        ),
        steps = steps,
        max_pairs = max_pairs,
        workspace = workspace
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contact_workspace_checks_native_face_and_krylov_allowances() {
        for invalid in [(0, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0)] {
            assert!(storage_bytes(invalid.0, invalid.1, invalid.2, invalid.3).is_err());
        }
        assert!(storage_bytes(8, 2, usize::MAX, 8).is_err());
        assert!(storage_bytes(usize::MAX, 2, 1, 8).is_err());
        assert!(storage_bytes(8, usize::MAX, 1, 8).is_err());
        assert!(storage_bytes(80, 2, 2, 32).unwrap() > storage_bytes(80, 2, 1, 16).unwrap());
        assert_eq!(
            storage_bytes(80, 2, 1, 64).unwrap(),
            storage_bytes(80, 2, 1, usize::MAX).unwrap()
        );
    }
}
