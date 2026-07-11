//! os_affinity capsule battery (moved from the capsule module to
//! keep it under the 300-line cap): structural refusal-or-proof on
//! every target, page audit, and hugepage advice.

use fs_substrate::os_affinity::{
    OsAffinityError, advise_hugepages, current_cpu, page_nodes, pin_current_thread,
};

#[test]
fn unsupported_targets_refuse_loudly_and_supported_ones_answer() {
    // On every target the API answers STRUCTURALLY: either a real
    // answer (Linux x86-64) or Unsupported — never a silent no-op.
    match current_cpu() {
        Ok(cpu) => {
            // If we can ask where we are, we can pin there and stay.
            pin_current_thread(&[cpu]).expect("pin to current cpu");
            assert_eq!(current_cpu().expect("still answers"), cpu);
        }
        Err(OsAffinityError::Unsupported(why)) => {
            assert!(!why.is_empty());
            assert!(matches!(
                pin_current_thread(&[0]),
                Err(OsAffinityError::Unsupported(_))
            ));
        }
        Err(e) => panic!("unexpected: {e}"),
    }
    // Argument validation is target-independent semantics on Linux;
    // elsewhere Unsupported wins (both are refusals, never no-ops).
    assert!(pin_current_thread(&[]).is_err());
}

#[test]
fn page_audit_reports_touched_pages_or_refuses() {
    let buf = vec![1u8; 1 << 20];
    match page_nodes(&buf, 4096) {
        Ok(nodes) => {
            // The buffer touches however many 4 KiB pages its byte range
            // spans, which depends on the allocation's alignment: an
            // unaligned 1 MiB buffer crosses one extra boundary (256 pages
            // when 4 KiB-aligned, 257 when not). aarch64 heap allocations
            // landed page-aligned here; x86-64's did not — compute the true
            // span from the actual start offset rather than assuming it.
            let start = buf.as_ptr() as usize;
            let expected = (start % 4096 + buf.len()).div_ceil(4096);
            assert_eq!(nodes.len(), expected);
            // Touched pages report a non-negative node id.
            assert!(nodes.iter().all(|&n| n >= 0), "touched pages have nodes");
        }
        Err(OsAffinityError::Unsupported(_)) => {}
        Err(e) => panic!("unexpected: {e}"),
    }
}

#[test]
fn hugepage_advice_reports_bytes_or_refuses() {
    let buf = vec![7u8; 8 << 20];
    match advise_hugepages(&buf, 4096) {
        Ok(advised) => {
            assert!(advised > 0 && advised <= buf.len());
            // Contents are untouched by advice.
            assert!(buf.iter().all(|&b| b == 7));
        }
        Err(OsAffinityError::Unsupported(_)) => {}
        Err(e) => panic!("unexpected: {e}"),
    }
    // Degenerate inputs are Ok(0) on supported targets, refusals elsewhere.
    match advise_hugepages(&[], 4096) {
        Ok(n) => assert_eq!(n, 0),
        Err(OsAffinityError::Unsupported(_)) => {}
        Err(e) => panic!("unexpected: {e}"),
    }
}
