//! Integration test asserting the Euler disc requirements traceability matrix
//! (bead `frankensim-euler-disc-emergent-flagship-t6314.8.13`).

use std::fs;
use std::path::Path;

#[test]
fn test_euler_disc_verification_matrix_is_valid_and_complete() {
    let manifest_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/e2e/campaigns/euler_disc/verification_matrix.toml");
    
    assert!(
        manifest_path.exists(),
        "Traceability matrix must exist at {}",
        manifest_path.display()
    );

    let content = fs::read_to_string(&manifest_path)
        .expect("Failed to read verification_matrix.toml");

    assert!(
        content.contains("org.frankensim.euler-disc.verification-matrix.v1"),
        "Matrix must declare exact schema version"
    );

    // Verify all 12 mandatory requirements are registered
    for i in 1..=12 {
        let req_id = format!("REQ-ED-{:03}", i);
        assert!(
            content.contains(&req_id),
            "Traceability matrix missing required requirement {}",
            req_id
        );
    }
}
