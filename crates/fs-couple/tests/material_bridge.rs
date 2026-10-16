//! Material strings on a fixed-nut / moving-bridge support.
//!
//! The integration regressions announced by fbebd8dd4 ("bind material
//! strings to fixed-nut moving-bridge dynamics") never landed. That
//! commit wrote this file as 66 bytes of tool error text, not Rust, so
//! every `cargo test -p fs-couple --test '*'` failed to compile. The
//! announced cases are listed in the ignored test below. They need
//! restoring from the remote lane's prepared validation bundle, and
//! until then they are NOT claimed as tested.

/// Placeholder for the lost regressions. Run it explicitly to see what is missing.
#[test]
#[ignore = "fbebd8dd4's material_bridge regressions were never committed; restore them from the validation bundle"]
fn material_bridge_regressions_are_missing() {
    panic!(
        "missing: direct shared-bridge dynamics; material/torsion retention; \
         stationary-gas substitutions; invalid frame/clock inputs; failed-step rollback"
    );
}
