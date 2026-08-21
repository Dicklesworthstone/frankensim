//! E6.2 lane probe: the canonical digest through THIS crate's native
//! build (vs the wasm lane and the main-workspace lane).
fn main() {
    let mut slot = fs_flyer_wasm::engine::EngineSlot::default();
    let init = slot.init(1903, 1.294, 11.0, 0, 0, 18.3, 120, false, false);
    assert!(init.starts_with("{\"ok\""), "{init}");
    for _ in 0..120 {
        slot.step(false, 0.0, 0.0);
    }
    println!("wasmcrate-native {}", slot.digest());
}
