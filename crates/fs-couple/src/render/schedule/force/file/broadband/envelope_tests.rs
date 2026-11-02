use super::*;

const EXAMPLE: &[u8] =
    include_bytes!("../../../../../../../../examples/modal_enveloped_tones.performance");

#[test]
fn enveloped_v1_file_keeps_controls_pressure_work_and_heat_across_callbacks() {
    let mut whole = BroadbandModalPerformance::from_bytes(EXAMPLE, 4800).unwrap();
    let mut split = BroadbandModalPerformance::from_bytes(EXAMPLE, 37).unwrap();
    assert_eq!(whole.envelope_count(), 2);
    assert_eq!(whole.continuous_force_tones(), 3);
    assert_eq!(whole.compiled_controls(), 2);
    let mut a = vec![0.0; 4800];
    let mut b = vec![0.0; 4800];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(37) {
        split.block(chunk).unwrap();
    }
    assert_eq!(a, b);
    assert!(a.iter().any(|value| value.abs() > 1e-7));
    assert_eq!(whole.forces, vec![vec![0.0]]);
    assert_eq!(whole.next_event, whole.events.len());
    for (a, b) in whole.observers.iter().zip(&split.observers) {
        assert_eq!(a.model().states(), b.model().states());
        assert_eq!(a.last_modal_frame(), b.last_modal_frame());
        assert!(a.last_modal_frame().viscous_dissipation_j > 0.0);
    }
    assert_eq!(whole.samples_rendered(), 4800);
    assert_eq!(
        whole.info().input_hash,
        fs_blake3::hash_domain("org.frankensim.modal-performance.broadband.v1", EXAMPLE)
    );
    assert!(whole.block(&mut [0.0]).is_err());
}
