use super::*;

const EXAMPLE: &[u8] =
    include_bytes!("../../../../../../../../examples/modal_swept_tones.performance");

#[test]
fn swept_v1_file_retains_callback_bits_and_ringdown_after_source_release() {
    let mut whole = BroadbandModalPerformance::from_bytes(EXAMPLE, 4800).unwrap();
    let mut split = BroadbandModalPerformance::from_bytes(EXAMPLE, 37).unwrap();
    assert_eq!(whole.frequency_sweep_count(), 3);
    assert_eq!(whole.envelope_count(), 3);
    assert_eq!(whole.continuous_force_tones(), 3);
    assert_eq!(whole.compiled_controls(), 2);
    assert!(whole.maximum_sweep_force_error_bound_n() > 0.0);
    assert!(whole.maximum_sweep_force_error_bound_n() <= 0.00001);
    let mut a = vec![0.0; 4800];
    let mut b = vec![0.0; 4800];
    whole.block(&mut a).unwrap();
    for chunk in b[..4000].chunks_mut(37) {
        split.block(chunk).unwrap();
    }
    let release_energy = split.observers[0].last_modal_frame().total_modal_energy_j;
    for chunk in b[4000..].chunks_mut(37) {
        split.block(chunk).unwrap();
        for observer in &split.observers {
            assert_eq!(observer.last_modal_frame().input_work_j, 0.0);
        }
    }
    assert_eq!(a, b);
    assert!(a.iter().any(|value| value.abs() > 1e-7));
    assert!(a[4000..].iter().any(|value| value.abs() > 1e-8));
    assert_eq!(whole.forces, vec![vec![0.0]]);
    assert_eq!(whole.next_event, whole.events.len());
    assert!(split.observers[0].last_modal_frame().total_modal_energy_j < release_energy);
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

#[test]
fn empty_sweep_section_retains_existing_envelope_performance_bits() {
    let original =
        include_str!("../../../../../../../../examples/modal_enveloped_tones.performance");
    let empty = original.replace("broadband 1", "frequency_sweeps 0\nbroadband 1");
    let mut a = BroadbandModalPerformance::from_bytes(original.as_bytes(), 4800).unwrap();
    let mut b = BroadbandModalPerformance::from_bytes(empty.as_bytes(), 4800).unwrap();
    assert_eq!(b.frequency_sweep_count(), 0);
    assert_eq!(b.maximum_sweep_force_error_bound_n(), 0.0);
    let mut pa = vec![0.0; 4800];
    let mut pb = vec![0.0; 4800];
    a.block(&mut pa).unwrap();
    b.block(&mut pb).unwrap();
    assert_eq!(pa, pb);
    assert_eq!(
        a.observers[0].last_modal_frame(),
        b.observers[0].last_modal_frame()
    );
}
