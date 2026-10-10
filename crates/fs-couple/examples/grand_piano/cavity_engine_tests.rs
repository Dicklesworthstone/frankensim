//! Reciprocal air storage through the existing struck piano and its clock.
//! These authored overlap matrices are mechanical integration fixtures, not
//! a geometry calibration or a measured piano cavity.
use super::*;
use fs_couple::vibroacoustic::CavityModes;

fn piano() -> Instrument {
    let course = super::super::geometry::demonstration_scale().unwrap()[48];
    Instrument::new(vec![course], &super::super::board::demonstration(),
        48_000, 4, 12, true).unwrap()
}

fn cavity(ports: usize, standing: bool, drag: f64, coupled: bool) -> CavityCoupling {
    let modes = if standing { 2 } else { 1 };
    let air = CavityModes {
        omegas: if standing { vec![0.0, std::f64::consts::TAU * 300.0] } else { vec![0.0] },
        lambdas: if standing { vec![0.3, 0.15] } else { vec![0.3] },
        interface: vec![vec![1.0]; modes], loss_factor: 0.0, rho0: 1.2, c0: 343.0,
    };
    let mut overlap = vec![0.0; ports * modes];
    if coupled {
        for port in 0..ports {
            overlap[port * modes] = 0.15 / (port + 1) as f64;
            if standing {
                let sign = if port % 2 == 0 { 1.0 } else { -1.0 };
                overlap[port * modes + 1] = sign * 0.25 / (port + 1) as f64;
            }
        }
    }
    let damping = if standing { vec![0.0, drag] } else { vec![0.0] };
    CavityCoupling::new(&air, ports, &overlap, &damping).unwrap()
}

fn exterior(ports: usize) -> radiation::Model {
    radiation::Model { ports, poles: vec![radiation::Pole {
        omega: std::f64::consts::TAU * 420.0, zeta: 0.25,
        coupling: (0..ports).map(|j| 250.0 / (j + 1) as f64).collect(),
    }] }
}

fn balanced(piano: &Instrument) {
    let defect = piano.accounting.input_work_j - piano.energy_j()
        - piano.accounting.dissipated_j();
    assert!(defect.abs() < 1e-7, "combined piano/air energy defect {defect:e} J");
}

#[test]
fn struck_piano_drives_cavity_storage_and_changes_the_actual_board_motion() {
    let mut loaded = piano();
    let mut bare = piano();
    loaded.configure_cavity(&cavity(loaded.bank.board_count, true, 30.0, true)).unwrap();
    assert!(loaded.has_cavity());
    assert_eq!(loaded.cavity_energy_j(), 0.0);
    assert_eq!(loaded.accounting.input_work_j, 0.0);
    loaded.note_on(69, 1.0).unwrap();
    bare.note_on(69, 1.0).unwrap();
    let mut peak_storage = 0.0_f64;
    let mut peak_output = 0.0_f64;
    let mut peak_difference = 0.0_f64;
    for _ in 0..2400 {
        let actual = loaded.step().unwrap();
        let original = bare.step().unwrap();
        peak_storage = peak_storage.max(loaded.cavity_energy_j());
        peak_output = peak_output.max(actual.abs());
        peak_difference = peak_difference.max((actual - original).abs());
    }
    assert!(peak_storage > 1e-15 && loaded.accounting.cavity_loss_j > 1e-15);
    assert!(peak_output > 1e-10 && peak_difference > 1e-12);
    assert_ne!(loaded.bank.q, bare.bank.q);
    assert_eq!(bare.cavity_energy_j(), 0.0);
    assert_eq!(bare.accounting.cavity_loss_j, 0.0);
    balanced(&loaded);
    balanced(&bare);
}

#[test]
fn uniform_compression_stores_energy_without_inventing_acoustic_damping() {
    let mut loaded = piano();
    let mut bare = piano();
    let model = cavity(loaded.bank.board_count, false, 0.0, true);
    assert_eq!(model.total_modes(), model.structural_modes());
    loaded.configure_cavity(&model).unwrap();
    loaded.note_on(69, 1.0).unwrap();
    bare.note_on(69, 1.0).unwrap();
    let mut peak_storage = 0.0_f64;
    for _ in 0..1200 {
        loaded.step().unwrap();
        bare.step().unwrap();
        peak_storage = peak_storage.max(loaded.cavity_energy_j());
    }
    assert!(peak_storage > 1e-15);
    assert_ne!(loaded.bank.q, bare.bank.q);
    assert_eq!(loaded.accounting.cavity_loss_j, 0.0);
    balanced(&loaded);
}

#[test]
fn cavity_and_exterior_radiation_share_contact_and_keep_separate_real_losses() {
    let mut loaded = piano();
    let mut exterior_only = piano();
    let ports = loaded.bank.board_count;
    let outside = exterior(ports);
    loaded.configure_radiation(&outside).unwrap();
    loaded.configure_cavity(&cavity(ports, true, 30.0, true)).unwrap();
    exterior_only.configure_radiation(&outside).unwrap();
    assert!(loaded.has_cavity() && loaded.has_radiation());
    assert_eq!(loaded.hammer_contact_count(), exterior_only.hammer_contact_count());
    loaded.note_on(69, 1.0).unwrap();
    exterior_only.note_on(69, 1.0).unwrap();
    let mut cavity_peak = 0.0_f64;
    let mut exterior_peak = 0.0_f64;
    for _ in 0..1800 {
        loaded.step().unwrap();
        exterior_only.step().unwrap();
        cavity_peak = cavity_peak.max(loaded.cavity_energy_j());
        exterior_peak = exterior_peak.max(loaded.radiation_energy_j());
    }
    assert!(cavity_peak > 1e-15 && exterior_peak > 1e-15);
    assert!(loaded.accounting.cavity_loss_j > 1e-15);
    assert!(loaded.accounting.radiation_loss_j > 1e-15);
    assert!(loaded.accounting.felt_loss_j > 0.0);
    assert_ne!(loaded.bank.q, exterior_only.bank.q);
    balanced(&loaded);
    balanced(&exterior_only);
}

#[test]
fn failed_output_frame_restores_both_air_states_and_future_board_traces() {
    let mut a = piano();
    let mut b = piano();
    let ports = a.bank.board_count;
    for instrument in [&mut a, &mut b] {
        instrument.configure_cavity(&cavity(ports, true, 30.0, true)).unwrap();
        instrument.configure_radiation(&exterior(ports)).unwrap();
        instrument.note_on(69, 1.0).unwrap();
    }
    let mut trace = vec![0.0; a.board_trace_len()];
    for _ in 0..480 {
        a.step_with_board_trace(&mut trace).unwrap();
        for substep in 0..b.substeps {
            b.mechanics_step().unwrap();
            assert_eq!(&trace[substep * ports..(substep + 1) * ports],
                &b.bank.v[b.bank.modes.len()..]);
        }
    }
    let energy = a.energy_j();
    let cavity_loss = a.accounting.cavity_loss_j;
    let exterior_loss = a.accounting.radiation_loss_j;
    let total_loss = a.accounting.dissipated_j();
    // Both air half-flows run before the invalid damper refuses the frame.
    // Replay below also detects internal air state that energy alone could miss.
    a.damper_drag_ns_m = f64::NAN;
    assert!(a.step_with_board_trace(&mut trace).is_err());
    assert_eq!(a.bank.q, b.bank.q);
    assert_eq!(a.bank.v, b.bank.v);
    assert_eq!(a.energy_j(), energy);
    assert_eq!(a.cavity_energy_j(), b.cavity_energy_j());
    assert_eq!(a.radiation_energy_j(), b.radiation_energy_j());
    assert_eq!(a.accounting.cavity_loss_j, cavity_loss);
    assert_eq!(a.accounting.radiation_loss_j, exterior_loss);
    assert_eq!(a.accounting.dissipated_j(), total_loss);
    a.damper_drag_ns_m = 0.4;
    let mut reference = trace.clone();
    for _ in 0..96 {
        assert_eq!(a.step_with_board_trace(&mut trace).unwrap(),
            b.step_with_board_trace(&mut reference).unwrap());
        assert_eq!(trace, reference);
        assert_eq!(a.bank.q, b.bank.q);
        assert_eq!(a.bank.v, b.bank.v);
        assert_eq!(a.cavity_energy_j(), b.cavity_energy_j());
        assert_eq!(a.radiation_energy_j(), b.radiation_energy_j());
    }
    balanced(&a);
}

#[test]
fn cavity_admission_is_cold_atomic_and_zero_overlap_keeps_the_original_note() {
    let mut a = piano();
    let mut b = piano();
    let ports = a.bank.board_count;
    assert!(a.configure_cavity(&cavity(ports - 1, true, 30.0, true)).is_err());
    assert!(!a.has_cavity());
    assert_eq!(a.energy_j(), 0.0);
    let zero = cavity(ports, true, 30.0, false);
    a.configure_cavity(&zero).unwrap();
    assert!(a.configure_cavity(&zero).is_err());
    a.note_on(69, 0.5).unwrap();
    b.note_on(69, 0.5).unwrap();
    assert!(b.configure_cavity(&zero).is_err());
    assert!(!b.has_cavity());
    for _ in 0..480 {
        assert_eq!(a.step().unwrap(), b.step().unwrap());
        assert_eq!(a.bank.q, b.bank.q);
        assert_eq!(a.bank.v, b.bank.v);
    }
    assert_eq!(a.cavity_energy_j(), 0.0);
    assert_eq!(a.accounting.cavity_loss_j, 0.0);
    assert_eq!(a.accounting.dissipated_j(), b.accounting.dissipated_j());
    let mut held = piano();
    held.silent_key_down(69).unwrap();
    assert!(held.configure_cavity(&zero).is_err());
    held.note_off(69).unwrap();
    held.configure_cavity(&zero).unwrap();
}
