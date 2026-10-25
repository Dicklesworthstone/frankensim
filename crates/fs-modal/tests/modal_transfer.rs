//! G1/G3 physical modal transfer: independent field, momentum and energy checks.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::ModePair;
use fs_modal::time::ModalAcousticState;
use fs_modal::transfer::{ModalTransferControls, ModalTransferError, transfer_modal_state};
use fs_sparse::{Coo, Csr};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(gate, arena,
            StreamKey { seed: 34, kernel_id: 34, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic))
    })
}
fn diagonal(values: &[f64]) -> Csr {
    let mut matrix = Coo::new(values.len(), values.len());
    for (i, &value) in values.iter().enumerate() { matrix.push(i, i, value); }
    matrix.assemble()
}
fn mode(lambda: f64, phi: &[f64]) -> ModePair {
    ModePair { lambda, phi: phi.to_vec(), residual: 0.0, interval: (lambda, lambda) }
}
fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState { displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v }
}
fn controls() -> ModalTransferControls {
    ModalTransferControls {
        maximum_dofs: 16,
        maximum_modes: 16,
        maximum_factor_bytes: 1 << 20,
        maximum_mass_orthogonality_error: 1e-11,
        maximum_relative_eigen_residual: 1e-11,
        maximum_displacement_residual_m_sqrt_kg: 1e-12,
        maximum_kinetic_residual_j: 1e-20,
        maximum_elastic_residual_j: 1e-20,
        maximum_abs_parameter_work_j: 1000.0,
        energy_tolerance_j: 1e-10,
    }
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 1e-11 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn sign_and_order_changes_preserve_physical_fields_and_canonical_momentum() {
    let mass = diagonal(&[4.0, 9.0]);
    let stiffness = diagonal(&[16.0, 81.0]);
    let old = [mode(4.0, &[0.5, 0.0]), mode(9.0, &[0.0, 1.0/3.0])];
    let new = [mode(9.0, &[0.0, -1.0/3.0]), mode(4.0, &[0.5, 0.0])];
    let states = [state(2.0, 4.0), state(-3.0, 5.0)];
    with_cx(&CancelGate::new(), |cx| {
        let transfer = transfer_modal_state(&mass, &stiffness, &old, &new,
            &states, &[0.0, 0.0], controls(), cx).unwrap();
        // The physical state is u=[1,-1], v=[2,5/3], p=Mv=[8,15].
        // Reordering the amplitudes by sorted frequency alone is incorrect.
        close(transfer.states[0].displacement_m_sqrt_kg, 3.0);
        close(transfer.states[0].velocity_m_sqrt_kg_per_s, -5.0);
        close(transfer.states[1].displacement_m_sqrt_kg, 2.0);
        close(transfer.states[1].velocity_m_sqrt_kg_per_s, 4.0);
        close(0.5 * transfer.states[1].displacement_m_sqrt_kg, 1.0);
        close(-transfer.states[0].displacement_m_sqrt_kg / 3.0, -1.0);
        close(4.0 * 0.5 * transfer.states[1].velocity_m_sqrt_kg_per_s, 8.0);
        close(-3.0 * transfer.states[0].velocity_m_sqrt_kg_per_s, 15.0);
        close(transfer.previous_modal_energy_j, 69.0);
        close(transfer.target_full_energy_j, 69.0);
        close(transfer.retained_modal_energy_j, 69.0);
        close(transfer.parameter_work_j, 0.0);
        assert!(transfer.displacement_residual_m_sqrt_kg < 1e-12);
        assert!(transfer.kinetic_residual_j < 1e-20);
        assert!(transfer.elastic_residual_j < 1e-20);
        let replay = transfer_modal_state(&mass, &stiffness, &old, &new,
            &states, &[0.0, 0.0], controls(), cx).unwrap();
        assert_eq!(transfer.states, replay.states);
        assert_eq!(transfer.target_full_energy_j.to_bits(), replay.target_full_energy_j.to_bits());
        assert_eq!(states, [state(2.0, 4.0), state(-3.0, 5.0)]);
    });
}

#[test]
fn degenerate_subspace_rotation_transfers_motion_without_mode_matching() {
    let mass = diagonal(&[1.0, 1.0, 1.0]);
    let stiffness = diagonal(&[25.0, 25.0, 81.0]);
    let old = [mode(25.0, &[1.0, 0.0, 0.0]), mode(25.0, &[0.0, 1.0, 0.0]),
        mode(81.0, &[0.0, 0.0, 1.0])];
    let new = [mode(81.0, &[0.0, 0.0, 1.0]), mode(25.0, &[0.6, 0.8, 0.0]),
        mode(25.0, &[-0.8, 0.6, 0.0])];
    let states = [state(2.0, 0.3), state(-1.0, -0.4), state(0.5, 0.8)];
    with_cx(&CancelGate::new(), |cx| {
        let transfer = transfer_modal_state(&mass, &stiffness, &old, &new,
            &states, &[0.0; 3], controls(), cx).unwrap();
        for (actual, expected) in transfer.states.iter()
            .zip([state(0.5, 0.8), state(0.4, -0.14), state(-2.2, -0.48)]) {
            close(actual.displacement_m_sqrt_kg, expected.displacement_m_sqrt_kg);
            close(actual.velocity_m_sqrt_kg_per_s, expected.velocity_m_sqrt_kg_per_s);
        }
        // Independent physical-coordinate storage:
        // (0.3²+0.4²+0.8² + 25*(2²+1²) + 81*0.5²)/2.
        close(transfer.previous_modal_energy_j, 73.07);
        close(transfer.retained_modal_energy_j, 73.07);
        close(transfer.parameter_work_j, 0.0);
        assert!(transfer.displacement_residual_m_sqrt_kg < 1e-12);
        assert!(transfer.projection_energy_residual_j.abs() < 1e-10);
    });
}

#[test]
fn shifted_equilibrium_and_changed_stiffness_report_independent_parameter_work() {
    let mass = diagonal(&[4.0, 9.0]);
    let stiffness = diagonal(&[100.0, 144.0]);
    let old = [mode(4.0, &[0.5, 0.0]), mode(9.0, &[0.0, 1.0/3.0])];
    let new = [mode(25.0, &[0.5, 0.0]), mode(16.0, &[0.0, 1.0/3.0])];
    let states = [state(2.0, 4.0), state(-3.0, 5.0)];
    with_cx(&CancelGate::new(), |cx| {
        let transfer = transfer_modal_state(&mass, &stiffness, &old, &new,
            &states, &[0.25, 0.5], controls(), cx).unwrap();
        // u_new-relative=[1.25,-0.5], while physical velocity is unchanged.
        close(transfer.states[0].displacement_m_sqrt_kg, 2.5);
        close(transfer.states[1].displacement_m_sqrt_kg, -1.5);
        close(transfer.states[0].velocity_m_sqrt_kg_per_s, 4.0);
        close(transfer.states[1].velocity_m_sqrt_kg_per_s, 5.0);
        // Kinetic=20.5; target elastic=(100*1.25²+144*0.5²)/2=96.125.
        close(transfer.previous_modal_energy_j, 69.0);
        close(transfer.target_full_energy_j, 116.625);
        close(transfer.retained_modal_energy_j, 116.625);
        close(transfer.parameter_work_j, 47.625);
        let mut tight = controls();
        tight.maximum_abs_parameter_work_j = 47.0;
        assert!(matches!(transfer_modal_state(&mass, &stiffness, &old, &new,
            &states, &[0.25, 0.5], tight, cx), Err(ModalTransferError::Limit {
                quantity: "absolute parameter work", ..
            })));
        assert_eq!(states, [state(2.0, 4.0), state(-3.0, 5.0)]);
    });
}

#[test]
fn truncation_requires_explicit_displacement_and_energy_allowances() {
    let mass = diagonal(&[1.0, 1.0]);
    let stiffness = diagonal(&[4.0, 9.0]);
    let old = [mode(4.0, &[1.0, 0.0]), mode(9.0, &[0.0, 1.0])];
    let states = [state(0.5, 3.0), state(2.0, -4.0)];
    with_cx(&CancelGate::new(), |cx| {
        assert!(matches!(transfer_modal_state(&mass, &stiffness, &old, &old[..1],
            &states, &[0.0; 2], controls(), cx), Err(ModalTransferError::Limit { .. })));
        let mut allowed = controls();
        allowed.maximum_displacement_residual_m_sqrt_kg = 2.0;
        allowed.maximum_kinetic_residual_j = 8.0;
        allowed.maximum_elastic_residual_j = 18.0;
        let transfer = transfer_modal_state(&mass, &stiffness, &old, &old[..1],
            &states, &[0.0; 2], allowed, cx).unwrap();
        assert_eq!(transfer.states, vec![state(0.5, 3.0)]);
        close(transfer.previous_modal_energy_j, 31.0);
        close(transfer.target_full_energy_j, 31.0);
        close(transfer.retained_modal_energy_j, 5.0);
        close(transfer.parameter_work_j, 0.0);
        close(transfer.displacement_residual_m_sqrt_kg, 2.0);
        close(transfer.kinetic_residual_j, 8.0);
        close(transfer.elastic_residual_j, 18.0);
        close(transfer.projection_energy_residual_j, 0.0);
        assert_eq!(states, [state(0.5, 3.0), state(2.0, -4.0)]);
    });
}

#[test]
fn malformed_basis_state_operator_budgets_and_cancellation_refuse() {
    let mass = diagonal(&[1.0, 1.0]);
    let stiffness = diagonal(&[4.0, 9.0]);
    let modes = [mode(4.0, &[1.0, 0.0]), mode(9.0, &[0.0, 1.0])];
    let states = [state(0.5, 3.0), state(2.0, -4.0)];
    with_cx(&CancelGate::new(), |cx| {
        let wrong = [mode(5.0, &[1.0, 0.0]), mode(9.0, &[0.0, 1.0])];
        assert!(matches!(transfer_modal_state(&mass, &stiffness, &modes, &wrong,
            &states, &[0.0; 2], controls(), cx), Err(ModalTransferError::Limit {
                quantity: "target eigen residual", ..
            })));
        let bad_basis = [mode(4.0, &[2.0, 0.0]), mode(9.0, &[0.0, 1.0])];
        assert!(transfer_modal_state(&mass, &stiffness, &bad_basis, &modes,
            &states, &[0.0; 2], controls(), cx).is_err());
        assert!(transfer_modal_state(&mass, &stiffness, &modes, &modes,
            &[state(f64::NAN, 0.0), states[1]], &[0.0; 2], controls(), cx).is_err());
        assert!(transfer_modal_state(&mass, &stiffness, &modes, &modes,
            &states, &[f64::INFINITY, 0.0], controls(), cx).is_err());
        assert!(transfer_modal_state(&diagonal(&[-1.0, 1.0]), &stiffness, &modes, &modes,
            &states, &[0.0; 2], controls(), cx).is_err());
        let mut asymmetric = Coo::new(2, 2);
        asymmetric.push(0, 0, 4.0);
        asymmetric.push(1, 1, 9.0);
        asymmetric.push(0, 1, 0.1);
        assert!(transfer_modal_state(&mass, &asymmetric.assemble(), &modes, &modes,
            &states, &[0.0; 2], controls(), cx).is_err());
        let mut limited = controls();
        limited.maximum_dofs = 1;
        assert!(matches!(transfer_modal_state(&mass, &stiffness, &modes, &modes,
            &states, &[0.0; 2], limited, cx), Err(ModalTransferError::Budget { .. })));
        assert!(transfer_modal_state(&mass, &stiffness, &modes, &modes,
            &states[..1], &[0.0; 2], controls(), cx).is_err());
    });
    let cancelled = CancelGate::new();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        assert!(matches!(transfer_modal_state(&mass, &stiffness, &modes, &modes,
            &states, &[0.0; 2], controls(), cx), Err(ModalTransferError::Cancelled)));
    });
    assert_eq!(states, [state(0.5, 3.0), state(2.0, -4.0)]);
}
