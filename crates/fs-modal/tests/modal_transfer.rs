//! G1/G3 physical modal transfer: independent field, momentum and energy checks.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::ModePair;
use fs_modal::time::ModalAcousticState;
use fs_modal::transfer::{
    ModalMassTransferControls, ModalTransferControls, ModalTransferError,
    transfer_modal_state, transfer_modal_state_with_mass_change,
};
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

fn mass_controls() -> ModalMassTransferControls {
    ModalMassTransferControls {
        transfer: controls(),
        maximum_relative_momentum_residual: 1e-12,
    }
}

#[test]
fn changed_scalar_inertia_preserves_displacement_and_momentum_with_signed_work() {
    let previous_mass = diagonal(&[4.0]);
    let stiffness = diagonal(&[16.0]);
    let old = [mode(4.0, &[0.5])];
    let states = [state(2.0, 6.0)];
    // Independently: u=1 m, p=12 kg m/s, E_old=26 J. The same
    // stiffness under a new inertia has E_new=p²/(2m)+ku²/2.
    with_cx(&CancelGate::new(), |cx| {
        for target_m in [1.0_f64, 4.0, 9.0] {
            let target_mass = diagonal(&[target_m]);
            let new = [mode(16.0 / target_m, &[1.0 / target_m.sqrt()])];
            let result = transfer_modal_state_with_mass_change(
                &previous_mass,
                &target_mass,
                &stiffness,
                &old,
                &new,
                &states,
                &[0.0],
                mass_controls(),
                cx,
            )
            .unwrap();
            let transfer = result.transfer;
            close(
                transfer.states[0].displacement_m_sqrt_kg / target_m.sqrt(),
                1.0,
            );
            close(
                transfer.states[0].velocity_m_sqrt_kg_per_s * target_m.sqrt(),
                12.0,
            );
            close(transfer.previous_modal_energy_j, 26.0);
            close(transfer.target_full_energy_j, 72.0 / target_m + 8.0);
            close(transfer.retained_modal_energy_j, 72.0 / target_m + 8.0);
            close(transfer.parameter_work_j, 72.0 / target_m - 18.0);
            assert!(result.relative_momentum_residual <= 1e-12);
            assert!(transfer.kinetic_residual_j < 1e-20);
        }
        let mut limited = mass_controls();
        limited.transfer.maximum_abs_parameter_work_j = 9.0;
        assert!(matches!(
            transfer_modal_state_with_mass_change(
                &previous_mass,
                &diagonal(&[9.0]),
                &stiffness,
                &old,
                &[mode(16.0 / 9.0, &[1.0 / 3.0])],
                &states,
                &[0.0],
                limited,
                cx
            ),
            Err(ModalTransferError::Limit {
                quantity: "absolute parameter work",
                ..
            })
        ));
    });
    assert_eq!(states, [state(2.0, 6.0)]);
}

#[test]
fn dense_changed_mass_and_rotated_degenerate_modes_preserve_physical_state() {
    let previous_mass = diagonal(&[1.0, 1.0]);
    let old = [mode(9.0, &[1.0, 0.0]), mode(16.0, &[0.0, 1.0])];
    let states = [state(0.4, 2.0), state(-0.3, -1.0)];
    let mut mass = Coo::new(2, 2);
    let mut stiffness = Coo::new(2, 2);
    for (i, j, value) in [(0, 0, 2.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 2.0)] {
        mass.push(i, j, value);
        stiffness.push(i, j, 25.0 * value);
    }
    let a = 1.0 / 6.0_f64.sqrt();
    let b = 1.0 / 2.0_f64.sqrt();
    // Rotate the entire mass-orthonormal repeated-eigenvalue subspace;
    // neither mode indices nor coordinate-wise mass scaling can do this.
    let new = [
        mode(25.0, &[0.6 * a + 0.8 * b, 0.6 * a - 0.8 * b]),
        mode(25.0, &[-0.8 * a + 0.6 * b, -0.8 * a - 0.6 * b]),
    ];
    with_cx(&CancelGate::new(), |cx| {
        let mass = mass.assemble();
        let stiffness = stiffness.assemble();
        let result = transfer_modal_state_with_mass_change(
            &previous_mass,
            &mass,
            &stiffness,
            &old,
            &new,
            &states,
            &[0.0; 2],
            mass_controls(),
            cx,
        )
        .unwrap();
        let transfer = &result.transfer;
        let mut u = [0.0; 2];
        let mut v = [0.0; 2];
        for (mode, state) in new.iter().zip(&transfer.states) {
            for i in 0..2 {
                u[i] += mode.phi[i] * state.displacement_m_sqrt_kg;
                v[i] += mode.phi[i] * state.velocity_m_sqrt_kg_per_s;
            }
        }
        close(u[0], 0.4);
        close(u[1], -0.3);
        close(v[0], 5.0 / 3.0);
        close(v[1], -4.0 / 3.0);
        close(2.0 * v[0] + v[1], 2.0);
        close(v[0] + 2.0 * v[1], -1.0);
        close(transfer.previous_modal_energy_j, 3.94);
        close(transfer.target_full_energy_j, 67.0 / 12.0);
        close(transfer.retained_modal_energy_j, 67.0 / 12.0);
        close(transfer.parameter_work_j, 67.0 / 12.0 - 3.94);
        assert!(result.relative_momentum_residual < 1e-12);
        let replay = transfer_modal_state_with_mass_change(
            &previous_mass,
            &mass,
            &stiffness,
            &old,
            &new,
            &states,
            &[0.0; 2],
            mass_controls(),
            cx,
        )
        .unwrap();
        assert_eq!(transfer.states, replay.transfer.states);
        assert_eq!(
            result.relative_momentum_residual.to_bits(),
            replay.relative_momentum_residual.to_bits()
        );
    });
}

#[test]
fn changed_mass_truncation_accounts_for_target_velocity_energy() {
    let previous_mass = diagonal(&[4.0, 9.0]);
    let mass = diagonal(&[1.0, 4.0]);
    let stiffness = diagonal(&[4.0, 36.0]);
    let old = [mode(4.0, &[0.5, 0.0]), mode(9.0, &[0.0, 1.0 / 3.0])];
    let new = [mode(4.0, &[1.0, 0.0])];
    let states = [state(2.0, 4.0), state(3.0, 6.0)];
    with_cx(&CancelGate::new(), |cx| {
        assert!(matches!(
            transfer_modal_state_with_mass_change(
                &previous_mass,
                &mass,
                &stiffness,
                &old,
                &new,
                &states,
                &[0.0; 2],
                mass_controls(),
                cx
            ),
            Err(ModalTransferError::Limit { .. })
        ));
        let mut permitted = mass_controls();
        permitted.transfer.maximum_displacement_residual_m_sqrt_kg = 2.0;
        permitted.transfer.maximum_kinetic_residual_j = 40.5;
        permitted.transfer.maximum_elastic_residual_j = 18.0;
        let transfer = transfer_modal_state_with_mass_change(
            &previous_mass,
            &mass,
            &stiffness,
            &old,
            &new,
            &states,
            &[0.0; 2],
            permitted,
            cx,
        )
        .unwrap()
        .transfer;
        assert_eq!(transfer.states, vec![state(1.0, 8.0)]);
        close(transfer.previous_modal_energy_j, 74.5);
        close(transfer.target_full_energy_j, 92.5);
        close(transfer.retained_modal_energy_j, 34.0);
        close(transfer.parameter_work_j, 18.0);
        close(transfer.displacement_residual_m_sqrt_kg, 2.0);
        close(transfer.kinetic_residual_j, 40.5);
        close(transfer.elastic_residual_j, 18.0);
        close(transfer.projection_energy_residual_j, 0.0);
    });
}

#[test]
fn changed_mass_refuses_bad_mass_space_source_normalization_and_solve_budget() {
    let mass = diagonal(&[4.0]);
    let stiffness = diagonal(&[16.0]);
    let old = [mode(4.0, &[0.5])];
    let new = [mode(16.0 / 9.0, &[1.0 / 3.0])];
    let states = [state(2.0, 6.0)];
    with_cx(&CancelGate::new(), |cx| {
        for target in [diagonal(&[-9.0]), diagonal(&[0.0]), diagonal(&[9.0, 9.0])] {
            assert!(
                transfer_modal_state_with_mass_change(
                    &mass,
                    &target,
                    &stiffness,
                    &old,
                    &new,
                    &states,
                    &[0.0],
                    mass_controls(),
                    cx
                )
                .is_err()
            );
        }
        // Source normalization belongs to the OLD mass, never the target mass.
        assert!(matches!(
            transfer_modal_state_with_mass_change(
                &mass,
                &diagonal(&[9.0]),
                &stiffness,
                &[mode(4.0, &[1.0 / 3.0])],
                &new,
                &states,
                &[0.0],
                mass_controls(),
                cx
            ),
            Err(ModalTransferError::Limit {
                quantity: "mass orthogonality",
                ..
            })
        ));
        for tolerance in [-1.0, 1.0, f64::NAN, f64::INFINITY] {
            let mut invalid = mass_controls();
            invalid.maximum_relative_momentum_residual = tolerance;
            assert!(matches!(
                transfer_modal_state_with_mass_change(
                    &mass,
                    &diagonal(&[9.0]),
                    &stiffness,
                    &old,
                    &new,
                    &states,
                    &[0.0],
                    invalid,
                    cx
                ),
                Err(ModalTransferError::Invalid(_))
            ));
        }
    });
    let cancelled = CancelGate::new();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        assert!(matches!(
            transfer_modal_state_with_mass_change(
                &mass,
                &diagonal(&[9.0]),
                &stiffness,
                &old,
                &new,
                &states,
                &[0.0],
                mass_controls(),
                cx
            ),
            Err(ModalTransferError::Cancelled)
        ));
    });
    assert_eq!(states, [state(2.0, 6.0)]);
}
