use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchAuthority, NormalPatchView, PartialSlipCheckpoint,
    PartialSlipError, PartialSlipInterface, PartialSlipKinematics, PartialSlipLaw,
    PartialSlipParameters, PartialSlipState, PartialSlipStateKind, TangentFrame,
};

fn parameters() -> PartialSlipParameters {
    PartialSlipParameters {
        static_mu: 0.8,
        kinetic_mu: 0.4,
        tangential_stiffness_n_per_m: 10_000.0,
        torsional_stiffness_nm_per_rad: 100.0,
        torsional_capacity_factor: 0.5,
        partial_slip_onset_fraction: 0.5,
        partial_slip_hardening_fraction: 0.4,
    }
}
fn law(p: PartialSlipParameters) -> PartialSlipLaw {
    PartialSlipLaw::new(
        "two-region/source-response",
        "synthetic/varying-coefficients",
        p,
    )
    .unwrap()
}
fn patch() -> NormalPatchView {
    NormalPatchView::new(
        "patch",
        "normal",
        "fixture",
        NormalPatchAuthority::SyntheticFixture,
        100.0,
        0.02,
        0.01,
        1.0e-4,
    )
    .unwrap()
}
fn interface() -> PartialSlipInterface {
    PartialSlipInterface::new(
        "A->B",
        "history",
        "fixture",
        NormalPatchAuthority::SyntheticFixture,
    )
    .unwrap()
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn key() -> GeneralizedWorkOwnership {
    GeneralizedWorkOwnership::new("patch", "interval", "x", "y", "spin").unwrap()
}
fn rest() -> PartialSlipKinematics {
    PartialSlipKinematics {
        creepage: [0.0; 2],
        rolling_speed_mps: 1.0,
        torsional_spin_rad_per_s: 0.0,
        contact_spin_rad_per_s: 0.0,
        dt_s: 0.1,
    }
}
fn checkpoint(state: PartialSlipState) -> PartialSlipCheckpoint {
    PartialSlipCheckpoint::new(patch(), interface(), law(parameters()), state).unwrap()
}
fn preload() -> PartialSlipState {
    PartialSlipState::from_region_checkpoint(
        [[1.0e-4, -2.0e-4]; 2],
        [3.0e-4; 2],
        [false; 2],
        Some(frame()),
        7,
    )
    .unwrap()
}
// Independent energy oracle from the public regional strains and supplied law.
fn energy(p: PartialSlipParameters, state: &PartialSlipState) -> f64 {
    let e = state.region_elastic_displacement_m();
    let twist = state.region_elastic_twist_rad();
    let h = p.partial_slip_hardening_fraction;
    [1.0 - h, h]
        .into_iter()
        .enumerate()
        .map(|(i, w)| {
            0.5 * w
                * (p.tangential_stiffness_n_per_m * (e[i][0] * e[i][0] + e[i][1] * e[i][1])
                    + p.torsional_stiffness_nm_per_rad * twist[i] * twist[i])
        })
        .sum()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3.0e-12 * a.abs().max(b.abs()).max(1.0e-5),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn stiffness_cycles_exchange_signed_work_without_inventing_heat_or_erasing_preload() {
    let old = checkpoint(preload());
    let snapshot = old.clone();
    let mut p = parameters();
    p.tangential_stiffness_n_per_m *= 2.0;
    p.torsional_stiffness_nm_per_rad *= 3.0;
    let changed = law(p)
        .advance_with_coefficients(&patch(), &interface(), frame(), rest(), &key(), &old)
        .unwrap();
    assert_eq!(old, snapshot);
    assert_eq!(
        changed.next_state.region_elastic_displacement_m(),
        old.state().region_elastic_displacement_m()
    );
    assert_eq!(
        changed.next_state.region_elastic_twist_rad(),
        old.state().region_elastic_twist_rad()
    );
    assert_eq!(
        changed.next_state.sliding_regions(),
        old.state().sliding_regions()
    );
    assert_eq!(changed.next_state.accepted_steps(), 8);
    close(changed.tangent_force_n[0], -2.0);
    close(changed.tangent_force_n[1], 4.0);
    close(changed.torsional_moment_nm, -0.09);
    let expected = 0.5 * 10_000.0 * (1.0e-8 + 4.0e-8) + 0.5 * 200.0 * 9.0e-8;
    close(changed.generalized_work.coefficient_work_j, expected);
    close(changed.storage_change_j, expected);
    assert_eq!(changed.generalized_work.work_into_interface_j, 0.0);
    assert_eq!(changed.dissipation.heat_j, 0.0);
    assert_eq!(changed.generalized_work.numerical_dissipation_j, 0.0);
    close(changed.generalized_work.closure_residual_j, 0.0);
    let restored = law(parameters())
        .advance_with_coefficients(
            &patch(),
            &interface(),
            frame(),
            rest(),
            &key(),
            &changed.checkpoint,
        )
        .unwrap();
    close(restored.generalized_work.coefficient_work_j, -expected);
    close(restored.stored_energy_j, energy(parameters(), old.state()));
    assert_eq!(restored.dissipation.heat_j, 0.0);
}

#[test]
fn mixed_motion_uses_original_return_map_and_closes_against_the_original_stored_energy() {
    let old = checkpoint(
        PartialSlipState::from_region_checkpoint(
            [[2.0e-4, -1.0e-4], [-3.0e-4, 4.0e-4]],
            [1.0e-3, -2.0e-3],
            [true, false],
            Some(frame()),
            9,
        )
        .unwrap(),
    );
    let mut p = parameters();
    p.tangential_stiffness_n_per_m *= 1.5;
    p.torsional_stiffness_nm_per_rad *= 0.7;
    p.static_mu = 0.55;
    p.kinetic_mu = 0.3;
    let law = law(p);
    let motion = PartialSlipKinematics {
        creepage: [0.04, -0.025],
        rolling_speed_mps: 0.5,
        torsional_spin_rad_per_s: 0.8,
        contact_spin_rad_per_s: 2.0,
        dt_s: 0.3,
    };
    let transported_frame = TangentFrame::new([0.3, 0.4, 0.5], [0.6, -0.2, 0.1]).unwrap();
    let direct = law
        .advance(
            &patch(),
            &interface(),
            transported_frame,
            motion,
            &key(),
            old.state(),
        )
        .unwrap();
    let changed = law
        .advance_with_coefficients(
            &patch(),
            &interface(),
            transported_frame,
            motion,
            &key(),
            &old,
        )
        .unwrap();
    assert_eq!(changed.wrench, direct.wrench);
    assert_eq!(changed.next_state, direct.next_state);
    assert_eq!(changed.dissipation, direct.dissipation);
    assert_eq!(
        changed.generalized_work.numerical_dissipation_j,
        direct.generalized_work.numerical_dissipation_j
    );
    assert!(changed.dissipation.heat_j > 0.0);
    assert!(changed.generalized_work.numerical_dissipation_j > 0.0);
    close(changed.stored_energy_j, energy(p, &changed.next_state));
    close(
        changed.storage_change_j,
        energy(p, &changed.next_state) - energy(parameters(), old.state()),
    );
    close(
        changed.generalized_work.coefficient_work_j,
        energy(p, old.state()) - energy(parameters(), old.state()),
    );
    let d = changed.stored_energy_j - energy(parameters(), old.state());
    close(
        d + changed.dissipation.heat_j + changed.generalized_work.numerical_dissipation_j,
        changed.generalized_work.work_into_interface_j
            + changed.generalized_work.coefficient_work_j,
    );
    close(changed.generalized_work.closure_residual_j, 0.0);
}

#[test]
fn decreased_friction_relaxes_retained_strain_with_no_fictitious_parameter_heat() {
    let old = checkpoint(
        PartialSlipState::from_region_checkpoint(
            [[0.001, 0.0]; 2],
            [0.0; 2],
            [false; 2],
            Some(frame()),
            3,
        )
        .unwrap(),
    );
    let mut p = parameters();
    p.static_mu = 0.02;
    p.kinetic_mu = 0.01;
    let changed = law(p)
        .advance_with_coefficients(&patch(), &interface(), frame(), rest(), &key(), &old)
        .unwrap();
    assert_eq!(changed.state, PartialSlipStateKind::GrossSlide);
    assert_eq!(changed.generalized_work.coefficient_work_j, 0.0);
    assert_eq!(changed.generalized_work.work_into_interface_j, 0.0);
    assert!(changed.dissipation.heat_j > 0.0);
    assert!(changed.generalized_work.numerical_dissipation_j > 0.0);
    close(
        changed.stored_energy_j
            + changed.dissipation.heat_j
            + changed.generalized_work.numerical_dissipation_j,
        energy(parameters(), old.state()),
    );
    assert!(changed.next_state.region_elastic_displacement_m()[0][0] < 0.001);
}

#[test]
fn unchanged_coefficients_are_bit_identical_to_the_original_update() {
    let old = checkpoint(preload());
    let law = law(parameters());
    let motion = PartialSlipKinematics {
        creepage: [0.13, -0.04],
        torsional_spin_rad_per_s: 1.5,
        ..rest()
    };
    assert_eq!(
        law.advance_with_coefficients(&patch(), &interface(), frame(), motion, &key(), &old)
            .unwrap(),
        law.advance(&patch(), &interface(), frame(), motion, &key(), old.state())
            .unwrap()
    );
}

#[test]
fn ordinary_restore_and_source_or_region_changes_remain_strict_and_atomic() {
    let old = checkpoint(preload());
    let snapshot = old.clone();
    let mut p = parameters();
    p.tangential_stiffness_n_per_m *= 2.0;
    assert!(matches!(
        law(p).restore_checkpoint(&patch(), &interface(), &old),
        Err(PartialSlipError::CheckpointIdentityMismatch {
            field: "law_parameters"
        })
    ));
    let foreign = PartialSlipLaw::new("two-region/source-response", "another-source", p).unwrap();
    assert!(matches!(
        foreign.advance_with_coefficients(&patch(), &interface(), frame(), rest(), &key(), &old),
        Err(PartialSlipError::CheckpointIdentityMismatch {
            field: "law_source_id"
        })
    ));
    p.partial_slip_hardening_fraction = 0.7;
    assert!(matches!(
        law(p).advance_with_coefficients(&patch(), &interface(), frame(), rest(), &key(), &old),
        Err(PartialSlipError::CheckpointIdentityMismatch {
            field: "region_stiffness_fraction"
        })
    ));
    assert_eq!(old, snapshot);
}

#[test]
fn changing_normal_load_reuses_the_original_history_and_work_balance() {
    let old = checkpoint(preload());
    let evolving_patch = NormalPatchView::new(
        "patch",
        "normal",
        "fixture",
        NormalPatchAuthority::SyntheticFixture,
        0.1,
        0.015,
        0.008,
        5.0e-5,
    )
    .unwrap();
    let mut p = parameters();
    p.tangential_stiffness_n_per_m *= 1.5;
    let step = law(p)
        .advance_with_coefficients(&evolving_patch, &interface(), frame(), rest(), &key(), &old)
        .unwrap();
    assert_eq!(step.checkpoint.normal_patch(), &evolving_patch);
    assert_eq!(
        step.next_state.accepted_steps(),
        old.state().accepted_steps() + 1
    );
    assert!(step.dissipation.heat_j > 0.0);
    close(step.generalized_work.work_into_interface_j, 0.0);
    close(
        step.stored_energy_j - energy(parameters(), old.state())
            + step.dissipation.heat_j
            + step.generalized_work.numerical_dissipation_j,
        step.generalized_work.coefficient_work_j,
    );
    assert_eq!(
        law(p)
            .restore_checkpoint(&evolving_patch, &interface(), &step.checkpoint)
            .unwrap(),
        step.next_state
    );
    let changed_identity = NormalPatchView::new(
        "another-patch",
        "normal",
        "fixture",
        NormalPatchAuthority::SyntheticFixture,
        0.1,
        0.015,
        0.008,
        5.0e-5,
    )
    .unwrap();
    assert!(
        law(p)
            .advance_with_coefficients(
                &changed_identity,
                &interface(),
                frame(),
                rest(),
                &key(),
                &old
            )
            .is_err()
    );
    for mutating_onset in [false, true] {
        let mut p = parameters();
        if mutating_onset {
            p.partial_slip_onset_fraction = 0.6;
        } else {
            p.torsional_capacity_factor = 0.6;
        }
        assert!(
            law(p)
                .advance_with_coefficients(&patch(), &interface(), frame(), rest(), &key(), &old)
                .is_err()
        );
    }
}
