use fs_contact::tangential::smooth::{
    SmoothAuthorityPolicy, SmoothRegularization, SmoothTangentialAdapter, SmoothTangentialError,
    SmoothTangentialQuery,
};
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchAuthority, NormalPatchView, PARTIAL_SLIP_MODEL_ID,
    PartialSlipInterface, PartialSlipKinematics, PartialSlipLaw, PartialSlipParameters,
    PartialSlipStateKind, TangentFrame,
};

const TOLERANCE: f64 = 2.0e-10;

fn close(left: f64, right: f64) {
    assert!(
        (left - right).abs() <= TOLERANCE * left.abs().max(right.abs()).max(1.0),
        "left={left:.17e}, right={right:.17e}"
    );
}

fn close_vec(left: [f64; 3], right: [f64; 3]) {
    for index in 0..3 {
        close(left[index], right[index]);
    }
}

fn patch(authority: NormalPatchAuthority) -> NormalPatchView {
    NormalPatchView::new(
        "patch-a",
        "normal-card-v1",
        "fixture/normal-patch",
        authority,
        100.0,
        0.02,
        0.01,
        1.0e-4,
    )
    .expect("patch")
}

fn interface(authority: NormalPatchAuthority) -> PartialSlipInterface {
    PartialSlipInterface::new(
        "body-a->support-b",
        "history-a",
        "fixture/interface",
        authority,
    )
    .expect("interface")
}

fn law() -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        "fixture/partial-slip-law",
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: 0.4,
            tangential_stiffness_n_per_m: 10_000.0,
            torsional_stiffness_nm_per_rad: 100.0,
            torsional_capacity_factor: 0.5,
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .expect("law")
}

fn adapter(scale: f64) -> SmoothTangentialAdapter {
    SmoothTangentialAdapter::new(
        "smooth-adapter-v1",
        "fixture/smooth-config",
        SmoothRegularization {
            creepage_scale: scale,
            torsional_spin_scale_rad_per_s: scale,
            tangent_probe_creepage: 1.0e-5,
            tangent_probe_spin_rad_per_s: 1.0e-5,
        },
        SmoothAuthorityPolicy::test_only(),
    )
    .expect("adapter")
}

fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).expect("frame")
}

fn query(version: u64, interval: &str, creepage: [f64; 2], spin: f64) -> SmoothTangentialQuery {
    SmoothTangentialQuery {
        query_id: format!("query-{interval}"),
        expected_state_version: version,
        frame: frame(),
        kinematics: PartialSlipKinematics {
            creepage,
            rolling_speed_mps: 10.0,
            torsional_spin_rad_per_s: spin,
            dt_s: 0.001,
        },
        work_ownership: GeneralizedWorkOwnership::new("patch-a", interval, "qx", "qy", "qspin")
            .expect("work key"),
    }
}

fn state(
    adapter: &SmoothTangentialAdapter,
) -> fs_contact::tangential::smooth::SmoothTangentialState {
    adapter
        .initial_state(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            8,
        )
        .expect("state")
}

#[test]
fn fixed_branch_tangent_agrees_with_independent_centered_difference() {
    let adapter = adapter(1.0e-4);
    let state = state(&adapter);
    let base = query(0, "fixed", [0.1, 0.02], 0.01);
    let tangent = adapter
        .fixed_branch_tangent(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state,
            &base,
        )
        .expect("stick tangent");
    assert_eq!(tangent.branch, PartialSlipStateKind::Sticking);

    let h = 1.0e-5;
    let mut plus = base.clone();
    plus.kinematics.creepage[0] += h;
    let mut minus = base;
    minus.kinematics.creepage[0] -= h;
    let plus = adapter
        .prepare(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state,
            &plus,
        )
        .expect("plus");
    let minus = adapter
        .prepare(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state,
            &minus,
        )
        .expect("minus");
    close(
        tangent.derivative[0][0],
        (plus.residual.tangent_force_n[0] - minus.residual.tangent_force_n[0]) / (2.0 * h),
    );
}

#[test]
fn regularization_refines_toward_the_delegated_law() {
    let coarse = adapter(0.1);
    let fine = adapter(1.0e-4);
    let raw = query(0, "regularization", [0.7, 0.0], 0.0);
    let direct = law()
        .advance(
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            frame(),
            raw.kinematics,
            &raw.work_ownership,
            &fs_tribo::partial_slip::PartialSlipState::zero(),
        )
        .expect("delegated reference");
    let coarse_force = coarse
        .prepare(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state(&coarse),
            &raw,
        )
        .expect("coarse")
        .residual
        .tangent_force_n[0];
    let fine_force = fine
        .prepare(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state(&fine),
            &raw,
        )
        .expect("fine")
        .residual
        .tangent_force_n[0];
    assert!(
        (fine_force - direct.tangent_force_n[0]).abs()
            < (coarse_force - direct.tangent_force_n[0]).abs()
    );
}

#[test]
fn transition_neighborhoods_are_explicit_and_branch_derivatives_refuse() {
    let adapter = adapter(1.0e-4);
    let state = state(&adapter);
    let inputs = [
        ([0.2, 0.0], PartialSlipStateKind::Sticking),
        ([0.6, 0.0], PartialSlipStateKind::PartialSlip),
        ([1.2, 0.0], PartialSlipStateKind::GrossSlide),
    ];
    for (creepage, branch) in inputs {
        assert_eq!(
            adapter
                .prepare(
                    &law(),
                    &patch(NormalPatchAuthority::SyntheticFixture),
                    &interface(NormalPatchAuthority::SyntheticFixture),
                    &state,
                    &query(0, "transition", creepage, 0.0),
                )
                .expect("transition response")
                .branch,
            branch
        );
    }
    // The declared map subtracts approximately one scale at this magnitude,
    // so this raw value regularizes exactly onto the 0.5-capacity onset.
    let at_onset = query(0, "onset", [0.4001, 0.0], 0.0);
    assert!(matches!(
        adapter.fixed_branch_tangent(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state,
            &at_onset,
        ),
        Err(SmoothTangentialError::NoDerivativeOnBranchChange { .. })
    ));
}

#[test]
fn rollback_retry_reversal_and_recontact_are_transactional() {
    let adapter = adapter(1.0e-4);
    let patch = patch(NormalPatchAuthority::SyntheticFixture);
    let interface = interface(NormalPatchAuthority::SyntheticFixture);
    let first_state = state(&adapter);
    let first = adapter
        .prepare(
            &law(),
            &patch,
            &interface,
            &first_state,
            &query(0, "a", [1.2, 0.0], 0.0),
        )
        .expect("first candidate");
    assert_eq!(adapter.rollback(&first).expect("rollback"), first_state);
    let retried = adapter
        .prepare(
            &law(),
            &patch,
            &interface,
            &first_state,
            &query(0, "a", [1.2, 0.0], 0.0),
        )
        .expect("retry after rollback");
    let committed = adapter.commit(&first_state, &retried).expect("commit");
    let reversed = adapter
        .prepare(
            &law(),
            &patch,
            &interface,
            &committed,
            &query(1, "b", [-1.2, 0.0], 0.0),
        )
        .expect("reversal");
    assert!(retried.residual.tangent_force_n[0] < 0.0);
    assert!(reversed.residual.tangent_force_n[0] > 0.0);
    let recontact = adapter
        .initial_state(&law(), &patch, &interface, 8)
        .expect("explicit fresh recontact state");
    assert_eq!(recontact.committed_version(), 0);
}

#[test]
fn so2_covariance_and_power_invariance_hold_for_the_adapter_receipt() {
    let adapter = adapter(1.0e-4);
    let state = state(&adapter);
    let angle = 0.37;
    let c = angle.cos();
    let s = angle.sin();
    let raw = [0.6, 0.2];
    let rotated_raw = [c * raw[0] + s * raw[1], -s * raw[0] + c * raw[1]];
    let base = adapter
        .prepare(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state,
            &query(0, "base", raw, 0.3),
        )
        .expect("base");
    let mut rotated = query(0, "rotated", rotated_raw, 0.3);
    rotated.frame = frame().rotated(angle).expect("rotated frame");
    let rotated = adapter
        .prepare(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            &state,
            &rotated,
        )
        .expect("rotated");
    close_vec(
        base.action_reaction.action_on_declared_body.force_n,
        rotated.action_reaction.action_on_declared_body.force_n,
    );
    close_vec(
        base.action_reaction.action_on_declared_body.torque_nm,
        rotated.action_reaction.action_on_declared_body.torque_nm,
    );
    close(
        base.residual.endpoint_relative_power_w,
        rotated.residual.endpoint_relative_power_w,
    );
}

#[test]
fn action_reaction_duplicate_work_and_stale_future_refuse() {
    let adapter = adapter(1.0e-4);
    let patch = patch(NormalPatchAuthority::SyntheticFixture);
    let interface = interface(NormalPatchAuthority::SyntheticFixture);
    let initial = state(&adapter);
    let receipt = adapter
        .prepare(
            &law(),
            &patch,
            &interface,
            &initial,
            &query(0, "unique", [0.6, 0.1], 0.2),
        )
        .expect("receipt");
    close_vec(
        receipt.action_reaction.reaction_on_counterbody.force_n,
        receipt
            .action_reaction
            .action_on_declared_body
            .force_n
            .map(|value| -value),
    );
    close_vec(
        receipt.action_reaction.reaction_on_counterbody.torque_nm,
        receipt
            .action_reaction
            .action_on_declared_body
            .torque_nm
            .map(|value| -value),
    );
    let committed = adapter.commit(&initial, &receipt).expect("commit");
    assert!(matches!(
        adapter.prepare(
            &law(),
            &patch,
            &interface,
            &committed,
            &query(1, "unique", [0.6, 0.1], 0.2)
        ),
        Err(SmoothTangentialError::DuplicateWorkKey)
    ));
    assert!(matches!(
        adapter.prepare(
            &law(),
            &patch,
            &interface,
            &committed,
            &query(0, "stale", [0.0, 0.0], 0.0)
        ),
        Err(SmoothTangentialError::StaleState { .. })
    ));
    assert!(matches!(
        adapter.prepare(
            &law(),
            &patch,
            &interface,
            &committed,
            &query(2, "future", [0.0, 0.0], 0.0)
        ),
        Err(SmoothTangentialError::FutureState { .. })
    ));
}

#[test]
fn authority_policy_and_checkpoint_identity_refuse_without_promotion() {
    let adapter = SmoothTangentialAdapter::new(
        "smooth-adapter-v1",
        "fixture/smooth-config",
        adapter(1.0e-4).regularization(),
        SmoothAuthorityPolicy {
            allow_caller_declared: true,
            allow_synthetic_fixture: false,
            allow_estimated: false,
        },
    )
    .expect("strict adapter");
    assert!(matches!(
        adapter.initial_state(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            2,
        ),
        Err(SmoothTangentialError::AuthorityRefused { .. })
    ));

    let accepting = adapter(1.0e-4);
    let checkpoint = state(&accepting).checkpoint();
    let altered = SmoothTangentialAdapter::new(
        "other-adapter",
        "fixture/smooth-config",
        accepting.regularization(),
        SmoothAuthorityPolicy::test_only(),
    )
    .expect("altered identity");
    assert!(matches!(
        altered.restore_checkpoint(
            &law(),
            &patch(NormalPatchAuthority::SyntheticFixture),
            &interface(NormalPatchAuthority::SyntheticFixture),
            checkpoint,
        ),
        Err(SmoothTangentialError::CheckpointIdentityMismatch {
            field: "adapter_id"
        })
    ));
}
