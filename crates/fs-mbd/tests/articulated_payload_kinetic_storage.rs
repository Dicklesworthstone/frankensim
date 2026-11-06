//! Inventory work remains observable beside a much larger unchanged spin.

use fs_ga::{Mat3, Se3, Twist, Vec3};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::motion::inventory::release::{
    PayloadReleaseInput, SeparationImpulse, release_payload, release_payloads,
};
use fs_mbd::articulated::motion::inventory::supported::release::{
    SupportedPayloadReleaseControls, SupportedPayloadSeparationControls,
    release_payload_with_impulse,
};
use fs_mbd::articulated::motion::inventory::{
    PayloadCaptureInput, PayloadControls, PayloadError, PayloadState, capture_payload,
    capture_payloads, supported,
};
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

fn sphere(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [mass, 0., 0., 0., mass, 0., 0., 0., mass],
        },
    )
    .unwrap()
}

fn body() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(2.),
    )])
    .unwrap()
}

fn velocity(speed: f64) -> Twist {
    Twist::new(Vec3::new(0., 0., 1e12), Vec3::new(0., 0., speed))
}

fn free_state(dofs: usize) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(Se3::identity(), velocity(0.)),
        joint_position: vec![0.; dofs],
        joint_velocity: vec![0.; dofs],
    }
}

fn payload(mass: f64, speed: f64) -> PayloadState {
    PayloadState {
        inertia: sphere(mass),
        world_from_body: Se3::identity(),
        twist_body: velocity(speed),
    }
}

fn controls(tolerances: &[f64], work: f64) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 1e-12,
        angular_momentum_tolerance_nms: 1e-9,
        linear_momentum_tolerance_ns: 1e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 1e-9,
        maximum_energy_exchange_j: work,
        maximum_constraint_work_j: 1e-9,
        minimum_reciprocal_condition: 1e-10,
    }
}

fn kick(impulse: f64) -> SeparationImpulse {
    SeparationImpulse {
        impulse_on_payload_link_ns: Vec3::new(0., 0., impulse),
        ..SeparationImpulse::default()
    }
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-9, "{actual} != {expected}");
}

#[test]
fn free_and_locked_capture_release_retain_real_storage_and_work() {
    for locked in [false, true] {
        let model = if locked {
            ArticulatedModel::new(vec![
                Link::new("base", None, Se3::identity(), JointModel::FIXED, sphere(1.)),
                Link::new(
                    "locked rotor",
                    Some(0),
                    Se3::identity(),
                    JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
                    sphere(1.),
                ),
            ])
            .unwrap()
        } else {
            body()
        };
        let rows = if locked {
            vec![JointConstraint {
                coordinate: ConstraintCoordinate::Angle,
                coefficients: vec![1.],
                position: 0.,
                position_tolerance: 1e-9,
                velocity_tolerance: 1e-9,
            }]
        } else {
            vec![]
        };
        let tolerances = vec![1e-9; model.dof_count()];
        let state = free_state(model.dof_count());
        let target = model.link_count() - 1;
        let captured = capture_payload(
            &model,
            &state,
            target,
            payload(2., 2.),
            &rows,
            controls(&tolerances, 2.),
            || false,
        )
        .unwrap();
        // The complete absolute energies are both about 2e24 J.
        assert_eq!(
            captured.balance.kinetic_before_j,
            captured.balance.kinetic_after_j
        );
        close(captured.balance.kinetic_energy_change_j, -2.);
        close(captured.removed_relative_energy_j, 2.);
        close(captured.balance.inventory_decomposition_j, 0.);
        close(captured.balance.energy_residual_j, 0.);
        assert!(matches!(
            capture_payload(
                &model,
                &state,
                target,
                payload(2., 2.),
                &rows,
                controls(&tolerances, 1.99),
                || false,
            ),
            Err(PayloadError::EnergyExchange { .. })
        ));
        let released = release_payload(
            &captured.model_after,
            &captured.state_after,
            &captured.attachment,
            Some(kick(2.)),
            &rows,
            controls(&tolerances, 2.),
            || false,
        )
        .unwrap();
        assert_eq!(
            released.balance.kinetic_before_j,
            released.balance.kinetic_after_j
        );
        close(released.balance.kinetic_energy_change_j, 2.);
        close(released.supplied_work_j, 2.);
        close(released.constraint_work_j, 0.);
        close(released.balance.inventory_decomposition_j, 0.);
        close(released.balance.energy_residual_j, 0.);
    }
}

#[test]
fn simultaneous_inventory_keeps_loss_and_cancelling_source_recoil() {
    let model = body();
    let state = free_state(0);
    let incoming = [
        PayloadCaptureInput {
            target_link: 0,
            payload: payload(1., 2.),
        },
        PayloadCaptureInput {
            target_link: 0,
            payload: payload(1., -2.),
        },
    ];
    let captured =
        capture_payloads(&model, &state, &incoming, &[], controls(&[], 4.), || false).unwrap();
    assert_eq!(
        captured.balance.kinetic_before_j,
        captured.balance.kinetic_after_j
    );
    close(captured.balance.kinetic_energy_change_j, -4.);
    close(captured.removed_relative_energy_j, 4.);
    close(captured.balance.energy_residual_j, 0.);
    let outgoing = [
        PayloadReleaseInput {
            attachment: &captured.attachments[1],
            impulse: Some(kick(1.)),
        },
        PayloadReleaseInput {
            attachment: &captured.attachments[0],
            impulse: Some(kick(-1.)),
        },
    ];
    let released = release_payloads(
        &captured.model_after,
        &captured.state_after,
        &outgoing,
        &[],
        controls(&[], 1.),
        || false,
    )
    .unwrap();
    assert_eq!(released.state_after, captured.state_after);
    assert_eq!(
        released.balance.kinetic_before_j,
        released.balance.kinetic_after_j
    );
    close(released.balance.kinetic_energy_change_j, 1.);
    close(released.supplied_work_j, 1.);
    close(released.balance.inventory_decomposition_j, 0.);
    close(released.balance.energy_residual_j, 0.);
}

#[test]
fn finite_inertia_assembly_defect_is_neither_capture_loss_nor_release_work() {
    let model = body();
    let state = free_state(0);
    let tiny_moment = 2.0_f64.powi(-54);
    let incoming = PayloadState {
        inertia: SpatialInertia::new(
            2.,
            Vec3::default(),
            Mat3 {
                m: [
                    tiny_moment,
                    0.,
                    0.,
                    0.,
                    tiny_moment,
                    0.,
                    0.,
                    0.,
                    tiny_moment,
                ],
            },
        )
        .unwrap(),
        ..payload(2., 0.)
    };
    // The mass increment is representable; its small spin inertia is lost when
    // added to the carrier's I=2. At this spin that storage is still material.
    let defect = 0.5 * tiny_moment * 1e12 * 1e12;
    assert!(matches!(
        capture_payload(&model, &state, 0, incoming, &[], controls(&[], 0.), || {
            false
        }),
        Err(PayloadError::EnergyResidual { .. })
    ));
    let mut loose = controls(&[], 0.);
    loose.energy_tolerance_j = 1e8;
    let captured = capture_payload(&model, &state, 0, incoming, &[], loose, || false).unwrap();
    assert_eq!(
        captured.balance.kinetic_before_j,
        captured.balance.kinetic_after_j
    );
    close(captured.removed_relative_energy_j, 0.);
    close(captured.balance.kinetic_energy_change_j, -defect);
    close(captured.balance.inventory_decomposition_j, -defect);
    close(captured.balance.energy_residual_j, -defect);
    assert!(matches!(
        release_payload(
            &captured.model_after,
            &captured.state_after,
            &captured.attachment,
            None,
            &[],
            controls(&[], 0.),
            || false,
        ),
        Err(PayloadError::EnergyResidual { .. })
    ));
    let released = release_payload(
        &captured.model_after,
        &captured.state_after,
        &captured.attachment,
        None,
        &[],
        loose,
        || false,
    )
    .unwrap();
    close(released.supplied_work_j, 0.);
    close(released.balance.kinetic_energy_change_j, defect);
    close(released.balance.inventory_decomposition_j, defect);
    close(released.balance.energy_residual_j, defect);
}

#[test]
fn fixed_and_driven_supports_keep_separation_work_beside_common_spin() {
    for speed in [0., 1e9] {
        let model = body();
        let state = PrescribedBaseMotionState {
            time_s: 0.,
            base: BaseState::prescribed(Se3::identity(), velocity(speed), Twist::zero()),
            joint_position: vec![],
            joint_velocity: vec![],
        };
        let captured = supported::capture_payload(
            &model,
            &state,
            0,
            payload(2., speed + 2.),
            supported::SupportedPayloadControls {
                angular_momentum_tolerance_nms: 1e-9,
                linear_momentum_tolerance_ns: 1e-9,
                joint_momentum_tolerances: &[],
                energy_tolerance_j: 1e-9,
                maximum_removed_relative_energy_j: 4.,
            },
            || false,
        )
        .unwrap();
        close(captured.removed_relative_energy_j, 4.);
        close(captured.balance.accounting_kinetic_energy_change_j, -4.);
        close(captured.balance.kinetic_energy_change_j, -4. * speed - 4.);
        close(captured.balance.support_work_j, -4. * speed);
        close(captured.balance.energy_residual_j, 0.);
        let released = release_payload_with_impulse(
            &captured.model_after,
            &captured.state_after,
            &captured.attachment,
            kick(2.),
            SupportedPayloadSeparationControls {
                balance: SupportedPayloadReleaseControls {
                    angular_momentum_tolerance_nms: 1e-9,
                    linear_momentum_tolerance_ns: 1e-9,
                    joint_momentum_tolerances: &[],
                    energy_tolerance_j: 1e-9,
                },
                maximum_separation_work_j: 1.,
            },
            || false,
        )
        .unwrap();
        assert_eq!(
            released.balance.kinetic_before_j,
            released.balance.kinetic_after_j
        );
        close(released.balance.accounting_kinetic_energy_change_j, 1.);
        close(released.balance.kinetic_energy_change_j, 2. * speed + 1.);
        close(released.supplied_work_j, 1.);
        close(released.balance.support_work_j, 2. * speed);
        close(released.balance.accounting_inventory_decomposition_j, 0.);
        close(released.balance.energy_residual_j, 0.);
    }
}
