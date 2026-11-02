//! G1/G3/G4: sourced adhesive islands, independent finite-inertia work and retry.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::spinning::SourceAdhesiveSlidingPatch;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::multiple::*;
use fs_mbd::tangential_friction::adhesive::{RigidAdhesiveContact, RigidAdhesiveControls};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::adhesive::spinning::contact::{AdhesiveWrenchConfig, AdhesiveWrenchMode};
use fs_tribo::partial_slip::TangentFrame;

const SPIN: f64 = 1000.0;
const STRENGTH: f64 = 2e6;
const IA: f64 = 1.2e-27;
const IB: f64 = 2.5e-27;
// Micron-body translation is much stiffer than spin at these nanometre pads.
// Resolve the physical interval without weakening the original support gates.
const RESOLVED_FRACTION: f64 = 1e-4;

fn close(actual: f64, expected: f64, absolute: f64) {
    common::close(actual, expected, absolute, 2e-8);
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn controls() -> RigidAdhesiveIslandControls {
    RigidAdhesiveIslandControls {
        mechanical: RigidAdhesiveControls {
            maximum_bodies: 2,
            geometry_tolerance_m: 1e-18,
            velocity_tolerance_m_s: 2e-12,
            angular_velocity_tolerance_rad_s: 2e-6,
            linear_momentum_tolerance_ns: 1e-30,
            angular_momentum_tolerance_nms: 1e-34,
            energy_tolerance_j: 5e-31,
            contact: AdhesiveWrenchConfig {
                maximum_static_iterations: 96,
                maximum_iterations: 256,
                maximum_evaluations: 1024,
                force_tolerance_n: 2e-18,
                moment_tolerance_nm: 2e-26,
                displacement_tolerance_m: 1e-25,
                energy_tolerance_j: 1e-31,
                relative_matrix_tolerance: 1e-13,
                maximum_force_n: 1e-6,
                maximum_spin_moment_nm: 1e-12,
            },
        },
        maximum_contacts: 4,
        maximum_matrix_terms: 144,
        maximum_sweeps: 512,
        maximum_evaluations: 65536,
    }
}

struct Fixture {
    bodies: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    radius: f64,
}
impl Fixture {
    fn new(radius: f64) -> Self {
        let properties = vec![
            MassProperties::new(2e-15, Vec3::ZERO, Vec3::new(0.8e-27, 1e-27, IA)).unwrap(),
            MassProperties::new(3e-15, Vec3::ZERO, Vec3::new(1.5e-27, 2e-27, IB)).unwrap(),
            MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.5, 0.7, 0.9)).unwrap(),
        ];
        let pose = Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap();
        let bodies = vec![
            RigidBodyState::new(pose, Vec3::ZERO, Vec3::new(0.0, 0.0, IA * SPIN)).unwrap(),
            RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO).unwrap(),
            // This energetic spectator must not dilute microscopic work gates
            // or count against the two genuinely participating bodies.
            RigidBodyState::new(
                Pose::new(Vec3::new(3.0, 2.0, 1.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(3e4, -2e4, 1e4),
                Vec3::new(4.0, 2.0, 1.0),
            )
            .unwrap(),
        ];
        Self {
            bodies,
            properties,
            radius,
        }
    }
    fn contacts<'a>(
        &self,
        source: SourceAdhesiveSlidingPatch<'a>,
    ) -> [RigidAdhesiveIslandContact<'a>; 2] {
        // Two distinct disks centred at +/-a on a common plane. Their edges
        // touch; they do not double-count the same solid-contact area.
        [1.0, -1.0].map(|sign| {
            let arm = Vec3::new(sign * self.radius, 0.0, 0.0);
            RigidAdhesiveIslandContact {
                patch_id: if sign > 0.0 { "right" } else { "left" },
                contact: RigidAdhesiveContact {
                    binding: RigidContactBinding {
                        body_a: 0,
                        body_b: 1,
                        arm_a_body: arm,
                        arm_b_body: arm,
                    },
                    frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
                    source,
                },
            }
        })
    }
}

fn reference(radius: f64, area: f64) -> (f64, f64, f64, f64) {
    // Independent uniform-disk integrals at |u|=a*omega, not pressure-weighted
    // coefficients or a reuse of the production junction evaluator.
    let force = STRENGTH * area * 8.0 / (3.0 * PI);
    let moment = STRENGTH * area * radius * 8.0 / (9.0 * PI);
    let torque = 2.0 * (radius * force + moment);
    let stop_time = SPIN / (torque * (1.0 / IA + 1.0 / IB));
    (force, moment, torque, stop_time)
}

#[test]
fn two_symmetric_junctions_transfer_exact_rim_torque_and_relative_energy_to_finite_bodies() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let pull_off =
        2.0 * PI * normal.effective_radius_m() * normal.normal_state().adhesion_energy_j_per_m2();
    for load in [0.0, -0.5 * pull_off] {
        let response = normal.evaluate_contact_load(load).unwrap();
        let source = shear
            .sliding_patch_dmt(normal.normal_state(), &response)
            .unwrap();
        let f = Fixture::new(response.contact_radius_m);
        let (force, moment, torque, stop_time) = reference(f.radius, response.contact_area_m2);
        let dt = RESOLVED_FRACTION * stop_time;
        let contacts = f.contacts(source);
        let prepared = PreparedRigidAdhesiveIsland::new(
            &f.bodies,
            &f.properties,
            &contacts,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(
            prepared.patch_ids().collect::<Vec<_>>(),
            vec!["left", "right"]
        );
        assert_eq!(prepared.body_indices(), &[0, 1]);
        assert!(prepared.detached_contacts().is_empty());
        assert!(core::ptr::eq(
            prepared.contacts()[0].contact.source.normal_state(),
            normal.normal_state()
        ));
        assert!(core::ptr::eq(
            prepared.contacts()[0].contact.source.shear_source(),
            &shear
        ));
        // Both patches measure exactly the same spin. The complete 6x6 Gram is
        // rank deficient and retains force-spin cross terms; no diagonal or
        // invented compliance can produce these actual shared-body rows.
        let mobility = prepared.scaled_mobility_per_kg();
        for column in 0..6 {
            close(mobility[2 * 6 + column], mobility[5 * 6 + column], 1.0);
        }
        assert!(mobility[1 * 6 + 2].abs() > 1e10);
        close(
            mobility[1 * 6 + 2],
            -f.radius * f.radius * (1.0 / IA + 1.0 / IB),
            1.0,
        );
        let candidate = prepared.advance(&CancelGate::new()).unwrap();
        let report = &candidate.report.mechanical;
        let omega_after = SPIN - torque * dt * (1.0 / IA + 1.0 / IB);
        let heat = torque * dt * 0.5 * (SPIN + omega_after);
        for (i, contact) in report.contacts.iter().enumerate() {
            assert_eq!(contact.contact.mode, AdhesiveWrenchMode::Sliding);
            close(contact.contact.force_n[0], 0.0, 2e-18);
            close(contact.contact.force_n[1], [force, -force][i], 3e-18);
            close(contact.contact.spin_moment_nm, -moment, 3e-26);
            close(contact.spin_after_rad_s, omega_after, 3e-6);
            close(
                contact.slip_velocity_after_m_s[1],
                [-1.0, 1.0][i] * f.radius * omega_after,
                3e-12,
            );
            close(contact.dissipated_energy_j, heat / 2.0, 1e-30);
            assert_eq!(contact.normal_velocity_before_m_s, 0.0);
            assert_eq!(contact.normal_velocity_after_m_s, 0.0);
        }
        close(
            candidate.bodies_after[0].angular_momentum_body().z,
            IA * SPIN - torque * dt,
            3e-33,
        );
        close(
            candidate.bodies_after[1].angular_momentum_body().z,
            torque * dt,
            3e-33,
        );
        assert!(candidate.bodies_after[1].angular_momentum_body().z > 0.0);
        close(report.kinetic_before_j, 0.5 * IA * SPIN * SPIN, 1e-33);
        close(
            report.kinetic_after_j - report.kinetic_before_j,
            -heat,
            1e-30,
        );
        close(report.midpoint_impulse_work_j, -heat, 1e-30);
        close(report.dissipated_energy_j, heat, 1e-30);
        assert!(report.total_energy_residual_j.abs() <= controls().mechanical.energy_tolerance_j);
        assert!(
            norm(report.angular_momentum_change_nms)
                <= controls().mechanical.angular_momentum_tolerance_nms
        );
        let mut reversed = contacts;
        reversed.reverse();
        let retry = PreparedRigidAdhesiveIsland::new(
            &f.bodies,
            &f.properties,
            &reversed,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
        .advance(&CancelGate::new())
        .unwrap();
        assert_eq!(retry, candidate);
        let mut published = f.bodies.clone();
        let accepted = apply_rigid_adhesive_island(
            &mut published,
            &f.properties,
            &reversed,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(accepted, candidate.report);
        assert_eq!(&published[..2], candidate.bodies_after);
        assert_eq!(published[2], f.bodies[2]);
        for i in 0..2 {
            assert_eq!(published[i].pose(), f.bodies[i].pose());
        }
    }
}

#[test]
fn detached_sources_keep_identity_without_rows_recoil_or_heat_even_when_all_detached() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let attached_response = normal.evaluate_contact_load(0.0).unwrap();
    let gap_response = normal.evaluate(-0.5 * common::Z).unwrap();
    assert!(gap_response.normal_force_n < 0.0 && gap_response.contact_area_m2 == 0.0);
    let attached = shear
        .sliding_patch_dmt(normal.normal_state(), &attached_response)
        .unwrap();
    let gap = shear
        .sliding_patch_dmt(normal.normal_state(), &gap_response)
        .unwrap();
    let f = Fixture::new(attached_response.contact_radius_m);
    let pair = f.contacts(attached);
    let detached = RigidAdhesiveIslandContact {
        patch_id: "gap",
        contact: RigidAdhesiveContact {
            binding: RigidContactBinding {
                body_a: 1,
                body_b: 2,
                arm_a_body: Vec3::ZERO,
                arm_b_body: Vec3::ZERO,
            },
            source: gap,
            ..pair[0].contact
        },
    };
    let dt = RESOLVED_FRACTION * reference(f.radius, attached_response.contact_area_m2).3;
    let baseline = PreparedRigidAdhesiveIsland::new(
        &f.bodies,
        &f.properties,
        &pair,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&CancelGate::new())
    .unwrap();
    let prepared = PreparedRigidAdhesiveIsland::new(
        &f.bodies,
        &f.properties,
        &[pair[0], detached, pair[1]],
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(prepared.contacts().len(), 2);
    assert_eq!(prepared.detached_contacts().len(), 1);
    assert_eq!(prepared.body_indices(), &[0, 1]);
    assert_eq!(prepared.scaled_mobility_per_kg().len(), 36);
    let candidate = prepared.advance(&CancelGate::new()).unwrap();
    assert_eq!(candidate.bodies_after, baseline.bodies_after);
    assert_eq!(candidate.report.mechanical, baseline.report.mechanical);
    assert_eq!(candidate.report.detached_patch_ids, vec!["gap"]);

    let prepared = PreparedRigidAdhesiveIsland::new(
        &f.bodies,
        &f.properties,
        &[detached],
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(prepared.body_indices().is_empty());
    assert!(prepared.scaled_mobility_per_kg().is_empty());
    assert!(prepared.scaled_compliance_m_per_n().is_empty());
    let candidate = prepared.advance(&CancelGate::new()).unwrap();
    assert!(candidate.bodies_after.is_empty());
    assert!(candidate.report.mechanical.contacts.is_empty());
    assert_eq!(candidate.report.mechanical.dissipated_energy_j, 0.0);
    assert_eq!(candidate.report.mechanical.midpoint_impulse_work_j, 0.0);
    let mut bodies = f.bodies.clone();
    assert_eq!(
        apply_rigid_adhesive_island(
            &mut bodies,
            &f.properties,
            &[detached],
            dt,
            controls(),
            &CancelGate::new()
        )
        .unwrap(),
        candidate.report
    );
    assert_eq!(bodies, f.bodies);
}

#[test]
fn coarse_stopping_events_late_momentum_refusal_and_cancellation_never_publish_partial_bodies() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = Fixture::new(response.contact_radius_m);
    let contacts = f.contacts(source);
    let stop_time = reference(f.radius, response.contact_area_m2).3;
    let mut bodies = f.bodies.clone();
    // With the original bounded support work, this larger interval does not
    // resolve the strongly coupled opposite shear forces. It must not publish
    // an unconverged body update merely because the requested time is finite.
    assert!(matches!(
        apply_rigid_adhesive_island(
            &mut bodies,
            &f.properties,
            &contacts,
            0.25 * stop_time,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveIslandError::Native(_))
    ));
    assert_eq!(bodies, f.bodies);
    // An axial single disk removes the translation-conditioning limitation,
    // isolating the actual endpoint-event refusal from numerical work budgets.
    let mut axial = contacts[0];
    axial.contact.binding.arm_a_body = Vec3::ZERO;
    axial.contact.binding.arm_b_body = Vec3::ZERO;
    let spin_moment = (2.0 / 3.0) * STRENGTH * response.contact_area_m2 * f.radius;
    let axial_stop_time = SPIN / (spin_moment * (1.0 / IA + 1.0 / IB));
    assert!(matches!(
        apply_rigid_adhesive_island(
            &mut bodies,
            &f.properties,
            &[axial],
            1.5 * axial_stop_time,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveIslandError::UnresolvedSlipEvent { .. })
    ));
    assert_eq!(bodies, f.bodies);
    let prepared = PreparedRigidAdhesiveIsland::new(
        &bodies,
        &f.properties,
        &contacts,
        RESOLVED_FRACTION * stop_time,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let reference = prepared.advance(&CancelGate::new()).unwrap();
    let gate = CancelGate::new();
    gate.request();
    assert!(prepared.advance(&gate).is_err());
    assert!(
        apply_rigid_adhesive_island(
            &mut bodies,
            &f.properties,
            &contacts,
            RESOLVED_FRACTION * stop_time,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    // Admit nearby geometry but retain the much tighter physical angular-
    // momentum gate. Separated force application points create real net torque.
    let mut separated = contacts;
    separated[0].contact.binding.arm_a_body = separated[0]
        .contact
        .binding
        .arm_a_body
        .add(Vec3::new(1e-10, 0.0, 0.0));
    let mut permissive = controls();
    permissive.mechanical.geometry_tolerance_m = 2e-10;
    assert!(matches!(
        apply_rigid_adhesive_island(
            &mut bodies,
            &f.properties,
            &separated,
            RESOLVED_FRACTION * stop_time,
            permissive,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveIslandError::Audit(_))
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(prepared.advance(&CancelGate::new()).unwrap(), reference);
}

#[test]
fn held_normal_identity_and_whole_island_budgets_are_admitted_before_publication() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = Fixture::new(response.contact_radius_m);
    let contacts = f.contacts(source);
    let dt = RESOLVED_FRACTION * reference(f.radius, response.contact_area_m2).3;
    let mut duplicate = contacts;
    duplicate[1].patch_id = duplicate[0].patch_id;
    assert!(
        PreparedRigidAdhesiveIsland::new(
            &f.bodies,
            &f.properties,
            &duplicate,
            dt,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut separating = f.bodies.clone();
    separating[0] = RigidBodyState::new(
        separating[0].pose(),
        Vec3::new(0.0, 0.0, 1e-20),
        separating[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        PreparedRigidAdhesiveIsland::new(
            &separating,
            &f.properties,
            &contacts,
            dt,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    for limited in [
        RigidAdhesiveIslandControls {
            maximum_contacts: 1,
            ..controls()
        },
        RigidAdhesiveIslandControls {
            maximum_matrix_terms: 35,
            ..controls()
        },
        RigidAdhesiveIslandControls {
            maximum_evaluations: 1,
            ..controls()
        },
    ] {
        let mut bodies = f.bodies.clone();
        assert!(
            apply_rigid_adhesive_island(
                &mut bodies,
                &f.properties,
                &contacts,
                dt,
                limited,
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
    }
}
