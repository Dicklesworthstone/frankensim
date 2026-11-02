//! G1: independent finite-gap, area-crossing and full native recoil references.
//! G3: tangent gauge and original body indexing. G4: atomic physical refusals.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::dmt::{SourceDmtSphereContact, coupled::*};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::normal_motion::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;

const LENGTH: f64 = 1e-6;
const MASS: f64 = 1e-15;
const INERTIA: f64 = MASS * LENGTH * LENGTH;
const ADHESION: f64 = 0.05;
const KINETIC: f64 = 2e6;

fn close(a: f64, b: f64, absolute: f64) {
    common::close(a, b, absolute, 3e-9);
}
fn close_vec(a: Vec3, b: Vec3, absolute: f64) {
    close(a.x, b.x, absolute);
    close(a.y, b.y, absolute);
    close(a.z, b.z, absolute);
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn controls() -> RigidDmtControls {
    let mut contact = common::config();
    contact.coordinate_absolute_tolerance_m = 1e-23;
    contact.coordinate_relative_tolerance = 2e-12;
    contact.force_absolute_tolerance_n = 1e-20;
    contact.force_relative_tolerance = 1e-10;
    RigidDmtControls {
        contact,
        maximum_bodies: 4,
        geometry_tolerance_m: 1e-20,
        displacement_tolerance_m: 2e-21,
        velocity_tolerance_m_s: 1e-12,
        angular_velocity_tolerance_rad_s: 1e-6,
        linear_momentum_tolerance_ns: 1e-29,
        angular_momentum_tolerance_nms: 1e-34,
        energy_tolerance_j: 2e-26,
    }
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    reaction_arms: [Vec3; 2],
    frame: TangentFrame,
}
impl Fixture {
    fn contact<'a>(
        &self,
        normal: &'a SourceDmtSphereContact,
        shear: &'a fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear,
    ) -> RigidDmtContact<'a> {
        RigidDmtContact {
            binding: self.binding,
            frame: self.frame,
            reference_approach_m: 0.0,
            normal,
            shear,
        }
    }
    fn set_a_velocity(&mut self, velocity: Vec3, spin: f64) {
        self.bodies[0] = RigidBodyState::new(
            self.bodies[0].pose(),
            velocity.scale(self.properties[0].mass()),
            Vec3::new(
                0.0,
                0.0,
                self.properties[0].principal_inertia_body().z * spin,
            ),
        )
        .unwrap();
    }
    fn approach(&self) -> f64 {
        let a = self.bodies[0]
            .pose()
            .point_world_from_body(self.binding.arm_a_body)
            .unwrap();
        let b = self.bodies[1]
            .pose()
            .point_world_from_body(self.binding.arm_b_body)
            .unwrap();
        a.sub(b).dot(vector(self.frame.normal()))
    }
}
fn fixture(approach: f64, reaction_arms: [Vec3; 2]) -> Fixture {
    let properties = [
        MassProperties::new(
            2.0 * MASS,
            Vec3::ZERO,
            Vec3::new(0.8, 1.1, 1.4).scale(INERTIA),
        )
        .unwrap(),
        MassProperties::new(
            3.0 * MASS,
            Vec3::ZERO,
            Vec3::new(1.5, 1.7, 2.0).scale(INERTIA),
        )
        .unwrap(),
    ];
    let poses =
        reaction_arms.map(|arm| Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap());
    let bodies = poses.map(|pose| RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO).unwrap());
    Fixture {
        bodies,
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: reaction_arms[0].add(Vec3::new(0.0, 0.0, 0.5 * approach)),
            arm_b_body: reaction_arms[1].sub(Vec3::new(0.0, 0.0, 0.5 * approach)),
        },
        reaction_arms,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn axial(approach: f64) -> Fixture {
    fixture(
        approach,
        [Vec3::new(0.0, 0.0, LENGTH), Vec3::new(0.0, 0.0, -LENGTH)],
    )
}

// Independent inverse-inertia bilinear form; no production Gram or impulses.
fn mobility(f: &Fixture) -> [[f64; 3]; 3] {
    let axes = [f.frame.normal(), f.frame.longitudinal(), f.frame.lateral()].map(vector);
    let mut result = [[0.0; 3]; 3];
    for body in 0..2 {
        let inertia = f.properties[body].principal_inertia_body();
        let torques = axes.map(|axis| f.reaction_arms[body].cross(axis));
        for i in 0..3 {
            for j in 0..3 {
                let t = torques[j];
                result[i][j] += axes[i].dot(axes[j]) / f.properties[body].mass()
                    + torques[i].dot(Vec3::new(t.x / inertia.x, t.y / inertia.y, t.z / inertia.z));
            }
        }
    }
    result
}
fn kinetic(bodies: [RigidBodyState; 2], properties: [MassProperties; 2]) -> f64 {
    bodies
        .into_iter()
        .zip(properties)
        .map(|(body, p)| body.kinetic_energy(p).unwrap())
        .sum()
}
fn check_native(f: &Fixture, candidate: &RigidDmtCandidate) {
    let impulse = candidate.report.force_impulse_on_a_world_ns;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let expected = impulse.scale(sign);
        close_vec(
            candidate.bodies_after[i].linear_momentum_world(),
            f.bodies[i].linear_momentum_world().add(expected),
            1e-29,
        );
        close_vec(
            candidate.bodies_after[i].angular_momentum_body(),
            f.bodies[i]
                .angular_momentum_body()
                .add(f.reaction_arms[i].cross(expected)),
            1e-34,
        );
        assert_eq!(candidate.bodies_after[i].pose(), f.bodies[i].pose());
    }
    let delta_k = kinetic(candidate.bodies_after, f.properties) - kinetic(f.bodies, f.properties);
    close(candidate.report.kinetic_energy_change_j, delta_k, 2e-29);
    close(
        delta_k + candidate.report.normal_storage_change_j + candidate.report.dissipated_energy_j,
        0.0,
        2e-26,
    );
    assert!(norm(candidate.report.audit.linear_momentum_change_ns) < 1e-29);
    assert!(norm(candidate.report.audit.angular_momentum_change_nms) < 1e-34);
}

#[test]
fn finite_gap_attraction_drives_actual_recoil_and_geometry_owned_continuation() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let f = axial(-0.5 * common::Z);
    let dt = 5e-10;
    let gate = CancelGate::new();
    let contact = f.contact(&normal, &shear);
    let mut prepared =
        PreparedRigidDmtPair::new(&f.bodies, &f.properties, contact, dt, controls(), &gate)
            .unwrap();
    let start = prepared.initialize([0.0; 2]).unwrap();
    let candidate = prepared.advance(start, &gate).unwrap();
    let x0 = prepared.approach_m();
    let w = 1.0 / f.properties[0].mass() + 1.0 / f.properties[1].mass();
    let c = 0.5 * dt * dt * w;
    let stiffness = 2.0 * PI * normal.effective_radius_m() * ADHESION / common::Z;
    // In the finite gap, N = -k(z+x): solve the linear discrete gradient exactly.
    let dx = c * stiffness * (common::Z + x0) / (1.0 - 0.5 * c * stiffness);
    let force = -stiffness * (common::Z + x0 + 0.5 * dx);
    let storage = force * dx;
    close(candidate.state_after.position_m()[0], x0 + dx, 2e-21);
    close(candidate.report.applied_mean_force_n[0], force, 1e-18);
    close(
        candidate.report.relative_velocity_m_s[1][0],
        -dt * w * force,
        1e-12,
    );
    close(candidate.report.normal_storage_change_j, storage, 1e-28);
    close(candidate.report.kinetic_energy_change_j, -storage, 2e-26);
    assert!(candidate.report.kinetic_energy_change_j > 0.0);
    assert_eq!(candidate.report.contact.path.mean_real_area_m2, 0.0);
    assert_eq!(candidate.report.dissipated_energy_j, 0.0);
    assert_eq!(candidate.report.applied_mean_force_n[1..], [0.0; 2]);
    assert_eq!(candidate.state_after.regime(), DmtCoupledRegime::Released);
    assert!(candidate.bodies_after[0].linear_momentum_world().z > 0.0);
    assert!(candidate.bodies_after[1].linear_momentum_world().z < 0.0);
    close(
        candidate.report.marker_displacement_world_m[0]
            .sub(candidate.report.marker_displacement_world_m[1])
            .z,
        dx,
        2e-21,
    );
    check_native(&f, &candidate);

    let mut bodies = f.bodies;
    let mut history = start;
    let report = apply_rigid_dmt_contact(
        &mut bodies,
        &f.properties,
        contact,
        &mut history,
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(history, candidate.state_after);
    assert_eq!(bodies, candidate.bodies_after);
    let accepted = bodies;
    let accepted_history = history;
    assert!(matches!(
        apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            dt,
            controls(),
            &gate
        ),
        Err(RigidDmtError::GeometryMismatch { .. })
    ));
    assert_eq!(bodies, accepted);
    assert_eq!(history, accepted_history);
    // Pure axial translation needs no rotational transport. The existing pose
    // owner consumes the actual COM travel before the next normal phase.
    for i in 0..2 {
        let pose = Pose::new(
            bodies[i]
                .pose()
                .position_world()
                .add(report.body_displacement_world_m[i]),
            bodies[i].pose().orientation(),
        )
        .unwrap();
        bodies[i] = RigidBodyState::new(
            pose,
            bodies[i].linear_momentum_world(),
            bodies[i].angular_momentum_body(),
        )
        .unwrap();
    }
    apply_rigid_dmt_contact(
        &mut bodies,
        &f.properties,
        contact,
        &mut history,
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert!(history.position_m()[0] > accepted_history.position_m()[0]);
}

fn potential(normal: &SourceDmtSphereContact, x: f64) -> f64 {
    let pull = 2.0 * PI * normal.effective_radius_m() * ADHESION;
    let hertz = normal.evaluate(common::Z).unwrap().repulsive_force_n / common::Z.powf(1.5);
    if x >= 0.0 {
        0.4 * hertz * x.powf(2.5) - pull * (x + 0.5 * common::Z)
    } else {
        -0.5 * pull * (x + common::Z).max(0.0).powi(2) / common::Z
    }
}
// Manufacture an actual midpoint path from independent normal potential,
// analytic positive-approach area and explicit inverse-inertia response.
fn crossing(
    normal: &SourceDmtSphereContact,
    from: f64,
    to: f64,
    dt: f64,
) -> (Fixture, [f64; 3], f64) {
    let mut f = fixture(
        from,
        [
            Vec3::new(0.2, 0.1, 0.4).scale(LENGTH),
            Vec3::new(-0.12, 0.15, -0.3).scale(LENGTH),
        ],
    );
    let from = f.approach();
    // This tangent direction gives zero RELATIVE torsion, while each finite
    // body can recoil in all angular components. All three offdiagonals survive.
    let iz_a = f.properties[0].principal_inertia_body().z;
    let iz_b = f.properties[1].principal_inertia_body().z;
    let dx = f.reaction_arms[0].x / iz_a + f.reaction_arms[1].x / iz_b;
    let dy = f.reaction_arms[0].y / iz_a + f.reaction_arms[1].y / iz_b;
    let direction = [dx / dx.hypot(dy), dy / dx.hypot(dy)];
    let area = PI * normal.effective_radius_m() * (to.max(0.0).powi(2) - from.max(0.0).powi(2))
        / (2.0 * (to - from));
    let forces = [
        (potential(normal, to) - potential(normal, from)) / (to - from),
        KINETIC * area * direction[0],
        KINETIC * area * direction[1],
    ];
    let delta = [
        to - from,
        0.1 * common::Z * direction[0],
        0.1 * common::Z * direction[1],
    ];
    let w = mobility(&f);
    let velocity: [f64; 3] = core::array::from_fn(|i| {
        delta[i] / dt + 0.5 * dt * (0..3).map(|j| w[i][j] * forces[j]).sum::<f64>()
    });
    f.set_a_velocity(Vec3::new(velocity[1], velocity[2], velocity[0]), 0.0);
    (f, forces, area)
}

#[test]
fn attachment_and_detachment_use_integrated_real_area_and_all_native_cross_responses() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    for (from, to) in [
        (-0.2 * common::Z, 0.2 * common::Z),
        (0.2 * common::Z, -0.2 * common::Z),
    ] {
        let (f, expected_force, area) = crossing(&normal, from, to, dt);
        let w = mobility(&f);
        assert!(w[0][1].abs() > 1e12 && w[0][2].abs() > 1e12 && w[1][2].abs() > 1e12);
        let mut prepared = PreparedRigidDmtPair::new(
            &f.bodies,
            &f.properties,
            f.contact(&normal, &shear),
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let start = prepared.initialize([0.0; 2]).unwrap();
        let candidate = prepared.advance(start, &CancelGate::new()).unwrap();
        for i in 0..3 {
            close(
                candidate.report.applied_mean_force_n[i],
                expected_force[i],
                1e-18,
            );
            for j in 0..3 {
                close(candidate.report.mobility_per_kg[i][j], w[i][j], 2.0);
            }
        }
        close(candidate.state_after.position_m()[0], to, 2e-21);
        close(candidate.report.contact.path.mean_real_area_m2, area, 1e-27);
        close(
            candidate.report.dissipated_energy_j,
            KINETIC * area * 0.1 * common::Z,
            1e-28,
        );
        assert_eq!(candidate.report.contact.regime, DmtCoupledRegime::Sliding);
        assert_eq!(
            candidate.state_after.regime(),
            if to > 0.0 {
                DmtCoupledRegime::Sliding
            } else {
                DmtCoupledRegime::Released
            }
        );
        close(candidate.report.relative_spin_rad_s[1], 0.0, 1e-6);
        if to > 0.0 {
            close(
                area,
                0.25 * candidate.report.contact.path.contact.contact_area_m2,
                1e-27,
            );
        } else {
            assert_eq!(candidate.report.contact.path.contact.contact_area_m2, 0.0);
            assert!(candidate.report.dissipated_energy_j > 0.0);
        }
        check_native(&f, &candidate);
    }
}

#[test]
fn static_junction_balances_normal_shear_recoil_without_creating_friction_heat() {
    let (normal, shear) = common::sources(ADHESION, 2e8, KINETIC);
    let dt = 5e-10;
    let mut f = fixture(
        0.05 * common::Z,
        [
            Vec3::new(0.2, 0.0, 0.4).scale(LENGTH),
            Vec3::new(-0.12, 0.0, -0.3).scale(LENGTH),
        ],
    );
    let to = 0.1 * common::Z;
    let from = f.approach();
    let n = (potential(&normal, to) - potential(&normal, from)) / (to - from);
    let w = mobility(&f);
    let tangent = -w[1][0] * n / w[1][1];
    let normal_rate = (to - from) / dt + 0.5 * dt * (w[0][0] * n + w[0][1] * tangent);
    f.set_a_velocity(Vec3::new(0.0, 0.0, normal_rate), 0.0);
    let mut prepared = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        f.contact(&normal, &shear),
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let start = prepared.initialize([0.0; 2]).unwrap();
    assert_eq!(start.regime(), DmtCoupledRegime::Sticking);
    let candidate = prepared.advance(start, &CancelGate::new()).unwrap();
    assert_eq!(candidate.report.contact.regime, DmtCoupledRegime::Sticking);
    close(candidate.state_after.position_m()[0], to, 2e-21);
    close(candidate.report.applied_mean_force_n[1], tangent, 1e-18);
    assert!(tangent.abs() > candidate.report.contact.path.kinetic_capacity_n);
    assert!(tangent.abs() < candidate.report.contact.path.static_capacity_n);
    assert_eq!(candidate.report.dissipated_energy_j, 0.0);
    close(candidate.report.relative_velocity_m_s[1][1], 0.0, 1e-12);
    check_native(&f, &candidate);
}

#[test]
fn tangent_gauge_and_body_array_permutation_preserve_the_physical_pair() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    let (f, _, _) = crossing(&normal, -0.2 * common::Z, 0.2 * common::Z, dt);
    let base_contact = f.contact(&normal, &shear);
    let mut base = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        base_contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let state = base.initialize([0.0; 2]).unwrap();
    let expected = base.advance(state, &CancelGate::new()).unwrap();
    let angle = 0.73_f64;
    let rotated = RigidDmtContact {
        frame: TangentFrame::new([0.0, 0.0, 1.0], [angle.cos(), angle.sin(), 0.0]).unwrap(),
        ..base_contact
    };
    let mut prepared = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        rotated,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let rotated_state = prepared.initialize([0.0; 2]).unwrap();
    let actual = prepared.advance(rotated_state, &CancelGate::new()).unwrap();
    close_vec(
        actual.report.force_impulse_on_a_world_ns,
        expected.report.force_impulse_on_a_world_ns,
        1e-27,
    );
    close(
        actual.report.dissipated_energy_j,
        expected.report.dissipated_energy_j,
        1e-28,
    );
    for i in 0..2 {
        close_vec(
            actual.bodies_after[i].linear_momentum_world(),
            expected.bodies_after[i].linear_momentum_world(),
            1e-27,
        );
        close_vec(
            actual.bodies_after[i].angular_momentum_body(),
            expected.bodies_after[i].angular_momentum_body(),
            1e-33,
        );
    }
    let mut bodies = [f.bodies[1], f.bodies[0]];
    let properties = [f.properties[1], f.properties[0]];
    let permuted = RigidDmtContact {
        binding: RigidContactBinding {
            body_a: 1,
            body_b: 0,
            ..base_contact.binding
        },
        ..base_contact
    };
    let mut history = state;
    let report = apply_rigid_dmt_contact(
        &mut bodies,
        &properties,
        permuted,
        &mut history,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, [expected.bodies_after[1], expected.bodies_after[0]]);
    assert_eq!(history, expected.state_after);
    assert_eq!(report, expected.report);
}

#[test]
fn torsion_and_unresolved_sliding_capture_refuse_but_gap_spin_remains_free() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    for x in [0.2 * common::Z, -0.5 * common::Z] {
        let mut f = axial(x);
        f.set_a_velocity(Vec3::ZERO, 1000.0);
        let contact = f.contact(&normal, &shear);
        let prepared = PreparedRigidDmtPair::new(
            &f.bodies,
            &f.properties,
            contact,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let start = prepared.initialize([0.0; 2]).unwrap();
        let mut history = start;
        let mut bodies = f.bodies;
        let result = apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            dt,
            controls(),
            &CancelGate::new(),
        );
        if x > 0.0 {
            assert!(matches!(
                result,
                Err(RigidDmtError::UnsupportedTorsion { .. })
            ));
            assert_eq!(bodies, f.bodies);
            assert_eq!(history, start);
        } else {
            let report = result.unwrap();
            for rate in report.relative_spin_rad_s {
                close(rate, 1000.0, 1e-9);
            }
            assert_eq!(report.dissipated_energy_j, 0.0);
            assert_eq!(
                bodies[0].angular_momentum_body(),
                f.bodies[0].angular_momentum_body()
            );
        }
    }
    let mut f = axial(0.2 * common::Z);
    f.set_a_velocity(Vec3::new(1e-8, 0.0, 0.0), 0.0);
    let contact = f.contact(&normal, &shear);
    let prepared = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let start = prepared.initialize([0.0; 2]).unwrap();
    let mut history = start;
    let mut bodies = f.bodies;
    assert_eq!(start.regime(), DmtCoupledRegime::Sliding);
    assert!(matches!(
        apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            dt,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidDmtError::UnresolvedSlipEvent { .. })
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, start);
}

#[test]
fn foreign_history_geometry_cancellation_and_budget_refusals_leave_the_pair_retryable() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let (foreign_normal, foreign_shear) = common::sources(0.04, 3e6, KINETIC);
    let f = axial(-0.5 * common::Z);
    let dt = 5e-10;
    let contact = f.contact(&normal, &shear);
    let prepared = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let start = prepared.initialize([0.0; 2]).unwrap();
    let foreign = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        f.contact(&foreign_normal, &foreign_shear),
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .initialize([0.0; 2])
    .unwrap();
    let mut bodies = f.bodies;
    let mut history = foreign;
    assert!(matches!(
        apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            dt,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidDmtError::Contact(DmtCoupledError::StateMismatch))
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, foreign);
    history = start;
    let wrong_datum = RigidDmtContact {
        reference_approach_m: common::Z,
        ..contact
    };
    assert!(matches!(
        apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            wrong_datum,
            &mut history,
            dt,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidDmtError::GeometryMismatch { .. })
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            dt,
            controls(),
            &cancelled
        )
        .is_err()
    );
    let mut short_budget = controls();
    short_budget.contact.maximum_evaluations = 1;
    assert!(matches!(
        apply_rigid_dmt_contact(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            dt,
            short_budget,
            &CancelGate::new()
        ),
        Err(RigidDmtError::Contact(DmtCoupledError::Budget))
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, start);
    let mut retry = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let expected = retry.advance(start, &CancelGate::new()).unwrap();
    let report = apply_rigid_dmt_contact(
        &mut bodies,
        &f.properties,
        contact,
        &mut history,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(history, expected.state_after);
    assert_eq!(report, expected.report);
}
