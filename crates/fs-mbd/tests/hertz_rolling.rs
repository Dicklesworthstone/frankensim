//! Independent mechanics for simultaneous Hertz shear/spin and sphere bulk loss.
//! Synthetic Kelvin--Voigt data test the explicitly decoupled first-order
//! normal-pressure loss / local Coulomb shear approximation, not calibration.

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::{
    RigidHertzContact, RigidHertzControls, apply_rigid_hertz,
};
use fs_mbd::tangential_friction::hertz_rolling::*;
use fs_mbd::tangential_friction::rolling::{RigidRollingContact, RigidRollingHistory};
use fs_mbd::{MassProperties, PointKinematics, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::sliding_patch::contact::*;
use fs_tribo::{
    ApplicabilityRange, FrictionLaw, HertzSpherePlane, InputAuthority, InterfaceMedium,
    InterfaceSystemRef,
};

const R: f64 = 0.1;
const L: f64 = 0.05;
const TAU: f64 = 0.01;
const MIXED_SLIP: [f64; 2] = [0.0025, -0.0017];
const MIXED_ROLL: [f64; 2] = [0.03, -0.02];
const MIXED_SPIN: f64 = 0.25;
const DT: f64 = 0.0002;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn vector(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn array(x: Vec3) -> [f64; 3] {
    [x.x, x.y, x.z]
}
fn norm(x: Vec3) -> f64 {
    x.x.hypot(x.y).hypot(x.z)
}
fn z_cross(x: Vec3) -> Vec3 {
    Vec3::new(-x.y, x.x, 0.0)
}

fn rolling_law(tau: f64, source: &str) -> RollingLossLaw {
    let young = 1e6;
    let nu = 0.25;
    RollingLossLaw::ViscoelasticSphere(
        ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            source,
            InputAuthority::SyntheticFixture,
            ViscoelasticSphereParameters {
                sphere_surface: DeformingSphereSurface::SurfaceA,
                radius_m: R,
                young_modulus_pa: young,
                poisson_ratio: nu,
                density_kg_m3: 1000.0,
                shear_viscosity_pa_s: tau * young / (2.0 * (1.0 + nu)),
                bulk_viscosity_pa_s: tau * young / (3.0 * (1.0 - 2.0 * nu)),
            },
            QuasistaticSphereLimits {
                maximum_indentation_radius_ratio: 0.1,
                maximum_speed_wave_ratio: 0.1,
                maximum_relaxation_ratio: 0.1,
            },
            RollingLossApplicability::new(
                ApplicabilityRange::new(280.0, 320.0).unwrap(),
                ApplicabilityRange::new(0.0, 100.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap(),
    )
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
    hertz: HertzWrenchContact,
    rolling: RollingLossLaw,
    patch: RollingPatchReceipt,
    interface: InterfaceSystemRef,
    work: RollingWorkOwnership,
    radius: f64,
}
impl Fixture {
    fn contact(&self) -> RigidHertzRollingContact<'_> {
        RigidHertzRollingContact {
            binding: self.binding,
            frame: self.frame,
            hertz: &self.hertz,
            patch: &self.patch,
            interface: &self.interface,
            rolling: &self.rolling,
            temperature_k: 300.0,
            excitation_frequency_hz: 10.0,
            work: &self.work,
        }
    }
}
fn fixture(tau: f64, slip: [f64; 2], roll: [f64; 2], spin: f64, mu: [f64; 2]) -> Fixture {
    let interface = InterfaceSystemRef::new(
        "synthetic sphere/finite support",
        "original joint rolling history",
        "synthetic interface source",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let normal = HertzSpherePlane {
        effective_radius: R,
        reduced_modulus: 1e6 / (1.0 - 0.25_f64.powi(2)),
    }
    .response(&interface, 0.001)
    .unwrap();
    let patch = RollingPatchReceipt::new(
        "joint original patch",
        "original Hertz",
        "synthetic original normal source",
        InputAuthority::SyntheticFixture,
        normal.normal_force_n,
        PI * normal.contact_radius_m.powi(2),
        PatchCurvature::Principal {
            first_per_m: 1.0 / R,
            second_per_m: 1.0 / R,
        },
    )
    .unwrap();
    let mass = 4.0 * PI * 1000.0 * R.powi(3) / 3.0;
    let inertia = 0.4 * mass * R * R;
    let properties = [
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap(),
        MassProperties::new(7.0, Vec3::ZERO, Vec3::new(0.08, 0.13, 0.18)).unwrap(),
    ];
    let omega = Vec3::new(roll[0], roll[1], spin);
    let v = Vec3::new(slip[0], slip[1], 0.0).sub(z_cross(omega).scale(R));
    Fixture {
        bodies: [
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, R), UnitQuaternion::IDENTITY).unwrap(),
                v.scale(mass),
                properties[0].angular_momentum_body_checked(omega).unwrap(),
            )
            .unwrap(),
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, -L), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .unwrap(),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: Vec3::new(0.0, 0.0, -R),
            arm_b_body: Vec3::new(0.0, 0.0, L),
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        hertz: HertzWrenchContact::new(HertzWrenchParameters {
            normal_force_n: normal.normal_force_n,
            contact_radius_m: normal.contact_radius_m,
            static_mu: mu[0],
            kinetic_mu: mu[1],
        })
        .unwrap(),
        rolling: rolling_law(tau, "synthetic original viscosity"),
        patch,
        interface,
        work: RollingWorkOwnership::new(
            "joint original patch",
            "joint interval",
            "relative tangent rolling vector",
            RollingLossChannel::RollingHysteresis,
        )
        .unwrap(),
        radius: normal.contact_radius_m,
    }
}
fn controls() -> RigidHertzRollingControls {
    RigidHertzRollingControls {
        mechanical: RigidHertzControls {
            maximum_bodies: 4,
            geometry_tolerance_m: 2e-12,
            velocity_tolerance_m_s: 2e-11,
            angular_velocity_tolerance_rad_s: 2e-10,
            linear_momentum_tolerance_ns: 2e-11,
            angular_momentum_tolerance_nms: 2e-11,
            energy_tolerance_j: 2e-11,
            contact: HertzWrenchConfig {
                maximum_static_iterations: 96,
                maximum_iterations: 256,
                maximum_evaluations: 512,
                force_tolerance_n: 2e-9,
                moment_tolerance_nm: 2e-10,
                displacement_tolerance_m: 2e-13,
                energy_tolerance_j: 2e-12,
                relative_matrix_tolerance: 1e-13,
                maximum_force_n: 20.0,
                maximum_spin_moment_nm: 2.0,
            },
        },
        maximum_rolling_moment_nm: 1.0,
        rolling_moment_tolerance_nm: 2e-11,
    }
}
fn advance(
    f: &Fixture,
    bodies: &mut [RigidBodyState],
    history: &mut RigidRollingHistory,
    dt: f64,
) -> RigidHertzRollingReport {
    apply_rigid_hertz_rolling(
        bodies,
        &f.properties,
        history,
        f.contact(),
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn points(f: &Fixture, bodies: &[RigidBodyState; 2]) -> [PointKinematics; 2] {
    [
        bodies[0]
            .point_kinematics(f.properties[0], f.binding.arm_a_body)
            .unwrap(),
        bodies[1]
            .point_kinematics(f.properties[1], f.binding.arm_b_body)
            .unwrap(),
    ]
}
fn relative_point(p: [PointKinematics; 2]) -> Vec3 {
    p[0].point_velocity_world.sub(p[1].point_velocity_world)
}
fn relative_omega(p: [PointKinematics; 2]) -> Vec3 {
    p[0].angular_velocity_world.sub(p[1].angular_velocity_world)
}
fn kinetic(bodies: &[RigidBodyState; 2], properties: &[MassProperties; 2]) -> f64 {
    bodies
        .iter()
        .zip(properties)
        .map(|(b, p)| b.kinetic_energy(*p).unwrap())
        .sum()
}

#[test]
fn planar_static_sliding_and_traction_induced_rolling_match_independent_equations() {
    let h = 0.002;
    for (slip, q0, static_mode) in [(0.0, 0.03, true), (0.1, 0.03, false), (0.1, 0.0, false)] {
        let f = fixture(TAU, [slip, 0.0], [0.0, q0], 0.0, [0.3, 0.2]);
        let ia = f.properties[0].principal_inertia_body().y;
        let ib = f.properties[1].principal_inertia_body().y;
        let alpha =
            1.0 / f.properties[0].mass() + 1.0 / f.properties[1].mass() + R * R / ia + L * L / ib;
        let beta = -R / ia + L / ib;
        let gamma = 1.0 / ia + 1.0 / ib;
        let c = TAU * f.patch.normal_load_n() * R;
        let denominator = 1.0 + 0.5 * h * gamma * c;
        let force = if static_mode {
            (beta * c * q0 / denominator - 2.0 * slip / h)
                / (alpha - 0.5 * h * beta * beta * c / denominator)
        } else {
            -0.2 * f.patch.normal_load_n()
        };
        let qmid = (q0 + 0.5 * h * beta * force) / denominator;
        let moment = -c * qmid;
        let umid = slip + 0.5 * h * (alpha * force + beta * moment);
        let heat_surface = -h * force * umid;
        let heat_bulk = c * h * qmid * qmid;
        let mut bodies = f.bodies;
        let mut history = RigidRollingHistory::new();
        let report = advance(&f, &mut bodies, &mut history, h);
        assert_eq!(
            report.hertz.mode,
            if static_mode {
                HertzWrenchMode::Sticking
            } else {
                HertzWrenchMode::Sliding
            }
        );
        close(report.hertz.force_n[0], force, 3e-9);
        close(report.hertz.force_n[1], 0.0, 2e-12);
        close(report.hertz.spin_moment_nm, 0.0, 2e-12);
        close(
            bodies[0].linear_momentum_world().x,
            f.bodies[0].linear_momentum_world().x + h * force,
            1e-11,
        );
        close(bodies[1].linear_momentum_world().x, -h * force, 1e-11);
        close(
            bodies[0].angular_momentum_body().y,
            f.bodies[0].angular_momentum_body().y + h * (-R * force + moment),
            2e-12,
        );
        close(
            bodies[1].angular_momentum_body().y,
            h * (-L * force - moment),
            2e-12,
        );
        close(
            relative_point(points(&f, &bodies)).x,
            slip + h * (alpha * force + beta * moment),
            1e-11,
        );
        close(
            relative_omega(points(&f, &bodies)).y,
            q0 + h * (beta * force + gamma * moment),
            2e-10,
        );
        close(report.surface_heat_j, heat_surface, 2e-12);
        close(report.bulk_heat_j, heat_bulk, 3e-15);
        close(history.state().cumulative_heat_j(), heat_bulk, 3e-15);
        if q0 == 0.0 {
            assert!(qmid > 0.0);
            assert!(report.bulk_heat_j > 1e-10);
            assert!(bodies[0].angular_momentum_body().y > 0.0);
        }
        if static_mode {
            assert_eq!(report.surface_heat_j, 0.0);
            close(relative_point(points(&f, &bodies)).x, 0.0, 2e-11);
        }
        close(
            kinetic(&bodies, &f.properties) - kinetic(&f.bodies, &f.properties)
                + heat_surface
                + heat_bulk,
            0.0,
            2e-12,
        );
    }
}

#[test]
fn zero_viscosity_recovers_the_original_hertz_pair_without_extra_heat() {
    let f = fixture(0.0, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let initial_history = RigidRollingHistory::new();
    let gate = CancelGate::new();
    let prepared = PreparedRigidHertzRollingPair::new(
        &f.bodies,
        &f.properties,
        &initial_history,
        f.contact(),
        DT,
        controls(),
        &gate,
    )
    .unwrap();
    let candidate = prepared.advance(&f.hertz, &gate).unwrap();
    let repeated = prepared.advance(&f.hertz, &gate).unwrap();
    assert_eq!(candidate.bodies_after, repeated.bodies_after);
    assert_eq!(candidate.rolling_history, repeated.rolling_history);
    assert_eq!(candidate.report, repeated.report);
    let warmer = prepared
        .advance_at_temperature(&f.hertz, 310.0, &gate)
        .unwrap();
    assert_eq!(warmer.bodies_after, candidate.bodies_after);
    assert_eq!(warmer.rolling_history, candidate.rolling_history);
    assert_eq!(warmer.report.rolling_temperature_k, 310.0);
    assert!(
        prepared
            .advance_at_temperature(&f.hertz, 330.0, &gate)
            .is_err()
    );
    assert_eq!(initial_history, RigidRollingHistory::new());
    let mut original = f.bodies;
    let expected = apply_rigid_hertz(
        &mut original,
        &f.properties,
        RigidHertzContact {
            binding: f.binding,
            frame: f.frame,
            law: &f.hertz,
        },
        DT,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap();
    let mut joint = f.bodies;
    let mut history = RigidRollingHistory::new();
    let actual = advance(&f, &mut joint, &mut history, DT);
    assert_eq!(candidate.bodies_after, joint);
    assert_eq!(candidate.rolling_history, history);
    assert_eq!(candidate.report, actual);
    for i in 0..2 {
        close_vec(
            joint[i].linear_momentum_world(),
            original[i].linear_momentum_world(),
            2e-12,
        );
        close_vec(
            joint[i].angular_momentum_body(),
            original[i].angular_momentum_body(),
            2e-13,
        );
    }
    close(actual.surface_heat_j, expected.dissipated_energy_j, 3e-13);
    close(
        actual.hertz.spin_moment_nm,
        expected.contact.spin_moment_nm,
        2e-10,
    );
    assert_eq!(actual.bulk_heat_j, 0.0);
    assert_eq!(history.state().cumulative_heat_j(), 0.0);
    assert_eq!(actual.rolling.storage_change_j, 0.0);
}

#[test]
fn zero_friction_is_an_unconstrained_viscous_couple_with_both_bodies_recoiling() {
    let f = fixture(TAU, [0.0; 2], MIXED_ROLL, MIXED_SPIN, [0.0; 2]);
    let h = 0.05;
    let ia = f.properties[0].principal_inertia_body().x;
    let c = TAU * f.patch.normal_load_n() * R;
    let mut expected_heat = 0.0;
    let mut impulse = [0.0; 2];
    let mut after = [0.0; 2];
    for (i, ib) in [0.08, 0.13].into_iter().enumerate() {
        let gamma = 1.0 / ia + 1.0 / ib;
        let midpoint = MIXED_ROLL[i] / (1.0 + 0.5 * h * c * gamma);
        impulse[i] = -h * c * midpoint;
        after[i] = MIXED_ROLL[i] + gamma * impulse[i];
        expected_heat += 0.5 * (MIXED_ROLL[i].powi(2) - after[i].powi(2)) / gamma;
    }
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let report = advance(&f, &mut bodies, &mut history, h);
    assert_eq!(report.hertz.mode, HertzWrenchMode::Frictionless);
    assert_eq!(report.surface_heat_j, 0.0);
    for i in 0..2 {
        assert_eq!(
            bodies[i].linear_momentum_world(),
            f.bodies[i].linear_momentum_world()
        );
    }
    let j = Vec3::new(impulse[0], impulse[1], 0.0);
    close_vec(
        bodies[0].angular_momentum_body(),
        f.bodies[0].angular_momentum_body().add(j),
        3e-17,
    );
    close_vec(bodies[1].angular_momentum_body(), j.scale(-1.0), 3e-17);
    let omega = relative_omega(points(&f, &bodies));
    close(omega.x, after[0], 3e-15);
    close(omega.y, after[1], 3e-15);
    close(omega.z, MIXED_SPIN, 3e-15);
    close(report.bulk_heat_j, expected_heat, 3e-17);
    assert!(norm(relative_point(points(&f, &bodies))) > 1e-6);
    close(
        kinetic(&bodies, &f.properties) - kinetic(&f.bodies, &f.properties) + expected_heat,
        0.0,
        3e-17,
    );
}

// Direct Hertz pressure integration, independent of all analytic shape factors
// and native support evaluations. z=sqrt(1-r²/a²) gives pressure-area weight
// (3N/2pi)*z² dz dtheta. Simpson in z, periodic trapezoidal rule in theta.
fn surface_reference(
    u: [f64; 2],
    omega: f64,
    normal: f64,
    radius: f64,
    mu: f64,
    h: f64,
    radial: usize,
) -> ([f64; 2], f64, f64) {
    let angular = 1024;
    let angles: Vec<_> = (0..angular)
        .map(|i| {
            let theta = 2.0 * PI * (i as f64 + 0.5) / angular as f64;
            (theta.cos(), theta.sin())
        })
        .collect();
    let mut sum = [0.0; 4];
    for j in 0..=radial {
        let z = j as f64 / radial as f64;
        let r = radius * (1.0 - z * z).sqrt();
        let weight = z
            * z
            * if j == 0 || j == radial {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
        let mut ring = [0.0; 4];
        for &(cos, sin) in &angles {
            let x = r * cos;
            let y = r * sin;
            let vx = u[0] - omega * y;
            let vy = u[1] + omega * x;
            let speed = vx.hypot(vy);
            assert!(speed > 0.0);
            ring[0] += vx / speed;
            ring[1] += vy / speed;
            ring[2] += (x * vy - y * vx) / speed;
            ring[3] += speed;
        }
        for k in 0..4 {
            sum[k] += weight * ring[k] / angular as f64;
        }
    }
    let factor = normal * mu / radial as f64;
    (
        [-factor * sum[0], -factor * sum[1]],
        -factor * sum[2],
        factor * sum[3] * h,
    )
}

#[test]
fn mixed_shear_spin_and_roll_close_actual_body_work_and_independent_pressure_heat() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let report = advance(&f, &mut bodies, &mut history, DT);
    let before = points(&f, &f.bodies);
    let after = points(&f, &bodies);
    let u = relative_point(before).add(relative_point(after)).scale(0.5);
    let omega = relative_omega(before).add(relative_omega(after)).scale(0.5);
    let impulse = bodies[0]
        .linear_momentum_world()
        .sub(f.bodies[0].linear_momentum_world());
    let angular_change = f.bodies[0].pose().orientation().rotate_body_to_world(
        bodies[0]
            .angular_momentum_body()
            .sub(f.bodies[0].angular_momentum_body()),
    );
    let couple = angular_change.sub(before[0].arm_world.cross(impulse));
    let c = TAU * f.patch.normal_load_n() * R;
    let bulk_heat = c * DT * (omega.x * omega.x + omega.y * omega.y);
    close_vec(
        Vec3::new(couple.x, couple.y, 0.0),
        Vec3::new(-c * DT * omega.x, -c * DT * omega.y, 0.0),
        2e-16,
    );
    let (force, moment, surface_heat) = surface_reference(
        [u.x, u.y],
        omega.z,
        f.patch.normal_load_n(),
        f.radius,
        0.2,
        DT,
        256,
    );
    let (coarse_force, coarse_moment, coarse_heat) = surface_reference(
        [u.x, u.y],
        omega.z,
        f.patch.normal_load_n(),
        f.radius,
        0.2,
        DT,
        128,
    );
    close(force[0], coarse_force[0], 2e-8);
    close(force[1], coarse_force[1], 2e-8);
    close(moment, coarse_moment, 2e-10);
    close(surface_heat, coarse_heat, 2e-12);
    close(impulse.x / DT, force[0], 3e-8);
    close(impulse.y / DT, force[1], 3e-8);
    close(couple.z / DT, moment, 3e-10);
    close(report.surface_heat_j, surface_heat, 3e-12);
    close(report.bulk_heat_j, bulk_heat, 3e-17);
    assert!(report.surface_heat_j > 1e-8);
    assert!(report.bulk_heat_j > 1e-10);
    assert_eq!(report.hertz.mode, HertzWrenchMode::Sliding);
    let work = impulse.dot(u) + couple.dot(omega);
    let delta_k = kinetic(&bodies, &f.properties) - kinetic(&f.bodies, &f.properties);
    close(delta_k, work, 2e-16);
    close(work + surface_heat + bulk_heat, 0.0, 3e-12);
    close(report.midpoint_impulse_work_j, work, 2e-16);
    close(history.state().cumulative_heat_j(), bulk_heat, 3e-17);
    close(report.rolling.dissipation.total_heat_j, bulk_heat, 3e-17);
    assert_eq!(report.rolling.storage_change_j, 0.0);
    assert_eq!(report.rolling.generalized_work.ownership, f.work);
    for i in 0..2 {
        assert_eq!(bodies[i].pose(), f.bodies[i].pose());
    }
    // Joint heat differs from a sequential Hertz then rolling postprocess: the
    // same accepted midpoint velocities drive BOTH physical losses.
    let scaled = [u.x, u.y, f.radius * omega.z, R * omega.x, R * omega.y];
    for i in 0..5 {
        close(
            report.actual_midpoint_coordinates_m[i],
            DT * scaled[i],
            2e-15,
        );
    }
}

fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [aw, ax, ay, az] = a.components();
    let [bw, bx, by, bz] = b.components();
    UnitQuaternion::new(
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    )
    .unwrap()
}
fn rotate_body(body: RigidBodyState, rotation: UnitQuaternion) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            rotation.rotate_body_to_world(body.pose().position_world()),
            compose(rotation, body.pose().orientation()),
        )
        .unwrap(),
        rotation.rotate_body_to_world(body.linear_momentum_world()),
        body.angular_momentum_body(),
    )
    .unwrap()
}

#[test]
fn world_rotation_tangent_gauge_and_original_body_indices_preserve_the_joint_solution() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let mut baseline = f.bodies;
    let mut baseline_history = RigidRollingHistory::new();
    let expected = advance(&f, &mut baseline, &mut baseline_history, DT);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.73).unwrap();
    let spectator = RigidBodyState::new(
        Pose::new(Vec3::new(9.0, -4.0, 7.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(1e18, -3e18, 2e18),
        Vec3::new(3e17, -1e17, 2e17),
    )
    .unwrap();
    let mut bodies = [
        spectator,
        rotate_body(f.bodies[1], rotation),
        rotate_body(f.bodies[0], rotation),
    ];
    let properties = [
        MassProperties::new(5.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap(),
        f.properties[1],
        f.properties[0],
    ];
    let frame = TangentFrame::new(
        array(rotation.rotate_body_to_world(vector(f.frame.normal()))),
        array(rotation.rotate_body_to_world(vector(f.frame.longitudinal()))),
    )
    .unwrap()
    .rotated(0.61)
    .unwrap();
    let contact = RigidHertzRollingContact {
        binding: RigidContactBinding {
            body_a: 2,
            body_b: 1,
            ..f.binding
        },
        frame,
        ..f.contact()
    };
    let mut history = RigidRollingHistory::new();
    let report = apply_rigid_hertz_rolling(
        &mut bodies,
        &properties,
        &mut history,
        contact,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    for (index, original) in [(2, baseline[0]), (1, baseline[1])] {
        let expected_body = rotate_body(original, rotation);
        close_vec(
            bodies[index].linear_momentum_world(),
            expected_body.linear_momentum_world(),
            3e-11,
        );
        close_vec(
            bodies[index].angular_momentum_body(),
            expected_body.angular_momentum_body(),
            3e-12,
        );
        assert_eq!(bodies[index].pose(), expected_body.pose());
    }
    assert_eq!(bodies[0], spectator);
    assert_eq!(report.deforming_body, 2);
    close(report.surface_heat_j, expected.surface_heat_j, 3e-12);
    close(report.bulk_heat_j, expected.bulk_heat_j, 3e-17);
    close(
        history.state().cumulative_heat_j(),
        baseline_history.state().cumulative_heat_j(),
        3e-17,
    );
    close(report.kinetic_before_j, expected.kinetic_before_j, 2e-16);
    close(report.kinetic_after_j, expected.kinetic_after_j, 2e-16);
}

#[test]
fn source_geometry_approximation_and_unresolved_event_refusals_preserve_accepted_states() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let invalid_normal = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 1.1 * f.patch.normal_load_n(),
        contact_radius_m: f.radius,
        static_mu: 0.3,
        kinetic_mu: 0.2,
    })
    .unwrap();
    for contact in [
        RigidHertzRollingContact {
            temperature_k: 330.0,
            ..f.contact()
        },
        RigidHertzRollingContact {
            excitation_frequency_hz: 101.0,
            ..f.contact()
        },
        RigidHertzRollingContact {
            hertz: &invalid_normal,
            ..f.contact()
        },
        RigidHertzRollingContact {
            binding: RigidContactBinding {
                arm_a_body: Vec3::new(0.0, 0.0, -1.1 * R),
                ..f.binding
            },
            ..f.contact()
        },
    ] {
        assert!(
            apply_rigid_hertz_rolling(
                &mut bodies,
                &f.properties,
                &mut history,
                contact,
                DT,
                controls(),
                &CancelGate::new(),
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
        assert_eq!(history, RigidRollingHistory::new());
    }
    let too_fast = fixture(TAU, MIXED_SLIP, [0.3, 0.0], MIXED_SPIN, [0.3, 0.2]);
    let mut fast_bodies = too_fast.bodies;
    assert!(
        apply_rigid_hertz_rolling(
            &mut fast_bodies,
            &too_fast.properties,
            &mut history,
            too_fast.contact(),
            DT,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(fast_bodies, too_fast.bodies);
    // Zero initial rolling and an admitted midpoint do not admit the faster
    // endpoint induced by friction. The original small-parameter box must be
    // checked there too; the same initial state passed at h=0.002 above.
    let driven = fixture(TAU, [0.1, 0.0], [0.0; 2], 0.0, [0.3, 0.2]);
    let mut driven_bodies = driven.bodies;
    assert!(
        apply_rigid_hertz_rolling(
            &mut driven_bodies,
            &driven.properties,
            &mut history,
            driven.contact(),
            0.01,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(driven_bodies, driven.bodies);
    let coarse = fixture(0.0, [0.01, 0.0], [0.0; 2], 0.0, [0.6, 0.2]);
    let mut coarse_bodies = coarse.bodies;
    assert!(
        apply_rigid_hertz_rolling(
            &mut coarse_bodies,
            &coarse.properties,
            &mut history,
            coarse.contact(),
            0.1,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(coarse_bodies, coarse.bodies);
    // Resolving the physical interval, with the same model and initial state,
    // admits the phase without inventing endpoint arrest or heat.
    let mut refined = coarse.bodies;
    advance(&coarse, &mut refined, &mut history, 0.0001);
    assert!(relative_point(points(&coarse, &refined)).x > 0.0);
    // A released normal patch has no surface or bulk loss. Its inactive Hertz
    // length scale stays positive without manufacturing a nonzero normal area.
    let released_patch = RollingPatchReceipt::new(
        "joint original patch",
        "original Hertz",
        "synthetic original normal source",
        InputAuthority::SyntheticFixture,
        0.0,
        0.0,
        PatchCurvature::Principal {
            first_per_m: 1.0 / R,
            second_per_m: 1.0 / R,
        },
    )
    .unwrap();
    let released_hertz = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 0.0,
        contact_radius_m: f.radius,
        static_mu: 0.3,
        kinetic_mu: 0.2,
    })
    .unwrap();
    let mut released = f.bodies;
    let released_report = apply_rigid_hertz_rolling(
        &mut released,
        &f.properties,
        &mut RigidRollingHistory::new(),
        RigidHertzRollingContact {
            hertz: &released_hertz,
            patch: &released_patch,
            ..f.contact()
        },
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(released, f.bodies);
    assert_eq!(released_report.surface_heat_j, 0.0);
    assert_eq!(released_report.bulk_heat_j, 0.0);
}

#[test]
fn loaded_native_checkpoint_late_audit_source_budget_and_cancel_have_exact_retry() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let mut accepted = f.bodies;
    let mut history = RigidRollingHistory::new();
    let first = advance(&f, &mut accepted, &mut history, DT);
    let saved_history = history.clone();
    assert!(history.checkpoint().is_some());
    let static_friction = FrictionLaw::Coulomb {
        static_mu: 0.3,
        kinetic_mu: 0.2,
    };
    let restored = RigidRollingHistory::from_checkpoint(
        RigidRollingContact {
            binding: f.binding,
            frame: f.frame,
            patch: &f.patch,
            interface: &f.interface,
            law: &f.rolling,
            static_friction: &static_friction,
            temperature_k: 300.0,
            excitation_frequency_hz: 10.0,
            work: &f.work,
        },
        history.checkpoint().unwrap().clone(),
    )
    .unwrap();
    assert_eq!(restored, history);
    let mut expected = accepted;
    let mut expected_history = restored;
    let expected_report = advance(&f, &mut expected, &mut expected_history, DT);

    let mut shifted = accepted;
    shifted[1] = RigidBodyState::new(
        Pose::new(
            shifted[1]
                .pose()
                .position_world()
                .add(Vec3::new(1e-8, 0.0, 0.0)),
            shifted[1].pose().orientation(),
        )
        .unwrap(),
        shifted[1].linear_momentum_world(),
        shifted[1].angular_momentum_body(),
    )
    .unwrap();
    let shifted_before = shifted;
    let mut late_controls = controls();
    late_controls.mechanical.geometry_tolerance_m = 1e-7;
    late_controls.mechanical.angular_momentum_tolerance_nms = 1e-20;
    assert!(
        apply_rigid_hertz_rolling(
            &mut shifted,
            &f.properties,
            &mut history,
            f.contact(),
            DT,
            late_controls,
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(shifted, shifted_before);
    assert_eq!(history, saved_history);

    let replacement = rolling_law(
        TAU,
        "different viscosity source cannot replace loaded history",
    );
    let mut bodies = accepted;
    assert!(
        apply_rigid_hertz_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            RigidHertzRollingContact {
                rolling: &replacement,
                ..f.contact()
            },
            DT,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(bodies, accepted);
    assert_eq!(history, saved_history);
    let limited = RigidHertzRollingControls {
        maximum_rolling_moment_nm: 1e-10,
        ..controls()
    };
    assert!(
        apply_rigid_hertz_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            DT,
            limited,
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(bodies, accepted);
    assert_eq!(history, saved_history);
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_hertz_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            DT,
            controls(),
            &gate,
        )
        .is_err()
    );
    assert_eq!(bodies, accepted);
    assert_eq!(history, saved_history);
    let report = advance(&f, &mut bodies, &mut history, DT);
    assert_eq!(bodies, expected);
    assert_eq!(history, expected_history);
    assert_eq!(report, expected_report);
    close(
        history.state().cumulative_heat_j(),
        first.bulk_heat_j + report.bulk_heat_j,
        2e-18,
    );
    assert_eq!(history.state().accepted_steps(), 2);
    assert!(history.state().cumulative_heat_j() < first.surface_heat_j + report.surface_heat_j);
}
