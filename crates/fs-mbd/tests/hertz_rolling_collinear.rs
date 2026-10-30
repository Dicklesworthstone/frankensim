//! Physical stop/restart and exact two-coordinate native body references.
//! These fixtures exercise held geometry and synthetic source material only.

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::RigidHertzControls;
use fs_mbd::tangential_friction::hertz_rolling::collinear::*;
use fs_mbd::tangential_friction::hertz_rolling::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::sliding_patch::contact::*;
use fs_tribo::{
    ApplicabilityRange, HertzSpherePlane, InputAuthority, InterfaceMedium, InterfaceSystemRef,
};

const R: f64 = 0.1;
const L: f64 = 0.06;
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}; tol={tol:e}");
}
fn close_vec(a: Vec3, b: Vec3, tol: f64) {
    close(a.x, b.x, tol);
    close(a.y, b.y, tol);
    close(a.z, b.z, tol);
}
fn controls() -> RigidHertzRollingCollinearControls {
    RigidHertzRollingCollinearControls {
        mechanical: RigidHertzRollingControls {
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
        },
        maximum_phases: 16,
        maximum_root_iterations: 128,
        time_tolerance_s: 1e-13,
    }
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
    history: RigidHertzRollingCollinearHistory,
    tau: f64,
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
    fn prepare(&self, dt: f64) -> PreparedRigidHertzRollingCollinear<'_> {
        PreparedRigidHertzRollingCollinear::new(
            &self.bodies,
            &self.properties,
            &self.history,
            self.contact(),
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
    }
    fn run(&self, dt: f64) -> RigidHertzRollingCollinearCandidate {
        self.prepare(dt)
            .advance(&self.hertz, &CancelGate::new())
            .unwrap()
    }
    fn response(&self) -> [[f64; 2]; 2] {
        let mass = self.properties[0].mass();
        let inertia = 0.4 * mass * R * R;
        let a = 1.0 / mass + 1.0 / 7.0 + R * R / inertia + L * L / 0.13;
        let b = R * (-R / inertia + L / 0.13);
        [[a, b], [b, R * R * (1.0 / inertia + 1.0 / 0.13)]]
    }
}
fn fixture(tau: f64, mu: [f64; 2], slip: f64, rolling_rate: f64) -> Fixture {
    let interface = InterfaceSystemRef::new(
        "sphere/finite plane",
        "pure shear and material rolling",
        "synthetic source",
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
        "collinear original patch",
        "Hertz normal",
        "synthetic normal source",
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
    let bodies = [
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, R), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(mass * (slip + R * rolling_rate), 0.0, 0.0),
            Vec3::new(0.0, inertia * rolling_rate, 0.0),
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, -L), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap(),
    ];
    Fixture {
        bodies,
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
        rolling: RollingLossLaw::ViscoelasticSphere(
            ViscoelasticSphereRollingCard::new(
                VISCOELASTIC_SPHERE_LAW_ID,
                "synthetic Kelvin-Voigt sphere",
                InputAuthority::SyntheticFixture,
                ViscoelasticSphereParameters {
                    sphere_surface: DeformingSphereSurface::SurfaceA,
                    radius_m: R,
                    young_modulus_pa: 1e6,
                    poisson_ratio: 0.25,
                    density_kg_m3: 1000.0,
                    shear_viscosity_pa_s: tau * 1e6 / 2.5,
                    bulk_viscosity_pa_s: tau * 1e6 / 1.5,
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
        ),
        patch,
        interface,
        work: RollingWorkOwnership::new(
            "collinear original patch",
            "trajectory interval",
            "relative tangent rolling vector",
            RollingLossChannel::RollingHysteresis,
        )
        .unwrap(),
        history: RigidHertzRollingCollinearHistory::default(),
        tau,
    }
}
fn audit(f: &Fixture, r: &RigidHertzRollingCollinearCandidate) {
    let before = f
        .bodies
        .iter()
        .zip(&f.properties)
        .map(|(b, p)| b.kinetic_energy(*p).unwrap())
        .sum::<f64>();
    let after = r
        .bodies_after
        .iter()
        .zip(&f.properties)
        .map(|(b, p)| b.kinetic_energy(*p).unwrap())
        .sum::<f64>();
    close(
        after - before + r.report.surface_heat_j + r.report.bulk_heat_j,
        0.0,
        2e-11,
    );
    close(before, r.report.kinetic_before_j, 2e-16);
    close(after, r.report.kinetic_after_j, 2e-16);
    close_vec(r.report.linear_momentum_change_ns, Vec3::ZERO, 2e-12);
    close_vec(r.report.angular_momentum_change_nms, Vec3::ZERO, 2e-12);
    close(
        r.history.rolling_history().state().cumulative_heat_j()
            - f.history.rolling_history().state().cumulative_heat_j(),
        r.report.bulk_heat_j,
        2e-16,
    );
    for i in 0..2 {
        assert_eq!(
            f.bodies[i].pose(),
            r.bodies_after[i].pose(),
            "this owner must not pretend to advance geometry"
        );
    }
}

#[test]
fn moving_slip_uses_kinetic_force_until_exact_stop_even_with_large_static_capacity() {
    let f = fixture(0.0, [0.8, 0.1], 0.004, 0.0);
    let w = f.response();
    let force = 0.1 * f.patch.normal_load_n();
    let stop = 0.004 / (w[0][0] * force);
    let r = f.run(0.2);
    audit(&f, &r);
    assert_eq!(r.report.phases.len(), 2);
    let first = &r.report.phases[0];
    assert_eq!(first.mode, RigidHertzRollingCollinearMode::SlidingPositive);
    assert!(first.reached_stop);
    close(first.duration_s, stop, 2e-15);
    close(first.maximum_shear_force_n, force, 2e-15);
    close(
        r.report.surface_heat_j,
        0.5 * 0.004_f64.powi(2) / w[0][0],
        2e-17,
    );
    close(r.report.bulk_heat_j, 0.0, 0.0);
    assert_eq!(
        r.report.phases[1].mode,
        RigidHertzRollingCollinearMode::Sticking
    );
    assert!(r.history.is_stopped());
    close(r.report.force_impulse_world_ns.x, -0.004 / w[0][0], 2e-15);
    close(
        first.rolling_rate_after_rad_s,
        (w[1][0] * (-0.004 / w[0][0])) / R,
        2e-14,
    );
}

#[test]
fn nonzero_tiny_slip_is_integrated_and_never_reclassified_as_initial_stick() {
    let f = fixture(0.0, [0.8, 0.1], 1e-20, 0.0);
    let r = f.run(0.001);
    assert_eq!(
        r.report.phases[0].mode,
        RigidHertzRollingCollinearMode::SlidingPositive
    );
    assert!(r.report.phases[0].duration_s > 0.0);
    let w = f.response();
    close(
        r.report.phases[0].duration_s,
        1e-20 / (w[0][0] * 0.1 * f.patch.normal_load_n()),
        1e-34,
    );
    assert!(r.report.surface_heat_j > 0.0);
    close(r.report.surface_heat_j, 0.5e-40 / w[0][0], 1e-53);
    assert!(r.history.is_stopped());
}

#[test]
fn coupled_stop_then_static_rolling_matches_independent_closed_form_work() {
    let u0 = 0.004;
    let r0 = 0.001;
    let f = fixture(0.01, [0.6, 0.1], u0, r0 / R);
    let dt = 0.4;
    let [[a, b], [_b, c]] = f.response();
    let gamma = f.tau * f.patch.normal_load_n() / R;
    let force = -0.1 * f.patch.normal_load_n();
    let lambda = gamma * c;
    let equilibrium = b * force / lambda;
    let difference = r0 - equilibrium;
    let integral = |t: f64| equilibrium * t + difference * (-(-lambda * t).exp_m1()) / lambda;
    let shear = |t: f64| u0 + a * force * t - b * gamma * integral(t);
    let mut lo = 0.0;
    let mut hi = 0.02;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if shear(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let stop = 0.5 * (lo + hi);
    let at_stop = equilibrium + difference * (-lambda * stop).exp();
    let remaining = dt - stop;
    let lambda_static = gamma * (c - b * b / a);
    let static_integral = at_stop * (-(-lambda_static * remaining).exp_m1()) / lambda_static;
    let expected_rolling_integral = integral(stop) + static_integral;
    let first_square = equilibrium * equilibrium * stop
        + 2.0 * equilibrium * difference * (-(-lambda * stop).exp_m1()) / lambda
        + difference * difference * (-(-2.0 * lambda * stop).exp_m1()) / (2.0 * lambda);
    let second_square =
        at_stop * at_stop * (-(-2.0 * lambda_static * remaining).exp_m1()) / (2.0 * lambda_static);
    let expected_shear_impulse = force * stop + (b * gamma / a) * static_integral;
    let r = f.run(dt);
    audit(&f, &r);
    assert_eq!(r.report.phases.len(), 2);
    close(r.report.phases[0].duration_s, stop, 3e-13);
    assert_eq!(
        r.report.phases[1].mode,
        RigidHertzRollingCollinearMode::Sticking
    );
    assert!(r.report.phases[1].bulk_heat_j > 0.0);
    assert_eq!(r.report.phases[1].surface_heat_j, 0.0);
    close(
        r.report.force_impulse_world_ns.x,
        expected_shear_impulse,
        3e-13,
    );
    close(
        r.report.rolling_angular_impulse_world_nms.y,
        -R * gamma * expected_rolling_integral,
        3e-14,
    );
    close(
        r.report.bulk_heat_j,
        gamma * (first_square + second_square),
        5e-14,
    );
    close(
        r.report.phases[1].rolling_rate_after_rad_s,
        at_stop * (-lambda_static * remaining).exp() / R,
        2e-12,
    );
}

#[test]
fn insufficient_static_capacity_causes_physical_reversal_at_the_stop() {
    let f = fixture(0.01, [0.0004, 0.0002], 0.0003, -0.09);
    let r = f.run(0.5);
    audit(&f, &r);
    assert!(r.report.phases.len() >= 2);
    assert!(r.report.phases[0].reached_stop);
    assert_eq!(
        r.report.phases[0].mode,
        RigidHertzRollingCollinearMode::SlidingPositive
    );
    assert_eq!(
        r.report.phases[1].mode,
        RigidHertzRollingCollinearMode::SlidingNegative
    );
    let [[a, b], _] = f.response();
    let static_required = (b / a)
        * (f.tau * f.patch.normal_load_n() / R)
        * (R * r.report.phases[1].rolling_rate_before_rad_s);
    assert!(static_required.abs() > 0.0004 * f.patch.normal_load_n());
    assert!(r.report.phases[1].shear_displacement_m < 0.0);
    assert!(r.report.phases[1].surface_heat_j > 0.0);

    // Continue through the outgoing negative phase's first POSITIVE-time zero.
    // Its initial zero is departure, not arrest. This exercises the analytic
    // extremum used to exclude that zero, then the eventual static branch.
    let [[a, b], [_b, c]] = f.response();
    let gamma = f.tau * f.patch.normal_load_n() / R;
    let lambda = gamma * c;
    let force = 0.0002 * f.patch.normal_load_n();
    let path = |u: f64, r0: f64, force: f64, t: f64| {
        let equilibrium = b * force / lambda;
        let difference = r0 - equilibrium;
        let integral = equilibrium * t + difference * (-(-lambda * t).exp_m1()) / lambda;
        (
            u + a * force * t - b * gamma * integral,
            equilibrium + difference * (-lambda * t).exp(),
        )
    };
    let mut lo = 0.0;
    let mut hi = 0.5;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if path(0.0003, -0.009, -force, mid).0 > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let first_stop = 0.5 * (lo + hi);
    let rolling_at_stop = path(0.0003, -0.009, -force, first_stop).1;
    lo = 0.1;
    hi = 5.0;
    assert!(path(0.0, rolling_at_stop, force, lo).0 < 0.0);
    assert!(path(0.0, rolling_at_stop, force, hi).0 > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if path(0.0, rolling_at_stop, force, mid).0 < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let outgoing_stop = 0.5 * (lo + hi);
    let complete = f.run(5.0);
    audit(&f, &complete);
    assert_eq!(complete.report.phases.len(), 3);
    close(complete.report.phases[0].duration_s, first_stop, 2e-11);
    close(complete.report.phases[1].duration_s, outgoing_stop, 2e-10);
    assert!(complete.report.phases[1].reached_stop);
    assert_eq!(
        complete.report.phases[1].mode,
        RigidHertzRollingCollinearMode::SlidingNegative
    );
    assert_eq!(
        complete.report.phases[2].mode,
        RigidHertzRollingCollinearMode::Sticking
    );
    assert!(complete.history.is_stopped());
    assert!(complete.report.phases[2].bulk_heat_j > 0.0);
    assert_eq!(complete.report.phases[2].surface_heat_j, 0.0);
}

#[test]
fn split_windows_retain_physical_stop_checkpoint_and_original_cumulative_heat() {
    let mut f = fixture(0.01, [0.6, 0.1], 0.004, 0.01);
    let whole = f.run(0.4);
    let first = f.run(0.1);
    assert!(first.history.is_stopped());
    let before = f.bodies;
    let history = f.history.clone();
    f.bodies = first.bodies_after;
    f.history = first.history.clone();
    let second = f.run(0.3);
    assert!(second.history.is_stopped());
    for i in 0..2 {
        close_vec(
            second.bodies_after[i].linear_momentum_world(),
            whole.bodies_after[i].linear_momentum_world(),
            3e-12,
        );
        close_vec(
            second.bodies_after[i].angular_momentum_body(),
            whole.bodies_after[i].angular_momentum_body(),
            3e-13,
        );
    }
    close(
        first.report.surface_heat_j + second.report.surface_heat_j,
        whole.report.surface_heat_j,
        3e-13,
    );
    close(
        first.report.bulk_heat_j + second.report.bulk_heat_j,
        whole.report.bulk_heat_j,
        3e-13,
    );
    close(
        second.history.rolling_history().state().cumulative_heat_j(),
        whole.history.rolling_history().state().cumulative_heat_j(),
        3e-13,
    );
    // A located stop binds the physical inertia as well as its momenta. These
    // tiny changes would leave the new slip below the loose velocity gate.
    let original_properties = f.properties;
    let old = f.properties[0];
    f.properties[0] = MassProperties::new(
        old.mass() * (1.0 + 1e-13),
        Vec3::ZERO,
        old.principal_inertia_body(),
    )
    .unwrap();
    assert!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.3,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    f.properties = original_properties;
    let old = f.properties[1];
    let mut inertia = old.principal_inertia_body();
    inertia.y *= 1.0 + 1e-12;
    f.properties[1] = MassProperties::new(old.mass(), Vec3::ZERO, inertia).unwrap();
    assert!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.3,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    f.properties = original_properties;
    f.bodies[0] = RigidBodyState::new(
        f.bodies[0].pose(),
        f.bodies[0]
            .linear_momentum_world()
            .add(Vec3::new(1e-15, 0.0, 0.0)),
        f.bodies[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.3,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    f.bodies = before;
    f.history = history;
    assert_eq!(whole, f.run(0.4));
}

#[test]
fn full_inertia_subspace_domain_budgets_and_cancellation_refuse_atomically() {
    let mut f = fixture(0.01, [0.6, 0.1], 0.004, 0.01);
    let before = f.bodies;
    let history = f.history.clone();
    let prepared = f.prepare(0.4);
    let good = prepared.advance(&f.hertz, &CancelGate::new()).unwrap();
    assert!(
        prepared
            .advance_at_temperature(&f.hertz, 500.0, &CancelGate::new())
            .is_err()
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(prepared.advance(&f.hertz, &cancelled).is_err());
    assert_eq!(
        good,
        prepared.advance(&f.hertz, &CancelGate::new()).unwrap()
    );
    let mut c = controls();
    c.maximum_phases = 1;
    assert!(matches!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.4,
            c,
            &CancelGate::new()
        )
        .unwrap()
        .advance(&f.hertz, &CancelGate::new()),
        Err(RigidHertzRollingCollinearError::PhaseBudget)
    ));
    c = controls();
    c.maximum_root_iterations = 1;
    assert!(matches!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.4,
            c,
            &CancelGate::new()
        )
        .unwrap()
        .advance(&f.hertz, &CancelGate::new()),
        Err(RigidHertzRollingCollinearError::RootBudget)
    ));
    assert_eq!(f.bodies, before);
    assert_eq!(f.history, history);
    let old = f.properties[1];
    f.properties[1] = MassProperties::from_principal_frame(
        old.mass(),
        Vec3::ZERO,
        old.principal_inertia_body(),
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.4).unwrap(),
    )
    .unwrap();
    assert!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.4,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    f.properties[1] = old;
    f.bodies[0] = RigidBodyState::new(
        f.bodies[0].pose(),
        f.bodies[0]
            .linear_momentum_world()
            .add(Vec3::new(0.0, 0.001, 0.0)),
        f.bodies[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        PreparedRigidHertzRollingCollinear::new(
            &f.bodies,
            &f.properties,
            &f.history,
            f.contact(),
            0.4,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );

    // A tiny circular patch must not hide a real excluded spin merely because
    // its work coordinate is a*omega. The original full pair admits this
    // source/geometry; only the declared two-coordinate restriction refuses.
    let mut tiny = fixture(0.0, [0.8, 0.1], 0.004, 0.0);
    let normal = HertzSpherePlane {
        effective_radius: R,
        reduced_modulus: 1e6 / (1.0 - 0.25_f64.powi(2)),
    }
    .response(&tiny.interface, 1e-39)
    .unwrap();
    tiny.patch = RollingPatchReceipt::new(
        "collinear original patch",
        "Hertz normal",
        "synthetic normal source",
        InputAuthority::SyntheticFixture,
        normal.normal_force_n,
        PI * normal.contact_radius_m.powi(2),
        PatchCurvature::Principal {
            first_per_m: 1.0 / R,
            second_per_m: 1.0 / R,
        },
    )
    .unwrap();
    tiny.hertz = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: normal.normal_force_n,
        contact_radius_m: normal.contact_radius_m,
        ..tiny.hertz.parameters()
    })
    .unwrap();
    tiny.bodies[0] = RigidBodyState::new(
        tiny.bodies[0].pose(),
        tiny.bodies[0].linear_momentum_world(),
        tiny.properties[0]
            .angular_momentum_body_checked(Vec3::new(0.0, 0.0, 1.0))
            .unwrap(),
    )
    .unwrap();
    PreparedRigidHertzRollingPair::new(
        &tiny.bodies,
        &tiny.properties,
        tiny.history.rolling_history(),
        tiny.contact(),
        0.01,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        PreparedRigidHertzRollingCollinear::new(
            &tiny.bodies,
            &tiny.properties,
            &tiny.history,
            tiny.contact(),
            0.01,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidHertzRollingCollinearError::Invalid(
            "initial motion leaves the invariant longitudinal/lateral-rolling subspace"
        ))
    ));
}

#[test]
fn zero_friction_and_zero_load_keep_real_rolling_motion_and_original_limits() {
    let mut f = fixture(0.01, [0.0, 0.0], 0.004, 0.02);
    let r = f.run(0.4);
    audit(&f, &r);
    assert_eq!(r.report.surface_heat_j, 0.0);
    assert!(r.report.bulk_heat_j > 0.0);
    assert!(!r.history.is_stopped());
    f.patch = RollingPatchReceipt::new(
        "collinear original patch",
        "Hertz normal",
        "synthetic normal source",
        InputAuthority::SyntheticFixture,
        0.0,
        0.0,
        PatchCurvature::Principal {
            first_per_m: 1.0 / R,
            second_per_m: 1.0 / R,
        },
    )
    .unwrap();
    f.hertz = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 0.0,
        ..f.hertz.parameters()
    })
    .unwrap();
    let released = f.run(0.4);
    assert_eq!(released.bodies_after, f.bodies);
    assert_eq!(released.report.surface_heat_j, 0.0);
    assert_eq!(released.report.bulk_heat_j, 0.0);
}
