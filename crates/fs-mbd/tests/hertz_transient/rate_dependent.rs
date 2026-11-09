//! G1/G3/G4: actual rate-law pressure traction changes the original native pair.
use super::*;
use fs_mbd::tangential_friction::hertz::rate_dependent::*;
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};

fn quadrature() -> HertzFrictionQuadrature {
    HertzFrictionQuadrature {
        radial_panels: 16,
        angular_panels: 16,
    }
}
fn patch(law: FrictionLaw) -> HertzFrictionPatch {
    HertzFrictionPatch::new(10.0, R, law).unwrap()
}
fn affine() -> HertzFrictionPatch {
    patch(FrictionLaw::VelocityDependent {
        static_mu: 0.9,
        mu_zero: 0.4,
        slope_per_speed: 0.8,
    })
}
fn contact<'a>(f: &Fixture, law: &'a HertzFrictionPatch) -> RigidHertzRateContact<'a> {
    RigidHertzRateContact {
        binding: f.binding,
        frame: f.frame,
        law,
    }
}
fn centre_pair(velocity: f64, spin: f64) -> Fixture {
    let mut f = fixture();
    f.binding.arm_a_body = Vec3::ZERO;
    f.binding.arm_b_body = Vec3::ZERO;
    let pose = Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap();
    f.bodies = [
        RigidBodyState::new(
            pose,
            Vec3::new(f.properties[0].mass() * velocity, 0.0, 0.0),
            Vec3::new(0.0, 0.0, f.properties[0].principal_inertia_body().z * spin),
        )
        .unwrap(),
        RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO).unwrap(),
    ];
    f
}
fn run(f: &Fixture, law: &HertzFrictionPatch, dt: f64) -> RigidHertzRateTransientCandidate {
    prepare_rigid_hertz_rate_transient(
        &f.bodies,
        &f.properties,
        contact(f, law),
        dt,
        controls(),
        quadrature(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap()
}

#[test]
fn affine_sliding_and_spin_recoil_follow_independent_exact_time_and_travel() {
    // Pressure moments give F=-N(mu0+b*u) and
    // M=-N(3*pi*mu0*a/16 + 2*b*a^2*omega/5). Each pure native mode
    // therefore obeys y'=-A-B*y, including its exact finite stopping event.
    let law = affine();
    for spin in [false, true] {
        let initial = if spin { 3.0 } else { 0.3 };
        let f = centre_pair(
            if spin { 0.0 } else { initial },
            if spin { initial } else { 0.0 },
        );
        let mobility: f64 = f
            .properties
            .iter()
            .map(|p| {
                if spin {
                    p.principal_inertia_body().z.recip()
                } else {
                    p.mass().recip()
                }
            })
            .sum();
        let a = mobility
            * 10.0
            * 0.4
            * if spin {
                3.0 * core::f64::consts::PI * R / 16.0
            } else {
                1.0
            };
        let b = mobility * 10.0 * 0.8 * if spin { 2.0 * R * R / 5.0 } else { 1.0 };
        let stop = (b * initial / a).ln_1p() / b;
        let mut stopped: Option<[RigidBodyState; 2]> = None;
        for dt in [0.4 * stop, 1.7 * stop, 4.0 * stop] {
            let moving = dt.min(stop);
            let decay = (-b * moving).exp();
            let endpoint = if dt >= stop {
                0.0
            } else {
                initial * decay + (a / b) * (-b * moving).exp_m1()
            };
            let travel = (initial + a / b) * -(-b * moving).exp_m1() / b - (a / b) * moving;
            let expected_heat = 0.5 * (initial - endpoint) * (initial + endpoint) / mobility;
            let result = run(&f, &law, dt);
            let mechanical = &result.report.mechanical;
            let path = &mechanical.transient;
            close(
                if spin {
                    mechanical.spin_rad_s[1]
                } else {
                    mechanical.slip_velocity_m_s[1][0]
                },
                endpoint,
                4e-9,
            );
            close(
                if spin {
                    path.spin_increment_rad
                } else {
                    path.slip_increment_m[0]
                },
                travel,
                if spin { 5e-9 } else { 2e-10 },
            );
            close(path.dissipated_energy_j, expected_heat, 8e-11);
            close(
                kinetic(result.bodies_after, f.properties) - kinetic(f.bodies, f.properties),
                -expected_heat,
                8e-11,
            );
            close(
                mechanical.audit.total_energy_residual_j,
                0.0,
                controls().energy_tolerance_j,
            );
            close_vec(
                mechanical.audit.linear_momentum_change_ns,
                Vec3::ZERO,
                controls().linear_momentum_tolerance_ns,
            );
            close_vec(
                mechanical.audit.angular_momentum_change_nms,
                Vec3::ZERO,
                controls().angular_momentum_tolerance_nms,
            );
            assert!(result.report.kinetic_coefficient_bounds[0] > 0.0);
            if dt > stop {
                let bracket = path.stopping_time_bracket_s.unwrap();
                assert!(stop >= bracket[0] - transient().stopping_time_tolerance_s);
                assert!(stop <= bracket[1] + transient().stopping_time_tolerance_s);
                if let Some(previous) = stopped {
                    for i in 0..2 {
                        close_vec(
                            result.bodies_after[i].linear_momentum_world(),
                            previous[i].linear_momentum_world(),
                            2e-11,
                        );
                        close_vec(
                            result.bodies_after[i].angular_momentum_body(),
                            previous[i].angular_momentum_body(),
                            2e-12,
                        );
                    }
                }
                stopped = Some(result.bodies_after);
            } else {
                assert!(path.stopping_time_bracket_s.is_none());
            }
        }
    }
}

#[test]
fn coulomb_reduction_and_repeated_curve_trials_retain_the_original_native_pair() {
    let f = fixture();
    let law = patch(FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.4,
    });
    let dt = 0.008;
    let original = prepare_rigid_hertz_transient(
        &f.bodies,
        &f.properties,
        f.contact(&super::law()),
        dt,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    let reduced = run(&f, &law, dt);
    assert_eq!(reduced.bodies_after, original.bodies_after);
    assert_eq!(reduced.report.mechanical, original.report);
    let prepared = PreparedRigidHertzRateTransient::new(
        &f.bodies,
        &f.properties,
        contact(&f, &law),
        dt,
        controls(),
        quadrature(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    for curve in [
        affine(),
        patch(FrictionLaw::Hyperbolic {
            static_mu: 0.8,
            kinetic_mu: 0.2,
            characteristic_speed: 0.1,
        }),
    ] {
        let retry = prepared.advance(&curve, &CancelGate::new()).unwrap();
        assert_eq!(retry, run(&f, &curve, dt));
        assert_ne!(retry.bodies_after, reduced.bodies_after);
    }
    for (normal, radius) in [(11.0, R), (10.0, 1.1 * R)] {
        let foreign = HertzFrictionPatch::new(normal, radius, law.friction_law().clone()).unwrap();
        assert!(prepared.advance(&foreign, &CancelGate::new()).is_err());
    }
    assert_eq!(prepared.advance(&law, &CancelGate::new()).unwrap(), reduced);
}

#[test]
fn velocity_weakening_can_feed_spin_while_total_native_energy_decreases() {
    let f = centre_pair(0.3, 0.1);
    let law = HertzFrictionPatch::new(
        2.0,
        R,
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.0,
            characteristic_speed: 0.3,
            viscous_per_speed: 0.0,
        },
    )
    .unwrap();
    // Independently differentiated small-spin pressure moment is positive
    // when mu(U)/U + mu'(U)<0. A componentwise opposing-work gate is wrong.
    let (mu, slope) = law
        .friction_law()
        .kinetic_coefficient_and_slope(0.3)
        .unwrap();
    let leading_moment = -2.0 * R * R * 0.1 * (mu / 0.3 + slope) / 5.0;
    assert!(leading_moment > 0.0);
    let dt = 0.0001;
    let result = run(&f, &law, dt);
    let m = &result.report.mechanical;
    assert!(m.spin_rad_s[1] > m.spin_rad_s[0]);
    assert!(m.spin_angular_impulse_world_nms.z > 0.0);
    assert!(m.slip_velocity_m_s[1][0] < m.slip_velocity_m_s[0][0]);
    assert!(m.dissipated_energy_j > 0.0);
    close(
        m.spin_angular_impulse_world_nms.z / dt,
        leading_moment,
        0.01 * leading_moment,
    );
    close(
        kinetic(result.bodies_after, f.properties) - kinetic(f.bodies, f.properties),
        -m.dissipated_energy_j,
        2e-10,
    );
}

#[test]
fn rate_law_refusals_and_publication_preserve_unrelated_bodies_and_exact_retry() {
    let f = fixture();
    let law = affine();
    let mut bodies = vec![f.bodies[0], f.bodies[1], f.bodies[0]];
    let properties = vec![f.properties[0], f.properties[1], f.properties[0]];
    let saved = bodies.clone();
    let dt = 0.008;
    let exhausted = HertzTransientControls {
        maximum_evaluations: 1,
        ..transient()
    };
    assert!(
        apply_rigid_hertz_rate_transient(
            &mut bodies,
            &properties,
            contact(&f, &law),
            dt,
            controls(),
            quadrature(),
            exhausted,
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_hertz_rate_transient(
            &mut bodies,
            &properties,
            contact(&f, &law),
            dt,
            controls(),
            quadrature(),
            transient(),
            &gate,
        )
        .is_err()
    );
    assert_eq!(bodies, saved);
    let expected = run(&f, &law, dt);
    let report = apply_rigid_hertz_rate_transient(
        &mut bodies,
        &properties,
        contact(&f, &law),
        dt,
        controls(),
        quadrature(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, expected.report);
    assert_eq!(bodies[..2], expected.bodies_after);
    assert_eq!(bodies[2], saved[2]);
    // The same actual normal-motion owner rejects an oblique arm which would
    // require an unsolved normal reaction to keep this pressure patch closed.
    let mut oblique = fixture();
    oblique.binding.arm_a_body.z = R;
    oblique.bodies[0] = RigidBodyState::new(
        Pose::new(
            oblique.binding.arm_a_body.scale(-1.0),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        f.bodies[0].linear_momentum_world(),
        f.bodies[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        prepare_rigid_hertz_rate_transient(
            &oblique.bodies,
            &oblique.properties,
            contact(&oblique, &law),
            dt,
            controls(),
            quadrature(),
            transient(),
            &CancelGate::new(),
        )
        .is_err()
    );
}
