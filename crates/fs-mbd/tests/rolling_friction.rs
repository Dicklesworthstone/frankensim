//! Original finite-body rolling against independent constrained mechanics.
//! All material/geometry values are synthetic; no experimental calibration.

use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::rolling::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::{
    ApplicabilityRange, FrictionLaw, HertzSpherePlane, InputAuthority, InterfaceMedium,
    InterfaceSystemRef,
};

const R: f64 = 0.1;
const D: f64 = 0.05;
const TAU: f64 = 0.01;
const THETA: f64 = 0.47;
const INITIAL_RATE: [f64; 2] = [0.03, -0.02];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tol={tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn z_cross(v: Vec3) -> Vec3 {
    Vec3::new(-v.y, v.x, 0.0)
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
    law: RollingLossLaw,
    patch: RollingPatchReceipt,
    interface: InterfaceSystemRef,
    friction: FrictionLaw,
    work: RollingWorkOwnership,
}
impl Fixture {
    fn contact(&self) -> RigidRollingContact<'_> {
        RigidRollingContact {
            binding: self.binding,
            frame: self.frame,
            patch: &self.patch,
            interface: &self.interface,
            law: &self.law,
            static_friction: &self.friction,
            temperature_k: 300.0,
            excitation_frequency_hz: 10.0,
            work: &self.work,
        }
    }
}
fn material(tau: f64, surface: DeformingSphereSurface) -> ViscoelasticSphereParameters {
    let young = 1e6;
    let nu = 0.25;
    ViscoelasticSphereParameters {
        sphere_surface: surface,
        radius_m: R,
        young_modulus_pa: young,
        poisson_ratio: nu,
        density_kg_m3: 1000.0,
        shear_viscosity_pa_s: tau * young / (2.0 * (1.0 + nu)),
        bulk_viscosity_pa_s: tau * young / (3.0 * (1.0 - 2.0 * nu)),
    }
}
fn rolling_law(tau: f64, surface: DeformingSphereSurface, source: &str) -> RollingLossLaw {
    RollingLossLaw::ViscoelasticSphere(
        ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            source,
            InputAuthority::SyntheticFixture,
            material(tau, surface),
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
fn fixture(tau: f64, rate: [f64; 2], indentation: f64) -> Fixture {
    let interface = InterfaceSystemRef::new(
        "fixture/sphere-A/finite-plane-B",
        "fixture/rolling-history",
        "fixture/original-interface",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let p = material(tau, DeformingSphereSurface::SurfaceA);
    let normal = HertzSpherePlane {
        effective_radius: R,
        reduced_modulus: p.young_modulus_pa / (1.0 - p.poisson_ratio * p.poisson_ratio),
    }
    .response(&interface, indentation)
    .unwrap();
    let patch = RollingPatchReceipt::new(
        "fixture/rolling-patch",
        "original/Hertz",
        "fixture/normal-source",
        InputAuthority::SyntheticFixture,
        normal.normal_force_n,
        core::f64::consts::PI * normal.contact_radius_m * normal.contact_radius_m,
        PatchCurvature::Principal {
            first_per_m: R.recip(),
            second_per_m: R.recip(),
        },
    )
    .unwrap();
    let mass = (4.0 * core::f64::consts::PI / 3.0) * p.density_kg_m3 * R.powi(3);
    let inertia = 0.4 * mass * R * R;
    let properties = [
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap(),
        MassProperties::new(7.0, Vec3::ZERO, Vec3::new(0.08, 0.13, 0.18)).unwrap(),
    ];
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), THETA).unwrap();
    let omega = rotation.rotate_body_to_world(Vec3::new(rate[0], rate[1], 0.0));
    let bodies = [
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, R), UnitQuaternion::IDENTITY).unwrap(),
            z_cross(omega).scale(-mass * R),
            omega.scale(inertia),
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, -D), rotation).unwrap(),
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
            arm_b_body: Vec3::new(0.0, 0.0, D),
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        law: rolling_law(tau, DeformingSphereSurface::SurfaceA, "fixture/viscosity"),
        patch,
        interface,
        friction: FrictionLaw::Coulomb {
            static_mu: 0.3,
            kinetic_mu: 0.2,
        },
        work: RollingWorkOwnership::new(
            "fixture/rolling-patch",
            "interval-1",
            "relative-tangent-rolling-vector",
            RollingLossChannel::RollingHysteresis,
        )
        .unwrap(),
    }
}
fn controls(step: f64) -> RigidRollingControls {
    RigidRollingControls {
        maximum_bodies: 2,
        maximum_substeps: 4096,
        maximum_step_s: step,
        geometry_tolerance_m: 1e-12,
        velocity_tolerance_m_s: 1e-11,
        angular_velocity_tolerance_rad_s: 1e-10,
        linear_momentum_tolerance_ns: 1e-11,
        angular_momentum_tolerance_nms: 1e-11,
        energy_tolerance_j: 1e-12,
    }
}
fn advance(
    f: &Fixture,
    bodies: &mut [RigidBodyState; 2],
    history: &mut RigidRollingHistory,
    dt: f64,
    c: RigidRollingControls,
) -> RigidRollingReport {
    apply_rigid_rolling(
        bodies,
        &f.properties,
        history,
        f.contact(),
        dt,
        c,
        &CancelGate::new(),
    )
    .unwrap()
}

// Independent force/torque elimination in each physical principal direction.
// For sphere inertia I_s and plane inertia I_p, relative no-slip acceleration
// is a F + b J M, a=1/m_s+1/m_p+R²/I_s+D²/I_p, b=R/I_s-D/I_p.
// Thus F=-(b/a) J M and rolling inverse inertia is 1/I_s+1/I_p-b²/a.
// No production response matrix, projection factor, or returned impulse enters.
fn reference(
    f: &Fixture,
    duration: f64,
    steps: Option<usize>,
    initial: [f64; 2],
) -> ([RigidBodyState; 2], f64, [f64; 2]) {
    let sphere_inertia = f.properties[0].principal_inertia_body().x;
    let plane_inertias = [
        f.properties[1].principal_inertia_body().x,
        f.properties[1].principal_inertia_body().y,
    ];
    let rotation = f.bodies[1].pose().orientation();
    let directions = [
        rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        rotation.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0)),
    ];
    // The proportional Kelvin-Voigt fixture has A=tau independently of nu.
    let damping = TAU * f.patch.normal_load_n() * R;
    let mut couple = Vec3::ZERO;
    let mut force = Vec3::ZERO;
    let mut heat = 0.0;
    let mut kappa = [0.0; 2];
    for i in 0..2 {
        let a = 1.0 / f.properties[0].mass()
            + 1.0 / f.properties[1].mass()
            + R * R / sphere_inertia
            + D * D / plane_inertias[i];
        let b = R / sphere_inertia - D / plane_inertias[i];
        kappa[i] = 1.0 / sphere_inertia + 1.0 / plane_inertias[i] - b * b / a;
        let rate = damping * kappa[i];
        let ratio = if let Some(n) = steps {
            let alpha = 0.5 * duration * rate / n as f64;
            ((1.0 - alpha) / (1.0 + alpha)).powi(n as i32)
        } else {
            (-duration * rate).exp()
        };
        let final_rate = initial[i] * ratio;
        let impulse = (final_rate - initial[i]) / kappa[i];
        couple = couple.add(directions[i].scale(impulse));
        force = force.add(z_cross(directions[i]).scale(-b * impulse / a));
        heat += 0.5 * (initial[i] * initial[i] - final_rate * final_rate) / kappa[i];
    }
    let sphere_angular = f.bodies[0]
        .angular_momentum_body()
        .add(z_cross(force).scale(-R))
        .add(couple);
    let plane_angular_world = z_cross(force).scale(-D).sub(couple);
    let bodies = [
        RigidBodyState::new(
            f.bodies[0].pose(),
            f.bodies[0].linear_momentum_world().add(force),
            sphere_angular,
        )
        .unwrap(),
        RigidBodyState::new(
            f.bodies[1].pose(),
            f.bodies[1].linear_momentum_world().sub(force),
            rotation.rotate_world_to_body(plane_angular_world),
        )
        .unwrap(),
    ];
    (bodies, heat, kappa)
}

#[test]
fn both_axes_slow_the_sphere_and_recoil_the_finite_support_with_original_heat() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let report = advance(&f, &mut bodies, &mut history, 0.2, controls(0.2));
    let (expected, heat, eigenvalues) = reference(&f, 0.2, Some(1), INITIAL_RATE);
    assert_eq!(report.substeps, 1);
    assert!(report.rolling_mobility[0][1].abs() > 1.0);
    for i in 0..2 {
        close_vec(
            bodies[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            2e-16,
        );
        close_vec(
            bodies[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            2e-17,
        );
        assert_eq!(bodies[i].pose(), f.bodies[i].pose());
    }
    let q = f.bodies[1].pose().orientation();
    for (i, direction) in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)]
        .into_iter()
        .enumerate()
    {
        let u = q.rotate_body_to_world(direction);
        close(
            report.rolling_mobility[0][0] * u.x + report.rolling_mobility[0][1] * u.y,
            eigenvalues[i] * u.x,
            5e-13,
        );
        close(
            report.rolling_mobility[1][0] * u.x + report.rolling_mobility[1][1] * u.y,
            eigenvalues[i] * u.y,
            5e-13,
        );
    }
    assert!(norm(bodies[1].linear_momentum_world()) > 1e-4);
    assert!(norm(bodies[1].angular_momentum_body()) > 1e-6);
    assert!(norm(bodies[0].linear_momentum_world()) < norm(f.bodies[0].linear_momentum_world()));
    close(report.dissipated_heat_j, heat, 2e-18);
    close(history.state().cumulative_heat_j(), heat, 2e-18);
    close(report.static_reaction_work_j, 0.0, 1e-19);
    close(
        report.kinetic_before_j - report.kinetic_after_j,
        heat,
        2e-18,
    );
    close(report.total_energy_residual_j, 0.0, 2e-18);
    close(report.work_energy_residual_j, 0.0, 2e-18);
    close_vec(report.linear_momentum_change_ns, Vec3::ZERO, 1e-16);
    close_vec(report.angular_momentum_change_nms, Vec3::ZERO, 1e-17);
    assert_eq!(report.deforming_body, 0);
    assert!(report.maximum_static_force_n < report.static_capacity_n);
    assert_eq!(report.final_step.generalized_work.ownership, f.work);
}

#[test]
fn unrelated_large_energy_cannot_hide_pair_work_or_change_original_body_indices() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let mut isolated = f.bodies;
    let mut isolated_history = RigidRollingHistory::new();
    let mut expected = advance(&f, &mut isolated, &mut isolated_history, 0.2, controls(0.2));
    let spectator = RigidBodyState::new(
        Pose::new(Vec3::new(11.0, -7.0, 5.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(1e18, -2e18, 3e18),
        Vec3::new(-3e17, 2e17, 1e17),
    )
    .unwrap();
    let spectator_properties =
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(3.0, 4.0, 5.0)).unwrap();
    let spectator_energy = spectator.kinetic_energy(spectator_properties).unwrap();
    assert!(spectator_energy.is_finite());
    assert!(expected.dissipated_heat_j > controls(0.2).energy_tolerance_j);
    // A whole-scene subtraction loses this real, independently verified heat.
    assert_eq!(
        spectator_energy + expected.kinetic_before_j,
        spectator_energy + expected.kinetic_after_j,
    );
    let mut scene = [spectator, f.bodies[1], f.bodies[0]];
    let accepted = scene;
    let properties = [spectator_properties, f.properties[1], f.properties[0]];
    let contact = RigidRollingContact {
        binding: RigidContactBinding {
            body_a: 2,
            body_b: 1,
            ..f.binding
        },
        ..f.contact()
    };
    let mut history = RigidRollingHistory::new();
    // The existing scene-size admission still applies before pair staging.
    assert!(
        apply_rigid_rolling(
            &mut scene,
            &properties,
            &mut history,
            contact,
            0.2,
            controls(0.2),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(scene, accepted);
    assert_eq!(history, RigidRollingHistory::new());
    let report = apply_rigid_rolling(
        &mut scene,
        &properties,
        &mut history,
        contact,
        0.2,
        RigidRollingControls {
            maximum_bodies: 3,
            ..controls(0.2)
        },
        &CancelGate::new(),
    )
    .unwrap();
    expected.deforming_body = 2;
    assert_eq!(report, expected);
    assert_eq!(history, isolated_history);
    assert_eq!(scene[2], isolated[0]);
    assert_eq!(scene[1], isolated[1]);
    assert_eq!(scene[0], spectator);
}

#[test]
fn anisotropic_two_axis_decay_converges_to_independent_matrix_exponential() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let (exact, exact_heat, _) = reference(&f, 2.0, None, INITIAL_RATE);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut bodies = f.bodies;
        let mut history = RigidRollingHistory::new();
        let report = advance(&f, &mut bodies, &mut history, 2.0, controls(2.0 / n as f64));
        assert_eq!(report.substeps, n);
        let (discrete, discrete_heat, _) = reference(&f, 2.0, Some(n), INITIAL_RATE);
        for i in 0..2 {
            close_vec(
                bodies[i].linear_momentum_world(),
                discrete[i].linear_momentum_world(),
                2e-16,
            );
            close_vec(
                bodies[i].angular_momentum_body(),
                discrete[i].angular_momentum_body(),
                2e-17,
            );
        }
        close(report.dissipated_heat_j, discrete_heat, 2e-18);
        let momentum_error = norm(
            bodies[0]
                .linear_momentum_world()
                .sub(exact[0].linear_momentum_world()),
        );
        let heat_error = (report.dissipated_heat_j - exact_heat).abs();
        assert!(momentum_error > 1e-10 && heat_error > 1e-12);
        errors.push((momentum_error, heat_error));
    }
    for pair in errors.windows(2) {
        assert!((3.95..4.05).contains(&(pair[0].0 / pair[1].0)));
        assert!((3.95..4.05).contains(&(pair[0].1 / pair[1].1)));
    }
}

#[test]
fn automatic_substeps_prevent_reversal_and_exact_zero_limits_do_not_snap() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let report = advance(&f, &mut bodies, &mut history, 20.0, controls(20.0));
    assert!(report.substeps > 2);
    let q = f.bodies[1].pose().orientation();
    let final_rate = q.rotate_world_to_body(Vec3::new(
        report.rolling_rate_after_rad_s[0],
        report.rolling_rate_after_rad_s[1],
        0.0,
    ));
    assert!(final_rate.x > 0.0 && final_rate.x < INITIAL_RATE[0]);
    assert!(final_rate.y < 0.0 && final_rate.y > INITIAL_RATE[1]);
    let (expected, heat, _) = reference(&f, 20.0, Some(report.substeps), INITIAL_RATE);
    close_vec(
        bodies[0].linear_momentum_world(),
        expected[0].linear_momentum_world(),
        2e-16,
    );
    close(report.dissipated_heat_j, heat, 2e-18);
    for fixture in [
        fixture(TAU, [0.0; 2], 0.001),
        fixture(0.0, INITIAL_RATE, 0.001),
        fixture(TAU, INITIAL_RATE, 0.0),
    ] {
        let mut body = fixture.bodies;
        let mut state = RigidRollingHistory::new();
        let mut bounded = controls(1e-6);
        bounded.maximum_substeps = 1;
        let zero = advance(&fixture, &mut body, &mut state, 1e6, bounded);
        assert_eq!(body, fixture.bodies);
        assert_eq!(zero.dissipated_heat_j, 0.0);
        assert_eq!(zero.traction_impulse_world_ns, Vec3::ZERO);
        assert_eq!(zero.rolling_angular_impulse_world_nms, Vec3::ZERO);
        assert_eq!(zero.substeps, 1);
        assert_eq!(state.state().cumulative_heat_j(), 0.0);
    }
}

fn compose(q: UnitQuaternion, r: UnitQuaternion) -> UnitQuaternion {
    let [a, b, c, d] = q.components();
    let [e, f, g, h] = r.components();
    UnitQuaternion::new(
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )
    .unwrap()
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

#[test]
fn full_tensor_body_coordinates_preserve_rolling_recoil_mobility_and_heat() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let mut principal_bodies = f.bodies;
    let mut principal_history = RigidRollingHistory::new();
    let expected = advance(
        &f,
        &mut principal_bodies,
        &mut principal_history,
        0.3,
        controls(0.1),
    );
    // Rotate each BODY coordinate system, retaining the same world geometry,
    // material, velocity, inertia tensor and physical initial momentum.
    // The support now has products of inertia in its geometry frame.
    let labels = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.81).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 0.7, 1.3), -0.63).unwrap(),
    ];
    let mut dense = fixture(TAU, INITIAL_RATE, 0.001);
    for k in 0..2 {
        let [w, x, y, z] = labels[k].components();
        let inverse = UnitQuaternion::new(w, -x, -y, -z).unwrap();
        dense.properties[k] = MassProperties::from_principal_frame(
            f.properties[k].mass(),
            Vec3::ZERO,
            f.properties[k].principal_inertia_body(),
            labels[k],
        )
        .unwrap();
        dense.bodies[k] = RigidBodyState::new(
            Pose::new(
                f.bodies[k].pose().position_world(),
                compose(f.bodies[k].pose().orientation(), inverse),
            )
            .unwrap(),
            f.bodies[k].linear_momentum_world(),
            labels[k].rotate_body_to_world(f.bodies[k].angular_momentum_body()),
        )
        .unwrap();
    }
    assert!(dense.properties[1].inertia_tensor_body().m[1].abs() > 0.01);
    dense.binding.arm_a_body = labels[0].rotate_body_to_world(f.binding.arm_a_body);
    dense.binding.arm_b_body = labels[1].rotate_body_to_world(f.binding.arm_b_body);
    let mut bodies = dense.bodies;
    let mut history = RigidRollingHistory::new();
    let actual = advance(&dense, &mut bodies, &mut history, 0.3, controls(0.1));
    for k in 0..2 {
        close_vec(
            bodies[k].linear_momentum_world(),
            principal_bodies[k].linear_momentum_world(),
            3e-16,
        );
        close_vec(
            bodies[k]
                .pose()
                .orientation()
                .rotate_body_to_world(bodies[k].angular_momentum_body()),
            principal_bodies[k]
                .pose()
                .orientation()
                .rotate_body_to_world(principal_bodies[k].angular_momentum_body()),
            3e-17,
        );
        close(
            actual.rolling_rate_after_rad_s[k],
            expected.rolling_rate_after_rad_s[k],
            3e-16,
        );
        for j in 0..2 {
            close(
                actual.rolling_mobility[k][j],
                expected.rolling_mobility[k][j],
                5e-12,
            );
        }
    }
    for i in 0..4 {
        for j in 0..4 {
            close(actual.mobility[i][j], expected.mobility[i][j], 5e-12);
        }
    }
    close_vec(
        actual.traction_impulse_world_ns,
        expected.traction_impulse_world_ns,
        3e-16,
    );
    close_vec(
        actual.rolling_angular_impulse_world_nms,
        expected.rolling_angular_impulse_world_nms,
        3e-17,
    );
    close(actual.dissipated_heat_j, expected.dissipated_heat_j, 3e-18);
    close(
        history.state().cumulative_heat_j(),
        principal_history.state().cumulative_heat_j(),
        3e-18,
    );
    close(actual.kinetic_after_j, expected.kinetic_after_j, 3e-18);
    close(actual.static_reaction_work_j, 0.0, 1e-19);
    close(actual.work_energy_residual_j, 0.0, 3e-18);
    close(actual.total_energy_residual_j, 0.0, 3e-18);
    close_vec(actual.linear_momentum_change_ns, Vec3::ZERO, 1e-16);
    close_vec(actual.angular_momentum_change_nms, Vec3::ZERO, 1e-17);
}

#[test]
fn world_rotation_surface_order_and_rate_reversal_preserve_physical_work() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let mut original = f.bodies;
    let original_report = advance(
        &f,
        &mut original,
        &mut RigidRollingHistory::new(),
        0.3,
        controls(0.1),
    );
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.73).unwrap();
    let mut rotated = fixture(TAU, INITIAL_RATE, 0.001);
    for body in &mut rotated.bodies {
        *body = RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(body.pose().position_world()),
                compose(q, body.pose().orientation()),
            )
            .unwrap(),
            q.rotate_body_to_world(body.linear_momentum_world()),
            body.angular_momentum_body(),
        )
        .unwrap();
    }
    rotated.frame = TangentFrame::new(
        array(q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
        array(q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
    )
    .unwrap();
    let mut actual = rotated.bodies;
    let report = advance(
        &rotated,
        &mut actual,
        &mut RigidRollingHistory::new(),
        0.3,
        controls(0.1),
    );
    for i in 0..2 {
        close_vec(
            actual[i].linear_momentum_world(),
            q.rotate_body_to_world(original[i].linear_momentum_world()),
            3e-16,
        );
        close_vec(
            actual[i].angular_momentum_body(),
            original[i].angular_momentum_body(),
            3e-17,
        );
    }
    close(
        report.dissipated_heat_j,
        original_report.dissipated_heat_j,
        3e-18,
    );
    let mut swapped = fixture(TAU, INITIAL_RATE, 0.001);
    swapped.bodies.swap(0, 1);
    swapped.properties.swap(0, 1);
    swapped.binding = RigidContactBinding {
        body_a: 0,
        body_b: 1,
        arm_a_body: f.binding.arm_b_body,
        arm_b_body: f.binding.arm_a_body,
    };
    swapped.frame = TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap();
    swapped.law = rolling_law(TAU, DeformingSphereSurface::SurfaceB, "fixture/viscosity");
    swapped.interface = InterfaceSystemRef::new(
        "fixture/finite-plane-A/sphere-B",
        "fixture/rolling-history",
        "fixture/original-interface",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let mut swapped_bodies = swapped.bodies;
    let swapped_report = advance(
        &swapped,
        &mut swapped_bodies,
        &mut RigidRollingHistory::new(),
        0.3,
        controls(0.1),
    );
    for i in 0..2 {
        close_vec(
            swapped_bodies[1 - i].linear_momentum_world(),
            original[i].linear_momentum_world(),
            3e-16,
        );
        close_vec(
            swapped_bodies[1 - i].angular_momentum_body(),
            original[i].angular_momentum_body(),
            3e-17,
        );
    }
    assert_eq!(swapped_report.deforming_body, 1);
    close(
        swapped_report.dissipated_heat_j,
        original_report.dissipated_heat_j,
        3e-18,
    );
    let reversed = fixture(TAU, INITIAL_RATE.map(|v| -v), 0.001);
    let mut reversed_bodies = reversed.bodies;
    let reverse_report = advance(
        &reversed,
        &mut reversed_bodies,
        &mut RigidRollingHistory::new(),
        0.3,
        controls(0.1),
    );
    for i in 0..2 {
        close_vec(
            reversed_bodies[i].linear_momentum_world(),
            original[i].linear_momentum_world().scale(-1.0),
            3e-16,
        );
        close_vec(
            reversed_bodies[i].angular_momentum_body(),
            original[i].angular_momentum_body().scale(-1.0),
            3e-17,
        );
    }
    close(
        reverse_report.dissipated_heat_j,
        original_report.dissipated_heat_j,
        3e-18,
    );
}

#[test]
fn source_geometry_no_slip_spin_and_static_traction_boundaries_refuse() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let weak = FrictionLaw::Coulomb {
        static_mu: 0.0,
        kinetic_mu: 0.0,
    };
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let contact = RigidRollingContact {
        static_friction: &weak,
        ..f.contact()
    };
    assert!(matches!(
        apply_rigid_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            contact,
            0.2,
            controls(0.2),
            &CancelGate::new()
        ),
        Err(RigidRollingError::StaticTraction { .. })
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, RigidRollingHistory::new());
    let mismatch = RollingPatchReceipt::new(
        f.patch.patch_id(),
        f.patch.normal_model_id(),
        f.patch.source_id(),
        f.patch.authority(),
        f.patch.normal_load_n() * 1.01,
        f.patch.contact_area_m2(),
        f.patch.curvature(),
    )
    .unwrap();
    for contact in [
        RigidRollingContact {
            temperature_k: 321.0,
            ..f.contact()
        },
        RigidRollingContact {
            patch: &mismatch,
            ..f.contact()
        },
        RigidRollingContact {
            binding: RigidContactBinding {
                arm_a_body: Vec3::new(0.0, 0.0, -R * 0.9),
                ..f.binding
            },
            ..f.contact()
        },
    ] {
        assert!(
            apply_rigid_rolling(
                &mut bodies,
                &f.properties,
                &mut history,
                contact,
                0.2,
                controls(0.2),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
        assert_eq!(history, RigidRollingHistory::new());
    }
    for (linear, angular) in [
        (
            f.bodies[0]
                .linear_momentum_world()
                .add(Vec3::new(1e-3, 0.0, 0.0)),
            f.bodies[0].angular_momentum_body(),
        ),
        (
            f.bodies[0].linear_momentum_world(),
            f.bodies[0]
                .angular_momentum_body()
                .add(Vec3::new(0.0, 0.0, 1e-3)),
        ),
    ] {
        let mut altered = f.bodies;
        altered[0] = RigidBodyState::new(altered[0].pose(), linear, angular).unwrap();
        let initial = altered;
        assert!(
            apply_rigid_rolling(
                &mut altered,
                &f.properties,
                &mut history,
                f.contact(),
                0.2,
                controls(0.2),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(altered, initial);
    }
}

#[test]
fn loaded_late_audit_budget_source_and_cancel_failures_preserve_exact_retry() {
    let f = fixture(TAU, INITIAL_RATE, 0.001);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    advance(&f, &mut bodies, &mut history, 0.2, controls(0.1));
    assert!(history.state().cumulative_heat_j() > 0.0);
    let restored =
        RigidRollingHistory::from_checkpoint(f.contact(), history.checkpoint().unwrap().clone())
            .unwrap();
    assert_eq!(restored, history);
    // A proximity-admitted point mismatch produces an independently measurable
    // angular impulse defect only AFTER candidate impulses/history were staged.
    bodies[1] = RigidBodyState::new(
        Pose::new(
            bodies[1]
                .pose()
                .position_world()
                .add(Vec3::new(1e-7, 0.0, 0.0)),
            bodies[1].pose().orientation(),
        )
        .unwrap(),
        bodies[1].linear_momentum_world(),
        bodies[1].angular_momentum_body(),
    )
    .unwrap();
    let accepted = bodies;
    let accepted_history = history.clone();
    let mut loose = controls(0.1);
    loose.geometry_tolerance_m = 1e-6;
    loose.angular_momentum_tolerance_nms = 1e-9;
    let mut strict = loose;
    strict.angular_momentum_tolerance_nms = 1e-22;
    assert!(matches!(
        apply_rigid_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            0.3,
            strict,
            &CancelGate::new()
        ),
        Err(RigidRollingError::Audit(_))
    ));
    assert_eq!(bodies, accepted);
    assert_eq!(history, accepted_history);
    let mut short = loose;
    short.maximum_substeps = 1;
    assert!(matches!(
        apply_rigid_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            20.0,
            short,
            &CancelGate::new()
        ),
        Err(RigidRollingError::SubstepBudget { .. })
    ));
    let replacement = rolling_law(
        TAU,
        DeformingSphereSurface::SurfaceA,
        "fixture/replacement-source",
    );
    assert!(matches!(
        apply_rigid_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            RigidRollingContact {
                law: &replacement,
                ..f.contact()
            },
            0.3,
            loose,
            &CancelGate::new()
        ),
        Err(RigidRollingError::Native(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_rolling(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            0.3,
            loose,
            &cancelled
        )
        .is_err()
    );
    assert_eq!(bodies, accepted);
    assert_eq!(history, accepted_history);
    let mut reference_bodies = accepted;
    let mut reference_history = accepted_history;
    let expected = advance(
        &f,
        &mut reference_bodies,
        &mut reference_history,
        0.3,
        loose,
    );
    let actual = advance(&f, &mut bodies, &mut history, 0.3, loose);
    assert_eq!(bodies, reference_bodies);
    assert_eq!(history, reference_history);
    assert_eq!(actual, expected);
}
