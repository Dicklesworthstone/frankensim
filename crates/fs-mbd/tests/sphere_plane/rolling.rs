//! Material-derived motion and actual deformation work; synthetic source data.
use fs_mbd::sphere_plane::rolling::{
    SphereRollingBinding, SphereRollingError, advance_viscoelastic_sphere_rolling,
};
use fs_mbd::sphere_plane::{
    SpherePlane, SpherePlaneConfig, SpherePlaneError, SpherePlaneMode, SpherePlaneState,
    advance_sphere_plane_with_support_friction,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::rolling_loss::viscoelastic::{
    DeformingSphereSurface, QuasistaticSphereLimits, ViscoelasticSphereParameters,
    ViscoelasticSphereRollingCard,
};
use fs_tribo::rolling_loss::{
    PatchCurvature, RollingLossApplicability, RollingLossChannel, RollingLossError,
    RollingPatchReceipt, RollingWorkOwnership,
};
use fs_tribo::{
    ApplicabilityRange, FrictionLaw, HertzSpherePlane, InputAuthority, InterfaceMedium,
    InterfaceSystemRef,
};

const SUPPORT: FrictionLaw = FrictionLaw::Coulomb {
    static_mu: 0.4,
    kinetic_mu: 0.2,
};
const RADIUS: f64 = 0.02;

#[derive(Clone)]
struct Fixture {
    scene: SpherePlane,
    card: ViscoelasticSphereRollingCard,
    patch: RollingPatchReceipt,
    interface: InterfaceSystemRef,
    ownership: RollingWorkOwnership,
}

impl Fixture {
    fn new(relaxation: f64) -> Self {
        let young = 2.0e7;
        let poisson = 0.3;
        let density = 1100.0;
        let mass = (4.0 * core::f64::consts::PI / 3.0) * RADIUS.powi(3) * density;
        let inertia = 0.4 * mass * RADIUS.powi(2);
        let parameters = ViscoelasticSphereParameters {
            sphere_surface: DeformingSphereSurface::SurfaceA,
            radius_m: RADIUS,
            young_modulus_pa: young,
            poisson_ratio: poisson,
            density_kg_m3: density,
            shear_viscosity_pa_s: relaxation * young / (2.0 * (1.0 + poisson)),
            bulk_viscosity_pa_s: relaxation * young / (3.0 * (1.0 - 2.0 * poisson)),
        };
        let card = ViscoelasticSphereRollingCard::new(
            "first-order material sphere",
            "synthetic/viscosities",
            InputAuthority::SyntheticFixture,
            parameters,
            QuasistaticSphereLimits {
                maximum_indentation_radius_ratio: 0.03,
                maximum_speed_wave_ratio: 0.01,
                maximum_relaxation_ratio: 0.05,
            },
            RollingLossApplicability::new(
                ApplicabilityRange::new(280.0, 320.0).unwrap(),
                ApplicabilityRange::new(0.0, 100.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
        let interface = InterfaceSystemRef::new(
            "sphere/plane",
            "rolling-material-history",
            "synthetic/sphere-material",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap();
        let normal = HertzSpherePlane {
            effective_radius: RADIUS,
            reduced_modulus: young / (1.0 - poisson * poisson),
        }
        .response(&interface, 2e-4)
        .unwrap();
        let patch = RollingPatchReceipt::new(
            "material-sphere-patch",
            "original-Hertz",
            "synthetic/normal",
            InputAuthority::SyntheticFixture,
            normal.normal_force_n,
            core::f64::consts::PI * normal.contact_radius_m.powi(2),
            PatchCurvature::Principal {
                first_per_m: 1.0 / RADIUS,
                second_per_m: 1.0 / RADIUS,
            },
        )
        .unwrap();
        let scene = SpherePlane {
            properties: MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia))
                .unwrap(),
            radius: RADIUS,
            point_world: Vec3::ZERO,
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            velocity_world: Vec3::ZERO,
            force_world: Vec3::new(0.0, 0.0, -normal.normal_force_n),
            restitution: 0.0,
            friction_coefficient: 0.4,
        };
        let ownership = RollingWorkOwnership::new(
            patch.patch_id(),
            "rolling-interval",
            "sphere-rolling-angle",
            RollingLossChannel::RollingHysteresis,
        )
        .unwrap();
        Self {
            scene,
            card,
            patch,
            interface,
            ownership,
        }
    }
    fn binding(&self) -> SphereRollingBinding<'_> {
        SphereRollingBinding {
            card: &self.card,
            patch: &self.patch,
            interface: &self.interface,
            ownership: &self.ownership,
            temperature_kelvin: 300.0,
            excitation_frequency_hz: 10.0,
        }
    }
    fn state(
        &self,
        relative_velocity: Vec3,
        spin: f64,
        attitude: UnitQuaternion,
    ) -> SpherePlaneState {
        let omega = self
            .scene
            .normal_world
            .cross(relative_velocity)
            .scale(RADIUS.recip())
            .add(self.scene.normal_world.scale(spin));
        SpherePlaneState::new(
            RigidBodyState::new(
                Pose::new(self.scene.normal_world.scale(RADIUS), attitude).unwrap(),
                self.scene
                    .velocity_world
                    .add(relative_velocity)
                    .scale(self.scene.properties.mass()),
                attitude.rotate_world_to_body(
                    omega.scale(self.scene.properties.principal_inertia_body().x),
                ),
            )
            .unwrap(),
            0.0,
        )
        .unwrap()
    }
    fn effective_mass(&self) -> f64 {
        1.4 * self.scene.properties.mass()
    }
    fn beta(&self) -> f64 {
        self.card.rolling_time_s() * self.patch.normal_load_n() / RADIUS
    }
    fn rate(&self) -> f64 {
        self.beta() / self.effective_mass()
    }
}

#[track_caller]
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1e-6),
        "{a:e} != {b:e}"
    );
}
#[track_caller]
fn vector_close(a: Vec3, b: Vec3) {
    // Relative vector accuracy must be invariant under the tested rotation;
    // a nearly cancelling component is not an independent magnitude scale.
    let scale = a
        .norm_squared()
        .sqrt()
        .max(b.norm_squared().sqrt())
        .max(1e-6);
    assert!(
        a.sub(b).norm_squared().sqrt() <= 2e-11 * scale,
        "{a:?} != {b:?}"
    );
}
fn velocity(f: &Fixture, state: SpherePlaneState) -> Vec3 {
    state
        .body()
        .center_of_mass_velocity_world(f.scene.properties)
        .unwrap()
}
fn omega(f: &Fixture, state: SpherePlaneState) -> Vec3 {
    let body = state.body();
    body.pose().orientation().rotate_body_to_world(
        f.scene
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())
            .unwrap(),
    )
}

#[test]
fn material_coast_down_has_exact_exponential_motion_distance_and_bulk_heat() {
    let f = Fixture::new(2e-4);
    let speed = 0.01;
    let before = f.state(Vec3::new(speed, 0.0, 0.0), 0.0, UnitQuaternion::IDENTITY);
    for scaled_time in [1e-8, 0.2, 2.0, 20.0] {
        let h = scaled_time / f.rate();
        let result = advance_viscoelastic_sphere_rolling(
            &f.scene,
            before,
            &SUPPORT,
            f.binding(),
            h,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        let decay = (-scaled_time).exp();
        let distance = speed * -(-scaled_time).exp_m1() / f.rate();
        let heat = 0.5 * f.effective_mass() * speed * speed * -(-2.0 * scaled_time).exp_m1();
        let after = result.motion.state_after;
        close(velocity(&f, after).x, speed * decay);
        assert!(velocity(&f, after).x > 0.0, "no finite-time viscous arrest");
        close(after.body().pose().position_world().x, distance);
        close(omega(&f, after).y, speed * decay / RADIUS);
        close(result.deformation_heat_j, heat);
        close(result.motion.dissipated_energy, heat);
        close(
            result.maximum_required_traction_n,
            f.scene.properties.mass() * f.rate() * speed,
        );
        close(
            result.rolling_angular_impulse_world.y,
            -f.card.rolling_time_s() * f.patch.normal_load_n() * distance,
        );
        close(result.motion.energy_balance_residual, 0.0);
        vector_close(result.motion.linear_momentum_residual, Vec3::ZERO);
        vector_close(result.motion.angular_momentum_residual, Vec3::ZERO);
        let rotation =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), distance / RADIUS).unwrap();
        vector_close(
            after
                .body()
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        );
        assert_eq!(after.mode(), SpherePlaneMode::Rolling);
        assert_eq!(result.binding.card.source_id(), f.card.source_id());
    }
}

#[test]
fn driven_turn_and_reversal_integrate_real_force_and_continuous_torque_work() {
    for terminal in [Vec3::new(0.0, 0.006, 0.0), Vec3::new(-0.005, 0.0, 0.0)] {
        let mut f = Fixture::new(2e-4);
        let tangent_force = terminal.scale(f.beta());
        f.scene.force_world = f.scene.force_world.add(tangent_force);
        let u0 = Vec3::new(0.01, 0.0, 0.0);
        let before = f.state(u0, 0.0, UnitQuaternion::IDENTITY);
        let h = 3.0 / f.rate();
        let result = advance_viscoelastic_sphere_rolling(
            &f.scene,
            before,
            &SUPPORT,
            f.binding(),
            h,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        let difference = u0.sub(terminal);
        let j1 = (1.0 - (-3.0_f64).exp()) / f.rate();
        let j2 = (1.0 - (-6.0_f64).exp()) / (2.0 * f.rate());
        let expected_velocity = terminal.add(difference.scale((-3.0_f64).exp()));
        let displacement = terminal.scale(h).add(difference.scale(j1));
        // Integrate |u_inf + (u0-u_inf) exp(-rate*t)|² directly.
        let heat = f.beta()
            * (terminal.dot(terminal) * h
                + 2.0 * terminal.dot(difference) * j1
                + difference.dot(difference) * j2);
        vector_close(velocity(&f, result.motion.state_after), expected_velocity);
        vector_close(
            result.motion.state_after.body().pose().position_world(),
            before.body().pose().position_world().add(displacement),
        );
        close(result.deformation_heat_j, heat);
        close(
            result.motion.external_force_work,
            tangent_force.dot(displacement),
        );
        close(result.motion.energy_balance_residual, 0.0);
        vector_close(result.motion.angular_momentum_residual, Vec3::ZERO);
        if terminal.x < 0.0 {
            assert!(velocity(&f, result.motion.state_after).x < 0.0);
        } else {
            assert!(velocity(&f, result.motion.state_after).y > 0.0);
        }
        assert!(result.deformation_heat_j > 0.0);
    }
}

#[test]
fn translating_support_rotated_world_and_split_continuation_preserve_work_and_spin() {
    let mut base = Fixture::new(2e-4);
    base.scene.force_world = base
        .scene
        .force_world
        .add(Vec3::new(0.0, 0.004 * base.beta(), 0.0));
    let initial = base.state(Vec3::new(0.01, 0.0, 0.0), 1.3, UnitQuaternion::IDENTITY);
    let h = 0.6;
    let controls = SpherePlaneConfig {
        max_attitude_step_s: 0.002,
        ..SpherePlaneConfig::default()
    };
    let whole = advance_viscoelastic_sphere_rolling(
        &base.scene,
        initial,
        &SUPPORT,
        base.binding(),
        h,
        controls,
        || false,
    )
    .unwrap();
    let first = advance_viscoelastic_sphere_rolling(
        &base.scene,
        initial,
        &SUPPORT,
        base.binding(),
        0.5 * h,
        controls,
        || false,
    )
    .unwrap();
    let second = advance_viscoelastic_sphere_rolling(
        &base.scene,
        first.motion.state_after,
        &SUPPORT,
        base.binding(),
        0.5 * h,
        controls,
        || false,
    )
    .unwrap();
    vector_close(
        velocity(&base, whole.motion.state_after),
        velocity(&base, second.motion.state_after),
    );
    vector_close(
        whole.motion.state_after.body().pose().position_world(),
        second.motion.state_after.body().pose().position_world(),
    );
    vector_close(
        omega(&base, whole.motion.state_after),
        omega(&base, second.motion.state_after),
    );
    close(
        whole.deformation_heat_j,
        first.deformation_heat_j + second.deformation_heat_j,
    );
    close(omega(&base, whole.motion.state_after).z, 1.3);
    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let rotate = |v| rotation.rotate_body_to_world(v);
    let boost = Vec3::new(0.003, -0.002, 0.004);
    let mut moved = base.clone();
    moved.scene.normal_world = rotate(base.scene.normal_world);
    moved.scene.force_world = rotate(base.scene.force_world);
    moved.scene.velocity_world = boost;
    let moved_initial = moved.state(rotate(Vec3::new(0.01, 0.0, 0.0)), 1.3, rotation);
    let result = advance_viscoelastic_sphere_rolling(
        &moved.scene,
        moved_initial,
        &SUPPORT,
        moved.binding(),
        h,
        controls,
        || false,
    )
    .unwrap();
    vector_close(
        velocity(&moved, result.motion.state_after),
        rotate(velocity(&base, whole.motion.state_after)).add(boost),
    );
    vector_close(
        omega(&moved, result.motion.state_after),
        rotate(omega(&base, whole.motion.state_after)),
    );
    vector_close(
        result.motion.state_after.body().pose().position_world(),
        rotate(whole.motion.state_after.body().pose().position_world()).add(boost.scale(h)),
    );
    vector_close(
        result.motion.contact_impulse_world,
        rotate(whole.motion.contact_impulse_world),
    );
    close(result.deformation_heat_j, whole.deformation_heat_j);
    close(
        result.motion.plane_work,
        result.motion.contact_impulse_world.dot(boost),
    );
    close(
        result.motion.external_force_work,
        whole.motion.external_force_work + moved.scene.force_world.dot(boost.scale(h)),
    );
    close(result.motion.energy_balance_residual, 0.0);
    vector_close(result.motion.angular_momentum_residual, Vec3::ZERO);
    let axis = Vec3::new(0.0, 1.0, 0.0);
    vector_close(
        result
            .motion
            .state_after
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(axis),
        rotate(
            whole
                .motion
                .state_after
                .body()
                .pose()
                .orientation()
                .rotate_body_to_world(axis),
        ),
    );
}

#[test]
fn elastic_and_nearly_elastic_limits_retain_original_driven_rolling() {
    for relaxation in [0.0, 1e-16] {
        let mut f = Fixture::new(relaxation);
        f.scene.force_world = f.scene.force_world.add(Vec3::new(1e-4, 2e-4, 0.0));
        let initial = f.state(Vec3::ZERO, 0.0, UnitQuaternion::IDENTITY);
        let h = 0.2;
        let result = advance_viscoelastic_sphere_rolling(
            &f.scene,
            initial,
            &SUPPORT,
            f.binding(),
            h,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        let original = advance_sphere_plane_with_support_friction(
            &f.scene,
            initial,
            &SUPPORT,
            h,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        vector_close(
            velocity(&f, result.motion.state_after),
            velocity(&f, original.state_after),
        );
        vector_close(
            result.motion.state_after.body().pose().position_world(),
            original.state_after.body().pose().position_world(),
        );
        vector_close(
            omega(&f, result.motion.state_after),
            omega(&f, original.state_after),
        );
        close(
            result.motion.external_force_work,
            original.external_force_work,
        );
        if relaxation == 0.0 {
            assert_eq!(result.deformation_heat_j, 0.0);
        } else {
            assert!(result.deformation_heat_j > 0.0);
        }
    }
    let mut released = Fixture::new(2e-4);
    released.scene.force_world = Vec3::ZERO;
    released.patch = RollingPatchReceipt::new(
        released.patch.patch_id(),
        released.patch.normal_model_id(),
        released.patch.source_id(),
        released.patch.authority(),
        0.0,
        0.0,
        released.patch.curvature(),
    )
    .unwrap();
    let initial = released.state(Vec3::new(0.01, 0.0, 0.0), 0.0, UnitQuaternion::IDENTITY);
    let result = advance_viscoelastic_sphere_rolling(
        &released.scene,
        initial,
        &SUPPORT,
        released.binding(),
        0.2,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(result.deformation_heat_j, 0.0);
    vector_close(
        velocity(&released, result.motion.state_after),
        velocity(&released, initial),
    );
    assert_eq!(result.motion.contact_impulse_world, Vec3::ZERO);
}

#[test]
fn no_slip_feasibility_source_identity_and_late_cancellation_refuse_atomically() {
    let mut f = Fixture::new(2e-4);
    let drive = 0.008 * f.beta();
    f.scene.force_world = f.scene.force_world.add(Vec3::new(drive, 0.0, 0.0));
    let initial = f.state(Vec3::ZERO, 0.0, UnitQuaternion::IDENTITY);
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6 * drive / f.patch.normal_load_n(),
        kinetic_mu: 0.0,
    };
    // Initial reaction is (2/7)*drive, but tends to the whole drive. The later
    // cone failure must refuse even though this phase starts feasibly at rest.
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            initial,
            &law,
            f.binding(),
            3.0 / f.rate(),
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::TractionExceeded { .. })
    ));
    let mut hot = f.binding();
    hot.temperature_kelvin = 330.0;
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            initial,
            &SUPPORT,
            hot,
            0.2,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::Constitutive(
            RollingLossError::OutsideApplicability { .. }
        ))
    ));
    let mut wrong = f.scene;
    let i = wrong.properties.principal_inertia_body().x * 0.9;
    wrong.properties =
        MassProperties::new(wrong.properties.mass(), Vec3::ZERO, Vec3::new(i, i, i)).unwrap();
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &wrong,
            initial,
            &SUPPORT,
            f.binding(),
            0.2,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::Dynamics(SpherePlaneError::Invalid(_)))
    ));
    let sliding = SpherePlaneState::new(
        RigidBodyState::new(initial.body().pose(), Vec3::new(1e-3, 0.0, 0.0), Vec3::ZERO).unwrap(),
        0.0,
    )
    .unwrap();
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            sliding,
            &SUPPORT,
            f.binding(),
            0.2,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::Dynamics(SpherePlaneError::Invalid(_)))
    ));
    // Unrelated normal spin must not raise the no-slip tolerance enough to
    // conceal resolved tangential motion or its missing surface-friction heat.
    let spinning_slip = SpherePlaneState::new(
        RigidBodyState::new(
            initial.body().pose(),
            sliding.body().linear_momentum_world(),
            Vec3::new(
                0.0,
                0.0,
                1e14 * f.scene.properties.principal_inertia_body().z,
            ),
        )
        .unwrap(),
        0.0,
    )
    .unwrap();
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            spinning_slip,
            &SUPPORT,
            f.binding(),
            0.2,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::Dynamics(SpherePlaneError::Invalid(_)))
    ));
    let rotating = f.state(Vec3::new(0.0, 0.005, 0.0), 1.0, UnitQuaternion::IDENTITY);
    let mut polls = 0;
    let accepted = advance_viscoelastic_sphere_rolling(
        &f.scene,
        rotating,
        &SUPPORT,
        f.binding(),
        0.2,
        SpherePlaneConfig::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let final_poll = polls;
    polls = 0;
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            rotating,
            &SUPPORT,
            f.binding(),
            0.2,
            SpherePlaneConfig::default(),
            || {
                polls += 1;
                polls == final_poll
            },
        ),
        Err(SphereRollingError::Dynamics(SpherePlaneError::Cancelled))
    ));
    let retry = advance_viscoelastic_sphere_rolling(
        &f.scene,
        rotating,
        &SUPPORT,
        f.binding(),
        0.2,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(retry.motion.state_after, accepted.motion.state_after);
    assert_eq!(retry.deformation_heat_j, accepted.deformation_heat_j);
    assert_eq!(rotating.time_s(), 0.0);
}

#[test]
fn changing_axis_attitude_refines_without_changing_analytic_heat_or_momenta() {
    let mut f = Fixture::new(2e-4);
    f.scene.force_world = f
        .scene
        .force_world
        .add(Vec3::new(0.0, 0.006 * f.beta(), 0.0));
    let before = f.state(Vec3::new(0.01, 0.0, 0.0), 1.3, UnitQuaternion::IDENTITY);
    let h = 0.8;
    let solve = |step| {
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            before,
            &SUPPORT,
            f.binding(),
            h,
            SpherePlaneConfig {
                max_attitude_step_s: step,
                ..SpherePlaneConfig::default()
            },
            || false,
        )
        .unwrap()
    };
    let reference = solve(0.0002);
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let expected = reference
        .motion
        .state_after
        .body()
        .pose()
        .orientation()
        .rotate_body_to_world(axis);
    let mut errors = Vec::new();
    for step in [0.04, 0.02, 0.01] {
        let result = solve(step);
        let actual = result
            .motion
            .state_after
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(axis);
        errors.push(actual.sub(expected).norm_squared().sqrt());
        vector_close(
            velocity(&f, result.motion.state_after),
            velocity(&f, reference.motion.state_after),
        );
        vector_close(
            omega(&f, result.motion.state_after),
            omega(&f, reference.motion.state_after),
        );
        close(result.deformation_heat_j, reference.deformation_heat_j);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
}
