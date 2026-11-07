//! Real prefix motion, work and explicit handoff at material-rolling breakaway.

use super::*;
use fs_mbd::sphere_plane::SpherePlaneEventKind;
use fs_mbd::sphere_plane::rolling::advance_viscoelastic_sphere_rolling_until_breakaway;

fn driven_fixture() -> (Fixture, FrictionLaw, f64) {
    let mut fixture = Fixture::new(2e-4);
    let drive = 0.008 * fixture.beta();
    fixture.scene.force_world = fixture.scene.force_world.add(Vec3::new(drive, 0.0, 0.0));
    let support = FrictionLaw::Coulomb {
        static_mu: 0.6 * drive / fixture.patch.normal_load_n(),
        kinetic_mu: 0.3 * drive / fixture.patch.normal_load_n(),
    };
    (fixture, support, drive)
}

#[test]
fn force_driven_breakaway_returns_the_exact_physical_prefix_even_for_a_long_request() {
    let (f, support, drive) = driven_fixture();
    let initial = f.state(Vec3::ZERO, 0.0, UnitQuaternion::IDENTITY);
    // Initial traction is -2F/7; the material drag makes it tend to -F.
    // At 0.6F capacity, exp(-rate*t)=(1-0.6)/(5/7)=0.56.
    let expected_time = -0.56_f64.ln() / f.rate();
    let requested = 1000.0 / f.rate();
    let result = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        initial,
        Vec3::ZERO,
        &support,
        f.binding(),
        requested,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    let event = result.breakaway.unwrap();
    close(result.accepted_duration_s, expected_time);
    close(result.remaining_duration_s, requested - expected_time);
    let h = result.accepted_duration_s;
    let terminal = drive / f.beta();
    let decay = (-f.rate() * h).exp();
    let j1 = -(-f.rate() * h).exp_m1() / f.rate();
    let j2 = -(-2.0 * f.rate() * h).exp_m1() / (2.0 * f.rate());
    let speed = terminal * (1.0 - decay);
    let distance = terminal * (h - j1);
    let heat = f.beta() * terminal * terminal * (h - 2.0 * j1 + j2);
    let rolling = result.rolling;
    let motion = rolling.motion;
    close(velocity(&f, motion.state_after).x, speed);
    close(omega(&f, motion.state_after).y, speed / RADIUS);
    close(
        motion.state_after.body().pose().position_world().x,
        distance,
    );
    close(rolling.deformation_heat_j, heat);
    close(motion.external_force_work, drive * distance);
    close(motion.external_torque_work, 0.0);
    close(motion.energy_balance_residual, 0.0);
    vector_close(motion.linear_momentum_residual, Vec3::ZERO);
    vector_close(motion.angular_momentum_residual, Vec3::ZERO);
    assert!(rolling.maximum_required_traction_n <= event.available_static_traction_n);
    vector_close(
        event.required_static_traction_world_n,
        Vec3::new(-0.6 * drive, 0.0, 0.0),
    );
    vector_close(
        event.rolling_moment_world_nm,
        Vec3::new(
            0.0,
            -f.card.rolling_time_s() * f.patch.normal_load_n() * speed,
            0.0,
        ),
    );
    assert!(rolling.deformation_heat_j > 0.0);
    assert_eq!(motion.state_after.mode(), SpherePlaneMode::Sliding);
    assert_eq!(motion.events.len(), 2);
    assert_eq!(motion.events[0].kind, SpherePlaneEventKind::Support);
    assert_eq!(motion.events[1].kind, SpherePlaneEventKind::Sliding);
    close(motion.events[1].time_s, expected_time);
    assert_eq!(motion.events[1].normal_impulse, 0.0);
    assert_eq!(motion.events[1].tangent_impulse_world, Vec3::ZERO);
    assert!(matches!(
        advance_viscoelastic_sphere_rolling(
            &f.scene,
            initial,
            &support,
            f.binding(),
            3.0 / f.rate(),
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::TractionExceeded { .. })
    ));
}

#[test]
fn combined_motor_breakaway_preserves_continuation_and_moving_frame_work() {
    let mut f = Fixture::new(2e-4);
    let drive = 0.008 * f.beta();
    let force = Vec3::new(drive, 0.0, 0.0);
    f.scene.force_world = f.scene.force_world.add(force);
    let inertial = Vec3::new(drive, 0.5 * drive, 0.0);
    // Newton/Euler balance at rest chooses traction (0,F/2). Its subsequent
    // two-axis path crosses the disk, not a componentwise force bound.
    let torque = f
        .scene
        .normal_world
        .cross(
            inertial
                .scale(f.effective_mass() / f.scene.properties.mass())
                .sub(force),
        )
        .scale(RADIUS)
        .add(
            f.scene
                .normal_world
                .scale(0.8 * f.scene.properties.principal_inertia_body().z),
        );
    let support = FrictionLaw::Coulomb {
        static_mu: 0.75 * drive / f.patch.normal_load_n(),
        kinetic_mu: 0.3 * drive / f.patch.normal_load_n(),
    };
    let initial = f.state(Vec3::ZERO, 0.6, UnitQuaternion::IDENTITY);
    let requested = 3.0 / f.rate();
    let expected_decay = (2.0 - 1.8125_f64.sqrt()) / 2.5;
    let expected_time = -expected_decay.ln() / f.rate();
    let controls = SpherePlaneConfig {
        max_attitude_step_s: 0.002,
        ..SpherePlaneConfig::default()
    };
    let whole = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        initial,
        torque,
        &support,
        f.binding(),
        requested,
        controls,
        || false,
    )
    .unwrap();
    close(whole.accepted_duration_s, expected_time);
    let split = 0.25 * expected_time;
    let first = advance_viscoelastic_sphere_rolling_with_torque(
        &f.scene,
        initial,
        torque,
        &support,
        f.binding(),
        split,
        controls,
        || false,
    )
    .unwrap();
    let second = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        first.motion.state_after,
        torque,
        &support,
        f.binding(),
        requested - split,
        controls,
        || false,
    )
    .unwrap();
    close(
        split + second.accepted_duration_s,
        whole.accepted_duration_s,
    );
    close(second.remaining_duration_s, whole.remaining_duration_s);
    vector_close(
        velocity(&f, second.rolling.motion.state_after),
        velocity(&f, whole.rolling.motion.state_after),
    );
    vector_close(
        omega(&f, second.rolling.motion.state_after),
        omega(&f, whole.rolling.motion.state_after),
    );
    vector_close(
        second
            .rolling
            .motion
            .state_after
            .body()
            .pose()
            .position_world(),
        whole
            .rolling
            .motion
            .state_after
            .body()
            .pose()
            .position_world(),
    );
    close(
        first.deformation_heat_j + second.rolling.deformation_heat_j,
        whole.rolling.deformation_heat_j,
    );
    close(
        first.motion.external_torque_work + second.rolling.motion.external_torque_work,
        whole.rolling.motion.external_torque_work,
    );
    close(
        first.motion.external_force_work + second.rolling.motion.external_force_work,
        whole.rolling.motion.external_force_work,
    );
    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let rotate = |value| rotation.rotate_body_to_world(value);
    let boost = Vec3::new(0.003, -0.002, 0.004);
    let mut moved = f.clone();
    moved.scene.normal_world = rotate(f.scene.normal_world);
    moved.scene.force_world = rotate(f.scene.force_world);
    moved.scene.velocity_world = boost;
    let result = advance_viscoelastic_sphere_rolling_until_breakaway(
        &moved.scene,
        moved.state(Vec3::ZERO, 0.6, rotation),
        rotate(torque),
        &support,
        moved.binding(),
        requested,
        controls,
        || false,
    )
    .unwrap();
    close(result.accepted_duration_s, whole.accepted_duration_s);
    vector_close(
        velocity(&moved, result.rolling.motion.state_after),
        rotate(velocity(&f, whole.rolling.motion.state_after)).add(boost),
    );
    vector_close(
        result.breakaway.unwrap().required_static_traction_world_n,
        rotate(whole.breakaway.unwrap().required_static_traction_world_n),
    );
    vector_close(
        result.breakaway.unwrap().rolling_moment_world_nm,
        rotate(whole.breakaway.unwrap().rolling_moment_world_nm),
    );
    close(
        result.rolling.deformation_heat_j,
        whole.rolling.deformation_heat_j,
    );
    close(
        result.rolling.motion.external_torque_work,
        whole.rolling.motion.external_torque_work,
    );
    close(
        result.rolling.motion.plane_work,
        result.rolling.motion.contact_impulse_world.dot(boost),
    );
    close(result.rolling.motion.energy_balance_residual, 0.0);
    vector_close(result.rolling.motion.angular_momentum_residual, Vec3::ZERO);
}

#[test]
fn immediate_breakaway_changes_mode_without_an_impulse_or_invented_motion() {
    let f = Fixture::new(2e-4);
    let initial = f.state(Vec3::ZERO, 0.0, UnitQuaternion::IDENTITY);
    let torque = Vec3::new(0.0, 0.8 * f.patch.normal_load_n() * RADIUS, 0.0);
    let result = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        initial,
        torque,
        &SUPPORT,
        f.binding(),
        0.1,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    let event = result.breakaway.unwrap();
    assert_eq!(result.accepted_duration_s, 0.0);
    assert_eq!(result.remaining_duration_s, 0.1);
    assert_eq!(result.rolling.motion.state_after.body(), initial.body());
    assert_eq!(result.rolling.motion.state_after.time_s(), initial.time_s());
    assert_eq!(
        result.rolling.motion.state_after.mode(),
        SpherePlaneMode::Sliding
    );
    assert_eq!(result.rolling.motion.contact_impulse_world, Vec3::ZERO);
    assert_eq!(
        result.rolling.motion.contact_angular_impulse_world,
        Vec3::ZERO
    );
    assert_eq!(
        result.rolling.motion.external_torque_impulse_world,
        Vec3::ZERO
    );
    assert_eq!(result.rolling.motion.external_torque_work, 0.0);
    assert_eq!(result.rolling.motion.external_force_work, 0.0);
    assert_eq!(result.rolling.motion.plane_work, 0.0);
    assert_eq!(result.rolling.deformation_heat_j, 0.0);
    assert!(
        event.required_static_traction_world_n.norm_squared().sqrt()
            > event.available_static_traction_n
    );
    let coasting = f.state(Vec3::new(0.005, 0.0, 0.0), 0.7, UnitQuaternion::IDENTITY);
    let result = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        coasting,
        Vec3::ZERO,
        &FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        f.binding(),
        0.1,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(result.accepted_duration_s, 0.0);
    assert_eq!(result.rolling.motion.state_after.body(), coasting.body());
    assert_eq!(result.rolling.deformation_heat_j, 0.0);
    assert!(result.rolling.response_before.dissipated_power_w > 0.0);
    assert!(result.breakaway.unwrap().rolling_moment_world_nm.y < 0.0);
    assert!(matches!(
        advance_viscoelastic_sphere_rolling_until_breakaway(
            &f.scene,
            initial,
            torque,
            &SUPPORT,
            f.binding(),
            0.1,
            SpherePlaneConfig {
                max_events: 1,
                ..SpherePlaneConfig::default()
            },
            || false,
        ),
        Err(SphereRollingError::Dynamics(SpherePlaneError::EventBudget))
    ));
}

#[test]
fn feasible_coast_and_elastic_motion_complete_without_a_false_breakaway() {
    for relaxation in [0.0, 2e-4] {
        let f = Fixture::new(relaxation);
        let initial = f.state(Vec3::new(0.005, 0.002, 0.0), 0.7, UnitQuaternion::IDENTITY);
        let result = advance_viscoelastic_sphere_rolling_until_breakaway(
            &f.scene,
            initial,
            Vec3::ZERO,
            &SUPPORT,
            f.binding(),
            0.2,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        let original = advance_viscoelastic_sphere_rolling(
            &f.scene,
            initial,
            &SUPPORT,
            f.binding(),
            0.2,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        assert!(result.breakaway.is_none());
        assert_eq!(result.accepted_duration_s, 0.2);
        assert_eq!(result.remaining_duration_s, 0.0);
        assert_eq!(
            result.rolling.motion.state_after,
            original.motion.state_after
        );
        assert_eq!(
            result.rolling.deformation_heat_j,
            original.deformation_heat_j
        );
        assert_eq!(
            result.rolling.motion.state_after.mode(),
            SpherePlaneMode::Rolling
        );
    }
}

#[test]
fn source_limits_and_final_cancellation_refuse_the_entire_uncommitted_prefix() {
    let (f, support, _) = driven_fixture();
    let initial = f.state(Vec3::ZERO, 0.0, UnitQuaternion::IDENTITY);
    let mut polls = 0;
    let accepted = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        initial,
        Vec3::ZERO,
        &support,
        f.binding(),
        3.0 / f.rate(),
        SpherePlaneConfig::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let last_poll = polls;
    polls = 0;
    assert!(matches!(
        advance_viscoelastic_sphere_rolling_until_breakaway(
            &f.scene,
            initial,
            Vec3::ZERO,
            &support,
            f.binding(),
            3.0 / f.rate(),
            SpherePlaneConfig::default(),
            || {
                polls += 1;
                polls == last_poll
            },
        ),
        Err(SphereRollingError::Dynamics(SpherePlaneError::Cancelled))
    ));
    let retry = advance_viscoelastic_sphere_rolling_until_breakaway(
        &f.scene,
        initial,
        Vec3::ZERO,
        &support,
        f.binding(),
        3.0 / f.rate(),
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(
        retry.rolling.motion.state_after,
        accepted.rolling.motion.state_after
    );
    assert_eq!(
        retry.rolling.deformation_heat_j,
        accepted.rolling.deformation_heat_j
    );
    assert_eq!(initial.time_s(), 0.0);
    let mut outside = f.clone();
    let drive = 0.2 * outside.beta();
    outside.scene.force_world = Vec3::new(drive, 0.0, -outside.patch.normal_load_n());
    let weak = FrictionLaw::Coulomb {
        static_mu: 0.6 * drive / outside.patch.normal_load_n(),
        kinetic_mu: 0.3 * drive / outside.patch.normal_load_n(),
    };
    // The original material's speed/relaxation domain ends before traction
    // breakaway. An event locator cannot authorize that invalid prefix.
    assert!(matches!(
        advance_viscoelastic_sphere_rolling_until_breakaway(
            &outside.scene,
            initial,
            Vec3::ZERO,
            &weak,
            outside.binding(),
            3.0 / outside.rate(),
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::Constitutive(
            RollingLossError::OutsideApplicability { .. }
        ))
    ));
    let mut hot = f.binding();
    hot.temperature_kelvin = 330.0;
    assert!(matches!(
        advance_viscoelastic_sphere_rolling_until_breakaway(
            &f.scene,
            initial,
            Vec3::new(0.0, f.patch.normal_load_n(), 0.0),
            &support,
            hot,
            0.1,
            SpherePlaneConfig::default(),
            || false,
        ),
        Err(SphereRollingError::Constitutive(
            RollingLossError::OutsideApplicability { .. }
        ))
    ));
}
