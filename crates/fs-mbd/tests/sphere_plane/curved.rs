//! G1/G3/G4: resolved changing-direction Coulomb slip on native spheres.
use super::*;

fn curved_scene() -> SpherePlane {
    SpherePlane {
        force_world: Vec3::new(0.0, -20.0, 2.0),
        // kappa=1/m+r²/I=7/4; mu*N=8/7, hence slip deceleration B=2.
        friction_coefficient: 2.0 / 35.0,
        ..scene()
    }
}

fn initial() -> SpherePlaneState {
    state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    )
}

fn controls(step: f64) -> SpherePlaneConfig {
    SpherePlaneConfig {
        max_sliding_step_s: step,
        max_attitude_step_s: 0.001,
        ..SpherePlaneConfig::default()
    }
}

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn slip(scene: &SpherePlane, body: RigidBodyState) -> Vec3 {
    velocity(body)
        .sub(scene.velocity_world)
        .add(spin(body).cross(scene.normal_world.scale(-RADIUS)))
}

#[test]
fn g1_curved_sliding_converges_to_independent_exact_two_axis_flow() {
    let scene = curved_scene();
    let before = initial();
    // Independent continuum solution: u'=a-B*u/|u|, a=(0,1), B=2.
    // With w=tan(theta/2), u_cross=w², u_drive=(w-w³)/2 and
    // t=2/3-w/2-w³/6. At w=1/2, T=19/48. Integrating these
    // expressions gives integral(u_cross)=233/960, integral(u_drive)=27/512.
    // Native COM velocity is (5/7)*(v0+a*t)+(2/7)*u; spin follows r×J.
    let duration = 19.0 / 48.0;
    let expected_v = Vec3::new(11.0 / 14.0, 0.0, 113.0 / 336.0);
    let expected_x = Vec3::new(1183.0 / 3360.0, RADIUS, 2291.0 / 32256.0);
    let expected_spin = Vec3::new(25.0 / 84.0, 0.0, -15.0 / 14.0);
    // Independently mu*N*integral(|u|dt) = (8/7)*(129/512).
    let expected_loss = 129.0 / 448.0;
    let mut errors = [[0.0; 4]; 3];
    for (rung, count) in [16, 32, 64].into_iter().enumerate() {
        let report = advance_sphere_plane(
            &scene,
            before,
            duration,
            controls(duration / count as f64),
            || false,
        )
        .unwrap();
        let body = report.state_after.body();
        assert_eq!(report.state_after.mode(), SpherePlaneMode::Sliding);
        assert!(report.sliding_substeps >= count);
        assert!(
            !report
                .events
                .iter()
                .any(|e| e.kind == SpherePlaneEventKind::Rolling)
        );
        balances(&scene, before, duration, &report);
        errors[rung] = [
            norm(velocity(body).sub(expected_v)),
            norm(body.pose().position_world().sub(expected_x)),
            norm(spin(body).sub(expected_spin)),
            (report.dissipated_energy - expected_loss).abs(),
        ];
    }
    for observable in 0..4 {
        for rung in 0..2 {
            let ratio = errors[rung][observable] / errors[rung + 1][observable];
            assert!((1.8..2.2).contains(&ratio), "continuum errors={errors:?}");
        }
    }
    assert!(errors[2][0] < 0.002);
    assert!(errors[2][1] < 0.0005);
    assert!(errors[2][2] < 0.01);
    assert!(errors[2][3] < 0.001);
}

#[test]
fn g1_curved_slip_arrest_converges_then_preserves_driven_rolling() {
    let scene = curved_scene();
    let before = initial();
    // The same exact solution reaches w=0 at t=2/3 s. Thereafter rolling
    // accelerates along z at 5/7 m/s², with no additional friction loss.
    let expected_x = Vec3::new(83.0 / 105.0, RADIUS, 8.0 / 21.0);
    let expected_v = Vec3::new(5.0 / 7.0, 0.0, 5.0 / 7.0);
    let mut errors = [[0.0; 2]; 3];
    for (rung, count) in [32, 64, 128].into_iter().enumerate() {
        let h = 1.0 / count as f64;
        let report = advance_sphere_plane(&scene, before, 1.0, controls(h), || false).unwrap();
        let body = report.state_after.body();
        let arrest = report
            .events
            .iter()
            .find(|e| e.kind == SpherePlaneEventKind::Rolling)
            .unwrap();
        assert!((arrest.time_s - 2.0 / 3.0).abs() <= 3.0 * h);
        assert_eq!(report.state_after.mode(), SpherePlaneMode::Rolling);
        close_vector(slip(&scene, body), Vec3::ZERO);
        close_vector(velocity(body), expected_v);
        balances(&scene, before, 1.0, &report);
        errors[rung] = [
            norm(body.pose().position_world().sub(expected_x)),
            (report.dissipated_energy - 1.0 / 3.0).abs(),
        ];
    }
    for observable in 0..2 {
        for rung in 0..2 {
            let ratio = errors[rung][observable] / errors[rung + 1][observable];
            assert!((1.7..2.3).contains(&ratio), "arrest errors={errors:?}");
        }
    }
}

#[test]
fn g1_coarse_turning_steps_remain_passive_and_supercritical_drive_cannot_stick() {
    let mut scene = curved_scene();
    scene.force_world = Vec3::new(-3.0, -20.0, 0.2);
    let before = initial();
    // A single 1 s endpoint projection would accelerate the INITIAL slip.
    // Resolving that rejected step gives a passive half-step followed by roll.
    let report = advance_sphere_plane(&scene, before, 1.0, controls(1.0), || false).unwrap();
    assert_eq!(report.sliding_substeps, 2);
    assert_eq!(report.state_after.mode(), SpherePlaneMode::Rolling);
    assert!(report.dissipated_energy > 0.0);
    balances(&scene, before, 1.0, &report);

    // A=3 exceeds B=2. Although a coarse implicit step can put u on the disk's
    // apex, transverse continuous slip cannot arrest under this drive.
    scene.force_world = Vec3::new(0.0, -20.0, 6.0);
    let before = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(0.01, 0.0, -1.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let report = advance_sphere_plane(&scene, before, 1.0, controls(0.1), || false).unwrap();
    assert_eq!(report.state_after.mode(), SpherePlaneMode::Sliding);
    assert!(
        !report
            .events
            .iter()
            .any(|e| e.kind == SpherePlaneEventKind::Rolling)
    );
    assert!(report.sliding_substeps > 10);
    assert!(slip(&scene, report.state_after.body()).z > 0.7);
    balances(&scene, before, 1.0, &report);
}

#[test]
fn g3_curved_motion_preserves_moving_plane_work_continuation_and_atomic_refusal() {
    let scene = curved_scene();
    let before = initial();
    let duration = 0.5;
    let config = controls(1.0 / 128.0);
    let whole = advance_sphere_plane(&scene, before, duration, config, || false).unwrap();
    let first = advance_sphere_plane(&scene, before, duration / 2.0, config, || false).unwrap();
    let second =
        advance_sphere_plane(&scene, first.state_after, duration / 2.0, config, || false).unwrap();
    close_vector(
        second.state_after.body().pose().position_world(),
        whole.state_after.body().pose().position_world(),
    );
    close_vector(
        velocity(second.state_after.body()),
        velocity(whole.state_after.body()),
    );
    close_vector(
        spin(second.state_after.body()),
        spin(whole.state_after.body()),
    );
    close(
        first.dissipated_energy + second.dissipated_energy,
        whole.dissipated_energy,
    );
    close_vector(
        first
            .contact_impulse_world
            .add(second.contact_impulse_world),
        whole.contact_impulse_world,
    );

    let boost = Vec3::new(0.3, 0.2, -0.15);
    let moving_scene = SpherePlane {
        velocity_world: boost,
        ..scene
    };
    let moving = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(1.0, 0.0, 0.0).add(boost),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let result = advance_sphere_plane(&moving_scene, moving, duration, config, || false).unwrap();
    close_vector(
        result
            .state_after
            .body()
            .pose()
            .position_world()
            .sub(whole.state_after.body().pose().position_world()),
        boost.scale(duration),
    );
    close_vector(
        velocity(result.state_after.body()).sub(velocity(whole.state_after.body())),
        boost,
    );
    close_vector(
        spin(result.state_after.body()),
        spin(whole.state_after.body()),
    );
    close(result.dissipated_energy, whole.dissipated_energy);
    close(result.plane_work, result.contact_impulse_world.dot(boost));
    close(
        result.external_force_work - whole.external_force_work,
        scene.force_world.dot(boost) * duration,
    );
    balances(&moving_scene, moving, duration, &result);

    assert!(matches!(
        advance_sphere_plane(
            &scene,
            before,
            duration,
            SpherePlaneConfig {
                max_sliding_substeps: 1,
                ..config
            },
            || false
        ),
        Err(SpherePlaneError::SlidingBudget)
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_sphere_plane(&scene, before, duration, config, || {
            polls += 1;
            polls == 20
        }),
        Err(SpherePlaneError::Cancelled)
    ));
    let retry = advance_sphere_plane(&scene, before, duration, config, || false).unwrap();
    assert_eq!(retry.state_after, whole.state_after);
}
