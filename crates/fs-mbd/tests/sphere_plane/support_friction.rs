//! G1/G3/G4: physical breakaway, kinetic slip and event-resolved recapture.
use super::*;
use fs_mbd::sphere_plane::advance_sphere_plane_with_support_friction;
use fs_tribo::FrictionLaw;

fn law() -> FrictionLaw {
    FrictionLaw::Coulomb {
        static_mu: 0.3,
        kinetic_mu: 0.1,
    }
}

fn supported(speed: f64) -> SpherePlaneState {
    state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(speed, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    )
}

#[test]
fn g1_static_threshold_and_kinetic_hysteresis_use_original_sphere_inertia() {
    // A solid sphere requires |F_t|/3.5 of static traction, not |F_t|.
    // N=20, mu_s=.3 therefore permit rolling through |F_t|=21 N;
    // an already sliding sphere feels only mu_k*N=2 N.
    let duration = 0.4;
    for direction in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.6, 0.0, 0.8)] {
        for force in [14.0, 21.0, 22.0] {
            let scene = SpherePlane {
                force_world: direction.scale(force).add(Vec3::new(0.0, -20.0, 0.0)),
                ..scene()
            };
            let before = supported(0.0);
            let report = advance_sphere_plane_with_support_friction(
                &scene,
                before,
                &law(),
                duration,
                SpherePlaneConfig::default(),
                || false,
            )
            .unwrap();
            let rolling = force <= 21.0;
            let acceleration = if rolling {
                force / 2.8
            } else {
                (force - 2.0) / MASS
            };
            assert_eq!(
                report.state_after.mode(),
                if rolling {
                    SpherePlaneMode::Rolling
                } else {
                    SpherePlaneMode::Sliding
                }
            );
            close_vector(
                velocity(report.state_after.body()),
                direction.scale(acceleration * duration),
            );
            close_vector(
                report.state_after.body().pose().position_world(),
                Vec3::new(0.0, RADIUS, 0.0)
                    .add(direction.scale(0.5 * acceleration * duration * duration)),
            );
            let slip_acceleration = if rolling { 0.0 } else { force / MASS - 3.5 };
            close(
                report.dissipated_energy,
                slip_acceleration * duration * duration,
            );
            balances(&scene, before, duration, &report);
        }
    }
    let scene = SpherePlane {
        force_world: Vec3::new(14.0, -20.0, 0.0),
        ..scene()
    };
    let before = supported(0.4);
    let report = advance_sphere_plane_with_support_friction(
        &scene,
        before,
        &law(),
        duration,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(report.state_after.mode(), SpherePlaneMode::Sliding);
    close(velocity(report.state_after.body()).x, 0.4 + 6.0 * duration);
    close(spin(report.state_after.body()).z, -5.0 * duration);
    close(
        report.dissipated_energy,
        2.0 * (0.4 * duration + 0.5 * 3.5 * duration * duration),
    );
    balances(&scene, before, duration, &report);
}

#[test]
fn g1_arrest_resticks_inside_static_band_and_reverses_above_it() {
    let duration = 0.4;
    let before = supported(1.4);
    for force in [-14.0, -24.0] {
        let scene = SpherePlane {
            force_world: Vec3::new(force, -20.0, 0.0),
            ..scene()
        };
        let report = advance_sphere_plane_with_support_friction(
            &scene,
            before,
            &law(),
            duration,
            SpherePlaneConfig::default(),
            || false,
        )
        .unwrap();
        let arrest = 1.4 / (3.5 - force / MASS);
        let remaining = duration - arrest;
        let v_arrest = 1.4 + (force - 2.0) / MASS * arrest;
        let x_arrest = 1.4 * arrest + 0.5 * (force - 2.0) / MASS * arrest * arrest;
        let holds = force == -14.0;
        let a_after = if holds {
            force / 2.8
        } else {
            (force + 2.0) / MASS
        };
        close(
            velocity(report.state_after.body()).x,
            v_arrest + a_after * remaining,
        );
        close(
            report.state_after.body().pose().position_world().x,
            x_arrest + v_arrest * remaining + 0.5 * a_after * remaining * remaining,
        );
        close(
            report
                .events
                .iter()
                .find(|e| e.kind == SpherePlaneEventKind::Rolling)
                .unwrap()
                .time_s,
            arrest,
        );
        assert_eq!(
            report.state_after.mode(),
            if holds {
                SpherePlaneMode::Rolling
            } else {
                SpherePlaneMode::Sliding
            }
        );
        let loss = 1.4 * arrest
            + if holds {
                0.0
            } else {
                (-force / MASS - 3.5) * remaining * remaining
            };
        close(report.dissipated_energy, loss);
        if !holds {
            close(
                report
                    .events
                    .iter()
                    .find(|e| e.kind == SpherePlaneEventKind::Sliding)
                    .unwrap()
                    .time_s,
                arrest,
            );
        }
        balances(&scene, before, duration, &report);
    }
    // Zero kinetic friction still permits force-driven slip arrest followed by
    // static recapture. Integrating the entire frictionless interval would miss it.
    let scene = SpherePlane {
        force_world: Vec3::new(-14.0, -20.0, 0.0),
        ..scene()
    };
    let report = advance_sphere_plane_with_support_friction(
        &scene,
        before,
        &FrictionLaw::Coulomb {
            static_mu: 0.3,
            kinetic_mu: 0.0,
        },
        duration,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(report.state_after.mode(), SpherePlaneMode::Rolling);
    close(velocity(report.state_after.body()).x, -1.0);
    close(spin(report.state_after.body()).z, 2.0);
    close(report.state_after.body().pose().position_world().x, 0.04);
    close(report.dissipated_energy, 0.0);
    balances(&scene, before, duration, &report);
}

#[test]
fn g1_support_law_does_not_replace_the_instantaneous_impact_law() {
    let scene = scene();
    let before = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(2.0, -1.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let original =
        advance_sphere_plane(&scene, before, 0.03, SpherePlaneConfig::default(), || false).unwrap();
    let report = advance_sphere_plane_with_support_friction(
        &scene,
        before,
        &FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.05,
        },
        0.03,
        SpherePlaneConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(report.state_after, original.state_after);
    assert_eq!(report.dissipated_energy, original.dissipated_energy);
    close(report.events[0].normal_impulse, 3.0);
    close(report.events[0].tangent_impulse_world.x, -0.6);
    close(velocity(report.state_after.body()).x, 1.7);
    close(spin(report.state_after.body()).z, -1.5);
    balances(&scene, before, 0.03, &report);
}

#[test]
fn g3_curved_slip_retains_kinetic_flow_without_coarse_static_capture() {
    // Transverse driving equals or exceeds the kinetic, but not static, limit.
    // At equality, the continuum invariant |u|+u_z is positive and the slip
    // approaches a nonzero limiting speed. Coarse disk capture is unphysical.
    let before = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(0.4, 0.0, -1.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    for force in [4.0, 8.0] {
        let scene = SpherePlane {
            force_world: Vec3::new(0.0, -20.0, force),
            friction_coefficient: 2.0 / 35.0,
            ..scene()
        };
        let support = FrictionLaw::Coulomb {
            static_mu: 0.3,
            kinetic_mu: scene.friction_coefficient,
        };
        for step in [0.4, 0.25, 0.025, 0.0125] {
            let config = SpherePlaneConfig {
                max_sliding_step_s: step,
                ..SpherePlaneConfig::default()
            };
            let original = advance_sphere_plane(&scene, before, 0.4, config, || false).unwrap();
            let report = advance_sphere_plane_with_support_friction(
                &scene,
                before,
                &support,
                0.4,
                config,
                || false,
            )
            .unwrap();
            assert_eq!(report.state_after, original.state_after);
            assert_eq!(report.state_after.mode(), SpherePlaneMode::Sliding);
            assert!(
                !report
                    .events
                    .iter()
                    .any(|e| e.kind == SpherePlaneEventKind::Rolling)
            );
            assert_eq!(report.dissipated_energy, original.dissipated_energy);
            balances(&scene, before, 0.4, &report);
        }
    }
}

#[test]
fn g3_g4_moving_support_continuation_and_refusals_preserve_physical_state() {
    let scene = SpherePlane {
        velocity_world: Vec3::new(0.7, 0.0, -0.2),
        force_world: Vec3::new(-14.0, -20.0, 0.0),
        ..scene()
    };
    let before = state(
        Vec3::new(0.0, RADIUS, 0.0),
        Vec3::new(2.1, 0.0, -0.2),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let config = SpherePlaneConfig::default();
    let one =
        advance_sphere_plane_with_support_friction(&scene, before, &law(), 0.4, config, || false)
            .unwrap();
    let first =
        advance_sphere_plane_with_support_friction(&scene, before, &law(), 0.1, config, || false)
            .unwrap();
    let second = advance_sphere_plane_with_support_friction(
        &scene,
        first.state_after,
        &law(),
        0.3,
        config,
        || false,
    )
    .unwrap();
    close_vector(
        one.state_after.body().pose().position_world(),
        second.state_after.body().pose().position_world(),
    );
    close_vector(
        velocity(one.state_after.body()),
        velocity(second.state_after.body()),
    );
    close_vector(
        spin(one.state_after.body()),
        spin(second.state_after.body()),
    );
    close(
        one.dissipated_energy,
        first.dissipated_energy + second.dissipated_energy,
    );
    close(one.plane_work, first.plane_work + second.plane_work);
    close(
        one.plane_work,
        one.contact_impulse_world.dot(scene.velocity_world),
    );
    balances(&scene, before, 0.4, &one);
    for invalid in [
        FrictionLaw::Coulomb {
            static_mu: 0.1,
            kinetic_mu: 0.2,
        },
        FrictionLaw::Coulomb {
            static_mu: f64::NAN,
            kinetic_mu: 0.1,
        },
        FrictionLaw::Coulomb {
            static_mu: 0.3,
            kinetic_mu: -0.1,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.3,
            kinetic_mu: 0.1,
            characteristic_speed: 0.1,
            viscous_per_speed: 0.0,
        },
    ] {
        assert!(
            advance_sphere_plane_with_support_friction(
                &scene,
                before,
                &invalid,
                0.4,
                config,
                || false
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        advance_sphere_plane_with_support_friction(&scene, before, &law(), 0.4, config, || {
            polls += 1;
            polls == 4
        }),
        Err(SpherePlaneError::Cancelled)
    ));
    assert!(matches!(
        advance_sphere_plane_with_support_friction(
            &scene,
            before,
            &law(),
            0.4,
            SpherePlaneConfig {
                max_events: 1,
                ..config
            },
            || false,
        ),
        Err(SpherePlaneError::EventBudget)
    ));
    let retry =
        advance_sphere_plane_with_support_friction(&scene, before, &law(), 0.4, config, || false)
            .unwrap();
    assert_eq!(retry.state_after, one.state_after);
}
