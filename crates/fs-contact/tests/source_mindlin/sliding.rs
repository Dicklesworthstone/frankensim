//! The same original bulk, interface and Hertz owners drive distributed spin/slip.

use core::f64::consts::PI;

use fs_contact::interface_binding::mindlin::sliding::{
    SourceHertzSlidingInput, SourceHertzSlidingPatch,
};

use super::*;

fn sliding(source: SourceMindlinContact) -> SourceHertzSlidingPatch {
    let checkpoint = source.initial_state();
    SourceHertzSlidingPatch::from_unloaded_source(source, &checkpoint).unwrap()
}

fn motion(source: &SourceHertzSlidingPatch) -> SourceHertzSlidingInput {
    SourceHertzSlidingInput {
        normal_force_n: source.normal_response().normal_force_n,
        slip_velocity_m_s: [0.03, -0.04],
        spin_rad_s: 12.0,
        dt_s: 1e-4,
        contact_temperature_k: Some(source.source().temperature_kelvin()),
        nominal_pressure_pa: None,
    }
}

#[test]
fn sourced_hertz_pressure_drives_pure_and_mixed_sliding_spin_with_original_receipts() {
    let original = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let source = sliding(original.clone());
    assert_eq!(source.source(), &original);
    assert_eq!(source.normal_source(), original.normal_source());
    assert_eq!(source.normal_response(), original.normal_response());
    assert_eq!(source.friction_source(), original.friction_source());
    let p = original.normal_response();
    let capacity = 0.6 * p.normal_force_n;
    let first_moment_radius = 3.0 * PI * p.patch_radius_m / 16.0;
    let input = motion(&source);
    let translate = source
        .evaluate(SourceHertzSlidingInput {
            spin_rad_s: 0.0,
            ..input
        })
        .unwrap();
    common::close(translate.wrench.force_n[0], -0.6 * capacity);
    common::close(translate.wrench.force_n[1], 0.8 * capacity);
    assert_eq!(translate.wrench.spin_moment_nm, 0.0);
    common::close(translate.dissipated_j, capacity * 0.05 * input.dt_s);
    let rotate = source
        .evaluate(SourceHertzSlidingInput {
            slip_velocity_m_s: [0.0; 2],
            ..input
        })
        .unwrap();
    assert_eq!(rotate.wrench.force_n, [0.0; 2]);
    common::close(
        rotate.wrench.spin_moment_nm,
        -capacity * first_moment_radius,
    );
    common::close(
        rotate.dissipated_j,
        capacity * first_moment_radius * input.spin_rad_s * input.dt_s,
    );
    let edge_speed = p.patch_radius_m * input.spin_rad_s;
    let mixed = source
        .evaluate(SourceHertzSlidingInput {
            slip_velocity_m_s: [0.6 * edge_speed, -0.8 * edge_speed],
            ..input
        })
        .unwrap();
    common::close(mixed.wrench.force_n[0], -0.6 * capacity * 9.0 * PI / 32.0);
    common::close(mixed.wrench.force_n[1], 0.8 * capacity * 9.0 * PI / 32.0);
    common::close(
        mixed.wrench.spin_moment_nm,
        -capacity * p.patch_radius_m * 9.0 * PI / 128.0,
    );
    common::close(
        mixed.dissipated_j,
        capacity * edge_speed * (45.0 * PI / 128.0) * input.dt_s,
    );
    assert_eq!(
        mixed.normal_model_identity,
        original.normal_source().model().identity()
    );
    assert_eq!(
        mixed.friction_state_identity,
        original.friction_source().identity()
    );
    assert_eq!(
        source.normal_response().authority,
        InputAuthority::SyntheticFixture
    );
}

#[test]
fn bulk_interface_and_normal_geometry_substitutions_change_actual_sliding_wrenches() {
    let soft = sliding(resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap());
    let stiff = sliding(resolve(normal(3e7, Some(0.6), Dims::NONE, false)).unwrap());
    let slippery = sliding(resolve(normal(1e7, Some(0.2), Dims::NONE, false)).unwrap());
    let original = normal(1e7, Some(0.6), Dims::NONE, false);
    let wider = sliding(
        resolve(
            SourcePointContact::new(
                original.model().clone(),
                PointContactGeometry::SpherePlane {
                    effective_radius_m: 0.4,
                },
                original.envelope(),
            )
            .unwrap(),
        )
        .unwrap(),
    );
    let pure_spin = |source: &SourceHertzSlidingPatch| {
        source
            .evaluate(SourceHertzSlidingInput {
                slip_velocity_m_s: [0.0; 2],
                ..motion(source)
            })
            .unwrap()
    };
    let a = pure_spin(&soft);
    let b = pure_spin(&stiff);
    let c = pure_spin(&slippery);
    let d = pure_spin(&wider);
    assert!(b.wrench.spin_moment_nm.abs() > a.wrench.spin_moment_nm.abs());
    common::close(c.wrench.spin_moment_nm, a.wrench.spin_moment_nm / 3.0);
    common::close(
        wider.normal_response().patch_radius_m,
        2.0 * soft.normal_response().patch_radius_m,
    );
    common::close(d.wrench.spin_moment_nm, 4.0 * a.wrench.spin_moment_nm);
    assert_ne!(a.normal_model_identity, b.normal_model_identity);
    assert_ne!(a.friction_state_identity, c.friction_state_identity);
}

#[test]
fn loaded_or_foreign_elastic_history_cannot_be_discarded_by_selecting_gross_sliding() {
    let source = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let initial = source.initial_state();
    let loaded = source
        .advance(
            &initial,
            &step(&source, 0.2 * source.law().gross_slip_displacement_m()),
        )
        .unwrap();
    let saved = loaded.state.clone();
    assert!(source.stored_energy_j(&saved) > 0.0);
    assert!(SourceHertzSlidingPatch::from_unloaded_source(source.clone(), &loaded.state).is_err());
    assert_eq!(loaded.state, saved);
    let other = resolve(normal_with_history(
        1e7,
        Some(0.6),
        Dims::NONE,
        false,
        "different specimen history",
    ))
    .unwrap();
    assert_eq!(other.law().parameters(), source.law().parameters());
    assert!(SourceHertzSlidingPatch::from_unloaded_source(other, &initial).is_err());
    let accepted = SourceHertzSlidingPatch::from_unloaded_source(source.clone(), &initial).unwrap();
    assert_eq!(initial, source.initial_state());
    let input = motion(&accepted);
    assert_eq!(
        accepted.evaluate(input).unwrap(),
        accepted.evaluate(input).unwrap()
    );
}

#[test]
fn fixed_source_and_original_normal_envelope_refuse_without_changing_the_binding() {
    let source = sliding(resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap());
    let saved = source.clone();
    let input = motion(&source);
    let expected = source.evaluate(input).unwrap();
    for bad in [
        SourceHertzSlidingInput {
            normal_force_n: input.normal_force_n * 1.01,
            ..input
        },
        SourceHertzSlidingInput {
            contact_temperature_k: None,
            ..input
        },
        SourceHertzSlidingInput {
            contact_temperature_k: Some(300.0),
            ..input
        },
        SourceHertzSlidingInput {
            nominal_pressure_pa: Some(1e5),
            ..input
        },
        SourceHertzSlidingInput {
            dt_s: 0.002,
            ..input
        },
        SourceHertzSlidingInput { dt_s: 0.0, ..input },
        SourceHertzSlidingInput {
            slip_velocity_m_s: [f64::NAN, 0.0],
            ..input
        },
        SourceHertzSlidingInput {
            spin_rad_s: f64::INFINITY,
            ..input
        },
    ] {
        assert!(source.evaluate(bad).is_err());
        assert_eq!(source, saved);
    }
    assert_eq!(source.evaluate(input).unwrap(), expected);
}

#[test]
fn explicitly_sourced_frictionless_patch_and_rest_do_not_invent_stiction_or_heat() {
    let source = sliding(resolve(normal(1e7, Some(0.0), Dims::NONE, false)).unwrap());
    let moving = source.evaluate(motion(&source)).unwrap();
    assert!(!moving.wrench.at_rest);
    assert_eq!(moving.wrench.force_n, [0.0; 2]);
    assert_eq!(moving.wrench.spin_moment_nm, 0.0);
    assert_eq!(moving.dissipated_j, 0.0);
    let friction = sliding(resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap());
    let rest = friction
        .evaluate(SourceHertzSlidingInput {
            slip_velocity_m_s: [0.0; 2],
            spin_rad_s: 0.0,
            ..motion(&friction)
        })
        .unwrap();
    assert!(rest.wrench.at_rest);
    assert_eq!(rest.wrench.force_n, [0.0; 2]);
    assert_eq!(rest.wrench.spin_moment_nm, 0.0);
    assert_eq!(rest.dissipated_j, 0.0);
}
