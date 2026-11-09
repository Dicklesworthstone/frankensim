//! G0/G1/G3: original JKR area, tensile rolling memory and independent work.
#[path = "source_point/common.rs"]
mod common;
#[path = "adhesive_rolling/events.rs"]
mod events;

#[path = "adhesive_rolling/temperature.rs"]
mod temperature;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::rolling::*;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_tribo::InputAuthority;

fn normal(label: &str) -> SourceJkrSphereContact {
    let a = common::solid(label, 1000.0, 1e7);
    let b = common::solid("rolling-B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}
fn controls() -> JkrRollingControls {
    JkrRollingControls {
        critical_rolling_displacement_m: 1e-8,
        maximum_critical_displacement_to_contact_radius: 0.01,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn original_jkr_zero_and_tensile_loads_have_nonzero_rolling_stiffness() {
    let normal = normal("rolling-A");
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let a0 = (4.5 * PI * w * r * r / e).cbrt();
    let fc = 1.5 * PI * w * r;
    let mut stiffness = Vec::new();
    for load in [0.0, -0.5 * fc, -fc] {
        let response = normal.evaluate_force_controlled(load).unwrap();
        let source =
            SourceJkrRollingContact::new(normal.clone(), response.indentation_m, controls())
                .unwrap();
        let factor = (response.contact_radius_m / a0).powf(1.5);
        let k = 4.0 * fc / r * factor;
        stiffness.push(source.law().stiffness_n_per_m());
        close(source.law().stiffness_n_per_m(), k);
        assert!(
            (source.normal_response().normal_force_n - response.normal_force_n).abs() < 1e-12 * fc
        );
        assert_eq!(
            source.normal_response().normal_state_identity(),
            normal.normal_state().identity()
        );
        assert_eq!(
            source.normal_response().authority,
            InputAuthority::SyntheticFixture
        );
        let angle = 0.25 * controls().critical_rolling_displacement_m / r;
        let step = source.advance(&source.initial_state(), angle).unwrap();
        close(step.rolling().endpoint_moment_n_m, -k * r * r * angle);
        close(
            step.rolling().storage_after_j,
            0.5 * k * (r * angle).powi(2),
        );
        close(
            step.rolling().applied_work_j,
            step.rolling().storage_after_j,
        );
        assert_eq!(step.rolling().irreversible_loss_j, 0.0);
        assert!(step.rolling().endpoint_moment_n_m < 0.0);
    }
    close(stiffness[2] / stiffness[0], 0.5);
}

#[test]
fn source_owned_yield_reversal_and_partition_preserve_distinct_storage_and_loss() {
    let normal = normal("rolling-A");
    let response = normal.evaluate_force_controlled(0.0).unwrap();
    let source = SourceJkrRollingContact::new(normal, response.indentation_m, controls()).unwrap();
    let theta = source.law().critical_rolling_angle_rad();
    let initial = source.initial_state();
    let full = source.advance(&initial, 3.0 * theta).unwrap();
    let k = source.law().stiffness_n_per_m();
    let xi = controls().critical_rolling_displacement_m;
    close(full.rolling().irreversible_loss_j, 2.0 * k * xi * xi);
    close(full.rolling().storage_after_j, 0.5 * k * xi * xi);
    close(
        -full.rolling().mean_moment_n_m * 3.0 * theta,
        full.rolling().applied_work_j,
    );
    let mut state = initial;
    let (mut work, mut loss) = (0.0, 0.0);
    for increment in [0.25, 0.75, 0.5, 1.5] {
        let step = source.advance(&state, increment * theta).unwrap();
        work += step.rolling().applied_work_j;
        loss += step.rolling().irreversible_loss_j;
        state = step.state();
    }
    close(state.elastic_shift_m(), full.state().elastic_shift_m());
    close(work, full.rolling().applied_work_j);
    close(loss, full.rolling().irreversible_loss_j);
    let unload = source.advance(&state, -theta).unwrap();
    assert_eq!(unload.rolling().irreversible_loss_j, 0.0);
    assert!(unload.rolling().applied_work_j < 0.0);
    close(
        source.stored_energy_j(&state).unwrap(),
        full.rolling().storage_after_j,
    );
    let reverse = source.advance(&state, -3.0 * theta).unwrap();
    close(reverse.state().elastic_shift_m(), -xi);
    close(reverse.rolling().irreversible_loss_j, k * xi * xi);
}

#[test]
fn changed_normal_state_or_source_cannot_discard_or_rebind_rolling_memory() {
    let original = normal("rolling-A");
    let indentation = original
        .evaluate_force_controlled(0.0)
        .unwrap()
        .indentation_m;
    let source = SourceJkrRollingContact::new(original.clone(), indentation, controls()).unwrap();
    let loaded = source
        .advance(
            &source.initial_state(),
            0.5 * source.law().critical_rolling_angle_rad(),
        )
        .unwrap()
        .state();
    let saved = loaded;
    let changed = [
        SourceJkrRollingContact::new(original.clone(), indentation * 1.1, controls()).unwrap(),
        SourceJkrRollingContact::new(normal("other-source-A"), indentation, controls()).unwrap(),
        SourceJkrRollingContact::new(
            original,
            indentation,
            JkrRollingControls {
                critical_rolling_displacement_m: 2e-8,
                ..controls()
            },
        )
        .unwrap(),
    ];
    for other in changed {
        assert!(matches!(
            other.advance(&loaded, 0.0),
            Err(SourceJkrRollingError::StateMismatch)
        ));
        assert!(other.stored_energy_j(&loaded).is_err());
    }
    assert!(source.advance(&loaded, f64::NAN).is_err());
    assert_eq!(loaded, saved);
    let no_op = source.advance(&loaded, 0.0).unwrap();
    assert_eq!(no_op.state(), loaded);
    assert_eq!(no_op.rolling().applied_work_j, 0.0);
    assert_eq!(no_op.rolling().irreversible_loss_j, 0.0);
}

#[test]
fn prescribed_normal_loading_retains_rolling_memory_and_accounts_for_both_stores() {
    let normal = normal("rolling-A");
    let initial_response = normal.evaluate_force_controlled(0.0).unwrap();
    let source =
        SourceJkrRollingContact::new(normal.clone(), initial_response.indentation_m, controls())
            .unwrap();
    let initial = source
        .advance(
            &source.initial_state(),
            2.5 * source.law().critical_rolling_angle_rad(),
        )
        .unwrap()
        .state();
    let target = 0.8 * normal.envelope().maximum_indentation_m;
    let step = source.advance_normal(&initial, target).unwrap();
    let next = step.source();
    let end = next.normal_response();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let r = normal.effective_radius_m();
    let fc = 1.5 * PI * w * r;
    let a0 = (4.5 * PI * w * r * r / e).cbrt();
    let primary_stiffness = 4.0 * fc / r * (end.contact_radius_m / a0).powf(1.5);
    close(next.law().stiffness_n_per_m(), primary_stiffness);
    close(
        step.rolling_parameter_work_j(),
        0.5 * (primary_stiffness - source.law().stiffness_n_per_m())
            * initial.elastic_shift_m().powi(2),
    );
    assert!(step.rolling_parameter_work_j() > 0.0);
    assert_eq!(step.state().elastic_shift_m(), initial.elastic_shift_m());
    assert_eq!(
        step.state().accumulated_angle_rad(),
        initial.accumulated_angle_rad()
    );
    assert_eq!(next.normal_source(), source.normal_source());
    assert_eq!(next.controls(), source.controls());
    assert_eq!(end.authority, InputAuthority::SyntheticFixture);

    // Independent full normal potential at a prescribed radius and indentation.
    let potential = |a: f64, delta: f64| {
        let hertz = a * a / r;
        e * a * (delta - hertz / 3.0).powi(2) + (4.0 / 45.0) * e * a * hertz * hertz
            - PI * w * a * a
    };
    let normal_change = potential(end.contact_radius_m, end.indentation_m)
        - potential(
            initial_response.contact_radius_m,
            initial_response.indentation_m,
        );
    close(step.normal_contact_work_j(), normal_change);
    close(
        step.applied_work_j(),
        normal_change + step.rolling_storage_after_j() - step.rolling_storage_before_j(),
    );
    assert!(matches!(
        next.advance(&initial, 0.0),
        Err(SourceJkrRollingError::StateMismatch)
    ));

    // The complete prior plastic offset survives; unloading does not manufacture
    // a fresh stress-free contact simply because the normal load changed.
    let unload = next
        .advance(&step.state(), -next.law().critical_rolling_angle_rad())
        .unwrap();
    close(unload.state().elastic_shift_m(), 0.0);
    assert_eq!(unload.rolling().irreversible_loss_j, 0.0);
    close(
        unload.rolling().applied_work_j,
        -step.rolling_storage_after_j(),
    );
    let return_step = next
        .advance_normal(&step.state(), initial_response.indentation_m)
        .unwrap();
    assert_eq!(return_step.source(), &source);
    assert_eq!(return_step.state(), initial);
    close(
        return_step.normal_contact_work_j(),
        -step.normal_contact_work_j(),
    );
    assert_eq!(
        return_step.rolling_parameter_work_j(),
        -step.rolling_parameter_work_j()
    );
    close(return_step.applied_work_j(), -step.applied_work_j());
}

#[test]
fn attached_tensile_area_continuation_preserves_work_under_subdivision() {
    let normal = normal("rolling-A");
    let delta0 = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let source = SourceJkrRollingContact::new(normal.clone(), delta0, controls()).unwrap();
    let state = source
        .advance(
            &source.initial_state(),
            -0.5 * source.law().critical_rolling_angle_rad(),
        )
        .unwrap()
        .state();
    let delta1 = normal
        .evaluate_force_controlled(0.5 * normal.force_controlled_pull_off_n())
        .unwrap()
        .indentation_m;
    let delta2 = normal
        .evaluate_force_controlled(normal.force_controlled_pull_off_n())
        .unwrap()
        .indentation_m;
    let direct = source.advance_normal(&state, delta2).unwrap();
    let first = source.advance_normal(&state, delta1).unwrap();
    let second = first
        .source()
        .advance_normal(&first.state(), delta2)
        .unwrap();
    assert_eq!(second.source(), direct.source());
    assert_eq!(second.state(), direct.state());
    close(
        first.normal_contact_work_j() + second.normal_contact_work_j(),
        direct.normal_contact_work_j(),
    );
    close(
        first.rolling_parameter_work_j() + second.rolling_parameter_work_j(),
        direct.rolling_parameter_work_j(),
    );
    close(
        first.applied_work_j() + second.applied_work_j(),
        direct.applied_work_j(),
    );
    assert!(direct.source().normal_response().normal_force_n < 0.0);
    close(
        direct.source().law().stiffness_n_per_m(),
        0.5 * source.law().stiffness_n_per_m(),
    );
    assert!(direct.rolling_parameter_work_j() < 0.0);
    let continued = direct.source().advance(&direct.state(), 0.0).unwrap();
    assert!(continued.rolling().endpoint_moment_n_m > 0.0);
    assert_eq!(continued.rolling().irreversible_loss_j, 0.0);
}

#[test]
fn normal_continuation_refuses_detachment_and_inapplicable_area_atomically() {
    let normal_source = normal("rolling-A");
    let delta = 0.8 * normal_source.envelope().maximum_indentation_m;
    let source = SourceJkrRollingContact::new(normal_source.clone(), delta, controls()).unwrap();
    let state = source
        .advance(
            &source.initial_state(),
            source.law().critical_rolling_angle_rad(),
        )
        .unwrap()
        .state();
    let saved = source.clone();
    for target in [
        f64::NAN,
        f64::INFINITY,
        1.01 * normal_source.envelope().maximum_indentation_m,
        1.01 * normal_source.minimum_attached_indentation_m(),
    ] {
        assert!(source.advance_normal(&state, target).is_err());
    }
    let wrong = SourceJkrRollingContact::new(normal("other-rolling-A"), delta, controls()).unwrap();
    assert!(matches!(
        wrong.advance_normal(&state, delta),
        Err(SourceJkrRollingError::StateMismatch)
    ));
    let fold = normal_source.minimum_attached_indentation_m();
    let radius = normal_source
        .evaluate_attached(fold)
        .unwrap()
        .contact_radius_m;
    let limited = SourceJkrRollingContact::new(
        normal_source,
        delta,
        JkrRollingControls {
            critical_rolling_displacement_m: 0.015 * radius,
            maximum_critical_displacement_to_contact_radius: 0.01,
        },
    )
    .unwrap();
    assert!(matches!(
        limited.advance_normal(&limited.initial_state(), fold),
        Err(SourceJkrRollingError::Rolling(
            fs_tribo::rolling_loss::adhesive::AdhesiveRollingError::OutsideSmallShift { .. }
        ))
    ));
    assert_eq!(source, saved);
    let no_op = source.advance_normal(&state, delta).unwrap();
    assert_eq!(no_op.source(), &source);
    assert_eq!(no_op.state(), state);
    assert_eq!(no_op.normal_contact_work_j(), 0.0);
    assert_eq!(no_op.rolling_parameter_work_j(), 0.0);
    assert_eq!(no_op.applied_work_j(), 0.0);
}
