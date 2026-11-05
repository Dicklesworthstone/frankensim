//! G0/G1/G3: original adhesive work, wrench projection and atomic receiver state.
#[path = "source_point/common.rs"]
mod common;
use core::f64::consts::PI;
use fs_contact::normal_patch::{
    NormalLawSource,
    adhesive::{JkrEnvelope, SourceJkrSphereContact},
    maugis::{MaugisEnvelope, SourceMaugisSphereContact},
    source_port::*,
};
use fs_exec::CancelGate;

fn jkr() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
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
fn maugis() -> SourceMaugisSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let radius = 0.01;
    let a0 =
        (PI * 1e-7 * radius * radius / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0)).cbrt();
    let d0 = a0 * a0 / radius;
    SourceMaugisSphereContact::new(
        normal,
        radius,
        2.0 * d0 / (PI * 0.2),
        MaugisEnvelope {
            maximum_indentation_m: 20.0 * d0,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}
fn geometry() -> SourceNormalPortGeometry {
    SourceNormalPortGeometry {
        normal: [0.0, 0.0, 1.0],
        moment_arm_m: [0.1, 0.2, 0.0],
    }
}
fn controls() -> SourceNormalPortControls {
    SourceNormalPortControls {
        maximum_step_s: 0.001,
        maximum_absolute_approach_rate_m_per_s: 1.0,
        maximum_force_n: 1e3,
        maximum_absolute_energy_j: 1.0,
        kinematic_absolute_tolerance_m: 1e-15,
        work_absolute_tolerance_j: 1e-20,
        work_relative_tolerance: 1e-9,
        maximum_steps: 100,
    }
}
fn owner(source: NormalLawSource, from: f64) -> SourceNormalPathPort {
    SourceNormalPathPort::new(
        source,
        geometry(),
        controls(),
        -from,
        0.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn step(from: f64, to: f64, time: f64) -> SourceNormalPortStep {
    SourceNormalPortStep {
        end_time_s: time,
        gap_m: -to,
        action: NormalPortBodyIncrement {
            translation_m: [0.0, 0.0, -(to - from)],
            rotation_rad: [0.0; 3],
        },
        reaction: NormalPortBodyIncrement::default(),
        converged: true,
    }
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 1e-9 * scale.max(1e-30), "{a:e} != {b:e}");
}

#[test]
fn positive_gap_and_tensile_force_keep_original_work_and_both_wrenches() {
    let source = jkr();
    let from = 0.25 * source.minimum_attached_indentation_m();
    let to = 0.1 * source.minimum_attached_indentation_m();
    let expected = source.attached_path_average(from, to).unwrap();
    let mut port = owner(source.clone().into(), from);
    let r = port
        .advance(step(from, to, 1e-4), &CancelGate::new_clock_free())
        .unwrap();
    assert!(r.from_approach_m < 0.0 && r.to_approach_m < 0.0 && r.mean_force_n < 0.0);
    assert_eq!(r.source_identity, source.identity());
    close(r.mean_force_n, expected, expected.abs());
    close(
        r.source_work_j,
        expected * (to - from),
        (expected * (to - from)).abs(),
    );
    close(r.work_on_bodies_j, -r.source_work_j, r.source_work_j.abs());
    assert!(r.work_on_bodies_j > 0.0);
    for i in 0..3 {
        assert_eq!(
            r.wrenches.action_force_n[i],
            -r.wrenches.reaction_force_n[i]
        );
        assert_eq!(
            r.wrenches.action_moment_n_m[i],
            -r.wrenches.reaction_moment_n_m[i]
        );
    }
    close(
        r.wrenches.action_moment_n_m[0],
        0.2 * expected,
        expected.abs(),
    );
    close(
        r.wrenches.action_moment_n_m[1],
        -0.1 * expected,
        expected.abs(),
    );
    assert_eq!(port.approach_m(), to);
    assert_eq!(port.accepted_steps(), 1);
}

#[test]
fn rotational_body_work_and_refinement_recover_the_same_conservative_storage() {
    let source = jkr();
    let mut full = owner(source.clone().into(), 0.0);
    let mut split = owner(source.into(), 0.0);
    let to = 1e-6;
    let r = full
        .advance(step(0.0, to, 2e-4), &CancelGate::new_clock_free())
        .unwrap();
    let mut total = 0.0;
    for (from, to, time) in [(0.0, 0.4e-6, 1e-4), (0.4e-6, 1e-6, 2e-4)] {
        let mut s = step(from, to, time);
        // Equal translation carries no relative work. The original moment
        // supplies ALL the normal work through a fixed-frame rotation.
        s.action.translation_m = [1.0, 2.0, 3.0];
        s.reaction.translation_m = s.action.translation_m;
        s.action.rotation_rad = [-(to - from) / geometry().moment_arm_m[1], 0.0, 0.0];
        total += split
            .advance(s, &CancelGate::new_clock_free())
            .unwrap()
            .work_on_bodies_j;
    }
    close(total, r.work_on_bodies_j, r.source_work_j.abs());
    close(
        full.stored_energy_j(),
        split.stored_energy_j(),
        full.stored_energy_j().abs(),
    );
    let back = split
        .advance(step(to, 0.0, 3e-4), &CancelGate::new_clock_free())
        .unwrap();
    close(total + back.work_on_bodies_j, 0.0, r.source_work_j.abs());
}

#[test]
fn unilateral_onset_retains_inactive_path_and_does_not_omit_hunt_crossley_loss() {
    let source = common::source(1e7, None);
    let from = -2e-5;
    let to = 1e-5;
    let active = source
        .conservative_path_average(0.0, to)
        .unwrap()
        .mean_force_n;
    let mut port = owner(source.into(), from);
    let r = port
        .advance(step(from, to, 2e-4), &CancelGate::new_clock_free())
        .unwrap();
    close(r.mean_force_n, active / 3.0, active);
    assert_eq!(r.stored_energy_before_j, 0.0);
    close(
        r.work_on_bodies_j,
        -r.stored_energy_after_j,
        r.stored_energy_after_j,
    );
    assert!(matches!(
        SourceNormalPathPort::new(
            common::source(1e7, Some(0.4)).into(),
            geometry(),
            controls(),
            0.0,
            0.0,
            &CancelGate::new_clock_free()
        ),
        Err(SourceNormalPortError::Invalid(_))
    ));
}

#[test]
fn admitted_maugis_noncontact_has_native_attraction_and_restricted_branch_refuses() {
    use fs_contact::normal_patch::maugis::{
        MaugisWorkBudget, continuation::MaugisContinuationBudget,
        mechanical::MaugisAttachedPotential,
    };
    let source = maugis();
    let limits = source.noncontact_limits().unwrap().unwrap();
    let width = limits.first_contact_indentation_m - limits.minimum_indentation_m;
    let from = limits.minimum_indentation_m + 0.4 * width;
    let to = limits.minimum_indentation_m + 0.6 * width;
    let potential = MaugisAttachedPotential::new(
        source.clone(),
        MaugisContinuationBudget {
            work: MaugisWorkBudget {
                absolute_tolerance_j: 1e-30,
                relative_tolerance: 1e-10,
                maximum_evaluations: 10000,
            },
            energy_absolute_tolerance_j: 1e-27,
            energy_relative_tolerance: 1e-9,
        },
    )
    .unwrap();
    assert!(
        SourceNormalPathPort::new(
            potential.clone().into(),
            geometry(),
            controls(),
            -from,
            0.0,
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut c = controls();
    c.work_absolute_tolerance_j = 1e-27;
    let mut port = SourceNormalPathPort::new(
        potential.with_noncontact_continuation().unwrap().into(),
        geometry(),
        c,
        -from,
        0.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let r = port
        .advance(step(from, to, 1e-4), &CancelGate::new_clock_free())
        .unwrap();
    assert!(r.mean_force_n < 0.0);
    assert_eq!(
        source.evaluate_noncontact(to).unwrap().contact_response(),
        None
    );
    close(r.work_on_bodies_j, -r.source_work_j, r.source_work_j.abs());
}

#[test]
fn body_motion_is_checked_before_callback_and_work_cannot_hide_in_large_storage() {
    let mut port = owner(jkr().into(), 0.0);
    let before = port.clone();
    let mut bad = step(0.0, 1e-6, 1e-4);
    bad.action = NormalPortBodyIncrement::default();
    let mut called = false;
    let rejected = port.advance_with_admission(bad, &CancelGate::new_clock_free(), |_| {
        called = true;
        Ok::<(), SourceNormalPortError>(())
    });
    assert!(matches!(
        rejected,
        Err(SourceNormalPortError::Kinematics { .. })
    ));
    assert!(!called);
    assert_eq!(port, before);
    // Even a permissive geometric tolerance cannot silently publish work
    // inconsistent with the actual receiver motion.
    let mut c = controls();
    c.kinematic_absolute_tolerance_m = 1e-3;
    let mut loose = SourceNormalPathPort::new(
        jkr().into(),
        geometry(),
        c,
        0.0,
        0.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let old = loose.clone();
    assert!(matches!(
        loose.advance(bad, &CancelGate::new_clock_free()),
        Err(SourceNormalPortError::Work {
            channel: "action/reaction projection",
            ..
        })
    ));
    assert_eq!(loose, old);
}

#[test]
fn late_receiver_cancellation_domain_and_replay_refusals_preserve_accepted_state() {
    let source = jkr();
    let mut port = owner(source.clone().into(), 0.0);
    let initial = port.clone();
    let s = step(0.0, 1e-6, 1e-4);
    let failed = port.advance_with_admission(s, &CancelGate::new_clock_free(), |_| {
        Err::<(), _>(SourceNormalPortError::Invalid("receiver refused"))
    });
    assert!(failed.is_err());
    assert_eq!(port, initial);
    let gate = CancelGate::new_clock_free();
    let cancelled = port.advance_with_admission(s, &gate, |_| {
        gate.request();
        Ok::<(), SourceNormalPortError>(())
    });
    assert_eq!(cancelled, Err(SourceNormalPortError::Cancelled));
    assert_eq!(port, initial);
    let outside = source.minimum_attached_indentation_m().next_down();
    assert!(
        port.advance(step(0.0, outside, 1e-4), &CancelGate::new_clock_free())
            .is_err()
    );
    assert_eq!(port, initial);
    port.advance(s, &CancelGate::new_clock_free()).unwrap();
    let accepted = port.clone();
    assert_eq!(
        port.advance(s, &CancelGate::new_clock_free()),
        Err(SourceNormalPortError::StaleTime)
    );
    assert_eq!(port, accepted);
}

#[test]
fn source_force_and_step_budgets_are_enforced_before_acceptance() {
    let source = jkr();
    let initial = source.evaluate_attached(0.0).unwrap();
    for restricted in [
        SourceNormalPortControls {
            maximum_force_n: 1.1 * initial.normal_force_n.abs(),
            ..controls()
        },
        SourceNormalPortControls {
            maximum_absolute_energy_j: (1.1 * initial.reversible_energy_j.abs()).max(1e-30),
            ..controls()
        },
        SourceNormalPortControls {
            maximum_absolute_approach_rate_m_per_s: 0.01,
            ..controls()
        },
    ] {
        let mut limited = SourceNormalPathPort::new(
            source.clone().into(),
            geometry(),
            restricted,
            0.0,
            0.0,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let old = limited.clone();
        assert!(
            limited
                .advance(step(0.0, 8e-6, 1e-4), &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(limited, old);
    }
    // The source's own stricter rate and duration bounds remain effective.
    for (end, time) in [(1e-5, 1e-5), (1e-5, 0.002)] {
        let mut c = controls();
        c.maximum_step_s = 0.01;
        c.maximum_absolute_approach_rate_m_per_s = 10.0;
        let mut limited = SourceNormalPathPort::new(
            common::source(1e7, None).into(),
            geometry(),
            c,
            0.0,
            0.0,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let old = limited.clone();
        assert!(matches!(
            limited.advance(step(0.0, end, time), &CancelGate::new_clock_free()),
            Err(SourceNormalPortError::Source(_))
        ));
        assert_eq!(limited, old);
    }
    let mut c = controls();
    c.maximum_steps = 1;
    let mut port = SourceNormalPathPort::new(
        jkr().into(),
        geometry(),
        c,
        0.0,
        0.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    port.advance(step(0.0, 1e-6, 1e-4), &CancelGate::new_clock_free())
        .unwrap();
    let accepted = port.clone();
    assert_eq!(
        port.advance(step(1e-6, 1e-6, 2e-4), &CancelGate::new_clock_free()),
        Err(SourceNormalPortError::StepBudget)
    );
    assert_eq!(port, accepted);
    let mut zero = geometry();
    zero.normal = [0.0; 3];
    assert!(
        SourceNormalPathPort::new(
            jkr().into(),
            zero,
            controls(),
            0.0,
            0.0,
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}
