//! G1/G3: finite-range attraction, loading stiffness, instabilities and work.
#[path = "source_point/common.rs"]
mod common;

use fs_contact::normal_patch::dmt::{fixture::*, *};
use fs_tribo::InputAuthority;

fn contact(length_scale: f64) -> SourceDmtSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 0.01 * length_scale)
        .normal_state()
        .clone();
    SourceDmtSphereContact::new(
        normal,
        1e-6 * length_scale,
        1e-9 * length_scale,
        DmtEnvelope {
            maximum_indentation_m: 1e-9 * length_scale,
            half_space_depth_m: 1e-3 * length_scale,
            layer_thickness_m: 1e-3 * length_scale,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}
fn fixture(stiffness_ratio: f64) -> DmtSpringFixture {
    let c = contact(1.0);
    let attraction_slope = -c.force_controlled_pull_off_n() / c.interaction_range_m();
    DmtSpringFixture::new(c, stiffness_ratio * attraction_slope).unwrap()
}
fn close(a: f64, b: f64, scale: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * scale,
        "{a:e} != {b:e}, error/scale = {:e}",
        (a - b).abs() / scale
    );
}
fn energy_at(f: &DmtSpringFixture, stage: f64, delta: f64) -> f64 {
    f.contact().evaluate(delta).unwrap().reversible_energy_j
        + 0.5 * f.stiffness_n_per_m() * (stage - delta).powi(2)
}

#[test]
fn source_and_equilibrium_survive_finite_compliance() {
    for ratio in [0.1, 0.75, 1.25, 10.0, 1e6] {
        let f = fixture(ratio);
        let z = f.contact().interaction_range_m();
        let pull = -f.contact().force_controlled_pull_off_n();
        for d in [0.0, 1e-5, 0.02, 0.3, 1.0] {
            let expected = f.contact().evaluate(d * z).unwrap();
            let stage = expected.indentation_m + expected.normal_force_n / f.stiffness_n_per_m();
            let actual = f.evaluate_contact(stage).unwrap();
            close(
                actual.contact.indentation_m,
                expected.indentation_m,
                z,
                2e-13,
            );
            close(
                actual.contact.normal_force_n,
                expected.normal_force_n,
                pull.max(expected.normal_force_n.abs()),
                2e-13,
            );
            close(
                actual.equilibrium_residual_n,
                0.0,
                expected.normal_force_n.abs()
                    + f.stiffness_n_per_m() * (stage.abs() + actual.contact.indentation_m.abs()),
                // Near rigid control, D-delta subtracts nearly equal
                // coordinates. Its rounding scale is K(|D|+|delta|), not F.
                16.0 * f64::EPSILON,
            );
            assert_eq!(
                actual.contact.normal_state_identity(),
                f.contact().normal_state().identity()
            );
            assert_eq!(actual.contact.authority, InputAuthority::SyntheticFixture);
            close(
                actual.total_stored_energy_j,
                energy_at(&f, stage, actual.contact.indentation_m),
                actual.total_stored_energy_j.abs().max(pull * z),
                2e-12,
            );
        }
    }
}

#[test]
fn attraction_gradient_selects_stable_gap_or_snap_hysteresis() {
    for ratio in [1.1, 2.0, 10.0] {
        let f = fixture(ratio);
        assert!(!f.has_snap_hysteresis());
        let z = f.contact().interaction_range_m();
        let pull = -f.contact().force_controlled_pull_off_n();
        let stage = 0.5 * (-z + f.touch_stage_m());
        let gap = f.evaluate_noncontact(stage).unwrap();
        close(gap.contact.indentation_m, -0.5 * z, z, 1e-14);
        close(gap.contact.normal_force_n, -0.5 * pull, pull, 1e-14);
        assert_eq!(gap.contact.regime, DmtRegime::AttractiveGap);
        assert_eq!(gap.contact.contact_area_m2, 0.0);
        let shear = gap.contact.friction_contact_state().unwrap();
        assert_eq!(shear.area_m2(), 0.0);
        assert!(gap.stage_tangent_n_per_m.unwrap() < 0.0);
        assert!(f.stiffness_n_per_m() + gap.contact.tangent_n_per_m.unwrap() > 0.0);
    }
    let f = fixture(0.25);
    let z = f.contact().interaction_range_m();
    let s = -f.contact().force_controlled_pull_off_n() / z;
    assert!(f.has_snap_hysteresis());
    // The algebraic tail solution at D=-2z is a maximum of total potential.
    let stage = -2.0 * z;
    let unstable_delta = (f.stiffness_n_per_m() * stage + s * z) / (f.stiffness_n_per_m() - s);
    let h = z * 1e-3;
    let curvature = (energy_at(&f, stage, unstable_delta + h)
        + energy_at(&f, stage, unstable_delta - h)
        - 2.0 * energy_at(&f, stage, unstable_delta))
        / (h * h);
    close(curvature, f.stiffness_n_per_m() - s, s, 1e-8);
    assert!(curvature < 0.0);
    assert_eq!(
        f.evaluate_noncontact(stage).unwrap().contact.regime,
        DmtRegime::Separated
    );
    assert!(matches!(
        f.evaluate_noncontact(-0.9 * z),
        Err(DmtFixtureError::NoNoncontactEquilibrium)
    ));
}

#[test]
fn soft_fixture_localizes_both_jumps_and_counts_release_once() {
    let f = fixture(0.25);
    let z = f.contact().interaction_range_m();
    let pull = -f.contact().force_controlled_pull_off_n();
    close(f.touch_stage_m(), -4.0 * z, z, 1e-15);
    let initial = f.detached_state(-6.0 * z).unwrap();
    let approach = f.advance_monotone(initial, 0.0).unwrap();
    let attach = approach.jump.unwrap();
    assert_eq!(attach.kind, DmtFixtureJumpKind::SnapIn);
    assert_eq!(attach.stage_displacement_m, -z);
    assert_eq!(attach.before_indentation_m, -z);
    assert!(attach.after_indentation_m > 0.0);
    let event = f.evaluate_contact(-z).unwrap();
    close(
        attach.released_energy_j,
        -event.total_stored_energy_j,
        pull * z,
        2e-14,
    );
    assert_eq!(approach, f.advance_monotone(initial, 0.0).unwrap());
    close(
        approach.prescribed_work_j,
        approach.response.total_stored_energy_j + approach.irreversible_work_j,
        pull * z,
        3e-13,
    );
    let retreat = f
        .advance_monotone(approach.response.state(), -6.0 * z)
        .unwrap();
    let detach = retreat.jump.unwrap();
    assert_eq!(detach.kind, DmtFixtureJumpKind::SnapOff);
    assert_eq!(detach.stage_displacement_m, f.touch_stage_m());
    assert_eq!(detach.before_indentation_m, 0.0);
    assert_eq!(detach.after_indentation_m, f.touch_stage_m());
    close(detach.released_energy_j, 1.5 * pull * z, pull * z, 1e-15);
    assert_eq!(retreat.response.contact.regime, DmtRegime::Separated);
    assert_eq!(retreat.response.total_stored_energy_j, 0.0);
    let again = f
        .advance_monotone(retreat.response.state(), -7.0 * z)
        .unwrap();
    assert_eq!(again.jump, None);
    assert_eq!(again.irreversible_work_j, 0.0);
    assert_eq!(again.prescribed_work_j, 0.0);
    close(
        approach.prescribed_work_j + retreat.prescribed_work_j,
        attach.released_energy_j + detach.released_energy_j,
        pull * z,
        3e-13,
    );
}

fn integrate_contact(f: &DmtSpringFixture, from: f64, to: f64) -> f64 {
    let n = 4096;
    let step = (to - from) / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let stage = if i == n { to } else { from + i as f64 * step };
        let force = f.evaluate_contact(stage).unwrap().contact.normal_force_n;
        sum += force
            * if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
    }
    sum * step / 3.0
}

#[test]
fn cycle_work_matches_independent_force_integral_and_stiff_cycle_is_reversible() {
    for ratio in [0.1, 0.5, 0.9] {
        let f = fixture(ratio);
        let z = f.contact().interaction_range_m();
        let scale = -f.contact().force_controlled_pull_off_n() * z;
        let initial = f.detached_state(2.0 * f.touch_stage_m()).unwrap();
        let approach = f.advance_monotone(initial, 0.0).unwrap();
        let retreat = f
            .advance_monotone(approach.response.state(), initial.stage_displacement_m())
            .unwrap();
        let integral =
            integrate_contact(&f, -z, 0.0) + integrate_contact(&f, 0.0, f.touch_stage_m());
        close(
            integral,
            approach.irreversible_work_j + retreat.irreversible_work_j,
            scale,
            2e-8,
        );
        close(
            integral,
            approach.prescribed_work_j + retreat.prescribed_work_j,
            scale,
            2e-8,
        );
    }
    let f = fixture(2.0);
    let z = f.contact().interaction_range_m();
    let initial = f.detached_state(-3.0 * z).unwrap();
    let mut state = initial;
    let mut work = 0.0;
    for d in [-1.0, -0.75, -0.5, 0.0, -0.5, -0.75, -1.0, -3.0] {
        let step = f.advance_monotone(state, d * z).unwrap();
        assert!(step.jump.is_none());
        assert_eq!(step.irreversible_work_j, 0.0);
        work += step.prescribed_work_j;
        state = step.response.state();
    }
    close(
        work,
        0.0,
        -f.contact().force_controlled_pull_off_n() * z,
        3e-13,
    );
    assert_eq!(state, initial);
}

#[test]
fn potential_gradient_and_condensed_tangent_match_the_stage_force() {
    for ratio in [0.25, 2.0, 100.0] {
        let f = fixture(ratio);
        let z = f.contact().interaction_range_m();
        let pull = -f.contact().force_controlled_pull_off_n();
        for stage in [0.5 * f.touch_stage_m(), 0.0, 0.5 * z] {
            let h = 1e-5 * z;
            let a = f.evaluate_contact(stage - h).unwrap();
            let b = f.evaluate_contact(stage + h).unwrap();
            let c = f.evaluate_contact(stage).unwrap();
            close(
                (b.total_stored_energy_j - a.total_stored_energy_j) / (2.0 * h),
                c.contact.normal_force_n,
                pull.max(c.contact.normal_force_n.abs()),
                4e-8,
            );
            close(
                (b.contact.normal_force_n - a.contact.normal_force_n) / (2.0 * h),
                c.stage_tangent_n_per_m.unwrap(),
                f.stiffness_n_per_m(),
                1e-8,
            );
            let port = c.stage_port([0.0, 0.0, 1.0], [0.1, 0.0, 0.0]).unwrap();
            assert_eq!(port.tangent_n_per_m, c.stage_tangent_n_per_m.unwrap());
            for k in 0..3 {
                assert_eq!(port.action_force_n[k], -port.reaction_force_n[k]);
                assert_eq!(port.action_moment_n_m[k], -port.reaction_moment_n_m[k]);
            }
            assert_eq!(port.dissipated_power_w, 0.0);
        }
    }
    let f = fixture(2.0);
    let z = f.contact().interaction_range_m();
    let center = f.evaluate_noncontact(-0.75 * z).unwrap();
    let left = f.evaluate_noncontact(-0.75001 * z).unwrap();
    let right = f.evaluate_noncontact(-0.74999 * z).unwrap();
    close(
        (right.total_stored_energy_j - left.total_stored_energy_j) / (0.00002 * z),
        center.contact.normal_force_n,
        center.contact.normal_force_n.abs(),
        1e-9,
    );
    close(
        center.stage_tangent_n_per_m.unwrap(),
        -f.stiffness_n_per_m(),
        f.stiffness_n_per_m(),
        1e-14,
    );
}

#[test]
fn length_scaling_refusals_and_immutable_source_history() {
    let f = fixture(0.5);
    let scaled = DmtSpringFixture::new(contact(4.0), 4.0 * f.stiffness_n_per_m()).unwrap();
    let z = f.contact().interaction_range_m();
    let initial = f.detached_state(-4.0 * z).unwrap();
    let a = f.advance_monotone(initial, 0.0).unwrap();
    let b = scaled
        .advance_monotone(scaled.detached_state(-16.0 * z).unwrap(), 0.0)
        .unwrap();
    close(
        b.response.contact.indentation_m,
        4.0 * a.response.contact.indentation_m,
        z,
        1e-13,
    );
    close(
        b.response.contact.normal_force_n,
        16.0 * a.response.contact.normal_force_n,
        a.response.contact.normal_force_n.abs(),
        1e-13,
    );
    close(
        b.irreversible_work_j,
        64.0 * a.irreversible_work_j,
        a.irreversible_work_j,
        1e-12,
    );
    assert!(matches!(
        scaled.advance_monotone(initial, 0.0),
        Err(DmtFixtureError::StateMismatch)
    ));
    let c = contact(1.0);
    let critical = -c.force_controlled_pull_off_n() / c.interaction_range_m();
    assert!(matches!(
        DmtSpringFixture::new(c.clone(), critical),
        Err(DmtFixtureError::NeutralAttraction)
    ));
    assert!(DmtSpringFixture::new(c.clone(), -1.0).is_err());
    assert!(DmtSpringFixture::new(c, f64::INFINITY).is_err());
    assert!(f.evaluate_contact(f.touch_stage_m().next_down()).is_err());
    assert!(f.evaluate_contact(f.maximum_stage_m().next_up()).is_err());
    assert!(f.advance_monotone(initial, f64::NAN).is_err());
    assert!(f.detached_state(-0.5 * z).is_err());
    let corner = f.evaluate_contact(f.touch_stage_m()).unwrap();
    assert!(corner.stage_tangent_n_per_m.is_none());
    assert!(corner.stage_port([0.0, 0.0, 1.0], [0.0; 3]).is_err());
}
