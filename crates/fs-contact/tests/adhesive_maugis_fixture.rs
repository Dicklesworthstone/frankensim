//! Original sourced Maugis contact under real loading-machine compliance.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::continuation::{
    MaugisContinuationBudget, MaugisDisplacementContinuation, MaugisPhase, MaugisSnapKind,
};
use fs_contact::normal_patch::maugis::fixture::{MaugisFixtureError, MaugisSpringFixture};
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, MaugisError, MaugisWorkBudget, SourceMaugisSphereContact,
};
use fs_tribo::InputAuthority;

fn model(lambda: f64) -> SourceMaugisSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let r: f64 = 0.01;
    let a0 = (PI * 1e-7 * r * r / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0)).cbrt();
    let d0 = a0 * a0 / r;
    SourceMaugisSphereContact::new(
        normal,
        r,
        2.0 * d0 / (PI * lambda),
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
fn scales(law: &SourceMaugisSphereContact) -> (f64, f64) {
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let r = law.effective_radius_m();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let a0 = (PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    (a0 * a0 / r, PI * w * r)
}
fn fixture(lambda: f64, stiffness: f64) -> MaugisSpringFixture {
    let law = model(lambda);
    let (d0, p0) = scales(&law);
    MaugisSpringFixture::new(law, stiffness * p0 / d0)
        .unwrap_or_else(|error| panic!("lambda={lambda}, stiffness={stiffness}: {error:?}"))
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scale={scale:e}"
    );
}
fn budget(fixture: &MaugisSpringFixture) -> MaugisContinuationBudget {
    let (d0, p0) = scales(fixture.contact());
    MaugisContinuationBudget {
        work: MaugisWorkBudget {
            absolute_tolerance_j: 1e-10 * p0 * d0,
            relative_tolerance: 1e-10,
            maximum_evaluations: 20000,
        },
        energy_absolute_tolerance_j: 1e-8 * p0 * d0,
        energy_relative_tolerance: 1e-8,
    }
}

#[test]
fn original_contact_equations_spring_energy_and_stage_derivatives_agree() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        for stiffness in [0.1, 1.0, 10.0, 1000.0] {
            let f = fixture(lambda, stiffness);
            let law = f.contact();
            let (d0, p0) = scales(law);
            // Above the force-controlled fold, hence stable for every K here.
            for d in [1.0, 3.0, 10.0] {
                let expected = law.evaluate_contact_energy(d * d0).unwrap();
                let stage =
                    expected.indentation_m() + expected.normal_force_n() / f.stiffness_n_per_m();
                let result = f.evaluate_contact(stage).unwrap();
                let actual = result.equilibrium();
                close(actual.indentation_m(), expected.indentation_m(), d0, 2e-12);
                close(
                    actual.normal_force_n(),
                    expected.normal_force_n(),
                    p0,
                    2e-11,
                );
                close(result.equilibrium_residual_n(), 0.0, p0, 2e-11);
                let spring = 0.5 * f.stiffness_n_per_m() * (stage - actual.indentation_m()).powi(2);
                close(result.spring_energy_j(), spring, p0 * d0, 2e-10);
                close(
                    result.total_potential_j(),
                    actual.elastic_energy_j() + actual.surface_energy_j() + spring,
                    p0 * d0,
                    2e-10,
                );
                assert_eq!(actual.model_identity(), law.identity());
                assert_eq!(actual.authority(), InputAuthority::SyntheticFixture);
                let h = 1e-5 * d0;
                let left = f.evaluate_contact(stage - h).unwrap();
                let right = f.evaluate_contact(stage + h).unwrap();
                close(
                    (right.total_potential_j() - left.total_potential_j()) / (2.0 * h),
                    actual.normal_force_n(),
                    p0,
                    3e-6,
                );
                close(
                    (right.equilibrium().normal_force_n() - left.equilibrium().normal_force_n())
                        / (2.0 * h),
                    result.stage_tangent_n_per_m().unwrap(),
                    p0 / d0,
                    2e-6,
                );
                let port = result
                    .stage_port([0.0, 0.0, 1.0], [0.01, 0.0, 0.0])
                    .unwrap();
                assert_eq!(port.action_force_n[2], actual.normal_force_n());
                assert_eq!(port.reaction_force_n[2], -actual.normal_force_n());
                assert_eq!(
                    port.tangent_n_per_m,
                    result.stage_tangent_n_per_m().unwrap()
                );
            }
        }
    }
}

#[test]
fn finite_loading_stiffness_moves_the_fold_between_control_limits() {
    let law = model(5.0);
    let (d0, p0) = scales(&law);
    let displacement = MaugisDisplacementContinuation::new(law.clone()).unwrap();
    let force = law
        .evaluate_contact_load(law.force_controlled_pull_off_n())
        .unwrap();
    let mut previous = force.indentation_m;
    for stiffness in [0.05, 0.5, 5.0, 50.0, 500.0] {
        let f = MaugisSpringFixture::new(law.clone(), stiffness * p0 / d0).unwrap();
        let fold = f.separation().before();
        let actual = fold.equilibrium();
        assert!(actual.indentation_m() < previous);
        previous = actual.indentation_m();
        assert!(actual.indentation_m() > law.minimum_attached_indentation_m());
        close(
            actual.tangent_n_per_m().unwrap(),
            -f.stiffness_n_per_m(),
            f.stiffness_n_per_m(),
            2e-9,
        );
        assert!(fold.stage_tangent_n_per_m().is_none());
        assert!(fold.stage_port([0.0, 0.0, 1.0], [0.0; 3]).is_err());
        assert!(
            f.evaluate_contact(f.minimum_contact_stage_m() - 1e-7 * d0)
                .is_err()
        );
    }
    let stiff = MaugisSpringFixture::new(law.clone(), 1e8 * p0 / d0).unwrap();
    close(
        stiff.separation().stage_displacement_m(),
        displacement.separation().indentation_m(),
        d0,
        2e-8,
    );
    let soft = MaugisSpringFixture::new(law.clone(), 1e-6 * p0 / d0).unwrap();
    close(
        soft.separation().before().equilibrium().normal_force_n(),
        law.force_controlled_pull_off_n(),
        p0,
        1e-9,
    );
}

#[test]
fn cancellation_conditioned_snap_landing_retains_force_and_refuses_unresolved_stiffness() {
    let f = fixture(1000.0, 1000.0);
    let (d0, p0) = scales(f.contact());
    let landing = f.onset().after();
    let contact = landing.equilibrium();
    let a = contact.contact_response().unwrap().contact_radius_m;
    let hertz_indentation = a * a / f.contact().effective_radius_m();
    // The original indentation subtracts the adhesive deflection from this
    // much larger Hertz term. K*|delta| alone misses that conditioning.
    assert!(hertz_indentation > 1000.0 * contact.indentation_m().abs());
    assert!(landing.equilibrium_residual_n().abs() <= landing.equilibrium_allowance_n());
    assert!(landing.equilibrium_allowance_n() < 1e-8 * p0);
    close(
        landing.stage_displacement_m(),
        contact.indentation_m() + contact.normal_force_n() / f.stiffness_n_per_m(),
        d0,
        4e-15,
    );
    assert_eq!(
        landing.equilibrium_residual_n(),
        contact.normal_force_n()
            - f.stiffness_n_per_m() * (landing.stage_displacement_m() - contact.indentation_m())
    );
    // At this stiffness the operand-derived allowance reaches the adhesion
    // load scale. It must refuse rather than admit arbitrary force mismatch.
    assert!(matches!(
        MaugisSpringFixture::new(f.contact().clone(), 1e13 * p0 / d0),
        Err(MaugisFixtureError::Contact(MaugisError::Unrepresentable {
            field: "resolved loading-spring force equilibrium"
        }))
    ));
}

#[test]
fn finite_spring_noncontact_disk_uses_deformed_geometry_and_no_shear_area() {
    for (lambda, stiffness) in [(0.02, 1.0), (0.2, 10.0), (0.8, 1000.0)] {
        let f = fixture(lambda, stiffness);
        let law = f.contact();
        let (d0, p0) = scales(law);
        let [minimum, touch] = f.noncontact_stage_interval_m().unwrap();
        let stage = 0.5 * minimum + 0.5 * touch;
        let result = f.evaluate_noncontact(stage).unwrap();
        let actual = result.equilibrium();
        let sigma = law.normal_state().adhesion_energy_j_per_m2() / law.interaction_range_m();
        let modulus = law.normal_state().elastic().reduced_modulus_pa();
        let radius = law.effective_radius_m();
        let c = actual.cohesive_radius_m();
        let a = 1.0 / (2.0 * radius) - PI * sigma / f.stiffness_n_per_m();
        let b = 4.0 * sigma / (PI * modulus);
        close(
            stage,
            a * c * c - b * c - law.interaction_range_m(),
            d0,
            2e-11,
        );
        close(actual.normal_force_n(), -PI * sigma * c * c, p0, 2e-12);
        close(
            actual.central_gap_m(),
            -actual.indentation_m() - 2.0 * sigma * c / modulus,
            d0,
            2e-11,
        );
        assert_eq!(actual.friction_contact_state().unwrap().area_m2(), 0.0);
        assert!(actual.central_gap_m() > 0.0);
        let fold = f.evaluate_noncontact(minimum).unwrap();
        close(
            fold.equilibrium().cohesive_radius_m(),
            b / (2.0 * a),
            c,
            2e-10,
        );
        assert!(fold.stage_tangent_n_per_m().is_none());
        let h = 1e-5 * (touch - minimum);
        let left = f.evaluate_noncontact(stage - h).unwrap();
        let right = f.evaluate_noncontact(stage + h).unwrap();
        close(
            (right.total_potential_j() - left.total_potential_j()) / (2.0 * h),
            actual.normal_force_n(),
            p0,
            2e-6,
        );
        close(
            (right.equilibrium().normal_force_n() - left.equilibrium().normal_force_n())
                / (2.0 * h),
            result.stage_tangent_n_per_m().unwrap(),
            p0 / d0,
            2e-6,
        );
        let disk = f.evaluate_noncontact(touch).unwrap();
        let contact = f.evaluate_contact(touch).unwrap();
        close(
            disk.total_potential_j(),
            contact.total_potential_j(),
            p0 * d0,
            2e-10,
        );
        close(
            disk.equilibrium().normal_force_n(),
            contact.equilibrium().normal_force_n(),
            p0,
            2e-11,
        );
        assert_eq!(
            contact
                .equilibrium()
                .friction_contact_state()
                .unwrap()
                .area_m2(),
            0.0
        );
        assert!(f.evaluate_noncontact(minimum - 1e-6 * d0).is_err());
        assert!(f.evaluate_noncontact(touch + 1e-6 * d0).is_err());
    }
}

#[test]
fn loading_unloading_jumps_release_contact_and_spring_storage_once() {
    for lambda in [0.02, 0.2, 2.0, 5.0] {
        for stiffness in [0.2, 10.0] {
            let f = fixture(lambda, stiffness);
            let (d0, p0) = scales(f.contact());
            let initial = f.separation().stage_displacement_m() - 10.0 * d0;
            let old = f.initialize(MaugisPhase::Detached, initial).unwrap();
            let copy = old;
            let load = f.advance_monotone(&old, 2.0 * d0, budget(&f)).unwrap();
            assert_eq!(load.jump.unwrap().kind(), MaugisSnapKind::AttractionOnset);
            assert_eq!(old, copy);
            assert_eq!(
                load,
                f.advance_monotone(&old, 2.0 * d0, budget(&f)).unwrap()
            );
            let unload = f.advance_monotone(&load.next, initial, budget(&f)).unwrap();
            assert_eq!(unload.jump.unwrap().kind(), MaugisSnapKind::Separation);
            assert_eq!(
                unload.next.response().equilibrium().phase(),
                MaugisPhase::Detached
            );
            close(
                load.work.work_j + unload.work.work_j,
                load.released_energy_j + unload.released_energy_j,
                p0 * d0,
                2e-7,
            );
            for jump in [load.jump.unwrap(), unload.jump.unwrap()] {
                close(
                    jump.released_energy_j(),
                    jump.before().total_potential_j() - jump.after().total_potential_j(),
                    p0 * d0,
                    1e-12,
                );
                assert!(jump.released_energy_j() > 0.0);
                assert_eq!(
                    jump.before().stage_displacement_m(),
                    jump.after().stage_displacement_m()
                );
            }
            let held = f
                .advance_monotone(&unload.next, initial, budget(&f))
                .unwrap();
            assert!(held.jump.is_none());
            assert_eq!(held.released_energy_j, 0.0);
            assert_eq!(held.work.work_j, 0.0);
            let farther = f
                .advance_monotone(&unload.next, initial - d0, budget(&f))
                .unwrap();
            assert!(farther.jump.is_none());
            assert_eq!(farther.work.work_j, 0.0);
        }
    }
}

#[test]
fn source_stiffness_domain_and_work_budget_refusals_preserve_history() {
    let f = fixture(5.0, 1.0);
    let (d0, _) = scales(f.contact());
    for k in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(MaugisSpringFixture::new(f.contact().clone(), k).is_err());
    }
    let old = f.initialize(MaugisPhase::Contact, 0.0).unwrap();
    let copy = old;
    let different = fixture(5.0, 2.0);
    assert!(matches!(
        different.advance_monotone(&old, d0, budget(&different)),
        Err(MaugisFixtureError::HistoryMismatch)
    ));
    assert!(f.evaluate_contact(f.maximum_stage_m() + d0).is_err());
    assert!(f.evaluate_contact(f64::NAN).is_err());
    assert!(f.evaluate_noncontact(0.0).is_err());
    assert!(f.evaluate_detached(0.0).is_err());
    let maximum = f
        .advance_monotone(&old, f.maximum_stage_m(), budget(&f))
        .unwrap();
    assert_eq!(
        maximum.next.response().stage_displacement_m(),
        f.maximum_stage_m()
    );
    assert!(maximum.energy_residual_j.abs() <= maximum.energy_tolerance_j);
    let mut short = budget(&f);
    short.work.maximum_evaluations = 5;
    short.work.absolute_tolerance_j *= 1e-10;
    short.work.relative_tolerance = 0.0;
    assert!(matches!(
        f.advance_monotone(&old, f.maximum_stage_m(), short),
        Err(MaugisFixtureError::Contact(MaugisError::WorkBudgetExceeded))
    ));
    assert_eq!(old, copy);
    assert!(f.advance_monotone(&old, d0, budget(&f)).is_ok());
}

#[test]
fn unresolved_contact_disk_crossover_refuses_without_silently_switching_branches() {
    let law = model(0.2);
    let limits = law.noncontact_limits().unwrap().unwrap();
    let c = limits.first_contact_cohesive_radius_m;
    let sigma = law.normal_state().adhesion_energy_j_per_m2() / law.interaction_range_m();
    let r = law.effective_radius_m();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    // At this stiffness, the disk's finite-K fold reaches zero-area touch.
    let critical = 2.0 * PI * sigma * c / (c / r - 4.0 * sigma / (PI * e));
    for factor in [1.0 - 1e-12, 1.0, 1.0 + 1e-12] {
        assert!(MaugisSpringFixture::new(law.clone(), critical * factor).is_err());
    }
}
