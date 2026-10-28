//! G1/G3: load control, finite machine stiffness, stability and physical cycle work.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::fixture::{JkrLoadingError, JkrSpringFixture};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrEnvelope, JkrJumpKind, SourceJkrSphereContact,
};

fn model() -> SourceJkrSphereContact {
    let a = common::solid("fixture-A", 1000.0, 1e7);
    let b = common::solid("fixture-B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.05).normal_state().clone();
    SourceJkrSphereContact::new(
        normal,
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
fn scales(contact: &SourceJkrSphereContact) -> (f64, f64, f64) {
    let h = -contact.minimum_attached_indentation_m() / 3.0;
    let f = -contact.force_controlled_pull_off_n() / 9.0;
    (h, f, f / h)
}
fn near(actual: f64, expected: f64, relative: f64, reference: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()).max(reference),
        "{actual:e} differs from {expected:e} (reference {reference:e})"
    );
}

#[test]
fn force_control_uses_stable_radius_and_preserves_original_contact_sources() {
    let law = model();
    let (h, f0, _) = scales(&law);
    let minimum = law.force_controlled_pull_off_n();
    let maximum = law
        .evaluate_attached(law.envelope().maximum_indentation_m)
        .unwrap()
        .normal_force_n;
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let r = law.effective_radius_m();
    for load in [
        minimum,
        minimum.next_up(),
        -8.0 * f0,
        0.0,
        0.5 * maximum,
        maximum,
    ] {
        let response = law.evaluate_force_controlled(load).unwrap();
        // Original dimensional load-radius relation, independently reconstructed.
        let radius = ((3.0 * r / (4.0 * e))
            * (load + 3.0 * PI * w * r + (6.0 * PI * w * r * (load - minimum)).sqrt()))
        .cbrt();
        near(response.contact_radius_m, radius, 4e-12, radius);
        near(response.normal_force_n, load, 2e-14, f0);
        near(
            response.indentation_m,
            radius * radius / r - (2.0 * PI * w * radius / e).sqrt(),
            5e-12,
            h,
        );
        assert!(response.tangent_n_per_m.unwrap() >= 0.0);
        assert_eq!(
            response.normal_state_identity(),
            law.normal_state().identity()
        );
    }
    assert_eq!(
        law.evaluate_force_controlled(minimum)
            .unwrap()
            .tangent_n_per_m,
        Some(0.0)
    );
    assert!(matches!(
        law.evaluate_force_controlled(minimum.next_down()),
        Err(JkrLoadingError::NoForceEquilibrium { .. })
    ));
    assert!(law.evaluate_force_controlled(maximum.next_up()).is_err());
    assert!(law.evaluate_force_controlled(f64::NAN).is_err());
}

#[test]
fn finite_machine_stiffness_sets_the_actual_stability_fold() {
    let law = model();
    let (h, f0, k0) = scales(&law);
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let r = law.effective_radius_m();
    let mut previous_force = -10.0 * f0;
    for ratio in [1e-3, 0.02, 1.0, 50.0, 1e4] {
        let k = ratio * k0;
        let fixture = JkrSpringFixture::new(law.clone(), k).unwrap();
        let critical = fixture.critical_contact();
        assert!(critical.normal_force_n >= -9.0 * f0 - 1e-14 * f0);
        assert!(critical.normal_force_n <= -5.0 * f0 + 1e-14 * f0);
        assert!(critical.normal_force_n > previous_force);
        previous_force = critical.normal_force_n;
        assert!(critical.indentation_m > law.minimum_attached_indentation_m());
        assert!(critical.indentation_m < 0.0);
        near(critical.tangent_n_per_m.unwrap(), -k, 5e-10, k);
        near(
            fixture.minimum_attached_stage_m(),
            critical.indentation_m + critical.normal_force_n / k,
            2e-14,
            h,
        );
        // Independent Hessian of U(a,delta)+K(D-delta)^2/2.
        let a = critical.contact_radius_m;
        let delta = critical.indentation_m;
        let aa = -4.0 * e * a / r * (delta - a * a / r) - 2.0 * PI * w;
        let ad = 2.0 * e * (delta - a * a / r);
        let dd = 2.0 * e * a + k;
        assert!(aa > 0.0);
        near(aa * dd, ad * ad, 2e-9, aa * dd + ad * ad);
        let fold = fixture
            .evaluate_attached(fixture.minimum_attached_stage_m())
            .unwrap();
        assert_eq!(fold.stage_tangent_n_per_m, None);
        assert!(fold.stage_port([0., 0., 1.], [0.; 3]).is_err());
        let above = fixture.evaluate_attached(0.0).unwrap().contact;
        let a = above.contact_radius_m;
        let delta = above.indentation_m;
        let aa = -4.0 * e * a / r * (delta - a * a / r) - 2.0 * PI * w;
        let ad = 2.0 * e * (delta - a * a / r);
        assert!(aa * (2.0 * e * a + k) - ad * ad > 0.0);
    }
}

#[test]
fn spring_force_energy_gradient_and_condensed_tangent_agree() {
    let law = model();
    let (h, f0, k0) = scales(&law);
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let r = law.effective_radius_m();
    for ratio in [0.3, 5.0] {
        let k = ratio * k0;
        let fixture = JkrSpringFixture::new(law.clone(), k).unwrap();
        for drive in [
            fixture.minimum_attached_stage_m() * 0.5,
            0.0,
            fixture.maximum_stage_m() * 0.2,
            fixture.maximum_stage_m() * 0.6,
        ] {
            let response = fixture.evaluate_attached(drive).unwrap();
            let contact = response.contact;
            near(
                contact.normal_force_n,
                k * (drive - contact.indentation_m),
                2e-12,
                f0,
            );
            near(response.equilibrium_residual_n, 0.0, 2e-12, f0);
            near(
                response.spring_energy_j,
                0.5 * k * (drive - contact.indentation_m).powi(2),
                2e-12,
                f0 * h,
            );
            let a = contact.contact_radius_m;
            let delta = contact.indentation_m;
            let normal_energy = e * a * delta * delta - 2.0 * e * a.powi(3) * delta / (3.0 * r)
                + e * a.powi(5) / (5.0 * r * r)
                - PI * w * a * a;
            near(
                response.total_stored_energy_j,
                normal_energy + 0.5 * k * (drive - delta).powi(2),
                4e-12,
                f0 * h,
            );
            let step = h * 1e-5;
            let plus = fixture.evaluate_attached(drive + step).unwrap();
            let minus = fixture.evaluate_attached(drive - step).unwrap();
            near(
                (plus.total_stored_energy_j - minus.total_stored_energy_j) / (2.0 * step),
                contact.normal_force_n,
                2e-6,
                f0,
            );
            near(
                (plus.contact.normal_force_n - minus.contact.normal_force_n) / (2.0 * step),
                response.stage_tangent_n_per_m.unwrap(),
                2e-6,
                k0,
            );
            let kc = contact.tangent_n_per_m.unwrap();
            near(
                response.stage_tangent_n_per_m.unwrap(),
                k * kc / (k + kc),
                2e-11,
                k0,
            );
            let port = response.stage_port([0., 0., 1.], [0.02, 0., 0.]).unwrap();
            near(port.action_force_n[2], contact.normal_force_n, 1e-15, f0);
            near(
                port.action_moment_n_m[1],
                -0.02 * contact.normal_force_n,
                1e-15,
                0.02 * f0,
            );
            assert_eq!(port.reaction_force_n, port.action_force_n.map(|f| -f));
            assert_eq!(
                port.tangent_n_per_m,
                response.stage_tangent_n_per_m.unwrap()
            );
            assert_eq!(port.irreversible_work_j, 0.0);
            assert_eq!(port.dissipated_power_w, 0.0);
        }
    }
}

#[test]
fn complete_drive_cycle_matches_independent_radius_work_and_event_energy() {
    let law = model();
    let (_, f0, k0) = scales(&law);
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let r = law.effective_radius_m();
    for ratio in [0.1, 1.0, 100.0] {
        let k = ratio * k0;
        let fixture = JkrSpringFixture::new(law.clone(), k).unwrap();
        let far = 1.2 * fixture.minimum_attached_stage_m();
        let initial = fixture.detached_state(far).unwrap();
        let loading = fixture
            .advance_monotone(initial, fixture.maximum_stage_m() * 0.4)
            .unwrap();
        let attach = loading.jump.unwrap();
        assert_eq!(attach.kind, JkrJumpKind::Attach);
        assert_eq!(attach.stage_displacement_m, 0.0);
        assert_eq!(attach.before_indentation_m, 0.0);
        assert!(attach.after_indentation_m > 0.0);
        assert!(attach.released_energy_j > 0.0);
        let unloading = fixture
            .advance_monotone(loading.response.state(), far)
            .unwrap();
        let detach = unloading.jump.unwrap();
        assert_eq!(detach.kind, JkrJumpKind::Detach);
        assert_eq!(
            detach.stage_displacement_m,
            fixture.minimum_attached_stage_m()
        );
        assert!(detach.before_indentation_m < 0.0);
        assert!(detach.after_indentation_m < detach.before_indentation_m);
        assert_eq!(unloading.response.contact.branch, JkrBranch::Detached);
        assert_eq!(unloading.response.total_stored_energy_j, 0.0);
        let work = loading.prescribed_work_j + unloading.prescribed_work_j;
        near(
            work,
            attach.released_energy_j + detach.released_energy_j,
            5e-14,
            work,
        );
        // Loading and unloading cancel above D=0. Integrate the original
        // dimensional F(a)*(delta'(a)+F'(a)/K) from attachment to separation.
        let from = fixture
            .evaluate_attached(0.0)
            .unwrap()
            .contact
            .contact_radius_m;
        let to = fixture.critical_contact().contact_radius_m;
        let count = 2048;
        let step = (to - from) / f64::from(count);
        let c = (2.0 * PI * w / e).sqrt();
        let integrand = |a: f64| {
            let force = 4.0 * e * a.powi(3) / (3.0 * r) - 2.0 * e * c * a.powf(1.5);
            let delta_prime = 2.0 * a / r - c / (2.0 * a.sqrt());
            let force_prime = 4.0 * e * a * a / r - 3.0 * e * c * a.sqrt();
            force * (delta_prime + force_prime / k)
        };
        let mut sum = integrand(from) + integrand(to);
        for i in 1..count {
            sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(from + f64::from(i) * step);
        }
        let integral = step * sum / 3.0;
        near(
            work,
            integral,
            3e-10,
            f0 * (-law.minimum_attached_indentation_m()),
        );
        let again = fixture
            .advance_monotone(unloading.response.state(), far * 1.1)
            .unwrap();
        assert_eq!(again.jump, None);
        assert_eq!(again.irreversible_work_j, 0.0);
        assert_eq!(again.prescribed_work_j, 0.0);
    }
}

#[test]
fn soft_and_stiff_fixture_limits_and_force_query_are_consistent() {
    let law = model();
    let (h, f0, k0) = scales(&law);
    let soft = JkrSpringFixture::new(law.clone(), 1e-3 * k0).unwrap();
    let stiff = JkrSpringFixture::new(law.clone(), 1e4 * k0).unwrap();
    near(soft.critical_contact().normal_force_n, -9.0 * f0, 2e-6, f0);
    near(stiff.critical_contact().normal_force_n, -5.0 * f0, 4e-4, f0);
    near(stiff.critical_contact().indentation_m, -3.0 * h, 1e-6, h);
    assert!(stiff.evaluate_attached(0.0).unwrap().contact.indentation_m < 0.002 * h);
    for fixture in [soft, stiff] {
        let direct = law.evaluate_force_controlled(-6.0 * f0).unwrap();
        let drive = direct.indentation_m + direct.normal_force_n / fixture.stiffness_n_per_m();
        let coupled = fixture.evaluate_attached(drive).unwrap();
        near(
            coupled.contact.normal_force_n,
            direct.normal_force_n,
            2e-10,
            f0,
        );
        near(
            coupled.contact.contact_radius_m,
            direct.contact_radius_m,
            2e-10,
            direct.contact_radius_m,
        );
    }
}

#[test]
fn invalid_loading_and_foreign_history_refuse_without_consuming_state() {
    let law = model();
    let (_, _, k0) = scales(&law);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e12 * k0] {
        assert!(JkrSpringFixture::new(law.clone(), bad).is_err());
    }
    let fixture = JkrSpringFixture::new(law.clone(), k0).unwrap();
    let other = JkrSpringFixture::new(law, 2.0 * k0).unwrap();
    let previous = fixture
        .detached_state(fixture.minimum_attached_stage_m() * 1.1)
        .unwrap();
    assert!(matches!(
        other.advance_monotone(previous, 0.0),
        Err(JkrLoadingError::StateMismatch)
    ));
    assert!(matches!(
        fixture.evaluate_attached(fixture.minimum_attached_stage_m().next_down()),
        Err(JkrLoadingError::NoFixtureEquilibrium { .. })
    ));
    assert!(matches!(
        fixture.advance_monotone(previous, fixture.maximum_stage_m().next_up()),
        Err(JkrLoadingError::StageOutsideEnvelope { .. })
    ));
    assert!(fixture.detached_state(1e-9).is_err());
    assert!(fixture.advance_monotone(previous, f64::NAN).is_err());
    assert_eq!(previous.branch(), JkrBranch::Detached);
    let accepted = fixture.advance_monotone(previous, 0.0).unwrap();
    assert_eq!(accepted.response.state().branch(), JkrBranch::Attached);
    let smooth = fixture
        .advance_monotone(accepted.response.state(), fixture.maximum_stage_m() * 0.01)
        .unwrap();
    assert_eq!(smooth.jump, None);
    assert_eq!(smooth.irreversible_work_j, 0.0);
}

#[test]
fn adjacent_stage_above_fold_resolves_attached_area_without_a_fake_smooth_fold() {
    let law = model();
    let (_, _, k0) = scales(&law);
    for ratio in [1e-3, 1.0, 1e3] {
        let k = ratio * k0;
        let fixture = JkrSpringFixture::new(law.clone(), k).unwrap();
        let fold = fixture.minimum_attached_stage_m();
        let above = fixture.evaluate_attached(fold.next_up()).unwrap();
        assert!(above.contact.contact_radius_m > fixture.critical_contact().contact_radius_m);
        assert!(above.contact.indentation_m >= fixture.critical_contact().indentation_m);
        assert!(above.stage_tangent_n_per_m.unwrap().is_finite());
        let scale = above.contact.normal_force_n.abs()
            + k * above.stage_displacement_m.abs()
            + k * above.contact.indentation_m.abs();
        near(
            above.equilibrium_residual_n,
            0.0,
            128.0 * f64::EPSILON,
            scale,
        );
        assert_eq!(
            fixture
                .evaluate_attached(fold)
                .unwrap()
                .stage_tangent_n_per_m,
            None
        );
    }
}
