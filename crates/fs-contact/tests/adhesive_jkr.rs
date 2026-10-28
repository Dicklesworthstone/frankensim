//! G1/G3: independent JKR reconstruction, energy, history and friction coupling.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::*;
use fs_tribo::{
    ContactFrame, FrictionLaw, InputAuthority, TangentialSlip, adhesive::AdhesiveFrictionLaw,
};

fn envelope() -> JkrEnvelope {
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
    }
}
fn source(work: f64, modulus_scale: f64, range: f64) -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7 * modulus_scale);
    let b = common::solid("B", 2000.0, 2e7 * modulus_scale);
    let normal = common::bind(&a, &b, None, work).normal_state().clone();
    SourceJkrSphereContact::new(normal, 0.01, range, envelope()).unwrap()
}
fn model() -> SourceJkrSphereContact {
    source(0.05, 1.0, 1e-9)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:e} differs from {expected:e}"
    );
}
fn quantities(law: &SourceJkrSphereContact) -> (f64, f64, f64, f64) {
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let r = law.effective_radius_m();
    let a_s = (PI * w * r * r / (8.0 * e)).cbrt();
    (e, w, r, a_s)
}

#[test]
fn physical_radius_force_area_and_potential_follow_original_jkr_equations() {
    let law = model();
    let (e, w, r, a_s) = quantities(&law);
    let c = (2.0 * PI * w / e).sqrt();
    for multiple in [1.01, 1.2, 2.0, 4.0] {
        let a = multiple * a_s;
        let delta = a * a / r - c * a.sqrt();
        let result = law.evaluate_attached(delta).unwrap();
        close(result.contact_radius_m, a, 3e-12);
        close(result.contact_area_m2, PI * a * a, 5e-12);
        close(
            result.repulsive_normal_force_n,
            4.0 * e * a.powi(3) / (3.0 * r),
            5e-12,
        );
        close(
            result.normal_force_n,
            4.0 * e * a.powi(3) / (3.0 * r) - (8.0 * PI * w * e * a.powi(3)).sqrt(),
            5e-12,
        );
        let potential = 8.0 * e * a.powi(5) / (15.0 * r * r)
            - 4.0 * e * c * a.powf(3.5) / (3.0 * r)
            + PI * w * a * a;
        close(result.reversible_energy_j, potential, 1e-10);
        let tangent =
            (4.0 * e * a * a / r - 3.0 * e * c * a.sqrt()) / (2.0 * a / r - c / (2.0 * a.sqrt()));
        close(result.tangent_n_per_m.unwrap(), tangent, 5e-11);
        assert_eq!(
            result.normal_state_identity(),
            law.normal_state().identity()
        );
        assert_eq!(result.authority, InputAuthority::SyntheticFixture);
    }
}

#[test]
fn energy_gradient_and_force_tangent_agree_on_tensile_and_compressive_branches() {
    let law = model();
    let critical = law.minimum_attached_indentation_m();
    for delta in [0.7 * critical, 0.0, 1e-6, 8e-6] {
        let h = 1e-5 * critical.abs();
        let center = law.evaluate_attached(delta).unwrap();
        let left = law.evaluate_attached(delta - h).unwrap();
        let right = law.evaluate_attached(delta + h).unwrap();
        close(
            (right.reversible_energy_j - left.reversible_energy_j) / (2.0 * h),
            center.normal_force_n,
            2e-8,
        );
        close(
            (right.normal_force_n - left.normal_force_n) / (2.0 * h),
            center.tangent_n_per_m.unwrap(),
            2e-8,
        );
    }
}

#[test]
fn displacement_spinodal_differs_from_force_controlled_pull_off() {
    let law = model();
    let (e, w, r, a_s) = quantities(&law);
    let minimum = law.minimum_attached_indentation_m();
    close(minimum, -3.0 * a_s * a_s / r, 3e-14);
    let terminal = law.evaluate_attached(minimum).unwrap();
    close(terminal.normal_force_n, -(5.0 / 6.0) * PI * w * r, 3e-14);
    assert_eq!(terminal.tangent_n_per_m, None);
    assert!(terminal.reversible_energy_j > 0.0);
    assert!(matches!(
        terminal.point_port([0.0, 0.0, 1.0], [0.0; 3]),
        Err(JkrError::SingularTangent)
    ));
    let adjacent = law.evaluate_attached(minimum.next_up()).unwrap();
    assert!(adjacent.tangent_n_per_m.unwrap().is_finite());
    assert!(adjacent.tangent_n_per_m.unwrap() < 0.0);
    assert!(matches!(
        law.evaluate_attached(minimum.next_down()),
        Err(JkrError::NoAttachedEquilibrium { .. })
    ));
    let a_load = (9.0 * PI * w * r * r / (8.0 * e)).cbrt();
    let delta_load = a_load * a_load / r - (2.0 * PI * w * a_load / e).sqrt();
    let load_minimum = law.evaluate_attached(delta_load).unwrap();
    close(load_minimum.normal_force_n, -1.5 * PI * w * r, 3e-14);
    close(
        law.force_controlled_pull_off_n(),
        load_minimum.normal_force_n,
        3e-14,
    );
    assert!(load_minimum.tangent_n_per_m.unwrap().abs() < 1e-9 * e * a_load);
    assert!(minimum < delta_load);
}

#[test]
fn branch_jumps_preserve_closed_cycle_energy_and_do_not_repeat() {
    let law = model();
    let minimum = law.minimum_attached_indentation_m();
    let detached = law.detached_state(2.0 * minimum).unwrap();
    let touch = law.advance_monotone(detached, 8e-6).unwrap();
    assert_eq!(touch.response.branch, JkrBranch::Attached);
    assert_eq!(touch.jump.unwrap().kind, JkrJumpKind::Attach);
    assert_eq!(touch.jump.unwrap().indentation_m, 0.0);
    close(
        touch.prescribed_work_j,
        touch.response.reversible_energy_j + touch.irreversible_work_j,
        2e-14,
    );
    let unload = law
        .advance_monotone(touch.response.state(), 2.0 * minimum)
        .unwrap();
    assert_eq!(unload.response.branch, JkrBranch::Detached);
    assert_eq!(unload.jump.unwrap().kind, JkrJumpKind::Detach);
    assert_eq!(unload.jump.unwrap().indentation_m, minimum);
    assert_eq!(unload.response.normal_force_n, 0.0);
    assert_eq!(unload.response.contact_area_m2, 0.0);
    assert_eq!(unload.response.reversible_energy_j, 0.0);
    let again = law
        .advance_monotone(unload.response.state(), 3.0 * minimum)
        .unwrap();
    assert_eq!(again.jump, None);
    assert_eq!(again.irreversible_work_j, 0.0);
    assert_eq!(again.prescribed_work_j, 0.0);

    // Independent quadrature of F(a) d(delta(a)), using x=sqrt(a/a_s),
    // avoids the singular endpoint tangent. The loading/unloading portions
    // above touch cancel, leaving the attached return from delta=0 to delta_s.
    let (e, w, r, a_s) = quantities(&law);
    let c = (2.0 * PI * w / e).sqrt();
    let x_touch = 4.0_f64.cbrt();
    let n = 2048;
    let h = (1.0 - x_touch) / f64::from(n);
    let integrand = |x: f64| {
        let a = a_s * x * x;
        let force = 4.0 * e * a.powi(3) / (3.0 * r) - (8.0 * PI * w * e * a.powi(3)).sqrt();
        force * (4.0 * a_s * a_s * x.powi(3) / r - c * a_s.sqrt())
    };
    let mut integral = integrand(x_touch) + integrand(1.0);
    for i in 1..n {
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(x_touch + f64::from(i) * h);
    }
    integral *= h / 3.0;
    let loss = touch.irreversible_work_j + unload.irreversible_work_j;
    assert!(loss > 0.0);
    close(
        touch.prescribed_work_j + unload.prescribed_work_j,
        loss,
        3e-13,
    );
    close(integral, loss, 2e-11);

    // A smooth reversal retraces its conservative branch and dissipates nothing.
    let a = law.evaluate_attached(4e-6).unwrap();
    let b = law.advance_monotone(a.state(), 6e-6).unwrap();
    let c = law.advance_monotone(b.response.state(), 4e-6).unwrap();
    assert_eq!(b.irreversible_work_j + c.irreversible_work_j, 0.0);
    assert_eq!(b.prescribed_work_j, -c.prescribed_work_j);
}

#[test]
fn zero_net_load_and_tensile_junctions_have_real_area_and_friction() {
    let law = model();
    let (e, w, r, _) = quantities(&law);
    let a0 = (9.0 * PI * w * r * r / (2.0 * e)).cbrt();
    let unloaded = law.evaluate_attached(a0 * a0 / (3.0 * r)).unwrap();
    assert!(unloaded.normal_force_n.abs() < 2e-14 * PI * w * r);
    close(unloaded.contact_radius_m, a0, 3e-14);
    let shear = AdhesiveFrictionLaw::new(
        FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        2e6,
        1e6,
    )
    .unwrap();
    let interface = law.normal_state().elastic().interface().interface();
    let slip = TangentialSlip::new(
        &ContactFrame::new([0.0, 0.0, 1.0]).unwrap(),
        [0.03, 0.04, 0.0],
    )
    .unwrap();
    for normal in [
        unloaded,
        law.evaluate_attached(0.5 * law.minimum_attached_indentation_m())
            .unwrap(),
    ] {
        let response = shear
            .evaluate(interface, normal.friction_contact_state().unwrap(), slip)
            .unwrap();
        assert!(response.dissipated_power_w > 0.0);
        close(
            response.dissipated_power_w,
            1e6 * normal.contact_area_m2 * 0.05,
            3e-14,
        );
        let port = normal.point_port([0.6, 0.8, 0.0], [0.2, 0.0, 0.0]).unwrap();
        close(port.action_force_n[0], 0.6 * normal.normal_force_n, 3e-14);
        close(
            port.action_moment_n_m[2],
            0.2 * port.action_force_n[1],
            3e-14,
        );
        for axis in 0..3 {
            assert_eq!(port.action_force_n[axis], -port.reaction_force_n[axis]);
        }
        assert_eq!(port.irreversible_work_j, 0.0);
    }
}

#[test]
fn source_scaling_hertz_limit_and_history_identity_remain_physical() {
    let base = model();
    let stronger = source(0.4, 1.0, 1e-9);
    let a = base.evaluate_attached(0.0).unwrap();
    let b = stronger.evaluate_attached(0.0).unwrap();
    close(b.contact_radius_m / a.contact_radius_m, 2.0, 5e-14);
    close(b.normal_force_n / a.normal_force_n, 8.0, 5e-14);
    close(b.reversible_energy_j / a.reversible_energy_j, 32.0, 8e-14);
    assert_ne!(base.identity(), stronger.identity());
    assert!(matches!(
        stronger.advance_monotone(a.state(), 0.0),
        Err(JkrError::StateMismatch)
    ));
    let stiffer = source(0.05, 8.0, 1e-9).evaluate_attached(0.0).unwrap();
    close(stiffer.contact_radius_m / a.contact_radius_m, 0.5, 5e-14);
    close(
        stiffer.reversible_energy_j / a.reversible_energy_j,
        0.25,
        8e-14,
    );
    let weak = source(1e-15, 1.0, 1e-30);
    let delta = 1e-5_f64;
    let weak_result = weak.evaluate_attached(delta).unwrap();
    let modulus = weak.normal_state().elastic().reduced_modulus_pa();
    let k = (4.0 / 3.0) * modulus * weak.effective_radius_m().sqrt();
    close(weak_result.normal_force_n, k * delta.powf(1.5), 1e-6);
    close(
        weak_result.reversible_energy_j,
        0.4 * k * delta.powf(2.5),
        1e-6,
    );
}

#[test]
fn approximation_bounds_and_singular_inputs_refuse_without_mutating_history() {
    let law = model();
    let normal = law.normal_state().clone();
    for (radius, range) in [
        (0.0, 1e-9),
        (-1.0, 1e-9),
        (0.01, 0.0),
        (f64::NAN, 1e-9),
        (0.01, 1.0),
    ] {
        assert!(SourceJkrSphereContact::new(normal.clone(), radius, range, envelope()).is_err());
    }
    for bad in [
        JkrEnvelope {
            layer_thickness_m: 1e-7,
            ..envelope()
        },
        JkrEnvelope {
            maximum_patch_to_radius: 1e-5,
            ..envelope()
        },
        JkrEnvelope {
            maximum_nominal_strain: 1e-5,
            ..envelope()
        },
        JkrEnvelope {
            maximum_repulsive_pressure_to_yield: 1e-8,
            ..envelope()
        },
        JkrEnvelope {
            maximum_process_zone_to_patch: 1e-20,
            ..envelope()
        },
        JkrEnvelope {
            minimum_tabor_parameter: 1.0,
            ..envelope()
        },
    ] {
        assert!(SourceJkrSphereContact::new(normal.clone(), 0.01, 1e-9, bad).is_err());
    }
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    let zero = common::bind(&a, &b, None, 0.0).normal_state().clone();
    assert!(SourceJkrSphereContact::new(zero, 0.01, 1e-9, envelope()).is_err());
    let start = law.detached_state(-1e-5).unwrap();
    for delta in [f64::NAN, f64::INFINITY, 1.1e-5] {
        assert!(law.advance_monotone(start, delta).is_err());
    }
    assert_eq!(start.branch(), JkrBranch::Detached);
    assert_eq!(start.indentation_m(), -1e-5);
    assert!(law.detached_state(1e-6).is_err());
    let r = law.evaluate_attached(0.0).unwrap();
    assert!(r.point_port([0.0, 0.0, 2.0], [0.0; 3]).is_err());
}
