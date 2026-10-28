//! G1/G3: low-Tabor contact, finite-gap attraction and conservative work.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::dmt::*;
use fs_tribo::{
    ContactFrame, FrictionLaw, InputAuthority, TangentialSlip, adhesive::AdhesiveFrictionLaw,
};

fn envelope() -> DmtEnvelope {
    DmtEnvelope {
        maximum_indentation_m: 1e-9,
        half_space_depth_m: 1e-3,
        layer_thickness_m: 1e-3,
        maximum_interaction_radius_to_radius: 0.1,
        maximum_interaction_radius_to_depth: 0.1,
        maximum_interaction_radius_to_layer: 0.1,
        maximum_nominal_strain: 0.1,
        maximum_repulsive_pressure_to_yield: 0.3,
        maximum_adhesive_stress_to_yield: 0.1,
        maximum_interaction_range_to_radius: 0.01,
        maximum_tabor_parameter: 0.1,
    }
}
fn source(work: f64, modulus_scale: f64, range: f64) -> SourceDmtSphereContact {
    let a = common::solid("A", 1000.0, 1e11 * modulus_scale);
    let b = common::solid("B", 2000.0, 2e11 * modulus_scale);
    let normal = common::bind(&a, &b, None, work).normal_state().clone();
    SourceDmtSphereContact::new(normal, 1e-6, range, envelope()).unwrap()
}
fn model() -> SourceDmtSphereContact {
    source(0.01, 1.0, 1e-9)
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}, error/scale={:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn original_dimensional_contact_and_rigid_gap_equations_agree() {
    let model = model();
    let r = model.effective_radius_m();
    let w = model.normal_state().adhesion_energy_j_per_m2();
    let e = model.normal_state().elastic().reduced_modulus_pa();
    let z = model.interaction_range_m();
    let pull = 2.0 * PI * r * w;
    for delta in [-2.0 * z, -z, -0.9 * z, -0.5 * z, 0.0, 0.01 * z, 0.1 * z, z] {
        let actual = model.evaluate(delta).unwrap();
        let (radius, repulsive, force, potential) = if delta > 0.0 {
            let a = (r * delta).sqrt();
            let hertz = 4.0 * e * a.powi(3) / (3.0 * r);
            (
                a,
                hertz,
                hertz - pull,
                8.0 * e * r.sqrt() * delta.powf(2.5) / 15.0 - pull * delta - 0.5 * pull * z,
            )
        } else if delta > -z {
            (
                0.0,
                0.0,
                -pull * (1.0 + delta / z),
                -PI * r * w * (delta + z).powi(2) / z,
            )
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };
        close(actual.contact_radius_m, radius, r, 1e-15);
        close(
            actual.repulsive_normal_force_n,
            repulsive,
            repulsive.max(pull),
            1e-14,
        );
        close(actual.normal_force_n, force, repulsive.max(pull), 1e-14);
        close(
            actual.reversible_energy_j,
            potential,
            potential.abs().max(pull * z),
            2e-14,
        );
        close(actual.contact_area_m2, PI * radius * radius, r * r, 1e-15);
        assert_eq!(
            actual.normal_state_identity(),
            model.normal_state().identity()
        );
        assert_eq!(actual.authority, InputAuthority::SyntheticFixture);
    }
}

#[test]
fn potential_gradient_and_force_tangent_match_each_smooth_branch() {
    let model = model();
    let z = model.interaction_range_m();
    for delta in [-2.0 * z, -0.6 * z, 0.1 * z, 0.8 * z] {
        let h = z * 1e-5;
        let minus = model.evaluate(delta - h).unwrap();
        let plus = model.evaluate(delta + h).unwrap();
        let center = model.evaluate(delta).unwrap();
        let force = (plus.reversible_energy_j - minus.reversible_energy_j) / (2.0 * h);
        let tangent = (plus.normal_force_n - minus.normal_force_n) / (2.0 * h);
        close(
            force,
            center.normal_force_n,
            center.normal_force_n.abs().max(1e-9),
            2e-8,
        );
        close(
            tangent,
            center.tangent_n_per_m.unwrap(),
            center.tangent_n_per_m.unwrap().abs().max(1.0),
            2e-8,
        );
    }
}

#[test]
fn touch_and_range_cutoff_are_continuous_but_keep_one_sided_tangents() {
    let model = model();
    let z = model.interaction_range_m();
    let pull = -model.force_controlled_pull_off_n();
    for corner in [-z, 0.0] {
        let at = model.evaluate(corner).unwrap();
        assert_eq!(at.tangent_n_per_m, None);
        assert_ne!(
            at.one_sided_tangents_n_per_m[0],
            at.one_sided_tangents_n_per_m[1]
        );
        assert!(matches!(
            at.point_port([0., 0., 1.], [0.; 3]),
            Err(DmtError::CornerTangent)
        ));
        for direction in [-1.0, 1.0] {
            let nearby = model.evaluate(corner + direction * z * 1e-8).unwrap();
            close(nearby.normal_force_n, at.normal_force_n, pull, 2e-8);
            close(
                nearby.reversible_energy_j,
                at.reversible_energy_j,
                pull * z,
                2e-8,
            );
        }
    }
    assert!(model.evaluate(-0.5 * z).unwrap().tangent_n_per_m.unwrap() < 0.0);
}

#[test]
fn force_controlled_zero_and_tensile_contacts_retain_real_shear_area() {
    let model = model();
    let normal = model.normal_state();
    let pull = -model.force_controlled_pull_off_n();
    let law = AdhesiveFrictionLaw::new(
        FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        3e6,
        2e6,
    )
    .unwrap();
    let slip =
        TangentialSlip::new(&ContactFrame::new([0., 0., 1.]).unwrap(), [0.1, 0., 0.]).unwrap();
    let r = model.effective_radius_m();
    let e = normal.elastic().reduced_modulus_pa();
    let a0 = (3.0 * r * pull / (4.0 * e)).cbrt();
    let zero = model.evaluate_contact_load(0.0).unwrap();
    close(zero.contact_radius_m, a0, a0, 1e-14);
    for load in [0.0, -0.7 * pull] {
        let response = model.evaluate_contact_load(load).unwrap();
        close(response.normal_force_n, load, pull, 1e-14);
        assert!(response.contact_area_m2 > 0.0);
        let shear = law
            .evaluate(
                normal.elastic().interface().interface(),
                response.friction_contact_state().unwrap(),
                slip,
            )
            .unwrap();
        close(
            shear.static_capacity_n,
            3e6 * response.contact_area_m2,
            shear.static_capacity_n,
            1e-14,
        );
        assert!(shear.dissipated_power_w > 0.0);
    }
    let gap = model.evaluate(-0.5 * model.interaction_range_m()).unwrap();
    assert_eq!(gap.regime, DmtRegime::AttractiveGap);
    assert!(gap.normal_force_n < 0.0);
    let shear = law
        .evaluate(
            normal.elastic().interface().interface(),
            gap.friction_contact_state().unwrap(),
            slip,
        )
        .unwrap();
    assert_eq!(shear.static_capacity_n, 0.0);
    assert_eq!(shear.dissipated_power_w, 0.0);
}

#[test]
fn independent_work_integral_and_closed_cycle_create_no_hysteresis() {
    let model = model();
    let z = model.interaction_range_m();
    let mut integrated = 0.0;
    for (from, to) in [(-2.0 * z, -z), (-z, 0.0), (0.0, z)] {
        let n = 8192;
        let h = (to - from) / n as f64;
        let mut sum = model.evaluate(from).unwrap().normal_force_n
            + model.evaluate(to).unwrap().normal_force_n;
        for i in 1..n {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 })
                * model.evaluate(from + i as f64 * h).unwrap().normal_force_n;
        }
        integrated += sum * h / 3.0;
    }
    let exact = model.reversible_work_between(-2.0 * z, z).unwrap();
    close(integrated, exact, exact.abs(), 2e-8);
    close(
        exact + model.reversible_work_between(z, -2.0 * z).unwrap(),
        0.0,
        exact.abs(),
        1e-15,
    );
    let port = model
        .evaluate(-0.4 * z)
        .unwrap()
        .point_port([0., 0., 1.], [1., 2., 3.])
        .unwrap();
    assert_eq!(port.irreversible_work_j, 0.0);
    assert_eq!(port.dissipated_power_w, 0.0);
    for i in 0..3 {
        assert_eq!(port.action_force_n[i], -port.reaction_force_n[i]);
        assert_eq!(port.action_moment_n_m[i], -port.reaction_moment_n_m[i]);
    }
    close(
        port.action_moment_n_m[0],
        2.0 * port.action_force_n[2],
        port.action_force_n[2].abs(),
        1e-15,
    );
}

#[test]
fn source_modulus_and_work_scaling_and_weak_adhesion_limit_are_physical() {
    let original = model();
    let softer = source(0.01, 0.5, 1e-9);
    let stronger = source(0.02, 1.0, 1e-9);
    let weak = source(1e-10, 1.0, 1e-9);
    let delta = 2e-10;
    let a = original.evaluate(delta).unwrap();
    let b = softer.evaluate(delta).unwrap();
    let c = stronger.evaluate(delta).unwrap();
    close(
        a.contact_radius_m,
        b.contact_radius_m,
        a.contact_radius_m,
        1e-15,
    );
    close(
        a.repulsive_normal_force_n,
        2.0 * b.repulsive_normal_force_n,
        a.repulsive_normal_force_n,
        1e-14,
    );
    close(
        stronger.force_controlled_pull_off_n(),
        2.0 * original.force_controlled_pull_off_n(),
        original.force_controlled_pull_off_n().abs(),
        1e-14,
    );
    close(
        c.repulsive_normal_force_n,
        a.repulsive_normal_force_n,
        a.repulsive_normal_force_n,
        1e-14,
    );
    let tiny = weak.evaluate(delta).unwrap();
    close(
        tiny.normal_force_n,
        tiny.repulsive_normal_force_n,
        tiny.repulsive_normal_force_n,
        1e-8,
    );
    assert_ne!(original.identity(), softer.identity());
    assert_ne!(original.identity(), stronger.identity());
}

#[test]
fn invalid_data_and_non_dmt_operating_regimes_refuse() {
    let model = model();
    let normal = model.normal_state().clone();
    assert!(SourceDmtSphereContact::new(normal.clone(), 1e-6, 1e-12, envelope()).is_err());
    assert!(SourceDmtSphereContact::new(normal.clone(), 1e-6, 1e-7, envelope()).is_err());
    assert!(SourceDmtSphereContact::new(normal.clone(), 1e-6, 0.0, envelope()).is_err());
    let mut bad = envelope();
    bad.maximum_tabor_parameter = 0.2;
    assert!(SourceDmtSphereContact::new(normal.clone(), 1e-6, 1e-9, bad).is_err());
    bad = envelope();
    bad.maximum_repulsive_pressure_to_yield = 1e-4;
    assert!(SourceDmtSphereContact::new(normal, 1e-6, 1e-9, bad).is_err());
    for delta in [f64::NAN, f64::INFINITY, 2e-9] {
        assert!(model.evaluate(delta).is_err());
    }
    assert!(matches!(
        model.evaluate_contact_load(1.01 * model.force_controlled_pull_off_n()),
        Err(DmtError::NoContactEquilibrium { .. })
    ));
    assert!(model.evaluate_contact_load(f64::NAN).is_err());
    let at_pull_off = model
        .evaluate_contact_load(model.force_controlled_pull_off_n())
        .unwrap();
    assert_eq!(at_pull_off.contact_area_m2, 0.0);
    assert_eq!(at_pull_off.regime, DmtRegime::Contact);
    assert!(
        model
            .evaluate(1e-10)
            .unwrap()
            .point_port([0.; 3], [0.; 3])
            .is_err()
    );
}
