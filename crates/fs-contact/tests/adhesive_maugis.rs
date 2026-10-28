//! Physical equations and original-source integration for finite-zone adhesion.
#[path = "source_point/common.rs"]
mod common;
use core::f64::consts::PI;
use fs_contact::normal_patch::maugis::*;
use fs_tribo::InputAuthority;

fn model(lambda: f64) -> SourceMaugisSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let r = 0.01;
    let e = normal.elastic().reduced_modulus_pa();
    let a0 = (PI * 1e-7 * r * r / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / r;
    let envelope = MaugisEnvelope {
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
    };
    SourceMaugisSphereContact::new(normal, r, 2.0 * d0 / (PI * lambda), envelope).unwrap()
}
fn scales(law: &SourceMaugisSphereContact) -> (f64, f64, f64) {
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let r = law.effective_radius_m();
    let a0 = (PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    (a0, a0 * a0 / r, PI * w * r)
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}
#[test]
fn original_maugis_closure_load_and_indentation_equations() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 20.0] {
        let law = model(lambda);
        let (a0, d0, p0) = scales(&law);
        close(law.maugis_parameter(), lambda, 2e-15);
        for d in [0.0, 1.0, 4.0, 16.0] {
            let s = law.evaluate_attached(d * d0).unwrap();
            let a = s.contact_radius_m / a0;
            let m = s.cohesive_radius_m / s.contact_radius_m;
            let u = (m * m - 1.0).sqrt();
            let theta = (1.0 / m).acos();
            let closure = 0.5 * lambda * a * a * ((m * m - 2.0) * theta + u)
                + (4.0 / 3.0) * lambda * lambda * a * (u * theta - m + 1.0);
            close(closure, 1.0, 2e-8);
            close(
                s.normal_force_n / p0,
                a.powi(3) - lambda * a * a * (u + m * m * theta),
                2e-8,
            );
            assert!((s.indentation_m / d0 - (a * a - (4.0 / 3.0) * lambda * a * u)).abs() < 2e-8);
            close(s.repulsive_normal_force_n, p0 * a.powi(3), 3e-14);
            close(s.contact_area_m2, PI * s.contact_radius_m.powi(2), 2e-15);
            assert!(s.cohesive_radius_m > s.contact_radius_m);
        }
    }
}
#[test]
fn resolves_dmt_and_jkr_limits_without_blending_them() {
    for (lambda, tol) in [(1e-7, 8e-4), (1000.0, 8e-4)] {
        let law = model(lambda);
        let (a0, d0, p0) = scales(&law);
        for p in [0.0, 1.0, 5.0] {
            let s = law.evaluate_contact_load(p * p0).unwrap();
            let a3 = if lambda < 1.0 {
                p + 2.0
            } else {
                p + 3.0 + (6.0 * p + 9.0).sqrt()
            };
            let a = a3.cbrt();
            close(s.contact_radius_m / a0, a, tol);
            let d = if lambda < 1.0 {
                a * a
            } else {
                a * a - (8.0 * a / 3.0).sqrt()
            };
            close(s.indentation_m / d0, d, tol);
            assert!((s.normal_force_n / p0 - p).abs() < 2e-12);
        }
        close(
            law.force_controlled_pull_off_n() / p0,
            if lambda < 1.0 { -2.0 } else { -1.5 },
            tol,
        );
    }
}
#[test]
fn branch_endpoints_keep_finite_attraction_and_refuse_detachment_inference() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 20.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let min = law.minimum_attached_indentation_m();
        let end = law.evaluate_attached(min).unwrap();
        assert!(end.normal_force_n < 0.0);
        assert!(end.endpoint.is_some());
        assert!(end.tangent_n_per_m.is_none());
        assert!(matches!(
            end.point_port([0.0, 0.0, 1.0], [0.0; 3]),
            Err(MaugisError::EndpointTangent)
        ));
        if lambda < 0.1 {
            assert_eq!(end.endpoint, Some(MaugisEndpoint::ZeroArea));
            assert_eq!(end.contact_area_m2, 0.0);
        }
        if lambda > 2.0 {
            assert_eq!(end.endpoint, Some(MaugisEndpoint::DisplacementFold));
            assert!(end.contact_area_m2 > 0.0);
        }
        assert!(matches!(
            law.evaluate_attached(min - 1e-6 * d0),
            Err(MaugisError::NoContactEquilibrium { .. })
        ));
        let pull = law.force_controlled_pull_off_n();
        assert!(pull >= -2.00001 * p0 && pull <= -1.49999 * p0);
        let load_end = law.evaluate_contact_load(pull).unwrap();
        assert!(load_end.indentation_m >= min);
        assert!(load_end.contact_area_m2 > 0.0);
        assert!(matches!(
            law.evaluate_contact_load(pull - 1e-6 * p0),
            Err(MaugisError::BelowPullOff { .. })
        ));
    }
}
#[test]
fn analytic_tangent_matches_independent_force_differences() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        let law = model(lambda);
        let (_, d0, _) = scales(&law);
        for delta in [0.0, 0.5, 3.0, 10.0] {
            let d = delta * d0;
            let h = 1e-5 * d0;
            let c = law.evaluate_attached(d).unwrap();
            let l = law.evaluate_attached(d - h).unwrap();
            let r = law.evaluate_attached(d + h).unwrap();
            close(
                c.tangent_n_per_m.unwrap(),
                (r.normal_force_n - l.normal_force_n) / (2.0 * h),
                5e-7,
            );
        }
    }
}
#[test]
fn radial_pressure_integrates_to_original_signed_force() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let s = law.evaluate_attached(2.0 * d0).unwrap();
        let n = 8192;
        let h = 0.5 * PI / n as f64;
        let mut integral = 0.0;
        for i in 0..=n {
            let theta = i as f64 * h;
            let r = s.contact_radius_m * theta.sin();
            let f = 2.0 * PI * s.pressure_pa(r).unwrap() * r * s.contact_radius_m * theta.cos();
            integral += if i == 0 || i == n {
                f
            } else if i % 2 == 0 {
                2.0 * f
            } else {
                4.0 * f
            };
        }
        integral *= h / 3.0;
        let sigma = law.normal_state().adhesion_energy_j_per_m2() / law.interaction_range_m();
        integral -= PI
            * sigma
            * (s.cohesive_radius_m - s.contact_radius_m)
            * (s.cohesive_radius_m + s.contact_radius_m);
        assert!(
            (integral - s.normal_force_n).abs() < 2e-7 * p0,
            "lambda={lambda}: {integral:e}, {:?}",
            s.normal_force_n
        );
        close(s.pressure_pa(s.contact_radius_m).unwrap(), -sigma, 2e-15);
        assert_eq!(s.pressure_pa(s.cohesive_radius_m * 1.01).unwrap(), 0.0);
    }
}
#[test]
fn branch_work_matches_independent_quadrature_and_reverses() {
    for lambda in [0.2, 1.0, 5.0] {
        let law = model(lambda);
        let (_, d0, p0) = scales(&law);
        let budget = MaugisWorkBudget {
            absolute_tolerance_j: 1e-11 * p0 * d0,
            relative_tolerance: 1e-10,
            maximum_evaluations: 10000,
        };
        let work = law.reversible_work_between(0.0, 6.0 * d0, budget).unwrap();
        let reverse = law.reversible_work_between(6.0 * d0, 0.0, budget).unwrap();
        assert_eq!(work.work_j, -reverse.work_j);
        let left = law.reversible_work_between(0.0, 2.0 * d0, budget).unwrap();
        let right = law
            .reversible_work_between(2.0 * d0, 6.0 * d0, budget)
            .unwrap();
        close(work.work_j, left.work_j + right.work_j, 3e-9);
        let n = 4096;
        let h = 6.0 * d0 / n as f64;
        let mut reference = 0.0;
        for i in 0..=n {
            let f = law.evaluate_attached(i as f64 * h).unwrap().normal_force_n;
            reference += if i == 0 || i == n { 0.5 * f } else { f };
        }
        close(work.work_j, reference * h, 2e-6);
        assert!(work.estimated_absolute_error_j < 1e-7 * p0 * d0);
        assert!(work.evaluations <= budget.maximum_evaluations);
        let same = law.reversible_work_between(d0, d0, budget).unwrap();
        assert_eq!(same.work_j, 0.0);
        // Distinct inputs must never collapse into an exact zero-work claim.
        match law.reversible_work_between(d0, f64::from_bits(d0.to_bits() + 1), budget) {
            Ok(v) => assert_ne!(v.work_j, 0.0),
            Err(MaugisError::Unrepresentable { .. } | MaugisError::WorkBudgetExceeded) => {}
            x => panic!("{x:?}"),
        }
        assert!(matches!(
            law.reversible_work_between(
                0.0,
                6.0 * d0,
                MaugisWorkBudget {
                    maximum_evaluations: 5,
                    ..budget
                }
            ),
            Err(MaugisError::WorkBudgetExceeded)
        ));
    }
}
#[test]
fn original_sources_real_area_signed_ports_and_applicability_are_preserved() {
    let law = model(1.0);
    let (_, d0, p0) = scales(&law);
    let s = law.evaluate_contact_load(-0.5 * p0).unwrap();
    assert_eq!(s.normal_state_identity(), law.normal_state().identity());
    assert_eq!(s.authority, InputAuthority::SyntheticFixture);
    assert!(s.contact_area_m2 > 0.0 && s.repulsive_normal_force_n > 0.0);
    s.friction_contact_state().unwrap();
    let port = s.point_port([0.0, 0.0, 1.0], [0.2, 0.0, 0.0]).unwrap();
    assert!(port.action_force_n[2] < 0.0);
    for i in 0..3 {
        assert_eq!(port.action_force_n[i], -port.reaction_force_n[i]);
        assert_eq!(port.action_moment_n_m[i], -port.reaction_moment_n_m[i]);
    }
    assert_eq!(port.irreversible_work_j, 0.0);
    assert_eq!(port.dissipated_power_w, 0.0);
    assert!(matches!(
        law.evaluate_attached(21.0 * d0),
        Err(MaugisError::OutsideApplicability { .. })
    ));
    assert!(matches!(
        law.evaluate_attached(f64::NAN),
        Err(MaugisError::InvalidInput { .. })
    ));
    assert!(s.pressure_pa(-1.0).is_err());
    let mut env = law.envelope();
    env.maximum_adhesive_stress_to_yield = 1e-15;
    assert!(matches!(
        SourceMaugisSphereContact::new(
            law.normal_state().clone(),
            law.effective_radius_m(),
            law.interaction_range_m(),
            env
        ),
        Err(MaugisError::OutsideApplicability { .. })
    ));
}
