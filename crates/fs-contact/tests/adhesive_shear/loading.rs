//! Real source strengths compete with reversible area loss and machine stability.
use super::*;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::{
    JkrMixedModeEnvelope, SourceJkrMixedMode,
    loading::{
        JkrShearBoundaryKind, JkrShearLoadingError, JkrShearLoadingLimits, JkrShearSpringFixture,
    },
};

fn mixed(strength: f64, length: f64) -> (SourceJkrMixedMode, SourceAdhesiveShear) {
    // At fixed E,W, the JKR similarity is R~l^3, a~l^2, delta/u~l,
    // F~l^3, K~l^2, U~l^4. Source shear strength ~1/l preserves T/(tau A).
    let normal = normal(
        Some((strength / length, 0.5 * strength / length)),
        Pressure::DIMS,
        false,
        false,
    );
    let shear =
        SourceAdhesiveShear::resolve(&normal, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let normal = SourceJkrSphereContact::new(
        normal,
        0.01 * length.powi(3),
        1e-9 * length,
        JkrEnvelope {
            maximum_indentation_m: 1e-5 * length,
            half_space_depth_m: 0.1 * length * length,
            layer_thickness_m: 0.1 * length * length,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap();
    (
        SourceJkrMixedMode::new(
            normal,
            JkrMixedModeEnvelope {
                maximum_tangential_displacement_to_radius: 0.1,
                maximum_mean_shear_stress_to_yield: 0.1,
            },
        )
        .unwrap(),
        shear,
    )
}
fn limits(length: f64) -> JkrShearLoadingLimits {
    JkrShearLoadingLimits {
        maximum_iterations: 80,
        maximum_absolute_stage_m: 1.0 * length,
    }
}
fn fixture(strength: f64, stiffness: f64) -> JkrShearSpringFixture {
    let (normal, shear) = mixed(strength, 1.0);
    JkrShearSpringFixture::new(normal, shear, 2e-6, stiffness, limits(1.0)).unwrap()
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scale={scale:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn source_strength_selects_slip_or_peeling_on_the_actual_reduced_area() {
    for (strength, kind) in [
        (1000.0, JkrShearBoundaryKind::GrossSlip),
        (3e6, JkrShearBoundaryKind::RadiusInstability),
    ] {
        let law = fixture(strength, 2000.0);
        let initial = law.initial_response();
        let boundary = law.first_boundary();
        assert_eq!(boundary.kind, kind);
        let response = boundary.response;
        let contact = response.contact;
        let handoff = law
            .advance_until_event(initial.state(), 1.5 * response.stage_displacement_m)
            .unwrap();
        assert_eq!(handoff.boundary, Some(kind));
        assert_eq!(handoff.response, response);
        assert!(handoff.untraversed_stage_displacement_m > 0.0);
        assert!(law.evaluate(1.5 * response.stage_displacement_m).is_err());
        assert!(contact.contact_area_m2 < initial.contact.contact_area_m2);
        let actual = law
            .shear_source()
            .law()
            .static_capacity_n(contact.friction_contact_state().unwrap())
            .unwrap();
        assert_eq!(response.static_capacity_n, actual);
        assert!(response.static_margin_n >= 0.0);
        close(
            contact.tangential_force_n[0],
            law.stiffness_n_per_m()
                * (response.stage_displacement_m - contact.tangential_displacement_m[0]),
            contact.tangential_force_n[0],
            3e-15,
        );
        assert_eq!(
            contact.normal_state_identity(),
            law.shear_source().normal_state_identity()
        );
        let held = law
            .shear_source()
            .bind_mixed_mode(law.contact().normal_source().normal_state(), &contact)
            .unwrap();
        assert_eq!(
            held.tangential().contact().area_m2(),
            contact.contact_area_m2
        );
        if kind == JkrShearBoundaryKind::GrossSlip {
            close(contact.tangential_force_n[0], actual, actual, 1e-14);
            assert!(response.fixture_radius_curvature_n_per_m > 0.0);
            assert!(response.stage_tangent_n_per_m.is_some());
            // Using the initial, larger area would permit an inadmissible load.
            assert!(actual < strength * initial.contact.contact_area_m2);
        } else {
            assert!(actual > contact.tangential_force_n[0]);
            close(
                response.fixture_radius_curvature_n_per_m,
                0.0,
                PI * 0.05,
                2e-12,
            );
            assert_eq!(response.stage_tangent_n_per_m, None);
        }
    }
}

#[test]
fn apparatus_stiffness_changes_the_peeling_fold_with_correct_soft_and_rigid_limits() {
    let mut previous_shift = 0.0;
    let mut previous_radius = f64::INFINITY;
    for stiffness in [0.1, 100.0, 2000.0, 1e5, 1e7] {
        let law = fixture(3e6, stiffness);
        let response = law.first_boundary().response;
        let contact = response.contact;
        let a = contact.contact_radius_m;
        let u = contact.tangential_displacement_m[0];
        let k = 8.0 * law.contact().reduced_shear_modulus_pa();
        let e = law
            .contact()
            .normal_source()
            .normal_state()
            .elastic()
            .reduced_modulus_pa();
        let r = law.contact().normal_source().effective_radius_m();
        let w = law
            .contact()
            .normal_source()
            .normal_state()
            .adhesion_energy_j_per_m2();
        // Independent Hessian of U(a,delta,u)+K(D-u)^2/2 at fixed delta,D.
        let haa = 4.0 * e * a * (a * a / r - law.indentation_m()) / r - 2.0 * PI * w;
        let hau = k * u;
        close(
            haa - hau * hau / (k * a + stiffness),
            response.fixture_radius_curvature_n_per_m,
            PI * w,
            1e-12,
        );
        assert!(u > previous_shift && a < previous_radius);
        previous_shift = u;
        previous_radius = a;
        if stiffness == 0.1 {
            let contact_tangent = contact.coupled_stiffness_n_per_m.unwrap()[1][1];
            close(contact_tangent, -stiffness, k * a, 2e-12);
        }
        if stiffness == 1e7 {
            let rigid = law
                .contact()
                .critical_tangential_displacement_m(law.indentation_m())
                .unwrap();
            close(u, rigid, rigid, 2e-7);
            assert!(u < rigid);
        }
    }
}

#[test]
fn finite_stage_equilibrium_energy_gradient_and_tangent_share_the_original_potential() {
    let law = fixture(3e6, 2000.0);
    let limit = law.first_boundary().response.stage_displacement_m;
    for fraction in [0.05, 0.4, 0.8, 0.95] {
        let stage = fraction * limit;
        let center = law.evaluate(stage).unwrap();
        assert!(center.fixture_radius_curvature_n_per_m > 0.0);
        let h = 1e-5 * limit;
        let left = law.evaluate(stage - h).unwrap();
        let right = law.evaluate(stage + h).unwrap();
        let force = center.contact.tangential_force_n[0];
        close(
            (right.total_stored_energy_j - left.total_stored_energy_j) / (2.0 * h),
            force,
            force,
            2e-8,
        );
        let stiffness = center.stage_tangent_n_per_m.unwrap();
        close(
            (right.contact.tangential_force_n[0] - left.contact.tangential_force_n[0]) / (2.0 * h),
            stiffness,
            law.stiffness_n_per_m(),
            2e-7,
        );
        close(
            center.spring_energy_j,
            force * force / (2.0 * law.stiffness_n_per_m()),
            center.spring_energy_j,
            3e-15,
        );
        assert_eq!(
            center.total_stored_energy_j,
            center.spring_energy_j + center.contact.reversible_energy_j
        );
        close(center.force_balance_residual_n, 0.0, force, 3e-14);
    }
    let adjacent = law.evaluate(limit.next_down()).unwrap();
    assert!(adjacent.fixture_radius_curvature_n_per_m > 0.0);
    assert!(adjacent.stage_tangent_n_per_m.is_some());
    assert!(law.evaluate(limit.next_up()).is_err());
    let tiny = law.evaluate(1e-30 * limit).unwrap();
    assert!(tiny.contact.tangential_force_n[0] > 0.0);
    close(
        tiny.force_balance_residual_n,
        0.0,
        tiny.contact.tangential_force_n[0],
        3e-14,
    );
}

#[test]
fn exact_work_survives_path_partition_unloading_reversal_and_event_handoff() {
    let law = fixture(3e6, 2000.0);
    let initial = law.initial_response();
    let limit = law.first_boundary().response.stage_displacement_m;
    let target = 0.85 * limit;
    let full = law.advance_until_event(initial.state(), target).unwrap();
    let first = law
        .advance_until_event(initial.state(), 0.3 * limit)
        .unwrap();
    let second = law
        .advance_until_event(first.response.state(), target)
        .unwrap();
    close(
        first.prescribed_work_j + second.prescribed_work_j,
        full.prescribed_work_j,
        full.prescribed_work_j.abs(),
        2e-15,
    );
    let n = 512;
    let mut integral = law.evaluate(0.0).unwrap().contact.tangential_force_n[0]
        + law.evaluate(target).unwrap().contact.tangential_force_n[0];
    for i in 1..n {
        integral += (if i % 2 == 0 { 2.0 } else { 4.0 })
            * law
                .evaluate(target * f64::from(i) / f64::from(n))
                .unwrap()
                .contact
                .tangential_force_n[0];
    }
    integral *= target / (3.0 * f64::from(n));
    close(
        integral,
        full.prescribed_work_j,
        full.prescribed_work_j.abs(),
        2e-9,
    );
    let reversed = law
        .advance_until_event(full.response.state(), -0.4 * limit)
        .unwrap();
    let returned = law
        .advance_until_event(reversed.response.state(), 0.0)
        .unwrap();
    close(
        full.prescribed_work_j + reversed.prescribed_work_j + returned.prescribed_work_j,
        0.0,
        full.prescribed_work_j.abs(),
        3e-15,
    );
    assert_eq!(returned.response, initial);
    let handoff = law
        .advance_until_event(initial.state(), 1.4 * limit)
        .unwrap();
    assert_eq!(
        handoff.boundary,
        Some(JkrShearBoundaryKind::RadiusInstability)
    );
    assert_eq!(handoff.response, law.first_boundary().response);
    close(
        handoff.untraversed_stage_displacement_m,
        0.4 * limit,
        limit,
        3e-15,
    );
    assert_eq!(
        handoff.prescribed_work_j,
        handoff.response.total_stored_energy_j - initial.total_stored_energy_j
    );
    assert_eq!(initial, law.initial_response());
    let reflected = law
        .advance_until_event(initial.state(), -1.4 * limit)
        .unwrap();
    assert_eq!(reflected.prescribed_work_j, handoff.prescribed_work_j);
    assert_eq!(
        reflected.response.contact.contact_area_m2,
        handoff.response.contact.contact_area_m2
    );
    assert_eq!(
        reflected.response.contact.tangential_force_n[0],
        -handoff.response.contact.tangential_force_n[0]
    );
    assert_eq!(
        reflected.untraversed_stage_displacement_m,
        -handoff.untraversed_stage_displacement_m
    );
}

#[test]
fn source_similarity_preserves_competing_boundary_and_machine_work() {
    for strength in [1000.0, 3e6] {
        let original = fixture(strength, 2000.0);
        let terminal = original.first_boundary();
        let stage = 0.65 * terminal.response.stage_displacement_m;
        let reference = original
            .advance_until_event(original.initial_response().state(), stage)
            .unwrap();
        for length in [0.5, 2.0, 5.0] {
            let (normal, shear) = mixed(strength, length);
            let scaled = JkrShearSpringFixture::new(
                normal,
                shear,
                2e-6 * length,
                2000.0 * length * length,
                limits(length),
            )
            .unwrap();
            let boundary = scaled.first_boundary();
            assert_eq!(boundary.kind, terminal.kind);
            close(
                boundary.response.stage_displacement_m,
                terminal.response.stage_displacement_m * length,
                terminal.response.stage_displacement_m * length,
                3e-14,
            );
            close(
                boundary.response.contact.contact_radius_m,
                terminal.response.contact.contact_radius_m * length * length,
                terminal.response.contact.contact_radius_m * length * length,
                3e-14,
            );
            let step = scaled
                .advance_until_event(scaled.initial_response().state(), stage * length)
                .unwrap();
            close(
                step.response.contact.tangential_force_n[0],
                reference.response.contact.tangential_force_n[0] * length.powi(3),
                reference.response.contact.tangential_force_n[0] * length.powi(3),
                3e-14,
            );
            close(
                step.prescribed_work_j,
                reference.prescribed_work_j * length.powi(4),
                reference.prescribed_work_j.abs() * length.powi(4),
                2e-11,
            );
        }
    }
}

#[test]
fn capacity_zero_source_mismatch_and_explicit_budgets_do_not_publish_a_trial() {
    let zero = fixture(0.0, 2000.0);
    assert_eq!(zero.first_boundary().kind, JkrShearBoundaryKind::GrossSlip);
    assert_eq!(zero.first_boundary().response.stage_displacement_m, 0.0);
    let stop = zero
        .advance_until_event(zero.initial_response().state(), 1e-6)
        .unwrap();
    assert_eq!(stop.prescribed_work_j, 0.0);
    assert_eq!(stop.untraversed_stage_displacement_m, 1e-6);
    let law = fixture(1000.0, 2000.0);
    let original = law.initial_response();
    for target in [f64::NAN, 2.0] {
        assert!(law.advance_until_event(original.state(), target).is_err());
        assert_eq!(law.initial_response(), original);
    }
    let other = fixture(1000.0, 3000.0);
    assert!(matches!(
        law.advance_until_event(other.initial_response().state(), 0.0),
        Err(JkrShearLoadingError::StateMismatch)
    ));
    let (normal, _) = mixed(1000.0, 1.0);
    let (_, foreign) = mixed(2000.0, 1.0);
    assert!(matches!(
        JkrShearSpringFixture::new(normal, foreign, 2e-6, 2000.0, limits(1.0)),
        Err(JkrShearLoadingError::StateMismatch)
    ));
    let (normal, shear) = mixed(3e6, 1.0);
    let before = normal.clone();
    let mut budget = limits(1.0);
    budget.maximum_iterations = 1;
    assert!(matches!(
        JkrShearSpringFixture::new(normal.clone(), shear, 2e-6, 2000.0, budget),
        Err(JkrShearLoadingError::NoConvergence { .. })
    ));
    assert_eq!(normal, before);
}
