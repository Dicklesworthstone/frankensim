//! G1/G3: original pressure fields, conservative finite patches and singular edges.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeEnvelope, SourceJkrMixedMode};
use fs_contact::normal_patch::{adhesive::*, maugis::*, traction::*};
use fs_tribo::InputAuthority;

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

fn mixed_jkr() -> SourceJkrMixedMode {
    SourceJkrMixedMode::new(
        jkr(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}

#[test]
fn mixed_mode_traction_uses_actual_peeled_area_and_normal_energy_release() {
    let source = mixed_jkr();
    let normal = source.normal_source();
    let delta = 1e-6;
    let critical = source.critical_tangential_displacement_m(delta).unwrap();
    let shift = [0.54 * critical, 0.72 * critical];
    let response = source
        .evaluate_displacement_controlled(delta, shift)
        .unwrap();
    let profile =
        AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, shift).unwrap();
    let pure = AdhesiveTractionProfile::jkr_attached(normal, delta).unwrap();
    assert!(profile.contact_radius_m() < 0.95 * pure.contact_radius_m());
    assert_eq!(profile.contact_radius_m(), response.contact_radius_m);
    assert_eq!(profile.model_identity(), source.identity());
    assert_ne!(profile.model_identity(), pure.model_identity());
    assert_eq!(
        profile.normal_state_identity(),
        normal.normal_state().identity()
    );
    assert_eq!(profile.authority(), InputAuthority::SyntheticFixture);

    let a = response.contact_radius_m;
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    // Independent elastic pressure at fixed approach and equilibrated radius.
    // Normal energy release is W minus the positive shear release; substituting
    // total W into a pure-normal field would overestimate its tensile rim.
    let peak = 2.0 * e * a / (PI * r);
    let q = e * (a * a / r - delta) / (PI * a);
    let pure_q = (2.0 * e * w / (PI * a)).sqrt();
    assert!(q < 0.9 * pure_q);
    for fraction in [0.0_f64, 0.3, 0.8, 0.99] {
        let s = (1.0 - fraction * fraction).sqrt();
        close(
            profile.pressure_pa(fraction * a).unwrap(),
            peak * s - q / s,
            peak + q / s,
            3e-14,
        );
    }
    let whole = profile.annulus(0.0, a).unwrap();
    let hertz = (4.0 / 3.0) * e * a.powi(3) / r;
    let tensile = -2.0 * e * a * (a * a / r - delta);
    close(whole.hertz_force_n, hertz, hertz, 3e-14);
    close(whole.adhesive_force_n, tensile, tensile.abs(), 3e-14);
    close(
        whole.normal_force_n,
        response.normal_force_n,
        hertz - tensile,
        3e-14,
    );
    close(
        whole.solid_contact_area_m2,
        response.contact_area_m2,
        response.contact_area_m2,
        3e-14,
    );
    close(
        PI * a * q * q / (2.0 * e),
        response.normal_energy_release_rate_j_m2,
        w,
        3e-14,
    );
    assert_eq!(
        profile.pressure_pa(a),
        Err(AdhesiveTractionError::SingularEdge)
    );
    assert_eq!(profile.pressure_pa(a.next_up()).unwrap(), 0.0);
}

#[test]
fn mixed_mode_force_traction_retains_tension_and_the_correct_peeling_fold() {
    let source = mixed_jkr();
    let pull_off = source.normal_source().force_controlled_pull_off_n();
    for load in [0.5 * pull_off, 0.0, -2.0 * pull_off] {
        let critical = source.critical_tangential_force_n(load).unwrap();
        for fraction in [0.0, 0.7, 1.0] {
            let shear = [fraction * critical, 0.0];
            let response = source.evaluate_force_controlled(load, shear).unwrap();
            let profile =
                AdhesiveTractionProfile::jkr_mixed_mode_force(&source, load, shear).unwrap();
            let a = profile.contact_radius_m();
            let whole = profile.annulus(0.0, a).unwrap();
            let scale = whole.hertz_force_n - whole.adhesive_force_n;
            close(whole.normal_force_n, load, scale, 1e-13);
            assert_eq!(a, response.contact_radius_m);
            assert_eq!(profile.normal_force_n(), response.normal_force_n);
            if fraction == 1.0 {
                assert!(response.marginal);
            }
        }
        assert!(matches!(
            AdhesiveTractionProfile::jkr_mixed_mode_force(&source, load, [critical.next_up(), 0.0]),
            Err(AdhesiveTractionError::MixedMode(_))
        ));
    }
    assert!(matches!(
        AdhesiveTractionProfile::jkr_mixed_mode_force(&source, pull_off.next_down(), [0.0; 2]),
        Err(AdhesiveTractionError::MixedMode(_))
    ));
}

#[test]
fn mixed_mode_spatial_transfer_recovers_original_force_and_affine_work() {
    use fs_contact::normal_patch::traction::sectors::RadialMomentBudget;
    use fs_exec::CancelGate;
    let source = mixed_jkr();
    let delta = 1e-6;
    let critical = source.critical_tangential_displacement_m(delta).unwrap();
    let profile =
        AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, [0.8 * critical, 0.0])
            .unwrap();
    let a = profile.contact_radius_m();
    let whole = profile.annulus(0.0, a).unwrap();
    let scale = whole.hertz_force_n - whole.adhesive_force_n;
    let radii = [0.0, 0.31 * a, 0.8 * a, a.next_down(), a, 1.2 * a];
    let mut forces = [f64::NAN; 5];
    profile
        .transfer_annular_forces(&radii, &mut forces, 1e-12 * scale)
        .unwrap();
    close(forces.iter().sum(), profile.normal_force_n(), scale, 1e-12);
    assert!(forces[3] < 0.0);
    assert_eq!(forces[4], 0.0);

    let exact = RadialMomentBudget {
        absolute_tolerance_n_m: 0.0,
        relative_tolerance: 0.0,
        maximum_evaluations: 0,
    };
    let angular_width = 0.83;
    let start_angle = 0.27_f64;
    let end_angle = start_angle + angular_width;
    let sector = profile
        .sector(0.0, a, start_angle, end_angle, exact, &CancelGate::new())
        .unwrap();
    // Exact radial first moments for Hertz and flat-punch pressures, before
    // the azimuth integral: pi*a^3*(peak/16 - q/4).
    let e = source
        .normal_source()
        .normal_state()
        .elastic()
        .reduced_modulus_pa();
    let r = source.normal_source().effective_radius_m();
    let peak = 2.0 * e * a / (PI * r);
    let q = e * (a * a / r - delta) / (PI * a);
    let radial = PI * a.powi(3) * (peak / 16.0 - q / 4.0);
    let first = [
        radial * (end_angle.sin() - start_angle.sin()),
        radial * (start_angle.cos() - end_angle.cos()),
    ];
    for (actual, expected) in sector.first_moment_n_m.into_iter().zip(first) {
        close(actual, expected, scale * a, 1e-13);
    }
    let (power, error) = sector.affine_normal_power(0.37, [2.0, -3.0]).unwrap();
    let expected = 0.37 * profile.normal_force_n() * angular_width / (2.0 * PI) + 2.0 * first[0]
        - 3.0 * first[1];
    close(power, expected, scale * (0.37 + 5.0 * a), 1e-13);
    assert_eq!(error, 0.0);
}

#[test]
fn mixed_mode_spatial_zero_shear_rotation_and_source_refusals() {
    let source = mixed_jkr();
    let delta = 1e-6;
    let pure = AdhesiveTractionProfile::jkr_attached(source.normal_source(), delta).unwrap();
    let zero =
        AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, [0.0; 2]).unwrap();
    assert_eq!(zero.contact_radius_m(), pure.contact_radius_m());
    for fraction in [0.0, 0.3, 0.9] {
        let radius = fraction * zero.contact_radius_m();
        let a = zero.pressure_pa(radius).unwrap();
        let b = pure.pressure_pa(radius).unwrap();
        close(a, b, a.abs().max(b.abs()), 1e-12);
    }
    let critical = source.critical_tangential_displacement_m(delta).unwrap();
    let displacement = 0.7 * critical;
    let first =
        AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, [displacement, 0.0])
            .unwrap();
    let rotated =
        AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, [0.0, -displacement])
            .unwrap();
    assert_eq!(first, rotated);
    let fold =
        AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, [critical, 0.0])
            .unwrap();
    assert!(fold.contact_radius_m() < first.contact_radius_m());
    for shift in [[critical.next_up(), 0.0], [f64::NAN, 0.0]] {
        assert!(matches!(
            AdhesiveTractionProfile::jkr_mixed_mode_displacement(&source, delta, shift),
            Err(AdhesiveTractionError::MixedMode(_))
        ));
    }
}

fn maugis(lambda: f64) -> SourceMaugisSphereContact {
    let a = common::solid("A", 1000.0, 1e11);
    let b = common::solid("B", 2000.0, 2e11);
    let normal = common::bind(&a, &b, None, 1e-7).normal_state().clone();
    let r = 0.01;
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

fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; error/scale={:e}",
        (a - b).abs() / scale
    );
}

#[test]
fn jkr_retains_tensile_rim_and_integrates_its_singularity_without_clipping() {
    let source = jkr();
    for delta in [source.minimum_attached_indentation_m(), 0.0, 1e-6, 8e-6] {
        let response = source.evaluate_attached(delta).unwrap();
        let profile = AdhesiveTractionProfile::jkr_attached(&source, delta).unwrap();
        let a = response.contact_radius_m;
        let e = source.normal_state().elastic().reduced_modulus_pa();
        let w = source.normal_state().adhesion_energy_j_per_m2();
        let peak = 2.0 * e * a / (PI * source.effective_radius_m());
        let tensile = (2.0 * e * w / (PI * a)).sqrt();
        for fraction in [0.0_f64, 0.3, 0.8, 0.99] {
            let s = (1.0 - fraction * fraction).sqrt();
            close(
                profile.pressure_pa(fraction * a).unwrap(),
                peak * s - tensile / s,
                peak + tensile / s,
                5e-14,
            );
        }
        assert_eq!(
            profile.pressure_pa(a),
            Err(AdhesiveTractionError::SingularEdge)
        );
        assert_eq!(profile.pressure_pa(a.next_up()).unwrap(), 0.0);
        let whole = profile.annulus(0.0, a).unwrap();
        close(
            whole.normal_force_n,
            response.normal_force_n,
            whole.hertz_force_n - whole.adhesive_force_n,
            5e-14,
        );
        close(
            whole.solid_contact_area_m2,
            response.contact_area_m2,
            response.contact_area_m2,
            5e-15,
        );
        let rim = profile.annulus((1.0 - 1e-10) * a, a).unwrap();
        assert!(rim.mean_pressure_pa.is_finite() && rim.mean_pressure_pa < 0.0);
        assert!(rim.adhesive_force_n < 0.0);
        let center = profile.annulus(0.0, a * 1e-10).unwrap();
        close(
            center.mean_pressure_pa,
            peak - tensile,
            peak + tensile,
            5e-15,
        );
        assert_eq!(profile.model_identity(), source.identity());
        assert_eq!(
            profile.normal_state_identity(),
            source.normal_state().identity()
        );
        assert_eq!(profile.authority(), InputAuthority::SyntheticFixture);
    }
    let profile = AdhesiveTractionProfile::jkr_attached(&source, 0.0).unwrap();
    assert!(profile.normal_force_n() < 0.0);
    assert!(
        profile
            .annulus(0.0, 0.5 * profile.contact_radius_m())
            .unwrap()
            .normal_force_n
            > 0.0
    );
}

#[test]
fn annular_transfer_preserves_force_and_virtual_work_under_refinement() {
    let profile = AdhesiveTractionProfile::jkr_attached(&jkr(), 0.0).unwrap();
    let a = profile.contact_radius_m();
    let scale = profile.annulus(0.0, a).unwrap().hertz_force_n;
    for count in [2, 17, 128] {
        let edges: Vec<_> = (0..=count).map(|i| a * (i as f64 / count as f64)).collect();
        let mut forces = vec![123.0; count];
        let residual = profile
            .transfer_annular_forces(&edges, &mut forces, 1e-12 * scale)
            .unwrap();
        close(residual, 0.0, scale, 1e-12);
        let velocity = -0.3;
        let power: f64 = forces.iter().map(|f| f * velocity).sum();
        close(
            power,
            profile.normal_force_n() * velocity,
            scale * velocity.abs(),
            1e-12,
        );
    }
    let whole = profile.annulus(0.0, a).unwrap();
    let left = profile.annulus(0.0, 0.83 * a).unwrap();
    let right = profile.annulus(0.83 * a, 1.4 * a).unwrap();
    close(
        left.normal_force_n + right.normal_force_n,
        whole.normal_force_n,
        scale,
        5e-14,
    );
    close(
        right.mean_pressure_pa * right.area_m2,
        right.normal_force_n,
        scale,
        5e-14,
    );
}

#[test]
fn native_source_binding_preserves_its_original_branch_admission() {
    use fs_contact::normal_patch::NormalLawSource;
    let original = jkr();
    let binding: NormalLawSource = original.clone().into();
    let profile = binding.adhesive_traction_profile(0.0).unwrap().unwrap();
    assert_eq!(profile.model_identity(), original.identity());
    assert!(
        binding
            .adhesive_traction_profile(original.minimum_attached_indentation_m().next_down())
            .is_err()
    );
    let source = maugis(0.2);
    let limits = source.noncontact_limits().unwrap().unwrap();
    let d = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let potential = mechanical::MaugisAttachedPotential::new(
        source,
        continuation::MaugisContinuationBudget {
            work: MaugisWorkBudget {
                absolute_tolerance_j: 1e-25,
                relative_tolerance: 1e-10,
                maximum_evaluations: 10000,
            },
            energy_absolute_tolerance_j: 1e-22,
            energy_relative_tolerance: 1e-9,
        },
    )
    .unwrap();
    let restricted: NormalLawSource = potential.clone().into();
    assert!(restricted.adhesive_traction_profile(d).is_err());
    let extended: NormalLawSource = potential.with_noncontact_continuation().unwrap().into();
    let profile = extended.adhesive_traction_profile(d).unwrap().unwrap();
    assert_eq!(profile.contact_radius_m(), 0.0);
    assert!(profile.normal_force_n() < 0.0);
}

#[test]
fn maugis_patch_loads_match_original_pressure_and_resultant_across_tabor_regimes() {
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        let source = maugis(lambda);
        for delta in [
            source.minimum_attached_indentation_m(),
            0.0,
            0.5 * source.envelope().maximum_indentation_m,
        ] {
            let equilibrium = source.evaluate_contact_energy(delta).unwrap();
            let response = equilibrium.contact_response().unwrap();
            let profile = AdhesiveTractionProfile::maugis(&source, &equilibrium).unwrap();
            let a = profile.contact_radius_m();
            let c = profile.interaction_radius_m();
            let sigma =
                source.normal_state().adhesion_energy_j_per_m2() / source.interaction_range_m();
            for fraction in [0.0, 0.2, 0.8, 1.0] {
                close(
                    profile.pressure_pa(fraction * a).unwrap(),
                    response.pressure_pa(fraction * a).unwrap(),
                    sigma + response.repulsive_normal_force_n / (PI * c * c),
                    2e-13,
                );
            }
            assert_eq!(profile.pressure_pa(c).unwrap(), -sigma);
            let whole = profile.annulus(0.0, c).unwrap();
            let scale = whole.hertz_force_n - whole.adhesive_force_n;
            close(
                whole.normal_force_n,
                equilibrium.normal_force_n(),
                scale,
                2e-12,
            );
            let edges: Vec<_> = (0..=97).map(|i| c * (f64::from(i) / 97.0)).collect();
            let mut forces = vec![0.0; 97];
            profile
                .transfer_annular_forces(&edges, &mut forces, 3e-12 * scale)
                .unwrap();
            if a > 0.0 {
                // A narrow interior ring converges to point pressure even when
                // the cohesive annulus is very thin or very wide.
                let r = 0.4 * a;
                let ring = profile.annulus(r, r + 1e-8 * a).unwrap();
                close(
                    ring.mean_pressure_pa,
                    response.pressure_pa(r + 0.5e-8 * a).unwrap(),
                    sigma + response.repulsive_normal_force_n / (PI * a * a),
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn noncontact_maugis_attraction_does_not_fabricate_solid_area() {
    let source = maugis(0.2);
    let limits = source.noncontact_limits().unwrap().unwrap();
    let delta = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let equilibrium = source.evaluate_noncontact(delta).unwrap();
    let profile = AdhesiveTractionProfile::maugis(&source, &equilibrium).unwrap();
    let c = profile.interaction_radius_m();
    let ring = profile.annulus(0.0, c).unwrap();
    assert_eq!(ring.solid_contact_area_m2, 0.0);
    assert_eq!(ring.hertz_force_n, 0.0);
    assert!(ring.normal_force_n < 0.0);
    close(
        ring.normal_force_n,
        equilibrium.normal_force_n(),
        equilibrium.normal_force_n().abs(),
        5e-15,
    );
    let detached = source
        .evaluate_detached(-2.0 * source.interaction_range_m())
        .unwrap();
    let zero = AdhesiveTractionProfile::maugis(&source, &detached).unwrap();
    assert_eq!(zero.pressure_pa(0.0).unwrap(), 0.0);
    assert_eq!(zero.annulus(0.0, c).unwrap().normal_force_n, 0.0);
}

#[test]
fn source_domain_and_partition_failures_leave_output_untouched() {
    let source = jkr();
    assert!(
        AdhesiveTractionProfile::jkr_attached(
            &source,
            source.minimum_attached_indentation_m().next_down()
        )
        .is_err()
    );
    let profile = AdhesiveTractionProfile::jkr_attached(&source, 0.0).unwrap();
    let a = profile.contact_radius_m();
    for edges in [
        [0.0, 0.5 * a, 0.9 * a],
        [0.0, a, 0.5 * a],
        [0.0, 0.5 * a, f64::NAN],
    ] {
        let mut output = [123.0, 456.0];
        assert!(
            profile
                .transfer_annular_forces(&edges, &mut output, 1e-12)
                .is_err()
        );
        assert_eq!(output, [123.0, 456.0]);
    }
    let mut output = [123.0, 456.0];
    assert!(
        profile
            .transfer_annular_forces(&[0.0, a, 2.0 * a], &mut output, -1.0)
            .is_err()
    );
    assert_eq!(output, [123.0, 456.0]);
    assert!(profile.pressure_pa(f64::NAN).is_err());
    assert!(profile.annulus(a, a).is_err());
    let first = maugis(0.2);
    let other = maugis(0.3);
    let equilibrium = first.evaluate_contact_energy(0.0).unwrap();
    assert_eq!(
        AdhesiveTractionProfile::maugis(&other, &equilibrium),
        Err(AdhesiveTractionError::SourceMismatch)
    );
}
