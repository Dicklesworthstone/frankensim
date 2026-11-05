//! Original-source force, torque and affine-work transfer into finite patches.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_contact::normal_patch::{
    adhesive::*,
    maugis::*,
    traction::{sectors::*, *},
};
use fs_exec::CancelGate;
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

fn exact() -> RadialMomentBudget {
    RadialMomentBudget {
        absolute_tolerance_n_m: 0.0,
        relative_tolerance: 0.0,
        maximum_evaluations: 0,
    }
}
fn budget(force: f64, radius: f64) -> PolarTransferBudget {
    PolarTransferBudget {
        radial_moment: RadialMomentBudget {
            absolute_tolerance_n_m: 1e-12 * force * radius,
            relative_tolerance: 1e-10,
            maximum_evaluations: MAX_RADIAL_MOMENT_EVALUATIONS,
        },
        force_absolute_tolerance_n: 1e-11 * force,
        moment_absolute_tolerance_n_m: 1e-11 * force * radius,
    }
}
fn close(x: f64, y: f64, scale: f64, tolerance: f64) {
    assert!((x - y).abs() <= tolerance * scale, "{x:e} != {y:e}");
}

#[test]
fn zero_resultant_jkr_sector_retains_its_couple_and_rotational_work() {
    let source = jkr();
    let e = source.normal_state().elastic().reduced_modulus_pa();
    let w = source.normal_state().adhesion_energy_j_per_m2();
    let r = source.effective_radius_m();
    let a0 = (4.5 * PI * w * r * r / e).cbrt();
    let profile = AdhesiveTractionProfile::jkr_attached(&source, a0 * a0 / (3.0 * r)).unwrap();
    let a = profile.contact_radius_m();
    let whole = profile.annulus(0.0, a).unwrap();
    close(whole.normal_force_n, 0.0, whole.hertz_force_n, 2e-11);
    let cell = profile
        .sector(0.0, a, 0.0, FRAC_PI_2, exact(), &CancelGate::new())
        .unwrap();
    // Full-disk Hertz and tensile centroids are 3πa/16 and πa/4.
    // At zero resultant their difference still produces a finite sector couple.
    let expected = -whole.hertz_force_n * a / 32.0;
    for moment in cell.first_moment_n_m {
        close(moment, expected, expected.abs(), 1e-10);
    }
    assert_eq!(cell.moment_error_estimate_n_m, 0.0);
    let (force, torque) = cell.wrench();
    assert_eq!(force, [0.0, 0.0, cell.normal_force_n]);
    assert_eq!(
        torque,
        [cell.first_moment_n_m[1], -cell.first_moment_n_m[0], 0.0]
    );
    let (power, error) = cell.affine_normal_power(0.0, [-3.0, 7.0]).unwrap();
    close(
        power,
        torque[0] * 7.0 + torque[1] * 3.0,
        expected.abs() * 10.0,
        5e-15,
    );
    assert!(power != 0.0);
    assert_eq!(error, 0.0);
}

#[test]
fn polar_refinement_conserves_source_force_moment_and_affine_power() {
    let source = jkr();
    let profile = AdhesiveTractionProfile::jkr_attached(&source, 0.0).unwrap();
    let a = profile.contact_radius_m();
    let whole = profile.annulus(0.0, a).unwrap();
    let scale = whole.hertz_force_n - whole.adhesive_force_n;
    let gate = CancelGate::new();
    for (rings, sectors) in [(2, 3), (9, 17), (16, 64)] {
        let radii: Vec<_> = (0..=rings)
            .map(|i| 1.2 * a * i as f64 / rings as f64)
            .collect();
        let angles: Vec<_> = (0..=sectors)
            .map(|i| TAU * (i as f64 / sectors as f64))
            .collect();
        let mut cells = vec![SectorAdhesiveTraction::default(); rings * sectors];
        let mut scratch = vec![PolarRingScratch::default(); rings];
        let report = profile
            .transfer_polar_grid(
                &radii,
                &angles,
                &mut cells,
                &mut scratch,
                budget(scale, a),
                &gate,
            )
            .unwrap();
        assert_eq!(report.model_identity, source.identity());
        assert_eq!(
            report.normal_state_identity,
            source.normal_state().identity()
        );
        assert_eq!(report.authority, InputAuthority::SyntheticFixture);
        assert_eq!(report.radial_evaluations, 0);
        let power: f64 = cells
            .iter()
            .map(|cell| cell.affine_normal_power(-0.3, [7000.0, -9000.0]).unwrap().0)
            .sum();
        close(
            power,
            -0.3 * profile.normal_force_n(),
            scale * (0.3 + 16000.0 * a),
            2e-12,
        );
        close(
            cells.iter().map(|c| c.solid_contact_area_m2).sum(),
            whole.solid_contact_area_m2,
            whole.solid_contact_area_m2,
            2e-12,
        );
        for component in 0..2 {
            close(
                cells.iter().map(|c| c.first_moment_n_m[component]).sum(),
                0.0,
                scale * a,
                2e-12,
            );
        }
    }
}

#[test]
fn maugis_radial_moments_add_and_one_radial_solve_serves_all_angles() {
    let gate = CancelGate::new();
    for lambda in [0.02, 0.2, 1.0, 5.0, 1000.0] {
        let source = maugis(lambda);
        let state = source.evaluate_contact_energy(0.0).unwrap();
        let profile = AdhesiveTractionProfile::maugis(&source, &state).unwrap();
        let a = profile.contact_radius_m();
        let c = profile.interaction_radius_m();
        let whole = profile.annulus(0.0, c).unwrap();
        let scale = whole.hertz_force_n - whole.adhesive_force_n;
        let limits = budget(scale, c);
        let quarter = profile
            .sector(0.0, c, 0.0, FRAC_PI_2, limits.radial_moment, &gate)
            .unwrap();
        let mut radii = vec![0.0, 0.45 * a, a];
        if c > a {
            radii.push(c);
        }
        let mut sum = [0.0; 2];
        for pair in radii.windows(2) {
            let cell = profile
                .sector(
                    pair[0],
                    pair[1],
                    0.0,
                    FRAC_PI_2,
                    limits.radial_moment,
                    &gate,
                )
                .unwrap();
            for k in 0..2 {
                sum[k] += cell.first_moment_n_m[k];
            }
        }
        for k in 0..2 {
            close(sum[k], quarter.first_moment_n_m[k], scale * c, 2e-10);
        }
        let mut previous_evaluations = None;
        for sectors in [4, 11] {
            let angles: Vec<_> = (0..=sectors)
                .map(|i| TAU * (i as f64 / sectors as f64))
                .collect();
            let mut cells = vec![SectorAdhesiveTraction::default(); (radii.len() - 1) * sectors];
            let mut scratch = vec![PolarRingScratch::default(); radii.len() - 1];
            let report = profile
                .transfer_polar_grid(&radii, &angles, &mut cells, &mut scratch, limits, &gate)
                .unwrap();
            assert!(report.radial_evaluations > 0);
            if let Some(previous) = previous_evaluations {
                assert_eq!(report.radial_evaluations, previous);
            }
            previous_evaluations = Some(report.radial_evaluations);
        }
        // Complete circles have identically zero first moment: no numerical
        // centroid should be requested, even with an exact zero-work budget.
        let circle = profile.sector(0.0, c, 0.0, TAU, exact(), &gate).unwrap();
        assert_eq!(circle.normal_force_n, whole.normal_force_n);
        assert_eq!(circle.first_moment_n_m, [0.0, 0.0]);
    }
}

#[test]
fn noncontact_disk_has_exact_wrench_without_inventing_solid_contact() {
    let source = maugis(0.2);
    let limits = source.noncontact_limits().unwrap().unwrap();
    let delta = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let state = source.evaluate_noncontact(delta).unwrap();
    let profile = AdhesiveTractionProfile::maugis(&source, &state).unwrap();
    let c = profile.interaction_radius_m();
    let quarter = profile
        .sector(0.0, c, 0.0, FRAC_PI_2, exact(), &CancelGate::new())
        .unwrap();
    assert_eq!(quarter.solid_contact_area_m2, 0.0);
    assert_eq!(quarter.hertz_force_n, 0.0);
    assert!(quarter.normal_force_n < 0.0);
    for moment in quarter.first_moment_n_m {
        close(
            moment,
            state.normal_force_n() * c / (3.0 * PI),
            state.normal_force_n().abs() * c,
            5e-15,
        );
    }
    let detached = source
        .evaluate_detached(-2.0 * source.interaction_range_m())
        .unwrap();
    let zero = AdhesiveTractionProfile::maugis(&source, &detached)
        .unwrap()
        .sector(0.0, c, 0.0, FRAC_PI_2, exact(), &CancelGate::new())
        .unwrap();
    assert_eq!(zero.normal_force_n, 0.0);
    assert_eq!(zero.first_moment_n_m, [0.0, 0.0]);
}

#[test]
fn numerical_budget_invalid_partition_and_cancellation_refuse_atomically() {
    let source = maugis(0.2);
    let state = source.evaluate_contact_energy(0.0).unwrap();
    let profile = AdhesiveTractionProfile::maugis(&source, &state).unwrap();
    let c = profile.interaction_radius_m();
    let scale = profile.annulus(0.0, c).unwrap().adhesive_force_n.abs();
    let sentinel = SectorAdhesiveTraction {
        normal_force_n: 123.0,
        first_moment_n_m: [456.0, 789.0],
        ..SectorAdhesiveTraction::default()
    };
    let mut output = [sentinel; 4];
    let mut scratch = [PolarRingScratch::default(); 2];
    let mut limits = budget(scale, c);
    limits.radial_moment.maximum_evaluations = 27;
    assert!(matches!(
        profile.transfer_polar_grid(
            &[0.0, 0.5 * c, c],
            &[0.0, PI, TAU],
            &mut output,
            &mut scratch,
            limits,
            &CancelGate::new()
        ),
        Err(SectorTractionError::MomentBudget { .. })
    ));
    assert_eq!(output, [sentinel; 4]);
    for radii in [
        [0.0, 0.5 * c, 0.9 * c],
        [0.0, 1.2 * c, c],
        [0.0, 0.5 * c, f64::NAN],
    ] {
        assert!(
            profile
                .transfer_polar_grid(
                    &radii,
                    &[0.0, PI, TAU],
                    &mut output,
                    &mut scratch,
                    budget(scale, c),
                    &CancelGate::new()
                )
                .is_err()
        );
        assert_eq!(output, [sentinel; 4]);
    }
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        profile.transfer_polar_grid(
            &[0.0, 0.5 * c, c],
            &[0.0, PI, TAU],
            &mut output,
            &mut scratch,
            budget(scale, c),
            &cancelled
        ),
        Err(SectorTractionError::Cancelled)
    );
    assert_eq!(output, [sentinel; 4]);
    let mut circle = [sentinel];
    let mut radial = [PolarRingScratch::default()];
    limits.radial_moment = exact();
    let report = profile
        .transfer_polar_grid(
            &[0.0, c],
            &[0.0, TAU],
            &mut circle,
            &mut radial,
            limits,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(report.radial_evaluations, 0);
    assert_eq!(circle[0].first_moment_n_m, [0.0, 0.0]);
}
