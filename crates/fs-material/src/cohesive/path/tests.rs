use super::*;
use crate::cohesive::{BilinearCohesive, CohesiveViscosity, mixed_mode::MixedModeParameters};

fn common() -> BilinearCohesive {
    BilinearCohesive::new(1e9, 2e5, 100.0, 2e9, 0.02).unwrap()
}

fn mixed() -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e8,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 1.6e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 400.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

fn close_vector(a: [f64; 3], b: [f64; 3], relative: f64) {
    let error = norm(core::array::from_fn(|i| a[i] - b[i])).unwrap();
    assert!(
        error <= relative * (1.0 + norm(a).unwrap().max(norm(b).unwrap())),
        "{a:?} != {b:?}; error={error:e}"
    );
}

fn at(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    core::array::from_fn(|i| (1.0 - s) * a[i] + s * b[i])
}

#[test]
fn elastic_mean_integrates_the_normal_switch_and_keeps_endpoint_tangents_separate() {
    let law = common();
    let old = law
        .bonded_point()
        .trial([-1e-4, 3e-5, -2e-5], 1e-9)
        .unwrap()
        .point;
    let end = [1e-4, -1e-5, 6e-5];
    let (step, mean) = old.trial_with_mean_traction(end, 1e-9).unwrap();
    close_vector(mean, [-2.5e4, 1e4, 2e4], 1e-14);
    assert_eq!(step, old.trial(end, 1e-9).unwrap());
    assert_ne!(mean, step.traction_pa);
    close(step.energy.work, -3.8, 1e-13);
    let (_, resting) = step.point.trial_with_mean_traction(end, 1e-9).unwrap();
    assert_eq!(resting, step.traction_pa);
}

#[test]
fn vector_mean_matches_independent_traction_quadrature_through_all_events() {
    for law in [common(), mixed()] {
        for (start, end) in [
            ([-1.5e-4, 3e-5, 4e-5], [5e-4, 1.4e-3, -7e-4]),
            ([3.5e-4, 2.5e-4, 1e-4], [-2e-4, 1.4e-3, -5e-4]),
            ([3e-4, 2.5e-4, 0.0], [1e-4, -1.8e-3, 8e-4]),
            ([4e-4, 2e-4, 1e-4], [8e-5, -4e-5, 2e-5]),
            ([6e-5, -3e-5, 2e-5], [-5e-4, 0.012, 0.004]),
            ([0.0, 0.014, 0.0], [-4e-4, -0.002, 0.001]),
        ] {
            let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
            let (step, mean) = old.trial_with_mean_traction(end, 1e-8).unwrap();
            let mut quadrature = [0.0; 3];
            let count = 32_768;
            for i in 0..count {
                let s = (f64::from(i) + 0.5) / f64::from(count);
                let traction = old.trial(at(start, end, s), 1e-8).unwrap().traction_pa;
                for j in 0..3 {
                    quadrature[j] += traction[j] / f64::from(count);
                }
            }
            close_vector(mean, quadrature, 2e-7);
            close(
                dot(mean, core::array::from_fn(|i| end[i] - start[i])),
                step.energy.work,
                2e-12,
            );
            assert_eq!(step, old.trial(end, 1e-8).unwrap());
        }
    }
}

#[test]
fn rotated_shear_paths_rotate_the_vector_mean_without_reselecting_toughness() {
    let rotate = |g: [f64; 3]| [g[0], 0.6 * g[1] - 0.8 * g[2], 0.8 * g[1] + 0.6 * g[2]];
    for law in [common(), mixed()] {
        let mut point = law.bonded_point();
        let mut rotated = law.bonded_point();
        for end in [
            [7e-5, -2e-5, 1e-5],
            [3e-4, 6e-4, -2e-4],
            [-1e-4, -8e-4, 5e-4],
            [9e-5, 3e-5, -4e-5],
            [7e-4, 0.012, 0.004],
            [-3e-4, -2e-4, 1e-4],
        ] {
            let (a, mean) = point.trial_with_mean_traction(end, 1e-8).unwrap();
            let (b, rotated_mean) = rotated.trial_with_mean_traction(rotate(end), 1e-8).unwrap();
            close_vector(rotate(mean), rotated_mean, 3e-12);
            close(a.energy.work, b.energy.work, 3e-12);
            close(a.point.damage(), b.point.damage(), 3e-12);
            point = a.point;
            rotated = b.point;
        }
        assert!(point.failed());
        assert!(rotated.failed());
    }
}

#[test]
fn nonproportional_means_and_work_are_additive_under_path_subdivision() {
    for law in [common(), mixed()] {
        for (start, end) in [
            ([7e-5, 3e-5, 0.0], [7e-4, 1.2e-3, -4e-4]),
            ([4e-4, 2e-4, 1e-4], [-4e-4, -0.003, 0.001]),
            ([-3e-4, 0.0, 0.0], [5e-4, 0.012, -0.004]),
        ] {
            let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
            let (whole, mean) = old.trial_with_mean_traction(end, 1e-8).unwrap();
            let mut point = old;
            let mut divided_mean = [0.0; 3];
            let mut work = 0.0;
            for i in 1..=64 {
                let (step, part_mean) = point
                    .trial_with_mean_traction(at(start, end, f64::from(i) / 64.0), 1e-8)
                    .unwrap();
                for j in 0..3 {
                    divided_mean[j] += part_mean[j] / 64.0;
                }
                work += step.energy.work;
                point = step.point;
            }
            close_vector(mean, divided_mean, 5e-12);
            close(whole.energy.work, work, 5e-12);
            close(whole.point.damage(), point.damage(), 5e-12);
            if let (Some(a), Some(b)) = (whole.point.mixed_history, point.mixed_history) {
                close(a.energy_shear_fraction, b.energy_shear_fraction, 5e-12);
            }
        }
    }
}

#[test]
fn complete_rupture_retains_work_when_the_initiation_fraction_is_below_machine_resolution() {
    for stiffness in [1e9, 1e24] {
        let law = BilinearCohesive::new(stiffness, 1.0, 1e-3, stiffness, 0.01).unwrap();
        let direction = [0.6, 0.48, -0.64];
        let end = direction.map(|x| x * 0.004);
        let (step, mean) = law
            .bonded_point()
            .trial_with_mean_traction(end, 1e-14)
            .unwrap();
        assert!(step.point.failed());
        assert_eq!(step.traction_pa, [0.0; 3]);
        close_vector(mean, direction.map(|x| 0.25 * x), 1e-13);
        close(dot(mean, end), 1e-3, 1e-14);
    }
}

#[test]
fn short_softening_intervals_retain_traction_close_to_complete_failure() {
    let law = common();
    let start = [law.failure - 1e-12, 0.0, 0.0];
    let end = [law.failure - 5e-13, 0.0, 0.0];
    let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
    let (step, mean) = old.trial_with_mean_traction(end, 1e-8).unwrap();
    let expected = 0.5 * old.traction_pa()[0] + 0.5 * step.traction_pa[0];
    close(mean[0], expected, 2e-11);
    assert!(mean[0] > 0.0);
    assert_eq!([mean[1], mean[2]], [0.0; 2]);
}

#[test]
fn timed_means_include_actual_bond_viscosity_and_preserve_subdivided_work() {
    let viscosity = CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap();
    for law in [
        common().with_viscosity(viscosity),
        mixed().with_viscosity(viscosity),
    ] {
        let start = [3e-4, 2e-4, 1e-4];
        let end = [-4e-4, -0.006, 0.002];
        let dt = 0.2;
        let old = law
            .bonded_point()
            .trial_over_time(start, dt, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        let (whole, mean) = old
            .trial_over_time_with_mean_traction(end, dt, 1e-8)
            .unwrap();
        assert!(whole.viscous_dissipated_j_m2 > 0.0);
        close(
            dot(mean, core::array::from_fn(|i| end[i] - start[i])),
            whole.total_work_j_m2(),
            3e-12,
        );
        let mut quadrature = [0.0; 3];
        for i in 0..16_384 {
            let s = (f64::from(i) + 0.5) / 16_384.0;
            let traction = old
                .trial_over_time(at(start, end, s), dt * s, 1e-8)
                .unwrap()
                .traction_pa;
            for j in 0..3 {
                quadrature[j] += traction[j] / 16_384.0;
            }
        }
        close_vector(mean, quadrature, 3e-7);
        let mut point = old;
        let mut divided = [0.0; 3];
        let mut work = 0.0;
        let mut viscous = 0.0;
        for i in 1..=64 {
            let (step, part) = point
                .trial_over_time_with_mean_traction(
                    at(start, end, f64::from(i) / 64.0),
                    dt / 64.0,
                    1e-8,
                )
                .unwrap();
            for j in 0..3 {
                divided[j] += part[j] / 64.0;
            }
            work += step.total_work_j_m2();
            viscous += step.viscous_dissipated_j_m2;
            point = step.equilibrium.point;
        }
        close_vector(mean, divided, 5e-12);
        close(work, whole.total_work_j_m2(), 5e-12);
        close(viscous, whole.viscous_dissipated_j_m2, 5e-12);
        assert!(point.failed());
    }
}

#[test]
fn mean_traction_refusals_leave_accepted_history_unchanged() {
    let law = mixed().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 0.1).unwrap());
    let point = law.bonded_point();
    let saved = point;
    assert!(
        point
            .trial_with_mean_traction([1e-4, 0.0, 0.0], 1e-8)
            .is_err()
    );
    for (jump, dt) in [
        ([1e-4, 0.0, 0.0], 0.0),
        ([0.01, 0.0, 0.0], 0.001),
        ([0.03, 0.0, 0.0], 1.0),
        ([f64::NAN, 0.0, 0.0], 1.0),
    ] {
        assert!(
            point
                .trial_over_time_with_mean_traction(jump, dt, 1e-8)
                .is_err()
        );
    }
    assert_eq!(point, saved);
}


// G0/G1: independently perturb the EXISTING physical mean-traction owner,
// including its retained history and exact viscosity loss. Endpoint tangents
// cannot satisfy this check through softening or a changing path direction.
fn check_path_tangent(
    old: CohesivePoint<'_>,
    end: [f64; 3],
    dt: f64,
    perturbation: f64,
) -> [[f64; 3]; 3] {
    let saved = old;
    let reference = old.trial_over_time_with_mean_traction(end, dt, 1e-8).unwrap();
    let (step, mean, tangent) = old
        .trial_over_time_with_mean_traction_tangent(end, dt, 1e-8)
        .unwrap();
    assert_eq!(step, reference.0);
    assert_eq!(mean, reference.1);
    let tangent = tangent.expect("smooth covered path must provide its derivative");
    let scale = tangent.iter().flatten().fold(1.0_f64, |a, b| a.max(b.abs()));
    for j in 0..3 {
        let mut plus = end;
        let mut minus = end;
        plus[j] += perturbation;
        minus[j] -= perturbation;
        let a = old
            .trial_over_time_with_mean_traction(plus, dt, 1e-8)
            .unwrap();
        let b = old
            .trial_over_time_with_mean_traction(minus, dt, 1e-8)
            .unwrap();
        for i in 0..3 {
            let difference = (a.1[i] - b.1[i]) / (2.0 * perturbation);
            assert!(
                (difference - tangent[i][j]).abs() <= 3e-7 * scale,
                "row {i}, column {j}: tangent={} difference={difference}, old={:?}, end={end:?}",
                tangent[i][j],
                old.jump_m(),
            );
        }
        // Derivative of actual path work uses BOTH the varying mean and the
        // varying endpoint jump; this also exercises the viscous path integral.
        let difference = (a.0.total_work_j_m2() - b.0.total_work_j_m2())
            / (2.0 * perturbation);
        let predicted = mean[j]
            + (0..3).map(|i| tangent[i][j] * (end[i] - old.jump[i])).sum::<f64>();
        let work_scale = 1.0 + norm(mean).unwrap()
            + scale * norm(core::array::from_fn(|i| end[i] - old.jump[i])).unwrap();
        assert!((difference - predicted).abs() <= 3e-7 * work_scale);
    }
    assert_eq!(old, saved);
    tangent
}

#[test]
fn path_tangents_follow_damage_reversal_rupture_and_viscous_integrity() {
    let viscosity = CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap();
    for law in [
        common(),
        common().with_viscosity(viscosity),
        mixed(),
        mixed().with_viscosity(viscosity),
    ] {
        for (start, end) in [
            ([3.5e-4, 2.5e-4, 1e-4], [-2e-4, 1.4e-3, -5e-4]),
            ([3e-4, 2.5e-4, 0.0], [1e-4, -1.8e-3, 8e-4]),
            ([4e-4, 2e-4, 1e-4], [8e-5, -4e-5, 2e-5]),
            ([4e-4, 2e-4, 1e-4], [4.01e-4, 2.002e-4, 1.003e-4]),
            ([3e-4, 4e-4, -1e-4], [7e-4, 0.012, 0.004]),
            ([0.0, 0.014, 0.0], [-4e-4, -0.002, 0.001]),
        ] {
            let old = law.bonded_point()
                .trial_over_time(start, 0.2, 1e-8).unwrap().equilibrium.point;
            let tangent = check_path_tangent(old, end, 0.2, 1e-9);
            if law.viscosity().is_none() {
                // Once the envelope is frozen, the intrinsic path force has a
                // symmetric derivative even during irreversible radial growth.
                let scale = tangent.iter().flatten().fold(1.0_f64, |a, b| a.max(b.abs()));
                for i in 0..3 {
                    for j in 0..3 {
                        assert!((tangent[i][j] - tangent[j][i]).abs() < 2e-12 * scale);
                    }
                }
            }
        }
    }
    let law = common();
    for (start, end) in [
        // Damage initiation exactly at the normal sign switch; the ensuing
        // softening segment is in compression, including its derivative.
        ([1e-4, 0.0, 0.0], [-1e-4, 4e-4, 0.0]),
        ([-1.5e-4, 3e-5, 4e-5], [5e-4, 1.4e-3, -7e-4]),
        ([6e-5, -3e-5, 2e-5], [-5e-4, 0.012, 0.004]),
    ] {
        let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
        check_path_tangent(old, end, 0.2, 1e-9);
    }
}

#[test]
fn held_elastic_path_has_half_stiffness_and_full_surviving_dashpot_slope() {
    let viscosity = CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap();
    for law in [common().with_viscosity(viscosity), mixed().with_viscosity(viscosity)] {
        let jump = [3e-5, -2e-5, 1e-5];
        let old = law.bonded_point().trial_over_time(jump, 0.2, 1e-8)
            .unwrap().equilibrium.point;
        let (_, _, tangent) = old
            .trial_over_time_with_mean_traction_tangent(jump, 0.2, 1e-8).unwrap();
        let tangent = tangent.unwrap();
        let stiffness = law.mixed_mode_parameters().map_or([law.stiffness; 3], |p| [
            p.normal_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
        ]);
        let eta = [viscosity.normal_pa_s_m(), viscosity.tangential_pa_s_m(),
            viscosity.tangential_pa_s_m()];
        for i in 0..3 {
            for j in 0..3 {
                close(tangent[i][j], if i == j { 0.5 * stiffness[i] + eta[i] / 0.2 }
                    else { 0.0 }, 2e-15);
            }
        }
        assert_eq!(old.jump_m(), jump);
        assert_eq!(old.damage(), 0.0);
    }
}

#[test]
fn full_rupture_derivative_retains_tiny_onset_and_short_softening_limits() {
    let unit = [0.6, 0.48, -0.64];
    let radius = 0.004;
    for stiffness in [1e9, 1e24] {
        let law = BilinearCohesive::new(stiffness, 1.0, 1e-3, stiffness, 0.01).unwrap();
        let (step, mean, tangent) = law.bonded_point()
            .trial_with_mean_traction_tangent(unit.map(|value| value * radius), 1e-14)
            .unwrap();
        assert!(step.point.failed());
        let tangent = tangent.unwrap();
        // The whole radial ramp has spent Gc: mean = Gc*g/|g|^2.
        for i in 0..3 {
            close(mean[i], 1e-3 * unit[i] / radius, 1e-12);
            for j in 0..3 {
                let expected = (1e-3 / (radius * radius))
                    * (if i == j { 1.0 } else { 0.0 } - 2.0 * unit[i] * unit[j]);
                close(tangent[i][j], expected, 2e-10);
            }
        }
    }
    let law = common();
    let start = [law.failure - 1e-12, 0.0, 0.0];
    let end = [law.failure - 5e-13, 0.0, 0.0];
    let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
    let (_, _, tangent) = old.trial_with_mean_traction_tangent(end, 1e-8).unwrap();
    close(tangent.unwrap()[0][0], -0.5 * law.softening_stiffness_pa_m(), 2e-12);
}


#[test]
fn first_mixed_mode_path_derivatives_follow_selection_reclosure_and_viscosity() {
    let viscosity = CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap();
    for exponent in [1.0, 1.7, 2.5] {
        let mut parameters = mixed().mixed_mode_parameters().unwrap();
        parameters.bk_exponent = exponent;
        for viscous in [false, true] {
            let law = BilinearCohesive::new_mixed_mode(parameters, 2e9, 0.02).unwrap();
            let law = if viscous { law.with_viscosity(viscosity) } else { law };
            for (start, end) in [
                ([0.0, 0.0, 0.0], [4e-4, 6e-4, 2e-4]),
                ([6e-5, 3e-5, -2e-5], [4e-4, 6e-4, -3e-4]),
                ([-1.5e-4, 3e-5, 4e-5], [5e-4, 1.4e-3, -7e-4]),
                ([6e-5, -3e-5, 2e-5], [-5e-4, 0.012, 0.004]),
                // The first quadratic onset is exactly at the normal switch.
                ([1e-4, 0.0, 0.0], [-1e-4, 8e-4, 0.0]),
                ([0.0, 0.0, 0.0], [6e-4, 0.0, 0.0]),
                ([-1e-4, 0.0, 0.0], [-1e-4, 0.001, 0.0003]),
            ] {
                let old = law.bonded_point()
                    .trial_over_time(start, 0.2, 1e-8).unwrap().equilibrium.point;
                assert!(old.mixed_mode_history().is_none());
                check_path_tangent(old, end, 0.2, 1e-9);
                let step = old.trial_over_time(end, 0.2, 1e-8).unwrap();
                assert!(step.equilibrium.point.mixed_mode_history().is_some());
                assert!(old.mixed_mode_history().is_none());
            }
        }
    }
}

#[test]
fn first_mixed_mode_selection_can_precede_outward_energy_radius_growth() {
    let mut parameters = mixed().mixed_mode_parameters().unwrap();
    parameters.tangential_strength_pa = 2.8e5;
    let start = [1.6e-4, 3.85e-4, 0.0];
    let end = [4.4e-4, -7e-4, 0.0];
    let (history, initiation) =
        super::super::mixed_mode::initiate(parameters.into(), start, end).unwrap();
    let stiffness = [
        parameters.normal_stiffness_pa_m,
        parameters.tangential_stiffness_pa_m,
        parameters.tangential_stiffness_pa_m,
    ];
    let delta = core::array::from_fn::<_, 3, _>(|i| end[i] - start[i]);
    let at_initiation = at(start, end, initiation);
    // The strength surface has been crossed while elastic energy still falls.
    // Softening starts at a later outward crossing of the selected energy radius.
    assert!((0..3).map(|i| stiffness[i] * at_initiation[i] * delta[i]).sum::<f64>() < 0.0);
    let weights = stiffness.map(det::sqrt);
    let softening = outward_crossing(
        core::array::from_fn(|i| weights[i] * start[i]),
        core::array::from_fn(|i| weights[i] * end[i]),
        history.onset_sqrt_energy,
        initiation,
    ).unwrap();
    assert!(softening > initiation);
    for viscous in [false, true] {
        let law = BilinearCohesive::new_mixed_mode(parameters, 2e9, 0.02).unwrap();
        let law = if viscous {
            law.with_viscosity(CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap())
        } else {
            law
        };
        let old = law.bonded_point()
            .trial_over_time(start, 0.2, 1e-8).unwrap().equilibrium.point;
        assert!(old.mixed_mode_history().is_none());
        check_path_tangent(old, end, 0.2, 1e-9);
    }
}

#[test]
fn first_mixed_mode_mean_tangent_rotates_with_the_shear_basis() {
    let rotation = [[1.0, 0.0, 0.0], [0.0, 0.6, -0.8], [0.0, 0.8, 0.6]];
    let rotate = |g: [f64; 3]| rotation.map(|row| dot(row, g));
    let start = [6e-5, 3e-5, -2e-5];
    let end = [4e-4, 6e-4, -3e-4];
    for law in [
        mixed(),
        mixed().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap()),
    ] {
        let old = law.bonded_point()
            .trial_over_time(start, 0.2, 1e-8).unwrap().equilibrium.point;
        let rotated = law.bonded_point()
            .trial_over_time(rotate(start), 0.2, 1e-8).unwrap().equilibrium.point;
        let (step, mean, tangent) = old
            .trial_over_time_with_mean_traction_tangent(end, 0.2, 1e-8).unwrap();
        let (rotated_step, rotated_mean, rotated_tangent) = rotated
            .trial_over_time_with_mean_traction_tangent(rotate(end), 0.2, 1e-8).unwrap();
        close_vector(rotate(mean), rotated_mean, 3e-12);
        close(step.total_work_j_m2(), rotated_step.total_work_j_m2(), 3e-12);
        let tangent = tangent.unwrap();
        let rotated_tangent = rotated_tangent.unwrap();
        let scale = tangent.iter().flatten().fold(1.0_f64, |a, b| a.max(b.abs()));
        for i in 0..3 {
            for j in 0..3 {
                let expected = (0..3).map(|a| (0..3).map(|b| {
                    rotation[i][a] * tangent[a][b] * rotation[j][b]
                }).sum::<f64>()).sum::<f64>();
                assert!((expected - rotated_tangent[i][j]).abs() <= 3e-12 * scale);
            }
        }
        assert!(old.mixed_mode_history().is_none());
        assert!(rotated.mixed_mode_history().is_none());
    }
}

#[test]
fn first_mixed_mode_full_rupture_retains_the_bk_toughness_gradient() {
    let end = [0.003, 0.008, -0.004];
    for exponent in [1.0, 1.7, 2.5] {
        let mut p = mixed().mixed_mode_parameters().unwrap();
        p.bk_exponent = exponent;
        let law = BilinearCohesive::new_mixed_mode(p, 2e9, 0.02).unwrap();
        let (step, mean, tangent) = law.bonded_point()
            .trial_with_mean_traction_tangent(end, 1e-8).unwrap();
        assert!(step.point.failed());
        let tangent = tangent.unwrap();
        let stiffness = [p.normal_stiffness_pa_m, p.tangential_stiffness_pa_m,
            p.tangential_stiffness_pa_m];
        let radius_squared = (0..3).map(|i| stiffness[i] * end[i] * end[i]).sum::<f64>();
        let beta = p.tangential_stiffness_pa_m
            * (end[1] * end[1] + end[2] * end[2]) / radius_squared;
        let difference = p.mode_ii_fracture_energy_j_m2 - p.mode_i_fracture_energy_j_m2;
        let energy = p.mode_i_fracture_energy_j_m2 + difference * beta.powf(exponent);
        let beta_gradient = core::array::from_fn::<_, 3, _>(|j| {
            2.0 * (if j == 0 { -beta } else { 1.0 - beta })
                * stiffness[j] * end[j] / radius_squared
        });
        let energy_gradient = beta_gradient.map(|value| {
            difference * exponent * beta.powf(exponent - 1.0) * value
        });
        close(step.energy.work, energy, 2e-12);
        // Independent closed form for a fully spent proportional ramp:
        // mean_i = Gc(beta) K_i g_i / (g dot K g).
        // Differentiating Gc is essential and need not give a symmetric matrix.
        for i in 0..3 {
            close(mean[i], energy * stiffness[i] * end[i] / radius_squared, 2e-12);
            for j in 0..3 {
                let expected = if i == j { energy * stiffness[i] / radius_squared } else { 0.0 }
                    + (stiffness[i] * end[i] / radius_squared) * energy_gradient[j]
                    - 2.0 * energy * (stiffness[i] * end[i] / radius_squared)
                        * (stiffness[j] * end[j] / radius_squared);
                close(tangent[i][j], expected, 2e-10);
            }
        }
        assert!((tangent[0][1] - tangent[1][0]).abs() > 1e3);
    }
}

#[test]
fn unsupported_path_derivatives_preserve_original_force_and_explicit_refusals() {
    let common = common();
    let mixed = mixed();
    let held = common.bonded_point().trial([4e-4, 0.0, 0.0], 1e-8).unwrap().point;
    for (old, end) in [
        (mixed.bonded_point(), [0.0, 4e-4, 1e-4]), // First BK path held at normal kink.
        (common.bonded_point(), [0.0, 4e-4, 1e-4]), // Entire path at contact kink.
        (held, held.jump_m()), // Held irreversible maximum.
    ] {
        let original = old.trial_with_mean_traction(end, 1e-8).unwrap();
        let (step, mean, tangent) = old.trial_with_mean_traction_tangent(end, 1e-8).unwrap();
        assert_eq!((step, mean), original);
        assert_eq!(tangent, None);
        let timed = old.trial_over_time_with_mean_traction_tangent(end, 0.2, 1e-8).unwrap();
        assert_eq!(timed.2, None);
    }
    // A held state strictly inside the retained history is smooth.
    let unloaded = held.trial([1e-4, 0.0, 0.0], 1e-8).unwrap().point;
    assert!(unloaded.trial_with_mean_traction_tangent(unloaded.jump_m(), 1e-8)
        .unwrap().2.is_some());
    let saved = held;
    assert_eq!(
        held.trial_with_mean_traction_tangent([0.03, 0.0, 0.0], 1e-8).unwrap_err(),
        held.trial_with_mean_traction([0.03, 0.0, 0.0], 1e-8).unwrap_err(),
    );
    assert!(held.trial_over_time_with_mean_traction_tangent([4e-4, 0.0, 0.0], 0.0, 1e-8).is_err());
    assert_eq!(held, saved);
}
