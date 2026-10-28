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
