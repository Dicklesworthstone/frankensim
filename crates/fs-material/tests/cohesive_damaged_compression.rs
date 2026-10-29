//! The damaged-area measure is checked against independent original-law samples.
use fs_material::cohesive::mixed_mode::MixedModeParameters;
use fs_material::cohesive::{BilinearCohesive, CohesivePoint, CohesiveViscosity};

fn laws() -> [BilinearCohesive; 2] {
    [
        BilinearCohesive::new(1e9, 2e5, 100.0, 2e9, 0.02).unwrap(),
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
        .unwrap(),
    ]
}
fn at(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    core::array::from_fn(|i| (1.0 - s) * a[i] + s * b[i])
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn quadrature(old: CohesivePoint<'_>, end: [f64; 3]) -> f64 {
    let n = 65_536;
    let mut value = 0.0;
    for i in 0..n {
        let g = at(old.jump_m(), end, (f64::from(i) + 0.5) / f64::from(n));
        let state = old.trial(g, 1e-8).unwrap().point;
        value += state.damage() * old.law().compression_stiffness_pa_m() * (-g[0]).max(0.0)
            / f64::from(n);
    }
    value
}

#[test]
fn original_damage_pressure_integral_handles_onset_reversal_sign_switch_and_failure() {
    for law in laws() {
        for (start, end) in [
            ([-2e-4, 0.0, 0.0], [-4e-4, 8e-4, 3e-4]),
            ([4e-4, 2e-4, 1e-4], [-4e-4, -0.003, 0.001]),
            ([-3e-4, 2e-4, 0.0], [4e-4, 0.012, -0.004]),
            ([-1e-4, 0.0006, -0.0003], [-3e-4, 0.0001, 0.0002]),
            ([-1e-4, 0.014, 0.0], [-4e-4, -0.003, 0.001]),
            ([-3e-4, 0.0, 0.0], [2e-4, 0.0, 0.0]),
            ([1e-4, 0.0, 0.0], [4e-4, 0.004, 0.0]),
        ] {
            let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
            let (step, mean, p) = old
                .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
                .unwrap();
            let original = old
                .trial_over_time_with_mean_traction(end, 0.01, 1e-8)
                .unwrap();
            assert_eq!((step, mean), original);
            assert!(p >= 0.0);
            close(p, quadrature(old, end), 2e-7);
        }
    }
}

#[test]
fn measure_is_additive_uses_current_damage_along_path_and_introduces_no_storage() {
    for law in laws() {
        let start = [-3e-4, 7e-5, 2e-5];
        let end = [-5e-4, 0.007, -0.003];
        let initial = law.bonded_point().trial(start, 1e-8).unwrap().point;
        let (whole, _, p) = initial
            .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
            .unwrap();
        let mut point = initial;
        let mut divided = 0.0;
        let mut work = 0.0;
        for i in 1..=64 {
            let (step, _, pressure) = point
                .trial_over_time_with_damage_weighted_compression(
                    at(start, end, f64::from(i) / 64.0),
                    0.01 / 64.0,
                    1e-8,
                )
                .unwrap();
            divided += pressure / 64.0;
            work += step.total_work_j_m2();
            point = step.equilibrium.point;
        }
        close(p, divided, 1e-11);
        close(work, whole.total_work_j_m2(), 1e-11);
        close(point.damage(), whole.equilibrium.point.damage(), 1e-12);
        close(
            point.stored_energy_j_m2(),
            whole.equilibrium.point.stored_energy_j_m2(),
            1e-12,
        );
        let endpoint_shortcut = whole.equilibrium.point.damage() * 2e9 * 0.0004;
        assert!(p < endpoint_shortcut * 0.98);
    }
}

#[test]
fn damaged_pressure_is_independent_of_bond_viscosity_and_resolves_tiny_damage_onset() {
    let law = laws()[0].clone();
    let viscous = law
        .clone()
        .with_viscosity(CohesiveViscosity::new(1e5, 2e5, 100.0).unwrap());
    let start = [-2e-4, 0.0, 0.0];
    let end = [-4e-4, 6e-4, 0.0];
    let original = law.bonded_point().trial(start, 1e-8).unwrap().point;
    let damped = viscous
        .bonded_point()
        .trial_over_time(start, 0.01, 1e-8)
        .unwrap()
        .equilibrium
        .point;
    let a = original
        .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
        .unwrap();
    let b = damped
        .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
        .unwrap();
    assert_eq!(a.2, b.2);
    assert!(b.0.viscous_dissipated_j_m2 > 0.0);
    assert_ne!(a.1, b.1);
    let onset = law.onset_separation_m();
    let old = law
        .bonded_point()
        .trial([-2e-4, onset, 0.0], 1e-8)
        .unwrap()
        .point;
    for fraction in [1e-10, 1e-8, 1e-6, 1e-4] {
        let end = [-2e-4, onset * (1.0 + fraction), 0.0];
        let (_, _, p) = old
            .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
            .unwrap();
        assert!(p > 0.0);
        let independent = quadrature(old, end);
        assert!(
            (p - independent).abs() < 2e-10 + 1e-6 * independent.abs(),
            "{fraction}: {p:e} != {independent:e}"
        );
    }
}
