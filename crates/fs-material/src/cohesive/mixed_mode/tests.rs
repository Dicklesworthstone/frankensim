use super::*;
fn parameters() -> MixedModeParameters {
    MixedModeParameters {
        normal_stiffness_pa_m: 1e9,
        tangential_stiffness_pa_m: 4e8,
        normal_strength_pa: 2e5,
        tangential_strength_pa: 1.6e5,
        mode_i_fracture_energy_j_m2: 100.0,
        mode_ii_fracture_energy_j_m2: 400.0,
        bk_exponent: 2.0,
    }
}
fn law() -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(parameters(), 2e9, 0.02).unwrap()
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

#[test]
fn independent_pure_modes_and_proportional_bk_mix_spend_their_own_toughness() {
    let law = law();
    let p = parameters();
    for direction in [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.6, -0.8],
        [0.6, 0.8, 0.0],
    ] {
        let stiffness = p.normal_stiffness_pa_m * direction[0] * direction[0]
            + p.tangential_stiffness_pa_m
                * (direction[1] * direction[1] + direction[2] * direction[2]);
        let beta = p.tangential_stiffness_pa_m
            * (direction[1] * direction[1] + direction[2] * direction[2])
            / stiffness;
        let fracture = 100.0 + 300.0 * beta * beta;
        let onset = 1.0
            / ((p.normal_stiffness_pa_m * direction[0] / p.normal_strength_pa).powi(2)
                + (p.tangential_stiffness_pa_m * direction[1] / p.tangential_strength_pa).powi(2)
                + (p.tangential_stiffness_pa_m * direction[2] / p.tangential_strength_pa).powi(2))
            .sqrt();
        let failure = 2.0 * fracture / (stiffness * onset);
        let mut point = law.bonded_point();
        let mut work = 0.0;
        for step in 1..=120 {
            let distance = failure * f64::from(step) / 100.0;
            let trial = point.trial(direction.map(|x| x * distance), 1e-8).unwrap();
            assert!(trial.point.damage() >= point.damage());
            assert!(trial.energy.dissipated >= 0.0);
            work += trial.energy.work;
            point = trial.point;
        }
        assert!(point.failed());
        close(work, fracture, 2e-12);
        close(point.dissipated_j_m2(), fracture, 2e-12);
        close(
            point.mixed_mode_history().unwrap().energy_shear_fraction,
            beta,
            1e-13,
        );
        assert_eq!(point.traction_pa(), [0.0; 3]);
    }
}

#[test]
fn first_initiation_sensitive_and_retained_history_tangents_match_differences() {
    let law = law();
    let virgin = law.bonded_point();
    let elastic = virgin.trial([6e-5, 3e-5, -2e-5], 1e-8).unwrap().point;
    let damaged = virgin.trial([0.0004, 0.0003, 0.0002], 1e-8).unwrap().point;
    for (point, end) in [
        (virgin, [0.0004, 0.0006, 0.0002]),
        (elastic, [0.0004, 0.0006, -0.0003]),
        (damaged, [0.0006, 0.001, 0.0004]),
        (damaged, [-0.0002, 0.001, -0.0003]),
        (damaged, [0.0001, 0.0001, -0.0001]),
    ] {
        let tangent = point.trial(end, 1e-8).unwrap().tangent_pa_m.unwrap();
        for j in 0..3 {
            let h = 2e-9;
            let mut plus = end;
            let mut minus = end;
            plus[j] += h;
            minus[j] -= h;
            let a = point.trial(plus, 1e-8).unwrap().traction_pa;
            let b = point.trial(minus, 1e-8).unwrap().traction_pa;
            for i in 0..3 {
                assert!(
                    ((a[i] - b[i]) / (2.0 * h) - tangent[i][j]).abs() < 2e-7 * 1e9,
                    "i={i},j={j},fd={},exact={}",
                    (a[i] - b[i]) / (2.0 * h),
                    tangent[i][j]
                );
            }
        }
    }
}

#[test]
fn nonproportional_work_is_exact_and_onset_mix_does_not_change_or_heal() {
    let law = law();
    let point = law
        .bonded_point()
        .trial([0.0004, 0.0002, 0.0001], 1e-8)
        .unwrap()
        .point;
    let end = [-0.0003, 0.0015, -0.0005];
    let delta = core::array::from_fn::<_, 3, _>(|i| end[i] - point.jump[i]);
    let mut independent = 0.0;
    for i in 0..20_000 {
        let s = (f64::from(i) + 0.5) / 20_000.0;
        let g = core::array::from_fn(|j| point.jump[j] + s * delta[j]);
        independent += dot(point.trial(g, 1e-8).unwrap().traction_pa, delta) / 20_000.0;
    }
    let trial = point.trial(end, 1e-8).unwrap();
    assert!((trial.energy.work - independent).abs() < 2e-5);
    assert_eq!(trial.point.mixed_mode_history(), point.mixed_mode_history());
    let unload = trial.point.trial([-0.0001, 0.0001, 0.0], 1e-8).unwrap();
    close(unload.point.damage(), trial.point.damage(), 1e-14);
    assert_eq!(unload.energy.dissipated, 0.0);
    let failed = unload.point.trial([0.0, 0.015, 0.0], 1e-8).unwrap().point;
    assert!(failed.failed());
    let closed = failed.trial([-0.001, 0.001, 0.0], 1e-8).unwrap();
    assert_eq!(closed.traction_pa, [-2e6, 0.0, 0.0]);
    assert_eq!(closed.energy.dissipated, 0.0);
}

#[test]
fn equal_mode_parameters_recover_original_law_on_nonproportional_paths() {
    let standard = BilinearCohesive::new(1e9, 1e6, 2000.0, 2e9, 0.02).unwrap();
    let mixed = BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 1e9,
            normal_strength_pa: 1e6,
            tangential_strength_pa: 1e6,
            mode_i_fracture_energy_j_m2: 2000.0,
            mode_ii_fracture_energy_j_m2: 2000.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap();
    let mut a = standard.bonded_point();
    let mut b = mixed.bonded_point();
    for end in [
        [0.0002, 0.0001, 0.0],
        [0.0015, 0.0009, 0.0003],
        [-0.001, 0.002, 0.0005],
        [0.0001, 0.0002, 0.0],
        [0.003, 0.001, 0.0005],
        [0.0, 0.006, 0.0],
        [-0.001, 0.001, 0.0],
    ] {
        let x = a.trial(end, 1e-8).unwrap();
        let y = b.trial(end, 1e-8).unwrap();
        for i in 0..3 {
            close(x.traction_pa[i], y.traction_pa[i], 3e-12);
        }
        for (u, v) in [
            (x.energy.work, y.energy.work),
            (x.energy.dissipated, y.energy.dissipated),
            (x.energy.stored_change, y.energy.stored_change),
            (x.point.damage(), y.point.damage()),
        ] {
            close(u, v, 2e-11);
        }
        for i in 0..3 {
            for j in 0..3 {
                close(
                    x.branch_tangent_pa_m[i][j],
                    y.branch_tangent_pa_m[i][j],
                    2e-10,
                );
            }
        }
        a = x.point;
        b = y.point;
    }
}

#[test]
fn straight_path_subdivision_keeps_first_mode_selection_work_and_damage() {
    let law = law();
    let start = law
        .bonded_point()
        .trial([0.00007, 0.00003, 0.0], 1e-8)
        .unwrap()
        .point;
    let end = [0.0007, 0.0012, -0.0004];
    let one = start.trial(end, 1e-8).unwrap();
    let mut point = start;
    let mut work = 0.0;
    let mut loss = 0.0;
    for i in 1..=80 {
        let s = f64::from(i) / 80.0;
        let g = core::array::from_fn(|j| start.jump[j] + s * (end[j] - start.jump[j]));
        let step = point.trial(g, 1e-8).unwrap();
        work += step.energy.work;
        loss += step.energy.dissipated;
        point = step.point;
    }
    close(work, one.energy.work, 1e-11);
    close(loss, one.energy.dissipated, 1e-11);
    close(point.damage(), one.point.damage(), 1e-12);
    close(
        point.mixed_mode_history().unwrap().energy_shear_fraction,
        one.point
            .mixed_mode_history()
            .unwrap()
            .energy_shear_fraction,
        1e-12,
    );
}

#[test]
fn compression_shear_and_invalid_data_keep_explicit_bounds() {
    let law = law();
    let compression = law
        .bonded_point()
        .trial([-0.01, 0.0, 0.0], 1e-8)
        .unwrap()
        .point;
    assert_eq!(compression.damage(), 0.0);
    assert!(compression.mixed_mode_history().is_none());
    let shear = compression.trial([-0.01, 0.006, 0.0], 1e-7).unwrap();
    assert!(shear.point.failed());
    close(shear.energy.dissipated, 400.0, 1e-12);
    let saved = shear.point;
    assert!(saved.trial([0.1, 0.0, 0.0], 1e-8).is_err());
    assert_eq!(saved, shear.point);
    let mut p = parameters();
    p.bk_exponent = 0.5;
    assert!(BilinearCohesive::new_mixed_mode(p, 2e9, 0.02).is_err());
    let mut p = parameters();
    p.tangential_strength_pa = -1.0;
    assert!(BilinearCohesive::new_mixed_mode(p, 2e9, 0.02).is_err());
    assert!(BilinearCohesive::new_mixed_mode(parameters(), 2e9, 0.002).is_err());
}
