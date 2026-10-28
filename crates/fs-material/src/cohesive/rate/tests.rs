use super::*;
use crate::cohesive::mixed_mode::MixedModeParameters;
fn law() -> BilinearCohesive {
    BilinearCohesive::new(1e9, 1e6, 2000.0, 2e9, 0.02)
        .unwrap()
        .with_viscosity(CohesiveViscosity::new(2e7, 3e7, 0.2).unwrap())
}
fn mixed() -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 2e9,
            normal_strength_pa: 1e6,
            tangential_strength_pa: 1.6e6,
            mode_i_fracture_energy_j_m2: 2000.0,
            mode_ii_fracture_energy_j_m2: 4500.0,
            bk_exponent: 1.5,
        },
        2e9,
        0.02,
    )
    .unwrap()
    .with_viscosity(CohesiveViscosity::new(2e7, 3e7, 0.2).unwrap())
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * (a.abs() + b.abs() + 1e-15),
        "{a} != {b}"
    );
}

#[test]
fn intact_kelvin_voigt_response_has_physical_rate_and_path_work() {
    let law = law();
    let p = law.bonded_point();
    let g = [0.0003, 0.0002, -0.0001];
    let slow = p.trial_over_time(g, 0.5, 1e-8).unwrap();
    let fast = p.trial_over_time(g, 0.25, 1e-8).unwrap();
    assert_eq!(slow.equilibrium, fast.equilibrium);
    let eta = [2e7, 3e7, 3e7];
    let mut expected_loss = 0.0;
    for i in 0..3 {
        close(slow.traction_pa[i], 1e9 * g[i] + eta[i] * g[i] / 0.5, 1e-14);
        expected_loss += eta[i] * g[i] * g[i] / 0.5;
    }
    close(slow.viscous_dissipated_j_m2, expected_loss, 1e-14);
    close(fast.viscous_dissipated_j_m2, 2.0 * expected_loss, 1e-14);
    close(
        slow.total_work_j_m2(),
        slow.equilibrium.energy.stored_change + expected_loss,
        1e-14,
    );
    assert_eq!(slow.equilibrium.energy.dissipated, 0.0);
}

#[test]
fn total_consistent_tangent_includes_damage_and_mixed_initiation_sensitivity() {
    for law in [law(), mixed()] {
        let virgin = law.bonded_point();
        let damaged = virgin
            .trial_over_time([0.002, 0.0003, 0.0], 0.2, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        for (p, g) in [
            (virgin, [0.0016, 0.0008, 0.0004]),
            (virgin, [-0.001, 0.0018, 0.0005]),
            (damaged, [0.0005, 0.0002, -0.0001]),
            (damaged, [0.0024, 0.001, 0.0005]),
        ] {
            let result = p.trial_over_time(g, 0.3, 1e-8).unwrap();
            let tangent = result.tangent_pa_m.unwrap();
            for j in 0..3 {
                let mut a = g;
                let mut b = g;
                a[j] += 1e-9;
                b[j] -= 1e-9;
                let ta = p.trial_over_time(a, 0.3, 1e-8).unwrap().traction_pa;
                let tb = p.trial_over_time(b, 0.3, 1e-8).unwrap().traction_pa;
                for i in 0..3 {
                    let fd = (ta[i] - tb[i]) / 2e-9;
                    assert!(
                        (fd - tangent[i][j]).abs() < 2e-6 * 1e9,
                        "{fd} != {}",
                        tangent[i][j]
                    );
                }
            }
        }
    }
}

#[test]
fn analytic_loss_matches_independent_nonproportional_rate_quadrature() {
    for law in [law(), mixed()] {
        let p = law
            .bonded_point()
            .trial_over_time([0.0013, 0.0004, 0.0], 0.2, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        for g in [
            [-0.001, 0.0027, 0.0006],
            [0.005, -0.001, 0.001],
            [-0.002, -0.0005, 0.0],
        ] {
            let dt = 0.4;
            let d = core::array::from_fn::<_, 3, _>(|i| g[i] - p.jump[i]);
            let result = p.trial_over_time(g, dt, 1e-8).unwrap();
            let mut quadrature = 0.0;
            let eta = law.viscosity().unwrap().coefficients();
            let factor = (0..3).map(|i| eta[i] * d[i] * d[i] / dt).sum::<f64>();
            let count = 20000;
            for n in 0..count {
                let s = (f64::from(n) + 0.5) / f64::from(count);
                let q = core::array::from_fn(|i| p.jump[i] + s * d[i]);
                let step = p.trial_equilibrium(q, 1e-8).unwrap();
                quadrature += factor * integrity(&step.point) / f64::from(count);
            }
            assert!(
                (result.viscous_dissipated_j_m2 - quadrature).abs() < 3e-5,
                "{} != {quadrature}",
                result.viscous_dissipated_j_m2
            );
            close(
                result.total_work_j_m2(),
                result.equilibrium.energy.stored_change
                    + result.equilibrium.energy.dissipated
                    + result.viscous_dissipated_j_m2,
                1e-12,
            );
        }
    }
}

#[test]
fn constant_rate_subdivision_preserves_total_fracture_and_viscous_work() {
    for law in [law(), mixed()] {
        let start = law.bonded_point();
        let end = [0.006, 0.0012, -0.0008];
        let dt = 0.8;
        let whole = start.trial_over_time(end, dt, 1e-8).unwrap();
        let mut point = start;
        let mut work = 0.0;
        let mut loss = 0.0;
        for i in 1..=80 {
            let step = point
                .trial_over_time(end.map(|v| v * f64::from(i) / 80.0), dt / 80.0, 1e-8)
                .unwrap();
            work += step.total_work_j_m2();
            loss += step.viscous_dissipated_j_m2;
            point = step.equilibrium.point;
        }
        close(work, whole.total_work_j_m2(), 2e-12);
        close(loss, whole.viscous_dissipated_j_m2, 2e-12);
        close(point.damage(), whole.equilibrium.point.damage(), 1e-12);
    }
}

#[test]
fn compression_is_viscous_only_while_bonds_survive_and_unloading_does_not_heal() {
    let law = law();
    let p = law.bonded_point();
    let compression = p.trial_over_time([-0.001, 0.0, 0.0], 0.2, 1e-8).unwrap();
    close(compression.traction_pa[0], -2.1e6, 1e-14);
    close(compression.viscous_dissipated_j_m2, 100.0, 1e-14);
    assert_eq!(compression.equilibrium.point.damage(), 0.0);
    let damaged = p
        .trial_over_time([0.002, 0.0, 0.0], 0.2, 1e-8)
        .unwrap()
        .equilibrium
        .point;
    let unload = damaged
        .trial_over_time([0.0002, 0.0, 0.0], 0.2, 1e-8)
        .unwrap();
    assert_eq!(unload.equilibrium.point.damage(), damaged.damage());
    assert!(unload.viscous_dissipated_j_m2 > 0.0);
    let failed = p
        .trial_over_time([0.006, 0.0, 0.0], 0.2, 1e-8)
        .unwrap()
        .equilibrium
        .point;
    let closed = failed
        .trial_over_time([-0.001, 0.001, -0.001], 0.2, 1e-8)
        .unwrap();
    assert_eq!(closed.traction_pa, [-2e6, 0.0, 0.0]);
    assert_eq!(closed.viscous_dissipated_j_m2, 0.0);
}

#[test]
fn tiny_softening_intervals_retain_positive_loss_and_refusals_are_atomic() {
    let law = law();
    let p = law.bonded_point();
    for fraction in [0.7, 0.999, 1.0 - 1e-8] {
        let a = law.failure_separation_m() * fraction;
        let point = p
            .trial_over_time([a, 0.0, 0.0], 0.2, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        let dx = (law.failure_separation_m() - a) * 1e-5;
        let end = [a + dx, 0.0, 0.0];
        let step = point.trial_over_time(end, 0.1, 1e-8).unwrap();
        assert!(step.viscous_dissipated_j_m2 > 0.0);
        let factor = 2e7 * (end[0] - a) * (end[0] - a) / 0.1;
        let bound = integrity(&point) * factor;
        assert!(step.viscous_dissipated_j_m2 <= bound * (1.0 + 1e-12));
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(p.trial_over_time([0.001, 0.0, 0.0], dt, 1e-8).is_err());
    }
    assert!(p.trial([0.001, 0.0, 0.0], 1e-8).is_err());
    assert!(p.trial_over_time([0.001, 0.0, 0.0], 1e-6, 1e-8).is_err());
    assert_eq!(p, law.bonded_point());
    assert!(CohesiveViscosity::new(-1.0, 1.0, 1.0).is_err());
    assert!(CohesiveViscosity::new(1.0, f64::NAN, 1.0).is_err());
}

#[test]
fn complete_rupture_ramp_retains_the_analytic_logarithmic_viscous_work() {
    // G0/G3: independent integral for a monotonic full-rupture ramp:
    // mean(1-d) = onset*failure*ln(failure/onset)/((failure-onset)*end).
    // Very small initiation fractions must not set a timestep-dependent
    // fracture threshold or erase the nonnegative logarithmic contribution.
    for stiffness in [1e9, 1e12, 1e24] {
        let law = BilinearCohesive::new(stiffness, 1.0, 1e-3, stiffness, 0.01)
            .unwrap()
            .with_viscosity(CohesiveViscosity::new(1e6, 1e6, 1.0).unwrap());
        let end = 0.004;
        let dt = 0.1;
        let onset = 1.0 / stiffness;
        let failure = 0.002;
        let mean = (onset / end) * (failure / (failure - onset)) * (failure / onset).ln();
        let expected = 1e6 * end * end / dt * mean;
        let result = law
            .bonded_point()
            .trial_over_time([end, 0.0, 0.0], dt, 1e-12)
            .unwrap();
        assert!(result.equilibrium.point.failed());
        assert_eq!(result.traction_pa, [0.0; 3]);
        assert!(
            (result.viscous_dissipated_j_m2 / expected - 1.0).abs() < 1e-12,
            "K={stiffness}: {} != {expected}",
            result.viscous_dissipated_j_m2
        );
        close(result.equilibrium.energy.dissipated, 1e-3, 1e-13);
    }
}

#[test]
fn viscous_reversal_through_compression_preserves_subdivided_work() {
    // G3: damage arrests during unloading and resumes only at the outward
    // retained-maximum crossing, possibly beyond the normal-contact switch.
    for law in [law(), mixed()] {
        let start = law
            .bonded_point()
            .trial_over_time([0.0012, 0.0004, 0.0], 0.2, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        for end in [[-0.002, -0.003, 0.0006], [0.006, -0.003, 0.0005]] {
            let whole = start.trial_over_time(end, 0.8, 1e-8).unwrap();
            let mut point = start;
            let mut loss = 0.0;
            let mut work = 0.0;
            for n in 1..=64 {
                let s = f64::from(n) / 64.0;
                let jump = core::array::from_fn(|i| (1.0 - s) * start.jump[i] + s * end[i]);
                let step = point.trial_over_time(jump, 0.8 / 64.0, 1e-8).unwrap();
                loss += step.viscous_dissipated_j_m2;
                work += step.total_work_j_m2();
                point = step.equilibrium.point;
            }
            close(loss, whole.viscous_dissipated_j_m2, 4e-12);
            close(work, whole.total_work_j_m2(), 4e-12);
            close(point.damage(), whole.equilibrium.point.damage(), 1e-13);
        }
    }
}
