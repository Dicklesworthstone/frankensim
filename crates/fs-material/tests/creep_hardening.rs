//! G0/G1: independent primary-creep solutions and exact secondary reductions.
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::uniaxial::UniaxialNortonJ2;
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{
    CreepState, CreepTolerances, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::plastic::J2Plasticity;
use fs_material::tensor::{rotate, rotation, von_mises};
use fs_material::{IsotropicElastic, SmallStrainLaw, Uniaxial, Voigt};

const E: f64 = 1.2e9;
const STRESS: f64 = 3e6;
const RATE: f64 = 2e-3;
const A: f64 = 8e-4;
fn elastic() -> IsotropicElastic {
    IsotropicElastic::new(E, 0.2, 0.05).unwrap()
}
fn norton(n: f64, dt: f64) -> NortonCreep {
    NortonCreep::new(
        elastic(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        500.0,
        dt,
    )
    .unwrap()
}
fn primary(n: f64, dt: f64, m: f64, f: f64) -> NortonCreep {
    norton(n, dt)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: m,
            residual_rate_fraction: f,
        })
        .unwrap()
}
fn tolerance() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 3e-6,
        energy_j_m3: 2e-7,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b} (tol {tol})");
}
fn dev(a: f64) -> Voigt {
    [a, -0.5 * a, -0.5 * a, 0.0, 0.0, 0.0]
}
// Positive root of d^2 + b*d - c = 0, without cancellation.
fn quadratic(b: f64, c: f64) -> f64 {
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[test]
fn primary_tensor_and_axial_returns_match_independent_quadratic_and_work() {
    let dt = 0.7;
    let law = primary(1.0, dt, 1.0, 0.0);
    let old = CreepState {
        creep_strain: dev(2e-4),
        equivalent_creep_strain: 9e-4,
    };
    let strain: Voigt = core::array::from_fn(|i| old.creep_strain[i] + dev(0.004)[i]);
    let h = 3.0 * law.elastic().lame().1;
    let q_trial = h * 0.004;
    let k = dt * RATE * A / STRESS;
    let d = quadratic(A + old.equivalent_creep_strain + k * h, k * q_trial);
    let q = q_trial - h * d;
    let result = law
        .advance_with_energy(&old.creep_strain, &old, &strain, tolerance())
        .unwrap();
    close(
        result.state.equivalent_creep_strain - old.equivalent_creep_strain,
        d,
        3e-17,
    );
    close(von_mises(&result.stress_pa), q, 3e-8);
    close(result.energy.creep_dissipation_j_m3, q * d, 1e-8);
    close(result.energy.physical_dissipation_j_m3, q * d, 1e-8);
    close(result.energy.residual_j_m3, 0.0, 2e-7);
    let radial = 1.0 / (1.0 + h * k / (A + old.equivalent_creep_strain + 2.0 * d));
    let action: f64 = (0..6).map(|j| result.tangent_pa[0][j] * dev(1.0)[j]).sum();
    close(action, 2.0 * law.elastic().lame().1 * radial, 2e-5);
    let axial = UniaxialNorton::new(law);
    let state = axial.initial_state();
    let d = quadratic(A + k * E, k * E * 0.004);
    let response = axial.response(0.004, &state).unwrap();
    close(response.state.creep_strain(), d, 3e-17);
    close(response.stress_pa, E * (0.004 - d), 3e-8);
    close(response.tangent_pa, E / (1.0 + E * k / (A + 2.0 * d)), 3e-5);
    axial
        .advance_with_energy(0.0, &state, 0.004, tolerance())
        .unwrap();
    let zero = axial.response(0.0, &state).unwrap();
    close(zero.tangent_pa, E / (1.0 + E * dt * RATE / STRESS), 3e-5);
    assert_eq!(zero.state, state);
}

#[test]
fn sub_stress_roundoff_creep_still_solves_its_history_and_tiny_hardening_remains_admissible() {
    let dt = 1.5e-17;
    let a = 1e-21;
    let law = norton(1.0, dt)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: a,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap();
    let h = 3.0 * elastic().lame().1;
    let strain = dev(1e6 / h);
    let result = law
        .advance_with_energy(&strain, &law.initial_state(), &strain, tolerance())
        .unwrap();
    let k = dt * RATE * a / STRESS;
    let expected = quadratic(a + k * h, k * 1e6);
    close(result.state.equivalent_creep_strain, expected, 2e-33);
    let rate = law
        .equivalent_creep_rate_at_strain_s(
            von_mises(&result.stress_pa),
            result.state.equivalent_creep_strain,
        )
        .unwrap();
    close(result.state.equivalent_creep_strain, dt * rate, 2e-33);
    // The positive effective stiffness is representable even though H/(3G+H)
    // is not. Its logarithm must reach the primary return without rounding to 0.
    let combined = |hardening| {
        NortonJ2::new(
            primary(1.0, 0.25, 1.0, 0.0),
            J2Plasticity::new(elastic(), 2e6, hardening).unwrap(),
        )
        .unwrap()
    };
    // At 1e-316, H/(3G+H) still rounds to zero while the stored hardening energy
    // H alpha^2 / 2 stays representable. At 1e-320 that energy underflowed as well, and the
    // energy account refuses a nonzero energy that underflows (plastic/energy.rs, product).
    let tiny_hardening = 1e-316;
    assert_eq!(
        tiny_hardening / (3.0 * elastic().lame().1 + tiny_hardening),
        0.0
    );
    let tiny = combined(tiny_hardening);
    let perfect = combined(0.0);
    let p = perfect
        .advance_with_energy(
            &[0.0; 6],
            &perfect.initial_state(),
            &dev(0.006),
            tolerance(),
        )
        .unwrap();
    let t = tiny
        .advance_with_energy(&[0.0; 6], &tiny.initial_state(), &dev(0.006), tolerance())
        .unwrap();
    close(
        t.state.creep.equivalent_creep_strain,
        p.state.creep.equivalent_creep_strain,
        3e-16,
    );
    close(t.state.plastic.alpha, p.state.plastic.alpha, 3e-16);
}

#[test]
fn constant_stress_primary_creep_converges_to_integrated_strain_hardening_law() {
    let sigma = 2e6;
    let r = RATE * (sigma / STRESS).powi(3);
    let total = 4.0;
    let exact = A * ((1.0 + 2.0 * r * total / A).sqrt() - 1.0);
    let elastic_strain = [sigma / E, -0.2 * sigma / E, -0.2 * sigma / E, 0.0, 0.0, 0.0];
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let dt = total / f64::from(steps);
        let law = primary(3.0, dt, 1.0, 0.0);
        let mut state = law.initial_state();
        let mut old_strain = elastic_strain;
        let mut alpha = 0.0;
        let mut previous_rate = r;
        for _ in 0..steps {
            let d = quadratic(A + alpha, dt * r * A);
            alpha += d;
            let strain = core::array::from_fn(|i| elastic_strain[i] + dev(alpha)[i]);
            let result = law
                .advance_with_energy(&old_strain, &state, &strain, tolerance())
                .unwrap();
            close(result.stress_pa[0], sigma, 5e-8);
            close(result.state.equivalent_creep_strain, alpha, 6e-17);
            assert!(d / dt < previous_rate);
            previous_rate = d / dt;
            state = result.state;
            old_strain = strain;
        }
        errors.push((alpha - exact).abs());
        let unloaded = law
            .advance_with_energy(&old_strain, &state, &state.creep_strain, tolerance())
            .unwrap();
        assert_eq!(unloaded.state, state);
        assert_eq!(unloaded.energy.creep_dissipation_j_m3, 0.0);
        assert!(unloaded.stress_pa.iter().all(|s| *s == 0.0));
    }
    for e in errors.windows(2) {
        assert!(e[0] / e[1] > 1.8 && e[0] / e[1] < 2.2, "{errors:?}");
    }
}

#[test]
fn nonproportional_primary_flow_has_exact_tangent_and_rotates_with_history() {
    let law = primary(3.2, 0.4, 0.8, 0.06);
    let first_strain = [0.002, -0.0005, 0.0001, 0.001, -0.0002, 0.0006];
    let old = law.update_state(&first_strain, &law.initial_state());
    let frozen = old.clone();
    let strain = [-0.001, 0.0018, -0.0004, 0.0002, 0.0012, -0.0008];
    let result = law
        .advance_with_energy(&first_strain, &old, &strain, tolerance())
        .unwrap();
    assert!(result.energy.creep_dissipation_j_m3 > 0.0);
    for column in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        let step = 1e-8;
        plus[column] += step;
        minus[column] -= step;
        let p = law.stress(&plus, &old);
        let m = law.stress(&minus, &old);
        for row in 0..6 {
            close(
                (p[row] - m[row]) / (2.0 * step),
                result.tangent_pa[row][column],
                0.4,
            );
        }
    }
    // rotation() takes a unit axis. [0.4, -0.3, 0.8] has length 0.943, and the matrix it gave
    // was not orthogonal: it moved the creep strain's trace to -1.2e-7, which the law refuses.
    let axis: [f64; 3] = [0.4, -0.3, 0.8];
    let length = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let rotation = rotation(axis.map(|a| a / length), 0.71);
    let rotated = CreepState {
        creep_strain: rotate(&old.creep_strain, &rotation),
        equivalent_creep_strain: old.equivalent_creep_strain,
    };
    let actual = law.stress(&rotate(&strain, &rotation), &rotated);
    let expected = rotate(&result.stress_pa, &rotation);
    for i in 0..6 {
        close(actual[i], expected[i], 5e-8);
    }
    assert_eq!(old, frozen);
    assert_eq!(law.stress(&strain, &old), result.stress_pa);
}

#[test]
fn primary_creep_and_yielding_solve_both_rates_simultaneously_including_perfect_plasticity() {
    let h = 3.0 * elastic().lame().1;
    let q_trial = h * 0.006;
    let yield_stress = 2e6;
    let k = 0.25 * RATE * A / STRESS;
    for hardening in [0.0, 1e8] {
        let law = NortonJ2::new(
            primary(1.0, 0.25, 1.0, 0.0),
            J2Plasticity::new(elastic(), yield_stress, hardening).unwrap(),
        )
        .unwrap();
        let effective = yield_stress + hardening / (h + hardening) * (q_trial - yield_stress);
        let modulus = h * hardening / (h + hardening);
        let dc = quadratic(A + k * modulus, k * effective);
        let q = effective - modulus * dc;
        let dp = (q_trial - yield_stress - h * dc) / (h + hardening);
        let old = law.initial_state();
        let result = law
            .advance_with_energy(&[0.0; 6], &old, &dev(0.006), tolerance())
            .unwrap();
        close(result.state.creep.equivalent_creep_strain, dc, 3e-17);
        close(result.state.plastic.alpha, dp, 3e-17);
        close(von_mises(&result.stress_pa), q, 4e-8);
        close(result.energy.creep_dissipation_j_m3, q * dc, 1e-8);
        close(
            result.energy.plastic_dissipation_j_m3,
            yield_stress * dp,
            1e-8,
        );
        let axial = UniaxialNortonJ2::new(law);
        let old = axial.initial_state();
        let result = axial.response(0.006, &old).unwrap();
        let step = 1e-8;
        let fd =
            (axial.stress(0.006 + step, &old) - axial.stress(0.006 - step, &old)) / (2.0 * step);
        close(result.tangent_pa, fd, 0.2);
        axial
            .advance_with_energy(0.0, &old, 0.006, tolerance())
            .unwrap();
    }
    // Secondary creep alone relaxes this trial below yield, whereas primary
    // hardening slows it sufficiently to activate plasticity in the same step.
    let c = primary(1.0, 0.25, 1.0, 0.0);
    let coupled =
        NortonJ2::new(c, J2Plasticity::new(elastic(), yield_stress, 1e8).unwrap()).unwrap();
    let result = coupled
        .advance_with_energy(
            &[0.0; 6],
            &coupled.initial_state(),
            &dev(0.0016),
            tolerance(),
        )
        .unwrap();
    assert!(result.state.plastic.alpha > 0.0);
}

#[test]
fn residual_secondary_rate_and_exact_disabled_reductions_are_explicit() {
    let law = primary(2.5, 0.5, 0.7, 0.15);
    let mut previous = RATE;
    for alpha in [0.0, A, 4.0 * A, 100.0 * A, 1e12 * A] {
        let expected = RATE * (0.15 + 0.85 * (1.0 + alpha / A).powf(-0.7));
        let actual = law
            .equivalent_creep_rate_at_strain_s(STRESS, alpha)
            .unwrap();
        close(actual, expected, 2e-17);
        assert!(actual <= previous + 2e-17 && actual >= 0.15 * RATE);
        previous = actual;
    }
    assert!(law.equivalent_creep_rate_at_strain_s(STRESS, -1.0).is_err());
    for n in [1.0, 3.0] {
        let plain = norton(n, 0.7);
        for (m, f) in [(0.0, 0.0), (2.0, 1.0)] {
            let reduced = primary(n, 0.7, m, f);
            for strain in [[0.0; 6], [0.003, -0.001, 0.0004, 0.001, 0.0, 0.0]] {
                let old = plain.initial_state();
                let expected = plain
                    .advance_with_energy(&[0.0; 6], &old, &strain, tolerance())
                    .unwrap();
                let actual = reduced
                    .advance_with_energy(&[0.0; 6], &old, &strain, tolerance())
                    .unwrap();
                assert_eq!(actual, expected);
            }
        }
    }
    for p in [
        StrainHardeningCreep {
            reference_strain: 0.0,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        },
        StrainHardeningCreep {
            reference_strain: A,
            exponent: -1.0,
            residual_rate_fraction: 0.0,
        },
        StrainHardeningCreep {
            reference_strain: A,
            exponent: 1.0,
            residual_rate_fraction: 1.01,
        },
        StrainHardeningCreep {
            reference_strain: A,
            exponent: f64::NAN,
            residual_rate_fraction: 0.0,
        },
    ] {
        assert!(norton(1.0, 0.25).with_strain_hardening(p).is_err());
    }
}
