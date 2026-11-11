//! Independent finite primary/recovery references on synthetic authored laws.
//! These cases verify the discrete mechanics, not measured material calibration.

use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::FiniteInelasticMaterial;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::creep::recovery::{
    FiniteRecoveringCreepState, FiniteRecoveringNortonCreep, FiniteRecoveringNortonJ2,
};
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{
    FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2Response, FiniteJ2State,
};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const A: f64 = 0.025;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn native(eta: f64, n: f64, dt: f64) -> FiniteNortonCreep {
    FiniteNortonCreep::new(
        MU,
        K,
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 3.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
        FiniteOverstressParameters {
            viscosity_pa_s: eta,
            stress_scale_pa: 2e6,
            exponent: n,
        },
        dt,
    )
    .unwrap()
}

fn primary(eta: f64, n: f64, dt: f64) -> FiniteNortonCreep {
    native(eta, n, dt)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap()
}

fn law(eta: f64, n: f64, dt: f64, k: f64) -> FiniteRecoveringNortonCreep {
    FiniteRecoveringNortonCreep::new(primary(eta, n, dt), k).unwrap()
}

fn coaxial(mean: f64, axial: f64) -> [f64; 9] {
    [
        (mean + axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
    ]
}

fn history(net: f64, alpha: f64, resistance: f64) -> FiniteRecoveringCreepState {
    FiniteRecoveringCreepState {
        native: FiniteJ2State {
            inverse_plastic_cauchy_green: [(-2.0 * net).exp(), net.exp(), net.exp(), 0.0, 0.0, 0.0],
            alpha,
            creep_alpha: 0.0,
        },
        resistance_strain: resistance,
    }
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}, difference={:.8e}",
        (actual - expected).abs()
    );
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn same_native(actual: &FiniteJ2Response<FiniteRecoveringCreepState>, expected: &FiniteJ2Response) {
    assert_eq!(actual.state.native, expected.state);
    assert_eq!(actual.first_piola_pa, expected.first_piola_pa);
    assert_eq!(actual.kirchhoff_pa, expected.kirchhoff_pa);
    assert_eq!(actual.tangent_pa, expected.tangent_pa);
    assert_eq!(
        actual.incremental_potential_j_m3,
        expected.incremental_potential_j_m3
    );
    assert_eq!(actual.flow_potential_j_m3, expected.flow_potential_j_m3);
    assert_eq!(actual.stored_energy_j_m3, expected.stored_energy_j_m3);
    assert_eq!(
        actual.plastic_dissipation_j_m3,
        expected.plastic_dissipation_j_m3
    );
    assert_eq!(actual.plastic_increment, expected.plastic_increment);
    assert_eq!(
        actual.creep_dissipation_j_m3,
        expected.creep_dissipation_j_m3
    );
    assert_eq!(actual.creep_increment, expected.creep_increment);
    assert_eq!(
        actual.equivalent_kirchhoff_pa,
        expected.equivalent_kirchhoff_pa
    );
    assert_eq!(actual.jacobian, expected.jacobian);
}

#[test]
fn zero_recovery_reduces_exactly_to_original_finite_primary_return() {
    let original = primary(20e6, 2.5, 0.3)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: 0.015,
            exponent: 0.8,
            residual_rate_fraction: 0.17,
        })
        .unwrap();
    let material = FiniteRecoveringNortonCreep::new(original.clone(), 0.0).unwrap();
    assert_eq!(material.norton(), &original);
    assert_eq!(material.recovery_rate_s(), 0.0);
    let mut state = material.initial_state();
    for f in [
        coaxial(0.01, 0.07),
        [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96],
        coaxial(-0.02, -0.12),
    ] {
        let old = state.clone();
        let expected = original.response(&f, &old.native).unwrap();
        let actual = material.response(&f, &old).unwrap();
        same_native(&actual, &expected);
        assert_eq!(actual.state.resistance_strain, actual.state.native.alpha);
        assert_eq!(old, state);
        assert_eq!(
            material.stress_and_potential(&f, &old).unwrap(),
            (actual.first_piola_pa, actual.incremental_potential_j_m3)
        );
        state = actual.state;
    }
}

#[test]
fn recovering_finite_return_matches_independent_quadratic_and_cubic_potential() {
    let (eta, dt, k, net, mean, axial) = (60e6, 0.4, 1.3, 0.04, 0.02, 0.18);
    let material = law(eta, 1.0, dt, k);
    let old = history(net, 0.08, 0.015);
    let f = coaxial(mean, axial);
    let r = material.response(&f, &old).unwrap();
    let modulus = eta / dt;
    let trial = 3.0 * MU * (axial - net);
    let divisor = 1.0 + k * dt;
    let a = A * divisor;
    // c^2 + (a+r_old+3mu*a*dt/eta)*c - a*dt*q_trial/eta = 0.
    let b = a + old.resistance_strain + a * 3.0 * MU / modulus;
    let drive = a * trial / modulus;
    let dc = 2.0 * drive / (b + (b * b + 4.0 * drive).sqrt());
    let q = trial - 3.0 * MU * dc;
    let stored = 0.5 * K * (3.0 * mean).powi(2) + q * q / (6.0 * MU);
    let potential =
        modulus * (0.5 * (1.0 + old.resistance_strain / a) * dc * dc + dc.powi(3) / (3.0 * a));
    near(r.plastic_increment, dc, 2e-13, 2e-11);
    near(r.state.native.alpha, old.native.alpha + dc, 2e-13, 2e-11);
    near(
        r.state.resistance_strain,
        (old.resistance_strain + dc) / divisor,
        2e-13,
        2e-11,
    );
    near(r.equivalent_kirchhoff_pa, q, 3e-5, 2e-11);
    near(
        material.endpoint_equivalent_stress_from(dc, &old).unwrap(),
        q,
        3e-5,
        2e-11,
    );
    near(r.stored_energy_j_m3, stored, 2e-5, 2e-11);
    near(r.inelastic_dissipation_j_m3(), q * dc, 2e-5, 2e-11);
    near(r.flow_potential_j_m3, potential, 2e-5, 2e-10);
    near(
        r.incremental_potential_j_m3,
        stored + potential,
        2e-5,
        2e-10,
    );
    assert_eq!(r.state.native.creep_alpha, 0.0);
    assert_eq!(r.creep_increment, 0.0);
    assert_eq!(r.creep_dissipation_j_m3, 0.0);
    let expected = history(net + dc, old.native.alpha + dc, 0.0);
    for (actual, expected) in r
        .state
        .native
        .inverse_plastic_cauchy_green
        .iter()
        .zip(expected.native.inverse_plastic_cauchy_green)
    {
        near(*actual, expected, 2e-13, 2e-11);
    }
    material.validate_committed(&f, &r.state).unwrap();
    near(
        material.stored_energy(&f, &r.state).unwrap(),
        stored,
        2e-5,
        2e-11,
    );
    assert_eq!(
        material
            .norton()
            .strain_hardening()
            .unwrap()
            .reference_strain,
        A
    );
    assert_eq!(old, history(net, 0.08, 0.015));
}

#[test]
fn sustained_kirchhoff_stress_preserves_stationary_resistance_and_permanent_flow() {
    let (eta, k, q) = (20e6, 1.2, 1e6);
    let unresisted_rate = q / eta;
    // Solve k*r*(a+r)=a*rate_Norton(q); the stationary physical creep rate is k*r.
    let equilibrium =
        2.0 * A * unresisted_rate / k / (A + (A * A + 4.0 * A * unresisted_rate / k).sqrt());
    let mut old = history(0.04, 0.08, equilibrium);
    let mut net = 0.04;
    let initial_alpha = old.native.alpha;
    let mut elapsed = 0.0;
    for dt in [0.1, 0.3, 0.2, 0.5] {
        let material = law(eta, 1.0, dt, k);
        let dc = dt * k * equilibrium;
        net += dc;
        elapsed += dt;
        let f = coaxial(0.02, net + q / (3.0 * MU));
        let step = material.response(&f, &old).unwrap();
        near(step.plastic_increment, dc, 2e-13, 2e-10);
        near(step.equivalent_kirchhoff_pa, q, 3e-5, 2e-10);
        near(step.state.resistance_strain, equilibrium, 2e-13, 2e-10);
        near(
            step.state.native.alpha,
            initial_alpha + elapsed * k * equilibrium,
            2e-13,
            2e-10,
        );
        near(step.inelastic_dissipation_j_m3(), q * dc, 2e-5, 2e-10);
        old = step.state;
    }
}

#[test]
fn stress_free_finite_set_recovers_only_resistance_then_creeps_faster_on_reload() {
    let material = law(30e6, 2.0, 0.5, 1.5);
    // This finite permanent set has exactly representable F*C_c^-1*F^T=I.
    // Its cumulative path can exceed its net strain after prior load reversal.
    let set = [2.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0];
    let original = FiniteRecoveringCreepState {
        native: FiniteJ2State {
            inverse_plastic_cauchy_green: [0.25, 4.0, 1.0, 0.0, 0.0, 0.0],
            alpha: 1.2,
            creep_alpha: 0.0,
        },
        resistance_strain: 0.3,
    };
    let mut state = original.clone();
    for i in 1..=8 {
        let old = state.clone();
        let step = material.response(&set, &old).unwrap();
        assert_eq!(step.state.native, original.native);
        assert_eq!(step.plastic_increment, 0.0);
        assert_eq!(step.inelastic_dissipation_j_m3(), 0.0);
        assert_eq!(step.stored_energy_j_m3, 0.0);
        assert_eq!(step.flow_potential_j_m3, 0.0);
        assert_eq!(step.first_piola_pa, [0.0; 9]);
        near(
            step.state.resistance_strain,
            0.3 / 1.75_f64.powi(i),
            1e-16,
            2e-15,
        );
        material.validate_committed(&set, &step.state).unwrap();
        assert_eq!(material.stored_energy(&set, &step.state).unwrap(), 0.0);
        assert_eq!(old, state);
        state = step.state;
    }
    let reload = multiply(&coaxial(0.0, 0.08), &set);
    let without_dwell = material.response(&reload, &original).unwrap();
    let after_dwell = material.response(&reload, &state).unwrap();
    assert!(after_dwell.plastic_increment > without_dwell.plastic_increment);
    assert!(after_dwell.equivalent_kirchhoff_pa < without_dwell.equivalent_kirchhoff_pa);
    assert!(after_dwell.state.native.alpha > original.native.alpha);
    assert!(after_dwell.state.resistance_strain > state.resistance_strain);
}

#[test]
fn nonlinear_recovery_has_objective_history_and_consistent_potential_and_tangent() {
    let original = native(20e6, 2.5, 0.3)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: 0.015,
            exponent: 0.8,
            residual_rate_fraction: 0.17,
        })
        .unwrap();
    let material = FiniteRecoveringNortonCreep::new(original, 1.8).unwrap();
    let old = material
        .response(&coaxial(0.01, 0.07), &material.initial_state())
        .unwrap()
        .state;
    assert!(old.native.alpha > old.resistance_strain);
    let f = [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96];
    let r = material.response(&f, &old).unwrap();
    let multiplier = 0.17 + 0.83 * (1.0 + r.state.resistance_strain / 0.015).powf(-0.8);
    let expected_increment =
        0.3 * 2e6 / 20e6 * (r.equivalent_kirchhoff_pa / 2e6).powf(2.5) * multiplier;
    near(r.plastic_increment, expected_increment, 2e-13, 2e-11);
    assert_eq!(
        material.stress_and_potential(&f, &old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    for column in 0..9 {
        let h = 2e-6;
        let (mut plus, mut minus) = (f, f);
        plus[column] += h;
        minus[column] -= h;
        let (pp, wp) = material.stress_and_potential(&plus, &old).unwrap();
        let (pm, wm) = material.stress_and_potential(&minus, &old).unwrap();
        near((wp - wm) / (2.0 * h), r.first_piola_pa[column], 0.25, 3e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * h),
                r.tangent_pa[row][column],
                2.0,
                3e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                2e-12,
            );
        }
    }
    let (s, c) = 1.15_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = material.response(&multiply(&rotation, &f), &old).unwrap();
    near(
        rotated.state.resistance_strain,
        r.state.resistance_strain,
        2e-13,
        2e-11,
    );
    near(
        rotated.state.native.alpha,
        r.state.native.alpha,
        2e-13,
        2e-11,
    );
    near(
        rotated.flow_potential_j_m3,
        r.flow_potential_j_m3,
        2e-5,
        2e-10,
    );
    near(
        rotated.inelastic_dissipation_j_m3(),
        r.inelastic_dissipation_j_m3(),
        2e-5,
        2e-11,
    );
    for (actual, expected) in rotated
        .first_piola_pa
        .iter()
        .zip(multiply(&rotation, &r.first_piola_pa))
    {
        near(*actual, expected, 3e-5, 2e-11);
    }
    for (actual, expected) in rotated
        .state
        .native
        .inverse_plastic_cauchy_green
        .iter()
        .zip(r.state.native.inverse_plastic_cauchy_green)
    {
        near(*actual, expected, 2e-13, 2e-11);
    }
    // Pure volumetric stretch changes pressure only, not isochoric creep or resistance.
    let dilated = material.response(&f.map(|x| 1.04 * x), &old).unwrap();
    near(dilated.plastic_increment, r.plastic_increment, 2e-13, 2e-11);
    near(
        dilated.state.resistance_strain,
        r.state.resistance_strain,
        2e-13,
        2e-11,
    );
}

#[test]
fn disabled_primary_keeps_secondary_mechanics_while_retaining_resistance_memory() {
    for (exponent, fraction) in [(0.0, 0.1), (1.0, 1.0)] {
        let secondary = native(30e6, 2.0, 0.4);
        let selected = secondary
            .clone()
            .with_strain_hardening(StrainHardeningCreep {
                reference_strain: A,
                exponent,
                residual_rate_fraction: fraction,
            })
            .unwrap();
        let material = FiniteRecoveringNortonCreep::new(selected, 1.5).unwrap();
        let old = history(0.04, 0.08, 0.025);
        let f = coaxial(0.02, 0.12);
        let expected = secondary.response(&f, &old.native).unwrap();
        let actual = material.response(&f, &old).unwrap();
        same_native(&actual, &expected);
        near(
            actual.state.resistance_strain,
            (old.resistance_strain + actual.plastic_increment) / 1.6,
            0.0,
            2e-15,
        );
    }
}

#[test]
fn invalid_recovery_and_history_refuse_without_changing_accepted_state() {
    assert!(FiniteRecoveringNortonCreep::new(native(30e6, 1.0, 0.5), 1.0).is_err());
    for rate in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(FiniteRecoveringNortonCreep::new(primary(30e6, 1.0, 0.5), rate).is_err());
    }
    assert!(FiniteRecoveringNortonCreep::new(primary(30e6, 1.0, 2.0), f64::MAX).is_err());
    let unrepresentable = primary(30e6, 1.0, 1.0)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: f64::MAX,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap();
    assert!(FiniteRecoveringNortonCreep::new(unrepresentable, 1.0).is_err());
    let material = law(30e6, 1.0, 0.5, 1.5);
    let accepted = history(0.04, 0.08, 0.025);
    for resistance in [-0.01, f64::NAN, f64::INFINITY, 0.1] {
        let invalid = FiniteRecoveringCreepState {
            resistance_strain: resistance,
            ..accepted.clone()
        };
        assert!(material.validate_history(&invalid).is_err());
        assert!(material.response(&ID, &invalid).is_err());
        assert!(material.stress_and_potential(&ID, &invalid).is_err());
        assert!(material.stored_energy(&ID, &invalid).is_err());
        assert!(material.validate_committed(&ID, &invalid).is_err());
        assert!(
            material
                .endpoint_equivalent_stress_from(0.001, &invalid)
                .is_err()
        );
    }
    let mut invalid = material.initial_state();
    invalid.resistance_strain = f64::from_bits(1);
    assert!(material.validate_history(&invalid).is_err());
    assert!(
        material
            .response(&[-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], &accepted)
            .is_err()
    );
    assert_eq!(accepted, history(0.04, 0.08, 0.025));
}

fn combined_native(eta: f64, n: f64, dt: f64, linear: f64, saturation: f64) -> FiniteNortonJ2 {
    let backbone = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: linear,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 3.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
    .with_voce_hardening(saturation, 12.0)
    .unwrap();
    FiniteNortonJ2::new(
        backbone,
        FiniteOverstressParameters {
            viscosity_pa_s: eta,
            stress_scale_pa: 2e6,
            exponent: n,
        },
        dt,
        3.0,
    )
    .unwrap()
}

fn combined_primary(eta: f64, n: f64, dt: f64, linear: f64, saturation: f64) -> FiniteNortonJ2 {
    combined_native(eta, n, dt, linear, saturation)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap()
}

fn combined_history(
    net: f64,
    plastic: f64,
    creep: f64,
    resistance: f64,
) -> FiniteRecoveringCreepState {
    let mut state = history(net, plastic, resistance);
    state.native.creep_alpha = creep;
    state
}

#[test]
fn zero_recovery_keeps_exact_joint_primary_flow_and_both_native_permanent_paths() {
    let original = combined_primary(20e6, 2.5, 0.3, H, 15e6)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: 0.015,
            exponent: 0.8,
            residual_rate_fraction: 0.17,
        })
        .unwrap();
    let material = FiniteRecoveringNortonJ2::new(original.clone(), 0.0).unwrap();
    assert_eq!(material.norton_j2(), &original);
    assert_eq!(material.recovery_rate_s(), 0.0);
    assert_eq!(material.duration_s(), 0.3);
    assert_eq!(FiniteInelasticMaterial::duration_s(&material), Some(0.3));
    assert_eq!(
        FiniteInelasticMaterial::initial_state(&material),
        material.initial_state()
    );
    let mut state = material.initial_state();
    for f in [
        coaxial(0.01, 0.07),
        [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96],
        coaxial(-0.02, -0.12),
    ] {
        let old = state.clone();
        let expected = original.response(&f, &old.native).unwrap();
        let actual = FiniteInelasticMaterial::response(&material, &f, &old).unwrap();
        same_native(&actual, &expected);
        assert_eq!(
            actual.state.resistance_strain,
            actual.state.native.creep_alpha
        );
        assert_eq!(old, state);
        assert_eq!(
            FiniteInelasticMaterial::stress_and_potential(&material, &f, &old).unwrap(),
            (actual.first_piola_pa, actual.incremental_potential_j_m3)
        );
        FiniteInelasticMaterial::validate_committed(&material, &f, &actual.state, 1e-4).unwrap();
        near(
            FiniteInelasticMaterial::stored_energy(&material, &f, &actual.state).unwrap(),
            actual.stored_energy_j_m3,
            2e-5,
            2e-11,
        );
        state = actual.state;
    }
    assert!(state.native.alpha > 0.0);
    assert!(state.native.creep_alpha > 0.0);
}

#[test]
fn recovered_joint_return_matches_manufactured_perfect_linear_and_voce_endpoints() {
    let (eta, dt, k, mean) = (40e6, 0.2, 1.7, 0.015);
    let (old_plastic, old_creep, dp) = (0.02_f64, 0.035_f64, 0.025_f64);
    let old_net = old_plastic + old_creep;
    let old = combined_history(old_net, old_plastic, old_creep, 0.009);
    let divisor = 1.0 + k * dt;
    let a = A * divisor;
    for (linear, saturation) in [(0.0, 0.0), (H, 0.0), (H, 15e6)] {
        let material =
            FiniteRecoveringNortonJ2::new(combined_primary(eta, 1.0, dt, linear, saturation), k)
                .unwrap();
        let plastic = old_plastic + dp;
        let saturation_fraction = 1.0 - (-12.0 * plastic).exp();
        let q = Y + linear * plastic + saturation * saturation_fraction;
        // Author a yielded endpoint, solve its independent quadratic creep
        // relation, then construct the compatible finite deformation.
        let shift = a + old.resistance_strain;
        let drive = a * dt * q / eta;
        let dc = 2.0 * drive / (shift + (shift * shift + 4.0 * drive).sqrt());
        let total = old_net + dp + dc;
        let f = coaxial(mean, total + q / (3.0 * MU));
        let r = material.response(&f, &old).unwrap();
        let stored = 0.5 * K * (3.0 * mean).powi(2)
            + q * q / (6.0 * MU)
            + 0.5 * linear * plastic * plastic
            + saturation * (plastic - saturation_fraction / 12.0);
        let creep_potential =
            eta / dt * (0.5 * (1.0 + old.resistance_strain / a) * dc * dc + dc.powi(3) / (3.0 * a));
        let flow = Y * dp + creep_potential;
        near(r.plastic_increment, dp, 2e-13, 2e-11);
        near(r.creep_increment, dc, 2e-13, 2e-11);
        near(r.state.native.alpha, plastic, 2e-13, 2e-11);
        near(r.state.native.creep_alpha, old_creep + dc, 2e-13, 2e-11);
        near(
            r.state.resistance_strain,
            (old.resistance_strain + dc) / divisor,
            2e-13,
            2e-11,
        );
        near(r.equivalent_kirchhoff_pa, q, 3e-5, 2e-11);
        near(r.stored_energy_j_m3, stored, 2e-5, 2e-11);
        near(r.plastic_dissipation_j_m3, Y * dp, 2e-5, 2e-11);
        near(r.creep_dissipation_j_m3, q * dc, 2e-5, 2e-11);
        near(r.flow_potential_j_m3, flow, 2e-5, 2e-10);
        near(r.incremental_potential_j_m3, stored + flow, 2e-5, 2e-10);
        let expected = combined_history(total, plastic, old_creep + dc, 0.0);
        for (actual, expected) in r
            .state
            .native
            .inverse_plastic_cauchy_green
            .iter()
            .zip(expected.native.inverse_plastic_cauchy_green)
        {
            near(*actual, expected, 2e-13, 2e-11);
        }
        material.validate_committed(&f, &r.state, 1e-4).unwrap();
        near(
            material.stored_energy(&f, &r.state).unwrap(),
            stored,
            2e-5,
            2e-11,
        );
        assert_eq!(
            material
                .norton_j2()
                .strain_hardening()
                .unwrap()
                .reference_strain,
            A
        );
    }
    assert_eq!(
        old,
        combined_history(old_net, old_plastic, old_creep, 0.009)
    );
}

#[test]
fn stationary_resistance_allows_subyield_creep_without_plastic_hardening() {
    let (eta, k, q) = (20e6, 1.2, 0.8e6);
    let equilibrium = 2.0 * A * q / (eta * k) / (A + (A * A + 4.0 * A * q / (eta * k)).sqrt());
    let mut old = combined_history(0.04, 0.02, 0.08, equilibrium);
    let mut net = 0.04;
    let mut elapsed = 0.0;
    for dt in [0.1, 0.3, 0.2, 0.5] {
        let material =
            FiniteRecoveringNortonJ2::new(combined_primary(eta, 1.0, dt, H, 0.0), k).unwrap();
        let dc = dt * k * equilibrium;
        net += dc;
        elapsed += dt;
        let f = coaxial(0.02, net + q / (3.0 * MU));
        let step = material.response(&f, &old).unwrap();
        assert_eq!(step.plastic_increment, 0.0);
        assert_eq!(step.plastic_dissipation_j_m3, 0.0);
        assert_eq!(step.state.native.alpha, 0.02);
        near(step.creep_increment, dc, 2e-13, 2e-10);
        near(
            step.state.native.creep_alpha,
            0.08 + elapsed * k * equilibrium,
            2e-13,
            2e-10,
        );
        near(step.state.resistance_strain, equilibrium, 2e-13, 2e-10);
        near(step.equivalent_kirchhoff_pa, q, 3e-5, 2e-10);
        near(step.creep_dissipation_j_m3, q * dc, 2e-5, 2e-10);
        old = step.state;
    }
}

#[test]
fn unloaded_joint_history_keeps_plastic_storage_and_permanent_set_during_recovery() {
    let material =
        FiniteRecoveringNortonJ2::new(combined_primary(30e6, 2.0, 0.5, H, 0.0), 1.5).unwrap();
    let set = [2.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0];
    let original = FiniteRecoveringCreepState {
        native: FiniteJ2State {
            inverse_plastic_cauchy_green: [0.25, 4.0, 1.0, 0.0, 0.0, 0.0],
            alpha: 0.2,
            creep_alpha: 1.0,
        },
        // Larger than plastic alpha: only permanent creep bounds resistance.
        resistance_strain: 0.3,
    };
    let stored = material.stored_energy(&set, &original).unwrap();
    near(stored, 0.5 * H * original.native.alpha.powi(2), 1e-8, 2e-14);
    let mut state = original.clone();
    for i in 1..=8 {
        let old = state.clone();
        let step = material.response(&set, &old).unwrap();
        assert_eq!(step.state.native, original.native);
        assert_eq!(step.plastic_increment, 0.0);
        assert_eq!(step.creep_increment, 0.0);
        assert_eq!(step.inelastic_dissipation_j_m3(), 0.0);
        assert_eq!(step.stored_energy_j_m3, stored);
        assert_eq!(step.flow_potential_j_m3, 0.0);
        assert_eq!(step.first_piola_pa, [0.0; 9]);
        near(
            step.state.resistance_strain,
            0.3 / 1.75_f64.powi(i),
            1e-16,
            2e-15,
        );
        material
            .validate_committed(&set, &step.state, 1e-4)
            .unwrap();
        assert_eq!(material.stored_energy(&set, &step.state).unwrap(), stored);
        assert_eq!(old, state);
        state = step.state;
    }
    let reload = multiply(&coaxial(0.0, 0.08), &set);
    let without_dwell = material.response(&reload, &original).unwrap();
    let after_dwell = material.response(&reload, &state).unwrap();
    assert!(after_dwell.creep_increment > without_dwell.creep_increment);
    assert!(after_dwell.plastic_increment < without_dwell.plastic_increment);
    assert!(after_dwell.equivalent_kirchhoff_pa < without_dwell.equivalent_kirchhoff_pa);
    assert!(after_dwell.state.native.creep_alpha > original.native.creep_alpha);
}

#[test]
fn nonlinear_joint_recovery_retains_objectivity_and_full_potential_tangent() {
    let native = combined_primary(20e6, 2.5, 0.3, H, 15e6)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: 0.015,
            exponent: 0.8,
            residual_rate_fraction: 0.17,
        })
        .unwrap();
    let material = FiniteRecoveringNortonJ2::new(native, 1.8).unwrap();
    let old = material
        .response(&coaxial(0.01, 0.07), &material.initial_state())
        .unwrap()
        .state;
    assert!(old.native.creep_alpha > old.resistance_strain);
    let f = [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96];
    let r = material.response(&f, &old).unwrap();
    assert!(r.plastic_increment > 0.0);
    let multiplier = 0.17 + 0.83 * (1.0 + r.state.resistance_strain / 0.015).powf(-0.8);
    let expected_increment =
        0.3 * 2e6 / 20e6 * (r.equivalent_kirchhoff_pa / 2e6).powf(2.5) * multiplier;
    near(r.creep_increment, expected_increment, 2e-13, 2e-11);
    let plastic = r.state.native.alpha;
    near(
        r.equivalent_kirchhoff_pa,
        Y + H * plastic + 15e6 * (1.0 - (-12.0 * plastic).exp()),
        3e-5,
        2e-11,
    );
    assert_eq!(
        material.stress_and_potential(&f, &old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    for column in 0..9 {
        let h = 2e-6;
        let (mut plus, mut minus) = (f, f);
        plus[column] += h;
        minus[column] -= h;
        let (pp, wp) = material.stress_and_potential(&plus, &old).unwrap();
        let (pm, wm) = material.stress_and_potential(&minus, &old).unwrap();
        near((wp - wm) / (2.0 * h), r.first_piola_pa[column], 0.25, 3e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * h),
                r.tangent_pa[row][column],
                2.0,
                3e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                2e-12,
            );
        }
    }
    let (s, c) = 1.15_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = material.response(&multiply(&rotation, &f), &old).unwrap();
    near(
        rotated.state.resistance_strain,
        r.state.resistance_strain,
        2e-13,
        2e-11,
    );
    near(
        rotated.state.native.alpha,
        r.state.native.alpha,
        2e-13,
        2e-11,
    );
    near(
        rotated.state.native.creep_alpha,
        r.state.native.creep_alpha,
        2e-13,
        2e-11,
    );
    near(
        rotated.flow_potential_j_m3,
        r.flow_potential_j_m3,
        2e-5,
        2e-10,
    );
    near(
        rotated.inelastic_dissipation_j_m3(),
        r.inelastic_dissipation_j_m3(),
        2e-5,
        2e-11,
    );
    for (actual, expected) in rotated
        .first_piola_pa
        .iter()
        .zip(multiply(&rotation, &r.first_piola_pa))
    {
        near(*actual, expected, 3e-5, 2e-11);
    }
    for (actual, expected) in rotated
        .state
        .native
        .inverse_plastic_cauchy_green
        .iter()
        .zip(r.state.native.inverse_plastic_cauchy_green)
    {
        near(*actual, expected, 2e-13, 2e-11);
    }
    let dilated = material.response(&f.map(|x| 1.04 * x), &old).unwrap();
    near(dilated.plastic_increment, r.plastic_increment, 2e-13, 2e-11);
    near(dilated.creep_increment, r.creep_increment, 2e-13, 2e-11);
    near(
        dilated.state.resistance_strain,
        r.state.resistance_strain,
        2e-13,
        2e-11,
    );
}

#[test]
fn joint_secondary_limit_and_history_admission_keep_resistance_separate_from_plasticity() {
    let secondary = combined_native(30e6, 2.0, 0.4, H, 0.0);
    assert!(FiniteRecoveringNortonJ2::new(secondary.clone(), 1.5).is_err());
    let old = combined_history(0.04, 0.02, 0.08, 0.025);
    let f = coaxial(0.02, 0.12);
    for (exponent, fraction) in [(0.0, 0.1), (1.0, 1.0)] {
        let selected = secondary
            .clone()
            .with_strain_hardening(StrainHardeningCreep {
                reference_strain: A,
                exponent,
                residual_rate_fraction: fraction,
            })
            .unwrap();
        let material = FiniteRecoveringNortonJ2::new(selected, 1.5).unwrap();
        let expected = secondary.response(&f, &old.native).unwrap();
        let actual = material.response(&f, &old).unwrap();
        same_native(&actual, &expected);
        assert!(actual.plastic_increment > 0.0);
        near(
            actual.state.resistance_strain,
            (old.resistance_strain + actual.creep_increment) / 1.6,
            0.0,
            2e-15,
        );
        for resistance in [-0.01, f64::NAN, f64::INFINITY, 0.1] {
            let invalid = FiniteRecoveringCreepState {
                resistance_strain: resistance,
                ..old.clone()
            };
            assert!(material.validate_history(&invalid).is_err());
            assert!(material.response(&f, &invalid).is_err());
            assert!(material.stress_and_potential(&f, &invalid).is_err());
            assert!(material.stored_energy(&f, &invalid).is_err());
            assert!(material.validate_committed(&f, &invalid, 1e-4).is_err());
        }
        // Plastic accumulation cannot excuse excessive creep resistance.
        let invalid = combined_history(0.04, 0.2, 0.01, 0.02);
        assert!(material.validate_history(&invalid).is_err());
        for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                FiniteInelasticMaterial::validate_committed(
                    &material,
                    &f,
                    &actual.state,
                    tolerance
                )
                .is_err()
            );
        }
    }
    assert_eq!(old, combined_history(0.04, 0.02, 0.08, 0.025));
}
