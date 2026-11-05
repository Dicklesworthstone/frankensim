//! Finite porous compression on the native logarithmic metric and return map.
//! Synthetic coefficients verify mathematics; these are not foam calibrations.

use fs_material::SmallStrainLaw;
use fs_material::elastic::IsotropicElastic;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::porous::PorousPlasticity;

const MU: f64 = 80.0;
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;

fn parameters() -> FinitePorousParameters {
    FinitePorousParameters {
        shear_modulus_pa: MU,
        bulk_modulus_pa: K,
        yield_stress_pa: Y,
        hardening_modulus_pa: H,
        pressure_sensitivity: A,
    }
}

fn limits() -> FinitePorousLimits {
    FinitePorousLimits {
        minimum_jacobian: 0.2,
        maximum_elastic_log_strain: 0.5,
        maximum_plastic_log_strain: 2.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}

fn law() -> FinitePorous {
    FinitePorous::new(parameters(), limits()).unwrap()
}

fn diagonal(x: f64) -> [f64; 9] {
    [x, 0.0, 0.0, 0.0, x, 0.0, 0.0, 0.0, x]
}

fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| a[3 * (i % 3) + i / 3])
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| (0..3).map(|k| a[3 * (i / 3) + k] * b[3 * k + i % 3]).sum())
}

fn expand(v: &[f64; 6]) -> [f64; 9] {
    [v[0], v[3], v[5], v[3], v[1], v[4], v[5], v[4], v[2]]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

fn near_array<const N: usize>(actual: &[f64; N], expected: &[f64; N], relative: f64) {
    let scale = expected.iter().fold(1.0_f64, |s, x| s.max(x.abs()));
    for i in 0..N {
        near(actual[i], expected[i], relative * scale);
    }
}

#[test]
fn hydrostatic_compression_retains_finite_volume_loss_storage_and_yield_heat() {
    let material = law();
    let initial = material.initial_state();
    let saved = initial.clone();
    let stretch = 0.8_f64;
    let beta = A / (1.0 + A * A / 9.0).sqrt();
    let alpha = (-K * beta * (3.0 * stretch.ln()) - Y) / (H + K * beta * beta);
    let mean_kirchhoff = -(Y + H * alpha) / beta;
    let plastic_j = (-beta * alpha).exp();
    let storage = 0.5 * mean_kirchhoff.powi(2) / K + 0.5 * H * alpha * alpha;
    let heat = Y * alpha;
    let f = diagonal(stretch);
    let response = material.response(&f, &initial).unwrap();
    assert_eq!(
        initial, saved,
        "a trial must not advance caller-owned history"
    );
    assert!(alpha > 0.5 && plastic_j < 0.55);
    near(response.state.alpha, alpha, 2e-12);
    near(response.plastic_increment, alpha, 2e-12);
    near(response.jacobian, stretch.powi(3), 2e-14);
    near(response.plastic_jacobian, plastic_j, 2e-12);
    near(
        material.plastic_jacobian(&response.state).unwrap(),
        plastic_j,
        2e-12,
    );
    near_array(&response.kirchhoff_pa, &diagonal(mean_kirchhoff), 2e-11);
    near_array(
        &response.first_piola_pa,
        &diagonal(mean_kirchhoff / stretch),
        2e-11,
    );
    near(response.equivalent_kirchhoff_pa, Y + H * alpha, 2e-10);
    near(response.yield_residual_pa, 0.0, 2e-10);
    near(response.stored_energy_j_m3, storage, 2e-11);
    near(response.plastic_dissipation_j_m3, heat, 2e-12);
    near(response.incremental_potential_j_m3, storage + heat, 3e-11);
    near(
        material.stored_energy(&f, &response.state).unwrap(),
        storage,
        3e-11,
    );
    material
        .validate_committed(&f, &response.state, 2e-9)
        .unwrap();
    let inverse_plastic_stretch_squared = plastic_j.powf(-2.0 / 3.0);
    near_array(
        &response.state.inverse_plastic_cauchy_green,
        &[
            inverse_plastic_stretch_squared,
            inverse_plastic_stretch_squared,
            inverse_plastic_stretch_squared,
            0.0,
            0.0,
            0.0,
        ],
        2e-11,
    );

    // Full unloading restores zero elastic stress but retains the plastic volume.
    let stress_free_stretch = (-beta * alpha / 3.0).exp();
    let unloaded = material
        .response(&diagonal(stress_free_stretch), &response.state)
        .unwrap();
    assert_eq!(unloaded.state, response.state);
    near(unloaded.plastic_increment, 0.0, 0.0);
    near(unloaded.plastic_dissipation_j_m3, 0.0, 0.0);
    near_array(&unloaded.first_piola_pa, &[0.0; 9], 2e-10);
    near(unloaded.stored_energy_j_m3, 0.5 * H * alpha * alpha, 2e-10);
    let reload_stretch = (stretch * stress_free_stretch).sqrt();
    let reloaded = material
        .response(&diagonal(reload_stretch), &unloaded.state)
        .unwrap();
    assert_eq!(reloaded.state, unloaded.state);
    near_array(
        &reloaded.kirchhoff_pa,
        &diagonal(0.5 * mean_kirchhoff),
        3e-10,
    );
}

#[test]
fn all_nine_tangent_columns_and_potential_gradient_match_after_noncoaxial_history() {
    let material = law();
    let compacted = material
        .response(&diagonal(0.8), &material.initial_state())
        .unwrap()
        .state;
    let prior_f = [0.89, 0.11, -0.02, 0.03, 0.93, 0.06, -0.01, 0.04, 0.84];
    let old = material.response(&prior_f, &compacted).unwrap().state;
    let saved = old.clone();
    let f = [1.03, -0.06, 0.04, 0.09, 0.91, 0.03, 0.02, -0.01, 0.87];
    let response = material.response(&f, &old).unwrap();
    assert!(response.plastic_increment > 0.1);
    let scale = response
        .tangent_pa
        .iter()
        .flatten()
        .fold(1.0_f64, |s, x| s.max(x.abs()));
    let h = 2e-6;
    for column in 0..9 {
        let mut plus = f;
        let mut minus = f;
        plus[column] += h;
        minus[column] -= h;
        let (p_plus, potential_plus) = material.stress_and_potential(&plus, &old).unwrap();
        let (p_minus, potential_minus) = material.stress_and_potential(&minus, &old).unwrap();
        for row in 0..9 {
            let derivative = (p_plus[row] - p_minus[row]) / (2.0 * h);
            near(response.tangent_pa[row][column], derivative, 2e-6 * scale);
            near(
                response.tangent_pa[row][column],
                response.tangent_pa[column][row],
                2e-8 * scale,
            );
        }
        near(
            (potential_plus - potential_minus) / (2.0 * h),
            response.first_piola_pa[column],
            2e-6,
        );
    }
    assert_eq!(old, saved);
}

#[test]
fn superposed_rotation_is_objective_and_accepted_holds_do_not_create_plastic_flow() {
    let material = law();
    let compacted = material
        .response(&diagonal(0.8), &material.initial_state())
        .unwrap()
        .state;
    let accepted_f = [0.89, 0.11, -0.02, 0.03, 0.93, 0.06, -0.01, 0.04, 0.84];
    let accepted = material.response(&accepted_f, &compacted).unwrap();
    let rotation = [
        -2.0 / 3.0,
        2.0 / 15.0,
        11.0 / 15.0,
        2.0 / 3.0,
        -1.0 / 3.0,
        2.0 / 3.0,
        1.0 / 3.0,
        14.0 / 15.0,
        2.0 / 15.0,
    ];
    let held = material.response(&accepted_f, &accepted.state).unwrap();
    let rotated_hold = material
        .response(&multiply(&rotation, &accepted_f), &accepted.state)
        .unwrap();
    for response in [&held, &rotated_hold] {
        assert_eq!(response.state, accepted.state);
        near(response.plastic_increment, 0.0, 0.0);
        near(response.plastic_dissipation_j_m3, 0.0, 0.0);
        near(
            response.stored_energy_j_m3,
            accepted.stored_energy_j_m3,
            2e-9,
        );
    }
    near_array(
        &rotated_hold.first_piola_pa,
        &multiply(&rotation, &held.first_piola_pa),
        2e-9,
    );

    let f = [1.03, -0.06, 0.04, 0.09, 0.91, 0.03, 0.02, -0.01, 0.87];
    let original = material.response(&f, &accepted.state).unwrap();
    let rotated = material
        .response(&multiply(&rotation, &f), &accepted.state)
        .unwrap();
    near_array(
        &rotated.first_piola_pa,
        &multiply(&rotation, &original.first_piola_pa),
        2e-9,
    );
    near_array(
        &rotated.kirchhoff_pa,
        &multiply(
            &multiply(&rotation, &original.kirchhoff_pa),
            &transpose(&rotation),
        ),
        2e-9,
    );
    near_array(
        &rotated.state.inverse_plastic_cauchy_green,
        &original.state.inverse_plastic_cauchy_green,
        2e-9,
    );
    near(rotated.state.alpha, original.state.alpha, 2e-10);
    near(rotated.plastic_jacobian, original.plastic_jacobian, 2e-10);
    near(
        rotated.incremental_potential_j_m3,
        original.incremental_potential_j_m3,
        2e-9,
    );
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let expected: f64 = (0..3)
                        .flat_map(|m| (0..3).map(move |n| (m, n)))
                        .map(|(m, n)| {
                            rotation[3 * i + m]
                                * rotation[3 * k + n]
                                * original.tangent_pa[3 * m + j][3 * n + l]
                        })
                        .sum();
                    near(rotated.tangent_pa[3 * i + j][3 * k + l], expected, 2e-7);
                }
            }
        }
    }
}

#[test]
fn zero_pressure_sensitivity_recovers_native_finite_j2_through_reversal() {
    let mut p = parameters();
    p.pressure_sensitivity = 0.0;
    let porous = FinitePorous::new(p, limits()).unwrap();
    let j2 = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    let mut porous_state = porous.initial_state();
    let mut j2_state = j2.initial_state();
    for f in [
        diagonal(0.9),
        [0.87, 0.1, 0.0, 0.0, 1.02, 0.0, 0.0, 0.0, 0.98],
        [1.04, -0.06, 0.02, 0.05, 0.95, 0.03, 0.0, 0.0, 0.91],
    ] {
        let actual = porous.response(&f, &porous_state).unwrap();
        let expected = j2.response(&f, &j2_state).unwrap();
        near_array(&actual.first_piola_pa, &expected.first_piola_pa, 3e-9);
        near_array(&actual.kirchhoff_pa, &expected.kirchhoff_pa, 3e-9);
        for row in 0..9 {
            near_array(&actual.tangent_pa[row], &expected.tangent_pa[row], 5e-8);
        }
        near_array(
            &actual.state.inverse_plastic_cauchy_green,
            &expected.state.inverse_plastic_cauchy_green,
            3e-9,
        );
        near(actual.state.alpha, expected.state.alpha, 2e-10);
        near(actual.stored_energy_j_m3, expected.stored_energy_j_m3, 3e-9);
        near(
            actual.incremental_potential_j_m3,
            expected.incremental_potential_j_m3,
            3e-9,
        );
        near(
            actual.plastic_dissipation_j_m3,
            expected.plastic_dissipation_j_m3,
            2e-9,
        );
        near(actual.plastic_jacobian, 1.0, 3e-10);
        porous_state = actual.state;
        j2_state = expected.state;
    }
}

#[test]
fn logarithmic_return_matches_small_strain_owner_and_preserves_stress_unit_scaling() {
    let material = law();
    let e = [-0.012_f64, -0.012, -0.004, 0.002, 0.0, 0.0];
    // Exact exponential of this symmetric logarithmic strain, independent of
    // the finite owner's eigensolver and metric reconstruction.
    let axial = e[0].exp();
    let f = [
        axial * e[3].cosh(),
        axial * e[3].sinh(),
        0.0,
        axial * e[3].sinh(),
        axial * e[3].cosh(),
        0.0,
        0.0,
        0.0,
        e[2].exp(),
    ];
    let small = PorousPlasticity::new(
        IsotropicElastic::new(
            9.0 * K * MU / (3.0 * K + MU),
            (3.0 * K - 2.0 * MU) / (2.0 * (3.0 * K + MU)),
            0.05,
        )
        .unwrap(),
        Y,
        H,
        A,
    )
    .unwrap();
    let reference = small.response(&e, &small.initial_state()).unwrap();
    let initial = material.initial_state();
    let actual = material.response(&f, &initial).unwrap();
    assert!(actual.plastic_increment > 0.0);
    near_array(&actual.kirchhoff_pa, &expand(&reference.stress_pa), 3e-10);
    near(actual.state.alpha, reference.state.alpha, 3e-12);
    near(
        actual.plastic_jacobian,
        reference.state.plastic_strain[..3]
            .iter()
            .sum::<f64>()
            .exp(),
        3e-12,
    );
    near(
        actual.stored_energy_j_m3,
        small
            .stored_energy(&e, &reference.state)
            .unwrap()
            .total_j_m3(),
        3e-11,
    );
    let frozen = material.frozen(&initial).unwrap();
    let (stress, potential) = frozen.stress_and_potential(&f).unwrap();
    near_array(&stress, &actual.first_piola_pa, 1e-12);
    near(potential, actual.incremental_potential_j_m3, 1e-12);
    near(frozen.energy_density(&f), potential, 1e-12);
    near_array(&frozen.first_piola(&f).unwrap(), &stress, 1e-12);
    for (row, expected) in frozen
        .consistent_tangent(&f)
        .unwrap()
        .iter()
        .zip(&actual.tangent_pa)
    {
        near_array(row, expected, 1e-12);
    }
    for factor in [1e-6, 1e6] {
        let mut p = parameters();
        p.shear_modulus_pa *= factor;
        p.bulk_modulus_pa *= factor;
        p.yield_stress_pa *= factor;
        p.hardening_modulus_pa *= factor;
        let scaled = FinitePorous::new(p, limits()).unwrap();
        let result = scaled.response(&f, &scaled.initial_state()).unwrap();
        near_array(
            &result.first_piola_pa.map(|x| x / factor),
            &actual.first_piola_pa,
            3e-9,
        );
        for row in 0..9 {
            near_array(
                &result.tangent_pa[row].map(|x| x / factor),
                &actual.tangent_pa[row],
                3e-9,
            );
        }
        near_array(
            &result.state.inverse_plastic_cauchy_green,
            &actual.state.inverse_plastic_cauchy_green,
            3e-10,
        );
        near(result.state.alpha, actual.state.alpha, 3e-11);
        near(
            result.stored_energy_j_m3 / factor,
            actual.stored_energy_j_m3,
            3e-10,
        );
        near(
            result.plastic_dissipation_j_m3 / factor,
            actual.plastic_dissipation_j_m3,
            3e-10,
        );
    }
}

#[test]
fn invalid_coefficients_metrics_bounds_and_committed_states_refuse_without_repairs() {
    for index in 0..6 {
        let mut p = parameters();
        match index {
            0 => p.shear_modulus_pa = 0.0,
            1 => p.bulk_modulus_pa = f64::INFINITY,
            2 => p.yield_stress_pa = 0.0,
            3 => p.hardening_modulus_pa = -1.0,
            4 => p.pressure_sensitivity = -0.1,
            _ => p.pressure_sensitivity = f64::NAN,
        }
        assert!(FinitePorous::new(p, limits()).is_err());
    }
    for index in 0..5 {
        let mut domain = limits();
        match index {
            0 => domain.minimum_jacobian = 0.0,
            1 => domain.minimum_jacobian = 1.0,
            2 => domain.maximum_elastic_log_strain = 0.0,
            3 => domain.maximum_plastic_log_strain = -1.0,
            _ => domain.maximum_accumulated_plastic_strain = f64::INFINITY,
        }
        assert!(FinitePorous::new(parameters(), domain).is_err());
    }
    let material = law();
    for index in 0..7 {
        let mut state = material.initial_state();
        match index {
            0 => state.alpha = -1.0,
            1 => state.alpha = f64::NAN,
            2 => {
                state.alpha = 1.0;
                state.inverse_plastic_cauchy_green[0] = 0.0;
            }
            3 => {
                state.alpha = 1.0;
                state.inverse_plastic_cauchy_green[0] = f64::NAN;
            }
            4 => state.inverse_plastic_cauchy_green = [2.0, 2.0, 2.0, 0.0, 0.0, 0.0],
            5 => state.alpha = 4.0,
            _ => {
                state.alpha = 1e-8;
                state.inverse_plastic_cauchy_green[0] = 2.0;
            }
        }
        assert!(material.validate_history(&state).is_err());
        assert!(material.frozen(&state).is_err());
    }
    let initial = material.initial_state();
    let saved = initial.clone();
    for f in [[0.0; 9], diagonal(0.1), diagonal(-1.0), [f64::NAN; 9]] {
        assert!(material.response(&f, &initial).is_err());
    }
    assert_eq!(initial, saved);
    let compressed = material.response(&diagonal(0.8), &initial).unwrap();
    assert!(
        material
            .validate_committed(&diagonal(1.0), &compressed.state, 1e-9)
            .is_err()
    );
    assert!(
        material
            .validate_committed(&diagonal(0.8), &compressed.state, f64::NAN)
            .is_err()
    );
    let mut p = parameters();
    p.pressure_sensitivity = 0.0;
    let j2_limit = FinitePorous::new(p, limits()).unwrap();
    assert!(j2_limit.validate_history(&compressed.state).is_err());
    for index in 0..3 {
        let mut domain = limits();
        match index {
            0 => domain.maximum_elastic_log_strain = 0.001,
            1 => domain.maximum_plastic_log_strain = 0.01,
            _ => domain.maximum_accumulated_plastic_strain = 0.01,
        }
        let bounded = FinitePorous::new(parameters(), domain).unwrap();
        assert!(
            bounded
                .response(&diagonal(0.8), &bounded.initial_state())
                .is_err()
        );
    }
}
