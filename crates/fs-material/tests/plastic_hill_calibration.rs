//! Independent sheet-flow checks for explicit Lankford calibration.
//! Synthetic r-values exercise the constitutive contract, not alloy validation.

use fs_material::IsotropicElastic;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};

const R: [f64; 3] = [1.8, 1.1, 2.2];
const SHEAR: [f64; 2] = [0.8, 1.3];
const E: f64 = 210e9;
const NU: f64 = 0.3;
const Y: f64 = 250e6;
const H: f64 = 1.2e9;

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn equivalent(c: Hill48Coefficients, s: &[f64; 6]) -> f64 {
    (c.f * (s[1] - s[2]).powi(2)
        + c.g * (s[2] - s[0]).powi(2)
        + c.h * (s[0] - s[1]).powi(2)
        + 2.0 * c.l * s[4].powi(2)
        + 2.0 * c.m * s[5].powi(2)
        + 2.0 * c.n * s[3].powi(2))
    .sqrt()
}

#[test]
fn actual_associated_flow_recovers_all_three_measured_r_values_in_both_load_signs() {
    let c = Hill48Coefficients::from_planar_r_values(R, SHEAR).unwrap();
    let law = Hill48::new(IsotropicElastic::new(E, NU, 0.04).unwrap(), Y, H, c).unwrap();
    let gamma = 0.004;
    let q = Y + H * gamma;
    let diagonal = core::f64::consts::FRAC_1_SQRT_2;
    for ((cosine, sine), expected_r) in [(1.0, 0.0), (diagonal, diagonal), (0.0, 1.0)]
        .into_iter()
        .zip(R)
    {
        for sign in [-1.0, 1.0] {
            let unit = [cosine * cosine, sine * sine, 0.0, cosine * sine, 0.0, 0.0];
            let axial = sign * q / equivalent(c, &unit);
            let stress = unit.map(|v| v * axial);
            // Independent gradient in TENSOR shear components, and an exact
            // backward-Euler endpoint constructed from stress compliance.
            let normal = [
                ((c.g + c.h) * stress[0] - c.h * stress[1]) / q,
                ((c.f + c.h) * stress[1] - c.h * stress[0]) / q,
                (-c.g * stress[0] - c.f * stress[1]) / q,
                c.n * stress[3] / q,
                0.0,
                0.0,
            ];
            let plastic = normal.map(|n| gamma * n);
            let trace = stress[0] + stress[1] + stress[2];
            let strain = core::array::from_fn(|i| {
                plastic[i]
                    + if i < 3 {
                        ((1.0 + NU) * stress[i] - NU * trace) / E
                    } else {
                        (1.0 + NU) * stress[i] / E
                    }
            });
            let response = law.response(&strain, &law.initial_state()).unwrap();
            near(response.state.alpha, gamma, 2e-15, 1e-12);
            for i in 0..6 {
                near(response.stress_pa[i], stress[i], 1e-4, 2e-12);
                near(response.state.plastic_strain[i], plastic[i], 3e-15, 2e-12);
            }
            let p = response.state.plastic_strain;
            let width = sine * sine * p[0] + cosine * cosine * p[1] - 2.0 * sine * cosine * p[3];
            near(width / p[2], expected_r, 2e-12, 2e-12);
            near(p[0] + p[1] + p[2], 0.0, 3e-18, 0.0);
        }
    }
}

#[test]
fn rolling_reference_yield_and_independent_transverse_shear_ratios_are_preserved() {
    let c = Hill48Coefficients::from_planar_r_values(R, SHEAR).unwrap();
    near(c.g + c.h, 1.0, 3e-16, 0.0);
    let law = Hill48::new(IsotropicElastic::new(E, NU, 0.04).unwrap(), Y, H, c).unwrap();
    near(
        law.equivalent_stress_pa(&[Y, 0.0, 0.0, 0.0, 0.0, 0.0])
            .unwrap(),
        Y,
        1e-6,
        2e-14,
    );
    for (component, ratio) in [4, 5].into_iter().zip(SHEAR) {
        let mut stress = [0.0; 6];
        stress[component] = ratio * Y / 3.0_f64.sqrt();
        near(law.equivalent_stress_pa(&stress).unwrap(), Y, 1e-6, 2e-14);
    }
    // Planar measurements do not determine L or M: changing only the supplied
    // transverse-shear data leaves every in-plane coefficient unchanged.
    let changed = Hill48Coefficients::from_planar_r_values(R, [1.4, 0.7]).unwrap();
    assert_eq!(
        (c.f, c.g, c.h, c.n),
        (changed.f, changed.g, changed.h, changed.n)
    );
    assert_ne!((c.l, c.m), (changed.l, changed.m));
}

#[test]
fn isotropic_data_recover_mises_without_losing_small_positive_r_values() {
    assert_eq!(
        Hill48Coefficients::from_planar_r_values([1.0; 3], [1.0; 2]).unwrap(),
        Hill48Coefficients::isotropic()
    );
    let tiny = Hill48Coefficients::from_planar_r_values([1e-16, 1.0, 1.0], [1.0; 2]).unwrap();
    assert!(tiny.h > 0.0);
    near(tiny.h / tiny.g, 1e-16, 0.0, 1e-14);
}

#[test]
fn incomplete_or_unrepresentable_planar_calibration_refuses() {
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for index in 0..3 {
            let mut planar = R;
            planar[index] = invalid;
            assert!(Hill48Coefficients::from_planar_r_values(planar, SHEAR).is_err());
        }
        for index in 0..2 {
            let mut shear = SHEAR;
            shear[index] = invalid;
            assert!(Hill48Coefficients::from_planar_r_values(R, shear).is_err());
        }
    }
    for unrepresentable in [1e-300, 1e300] {
        assert!(Hill48Coefficients::from_planar_r_values(R, [unrepresentable, 1.0]).is_err());
    }
    assert!(Hill48Coefficients::from_planar_r_values([1e-300, 1.0, 1e100], [1.0; 2]).is_err());
    assert!(Hill48Coefficients::from_planar_r_values([1.0, 1e-300, 1.0], [1.0; 2]).is_err());
}
