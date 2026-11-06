//! G0/G1/G3: independent Hill returns, finite stress scaling and lazy tangents.
//! These numerical coefficients are synthetic, not a measured alloy calibration.

use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients, HillState};
use fs_material::tensor::{contract, rotate, rotation, trace};
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};

const E: f64 = 210_000.0;
const NU: f64 = 0.29;
const Y: f64 = 180.0;
const HARDENING: f64 = 2400.0;
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn coefficients() -> Hill48Coefficients {
    Hill48Coefficients {
        f: 0.8,
        g: 1.2,
        h: 0.4,
        l: 1.7,
        m: 0.9,
        n: 2.1,
    }
}
fn law(hardening: f64) -> Hill48 {
    Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        Y,
        hardening,
        coefficients(),
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        actual.is_finite()
            && (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn finite_equivalent_stress_avoids_overflow_in_deviatoric_coordinates() {
    let law = law(HARDENING);
    for sign in [-1.0, 1.0] {
        let stress = [sign * 1e308, sign * 1e308, 0.0, 0.0, 0.0, 0.0];
        // The sum of normal differences overflows before multiplication by
        // 1/sqrt(6), but the original Hill quadratic has a finite square root.
        near(
            law.equivalent_stress_pa(&stress).unwrap() / 1e308,
            (coefficients().f + coefficients().g).sqrt(),
            2e-14,
            2e-14,
        );
    }
    assert_eq!(
        law.equivalent_stress_pa(&[1e308, 1e308, 1e308, 0.0, 0.0, 0.0])
            .unwrap(),
        0.0
    );
    let shear = Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        Y,
        HARDENING,
        Hill48Coefficients {
            n: 0.25,
            ..coefficients()
        },
    )
    .unwrap();
    // Mandel shear sqrt(2)*sigma_xy itself overflows, although q is finite.
    near(
        shear
            .equivalent_stress_pa(&[0.0, 0.0, 0.0, 1.5e308, 0.0, 0.0])
            .unwrap()
            / 1.5e308,
        0.5_f64.sqrt(),
        2e-14,
        2e-14,
    );
    assert!(
        law.equivalent_stress_pa(&[0.0, 0.0, 0.0, 1.5e308, 0.0, 0.0])
            .is_err()
    );
    assert!(law.equivalent_stress_pa(&[f64::NAN; 6]).is_err());
    assert!(law.equivalent_stress_pa(&[f64::INFINITY; 6]).is_err());
}

#[test]
fn independent_noncoaxial_endpoint_recovers_flow_and_work_conjugacy() {
    let law = law(HARDENING);
    let c = coefficients();
    let gamma = 0.0031;
    let direction: Voigt = [-0.7, 0.4, 0.3, 0.2, -0.3, 0.5];
    // Independent six-component quadratic, without the native modal return.
    let phi = (c.f * (direction[1] - direction[2]).powi(2)
        + c.g * (direction[2] - direction[0]).powi(2)
        + c.h * (direction[0] - direction[1]).powi(2)
        + 2.0 * c.n * direction[3].powi(2)
        + 2.0 * c.l * direction[4].powi(2)
        + 2.0 * c.m * direction[5].powi(2))
    .sqrt();
    let q = Y + HARDENING * gamma;
    let stress = direction.map(|value| value * q / phi);
    // Tensor shear already occurs twice in sigma:epsilon. No engineering
    // factor of two belongs in the plastic tensor's shear flow components.
    let flow = [
        (c.g * (stress[0] - stress[2]) + c.h * (stress[0] - stress[1])) / q,
        (c.f * (stress[1] - stress[2]) + c.h * (stress[1] - stress[0])) / q,
        (c.g * (stress[2] - stress[0]) + c.f * (stress[2] - stress[1])) / q,
        c.n * stress[3] / q,
        c.l * stress[4] / q,
        c.m * stress[5] / q,
    ];
    let plastic = flow.map(|value| gamma * value);
    let elastic: Voigt = core::array::from_fn(|i| {
        if i < 3 {
            ((1.0 + NU) * stress[i] - NU * trace(&stress)) / E
        } else {
            (1.0 + NU) * stress[i] / E
        }
    });
    let strain = core::array::from_fn(|i| elastic[i] + plastic[i]);
    let response = law.response(&strain, &law.initial_state()).unwrap();
    for i in 0..6 {
        near(response.stress_pa[i], stress[i], 2e-10, 2e-12);
        near(response.state.plastic_strain[i], plastic[i], 3e-16, 2e-12);
    }
    near(response.plastic_increment, gamma, 3e-16, 2e-12);
    near(response.state.alpha, gamma, 3e-16, 2e-12);
    near(contract(&stress, &plastic), q * gamma, 2e-14, 2e-12);
}

#[test]
fn stress_and_history_only_paths_match_full_response_in_original_local_axes() {
    let local = law(HARDENING);
    let rotation = rotation([1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0], 0.73);
    for axes in [IDENTITY, rotation] {
        let oriented = OrientedHill48::new(local.clone(), axes).unwrap();
        let mut state = oriented.initial_state();
        for material_strain in [
            [0.0001, -0.00003, 0.00002, 0.0001, 0.0, 0.0],
            [-0.005, 0.001, 0.002, 0.001, -0.0002, 0.0004],
            [0.004, -0.003, -0.001, -0.001, 0.0004, -0.0008],
        ] {
            let strain = rotate(&material_strain, &axes);
            let before = state.clone();
            let full = oriented.response(&strain, &state).unwrap();
            assert_eq!(oriented.stress(&strain, &state), full.stress_pa);
            assert_eq!(oriented.tangent(&strain, &state), full.tangent_pa);
            assert_eq!(oriented.update_state(&strain, &state), full.state);
            assert_eq!(state, before);
            let native = local.response(&material_strain, &state).unwrap();
            assert_eq!(local.stress(&material_strain, &state), native.stress_pa);
            assert_eq!(local.update_state(&material_strain, &state), native.state);
            for i in 0..6 {
                near(
                    full.state.plastic_strain[i],
                    native.state.plastic_strain[i],
                    3e-16,
                    2e-12,
                );
            }
            near(full.state.alpha, native.state.alpha, 3e-16, 2e-12);
            // Inspecting the yielded endpoint still uses the existing roundoff
            // envelope and must not invent another plastic increment.
            assert_eq!(oriented.update_state(&strain, &full.state), full.state);
            state = full.state;
        }
        for bad in [[f64::NAN; 6], [0.1; 6]] {
            let before = state.clone();
            assert!(oriented.response(&bad, &state).is_err());
            assert!(
                oriented
                    .stress(&bad, &state)
                    .iter()
                    .all(|value| value.is_nan())
            );
            assert_eq!(oriented.update_state(&bad, &state), before);
        }
    }
}

#[test]
fn exhausted_accumulation_refuses_on_history_only_paths_without_losing_work() {
    let local = law(0.0);
    let oriented = OrientedHill48::new(local.clone(), IDENTITY).unwrap();
    let history = HillState {
        plastic_strain: [0.0; 6],
        alpha: 1e20,
    };
    let original = history.clone();
    let strain = [0.0, 0.0, 0.0, 0.01, 0.0, 0.0];
    assert!(local.response(&strain, &history).is_err());
    assert!(
        local
            .stress(&strain, &history)
            .iter()
            .all(|value| value.is_nan())
    );
    assert_eq!(local.update_state(&strain, &history), original);
    assert!(oriented.response(&strain, &history).is_err());
    assert!(
        oriented
            .stress(&strain, &history)
            .iter()
            .all(|value| value.is_nan())
    );
    assert_eq!(oriented.update_state(&strain, &history), original);
    assert_eq!(history, original);
    // That same history remains valid for an elastic hold.
    assert_eq!(oriented.update_state(&[0.0; 6], &history), history);
}
