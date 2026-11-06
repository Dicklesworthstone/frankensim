//! G0/G2/G3: independent directional coupon solutions, J2 recovery, tensor
//! derivatives, hydrostatic invariance and work accounting. Coefficients are
//! synthetic; these tests do not establish a real rolled alloy calibration.

use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients, HillError, HillState};
use fs_material::tensor::{contract, trace};
use fs_material::{SmallStrainLaw, Voigt};

const E: f64 = 210e9;
const NU: f64 = 0.3;
const Y: f64 = 250e6;
const HARDENING: f64 = 1.2e9;

fn coefficients() -> Hill48Coefficients {
    Hill48Coefficients {
        f: 0.25,
        g: 0.9,
        h: 0.45,
        l: 1.4,
        m: 2.1,
        n: 1.8,
    }
}

fn material(c: Hill48Coefficients, hardening: f64) -> Hill48 {
    Hill48::new(IsotropicElastic::new(E, NU, 0.05).unwrap(), Y, hardening, c).unwrap()
}

fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "actual={a:.16e}, expected={b:.16e}, difference={:.8e}",
        (a - b).abs()
    );
}

fn vector(a: &Voigt, b: &Voigt, absolute: f64, relative: f64) {
    for i in 0..6 {
        near(a[i], b[i], absolute, relative);
    }
}

/// Independent definition, without using the law's modes or equivalent stress.
fn hill_equivalent(c: Hill48Coefficients, s: &Voigt) -> f64 {
    (c.f * (s[1] - s[2]).powi(2)
        + c.g * (s[2] - s[0]).powi(2)
        + c.h * (s[0] - s[1]).powi(2)
        + 2.0 * c.n * s[3].powi(2)
        + 2.0 * c.l * s[4].powi(2)
        + 2.0 * c.m * s[5].powi(2))
    .sqrt()
}

#[test]
fn all_unit_yield_ratios_recover_j2_on_nonproportional_load_paths() {
    let c = Hill48Coefficients::from_yield_ratios([1.0; 6]).unwrap();
    assert_eq!(c, Hill48Coefficients::isotropic());
    for hardening in [0.0, HARDENING] {
        let hill = material(c, hardening);
        let j2 = J2Plasticity::new(hill.elastic().clone(), Y, hardening).unwrap();
        let mut hs = hill.initial_state();
        let mut js = j2.initial_state();
        for strain in [
            [0.0001, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.005, -0.001, 0.0004, 0.001, 0.0002, -0.0008],
            [0.002, 0.003, -0.001, -0.002, 0.0005, 0.001],
            [-0.004, 0.001, 0.001, 0.0, -0.001, 0.0015],
            [0.0; 6],
        ] {
            let r = hill.response(&strain, &hs).unwrap();
            vector(&r.stress_pa, &j2.stress(&strain, &js), 1e-4, 2e-12);
            let tangent = j2.tangent(&strain, &js);
            for (a, b) in r.tangent_pa.iter().zip(tangent.iter()) {
                vector(a, b, 0.03, 3e-12);
            }
            js = j2.update_state(&strain, &js);
            vector(&r.state.plastic_strain, &js.plastic_strain, 2e-16, 2e-12);
            near(r.state.alpha, js.alpha, 2e-16, 2e-12);
            hs = r.state;
        }
    }
}

#[test]
fn directional_tension_compression_and_shear_match_closed_form_coupons() {
    let c = coefficients();
    let law = material(c, HARDENING);
    let alpha = 0.003;
    let radius = Y + HARDENING * alpha;
    // Pure axial stress has q=k|s| and the associated normal plastic strain
    // direction is the gradient of the ORIGINAL Hill quadratic form.
    for axis in 0..3 {
        let coefficients = match axis {
            0 => [c.g + c.h, -c.h, -c.g],
            1 => [-c.h, c.f + c.h, -c.f],
            _ => [-c.g, -c.f, c.f + c.g],
        };
        let k = coefficients[axis].sqrt();
        for sign in [-1.0, 1.0] {
            let sigma = sign * radius / k;
            let expected_plastic: Voigt = core::array::from_fn(|i| {
                if i < 3 {
                    sign * alpha * coefficients[i] / k
                } else {
                    0.0
                }
            });
            let strain = core::array::from_fn(|i| {
                expected_plastic[i]
                    + if i == axis {
                        sigma / E
                    } else if i < 3 {
                        -NU * sigma / E
                    } else {
                        0.0
                    }
            });
            let r = law.response(&strain, &law.initial_state()).unwrap();
            vector(
                &r.stress_pa,
                &core::array::from_fn(|i| if i == axis { sigma } else { 0.0 }),
                3e-5,
                2e-13,
            );
            vector(&r.state.plastic_strain, &expected_plastic, 3e-16, 2e-13);
            near(r.state.alpha, alpha, 3e-16, 2e-13);
            near(trace(&r.state.plastic_strain), 0.0, 3e-18, 0.0);
            near(hill_equivalent(c, &r.stress_pa), radius, 3e-5, 2e-13);
        }
    }
    let twice_mu = E / (1.0 + NU);
    for (component, coefficient) in [(3, c.n), (4, c.l), (5, c.m)] {
        let sigma = radius / (2.0 * coefficient).sqrt();
        let plastic_shear = alpha * (0.5 * coefficient).sqrt();
        let mut strain = [0.0; 6];
        strain[component] = sigma / twice_mu + plastic_shear;
        let r = law.response(&strain, &law.initial_state()).unwrap();
        near(r.stress_pa[component], sigma, 3e-5, 2e-13);
        near(
            r.state.plastic_strain[component],
            plastic_shear,
            3e-16,
            2e-13,
        );
        near(r.state.alpha, alpha, 3e-16, 2e-13);
        near(hill_equivalent(c, &r.stress_pa), radius, 3e-5, 2e-13);
    }
}

#[test]
fn directional_ratios_have_explicit_shear_normalization() {
    let ratios = [1.0, 1.15, 0.9, 1.25, 0.8, 1.1];
    let law = material(Hill48Coefficients::from_yield_ratios(ratios).unwrap(), 0.0);
    for component in 0..6 {
        let mut stress = [0.0; 6];
        stress[component] =
            Y * ratios[component] / if component < 3 { 1.0 } else { 3.0_f64.sqrt() };
        near(law.equivalent_stress_pa(&stress).unwrap(), Y, 2e-5, 2e-14);
    }
}

#[test]
fn sheet_r_fit_predicts_directional_yields_and_keeps_transverse_shear_independent() {
    let [r0, r45, r90] = [1.8, 1.1, 2.4];
    let shear_ratios = [0.75, 1.4];
    let law = material(
        Hill48Coefficients::from_planar_r_values([r0, r45, r90], shear_ratios).unwrap(),
        0.0,
    );
    // Closed-form predictions written directly in terms of coupon r-values,
    // not reconstructed from the production Hill coefficients.
    let transverse = Y * (r90 * (1.0 + r0) / (r0 * (1.0 + r90))).sqrt();
    let thickness = Y * (r90 * (1.0 + r0) / (r0 + r90)).sqrt();
    let in_plane_shear = Y * ((1.0 + r0) * r90 / ((2.0 * r45 + 1.0) * (r0 + r90))).sqrt();
    let expected = [
        Y,
        transverse,
        thickness,
        in_plane_shear,
        Y * shear_ratios[0] / 3.0_f64.sqrt(),
        Y * shear_ratios[1] / 3.0_f64.sqrt(),
    ];
    for axis in 0..6 {
        let mut stress = [0.0; 6];
        stress[axis] = expected[axis];
        near(law.equivalent_stress_pa(&stress).unwrap(), Y, 3e-5, 2e-13);
    }
    let diagonal = Y * (2.0 * r90 * (1.0 + r0) / ((1.0 + r45) * (r0 + r90))).sqrt();
    let stress = [
        0.5 * diagonal,
        0.5 * diagonal,
        0.0,
        0.5 * diagonal,
        0.0,
        0.0,
    ];
    near(law.equivalent_stress_pa(&stress).unwrap(), Y, 3e-5, 2e-13);
    assert!((transverse / Y - 1.0).abs() > 0.01);
    assert!((diagonal / Y - 1.0).abs() > 0.1);
}

#[test]
fn anisotropic_tangent_is_the_derivative_at_fixed_accepted_history() {
    let law = material(coefficients(), HARDENING);
    let old = law
        .response(
            &[0.004, 0.001, -0.0005, 0.001, -0.0002, 0.0],
            &law.initial_state(),
        )
        .unwrap()
        .state;
    for (strain, state) in [
        (
            [0.0001, -0.00003, 0.00002, 0.00004, 0.0, 0.00001],
            law.initial_state(),
        ),
        ([-0.003, 0.006, 0.001, -0.001, 0.002, -0.0004], old),
    ] {
        let expected = law.response(&strain, &state).unwrap();
        for column in 0..6 {
            let mut plus = strain;
            let mut minus = strain;
            plus[column] += 1e-8;
            minus[column] -= 1e-8;
            let a = law.response(&plus, &state).unwrap();
            let b = law.response(&minus, &state).unwrap();
            for row in 0..6 {
                near(
                    (a.stress_pa[row] - b.stress_pa[row]) / 2e-8,
                    expected.tangent_pa[row][column],
                    30.0,
                    3e-7,
                );
                let wr = if row < 3 { 1.0 } else { 2.0 };
                let wc = if column < 3 { 1.0 } else { 2.0 };
                near(
                    wr * expected.tangent_pa[row][column],
                    wc * expected.tangent_pa[column][row],
                    1e-4,
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn accepted_endpoints_do_not_generate_plastic_steps_from_roundoff() {
    for hardening in [0.0, HARDENING] {
        let law = material(coefficients(), hardening);
        let mut state = law.initial_state();
        for strain in [
            [0.005, -0.001, 0.0004, 0.001, 0.0002, -0.0008],
            [0.002, 0.003, -0.001, -0.002, 0.0005, 0.001],
            [-0.004, 0.001, 0.001, 0.0, -0.001, 0.0015],
        ] {
            let accepted = law.response(&strain, &state).unwrap();
            for _ in 0..3 {
                let repeated = law.response(&strain, &accepted.state).unwrap();
                assert_eq!(repeated.state, accepted.state);
                assert_eq!(repeated.plastic_increment, 0.0);
                vector(&repeated.stress_pa, &accepted.stress_pa, 1e-4, 2e-13);
            }
            state = accepted.state;
        }
    }
}

#[test]
fn hydrostatic_pressure_preserves_yielding_and_flow_is_isochoric() {
    let law = material(coefficients(), HARDENING);
    let pure = law
        .response(
            &[-0.008, -0.008, -0.008, 0.0, 0.0, 0.0],
            &law.initial_state(),
        )
        .unwrap();
    assert_eq!(pure.state, law.initial_state());
    assert_eq!(pure.plastic_increment, 0.0);
    let strain = [0.003, -0.001, 0.0005, 0.001, 0.0005, 0.0];
    let shifted = core::array::from_fn(|i| strain[i] - if i < 3 { 0.002 } else { 0.0 });
    let a = law.response(&strain, &law.initial_state()).unwrap();
    let b = law.response(&shifted, &law.initial_state()).unwrap();
    vector(
        &a.state.plastic_strain,
        &b.state.plastic_strain,
        3e-16,
        2e-13,
    );
    near(a.state.alpha, b.state.alpha, 3e-16, 2e-13);
    for i in 0..6 {
        near(
            a.stress_pa[i] - b.stress_pa[i],
            if i < 3 {
                E * 0.002 / (1.0 - 2.0 * NU)
            } else {
                0.0
            },
            3e-5,
            2e-13,
        );
    }
    near(trace(&a.state.plastic_strain), 0.0, 3e-18, 0.0);
}

#[test]
fn equivalent_stress_preserves_small_deviatoric_differences_under_large_pressure() {
    let law = material(Hill48Coefficients::isotropic(), HARDENING);
    for magnitude in [1e-200_f64, 1e-9, 1.0, 1e9, 1e200, 1e300] {
        for sign in [-1.0, 1.0] {
            let x = sign * magnitude;
            let previous = x.next_down();
            // For two equal principal stresses, Mises q is exactly the
            // difference from the third. These neighboring inputs subtract
            // exactly; dividing each by their common scale first does not.
            let expected = (x - previous).abs();
            let actual = law
                .equivalent_stress_pa(&[x, previous, x, 0.0, 0.0, 0.0])
                .unwrap();
            near(actual, expected, 0.0, 8.0 * f64::EPSILON);
            assert_eq!(
                law.equivalent_stress_pa(&[x, x, x, 0.0, 0.0, 0.0]).unwrap(),
                0.0,
            );
        }
    }
}

#[test]
fn compression_retains_permanent_set_and_accounts_for_work() {
    let law = material(coefficients(), HARDENING);
    let initial = law.initial_state();
    let tolerances = EnergyTolerances {
        energy_j_m3: 1e-6,
        yield_pa: 1e-4,
    };
    let strain = [-0.006, 0.001, 0.002, 0.001, -0.0003, 0.0006];
    let first = law
        .advance_with_energy(&[0.0; 6], &initial, &strain, tolerances)
        .unwrap();
    assert!(first.response.state.alpha > 0.0);
    near(
        first.energy.physical_dissipation_j_m3,
        Y * first.response.state.alpha,
        1e-8,
        2e-14,
    );
    near(
        first.energy.final_state.hardening_j_m3,
        0.5 * HARDENING * first.response.state.alpha.powi(2),
        1e-8,
        2e-14,
    );
    assert!(first.energy.numerical_dissipation_j_m3 > 0.0);
    assert!(first.energy.residual_j_m3.abs() < 1e-6);
    let set = first.response.state.plastic_strain;
    let unload = law
        .advance_with_energy(&strain, &first.response.state, &set, tolerances)
        .unwrap();
    vector(&unload.response.stress_pa, &[0.0; 6], 1e-5, 0.0);
    assert_eq!(unload.response.state, first.response.state);
    assert_eq!(unload.energy.physical_dissipation_j_m3, 0.0);
    assert!(contract(&set, &set) > 0.0);
    assert_eq!(initial, HillState::default());
}

#[test]
fn signed_normal_coefficients_are_admitted_by_convexity_not_component_sign() {
    let c = Hill48Coefficients {
        f: -0.2,
        g: 0.8,
        h: 1.2,
        ..coefficients()
    };
    let law = material(c, HARDENING);
    let strain = [0.005, -0.001, 0.001, 0.001, 0.0002, -0.0001];
    let r = law.response(&strain, &law.initial_state()).unwrap();
    near(
        hill_equivalent(c, &r.stress_pa),
        Y + HARDENING * r.state.alpha,
        4e-5,
        2e-13,
    );
    for invalid in [
        Hill48Coefficients { f: -2.0, ..c },
        Hill48Coefficients {
            f: 0.0,
            g: 0.0,
            h: 1.0,
            ..c
        },
        Hill48Coefficients { n: 0.0, ..c },
        Hill48Coefficients { f: f64::NAN, ..c },
    ] {
        assert!(Hill48::new(law.elastic().clone(), Y, HARDENING, invalid).is_err());
    }
}

#[test]
fn stress_unit_rescaling_does_not_change_plastic_flow() {
    let strain = [0.005, -0.001, 0.0005, 0.001, 0.0004, -0.0008];
    let reference = material(coefficients(), HARDENING)
        .response(&strain, &HillState::default())
        .unwrap();
    for scale in [1e-120, 1e120] {
        let law = Hill48::new(
            IsotropicElastic::new(E * scale, NU, 0.05).unwrap(),
            Y * scale,
            HARDENING * scale,
            coefficients(),
        )
        .unwrap();
        let r = law.response(&strain, &law.initial_state()).unwrap();
        vector(
            &r.stress_pa.map(|x| x / scale),
            &reference.stress_pa,
            1e-4,
            3e-12,
        );
        vector(
            &r.state.plastic_strain,
            &reference.state.plastic_strain,
            3e-16,
            3e-12,
        );
        near(r.state.alpha, reference.state.alpha, 3e-16, 3e-12);
        for i in 0..6 {
            vector(
                &r.tangent_pa[i].map(|x| x / scale),
                &reference.tangent_pa[i],
                0.01,
                3e-12,
            );
        }
    }
}

#[test]
fn invalid_history_and_domain_exits_refuse_without_mutating_state() {
    let law = material(coefficients(), HARDENING);
    for state in [
        HillState {
            alpha: -0.001,
            ..HillState::default()
        },
        HillState {
            plastic_strain: [0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
            alpha: 1.0,
        },
        HillState {
            plastic_strain: [0.001, -0.001, 0.0, 0.0, 0.0, 0.0],
            alpha: 0.0,
        },
    ] {
        let saved = state.clone();
        assert!(matches!(
            law.response(&[0.0; 6], &state),
            Err(HillError::State(_))
        ));
        assert_eq!(state, saved);
    }
    for strain in [[f64::NAN; 6], [0.1; 6]] {
        assert!(law.response(&strain, &law.initial_state()).is_err());
    }
    // Native infallible traits expose a failed constitutive candidate to the
    // checked solver through nonfinite stress/tangent, preserving old history.
    assert!(
        law.stress(&[0.1; 6], &law.initial_state())
            .iter()
            .all(|v| v.is_nan())
    );
    assert_eq!(
        law.update_state(&[0.1; 6], &law.initial_state()),
        law.initial_state()
    );
}
