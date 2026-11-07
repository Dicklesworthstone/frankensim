//! Temperature-dependent table flow on the original finite thermal owner.
//! Synthetic curves check equations and retained histories, not alloy calibration.

use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters, FiniteJ2ThermalTabulatedHardening, FiniteNortonThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const KNOTS: [(f64, f64); 5] = [
    (0.0, 2e6),
    (0.04, 5e6),
    (0.12, 6e6),
    (0.3, 6e6),
    (0.8, 20e6),
];
const SLOPES: [f64; 5] = [-4_000.0, -10_000.0, -14_000.0, -14_000.0, -35_000.0];

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.1,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 2.0,
        maximum_accumulated_plastic_strain: 0.8,
    }
}
fn parameters() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
        reference: FiniteJ2Parameters {
            shear_modulus_pa: 80e6,
            bulk_modulus_pa: 160e6,
            yield_stress_pa: 2e6,
            hardening_modulus_pa: 0.0,
        },
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: -20_000.0,
            bulk_modulus_pa_k: -30_000.0,
            yield_stress_pa_k: SLOPES[0],
            hardening_modulus_pa_k: 0.0,
        },
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}
fn fixed_elasticity() -> FiniteJ2ThermalParameters {
    let mut p = parameters();
    p.slopes.shear_modulus_pa_k = 0.0;
    p.slopes.bulk_modulus_pa_k = 0.0;
    p.logarithmic_expansion_per_k = 0.0;
    p.logarithmic_expansion_slope_per_k2 = 0.0;
    p
}
fn curve(ordinates: [f64; 5]) -> TabulatedHardening {
    TabulatedHardening::new(
        KNOTS
            .iter()
            .zip(ordinates)
            .map(|(&(plastic_strain, _), yield_stress_pa)| HardeningPoint {
                plastic_strain,
                yield_stress_pa,
            })
            .collect(),
    )
    .unwrap()
}
fn reference() -> TabulatedHardening {
    curve(KNOTS.map(|p| p.1))
}
fn hardening() -> FiniteJ2ThermalTabulatedHardening {
    FiniteJ2ThermalTabulatedHardening::new(reference(), SLOPES.to_vec()).unwrap()
}
fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(parameters(), limits())
        .unwrap()
        .with_tabulated_hardening(hardening())
        .unwrap()
}
fn coaxial(e: f64) -> [f64; 9] {
    [
        e.exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * e).exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * e).exp(),
    ]
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}
fn ordinates(temperature: f64) -> [f64; 5] {
    core::array::from_fn(|i| SLOPES[i].mul_add(temperature - 300.0, KNOTS[i].1))
}
fn value(alpha: f64, y: [f64; 5]) -> f64 {
    for index in 0..4 {
        if alpha <= KNOTS[index + 1].0 {
            return y[index]
                + (y[index + 1] - y[index])
                    * ((alpha - KNOTS[index].0) / (KNOTS[index + 1].0 - KNOTS[index].0));
        }
    }
    panic!("independent oracle outside support");
}
fn primitive(alpha: f64, y: [f64; 5]) -> f64 {
    let mut result = 0.0;
    for index in 0..4 {
        let (lower, upper) = (KNOTS[index].0, KNOTS[index + 1].0);
        if alpha <= lower {
            break;
        }
        let dx = alpha.min(upper) - lower;
        result +=
            dx * ((y[index] - y[0]) + 0.5 * dx * ((y[index + 1] - y[index]) / (upper - lower)));
    }
    result
}
fn root(mut lower: f64, mut upper: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lower) <= 0.0 && f(upper) >= 0.0);
    for _ in 0..120 {
        let center = 0.5 * lower + 0.5 * upper;
        if f(center) < 0.0 {
            lower = center;
        } else {
            upper = center;
        }
    }
    0.5 * lower + 0.5 * upper
}
fn kinetics() -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 40e6,
        stress_scale_pa: 3e6,
        exponent: 2.0,
        activation_temperature_k: 600.0,
    }
}

#[test]
fn constant_table_recovers_the_original_finite_owner_and_its_physical_heat() {
    let mut p = fixed_elasticity();
    p.slopes.yield_stress_pa_k = 0.0;
    let hardening = FiniteJ2ThermalTabulatedHardening::new(reference(), vec![0.0; 5]).unwrap();
    let family = AffineThermalFiniteJ2::new(p, limits())
        .unwrap()
        .with_tabulated_hardening(hardening.clone())
        .unwrap();
    assert_eq!(family.tabulated_hardening(), Some(&hardening));
    assert_eq!(hardening.reference(), &reference());
    assert_eq!(hardening.yield_stress_slopes_pa_k(), &[0.0; 5]);
    let native = FiniteJ2::new(p.reference, limits())
        .unwrap()
        .with_tabulated_hardening(reference())
        .unwrap();
    let f = [1.2, 0.12, 0.01, 0.0, 0.94, 0.04, 0.02, 0.0, 0.95];
    let old = family.initial_state();
    let expected = native.response(&f, &old).unwrap();
    let step = family
        .checked_increment(&ID, &old, 300.0, &f, 410.0, 1e-3)
        .unwrap();
    assert_eq!(step.mechanical, expected);
    assert_eq!(step.coefficient_storage_change_j_m3, 0.0);
    assert_eq!(step.thermoelastic_storage_change_j_m3, 0.0);
    assert_eq!(step.reversible_heat_j_m3, 0.0);
    assert_eq!(step.caloric_heat_j_m3(), expected.plastic_dissipation_j_m3);
    near(step.balance_error_j_m3, 0.0, 3e-7, 0.0);
}

#[test]
fn exact_signed_table_primitive_supplies_entropy_without_inventing_capacity() {
    let base = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    let family = family();
    let history = family
        .at_temperature(300.0)
        .unwrap()
        .response(&coaxial(0.17), &family.initial_state())
        .unwrap()
        .state;
    let f = [1.12, 0.1, 0.02, -0.03, 0.95, 0.05, 0.01, 0.0, 1.01];
    let temperature = 370.0;
    let baseline = base.thermodynamics(&f, &history, temperature).unwrap();
    let actual = family.thermodynamics(&f, &history, temperature).unwrap();
    let stored = primitive(history.alpha, ordinates(temperature));
    let derivative = primitive(history.alpha, SLOPES);
    near(
        actual.free_energy_j_m3 - baseline.free_energy_j_m3,
        stored,
        2e-8,
        2e-13,
    );
    near(
        actual.entropy_j_m3_k - baseline.entropy_j_m3_k,
        -derivative,
        1e-10,
        2e-13,
    );
    near(
        actual.internal_energy_j_m3 - baseline.internal_energy_j_m3,
        stored - temperature * derivative,
        3e-8,
        2e-13,
    );
    assert_eq!(actual.heat_capacity_j_m3_k, baseline.heat_capacity_j_m3_k);
    assert_eq!(actual.first_piola_pa, baseline.first_piola_pa);
    assert_eq!(
        actual.first_piola_temperature_derivative_fixed_history_pa_k,
        baseline.first_piola_temperature_derivative_fixed_history_pa_k
    );
    assert_eq!(
        family.minimum_mechanical_heat_capacity_j_m3_k(),
        base.minimum_mechanical_heat_capacity_j_m3_k()
    );
    let dt = 0.05;
    let plus = family
        .thermodynamics(&f, &history, temperature + dt)
        .unwrap();
    let minus = family
        .thermodynamics(&f, &history, temperature - dt)
        .unwrap();
    near(
        -(plus.free_energy_j_m3 - minus.free_energy_j_m3) / (2.0 * dt),
        actual.entropy_j_m3_k,
        0.02,
        3e-7,
    );
    near(
        (plus.internal_energy_j_m3 - minus.internal_energy_j_m3) / (2.0 * dt),
        actual.heat_capacity_j_m3_k,
        0.05,
        3e-6,
    );
}

#[test]
fn a_thermal_yield_plateau_keeps_entropy_and_a_closed_unloaded_cycle_keeps_history() {
    let mut p = fixed_elasticity();
    p.slopes.yield_stress_pa_k = 0.0;
    let slopes = KNOTS.map(|(_, y)| -(y - 2e6) / 200.0);
    let hardening = FiniteJ2ThermalTabulatedHardening::new(reference(), slopes.to_vec()).unwrap();
    let family = AffineThermalFiniteJ2::new(p, limits())
        .unwrap()
        .with_tabulated_hardening(hardening)
        .unwrap();
    let old = family
        .at_temperature(300.0)
        .unwrap()
        .response(&coaxial(0.18), &family.initial_state())
        .unwrap()
        .state;
    let f = coaxial(old.alpha);
    let endpoint = family.thermodynamics(&f, &old, 500.0).unwrap();
    near(
        endpoint.entropy_j_m3_k,
        -primitive(old.alpha, slopes),
        1e-10,
        2e-13,
    );
    assert!(endpoint.entropy_j_m3_k > 0.0);
    near(endpoint.free_energy_j_m3, 0.0, 2e-8, 0.0);
    near(
        endpoint.internal_energy_j_m3,
        500.0 * endpoint.entropy_j_m3_k,
        2e-8,
        2e-13,
    );
    for (old_temperature, new_temperature) in [(300.0, 500.0), (500.0, 300.0)] {
        let step = family
            .checked_increment(&f, &old, old_temperature, &f, new_temperature, 1e-3)
            .unwrap();
        assert_eq!(step.mechanical.state, old);
        assert_eq!(step.mechanical.plastic_increment, 0.0);
        let change = (new_temperature - old_temperature) * primitive(old.alpha, slopes);
        near(step.coefficient_storage_change_j_m3, change, 2e-8, 2e-13);
        near(step.thermoelastic_storage_change_j_m3, -change, 2e-8, 2e-13);
        near(step.internal_energy_change_j_m3(), 0.0, 5e-8, 0.0);
        near(step.caloric_heat_j_m3(), 0.0, 5e-8, 0.0);
        near(step.balance_error_j_m3, 0.0, 5e-8, 0.0);
    }
}

#[test]
fn heating_at_held_deformation_advances_yielding_and_cooling_does_not_reset_plasticity() {
    let p = fixed_elasticity();
    let family = AffineThermalFiniteJ2::new(p, limits())
        .unwrap()
        .with_tabulated_hardening(hardening())
        .unwrap();
    let f = coaxial(0.14);
    let old = family
        .at_temperature(300.0)
        .unwrap()
        .response(&f, &family.initial_state())
        .unwrap()
        .state;
    let retained = old.clone();
    assert!(family.validate_committed(&f, &old, 450.0, 1e-3).is_err());
    let expected = root(old.alpha, 0.14, |alpha| {
        3.0 * p.reference.shear_modulus_pa * alpha + value(alpha, ordinates(450.0))
            - 3.0 * p.reference.shear_modulus_pa * 0.14
    });
    let heated = family
        .checked_increment(&f, &old, 300.0, &f, 450.0, 1e-3)
        .unwrap();
    near(heated.mechanical.state.alpha, expected, 3e-14, 5e-12);
    assert!(heated.mechanical.state.alpha > old.alpha);
    near(
        heated.mechanical.equivalent_kirchhoff_pa,
        value(expected, ordinates(450.0)),
        4e-6,
        5e-12,
    );
    near(
        heated.mechanical.plastic_dissipation_j_m3,
        ordinates(450.0)[0] * (expected - old.alpha),
        3e-7,
        5e-12,
    );
    near(heated.mechanical_work_j_m3, 0.0, 1e-12, 0.0);
    near(heated.balance_error_j_m3, 0.0, 3e-7, 0.0);
    assert_eq!(old, retained);
    let cooled = family
        .checked_increment(&f, &heated.mechanical.state, 450.0, &f, 300.0, 1e-3)
        .unwrap();
    assert_eq!(cooled.mechanical.state, heated.mechanical.state);
    assert_eq!(cooled.mechanical.plastic_dissipation_j_m3, 0.0);
    near(cooled.caloric_heat_j_m3(), 0.0, 3e-7, 0.0);
}

#[test]
fn noncoaxial_hot_return_keeps_all_tensor_derivatives_and_objectivity() {
    let family = family();
    let old = family
        .at_temperature(300.0)
        .unwrap()
        .response(&coaxial(0.07), &family.initial_state())
        .unwrap()
        .state;
    let law = family.at_temperature(410.0).unwrap();
    let f = [1.3, 0.15, 0.03, 0.02, 0.9, 0.08, 0.01, -0.03, 0.93];
    let r = law.response(&f, &old).unwrap();
    assert!(r.plastic_increment > 0.0);
    for column in 0..9 {
        let h = 1e-6;
        let mut plus = f;
        let mut minus = f;
        plus[column] += h;
        minus[column] -= h;
        let (pp, wp) = law.stress_and_potential(&plus, &old).unwrap();
        let (pm, wm) = law.stress_and_potential(&minus, &old).unwrap();
        near((wp - wm) / (2.0 * h), r.first_piola_pa[column], 0.3, 5e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * h),
                r.tangent_pa[row][column],
                3.0,
                5e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                3e-12,
            );
        }
    }
    let multiply = |a: &[f64; 9], b: &[f64; 9]| -> [f64; 9] {
        core::array::from_fn(|ij| {
            (0..3)
                .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
                .sum()
        })
    };
    let (s, c) = 0.67_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let turned = law.response(&multiply(&rotation, &f), &old).unwrap();
    for (a, b) in turned
        .first_piola_pa
        .into_iter()
        .zip(multiply(&rotation, &r.first_piola_pa))
    {
        near(a, b, 5e-6, 3e-12);
    }
    near(turned.state.alpha, r.state.alpha, 3e-14, 3e-12);
    near(turned.stored_energy_j_m3, r.stored_energy_j_m3, 3e-7, 3e-12);
}

#[test]
fn hot_tables_compose_with_native_perzyna_secondary_and_primary_creep() {
    let p = fixed_elasticity();
    let base = AffineThermalFiniteJ2::new(p, limits()).unwrap();
    let temperature = 410.0;
    let duration = 0.12;
    let f = coaxial(0.2);
    let primary = StrainHardeningCreep {
        reference_strain: 0.015,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    };
    for mechanism in 0..3 {
        let selected = base.clone().with_tabulated_hardening(hardening()).unwrap();
        let (family, opposite) = if mechanism == 0 {
            (
                selected.with_overstress(kinetics()).unwrap(),
                base.clone()
                    .with_overstress(kinetics())
                    .unwrap()
                    .with_tabulated_hardening(hardening())
                    .unwrap(),
            )
        } else if mechanism == 1 {
            (
                selected.with_norton_creep(kinetics(), 2.0).unwrap(),
                base.clone()
                    .with_norton_creep(kinetics(), 2.0)
                    .unwrap()
                    .with_tabulated_hardening(hardening())
                    .unwrap(),
            )
        } else {
            (
                selected
                    .with_norton_creep(kinetics(), 2.0)
                    .unwrap()
                    .with_norton_strain_hardening(primary)
                    .unwrap(),
                base.clone()
                    .with_norton_creep(kinetics(), 2.0)
                    .unwrap()
                    .with_norton_strain_hardening(primary)
                    .unwrap()
                    .with_tabulated_hardening(hardening())
                    .unwrap(),
            )
        };
        assert_eq!(family, opposite);
        let old = family.initial_state();
        assert!(family.at_temperature(temperature).is_err());
        assert!(
            family
                .at_temperature_for_inspection(temperature)
                .unwrap()
                .response(&f, &old)
                .is_err()
        );
        let native = FiniteJ2::new(
            FiniteJ2Parameters {
                yield_stress_pa: ordinates(temperature)[0],
                ..p.reference
            },
            limits(),
        )
        .unwrap()
        .with_tabulated_hardening(curve(ordinates(temperature)))
        .unwrap();
        let kinetic = kinetics();
        let rate = FiniteOverstressParameters {
            viscosity_pa_s: kinetic.reference_viscosity_pa_s
                * (kinetic.activation_temperature_k * (1.0 / temperature - 1.0 / 300.0)).exp(),
            stress_scale_pa: kinetic.stress_scale_pa,
            exponent: kinetic.exponent,
        };
        let expected = if mechanism == 0 {
            FiniteOverstressJ2::new(native, rate, duration)
                .unwrap()
                .response(&f, &old)
                .unwrap()
        } else {
            let mut combined = FiniteNortonJ2::new(native, rate, duration, 2.0).unwrap();
            if mechanism == 2 {
                combined = combined.with_strain_hardening(primary).unwrap();
            }
            combined.response(&f, &old).unwrap()
        };
        let step = family
            .checked_increment_for_step(&ID, &old, 300.0, &f, temperature, 1e-3, duration)
            .unwrap();
        let actual = &step.mechanical;
        assert!(actual.plastic_increment > 0.0);
        if mechanism != 0 {
            assert!(actual.state.creep_alpha > 0.0);
        }
        near(actual.state.alpha, expected.state.alpha, 3e-14, 5e-12);
        near(
            actual.state.creep_alpha,
            expected.state.creep_alpha,
            3e-14,
            5e-12,
        );
        near(
            actual.inelastic_dissipation_j_m3(),
            expected.inelastic_dissipation_j_m3(),
            3e-7,
            5e-12,
        );
        near(
            actual.flow_potential_j_m3,
            expected.flow_potential_j_m3,
            3e-7,
            5e-12,
        );
        for row in 0..9 {
            near(
                actual.first_piola_pa[row],
                expected.first_piola_pa[row],
                4e-6,
                5e-12,
            );
            for column in 0..9 {
                near(
                    actual.tangent_pa[row][column],
                    expected.tangent_pa[row][column],
                    3e-5,
                    5e-12,
                );
            }
        }
        near(
            step.reversible_heat_j_m3,
            temperature * primitive(actual.state.alpha, SLOPES),
            3e-7,
            5e-12,
        );
        near(step.balance_error_j_m3, 0.0, 4e-7, 0.0);
    }
}

#[test]
fn signed_storage_derivatives_cross_zero_without_losing_the_original_curve() {
    let slopes = [0.0, 2000.0, -2000.0, 1000.0, -1000.0];
    let hardening = FiniteJ2ThermalTabulatedHardening::new(reference(), slopes.to_vec()).unwrap();
    for alpha in [0.0, 0.01, 0.04, 0.08, 0.12, 0.24, 0.3, 0.5, 0.8] {
        near(
            hardening
                .hardening_energy_temperature_derivative_j_m3_k(alpha)
                .unwrap(),
            primitive(alpha, slopes),
            1e-12,
            5e-14,
        );
    }
    assert!(
        hardening
            .hardening_energy_temperature_derivative_j_m3_k(-1e-6)
            .is_err()
    );
    assert!(
        hardening
            .hardening_energy_temperature_derivative_j_m3_k(0.800001)
            .is_err()
    );
    assert!(
        hardening
            .hardening_energy_temperature_derivative_j_m3_k(f64::NAN)
            .is_err()
    );
    assert_eq!(hardening.reference(), &reference());
}

#[test]
fn admission_covers_all_temperature_endpoints_and_refuses_conflicting_histories() {
    let base = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    for slopes in [vec![0.0; 4], vec![f64::NAN; 5]] {
        assert!(FiniteJ2ThermalTabulatedHardening::new(reference(), slopes).is_err());
    }
    for slopes in [
        [-4_000.0, -30_000.0, -14_000.0, -14_000.0, -35_000.0],
        [-4_000.0, -10_000.0, -18_000.0, -18_000.0, -35_000.0],
        [-4_001.0, -10_000.0, -14_000.0, -14_000.0, -35_000.0],
    ] {
        let table = FiniteJ2ThermalTabulatedHardening::new(reference(), slopes.to_vec()).unwrap();
        assert!(base.clone().with_tabulated_hardening(table).is_err());
    }
    let selected = base.clone().with_tabulated_hardening(hardening()).unwrap();
    assert_eq!(
        selected
            .clone()
            .with_tabulated_hardening(hardening())
            .unwrap(),
        selected
    );
    assert!(selected.clone().with_voce_hardening(0.0, 1.0, 0.0).is_err());
    assert!(
        base.clone()
            .with_voce_hardening(0.0, 1.0, 0.0)
            .unwrap()
            .with_tabulated_hardening(hardening())
            .is_err()
    );
    let replacement = FiniteJ2ThermalTabulatedHardening::new(
        curve([2e6, 5.1e6, 6e6, 6e6, 20e6]),
        SLOPES.to_vec(),
    )
    .unwrap();
    assert!(
        selected
            .clone()
            .with_tabulated_hardening(replacement)
            .is_err()
    );
    let mut p = parameters();
    p.reference.hardening_modulus_pa = 1e6;
    assert!(
        AffineThermalFiniteJ2::new(p, limits())
            .unwrap()
            .with_tabulated_hardening(hardening())
            .is_err()
    );
    p = parameters();
    p.reference_temperature_k = 280.0;
    p.slopes.hardening_modulus_pa_k = 1.0;
    assert!(
        AffineThermalFiniteJ2::new(p, limits())
            .unwrap()
            .with_tabulated_hardening(hardening())
            .is_err()
    );
    assert!(
        AffineThermalFiniteJ2::new(
            parameters(),
            FiniteJ2Limits {
                maximum_accumulated_plastic_strain: 0.800001,
                ..limits()
            }
        )
        .unwrap()
        .with_tabulated_hardening(hardening())
        .is_err()
    );
    let p = parameters();
    let creep = AffineThermalFiniteJ2::new_norton(
        FiniteNortonThermalParameters {
            reference_temperature_k: p.reference_temperature_k,
            temperature_range_k: p.temperature_range_k,
            shear_modulus_pa: p.reference.shear_modulus_pa,
            bulk_modulus_pa: p.reference.bulk_modulus_pa,
            shear_modulus_slope_pa_k: p.slopes.shear_modulus_pa_k,
            bulk_modulus_slope_pa_k: p.slopes.bulk_modulus_pa_k,
            logarithmic_expansion_per_k: p.logarithmic_expansion_per_k,
            logarithmic_expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
        },
        limits(),
        kinetics(),
    )
    .unwrap();
    assert!(creep.with_tabulated_hardening(hardening()).is_err());
    assert!(selected.at_temperature(500.01).is_err());
    let old = selected.initial_state();
    assert!(
        selected
            .at_temperature(500.0)
            .unwrap()
            .response(&coaxial(1.1), &old)
            .is_err()
    );
    assert_eq!(old, selected.initial_state());
}
