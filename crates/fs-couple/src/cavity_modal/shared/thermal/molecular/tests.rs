//! G1 independent constitutive/time-energy oracles; G3 diffusion refinement.
//! Synthetic reservoir rates test the equations, not measured gas accuracy.
use super::*;
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

fn material(reservoirs: bool, split: bool) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            molar_mass: 0.0289644,
            gamma: 1.25,
            sutherland_beta: 1.458e-6,
            sutherland_s: 110.4,
            conductivity: ConductivityModel::Eucken,
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let ids: Vec<_> = if split { vec!["A", "B"] } else { vec!["A"] }
        .into_iter()
        .map(|s| SpeciesId::new(s).unwrap())
        .collect();
    let species = ids
        .iter()
        .map(|id| SpeciesRelaxation {
            species: id.clone(),
            processes: if reservoirs {
                vec![RelaxationProcess {
                    reservoir: InternalReservoir::Rotational {
                        degrees_of_freedom: 2,
                    },
                    rate: RelaxationRateLaw::PartnerIndependent {
                        p_tau_ref_pa_s: 101325.0 * 0.025,
                        t_ref_k: 293.15,
                        exponent: 0.0,
                        t_min_k: 200.0,
                        t_max_k: 400.0,
                    },
                }]
            } else {
                vec![]
            },
        })
        .collect();
    let components = ids
        .iter()
        .map(|id| (id.clone(), 1.0 / ids.len() as f64))
        .collect();
    RelaxationAbsorption::new(species)
        .unwrap()
        .evaluate(&RelaxingMedium::new(gas, components).unwrap())
        .unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "{actual:e} != {expected:e}"
    );
}

fn complex_close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn one_thermal_pole_and_reservoir_match_independent_quadratic_and_limits() {
    let law = material(true, false);
    let gas = law.gas();
    let modes = [ThermalMode {
        weight: 1.0,
        rate_per_s: 35.0,
    }];
    let result = realize(&law, &modes, [0.001, 1e5], 2, &mut || false).unwrap();
    let a = law.frozen_cv_j_kg_k() + gas.specific_gas_constant;
    let r = 40.0;
    let d = gas.specific_heat_cp * modes[0].rate_per_s;
    for omega in [0.001, 0.1, 3.0, 17.0, 40.0, 90.0, 1000.0, 1e5] {
        let s = C64::new(0.0, omega);
        let numerator =
            (s * s).scale(a) + s.scale(gas.specific_heat_cp * r + d) + C64::from_re(d * r);
        let denominator = (s * s).scale(a - gas.specific_gas_constant)
            + s.scale((gas.specific_heat_cp - gas.specific_gas_constant) * r + d)
            + C64::from_re(d * r);
        let exact = (numerator / denominator).scale(gas.pressure);
        let actual = result.bulk_at(gas.pressure, omega);
        complex_close(actual, exact, 3e-12);
        close(actual.im, exact.im, 3e-12);
        assert!(actual.im > 0.0);
    }
    close(result.bulk_at(gas.pressure, 0.0).re, gas.pressure, 1e-15);
    close(
        result.frozen_bulk_modulus_pa,
        gas.pressure * a / law.frozen_cv_j_kg_k(),
        3e-12,
    );
}

// Independent complex modified-Bessel power series. q=z^2/4; return
// 2 I1(z)/(z I0(z)) without a complex square root or either owner kernel.
fn thermal_bessel(z_squared: C64) -> C64 {
    let q = z_squared.scale(0.25);
    let (mut i0, mut twice_i1_over_z) = (C64::ONE, C64::ONE);
    let (mut term0, mut term1) = (C64::ONE, C64::ONE);
    for k in 1..100 {
        term0 = (term0 * q).scale(1.0 / (k * k) as f64);
        term1 = (term1 * q).scale(1.0 / (k * (k + 1)) as f64);
        i0 = i0 + term0;
        twice_i1_over_z = twice_i1_over_z + term1;
    }
    twice_i1_over_z / i0
}

#[test]
fn joint_cylindrical_law_converges_to_independent_complex_bessel_heat_diffusion() {
    let law = material(true, false);
    let gas = law.gas();
    let radius = 0.001;
    let mut previous = f64::INFINITY;
    for poles in [2, 4, 8] {
        let wall = fs_phs::ZkDiffusive::new(
            radius,
            gas.density,
            gas.dynamic_viscosity,
            gas.gamma,
            gas.prandtl,
            poles,
            [10.0, 5000.0],
            0.95,
        )
        .unwrap();
        let modes: Vec<_> = wall
            .shunt_terms()
            .iter()
            .map(|term| ThermalMode {
                weight: term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0)),
                rate_per_s: term.rate_per_s,
            })
            .collect();
        let result = realize(&law, &modes, [10.0, 5000.0], 2 * poles, &mut || false).unwrap();
        let mut error = 0.0_f64;
        for omega in [10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 5000.0] {
            let cp = law.caloric_response(omega).unwrap().cp_j_kg_k.conj();
            let z2 = (C64::new(0.0, omega) * cp)
                .scale(gas.density * radius * radius / gas.thermal_conductivity);
            let f = thermal_bessel(z2);
            let beta_p = C64::ONE - ((C64::ONE - f) / cp).scale(gas.specific_gas_constant);
            let exact = beta_p.recip().scale(gas.pressure);
            let actual = result.bulk_at(gas.pressure, omega);
            error = error.max((actual - exact).abs() / exact.abs());
            error = error.max((actual.im - exact.im).abs() / exact.im);
        }
        assert!(
            error < previous,
            "poles={poles}, error={error:e}, previous={previous:e}"
        );
        previous = error;
    }
    assert!(previous < 2e-5, "final complex/loss error={previous:e}");
}

#[test]
fn reservoir_splitting_and_nonrelaxing_gas_preserve_the_original_physics() {
    let single = material(true, false);
    let split = material(true, true);
    let none = material(false, false);
    let modes = [
        ThermalMode {
            weight: 0.7,
            rate_per_s: 13.0,
        },
        ThermalMode {
            weight: 0.3,
            rate_per_s: 170.0,
        },
    ];
    let a = realize(&single, &modes, [0.1, 1e4], 4, &mut || false).unwrap();
    let b = realize(&split, &modes, [0.1, 1e4], 6, &mut || false).unwrap();
    let c = realize(&none, &modes, [0.1, 1e4], 2, &mut || false).unwrap();
    for omega in [0.1, 3.0, 30.0, 100.0, 1000.0, 1e4] {
        complex_close(
            a.bulk_at(single.gas().pressure, omega),
            b.bulk_at(split.gas().pressure, omega),
            5e-12,
        );
        let s = C64::new(0.0, omega);
        let heat = modes.iter().fold(C64::ZERO, |sum, mode| {
            sum + (C64::from_re(mode.rate_per_s) / (s + C64::from_re(mode.rate_per_s)))
                .scale(mode.weight)
        });
        let gas = none.gas();
        let exact = (C64::ONE + heat.scale(gas.gamma - 1.0))
            .recip()
            .scale(gas.pressure * gas.gamma);
        complex_close(c.bulk_at(gas.pressure, omega), exact, 3e-12);
    }
}

#[test]
fn vanishing_wall_exchange_recovers_the_material_molecular_bulk_modulus() {
    let law = material(true, false);
    let modes = [ThermalMode {
        weight: 1.0,
        rate_per_s: 1e-7,
    }];
    let result = realize(&law, &modes, [0.1, 1000.0], 2, &mut || false).unwrap();
    for omega in [0.1, 3.0, 40.0, 1000.0] {
        let expected = law.caloric_response(omega).unwrap().bulk_modulus_pa.conj();
        complex_close(result.bulk_at(law.gas().pressure, omega), expected, 2e-6);
    }
}

#[test]
fn step_compression_matches_original_temperature_ode_storage_and_heat() {
    let law = material(true, false);
    let gas = law.gas();
    let lambda = 35.0;
    let result = realize(
        &law,
        &[ThermalMode {
            weight: 1.0,
            rate_per_s: lambda,
        }],
        [0.01, 1e4],
        2,
        &mut || false,
    )
    .unwrap();
    let capacity = law.processes()[0].heat_capacity_j_kg_k;
    let frozen_cv = law.frozen_cv_j_kg_k();
    let exchange = capacity / law.processes()[0].relaxation_time_s;
    let conduction = gas.specific_heat_cp * lambda;
    let strain = 1e-5;
    // Primary temperature, original molecular temperature, integrated heat
    // density. Following the instantaneous strain, molecular temperature is
    // still unchanged and only translational/frozen energy has responded.
    let mut x = [
        -gas.specific_gas_constant * gas.temperature * strain / frozen_cv,
        0.0,
        0.0,
    ];
    let derivative = |x: &[f64; 3]| {
        [
            (-(conduction + exchange) * x[0] + exchange * x[1]) / frozen_cv,
            exchange * (x[0] - x[1]) / capacity,
            gas.density / gas.temperature
                * (conduction * x[0] * x[0] + exchange * (x[0] - x[1]).powi(2)),
        ]
    };
    let dt = 1e-5;
    for tick in 1..=10000 {
        let stage = |d: &[f64; 3], factor: f64| {
            [
                x[0] + factor * dt * d[0],
                x[1] + factor * dt * d[1],
                x[2] + factor * dt * d[2],
            ]
        };
        let a = derivative(&x);
        let b = derivative(&stage(&a, 0.5));
        let c = derivative(&stage(&b, 0.5));
        let d = derivative(&stage(&c, 1.0));
        for i in 0..3 {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
        if tick % 1000 != 0 {
            continue;
        }
        let time = tick as f64 * dt;
        let mut modulus = gas.pressure;
        let mut storage = 0.5 * gas.pressure * strain * strain;
        let mut heat = 0.0;
        for arm in &result.arms {
            let decay = det::exp(-time / arm.relaxation_time_s);
            modulus += arm.increment_pa * decay;
            storage += 0.5 * arm.increment_pa * strain * strain * decay * decay;
            heat += 0.5 * arm.increment_pa * strain * strain * (1.0 - decay * decay);
        }
        let pressure = -gas.pressure * strain + gas.density * gas.specific_gas_constant * x[0];
        let availability = 0.5 * gas.pressure * strain * strain
            + gas.density / (2.0 * gas.temperature)
                * (frozen_cv * x[0] * x[0] + capacity * x[1] * x[1]);
        close(pressure, -strain * modulus, 2e-11);
        close(availability, storage, 2e-11);
        close(x[2], heat, 2e-10);
        close(
            storage + heat,
            0.5 * result.frozen_bulk_modulus_pa * strain * strain,
            3e-12,
        );
    }
}

#[test]
fn admission_bounds_joint_coordinates_and_polls_before_partial_publication() {
    let law = material(true, false);
    let modes = [ThermalMode {
        weight: 1.0,
        rate_per_s: 35.0,
    }];
    assert!(realize(&law, &modes, [1.0, 1000.0], 1, &mut || false).is_err());
    assert!(
        realize(
            &law,
            &modes,
            [1.0, 1000.0],
            MAX_COORDINATES + 1,
            &mut || false
        )
        .is_err()
    );
    assert!(
        realize(
            &law,
            &[ThermalMode {
                weight: 0.9,
                rate_per_s: 35.0
            }],
            [1.0, 1000.0],
            2,
            &mut || false
        )
        .is_err()
    );
    assert!(realize(&law, &modes, [0.0, 1000.0], 2, &mut || false).is_err());
    let mut polls = 0;
    let result = realize(&law, &modes, [1.0, 1000.0], 2, &mut || {
        polls += 1;
        polls > 4
    });
    assert_eq!(
        result.unwrap_err(),
        "molecular cavity thermal preparation cancelled"
    );
}
