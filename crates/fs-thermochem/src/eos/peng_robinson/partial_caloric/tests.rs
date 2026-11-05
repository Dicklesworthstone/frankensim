use super::super::{PengRobinsonComponentV1, PengRobinsonPhaseRootV1, PengRobinsonPhaseV1};
use super::*;
use fs_qty::{MolarMass, Pressure};

fn sources() -> [PengRobinsonComponentV1; 2] {
    [
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
            molar_mass: MolarMass::new(0.028),
        },
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(150.687),
            critical_pressure: Pressure::new(4.863e6),
            acentric_factor: -0.00219,
            molar_mass: MolarMass::new(0.004),
        },
    ]
}
fn binary() -> PengRobinsonMixtureV1 {
    PengRobinsonMixtureV1::new(
        sources().to_vec(),
        vec![0.75, 0.25],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap()
}
const RHO: f64 = 47.889_003_167_633_395;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}, relative={relative:e}"
    );
}
fn euler(model: &PengRobinsonMixtureV1, t: f64, rho: f64) {
    let evaluated = model
        .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
        .unwrap();
    assert!(core::ptr::eq(evaluated.mixture(), model));
    assert_eq!(
        *evaluated.state(),
        model
            .evaluate_density_caloric(Temperature::new(t), Density::new(rho))
            .unwrap()
    );
    let h = model
        .mole_fractions()
        .iter()
        .zip(evaluated.components())
        .map(|(x, value)| x * value.residual_enthalpy_j_mol())
        .sum::<f64>();
    let v = model
        .mole_fractions()
        .iter()
        .zip(evaluated.components())
        .map(|(x, value)| x * value.molar_volume_m3_mol())
        .sum::<f64>();
    close(
        h,
        evaluated.state().residual_enthalpy_j_mol(),
        2048.0 * f64::EPSILON,
    );
    close(
        v,
        evaluated.state().molar_volume_m3_per_mol(),
        2048.0 * f64::EPSILON,
    );
}

#[test]
fn partial_calorics_match_independent_extensive_derivatives() {
    // Independent 90-digit Decimal evaluation: perturb one amount by +/-1e-20
    // mol, renormalize x_i=n_i/sum(n), solve the dimensional PR pressure equation
    // at unchanged T,p, then differentiate total N*H^R and total N*v. Source
    // constants and the initial density are exact f64 inputs. This computes
    // the defining extensive derivatives, not the analytic production formula.
    let model = binary();
    for (t, rho, expected) in [
        (
            300.0,
            RHO,
            [
                (-891.758_133_724_153_1, 4.506_224_758_994_829_4e-4),
                (-354.790_656_400_891_84, 4.857_151_941_815_638e-4),
            ],
        ),
        // One original component's alpha base has changed sign here.
        (
            2200.0,
            50.0,
            [
                (1137.387_737_574_336_3, 4.416_830_191_014_400_8e-4),
                (920.326_532_657_414_1, 4.349_509_426_956_799_4e-4),
            ],
        ),
    ] {
        let evaluated = model
            .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
            .unwrap();
        for (actual, (h, v)) in evaluated.components().iter().zip(expected) {
            close(actual.residual_enthalpy_j_mol(), h, 2e-12);
            close(actual.molar_volume_m3_mol(), v, 2e-12);
        }
        assert!(
            (evaluated.components()[0].residual_enthalpy_j_mol()
                - evaluated.state().residual_enthalpy_j_mol())
            .abs()
                > 50.0
        );
        euler(&model, t, rho);
    }
}

fn nearby_root(
    model: &PengRobinsonMixtureV1,
    t: f64,
    p: f64,
    volume: f64,
) -> PengRobinsonPhaseRootV1 {
    let state = model
        .evaluate_pt(Temperature::new(t), Pressure::new(p))
        .unwrap();
    state
        .phase_roots()
        .iter()
        .min_by(|a, b| {
            (a.molar_volume_m3_per_mol() - volume)
                .abs()
                .total_cmp(&(b.molar_volume_m3_per_mol() - volume).abs())
        })
        .unwrap()
        .clone()
}
fn derivative(mut evaluate: impl FnMut(f64) -> f64, step: f64) -> f64 {
    (evaluate(-2.0 * step) - 8.0 * evaluate(-step) + 8.0 * evaluate(step) - evaluate(2.0 * step))
        / (12.0 * step)
}

#[test]
fn partial_responses_satisfy_gibbs_helmholtz_and_fixed_pressure_amount_derivatives() {
    let model = binary();
    for (t, rho) in [(300.0, RHO), (360.0, 80.0), (2200.0, 50.0), (4000.0, 30.0)] {
        let evaluated = model
            .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
            .unwrap();
        let p = evaluated.state().pressure().value();
        let volume = evaluated.state().molar_volume_m3_per_mol();
        for (i, expected) in evaluated.components().iter().enumerate() {
            // Independent original fugacity path: hbar_i^R=-R*T²*d ln(phi_i)/dT|p.
            let dh = derivative(
                |dt| nearby_root(&model, t + dt, p, volume).log_fugacity_coefficients()[i],
                t * 1e-4,
            );
            close(expected.residual_enthalpy_j_mol(), -R * t * t * dh, 2e-7);
            let dv = derivative(
                |dp| nearby_root(&model, t, p + dp, volume).log_fugacity_coefficients()[i],
                p * 1e-4,
            );
            close(expected.molar_volume_m3_mol(), R * t * (1.0 / p + dv), 2e-7);
            for (property, target) in [
                (0, expected.residual_enthalpy_j_mol()),
                (1, expected.molar_volume_m3_mol()),
            ] {
                let response = derivative(
                    |dn| {
                        let mut amounts = model.mole_fractions().to_vec();
                        amounts[i] += dn;
                        let total = amounts.iter().sum::<f64>();
                        let changed = PengRobinsonMixtureV1::new(
                            model.components().to_vec(),
                            amounts.iter().map(|n| n / total).collect(),
                            model.binary_interactions().to_vec(),
                        )
                        .unwrap();
                        let root = nearby_root(&changed, t, p, volume);
                        total
                            * if property == 0 {
                                root.residual_enthalpy_j_mol()
                            } else {
                                root.molar_volume_m3_per_mol()
                            }
                    },
                    1e-4,
                );
                close(response, target, 2e-7);
            }
        }
        euler(&model, t, rho);
    }
}

#[test]
fn partial_calorics_keep_pure_limits_signed_volumes_and_original_order() {
    let pure = PengRobinsonMixtureV1::new(vec![sources()[0]], vec![1.0], vec![]).unwrap();
    for (t, rho) in [(300.0, RHO), (300.0, 1e-15), (4000.0, 30.0)] {
        let evaluated = pure
            .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
            .unwrap();
        assert_eq!(evaluated.components().len(), 1);
        assert_eq!(
            evaluated.components()[0]
                .residual_enthalpy_j_mol()
                .to_bits(),
            evaluated.state().residual_enthalpy_j_mol().to_bits()
        );
        assert_eq!(
            evaluated.components()[0].molar_volume_m3_mol().to_bits(),
            evaluated.state().molar_volume_m3_per_mol().to_bits()
        );
    }
    let model = binary();
    let signed = model
        .evaluate_density_partial_calorics(Temperature::new(100.0), Density::new(40.5))
        .unwrap();
    assert!(signed.components()[0].molar_volume_m3_mol() < 0.0);
    assert!(signed.components()[1].molar_volume_m3_mol() > 0.0);
    euler(&model, 100.0, 40.5);
    let reversed = PengRobinsonMixtureV1::new(
        sources().into_iter().rev().collect(),
        vec![0.25, 0.75],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    for (t, rho) in [(300.0, RHO), (2200.0, 50.0)] {
        let left = model
            .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
            .unwrap();
        let right = reversed
            .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
            .unwrap();
        for (left, right) in left
            .components()
            .iter()
            .zip(right.components().iter().rev())
        {
            close(
                left.residual_enthalpy_j_mol(),
                right.residual_enthalpy_j_mol(),
                2e-13,
            );
            close(
                left.molar_volume_m3_mol(),
                right.molar_volume_m3_mol(),
                2e-13,
            );
        }
    }
}

#[test]
fn dilute_partial_enthalpy_retains_each_species_second_virial() {
    let model = binary();
    let t = 300.0;
    let mixed = model.mixing(t).unwrap();
    let mut derivative_rows = [0.0; 2];
    model
        .mixing_caloric(t, |record| {
            if let MixingRecord::Pair {
                row,
                column,
                derivative,
                ..
            } = record
            {
                derivative_rows[row] += model.mole_fractions()[column] * derivative;
            }
        })
        .unwrap();
    for rho in [1e-6, 1e-12, 1e-18] {
        let evaluated = model
            .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
            .unwrap();
        let v = evaluated.state().molar_volume_m3_per_mol();
        for (i, actual) in evaluated.components().iter().enumerate() {
            let attraction_row = (0..2)
                .map(|j| model.mole_fractions()[j] * mixed.aij[2 * i + j])
                .sum::<f64>();
            // Fixed-p second-virial limit from -R*T²*d ln(phi_i)/dT|p.
            let expected = (R * t * mixed.bi[i]
                + t * (2.0 * derivative_rows[i] - mixed.caloric.da_mix_dt)
                - 2.0 * (2.0 * attraction_row - mixed.caloric.a_mix))
                / v;
            assert_ne!(actual.residual_enthalpy_j_mol(), 0.0);
            close(actual.residual_enthalpy_j_mol(), expected, 2e-8);
        }
        euler(&model, t, rho);
    }
}

#[test]
fn partial_calorics_keep_the_original_refusals_and_component_bound() {
    let model = binary();
    for (t, rho) in [
        (300.0, 0.0),
        (300.0, -1.0),
        (300.0, f64::NAN),
        (300.0, f64::INFINITY),
        (0.0, RHO),
        (-1.0, RHO),
        (f64::NAN, RHO),
        (f64::INFINITY, RHO),
        (300.0, 2000.0),
    ] {
        assert_eq!(
            model
                .evaluate_density_partial_calorics(Temperature::new(t), Density::new(rho))
                .unwrap_err(),
            model
                .evaluate_density_caloric(Temperature::new(t), Density::new(rho))
                .unwrap_err()
        );
    }
    let pure = PengRobinsonMixtureV1::new(vec![sources()[0]], vec![1.0], vec![]).unwrap();
    let phase = pure
        .evaluate_pt(Temperature::new(150.0), Pressure::new(1e6))
        .unwrap();
    let unstable = phase.phase_roots()[1].mass_density_kg_per_m3();
    assert_eq!(
        pure.evaluate_density_partial_calorics(Temperature::new(150.0), Density::new(unstable))
            .unwrap_err(),
        PengRobinsonErrorV1::MechanicallyUnstable
    );
    for branch in [PengRobinsonPhaseV1::Liquid, PengRobinsonPhaseV1::Vapor] {
        euler(&pure, 150.0, phase.select(branch).mass_density_kg_per_m3());
    }
    let zero = PengRobinsonMixtureV1::new(sources().to_vec(), vec![1.0, 0.0], vec![]).unwrap();
    let zero = zero
        .evaluate_density_partial_calorics(Temperature::new(300.0), Density::new(RHO))
        .unwrap();
    assert_eq!(zero.components().len(), 2);
    assert!(zero.components()[1].residual_enthalpy_j_mol().is_finite());
    let maximum = PengRobinsonMixtureV1::new(
        vec![sources()[0]; MAX_PENG_ROBINSON_COMPONENTS_V1],
        vec![1.0 / MAX_PENG_ROBINSON_COMPONENTS_V1 as f64; MAX_PENG_ROBINSON_COMPONENTS_V1],
        vec![],
    )
    .unwrap();
    assert_eq!(
        maximum
            .evaluate_density_partial_calorics(Temperature::new(300.0), Density::new(RHO))
            .unwrap()
            .components()
            .len(),
        MAX_PENG_ROBINSON_COMPONENTS_V1
    );
    euler(&maximum, 300.0, RHO);
}
