//! Independent original-source entropy-pair and inviscid face regressions.

use super::EntropyFaceReport;
use crate::multicomponent::tests::{bank, source};
use crate::multicomponent::{CaloricDomain, Conserved, Error, MixtureGas};
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

// Rows are [Cp/R constant, Cp/R linear coefficient, M/R]. The oracle
// independently integrates these polynomials instead of calling the flux mean.
const ROWS: [[f64; 3]; 2] = [[3.5, 0.2, 1.0], [4.0, 0.4, 2.0]];
const T0: f64 = 0.001;

fn polynomial_gas(formation: f64) -> MixtureGas {
    MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.2, formation),
            source("B", 2.0 * R, 4.0, 0.4, -2.0 * formation),
        ]),
        CaloricDomain::new(Temperature::new(T0), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

fn state(gas: &MixtureGas, rho: &[f64], velocity: f64, t: f64) -> Conserved {
    let concentration: f64 = rho
        .iter()
        .zip(gas.source().models())
        .map(|(density, model)| density / model.molar_mass().value())
        .sum();
    gas.conserved(rho, velocity, R * t * concentration).unwrap()
}

fn values(q: &Conserved) -> Vec<f64> {
    let mut result = q.partial_densities.clone();
    result.extend([q.momentum, q.energy]);
    result
}

fn evaluate(
    gas: &MixtureGas,
    left: &Conserved,
    right: &Conserved,
    dissipative: bool,
) -> (Vec<f64>, EntropyFaceReport) {
    let mut out = vec![0.0; gas.species_count() + 2];
    let report = if dissipative {
        gas.entropy_stable_flux(left.as_cell(), right.as_cell(), &mut out, || false)
    } else {
        gas.entropy_conservative_flux(left.as_cell(), right.as_cell(), &mut out, || false)
    }
    .unwrap();
    (out, report)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-11 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

// Mathematical entropy is minus physical mixture entropy. Each entry follows
// a direct derivative with respect to [rho_A,rho_B,momentum,referenced energy].
fn entropy_variables(rho: &[f64], u: f64, t: f64) -> [f64; 4] {
    let mut v = [0.0; 4];
    for k in 0..2 {
        let [a, b, mass] = ROWS[k];
        let e = ((a - 1.0) * (t - T0) + 0.5 * b * (t * t - T0 * T0)) / mass;
        let s = ((a - 1.0) * t.ln() + b * t - (rho[k] / mass).ln()) / mass;
        v[k] = e / t + 1.0 / mass - s - 0.5 * u * u / t;
    }
    v[2] = u / t;
    v[3] = -1.0 / t;
    v
}

fn potential(rho: &[f64], u: f64) -> f64 {
    u * (rho[0] / ROWS[0][2] + rho[1] / ROWS[1][2])
}

fn mean_internal(row: [f64; 3], left: f64, right: f64) -> f64 {
    let [a, b, mass] = row;
    let divided_t = if left == right {
        left
    } else {
        left * right * ((right - left) / left).ln_1p() / (right - left)
    };
    ((a - 1.0) * divided_t + 0.5 * b * left * right - (a - 1.0) * T0 - 0.5 * b * T0 * T0) / mass
}

#[test]
fn g0_variable_source_entropy_pair_and_dissipation_match_independent_derivatives() {
    let gas = polynomial_gas(0.0);
    let left = state(&gas, &[0.8, 0.2], -0.3, 2.0);
    let right = state(&gas, &[0.3, 0.7], 1.2, 7.0);
    let pl = gas.primitive(left.as_cell(), || false).unwrap();
    let pr = gas.primitive(right.as_cell(), || false).unwrap();
    let vl = entropy_variables(&left.partial_densities, pl.velocity, pl.temperature);
    let vr = entropy_variables(&right.partial_densities, pr.velocity, pr.temperature);
    let ql = values(&left);
    let qr = values(&right);
    let jump_potential = potential(&right.partial_densities, pr.velocity)
        - potential(&left.partial_densities, pl.velocity);
    let (ec, ec_report) = evaluate(&gas, &left, &right, false);
    for (actual, expected) in ec.iter().zip([
        0.22939772576023482,
        0.1796030100332838,
        2.5951614422181946,
        4.539804789164915,
    ]) {
        near(*actual, expected);
    }
    let ec_residual = jump_potential - (0..4).map(|k| (vr[k] - vl[k]) * ec[k]).sum::<f64>();
    near(ec_residual, 0.0);
    near(ec_report.entropy_conservative_residual_w_m2_k, ec_residual);
    near(ec_report.entropy_production_w_m2_k, ec_residual);
    assert_eq!(ec_report.dissipation_entropy_w_m2_k, 0.0);

    let (es, report) = evaluate(&gas, &left, &right, true);
    let speed = (pl.velocity.abs() + pl.sound_speed).max(pr.velocity.abs() + pr.sound_speed);
    let convexity: f64 = (0..4).map(|k| (vr[k] - vl[k]) * (qr[k] - ql[k])).sum();
    let dissipation = 0.5 * speed * convexity;
    assert!(convexity > 0.0 && dissipation > 0.0);
    near(report.signal_speed, speed);
    near(dissipation, 11.711975129142397);
    near(report.dissipation_entropy_w_m2_k, dissipation);
    let production = jump_potential - (0..4).map(|k| (vr[k] - vl[k]) * es[k]).sum::<f64>();
    near(report.entropy_production_w_m2_k, production);
    near(production, ec_residual + dissipation);
    for k in 0..4 {
        near(es[k], ec[k] - 0.5 * speed * (qr[k] - ql[k]));
    }

    // Reversing trace order preserves the symmetric conservative flux, and
    // reverses the dissipative correction, not the coordinate orientation.
    let (reversed_ec, _) = evaluate(&gas, &right, &left, false);
    let (reversed_es, reversed_report) = evaluate(&gas, &right, &left, true);
    for k in 0..4 {
        near(reversed_ec[k], ec[k]);
        near(es[k] + reversed_es[k], 2.0 * ec[k]);
    }
    near(reversed_report.dissipation_entropy_w_m2_k, dissipation);
}

fn piecewise_gas(formation: f64) -> MixtureGas {
    use fs_evidence::ValidityDomain;
    use fs_thermochem::{
        MolarEnergyQuantityV1, MolarThermalQuantityV1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    };
    let lower = source("A", R, 3.5, 0.2, formation);
    let upper = source("A", R, 5.5, 0.1, formation);
    let mut first = lower.regions()[0].card().clone();
    first.validity = ValidityDomain::unconstrained().with("T", T0, 3.0);
    let mut second = upper.regions()[0].card().clone();
    second.validity = ValidityDomain::unconstrained().with("T", 3.0, 50.0);
    let continuous = Nasa9StandardStateModelV1::new(
        lower.species().clone(),
        lower.molar_mass(),
        lower.convention().clone(),
        vec![
            Nasa9RegionV1::from_card(first).unwrap(),
            Nasa9RegionV1::from_card(second).unwrap(),
        ],
    )
    .unwrap()
    .with_caloric_continuity(
        Temperature::new(2.0),
        MolarEnergyQuantityV1::new(1e6),
        MolarThermalQuantityV1::new(1e6),
    )
    .unwrap();
    MixtureGas::new(
        bank(vec![
            continuous,
            source("B", 2.0 * R, 4.0, 0.4, -2.0 * formation),
        ]),
        CaloricDomain::new(Temperature::new(T0), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_piecewise_original_internal_energy_and_large_formation_gauges_are_retained() {
    let mut reference = None;
    for formation in [0.0, 1e18] {
        let gas = piecewise_gas(formation);
        let left = state(&gas, &[0.8, 0.2], -0.3, 2.0);
        let right = state(&gas, &[0.3, 0.7], 1.2, 7.0);
        let (ec, ec_report) = evaluate(&gas, &left, &right, false);
        let (es, es_report) = evaluate(&gas, &left, &right, true);

        // Independently integrate e_ref(T) d(1/T) through the actual Cp seam.
        // The constants make e_ref continuous and retain its original Tmin datum.
        let c1 = 2.5 * T0 + 0.1 * T0 * T0;
        let c2 = c1 + 2.0 * 3.0 + 0.5 * (0.1 - 0.2) * 9.0;
        let interval = |a: f64, b: f64, c: f64, lo: f64, hi: f64| {
            -(a - 1.0) * (hi / lo).ln() - 0.5 * b * (hi - lo) - c * (1.0 / hi - 1.0 / lo)
        };
        let e_a = (interval(3.5, 0.2, c1, 2.0, 3.0) + interval(5.5, 0.1, c2, 3.0, 7.0))
            / (1.0 / 7.0 - 1.0 / 2.0);
        let e_b = mean_internal(ROWS[1], 2.0, 7.0);
        let mean_u = 0.45;
        let pressure = (0.55 + 0.45 / 2.0) / (0.5 / 2.0 + 0.5 / 7.0);
        let expected_energy =
            e_a * ec[0] + e_b * ec[1] + 0.5 * (-0.3) * 1.2 * (ec[0] + ec[1]) + mean_u * pressure;
        near(ec[3], expected_energy);
        near(ec[3], 4.820786726889179);
        near(ec_report.entropy_conservative_residual_w_m2_k, 0.0);
        assert!(es_report.entropy_production_w_m2_k > 0.0);
        // Reattach the original per-species datum and check that numerical
        // dissipation transports physical total energy in the same gauge.
        let original = |f: &[f64]| {
            f[3] + f[0] * gas.energy_references()[0] + f[1] * gas.energy_references()[1]
        };
        let original_jump = right.energy - left.energy
            + (formation + c1) * (right.partial_densities[0] - left.partial_densities[0])
            + (-formation + 0.0015001) * (right.partial_densities[1] - left.partial_densities[1]);
        near(
            original(&es),
            original(&ec) - 0.5 * es_report.signal_speed * original_jump,
        );

        if let Some((original_left, original_right, original_ec, original_es, ce, se)) = &reference
        {
            assert_eq!(&left, original_left);
            assert_eq!(&right, original_right);
            assert_eq!(&ec, original_ec);
            assert_eq!(&es, original_es);
            assert_eq!(&ec_report, ce);
            assert_eq!(&es_report, se);
        } else {
            reference = Some((left, right, ec, es, ec_report, es_report));
        }
    }
}

#[test]
fn g0_equal_and_adjacent_temperatures_keep_consistency_and_the_source_energy_mean() {
    let mut gauge_reference = None;
    for formation in [0.0, 1e18] {
        let gas = polynomial_gas(formation);
        for velocity in [0.0, 0.7] {
            let uniform = state(&gas, &[0.8, 0.2], velocity, 2.0);
            let p = gas.primitive(uniform.as_cell(), || false).unwrap();
            let physical = [
                uniform.partial_densities[0] * p.velocity,
                uniform.partial_densities[1] * p.velocity,
                uniform.momentum * p.velocity + p.pressure,
                (uniform.energy + p.pressure) * p.velocity,
            ];
            for dissipative in [false, true] {
                let (out, report) = evaluate(&gas, &uniform, &uniform, dissipative);
                assert_eq!(out, physical);
                assert_eq!(report.entropy_conservative_residual_w_m2_k, 0.0);
                assert_eq!(report.dissipation_entropy_w_m2_k, 0.0);
                assert_eq!(report.entropy_production_w_m2_k, 0.0);
            }
        }
        let stationary_left = state(&gas, &[0.8, 0.2], 0.0, 2.0);
        let stationary_right = state(&gas, &[0.3, 0.7], 0.0, 7.0);
        let (stationary, stationary_report) =
            evaluate(&gas, &stationary_left, &stationary_right, false);
        assert_eq!([stationary[0], stationary[1], stationary[3]], [0.0; 3]);
        assert!(stationary[2] > 0.0);
        assert_eq!(stationary_report.entropy_production_w_m2_k, 0.0);

        let mut weak_results = Vec::new();
        for (rho_right, t_right) in [
            ([0.3, 0.7], 2.0),
            ([0.3, 0.7], 2.0_f64.next_up()),
            ([0.8_f64.next_up(), 0.2_f64.next_up()], 2.0_f64.next_up()),
        ] {
            let left = state(&gas, &[0.8, 0.2], 0.4, 2.0);
            let right = state(&gas, &rho_right, 0.4, t_right);
            // Source-owned inversion with explicit hints retains the actual
            // adjacent temperatures before entering the native flux interface.
            let pl = gas
                .decode_cell(left.as_cell(), Some(2.0), &mut || false)
                .unwrap();
            let pr = gas
                .decode_cell(right.as_cell(), Some(t_right), &mut || false)
                .unwrap();
            assert_eq!(pl.temperature, 2.0);
            assert_eq!(pr.temperature, t_right);
            let ql = values(&left);
            let qr = values(&right);
            let mut actual = [0.0; 4];
            let report = gas
                .interface_entropy(&ql, pl, &qr, pr, &mut actual, &mut || false)
                .unwrap();
            let mut mass = [0.0; 2];
            for k in 0..2 {
                let rho_left = left.partial_densities[k];
                let delta = rho_right[k] - rho_left;
                mass[k] = 0.4 * delta / (delta / rho_left).ln_1p();
            }
            let pressure = (0.5 * (0.8 + rho_right[0]) + 0.25 * (0.2 + rho_right[1]))
                / (0.5 / 2.0 + 0.5 / t_right);
            let ec = [
                mass[0],
                mass[1],
                0.4 * (mass[0] + mass[1]) + pressure,
                mean_internal(ROWS[0], 2.0, t_right) * mass[0]
                    + mean_internal(ROWS[1], 2.0, t_right) * mass[1]
                    + 0.5 * 0.4 * 0.4 * (mass[0] + mass[1])
                    + 0.4 * pressure,
            ];
            for k in 0..4 {
                near(
                    actual[k],
                    ec[k] - 0.5 * report.signal_speed * (qr[k] - ql[k]),
                );
            }
            weak_results.push(actual);
        }
        if let Some(previous) = &gauge_reference {
            assert_eq!(&weak_results, previous);
        } else {
            gauge_reference = Some(weak_results);
        }
    }
}

#[test]
fn g4_common_absence_is_preserved_and_unsupported_support_or_cancellation_refuses() {
    let gas = MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.2, 0.0),
            source("B", 2.0 * R, 4.0, 0.4, 0.0),
            source("C", R, 5.0, 0.0, 0.0),
        ]),
        CaloricDomain::new(Temperature::new(T0), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let left = state(&gas, &[0.8, 0.2, 0.0], -0.3, 2.0);
    let right = state(&gas, &[0.3, 0.7, 0.0], 1.2, 7.0);
    let binary = polynomial_gas(0.0);
    let binary_left = state(&binary, &[0.8, 0.2], -0.3, 2.0);
    let binary_right = state(&binary, &[0.3, 0.7], 1.2, 7.0);
    for dissipative in [false, true] {
        let (out, report) = evaluate(&gas, &left, &right, dissipative);
        let (expected, expected_report) =
            evaluate(&binary, &binary_left, &binary_right, dissipative);
        assert_eq!(out[2], 0.0);
        for (actual, expected) in [out[0], out[1], out[3], out[4]].into_iter().zip(expected) {
            near(actual, expected);
        }
        near(
            report.entropy_production_w_m2_k,
            expected_report.entropy_production_w_m2_k,
        );
    }

    let absent = state(&gas, &[1.0, 0.0, 0.0], 0.0, 2.0);
    let mut out = [137.0; 5];
    let error = gas
        .entropy_stable_flux(absent.as_cell(), right.as_cell(), &mut out, || false)
        .unwrap_err();
    assert!(matches!(
        error,
        Error::Invalid("entropy mixture flux requires common positive species support")
    ));
    let error = gas
        .entropy_conservative_flux(left.as_cell(), right.as_cell(), &mut out, || true)
        .unwrap_err();
    assert_eq!(error, Error::Cancelled);
}
