use super::*;
use crate::{MolarThermalQuantityV1, SpeciesId};

fn model(
    a: &[(usize, usize, f64)],
    b: &[(usize, usize, f64)],
    cp: &[[f64; 2]],
    datum: f64,
) -> MulticomponentLiquidCaloricModelV1 {
    let n = cp.len();
    let ids: Vec<_> = (0..n)
        .map(|i| SpeciesId::new(format!("T{i:03}")).unwrap())
        .collect();
    let matrix = |pairs: &[(usize, usize, f64)]| {
        let mut values = vec![0.0; n * n];
        for &(i, j, value) in pairs {
            values[i * n + j] = value;
            values[j * n + i] = value;
        }
        values
    };
    let potential = MulticomponentRegularSolutionV1::new(
        ids.clone(),
        matrix(a),
        matrix(b),
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    MulticomponentLiquidCaloricModelV1::new(
        potential,
        ids.into_iter()
            .zip(cp)
            .map(|(id, &[intercept, slope])| {
                PureLiquidHeatCapacityV1::new(
                    id,
                    Temperature::new(datum),
                    [100.0, 1000.0_f64.max(datum)]
                        .into_iter()
                        .map(|t| {
                            (
                                Temperature::new(t),
                                MolarThermalQuantityV1::new(intercept + slope * t),
                            )
                        })
                        .collect(),
                )
                .unwrap()
            })
            .collect(),
    )
    .unwrap()
}
fn controls() -> MulticomponentLiquidTransportSolveV1 {
    MulticomponentLiquidTransportSolveV1 {
        temperature_bracket: [Temperature::new(200.0), Temperature::new(800.0)],
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
fn state(
    model: &MulticomponentLiquidCaloricModelV1,
    t: f64,
    amounts: &[f64],
) -> MulticomponentCaloricStateV1 {
    model
        .homogeneous(Temperature::new(t), amounts, || false)
        .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn totals(
    t: f64,
    n: &[f64],
    cp: &[[f64; 2]],
    a: &[(usize, usize, f64)],
    b: &[(usize, usize, f64)],
) -> [f64; 2] {
    let total: f64 = n.iter().sum();
    let h: f64 = n
        .iter()
        .zip(cp)
        .map(|(&n, &[a, b])| n * (a * (t - 300.0) + 0.5 * b * (t * t - 300.0 * 300.0)))
        .sum();
    let s: f64 = n
        .iter()
        .zip(cp)
        .map(|(&n, &[a, b])| n * (a * (t / 300.0).ln() + b * (t - 300.0)))
        .sum();
    let pair = |pairs: &[(usize, usize, f64)]| -> f64 {
        pairs
            .iter()
            .map(|&(i, j, coefficient)| coefficient * n[i] * n[j] / total)
            .sum()
    };
    let ideal: f64 = n.iter().map(|&v| v * (v / total).ln()).sum();
    [h + R * pair(b), s - R * (ideal + pair(a))]
}

#[test]
fn signed_species_exchange_recovers_original_variable_cp_temperature_and_entropy() {
    let cp = [[20.0, 0.02], [40.0, 0.04], [60.0, 0.01]];
    let a = [(0, 1, 0.4), (0, 2, -0.2), (1, 2, 0.1)];
    let b = [(0, 1, 90.0), (0, 2, -20.0), (1, 2, 30.0)];
    let model = model(&a, &b, &cp, 300.0);
    let old = [0.7, 1.2, 0.6];
    let new = [0.9, 0.8, 0.5];
    let initial = state(&model, 420.0, &old);
    let final_state = state(&model, 460.0, &new);
    let before = totals(420.0, &old, &cp, &a, &b);
    let after = totals(460.0, &new, &cp, &a, &b);
    let change = model
        .homogeneous_transport_change(&initial, &final_state, || false)
        .unwrap();
    near(change.enthalpy_change_j, after[0] - before[0], 2e-13, 1e-9);
    near(
        change.entropy_change_j_per_k,
        after[1] - before[1],
        2e-13,
        1e-11,
    );
    let result = model
        .homogeneous_transport_candidate(&initial, &new, after[0] - before[0], controls(), || false)
        .unwrap();
    near(result.final_state().temperature().value(), 460.0, 0.0, 1e-9);
    near(
        result.change().entropy_change_j_per_k,
        change.entropy_change_j_per_k,
        2e-12,
        1e-10,
    );
    assert_eq!(result.final_state().component_moles(), &new);
    near(
        result.energy_residual_j(),
        0.0,
        0.0,
        controls().energy_tolerance_j,
    );

    // Exporting species can decrease a cell's entropy; only a closed spatial
    // transaction may impose a nonnegative total production condition.
    let depleted: Vec<_> = old.iter().map(|n| 0.5 * n).collect();
    let exported = model
        .homogeneous_transport_candidate(
            &initial,
            &depleted,
            -0.5 * initial.enthalpy_j(),
            controls(),
            || false,
        )
        .unwrap();
    near(
        exported.final_state().temperature().value(),
        420.0,
        0.0,
        1e-9,
    );
    assert!(exported.change().entropy_change_j_per_k < 0.0);
}

#[test]
fn caloric_increments_survive_large_datums_and_quadratic_mixing() {
    let huge = model(&[], &[], &[[20.0, 0.0], [30.0, 0.0]], 1e15);
    let initial = state(&huge, 400.0, &[2.0, 3.0]);
    let final_state = state(&huge, 400.0 + 1e-9, &[2.0, 3.0]);
    assert_eq!(initial.enthalpy_j(), final_state.enthalpy_j());
    let change = huge
        .homogeneous_transport_change(&initial, &final_state, || false)
        .unwrap();
    near(
        change.enthalpy_change_j,
        130.0 * (final_state.temperature().value() - initial.temperature().value()),
        1e-14,
        0.0,
    );
    assert!(change.entropy_change_j_per_k > 0.0);

    let ideal = model(&[], &[], &[[20.0, 0.0]; 2], 300.0);
    let e = 2.0_f64.powi(-30);
    let from = state(&ideal, 300.0, &[1.0, 1.0]);
    let to = state(&ideal, 300.0, &[1.0 + e, 1.0 - e]);
    let change = ideal
        .homogeneous_transport_change(&from, &to, || false)
        .unwrap();
    near(change.entropy_change_j_per_k, -R * e * e, 2e-14, 0.0);
}

#[test]
fn finite_face_divided_enthalpy_matches_affine_cp_and_gibbs_helmholtz_identity() {
    let cp = [[20.0, 0.02], [30.0, 0.04], [60.0, 0.01]];
    let a = [(0, 1, 0.3), (0, 2, -0.2), (1, 2, 0.1)];
    let b = [(0, 1, 90.0), (0, 2, -50.0), (1, 2, 20.0)];
    let model = model(&a, &b, &cp, 300.0);
    let temperatures = [320.0_f64, 610.0];
    let x = [[0.2, 0.3, 0.5], [0.4, 0.35, 0.25]];
    let points = [0, 1].map(|i| MulticomponentLiquidTransportPointV1 {
        temperature: Temperature::new(temperatures[i]),
        mole_fractions: &x[i],
    });
    let face = model.transport_efforts(points, || false).unwrap();
    let middle: Vec<_> = (0..3).map(|i| 0.5 * (x[0][i] + x[1][i])).collect();
    let t_star = temperatures[0] * temperatures[1] * (temperatures[1] / temperatures[0]).ln()
        / (temperatures[1] - temperatures[0]);
    let endpoint_b = |y: &[f64]| -> f64 { b.iter().map(|&(i, j, v)| v * y[i] * y[j]).sum() };
    let beta = 1.0 / temperatures[1] - 1.0 / temperatures[0];
    let mut gauges = Vec::new();
    for i in 0..3 {
        let row: f64 = b
            .iter()
            .map(|&(j, k, v)| {
                if j == i {
                    v * middle[k]
                } else if k == i {
                    v * middle[j]
                } else {
                    0.0
                }
            })
            .sum();
        let expected_h = cp[i][0] * (t_star - 300.0)
            + 0.5 * cp[i][1] * (temperatures[0] * temperatures[1] - 300.0 * 300.0)
            + R * (row - 0.5 * (endpoint_b(&x[0]) + endpoint_b(&x[1])));
        near(
            face.partial_enthalpies_j_per_mol[i],
            expected_h,
            4e-13,
            1e-9,
        );
        let mu_over_t = |side: usize| {
            let t = temperatures[side];
            let pure_h = cp[i][0] * (t - 300.0) + 0.5 * cp[i][1] * (t * t - 90000.0);
            let pure_s = cp[i][0] * (t / 300.0).ln() + cp[i][1] * (t - 300.0);
            let state = model
                .thermodynamics()
                .evaluate(Temperature::new(t), &x[side], || false)
                .unwrap();
            pure_h / t - pure_s + R * state.reduced_chemical_potentials()[i].unwrap()
        };
        gauges.push(
            mu_over_t(1)
                - mu_over_t(0)
                - face.chemical_efforts_j_per_mol_k[i]
                - face.partial_enthalpies_j_per_mol[i] * beta,
        );
    }
    near(gauges[0], gauges[1], 0.0, 2e-12);
    near(gauges[0], gauges[2], 0.0, 2e-12);
    let reverse = model
        .transport_efforts([points[1], points[0]], || false)
        .unwrap();
    for i in 0..3 {
        assert_eq!(
            reverse.chemical_efforts_j_per_mol_k[i],
            -face.chemical_efforts_j_per_mol_k[i]
        );
        near(
            reverse.partial_enthalpies_j_per_mol[i],
            face.partial_enthalpies_j_per_mol[i],
            2e-14,
            1e-10,
        );
    }
    let uniform = model
        .transport_efforts(
            [
                points[0],
                MulticomponentLiquidTransportPointV1 {
                    mole_fractions: &x[0],
                    ..points[1]
                },
            ],
            || false,
        )
        .unwrap();
    assert_eq!(uniform.chemical_efforts_j_per_mol_k, vec![0.0; 3]);
}

#[test]
fn abundant_carrier_effort_retains_independent_positive_trace() {
    let model = model(&[], &[], &[[20.0, 0.0]; 2], 300.0);
    let x = [[1.0, 1e-100], [1.0, 2e-100]];
    let points = [0, 1].map(|i| MulticomponentLiquidTransportPointV1 {
        temperature: Temperature::new(400.0),
        mole_fractions: &x[i],
    });
    let face = model.transport_efforts(points, || false).unwrap();
    near(
        face.chemical_efforts_j_per_mol_k[0],
        -R * 1e-100,
        1e-14,
        0.0,
    );
    near(
        face.chemical_efforts_j_per_mol_k[1],
        R * 2.0_f64.ln(),
        1e-14,
        0.0,
    );
    let before = state(&model, 300.0, &x[0]);
    let after = state(&model, 300.0, &x[1]);
    let change = model
        .homogeneous_transport_change(&before, &after, || false)
        .unwrap();
    let expected = -R * 1e-100 * (1e-100_f64.ln() + 2.0 * 2.0_f64.ln() - 1.0);
    near(change.entropy_change_j_per_k, expected, 2e-14, 0.0);
}

#[test]
fn finite_face_enthalpy_contrasts_cancel_a_huge_common_pure_datum_before_rounding() {
    let b = [(0, 1, 3.0), (0, 2, -7.0), (1, 2, 5.0)];
    let x = [[0.2, 0.3, 0.5], [0.4, 0.35, 0.25]];
    let mean: Vec<_> = (0..3).map(|i| 0.5 * (x[0][i] + x[1][i])).collect();
    let row = |i| -> f64 {
        b.iter()
            .map(|&(j, k, coefficient)| {
                if j == i {
                    coefficient * mean[k]
                } else if k == i {
                    coefficient * mean[j]
                } else {
                    0.0
                }
            })
            .sum()
    };
    for datum in [300.0, 1e15] {
        let model = model(&[], &b, &[[80.0, 0.0]; 3], datum);
        let points = [0, 1].map(|i| MulticomponentLiquidTransportPointV1 {
            temperature: Temperature::new([320.0, 610.0][i]),
            mole_fractions: &x[i],
        });
        let face = model.transport_efforts(points, || false).unwrap();
        for i in 0..3 {
            near(
                face.partial_enthalpy_contrasts_j_per_mol[i],
                R * (row(i) - row(0)),
                3e-13,
                1e-12,
            );
        }
    }
}

#[test]
fn global_bracket_admission_checks_all_tangent_directions_and_both_temperatures() {
    let cp = [[20.0, 0.0]; 3];
    let good = model(&[(0, 1, 0.4)], &[(0, 2, 30.0)], &cp, 300.0);
    assert!(
        good.admit_transport_temperature_interval(controls().temperature_bracket, || false)
            .unwrap()
            > 0.0
    );
    // This original law is locally stable at uniform x, but not globally
    // convex: at x=[.25,.25,.5], v=[1,1,-2] has negative curvature.
    let attractive = model(&[(0, 1, -8.5)], &[], &cp, 300.0);
    assert!(matches!(
        attractive.admit_transport_temperature_interval(controls().temperature_bracket, || false),
        Err(MulticomponentCaloricErrorV1::TransportConvexity { .. })
    ));
    let cold_spinodal = model(&[], &[(0, 1, 500.0)], &cp[..2], 300.0);
    assert!(
        cold_spinodal
            .admit_transport_temperature_interval(
                [Temperature::new(400.0), Temperature::new(800.0)],
                || false,
            )
            .is_ok()
    );
    assert!(matches!(
        cold_spinodal
            .admit_transport_temperature_interval(controls().temperature_bracket, || false),
        Err(MulticomponentCaloricErrorV1::TransportConvexity { .. })
    ));
}

#[test]
fn signed_transport_refuses_changed_owner_absence_unsupported_heat_and_final_cancel() {
    let model = model(&[], &[], &[[20.0, 0.0]; 2], 300.0);
    let initial = state(&model, 400.0, &[0.3, 0.7]);
    let input_copy = initial.clone();
    assert!(
        model
            .homogeneous_transport_candidate(&initial, &[0.0, 1.0], 0.0, controls(), || false)
            .is_err()
    );
    assert!(matches!(
        model.homogeneous_transport_candidate(&initial, &[0.3, 0.7], 1e9, controls(), || false),
        Err(MulticomponentCaloricErrorV1::Caloric(
            BinaryCaloricErrorV1::HeatOutsideBracket { .. }
        ))
    ));
    let wrong = super::tests::model(&[], &[], &[[30.0, 0.0]; 2], 300.0);
    assert!(
        wrong
            .homogeneous_transport_candidate(&initial, &[0.3, 0.7], 0.0, controls(), || false)
            .is_err()
    );
    let mut count = 0;
    let result = model
        .homogeneous_transport_candidate(&initial, &[0.4, 0.6], 40.0, controls(), || {
            count += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    let cancelled =
        model.homogeneous_transport_candidate(&initial, &[0.4, 0.6], 40.0, controls(), || {
            polls += 1;
            polls == count
        });
    assert_eq!(cancelled, Err(MulticomponentCaloricErrorV1::Cancelled));
    assert_eq!(initial, input_copy);
    assert_eq!(
        result,
        model
            .homogeneous_transport_candidate(&initial, &[0.4, 0.6], 40.0, controls(), || false)
            .unwrap()
    );
}
