use super::*;

/// Construct the complete original friction rows independently of any
/// reduced matrix or internal normalization.
fn original_forces(x: &[f64], molar: &[f64], d: &[f64], concentration: f64) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|i| {
            -sum((0..n)
                .filter(|&j| j != i)
                .map(|j| (x[j] * molar[i] - x[i] * molar[j]) / d[i * n + j]))
                / concentration
        })
        .collect()
}

#[test]
fn g1_sparse_compositions_recover_manufactured_fluxes_in_both_frames() {
    for (n, active) in [
        (3, vec![(1, 1.0_f64)]),
        (9, vec![(1, 0.25), (7, 0.75)]),
        (9, vec![(0, 0.125), (3, 0.375), (8, 0.5)]),
        (128, vec![(17, 1.0)]),
        (128, vec![(1, 0.25), (126, 0.75)]),
        (128, vec![(0, 0.125), (17, 0.125), (63, 0.25), (127, 0.5)]),
    ] {
        let masses: Vec<f64> = (0..n).map(|i| 0.002 + 0.003 * (i % 17) as f64).collect();
        let owner = model(&masses);
        let mut x = vec![0.0; n];
        for &(i, value) in &active {
            x[i] = value;
        }
        let gauge = active.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap().0;
        let mut d = matrix(n, 0.0);
        for i in 0..n {
            for j in i + 1..n {
                let value = 1e-5 * (1 + (i + 3 * j) % 11) as f64;
                d[i * n + j] = value;
                d[j * n + i] = value;
            }
        }
        let concentration = 2.5;
        let mut work = owner.workspace().unwrap();
        for mass_frame in [true, false] {
            let mut primary: Vec<f64> = (0..n)
                .map(|i| {
                    let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                    sign * (1 + i % 7) as f64 * 1e-9
                })
                .collect();
            primary[gauge] = -sum((0..n).filter(|&i| i != gauge).map(|i| primary[i]));
            let molar: Vec<f64> = primary
                .iter()
                .zip(&masses)
                .map(|(&value, &mass)| if mass_frame { value / mass } else { value })
                .collect();
            let forces = original_forces(&x, &molar, &d, concentration);
            let actual = if mass_frame {
                owner
                    .evaluate_driving_forces(
                        &x,
                        &forces,
                        AmountConcentration::new(concentration),
                        &d,
                        MaxwellStefanSolveV1::default(),
                        &mut work,
                        || false,
                    )
                    .unwrap()
                    .mass_fluxes()
                    .to_vec()
            } else {
                owner
                    .evaluate_molar_driving_forces(
                        &x,
                        &forces,
                        AmountConcentration::new(concentration),
                        &d,
                        MaxwellStefanMolarSolveV1::default(),
                        &mut work,
                        || false,
                    )
                    .unwrap()
                    .molar_fluxes()
                    .to_vec()
            };
            for (i, (&got, &expected)) in actual.iter().zip(&primary).enumerate() {
                assert!(
                    (got - expected).abs() <= 2e-11 * expected.abs().max(1e-9),
                    "n={n}, species={i}, mass_frame={mass_frame}, got={got:e}, expected={expected:e}",
                );
                if x[i] == 0.0 {
                    assert_ne!(
                        got, 0.0,
                        "an absent species still has its authored driving force"
                    );
                }
            }
        }
    }
}

#[test]
fn g3_absent_catalog_extension_preserves_the_binary_response() {
    let n = MAX_MAXWELL_STEFAN_SPECIES_V1;
    let mut masses = vec![0.044; n];
    masses[1] = 0.002;
    masses[n - 2] = 0.028;
    let owner = model(&masses);
    let mut work = owner.workspace().unwrap();
    let mut x = vec![0.0; n];
    x[1] = 0.25;
    x[n - 2] = 0.75;
    let mut driving = vec![0.0; n];
    driving[1] = 0.2;
    driving[n - 2] = -0.2;
    let c = 3.0;
    let expected =
        -c * 0.5 * 0.2 * masses[1] * masses[n - 2] / (x[1] * masses[1] + x[n - 2] * masses[n - 2]);
    for inactive_pair in [f64::from_bits(1), 1e-300, 1e200] {
        let mut d = matrix(n, inactive_pair);
        for i in 0..n {
            for j in i + 1..n {
                if x[i] > 0.0 || x[j] > 0.0 {
                    d[i * n + j] = 0.5;
                    d[j * n + i] = 0.5;
                }
            }
        }
        let result = owner
            .evaluate_driving_forces(
                &x,
                &driving,
                AmountConcentration::new(c),
                &d,
                MaxwellStefanSolveV1::default(),
                &mut work,
                || false,
            )
            .unwrap();
        near(result.mass_fluxes()[1], expected, 5e-14);
        near(result.mass_fluxes()[n - 2], -expected, 5e-14);
        for i in 0..n {
            if x[i] == 0.0 {
                assert_eq!(result.mass_fluxes()[i], 0.0);
            }
        }
        let molar = owner
            .evaluate_molar_driving_forces(
                &x,
                &driving,
                AmountConcentration::new(c),
                &d,
                MaxwellStefanMolarSolveV1::default(),
                &mut work,
                || false,
            )
            .unwrap();
        near(molar.molar_fluxes()[1], -c * 0.5 * 0.2, 5e-14);
        near(molar.molar_fluxes()[n - 2], c * 0.5 * 0.2, 5e-14);
    }
}

#[test]
fn g1_positive_trace_species_retain_their_finite_friction() {
    let owner = model(&[1.0; 3]);
    let mut work = owner.workspace().unwrap();
    for trace in [1e-100, 1e-200, 1e-300] {
        let x = [0.4, 0.6, trace];
        let d = [0.0, 1.0, trace, 1.0, 0.0, trace, trace, trace, 0.0];
        let result = owner
            .evaluate_driving_forces(
                &x,
                &[0.1, -0.1, 0.0],
                AmountConcentration::new(1.0),
                &d,
                MaxwellStefanSolveV1::default(),
                &mut work,
                || false,
            )
            .unwrap();
        // The third original row and zero-sum frame imply z2=0. The
        // trace/D ratio is one, so it doubles each carrier's friction:
        // z0=-0.1/2, whereas dropping this positive species gives -0.1.
        near(result.mass_fluxes()[0], -0.05, 3e-13);
        near(result.mass_fluxes()[1], 0.05, 3e-13);
        assert!(result.mass_fluxes()[2].abs() < 1e-310);
    }
}

#[test]
fn g1_binary_closed_form_retains_pure_limits_and_extreme_mass_ratios() {
    for masses in [[0.002, 0.028], [1e-200, 1.0], [1.0, 1e-200]] {
        let owner = model(&masses);
        let mut work = owner.workspace().unwrap();
        let d = matrix(2, 0.5);
        for x in [[0.0, 1.0], [0.25, 0.75], [1.0, 0.0]] {
            let mean = x[0] * masses[0] + x[1] * masses[1];
            let expected = -0.5 * 0.1 * masses[0] * masses[1] / mean;
            let mass = owner
                .evaluate_driving_forces(
                    &x,
                    &[0.1, -0.1],
                    AmountConcentration::new(1.0),
                    &d,
                    MaxwellStefanSolveV1::default(),
                    &mut work,
                    || false,
                )
                .unwrap();
            for (got, expected) in mass.mass_fluxes().iter().zip([expected, -expected]) {
                assert!((got / expected - 1.0).abs() < 5e-14);
            }
            let molar = owner
                .evaluate_molar_driving_forces(
                    &x,
                    &[0.1, -0.1],
                    AmountConcentration::new(1.0),
                    &d,
                    MaxwellStefanMolarSolveV1::default(),
                    &mut work,
                    || false,
                )
                .unwrap();
            near(molar.molar_fluxes()[0], -0.05, 3e-14);
            near(molar.molar_fluxes()[1], 0.05, 3e-14);
        }
    }
}

#[test]
fn g0_composition_changes_and_cancelled_solves_reuse_original_storage() {
    let owner = model(&[0.002, 0.016, 0.028, 0.044, 0.018, 0.030, 0.032, 0.040]);
    let mut work = owner.workspace().unwrap();
    let pointer = work.mass_flux.as_ptr();
    let bytes = work.workspace_bytes();
    let d = matrix(8, 0.75);
    let driving = [0.1, -0.2, 0.05, 0.0, 0.05, 0.0, 0.0, 0.0];
    for x in [
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.25, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0],
        [0.125; 8],
        [0.0, 0.0, 0.125, 0.0, 0.375, 0.0, 0.0, 0.5],
    ] {
        let mut fresh = owner.workspace().unwrap();
        let expected = owner
            .evaluate_driving_forces(
                &x,
                &driving,
                AmountConcentration::new(1.0),
                &d,
                MaxwellStefanSolveV1::default(),
                &mut fresh,
                || false,
            )
            .unwrap()
            .mass_fluxes()
            .to_vec();
        for limit in [0, 12, 23, 27] {
            let mut polls = 0;
            assert_eq!(
                owner
                    .evaluate_driving_forces(
                        &x,
                        &driving,
                        AmountConcentration::new(1.0),
                        &d,
                        MaxwellStefanSolveV1::default(),
                        &mut work,
                        || {
                            polls += 1;
                            polls > limit
                        },
                    )
                    .unwrap_err(),
                MaxwellStefanErrorV1::Cancelled,
            );
            let result = owner
                .evaluate_driving_forces(
                    &x,
                    &driving,
                    AmountConcentration::new(1.0),
                    &d,
                    MaxwellStefanSolveV1::default(),
                    &mut work,
                    || false,
                )
                .unwrap();
            assert_eq!(result.mass_fluxes().as_ptr(), pointer);
            assert!(
                result
                    .mass_fluxes()
                    .iter()
                    .zip(&expected)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
            assert_eq!(work.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g0_invalid_absent_pair_data_still_refuses() {
    let owner = model(&[0.002, 0.028, 0.044, 0.018]);
    let mut work = owner.workspace().unwrap();
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut d = matrix(4, 0.5);
        d[2 * 4 + 3] = invalid;
        d[3 * 4 + 2] = invalid;
        assert!(matches!(
            owner.evaluate_driving_forces(
                &[0.25, 0.75, 0.0, 0.0],
                &[0.1, -0.1, 0.0, 0.0],
                AmountConcentration::new(1.0),
                &d,
                MaxwellStefanSolveV1::default(),
                &mut work,
                || false,
            ),
            Err(MaxwellStefanErrorV1::Invalid(_))
        ));
    }
}
