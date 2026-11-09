use super::super::{
    AmountConcentration, Frame, MaxwellStefanMolarSolveV1, MaxwellStefanSolveV1, MaxwellStefanV1,
    MaxwellStefanWorkspaceV1, MolarMass, SpeciesId,
};
use super::*;

fn model(masses: &[f64]) -> MaxwellStefanV1 {
    MaxwellStefanV1::new(
        masses
            .iter()
            .enumerate()
            .map(|(i, &mass)| {
                (
                    SpeciesId::new(&format!("F{i:03}")).unwrap(),
                    MolarMass::new(mass),
                )
            })
            .collect(),
    )
    .unwrap()
}

fn pair_data(n: usize) -> Vec<f64> {
    let mut values = vec![0.0; n * n];
    for i in 0..n {
        for j in i + 1..n {
            let value = 1e-5 * (1 + (3 * i + 7 * j) % 17) as f64;
            values[i * n + j] = value;
            values[j * n + i] = value;
        }
    }
    values
}

fn evaluate(
    owner: &MaxwellStefanV1,
    x: &[f64],
    driving: &[f64],
    d: &[f64],
    frame: Frame,
    work: &mut MaxwellStefanWorkspaceV1,
    cancelled: impl FnMut() -> bool,
) -> Result<(Vec<f64>, Vec<f64>), MaxwellStefanErrorV1> {
    match frame {
        Frame::MassAveraged => {
            let result = owner.evaluate_driving_forces(
                x,
                driving,
                AmountConcentration::new(2.5),
                d,
                MaxwellStefanSolveV1::default(),
                work,
                cancelled,
            )?;
            Ok((
                result.mass_fluxes().to_vec(),
                result.molar_fluxes().to_vec(),
            ))
        }
        Frame::MolarAveraged => {
            let result = owner.evaluate_molar_driving_forces(
                x,
                driving,
                AmountConcentration::new(2.5),
                d,
                MaxwellStefanMolarSolveV1::default(),
                work,
                cancelled,
            )?;
            Ok((
                result.mass_fluxes().to_vec(),
                result.molar_fluxes().to_vec(),
            ))
        }
    }
}

fn same_bits(actual: &(Vec<f64>, Vec<f64>), expected: &(Vec<f64>, Vec<f64>)) {
    for (a, b) in actual
        .0
        .iter()
        .chain(&actual.1)
        .zip(expected.0.iter().chain(&expected.1))
    {
        assert_eq!(a.to_bits(), b.to_bits(), "{a:e} != {b:e}");
    }
}

#[test]
fn repeated_spatial_directions_recover_original_fluxes_and_fresh_solve_bits() {
    for x in [
        vec![0.25, 0.25, 0.5],
        vec![0.5, 0.0, 1e-250, 0.125, 0.125, 0.0, 0.0, 0.25],
        vec![1.0 / 32.0; 32],
        vec![1.0 / 128.0; 128],
    ] {
        let n = x.len();
        let masses: Vec<_> = (0..n).map(|i| 0.005 + 0.002 * (i % 13) as f64).collect();
        let owner = model(&masses);
        let d = pair_data(n);
        let gauge = x
            .iter()
            .position(|&v| v == x.iter().copied().fold(0.0, f64::max))
            .unwrap();
        let active = x.iter().filter(|&&v| v > 0.0).count();
        for frame in [Frame::MassAveraged, Frame::MolarAveraged] {
            let mut work = owner.workspace().unwrap();
            let bytes = work.workspace_bytes();
            for direction in 0..4 {
                let mut primary: Vec<_> = (0..n)
                    .map(|i| (((i + direction) % 7) as f64 - 3.0) * 1e-6)
                    .collect();
                primary[gauge] = 0.0;
                primary[gauge] = -sum(primary.iter().copied());
                let molar: Vec<_> = primary
                    .iter()
                    .zip(&masses)
                    .map(|(&value, &mass)| match frame {
                        Frame::MassAveraged => value / mass,
                        Frame::MolarAveraged => value,
                    })
                    .collect();
                // Manufacture forces from complete original pair friction,
                // independently of the reduced matrix and normalized solver.
                let driving: Vec<_> = (0..n)
                    .map(|i| {
                        -sum((0..n)
                            .filter(|&j| j != i)
                            .map(|j| (x[j] * molar[i] - x[i] * molar[j]) / d[i * n + j]))
                            / 2.5
                    })
                    .collect();
                let mut warm_polls = 0;
                let result = evaluate(&owner, &x, &driving, &d, frame, &mut work, || {
                    warm_polls += 1;
                    false
                })
                .unwrap();
                let mut cold_polls = 0;
                let cold = evaluate(
                    &owner,
                    &x,
                    &driving,
                    &d,
                    frame,
                    &mut owner.workspace().unwrap(),
                    || {
                        cold_polls += 1;
                        false
                    },
                )
                .unwrap();
                same_bits(&result, &cold);
                let actual = match frame {
                    Frame::MassAveraged => &result.0,
                    Frame::MolarAveraged => &result.1,
                };
                for (&value, &expected) in actual.iter().zip(&primary) {
                    assert!((value - expected).abs() <= 3e-10 * expected.abs() + 2e-18);
                }
                assert_eq!(work.factorization.dimension, active);
                assert_eq!(work.workspace_bytes(), bytes);
                if direction > 0 {
                    assert!(
                        warm_polls < cold_polls,
                        "repeated direction must skip elimination"
                    );
                }
            }
        }
    }
}

#[test]
fn changed_composition_coefficients_owner_and_frame_refresh_the_factorization() {
    let owner = model(&[0.005, 0.02, 0.04, 0.08, 0.12]);
    let other = model(&[0.01, 0.03, 0.06, 0.09, 0.1]);
    let mut work = owner.workspace().unwrap();
    let d = pair_data(5);
    let mut changed_d = d.clone();
    changed_d[1] *= 1.25;
    changed_d[5] = changed_d[1];
    let driving = [0.2, -0.1, 0.05, -0.07, -0.08];
    let original = [0.25, 0.125, 0.125, 0.25, 0.25];
    let changed = [0.125, 0.125, 0.5, 0.125, 0.125];
    let sparse = [0.5, 0.0, 0.25, 0.0, 0.25];
    for (source, x, pairs, frame) in [
        (&owner, &original, &d, Frame::MassAveraged),
        (&owner, &changed, &d, Frame::MassAveraged),
        (&owner, &changed, &changed_d, Frame::MassAveraged),
        (&owner, &changed, &changed_d, Frame::MolarAveraged),
        (&other, &changed, &changed_d, Frame::MassAveraged),
        (&other, &sparse, &changed_d, Frame::MassAveraged),
        (&owner, &original, &d, Frame::MassAveraged),
    ] {
        let result = evaluate(source, x, &driving, pairs, frame, &mut work, || false).unwrap();
        let cold = evaluate(
            source,
            x,
            &driving,
            pairs,
            frame,
            &mut source.workspace().unwrap(),
            || false,
        )
        .unwrap();
        same_bits(&result, &cold);
        assert_eq!(
            work.factorization.dimension,
            x.iter().filter(|&&v| v > 0.0).count()
        );
    }
}

#[test]
fn reused_factors_still_admit_current_sources_and_refresh_the_physical_scale() {
    let owner = model(&[0.005, 0.02, 0.04, 0.08, 0.12]);
    let x = [0.5, 0.25, 0.25, 0.0, 0.0];
    let driving = [0.2, -0.1, 0.05, -0.07, -0.08];
    let d = pair_data(5);
    let mut work = owner.workspace().unwrap();
    let original = evaluate(
        &owner,
        &x,
        &driving,
        &d,
        Frame::MassAveraged,
        &mut work,
        || false,
    )
    .unwrap();
    let key = work.factorization.source.clone();
    let factors = work.factorization.factors.clone();

    // A change only between two absent species cannot alter the physical
    // friction matrix, but the entire source table must still be admitted.
    let mut absent_pair = d.clone();
    absent_pair[3 * 5 + 4] = 1e-300;
    absent_pair[4 * 5 + 3] = 1e-300;
    let unchanged = evaluate(
        &owner,
        &x,
        &driving,
        &absent_pair,
        Frame::MassAveraged,
        &mut work,
        || false,
    )
    .unwrap();
    same_bits(&unchanged, &original);
    assert_eq!(work.factorization.source, key);
    assert_eq!(work.factorization.factors, factors);
    absent_pair[4 * 5 + 3] = 0.0;
    assert!(matches!(
        evaluate(
            &owner,
            &x,
            &driving,
            &absent_pair,
            Frame::MassAveraged,
            &mut work,
            || false
        ),
        Err(MaxwellStefanErrorV1::Invalid(_))
    ));

    // Uniform source scaling leaves D_min/D_ij exactly unchanged while the
    // published physical flux must use the newly admitted diffusivity scale.
    let doubled: Vec<_> = d.iter().map(|value| 2.0 * value).collect();
    let scaled = evaluate(
        &owner,
        &x,
        &driving,
        &doubled,
        Frame::MassAveraged,
        &mut work,
        || false,
    )
    .unwrap();
    for (&new, &old) in scaled
        .0
        .iter()
        .chain(&scaled.1)
        .zip(original.0.iter().chain(&original.1))
    {
        assert_eq!(new.to_bits(), (2.0 * old).to_bits());
    }
    assert_eq!(work.factorization.source, key);
    assert_eq!(work.factorization.factors, factors);
}

fn pivoted_matrix(n: usize) -> Vec<f64> {
    (0..n * n)
        .map(|index| {
            let i = (index / n + 1) % n;
            let j = index % n;
            if i == j {
                n as f64 + 1.0 + i as f64 / 8.0
            } else {
                ((i + 2 * j) % 7) as f64 / 16.0
            }
        })
        .collect()
}

fn right_for(matrix: &[f64], values: &[f64]) -> Vec<f64> {
    let n = values.len();
    (0..n)
        .map(|i| sum((0..n).map(|j| matrix[i * n + j] * values[j])))
        .collect()
}

#[test]
fn chronological_elimination_replays_late_pivots_for_multiple_right_hand_sides() {
    for n in 3..=16 {
        let a = pivoted_matrix(n);
        let mut factors = Factorization::new(n).unwrap();
        let mut scales = vec![0.0; n];
        for direction in 0..3 {
            let expected: Vec<_> = (0..n)
                .map(|i| (i + direction) as f64 / 7.0 - 1.25)
                .collect();
            let mut b = right_for(&a, &expected);
            factors
                .solve(&a, &mut b, &mut scales, n, &mut || false)
                .unwrap();
            for (&actual, &expected) in b.iter().zip(&expected) {
                assert!((actual - expected).abs() < 2e-14);
            }
            assert!(
                factors
                    .pivots
                    .iter()
                    .enumerate()
                    .skip(1)
                    .any(|(i, &p)| i != p)
            );
        }
    }
}

#[test]
fn cancelled_refactorizations_and_rhs_solves_leave_reusable_storage() {
    let n = 5;
    let a = pivoted_matrix(n);
    let expected = [0.25, -0.5, 0.75, -1.0, 1.25];
    let original_rhs = right_for(&a, &expected);
    let mut scales = vec![0.0; n];
    let mut factors = Factorization::new(n).unwrap();
    let mut solved = original_rhs.clone();
    let mut total_polls = 0;
    factors
        .solve(&a, &mut solved, &mut scales, n, &mut || {
            total_polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=total_polls {
        let mut interrupted = Factorization::new(n).unwrap();
        let mut rhs = original_rhs.clone();
        let mut polls = 0;
        assert_eq!(
            interrupted.solve(&a, &mut rhs, &mut scales, n, &mut || {
                polls += 1;
                polls == stop
            }),
            Err(MaxwellStefanErrorV1::Cancelled)
        );
        rhs.copy_from_slice(&original_rhs);
        interrupted
            .solve(&a, &mut rhs, &mut scales, n, &mut || false)
            .unwrap();
        assert_eq!(
            rhs.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            solved.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }
    let mut changed = a.clone();
    changed[0] += 0.125;
    let mut polls = 0;
    assert_eq!(
        factors.solve(
            &changed,
            &mut original_rhs.clone(),
            &mut scales,
            n,
            &mut || {
                polls += 1;
                polls == 3
            }
        ),
        Err(MaxwellStefanErrorV1::Cancelled)
    );
    assert_eq!(factors.dimension, 0);
    let mut rhs = original_rhs.clone();
    factors
        .solve(&a, &mut rhs, &mut scales, n, &mut || false)
        .unwrap();
    assert_eq!(rhs, solved);

    // Once factors are complete, cancellation only disturbs RHS scratch.
    let mut warm_polls = 0;
    factors
        .solve(&a, &mut original_rhs.clone(), &mut scales, n, &mut || {
            warm_polls += 1;
            false
        })
        .unwrap();
    assert!(warm_polls < total_polls);
    for stop in 1..=warm_polls {
        let mut rhs = original_rhs.clone();
        let mut polls = 0;
        assert_eq!(
            factors.solve(&a, &mut rhs, &mut scales, n, &mut || {
                polls += 1;
                polls == stop
            }),
            Err(MaxwellStefanErrorV1::Cancelled)
        );
        assert_eq!(factors.dimension, n);
        rhs.copy_from_slice(&original_rhs);
        factors
            .solve(&a, &mut rhs, &mut scales, n, &mut || false)
            .unwrap();
        assert_eq!(rhs, solved);
    }

    let singular = vec![0.0; n * n];
    assert_eq!(
        factors.solve(
            &singular,
            &mut original_rhs.clone(),
            &mut scales,
            n,
            &mut || false
        ),
        Err(MaxwellStefanErrorV1::SingularSystem)
    );
    assert_eq!(factors.dimension, 0);
    let mut rhs = original_rhs;
    factors
        .solve(&a, &mut rhs, &mut scales, n, &mut || false)
        .unwrap();
    assert_eq!(rhs, solved);
}

#[test]
fn workspace_budgets_include_all_factor_storage_without_growth() {
    let n = 32;
    let owner = model(&vec![0.02; n]);
    let mut work = owner.workspace().unwrap();
    let expected =
        (4 * n * n + 5 * n) * core::mem::size_of::<f64>() + n * core::mem::size_of::<usize>();
    assert_eq!(work.workspace_bytes(), expected);
    let source = work.factorization.source.as_ptr();
    let factors = work.factorization.factors.as_ptr();
    let pivots = work.factorization.pivots.as_ptr();
    let x = vec![1.0 / n as f64; n];
    let d = pair_data(n);
    for direction in 0..3 {
        let mut driving = vec![0.0; n];
        driving[direction] = 0.1;
        driving[n - 1] = -0.1;
        evaluate(
            &owner,
            &x,
            &driving,
            &d,
            Frame::MassAveraged,
            &mut work,
            || false,
        )
        .unwrap();
        assert_eq!(work.workspace_bytes(), expected);
        assert_eq!(work.factorization.source.as_ptr(), source);
        assert_eq!(work.factorization.factors.as_ptr(), factors);
        assert_eq!(work.factorization.pivots.as_ptr(), pivots);
    }
}
