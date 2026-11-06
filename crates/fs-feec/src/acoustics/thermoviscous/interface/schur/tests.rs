use super::*;

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-25),
        "{actual:?} != {expected:?}"
    );
}

/// G1: two finite-volume shear or heat coordinates in perfect contact. The
/// complex impedances retain inertia/storage and independent positive losses.
#[test]
fn two_original_impedances_share_velocity_and_exchange_equal_opposite_work() {
    let impedances = [C64::new(2.0, 7.0), C64::new(5.0, 3.0)];
    let factors = impedances.map(|a| lu_complex_with_cancel(&[a], 1, &mut || false).unwrap());
    let left_load = [2e-6];
    let right_load = [-2e-6];
    let loads = [&left_load[..], &right_load[..]];
    let mut solve = |side: usize, rhs: &[C64], cancelled: &mut _| {
        let mut x = rhs.to_vec();
        factors[side].solve_with_cancel(&mut x, cancelled)?;
        Ok(x)
    };
    let mut cancelled = || false;
    let system =
        ConformingSchur::prepare(loads, [1, 1], 1, 1e-12, &mut solve, &mut cancelled).unwrap();
    let force = C64::new(0.3, -0.7);
    let (x, traction) = system
        .response(
            loads,
            1,
            [1, 1],
            [1, 1],
            [&[force], &[C64::ZERO]],
            &mut solve,
            &mut cancelled,
        )
        .unwrap();
    let exact = force / (impedances[0] + impedances[1]);
    near(x[0][0], exact, 3e-14);
    near(x[1][0], exact, 3e-14);
    near(
        traction[0].scale(left_load[0]),
        C64::ZERO - impedances[1] * exact,
        3e-14,
    );
    let work = std::array::from_fn::<_, 2, _>(|side| {
        0.5 * (traction[0].conj() * dual(loads[side], 1, 0, &x[side])).re
    });
    assert!(work[0] < 0.0 && work[1] > 0.0);
    assert!((work[0] + work[1]).abs() < 1e-14 * work[1]);
    let heat = 0.5 * (impedances[0].re + impedances[1].re) * exact.norm_sq();
    assert!((0.5 * (force.conj() * exact).re - heat).abs() < 1e-14 * heat);
}

/// G1/G3: nonsymmetric original momentum/thermal blocks versus an independent
/// primal solve, under 16 orders of interface-coordinate rescaling. This is a
/// mathematical Schur test; it does not replace the native tetrahedral tests.
#[test]
fn full_vector_and_thermal_connection_matches_primal_under_basis_scaling() {
    let n = 4;
    let matrices: [Vec<C64>; 2] = std::array::from_fn(|side| {
        let mut a = vec![C64::ZERO; n * n];
        for i in 0..n {
            a[i * n + i] = C64::new(1.0 + i as f64 + 2.0 * side as f64, 0.5 + 0.2 * i as f64);
        }
        a[3] = C64::new(0.4, 0.0);
        a[12] = C64::new(-0.4, 0.0);
        a[1] = C64::new(0.2, 0.0);
        a[4] = C64::new(0.2, 0.0);
        a
    });
    let factors = matrices
        .each_ref()
        .map(|a| lu_complex_with_cancel(a, n, &mut || false).unwrap());
    let rhs = [
        vec![
            C64::new(0.7, 0.2),
            C64::new(-0.3, 0.1),
            C64::new(0.2, -0.4),
            C64::new(0.8, 0.1),
        ],
        vec![
            C64::new(0.1, -0.2),
            C64::ZERO,
            C64::new(-0.2, 0.3),
            C64::new(-0.1, 0.5),
        ],
    ];
    let primal: Vec<_> = matrices[0]
        .iter()
        .zip(&matrices[1])
        .map(|(&a, &b)| a + b)
        .collect();
    let mut exact: Vec<_> = rhs[0].iter().zip(&rhs[1]).map(|(&a, &b)| a + b).collect();
    lu_complex_with_cancel(&primal, n, &mut || false)
        .unwrap()
        .solve_with_cancel(&mut exact, &mut || false)
        .unwrap();
    for basis in [[1.0; 4], [1e-8, 1e3, 5.0, 1e8]] {
        let mut left = vec![0.0; n * n];
        for i in 0..n {
            left[i * n + i] = 2e-6 * basis[i];
        }
        left[1] = 0.3e-6 * basis[1];
        left[4] = 0.3e-6 * basis[0];
        let right: Vec<_> = left.iter().map(|v| -v).collect();
        let loads = [&left[..], &right[..]];
        let mut solve = |side: usize, rhs: &[C64], cancelled: &mut _| {
            let mut x = rhs.to_vec();
            factors[side].solve_with_cancel(&mut x, cancelled)?;
            Ok(x)
        };
        let mut cancelled = || false;
        let system =
            ConformingSchur::prepare(loads, [n, n], n, 1e-12, &mut solve, &mut cancelled).unwrap();
        let (fields, reaction) = system
            .response(
                loads,
                3,
                [3, 3],
                [3, 3],
                [&rhs[0], &rhs[1]],
                &mut solve,
                &mut cancelled,
            )
            .unwrap();
        let mut interface_work = [0.0; 2];
        for side in 0..2 {
            for (&a, &b) in fields[side].iter().zip(&exact) {
                near(a, b, 2e-13);
            }
            for column in 0..n {
                interface_work[side] += 0.5
                    * (reaction[column].conj() * dual(loads[side], n, column, &fields[side])).re;
            }
        }
        assert!((interface_work[0] + interface_work[1]).abs() < 1e-13);
    }
}

/// G4: a trace fixed on one side arrests the other without demanding exact
/// zero relative to its own roundoff. The unrelated field remains live.
#[test]
fn one_sided_stationary_constraint_uses_physical_free_field_scale() {
    let impedances = [C64::new(3.0, 8.0), C64::new(0.07, 2.1)];
    let factors = impedances.map(|a| lu_complex_with_cancel(&[a], 1, &mut || false).unwrap());
    let loads = [&[0.0][..], &[-1.3e-7][..]];
    let mut solve = |side: usize, rhs: &[C64], cancelled: &mut _| {
        let mut x = rhs.to_vec();
        factors[side].solve_with_cancel(&mut x, cancelled)?;
        Ok(x)
    };
    let mut cancelled = || false;
    let system =
        ConformingSchur::prepare(loads, [1, 1], 1, 1e-12, &mut solve, &mut cancelled).unwrap();
    let (fields, reaction) = system
        .response(
            loads,
            1,
            [1, 1],
            [1, 1],
            [&[C64::new(0.2, 0.3)], &[C64::new(0.7, -0.9)]],
            &mut solve,
            &mut cancelled,
        )
        .unwrap();
    near(fields[0][0], C64::new(0.2, 0.3) / impedances[0], 2e-14);
    assert!(fields[1][0].abs() < 1e-14);
    near(reaction[0].scale(1.3e-7), C64::new(0.7, -0.9), 2e-14);
}

#[test]
fn cancellation_refuses_the_candidate_and_prepared_connection_remains_reusable() {
    let loads = [&[1.0][..], &[-1.0][..]];
    let mut solve = |_side: usize, rhs: &[C64], cancelled: &mut _| {
        poll(cancelled)?;
        Ok(rhs.to_vec())
    };
    let mut cancelled = || false;
    let system =
        ConformingSchur::prepare(loads, [1, 1], 1, 1e-12, &mut solve, &mut cancelled).unwrap();
    assert!(
        system
            .response(
                loads,
                1,
                [1, 1],
                [1, 1],
                [&[C64::new(1.0, 0.0)], &[C64::ZERO]],
                &mut |_side, rhs, cancelled: &mut fn() -> bool| {
                    poll(cancelled)?;
                    Ok(rhs.to_vec())
                },
                &mut ((|| true) as fn() -> bool)
            )
            .is_err()
    );
    let (fields, _) = system
        .response(
            loads,
            1,
            [1, 1],
            [1, 1],
            [&[C64::new(1.0, 0.0)], &[C64::ZERO]],
            &mut solve,
            &mut cancelled,
        )
        .unwrap();
    assert_eq!(fields[0][0], C64::new(0.5, 0.0));
    assert_eq!(fields[1][0], C64::new(0.5, 0.0));
}
