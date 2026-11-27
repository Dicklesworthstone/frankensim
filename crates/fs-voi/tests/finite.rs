//! Exact finite-menu information values and refusal boundaries.
use fs_voi::finite::{
    FINITE_EVSI_VERSION, FiniteDecision, FiniteError, FiniteTest, MASS_TOLERANCE, MAX_AXIS, evsi,
    value_now,
};

#[test]
fn two_actions_two_hypotheses_by_hand() {
    let d = FiniteDecision {
        prior: &[0.5, 0.5],
        payoff: &[&[6., 0.], &[0., 4.]],
    };
    let test = FiniteTest {
        likelihood: &[&[0.75, 0.25], &[0.25, 0.75]],
        cost: 1.,
    };
    assert_eq!(value_now(&d), Ok((0, 3.)));
    let result = evsi(&d, &test).unwrap();
    // Joint best utilities are 2.25 and 1.5. Perfect information yields 5.
    assert_eq!(result.evsi, 0.75);
    assert_eq!(result.evpi, 2.);
    assert_eq!(result.net_value, -0.25);
    assert_eq!(result.posterior_policy, [Some(0), Some(1)]);
    assert_eq!(result.outcome_probability, [0.5, 0.5]);
    assert_eq!(result, evsi(&d, &test).unwrap());
    assert_eq!(FINITE_EVSI_VERSION, 1);
}

#[test]
fn action_consequence_and_impossible_outcomes() {
    let d = FiniteDecision {
        prior: &[0.6, 0.4],
        payoff: &[&[-1., -2.], &[-1., -2.], &[-2., -5.]],
    };
    let test = FiniteTest {
        likelihood: &[&[0.1, 0.9, 0.], &[0.6, 0.4, 0.]],
        cost: 0.,
    };
    let r = evsi(&d, &test).unwrap();
    assert_eq!((r.evsi, r.evpi), (0., 0.));
    assert_eq!(r.posterior_policy, [Some(0), Some(0), None]);

    let d = FiniteDecision {
        prior: &[0.9, 0.1],
        payoff: &[&[2., 0.], &[0., 2.]],
    };
    let r = evsi(&d, &test).unwrap();
    assert_eq!(r.evsi.to_bits(), 0);
    assert!(r.evpi > 0.);
    assert_eq!(r.posterior_policy, [Some(0), Some(0), None]);
    let noisy = FiniteTest {
        likelihood: &[&[0.3, 0.7], &[0.3, 0.7]],
        cost: 0.,
    };
    assert_eq!(evsi(&d, &noisy).unwrap().evsi.to_bits(), 0);
}

#[test]
fn perfect_information_permutations_and_dominated_actions() {
    let d = FiniteDecision {
        prior: &[0.25, 0.75],
        payoff: &[&[3., -1.], &[0., 2.]],
    };
    let test = FiniteTest {
        likelihood: &[&[0., 1.], &[1., 0.]],
        cost: 0.,
    };
    let r = evsi(&d, &test).unwrap();
    assert_eq!(r.evsi, 0.75);
    assert_eq!(r.evsi.to_bits(), r.evpi.to_bits());
    assert_eq!(r.posterior_policy, [Some(1), Some(0)]);
    let d = FiniteDecision {
        prior: &[0.75, 0.25],
        payoff: &[&[-9., -9.], &[2., 0.], &[-1., 3.]],
    };
    let permuted = evsi(&d, &test).unwrap();
    assert_eq!(r.evsi, permuted.evsi);
    assert_eq!(permuted.posterior_policy, [Some(2), Some(1)]);
}

#[test]
fn agrees_with_integer_joint_payoffs_across_menus() {
    // The prior and each likelihood have denominator eight. Compute every
    // expectation independently in integer numerators with denominator 64.
    let mut seed = 1_u32;
    let mut draw = |limit| {
        seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        (seed >> 8) as i32 % limit
    };
    for _ in 0..512 {
        let h0 = draw(9);
        let h1 = draw(9 - h0);
        let p = [h0, h1, 8 - h0 - h1];
        let u = [[(); 3]; 4].map(|row| row.map(|()| draw(33) - 16));
        let l = [(); 3].map(|()| {
            let v = draw(9);
            [v, 8 - v]
        });
        let prior = p.map(|v| f64::from(v) / 8.);
        let payoff = u.map(|row| row.map(f64::from));
        let likelihood = l.map(|row| row.map(|v| f64::from(v) / 8.));
        let d = FiniteDecision {
            prior: &prior,
            payoff: &payoff.each_ref().map(|r| r.as_slice()),
        };
        let t = FiniteTest {
            likelihood: &likelihood.each_ref().map(|r| r.as_slice()),
            cost: 0.0,
        };
        let now = u
            .iter()
            .map(|row| (0..3).map(|h| p[h] * row[h]).sum::<i32>())
            .max()
            .unwrap();
        let perfect: i32 = (0..3)
            .map(|h| p[h] * u.iter().map(|row| row[h]).max().unwrap())
            .sum();
        let mut after = 0;
        let r = evsi(&d, &t).unwrap();
        for y in 0..2 {
            let score = u.map(|row| (0..3).map(|h| p[h] * l[h][y] * row[h]).sum::<i32>());
            let best = *score.iter().max().unwrap();
            let mass: i32 = (0..3).map(|h| p[h] * l[h][y]).sum();
            assert_eq!(r.outcome_probability[y], f64::from(mass) / 64.);
            if mass == 0 {
                assert_eq!(r.posterior_policy[y], None);
            } else {
                assert_eq!(r.posterior_policy[y], score.iter().position(|x| *x == best));
            }
            after += best;
        }
        assert!((r.evsi - f64::from(after - 8 * now) / 64.).abs() < 1e-13);
        assert_eq!(r.evpi, f64::from(perfect - now) / 8.);
        assert!(r.evsi >= 0. && r.evsi <= r.evpi);
    }
}

#[test]
fn validates_all_rows_and_axis_bounds() {
    let decide = |p: &[f64]| {
        value_now(&FiniteDecision {
            prior: p,
            payoff: &[&[1., 2.]],
        })
    };
    assert_eq!(decide(&[]), Err(FiniteError::InvalidShape));
    for p in [
        [0., 0.],
        [-0.1, 1.],
        [f64::NAN, 1.],
        [0.5, 0.5 + 2. * MASS_TOLERANCE],
    ] {
        assert_eq!(decide(&p), Err(FiniteError::InvalidPrior));
    }
    assert_eq!(
        value_now(&FiniteDecision {
            prior: &[1.],
            payoff: &[&[f64::INFINITY]]
        }),
        Err(FiniteError::NonFinitePayoff)
    );
    let too_wide = vec![0.; MAX_AXIS + 1];
    assert_eq!(
        value_now(&FiniteDecision {
            prior: &too_wide,
            payoff: &[&too_wide]
        }),
        Err(FiniteError::InvalidShape)
    );
    let d = FiniteDecision {
        prior: &[1., 0.],
        payoff: &[&[1., 2.]],
    };
    for (rows, error) in [
        (&[][..], FiniteError::InvalidShape),
        (&[&[1.][..], &[0.5, 0.5]][..], FiniteError::InvalidShape),
        (&[&[1.][..], &[0.8]][..], FiniteError::InvalidLikelihood(1)),
    ] {
        assert_eq!(
            evsi(
                &d,
                &FiniteTest {
                    likelihood: rows,
                    cost: 0.
                }
            ),
            Err(error)
        );
    }
    for cost in [-1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            evsi(
                &d,
                &FiniteTest {
                    likelihood: &[&[1.], &[1.]],
                    cost
                }
            ),
            Err(FiniteError::InvalidCost)
        );
    }
    assert_eq!(
        value_now(&FiniteDecision {
            prior: &[1. - MASS_TOLERANCE / 2.],
            payoff: &[&[2.]]
        }),
        Ok((0, 2.))
    );
}

#[test]
fn rare_observation_ignores_zero_prior_and_keeps_its_policy() {
    let tiny = f64::from_bits(2);
    let d = FiniteDecision {
        prior: &[0., 0.5, 0.5],
        payoff: &[&[0., 4., 0.], &[100., 0., 2.]],
    };
    let t = FiniteTest {
        likelihood: &[&[1., 0.], &[0., 1.], &[tiny, 1.]],
        cost: 0.,
    };
    let r = evsi(&d, &t).unwrap();
    assert_eq!(r.outcome_probability[0].to_bits(), 1);
    assert_eq!(r.posterior_policy[0], Some(1));
    assert_eq!(r.evsi.to_bits(), 2);
}
