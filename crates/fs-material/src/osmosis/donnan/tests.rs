use super::*;

const TEMPERATURE: f64 = 298.15;
const RT: f64 = R * TEMPERATURE;

fn problem(totals: &[f64]) -> DonnanProblem<'_> {
    DonnanProblem {
        temperature_k: TEMPERATURE,
        volumes_m3: [1.0, 1.0],
        mobile_totals_mol: totals,
        fixed_charge_mol: [-1.0, 0.0],
        trapped_particles_mol: [0.0, 0.0],
    }
}
fn solve(law: &IdealDonnan, input: DonnanProblem<'_>) -> (Vec<f64>, Vec<f64>, DonnanReport) {
    let mut left = vec![0.0; law.valences().len()];
    let mut right = left.clone();
    let report = law
        .equilibrate_into(
            input,
            DonnanControls::default(),
            &mut left,
            &mut right,
            || false,
        )
        .unwrap();
    (left, right, report)
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn psi(report: DonnanReport) -> f64 {
    match report.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            ..
        } => dimensionless_left_minus_right,
        other => panic!("expected finite potential, received {other:?}"),
    }
}

#[test]
fn g1_finite_bath_reference_has_exact_partition_potential_and_osmotic_force() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let (left, right, report) = solve(&law, problem(&[2.0, 1.0]));
    for (actual, expected) in
        left.iter()
            .chain(&right)
            .zip([4.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0])
    {
        near(*actual, expected, 2e-13);
    }
    near(psi(report), -2.0_f64.ln(), 2e-13);
    if let DonnanPotential::Finite {
        voltage_left_minus_right_v,
        ..
    } = report.potential
    {
        near(
            voltage_left_minus_right_v,
            -RT * 2.0_f64.ln() / FARADAY_C_PER_MOL,
            1e-14,
        );
    }
    near(
        report.osmotic_pressure_pa[0] - report.osmotic_pressure_pa[1],
        RT / 3.0,
        2e-9,
    );
    near(
        report.solvent_driving_pressure_pa([0.0, 0.0]).unwrap(),
        -RT / 3.0,
        2e-9,
    );
    near(
        report.solvent_driving_pressure_pa([RT / 3.0, 0.0]).unwrap(),
        0.0,
        2e-9,
    );
    assert!(report.charge_residual_mol.iter().all(|q| q.abs() < 2e-13));
}

#[test]
fn g1_uncharged_compartments_equalize_every_ionic_concentration() {
    let law = IdealDonnan::new(vec![1, 2, -1]).unwrap();
    let input = DonnanProblem {
        volumes_m3: [1.0, 3.0],
        fixed_charge_mol: [0.0, 0.0],
        ..problem(&[2.0, 1.0, 4.0])
    };
    let (left, right, report) = solve(&law, input);
    near(psi(report), 0.0, 1e-14);
    for i in 0..3 {
        near(left[i], input.mobile_totals_mol[i] / 4.0, 2e-14);
        near(left[i], right[i] / 3.0, 2e-14);
    }
    near(
        report.osmotic_pressure_pa[0],
        report.osmotic_pressure_pa[1],
        2e-10,
    );
    near(
        report.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        0.0,
        1e-20,
    );
}

#[test]
fn g3_multivalent_mixture_preserves_charge_species_and_electrochemical_potentials() {
    let law = IdealDonnan::new(vec![1, -1, 2, -2, 1]).unwrap();
    let input = DonnanProblem {
        volumes_m3: [0.7, 1.3],
        fixed_charge_mol: [-2.0, -1.0],
        ..problem(&[3.0, 2.0, 1.0, 0.5, 1.0])
    };
    let (left, right, report) = solve(&law, input);
    let potential = psi(report);
    for (i, &z) in law.valences().iter().enumerate() {
        near(left[i] + right[i], input.mobile_totals_mol[i], 2e-15);
        let chemical_difference = (left[i] / 0.7).ln() - (right[i] / 1.3).ln();
        near(chemical_difference + f64::from(z) * potential, 0.0, 3e-13);
    }
    near(left[0], 3.0 * left[4], 1e-14);
    near(right[0], 3.0 * right[4], 1e-14);
    assert!(report.charge_residual_mol.iter().all(|q| q.abs() < 2e-12));
    assert!(report.osmotic_pressure_drop_slope_pa_per_m3.unwrap() >= 0.0);
}

#[test]
fn g3_reversing_compartments_reverses_potential_and_swaps_partition() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let input = DonnanProblem {
        volumes_m3: [0.6, 1.4],
        trapped_particles_mol: [0.1, 0.03],
        ..problem(&[2.0, 1.0])
    };
    let (left, right, a) = solve(&law, input);
    let reversed = DonnanProblem {
        volumes_m3: [1.4, 0.6],
        fixed_charge_mol: [0.0, -1.0],
        trapped_particles_mol: [0.03, 0.1],
        ..input
    };
    let (rleft, rright, b) = solve(&law, reversed);
    for i in 0..2 {
        near(left[i], rright[i], 3e-13);
        near(right[i], rleft[i], 3e-13);
    }
    near(psi(a), -psi(b), 3e-13);
    near(a.osmotic_pressure_pa[0], b.osmotic_pressure_pa[1], 2e-9);
    near(
        a.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        b.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        3e-9,
    );
}

#[test]
fn g1_charge_and_osmotic_particle_counts_are_independent() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let (left, right, a) = solve(&law, problem(&[2.0, 1.0]));
    let (bleft, bright, b) = solve(
        &law,
        DonnanProblem {
            trapped_particles_mol: [0.02, 0.03],
            ..problem(&[2.0, 1.0])
        },
    );
    assert_eq!((left, right), (bleft, bright));
    near(
        b.osmotic_pressure_pa[0] - a.osmotic_pressure_pa[0],
        0.02 * RT,
        1e-9,
    );
    near(
        b.osmotic_pressure_pa[1] - a.osmotic_pressure_pa[1],
        0.03 * RT,
        1e-9,
    );
    near(
        b.osmotic_pressure_drop_slope_pa_per_m3.unwrap()
            - a.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        0.05 * RT,
        1e-9,
    );
}

#[test]
fn g1_extensive_scaling_preserves_concentrations_and_potential() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let (left, right, a) = solve(&law, problem(&[2.0, 1.0]));
    let scaled = DonnanProblem {
        volumes_m3: [7.0, 7.0],
        fixed_charge_mol: [-7.0, 0.0],
        ..problem(&[14.0, 7.0])
    };
    let (bleft, bright, b) = solve(&law, scaled);
    for i in 0..2 {
        near(bleft[i], 7.0 * left[i], 3e-12);
        near(bright[i], 7.0 * right[i], 3e-12);
    }
    near(psi(a), psi(b), 3e-13);
    near(a.osmotic_pressure_pa[0], b.osmotic_pressure_pa[0], 2e-9);
    near(
        a.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        7.0 * b.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        3e-9,
    );
}

#[test]
fn g1_counterion_only_limits_and_pure_solvent_need_no_invented_salt() {
    for z in [-1, 1] {
        let law = IdealDonnan::new(vec![z]).unwrap();
        let input = DonnanProblem {
            fixed_charge_mol: [-2.0 * f64::from(z), 0.0],
            ..problem(&[2.0])
        };
        let (left, right, report) = solve(&law, input);
        assert_eq!(
            (left.as_slice(), right.as_slice()),
            (&[2.0][..], &[0.0][..])
        );
        assert_eq!(
            report.potential,
            if z > 0 {
                DonnanPotential::NegativeInfinite
            } else {
                DonnanPotential::PositiveInfinite
            }
        );
        assert_eq!(report.charge_residual_mol, [0.0, 0.0]);
        near(report.osmotic_pressure_pa[0], 2.0 * RT, 1e-12);
        near(
            report.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
            2.0 * RT,
            1e-12,
        );
    }
    let empty = IdealDonnan::new(vec![]).unwrap();
    let (_, _, report) = solve(
        &empty,
        DonnanProblem {
            fixed_charge_mol: [0.0, 0.0],
            volumes_m3: [1.0, 2.0],
            trapped_particles_mol: [0.1, 0.2],
            ..problem(&[])
        },
    );
    assert_eq!(report.potential, DonnanPotential::Indeterminate);
    assert_eq!(report.iterations, 0);
    near(
        report.osmotic_pressure_pa[0],
        report.osmotic_pressure_pa[1],
        1e-12,
    );
}

#[test]
fn g0_non_neutral_and_locally_infeasible_inventories_refuse_atomically() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let mut left = [41.0, 42.0];
    let mut right = [43.0, 44.0];
    let unmatched = law.equilibrate_into(
        problem(&[1.0, 1.0]),
        DonnanControls::default(),
        &mut left,
        &mut right,
        || false,
    );
    assert!(matches!(
        unmatched,
        Err(DonnanError::NotElectroneutral { .. })
    ));
    let impossible = DonnanProblem {
        fixed_charge_mol: [2.0, -2.0],
        ..problem(&[1.0, 1.0])
    };
    assert_eq!(
        law.equilibrate_into(
            impossible,
            DonnanControls::default(),
            &mut left,
            &mut right,
            || false
        )
        .unwrap_err(),
        DonnanError::InfeasibleCharge
    );
    assert_eq!((left, right), ([41.0, 42.0], [43.0, 44.0]));
}

#[test]
fn g0_domain_layout_and_nonfinite_results_are_checked_before_publication() {
    assert!(IdealDonnan::new(vec![0]).is_err());
    assert!(IdealDonnan::new(vec![1; MAX_DONNAN_SPECIES + 1]).is_err());
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let mut left = [41.0, 42.0];
    let mut right = [43.0, 44.0];
    for input in [
        DonnanProblem {
            temperature_k: 0.0,
            ..problem(&[2.0, 1.0])
        },
        DonnanProblem {
            volumes_m3: [0.0, 1.0],
            ..problem(&[2.0, 1.0])
        },
        DonnanProblem {
            fixed_charge_mol: [f64::NAN, 0.0],
            ..problem(&[2.0, 1.0])
        },
        DonnanProblem {
            trapped_particles_mol: [-1.0, 0.0],
            ..problem(&[2.0, 1.0])
        },
        problem(&[-2.0, 1.0]),
        problem(&[f64::INFINITY, 1.0]),
        problem(&[2.0]),
        DonnanProblem {
            volumes_m3: [1e-310, 1e-310],
            ..problem(&[2.0, 1.0])
        },
    ] {
        assert!(
            law.equilibrate_into(
                input,
                DonnanControls::default(),
                &mut left,
                &mut right,
                || false
            )
            .is_err()
        );
        assert_eq!((left, right), ([41.0, 42.0], [43.0, 44.0]));
    }
}

#[test]
fn g0_work_exhaustion_and_late_cancellation_leave_both_outputs_unchanged() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let input = problem(&[2.0, 1.0]);
    let (mut left, mut right) = ([41.0, 42.0], [43.0, 44.0]);
    let budget = DonnanControls {
        maximum_iterations: 1,
        ..DonnanControls::default()
    };
    assert_eq!(
        law.equilibrate_into(input, budget, &mut left, &mut right, || false)
            .unwrap_err(),
        DonnanError::NonConvergence
    );
    assert_eq!((left, right), ([41.0, 42.0], [43.0, 44.0]));
    let mut polls_to_completion = 0;
    law.equilibrate_into(
        input,
        DonnanControls::default(),
        &mut left,
        &mut right,
        || {
            polls_to_completion += 1;
            false
        },
    )
    .unwrap();
    (left, right) = ([41.0, 42.0], [43.0, 44.0]);
    let mut polls = 0;
    let result = law.equilibrate_into(
        input,
        DonnanControls::default(),
        &mut left,
        &mut right,
        || {
            polls += 1;
            polls == polls_to_completion
        },
    );
    assert_eq!(result.unwrap_err(), DonnanError::Cancelled);
    assert_eq!((left, right), ([41.0, 42.0], [43.0, 44.0]));
}

#[test]
fn g1_logistic_retains_tiny_representable_minority_without_a_floor() {
    for odds in [-1000.0, 1000.0] {
        let (left, right) = split_amount(1e300, odds);
        let expected = (1e300_f64.ln() - 1000.0).exp();
        let minority = left.min(right);
        assert!(minority > 0.0 && minority.is_finite());
        near(minority / expected, 1.0, 2e-12);
        assert_eq!(left + right, 1e300);
    }
    assert_eq!(split_amount(0.0, 1000.0), (0.0, 0.0));
}

#[test]
fn g1_equilibrium_minimizes_ideal_mixing_energy_and_has_positive_analytic_water_stiffness() {
    let law = IdealDonnan::new(vec![1, -1]).unwrap();
    let (left, right, report) = solve(&law, problem(&[2.0, 1.0]));
    let energy = |values: &[f64]| {
        values
            .iter()
            .map(|&n| if n == 0.0 { 0.0 } else { n * (n.ln() - 1.0) })
            .sum::<f64>()
    };
    assert!(energy(&left) + energy(&right) < energy(&[1.0, 0.0]) + energy(&[1.0, 1.0]));
    let analytic = report.osmotic_pressure_drop_slope_pa_per_m3.unwrap();
    near(analytic, RT * 17.0 / 27.0, 2e-9);
    let h = 1e-5;
    let (_, _, plus) = solve(
        &law,
        DonnanProblem {
            volumes_m3: [1.0 - h, 1.0 + h],
            ..problem(&[2.0, 1.0])
        },
    );
    let (_, _, minus) = solve(
        &law,
        DonnanProblem {
            volumes_m3: [1.0 + h, 1.0 - h],
            ..problem(&[2.0, 1.0])
        },
    );
    let numeric = ((plus.osmotic_pressure_pa[0] - plus.osmotic_pressure_pa[1])
        - (minus.osmotic_pressure_pa[0] - minus.osmotic_pressure_pa[1]))
        / (2.0 * h);
    near(analytic, numeric, 1e-5);
}
