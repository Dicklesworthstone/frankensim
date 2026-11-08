//! G0/G1/G3/G4 references for finite, activity-consistent ionic equilibrium.
use super::*;
use crate::osmosis::{donnan::*, nonideal::VirialSolutionModel};
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T;

fn model(caps: Vec<f64>, matrix: Vec<f64>) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, 18e-6, 1.0, caps, matrix).unwrap())
}

fn binary(b: f64, caps: [f64; 2]) -> VirialDonnan {
    VirialDonnan::new(
        vec![1, -1],
        model(caps.to_vec(), vec![b, 0.0, 0.0, b]),
        2e-12,
    )
    .unwrap()
}

fn input(totals: &[f64]) -> DonnanProblem<'_> {
    DonnanProblem {
        temperature_k: T,
        volumes_m3: [1.0, 1.0],
        mobile_totals_mol: totals,
        fixed_charge_mol: [-1.0, 0.0],
        trapped_particles_mol: [0.0, 0.0],
    }
}

fn solve(law: &mut VirialDonnan, problem: DonnanProblem<'_>) -> (Vec<f64>, Vec<f64>, DonnanReport) {
    let mut left = vec![0.0; problem.mobile_totals_mol.len()];
    let mut right = left.clone();
    let report = law
        .equilibrate_into(
            problem,
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

fn potential(report: &DonnanReport) -> f64 {
    match report.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            ..
        } => dimensionless_left_minus_right,
        other => panic!("expected finite potential, got {other:?}"),
    }
}

// Independent one-dimensional constrained free-energy minimum. Electroneutrality
// gives nL=[1+t,t], nR=[1-t,1-t], eliminating both charge and voltage unknowns.
fn binary_root(b: f64, volumes: [f64; 2]) -> f64 {
    let [vl, vr] = volumes;
    let f = |t: f64| {
        ((1.0 + t) / vl).ln() + (t / vl).ln() - 2.0 * ((1.0 - t) / vr).ln()
            + b * ((1.0 + 2.0 * t) / vl - 2.0 * (1.0 - t) / vr)
    };
    let (mut lo, mut hi) = (1e-12, 1.0 - 1e-12);
    assert!(f(lo) < 0.0 && f(hi) > 0.0);
    for _ in 0..80 {
        let mid = lo + 0.5 * (hi - lo);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo + 0.5 * (hi - lo)
}

fn binary_pressure(b: f64, volume: f64, n: [f64; 2]) -> f64 {
    let [p, m] = n.map(|x| x / volume);
    RT * (p + m + 0.5 * b * (p * p + m * m))
}

fn binary_energy(b: f64, volumes: [f64; 2], t: f64) -> f64 {
    [[1.0 + t, t], [1.0 - t, 1.0 - t]]
        .into_iter()
        .zip(volumes)
        .map(|(amounts, v)| {
            amounts
                .into_iter()
                .map(|n| n * ((n / v).ln() - 1.0) + 0.5 * b * n * n / v)
                .sum::<f64>()
        })
        .sum::<f64>()
        * RT
}

#[test]
fn g1_zero_interactions_recover_original_ideal_partition_with_absent_species() {
    let valences = vec![1, -1, 2, -2, 1, -1];
    let totals = [3.0, 2.0, 1.0, 0.5, 1.0, 0.0];
    let ideal = IdealDonnan::new(valences.clone()).unwrap();
    let mut virial =
        VirialDonnan::new(valences, model(vec![20.0; 6], vec![0.0; 36]), 2e-12).unwrap();
    let problem = DonnanProblem {
        volumes_m3: [0.7, 1.3],
        fixed_charge_mol: [-2.0, -1.0],
        trapped_particles_mol: [0.1, 0.04],
        ..input(&totals)
    };
    let (mut ideal_left, mut ideal_right) = ([0.0; 6], [0.0; 6]);
    let reference = ideal
        .equilibrate_into(
            problem,
            DonnanControls::default(),
            &mut ideal_left,
            &mut ideal_right,
            || false,
        )
        .unwrap();
    let (left, right, report) = solve(&mut virial, problem);
    for i in 0..6 {
        near(left[i], ideal_left[i], 3e-12);
        near(right[i], ideal_right[i], 3e-12);
    }
    assert_eq!((left[5], right[5]), (0.0, 0.0));
    near(potential(&report), potential(&reference), 3e-12);
    for side in 0..2 {
        near(
            report.osmotic_pressure_pa[side],
            reference.osmotic_pressure_pa[side],
            3e-8,
        );
    }
}

#[test]
fn g1_original_activities_change_partition_voltage_pressure_and_free_energy() {
    let b = 0.8;
    let mut law = binary(b, [10.0; 2]);
    let problem = input(&[2.0, 1.0]);
    let (left, right, report) = solve(&mut law, problem);
    let t = binary_root(b, problem.volumes_m3);
    assert!((t - 1.0 / 3.0).abs() > 0.01);
    for (actual, expected) in left
        .iter()
        .chain(&right)
        .zip([1.0 + t, t, 1.0 - t, 1.0 - t])
    {
        near(*actual, expected, 2e-11);
    }
    let psi = -((1.0 + t) / (1.0 - t)).ln() - 2.0 * b * t;
    near(potential(&report), psi, 3e-11);
    near(
        report.osmotic_pressure_pa[0],
        binary_pressure(b, 1.0, [1.0 + t, t]),
        2e-7,
    );
    near(
        report.osmotic_pressure_pa[1],
        binary_pressure(b, 1.0, [1.0 - t; 2]),
        2e-7,
    );
    let energy = binary_energy(b, problem.volumes_m3, left[1]);
    for alternative in [t - 0.01, t + 0.01, 1.0 / 3.0] {
        assert!(energy < binary_energy(b, problem.volumes_m3, alternative));
    }
    assert!(report.charge_residual_mol.iter().all(|q| q.abs() < 2e-11));
}

#[test]
fn g3_mixed_valence_cross_interactions_obey_each_original_electrochemical_equation() {
    let valences = vec![1, -1, 2, -2, 1];
    let totals = [3.0, 2.0, 1.0, 0.5, 1.0];
    let mut matrix = vec![0.02; 25];
    for i in 0..5 {
        matrix[5 * i + i] = 0.15 + 0.03 * i as f64;
    }
    let mut law = VirialDonnan::new(
        valences.clone(),
        model(vec![20.0; 5], matrix.clone()),
        2e-12,
    )
    .unwrap();
    let problem = DonnanProblem {
        volumes_m3: [0.7, 1.3],
        fixed_charge_mol: [-2.0, -1.0],
        trapped_particles_mol: [0.1, 0.04],
        ..input(&totals)
    };
    let (left, right, report) = solve(&mut law, problem);
    let psi = potential(&report);
    let cl: Vec<f64> = left.iter().map(|n| n / 0.7).collect();
    let cr: Vec<f64> = right.iter().map(|n| n / 1.3).collect();
    for i in 0..5 {
        let excess: f64 = (0..5).map(|j| matrix[5 * i + j] * (cl[j] - cr[j])).sum();
        near(
            cl[i].ln() - cr[i].ln() + excess + f64::from(valences[i]) * psi,
            0.0,
            2e-10,
        );
        near(left[i] + right[i], totals[i], 2e-15);
    }
    for (side, c) in [&cl, &cr].into_iter().enumerate() {
        let quadratic: f64 = (0..5)
            .map(|i| (0..5).map(|j| c[i] * matrix[5 * i + j] * c[j]).sum::<f64>())
            .sum();
        let pressure = RT
            * (c.iter().sum::<f64>()
                + 0.5 * quadratic
                + problem.trapped_particles_mol[side] / problem.volumes_m3[side]);
        near(report.osmotic_pressure_pa[side], pressure, 2e-7);
    }
    let reversed = DonnanProblem {
        volumes_m3: [1.3, 0.7],
        fixed_charge_mol: [-1.0, -2.0],
        trapped_particles_mol: [0.04, 0.1],
        ..problem
    };
    let (reverse_left, reverse_right, reverse) = solve(&mut law, reversed);
    for i in 0..5 {
        near(reverse_left[i], right[i], 2e-10);
        near(reverse_right[i], left[i], 2e-10);
    }
    near(potential(&reverse), -psi, 2e-10);
}

#[test]
fn g1_reduced_osmotic_stiffness_matches_independent_partition_sensitivity() {
    let b = 0.8;
    let mut law = binary(b, [10.0; 2]);
    let problem = DonnanProblem {
        volumes_m3: [0.8, 1.2],
        trapped_particles_mol: [0.1, 0.04],
        ..input(&[2.0, 1.0])
    };
    let (_, _, report) = solve(&mut law, problem);
    let pressure_drop = |x: f64| {
        let volumes = [0.8 - x, 1.2 + x];
        let t = binary_root(b, volumes);
        binary_pressure(b, volumes[0], [1.0 + t, t]) - binary_pressure(b, volumes[1], [1.0 - t; 2])
            + RT * (0.1 / volumes[0] - 0.04 / volumes[1])
    };
    let h = 1e-5;
    let reference = (pressure_drop(h) - pressure_drop(-h)) / (2.0 * h);
    let actual = report.osmotic_pressure_drop_slope_pa_per_m3.unwrap();
    assert!(actual > 0.0);
    near(actual, reference, 2e-6 * reference.abs());
}

#[test]
fn g0_supported_nonideal_root_is_found_when_ideal_seed_exceeds_original_cap() {
    let mut law = binary(0.8, [1.32, 1.0]);
    let (left, right, _) = solve(&mut law, input(&[2.0, 1.0]));
    let expected = binary_root(0.8, [1.0, 1.0]);
    assert!(4.0 / 3.0 > law.model().maximum_concentrations_mol_m3()[0]);
    near(left[0], 1.0 + expected, 2e-11);
    near(left[1], expected, 2e-11);
    assert!(left[0] <= 1.32 && right[0] <= 1.32);
}

#[test]
fn g0_counterion_only_limit_keeps_exact_zeros_and_original_nonideal_pressure() {
    let mut law = binary(0.8, [10.0; 2]);
    let problem = DonnanProblem {
        trapped_particles_mol: [0.2, 0.03],
        ..input(&[1.0, 0.0])
    };
    let (left, right, report) = solve(&mut law, problem);
    assert_eq!(left, [1.0, 0.0]);
    assert_eq!(right, [0.0, 0.0]);
    assert_eq!(report.potential, DonnanPotential::NegativeInfinite);
    near(report.osmotic_pressure_pa[0], RT * 1.6, 2e-9);
    near(report.osmotic_pressure_pa[1], RT * 0.03, 2e-10);
    near(
        report.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
        RT * 2.03,
        2e-8,
    );
}

#[test]
fn g4_domain_work_and_cancellation_refusals_leave_outputs_unchanged_and_reusable() {
    let problem = input(&[2.0, 1.0]);
    let mut left = [17.0, 19.0];
    let mut right = [-23.0, -29.0];
    let mut unsupported = binary(0.8, [0.9, 0.9]);
    assert!(
        unsupported
            .equilibrate_into(
                problem,
                DonnanControls::default(),
                &mut left,
                &mut right,
                || false,
            )
            .is_err()
    );
    assert_eq!((left, right), ([17.0, 19.0], [-23.0, -29.0]));
    let mut law = binary(0.8, [10.0; 2]);
    assert!(
        law.equilibrate_into(
            problem,
            DonnanControls::default(),
            &mut left,
            &mut right,
            || true,
        )
        .is_err()
    );
    assert_eq!((left, right), ([17.0, 19.0], [-23.0, -29.0]));
    assert!(
        law.equilibrate_into(
            problem,
            DonnanControls {
                maximum_iterations: 1,
                ..Default::default()
            },
            &mut left,
            &mut right,
            || false,
        )
        .is_err()
    );
    assert_eq!((left, right), ([17.0, 19.0], [-23.0, -29.0]));
    law.equilibrate_into(
        problem,
        DonnanControls::default(),
        &mut left,
        &mut right,
        || false,
    )
    .unwrap();
    near(left[1], binary_root(0.8, [1.0, 1.0]), 2e-11);
}

#[test]
fn g3_volume_and_amount_scaling_preserve_trace_ions_and_original_pressure_sensitivity() {
    let totals = [2.0, 1.0, 1e-24, 1e-24];
    let mut reference: Option<(Vec<f64>, Vec<f64>, DonnanReport)> = None;
    for (volume_scale, amount_scale) in [(1.0, 1.0), (1e-6, 1e-4), (1e6, 1e-4)] {
        let concentration_scale = amount_scale / volume_scale;
        let mut matrix = vec![0.0; 16];
        for i in 0..4 {
            matrix[4 * i + i] = 0.8 / concentration_scale;
        }
        let mut law = VirialDonnan::new(
            vec![1, -1, 2, -2],
            model(vec![10.0 * concentration_scale; 4], matrix),
            2e-12,
        )
        .unwrap();
        let scaled_totals = totals.map(|n| n * amount_scale);
        let problem = DonnanProblem {
            volumes_m3: [0.8 * volume_scale, 1.2 * volume_scale],
            fixed_charge_mol: [-amount_scale, 0.0],
            trapped_particles_mol: [0.1 * amount_scale, 0.02 * amount_scale],
            ..input(&scaled_totals)
        };
        let (left, right, report) = solve(&mut law, problem);
        assert!(left[2] > 0.0 && right[2] > 0.0 && left[3] > 0.0 && right[3] > 0.0);
        if let Some((ref_left, ref_right, ref_report)) = &reference {
            for i in 0..4 {
                near(left[i] / amount_scale, ref_left[i], 5e-11 * totals[i]);
                near(right[i] / amount_scale, ref_right[i], 5e-11 * totals[i]);
            }
            near(potential(&report), potential(ref_report), 5e-11);
            for side in 0..2 {
                near(
                    report.osmotic_pressure_pa[side] / concentration_scale,
                    ref_report.osmotic_pressure_pa[side],
                    5e-7,
                );
            }
            let slope = report.osmotic_pressure_drop_slope_pa_per_m3.unwrap() * volume_scale
                / concentration_scale;
            near(
                slope,
                ref_report.osmotic_pressure_drop_slope_pa_per_m3.unwrap(),
                1e-6,
            );
        } else {
            reference = Some((left, right, report));
        }
    }
}
